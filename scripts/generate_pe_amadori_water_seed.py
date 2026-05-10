#!/usr/bin/env python3
"""
Generate an explicit-water cluster seed for the PE-Amadori 1,2-proton shift.

Diagnostic findings (2026-04-20) on the dry pe_amadori target:
- The Amadori rearrangement is two coupled 1,3-H shifts:
    O0–H17 → C2 (enolisation: O loses H, C gains H)
    C1–H18 → N3 (proton hop on backbone)
- Donor–acceptor distances are 2.39 Å (O0···C2) and 2.52 Å (C1···N3),
  giving 3-membered transition states with prohibitive strain (~70 kcal/mol).
- Inserting water as a Grotthuss relay opens 5-membered TS rings, lowering
  the barrier by ~50 kcal/mol toward the literature value (Ea ≈ 19.8 kcal/mol).

This script:
1. Reads dry pe_amadori reactant/product geometries (33 atoms each).
2. Places 2 explicit waters bridging the two proton-transfer channels.
3. Relaxes water positions with xTB+ALPB(water) while the solute is frozen.
4. Builds an analogous product geometry with waters in identical relative
   positions (so DE-GSM has a clean R/P pair with the same water molecules).
5. Writes `reactant.xyz`, `product.xyz`, and `run_xtb.sh` into
   `data/geometries/xtb_inputs/pe_amadori_water/`.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "data" / "geometries" / "xtb_inputs" / "pe_amadori"
DEST_DIR = ROOT / "data" / "geometries" / "xtb_inputs" / "pe_amadori_water"

# Donor / acceptor indices (0-based), confirmed by geometric analysis:
#   O0 (donor)  -> C2 (acceptor)   transfer of H17
#   C1 (donor)  -> N3 (acceptor)   transfer of H18
DONOR_O, ACCEPTOR_C = 0, 2
DONOR_C, ACCEPTOR_N = 1, 3

# Water geometry parameters
OH_BOND = 0.97  # Å
HOH_ANGLE_DEG = 104.5
WATER_OFFSET = 1.6  # Å above the donor-acceptor midpoint (perpendicular to plane)


def read_xyz(path: Path) -> Tuple[List[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    n = int(lines[0].strip())
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return syms, np.array(coords)


def write_xyz(path: Path, syms: List[str], coords: np.ndarray, comment: str = "") -> None:
    lines = [str(len(syms)), comment]
    for s, (x, y, z) in zip(syms, coords):
        lines.append(f"{s:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_water(o_pos: np.ndarray, donor_pos: np.ndarray,
               acceptor_pos: np.ndarray) -> np.ndarray:
    """Build a 3-atom water (O, H_donor_side, H_acceptor_side).

    H pointing toward acceptor: hydrogen-bond donor of water.
    H pointing toward donor: free hydrogen (will accept H-bond from donor).

    The geometry is approximate; xTB relaxation refines it.
    """
    # Bisector vector points from O toward midpoint(donor, acceptor)
    midpoint = (donor_pos + acceptor_pos) / 2.0
    bisector = midpoint - o_pos
    bisector /= np.linalg.norm(bisector)

    # Perpendicular axis (in the plane defined by O, donor, acceptor)
    da_vec = acceptor_pos - donor_pos
    da_vec /= np.linalg.norm(da_vec)
    perp = da_vec - np.dot(da_vec, bisector) * bisector
    perp /= np.linalg.norm(perp)

    half_angle = np.deg2rad(HOH_ANGLE_DEG / 2.0)
    # H pointing toward acceptor side
    h_acc_dir = np.cos(half_angle) * bisector + np.sin(half_angle) * perp
    h_acc = o_pos + OH_BOND * h_acc_dir
    # H pointing toward donor side (opposite perp component)
    h_don_dir = np.cos(half_angle) * bisector - np.sin(half_angle) * perp
    h_don = o_pos + OH_BOND * h_don_dir

    return np.array([o_pos, h_acc, h_don])


def place_water_bridge(donor: np.ndarray, acceptor: np.ndarray,
                       reference_normal: np.ndarray) -> np.ndarray:
    """Place a water O at midpoint of donor-acceptor + offset along normal.

    Returns 3 atoms: O, H_acceptor_side, H_donor_side.
    """
    midpoint = (donor + acceptor) / 2.0
    o_pos = midpoint + WATER_OFFSET * reference_normal
    return make_water(o_pos, donor, acceptor)


def molecular_normal(coords: np.ndarray) -> np.ndarray:
    """Estimate a molecular plane normal using SVD of centered heavy-atom coords."""
    centered = coords - coords.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
    # Smallest singular vector = plane normal
    return vh[-1] / np.linalg.norm(vh[-1])


def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute optimal rotation+translation that aligns ``mobile`` onto ``target``.

    Returns (rotation_matrix, mobile_centroid, target_centroid). The aligned
    coordinates are: ``(mobile - mobile_centroid) @ R + target_centroid``.
    Both inputs must have identical shape ``(N, 3)``.
    """
    if mobile.shape != target.shape:
        raise ValueError("kabsch_align: shape mismatch")
    mc = mobile.mean(axis=0)
    tc = target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _s, vh = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vh.T @ u.T))
    rot = vh.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, mc, tc


def build_water_augmented_geometry(syms: List[str], coords: np.ndarray) -> Tuple[List[str], np.ndarray]:
    """Insert 2 waters bridging the two proton-transfer channels."""
    # Use a normal computed from the reactive heavy atoms to place water on the
    # accessible face. We pick the normal to the plane (O0, C1, C2, N3).
    reactive_atoms = coords[[DONOR_O, DONOR_C, ACCEPTOR_C, ACCEPTOR_N]]
    normal = molecular_normal(reactive_atoms)
    # Ensure it points away from N3's other substituents (heuristic)
    centroid = coords.mean(axis=0)
    if np.dot(normal, reactive_atoms.mean(axis=0) - centroid) < 0:
        normal = -normal

    # Water 1: bridges O0 (donor) -> C2 (acceptor)  for H17
    water1 = place_water_bridge(coords[DONOR_O], coords[ACCEPTOR_C], normal)
    # Water 2: bridges C1 (donor) -> N3 (acceptor)  for H18
    water2 = place_water_bridge(coords[DONOR_C], coords[ACCEPTOR_N], normal)

    new_syms = list(syms) + ["O", "H", "H", "O", "H", "H"]
    new_coords = np.vstack([coords, water1, water2])
    return new_syms, new_coords


def relax_waters(geom_path: Path, n_solute: int, work_dir: Path) -> Path:
    """Run xTB optimisation with the solute frozen; return path to xtbopt.xyz."""
    constrain = work_dir / "constrain.inp"
    # Atoms 1-N (1-indexed) frozen
    constrain.write_text(
        "$constrain\n"
        "  force constant=1.0\n"
        f"  atoms: 1-{n_solute}\n"
        "$end\n",
        encoding="utf-8",
    )
    # Copy geometry into work_dir
    local_geom = work_dir / geom_path.name
    if local_geom.resolve() != geom_path.resolve():
        shutil.copy(geom_path, local_geom)
    print(f"  [xTB] Relaxing {geom_path.name} (solute atoms 1-{n_solute} frozen)...")
    cmd = [
        "xtb", local_geom.name,
        "--opt", "tight",
        "--input", "constrain.inp",
        "--gfn", "2",
        "--alpb", "water",
        "--chrg", "0",
    ]
    result = subprocess.run(
        cmd, cwd=work_dir, capture_output=True, text=True, timeout=600
    )
    if result.returncode != 0:
        print("  [xTB] FAILED. Stderr:")
        print(result.stderr[-2000:])
        raise RuntimeError(f"xTB optimisation failed for {geom_path.name}")
    out = work_dir / "xtbopt.xyz"
    if not out.exists():
        raise RuntimeError(f"xTB did not produce xtbopt.xyz for {geom_path.name}")
    return out


def cleanup(work_dir: Path) -> None:
    for name in [
        "xtbopt.log", "xtbtopo.mol", "charges", "wbo", "xtbrestart",
        "constrain.inp", "xtbopt.xyz", ".xtbtmpmol",
        "reactant_initial.xyz", "product_initial.xyz",
    ]:
        f = work_dir / name
        if f.exists():
            f.unlink()


def main() -> int:
    print("=== PE-Amadori explicit water seed generator ===")
    if not SRC_DIR.exists():
        print(f"ERROR: source directory not found: {SRC_DIR}")
        return 1
    react_src = SRC_DIR / "reactant.xyz"
    prod_src = SRC_DIR / "product.xyz"
    if not react_src.exists() or not prod_src.exists():
        print(f"ERROR: missing reactant.xyz or product.xyz in {SRC_DIR}")
        return 1

    DEST_DIR.mkdir(parents=True, exist_ok=True)

    # Read dry geometries (33 atoms each)
    syms_r, coords_r = read_xyz(react_src)
    syms_p, coords_p_raw = read_xyz(prod_src)
    n_solute = len(syms_r)
    if len(syms_p) != n_solute:
        print(f"ERROR: atom count mismatch R={n_solute} P={len(syms_p)}")
        return 1
    print(f"Loaded dry pe_amadori: {n_solute} atoms each (R, P)")

    # Align dry product onto dry reactant frame using a Kabsch fit on heavy
    # atoms. The dry product geometry has arbitrary absolute coordinates from
    # RDKit; we want both endpoints in the same frame so that the transplanted
    # waters do not collide with shifted heavy atoms.
    heavy_idx = [i for i, s in enumerate(syms_r) if s != "H"]
    rot, p_centroid, r_centroid = kabsch_align(
        coords_p_raw[heavy_idx], coords_r[heavy_idx]
    )
    coords_p = (coords_p_raw - p_centroid) @ rot + r_centroid
    rmsd_after = float(np.sqrt(((coords_p[heavy_idx] - coords_r[heavy_idx]) ** 2).sum() / len(heavy_idx)))
    print(f"Aligned dry product to dry reactant (heavy-atom RMSD={rmsd_after:.3f} Å)")

    # Build hydrated REACTANT geometry (place 2 waters using R coordinates).
    print("Placing 2 explicit waters bridging proton-transfer channels...")
    syms_rw, coords_rw = build_water_augmented_geometry(syms_r, coords_r)

    # Write initial (un-relaxed) reactant guess
    react_init = DEST_DIR / "reactant_initial.xyz"
    write_xyz(react_init, syms_rw, coords_rw, comment="pe_amadori + 2 H2O (initial)")

    # Relax reactant waters with the solute frozen (xTB+ALPB)
    print()
    print("Relaxing reactant water positions (solute frozen)...")
    react_opt = relax_waters(react_init, n_solute, DEST_DIR)
    syms_rr, coords_rr = read_xyz(react_opt)
    write_xyz(DEST_DIR / "reactant.xyz", syms_rr, coords_rr,
              comment="pe_amadori + 2 H2O (waters relaxed, solute frozen)")

    # Build hydrated PRODUCT by transplanting the relaxed waters from reactant.
    # CRITICAL: the waters MUST be byte-identical between R and P. If we let
    # xTB relax the product waters independently, they drift to a different
    # local minimum (per-atom diffs up to 5 Å observed 2026-04-20), and DE-GSM
    # then has to evolve water reorganisation along the string in addition to
    # the proton transfers — yielding a "TS" that is not a true saddle.
    coords_pw = np.vstack([coords_p, coords_rr[n_solute:]])
    syms_pw = list(syms_p) + list(syms_rr[n_solute:])

    # Sanity check: no severe clashes after the transplant. Heavy R↔P atom
    # displacement is bounded by Kabsch alignment but we still want to verify.
    coords_pw = np.asarray(coords_pw)
    min_d = float("inf")
    min_pair = None
    for i in range(len(coords_pw)):
        for j in range(i + 1, len(coords_pw)):
            # Skip water O–H bonds (within each water, atoms are consecutive)
            if i >= n_solute and j >= n_solute and (j - i) <= 2 and (i - n_solute) // 3 == (j - n_solute) // 3:
                continue
            d = float(np.linalg.norm(coords_pw[i] - coords_pw[j]))
            if d < min_d:
                min_d = d
                min_pair = (i, syms_pw[i], j, syms_pw[j])
    print(f"Product transplant min nonbonded distance: {min_d:.3f} Å at {min_pair}")
    if min_d < 0.85:
        print(f"WARNING: clash detected (<0.85 Å); GSM may struggle.")

    write_xyz(DEST_DIR / "product.xyz", syms_pw, coords_pw,
              comment="pe_amadori product (aligned) + 2 H2O (waters identical to R, no relax)")
    syms_pr, coords_pr = syms_pw, coords_pw

    # Write run_xtb.sh for downstream pipeline
    run_sh = DEST_DIR / "run_xtb.sh"
    run_sh.write_text(
        "#!/bin/bash\n"
        "# DE-GSM path search for PE-Amadori with explicit waters\n"
        "xtb reactant.xyz --path product.xyz --gfn 2 --alpb water --chrg 0\n",
        encoding="utf-8",
    )
    run_sh.chmod(0o755)

    cleanup(DEST_DIR)
    print()
    print(f"✓ Wrote hydrated geometries to {DEST_DIR}")
    print(f"  reactant.xyz: {len(syms_rr)} atoms")
    print(f"  product.xyz : {len(syms_pr)} atoms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
