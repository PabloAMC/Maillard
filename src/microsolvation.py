"""Generalised microsolvation utilities for transition-state search.

Provides:
- ``place_proton_shuttle_water``: place a single explicit water bridging a
  donor → H → acceptor proton transfer (turns a 3-membered TS into a
  5-membered ring, lowering the barrier ~50 kcal/mol for 1,3-H shifts).
- ``kabsch_align``: heavy-atom RMSD alignment of two geometries.
- ``relax_solvent``: xTB+ALPB optimisation with a frozen solute mask.
- ``build_hydrated_endpoints``: ties everything together for an elementary
  step (donor, acceptor, h_atom indices) → hydrated reactant/product
  geometries with byte-identical water positions.

Design notes (proven on the 2026-04-20 pe_amadori work):
- Waters MUST occupy identical positions in R and P, otherwise DE-GSM tries
  to evolve water reorganisation along the string and the resulting "TS"
  has zero imaginary frequencies.
- The dry product geometry is first Kabsch-aligned onto the dry reactant
  heavy-atom frame so transplanted waters do not collide with shifted heavy
  atoms.
- Optional CREST/QCG cluster-continuum sampling exists in ``src.solvation``
  but is overkill for single proton-shuttle waters; the analytic placement
  here is faster and reproducible.
"""
from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Ideal isolated-water geometry parameters (TIP3P-ish)
OH_BOND_A = 0.97
HOH_ANGLE_DEG = 104.5
DEFAULT_WATER_OFFSET_A = 1.6
CLASH_THRESHOLD_A = 0.85


# ---------------------------------------------------------------------------
# XYZ helpers
# ---------------------------------------------------------------------------

def read_xyz(path: Path) -> Tuple[List[str], np.ndarray]:
    """Read a single-frame XYZ file. Returns (symbols, coords[N,3])."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    n = int(lines[0].strip())
    syms: List[str] = []
    coords: List[List[float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return syms, np.asarray(coords, dtype=float)


def write_xyz(path: Path, syms: Sequence[str], coords: np.ndarray, comment: str = "") -> None:
    lines = [str(len(syms)), comment]
    for s, (x, y, z) in zip(syms, coords):
        lines.append(f"{s:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------

def kabsch_align(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Optimal rotation+translation aligning ``mobile`` onto ``target``.

    Returns ``(R, mobile_centroid, target_centroid)`` such that the aligned
    coordinates are ``(mobile - mobile_centroid) @ R + target_centroid``.
    """
    if mobile.shape != target.shape:
        raise ValueError("kabsch_align: shape mismatch")
    mc = mobile.mean(axis=0)
    tc = target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _s, vh = np.linalg.svd(h)
    d = float(np.sign(np.linalg.det(vh.T @ u.T)))
    # We want ``aligned = (mobile - mc) @ rot + tc`` to match ``target``.
    # The Kabsch optimal rotation (col-vec convention) is
    # ``R = V @ diag(1,1,d) @ U.T``; in row-vector form we need ``R.T``.
    rot = u @ np.diag([1.0, 1.0, d]) @ vh
    return rot, mc, tc


def molecular_normal(coords: np.ndarray) -> np.ndarray:
    """Estimate a plane normal via SVD of centered coords."""
    centered = coords - coords.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centered, full_matrices=False)
    n = vh[-1]
    return n / np.linalg.norm(n)


def make_water(o_pos: np.ndarray, donor_pos: np.ndarray, acceptor_pos: np.ndarray) -> np.ndarray:
    """Build a 3-atom water (O, H_acceptor_side, H_donor_side).

    The water is oriented so its hydrogen-bond donor H points toward the
    proton acceptor, and its free H points back toward the proton donor.
    Geometry is approximate; xTB relaxation refines it.
    """
    midpoint = 0.5 * (donor_pos + acceptor_pos)
    bisector = midpoint - o_pos
    bnorm = float(np.linalg.norm(bisector))
    da_vec = acceptor_pos - donor_pos
    da_vec /= np.linalg.norm(da_vec)
    if bnorm < 1e-6:
        # Degenerate: O lies on the donor-acceptor midpoint. Pick an
        # arbitrary axis perpendicular to the donor-acceptor line as the
        # H-H bisector. The water plane will be perpendicular to da_vec.
        ref = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(ref, da_vec)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        bisector = ref - np.dot(ref, da_vec) * da_vec
        bisector /= np.linalg.norm(bisector)
    else:
        bisector /= bnorm
    perp = da_vec - np.dot(da_vec, bisector) * bisector
    perp /= np.linalg.norm(perp)

    half_angle = np.deg2rad(HOH_ANGLE_DEG / 2.0)
    h_acc_dir = np.cos(half_angle) * bisector + np.sin(half_angle) * perp
    h_don_dir = np.cos(half_angle) * bisector - np.sin(half_angle) * perp
    return np.array([
        o_pos,
        o_pos + OH_BOND_A * h_acc_dir,
        o_pos + OH_BOND_A * h_don_dir,
    ])


def place_proton_shuttle_water(
    coords: np.ndarray,
    donor_idx: int,
    acceptor_idx: int,
    *,
    reference_normal: np.ndarray | None = None,
    offset: float = DEFAULT_WATER_OFFSET_A,
) -> np.ndarray:
    """Place 1 explicit water bridging donor → acceptor (5-mem TS topology).

    The water O is placed ``offset`` Å above the donor-acceptor midpoint along
    ``reference_normal`` (computed from the molecular plane if not provided).

    Returns ``(3, 3)`` array with [O, H_accept_side, H_donor_side] coords.
    """
    donor = coords[donor_idx]
    acceptor = coords[acceptor_idx]
    if reference_normal is None:
        reference_normal = molecular_normal(coords)
    midpoint = 0.5 * (donor + acceptor)
    o_pos = midpoint + offset * reference_normal
    return make_water(o_pos, donor, acceptor)


# ---------------------------------------------------------------------------
# Clash detection
# ---------------------------------------------------------------------------

@dataclass
class ClashReport:
    min_distance_a: float
    pair: Tuple[int, str, int, str] | None
    has_clash: bool


def check_no_clash(
    syms: Sequence[str],
    coords: np.ndarray,
    n_solute: int,
    *,
    threshold: float = CLASH_THRESHOLD_A,
) -> ClashReport:
    """Compute the minimum nonbonded distance, skipping intra-water O–H bonds."""
    n = len(coords)
    min_d = float("inf")
    pair = None
    for i in range(n):
        for j in range(i + 1, n):
            # Skip intra-water bonds (consecutive triplets after solute block)
            if i >= n_solute and j >= n_solute:
                if (i - n_solute) // 3 == (j - n_solute) // 3:
                    continue
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if d < min_d:
                min_d = d
                pair = (i, syms[i], j, syms[j])
    return ClashReport(min_distance_a=min_d, pair=pair, has_clash=min_d < threshold)


# ---------------------------------------------------------------------------
# xTB constrained optimisation (waters relaxed, solute frozen)
# ---------------------------------------------------------------------------

def relax_solvent(
    geom_path: Path,
    n_solute: int,
    work_dir: Path,
    *,
    charge: int = 0,
    solvent: str = "water",
    timeout_s: int = 600,
) -> Path:
    """Run xTB --opt with the first ``n_solute`` atoms frozen.

    Returns path to the optimised xyz (in ``work_dir/xtbopt.xyz``).
    Raises ``RuntimeError`` if xTB fails or does not produce output.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    constrain = work_dir / "constrain.inp"
    constrain.write_text(
        "$constrain\n"
        "  force constant=1.0\n"
        f"  atoms: 1-{n_solute}\n"
        "$end\n",
        encoding="utf-8",
    )
    local_geom = work_dir / Path(geom_path).name
    if Path(geom_path).resolve() != local_geom.resolve():
        shutil.copy(geom_path, local_geom)
    logger.info(
        "[microsolvation] xTB constrained opt on %s (atoms 1-%d frozen)",
        local_geom.name, n_solute,
    )
    cmd = [
        "xtb", local_geom.name,
        "--opt", "tight",
        "--input", "constrain.inp",
        "--gfn", "2",
        "--alpb", solvent,
        "--chrg", str(charge),
    ]
    result = subprocess.run(
        cmd, cwd=work_dir, capture_output=True, text=True, timeout=timeout_s,
    )
    if result.returncode != 0:
        logger.error("[microsolvation] xTB failed:\n%s", result.stderr[-2000:])
        raise RuntimeError(f"xTB constrained opt failed for {local_geom.name}")
    out = work_dir / "xtbopt.xyz"
    if not out.exists():
        raise RuntimeError(f"xTB did not produce xtbopt.xyz for {local_geom.name}")
    return out


# ---------------------------------------------------------------------------
# High-level: build hydrated R/P endpoints for an elementary proton-shuttle step
# ---------------------------------------------------------------------------

@dataclass
class HydratedEndpoints:
    reactant_xyz_path: Path
    product_xyz_path: Path
    n_solute: int
    n_total: int
    heavy_rmsd_after_align: float
    clash_report_reactant: ClashReport
    clash_report_product: ClashReport


def build_hydrated_endpoints(
    dry_reactant_xyz: Path,
    dry_product_xyz: Path,
    work_dir: Path,
    *,
    proton_shuttles: Sequence[Tuple[int, int]],
    charge: int = 0,
    solvent: str = "water",
    offset: float = DEFAULT_WATER_OFFSET_A,
    relax_product_waters: bool = True,
) -> HydratedEndpoints:
    """Build hydrated R/P endpoints with ``len(proton_shuttles)`` waters.

    Each entry of ``proton_shuttles`` is ``(donor_idx, acceptor_idx)`` (0-based,
    matching atom order in the dry XYZ files). One explicit water is added
    per shuttle, bridging donor and acceptor.

    Steps:
      1. Read dry R, P (must have identical atom count and order).
      2. Kabsch-align P onto R using heavy atoms.
      3. Place ``len(proton_shuttles)`` waters near R coordinates.
      4. Constrained xTB opt: relax waters in R, solute frozen.
      5. Transplant the relaxed waters into the aligned P.
      6. (Optional, default True) Constrained xTB opt: relax waters in P,
         solute frozen. Recommended whenever P heavy atoms differ
         appreciably from R (RMSD >~ 0.5 Å); otherwise the unrelaxed water
         can leave P higher in energy than the TS, causing DE-GSM to place
         the peak at the endpoint.
      7. Sanity-check clashes in both R and P.
      8. Write reactant.xyz, product.xyz to ``work_dir``.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    syms_r, coords_r = read_xyz(dry_reactant_xyz)
    syms_p, coords_p_raw = read_xyz(dry_product_xyz)
    n_solute = len(syms_r)
    if len(syms_p) != n_solute:
        raise ValueError(
            f"reactant/product atom count mismatch: R={n_solute} P={len(syms_p)}"
        )

    # Kabsch align P onto R
    heavy_idx = [i for i, s in enumerate(syms_r) if s != "H"]
    rot, p_centroid, r_centroid = kabsch_align(coords_p_raw[heavy_idx], coords_r[heavy_idx])
    coords_p = (coords_p_raw - p_centroid) @ rot + r_centroid
    rmsd = float(np.sqrt(((coords_p[heavy_idx] - coords_r[heavy_idx]) ** 2).sum() / len(heavy_idx)))
    logger.info("[microsolvation] heavy-atom RMSD after Kabsch align: %.3f Å", rmsd)

    # Place waters from R reference
    reactive_atoms_idx: List[int] = []
    for d, a in proton_shuttles:
        reactive_atoms_idx.extend([d, a])
    normal = molecular_normal(coords_r[reactive_atoms_idx])
    centroid_r = coords_r.mean(axis=0)
    if np.dot(normal, coords_r[reactive_atoms_idx].mean(axis=0) - centroid_r) < 0:
        normal = -normal

    waters: List[np.ndarray] = []
    new_syms: List[str] = list(syms_r)
    for donor, acceptor in proton_shuttles:
        w = place_proton_shuttle_water(
            coords_r, donor, acceptor, reference_normal=normal, offset=offset,
        )
        waters.append(w)
        new_syms.extend(["O", "H", "H"])
    coords_rw = np.vstack([coords_r] + waters)

    init_path = work_dir / "reactant_initial.xyz"
    write_xyz(init_path, new_syms, coords_rw, comment="hydrated reactant (initial)")

    # Relax waters in R
    opt_path = relax_solvent(init_path, n_solute, work_dir, charge=charge, solvent=solvent)
    syms_rr, coords_rr = read_xyz(opt_path)
    react_path = work_dir / "reactant.xyz"
    write_xyz(react_path, syms_rr, coords_rr, comment="hydrated reactant (waters relaxed)")

    # Transplant relaxed waters into aligned P
    coords_pw = np.vstack([coords_p, coords_rr[n_solute:]])
    syms_pw: List[str] = list(syms_p) + list(syms_rr[n_solute:])
    prod_initial = work_dir / "product_initial.xyz"
    write_xyz(prod_initial, syms_pw, coords_pw,
              comment="hydrated product (waters from R, before P relax)")

    if relax_product_waters:
        prod_opt = relax_solvent(
            prod_initial, n_solute, work_dir,
            charge=charge, solvent=solvent,
        )
        syms_pw_r, coords_pw_r = read_xyz(prod_opt)
        # Keep solute coordinates from the aligned P (xTB freeze should
        # leave them invariant, but enforce it explicitly to be safe).
        coords_pw_r[:n_solute] = coords_p
        coords_pw = coords_pw_r
        syms_pw = list(syms_pw_r)
    prod_path = work_dir / "product.xyz"
    write_xyz(prod_path, syms_pw, coords_pw,
              comment="hydrated product" + (" (waters relaxed)" if relax_product_waters else " (waters from R)"))

    # Sanity checks
    clash_r = check_no_clash(syms_rr, coords_rr, n_solute)
    clash_p = check_no_clash(syms_pw, coords_pw, n_solute)
    if clash_r.has_clash:
        logger.warning("[microsolvation] reactant has clash %.3f Å at %s",
                       clash_r.min_distance_a, clash_r.pair)
    if clash_p.has_clash:
        logger.warning("[microsolvation] product has clash %.3f Å at %s",
                       clash_p.min_distance_a, clash_p.pair)

    # Clean intermediate xTB scratch
    for name in ("xtbopt.log", "xtbtopo.mol", "charges", "wbo", "xtbrestart",
                 "constrain.inp", "xtbopt.xyz", ".xtbtmpmol",
                 "reactant_initial.xyz", "product_initial.xyz"):
        f = work_dir / name
        if f.exists():
            f.unlink()

    return HydratedEndpoints(
        reactant_xyz_path=react_path,
        product_xyz_path=prod_path,
        n_solute=n_solute,
        n_total=len(syms_pw),
        heavy_rmsd_after_align=rmsd,
        clash_report_reactant=clash_r,
        clash_report_product=clash_p,
    )
