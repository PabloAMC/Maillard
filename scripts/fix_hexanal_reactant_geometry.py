#!/usr/bin/env python3
"""Regenerate the hexanal_radical_quench reactant geometry.

The historical ``reactant.xyz`` ships with the [SH] thiyl radical literally
sitting on top of the hexanal backbone (S–C5 = 0.70 Å, S–O = 0.62 Å), an
artefact of the ``source.kind == "manual"`` branch in
``scripts/generators/integrate_pubchem_geometries.py::_refresh_reactant``,
which never re-aligns the manual fragment relative to the host molecule.

This module rebuilds the fragment in a chemically sensible HAT pre-reactive
geometry:

* place the migrating ``H`` along the carbonyl-oxygen lone-pair direction
  (the C→O axis), at d(O,H) = 1.5 Å (well above the bonded 0.98 Å so the
  bond is "forming" but not yet bonded);
* place the sulphur further along the same axis, at d(O,S) = 2.8 Å, which
  preserves d(S,H) ≈ 1.30 Å (close to the equilibrium S–H length 1.34 Å).

The hexanal heavy-atom skeleton is left untouched.  A backup of the
previous file is written to
``results/computational_gap_refinement/pubchem_geometry_refresh/hexanal_radical_quench/``
so the change is traceable and reversible.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TARGET_DIR = ROOT / "data" / "geometries" / "xtb_inputs" / "hexanal_radical_quench"
BACKUP_DIR = (
    ROOT
    / "results"
    / "computational_gap_refinement"
    / "pubchem_geometry_refresh"
    / "hexanal_radical_quench"
)

# Indices in the canonical reactant.xyz layout (verified against the manifest):
#   0..5   carbon chain (C5 = aldehyde carbon)
#   6      aldehyde O
#   7..18  hydrogens of hexanal
#   19     S (manual fragment)
#   20     H (manual fragment, the migrating proton)
ALDEHYDE_C_INDEX = 5
ALDEHYDE_O_INDEX = 6
SULFUR_INDEX = 19
MIGRATING_H_INDEX = 20

# Pre-reactive HAT contact distances (Å).
TARGET_OH_DISTANCE = 1.50
TARGET_OS_DISTANCE = 2.80


def _read_xyz(path: Path) -> Tuple[List[str], np.ndarray, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1]
    symbols: List[str] = []
    coords: List[List[float]] = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.asarray(coords, dtype=float), comment


def _write_xyz(path: Path, symbols: List[str], coords: np.ndarray, comment: str) -> None:
    out = [f"{len(symbols)}", comment]
    for sym, (x, y, z) in zip(symbols, coords):
        out.append(f"{sym:<2s}  {x:>12.6f}  {y:>12.6f}  {z:>12.6f}")
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _validate_layout(symbols: List[str]) -> None:
    expected = {
        ALDEHYDE_C_INDEX: "C",
        ALDEHYDE_O_INDEX: "O",
        SULFUR_INDEX: "S",
        MIGRATING_H_INDEX: "H",
    }
    for idx, sym in expected.items():
        if symbols[idx] != sym:
            raise RuntimeError(
                f"Unexpected layout: atom {idx} is {symbols[idx]!r}, expected {sym!r}. "
                "Refusing to mutate this geometry."
            )


def _reposition_sh(coords: np.ndarray) -> np.ndarray:
    new_coords = coords.copy()
    c = new_coords[ALDEHYDE_C_INDEX]
    o = new_coords[ALDEHYDE_O_INDEX]
    axis = o - c
    axis /= np.linalg.norm(axis)
    new_coords[MIGRATING_H_INDEX] = o + TARGET_OH_DISTANCE * axis
    new_coords[SULFUR_INDEX] = o + TARGET_OS_DISTANCE * axis
    return new_coords


def _check_no_clashes(symbols: List[str], coords: np.ndarray) -> None:
    """Fail loudly if any heavy/heavy distance is below 1.4 Å."""
    movers = (SULFUR_INDEX, MIGRATING_H_INDEX)
    for mover in movers:
        for j, sym in enumerate(symbols):
            if j == mover:
                continue
            if j in movers:
                continue  # S–H pair is checked separately
            dist = float(np.linalg.norm(coords[mover] - coords[j]))
            min_allowed = 1.4 if sym != "H" else 1.5
            if dist < min_allowed:
                raise RuntimeError(
                    f"Clash after repositioning: atom {mover} ({symbols[mover]}) is "
                    f"{dist:.3f} Å from atom {j} ({sym}); minimum allowed {min_allowed} Å."
                )
    sh = float(np.linalg.norm(coords[SULFUR_INDEX] - coords[MIGRATING_H_INDEX]))
    if not 1.20 <= sh <= 1.45:
        raise RuntimeError(f"Unexpected S–H distance after repositioning: {sh:.3f} Å.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the new geometry diagnostics but do not write files.")
    args = parser.parse_args()

    reactant_path = TARGET_DIR / "reactant.xyz"
    symbols, coords, comment = _read_xyz(reactant_path)
    _validate_layout(symbols)

    o = coords[ALDEHYDE_O_INDEX]
    print("Before:")
    print(f"  d(O, H_mig) = {np.linalg.norm(o - coords[MIGRATING_H_INDEX]):.3f} Å")
    print(f"  d(O, S)     = {np.linalg.norm(o - coords[SULFUR_INDEX]):.3f} Å")
    print(f"  d(S, C5)    = {np.linalg.norm(coords[SULFUR_INDEX] - coords[ALDEHYDE_C_INDEX]):.3f} Å")

    new_coords = _reposition_sh(coords)
    _check_no_clashes(symbols, new_coords)

    o_new = new_coords[ALDEHYDE_O_INDEX]
    print("After:")
    print(f"  d(O, H_mig) = {np.linalg.norm(o_new - new_coords[MIGRATING_H_INDEX]):.3f} Å")
    print(f"  d(O, S)     = {np.linalg.norm(o_new - new_coords[SULFUR_INDEX]):.3f} Å")
    print(f"  d(S, C5)    = {np.linalg.norm(new_coords[SULFUR_INDEX] - new_coords[ALDEHYDE_C_INDEX]):.3f} Å")
    print(f"  d(S, H_mig) = {np.linalg.norm(new_coords[SULFUR_INDEX] - new_coords[MIGRATING_H_INDEX]):.3f} Å")

    if args.dry_run:
        print("\n--dry-run: no files written.")
        return 0

    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    backup_path = BACKUP_DIR / "reactant.pre_sh_relocation.xyz"
    shutil.copy2(reactant_path, backup_path)
    print(f"\nBacked up old geometry → {backup_path.relative_to(ROOT)}")

    new_comment = (
        "hexanal + [SH] HAT pre-reactive geometry; "
        "SH placed along C=O axis (d(O,H)=1.5 Å, d(O,S)=2.8 Å)"
    )
    _write_xyz(reactant_path, symbols, new_coords, new_comment)
    print(f"Wrote new geometry → {reactant_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
