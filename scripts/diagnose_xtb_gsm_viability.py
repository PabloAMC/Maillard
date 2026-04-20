#!/usr/bin/env python3
"""Diagnose whether xTB captures a usable barrier before integrating pyGSM.

This script evaluates a small panel of reactions with:
1. A 1D relaxed scan along the manifest forming bond.
2. A 2D constrained scan along the forming bond plus the next-largest
   pairwise distance change between reactant and product.

The goal is not to produce production-grade barriers, but to answer a gating
question: does xTB show a positive saddle-like region on a coarse reaction
surface at all?
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from ase.constraints import FixBondLength
from ase.optimize import BFGS

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src._native_io import suppress_native_output as _suppress_native_output
from src.xtb_backend import (
    _xyz_to_ase_atoms,
    get_xtb_ase_calculator,
    run_xtb_relaxed_scan,
)

DEFAULT_TARGETS = [
    "hexanal_radical_quench",
    "lysinoalanine_crosslink",
    "aa_ring_open_dicarbonyl",
]
EV_TO_KCAL_MOL = 23.0605
MIN_SCAN_DISTANCE = 1.2


def _load_manifest(path: Path) -> Dict[str, Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    jobs = payload.get("jobs", [])
    return {str(job.get("target_id")): job for job in jobs}


def _pair_metadata(atoms_r, atoms_p, i_atom: int, j_atom: int) -> Dict[str, Any]:
    d_reactant = float(atoms_r.get_distance(i_atom, j_atom))
    d_product = float(atoms_p.get_distance(i_atom, j_atom))
    symbols = [atoms_r[i_atom].symbol, atoms_r[j_atom].symbol]
    return {
        "pair": [int(i_atom), int(j_atom)],
        "symbols": symbols,
        "reactant_distance_angstrom": d_reactant,
        "product_distance_angstrom": d_product,
        "delta_angstrom": d_product - d_reactant,
        "abs_delta_angstrom": abs(d_product - d_reactant),
        "heavy_heavy": all(symbol != "H" for symbol in symbols),
    }


def _pick_second_coordinate(atoms_r, atoms_p, pair1: Tuple[int, int]) -> Optional[Dict[str, Any]]:
    candidates: List[Dict[str, Any]] = []
    n_atoms = len(atoms_r)
    pair1_set = {int(pair1[0]), int(pair1[1])}
    for i_atom in range(n_atoms):
        for j_atom in range(i_atom + 1, n_atoms):
            if {i_atom, j_atom} == pair1_set:
                continue
            meta = _pair_metadata(atoms_r, atoms_p, i_atom, j_atom)
            if meta["abs_delta_angstrom"] < 0.5:
                continue
            score = meta["abs_delta_angstrom"]
            if min(meta["reactant_distance_angstrom"], meta["product_distance_angstrom"]) < 2.2:
                score += 0.5
            if meta["heavy_heavy"]:
                score += 0.25
            meta["selection_score"] = score
            candidates.append(meta)
    if not candidates:
        return None
    candidates.sort(
        key=lambda item: (
            item["selection_score"],
            item["abs_delta_angstrom"],
            item["heavy_heavy"],
        ),
        reverse=True,
    )
    return candidates[0]


def _scan_axis(dist_reactant: float, dist_product: float, n_points: int) -> List[float]:
    dist_stop = max(float(dist_product), MIN_SCAN_DISTANCE)
    if n_points <= 1:
        return [float(dist_reactant)]
    return np.linspace(float(dist_reactant), dist_stop, int(n_points)).tolist()


def _run_2d_scan(
    reactant_xyz: str,
    charge: int,
    spin: int,
    pair1: Tuple[int, int],
    pair2: Tuple[int, int],
    atoms_p,
    *,
    grid_points: int = 6,
    fmax: float = 0.25,
    max_steps: int = 15,
) -> Dict[str, Any]:
    atoms_r = _xyz_to_ase_atoms(reactant_xyz)
    axis1 = _scan_axis(
        atoms_r.get_distance(pair1[0], pair1[1]),
        atoms_p.get_distance(pair1[0], pair1[1]),
        grid_points,
    )
    axis2 = _scan_axis(
        atoms_r.get_distance(pair2[0], pair2[1]),
        atoms_p.get_distance(pair2[0], pair2[1]),
        grid_points,
    )
    calc = get_xtb_ase_calculator(charge=charge, spin=spin, solvent="water")
    energies: List[List[Optional[float]]] = []
    failures = 0
    for d1 in axis1:
        row: List[Optional[float]] = []
        for d2 in axis2:
            try:
                atoms = atoms_r.copy()
                atoms.calc = calc
                atoms.set_distance(pair1[0], pair1[1], float(d1), fix=0)
                atoms.set_distance(pair2[0], pair2[1], float(d2), fix=0)
                atoms.set_constraint([
                    FixBondLength(pair1[0], pair1[1]),
                    FixBondLength(pair2[0], pair2[1]),
                ])
                with _suppress_native_output():
                    opt = BFGS(atoms, logfile=None)
                    opt.run(fmax=fmax, steps=max_steps)
                    row.append(float(atoms.get_potential_energy()))
            except Exception:
                failures += 1
                row.append(None)
        energies.append(row)

    start_energy = energies[0][0] if energies and energies[0] else None
    rel_grid: List[List[Optional[float]]] = []
    if start_energy is not None:
        for row in energies:
            rel_grid.append([
                None if energy is None else float(energy - start_energy)
                for energy in row
            ])
    else:
        rel_grid = [[None for _ in row] for row in energies]

    barrier_ev = _minimax_barrier(rel_grid)
    valid_rel = [value for row in rel_grid for value in row if value is not None]
    global_max_ev = max(valid_rel) if valid_rel else None
    return {
        "grid_points": int(grid_points),
        "axis1": axis1,
        "axis2": axis2,
        "energies_ev": energies,
        "relative_energies_ev": rel_grid,
        "start_energy_ev": start_energy,
        "failed_points": int(failures),
        "barrier_2d_ev": barrier_ev,
        "barrier_2d_kcal_mol": None if barrier_ev is None else barrier_ev * EV_TO_KCAL_MOL,
        "grid_global_max_ev": global_max_ev,
        "grid_global_max_kcal_mol": None if global_max_ev is None else global_max_ev * EV_TO_KCAL_MOL,
    }


def _minimax_barrier(rel_grid: List[List[Optional[float]]]) -> Optional[float]:
    if not rel_grid or not rel_grid[0] or rel_grid[0][0] is None:
        return None
    n_rows = len(rel_grid)
    n_cols = len(rel_grid[0])
    inf = float("inf")
    cost = [[inf for _ in range(n_cols)] for _ in range(n_rows)]
    cost[0][0] = max(0.0, float(rel_grid[0][0]))
    for i_row in range(n_rows):
        for j_col in range(n_cols):
            current = rel_grid[i_row][j_col]
            if current is None:
                continue
            if i_row == 0 and j_col == 0:
                continue
            candidates: List[float] = []
            if i_row > 0 and cost[i_row - 1][j_col] != inf:
                candidates.append(max(cost[i_row - 1][j_col], float(current)))
            if j_col > 0 and cost[i_row][j_col - 1] != inf:
                candidates.append(max(cost[i_row][j_col - 1], float(current)))
            if candidates:
                cost[i_row][j_col] = min(candidates)
    final_cost = cost[-1][-1]
    if final_cost == inf:
        return None
    return float(final_cost)


def _summarize_1d(scan_result: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not scan_result:
        return {
            "available": False,
            "barrier_1d_kcal_mol": None,
            "peak_index": None,
            "peak_distance": None,
            "peak_is_interior": False,
            "scan_point_count": 0,
        }
    points = list(scan_result.get("scan_points") or [])
    if not points:
        return {
            "available": False,
            "barrier_1d_kcal_mol": None,
            "peak_index": None,
            "peak_distance": None,
            "peak_is_interior": False,
            "scan_point_count": 0,
        }
    energies = [float(point["energy_ev"]) for point in points]
    peak_index = int(max(range(len(energies)), key=lambda idx: energies[idx]))
    barrier_ev = energies[peak_index] - energies[0]
    return {
        "available": True,
        "barrier_1d_kcal_mol": float(barrier_ev * EV_TO_KCAL_MOL),
        "peak_index": peak_index,
        "peak_distance": float(points[peak_index]["distance"]),
        "peak_is_interior": 0 < peak_index < (len(points) - 1),
        "scan_point_count": len(points),
        "scan_points": points,
    }


def _evaluate_target(job: Dict[str, Any]) -> Dict[str, Any]:
    reactant_path = ROOT / str(job["reactant_path"])
    product_path = ROOT / str(job["product_path"])
    reactant_xyz = reactant_path.read_text(encoding="utf-8")
    product_xyz = product_path.read_text(encoding="utf-8")
    charge = int(job.get("charge", 0) or 0)
    spin = int(job.get("spin", 0) or 0)
    pair1_list = list(job.get("forming_bond_atoms") or [])
    if len(pair1_list) != 2:
        raise RuntimeError(f"{job.get('target_id')}: missing forming_bond_atoms")
    pair1 = (int(pair1_list[0]), int(pair1_list[1]))
    atoms_r = _xyz_to_ase_atoms(reactant_xyz)
    atoms_p = _xyz_to_ase_atoms(product_xyz)
    pair1_meta = _pair_metadata(atoms_r, atoms_p, pair1[0], pair1[1])
    pair2_meta = _pick_second_coordinate(atoms_r, atoms_p, pair1)

    with _suppress_native_output():
        scan_1d = run_xtb_relaxed_scan(
            reactant_xyz,
            product_xyz,
            charge=charge,
            spin=spin,
            forming_bond=pair1,
            n_points=24,
            fmax=0.20,
            max_steps_per_point=20,
        )
    summary_1d = _summarize_1d(scan_1d)

    scan_2d: Dict[str, Any]
    if pair2_meta is None:
        scan_2d = {
            "available": False,
            "reason": "no_second_coordinate",
            "barrier_2d_kcal_mol": None,
            "grid_global_max_kcal_mol": None,
        }
    else:
        pair2 = tuple(pair2_meta["pair"])
        scan_2d = _run_2d_scan(
            reactant_xyz,
            charge,
            spin,
            pair1,
            pair2,
            atoms_p,
            grid_points=6,
            fmax=0.25,
            max_steps=15,
        )
        scan_2d["available"] = True

    support_1d = bool(
        summary_1d.get("available")
        and summary_1d.get("peak_is_interior")
        and (summary_1d.get("barrier_1d_kcal_mol") or 0.0) >= 3.0
    )
    support_2d = bool((scan_2d.get("barrier_2d_kcal_mol") or 0.0) >= 3.0)
    return {
        "target_id": job.get("target_id"),
        "charge": charge,
        "spin": spin,
        "pair1": pair1_meta,
        "pair2": pair2_meta,
        "scan_1d": summary_1d,
        "scan_2d": scan_2d,
        "supports_gsm": bool(support_1d or support_2d),
        "support_reason": {
            "scan_1d": support_1d,
            "scan_2d": support_2d,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="results/computational_gap_refinement/computational_gap_dft_job_manifest.json",
    )
    parser.add_argument(
        "--output",
        default="results/validation/xtb_gsm_viability_2026_04_19.json",
    )
    parser.add_argument("--targets", nargs="*", default=DEFAULT_TARGETS)
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.WARNING)
    logging.getLogger("tblite").setLevel(logging.ERROR)

    manifest = _load_manifest(ROOT / args.manifest)
    results: List[Dict[str, Any]] = []
    xtb_logger = logging.getLogger("src.xtb_backend")
    xtb_logger.setLevel(logging.ERROR)

    for target_id in args.targets:
        job = manifest.get(str(target_id))
        if job is None:
            raise KeyError(f"Target not found in manifest: {target_id}")
        result = _evaluate_target(job)
        results.append(result)
        barrier_1d = result["scan_1d"].get("barrier_1d_kcal_mol")
        barrier_2d = result["scan_2d"].get("barrier_2d_kcal_mol")
        print(
            f"{target_id}: pair1={result['pair1']['pair']} pair2="
            f"{None if result['pair2'] is None else result['pair2']['pair']} "
            f"barrier_1d={None if barrier_1d is None else round(barrier_1d, 2)} "
            f"barrier_2d={None if barrier_2d is None else round(barrier_2d, 2)} "
            f"supports_gsm={result['supports_gsm']}"
        )

    support_count = sum(1 for item in results if item.get("supports_gsm"))
    recommendation = "GSM_GO" if support_count >= 2 else "GSM_NO_GO"
    payload = {
        "targets": args.targets,
        "support_count": support_count,
        "recommendation": recommendation,
        "results": results,
    }
    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"ARTIFACT {output_path}")
    print(recommendation)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())