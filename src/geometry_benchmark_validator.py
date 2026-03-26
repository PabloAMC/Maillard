from __future__ import annotations

import io
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping

import numpy as np

from src.geometry_benchmark import GeometryBenchmarkEntry, load_geometry_benchmark_entries
from src.mlp_backend_adapters import build_candidate_adapter
from src.mlp_adoption_contract import load_mlp_candidates
from src.mlp_optimizer import compute_geometry_drift_metrics


def _calculate_rmsd(reference_xyz: str, test_xyz: str) -> float:
    try:
        from ase.build import minimize_rotation_and_translation
        from ase.io import read
    except ImportError:
        raise RuntimeError("ASE is required to compute RMSD for the geometry benchmark")

    with io.StringIO(reference_xyz.strip()) as ref_handle:
        ref_atoms = read(ref_handle, format="xyz")
    with io.StringIO(test_xyz.strip()) as test_handle:
        test_atoms = read(test_handle, format="xyz")

    if isinstance(ref_atoms, list):
        ref_atoms = ref_atoms[-1]
    if isinstance(test_atoms, list):
        test_atoms = test_atoms[-1]
    if len(ref_atoms) != len(test_atoms):
        raise ValueError("Geometry benchmark structures must preserve atom count")

    minimize_rotation_and_translation(ref_atoms, test_atoms)
    delta = ref_atoms.get_positions() - test_atoms.get_positions()
    return float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))


def _assess_geometry_candidate(entries: List[GeometryBenchmarkEntry], candidate) -> Dict[str, Any]:
    adapter = build_candidate_adapter(candidate)
    availability = adapter.probe_availability()
    if not availability.backend_available:
        return {
            "available": False,
            "backend_available": False,
            "reason": availability.reason,
        }
    if not availability.available:
        return {
            "available": False,
            "backend_available": True,
            "reason": availability.reason,
        }

    ground_state_entries = [entry for entry in entries if entry.benchmark_kind == "ground_state"]
    rmsd_rows: List[Dict[str, Any]] = []
    for entry in ground_state_entries:
        optimized_xyz = adapter.optimize_geometry(entry.xyz, fmax=0.05, max_steps=100)
        metrics = compute_geometry_drift_metrics(entry.xyz, optimized_xyz)
        rmsd_rows.append(
            {
                "benchmark_id": entry.benchmark_id,
                "chemistry_family": entry.chemistry_family,
                "rmsd_angstrom": float(metrics["rmsd_angstrom"] or 0.0),
                "max_atom_displacement_angstrom": float(metrics["max_atom_displacement_angstrom"] or 0.0),
                "hetero_atom_rmsd_angstrom": None if metrics["hetero_atom_rmsd_angstrom"] is None else float(metrics["hetero_atom_rmsd_angstrom"]),
                "sulfur_local_rmsd_angstrom": None if metrics["sulfur_local_rmsd_angstrom"] is None else float(metrics["sulfur_local_rmsd_angstrom"]),
                "sulfur_neighbor_max_delta_angstrom": None if metrics["sulfur_neighbor_max_delta_angstrom"] is None else float(metrics["sulfur_neighbor_max_delta_angstrom"]),
                "sulfur_bond_max_delta_angstrom": None if metrics["sulfur_bond_max_delta_angstrom"] is None else float(metrics["sulfur_bond_max_delta_angstrom"]),
                "sulfur_angle_max_delta_degrees": None if metrics["sulfur_angle_max_delta_degrees"] is None else float(metrics["sulfur_angle_max_delta_degrees"]),
            }
        )

    if not rmsd_rows:
        return {
            "available": False,
            "backend_available": True,
            "reason": "no_ground_state_entries",
        }

    rmsd_values = [float(row["rmsd_angstrom"]) for row in rmsd_rows]
    max_atom_drifts = [float(row["max_atom_displacement_angstrom"]) for row in rmsd_rows]
    hetero_rmsd_values = [float(row["hetero_atom_rmsd_angstrom"]) for row in rmsd_rows if row["hetero_atom_rmsd_angstrom"] is not None]
    sulfur_local_rmsd_values = [float(row["sulfur_local_rmsd_angstrom"]) for row in rmsd_rows if row["sulfur_local_rmsd_angstrom"] is not None]
    sulfur_neighbor_delta_values = [float(row["sulfur_neighbor_max_delta_angstrom"]) for row in rmsd_rows if row["sulfur_neighbor_max_delta_angstrom"] is not None]
    sulfur_bond_delta_values = [float(row["sulfur_bond_max_delta_angstrom"]) for row in rmsd_rows if row["sulfur_bond_max_delta_angstrom"] is not None]
    sulfur_angle_delta_values = [float(row["sulfur_angle_max_delta_degrees"]) for row in rmsd_rows if row["sulfur_angle_max_delta_degrees"] is not None]
    return {
        "available": True,
        "backend_available": True,
        "reason": "evaluated",
        "benchmark_entry_count": len(ground_state_entries),
        "evaluated_entry_count": len(rmsd_rows),
        "max_rmsd_angstrom": max(rmsd_values),
        "mean_rmsd_angstrom": sum(rmsd_values) / len(rmsd_values),
        "max_atom_displacement_angstrom": max(max_atom_drifts),
        "max_hetero_atom_rmsd_angstrom": max(hetero_rmsd_values) if hetero_rmsd_values else None,
        "max_sulfur_local_rmsd_angstrom": max(sulfur_local_rmsd_values) if sulfur_local_rmsd_values else None,
        "max_sulfur_neighbor_delta_angstrom": max(sulfur_neighbor_delta_values) if sulfur_neighbor_delta_values else None,
        "max_sulfur_bond_delta_angstrom": max(sulfur_bond_delta_values) if sulfur_bond_delta_values else None,
        "max_sulfur_angle_delta_degrees": max(sulfur_angle_delta_values) if sulfur_angle_delta_values else None,
        "entry_results": rmsd_rows,
    }


def build_geometry_assessment_artifact() -> Dict[str, Any]:
    entries = load_geometry_benchmark_entries()
    candidates = load_mlp_candidates()
    candidate_rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        if candidate.proposed_role == "geom_preopt":
            row = _assess_geometry_candidate(entries, candidate)
        else:
            row = {
                "available": False,
                "backend_available": False,
                "reason": "role_not_in_geometry_lane",
            }
        row.update(
            {
                "candidate_id": candidate.candidate_id,
                "model_family": candidate.model_family,
                "model_name": candidate.model_name,
                "proposed_role": candidate.proposed_role,
            }
        )
        candidate_rows.append(row)

    return {
        "summary": {
            "benchmark_id": "mlp_geometry_benchmark_v1",
            "candidate_count": len(candidate_rows),
            "ground_state_entry_count": sum(1 for entry in entries if entry.benchmark_kind == "ground_state"),
            "ts_seed_entry_count": sum(1 for entry in entries if entry.benchmark_kind == "ts_seed"),
        },
        "candidate_assessments": candidate_rows,
    }


def render_geometry_assessment_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Geometry Preoptimization Assessment",
        "",
        "| Candidate | Role | Available | Backend | Max RMSD (Å) | Max Atom Drift (Å) | Max Sulfur Bond Delta (Å) | Max Sulfur Angle Delta (deg) | Reason |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in payload.get("candidate_assessments", []):
        max_rmsd = row.get("max_rmsd_angstrom")
        max_atom_drift = row.get("max_atom_displacement_angstrom")
        sulfur_bond_delta = row.get("max_sulfur_bond_delta_angstrom")
        sulfur_angle_delta = row.get("max_sulfur_angle_delta_degrees")
        lines.append(
            f"| {row.get('candidate_id', 'unknown')} | {row.get('proposed_role', 'unknown')} | "
            f"{'yes' if row.get('available', False) else 'no'} | {'yes' if row.get('backend_available', False) else 'no'} | "
            f"{'' if max_rmsd is None else f'{float(max_rmsd):.3f}'} | {'' if max_atom_drift is None else f'{float(max_atom_drift):.3f}'} | "
            f"{'' if sulfur_bond_delta is None else f'{float(sulfur_bond_delta):.3f}'} | {'' if sulfur_angle_delta is None else f'{float(sulfur_angle_delta):.2f}'} | {row.get('reason', 'unknown')} |"
        )
    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark id: {summary.get('benchmark_id', 'unknown')}",
            f"Candidates assessed: {int(summary.get('candidate_count', 0))}",
            f"Ground-state entries: {int(summary.get('ground_state_entry_count', 0))}",
            f"TS-seed entries: {int(summary.get('ts_seed_entry_count', 0))}",
        ]
    )
    return "\n".join(lines) + "\n"