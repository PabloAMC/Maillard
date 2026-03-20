from __future__ import annotations

import io
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping

import numpy as np

from src.geometry_benchmark import GeometryBenchmarkEntry, load_geometry_benchmark_entries
from src.mlp_adoption_contract import load_mlp_candidates


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


def _assess_mace_geometry_candidate(entries: List[GeometryBenchmarkEntry], model_name: str) -> Dict[str, Any]:
    try:
        from src.mlp_optimizer import MLPOptimizer
    except ImportError:
        return {
            "available": False,
            "backend_available": False,
            "reason": "mace_backend_unavailable",
        }

    optimizer = MLPOptimizer(model_name=model_name, device="cpu")
    ground_state_entries = [entry for entry in entries if entry.benchmark_kind == "ground_state"]
    rmsd_rows: List[Dict[str, Any]] = []
    for entry in ground_state_entries:
        optimized_xyz = optimizer.optimize_geometry(entry.xyz, fmax=0.05, max_steps=100, drift_threshold=2.0)
        rmsd = _calculate_rmsd(entry.xyz, optimized_xyz)
        rmsd_rows.append(
            {
                "benchmark_id": entry.benchmark_id,
                "chemistry_family": entry.chemistry_family,
                "rmsd_angstrom": float(rmsd),
            }
        )

    if not rmsd_rows:
        return {
            "available": False,
            "backend_available": True,
            "reason": "no_ground_state_entries",
        }

    rmsd_values = [float(row["rmsd_angstrom"]) for row in rmsd_rows]
    return {
        "available": True,
        "backend_available": True,
        "reason": "evaluated",
        "benchmark_entry_count": len(ground_state_entries),
        "evaluated_entry_count": len(rmsd_rows),
        "max_rmsd_angstrom": max(rmsd_values),
        "mean_rmsd_angstrom": sum(rmsd_values) / len(rmsd_values),
        "entry_results": rmsd_rows,
    }


def build_geometry_assessment_artifact() -> Dict[str, Any]:
    entries = load_geometry_benchmark_entries()
    candidates = load_mlp_candidates()
    candidate_rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        if candidate.proposed_role == "geom_preopt" and candidate.model_family == "mace_mp":
            row = _assess_mace_geometry_candidate(entries, candidate.model_name)
        else:
            row = {
                "available": False,
                "backend_available": False,
                "reason": "no_candidate_backend_adapter",
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
            "benchmark_id": "p4_geometry_benchmark_v1",
            "candidate_count": len(candidate_rows),
            "ground_state_entry_count": sum(1 for entry in entries if entry.benchmark_kind == "ground_state"),
            "ts_seed_entry_count": sum(1 for entry in entries if entry.benchmark_kind == "ts_seed"),
        },
        "candidate_assessments": candidate_rows,
    }


def render_geometry_assessment_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# P4 Geometry Assessment",
        "",
        "| Candidate | Role | Available | Backend | Max RMSD (Å) | Mean RMSD (Å) | Reason |",
        "| --- | --- | --- | --- | ---: | ---: | --- |",
    ]
    for row in payload.get("candidate_assessments", []):
        max_rmsd = row.get("max_rmsd_angstrom")
        mean_rmsd = row.get("mean_rmsd_angstrom")
        lines.append(
            f"| {row.get('candidate_id', 'unknown')} | {row.get('proposed_role', 'unknown')} | "
            f"{'yes' if row.get('available', False) else 'no'} | {'yes' if row.get('backend_available', False) else 'no'} | "
            f"{'' if max_rmsd is None else f'{float(max_rmsd):.3f}'} | {'' if mean_rmsd is None else f'{float(mean_rmsd):.3f}'} | {row.get('reason', 'unknown')} |"
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