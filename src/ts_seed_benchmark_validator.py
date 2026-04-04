from __future__ import annotations

import io
from typing import Any, Dict, List, Mapping

import numpy as np

from src.mlp_adoption_contract import load_mlp_candidates
from src.mlp_backend_adapters import build_candidate_adapter
from src.ts_seed_benchmark import TSSeedBenchmarkEntry, load_ts_seed_benchmark_entries


def _calculate_rmsd(reference_xyz: str, test_xyz: str) -> float:
    try:
        from ase.build import minimize_rotation_and_translation
        from ase.io import read
    except ImportError:
        raise RuntimeError("ASE is required to compute RMSD for the TS-seed benchmark")

    with io.StringIO(reference_xyz.strip()) as ref_handle:
        ref_atoms = read(ref_handle, format="xyz")
    with io.StringIO(test_xyz.strip()) as test_handle:
        test_atoms = read(test_handle, format="xyz")

    if isinstance(ref_atoms, list):
        ref_atoms = ref_atoms[-1]
    if isinstance(test_atoms, list):
        test_atoms = test_atoms[-1]
    if len(ref_atoms) != len(test_atoms):
        raise ValueError("TS-seed benchmark structures must preserve atom count")

    minimize_rotation_and_translation(ref_atoms, test_atoms)
    delta = ref_atoms.get_positions() - test_atoms.get_positions()
    return float(np.sqrt(np.mean(np.sum(delta**2, axis=1))))


def _assess_ts_seed_candidate(entries: List[TSSeedBenchmarkEntry], candidate) -> Dict[str, Any]:
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

    rmsd_rows: List[Dict[str, Any]] = []
    for entry in entries:
        recovered_xyz = adapter.prepare_ts_seed(entry.challenged_seed_xyz, fmax=0.05, max_steps=100)
        rmsd = _calculate_rmsd(entry.reference_xyz, recovered_xyz)
        rmsd_rows.append(
            {
                "benchmark_id": entry.benchmark_id,
                "chemistry_family": entry.chemistry_family,
                "seed_rmsd_angstrom": float(rmsd),
            }
        )

    if not rmsd_rows:
        return {
            "available": False,
            "backend_available": True,
            "reason": "no_ts_seed_entries",
        }

    rmsd_values = [float(row["seed_rmsd_angstrom"]) for row in rmsd_rows]
    return {
        "available": True,
        "backend_available": True,
        "reason": "evaluated",
        "benchmark_entry_count": len(entries),
        "evaluated_entry_count": len(rmsd_rows),
        "max_seed_rmsd_angstrom": max(rmsd_values),
        "mean_seed_rmsd_angstrom": sum(rmsd_values) / len(rmsd_values),
        "entry_results": rmsd_rows,
    }


def build_ts_seed_assessment_artifact() -> Dict[str, Any]:
    entries = load_ts_seed_benchmark_entries()
    candidates = load_mlp_candidates()
    candidate_rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        if candidate.proposed_role == "ts_initialization":
            row = _assess_ts_seed_candidate(entries, candidate)
        else:
            row = {
                "available": False,
                "backend_available": False,
                "reason": "role_not_in_ts_seed_lane",
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
            "benchmark_id": "mlp_ts_seed_benchmark_v1",
            "candidate_count": len(candidate_rows),
            "entry_count": len(entries),
        },
        "candidate_assessments": candidate_rows,
    }


def render_ts_seed_assessment_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# TS Seed Recovery Assessment",
        "",
        "| Candidate | Role | Available | Backend | Max Seed RMSD (A) | Mean Seed RMSD (A) | Reason |",
        "| --- | --- | --- | --- | ---: | ---: | --- |",
    ]
    for row in payload.get("candidate_assessments", []):
        max_rmsd = row.get("max_seed_rmsd_angstrom")
        mean_rmsd = row.get("mean_seed_rmsd_angstrom")
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
            f"Entries: {int(summary.get('entry_count', 0))}",
        ]
    )
    return "\n".join(lines) + "\n"