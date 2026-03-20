#!/usr/bin/env python3
"""
generate_matrix_target_status.py — P0.4 artifact script.

Iterates all matrix benchmarks, computes MatrixBenchmarkEvidence, and emits
a JSON artifact at results/matrix_target_status.json.

For each compound in the ranking contract, classifies:
  - closed      : compound is externally or internally benchmarked
  - directional : class anchor exists (class_anchored or directional_transferred)
  - open        : only heuristic or missing calibration data

Also emits a per-benchmark promotion_ready flag based on the rule:
  ≥ 2 compounds are closed/directional AND ranking_contract_status != order_mismatch.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    evaluate_benchmark,
    summarize_evaluation,
    assess_matrix_benchmark_evidence,
    get_benchmark_files,
    get_matrix_ranking_contract,
    load_benchmark,
    get_benchmark_metadata,
)
from src.matrix_calibration_registry import describe_matrix_calibration
from src.benchmark_validation import determine_matrix_process_state


_CLOSED_STRENGTHS = {
    "externally_benchmarked",
    "internally_benchmarked",
    "literature_anchored",
    "conditional_literature_anchored",
    "compound_specific",
    "compound_specific_process_state",
    "process_state_mismatch",
}
_DIRECTIONAL_STRENGTHS = {
    "class_anchored",
    "directional_transferred",
}


def _resolve_compound_status(strength: str) -> str:
    if strength in _CLOSED_STRENGTHS:
        return "closed"
    if strength in _DIRECTIONAL_STRENGTHS:
        return "directional"
    return "open"


def _process_benchmark(bench_file: Path) -> dict | None:
    bench = load_benchmark(bench_file)
    metadata = get_benchmark_metadata(bench)
    if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
        return None

    protein_type = str(bench.get("protein_type", "free"))
    conditions = bench.get("conditions", {})
    process_state = determine_matrix_process_state(
        temperature_celsius=float(conditions.get("temp_C", 100.0)),
        time_minutes=float(conditions.get("time_min", 45.0)),
    )

    contract = get_matrix_ranking_contract(bench)
    observable_targets = contract.get("observable_targets", [])
    adverse_markers = contract.get("adverse_markers", [])
    all_compounds = [str(item.get("name", "")) for item in observable_targets] + [
        str(m) for m in adverse_markers
    ]

    compound_details = []
    closed_count = 0
    directional_count = 0
    for compound_name in all_compounds:
        if not compound_name:
            continue
        cal = describe_matrix_calibration(
            compound_name,
            protein_type=protein_type,
            process_state=process_state,
        )
        strength = str(cal.get("calibration_evidence_strength", "heuristic"))
        status = _resolve_compound_status(strength)
        if status == "closed":
            closed_count += 1
        elif status == "directional":
            directional_count += 1

        role = "adverse_marker" if compound_name in adverse_markers else "desirable_marker"
        compound_details.append(
            {
                "compound": compound_name,
                "role": role,
                "status": status,
                "calibration_evidence_strength": strength,
                "calibration_source": str(cal.get("calibration_source", "")),
            }
        )

    # P0.4 promotion rule: ≥2 anchored (closed or directional) AND no hard order_mismatch
    anchored_count = closed_count + directional_count
    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type=protein_type)
    promotion_ready = (
        anchored_count >= 2
        and summary.ranking_contract_status not in {"order_mismatch"}
    )

    return {
        "benchmark_id": str(bench.get("benchmark_id", bench_file.stem)),
        "bench_file": str(bench_file),
        "protein_type": protein_type,
        "process_state": process_state,
        "execution_path": metadata.execution_path,
        "ranking_contract_status": summary.ranking_contract_status,
        "compound_count": len(compound_details),
        "closed_count": closed_count,
        "directional_count": directional_count,
        "open_count": len(compound_details) - closed_count - directional_count,
        "promotion_ready": promotion_ready,
        "promotion_rule": ">=2 anchored (closed|directional) AND ranking_contract_status != order_mismatch",
        "compounds": compound_details,
    }


def main() -> int:
    output_dir = ROOT / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "matrix_target_status.json"

    bench_files = get_benchmark_files()
    rows = []
    for bench_file in bench_files:
        try:
            row = _process_benchmark(bench_file)
        except Exception as exc:
            print(f"[WARN] Skipping {bench_file.name}: {exc}", file=sys.stderr)
            continue
        if row is not None:
            rows.append(row)
            ready = "✓" if row["promotion_ready"] else "✗"
            print(
                f"  [{ready}] {row['benchmark_id']:60s}  "
                f"closed={row['closed_count']} directional={row['directional_count']} "
                f"open={row['open_count']}  ranking={row['ranking_contract_status']}"
            )

    artifact = {
        "schema_version": "1.0",
        "description": "Matrix target observability status — P0.4",
        "promotion_rule": ">=2 anchored AND ranking_contract_status != order_mismatch",
        "benchmarks": rows,
        "summary": {
            "total_benchmarks": len(rows),
            "promotion_ready": sum(1 for r in rows if r["promotion_ready"]),
        },
    }
    output_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nArtifact written to {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
