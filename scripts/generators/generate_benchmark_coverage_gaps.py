#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_benchmark_evidence_audit, get_benchmark_files, get_benchmark_metadata, load_benchmark

INTAKE_REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"


TARGET_PROTEINS = ["free", "pea_iso", "soy_iso", "pea_conc", "soy_conc", "myco"]
TARGET_PROCESS_STATES = ["ambient_slurry", "aqueous_pre_extrusion_model", "heated_matrix"]


def _build_rows() -> list[dict[str, object]]:
    bench_files = get_benchmark_files()
    benches = [load_benchmark(path) for path in bench_files]
    evidence = build_matrix_benchmark_evidence_audit(bench_files)

    rows: list[dict[str, object]] = []

    for protein in TARGET_PROTEINS:
        matching = [bench for bench in benches if str(bench.get("protein_type", "free")) == protein]
        rows.append({
            "dimension": "protein_system",
            "category": protein,
            "benchmark_count": len(matching),
            "status": "covered" if matching else "gap",
            "note": "present in current benchmark set" if matching else "no benchmark currently covers this protein system",
        })

    observed_states = {}
    for bench in benches:
        state = str((bench.get("process_metadata") or {}).get("state", "unspecified"))
        observed_states[state] = observed_states.get(state, 0) + 1
    for state in TARGET_PROCESS_STATES:
        count = observed_states.get(state, 0)
        rows.append({
            "dimension": "process_state",
            "category": state,
            "benchmark_count": count,
            "status": "covered" if count else "gap",
            "note": "present in current benchmark set" if count else "no benchmark currently covers this process state",
        })

    free_precursor = sum(1 for bench in benches if get_benchmark_metadata(bench).execution_path == "free_precursor")
    matrix_only = sum(1 for bench in benches if get_benchmark_metadata(bench).execution_path == "matrix_only")
    matrix_precursor = sum(1 for bench in benches if get_benchmark_metadata(bench).execution_path == "matrix_precursor_augmented")
    rows.extend([
        {
            "dimension": "execution_path",
            "category": "free_precursor",
            "benchmark_count": free_precursor,
            "status": "covered" if free_precursor else "gap",
            "note": "free-AA sulfur and safety envelope",
        },
        {
            "dimension": "execution_path",
            "category": "matrix_only",
            "benchmark_count": matrix_only,
            "status": "covered" if matrix_only else "gap",
            "note": "matrix off-flavour intake/headspace envelope",
        },
        {
            "dimension": "execution_path",
            "category": "matrix_precursor_augmented",
            "benchmark_count": matrix_precursor,
            "status": "covered" if matrix_precursor else "gap",
            "note": "matrix meaty-positive reproducibility harness envelope",
        },
    ])

    external_matrix_meaty = sum(
        1
        for row in evidence
        if row.target_profile in {"meaty_positive", "mixed"}
        and row.external_data_status == "external_quantitative"
    )
    rows.append({
        "dimension": "scientific_gap",
        "category": "external_matrix_meaty_positive",
        "benchmark_count": external_matrix_meaty,
        "status": "covered" if external_matrix_meaty else "gap",
        "note": "wet-lab quantitative meaty-positive matrix benchmarks are the main blocker for broad alt-protein validation",
    })

    intake = json.loads(INTAKE_REGISTRY_PATH.read_text(encoding="utf-8")) if INTAKE_REGISTRY_PATH.exists() else {}
    for gap in intake.get("structural_gaps", []):
        gap_id = str(gap.get("id", "unknown"))
        if gap_id not in {
            "ppi_meaty_positive_matrix_benchmark",
            "spi_meaty_positive_matrix_benchmark",
            "meaty_off_flavour_safety_tradeoff_panel",
        }:
            continue
        near_misses = ", ".join(
            str(item.get("entry_id", "unknown"))
            for item in gap.get("near_miss_candidates", [])
            if isinstance(item, dict)
        ) or "none"
        rows.append({
            "dimension": "structural_gap",
            "category": gap_id,
            "benchmark_count": 0,
            "status": "gap",
            "note": f"{str(gap.get('why', ''))} Near misses: {near_misses}.",
            "requires_primary_data": bool(gap.get("requires_primary_data", True)),
            "closure_outcome": str(gap.get("closure_outcome", "unknown")),
            "evidence_state": str(gap.get("evidence_state", "unknown")),
        })
    return rows


def _render_markdown(rows: list[dict[str, object]]) -> str:
    lines = [
        "# Benchmark Coverage Gaps",
        "",
        "| Dimension | Category | Benchmarks | Status | Note |",
        "| --- | --- | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['dimension']} | {row['category']} | {row['benchmark_count']} | {row['status']} | {row['note']} |"
        )
    lines.extend([
        "",
        f"Rows: {len(rows)}",
        f"Gaps: {sum(1 for row in rows if row['status'] == 'gap')}",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    rows = _build_rows()
    markdown = _render_markdown(rows)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    markdown_path = output_dir / "benchmark_coverage_gaps.md"
    json_path = output_dir / "benchmark_coverage_gaps.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
