from __future__ import annotations

import math
import json
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional

from src.benchmark_validation import (
    DEFAULT_BENCHMARK_THRESHOLDS,
    DEFAULT_TARGET_TAG,
    evaluate_benchmark,
    summarize_benchmarks,
)
from src.family_ingestion_plan import load_family_ingestion_plan
from src.literature_family_registry import build_family_payload_coverage_artifact
from src.literature_family_registry import resolve_family_descriptor
from src.matrix_targets import get_compound_panel_entry
from src.validation_contract import BenchmarkThresholds


def _load_benchmark_payload(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _base_family_row(plan_row: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "slr_family": str(plan_row.get("slr_family", "99")).zfill(2),
        "chemistry_family": str(plan_row.get("family_id", "unknown")),
        "display_name": str(plan_row.get("display_name", plan_row.get("family_id", "unknown"))),
        "strategic_posture": str(plan_row.get("strategic_posture", "unknown")),
        "benchmark_count": 0,
        "strict_ready_count": 0,
        "supported_count": 0,
        "quantitative_point_count": 0,
        "benchmark_ids": set(),
        "quantitative_benchmark_ids": set(),
        "execution_paths": defaultdict(int),
        "compound_ratios": [],
        "compound_abs_log10_errors": [],
        "primary_payload_count": 0,
        "supporting_payload_count": 0,
        "has_runtime_support": False,
    }


def _ensure_family_row(rows: Dict[str, Dict[str, Any]], family_id: str) -> Dict[str, Any]:
    if family_id in rows:
        return rows[family_id]
    descriptor = resolve_family_descriptor(family_id)
    row = _base_family_row(
        {
            "slr_family": str(descriptor.get("slr_family", "99")).zfill(2),
            "family_id": family_id,
            "display_name": str(descriptor.get("display_name", family_id)),
            "strategic_posture": str(descriptor.get("strategic_posture", "unknown")),
        }
    )
    rows[family_id] = row
    return row


def build_family_validation_overview_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> Dict[str, Any]:
    plan = load_family_ingestion_plan()
    payload_coverage = build_family_payload_coverage_artifact()
    coverage_by_family = {
        str(row.get("family_id", "unknown")): dict(row)
        for row in payload_coverage.get("families", [])
    }
    family_rows: Dict[str, Dict[str, Any]] = {
        str(row.get("family_id", "unknown")): {
            **_base_family_row(dict(row)),
            "primary_payload_count": int(coverage_by_family.get(str(row.get("family_id", "unknown")), {}).get("total_primary_payload_count", 0)),
            "supporting_payload_count": int(coverage_by_family.get(str(row.get("family_id", "unknown")), {}).get("total_supporting_payload_count", 0)),
            "has_runtime_support": bool(coverage_by_family.get(str(row.get("family_id", "unknown")), {}).get("has_runtime_support", False)),
        }
        for row in plan.get("families", [])
        if str(row.get("family_id", "")).strip()
    }

    summaries = summarize_benchmarks(benchmark_files, target_tag=target_tag, thresholds=thresholds)
    quantitative_points: List[Dict[str, Any]] = []

    for summary in summaries:
        bench_payload = _load_benchmark_payload(summary.bench_file)
        source_origin = str((bench_payload.get("source_metadata", {}) or {}).get("origin", "unknown"))
        for family_id in summary.chemistry_families:
            row = _ensure_family_row(family_rows, family_id)
            row["benchmark_count"] += 1
            row["strict_ready_count"] += int(bool(summary.strict_ready))
            row["supported_count"] += int(bool(summary.supported))
            row["benchmark_ids"].add(summary.benchmark_id)
            row["execution_paths"][summary.execution_path] += 1

        if not summary.supported or summary.matched_compounds <= 0:
            continue

        evaluation = evaluate_benchmark(summary.bench_file, target_tag=target_tag)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            if comparison.measured_ppb <= 0.0 or comparison.predicted_ppb <= 0.0:
                continue

            panel_entry = get_compound_panel_entry(comparison.compound) or get_compound_panel_entry(comparison.matched_name)
            chemistry_family = str((panel_entry or {}).get("chemistry_family", "")).strip()
            if not chemistry_family:
                continue

            descriptor = resolve_family_descriptor(chemistry_family)
            row = _ensure_family_row(family_rows, chemistry_family)
            abs_log10_error = abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb))
            row["quantitative_point_count"] += 1
            row["quantitative_benchmark_ids"].add(summary.benchmark_id)
            row["compound_ratios"].append(float(comparison.ratio))
            row["compound_abs_log10_errors"].append(abs_log10_error)

            quantitative_points.append(
                {
                    "benchmark_id": summary.benchmark_id,
                    "compound": comparison.compound,
                    "matched_name": comparison.matched_name,
                    "chemistry_family": chemistry_family,
                    "slr_family": str(descriptor.get("slr_family", row["slr_family"])).zfill(2),
                    "display_name": str(descriptor.get("display_name", row["display_name"])),
                    "strategic_posture": str(descriptor.get("strategic_posture", row["strategic_posture"])),
                    "execution_path": summary.execution_path,
                    "strict_ready": bool(summary.strict_ready),
                    "reference_signal_origin": str(summary.reference_signal_origin),
                    "source_origin": source_origin,
                    "measured_ppb": float(comparison.measured_ppb),
                    "predicted_ppb": float(comparison.predicted_ppb),
                    "uncertainty_pct": comparison.uncertainty_pct,
                    "compound_ratio": float(comparison.ratio),
                    "abs_log10_error": abs_log10_error,
                }
            )

    rendered_rows: List[Dict[str, Any]] = []
    for row in sorted(family_rows.values(), key=lambda item: (item["slr_family"], item["chemistry_family"])):
        rendered_rows.append(
            {
                "slr_family": row["slr_family"],
                "chemistry_family": row["chemistry_family"],
                "display_name": row["display_name"],
                "strategic_posture": row["strategic_posture"],
                "benchmark_count": int(row["benchmark_count"]),
                "strict_ready_count": int(row["strict_ready_count"]),
                "supported_count": int(row["supported_count"]),
                "quantitative_point_count": int(row["quantitative_point_count"]),
                "quantitative_benchmark_count": len(row["quantitative_benchmark_ids"]),
                "primary_payload_count": int(row["primary_payload_count"]),
                "supporting_payload_count": int(row["supporting_payload_count"]),
                "has_runtime_support": bool(row["has_runtime_support"]),
                "benchmark_ids": sorted(row["benchmark_ids"]),
                "execution_paths": dict(sorted(row["execution_paths"].items())),
                "median_compound_ratio": median(row["compound_ratios"]) if row["compound_ratios"] else None,
                "max_compound_ratio": max(row["compound_ratios"]) if row["compound_ratios"] else None,
                "mean_abs_log10_error": (
                    sum(row["compound_abs_log10_errors"]) / len(row["compound_abs_log10_errors"])
                    if row["compound_abs_log10_errors"]
                    else None
                ),
                "has_quantitative_parity": bool(row["compound_ratios"]),
                "integration_status": (
                    "quantitative_parity"
                    if row["compound_ratios"]
                    else "benchmark_linked"
                    if row["benchmark_count"] > 0
                    else "runtime_integrated"
                    if row["has_runtime_support"]
                    else "not_integrated"
                ),
            }
        )

    quantitative_families = [row for row in rendered_rows if row["has_quantitative_parity"]]
    benchmark_backed_families = [row for row in rendered_rows if row["benchmark_count"] > 0]
    integrated_families = [row for row in rendered_rows if row["has_runtime_support"]]

    return {
        "summary": {
            "family_count": len(rendered_rows),
            "integrated_family_count": len(integrated_families),
            "benchmark_backed_family_count": len(benchmark_backed_families),
            "quantitative_family_count": len(quantitative_families),
            "quantitative_point_count": len(quantitative_points),
            "families_without_quantitative_parity": [
                row["chemistry_family"] for row in rendered_rows if not row["has_quantitative_parity"]
            ],
            "integrated_families_without_quantitative_parity": [
                row["chemistry_family"] for row in integrated_families if not row["has_quantitative_parity"]
            ],
        },
        "families": rendered_rows,
        "integrated_families": integrated_families,
        "quantitative_families": quantitative_families,
        "quantitative_points": quantitative_points,
    }


def render_family_validation_overview_markdown(payload: Dict[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_family_validation_figures.py". Manual edits will be overwritten. -->',
        "",
        "# Family Validation Overview",
        "",
        "This artifact answers the product question directly: which chemistry families already have experiment-linked predictive closure, and which remain calibration-only or directional lanes.",
        "",
        f"- Families tracked: {int(summary.get('family_count', 0))}",
        f"- Families with landed runtime integration: {int(summary.get('integrated_family_count', 0))}",
        f"- Families with benchmark-linked experimental support: {int(summary.get('benchmark_backed_family_count', 0))}",
        f"- Families with compound-level quantitative parity points: {int(summary.get('quantitative_family_count', 0))}",
        f"- Quantitative compound points plotted: {int(summary.get('quantitative_point_count', 0))}",
        "",
        "How to read the PNG:",
        "",
        "- single panel: measured vs predicted compound parity, colored by chemistry family and shaped by execution path",
        "- only families with executable numeric benchmarks and matched measured compounds can appear in this scatter",
        "- integrated support or upstream lanes without direct measured compounds remain visible in the table below even when they have zero parity points",
        "",
        "| SLR | Family | Posture | Integrated | Benchmarks | Strict Ready | Quantitative Points | Median Ratio | Mean |log10 error| |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in payload.get("families", []):
        median_ratio = row.get("median_compound_ratio")
        mean_abs_log10_error = row.get("mean_abs_log10_error")
        lines.append(
            f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {row.get('strategic_posture', 'unknown')} | {bool(row.get('has_runtime_support', False))} | {int(row.get('benchmark_count', 0))} | {int(row.get('strict_ready_count', 0))} | {int(row.get('quantitative_point_count', 0))} | {median_ratio:.3f} | {mean_abs_log10_error:.3f} |"
            if median_ratio is not None and mean_abs_log10_error is not None
            else f"| {row.get('slr_family', '99')} | {row.get('display_name', 'unknown')} | {row.get('strategic_posture', 'unknown')} | {bool(row.get('has_runtime_support', False))} | {int(row.get('benchmark_count', 0))} | {int(row.get('strict_ready_count', 0))} | {int(row.get('quantitative_point_count', 0))} | - | - |"
        )
    return "\n".join(lines) + "\n"