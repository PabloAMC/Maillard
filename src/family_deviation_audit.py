from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional

from src.benchmark_labels import benchmark_label
from src.benchmark_validation import DEFAULT_BENCHMARK_THRESHOLDS, DEFAULT_TARGET_TAG
from src.family_validation_overview import build_family_validation_overview_artifact
from src.validation_contract import BenchmarkThresholds


def _p90(values: List[float]) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    idx = int(round(0.9 * (len(ordered) - 1)))
    return float(ordered[idx])


def _trimmed_mean(values: List[float], trim_fraction: float = 0.1) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    trim = int(len(ordered) * trim_fraction)
    if trim == 0 or len(ordered) < 2 * trim + 1:
        return float(sum(ordered) / len(ordered))
    trimmed = ordered[trim:-trim]
    return float(sum(trimmed) / len(trimmed))


def _classify_root_cause(point: Dict[str, Any]) -> str:
    source_origin = str(point.get("source_origin", "unknown"))
    signal_origin = str(point.get("reference_signal_origin", "measured_volatiles"))
    ratio = float(point.get("compound_ratio", 0.0))
    abs_log10_error = float(point.get("abs_log10_error", 0.0))

    if source_origin == "internal_reproducibility_candidate" and signal_origin == "reference_volatiles":
        if ratio >= 5.0 or abs_log10_error >= 0.5:
            return "internal_synthetic_reference_drift"
        return "internal_synthetic_reference_match"
    if signal_origin == "reference_volatiles":
        return "non_experimental_reference_comparator"
    if ratio >= 5.0 or abs_log10_error >= 0.5:
        return "model_or_mapping_mismatch"
    return "within_expected_error"


def _recommended_action(root_cause: str) -> str:
    if root_cause == "internal_synthetic_reference_drift":
        return "refresh_internal_reference_payload_or_remove_from_external_accuracy_claims"
    if root_cause == "non_experimental_reference_comparator":
        return "keep_as_regression_lane_not_as_external_accuracy_evidence"
    if root_cause == "model_or_mapping_mismatch":
        return "triage_mapping_projection_calibration"
    return "no_action_required"


def build_family_deviation_audit_artifact(
    benchmark_files: Optional[Iterable[str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
    high_ratio_threshold: float = 5.0,
    high_log10_error_threshold: float = 0.5,
    top_n_outliers: int = 20,
) -> Dict[str, Any]:
    payload = build_family_validation_overview_artifact(
        benchmark_files=benchmark_files,
        target_tag=target_tag,
        thresholds=thresholds,
    )
    points = [dict(row) for row in payload.get("quantitative_points", [])]

    root_cause_counts: Dict[str, int] = defaultdict(int)
    for row in points:
        root_cause = _classify_root_cause(row)
        row["root_cause"] = root_cause
        row["recommended_action"] = _recommended_action(root_cause)
        root_cause_counts[root_cause] += 1

    experimental_points = [
        row
        for row in points
        if str(row.get("reference_signal_origin", "")) == "measured_volatiles"
        and str(row.get("source_origin", "")) != "internal_reproducibility_candidate"
    ]
    by_ratio = sorted(points, key=lambda row: float(row.get("compound_ratio", 0.0)), reverse=True)
    by_log10 = sorted(points, key=lambda row: float(row.get("abs_log10_error", 0.0)), reverse=True)
    experimental_by_ratio = sorted(experimental_points, key=lambda row: float(row.get("compound_ratio", 0.0)), reverse=True)
    experimental_by_log10 = sorted(experimental_points, key=lambda row: float(row.get("abs_log10_error", 0.0)), reverse=True)
    experimental_high_ratio_count = sum(
        1 for row in experimental_points if float(row.get("compound_ratio", 0.0)) >= high_ratio_threshold
    )

    family_rows: List[Dict[str, Any]] = []
    for family in payload.get("families", []):
        family_id = str(family.get("chemistry_family", "unknown"))
        family_points = [row for row in points if str(row.get("chemistry_family", "")) == family_id]
        ratios = [float(row.get("compound_ratio", 0.0)) for row in family_points]
        abs_errors = [float(row.get("abs_log10_error", 0.0)) for row in family_points]
        high_ratio_count = sum(1 for value in ratios if value >= high_ratio_threshold)
        high_log_count = sum(1 for value in abs_errors if value >= high_log10_error_threshold)
        dominant = sorted(family_points, key=lambda row: float(row.get("compound_ratio", 0.0)), reverse=True)[:3]
        family_rows.append(
            {
                "slr_family": str(family.get("slr_family", "99")),
                "chemistry_family": family_id,
                "display_name": str(family.get("display_name", family_id)),
                "strategic_posture": str(family.get("strategic_posture", "unknown")),
                "quantitative_point_count": int(family.get("quantitative_point_count", 0)),
                "benchmark_count": int(family.get("benchmark_count", 0)),
                "median_ratio": family.get("median_compound_ratio"),
                "p90_ratio": _p90(ratios),
                "max_ratio": family.get("max_compound_ratio"),
                "mean_abs_log10_error": family.get("mean_abs_log10_error"),
                "trimmed_mean_abs_log10_error": _trimmed_mean(abs_errors),
                "high_ratio_count": high_ratio_count,
                "high_log10_error_count": high_log_count,
                "dominant_outliers": [
                    {
                        "benchmark_id": row.get("benchmark_id"),
                        "compound": row.get("compound"),
                        "execution_path": row.get("execution_path"),
                        "compound_ratio": row.get("compound_ratio"),
                        "abs_log10_error": row.get("abs_log10_error"),
                    }
                    for row in dominant
                ],
            }
        )

    ranked_families = sorted(
        family_rows,
        key=lambda row: (
            float(row.get("max_ratio") or 0.0),
            float(row.get("mean_abs_log10_error") or 0.0),
        ),
        reverse=True,
    )

    fix_queue_map: Dict[str, Dict[str, Any]] = {}
    for row in experimental_by_ratio:
        if float(row.get("compound_ratio", 0.0)) < high_ratio_threshold:
            continue
        ticket_key = "{benchmark}|{cause}".format(
            benchmark=row.get("benchmark_id", "unknown"),
            cause=row.get("root_cause", "unknown"),
        )
        bucket = fix_queue_map.setdefault(
            ticket_key,
            {
                "ticket_id": f"FD-{len(fix_queue_map) + 1:03d}",
                "benchmark_id": row.get("benchmark_id", "unknown"),
                "root_cause": row.get("root_cause", "unknown"),
                "recommended_action": row.get("recommended_action", "no_action_required"),
                "execution_path": row.get("execution_path", "unknown"),
                "compounds": [],
                "max_ratio_in_cluster": 0.0,
            },
        )
        bucket["compounds"].append(str(row.get("compound", "unknown")))
        bucket["max_ratio_in_cluster"] = max(
            float(bucket.get("max_ratio_in_cluster", 0.0)),
            float(row.get("compound_ratio", 0.0)),
        )

    return {
        "summary": {
            "quantitative_point_count": len(points),
            "family_count": len(payload.get("families", [])),
            "families_with_quantitative_points": sum(1 for row in family_rows if int(row.get("quantitative_point_count", 0)) > 0),
            "high_ratio_threshold": float(high_ratio_threshold),
            "high_log10_error_threshold": float(high_log10_error_threshold),
            "high_ratio_point_count": sum(1 for row in experimental_points if float(row.get("compound_ratio", 0.0)) >= high_ratio_threshold),
            "high_log10_error_point_count": sum(1 for row in experimental_points if float(row.get("abs_log10_error", 0.0)) >= high_log10_error_threshold),
            "root_cause_counts": dict(sorted(root_cause_counts.items())),
            "experimental_quantitative_point_count": len(experimental_points),
            "experimental_high_ratio_point_count": experimental_high_ratio_count,
            "max_observed_ratio": float(experimental_by_ratio[0].get("compound_ratio", 0.0)) if experimental_by_ratio else None,
            "worst_point": experimental_by_ratio[0] if experimental_by_ratio else None,
        },
        "family_deviation_table": ranked_families,
        "top_outliers_by_ratio": by_ratio[:top_n_outliers],
        "top_outliers_by_abs_log10_error": by_log10[:top_n_outliers],
        "experimental_top_outliers_by_ratio": experimental_by_ratio[:top_n_outliers],
        "experimental_top_outliers_by_abs_log10_error": experimental_by_log10[:top_n_outliers],
        "fix_queue": list(fix_queue_map.values()),
    }


def render_family_deviation_audit_markdown(payload: Dict[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Family Deviation Audit",
        "",
        "This artifact identifies which benchmark points dominate family-level error tails and where targeted closure work should start.",
        "",
        f"- Quantitative points analyzed: {int(summary.get('quantitative_point_count', 0))}",
        f"- Families tracked: {int(summary.get('family_count', 0))}",
        f"- Families with quantitative points: {int(summary.get('families_with_quantitative_points', 0))}",
        f"- High-ratio threshold: {float(summary.get('high_ratio_threshold', 0.0)):.2f}",
        f"- High-ratio points: {int(summary.get('high_ratio_point_count', 0))}",
        f"- High-log10-error threshold: {float(summary.get('high_log10_error_threshold', 0.0)):.2f}",
        f"- High-log10-error points: {int(summary.get('high_log10_error_point_count', 0))}",
        f"- Experimental points (measured-volatiles only): {int(summary.get('experimental_quantitative_point_count', 0))}",
        f"- Experimental high-ratio points: {int(summary.get('experimental_high_ratio_point_count', 0))}",
        "",
        "## Root-Cause Distribution",
        "",
    ]
    for cause, count in (summary.get("root_cause_counts", {}) or {}).items():
        lines.append(f"- {cause}: {int(count)}")

    lines.extend([
        "",
        "## Family Deviation Table",
        "",
        "| SLR | Family | Benchmarks | Quant Points | Median Ratio | P90 Ratio | Max Ratio | Mean |log10 error| | Trimmed Mean |log10 error| | High Ratio Pts |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in payload.get("family_deviation_table", []):
        median_ratio = row.get("median_ratio")
        p90_ratio = row.get("p90_ratio")
        max_ratio = row.get("max_ratio")
        mean_error = row.get("mean_abs_log10_error")
        trimmed_error = row.get("trimmed_mean_abs_log10_error")
        lines.append(
            "| {slr} | {name} | {bench} | {pts} | {med} | {p90} | {mx} | {mean} | {trim} | {high} |".format(
                slr=row.get("slr_family", "99"),
                name=row.get("display_name", "unknown"),
                bench=int(row.get("benchmark_count", 0)),
                pts=int(row.get("quantitative_point_count", 0)),
                med=f"{float(median_ratio):.3f}" if median_ratio is not None else "-",
                p90=f"{float(p90_ratio):.3f}" if p90_ratio is not None else "-",
                mx=f"{float(max_ratio):.3f}" if max_ratio is not None else "-",
                mean=f"{float(mean_error):.3f}" if mean_error is not None else "-",
                trim=f"{float(trimmed_error):.3f}" if trimmed_error is not None else "-",
                high=int(row.get("high_ratio_count", 0)),
            )
        )

    lines.extend([
        "",
        "## Fix Queue",
        "",
        "| Ticket | Benchmark | Root Cause | Action | Max Ratio | Compounds |",
        "| --- | --- | --- | --- | ---: | --- |",
    ])
    for row in payload.get("fix_queue", []):
        lines.append(
            f"| {row.get('ticket_id', 'FD-000')} | {benchmark_label(str(row.get('benchmark_id', 'unknown')))} | {row.get('root_cause', 'unknown')} | {row.get('recommended_action', 'none')} | {float(row.get('max_ratio_in_cluster', 0.0)):.3f} | {', '.join(sorted(set(row.get('compounds', []))))} |"
        )

    lines.extend([
        "",
        "## Top Outliers by Ratio",
        "",
        "| Benchmark | Family | Compound | Lane | Ratio | |log10 error| |",
        "| --- | --- | --- | --- | ---: | ---: |",
    ])
    for row in payload.get("top_outliers_by_ratio", [])[:20]:
        lines.append(
            f"| {benchmark_label(str(row.get('benchmark_id', 'unknown')))} | {row.get('display_name', row.get('chemistry_family', 'unknown'))} | {row.get('compound', 'unknown')} | {row.get('execution_path', 'unknown')} | {float(row.get('compound_ratio', 0.0)):.3f} | {float(row.get('abs_log10_error', 0.0)):.3f} |"
        )
    return "\n".join(lines) + "\n"