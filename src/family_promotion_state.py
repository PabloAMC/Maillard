from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, Iterable, Mapping, Optional

from src.benchmark_validation import evaluate_benchmark, get_benchmark_files, load_benchmark, summarize_benchmarks
from src.family_ingestion_plan import load_family_ingestion_plan
from src.family_validation_overview import build_family_validation_overview_artifact
from src.literature_family_registry import build_family_payload_coverage_artifact, resolve_family_descriptor
from src.matrix_targets import get_compound_panel_entry


def _normalize_name(name: str) -> str:
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())


def _benchmark_compound_names(bench: Mapping[str, Any], summary: Any) -> list[str]:
    names: list[str] = []
    measured = bench.get("measured_volatiles", {}) or {}
    if isinstance(measured, dict):
        names.extend(str(name) for name in measured.keys())
    reference = bench.get("reference_volatiles", {}) or {}
    if isinstance(reference, dict):
        names.extend(str(name) for name in reference.keys())
    names.extend(str(name) for name in getattr(summary, "ranked_observable_targets", []))
    names.extend(str(name) for name in getattr(summary, "adverse_markers", []))
    deduped: list[str] = []
    seen: set[str] = set()
    for name in names:
        normalized = _normalize_name(name)
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(str(name))
    return deduped


def _family_support_benchmark_rows(
    benchmark_files: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    summaries = summarize_benchmarks(benchmark_files)
    rows: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
        "support_benchmark_ids": set(),
        "support_quantitative_benchmark_ids": set(),
        "support_quantitative_point_count": 0,
        "explicit_benchmark_uncertainty": False,
    })

    for summary in summaries:
        bench = load_benchmark(summary.bench_file)
        measured = bench.get("measured_volatiles", {}) or {}
        support_families_seen: set[str] = set()
        for compound_name in _benchmark_compound_names(bench, summary):
            panel_entry = get_compound_panel_entry(compound_name) or {}
            for family_ref in panel_entry.get("supporting_families", []) or []:
                descriptor = resolve_family_descriptor(str(family_ref))
                family_id = str(descriptor.get("family_id", "")).strip()
                if family_id:
                    support_families_seen.add(family_id)
        for family_id in support_families_seen:
            rows[family_id]["support_benchmark_ids"].add(summary.benchmark_id)

        if not summary.supported or summary.matched_compounds <= 0:
            continue
        evaluation = evaluate_benchmark(summary.bench_file)
        for comparison in evaluation.comparisons:
            if comparison.matched_name is None:
                continue
            if comparison.measured_ppb <= 0.0 or comparison.predicted_ppb <= 0.0:
                continue
            panel_entry = get_compound_panel_entry(comparison.compound) or get_compound_panel_entry(comparison.matched_name) or {}
            supporting_families = panel_entry.get("supporting_families", []) or []
            if not supporting_families:
                continue
            uncertainty = (measured.get(comparison.compound) or measured.get(comparison.matched_name) or {}).get("uncertainty_pct")
            for family_ref in supporting_families:
                descriptor = resolve_family_descriptor(str(family_ref))
                family_id = str(descriptor.get("family_id", "")).strip()
                if not family_id:
                    continue
                rows[family_id]["support_quantitative_benchmark_ids"].add(summary.benchmark_id)
                rows[family_id]["support_quantitative_point_count"] += 1
                if uncertainty is not None:
                    rows[family_id]["explicit_benchmark_uncertainty"] = True

    rendered: Dict[str, Dict[str, Any]] = {}
    for family_id, row in rows.items():
        rendered[family_id] = {
            "support_benchmark_count": len(row["support_benchmark_ids"]),
            "support_quantitative_benchmark_count": len(row["support_quantitative_benchmark_ids"]),
            "support_quantitative_point_count": int(row["support_quantitative_point_count"]),
            "support_benchmark_ids": sorted(row["support_benchmark_ids"]),
            "explicit_benchmark_uncertainty": bool(row["explicit_benchmark_uncertainty"]),
        }
    return rendered


def _promotion_state(row: Mapping[str, Any]) -> str:
    if bool(row.get("has_quantitative_parity", False)) and int(row.get("benchmark_count", 0)) >= 1:
        return "near_quantitative"
    if (
        int(row.get("support_benchmark_count", 0)) >= 1
        and bool(row.get("explicit_uncertainty_bounds", False))
        and "benchmark_ready" in str(row.get("benchmarkability_status", ""))
    ):
        return "benchmark_linked_support"
    if int(row.get("primary_payload_count", 0)) >= 1 and bool(row.get("explicit_uncertainty_bounds", False)):
        return "calibration_grade"
    if int(row.get("supporting_payload_count", 0)) >= 1 or int(row.get("primary_payload_count", 0)) >= 1:
        return "directional_only"
    return "structural_gap"


def _promotion_rationale(row: Mapping[str, Any]) -> str:
    state = str(row.get("promotion_state", "unknown"))
    if state == "near_quantitative":
        return "compound-level quantitative parity exists on primary benchmark-linked observables"
    if state == "benchmark_linked_support":
        return "benchmark-linked observables directly constrain this family through supporting-family compounds with explicit uncertainty"
    if state == "calibration_grade":
        return "machine-readable runtime payloads exist with explicit bounds, but direct benchmark-linked closure is still missing"
    if state == "directional_only":
        return "runtime support exists, but it remains transferred or gap-limited rather than benchmark-linked"
    return "no machine-readable runtime or benchmark-linked support yet"


def build_family_promotion_state_artifact(
    benchmark_files: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    family_plan = {
        str(row.get("family_id", "unknown")): row
        for row in load_family_ingestion_plan().get("families", [])
    }
    overview = build_family_validation_overview_artifact(benchmark_files)
    coverage = build_family_payload_coverage_artifact()
    support_rows = _family_support_benchmark_rows(benchmark_files or get_benchmark_files())
    coverage_rows = {
        str(row.get("family_id", "unknown")): row
        for row in coverage.get("families", [])
    }
    families = []
    for row in overview.get("families", []):
        family_id = str(row.get("chemistry_family", "unknown"))
        coverage_row = coverage_rows.get(family_id, {})
        support_row = support_rows.get(family_id, {})
        plan_row = family_plan.get(family_id, {})
        explicit_uncertainty_bounds = bool(
            support_row.get("explicit_benchmark_uncertainty", False)
            or int(coverage_row.get("total_primary_payload_count", 0)) >= 1
        )
        merged = {
            **row,
            "benchmarkability_status": str(plan_row.get("benchmarkability_status", "unknown")),
            "strategic_posture": str(plan_row.get("strategic_posture", row.get("strategic_posture", "unknown"))),
            "primary_payload_count": int(coverage_row.get("total_primary_payload_count", 0)),
            "supporting_payload_count": int(coverage_row.get("total_supporting_payload_count", 0)),
            "support_benchmark_count": int(support_row.get("support_benchmark_count", 0)),
            "support_quantitative_benchmark_count": int(support_row.get("support_quantitative_benchmark_count", 0)),
            "support_quantitative_point_count": int(support_row.get("support_quantitative_point_count", 0)),
            "support_benchmark_ids": list(support_row.get("support_benchmark_ids", [])),
            "explicit_uncertainty_bounds": explicit_uncertainty_bounds,
        }
        merged["promotion_state"] = _promotion_state(merged)
        merged["promotion_rationale"] = _promotion_rationale(merged)
        families.append(merged)

    families.sort(key=lambda item: (str(item.get("slr_family", "99")), str(item.get("chemistry_family", "unknown"))))
    return {
        "schema_version": "1.0",
        "description": "Family-level promotion state defining whether each chemistry family is near-quantitative, benchmark-linked support, calibration-grade, directional-only, or still a structural gap.",
        "promotion_rule": {
            "near_quantitative": "primary benchmark-linked family with compound-level quantitative parity",
            "benchmark_linked_support": "non-primary family constrained by benchmark-linked supporting compounds with explicit uncertainty",
            "calibration_grade": "machine-readable primary payload support exists with explicit bounds but no benchmark-linked support yet",
            "directional_only": "runtime support exists but remains transferred or gap-limited",
            "structural_gap": "no machine-readable support or benchmark-linked closure exists",
        },
        "families": families,
        "summary": {
            "family_count": len(families),
            "state_counts": dict(sorted({state: sum(1 for row in families if row.get("promotion_state") == state) for state in {str(row.get("promotion_state", "unknown")) for row in families}}.items())),
            "promoted_non_quantitative_families": [
                row["chemistry_family"] for row in families
                if row.get("promotion_state") in {"benchmark_linked_support", "calibration_grade"}
                and not bool(row.get("has_quantitative_parity", False))
            ],
        },
    }


def render_family_promotion_state_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Promotion State",
        "",
        "| SLR | Family | Promotion State | Benchmarks | Support Benchmarks | Quantitative Points | Support Quant Points | Primary Payloads | Explicit Bounds | Rationale |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('slr_family', '99')} | {row.get('chemistry_family', 'unknown')} | {row.get('promotion_state', 'unknown')} | {int(row.get('benchmark_count', 0))} | {int(row.get('support_benchmark_count', 0))} | {int(row.get('quantitative_point_count', 0))} | {int(row.get('support_quantitative_point_count', 0))} | {int(row.get('primary_payload_count', 0))} | {'yes' if row.get('explicit_uncertainty_bounds', False) else 'no'} | {row.get('promotion_rationale', 'unknown')} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Families summarized: {int(summary.get('family_count', 0))}",
        f"Promoted non-quantitative families: {', '.join(str(item) for item in summary.get('promoted_non_quantitative_families', [])) or 'none'}",
    ])
    return "\n".join(lines) + "\n"