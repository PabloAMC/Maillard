from __future__ import annotations

from typing import Any, Dict, List, Mapping

from src.matrix_family_coverage import build_matrix_family_coverage_artifact


_SCOPE_PRIORITY_RANK = {
    "active_matrix_priority": 4,
    "bounded_next_candidate": 3,
    "scope_gap_to_rank": 2,
    "hold_current_posture": 1,
    "maintain_reference_core": 0,
    "defer": 0,
}

_IMPORTANCE_RANK = {
    "critical": 3,
    "high": 2,
    "medium": 1,
    "low": 0,
}

_CLOSABILITY_RANK = {
    "high": 3,
    "medium": 2,
    "low": 1,
    "hold": 0,
}


def _closability_band(row: Mapping[str, Any]) -> str:
    expansion_status = str(row.get("expansion_status", "unknown"))
    evidence_surface = str(row.get("evidence_surface", "unknown"))

    if expansion_status == "promote_primary_benchmark":
        return "high"
    if expansion_status == "bounded_expansion_candidate":
        return "medium"
    if expansion_status == "blocked_on_family_specific_evidence":
        return "medium" if evidence_surface == "generic_fat_fraction_plus_generic_lipid_chemistry" else "low"
    if expansion_status == "blocked_on_runtime_prior_and_benchmark":
        return "low"
    return "hold"


def _priority_score(row: Mapping[str, Any]) -> int:
    scope_priority = str(row.get("scope_priority", "defer"))
    importance_tier = str(row.get("importance_tier", "low"))
    closability = _closability_band(row)
    return (
        100 * _SCOPE_PRIORITY_RANK.get(scope_priority, 0)
        + 10 * _IMPORTANCE_RANK.get(importance_tier, 0)
        + _CLOSABILITY_RANK.get(closability, 0)
    )


def _evidence_landing(row: Mapping[str, Any]) -> str:
    matrix_family = str(row.get("matrix_family", "unknown"))
    expansion_status = str(row.get("expansion_status", "unknown"))
    evidence_surface = str(row.get("evidence_surface", "unknown"))

    if expansion_status == "promote_primary_benchmark":
        return "external_benchmark_plus_calibration_review"
    if expansion_status == "bounded_expansion_candidate" and evidence_surface == "bounded_calibration_prior":
        return "bounded_calibration_prior"
    if matrix_family == "soy_hydrolysate":
        return "qualitative_intake_only"
    if matrix_family == "extrusion_heavy_systems":
        return "process_regime_intake_and_warnings"
    if matrix_family == "coconut_oil_co_matrix":
        return "family_specific_calibration_or_tradeoff_benchmark"
    if expansion_status == "blocked_on_runtime_prior_and_benchmark":
        return "runtime_prior_plus_first_benchmark_design"
    return "maintain_current_evidence_surface"


def build_matrix_family_priority_ranking_artifact() -> Dict[str, Any]:
    coverage = build_matrix_family_coverage_artifact()
    ranked_rows: List[Dict[str, Any]] = []
    for row in coverage.get("families", []):
        if str(row.get("expansion_status", "")) == "reference_core":
            continue
        closability = _closability_band(row)
        ranked_rows.append(
            {
                "matrix_family": row.get("matrix_family", "unknown"),
                "category": row.get("category", "unknown"),
                "importance_tier": row.get("importance_tier", "unknown"),
                "runtime_posture": row.get("runtime_posture", "unknown"),
                "evidence_surface": row.get("evidence_surface", "unknown"),
                "support_class": row.get("support_class", "unknown"),
                "expansion_status": row.get("expansion_status", "unknown"),
                "scope_priority": row.get("scope_priority", "unknown"),
                "closability_band": closability,
                "evidence_landing": _evidence_landing(row),
                "priority_score": _priority_score(row),
                "primary_blocker": row.get("primary_blocker", "none"),
                "next_best_action": row.get("next_best_action", "unknown"),
            }
        )

    ranked_rows.sort(
        key=lambda row: (
            -int(row.get("priority_score", 0)),
            str(row.get("matrix_family", "unknown")),
        )
    )
    for index, row in enumerate(ranked_rows, start=1):
        row["priority_rank"] = index

    return {
        "summary": {
            "family_count": len(ranked_rows),
            "active_matrix_priorities": [
                row["matrix_family"]
                for row in ranked_rows
                if row["scope_priority"] == "active_matrix_priority"
            ],
            "bounded_next_candidates": [
                row["matrix_family"]
                for row in ranked_rows
                if row["scope_priority"] == "bounded_next_candidate"
            ],
            "scope_gap_priorities": [
                row["matrix_family"]
                for row in ranked_rows
                if row["scope_priority"] == "scope_gap_to_rank"
            ],
            "policy": "rank_matrix_family_scope_by_product_impact_and_closability_without_overclaiming_broad_pbma_coverage",
        },
        "families": ranked_rows,
    }


def render_matrix_family_priority_ranking_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Family Priority Ranking",
        "",
        "| Rank | Matrix Family | Scope Priority | Impact | Closability | Evidence Landing | Runtime Posture | Primary Blocker | Next Best Action |",
        "| ---: | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {int(row.get('priority_rank', 0))} | {row.get('matrix_family', 'unknown')} | {row.get('scope_priority', 'unknown')} | "
            f"{row.get('importance_tier', 'unknown')} | {row.get('closability_band', 'unknown')} | {row.get('evidence_landing', 'unknown')} | {row.get('runtime_posture', 'unknown')} | "
            f"{row.get('primary_blocker', 'none')} | {row.get('next_best_action', 'unknown')} |"
        )

    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Families ranked: {int(summary.get('family_count', 0))}",
        f"Active matrix priorities: {', '.join(str(item) for item in summary.get('active_matrix_priorities', [])) or 'none'}",
        f"Bounded next candidates: {', '.join(str(item) for item in summary.get('bounded_next_candidates', [])) or 'none'}",
        f"Scope-gap priorities: {', '.join(str(item) for item in summary.get('scope_gap_priorities', [])) or 'none'}",
        f"Policy: {summary.get('policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"