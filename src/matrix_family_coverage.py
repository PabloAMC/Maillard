from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src.artifact_io import load_json_mapping, repo_root


DEFAULT_MATRIX_FAMILY_COVERAGE_REGISTRY = repo_root() / "data" / "lit" / "matrix_family_coverage_registry.json"


def load_matrix_family_coverage_registry(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_MATRIX_FAMILY_COVERAGE_REGISTRY
    return load_json_mapping(path)


def _support_class(row: Mapping[str, Any]) -> str:
    runtime_posture = str(row.get("runtime_posture", "unknown"))
    if runtime_posture in {"quantitative_core", "directional_matrix", "qualitative_intake_only"}:
        return "explicit_supported"
    if runtime_posture == "indirect_generic_support":
        return "indirect_only"
    if runtime_posture == "open_gap":
        return "open_gap"
    return "unknown"


def _expansion_status(row: Mapping[str, Any]) -> str:
    matrix_family = str(row.get("matrix_family", "unknown"))
    runtime_posture = str(row.get("runtime_posture", "unknown"))
    evidence_surface = str(row.get("evidence_surface", "unknown"))
    if runtime_posture == "quantitative_core":
        return "reference_core"
    if matrix_family in {"pea_isolate", "soy_isolate"}:
        return "promote_primary_benchmark"
    if runtime_posture == "qualitative_intake_only":
        return "hold_intake_only"
    if matrix_family == "mycoprotein" and evidence_surface == "bounded_calibration_prior":
        return "bounded_expansion_candidate"
    if matrix_family == "extrusion_heavy_systems":
        return "hold_process_regime_only"
    if runtime_posture == "indirect_generic_support":
        return "blocked_on_family_specific_evidence"
    if runtime_posture == "open_gap":
        return "blocked_on_runtime_prior_and_benchmark"
    return "unknown"


def _primary_blocker(row: Mapping[str, Any]) -> str:
    unsupported = row.get("what_is_not_supported", [])
    if isinstance(unsupported, list) and unsupported:
        return str(unsupported[0])
    return "none"


def _scope_priority(row: Mapping[str, Any]) -> str:
    expansion_status = str(row.get("expansion_status", "unknown"))
    importance_tier = str(row.get("importance_tier", "unknown"))

    if expansion_status == "reference_core":
        return "maintain_reference_core"
    if expansion_status == "promote_primary_benchmark":
        return "active_matrix_priority"
    if expansion_status == "bounded_expansion_candidate":
        return "bounded_next_candidate"
    if expansion_status in {"blocked_on_family_specific_evidence", "blocked_on_runtime_prior_and_benchmark"} and importance_tier in {"critical", "high"}:
        return "scope_gap_to_rank"
    if expansion_status in {"hold_intake_only", "hold_process_regime_only"}:
        return "hold_current_posture"
    return "defer"


def build_matrix_family_coverage_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    payload = load_matrix_family_coverage_registry(file_path)
    families = []
    for raw_row in payload.get("families", []):
        row = dict(raw_row)
        row["support_class"] = _support_class(row)
        row["expansion_status"] = _expansion_status(row)
        row["primary_blocker"] = _primary_blocker(row)
        row["scope_priority"] = _scope_priority(row)
        row["artifact_count"] = len(row.get("artifacts", []))
        families.append(row)

    posture_counts: Dict[str, int] = {}
    category_counts: Dict[str, int] = {}
    support_class_counts: Dict[str, int] = {}
    for row in families:
        posture = str(row.get("runtime_posture", "unknown"))
        category = str(row.get("category", "unknown"))
        support_class = str(row.get("support_class", "unknown"))
        posture_counts[posture] = posture_counts.get(posture, 0) + 1
        category_counts[category] = category_counts.get(category, 0) + 1
        support_class_counts[support_class] = support_class_counts.get(support_class, 0) + 1

    explicit_supported = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) in {"quantitative_core", "directional_matrix", "qualitative_intake_only"}
    ]
    indirect_only = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) == "indirect_generic_support"
    ]
    open_gaps = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) == "open_gap"
    ]
    bounded_expansion_candidates = [
        row["matrix_family"]
        for row in families
        if str(row.get("expansion_status", "")) == "bounded_expansion_candidate"
    ]
    scope_hold_families = [
        row["matrix_family"]
        for row in families
        if str(row.get("expansion_status", "")) in {"hold_intake_only", "hold_process_regime_only"}
    ]
    evidence_blocked_families = [
        row["matrix_family"]
        for row in families
        if str(row.get("expansion_status", "")).startswith("blocked_on_")
    ]
    active_scope_priorities = [
        row["matrix_family"]
        for row in families
        if str(row.get("scope_priority", "")) in {"active_matrix_priority", "bounded_next_candidate"}
    ]
    scope_gap_priorities = [
        row["matrix_family"]
        for row in families
        if str(row.get("scope_priority", "")) == "scope_gap_to_rank"
    ]

    return {
        "summary": {
            "family_count": len(families),
            "posture_counts": dict(sorted(posture_counts.items())),
            "category_counts": dict(sorted(category_counts.items())),
            "support_class_counts": dict(sorted(support_class_counts.items())),
            "explicit_supported_families": explicit_supported,
            "indirect_only_families": indirect_only,
            "open_gap_families": open_gaps,
            "bounded_expansion_candidates": bounded_expansion_candidates,
            "scope_hold_families": scope_hold_families,
            "evidence_blocked_families": evidence_blocked_families,
            "active_scope_priorities": active_scope_priorities,
            "scope_gap_priorities": scope_gap_priorities,
            "policy": "matrix_family_scope_must_distinguish_explicit_support_from_generic_indirect_support",
            "expansion_policy": "do_not_broaden_matrix_scope_beyond_bounded_candidates_until_the_next_family_has_runtime_evidence_and_a_named_benchmark_or_calibration_landing",
        },
        "families": families,
    }


def render_matrix_family_coverage_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Family Coverage",
        "",
        "| Matrix Family | Category | Runtime Posture | Evidence Surface | Importance | Next Best Action |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {row.get('category', 'unknown')} | {row.get('runtime_posture', 'unknown')} | "
            f"{row.get('evidence_surface', 'unknown')} | {row.get('importance_tier', 'unknown')} | {row.get('next_best_action', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Expansion Gates",
        "",
        "| Matrix Family | Support Class | Expansion Status | Scope Priority | Primary Blocker | Artifacts |",
        "| --- | --- | --- | --- | --- | ---: |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {row.get('support_class', 'unknown')} | {row.get('expansion_status', 'unknown')} | {row.get('scope_priority', 'unknown')} | "
            f"{row.get('primary_blocker', 'none')} | {int(row.get('artifact_count', 0))} |"
        )

    lines.extend([
        "",
        "## Scope Notes",
        "",
        "| Matrix Family | Supported Today | Not Supported Today |",
        "| --- | --- | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {'; '.join(str(item) for item in row.get('what_is_supported', [])) or 'none'} | "
            f"{'; '.join(str(item) for item in row.get('what_is_not_supported', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Matrix families tracked: {int(summary.get('family_count', 0))}",
            f"Explicitly supported families: {', '.join(str(item) for item in summary.get('explicit_supported_families', [])) or 'none'}",
            f"Indirect-only families: {', '.join(str(item) for item in summary.get('indirect_only_families', [])) or 'none'}",
            f"Open-gap families: {', '.join(str(item) for item in summary.get('open_gap_families', [])) or 'none'}",
            f"Bounded expansion candidates: {', '.join(str(item) for item in summary.get('bounded_expansion_candidates', [])) or 'none'}",
            f"Scope-hold families: {', '.join(str(item) for item in summary.get('scope_hold_families', [])) or 'none'}",
            f"Evidence-blocked families: {', '.join(str(item) for item in summary.get('evidence_blocked_families', [])) or 'none'}",
            f"Active scope priorities: {', '.join(str(item) for item in summary.get('active_scope_priorities', [])) or 'none'}",
            f"Scope-gap priorities: {', '.join(str(item) for item in summary.get('scope_gap_priorities', [])) or 'none'}",
            f"Policy: {summary.get('policy', 'unknown')}",
            f"Expansion policy: {summary.get('expansion_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"