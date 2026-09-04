from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from src import data_paths

_IMPACT_BY_POSTURE = {
    "first_class_runtime_lane": 1.00,
    "immediate_expansion_lane": 0.92,
    "upstream_pretreatment_lane": 0.84,
    "high_value_support_lane": 0.82,
    "guardrail_lane": 0.74,
    "matrix_scope_lane": 0.68,
    "failure_mode_lane": 0.64,
    "upstream_precursor_sink": 0.60,
    "bounded_upstream_source": 0.58,
    "trapping_burden_modifier": 0.54,
    "first_class_core": 0.20,
}

_ACTION_TEMPLATES = {
    "02": "split adverse lipid markers from carbonyl-competition priors",
    "03": "promote thiamine from availability cue to calibrated sulfur-support modifier",
    "04": "encode nucleotide and ribose enrichment as explicit upstream state variables",
    "05": "bound peptide and glutathione support with reusable sulfur-lane priors",
    "06": "extend matrix-family scope without over-claiming transferred closure",
    "07": "land donor identity as a benchmark-facing formulation variable",
    "10": "bind fermentation pretreatment state to precursor release and pH-shift effects",
    "11": "wire lipid-Maillard competition constants into the runtime competition surface",
    "12": "expand safety damage proxies with benchmark-backed guardrail payloads",
    "14": "bound ascorbic-acid dicarbonyl source terms before they leak into safety claims",
}

_RUNTIME_QUEUE_POSTURE_PRIORITY = {
    "immediate_expansion_lane": 0,
    "first_class_runtime_lane": 1,
    "high_value_support_lane": 2,
    "guardrail_lane": 3,
    "upstream_pretreatment_lane": 4,
    "matrix_scope_lane": 5,
    "failure_mode_lane": 6,
    "upstream_precursor_sink": 7,
    "bounded_upstream_source": 8,
    "trapping_burden_modifier": 9,
    "first_class_core": 10,
}


def load_family_ingestion_plan(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else data_paths.FAMILY_INGESTION_PLAN
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_deep_research_backlog(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else data_paths.DEEP_RESEARCH_BACKLOG
    if not path.exists():
        return {"summary": {}, "items": []}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


def _normalize_source_file(path_value: Any) -> str:
    return Path(str(path_value or "")).name


def _select_backlog_items_for_family(family_row: Mapping[str, Any], backlog_items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source_file = _normalize_source_file(family_row.get("source_slr_file", ""))
    if not source_file:
        return []
    return [
        item
        for item in backlog_items
        if item.get("status") == "BACKLOG" and source_file in {str(file_name) for file_name in item.get("files", [])}
    ]


def _select_all_items_for_family(family_row: Mapping[str, Any], backlog_items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    source_file = _normalize_source_file(family_row.get("source_slr_file", ""))
    if not source_file:
        return []
    return [
        item
        for item in backlog_items
        if source_file in {str(file_name) for file_name in item.get("files", [])}
    ]


def _compute_impact_score(family_row: Mapping[str, Any], backlog_items: list[dict[str, Any]]) -> float:
    posture = str(family_row.get("strategic_posture", "unknown"))
    impact = float(_IMPACT_BY_POSTURE.get(posture, 0.50))
    if any("benchmark_payload" == str(payload) for payload in family_row.get("preferred_payload_types", [])):
        impact += 0.04
    if any(int(item.get("score_value", 0)) >= 8 for item in backlog_items):
        impact += 0.04
    if str(family_row.get("lane_concept", "")) == "computational_closure_candidate":
        impact += 0.03
    return _clamp(impact, 0.10, 1.00)


def _compute_ease_score(family_row: Mapping[str, Any], backlog_items: list[dict[str, Any]]) -> float:
    wave = int(family_row.get("implementation_wave", 99) or 99)
    dependency_count = len(list(family_row.get("depends_on", [])))
    target_module_count = len(list(family_row.get("target_runtime_modules", [])))
    ease = 0.98
    ease -= 0.10 * max(0, wave - 1)
    ease -= 0.05 * dependency_count
    ease -= 0.03 * max(0, target_module_count - 2)
    if any("benchmark_payload" == str(payload) for payload in family_row.get("preferred_payload_types", [])):
        ease += 0.05
    if str(family_row.get("runtime_concept", "")).endswith("lane"):
        ease += 0.02
    if not backlog_items:
        ease -= 0.10
    return _clamp(ease, 0.20, 1.00)


def _compute_evidence_score(backlog_items: list[dict[str, Any]]) -> float:
    if not backlog_items:
        return 0.0
    high_confidence = sum(1 for item in backlog_items if int(item.get("score_value", 0)) >= 8)
    evidence = 0.35
    evidence += 0.40 * min(high_confidence, 3) / 3.0
    evidence += 0.25 * min(len(backlog_items), 4) / 4.0
    return _clamp(evidence, 0.0, 1.0)


def _build_recommended_slice(family_row: Mapping[str, Any], backlog_items: list[dict[str, Any]]) -> Dict[str, Any]:
    slr_family = str(family_row.get("slr_family", "unknown"))
    citations = [str(item.get("citation", "")) for item in backlog_items[:3] if str(item.get("citation", "")).strip()]
    return {
        "focus": _ACTION_TEMPLATES.get(slr_family, str((family_row.get("next_curation_actions") or ["continue machine-readable runtime landing"])[0])),
        "target_modules": list(family_row.get("target_runtime_modules", []))[:3],
        "candidate_citations": citations,
        "preferred_payload_types": list(family_row.get("preferred_payload_types", [])),
    }


def build_family_ingestion_plan_artifact(
    file_path: Optional[Path | str] = None,
    backlog_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    payload = load_family_ingestion_plan(file_path)
    deep_research_backlog = load_deep_research_backlog(backlog_path)
    families = [dict(row) for row in payload.get("families", [])]
    backlog_items = [dict(row) for row in deep_research_backlog.get("items", [])]

    payload_type_counts: Dict[str, int] = {}
    posture_counts: Dict[str, int] = {}
    wave_counts: Dict[str, int] = {}
    mapped_scope_count = 0
    unmapped_scope_families = []
    backlog_family_count = 0
    total_family_backlog_citations = 0

    families.sort(key=lambda row: (int(row.get("implementation_wave", 99)), int(row.get("order_in_wave", 99)), str(row.get("slr_family", "99"))))

    for row in families:
        posture = str(row.get("strategic_posture", "unknown"))
        posture_counts[posture] = posture_counts.get(posture, 0) + 1

        wave = str(row.get("implementation_wave", "unknown"))
        wave_counts[wave] = wave_counts.get(wave, 0) + 1

        if str(row.get("scope_family_id", "")).strip():
            mapped_scope_count += 1
        else:
            unmapped_scope_families.append(str(row.get("slr_family", "unknown")))

        for payload_type in row.get("preferred_payload_types", []):
            key = str(payload_type)
            payload_type_counts[key] = payload_type_counts.get(key, 0) + 1

        family_backlog_items = _select_backlog_items_for_family(row, backlog_items)
        family_all_items = _select_all_items_for_family(row, backlog_items)
        backlog_count = len(family_backlog_items)
        backlog_8_of_8_count = sum(1 for item in family_all_items if int(item.get("score_value", 0)) >= 8)
        backlog_occurrence_count = sum(int(item.get("occurrence_count", 0)) for item in family_all_items)
        if backlog_count > 0:
            backlog_family_count += 1
            total_family_backlog_citations += backlog_count

        impact_score = _compute_impact_score(row, family_backlog_items)
        ease_score = _compute_ease_score(row, family_backlog_items)
        evidence_score = _compute_evidence_score(family_backlog_items)
        priority_score = round(100.0 * (0.45 * impact_score + 0.30 * ease_score + 0.25 * evidence_score), 1)

        row["deep_research_backlog"] = {
            "citation_count": backlog_count,
            "high_confidence_count": backlog_8_of_8_count,
            "occurrence_count": backlog_occurrence_count,
            "top_citations": [str(item.get("citation", "")) for item in family_backlog_items[:3]],
        }
        row["execution_priority"] = {
            "impact_score": round(impact_score, 3),
            "ease_score": round(ease_score, 3),
            "evidence_score": round(evidence_score, 3),
            "priority_score": float(priority_score),
            "recommended_slice": _build_recommended_slice(row, family_backlog_items),
        }

    first_wave = [
        str(row.get("slr_family", "unknown"))
        for row in families
        if int(row.get("implementation_wave", 99)) == 1
    ]
    active_sequence = [str(row.get("slr_family", "unknown")) for row in families]
    prioritized_families = sorted(
        [
            row
            for row in families
            if int(row.get("implementation_wave", 99)) > 0
            and int(row.get("deep_research_backlog", {}).get("citation_count", 0)) > 0
        ],
        key=lambda row: (
            int(_RUNTIME_QUEUE_POSTURE_PRIORITY.get(str(row.get("strategic_posture", "unknown")), 99)),
            int(row.get("implementation_wave", 99)),
            int(row.get("order_in_wave", 99)),
            -float(row.get("execution_priority", {}).get("priority_score", 0.0)),
        ),
    )
    recommended_runtime_queue = [str(row.get("slr_family", "unknown")) for row in prioritized_families[:5]]
    next_family = prioritized_families[0] if prioritized_families else None

    return {
        "source": str((Path(file_path) if file_path is not None else data_paths.FAMILY_INGESTION_PLAN).resolve().relative_to(data_paths.REPO_ROOT).as_posix()),
        "deep_research_source": str((Path(backlog_path) if backlog_path is not None else data_paths.DEEP_RESEARCH_BACKLOG).resolve().relative_to(data_paths.REPO_ROOT).as_posix()),
        "template_source": str(payload.get("template_source", "")),
        "summary": {
            "family_count": len(families),
            "payload_type_counts": dict(sorted(payload_type_counts.items())),
            "strategic_posture_counts": dict(sorted(posture_counts.items())),
            "implementation_wave_counts": dict(sorted(wave_counts.items())),
            "recommended_first_wave": first_wave,
            "active_build_sequence": active_sequence,
            "backlog_family_count": backlog_family_count,
            "backlog_citation_count": total_family_backlog_citations,
            "recommended_runtime_queue": recommended_runtime_queue,
            "recommended_next_family": None if next_family is None else {
                "slr_family": str(next_family.get("slr_family", "unknown")),
                "family_id": str(next_family.get("family_id", "unknown")),
                "display_name": str(next_family.get("display_name", "unknown")),
                "priority_score": float(next_family.get("execution_priority", {}).get("priority_score", 0.0)),
                "recommended_slice": dict(next_family.get("execution_priority", {}).get("recommended_slice", {})),
            },
            "mapped_scope_family_count": mapped_scope_count,
            "unmapped_scope_families": unmapped_scope_families,
            "state_of_art_policy": "extend_the_quantitative_core_by_explicit_family_lanes_with_machine_readable_payloads_not_narrative_only_docs",
            "identifier_policy": "scope_family_id_uses_the_same_canonical_chemistry_family_id_as_payloads_and_validation_outputs",
            "axis_policy": "chemistry_family_and_matrix_family_are_separate_axes_and_should_not_share_identifier_spaces",
        },
        "families": families,
    }


def render_family_ingestion_plan_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Ingestion Plan",
        "",
        "| SLR | Family | Strategic Posture | Runtime Concept | Payload Types | Wave | Next Build Action |",
        "| --- | --- | --- | --- | --- | ---: | --- |",
    ]
    for row in payload.get("families", []):
        next_action = "; ".join(str(item) for item in row.get("next_curation_actions", [])[:1]) or "none"
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {row.get('family_id', 'unknown')} | {row.get('strategic_posture', 'unknown')} | {row.get('runtime_concept', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('preferred_payload_types', [])) or 'none'} | {row.get('implementation_wave', 'unknown')} | {next_action} |"
        )

    lines.extend([
        "",
        "## Deep Research Priority Surface",
        "",
        "| Rank | SLR | Family | Backlog Citations | 8/8 Citations | Impact | Ease | Priority | Next Slice |",
        "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ])
    ranked_rows = sorted(
        [row for row in payload.get("families", []) if int(row.get("deep_research_backlog", {}).get("citation_count", 0)) > 0],
        key=lambda row: -float(row.get("execution_priority", {}).get("priority_score", 0.0)),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        backlog = row.get("deep_research_backlog", {})
        priority = row.get("execution_priority", {})
        next_slice = str(priority.get("recommended_slice", {}).get("focus", "continue runtime landing"))
        lines.append(
            f"| {rank} | {row.get('slr_family', 'unknown')} | {row.get('display_name', row.get('family_id', 'unknown'))} | "
            f"{int(backlog.get('citation_count', 0))} | {int(backlog.get('high_confidence_count', 0))} | "
            f"{float(priority.get('impact_score', 0.0)):.2f} | {float(priority.get('ease_score', 0.0)):.2f} | "
            f"{float(priority.get('priority_score', 0.0)):.1f} | {next_slice} |"
        )

    lines.extend([
        "",
        "## Target Modules",
        "",
        "| SLR | Target Modules | Target Observables or State Variables |",
        "| --- | --- | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {'; '.join(str(item) for item in row.get('target_runtime_modules', [])) or 'none'} | "
            f"{'; '.join(str(item) for item in row.get('target_compounds_or_state_variables', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    next_family = summary.get("recommended_next_family") or {}
    next_slice = next_family.get("recommended_slice", {}) if isinstance(next_family, Mapping) else {}
    lines.extend(
        [
            "",
            f"Families tracked: {int(summary.get('family_count', 0))}",
            f"Recommended first wave: {', '.join(str(item) for item in summary.get('recommended_first_wave', [])) or 'none'}",
            f"Active build sequence: {', '.join(str(item) for item in summary.get('active_build_sequence', [])) or 'none'}",
            f"Backlog-bearing families: {int(summary.get('backlog_family_count', 0))}",
            f"Backlog citations mapped to family plan: {int(summary.get('backlog_citation_count', 0))}",
            f"Recommended runtime queue: {', '.join(str(item) for item in summary.get('recommended_runtime_queue', [])) or 'none'}",
            (
                f"Recommended next slice: {next_family.get('slr_family', 'none')} {next_family.get('display_name', '')} "
                f"-> {next_slice.get('focus', 'none')}"
            ).strip(),
            f"Mapped scope families: {int(summary.get('mapped_scope_family_count', 0))}",
            f"Unmapped scope families: {', '.join(str(item) for item in summary.get('unmapped_scope_families', [])) or 'none'}",
            f"Policy: {summary.get('state_of_art_policy', 'unknown')}",
            f"Identifier policy: {summary.get('identifier_policy', 'unknown')}",
            f"Axis policy: {summary.get('axis_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"