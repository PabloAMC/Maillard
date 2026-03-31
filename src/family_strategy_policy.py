from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping

from src.chemistry_family_scope import load_chemistry_family_scope_registry
from src.family_ingestion_plan import load_family_ingestion_plan
from src.artifact_io import repo_root


def _sorted_scope_rows(payload: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows = [dict(row) for row in payload.get("families", []) if isinstance(row, Mapping)]
    rows.sort(key=lambda row: str(row.get("slr_family", "99")))
    return rows


def _sorted_plan_rows(payload: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows = [dict(row) for row in payload.get("families", []) if isinstance(row, Mapping)]
    rows.sort(key=lambda row: (int(row.get("implementation_wave", 99)), int(row.get("order_in_wave", 99)), str(row.get("slr_family", "99"))))
    return rows


def build_family_strategy_policy_artifact() -> Dict[str, Any]:
    scope = load_chemistry_family_scope_registry()
    plan = load_family_ingestion_plan()

    scope_rows = _sorted_scope_rows(scope)
    plan_rows = _sorted_plan_rows(plan)
    plan_by_family = {str(row.get("family_id", "")).strip(): row for row in plan_rows}

    first_class_core = [
        str(row.get("family_id", "")).strip()
        for row in scope_rows
        if str(row.get("current_status", "")).strip() == "first_class_core"
    ]
    bounded_lanes = [
        str(row.get("family_id", "")).strip()
        for row in scope_rows
        if str(row.get("current_status", "")).strip() == "bounded_lane"
    ]
    partially_encoded_high_priority = [
        str(row.get("family_id", "")).strip()
        for row in scope_rows
        if str(row.get("current_status", "")).strip() == "partially_encoded_high_priority"
    ]
    open_gaps = [
        str(row.get("family_id", "")).strip()
        for row in scope_rows
        if str(row.get("current_status", "")).strip() == "open_gap"
    ]

    next_family = next(iter(partially_encoded_high_priority), "")
    next_plan = plan_by_family.get(next_family, {})
    next_scope = next((row for row in scope_rows if str(row.get("family_id", "")).strip() == next_family), {})

    family_lane_positions = []
    for row in scope_rows:
        family_id = str(row.get("family_id", "")).strip()
        plan_row = plan_by_family.get(family_id, {})
        family_lane_positions.append(
            {
                "slr_family": str(row.get("slr_family", "")).zfill(2),
                "family_id": family_id,
                "display_name": str(row.get("display_name", "unknown")),
                "current_status": str(row.get("current_status", "unknown")),
                "priority": str(row.get("priority", "unknown")),
                "classification": (
                    "first_class_core"
                    if family_id in first_class_core
                    else "partially_encoded_high_priority"
                    if family_id in partially_encoded_high_priority
                    else "bounded_lane"
                    if family_id in bounded_lanes
                    else "open_gap"
                ),
                "preferred_payload_types": list(plan_row.get("preferred_payload_types", [])),
                "implementation_wave": int(plan_row.get("implementation_wave", 99)) if plan_row else 99,
            }
        )

    return {
        "source": {
            "chemistry_family_scope_registry": "data/lit/chemistry_family_scope_registry.json",
            "family_ingestion_plan_registry": "data/lit/family_ingestion_plan.json",
        },
        "summary": {
            "quantitative_trunk_family": "amino_acid_sugar_core",
            "benchmarked_foundation_statement": "keep_amino_acid_plus_sugar_core_with_strecker_and_sulfur_scoring_surfaces_as_the_benchmarked_foundation",
            "default_next_expansion_family": next_family,
            "default_next_expansion_reason": str(next_scope.get("why_it_matters", "")),
            "default_next_expansion_wave": int(next_plan.get("implementation_wave", 99)) if next_plan else 99,
            "family_lane_classification": {
                "first_class_core": first_class_core,
                "partially_encoded_high_priority": partially_encoded_high_priority,
                "bounded_lanes": bounded_lanes,
                "open_gaps": open_gaps,
            },
            "shared_ingestion_contract": {
                "machine_readable_only": True,
                "narrative_only_workflow_allowed": False,
                "required_surfaces": [
                    "benchmark_intake_registry",
                    "runtime_payloads",
                    "process_gap_registry",
                ],
                "policy": "new_families_use_the_same_machine_readable_ingestion_contract_as_the_core_and_do_not_create_a_parallel_markdown_only_workflow",
            },
            "lipid_crosstalk_dual_lane_policy": {
                "family_id": "lipid_oxidation_and_carbonylic_crosstalk",
                "observable_lane_payloads": ["benchmark_payload", "retention_payload"],
                "competition_lane_payloads": ["computational_prior", "structural_gap_entry"],
                "policy": "treat_lipid_oxidation_as_a_dual_lane_of_observable_adverse_markers_plus_carbonyl_competition_and_crosstalk_priors",
            },
            "compute_policy": {
                "mlp_status": "infrastructure_not_active_product_objective",
                "dft_policy": "selective_dft_is_reserved_for_benchmark_visible_sulfur_carbonyl_transfer_and_ts_sensitive_gaps_after_cheap_first_screening",
                "mlp_policy": "mlps_remain_bounded_offline_accelerators_until_local_geometry_or_ts_benchmarks_accept_them_on_maillard_relevant_systems",
            },
        },
        "family_lane_positions": family_lane_positions,
    }


def render_family_strategy_policy_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    classification = summary.get("family_lane_classification", {})
    ingestion_contract = summary.get("shared_ingestion_contract", {})
    dual_lane = summary.get("lipid_crosstalk_dual_lane_policy", {})
    compute_policy = summary.get("compute_policy", {})

    lines = [
        "# Family Strategy Policy",
        "",
        f"Quantitative trunk: {summary.get('quantitative_trunk_family', 'unknown')}",
        f"Benchmarked foundation: {summary.get('benchmarked_foundation_statement', 'unknown')}",
        f"Default next expansion family: {summary.get('default_next_expansion_family', 'unknown')}",
        f"Default next expansion reason: {summary.get('default_next_expansion_reason', 'unknown')}",
        "",
        "## Family Classification",
        "",
        f"First-class core: {', '.join(str(item) for item in classification.get('first_class_core', [])) or 'none'}",
        f"High-priority partial lanes: {', '.join(str(item) for item in classification.get('partially_encoded_high_priority', [])) or 'none'}",
        f"Bounded lanes: {', '.join(str(item) for item in classification.get('bounded_lanes', [])) or 'none'}",
        f"Open gaps: {', '.join(str(item) for item in classification.get('open_gaps', [])) or 'none'}",
        "",
        "## Ingestion Contract",
        "",
        f"Machine-readable only: {ingestion_contract.get('machine_readable_only', False)}",
        f"Narrative-only workflow allowed: {ingestion_contract.get('narrative_only_workflow_allowed', True)}",
        f"Required surfaces: {', '.join(str(item) for item in ingestion_contract.get('required_surfaces', [])) or 'none'}",
        f"Policy: {ingestion_contract.get('policy', 'unknown')}",
        "",
        "## Lipid Dual Lane",
        "",
        f"Family: {dual_lane.get('family_id', 'unknown')}",
        f"Observable lane payloads: {', '.join(str(item) for item in dual_lane.get('observable_lane_payloads', [])) or 'none'}",
        f"Competition lane payloads: {', '.join(str(item) for item in dual_lane.get('competition_lane_payloads', [])) or 'none'}",
        f"Policy: {dual_lane.get('policy', 'unknown')}",
        "",
        "## Compute Policy",
        "",
        f"MLP status: {compute_policy.get('mlp_status', 'unknown')}",
        f"DFT policy: {compute_policy.get('dft_policy', 'unknown')}",
        f"MLP policy: {compute_policy.get('mlp_policy', 'unknown')}",
        "",
        "## Lane Positions",
        "",
        "| SLR | Family | Classification | Priority | Payload Types | Wave |",
        "| --- | --- | --- | --- | --- | ---: |",
    ]
    for row in payload.get("family_lane_positions", []):
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {row.get('family_id', 'unknown')} | {row.get('classification', 'unknown')} | {row.get('priority', 'unknown')} | {', '.join(str(item) for item in row.get('preferred_payload_types', [])) or 'none'} | {row.get('implementation_wave', 'unknown')} |"
        )
    return "\n".join(lines) + "\n"