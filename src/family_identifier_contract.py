from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, Optional

from src import data_paths
from src.chemistry_family_scope import load_chemistry_family_scope_registry
from src.family_ingestion_plan import load_family_ingestion_plan
from src.literature_family_registry import (
    iter_benchmark_intake_entries,
    iter_computational_prior_entries,
    iter_flavor_reference_entries,
    iter_matrix_decision_panel_entries,
    iter_process_gap_entries,
    iter_retention_reference_entries,
)
from src.matrix_family_coverage import load_matrix_family_coverage_registry


def _sorted_plan_families(plan: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rows = [dict(row) for row in plan.get("families", []) if isinstance(row, Mapping)]
    rows.sort(key=lambda row: (int(row.get("implementation_wave", 99)), int(row.get("order_in_wave", 99)), str(row.get("slr_family", "99"))))
    return rows


def _nonempty_strings(values: Iterable[Any]) -> List[str]:
    return [str(value).strip() for value in values if str(value).strip()]


def build_family_identifier_contract_artifact() -> Dict[str, Any]:
    plan = load_family_ingestion_plan()
    scope = load_chemistry_family_scope_registry()
    matrix_scope = load_matrix_family_coverage_registry()

    plan_rows = _sorted_plan_families(plan)
    plan_family_ids = [str(row.get("family_id", "")).strip() for row in plan_rows]
    plan_slr_by_family = {str(row.get("family_id", "")).strip(): str(row.get("slr_family", "")).zfill(2) for row in plan_rows}

    scope_rows = [dict(row) for row in scope.get("families", []) if isinstance(row, Mapping)]
    scope_family_ids = [str(row.get("family_id", "")).strip() for row in scope_rows if str(row.get("family_id", "")).strip()]

    matrix_family_ids = [
        str(row.get("matrix_family", "")).strip()
        for row in matrix_scope.get("families", [])
        if isinstance(row, Mapping) and str(row.get("matrix_family", "")).strip()
    ]

    payload_groups = {
        "benchmark_intake": list(iter_benchmark_intake_entries()),
        "flavor_reference_payload": list(iter_flavor_reference_entries()),
        "retention_payload": list(iter_retention_reference_entries()),
        "computational_prior": list(iter_computational_prior_entries()),
        "structural_gap_entry": list(iter_process_gap_entries()),
        "decision_panel_entry": list(iter_matrix_decision_panel_entries()),
    }

    payload_usage: Dict[str, Dict[str, Any]] = {}
    payload_unknown_family_ids: List[str] = []
    for role, rows in payload_groups.items():
        family_ids = sorted({str(row.get("chemistry_family", "")).strip() for row in rows if str(row.get("chemistry_family", "")).strip()})
        missing_metadata = sum(1 for row in rows if not str(row.get("chemistry_family", "")).strip())
        payload_usage[role] = {
            "row_count": len(rows),
            "family_ids": family_ids,
            "missing_chemistry_family_count": int(missing_metadata),
        }
        for family_id in family_ids:
            if family_id not in plan_family_ids and family_id not in payload_unknown_family_ids:
                payload_unknown_family_ids.append(family_id)

    overlap = sorted(set(plan_family_ids).intersection(matrix_family_ids))
    plan_missing_from_scope = [family_id for family_id in plan_family_ids if family_id not in scope_family_ids]
    scope_missing_from_plan = [family_id for family_id in scope_family_ids if family_id not in plan_family_ids]

    family_rows = []
    for row in plan_rows:
        family_id = str(row.get("family_id", "")).strip()
        scope_match = next((entry for entry in scope_rows if str(entry.get("family_id", "")).strip() == family_id), {})
        family_rows.append(
            {
                "slr_family": plan_slr_by_family.get(family_id, "unknown"),
                "family_id": family_id,
                "display_name": str(row.get("display_name", "unknown")),
                "scope_registry_status": str(scope_match.get("current_status", "missing")),
                "scope_registry_priority": str(scope_match.get("priority", "missing")),
                "has_scope_registry_entry": bool(scope_match),
                "uses_canonical_scope_key": str(row.get("scope_family_id", "")).strip() == family_id,
            }
        )

    return {
        "source": {
            "family_ingestion_plan": data_paths.rel(data_paths.FAMILY_INGESTION_PLAN),
            "chemistry_family_scope_registry": data_paths.rel(data_paths.CHEMISTRY_FAMILY_SCOPE_REGISTRY),
            "matrix_family_coverage_registry": data_paths.rel(data_paths.MATRIX_FAMILY_COVERAGE_REGISTRY),
        },
        "summary": {
            "canonical_chemistry_family_count": len(plan_family_ids),
            "scope_registry_family_count": len(scope_family_ids),
            "matrix_family_count": len(matrix_family_ids),
            "plan_families_missing_from_scope": plan_missing_from_scope,
            "scope_families_missing_from_plan": scope_missing_from_plan,
            "payload_unknown_family_ids": sorted(payload_unknown_family_ids),
            "chemistry_matrix_axis_overlap": overlap,
            "scope_registry_covers_plan": not plan_missing_from_scope,
            "payload_surfaces_use_only_canonical_families": not payload_unknown_family_ids,
            "chemistry_and_matrix_axes_are_disjoint": not overlap,
            "identifier_policy": "chemistry_family_ids_are_canonical_across_scope_plan_payloads_and_validation_outputs",
            "axis_policy": "chemistry_family_and_matrix_family_remain_separate_axes_even_when_a_payload_references_both",
        },
        "families": family_rows,
        "payload_usage": payload_usage,
    }


def render_family_identifier_contract_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Identifier Contract",
        "",
        "| SLR | Family | In Scope Registry | Scope Status | Canonical Scope Key |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {row.get('family_id', 'unknown')} | {row.get('has_scope_registry_entry', False)} | {row.get('scope_registry_status', 'missing')} | {row.get('uses_canonical_scope_key', False)} |"
        )

    lines.extend([
        "",
        "## Payload Usage",
        "",
        "| Payload Role | Rows | Family IDs | Missing Chemistry Family |",
        "| --- | ---: | --- | ---: |",
    ])
    for role, row in payload.get("payload_usage", {}).items():
        lines.append(
            f"| {role} | {int(row.get('row_count', 0))} | {', '.join(str(item) for item in row.get('family_ids', [])) or 'none'} | {int(row.get('missing_chemistry_family_count', 0))} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Plan families missing from scope: {', '.join(str(item) for item in summary.get('plan_families_missing_from_scope', [])) or 'none'}",
            f"Scope families missing from plan: {', '.join(str(item) for item in summary.get('scope_families_missing_from_plan', [])) or 'none'}",
            f"Payload unknown family ids: {', '.join(str(item) for item in summary.get('payload_unknown_family_ids', [])) or 'none'}",
            f"Chemistry/matrix axis overlap: {', '.join(str(item) for item in summary.get('chemistry_matrix_axis_overlap', [])) or 'none'}",
            f"Identifier policy: {summary.get('identifier_policy', 'unknown')}",
            f"Axis policy: {summary.get('axis_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"