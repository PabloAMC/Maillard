from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src.artifact_io import load_json_mapping, repo_root, resolve_optional_path


DEFAULT_EXTRUSION_EXTERNAL_CLOSURE_CONTRACT = repo_root() / "data" / "protocols" / "extrusion_external_closure_contract.json"
PROCESS_STATE_CALIBRATIONS = repo_root() / "data" / "lit" / "process_state_calibrations.json"


def load_extrusion_external_closure_contract(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    return load_json_mapping(resolve_optional_path(file_path, DEFAULT_EXTRUSION_EXTERNAL_CLOSURE_CONTRACT))


def _process_state_lookup() -> Dict[str, Mapping[str, Any]]:
    payload = load_json_mapping(PROCESS_STATE_CALIBRATIONS)
    return {
        str(row.get("id", "")).strip(): row
        for row in payload.get("entries", [])
        if str(row.get("id", "")).strip()
    }


def _requirement_passed(requirement_key: str, markers: List[Mapping[str, Any]]) -> bool:
    if requirement_key == "direct_crosslink_marker_external_quantified":
        return any(
            bool(marker.get("repo_direct_anchor_available", False))
            for marker in markers
            if str(marker.get("closure_role", "")) == "direct_crosslink_marker"
        )
    if requirement_key == "lysine_damage_marker_external_quantified":
        return any(
            bool(marker.get("repo_direct_anchor_available", False))
            for marker in markers
            if str(marker.get("closure_role", "")) == "lysine_damage_marker"
        )
    return False


def build_extrusion_external_closure_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    contract = load_extrusion_external_closure_contract(file_path)
    process_lookup = _process_state_lookup()
    promotion_requirements = list(contract.get("promotion_requirements", []))

    backfill_rows = []
    backfill_by_matrix: Dict[str, Dict[str, List[str]]] = {}
    for row in contract.get("literature_backfill", []):
        anchor_id = str(row.get("anchor_id", "")).strip()
        source_row = process_lookup.get(anchor_id, {})
        available = bool(source_row)
        matrix_scope = [str(item) for item in row.get("matrix_scope", [])]
        rendered = {
            "anchor_id": anchor_id,
            "matrix_scope": matrix_scope,
            "support_level": str(row.get("support_level", "unknown")),
            "source_kind": str(row.get("source_kind", "unknown")),
            "available_in_repo": available,
            "validated_status": str(source_row.get("validated_status", "missing")) if available else "missing",
            "provenance_tier": str(source_row.get("provenance_tier", "missing")) if available else "missing",
            "what_it_supports": list(row.get("what_it_supports", [])),
            "why_not_sufficient": str(row.get("why_not_sufficient", "unknown")),
            "source_artifact": str(row.get("source_artifact", "unknown")),
        }
        backfill_rows.append(rendered)
        for matrix in matrix_scope:
            buckets = backfill_by_matrix.setdefault(matrix, {"supportive": [], "contextual": []})
            if str(row.get("support_level", "")).startswith("supportive"):
                buckets["supportive"].append(anchor_id)
            else:
                buckets["contextual"].append(anchor_id)

    matrix_rows = []
    for row in contract.get("matrices", []):
        matrix = str(row.get("matrix", "unknown"))
        markers = list(row.get("direct_marker_panel", []))
        requirement_rows = []
        for requirement in promotion_requirements:
            key = str(requirement.get("key", ""))
            passed = _requirement_passed(key, markers)
            requirement_rows.append(
                {
                    "key": key,
                    "label": str(requirement.get("label", key)),
                    "passed": passed,
                    "description": str(requirement.get("description", "")),
                }
            )

        direct_markers_missing = [
            str(marker.get("display_name", marker.get("marker_id", "unknown")))
            for marker in markers
            if not bool(marker.get("repo_direct_anchor_available", False)) and bool(marker.get("required", True))
        ]
        helper = backfill_by_matrix.get(matrix, {"supportive": [], "contextual": []})
        matrix_rows.append(
            {
                "matrix": matrix,
                "display_name": str(row.get("display_name", matrix)),
                "process_state": str(contract.get("process_state", "unknown")),
                "temperature_window_c": list(row.get("temperature_window_c", [])),
                "direct_closure_ready": len(direct_markers_missing) == 0,
                "direct_markers_missing": direct_markers_missing,
                "direct_marker_panel": markers,
                "supportive_anchor_ids": helper.get("supportive", []),
                "contextual_anchor_ids": helper.get("contextual", []),
                "companion_contract_id": str(row.get("companion_contract_id", "unknown")),
                "promotion_requirements": requirement_rows,
                "closed_requirements": [item["key"] for item in requirement_rows if bool(item["passed"])],
                "remaining_requirements": [item["key"] for item in requirement_rows if not bool(item["passed"])],
                "next_best_action": "land_external_extrusion_damage_package",
                "closure_state": "blocked_on_external_package" if direct_markers_missing else "ready_for_promotion_review",
            }
        )

    direct_closure_ready_count = sum(1 for row in matrix_rows if bool(row.get("direct_closure_ready", False)))
    supportive_anchor_count = sum(
        len(row.get("supportive_anchor_ids", [])) + len(row.get("contextual_anchor_ids", []))
        for row in matrix_rows
    )

    return {
        "summary": {
            "protocol_id": str(contract.get("protocol_id", "unknown")),
            "matrix_count": len(matrix_rows),
            "tracked_requirement_count": len(promotion_requirements),
            "requirements_closed_today": 0,
            "direct_closure_ready_matrices": direct_closure_ready_count,
            "supportive_anchor_count": supportive_anchor_count,
            "root_blocker": "external_quantitative_direct_damage_markers_under_extrusion",
            "next_best_action": "land_external_extrusion_damage_package",
            "selected_next_family_sprint": "dha_lysinoalanine_external_package",
            "selection_rationale": "After the shared extrusion blocker is made explicit, DHA/Lysinoalanine remains the next chemistry sprint because it is the highest-impact lane still blocked primarily by missing direct extrusion marker packages rather than by internal visibility or governance gaps.",
        },
        "contract": contract,
        "matrices": matrix_rows,
        "literature_backfill": backfill_rows,
    }


def render_extrusion_external_closure_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Extrusion External Closure",
        "",
        f"Protocol id: {summary.get('protocol_id', 'unknown')}",
        f"Root blocker: {summary.get('root_blocker', 'unknown')}",
        f"Next best action: {summary.get('next_best_action', 'unknown')}",
        f"Selected next family sprint: {summary.get('selected_next_family_sprint', 'unknown')}",
        f"Selection rationale: {summary.get('selection_rationale', 'unknown')}",
        "",
        "## Matrix Status",
        "",
        "| Matrix | Process State | Temp Window C | Direct Closure Ready | Missing Direct Markers | Supportive Anchors | Contextual Anchors | Remaining Requirements |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("matrices", []):
        temp_window = row.get("temperature_window_c", [])
        lines.append(
            f"| {row.get('display_name', 'unknown')} | {row.get('process_state', 'unknown')} | {', '.join(str(item) for item in temp_window) or 'none'} | {'yes' if row.get('direct_closure_ready', False) else 'no'} | {', '.join(str(item) for item in row.get('direct_markers_missing', [])) or 'none'} | {', '.join(str(item) for item in row.get('supportive_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('contextual_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_requirements', [])) or 'none'} |"
        )

    lines.extend([
        "",
        "## Literature Backfill",
        "",
        "| Anchor | Matrix Scope | Support Level | Available In Repo | Supports | Why Not Sufficient |",
        "| --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload.get("literature_backfill", []):
        lines.append(
            f"| {row.get('anchor_id', 'unknown')} | {', '.join(str(item) for item in row.get('matrix_scope', [])) or 'none'} | {row.get('support_level', 'unknown')} | {'yes' if row.get('available_in_repo', False) else 'no'} | {', '.join(str(item) for item in row.get('what_it_supports', [])) or 'none'} | {row.get('why_not_sufficient', 'unknown')} |"
        )

    lines.extend([
        "",
        f"Matrices tracked: {int(summary.get('matrix_count', 0))}",
        f"Supportive/contextual anchors encoded: {int(summary.get('supportive_anchor_count', 0))}",
        f"Direct-closure-ready matrices: {int(summary.get('direct_closure_ready_matrices', 0))}",
    ])
    return "\n".join(lines) + "\n"