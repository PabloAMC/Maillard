from __future__ import annotations

from typing import Any, Dict, Mapping

from src.dha_lysinoalanine_external_package import build_dha_lysinoalanine_external_package_artifact
from src.extrusion_external_closure import build_extrusion_external_closure_artifact
from src.hexanal_nonanal_calibration import build_hexanal_nonanal_calibration_artifact
from src.pea_soy_external_evidence import build_pea_soy_external_evidence_artifact
from src.pea_soy_mixed_external_package import build_pea_soy_mixed_external_package_artifact


def _format_measurement(value: Any) -> str:
    if value is None:
        return ""
    return f"{float(value):.6g}"


def build_objective_progress_artifact() -> Dict[str, Any]:
    external_payload = build_pea_soy_external_evidence_artifact()
    calibration_payload = build_hexanal_nonanal_calibration_artifact()
    extrusion_payload = build_extrusion_external_closure_artifact()
    dha_package_payload = build_dha_lysinoalanine_external_package_artifact()
    mixed_package_payload = build_pea_soy_mixed_external_package_artifact()

    objectives = [
        {
            "objective_id": "external_mixed_meaty_positive_package",
            "label": "External mixed meaty-positive package",
            "target_count": int(external_payload.get("summary", {}).get("lane_count", 0)),
            "closed_count": int(external_payload.get("summary", {}).get("external_mixed_meaty_positive_ready", 0)),
            "status": "blocked_on_external_data",
            "resolved_in_last_step": [
                "explicit_required_external_package_for_pea_and_soy",
                "promotion_delta_if_package_lands_today",
            ],
            "prediction_effect": "No promotion-ready lane is unlocked yet; the repo now states exactly which package would move the decision gate.",
        },
        {
            "objective_id": "hexanal_nonanal_ambiguity",
            "label": "Hexanal/Nonanal ambiguity",
            "target_count": int(calibration_payload.get("summary", {}).get("marker_count", 0)),
            "closed_count": int(calibration_payload.get("summary", {}).get("closed_marker_count", 0)),
            "status": "closed_internal_calibration_route",
            "resolved_in_last_step": [
                "prediction_validation_chain_exposed",
                "closed_marker_counts_visible",
            ],
            "prediction_effect": "All tracked Hexanal/Nonanal markers are within the accepted internal ratio band, so adverse-marker drift is calibration-closed without upgrading promotion posture.",
        },
        {
            "objective_id": "extrusion_direct_damage_package",
            "label": "Extrusion direct damage closure package",
            "target_count": int(extrusion_payload.get("summary", {}).get("tracked_requirement_count", 0)),
            "closed_count": int(extrusion_payload.get("summary", {}).get("requirements_closed_today", 0)),
            "status": "blocked_on_external_data",
            "resolved_in_last_step": [
                "shared_extrusion_blocker_encoded_as_contract_artifact",
                "dha_lysinoalanine_external_package_specified",
            ],
            "prediction_effect": "Extrusion reporting now states the exact direct-marker package still missing for DHA/LAL and lysine-damage claims, and names the external measurement bundle required to move closure while keeping thiamine retention and soy thermal-history anchors explicitly partial.",
        },
    ]

    return {
        "summary": {
            "objective_count": len(objectives),
            "closed_objective_count": sum(1 for row in objectives if int(row.get("closed_count", 0)) >= int(row.get("target_count", 0)) and int(row.get("target_count", 0)) > 0),
            "policy": "objective_progress_must_state_what_changed_without_turning_internal_closure_into_external_promotion",
        },
        "objectives": objectives,
        "hexanal_nonanal_prediction_change": list(calibration_payload.get("prediction_change_cascade", [])),
        "pea_soy_external_delta": list(external_payload.get("lanes", [])),
        "pea_soy_mixed_external_package_delta": list(mixed_package_payload.get("matrices", [])),
        "extrusion_external_delta": list(extrusion_payload.get("matrices", [])),
        "dha_lysinoalanine_external_package_delta": list(dha_package_payload.get("matrices", [])),
        "selected_next_family_sprint": str(dha_package_payload.get("summary", {}).get("sprint_name", "unknown")),
    }


def render_objective_progress_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Objective Progress",
        "",
        "| Objective | Closed / Target | Status | Resolved In Last Step | Prediction Effect |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in payload.get("objectives", []):
        lines.append(
            f"| {row.get('label', 'unknown')} | {int(row.get('closed_count', 0))} / {int(row.get('target_count', 0))} | {row.get('status', 'unknown')} | {', '.join(str(item) for item in row.get('resolved_in_last_step', [])) or 'none'} | {row.get('prediction_effect', 'unknown')} |"
        )

    lines.extend(
        [
            "",
            "## Hexanal Nonanal Prediction Change",
            "",
            "| Protein | Compound | Internal2026 ppb | ProtocolPilot2026 ppb | Ratio | Closure State |",
            "| --- | --- | ---: | ---: | ---: | --- |",
        ]
    )
    for row in payload.get("hexanal_nonanal_prediction_change", []):
        ratio_value = row.get("ratio")
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {row.get('compound', 'unknown')} | {_format_measurement(row.get('internal2026_ppb'))} | {_format_measurement(row.get('protocolpilot2026_ppb'))} | {'' if ratio_value is None else f'{float(ratio_value):.3f}'} | {row.get('closure_state', 'unknown')} |"
        )

    lines.extend(
        [
            "",
            "## Pea Soy Promotion Delta",
            "",
            "| Protein | Mixed External Present Today | Would Close If Package Lands | Still Remaining After Package |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in payload.get("pea_soy_external_delta", []):
        effect = row.get("package_effect_if_landed_today", {})
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {'yes' if row.get('mixed_meaty_positive_external_present', False) else 'no'} | {', '.join(str(item) for item in effect.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in effect.get('remaining_after_package', [])) or 'none'} |"
        )

    lines.extend(
        [
            "",
            "## Pea Soy Mixed External Package Delta",
            "",
            "| Matrix | Status | Benchmark Candidate | Current Anchor | Desirable Targets | Adverse Targets | Would Close | Remaining After Package |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload.get("pea_soy_mixed_external_package_delta", []):
        lines.append(
            f"| {row.get('matrix', 'unknown')} | {row.get('package_status', 'unknown')} | {row.get('benchmark_candidate', 'unknown')} | {row.get('current_external_anchor', 'unknown')} | {', '.join(str(item) for item in row.get('desirable_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('adverse_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_after_package', [])) or 'none'} |"
        )

    lines.extend(
        [
            "",
            "## Extrusion Direct Closure Delta",
            "",
            "| Matrix | Direct Closure Ready | Missing Direct Markers | Supportive Anchors | Contextual Anchors | Remaining Requirements |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload.get("extrusion_external_delta", []):
        lines.append(
            f"| {row.get('matrix', 'unknown')} | {'yes' if row.get('direct_closure_ready', False) else 'no'} | {', '.join(str(item) for item in row.get('direct_markers_missing', [])) or 'none'} | {', '.join(str(item) for item in row.get('supportive_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('contextual_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_requirements', [])) or 'none'} |"
        )

    lines.extend(
        [
            "",
            "## DHA Lysinoalanine External Package Delta",
            "",
            "| Matrix | Status | Direct Damage Targets | Paired Meaty Targets | Metadata Required | Would Close | Remaining After Package |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in payload.get("dha_lysinoalanine_external_package_delta", []):
        lines.append(
            f"| {row.get('matrix', 'unknown')} | {row.get('package_status', 'unknown')} | {', '.join(str(item) for item in row.get('direct_damage_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('paired_meaty_positive_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('extrusion_metadata_fields', [])) or 'none'} | {', '.join(str(item) for item in row.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_after_package', [])) or 'none'} |"
        )

    lines.append("")
    lines.append(f"Selected next family sprint: {payload.get('selected_next_family_sprint', 'unknown')}")
    lines.append(f"Policy: {payload.get('summary', {}).get('policy', 'unknown')}")
    return "\n".join(lines) + "\n"