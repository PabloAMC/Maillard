from __future__ import annotations

from typing import Any, Dict, Mapping

from src.hexanal_nonanal_calibration import build_hexanal_nonanal_calibration_artifact
from src.pea_soy_external_evidence import build_pea_soy_external_evidence_artifact


def build_hexanal_nonanal_resolution_artifact() -> Dict[str, Any]:
    calibration = build_hexanal_nonanal_calibration_artifact()
    external = build_pea_soy_external_evidence_artifact()
    external_by_protein = {str(row.get("protein_type", "unknown")): row for row in external.get("lanes", [])}

    resolution_rows = []
    for row in calibration.get("prediction_change_cascade", []):
        protein_type = str(row.get("protein_type", "unknown"))
        external_row = external_by_protein.get(protein_type, {})
        resolution_rows.append(
            {
                "protein_type": protein_type,
                "compound": str(row.get("compound", "unknown")),
                "ratio": row.get("ratio"),
                "ambiguity_status": "reduced_by_internal_calibration_closure" if str(row.get("closure_state", "")) == "calibration_closed" else "still_open",
                "promotion_effect": "no_external_promotion_upgrade",
                "remaining_external_blocker": str(external_row.get("next_best_action", "unknown")),
            }
        )

    return {
        "summary": {
            "marker_count": len(resolution_rows),
            "reduced_marker_count": sum(1 for row in resolution_rows if row.get("ambiguity_status") == "reduced_by_internal_calibration_closure"),
            "policy": "hexanal_nonanal_ambiguity_can_be_reduced_by_internal_calibration_closure_without_relabeling_the_lane_as_externally_promotion_ready",
        },
        "markers": resolution_rows,
    }


def render_hexanal_nonanal_resolution_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Hexanal Nonanal Resolution",
        "",
        "| Protein | Compound | Ratio | Ambiguity Status | Promotion Effect | Remaining External Blocker |",
        "| --- | --- | ---: | --- | --- | --- |",
    ]
    for row in payload.get("markers", []):
        ratio_value = row.get("ratio")
        lines.append(
            f"| {row.get('protein_type', 'unknown')} | {row.get('compound', 'unknown')} | {'' if ratio_value is None else f'{float(ratio_value):.3f}'} | {row.get('ambiguity_status', 'unknown')} | {row.get('promotion_effect', 'unknown')} | {row.get('remaining_external_blocker', 'unknown')} |"
        )
    lines.append("")
    lines.append(f"Policy: {payload.get('summary', {}).get('policy', 'unknown')}")
    return "\n".join(lines) + "\n"