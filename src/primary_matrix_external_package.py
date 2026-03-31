from __future__ import annotations

from typing import Any, Dict, Mapping

from src.matrix_primary_benchmark_campaign import build_matrix_primary_benchmark_campaign_artifact
from src.pea_soy_external_evidence import build_pea_soy_external_evidence_artifact


def build_primary_matrix_external_package_artifact() -> Dict[str, Any]:
    campaign = build_matrix_primary_benchmark_campaign_artifact()
    external = build_pea_soy_external_evidence_artifact()

    selected_matrix = str(campaign.get("summary", {}).get("selected_matrix", "pea_iso"))
    selected_arm = next(
        (row for row in campaign.get("arms", []) if str(row.get("matrix", "")) == selected_matrix),
        campaign.get("arms", [])[0] if campaign.get("arms", []) else {},
    )
    external_row = next(
        (row for row in external.get("lanes", []) if str(row.get("protein_type", "")) == selected_matrix),
        {},
    )

    package = {
        "selected_matrix": selected_matrix,
        "benchmark_candidate": str(selected_arm.get("benchmark_id", "unknown")),
        "temperature_c": selected_arm.get("protocol_temperature_c"),
        "ph": selected_arm.get("protocol_ph"),
        "time_points_min": list(selected_arm.get("protocol_time_points_min", [])),
        "required_replicates": int(selected_arm.get("required_replicates", 0) or 0),
        "desirable_targets": list(selected_arm.get("required_desirable_targets", [])),
        "adverse_targets": list(selected_arm.get("required_adverse_targets", [])),
        "companion_assays": list(selected_arm.get("companion_assays", [])),
        "current_external_anchor": str(external_row.get("external_benchmark_id", "unknown")),
        "current_external_target_profile": str(external_row.get("external_target_profile", "unknown")),
        "would_close_requirements": list(selected_arm.get("would_close_requirements", [])),
        "remaining_after_package": list(selected_arm.get("remaining_after_protocol", [])),
        "promotion_blocker_today": str(selected_arm.get("promotion_blocker", "unknown")),
        "status": "specified_not_yet_measured",
        "policy": "landing_a_primary_external_package_means_specifying_the_exact_measurement_bundle_needed_for_one_prioritized_matrix_lane_not_claiming_that_external_closure_already_exists_or_replacing_the_dual_pea_soy_package",
    }

    return {
        "summary": {
            "selected_matrix": package["selected_matrix"],
            "status": package["status"],
            "benchmark_candidate": package["benchmark_candidate"],
            "policy": package["policy"],
        },
        "package": package,
    }


def render_primary_matrix_external_package_markdown(payload: Mapping[str, Any]) -> str:
    package = payload.get("package", {})
    lines = [
        "# Primary Matrix External Package",
        "",
        f"Selected matrix: {package.get('selected_matrix', 'unknown')}",
        f"Status: {package.get('status', 'unknown')}",
        f"Benchmark candidate: {package.get('benchmark_candidate', 'unknown')}",
        f"Current external anchor: {package.get('current_external_anchor', 'unknown')}",
        f"Current external target profile: {package.get('current_external_target_profile', 'unknown')}",
        f"Promotion blocker today: {package.get('promotion_blocker_today', 'unknown')}",
        "",
        "## Package Specification",
        "",
        f"Temperature C: {'' if package.get('temperature_c') is None else f'{float(package.get('temperature_c')):.1f}'}",
        f"pH: {'' if package.get('ph') is None else f'{float(package.get('ph')):.1f}'}",
        f"Time points min: {', '.join(str(item) for item in package.get('time_points_min', [])) or 'none'}",
        f"Required replicates: {int(package.get('required_replicates', 0))}",
        f"Desirable targets: {', '.join(str(item) for item in package.get('desirable_targets', [])) or 'none'}",
        f"Adverse targets: {', '.join(str(item) for item in package.get('adverse_targets', [])) or 'none'}",
        f"Companion assays: {', '.join(str(item) for item in package.get('companion_assays', [])) or 'none'}",
        "",
        "## Decision Delta",
        "",
        f"Would close requirements: {', '.join(str(item) for item in package.get('would_close_requirements', [])) or 'none'}",
        f"Still remaining after package: {', '.join(str(item) for item in package.get('remaining_after_package', [])) or 'none'}",
        f"Policy: {package.get('policy', 'unknown')}",
    ]
    return "\n".join(lines) + "\n"