from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from src import data_paths
from src.artifact_io import load_json_mapping, resolve_optional_path
from src.primary_benchmark_campaign import build_matrix_primary_benchmark_campaign
from src.pea_soy_external_evidence import build_pea_soy_external_evidence_artifact


DEFAULT_PEA_SOY_MIXED_EXTERNAL_PACKAGE_CONTRACT = data_paths.PEA_SOY_MIXED_EXTERNAL_PACKAGE_CONTRACT


def load_pea_soy_mixed_external_package_contract(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    return load_json_mapping(resolve_optional_path(file_path, DEFAULT_PEA_SOY_MIXED_EXTERNAL_PACKAGE_CONTRACT))


def build_pea_soy_mixed_external_package_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    contract = load_pea_soy_mixed_external_package_contract(file_path)
    campaign = build_matrix_primary_benchmark_campaign()
    external = build_pea_soy_external_evidence_artifact()

    campaign_by_matrix = {str(row.get("matrix", "")): row for row in campaign.get("arms", [])}
    external_by_matrix = {str(row.get("protein_type", "")): row for row in external.get("lanes", [])}
    matrix_rows = []

    for row in contract.get("matrices", []):
        matrix = str(row.get("matrix", "unknown"))
        campaign_row = campaign_by_matrix.get(matrix, {})
        external_row = external_by_matrix.get(matrix, {})
        matrix_rows.append(
            {
                "matrix": matrix,
                "display_name": str(row.get("display_name", matrix)),
                "benchmark_candidate": str(campaign_row.get("benchmark_id", "unknown")),
                "package_status": str(contract.get("package_status", "unknown")),
                "temperature_c": campaign_row.get("protocol_temperature_c"),
                "ph": campaign_row.get("protocol_ph"),
                "time_points_min": list(campaign_row.get("protocol_time_points_min", [])),
                "required_replicates": int(contract.get("required_replicates", 0) or 0),
                "desirable_targets": list(campaign_row.get("required_desirable_targets", contract.get("global_desirable_targets", []))),
                "adverse_targets": list(campaign_row.get("required_adverse_targets", contract.get("global_adverse_targets", []))),
                "companion_assays": list(campaign_row.get("companion_assays", contract.get("global_companion_assays", []))),
                "current_external_anchor": str(external_row.get("external_benchmark_id", "unknown")),
                "current_external_target_profile": str(external_row.get("external_target_profile", "unknown")),
                "mixed_meaty_positive_external_present": bool(external_row.get("mixed_meaty_positive_external_present", False)),
                "would_close_requirements": list(contract.get("would_close_requirements", [])),
                "remaining_after_package": list(contract.get("remaining_after_package", [])),
                "promotion_blocker_today": str(campaign_row.get("promotion_blocker", "unknown")),
                "selection_reason": str(row.get("selection_reason", "unknown")),
            }
        )

    return {
        "summary": {
            "package_id": str(contract.get("package_id", "unknown")),
            "sprint_name": str(contract.get("sprint_name", "unknown")),
            "objective_id": str(contract.get("objective_id", "unknown")),
            "package_status": str(contract.get("package_status", "unknown")),
            "matrix_count": len(matrix_rows),
            "root_blocker": str(contract.get("root_blocker", "unknown")),
            "policy": "specifying_the_pea_soy_mixed_external_package_means_naming_the_exact_external_measurement_bundle_for_both_matrices_without_claiming_external_closure_today",
        },
        "contract": contract,
        "matrices": matrix_rows,
    }


def render_pea_soy_mixed_external_package_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    contract = payload.get("contract", {})
    lines = [
        "# Pea Soy Mixed External Package",
        "",
        f"Package id: {summary.get('package_id', 'unknown')}",
        f"Sprint name: {summary.get('sprint_name', 'unknown')}",
        f"Objective id: {summary.get('objective_id', 'unknown')}",
        f"Status: {summary.get('package_status', 'unknown')}",
        f"Root blocker: {summary.get('root_blocker', 'unknown')}",
        f"Purpose: {contract.get('purpose', 'unknown')}",
        "",
        "## Global Measurement Bundle",
        "",
        f"Required replicates: {int(contract.get('required_replicates', 0) or 0)}",
        f"Desirable targets: {', '.join(str(item) for item in contract.get('global_desirable_targets', [])) or 'none'}",
        f"Adverse targets: {', '.join(str(item) for item in contract.get('global_adverse_targets', [])) or 'none'}",
        f"Companion assays: {', '.join(str(item) for item in contract.get('global_companion_assays', [])) or 'none'}",
        f"Would close requirements: {', '.join(str(item) for item in contract.get('would_close_requirements', [])) or 'none'}",
        f"Remaining after package: {', '.join(str(item) for item in contract.get('remaining_after_package', [])) or 'none'}",
        "",
        "## Matrix Specification",
        "",
        "| Matrix | Benchmark Candidate | Temp C | pH | Time Points | Current Anchor | Current Profile | Desirable Targets | Adverse Targets | Companion Assays | Would Close | Remaining After Package |",
        "| --- | --- | ---: | ---: | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("matrices", []):
        temperature = row.get("temperature_c")
        ph = row.get("ph")
        lines.append(
            f"| {row.get('display_name', 'unknown')} | {row.get('benchmark_candidate', 'unknown')} | {'' if temperature is None else f'{float(temperature):.1f}'} | {'' if ph is None else f'{float(ph):.1f}'} | {', '.join(str(item) for item in row.get('time_points_min', [])) or 'none'} | {row.get('current_external_anchor', 'unknown')} | {row.get('current_external_target_profile', 'unknown')} | {', '.join(str(item) for item in row.get('desirable_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('adverse_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('companion_assays', [])) or 'none'} | {', '.join(str(item) for item in row.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_after_package', [])) or 'none'} |"
        )

    lines.extend([
        "",
        "## Matrix Rationale",
        "",
    ])
    for row in payload.get("matrices", []):
        lines.append(f"- {row.get('display_name', 'unknown')}: {row.get('selection_reason', 'unknown')}")

    lines.extend([
        "",
        f"Matrices tracked: {int(summary.get('matrix_count', 0))}",
        f"Policy: {summary.get('policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"