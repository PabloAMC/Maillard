from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src import data_paths
from src.artifact_io import load_json_mapping, resolve_optional_path
from src.extrusion_external_closure import build_extrusion_external_closure_artifact


DEFAULT_DHA_LYSINOALANINE_EXTERNAL_PACKAGE_CONTRACT = data_paths.DHA_LYSINOALANINE_EXTERNAL_PACKAGE_CONTRACT


def load_dha_lysinoalanine_external_package_contract(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    return load_json_mapping(resolve_optional_path(file_path, DEFAULT_DHA_LYSINOALANINE_EXTERNAL_PACKAGE_CONTRACT))


def _unique_preserving_order(items: List[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def build_dha_lysinoalanine_external_package_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    contract = load_dha_lysinoalanine_external_package_contract(file_path)
    extrusion_payload = build_extrusion_external_closure_artifact()
    extrusion_by_matrix = {str(row.get("matrix", "")): row for row in extrusion_payload.get("matrices", [])}

    global_targets = [str(item) for item in contract.get("paired_meaty_positive_targets", [])]
    metadata_fields = [str(item) for item in contract.get("extrusion_metadata_fields", [])]
    would_close_requirements = [str(item) for item in contract.get("would_close_requirements", [])]
    remaining_after_package = [str(item) for item in contract.get("remaining_after_package", [])]
    matrix_rows = []

    for row in contract.get("matrices", []):
        matrix = str(row.get("matrix", "unknown"))
        extrusion_row = extrusion_by_matrix.get(matrix, {})
        marker_lookup = {
            str(marker.get("marker_id", "")): marker
            for marker in extrusion_row.get("direct_marker_panel", [])
            if str(marker.get("marker_id", ""))
        }
        direct_damage_targets = [str(item) for item in row.get("direct_damage_targets", [])]
        marker_rows = []
        preferred_assays: List[str] = []
        display_targets: List[str] = []
        for target in direct_damage_targets:
            marker = marker_lookup.get(target, {})
            assays = [str(item) for item in marker.get("preferred_assays", [])]
            preferred_assays.extend(assays)
            display_targets.append(str(marker.get("display_name", target)))
            marker_rows.append(
                {
                    "marker_id": target,
                    "display_name": str(marker.get("display_name", target)),
                    "preferred_assays": assays,
                    "missing_reason": str(marker.get("missing_reason", "unknown")),
                    "repo_direct_anchor_available": bool(marker.get("repo_direct_anchor_available", False)),
                }
            )

        matrix_rows.append(
            {
                "matrix": matrix,
                "display_name": str(row.get("display_name", matrix)),
                "benchmark_candidate": str(row.get("benchmark_candidate", "unknown")),
                "package_status": str(contract.get("package_status", "unknown")),
                "temperature_window_c": list(row.get("temperature_window_c", [])),
                "required_replicates": int(contract.get("required_replicates", 0) or 0),
                "extrusion_metadata_fields": metadata_fields,
                "paired_meaty_positive_targets": global_targets,
                "direct_damage_targets": display_targets,
                "direct_damage_target_rows": marker_rows,
                "preferred_assays": _unique_preserving_order(preferred_assays),
                "supportive_anchor_ids": list(extrusion_row.get("supportive_anchor_ids", [])),
                "contextual_anchor_ids": list(extrusion_row.get("contextual_anchor_ids", [])),
                "would_close_requirements": would_close_requirements,
                "remaining_after_package": remaining_after_package,
                "selection_reason": str(row.get("selection_reason", "unknown")),
                "current_closure_state": str(extrusion_row.get("closure_state", "unknown")),
            }
        )

    return {
        "summary": {
            "package_id": str(contract.get("package_id", "unknown")),
            "sprint_name": str(contract.get("sprint_name", "unknown")),
            "objective_id": str(contract.get("objective_id", "unknown")),
            "package_status": str(contract.get("package_status", "unknown")),
            "matrix_count": len(matrix_rows),
            "tracked_requirement_count": int(extrusion_payload.get("summary", {}).get("tracked_requirement_count", 0)),
            "root_blocker": str(contract.get("root_blocker", "unknown")),
            "policy": "specifying_the_dha_lysinoalanine_external_package_means_naming_the_exact_external_measurement_bundle_needed_without_claiming_that_direct_closure_exists_today",
        },
        "contract": contract,
        "matrices": matrix_rows,
    }


def render_dha_lysinoalanine_external_package_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    contract = payload.get("contract", {})
    lines = [
        "# DHA Lysinoalanine External Package",
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
        f"Extrusion metadata fields: {', '.join(str(item) for item in contract.get('extrusion_metadata_fields', [])) or 'none'}",
        f"Paired meaty-positive targets: {', '.join(str(item) for item in contract.get('paired_meaty_positive_targets', [])) or 'none'}",
        f"Would close requirements: {', '.join(str(item) for item in contract.get('would_close_requirements', [])) or 'none'}",
        f"Remaining after package: {', '.join(str(item) for item in contract.get('remaining_after_package', [])) or 'none'}",
        "",
        "## Matrix Specification",
        "",
        "| Matrix | Benchmark Candidate | Temp Window C | Direct Damage Targets | Preferred Assays | Supportive Anchors | Contextual Anchors | Would Close | Remaining After Package |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("matrices", []):
        lines.append(
            f"| {row.get('display_name', 'unknown')} | {row.get('benchmark_candidate', 'unknown')} | {', '.join(str(item) for item in row.get('temperature_window_c', [])) or 'none'} | {', '.join(str(item) for item in row.get('direct_damage_targets', [])) or 'none'} | {', '.join(str(item) for item in row.get('preferred_assays', [])) or 'none'} | {', '.join(str(item) for item in row.get('supportive_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('contextual_anchor_ids', [])) or 'none'} | {', '.join(str(item) for item in row.get('would_close_requirements', [])) or 'none'} | {', '.join(str(item) for item in row.get('remaining_after_package', [])) or 'none'} |"
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
        f"Tracked promotion requirements: {int(summary.get('tracked_requirement_count', 0))}",
        f"Policy: {summary.get('policy', 'unknown')}",
    ])
    return "\n".join(lines) + "\n"