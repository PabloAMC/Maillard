from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from src.benchmark_validation import (
    build_matrix_observable_closure_audit,
    build_matrix_promotion_contract_artifact,
)
from src.hexanal_nonanal_calibration import build_hexanal_nonanal_calibration_artifact


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_PRIMARY_BENCHMARK_CONTRACT = _repo_root() / "data" / "protocols" / "ppi_spi_primary_benchmark_contract.json"


def load_primary_benchmark_contract(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_PRIMARY_BENCHMARK_CONTRACT
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _candidate_benchmark_id(matrix: str) -> str:
    mapping = {
        "pea_iso": "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
        "soy_iso": "soy_isolate_ribose_cysteine_100C_45min_Internal2026",
    }
    return mapping.get(str(matrix), f"{matrix}_mixed_matrix_candidate")


def build_matrix_primary_benchmark_campaign_artifact(
    contract_file: Optional[Path | str] = None,
) -> Dict[str, Any]:
    contract = load_primary_benchmark_contract(contract_file)
    promotion_payload = build_matrix_promotion_contract_artifact()
    closure_payload = build_matrix_observable_closure_audit()
    calibration_payload = build_hexanal_nonanal_calibration_artifact()

    promotion_lookup = {
        str(row.get("benchmark_id", "")): row
        for row in promotion_payload.get("benchmarks", [])
        if str(row.get("benchmark_id", "")).strip()
    }
    closure_lookup = {
        str(row.get("benchmark_id", "")): row
        for row in closure_payload.get("benchmarks", [])
        if str(row.get("benchmark_id", "")).strip()
    }
    calibration_lookup = {
        str(row.get("protein_type", "")): row
        for row in calibration_payload.get("lanes", [])
        if str(row.get("protein_type", "")).strip()
    }

    arm_rows = []
    for matrix_row in contract.get("matrices", []):
        matrix = str(matrix_row.get("matrix", "unknown"))
        benchmark_id = _candidate_benchmark_id(matrix)
        promotion_row = promotion_lookup.get(benchmark_id, {})
        closure_row = closure_lookup.get(benchmark_id, {})
        calibration_row = calibration_lookup.get(matrix, {})
        requirement_rows = list(promotion_row.get("requirements", []))
        unmet_requirements = [row for row in requirement_rows if not bool(row.get("passed", False))]
        compound_rows = list(closure_row.get("compounds", []))
        mechanistic_blockers = [
            str(row.get("compound", "unknown"))
            for row in compound_rows
            if str(row.get("closure_action", "")) == "mechanistic_blocker"
        ]
        evidence_or_calibration_blockers = [
            str(row.get("compound", "unknown"))
            for row in compound_rows
            if str(row.get("closure_action", "")) == "evidence_or_calibration_blocker"
        ]
        transfer_ready = [
            str(row.get("compound", "unknown"))
            for row in compound_rows
            if str(row.get("closure_action", "")) == "class_level_transfer_acceptable"
        ]
        arm_rows.append(
            {
                "matrix": matrix,
                "benchmark_id": benchmark_id,
                "process_state": str(promotion_row.get("process_state", matrix_row.get("format", "unknown"))),
                "target_profile": str(promotion_row.get("target_profile", "mixed")),
                "promotion_blocker": str(promotion_row.get("promotion_blocker", "unknown")),
                "protocol_temperature_c": float(matrix_row.get("temp_C", 0.0) or 0.0),
                "protocol_ph": float(matrix_row.get("ph", 0.0) or 0.0),
                "protocol_time_points_min": list(matrix_row.get("time_points_min", [])),
                "required_replicates": int(contract.get("analytical_requirements", {}).get("replicates", 0) or 0),
                "required_desirable_targets": list(contract.get("required_panel", {}).get("desirable", [])),
                "required_adverse_targets": list(contract.get("required_panel", {}).get("adverse", [])),
                "required_safety_targets": list(contract.get("required_panel", {}).get("safety", [])),
                "companion_assays": list(contract.get("companion_assays", [])),
                "calibration_closure_action": str(calibration_row.get("closure_action", "missing_data")),
                "calibration_passed": bool(calibration_row.get("passed", False)),
                "hexanal_ratio": calibration_row.get("compounds", {}).get("Hexanal", {}).get("ratio"),
                "nonanal_ratio": calibration_row.get("compounds", {}).get("Nonanal", {}).get("ratio"),
                "unmet_promotion_requirements": unmet_requirements,
                "mechanistic_blockers": mechanistic_blockers,
                "evidence_or_calibration_blockers": evidence_or_calibration_blockers,
                "transfer_ready_targets": transfer_ready,
                "would_close_requirements": [
                    row["key"]
                    for row in unmet_requirements
                    if row["key"] in {
                        "comparator_is_measured_volatiles",
                        "external_quantitative_origin",
                        "minimum_quantitative_closed_targets",
                    }
                ],
                "remaining_after_protocol": [
                    row["key"]
                    for row in unmet_requirements
                    if row["key"] == "no_internal_or_directional_dependencies"
                ],
            }
        )

    selected_target = promotion_payload.get("selected_promotion_target", {}) or {}
    selected_matrix = str(selected_target.get("protein_type", "pea_iso"))
    selected_arm = next((row for row in arm_rows if row["matrix"] == selected_matrix), arm_rows[0] if arm_rows else None)
    fallback_arm = next((row for row in arm_rows if selected_arm is not None and row["matrix"] != selected_arm["matrix"]), None)

    return {
        "summary": {
            "protocol_id": str(contract.get("protocol_id", "unknown")),
            "selected_matrix": selected_arm.get("matrix", "unknown") if selected_arm else "unknown",
            "selected_benchmark_id": selected_arm.get("benchmark_id", "unknown") if selected_arm else "unknown",
            "fallback_matrix": fallback_arm.get("matrix", "none") if fallback_arm else "none",
            "arm_count": len(arm_rows),
            "primary_data_blocker": "external_quantitative_measured_volatiles_for_mixed_matrix_lane",
            "next_best_action": "run_primary_benchmark_protocol_and_land_results_as_benchmark_json_plus_process_state_calibration",
        },
        "contract": contract,
        "selected_target": selected_target,
        "arms": arm_rows,
    }


def render_matrix_primary_benchmark_campaign_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Matrix Primary Benchmark Campaign",
        "",
        f"Protocol id: {summary.get('protocol_id', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Selected benchmark id: {summary.get('selected_benchmark_id', 'unknown')}",
        f"Fallback matrix: {summary.get('fallback_matrix', 'none')}",
        f"Primary data blocker: {summary.get('primary_data_blocker', 'unknown')}",
        f"Next best action: {summary.get('next_best_action', 'unknown')}",
        "",
        "## Campaign Arms",
        "",
        "| Matrix | Benchmark | Temp C | pH | Time Points | Calibration Route | Promotion Blocker | Would Close | Remaining After Protocol |",
        "| --- | --- | ---: | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("arms", []):
        lines.append(
            f"| {row.get('matrix', 'unknown')} | {row.get('benchmark_id', 'unknown')} | {float(row.get('protocol_temperature_c', 0.0)):.1f} | "
            f"{float(row.get('protocol_ph', 0.0)):.1f} | {', '.join(str(item) for item in row.get('protocol_time_points_min', [])) or 'none'} | "
            f"{row.get('calibration_closure_action', 'unknown')} | "
            f"{row.get('promotion_blocker', 'unknown')} | {', '.join(str(item) for item in row.get('would_close_requirements', [])) or 'none'} | "
            f"{', '.join(str(item) for item in row.get('remaining_after_protocol', [])) or 'none'} |"
        )

    lines.extend([
        "",
        "## Required Panel",
        "",
        "| Matrix | Transfer-Ready Targets | Evidence/Calibration Blockers | Mechanistic Blockers | Hexanal Ratio | Nonanal Ratio | Companion Assays | Replicates |",
        "| --- | --- | --- | --- | ---: | ---: | --- | ---: |",
    ])
    for row in payload.get("arms", []):
        lines.append(
            f"| {row.get('matrix', 'unknown')} | {', '.join(str(item) for item in row.get('transfer_ready_targets', [])) or 'none'} | "
            f"{', '.join(str(item) for item in row.get('evidence_or_calibration_blockers', [])) or 'none'} | "
            f"{', '.join(str(item) for item in row.get('mechanistic_blockers', [])) or 'none'} | {'' if row.get('hexanal_ratio') is None else f'{float(row.get('hexanal_ratio')):.3f}'} | {'' if row.get('nonanal_ratio') is None else f'{float(row.get('nonanal_ratio')):.3f}'} | "
            f"{', '.join(str(item) for item in row.get('companion_assays', [])) or 'none'} | {int(row.get('required_replicates', 0))} |"
        )

    lines.extend([
        "",
        "## Promotion Delta",
        "",
        "| Matrix | Requirement | Passed Today | Detail |",
        "| --- | --- | --- | --- |",
    ])
    for row in payload.get("arms", []):
        for requirement in row.get("unmet_promotion_requirements", []):
            lines.append(
                f"| {row.get('matrix', 'unknown')} | {requirement.get('label', 'unknown')} | {requirement.get('passed', False)} | {requirement.get('detail', 'unknown')} |"
            )
    return "\n".join(lines) + "\n"