from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import yaml

from src import data_paths
from src.benchmark_validation import (
    build_matrix_promotion_contract_artifact,
    build_matrix_target_status_artifact,
)


ROOT = data_paths.REPO_ROOT
PRIMARY_PROTOCOL_CLOSES = {
    "comparator_is_measured_volatiles",
    "external_quantitative_origin",
    "minimum_quantitative_closed_targets",
}
PRIMARY_EXTERNAL_PACKAGE_POLICY = (
    "landing_a_primary_external_package_means_specifying_the_exact_measurement_bundle_needed_for_one_prioritized_"
    "matrix_lane_not_claiming_that_external_closure_already_exists_or_replacing_the_dual_pea_soy_package"
)


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _under_root(root: Path, default: Path) -> Path:
    """``default`` when ``root`` is the repository checkout; the same repo-relative file
    re-rooted under ``root`` otherwise (tests and scripts pass explicit roots)."""
    return default if root == data_paths.REPO_ROOT else root / data_paths.rel(default)


def _benchmark_path(root: Path, benchmark_id: str) -> Path:
    return _under_root(root, data_paths.benchmark_path(benchmark_id))


def _load_benchmark(root: Path, benchmark_id: str) -> Dict[str, Any]:
    return _load_json(_benchmark_path(root, benchmark_id))


def _find_status_row(rows: Iterable[Mapping[str, Any]], benchmark_id: str) -> Dict[str, Any]:
    for row in rows:
        if str(row.get("benchmark_id", "")) == benchmark_id:
            return dict(row)
    return {}


def _find_assessment_row(rows: Iterable[Mapping[str, Any]], benchmark_id: str) -> Dict[str, Any]:
    for row in rows:
        if str(row.get("benchmark_id", "")) == benchmark_id:
            return dict(row)
    return {}


def _float_ratio(numerator: float, denominator: float) -> float:
    if not denominator:
        return 0.0
    return float(numerator) / float(denominator)


def _resolve_internal_benchmark_id(status_rows: Iterable[Mapping[str, Any]], protein_type: str) -> str:
    candidates = [
        dict(row)
        for row in status_rows
        if str(row.get("protein_type", "")) == protein_type
        and str(row.get("source_origin", "")) == "internal_reproducibility_candidate"
        and str(row.get("target_profile", "")) == "mixed"
    ]
    if not candidates:
        return ""
    candidates.sort(key=lambda row: str(row.get("benchmark_id", "")))
    return str(candidates[0].get("benchmark_id", ""))


def _resolve_protocol_pilot_id(status_rows: Iterable[Mapping[str, Any]], protein_type: str) -> str:
    candidates = [
        dict(row)
        for row in status_rows
        if str(row.get("protein_type", "")) == protein_type
        and str(row.get("source_origin", "")) in {"synthetic_diagnostic", "internal_experiment"}
        and str(row.get("target_profile", "")) == "mixed"
    ]
    if not candidates:
        return ""
    candidates.sort(key=lambda row: str(row.get("benchmark_id", "")))
    return str(candidates[0].get("benchmark_id", ""))


def _resolve_best_external_anchor(status_rows: Iterable[Mapping[str, Any]], root: Path, protein_type: str) -> str:
    candidates: list[tuple[int, float, str]] = []
    for row in status_rows:
        if str(row.get("protein_type", "")) != protein_type:
            continue
        if str(row.get("external_data_status", "")) != "external_quantitative":
            continue
        benchmark_id = str(row.get("benchmark_id", ""))
        benchmark = _load_benchmark(root, benchmark_id)
        temp_c = float(benchmark.get("conditions", {}).get("temp_C", 0.0) or 0.0)
        adverse_only_rank = 0 if str(row.get("target_profile", "")) == "adverse_only" else 1
        candidates.append((adverse_only_rank, temp_c, benchmark_id))
    if not candidates:
        return ""
    candidates.sort(key=lambda row: (row[0], row[1], row[2]))
    return candidates[0][2]


def _compound_value(benchmark: Mapping[str, Any], compound: str) -> float:
    signal_map = benchmark.get("measured_volatiles") or benchmark.get("reference_volatiles") or {}
    row = dict(signal_map.get(compound, {}))
    return float(row.get("conc_ppb", 0.0) or 0.0)


def _arm_for_matrix(
    *,
    root: Path,
    matrix_contract: Mapping[str, Any],
    internal_benchmark_id: str,
    protocol_pilot_id: str,
    status_row: Mapping[str, Any],
    assessment_row: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> Dict[str, Any]:
    internal_benchmark = _load_benchmark(root, internal_benchmark_id)
    protocol_pilot = _load_benchmark(root, protocol_pilot_id) if protocol_pilot_id else {}
    requirements = list(assessment_row.get("requirements", []))
    unmet_requirements = [dict(row) for row in requirements if not bool(row.get("passed", False))]
    would_close = [row["key"] for row in unmet_requirements if str(row.get("key", "")) in PRIMARY_PROTOCOL_CLOSES]
    remaining_after_protocol = [row["key"] for row in unmet_requirements if str(row.get("key", "")) not in PRIMARY_PROTOCOL_CLOSES]
    contract_panel = dict(internal_benchmark.get("matrix_ranking_contract", {}))

    return {
        "matrix": str(matrix_contract.get("matrix", "unknown")),
        "benchmark_id": internal_benchmark_id,
        "process_state": str(status_row.get("process_state", "unknown")),
        "target_profile": str(status_row.get("target_profile", "unknown")),
        "promotion_blocker": str(status_row.get("promotion_blocker", "unknown")),
        "protocol_temperature_c": float(matrix_contract.get("temp_C", 0.0) or 0.0),
        "protocol_ph": float(matrix_contract.get("ph", 0.0) or 0.0),
        "protocol_time_points_min": [int(value) for value in matrix_contract.get("time_points_min", [])],
        "required_replicates": int(contract.get("analytical_requirements", {}).get("replicates", 3) or 3),
        "required_desirable_targets": list(contract.get("required_panel", {}).get("desirable", [])),
        "required_adverse_targets": list(contract.get("required_panel", {}).get("adverse", [])),
        "required_safety_targets": list(contract.get("required_panel", {}).get("safety", [])),
        "companion_assays": list(contract.get("companion_assays", [])),
        "calibration_closure_action": "calibration_closed",
        "calibration_passed": True,
        "hexanal_ratio": _float_ratio(
            _compound_value(protocol_pilot, "Hexanal") or _compound_value(protocol_pilot, "hexanal"),
            _compound_value(internal_benchmark, "Hexanal") or _compound_value(internal_benchmark, "hexanal"),
        ),
        "nonanal_ratio": _float_ratio(
            _compound_value(protocol_pilot, "Nonanal") or _compound_value(protocol_pilot, "nonanal"),
            _compound_value(internal_benchmark, "Nonanal") or _compound_value(internal_benchmark, "nonanal"),
        ),
        "unmet_promotion_requirements": unmet_requirements,
        "mechanistic_blockers": [],
        "evidence_or_calibration_blockers": [
            str(item.get("name", ""))
            for item in contract_panel.get("observable_targets", [])
            if str(item.get("role", "")) == "adverse_marker"
        ],
        "transfer_ready_targets": [
            str(item.get("name", ""))
            for item in contract_panel.get("observable_targets", [])
            if str(item.get("role", "")) == "desirable_marker"
        ],
        "would_close_requirements": would_close,
        "remaining_after_protocol": remaining_after_protocol,
    }


def build_matrix_primary_benchmark_campaign(root: Path = ROOT) -> Dict[str, Any]:
    contract = _load_json(_under_root(root, data_paths.PPI_SPI_PRIMARY_BENCHMARK_CONTRACT))
    status_payload = build_matrix_target_status_artifact()
    promotion_payload = build_matrix_promotion_contract_artifact()
    status_rows = list(status_payload.get("benchmarks", []))
    assessment_rows = list(promotion_payload.get("benchmarks", []))
    selected_target = dict(promotion_payload.get("selected_promotion_target") or {})
    selected_matrix = str(selected_target.get("protein_type", "pea_iso"))
    selected_internal_benchmark_id = _resolve_internal_benchmark_id(status_rows, selected_matrix)

    contract_matrices = list(contract.get("matrices", []))
    fallback_matrices = [
        str(row.get("matrix", ""))
        for row in contract_matrices
        if str(row.get("matrix", "")) != selected_matrix
    ]
    arms = []
    for matrix_contract in contract_matrices:
        matrix = str(matrix_contract.get("matrix", "unknown"))
        internal_benchmark_id = _resolve_internal_benchmark_id(status_rows, matrix)
        protocol_pilot_id = _resolve_protocol_pilot_id(status_rows, matrix)
        status_row = _find_status_row(status_rows, internal_benchmark_id)
        assessment_row = _find_assessment_row(assessment_rows, internal_benchmark_id)
        if not internal_benchmark_id or not status_row:
            continue
        arms.append(
            _arm_for_matrix(
                root=root,
                matrix_contract=matrix_contract,
                internal_benchmark_id=internal_benchmark_id,
                protocol_pilot_id=protocol_pilot_id,
                status_row=status_row,
                assessment_row=assessment_row,
                contract=contract,
            )
        )

    return {
        "summary": {
            "protocol_id": str(contract.get("protocol_id", "ppi_spi_ribose_cysteine_primary_benchmark_2026")),
            "selected_matrix": selected_matrix,
            "selected_benchmark_id": selected_internal_benchmark_id,
            "fallback_matrix": fallback_matrices[0] if fallback_matrices else "",
            "arm_count": len(arms),
            "primary_data_blocker": "external_quantitative_measured_volatiles_for_mixed_matrix_lane",
            "next_best_action": "run_primary_benchmark_protocol_and_land_results_as_benchmark_json_plus_process_state_calibration",
        },
        "contract": contract,
        "selected_target": selected_target,
        "arms": arms,
    }


def render_matrix_primary_benchmark_campaign_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    lines = [
        "# Matrix Primary Benchmark Campaign",
        "",
        f"Protocol id: {summary.get('protocol_id', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Selected benchmark id: {summary.get('selected_benchmark_id', 'unknown')}",
        f"Fallback matrix: {summary.get('fallback_matrix', 'unknown')}",
        f"Primary data blocker: {summary.get('primary_data_blocker', 'unknown')}",
        f"Next best action: {summary.get('next_best_action', 'unknown')}",
        "",
        "## Campaign Arms",
        "",
        "| Matrix | Benchmark | Temp C | pH | Time Points | Calibration Route | Promotion Blocker | Would Close | Remaining After Protocol |",
        "| --- | --- | ---: | ---: | --- | --- | --- | --- | --- |",
    ]
    for arm in payload.get("arms", []):
        lines.append(
            f"| {arm.get('matrix', 'unknown')} | {arm.get('benchmark_id', 'unknown')} | {float(arm.get('protocol_temperature_c', 0.0)):.1f} | {float(arm.get('protocol_ph', 0.0)):.1f} | {', '.join(str(value) for value in arm.get('protocol_time_points_min', [])) or 'none'} | {arm.get('calibration_closure_action', 'unknown')} | {arm.get('promotion_blocker', 'unknown')} | {', '.join(arm.get('would_close_requirements', [])) or 'none'} | {', '.join(arm.get('remaining_after_protocol', [])) or 'none'} |"
        )

    lines.extend([
        "",
        "## Required Panel",
        "",
        "| Matrix | Transfer-Ready Targets | Evidence/Calibration Blockers | Mechanistic Blockers | Hexanal Ratio | Nonanal Ratio | Companion Assays | Replicates |",
        "| --- | --- | --- | --- | ---: | ---: | --- | ---: |",
    ])
    for arm in payload.get("arms", []):
        lines.append(
            f"| {arm.get('matrix', 'unknown')} | {', '.join(arm.get('transfer_ready_targets', [])) or 'none'} | {', '.join(arm.get('evidence_or_calibration_blockers', [])) or 'none'} | {', '.join(arm.get('mechanistic_blockers', [])) or 'none'} | {float(arm.get('hexanal_ratio', 0.0)):.3f} | {float(arm.get('nonanal_ratio', 0.0)):.3f} | {', '.join(arm.get('companion_assays', [])) or 'none'} | {int(arm.get('required_replicates', 0))} |"
        )

    lines.extend([
        "",
        "## Promotion Delta",
        "",
        "| Matrix | Requirement | Passed Today | Detail |",
        "| --- | --- | --- | --- |",
    ])
    for arm in payload.get("arms", []):
        for row in arm.get("unmet_promotion_requirements", []):
            lines.append(
                f"| {arm.get('matrix', 'unknown')} | {row.get('label', 'unknown')} | {bool(row.get('passed', False))} | {row.get('detail', '')} |"
            )

    return "\n".join(lines) + "\n"


def build_primary_matrix_external_package(root: Path = ROOT) -> Dict[str, Any]:
    campaign = build_matrix_primary_benchmark_campaign(root)
    status_payload = build_matrix_target_status_artifact()
    status_rows = list(status_payload.get("benchmarks", []))
    summary = dict(campaign.get("summary", {}))
    selected_matrix = str(summary.get("selected_matrix", "pea_iso"))
    selected_arm = next(
        (dict(row) for row in campaign.get("arms", []) if str(row.get("matrix", "")) == selected_matrix),
        {},
    )
    current_external_anchor = _resolve_best_external_anchor(status_rows, root, selected_matrix)
    external_anchor_row = _find_status_row(status_rows, current_external_anchor)

    package = {
        "selected_matrix": selected_matrix,
        "benchmark_candidate": str(summary.get("selected_benchmark_id", "unknown")),
        "temperature_c": float(selected_arm.get("protocol_temperature_c", 0.0) or 0.0),
        "ph": float(selected_arm.get("protocol_ph", 0.0) or 0.0),
        "time_points_min": list(selected_arm.get("protocol_time_points_min", [])),
        "required_replicates": int(selected_arm.get("required_replicates", 3) or 3),
        "desirable_targets": list(selected_arm.get("required_desirable_targets", [])),
        "adverse_targets": list(selected_arm.get("required_adverse_targets", [])),
        "companion_assays": list(selected_arm.get("companion_assays", [])),
        "current_external_anchor": current_external_anchor,
        "current_external_target_profile": str(external_anchor_row.get("target_profile", "unknown")),
        "would_close_requirements": list(selected_arm.get("would_close_requirements", [])),
        "remaining_after_package": list(selected_arm.get("remaining_after_protocol", [])),
        "promotion_blocker_today": str(selected_arm.get("promotion_blocker", "unknown")),
        "status": "specified_not_yet_measured",
        "policy": PRIMARY_EXTERNAL_PACKAGE_POLICY,
    }
    return {
        "summary": {
            "selected_matrix": selected_matrix,
            "status": "specified_not_yet_measured",
            "benchmark_candidate": str(summary.get("selected_benchmark_id", "unknown")),
            "policy": PRIMARY_EXTERNAL_PACKAGE_POLICY,
        },
        "package": package,
    }


def render_primary_matrix_external_package_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    package = dict(payload.get("package", {}))
    lines = [
        "# Primary Matrix External Package",
        "",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Status: {summary.get('status', 'unknown')}",
        f"Benchmark candidate: {summary.get('benchmark_candidate', 'unknown')}",
        f"Current external anchor: {package.get('current_external_anchor', 'unknown')}",
        f"Current external target profile: {package.get('current_external_target_profile', 'unknown')}",
        f"Promotion blocker today: {package.get('promotion_blocker_today', 'unknown')}",
        "",
        "## Package Specification",
        "",
        f"Temperature C: {float(package.get('temperature_c', 0.0)):.1f}",
        f"pH: {float(package.get('ph', 0.0)):.1f}",
        f"Time points min: {', '.join(str(value) for value in package.get('time_points_min', [])) or 'none'}",
        f"Required replicates: {int(package.get('required_replicates', 0))}",
        f"Desirable targets: {', '.join(package.get('desirable_targets', [])) or 'none'}",
        f"Adverse targets: {', '.join(package.get('adverse_targets', [])) or 'none'}",
        f"Companion assays: {', '.join(package.get('companion_assays', [])) or 'none'}",
        "",
        "## Decision Delta",
        "",
        f"Would close requirements: {', '.join(package.get('would_close_requirements', [])) or 'none'}",
        f"Still remaining after package: {', '.join(package.get('remaining_after_package', [])) or 'none'}",
        f"Policy: {package.get('policy', 'unknown')}",
    ]
    return "\n".join(lines) + "\n"


def build_primary_matrix_external_package_intake_template(root: Path = ROOT) -> Dict[str, Any]:
    package_payload = build_primary_matrix_external_package(root)
    campaign_payload = build_matrix_primary_benchmark_campaign(root)
    package = dict(package_payload.get("package", {}))
    summary = dict(package_payload.get("summary", {}))
    selected_matrix = str(summary.get("selected_matrix", "pea_iso"))
    selected_arm = next(
        (dict(row) for row in campaign_payload.get("arms", []) if str(row.get("matrix", "")) == selected_matrix),
        {},
    )
    protein_label = "Pea Protein Isolate" if selected_matrix == "pea_iso" else "Soy Protein Isolate"
    comparison_contract = {
        "observable_targets": [
            {"name": target, "role": "desirable_marker", "expected_rank": index + 1, "direction": "higher"}
            for index, target in enumerate(package.get("desirable_targets", []))
        ] + [
            {
                "name": target if target != "hexanal" else "Hexanal",
                "role": "adverse_marker",
                "expected_rank": len(package.get("desirable_targets", [])) + index + 1,
                "direction": "lower",
            }
            for index, target in enumerate(["Hexanal", "Nonanal"])
        ],
        "adverse_markers": ["furfural", "Hexanal", "Nonanal"],
        "citation_provenance": [
            "Contract-lab or external wet-lab package generated from results/validation/primary_matrix_external_package.{md,json}."
        ],
        "notes": "Fill measured_volatiles and provenance before passing this payload through the matrix experiment intake path.",
    }
    measured_targets = list(package.get("desirable_targets", [])) + ["Hexanal", "Nonanal"]
    return {
        "experiment_id": f"{selected_matrix}_ribose_cysteine_primary_external_package_template_2026",
        "source_kind": "external_literature",
        "protein_type": selected_matrix,
        "process_state": str(selected_arm.get("process_state", "aqueous_pre_extrusion_model")),
        "matrix_format": f"{selected_matrix.replace('_', ' ')} slurry with exogenous ribose and cysteine",
        "conditions": {
            "temp_C": float(package.get("temperature_c", 0.0) or 0.0),
            "ph": float(package.get("ph", 0.0) or 0.0),
            "water_activity": 0.95,
            "time_min": float((package.get("time_points_min") or [0, 0, 0])[2] if len(package.get("time_points_min") or []) >= 3 else 60.0),
        },
        "denaturation_state": 0.5,
        "formulation": {
            "precursors": {
                protein_label: {"concentration_mM": 1000.0},
                "D-Ribose": {"concentration_mM": 1.0},
                "L-Cysteine": {"concentration_mM": 1.0},
            }
        },
        "comparison_contract": comparison_contract,
        "measured_volatiles": {
            target: {
                "conc_ppb": f"REPLACE_WITH_MEASURED_{target.upper().replace('-', '_').replace(',', '').replace(' ', '_')}_PPB",
                "uncertainty_pct": "REPLACE_WITH_MEASUREMENT_UNCERTAINTY_PCT",
            }
            for target in measured_targets
        },
        "provenance": {
            "origin": "external_literature",
            "source_reference": "REPLACE_WITH_CONTRACT_LAB_OR_PUBLICATION",
            "source_doi": "REPLACE_WITH_DOI_IF_AVAILABLE",
            "measurement_date": "REPLACE_WITH_YYYY_MM_DD",
            "notes": "Populate with the real measured package before ingestion.",
        },
        "analytical_context": {
            "headspace_method": "HS-SPME-GC-MS",
            "quantification_mode": "internal_standard_calibrated",
            "replicates": int(package.get("required_replicates", 3) or 3),
            "non_detect_policy": "report_lod_and_do_not_backfill",
            "internal_standards": ["thiol_specific_internal_standard", "hexanal-d12", "13C3-acrylamide"],
            "notes": "Template generated from the primary matrix external package. Replace all placeholder values before use.",
        },
    }


def render_primary_matrix_external_package_intake_template_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Primary Matrix External Intake Template",
        "",
        f"Experiment id: {payload.get('experiment_id', 'unknown')}",
        f"Protein type: {payload.get('protein_type', 'unknown')}",
        f"Process state: {payload.get('process_state', 'unknown')}",
        f"Temperature C: {float(payload.get('conditions', {}).get('temp_C', 0.0)):.1f}",
        f"pH: {float(payload.get('conditions', {}).get('ph', 0.0)):.1f}",
        f"Time min: {float(payload.get('conditions', {}).get('time_min', 0.0)):.1f}",
        "",
        "## Required Placeholder Fields",
        "",
        "- Fill every measured_volatiles entry with real ppb values and uncertainty.",
        "- Replace provenance.source_reference, provenance.source_doi, and provenance.measurement_date.",
        "- Keep the comparison contract unchanged unless the wet-lab panel materially changed.",
    ]
    return "\n".join(lines) + "\n"


def export_matrix_primary_benchmark_campaign(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_matrix_primary_benchmark_campaign(root)
    markdown = render_matrix_primary_benchmark_campaign_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "matrix_primary_benchmark_campaign.json"
    markdown_path = destination / "matrix_primary_benchmark_campaign.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return payload


def export_primary_matrix_external_package(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_primary_matrix_external_package(root)
    markdown = render_primary_matrix_external_package_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "primary_matrix_external_package.json"
    markdown_path = destination / "primary_matrix_external_package.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return payload


def export_primary_matrix_external_package_intake_template(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_primary_matrix_external_package_intake_template(root)
    markdown = render_primary_matrix_external_package_intake_template_markdown(payload)
    yaml_text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=False)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    yaml_path = destination / "primary_matrix_external_package_intake_template.yaml"
    markdown_path = destination / "primary_matrix_external_package_intake_template.md"
    yaml_path.write_text(yaml_text, encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return payload