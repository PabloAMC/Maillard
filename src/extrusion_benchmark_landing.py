from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import yaml

from src import data_paths
from src import data_access
from src.doe_generator import build_extrusion_benchmark_protocol


ROOT = data_paths.REPO_ROOT


def _under_root(root: Path, default: Path) -> Path:
    """``default`` when ``root`` is the repository checkout; the same repo-relative file
    re-rooted under ``root`` otherwise (tests and scripts pass explicit roots)."""
    return default if root == data_paths.REPO_ROOT else root / data_paths.rel(default)


def _load_json(path: Path) -> Dict[str, Any]:
    # Strict since 2026-09-01 (used to return {} for a missing registry).
    return data_access.load_json(path)


def _find_row(rows: Iterable[Mapping[str, Any]], key: str, value: str) -> Dict[str, Any]:
    for row in rows:
        if str(row.get(key, "")) == value:
            return dict(row)
    return {}


def _required_tradeoff_panel(contract: Mapping[str, Any]) -> Dict[str, Any]:
    required_panel = dict(contract.get("required_panel", {}))
    desirable = [str(item) for item in required_panel.get("desirable", [])]
    adverse = [str(item) for item in required_panel.get("adverse", [])]
    safety = [str(item) for item in required_panel.get("safety", [])]
    return {
        "desirable_targets": desirable,
        "adverse_targets": adverse,
        "safety_targets": safety,
        "same_run_minimum": [
            "2-methyl-3-furanthiol",
            "2-furfurylthiol",
            "bis(2-methyl-3-furyl) disulfide",
            "2,5-dimethylpyrazine",
            "Hexanal",
            "2-pentylfuran",
            "furfural",
            "acrylamide",
        ],
    }


def build_extrusion_external_closure_package(root: Path = ROOT) -> Dict[str, Any]:
    closure_contract = _load_json(_under_root(root, data_paths.EXTRUSION_EXTERNAL_CLOSURE_CONTRACT))
    primary_contract = _load_json(_under_root(root, data_paths.PPI_SPI_PRIMARY_BENCHMARK_CONTRACT))
    protocol = build_extrusion_benchmark_protocol(root)
    selected_matrix = str(protocol.get("selected_protein_type", "soy_iso"))
    matrix_row = _find_row(closure_contract.get("matrices", []), "matrix", selected_matrix)
    tradeoff_panel = _required_tradeoff_panel(primary_contract)
    required_metadata = list(primary_contract.get("analytical_requirements", {}).get("required_metadata", []))
    required_metadata.extend(
        [
            "extruder model and screw configuration",
            "screw speed rpm",
            "feed rate kg per h",
            "barrel temperature profile C",
            "die exit temperature C",
            "mean residence time seconds or equivalent",
            "extrusion moisture wt pct",
        ]
    )
    supportive_ids = [str(item) for item in matrix_row.get("supportive_anchor_ids", [])]
    contextual_ids = [str(item) for item in matrix_row.get("contextual_anchor_ids", [])]
    process_arms = []
    for arm in protocol.get("process_arms", []):
        process_arms.append(
            {
                "arm_id": str(arm.get("arm_id", "unknown")),
                "sme_kj_per_kg": float(arm.get("sme_kj_per_kg", 0.0) or 0.0),
                "barrel_temperature_celsius": float(arm.get("barrel_temperature_celsius", 0.0) or 0.0),
                "effective_temperature_celsius": float(arm.get("effective_temperature_celsius", 0.0) or 0.0),
                "moisture_regime": str(arm.get("moisture_regime", "unknown")),
                "water_activity": float(arm.get("water_activity", 0.0) or 0.0),
                "damage_predictions": dict(arm.get("predicted_damage_load", {})),
            }
        )

    return {
        "summary": {
            "contract_id": str(closure_contract.get("protocol_id", "extrusion_external_closure_2026")),
            "status": "specified_not_yet_measured",
            "selected_matrix": selected_matrix,
            "process_state": str(closure_contract.get("process_state", "extrusion_structured")),
            "arm_count": len(process_arms),
            "next_best_action": "run_spi_extrusion_arms_and_land_direct_damage_plus_tradeoff_panel_as_external_package",
        },
        "promotion_requirements": list(closure_contract.get("promotion_requirements", [])),
        "selected_matrix_contract": matrix_row,
        "same_run_tradeoff_panel": tradeoff_panel,
        "required_metadata": required_metadata,
        "process_arms": process_arms,
        "supportive_anchor_ids": supportive_ids,
        "contextual_anchor_ids": contextual_ids,
    }


def render_extrusion_external_closure_package_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    matrix_row = dict(payload.get("selected_matrix_contract", {}))
    tradeoff_panel = dict(payload.get("same_run_tradeoff_panel", {}))
    lines = [
        "# Extrusion External Closure Package",
        "",
        f"Contract id: {summary.get('contract_id', 'unknown')}",
        f"Status: {summary.get('status', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Process state: {summary.get('process_state', 'unknown')}",
        f"Next best action: {summary.get('next_best_action', 'unknown')}",
        "",
        "## Direct Marker Panel",
        "",
        "| Marker | Closure Role | Preferred Assays | Required | Repo Direct Anchor Available |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in matrix_row.get("direct_marker_panel", []):
        lines.append(
            f"| {row.get('display_name', 'unknown')} | {row.get('closure_role', 'unknown')} | {', '.join(row.get('preferred_assays', [])) or 'none'} | {bool(row.get('required', False))} | {bool(row.get('repo_direct_anchor_available', False))} |"
        )
    lines.extend(
        [
            "",
            "## Same-Run Tradeoff Panel",
            "",
            f"Desirable targets: {', '.join(tradeoff_panel.get('desirable_targets', [])) or 'none'}",
            f"Adverse targets: {', '.join(tradeoff_panel.get('adverse_targets', [])) or 'none'}",
            f"Safety targets: {', '.join(tradeoff_panel.get('safety_targets', [])) or 'none'}",
            f"Minimum same-run panel: {', '.join(tradeoff_panel.get('same_run_minimum', [])) or 'none'}",
            "",
            "## Process Arms",
            "",
            "| Arm | SME kJ/kg | Barrel Temp C | Effective Temp C | Moisture Regime | Water Activity |",
            "| --- | ---: | ---: | ---: | --- | ---: |",
        ]
    )
    for arm in payload.get("process_arms", []):
        lines.append(
            f"| {arm.get('arm_id', 'unknown')} | {float(arm.get('sme_kj_per_kg', 0.0)):.0f} | {float(arm.get('barrel_temperature_celsius', 0.0)):.1f} | {float(arm.get('effective_temperature_celsius', 0.0)):.1f} | {arm.get('moisture_regime', 'unknown')} | {float(arm.get('water_activity', 0.0)):.2f} |"
        )
    lines.extend(["", "## Required Metadata", ""])
    for item in payload.get("required_metadata", []):
        lines.append(f"- {item}")
    lines.extend(["", "## Context Anchors", ""])
    lines.append(f"Supportive anchors: {', '.join(payload.get('supportive_anchor_ids', [])) or 'none'}")
    lines.append(f"Contextual anchors: {', '.join(payload.get('contextual_anchor_ids', [])) or 'none'}")
    return "\n".join(lines) + "\n"


def build_extrusion_external_closure_workbook(root: Path = ROOT) -> Dict[str, Any]:
    package = build_extrusion_external_closure_package(root)
    protocol = build_extrusion_benchmark_protocol(root)
    hme_anchor = dict(protocol.get("closest_repo_backed_hme_anchor", {}))
    aligned_benchmark_id = "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    experiments = []
    for arm in package.get("process_arms", []):
        experiments.append(
            {
                "experiment_id": f"{package['summary']['selected_matrix']}_{arm['arm_id']}_external_closure_template_2026",
                "source_kind": "external_literature",
                "process_state": package["summary"]["process_state"],
                "protein_type": package["summary"]["selected_matrix"],
                "arm_id": arm["arm_id"],
                "benchmark_alignment": {
                    "target_benchmark_id": aligned_benchmark_id,
                    "notes": "Use the existing measured soy protocol-pilot benchmark as the closest aqueous baseline comparator for support-delta review.",
                },
                "target_process_settings": {
                    "sme_kj_per_kg": arm["sme_kj_per_kg"],
                    "barrel_temperature_celsius": arm["barrel_temperature_celsius"],
                    "effective_temperature_celsius": arm["effective_temperature_celsius"],
                    "extrusion_moisture_wt_pct": float(hme_anchor.get("extrusion_moisture_wt_pct", 0.0) or 0.0),
                    "screw_speed_rpm": float(hme_anchor.get("screw_speed_rpm", 0.0) or 0.0),
                    "feed_rate_kg_per_h": float(hme_anchor.get("feed_rate_kg_per_h", 0.0) or 0.0),
                    "die_exit_temp_C": float(hme_anchor.get("die_exit_temp_C", 0.0) or 0.0),
                },
                "measured_damage_markers": {
                    "reactive_lysine_fraction": "REPLACE_WITH_MEASURED_REACTIVE_LYSINE_FRACTION",
                    "furosine_mg_per_kg": "REPLACE_WITH_MEASURED_FUROSINE_MG_PER_KG",
                    "lysinoalanine_mg_per_kg": "REPLACE_WITH_MEASURED_LYSINOALANINE_MG_PER_KG",
                },
                "measured_volatiles_ppb": {
                    name: f"REPLACE_WITH_MEASURED_{name.upper().replace('-', '_').replace(',', '').replace(' ', '_')}_PPB"
                    for name in package.get("same_run_tradeoff_panel", {}).get("same_run_minimum", [])
                },
                "required_metadata": list(package.get("required_metadata", [])),
                "provenance": {
                    "origin": "REPLACE_WITH_CONTRACT_LAB_OR_PUBLICATION",
                    "source_reference": "REPLACE_WITH_TRACEABLE_CITATION_OR_REPORT",
                    "source_doi": "REPLACE_WITH_DOI_IF_AVAILABLE",
                    "measurement_date": "REPLACE_WITH_YYYY_MM_DD",
                    "notes": "Keep every numeric value traceable. Do not infer missing direct-damage markers.",
                },
            }
        )
    return {
        "package_id": "spi_extrusion_external_closure_workbook_2026",
        "status": "specified_not_yet_measured",
        "source_contract_id": package["summary"]["contract_id"],
        "selected_matrix": package["summary"]["selected_matrix"],
        "experiments": experiments,
    }


def render_extrusion_external_closure_workbook_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Extrusion External Closure Workbook",
        "",
        f"Package id: {payload.get('package_id', 'unknown')}",
        f"Status: {payload.get('status', 'unknown')}",
        f"Selected matrix: {payload.get('selected_matrix', 'unknown')}",
        "",
        "## Experiments",
        "",
    ]
    for experiment in payload.get("experiments", []):
        lines.extend(
            [
                f"- {experiment.get('experiment_id', 'unknown')}: arm {experiment.get('arm_id', 'unknown')}, SME {float(experiment.get('target_process_settings', {}).get('sme_kj_per_kg', 0.0)):.0f} kJ/kg",
                f"  Requires direct markers: {', '.join(experiment.get('measured_damage_markers', {}).keys())}",
                f"  Requires same-run volatiles: {', '.join(experiment.get('measured_volatiles_ppb', {}).keys())}",
            ]
        )
    return "\n".join(lines) + "\n"


def _prior_ids() -> Dict[str, str]:
    return {
        "disulfide_context": "raman_sds_extrusion_disulfide_severity",
        "mft_disulfide_trapping": "acs_jafc_3c02618_mft_disulfide_trapping_v1",
        "protein_binding_hierarchy": "acs_jafc_0c01925_protein_binding_hierarchy_v1",
    }


def build_extrusion_disulfide_follow_on_package(root: Path = ROOT) -> Dict[str, Any]:
    protocol = build_extrusion_benchmark_protocol(root)
    primary_contract = _load_json(_under_root(root, data_paths.PPI_SPI_PRIMARY_BENCHMARK_CONTRACT))
    process_payload = _load_json(_under_root(root, data_paths.PROCESS_STATE_CALIBRATIONS))
    prior_payload = _load_json(_under_root(root, data_paths.COMPUTATIONAL_PRIORS))
    ids = _prior_ids()
    process_entry = _find_row(process_payload.get("entries", []), "id", ids["disulfide_context"])
    prior_rows = list(prior_payload.get("retention_binding_priors", []))
    trapping_prior = _find_row(prior_rows, "id", ids["mft_disulfide_trapping"])
    binding_prior = _find_row(prior_rows, "id", ids["protein_binding_hierarchy"])
    same_run_volatiles = [
        "2-methyl-3-furanthiol",
        "2-furfurylthiol",
        "bis(2-methyl-3-furyl) disulfide",
        "2,5-dimethylpyrazine",
        "Hexanal",
        "2-pentylfuran",
    ]
    return {
        "summary": {
            "follow_on_id": "extrusion_5_8_disulfide_retention_follow_on_2026",
            "status": "trigger_ready_pending_first_measured_spi_extrusion_panel",
            "selected_matrix": str(protocol.get("selected_protein_type", "soy_iso")),
            "trigger_protocol_id": str(protocol.get("protocol_id", "spi_extrusion_mvp_benchmark_2026")),
            "arm_ids": [str(arm.get("arm_id", "unknown")) for arm in protocol.get("process_arms", [])],
            "objective": "fit free_SH_to_disulfide_severity_against_sulfur_retention_without_claiming_absolute_coefficients_before_the_first_measured_panel_lands",
        },
        "supporting_runtime_sources": {
            "process_state_calibration_id": str(process_entry.get("id", ids["disulfide_context"])),
            "retention_binding_prior_ids": [
                str(trapping_prior.get("id", ids["mft_disulfide_trapping"])),
                str(binding_prior.get("id", ids["protein_binding_hierarchy"])),
            ],
        },
        "required_same_run_observables": {
            "process_state_assays": [
                "Ellman free SH",
                "OPA free amino groups",
                "DSC or equivalent denaturation proxy",
                "furosine",
                "post-extrusion pH",
            ],
            "volatile_panel": same_run_volatiles,
            "feed_reference_assays": [
                "pre-extrusion Ellman free SH",
                "pre-extrusion OPA free amino groups",
            ],
        },
        "derived_metrics": [
            {
                "metric_id": "free_sh_retention_fraction",
                "definition": "post_extrusion_free_sh / pre_extrusion_free_sh",
                "reason": "Tracks how much sulfur accessibility survives extrusion severity.",
            },
            {
                "metric_id": "disulfide_pressure_proxy",
                "definition": "1 - free_sh_retention_fraction",
                "reason": "Bounded severity proxy aligned with the runtime disulfide context rather than an absolute fitted coefficient.",
            },
            {
                "metric_id": "sulfur_to_pyrazine_retention_ratio",
                "definition": "mft_recovery_fraction / pyrazine_recovery_fraction",
                "reason": "Separates sulfur trapping from general volatile release or stripping.",
            },
            {
                "metric_id": "furyl_disulfide_to_mft_ratio",
                "definition": "bis_2_methyl_3_furyl_disulfide_ppb / max(mft_ppb, epsilon)",
                "reason": "Connects sulfur oxidation/trapping pressure to the same arm's sulfur-positive yield.",
            },
            {
                "metric_id": "furosine_damage_gradient",
                "definition": "furosine_high_sme - furosine_low_sme",
                "reason": "Anchors the retention follow-on against protein damage severity, not just volatile movement.",
            },
        ],
        "runtime_context": {
            "disulfide_context_support": list(process_entry.get("what_it_supports", [])),
            "trapping_prior_notes": str(trapping_prior.get("notes", "")),
            "binding_prior_notes": str(binding_prior.get("notes", "")),
            "companion_tradeoff_panel": list(_required_tradeoff_panel(primary_contract).get("same_run_minimum", [])),
        },
    }


def render_extrusion_disulfide_follow_on_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    observables = dict(payload.get("required_same_run_observables", {}))
    lines = [
        "# Extrusion 5.8 Follow-On Package",
        "",
        f"Follow-on id: {summary.get('follow_on_id', 'unknown')}",
        f"Status: {summary.get('status', 'unknown')}",
        f"Selected matrix: {summary.get('selected_matrix', 'unknown')}",
        f"Trigger protocol: {summary.get('trigger_protocol_id', 'unknown')}",
        f"Arms: {', '.join(summary.get('arm_ids', [])) or 'none'}",
        "",
        "## Required Same-Run Observables",
        "",
        f"Process-state assays: {', '.join(observables.get('process_state_assays', [])) or 'none'}",
        f"Volatile panel: {', '.join(observables.get('volatile_panel', [])) or 'none'}",
        f"Feed reference assays: {', '.join(observables.get('feed_reference_assays', [])) or 'none'}",
        "",
        "## Derived Metrics",
        "",
    ]
    for row in payload.get("derived_metrics", []):
        lines.append(f"- {row.get('metric_id', 'unknown')}: {row.get('definition', '')}")
    lines.extend(["", "## Runtime Sources", ""])
    sources = dict(payload.get("supporting_runtime_sources", {}))
    lines.append(f"Process-state calibration: {sources.get('process_state_calibration_id', 'unknown')}")
    lines.append(f"Retention priors: {', '.join(sources.get('retention_binding_prior_ids', [])) or 'none'}")
    return "\n".join(lines) + "\n"


def build_extrusion_disulfide_follow_on_workbook(root: Path = ROOT) -> Dict[str, Any]:
    follow_on = build_extrusion_disulfide_follow_on_package(root)
    protocol = build_extrusion_benchmark_protocol(root)
    reference_benchmark_id = "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    experiments = []
    for arm in protocol.get("process_arms", []):
        experiments.append(
            {
                "experiment_id": f"{protocol.get('selected_protein_type', 'soy_iso')}_{arm.get('arm_id', 'unknown')}_5_8_follow_on_template_2026",
                "arm_id": str(arm.get("arm_id", "unknown")),
                "reference_benchmark_id": reference_benchmark_id,
                "target_process_settings": {
                    "sme_kj_per_kg": float(arm.get("sme_kj_per_kg", 0.0) or 0.0),
                    "effective_temperature_celsius": float(arm.get("effective_temperature_celsius", 0.0) or 0.0),
                },
                "feed_reference_assays": {
                    "pre_extrusion_free_sh_umol_per_g": "REPLACE_WITH_FEED_FREE_SH_UMOL_PER_G",
                    "pre_extrusion_free_amino_groups_umol_per_g": "REPLACE_WITH_FEED_FREE_AMINO_UMOL_PER_G",
                },
                "post_extrusion_process_state_assays": {
                    "ellman_free_sh_umol_per_g": "REPLACE_WITH_POST_EXTRUSION_FREE_SH_UMOL_PER_G",
                    "opa_free_amino_groups_umol_per_g": "REPLACE_WITH_POST_EXTRUSION_FREE_AMINO_UMOL_PER_G",
                    "furosine_mg_per_kg": "REPLACE_WITH_POST_EXTRUSION_FUROSINE_MG_PER_KG",
                    "post_extrusion_ph": "REPLACE_WITH_POST_EXTRUSION_PH",
                },
                "measured_volatiles_ppb": {
                    name: f"REPLACE_WITH_MEASURED_{name.upper().replace('-', '_').replace(',', '').replace(' ', '_')}_PPB"
                    for name in follow_on.get("required_same_run_observables", {}).get("volatile_panel", [])
                },
            }
        )
    return {
        "package_id": "spi_extrusion_5_8_follow_on_workbook_2026",
        "status": "trigger_ready_pending_first_measured_spi_extrusion_panel",
        "trigger_protocol_id": follow_on.get("summary", {}).get("trigger_protocol_id", "spi_extrusion_mvp_benchmark_2026"),
        "selected_matrix": follow_on.get("summary", {}).get("selected_matrix", "soy_iso"),
        "derived_metrics": list(follow_on.get("derived_metrics", [])),
        "experiments": experiments,
    }


def render_extrusion_disulfide_follow_on_workbook_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Extrusion 5.8 Follow-On Workbook",
        "",
        f"Package id: {payload.get('package_id', 'unknown')}",
        f"Status: {payload.get('status', 'unknown')}",
        f"Selected matrix: {payload.get('selected_matrix', 'unknown')}",
        "",
        "## Experiments",
        "",
    ]
    for experiment in payload.get("experiments", []):
        lines.extend(
            [
                f"- {experiment.get('experiment_id', 'unknown')}: arm {experiment.get('arm_id', 'unknown')}",
                f"  Feed reference assays: {', '.join(experiment.get('feed_reference_assays', {}).keys())}",
                f"  Post-extrusion assays: {', '.join(experiment.get('post_extrusion_process_state_assays', {}).keys())}",
                f"  Volatiles: {', '.join(experiment.get('measured_volatiles_ppb', {}).keys())}",
            ]
        )
    return "\n".join(lines) + "\n"


def _write_yaml(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(yaml.safe_dump(dict(payload), sort_keys=False, allow_unicode=False), encoding="utf-8")


def export_extrusion_external_closure_package(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_external_closure_package(root)
    markdown = render_extrusion_external_closure_package_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "extrusion_external_closure_package.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (destination / "extrusion_external_closure_package.md").write_text(markdown, encoding="utf-8")
    return payload


def export_extrusion_external_closure_workbook(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_external_closure_workbook(root)
    markdown = render_extrusion_external_closure_workbook_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    _write_yaml(destination / "extrusion_external_closure_workbook.yaml", payload)
    (destination / "extrusion_external_closure_workbook.md").write_text(markdown, encoding="utf-8")
    return payload


def export_extrusion_disulfide_follow_on_package(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_disulfide_follow_on_package(root)
    markdown = render_extrusion_disulfide_follow_on_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "extrusion_disulfide_follow_on_package.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (destination / "extrusion_disulfide_follow_on_package.md").write_text(markdown, encoding="utf-8")
    return payload


def export_extrusion_disulfide_follow_on_workbook(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_disulfide_follow_on_workbook(root)
    markdown = render_extrusion_disulfide_follow_on_workbook_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    _write_yaml(destination / "extrusion_disulfide_follow_on_workbook.yaml", payload)
    (destination / "extrusion_disulfide_follow_on_workbook.md").write_text(markdown, encoding="utf-8")
    return payload