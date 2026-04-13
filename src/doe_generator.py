"""
doe_generator.py

S12.3: Model-Guided Active Learning (DoE Loop)
This module formalizes "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
When the pipeline identifies unresolvable literature gaps in the process registry, 
it generates a formal lab protocol optimized for calibration gain.
"""

import json
import os
import logging
from pathlib import Path
from typing import Any, Dict, Mapping

from src.extrusion import build_extrusion_process_profile

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]

PRIMARY_PROTOCOL_CONTRACT_PATH = ROOT / "data" / "protocols" / "ppi_spi_primary_benchmark_contract.json"
SAFETY_REFERENCE_PAYLOADS_PATH = ROOT / "data" / "lit" / "safety_reference_payloads.json"
BENCHMARK_INTAKE_REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

# Basic templated DoE protocols based on the gap type
DOE_TEMPLATES = {
    "missing_absolute_anchor": {
        "method": "SIDA GC-MS/O",
        "factors": ["Temperature (95C, 120C)", "Time (10m, 30m)"],
        "instructions": "Prepare standard aqueous target matrix. Add isotopically labeled internal standards (e.g. 13C-MFT, d3-methional) before extraction. Heat uniformly."
    },
    "blocking_benchmark_gap": {
        "method": "Multi-factorial Quantitative Headspace SIDA",
        "factors": ["Precursor Load (1x, 5x)", "Temperature (90C, 130C)", "Matrix (SPI, PPI)"],
        "instructions": "Standard PBMA formulation baseline. Use Safe+SPME extraction to capture both highly volatile (H2S) and semi-volatile (pyrazines) simultaneously with SIDA quantitation."
    },
    "missing_positive_flavor_anchor": {
        "method": "Targeted GC-MS (Furanone band)",
        "factors": ["Water Activity", "pH (5.5, 6.5)"],
        "instructions": "Focus resolution on HEMF and DMHF specifically. Ensure SPME fiber captures polar furanones effectively without early saturation by aldehydes."
    },
    "missing_kinetic_dataset": {
        "method": "Time-course LC-MS/MS",
        "factors": ["Time (0, 5, 10, 20, 60 min)", "Temp (80C, 100C, 120C)"],
        "instructions": "Sample dynamically during the heating phase to capture pseudo-first-order formation rates. Quench immediately in ice bath."
    },
    "missing_process_state_bundle": {
        "method": "Simultaneous Ellman/OPA and DSC",
        "factors": ["Heating time", "Moisture Content (Extrusion regime)"],
        "instructions": "Run Differential Scanning Calorimetry alongside Ellman's reagent for free thiols and OPA for primary amines on the exact same homogenized samples to eliminate lot-to-lot variance."
    }
}

def generate_active_learning_requests(gap_registry_path: str) -> dict:
    """Read the structural gaps and generate explicit DoE requests."""
    
    if not os.path.exists(gap_registry_path):
        logger.warning(f"Gap registry not found at {gap_registry_path}. Cannot generate DoE.")
        return {"active_learning_requests": []}
        
    with open(gap_registry_path, "r", encoding="utf-8") as f:
        gaps = json.load(f)
        
    requests = []
    
    for entry in gaps.get("entries", []):
        if entry.get("wet_lab_requirement") != "required":
            continue
            
        gap_type = entry.get("gap_type", "unknown")
        template = DOE_TEMPLATES.get(gap_type, {
            "method": "Standard GC-MS Profiling",
            "factors": ["Vary primary bottleneck"],
            "instructions": "Perform standard empirical screen to establish directional bounds."
        })
        
        request = {
            "request_id": f"doe_request_{entry['gap_id']}",
            "priority": "HIGH" if "blocking" in gap_type else "MEDIUM",
            "target_chemistry_family": entry.get("chemistry_family"),
            "observables": entry.get("observable_panel_tags", []),
            "justification": entry.get("why_literature_cannot_close_it"),
            "suggested_next_step": entry.get("cheapest_next_step"),
            "experimental_design": template
        }
        
        requests.append(request)
        
    return {"active_learning_requests": requests}


def export_doe_requests(gap_registry_path: str, output_path: str):
    """Generates and writes the DoE requests to a JSON artifact."""
    payload = generate_active_learning_requests(gap_registry_path)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        
    logger.info(f"Exported {len(payload['active_learning_requests'])} Active Learning DoE requests to {output_path}")
    return payload


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _resolve_acrylamide_reference(root: Path) -> Dict[str, Any]:
    payload = _load_json(root / "data" / "lit" / "safety_reference_payloads.json")
    for entry in payload.get("entries", []):
        if str(entry.get("id", "")) == "choi_2024_pea_acrylamide_asparagine":
            return dict(entry)
    return {}


def _resolve_hme_control_anchor(root: Path) -> Dict[str, Any]:
    payload = _load_json(root / "data" / "lit" / "benchmark_intake_registry.json")
    for entry in payload.get("eligible_references", []):
        if str(entry.get("id", "")) == "pmc_2026_hme_hexanal_baseline":
            return dict(entry)
    return {}


def build_extrusion_benchmark_protocol(root: Path = ROOT) -> Dict[str, Any]:
    contract = _load_json(root / "data" / "protocols" / "ppi_spi_primary_benchmark_contract.json")
    acrylamide_reference = _resolve_acrylamide_reference(root)
    hme_control_anchor = _resolve_hme_control_anchor(root)
    hme_key_values = dict(hme_control_anchor.get("key_values", {}))

    barrel_temperature_celsius = 145.0
    water_activity = 0.75
    moisture_regime = "hme"
    sme_levels = [120.0, 180.0]
    process_arms = []
    for sme in sme_levels:
        profile = build_extrusion_process_profile(
            base_temperature_celsius=barrel_temperature_celsius,
            water_activity=water_activity,
            protein_type="soy_iso",
            sme_kj_per_kg=sme,
            moisture_regime=moisture_regime,
        )
        process_arms.append(
            {
                "arm_id": f"spi_hme_{int(sme)}_kj_per_kg",
                "protein_type": "soy_iso",
                "sme_kj_per_kg": sme,
                "barrel_temperature_celsius": barrel_temperature_celsius,
                "moisture_regime": moisture_regime,
                "water_activity": water_activity,
                "effective_temperature_celsius": float(profile.get("effective_temperature_celsius", barrel_temperature_celsius)),
                "predicted_damage_load": dict(profile.get("total_damage_load", {})),
                "zone_profile": list(profile.get("zone_profile", [])),
            }
        )

    required_panel = {
        "core": ["2-methyl-3-furanthiol", "hexanal", "furosine"],
        "recommended_same_run_extensions": [
            "2-furfurylthiol",
            "2-pentylfuran",
            "2,5-dimethylpyrazine",
            "acrylamide",
        ],
    }
    return {
        "protocol_id": "spi_extrusion_mvp_benchmark_2026",
        "protocol_status": "review_ready_missing_locked_lab_specs",
        "selected_protein_type": "soy_iso",
        "selection_rationale": [
            "The existing primary benchmark contract already includes a high-temperature SPI arm with same-run safety coverage.",
            "Soy isolate has the strongest current repo-backed path for combining meaty sulfur targets, off-note markers, and safety instrumentation.",
            "Using SPI first keeps the extrusion benchmark aligned with the current acrylamide and furosine evidence surface rather than opening a second matrix-selection decision.",
            "The closest structured HME control anchor already in the repo supplies practical moisture, screw-speed, feed-rate, and barrel-profile defaults for the first lab translation.",
        ],
        "design_targets": {
            "barrel_temperature_celsius": barrel_temperature_celsius,
            "sme_levels_kj_per_kg": sme_levels,
            "moisture_regime": moisture_regime,
            "water_activity": water_activity,
            "required_panel": required_panel,
            "shared_precursors": dict(contract.get("precursor_addition", {})),
        },
        "process_arms": process_arms,
        "closest_repo_backed_hme_anchor": {
            "citation": str(hme_control_anchor.get("citation", "Li et al. (2026), Foods 15(5):912")),
            "matrix_family": str(hme_control_anchor.get("matrix_family", "spi_wheat_gluten_hme")),
            "dry_basis_composition": {
                "spi_fraction": float(hme_key_values.get("spi_dry_basis_fraction", 0.0) or 0.0),
                "wheat_gluten_fraction": float(hme_key_values.get("wheat_gluten_dry_basis_fraction", 0.0) or 0.0),
            },
            "extrusion_moisture_wt_pct": float(hme_key_values.get("extrusion_moisture_wt_pct", 0.0) or 0.0),
            "screw_speed_rpm": float(hme_key_values.get("screw_speed_rpm", 0.0) or 0.0),
            "feed_rate_kg_per_h": float(hme_key_values.get("feed_rate_kg_per_h", 0.0) or 0.0),
            "barrel_temp_profile_C": list(hme_key_values.get("barrel_temp_profile_C", [])),
            "die_exit_temp_C": float(hme_key_values.get("die_exit_temp_C", 0.0) or 0.0),
            "replicates": int(hme_key_values.get("replicates", 3) or 3),
            "control_off_note_anchors_ug_per_kg": {
                "hexanal": float(hme_key_values.get("hexanal_control_ug_per_kg", 0.0) or 0.0),
                "heptanal": float(hme_key_values.get("heptanal_control_ug_per_kg", 0.0) or 0.0),
                "nonanal": float(hme_key_values.get("nonanal_control_ug_per_kg", 0.0) or 0.0),
                "1-hexanol": float(hme_key_values.get("1-hexanol_control_ug_per_kg", 0.0) or 0.0),
                "2-pentylfuran": float(hme_key_values.get("2-pentylfuran_control_ug_per_kg", 0.0) or 0.0),
            },
            "translation_note": "This anchor is a mixed SPI:wheat-gluten HME control rather than a pure-SPI benchmark, so it narrows SOP defaults but does not close the extruder-model or isotopic-spike fields.",
        },
        "companion_assays": [
            "Ellman free SH",
            "OPA free amino groups",
            "DSC or equivalent denaturation proxy",
        ],
        "repo_backed_lab_specs": {
            "headspace_method": {
                "platform": "HS-SPME-GC-MS",
                "fiber": "DVB/CAR/PDMS",
                "equilibration_temperature_celsius": 40.0,
                "equilibration_time_minutes": 30.0,
                "extraction_time_minutes": 10.0,
                "thiol_internal_standard_identities": ["[2H2]-MFT", "[2H2]-FFT"],
                "aldehyde_internal_standard_identity": "hexanal-d12",
            },
            "safety_method": {
                "platform": str(acrylamide_reference.get("method", {}).get("instrument", "LC-MS/MS")),
                "internal_standard_identity": str(acrylamide_reference.get("method", {}).get("internal_standard", "13C3-acrylamide")),
                "replicates": int(acrylamide_reference.get("method", {}).get("replicates", 3) or 3),
            },
        },
        "missing_locked_lab_specs": [
            {
                "field": "extruder_model",
                "status": "missing_from_repo",
                "why_it_matters": "The protocol can specify barrel temperature and SME targets, but the workspace does not yet canonize a specific twin-screw platform or screw configuration.",
            },
            {
                "field": "thiol_internal_standard_concentrations",
                "status": "missing_from_repo",
                "why_it_matters": "The repo identifies [2H2]-MFT and [2H2]-FFT as the preferred thiol standards, but it does not lock the spike concentration for routine lab execution.",
            },
            {
                "field": "hexanal_internal_standard_concentration",
                "status": "missing_from_repo",
                "why_it_matters": "The repo names hexanal-d12, but it does not yet declare the exact working concentration for the shared headspace method.",
            },
            {
                "field": "acrylamide_internal_standard_concentration",
                "status": "missing_from_repo",
                "why_it_matters": "The repo names 13C3-acrylamide and the LC-MS/MS platform, but it does not pin the final spike concentration for the extrusion benchmark run.",
            },
        ],
        "notes": [
            "This protocol is suitable as a wet-lab request artifact today, but not yet as a fully locked SOP.",
            "The first extrusion campaign should keep feed rate and screw speed constant while varying SME through the lab's chosen screw profile or energy setting.",
            "The same-run Ellman and furosine measurements are intentional so item 5.8 can be revisited with real data rather than solver-only assumptions.",
        ],
    }


def render_extrusion_benchmark_protocol_markdown(payload: Mapping[str, Any]) -> str:
    design_targets = dict(payload.get("design_targets", {}))
    headspace = dict(payload.get("repo_backed_lab_specs", {}).get("headspace_method", {}))
    safety = dict(payload.get("repo_backed_lab_specs", {}).get("safety_method", {}))
    hme_anchor = dict(payload.get("closest_repo_backed_hme_anchor", {}))
    lines = [
        "# Extrusion Benchmark Protocol",
        "",
        f"Protocol status: {payload.get('protocol_status', 'unknown')}",
        f"Selected protein type: {payload.get('selected_protein_type', 'unknown')}",
        "",
        "## Why This Design",
        "",
    ]
    for reason in payload.get("selection_rationale", []):
        lines.append(f"- {reason}")

    lines.extend([
        "",
        "## Minimum Viable Design",
        "",
        f"Barrel temperature C: {float(design_targets.get('barrel_temperature_celsius', 0.0)):.1f}",
        f"SME levels kJ/kg: {', '.join(str(int(value)) for value in design_targets.get('sme_levels_kj_per_kg', [])) or 'none'}",
        f"Moisture regime: {design_targets.get('moisture_regime', 'unknown')}",
        f"Water activity: {float(design_targets.get('water_activity', 0.0)):.2f}",
        f"Shared precursors: {', '.join(f'{key}={value}' for key, value in design_targets.get('shared_precursors', {}).items()) or 'none'}",
        f"Core panel: {', '.join(design_targets.get('required_panel', {}).get('core', [])) or 'none'}",
        f"Recommended extensions: {', '.join(design_targets.get('required_panel', {}).get('recommended_same_run_extensions', [])) or 'none'}",
        "",
        "| Arm | SME kJ/kg | Effective Temp C | Furosine mg/kg | LAL mg/kg |",
        "| --- | ---: | ---: | ---: | ---: |",
    ])
    for row in payload.get("process_arms", []):
        damage = dict(row.get("predicted_damage_load", {}))
        lines.append(
            f"| {row.get('arm_id', 'unknown')} | {float(row.get('sme_kj_per_kg', 0.0)):.0f} | {float(row.get('effective_temperature_celsius', 0.0)):.1f} | {float(damage.get('furosine_mg_per_kg', 0.0)):.1f} | {float(damage.get('lal_mg_per_kg', 0.0)):.1f} |"
        )

    lines.extend([
        "",
        "## Repo-Backed Lab Specs",
        "",
        f"Headspace platform: {headspace.get('platform', 'unknown')}",
        f"Fiber: {headspace.get('fiber', 'unknown')}",
        f"Headspace timing: {float(headspace.get('equilibration_temperature_celsius', 0.0)):.0f} C equilibration, {float(headspace.get('equilibration_time_minutes', 0.0)):.0f} min equilibration, {float(headspace.get('extraction_time_minutes', 0.0)):.0f} min extraction",
        f"Thiol internal standards: {', '.join(headspace.get('thiol_internal_standard_identities', [])) or 'none'}",
        f"Aldehyde internal standard: {headspace.get('aldehyde_internal_standard_identity', 'unknown')}",
        f"Safety platform: {safety.get('platform', 'unknown')}",
        f"Safety internal standard: {safety.get('internal_standard_identity', 'unknown')}",
        "",
        "## Closest Repo-Backed HME Anchor",
        "",
        f"Citation: {hme_anchor.get('citation', 'unknown')}",
        f"Matrix family: {hme_anchor.get('matrix_family', 'unknown')}",
        f"Dry-basis composition: SPI={float(hme_anchor.get('dry_basis_composition', {}).get('spi_fraction', 0.0)):.2f}, wheat gluten={float(hme_anchor.get('dry_basis_composition', {}).get('wheat_gluten_fraction', 0.0)):.2f}",
        f"Control extrusion settings: {float(hme_anchor.get('extrusion_moisture_wt_pct', 0.0)):.0f}% moisture, {float(hme_anchor.get('screw_speed_rpm', 0.0)):.0f} rpm, {float(hme_anchor.get('feed_rate_kg_per_h', 0.0)):.1f} kg/h",
        f"Barrel profile C: {', '.join(str(int(value)) for value in hme_anchor.get('barrel_temp_profile_C', [])) or 'none'}",
        f"Control off-note anchors ug/kg: {', '.join(f'{name}={value:.2f}' for name, value in hme_anchor.get('control_off_note_anchors_ug_per_kg', {}).items()) or 'none'}",
        f"Translation note: {hme_anchor.get('translation_note', 'none')}",
        "",
        "## Missing Locked Lab Specs",
        "",
    ])
    for row in payload.get("missing_locked_lab_specs", []):
        lines.append(
            f"- {row.get('field', 'unknown')}: {row.get('why_it_matters', '')}"
        )

    lines.extend([
        "",
        "## Notes",
        "",
    ])
    for note in payload.get("notes", []):
        lines.append(f"- {note}")

    return "\n".join(lines) + "\n"


def build_extrusion_sop_lock_register(root: Path = ROOT) -> Dict[str, Any]:
    protocol = build_extrusion_benchmark_protocol(root)
    anchor = dict(protocol.get("closest_repo_backed_hme_anchor", {}))
    headspace = dict(protocol.get("repo_backed_lab_specs", {}).get("headspace_method", {}))
    safety = dict(protocol.get("repo_backed_lab_specs", {}).get("safety_method", {}))
    contract = _load_json(root / "data" / "protocols" / "ppi_spi_primary_benchmark_contract.json")
    return {
        "summary": {
            "status": "partially_locked_repo_backed",
            "full_sop_lock_available": False,
            "reason": "The repo canonizes the process envelope, platform type, and analytical identities, but not an exact extruder brand/model or final isotope spike concentrations.",
        },
        "locked_fields": {
            "extruder_platform": "twin_screw",
            "moisture_regime": str(protocol.get("design_targets", {}).get("moisture_regime", "hme")),
            "screw_speed_rpm": float(anchor.get("screw_speed_rpm", 0.0) or 0.0),
            "feed_rate_kg_per_h": float(anchor.get("feed_rate_kg_per_h", 0.0) or 0.0),
            "barrel_temp_profile_C": list(anchor.get("barrel_temp_profile_C", [])),
            "die_exit_temp_C": float(anchor.get("die_exit_temp_C", 0.0) or 0.0),
            "headspace_platform": str(headspace.get("platform", "HS-SPME-GC-MS")),
            "headspace_fiber": str(headspace.get("fiber", "DVB/CAR/PDMS")),
            "thiol_internal_standard_identities": list(headspace.get("thiol_internal_standard_identities", [])),
            "aldehyde_internal_standard_identity": str(headspace.get("aldehyde_internal_standard_identity", "hexanal-d12")),
            "safety_platform": str(safety.get("platform", "UHPLC Nexera plus TQMS8040")),
            "safety_internal_standard_identity": str(safety.get("internal_standard_identity", "13C3-acrylamide")),
            "thiol_quantification_route": str(contract.get("analytical_requirements", {}).get("thiol_quantification", "thiol-specific derivatization with an internal standard suitable for MFT and FFT")),
        },
        "unresolved_fields": list(protocol.get("missing_locked_lab_specs", [])),
    }


def render_extrusion_sop_lock_register_markdown(payload: Mapping[str, Any]) -> str:
    summary = dict(payload.get("summary", {}))
    locked = dict(payload.get("locked_fields", {}))
    lines = [
        "# Extrusion SOP Lock Register",
        "",
        f"Status: {summary.get('status', 'unknown')}",
        f"Full SOP lock available: {summary.get('full_sop_lock_available', False)}",
        f"Reason: {summary.get('reason', 'unknown')}",
        "",
        "## Locked Fields",
        "",
    ]
    for key, value in locked.items():
        if isinstance(value, list):
            rendered = ", ".join(str(item) for item in value) if value else "none"
        else:
            rendered = str(value)
        lines.append(f"- {key}: {rendered}")
    lines.extend([
        "",
        "## Unresolved Fields",
        "",
    ])
    for row in payload.get("unresolved_fields", []):
        lines.append(f"- {row.get('field', 'unknown')}: {row.get('why_it_matters', '')}")
    return "\n".join(lines) + "\n"


def export_extrusion_sop_lock_register(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_sop_lock_register(root)
    markdown = render_extrusion_sop_lock_register_markdown(payload)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "extrusion_sop_lock_register.json"
    markdown_path = destination / "extrusion_sop_lock_register.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return payload


def export_extrusion_benchmark_protocol(output_dir: str, *, root: Path = ROOT) -> Dict[str, Any]:
    payload = build_extrusion_benchmark_protocol(root)
    markdown = render_extrusion_benchmark_protocol_markdown(payload)

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "extrusion_benchmark_protocol.json"
    markdown_path = destination / "extrusion_benchmark_protocol.md"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    logger.info("Exported extrusion benchmark protocol to %s and %s", json_path, markdown_path)
    return payload
