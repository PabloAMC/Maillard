from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src.artifact_io import load_optional_json_mapping
from src.family_ingestion_plan import load_family_ingestion_plan
from src.literature_family_registry import (
    iter_computational_prior_entries as iter_family_computational_prior_entries,
    iter_flavor_reference_entries as iter_family_flavor_reference_entries,
    iter_retention_reference_entries as iter_family_retention_reference_entries,
)
from src.matrix_prior_registry import summarize_family_prior_bundle
from src.matrix_targets import get_compound_panel_entry, iter_target_panel_entries
from src.projection_metadata import ProjectionMetadataMap


ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"


def _load_json_payload(file_name: str) -> dict[str, Any]:
    return load_optional_json_mapping(DATA_LIT_DIR / file_name)


FLAVOR_REFERENCE_PAYLOADS = _load_json_payload("flavor_reference_payloads.json")
RETENTION_REFERENCE_PAYLOADS = _load_json_payload("retention_reference_payloads.json")
COMPUTATIONAL_PRIORS_PAYLOAD = _load_json_payload("computational_priors.json")
try:
    FAMILY_INGESTION_PLAN_PAYLOAD = load_family_ingestion_plan()
except FileNotFoundError:
    FAMILY_INGESTION_PLAN_PAYLOAD = {}


_FAMILY_ENTRIES_BY_SLR = {
    str(entry.get("slr_family", "")).zfill(2): dict(entry)
    for entry in FAMILY_INGESTION_PLAN_PAYLOAD.get("families", [])
    if isinstance(entry, Mapping)
}

_RUNTIME_LANE_PRIOR_REGISTRY = {
    "pyrazine_control": {
        "section_name": "pyrazine_control_priors",
        "default_id": "slr_pyrazine_control_surface_v1",
        "slr_family": "01",
    },
    "strecker_crosstalk": {
        "section_name": "strecker_crosstalk_priors",
        "default_id": "lincoln_2025_polyphenol_crosstalk_v1",
        "slr_family": "02",
    },
    "furanone_support": {
        "section_name": "furanone_priors",
        "default_id": "blank_fay_1996_furanone_expectation_v1",
        "slr_family": "01",
    },
    "thiamine_fragmentation": {
        "section_name": "thiamine_pathway_priors",
        "default_id": "cerny_2007_thiamine_split_v1",
        "slr_family": "03",
    },
}


from src.text_utils import normalize_name_spaced as _normalize_name


def _sigmoid(value: float, center: float, width: float) -> float:
    width = max(float(width), 1.0e-6)
    exponent = max(-60.0, min(60.0, -(float(value) - center) / width))
    return 1.0 / (1.0 + math.exp(exponent))


def _log_floor_decay(value: Optional[float], *, reference: float, slope: float, floor: float) -> float:
    if value is None or value <= 0.0:
        return 1.0
    factor = 1.0 - slope * math.log1p(float(value) / max(reference, 1.0e-6))
    return max(float(floor), min(1.0, float(factor)))


def _family_runtime_descriptor(slr_family: str) -> Dict[str, Any]:
    row = _FAMILY_ENTRIES_BY_SLR.get(str(slr_family).zfill(2), {})
    if not row:
        return {"slr_family": str(slr_family).zfill(2)}
    return {
        "slr_family": str(row.get("slr_family", slr_family)).zfill(2),
        "family_id": str(row.get("family_id", "unknown")),
        "display_name": str(row.get("display_name", "unknown")),
        "strategic_posture": str(row.get("strategic_posture", "unknown")),
        "runtime_concept": str(row.get("runtime_concept", "unknown")),
        "implementation_wave": int(row.get("implementation_wave", 99)),
        "order_in_wave": int(row.get("order_in_wave", 99)),
    }


def _family_sort_key(slr_family: str) -> tuple[int, int, str]:
    descriptor = _family_runtime_descriptor(slr_family)
    return (
        int(descriptor.get("implementation_wave", 99)),
        int(descriptor.get("order_in_wave", 99)),
        str(descriptor.get("slr_family", slr_family)),
    )


def _entry_matches_compound_name(entry: Mapping[str, Any], compound_name: str) -> bool:
    requested = _normalize_name(compound_name)
    candidate_values: List[str] = []
    for key in ["compound", "display_name", "canonical_name", "name"]:
        value = str(entry.get(key, "")).strip()
        if value:
            candidate_values.append(value)
    aliases = entry.get("aliases", []) or []
    if isinstance(aliases, list):
        candidate_values.extend(str(value).strip() for value in aliases if str(value).strip())
    for key in ["polyphenol_examples", "required_sugars", "target_compounds_or_state_variables"]:
        values = entry.get(key, []) or []
        if isinstance(values, list):
            candidate_values.extend(str(value).strip() for value in values if str(value).strip())
    return any(_normalize_name(value) == requested for value in candidate_values)


def _process_state_matches(entry: Mapping[str, Any], process_state: Optional[str]) -> bool:
    if process_state is None:
        return True
    requested = _normalize_name(process_state)
    values = entry.get("process_state_scope", entry.get("process_state_applicability", [])) or []
    if not isinstance(values, list) or not values:
        return True
    return requested in {_normalize_name(value) for value in values}


def _protein_type_matches_retention_matrix(matrix_family: str, protein_type: Optional[str]) -> bool:
    if protein_type is None:
        return True
    normalized_protein = _normalize_name(protein_type)
    normalized_matrix = _normalize_name(matrix_family)
    if not normalized_protein:
        return True
    if normalized_protein.startswith("pea"):
        return normalized_matrix.startswith("pea")
    if normalized_protein.startswith("soy"):
        return normalized_matrix.startswith("soy")
    if normalized_protein.startswith("myco"):
        return "myco" in normalized_matrix or "fung" in normalized_matrix
    return normalized_protein in normalized_matrix


def _resolve_flavor_reference_pipeline_role(entry: Mapping[str, Any]) -> str:
    explicit = str(entry.get("pipeline_role", "")).strip().lower()
    if explicit:
        return explicit
    benchmark_role = str(entry.get("benchmark_role", "")).strip().lower()
    return _DEFAULT_FLAVOR_PIPELINE_ROLE_BY_BENCHMARK.get(benchmark_role, "reference_only")


def query_family_runtime_priors(
    *,
    runtime_lane: Optional[str] = None,
    family: Optional[str] = None,
    protein_type: Optional[str] = None,
    process_state: Optional[str] = None,
    compound_name: Optional[str] = None,
    entry_id: Optional[str] = None,
) -> List[Dict[str, Any]]:
    lane_metadata = _RUNTIME_LANE_PRIOR_REGISTRY.get(str(runtime_lane or "").strip(), {}) if runtime_lane is not None else {}
    section_name = str(lane_metadata.get("section_name", "")).strip()
    query_family = family
    if query_family is None and lane_metadata.get("slr_family"):
        query_family = str(lane_metadata.get("slr_family", ""))

    rows: List[Dict[str, Any]] = []
    for entry in iter_family_computational_prior_entries(family=query_family, protein_type=protein_type):
        row = dict(entry)
        if section_name and str(row.get("_section_name", "")).strip() != section_name:
            continue
        if entry_id is not None and str(row.get("id", "")).strip() != str(entry_id).strip():
            continue
        if compound_name is not None and not _entry_matches_compound_name(row, compound_name):
            continue
        if not _process_state_matches(row, process_state):
            continue
        row["runtime_lane"] = str(runtime_lane or row.get("runtime_lane", "generic_prior_query"))
        row["prior_section"] = str(row.get("_section_name", section_name or "unknown"))
        if "family" not in row:
            descriptor = row.get("family_descriptor", {}) if isinstance(row.get("family_descriptor", {}), Mapping) else {}
            row["family"] = {
                "slr_family": str(descriptor.get("slr_family", row.get("slr_family_source", ""))).zfill(2),
                "family_id": str(descriptor.get("family_id", row.get("chemistry_family", "unknown"))),
                "display_name": str(descriptor.get("display_name", row.get("chemistry_family", "unknown"))),
                "strategic_posture": str(descriptor.get("strategic_posture", "unknown")),
            }
        rows.append(row)
    rows.sort(key=lambda row: (_family_sort_key(str(row.get("slr_family_source", row.get("family", {}).get("slr_family", "99")))), str(row.get("id", "unknown"))))
    return rows


def get_family_runtime_prior(
    *,
    runtime_lane: str,
    entry_id: Optional[str] = None,
) -> dict[str, Any]:
    lane_metadata = _RUNTIME_LANE_PRIOR_REGISTRY.get(str(runtime_lane).strip(), {})
    requested_id = str(entry_id or lane_metadata.get("default_id", "")).strip() or None
    rows = query_family_runtime_priors(runtime_lane=runtime_lane, entry_id=requested_id)
    return rows[0] if rows else {}


def get_pyrazine_control_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="pyrazine_control")


def get_furanone_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="furanone_support")


def get_thiamine_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="thiamine_fragmentation")


def get_strecker_crosstalk_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="strecker_crosstalk")


_DEFAULT_FLAVOR_PIPELINE_ROLE_BY_BENCHMARK = {
    "reference_anchor": "secondary_marker",
    "directional_comparison_anchor": "diagnostic_marker",
    "pbma_counterexample": "reference_only",
    "low_confidence_mechanistic_anchor": "optimization_constraint",
}


def query_flavor_reference_entries(
    *,
    family: Optional[str] = None,
    compound_name: Optional[str] = None,
    pipeline_role: Optional[str] = None,
    entry_id: Optional[str] = None,
    payload_section: Optional[str] = None,
    scoring_only: bool = False,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for entry in iter_family_flavor_reference_entries(family=family):
        row = dict(entry)
        row["section"] = str(row.get("_section_name", "unknown"))
        row["pipeline_role"] = _resolve_flavor_reference_pipeline_role(row)
        if entry_id is not None and str(row.get("id", "")).strip() != str(entry_id).strip():
            continue
        if payload_section is not None and str(row.get("section", "")).strip() != str(payload_section).strip():
            continue
        if compound_name is not None and not _entry_matches_compound_name(row, compound_name):
            continue
        if pipeline_role is not None and str(row.get("pipeline_role", "")).strip().lower() != str(pipeline_role).strip().lower():
            continue
        if scoring_only and not _role_supports_runtime_scoring(str(row.get("pipeline_role", ""))):
            continue
        rows.append(row)
    rows.sort(key=lambda row: (str(row.get("pipeline_role", "unknown")), str(row.get("compound", "unknown"))))
    return rows


def _iter_flavor_entries() -> Iterable[tuple[str, Dict[str, Any]]]:
    for entry in query_flavor_reference_entries():
        yield str(entry.get("section", "unknown")), dict(entry)


def _flavor_entry_by_id(entry_id: str) -> Dict[str, Any]:
    rows = query_flavor_reference_entries(entry_id=entry_id)
    return rows[0] if rows else {}


def get_flavor_reference_pipeline_role(entry_id: str) -> str:
    entry = _flavor_entry_by_id(entry_id)
    return _resolve_flavor_reference_pipeline_role(entry)


def build_flavor_reference_policy_summary() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for section_name, entry in _iter_flavor_entries():
        rows.append(
            {
                "id": str(entry.get("id", "unknown")),
                "section": section_name,
                "compound": str(entry.get("compound", "unknown")),
                "benchmark_role": str(entry.get("benchmark_role", "unknown")),
                "pipeline_role": str(entry.get("pipeline_role", "reference_only")),
                "target_direction": str(entry.get("target_direction", "")),
                "source_citation": str(entry.get("source_citation", "unknown")),
            }
        )
    rows.sort(key=lambda row: (row["pipeline_role"], row["compound"]))
    return rows


def _role_supports_runtime_scoring(role: str) -> bool:
    return str(role).strip().lower() in {
        "primary_target",
        "secondary_marker",
        "diagnostic_marker",
        "optimization_constraint",
    }


def query_retention_reference_entries(
    *,
    family: Optional[str] = None,
    compound_name: Optional[str] = None,
    protein_type: Optional[str] = None,
    process_state: Optional[str] = None,
    surrogate_type: Optional[str] = None,
    entry_id: Optional[str] = None,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for entry in iter_family_retention_reference_entries(family=family):
        row = dict(entry)
        row["section"] = str(row.get("_section_name", "unknown"))
        row["matrix_family"] = str(row.get("_matrix_family", "unknown"))
        surrogate = row.get("runtime_surrogate", {}) if isinstance(row.get("runtime_surrogate", {}), Mapping) else {}
        row["runtime_surrogate_type"] = str(surrogate.get("type", "")) if surrogate else ""
        if entry_id is not None and str(row.get("id", "")).strip() != str(entry_id).strip():
            continue
        if compound_name is not None and not _entry_matches_compound_name(row, compound_name):
            continue
        if protein_type is not None and not _protein_type_matches_retention_matrix(str(row.get("matrix_family", "")), protein_type):
            continue
        if not _process_state_matches(row, process_state):
            continue
        if surrogate_type is not None and str(row.get("runtime_surrogate_type", "")).strip() != str(surrogate_type).strip():
            continue
        rows.append(row)
    rows.sort(key=lambda row: (str(row.get("section", "unknown")), str(row.get("matrix_family", "unknown")), str(row.get("id", "unknown"))))
    return rows


def _iter_retention_entries() -> Iterable[tuple[str, str, Dict[str, Any]]]:
    for entry in query_retention_reference_entries():
        yield str(entry.get("section", "unknown")), str(entry.get("matrix_family", "unknown")), dict(entry)


def _retention_entry_by_id(entry_id: str) -> Dict[str, Any]:
    rows = query_retention_reference_entries(entry_id=entry_id)
    return rows[0] if rows else {}


def _retention_source_label(entry: Mapping[str, Any]) -> str:
    return str(entry.get("source_citation") or entry.get("id") or "unknown retention reference")


def _runtime_volatile_family(normalized_name: str) -> str:
    if any(token in normalized_name for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
        return "sulfur"
    if "pyrazine" in normalized_name:
        return "pyrazine"
    if any(token in normalized_name for token in ["furan", "furfural"]):
        return "furan"
    if any(token in normalized_name for token in ["ol", "alcohol"]):
        return "alcohol"
    if any(token in normalized_name for token in ["anal", "enal", "aldehyde"]):
        return "aldehyde"
    return "other"


def _build_extrusion_runtime_surrogate(
    compound_name: str,
    *,
    protein_type: Optional[str],
    temperature_celsius: float,
    time_minutes: Optional[float],
    water_activity: Optional[float],
    process_state: Optional[str],
) -> Dict[str, Any]:
    if process_state not in {"aqueous_pre_extrusion_model", "extrusion_structured"}:
        return {}

    protein = str(protein_type or "")
    if not protein.startswith(("pea", "soy", "myco")):
        return {}

    normalized = _normalize_name(compound_name)
    family = _runtime_volatile_family(normalized)
    aw = 0.55 if water_activity is None else max(0.05, min(1.0, float(water_activity)))
    dryness = max(0.0, min(1.0, (0.65 - aw) / 0.35))
    thermal_drive = _sigmoid(float(temperature_celsius), 150.0, 9.0)
    residence_factor = 1.0 if time_minutes is None else max(0.88, min(1.15, 1.0 + 0.06 * math.log1p(max(float(time_minutes), 0.0) / 3.0)))

    if family in {"aldehyde", "alcohol"}:
        moisture_factor = max(0.42, 1.0 - 0.34 * dryness)
        structure_factor = max(0.55, 1.0 - 0.18 * thermal_drive)
    elif family == "furan":
        moisture_factor = max(0.50, 1.0 - 0.22 * dryness)
        structure_factor = max(0.62, 1.0 - 0.12 * thermal_drive)
    elif family in {"sulfur", "pyrazine"}:
        moisture_factor = 1.0 + 0.12 * (1.0 - dryness)
        structure_factor = max(0.68, 1.0 - 0.10 * dryness * thermal_drive)
    else:
        moisture_factor = max(0.55, 1.0 - 0.18 * dryness)
        structure_factor = max(0.70, 1.0 - 0.08 * thermal_drive)

    if process_state == "aqueous_pre_extrusion_model":
        state_factor = 1.0 if family in {"sulfur", "pyrazine"} else 0.92
        retention_mode = "extrusion_moisture_redistribution"
    else:
        state_factor = 0.82 if family in {"sulfur", "pyrazine"} else 0.72
        retention_mode = "extrusion_structured_entrapment"

    surrogate_factor = max(0.18, min(1.35, moisture_factor * structure_factor * state_factor * residence_factor))
    return {
        "retention_runtime_mode": retention_mode,
        "retention_reference_sources": [
            "Extrusion severity surrogate derived from shared matrix accessibility/process-state assumptions",
        ],
        "extrusion_moisture_factor": float(moisture_factor),
        "extrusion_structure_factor": float(structure_factor * state_factor),
        "dynamic_retention_factor": float(surrogate_factor),
    }


def get_retention_ph_release_profile(
    compound_name: str,
    *,
    protein_type: Optional[str],
) -> Dict[str, Any]:
    for entry in query_retention_reference_entries(
        compound_name=compound_name,
        protein_type=protein_type,
        surrogate_type="ph_release_modifier",
    ):
        surrogate = entry.get("runtime_surrogate", {})
        if isinstance(surrogate, Mapping):
            return {
                **dict(surrogate),
                "source": _retention_source_label(entry),
                "entry_id": str(entry.get("id", "unknown")),
            }
    return {}


def describe_retention_runtime(
    compound_name: str,
    *,
    protein_type: Optional[str],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
    water_activity: Optional[float] = None,
    process_state: Optional[str],
) -> Dict[str, Any]:
    normalized = _normalize_name(compound_name)
    protein = str(protein_type or "")
    temperature = float(temperature_celsius if temperature_celsius is not None else 25.0)
    duration = None if time_minutes is None else float(time_minutes)

    release_factor = 1.0
    temporal_attenuation_factor = 1.0
    dynamic_retention_factor = 1.0
    retention_mode = "static_class_profile"
    sources: List[str] = []
    extrusion_moisture_factor = 1.0
    extrusion_structure_factor = 1.0

    if protein.startswith("pea") and normalized == "hexanal":
        retention_mode = "direct_binding_plus_ph_release_reference"
        sources.extend(
            filter(
                None,
                [
                    _retention_source_label(_retention_entry_by_id("jafc_3c05991_ppi_hexanal_binding")),
                    _retention_source_label(_retention_entry_by_id("karolkowski_2021_ppi_hexanal_ph_release")),
                ],
            )
        )
    elif protein.startswith("pea") and normalized in {"2 pentylfuran", "2 pentyl furan"}:
        retention_mode = "ph_release_reference"
        sources.append(_retention_source_label(_retention_entry_by_id("karolkowski_2021_ppi_2_pentylfuran_native_panel")))

    if protein.startswith("soy") and normalized == "hexanal":
        xu_hexanal = _retention_entry_by_id("xu_2023_spi_hexanal_temporal_profile")
        xu_surrogate = xu_hexanal.get("runtime_surrogate", {}) if isinstance(xu_hexanal, Mapping) else {}
        temperature_gate = float(xu_surrogate.get("temperature_gate_celsius", 90.0) or 90.0)
        reference_minutes = float(xu_surrogate.get("reference_minutes", 8.0) or 8.0)
        log_slope = float(xu_surrogate.get("log_slope", 0.18) or 0.18)
        floor = float(xu_surrogate.get("floor", 0.62) or 0.62)
        release_factor = 1.0 + 0.55 * _sigmoid(temperature, 58.0, 10.0)
        temporal_attenuation_factor = _log_floor_decay(duration, reference=reference_minutes, slope=log_slope, floor=floor) if temperature >= temperature_gate else 1.0
        dynamic_retention_factor = max(0.85, min(1.65, release_factor * temporal_attenuation_factor))
        retention_mode = "reversible_release_plus_temporal_attenuation"
        sources.extend([
            _retention_source_label(_retention_entry_by_id("ince_2024_glycinin_hexanal_binding")),
            _retention_source_label(_retention_entry_by_id("jafc_3c05991_ppi_hexanal_binding")),
        ])
        if temperature >= temperature_gate:
            sources.append(_retention_source_label(xu_hexanal))
    elif protein.startswith("soy") and normalized in {"2 pentylfuran", "2 pentyl furan"}:
        xu_furan = _retention_entry_by_id("xu_2023_spi_2_pentylfuran_temporal_profile")
        xu_surrogate = xu_furan.get("runtime_surrogate", {}) if isinstance(xu_furan, Mapping) else {}
        temperature_gate = float(xu_surrogate.get("temperature_gate_celsius", 90.0) or 90.0)
        reference_minutes = float(xu_surrogate.get("reference_minutes", 8.0) or 8.0)
        log_slope = float(xu_surrogate.get("log_slope", 0.24) or 0.24)
        floor = float(xu_surrogate.get("floor", 0.55) or 0.55)
        temporal_attenuation_factor = _log_floor_decay(duration, reference=reference_minutes, slope=log_slope, floor=floor) if temperature >= temperature_gate else 1.0
        dynamic_retention_factor = temporal_attenuation_factor
        retention_mode = "temporal_attenuation"
        sources.extend([
            _retention_source_label(xu_furan),
            _retention_source_label(_retention_entry_by_id("shu_2024_heated_spi_2_pentylfuran_censored")),
        ])
    elif protein.startswith("soy") and any(token in normalized for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
        retention_mode = "sulfur_proxy_reference"
        sources.append(_retention_source_label(_retention_entry_by_id("zhang_2026_spi_lenthionine_retention")))

    if not sources and process_state == "heated_matrix" and protein.startswith("soy"):
        sources.append("heated soy process-state carryover")

    extrusion_surrogate = _build_extrusion_runtime_surrogate(
        compound_name,
        protein_type=protein_type,
        temperature_celsius=temperature,
        time_minutes=duration,
        water_activity=water_activity,
        process_state=process_state,
    )
    if extrusion_surrogate:
        extrusion_moisture_factor = float(extrusion_surrogate.get("extrusion_moisture_factor", 1.0))
        extrusion_structure_factor = float(extrusion_surrogate.get("extrusion_structure_factor", 1.0))
        dynamic_retention_factor *= float(extrusion_surrogate.get("dynamic_retention_factor", 1.0))
        dynamic_retention_factor = max(0.05, min(1.65, dynamic_retention_factor))
        extrusion_sources = extrusion_surrogate.get("retention_reference_sources", [])
        if extrusion_sources:
            sources.extend(str(item) for item in extrusion_sources)
        extrusion_mode = str(extrusion_surrogate.get("retention_runtime_mode", "extrusion_surrogate"))
        retention_mode = extrusion_mode if retention_mode == "static_class_profile" else f"{retention_mode}+{extrusion_mode}"

    return {
        "retention_runtime_mode": retention_mode,
        "retention_reference_sources": sources,
        "reversible_release_factor": float(release_factor),
        "temporal_attenuation_factor": float(temporal_attenuation_factor),
        "dynamic_retention_factor": float(dynamic_retention_factor),
        "extrusion_moisture_factor": float(extrusion_moisture_factor),
        "extrusion_structure_factor": float(extrusion_structure_factor),
    }


def _projection_rows_to_signal_map(projection_metadata: ProjectionMetadataMap) -> Dict[str, float]:
    signals: Dict[str, float] = {}
    for row in projection_metadata.values():
        compound = str(row.get("compound", "")).strip()
        if not compound:
            continue
        normalized = _normalize_name(compound)
        signals[normalized] = max(signals.get(normalized, 0.0), float(row.get("observable_ppb", 0.0)))
    return signals


def _compound_signal(signal_map: Dict[str, float], names: Iterable[str]) -> float:
    return max(float(signal_map.get(_normalize_name(name), 0.0)) for name in names)


def _pyrazine_signal(signal_map: Dict[str, float]) -> float:
    total = 0.0
    seen: set[str] = set()
    for normalized_name, value in signal_map.items():
        if "pyrazine" not in normalized_name or normalized_name in seen:
            continue
        total += float(value)
        seen.add(normalized_name)
    return total


def _normalize_intervention_tokens(interventions: Optional[List[Any]]) -> List[str]:
    tokens: List[str] = []
    for intervention in interventions or []:
        if isinstance(intervention, Mapping):
            name = intervention.get("name")
            if name is not None:
                tokens.append(_normalize_name(str(name)))
            else:
                for key, value in intervention.items():
                    tokens.append(_normalize_name(str(key)))
                    if isinstance(value, str) and value.strip():
                        tokens.append(_normalize_name(value))
        else:
            tokens.append(_normalize_name(str(intervention)))
    return tokens


def _build_family_lane_payload(
    slr_family: str,
    *,
    active: bool,
    summary: str,
    metrics: Mapping[str, Any],
) -> Dict[str, Any]:
    return {
        **_family_runtime_descriptor(slr_family),
        "active": bool(active),
        "summary": str(summary),
        **dict(metrics),
    }


def _lipid_benchmark_target_panel() -> List[Dict[str, Any]]:
    benchmarkable_states = {
        "externally_benchmarked",
        "internally_benchmarked",
        "conditional_calibration",
    }
    rows = []
    for row in iter_target_panel_entries(
        target_class="adverse_lipid_markers",
        chemistry_family="lipid_oxidation_and_carbonylic_crosstalk",
        observable_kind="volatile",
    ):
        if str(row.get("evidence_state", "still_missing")) not in benchmarkable_states:
            continue
        rows.append(dict(row))
    return rows


def _build_lipid_crosstalk_lane(
    *,
    signal_map: Dict[str, float],
    normalized_sugars: List[str],
    lipids: List[str],
    polyphenol_active: bool,
    protein_type: Optional[str],
) -> Dict[str, Any]:
    lipid_marker_signals = {
        "hexanal": _compound_signal(signal_map, ["hexanal"]),
        "2_pentylfuran": _compound_signal(signal_map, ["2-pentylfuran", "2 pentylfuran"]),
        "1_octen_3_ol": _compound_signal(signal_map, ["1-octen-3-ol", "1 octen 3 ol"]),
        "nonanal": _compound_signal(signal_map, ["nonanal"]),
        "2_4_decadienal": _compound_signal(
            signal_map,
            ["2,4-decadienal", "2 4 decadienal", "e,e-2,4-decadienal", "e e 2 4 decadienal"],
        ),
    }
    lipid_marker_signal_ppb = sum(float(value) for value in lipid_marker_signals.values())
    marker_count = sum(1 for value in lipid_marker_signals.values() if float(value) > 0.0)
    lipid_input_count = len(lipids)
    donor_pressure = 1.0 if any("ribose" in sugar or "xylose" in sugar or "fructose" in sugar for sugar in normalized_sugars) else 0.55
    carbonyl_competition_factor = min(1.75, 0.18 * lipid_input_count + donor_pressure * lipid_marker_signal_ppb / 120.0)
    retention_rows = query_retention_reference_entries(
        family="02",
        protein_type=protein_type,
        process_state="heated_matrix",
    )
    crosstalk_prior_rows = query_family_runtime_priors(
        runtime_lane="strecker_crosstalk",
        family="02",
        process_state="heated_matrix",
    )
    benchmark_ready_targets = _lipid_benchmark_target_panel()
    crosstalk_prior_active = bool(crosstalk_prior_rows) and bool(polyphenol_active or lipid_input_count or lipid_marker_signal_ppb > 0.0)
    strecker_suppression_factor = min(
        1.0,
        (carbonyl_competition_factor / 1.75) * (1.15 if polyphenol_active and crosstalk_prior_active else 1.0),
    )
    maillard_closure_pressure = min(2.0, carbonyl_competition_factor * (1.0 + 0.25 * float(polyphenol_active and crosstalk_prior_active)))
    active = bool(lipid_input_count or lipid_marker_signal_ppb > 0.0)
    summary = "No explicit lipid-driven crosstalk lane active."
    if active:
        summary = (
            "Lipid-derived adverse markers are present or lipid precursors are active, so the lane should be split into adverse-marker retention and carbonyl-competition support rather than treated as a generic off-note penalty."
        )
        if polyphenol_active:
            summary = (
                "Lipid-derived adverse markers coexist with polyphenol chemistry, so oxidative retention and Strecker-suppression crosstalk should be treated as coupled sub-lanes."
            )
    return _build_family_lane_payload(
        "02",
        active=active,
        summary=summary,
        metrics={
            "lipid_marker_signals_ppb": lipid_marker_signals,
            "lipid_marker_signal_ppb": float(lipid_marker_signal_ppb),
            "lipid_marker_count": int(marker_count),
            "lipid_input_count": int(lipid_input_count),
            "polyphenol_crosstalk_active": bool(polyphenol_active and active),
            "carbonyl_competition_factor": float(carbonyl_competition_factor),
            "benchmark_ready_targets": [str(row.get("display_name", row.get("canonical_name", "unknown"))) for row in benchmark_ready_targets],
            "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
            "competition_prior_ids": [str(row.get("id", "unknown")) for row in crosstalk_prior_rows],
            "strecker_suppression_factor": float(strecker_suppression_factor),
            "maillard_closure_pressure": float(maillard_closure_pressure),
            "runtime_sub_lanes": {
                "adverse_marker_generation_and_retention": {
                    "active": bool(active),
                    "benchmark_ready_target_count": len(benchmark_ready_targets),
                    "benchmark_ready_targets": [str(row.get("display_name", row.get("canonical_name", "unknown"))) for row in benchmark_ready_targets],
                    "retention_reference_count": len(retention_rows),
                    "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
                },
                "carbonyl_competition_and_crosstalk": {
                    "active": bool(crosstalk_prior_active),
                    "competition_prior_count": len(crosstalk_prior_rows),
                    "competition_prior_ids": [str(row.get("id", "unknown")) for row in crosstalk_prior_rows],
                    "strecker_suppression_factor": float(strecker_suppression_factor),
                    "maillard_closure_pressure": float(maillard_closure_pressure),
                },
            },
        },
    )


def _build_state_marker_payload(
    marker_name: str,
    *,
    active: bool,
    state_value: Any,
    state_value_summary: str,
    influence_mode: str,
    family_lane: Mapping[str, Any],
    summary: str,
    extras: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    panel_entry = get_compound_panel_entry(marker_name) or {}
    payload = {
        "marker_id": str(panel_entry.get("canonical_name", _normalize_name(marker_name))),
        "display_name": str(panel_entry.get("display_name", marker_name)),
        "target_class": str(panel_entry.get("target_class", "pretreatment_state_markers")),
        "panel_role": str(panel_entry.get("panel_role", "report_only")),
        "observable_kind": str(panel_entry.get("observable_kind", "state_variable")),
        "modeling_regimes": list(panel_entry.get("modeling_regimes", [])),
        "chemistry_family": panel_entry.get("chemistry_family"),
        "supporting_families": list(panel_entry.get("supporting_families", [])),
        "observable_panel_tags": list(panel_entry.get("observable_panel_tags", [])),
        "active": bool(active),
        "state_value": state_value,
        "state_value_summary": state_value_summary,
        "influence_mode": influence_mode,
        "summary": summary,
        "family_lane": {
            "slr_family": str(family_lane.get("slr_family", "")),
            "display_name": str(family_lane.get("display_name", "")),
            "active": bool(family_lane.get("active", False)),
        },
    }
    if extras:
        payload.update(extras)
    return payload


def _build_family_state_markers(
    *,
    thiamine_metadata: Mapping[str, Any],
    thiamine_fraction_estimate: float,
    thiamine_mode: str,
    thiamine_support_lane: Mapping[str, Any],
    nucleotide_support_lane: Mapping[str, Any],
    fermentation_pretreatment_lane: Mapping[str, Any],
    caramelization_lane: Mapping[str, Any],
) -> List[Dict[str, Any]]:
    markers = [
        _build_state_marker_payload(
            "thiamine_availability",
            active=bool(thiamine_support_lane.get("active", False)),
            state_value=bool(thiamine_metadata.get("available", False)),
            state_value_summary=(
                f"available={bool(thiamine_metadata.get('available', False))}, source={thiamine_metadata.get('source', 'unknown')}, mode={thiamine_mode}, mft_fraction_estimate={thiamine_fraction_estimate:.2f}"
            ),
            influence_mode="upstream_state_only",
            family_lane=thiamine_support_lane,
            summary=str(thiamine_support_lane.get("summary", "")),
            extras={
                "source": str(thiamine_metadata.get("source", "unknown")),
                "explicit": bool(thiamine_metadata.get("explicit", False)),
                "inferred_from_inputs": bool(thiamine_metadata.get("inferred_from_inputs", False)),
            },
        ),
        _build_state_marker_payload(
            "nucleotide_enrichment",
            active=bool(nucleotide_support_lane.get("active", False)),
            state_value=float(nucleotide_support_lane.get("nucleotide_support_score", 0.0)),
            state_value_summary=(
                f"support_score={float(nucleotide_support_lane.get('nucleotide_support_score', 0.0)):.2f}, nucleotide_active={bool(nucleotide_support_lane.get('nucleotide_support_active', False))}, ribose_delivery_active={bool(nucleotide_support_lane.get('ribose_delivery_active', False))}"
            ),
            influence_mode="upstream_state_only",
            family_lane=nucleotide_support_lane,
            summary=str(nucleotide_support_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "free_amino_acid_enrichment",
            active=bool(fermentation_pretreatment_lane.get("precursor_release_active", False)),
            state_value=bool(fermentation_pretreatment_lane.get("precursor_release_active", False)),
            state_value_summary=(
                f"precursor_release_active={bool(fermentation_pretreatment_lane.get('precursor_release_active', False))}, pretreatment_support_score={float(fermentation_pretreatment_lane.get('pretreatment_support_score', 0.0)):.2f}"
            ),
            influence_mode="upstream_state_only",
            family_lane=fermentation_pretreatment_lane,
            summary=str(fermentation_pretreatment_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "pretreatment_ph_shift",
            active=bool(fermentation_pretreatment_lane.get("pretreatment_pH_shift_active", False)),
            state_value=bool(fermentation_pretreatment_lane.get("pretreatment_pH_shift_active", False)),
            state_value_summary=(
                f"pretreatment_pH_shift_active={bool(fermentation_pretreatment_lane.get('pretreatment_pH_shift_active', False))}, fermentation_cues_active={bool(fermentation_pretreatment_lane.get('fermentation_cues_active', False))}"
            ),
            influence_mode="upstream_state_only",
            family_lane=fermentation_pretreatment_lane,
            summary=str(fermentation_pretreatment_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "caramelization_severity",
            active=bool(caramelization_lane.get("active", False)),
            state_value=float(caramelization_lane.get("severity_signal_ppb", 0.0)),
            state_value_summary=(
                f"severity_signal_ppb={float(caramelization_lane.get('severity_signal_ppb', 0.0)):.2f}, severity_penalty_factor={float(caramelization_lane.get('severity_penalty_factor', 0.0)):.2f}"
            ),
            influence_mode="upstream_state_plus_marker_panel",
            family_lane=caramelization_lane,
            summary=str(caramelization_lane.get("summary", "")),
            extras={
                "supportive_furanone_context": float(caramelization_lane.get("caramelization_support_score", 0.0)),
            },
        ),
    ]
    return markers


def _build_core_maillard_lane(
    *,
    normalized_sugars: List[str],
    normalized_amino: List[str],
    sulfur_signal: float,
    strecker_balance_score: float,
) -> Dict[str, Any]:
    active = bool(normalized_sugars and normalized_amino)
    core_support_score = min(1.0, 0.45 * float(bool(normalized_sugars)) + 0.35 * float(bool(normalized_amino)) + 0.20 * min(max(strecker_balance_score, 0.0), 1.0))
    summary = "No explicit amino acid-sugar core lane active."
    if active:
        summary = "The amino acid-sugar core lane is active and remains the quantitative trunk for the current formulation."
        if sulfur_signal > 0.0:
            summary = "The amino acid-sugar core lane is active and already supports sulfur-positive routing signals."
    return _build_family_lane_payload(
        "01",
        active=active,
        summary=summary,
        metrics={
            "core_support_score": float(core_support_score),
            "sugar_count": int(len(normalized_sugars)),
            "amino_count": int(len(normalized_amino)),
        },
    )

def _build_positive_adverse_crosstalk_lane(
    *,
    sulfur_signal: float,
    strecker_balance_score: float,
    lipid_crosstalk_lane: Mapping[str, Any],
) -> Dict[str, Any]:
    """
    [Reframe: Family 11 — Positive vs. adverse Maillard/lipid crosstalk closure]

    Hexanal suppression is NOT a standalone heuristic tweak.  It only makes sense
    as part of the broader competitive-flux landscape described by Family 11:
    the same carbonyl pool that drives positive meaty/sulfury Maillard character
    simultaneously competes with lipid-oxidation pathways for available precursors.

    This lane quantifies that competition as a kinetic proxy:
      positive_signal  = sulfur_signal + 10 * strecker_balance_score
      suppression_factor = clipped log10 ratio of positive vs. lipid signal

    Any future compound-specific Kd estimation for hexanal binding should be
    routed here (Family 11) as an upgrade to this heuristic, NOT added to the
    retention pipeline as a separate silo.  Retention/Kd is part of the
    computational closure ladder, not a standalone sprint.

    Source: data/Gemini_Deep_Research/11_maillard_lipid_crosstalk.md
    """
    lipid_marker_signal_ppb = float(lipid_crosstalk_lane.get("lipid_marker_signal_ppb", 0.0))
    positive_signal = sulfur_signal + 10.0 * strecker_balance_score
    active = bool(positive_signal > 1.0 and lipid_marker_signal_ppb > 1.0)
    hexanal_suppression_factor = 0.0
    if active:
        hexanal_suppression_factor = min(1.0, 0.25 * math.log10(1.0 + positive_signal / max(lipid_marker_signal_ppb, 1.0)))
    
    summary = "No positive-vs-adverse crosstalk active."
    if active:
        summary = "Positive Maillard pathways outcompete oxidative pathways, suppressing hexanal/lipid marker accumulation."
        
    return _build_family_lane_payload(
        "11",
        active=active,
        summary=summary,
        metrics={
            "positive_maillard_signal": float(positive_signal),
            "hexanal_suppression_factor": float(hexanal_suppression_factor),
        },
    )


def _build_quinone_cys_sink_lane(
    *,
    normalized_amino: List[str],
    polyphenol_active: bool,
) -> Dict[str, Any]:
    """
    [Reframe: Family 13 — Polyphenol/amino-acid capping as an upstream precursor sink]

    Polyphenol interactions are NOT primarily a safety or antioxidant topic in
    this framework.  They represent a quantitative upstream sink that consumes
    free amino acids (especially Cys > Lys) via quinone conjugation BEFORE those
    amino acids can enter core Maillard pathways.  The priority ordering
    (Cys >> Lys) is dictated by relative nucleophilicity, not by nutritional concern.

    This lane computes a quinone budget that explicitly reduces available
    precursor pool for downstream flavor formation.

    Source: data/Gemini_Deep_Research/13_polyphenol_amino_capping.md
    """
    cys_count = sum(1 for a in normalized_amino if "cysteine" in a or "cys" in a)
    lys_count = sum(1 for a in normalized_amino if "lysine" in a or "lys" in a)
    active = bool(polyphenol_active and (cys_count > 0 or lys_count > 0))
    quinone_budget_consumed = 0.0
    cys_penalty = 0.0
    
    if active:
        if cys_count > 0:
            cys_penalty = 0.8  # Quinone consumes Cys heavily over Lys
            quinone_budget_consumed = 1.0
        elif lys_count > 0:
            quinone_budget_consumed = 0.4
            
    summary = "No quinone-driven precursor sink active."
    if active:
        summary = "Quinones actively deplete available Cysteine/Lysine, throttling downstream meaty/savory Maillard pathways."
        
    return _build_family_lane_payload(
        "13",
        active=active,
        summary=summary,
        metrics={
            "quinone_budget_consumed": float(quinone_budget_consumed),
            "cysteine_depletion_penalty": float(cys_penalty),
        },
    )


def _build_amino_dicarbonyl_source_lane(
    *,
    normalized_amino: List[str],
    core_support_score: float,
) -> Dict[str, Any]:
    """
    [Reframe: Family 14 — Ascorbic acid / AA degradation as a dicarbonyl source term]

    Ascorbic acid and related amino acid degradation pathways contribute to the
    dicarbonyl pool independently of the main sugar fragmentation route.  This
    lane adds that source term to the flavor formation budget rather than treating
    AA as a passive participant.  This is NOT a safety topic; it is a precursor
    routing issue within the computational closure ladder.

    Source: data/Gemini_Deep_Research/14_ascorbic_acid_maillard.md
    """
    active = len(normalized_amino) > 0 and core_support_score > 0.0
    aa_dicarbonyl_flux = 0.0
    if active:
        aa_dicarbonyl_flux = min(1.0, 0.15 * len(normalized_amino) * core_support_score)
        
    summary = "No measurable AA-driven dicarbonyl sourcing."
    if active:
        summary = "Amino acid degradation acts as an explicit source term for reactive dicarbonyls, bypassing standard sugar fragmentation bottlenecks."
        
    return _build_family_lane_payload(
        "14",
        active=active,
        summary=summary,
        metrics={
            "aa_dicarbonyl_flux": float(aa_dicarbonyl_flux),
        },
    )


def _build_sugar_diversion_lane(
    *,
    normalized_sugars: List[str],
    lipids: List[str],
) -> Dict[str, Any]:
    """
    [Reframe: Family 15 — Phospholipid/amine interfacial sugar diversion]

    Phosphatidylethanolamine (PE) acts as one of the most reactive primary amines
    in the early Maillard reaction, especially in water/oil interfaces during
    extrusion.  This lane models the competitive diversion of reducing sugars
    toward PE rather than free amino acids.  It is NOT a lipid oxidation story;
    it is a precursor availability and interfacial reaction acceleration story
    that belongs to the computational closure ladder for Family 15.

    Source: data/Gemini_Deep_Research/15_phospholipid_amine_maillard.md
    """
    pe_active = any("phosphatidylethanolamine" in str(l).lower() or "pe" in str(l).lower().split() for l in lipids)
    active = bool(pe_active and normalized_sugars)
    diversion_factor = 0.0
    
    if active:
        diversion_factor = min(0.8, 0.4 * len(normalized_sugars))
        
    summary = "No PE-driven sugar diversion active."
    if active:
        summary = "Phosphatidylethanolamine acts as a highly reactive primary amine, aggressively diverting early reducing sugars away from standard flavor pathways via interfacial acceleration."
        
    return _build_family_lane_payload(
        "15",
        active=active,
        summary=summary,
        metrics={
            "pe_sugar_diversion_factor": float(diversion_factor),
        },
    )


def _build_trapping_burden_modifier_lane(
    *,
    normalized_additives: List[str],
    protein_label: str,
) -> Dict[str, Any]:
    """
    [Reframe: Family 16 — Melanoidin polymerisation as a macromolecular trapping modifier]

    This lane expresses the generic volatile-trapping burden imposed by high-MW
    melanoidin networks, cross-linked proteins, and polysaccharide matrices.  It
    is intentionally framed as a MODIFIER (not an absolute release curve) because
    compound-specific retention is part of the computational closure ladder and
    must be earned by Kd evidence, not assumed from macromolecular heuristics.

    Safety implications of AGE cross-linking belong to Family 12
    (protein_damage_markers), not here.  This lane only governs flavor release
    attenuation from the polymerised bulk phase.

    Source: data/Gemini_Deep_Research/16_melanoidin_polymerization.md
    """
    macromolecular_burden = 0.0
    if protein_label != "free":
        macromolecular_burden += 0.40
    if any("starch" in a or "hydrocolloid" in a or "fiber" in a for a in normalized_additives):
        macromolecular_burden += 0.35
        
    active = float(macromolecular_burden) > 0.0
    summary = "Matrix entrapment modifying burden is negligible."
    if active:
        summary = "Matrix macromolecules impose a generic trapping burden on volatile release, strictly attenuating absolute yield without pretending to be a compound-specific release curve."
        
    return _build_family_lane_payload(
        "16",
        active=active,
        summary=summary,
        metrics={
            "trapping_burden_modifier": min(1.0, macromolecular_burden),
        },
    )


def _build_upstream_interception_layer(
    family_lane_summary: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """
    [Priority 2B Architecture]
    Groups families 13, 15, and 16 into one explicit upstream runtime layer
    that modifies effective precursor availability before the main flavor solver.
    
    Exposes unified scientist-readable state variables:
      - effective_cysteine_fraction
      - effective_sugar_fraction
      - volatile_retention_factor
    """
    f13 = family_lane_summary.get("13", {})
    f15 = family_lane_summary.get("15", {})
    f16 = family_lane_summary.get("16", {})
    
    f13_metrics = f13.get("metrics", {}) or {}
    f15_metrics = f15.get("metrics", {}) or {}
    f16_metrics = f16.get("metrics", {}) or {}
    
    active = any([f13.get("active", False), f15.get("active", False), f16.get("active", False)])
    
    return {
        "active": active,
        "interceptors": {
            "cysteine_depletion": float(f13_metrics.get("cysteine_depletion_penalty", 0.0)),
            "quinone_load": float(f13_metrics.get("quinone_budget_consumed", 0.0)),
            "sugar_diversion": float(f15_metrics.get("pe_sugar_diversion_factor", 0.0)),
            "trapping_burden": float(f16_metrics.get("trapping_burden_modifier", 0.0)),
        },
        "state_variables": {
            "effective_cysteine_fraction": round(1.0 - float(f13_metrics.get("cysteine_depletion_penalty", 0.0)), 4),
            "effective_sugar_fraction": round(1.0 - float(f15_metrics.get("pe_sugar_diversion_factor", 0.0)), 4),
            "volatile_retention_factor": round(1.0 - float(f16_metrics.get("trapping_burden_modifier", 0.0)), 4),
        },
        "summary": "Unified layer for Families 13 (Polyphenol Sink), 15 (PE Diversion), and 16 (Trapping Burden)."
    }



def _classify_donor_family(normalized_sugar: str) -> str:
    if any(token in normalized_sugar for token in ["phosphate", "ribose 5 phosphate", "glucose 6 phosphate", "fructose 6 phosphate"]):
        return "phosphorylated"
    if any(token in normalized_sugar for token in ["ribose", "xylose", "arabinose"]):
        return "pentose"
    if "fructose" in normalized_sugar:
        return "fructose"
    if "glucose" in normalized_sugar:
        return "glucose"
    return "other_reducing_sugar"


def _build_donor_hierarchy_lane(normalized_sugars: List[str]) -> Dict[str, Any]:
    donor_counts: Dict[str, int] = {}
    donor_inputs: Dict[str, List[str]] = {}
    donor_strength = {
        "phosphorylated": 1.0,
        "pentose": 0.95,
        "fructose": 0.75,
        "glucose": 0.50,
        "other_reducing_sugar": 0.35,
    }
    for sugar in normalized_sugars:
        donor_class = _classify_donor_family(sugar)
        donor_counts[donor_class] = donor_counts.get(donor_class, 0) + 1
        donor_inputs.setdefault(donor_class, []).append(sugar)

    dominant_donor_class = "none"
    if donor_counts:
        dominant_donor_class = max(
            donor_counts,
            key=lambda donor_class: (donor_strength.get(donor_class, 0.0), donor_counts[donor_class], donor_class),
        )

    active = bool(normalized_sugars)
    donor_hierarchy_score = donor_strength.get(dominant_donor_class, 0.0)
    mixed_donor_system = len(donor_counts) > 1
    summary = "No explicit carbonyl donor hierarchy lane active."
    if active:
        summary = (
            f"Dominant donor class is {dominant_donor_class}, so sugar identity should be treated as a first-class routing variable rather than generic sugar loading."
        )
    return _build_family_lane_payload(
        "07",
        active=active,
        summary=summary,
        metrics={
            "dominant_donor_class": dominant_donor_class,
            "donor_class_counts": donor_counts,
            "donor_inputs": donor_inputs,
            "donor_hierarchy_score": float(donor_hierarchy_score),
            "mixed_donor_system": bool(mixed_donor_system),
            "supports_fast_sulfur_routing": dominant_donor_class in {"pentose", "phosphorylated"},
        },
    )


def _build_thiamine_support_lane(
    *,
    thiamine_active: bool,
    thiamine_fraction_estimate: float,
    thiamine_mode: str,
) -> Dict[str, Any]:
    summary = "No explicit thiamine support lane active."
    if thiamine_active:
        summary = "Thiamine support is active, so sulfur routing should be interpreted as augmented rather than purely core-derived."
    return _build_family_lane_payload(
        "03",
        active=thiamine_active,
        summary=summary,
        metrics={
            "thiamine_support_score": float(thiamine_fraction_estimate),
            "thiamine_mode": thiamine_mode,
        },
    )


def _build_nucleotide_support_lane(
    *,
    normalized_sugars: List[str],
    normalized_additives: List[str],
    normalized_interventions: List[str],
) -> Dict[str, Any]:
    support_pool = normalized_sugars + normalized_additives + normalized_interventions
    nucleotide_tokens = ["imp", "gmp", "inosinate", "guanylate", "yeast extract"]
    nucleotide_active = any(any(token in value for token in nucleotide_tokens) for value in support_pool)
    ribose_delivery_active = any(any(token in value for token in ["ribose", "ribose 5 phosphate"]) for value in normalized_sugars + normalized_additives)
    nucleotide_support_score = min(1.0, 0.55 * float(nucleotide_active) + 0.45 * float(ribose_delivery_active))
    summary = "No explicit nucleotide support lane active."
    if nucleotide_active or ribose_delivery_active:
        summary = "Nucleotide or ribose-delivery support is active, so umami amplification should be treated as an upstream support lane."
    return _build_family_lane_payload(
        "04",
        active=bool(nucleotide_active or ribose_delivery_active),
        summary=summary,
        metrics={
            "nucleotide_support_active": bool(nucleotide_active),
            "ribose_delivery_active": bool(ribose_delivery_active),
            "nucleotide_support_score": float(nucleotide_support_score),
        },
    )


def _build_sulfur_peptide_support_lane(
    *,
    normalized_additives: List[str],
    normalized_amino: List[str],
    normalized_interventions: List[str],
) -> Dict[str, Any]:
    support_pool = normalized_additives + normalized_amino + normalized_interventions
    glutathione_active = any("glutathione" in value for value in support_pool)
    peptide_support_active = any(any(token in value for token in ["peptide", "hydrolysate", "autolysate"]) for value in support_pool)
    sulfur_peptide_support_score = min(1.0, 0.60 * float(glutathione_active) + 0.40 * float(peptide_support_active))
    summary = "No explicit glutathione or peptide support lane active."
    if glutathione_active or peptide_support_active:
        summary = "Glutathione or peptide support is active, so sulfur intensity should be interpreted as partly matrix-supported rather than purely free-cysteine-driven."
    return _build_family_lane_payload(
        "05",
        active=bool(glutathione_active or peptide_support_active),
        summary=summary,
        metrics={
            "glutathione_active": bool(glutathione_active),
            "peptide_support_active": bool(peptide_support_active),
            "sulfur_peptide_support_score": float(sulfur_peptide_support_score),
        },
    )


def _build_matrix_scope_lane(*, protein_label: str) -> Dict[str, Any]:
    normalized_protein = _normalize_name(protein_label)
    alternative_matrix_active = bool(normalized_protein and normalized_protein not in {"free", "pea_iso", "pea_conc", "soy_iso", "soy_conc"})
    matrix_uncertainty_factor = 0.75 if alternative_matrix_active else 0.0
    summary = "No alternative-matrix scope lane active."
    if alternative_matrix_active:
        summary = "Alternative protein matrix scope is active, so matrix-transfer uncertainty should be treated as a first-class recommendation constraint."
    return _build_family_lane_payload(
        "06",
        active=alternative_matrix_active,
        summary=summary,
        metrics={
            "matrix_scope_active": bool(alternative_matrix_active),
            "matrix_uncertainty_factor": float(matrix_uncertainty_factor),
            "protein_type": protein_label,
        },
    )


def _build_fermentation_pretreatment_lane(
    *,
    normalized_additives: List[str],
    normalized_amino: List[str],
    normalized_interventions: List[str],
    pH: Optional[float],
    thiamine_active: bool,
) -> Dict[str, Any]:
    pretreatment_pool = normalized_additives + normalized_amino + normalized_interventions
    fermentation_active = any(
        any(token in value for token in ["ferment", "koji", "miso", "starter culture", "culture", "yeast extract"])
        for value in pretreatment_pool
    )
    precursor_release_active = any(
        any(token in value for token in ["hydrolysate", "hydroly", "protease", "peptide", "autolysate"])
        for value in pretreatment_pool
    )
    off_note_cleanup_active = any(
        any(token in value for token in ["yeast fermentation", "yeast extract", "koji", "ferment"])
        for value in pretreatment_pool
    )
    nucleotide_support_active = any(
        any(token in value for token in ["imp", "gmp", "inosinate", "guanylate", "yeast extract"])
        for value in pretreatment_pool
    )
    pretreatment_pH_shift_active = bool((fermentation_active or off_note_cleanup_active) and pH is not None and float(pH) < 6.2)
    active = bool(
        fermentation_active
        or precursor_release_active
        or off_note_cleanup_active
        or nucleotide_support_active
    )
    pretreatment_support_score = min(
        1.0,
        0.35 * float(precursor_release_active)
        + 0.30 * float(off_note_cleanup_active)
        + 0.20 * float(nucleotide_support_active)
        + 0.15 * float(thiamine_active),
    )
    summary = "No explicit fermentation pretreatment lane active."
    if active:
        summary = (
            "Pretreatment cues indicate upstream fermentation or hydrolysis support, so precursor loading and off-note burden should be interpreted before thermal chemistry."
        )
    return _build_family_lane_payload(
        "10",
        active=active,
        summary=summary,
        metrics={
            "fermentation_cues_active": bool(fermentation_active),
            "precursor_release_active": bool(precursor_release_active),
            "off_note_cleanup_active": bool(off_note_cleanup_active),
            "nucleotide_support_active": bool(nucleotide_support_active),
            "pretreatment_pH_shift_active": bool(pretreatment_pH_shift_active),
            "pretreatment_support_score": float(pretreatment_support_score),
        },
    )


def _build_caramelization_lane(
    *,
    signal_map: Dict[str, float],
    furanone_expected: List[str],
    furanone_observed: List[str],
) -> Dict[str, Any]:
    severity_marker_signals = {
        "furfural": _compound_signal(signal_map, ["furfural"]),
        "hmf": _compound_signal(signal_map, ["5-hydroxymethylfurfural", "5-hydroxymethylfurfural (hmf)", "hmf"]),
        "2_acetylfuran": _compound_signal(signal_map, ["2-acetylfuran", "2 acetylfuran"]),
    }
    severity_signal_ppb = sum(float(value) for value in severity_marker_signals.values())
    supportive_furanones_observed = len(furanone_observed)
    active = bool(severity_signal_ppb > 0.0 or furanone_expected)
    caramelization_support_score = min(1.0, 0.40 * float(bool(furanone_expected)) + 0.30 * min(supportive_furanones_observed, 2) / 2.0)
    severity_penalty_factor = min(1.5, severity_signal_ppb / 60.0)
    summary = "No explicit carbohydrate pyrolysis or caramelization lane active."
    if active:
        summary = "Carbohydrate pyrolysis or caramelization markers are active, so helpful browning support must be separated from over-furan drift."
    return _build_family_lane_payload(
        "09",
        active=active,
        summary=summary,
        metrics={
            "severity_marker_signals_ppb": severity_marker_signals,
            "severity_signal_ppb": float(severity_signal_ppb),
            "caramelization_support_score": float(caramelization_support_score),
            "severity_penalty_factor": float(severity_penalty_factor),
        },
    )


def _build_off_note_guardrail_lane(
    *,
    normalized_additives: List[str],
    normalized_interventions: List[str],
    lipid_crosstalk_lane: Mapping[str, Any],
    polyphenol_active: bool,
    protein_label: str,
) -> Dict[str, Any]:
    guardrail_pool = normalized_additives + normalized_interventions
    antioxidant_guardrail_active = any(
        any(token in value for token in ["rosemary", "catechin", "tannic acid", "green tea extract", "grape seed extract", "polyphenol"])
        for value in guardrail_pool
    )
    calcium_guardrail_active = any("calcium carbonate" in value for value in guardrail_pool)
    lipid_marker_signal_ppb = float(lipid_crosstalk_lane.get("lipid_marker_signal_ppb", 0.0))
    dicarbonyl_trapping_factor = 1.0 if polyphenol_active else (0.4 if antioxidant_guardrail_active else 0.0)
    amino_group_blocking_factor = min(1.5, lipid_marker_signal_ppb / 120.0) if protein_label != "free" else 0.0
    suppression_pressure_active = bool(
        lipid_marker_signal_ppb > 0.0
        or antioxidant_guardrail_active
        or calcium_guardrail_active
        or polyphenol_active
    )
    summary = "No explicit off-note guardrail lane active."
    if suppression_pressure_active:
        summary = (
            "Plant-matrix guardrails are active, so off-note suppression and amino-group blocking risk should constrain optimistic flavor claims."
        )
    return _build_family_lane_payload(
        "08",
        active=suppression_pressure_active,
        summary=summary,
        metrics={
            "suppression_pressure_active": bool(suppression_pressure_active),
            "antioxidant_guardrail_active": bool(antioxidant_guardrail_active),
            "polyphenol_guardrail_active": bool(polyphenol_active),
            "acrylamide_guardrail_active": bool(calcium_guardrail_active),
            "dicarbonyl_trapping_factor": float(dicarbonyl_trapping_factor),
            "amino_group_blocking_factor": float(amino_group_blocking_factor),
        },
    )


def _build_family_lane_adjustments(family_lane_summary: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    per_lane: Dict[str, Dict[str, float]] = {}

    core_lane = family_lane_summary.get("01", {})
    per_lane["01"] = {
        "target_score_delta": 0.08 * float(core_lane.get("core_support_score", 0.0)) if core_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.0,
    }

    lipid_lane = family_lane_summary.get("02", {})
    per_lane["02"] = {
        "target_score_delta": -0.06 * min(1.5, float(lipid_lane.get("lipid_marker_signal_ppb", 0.0)) / 100.0) if lipid_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.14 * float(lipid_lane.get("maillard_closure_pressure", 0.0)) if lipid_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.18 * min(1.5, float(lipid_lane.get("lipid_marker_signal_ppb", 0.0)) / 100.0) if lipid_lane.get("active") else 0.0,
    }

    thiamine_lane = family_lane_summary.get("03", {})
    per_lane["03"] = {
        "target_score_delta": 0.18 * float(thiamine_lane.get("thiamine_support_score", 0.0)) if thiamine_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.0,
    }

    nucleotide_lane = family_lane_summary.get("04", {})
    per_lane["04"] = {
        "target_score_delta": 0.16 * float(nucleotide_lane.get("nucleotide_support_score", 0.0)) if nucleotide_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.03 * float(nucleotide_lane.get("nucleotide_support_active", False)) if nucleotide_lane.get("active") else 0.0,
    }

    peptide_lane = family_lane_summary.get("05", {})
    per_lane["05"] = {
        "target_score_delta": 0.12 * float(peptide_lane.get("sulfur_peptide_support_score", 0.0)) if peptide_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.0,
    }

    matrix_lane = family_lane_summary.get("06", {})
    per_lane["06"] = {
        "target_score_delta": -0.12 * float(matrix_lane.get("matrix_uncertainty_factor", 0.0)) if matrix_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.05 * float(matrix_lane.get("matrix_uncertainty_factor", 0.0)) if matrix_lane.get("active") else 0.0,
    }

    donor_lane = family_lane_summary.get("07", {})
    donor_score = float(donor_lane.get("donor_hierarchy_score", 0.0))
    supports_fast = bool(donor_lane.get("supports_fast_sulfur_routing", False))
    per_lane["07"] = {
        "target_score_delta": (0.14 * donor_score) if donor_lane.get("active") else 0.0,
        "off_flavour_risk_delta": (-0.04 * donor_score) if donor_lane.get("active") and supports_fast else 0.02 * max(0.0, 0.6 - donor_score),
    }

    guardrail_lane = family_lane_summary.get("08", {})
    per_lane["08"] = {
        "target_score_delta": -0.08 * float(guardrail_lane.get("amino_group_blocking_factor", 0.0)) if guardrail_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.18 * float(guardrail_lane.get("dicarbonyl_trapping_factor", 0.0)) + 0.08 * float(guardrail_lane.get("amino_group_blocking_factor", 0.0)) if guardrail_lane.get("active") else 0.0,
    }

    caramelization_lane = family_lane_summary.get("09", {})
    per_lane["09"] = {
        "target_score_delta": 0.08 * float(caramelization_lane.get("caramelization_support_score", 0.0)) - 0.10 * float(caramelization_lane.get("severity_penalty_factor", 0.0)) if caramelization_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.06 * float(caramelization_lane.get("severity_penalty_factor", 0.0)) if caramelization_lane.get("active") else 0.0,
    }

    fermentation_lane = family_lane_summary.get("10", {})
    off_note_cleanup_bonus = 0.08 if bool(fermentation_lane.get("off_note_cleanup_active", False)) else 0.0
    per_lane["10"] = {
        "target_score_delta": 0.16 * float(fermentation_lane.get("pretreatment_support_score", 0.0)) if fermentation_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -off_note_cleanup_bonus if fermentation_lane.get("active") else 0.0,
    }

    family_11_lane = family_lane_summary.get("11", {})
    per_lane["11"] = {
        "target_score_delta": 0.05 * float(family_11_lane.get("hexanal_suppression_factor", 0.0)) if family_11_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.15 * float(family_11_lane.get("hexanal_suppression_factor", 0.0)) if family_11_lane.get("active") else 0.0,
    }

    family_13_lane = family_lane_summary.get("13", {})
    per_lane["13"] = {
        "target_score_delta": -0.12 * float(family_13_lane.get("cysteine_depletion_penalty", 0.0)) if family_13_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.05 * float(family_13_lane.get("quinone_budget_consumed", 0.0)) if family_13_lane.get("active") else 0.0,
    }

    family_14_lane = family_lane_summary.get("14", {})
    per_lane["14"] = {
        "target_score_delta": 0.15 * float(family_14_lane.get("aa_dicarbonyl_flux", 0.0)) if family_14_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.05 * float(family_14_lane.get("aa_dicarbonyl_flux", 0.0)) if family_14_lane.get("active") else 0.0,
    }

    family_15_lane = family_lane_summary.get("15", {})
    per_lane["15"] = {
        "target_score_delta": -0.18 * float(family_15_lane.get("pe_sugar_diversion_factor", 0.0)) if family_15_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.08 * float(family_15_lane.get("pe_sugar_diversion_factor", 0.0)) if family_15_lane.get("active") else 0.0,
    }

    family_16_lane = family_lane_summary.get("16", {})
    per_lane["16"] = {
        "target_score_delta": -0.20 * float(family_16_lane.get("trapping_burden_modifier", 0.0)) if family_16_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.10 * float(family_16_lane.get("trapping_burden_modifier", 0.0)) if family_16_lane.get("active") else 0.0,
    }

    total_target_score_delta = sum(float(row.get("target_score_delta", 0.0)) for row in per_lane.values())
    total_maillard_closure_delta = sum(float(row.get("maillard_closure_delta", 0.0)) for row in per_lane.values())
    total_off_flavour_risk_delta = sum(float(row.get("off_flavour_risk_delta", 0.0)) for row in per_lane.values())
    return {
        "per_lane": per_lane,
        "target_score_delta": float(total_target_score_delta),
        "maillard_closure_delta": float(total_maillard_closure_delta),
        "off_flavour_risk_delta": float(total_off_flavour_risk_delta),
    }


def _reference_panel_value(payload_section: str, entry_id: str, panel_key: str) -> float:
    for entry in query_flavor_reference_entries(entry_id=entry_id, payload_section=payload_section):
        values = entry.get("numeric_band_or_point", {}).get("values", {})
        if panel_key in values:
            return float(values[panel_key])
    return 0.0


def _resolve_thiamine_availability(
    thiamine_availability: Any,
    *,
    normalized_sugars: List[str],
    normalized_amino: List[str],
    normalized_additives: List[str],
    protein_label: str,
) -> Dict[str, Any]:
    inferred_from_inputs = any(
        "thiamine" in value or "vitamin b1" in value
        for value in normalized_additives + normalized_amino + normalized_sugars
    )

    if isinstance(thiamine_availability, Mapping):
        available = bool(thiamine_availability.get("available", False))
        source = str(
            thiamine_availability.get(
                "source",
                "explicit_formulation_metadata" if available else "explicitly_disabled",
            )
        )
        explicit = True
    elif isinstance(thiamine_availability, bool):
        available = thiamine_availability
        source = "explicit_formulation_metadata" if available else "explicitly_disabled"
        explicit = True
    elif isinstance(thiamine_availability, str) and thiamine_availability.strip():
        normalized = thiamine_availability.strip().lower().replace("-", "_").replace(" ", "_")
        positive_tokens = {"present", "available", "fortified", "explicit_additive", "pbma_fortified", "added_thiamine"}
        negative_tokens = {"absent", "unavailable", "disabled", "native_matrix_default", "benchmark_native_default"}
        available = normalized in positive_tokens
        if normalized in negative_tokens:
            available = False
        source = normalized
        explicit = True
    else:
        available = inferred_from_inputs
        explicit = False
        if inferred_from_inputs:
            source = "ingredient_list_inference"
        elif protein_label.startswith(("pea", "soy")):
            source = "native_matrix_default_inactive"
        else:
            source = "no_thiamine_evidence"

    return {
        "available": bool(available),
        "source": source,
        "explicit": bool(explicit),
        "inferred_from_inputs": bool(inferred_from_inputs),
    }


def build_family_upstream_contract(
    *,
    sugars: Optional[List[str]] = None,
    amino_acids: Optional[List[str]] = None,
    additives: Optional[List[str]] = None,
    interventions: Optional[List[Any]] = None,
    protein_type: Optional[str] = None,
    pH: Optional[float] = None,
    thiamine_availability: Any = None,
    molar_ratios: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    sugars = sugars or []
    amino_acids = amino_acids or []
    additives = additives or []
    normalized_sugars = [_normalize_name(value) for value in sugars]
    normalized_additives = [_normalize_name(value) for value in additives]
    normalized_amino = [_normalize_name(value) for value in amino_acids]
    normalized_interventions = _normalize_intervention_tokens(interventions)
    protein_label = str(protein_type or "free")

    thiamine_metadata = _resolve_thiamine_availability(
        thiamine_availability,
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        protein_label=protein_label,
    )
    thiamine_active = bool(thiamine_metadata["available"])

    thiamine_priors = get_thiamine_priors()
    pentose_active = any(token in sugar for sugar in normalized_sugars for token in ["ribose", "xylose", "arabinose"])
    thiamine_fraction_estimate = 0.0
    thiamine_mode = "inactive"
    if thiamine_active:
        if pentose_active:
            thiamine_fraction_estimate = float(thiamine_priors.get("mixed_pentose_fraction", 0.5) or 0.5)
            thiamine_mode = "mixed_thiamine_plus_pentose"
        else:
            thiamine_fraction_estimate = float(thiamine_priors.get("thiamine_only_fraction", 1.0) or 1.0)
            thiamine_mode = "thiamine_only"

    donor_hierarchy_lane = _build_donor_hierarchy_lane(normalized_sugars)
    fermentation_pretreatment_lane = _build_fermentation_pretreatment_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_interventions=normalized_interventions,
        pH=pH,
        thiamine_active=thiamine_active,
    )
    thiamine_support_lane = _build_thiamine_support_lane(
        thiamine_active=thiamine_active,
        thiamine_fraction_estimate=thiamine_fraction_estimate,
        thiamine_mode=thiamine_mode,
    )
    nucleotide_support_lane = _build_nucleotide_support_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
    )
    sulfur_peptide_support_lane = _build_sulfur_peptide_support_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_interventions=normalized_interventions,
    )
    matrix_scope_lane = _build_matrix_scope_lane(protein_label=protein_label)

    donor_weight_by_class = {
        "phosphorylated": 1.35,
        "pentose": 1.20,
        "fructose": 1.05,
        "glucose": 0.92,
        "other_reducing_sugar": 0.88,
    }
    effective_molar_ratios = {
        str(key): float(value)
        for key, value in (molar_ratios or {}).items()
        if value is not None
    }
    donor_pool_factors: Dict[str, float] = {}
    ratio_key_lookup = {_normalize_name(key): key for key in effective_molar_ratios}
    dominant_donor_class = str(donor_hierarchy_lane.get("dominant_donor_class", "none"))
    mixed_donor_system = bool(donor_hierarchy_lane.get("mixed_donor_system", False))
    for sugar in normalized_sugars:
        donor_class = _classify_donor_family(sugar)
        factor = float(donor_weight_by_class.get(donor_class, 1.0))
        if mixed_donor_system and dominant_donor_class in {"phosphorylated", "pentose"} and donor_class != dominant_donor_class:
            factor *= 0.85
        donor_pool_factors[sugar] = float(factor)
        ratio_key = ratio_key_lookup.get(sugar)
        if ratio_key is not None:
            effective_molar_ratios[ratio_key] = float(effective_molar_ratios[ratio_key]) * factor

    effective_pH = None if pH is None else float(pH)
    if fermentation_pretreatment_lane.get("pretreatment_pH_shift_active") and effective_pH is not None:
        effective_pH = max(4.8, effective_pH - 0.35)

    added_precursors: List[str] = []
    added_precursor_ratios: Dict[str, float] = {}
    support_pool = normalized_sugars + normalized_additives + normalized_amino
    if thiamine_active and not any("thiamine" in value or "vitamin b1" in value for value in support_pool):
        added_precursors.append("thiamine")
        added_precursor_ratios["thiamine"] = float(max(0.12, 0.25 * max(thiamine_fraction_estimate, 0.5)))

    return {
        "effective_molar_ratios": effective_molar_ratios,
        "effective_pH": effective_pH,
        "pretreatment_active": bool(fermentation_pretreatment_lane.get("active", False)),
        "pretreatment_interventions": [
            token
            for token in normalized_interventions
            if token in {"yeast fermentation", "protease hydrolysis", "yeast_fermentation", "protease_hydrolysis"}
        ],
        "donor_pool_factors": donor_pool_factors,
        "donor_limited": bool(donor_hierarchy_lane.get("active", False) and dominant_donor_class in {"glucose", "other_reducing_sugar", "none"}),
        "dominant_donor_class": dominant_donor_class,
        "supports_fast_sulfur_routing": bool(donor_hierarchy_lane.get("supports_fast_sulfur_routing", False)),
        "added_precursors": added_precursors,
        "added_precursor_ratios": added_precursor_ratios,
        "thiamine_metadata": thiamine_metadata,
        "thiamine_fraction_estimate": float(thiamine_fraction_estimate),
        "thiamine_mode": thiamine_mode,
        "family_lanes": {
            "03": thiamine_support_lane,
            "04": nucleotide_support_lane,
            "05": sulfur_peptide_support_lane,
            "06": matrix_scope_lane,
            "07": donor_hierarchy_lane,
            "10": fermentation_pretreatment_lane,
        },
        "summary": {
            "donor_identity_first_class": bool(donor_hierarchy_lane.get("active", False)),
            "fermentation_pretreatment_active": bool(fermentation_pretreatment_lane.get("active", False)),
            "nucleotide_support_active": bool(nucleotide_support_lane.get("active", False)),
            "sulfur_peptide_support_active": bool(sulfur_peptide_support_lane.get("active", False)),
            "matrix_scope_active": bool(matrix_scope_lane.get("active", False)),
        },
    }


def build_flavor_axis_summary(
    *,
    projection_metadata: ProjectionMetadataMap,
    sugars: Optional[List[str]] = None,
    amino_acids: Optional[List[str]] = None,
    additives: Optional[List[str]] = None,
    lipids: Optional[List[str]] = None,
    interventions: Optional[List[Any]] = None,
    protein_type: Optional[str] = None,
    pH: Optional[float] = None,
    thiamine_availability: Any = None,
    molar_ratios: Optional[Mapping[str, Any]] = None,
    family_upstream_contract: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    signal_map = _projection_rows_to_signal_map(projection_metadata)
    sugars = sugars or []
    amino_acids = amino_acids or []
    additives = additives or []
    lipids = lipids or []
    protein_label = str(protein_type or "free")

    mft_signal = _compound_signal(signal_map, ["2-Methyl-3-furanthiol (MFT)", "2-methyl-3-furanthiol", "MFT"])
    fft_signal = _compound_signal(signal_map, ["2-Furfurylthiol (FFT)", "2-furfurylthiol", "FFT"])
    sulfur_signal = max(mft_signal, fft_signal)
    methional_signal = _compound_signal(signal_map, ["Methional (3-(methylthio)propanal)", "methional"])
    two_methylbutanal_signal = _compound_signal(signal_map, ["2-Methylbutanal", "2-methylbutanal"])
    three_methylbutanal_signal = _compound_signal(signal_map, ["3-Methylbutanal", "3-methylbutanal"])
    phenylacetaldehyde_signal = _compound_signal(signal_map, ["Phenylacetaldehyde", "phenylacetaldehyde"])
    benzaldehyde_signal = _compound_signal(signal_map, ["Benzaldehyde", "benzaldehyde"])

    two_methylbutanal_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_2_methylbutanal_panel",
        "regular_ground_beef",
    ) or 23.51
    methional_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_methional_panel",
        "regular_ground_beef",
    ) or 4.66
    three_methylbutanal_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_3_methylbutanal_panel",
        "regular_ground_beef",
    ) or 16.91
    phenylacetaldehyde_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_phenylacetaldehyde_panel",
        "regular_ground_beef",
    ) or 11.85
    benzaldehyde_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_benzaldehyde_panel",
        "regular_ground_beef",
    ) or 23.65

    weighted_strecker = (
        (2.0 * min(two_methylbutanal_signal / max(two_methylbutanal_target, 1.0e-6), 1.5) if _role_supports_runtime_scoring(get_flavor_reference_pipeline_role("hernandez_2023_2_methylbutanal_panel")) else 0.0)
        + (1.0 * min(methional_signal / max(methional_target, 1.0e-6), 1.5) if _role_supports_runtime_scoring(get_flavor_reference_pipeline_role("hernandez_2023_methional_panel")) else 0.0)
        + (0.6 * min(phenylacetaldehyde_signal / max(phenylacetaldehyde_target, 1.0e-6), 1.5) if _role_supports_runtime_scoring(get_flavor_reference_pipeline_role("hernandez_2023_phenylacetaldehyde_panel")) else 0.0)
    )
    strecker_balance_score = weighted_strecker / 3.6
    strecker_gap_penalty = 0.0
    if sulfur_signal > 0.25:
        strecker_gap_penalty = max(0.0, (0.4 - strecker_balance_score) * 2.5)

    pyrazine_priors = get_pyrazine_control_priors()
    normalized_sugars = [_normalize_name(value) for value in sugars]
    normalized_additives = [_normalize_name(value) for value in additives]
    normalized_amino = [_normalize_name(value) for value in amino_acids]
    normalized_interventions = _normalize_intervention_tokens(interventions)

    upstream_contract = dict(
        family_upstream_contract
        or build_family_upstream_contract(
            sugars=sugars,
            amino_acids=amino_acids,
            additives=additives,
            interventions=interventions,
            protein_type=protein_type,
            pH=pH,
            thiamine_availability=thiamine_availability,
            molar_ratios=molar_ratios,
        )
    )

    pyrazine_signal = _pyrazine_signal(signal_map)
    ph_factor = _sigmoid(float(pH if pH is not None else 6.0), 6.7, 0.85)

    sugar_bias = 0.0
    if any(token in sugar for sugar in normalized_sugars for token in ["fructose"]):
        sugar_bias += 0.55
    if any(token in sugar for sugar in normalized_sugars for token in ["xylose", "arabinose", "ribose"]):
        sugar_bias += 0.35
    if any(token in sugar for sugar in normalized_sugars for token in ["glucose"]):
        sugar_bias += 0.18

    peptide_bias = 0.0
    if protein_label != "free":
        peptide_bias += 0.10
    if any(token in value for value in normalized_additives + normalized_amino for token in ["peptide", "hydrolysate", "glutathione"]):
        peptide_bias += 0.25

    pyrazine_propensity = min(
        float(pyrazine_priors.get("propensity_cap", 1.75) or 1.75),
        ph_factor + sugar_bias + peptide_bias,
    )
    pyrazine_burden = pyrazine_signal * (1.0 + 0.5 * pyrazine_propensity)
    pyrazine_penalty = 0.0
    if pyrazine_burden > 20.0:
        pyrazine_penalty = math.log10(1.0 + pyrazine_burden / 20.0) * (0.55 + pyrazine_propensity)

    furanone_priors = get_furanone_priors()
    pentose_active = any(token in sugar for sugar in normalized_sugars for token in ["ribose", "xylose", "arabinose"])
    alanine_active = any("alanine" in value for value in normalized_additives + normalized_amino)
    glycine_active = any("glycine" in value for value in normalized_additives + normalized_amino)
    endogenous_proxy = protein_label != "free"
    furanone_expected: List[str] = []
    if pentose_active and (alanine_active or endogenous_proxy):
        furanone_expected.append("HEMF")
    if pentose_active and (glycine_active or alanine_active or endogenous_proxy):
        furanone_expected.append("DMHF")
    furanone_observed = [
        compound
        for compound in ["HEMF", "DMHF"]
        if signal_map.get(_normalize_name(compound), 0.0) > 0.0
    ]
    missing_furanones = [compound for compound in furanone_expected if compound not in furanone_observed]
    furanone_support_score = 1.0
    furanone_penalty = 0.0
    if furanone_expected:
        furanone_support_score = len(furanone_observed) / len(furanone_expected)
        confidence_scale = {
            "low": 0.35,
            "medium": 0.60,
            "high": 0.90,
        }.get(str(furanone_priors.get("confidence_tier", "low")).lower(), 0.35)
        furanone_penalty = (1.0 - furanone_support_score) * confidence_scale

    thiamine_priors = get_thiamine_priors()
    thiamine_metadata = dict(upstream_contract.get("thiamine_metadata", {})) or _resolve_thiamine_availability(
        thiamine_availability,
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        protein_label=protein_label,
    )
    thiamine_active = bool(thiamine_metadata["available"])
    thiamine_fraction_estimate = float(upstream_contract.get("thiamine_fraction_estimate", 0.0) or 0.0)
    thiamine_mode = str(upstream_contract.get("thiamine_mode", "inactive"))
    if thiamine_active and thiamine_fraction_estimate <= 0.0:
        if pentose_active:
            thiamine_fraction_estimate = float(thiamine_priors.get("mixed_pentose_fraction", 0.5) or 0.5)
            thiamine_mode = "mixed_thiamine_plus_pentose"
        else:
            thiamine_fraction_estimate = float(thiamine_priors.get("thiamine_only_fraction", 1.0) or 1.0)
            thiamine_mode = "thiamine_only"

    strecker_crosstalk_priors = get_strecker_crosstalk_priors()
    polyphenol_tokens = [
        str(token).strip().lower()
        for token in strecker_crosstalk_priors.get("polyphenol_examples", [])
        if str(token).strip()
    ] or ["catechin", "tannic acid", "green tea extract", "grape seed extract", "polyphenol"]
    required_sugars = [
        str(token).strip().lower()
        for token in strecker_crosstalk_priors.get("required_sugars", [])
        if str(token).strip()
    ] or ["glucose"]
    polyphenol_active = any(
        token in value
        for value in normalized_additives
        for token in polyphenol_tokens
    )
    crosstalk_active = bool(
        polyphenol_active
        and lipids
        and any(any(required in sugar for required in required_sugars) for sugar in normalized_sugars)
    )
    lincoln_crosstalk_prior = {
        "active": crosstalk_active,
        "source": str(strecker_crosstalk_priors.get("source", "Lincoln et al. (2025), Sustainable Food Proteins")),
        "confidence_tier": str(strecker_crosstalk_priors.get("confidence_tier", "medium_low")),
        "effect_direction": str(strecker_crosstalk_priors.get("effect_direction", "suppress_strecker_and_moderate_oxidative_crosstalk")),
        "summary": (
            str(strecker_crosstalk_priors.get("summary", "Glucose plus lipid plus polyphenol systems can suppress Strecker aldehydes through alpha-dicarbonyl competition."))
            if crosstalk_active
            else "Lincoln 2025 qualitative prior inactive for this formulation context."
        ),
    }

    lipid_crosstalk_lane = _build_lipid_crosstalk_lane(
        signal_map=signal_map,
        normalized_sugars=normalized_sugars,
        lipids=lipids,
        polyphenol_active=polyphenol_active,
        protein_type=protein_label,
    )
    donor_hierarchy_lane = _build_donor_hierarchy_lane(normalized_sugars)
    fermentation_pretreatment_lane = _build_fermentation_pretreatment_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_interventions=normalized_interventions,
        pH=pH,
        thiamine_active=thiamine_active,
    )
    off_note_guardrail_lane = _build_off_note_guardrail_lane(
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        lipid_crosstalk_lane=lipid_crosstalk_lane,
        polyphenol_active=polyphenol_active,
        protein_label=protein_label,
    )
    core_maillard_lane = _build_core_maillard_lane(
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        sulfur_signal=sulfur_signal,
        strecker_balance_score=strecker_balance_score,
    )
    thiamine_support_lane = _build_thiamine_support_lane(
        thiamine_active=thiamine_active,
        thiamine_fraction_estimate=thiamine_fraction_estimate,
        thiamine_mode=thiamine_mode,
    )
    nucleotide_support_lane = _build_nucleotide_support_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
    )
    sulfur_peptide_support_lane = _build_sulfur_peptide_support_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_interventions=normalized_interventions,
    )
    matrix_scope_lane = _build_matrix_scope_lane(protein_label=protein_label)
    caramelization_lane = _build_caramelization_lane(
        signal_map=signal_map,
        furanone_expected=furanone_expected,
        furanone_observed=furanone_observed,
    )
    family_11_lane = _build_positive_adverse_crosstalk_lane(
        sulfur_signal=sulfur_signal,
        strecker_balance_score=strecker_balance_score,
        lipid_crosstalk_lane=lipid_crosstalk_lane,
    )
    family_13_lane = _build_quinone_cys_sink_lane(
        normalized_amino=normalized_amino,
        polyphenol_active=polyphenol_active,
    )
    family_14_lane = _build_amino_dicarbonyl_source_lane(
        normalized_amino=normalized_amino,
        core_support_score=float(core_maillard_lane.get("core_support_score", 0.0)),
    )
    family_15_lane = _build_sugar_diversion_lane(
        normalized_sugars=normalized_sugars,
        lipids=lipids,
    )
    family_16_lane = _build_trapping_burden_modifier_lane(
        normalized_additives=normalized_additives,
        protein_label=protein_label,
    )
    family_lane_summary = {
        row["slr_family"]: row
        for row in [
            core_maillard_lane,
            lipid_crosstalk_lane,
            thiamine_support_lane,
            nucleotide_support_lane,
            sulfur_peptide_support_lane,
            matrix_scope_lane,
            donor_hierarchy_lane,
            caramelization_lane,
            fermentation_pretreatment_lane,
            off_note_guardrail_lane,
            family_11_lane,
            family_13_lane,
            family_14_lane,
            family_15_lane,
            family_16_lane,
        ]
    }
    family_lane_adjustments = _build_family_lane_adjustments(family_lane_summary)
    family_prior_bundle = summarize_family_prior_bundle(
        protein_type=protein_label,
        process_state="heated_matrix",
    )
    family_state_markers = _build_family_state_markers(
        thiamine_metadata=thiamine_metadata,
        thiamine_fraction_estimate=thiamine_fraction_estimate,
        thiamine_mode=thiamine_mode,
        thiamine_support_lane=thiamine_support_lane,
        nucleotide_support_lane=nucleotide_support_lane,
        fermentation_pretreatment_lane=fermentation_pretreatment_lane,
        caramelization_lane=caramelization_lane,
    )
    active_family_lanes = [
        slr_family
        for slr_family, row in sorted(family_lane_summary.items(), key=lambda item: _family_sort_key(item[0]))
        if bool(row.get("active", False))
    ]
    active_family_lane_names = [
        str(family_lane_summary[slr_family].get("display_name", slr_family))
        for slr_family in active_family_lanes
    ]

    return {
        "strecker_balance_score": float(strecker_balance_score),
        "strecker_gap_penalty": float(strecker_gap_penalty),
        "strecker_signals_ppb": {
            "methional": float(methional_signal),
            "2_methylbutanal": float(two_methylbutanal_signal),
            "3_methylbutanal": float(three_methylbutanal_signal),
            "phenylacetaldehyde": float(phenylacetaldehyde_signal),
            "benzaldehyde": float(benzaldehyde_signal),
        },
        "strecker_reference_targets_ppb": {
            "methional": float(methional_target),
            "2_methylbutanal": float(two_methylbutanal_target),
            "3_methylbutanal": float(three_methylbutanal_target),
            "phenylacetaldehyde": float(phenylacetaldehyde_target),
            "benzaldehyde": float(benzaldehyde_target),
        },
        "flavor_reference_policy": build_flavor_reference_policy_summary(),
        "pyrazine_signal_ppb": float(pyrazine_signal),
        "pyrazine_propensity": float(pyrazine_propensity),
        "pyrazine_burden": float(pyrazine_burden),
        "pyrazine_penalty": float(pyrazine_penalty),
        "pyrazine_drivers": {
            "pH_factor": float(ph_factor),
            "sugar_bias": float(sugar_bias),
            "peptide_bias": float(peptide_bias),
        },
        "furanone_expected": furanone_expected,
        "furanone_observed": furanone_observed,
        "furanone_support_score": float(furanone_support_score),
        "furanone_penalty": float(furanone_penalty),
        "furanone_missing": missing_furanones,
        "furanone_confidence": str(furanone_priors.get("confidence_tier", "low")),
        "thiamine_pathway_active": bool(thiamine_active),
        "thiamine_availability_source": str(thiamine_metadata["source"]),
        "thiamine_availability_explicit": bool(thiamine_metadata["explicit"]),
        "thiamine_inferred_from_inputs": bool(thiamine_metadata["inferred_from_inputs"]),
        "thiamine_mft_fraction_estimate": float(thiamine_fraction_estimate),
        "thiamine_provenance_mode": thiamine_mode,
        "lincoln_crosstalk_prior": lincoln_crosstalk_prior,
        "family_upstream_contract": upstream_contract,
        "family_lane_summary": family_lane_summary,
        "family_lane_adjustments": family_lane_adjustments,
        "family_prior_bundle": family_prior_bundle,
        "family_state_markers": family_state_markers,
        "active_family_lanes": active_family_lanes,
        "active_family_lane_names": active_family_lane_names,
    }