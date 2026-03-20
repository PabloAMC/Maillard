from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src.projection_metadata import ProjectionMetadataMap


ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"


def _load_json_payload(file_name: str) -> dict[str, Any]:
    payload_path = DATA_LIT_DIR / file_name
    if not payload_path.exists():
        return {}
    with open(payload_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


FLAVOR_REFERENCE_PAYLOADS = _load_json_payload("flavor_reference_payloads.json")
RETENTION_REFERENCE_PAYLOADS = _load_json_payload("retention_reference_payloads.json")
COMPUTATIONAL_PRIORS_PAYLOAD = _load_json_payload("computational_priors.json")


def _normalize_name(name: str) -> str:
    normalized = str(name).lower().replace("_", " ").replace("-", " ")
    return " ".join(normalized.split())


def _sigmoid(value: float, center: float, width: float) -> float:
    width = max(float(width), 1.0e-6)
    exponent = max(-60.0, min(60.0, -(float(value) - center) / width))
    return 1.0 / (1.0 + math.exp(exponent))


def _log_floor_decay(value: Optional[float], *, reference: float, slope: float, floor: float) -> float:
    if value is None or value <= 0.0:
        return 1.0
    factor = 1.0 - slope * math.log1p(float(value) / max(reference, 1.0e-6))
    return max(float(floor), min(1.0, float(factor)))


def get_pyrazine_control_priors() -> dict[str, Any]:
    for entry in COMPUTATIONAL_PRIORS_PAYLOAD.get("pyrazine_control_priors", []):
        if str(entry.get("id", "")).strip() == "slr_pyrazine_control_surface_v1":
            return dict(entry)
    return {}


def get_furanone_priors() -> dict[str, Any]:
    for entry in COMPUTATIONAL_PRIORS_PAYLOAD.get("furanone_priors", []):
        if str(entry.get("id", "")).strip() == "blank_fay_1996_furanone_expectation_v1":
            return dict(entry)
    return {}


def get_thiamine_priors() -> dict[str, Any]:
    for entry in COMPUTATIONAL_PRIORS_PAYLOAD.get("thiamine_pathway_priors", []):
        if str(entry.get("id", "")).strip() == "cerny_2007_thiamine_split_v1":
            return dict(entry)
    return {}


def get_strecker_crosstalk_priors() -> dict[str, Any]:
    for entry in COMPUTATIONAL_PRIORS_PAYLOAD.get("strecker_crosstalk_priors", []):
        if str(entry.get("id", "")).strip() == "lincoln_2025_polyphenol_crosstalk_v1":
            return dict(entry)
    return {}


_DEFAULT_FLAVOR_PIPELINE_ROLE_BY_BENCHMARK = {
    "reference_anchor": "secondary_marker",
    "directional_comparison_anchor": "diagnostic_marker",
    "pbma_counterexample": "reference_only",
    "low_confidence_mechanistic_anchor": "optimization_constraint",
}


def _iter_flavor_entries() -> Iterable[tuple[str, Dict[str, Any]]]:
    for section_name, entries in FLAVOR_REFERENCE_PAYLOADS.items():
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if isinstance(entry, Mapping):
                yield str(section_name), dict(entry)


def _flavor_entry_by_id(entry_id: str) -> Dict[str, Any]:
    requested = str(entry_id).strip()
    for _, entry in _iter_flavor_entries():
        if str(entry.get("id", "")).strip() == requested:
            return entry
    return {}


def get_flavor_reference_pipeline_role(entry_id: str) -> str:
    entry = _flavor_entry_by_id(entry_id)
    explicit = str(entry.get("pipeline_role", "")).strip().lower()
    if explicit:
        return explicit
    benchmark_role = str(entry.get("benchmark_role", "")).strip().lower()
    return _DEFAULT_FLAVOR_PIPELINE_ROLE_BY_BENCHMARK.get(benchmark_role, "reference_only")


def build_flavor_reference_policy_summary() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for section_name, entry in _iter_flavor_entries():
        rows.append(
            {
                "id": str(entry.get("id", "unknown")),
                "section": section_name,
                "compound": str(entry.get("compound", "unknown")),
                "benchmark_role": str(entry.get("benchmark_role", "unknown")),
                "pipeline_role": get_flavor_reference_pipeline_role(str(entry.get("id", "unknown"))),
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


def _iter_retention_entries() -> Iterable[tuple[str, str, Dict[str, Any]]]:
    for section_name, section_payload in RETENTION_REFERENCE_PAYLOADS.items():
        if not isinstance(section_payload, Mapping):
            continue
        for matrix_family, entries in section_payload.items():
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if isinstance(entry, Mapping):
                    yield str(section_name), str(matrix_family), dict(entry)


def _retention_entry_by_id(entry_id: str) -> Dict[str, Any]:
    requested = str(entry_id).strip()
    for _, _, entry in _iter_retention_entries():
        if str(entry.get("id", "")).strip() == requested:
            return entry
    return {}


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
    normalized = _normalize_name(compound_name)
    protein = str(protein_type or "")

    candidate_ids: List[str] = []
    if protein.startswith("pea") and normalized == "hexanal":
        candidate_ids.append("karolkowski_2021_ppi_hexanal_ph_release")
    elif protein.startswith("pea") and normalized in {"2 pentylfuran", "2 pentyl furan"}:
        candidate_ids.append("karolkowski_2021_ppi_2_pentylfuran_native_panel")

    for entry_id in candidate_ids:
        entry = _retention_entry_by_id(entry_id)
        surrogate = entry.get("runtime_surrogate", {})
        if isinstance(surrogate, Mapping) and surrogate.get("type") == "ph_release_modifier":
            return {
                **dict(surrogate),
                "source": _retention_source_label(entry),
                "entry_id": str(entry.get("id", entry_id)),
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


def _reference_panel_value(payload_section: str, entry_id: str, panel_key: str) -> float:
    for entry in FLAVOR_REFERENCE_PAYLOADS.get(payload_section, []):
        if str(entry.get("id", "")) != entry_id:
            continue
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


def build_flavor_axis_summary(
    *,
    projection_metadata: ProjectionMetadataMap,
    sugars: Optional[List[str]] = None,
    amino_acids: Optional[List[str]] = None,
    additives: Optional[List[str]] = None,
    lipids: Optional[List[str]] = None,
    protein_type: Optional[str] = None,
    pH: Optional[float] = None,
    thiamine_availability: Any = None,
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
    thiamine_metadata = _resolve_thiamine_availability(
        thiamine_availability,
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        protein_label=protein_label,
    )
    thiamine_active = bool(thiamine_metadata["available"])
    thiamine_fraction_estimate = 0.0
    thiamine_mode = "inactive"
    if thiamine_active:
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
    }