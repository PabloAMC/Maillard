from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

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


def describe_retention_runtime(
    compound_name: str,
    *,
    protein_type: Optional[str],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
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

    if protein.startswith("pea") and normalized == "hexanal":
        retention_mode = "direct_binding_anchor"
        sources.append("JAFC 10.1021/acs.jafc.3c05991 direct PPI hexanal binding")

    if protein.startswith("soy") and normalized == "hexanal":
        release_factor = 1.0 + 0.55 * _sigmoid(temperature, 58.0, 10.0)
        temporal_attenuation_factor = _log_floor_decay(duration, reference=8.0, slope=0.18, floor=0.62) if temperature >= 90.0 else 1.0
        dynamic_retention_factor = max(0.85, min(1.65, release_factor * temporal_attenuation_factor))
        retention_mode = "reversible_release_plus_temporal_attenuation"
        sources.extend([
            "Ince 2024 reversible non-covalent soy hexanal binding",
            "JAFC 10.1021/acs.jafc.3c05991 direct hexanal binding",
        ])
        if temperature >= 90.0:
            sources.append("Xu 2023 heated SPI temporal off-flavour attenuation prior")
    elif protein.startswith("soy") and normalized in {"2 pentylfuran", "2 pentyl furan"}:
        temporal_attenuation_factor = _log_floor_decay(duration, reference=8.0, slope=0.24, floor=0.55) if temperature >= 90.0 else 1.0
        dynamic_retention_factor = temporal_attenuation_factor
        retention_mode = "temporal_attenuation"
        sources.extend([
            "Xu 2023 heated SPI temporal off-flavour attenuation prior",
            "Shu 2024 heated SPI 2-pentylfuran below-detection endpoint",
        ])
    elif protein.startswith("soy") and any(token in normalized for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
        retention_mode = "sulfur_proxy_reference"
        sources.append("Zhang 2026 SPI sulfur proxy retention prior")

    if not sources and process_state == "heated_matrix" and protein.startswith("soy"):
        sources.append("heated soy process-state carryover")

    return {
        "retention_runtime_mode": retention_mode,
        "retention_reference_sources": sources,
        "reversible_release_factor": float(release_factor),
        "temporal_attenuation_factor": float(temporal_attenuation_factor),
        "dynamic_retention_factor": float(dynamic_retention_factor),
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


def build_flavor_axis_summary(
    *,
    projection_metadata: ProjectionMetadataMap,
    sugars: Optional[List[str]] = None,
    amino_acids: Optional[List[str]] = None,
    additives: Optional[List[str]] = None,
    lipids: Optional[List[str]] = None,
    protein_type: Optional[str] = None,
    pH: Optional[float] = None,
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
    phenylacetaldehyde_signal = _compound_signal(signal_map, ["Phenylacetaldehyde", "phenylacetaldehyde"])

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
    phenylacetaldehyde_target = _reference_panel_value(
        "strecker_reference_anchors",
        "hernandez_2023_phenylacetaldehyde_panel",
        "regular_ground_beef",
    ) or 11.85

    weighted_strecker = (
        2.0 * min(two_methylbutanal_signal / max(two_methylbutanal_target, 1.0e-6), 1.5)
        + 1.0 * min(methional_signal / max(methional_target, 1.0e-6), 1.5)
        + 0.6 * min(phenylacetaldehyde_signal / max(phenylacetaldehyde_target, 1.0e-6), 1.5)
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

    thiamine_priors = get_thiamine_priors()
    thiamine_active = any("thiamine" in value for value in normalized_additives + normalized_amino + normalized_sugars)
    thiamine_fraction_estimate = 0.0
    thiamine_mode = "inactive"
    if thiamine_active:
        if pentose_active:
            thiamine_fraction_estimate = float(thiamine_priors.get("mixed_pentose_fraction", 0.5) or 0.5)
            thiamine_mode = "mixed_thiamine_plus_pentose"
        else:
            thiamine_fraction_estimate = float(thiamine_priors.get("thiamine_only_fraction", 1.0) or 1.0)
            thiamine_mode = "thiamine_only"

    polyphenol_active = any(
        token in value
        for value in normalized_additives
        for token in ["catechin", "tannic", "green tea", "grape seed", "polyphenol"]
    )
    lincoln_crosstalk_prior = {
        "active": bool(polyphenol_active and lipids and any("glucose" in sugar for sugar in normalized_sugars)),
        "summary": (
            "Lincoln 2025 qualitative prior active: glucose plus lipids plus polyphenols can suppress Strecker aldehydes via dicarbonyl competition."
            if polyphenol_active and lipids and any("glucose" in sugar for sugar in normalized_sugars)
            else "Lincoln 2025 qualitative prior inactive for this formulation context."
        ),
    }

    return {
        "strecker_balance_score": float(strecker_balance_score),
        "strecker_gap_penalty": float(strecker_gap_penalty),
        "strecker_signals_ppb": {
            "methional": float(methional_signal),
            "2_methylbutanal": float(two_methylbutanal_signal),
            "phenylacetaldehyde": float(phenylacetaldehyde_signal),
        },
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
        "furanone_confidence": str(furanone_priors.get("confidence_tier", "low")),
        "thiamine_pathway_active": bool(thiamine_active),
        "thiamine_mft_fraction_estimate": float(thiamine_fraction_estimate),
        "thiamine_provenance_mode": thiamine_mode,
        "lincoln_crosstalk_prior": lincoln_crosstalk_prior,
    }