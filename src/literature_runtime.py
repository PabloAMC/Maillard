from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from src import data_paths
from src import data_access
from src.family_ingestion_plan import load_family_ingestion_plan
from src.kokumi_scoring import build_kokumi_support_profile
from src.literature_family_registry import (
    iter_benchmark_intake_entries,
    iter_computational_prior_entries as iter_family_computational_prior_entries,
    iter_flavor_reference_entries as iter_family_flavor_reference_entries,
    iter_process_gap_entries as iter_family_process_gap_entries,
    iter_retention_reference_entries as iter_family_retention_reference_entries,
)
from src.matrix_prior_registry import summarize_family_prior_bundle
from src.matrix_targets import get_compound_panel_entry, iter_target_panel_entries
from src.lipid_oxidation import summarize_lipid_runtime_split
from src.projection_metadata import ProjectionMetadataMap
from src.safety import (
    MG_PER_KG_TO_PPB,
    get_safety_reference_payload,
    get_safety_reference_range,
    mg_per_100g_protein_to_ppb,
    predict_acrylamide,
    predict_cel,
    predict_cml,
    predict_furosine,
)


def _load_json_payload(payload_path: Path) -> dict[str, Any]:
    # Strict since 2026-09-01: a missing registry used to load as {} and the runtime
    # carried on without it.
    return data_access.load_json(payload_path)


FLAVOR_REFERENCE_PAYLOADS = _load_json_payload(data_paths.FLAVOR_REFERENCE_PAYLOADS)
RETENTION_REFERENCE_PAYLOADS = _load_json_payload(data_paths.RETENTION_REFERENCE_PAYLOADS)
COMPUTATIONAL_PRIORS_PAYLOAD = _load_json_payload(data_paths.COMPUTATIONAL_PRIORS)
PROCESS_STATE_CALIBRATIONS_PAYLOAD = _load_json_payload(data_paths.PROCESS_STATE_CALIBRATIONS)
PROTEIN_SOURCE_REGISTRY_PAYLOAD = _load_json_payload(data_paths.PROTEIN_SOURCE_REGISTRY)
FAMILY_INGESTION_PLAN_PAYLOAD = load_family_ingestion_plan()


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
    "nucleotide_pathway": {
        "section_name": "nucleotide_pathway_priors",
        "default_id": "matoba_1988_nucleotide_hydrolysis_v1",
        "slr_family": "04",
    },
    "sulfur_peptide_support": {
        "section_name": "sulfur_peptide_priors",
        "default_id": "wang_xu_glutathione_peptide_support_v1",
        "slr_family": "05",
    },
}

_PRIMARY_CONTRACT_TAGS_BY_FAMILY = {
    "03": {"core_meaty", "mft", "extrusion_survival", "beef_realistic"},
    "04": {"nucleotide_support", "ribose_delivery", "umami_support", "euc_anchor", "rate_limiting_step", "process_state_calibration", "sensory_threshold"},
    "05": {"core_meaty", "matrix_intake"},
}


_PROTEIN_SOURCE_ALIAS_MAP = {
    "free": "",
    "myco": "mycoprotein",
    "pea_iso": "pea_isolate",
    "pea_conc": "pea_concentrate",
    "soy_iso": "soy_isolate",
    "soy_conc": "soy_concentrate",
}


_PROTEIN_SOURCE_PROFILES = {
    str(entry.get("source_id", "")).strip().lower(): dict(entry)
    for entry in PROTEIN_SOURCE_REGISTRY_PAYLOAD.get("sources", [])
    if isinstance(entry, Mapping) and str(entry.get("source_id", "")).strip()
}

# 2026-08-27 (Wave T3, finding T1-01). protein_source_registry.json is a MOCK. Its
# own description has always said so; nothing read that sentence, so the mock status
# never reached either stderr or any output payload while the numbers underneath it
# drove `matrix_uncertainty_factor` (see `_alternative_matrix_lane` below) and the
# meaty-potential score. Warned at load, in the shape of the family-12 molar-ratio
# unit warning in `_resolve_concentration_unit`: state the defect, name what depends
# on it, rescale nothing, invent nothing.
PROTEIN_SOURCE_REGISTRY_UNSOURCED = (
    str(PROTEIN_SOURCE_REGISTRY_PAYLOAD.get("source_status", "")).strip()
    == "no_verifiable_source"
)
if PROTEIN_SOURCE_REGISTRY_UNSOURCED:
    warnings.warn(
        f"src.literature_runtime: {data_paths.rel(data_paths.PROTEIN_SOURCE_REGISTRY)} declares "
        "source_status='no_verifiable_source' -- every value in it is a MOCKED "
        "placeholder whose only declared upstream is the LLM digest "
        "data/Gemini_Deep_Research/06_alternative_proteins.md, which is not provenance. "
        "It is nonetheless LIVE: hydrolysate_observability_bias, off_note_penalty and "
        "lox_activity_flag enter matrix_uncertainty_factor directly, and "
        "meaty_potential_multiplier drives the meaty-potential score. The plant-source "
        "DIFFERENTIATION these values encode is not evidence. Values are NOT substituted.",
        RuntimeWarning,
        stacklevel=2,
    )

#: Surfaced verbatim on the family-06 lane payload so the mock status travels with the
#: number it contaminates instead of living only in a warning nobody sees.
PROTEIN_SOURCE_PROVENANCE = {
    "registry": data_paths.rel(data_paths.PROTEIN_SOURCE_REGISTRY),
    "source_status": str(PROTEIN_SOURCE_REGISTRY_PAYLOAD.get("source_status", "") or "unstated"),
    "value_basis": str(PROTEIN_SOURCE_REGISTRY_PAYLOAD.get("value_basis", "") or "unstated"),
    "declared_upstream": "data/Gemini_Deep_Research/06_alternative_proteins.md (LLM digest)",
    "affects": [
        "matrix_uncertainty_factor",
        "matrix_source_support_score",
        "process_state_transfer_confidence",
    ],
    "warning": (
        "MOCKED VALUES. Source-to-source differences in these outputs are not evidence; "
        "they reproduce an ordering someone wrote down. Wave T3 (2026-08-27), finding T1-01."
    ),
}


def _normalize_name(name: str) -> str:
    normalized = str(name).lower().replace("_", " ").replace("-", " ")
    return " ".join(normalized.split())


#: Tokens whose legitimate matches are morphological variants ("ferment" ->
#: "fermentation"/"fermented", "culture" -> "cultured"). They match at a word
#: START; every other token must match a whole word. Plain substring matching
#: used to fire on "important"/"imported" for "imp" and on "asnase"
#: (asparaginase) for "asn".
_STEM_MATCH_TOKENS = {"ferment", "culture", "hydroly", "cultur"}


def _token_matches(value: str, token: str) -> bool:
    """Word-boundary token matching over an already-normalized name."""
    target = str(token).strip().lower().replace("_", " ").replace("-", " ")
    if not target:
        return False
    text = str(value).strip().lower()
    if " " in target:
        return target in text
    words = text.split()
    if target in _STEM_MATCH_TOKENS:
        return any(word.startswith(target) for word in words)
    return any(word == target or word == target + "s" for word in words)


def _any_token_matches(value: str, tokens: Sequence[str]) -> bool:
    return any(_token_matches(value, token) for token in tokens)


#: CANONICAL UNIT for every `molar_ratios` / `effective_molar_ratios` mapping in
#: the repo: millimolar (mM). `src.safety.predict_acrylamide` and the AGE
#: predictors consume these values as mM, and three different authoring
#: conventions were found in-repo (mM, mol/L, and unitless 1:1:1-style ratios),
#: which move the predicted acrylamide by ~3300x. Payloads may declare their own
#: unit with a `concentration_unit` entry; values are NEVER silently rescaled.
MOLAR_RATIO_CANONICAL_UNIT = "mM"

#: Reserved keys inside a molar_ratios mapping that carry metadata, not a
#: concentration.
MOLAR_RATIO_METADATA_KEYS = {"concentration_unit", "__concentration_unit__", "units", "unit"}

#: Heuristic: if every declared concentration is <= this and no unit was
#: declared, the mapping was almost certainly authored as unitless ratios.
UNITLESS_RATIO_SUSPICION_MAX = 5.0


def _extract_declared_concentration_unit(mapping: Optional[Mapping[str, Any]]) -> str:
    for key, value in (mapping or {}).items():
        if str(key).strip().lower() in MOLAR_RATIO_METADATA_KEYS and isinstance(value, str):
            return str(value).strip()
    return ""


def _numeric_ratio_items(mapping: Optional[Mapping[str, Any]]) -> Dict[str, float]:
    numeric: Dict[str, float] = {}
    for key, value in (mapping or {}).items():
        if str(key).strip().lower() in MOLAR_RATIO_METADATA_KEYS:
            continue
        try:
            numeric[str(key)] = float(value)
        except (TypeError, ValueError):
            continue
    return numeric


def _assess_concentration_unit(
    mapping: Optional[Mapping[str, Any]],
    *,
    context: str,
) -> Dict[str, Any]:
    """Resolve the declared unit of a molar_ratios mapping and warn on ambiguity.

    Returns the declared unit (empty when none), the unit actually assumed, and
    whether the values look like unitless ratios. Nothing is rescaled: the caller
    keeps the numbers it was given, it just learns that they are suspect.
    """
    declared = _extract_declared_concentration_unit(mapping)
    numeric = _numeric_ratio_items(mapping)
    values = [value for value in numeric.values() if value > 0.0]
    looks_unitless = bool(values) and not declared and max(values) <= UNITLESS_RATIO_SUSPICION_MAX
    if looks_unitless:
        # The message deliberately names only the KEYS, not the values, so the
        # default warning filter collapses a whole optimizer sweep into one
        # emission per distinct ingredient set instead of one per trial.
        warnings.warn(
            f"{context}: molar_ratios for {sorted(numeric)} carry NO declared "
            f"concentration_unit and every value is <= {UNITLESS_RATIO_SUSPICION_MAX}; they are "
            f"being consumed as {MOLAR_RATIO_CANONICAL_UNIT} (the canonical unit) but look like "
            "unitless ratios. Acrylamide and AGE predictions scale with these values - declare "
            "concentration_unit in the payload. Values are NOT rescaled.",
            RuntimeWarning,
            stacklevel=3,
        )
    if declared and declared.strip().lower() not in {"mm", "millimolar"}:
        warnings.warn(
            f"{context}: molar_ratios declare concentration_unit={declared!r}, but the safety "
            f"lane consumes {MOLAR_RATIO_CANONICAL_UNIT}. Values are NOT rescaled; convert "
            "upstream.",
            RuntimeWarning,
            stacklevel=3,
        )
    return {
        "declared_concentration_unit": declared,
        "assumed_concentration_unit": MOLAR_RATIO_CANONICAL_UNIT,
        "looks_like_unitless_ratios": bool(looks_unitless),
    }


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


def _iter_process_state_calibration_entries() -> Iterable[Dict[str, Any]]:
    entries = PROCESS_STATE_CALIBRATIONS_PAYLOAD.get("entries", []) if isinstance(PROCESS_STATE_CALIBRATIONS_PAYLOAD, Mapping) else []
    for entry in entries:
        if isinstance(entry, Mapping):
            yield dict(entry)


def _process_state_calibration_entry(entry_id: str) -> Dict[str, Any]:
    for entry in _iter_process_state_calibration_entries():
        if str(entry.get("id", "")).strip() == str(entry_id).strip():
            return entry
    return {}


def _computational_prior_entry(entry_id: str) -> Dict[str, Any]:
    rows = query_family_runtime_priors(entry_id=entry_id)
    return rows[0] if rows else {}


def _runtime_source_label(entry: Mapping[str, Any]) -> str:
    return str(
        entry.get("source_citation")
        or entry.get("source")
        or entry.get("citation")
        or entry.get("id")
        or "unknown runtime reference"
    )


def _resolve_protein_source_id(protein_label: str) -> str:
    normalized = str(protein_label or "").strip().lower()
    if not normalized:
        return ""
    return _PROTEIN_SOURCE_ALIAS_MAP.get(normalized, normalized)


def _protein_source_profile(source_id: str) -> Dict[str, Any]:
    normalized = str(source_id or "").strip().lower()
    if not normalized:
        return {}
    return dict(_PROTEIN_SOURCE_PROFILES.get(normalized, {}))


def _interpolate_profile(points: Iterable[float], values: Iterable[float], query_value: Optional[float]) -> Optional[float]:
    if query_value is None:
        return None
    pairs = sorted((float(point), float(value)) for point, value in zip(points, values))
    if not pairs:
        return None
    query = float(query_value)
    if query <= pairs[0][0]:
        return float(pairs[0][1])
    for (left_point, left_value), (right_point, right_value) in zip(pairs, pairs[1:]):
        if query <= right_point:
            span = max(right_point - left_point, 1.0e-6)
            weight = (query - left_point) / span
            return float(left_value + (right_value - left_value) * weight)
    return float(pairs[-1][1])


def _interpolate_numeric_mapping(values_by_point: Mapping[str, Any], query_value: Optional[float]) -> Optional[float]:
    pairs: List[tuple[float, float]] = []
    for point, value in values_by_point.items():
        try:
            pairs.append((float(point), float(value)))
        except (TypeError, ValueError):
            continue
    if not pairs:
        return None
    return _interpolate_profile([point for point, _ in pairs], [value for _, value in pairs], query_value)


def _estimate_thiamine_extrusion_survival(
    thiamine_priors: Mapping[str, Any],
    *,
    process_state: Optional[str],
    temperature_celsius: Optional[float],
    water_activity: Optional[float],
) -> Dict[str, Any]:
    calibration_entry = _process_state_calibration_entry(str(thiamine_priors.get("extrusion_retention_reference_id", "")))
    if process_state not in {"aqueous_pre_extrusion_model", "extrusion_structured"}:
        return {
            "extrusion_survival_factor": 1.0,
            "extrusion_reference_id": str(calibration_entry.get("id", "")),
            "extrusion_reference_source": str(calibration_entry.get("source_citation", "")),
            "extrusion_reference_active": False,
        }

    anchors = calibration_entry.get("numeric_anchors", {}) if isinstance(calibration_entry.get("numeric_anchors", {}), Mapping) else {}
    baseline = _interpolate_profile(
        [140.0, 160.0, 180.0],
        [
            float(anchors.get("retention_pct_140c", 68.0)) / 100.0,
            float(anchors.get("retention_pct_160c", 31.0)) / 100.0,
            float(anchors.get("retention_pct_180c", 8.0)) / 100.0,
        ],
        temperature_celsius if temperature_celsius is not None else 160.0,
    )
    survival = float(baseline if baseline is not None else 1.0)
    aw = None if water_activity is None else max(0.0, min(1.0, float(water_activity)))
    if aw is not None:
        moisture_penalty = 1.0 - 0.5 * max(0.0, min(1.0, (aw - 0.65) / 0.10))
        survival *= moisture_penalty
    survival = max(0.04, min(1.0, survival))
    return {
        "extrusion_survival_factor": float(survival),
        "extrusion_reference_id": str(calibration_entry.get("id", "")),
        "extrusion_reference_source": str(calibration_entry.get("source_citation", "")),
        "extrusion_reference_active": True,
        "extrusion_aw_surrogate": aw,
    }


def _build_thiamine_calibration_context(
    *,
    thiamine_priors: Mapping[str, Any],
    baseline_fraction: float,
    thiamine_mode: str,
    molar_ratios: Optional[Mapping[str, Any]],
    pH: Optional[float],
    process_state: Optional[str],
    temperature_celsius: Optional[float],
    water_activity: Optional[float],
    sulfur_partner_active: bool,
) -> Dict[str, Any]:
    benchmark_rows = query_benchmark_intake_entries(family="03", primary_only=True)
    aw_support_prior = _computational_prior_entry("arabshahi_1988_aw_dependent_thiamine_ea_v1")
    alone_points = thiamine_priors.get("thiamine_alone_ph_points", []) or [4.0, 5.0, 6.0, 7.0, 8.0]
    alone_yields = thiamine_priors.get("thiamine_alone_mft_yield_ug_per_g", []) or [0.02, 0.14, 0.58, 0.42, 0.19]
    reference_ph = float(thiamine_priors.get("reference_ph", 6.0) or 6.0)
    reference_yield = float(thiamine_priors.get("reference_yield_ug_per_g", 0.58) or 0.58)
    alone_yield = _interpolate_profile(alone_points, alone_yields, pH if pH is not None else reference_ph)
    pH_factor = 1.0 if alone_yield is None else max(0.03, min(1.0, float(alone_yield) / max(reference_yield, 1.0e-6)))
    yield_mode = "thiamine_alone_reference"
    reference_yield_value = float(alone_yield if alone_yield is not None else reference_yield)

    mixed_system_active = thiamine_mode == "mixed_thiamine_plus_pentose" and sulfur_partner_active
    mixed_peak = float(thiamine_priors.get("mixed_system_synergy_factor_peak", 4.3) or 4.3)
    synergy_floor = float(thiamine_priors.get("beef_realistic_synergy_factor", 2.64) or 2.64)
    synergy_factor = 1.0
    precursor_loading_mM = 0.0
    for key, value in (molar_ratios or {}).items():
        normalized_key = _normalize_name(str(key))
        try:
            numeric_value = max(0.0, float(value))
        except (TypeError, ValueError):
            continue
        if "thiamine" in normalized_key or "vitamin b1" in normalized_key:
            precursor_loading_mM += numeric_value
        elif any(token in normalized_key for token in ["cysteine", "ribose", "xylose", "arabinose"]):
            precursor_loading_mM += numeric_value
    if mixed_system_active:
        if pH is None:
            synergy_factor = synergy_floor
        else:
            peak_weight = math.exp(-((float(pH) - 5.75) / 1.1) ** 2)
            synergy_factor = synergy_floor + (mixed_peak - synergy_floor) * peak_weight
        optimal_ph = float(thiamine_priors.get("mixed_system_optimal_ph", 5.5) or 5.5)
        reference_mixed_ph = float(thiamine_priors.get("mixed_system_reference_ph", 6.0) or 6.0)
        optimal_yield = float(thiamine_priors.get("mixed_system_optimal_yield_ug_per_g", 3.11) or 3.11)
        reference_mixed_yield = float(thiamine_priors.get("mixed_system_reference_yield_ug_per_g", 2.47) or 2.47)
        if pH is not None and optimal_ph <= float(pH) <= reference_mixed_ph:
            mixed_factor = _interpolate_profile(
                [optimal_ph, reference_mixed_ph],
                [optimal_yield / max(reference_mixed_yield, 1.0e-6), 1.0],
                float(pH),
            )
            if mixed_factor is not None:
                pH_factor = max(0.05, min(1.3, float(mixed_factor)))
                yield_mode = "mixed_system_optimal_window"
                reference_yield_value = float(_interpolate_profile([optimal_ph, reference_mixed_ph], [optimal_yield, reference_mixed_yield], float(pH)) or reference_mixed_yield)

    reference_loading_mM = float(thiamine_priors.get("thiamine_only_reference_total_precursor_mM", 10.0) or 10.0)
    if mixed_system_active:
        if yield_mode == "mixed_system_optimal_window":
            reference_loading_mM = float(thiamine_priors.get("mixed_system_optimal_total_precursor_mM", 30.0) or 30.0)
        else:
            reference_loading_mM = float(thiamine_priors.get("mixed_system_reference_total_precursor_mM", 30.0) or 30.0)
    precursor_loading_mM = float(precursor_loading_mM or reference_loading_mM)
    beef_reference_yield = float(thiamine_priors.get("beef_reference_yield_ug_per_g", 3.14) or 3.14)
    beef_reference_loading = float(thiamine_priors.get("beef_reference_total_precursor_mM", 2.1) or 2.1)
    dilute_loading_uplift_factor = 1.0
    if mixed_system_active:
        dilute_loading_uplift_ceiling = float(
            thiamine_priors.get("mixed_system_dilute_loading_uplift_ceiling", 1.0) or 1.0
        )
        if dilute_loading_uplift_ceiling > 1.0 and reference_loading_mM > beef_reference_loading:
            dilute_loading_uplift_factor = float(
                _interpolate_profile(
                    [beef_reference_loading, reference_loading_mM],
                    [dilute_loading_uplift_ceiling, 1.0],
                    max(beef_reference_loading, min(reference_loading_mM, precursor_loading_mM)),
                )
                or 1.0
            )
    baseline_efficiency = beef_reference_yield / max(beef_reference_loading, 1.0e-6)
    current_efficiency = float(reference_yield_value) / max(precursor_loading_mM, 1.0e-6)
    observable_efficiency_factor = current_efficiency / max(baseline_efficiency, 1.0e-6)
    observable_efficiency_factor = max(
        float(thiamine_priors.get("observable_efficiency_floor", 0.03) or 0.03),
        min(
            float(thiamine_priors.get("observable_efficiency_ceiling", 1.1) or 1.1),
            float(observable_efficiency_factor),
        ),
    )

    extrusion_context = _estimate_thiamine_extrusion_survival(
        thiamine_priors,
        process_state=process_state,
        temperature_celsius=temperature_celsius,
        water_activity=water_activity,
    )
    aw = None if water_activity is None else max(0.05, min(0.98, float(water_activity)))
    aw_support_active = False
    aw_dependent_ea = 0.0
    aw_modulation_factor = 1.0
    aw_transfer_weight = 1.0
    aw_reference = float(aw_support_prior.get("reference_aw", 0.65) or 0.65)
    aw_reference_ea = _interpolate_numeric_mapping(
        aw_support_prior.get("starch_ea_kcal_per_mol_by_aw", {}) if isinstance(aw_support_prior.get("starch_ea_kcal_per_mol_by_aw", {}), Mapping) else {},
        aw_reference,
    )
    if aw is not None:
        aw_dependent_ea_value = _interpolate_numeric_mapping(
            aw_support_prior.get("starch_ea_kcal_per_mol_by_aw", {}) if isinstance(aw_support_prior.get("starch_ea_kcal_per_mol_by_aw", {}), Mapping) else {},
            aw,
        )
        if aw_dependent_ea_value is not None and aw_reference_ea is not None:
            aw_support_active = bool(aw_support_prior)
            aw_dependent_ea = float(aw_dependent_ea_value)
            aw_modulation_factor = max(0.75, min(1.15, aw_dependent_ea / max(float(aw_reference_ea), 1.0e-6)))
            aw_transfer_weight = _thiamine_aw_transfer_weight(process_state)
            aw_modulation_factor = 1.0 - aw_transfer_weight * (1.0 - aw_modulation_factor)
    synergy_scaling = 1.0
    if mixed_system_active:
        synergy_scaling = 0.85 + 0.15 * min(1.0, synergy_factor / max(mixed_peak, 1.0e-6))
    effective_fraction = max(
        0.0,
        min(
            1.0,
            float(baseline_fraction)
            * float(pH_factor)
            * float(synergy_scaling)
            * float(dilute_loading_uplift_factor)
            * float(aw_modulation_factor)
            * float(extrusion_context.get("extrusion_survival_factor", 1.0)),
        ),
    )
    return {
        "baseline_fraction": float(baseline_fraction),
        "effective_fraction": float(effective_fraction),
        "pH_factor": float(pH_factor),
        "reference_yield_ug_per_g": float(reference_yield_value),
        "yield_mode": yield_mode,
        "synergy_factor": float(synergy_factor),
        "synergy_scaling": float(synergy_scaling),
        "precursor_loading_mM": float(precursor_loading_mM),
        "reference_loading_mM": float(reference_loading_mM),
        "dilute_loading_uplift_factor": float(dilute_loading_uplift_factor),
        "observable_efficiency_factor": float(observable_efficiency_factor),
        "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in benchmark_rows],
        "benchmark_anchor_citations": [str(row.get("citation", "unknown")) for row in benchmark_rows],
        "aw_reference_active": bool(aw_support_active),
        "aw_reference_id": str(aw_support_prior.get("id", "")) if aw_support_active else "",
        "aw_reference_source": str(aw_support_prior.get("source", "")) if aw_support_active else "",
        "aw_reference_matrix": str(aw_support_prior.get("reference_conditions", {}).get("matrix", "")) if aw_support_active and isinstance(aw_support_prior.get("reference_conditions", {}), Mapping) else "",
        "aw_reference_value": float(aw) if aw is not None else 0.0,
        "aw_reference_ea_kcal_per_mol": float(aw_dependent_ea),
        "aw_modulation_factor": float(aw_modulation_factor),
        "aw_transfer_weight": float(aw_transfer_weight),
        **extrusion_context,
    }


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

def _thiamine_aw_transfer_weight(process_state: Optional[str]) -> float:
    normalized = _normalize_name(process_state or "")
    if normalized in {_normalize_name("heated_matrix"), _normalize_name("extrusion_structured")}:
        return 1.0
    if normalized == _normalize_name("aqueous_pre_extrusion_model"):
        return 0.9
    if normalized == _normalize_name("intermediate_matrix"):
        return 0.8
    if normalized == _normalize_name("ambient_slurry"):
        return 0.7
    return 1.0


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
    default_id = str(lane_metadata.get("default_id", "")).strip()
    rows.sort(
        key=lambda row: (
            _family_sort_key(str(row.get("slr_family_source", row.get("family", {}).get("slr_family", "99")))),
            0 if default_id and str(row.get("id", "")).strip() == default_id else 1,
            str(row.get("id", "unknown")),
        )
    )
    return rows


def _matches_primary_contract_surface(row: Mapping[str, Any], target_slr: str) -> bool:
    required_tags = _PRIMARY_CONTRACT_TAGS_BY_FAMILY.get(str(target_slr).zfill(2), set())
    if not required_tags:
        return True
    row_tags = {
        str(tag).strip()
        for tag in (row.get("observable_panel_tags", []) or [])
        if str(tag).strip()
    }
    return bool(row_tags.intersection(required_tags))


def get_family_runtime_prior(
    *,
    runtime_lane: str,
    entry_id: Optional[str] = None,
) -> dict[str, Any]:
    lane_metadata = _RUNTIME_LANE_PRIOR_REGISTRY.get(str(runtime_lane).strip(), {})
    requested_id = str(entry_id or lane_metadata.get("default_id", "")).strip() or None
    rows = query_family_runtime_priors(runtime_lane=runtime_lane, entry_id=requested_id)
    return rows[0] if rows else {}


def query_benchmark_intake_entries(
    *,
    family: Optional[str] = None,
    entry_id: Optional[str] = None,
    matrix_family: Optional[str] = None,
    kind: Optional[str] = None,
    status: Optional[str] = None,
    process_state: Optional[str] = None,
    primary_only: bool = False,
) -> List[Dict[str, Any]]:
    descriptor = _family_runtime_descriptor(str(family or "")) if family is not None else {}
    target_family_id = str(descriptor.get("family_id", family or "")).strip()
    target_slr = str(descriptor.get("slr_family", family or "")).zfill(2) if family is not None else ""
    rows: List[Dict[str, Any]] = []
    for entry in iter_benchmark_intake_entries(family=family):
        row = dict(entry)
        if row.get("payload_role") == "retention_payload":
            continue
        if entry_id is not None and str(row.get("id", "")).strip() != str(entry_id).strip():
            continue
        if primary_only:
            is_primary = False
            if target_family_id and str(row.get("chemistry_family", "")).strip() == target_family_id:
                is_primary = True
            if target_slr and str(row.get("slr_family_source", "")).zfill(2) == target_slr:
                is_primary = True
            if not is_primary:
                continue
            if str(row.get("kind", "")).strip() == "directional_prior":
                continue
            if target_slr and not _matches_primary_contract_surface(row, target_slr):
                continue
        if matrix_family is not None and _normalize_name(str(row.get("matrix_family", ""))) != _normalize_name(matrix_family):
            continue
        if kind is not None and str(row.get("kind", "")).strip() != str(kind).strip():
            continue
        if status is not None and str(row.get("status", "")).strip() != str(status).strip():
            continue
        if not _process_state_matches(row, process_state):
            continue
        rows.append(row)
    rows.sort(
        key=lambda row: (
            _family_sort_key(str(row.get("slr_family_source", row.get("family_descriptor", {}).get("slr_family", "99")))),
            str(row.get("id", "unknown")),
        )
    )
    return rows


def query_dft_kinetic_priors(
    *,
    family: Optional[str] = None,
    reaction_key: Optional[str] = None,
    active_arrhenius_key: Optional[str] = None,
    process_state: Optional[str] = None,
) -> List[Dict[str, Any]]:
    descriptor = _family_runtime_descriptor(str(family or "")) if family is not None else {}
    target_family_id = str(descriptor.get("family_id", family or "")).strip()
    target_slr = str(descriptor.get("slr_family", family or "")).zfill(2) if family is not None else ""
    dft_payload = COMPUTATIONAL_PRIORS_PAYLOAD.get("dft_kinetic_priors", {})
    entries = dft_payload.get("entries", []) if isinstance(dft_payload, Mapping) else []
    rows: List[Dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, Mapping):
            continue
        row = dict(entry)
        row.setdefault("id", str(row.get("reaction_key", "unknown")))
        slr_family = str(row.get("slr_family", "")).zfill(2)
        family_descriptor = _family_runtime_descriptor(slr_family)
        if family_descriptor:
            row.setdefault("chemistry_family", str(family_descriptor.get("family_id", "unknown")))
            row["family_descriptor"] = {
                "slr_family": str(family_descriptor.get("slr_family", slr_family)).zfill(2),
                "family_id": str(family_descriptor.get("family_id", row.get("chemistry_family", "unknown"))),
                "display_name": str(family_descriptor.get("display_name", row.get("chemistry_family", "unknown"))),
                "strategic_posture": str(family_descriptor.get("strategic_posture", "unknown")),
            }
        if target_family_id and str(row.get("chemistry_family", "")).strip() != target_family_id and slr_family != target_slr:
            continue
        if reaction_key is not None and str(row.get("reaction_key", "")).strip() != str(reaction_key).strip():
            continue
        if active_arrhenius_key is not None and str(row.get("active_arrhenius_key", "")).strip() != str(active_arrhenius_key).strip():
            continue
        if not _process_state_matches(row, process_state):
            continue
        row["_section_name"] = "dft_kinetic_priors"
        rows.append(row)
    rows.sort(key=lambda row: (_family_sort_key(str(row.get("slr_family", "99"))), str(row.get("reaction_key", "unknown"))))
    return rows


def get_pyrazine_control_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="pyrazine_control")


def get_furanone_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="furanone_support")


def get_thiamine_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="thiamine_fragmentation")


def get_nucleotide_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="nucleotide_pathway")


def get_sulfur_peptide_priors() -> dict[str, Any]:
    return get_family_runtime_prior(runtime_lane="sulfur_peptide_support")


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
    runtime_prior_ids: List[str] = []
    process_state_calibration_ids: List[str] = []
    sulfur_binding_factor = 1.0
    partition_temperature_factor = 1.0

    partition_calibration = _process_state_calibration_entry("acs_jafc_3c05991_ppi_spi_partitioning")
    sulfur_binding_prior = _computational_prior_entry("acs_jafc_3c02618_mft_disulfide_trapping_v1")
    protein_binding_prior = _computational_prior_entry("acs_jafc_0c01925_protein_binding_hierarchy_v1")
    pea_sh_crosscheck = _process_state_calibration_entry("malia_2025_pea_free_sh_crosscheck")
    extrusion_disulfide_calibration = _process_state_calibration_entry("raman_sds_extrusion_disulfide_severity")

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
        if partition_calibration:
            partition_temperature_factor = 0.96 + 0.14 * _sigmoid(temperature, 52.0, 11.0)
            dynamic_retention_factor *= partition_temperature_factor
            sources.append(_runtime_source_label(partition_calibration))
            process_state_calibration_ids.append(str(partition_calibration.get("id", "unknown")))
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
        if partition_calibration:
            partition_temperature_factor = 0.94 + 0.18 * _sigmoid(temperature, 58.0, 10.0)
            dynamic_retention_factor *= partition_temperature_factor
            sources.append(_runtime_source_label(partition_calibration))
            process_state_calibration_ids.append(str(partition_calibration.get("id", "unknown")))
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

    sulfur_family = any(token in normalized for token in ["mft", "fft", "thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"])
    pyrazine_family = "pyrazine" in normalized
    if protein.startswith(("pea", "soy")) and (sulfur_family or pyrazine_family):
        thermal_binding_drive = _sigmoid(temperature, 102.0, 14.0)
        residence_drive = 0.0 if duration is None else min(1.0, 0.28 * math.log1p(max(duration, 0.0)))
        dryness = 0.0 if water_activity is None else max(0.0, min(1.0, (0.62 - max(0.05, min(0.98, float(water_activity)))) / 0.30))
        process_drive = 1.0 if process_state == "extrusion_structured" else 0.72 if process_state == "aqueous_pre_extrusion_model" else 0.46 if process_state == "heated_matrix" else 0.0
        disulfide_pressure = min(1.0, 0.60 * thermal_binding_drive + 0.22 * residence_drive + 0.18 * dryness)
        binding_pressure = min(1.0, disulfide_pressure * process_drive)
        hierarchy_bias = 0.22 if sulfur_family else 0.10
        sulfur_binding_factor = max(0.22, 1.0 - binding_pressure * (0.44 + hierarchy_bias))
        if protein.startswith("pea") and pea_sh_crosscheck:
            sulfur_binding_factor = min(1.0, sulfur_binding_factor + 0.06)
            sources.append(_runtime_source_label(pea_sh_crosscheck))
            process_state_calibration_ids.append(str(pea_sh_crosscheck.get("id", "unknown")))
        if process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"} and extrusion_disulfide_calibration:
            sulfur_binding_factor = max(0.18, sulfur_binding_factor - 0.08 * process_drive)
            sources.append(_runtime_source_label(extrusion_disulfide_calibration))
            process_state_calibration_ids.append(str(extrusion_disulfide_calibration.get("id", "unknown")))
        if sulfur_binding_prior:
            sources.append(_runtime_source_label(sulfur_binding_prior))
            runtime_prior_ids.append(str(sulfur_binding_prior.get("id", "unknown")))
        if protein_binding_prior:
            sources.append(_runtime_source_label(protein_binding_prior))
            runtime_prior_ids.append(str(protein_binding_prior.get("id", "unknown")))
        dynamic_retention_factor *= sulfur_binding_factor
        if sulfur_family:
            retention_mode = "sulfur_binding_prior" if retention_mode == "static_class_profile" else f"{retention_mode}+sulfur_binding_prior"
        elif pyrazine_family:
            retention_mode = "pyrazine_binding_prior" if retention_mode == "static_class_profile" else f"{retention_mode}+pyrazine_binding_prior"

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
        "partition_temperature_factor": float(partition_temperature_factor),
        "sulfur_binding_factor": float(sulfur_binding_factor),
        "runtime_prior_ids": sorted(set(runtime_prior_ids)),
        "process_state_calibration_ids": sorted(set(process_state_calibration_ids)),
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


def _extract_benchmark_marker_targets(
    benchmark_row: Mapping[str, Any],
    *,
    key_name: str,
) -> Dict[str, float]:
    key_values = benchmark_row.get("key_values", {}) if isinstance(benchmark_row.get("key_values", {}), Mapping) else {}
    raw_targets = key_values.get(key_name, {}) if isinstance(key_values, Mapping) else {}
    if not isinstance(raw_targets, Mapping):
        return {}
    return {str(compound): float(value) for compound, value in raw_targets.items()}


def _build_benchmark_anchor_context(benchmark_row: Mapping[str, Any]) -> Dict[str, Any]:
    key_values = benchmark_row.get("key_values", {}) if isinstance(benchmark_row.get("key_values", {}), Mapping) else {}
    context: Dict[str, Any] = {
        "id": str(benchmark_row.get("id", "unknown")),
        "citation": str(benchmark_row.get("citation", "unknown")),
        "matrix_family": str(benchmark_row.get("matrix_family", "unknown")),
        "status": str(benchmark_row.get("status", "unknown")),
        "volatile_output_mode": str(key_values.get("volatile_output_mode", "unknown")),
    }
    if "uht_temp_C" in key_values:
        context["temp_C"] = float(key_values.get("uht_temp_C", 0.0))
    elif "spi_temp_C" in key_values:
        context["temp_C"] = float(key_values.get("spi_temp_C", 0.0))
    if "uht_hold_seconds" in key_values:
        context["hold_seconds"] = float(key_values.get("uht_hold_seconds", 0.0))
    elif "heating_time_h" in key_values:
        context["heating_time_h"] = float(key_values.get("heating_time_h", 0.0))
    return context


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
    lipid_input_count = len(lipids)
    donor_pressure = 1.0 if any("ribose" in sugar or "xylose" in sugar or "fructose" in sugar for sugar in normalized_sugars) else 0.55
    split = summarize_lipid_runtime_split(
        {
            "Hexanal": lipid_marker_signals["hexanal"],
            "2-Pentylfuran": lipid_marker_signals["2_pentylfuran"],
            "1-Octen-3-ol": lipid_marker_signals["1_octen_3_ol"],
            "Nonanal": lipid_marker_signals["nonanal"],
            "2,4-Decadienal": lipid_marker_signals["2_4_decadienal"],
        },
        lipid_input_count=lipid_input_count,
        donor_pressure=donor_pressure,
        polyphenol_active=polyphenol_active,
    )
    lipid_marker_signal_ppb = float(split.get("lipid_marker_signal_ppb", 0.0))
    marker_count = int(split.get("lipid_marker_count", 0))
    carbonyl_competition_factor = float(split.get("carbonyl_competition_factor", 0.0))
    retention_rows = query_retention_reference_entries(
        family="02",
        protein_type=protein_type,
        process_state="heated_matrix",
    )
    benchmark_rows = query_benchmark_intake_entries(
        family="02",
        process_state="heated_matrix",
        primary_only=True,
    )
    benchmark_rows = [
        row
        for row in benchmark_rows
        if isinstance(row.get("key_values", {}), Mapping)
        and isinstance(row.get("key_values", {}).get("tracked_uht_markers_ug_per_l", {}), Mapping)
    ]
    primary_benchmark_row = next(
        (
            row for row in benchmark_rows
            if str(row.get("id", "")).strip() == "trikusuma_2019"
        ),
        benchmark_rows[0] if benchmark_rows else {},
    )
    primary_benchmark_ids = [str(primary_benchmark_row.get("id", "unknown"))] if primary_benchmark_row else []
    benchmark_marker_targets = _extract_benchmark_marker_targets(
        primary_benchmark_row,
        key_name="tracked_uht_markers_ug_per_l",
    )
    crosstalk_prior_rows = query_family_runtime_priors(
        runtime_lane="strecker_crosstalk",
        family="02",
        process_state="heated_matrix",
    )
    benchmark_ready_targets = _lipid_benchmark_target_panel()
    crosstalk_prior_active = bool(crosstalk_prior_rows) and bool(polyphenol_active or lipid_input_count or lipid_marker_signal_ppb > 0.0)
    strecker_suppression_factor = float(split.get("strecker_suppression_factor", 0.0))
    if polyphenol_active and crosstalk_prior_active and not bool(split.get("polyphenol_crosstalk_active", False)):
        strecker_suppression_factor = min(1.0, strecker_suppression_factor * 1.15)
    maillard_closure_pressure = float(split.get("maillard_closure_pressure", 0.0))
    active = bool(lipid_input_count or lipid_marker_signal_ppb > 0.0)
    summary = "No explicit lipid-driven crosstalk lane active."
    if active:
        summary = (
            "Lipid-derived adverse markers are present or lipid precursors are active, so the lane should be split into adverse-marker retention and carbonyl-competition support rather than treated as a generic off-note penalty."
        )
        if primary_benchmark_row:
            summary = (
                "Lipid-derived adverse markers are present or lipid precursors are active, so the lane is benchmark-facing against the Trikusuma 2019 heated pea panel and must stay split into adverse-marker retention and carbonyl-competition support."
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
            "dominant_marker": str(split.get("dominant_marker", "none")),
            "polyphenol_crosstalk_active": bool(split.get("polyphenol_crosstalk_active", False)),
            "carbonyl_competition_factor": float(carbonyl_competition_factor),
            "benchmark_ready_targets": [str(row.get("display_name", row.get("canonical_name", "unknown"))) for row in benchmark_ready_targets],
            "benchmark_anchor_ids": primary_benchmark_ids,
            "primary_benchmark_id": str(primary_benchmark_row.get("id", "unknown")) if primary_benchmark_row else "",
            "primary_benchmark_citation": str(primary_benchmark_row.get("citation", "")) if primary_benchmark_row else "",
            "benchmark_marker_targets_ug_per_l": benchmark_marker_targets,
            "benchmark_anchor_context": _build_benchmark_anchor_context(primary_benchmark_row) if primary_benchmark_row else {},
            "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
            "competition_prior_ids": [str(row.get("id", "unknown")) for row in crosstalk_prior_rows],
            "strecker_suppression_factor": float(strecker_suppression_factor),
            "maillard_closure_pressure": float(maillard_closure_pressure),
            "runtime_sub_lanes": {
                "adverse_marker_generation_and_retention": {
                    **dict(split.get("runtime_sub_lanes", {}).get("adverse_marker_generation_and_retention", {})),
                    "benchmark_ready_target_count": len(benchmark_ready_targets),
                    "benchmark_ready_targets": [str(row.get("display_name", row.get("canonical_name", "unknown"))) for row in benchmark_ready_targets],
                    "benchmark_anchor_count": len(primary_benchmark_ids),
                    "benchmark_anchor_ids": primary_benchmark_ids,
                    "benchmark_marker_targets_ug_per_l": benchmark_marker_targets,
                    "benchmark_anchor_context": _build_benchmark_anchor_context(primary_benchmark_row) if primary_benchmark_row else {},
                    "retention_reference_count": len(retention_rows),
                    "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
                },
                "carbonyl_competition_and_crosstalk": {
                    **dict(split.get("runtime_sub_lanes", {}).get("carbonyl_competition_and_crosstalk", {})),
                    "active": bool(crosstalk_prior_active),
                    "competition_prior_count": len(crosstalk_prior_rows),
                    "competition_prior_ids": [str(row.get("id", "unknown")) for row in crosstalk_prior_rows],
                    "strecker_suppression_factor": float(strecker_suppression_factor),
                    "maillard_closure_pressure": float(maillard_closure_pressure),
                },
            },
        },
    )


def _build_lipid_maillard_competition_lane(
    *,
    signal_map: Dict[str, float],
    normalized_sugars: List[str],
    normalized_amino: List[str],
    lipids: List[str],
    protein_type: Optional[str],
    lipid_crosstalk_lane: Mapping[str, Any],
) -> Dict[str, Any]:
    hexanal_signal = float(dict(lipid_crosstalk_lane.get("lipid_marker_signals_ppb", {})).get("hexanal", 0.0))
    total_lipid_signal = float(lipid_crosstalk_lane.get("lipid_marker_signal_ppb", 0.0))
    mft_signal = _compound_signal(signal_map, ["2-Methyl-3-furanthiol (MFT)", "2-methyl-3-furanthiol", "MFT"])
    fft_signal = _compound_signal(signal_map, ["2-Furfurylthiol (FFT)", "2-furfurylthiol", "FFT"])
    methional_signal = _compound_signal(signal_map, ["Methional (3-(methylthio)propanal)", "methional"])
    sulfur_partner_signal = max(mft_signal, fft_signal, methional_signal)
    cysteine_partner_active = any("cysteine" in value for value in normalized_amino)
    fast_donor_active = any(any(token in value for token in ["ribose", "xylose", "arabinose"]) for value in normalized_sugars)
    lipid_context_active = bool(lipids or lipid_crosstalk_lane.get("active") or total_lipid_signal > 0.0)
    competition_window_active = bool(lipid_context_active and (cysteine_partner_active or sulfur_partner_signal > 0.0))

    kinetic_prior = next(
        iter(
            query_dft_kinetic_priors(
                family="11",
                reaction_key="hexanal_radical_quench",
                process_state="heated_matrix",
            )
        ),
        {},
    )
    barrier_kj_mol = float(kinetic_prior.get("barrier_kj_mol", 31.72) or 31.72)
    quench_reactivity_factor = max(0.25, min(1.0, 1.05 - max(0.0, barrier_kj_mol - 25.0) / 45.0))

    benchmark_rows = query_benchmark_intake_entries(
        family="11",
        primary_only=True,
        kind="quantitative_benchmark",
    )
    hexanal_baseline_anchors_ug_per_kg: Dict[str, float] = {}
    suppression_anchor_fraction = 0.288
    for row in benchmark_rows:
        key_values = row.get("key_values", {}) if isinstance(row.get("key_values", {}), Mapping) else {}
        if "hexanal_control_ug_per_kg" in key_values:
            hexanal_baseline_anchors_ug_per_kg[str(row.get("id", "unknown"))] = float(key_values.get("hexanal_control_ug_per_kg", 0.0))
            suppression_anchor_fraction = float(key_values.get("hexanal_suppression_80min_pct", 28.8)) / 100.0
        elif "hexanal_raw_pea_ug_per_kg" in key_values:
            hexanal_baseline_anchors_ug_per_kg[str(row.get("id", "unknown"))] = float(key_values.get("hexanal_raw_pea_ug_per_kg", 0.0))

    protein_label = str(protein_type or "free")
    selected_anchor_id = ""
    if protein_label.startswith("pea") and "acs_2020_raw_pea_hexanal_baseline" in hexanal_baseline_anchors_ug_per_kg:
        selected_anchor_id = "acs_2020_raw_pea_hexanal_baseline"
    elif protein_label.startswith("soy") and "pmc_2026_hme_hexanal_baseline" in hexanal_baseline_anchors_ug_per_kg:
        selected_anchor_id = "pmc_2026_hme_hexanal_baseline"
    elif hexanal_baseline_anchors_ug_per_kg:
        selected_anchor_id = next(iter(hexanal_baseline_anchors_ug_per_kg))

    selected_hexanal_anchor = float(hexanal_baseline_anchors_ug_per_kg.get(selected_anchor_id, 0.0))
    proxy_scaled_anchor = max(25.0, selected_hexanal_anchor * 0.25) if selected_hexanal_anchor > 0.0 else 120.0
    hexanal_pressure = min(1.0, max(hexanal_signal, total_lipid_signal) / proxy_scaled_anchor) if competition_window_active else 0.0
    radical_partner_pressure = min(
        1.4,
        0.55 * float(cysteine_partner_active)
        + 0.25 * float(fast_donor_active)
        + 0.35 * min(1.5, sulfur_partner_signal / 12.0),
    )
    suppression_cap = min(0.45, suppression_anchor_fraction + (0.08 if fast_donor_active else 0.0))
    hexanal_suppression_fraction = min(
        0.45,
        suppression_cap * quench_reactivity_factor * radical_partner_pressure * max(0.35, hexanal_pressure),
    ) if competition_window_active else 0.0
    predicted_quenched_hexanal_ppb = max(0.0, hexanal_signal * (1.0 - hexanal_suppression_fraction))

    summary = "No explicit Maillard/lipid competition lane active."
    if competition_window_active:
        summary = (
            "Hexanal-bearing lipid chemistry now overlaps with sulfur-positive Maillard routing, so the Family 11 competition lane should expose explicit hexanal/MFT suppression priors rather than inheriting a generic lipid penalty."
        )

    kinetic_prior_ids = [str(kinetic_prior.get("id", kinetic_prior.get("reaction_key", "")))] if kinetic_prior else []
    return _build_family_lane_payload(
        "11",
        active=competition_window_active,
        summary=summary,
        metrics={
            "competition_window_active": bool(competition_window_active),
            "hexanal_signal_ppb": float(hexanal_signal),
            "sulfur_partner_signal_ppb": float(sulfur_partner_signal),
            "mft_signal_ppb": float(mft_signal),
            "fft_signal_ppb": float(fft_signal),
            "cysteine_partner_active": bool(cysteine_partner_active),
            "fast_donor_active": bool(fast_donor_active),
            "kinetic_prior_ids": kinetic_prior_ids,
            "kinetic_prior": {
                "reaction_key": str(kinetic_prior.get("reaction_key", "")),
                "active_arrhenius_key": str(kinetic_prior.get("active_arrhenius_key", "")),
                "barrier_kj_mol": float(barrier_kj_mol),
                "honest_label": str(kinetic_prior.get("honest_label", "")),
            },
            "upstream_crosstalk_prior_ids": [str(item) for item in lipid_crosstalk_lane.get("competition_prior_ids", [])],
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in benchmark_rows],
            "selected_benchmark_anchor_id": selected_anchor_id,
            "hexanal_baseline_anchors_ug_per_kg": hexanal_baseline_anchors_ug_per_kg,
            "suppression_anchor_fraction": float(suppression_anchor_fraction),
            "radical_quench_reactivity_factor": float(quench_reactivity_factor),
            "hexanal_pressure": float(hexanal_pressure),
            "radical_partner_pressure": float(radical_partner_pressure),
            "hexanal_suppression_fraction": float(hexanal_suppression_fraction),
            "predicted_quenched_hexanal_ppb": float(predicted_quenched_hexanal_ppb),
        },
    )


def _sum_effective_molar_ratios(
    effective_molar_ratios: Mapping[str, Any],
    *,
    tokens: List[str],
) -> float:
    """Sum the mM loading of every ingredient matching one of `tokens`.

    Word-boundary matched since 2026-08-27: substring matching parsed "asnase"
    (asparaginase, an acrylamide MITIGATION enzyme) as asparagine, i.e. as an
    acrylamide precursor at the enzyme's own concentration.
    Values are consumed as MOLAR_RATIO_CANONICAL_UNIT (mM).
    """
    total = 0.0
    requested = [str(token).strip().lower() for token in tokens if str(token).strip()]
    for raw_name, raw_value in _numeric_ratio_items(effective_molar_ratios).items():
        normalized_name = _normalize_name(str(raw_name))
        if _any_token_matches(normalized_name, requested):
            total += raw_value
    return float(total)


def _build_protein_damage_markers_lane(
    *,
    normalized_sugars: List[str],
    normalized_amino: List[str],
    normalized_additives: List[str],
    protein_type: Optional[str],
    process_state: Optional[str],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
    water_activity: Optional[float],
    effective_molar_ratios: Mapping[str, Any],
    pH: Optional[float] = None,
    concentration_unit: Optional[str] = None,
) -> Dict[str, Any]:
    temperature = float(temperature_celsius if temperature_celsius is not None else 25.0)
    duration = float(time_minutes if time_minutes is not None else 0.0)
    protein_label = str(protein_type or "free")
    process_label = str(process_state or "heated_matrix")

    unit_assessment = _assess_concentration_unit(
        dict(effective_molar_ratios or {}, **({"concentration_unit": concentration_unit} if concentration_unit else {})),
        context="protein_damage_markers lane (family 12)",
    )

    reducing_sugar_active = any(
        _any_token_matches(sugar, ["ribose", "glucose", "fructose", "xylose", "maltose", "sugar"])
        for sugar in normalized_sugars
    )
    polyphenol_factor = 0.20 if any("polyphenol" in value or "grape seed" in value or "catechin" in value for value in normalized_additives) else 0.0
    soy_isoflavone_prior = _computational_prior_entry("nakagawa_2004_isoflavone_dicarbonyl_sink_v1") if protein_label.startswith("soy") else {}
    soy_isoflavone_factor = 0.0
    if soy_isoflavone_prior:
        soy_isoflavone_factor = float(soy_isoflavone_prior.get("age_inhibition_fraction_at_100uM", 0.0) or 0.0)
        soy_isoflavone_factor = max(0.0, min(1.0, 0.45 * soy_isoflavone_factor))
    explicit_sugar_mM = _sum_effective_molar_ratios(
        effective_molar_ratios,
        tokens=["ribose", "glucose", "fructose", "xylose", "maltose", "sugar"],
    )
    explicit_lysine_mM = _sum_effective_molar_ratios(
        effective_molar_ratios,
        tokens=["lysine", "lys"],
    )
    explicit_asparagine_mM = _sum_effective_molar_ratios(
        effective_molar_ratios,
        tokens=["asparagine", "asn"],
    )

    protein_bound_lysine_proxy = {
        "free": 0.0,
        "pea_iso": 0.45,
        "pea_conc": 0.38,
        "soy_iso": 0.62,
        "soy_conc": 0.55,
        "myco": 0.30,
    }.get(protein_label, 0.35)
    reducing_sugar_mM = explicit_sugar_mM if explicit_sugar_mM > 0.0 else (0.90 if reducing_sugar_active else 0.0)
    lysine_mM = explicit_lysine_mM if explicit_lysine_mM > 0.0 else protein_bound_lysine_proxy
    asparagine_mM = explicit_asparagine_mM

    benchmark_rows = query_benchmark_intake_entries(
        family="12",
        primary_only=True,
        kind="quantitative_benchmark",
    )
    benchmark_by_id = {str(row.get("id", "unknown")): dict(row) for row in benchmark_rows}
    benchmark_anchor_ids = [str(row.get("id", "unknown")) for row in benchmark_rows]
    benchmark_anchor_citations = [str(row.get("citation", "unknown")) for row in benchmark_rows]

    foods_benchmark = benchmark_by_id.get("foods_2023_cml_cel_proxy_benchmark", {})
    lysine_loss_benchmark = benchmark_by_id.get("acs_2022_pba_lysine_loss_benchmark", {})
    furosine_benchmark = benchmark_by_id.get("ramirez_jimenez_2000_furosine_crossover_benchmark", {})
    acrylamide_benchmark = benchmark_by_id.get("acs_ref3_spi_acrylamide_fast_kinetics", {})

    # UNIT ALIGNMENT (2026-08-27). src.safety's predictors return ppb (ug/kg of
    # food) - see the UNITS contract there. The intake-registry anchors are
    # authored in mg/kg and mg per 100 g protein, and were previously divided
    # into ppb predictions directly, which is a 1000x (CML/CEL) and 2000x
    # (furosine) unit error. Every anchor is converted to ppb here, once.
    # The furosine ratio in particular used to be min()-capped at 2.0 for every
    # heated formulation; it is a real, varying quantity again.
    representative_cml_ppb = MG_PER_KG_TO_PPB * float(
        (foods_benchmark.get("key_values", {}) or {}).get("representative_cml_mg_per_kg", 32.0) or 32.0
    )
    representative_cel_ppb = MG_PER_KG_TO_PPB * float(
        (foods_benchmark.get("key_values", {}) or {}).get("representative_cel_mg_per_kg", 55.0) or 55.0
    )
    furosine_peak_ppb = mg_per_100g_protein_to_ppb(
        float((furosine_benchmark.get("key_values", {}) or {}).get("peak_furosine_mg_per_100g_protein", 8.7) or 8.7)
    )
    # ug/kg IS ppb - no conversion, stated so the asymmetry is not read as an oversight.
    acrylamide_anchor_ppb = float(
        (acrylamide_benchmark.get("key_values", {}) or {}).get("acrylamide_end_ug_per_kg", 62.62) or 62.62
    )
    lysine_loss_anchor = float(
        (lysine_loss_benchmark.get("key_values", {}) or {}).get(
            "spi_lysine_loss_pct" if protein_label.startswith("soy") else "ppi_lysine_loss_pct",
            31.4 if protein_label.startswith("soy") else 15.4,
        )
        or (31.4 if protein_label.startswith("soy") else 15.4)
    )

    thermal_damage_window = _sigmoid(temperature, 118.0, 10.0)
    extrusion_damage_window = _sigmoid(temperature, 138.0, 7.5) if process_label in {"extrusion_structured", "hme", "lme"} else 0.0
    residence_scale = max(0.40, min(1.35, 0.80 + 0.12 * math.log1p(max(duration, 0.0))))
    lysine_loss_pct = lysine_loss_anchor * max(0.15, thermal_damage_window) * residence_scale

    # pH plumbing (2026-08-27): this used to read
    # `pH=float(6.0 if process_state is None else 6.0)` - a tautology that pinned
    # every runtime acrylamide prediction to pH 6.0 regardless of the
    # formulation. The formulation pH is now used, falling back to 6.0 only when
    # the caller genuinely has none.
    acrylamide_pH = float(pH) if pH is not None else 6.0
    predicted_acrylamide_ppb = 0.0
    acrylamide_assessment_status = "not_assessed"
    acrylamide_not_assessed_reason = ""
    if asparagine_mM <= 0.0:
        acrylamide_not_assessed_reason = "no_asparagine_in_formulation"
    elif reducing_sugar_mM <= 0.0:
        acrylamide_not_assessed_reason = "no_reducing_sugar_in_formulation"
    elif temperature < 105.0:
        acrylamide_not_assessed_reason = "below_acrylamide_lane_temperature_floor_105C"
    if asparagine_mM > 0.0 and reducing_sugar_mM > 0.0 and temperature >= 105.0:
        acrylamide_assessment_status = "assessed"
        predicted_acrylamide_ppb = float(
            predict_acrylamide(
                asparagine_mM=asparagine_mM,
                reducing_sugar_mM=reducing_sugar_mM,
                temp_C=temperature,
                time_min=max(duration, 1.0),
                pH=acrylamide_pH,
                water_activity=water_activity,
                effective_temp_c=temperature if process_label in {"extrusion_structured", "hme", "lme"} else None,
            ).acrylamide_ppb
        )

    predicted_cml_proxy = 0.0
    predicted_cel_proxy = 0.0
    predicted_furosine_proxy = 0.0
    if lysine_mM > 0.0 and reducing_sugar_mM > 0.0 and temperature >= 95.0:
        predicted_cml_proxy = float(
            predict_cml(
                lysine_mM=lysine_mM,
                reducing_sugar_mM=reducing_sugar_mM,
                temp_C=temperature,
                time_min=max(duration, 1.0),
                water_activity=water_activity,
                effective_temp_c=temperature if process_label in {"extrusion_structured", "hme", "lme"} else None,
                polyphenol_factor=polyphenol_factor,
                soy_isoflavone_factor=soy_isoflavone_factor,
            )
        )
        predicted_cel_proxy = float(
            predict_cel(
                lysine_mM=lysine_mM,
                reducing_sugar_mM=reducing_sugar_mM,
                temp_C=temperature,
                time_min=max(duration, 1.0),
                water_activity=water_activity,
                effective_temp_c=temperature if process_label in {"extrusion_structured", "hme", "lme"} else None,
                polyphenol_factor=polyphenol_factor,
                soy_isoflavone_factor=soy_isoflavone_factor,
            )
        )
        predicted_furosine_proxy = float(
            predict_furosine(
                temp_C=temperature,
                time_min=max(duration, 1.0),
                lysine_mM=lysine_mM,
                reducing_sugar_mM=reducing_sugar_mM,
                protein_type=protein_label,
                water_activity=water_activity,
                effective_temp_c=temperature if process_label in {"extrusion_structured", "hme", "lme"} else None,
            )
        )

    acrylamide_ratio = predicted_acrylamide_ppb / max(acrylamide_anchor_ppb, 1.0e-6) if predicted_acrylamide_ppb > 0.0 else 0.0
    cml_ratio = predicted_cml_proxy / max(representative_cml_ppb, 1.0e-6) if predicted_cml_proxy > 0.0 else 0.0
    cel_ratio = predicted_cel_proxy / max(representative_cel_ppb, 1.0e-6) if predicted_cel_proxy > 0.0 else 0.0
    furosine_ratio = predicted_furosine_proxy / max(furosine_peak_ppb, 1.0e-6) if predicted_furosine_proxy > 0.0 else 0.0
    lysine_loss_ratio = lysine_loss_pct / max(lysine_loss_anchor, 1.0e-6) if lysine_loss_pct > 0.0 else 0.0
    damage_burden_score = min(
        1.6,
        0.22 * min(2.0, acrylamide_ratio)
        + 0.18 * min(2.0, cml_ratio)
        + 0.24 * min(2.0, cel_ratio)
        + 0.16 * min(2.0, furosine_ratio)
        + 0.20 * min(2.0, lysine_loss_ratio),
    )

    selected_benchmark_anchor_id = ""
    anchor_pressure = {
        "acs_ref3_spi_acrylamide_fast_kinetics": acrylamide_ratio,
        "foods_2023_cml_cel_proxy_benchmark": max(cml_ratio, cel_ratio),
        "ramirez_jimenez_2000_furosine_crossover_benchmark": furosine_ratio,
        "acs_2022_pba_lysine_loss_benchmark": lysine_loss_ratio,
    }
    if any(value > 0.0 for value in anchor_pressure.values()):
        selected_benchmark_anchor_id = max(anchor_pressure.items(), key=lambda item: item[1])[0]

    foods_reference = get_safety_reference_payload("foods_2023_pba_cml_cel_benchmark") or {}
    foods_cml_range = get_safety_reference_range("commercial_pbma_burgers", "foods_2023_pba_cml_cel_benchmark") or {}
    pmc_cml_range = get_safety_reference_range("commercial_pba_products", "pmc_2024_pba_cml_cel_ranges") or {}
    furosine_reference = get_safety_reference_payload("ramirez_jimenez_2000_furosine_crossover") or {}

    active = bool(
        (temperature_celsius is not None and temperature >= 95.0 and reducing_sugar_mM > 0.0 and (lysine_mM > 0.0 or asparagine_mM > 0.0 or protein_label != "free"))
        or predicted_acrylamide_ppb > 0.0
        or predicted_cml_proxy > 0.0
        or predicted_cel_proxy > 0.0
        or predicted_furosine_proxy > 0.0
    )
    summary = "No explicit protein-damage guardrail lane active."
    if active:
        summary = (
            "Protein damage markers now travel as a first-class guardrail lane, so reactive lysine loss, AGE proxies, and furosine crossover remain visible alongside aroma-positive lanes instead of staying buried inside aggregate safety output."
        )

    return _build_family_lane_payload(
        "12",
        active=active,
        summary=summary,
        metrics={
            "damage_guardrail_active": bool(active),
            "reducing_sugar_mM": float(reducing_sugar_mM),
            "reactive_lysine_mM": float(lysine_mM),
            "asparagine_mM": float(asparagine_mM),
            "reactive_lysine_loss_pct": float(lysine_loss_pct),
            "predicted_acrylamide_ppb": float(predicted_acrylamide_ppb),
            "acrylamide_assessment_status": acrylamide_assessment_status,
            "acrylamide_not_assessed_reason": acrylamide_not_assessed_reason,
            "acrylamide_pH_used": float(acrylamide_pH),
            "predicted_cml_proxy": float(predicted_cml_proxy),
            "predicted_cel_proxy": float(predicted_cel_proxy),
            "predicted_furosine_proxy": float(predicted_furosine_proxy),
            # Units are declared, not implied: every predicted_* value above and
            # every *_anchor below is ppb (ug/kg of food).
            "damage_marker_units": "ppb (ug/kg food)",
            "concentration_input_unit": MOLAR_RATIO_CANONICAL_UNIT,
            "declared_concentration_unit": unit_assessment["declared_concentration_unit"],
            "concentration_unit_warning": bool(unit_assessment["looks_like_unitless_ratios"]),
            "acrylamide_anchor_ppb": float(acrylamide_anchor_ppb),
            "cml_anchor_ppb": float(representative_cml_ppb),
            "cel_anchor_ppb": float(representative_cel_ppb),
            "furosine_anchor_ppb": float(furosine_peak_ppb),
            "furosine_anchor_conversion": (
                "8.7 mg/100 g protein x 20% protein x 10 -> 17.4 mg/kg x 1000 -> 17400 ppb"
            ),
            "benchmark_anchor_ids": benchmark_anchor_ids,
            "benchmark_anchor_citations": benchmark_anchor_citations,
            "selected_benchmark_anchor_id": selected_benchmark_anchor_id,
            "safety_reference_ids": [
                "foods_2023_pba_cml_cel_benchmark",
                "pmc_2024_pba_cml_cel_ranges",
                "ramirez_jimenez_2000_furosine_crossover",
            ],
            "runtime_prior_ids": [str(soy_isoflavone_prior.get("id", ""))] if soy_isoflavone_prior else [],
            "safety_reference_citations": [
                str(foods_reference.get("source_citation", "")),
                str((get_safety_reference_payload("pmc_2024_pba_cml_cel_ranges") or {}).get("source_citation", "")),
                str(furosine_reference.get("source_citation", "")),
            ],
            "runtime_prior_citations": [str(soy_isoflavone_prior.get("source", ""))] if soy_isoflavone_prior else [],
            "age_reference_ranges_mg_per_kg": {
                "foods_2023_cml": dict(foods_cml_range),
                "pmc_2024_cml_cel": dict(pmc_cml_range),
            },
            "furosine_reference_range": dict(get_safety_reference_range("mild_legume_extrudate", "ramirez_jimenez_2000_furosine_crossover") or {}),
            "acrylamide_ratio_vs_anchor": float(acrylamide_ratio),
            "cml_ratio_vs_foods_anchor": float(cml_ratio),
            "cel_ratio_vs_foods_anchor": float(cel_ratio),
            "furosine_ratio_vs_crossover": float(furosine_ratio),
            "lysine_loss_ratio_vs_anchor": float(lysine_loss_ratio),
            "damage_burden_score": float(damage_burden_score),
            "thermal_damage_window": float(thermal_damage_window),
            "extrusion_damage_window": float(extrusion_damage_window),
            "soy_isoflavone_factor": float(soy_isoflavone_factor),
        },
    )


def _build_polyphenol_amino_capping_lane(
    *,
    normalized_additives: List[str],
    normalized_amino: List[str],
    normalized_support_cues: List[str],
    normalized_interventions: List[str],
    protein_label: str,
    pH: Optional[float],
    temperature_celsius: Optional[float],
    process_state: Optional[str],
) -> Dict[str, Any]:
    strecker_crosstalk_priors = get_strecker_crosstalk_priors()
    polyphenol_tokens = [
        str(token).strip().lower()
        for token in strecker_crosstalk_priors.get("polyphenol_examples", [])
        if str(token).strip()
    ]
    polyphenol_tokens.extend([
        "polyphenol",
        "chlorogenic acid",
        "cga",
        "caffeic acid",
        "egcg",
        "epigallocatechin gallate",
    ])
    support_pool = normalized_additives + normalized_support_cues + normalized_interventions
    polyphenol_active = any(
        token in value
        for value in support_pool
        for token in polyphenol_tokens
    )
    cysteine_target_present = any("cysteine" in value for value in normalized_amino + normalized_additives)
    lysine_target_present = any("lysine" in value for value in normalized_amino + normalized_additives) or protein_label != "free"

    kinetic_prior = next(
        iter(
            query_dft_kinetic_priors(
                family="13",
                reaction_key="quinone_cys_michael",
                process_state=process_state,
            )
        ),
        {},
    )
    barrier_kj_mol = float(kinetic_prior.get("barrier_kj_mol", 74.65) or 74.65)
    pH_window = _sigmoid(6.2 if pH is None else float(pH), 6.4, 0.75)
    reference_temperature_celsius = 105.0 if temperature_celsius is None else float(temperature_celsius)
    thermal_drive = _sigmoid(reference_temperature_celsius, 85.0, 14.0)
    kinetic_feasibility = max(0.20, min(1.0, 1.10 - max(0.0, barrier_kj_mol - 55.0) / 65.0))
    quinone_budget = float(polyphenol_active) * thermal_drive * (0.45 + 0.55 * pH_window) * kinetic_feasibility
    thiol_preference_ratio_vs_amine = 1.0e5
    cysteine_depletion_factor = 0.0
    lysine_depletion_factor = 0.0
    if polyphenol_active:
        cysteine_depletion_factor = min(0.80, quinone_budget * (0.58 if cysteine_target_present else 0.18))
        lysine_channel_weight = 0.0
        if lysine_target_present:
            lysine_channel_weight = 0.18 if cysteine_target_present else 0.55
        lysine_depletion_factor = min(0.28, quinone_budget * lysine_channel_weight)

    active = bool(polyphenol_active and quinone_budget > 0.05)
    summary = "No explicit polyphenol-amino capping lane active."
    if active:
        summary = (
            "Polyphenol chemistry is active and now acts as an upstream precursor sink, so quinone capture can deplete sulfur-positive Maillard partners before they reach the benchmark-facing aroma lanes."
        )
        if kinetic_prior:
            summary = (
                "Polyphenol chemistry is active and now acts as an upstream precursor sink, bounded by the Family 13 quinone-thiol kinetic prior rather than by a generic antioxidant keyword."
            )

    return _build_family_lane_payload(
        "13",
        active=active,
        summary=summary,
        metrics={
            "polyphenol_active": bool(polyphenol_active),
            "cysteine_target_present": bool(cysteine_target_present),
            "lysine_target_present": bool(lysine_target_present),
            "quinone_budget": float(quinone_budget),
            "quinone_generation_factor": float((0.45 + 0.55 * pH_window) * thermal_drive) if polyphenol_active else 0.0,
            "cysteine_depletion_factor": float(cysteine_depletion_factor),
            "lysine_depletion_factor": float(lysine_depletion_factor),
            "thiol_preference_ratio_vs_amine": float(thiol_preference_ratio_vs_amine),
            "kinetic_prior_ids": [str(kinetic_prior.get("id", kinetic_prior.get("reaction_key", "")))] if kinetic_prior else [],
            "kinetic_prior": {
                "reaction_key": str(kinetic_prior.get("reaction_key", "")),
                "active_arrhenius_key": str(kinetic_prior.get("active_arrhenius_key", "")),
                "barrier_kj_mol": float(barrier_kj_mol),
                "honest_label": str(kinetic_prior.get("honest_label", "")),
            },
        },
    )


def _build_ascorbic_acid_maillard_lane(
    *,
    normalized_additives: List[str],
    normalized_support_cues: List[str],
    normalized_interventions: List[str],
    normalized_amino: List[str],
    protein_label: str,
    pH: Optional[float],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
    water_activity: Optional[float],
    process_state: Optional[str],
) -> Dict[str, Any]:
    support_pool = normalized_additives + normalized_support_cues + normalized_interventions
    ascorbic_tokens = [
        "ascorbic acid",
        "ascorbate",
        "vitamin c",
        "sodium ascorbate",
        "calcium ascorbate",
    ]
    ascorbic_acid_active = any(
        token in value
        for value in support_pool
        for token in ascorbic_tokens
    )

    kinetic_prior = next(
        iter(
            query_dft_kinetic_priors(
                family="14",
                reaction_key="aa_ring_open_dicarbonyl",
                process_state=process_state,
            )
        ),
        {},
    )
    if pH is None:
        effective_ea_kj_mol = 31.70
    else:
        effective_ea_kj_mol = _interpolate_profile(
            [5.0, 7.0, 9.5],
            [15.77, 31.70, 47.53],
            float(pH),
        ) or 31.70
    aw = None if water_activity is None else max(0.05, min(0.98, float(water_activity)))
    water_activity_modulation_factor = 0.70
    if aw is not None:
        aw_peak_distance = (aw - 0.64) / 0.24
        water_activity_modulation_factor = max(0.18, 1.0 - aw_peak_distance * aw_peak_distance)
    thermal_drive = _sigmoid(25.0 if temperature_celsius is None else float(temperature_celsius), 122.0, 13.0)
    residence_factor = min(1.20, 0.45 + 0.18 * math.log1p(max(0.0, float(time_minutes or 0.0))))
    pH_reactivity_factor = max(0.35, min(1.55, 31.70 / max(float(effective_ea_kj_mol), 1.0e-6)))
    dicarbonyl_source_pressure = float(ascorbic_acid_active) * thermal_drive * residence_factor * water_activity_modulation_factor * pH_reactivity_factor
    ascorbic_acid_loss = min(1.0, dicarbonyl_source_pressure * 0.82)
    crosslink_partner_present = any(
        any(token in value for token in ["lysine", "arginine", "histidine"])
        for value in normalized_amino + normalized_additives
    ) or protein_label != "free"
    pentosidine_load = min(
        1.0,
        dicarbonyl_source_pressure * (0.72 if crosslink_partner_present else 0.25),
    )

    active = bool(ascorbic_acid_active and dicarbonyl_source_pressure > 0.05)
    summary = "No explicit ascorbic-acid Maillard lane active."
    if active:
        summary = (
            "Ascorbic acid is active as a bounded upstream dicarbonyl source, so thermal severity can raise cross-link pressure and safety burden even when the main sugar route looks unchanged."
        )
        if kinetic_prior:
            summary = (
                "Ascorbic acid is active as a bounded upstream dicarbonyl source, anchored to the Family 14 HCW kinetic prior instead of being silently folded into the core sugar lane."
            )

    return _build_family_lane_payload(
        "14",
        active=active,
        summary=summary,
        metrics={
            "ascorbic_acid_active": bool(ascorbic_acid_active),
            "crosslink_partner_present": bool(crosslink_partner_present),
            "effective_ea_kj_mol": float(effective_ea_kj_mol),
            "water_activity_modulation_factor": float(water_activity_modulation_factor),
            "thermal_drive": float(thermal_drive),
            "residence_factor": float(residence_factor),
            "pH_reactivity_factor": float(pH_reactivity_factor),
            "ascorbic_acid_loss": float(ascorbic_acid_loss),
            "dicarbonyl_source_pressure": float(dicarbonyl_source_pressure),
            "pentosidine_load": float(pentosidine_load),
            "kinetic_prior_ids": [str(kinetic_prior.get("id", kinetic_prior.get("reaction_key", "")))] if kinetic_prior else [],
            "kinetic_prior": {
                "reaction_key": str(kinetic_prior.get("reaction_key", "")),
                "active_arrhenius_key": str(kinetic_prior.get("active_arrhenius_key", "")),
                "barrier_kj_mol": float(kinetic_prior.get("barrier_kj_mol", 31.70) or 31.70),
                "honest_label": str(kinetic_prior.get("honest_label", "")),
            },
        },
    )


def _build_phospholipid_amine_sink_lane(
    *,
    normalized_sugars: List[str],
    normalized_additives: List[str],
    normalized_lipids: List[str],
    normalized_support_cues: List[str],
    pH: Optional[float],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
    process_state: Optional[str],
) -> Dict[str, Any]:
    support_pool = normalized_additives + normalized_lipids + normalized_support_cues
    explicit_pe_active = any(
        any(token in value for token in ["phosphatidylethanolamine", "pe ", " pe", "soy lecithin", "lecithin", "phospholipid", "phosphatidyl"])
        for value in support_pool
    )
    sugar_present = bool(normalized_sugars)

    schiff_prior = next(
        iter(
            query_dft_kinetic_priors(
                family="15",
                reaction_key="pe_schiff_base",
                process_state=process_state,
            )
        ),
        {},
    )
    amadori_prior = next(
        iter(
            query_dft_kinetic_priors(
                family="15",
                reaction_key="pe_amadori",
                process_state=process_state,
            )
        ),
        {},
    )

    reference_temperature_celsius = 25.0 if temperature_celsius is None else float(temperature_celsius)
    thermal_drive = _sigmoid(reference_temperature_celsius, 118.0, 16.0)
    residence_factor = min(1.15, 0.42 + 0.18 * math.log1p(max(0.0, float(time_minutes or 0.0))))
    pH_window = _sigmoid(6.2 if pH is None else float(pH), 6.0, 0.85)
    schiff_feasibility = max(0.18, min(1.0, 1.08 - max(0.0, float(schiff_prior.get("barrier_kj_mol", 92.9) or 92.9) - 75.0) / 80.0))
    amadori_feasibility = max(0.24, min(1.0, 1.10 - max(0.0, float(amadori_prior.get("barrier_kj_mol", 82.9) or 82.9) - 68.0) / 75.0))
    pe_cue_strength = 1.0 if any(any(token in value for token in ["phosphatidylethanolamine", "pe ", " pe"]) for value in support_pool) else 0.72
    glycation_pressure = float(explicit_pe_active) * pe_cue_strength * thermal_drive * residence_factor * (0.35 + 0.65 * pH_window)
    pe_glycation_fraction = min(0.45, glycation_pressure * (0.22 + 0.18 * schiff_feasibility + 0.22 * amadori_feasibility) * float(sugar_present))
    sugar_sink_fraction = min(0.30, glycation_pressure * (0.18 + 0.24 * amadori_feasibility) * float(sugar_present))
    available_sugar_retention_factor = max(0.70, 1.0 - sugar_sink_fraction)

    active = bool(explicit_pe_active and sugar_present and sugar_sink_fraction > 0.02)
    summary = "No explicit phospholipid-amine sugar-sink lane active."
    if active:
        summary = (
            "Phospholipid-amine chemistry is active as a bounded upstream sugar sink, so part of the reducing-sugar pool should be treated as diverted into PE glycation before it reaches the core amino acid lane."
        )
        if schiff_prior and amadori_prior:
            summary = (
                "Phospholipid-amine chemistry is active as a bounded upstream sugar sink, anchored to the Family 15 PE Schiff-base and PE Amadori priors instead of being silently absorbed into the generic lipid background."
            )

    return _build_family_lane_payload(
        "15",
        active=active,
        summary=summary,
        metrics={
            "phospholipid_active": bool(explicit_pe_active),
            "pe_glycation_fraction": float(pe_glycation_fraction),
            "sugar_sink_fraction": float(sugar_sink_fraction),
            "available_sugar_retention_factor": float(available_sugar_retention_factor),
            "thermal_drive": float(thermal_drive),
            "residence_factor": float(residence_factor),
            "pH_window": float(pH_window),
            "kinetic_prior_ids": [
                str(schiff_prior.get("id", schiff_prior.get("reaction_key", ""))),
                str(amadori_prior.get("id", amadori_prior.get("reaction_key", ""))),
            ] if schiff_prior or amadori_prior else [],
            "kinetic_priors": {
                "pe_schiff_base": {
                    "reaction_key": str(schiff_prior.get("reaction_key", "")),
                    "active_arrhenius_key": str(schiff_prior.get("active_arrhenius_key", "")),
                    "barrier_kj_mol": float(schiff_prior.get("barrier_kj_mol", 92.9) or 92.9),
                },
                "pe_amadori": {
                    "reaction_key": str(amadori_prior.get("reaction_key", "")),
                    "active_arrhenius_key": str(amadori_prior.get("active_arrhenius_key", "")),
                    "barrier_kj_mol": float(amadori_prior.get("barrier_kj_mol", 82.9) or 82.9),
                },
            },
        },
    )


def _build_melanoidin_polymerization_lane(
    *,
    normalized_sugars: List[str],
    normalized_additives: List[str],
    normalized_amino: List[str],
    protein_label: str,
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
    water_activity: Optional[float],
    process_state: Optional[str],
) -> Dict[str, Any]:
    browning_core_active = bool(normalized_sugars and (normalized_amino or protein_label != "free"))
    support_pool = normalized_additives + normalized_amino
    sulfur_partner_active = any(any(token in value for token in ["cysteine", "methionine", "glutathione"]) for value in support_pool)
    gum_arabic_active = any(any(token in value for token in ["gum arabic", "arabic gum"]) for value in normalized_additives)
    hydrolysate_support_active = any(any(token in value for token in ["hydrolysate", "peptide", "autolysate"]) for value in support_pool)
    thermal_drive = _sigmoid(25.0 if temperature_celsius is None else float(temperature_celsius), 118.0, 15.0)
    residence_factor = 0.0 if time_minutes is None else min(1.25, 0.40 + 0.19 * math.log1p(max(0.0, float(time_minutes or 0.0))))
    aw = None if water_activity is None else max(0.05, min(0.98, float(water_activity)))
    browning_aw_factor = 0.82
    if aw is not None:
        aw_distance = (aw - 0.58) / 0.28
        browning_aw_factor = max(0.22, 1.0 - aw_distance * aw_distance)
    process_bias = 1.10 if process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"} else 0.92 if process_state == "heated_matrix" else 0.55

    retention_rows = query_retention_reference_entries(
        family="16",
        compound_name="2-Furfurylthiol",
        process_state=process_state,
    )
    architecture_row = next(
        (
            row
            for row in query_benchmark_intake_entries(family="16", primary_only=True)
            if str(row.get("id", "")).strip() == "jafc_2019_ref21_pea_gum_arabic_architecture_anchor"
        ),
        {},
    )
    architecture_calibration = _process_state_calibration_entry("jafc_2019_ref21_pea_gum_arabic_architecture_state")
    architecture_prior_rows = query_family_runtime_priors(
        family="16",
        entry_id="jafc_2019_ref21_pea_gum_arabic_architecture_v1",
    )
    architecture_prior = architecture_prior_rows[0] if architecture_prior_rows else {}
    fft_fold_reduction_anchor = max(
        [
            float(row.get("numeric_reference", {}).get("value", 0.0))
            for row in retention_rows
            if isinstance(row.get("numeric_reference", {}), Mapping)
            and str(row.get("numeric_reference", {}).get("units", "")).strip() == "fold_reduction"
        ]
        or [0.0]
    )
    melanoidin_mass = min(1.40, float(browning_core_active) * thermal_drive * residence_factor * browning_aw_factor * process_bias)
    thiol_scavenging_factor = min(
        0.92,
        melanoidin_mass * (0.24 + 0.34 * min(1.0, fft_fold_reduction_anchor / 16.0)) * (1.15 if sulfur_partner_active else 0.85),
    )
    matrix_trapping_pressure = min(1.0, melanoidin_mass * (0.55 + 0.45 * float(sulfur_partner_active)))
    architecture_anchors = architecture_calibration.get("numeric_anchors", {}) if isinstance(architecture_calibration.get("numeric_anchors", {}), Mapping) else {}
    native_mw_kda = float(architecture_anchors.get("native_mw_kda", 0.0) or 0.0)
    conjugated_mw_range = architecture_anchors.get("conjugated_mw_kda_range", [])
    native_rg_nm = float(architecture_anchors.get("native_radius_of_gyration_nm", 0.0) or 0.0)
    conjugated_rg_range = architecture_anchors.get("conjugated_radius_of_gyration_nm_range", [])
    conjugated_mw_midpoint = 0.0
    if isinstance(conjugated_mw_range, list) and len(conjugated_mw_range) == 2:
        conjugated_mw_midpoint = 0.5 * (float(conjugated_mw_range[0]) + float(conjugated_mw_range[1]))
    conjugated_rg_midpoint = 0.0
    if isinstance(conjugated_rg_range, list) and len(conjugated_rg_range) == 2:
        conjugated_rg_midpoint = 0.5 * (float(conjugated_rg_range[0]) + float(conjugated_rg_range[1]))
    native_mark_houwink = float(architecture_anchors.get("native_mark_houwink_exponent", 0.0) or 0.0)
    conjugated_mark_houwink = float(architecture_anchors.get("conjugated_mark_houwink_exponent", 0.0) or 0.0)
    mark_houwink_drop_fraction = 0.0
    if native_mark_houwink > 0.0 and conjugated_mark_houwink > 0.0:
        mark_houwink_drop_fraction = max(0.0, (native_mark_houwink - conjugated_mark_houwink) / native_mark_houwink)
    mw_growth_fraction = 0.0
    if native_mw_kda > 0.0 and conjugated_mw_midpoint > 0.0:
        mw_growth_fraction = max(0.0, (conjugated_mw_midpoint - native_mw_kda) / native_mw_kda)
    rg_growth_fraction = 0.0
    if native_rg_nm > 0.0 and conjugated_rg_midpoint > 0.0:
        rg_growth_fraction = max(0.0, (conjugated_rg_midpoint - native_rg_nm) / native_rg_nm)
    architecture_shift_active = bool(architecture_calibration and architecture_row and gum_arabic_active and hydrolysate_support_active)
    architecture_shift_score = min(
        1.0,
        1.8 * mw_growth_fraction + 2.4 * rg_growth_fraction + 1.5 * mark_houwink_drop_fraction,
    ) if architecture_shift_active else 0.0
    if architecture_shift_active:
        thiol_scavenging_factor = min(0.92, thiol_scavenging_factor + 0.10 * architecture_shift_score)
        matrix_trapping_pressure = min(1.0, matrix_trapping_pressure + 0.12 * architecture_shift_score)

    active = bool(time_minutes is not None and browning_core_active and melanoidin_mass > 0.10)
    summary = "No explicit melanoidin-polymerization trapping lane active."
    if active:
        summary = (
            "Browning burden is high enough that melanoidin build-up should be treated as a bounded trapping lane, reducing sulfur-positive observability rather than being interpreted as uniformly helpful browning."
        )
        if retention_rows:
            summary = (
                "Browning burden is high enough that melanoidin build-up should be treated as a bounded trapping lane, anchored to the Family 16 FFT trapping reference instead of a generic dark-color heuristic."
            )
        if architecture_shift_active:
            summary = (
                "Browning burden is high enough that melanoidin build-up should be treated as a bounded trapping lane, and gum arabic plus hydrolysate cues now also activate the Ref. 21 architecture anchor for higher polymer-growth pressure."
            )

    return _build_family_lane_payload(
        "16",
        active=active,
        summary=summary,
        metrics={
            "melanoidin_mass": float(melanoidin_mass),
            "thiol_scavenging_factor": float(thiol_scavenging_factor),
            "matrix_trapping_pressure": float(matrix_trapping_pressure),
            "thermal_drive": float(thermal_drive),
            "residence_factor": float(residence_factor),
            "browning_aw_factor": float(browning_aw_factor),
            "sulfur_partner_active": bool(sulfur_partner_active),
            "gum_arabic_active": bool(gum_arabic_active),
            "hydrolysate_support_active": bool(hydrolysate_support_active),
            "architecture_shift_active": bool(architecture_shift_active),
            "architecture_shift_score": float(architecture_shift_score),
            "mw_growth_fraction": float(mw_growth_fraction),
            "radius_of_gyration_growth_fraction": float(rg_growth_fraction),
            "mark_houwink_drop_fraction": float(mark_houwink_drop_fraction),
            "fft_fold_reduction_anchor": float(fft_fold_reduction_anchor),
            "architecture_reference_ids": [str(architecture_row.get("id", "unknown"))] if architecture_row else [],
            "architecture_reference_citations": [str(architecture_row.get("citation", "unknown"))] if architecture_row else [],
            "architecture_prior_ids": [str(architecture_prior.get("id", "unknown"))] if architecture_prior else [],
            "process_state_calibration_ids": [str(architecture_calibration.get("id", "unknown"))] if architecture_calibration else [],
            "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
            "retention_reference_citations": [str(row.get("source_citation", "unknown")) for row in retention_rows],
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
    protein_damage_lane: Mapping[str, Any],
    polyphenol_capping_lane: Mapping[str, Any],
    ascorbic_acid_lane: Mapping[str, Any],
    phospholipid_amine_lane: Mapping[str, Any],
    melanoidin_polymerization_lane: Mapping[str, Any],
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
            "protein_damage_burden",
            active=bool(protein_damage_lane.get("active", False)),
            state_value=float(protein_damage_lane.get("damage_burden_score", 0.0)),
            state_value_summary=(
                f"damage_burden_score={float(protein_damage_lane.get('damage_burden_score', 0.0)):.2f}, reactive_lysine_loss_pct={float(protein_damage_lane.get('reactive_lysine_loss_pct', 0.0)):.1f}, selected_anchor={protein_damage_lane.get('selected_benchmark_anchor_id', '') or 'none'}"
            ),
            influence_mode="guardrail_lane",
            family_lane=protein_damage_lane,
            summary=str(protein_damage_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "polyphenol_precursor_sink",
            active=bool(polyphenol_capping_lane.get("active", False)),
            state_value=float(polyphenol_capping_lane.get("cysteine_depletion_factor", 0.0)),
            state_value_summary=(
                f"quinone_budget={float(polyphenol_capping_lane.get('quinone_budget', 0.0)):.2f}, cysteine_depletion_factor={float(polyphenol_capping_lane.get('cysteine_depletion_factor', 0.0)):.2f}, lysine_depletion_factor={float(polyphenol_capping_lane.get('lysine_depletion_factor', 0.0)):.2f}"
            ),
            influence_mode="upstream_precursor_sink",
            family_lane=polyphenol_capping_lane,
            summary=str(polyphenol_capping_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "ascorbic_dicarbonyl_source",
            active=bool(ascorbic_acid_lane.get("active", False)),
            state_value=float(ascorbic_acid_lane.get("dicarbonyl_source_pressure", 0.0)),
            state_value_summary=(
                f"dicarbonyl_source_pressure={float(ascorbic_acid_lane.get('dicarbonyl_source_pressure', 0.0)):.2f}, ascorbic_acid_loss={float(ascorbic_acid_lane.get('ascorbic_acid_loss', 0.0)):.2f}, pentosidine_load={float(ascorbic_acid_lane.get('pentosidine_load', 0.0)):.2f}"
            ),
            influence_mode="bounded_upstream_source",
            family_lane=ascorbic_acid_lane,
            summary=str(ascorbic_acid_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "pe_glycation_fraction",
            active=bool(phospholipid_amine_lane.get("active", False)),
            state_value=float(phospholipid_amine_lane.get("pe_glycation_fraction", 0.0)),
            state_value_summary=(
                f"pe_glycation_fraction={float(phospholipid_amine_lane.get('pe_glycation_fraction', 0.0)):.2f}, sugar_sink_fraction={float(phospholipid_amine_lane.get('sugar_sink_fraction', 0.0)):.2f}, sugar_retention_factor={float(phospholipid_amine_lane.get('available_sugar_retention_factor', 1.0)):.2f}"
            ),
            influence_mode="upstream_precursor_sink",
            family_lane=phospholipid_amine_lane,
            summary=str(phospholipid_amine_lane.get("summary", "")),
        ),
        _build_state_marker_payload(
            "melanoidin_trapping_burden",
            active=bool(melanoidin_polymerization_lane.get("active", False)),
            state_value=float(melanoidin_polymerization_lane.get("matrix_trapping_pressure", 0.0)),
            state_value_summary=(
                f"melanoidin_mass={float(melanoidin_polymerization_lane.get('melanoidin_mass', 0.0)):.2f}, thiol_scavenging_factor={float(melanoidin_polymerization_lane.get('thiol_scavenging_factor', 0.0)):.2f}, fft_fold_reduction_anchor={float(melanoidin_polymerization_lane.get('fft_fold_reduction_anchor', 0.0)):.1f}"
            ),
            influence_mode="bounded_trapping_lane",
            family_lane=melanoidin_polymerization_lane,
            summary=str(melanoidin_polymerization_lane.get("summary", "")),
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


def _build_donor_hierarchy_lane(
    normalized_sugars: List[str],
    *,
    normalized_amino: Optional[List[str]] = None,
    normalized_additives: Optional[List[str]] = None,
    normalized_interventions: Optional[List[str]] = None,
    pH: Optional[float] = None,
) -> Dict[str, Any]:
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

    cue_pool = list(normalized_additives or []) + list(normalized_interventions or [])
    peptide_intensification_active = any(
        any(token in value for token in ["hydrolysate", "peptide", "yeast extract", "protease hydrolysis", "autolysate"])
        for value in cue_pool
    )
    normalized_amino = list(normalized_amino or [])
    prior_rows = query_family_runtime_priors(runtime_lane="pyrazine_control", family="07", process_state="heated_matrix")
    pyrazine_prior = dict(prior_rows[0]) if prior_rows else get_pyrazine_control_priors()
    anchor_rows = query_flavor_reference_entries(family="07", payload_section="pyrazine_reference_anchors")
    rhamnose_hdmf_prior = _computational_prior_entry("blank_1997_rhamnose_proline_hdmf_uplift_v1")
    rhamnose_active = any("rhamnose" in sugar for sugar in normalized_sugars)
    proline_active = any("proline" in value for value in normalized_amino)
    rhamnose_hdmf_reference_active = bool(rhamnose_active and proline_active and rhamnose_hdmf_prior)
    rhamnose_hdmf_molar_yield = float(rhamnose_hdmf_prior.get("hdmf_molar_yield_fraction_lower_bound", 0.0) or 0.0)
    rhamnose_hdmf_reference_score = 0.0
    if rhamnose_hdmf_reference_active:
        rhamnose_hdmf_reference_score = min(1.0, rhamnose_hdmf_molar_yield / 0.5)

    active = bool(normalized_sugars)
    donor_hierarchy_score = donor_strength.get(dominant_donor_class, 0.0)
    mixed_donor_system = len(donor_counts) > 1
    propensity_cap = float(pyrazine_prior.get("propensity_cap", 1.75) or 1.75)
    effective_ph = 6.2 if pH is None else float(pH)
    alkaline_pyrazine_factor = _sigmoid(effective_ph, 6.8, 0.7)
    low_ph_suppression_factor = 1.0 - alkaline_pyrazine_factor
    fructose_bias = 0.08 if "fructose" in donor_counts else 0.0
    if dominant_donor_class == "phosphorylated":
        fructose_bias = 0.04
    donor_hierarchy_score = min(
        1.0,
        donor_hierarchy_score * (0.92 + 0.12 * low_ph_suppression_factor),
    )
    pyrazine_propensity = donor_strength.get(dominant_donor_class, 0.0) * (0.70 + 0.85 * alkaline_pyrazine_factor + fructose_bias)
    if mixed_donor_system and dominant_donor_class in {"phosphorylated", "pentose"}:
        pyrazine_propensity += 0.08
    if peptide_intensification_active:
        pyrazine_propensity *= 1.12
    pyrazine_propensity = min(propensity_cap, pyrazine_propensity)
    pyrazine_pressure_score = min(1.0, pyrazine_propensity / max(propensity_cap, 1.0e-6))

    summary = "No explicit carbonyl donor hierarchy lane active."
    if active:
        summary = (
            f"Dominant donor class is {dominant_donor_class}, so sugar identity should be treated as a first-class routing variable rather than generic sugar loading."
        )
        if anchor_rows and prior_rows:
            summary = (
                f"Dominant donor class is {dominant_donor_class}, and donor identity is now constrained by pyrazine-control references plus a bounded prior instead of a fixed sugar heuristic."
            )
        if alkaline_pyrazine_factor > 0.55 and peptide_intensification_active:
            summary = (
                f"Dominant donor class is {dominant_donor_class}, but alkaline plus peptide-rich conditions move the lane into a pyrazine-pressure regime that must stay bounded by the literature control surface."
            )
        elif rhamnose_hdmf_reference_active:
            summary = (
                "Rhamnose and proline are both active, so the lane now exposes a bounded HDMF-support reference rather than treating this donor pair like a generic low-priority reducing sugar."
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
            "prior_ids": [str(row.get("id", "unknown")) for row in prior_rows] or ([str(pyrazine_prior.get("id", "unknown"))] if pyrazine_prior else []),
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in anchor_rows],
            "benchmark_anchor_citations": [str(row.get("source_citation", "unknown")) for row in anchor_rows],
            "peptide_intensification_active": bool(peptide_intensification_active),
            "rhamnose_hdmf_reference_active": bool(rhamnose_hdmf_reference_active),
            "rhamnose_hdmf_prior_ids": [str(rhamnose_hdmf_prior.get("id", "unknown"))] if rhamnose_hdmf_reference_active else [],
            "rhamnose_hdmf_reference_citations": [str(rhamnose_hdmf_prior.get("source", "unknown"))] if rhamnose_hdmf_reference_active else [],
            "rhamnose_hdmf_reference_score": float(rhamnose_hdmf_reference_score),
            "rhamnose_hdmf_molar_yield_fraction_lower_bound": float(rhamnose_hdmf_molar_yield),
            "rhamnose_hdmf_odt_ug_per_l": float(rhamnose_hdmf_prior.get("hdmf_odt_ug_per_l", 0.0) or 0.0),
            "rhamnose_hdmf_oav_lower_bound": float(rhamnose_hdmf_prior.get("hdmf_oav_lower_bound", 0.0) or 0.0),
            "alkaline_pyrazine_factor": float(alkaline_pyrazine_factor),
            "low_ph_suppression_factor": float(low_ph_suppression_factor),
            "pyrazine_pressure_score": float(pyrazine_pressure_score),
        },
    )


def _build_thiamine_support_lane(
    *,
    thiamine_active: bool,
    thiamine_fraction_estimate: float,
    thiamine_mode: str,
    thiamine_calibration: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    calibration = dict(thiamine_calibration or {})
    summary = "No explicit thiamine support lane active."
    if thiamine_active:
        summary = "Thiamine support is active, so sulfur routing should be interpreted as augmented rather than purely core-derived."
        if float(calibration.get("extrusion_survival_factor", 1.0)) < 0.95:
            summary = "Thiamine support is active, but the De Leyn 2019 extrusion-retention anchor sharply attenuates pre-extrusion donor survival under this process state."
        elif bool(calibration.get("aw_reference_active")) and float(calibration.get("aw_modulation_factor", 1.0)) < 0.97:
            summary = "Thiamine support is active, but the Arabshahi 1988 aw-dependent anchor attenuates support under wetter starch-like conditions relative to the mid-aw reference state."
        elif str(calibration.get("yield_mode", "")) == "mixed_system_optimal_window":
            summary = "Thiamine support is active and sits inside the Cerny 2008 mixed-system pH window, so MFT amplification should be treated as benchmark-calibrated rather than availability-only."
    return _build_family_lane_payload(
        "03",
        active=thiamine_active,
        summary=summary,
        metrics={
            "thiamine_support_score": float(thiamine_fraction_estimate),
            "thiamine_mode": thiamine_mode,
            "baseline_thiamine_fraction": float(calibration.get("baseline_fraction", thiamine_fraction_estimate)),
            "thiamine_pH_factor": float(calibration.get("pH_factor", 1.0)),
            "thiamine_reference_yield_ug_per_g": float(calibration.get("reference_yield_ug_per_g", 0.0)),
            "thiamine_reference_yield_mode": str(calibration.get("yield_mode", "")),
            "thiamine_synergy_factor": float(calibration.get("synergy_factor", 1.0)),
            "extrusion_survival_factor": float(calibration.get("extrusion_survival_factor", 1.0)),
            "extrusion_reference_id": str(calibration.get("extrusion_reference_id", "")),
            "extrusion_reference_source": str(calibration.get("extrusion_reference_source", "")),
            "thiamine_aw_reference_active": bool(calibration.get("aw_reference_active", False)),
            "thiamine_aw_reference_id": str(calibration.get("aw_reference_id", "")),
            "thiamine_aw_reference_source": str(calibration.get("aw_reference_source", "")),
            "thiamine_aw_reference_matrix": str(calibration.get("aw_reference_matrix", "")),
            "thiamine_aw_reference_value": float(calibration.get("aw_reference_value", 0.0)),
            "thiamine_aw_reference_ea_kcal_per_mol": float(calibration.get("aw_reference_ea_kcal_per_mol", 0.0)),
            "thiamine_aw_modulation_factor": float(calibration.get("aw_modulation_factor", 1.0)),
            "thiamine_aw_transfer_weight": float(calibration.get("aw_transfer_weight", 1.0)),
            "benchmark_anchor_ids": list(calibration.get("benchmark_anchor_ids", [])),
            "benchmark_anchor_citations": list(calibration.get("benchmark_anchor_citations", [])),
        },
    )


def _estimate_nucleotide_half_life_hours(
    nucleotide_priors: Mapping[str, Any],
    *,
    nucleotide_key: str,
    pH: Optional[float],
    temperature_celsius: Optional[float],
) -> Optional[float]:
    temperature_profile = nucleotide_priors.get(f"{nucleotide_key}_half_life_hours_at_ph7_by_temp_c", {})
    ph_profile = nucleotide_priors.get(f"{nucleotide_key}_half_life_hours_at_100c_by_ph", {})
    if not isinstance(temperature_profile, Mapping) or not isinstance(ph_profile, Mapping):
        return None
    base_half_life = _interpolate_numeric_mapping(temperature_profile, temperature_celsius if temperature_celsius is not None else 100.0)
    if base_half_life is None:
        return None
    ph_reference = _interpolate_numeric_mapping(ph_profile, 7.0) or 1.0
    ph_half_life = _interpolate_numeric_mapping(ph_profile, pH if pH is not None else 7.0)
    if ph_half_life is not None:
        base_half_life *= float(ph_half_life) / max(float(ph_reference), 1.0e-6)
    return float(base_half_life)


def _build_nucleotide_calibration_context(
    *,
    nucleotide_priors: Mapping[str, Any],
    normalized_sugars: List[str],
    normalized_additives: List[str],
    normalized_interventions: List[str],
    pH: Optional[float],
    process_state: Optional[str],
    temperature_celsius: Optional[float],
    time_minutes: Optional[float],
) -> Dict[str, Any]:
    benchmark_rows = query_benchmark_intake_entries(family="04", primary_only=True)
    soladoye_row = next(
        (row for row in benchmark_rows if str(row.get("id", "")).strip() == "soladoye_2020_sous_vide_euc_anchor"),
        {},
    )
    soladoye_prior = get_family_runtime_prior(
        runtime_lane="nucleotide_pathway",
        entry_id="soladoye_2020_low_temp_euc_window_v1",
    )
    yeast_extract_grade_prior = get_family_runtime_prior(
        runtime_lane="nucleotide_pathway",
        entry_id="ahlberg_2021_yeast_extract_nucleotide_grade_window_v1",
    )
    mushroom_profile_prior = get_family_runtime_prior(
        runtime_lane="nucleotide_pathway",
        entry_id="cui_2022_mushroom_gmp_euc_window_v1",
    )
    support_pool = normalized_sugars + normalized_additives + normalized_interventions
    # Word-boundary matched (2026-08-27): "imp" as a substring fired on
    # "important"/"imported"/"improver" in free-text cue lists.
    imp_active = any(_any_token_matches(value, ["imp", "inosinate"]) for value in support_pool)
    gmp_active = any(_any_token_matches(value, ["gmp", "guanylate"]) for value in support_pool)
    yeast_extract_active = any("yeast extract" in value for value in support_pool)
    mushroom_species = ""
    for candidate in ["shiitake", "porcini", "enoki", "oyster"]:
        if any(candidate in value for value in normalized_additives + normalized_interventions):
            mushroom_species = candidate
            break
    explicit_ribose_active = any(
        any(token in value for token in ["ribose", "ribose 5 phosphate", "ribose-5-phosphate", "r5p"])
        for value in normalized_sugars + normalized_additives
    )

    imp_share = 1.0 if imp_active else 0.0
    gmp_share = 1.0 if gmp_active else 0.0
    if yeast_extract_active and imp_share + gmp_share <= 0.0:
        imp_share = float(nucleotide_priors.get("yeast_extract_imp_share", 0.6) or 0.6)
        gmp_share = float(nucleotide_priors.get("yeast_extract_gmp_share", 0.4) or 0.4)
    share_total = imp_share + gmp_share

    default_duration = None
    if time_minutes is not None:
        default_duration = float(time_minutes)
    elif process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"}:
        default_duration = float(nucleotide_priors.get("default_extrusion_residence_minutes", 2.0) or 2.0)

    imp_survival = 1.0
    gmp_survival = 1.0
    if share_total > 0.0 and default_duration is not None and temperature_celsius is not None:
        duration_hours = max(float(default_duration), 0.0) / 60.0
        imp_half_life = _estimate_nucleotide_half_life_hours(
            nucleotide_priors,
            nucleotide_key="imp",
            pH=pH,
            temperature_celsius=temperature_celsius,
        )
        gmp_half_life = _estimate_nucleotide_half_life_hours(
            nucleotide_priors,
            nucleotide_key="gmp",
            pH=pH,
            temperature_celsius=temperature_celsius,
        )
        if imp_half_life is not None:
            imp_survival = max(0.0, min(1.0, 0.5 ** (duration_hours / max(float(imp_half_life), 1.0e-6))))
        if gmp_half_life is not None:
            gmp_survival = max(0.0, min(1.0, 0.5 ** (duration_hours / max(float(gmp_half_life), 1.0e-6))))

    nucleotide_survival_factor = 0.0
    if share_total > 0.0:
        nucleotide_survival_factor = (imp_share * imp_survival + gmp_share * gmp_survival) / share_total

    hydrolysis_fraction = max(0.0, 1.0 - nucleotide_survival_factor) if share_total > 0.0 else 0.0
    max_ribose_delivery = float(nucleotide_priors.get("max_ribose_delivery_fraction", 0.4) or 0.4)
    ribose_delivery_factor = 1.0 if explicit_ribose_active else min(
        max_ribose_delivery,
        hydrolysis_fraction * float(nucleotide_priors.get("ribose_release_efficiency_multiplier", 2.3) or 2.3),
    )
    umami_support_factor = float(nucleotide_survival_factor if share_total > 0.0 else 0.0)
    nucleotide_support_score = min(
        1.0,
        float(nucleotide_priors.get("umami_weight", 0.6) or 0.6) * umami_support_factor
        + float(nucleotide_priors.get("ribose_weight", 0.4) or 0.4) * ribose_delivery_factor,
    )
    ribose_shift_active = bool(
        share_total > 0.0
        and (
            ribose_delivery_factor > umami_support_factor
            or (hydrolysis_fraction >= 0.35 and ribose_delivery_factor >= 0.25)
        )
    )
    soladoye_key_values = soladoye_row.get("key_values", {}) if isinstance(soladoye_row.get("key_values", {}), Mapping) else {}
    soladoye_euc_window = soladoye_key_values.get("euc_percent_msg_by_condition", {}) if isinstance(soladoye_key_values.get("euc_percent_msg_by_condition", {}), Mapping) else {}
    yeast_extract_grade_reference_active = bool(yeast_extract_active and yeast_extract_grade_prior)
    mushroom_profiles = mushroom_profile_prior.get("species_profiles", {}) if isinstance(mushroom_profile_prior.get("species_profiles", {}), Mapping) else {}
    mushroom_profile = mushroom_profiles.get(mushroom_species, {}) if isinstance(mushroom_profiles.get(mushroom_species, {}), Mapping) else {}
    mushroom_reference_active = bool(mushroom_species and mushroom_profile)
    low_temp_euc_window_active = bool(
        soladoye_row
        and process_state == "heated_matrix"
        and temperature_celsius is not None
        and float(temperature_celsius) <= 75.0
        and default_duration is not None
        and float(default_duration) >= 120.0
        and share_total > 0.0
    )

    return {
        "active": bool(share_total > 0.0 or explicit_ribose_active),
        "imp_active": bool(imp_active),
        "gmp_active": bool(gmp_active),
        "yeast_extract_active": bool(yeast_extract_active),
        "yeast_extract_grade_reference_active": bool(yeast_extract_grade_reference_active),
        "yeast_extract_grade_prior_ids": [str(yeast_extract_grade_prior.get("id", "unknown"))] if yeast_extract_grade_reference_active else [],
        "yeast_extract_grade_reference_citations": [str(yeast_extract_grade_prior.get("source", "unknown"))] if yeast_extract_grade_reference_active else [],
        "yeast_extract_standard_imp_mg_per_100g_dw_range": list(yeast_extract_grade_prior.get("standard_autolysate_imp_mg_per_100g_dw_range", [])) if yeast_extract_grade_reference_active else [],
        "yeast_extract_high_nucleotide_imp_mg_per_100g_dw_range": list(yeast_extract_grade_prior.get("high_nucleotide_imp_mg_per_100g_dw_range", [])) if yeast_extract_grade_reference_active else [],
        "yeast_extract_high_nucleotide_gmp_mg_per_100g_dw_range": list(yeast_extract_grade_prior.get("high_nucleotide_gmp_mg_per_100g_dw_range", [])) if yeast_extract_grade_reference_active else [],
        "yeast_extract_amp_deaminase_euc_uplift_factor_vs_standard_range": list(yeast_extract_grade_prior.get("amp_deaminase_euc_uplift_factor_vs_standard_range", [])) if yeast_extract_grade_reference_active else [],
        "mushroom_reference_active": bool(mushroom_reference_active),
        "mushroom_prior_ids": [str(mushroom_profile_prior.get("id", "unknown"))] if mushroom_reference_active else [],
        "mushroom_reference_citations": [str(mushroom_profile_prior.get("source", "unknown"))] if mushroom_reference_active else [],
        "mushroom_selected_species": str(mushroom_species),
        "mushroom_gmp_mg_per_100g_dw": float(mushroom_profile.get("gmp_mg_per_100g_dw", 0.0) or 0.0),
        "mushroom_nucleotide_only_euc_g_msg_per_100g": float(mushroom_profile.get("nucleotide_only_euc_g_msg_per_100g", 0.0) or 0.0),
        "mushroom_total_euc_g_msg_per_100g": float(mushroom_profile.get("total_euc_g_msg_per_100g", 0.0) or 0.0),
        "explicit_ribose_active": bool(explicit_ribose_active),
        "nucleotide_support_active": bool(share_total > 0.0),
        "ribose_delivery_active": bool(explicit_ribose_active or ribose_delivery_factor > 0.0),
        "nucleotide_survival_factor": float(nucleotide_survival_factor),
        "imp_survival_factor": float(imp_survival if share_total > 0.0 else 0.0),
        "gmp_survival_factor": float(gmp_survival if share_total > 0.0 else 0.0),
        "umami_support_factor": float(umami_support_factor),
        "ribose_delivery_factor": float(ribose_delivery_factor),
        "hydrolysis_fraction": float(hydrolysis_fraction),
        "nucleotide_support_score": float(nucleotide_support_score),
        "ribose_shift_active": ribose_shift_active,
        "umami_reference_mode": (
            "preserved_nucleotide_pool"
            if share_total > 0.0 and not ribose_shift_active
            else "hydrolyzing_nucleotide_pool"
            if share_total > 0.0
            else "explicit_ribose_delivery"
            if explicit_ribose_active
            else "inactive"
        ),
        "imp_fold_reduction_at_3mm": float(nucleotide_priors.get("imp_fold_reduction_at_3mm", 45.2) or 45.2),
        "gmp_fold_reduction_at_3mm": float(nucleotide_priors.get("gmp_fold_reduction_at_3mm", 29.8) or 29.8),
        "umami_synergy_constant": float(nucleotide_priors.get("umami_synergy_constant", 1218.0) or 1218.0),
        "low_temp_euc_window_active": bool(low_temp_euc_window_active),
        "soladoye_reference_ids": [str(soladoye_row.get("id", "unknown"))] if soladoye_row else [],
        "soladoye_reference_citations": [str(soladoye_row.get("citation", "unknown"))] if soladoye_row else [],
        "soladoye_prior_ids": [str(soladoye_prior.get("id", "unknown"))] if soladoye_prior else [],
        "soladoye_raw_euc_percent_msg": float(soladoye_key_values.get("raw_euc_percent_msg", 0.0) or 0.0),
        "soladoye_euc_percent_msg_by_condition": {
            str(key): float(value)
            for key, value in soladoye_euc_window.items()
            if value is not None
        },
        "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in benchmark_rows],
        "benchmark_anchor_citations": [str(row.get("citation", "unknown")) for row in benchmark_rows],
    }


def _build_nucleotide_support_lane(
    *,
    normalized_sugars: List[str],
    normalized_additives: List[str],
    normalized_interventions: List[str],
    nucleotide_calibration: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    support_pool = normalized_sugars + normalized_additives + normalized_interventions
    nucleotide_tokens = ["imp", "gmp", "inosinate", "guanylate", "yeast extract"]
    nucleotide_active = any(_any_token_matches(value, nucleotide_tokens) for value in support_pool)
    ribose_delivery_active = any(_any_token_matches(value, ["ribose", "ribose 5 phosphate", "ribose-5-phosphate", "r5p"]) for value in normalized_sugars + normalized_additives)
    calibration = dict(nucleotide_calibration or {})
    active = bool(calibration.get("active", nucleotide_active or ribose_delivery_active))
    nucleotide_support_score = float(
        calibration.get(
            "nucleotide_support_score",
            min(1.0, 0.55 * float(nucleotide_active) + 0.45 * float(ribose_delivery_active)),
        )
    )
    summary = "No explicit nucleotide support lane active."
    if active:
        if bool(calibration.get("ribose_shift_active", False)):
            summary = "Thermal severity is shifting nucleotide support away from preserved IMP/GMP and toward ribose delivery for downstream Maillard chemistry."
        elif bool(calibration.get("low_temp_euc_window_active", False)):
            summary = "Low-temperature heating remains inside the Soladoye EUC-reference window, so preserved nucleotide support should stay visible instead of being treated as an automatic ribose-shift case."
        elif bool(calibration.get("yeast_extract_grade_reference_active", False)) or bool(calibration.get("mushroom_reference_active", False)):
            summary = "Nucleotide support is active, and the lane now carries ingredient-source references for yeast-extract grade and mushroom GMP content so formulation-facing EUC context stays explicit."
        elif bool(calibration.get("explicit_ribose_active", False)) and not bool(calibration.get("nucleotide_support_active", nucleotide_active)):
            summary = "Explicit ribose delivery is active, so Family 04 is supporting downstream chemistry without preserved nucleotide-driven umami support."
        else:
            summary = "Nucleotide support is active, so upstream umami amplification and bounded ribose delivery should both be carried into prediction context."
    return _build_family_lane_payload(
        "04",
        active=active,
        summary=summary,
        metrics={
            "nucleotide_support_active": bool(calibration.get("nucleotide_support_active", nucleotide_active)),
            "ribose_delivery_active": bool(calibration.get("ribose_delivery_active", ribose_delivery_active)),
            "nucleotide_support_score": float(nucleotide_support_score),
            "explicit_ribose_active": bool(calibration.get("explicit_ribose_active", ribose_delivery_active)),
            "yeast_extract_grade_reference_active": bool(calibration.get("yeast_extract_grade_reference_active", False)),
            "yeast_extract_grade_prior_ids": list(calibration.get("yeast_extract_grade_prior_ids", [])),
            "yeast_extract_grade_reference_citations": list(calibration.get("yeast_extract_grade_reference_citations", [])),
            "yeast_extract_standard_imp_mg_per_100g_dw_range": list(calibration.get("yeast_extract_standard_imp_mg_per_100g_dw_range", [])),
            "yeast_extract_high_nucleotide_imp_mg_per_100g_dw_range": list(calibration.get("yeast_extract_high_nucleotide_imp_mg_per_100g_dw_range", [])),
            "yeast_extract_high_nucleotide_gmp_mg_per_100g_dw_range": list(calibration.get("yeast_extract_high_nucleotide_gmp_mg_per_100g_dw_range", [])),
            "yeast_extract_amp_deaminase_euc_uplift_factor_vs_standard_range": list(calibration.get("yeast_extract_amp_deaminase_euc_uplift_factor_vs_standard_range", [])),
            "mushroom_reference_active": bool(calibration.get("mushroom_reference_active", False)),
            "mushroom_prior_ids": list(calibration.get("mushroom_prior_ids", [])),
            "mushroom_reference_citations": list(calibration.get("mushroom_reference_citations", [])),
            "mushroom_selected_species": str(calibration.get("mushroom_selected_species", "")),
            "mushroom_gmp_mg_per_100g_dw": float(calibration.get("mushroom_gmp_mg_per_100g_dw", 0.0)),
            "mushroom_nucleotide_only_euc_g_msg_per_100g": float(calibration.get("mushroom_nucleotide_only_euc_g_msg_per_100g", 0.0)),
            "mushroom_total_euc_g_msg_per_100g": float(calibration.get("mushroom_total_euc_g_msg_per_100g", 0.0)),
            "nucleotide_survival_factor": float(calibration.get("nucleotide_survival_factor", float(nucleotide_active))),
            "imp_survival_factor": float(calibration.get("imp_survival_factor", 0.0)),
            "gmp_survival_factor": float(calibration.get("gmp_survival_factor", 0.0)),
            "umami_support_factor": float(calibration.get("umami_support_factor", float(nucleotide_active))),
            "ribose_delivery_factor": float(calibration.get("ribose_delivery_factor", float(ribose_delivery_active))),
            "hydrolysis_fraction": float(calibration.get("hydrolysis_fraction", 0.0)),
            "ribose_shift_active": bool(calibration.get("ribose_shift_active", False)),
            "umami_reference_mode": str(calibration.get("umami_reference_mode", "inactive")),
            "imp_fold_reduction_at_3mm": float(calibration.get("imp_fold_reduction_at_3mm", 0.0)),
            "gmp_fold_reduction_at_3mm": float(calibration.get("gmp_fold_reduction_at_3mm", 0.0)),
            "umami_synergy_constant": float(calibration.get("umami_synergy_constant", 0.0)),
            "low_temp_euc_window_active": bool(calibration.get("low_temp_euc_window_active", False)),
            "soladoye_reference_ids": list(calibration.get("soladoye_reference_ids", [])),
            "soladoye_reference_citations": list(calibration.get("soladoye_reference_citations", [])),
            "soladoye_prior_ids": list(calibration.get("soladoye_prior_ids", [])),
            "soladoye_raw_euc_percent_msg": float(calibration.get("soladoye_raw_euc_percent_msg", 0.0)),
            "soladoye_euc_percent_msg_by_condition": dict(calibration.get("soladoye_euc_percent_msg_by_condition", {})),
            "benchmark_anchor_ids": list(calibration.get("benchmark_anchor_ids", [])),
            "benchmark_anchor_citations": list(calibration.get("benchmark_anchor_citations", [])),
        },
    )


def _build_sulfur_peptide_support_lane(
    *,
    normalized_additives: List[str],
    normalized_amino: List[str],
    normalized_support_cues: List[str],
    normalized_interventions: List[str],
    protein_label: str,
    degree_of_hydrolysis: Optional[float] = None,
) -> Dict[str, Any]:
    support_pool = normalized_additives + normalized_amino + normalized_support_cues + normalized_interventions
    glutathione_active = any("glutathione" in value for value in support_pool)
    hydrolysate_active = any("hydrolysate" in value for value in support_pool)
    autolysate_active = any("autolysate" in value for value in support_pool)
    hydrolysis_fraction = None if degree_of_hydrolysis is None else max(0.0, min(1.0, float(degree_of_hydrolysis)))
    peptide_support_active = any(any(token in value for token in ["peptide", "hydrolysate", "autolysate"]) for value in support_pool)
    if hydrolysis_fraction is not None and hydrolysis_fraction < 0.99:
        peptide_support_active = True
    if hydrolysis_fraction is not None and hydrolysis_fraction >= 0.55:
        hydrolysate_active = True
    sulfur_peptide_priors = get_sulfur_peptide_priors()
    kokumi_priors = get_family_runtime_prior(
        runtime_lane="sulfur_peptide_support",
        entry_id="ohsu_2025_kokumi_casr_support_v1",
    )
    benchmark_rows = query_benchmark_intake_entries(
        family="05",
        primary_only=True,
        process_state="heated_matrix",
    )
    kokumi_reference_rows = query_benchmark_intake_entries(
        entry_id="ohsu_2025_kokumi_casr_anchor",
    )
    retention_rows = query_retention_reference_entries(
        family="05",
        protein_type=protein_label,
        process_state="heated_matrix",
    )
    glutathione_support_ratio = 0.5 * (
        float(sulfur_peptide_priors.get("gsh_mft_ratio_vs_free_cysteine", 2.25) or 2.25)
        + float(sulfur_peptide_priors.get("gsh_fft_ratio_vs_free_cysteine", 2.40) or 2.40)
    )
    generic_peptide_ratio = float(sulfur_peptide_priors.get("generic_peptide_ratio_vs_free_cysteine", 1.41) or 1.41)
    hydrolysate_ratio = float(sulfur_peptide_priors.get("hydrolysate_ratio_vs_free_cysteine", 1.60) or 1.60)
    peak_sequence_ratio = float(sulfur_peptide_priors.get("peak_sequence_ratio_vs_free_cysteine", 2.68) or 2.68)
    buffered_release_bonus = float(sulfur_peptide_priors.get("gsh_buffered_release_bonus", 0.08) or 0.08)
    pyrazine_tradeoff_ratio = float(sulfur_peptide_priors.get("pyrazine_tradeoff_ratio_vs_free_cysteine", 0.75) or 0.75)
    selected_peptide_ratio = 1.0
    peptide_mode = "inactive"
    if peptide_support_active:
        if hydrolysate_active or autolysate_active:
            selected_peptide_ratio = hydrolysate_ratio
            peptide_mode = "hydrolysate_supported"
        else:
            selected_peptide_ratio = generic_peptide_ratio
            peptide_mode = "generic_peptide"
    if hydrolysis_fraction is not None and peptide_support_active and not glutathione_active:
        if hydrolysate_active or autolysate_active:
            selected_peptide_ratio = 1.0 + (hydrolysate_ratio - 1.0) * hydrolysis_fraction
        else:
            selected_peptide_ratio = 1.0 + (generic_peptide_ratio - 1.0) * hydrolysis_fraction

    peptide_accessibility_factor = 1.0
    if hydrolysis_fraction is not None:
        peptide_accessibility_factor = 0.30 + 0.70 * hydrolysis_fraction
    elif hydrolysate_active or autolysate_active:
        peptide_accessibility_factor = 0.92
    elif peptide_support_active:
        peptide_accessibility_factor = 0.78

    selected_support_ratio = max(
        glutathione_support_ratio if glutathione_active else 1.0,
        selected_peptide_ratio,
    )
    free_cysteine_equivalent_factor = selected_support_ratio * peptide_accessibility_factor
    sulfur_peptide_support_score = 0.0
    if glutathione_active or peptide_support_active:
        sulfur_peptide_support_score = max(
            0.0,
            min(
                1.0,
                (selected_support_ratio - 1.0) / max(peak_sequence_ratio - 1.0, 1.0e-6),
            ),
        )
        if glutathione_active and peptide_support_active:
            sulfur_peptide_support_score = min(1.0, sulfur_peptide_support_score + buffered_release_bonus)
    sulfur_proxy_retention_factor = max(
        [
            float(row.get("numeric_reference", {}).get("value", 0.0)) / 100.0
            for row in retention_rows
            if isinstance(row.get("numeric_reference", {}), Mapping)
            and str(row.get("numeric_reference", {}).get("units", "")).strip() == "percent_retained"
        ]
        or [0.0]
    )
    gamma_glutamyl_peptide_active = bool(hydrolysate_active or autolysate_active)
    kokumi_profile = build_kokumi_support_profile(
        glutathione_active=glutathione_active,
        gamma_glutamyl_peptide_active=gamma_glutamyl_peptide_active,
        peptide_accessibility_factor=peptide_accessibility_factor,
        kokumi_priors=kokumi_priors,
    )
    kokumi_support_active = bool(kokumi_profile.get("kokumi_support_active", False))
    summary = "No explicit glutathione or peptide support lane active."
    if glutathione_active or peptide_support_active:
        summary = "Glutathione or peptide support is active, so sulfur intensity should be interpreted as a bounded GSH or hydrolysate-assisted lane rather than a purely free-cysteine signal."
        if benchmark_rows:
            summary = "Glutathione or peptide support is active, and the lane is anchored to the existing soy-hydrolysate intake record while using bounded GSH and peptide uplift priors for sulfur-positive interpretation."
    if kokumi_support_active and kokumi_reference_rows:
        summary = (
            f"{summary} The same Family 05 cues also activate a bounded kokumi-support lane tied to the Ohsu CaSR reference, so mouthfulness support should be reported alongside volatile sulfur output."
        )
    prior_ids = [str(sulfur_peptide_priors.get("id", ""))] if sulfur_peptide_priors else []
    kokumi_prior_ids = [str(kokumi_priors.get("id", ""))] if kokumi_priors else []
    return _build_family_lane_payload(
        "05",
        active=bool(glutathione_active or peptide_support_active),
        summary=summary,
        metrics={
            "glutathione_active": bool(glutathione_active),
            "hydrolysate_active": bool(hydrolysate_active),
            "autolysate_active": bool(autolysate_active),
            "peptide_support_active": bool(peptide_support_active),
            "sulfur_peptide_support_score": float(sulfur_peptide_support_score),
            "glutathione_support_ratio": float(glutathione_support_ratio if glutathione_active else 1.0),
            "selected_peptide_ratio": float(selected_peptide_ratio),
            "peak_sequence_ratio": float(peak_sequence_ratio),
            "peptide_mode": peptide_mode,
            "degree_of_hydrolysis": hydrolysis_fraction,
            "peptide_accessibility_factor": float(peptide_accessibility_factor),
            "free_cysteine_equivalent_factor": float(free_cysteine_equivalent_factor),
            "pyrazine_tradeoff_ratio_vs_free_cysteine": float(pyrazine_tradeoff_ratio),
            "sulfur_proxy_retention_factor": float(sulfur_proxy_retention_factor),
            "gamma_glutamyl_peptide_active": bool(gamma_glutamyl_peptide_active),
            **kokumi_profile,
            "kokumi_prior_ids": kokumi_prior_ids,
            "kokumi_reference_ids": [str(row.get("id", "unknown")) for row in kokumi_reference_rows],
            "kokumi_reference_citations": [str(row.get("citation", "unknown")) for row in kokumi_reference_rows],
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in benchmark_rows],
            "benchmark_anchor_citations": [str(row.get("citation", "unknown")) for row in benchmark_rows],
            "retention_reference_ids": [str(row.get("id", "unknown")) for row in retention_rows],
            "prior_ids": prior_ids,
        },
    )


def _build_matrix_scope_lane(*, protein_label: str) -> Dict[str, Any]:
    normalized_protein = _normalize_name(protein_label)
    source_id = _resolve_protein_source_id(protein_label)
    source_profile = _protein_source_profile(source_id)
    alternative_matrix_active = bool(normalized_protein and normalized_protein not in {"free", "pea iso", "pea conc", "soy iso", "soy conc"})
    quantitative_benchmarks = query_benchmark_intake_entries(
        family="06",
        primary_only=True,
        kind="quantitative_benchmark",
        process_state="heated_matrix",
    )
    process_state_rows = query_benchmark_intake_entries(
        family="06",
        primary_only=True,
        kind="calibration_reference",
    )
    exact_benchmark_row = next(
        (
            row for row in quantitative_benchmarks
            if source_id and source_id in str(row.get("matrix_family", "")).strip().lower()
        ),
        {},
    )
    transferable_benchmarks = [
        row
        for row in quantitative_benchmarks
        if "benchmark_validation" in {
            str(module).strip() for module in row.get("target_modules", []) if str(module).strip()
        }
    ]
    if exact_benchmark_row and not any(
        str(row.get("id", "")).strip() == str(exact_benchmark_row.get("id", "")).strip()
        for row in transferable_benchmarks
    ):
        transferable_benchmarks = [exact_benchmark_row, *transferable_benchmarks]
    structural_gaps = [dict(entry) for entry in iter_family_process_gap_entries(family="06")]
    meaty_potential_multiplier = float(source_profile.get("meaty_potential_multiplier", 0.0) or 0.0)
    hydrolysate_observability_bias = float(source_profile.get("hydrolysate_observability_bias", 0.0) or 0.0)
    off_note_penalty = float(source_profile.get("off_note_penalty", 0.0) or 0.0)
    lox_activity_flag = bool(source_profile.get("lox_activity_flag", False))
    matrix_source_support_score = max(
        0.0,
        min(1.0, meaty_potential_multiplier * (1.0 - 0.5 * off_note_penalty)),
    )
    nearest_benchmark_row: Dict[str, Any] = {}
    if not exact_benchmark_row and transferable_benchmarks and source_profile:
        nearest_benchmark_row = min(
            transferable_benchmarks,
            key=lambda row: abs(
                float(
                    _protein_source_profile(
                        "wheat_gluten" if "wheat_gluten" in str(row.get("matrix_family", "")).strip().lower() else "soy_isolate"
                    ).get("meaty_potential_multiplier", 0.0)
                ) - meaty_potential_multiplier
            ),
        )
    selected_benchmark_row = exact_benchmark_row if exact_benchmark_row else nearest_benchmark_row
    if not selected_benchmark_row and transferable_benchmarks:
        selected_benchmark_row = transferable_benchmarks[0]
    if not selected_benchmark_row and quantitative_benchmarks:
        selected_benchmark_row = quantitative_benchmarks[0]
    benchmark_transfer_mode = "inactive"
    if selected_benchmark_row:
        benchmark_transfer_mode = "exact_source_benchmark" if exact_benchmark_row else "nearest_source_transfer"
    matrix_uncertainty_factor = 0.0
    if alternative_matrix_active:
        matrix_uncertainty_factor = 0.20
        matrix_uncertainty_factor += 0.45 * max(0.0, 1.0 - hydrolysate_observability_bias)
        matrix_uncertainty_factor += 0.20 * off_note_penalty
        matrix_uncertainty_factor += 0.10 * float(lox_activity_flag)
        matrix_uncertainty_factor += 0.10 * float(bool(structural_gaps))
        if benchmark_transfer_mode == "exact_source_benchmark":
            matrix_uncertainty_factor -= 0.10
        elif benchmark_transfer_mode == "nearest_source_transfer":
            matrix_uncertainty_factor -= 0.03
        if process_state_rows:
            matrix_uncertainty_factor -= 0.05
        matrix_uncertainty_factor = max(0.20, min(0.85, matrix_uncertainty_factor))
    process_state_transfer_confidence = 1.0 - matrix_uncertainty_factor if alternative_matrix_active else 0.0
    summary = "No alternative-matrix scope lane active."
    if alternative_matrix_active:
        summary = "Alternative protein matrix scope is active, so source-specific meaty potential and process-state transfer uncertainty should be treated as first-class recommendation constraints."
        if benchmark_transfer_mode == "exact_source_benchmark":
            summary = "Alternative protein matrix scope is active with a direct source-matched benchmark anchor, but the unresolved same-experiment SH plus OPA plus DSC gap still bounds transfer confidence."
    return _build_family_lane_payload(
        "06",
        active=alternative_matrix_active,
        summary=summary,
        metrics={
            "matrix_scope_active": bool(alternative_matrix_active),
            "matrix_uncertainty_factor": float(matrix_uncertainty_factor),
            "protein_type": protein_label,
            "source_id": source_id,
            "source_profile_available": bool(source_profile),
            "meaty_potential_multiplier": float(meaty_potential_multiplier),
            "hydrolysate_observability_bias": float(hydrolysate_observability_bias),
            "off_note_penalty": float(off_note_penalty),
            "lox_activity_flag": bool(lox_activity_flag),
            "matrix_source_support_score": float(matrix_source_support_score),
            "process_state_transfer_confidence": float(process_state_transfer_confidence),
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in transferable_benchmarks],
            "matrix_reference_anchor_ids": [str(row.get("id", "unknown")) for row in quantitative_benchmarks],
            "selected_benchmark_anchor_id": str(selected_benchmark_row.get("id", "")) if selected_benchmark_row else "",
            "selected_benchmark_context": _build_benchmark_anchor_context(selected_benchmark_row) if selected_benchmark_row else {},
            "benchmark_transfer_mode": benchmark_transfer_mode,
            "process_state_anchor_ids": [str(row.get("id", "unknown")) for row in process_state_rows],
            "structural_gap_ids": [str(row.get("gap_id", "unknown")) for row in structural_gaps],
            # 2026-08-27 (Wave T3, T1-01): the five source-profile numbers above are
            # MOCKED. This block ships the registry's own admission alongside the
            # numbers it contaminates, so a consumer of matrix_uncertainty_factor
            # cannot read source differentiation as evidence.
            "protein_source_provenance": dict(PROTEIN_SOURCE_PROVENANCE),
            "protein_source_profile_unsourced": bool(PROTEIN_SOURCE_REGISTRY_UNSOURCED),
        },
    )


def _build_fermentation_pretreatment_lane(
    *,
    normalized_additives: List[str],
    normalized_amino: List[str],
    normalized_support_cues: List[str],
    normalized_interventions: List[str],
    pH: Optional[float],
    thiamine_active: bool,
) -> Dict[str, Any]:
    pretreatment_pool = normalized_additives + normalized_amino + normalized_support_cues + normalized_interventions
    # Word-boundary matched (2026-08-27). "ferment"/"culture" keep matching
    # their morphological variants (fermentation, fermented, cultured) via the
    # stem-token list; "miso"/"koji" no longer fire on longer words.
    fermentation_active = any(
        _any_token_matches(value, ["ferment", "koji", "miso", "starter culture", "culture", "yeast extract"])
        for value in pretreatment_pool
    )
    precursor_release_active = any(
        _any_token_matches(value, ["hydrolysate", "hydroly", "protease", "peptide", "autolysate"])
        for value in pretreatment_pool
    )
    off_note_cleanup_active = any(
        _any_token_matches(value, ["yeast fermentation", "yeast extract", "koji", "ferment"])
        for value in pretreatment_pool
    )
    nucleotide_support_active = any(
        _any_token_matches(value, ["imp", "gmp", "inosinate", "guanylate", "yeast extract"])
        for value in pretreatment_pool
    )
    benchmark_rows = query_benchmark_intake_entries(family="10")
    prior_rows = query_family_runtime_priors(family="10", process_state="heated_matrix")
    flavor_rows = query_flavor_reference_entries(family="10", payload_section="sulfur_reference_anchors")
    hydrolysate_cue_active = any(
        _any_token_matches(value, ["hydrolysate", "protease", "peptide", "autolysate"])
        for value in pretreatment_pool
    )
    explicit_nucleotide_cue = any(
        _any_token_matches(value, ["imp", "gmp", "inosinate", "guanylate"])
        for value in pretreatment_pool
    )
    pretreatment_pH_shift_active = bool((fermentation_active or off_note_cleanup_active) and pH is not None and float(pH) < 6.2)
    active = bool(
        fermentation_active
        or precursor_release_active
        or off_note_cleanup_active
        or nucleotide_support_active
    )
    nishimura_row = next((row for row in benchmark_rows if str(row.get("id", "")) == "nishimura_abe_2024"), {})
    matoba_row = next((row for row in benchmark_rows if str(row.get("id", "")) == "matoba_1988_nucleotide_hydrolysis"), {})
    selected_benchmark_row: Dict[str, Any] = {}
    if hydrolysate_cue_active and nishimura_row:
        selected_benchmark_row = nishimura_row
    elif (nucleotide_support_active or explicit_nucleotide_cue) and matoba_row:
        selected_benchmark_row = matoba_row
    elif nishimura_row:
        selected_benchmark_row = nishimura_row
    elif matoba_row:
        selected_benchmark_row = matoba_row
    elif benchmark_rows:
        selected_benchmark_row = dict(benchmark_rows[0])

    free_amino_acid_enrichment_score = min(
        1.0,
        0.55 * float(precursor_release_active)
        + 0.25 * float(hydrolysate_cue_active)
        + 0.20 * float(bool(nishimura_row)),
    )
    fermentation_cleanup_score = min(
        1.0,
        0.45 * float(off_note_cleanup_active)
        + 0.25 * float(fermentation_active)
        + 0.15 * float(pretreatment_pH_shift_active)
        + 0.15 * float(bool(flavor_rows)),
    )
    nucleotide_enrichment_score = min(
        1.0,
        0.50 * float(nucleotide_support_active)
        + 0.20 * float(explicit_nucleotide_cue)
        + 0.30 * float(bool(matoba_row)),
    )
    evidence_support_score = min(
        1.0,
        0.40 * min(len(benchmark_rows), 2) / 2.0
        + 0.35 * min(len(prior_rows), 3) / 3.0
        + 0.25 * min(len(flavor_rows), 2) / 2.0,
    )
    pretreatment_support_score = min(
        1.0,
        (
            0.38 * free_amino_acid_enrichment_score
            + 0.28 * fermentation_cleanup_score
            + 0.22 * nucleotide_enrichment_score
            + 0.12 * float(thiamine_active)
        )
        * (0.80 + 0.20 * evidence_support_score),
    )
    summary = "No explicit fermentation pretreatment lane active."
    if active:
        summary = (
            "Pretreatment cues indicate upstream fermentation or hydrolysis support, so precursor loading and off-note burden should be interpreted before thermal chemistry."
        )
        if selected_benchmark_row:
            summary = (
                "Pretreatment cues indicate upstream fermentation or hydrolysis support, and the lane is now bounded by literature anchors for hydrolysate chemistry or nucleotide release instead of keyword matching alone."
            )
    ordered_benchmark_rows: List[Dict[str, Any]] = []
    if selected_benchmark_row:
        ordered_benchmark_rows.append(dict(selected_benchmark_row))
    if matoba_row and str(matoba_row.get("id", "")) != str(selected_benchmark_row.get("id", "")):
        ordered_benchmark_rows.append(dict(matoba_row))
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
            "free_amino_acid_enrichment_score": float(free_amino_acid_enrichment_score),
            "fermentation_cleanup_score": float(fermentation_cleanup_score),
            "nucleotide_enrichment_score": float(nucleotide_enrichment_score),
            "evidence_support_score": float(evidence_support_score),
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in ordered_benchmark_rows],
            "benchmark_anchor_citations": [str(row.get("citation", "unknown")) for row in ordered_benchmark_rows],
            "selected_benchmark_anchor_id": str(selected_benchmark_row.get("id", "")) if selected_benchmark_row else "",
            "selected_benchmark_context": _build_benchmark_anchor_context(selected_benchmark_row) if selected_benchmark_row else {},
            "prior_ids": [str(row.get("id", "unknown")) for row in prior_rows],
            "flavor_anchor_ids": [str(row.get("id", "unknown")) for row in flavor_rows],
            "flavor_anchor_citations": [str(row.get("source_citation", "unknown")) for row in flavor_rows],
        },
    )


def _build_caramelization_lane(
    *,
    signal_map: Dict[str, float],
    furanone_expected: List[str],
    furanone_observed: List[str],
    pH: Optional[float],
) -> Dict[str, Any]:
    prior_rows = query_family_runtime_priors(runtime_lane="furanone_support", family="09", process_state="heated_matrix")
    default_prior = get_furanone_priors()
    mgo_hdmf_prior = _computational_prior_entry("brands_2002_mgo_hdmf_c3_route_v1")
    carbonyl_anchor_rows = query_flavor_reference_entries(family="09", payload_section="carbonyl_reference_anchors")
    furanone_anchor_rows = query_flavor_reference_entries(family="09", payload_section="furanone_reference_anchors")
    severity_targets = iter_target_panel_entries(
        chemistry_family="carbohydrate_pyrolysis_and_caramelization",
        target_class="severity_markers",
    )
    severity_marker_signals = {
        "furfural": _compound_signal(signal_map, ["furfural"]),
        "hmf": _compound_signal(signal_map, ["5-hydroxymethylfurfural", "5-hydroxymethylfurfural (hmf)", "hmf"]),
        "2_acetylfuran": _compound_signal(signal_map, ["2-acetylfuran", "2 acetylfuran"]),
    }
    furfural_anchor = next(
        (row for row in carbonyl_anchor_rows if str(row.get("id", "")) == "hernandez_2023_furfural_ratio_anchor"),
        {},
    )
    furfural_anchor_payload = furfural_anchor.get("numeric_band_or_point", {}) if isinstance(furfural_anchor.get("numeric_band_or_point", {}), Mapping) else {}
    furfural_meat_reference_ppb = float(furfural_anchor_payload.get("meat_reference", 0.0) or 0.0)
    furfural_reference_ratio = 0.0
    if furfural_meat_reference_ppb > 0.0 and severity_marker_signals["furfural"] > 0.0:
        furfural_reference_ratio = float(severity_marker_signals["furfural"]) / furfural_meat_reference_ppb
    severity_signal_ppb = sum(float(value) for value in severity_marker_signals.values())
    supportive_furanones_observed = len(furanone_observed)
    active = bool(severity_signal_ppb > 0.0 or furanone_expected)
    fragmentation_bias = _sigmoid(float(pH if pH is not None else 6.0), 6.55, 0.35)
    severity_bias = min(1.0, severity_signal_ppb / 40.0)
    mgo_hdmf_fragmentation_context_score = min(1.0, 0.65 * fragmentation_bias + 0.35 * severity_bias)
    mgo_hdmf_reference_active = bool(mgo_hdmf_prior and active and mgo_hdmf_fragmentation_context_score >= 0.45)
    mgo_hdmf_reference_score = float(mgo_hdmf_fragmentation_context_score if mgo_hdmf_reference_active else 0.0)
    caramelization_support_score = min(
        1.0,
        0.30 * float(bool(furanone_expected))
        + 0.20 * min(supportive_furanones_observed, 2) / 2.0
        + 0.25 * min(len(prior_rows), 1)
        + 0.25 * min(len(furanone_anchor_rows), 2) / 2.0,
    )
    severity_penalty_factor = min(1.5, max(severity_signal_ppb / 60.0, furfural_reference_ratio / 3.0))
    summary = "No explicit carbohydrate pyrolysis or caramelization lane active."
    if active:
        summary = "Carbohydrate pyrolysis or caramelization markers are active, so helpful browning support must be separated from over-furan drift."
        if prior_rows or carbonyl_anchor_rows or furanone_anchor_rows:
            summary = "Carbohydrate pyrolysis or caramelization markers are active, and the lane now carries explicit furanone and carbonyl anchors so browning support is separated from benchmark-visible over-furan drift."
        if mgo_hdmf_reference_active:
            # CLAIM CORRECTED 2026-08-29 (Wave Q1). This string used to say the
            # lane "carries an explicit methylglyoxal-to-HDMF C3 anchor". Two
            # things were wrong with that. (1) ANCHOR overstates what this lane
            # holds: the only thing behind this branch is the computational prior
            # `brands_2002_mgo_hdmf_c3_route_v1`, which the repo itself flags
            # `source_anchor_status: unanchored_no_doi_field` -- no record in any
            # data/lit/*.json resolves it to a paper. (2) THE LANE does not carry
            # a reaction edge at all; this module is the family-lane NARRATIVE
            # layer and is not wired to the kinetic core. A real MGO -> DMHF edge
            # DOES now exist, but it lives in the core (`r_mgo_dmhf`, 2 MGO ->
            # DMHF, rate key `k_mgo_dmhf`, added by B7 from Wang & Ho 2008 JAFC
            # 56:7405-7409 Fig. 1), and its own level is flagged
            # `digitised_from_a_bar_chart` / `single_temperature_no_ea_licensed`.
            # The sentence now says which of the two the reader is looking at.
            summary = (
                "Carbohydrate pyrolysis or caramelization markers are active, and "
                "an unanchored methylglyoxal-to-HDMF C3 prior is in scope, so "
                "fragmentation-heavy browning is not treated as amino-acid-dependent "
                "only. This lane carries the PRIOR, not a rate: the prior resolves to "
                "no DOI, and the model's actual C3+C3 edge is the kinetic core's "
                "`r_mgo_dmhf` (Wang & Ho 2008), whose level is a digitised bar-chart "
                "value at a single temperature."
            )
    return _build_family_lane_payload(
        "09",
        active=active,
        summary=summary,
        metrics={
            "severity_marker_signals_ppb": severity_marker_signals,
            "severity_signal_ppb": float(severity_signal_ppb),
            "caramelization_support_score": float(caramelization_support_score),
            "severity_penalty_factor": float(severity_penalty_factor),
            "furanone_support_active": bool(furanone_expected or furanone_observed),
            "prior_ids": [str(default_prior.get("id", "unknown"))] if default_prior else [],
            "mgo_hdmf_reference_active": bool(mgo_hdmf_reference_active),
            "mgo_hdmf_prior_ids": [str(mgo_hdmf_prior.get("id", "unknown"))] if mgo_hdmf_reference_active else [],
            "mgo_hdmf_reference_citations": [str(mgo_hdmf_prior.get("source", "unknown"))] if mgo_hdmf_reference_active else [],
            "mgo_hdmf_fragmentation_context_score": float(mgo_hdmf_fragmentation_context_score),
            "mgo_hdmf_reference_score": float(mgo_hdmf_reference_score),
            "mgo_hdmf_reference_yield_ug_per_g": float(mgo_hdmf_prior.get("hdmf_reference_yield_ug_per_g", 0.0) or 0.0),
            "mgo_hdmf_reference_mgo_mM": float(mgo_hdmf_prior.get("reference_conditions", {}).get("mgo_mM", 0.0) or 0.0) if isinstance(mgo_hdmf_prior.get("reference_conditions", {}), Mapping) else 0.0,
            "mgo_hdmf_route_independent_of_amino_acids": bool(mgo_hdmf_prior.get("route_independent_of_amino_acids", False)),
            "carbonyl_anchor_ids": [str(row.get("id", "unknown")) for row in carbonyl_anchor_rows],
            "furanone_anchor_ids": [str(row.get("id", "unknown")) for row in furanone_anchor_rows],
            "severity_target_ids": [str(row.get("canonical_name", "unknown")) for row in severity_targets],
            "furfural_meat_reference_ppb": float(furfural_meat_reference_ppb),
            "furfural_reference_ratio": float(furfural_reference_ratio),
        },
    )


def _build_off_note_guardrail_lane(
    *,
    signal_map: Dict[str, float],
    normalized_additives: List[str],
    normalized_interventions: List[str],
    normalized_sugars: List[str],
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
    cyclodextrin_guardrail_active = any("cyclodextrin" in value for value in guardrail_pool)
    lipid_marker_signal_ppb = float(lipid_crosstalk_lane.get("lipid_marker_signal_ppb", 0.0))
    benchmark_rows = query_benchmark_intake_entries(family="08", primary_only=True)
    crosstalk_prior_rows = query_family_runtime_priors(runtime_lane="strecker_crosstalk", family="08", process_state="heated_matrix")
    crosstalk_prior = dict(crosstalk_prior_rows[0]) if crosstalk_prior_rows else get_strecker_crosstalk_priors()
    cyclodextrin_prior = _computational_prior_entry("bhandari_1998_beta_cd_aldehyde_binding_v1")
    cyclodextrin_observed_signals = {
        "hexanal": _compound_signal(signal_map, ["hexanal"]),
        "nonanal": _compound_signal(signal_map, ["nonanal"]),
        "e_2_nonenal": _compound_signal(signal_map, ["(E)-2-Nonenal", "e-2-nonenal", "2-nonenal", "2 nonenal"]),
        "2_pentylfuran": _compound_signal(signal_map, ["2-pentylfuran", "2 pentylfuran"]),
        "1_octen_3_ol": _compound_signal(signal_map, ["1-octen-3-ol", "1 octen 3 ol"]),
    }
    cyclodextrin_target_compounds = [
        compound
        for compound, value in cyclodextrin_observed_signals.items()
        if value > 0.0
        and compound in (cyclodextrin_prior.get("oav_reduction_factor_at_1pct_w_w_by_compound", {}) or {})
    ]
    cyclodextrin_weighted_oav_reduction_factor = 0.0
    cyclodextrin_weighted_headspace_reduction_fraction = 0.0
    cyclodextrin_reference_active = bool(cyclodextrin_guardrail_active and cyclodextrin_prior and cyclodextrin_target_compounds)
    if cyclodextrin_reference_active:
        total_signal = sum(cyclodextrin_observed_signals[compound] for compound in cyclodextrin_target_compounds)
        if total_signal > 0.0:
            oav_map = cyclodextrin_prior.get("oav_reduction_factor_at_1pct_w_w_by_compound", {}) or {}
            headspace_map = cyclodextrin_prior.get("headspace_reduction_fraction_at_1pct_w_w_by_compound", {}) or {}
            cyclodextrin_weighted_oav_reduction_factor = sum(
                cyclodextrin_observed_signals[compound] * float(oav_map.get(compound, 0.0) or 0.0)
                for compound in cyclodextrin_target_compounds
            ) / total_signal
            cyclodextrin_weighted_headspace_reduction_fraction = sum(
                cyclodextrin_observed_signals[compound] * float(headspace_map.get(compound, 0.0) or 0.0)
                for compound in cyclodextrin_target_compounds
            ) / total_signal
    required_sugars = [_normalize_name(str(item)) for item in crosstalk_prior.get("required_sugars", []) or []]
    sugar_requirement_met = not required_sugars or any(
        any(required in sugar for required in required_sugars)
        for sugar in normalized_sugars
    )
    polyphenol_crosstalk_active = bool(
        polyphenol_active
        and antioxidant_guardrail_active
        and sugar_requirement_met
        and bool(lipid_crosstalk_lane.get("active", False))
    )
    dicarbonyl_trapping_factor = 0.0
    if antioxidant_guardrail_active:
        dicarbonyl_trapping_factor = 0.35
    if polyphenol_active:
        dicarbonyl_trapping_factor = max(dicarbonyl_trapping_factor, 0.55)
    if polyphenol_crosstalk_active:
        dicarbonyl_trapping_factor = max(dicarbonyl_trapping_factor, 0.85)
    amino_group_blocking_factor = 0.0
    if protein_label != "free":
        amino_group_blocking_factor = min(
            1.0,
            0.18 * float(antioxidant_guardrail_active)
            + 0.42 * float(polyphenol_crosstalk_active)
            + 0.30 * min(1.0, lipid_marker_signal_ppb / 120.0),
        )
    safety_reference = get_safety_reference_payload("squeo_2023_pbpi_acrylamide") or {}
    safety_range = get_safety_reference_range("soy_wet_extraction", "squeo_2023_pbpi_acrylamide") or {}
    acrylamide_guardrail_factor = 0.25 * float(calcium_guardrail_active) + 0.10 * float(bool(safety_reference))
    suppression_pressure_score = min(
        1.0,
        0.45 * float(dicarbonyl_trapping_factor)
        + 0.35 * float(amino_group_blocking_factor)
        + 0.20 * float(acrylamide_guardrail_factor),
    )
    if cyclodextrin_reference_active:
        suppression_pressure_score = min(
            1.0,
            suppression_pressure_score + 0.18 * min(1.0, cyclodextrin_weighted_headspace_reduction_fraction / 0.74),
        )
    suppression_pressure_active = bool(
        lipid_marker_signal_ppb > 0.0
        or antioxidant_guardrail_active
        or calcium_guardrail_active
        or polyphenol_active
        or cyclodextrin_guardrail_active
    )
    summary = "No explicit off-note guardrail lane active."
    if suppression_pressure_active:
        summary = (
            "Plant-matrix guardrails are active, so off-note suppression and amino-group blocking risk should constrain optimistic flavor claims."
        )
        if benchmark_rows and safety_reference:
            summary = (
                "Plant-matrix guardrails are active and now carry explicit safety provenance from Squeo 2023 plus polyphenol crosstalk support, so suppression claims stay bounded by reference data instead of cue detection alone."
            )
        if cyclodextrin_reference_active:
            summary = (
                "Plant-matrix guardrails are active and now also expose bounded beta-cyclodextrin sequestration support, so aldehyde suppression is treated as partial headspace reduction rather than complete cleanup."
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
            "cyclodextrin_guardrail_active": bool(cyclodextrin_guardrail_active),
            "dicarbonyl_trapping_factor": float(dicarbonyl_trapping_factor),
            "amino_group_blocking_factor": float(amino_group_blocking_factor),
            "suppression_pressure_score": float(suppression_pressure_score),
            "polyphenol_crosstalk_active": bool(polyphenol_crosstalk_active),
            "sugar_requirement_met": bool(sugar_requirement_met),
            "benchmark_anchor_ids": [str(row.get("id", "unknown")) for row in benchmark_rows],
            "benchmark_anchor_citations": [str(row.get("citation", "unknown")) for row in benchmark_rows],
            "crosstalk_prior_ids": [str(row.get("id", "unknown")) for row in crosstalk_prior_rows] or ([str(crosstalk_prior.get("id", "unknown"))] if crosstalk_prior else []),
            "cyclodextrin_prior_ids": [str(cyclodextrin_prior.get("id", "unknown"))] if cyclodextrin_reference_active else [],
            "cyclodextrin_reference_active": bool(cyclodextrin_reference_active),
            "cyclodextrin_target_compounds": cyclodextrin_target_compounds,
            "cyclodextrin_weighted_oav_reduction_factor": float(cyclodextrin_weighted_oav_reduction_factor),
            "cyclodextrin_weighted_headspace_reduction_fraction": float(cyclodextrin_weighted_headspace_reduction_fraction),
            "cyclodextrin_reference_loading_wt_pct": float(cyclodextrin_prior.get("reference_conditions", {}).get("beta_cyclodextrin_loading_wt_pct", 0.0) or 0.0),
            "safety_reference_ids": [str(safety_reference.get("id", "unknown"))] if safety_reference else [],
            "safety_reference_citations": [str(safety_reference.get("source_citation", "unknown"))] if safety_reference else [],
            "acrylamide_reference_mean_ug_per_kg": float(safety_range.get("mean", safety_range.get("point", 0.0)) or 0.0),
            "acrylamide_reference_range_ug_per_kg": {
                "min": float(safety_range.get("min", 0.0) or 0.0),
                "max": float(safety_range.get("max", 0.0) or 0.0),
            },
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
        "target_score_delta": (
            0.06 * float(matrix_lane.get("matrix_source_support_score", 0.0))
            - 0.12 * float(matrix_lane.get("matrix_uncertainty_factor", 0.0))
        ) if matrix_lane.get("active") else 0.0,
        "off_flavour_risk_delta": (
            0.10 * float(matrix_lane.get("off_note_penalty", 0.0))
            + 0.05 * float(matrix_lane.get("matrix_uncertainty_factor", 0.0))
        ) if matrix_lane.get("active") else 0.0,
    }

    donor_lane = family_lane_summary.get("07", {})
    donor_score = float(donor_lane.get("donor_hierarchy_score", 0.0))
    supports_fast = bool(donor_lane.get("supports_fast_sulfur_routing", False))
    pyrazine_pressure = float(donor_lane.get("pyrazine_pressure_score", 0.0))
    per_lane["07"] = {
        "target_score_delta": (0.14 * donor_score - 0.05 * pyrazine_pressure) if donor_lane.get("active") else 0.0,
        "off_flavour_risk_delta": (0.10 * pyrazine_pressure - 0.06 * donor_score) if donor_lane.get("active") and supports_fast else 0.02 * max(0.0, 0.6 - donor_score),
    }

    guardrail_lane = family_lane_summary.get("08", {})
    per_lane["08"] = {
        "target_score_delta": (
            -0.08 * float(guardrail_lane.get("amino_group_blocking_factor", 0.0))
            - 0.04 * float(guardrail_lane.get("suppression_pressure_score", 0.0))
        ) if guardrail_lane.get("active") else 0.0,
        "off_flavour_risk_delta": (
            0.14 * float(guardrail_lane.get("dicarbonyl_trapping_factor", 0.0))
            + 0.08 * float(guardrail_lane.get("amino_group_blocking_factor", 0.0))
            + 0.06 * float(guardrail_lane.get("suppression_pressure_score", 0.0))
        ) if guardrail_lane.get("active") else 0.0,
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

    lipid_maillard_lane = family_lane_summary.get("11", {})
    hexanal_suppression_fraction = float(lipid_maillard_lane.get("hexanal_suppression_fraction", 0.0))
    per_lane["11"] = {
        "target_score_delta": 0.10 * hexanal_suppression_fraction if lipid_maillard_lane.get("active") else 0.0,
        "maillard_closure_delta": 0.06 * hexanal_suppression_fraction if lipid_maillard_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.18 * hexanal_suppression_fraction if lipid_maillard_lane.get("active") else 0.0,
    }

    protein_damage_lane = family_lane_summary.get("12", {})
    damage_burden = float(protein_damage_lane.get("damage_burden_score", 0.0))
    per_lane["12"] = {
        "target_score_delta": -0.14 * damage_burden if protein_damage_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.08 * damage_burden if protein_damage_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.10 * damage_burden if protein_damage_lane.get("active") else 0.0,
    }

    polyphenol_capping_lane = family_lane_summary.get("13", {})
    cysteine_depletion = float(polyphenol_capping_lane.get("cysteine_depletion_factor", 0.0))
    lysine_depletion = float(polyphenol_capping_lane.get("lysine_depletion_factor", 0.0))
    per_lane["13"] = {
        "target_score_delta": (-0.18 * cysteine_depletion - 0.05 * lysine_depletion) if polyphenol_capping_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.10 * float(polyphenol_capping_lane.get("quinone_budget", 0.0)) if polyphenol_capping_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.0,
    }

    ascorbic_acid_lane = family_lane_summary.get("14", {})
    dicarbonyl_pressure = float(ascorbic_acid_lane.get("dicarbonyl_source_pressure", 0.0))
    pentosidine_load = float(ascorbic_acid_lane.get("pentosidine_load", 0.0))
    per_lane["14"] = {
        "target_score_delta": (-0.08 * dicarbonyl_pressure - 0.10 * pentosidine_load) if ascorbic_acid_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.06 * float(ascorbic_acid_lane.get("ascorbic_acid_loss", 0.0)) if ascorbic_acid_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.08 * dicarbonyl_pressure if ascorbic_acid_lane.get("active") else 0.0,
    }

    phospholipid_lane = family_lane_summary.get("15", {})
    sugar_sink_fraction = float(phospholipid_lane.get("sugar_sink_fraction", 0.0))
    per_lane["15"] = {
        "target_score_delta": -0.10 * sugar_sink_fraction if phospholipid_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.12 * float(phospholipid_lane.get("pe_glycation_fraction", 0.0)) if phospholipid_lane.get("active") else 0.0,
        "off_flavour_risk_delta": -0.02 * sugar_sink_fraction if phospholipid_lane.get("active") else 0.0,
    }

    melanoidin_lane = family_lane_summary.get("16", {})
    trapping_pressure = float(melanoidin_lane.get("matrix_trapping_pressure", 0.0))
    per_lane["16"] = {
        "target_score_delta": -0.12 * trapping_pressure if melanoidin_lane.get("active") else 0.0,
        "maillard_closure_delta": -0.10 * float(melanoidin_lane.get("thiol_scavenging_factor", 0.0)) if melanoidin_lane.get("active") else 0.0,
        "off_flavour_risk_delta": 0.12 * trapping_pressure if melanoidin_lane.get("active") else 0.0,
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
    lipids: Optional[List[str]] = None,
    support_cues: Optional[List[str]] = None,
    interventions: Optional[List[Any]] = None,
    protein_type: Optional[str] = None,
    pH: Optional[float] = None,
    thiamine_availability: Any = None,
    molar_ratios: Optional[Mapping[str, Any]] = None,
    process_state: Optional[str] = None,
    temperature_celsius: Optional[float] = None,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    degree_of_hydrolysis: Optional[float] = None,
) -> Dict[str, Any]:
    sugars = sugars or []
    amino_acids = amino_acids or []
    additives = additives or []
    lipids = lipids or []
    normalized_sugars = [_normalize_name(value) for value in sugars]
    normalized_additives = [_normalize_name(value) for value in additives]
    normalized_lipids = [_normalize_name(value) for value in lipids]
    normalized_amino = [_normalize_name(value) for value in amino_acids]
    normalized_support_cues = [_normalize_name(value) for value in (support_cues or [])]
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
    sulfur_partner_active = any("cysteine" in value for value in normalized_amino + normalized_additives)
    thiamine_fraction_estimate = 0.0
    base_thiamine_fraction = 0.0
    thiamine_mode = "inactive"
    if thiamine_active:
        if pentose_active:
            base_thiamine_fraction = float(thiamine_priors.get("mixed_pentose_fraction", 0.5) or 0.5)
            thiamine_mode = "mixed_thiamine_plus_pentose"
        else:
            base_thiamine_fraction = float(thiamine_priors.get("thiamine_only_fraction", 1.0) or 1.0)
            thiamine_mode = "thiamine_only"

    donor_hierarchy_lane = _build_donor_hierarchy_lane(
        normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        pH=pH,
    )
    fermentation_pretreatment_lane = _build_fermentation_pretreatment_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        pH=pH,
        thiamine_active=thiamine_active,
    )
    sulfur_peptide_support_lane = _build_sulfur_peptide_support_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        protein_label=protein_label,
        degree_of_hydrolysis=degree_of_hydrolysis,
    )
    matrix_scope_lane = _build_matrix_scope_lane(protein_label=protein_label)
    polyphenol_capping_lane = _build_polyphenol_amino_capping_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        protein_label=protein_label,
        pH=pH,
        temperature_celsius=temperature_celsius,
        process_state=process_state,
    )
    ascorbic_acid_lane = _build_ascorbic_acid_maillard_lane(
        normalized_additives=normalized_additives,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        normalized_amino=normalized_amino,
        protein_label=protein_label,
        pH=pH,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
    )

    donor_weight_by_class = {
        "phosphorylated": 1.35,
        "pentose": 1.20,
        "fructose": 1.05,
        "glucose": 0.92,
        "other_reducing_sugar": 0.88,
    }
    # `concentration_unit` passthrough (2026-08-27): payloads may declare the
    # unit of their molar_ratios explicitly; metadata keys are separated from
    # concentrations here so a declared unit cannot crash the float() coercion.
    declared_concentration_unit = _extract_declared_concentration_unit(molar_ratios)
    base_molar_ratios = _numeric_ratio_items(molar_ratios)
    effective_molar_ratios = dict(base_molar_ratios)
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
    cysteine_ratio_key = ratio_key_lookup.get("cysteine")
    if cysteine_ratio_key is not None and polyphenol_capping_lane.get("active"):
        effective_molar_ratios[cysteine_ratio_key] = float(effective_molar_ratios[cysteine_ratio_key]) * max(
            0.05,
            1.0 - float(polyphenol_capping_lane.get("cysteine_depletion_factor", 0.0)),
        )
    lysine_ratio_key = ratio_key_lookup.get("lysine")
    if lysine_ratio_key is not None and polyphenol_capping_lane.get("active"):
        effective_molar_ratios[lysine_ratio_key] = float(effective_molar_ratios[lysine_ratio_key]) * max(
            0.10,
            1.0 - float(polyphenol_capping_lane.get("lysine_depletion_factor", 0.0)),
        )
    phospholipid_amine_lane = _build_phospholipid_amine_sink_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_lipids=normalized_lipids,
        normalized_support_cues=normalized_support_cues,
        pH=pH,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        process_state=process_state,
    )
    sugar_retention_factor = float(phospholipid_amine_lane.get("available_sugar_retention_factor", 1.0) or 1.0)
    for sugar in normalized_sugars:
        ratio_key = ratio_key_lookup.get(sugar)
        if ratio_key is None or not phospholipid_amine_lane.get("active"):
            continue
        donor_class = _classify_donor_family(sugar)
        donor_sink_weight = {
            "phosphorylated": 1.05,
            "pentose": 1.00,
            "fructose": 0.96,
            "glucose": 0.90,
            "other_reducing_sugar": 0.86,
        }.get(donor_class, 0.90)
        sink_adjusted_ratio = float(effective_molar_ratios[ratio_key]) * max(
            0.65,
            1.0 - (1.0 - sugar_retention_factor) * donor_sink_weight,
        )
        base_ratio = float(base_molar_ratios.get(ratio_key, sink_adjusted_ratio))
        effective_molar_ratios[ratio_key] = min(
            sink_adjusted_ratio,
            base_ratio * max(0.65, sugar_retention_factor),
        )
    melanoidin_polymerization_lane = _build_melanoidin_polymerization_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        protein_label=protein_label,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
    )

    effective_pH = None if pH is None else float(pH)
    if fermentation_pretreatment_lane.get("pretreatment_pH_shift_active") and effective_pH is not None:
        effective_pH = max(4.8, effective_pH - 0.35)

    thiamine_calibration: Dict[str, Any] = {}
    nucleotide_priors = get_nucleotide_priors()
    nucleotide_calibration = _build_nucleotide_calibration_context(
        nucleotide_priors=nucleotide_priors,
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        pH=effective_pH,
        process_state=process_state,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
    )
    if thiamine_active:
        thiamine_calibration = _build_thiamine_calibration_context(
            thiamine_priors=thiamine_priors,
            baseline_fraction=base_thiamine_fraction,
            thiamine_mode=thiamine_mode,
            molar_ratios=effective_molar_ratios,
            pH=effective_pH,
            process_state=process_state,
            temperature_celsius=temperature_celsius,
            water_activity=water_activity,
            sulfur_partner_active=sulfur_partner_active,
        )
        thiamine_fraction_estimate = float(thiamine_calibration.get("effective_fraction", base_thiamine_fraction))

    thiamine_support_lane = _build_thiamine_support_lane(
        thiamine_active=thiamine_active,
        thiamine_fraction_estimate=thiamine_fraction_estimate,
        thiamine_mode=thiamine_mode,
        thiamine_calibration=thiamine_calibration,
    )
    nucleotide_support_lane = _build_nucleotide_support_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        nucleotide_calibration=nucleotide_calibration,
    )

    thiamine_ratio_key = ratio_key_lookup.get("thiamine")
    if thiamine_ratio_key is not None and thiamine_support_lane.get("active"):
        baseline_fraction = float(
            thiamine_calibration.get("baseline_fraction", thiamine_fraction_estimate)
            or thiamine_fraction_estimate
            or 1.0
        )
        effective_fraction = float(
            thiamine_calibration.get("effective_fraction", thiamine_fraction_estimate)
            or thiamine_fraction_estimate
        )
        thiamine_reactivity_factor = max(
            0.05,
            min(1.35, effective_fraction / max(baseline_fraction, 1.0e-6)),
        )
        effective_molar_ratios[thiamine_ratio_key] = (
            float(effective_molar_ratios[thiamine_ratio_key]) * thiamine_reactivity_factor
        )

    amino_accessibility_factor = float(sulfur_peptide_support_lane.get("peptide_accessibility_factor", 1.0) or 1.0)
    cysteine_accessibility_factor = float(sulfur_peptide_support_lane.get("free_cysteine_equivalent_factor", 1.0) or 1.0)
    melanoidin_thiol_retention_factor = 1.0
    melanoidin_supported_context = bool(
        protein_label != "free"
        or thiamine_support_lane.get("active")
        or sulfur_peptide_support_lane.get("active")
    )
    if melanoidin_polymerization_lane.get("active") and melanoidin_supported_context:
        melanoidin_thiol_retention_factor = max(
            0.18,
            1.0 - 0.65 * float(melanoidin_polymerization_lane.get("thiol_scavenging_factor", 0.0)),
        )

    if cysteine_ratio_key is not None:
        effective_molar_ratios[cysteine_ratio_key] = float(effective_molar_ratios[cysteine_ratio_key]) * cysteine_accessibility_factor * melanoidin_thiol_retention_factor
    if lysine_ratio_key is not None and degree_of_hydrolysis is not None:
        effective_molar_ratios[lysine_ratio_key] = float(effective_molar_ratios[lysine_ratio_key]) * amino_accessibility_factor

    added_precursors: List[str] = []
    added_precursor_ratios: Dict[str, float] = {}
    support_pool = normalized_sugars + normalized_additives + normalized_amino
    if thiamine_active and not any("thiamine" in value or "vitamin b1" in value for value in support_pool):
        added_precursors.append("thiamine")
        added_precursor_ratios["thiamine"] = float(max(0.01, 0.25 * max(thiamine_fraction_estimate, 0.0)))
    if (
        not pentose_active
        and not bool(nucleotide_calibration.get("explicit_ribose_active", False))
        and float(nucleotide_calibration.get("ribose_delivery_factor", 0.0)) >= float(nucleotide_priors.get("ribose_addition_threshold", 0.08) or 0.08)
    ):
        added_precursors.append("ribose")
        added_precursor_ratios["ribose"] = float(
            max(
                0.02,
                float(nucleotide_priors.get("ribose_addition_ratio_scale", 0.35) or 0.35)
                * float(nucleotide_calibration.get("ribose_delivery_factor", 0.0)),
            )
        )

    return {
        "effective_molar_ratios": effective_molar_ratios,
        # Declared unit of the mapping above; empty string means "undeclared, and
        # therefore assumed to be MOLAR_RATIO_CANONICAL_UNIT (mM)".
        "concentration_unit": declared_concentration_unit,
        "effective_pH": effective_pH,
        "pretreatment_active": bool(fermentation_pretreatment_lane.get("active", False)),
        "pretreatment_interventions": [
            token
            for token in normalized_interventions
            if token in {"yeast fermentation", "protease hydrolysis", "yeast_fermentation", "protease_hydrolysis"}
        ],
        "support_cues": normalized_support_cues,
        "donor_pool_factors": donor_pool_factors,
        "donor_limited": bool(donor_hierarchy_lane.get("active", False) and dominant_donor_class in {"glucose", "other_reducing_sugar", "none"}),
        "dominant_donor_class": dominant_donor_class,
        "supports_fast_sulfur_routing": bool(donor_hierarchy_lane.get("supports_fast_sulfur_routing", False)),
        "added_precursors": added_precursors,
        "added_precursor_ratios": added_precursor_ratios,
        "thiamine_metadata": thiamine_metadata,
        "thiamine_fraction_baseline": float(base_thiamine_fraction),
        "thiamine_fraction_estimate": float(thiamine_fraction_estimate),
        "thiamine_mode": thiamine_mode,
        "thiamine_calibration": thiamine_calibration,
        "nucleotide_calibration": nucleotide_calibration,
        "degree_of_hydrolysis": None if degree_of_hydrolysis is None else float(max(0.0, min(1.0, degree_of_hydrolysis))),
        "family_lanes": {
            "03": thiamine_support_lane,
            "04": nucleotide_support_lane,
            "05": sulfur_peptide_support_lane,
            "06": matrix_scope_lane,
            "07": donor_hierarchy_lane,
            "10": fermentation_pretreatment_lane,
            "13": polyphenol_capping_lane,
            "14": ascorbic_acid_lane,
            "15": phospholipid_amine_lane,
            "16": melanoidin_polymerization_lane,
        },
        "summary": {
            "donor_identity_first_class": bool(donor_hierarchy_lane.get("active", False)),
            "fermentation_pretreatment_active": bool(fermentation_pretreatment_lane.get("active", False)),
            "nucleotide_support_active": bool(nucleotide_support_lane.get("active", False)),
            "sulfur_peptide_support_active": bool(sulfur_peptide_support_lane.get("active", False)),
            "matrix_scope_active": bool(matrix_scope_lane.get("active", False)),
            "polyphenol_precursor_sink_active": bool(polyphenol_capping_lane.get("active", False)),
            "ascorbic_dicarbonyl_source_active": bool(ascorbic_acid_lane.get("active", False)),
            "phospholipid_sugar_sink_active": bool(phospholipid_amine_lane.get("active", False)),
            "melanoidin_trapping_active": bool(melanoidin_polymerization_lane.get("active", False)),
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
    process_state: Optional[str] = None,
    temperature_celsius: Optional[float] = None,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    degree_of_hydrolysis: Optional[float] = None,
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
            lipids=lipids,
            interventions=interventions,
            protein_type=protein_type,
            pH=pH,
            thiamine_availability=thiamine_availability,
            molar_ratios=molar_ratios,
            process_state=process_state,
            temperature_celsius=temperature_celsius,
            time_minutes=time_minutes,
            water_activity=water_activity,
            degree_of_hydrolysis=degree_of_hydrolysis,
        )
    )
    normalized_support_cues = [_normalize_name(value) for value in (upstream_contract.get("support_cues", []) or [])]

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
    thiamine_calibration = dict(upstream_contract.get("thiamine_calibration", {}))
    thiamine_fraction_estimate = float(upstream_contract.get("thiamine_fraction_estimate", 0.0) or 0.0)
    thiamine_mode = str(upstream_contract.get("thiamine_mode", "inactive"))
    if thiamine_active and thiamine_fraction_estimate <= 0.0:
        if pentose_active:
            thiamine_fraction_estimate = float(thiamine_priors.get("mixed_pentose_fraction", 0.5) or 0.5)
            thiamine_mode = "mixed_thiamine_plus_pentose"
        else:
            thiamine_fraction_estimate = float(thiamine_priors.get("thiamine_only_fraction", 1.0) or 1.0)
            thiamine_mode = "thiamine_only"
    if thiamine_active and not thiamine_calibration:
        thiamine_calibration = _build_thiamine_calibration_context(
            thiamine_priors=thiamine_priors,
            baseline_fraction=float(upstream_contract.get("thiamine_fraction_baseline", thiamine_fraction_estimate) or thiamine_fraction_estimate),
            thiamine_mode=thiamine_mode,
            pH=pH,
            process_state=process_state,
            temperature_celsius=temperature_celsius,
            water_activity=water_activity,
            sulfur_partner_active=any("cysteine" in value for value in normalized_amino + normalized_additives),
        )
        thiamine_fraction_estimate = float(thiamine_calibration.get("effective_fraction", thiamine_fraction_estimate))

    nucleotide_priors = get_nucleotide_priors()
    nucleotide_calibration = dict(upstream_contract.get("nucleotide_calibration", {}))
    if (not nucleotide_calibration) and any(
        _any_token_matches(value, ["imp", "gmp", "inosinate", "guanylate", "yeast extract", "ribose", "ribose 5 phosphate", "ribose-5-phosphate", "r5p"])
        for value in normalized_sugars + normalized_additives + normalized_interventions
    ):
        nucleotide_calibration = _build_nucleotide_calibration_context(
            nucleotide_priors=nucleotide_priors,
            normalized_sugars=normalized_sugars,
            normalized_additives=normalized_additives,
            normalized_interventions=normalized_interventions,
            pH=float(upstream_contract.get("effective_pH", pH)) if upstream_contract.get("effective_pH", pH) is not None else None,
            process_state=process_state,
            temperature_celsius=temperature_celsius,
            time_minutes=time_minutes,
        )

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
    lipid_maillard_competition_lane = _build_lipid_maillard_competition_lane(
        signal_map=signal_map,
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        lipids=lipids,
        protein_type=protein_label,
        lipid_crosstalk_lane=lipid_crosstalk_lane,
    )
    donor_hierarchy_lane = _build_donor_hierarchy_lane(
        normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        pH=pH,
    )
    fermentation_pretreatment_lane = _build_fermentation_pretreatment_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        pH=pH,
        thiamine_active=thiamine_active,
    )
    off_note_guardrail_lane = _build_off_note_guardrail_lane(
        signal_map=signal_map,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        normalized_sugars=normalized_sugars,
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
        thiamine_calibration=thiamine_calibration,
    )
    nucleotide_support_lane = _build_nucleotide_support_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_interventions=normalized_interventions,
        nucleotide_calibration=nucleotide_calibration,
    )
    sulfur_peptide_support_lane = _build_sulfur_peptide_support_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        protein_label=protein_label,
    )
    matrix_scope_lane = _build_matrix_scope_lane(protein_label=protein_label)
    caramelization_lane = _build_caramelization_lane(
        signal_map=signal_map,
        furanone_expected=furanone_expected,
        furanone_observed=furanone_observed,
        pH=pH,
    )
    protein_damage_lane = _build_protein_damage_markers_lane(
        normalized_sugars=normalized_sugars,
        normalized_amino=normalized_amino,
        normalized_additives=normalized_additives,
        protein_type=protein_label,
        process_state=process_state,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        effective_molar_ratios=upstream_contract.get("effective_molar_ratios", {}) or {},
        pH=pH,
        concentration_unit=str(upstream_contract.get("concentration_unit", "") or "") or None,
    )
    polyphenol_capping_lane = dict(upstream_contract.get("family_lanes", {}).get("13", {})) or _build_polyphenol_amino_capping_lane(
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        protein_label=protein_label,
        pH=pH,
        temperature_celsius=temperature_celsius,
        process_state=process_state,
    )
    ascorbic_acid_lane = dict(upstream_contract.get("family_lanes", {}).get("14", {})) or _build_ascorbic_acid_maillard_lane(
        normalized_additives=normalized_additives,
        normalized_support_cues=normalized_support_cues,
        normalized_interventions=normalized_interventions,
        normalized_amino=normalized_amino,
        protein_label=protein_label,
        pH=pH,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
    )
    phospholipid_amine_lane = dict(upstream_contract.get("family_lanes", {}).get("15", {})) or _build_phospholipid_amine_sink_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_lipids=[_normalize_name(value) for value in lipids],
        normalized_support_cues=normalized_support_cues,
        pH=pH,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        process_state=process_state,
    )
    melanoidin_polymerization_lane = dict(upstream_contract.get("family_lanes", {}).get("16", {})) or _build_melanoidin_polymerization_lane(
        normalized_sugars=normalized_sugars,
        normalized_additives=normalized_additives,
        normalized_amino=normalized_amino,
        protein_label=protein_label,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
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
            lipid_maillard_competition_lane,
            protein_damage_lane,
            polyphenol_capping_lane,
            ascorbic_acid_lane,
            phospholipid_amine_lane,
            melanoidin_polymerization_lane,
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
        protein_damage_lane=protein_damage_lane,
        polyphenol_capping_lane=polyphenol_capping_lane,
        ascorbic_acid_lane=ascorbic_acid_lane,
        phospholipid_amine_lane=phospholipid_amine_lane,
        melanoidin_polymerization_lane=melanoidin_polymerization_lane,
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
        "kokumi_support_active": bool(sulfur_peptide_support_lane.get("kokumi_support_active", False)),
        "kokumi_support_signal": float(sulfur_peptide_support_lane.get("kokumi_signal_score", 0.0)),
        "kokumi_signal_mode": str(sulfur_peptide_support_lane.get("kokumi_signal_mode", "inactive")),
        "lincoln_crosstalk_prior": lincoln_crosstalk_prior,
        "family_upstream_contract": upstream_contract,
        "family_lane_summary": family_lane_summary,
        "family_lane_adjustments": family_lane_adjustments,
        "family_prior_bundle": family_prior_bundle,
        "family_state_markers": family_state_markers,
        "active_family_lanes": active_family_lanes,
        "active_family_lane_names": active_family_lane_names,
    }