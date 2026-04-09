from __future__ import annotations

"""
Correction factors for amino acid reactivity in protein matrices
vs. free amino acid model systems.

These are empirical correction factors derived from literature comparing
model system predictions to protein matrix experiments. They should be
recalibrated as you generate your own experimental data.

Current legume-matrix anchoring inside the repo:
- pea: Asen 2022 DSC/DTNB + Li 2025 Ellman SH envelope, with older repo
    literature synthesis retained as a conservative interpolation layer
- soy: repo literature synthesis around soy glycinin/beta-conglycinin burial,
  sulfur limitation, extrusion-driven accessibility changes, and soy
  protein-polysaccharide conjugate trapping documented in
    docs/research/archives/Maillard_Plant_based.md and
    docs/research/archives/Elicit - Maillard Pathways in Plant-Based Cooking - Report.md
"""

import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from src.matrix_prior_registry import (
    get_accessibility_window_entry,
    get_denaturation_heuristic_entry,
    get_matrix_correction_entry,
    get_volatile_class_profile_entry,
    summarize_matrix_prior_bundle,
)
from src.literature_runtime import describe_retention_runtime


ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"


def _load_json_payload(file_name: str) -> dict:
    payload_path = DATA_LIT_DIR / file_name
    with open(payload_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


PROCESS_STATE_CALIBRATION_PAYLOAD = _load_json_payload("process_state_calibrations.json")
PROTEIN_SOURCE_REGISTRY_PAYLOAD = _load_json_payload("protein_source_registry.json")

@dataclass(frozen=True)
class ProteinSourceProfile:
    source_id: str
    meaty_potential_multiplier: float
    hydrolysate_observability_bias: float
    off_note_penalty: float
    lox_activity_flag: bool
    methoxypyrazine_ceiling: float
    aa_composition: dict[str, float]

def _build_protein_source_profiles() -> dict[str, ProteinSourceProfile]:
    profiles = {}
    for entry in PROTEIN_SOURCE_REGISTRY_PAYLOAD.get("sources", []):
        profiles[entry["source_id"]] = ProteinSourceProfile(
            source_id=entry["source_id"],
            meaty_potential_multiplier=float(entry.get("meaty_potential_multiplier", 1.0)),
            hydrolysate_observability_bias=float(entry.get("hydrolysate_observability_bias", 1.0)),
            off_note_penalty=float(entry.get("off_note_penalty", 0.0)),
            lox_activity_flag=bool(entry.get("lox_activity_flag", False)),
            methoxypyrazine_ceiling=float(entry.get("methoxypyrazine_ceiling", 1.0)),
            aa_composition={k: float(v) for k, v in entry.get("aa_composition", {}).items()},
        )
    return profiles

PROTEIN_SOURCE_PROFILES = _build_protein_source_profiles()

def get_protein_source_profile(source_id: str) -> ProteinSourceProfile | None:
    if not source_id:
        return None
    return PROTEIN_SOURCE_PROFILES.get(str(source_id).strip().lower())



def get_process_state_calibration_payload(protein_type: ProteinType | str) -> list[dict]:
    normalized = protein_type.value if isinstance(protein_type, ProteinType) else str(protein_type)
    return [
        entry
        for entry in PROCESS_STATE_CALIBRATION_PAYLOAD.get("entries", [])
        if str(entry.get("protein_type", "")) == normalized
    ]


class ProteinType(Enum):
    FREE_AMINO_ACID = "free"
    PEA_CONCENTRATE = "pea_conc"    # ~60% protein, fibrous matrix
    PEA_ISOLATE = "pea_iso"         # ~85% protein
    SOY_CONCENTRATE = "soy_conc"
    SOY_ISOLATE = "soy_iso"
    MYCOPROTEIN = "myco"

@dataclass
class MatrixCorrection:
    protein_type: ProteinType
    lysine_accessibility: float     # fraction of total lysine reactive
    cysteine_accessibility: float   # usually lower due to disulfide burial
    lysine_accessibility_native: float
    lysine_accessibility_denatured: float
    cysteine_accessibility_native: float
    cysteine_accessibility_denatured: float
    volatile_retention: float       # fraction escaping matrix (rest is bound)
    volatile_retention_native: float
    volatile_retention_denatured: float
    source: str


@dataclass(frozen=True)
class EffectiveMatrixCorrection:
    protein_type: ProteinType
    lysine_accessibility: float
    cysteine_accessibility: float
    volatile_retention: float
    source: str


@dataclass(frozen=True)
class AccessibilityLiteratureWindow:
    protein_type: ProteinType
    lysine_min: float
    lysine_max: float
    cysteine_min: float
    cysteine_max: float
    source: str


@dataclass(frozen=True)
class DenaturationHeuristic:
    protein_type: ProteinType
    midpoint_celsius: float
    width_celsius: float
    time_gain_celsius: float
    time_reference_minutes: float
    acidic_ph_gain_celsius: float
    reference_ph: float
    source: str


@dataclass(frozen=True)
class VolatileClassRetentionProfile:
    protein_type: ProteinType
    native_factors: dict[str, float]
    denatured_factors: dict[str, float]
    source: str


@dataclass(frozen=True)
class AccessibilityState:
    """
    Canonical accessibility state for a protein matrix run.

    Maps a denaturation float (0–1) plus protein_type to a named profile
    so that the pipeline and reports can reason about *why* accessibility
    constrains predictions, not just *how much*.

    Profiles:
      protein_embedded  — 0.0–0.2: most reactive sites buried; chemistry
                          limited by physical accessibility, not kinetics.
      peptide_bound     — 0.2–0.6: partial denaturation; mixed burial/exposure.
      partially_opened  — 0.6–0.9: processing has exposed most reactive sites.
      free_like         — 0.9–1.0: fully denatured; behaves like free precursor.
    """
    protein_type: str
    denaturation_state: float
    profile: str               # one of the four canonical labels above
    dominant_source: str       # 'denaturation_state_arg' | 'estimated_from_conditions'
    accessibility_warning: bool  # True when accessibility dominates uncertainty


def classify_accessibility_state(
    protein_type: str,
    denaturation_state: float,
    *,
    dominant_source: str = "denaturation_state_arg",
) -> "AccessibilityState":
    """Return the canonical AccessibilityState for a given protein_type and denaturation float."""
    d = max(0.0, min(1.0, float(denaturation_state)))
    pt_lower = str(protein_type).strip().lower()

    # Free systems are always fully accessible
    if pt_lower in {"free", "free_amino_acid", ""}:
        return AccessibilityState(
            protein_type=protein_type,
            denaturation_state=1.0,
            profile="free_like",
            dominant_source="protein_type_is_free",
            accessibility_warning=False,
        )

    if d < 0.2:
        profile = "protein_embedded"
        warning = True
    elif d < 0.6:
        profile = "peptide_bound"
        warning = True
    elif d < 0.9:
        profile = "partially_opened"
        warning = False
    else:
        profile = "free_like"
        warning = False

    return AccessibilityState(
        protein_type=protein_type,
        denaturation_state=d,
        profile=profile,
        dominant_source=dominant_source,
        accessibility_warning=warning,
    )



_LYSINE_IDENTIFIERS = {
    "lysine",
    "l-lysine",
    "nccccc(n)c(=o)o",
}

_CYSTEINE_IDENTIFIERS = {
    "cysteine",
    "l-cysteine",
    "nc(cs)c(=o)o",
}


def _build_accessibility_literature_windows() -> dict[ProteinType, AccessibilityLiteratureWindow]:
    windows: dict[ProteinType, AccessibilityLiteratureWindow] = {}
    for protein_type in ProteinType:
        entry = get_accessibility_window_entry(protein_type.value)
        if entry is None:
            continue
        windows[protein_type] = AccessibilityLiteratureWindow(
            protein_type=protein_type,
            lysine_min=float(entry["lysine_min"]),
            lysine_max=float(entry["lysine_max"]),
            cysteine_min=float(entry["cysteine_min"]),
            cysteine_max=float(entry["cysteine_max"]),
            source=str(entry["source"]),
        )
    return windows


def _build_denaturation_heuristics() -> dict[ProteinType, DenaturationHeuristic]:
    heuristics: dict[ProteinType, DenaturationHeuristic] = {}
    for protein_type in ProteinType:
        if protein_type == ProteinType.FREE_AMINO_ACID:
            continue
        entry = get_denaturation_heuristic_entry(protein_type.value)
        if entry is None:
            continue
        heuristics[protein_type] = DenaturationHeuristic(
            protein_type=protein_type,
            midpoint_celsius=float(entry["midpoint_celsius"]),
            width_celsius=float(entry["width_celsius"]),
            time_gain_celsius=float(entry["time_gain_celsius"]),
            time_reference_minutes=float(entry["time_reference_minutes"]),
            acidic_ph_gain_celsius=float(entry["acidic_ph_gain_celsius"]),
            reference_ph=float(entry["reference_ph"]),
            source=str(entry["source"]),
        )
    return heuristics


def _build_volatile_class_retention_profiles() -> dict[ProteinType, VolatileClassRetentionProfile]:
    profiles: dict[ProteinType, VolatileClassRetentionProfile] = {}
    for protein_type in ProteinType:
        if protein_type == ProteinType.FREE_AMINO_ACID:
            continue
        entry = get_volatile_class_profile_entry(protein_type.value)
        if entry is None:
            continue
        profiles[protein_type] = VolatileClassRetentionProfile(
            protein_type=protein_type,
            native_factors={key: float(value) for key, value in entry.get("native_factors", {}).items()},
            denatured_factors={key: float(value) for key, value in entry.get("denatured_factors", {}).items()},
            source=str(entry["source"]),
        )
    return profiles


def _build_matrix_corrections() -> dict[ProteinType, MatrixCorrection]:
    corrections: dict[ProteinType, MatrixCorrection] = {}
    for protein_type in ProteinType:
        entry = get_matrix_correction_entry(protein_type.value)
        if entry is None:
            continue
        corrections[protein_type] = MatrixCorrection(
            protein_type=protein_type,
            lysine_accessibility=float(entry["lysine_accessibility_mid"]),
            cysteine_accessibility=float(entry["cysteine_accessibility_mid"]),
            lysine_accessibility_native=float(entry["lysine_accessibility_native"]),
            lysine_accessibility_denatured=float(entry["lysine_accessibility_denatured"]),
            cysteine_accessibility_native=float(entry["cysteine_accessibility_native"]),
            cysteine_accessibility_denatured=float(entry["cysteine_accessibility_denatured"]),
            volatile_retention=float(entry["volatile_retention_mid"]),
            volatile_retention_native=float(entry["volatile_retention_native"]),
            volatile_retention_denatured=float(entry["volatile_retention_denatured"]),
            source=str(entry["source"]),
        )
    return corrections


ACCESSIBILITY_LITERATURE_WINDOWS = _build_accessibility_literature_windows()
DENATURATION_HEURISTICS = _build_denaturation_heuristics()
VOLATILE_CLASS_RETENTION_PROFILES = _build_volatile_class_retention_profiles()
MATRIX_CORRECTIONS = _build_matrix_corrections()


def clamp_denaturation_state(denaturation_state: float) -> float:
    return max(0.0, min(1.0, float(denaturation_state)))


def _coerce_protein_type(protein_type: ProteinType | str) -> ProteinType:
    if isinstance(protein_type, ProteinType):
        return protein_type
    return ProteinType(protein_type)


def estimate_denaturation_state(
    protein_type: ProteinType | str,
    temperature_celsius: float,
    time_minutes: Optional[float] = None,
    pH: Optional[float] = None,
) -> float:
    p_type = _coerce_protein_type(protein_type)
    if p_type == ProteinType.FREE_AMINO_ACID:
        return 1.0

    heuristic = DENATURATION_HEURISTICS.get(p_type)
    if heuristic is None:
        return 0.5

    duration = max(0.0, float(time_minutes if time_minutes is not None else 0.0))
    acidity_shift = 0.0
    if pH is not None:
        acidity_shift = heuristic.acidic_ph_gain_celsius * max(0.0, heuristic.reference_ph - float(pH))

    time_shift = heuristic.time_gain_celsius * math.log1p(duration / heuristic.time_reference_minutes)
    effective_midpoint = heuristic.midpoint_celsius - time_shift - acidity_shift
    exponent = max(-60.0, min(60.0, -(float(temperature_celsius) - effective_midpoint) / heuristic.width_celsius))
    return clamp_denaturation_state(1.0 / (1.0 + math.exp(exponent)))


def resolve_effective_denaturation_state(
    protein_type: ProteinType | str,
    temperature_celsius: float,
    time_minutes: Optional[float] = None,
    pH: Optional[float] = None,
    explicit_denaturation_state: Optional[float] = None,
) -> float:
    if explicit_denaturation_state is not None:
        return clamp_denaturation_state(explicit_denaturation_state)
    return estimate_denaturation_state(
        protein_type,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        pH=pH,
    )


def get_accessibility_literature_window(
    protein_type: ProteinType,
) -> AccessibilityLiteratureWindow | None:
    return ACCESSIBILITY_LITERATURE_WINDOWS.get(protein_type)


def classify_volatile_matrix_family(name: str, smiles: Optional[str] = None) -> str:
    normalized = name.strip().lower()
    if any(token in normalized for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
        return "sulfur"
    if "pyrazine" in normalized:
        return "pyrazine"
    if any(token in normalized for token in ["furan", "furfural"]):
        return "furan"
    if any(token in normalized for token in ["ol", "alcohol"]):
        return "alcohol"
    if any(token in normalized for token in ["anal", "enal", "aldehyde"]):
        return "aldehyde"
    if smiles:
        smi = smiles.lower()
        if "s" in smi:
            return "sulfur"
        if "n" in smi and "c1" in smi:
            return "pyrazine"
        if "o" in smi and "c1" in smi:
            return "furan"
    return "other"


def get_volatile_class_retention_factor(
    name: str,
    protein_type: ProteinType | str,
    denaturation_state: float,
    smiles: Optional[str] = None,
) -> float:
    p_type = _coerce_protein_type(protein_type)
    if p_type == ProteinType.FREE_AMINO_ACID:
        return 1.0
    profile = VOLATILE_CLASS_RETENTION_PROFILES.get(p_type)
    if profile is None:
        return 1.0
    family = classify_volatile_matrix_family(name, smiles=smiles)
    state = clamp_denaturation_state(denaturation_state)
    native = profile.native_factors.get(family, profile.native_factors.get("other", 1.0))
    denatured = profile.denatured_factors.get(family, profile.denatured_factors.get("other", 1.0))
    return native + state * (denatured - native)


def resolve_compound_matrix_retention(
    name: str,
    protein_type: ProteinType | str,
    denaturation_state: float = 0.5,
    smiles: Optional[str] = None,
    temperature_celsius: Optional[float] = None,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    process_state: Optional[str] = None,
    protein_source: Optional[str] = None,
) -> float:
    p_type = _coerce_protein_type(protein_type)
    base = resolve_matrix_correction(p_type, denaturation_state).volatile_retention
    class_factor = get_volatile_class_retention_factor(
        name,
        protein_type=p_type,
        denaturation_state=denaturation_state,
        smiles=smiles,
    )
    runtime = describe_retention_runtime(
        name,
        protein_type=p_type.value,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
    )
    dynamic_factor = float(runtime.get("dynamic_retention_factor", 1.0))
    
    source_factor = 1.0
    profile = get_protein_source_profile(protein_source)
    if profile:
        source_factor = profile.meaty_potential_multiplier
        if "methoxypyrazine" in name.lower():
            source_factor = min(source_factor, profile.methoxypyrazine_ceiling)
        compound_family = classify_volatile_matrix_family(name, smiles=smiles)
        if compound_family in ["aldehyde", "alcohol"] and profile.lox_activity_flag:
            source_factor *= (1.0 + profile.off_note_penalty)

    return max(0.01, min(1.0, base * class_factor * dynamic_factor * source_factor))


def describe_compound_matrix_retention(
    name: str,
    protein_type: ProteinType | str,
    denaturation_state: float = 0.5,
    smiles: Optional[str] = None,
    temperature_celsius: Optional[float] = None,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    process_state: Optional[str] = None,
    protein_source: Optional[str] = None,
) -> dict[str, object]:
    p_type = _coerce_protein_type(protein_type)
    base = resolve_matrix_correction(p_type, denaturation_state).volatile_retention
    class_factor = get_volatile_class_retention_factor(
        name,
        protein_type=p_type,
        denaturation_state=denaturation_state,
        smiles=smiles,
    )
    runtime = describe_retention_runtime(
        name,
        protein_type=p_type.value,
        temperature_celsius=temperature_celsius,
        time_minutes=time_minutes,
        water_activity=water_activity,
        process_state=process_state,
    )
    dynamic_factor = float(runtime.get("dynamic_retention_factor", 1.0))
    
    source_factor = 1.0
    profile = get_protein_source_profile(protein_source)
    if profile:
        source_factor = profile.meaty_potential_multiplier
        if "methoxypyrazine" in name.lower():
            source_factor = min(source_factor, profile.methoxypyrazine_ceiling)
        compound_family = classify_volatile_matrix_family(name, smiles=smiles)
        if compound_family in ["aldehyde", "alcohol"] and profile.lox_activity_flag:
            source_factor *= (1.0 + profile.off_note_penalty)

    matrix_factor = max(0.01, min(1.0, base * class_factor * dynamic_factor * source_factor))
    return {
        "base_matrix_factor": float(base),
        "class_matrix_factor": float(class_factor),
        "dynamic_retention_factor": float(dynamic_factor),
        "source_heuristic_factor": float(source_factor),
        "matrix_factor": float(matrix_factor),
        **runtime,
    }


def build_matrix_explainability(
    *,
    protein_type: ProteinType | str,
    effective_denaturation_state: float,
    temperature_celsius: float,
    time_minutes: Optional[float],
    pH: Optional[float],
    dominant_source: Optional[str] = None,
) -> dict[str, object]:
    p_type = _coerce_protein_type(protein_type)
    effective = resolve_matrix_correction(p_type, effective_denaturation_state)
    prior_summary = summarize_matrix_prior_bundle(p_type.value)
    state_source = dominant_source or "denaturation_state_arg"
    acc_state = classify_accessibility_state(
        p_type.value,
        effective_denaturation_state,
        dominant_source=state_source,
    )
    if p_type == ProteinType.FREE_AMINO_ACID:
        denaturation_source = "explicit/free"
    elif state_source == "estimated_from_conditions":
        denaturation_source = DENATURATION_HEURISTICS.get(p_type).source if p_type in DENATURATION_HEURISTICS else "estimated_from_conditions"
    else:
        denaturation_source = "explicit_override"
    uncertainty_postures = sorted(
        {
            str(item.get("uncertainty_posture", "unknown"))
            for item in prior_summary.values()
            if isinstance(item, dict) and item.get("uncertainty_posture")
        }
    )
    process_state_applicability = sorted(
        {
            str(state)
            for item in prior_summary.values()
            if isinstance(item, dict)
            for state in item.get("process_state_applicability", [])
        }
    )
    return {
        "protein_type": p_type.value,
        "effective_denaturation_state": float(effective_denaturation_state),
        "temperature_celsius": float(temperature_celsius),
        "time_minutes": None if time_minutes is None else float(time_minutes),
        "pH": None if pH is None else float(pH),
        "lysine_accessibility": float(effective.lysine_accessibility),
        "cysteine_accessibility": float(effective.cysteine_accessibility),
        "bulk_volatile_retention": float(effective.volatile_retention),
        "literature_window": (
            {
                "lysine_min": window.lysine_min,
                "lysine_max": window.lysine_max,
                "cysteine_min": window.cysteine_min,
                "cysteine_max": window.cysteine_max,
                "source": window.source,
            }
            if (window := get_accessibility_literature_window(p_type)) is not None
            else None
        ),
        "denaturation_source": denaturation_source,
        "prior_summary": prior_summary,
        "matrix_prior_uncertainty_postures": uncertainty_postures,
        "matrix_prior_process_state_applicability": process_state_applicability,
        # P1: canonical accessibility state
        "accessibility_profile": acc_state.profile,
        "accessibility_warning": acc_state.accessibility_warning,
        "accessibility_dominant_source": acc_state.dominant_source,
    }


def resolve_matrix_correction(
    protein_type: ProteinType,
    denaturation_state: float = 0.5,
) -> EffectiveMatrixCorrection:
    corr = MATRIX_CORRECTIONS.get(
        protein_type, MATRIX_CORRECTIONS[ProteinType.FREE_AMINO_ACID]
    )
    state = clamp_denaturation_state(denaturation_state)
    lys_factor = corr.lysine_accessibility_native + state * (
        corr.lysine_accessibility_denatured - corr.lysine_accessibility_native
    )
    cys_factor = corr.cysteine_accessibility_native + state * (
        corr.cysteine_accessibility_denatured - corr.cysteine_accessibility_native
    )
    volatile_retention = corr.volatile_retention_native + state * (
        corr.volatile_retention_denatured - corr.volatile_retention_native
    )
    return EffectiveMatrixCorrection(
        protein_type=protein_type,
        lysine_accessibility=lys_factor,
        cysteine_accessibility=cys_factor,
        volatile_retention=volatile_retention,
        source=corr.source,
    )


def _classify_reactive_aa_key(key: str) -> str | None:
    normalized = key.strip().lower()
    if normalized in _LYSINE_IDENTIFIERS:
        return "lysine"
    if normalized in _CYSTEINE_IDENTIFIERS:
        return "cysteine"
    return None

def apply_matrix_correction(
    predicted_concentrations: dict[str, float],
    reactive_amino_acids: dict[str, float],
    protein_type: ProteinType,
    denaturation_state: float = 0.5,  # 0=native, 1=fully denatured
    protein_source: Optional[str] = None,
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Scale predicted volatile concentrations and reactive AA concentrations
    by matrix accessibility factors and specific protein source heuristics.
    """
    corr = resolve_matrix_correction(
        protein_type,
        denaturation_state=denaturation_state,
    )
    
    profile = get_protein_source_profile(protein_source)

    corrected_aa = {}
    for aa, conc in reactive_amino_acids.items():
        aa_class = _classify_reactive_aa_key(aa)
        
        aa_scaler = 1.0
        if profile:
            if aa_class and aa_class in profile.aa_composition:
                aa_scaler = profile.aa_composition[aa_class]
            elif "methionine" in profile.aa_composition and "methionine" in aa.lower():
                aa_scaler = profile.aa_composition["methionine"]

        if aa_class == "lysine":
            corrected_aa[aa] = conc * corr.lysine_accessibility * aa_scaler
        elif aa_class == "cysteine":
            corrected_aa[aa] = conc * corr.cysteine_accessibility * aa_scaler
        else:
            corrected_aa[aa] = conc * ((corr.lysine_accessibility + corr.cysteine_accessibility) / 2.0) * aa_scaler

    corrected_volatiles = {}
    for compound, conc in predicted_concentrations.items():
        base_retention = corr.volatile_retention
        
        if profile:
            # Base capacity shift for desirable meatiness
            scaler = profile.meaty_potential_multiplier
            
            # Heuristic: methoxypyrazine non-correctable
            if "methoxypyrazine" in compound.lower():
                scaler = min(scaler, profile.methoxypyrazine_ceiling)

            # Heuristic: Off-note penalty + LOX flag
            compound_family = classify_volatile_matrix_family(compound)
            if compound_family in ["aldehyde", "alcohol"] and profile.lox_activity_flag:
                scaler *= (1.0 + profile.off_note_penalty)

            base_retention *= scaler

        corrected_volatiles[compound] = conc * base_retention

    return corrected_volatiles, corrected_aa

