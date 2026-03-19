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
  docs/Maillard_Plant_based.md and
  docs/Elicit - Maillard Pathways in Plant-Based Cooking - Report.md
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


ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"


def _load_json_payload(file_name: str) -> dict:
    payload_path = DATA_LIT_DIR / file_name
    with open(payload_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


PROCESS_STATE_CALIBRATION_PAYLOAD = _load_json_payload("process_state_calibrations.json")


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
) -> float:
    p_type = _coerce_protein_type(protein_type)
    base = resolve_matrix_correction(p_type, denaturation_state).volatile_retention
    class_factor = get_volatile_class_retention_factor(
        name,
        protein_type=p_type,
        denaturation_state=denaturation_state,
        smiles=smiles,
    )
    return max(0.01, min(1.0, base * class_factor))


def build_matrix_explainability(
    *,
    protein_type: ProteinType | str,
    effective_denaturation_state: float,
    temperature_celsius: float,
    time_minutes: Optional[float],
    pH: Optional[float],
) -> dict[str, object]:
    p_type = _coerce_protein_type(protein_type)
    effective = resolve_matrix_correction(p_type, effective_denaturation_state)
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
        "denaturation_source": DENATURATION_HEURISTICS.get(p_type).source if p_type in DENATURATION_HEURISTICS else "explicit/free",
        "prior_summary": summarize_matrix_prior_bundle(p_type.value),
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
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Scale predicted volatile concentrations and reactive AA concentrations
    by matrix accessibility factors.

    denaturation_state: extrusion/heating increases accessibility.
    1.0 (fully denatured) approaches the free AA model system.
    """
    corr = resolve_matrix_correction(
        protein_type,
        denaturation_state=denaturation_state,
    )

    corrected_aa = {}
    for aa, conc in reactive_amino_acids.items():
        aa_class = _classify_reactive_aa_key(aa)
        if aa_class == "lysine":
            corrected_aa[aa] = conc * corr.lysine_accessibility
        elif aa_class == "cysteine":
            corrected_aa[aa] = conc * corr.cysteine_accessibility
        else:
            corrected_aa[aa] = conc * (corr.lysine_accessibility + corr.cysteine_accessibility) / 2.0

    corrected_volatiles = {
        compound: conc * corr.volatile_retention
        for compound, conc in predicted_concentrations.items()
    }
    return corrected_volatiles, corrected_aa
