from __future__ import annotations

"""
Correction factors for amino acid reactivity in protein matrices
vs. free amino acid model systems.

These are empirical correction factors derived from literature comparing
model system predictions to protein matrix experiments. They should be
recalibrated as you generate your own experimental data.

Current legume-matrix anchoring inside the repo:
- pea: Asen 2022 DSC/DTNB + Malia 2025 Ellman SH envelope, with older repo
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


ACCESSIBILITY_LITERATURE_WINDOWS = {
    ProteinType.PEA_ISOLATE: AccessibilityLiteratureWindow(
        protein_type=ProteinType.PEA_ISOLATE,
        lysine_min=0.30,
        lysine_max=0.45,
        cysteine_min=0.00,
        cysteine_max=0.08,
        source=(
            "Asen 2022 DSC/DTNB + Malia 2025 Ellman SH envelope, retained as a "
            "conservative pea-isolate interpolation while exact benchmark-condition "
            "values remain unmeasured"
        ),
    ),
    ProteinType.SOY_ISOLATE: AccessibilityLiteratureWindow(
        protein_type=ProteinType.SOY_ISOLATE,
        lysine_min=0.36,
        lysine_max=0.50,
        cysteine_min=0.10,
        cysteine_max=0.24,
        source=(
            "Repo soy literature synthesis from Kutzli 2020 / Naik 2021: glycinin/"
            "beta-conglycinin burial, partial reopening after processing, sulfur still"
            " substantially less accessible than free-AA systems"
        ),
    ),
}


DENATURATION_HEURISTICS = {
    ProteinType.PEA_ISOLATE: DenaturationHeuristic(
        protein_type=ProteinType.PEA_ISOLATE,
        midpoint_celsius=82.0,
        width_celsius=14.0,
        time_gain_celsius=6.0,
        time_reference_minutes=10.0,
        acidic_ph_gain_celsius=2.5,
        reference_ph=6.0,
        source=(
            "Asen 2022 DSC/DTNB thermal window + Malia 2025 free-SH response; "
            "calibrated so pea isolate stays mostly native at 40C but opens "
            "progressively by 90-140C"
        ),
    ),
    ProteinType.PEA_CONCENTRATE: DenaturationHeuristic(
        protein_type=ProteinType.PEA_CONCENTRATE,
        midpoint_celsius=84.0,
        width_celsius=14.0,
        time_gain_celsius=5.5,
        time_reference_minutes=10.0,
        acidic_ph_gain_celsius=2.0,
        reference_ph=6.0,
        source=(
            "Pea concentrate inherits the Asen 2022 / Malia 2025 pea thermal envelope "
            "with a modest fiber-burial penalty"
        ),
    ),
    ProteinType.SOY_ISOLATE: DenaturationHeuristic(
        protein_type=ProteinType.SOY_ISOLATE,
        midpoint_celsius=84.0,
        width_celsius=14.0,
        time_gain_celsius=5.0,
        time_reference_minutes=10.0,
        acidic_ph_gain_celsius=2.0,
        reference_ph=6.0,
        source="Repo soy literature synthesis (glycinin/beta-conglycinin burial, partial reopening under heat/process)",
    ),
    ProteinType.SOY_CONCENTRATE: DenaturationHeuristic(
        protein_type=ProteinType.SOY_CONCENTRATE,
        midpoint_celsius=86.0,
        width_celsius=14.0,
        time_gain_celsius=4.5,
        time_reference_minutes=10.0,
        acidic_ph_gain_celsius=1.8,
        reference_ph=6.0,
        source="Soy concentrate inherits the soy-isolate thermal envelope with added concentrate burial",
    ),
    ProteinType.MYCOPROTEIN: DenaturationHeuristic(
        protein_type=ProteinType.MYCOPROTEIN,
        midpoint_celsius=78.0,
        width_celsius=12.0,
        time_gain_celsius=5.0,
        time_reference_minutes=10.0,
        acidic_ph_gain_celsius=1.5,
        reference_ph=6.0,
        source="Conservative mycoprotein thermal-unfolding placeholder pending dedicated calibration",
    ),
}


VOLATILE_CLASS_RETENTION_PROFILES = {
    ProteinType.PEA_ISOLATE: VolatileClassRetentionProfile(
        protein_type=ProteinType.PEA_ISOLATE,
        native_factors={
            "aldehyde": 0.82,
            "furan": 0.92,
            "alcohol": 1.05,
            "sulfur": 1.08,
            "pyrazine": 1.02,
            "other": 1.0,
        },
        denatured_factors={
            "aldehyde": 0.92,
            "furan": 0.97,
            "alcohol": 1.02,
            "sulfur": 1.04,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        source="Conservative class-aware trapping overlay: aldehydes bind most strongly, sulfur volatiles and alcohols escape somewhat more readily",
    ),
    ProteinType.PEA_CONCENTRATE: VolatileClassRetentionProfile(
        protein_type=ProteinType.PEA_CONCENTRATE,
        native_factors={
            "aldehyde": 0.78,
            "furan": 0.88,
            "alcohol": 1.03,
            "sulfur": 1.05,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        denatured_factors={
            "aldehyde": 0.88,
            "furan": 0.94,
            "alcohol": 1.01,
            "sulfur": 1.03,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        source="Pea concentrate inherits stronger aldehyde/furan trapping from the pea family plus a concentrate penalty",
    ),
    ProteinType.SOY_ISOLATE: VolatileClassRetentionProfile(
        protein_type=ProteinType.SOY_ISOLATE,
        native_factors={
            "aldehyde": 0.88,
            "furan": 0.96,
            "alcohol": 1.04,
            "sulfur": 1.06,
            "pyrazine": 1.01,
            "other": 1.0,
        },
        denatured_factors={
            "aldehyde": 0.95,
            "furan": 0.99,
            "alcohol": 1.02,
            "sulfur": 1.03,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        source="Conservative soy overlay: aldehydes still bind more strongly than sulfur/alcohol classes, but less severely than pea",
    ),
    ProteinType.SOY_CONCENTRATE: VolatileClassRetentionProfile(
        protein_type=ProteinType.SOY_CONCENTRATE,
        native_factors={
            "aldehyde": 0.84,
            "furan": 0.93,
            "alcohol": 1.03,
            "sulfur": 1.04,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        denatured_factors={
            "aldehyde": 0.91,
            "furan": 0.97,
            "alcohol": 1.01,
            "sulfur": 1.02,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        source="Soy concentrate inherits the soy family overlay with added burial",
    ),
    ProteinType.MYCOPROTEIN: VolatileClassRetentionProfile(
        protein_type=ProteinType.MYCOPROTEIN,
        native_factors={
            "aldehyde": 0.92,
            "furan": 0.98,
            "alcohol": 1.02,
            "sulfur": 1.03,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        denatured_factors={
            "aldehyde": 0.97,
            "furan": 1.0,
            "alcohol": 1.01,
            "sulfur": 1.01,
            "pyrazine": 1.0,
            "other": 1.0,
        },
        source="Conservative mycoprotein class overlay pending dedicated calibration",
    ),
}

MATRIX_CORRECTIONS = {
    ProteinType.FREE_AMINO_ACID: MatrixCorrection(
        protein_type=ProteinType.FREE_AMINO_ACID,
        lysine_accessibility=1.0,
        cysteine_accessibility=1.0,
        lysine_accessibility_native=1.0,
        lysine_accessibility_denatured=1.0,
        cysteine_accessibility_native=1.0,
        cysteine_accessibility_denatured=1.0,
        volatile_retention=1.0,
        volatile_retention_native=1.0,
        volatile_retention_denatured=1.0,
        source="model system, no correction"
    ),
    ProteinType.PEA_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.PEA_ISOLATE,
        # TNBS/OPA envelope in the repo suggests commercial pea isolates sit well
        # below free-AA accessibility even after wet-heat treatment.
        lysine_accessibility=0.38,
        cysteine_accessibility=0.04,
        lysine_accessibility_native=0.32,
        lysine_accessibility_denatured=0.44,
        cysteine_accessibility_native=0.01,
        cysteine_accessibility_denatured=0.07,
        volatile_retention=0.50,
        volatile_retention_native=0.42,
        volatile_retention_denatured=0.58,
        source=(
            "Asen 2022 DSC/DTNB + Malia 2025 Ellman SH envelope; values remain a "
            "conservative PPI interpolation because the exact 95C pH 5.5 benchmark "
            "condition still lacks direct wet-lab measurement"
        )
    ),
    ProteinType.PEA_CONCENTRATE: MatrixCorrection(
        protein_type=ProteinType.PEA_CONCENTRATE,
        lysine_accessibility=0.28,
        cysteine_accessibility=0.03,
        lysine_accessibility_native=0.22,
        lysine_accessibility_denatured=0.34,
        cysteine_accessibility_native=0.005,
        cysteine_accessibility_denatured=0.055,
        volatile_retention=0.35,
        volatile_retention_native=0.27,
        volatile_retention_denatured=0.43,
        source="Pea-isolate accessibility envelope with added fiber burial penalty"
    ),
    ProteinType.SOY_CONCENTRATE: MatrixCorrection(
        protein_type=ProteinType.SOY_CONCENTRATE,
        lysine_accessibility=0.34,
        cysteine_accessibility=0.15,
        lysine_accessibility_native=0.30,
        lysine_accessibility_denatured=0.38,
        cysteine_accessibility_native=0.11,
        cysteine_accessibility_denatured=0.19,
        volatile_retention=0.45,
        volatile_retention_native=0.37,
        volatile_retention_denatured=0.53,
        source=(
            "Soy concentrate envelope derived from soy-isolate globulin burial "
            "(glycinin/beta-conglycinin) with added fiber/concentrate penalty; "
            "repo anchors: Kutzli 2020 / Naik 2021 synthesis plus soy-polysaccharide "
            "conjugation trapping notes in docs/Maillard_Plant_based.md"
        )
    ),
    ProteinType.SOY_ISOLATE: MatrixCorrection(
        protein_type=ProteinType.SOY_ISOLATE,
        # Soy is kept above pea because the repo literature synthesis treats soy
        # isolates as lysine-rich globulins whose extrusion state can expose more
        # reactive sites, while still remaining well below free-AA behavior and
        # sulfur-limited relative to true meat-like precursor systems.
        lysine_accessibility=0.43,
        cysteine_accessibility=0.17,
        lysine_accessibility_native=0.38,
        lysine_accessibility_denatured=0.48,
        cysteine_accessibility_native=0.14,
        cysteine_accessibility_denatured=0.20,
        volatile_retention=0.55,
        volatile_retention_native=0.47,
        volatile_retention_denatured=0.63,
        source=(
            "Soy isolate envelope anchored to repo literature on glycinin/beta-conglycinin "
            "burial, lysine-rich but sulfur-poor plant-protein composition, and extrusion- "
            "altered reactive-site accessibility; see docs/Maillard_Plant_based.md and "
            "docs/Elicit - Maillard Pathways in Plant-Based Cooking - Report.md "
            "(Kutzli 2020 / Naik 2021 synthesis)"
        )
    ),
    ProteinType.MYCOPROTEIN: MatrixCorrection(
        protein_type=ProteinType.MYCOPROTEIN,
        lysine_accessibility=0.50,
        cysteine_accessibility=0.24,
        lysine_accessibility_native=0.42,
        lysine_accessibility_denatured=0.58,
        cysteine_accessibility_native=0.18,
        cysteine_accessibility_denatured=0.30,
        volatile_retention=0.60,
        volatile_retention_native=0.54,
        volatile_retention_denatured=0.66,
        source="Conservative mycoprotein estimate pending dedicated calibration"
    ),
}


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
