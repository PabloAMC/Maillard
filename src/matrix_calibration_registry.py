from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass(frozen=True)
class MatrixCalibrationRecord:
    protein_type: str
    process_state: str
    compound: str
    observable_factor: float
    evidence_strength: str
    source: str
    fallback_mode: str
    notes: str = ""


@dataclass(frozen=True)
class MatrixRuntimeCompositionRule:
    protein_type: str
    compound: str
    mode: str
    active_process_states: tuple[str, ...]
    source: str
    notes: str = ""


_MATRIX_CALIBRATION_RECORDS = (
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        observable_factor=1.0,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
        notes="Reference compound for the pea matrix-only intake/headspace lane.",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="2-pentylfuran",
        observable_factor=1.0,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
        notes="Reference furan marker for the pea matrix-only intake/headspace lane.",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="1-hexanol",
        observable_factor=1.0,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="ambient_slurry",
        compound="nonanal",
        observable_factor=1.0,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 pea isolate ambient slurry baseline",
        fallback_mode="compound_specific",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="hexanal",
        observable_factor=0.22877612093571738,
        evidence_strength="conditional_literature_anchored",
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT aldehyde anchor carried onto the matrix-only oxidation/headspace lane. This is a process-state-specific observable correction, not a global oxidation law.",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="2-pentylfuran",
        observable_factor=0.019473307397293472,
        evidence_strength="conditional_literature_anchored",
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT furan anchor carried onto the matrix-only oxidation/headspace lane.",
    ),
    MatrixCalibrationRecord(
        protein_type="pea_iso",
        process_state="heated_matrix",
        compound="nonanal",
        observable_factor=0.009595650239086601,
        evidence_strength="conditional_literature_anchored",
        source="Trikusuma 2019 UHT pea beverage heated headspace anchor",
        fallback_mode="compound_specific_process_state",
        notes="Heated pea UHT nonanal anchor carried onto the matrix-only oxidation/headspace lane.",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="hexanal",
        observable_factor=0.453 / 0.205,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
        notes="Soy release ratio carried relative to the pea reference intake lane.",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="2-pentylfuran",
        observable_factor=2.972 / 0.502,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="1-hexanol",
        observable_factor=0.143 / 0.063,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="ambient_slurry",
        compound="nonanal",
        observable_factor=0.160 / 0.150,
        evidence_strength="literature_anchored",
        source="Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
        fallback_mode="compound_specific",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="heated_matrix",
        compound="hexanal",
        observable_factor=(0.453 / 0.205) * (1.0 - 0.7060),
        evidence_strength="conditional_literature_anchored",
        source="Shu 2024 heated soy off-flavour attenuation carried onto the Pratap-Singh soy ambient baseline",
        fallback_mode="compound_specific_process_state",
        notes="High-severity soy treatment prior for heated matrix states. Useful for reliability and directional accuracy, but not a meaty benchmark anchor.",
    ),
    MatrixCalibrationRecord(
        protein_type="soy_iso",
        process_state="heated_matrix",
        compound="2-pentylfuran",
        observable_factor=(2.972 / 0.502) * 0.03,
        evidence_strength="conditional_literature_anchored",
        source="Shu 2024 heated soy 2-pentylfuran below-detection carryover onto the Pratap-Singh soy ambient baseline",
        fallback_mode="compound_specific_process_state",
        notes="Below-detection is carried as a small non-zero censoring surrogate so heated soy ranking stays numerically stable while honoring the reported severe attenuation.",
    ),
)


_MATRIX_RUNTIME_COMPOSITION_RULES = (
    MatrixRuntimeCompositionRule(
        protein_type="soy_iso",
        compound="hexanal",
        mode="compose_dynamic_retention",
        active_process_states=("intermediate_matrix", "heated_matrix"),
        source="Ince 2024 reversible soy hexanal binding plus Xu 2023 thermal attenuation prior",
        notes="Ambient slurry remains frozen to preserve the historical Pratap-Singh benchmark calibration.",
    ),
)


def _normalize_compound(name: str) -> str:
    return str(name).strip().lower()


def determine_matrix_process_state(*, temperature_celsius: float, time_minutes: float) -> str:
    if temperature_celsius <= 55.0 and time_minutes <= 30.0:
        return "ambient_slurry"
    if temperature_celsius >= 110.0 or time_minutes >= 90.0:
        return "heated_matrix"
    return "intermediate_matrix"


def get_matrix_calibration_record(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Optional[MatrixCalibrationRecord]:
    if not protein_type:
        return None
    normalized = _normalize_compound(compound)
    requested_state = process_state or "ambient_slurry"

    for record in _MATRIX_CALIBRATION_RECORDS:
        if record.protein_type == protein_type and record.process_state == requested_state and _normalize_compound(record.compound) == normalized:
            return record

    for record in _MATRIX_CALIBRATION_RECORDS:
        if record.protein_type == protein_type and _normalize_compound(record.compound) == normalized:
            return MatrixCalibrationRecord(
                protein_type=record.protein_type,
                process_state=requested_state,
                compound=record.compound,
                observable_factor=record.observable_factor,
                evidence_strength="process_state_mismatch",
                source=record.source,
                fallback_mode="nearest_process_state",
                notes=f"Requested process state '{requested_state}' falls back to '{record.process_state}'.",
            )
    return None


def get_matrix_runtime_composition_policy(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Dict[str, str]:
    if not protein_type:
        return {
            "mode": "static_observable_calibration",
            "source": "none",
            "notes": "No protein-type-specific runtime composition policy is registered.",
        }

    normalized = _normalize_compound(compound)
    requested_state = process_state or "ambient_slurry"
    for rule in _MATRIX_RUNTIME_COMPOSITION_RULES:
        if rule.protein_type != protein_type:
            continue
        if _normalize_compound(rule.compound) != normalized:
            continue
        if requested_state in rule.active_process_states:
            return {
                "mode": rule.mode,
                "source": rule.source,
                "notes": rule.notes,
            }

    return {
        "mode": "static_observable_calibration",
        "source": "historical_calibration_default",
        "notes": "Observable calibration is used as-is for this compound/process-state pair.",
    }


def describe_matrix_calibration(
    compound: str,
    *,
    protein_type: Optional[str],
    process_state: Optional[str],
) -> Dict[str, object]:
    record = get_matrix_calibration_record(
        compound,
        protein_type=protein_type,
        process_state=process_state,
    )
    if record is None:
        return {
            "calibration_source": "class_fallback",
            "calibration_process_state": process_state or "unknown",
            "calibration_evidence_strength": "heuristic",
            "calibration_fallback_mode": "class_level",
            "calibration_observable_factor": None,
            "calibration_notes": "No compound-specific matrix calibration is registered for this compound/process state.",
        }
    return {
        "calibration_source": record.source,
        "calibration_process_state": record.process_state,
        "calibration_evidence_strength": record.evidence_strength,
        "calibration_fallback_mode": record.fallback_mode,
        "calibration_observable_factor": float(record.observable_factor),
        "calibration_notes": record.notes,
    }
