"""Lane E (sprint 2026-05-10b): scope-check gate.

Verifies that the calibration convex hull is honored: in-scope formulations
return ``in_scope=True`` and out-of-scope formulations (new protein types,
non-calibrated process states, out-of-range temperature/time/pH) are flagged
with a structured reason list.
"""
from __future__ import annotations

from src.matrix_calibration_registry import (
    ScopeAssessment,
    is_formulation_in_calibration_scope,
)


def test_pea_iso_heated_within_envelope_is_in_scope():
    assessment = is_formulation_in_calibration_scope(
        protein_type="pea_iso",
        process_state="heated_matrix",
        temperature_celsius=105.0,
        time_minutes=45.0,
        pH=5.5,
    )
    assert isinstance(assessment, ScopeAssessment)
    assert assessment.in_scope is True
    assert assessment.reasons == ()


def test_wheat_gluten_is_out_of_scope_with_explicit_reason():
    assessment = is_formulation_in_calibration_scope(
        protein_type="wheat_gluten",
        process_state="heated_matrix",
        temperature_celsius=150.0,
        time_minutes=30.0,
        pH=6.5,
    )
    assert assessment.in_scope is False
    assert any("wheat_gluten" in reason for reason in assessment.reasons)


def test_extrusion_state_is_out_of_scope_for_pea_iso():
    assessment = is_formulation_in_calibration_scope(
        protein_type="pea_iso",
        process_state="extrusion_structured",
        temperature_celsius=160.0,
        time_minutes=2.0,
        pH=6.0,
    )
    assert assessment.in_scope is False
    nearest = assessment.nearest_calibrated
    assert nearest.get("protein_type") == "pea_iso"
    assert nearest.get("process_state") in {"ambient_slurry", "heated_matrix"}


def test_temperature_outside_envelope_demotes_in_scope_to_false():
    assessment = is_formulation_in_calibration_scope(
        protein_type="pea_iso",
        process_state="heated_matrix",
        temperature_celsius=200.0,  # well above 145°C calibrated upper bound
        time_minutes=30.0,
        pH=6.0,
    )
    assert assessment.in_scope is False
    assert any("temperature" in reason for reason in assessment.reasons)


def test_missing_protein_type_is_out_of_scope():
    assessment = is_formulation_in_calibration_scope(
        protein_type=None,
        process_state="heated_matrix",
    )
    assert assessment.in_scope is False
    assert "protein_type is missing" in assessment.reasons


def test_to_dict_round_trip_is_json_serializable():
    assessment = is_formulation_in_calibration_scope(
        protein_type="myco",
        process_state="heated_matrix",
    )
    payload = assessment.to_dict()
    assert payload["in_scope"] is False
    assert isinstance(payload["reasons"], list)
    assert isinstance(payload["nearest_calibrated"], dict)
