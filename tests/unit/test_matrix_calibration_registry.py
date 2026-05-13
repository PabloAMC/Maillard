from __future__ import annotations

import json

from src.matrix_calibration_registry import describe_matrix_calibration, get_matrix_calibration_record


def test_runtime_multiplier_scales_compound_specific_record(monkeypatch):
    monkeypatch.setenv("MAILLARD_MATRIX_CALIBRATION_MULTIPLIERS", json.dumps({"soy_iso": 1.5}))

    record = get_matrix_calibration_record(
        "Hexanal",
        protein_type="soy_iso",
        process_state="ambient_slurry",
    )

    assert record is not None
    assert record.observable_factor == 1.5 * (0.453 / 0.205)


def test_runtime_multiplier_scales_class_anchor(monkeypatch):
    monkeypatch.setenv("MAILLARD_MATRIX_CALIBRATION_MULTIPLIERS", json.dumps({"soy_iso": 1.5}))

    calibration = describe_matrix_calibration(
        "2-Methyl-3-furanthiol (MFT)",
        protein_type="soy_iso",
        process_state="ambient_slurry",
    )

    assert calibration["calibration_observable_factor"] == 1.5
    assert calibration["calibration_observable_multiplier"] == 1.5
    assert calibration["calibration_observable_multiplier_source"] == "runtime_override"