from __future__ import annotations

import json

import pytest

from src.matrix_calibration_registry import describe_matrix_calibration, get_matrix_calibration_record


def test_runtime_multiplier_scales_compound_specific_record(monkeypatch):
    monkeypatch.setenv("MAILLARD_MATRIX_CALIBRATION_MULTIPLIERS", json.dumps({"soy_iso": 1.5}))

    record = get_matrix_calibration_record(
        "Hexanal",
        protein_type="soy_iso",
        process_state="ambient_slurry",
    )

    assert record is not None
    # RE-PINNED 2026-08-27 (Wave O refit to content-corrected anchors, owner-approved).
    # The soy ambient hexanal factor was refitted 0.453/0.205 = 2.2097561 -> 9.54007 against
    # the paper's verified 1621.71 ppb (one shared scale of 4.317249x across both ambient
    # lanes; results/validation/matrix_observability_refit_pratap_singh.{json,md}).
    # THIS TEST IS ABOUT THE MULTIPLIER, NOT THE CONSTANT, so it is expressed as
    # "the multiplier scales whatever the registry ships" rather than as a literal product.
    # That keeps it green through any future refit while still failing if the multiplier
    # stops being applied -- which is the behaviour under test.
    shipped = get_matrix_calibration_record(
        "Hexanal",
        protein_type="pea_iso",  # any lane the soy_iso multiplier does NOT touch
        process_state="ambient_slurry",
    )
    assert shipped is not None
    # RE-PINNED 2026-08-28 (Wave Y): the constants moved sides. Wave O's shared ambient
    # hexanal scale now lives in MATRIX_BENCHMARK_BASE_MARKER_YIELDS['Hexanal'] and every
    # hexanal observability constant was divided by it, so soy ambient reads 0.453/0.205 and
    # pea ambient reads 1.0. Record: results/validation/matrix_marker_yield_rederivation.md.
    #
    # The comment above already says this test is about the MULTIPLIER rather than the
    # constant, so the assertions are now written that way instead of being re-pinned to new
    # literals -- which makes them survive the next relocation as well as this one.
    assert record.observable_factor == pytest.approx(1.5 * (0.453 / 0.205))
    # and the multiplier really is a multiplier, not a replacement
    assert record.observable_factor != pytest.approx(0.453 / 0.205)
    assert shipped.observable_factor == pytest.approx(1.0), (
        "the pea lane must NOT be scaled by a soy_iso multiplier"
    )


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