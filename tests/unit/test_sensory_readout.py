"""Tests for src.sensory_readout (S24 Tier A.2)."""
from __future__ import annotations

from src.pipeline import FormulationResult, UncertaintyEnvelope
from src.sensory_readout import (
    AXIS_MEATY,
    AXIS_OFF_NOTE,
    AXIS_SAFETY,
    AXIS_UNCLASSIFIED,
    CompoundOAV,
    build_sensory_readout,
    classify_axis,
    compute_oav,
    roll_up_axes,
)


# ---------------------------------------------------------------------------
# compute_oav
# ---------------------------------------------------------------------------

def test_compute_oav_basic_ratio():
    assert compute_oav(2.0, 0.5) == 4.0
    assert compute_oav(0.0001, 0.0001) == 1.0


def test_compute_oav_returns_none_when_odt_missing_or_invalid():
    assert compute_oav(1.0, None) is None
    assert compute_oav(1.0, 0.0) is None
    assert compute_oav(1.0, -1.0) is None
    assert compute_oav(1.0, float("nan")) is None


def test_compute_oav_returns_none_when_predicted_invalid():
    assert compute_oav(None, 0.5) is None
    assert compute_oav(float("nan"), 0.5) is None


def test_compute_oav_clamps_negative_predictions_to_zero():
    assert compute_oav(-2.0, 1.0) == 0.0


# ---------------------------------------------------------------------------
# classify_axis
# ---------------------------------------------------------------------------

def test_classify_axis_meaty_off_note_safety():
    assert classify_axis("2-Methyl-3-furanthiol (MFT)") == AXIS_MEATY
    assert classify_axis("Hexanal") == AXIS_OFF_NOTE
    assert classify_axis("Acrylamide") == AXIS_SAFETY
    assert classify_axis("benzaldehyde") == AXIS_UNCLASSIFIED


def test_classify_axis_safety_wins_over_meaty_when_overlap():
    # A hypothetical name that contains both keywords; safety must win to
    # match the existing _suggest_template ordering in experiment_value.
    assert classify_axis("acrylamide-furanthiol adduct") == AXIS_SAFETY


# ---------------------------------------------------------------------------
# roll_up_axes
# ---------------------------------------------------------------------------

def _row(compound: str, axis: str, oav: float | None, *, above: bool | None = None):
    return CompoundOAV(
        compound=compound,
        axis=axis,
        odour_threshold_ug_per_kg=1.0 if oav is not None else None,
        predicted_ppb=oav if oav is not None else 0.0,
        predicted_p5=None,
        predicted_p50=None,
        predicted_p95=None,
        oav=oav,
        oav_p5=None,
        oav_p95=None,
        above_threshold=(above if above is not None else (oav is not None and oav >= 1.0)),
    )


def test_roll_up_axes_picks_top_contributor_and_counts_above_threshold():
    rows = [
        _row("MFT", AXIS_MEATY, 12.0),
        _row("furfurylthiol", AXIS_MEATY, 0.5),
        _row("hexanal", AXIS_OFF_NOTE, 1.5),
        _row("acrylamide", AXIS_SAFETY, 0.2),
        _row("mystery_no_odt", AXIS_MEATY, None),
    ]
    axes = roll_up_axes(rows)

    assert axes[AXIS_MEATY].compounds_with_odt == 2
    assert axes[AXIS_MEATY].compounds_without_odt == 1
    assert axes[AXIS_MEATY].above_threshold_count == 1
    assert axes[AXIS_MEATY].max_oav == 12.0
    assert axes[AXIS_MEATY].top_contributor == "MFT"

    assert axes[AXIS_OFF_NOTE].above_threshold_count == 1
    assert axes[AXIS_OFF_NOTE].top_contributor == "hexanal"

    assert axes[AXIS_SAFETY].above_threshold_count == 0
    assert axes[AXIS_SAFETY].max_oav == 0.2


def test_roll_up_axes_handles_axis_with_no_compounds():
    axes = roll_up_axes([_row("MFT", AXIS_MEATY, 5.0)])
    assert axes[AXIS_OFF_NOTE].compounds_with_odt == 0
    assert axes[AXIS_OFF_NOTE].max_oav is None
    assert axes[AXIS_OFF_NOTE].top_contributor is None


# ---------------------------------------------------------------------------
# build_sensory_readout
# ---------------------------------------------------------------------------

def _make_envelope(compound: str, p50: float) -> UncertaintyEnvelope:
    return UncertaintyEnvelope(
        compound=compound,
        predicted_ppb=p50,
        predicted_p5=p50 / 3.0,
        predicted_p50=p50,
        predicted_p95=p50 * 3.0,
        ci_level_pct=90,
        support_count=1,
        envelope_source="test",
    )


def test_build_sensory_readout_against_curated_specs():
    # MFT + Hexanal + Acrylamide all have curated ODTs in
    # data/species/{desirable_targets,off_flavour_targets}.yml.
    predicted = {
        "2-Methyl-3-furanthiol (MFT)": 0.5,  # ODT 1e-4 -> OAV ~5000
        "Hexanal": 0.1,                       # ODT > 0 -> sub-threshold
        "Acrylamide": 0.0,                    # OAV exactly 0
        "mystery_compound_xyz": 4.2,          # no ODT, excluded from rollups
    }
    envelopes = {name: _make_envelope(name, ppb) for name, ppb in predicted.items()}
    result = FormulationResult(
        name="fixture",
        target_score=0.0,
        off_flavour_risk=0.0,
        predicted_ppb=predicted,
        uncertainty_envelopes=envelopes,
    )

    readout = build_sensory_readout(result)

    assert len(readout.per_compound) == 4
    by_name = {row.compound: row for row in readout.per_compound}

    mft = by_name["2-Methyl-3-furanthiol (MFT)"]
    assert mft.axis == AXIS_MEATY
    assert mft.odour_threshold_ug_per_kg is not None
    assert mft.oav is not None and mft.oav > 1.0
    assert mft.above_threshold is True
    assert mft.oav_p5 is not None and mft.oav_p95 is not None
    assert mft.oav_p5 <= mft.oav <= mft.oav_p95

    hexanal = by_name["Hexanal"]
    assert hexanal.axis == AXIS_OFF_NOTE
    assert hexanal.oav is not None

    mystery = by_name["mystery_compound_xyz"]
    assert mystery.oav is None
    assert mystery.above_threshold is False

    assert readout.axes[AXIS_MEATY].above_threshold_count >= 1
    assert readout.axes[AXIS_MEATY].top_contributor == "2-Methyl-3-furanthiol (MFT)"

    assert "meaty" in readout.headline
    # Mystery compound triggers the "no curated ODT" note.
    assert any("no curated odour threshold" in n for n in readout.notes)


def test_build_sensory_readout_handles_empty_predictions():
    result = FormulationResult(name="empty", target_score=0.0, off_flavour_risk=0.0)
    readout = build_sensory_readout(result)
    assert readout.per_compound == []
    assert readout.headline == "no predicted compounds"
    assert any("empty" in n for n in readout.notes)
