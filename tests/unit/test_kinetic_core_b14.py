"""Wave B14 (2026-09-07): a declared flat water-activity term on the acrylamide lane, inside its measured window."""
from __future__ import annotations

import pytest

from src.kinetic_core import acrylamide_conditions as ac
from src.kinetic_core import engine
from src.kinetic_core.engine import (
    ACRYLAMIDE, AW_TERM_LANES, AW_TERM_WINDOWS, CoreDraw, FormulationSpec, ProcessSpec, ThermalProgram,
    axis_refusal, declare_envelope, predict,
)

ACR = "Acrylamide"


def _spec(aw=None, name="pot", temp=160.0):
    return FormulationSpec(name=name, precursors={"L-Asparagine": 10.0, "D-Glucose": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(temp, 20.0), ph=6.0, water_activity=aw))


def test_the_band_and_window_are_the_source_table():
    assert ac.AW_WINDOW == (0.88, 0.99)
    assert ac.REFERENCE_AW == 0.92
    assert ac.KF_TABLE[0.92] == (3.57, 1.38)          # the lane's shipped k_int1_acr, Table 2 a_w 0.92 column
    assert ac.AW_SCALE_BAND == (0.41, 1.39)           # 1.45/3.57 .. (3.57+1.38)/3.57
    assert ac.aw_multiplier(None) == 1.0 and ac.aw_multiplier(0.5) == 1.0 and ac.aw_multiplier(0.95) == 1.0
    assert ac.aw_multiplier(0.95, 0.7) == 0.7 and ac.aw_multiplier(0.5, 0.7) == 1.0   # the draw acts only inside


def test_declarations_name_the_window_state():
    assert ac.declarations(_spec().process) == []
    inside = ac.declarations(_spec(0.90).process)
    assert len(inside) == 1 and "inside the measured window" in inside[0] and "FLAT" in inside[0]
    outside = ac.declarations(_spec(0.50).process)
    assert len(outside) == 1 and "OUTSIDE the measured window" in outside[0]


def test_a_prediction_at_the_centre_is_unchanged_and_the_draw_moves_it():
    centre = predict(_spec(), [ACR]).concentrations_ug_per_l[ACR]
    inside = predict(_spec(0.95), [ACR]).concentrations_ug_per_l[ACR]
    assert inside == pytest.approx(centre, rel=1e-9)            # the declared multiplier is exactly 1
    low = predict(_spec(0.95), [ACR], draw=CoreDraw(acrylamide_aw_scale=0.5)).concentrations_ug_per_l[ACR]
    assert 0 < low < inside
    outside = predict(_spec(0.5), [ACR], draw=CoreDraw(acrylamide_aw_scale=0.5)).concentrations_ug_per_l[ACR]
    assert outside == pytest.approx(centre, rel=1e-9)           # no term outside the window, draw or not


def test_the_axis_is_answered_inside_the_window_and_refused_across_it():
    assert ACRYLAMIDE in AW_TERM_LANES and AW_TERM_WINDOWS[ACRYLAMIDE] == ac.AW_WINDOW
    a, b = _spec(0.88, "a"), _spec(0.99, "b")
    da, db = declare_envelope(a, [ACR]), declare_envelope(b, [ACR])
    assert axis_refusal(a, b, da, db) is None
    c = _spec(0.5, "c")
    dc = declare_envelope(c, [ACR])
    reason = axis_refusal(a, c, da, dc)
    assert reason is not None and "measured only inside 0.88-0.99" in reason
    # the envelope declaration prints the term, not the old "metadata only" line
    assert any("WATER ACTIVITY (B14)" in w for w in da.warnings)
    assert not any("METADATA ONLY" in w for w in da.warnings)


def test_the_compare_of_the_two_window_arms_is_flat():
    out = engine.compare(_spec(0.88, "a"), _spec(0.99, "b"), [ACR])
    assert out["comparable"] is True
    row = next(r for r in out["ratios"]["rows"] if r["compound"] == ACR)
    assert row["ratio_a_over_b"] == pytest.approx(1.0, rel=1e-9) and row["direction"] == "equal"
