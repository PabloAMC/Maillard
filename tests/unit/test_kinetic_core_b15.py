"""Wave B15 (2026-09-07): pH and the dry-side a_w shape on the acrylamide lane, declared from De Vleeschouwer 2006/2007."""
from __future__ import annotations

import pytest

from src.kinetic_core import acrylamide_conditions as ac
from src.kinetic_core import engine
from src.kinetic_core.engine import (
    ACRYLAMIDE, NO_PH_TERM_LANES, PH_TERM_WINDOWS, CoreDraw, FormulationSpec, ProcessSpec, ThermalProgram,
    axis_refusal, declare_envelope, predict,
)

ACR = "Acrylamide"


def _spec(ph=6.8, aw=None, name="pot", temp=160.0):
    return FormulationSpec(name=name, precursors={"L-Asparagine": 10.0, "D-Glucose": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(temp, 20.0), ph=ph, water_activity=aw))


def test_the_ph_slopes_are_de_vleeschouwer_2006_in_decades():
    import math
    # phosphate Table 1: the paper's "log-linear" slopes are NATURAL-log slopes; two-point checks
    assert math.log(ac.KF_PH_2006[8.0][0] / ac.KF_PH_2006[4.0][0]) / 4.0 == pytest.approx(ac.LN_SLOPE_FORMATION, abs=0.01)
    assert math.log(ac.KE_PH_2006[8.0][0] / ac.KE_PH_2006[4.0][0]) / 4.0 == pytest.approx(ac.LN_SLOPE_ELIMINATION, abs=0.01)
    assert ac.PH_EXPONENT_FORMATION == pytest.approx(0.5414 / math.log(10), abs=1e-3)   # 0.235 decades per unit
    assert ac.PH_EXPONENT_ELIMINATION == pytest.approx(0.3442 / math.log(10), abs=1e-3)  # 0.149
    assert ac.PH_EXPONENT_FORMATION_BAND == (0.114, 0.281) and ac.PH_EXPONENT_ELIMINATION_BAND == (0.116, 0.155)
    assert ac.ph_factor(6.8, ac.PH_EXPONENT_FORMATION) == 1.0 and ac.ph_factor(None, 0.5) == 1.0
    assert ac.ph_factor(4.0, ac.PH_EXPONENT_FORMATION) == pytest.approx(10 ** (0.235 * (4.0 - 6.8)))
    assert ac.ph_factor(9.0, 0.5) == ac.ph_factor(8.0, 0.5)             # held at the window edge


def test_the_elimination_shape_is_the_2007_table_joined_to_the_2008_flat_window():
    assert ac.AW_ELIMINATION_TABLE == ((0.34, 0.76), (0.59, 0.6), (0.73, 0.33), (0.82, 0.37), (0.88, 1.0))
    assert ac.aw_elimination_multiplier(None) == 1.0 and ac.aw_elimination_multiplier(0.95) == 1.0
    assert ac.aw_elimination_multiplier(0.88) == 1.0                     # the 2008 window is flat
    assert ac.aw_elimination_multiplier(0.82) == pytest.approx(0.37)
    assert ac.aw_elimination_multiplier(0.34) == pytest.approx(0.76)     # the window floor is the first point
    assert ac.aw_elimination_multiplier(0.30) == 1.0                     # outside the window: no term
    assert ac.aw_elimination_multiplier(0.82, 0.0) == 1.0                # the deficit scale at 0 = no effect
    assert ac.aw_elimination_multiplier(0.82, 1.2) == pytest.approx(1 - 1.2 * 0.63)


def test_apply_scales_only_the_named_steps():
    from src.kinetic_core.engine import core_parameters
    base = core_parameters(ACRYLAMIDE)
    out, decl = ac.apply(base, _spec(ph=4.0, aw=0.82).process)
    for key, p in base.items():
        if key in ("k_asn_glc", "k_acr_dp", "k_int1_acr"):
            continue
        assert out[key] is p, key
    assert out["k_asn_glc"].k_ref == pytest.approx(base["k_asn_glc"].k_ref * ac.ph_factor(4.0, ac.PH_EXPONENT_FORMATION))
    assert out["k_acr_dp"].k_ref == pytest.approx(base["k_acr_dp"].k_ref * ac.ph_factor(4.0, ac.PH_EXPONENT_ELIMINATION) * 0.37)
    assert len(decl) == 2 and "pH (B15)" in decl[1] and "inside the measured window" in decl[0]
    same, decl0 = ac.apply(base, _spec().process)
    assert all(same[k] is base[k] for k in base) and decl0 == []


def test_a_prediction_at_the_references_is_unchanged_and_moves_off_them():
    centre = predict(_spec(), [ACR]).concentrations_ug_per_l[ACR]
    acid = predict(_spec(ph=4.0), [ACR]).concentrations_ug_per_l[ACR]
    basic = predict(_spec(ph=8.0), [ACR]).concentrations_ug_per_l[ACR]
    assert 0 < acid < centre < basic                                     # formation rises 0.54 dex per pH unit
    held = predict(_spec(ph=9.0), [ACR]).concentrations_ug_per_l[ACR]
    assert held == pytest.approx(basic, rel=1e-9)                        # held at pH 8
    drawn = predict(_spec(ph=4.0), [ACR], draw=CoreDraw(acrylamide_ph_exponent_formation=0.114)).concentrations_ug_per_l[ACR]  # the band floor
    assert acid < drawn < centre
    # the dry-side elimination shape: slower elimination at a_w 0.82 than at 0.92 -> more acrylamide
    a82 = predict(_spec(aw=0.82), [ACR]).concentrations_ug_per_l[ACR]
    a92 = predict(_spec(aw=0.92), [ACR]).concentrations_ug_per_l[ACR]
    assert a82 >= a92
    off = predict(_spec(aw=0.82), [ACR], draw=CoreDraw(acrylamide_aw_elimination_scale=0.0)).concentrations_ug_per_l[ACR]
    assert off == pytest.approx(a92, rel=1e-9)


def test_ph_comparisons_are_answered_inside_the_window_and_refused_across_it():
    assert ACRYLAMIDE not in NO_PH_TERM_LANES and PH_TERM_WINDOWS[ACRYLAMIDE] == (4.0, 8.0)
    a, b, c = _spec(ph=4.0, name="a"), _spec(ph=8.0, name="b"), _spec(ph=9.0, name="c")
    da, db, dc = (declare_envelope(s, [ACR]) for s in (a, b, c))
    assert axis_refusal(a, b, da, db) is None
    reason = axis_refusal(a, c, da, dc)
    assert reason is not None and "measured only inside pH 4-8" in reason
    out = engine.compare(a, b, [ACR])
    row = next(r for r in out["ratios"]["rows"] if r["compound"] == ACR)
    assert row["ratio_a_over_b"] < 1.0
