"""Wave B12 (2026-09-07): water activity and pH on the trunk lane, declared from measured ratios."""
from __future__ import annotations

import pytest

from src.kinetic_core import engine, trunk_conditions as tc
from src.kinetic_core.engine import (
    AW_TERM_LANES, NO_PH_TERM_LANES, TRUNK, FormulationSpec, ProcessSpec, ThermalProgram,
    axis_refusal, core_parameters, declare_envelope, predict,
)


def _spec(ph=6.8, aw=None, name="pot"):
    return FormulationSpec(name=name, precursors={"D-Glucose": 100.0, "Glycine": 100.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(120.0, 60.0), ph=ph, water_activity=aw))


def test_the_factors_are_exactly_one_at_the_references():
    for aw in (None, 0.98, 1.0):
        assert tc.aw_multiplier(aw) == 1.0
    assert tc.ph_factor(6.8) == 1.0 and tc.ph_factor(None) == 1.0
    m, f, w = tc.factors(_spec().process)
    assert (m, f, w) == (1.0, 1.0, [])


def test_the_aw_table_is_pereyra_gonzales_normalised_to_solution():
    # Table 1 ratios at 50 and 60 C, averaged
    assert tc.aw_multiplier(0.43) == pytest.approx(3.73)
    assert tc.aw_multiplier(0.69) == pytest.approx(3.68)
    assert tc.aw_multiplier(0.85) == pytest.approx(2.58)
    assert tc.aw_multiplier(0.33) == pytest.approx(3.05)
    assert tc.aw_multiplier(0.20) == tc.aw_multiplier(0.33)          # held below the table
    assert 1.0 < tc.aw_multiplier(0.90) < tc.aw_multiplier(0.85)     # interpolated
    lo, hi = tc.aw_band(0.6)
    assert lo == 1.0 and hi == pytest.approx(1.0 + (tc.aw_multiplier(0.6) - 1.0) * 1.2)


def test_the_ph_exponent_is_the_martins_2003_contrast():
    # k(6.8)/k(5.5) per step at 100 and 120 C, decades per pH unit over 1.3 units
    import math
    ratios = (0.57 / 0.19, 1.56 / 0.10, 1.55 / 0.18, 8.89 / 1.11, 6.29 / 0.86, 8.62 / 0.88)
    exps = [math.log10(r) / 1.3 for r in ratios]
    assert min(exps) == pytest.approx(tc.PH_EXPONENT_BAND[0], abs=0.01)
    assert max(exps) == pytest.approx(tc.PH_EXPONENT_BAND[1], abs=0.01)
    assert sum(exps) / len(exps) == pytest.approx(tc.PH_EXPONENT_DECADES_PER_UNIT, abs=0.01)
    assert tc.ph_factor(5.5) == pytest.approx(10 ** (-0.69 * 1.3))


def test_apply_scales_only_the_named_steps():
    base = core_parameters(TRUNK)
    out, warns = tc.apply(base, _spec(ph=5.5, aw=0.6).process)
    for key in tc.AW_STEPS:
        assert out[key].k_ref == pytest.approx(base[key].k_ref * tc.aw_multiplier(0.6))
    for key in tc.PH_STEPS:
        assert out[key].k_ref == pytest.approx(base[key].k_ref * tc.ph_factor(5.5))
    # B18 (2026-09-08): the pyrazine step's own pH term scales its two Strecker steps, reference 6.8
    from src.kinetic_core.parameters_pyrazine import PYRAZINE_PH_STEPS, pyrazine_ph_factor
    for key in PYRAZINE_PH_STEPS:
        assert out[key].k_ref == pytest.approx(base[key].k_ref * pyrazine_ph_factor(5.5))
    # B22 (2026-09-09): methionine's Strecker steps are glycine's times a ratio and take the same term
    from src.kinetic_core.parameters_methionine import METHIONINE_PH_STEPS
    for key in METHIONINE_PH_STEPS:
        assert out[key].k_ref == pytest.approx(base[key].k_ref * pyrazine_ph_factor(5.5))
    from src.kinetic_core.parameters_proline import PROLINE_PH_STEPS
    for key in PROLINE_PH_STEPS:
        assert out[key].k_ref == pytest.approx(base[key].k_ref * pyrazine_ph_factor(5.5))
    # B41 (2026-09-11): the formic-acid exit from 3-DG carries its own pH term (Martins 2003 k6), declared
    # in tc.THREE_DEOXY_EXIT_PH; it scales by 10^(exponent * (pH - 6.8)) and nothing else new moves.
    for key, (exponent, _band, _src) in tc.THREE_DEOXY_EXIT_PH.items():
        assert out[key].k_ref == pytest.approx(base[key].k_ref * 10 ** (exponent * (5.5 - 6.8)))
    untouched = [k for k in base if k not in tc.AW_STEPS + tc.PH_STEPS + PYRAZINE_PH_STEPS + METHIONINE_PH_STEPS + PROLINE_PH_STEPS + tuple(tc.THREE_DEOXY_EXIT_PH)]
    assert untouched
    for key in untouched:
        assert out[key] == base[key], key
    assert any("WATER ACTIVITY TERM" in w for w in warns) and any("pH TERM" in w for w in warns)
    assert any("PYRAZINE pH TERM" in w for w in warns)
    same, none = tc.apply(base, _spec().process)
    assert same == dict(base) and none == []


def test_a_trunk_prediction_at_the_references_is_unchanged_and_moves_off_them():
    ref = predict(_spec(), ["5-HMF"])
    again = predict(_spec(), ["5-HMF"])
    assert ref.require("5-HMF") == again.require("5-HMF")
    assert not any("TERM (B12)" in w for w in ref.declaration.warnings)
    low_aw = predict(_spec(aw=0.6), ["5-HMF"])
    assert low_aw.require("5-HMF") != ref.require("5-HMF")
    assert any("WATER ACTIVITY TERM (B12)" in w for w in low_aw.declaration.warnings)
    acid = predict(_spec(ph=5.5), ["5-HMF"])
    assert any("pH TERM (B12)" in w for w in acid.declaration.warnings)
    assert acid.run_metadata.get("condition_terms")


def test_the_axes_are_now_answerable_on_the_trunk_and_still_refused_elsewhere():
    assert TRUNK not in NO_PH_TERM_LANES and TRUNK in AW_TERM_LANES
    a, b = _spec(aw=0.3, name="a"), _spec(aw=0.6, name="b")
    da, db = declare_envelope(a, ["5-HMF"]), declare_envelope(b, ["5-HMF"])
    assert axis_refusal(a, b, da, db) is None
    c, d = _spec(ph=5.5, name="c"), _spec(ph=6.8, name="d")
    assert axis_refusal(c, d, declare_envelope(c, ["5-HMF"]), declare_envelope(d, ["5-HMF"])) is None
    acr = lambda aw, n: FormulationSpec(name=n, precursors={"L-Asparagine": 10.0, "D-Glucose": 10.0},  # noqa: E731
                                        process=ProcessSpec(thermal=ThermalProgram.isothermal(150.0, 20.0), ph=6.0, water_activity=aw))
    e, f = acr(0.3, "e"), acr(0.6, "f")
    reason = axis_refusal(e, f, declare_envelope(e, ["Acrylamide"]), declare_envelope(f, ["Acrylamide"]))
    assert reason and "WATER ACTIVITY" in reason
