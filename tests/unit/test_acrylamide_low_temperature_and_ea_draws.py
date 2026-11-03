"""CF-039 / CF-106 (2026-11-03): the acrylamide lane refuses programmes peaking inside 60-120 C, and the
envelope draws the two measured Ea (k_int1_acr across laboratories, k_asn_glc from its CI) with k(160 C) held."""
from __future__ import annotations

import pytest

from src.kinetic_core import acrylamide_conditions as ac
from src.kinetic_core.engine import (
    CoreDraw, FormulationSpec, ProcessSpec, ThermalProgram, declare_envelope, predict,
)
from src.kinetic_core.parameters_acrylamide import MEASURED_ACRYLAMIDE
from src.kinetic_core.uncertainty import SECOND_LABORATORY_BAND, core_priors

ACR = "Acrylamide"


def _spec(temp: float, minutes: float = 20.0) -> FormulationSpec:
    return FormulationSpec(name="pot", precursors={"L-Asparagine": 10.0, "D-Glucose": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(temp, minutes), ph=6.8))


def _refused(temp: float) -> bool:
    return any("CF-106" in r for r in declare_envelope(_spec(temp), [ACR]).reasons)


def test_the_lane_is_refused_only_inside_the_low_temperature_interval():
    assert ac.LOW_T_REFUSAL_C == (60.0, 120.0)
    assert _refused(90.0) and _refused(119.0)
    assert not _refused(120.0) and not _refused(160.0)     # the fit window is answered
    assert not _refused(25.0) and not _refused(60.0)       # the near-zero answer at <= 60 C stands


def test_a_programme_that_reaches_the_fit_window_is_answered():
    ramp = FormulationSpec(name="ramp", precursors={"L-Asparagine": 10.0, "D-Glucose": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram(((30.0, 100.0), (5.0, 170.0))), ph=6.8))
    assert not any("CF-106" in r for r in declare_envelope(ramp, [ACR]).reasons)


def test_the_two_ea_priors_are_declared_and_centred_on_the_shipped_constants():
    priors = {p.key: p for p in core_priors()}
    int1 = priors["acrylamide.k_int1_acr.ea_kj_mol"]
    assert int1.distribution == "uniform" and int1.band == ac.EA_INT1_BAND_KJ_MOL == (95.0, 189.0)
    assert int1.centre == MEASURED_ACRYLAMIDE["k_int1_acr"].ea_kj_mol and int1.reason == SECOND_LABORATORY_BAND
    asn = priors["acrylamide.k_asn_glc.ea_kj_mol"]
    assert asn.distribution == "normal" and asn.centre == MEASURED_ACRYLAMIDE["k_asn_glc"].ea_kj_mol
    assert asn.sigma == pytest.approx(MEASURED_ACRYLAMIDE["k_asn_glc"].ea_ci95_kj_mol / 1.96)


def test_an_ea_draw_holds_k_at_160_c_and_moves_the_temperature_law():
    def acr(temp, **draw):
        return predict(_spec(temp), [ACR], draw=CoreDraw(**draw) if draw else None).concentrations_ug_per_l[ACR]

    at_ref = acr(160.0)
    assert acr(160.0, acrylamide_ea_int1=95.0) == pytest.approx(at_ref, rel=1e-6)
    assert acr(160.0, acrylamide_ea_int1=189.0) == pytest.approx(at_ref, rel=1e-6)
    # the yield need not be monotone in Ea (a faster formation step also exposes acrylamide to elimination
    # sooner), so the temperature law is checked on the rate itself: above T_ref a higher barrier is faster
    k = lambda ea: ac.apply(MEASURED_ACRYLAMIDE, _spec(190.0).process, ea_overrides={"k_int1_acr": ea})[0]["k_int1_acr"]
    t_ref, t_hot = 160.0 + 273.15, 190.0 + 273.15
    assert k(95.0).k_at(t_ref) == pytest.approx(MEASURED_ACRYLAMIDE["k_int1_acr"].k_at(t_ref), rel=1e-12)
    assert k(95.0).k_at(t_hot) < MEASURED_ACRYLAMIDE["k_int1_acr"].k_at(t_hot) < k(189.0).k_at(t_hot)
    assert acr(190.0, acrylamide_ea_int1=95.0) != pytest.approx(acr(190.0), rel=1e-3)   # and the draw does move it
    # the centre (no draw) is the shipped prediction: CoreDraw() with None Ea changes nothing
    assert acr(190.0) == pytest.approx(predict(_spec(190.0), [ACR], draw=CoreDraw()).concentrations_ug_per_l[ACR], rel=1e-9)
