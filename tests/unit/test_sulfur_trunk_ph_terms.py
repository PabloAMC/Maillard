"""sulfur.TRUNK_PH_TERMS (sub-model 1, 2026-10-29): the trunk's measured pH terms inside the sulfur lane, off by default."""
from __future__ import annotations

import pytest

from src.kinetic_core import engine, sulfur, trunk_conditions


def _k(ph_label: float) -> dict:
    return sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), 373.15, ph_label, ph_label=ph_label)


def test_off_by_default_and_then_the_lane_is_ph_blind_on_the_amadori_steps():
    assert sulfur.TRUNK_PH_TERMS is False
    assert _k(5.5)["r_ama_odg"] == _k(6.8)["r_ama_odg"]


def test_on_it_scales_each_step_by_its_measured_slope_on_the_label_scale(monkeypatch):
    monkeypatch.setattr(sulfur, "TRUNK_PH_TERMS", True)
    slope = trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT
    exit_slope = trunk_conditions.THREE_DEOXY_EXIT_PH["k_tdg_fa"][0]
    at_ref, low = _k(trunk_conditions.REFERENCE_PH), _k(5.5)
    for step in ("r_ama_odg", "r_ama_tdg", "r_ama_mgo"):
        assert low[step] / at_ref[step] == pytest.approx(10 ** (slope * (5.5 - 6.8)))
    assert low["r_tdg_fa"] / at_ref["r_tdg_fa"] == pytest.approx(10 ** (exit_slope * (5.5 - 6.8)))
    assert low["r_odg_aa"] == at_ref["r_odg_aa"]   # no measured pH term on this step: untouched


def test_a_fit_overrides_a_slope_by_key_and_nothing_else(monkeypatch):
    monkeypatch.setattr(sulfur, "TRUNK_PH_TERMS", True)
    monkeypatch.setattr(sulfur, "TRUNK_PH_EXPONENT_OVERRIDE", {"k_ama_odg": 0.2, "not_a_step": 9.0})
    at_ref, low = _k(trunk_conditions.REFERENCE_PH), _k(5.5)
    assert low["r_ama_odg"] / at_ref["r_ama_odg"] == pytest.approx(10 ** (0.2 * (5.5 - 6.8)))
    assert low["r_ama_tdg"] / at_ref["r_ama_tdg"] == pytest.approx(10 ** (trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT * (5.5 - 6.8)))


def test_trunk_solutes_enter_the_balance_only_when_switched_on(monkeypatch):
    from src.kinetic_core import ph_state

    drift = engine.core_ph_drift()
    spec = ph_state.BufferSpec(kind="none", declared=True)
    state = {"Gly": 200.0, "AMA": 10.0, "SB": 1.0}
    assert ph_state.titratable_inventory(state, drift, spec) == ()
    monkeypatch.setattr(ph_state, "TRUNK_SOLUTES_IN_BALANCE", True)
    groups = {}
    for g, c in ph_state.titratable_inventory(state, drift, spec):
        groups[g.key] = groups.get(g.key, 0.0) + c
    assert groups["glycine"] == 200.0 and groups["arp"] == 11.0
    # 200 mM glycine buffers like ~0.1 M phosphate near its carboxyl pKa (the review's sizing): capacity > 50 mM/pH
    assert ph_state.buffer_capacity_mmol_per_ph(state, 298.15, 2.4, drift, spec) > 50.0


def test_a_clamped_run_hands_the_tables_its_label_scale_ph(monkeypatch):
    seen = []
    real = sulfur.sulfur_rate_constants_at

    def spy(parameters, t, ph, **kw):
        seen.append(kw.get("ph_label"))
        return real(parameters, t, ph, **kw)

    monkeypatch.setattr(sulfur, "sulfur_rate_constants_at", spy)
    from src.kinetic_core.ph_state import BufferSpec

    sulfur.integrate_sulfur(engine.core_parameters(engine.SULFUR), 363.15, {"AMA": 100.0}, [0.0, 60.0], ph=7.0,
                            buffer_spec=BufferSpec(kind="clamped", declared=True), clamp_ph_label=7.6)
    assert seen and set(seen) == {7.6}
