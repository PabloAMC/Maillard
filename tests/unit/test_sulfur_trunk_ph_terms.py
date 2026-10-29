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
