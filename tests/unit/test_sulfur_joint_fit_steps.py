"""The joint fit's switchable steps (2026-10-30): MGO + cysteine Strecker (and its B59a H2S share), the hexose -> C5
entry, and sub-model 3's three MGO steps -- zero by default, split by the registered share, and element-balanced."""
from __future__ import annotations

import dataclasses

import pytest

from src.kinetic_core import engine, parameters_sulfur as ps, sulfur

STEPS = ("r_mgo_cys_strecker", "r_mgo_cys_h2s", "r_odg_dpo", "r_nf_decay_oh_mgo", "r_tdp_mgo", "r_pent_mgo",
         "r_pent_amine", "r_pama_dpo", "r_pama_tdp", "r_pama_dpo_n", "r_pama_tdp_n", "r_dpo_aa")


def _k(**kw) -> dict:
    return sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), 393.15, 6.0, ph_label=6.0, **kw)


def test_every_new_step_is_zero_by_default_and_listed_inactive():
    k = _k()
    assert all(k[s] == 0.0 for s in STEPS)
    assert set(STEPS) <= set(sulfur.inactive_reactions())


def test_the_mgo_strecker_splits_by_the_registered_h2s_share(monkeypatch):
    monkeypatch.setattr(ps, "MGO_CYSTEINE_STRECKER", ps.MgoCysteineStrecker(k_l_per_mmol_min_100c_ph5=1e-3))
    monkeypatch.setattr(ps, "STRECKER_SULFUR", dataclasses.replace(ps.STRECKER_SULFUR, enabled=True, h2s_share=0.08))
    k = _k()
    total = k["r_mgo_cys_strecker"] + k["r_mgo_cys_h2s"]
    assert total > 0 and k["r_mgo_cys_h2s"] / total == pytest.approx(0.08)
    monkeypatch.setattr(ps, "STRECKER_SULFUR", dataclasses.replace(ps.STRECKER_SULFUR, enabled=False))
    assert _k()["r_mgo_cys_h2s"] == 0.0


def test_the_hexose_c5_entry_scales_with_its_constant(monkeypatch):
    monkeypatch.setattr(ps, "HEXOSE_C5_ENTRY", ps.HexoseC5Entry(k_per_min_100c=2e-3, ea_kj_mol=100.0))
    assert _k()["r_odg_dpo"] > 2e-3          # 120 C > the 100 C anchor
    assert "r_odg_dpo" not in sulfur.inactive_reactions()


def test_the_pentose_amine_entry_runs_at_arps_rates_with_its_ph_slopes(monkeypatch):
    monkeypatch.setattr(ps, "PENTOSE_AMINE_ENTRY", ps.PentoseAmineEntry(
        k_entry_l_per_mmol_min_100c_ph6=1e-4, entry_ph_slope=0.5, k_dpo_aa_per_min_100c_ph6=1e-2, aa_ph_slope=0.3))
    k6, k8 = _k(), sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), 393.15, 8.0, ph_label=8.0)
    entry = ps.PENTOSE_AMINE_ENTRY
    assert k6["r_pama_dpo"] == pytest.approx(entry.pama_dpo_factor * (k6["r_arp_dpo"] + k6["r_arp_dpo_thermal"]))
    assert k6["r_pama_tdp"] == pytest.approx(entry.pama_tdp_factor * (k6["r_arp_tdp"] + k6["r_arp_tdp_thermal"]))
    assert k6["r_pama_dpo_n"] == 0.0     # amine_return 1: fully catalytic
    assert k8["r_pent_amine"] / k6["r_pent_amine"] == pytest.approx(10.0 ** (0.5 * 2.0), rel=1e-6)
    assert k8["r_dpo_aa"] / k6["r_dpo_aa"] == pytest.approx(10.0 ** (0.3 * 2.0), rel=1e-6)
    assert not {"r_pent_amine", "r_pama_dpo", "r_pama_tdp", "r_dpo_aa"} & set(sulfur.inactive_reactions())


def test_pama_runs_at_desclaux_measured_constants(monkeypatch):
    """The frozen factors put PAMA's enolisations on Desclaux 2006's glycine-xylose constants (pH 6, 100 C)."""
    monkeypatch.setattr(ps, "PENTOSE_AMINE_ENTRY", ps.PentoseAmineEntry(k_entry_l_per_mmol_min_100c_ph6=1e-4))
    k = sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), 373.15, 6.0, ph_label=6.0)
    assert k["r_pama_dpo"] == pytest.approx(2.79e-2, rel=1e-3)
    assert k["r_pama_tdp"] == pytest.approx(5.50e-4, rel=1e-3)


def test_the_amine_return_splits_pama_breakdown(monkeypatch):
    monkeypatch.setattr(ps, "PENTOSE_AMINE_ENTRY", ps.PentoseAmineEntry(k_entry_l_per_mmol_min_100c_ph6=1e-4,
                                                                         amine_return=0.25))
    k = _k()
    assert k["r_pama_dpo"] / (k["r_pama_dpo"] + k["r_pama_dpo_n"]) == pytest.approx(0.25)
    assert k["r_pama_tdp_n"] / (k["r_pama_tdp"] + k["r_pama_tdp_n"]) == pytest.approx(0.75)


@pytest.mark.parametrize("key", ["r_mgo_cys_strecker", "r_mgo_cys_h2s", "r_odg_dpo", "r_pent_amine", "r_pama_dpo",
                                 "r_pama_tdp", "r_pama_dpo_n", "r_pama_tdp_n", "r_dpo_aa"])
def test_the_new_steps_conserve_carbon_nitrogen_and_sulfur(key):
    from src.kinetic_core.sulfur import FULL_REACTIONS

    reaction = next(r for r in FULL_REACTIONS if r.key == key)
    for element in ("carbon", "nitrogen", "sulfur"):
        before = sum(total_element_of(k, element) * n for k, n in reaction.reactants.items())
        after = sum(total_element_of(k, element) * n for k, n in reaction.products.items())
        assert before == pytest.approx(after), (key, element)


def total_element_of(species: str, element: str) -> float:
    import numpy as np
    from src.kinetic_core.sulfur import SULFUR_INDEX, N_SULFUR_STATE, total_element

    y = np.zeros(N_SULFUR_STATE)
    y[SULFUR_INDEX[species]] = 1.0
    return float(total_element(y, element))
