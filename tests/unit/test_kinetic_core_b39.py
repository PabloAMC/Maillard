"""
WAVES B39 / B40 / B41 (2026-09-11): the fed 3-deoxyglucosone triangle and the formic-acid exit's pH term.

B39 fitted the triangle on Mittelmaier 2011's fed pots and Zhang 2021's ratios and DID NOT SHIP (the
fed peak came three times early). B40 added a pH term to both 3-DG exits and DID NOT SHIP (Leitzen's
methylglyoxal row went 1.28x -> 33x). B41 shipped the formic-acid exit's term with the refit. These
tests hold what shipped: the structure, the literals against the B41 report, the term's scope, the
fed pot's behaviour, and the three ship-rule records.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from src import data_paths
from src.kinetic_core import network, trunk_conditions as TC, uncertainty as U
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, TRUNK, core_parameters, predict
from src.kinetic_core.parameters_dicarbonyl import FED_3DEOXY_COORDINATES, FED_3DEOXY_PARAMETERS, FROZEN_B39, SHIPPED_B39, with_fed_3deoxy

VAL = data_paths.VALIDATION_DIR


def test_the_triangle_is_on_the_trunk_and_balances():
    assert network.FED_3DEOXY_REACTION_KEYS == ("r_ddg_tdg", "r_ddg_dgal", "r_dgal_ddg")
    keys = {r.key: r for r in network.TRUNK_REACTIONS}
    assert keys["r_ddg_tdg"].reactants == {"DDG": 1} and keys["r_ddg_tdg"].products == {"TDG": 1}
    assert keys["r_ddg_dgal"].products == {"DGAL": 1} and keys["r_dgal_ddg"].reactants == {"DGAL": 1}
    network.validate_balance(network.TRUNK_REACTIONS)


def test_the_frozen_literals_are_the_b41_optimum_and_the_engine_reads_them():
    assert SHIPPED_B39
    report = json.loads((VAL / "kinetic_core_b41_fit_report.json").read_text())
    for k in FED_3DEOXY_COORDINATES:
        assert FROZEN_B39[k] == pytest.approx(report["frozen_parameters"]["fed_3deoxy"][k], abs=1e-9), k
    p = core_parameters(TRUNK)
    assert p["k_tdg_ddg"].k_ref == pytest.approx(10.0 ** FROZEN_B39["log10_k_tdg_ddg_100C"])
    assert p["k_ddg_tdg"].k_ref == pytest.approx(10.0 ** FROZEN_B39["log10_k_ddg_tdg_100C"])
    assert p["k_ddg_hmf"].k_ref == pytest.approx(10.0 ** FROZEN_B39["log10_k_ddg_hmf_100C"])
    assert p["k_dgal_ddg"].ea_kj_mol == pytest.approx(36.9) and "barrier_declared_equal_to_forward" in p["k_dgal_ddg"].flags
    assert all("fitted_wave_b39" in FED_3DEOXY_PARAMETERS[k].flags for k in ("k_tdg_ddg", "k_ddg_tdg", "k_ddg_dgal", "k_dgal_ddg", "k_ddg_hmf"))


def test_the_ph_term_is_on_the_formic_acid_exit_only_and_the_rejected_one_is_on_record():
    assert TC.THREE_DEOXY_EXIT_PH_TERM is True
    assert set(TC.THREE_DEOXY_EXIT_PH) == {"k_tdg_fa"}
    assert set(TC.THREE_DEOXY_EXIT_PH_REJECTED_B40) == {"k_tdg_mgo"}
    exp, band, _ = TC.THREE_DEOXY_EXIT_PH["k_tdg_fa"]
    assert band[0] <= exp <= band[1]
    # at the reference pH the factor is exactly one; at pH 5 it is below one, and the declaration says so
    from types import SimpleNamespace
    base = core_parameters(TRUNK)
    at68, _ = TC.apply(dict(base), SimpleNamespace(ph=6.8, water_activity=None))
    at5, w = TC.apply(dict(base), SimpleNamespace(ph=5.0, water_activity=None))
    assert at68["k_tdg_fa"].k_ref == pytest.approx(base["k_tdg_fa"].k_ref)
    assert at5["k_tdg_fa"].k_ref == pytest.approx(base["k_tdg_fa"].k_ref * 10 ** (exp * (5.0 - 6.8)))
    assert at5["k_tdg_mgo"].k_ref == pytest.approx(base["k_tdg_mgo"].k_ref)
    assert any("3-DG EXIT (B40)" in x for x in w)


def test_the_fed_pot_peaks_inside_the_papers_bracket_and_keeps_its_pool():
    """Mittelmaier 2011: 3,4-DGE from 200 uM 3-DG peaks at 30 min (grid 20/30/60) at 26.7 uM."""
    best_t, best_c = None, -1.0
    for t in (10.0, 15.0, 20.0, 23.0, 26.0, 30.0, 40.0, 60.0):
        spec = FormulationSpec(name="fed", precursors={"3-deoxyglucosone": 0.2},
                               process=ProcessSpec(thermal=ThermalProgram.isothermal(120.0, t), ph=5.0, water_activity=0.99, matrix="water"))
        c = predict(spec, ["3,4-dideoxyglucosone", "3-deoxygalactosone"]).concentrations_ug_per_l["3,4-dideoxyglucosone"] / 144.13
        if c > best_c:
            best_t, best_c = t, c
    assert 20.0 <= best_t <= 60.0, best_t
    assert 26.7 / 2 <= best_c <= 26.7 * 2, best_c
    spec = FormulationSpec(name="fed", precursors={"3-deoxyglucosone": 0.2},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(120.0, 60.0), ph=5.0, water_activity=0.99, matrix="water"))
    c = predict(spec, ["3-deoxygalactosone"]).concentrations_ug_per_l["3-deoxygalactosone"]
    assert c > 0.0   # the epimer exists now


def test_the_inert_before_is_still_expressible():
    off = with_fed_3deoxy({"log10_k_tdg_ddg_100C": None, "log10_k_ddg_tdg_100C": -30.0, "log10_k_ddg_dgal_100C": -30.0,
                           "log10_k_dgal_ddg_100C": -30.0, "log10_k_ddg_hmf_100C": None})
    assert off["k_ddg_tdg"].k_ref == 0.0 and off["k_dgal_ddg"].k_ref == 0.0


def test_the_envelope_carries_the_five_rows_and_retired_the_printed_band():
    keys = {p.key for p in U.core_priors()}
    for name in FED_3DEOXY_COORDINATES:
        assert f"b39.{name}" in keys
    assert not any(k.startswith("b34.k_tdg_ddg.") for k in keys)


@pytest.mark.parametrize("wave,verdict", [("b39", "DO NOT SHIP"), ("b40", "DO NOT SHIP"), ("b41", "SHIP")])
def test_the_three_ship_rules_say_what_the_waves_say(wave, verdict):
    art = json.loads((VAL / f"kinetic_core_{wave}_ship_rule.json").read_text())
    assert art["verdict"] == verdict
    before = VAL / f"_{wave}_baseline" / f"core_panel_scores_before_{wave}.json"
    after = VAL / f"_{wave}_baseline" / f"core_panel_scores_after_{wave}.json"
    assert before.exists() and after.exists()
