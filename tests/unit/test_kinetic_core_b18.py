"""Wave B18 (2026-09-08): the pyrazine step on the trunk lane, trunk-only (kinetic_core_b18_prereg.md)."""
from __future__ import annotations

import json
import math

import pytest

from src import data_paths
from src.kinetic_core import network, species
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, predict
from src.kinetic_core.parameters_pyrazine import (
    FROZEN_B18, K_COND, PYRAZINE_KEYS, PYRAZINE_PARAMETERS, PYRAZINE_PH_STEPS, PYRAZINE_REFERENCE_PH,
    pyrazine_ph_factor, with_fitted_pyrazine,
)

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b18_fit_report.json"


def _pot(precursors=None, t_c=120.0, minutes=60.0, ph=6.8):
    return FormulationSpec(name="pot", precursors=precursors or {"D-Glucose": 100.0, "Glycine": 100.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes), ph=ph))


def test_the_five_steps_run_on_the_trunk_only_and_balance():
    assert set(network.PYRAZINE_REACTION_KEYS) == {"r_go_ak", "r_mgo_ak", "r_akg_pz", "r_akm_dmp", "r_ak_mpz"}
    assert not set(network.PYRAZINE_REACTION_KEYS) & set(network.REACTION_KEYS)
    network.validate_balance(network.TRUNK_REACTIONS)
    from src.kinetic_core.sulfur import FULL_REACTION_KEYS
    assert not set(network.PYRAZINE_REACTION_KEYS) & set(FULL_REACTION_KEYS)
    cond = [r for r in network.PYRAZINE_REACTIONS if r.parameter_key == "k_cond"]
    assert len(cond) == 3 and all(r.order == 2 for r in cond)


def test_the_species_are_appended_after_every_existing_one_and_are_trunk_only():
    keys = list(species.SPECIES_KEYS)
    assert keys[-13:-8] == ["PZ", "DMP", "MPZ", "AKG", "AKM"]      # B20 appended LYSP / FLP / CML / CEL after them; B22 MET / MTAL / MSH / DMDS
    assert set(("PZ", "DMP", "MPZ", "AKG", "AKM")) <= set(species.TRUNK_ONLY_KEYS)
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    assert not {"PZ", "DMP", "MPZ", "AKG", "AKM"} & set(SULFUR_STATE_KEYS)


def test_the_frozen_literals_match_the_fit_report():
    if not FIT_REPORT.exists():
        pytest.skip("B18 fit report not generated yet")
    reported = json.loads(FIT_REPORT.read_text())["frozen_parameters"]["pyrazine"]
    for key, value in FROZEN_B18.items():
        assert value == pytest.approx(float(reported[key]), rel=1e-12, abs=1e-12), key


def test_exactly_two_constants_are_fitted_and_the_condensation_is_declared():
    fitted = [k for k, p in PYRAZINE_PARAMETERS.items() if p.evidence_class == "derived_from_fit_data"]
    assert fitted == ["k_go_ak", "k_mgo_ak"]
    assert PYRAZINE_PARAMETERS["k_cond"] is K_COND and K_COND.evidence_class == "bounded_from_a_timescale_bracket"
    assert K_COND.ea_kj_mol == 0.0 and "declared_fast" in K_COND.flags
    assert tuple(PYRAZINE_KEYS) == ("k_go_ak", "k_mgo_ak", "k_cond") and PYRAZINE_PH_STEPS == ("k_go_ak", "k_mgo_ak")
    block = with_fitted_pyrazine(-7.0, 100.0, -8.0, 110.0)
    assert block["k_go_ak"].k_ref == pytest.approx(1e-7) and block["k_mgo_ak"].ea_kj_mol == 110.0


def test_the_ph_term_is_one_at_the_reference_and_falls_toward_acid():
    assert pyrazine_ph_factor(None) == 1.0 and pyrazine_ph_factor(PYRAZINE_REFERENCE_PH) == pytest.approx(1.0)
    assert pyrazine_ph_factor(9.0) > pyrazine_ph_factor(8.0) > pyrazine_ph_factor(7.0) > pyrazine_ph_factor(5.0)
    # two slopes: the acid side is steeper than the alkaline side (Leahy's ratios)
    hi = math.log10(pyrazine_ph_factor(9.0) / pyrazine_ph_factor(8.0))
    lo = math.log10(pyrazine_ph_factor(7.0) / pyrazine_ph_factor(6.0))
    assert lo > hi > 0


def test_a_trunk_pot_answers_the_three_pyrazines_with_the_caveats():
    run = predict(_pot(), ["2,5-dimethylpyrazine", "pyrazine", "methylpyrazine", "5-HMF"])
    assert run.answered
    dmp, pz, mpz = (run.require(c) for c in ("2,5-dimethylpyrazine", "pyrazine", "methylpyrazine"))
    # Until B21 the trunk made almost no glyoxal in water and the order was DMP > MPZ > PZ; B21's aqueous
    # glucosone route (2026-09-09) supplies glyoxal and the order is PZ > MPZ > DMP, Xia 2022's direction.
    assert pz > mpz > dmp > 0
    assert sum("PYRAZINES (B18)" in w for w in run.declaration.warnings) == 2
    acid = predict(_pot(ph=5.0), ["2,5-dimethylpyrazine"]).require("2,5-dimethylpyrazine")
    assert acid < dmp


def test_a_sulfur_lane_request_for_a_pyrazine_is_refused_by_name():
    run = predict(_pot(precursors={"D-Ribose": 100.0, "L-Cysteine": 33.0}, t_c=145.0, minutes=20.0, ph=5.0),
                  ["2-methyl-3-furanthiol", "2,5-dimethylpyrazine"])
    assert not run.answered
    assert any("PYRAZINE TARGETS" in r and "trunk" in r for r in run.declaration.reasons), run.declaration.reasons


def test_the_ship_rule_record_agrees_with_the_report():
    ship = data_paths.VALIDATION_DIR / "kinetic_core_b18_ship_rule.json"
    if not ship.exists():
        pytest.skip("B18 ship rule not generated yet")
    payload = json.loads(ship.read_text())
    assert payload["verdict"] == "SHIP" and payload["T1"]["pass"] and payload["T2"]["pass"] and payload["T5"]["pass"]
    assert payload["T2"]["worst_dex"] < 0.05
    assert payload["frozen_parameters"]["pyrazine"]["log10_k_go_ak_100C"] == pytest.approx(FROZEN_B18["log10_k_go_ak_100C"])
