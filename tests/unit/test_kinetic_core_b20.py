"""Wave B20 (2026-09-09): the glycation arm on protein-bound lysine, trunk-only (kinetic_core_b20_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import network, species
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, predict
from src.kinetic_core.parameters_glycation import (
    EA_CML_LOSS_KJ_MOL, EA_FLP_CML_KJ_MOL, FROZEN_B20, GLYCATION_COORDINATES, GLYCATION_KEYS, GLYCATION_PARAMETERS,
    available_fraction, with_fitted_glycation,
)

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b20_fit_report.json"
SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b20_ship_rule.json"


def _pot(matrix=None, g_per_l=None, t_c=120.0, minutes=30.0):
    return FormulationSpec(name="pot", precursors={"D-Glucose": 150.0, "Glycine": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes), ph=6.8, matrix=matrix or "water",
                                               protein_g_per_l=g_per_l))


def test_the_five_steps_run_on_the_trunk_only_and_balance():
    assert set(network.GLYCATION_REACTION_KEYS) == {"r_glc_lysp", "r_flp_cml", "r_flp_cel", "r_flp_decay", "r_cml_loss"}
    network.validate_balance(network.TRUNK_REACTIONS)
    from src.kinetic_core.sulfur import FULL_REACTION_KEYS
    assert not set(network.GLYCATION_REACTION_KEYS) & set(FULL_REACTION_KEYS)
    by_key = {r.key: r for r in network.GLYCATION_REACTIONS}
    assert by_key["r_glc_lysp"].order == 2 and by_key["r_flp_decay"].products == {"TDG": 1, "LYSP": 1}


def test_the_species_are_appended_last_and_are_trunk_only():
    keys = list(species.SPECIES_KEYS)
    # REWRITTEN 2026-09-10. This asserted a NEGATIVE SLICE, and every later wave shifted it:
    # five wave tests broke at once when B24b appended two species. What a wave actually needs
    # is that its own species come AFTER everything that existed before it -- an ORDERING, not
    # a position -- and an ordering survives any number of later appends.
    _b20 = ["LYSP", "FLP", "CML", "CEL"]
    assert [k for k in keys if k in _b20] == _b20
    assert min(species.INDEX[k] for k in _b20) > species.INDEX["MPZ"]
    assert {"LYSP", "FLP", "CML", "CEL"} <= set(species.TRUNK_ONLY_KEYS)
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    assert not {"LYSP", "FLP", "CML", "CEL"} & set(SULFUR_STATE_KEYS)


def test_the_frozen_literals_match_the_fit_report():
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    for key in GLYCATION_COORDINATES:
        assert FROZEN_B20[key] == pytest.approx(report["frozen_parameters"]["glycation"][key], abs=1e-9), key
    assert set(GLYCATION_KEYS) == {"k_glyc", "k_flp_cml", "k_flp_cel", "k_flp_decay", "k_cml_loss"}


def test_the_barriers_are_declared_not_fitted():
    assert GLYCATION_PARAMETERS["k_flp_cml"].ea_kj_mol == EA_FLP_CML_KJ_MOL == 113.0
    assert GLYCATION_PARAMETERS["k_cml_loss"].ea_kj_mol == EA_CML_LOSS_KJ_MOL == 0.0
    for p in GLYCATION_PARAMETERS.values():
        assert "barrier_declared_not_fitted" in p.flags and p.evidence_class == "derived_from_fit_data"
    on = with_fitted_glycation(-4.0, -3.0, -3.5, -2.0, -1.0)
    assert on["k_glyc"].k_ref == pytest.approx(1e-4) and on["k_glyc"].order == 2
    assert available_fraction((0.4, 1.0)) == pytest.approx(0.7)


def test_a_loaded_pot_answers_cml_and_cel_with_the_caveat_and_an_unloaded_one_is_refused():
    run = predict(_pot("soy_isolate", 30.0), ["CML", "CEL", "fructosyl-lysine"])
    assert run.answered
    c = run.concentrations_ug_per_l
    assert c["CML"] > 0 and c["CEL"] > 0 and c["fructosyl-lysine"] > c["CML"]
    assert any(w.startswith("GLYCATION (B20)") for w in run.declaration.warnings)
    refused = predict(_pot(), ["CML"])
    assert not refused.answered
    assert any("GLYCATION TARGETS need a protein" in r for r in refused.declaration.reasons)


def test_without_a_loading_the_arm_is_inert_so_every_earlier_pot_reproduces():
    a = predict(_pot(), ["HMF", "furfural"])
    b = predict(FormulationSpec(name="pot", precursors={"D-Glucose": 150.0, "Glycine": 10.0},
                                process=ProcessSpec(thermal=ThermalProgram.isothermal(120.0, 30.0), ph=6.8)), ["HMF", "furfural"])
    assert a.concentrations_ug_per_l == b.concentrations_ug_per_l


def test_a_sulfur_lane_request_for_cml_is_refused_by_name():
    spec = FormulationSpec(name="pot", precursors={"D-Ribose": 10.0, "L-Cysteine": 10.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(145.0, 20.0), ph=5.0, matrix="soy_isolate", protein_g_per_l=30.0))
    run = predict(spec, ["CML"])
    assert not run.answered and any("GLYCATION TARGETS" in r for r in run.declaration.reasons)


def test_the_ship_rule_record_agrees_with_the_report():
    ship = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    assert ship["frozen"] == report["frozen_parameters"]["glycation"]
    assert ship["verdict"] in ("SHIP", "DO NOT SHIP")
