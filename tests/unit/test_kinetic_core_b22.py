"""Wave B22 (2026-09-09): the methionine chain on the sugar path, pre-registered, run and REFUSED (kinetic_core_b22_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import network, operative_parameters, species
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, b1_fitted, predict
from src.kinetic_core.parameters_methionine import (
    FROZEN_B22, INERT_B22, METHIONINE_COORDINATES, METHIONINE_KEYS, METHIONINE_SHIPPED, with_fitted_methionine,
)

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b22_fit_report.json"
SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b22_ship_rule.json"


def _pot(precursors, t_c=120.0, minutes=10.0, ph=6.2):
    return FormulationSpec(name="pot", precursors=precursors, process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes), ph=ph))


def test_the_four_steps_exist_balance_and_are_trunk_only():
    assert set(network.METHIONINE_REACTION_KEYS) == {
        "r_go_met", "r_mgo_met", "r_mtal_msh", "r_msh_dmds",
        "r_marp_mtal", "r_marp_loss",           # wave B22b, 2026-09-10
    }
    network.validate_balance(network.TRUNK_REACTIONS)
    # REWRITTEN 2026-09-10. This asserted a NEGATIVE SLICE, and every later wave shifted it:
    # five wave tests broke at once when B24b appended two species. What a wave actually needs
    # is that its own species come AFTER everything that existed before it -- an ORDERING, not
    # a position -- and an ordering survives any number of later appends.
    _b22 = ["MET", "MTAL", "MSH", "DMDS"]
    _keys = list(species.SPECIES_KEYS)
    assert [k for k in _keys if k in _b22] == _b22
    assert min(species.INDEX[k] for k in _b22) > species.INDEX["CEL"]
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    assert not {"MET", "MTAL", "MSH", "DMDS"} & set(SULFUR_STATE_KEYS)


def test_the_record_matches_the_report_and_nothing_is_installed():
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    for key in METHIONINE_COORDINATES:
        assert FROZEN_B22[key] == pytest.approx(report["frozen_parameters"]["methionine"][key], abs=1e-9), key
    ship = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    assert ship["verdict"] == "DO NOT SHIP" and METHIONINE_SHIPPED is False
    p = operative_parameters(b1_fitted())
    for key in METHIONINE_KEYS:
        assert p[key].k_ref == 0.0, key       # the steps exist and carry no flux
    on = with_fitted_methionine(*[FROZEN_B22[k] for k in METHIONINE_COORDINATES])
    assert on["k_go_met"].k_ref > 0 and on["k_go_met"].order == 2
    assert with_fitted_methionine(*[INERT_B22[k] for k in METHIONINE_COORDINATES])["k_mtal_msh"].k_ref == 0.0


def test_methional_is_refused_with_the_verdict_and_methionine_is_charged_as_glycine():
    run = predict(_pot({"L-Methionine": 0.268, "D-Fructose": 111.0, "D-Glucose": 83.0}), ["methional"])
    assert not run.answered
    assert any("did not ship" in r for r in run.declaration.reasons)
    hmf = predict(_pot({"L-Methionine": 200.0, "D-Glucose": 200.0}, ph=7.5), ["5-HMF"])
    assert hmf.answered and any("charged as GLYCINE" in w for w in hmf.declaration.warnings)
    gly = predict(_pot({"Glycine": 200.0, "D-Glucose": 200.0}, ph=7.5), ["5-HMF"])
    assert hmf.concentrations_ug_per_l["5-HMF"] == pytest.approx(gly.concentrations_ug_per_l["5-HMF"], rel=1e-6)


def test_dimethyl_trisulfide_is_refused_with_its_reason():
    run = predict(_pot({"L-Methionine": 1.0, "D-Glucose": 100.0}), ["dimethyl trisulfide"])
    assert not run.answered and any("hydrogen sulfide" in r for r in run.declaration.reasons)
