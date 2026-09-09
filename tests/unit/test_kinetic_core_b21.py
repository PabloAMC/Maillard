"""Wave B21 (2026-09-09): the aqueous glucosone route to glyoxal (kinetic_core_b21_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import network, operative_parameters
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, b1_fitted, predict
from src.kinetic_core.parameters_dicarbonyl import (
    AQUEOUS_GLYOXAL_COORDINATES, AQUEOUS_GLYOXAL_KEYS, AQUEOUS_GLYOXAL_PARAMETERS, DICARBONYL_PARAMETERS, EA_AMA_G_KJ_MOL,
    EA_G_GO_AQUEOUS_KJ_MOL, FROZEN_B21, GLASS_K_G_GO, with_aqueous_glyoxal,
)

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b21_fit_report.json"
SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b21_ship_rule.json"


def test_the_step_is_on_the_trunk_only_and_balances():
    assert network.AQUEOUS_GLYOXAL_REACTION_KEYS == ("r_ama_g",)
    r = next(r for r in network.TRUNK_REACTIONS if r.key == "r_ama_g")
    assert r.reactants == {"AMA": 1} and r.products == {"G": 1, "Gly": 1} and r.parameter_key == "k_ama_g"
    network.validate_balance(network.TRUNK_REACTIONS)
    from src.kinetic_core.sulfur import FULL_REACTION_KEYS
    assert "r_ama_g" not in FULL_REACTION_KEYS


def test_the_frozen_literals_match_the_fit_report_and_the_operative_set_reads_them():
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    for key in AQUEOUS_GLYOXAL_COORDINATES:
        assert FROZEN_B21[key] == pytest.approx(report["frozen_parameters"]["aqueous_glyoxal"][key], abs=1e-9), key
    p = operative_parameters(b1_fitted())
    assert p["k_ama_g"].k_ref == pytest.approx(10.0 ** FROZEN_B21["log10_k_ama_g_100C"])
    assert p["k_g_go"].k_ref == pytest.approx(10.0 ** FROZEN_B21["log10_k_g_go_aqueous_100C"])
    assert p["k_g_go"].ea_kj_mol == EA_G_GO_AQUEOUS_KJ_MOL and p["k_ama_g"].ea_kj_mol == EA_AMA_G_KJ_MOL


def test_the_glass_value_stays_as_the_record():
    assert DICARBONYL_PARAMETERS["k_g_go"] is GLASS_K_G_GO
    assert GLASS_K_G_GO.ea_kj_mol == pytest.approx(93.8)
    assert set(AQUEOUS_GLYOXAL_KEYS) == {"k_ama_g", "k_g_go"}
    on = with_aqueous_glyoxal(-2.0, -1.0)
    assert on["k_ama_g"].k_ref == pytest.approx(1e-2) and on["k_g_go"].k_ref == pytest.approx(1e-1)
    for q in AQUEOUS_GLYOXAL_PARAMETERS.values():
        assert "fitted_wave_b21" in q.flags


def test_a_sugar_amine_pot_now_makes_glyoxal_and_says_where_it_came_from():
    spec = FormulationSpec(name="pot", precursors={"D-Glucose": 100.0, "Glycine": 30.0},
                           process=ProcessSpec(thermal=ThermalProgram.isothermal(130.0, 21.0), ph=7.0))
    run = predict(spec, ["glyoxal", "2,5-dimethylpyrazine"])
    assert run.answered and run.concentrations_ug_per_l["glyoxal"] > 1000.0      # about 0.4 mmol/L = 23 mg/L
    assert any(w.startswith("GLYOXAL SUPPLY (B21)") for w in run.declaration.warnings)


def test_the_ship_rule_record_agrees_with_the_report():
    ship = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    assert ship["frozen"] == report["frozen_parameters"]["aqueous_glyoxal"]
    assert ship["verdict"] in ("SHIP", "DO NOT SHIP")
