"""Wave B24 (2026-09-09): 2-acetyl-1-pyrroline from proline, pre-registered, run and REFUSED (kinetic_core_b24_prereg.md)."""
from __future__ import annotations

import json

import pytest

from src import data_paths
from src.kinetic_core import network, operative_parameters, species
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, b1_fitted, predict
from src.kinetic_core.parameters_proline import FROZEN_B24, PROLINE_COORDINATES, PROLINE_KEYS, PROLINE_SHIPPED, with_fitted_proline

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b24_fit_report.json"
SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b24_ship_rule.json"


def _pot(precursors, t_c=100.0, minutes=30.0, ph=7.0):
    return FormulationSpec(name="pot", precursors=precursors, process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes), ph=ph))


def test_the_two_steps_exist_balance_and_are_trunk_only():
    assert set(network.PROLINE_REACTION_KEYS) == {"r_mgo_pro", "r_pyrl_ap"}
    network.validate_balance(network.TRUNK_REACTIONS)
    assert list(species.SPECIES_KEYS)[-3:] == ["PRO", "PYRL", "AP"]
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    assert not {"PRO", "PYRL", "AP"} & set(SULFUR_STATE_KEYS)


def test_the_record_matches_the_report_and_nothing_is_installed():
    report = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    for key in PROLINE_COORDINATES:
        assert FROZEN_B24[key] == pytest.approx(report["frozen_parameters"]["proline"][key], abs=1e-9), key
    ship = json.loads(SHIP_RULE.read_text(encoding="utf-8"))
    assert ship["verdict"] == "DO NOT SHIP" and PROLINE_SHIPPED is False
    p = operative_parameters(b1_fitted())
    for key in PROLINE_KEYS:
        assert p[key].k_ref == 0.0, key
    on = with_fitted_proline(*[FROZEN_B24[k] for k in PROLINE_COORDINATES])
    assert on["k_pyrl_ap"].k_ref == pytest.approx(10.0 ** FROZEN_B24["log10_k_pyrl_ap_100C"]) and on["k_pyrl_ap"].order == 2


def test_the_target_is_refused_with_the_verdict_and_proline_is_charged_as_glycine():
    run = predict(_pot({"L-Proline": 400.0, "methylglyoxal": 4.0}), ["2-acetyl-1-pyrroline"])
    assert not run.answered and any("did not ship" in r for r in run.declaration.reasons)
    hmf = predict(_pot({"L-Proline": 100.0, "D-Glucose": 100.0}, t_c=120.0, minutes=60.0), ["5-HMF"])
    gly = predict(_pot({"Glycine": 100.0, "D-Glucose": 100.0}, t_c=120.0, minutes=60.0), ["5-HMF"])
    assert hmf.answered and any("charged as GLYCINE" in w for w in hmf.declaration.warnings)
    assert hmf.concentrations_ug_per_l["5-HMF"] == pytest.approx(gly.concentrations_ug_per_l["5-HMF"], rel=1e-6)


def test_the_registry_now_carries_the_compound():
    from src import compound_keys
    key = compound_keys.resolve("2-acetyl-1-pyrroline")
    assert key is not None and key.id == "2_acetyl_1_pyrroline"
