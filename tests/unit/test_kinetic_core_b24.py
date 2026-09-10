"""Wave B24 (2026-09-09): 2-acetyl-1-pyrroline from proline, pre-registered, run and REFUSED (kinetic_core_b24_prereg.md)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import network, operative_parameters, species
from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, b1_fitted, predict
from src.kinetic_core.parameters_proline import FROZEN_B24, PROLINE_COORDINATES, PROLINE_KEYS, PROLINE_SHIPPED, with_fitted_proline

FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b24_fit_report.json"
SHIP_RULE = data_paths.VALIDATION_DIR / "kinetic_core_b24_ship_rule.json"


def _pot(precursors, t_c=100.0, minutes=30.0, ph=7.0):
    return FormulationSpec(name="pot", precursors=precursors, process=ProcessSpec(thermal=ThermalProgram.isothermal(t_c, minutes), ph=ph))


REPO = Path(__file__).resolve().parents[2]


def test_the_steps_exist_balance_and_are_trunk_only():
    """WIDENED BY WAVE B24b (2026-09-10): two steps became four, and two species became five."""
    assert set(network.PROLINE_REACTION_KEYS) == {
        "r_mgo_pro", "r_pyrl_ap", "r_pyrl_ha_athp", "r_pyrl_loss"}
    network.validate_balance(network.TRUNK_REACTIONS)
    # 2026-09-10: an ORDERING, not a position. B22b appended two more species after these and a
    # negative slice would break again on the next wave, which is what happened to five other
    # wave tests the same day.
    _mine = ["PRO", "PYRL", "AP", "ACETOL", "ATHP"]
    assert [k for k in species.SPECIES_KEYS if k in _mine] == _mine
    assert min(species.INDEX[k] for k in _mine) > species.INDEX["DMDS"]
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    assert not {"PRO", "PYRL", "AP", "ACETOL", "ATHP"} & set(SULFUR_STATE_KEYS)


def test_b24b_is_inert_and_its_loss_is_refuted_not_merely_unfitted():
    """
    The finding B24b exists for. Its 1-pyrroline loss ran to a bound, and a bound normally means
    "widen it". It was tested instead: a FIRST-ORDER loss moves the source's two experiments
    together at every rate, while the source demands they differ by 428x per pyrroline. The
    structure is refuted, so nothing is installed and the band is not widened.
    """
    from src.kinetic_core.parameters_proline import (
        B24B_SHIPPED, FROZEN_B24B, PROLINE_B24B_COORDINATES, PROLINE_PARAMETERS)

    assert B24B_SHIPPED is False and FROZEN_B24B == {}
    for key in ("k_ha_athp", "k_pyrl_loss"):
        assert PROLINE_PARAMETERS[key].k_ref == 0.0, key
    assert set(PROLINE_B24B_COORDINATES) == {"log10_k_ha_athp_100C", "log10_k_pyrl_loss_100C"}
    ship = json.loads((REPO / "results/validation/kinetic_core_b24b_ship_rule.json").read_text())
    assert ship["verdict"] == "DO NOT SHIP"
    # the ordering of the switch DID come out right; that is the half that worked
    assert ship["T1"]["ordering_increasing"] is True
    assert ship["T1"]["pass"] is False          # and the magnitude did not
    assert ship["T3"]["pass"] is True           # B24's fed rows survived the new sinks


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
