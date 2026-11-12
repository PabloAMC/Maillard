"""Toy-input tests of the step-5a display check (scripts/generators/build_step5a_display_check.py): the mapping from a
keyed arm to an engine run, the rival's prediction, and the display rule. No integration, no keyed data."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location(
    "build_step5a_display_check", ROOT / "scripts" / "generators" / "build_step5a_display_check.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def _arm(precursors, ph=6.0, matrix="buffer", buffer="none", temp=100, minutes=60, **kw):
    return {"lever_value": kw.get("lever_value"), "value": kw.get("value", 1.0),
            "conditions": {"ph": ph, "ph_basis": "initial", "temp_C": temp, "time_min": minutes,
                           "precursors": precursors, "buffer": buffer, "matrix": matrix}}


def test_parse_amount_forms():
    assert S.parse_amount("4.0 mmol in 40 mL") == ("molar", 100.0)
    assert S.parse_amount("3.5 g/100 mL") == ("g_per_l", 35.0)
    assert S.parse_amount("0.25%") == ("g_per_l", 2.5)
    assert S.parse_amount("1 M") == ("molar", 1000.0)
    assert S.parse_amount("50 mM") == ("molar", 50.0)
    assert S.parse_amount("0.4 g (amino acid:sugar 1:2 w/w)") == ("mass", 0.4)
    assert S.parse_amount("protein-to-sugar 1:0.5") == ("per_protein", 0.5)
    assert S.parse_amount("unstated (named in title)") is None


def test_pot_total_rules_in_order():
    assert S.pot_total_g({"xylose": "2.0 g", "water": "8 g (total 200 g)"}, "none") == 200.0
    assert S.pot_total_g({"hydrolysate": "1.0 g", "xylose": "0.3 g", "total solids": "10 %"}, "") == pytest.approx(13.0)
    assert S.pot_total_g({"xylose": "1.0 g", "protein": "2.0 g"}, "none (distilled water, 10% w/v solids)") == 30.0
    assert S.pot_total_g({"xylose": "0.1 g"}, "none (2 mL distilled water)") == pytest.approx(2.1)
    assert S.pot_total_g({"hydrolysate": "30 g", "xylose": "0.4 g"}, "none") == pytest.approx(30.4)


def test_pure_solution_maps_to_millimolar():
    spec, why = S.arm_spec(_arm({"D-xylose": "1.5 g/100 mL", "L-cysteine": "10 mM"}), "pH", "front_door")
    assert why == ""
    assert spec.precursors["xylose"] == pytest.approx(1000 * 15.0 / 150.13)
    assert spec.precursors["cysteine"] == 10.0
    assert spec.process.ph == 6.0 and spec.process.buffer.kind == "none" and spec.process.buffer.declared


def test_hydrolysate_is_refused_at_the_front_door_and_left_out_in_the_secondary_reading():
    arm = _arm({"pork hydrolysate": "40 mL", "xylose": "3.5 g/100 mL"}, matrix="other (pork hydrolysate)")
    assert S.arm_spec(arm, "precursor dose", "front_door")[1].startswith("matrix material not on file")
    spec, why = S.arm_spec(arm, "precursor dose", "added_precursors")
    assert why == "" and set(spec.precursors) == {"xylose"}
    clamped, _ = S.arm_spec(arm, "precursor dose", "added_precursors_ph_clamped")
    assert clamped.process.buffer.kind == "clamped"
    # the matrix field alone refuses too, when no material is keyed among the precursors
    assert S.arm_spec(_arm({"xylose": "1 M"}, matrix="meat"), "pH", "front_door")[0] is None


def test_missing_inputs_are_never_filled_in():
    assert S.arm_spec(_arm({"xylose": "1 M"}, ph=None), "time", "front_door")[1] == "pH not printed"
    assert S.arm_spec(_arm({"xylose": "1 M"}, temp=None), "pH", "front_door")[1] == "temperature or time not printed"
    assert S.arm_spec(_arm({"xylose": "1 M", "galactose": "1 M"}), "pH", "front_door")[1] == "unmapped precursor: galactose"
    assert S.arm_spec(_arm({"glycine": "unstated"}), "pH", "front_door")[1].startswith("charge not derivable: glycine")
    # a numeric lever's own value fills its condition (the lever is the pH)
    spec, _ = S.arm_spec(_arm({"xylose": "1 M"}, ph=None, lever_value=5.0), "pH", "front_door")
    assert spec.process.ph == 5.0


def test_a_lever_the_engine_has_no_term_for_is_not_simulable():
    item = {"compound": "Furfural", "lever_class": "sugar identity", "contrast": "toy-1",
            "arms": [_arm({"ribose": "0.1 M", "cysteine": "0.1 M"}), _arm({"xylose": "0.1 M", "cysteine": "0.1 M"})]}
    assert S.classify(item, "front_door")[0].startswith("the engine has no term for this lever")
    assert S.classify({**item, "compound": "Benzaldehyde"}, "front_door")[0] == "no engine target for this odorant"


def test_class_total_is_not_the_unsubstituted_pyrazine():
    assert S.target_alias("Pyrazines") is None
    assert S.target_alias("Unsubstituted pyrazine") == "pyrazine"


def test_rival_prediction_numeric_and_categorical():
    cells = {("X", "pH", "b"): {"fit": {"estimate": -0.4}},
             ("X", "sugar identity", "all"): {"fit": {"reference": "glucose", "levels": {
                 "xylose": {"estimate": 0.5}, "ribose": {"estimate": 0.8}, "maltose": {"reading": "not connected"}}}}}
    num = {"cell": ["X", "pH", "b"], "kind": "numeric"}
    cat = {"cell": ["X", "sugar identity", "all"], "kind": "categorical", "level": "ribose", "against": "xylose"}
    assert S.rival_prediction(cells, num) == -0.4
    assert S.rival_prediction(cells, cat) == pytest.approx(0.3)
    assert S.rival_prediction(cells, {**cat, "against": "glucose"}) == 0.8
    assert S.rival_prediction(cells, {**cat, "against": "maltose"}) is None
    assert S.rival_prediction(cells, {**num, "cell": ["Y", "pH", "b"]}) is None


def _row(lab, m, engine, table=None, censored=0):
    return {"lab": lab, "m": m, "span_steps": 1.0, "engine": engine, "table": table, "censored_arms": censored,
            "err_engine": abs(engine - m), "err_no_change": abs(m), "err_table": abs((table or 0.0) - m)}


def test_display_needs_both_rivals_beaten_in_every_laboratory():
    rng = np.random.default_rng(0)
    good = [_row("a", 1.0, 0.9, 0.5), _row("a", 0.8, 0.8, 0.5), _row("b", 1.2, 1.1, 0.5), _row("b", 0.9, 1.0, 0.5)]
    out = S.score_rows(good, rng)
    assert out["displayed"] and out["engine_minus_no_change"]["reading"] == "beats"
    assert out["engine_minus_table"]["ci90"][1] < 0
    # one laboratory where the table is closer: the interval reaches zero or beyond, and the cell is not displayed
    mixed = good[:2] + [_row("b", 0.5, 0.9, 0.5), _row("b", 0.5, 1.0, 0.5)]
    out = S.score_rows(mixed, rng)
    assert not out["displayed"] and out["engine_minus_table"]["reading"] != "beats"
    # one laboratory is never judged
    out = S.score_rows(good[:2], rng)
    assert not out["displayed"] and out["engine_minus_no_change"]["reading"].startswith("not judged")


def test_engine_worse_than_no_change_loses():
    rows = [_row("a", 0.1, -1.0), _row("b", -0.1, 1.2), _row("a", 0.2, -0.9), _row("b", 0.1, 1.0)]
    out = S.score_rows(rows, np.random.default_rng(0))
    assert out["engine_minus_no_change"]["reading"] == "loses" and not out["displayed"]


def test_censored_rows_enter_the_sign_score_only():
    rows = [_row("a", 1.0, 0.9), _row("b", 1.0, 0.8), _row("b", 2.0, -1.0, censored=1)]
    out = S.score_rows(rows, np.random.default_rng(0))
    assert out["contrasts"] == 2 and out["contrasts_with_a_censored_arm"] == 1
    assert out["sign_accuracy"]["contrasts"] == 3 and out["sign_accuracy"]["engine"] == pytest.approx(2 / 3)
    assert out["sign_accuracy"]["table"] == 0.0        # no table estimate: zero, a miss


def test_wls_slope_matches_two_point_difference():
    assert S.wls_slope([0.0, 2.0], [1.0, 2.0], [0.1, 0.3]) == pytest.approx(0.5)
