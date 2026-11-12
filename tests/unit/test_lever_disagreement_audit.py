"""Toy-input tests of the disagreement audit (scripts/generators/build_lever_disagreement_audit.py) and of the
precursor banding it relies on (PREREG amendment A5). No engine, no keyed data."""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location(
    "build_lever_disagreement_audit", ROOT / "scripts" / "generators" / "build_lever_disagreement_audit.py")
DA = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(DA)
LT = DA.LT


def _it(lab, y, v=0.0025, **kw):
    base = {"lab": lab, "source": lab + "_paper", "y": y, "v": v, "compound": "X", "lever_class": "pH", "band": "b",
            "quantification": "absolute", "matrix": "model system", "buffered": True, "ph_basis": "initial",
            "digitised": False, "temp_C": 100.0, "time_min": 60.0, "ph": 6.0, "lever_mid": 6.0}
    base.update(kw)
    return base


def test_classes_first_match():
    assert DA.classify([_it("a", 0.5), _it("b", -0.5)]) == "conflict between laboratories"
    assert DA.classify([_it("a", 0.5), _it("a", -0.5), _it("b", 0.4)]) == "sign changes inside one laboratory"
    assert DA.classify([_it("a", 0.5), _it("b", 0.05)]) == "direction agreed"
    assert DA.classify([_it("a", 0.05), _it("b", -0.05)]) == "no clear effect"


def test_pot_sd_floor_keeps_a_tiny_printed_error_from_making_a_conflict():
    # 0.15 dex with a printed SD near zero is still inside two floored SDs (0.2 dex)
    assert DA.classify([_it("a", 0.15, v=1e-8), _it("b", -0.15, v=1e-8)]) == "no clear effect"


def test_stratified_gap_sees_a_condition_that_separates_laboratories():
    rows = [_it("a", 0.5), _it("b", 0.5), _it("c", -0.5, quantification="peak area or relative"),
            _it("d", -0.5, quantification="peak area or relative")]
    g = DA.stratified_gap({("X", "pH", "b"): rows},
                          lambda p, q: None if p["lab"] == q["lab"] else DA._match("quantification", p, q))
    assert g["gap"] == pytest.approx(1.0) and g["pairs_true"] == 2 and g["pairs_false"] == 4


def test_match_is_none_when_a_condition_is_not_recorded():
    assert DA._match("matrix", _it("a", 0.1, matrix=None), _it("b", 0.1)) is None
    assert DA._match("temperature within 10 C", _it("a", 0.1, temp_C=100.0), _it("b", 0.1, temp_C=125.0)) is False


def test_variance_split_recovers_a_pure_laboratory_effect():
    rows = [_it("a", 0.5), _it("a", 0.5), _it("b", -0.5), _it("b", -0.5)]
    v = DA.variance_split({("X", "pH", "b"): rows})
    assert v["sd_within_lab"] == pytest.approx(0.0) and v["share_between_lab"] == pytest.approx(1.0)


def _c(lever, *precursor_sets):
    return ({"lever": lever}, [{"conditions": {"precursors": p}} for p in precursor_sets])


@pytest.mark.parametrize("lever, sets, expected", [
    ("cysteine concentration", ({"cysteine": "1 mM", "ribose": "10 mM"}, {"cysteine": "5 mM", "ribose": "10 mM"}),
     "sulfur donor varied"),
    ("ribose added", ({"thiamine": "30 umol", "ribose": "0"}, {"thiamine": "30 umol", "ribose": "670 umol"}),
     "sugar varied"),
    ("lysine concentration", ({"lysine": "0.2 M"}, {"lysine": "0.4 M"}), "amino acid varied"),
    ("2,4-decadienal amount", ({"2,4-decadienal": "0"}, {"2,4-decadienal": "50 umol"}),
     "lipid oxidation product varied"),
    ("yeast extract added", ({"yeast extract": "0", "water": "200"}, {"yeast extract": "2 g", "water": "198"}),
     "other precursor varied"),
    ("glucose dose (precursors not itemised)", ({}, {}), "sugar varied"),
])
def test_precursor_class_reads_the_varied_precursor(lever, sets, expected):
    c, arms = _c(lever, *sets)
    assert LT.precursor_class(c, arms) == expected


@pytest.mark.parametrize("matrix, expected", [
    ("buffer", "model system"), ("water", "model system"), ("other (acetonitrile/buffer 2:1)", "model system"),
    ("other (Tenebrio molitor mealworm protein powder in water)", "food or protein matrix"),
    ("plant protein (flaxseed protein hydrolysate)", "food or protein matrix"),
    ("buffer + lard (other: lard in phosphate buffer)", "model system"),
])
def test_matrix_class_reads_the_leading_category(matrix, expected):
    c = {"arms": [{"conditions": {"matrix": matrix}}]}
    assert DA.moderators(c)["matrix"] == expected


def test_ph_adjusted_without_buffer_is_not_buffered():
    arm = lambda b: {"arms": [{"conditions": {"buffer": b}}]}
    assert DA.moderators(arm("none (pH adjusted with 3 mol/L NaOH)"))["buffered"] is False
    assert DA.moderators(arm("0.2 M phosphate"))["buffered"] is True
    assert DA.moderators(arm("none stated"))["buffered"] is None
