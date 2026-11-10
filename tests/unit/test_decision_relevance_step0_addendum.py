"""Toy-input tests of the step-0 addendum (scripts/generators/build_decision_relevance_step0_addendum.py): lever
classification, sign accuracy and effect binning. No engine runs."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "build_decision_relevance_step0_addendum",
    ROOT / "scripts" / "generators" / "build_decision_relevance_step0_addendum.py")
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)


def _spec_of(ph=6.0, temp=120.0, time=60.0, prec=None, buffer="phosphate", matrix="water"):
    thermal = SimpleNamespace(segments=((time, temp),))
    proc = SimpleNamespace(ph=ph, thermal=thermal, buffer=buffer, matrix=matrix)
    return SimpleNamespace(precursors=prec or {"cysteine": 10.0, "ribose": 10.0}, process=proc)


def test_lever_class_single_field():
    assert A.lever_class(_spec_of(ph=5.0), _spec_of(ph=7.0))[0] == "pH"
    assert A.lever_class(_spec_of(time=30), _spec_of(time=60))[0] == "time"
    assert A.lever_class(_spec_of(temp=100), _spec_of(temp=140))[0] == "temperature"
    assert A.lever_class(_spec_of(prec={"cysteine": 5, "ribose": 10}),
                         _spec_of(prec={"cysteine": 10, "ribose": 10}))[0] == "precursor dose"
    assert A.lever_class(_spec_of(prec={"cysteine": 10, "ribose": 10}),
                         _spec_of(prec={"cysteine": 10, "glucose": 10}))[0] == "precursor identity"


def test_lever_class_mixed_and_none():
    cls, fields = A.lever_class(_spec_of(ph=5.0, time=30), _spec_of(ph=7.0, time=60))
    assert cls == "mixed" and fields == ["pH", "time"]
    assert A.lever_class(_spec_of(), _spec_of())[0] == "none detected"


def test_sign_accuracy_is_orientation_invariant():
    rows = [{"p": 0.5, "m": 0.2}, {"p": -0.4, "m": 0.3}, {"p": 0.6, "m": 0.4}]
    flipped = [{"p": -r["p"], "m": -r["m"]} for r in rows]
    assert A.sign_accuracy(rows) == A.sign_accuracy(flipped)
    assert A.sign_accuracy(rows)[f"|p| >= {A.SIGN_CUT_P}"]["right"] == 2


def test_effect_bins_cover_the_line():
    for x in (0.0, 0.049, 0.05, 0.29, 0.3, 0.99, 1.0, 7.0):
        assert A.effect_bin({"p": x}) != "?"
        assert A.effect_bin({"p": -x}) == A.effect_bin({"p": x})
