"""Toy-input tests of the PREREG A6 field derivations (scripts/generators/build_lever_field_readings.py). No engine,
no keyed data."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.generators import build_lever_field_readings as F  # noqa: E402


def test_control_from_buffer_text():
    assert F.control_from_text("0.2 M phosphate") == "buffered"
    assert F.control_from_text("none (pH adjusted with 3 mol/L NaOH)") == "adjusted_unbuffered"
    assert F.control_from_text("none (distilled water, 10% w/v solids)") == "unadjusted"
    assert F.control_from_text("none") == "unstated"          # a bare 'none' needs a reader check
    assert F.control_from_text("unstated") == "unstated"


def test_check_overrides_derived_control():
    arms = [{"conditions": {"buffer": "none"}}] * 2
    assert F.contrast_control(arms, None, "s-1") == ("unstated", "derived")
    chk = {"fields": ["F1"], "ph_control": "adjusted_unbuffered", "ph_control_overrides": {"s-2": "buffered"}}
    assert F.contrast_control(arms, chk, "s-1") == ("adjusted_unbuffered", "check")
    assert F.contrast_control(arms, chk, "s-2") == ("buffered", "check")
    assert F.contrast_control(arms, None, "s-3", keyed="adjusted_unbuffered") == ("adjusted_unbuffered", "keyed")


def test_final_ph_from_the_template_field():
    assert F.arm_final_ph({"conditions": {"ph": 7.0, "ph_basis": "initial", "final_ph": 6.2}}, None) == 6.2
    assert F.arm_final_ph({"conditions": {"ph": 5.5, "ph_basis": "final"}}, None) == 5.5
    assert F.arm_final_ph({"conditions": {"ph": 7.0, "ph_basis": "initial"}}, None) is None


def test_same_reactant_matches_the_level_not_its_neighbours():
    assert F._same_reactant("glycine", "glycine (2 mmol)")
    assert F._same_reactant("GSH", "glutathione")
    assert F._same_reactant("N-acetyl-L-cysteine", "n-acetylcysteine")
    assert F._same_reactant("glutamic acid", "free amino acids (lysine + glutamic acid)")
    assert not F._same_reactant("ribose 5-phosphate", "phosphate")
    assert not F._same_reactant("glycine", "gly-gly")
    assert not F._same_reactant("cysteine", "n-acetylcysteine")


def test_base_mixture_drops_the_level_and_zero_doses():
    none_arm = {"lever_value": "none", "conditions": {"precursors": {"ribose": "0.1 M", "cysteine": "0 M"}}}
    cys_arm = {"lever_value": "cysteine", "conditions": {"precursors": {"ribose": "0.1 M", "cysteine": "0.1 M"}}}
    assert F.arm_base(none_arm, None)[0] == F.arm_base(cys_arm, None)[0] == frozenset({"carbonyl"})


def test_final_ph_slope_refuses_a_run_with_a_narrow_final_spread():
    arms = [{"value": 10.0, "censor": "none"}, {"value": 20.0, "censor": "none"}, {"value": 5.0, "censor": "none"}]
    # final pH 3.2, 3.4, 5.0: the longest monotone run is 3.2 -> 3.4 (0.2 unit): no slope
    assert F.final_ph_slope({}, arms, [3.2, 3.4, 5.0]) is None
    ok = F.final_ph_slope({}, arms[:2], [5.5, 6.8])
    assert ok is not None and ok[0] > 0
