"""Toy-input tests of the lever-table model (scripts/generators/build_lever_table.py): ladders, censoring, the
exact Gaussian fit, the network of categorical levels and the reading rule. No engine, no keyed data."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("build_lever_table",
                                               ROOT / "scripts" / "generators" / "build_lever_table.py")
LT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(LT)


def _arm(x, value, censor="none", limit=None, ph_basis="initial", temp=140, time=60):
    return {"lever_value": x, "value": value, "censor": censor, "limit": limit,
            "conditions": {"ph": x, "ph_basis": ph_basis, "temp_C": temp, "time_min": time}}


def test_two_arm_ph_contrast_is_the_normalised_log_ratio():
    o = LT.numeric_observation({"lever_class": "pH"}, [_arm(5.0, 10.0), _arm(7.0, 1.0)])
    assert o["y"] == pytest.approx(-0.5)            # -1 dex over 2 pH units
    assert o["band"] == "pH 6-7.5, initial or label" and not o["non_monotone"]   # midpoint 6.0 (PREREG sec. 4)


def test_ladder_is_one_slope_and_a_peak_is_flagged():
    up = LT.numeric_observation({"lever_class": "pH"}, [_arm(4.0, 1.0), _arm(5.0, 10.0), _arm(6.0, 100.0)])
    assert up["y"] == pytest.approx(1.0) and up["arms_used"] == 3
    peak = LT.numeric_observation({"lever_class": "pH"},
                                  [_arm(4.0, 1.0), _arm(5.0, 10.0), _arm(6.0, 100.0), _arm(7.0, 10.0)])
    assert peak["non_monotone"] and peak["arms_used"] == 3 and peak["y"] == pytest.approx(1.0)


def test_censoring_rules():
    arms = [_arm(5.0, 8.0), _arm(6.0, None, censor="trace"), _arm(7.0, None, censor="nd_no_limit")]
    LT._trace_fill(arms)
    assert arms[1]["censor"] == "below_value" and arms[1]["limit"] == 8.0
    assert LT.arm_log_and_var(arms[2], False)[0] is None
    ly, var, _ = LT.arm_log_and_var({"value": 10.0, "sd": 1.0, "n": 4}, False)
    assert ly == pytest.approx(1.0) and var == pytest.approx((0.1 / math.log(10)) ** 2 / 4)


def test_numeric_fit_recovers_a_common_effect_across_labs():
    rng = np.random.default_rng(3)
    rows = [{"y": -0.4 + rng.normal(0, 0.05), "v": 0.003, "lab": f"L{i}", "source": f"s{i}"} for i in range(6)]
    s = LT.numeric_fit(rows, LT.PRIOR_SCALE)
    assert s["ci90"][0] < -0.4 < s["ci90"][1] and s["ci90"][1] < 0
    assert s["p_same_sign_new_lab"] > 0.9
    assert LT.reading(6, s) == "answered by the literature"
    assert LT.reading(2, s) == "too few laboratories"


def test_heterogeneous_labs_make_the_direction_unclear():
    rows = [{"y": y, "v": 0.003, "lab": f"L{i}", "source": f"s{i}"} for i, y in enumerate((0.8, -0.7, 0.6, -0.5))]
    s = LT.numeric_fit(rows, LT.PRIOR_SCALE)
    assert LT.reading(4, s) == "direction unclear"


def test_categorical_network_connects_through_a_shared_level():
    # study 1: glucose vs ribose; study 2: ribose vs xylose (no glucose arm) -> xylose reaches glucose via ribose
    rows = [{"pairs": [{"level": "ribose", "against": "glucose", "y": 1.0, "v": 0.01}], "shared_var": 0.005,
             "lab": "A", "source": "s1"},
            {"pairs": [{"level": "xylose", "against": "ribose", "y": 0.2, "v": 0.01}], "shared_var": 0.005,
             "lab": "B", "source": "s2"},
            {"pairs": [{"level": "fructose", "against": "rhamnose", "y": 0.3, "v": 0.01}], "shared_var": 0.005,
             "lab": "C", "source": "s3"}]
    f = LT.categorical_fit(rows, "glucose", LT.PRIOR_SCALE)
    assert f["levels"]["xylose"]["estimate"] == pytest.approx(1.2, abs=0.15)
    assert f["levels"]["fructose"] == {"reading": "not connected"}


@pytest.mark.parametrize("scale", [0.25, 0.5, 1.0])
def test_the_grid_prior_is_the_half_normal_it_declares(scale):
    """GRID is log-spaced: each point's weight must be the half-normal density times its cell width, else the grid
    applies ~ half-normal / s and piles prior mass near zero (sparse cells' intervals came out too narrow)."""
    w = np.exp(LT._grid_log_prior(scale))
    w /= w.sum()
    assert (w * LT.GRID).sum() == pytest.approx(scale * math.sqrt(2.0 / math.pi), rel=0.03)
    true_below = math.erf(0.1 / scale / math.sqrt(2.0))           # P(|N(0, scale)| < 0.1)
    assert w[LT.GRID < 0.1].sum() == pytest.approx(true_below, abs=0.03)


# ------------------------------------------------------------- sensitivity readings and fragility (RULES_V2_SPEC)
def test_a_flat_ladder_that_wiggles_has_no_turning_point_and_a_real_peak_has_one():
    flat = LT.numeric_observation({"lever_class": "pH"}, [_arm(4.0, 3.9), _arm(5.0, 3.6), _arm(6.0, 4.1), _arm(7.0, 3.7)])
    peak = LT.numeric_observation({"lever_class": "pH"}, [_arm(4.0, 1.0), _arm(5.0, 10.0), _arm(6.0, 1.0)])
    assert flat["non_monotone"] and not flat["turning_point"]
    assert peak["non_monotone"] and peak["turning_point"]


def test_whole_ladder_slope_of_a_symmetric_peak_is_zero_while_the_longest_run_is_not():
    peak = LT.numeric_observation({"lever_class": "pH"}, [_arm(4.0, 1.0), _arm(5.0, 10.0), _arm(6.0, 1.0)])
    assert abs(peak["y"]) == pytest.approx(1.0)
    assert peak["alt"]["whole_ladder"][0] == pytest.approx(0.0, abs=1e-9)


def test_without_censored_arms_a_two_arm_contrast_has_no_alternative_and_a_ladder_keeps_its_measured_arms():
    two = LT.numeric_observation({"lever_class": "pH"}, [_arm(5.0, 10.0), _arm(7.0, None, censor="nd_with_limit", limit=0.1)])
    assert two["y"] == pytest.approx(-1.0) and two["censored_arms"] == 1 and two["alt"]["no_censored"] is None
    lad = LT.numeric_observation({"lever_class": "pH"}, [_arm(5.0, 100.0), _arm(6.0, 10.0),
                                                          _arm(7.0, None, censor="nd_with_limit", limit=0.001)])
    assert lad["alt"]["no_censored"][0] == pytest.approx(-1.0) and lad["y"] < -1.5


def test_buffer_state_reads_ph_adjusted_without_buffer_as_unbuffered():
    arms = lambda b: [{"conditions": {"buffer": b}}]
    assert LT.buffer_state(arms("none (pH adjusted with 3 mol/L NaOH)")) == "unbuffered"
    assert LT.buffer_state(arms("0.2 M phosphate")) == "buffered"
    assert LT.buffer_state(arms("unstated")) == "unstated"


def test_a_cell_answered_only_through_detection_limits_is_flagged_fragile():
    def row(lab, censored):
        arms = [_arm(5.0, 100.0), _arm(6.0, 50.0),
                _arm(7.0, None, censor="nd_with_limit", limit=0.01) if censored else _arm(7.0, 45.0)]
        o = LT.numeric_observation({"lever_class": "pH"}, arms)
        return {**o, "lab": lab, "source": lab}
    rows = [row("a", True), row("b", True), row("c", True)]
    fr = LT.numeric_fragility(rows, "pH", "answered by the literature")
    assert fr["contrasts_with_a_censored_arm"] == 3
    assert "without censored arms" in fr["variants"]
    assert abs(fr["variants"]["without censored arms"]["estimate"]) < abs(sum(r["y"] for r in rows) / 3)
