"""Toy-input tests of step 0 of the decision-relevance pre-registration (scripts/generators/build_decision_relevance_step0.py):
the stage-0c estimator (pair-graph residuals, the floor, G_min, LOSO) and the stage-0a oracle. No engine, no data."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "build_decision_relevance_step0", ROOT / "scripts" / "generators" / "build_decision_relevance_step0.py")
S0 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S0)


def _toy_pairs(n_studies, sigma_e, rng, arms=4):
    """Arms of each study get residuals r_i = u + e_i, e_i ~ N(0, sigma_e^2); pairs carry e_ab = r_a - r_b."""
    pairs = []
    for s in range(n_studies):
        u = rng.normal(0, 1.0)
        r = u + rng.normal(0, sigma_e, size=arms)
        for i in range(arms):
            for j in range(i + 1, arms):
                pairs.append({"study": f"s{s}", "compound": "MFT", "lane": "sulfur", "a": f"s{s}a{i}",
                              "b": f"s{s}a{j}", "e": float(r[i] - r[j]), "oos": True})
    return pairs


def test_pair_graph_recovers_centred_arm_residuals():
    r = {"a": 0.3, "b": -0.1, "c": 0.5}
    pairs = [(x, y, r[x] - r[y]) for x, y in (("a", "b"), ("a", "c"), ("b", "c"))]
    got = S0._lstsq_arm_residuals(pairs)
    mean = sum(r.values()) / 3
    for k in r:
        assert got[k] == pytest.approx(r[k] - mean, abs=1e-12)


def test_sigma_r_recovers_the_ratio_sd_and_ignores_the_study_offset():
    rng = np.random.default_rng(1)
    sigma_e = 0.3                                     # per-arm independent part: ratio SD = sqrt 2 x 0.3
    pairs = _toy_pairs(40, sigma_e, rng, arms=5)
    groups = S0.groups_from_pairs(pairs, {})
    for g in groups:
        g["measurement_var"] = [0.0] * g["k"]         # no measurement part in the toy
    est = S0.sigma_r_estimate(groups)
    assert est["status"] == "ESTIMATED"
    assert est["sigma_r"] == pytest.approx(math.sqrt(2) * sigma_e, rel=0.12)
    assert est["t_interval_90"][0] <= est["sigma_r"] <= est["t_interval_90"][1]
    assert est["loso_min"] <= est["sigma_r"] <= est["loso_max"]


def test_floor_is_the_measurement_ratio_sd_never_zero():
    rng = np.random.default_rng(2)
    pairs = _toy_pairs(6, 1e-6, rng)                  # no discrepancy at all
    groups = S0.groups_from_pairs(pairs, {})
    est = S0.sigma_r_estimate(groups)
    assert est["floored"] is True
    assert est["sigma_r"] == pytest.approx(math.sqrt(2.0 * est["sigma_m2"]))
    assert est["sigma_r"] > 0


def test_fewer_than_g_min_studies_is_undecided():
    rng = np.random.default_rng(3)
    groups = S0.groups_from_pairs(_toy_pairs(S0.G_MIN - 1, 0.3, rng), {})
    assert S0.sigma_r_estimate(groups)["status"].startswith("UNDECIDED")
    s0c = {"primary": {"sulfur": S0.sigma_r_estimate(groups)}, "pooled": S0.sigma_r_estimate(groups)}
    est, rule = S0.sigma_r_for("sulfur", s0c)
    assert est is None and rule.startswith("UNDECIDED")


def test_lane_rule_falls_back_to_the_pooled_sigma_r():
    rng = np.random.default_rng(4)
    big = S0.sigma_r_estimate(S0.groups_from_pairs(_toy_pairs(8, 0.3, rng), {}))
    small = S0.sigma_r_estimate(S0.groups_from_pairs(_toy_pairs(2, 0.3, rng), {}))
    est, rule = S0.sigma_r_for("lipid", {"primary": {"lipid": small, "sulfur": big}, "pooled": big})
    assert est is big and rule.startswith("pooled")


def test_oracle_ceiling_large_levers_eligible_small_levers_omega_bound():
    rng = np.random.default_rng(5)
    big = np.concatenate([rng.normal(2.5, 0.2, 500), rng.normal(-2.5, 0.2, 500)])
    small = rng.normal(0.0, 0.3, 1000)
    assert S0.oracle(big, 0.59, np.random.default_rng(0))["eligible"] is True
    assert S0.oracle(small, 0.59, np.random.default_rng(0))["eligible"] is False
    # a perfect design with no discrepancy left shows every worthwhile lever correctly
    o = S0.oracle(small, 1e-6, np.random.default_rng(0))
    assert o["P_wrong"] == pytest.approx(0.0, abs=1e-9)


def test_sigma_star_is_where_the_bar_is_last_met():
    L = np.full(400, 1.0)
    s = S0.sigma_star(L, seed=1)
    assert s is not None and 0.2 < s < 1.0 / S0.Z90 + 1e-9


def test_lever_ratio_of_ratios():
    rec = {"values": {"A": {"MFT": 1.0, "FFT": 1.0}, "B": {"MFT": 10.0, "FFT": 1.0}}}
    assert S0.lever(rec, "A", "B", "MFT/FFT") == pytest.approx(1.0)
    assert S0.lever(rec, "A", "B", "MFT") == pytest.approx(1.0)
    assert S0.lever({"values": {"A": {"MFT": 0.0}, "B": {"MFT": 1.0}}}, "A", "B", "MFT") is None
