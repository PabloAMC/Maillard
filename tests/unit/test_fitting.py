"""src/kinetic_core/fitting.py on a synthetic problem with a known answer, and its constants against B8's."""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import fitting

ROOT = Path(__file__).resolve().parents[2]
SIGMA = 0.02  # dex, the synthetic measurement error
TRUE = {"A": (1.0, -2.0), "B": (1.0, -2.0)}  # (log10 amplitude, log10 rate) per group


def _rows(truth=TRUE, seed=7):
    rng = np.random.default_rng(seed)
    rows = []
    for group, (a, lk) in truth.items():
        for t in np.linspace(5.0, 200.0, 12):
            y = 10.0 ** a * math.exp(-(10.0 ** lk) * t) * 10.0 ** (SIGMA * rng.standard_normal())
            rows.append((group, float(t), float(y)))
    return tuple(rows)


def _residuals(theta, rows):
    a, lk, _ignored = theta
    # a small, rough "integrator noise" term: the reason the finite-difference step must not be tiny
    rough = 2e-4 * math.sin(3e5 * (a + lk))
    return np.array([(a - (10.0 ** lk) * t / math.log(10.0) - math.log10(y) + rough) / SIGMA for _, t, y in rows])


def _problem(rows=None, upper=(3.0, 0.0, 135.0)):
    return fitting.Problem(keys=("log10_A", "log10_k", "ea_ignored"), kinds=("log10k", "log10k", "Ea"),
                           lower=np.array([-3.0, -5.0, 40.0]), upper=np.array(upper),
                           rows=_rows() if rows is None else rows, residuals=_residuals, group_of=lambda r: r[0])


def test_constants_match_b8():
    laplace_src = (ROOT / "scripts/generators/generate_kinetic_core_b8_laplace.py").read_text()
    fit_src = (ROOT / "scripts/generators/generate_kinetic_core_b8_fit.py").read_text()
    assert eval(re.search(r"^THRESHOLDS = (\{.*\})", laplace_src, re.M).group(1)) == fitting.THRESHOLDS
    assert eval(re.search(r"^STEPS = (\{.*\})", laplace_src, re.M).group(1)) == fitting.STEPS
    assert float(re.search(r"^NULL_WEIGHT_TOL = (\S+)", laplace_src, re.M).group(1)) == fitting.NULL_WEIGHT_TOL
    assert "diff_step=3e-2, ftol=1e-6, xtol=1e-8" in fit_src
    s = fitting.LEAST_SQUARES_SETTINGS
    assert (s["ftol"], s["xtol"], s["method"], s["x_scale"]) == (1e-6, 1e-8, "trf", "jac")
    assert "diff_step" not in s  # replaced by absolute per-kind steps (see the module docstring)


def test_scipy_relative_step_fails_near_zero_where_absolute_steps_do_not():
    """The bug class found 2026-10-29: B8's relative diff_step stalls a fit whose coordinate starts at 0."""
    from scipy.optimize import least_squares

    p = _problem()
    start = np.array([0.0, -3.5, 80.0])
    scipy_relative = least_squares(p.evaluate, start, bounds=(p.lower, p.upper), method="trf", x_scale="jac",
                                   diff_step=fitting.B8_RELATIVE_DIFF_STEP, ftol=1e-6, xtol=1e-8)
    ours = fitting.fit(p, start)
    assert scipy_relative.cost > 100 * ours.cost


def test_fit_recovers_truth_on_a_rough_objective_from_a_far_start():
    p = _problem()
    f = fitting.fit(p, np.array([0.0, -3.5, 80.0]))
    lap = fitting.laplace(p, f.theta)
    assert abs(f.theta[0] - 1.0) < 3 * lap["sigma"]["log10_A"]
    assert abs(f.theta[1] + 2.0) < 3 * lap["sigma"]["log10_k"]
    assert 0.5 < lap["chi2_red"] < 2.0


def test_laplace_sigma_matches_linear_regression():
    p = _problem()
    f = fitting.fit(p, np.array([1.0, -2.0, 80.0]))
    lap = fitting.laplace(p, f.theta)
    t = np.array([r[1] for r in p.rows])
    x = np.column_stack([np.ones_like(t), -t / math.log(10.0)])  # log10 y = a - b t / ln10, b = 10**lk
    cov_ab = lap["chi2_red"] * SIGMA ** 2 * np.linalg.inv(x.T @ x)  # in parameter units (dex)
    sigma_a = math.sqrt(cov_ab[0, 0])
    assert lap["sigma"]["log10_A"] == pytest.approx(sigma_a, rel=0.15)


def test_unconstrained_direction_gets_infinite_sigma_and_is_unidentified():
    p = _problem()
    f = fitting.fit(p, np.array([1.0, -2.0, 80.0]))
    lap = fitting.laplace(p, f.theta)
    assert math.isinf(lap["sigma"]["ea_ignored"]) and not lap["identified"]["ea_ignored"]
    assert lap["identified"]["log10_A"] and lap["identified"]["log10_k"]


def test_bound_rule_is_b8s():
    p = _problem()
    flags = fitting.identification(p, theta=[3.0 - 1e-4, -2.0, 60.0], sigma=[0.7, 5.0, 10.0])
    assert flags["at_bound"]["log10_A"] and not flags["identified"]["log10_A"]  # on bound, sigma 0.7 > band 6 / 10
    assert not flags["identified"]["log10_k"]  # sigma above the log10k threshold
    flags = fitting.identification(p, theta=[3.0 - 1e-4, -2.0, 40.0 + 1e-3], sigma=[0.05, 0.1, 50.0])
    assert flags["at_bound"]["log10_A"] and flags["identified"]["log10_A"]  # on bound but tight
    assert flags["at_bound"]["ea_ignored"] and not flags["identified"]["ea_ignored"]  # 0.003 kJ above the floor


def test_a_fit_pinned_on_a_bound_is_flagged():
    p = _problem(upper=(3.0, -2.2, 135.0))  # the true log10 k (-2.0) lies outside the band
    f = fitting.fit(p, np.array([1.0, -3.0, 80.0]))
    assert fitting.laplace(p, f.theta)["at_bound"]["log10_k"]


def test_two_starts_converge():
    out = fitting.two_starts(_problem(), np.array([0.0, -3.5, 60.0]), np.array([2.0, -1.0, 120.0]))
    assert out["converged"] and out["unidentified"] == ["ea_ignored"]


def test_leave_one_group_out_finds_tension_only_where_it_exists():
    p = _problem()
    full = fitting.fit(p, np.array([1.0, -2.0, 80.0]))
    out = fitting.leave_one_group_out(p, full.theta, fitting.laplace(p, full.theta))
    assert all(not g["beyond_2_sigma"] and not g["predictive"]["tension"] for g in out.values())
    tense = _problem(rows=_rows({"A": (1.0, -2.0), "B": (1.0, -1.6)}))
    full = fitting.fit(tense, np.array([1.0, -2.0, 80.0]))
    out = fitting.leave_one_group_out(tense, full.theta, fitting.laplace(tense, full.theta))
    assert out["A"]["predictive"]["tension"] and out["B"]["predictive"]["tension"]
    # the parameter shift in the full fit's sigmas can miss it: the tension inflates the full chi-square
    assert out["A"]["predictive"]["rms_held_out"] > 3 * out["A"]["predictive"]["rms_in_sample"]


def test_failed_evaluation_is_a_bad_point_not_a_crash():
    def boom(theta, rows):
        raise RuntimeError("integration failed")

    p = fitting.Problem(keys=("x",), kinds=("log10k",), lower=np.array([0.0]), upper=np.array([1.0]),
                        rows=(1, 2, 3), residuals=boom)
    assert np.all(p.evaluate(np.array([0.5])) == p.failed_residual)


def test_cost_drop_threshold():
    assert fitting.cost_drop_is_significant(0.080, 0.0412)  # 0.5 x 3.84 x 0.0412 = 0.0791
    assert not fitting.cost_drop_is_significant(0.078, 0.0412)
