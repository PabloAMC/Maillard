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


def test_parallel_jacobians_change_nothing(monkeypatch):
    """Strategy M2: the Jacobian's points evaluated by forked workers give the same fit, bit for bit, with the
    same evaluation count, and the same Laplace (one coordinate pinned at its upper bound: the one-sided branch)."""
    p = _problem(upper=(3.0, -2.0, 135.0))
    start = np.array([0.5, -2.5, 80.0])
    runs = {}
    for workers in ("1", "3"):
        monkeypatch.setenv("MAILLARD_WORKERS", workers)
        f = fitting.fit(p, start, budget=400)
        runs[workers] = (f, fitting.laplace(p, f.theta))
    (f1, l1), (f3, l3) = runs["1"], runs["3"]
    assert np.array_equal(f1.theta, f3.theta) and f1.evaluations == f3.evaluations
    assert np.array_equal(f1.residuals, f3.residuals)
    assert l1["sigma"] == l3["sigma"] and l1["chi2_red"] == l3["chi2_red"]


# ---------------------------------------------------------------------------------------------------------------------
# The measurement model (M6) and the diagnostics (M7)
# ---------------------------------------------------------------------------------------------------------------------

def _mm(tau_semi=float("inf")):
    return fitting.MeasurementModel(series_of=lambda row: row["series"], sigma_of=lambda row: row["sigma"],
                                    method_of=lambda row: row["method"],
                                    scale_tau={**fitting.SCALE_TAU, "semi_quantitative": tau_semi})


def _series_rows(n=6, method="semi_quantitative", name="s", sigmas=None):
    sigmas = sigmas if sigmas is not None else [0.05 + 0.02 * i for i in range(n)]
    return [{"series": name, "sigma": s, "method": method} for s in sigmas]


@pytest.mark.parametrize("tau", [0.3, 2.0])
def test_whitening_is_the_exact_marginal_likelihood(tau):
    """z.z equals r' (D + tau^2 11')^-1 r for a series with a shared, integrated-out offset."""
    rows = _series_rows()
    sig = np.array([r["sigma"] for r in rows])
    r = np.random.default_rng(1).normal(size=len(rows)) * 0.1
    z = _mm(tau).whiten(r, rows)
    cov = np.diag(sig ** 2) + tau ** 2 * np.ones((len(rows), len(rows)))
    assert z @ z == pytest.approx(float(r @ np.linalg.solve(cov, r)), rel=1e-10)


def test_whitened_series_have_identity_covariance():
    rng = np.random.default_rng(3)
    rows = _series_rows(n=4)
    sig = np.array([r["sigma"] for r in rows])
    tau = 0.5
    draws = np.array([_mm(tau).whiten(rng.normal(0.0, tau) + rng.normal(size=4) * sig, rows) for _ in range(40000)])
    assert np.allclose(np.cov(draws.T), np.eye(4), atol=0.03)


def test_a_scale_free_series_sees_only_its_shape():
    rows = _series_rows() + _series_rows(name="abs", method="sida")
    r = np.random.default_rng(2).normal(size=len(rows)) * 0.1
    shifted = r.copy()
    shifted[:6] += 0.7  # an unknown response factor of 10^0.7 on the semi-quantitative series only
    mm = _mm()
    assert np.allclose(mm.whiten(r, rows), mm.whiten(shifted, rows))
    assert mm.nuisance_dof(rows) == 1
    shifted[6:] += 0.7  # the same shift on the absolute series is seen
    assert not np.allclose(mm.whiten(r, rows)[6:], mm.whiten(shifted, rows)[6:])


def test_series_offsets_are_the_posterior_mean_of_the_scale():
    rows = _series_rows() + _series_rows(name="abs", method="sida")
    sig = np.array([r["sigma"] for r in rows[:6]])
    r = np.random.default_rng(4).normal(size=len(rows)) * 0.1 + 0.7
    weighted_mean = float(np.sum(r[:6] / sig ** 2) / np.sum(1 / sig ** 2))
    free = _mm().offsets(r, rows)
    assert free["s"] == pytest.approx(weighted_mean) and free["abs"] == 0.0
    shrunk = _mm(0.3).offsets(r, rows)["s"]
    assert 0.0 < shrunk < weighted_mean
    # the residual a reader judges a scale-free series by has no weighted mean left
    assert float(np.sum((r[:6] - free["s"]) / sig ** 2)) == pytest.approx(0.0, abs=1e-9)


def test_a_series_cannot_mix_method_classes():
    rows = _series_rows(n=2) + [{"series": "s", "sigma": 0.1, "method": "sida"}]
    with pytest.raises(ValueError, match="mixes method classes"):
        _mm().whiten(np.zeros(3), rows)


def test_a_fit_recovers_the_rate_from_a_series_with_an_unknown_scale():
    """Two series of one decay: an absolute one and a semi-quantitative one read 5x high. With the measurement model
    the rate is recovered; treating the semi-quantitative series as absolute biases it."""
    times = np.linspace(10.0, 200.0, 8)
    true_a, true_lk = 1.0, -2.0
    rows = []
    for name, method, factor in (("abs", "sida", 1.0), ("semi", "semi_quantitative", 5.0)):
        for t in times:
            y = factor * 10.0 ** true_a * math.exp(-(10.0 ** true_lk) * t)
            rows.append({"series": name, "sigma": 0.03, "method": method, "t": float(t), "y": y})
    rows = tuple(rows)

    def raw(theta, rows_):
        a, lk = theta
        return np.array([a - (10.0 ** lk) * row["t"] / math.log(10.0) - math.log10(row["y"]) for row in rows_])

    def problem(mm):
        return fitting.Problem(keys=("log10_A", "log10_k"), kinds=("log10k", "log10k"), lower=np.array([-2.0, -5.0]),
                               upper=np.array([3.0, 0.0]), rows=rows,
                               residuals=lambda theta, rows_: mm.whiten(raw(theta, rows_), rows_),
                               nuisance_dof=mm.nuisance_dof)

    honest = fitting.fit(problem(_mm()), np.array([0.5, -2.5]), budget=600)
    naive = fitting.fit(problem(_mm(tau_semi=0.0)), np.array([0.5, -2.5]), budget=600)
    assert honest.theta[1] == pytest.approx(true_lk, abs=1e-3)
    assert honest.theta[0] == pytest.approx(true_a, abs=1e-3)
    assert abs(naive.theta[0] - true_a) > 0.2  # the absolute level is dragged toward the 5x series
    lap = fitting.laplace(problem(_mm()), honest.theta)
    assert np.isfinite(lap["chi2_red"])  # 16 rows - 2 constants - 1 scale-free series


def test_duplicates_are_counted_once():
    rows = [{"ds": "zhai2023_table1", "v": 1}, {"ds": "kang2026_s4", "v": 1}, {"ds": "other", "v": 2}]
    kept, dropped = fitting.count_once(rows, lambda r: r["ds"], {"kang2026_s4": "zhai2023_table1"})
    assert [r["ds"] for r in kept] == ["zhai2023_table1", "other"]
    assert [r["ds"] for r in dropped] == ["kang2026_s4"]
    kept, dropped = fitting.count_once(rows[1:], lambda r: r["ds"], {"kang2026_s4": "zhai2023_table1"})
    assert [r["ds"] for r in kept] == ["kang2026_s4", "other"] and not dropped  # the copy stands in for its original


def test_diagnostics_flag_a_biased_laboratory_and_a_trend():
    rows = tuple({"lab": lab, "T": t} for lab in ("A", "B") for t in (90.0, 100.0, 110.0, 120.0, 130.0, 140.0))

    def residuals(theta, rows_):
        return np.array([(3.0 if row["lab"] == "B" else 0.0) + 0.05 * (row["T"] - 115.0) for row in rows_])

    p = fitting.Problem(keys=("x",), kinds=("log10k",), lower=np.array([-1.0]), upper=np.array([1.0]), rows=rows,
                        residuals=residuals)
    report = fitting.diagnostics(p, np.array([0.0]), {"lab": lambda r: r["lab"], "T": lambda r: r["T"]})
    assert any(f.startswith("lab=B") for f in report["flags"])
    assert not any(f.startswith("lab=A") for f in report["flags"])
    assert report["trends"]["T"]["rho"] > 0
    assert "| B | 6 |" in fitting.diagnostics_markdown(report)


def test_series_noise_is_estimated_where_the_stated_error_is_too_small():
    """A series whose stated error (0.001 dex, a vector read) is far below its real scatter (0.1 dex) gets an extra
    noise near 0.1; a series whose stated error is honest gets none or little; a two-point series gets none."""
    rng = np.random.default_rng(5)
    rows = ([{"series": "tiny", "sigma": 0.001, "method": "calibrated"} for _ in range(40)]
            + [{"series": "honest", "sigma": 0.1, "method": "calibrated"} for _ in range(40)]
            + [{"series": "pair", "sigma": 0.001, "method": "calibrated"} for _ in range(2)])
    r = np.concatenate([rng.normal(0, 0.1, 40), rng.normal(0, 0.1, 40), rng.normal(0, 0.1, 2)])
    mm = _mm()
    noise = fitting.estimate_series_noise(mm, r, rows)
    assert noise["tiny"] == pytest.approx(0.1, rel=0.25)
    assert noise.get("honest", 0.0) < 0.05
    assert "pair" not in noise
    with_noise = fitting.MeasurementModel(series_of=mm.series_of, sigma_of=mm.sigma_of, method_of=mm.method_of,
                                          extra_noise=noise)
    z = with_noise.whiten(r, rows)
    assert float(np.mean(z[:40] ** 2)) == pytest.approx(1.0, rel=0.35)   # the tiny series is no longer 10^4 x over-weighted


def test_a_scale_free_series_noise_ignores_its_offset():
    rows = [{"series": "semi", "sigma": 0.01, "method": "semi_quantitative"} for _ in range(30)]
    r = 0.7 + np.random.default_rng(6).normal(0, 0.05, 30)   # a response factor of 10^0.7 plus 0.05 dex of scatter
    noise = fitting.estimate_series_noise(_mm(), r, rows)
    assert noise["semi"] == pytest.approx(0.05, rel=0.3)


def test_fit_with_series_noise_recovers_the_rate_when_one_series_states_a_tiny_error():
    """Two series of one decay, both with 0.05 dex real scatter; one states 0.0005 dex (a vector read). Without the
    noise term that series dominates; with it the rate is recovered and both series carry their real scatter."""
    rng = np.random.default_rng(9)
    times = np.linspace(10.0, 200.0, 12)
    rows = []
    for name, stated in (("vector", 0.0005), ("table", 0.05)):
        for t in times:
            y = 10.0 * math.exp(-0.01 * t) * 10.0 ** rng.normal(0, 0.05)
            rows.append({"series": name, "sigma": stated, "method": "calibrated", "t": float(t), "y": y})
    rows = tuple(rows)

    def raw(theta, rows_):
        a, lk = theta
        return np.array([a - (10.0 ** lk) * r["t"] / math.log(10.0) - math.log10(r["y"]) for r in rows_])

    def make(mm):
        return fitting.Problem(keys=("log10_A", "log10_k"), kinds=("log10k", "log10k"), lower=np.array([-2.0, -5.0]),
                               upper=np.array([3.0, 0.0]), rows=rows,
                               residuals=lambda th, rr: mm.whiten(raw(th, rr), rr), nuisance_dof=mm.nuisance_dof)

    f_, mm = fitting.fit_with_series_noise(make, _mm(), raw, np.array([0.5, -2.5]), budget=600)
    assert f_.theta[1] == pytest.approx(-2.0, abs=0.05)
    assert mm.extra_noise["vector"] == pytest.approx(0.05, rel=0.5)
    # the honest series' TOTAL error stays near its true 0.05 (12 points: the moment estimate fluctuates ~40 %)
    assert math.sqrt(0.05 ** 2 + mm.extra_noise.get("table", 0.0) ** 2) < 0.08


def test_evaluate_many_is_used_for_batches_and_must_match_the_single_evaluation():
    # sub-model 1 (2026-10-29): a flattened batch evaluation feeds the Jacobian; non-finite values become the
    # failed residual exactly as in the single-point path
    calls = []

    def residuals(theta, rows):
        return np.array([theta[0] - 1.0, theta[1] + 2.0])

    def many(points, rows):
        calls.append(len(points))
        out = [residuals(x, rows) for x in points]
        out[-1] = np.array([np.nan, 0.0])
        return out

    problem = fitting.Problem(keys=("a", "b"), kinds=("ph_slope", "log10_offset"), lower=np.array([0.0, -3.0]),
                              upper=np.array([2.0, 3.0]), rows=(1, 2), residuals=residuals, evaluate_many=many)
    values = fitting._evaluate_all(problem, [np.array([0.5, 0.0]), np.array([1.0, 1.0])], None)
    assert calls == [2] and np.allclose(values[0], [-0.5, 2.0]) and values[1][0] == problem.failed_residual
    assert fitting._step("ph_slope") == 0.05 and fitting._threshold("log10_offset") == 0.5


# ---------------------------------------------------------------------------------------------------------------------
# The acceptance protocol's measurement pieces (2026-10-30, STRATEGY sec. 9): the between-paper error of an absolute
# method class, the log predictive density and the interval coverage
# ---------------------------------------------------------------------------------------------------------------------

def _abs_rows(n_papers=40, per=5, sigma=0.05):
    return [{"series": f"p{g}|pot{i}", "paper": f"p{g}", "sigma": sigma, "method": "sida"}
            for g in range(n_papers) for i in range(per)]


def _between_model(tau=0.0):
    return fitting.MeasurementModel(series_of=lambda row: row["series"], sigma_of=lambda row: row["sigma"],
                                    method_of=lambda row: row["method"], level_of=lambda row: row["paper"],
                                    between_tau={"sida": tau} if tau else {})


def test_the_between_paper_error_is_recovered_from_paper_offsets():
    rng = np.random.default_rng(7)
    rows = _abs_rows()
    offsets = {f"p{g}": rng.normal(0.0, 0.4) for g in range(40)}
    r = np.array([offsets[row["paper"]] + rng.normal(0.0, row["sigma"]) for row in rows])
    tau = fitting.estimate_between_tau(_between_model(), r, rows)["sida"]
    assert tau == pytest.approx(0.4, rel=0.25)


def test_no_between_term_leaves_absolute_rows_exactly_as_before():
    rows = _abs_rows(n_papers=3)
    r = np.random.default_rng(1).normal(size=len(rows)) * 0.1
    plain = fitting.MeasurementModel(series_of=lambda row: row["series"], sigma_of=lambda row: row["sigma"],
                                     method_of=lambda row: row["method"])
    assert np.allclose(plain.whiten(r, rows), _between_model().whiten(r, rows))


def test_log_density_is_the_exact_gaussian_with_a_shared_paper_offset():
    from scipy.stats import multivariate_normal

    rows = _abs_rows(n_papers=1, per=4)
    r = np.random.default_rng(3).normal(size=4) * 0.2
    tau = 0.3
    cov = np.diag([row["sigma"] ** 2 for row in rows]) + tau ** 2 * np.ones((4, 4))
    assert _between_model(tau).log_density(r, rows) == pytest.approx(
        float(multivariate_normal(np.zeros(4), cov).logpdf(r)), rel=1e-9)


def test_ninety_percent_intervals_cover_about_ninety_percent():
    rng = np.random.default_rng(11)
    rows = _abs_rows(n_papers=200, per=3)
    offsets = {f"p{g}": rng.normal(0.0, 0.3) for g in range(200)}
    r = np.array([offsets[row["paper"]] + rng.normal(0.0, row["sigma"]) for row in rows])
    inside, scored = _between_model(0.3).coverage(r, rows, level=0.90)
    assert scored == len(rows) and 0.86 <= inside / scored <= 0.94


def test_compact_json_round_trips_and_writes_a_row_per_line():
    import json

    from src.provenance import dumps_compact

    payload = {"rows": [{"pot": f"p{i}", "measured": 1.5, "predicted": None, "flags": [1, 2]} for i in range(50)],
               "meta": {"n": 50}}
    text = dumps_compact(payload)
    assert json.loads(text) == payload
    assert len(text.splitlines()) <= 50 + 6
