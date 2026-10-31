"""The two-level measurement model (plan step B7a, 2026-10-31; results/validation/kinetic_core_dpo_partition_prereg.md
secs. 2, 4, 7): species series nested in a figure share the figure's level and deviate from it by a response factor,
Sigma_G = D + tauG^2 11' + sum_s tauW^2 1_s 1_s'. Checked against brute-force linear algebra, against the one-level
model in its limits, and on the case it exists for (a ratio error the shared level cannot absorb)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from src.kinetic_core import fitting

INF = float("inf")


def _rows(n_groups=3, species=("TDP", "DPO"), per=(4, 3), method="semi_quantitative", seed=0, paper="p"):
    """Outer groups paper|figure, each with one series per species (the level groups, keyed by series)."""
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_groups):
        for sp, n in zip(species, per):
            for _ in range(n):
                rows.append({"group": f"{paper}{g}|fig{g}", "series": f"{paper}{g}|fig{g}|{sp}",
                             "sigma": float(rng.uniform(0.03, 0.2)), "method": method})
    return rows


def _model(tau_g=INF, tau_w=None, method="semi_quantitative", grouped=True, **kw):
    scale = dict(fitting.SCALE_TAU)
    between = {}
    if method == "sida":
        between = {"sida": tau_g} if tau_g else {}
    else:
        scale[method] = tau_g
    return fitting.MeasurementModel(
        series_of=lambda row: row["series"], sigma_of=lambda row: row["sigma"], method_of=lambda row: row["method"],
        scale_tau=scale, between_tau=between,
        group_of=(lambda row: row["group"]) if grouped else None,
        within_tau={method: tau_w} if tau_w is not None else {}, **kw)


def _blocks(rows):
    out = {}
    for i, row in enumerate(rows):
        out.setdefault(row["group"], []).append(i)
    return out


def _sigma_g(rows, idx, tau_g, tau_w, include_outer=True):
    s = np.array([rows[i]["sigma"] for i in idx])
    cov = np.diag(s ** 2)
    labels = [rows[i]["series"] for i in idx]
    for key in set(labels):
        m = np.array([lab == key for lab in labels], dtype=float)
        cov += tau_w ** 2 * np.outer(m, m)
    if include_outer:
        cov += tau_g ** 2 * np.ones((len(idx), len(idx)))
    return cov


# --- 1. no grouping: bit-identical --------------------------------------------------------------------------------

def test_without_a_grouping_every_output_is_bit_identical():
    rng = np.random.default_rng(1)
    rows = (_rows(method="semi_quantitative") + _rows(method="sida", seed=2, paper="q")
            + [{"group": "x", "series": "cal", "sigma": 0.1, "method": "calibrated"} for _ in range(4)])
    for row in rows[-4:]:
        row["series"] = "cal"
    r = rng.normal(size=len(rows)) * 0.3
    r[5] = np.nan
    plain = fitting.MeasurementModel(series_of=lambda row: row["series"], sigma_of=lambda row: row["sigma"],
                                     method_of=lambda row: row["method"], between_tau={"sida": 0.3})
    variants = [
        fitting.MeasurementModel(series_of=plain.series_of, sigma_of=plain.sigma_of, method_of=plain.method_of,
                                 between_tau={"sida": 0.3}, within_tau={"semi_quantitative": 0.2}),     # no group_of
        fitting.MeasurementModel(series_of=plain.series_of, sigma_of=plain.sigma_of, method_of=plain.method_of,
                                 between_tau={"sida": 0.3}, group_of=lambda row: row["group"]),       # no within_tau
        fitting.MeasurementModel(series_of=plain.series_of, sigma_of=plain.sigma_of, method_of=plain.method_of,
                                 between_tau={"sida": 0.3}, group_of=lambda row: row["group"],
                                 within_tau={"relative": 0.2}),                                     # method unused
    ]
    clean = np.where(np.isfinite(r), r, 0.0)
    for v in variants:
        assert np.array_equal(plain.whiten(r, rows), v.whiten(r, rows), equal_nan=True)
        assert plain.offsets(r, rows) == v.offsets(r, rows)
        assert list(plain.offsets(r, rows)) == list(v.offsets(r, rows))
        assert plain.nuisance_dof(rows) == v.nuisance_dof(rows)
        assert plain.log_density(clean, rows) == v.log_density(clean, rows)
        assert plain.coverage(r, rows) == v.coverage(r, rows)
        assert fitting.estimate_series_noise(plain, clean, rows) == fitting.estimate_series_noise(v, clean, rows)
        assert fitting.estimate_between_tau(plain, clean, rows) == fitting.estimate_between_tau(v, clean, rows)
        assert v.group_offsets(r, rows) == {}


# --- 2. one level group per outer group, tauW -> 0: the one-level model -------------------------------------------

@pytest.mark.parametrize("tau_g,method", [(INF, "semi_quantitative"), (0.3, "semi_quantitative"), (0.3, "sida"),
                                          (0.0, "sida")])
def test_one_level_group_per_outer_group_with_a_vanishing_within_term_is_the_one_level_model(tau_g, method):
    rows = _rows(n_groups=4, species=("TDP",), per=(5,), method=method)
    r = np.random.default_rng(3).normal(size=len(rows)) * 0.3 + 0.2
    one = _model(tau_g, None, method, grouped=False)
    two = _model(tau_g, 1e-9, method)
    assert np.allclose(one.whiten(r, rows), two.whiten(r, rows), atol=1e-6, rtol=0)
    o1, o2 = one.offsets(r, rows), two.offsets(r, rows)
    assert o1.keys() == o2.keys() and all(abs(o1[k] - o2[k]) < 1e-6 for k in o1)
    assert one.log_density(r, rows) == pytest.approx(two.log_density(r, rows), abs=1e-6)
    assert one.nuisance_dof(rows) == two.nuisance_dof(rows)
    assert one.coverage(r, rows) == two.coverage(r, rows)


# --- 3. brute force: the quadratic form, the whitening's covariance, the log density -----------------------------

@pytest.mark.parametrize("tau_g", [0.0, 0.4, 2.0])
def test_finite_outer_tau_is_the_exact_gaussian(tau_g):
    from scipy.stats import multivariate_normal

    rows = _rows(n_groups=3, species=("TDP", "DPO", "FUR"), per=(3, 2, 4), method="sida")
    r = np.random.default_rng(4).normal(size=len(rows)) * 0.4
    tau_w = 0.25
    mm = _model(tau_g, tau_w, "sida")
    z = mm.whiten(r, rows)
    quad = logpdf = 0.0
    for idx in _blocks(rows).values():
        cov = _sigma_g(rows, idx, tau_g, tau_w)
        quad += float(r[idx] @ np.linalg.solve(cov, r[idx]))
        logpdf += float(multivariate_normal(np.zeros(len(idx)), cov).logpdf(r[idx]))
        # the whitening is linear and Cov z = I exactly
        w = np.array([mm.whiten(np.eye(len(rows))[k], rows)[idx] for k in idx]).T
        assert np.allclose(w @ cov @ w.T, np.eye(len(idx)), atol=1e-10)
    assert float(z @ z) == pytest.approx(quad, rel=1e-10)
    assert mm.log_density(r, rows) == pytest.approx(logpdf, rel=1e-10)
    assert mm.nuisance_dof(rows) == 0


def test_free_outer_level_is_projected_out_in_the_within_metric():
    rows = _rows(n_groups=3, species=("TDP", "DPO", "FUR"), per=(3, 2, 4))
    r = np.random.default_rng(5).normal(size=len(rows)) * 0.4
    tau_w = 0.25
    mm = _model(INF, tau_w)
    z = mm.whiten(r, rows)
    quad = logd = 0.0
    for idx in _blocks(rows).values():
        cov0 = _sigma_g(rows, idx, 0.0, tau_w, include_outer=False)
        one = np.ones(len(idx))
        a, b = np.linalg.solve(cov0, r[idx]), np.linalg.solve(cov0, one)
        q = float(r[idx] @ a) - float(one @ a) ** 2 / float(one @ b)
        quad += q
        logd += -0.5 * (q + float(np.linalg.slogdet(cov0)[1]) + math.log(float(one @ b))
                        + (len(idx) - 1) * math.log(2.0 * math.pi))
        # the whitened vector: Cov z = the projector off the free level, and a shift of the whole group is invisible
        w = np.array([mm.whiten(np.eye(len(rows))[k], rows)[idx] for k in idx]).T
        proj = w @ cov0 @ w.T
        assert np.allclose(proj @ proj, proj, atol=1e-10) and np.trace(proj) == pytest.approx(len(idx) - 1)
        shifted = r.copy()
        shifted[idx] += 0.9
        assert np.allclose(mm.whiten(shifted, rows), z, atol=1e-12)
    assert float(z @ z) == pytest.approx(quad, rel=1e-10)
    assert mm.log_density(r, rows) == pytest.approx(logd, rel=1e-10)
    assert mm.nuisance_dof(rows) == 3   # one per free OUTER group, not per species series


# --- 4. offsets: the BLUP from Henderson's mixed-model equations --------------------------------------------------

@pytest.mark.parametrize("tau_g", [INF, 0.5, 0.0])
def test_offsets_are_the_posterior_means_of_the_mixed_model(tau_g):
    rows = _rows(n_groups=2, species=("TDP", "DPO"), per=(4, 3), method="sida" if tau_g == 0.0 else "semi_quantitative")
    r = np.random.default_rng(6).normal(size=len(rows)) * 0.4 + 0.3
    r[2] = np.nan   # a failed row: left out of the posterior
    tau_w = 0.2
    method = rows[0]["method"]
    mm = _model(tau_g, tau_w, method)
    got, outer = mm.offsets(r, rows), mm.group_offsets(r, rows)
    for g, idx in _blocks(rows).items():
        ok = [i for i in idx if np.isfinite(r[i])]
        keys = sorted({rows[i]["series"] for i in idx})
        cols = [[rows[i]["series"] == k for i in ok] for k in keys]
        prec = [1.0 / tau_w ** 2] * len(keys)
        if tau_g != 0.0:   # an outer level with a flat (tau inf) or N(0, tauG^2) prior; tauG 0: no outer column
            cols, prec = [np.ones(len(ok))] + cols, [0.0 if math.isinf(tau_g) else 1.0 / tau_g ** 2] + prec
        x = np.column_stack(cols).astype(float)
        dinv = np.diag([1.0 / rows[i]["sigma"] ** 2 for i in ok])
        theta = np.linalg.solve(x.T @ dinv @ x + np.diag(prec), x.T @ dinv @ r[ok])
        if tau_g == 0.0:
            theta = np.concatenate([[0.0], theta])
        assert outer[g] == pytest.approx(theta[0], abs=1e-9)
        for j, k in enumerate(keys):
            assert got[k] == pytest.approx(theta[0] + theta[1 + j], abs=1e-9)


def test_a_failed_row_spoils_its_outer_group_only():
    rows = _rows(n_groups=2)
    r = np.random.default_rng(7).normal(size=len(rows)) * 0.3
    r[1] = np.nan
    z = _model(INF, 0.2).whiten(r, rows)
    first = _blocks(rows)["p0|fig0"]
    assert np.all(np.isnan(z[first])) and np.all(np.isfinite(np.delete(z, first)))


# --- 5. the purpose: a shared level cannot absorb a ratio error ---------------------------------------------------

def test_the_shared_figure_level_cannot_absorb_a_ratio_error():
    """TDP read 10x high and DPO 10x low in one figure: one free scale per species hides it completely; one shared
    figure level with a 0.2 dex response-factor term leaves (2 dex)^2 / (2 x 0.2^2) ~ 50 in the cost."""
    rows = _rows(n_groups=1, species=("TDP", "DPO"), per=(6, 6))
    for row in rows:
        row["sigma"] = 0.05
    noise = np.random.default_rng(8).normal(size=12) * 0.01
    r = np.array([1.0 if row["series"].endswith("TDP") else -1.0 for row in rows]) + noise
    one = _model(INF, None, grouped=False)
    two = _model(INF, 0.2)
    cost_one = float(np.sum(one.whiten(r, rows) ** 2))
    cost_two = float(np.sum(two.whiten(r, rows) ** 2))
    assert cost_one < 1.0
    assert cost_two > 40.0 and cost_two > 50 * cost_one
    # a COMMON error of the figure (both species 10x high) is still a free level and costs nothing
    common = np.ones(12) + noise
    assert float(np.sum(two.whiten(common, rows) ** 2)) < 1.0
    assert two.nuisance_dof(rows) == 1 and one.nuisance_dof(rows) == 2


# --- 6. the noise loop keeps the grouping -------------------------------------------------------------------------

def test_fit_with_series_noise_keeps_the_two_level_fields_across_rounds():
    rng = np.random.default_rng(9)
    times = np.linspace(10.0, 200.0, 8)
    rows = []
    for g in range(2):
        for sp, factor in (("TDP", 3.0), ("DPO", 3.0 * 10 ** 0.1)):
            for t in times:
                y = factor * 10.0 * math.exp(-0.01 * t) * 10.0 ** rng.normal(0, 0.05)
                rows.append({"group": f"p{g}", "series": f"p{g}|{sp}", "sigma": 0.005, "method": "semi_quantitative",
                             "t": float(t), "y": y})
    rows = tuple(rows)

    def raw(theta, rows_):
        a, lk = theta
        return np.array([a - (10.0 ** lk) * r["t"] / math.log(10.0) - math.log10(r["y"]) for r in rows_])

    seen = []

    def make(mm):
        seen.append(mm)
        return fitting.Problem(keys=("log10_A", "log10_k"), kinds=("log10k", "log10k"), lower=np.array([-2.0, -5.0]),
                               upper=np.array([3.0, 0.0]), rows=rows,
                               residuals=lambda th, rr: mm.whiten(raw(th, rr), rr), nuisance_dof=mm.nuisance_dof)

    start = _model(INF, 0.2)
    f_, final = fitting.fit_with_series_noise(make, start, raw, np.array([0.5, -2.5]), budget=400, rounds=3, tol=-1.0)
    assert len(seen) == 3
    for mm in seen[1:] + [final]:
        assert mm.group_of is start.group_of and dict(mm.within_tau) == {"semi_quantitative": 0.2}
    assert final.extra_noise and f_.theta[1] == pytest.approx(-2.0, abs=0.05)


# --- validation and the between-paper estimator ------------------------------------------------------------------

def test_an_outer_group_cannot_mix_method_classes_or_split_a_level_group():
    rows = _rows(n_groups=1)
    rows[0] = dict(rows[0], method="sida")
    rows[0]["series"] = "lone"
    with pytest.raises(ValueError, match="mixes method classes"):
        _model(INF, 0.2).whiten(np.zeros(len(rows)), rows)
    rows = _rows(n_groups=2)
    rows[0] = dict(rows[0], group="p1|fig1")   # one TDP row of figure 0 filed under figure 1
    with pytest.raises(ValueError, match="spans outer groups"):
        _model(INF, 0.2).whiten(np.zeros(len(rows)), rows)
    for bad in (0.0, INF, -0.1):
        with pytest.raises(ValueError, match="finite and > 0"):
            _model(INF, bad)


def test_the_between_paper_estimate_is_the_outer_tau_under_two_levels():
    rng = np.random.default_rng(10)
    rows = _rows(n_groups=300, species=("TDP", "DPO"), per=(3, 3), method="sida")
    mu = {f"p{g}|fig{g}": rng.normal(0.0, 0.4) for g in range(300)}
    dev = {row["series"]: rng.normal(0.0, 0.2) for row in rows}
    r = np.array([mu[row["group"]] + dev[row["series"]] + rng.normal(0.0, row["sigma"]) for row in rows])
    two = _model(0.0, 0.2, "sida")
    one = _model(0.0, None, "sida", grouped=False)
    assert fitting.estimate_between_tau(two, r, rows)["sida"] == pytest.approx(0.4, rel=0.15)
    # the one-level estimate over species series carries both terms, sqrt(0.4^2 + 0.2^2) ~ 0.45
    assert fitting.estimate_between_tau(one, r, rows)["sida"] == pytest.approx(math.hypot(0.4, 0.2), rel=0.15)


def test_ninety_percent_intervals_cover_about_ninety_percent_under_two_levels():
    rng = np.random.default_rng(11)
    for tau_g, method in ((0.3, "sida"), (INF, "semi_quantitative")):
        rows = _rows(n_groups=200, species=("TDP", "DPO"), per=(3, 2), method=method)
        mu = {f"p{g}|fig{g}": rng.normal(0.0, 0.3 if math.isfinite(tau_g) else 2.0) for g in range(200)}
        dev = {row["series"]: rng.normal(0.0, 0.2) for row in rows}
        r = np.array([mu[row["group"]] + dev[row["series"]] + rng.normal(0.0, row["sigma"]) for row in rows])
        inside, scored = _model(tau_g, 0.2, method).coverage(r, rows, level=0.90)
        assert scored == len(rows) and 0.87 <= inside / scored <= 0.93
