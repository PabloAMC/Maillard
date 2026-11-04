"""H0's scoring step (results/validation/kinetic_core_between_lab_tau_h0_prereg.md): the estimator's four identity
tests of sec. 2, the vectorised predictive against fitting.MeasurementModel, the St quadrature against its Gaussian
limit, and the blinding precondition (--score refuses to run without the ten pair records, each with its freeze
commit). Synthetic rows only."""
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))

import fit_thiol_lolo as lolo  # noqa: E402
import fit_thiol_lolo_pairs as pair_runner  # noqa: E402
import score_thiol_lolo_nested_tau as h0  # noqa: E402
from src.kinetic_core import fitting  # noqa: E402


def _row(paper, species, method, i, lab="lab_a", sigma=0.05, table="thiol"):
    return SimpleNamespace(ledger_id=f"{paper}_t{i % 2}", species=species, method=method, table=table, t_min=float(i),
                           sigma=sigma, series=f"{paper}|{species}|{method}", laboratory=lab, pot=f"{paper}_p{i}",
                           value=1.0, unit="mmol/L")


def _abs_rows(lab, n_groups, rows_per_group, method="sida", tag=""):
    rows = []
    for g in range(n_groups):
        for k in range(rows_per_group if not isinstance(rows_per_group, list) else rows_per_group[g]):
            rows.append(_row(f"p{tag}{lab[-1]}{g}{method[0]}", f"S{g}", method, k, lab=lab,
                             sigma=0.03 + 0.01 * ((g + k) % 4)))
    return rows


def _model_and_s2(rows, noise=0.1):
    model = lolo._model({fitting_series(r): noise for r in rows}, level=lolo.level_of)
    s2 = np.array([model.sigma(r) ** 2 for r in rows])
    return model, s2


def fitting_series(r):
    return lolo.fj.series_of(r)


def _residuals(rows, rng, tau=0.7):
    lay = h0.Layout.of_rows(rows)
    level = tau * rng.standard_normal(lay.G)
    return level[lay.gid] + 0.1 * rng.standard_normal(len(rows))


# ---- sec. 2 identity tests -------------------------------------------------------------------------------------

def test_identity_1_one_lab_one_method_equals_the_shipped_estimator():
    rng = np.random.default_rng(1)
    rows = _abs_rows("lab_a", 6, [1, 2, 3, 4, 1, 5])
    model, s2 = _model_and_s2(rows)
    r = _residuals(rows, rng)
    shipped = fitting.estimate_between_tau(model, r, rows)["sida"]
    est = h0.estimate([h0.Block("lab_a", h0.Layout.of_rows(rows), s2, r[None, :])], 1, omega2_fixed=0.0)
    assert shipped > 0.0
    assert abs(float(est["ABS"]["tau"][0]) - shipped) <= 1e-12 * max(1.0, shipped)
    assert est["ABS"]["status"] == "single-lab"   # one laboratory: flagged, the estimate itself unchanged


def test_identity_2_pooled_dbar_is_the_group_count_weighted_mean_of_the_methods():
    rng = np.random.default_rng(2)
    rows = _abs_rows("lab_a", 4, [2, 3, 1, 4], "sida") + _abs_rows("lab_a", 3, [3, 1, 2], "calibrated", tag="c")
    _, s2 = _model_and_s2(rows)
    r = _residuals(rows, rng)
    lay = h0.Layout.of_rows(rows)
    dbar = {}
    for method in ("sida", "calibrated"):
        d = []
        for g in range(lay.G):
            if lay.method_g[g] != method:
                continue
            idx = lay.gid == g
            w = 1.0 / s2[idx]
            m = float(np.sum(w * r[idx]) / np.sum(w))
            d.append(m * m - 1.0 / float(np.sum(w)))
        dbar[method] = (len(d), float(np.mean(d)))
    expected = sum(n * v for n, v in dbar.values()) / sum(n for n, _ in dbar.values())
    est = h0.estimate([h0.Block("lab_a", lay, s2, r[None, :])], 1, omega2_fixed=0.0)
    assert abs(float(est["ABS"]["tau2_raw"][0]) - expected) <= 1e-12 * max(1.0, abs(expected))


def test_identity_3_within_group_moment_has_expectation_c_omega2():
    rng = np.random.default_rng(3)
    rows = _abs_rows("lab_a", 3, [2, 3, 5])
    lay = h0.Layout.of_rows(rows)
    s2 = np.array([0.02, 0.05, 0.01, 0.03, 0.08, 0.02, 0.04, 0.06, 0.01, 0.05])
    omega2 = 0.09
    R = 100_000
    level = 3.0 * rng.standard_normal((R, lay.G))   # an arbitrary level: it must cancel
    r = level[:, lay.gid] + np.sqrt(s2 + omega2)[None, :] * rng.standard_normal((R, lay.n))
    _, e, c, sel = h0.within_moments(lay, s2, r)
    assert sel.all()
    mean, se = e.mean(axis=0), e.std(axis=0, ddof=1) / math.sqrt(R)
    assert np.all(np.abs(mean - c * omega2) <= 3.0 * se), (mean - c * omega2) / se


def test_identity_4_equal_lab_sizes_make_lab_weights_equal_group_weights():
    rng = np.random.default_rng(4)
    blocks, pooled = [], []
    for lab in ("lab_a", "lab_b", "lab_c"):
        rows = _abs_rows(lab, 4, [1, 3, 2, 2])
        _, s2 = _model_and_s2(rows)
        r = _residuals(rows, rng, tau=0.9)
        blocks.append(h0.Block(lab, h0.Layout.of_rows(rows), s2, r[None, :]))
        pooled.append((rows, r))
    rows_all = [x for rows, _ in pooled for x in rows]
    r_all = np.concatenate([r for _, r in pooled])
    model, _ = _model_and_s2(rows_all)
    shipped = fitting.estimate_between_tau(model, r_all, rows_all)["sida"]
    est = h0.estimate(blocks, 1, omega2_fixed=0.0)
    assert est["ABS"]["status"] == "normal" and shipped > 0.0
    assert abs(float(est["ABS"]["tau"][0]) - shipped) <= 1e-12 * max(1.0, shipped)


# ---- the vectorised predictive equals MeasurementModel; St reduces to the Gaussian ---------------------------

def test_vectorised_arm_matches_measurement_model_and_row_intervals():
    rng = np.random.default_rng(5)
    rows = (_abs_rows("lab_a", 3, [1, 3, 2])
            + [_row("figA", "DPO", "semi_quantitative", k, table="deoxyosone_submodel") for k in range(4)]
            + [_row("zhou2023", "MFT", "relative", k) for k in range(5)]
            + [_row("zhou2024", "FFT", "relative", 0)])
    _, s2 = _model_and_s2(rows)
    r = rng.standard_normal(len(rows))
    omega2, taus = 0.04, {"sida": 0.5, "calibrated": 0.5, "semi_figure": 0.45}
    arm = h0.gaussian_arm(h0.Layout.of_rows(rows), s2, r[None, :], omega2, taus)
    noise = {lolo.fj.series_of(x): math.sqrt(v - x.sigma ** 2 + omega2) for x, v in zip(rows, s2)}
    pm = lolo._model(noise, {"sida": 0.5, "calibrated": 0.5}, lolo.level_of)
    import dataclasses
    pm = dataclasses.replace(pm, scale_tau={**lolo.sm3.SCALE_TAU, "semi_figure": 0.45})
    assert abs(float(arm["lpd"][0]) - pm.log_density(r, rows)) <= 1e-9 * abs(pm.log_density(r, rows))
    for p, c, v in zip(lolo.row_intervals(pm, r, rows), arm["centred"][0], arm["var"][0]):
        if p["var"] is None:
            assert not np.isfinite(v)
        else:
            assert abs(p["var"] - v) <= 1e-12 and abs(p["centred"] - c) <= 1e-12


def test_st_reduces_to_the_gaussian_at_large_df():
    rng = np.random.default_rng(6)
    r, st2 = rng.standard_normal(4) * 0.5, np.array([0.01, 0.04, 0.02, 0.03])
    tau = 0.6
    W, B, A = np.sum(1 / st2), np.sum(r / st2), np.sum(r * r / st2)
    gauss = -0.5 * (A - tau ** 2 * B * B / (1 + tau ** 2 * W) + np.sum(np.log(st2)) + math.log1p(tau ** 2 * W)
                    + 4 * math.log(2 * math.pi))
    assert abs(h0.st_group_logdensity(r, st2, tau, 1e7) - gauss) < 1e-5
    from scipy.stats import norm
    x, sd = np.array([0.3, 1.0, 2.0]), np.array([0.2, 0.5, 1.0])
    for scale in (0.1, 0.8):
        got = h0.st_cdf_abs(x, sd, np.full(3, scale), 1e7, nodes=200)
        want = 2 * norm.cdf(x / np.sqrt(sd ** 2 + scale ** 2)) - 1
        assert np.allclose(got, want, atol=1e-6)


# ---- blinding: --score refuses without the ten pair records ------------------------------------------------

def _pair_record(a, b, commit="abc1234"):
    return {"artifact": f"{pair_runner.STEM}_{a}__{b}", "labs": [a, b], "freeze_commit": commit, "complete": True}


def test_score_refuses_without_any_pair_record(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "refuses" in str(exc.value.code)


def test_score_refuses_with_nine_of_ten_pair_records(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    for a, b in pair_runner.pairs()[:9]:
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b)))
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "missing" in str(exc.value.code)


def test_score_refuses_when_a_pair_record_has_no_freeze_commit(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    ps = pair_runner.pairs()
    for i, (a, b) in enumerate(ps):
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b, commit="" if i == 3 else "abc1234")))
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "freeze_commit" in str(exc.value.code)


def test_the_ten_complete_pair_records_pass_the_precondition(tmp_path, monkeypatch):
    for a, b in pair_runner.pairs():
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b)))
    assert len(h0.check_pair_records(tmp_path)) == 10


def test_the_sign_flip_and_the_exact_bootstrap_are_exact():
    counts, w = h0.multisets(5)
    assert counts.shape == (126, 5) and abs(w.sum() - 1.0) < 1e-12
    assert h0.sign_patterns(5).shape == (32, 5)
    d = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    res = h0.d3_decision(d[None, :], -d[None, :])
    assert res["p_better"][0] == 1 / 32 and res["lpd_pass"][0] and not res["sign_trigger"][0]
