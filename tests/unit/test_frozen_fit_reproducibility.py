"""scripts/generators/build_frozen_fit_reproducibility.py: the audit on a toy linear problem with a known optimum."""
from __future__ import annotations

import numpy as np
import pytest

from scripts.generators import build_frozen_fit_reproducibility as F

A = np.array([[1.0, 0.8], [0.0, 1.0], [1.0, 1.0], [0.5, -1.0]]) * 4.0
B = np.array([0.3, -0.2, 0.1, 0.4]) * 4.0


def _resid(shift: float = 0.0):
    """r = A x - (b + shift): linear, so S is an exact quadratic and one Gauss-Newton step reaches its optimum."""
    return lambda x: A @ np.asarray(x, dtype=float) - (B + shift)


def _spec(x0, shift=0.0, s_recorded=None, rows_then=None, lower=(-10.0, -10.0), upper=(10.0, 10.0)) -> F.FitSpec:
    res = _resid(shift)
    ids = [f"row{i}" for i in range(4)]
    if s_recorded is None:
        r = res(x0)
        s_recorded = float(r @ r)
    return F.FitSpec("toy", ("log10_a", "log10_b"), ("log10k", "log10k"), np.array(lower), np.array(upper),
                     np.asarray(x0, dtype=float), res, {i: {"target": float(b), "sigma": 1.0} for i, b in zip(ids, B + shift)},
                     rows_then, s_recorded, 4, 1.0, coordinate_keys={"a.log10_k_ref": ["log10_a"],
                                                                     "b.log10_k_ref": ["log10_b"]})


def _optimum(shift: float = 0.0) -> np.ndarray:
    return np.linalg.lstsq(A, B + shift, rcond=None)[0]


def test_frozen_optimum_reproduces():
    x_hat = _optimum()
    rec = F.audit_fit(_spec(x_hat), ["a.log10_k_ref", "b.log10_k_ref"])
    assert rec["verdict"]["class"] == "REPRODUCES"
    assert rec["near_optimality"]["one_step_decrease"] == pytest.approx(0.0, abs=1e-8)
    assert rec["near_optimality"]["gradient_norm"] == pytest.approx(0.0, abs=1e-8)


def test_one_gauss_newton_step_reaches_the_optimum_of_a_linear_problem():
    x0 = np.array([0.5, -0.5])
    r0 = _resid()(x0)
    jac, _ = F.forward_jacobian(_resid(), x0, r0, np.array([-10.0, -10.0]), np.array([10.0, 10.0]),
                                ("log10k", "log10k"))
    assert jac == pytest.approx(A, abs=1e-9)
    d = F.gauss_newton_step(r0, jac, x0, np.array([-10.0, -10.0]), np.array([10.0, 10.0]))
    assert x0 + d == pytest.approx(_optimum(), abs=1e-9)
    opt = F.near_optimality(_resid(), x0, r0, jac, np.array([-10.0, -10.0]), np.array([10.0, 10.0]), 1.0,
                            {"a": [0]})
    s_min = float(np.sum((A @ _optimum() - B) ** 2))
    assert opt["one_step_decrease"] == pytest.approx(float(r0 @ r0) - s_min, rel=1e-9)
    # the quadratic marginal of coordinate a: (x0 - x_hat)^2 / [(A'A)^-1]_00
    var = np.linalg.inv(A.T @ A)[0, 0]
    assert opt["coordinates"]["a"]["delta_quadratic"] == pytest.approx((x0[0] - _optimum()[0]) ** 2 / var, rel=1e-6)


def test_small_drift_is_benign_and_large_drift_is_material():
    x_hat = _optimum()
    r_frozen = _resid()(x_hat)
    s_rec = float(r_frozen @ r_frozen)
    # the data moved a little: S changes beyond 1 %, the frozen point is still inside today's 95 % region
    small = F.audit_fit(_spec(x_hat, shift=0.05, s_recorded=s_rec), ["a.log10_k_ref", "b.log10_k_ref"])
    assert not small["verdict"]["reproduced"]
    assert small["verdict"]["class"] == "DRIFTED-BENIGN"
    # the data moved a lot: today's optimum is far from the frozen point
    big = F.audit_fit(_spec(x_hat, shift=np.array([2.0, -2.0, 1.5, 2.0]), s_recorded=s_rec),
                      ["a.log10_k_ref", "b.log10_k_ref"])
    assert big["verdict"]["class"] == "DRIFTED-MATERIAL"
    assert big["verdict"]["coordinates_outside_local_95"]


def test_the_tolerance_does_not_scale_with_todays_misfit():
    """s2 is the freeze's chi2_red: a drift that inflates today's misfit must not widen its own tolerance."""
    x_hat = _optimum()
    spec = _spec(x_hat, shift=np.array([2.0, -2.0, 1.5, 2.0]), s_recorded=1.0)
    rec = F.audit_fit(spec, ["a.log10_k_ref"])
    assert rec["near_optimality"]["s2"] == 1.0
    assert rec["chi2_red_today"] > 1.0


def test_row_changes_classify_withdrawn_corrected_and_model_moves():
    then = {"a": {"r": 1.0, "target": 2.0, "sigma": 0.1}, "b": {"r": 0.5, "target": 3.0, "sigma": 0.1},
            "c": {"r": -0.2, "target": 1.0, "sigma": 0.1}, "d": {"r": 0.1, "target": 1.0, "sigma": 0.1}}
    now = {"b": {"r": 0.9, "target": 3.3, "sigma": 0.1}, "c": {"r": 2.0, "target": 1.0, "sigma": 0.1},
           "d": {"r": 0.1, "target": 1.0, "sigma": 0.1}, "e": {"r": 1.5, "target": 4.0, "sigma": 0.2}}
    rc = F.row_changes(then, now)
    kinds = {q["id"]: q["kind"] for q in rc["all"]}
    assert kinds == {"a": "withdrawn", "b": "corrected", "c": "model change", "d": "unchanged", "e": "added"}
    assert rc["all"][0]["id"] == "c"  # the largest change in r^2 first


def test_bounds_keep_the_step_inside_the_box():
    x0 = np.array([0.0, 0.0])
    r0 = _resid()(x0)
    lower, upper = np.array([-0.01, -10.0]), np.array([0.01, 10.0])
    jac, _ = F.forward_jacobian(_resid(), x0, r0, lower, upper, ("log10k", "log10k"))
    d = F.gauss_newton_step(r0, jac, x0, lower, upper)
    assert lower[0] - 1e-12 <= d[0] <= upper[0] + 1e-12


def test_markdown_renders_and_check_round_trips(tmp_path, monkeypatch):
    rec = F.audit_fit(_spec(_optimum()), ["a.log10_k_ref"])
    payload = {"artifact": "frozen_fit_reproducibility", "fits": {"B20": rec}, "method": F.METHOD}
    payload["summary"] = F.summary(payload)
    md = F.markdown(payload)
    assert "REPRODUCES" in md and "### B20" in md
    monkeypatch.setattr(F, "OUT_JSON", tmp_path / "x.json")
    monkeypatch.setattr(F, "OUT_MD", tmp_path / "x.md")
    monkeypatch.setattr(F, "_inputs", lambda fits: [])
    F.write_payload(payload)
    assert F.main(["--check"]) == 0
    (tmp_path / "x.md").write_text("stale\n", encoding="utf-8")
    assert F.main(["--check"]) == 1


def test_a_failed_quadratic_model_cannot_certify_benign():
    """Strongly curved r: the Gauss-Newton step overshoots, the trust ratio falls below 1/4, and a drift whose
    joint and marginal tests pass on the quadratic model is still not certified near-optimal."""
    def res(x):
        x = np.asarray(x, dtype=float)
        return np.array([x[0] ** 3 - 1.0, 0.1 * x[0]])

    x0 = np.array([0.3])
    r0 = res(x0)
    jac, _ = F.forward_jacobian(res, x0, r0, np.array([-10.0]), np.array([10.0]), ("log10k",))
    opt = F.near_optimality(res, x0, r0, jac, np.array([-10.0]), np.array([10.0]), 100.0, {"a": [0]},
                            backtracks=0)
    rho = F.trust_ratio(opt)
    assert rho is not None and rho < F.MIN_TRUST_RATIO
    v = F.classify(0.5 * float(r0 @ r0), float(r0 @ r0), 2, 2, 0.01, opt, 1)
    assert v["class"] == "DRIFTED-MATERIAL" and not v["local_model_valid"]
