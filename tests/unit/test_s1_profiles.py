"""scripts/generators/build_s1_profiles.py: the profile likelihood on a toy quadratic with a known interval."""
from __future__ import annotations

import math

import numpy as np
import pytest

from scripts.generators import build_s1_profiles as S


def _toy(correlated: bool = True, lower=(-10.0, -10.0), upper=(10.0, 10.0)) -> S.FitSpec:
    """r = A x - b, linear: S(x) = |A x - b|^2 is an exact quadratic. Its profile in x0 (x1 re-optimised) is
    Delta(psi) = (psi - x0_hat)^2 / Var(x0), Var = [(A^T A)^-1]_00, so the 95 % LR interval is
    x0_hat +/- sqrt(3.841 Var)."""
    a = np.array([[1.0, 0.8 if correlated else 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, -1.0]]) * 4.0
    b = np.array([0.3, -0.2, 0.1, 0.4]) * 4.0
    return S.FitSpec("toy", ("log10_a", "log10_b"), ("log10k", "log10k"), np.array(lower), np.array(upper),
                     np.array([0.0, 0.0]), lambda x: a @ np.asarray(x, dtype=float) - b, n_rows=4,
                     central_polish=1, meta={"A": a, "b": b})


def _exact(spec: S.FitSpec):
    a, b = spec.meta["A"], spec.meta["b"]
    x_hat = np.linalg.lstsq(a, b, rcond=None)[0]
    var = np.linalg.inv(a.T @ a)[0, 0]
    return x_hat, var


def test_profile_of_a_quadratic_recovers_its_known_interval():
    spec = _toy()
    x_hat, var = _exact(spec)
    rec = S.profile_fit(spec, [0], {0: 0.2 * math.sqrt(var)})
    assert rec["coordinates"]["log10_a"]["x_hat"] == pytest.approx(x_hat[0], abs=1e-4)
    raw = rec["coordinates"]["log10_a"]["intervals"]["raw"]
    half = math.sqrt(S.CHI2_1_95 * var)
    assert raw["closed_lower"] and raw["closed_upper"]
    assert raw["lower"] == pytest.approx(x_hat[0] - half, abs=2e-3)
    assert raw["upper"] == pytest.approx(x_hat[0] + half, abs=2e-3)
    # the profile is not the slice: with x1 held at its optimum the interval would be narrower
    slice_half = math.sqrt(S.CHI2_1_95 / (spec.meta["A"][:, 0] @ spec.meta["A"][:, 0]))
    assert half > 1.1 * slice_half
    # every walked point's Delta is the exact profile
    for p in rec["coordinates"]["log10_a"]["low_side"] + rec["coordinates"]["log10_a"]["high_side"]:
        assert p["delta"] == pytest.approx((p["value"] - x_hat[0]) ** 2 / var, rel=1e-3, abs=1e-5)


def test_a_bound_inside_the_interval_leaves_it_open_on_that_side():
    spec0 = _toy()
    x_hat, var = _exact(spec0)
    half = math.sqrt(S.CHI2_1_95 * var)
    spec = _toy(lower=(x_hat[0] - 0.5 * half, -10.0))
    rec = S.profile_fit(spec, [0], {0: 0.2 * math.sqrt(var)})
    raw = rec["coordinates"]["log10_a"]["intervals"]["raw"]
    assert not raw["closed_lower"] and raw["closed_upper"]
    assert raw["lower"] == pytest.approx(x_hat[0] - 0.5 * half)


def test_scaled_threshold_inflates_by_the_reduced_chi_square_only_above_one():
    t = S.thresholds(0.5, 4)
    assert t["scaled"] == pytest.approx(S.CHI2_1_95) and t["F"] is None
    t = S.thresholds(2.0, 4)
    assert t["scaled"] == pytest.approx(2.0 * S.CHI2_1_95)
    assert t["F"] == pytest.approx(2.0 * 7.708647, rel=1e-5)


def test_the_28_coordinates_are_the_s1_rows_missing_a_profile():
    rows = S.tier_rows()
    assert len(rows) == 28
    assert sum(S.GROUP_OF[c] == "sulfur" for c in rows) == 11


def test_no_band_on_record_is_never_a_search_bound():
    row = {"band": {"source": "x", "lo": -4.0, "hi": 0.0}, "band_is_search_bound": True,
           "identification": {"ci95_halfwidth": 0.3}}
    band, note = S.s4_band(row)
    assert band is None and "never a band" in note
    row = {"band": {"source": "envelope prior p", "sigma": 0.3 / 1.96, "distribution": "normal_log10"},
           "band_is_search_bound": False, "identification": {"ci95_halfwidth": 0.3}}
    band, note = S.s4_band(row)
    assert band is None and "Laplace sigma" in note
