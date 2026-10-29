"""
src/kinetic_core/fitting.py -- ONE TESTED WAY TO FIT, CHECK CONVERGENCE AND READ IDENTIFICATION (2026-10-29).

WHY THIS EXISTS
===============
Eighteen generators call ``scipy.optimize.least_squares`` on their own and forty-nine carry their own
identification logic. On 2026-10-29 the route calibration re-implemented both and got both wrong: it left the
optimiser's finite-difference step at its default (~1e-8, below the integrator's noise at rtol 1e-7), so the
Jacobian was noise and a "converged" fit was not at an optimum; and it tested bounds with an absolute tolerance
of 1e-6, so constants sitting on their bounds were reported identified. Both were settled long ago in B8
(``generate_kinetic_core_b8_fit.py``, ``generate_kinetic_core_b8_laplace.py``). This module carries B8's
settings and rules once, with tests, for every new wave (tasks/data_restructure_plan.md §7, "HOW WE WORK" item 3).

WHAT IT DOES
============
A ``Problem`` is a residual function over rows plus named, bounded coordinates of a declared kind. ``fit`` runs
the codebase's optimiser settings with an evaluation budget; ``two_starts`` fits from two starts and applies the
convergence criterion (costs within 1 %, every identified constant within its Laplace sigma: a status flag is not
convergence); ``laplace`` builds B8's central-difference Jacobian, its eigenbasis covariance (unconstrained
directions get an infinite sigma, not zero) and B8's identification rule; ``leave_one_group_out`` refits without
each group of rows (a source, a laboratory) from the full optimum and reports every constant's shift in sigmas
and every bound it hits; ``cost_drop_is_significant`` is the chi-square rule the envelope and the pre-registrations
use. Nothing here knows chemistry.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from . import parallel

#: B8's optimiser settings (generate_kinetic_core_b8_fit.py, fit_member): trust-region reflective, Jacobian
#: scaling, tolerances. B8 also passes ``diff_step = 3e-2``, which scipy reads as a RELATIVE step: 3 % of the
#: coordinate's value. At a coordinate near zero (a log10 k of -0.02) that step collapses to machine epsilon and the
#: Jacobian column is the integrator's noise; across units it is inconsistent (3 kJ/mol on an Ea of 100, 0.018 dex
#: on a log10 k of -0.6). Found 2026-10-29 by tests/unit/test_fitting.py. This module does not use it: ``fit``
#: supplies its own forward-difference Jacobian with an ABSOLUTE step per coordinate kind (B8's Laplace STEPS).
LEAST_SQUARES_SETTINGS: Dict[str, Any] = {"method": "trf", "x_scale": "jac", "ftol": 1e-6, "xtol": 1e-8}
B8_RELATIVE_DIFF_STEP = 3e-2  # recorded for the comparison test only
#: B8's Laplace (generate_kinetic_core_b8_laplace.py): identification thresholds and central-difference steps by
#: coordinate kind, and the null-space weight above which a sigma is unbounded. tests/unit/test_fitting.py
#: checks these equal B8's.
THRESHOLDS: Dict[str, float] = {"Ea": 60.0, "log10k": 3.0, "acid_yield": 0.5, "pka": 1.5}
STEPS: Dict[str, float] = {"Ea": 2.0, "log10k": 0.05, "acid_yield": 0.01, "pka": 0.1}
#: Kinds B8 did not have (sub-model 1, 2026-10-29; design review S3), kept apart so the table above stays B8's: a pH
#: slope in decades per pH unit (identified within half the width of Martins 2003's measured band, 0.37-0.92) and a
#: log10 offset of a rate constant from its tabulated value.
EXTRA_THRESHOLDS: Dict[str, float] = {"ph_slope": 0.25, "log10_offset": 0.5}
EXTRA_STEPS: Dict[str, float] = {"ph_slope": 0.05, "log10_offset": 0.05}


def _threshold(kind: str) -> float:
    return THRESHOLDS[kind] if kind in THRESHOLDS else EXTRA_THRESHOLDS[kind]


def _step(kind: str) -> float:
    return STEPS[kind] if kind in STEPS else EXTRA_STEPS[kind]


NULL_WEIGHT_TOL = 1e-6
#: On a bound when within this fraction of the band width of it (B8).
BOUND_FRACTION = 1e-3
#: The convergence criterion (route calibration, run log 1): two starts agree in cost to this relative tolerance.
COST_AGREEMENT = 0.01


@dataclass(frozen=True)
class Problem:
    """Named, bounded coordinates and a residual function of (theta, rows) returning one residual per row."""

    keys: Tuple[str, ...]
    kinds: Tuple[str, ...]
    lower: np.ndarray
    upper: np.ndarray
    rows: Tuple[Any, ...]
    residuals: Callable[[np.ndarray, Tuple[Any, ...]], np.ndarray]
    group_of: Callable[[Any], str] = field(default=lambda row: "all")
    #: what a failed evaluation (an integration error) returns per row: a bad point, not a crash of the fit
    failed_residual: float = 25.0
    #: degrees of freedom the residuals have already spent on nuisance parameters for these rows (a
    #: ``MeasurementModel``'s scale-free series: one each); subtracted wherever a reduced chi-square is formed
    nuisance_dof: Callable[[Tuple[Any, ...]], int] = field(default=lambda rows: 0)
    #: optional: residuals at many points in one call (a caller that can flatten points x work items over the
    #: workers; sub-model 1's fit runs every pot of every Jacobian point in one pool). Must equal [residuals(x, rows)].
    evaluate_many: Optional[Callable[[Sequence[np.ndarray], Tuple[Any, ...]], List[np.ndarray]]] = None

    def __post_init__(self):
        assert len(self.keys) == len(self.kinds) == len(self.lower) == len(self.upper), "keys, kinds and bounds must align"
        assert all(k in THRESHOLDS or k in EXTRA_THRESHOLDS for k in self.kinds), f"unknown coordinate kind in {self.kinds}"
        assert np.all(np.asarray(self.lower) < np.asarray(self.upper)), "every band must have positive width"

    def evaluate(self, theta: np.ndarray, rows: Optional[Tuple[Any, ...]] = None) -> np.ndarray:
        rows = self.rows if rows is None else rows
        try:
            r = np.asarray(self.residuals(np.asarray(theta, dtype=float), rows), dtype=float)
        except Exception:
            return np.full(len(rows), self.failed_residual)
        return np.where(np.isfinite(r), r, self.failed_residual)


@dataclass
class Fit:
    theta: np.ndarray
    cost: float
    start: np.ndarray
    start_cost: float
    status: int
    evaluations: int
    budget: int
    seconds: float
    residuals: np.ndarray

    @property
    def budget_exhausted(self) -> bool:
        return self.status == 0

    def as_dict(self, keys: Sequence[str]) -> Dict[str, Any]:
        return {"theta": dict(zip(keys, map(float, self.theta))), "cost": self.cost, "start_cost": self.start_cost,
                "start": dict(zip(keys, map(float, self.start))), "status": self.status,
                "budget_exhausted": self.budget_exhausted, "evaluations": self.evaluations, "budget": self.budget,
                "seconds": self.seconds}


def fit(problem: Problem, start: np.ndarray, budget: int = 1200, rows: Optional[Tuple[Any, ...]] = None,
        upper: Optional[np.ndarray] = None) -> Fit:
    """One least-squares fit with the codebase's settings. ``budget`` counts every evaluation, Jacobian columns
    included (least_squares' max_nfev does not count them)."""
    from scipy.optimize import least_squares

    upper = np.asarray(problem.upper if upper is None else upper, dtype=float)
    lower = np.asarray(problem.lower, dtype=float)
    start = np.clip(np.asarray(start, dtype=float), lower, upper)
    count = {"n": 0}

    def f(theta):
        count["n"] += 1
        return problem.evaluate(theta, rows)

    def f_many(points):
        count["n"] += len(points)  # counted here: a forked worker's increments do not come back
        return _evaluate_all(problem, points, rows)

    def jac(theta):
        return _forward_jacobian(problem, theta, lower, upper, f_many)

    t0 = time.time()
    r0 = f(start)
    result = least_squares(f, start, jac=jac, bounds=(lower, upper),
                           max_nfev=max(1, budget // (len(problem.keys) + 1)), **LEAST_SQUARES_SETTINGS)
    r = f(result.x)
    return Fit(theta=np.asarray(result.x, dtype=float), cost=float(0.5 * r @ r), start=start,
               start_cost=float(0.5 * r0 @ r0), status=int(result.status), evaluations=count["n"], budget=budget,
               seconds=round(time.time() - t0, 1), residuals=r)


def _evaluate_all(problem: Problem, points: Sequence[np.ndarray], rows) -> List[np.ndarray]:
    """``problem.evaluate`` at every point, in parallel processes when ``MAILLARD_WORKERS`` > 1
    (2026-10-29, strategy M2). The points of a finite-difference Jacobian are equal-sized, independent
    objective evaluations: the coarsest grain, so it scales where pot-level parallelism stalls on the
    slowest pot. Bit-identical to the serial loop (tests/unit/test_fitting.py)."""
    if problem.evaluate_many is not None:
        rows_ = problem.rows if rows is None else rows
        try:
            values = problem.evaluate_many([np.asarray(x, dtype=float) for x in points], rows_)
        except Exception:
            values = None
        if values is not None:
            return [np.where(np.isfinite(r), r, problem.failed_residual) for r in map(np.asarray, values)]
    return parallel.map_ordered(lambda x: problem.evaluate(x, rows), points)


def _forward_jacobian(problem: Problem, theta, lower, upper, f_many) -> np.ndarray:
    """Forward differences with B8's absolute step per coordinate kind, stepping inward at an upper bound.
    ``f_many`` evaluates a list of points (the centre first)."""
    theta = np.asarray(theta, dtype=float)
    steps, points = [], [theta]
    for j, kind in enumerate(problem.kinds):
        h = _step(kind) if theta[j] + _step(kind) <= upper[j] else -_step(kind)
        x = theta.copy()
        x[j] += h
        steps.append(h)
        points.append(x)
    values = f_many(points)
    r0 = values[0]
    jac = np.zeros((len(r0), len(theta)))
    for j, (h, r) in enumerate(zip(steps, values[1:])):
        jac[:, j] = (r - r0) / h
    return jac


def laplace(problem: Problem, theta: np.ndarray, rows: Optional[Tuple[Any, ...]] = None,
            upper: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """B8's Laplace at ``theta``: a central-difference Jacobian with a step per coordinate kind (one-sided at a
    bound), JtJ inverted on its eigenbasis, scaled by the reduced chi-square; a coordinate with weight in the
    dropped subspace gets an infinite sigma; B8's identification rule."""
    upper = np.asarray(problem.upper if upper is None else upper, dtype=float)
    lower = np.asarray(problem.lower, dtype=float)
    theta = np.asarray(theta, dtype=float)
    # Every point the differences need, evaluated in one batch (parallel when asked): the centre, then per
    # coordinate the step up and / or down that its bounds allow.
    points, plan = [theta], []
    for j, kind in enumerate(problem.kinds):
        h = _step(kind)
        up, down = theta[j] + h <= upper[j], theta[j] - h >= lower[j]
        xp, xm = theta.copy(), theta.copy()
        xp[j] += h
        xm[j] -= h
        plan.append((h, len(points) if up else None, len(points) + up if (down or not up) else None))
        points += ([xp] if up else []) + ([xm] if (down or not up) else [])
    values = _evaluate_all(problem, points, rows)
    r0 = values[0]
    jac = np.zeros((len(r0), len(theta)))
    for j, (h, i_up, i_down) in enumerate(plan):
        if i_up is not None and i_down is not None:
            jac[:, j] = (values[i_up] - values[i_down]) / (2 * h)
        elif i_up is not None:
            jac[:, j] = (values[i_up] - r0) / h
        else:
            jac[:, j] = (r0 - values[i_down]) / h
    chi2_red = float(r0 @ r0) / max(len(r0) - len(theta) - problem.nuisance_dof(problem.rows if rows is None else rows), 1)
    jtj = jac.T @ jac
    w, v = np.linalg.eigh(0.5 * (jtj + jtj.T))
    keep = w > 1e-12 * float(np.max(w)) if np.max(w) > 0 else np.zeros_like(w, dtype=bool)
    cov = (v[:, keep] / w[keep]) @ v[:, keep].T * chi2_red
    sigma = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    sigma = np.where(np.sum(v[:, ~keep] ** 2, axis=1) > NULL_WEIGHT_TOL, np.inf, sigma)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(sigma, sigma)
    flags = identification(problem, theta, sigma, upper)
    return {"chi2_red": chi2_red, "sigma_common": float(np.sqrt(chi2_red)),
            "sigma": dict(zip(problem.keys, map(float, sigma))), **flags,
            "correlation": np.where(np.isfinite(corr), corr, np.nan).round(3).tolist()}


def identification(problem: Problem, theta, sigma, upper: Optional[np.ndarray] = None) -> Dict[str, Dict[str, bool]]:
    """B8's rule: ON a bound within BOUND_FRACTION of the band width; identified when sigma is within the kind's
    threshold, unless on a bound with sigma above a tenth of the band (a one-sided direction a Gaussian
    misdescribes)."""
    upper = np.asarray(problem.upper if upper is None else upper, dtype=float)
    at_bound, identified = {}, {}
    for key, kind, t, s, lo, hi in zip(problem.keys, problem.kinds, theta, sigma, problem.lower, upper):
        width = hi - lo
        on = (t - lo) <= BOUND_FRACTION * width or (hi - t) <= BOUND_FRACTION * width
        ok = bool(np.isfinite(s) and s <= _threshold(kind))
        at_bound[key], identified[key] = bool(on), bool(ok and not (on and s > 0.1 * width))
    return {"at_bound": at_bound, "identified": identified}


def two_starts(problem: Problem, start_a: np.ndarray, start_b: np.ndarray, budget: int = 1200) -> Dict[str, Any]:
    """Fit from two starts; converged when the costs agree within COST_AGREEMENT and every identified constant
    agrees within its Laplace sigma (at the better optimum). Unidentified constants are listed, not compared."""
    fits = [fit(problem, start_a, budget), fit(problem, start_b, budget)]
    best = min(fits, key=lambda f_: f_.cost)
    lap = laplace(problem, best.theta)
    sig = np.array([lap["sigma"][k] for k in problem.keys])
    ident = np.array([lap["identified"][k] for k in problem.keys])
    gap = np.abs(fits[0].theta - fits[1].theta)
    disagreeing = [k for k, g, s, ok in zip(problem.keys, gap, sig, ident) if ok and g > s]
    cost_ok = abs(fits[0].cost - fits[1].cost) <= COST_AGREEMENT * max(min(fits[0].cost, fits[1].cost), 1e-12)
    return {"fits": fits, "best": best, "laplace": lap, "converged": bool(cost_ok and not disagreeing),
            "cost_agreement": cost_ok, "disagreeing": disagreeing,
            "unidentified": [k for k, ok in zip(problem.keys, ident) if not ok]}


def leave_one_group_out(problem: Problem, theta_full: np.ndarray, laplace_full: Dict[str, Any],
                        budget: int = 1200, groups: Optional[Sequence[str]] = None) -> Dict[str, Dict[str, Any]]:
    """Refit without each group's rows, starting from the full optimum; each constant's shift in the full fit's
    sigmas, the constants moving beyond 2 sigma, and the bounds hit. A bound hit is a loss of identification, not
    a tension, and is reported as such."""
    labels = sorted({problem.group_of(r) for r in problem.rows}) if groups is None else list(groups)
    sig = np.array([laplace_full["sigma"][k] for k in problem.keys])
    out = {}
    for g in labels:
        rows = tuple(r for r in problem.rows if problem.group_of(r) != g)
        held = tuple(r for r in problem.rows if problem.group_of(r) == g)
        f_ = fit(problem, theta_full, budget, rows=rows)
        shift = (f_.theta - np.asarray(theta_full)) / sig
        flags = identification(problem, f_.theta, sig)
        out[g] = {"fit": f_, "n_rows": len(rows),
                  "shift_sigma": dict(zip(problem.keys, map(float, shift))),
                  "beyond_2_sigma": [k for k, s in zip(problem.keys, shift) if np.isfinite(s) and abs(s) > 2.0],
                  "on_bound": [k for k, on in flags["at_bound"].items() if on],
                  "predictive": predictive_check(problem, f_, held, rows)}
    return out


def predictive_check(problem: Problem, f_: Fit, held: Tuple[Any, ...], kept: Tuple[Any, ...],
                     level: float = 0.95) -> Dict[str, Any]:
    """The held-out group's residuals under the model fitted without it, against the fit's own in-sample scatter:
    F = mean square held out / reduced chi-square of the leave-out fit, on (n_held, n_kept - k) degrees of freedom.
    Tension when F exceeds its ``level`` quantile. Unlike a parameter shift in the full fit's sigmas, this does
    not hide a tension that inflates the full fit's chi-square. Points within one series are not independent, so the
    test is anti-conservative; read it with that in mind."""
    from scipy.stats import f as f_dist

    r_held = problem.evaluate(f_.theta, held)
    dof = max(len(kept) - len(problem.keys) - problem.nuisance_dof(kept), 1)
    chi2_red = float(f_.residuals @ f_.residuals) / dof
    n_held = max(len(held) - problem.nuisance_dof(held), 1)
    ms_held = float(r_held @ r_held) / n_held
    stat = ms_held / chi2_red if chi2_red > 0 else np.inf
    return {"rms_held_out": float(np.sqrt(ms_held)), "rms_in_sample": float(np.sqrt(chi2_red)), "F": stat,
            "p": float(f_dist.sf(stat, n_held, dof)),
            "tension": bool(stat > f_dist.ppf(level, n_held, dof))}


def cost_drop_is_significant(cost_drop: float, chi2_red: float, df: int = 1, level: float = 0.95) -> bool:
    """A cost (0.5 x sum of squared residuals) lower by more than 0.5 x chi-square(df) quantile x reduced
    chi-square: the threshold the envelope and the pre-registrations use."""
    from scipy.stats import chi2

    return bool(cost_drop > 0.5 * float(chi2.ppf(level, df)) * chi2_red)


def summary_lines(problem: Problem, f_: Fit, lap: Dict[str, Any]) -> List[str]:
    """The numbers a run log quotes, printed from the result (prose only interprets them)."""
    lines = [f"cost {f_.start_cost:.4g} -> {f_.cost:.4g}; status {f_.status}"
             f"{' (budget exhausted)' if f_.budget_exhausted else ''}; {f_.evaluations} evaluations; {f_.seconds} s; "
             f"reduced chi-square {lap['chi2_red']:.4g}"]
    for k in problem.keys:
        s = lap["sigma"][k]
        tag = "on bound" if lap["at_bound"][k] else ("identified" if lap["identified"][k] else "UNIDENTIFIED")
        lines.append(f"  {k:34s} {f_.theta[problem.keys.index(k)]:10.4g} +- {s:.3g}  [{tag}]")
    return lines


# =====================================================================================================================
# THE MEASUREMENT MODEL (strategy M6, 2026-10-29) AND THE DIAGNOSTICS (M7)
# =====================================================================================================================
# A measured series is the truth times a calibration factor the paper may not pin down: semi-quantitative headspace
# work (response factor taken as 1) is off by an unknown factor that is SHARED by the whole series. In log space that
# is an additive offset s ~ N(0, tau^2) per series. Integrating it out exactly gives, per series, the covariance
# D + tau^2 11^T (D = diag(sigma_i^2)), and least squares on WHITENED residuals is then the exact marginal likelihood:
#     u = r / sigma;   z = u - c w^ (w^ . u),   w = 1 / sigma,   w^ = w / |w|,   c = 1 - 1 / sqrt(1 + tau^2 |w|^2)
# (Cov z = I; tests/unit/test_fitting.py checks it on simulated series.) The two limits need no number:
#   tau = 0    absolute data (isotope dilution, external calibration): z = u, the level counts;
#   tau = inf  a series whose level is unknown: c = 1, the weighted mean is projected out and only the series' SHAPE
#              (its within-study contrasts) counts, at the cost of one degree of freedom (``nuisance_dof``).
# A finite tau for calibrated data from different laboratories is a choice no source supplies; it belongs in the
# choices register as a swept cutoff, not here.
SCALE_TAU: Dict[str, float] = {
    "sida": 0.0,
    "calibrated": 0.0,
    "semi_quantitative": float("inf"),
    "relative": float("inf"),  # a series reported only relative to its own reference (percent of maximum, ...)
}


@dataclass(frozen=True)
class MeasurementModel:
    """How rows are measured: which series each belongs to, its error (dex) and its method class (a ``SCALE_TAU``
    key). Every row of one series must share a method class."""

    series_of: Callable[[Any], str]
    sigma_of: Callable[[Any], float]
    method_of: Callable[[Any], str]
    scale_tau: Mapping[str, float] = field(default_factory=lambda: dict(SCALE_TAU))
    #: extra noise per series (same units as sigma), added in quadrature to each row's own error; estimated from the
    #: data by ``estimate_series_noise`` (2026-10-29: rows with no printed error, or a vector-read error far below the
    #: series' real scatter, must not be given an invented floor)
    extra_noise: Mapping[str, float] = field(default_factory=dict)

    def sigma(self, row: Any) -> float:
        """The row's error with its series' extra noise added in quadrature."""
        s, e = float(self.sigma_of(row)), float(self.extra_noise.get(self.series_of(row), 0.0))
        return math.sqrt(s * s + e * e)

    def _series(self, rows: Sequence[Any]) -> Dict[str, List[int]]:
        out: Dict[str, List[int]] = {}
        for i, row in enumerate(rows):
            out.setdefault(self.series_of(row), []).append(i)
        for name, idx in out.items():
            methods = {self.method_of(rows[i]) for i in idx}
            if len(methods) != 1:
                raise ValueError(f"series {name!r} mixes method classes {sorted(methods)}")
            (method,) = methods
            if method not in self.scale_tau:
                raise KeyError(f"series {name!r}: unknown method class {method!r}")
        return out

    def whiten(self, log_residuals: np.ndarray, rows: Sequence[Any]) -> np.ndarray:
        """Residuals in dex -> whitened residuals, each series' calibration offset integrated out."""
        r = np.asarray(log_residuals, dtype=float)
        sigma = np.array([self.sigma(row) for row in rows])
        z = r / sigma
        for idx in self._series(rows).values():
            tau = float(self.scale_tau[self.method_of(rows[idx[0]])])
            if tau == 0.0:
                continue
            w = 1.0 / sigma[idx]
            norm2 = float(w @ w)
            c = 1.0 if math.isinf(tau) else 1.0 - 1.0 / math.sqrt(1.0 + tau * tau * norm2)
            w_hat = w / math.sqrt(norm2)
            z[idx] = z[idx] - c * w_hat * float(w_hat @ z[idx])
        return z

    def nuisance_dof(self, rows: Sequence[Any]) -> int:
        """One degree of freedom per scale-free series present in ``rows``."""
        return sum(1 for idx in self._series(rows).values()
                   if math.isinf(float(self.scale_tau[self.method_of(rows[idx[0]])])))


def estimate_series_noise(model: "MeasurementModel", raw_residuals: np.ndarray, rows: Sequence[Any],
                          n_constants: int = 0) -> Dict[str, float]:
    """Each series' extra noise from its own scatter around the current fit, by moments:
        s_k^2 = max(0, sum_i (r_i - m_k)^2 * n/(n - 1 - p_k) / n  -  mean_i sigma_i^2)
    where r are the raw residuals (dex or pH), m_k the series' mean residual if its calibration is free (tau > 0; the
    shape is all that is fitted) and 0 otherwise, and p_k the constants charged to the series (``n_constants`` spread
    over all rows). A series of one or two rows gets none (its scatter cannot be read). Iterate with the fit
    (``fit_with_series_noise``): the residuals move as the weights change."""
    r = np.asarray(raw_residuals, dtype=float)
    out: Dict[str, float] = {}
    n_total = max(len(rows), 1)
    for name, idx in model._series(rows).items():
        n = len(idx)
        tau = float(model.scale_tau[model.method_of(rows[idx[0]])])
        if n < 3:
            continue
        ri = r[idx]
        centred = ri - (float(np.mean(ri)) if tau > 0 else 0.0)
        lost = (1 if tau > 0 else 0) + n_constants * n / n_total
        dof = max(n - lost, 1.0)
        s2 = float(np.sum(centred ** 2)) / dof - float(np.mean([float(model.sigma_of(rows[i])) ** 2 for i in idx]))
        if s2 > 0:
            out[name] = math.sqrt(s2)
    return out


def fit_with_series_noise(make_problem: Callable[["MeasurementModel"], "Problem"], model: "MeasurementModel",
                          raw_residuals: Callable[[np.ndarray, Tuple[Any, ...]], np.ndarray], start: np.ndarray,
                          budget: int = 1200, rounds: int = 4, tol: float = 0.05,
                          estimator: Optional[Callable[..., Dict[str, float]]] = None) -> Tuple["Fit", "MeasurementModel"]:
    """Fit, estimate each series' extra noise from the residuals, refit with it, until the noise estimates move by
    less than ``tol`` (relative) or ``rounds`` is reached. ``make_problem(model)`` builds the Problem whose residuals
    are ``model.whiten(raw_residuals(theta, rows), rows)``. Returns the last fit and the model carrying its noise."""
    theta = np.asarray(start, dtype=float)
    noise: Dict[str, float] = dict(model.extra_noise)
    f_ = None
    for _ in range(max(1, rounds)):
        current = MeasurementModel(series_of=model.series_of, sigma_of=model.sigma_of, method_of=model.method_of,
                                   scale_tau=model.scale_tau, extra_noise=noise)
        problem = make_problem(current)
        f_ = fit(problem, theta, budget=budget)
        theta = f_.theta
        new = (estimator or estimate_series_noise)(current, raw_residuals(theta, problem.rows), problem.rows,
                                                   len(problem.keys))
        moved = max((abs(new.get(k, 0.0) - noise.get(k, 0.0)) / max(new.get(k, 0.0), noise.get(k, 0.0), 1e-12)
                     for k in set(new) | set(noise)), default=0.0)
        noise = new
        if moved < tol:
            break
    final = MeasurementModel(series_of=model.series_of, sigma_of=model.sigma_of, method_of=model.method_of,
                             scale_tau=model.scale_tau, extra_noise=noise)
    return f_, final


def count_once(rows: Sequence[Any], dataset_of: Callable[[Any], str], copy_of: Mapping[str, str]
               ) -> Tuple[Tuple[Any, ...], Tuple[Any, ...]]:
    """Duplicates counted once: drop the rows of a data set that is a confirmed copy (``copy_of[copy] = original``)
    when the original's rows are present; a copy whose original is absent stands in for it. Returns (kept,
    dropped). The copies come from the conflicts register's confirmed duplicates."""
    present = {dataset_of(row) for row in rows}
    kept, dropped = [], []
    for row in rows:
        original = copy_of.get(dataset_of(row))
        (dropped if original is not None and original in present else kept).append(row)
    return tuple(kept), tuple(dropped)


#: Reporting thresholds of the diagnostics (they flag, they never decide): a group's mean whitened residual more
#: than this many standard errors from zero; a rank correlation with a condition below this p.
DIAGNOSTIC_MEAN_SE = 2.5
DIAGNOSTIC_TREND_P = 0.01


def diagnostics(problem: Problem, theta: np.ndarray, attributes: Mapping[str, Callable[[Any], Any]],
                rows: Optional[Tuple[Any, ...]] = None) -> Dict[str, Any]:
    """The residuals at ``theta`` broken down by each named row attribute (series, laboratory, species, temperature,
    pH, time ...): per group the count, mean, RMS and share positive; per numeric attribute the rank correlation of
    residual with it. Flags the patterns the misfit protocol calls structural (docs/guides/OPERATING_PROCEDURE.md
    sec. 4): a group whose residuals share a sign beyond chance, and a residual trend with a condition."""
    from scipy.stats import spearmanr

    rows = problem.rows if rows is None else rows
    r = problem.evaluate(theta, rows)
    out: Dict[str, Any] = {"n": len(r), "rms": float(np.sqrt(np.mean(r ** 2))) if len(r) else float("nan"),
                           "by": {}, "trends": {}, "flags": []}
    for name, get in attributes.items():
        values = [get(row) for row in rows]
        groups: Dict[Any, List[int]] = {}
        for i, v in enumerate(values):
            groups.setdefault(v, []).append(i)
        table = {}
        for v, idx in sorted(groups.items(), key=lambda kv: str(kv[0])):
            g = r[idx]
            mean = float(np.mean(g))
            table[str(v)] = {"n": len(idx), "mean": mean, "rms": float(np.sqrt(np.mean(g ** 2))),
                             "share_positive": float(np.mean(g > 0))}
            if len(idx) >= 3 and abs(mean) * math.sqrt(len(idx)) > DIAGNOSTIC_MEAN_SE:
                out["flags"].append(f"{name}={v}: mean residual {mean:+.2f} over {len(idx)} rows (one sign)")
        out["by"][name] = table
        numeric = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)]
        if len(numeric) == len(values) and len(set(numeric)) >= 3:
            rho, p = spearmanr(numeric, r)
            out["trends"][name] = {"rho": float(rho), "p": float(p)}
            if p < DIAGNOSTIC_TREND_P:
                out["flags"].append(f"residual trends with {name}: rank correlation {rho:+.2f} (p = {p:.2g})")
    return out


def diagnostics_markdown(report: Dict[str, Any], title: str = "Residual diagnostics") -> str:
    """``diagnostics``' report as markdown: the flags first, then one table per attribute."""
    lines = [f"## {title}", "", f"{report['n']} rows, RMS whitened residual {report['rms']:.3g}.", ""]
    lines += (["**Flags**", ""] + [f"- {f}" for f in report["flags"]] + [""]) if report["flags"] else ["No flags.", ""]
    for name, table in report["by"].items():
        trend = report["trends"].get(name)
        lines.append(f"### By {name}" + (f" (rank correlation {trend['rho']:+.2f}, p = {trend['p']:.2g})" if trend else ""))
        lines += ["", "| group | n | mean | RMS | share > 0 |", "|---|---|---|---|---|"]
        lines += [f"| {g} | {t['n']} | {t['mean']:+.2f} | {t['rms']:.2f} | {t['share_positive']:.2f} |"
                  for g, t in table.items()]
        lines.append("")
    return "\n".join(lines)
