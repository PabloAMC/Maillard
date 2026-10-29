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

import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

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

    def __post_init__(self):
        assert len(self.keys) == len(self.kinds) == len(self.lower) == len(self.upper), "keys, kinds and bounds must align"
        assert all(k in THRESHOLDS for k in self.kinds), f"unknown coordinate kind in {self.kinds}"
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

    def jac(theta):
        return _forward_jacobian(problem, theta, lower, upper, f)

    t0 = time.time()
    r0 = f(start)
    result = least_squares(f, start, jac=jac, bounds=(lower, upper),
                           max_nfev=max(1, budget // (len(problem.keys) + 1)), **LEAST_SQUARES_SETTINGS)
    r = f(result.x)
    return Fit(theta=np.asarray(result.x, dtype=float), cost=float(0.5 * r @ r), start=start,
               start_cost=float(0.5 * r0 @ r0), status=int(result.status), evaluations=count["n"], budget=budget,
               seconds=round(time.time() - t0, 1), residuals=r)


def _forward_jacobian(problem: Problem, theta, lower, upper, f) -> np.ndarray:
    """Forward differences with B8's absolute step per coordinate kind, stepping inward at an upper bound."""
    theta = np.asarray(theta, dtype=float)
    r0 = f(theta)
    jac = np.zeros((len(r0), len(theta)))
    for j, kind in enumerate(problem.kinds):
        h = STEPS[kind] if theta[j] + STEPS[kind] <= upper[j] else -STEPS[kind]
        x = theta.copy()
        x[j] += h
        jac[:, j] = (f(x) - r0) / h
    return jac


def laplace(problem: Problem, theta: np.ndarray, rows: Optional[Tuple[Any, ...]] = None,
            upper: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """B8's Laplace at ``theta``: a central-difference Jacobian with a step per coordinate kind (one-sided at a
    bound), JtJ inverted on its eigenbasis, scaled by the reduced chi-square; a coordinate with weight in the
    dropped subspace gets an infinite sigma; B8's identification rule."""
    upper = np.asarray(problem.upper if upper is None else upper, dtype=float)
    lower = np.asarray(problem.lower, dtype=float)
    theta = np.asarray(theta, dtype=float)
    r0 = problem.evaluate(theta, rows)
    jac = np.zeros((len(r0), len(theta)))
    for j, kind in enumerate(problem.kinds):
        h = STEPS[kind]
        up, down = theta[j] + h <= upper[j], theta[j] - h >= lower[j]
        xp, xm = theta.copy(), theta.copy()
        xp[j] += h
        xm[j] -= h
        if up and down:
            jac[:, j] = (problem.evaluate(xp, rows) - problem.evaluate(xm, rows)) / (2 * h)
        elif up:
            jac[:, j] = (problem.evaluate(xp, rows) - r0) / h
        else:
            jac[:, j] = (r0 - problem.evaluate(xm, rows)) / h
    chi2_red = float(r0 @ r0) / max(len(r0) - len(theta), 1)
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
        ok = bool(np.isfinite(s) and s <= THRESHOLDS[kind])
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
    dof = max(len(kept) - len(problem.keys), 1)
    chi2_red = float(f_.residuals @ f_.residuals) / dof
    ms_held = float(r_held @ r_held) / max(len(held), 1)
    stat = ms_held / chi2_red if chi2_red > 0 else np.inf
    return {"rms_held_out": float(np.sqrt(ms_held)), "rms_in_sample": float(np.sqrt(chi2_red)), "F": stat,
            "p": float(f_dist.sf(stat, len(held), dof)),
            "tension": bool(stat > f_dist.ppf(level, len(held), dof))}


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
