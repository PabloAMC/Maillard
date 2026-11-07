"""
src/kinetic_core/predictive_densities.py -- THE MEASUREMENT MODEL AND THE PREDICTIVE DENSITY OF A UNIT (M1.1, 2026-11-05).

WHY THIS EXISTS
===============
M1 (EX-061) reached an artefact-driven verdict (results/validation/kinetic_core_m1_results_review.md): its floors
sat at measured values, it applied H0's between-laboratory tau to a C/C0 series, and its floor sensitivity used a
linearised density where the rival's prediction is strongly non-linear (exact lppd about -746 against the
linearised -175). The next three preregistered steps (M1b, M1c, the joint M1c + M2 refit) specify the fix: an
exact or adequately sampled predictive density on every form, floors from the assay, tau only where the value basis
keeps an offset, a tie band from replicate noise, an adequacy gate, Holm across guard groups, required agreement
over forms, and an exact per-species pH discrepancy. This module carries those pieces once, with tests
(tests/unit/test_predictive_densities.py), so each step's scorer calls them rather than re-implementing them.

WHAT IT IS
==========
Pure functions of arrays. Nothing here runs the engine: where an engine evaluation is needed (a likelihood at a
parameter draw, a charge-balance re-solve) the caller supplies the values or a callback.

UNITS: dex AND ln (read this before calling anything)
=====================================================
- "dex" means log10, "ln" the natural log. ``LN10 = ln 10``; x_ln = x_dex * LN10.
- Every likelihood here takes y (data), f (prediction) and cov in ONE log basis and returns the log density PER UNIT
  OF THAT BASIS. A density of n values in dex differs from the same density in ln by n * ln(LN10) nats:
  log p_ln = log p_dex - n ln(LN10) (``log_density_dex_to_ln``). Two rivals scored on the same unit in the same basis
  difference to the same d_u whichever basis is used; mixing bases inside a difference is an error.
- ``series_covariance`` takes omega and tau in dex (the units H0 reported them in) and converts them to the basis
  you name; it never guesses.
- Log densities are natural-log (nats) throughout; "lppd" is the log of the predictive density of the unit's data.

THE PREDICTIVE DENSITY
======================
lppd_u = log integral L_u(theta) p(theta) dtheta. Every estimator below takes the LOG INTEGRAND of a draw,
log L_u(theta_s) + log p(theta_s) (p normalised: Laplace Gaussian, profile density, uniform band), and the log
density of the proposal (or quadrature Gaussian) that produced theta_s:
- ``importance_lppd``: log (1/S) sum_s exp(log L + log p - log q), raw or Pareto-smoothed (PSIS, Vehtari et al.,
  JMLR 2024, through ArviZ's ``psislw``: Zhang & Stephens' estimator of the generalised Pareto tail with the weakly
  informative prior), with k-hat, the effective sample size and the acceptance rule (k-hat < 0.7, ESS >= max(100, a
  caller-given precision rule)).
- ``quadrature_lppd``: adaptive Gauss-Hermite on <= 2 integrated coordinates, centred where the caller puts it
  (the mode of L x p, with the inverse Hessian as its covariance), refusing a third dimension.
- ``plug_in_lppd`` and ``linearised_lppd``: the diagnostics (and M1's old primary).
ArviZ is wrapped here (``psis``); callers never import it.

THE REVISED DECISION RULE (sec. 7 item 5 and "The decision statistic (M1.1)", adopted 2026-11-06, PR #102)
=========================================================================================================
The last section implements the adopted rule as worded in docs/guides/STATUS_AND_PLAN.md (the prototype is
scripts/diagnostics/judge_calibration.py ``--revised``; where they differ, the plan's text is implemented):
- ``NoiseTerm`` / ``UnitNoiseModel``: the declared noise of a unit at the level its value varies (item 5). A bound is
  refused as an SD, a term estimated on another observable is refused, an instrument floor is never the whole noise,
  and the unit's ROLE follows from what was declared: "credit" (enters T), "falsification_only" (the adequacy gate
  only) or "not_decisive".
- ``conditioned_lppd`` / ``score_from_draws``: p(y_dec | y_spent) from draws, numerator and denominator each through
  the shared estimator (``importance_lppd``), experimental-fact forms mixed at their declared weights.
- ``build_null_pool`` / ``null_from_pool``: the parametric bootstrap of T under one rival as truth, truths outside the
  scoring ensemble, re-weighted to that rival's spent-conditioned posterior.
- ``judge_pair``: the verdict (``RevisedVerdict``) from per-form densities, adequacy p's and the two nulls;
  ``judge_from_draws`` runs all of it from draws.
"""
from __future__ import annotations

import itertools
import math
import warnings
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

# ----------------------------------------------------------------------------------------------------------------
# Constants and the dex <-> ln conversion
# ----------------------------------------------------------------------------------------------------------------

LN10 = math.log(10.0)
#: H0 (EX-059): the within-laboratory model discrepancy, a fractional error per predicted value, in dex. The ship
#: value is the WORSE of H0's two sets, 0.419 (warm; H0 prereg sec. 3.6 / D2, amendment A1 item 3): the plan, the H0
#: application and the M1b / M1c / M3 preregs use it. The reset set's 0.390 is a named sensitivity.
OMEGA_DEX = 0.419
OMEGA_DEX_RESET = 0.390
OMEGA_LN = OMEGA_DEX * LN10
#: H0's tau_ABS(all), dex: the between-laboratory offset of an ABSOLUTE level. The caller decides per series.
TAU_ABS_DEX = 0.964
#: the one-sided 95 % normal quantile: b_u = Z_90 x sd_rep is the 90 % half-width (the plan's interval convention)
Z_90 = 1.6448536269514722
#: the acceptance rule for an importance-sampled density (preregs M1b sec. 5.1, M1c + M2 sec. 7.1)
K_HAT_MAX = 0.7
ESS_MIN = 100.0
#: a floor this close (relative) to a data value is refused: M1's floors were the 30- and 50-min data points
#: (0.80 mM against a measured 0.8016, 1.22 against 1.22, 45 ug/L against 45)
FLOOR_DATA_RTOL = 0.01
#: the value bases a series may declare; only an ABSOLUTE level keeps a between-laboratory calibration offset
VALUE_BASES = ("absolute", "ratio", "c_over_c0", "first_point_normalised", "contrast", "share")
BASES = ("dex", "ln")
#: pH scales (M1b sec. 3.3, V5): pH after cooling, the label's 25 C value, or in situ at the process temperature
PH_SCALES = ("cooled", "label_25C", "in_situ")


class MeasurementModelRefusal(ValueError):
    """A measurement-model choice the design rules forbid (tau on a relative series, a floor at a datum, a mixed pH
    scale, a third quadrature dimension, ...). Raised, never silently corrected."""


def dex_to_ln(x):
    """A log10 quantity (a residual, an sd, omega) in natural-log units: x * ln 10. Variances scale by (ln 10)^2."""
    return np.asarray(x, dtype=float) * LN10 if np.ndim(x) else float(x) * LN10


def ln_to_dex(x):
    """A natural-log quantity in log10 units: x / ln 10."""
    return np.asarray(x, dtype=float) / LN10 if np.ndim(x) else float(x) / LN10


def log_density_dex_to_ln(log_density: float, n_values: int) -> float:
    """A log density of n values written in log10 (per dex^n) as a density per ln^n: subtract n ln(ln 10)."""
    return float(log_density) - int(n_values) * math.log(LN10)


def log_density_ln_to_dex(log_density: float, n_values: int) -> float:
    """The inverse of ``log_density_dex_to_ln``."""
    return float(log_density) + int(n_values) * math.log(LN10)


def _check_basis(basis: str) -> float:
    if basis not in BASES:
        raise MeasurementModelRefusal(f"basis must be one of {BASES}, got {basis!r}")
    return 1.0 if basis == "dex" else LN10


def _logsumexp(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return float("-inf")
    m = float(np.max(finite))
    return m + math.log(float(np.sum(np.exp(finite - m))))


# ----------------------------------------------------------------------------------------------------------------
# Likelihood pieces
# ----------------------------------------------------------------------------------------------------------------

def series_covariance(measurement_sd: Sequence[float], value_basis: str, basis: str = "dex",
                      omega_dex: Optional[float] = OMEGA_DEX, tau_dex: Optional[float] = None,
                      design: Optional[np.ndarray] = None, legacy_tau_override: bool = False) -> np.ndarray:
    """The covariance of a series' values, in ``basis`` units squared.

    ``measurement_sd`` is the sd of each measured POINT, already in ``basis`` units. ``design`` (values x points)
    maps points to the scored values (identity by default; a ratio 360/30 is the row [-1, +1, 0], and two ratios
    sharing the 360 min point are correlated through it). Each predicted point carries the model discrepancy omega
    independently (H0: omega per printed value), so

        Sigma = A diag(measurement_sd^2 + omega^2) A^T + tau^2 1 1^T

    with omega and tau given in dex and converted to ``basis``. tau is the shared offset of an absolute level
    (the level integrated out under N(0, tau^2)); it is refused unless ``value_basis == "absolute"``: in a ratio, a
    C/C0 or a normalised series a multiplicative calibration offset cancels, and putting tau there hides the misfit
    (M1's error). ``legacy_tau_override`` exists ONLY to reproduce a frozen record that made that error (M1's
    decisive comparison); a new score never sets it. ``omega_dex=None`` (or 0) leaves omega out (M1's old primary)."""
    scale = _check_basis(basis)
    if value_basis not in VALUE_BASES:
        raise MeasurementModelRefusal(f"value_basis must be one of {VALUE_BASES}, got {value_basis!r}")
    sd = np.asarray(measurement_sd, dtype=float)
    if sd.ndim != 1 or not np.all(np.isfinite(sd)) or np.any(sd < 0):
        raise MeasurementModelRefusal("measurement_sd must be a finite, non-negative vector (one per point)")
    a = np.eye(len(sd)) if design is None else np.asarray(design, dtype=float)
    if a.ndim != 2 or a.shape[1] != len(sd):
        raise MeasurementModelRefusal(f"design must be (values x {len(sd)} points), got {a.shape}")
    om = 0.0 if not omega_dex else float(omega_dex) * scale
    cov = a @ np.diag(sd ** 2 + om ** 2) @ a.T
    if tau_dex:
        if value_basis != "absolute" and not legacy_tau_override:
            raise MeasurementModelRefusal(
                f"tau refused on a {value_basis!r} series: an offset cancels in a relative value basis "
                "(lessons 2026-11-05); only an 'absolute' level carries a between-laboratory tau")
        t = float(tau_dex) * scale
        cov = cov + t ** 2 * np.ones((a.shape[0], a.shape[0]))
    return 0.5 * (cov + cov.T)


def gaussian_log_density(y: Sequence[float], f: Sequence[float], cov: np.ndarray) -> float:
    """log N(y; f, cov), nats, per unit of the basis y, f and cov share. A non-finite residual (a failed run, a
    zero prediction on a log scale) gives -inf. Equal to M1's ``log_normal``."""
    r = np.asarray(y, dtype=float) - np.asarray(f, dtype=float)
    if not np.all(np.isfinite(r)):
        return float("-inf")
    cov = np.atleast_2d(np.asarray(cov, dtype=float))
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        raise MeasurementModelRefusal("the covariance is not positive definite")
    return float(-0.5 * (r @ np.linalg.solve(cov, r) + logdet + len(r) * math.log(2.0 * math.pi)))


def assay_floor(lod: Optional[float] = None, loq: Optional[float] = None, use: str = "lod") -> float:
    """The detection floor from the ASSAY (a limit of detection or quantification, a calibration range's lower end),
    never from a measured point. Returns the chosen limit (same units as given)."""
    if use not in ("lod", "loq"):
        raise MeasurementModelRefusal(f"use must be 'lod' or 'loq', got {use!r}")
    v = lod if use == "lod" else loq
    if v is None or not (math.isfinite(float(v)) and float(v) > 0):
        raise MeasurementModelRefusal(f"no positive {use.upper()} given: a floor comes from the assay")
    return float(v)


def check_floor_not_data(floor: float, data_values: Iterable[float], rtol: float = FLOOR_DATA_RTOL) -> None:
    """Refuse a floor that equals a data value (relative ``rtol``, linear scale): a floor at a measured point
    flattens a rival onto that datum with zero width (M1's primary floors)."""
    fl = float(floor)
    for v in data_values:
        v = float(v)
        if math.isfinite(v) and abs(v - fl) <= rtol * max(abs(v), abs(fl)):
            raise MeasurementModelRefusal(
                f"floor {fl:g} equals the data value {v:g} (within {rtol:g} relative): a floor comes from the assay "
                "(LOD / LOQ), never from a measured point")


def censored_gaussian_log_density(y: Sequence[float], f: Sequence[float], cov: np.ndarray,
                                  below_floor: Sequence[bool], log_floor: Sequence[float]) -> float:
    """The log likelihood of a series in which some values are reported below the assay's floor.

    Observed values (``below_floor`` False) enter as a Gaussian density at y; censored values enter as the
    probability that the value lies below its floor, conditional on the observed ones:

        log p = log N(y_o; f_o, S_oo) + log P(Y_c < log_floor_c | y_o)

    where Y_c | y_o is Gaussian with the usual conditional mean and covariance (one censored value: the normal cdf;
    several: scipy's multivariate normal cdf). y, f, cov and ``log_floor`` are in one log basis (the floor is the log
    of the assay limit in that basis); y at a censored index is ignored. Refuses a censored value whose floor is not
    finite."""
    from scipy.stats import multivariate_normal, norm

    y, f = np.asarray(y, dtype=float), np.asarray(f, dtype=float)
    cov = np.atleast_2d(np.asarray(cov, dtype=float))
    c = np.asarray(below_floor, dtype=bool)
    lf = np.broadcast_to(np.asarray(log_floor, dtype=float), y.shape)
    if np.any(c & ~np.isfinite(lf)):
        raise MeasurementModelRefusal("a censored value needs a finite floor")
    o = ~c
    if not np.all(np.isfinite(f)):
        return float("-inf")
    out = 0.0
    mu_c, s_cc = f[c], cov[np.ix_(c, c)]
    if o.any():
        out += gaussian_log_density(y[o], f[o], cov[np.ix_(o, o)])
        if c.any():
            s_co, s_oo = cov[np.ix_(c, o)], cov[np.ix_(o, o)]
            k = np.linalg.solve(s_oo, s_co.T).T
            mu_c = mu_c + k @ (y[o] - f[o])
            s_cc = s_cc - k @ s_co.T
    if c.any():
        if c.sum() == 1:
            out += float(norm.logcdf(lf[c][0], loc=mu_c[0], scale=math.sqrt(float(s_cc[0, 0]))))
        else:
            p = float(multivariate_normal(mean=mu_c, cov=0.5 * (s_cc + s_cc.T), allow_singular=False).cdf(lf[c]))
            out += math.log(p) if p > 0 else float("-inf")
    return float(out)


def legacy_floored_prediction(f_log: Sequence[float], log_floor: float) -> Tuple[np.ndarray, np.ndarray]:
    """M1's withdrawn form (re-review B1a): a prediction below the floor scored AT the floor (a zero-width line at
    the floor). Returns the floored predictions and the mask of floored points. Here only to reproduce M1's record;
    a new score uses ``censored_gaussian_log_density`` on DATA reported below an assay floor."""
    f = np.asarray(f_log, dtype=float)
    mask = np.isfinite(f) & (f < float(log_floor))
    return np.where(mask, float(log_floor), f), mask


# ----------------------------------------------------------------------------------------------------------------
# The predictive density of a unit
# ----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class DensityEstimate:
    """One unit's lppd under one rival and one form, with the estimator's diagnostics."""

    lppd: float
    method: str                     # "plug_in" | "linearised" | "importance" | "psis" | "quadrature"
    accepted: bool
    reason: str = ""
    k_hat: Optional[float] = None
    ess: Optional[float] = None
    ess_required: Optional[float] = None
    mc_se: Optional[float] = None
    n_draws: Optional[int] = None
    n_nonfinite: int = 0
    extra: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        out = {k: getattr(self, k) for k in ("lppd", "method", "accepted", "reason", "k_hat", "ess", "ess_required",
                                              "mc_se", "n_draws", "n_nonfinite")}
        out.update(self.extra)
        return out


def plug_in_lppd(y: Sequence[float], f: Sequence[float], cov: np.ndarray) -> DensityEstimate:
    """The density at the point prediction (no parameter uncertainty): a diagnostic and a required-agreement form."""
    return DensityEstimate(gaussian_log_density(y, f, cov), "plug_in", True)


def linearised_lppd(y: Sequence[float], f: Sequence[float], cov: np.ndarray, jac: np.ndarray,
                    sigma_theta: np.ndarray) -> DensityEstimate:
    """The linearised (delta-method) predictive: N(y; f, cov + J Sigma_theta J^T). J is d f / d theta (values x
    coordinates, in the basis of f per coordinate unit), Sigma_theta the coordinates' covariance. Exact only where f is
    linear in theta over the posterior's width; M1's floor sensitivity showed how far that fails (-175 vs -746)."""
    j = np.atleast_2d(np.asarray(jac, dtype=float))
    extra = j @ np.asarray(sigma_theta, dtype=float) @ j.T
    return DensityEstimate(gaussian_log_density(y, f, np.asarray(cov, dtype=float) + extra), "linearised", True)


def psis(log_ratios: Sequence[float], r_eff: float = 1.0) -> Tuple[np.ndarray, float]:
    """Pareto-smoothed importance sampling (Vehtari et al., JMLR 2024) through ArviZ 0.23's ``psislw``: the largest
    M = ceil(min(S / 5, 3 sqrt(S / r_eff))) log ratios are replaced by the expected order statistics of a generalised
    Pareto fitted to them (Zhang & Stephens with the weakly informative prior on k), truncated at the largest raw
    ratio. Returns the smoothed, UNNORMALISED log ratios (on the input's scale, so logsumexp of them estimates the
    integral) and k-hat. k-hat is inf when the tail has <= 4 draws (too few for the fit)."""
    lr = np.asarray(log_ratios, dtype=float)
    if lr.ndim != 1:
        raise MeasurementModelRefusal("psis takes one vector of log ratios")
    if np.any(np.isnan(lr)):
        raise MeasurementModelRefusal("a NaN log ratio: a failed engine run must be decided by the caller")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # ArviZ's import-time FutureWarning about its 1.x refactor
        import arviz

        lw, k = arviz.psislw(lr.copy(), reff=float(r_eff), normalize=False)
    return np.asarray(lw, dtype=float), float(np.asarray(k))


def effective_sample_size(log_weights: Sequence[float]) -> float:
    """(sum w)^2 / sum w^2 for importance weights given as logs (any common scale)."""
    lw = np.asarray(log_weights, dtype=float)
    lw = lw[np.isfinite(lw)]
    if lw.size == 0:
        return 0.0
    w = np.exp(lw - np.max(lw))
    return float(np.sum(w) ** 2 / np.sum(w ** 2))


def tie_band_precision_ess(sd_rep: float, fraction: float = 0.15) -> float:
    """The preregs' precision rule tying the MC error to the tie band: SE_MC(lppd) <= fraction x sd_rep / sqrt 2 with
    SE_MC ~ 1 / sqrt(ESS), i.e. ESS >= 2 / (fraction x sd_rep)^2. Feed it to ``importance_lppd(precision_ess=...)``."""
    if not (sd_rep > 0 and math.isfinite(sd_rep)):
        return float("inf")
    return 2.0 / (fraction * float(sd_rep)) ** 2


def importance_lppd(log_integrand: Sequence[float], log_proposal: Sequence[float], smooth: bool = True,
                    precision_ess: float = 0.0, k_max: float = K_HAT_MAX, ess_min: float = ESS_MIN,
                    r_eff: float = 1.0) -> DensityEstimate:
    """lppd = log (1/S) sum_s exp(r_s), r_s = log_integrand_s - log_proposal_s, theta_s drawn from the proposal q.

    ``log_integrand`` is log L_u(theta_s) + log p(theta_s) with p NORMALISED (then the estimate is the predictive
    density itself, not a ratio of normalisers). ``smooth`` applies PSIS to the r_s (``psis``) and reports k-hat;
    ``smooth=False`` is the raw estimator (k-hat still reported). Draws with r_s = -inf (zero likelihood, a failed
    run the caller scored as impossible) count in S and are left out of the tail fit. The MC SE is the delta-method
    SE of the log of the mean weight. Accepted iff k-hat < ``k_max`` and ESS >= max(``ess_min``, ``precision_ess``)
    (the preregs' rule); otherwise the estimate is returned with ``accepted`` False and the reason."""
    li, lq = np.asarray(log_integrand, dtype=float), np.asarray(log_proposal, dtype=float)
    if li.shape != lq.shape or li.ndim != 1:
        raise MeasurementModelRefusal("log_integrand and log_proposal must be equal-length vectors")
    if np.any(np.isnan(li)) or np.any(np.isnan(lq)) or np.any(~np.isfinite(lq)):
        raise MeasurementModelRefusal("NaN log integrand or a non-finite proposal density: decide failed runs first")
    s = li.size
    r = li - lq
    fin = np.isfinite(r)
    n_bad = int(s - fin.sum())
    need = max(float(ess_min), float(precision_ess))
    if fin.sum() == 0:
        return DensityEstimate(float("-inf"), "psis" if smooth else "importance", False, "every draw has zero weight",
                               n_draws=s, n_nonfinite=n_bad, ess=0.0, ess_required=need)
    lw_fin, k = psis(r[fin], r_eff)
    lw = np.full(s, -np.inf)
    lw[fin] = lw_fin if smooth else r[fin]
    lppd = _logsumexp(lw) - math.log(s)
    ess = effective_sample_size(lw)
    w = np.exp(lw - np.max(lw[fin]))
    mean = float(np.mean(w))
    se = float(np.std(w, ddof=1) / (math.sqrt(s) * mean)) if s > 1 and mean > 0 else float("nan")
    reasons = []
    if not (k < k_max):
        reasons.append(f"k-hat {k:.3g} >= {k_max}")
    if not (ess >= need):
        reasons.append(f"ESS {ess:.1f} < {need:.1f}")
    return DensityEstimate(float(lppd), "psis" if smooth else "importance", not reasons, "; ".join(reasons),
                           k_hat=float(k), ess=float(ess), ess_required=need, mc_se=se, n_draws=s, n_nonfinite=n_bad)


def multivariate_t_draws(mean: Sequence[float], scale_cov: np.ndarray, n: int, rng: np.random.Generator,
                         df: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """The preregs' unit-tuned proposal: a multivariate Student-t (5 df) at ``mean`` with scale matrix ``scale_cov``
    (the inverse Hessian of -log(L x p) at the mode). Returns (draws n x d, log q at each draw)."""
    from scipy.stats import multivariate_t

    dist = multivariate_t(loc=np.atleast_1d(np.asarray(mean, dtype=float)),
                          shape=np.atleast_2d(np.asarray(scale_cov, dtype=float)), df=df)
    x = np.atleast_2d(dist.rvs(size=n, random_state=rng))
    x = x.reshape(n, -1)
    return x, np.atleast_1d(dist.logpdf(x)).astype(float)


def gaussian_parameter_log_density(theta: np.ndarray, mean: Sequence[float], cov: np.ndarray) -> np.ndarray:
    """log N(theta_s; mean, cov) for each row of ``theta`` (the Laplace density of identified coordinates)."""
    from scipy.stats import multivariate_normal

    x = np.atleast_2d(np.asarray(theta, dtype=float))
    m = np.atleast_1d(np.asarray(mean, dtype=float))
    return np.atleast_1d(multivariate_normal(mean=m, cov=np.atleast_2d(cov)).logpdf(x.reshape(-1, m.size)))


def quadrature_lppd(log_integrand: Callable[[np.ndarray], np.ndarray], center: Sequence[float], cov: np.ndarray,
                    n_nodes: int = 40) -> DensityEstimate:
    """Adaptive Gauss-Hermite: lppd = log integral exp(g(theta)) dtheta with g = log L + log p, for <= 2 integrated
    coordinates. The nodes are those of N(center, cov) (put ``center`` at the mode of g and ``cov`` at the inverse of
    its negative Hessian: then the rule is exact for a Gaussian integrand and needs few nodes); each node carries
    exp(g - log N(.; center, cov)). ``log_integrand`` takes an (m x d) array of points and returns m values. Refuses
    d > 2 (the preregs' fallback is for the one or two dominant coordinates only)."""
    c = np.atleast_1d(np.asarray(center, dtype=float))
    d = c.size
    if d > 2:
        raise MeasurementModelRefusal(f"quadrature on {d} coordinates refused: the fallback covers <= 2")
    cv = np.atleast_2d(np.asarray(cov, dtype=float))
    chol = np.linalg.cholesky(cv)
    z, wz = np.polynomial.hermite_e.hermegauss(int(n_nodes))     # probabilists' Hermite: weight exp(-z^2 / 2)
    wz = wz / math.sqrt(2.0 * math.pi)
    if d == 1:
        zz, ww = z[:, None], wz
    else:
        g1, g2 = np.meshgrid(z, z, indexing="ij")
        zz = np.column_stack([g1.ravel(), g2.ravel()])
        ww = np.outer(wz, wz).ravel()
    pts = c + zz @ chol.T
    g = np.asarray(log_integrand(pts), dtype=float).reshape(-1)
    if g.size != pts.shape[0]:
        raise MeasurementModelRefusal("log_integrand must return one value per point")
    if np.any(np.isnan(g)):
        raise MeasurementModelRefusal("NaN log integrand at a quadrature node")
    log_n = (-0.5 * np.sum(zz ** 2, axis=1) - 0.5 * d * math.log(2.0 * math.pi)
             - float(np.sum(np.log(np.diag(chol)))))
    lppd = _logsumexp(np.log(ww) + g - log_n)
    return DensityEstimate(float(lppd), "quadrature", bool(np.isfinite(lppd)),
                           "" if np.isfinite(lppd) else "zero integrand at every node", n_draws=int(pts.shape[0]),
                           extra={"dimensions": d, "nodes_per_dimension": int(n_nodes)})


# ----------------------------------------------------------------------------------------------------------------
# The decision helpers
# ----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class TieBand:
    sd_rep: float
    band: float
    n_replicates: int


def gaussian_replicates(y: Sequence[float], measurement_cov: np.ndarray, n: int, rng: np.random.Generator
                        ) -> np.ndarray:
    """R replicate data vectors y + e, e ~ N(0, measurement_cov): MEASUREMENT noise only (never omega)."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    return rng.multivariate_normal(y, np.atleast_2d(measurement_cov), size=int(n), check_valid="raise")


def tie_band(d_of_y: Callable[[np.ndarray], float], replicates: np.ndarray, z: float = Z_90) -> TieBand:
    """sd_rep(d_u): the spread of d_u = lppd_u(a) - lppd_u(b) over replicate data (``gaussian_replicates``, or the
    caller's own redraw rule, e.g. M1b's censoring-edge redraw), each evaluated at the FIXED predictive distributions
    by ``d_of_y``; b_u = z x sd_rep (z_0.95: the 90 % half-width). Non-finite replicate values are refused."""
    reps = np.atleast_2d(np.asarray(replicates, dtype=float))
    d = np.array([float(d_of_y(row)) for row in reps])
    if not np.all(np.isfinite(d)):
        raise MeasurementModelRefusal("a replicate gave a non-finite d_u")
    sd = float(np.std(d, ddof=1)) if d.size > 1 else float("nan")
    return TieBand(sd, z * sd, int(d.size))


def posterior_predictive_p(y: Sequence[float], mean: Sequence[float], cov: np.ndarray) -> float:
    """The adequacy p of a rival on a unit: the chi^2 tail of the Mahalanobis discrepancy of y from the predictive
    mean under the predictive covariance (the unit's Sigma with omega, plus the predictive spread). For one value this
    is the two-sided normal p."""
    from scipy.stats import chi2

    r = np.atleast_1d(np.asarray(y, dtype=float) - np.asarray(mean, dtype=float))
    if not np.all(np.isfinite(r)):
        return 0.0
    m2 = float(r @ np.linalg.solve(np.atleast_2d(cov), r))
    return float(chi2.sf(m2, r.size))


def unit_reading(d: float, band: float, p_a: float, p_b: float, mean_a: Sequence[float], mean_b: Sequence[float],
                 discrepancy_scale: float, p_min: float = 0.05) -> str:
    """A unit's reading: "a", "b", "tie" or "both inadequate".

    |d| <= band (or d not finite): "tie". Otherwise the side d points to (d = lppd(a) - lppd(b) > band: a) is
    credited only if that rival is ADEQUATE there (posterior-predictive p > ``p_min``) or the rivals' predictive means
    differ by more than ``discrepancy_scale`` (omega, or omega sqrt(c) / delta-pH, on the unit's scale; for a vector
    unit the largest component difference); otherwise "both inadequate"."""
    if not (math.isfinite(d) and math.isfinite(band)) or abs(d) <= band:
        return "tie"
    side, p = ("a", p_a) if d > 0 else ("b", p_b)
    gap = float(np.max(np.abs(np.atleast_1d(np.asarray(mean_a, float) - np.asarray(mean_b, float)))))
    if p > p_min or gap > discrepancy_scale:
        return side
    return "both inadequate"


def holm(p_values: Mapping[str, float], alpha: float = 0.05) -> Dict[str, Dict[str, Any]]:
    """Holm's step-down across groups at family-wise ``alpha``: sort p ascending, reject p_(i) while
    p_(i) <= alpha / (m - i + 1) (i from 1), stop at the first failure. Returns per group the Holm-adjusted p
    (monotone, capped at 1) and whether it is rejected."""
    items = sorted(p_values.items(), key=lambda kv: (kv[1], kv[0]))
    m = len(items)
    out: Dict[str, Dict[str, Any]] = {}
    running, still = 0.0, True
    for i, (g, p) in enumerate(items):
        if not (0.0 <= p <= 1.0):
            raise MeasurementModelRefusal(f"p for {g!r} is not a probability: {p!r}")
        running = max(running, min(1.0, (m - i) * p))
        still = still and p <= alpha / (m - i)
        out[g] = {"p": float(p), "p_holm": float(running), "rejected": bool(still)}
    return out


def required_agreement(d_by_form: Mapping[str, float], band_by_form: Mapping[str, float]) -> Dict[str, Any]:
    """A contrast's side over its required-agreement forms: "a" only if EVERY form has d > its own tie band, "b"
    only if every form has d < -band; otherwise "tie", with the forms that broke it. Every form needs a band (pass 0
    explicitly for a sign-only form)."""
    missing = sorted(set(d_by_form) - set(band_by_form))
    if missing:
        raise MeasurementModelRefusal(f"no tie band for forms {missing}")
    if not d_by_form:
        raise MeasurementModelRefusal("no forms")
    sides = {}
    for form, d in d_by_form.items():
        b = float(band_by_form[form])
        sides[form] = "a" if (math.isfinite(d) and d > b) else "b" if (math.isfinite(d) and d < -b) else "tie"
    vals = set(sides.values())
    side = vals.pop() if len(vals) == 1 and "tie" not in vals else "tie"
    majority = "a" if sum(v == "a" for v in sides.values()) >= sum(v == "b" for v in sides.values()) else "b"
    breaking = [] if side != "tie" else sorted(f for f, v in sides.items() if v != majority)
    return {"side": side, "per_form": sides, "breaking_forms": breaking}


# ----------------------------------------------------------------------------------------------------------------
# pH: the exact per-species discrepancy and the scale guard
# ----------------------------------------------------------------------------------------------------------------

def require_ph_scale(*scales: str) -> str:
    """Every pH in a comparison declares its scale ("cooled", "label_25C", "in_situ"); a mixed-scale comparison is
    refused (M1b V5). Returns the common scale."""
    if not scales:
        raise MeasurementModelRefusal("no pH scale declared")
    bad = [s for s in scales if s not in PH_SCALES]
    if bad:
        raise MeasurementModelRefusal(f"undeclared or unknown pH scale {bad}; declare one of {PH_SCALES}")
    if len(set(scales)) > 1:
        raise MeasurementModelRefusal(f"mixed-scale pH comparison refused: {sorted(set(scales))}")
    return scales[0]


@dataclass(frozen=True)
class PhReading:
    value: float
    scale: str
    temperature_c: Optional[float] = None


def ph_difference(predicted: PhReading, measured: PhReading) -> float:
    """predicted - measured, after ``require_ph_scale`` on both."""
    require_ph_scale(predicted.scale, measured.scale)
    return float(predicted.value) - float(measured.value)


def ph_species_discrepancy(resolve: Callable[[Mapping[str, float]], float], produced: Mapping[str, float],
                           acids: Iterable[str], omega_dex: float = OMEGA_DEX, scale: str = "cooled",
                           per_species: bool = False) -> Dict[str, Any]:
    """The model-discrepancy delta-pH of a predicted pH, by exact balance re-solves (M1b sec. 3.3, v3).

    ``produced`` holds the GROSS amount each produced titratable species was formed in (never a net change, which can
    be negative: refused). ``acids`` names the species that lower the pH as their amount rises; the rest raise it.
    Each species is scaled independently by 10^(+-omega); pH being monotone in each amount, the two corners (every
    acid up and every base down; the reverse) bound every combination, so

        delta-pH = | pH(acids x 10^omega, bases x 10^-omega) - pH(acids x 10^-omega, bases x 10^omega) | / 2

    ``resolve`` is the caller's charge-balance re-solve: amounts -> pH on ``scale`` (declared, checked). Returns the
    central pH, both corners, delta-pH, whether the corners bracket the central pH (the monotonicity the bound needs;
    flagged, not assumed), and with ``per_species`` each species' own half-spread (2 more solves per species)."""
    require_ph_scale(scale)
    amounts = {k: float(v) for k, v in produced.items()}
    if not amounts:
        raise MeasurementModelRefusal("no produced species")
    bad = [k for k, v in amounts.items() if not (math.isfinite(v) and v >= 0)]
    if bad:
        raise MeasurementModelRefusal(f"gross produced amounts must be finite and >= 0 (not net): {bad}")
    acid = set(acids)
    unknown = sorted(acid - set(amounts))
    if unknown:
        raise MeasurementModelRefusal(f"acids not among the produced species: {unknown}")
    up, down = 10.0 ** float(omega_dex), 10.0 ** (-float(omega_dex))
    central = float(resolve(dict(amounts)))
    acid_up = float(resolve({k: v * (up if k in acid else down) for k, v in amounts.items()}))
    base_up = float(resolve({k: v * (down if k in acid else up) for k, v in amounts.items()}))
    out: Dict[str, Any] = {"scale": scale, "omega_dex": float(omega_dex), "ph_central": central,
                           "ph_acids_up": acid_up, "ph_bases_up": base_up,
                           "delta_ph": 0.5 * abs(base_up - acid_up),
                           "corners_bracket_central": bool(acid_up <= central <= base_up)}
    if per_species:
        each = {}
        for k in amounts:
            hi = float(resolve({**amounts, k: amounts[k] * up}))
            lo = float(resolve({**amounts, k: amounts[k] * down}))
            each[k] = 0.5 * abs(hi - lo)
        out["per_species_half_spread"] = each
    return out


# ----------------------------------------------------------------------------------------------------------------
# The revised decision rule (sec. 7 item 5 and M1.1, adopted 2026-11-06, PR #102)
# ----------------------------------------------------------------------------------------------------------------
#
# DECLARED READINGS (where the plan's wording needs a reading to be computed)
# - T is the JOINT log predictive density of the credit units given the spent rows, rival a minus rival b:
#   T = log p_a(y_credit | y_spent) - log p_b(y_credit | y_spent) (the prototype's reading of "the pooled log predictive
#   density difference over the decisive units ... p(y_dec | y_spent)"). T_u is unit u alone given the spent rows, T_-u
#   the joint without u.
# - Mirror rule: b is credited when T < the 5th percentile of T under a as truth (i.e. -T exceeds the 95th percentile
#   of -T under a).
# - "Resolved" applies to every density the verdict reads: each rival's numerator for T, T_u and T_-u and the shared
#   denominator p(y_spent), in every scored form and experimental-fact component. An unresolved OBSERVED density makes
#   the pair "not computable". An unresolved bootstrap REPLICATE (any of its densities, in any form) is re-estimated
#   once when the caller supplies a retry draw set (``build_null_pool(retry_draws=...)``); one still unresolved is
#   never dropped or read as estimated but bounded (PR #105): in testing a's credit its T, T_u and T_-u are +inf in
#   BOTH nulls, in testing b's -inf, and a credit stands only if it survives the bound. +inf raises null b's upper
#   thresholds (harder for a) and null a's lower ones (easier for b to qualify too, which blocks a): the bound always
#   works against the credit tested. The crossing is strict in the null's weighted tail (a: W(T_null >= T) < alpha),
#   so a null whose weighted unresolved fraction reaches alpha cannot credit.
# - A credit needs the power check (``misspecified_false_credit``) supplied and passing: without it no credit.
# - ">= 1000 replicates" is the null's effective count: the Kish ESS of its weights (the replicate count when the
#   replicates are unweighted). A re-weighted null also needs the weights' PSIS k-hat < 0.7.
# - Adequacy family: every unit whose role is "credit" or "falsification_only" (an absolute level that cannot credit
#   still tests adequacy, as in the prototype); posterior-predictive p of each under p(y_u | y_spent) (the plan's
#   text: conditioned on the spent rows only), Holm across them, adequate iff every Holm-adjusted p >= 0.05.
#   CAVEAT: the judge calibration's 0.70 - 0.90 power (PR #99 / #102) used the prototype's adequacy, conditioned on
#   the spent AND decisive rows; this spent-only adequacy is uncalibrated except in synthetic worlds.
# - Integration over an estimate's sampling uncertainty: a scaled inverse chi-square on the variance (estimate s, df
#   nu), by a K-point quantile quadrature (K = 4, the prototype's); the bootstrap draws the true variance from it and
#   re-draws the estimate a replicate would print (s_rep^2 = sigma^2 chi^2_nu / nu).

NOMINAL_ALPHA = 0.05
NULL_REPLICATES_MIN = 1000
NULL_INPUTS = ("decisive_data", "conditioning_rows", "forced_inputs", "replicate_sds", "discrepancy_component",
               "form_index")
NOISE_TERM_KINDS = (
    "run_sd",             # the printed SD of independent runs (n - 1 df)
    "reproducibility",    # the same lab class's measured run-to-run reproducibility of the same observable
    "within_run_sd",      # a printed within-run replicate SD, or a same-observable within-run repeat
    "within_run_prior",   # a within-run per-point component integrated over a prior declared with items 1 - 5
    "instrument_floor",   # an instrument's repeatability or a reading error: a floor, never the whole noise
    "discrepancy",        # a model-discrepancy variance component (same observable and lab class, unfitted rows)
    "between_lab",        # the measured between-laboratory spread of the same observable and value basis
    "lane_within_lab",    # the lane's within-laboratory term (S1 - S4 preamble): falsification only
)
NOISE_STRUCTURES = ("per_point", "run_offset")
VALUE_LEVELS = ("within_run", "absolute", "cross_lab")
UNIT_ROLES = ("credit", "falsification_only", "not_decisive")
FORM_KINDS = ("experimental_fact", "modelling_choice", "omega_hedge")
QUAD_K = 4
MAX_QUAD_NODES = 4096
VERDICT_OUTCOMES = ("credited", "neither_credited_both_qualify", "undecided", "carried_sole_survivor",
                    "both_falsified", "not_computable", "not_decisive")
_RUN_LEVEL_KINDS = ("run_sd", "reproducibility")
_PER_POINT_NOISE_KINDS = ("run_sd", "reproducibility", "within_run_sd", "within_run_prior", "discrepancy")
_DF_REQUIRED = ("run_sd", "within_run_sd", "within_run_prior", "discrepancy")


def invchi2_quadrature(nu: float, k: int = QUAD_K) -> np.ndarray:
    """K variance multipliers at the (i - 1/2) / K quantiles of the scaled inverse chi-square posterior of sigma^2
    given an estimate s with nu df: sigma^2 = s^2 x multiplier (the prototype's ``invchi2_quad``)."""
    from scipy.stats import chi2

    u = (np.arange(int(k)) + 0.5) / int(k)
    return float(nu) / chi2.ppf(1.0 - u, float(nu))


def weighted_quantile(x: np.ndarray, w: Optional[np.ndarray], q: float) -> Tuple[float, float, float]:
    """(quantile, its Monte Carlo SE, ESS) of a weighted sample (the prototype's estimator, kept identical): the
    smallest x whose cumulative weight reaches q; SE = sqrt(q (1 - q) / ESS) / f(quantile), f a weighted Gaussian kernel
    density (Silverman bandwidth)."""
    x = np.asarray(x, dtype=float)
    w = np.ones_like(x) if w is None else np.asarray(w, dtype=float)
    ok = np.isfinite(x) & (w > 0)
    x, w = x[ok], w[ok]
    if x.size == 0:
        return float("nan"), float("nan"), 0.0
    w = w / w.sum()
    o = np.argsort(x)
    xs, cw = x[o], np.cumsum(w[o])
    qv = float(xs[min(np.searchsorted(cw, q), xs.size - 1)])
    ess = float(1.0 / np.sum(w ** 2))
    mu = float(w @ x)
    sd = float(np.sqrt(max(w @ (x - mu) ** 2, 1e-300)))
    h = 1.06 * sd * max(ess, 2.0) ** (-0.2)
    dens = float(w @ np.exp(-0.5 * ((x - qv) / h) ** 2)) / (h * math.sqrt(2 * math.pi)) if h > 0 else float("nan")
    se = math.sqrt(q * (1 - q) / max(ess, 1.0)) / dens if dens > 0 else float("inf")
    return qv, float(se), ess


_TAIL_TOL = 1e-12


def null_tail_test(x: Any, w: Optional[Any], t_obs: float, alpha: float, side: str) -> Dict[str, Any]:
    """The crossing test of an observed statistic against a weighted null that KEEPS +-inf (a bounded, unresolved
    replicate), unlike ``weighted_quantile``, which drops them.

    side "upper" (a's credit): crosses iff the null's weighted tail W(x >= t_obs) < alpha, i.e. t_obs > threshold with
    threshold = min{x_i : W(x > x_i) < alpha}; side "lower" (b's credit) the mirror (W(x <= t_obs) < alpha). The
    inequality is strict (the prototype's quantile crosses at a tail of exactly alpha), so an infinite mass >= alpha
    puts the threshold at +-inf and nothing crosses. Returns the threshold, its Monte Carlo SE (kernel density at the
    threshold over the finite part, Silverman bandwidth), the weights' ESS, the observed tail weight and its MC SE
    sqrt(p (1 - p) / ESS), the infinite mass, and ``crosses``. NaN in a weighted replicate is refused."""
    if side not in ("upper", "lower"):
        raise MeasurementModelRefusal("side must be 'upper' or 'lower'")
    x = np.asarray(x, dtype=float).ravel()
    w = np.ones_like(x) if w is None else np.asarray(w, dtype=float).ravel()
    if w.shape != x.shape:
        raise MeasurementModelRefusal("null weights and replicates differ in length")
    live = w > 0
    x, w = x[live], w[live]
    if x.size == 0:
        return {"threshold": float("nan"), "threshold_mc_se": float("nan"), "ess": 0.0, "tail_weight": float("nan"),
                "tail_weight_mc_se": float("nan"), "infinite_mass": float("nan"), "crosses": False}
    if np.any(np.isnan(x)):
        raise MeasurementModelRefusal("a null replicate is NaN: bound it (+-inf) or resolve it, never drop it")
    if not math.isfinite(float(t_obs)):
        raise MeasurementModelRefusal(f"the observed statistic is not finite ({t_obs!r}): it crosses nothing")
    w = w / w.sum()
    ess = float(1.0 / np.sum(w ** 2))
    s = 1.0 if side == "upper" else -1.0
    z, t = s * x, s * float(t_obs)                       # the upper-tail problem on z
    uv, inv = np.unique(z, return_inverse=True)
    wu = np.bincount(inv, weights=w, minlength=uv.size)
    above = np.concatenate([np.cumsum(wu[::-1])[::-1][1:], [0.0]])     # W(z > uv_j)
    j = int(np.argmax(above < alpha - _TAIL_TOL))
    thr_z = float(uv[j])
    tail = min(1.0, max(0.0, float(np.sum(w[z >= t]))))
    inf_mass = float(np.sum(w[np.isinf(z)]))
    crosses = bool(tail < alpha - _TAIL_TOL)
    fin = np.isfinite(z)
    se = float("nan")
    if math.isfinite(thr_z) and fin.sum() >= 2:
        wf = w[fin] / w[fin].sum()
        zf = z[fin]
        mu = float(wf @ zf)
        sd = float(np.sqrt(max(wf @ (zf - mu) ** 2, 1e-300)))
        h = 1.06 * sd * max(ess, 2.0) ** (-0.2)
        dens = (1.0 - inf_mass) * float(wf @ np.exp(-0.5 * ((zf - thr_z) / h) ** 2)) / (h * math.sqrt(2 * math.pi))
        q = 1.0 - alpha
        se = math.sqrt(q * (1 - q) / max(ess, 1.0)) / dens if dens > 0 else float("inf")
    return {"threshold": s * thr_z, "threshold_mc_se": se, "ess": ess, "tail_weight": tail,
            "tail_weight_mc_se": math.sqrt(tail * (1 - tail) / max(ess, 1.0)), "infinite_mass": inf_mass,
            "crosses": crosses}


def require_precision_ess(precision_ess: Any) -> float:
    """The density estimator's precision requirement, declared by the caller (e.g. ``tie_band_precision_ess``): no
    default, and 0 is refused (a silent 0 would accept any ESS >= 100 whatever the tie band needs)."""
    if precision_ess is None:
        raise MeasurementModelRefusal("declare precision_ess (tie_band_precision_ess); there is no default")
    v = float(precision_ess)
    if not (math.isfinite(v) and v > 0):
        raise MeasurementModelRefusal(f"precision_ess must be finite and > 0 (got {precision_ess!r}); a default of 0 "
                                      "is refused")
    return v


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    """The Wilson 95 % interval of a binomial rate k / n."""
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def misspecified_false_credit_gate(k: int, n: int, alpha: float = NOMINAL_ALPHA) -> Dict[str, Any]:
    """M1.1's power-check condition: false credit under the generator that adds the same-observable reproducibility
    term (else the lane's within-laboratory term); if the lower end of its Wilson 95 % interval exceeds ``alpha`` the
    case cannot credit.

    The check must be able to fail: an n at which an observed false-credit rate of 2 alpha (k = floor(2 alpha n))
    would still leave the Wilson lower end <= alpha is refused (``n_adequate`` False, ``passes`` False), so (0, 1)
    cannot pass the gate. ``min_n`` is the smallest n from which every n is adequate (100 at alpha = 0.05)."""
    k, n = int(k), int(n)
    lo, hi = wilson_interval(k, n)
    adequate = power_check_n_adequate(n, alpha)
    return {"k": k, "n": n, "rate": (k / n) if n else float("nan"), "wilson95": [lo, hi], "n_adequate": adequate,
            "min_n": power_check_min_n(alpha), "passes": bool(adequate and not (lo > alpha))}


def power_check_n_adequate(n: int, alpha: float = NOMINAL_ALPHA) -> bool:
    """Whether a power check of n sets could fail the gate at a false-credit rate of 2 alpha."""
    n = int(n)
    return bool(n > 0 and wilson_interval(int(math.floor(2.0 * alpha * n + 1e-9)), n)[0] > alpha)


def power_check_min_n(alpha: float = NOMINAL_ALPHA, n_max: int = 100000) -> int:
    """The smallest n from which every larger n is adequate (``power_check_n_adequate``)."""
    last_bad = 0
    for n in range(1, n_max + 1):
        if not power_check_n_adequate(n, alpha):
            last_bad = n
        elif n - last_bad > 200:
            break
    return last_bad + 1


# ---- the declared noise of a unit (item 5) ----------------------------------------------------------------------

@dataclass(frozen=True)
class NoiseTerm:
    """One declared noise term of a unit, in the unit's basis.

    ``sd`` is the estimate (a scalar, or one value per measured point); ``df`` its estimation degrees of freedom (the
    term is then integrated over its sampling uncertainty and re-drawn in the bootstrap); ``observable`` the observable
    it was ESTIMATED on (checked against the unit's: no transfer); ``is_bound`` marks a bound (an LOD, a worst case, a
    one-sided limit), which is refused: "no bound is used as a standard deviation"."""

    kind: str
    sd: Any
    observable: str
    structure: str = "per_point"
    df: Optional[float] = None
    is_bound: bool = False
    source: str = ""

    def __post_init__(self):
        if self.kind not in NOISE_TERM_KINDS:
            raise MeasurementModelRefusal(f"noise term kind must be one of {NOISE_TERM_KINDS}, got {self.kind!r}")
        if self.is_bound:
            raise MeasurementModelRefusal(
                f"{self.kind} term {self.source or self.observable!r} is a bound: no bound is used as a standard "
                "deviation (sec. 7 item 5)")
        if self.structure not in NOISE_STRUCTURES:
            raise MeasurementModelRefusal(f"structure must be one of {NOISE_STRUCTURES}, got {self.structure!r}")
        sd = np.atleast_1d(np.asarray(self.sd, dtype=float))
        if sd.ndim != 1 or not np.all(np.isfinite(sd)) or np.any(sd <= 0):
            raise MeasurementModelRefusal(f"{self.kind}: sd must be finite and > 0")
        if self.structure == "run_offset" and sd.size != 1:
            raise MeasurementModelRefusal(f"{self.kind}: a run offset has one sd")
        if self.df is not None and not (math.isfinite(float(self.df)) and float(self.df) > 0):
            raise MeasurementModelRefusal(f"{self.kind}: df must be > 0")
        if self.kind in _DF_REQUIRED and self.df is None:
            raise MeasurementModelRefusal(
                f"{self.kind}: declare its estimation df (it is integrated over its sampling / estimation uncertainty)")
        if self.kind == "instrument_floor" and (self.structure != "per_point" or self.df is not None):
            raise MeasurementModelRefusal("an instrument floor is a per-point floor with no df")
        if self.kind == "between_lab" and self.structure != "run_offset":
            raise MeasurementModelRefusal("a between-laboratory spread is a shared (run_offset) term")
        if self.kind in ("within_run_sd", "within_run_prior") and self.structure != "per_point":
            raise MeasurementModelRefusal(f"{self.kind} is a per-point term")

    @property
    def sd_array(self) -> np.ndarray:
        return np.atleast_1d(np.asarray(self.sd, dtype=float))

    def as_dict(self) -> Dict[str, Any]:
        sd = self.sd_array
        return {"kind": self.kind, "sd": float(sd[0]) if sd.size == 1 else sd.tolist(), "observable": self.observable,
                "structure": self.structure, "df": self.df, "source": self.source}


def printed_run_sd(sd: Any, observable: str, n: Optional[int] = None, replicate_level: Optional[str] = None,
                   source: str = "") -> NoiseTerm:
    """A printed SD as item 5 reads it: the SD of INDEPENDENT RUNS (``replicate_level == "independent_runs"``) is a
    run_sd term with n - 1 df, an unstated n taken as 2; an SD whose replicate level is unstated (or is a repeat of
    the measurement) counts as repeatability, a floor (an instrument_floor term)."""
    if replicate_level != "independent_runs":
        return NoiseTerm("instrument_floor", sd, observable, source=source or "printed SD, replicate level unstated")
    n_eff = 2 if n is None else int(n)
    if n_eff < 2:
        raise MeasurementModelRefusal("an SD of independent runs needs n >= 2")
    return NoiseTerm("run_sd", sd, observable, df=float(n_eff - 1), source=source)


@dataclass(frozen=True)
class UnitNoiseModel:
    """The declared noise model of one unit (sec. 7 item 5), shared by every rival.

    ``level`` is the level at which the compared value varies: "within_run" (a within-run ratio, shape or sign of a
    within-study contrast: a run-level term cancels), "absolute" (a level within one laboratory) or "cross_lab" (a
    comparison across laboratories). ``design`` (values x points) maps measured points to the scored values (identity
    by default). The covariance at variance node v is A diag(max(v_pp, floor^2)) A^T + v_off (A 1)(A 1)^T.

    Construction refuses: a term estimated on another observable (except the lane's within-laboratory term, which
    makes the unit falsification-only); an instrument floor as the whole noise of a reaction-system observable; a
    run-level (run_offset) term on a within-run value. ``role`` is then "credit", "falsification_only" (an absolute
    level without both a run-to-run term and a discrepancy estimate, or any unit carrying the lane term) or
    "not_decisive" (a cross-laboratory comparison with no measured between-laboratory spread)."""

    unit: str
    observable: str
    level: str
    terms: Tuple[NoiseTerm, ...]
    design: Optional[Any] = None
    reaction_system: bool = True

    def __post_init__(self):
        object.__setattr__(self, "terms", tuple(self.terms))
        if self.level not in VALUE_LEVELS:
            raise MeasurementModelRefusal(f"level must be one of {VALUE_LEVELS}, got {self.level!r}")
        if not self.terms:
            raise MeasurementModelRefusal(f"unit {self.unit!r}: no noise declared (zero is not assumed)")
        for t in self.terms:
            if not isinstance(t, NoiseTerm):
                raise MeasurementModelRefusal("terms must be NoiseTerm records")
            if t.kind != "lane_within_lab" and t.observable != self.observable:
                raise MeasurementModelRefusal(
                    f"unit {self.unit!r}: the {t.kind} term was estimated on {t.observable!r}, not on "
                    f"{self.observable!r}: no term estimated on another observable is transferred (sec. 7 item 5)")
        kinds = {t.kind for t in self.terms}
        if self.reaction_system and kinds <= {"instrument_floor"}:
            raise MeasurementModelRefusal(
                f"unit {self.unit!r}: an instrument's repeatability or a reading error is a floor, never the whole "
                "noise of a reaction-system observable (declare a run-to-run term, or for a within-run value a "
                "within_run_prior component bounded below by the floor)")
        if self.level == "within_run":
            if any(t.structure == "run_offset" for t in self.terms):
                raise MeasurementModelRefusal(f"unit {self.unit!r}: a run-level term cancels in a within-run value")
            if self.reaction_system and not (kinds & set(_PER_POINT_NOISE_KINDS)):
                raise MeasurementModelRefusal(
                    f"unit {self.unit!r}: a within-run value needs a per-point term beyond the floor")
        if self.design is not None:
            a = np.asarray(self.design, dtype=float)
            if a.ndim != 2:
                raise MeasurementModelRefusal("design must be (values x points)")
        n_nodes = int(np.prod([QUAD_K for t in self._integrated_terms()])) if self._integrated_terms() else 1
        if n_nodes > MAX_QUAD_NODES:
            raise MeasurementModelRefusal(f"unit {self.unit!r}: {n_nodes} quadrature nodes (> {MAX_QUAD_NODES})")

    # -- declaration ---------------------------------------------------------------------------------------------
    @property
    def role(self) -> str:
        kinds = {t.kind for t in self.terms}
        if "lane_within_lab" in kinds:
            return "falsification_only"
        if self.level == "within_run":
            return "credit"
        if self.level == "cross_lab" and "between_lab" not in kinds:
            return "not_decisive"
        if (kinds & set(_RUN_LEVEL_KINDS)) and "discrepancy" in kinds:
            return "credit"
        return "falsification_only"

    @property
    def redraws(self) -> FrozenSet[str]:
        """The null inputs this unit's noise brings: replicate SDs (printed run / within-run SDs) and the discrepancy
        component."""
        out = set()
        for t in self.terms:
            if t.kind in ("run_sd", "within_run_sd"):
                out.add("replicate_sds")
            if t.kind == "discrepancy":
                out.add("discrepancy_component")
        return frozenset(out)

    def as_dict(self) -> Dict[str, Any]:
        return {"unit": self.unit, "observable": self.observable, "level": self.level, "role": self.role,
                "terms": [t.as_dict() for t in self.terms],
                "design": None if self.design is None else np.asarray(self.design, float).tolist()}

    # -- the covariance --------------------------------------------------------------------------------------------
    def _integrated_terms(self) -> List[NoiseTerm]:
        return [t for t in self.terms if t.kind != "instrument_floor" and t.df is not None]

    def _matrix(self, n_values: int) -> np.ndarray:
        if self.design is None:
            return np.eye(int(n_values))
        a = np.asarray(self.design, dtype=float)
        if a.shape[0] != int(n_values):
            raise MeasurementModelRefusal(f"unit {self.unit!r}: design has {a.shape[0]} values, data {n_values}")
        return a

    def _point_vec(self, sd: np.ndarray, n_points: int) -> np.ndarray:
        if sd.size == 1:
            return np.full(n_points, float(sd[0]))
        if sd.size != n_points:
            raise MeasurementModelRefusal(f"unit {self.unit!r}: {sd.size} per-point sds for {n_points} points")
        return sd

    def variance_nodes(self, n_points: int) -> List[Tuple[np.ndarray, float]]:
        """The equal-weight quadrature nodes (per-point variance vector, run-offset variance) of the declared noise,
        each integrated term at its K inverse-chi-square quantiles, the floor applied per point."""
        floor2 = np.zeros(n_points)
        fixed_pp, fixed_off = np.zeros(n_points), 0.0
        integ: List[Tuple[NoiseTerm, np.ndarray]] = []
        for t in self.terms:
            v = self._point_vec(t.sd_array, n_points) ** 2 if t.structure == "per_point" else float(t.sd_array[0]) ** 2
            if t.kind == "instrument_floor":
                floor2 = np.maximum(floor2, v)
            elif t.df is None:
                if t.structure == "per_point":
                    fixed_pp = fixed_pp + v
                else:
                    fixed_off += v
            else:
                integ.append((t, v))
        mults = [invchi2_quadrature(t.df) for t, _ in integ]
        nodes = []
        for combo in (itertools.product(*mults) if mults else [()]):
            pp, off = fixed_pp.copy(), fixed_off
            for (t, v), m in zip(integ, combo):
                if t.structure == "per_point":
                    pp = pp + v * m
                else:
                    off += v * m
            nodes.append((np.maximum(pp, floor2), off))
        return nodes

    def log_lik(self, y: Any, f: Any) -> np.ndarray:
        """log p(y | f) under the declared noise (nats, the unit's basis), integrated over every term's estimation
        uncertainty. ``f`` is (S x values) predictions; ``y`` one data vector (returns S values) or (R x values)
        (returns R x S). A non-finite prediction gives -inf."""
        from scipy.linalg import solve_triangular

        f = np.atleast_2d(np.asarray(f, dtype=float))
        y = np.asarray(y, dtype=float)
        single = y.ndim == 1
        Y = np.atleast_2d(y)
        n = f.shape[1]
        if Y.shape[1] != n:
            raise MeasurementModelRefusal(f"unit {self.unit!r}: data have {Y.shape[1]} values, predictions {n}")
        a = self._matrix(n)
        a1 = a.sum(axis=1)
        r = Y[:, None, :] - f[None, :, :]
        bad = ~np.all(np.isfinite(r), axis=2)
        r = np.where(np.isfinite(r), r, 0.0).reshape(-1, n).T
        outs = []
        for pp, off in self.variance_nodes(a.shape[1]):
            cov = (a * pp[None, :]) @ a.T + off * np.outer(a1, a1)
            chol = np.linalg.cholesky(0.5 * (cov + cov.T))
            z = solve_triangular(chol, r, lower=True)
            q = np.sum(z ** 2, axis=0).reshape(Y.shape[0], f.shape[0])
            outs.append(-0.5 * q - float(np.sum(np.log(np.diag(chol)))) - 0.5 * n * math.log(2.0 * math.pi))
        st = np.stack(outs)
        m = np.max(st, axis=0)
        out = m + np.log(np.mean(np.exp(st - m), axis=0))
        out[bad] = -np.inf
        return out[0] if single else out

    def _draw_noise(self, n_rows: int, n_values: int, rng: np.random.Generator,
                    redraw_estimates: bool) -> Tuple[np.ndarray, Optional["UnitNoiseModel"]]:
        a = self._matrix(n_values)
        n_points = a.shape[1]
        floor2 = np.zeros(n_points)
        pp = np.zeros((n_rows, n_points))
        off = np.zeros(n_rows)
        new_terms = []
        for t in self.terms:
            v = self._point_vec(t.sd_array, n_points) ** 2 if t.structure == "per_point" else float(t.sd_array[0]) ** 2
            if t.kind == "instrument_floor":
                floor2 = np.maximum(floor2, v)
                new_terms.append(t)
                continue
            if t.df is None:
                mult = np.ones(n_rows)
                new_terms.append(t)
            else:
                mult = float(t.df) / rng.chisquare(float(t.df), n_rows)       # sigma^2 / s^2 from its posterior
                if redraw_estimates:
                    if n_rows != 1:
                        raise MeasurementModelRefusal("a bootstrap replicate is drawn one at a time")
                    est = float(mult[0]) * rng.chisquare(float(t.df)) / float(t.df)   # the estimate it prints
                    sd_new = math.sqrt(est) * t.sd_array
                    new_terms.append(replace(t, sd=float(sd_new[0]) if sd_new.size == 1 else sd_new.tolist()))
                else:
                    new_terms.append(t)
            if t.structure == "per_point":
                pp += mult[:, None] * np.atleast_1d(v)[None, :]
            else:
                off += mult * v
        pp = np.maximum(pp, floor2[None, :])
        e = rng.standard_normal((n_rows, n_points)) * np.sqrt(pp) + rng.standard_normal(n_rows)[:, None] * np.sqrt(off)[:, None]
        noise = e @ a.T
        return noise, (replace(self, terms=tuple(new_terms)) if redraw_estimates else None)

    def sample(self, f: Any, rng: np.random.Generator) -> np.ndarray:
        """Synthetic data at predictions ``f`` (R x values): each row's variances drawn from their estimation
        posteriors, then the noise (the predictive the scoring model integrates)."""
        f = np.atleast_2d(np.asarray(f, dtype=float))
        noise, _ = self._draw_noise(f.shape[0], f.shape[1], rng, False)
        return f + noise

    def bootstrap_draw(self, f: Any, rng: np.random.Generator) -> Tuple[np.ndarray, "UnitNoiseModel"]:
        """One bootstrap replicate at the truth ``f`` (values,): the data, and the noise model the replicate is scored
        under (every estimated term re-drawn: the replicate SDs and the discrepancy component)."""
        f = np.atleast_1d(np.asarray(f, dtype=float))
        noise, model = self._draw_noise(1, f.size, rng, True)
        return f + noise[0], model


# ---- forms ---------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FormDeclaration:
    """A form, declared with items 1 - 5 (M1.1): "experimental_fact" (an uncertain experimental fact: a discrete
    nuisance coordinate mixed in the density at ``weight``, equal unless ``weight_source`` cites one),
    "modelling_choice" (an unsourced modelling choice: keeps required agreement, each scored against its own null) or
    "omega_hedge" (a form that only hedged omega: dropped)."""

    name: str
    kind: str
    weight: Optional[float] = None
    weight_source: str = ""

    def __post_init__(self):
        if self.kind not in FORM_KINDS:
            raise MeasurementModelRefusal(f"form kind must be one of {FORM_KINDS}, got {self.kind!r}")
        if self.weight is not None and not (math.isfinite(float(self.weight)) and float(self.weight) > 0):
            raise MeasurementModelRefusal(f"form {self.name!r}: weight must be > 0")


BASE_FORM = "base"


def classify_forms(forms: Sequence[FormDeclaration]) -> Dict[str, Any]:
    """The forms as M1.1 classifies them: the experimental facts' normalised weights (equal unless sourced; unequal
    unsourced weights refused), the scored (required-agreement) forms, and the dropped omega hedges. With no
    experimental fact the mixture is the single component "base"; with no modelling choice the one scored form is
    "base"."""
    names = [f.name for f in forms]
    if len(set(names)) != len(names):
        raise MeasurementModelRefusal(f"duplicate form names {names}")
    facts = [f for f in forms if f.kind == "experimental_fact"]
    choices = [f.name for f in forms if f.kind == "modelling_choice"]
    dropped = [f.name for f in forms if f.kind == "omega_hedge"]
    if facts:
        raw = [1.0 if f.weight is None else float(f.weight) for f in facts]
        unsourced = [f.name for f in facts if not f.weight_source]
        if unsourced and len(set(raw)) > 1:
            raise MeasurementModelRefusal(
                f"unequal weights on experimental-fact forms without a source ({unsourced}): weights are equal unless "
                "sourced")
        tot = sum(raw)
        weights = {f.name: w / tot for f, w in zip(facts, raw)}
    else:
        weights = {BASE_FORM: 1.0}
    return {"fact_weights": weights, "scored_forms": choices or [BASE_FORM], "dropped_forms": dropped,
            "declared": [{"name": f.name, "kind": f.kind, "weight": f.weight, "weight_source": f.weight_source}
                         for f in forms]}


# ---- the measurement design of a case ------------------------------------------------------------------------------

@dataclass(frozen=True)
class MeasurementDesign:
    """Every scored unit's declared noise (decisive and spent), which units are spent (conditioning rows: they never
    enter T), the forms, per-experimental-fact noise overrides, and whether the case has forced inputs (measured
    courses that force the model: the null must redraw them)."""

    noise: Mapping[str, UnitNoiseModel]
    spent_units: Tuple[str, ...] = ()
    forms: Tuple[FormDeclaration, ...] = ()
    noise_by_fact_form: Mapping[str, Mapping[str, UnitNoiseModel]] = field(default_factory=dict)
    forced_inputs_present: bool = False

    def __post_init__(self):
        object.__setattr__(self, "spent_units", tuple(self.spent_units))
        object.__setattr__(self, "forms", tuple(self.forms))
        missing = sorted(set(self.spent_units) - set(self.noise))
        if missing:
            raise MeasurementModelRefusal(f"spent units without a noise model: {missing}")
        for u, m in self.noise.items():
            if m.unit != u:
                raise MeasurementModelRefusal(f"noise model keyed {u!r} declares unit {m.unit!r}")
        cls = classify_forms(self.forms)
        for k, over in self.noise_by_fact_form.items():
            if k not in cls["fact_weights"]:
                raise MeasurementModelRefusal(f"noise override for unknown experimental-fact form {k!r}")
            for u, m in over.items():
                base = self.noise.get(u)
                if base is None:
                    raise MeasurementModelRefusal(f"override for unknown unit {u!r}")
                sig = [(t.kind, t.structure, t.df is not None) for t in base.terms if t.kind != "instrument_floor"]
                sig2 = [(t.kind, t.structure, t.df is not None) for t in m.terms if t.kind != "instrument_floor"]
                if sig != sig2 or m.level != base.level or m.observable != base.observable:
                    raise MeasurementModelRefusal(
                        f"form {k!r} overrides unit {u!r} with a different structure: an experimental-fact form may "
                        "change fixed terms and floors only")

    @property
    def form_classes(self) -> Dict[str, Any]:
        return classify_forms(self.forms)

    def decisive_units(self) -> List[str]:
        return [u for u in self.noise if u not in self.spent_units]

    def units_by_role(self) -> Dict[str, List[str]]:
        out: Dict[str, List[str]] = {r: [] for r in UNIT_ROLES}
        for u in self.decisive_units():
            out[self.noise[u].role].append(u)
        return out

    def noise_for(self, fact_form: str, unit: str) -> UnitNoiseModel:
        return self.noise_by_fact_form.get(fact_form, {}).get(unit, self.noise[unit])

    def case_inputs(self) -> Dict[str, bool]:
        """Which of the null's listed inputs this case has (M1.1: the null must redraw every one present)."""
        redraws = set()
        for m in self.noise.values():
            redraws |= m.redraws
        return {"decisive_data": True, "conditioning_rows": bool(self.spent_units),
                "forced_inputs": bool(self.forced_inputs_present), "replicate_sds": "replicate_sds" in redraws,
                "discrepancy_component": "discrepancy_component" in redraws,
                "form_index": len(self.form_classes["fact_weights"]) > 1}


@dataclass(frozen=True)
class RivalDraws:
    """One rival's draws under one experimental-fact form: ``predictions[unit]`` is (S x values) for every scored unit
    (decisive and spent); ``log_prior_ratio`` log p(theta_s) / q(theta_s) when the draws are not from the prior."""

    predictions: Mapping[str, Any]
    log_prior_ratio: Optional[Any] = None

    @property
    def n(self) -> int:
        return int(np.atleast_2d(np.asarray(next(iter(self.predictions.values())))).shape[0])

    def pred(self, unit: str) -> np.ndarray:
        if unit not in self.predictions:
            raise MeasurementModelRefusal(f"no predictions for unit {unit!r}")
        p = np.asarray(self.predictions[unit], dtype=float)
        return p.reshape(self.n, -1)

    def lpr(self) -> np.ndarray:
        return np.zeros(self.n) if self.log_prior_ratio is None else np.asarray(self.log_prior_ratio, dtype=float)


#: A rival's ensemble per scored form: {scored form: {experimental-fact form: RivalDraws}}
RivalEnsemble = Mapping[str, Mapping[str, RivalDraws]]


# ---- the spent-conditioned predictive density -----------------------------------------------------------------------

def conditioned_lppd(log_lik_decisive: Any, log_lik_spent: Optional[Any] = None,
                     log_prior_ratio: Optional[Any] = None, *, precision_ess: float) -> DensityEstimate:
    """log p(y_dec | y_spent) from draws theta_s ~ q: log E_q[r L_dec L_spent] - log E_q[r L_spent] (r = p / q), the
    numerator and the denominator each through ``importance_lppd`` (PSIS) and each accepted only under the shared rule
    (k-hat < 0.7, ESS >= max(100, ``precision_ess``)). With no spent rows and draws from the prior the denominator is
    exactly 0. ``precision_ess`` has no default (``require_precision_ess``)."""
    precision_ess = require_precision_ess(precision_ess)
    ld = np.asarray(log_lik_decisive, dtype=float)
    s = ld.size
    lpr = np.zeros(s) if log_prior_ratio is None else np.asarray(log_prior_ratio, dtype=float)
    ls = np.zeros(s) if log_lik_spent is None else np.asarray(log_lik_spent, dtype=float)
    num = importance_lppd(lpr + ls + ld, np.zeros(s), precision_ess=precision_ess)
    if log_lik_spent is None and log_prior_ratio is None:
        den = DensityEstimate(0.0, "exact", True, "no spent rows, draws from the prior", n_draws=s)
    else:
        den = importance_lppd(lpr + ls, np.zeros(s), precision_ess=precision_ess)
    reasons = [f"numerator: {num.reason}"] * (not num.accepted) + [f"denominator: {den.reason}"] * (not den.accepted)
    return DensityEstimate(float(num.lppd - den.lppd), "psis_conditioned", bool(num.accepted and den.accepted),
                           "; ".join(reasons), k_hat=num.k_hat, ess=num.ess, ess_required=num.ess_required,
                           n_draws=s, n_nonfinite=num.n_nonfinite,
                           extra={"numerator": num.as_dict(), "denominator": den.as_dict()})


def _mix(per_fact: Mapping[str, Tuple[DensityEstimate, DensityEstimate]], weights: Mapping[str, float]
         ) -> DensityEstimate:
    """Mix (numerator, denominator) estimates over experimental-fact forms: log sum_k w_k N_k - log sum_k w_k D_k."""
    keys = list(per_fact)
    lw = np.array([math.log(weights[k]) for k in keys])
    num = _logsumexp(lw + np.array([per_fact[k][0].lppd for k in keys]))
    den = _logsumexp(lw + np.array([per_fact[k][1].lppd for k in keys]))
    bad = [f"{k} numerator: {per_fact[k][0].reason}" for k in keys if not per_fact[k][0].accepted]
    bad += [f"{k} denominator: {per_fact[k][1].reason}" for k in keys if not per_fact[k][1].accepted]
    ks = [e.k_hat for k in keys for e in per_fact[k] if e.k_hat is not None]
    es = [e.ess for k in keys for e in per_fact[k] if e.ess is not None and e.method != "exact"]
    return DensityEstimate(float(num - den), "psis_conditioned", not bad, "; ".join(bad),
                           k_hat=max(ks) if ks else None, ess=min(es) if es else None,
                           extra={"components": {k: {"numerator": per_fact[k][0].as_dict(),
                                                     "denominator": per_fact[k][1].as_dict()} for k in keys}})


def _score_form(y: Mapping[str, Any], design: MeasurementDesign, rival: Mapping[str, RivalDraws],
                weights: Mapping[str, float], credit: Sequence[str], adequacy: Sequence[str],
                noise: Optional[Mapping[Tuple[str, str], UnitNoiseModel]] = None,
                precision_ess: Optional[float] = None, rng: Optional[np.random.Generator] = None,
                ppc_reps: int = 0) -> Dict[str, Any]:
    """One rival, one scored form: the densities of T's keys and (with ``ppc_reps``) the adequacy p's with their
    Monte Carlo SE."""
    precision_ess = require_precision_ess(precision_ess)
    def nm(k, u):
        return noise[(k, u)] if noise is not None and (k, u) in noise else design.noise_for(k, u)

    missing = sorted(set(rival) ^ set(weights))
    if missing:
        raise MeasurementModelRefusal(f"the rival's experimental-fact components {sorted(rival)} differ from the "
                                      f"declared {sorted(weights)}")
    comp: Dict[str, Dict[str, Tuple[DensityEstimate, DensityEstimate]]] = {}
    keys = ["full"] + ([f"-{u}" for u in credit] + list(credit) if len(credit) > 1 else [])
    for key in keys:
        comp[key] = {}
    logpost, lls = {}, {}
    for k, dr in rival.items():
        lpr = dr.lpr()
        sp = np.zeros(dr.n)
        for u in design.spent_units:
            sp = sp + nm(k, u).log_lik(y[u], dr.pred(u))
        lu = {u: nm(k, u).log_lik(y[u], dr.pred(u)) for u in credit}
        exact_den = not design.spent_units and dr.log_prior_ratio is None
        base = lpr + sp
        den = (DensityEstimate(0.0, "exact", True, n_draws=dr.n) if exact_den
               else importance_lppd(base, np.zeros(dr.n), precision_ess=precision_ess))
        tot = base + sum(lu.values())
        comp["full"][k] = (importance_lppd(tot, np.zeros(dr.n), precision_ess=precision_ess), den)
        if len(credit) > 1:
            for u in credit:
                comp[f"-{u}"][k] = (importance_lppd(tot - lu[u], np.zeros(dr.n), precision_ess=precision_ess), den)
                comp[u][k] = (importance_lppd(base + lu[u], np.zeros(dr.n), precision_ess=precision_ess), den)
        logpost[k] = math.log(weights[k]) - math.log(dr.n) + base
        lls[k] = lu
    lppd = {key: _mix(comp[key], weights) for key in keys}
    out: Dict[str, Any] = {"lppd": lppd}
    if ppc_reps and adequacy:
        if rng is None:
            raise MeasurementModelRefusal("adequacy needs an rng")
        ks = list(rival)
        flat = np.concatenate([logpost[k] for k in ks])
        w = np.exp(flat - np.max(flat[np.isfinite(flat)])) if np.any(np.isfinite(flat)) else np.zeros_like(flat)
        w = np.where(np.isfinite(flat), w, 0.0)
        out["adequacy_p"], out["adequacy_mc_se"] = {}, {}
        if w.sum() <= 0:
            out["adequacy_p"] = {u: 0.0 for u in adequacy}
            out["adequacy_mc_se"] = {u: 0.0 for u in adequacy}
            return out
        w = w / w.sum()
        offs = np.cumsum([0] + [rival[k].n for k in ks])
        with np.errstate(divide="ignore"):
            lw_all = np.log(w)
        lw_k = {k: lw_all[offs[i]:offs[i + 1]] for i, k in enumerate(ks)}
        for u in adequacy:
            def neg_log_pred(Y):
                parts = [lw_k[k][None, :] + nm(k, u).log_lik(Y, rival[k].pred(u)) for k in ks]
                st = np.concatenate(parts, axis=1)
                m = np.max(st, axis=1, keepdims=True)
                m = np.where(np.isfinite(m), m, 0.0)
                return -(m[:, 0] + np.log(np.sum(np.exp(st - m), axis=1)))
            d_obs = float(neg_log_pred(np.atleast_2d(np.asarray(y[u], float)))[0])
            pick = rng.choice(w.size, size=int(ppc_reps), p=w)
            reps = np.empty((int(ppc_reps), np.asarray(y[u]).size))
            for i, k in enumerate(ks):
                sel = (pick >= offs[i]) & (pick < offs[i + 1])
                if sel.any():
                    reps[sel] = nm(k, u).sample(rival[k].pred(u)[pick[sel] - offs[i]], rng)
            d_rep = neg_log_pred(reps)
            p_u = float((1 + np.sum(d_rep >= d_obs)) / (int(ppc_reps) + 1))
            out["adequacy_p"][u] = p_u
            out["adequacy_mc_se"][u] = math.sqrt(p_u * (1.0 - p_u) / int(ppc_reps))
    return out


@dataclass(frozen=True)
class FormObservation:
    """One scored form's evidence on the observed data: ``lppd[rival][key]`` the density of each key of T ("full",
    "-u" and "u" per credit unit when there are several) and ``adequacy_p[rival][unit]`` each adequacy unit's
    posterior-predictive p under p(y_u | y_spent), with ``adequacy_mc_se`` its Monte Carlo SE (sqrt(p (1 - p) / R) for
    R predictive replicates; None when the caller's estimator did not report it). Rivals are keyed "a" and "b"."""

    lppd: Mapping[str, Mapping[str, DensityEstimate]]
    adequacy_p: Mapping[str, Mapping[str, float]]
    adequacy_mc_se: Optional[Mapping[str, Mapping[str, float]]] = None

    def T(self) -> Dict[str, float]:
        return {k: float(self.lppd["a"][k].lppd - self.lppd["b"][k].lppd) for k in self.lppd["a"]}


def score_from_draws(y: Mapping[str, Any], design: MeasurementDesign, draws_a: RivalEnsemble, draws_b: RivalEnsemble,
                     rng: np.random.Generator, ppc_reps: int = 400, *, precision_ess: float
                     ) -> Dict[str, FormObservation]:
    """The observed data's evidence per scored form: every density and every adequacy p the verdict reads."""
    cls = design.form_classes
    roles = design.units_by_role()
    credit, adequacy = roles["credit"], roles["credit"] + roles["falsification_only"]
    out = {}
    for f in cls["scored_forms"]:
        for name, d in (("a", draws_a), ("b", draws_b)):
            if f not in d:
                raise MeasurementModelRefusal(f"rival {name} has no draws for scored form {f!r}")
        sa = _score_form(y, design, draws_a[f], cls["fact_weights"], credit, adequacy, None, precision_ess, rng,
                         ppc_reps)
        sb = _score_form(y, design, draws_b[f], cls["fact_weights"], credit, adequacy, None, precision_ess, rng,
                         ppc_reps)
        out[f] = FormObservation({"a": sa["lppd"], "b": sb["lppd"]},
                                 {"a": sa.get("adequacy_p", {}), "b": sb.get("adequacy_p", {})},
                                 {"a": sa.get("adequacy_mc_se", {}), "b": sb.get("adequacy_mc_se", {})})
    return out


# ---- the null: a parametric bootstrap under one rival as truth ------------------------------------------------------

@dataclass(frozen=True)
class NullDistribution:
    """The null of T (and of every key the verdict reads) under one rival as truth: ``stats[form][key]`` replicate
    arrays, optional ``weights`` (a re-weighted pool), ``resolved`` per replicate, the inputs it ``redrew`` (from
    NULL_INPUTS) and whether its truths came from OUTSIDE the scoring ensemble."""

    truth: str
    stats: Mapping[str, Mapping[str, Any]]
    redrawn: FrozenSet[str]
    truths_outside_scoring_ensemble: bool
    resolved: Optional[Any] = None
    weights: Optional[Any] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_replicates(self) -> int:
        first = next(iter(next(iter(self.stats.values())).values()))
        return int(np.asarray(first).size)

    @property
    def n_effective(self) -> float:
        if self.weights is None:
            return float(self.n_replicates)
        w = np.asarray(self.weights, dtype=float)
        return float(w.sum() ** 2 / np.sum(w ** 2)) if w.sum() > 0 else 0.0

    def quantile(self, form: str, key: str, q: float) -> Tuple[float, float, float]:
        """The prototype's quantile of the RESOLVED, finite replicates (unresolved ones dropped): a diagnostic only;
        the verdict reads ``bounded`` through ``null_tail_test``."""
        x = np.asarray(self.stats[form][key], float)
        w = np.ones(x.size) if self.weights is None else np.asarray(self.weights, float)
        return weighted_quantile(x, np.where(self.resolved_mask(), w, 0.0), q)

    def resolved_mask(self) -> np.ndarray:
        n = self.n_replicates
        return np.ones(n, bool) if self.resolved is None else np.asarray(self.resolved, bool).reshape(n)

    @property
    def unresolved_fraction(self) -> float:
        """The weighted fraction of replicates with an unresolved density (the share the bound sets to +-inf)."""
        n = self.n_replicates
        w = np.ones(n) if self.weights is None else np.asarray(self.weights, float)
        tot = float(w.sum())
        return float(np.sum(w[~self.resolved_mask()]) / tot) if tot > 0 else float("nan")

    def bounded(self, form: str, key: str, sign: float) -> np.ndarray:
        """The replicates of ``key`` with every unresolved (or non-finite) one set to ``sign`` x inf: never dropped,
        never read as estimated (M1.1, PR #105)."""
        x = np.asarray(self.stats[form][key], float)
        return np.where(self.resolved_mask() & np.isfinite(x), x, math.copysign(math.inf, sign))


@dataclass(frozen=True)
class NullPool:
    """Bootstrap replicates under one rival as truth with truths drawn from its prior (x form weights); a synthetic
    or observed set's null is the pool re-weighted by p(y_spent | truth) (``null_from_pool``), the prototype's pool."""

    truth: str
    stats: Mapping[str, Mapping[str, np.ndarray]]
    resolved: np.ndarray
    redrawn: FrozenSet[str]
    truths_outside_scoring_ensemble: bool
    truth_form: np.ndarray
    truth_index: np.ndarray
    truths: Mapping[str, RivalDraws]
    design: MeasurementDesign
    extra: Dict[str, Any] = field(default_factory=dict)


def _same_draws(x: RivalDraws, y: RivalDraws) -> bool:
    """Whether truth draws ``x`` overlap the scoring draws ``y``: a shared buffer, an identical array, or ANY draw of
    ``x`` whose predictions on every common unit equal one draw of ``y`` (truths copied row by row, or a subset).
    Rival draws whose predictions do not depend on theta at all are flagged too (conservative: they have no truth
    outside the ensemble to draw)."""
    common = [u for u in x.predictions if u in y.predictions]
    for u in common:
        a, b = np.asarray(x.predictions[u]), np.asarray(y.predictions[u])
        if np.shares_memory(a, b) or (a.shape == b.shape and np.array_equal(a, b)):
            return True
    cols = [u for u in common if x.pred(u).shape[1] == y.pred(u).shape[1]]
    if not cols:
        return False

    def rows(d: RivalDraws) -> List[bytes]:
        m = np.ascontiguousarray(np.hstack([d.pred(u) for u in cols]) + 0.0, dtype=np.float64)   # -0.0 -> 0.0
        return [r.tobytes() for r in m]

    seen = set(rows(y))
    return any(r in seen for r in rows(x))


def build_null_pool(truth: str, truths: Mapping[str, RivalDraws], design: MeasurementDesign, draws_a: RivalEnsemble,
                    draws_b: RivalEnsemble, n: int, rng: np.random.Generator, *, precision_ess: float,
                    redraw_forced: Optional[Callable[[np.random.Generator, str, int],
                                                     Tuple[RivalEnsemble, RivalEnsemble, Mapping[str, np.ndarray]]]]
                    = None, retry_draws: Optional[Tuple[RivalEnsemble, RivalEnsemble]] = None) -> NullPool:
    """The parametric bootstrap of T under rival ``truth`` ("a" or "b") as truth.

    ``truths`` holds the truth rival's HELD-OUT draws per experimental-fact form (outside the scoring ensemble:
    checked against both ensembles). Each replicate draws the form index at the declared weights, a truth draw from
    that form (by its prior ratio when not from the prior), the decisive data and the conditioning rows at that truth
    under the form's declared noise, with every estimated term re-drawn (the replicate SDs and the discrepancy
    component); its T (every key) is scored exactly as the observed data. Forced inputs are redrawn only through
    ``redraw_forced(rng, form, index) -> (draws_a, draws_b, truth predictions {unit: values})``; without it a case that
    declares forced inputs gets a pool that does not redraw them (an incomplete null: the verdict is "not
    computable").

    Re-estimation (M1.1, PR #105): a replicate with any unresolved density is re-scored ONCE, on the same replicate
    data and re-drawn noise, with ``retry_draws = (draws_a, draws_b)``: a larger, independently drawn ensemble per rival
    (also checked to lie outside the truths). Limitation: the draws are the caller's, so without ``retry_draws`` there
    is no re-estimation, and a replicate whose forced inputs were redrawn (``redraw_forced``) is not retried (its
    ensembles are the replicate's own). A replicate still unresolved stays ``resolved = False`` and the verdict bounds
    it (``judge_pair``). ``extra`` records the counts. Not enforced: ``retry_draws`` is checked against the truths but
    not against the scoring draws, and "larger" is the caller's responsibility (its size is not checked)."""
    precision_ess = require_precision_ess(precision_ess)
    if truth not in ("a", "b"):
        raise MeasurementModelRefusal("truth must be 'a' or 'b'")
    cls = design.form_classes
    weights = cls["fact_weights"]
    if set(truths) != set(weights):
        raise MeasurementModelRefusal(f"truth draws {sorted(truths)} differ from the experimental-fact forms "
                                      f"{sorted(weights)}")
    ensembles = (draws_a, draws_b) + (tuple(retry_draws) if retry_draws is not None else ())
    outside = not any(_same_draws(truths[k], ens[f][k]) for ens in ensembles
                      for f in cls["scored_forms"] for k in weights if f in ens and k in ens[f])
    roles = design.units_by_role()
    credit = roles["credit"]
    units = design.decisive_units() + list(design.spent_units)
    ks = list(weights)
    pk = np.array([weights[k] for k in ks])
    stats: Dict[str, Dict[str, List[float]]] = {}
    resolved = np.zeros(int(n), bool)
    tf, ti = np.empty(int(n), dtype=object), np.zeros(int(n), int)
    n_retried = n_rescued = 0

    def score_replicate(y, rep_noise, da, db):
        vals: Dict[str, Dict[str, float]] = {}
        ok = True
        for f in cls["scored_forms"]:
            sa = _score_form(y, design, da[f], weights, credit, [], rep_noise, precision_ess)
            sb = _score_form(y, design, db[f], weights, credit, [], rep_noise, precision_ess)
            vals[f] = {}
            for key in sa["lppd"]:
                vals[f][key] = float(sa["lppd"][key].lppd - sb["lppd"][key].lppd)
                ok = ok and sa["lppd"][key].accepted and sb["lppd"][key].accepted
        return vals, ok

    for i in range(int(n)):
        k = ks[int(rng.choice(len(ks), p=pk))]
        td = truths[k]
        tw = td.lpr()
        p = np.exp(tw - np.max(tw))
        j = int(rng.choice(td.n, p=p / p.sum()))
        tf[i], ti[i] = k, j
        da, db = draws_a, draws_b
        tpred = {u: td.pred(u)[j] for u in units}
        if redraw_forced is not None and design.forced_inputs_present:
            da, db, tp = redraw_forced(rng, k, j)
            tpred.update({u: np.asarray(v, float) for u, v in tp.items()})
        y, rep_noise = {}, {}
        for u in units:
            yu, m = design.noise_for(k, u).bootstrap_draw(tpred[u], rng)
            y[u] = yu
            for kk in ks:          # the replicate's re-drawn estimates hold for every experimental-fact component
                base = design.noise_for(kk, u)
                rep_noise[(kk, u)] = replace(base, terms=tuple(_align_terms(base, m)))
        vals, ok = score_replicate(y, rep_noise, da, db)
        forced_here = redraw_forced is not None and design.forced_inputs_present
        if not ok and retry_draws is not None and not forced_here:
            n_retried += 1
            vals2, ok2 = score_replicate(y, rep_noise, retry_draws[0], retry_draws[1])
            if ok2:
                vals, ok, n_rescued = vals2, True, n_rescued + 1
        for f, d in vals.items():
            sf = stats.setdefault(f, {})
            for key, v in d.items():
                sf.setdefault(key, []).append(v)
        resolved[i] = ok
    redrawn = {"decisive_data", "replicate_sds", "discrepancy_component", "form_index"}
    if design.spent_units:
        redrawn.add("conditioning_rows")
    if design.forced_inputs_present and redraw_forced is not None:
        redrawn.add("forced_inputs")
    extra = {"reestimation": "retry_draws" if retry_draws is not None else "none (no retry draw set supplied)",
             "n_retried": int(n_retried), "n_rescued_by_retry": int(n_rescued),
             "n_unresolved_after_retry": int(np.sum(~resolved))}
    return NullPool(truth, {f: {k: np.asarray(v) for k, v in d.items()} for f, d in stats.items()}, resolved,
                    frozenset(redrawn), bool(outside), tf, ti, truths, design, extra)


def _align_terms(base: UnitNoiseModel, redrawn: UnitNoiseModel) -> List[NoiseTerm]:
    """The re-drawn model's terms matched to ``base``'s by position among the estimated (df-carrying) terms."""
    est = [t for t in redrawn.terms if t.kind != "instrument_floor" and t.df is not None]
    out, i = [], 0
    for t in base.terms:
        if t.kind != "instrument_floor" and t.df is not None:
            out.append(replace(t, sd=est[i].sd))
            i += 1
        else:
            out.append(t)
    return out


def null_from_pool(pool: NullPool, y_spent: Optional[Mapping[str, Any]] = None) -> NullDistribution:
    """The null for one data set: the pool re-weighted by p(y_spent | truth_j) under the declared noise (the truth
    rival's spent-conditioned posterior). With no spent rows the replicates are unweighted. Records the weights' ESS
    and PSIS k-hat."""
    d = pool.design
    extra: Dict[str, Any] = {"pool_size": int(pool.resolved.size),
                             "pool_resolved_fraction": float(np.mean(pool.resolved)) if pool.resolved.size else 0.0,
                             **pool.extra}
    if not d.spent_units:
        return NullDistribution(pool.truth, pool.stats, pool.redrawn, pool.truths_outside_scoring_ensemble,
                                pool.resolved, None, extra)
    if y_spent is None:
        raise MeasurementModelRefusal("a case with spent rows needs y_spent to condition its null")
    lw = np.zeros(pool.resolved.size)
    for k in set(pool.truth_form.tolist()):
        sel = pool.truth_form == k
        td = pool.truths[k]
        for u in d.spent_units:
            lw[sel] += d.noise_for(k, u).log_lik(y_spent[u], td.pred(u)[pool.truth_index[sel]])
    fin = np.isfinite(lw)
    w = np.zeros_like(lw)
    if fin.any():
        w[fin] = np.exp(lw[fin] - np.max(lw[fin]))
        extra["weights_k_hat"] = float(psis(lw[fin] - np.max(lw[fin]))[1]) if fin.sum() >= 10 else float("inf")
    else:
        extra["weights_k_hat"] = float("inf")
    return NullDistribution(pool.truth, pool.stats, pool.redrawn, pool.truths_outside_scoring_ensemble,
                            pool.resolved, w, extra)


# ---- the verdict ---------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class RevisedVerdict:
    """The revised M1.1 verdict of one pair, with every intermediate quantity a record needs."""

    outcome: str
    credited: Optional[str]
    carried: Tuple[str, ...]
    reasons: Tuple[str, ...]
    details: Dict[str, Any]

    def as_dict(self) -> Dict[str, Any]:
        return {"outcome": self.outcome, "credited": self.credited, "carried": list(self.carried),
                "reasons": list(self.reasons), **self.details}


def _null_problems(null: Optional[NullDistribution], name: str, forms: Sequence[str], keys: Sequence[str],
                   case_inputs: Mapping[str, bool], min_replicates: int) -> List[str]:
    if null is None:
        return [f"no null under {name} as truth"]
    out = []
    need = {k for k, v in case_inputs.items() if v}
    missing = sorted(need - set(null.redrawn))
    if missing:
        out.append(f"null under {name} does not redraw {missing}: not the rule's null")
    if not null.truths_outside_scoring_ensemble:
        out.append(f"null under {name}: truths not drawn outside the scoring ensemble")
    for f in forms:
        if f not in null.stats:
            out.append(f"null under {name} has no replicates for form {f!r}")
            continue
        miss = sorted(set(keys) - set(null.stats[f]))
        if miss:
            out.append(f"null under {name}, form {f!r}: no replicates of {miss}")
    if null.n_effective < min_replicates:
        out.append(f"null under {name}: {null.n_effective:.0f} effective replicates < {min_replicates}")
    k_w = null.extra.get("weights_k_hat")
    if null.weights is not None and not (k_w is not None and k_w < K_HAT_MAX):
        out.append(f"null under {name}: re-weighting k-hat {k_w} >= {K_HAT_MAX}")
    if null.resolved is not None and np.asarray(null.resolved).size != null.n_replicates:
        out.append(f"null under {name}: {np.asarray(null.resolved).size} resolution flags for {null.n_replicates} "
                   "replicates")
    # an unresolved replicate no longer makes the null incomplete: the verdict bounds it (PR #105)
    return out


def judge_pair(observed: Mapping[str, FormObservation], null_a: Optional[NullDistribution],
               null_b: Optional[NullDistribution], units: Mapping[str, UnitNoiseModel],
               case_inputs: Mapping[str, bool], forms: Sequence[FormDeclaration] = (),
               falsified: Optional[Mapping[str, bool]] = None, incumbent: Optional[str] = None,
               misspecified_false_credit: Optional[Tuple[int, int]] = None, alpha: float = NOMINAL_ALPHA,
               min_replicates: int = NULL_REPLICATES_MIN, names: Tuple[str, str] = ("a", "b")) -> RevisedVerdict:
    """The revised M1.1 decision between rivals a and b (sec. 7 item 5 and "The decision statistic (M1.1)").

    ``observed[form]``: each scored form's densities and adequacy p's (``FormObservation``), from
    ``score_from_draws`` or any estimator the caller ran (each density accepted or not by the shared rule).
    ``null_a`` / ``null_b``: the null of every key of T under a / b as truth. ``units``: the declared noise models of the
    decisive units (decide each unit's role). ``case_inputs``: for EVERY input in NULL_INPUTS, whether this case has it
    (the null must redraw each present one). ``falsified``: S3 falsification of a / b. ``misspecified_false_credit``:
    the power check's (k, n) false credit under the reproducibility generator: REQUIRED for a credit (without it, or
    failing, a qualifying rival is not credited).

    Order: S3 falsification first ("both_falsified", or "carried_sole_survivor": the survivor is carried, not
    credited), so a falsified rival is never returned as carried; forms and roles (no credit unit: "not_decisive");
    resolution of the observed densities and the null's completeness ("not computable", carried undecided); then a
    qualifies iff in EVERY scored form T crosses its null's upper 5 % tail under b (strictly: W(T_null >= T) < alpha),
    a is adequate on every adequacy unit after Holm, and every leave-one-unit-out reading stands; b by the mirror rule
    under a. Unresolved bootstrap replicates are bounded: a is credited iff a qualifies and b does not with them at
    +inf in both nulls; b iff b qualifies and a does not with them at -inf. Both qualify (each under the bound against
    it): "neither_credited_both_qualify"; one credit: "credited" only with the power-check gate supplied and passing;
    otherwise "undecided". Each null's weighted unresolved fraction is reported."""
    a_name, b_name = names
    details: Dict[str, Any] = {"rivals": {"a": a_name, "b": b_name}, "alpha": alpha}
    missing_inputs = sorted(set(NULL_INPUTS) - set(case_inputs))
    if missing_inputs:
        raise MeasurementModelRefusal(f"case_inputs must state every listed null input; missing {missing_inputs}")
    if not case_inputs.get("decisive_data"):
        raise MeasurementModelRefusal("every case has decisive data")
    cls = classify_forms(forms)
    details["forms"] = cls
    if sorted(observed) != sorted(cls["scored_forms"]):
        raise MeasurementModelRefusal(f"observed forms {sorted(observed)} differ from the scored forms "
                                      f"{cls['scored_forms']} (omega hedges are dropped, facts are mixed)")
    roles: Dict[str, List[str]] = {r: [] for r in UNIT_ROLES}
    for u, m in units.items():
        roles[m.role].append(u)
    credit, adequacy = roles["credit"], roles["credit"] + roles["falsification_only"]
    details.update({"unit_roles": roles, "noise": {u: m.as_dict() for u, m in units.items()},
                    "case_inputs": dict(case_inputs)})
    fals = {"a": bool((falsified or {}).get("a")), "b": bool((falsified or {}).get("b"))}
    details["falsified"] = fals
    if incumbent is not None and incumbent not in ("a", "b"):
        raise MeasurementModelRefusal("incumbent must be 'a', 'b' or None")
    details["incumbent"] = incumbent

    reestimation: List[str] = []          # filled once the nulls are read: every later verdict records it

    def verdict(outcome, credited=None, carried=(), reasons=()):
        return RevisedVerdict(outcome, credited, tuple(carried), tuple(reasons) + tuple(reestimation), details)

    # S3 falsification first (D1): a falsified rival is never returned as carried, whatever the units' roles
    if fals["a"] and fals["b"]:
        return verdict("both_falsified", reasons=["both rivals falsified under S3: neither credited nor a null"])
    if fals["a"] or fals["b"]:
        surv, dead = ("b", "a") if fals["a"] else ("a", "b")
        rs = [f"{dead} falsified under S3: neither credited nor used as a null; {surv} has no unfalsified rival to "
              "serve as its null and is carried, not credited"]
        if incumbent == dead:
            rs.append("the falsified rival is the incumbent: S3's refusal clause applies to the outputs that depended "
                      "on it")
        details["s3_refusal_applies"] = incumbent == dead
        return verdict("carried_sole_survivor", carried=(surv,), reasons=rs)
    if not credit:
        return verdict("not_decisive", carried=("a", "b"),
                       reasons=["no unit can credit under item 5: the data cannot decide"])

    keys = ["full"] + ([f"-{u}" for u in credit] + list(credit) if len(credit) > 1 else [])
    problems = []
    dens_record: Dict[str, Any] = {}
    for f, ob in observed.items():
        for r in ("a", "b"):
            for key in keys:
                if key not in ob.lppd.get(r, {}):
                    problems.append(f"form {f!r}, rival {r}: no density for {key!r}")
                    continue
                e = ob.lppd[r][key]
                if not e.accepted:
                    problems.append(f"form {f!r}, rival {r}, {key}: unresolved density ({e.reason})")
            miss = sorted(set(adequacy) - set(ob.adequacy_p.get(r, {})))
            if miss:
                problems.append(f"form {f!r}, rival {r}: no adequacy p for {miss}")
        dens_record[f] = {r: {k: ob.lppd[r][k].as_dict() for k in keys if k in ob.lppd.get(r, {})} for r in ("a", "b")}
    details["densities"] = dens_record
    for f, ob in observed.items():                 # a non-finite T would cross both sides (no weight is >= NaN)
        if not any("unresolved density" in pr or "no density" in pr for pr in problems if f"form {f!r}" in pr):
            bad = sorted(k for k, v in ob.T().items() if k in keys and not math.isfinite(v))
            if bad:
                problems.append(f"form {f!r}: non-finite observed T for {bad} (never read as a crossing)")
    problems += _null_problems(null_a, "a", cls["scored_forms"], keys, case_inputs, min_replicates)
    problems += _null_problems(null_b, "b", cls["scored_forms"], keys, case_inputs, min_replicates)
    for nm_, nl in (("a", null_a), ("b", null_b)):
        if nl is not None:
            details[f"null_{nm_}"] = {"n_replicates": nl.n_replicates, "n_effective": nl.n_effective,
                                      "redrawn": sorted(nl.redrawn),
                                      "truths_outside_scoring_ensemble": nl.truths_outside_scoring_ensemble,
                                      **nl.extra}
            if nl.resolved is None or np.asarray(nl.resolved).size == nl.n_replicates:
                details[f"null_{nm_}"].update({
                    "unresolved_fraction": nl.unresolved_fraction,
                    "unresolved_count": int(np.sum(~nl.resolved_mask()))})
            ex = nl.extra
            if str(ex.get("reestimation", "none")).startswith("retry_draws"):
                reestimation.append(f"null under {nm_}: re-estimation: {int(ex.get('n_retried', 0))} retried "
                                    f"({int(ex.get('n_rescued_by_retry', 0))} rescued)")
            else:
                reestimation.append(f"null under {nm_}: re-estimation: none")
    if problems:
        return verdict("not_computable", carried=("a", "b"),
                       reasons=["a verdict needs resolved densities and a complete null"] + problems)

    def weights(nl):
        return None if nl.weights is None else np.asarray(nl.weights, float)

    def conditions(sign: float) -> Tuple[Dict[str, Any], Dict[str, bool]]:
        """Every scored form's conditions with each null's unresolved replicates at ``sign`` x inf."""
        per_form: Dict[str, Any] = {}
        qual = {"a": True, "b": True}
        for f in cls["scored_forms"]:
            ob = observed[f]
            T = ob.T()
            rec: Dict[str, Any] = {"T": T["full"], "unit_contributions": {u: T[u] for u in credit}
                                   if len(credit) > 1 else {credit[0]: T["full"]}}
            for side, null, q, tail in (("a", null_b, 1.0 - alpha, "upper"), ("b", null_a, alpha, "lower")):
                def test(key, t):
                    return null_tail_test(null.bounded(f, key, sign), weights(null), t, alpha, tail)
                full = test("full", T["full"])
                crosses = full["crosses"]
                hp = holm({u: float(ob.adequacy_p[side][u]) for u in adequacy}, alpha)
                adequate = all(v["p_holm"] >= alpha for v in hp.values())
                se_ad = (ob.adequacy_mc_se or {}).get(side, {})
                loo: Dict[str, Any] = {}
                loo_ok = True
                for u in credit:
                    if len(credit) == 1:
                        loo[u] = {"rests_on": True, "T_u": T["full"], "threshold_u": full["threshold"],
                                  "alone_credits": crosses}
                        loo_ok = loo_ok and crosses
                        continue
                    tm, tu = test(f"-{u}", T[f"-{u}"]), test(u, T[u])
                    rests = not tm["crosses"]
                    alone = tu["crosses"]
                    loo[u] = {"T_minus_u": T[f"-{u}"], "threshold_minus_u": tm["threshold"], "rests_on": bool(rests),
                              "T_u": T[u], "threshold_u": tu["threshold"], "alone_credits": bool(alone)}
                    if rests and not alone:
                        loo_ok = False
                ok = crosses and adequate and loo_ok
                rec[side] = {"threshold": full["threshold"], "threshold_mc_se": full["threshold_mc_se"],
                             "threshold_ess": full["ess"], "threshold_quantile": q,
                             "null_tail_weight": full["tail_weight"], "null_tail_weight_mc_se": full["tail_weight_mc_se"],
                             "null_infinite_mass": full["infinite_mass"], "crosses": crosses,
                             "adequacy": {u: {"p": v["p"], "p_holm": v["p_holm"], "p_mc_se": se_ad.get(u)}
                                          for u, v in hp.items()},
                             "adequate": adequate, "leave_one_unit_out": loo, "loo_stands": loo_ok, "qualifies": ok}
                qual[side] = qual[side] and ok
            per_form[f] = rec
        return per_form, qual

    # testing a's credit: unresolved at +inf in both nulls; testing b's: at -inf (identical when none is unresolved)
    pf_a, q_plus = conditions(+1.0)
    pf_b, q_minus = conditions(-1.0)
    per_form = {f: {"T": pf_a[f]["T"], "unit_contributions": pf_a[f]["unit_contributions"], "a": pf_a[f]["a"],
                    "b": pf_b[f]["b"]} for f in cls["scored_forms"]}
    details["per_form"] = per_form
    details["qualifies"] = {"a": q_plus["a"], "b": q_minus["b"]}
    details["qualifies_under_bound"] = {"testing_a_unresolved_plus_inf": dict(q_plus),
                                        "testing_b_unresolved_minus_inf": dict(q_minus)}
    details["bound_differs"] = q_plus != q_minus
    credit_ok = {"a": q_plus["a"] and not q_plus["b"], "b": q_minus["b"] and not q_minus["a"]}
    gate = None
    if misspecified_false_credit is not None:
        gate = misspecified_false_credit_gate(*misspecified_false_credit, alpha=alpha)
    details["power_check_gate"] = gate if gate is not None else "not supplied"
    if credit_ok["a"] or credit_ok["b"]:
        side = "a" if credit_ok["a"] else "b"
        if gate is None:
            return verdict("undecided", carried=("a", "b"),
                           reasons=[f"{side} qualifies but the case cannot credit: a credit requires the power check "
                                    "(false credit under the reproducibility generator, misspecified_false_credit) "
                                    "supplied and passing; it was not supplied"])
        if not gate["n_adequate"]:
            return verdict("undecided", carried=("a", "b"),
                           reasons=[f"{side} qualifies but the case cannot credit: the power check's n = {gate['n']} "
                                    f"could not fail the gate at a false-credit rate of 2 x {alpha} (needs n >= "
                                    f"{gate['min_n']})"])
        if not gate["passes"]:
            return verdict("undecided", carried=("a", "b"),
                           reasons=[f"{side} qualifies but the case cannot credit: the power check's false credit "
                                    f"under the reproducibility generator has Wilson lower end "
                                    f"{gate['wilson95'][0]:.3f} > {alpha}"])
        return verdict("credited", credited=side,
                       reasons=[f"{side} credited in every scored form ({', '.join(cls['scored_forms'])})"])
    if q_plus["a"] and q_minus["b"]:
        return verdict("neither_credited_both_qualify", carried=("a", "b"),
                       reasons=["both rivals meet their conditions: neither is credited"])
    rs = ["neither rival meets its conditions: both carried (S3)"]
    for side, fav in (("a", q_minus), ("b", q_plus)):
        if fav[side] and not fav["b" if side == "a" else "a"]:
            rs.append(f"{side} would be credited were the unresolved replicates read in its favour: the credit does "
                      "not survive the bound")
    return verdict("undecided", carried=("a", "b"), reasons=rs)


def judge_from_draws(y: Mapping[str, Any], design: MeasurementDesign, draws_a: RivalEnsemble, draws_b: RivalEnsemble,
                     pool_a: Optional[NullPool], pool_b: Optional[NullPool], rng: np.random.Generator,
                     ppc_reps: int = 400, *, precision_ess: float, **judge_kwargs) -> RevisedVerdict:
    """``score_from_draws`` + ``null_from_pool`` for both rivals + ``judge_pair``. ``pool_a`` / ``pool_b`` are the
    bootstrap pools under a / b as truth (``build_null_pool``; reuse them across synthetic sets in a power check)."""
    obs = score_from_draws(y, design, draws_a, draws_b, rng, ppc_reps, precision_ess=precision_ess)
    y_spent = {u: y[u] for u in design.spent_units}
    na = null_from_pool(pool_a, y_spent) if pool_a is not None else None
    nb = null_from_pool(pool_b, y_spent) if pool_b is not None else None
    units = {u: design.noise[u] for u in design.decisive_units()}
    return judge_pair(obs, na, nb, units, design.case_inputs(), design.forms, **judge_kwargs)
