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
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

# ----------------------------------------------------------------------------------------------------------------
# Constants and the dex <-> ln conversion
# ----------------------------------------------------------------------------------------------------------------

LN10 = math.log(10.0)
#: H0 (EX-059, reset set): the within-laboratory model discrepancy, a fractional error per predicted value, in dex.
OMEGA_DEX = 0.390
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
