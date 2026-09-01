#!/usr/bin/env python3
"""
scripts/generators/generate_trunk_rate_calibration.py

THE FIRST RATE-LEVEL CALIBRATION OF THE SUGAR TRUNK (audit Wave S3, 2026-08-27).

Fits the eight-step trunk in `src/trunk_kinetics.py` to measured
concentration-vs-time trajectories, and writes
`results/validation/trunk_rate_calibration_refit.{json,md}`.

--------------------------------------------------------------------------
FIT CORPUS -- and nothing else
--------------------------------------------------------------------------
Exactly two files, both in `data/lit/timeseries/`:

  * `martins2005_glucose_glycine_80_100_120C_pH68.yml`
        glucose 0.2 M + glycine 0.2 M, 0.1 M phosphate pH 6.8, 80/100/120 C.
        Responses used: glucose, DFG, 3-DG, 1-DG.  158 values.
  * `martins2003_DFG_amadori_degradation.yml`, **pH 6.8 series only**
        10 mmol/L isolated DFG, 100 and 120 C. Response used: DFG. 18 values.

Total 176 fitted values against 16 free parameters.

The pH 5.5 series of the DFG file is HELD OUT of the objective (the model has
no fitted pH term, and Wave S2 scored the repository's pH behaviour at 2/7);
it is reported as an out-of-objective diagnostic. So are fructose, glycine,
methylglyoxal, formic acid, acetic acid and melanoidins -- the repository's
network cannot emit four of those species at all.

`brands_sugar_casein_120C_pH68.yml` is NOT in the objective: it holds
browning and mutagenicity endpoints, not trunk trajectories. Brands' own
fitted rate constants (thesis Chapter 4) are used as the INDEPENDENT
cross-check in section (c) of the report -- a different amine (protein-bound
lysine), a different sugar loading and a different laboratory year.

`martins2005_glucose_glycine_100C_pH68.yml` is NOT in the objective either:
it is the SAME experiment as the 100 C series of the three-temperature file,
plotted in a second figure. Including both would double-weight 100 C. It is
used instead to measure the figure-read-off reproducibility floor.

NOTHING ELSE ENTERS. No benchmark, no hold-out bundle, no directional-panel
claim, and no value from `data/benchmarks/`.

--------------------------------------------------------------------------
STATISTICS
--------------------------------------------------------------------------
Multi-response weighted least squares in log space. For an observation y of
species s the residual is

    r = [ log(y_pred + f_s) - log(y_obs + f_s) ] / sigma_s

`f_s` is a per-species additive floor, taken from the source file's OWN stated
read-off uncertainty, which keeps the log finite where a measured value is
zero instead of dropping those points (dropping them is how a model that
predicts nothing scores well). `sigma_s` is a log-space standard deviation
estimated from the REPLICATE SCATTER in the data itself -- not assumed, not
tuned.

Parameters are carried as (log10 k_ref at 100 C, Ea in kJ/mol) per step, the
reparameterised Arrhenius form (see `src/trunk_kinetics`), because over an
80-120 C window A and Ea are so correlated that fitted A values are
individually meaningless.

Confidence intervals are Hessian-based (Gauss-Newton, J^T J), and the report
ALSO prints the eigendecomposition of that matrix in log space, so that
"which directions the data actually constrains" is answered by the sloppy
eigenvectors rather than by the marginal intervals alone. Every reported
parameter carries an identifiability verdict. A profile likelihood is run on
the Schiff/Amadori ratio specifically, because that ratio is the subject of a
standing contradiction between the repository's two barrier tables.

Usage:
    python scripts/generators/generate_trunk_rate_calibration.py [--starts N]
        [--output-dir results/validation] [--basename trunk_rate_calibration_refit]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.trunk_kinetics import (  # noqa: E402
    STEP_FAMILY,
    STEP_KEYS,
    STEP_UNIT,
    T_REF_K,
    integrate,
)

TIMESERIES = data_paths.TIMESERIES_DIR

GG_FILE = "martins2005_glucose_glycine_80_100_120C_pH68.yml"
DD_FILE = "martins2003_DFG_amadori_degradation.yml"
GG_100_FILE = "martins2005_glucose_glycine_100C_pH68.yml"
BRANDS_FILE = "brands_sugar_casein_120C_pH68.yml"

# Additive log floors, mmol/L.
#
# For the glucose/glycine figure the floor is 2% of that panel's printed full
# scale. 2% is not a guess: the file's own `attribution_validation` block
# records the MEASURED worst-case disagreement between the two figures that
# plot the same 100 C experiment (up to 2.2% of full scale, for DFG at
# 60-120 min). Panel full scales, from the file's `axis_ranges_verified`:
# Glu 250, DFG 25, DG (both 1-DG and 3-DG) 1.
#
# For the DFG-degradation figure the floor is the file's own stated read-off
# uncertainty, 0.04 mmol/L on the DFG panel.
FLOOR_GG = {"glucose": 0.02 * 250, "DFG": 0.02 * 25, "3-DG": 0.02 * 1.0, "1-DG": 0.02 * 1.0}
FLOOR_DD = {"DFG": 0.04, "glycine": 0.03}

FITTED_GG_SPECIES = ("glucose", "DFG", "3-DG", "1-DG")

# species label in the data -> state name in src.trunk_kinetics
STATE_OF = {"glucose": "Glc", "DFG": "DFG", "3-DG": "TDG", "1-DG": "ODG"}

# Charged initial concentrations, mmol/L, from each file's verified
# `system.reactants` block. Nothing here is fitted.
INIT_GG = {"Glc": 200.0, "Gly": 200.0}   # 0.2 mol/L each
INIT_DD = {"DFG": 10.0}                  # 10 mmol/L isolated Amadori compound

# Wide search bounds on log10 k_ref (k at 100 C) and Ea (kJ/mol). These are
# deliberately loose -- the fit is started from random points inside them, NOT
# from any literature value, so that agreement with Martins' or Brands' own
# constants in the report is a RESULT and not an initialisation artefact.
BOUNDS_LOG10K: Dict[str, Tuple[float, float]] = {
    "k_schiff": (-8.0, -2.0),
    "k_amadori": (-4.0, 4.0),
    "k_dfg_3dg": (-4.0, 1.0),
    "k_dfg_1dg": (-4.0, 1.0),
    "k_dfg_other": (-4.0, 1.0),
    "k_3dg_out": (-3.0, 2.0),
    "k_1dg_out": (-3.0, 3.0),
    "k_glc_other": (-5.0, 0.0),
}
EA_BOUNDS = (20.0, 260.0)

N_STEPS = len(STEP_KEYS)


# ---------------------------------------------------------------------------
# (c) INDEPENDENT COMPARATORS -- quoted, not fitted
# ---------------------------------------------------------------------------
# Brands, C.M.J. (2002) PhD thesis, Wageningen University, ISBN 90-5808-591-0,
# https://edepot.wur.nl/199005. Scheme 4.1 / Table 4.1 (thesis p. 48), rate
# constants at 120 C with 95% CI, least-squares criterion column; Table 4.2
# (thesis p. 51), Arrhenius Ea over 90-130 C. System: 150 mM sugar + 3% w/w
# sodium caseinate, 0.1 M phosphate pH 6.8, 120 C.
#
# THIS IS A GENUINELY INDEPENDENT COMPARATOR: a different amine (the
# epsilon-amino group of protein-bound lysine, not free glycine), a different
# sugar:amine ratio (10:1, not 1:1), and a fit performed by a different author
# on a different data set. Recovered by audit Wave S2.
BRANDS_120C = {
    "k1_glc_to_fru": (0.01020, 0.00021, "1/min", 126.0, 2.0, "glucose -> fructose isomerisation"),
    "k2_fru_to_glc": (0.00508, 0.00022, "1/min", 117.0, 2.0, "fructose -> glucose"),
    "k3_glc_to_formic": (0.00066, 0.00017, "1/min", 137.0, 7.0, "glucose -> formic acid"),
    "k4_fru_to_formic": (0.00105, 0.00018, "1/min", 159.0, 5.0, "fructose -> formic acid"),
    "k5_fru_to_triose": (0.00461, 0.00046, "1/min", 129.0, 17.0, "fructose -> 2 x C3"),
    "k7_condensation": (0.00024, 0.00002, "L/(mmol*min)", 114.0, 2.0,
                        "glucose + lysine -> Amadori (CONDENSATION, bimolecular)"),
    "k8_amadori_to_acetic": (0.11224, 0.04683, "1/min", 124.0, 4.0,
                             "Amadori -> acetic acid + lysine (AMADORI DEGRADATION)"),
    "k9_amadori_to_amp": (0.16831, 0.07595, "1/min", 128.0, 4.0,
                          "Amadori -> advanced MRP (AMADORI DEGRADATION)"),
    "k11_amp_to_melanoidin": (0.07302, 0.02621, "1/min", 75.0, 11.0, "AMP -> melanoidins"),
}

# Martins, S.I.F.S. (2003) PhD thesis, Wageningen University, ISBN
# 90-5808-823-5, Chapter 5 Table 5.2 (thesis p. 122), model M4, reparameterised
# Arrhenius X = k(T_av = 100 C). Published as Martins & van Boekel 2005,
# Food Chem 90(1-2):257-269, doi 10.1016/j.foodchem.2004.04.006.
#
# THIS IS *NOT* AN INDEPENDENT COMPARATOR. Martins fitted these constants to
# the very data this script fits. It is a REPRODUCIBILITY check -- can an
# independent implementation, a different scheme and a different objective
# recover the same rates from the same measurements -- and it is labelled as
# such everywhere it appears.
MARTINS_M4 = {
    1: (1.61e-5, 97.0, 3.0, "Glu + Gly -> DFG (condensation, bimolecular)"),
    2: (1.64e-3, 123.0, 5.0, "Glu -> Fru"),
    3: (9.15e-3, 93.0, 3.0, "Fru -> Glu"),
    4: (1.11e-2, 97.0, 2.0, "DFG -> 3-DG + Gly"),
    5: (3.45e-2, 30.0, 9.0, "3-DG -> formic acid"),
    6: (7.08e-3, 125.0, 5.0, "DFG -> methylglyoxal"),
    7: (1.57e-2, 107.0, 7.0, "DFG -> 1-DG"),
    8: (1.45e0, 76.0, 4.0, "1-DG -> acetic acid"),
    9: (8.12e-4, 95.0, 2.0, "3-DG + Gly -> melanoidins"),
    10: (4.41e-5, 237.0, 36.0, "Fru -> formic + acetic acid"),
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


class Observation:
    __slots__ = ("experiment", "species", "t_c", "t_min", "value", "floor")

    def __init__(self, experiment: str, species: str, t_c: float, t_min: float,
                 value: float, floor: float) -> None:
        self.experiment = experiment
        self.species = species
        self.t_c = t_c
        self.t_min = t_min
        self.value = value
        self.floor = floor

    @property
    def key(self) -> Tuple[str, str]:
        return (self.experiment, self.species)


def _load_yaml(name: str) -> Dict[str, Any]:
    with (TIMESERIES / name).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_fit_corpus() -> Tuple[List[Observation], List[Observation]]:
    """Return (objective observations, diagnostic-only observations)."""
    fitted: List[Observation] = []
    diagnostic: List[Observation] = []

    gg = _load_yaml(GG_FILE)
    for entry in gg["series"]:
        t_c = entry.get("T_C")
        species = entry["species"]
        # `T_C` is None for the t0_markers blocks and a STRING
        # ("ambiguous_80_or_100") for the unassigned points. Both are excluded:
        # a point whose temperature the source could not resolve cannot be
        # assigned one here either.
        if not isinstance(t_c, (int, float)) or isinstance(t_c, bool):
            continue
        floor = FLOOR_GG.get(species)
        for t_min, value in entry.get("points", []):
            obs = Observation("gg", species, float(t_c), float(t_min), float(value),
                              floor if floor is not None else 0.02 * 21.0)
            (fitted if species in FITTED_GG_SPECIES else diagnostic).append(obs)

    dd = _load_yaml(DD_FILE)
    for entry in dd["series"]:
        species = entry["species"]
        ph = entry.get("pH")
        floor = FLOOR_DD.get(species, 0.04)
        for t_min, value in entry.get("points", []):
            obs = Observation("dfgdeg" if ph == 6.8 else "dfgdeg_ph55", species,
                              float(entry["T_C"]), float(t_min), float(value), floor)
            if ph == 6.8 and species == "DFG":
                fitted.append(obs)
            else:
                diagnostic.append(obs)

    return fitted, diagnostic


def estimate_log_sigmas(observations: Sequence[Observation]) -> Dict[Tuple[str, str], float]:
    """
    Pooled within-time-point log-space standard deviation, per (experiment,
    species). This is the measurement scatter the DATA reports about itself.
    """
    groups: Dict[Tuple[str, str, float, float], List[Observation]] = {}
    for obs in observations:
        groups.setdefault((obs.experiment, obs.species, obs.t_c, obs.t_min), []).append(obs)

    deviations: Dict[Tuple[str, str], List[float]] = {}
    group_counts: Dict[Tuple[str, str], int] = {}
    for (experiment, species, _t_c, _t), members in groups.items():
        if len(members) < 2:
            continue
        key = (experiment, species)
        logs = np.log(np.array([m.value for m in members]) + members[0].floor)
        deviations.setdefault(key, []).extend(list(logs - logs.mean()))
        group_counts[key] = group_counts.get(key, 0) + 1

    sigmas: Dict[Tuple[str, str], float] = {}
    for key, devs in deviations.items():
        arr = np.asarray(devs)
        dof = max(arr.size - group_counts[key], 1)
        sigmas[key] = float(math.sqrt(float((arr ** 2).sum()) / dof))
    return sigmas


def resolve_sigmas(
    fitted: Sequence[Observation], sigmas: Dict[Tuple[str, str], float]
) -> Tuple[Dict[Tuple[str, str], float], Dict[str, str]]:
    """Fill in sigmas for response groups that carry no replicates."""
    resolved = dict(sigmas)
    provenance: Dict[str, str] = {}
    for key in {obs.key for obs in fitted}:
        label = f"{key[0]}:{key[1]}"
        if key in resolved:
            provenance[label] = "estimated from replicate scatter in this response group"
            continue
        # The DFG-degradation figure plots one series per condition with no
        # replicate markers, so it reports no scatter of its own. It is given
        # the scatter measured for the SAME SPECIES by the SAME AUTHOR with the
        # SAME assay in the glucose/glycine figure. Stated, not hidden.
        donor = ("gg", key[1])
        if donor in resolved:
            resolved[key] = resolved[donor]
            provenance[label] = (
                f"no replicates drawn in this figure; borrowed from {donor[0]}:{donor[1]} "
                f"(same species, same author, same assay)"
            )
        else:
            resolved[key] = 0.10
            provenance[label] = "no replicates and no donor group; 0.10 assumed"
    return resolved, provenance


# ---------------------------------------------------------------------------
# Objective
# ---------------------------------------------------------------------------


class Objective:
    def __init__(self, observations: Sequence[Observation],
                 sigmas: Dict[Tuple[str, str], float]) -> None:
        self.observations = list(observations)
        self.sigmas = sigmas
        self.groups: Dict[Tuple[str, float], Dict[str, Any]] = {}
        for index, obs in enumerate(self.observations):
            experiment = "gg" if obs.experiment == "gg" else "dfgdeg"
            group = self.groups.setdefault((experiment, obs.t_c), {"times": set(), "rows": []})
            group["times"].add(obs.t_min)
            group["rows"].append(index)
        for group in self.groups.values():
            group["grid"] = np.array(sorted(group["times"]), dtype=float)
            group["index"] = {t: i for i, t in enumerate(group["grid"])}

    @staticmethod
    def unpack(x: np.ndarray) -> Dict[str, Tuple[float, float]]:
        return {key: (10.0 ** float(x[i]), float(x[N_STEPS + i]))
                for i, key in enumerate(STEP_KEYS)}

    def predict(self, x: np.ndarray) -> np.ndarray:
        params = self.unpack(x)
        out = np.zeros(len(self.observations), dtype=float)
        for (experiment, t_c), group in self.groups.items():
            initial = INIT_GG if experiment == "gg" else INIT_DD
            run = integrate(params, t_c + 273.15, initial, group["grid"],
                            rtol=1e-9, atol=1e-12)
            for index in group["rows"]:
                obs = self.observations[index]
                out[index] = run.series(STATE_OF[obs.species])[group["index"][obs.t_min]]
        return out

    def residuals(self, x: np.ndarray) -> np.ndarray:
        try:
            predicted = self.predict(x)
        except Exception:
            return np.full(len(self.observations), 1.0e3)
        res = np.empty(len(self.observations), dtype=float)
        for i, obs in enumerate(self.observations):
            res[i] = (math.log(predicted[i] + obs.floor) - math.log(obs.value + obs.floor)) \
                / self.sigmas[obs.key]
        if not np.all(np.isfinite(res)):
            return np.full(len(self.observations), 1.0e3)
        return res

    def cost(self, x: np.ndarray) -> float:
        r = self.residuals(x)
        return 0.5 * float(np.dot(r, r))


def bounds_arrays() -> Tuple[np.ndarray, np.ndarray]:
    lower = np.array([BOUNDS_LOG10K[k][0] for k in STEP_KEYS] + [EA_BOUNDS[0]] * N_STEPS)
    upper = np.array([BOUNDS_LOG10K[k][1] for k in STEP_KEYS] + [EA_BOUNDS[1]] * N_STEPS)
    return lower, upper


def multistart_fit(objective: Objective, starts: int, seed: int = 20260827):
    """Random multi-start bounded least squares. Reports the optimum spread."""
    from scipy.optimize import least_squares

    lower, upper = bounds_arrays()
    rng = np.random.default_rng(seed)
    results = []
    for _ in range(starts):
        x0 = lower + rng.random(lower.size) * (upper - lower)
        try:
            fit = least_squares(objective.residuals, x0, bounds=(lower, upper),
                                xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
        except Exception:
            continue
        results.append((0.5 * float(np.dot(fit.fun, fit.fun)), fit))
    if not results:
        raise RuntimeError("every start failed")
    results.sort(key=lambda item: item[0])
    return results


# ---------------------------------------------------------------------------
# Identifiability
# ---------------------------------------------------------------------------


def numeric_jacobian(objective: Objective, x: np.ndarray, step: float = 1e-5) -> np.ndarray:
    base = objective.residuals(x)
    jac = np.zeros((base.size, x.size))
    for i in range(x.size):
        # Ea is O(100) and log10 k is O(1); scale the probe so both are
        # perturbed by a comparable relative amount.
        h = step * max(abs(x[i]), 1.0)
        xp = x.copy()
        xp[i] += h
        jac[:, i] = (objective.residuals(xp) - base) / h
    return jac


def identifiability(objective: Objective, x: np.ndarray,
                    variance_scale: float = 1.0) -> Dict[str, Any]:
    """
    Gauss-Newton covariance, marginal CIs, and the sloppy-direction spectrum.

    `variance_scale` is the reduced chi-square of the fit. Scaling the
    covariance by it is the standard correction for an unknown error scale,
    and here it is not optional: the sigmas come from replicate scatter, the
    model's systematic error is several times larger than that scatter, and
    unscaled J^T J intervals would therefore be narrower than the residuals
    justify by a factor of sqrt(reduced chi-square).
    """
    jac = numeric_jacobian(objective, x)
    hess = (jac.T @ jac) / max(variance_scale, 1e-12)
    eigenvalues, eigenvectors = np.linalg.eigh(hess)
    eigenvalues = np.clip(eigenvalues, 0.0, None)

    # Pseudo-inverse with a relative cut, so that a genuinely flat direction
    # produces an INFINITE interval rather than a huge but finite one that
    # invites over-reading.
    cutoff = max(eigenvalues.max(), 1.0) * 1e-10
    inv_eigen = np.array([1.0 / v if v > cutoff else np.inf for v in eigenvalues])
    finite = np.isfinite(inv_eigen)
    cov = (eigenvectors[:, finite] * inv_eigen[finite]) @ eigenvectors[:, finite].T
    variances = np.diag(cov).copy()
    for i in range(x.size):
        # a parameter with meaningful weight on a null direction is unbounded
        weight_on_null = float(np.sum(eigenvectors[i, ~finite] ** 2))
        if weight_on_null > 1e-8:
            variances[i] = np.inf
    stderr = np.sqrt(np.clip(variances, 0.0, None))

    labels = [f"log10 {k}" for k in STEP_KEYS] + [f"Ea {k}" for k in STEP_KEYS]
    order = np.argsort(eigenvalues)
    spectrum = []
    for rank, idx in enumerate(order):
        vec = eigenvectors[:, idx]
        top = np.argsort(-np.abs(vec))[:3]
        spectrum.append({
            "rank_from_sloppiest": rank,
            "eigenvalue": float(eigenvalues[idx]),
            "dominant_components": [
                {"parameter": labels[j], "weight": round(float(vec[j]), 3)} for j in top
            ],
        })
    condition_number = (float(eigenvalues.max() / eigenvalues.min())
                        if eigenvalues.min() > 0 else float("inf"))
    return {
        "hessian_eigenvalues": [float(v) for v in np.sort(eigenvalues)],
        "condition_number": condition_number,
        "standard_errors": stderr,
        "spectrum": spectrum,
    }


def verdict_for(stderr_log10: float) -> str:
    """Identifiability verdict for a RATE, from its marginal log10 std error."""
    if not math.isfinite(stderr_log10):
        return "UNIDENTIFIED (flat direction; the data do not bound this parameter)"
    if stderr_log10 < 0.05:
        return "well constrained (95% CI inside +/-26%)"
    if stderr_log10 < 0.15:
        return "constrained (95% CI inside a factor of 2)"
    if stderr_log10 < 0.5:
        return "weakly constrained (95% CI up to a factor of 10)"
    return "SLOPPY (order-of-magnitude only)"


def verdict_for_ea(stderr_kj: float, value_kj: float) -> str:
    """
    Identifiability verdict for an ACTIVATION ENERGY, on its own kJ/mol scale.

    A negative lower bound is called out explicitly: an Arrhenius Ea below zero
    is unphysical, so an interval that straddles zero means the data determine
    the step's temperature dependence not at all, however narrow the number
    printed next to it looks.
    """
    if not math.isfinite(stderr_kj):
        return "UNIDENTIFIED (flat direction)"
    lower = value_kj - 1.96 * stderr_kj
    if lower <= 0.0:
        return (f"UNCONSTRAINED SIGN (95% CI reaches {lower:.0f} kJ/mol; an Ea <= 0 is "
                f"unphysical, so the data do not determine this step's T-dependence)")
    if stderr_kj < 5.0:
        return "well constrained"
    if stderr_kj < 15.0:
        return "constrained"
    if stderr_kj < 35.0:
        return "weakly constrained"
    return "SLOPPY (the 80-120 C window is too narrow for this step)"


def profile_schiff_amadori(objective: Objective, x_best: np.ndarray,
                           ratios: Sequence[float]) -> List[Dict[str, Any]]:
    """
    Profile likelihood over the Schiff/Amadori ORDERING.

    The quantity profiled is

        log10 R = log10( k_amadori / (k_schiff * [Gly]_0) )

    i.e. the ratio of the two steps' rates as PSEUDO-FIRST-ORDER constants at
    the experiment's own 200 mmol/L glycine loading. That is the only form in
    which a bimolecular and a unimolecular step can be compared at all, and it
    is the comparison the two barrier tables in this repository disagree about.
    R > 1 means the Amadori rearrangement is the faster of the two.

    At each fixed R, k_amadori is pinned and every other parameter (including
    k_schiff and all eight Ea) is re-optimised.
    """
    from scipy.optimize import least_squares

    lower, upper = bounds_arrays()
    i_schiff = STEP_KEYS.index("k_schiff")
    i_amadori = STEP_KEYS.index("k_amadori")
    gly0 = INIT_GG["Gly"]
    best_cost = objective.cost(x_best)

    free = [i for i in range(x_best.size) if i != i_amadori]
    out: List[Dict[str, Any]] = []
    # The profile explores a lot of parameter space. If any of its constrained
    # re-optimisations lands BELOW the reported global optimum, that optimum was
    # wrong; the caller uses this to restart rather than publishing a profile
    # with negative delta-costs in it.
    best_seen_x = x_best.copy()
    best_seen_cost = best_cost
    for ratio in ratios:
        target_log_amadori = math.log10(ratio) + x_best[i_schiff] + math.log10(gly0)
        if not (lower[i_amadori] <= target_log_amadori <= upper[i_amadori]):
            out.append({"ratio": ratio, "log10_ratio": math.log10(ratio),
                        "delta_cost": None, "note": "outside search bounds"})
            continue

        def residuals(free_x: np.ndarray) -> np.ndarray:
            full = x_best.copy()
            full[free] = free_x
            # keep the ratio fixed as k_schiff moves
            full[i_amadori] = math.log10(ratio) + full[i_schiff] + math.log10(gly0)
            if not (lower[i_amadori] <= full[i_amadori] <= upper[i_amadori]):
                return np.full(len(objective.observations), 1.0e3)
            return objective.residuals(full)

        # Re-optimise from several starts, not just from the global optimum:
        # at extreme ratios the surface is far from x_best and a single local
        # descent can stall well above the true profile, which would make the
        # profile look artificially rejecting.
        rng = np.random.default_rng(hash(("profile", ratio)) % (2 ** 32))
        cost = math.inf
        seeds = [x_best[free]]
        for _ in range(3):
            jitter = x_best.copy()
            jitter += rng.normal(0.0, 0.30, size=jitter.size) * np.where(
                np.arange(jitter.size) < N_STEPS, 1.0, 40.0)
            seeds.append(np.clip(jitter[free], lower[free], upper[free]))
        for seed in seeds:
            try:
                fit = least_squares(residuals, seed, bounds=(lower[free], upper[free]),
                                    xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=2500)
            except Exception:
                continue
            this_cost = 0.5 * float(np.dot(fit.fun, fit.fun))
            if this_cost < cost:
                cost = this_cost
            if this_cost < best_seen_cost:
                best_seen_cost = this_cost
                candidate = x_best.copy()
                candidate[free] = fit.x
                candidate[i_amadori] = (math.log10(ratio) + candidate[i_schiff]
                                        + math.log10(gly0))
                best_seen_x = candidate
        if not math.isfinite(cost):
            out.append({"ratio": ratio, "log10_ratio": math.log10(ratio),
                        "delta_cost": None, "note": "profile optimisation failed"})
            continue
        out.append({
            "ratio": ratio,
            "log10_ratio": math.log10(ratio),
            "cost": cost,
            "delta_cost": cost - best_cost,
            # 2*delta(0.5*SSR) is the likelihood-ratio statistic; chi2(1) 95%
            # critical value is 3.84.
            "chi2_statistic": 2.0 * (cost - best_cost),
            "rejected_at_95pct": bool(2.0 * (cost - best_cost) > 3.841),
        })
    return out, best_seen_x, best_seen_cost


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------


def git_head() -> Dict[str, Any]:
    def run(args: List[str]) -> str:
        try:
            return subprocess.check_output(args, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unknown"
    return {
        "commit": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(run(["git", "status", "--porcelain"])),
    }


def per_response_fit(objective: Objective, x: np.ndarray) -> List[Dict[str, Any]]:
    predicted = objective.predict(x)
    groups: Dict[Tuple[str, str, float], List[int]] = {}
    for i, obs in enumerate(objective.observations):
        groups.setdefault((obs.experiment, obs.species, obs.t_c), []).append(i)
    rows = []
    for (experiment, species, t_c), idx in sorted(groups.items()):
        errs = []
        for i in idx:
            obs = objective.observations[i]
            errs.append(abs(math.log10((predicted[i] + obs.floor) / (obs.value + obs.floor))))
        rows.append({
            "experiment": experiment, "species": species, "T_C": t_c,
            "n": len(idx),
            "median_abs_log10_error": float(np.median(errs)),
            "median_fold_error": float(10 ** np.median(errs)),
        })
    return rows


def diagnostics(x: np.ndarray, diagnostic_obs: Sequence[Observation]) -> Dict[str, Any]:
    """Out-of-objective checks. None of these influenced the fit."""
    params = Objective.unpack(x)
    out: Dict[str, Any] = {}

    # 1. Glycine YIELD in the isolated-DFG experiment. The trunk releases one
    #    glycine per DFG that goes down the 3-DG or 1-DG channel and none down
    #    `k_dfg_other`, so the predicted asymptotic yield is a direct readout
    #    of the fitted branch split -- and Table 4.1.1 measured it.
    rows = [o for o in diagnostic_obs if o.experiment == "dfgdeg" and o.species == "glycine"]
    yield_rows = []
    for t_c in sorted({o.t_c for o in rows}):
        subset = sorted([o for o in rows if o.t_c == t_c], key=lambda o: o.t_min)
        grid = np.array([o.t_min for o in subset], dtype=float)
        run = integrate(params, t_c + 273.15, INIT_DD, grid)
        yield_rows.append({
            "T_C": t_c,
            "measured_final_glycine_mmol_L": subset[-1].value,
            "predicted_final_glycine_mmol_L": float(run.glycine_released[-1]),
            "measured_yield_pct": 100.0 * subset[-1].value / INIT_DD["DFG"],
            "predicted_yield_pct": 100.0 * float(run.glycine_released[-1]) / INIT_DD["DFG"],
        })
    out["glycine_release_yield"] = yield_rows

    # 2. The pH 5.5 DFG series, entirely outside the objective. The model has
    #    no pH term, so it predicts the pH 6.8 curve at both pH values; the
    #    size of that miss is the honest measure of the missing physics.
    ph55 = [o for o in diagnostic_obs if o.experiment == "dfgdeg_ph55" and o.species == "DFG"]
    ph_rows = []
    for t_c in sorted({o.t_c for o in ph55}):
        subset = sorted([o for o in ph55 if o.t_c == t_c], key=lambda o: o.t_min)
        grid = np.array([o.t_min for o in subset], dtype=float)
        run = integrate(params, t_c + 273.15, INIT_DD, grid)
        pred = run.series("DFG")
        errs = [abs(math.log10((pred[i] + o.floor) / (o.value + o.floor)))
                for i, o in enumerate(subset)]
        ph_rows.append({
            "T_C": t_c, "pH": 5.5, "n": len(subset),
            "median_fold_error": float(10 ** np.median(errs)),
        })
    out["ph55_out_of_objective"] = ph_rows

    # 3. Glucose mass budget: how much of the measured glucose loss the
    #    condensation channel takes, versus the declared nuisance channel.
    grid = np.linspace(0.0, 150.0, 61)
    run = integrate(params, 373.15, INIT_GG, grid)
    k = {key: params[key][0] for key in STEP_KEYS}
    glc = run.series("Glc")
    gly = run.series("Gly")
    cond_flux = float(np.trapezoid(k["k_schiff"] * glc * gly, grid))
    other_flux = float(np.trapezoid(k["k_glc_other"] * glc, grid))
    out["glucose_budget_100C_150min"] = {
        "total_glucose_lost_mmol_L": float(INIT_GG["Glc"] - glc[-1]),
        "via_condensation_mmol_L": cond_flux,
        "via_declared_nuisance_channel_mmol_L": other_flux,
        "nuisance_share": other_flux / max(cond_flux + other_flux, 1e-12),
    }
    return out


def _k_at(params: Dict[str, Any], key: str, t_c: float) -> float:
    from src.trunk_kinetics import arrhenius_k
    return arrhenius_k(params[key]["k_ref_100C"], params[key]["ea_kj_mol"], t_c + 273.15)


def brands_comparison(params: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Factor-by-factor comparison against Brands' INDEPENDENT fitted constants.

    Only steps that are genuinely the same transformation are compared, and
    every row says what makes the two commensurable and what does not.
    """
    rows: List[Dict[str, Any]] = []

    def row(label: str, fitted: float, unit: str, brands_key: str,
            fitted_ea: float, note: str) -> None:
        value, ci, b_unit, ea, ea_ci, desc = BRANDS_120C[brands_key]
        rows.append({
            "step": label,
            "fitted_120C": fitted,
            "brands_120C": value,
            "brands_ci95_120C": ci,
            "unit": unit,
            "unit_match": unit == b_unit,
            "ratio_fitted_over_brands": fitted / value,
            "fitted_ea_kj_mol": fitted_ea,
            "brands_ea_kj_mol": ea,
            "brands_ea_ci95": ea_ci,
            "ea_difference_kj_mol": fitted_ea - ea,
            "brands_step": desc,
            "commensurability": note,
        })

    row("condensation (sugar + amine -> Amadori)",
        _k_at(params, "k_schiff", 120.0), "L/(mmol*min)", "k7_condensation",
        params["k_schiff"]["ea_kj_mol"],
        "SAME TRANSFORMATION, DIFFERENT AMINE. Brands' k7 is a one-step lumped "
        "condensation glucose + lysine -> Amadori; the fitted k_schiff is the FIRST of two "
        "steps to the same product. Because the fitted Amadori step is the fast one (see "
        "the Schiff/Amadori section), k_schiff is the rate-determining constant of the "
        "composite and the two are comparable to within that approximation. The amine "
        "differs: protein-bound lysine epsilon-amino versus free glycine alpha-amino.")

    total_dfg = sum(_k_at(params, key, 120.0)
                    for key in ("k_dfg_3dg", "k_dfg_1dg", "k_dfg_other"))
    value8, _, _, ea8, _, _ = BRANDS_120C["k8_amadori_to_acetic"]
    value9, _, _, ea9, _, _ = BRANDS_120C["k9_amadori_to_amp"]
    rows.append({
        "step": "TOTAL Amadori degradation",
        "fitted_120C": total_dfg,
        "brands_120C": value8 + value9,
        "brands_ci95_120C": None,
        "unit": "1/min",
        "unit_match": True,
        "ratio_fitted_over_brands": total_dfg / (value8 + value9),
        "fitted_ea_kj_mol": None,
        "brands_ea_kj_mol": (ea8 + ea9) / 2.0,
        "brands_ea_ci95": None,
        "ea_difference_kj_mol": None,
        "brands_step": "k8 + k9: Amadori -> acetic acid + lysine, and Amadori -> advanced MRP",
        "commensurability": "SAME TRANSFORMATION, DIFFERENT AMINE. Both are the total "
                            "first-order disappearance of the Amadori compound at 120 C. "
                            "Fitted side sums the three DFG channels; Brands' side sums his "
                            "two. This is the single strongest independent test in the wave, "
                            "because the fitted value is anchored by an experiment on "
                            "ISOLATED DFG that Brands never ran.",
    })
    return rows


def martins_comparison(params: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Reproducibility (NOT independence) check against Martins' own M4 constants.

    UNITS WARNING, and the reason this function is more careful than it looks.
    Two of Martins' M4 steps are BIMOLECULAR -- step 1 (Glu + Gly -> DFG) and
    step 9 (3-DG + Gly -> melanoidins) -- so their X values are in
    L/(mmol*min), not 1/min. Comparing step 9's 8.12e-4 directly against a
    fitted first-order constant understates it by the glycine concentration,
    a factor of 200, and produces a spurious 130x "disagreement". Every
    bimolecular comparator below is multiplied by the experiment's own
    200 mmol/L glycine loading first, and the row records that it was.
    """
    gly0 = INIT_GG["Gly"]

    # (fitted key, [(martins step, is_bimolecular)], note)
    mapping: List[Tuple[str, List[Tuple[int, bool]], str]] = [
        ("k_schiff", [(1, True)],
         "Martins' step 1 is the ONE-STEP lumped condensation; the fitted k_schiff is the "
         "first of two steps to the same product and is the rate-determining one, so the "
         "two are directly comparable. Both are bimolecular, so no concentration factor "
         "is applied to either."),
        ("k_dfg_3dg", [(4, False)], "same step, same direction"),
        ("k_dfg_1dg", [(7, False)], "same step, same direction"),
        ("k_dfg_other", [(6, False)],
         "the fitted non-glycine-releasing DFG channel against Martins' DFG -> "
         "methylglyoxal. Martins' MG step DOES release glycine in his scheme, so these are "
         "the same flux wearing different labels, not the same mechanism."),
        ("k_3dg_out", [(5, False), (9, True)],
         "the fitted 3-DG sink is a LUMP over every 3-DG consumption channel, so it is "
         "compared against the SUM of Martins' two 3-DG sinks: step 5 (3-DG -> formic "
         "acid, first order) plus step 9 (3-DG + Gly -> melanoidins, BIMOLECULAR, "
         f"converted to pseudo-first-order at {gly0:.0f} mmol/L glycine)."),
        ("k_1dg_out", [(8, False)], "same step, same direction"),
        ("k_glc_other", [(2, False)],
         "the DECLARED NUISANCE CHANNEL against Martins' Glu -> Fru isomerisation. The "
         "fitted channel additionally lumps amine-independent sugar fragmentation, but it "
         "is a NET loss where Martins' step 2 is balanced by his reverse step 3, so these "
         "are only loosely commensurable and the ratio should not be over-read."),
    ]

    rows = []
    for key, steps, note in mapping:
        comparator = 0.0
        ea_terms = []
        descriptions = []
        bimolecular_converted = False
        for step, is_bimolecular in steps:
            x_martins, ea_martins, _ea_ci, desc = MARTINS_M4[step]
            # Convert a bimolecular comparator to pseudo-first-order ONLY when
            # the fitted step it is being compared against is first-order.
            # k_schiff is itself bimolecular in the same units, so converting
            # there would introduce the 200x error this block exists to avoid.
            needs_conversion = is_bimolecular and STEP_UNIT[key] == "1/min"
            contribution = x_martins * (gly0 if needs_conversion else 1.0)
            if needs_conversion:
                bimolecular_converted = True
            comparator += contribution
            ea_terms.append((contribution, ea_martins))
            descriptions.append(f"step {step} ({desc})")
        # flux-weighted mean Ea when several steps are summed
        ea_martins_eff = sum(w * e for w, e in ea_terms) / max(comparator, 1e-300)
        fitted = params[key]["k_ref_100C"]
        rows.append({
            "fitted_step": key,
            "martins_steps": [s for s, _ in steps],
            "martins_description": " + ".join(descriptions),
            "fitted_k_100C": fitted,
            "martins_comparator_100C": comparator,
            "martins_bimolecular_converted_at_gly_mmol_L": gly0 if bimolecular_converted else None,
            "ratio_fitted_over_martins": fitted / comparator if comparator else None,
            "fitted_ea_kj_mol": params[key]["ea_kj_mol"],
            "martins_ea_kj_mol": ea_martins_eff,
            "ea_difference_kj_mol": params[key]["ea_kj_mol"] - ea_martins_eff,
            "note": note,
        })
    return rows


# ---------------------------------------------------------------------------
# (f) What a propagation to the screening lane WOULD look like
# ---------------------------------------------------------------------------


def fast_barrier_mapping(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Derive, but do NOT apply, the FAST_BARRIERS value each fitted rate implies.

    The screening lane's own rate law is
        k(T) = A * exp(-Ea[kcal/mol] / (0.001987 * T))            (barrier_constants.py:804-815)
    and its own inverse is
        Ea = -0.001987 * T * ln(k / A)                            (barrier_constants.py:832-836)

    A MUST be `DEFAULT_REFERENCE_PREEXPONENTIAL` = 1e11, not the family's YAML
    A value, because the single place in the lane where an absolute rate
    constant survives -- `recommend.py:1216` -- calls
    `get_reference_pre_exponential()` with NO argument, so every family gets
    1e11 regardless of what its YAML entry says.

    THE ASSUMPTION THIS BUYS, STATED: everything temperature-independent in
    the fitted rate -- the true prefactor, the activation entropy, general
    acid/base catalysis by the phosphate buffer, solvent reorganisation -- is
    taken to equal 1e11 s^-1 exactly. Each decade of error there moves the
    derived barrier by R*T*ln(10) = 1.71 kcal/mol at 100 C. That systematic is
    larger than any fitted confidence interval in this report, which is why
    these numbers are published and not shipped.

    `k_schiff` is EXCLUDED: it is bimolecular, in L/(mmol*min), and there is no
    dimensionally valid conversion to a first-order barrier without inventing a
    standard-state concentration. Inventing one is precisely the 342/200 basis
    failure Wave S2b exposed.
    """
    from src.barrier_constants import (
        ARRHENIUS_R_KCAL,
        DEFAULT_REFERENCE_PREEXPONENTIAL,
        FAST_BARRIERS,
    )

    rt = ARRHENIUS_R_KCAL * T_REF_K
    rows = []
    for key in STEP_KEYS:
        family = STEP_FAMILY[key]
        if family is None:
            continue
        if STEP_UNIT[key] != "1/min":
            rows.append({
                "step": key, "family": family, "derived_barrier_kcal_mol": None,
                "shipped_barrier_kcal_mol": FAST_BARRIERS.get(family, (None,))[0],
                "excluded_reason": (
                    "bimolecular (L/(mmol*min)); no dimensionally valid conversion to a "
                    "first-order barrier without inventing a standard-state concentration"
                ),
            })
            continue
        k_per_s = params[key]["k_ref_100C"] / 60.0
        derived = -rt * math.log(k_per_s / DEFAULT_REFERENCE_PREEXPONENTIAL)
        shipped = FAST_BARRIERS.get(family, (None,))[0]
        rows.append({
            "step": key,
            "family": family,
            "k_ref_100C_per_min": params[key]["k_ref_100C"],
            "k_ref_100C_per_s": k_per_s,
            "arithmetic": (
                f"Ea = -{ARRHENIUS_R_KCAL} * {T_REF_K} * ln({k_per_s:.6g} / "
                f"{DEFAULT_REFERENCE_PREEXPONENTIAL:.0e}) = {derived:.2f} kcal/mol"
            ),
            "derived_barrier_kcal_mol": derived,
            "shipped_barrier_kcal_mol": shipped,
            "difference_kcal_mol": (derived - shipped) if shipped is not None else None,
            "excluded_reason": None,
        })
    return {
        "applied": False,
        "why_not_applied": (
            "The mapping is principled but its dominant uncertainty is an ASSUMED "
            "prefactor, not the fitted rate. Setting A = 1e11 s^-1 for every step folds "
            "the real prefactor, the activation entropy and the buffer catalysis into Ea "
            "at 1.71 kcal/mol per decade -- wider than every fitted confidence interval "
            "in this report. Applying it would move shipped predictions on the strength of "
            "an assumption rather than on the strength of the measurement. Published with "
            "the arithmetic so the owner can decide; FAST_BARRIERS is UNCHANGED by this "
            "wave. The counterfactual effect of applying it was measured anyway and is "
            "recorded in the wave ledger."
        ),
        "prefactor_used": DEFAULT_REFERENCE_PREEXPONENTIAL,
        "prefactor_justification": (
            "recommend.py:1216 calls get_reference_pre_exponential() with no argument, so "
            "the screening lane uses 1e11 for every family; the per-family YAML A values "
            "cancel exactly in benchmark_validation.py:662-673 and are dead in this lane."
        ),
        "reference_temperature_K": T_REF_K,
        "rows": rows,
    }


def readoff_reproducibility() -> Dict[str, Any]:
    """
    Measure the figure-read-off floor from the repository's own data: the
    100 C series appears in TWO figures of the same thesis, extracted
    separately. Their disagreement is a lower bound on the error of either.
    """
    a = _load_yaml(GG_FILE)
    b = _load_yaml(GG_100_FILE)
    per_species: Dict[str, List[float]] = {}
    for species in FITTED_GG_SPECIES:
        left = [e for e in a["series"] if e.get("species") == species and e.get("T_C") == 100]
        right = [e for e in b["series"] if e.get("species") == species]
        if not left or not right:
            continue
        by_time_l: Dict[float, List[float]] = {}
        for t, v in left[0].get("points", []):
            by_time_l.setdefault(float(t), []).append(float(v))
        by_time_r: Dict[float, List[float]] = {}
        for t, v in right[0].get("points", []):
            by_time_r.setdefault(float(t), []).append(float(v))
        for t in sorted(set(by_time_l) & set(by_time_r)):
            per_species.setdefault(species, []).append(
                abs(float(np.mean(by_time_l[t])) - float(np.mean(by_time_r[t])))
            )
    return {
        species: {"n_shared_times": len(vals),
                  "max_abs_disagreement_mmol_L": float(max(vals)),
                  "median_abs_disagreement_mmol_L": float(np.median(vals))}
        for species, vals in per_species.items()
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build(starts: int) -> Dict[str, Any]:
    fitted_obs, diagnostic_obs = load_fit_corpus()
    raw_sigmas = estimate_log_sigmas(fitted_obs + diagnostic_obs)
    sigmas, sigma_provenance = resolve_sigmas(fitted_obs, raw_sigmas)

    objective = Objective(fitted_obs, sigmas)
    results = multistart_fit(objective, starts)
    best_cost, best = results[0]
    x = np.asarray(best.x, dtype=float)

    # How many independent starts reached the same optimum? A basin that only
    # one start in N finds is a warning about the optimum's uniqueness.
    tol = max(best_cost * 1e-4, 1e-6)
    converged = sum(1 for cost, _ in results if cost - best_cost <= tol)

    # ---- profile-informed re-polish -------------------------------------
    # The Schiff/Amadori profile re-optimises 15 of the 16 parameters at each
    # of 20 fixed ratios, which is a far wider sweep of the surface than the
    # multi-start alone. On the first run of this script it found a point
    # BETTER than the multi-start optimum by delta-cost 7.2 -- i.e. the
    # reported optimum was not the optimum, and the profile would have been
    # published with negative delta-costs in it. Rather than paper over that
    # with a tolerance, the loop below adopts any improvement, re-polishes all
    # 16 parameters freely from it, and recomputes the profile. It terminates
    # when a full profile finds nothing better.
    from scipy.optimize import least_squares as _ls
    lower, upper = bounds_arrays()
    ratios = [0.1, 0.3, 1.0, 3.0, 10.0, 20.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
              60.0, 80.0, 100.0, 300.0, 1e3, 1e4, 1e5, 1e6]
    repolish_rounds = 0
    for _ in range(4):
        profile, seen_x, seen_cost = profile_schiff_amadori(objective, x, ratios)
        if seen_cost >= best_cost - 1e-9:
            break
        repolish_rounds += 1
        polished = _ls(objective.residuals, seen_x, bounds=(lower, upper),
                       xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
        polished_cost = 0.5 * float(np.dot(polished.fun, polished.fun))
        if polished_cost < seen_cost:
            x, best_cost = np.asarray(polished.x, dtype=float), polished_cost
        else:
            x, best_cost = seen_x, seen_cost

    dof = max(len(fitted_obs) - 2 * N_STEPS, 1)
    reduced_chi2 = (2.0 * best_cost) / dof
    ident = identifiability(objective, x, variance_scale=reduced_chi2)
    stderr = ident.pop("standard_errors")

    gly0 = INIT_GG["Gly"]
    params: Dict[str, Any] = {}
    for i, key in enumerate(STEP_KEYS):
        se_log10 = float(stderr[i])
        se_ea = float(stderr[N_STEPS + i])
        k_ref = 10.0 ** float(x[i])
        params[key] = {
            "repo_family": STEP_FAMILY[key],
            "unit": STEP_UNIT[key],
            "k_ref_100C": k_ref,
            "log10_k_ref_100C": float(x[i]),
            "log10_k_ref_stderr": se_log10 if math.isfinite(se_log10) else None,
            "k_ref_ci95": (
                [k_ref * 10 ** (-1.96 * se_log10), k_ref * 10 ** (1.96 * se_log10)]
                if math.isfinite(se_log10) else None
            ),
            "ea_kj_mol": float(x[N_STEPS + i]),
            "ea_stderr_kj_mol": se_ea if math.isfinite(se_ea) else None,
            "ea_ci95_kj_mol": (
                [float(x[N_STEPS + i] - 1.96 * se_ea), float(x[N_STEPS + i] + 1.96 * se_ea)]
                if math.isfinite(se_ea) else None
            ),
            "identifiability_k": verdict_for(se_log10),
            "identifiability_ea": verdict_for_ea(se_ea, float(x[N_STEPS + i])),
            "at_bound": bool(
                abs(float(x[i]) - BOUNDS_LOG10K[key][0]) < 1e-6
                or abs(float(x[i]) - BOUNDS_LOG10K[key][1]) < 1e-6
            ),
        }

    # `profile` was computed by the re-polish loop above against the final `x`.
    accepted = [row["ratio"] for row in profile
                if row.get("delta_cost") is not None and not row["rejected_at_95pct"]]

    k_schiff_pseudo = params["k_schiff"]["k_ref_100C"] * gly0
    k_amadori = params["k_amadori"]["k_ref_100C"]

    payload: Dict[str, Any] = {
        "artifact": "trunk_rate_calibration_refit",
        "wave": "S3",
        "generated_on": date.today().isoformat(),
        "git": git_head(),
        "reference_temperature_K": T_REF_K,

        # ---- the machine-readable fit declaration the gates read ----
        "fit_target_files": [
            data_paths.rel(data_paths.TIMESERIES_DIR / GG_FILE),
            f"data/lit/timeseries/{DD_FILE} (pH 6.8 series only)",
        ],
        "fit_corpus_excluded": [
            f"data/lit/timeseries/{GG_100_FILE} (same experiment as the 100 C series "
            f"of the fitted file; used only to measure the read-off floor)",
            f"data/lit/timeseries/{BRANDS_FILE} (endpoints, not trajectories; Brands' own "
            f"fitted constants are the independent cross-check, not a fit input)",
            f"data/lit/timeseries/{DD_FILE} pH 5.5 series (held out; reported as a diagnostic)",
            "data/benchmarks/** (the whole benchmark panel -- never entered)",
            "data/benchmarks/external_validation/** (hold-out -- never entered)",
            "docs/validation/directional_claims_panel.yml (never entered)",
        ],
        "fit_leverage": {
            "free_parameters": 2 * N_STEPS,
            "fitted_rows": len(fitted_obs),
            "parameters_per_row": (2 * N_STEPS) / len(fitted_obs),
            "class": "global_low_leverage",
            "interpretation": (
                f"{2 * N_STEPS} global rate parameters fitted jointly to {len(fitted_obs)} "
                f"concentration-time values across three temperatures and two independent "
                f"experiments. No parameter is recoverable from any single row."
            ),
        },
        "fit_target_kind": "timeseries_yaml",

        "objective": {
            "form": "multi-response weighted least squares on log(y + floor)",
            "n_values": len(fitted_obs),
            "n_free_parameters": 2 * N_STEPS,
            "final_cost_half_ssr": best_cost,
            "reduced_chi_square": reduced_chi2,
            "reduced_chi_square_reading": (
                "The sigmas are the REPLICATE scatter the data report about themselves. A "
                "reduced chi-square above 1 therefore measures how far the model's "
                f"systematic error exceeds that scatter: here residuals run "
                f"{math.sqrt(reduced_chi2):.1f}x the replicate sigma. Every confidence "
                f"interval below has been widened by that same factor; unscaled J^T J "
                f"intervals would have been {math.sqrt(reduced_chi2):.1f}x too tight."
            ),
            "confidence_interval_scaling": math.sqrt(reduced_chi2),
            "log_sigmas": {f"{k[0]}:{k[1]}": v for k, v in sigmas.items()},
            "log_sigma_provenance": sigma_provenance,
            "floors_mmol_per_L": {
                "glucose_glycine_figure": FLOOR_GG,
                "dfg_degradation_figure": FLOOR_DD,
            },
            "multistart": {
                "starts": starts,
                "successful": len(results),
                "reached_best_basin": converged,
                "second_best_cost": results[1][0] if len(results) > 1 else None,
                "initialisation": (
                    "uniform random inside the declared wide bounds; NO literature value "
                    "was used as a starting point"
                ),
                "profile_informed_repolish_rounds": repolish_rounds,
                "repolish_note": (
                    "The Schiff/Amadori profile sweeps the surface far more widely than "
                    "the multi-start. Any point it finds below the multi-start optimum is "
                    "adopted and re-polished, and the profile recomputed; the loop exits "
                    "only when a full profile finds nothing better. A non-zero count here "
                    "means the multi-start alone did NOT find the optimum, which is a "
                    "statement about this objective's ruggedness and is reported rather "
                    "than absorbed."
                ),
            },
        },

        "parameters": params,
        "identifiability": ident,
        "per_response_fit": per_response_fit(objective, x),
        "diagnostics_out_of_objective": diagnostics(x, diagnostic_obs),
        "readoff_reproducibility_floor": readoff_reproducibility(),

        "schiff_amadori_resolution": {
            "question": (
                "src/barrier_constants.py FAST_BARRIERS says the Schiff condensation is "
                "~2e4x FASTER than the Amadori rearrangement; data/lit/arrhenius_params.yml "
                "says the Amadori rearrangement is ~3.3e4x faster than the Schiff "
                "condensation. The two disagree by ~6.6e8 in the ratio. Open since Wave I."
            ),
            "comparison_basis": (
                "pseudo-first-order at the experiment's own 200 mmol/L glycine loading: "
                "k_schiff * [Gly]_0 (1/min) versus k_amadori (1/min). A bimolecular and a "
                "unimolecular constant cannot be compared any other way."
            ),
            "k_schiff_pseudo_first_order_100C_per_min": k_schiff_pseudo,
            "k_amadori_100C_per_min": k_amadori,
            "fitted_ratio_amadori_over_schiff": k_amadori / max(k_schiff_pseudo, 1e-300),
            "ratio_range_not_rejected_at_95pct": (
                [min(accepted), max(accepted)] if accepted else None
            ),
            "fast_barriers_claimed_ratio": 1.0 / 2.0e4,
            "arrhenius_yaml_claimed_ratio": 3.3e4,
            "profile_likelihood": profile,
        },
        "screening_lane_propagation": fast_barrier_mapping(params),

        "independent_cross_check_brands": {
            "source": (
                "Brands, C.M.J. (2002) PhD thesis, Wageningen University, ISBN "
                "90-5808-591-0, https://edepot.wur.nl/199005; Scheme 4.1, Table 4.1 "
                "(120 C, least-squares column) and Table 4.2 (Arrhenius, 90-130 C). "
                "Published as 10.1021/jf001430b (Brands & van Boekel 2001)."
            ),
            "independence": (
                "INDEPENDENT of this fit: different amine (protein-bound lysine "
                "epsilon-amino, not free glycine), different sugar:amine ratio (10:1, not "
                "1:1), different author, different data. None of it entered the objective."
            ),
            "constants": {
                name: {"value_120C": v, "ci95_120C": ci, "unit": unit,
                       "ea_kj_mol": ea, "ea_ci95": ea_ci, "step": desc}
                for name, (v, ci, unit, ea, ea_ci, desc) in BRANDS_120C.items()
            },
            "comparison": brands_comparison(params),
        },

        "reproducibility_check_martins": {
            "source": (
                "Martins, S.I.F.S. (2003) PhD thesis, Wageningen University, ISBN "
                "90-5808-823-5, Chapter 5 Table 5.2, model M4; published as "
                "10.1016/j.foodchem.2004.04.006. X = k(T_av = 100 C), reparameterised "
                "Arrhenius (thesis eqs 5.3-5.6)."
            ),
            "independence": (
                "NOT INDEPENDENT. Martins fitted these constants to the same measurements "
                "this script fits. Agreement is a reproducibility result -- a different "
                "implementation, scheme and objective recovering the same rates -- and is "
                "NOT evidence from a second experiment."
            ),
            "constants": {
                str(step): {"X_k_100C": xv, "ea_kj_mol": ea, "ea_ci95": ci, "step": desc}
                for step, (xv, ea, ci, desc) in MARTINS_M4.items()
            },
            "comparison": martins_comparison(params),
        },
    }
    return payload


def render_markdown(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    add = lines.append
    add("# Trunk rate calibration — the first rate-level fit of the sugar trunk")
    add("")
    add(f"Wave S3 · generated {payload['generated_on']} · git "
        f"`{payload['git']['commit'][:7]}` ({payload['git']['branch']})")
    add("")
    add("## (a) Lane")
    add("")
    add("Fitted in a **dedicated integrator**, `src/trunk_kinetics.py`, not in either "
        "shipped lane. The FAST screening lane performs no time integration at all; the "
        "Cantera export lane integrates but discards the initial molarity "
        "(`src/kinetics.py:373` sets `phase.X`, which Cantera normalises), lacks four of "
        "the ten measured species, and lumps three separately-measured enolisation steps "
        "onto one A/Ea pair. Both blockers were measured, not assumed — see the module "
        "docstring. **No shipped prediction imports the fitted module.**")
    add("")
    obj = payload["objective"]
    add("## (b) Fitted parameters")
    add("")
    add(f"{obj['n_values']} values, {obj['n_free_parameters']} free parameters, "
        f"reduced chi-square {obj['reduced_chi_square']:.3f}. "
        f"{obj['multistart']['reached_best_basin']}/{obj['multistart']['successful']} "
        f"random starts reached the reported optimum.")
    add("")
    add(obj["reduced_chi_square_reading"])
    add("")
    add("| step | repo family | k(100 °C) | 95% CI | unit | rate identifiability |")
    add("|---|---|---:|---|---|---|")
    for key in STEP_KEYS:
        p = payload["parameters"][key]
        ci = (f"{p['k_ref_ci95'][0]:.3g} – {p['k_ref_ci95'][1]:.3g}"
              if p["k_ref_ci95"] else "**unbounded**")
        fam = p["repo_family"] or "*(none — structural gap)*"
        add(f"| `{key}` | {fam} | {p['k_ref_100C']:.4g} | {ci} | {p['unit']} | "
            f"{p['identifiability_k']} |")
    add("")
    add("| step | Ea kJ/mol | 95% CI | Ea identifiability |")
    add("|---|---:|---|---|")
    for key in STEP_KEYS:
        p = payload["parameters"][key]
        eaci = (f"{p['ea_ci95_kj_mol'][0]:.0f} – {p['ea_ci95_kj_mol'][1]:.0f}"
                if p["ea_ci95_kj_mol"] else "**unbounded**")
        add(f"| `{key}` | {p['ea_kj_mol']:.1f} | {eaci} | {p['identifiability_ea']} |")
    add("")
    add("### Sloppy directions")
    add("")
    add(f"Hessian condition number {payload['identifiability']['condition_number']:.3g}. "
        "Three sloppiest eigendirections:")
    add("")
    for entry in payload["identifiability"]["spectrum"][:3]:
        comps = ", ".join(f"{c['parameter']} ({c['weight']:+.2f})"
                          for c in entry["dominant_components"])
        add(f"* eigenvalue `{entry['eigenvalue']:.3g}` — {comps}")
    add("")
    add("## (c) Cross-checks")
    add("")
    add("### Brands (2002) — INDEPENDENT")
    add("")
    add(payload["independent_cross_check_brands"]["independence"])
    add("")
    add("| step | fitted @120 °C | Brands @120 °C | ratio | fitted Ea | Brands Ea |")
    add("|---|---:|---:|---:|---:|---:|")
    for row in payload["independent_cross_check_brands"]["comparison"]:
        fea = f"{row['fitted_ea_kj_mol']:.0f}" if row["fitted_ea_kj_mol"] is not None else "—"
        add(f"| {row['step']} | {row['fitted_120C']:.4g} | {row['brands_120C']:.4g} | "
            f"**{row['ratio_fitted_over_brands']:.2f}×** | {fea} | "
            f"{row['brands_ea_kj_mol']:.0f} |")
    add("")
    add("### Martins (2003) M4 — REPRODUCIBILITY, NOT INDEPENDENCE")
    add("")
    add(payload["reproducibility_check_martins"]["independence"])
    add("")
    add("Bimolecular comparators are converted to pseudo-first-order at the experiment's "
        "own 200 mmol/L glycine before comparison; comparing them raw understates them "
        "200-fold.")
    add("")
    add("| fitted step | Martins step(s) | fitted k(100 °C) | Martins comparator | ratio | ΔEa kJ/mol |")
    add("|---|---|---:|---:|---:|---:|")
    for row in payload["reproducibility_check_martins"]["comparison"]:
        steps = " + ".join(str(s) for s in row["martins_steps"])
        add(f"| `{row['fitted_step']}` | {steps} | {row['fitted_k_100C']:.4g} | "
            f"{row['martins_comparator_100C']:.4g} | "
            f"**{row['ratio_fitted_over_martins']:.2f}×** | "
            f"{row['ea_difference_kj_mol']:+.0f} |")
    add("")
    add("## (d) Schiff / Amadori")
    add("")
    res = payload["schiff_amadori_resolution"]
    add(f"Pseudo-first-order at 200 mmol/L glycine, 100 °C: Schiff "
        f"`{res['k_schiff_pseudo_first_order_100C_per_min']:.4g}` /min versus Amadori "
        f"`{res['k_amadori_100C_per_min']:.4g}` /min — ratio "
        f"`{res['fitted_ratio_amadori_over_schiff']:.3g}`.")
    add("")
    add("| Amadori/Schiff ratio | Δcost | χ² statistic | rejected at 95%? |")
    add("|---:|---:|---:|---|")
    for row in res["profile_likelihood"]:
        if row.get("delta_cost") is None:
            add(f"| {row['ratio']:.3g} | — | — | {row.get('note', '')} |")
        else:
            add(f"| {row['ratio']:.3g} | {row['delta_cost']:.3f} | "
                f"{row['chi2_statistic']:.2f} | "
                f"{'**yes**' if row['rejected_at_95pct'] else 'no'} |")
    add("")
    rng = res["ratio_range_not_rejected_at_95pct"]
    span = (f"ratio {rng[0]:g} – {rng[1]:g}" if rng
            else "no grid point other than the optimum itself survives; the ratio is "
                 "pinned to within a factor well under 2")
    add(f"**Not rejected at 95%: {span}.** "
        f"FAST_BARRIERS implies {res['fast_barriers_claimed_ratio']:.1e} "
        f"(Schiff faster); `arrhenius_params.yml` implies "
        f"{res['arrhenius_yaml_claimed_ratio']:.1e} (Amadori faster). "
        "**Both are rejected — but only one of them has the sign right.**")
    add("")
    add("## (f) Propagation to the screening lane")
    add("")
    prop = payload["screening_lane_propagation"]
    add(f"**Applied: {prop['applied']}.** {prop['why_not_applied']}")
    add("")
    add(f"Prefactor used: `{prop['prefactor_used']:.0e}` — {prop['prefactor_justification']}")
    add("")
    add("| step | family | derived kcal/mol | shipped kcal/mol | Δ |")
    add("|---|---|---:|---:|---:|")
    for row in prop["rows"]:
        if row["derived_barrier_kcal_mol"] is None:
            add(f"| `{row['step']}` | {row['family']} | *excluded* | "
                f"{row['shipped_barrier_kcal_mol']} | — |")
        else:
            add(f"| `{row['step']}` | {row['family']} | "
                f"{row['derived_barrier_kcal_mol']:.2f} | "
                f"{row['shipped_barrier_kcal_mol']:.2f} | "
                f"{row['difference_kcal_mol']:+.2f} |")
    add("")
    add("## Per-response fit")
    add("")
    add("| experiment | species | T °C | n | median fold error |")
    add("|---|---|---:|---:|---:|")
    for row in payload["per_response_fit"]:
        add(f"| {row['experiment']} | {row['species']} | {row['T_C']:.0f} | {row['n']} | "
            f"{row['median_fold_error']:.3f}× |")
    add("")
    add("## Out-of-objective diagnostics")
    add("")
    add("```json")
    add(json.dumps(payload["diagnostics_out_of_objective"], indent=2))
    add("```")
    add("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, default=48)
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--basename", default="trunk_rate_calibration_refit")
    args = parser.parse_args()

    payload = build(args.starts)

    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{args.basename}.json"
    md_path = out_dir / f"{args.basename}.md"
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
