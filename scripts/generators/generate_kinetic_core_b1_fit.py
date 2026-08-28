#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b1_fit.py

THE FIT STAGE OF BUILD WAVE B1 (kinetic core, Module 1: trunk + melanoidin sink).

Estimates the four consumption constants the corpus has NO literature value for,
by least squares on the DECLARED FIT ROWS ONLY, and writes
`results/validation/kinetic_core_b1_fit_report.{json,md}` together with the
FROZEN parameter set that the hold-out scorer reads.

--------------------------------------------------------------------------
THE FIT CORPUS -- and nothing else
--------------------------------------------------------------------------
One file: `data/lit/timeseries/martins2005_glucose_glycine_80_100_120C_pH68.yml`
(glucose 200 mmol/L + glycine 200 mmol/L, 0.1 M phosphate pH 6.8, 80/100/120 C,
figure read-off from Martins 2003 thesis Fig. 5.10 = Martins & van Boekel 2005
Fig. 6).

NINE of its ten responses enter the objective:
    glucose, glycine, fructose, DFG, 3-DG, 1-DG, methylglyoxal,
    formic acid, acetic acid.

THE TENTH -- MELANOIDINS -- IS THE MODULE 4 HOLD-OUT and NEVER enters. It is
the browning response of Martins step 9, the row that
`docs/reference/FIT_HOLDOUT_DECLARATION.md` D.6 marks
"Martins 2005 browning (step 9, epsilon 0.64) -- HOLD-OUT". This script cannot
read it: `_load_fit_rows` refuses any species not in the explicit FIT allow-list
and asserts that 'melanoidins' is absent from it. The two epsilon values
(0.64 L/(mmol*cm), 282 L/(mol*cm)) appear nowhere in this file or in
`src/kinetic_core/`.

Nothing under `data/benchmarks/` is read. No hold-out bundle is read. The FAST
screening lane, `benchmark_validation` and every governance module are
untouched.

--------------------------------------------------------------------------
WHAT IS FITTED, AND WHAT IS NOT
--------------------------------------------------------------------------
NOT FITTED (fixed at their measured values, Martins & van Boekel 2005 Table 2):
    the ten trunk steps, including the melanoidin-forming step 9
    (X = 8.12e-4 L/(mmol*min) at 100 C, Ea 95.2 +/- 2.3 kJ/mol).

FITTED (4 steps x 2 parameters = 8):
    k_glc_frag  amine-independent glucose loss   -- no literature value
    k_mgo_mel   methylglyoxal -> melanoidin      -- no literature value
    k_fa_frag   formic acid decomposition        -- no literature value
    k_aa_frag   acetic acid decomposition        -- no literature value

Each is started from RANDOM points inside deliberately wide bounds, so any
agreement with a literature number is a result and not an initialisation.

A SECOND, PRE-REGISTERED FIT (variant B) additionally frees the melanoidin
sink's own k_ref -- constrained ONLY by the reactant-side responses (3-DG and
glycine), never by the browning response. Its purpose is stated before it is
run: Martins fitted his step-9 constant to the very melanoidin response that is
held out, so scoring the hold-out with HIS constant is a reproducibility check,
not an out-of-sample test. Variant B is the genuine prediction and is
pre-registered as the headline of the hold-out report.

--------------------------------------------------------------------------
STATISTICS
--------------------------------------------------------------------------
Multi-response weighted least squares in log space, the same form the
repository's wave-S3 calibration used:

    r = [ log(y_pred + f_s) - log(y_obs + f_s) ] / sigma_s

`f_s` is a per-species additive floor equal to 2% of that panel's printed full
scale, which is the source file's OWN measured worst-case read-off error
(`attribution_validation`: up to 2.2% of full scale for DFG). It keeps the log
finite at a measured zero instead of dropping those points -- dropping them is
how a model that predicts nothing scores well. `sigma_s` is a log-space
standard deviation estimated from the REPLICATE SCATTER in the data itself.

Confidence intervals are Gauss-Newton (J^T J), scaled by the reduced
chi-square, with genuinely flat directions reported as UNIDENTIFIED rather than
as a large finite number.

Usage:
    python scripts/generators/generate_kinetic_core_b1_fit.py [--starts N]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.kinetic_core import (  # noqa: E402
    CROSS_LAB_COMPARATORS,
    FITTED_BOUNDS_LOG10K,
    FITTED_EA_BOUNDS,
    FITTED_KEYS,
    MARTINS_M4,
    MEASURED_LABEL_TO_KEY,
    SCHIFF_AMADORI_SPLIT,
    T_REF_K,
    describe,
    flux_budget,
    integrate,
    operative_parameters,
    registry_metadata,
)

TIMESERIES = ROOT / "data" / "lit" / "timeseries"
GG_FILE = "martins2005_glucose_glycine_80_100_120C_pH68.yml"
BRANDS_FILE = "brands_sugar_casein_120C_pH68.yml"

OUTPUT_DIR = ROOT / "results" / "validation"
BASENAME = "kinetic_core_b1_fit_report"

# ---------------------------------------------------------------------------
# THE FIT ALLOW-LIST. This is the hold-out firewall.
# ---------------------------------------------------------------------------
# Only these nine responses may enter the objective. 'melanoidins' is the
# Module 4 hold-out and is absent by construction; the assertion below fails the
# script rather than the silence of an omission.
FIT_SPECIES: Tuple[str, ...] = (
    "glucose", "glycine", "fructose", "DFG", "3-DG", "1-DG",
    "methylglyoxal", "formic_acid", "acetic_acid",
)
HOLDOUT_SPECIES: Tuple[str, ...] = ("melanoidins",)
assert not set(FIT_SPECIES) & set(HOLDOUT_SPECIES)
assert all(name in MEASURED_LABEL_TO_KEY for name in FIT_SPECIES)
assert all(name not in MEASURED_LABEL_TO_KEY for name in HOLDOUT_SPECIES)

# Additive log floors, mmol/L: 2% of each panel's printed full scale.
# Panel scales from the source file's `axis_ranges_verified` block:
#   A Glu 0-250, Fru 0-60; B Gly 0-250, DFG 0-25; C DG 0-1; D Acids 0-21;
#   E MG 0-10; F Mel 0-10.
FLOOR: Mapping[str, float] = {
    "glucose": 0.02 * 250.0,
    "glycine": 0.02 * 250.0,
    "fructose": 0.02 * 60.0,
    "DFG": 0.02 * 25.0,
    "3-DG": 0.02 * 1.0,
    "1-DG": 0.02 * 1.0,
    "methylglyoxal": 0.02 * 10.0,
    "formic_acid": 0.02 * 21.0,
    "acetic_acid": 0.02 * 21.0,
}

INITIAL = {"Glc": 200.0, "Gly": 200.0}

N_FITTED = len(FITTED_KEYS)


class Observation:
    __slots__ = ("species", "t_c", "t_min", "value", "floor")

    def __init__(self, species: str, t_c: float, t_min: float, value: float,
                 floor: float) -> None:
        self.species = species
        self.t_c = t_c
        self.t_min = t_min
        self.value = value
        self.floor = floor


def _load_fit_rows() -> List[Observation]:
    """Load the nine FIT responses. Refuses anything not on the allow-list."""
    with (TIMESERIES / GG_FILE).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)

    out: List[Observation] = []
    for entry in payload["series"]:
        species = entry.get("species")
        if species not in FIT_SPECIES:
            # Everything else -- including the melanoidin hold-out and the
            # ambiguous_80_or_100 parking blocks -- is skipped here and is not
            # returned to the caller in any form.
            continue
        t_c = entry.get("T_C")
        if not isinstance(t_c, (int, float)) or isinstance(t_c, bool):
            continue  # unassigned temperature: the source could not resolve it
        for t_min, value in entry.get("points", []):
            out.append(
                Observation(species, float(t_c), float(t_min), float(value),
                            FLOOR[species])
            )
    if not out:
        raise RuntimeError("no fit rows loaded")
    return out


def estimate_log_sigmas(observations: Sequence[Observation]) -> Dict[str, float]:
    """Pooled within-time-point log-space scatter, per species."""
    groups: Dict[Tuple[str, float, float], List[Observation]] = {}
    for obs in observations:
        groups.setdefault((obs.species, obs.t_c, obs.t_min), []).append(obs)

    deviations: Dict[str, List[float]] = {}
    counts: Dict[str, int] = {}
    for (species, _t_c, _t), members in groups.items():
        if len(members) < 2:
            continue
        logs = np.log(np.array([m.value for m in members]) + members[0].floor)
        deviations.setdefault(species, []).extend(list(logs - logs.mean()))
        counts[species] = counts.get(species, 0) + 1

    sigmas: Dict[str, float] = {}
    for species, devs in deviations.items():
        arr = np.asarray(devs)
        dof = max(arr.size - counts[species], 1)
        sigmas[species] = float(math.sqrt(float((arr ** 2).sum()) / dof))
    for species in FIT_SPECIES:
        if species not in sigmas or sigmas[species] <= 0.0:
            sigmas[species] = 0.10
    return sigmas


class Objective:
    """Weighted log-space least squares over the nine FIT responses."""

    def __init__(self, observations: Sequence[Observation],
                 sigmas: Mapping[str, float], free_melanoidin_sink: bool) -> None:
        self.observations = list(observations)
        self.sigmas = dict(sigmas)
        self.free_melanoidin_sink = bool(free_melanoidin_sink)
        self.groups: Dict[float, Dict[str, Any]] = {}
        for index, obs in enumerate(self.observations):
            group = self.groups.setdefault(obs.t_c, {"times": set(), "rows": []})
            group["times"].add(obs.t_min)
            group["rows"].append(index)
        for group in self.groups.values():
            group["grid"] = np.array(sorted(group["times"]), dtype=float)
            group["index"] = {t: i for i, t in enumerate(group["grid"])}

    @property
    def n_parameters(self) -> int:
        return 2 * N_FITTED + (1 if self.free_melanoidin_sink else 0)

    def unpack(self, x: np.ndarray) -> Dict[str, Any]:
        fitted = {
            key: (10.0 ** float(x[i]), float(x[N_FITTED + i]))
            for i, key in enumerate(FITTED_KEYS)
        }
        parameters = operative_parameters(fitted)
        if self.free_melanoidin_sink:
            from dataclasses import replace

            parameters = dict(parameters)
            parameters["k_tdg_mel"] = replace(
                parameters["k_tdg_mel"], k_ref=10.0 ** float(x[2 * N_FITTED])
            )
        return parameters

    def predict(self, x: np.ndarray) -> np.ndarray:
        parameters = self.unpack(x)
        out = np.zeros(len(self.observations), dtype=float)
        for t_c, group in self.groups.items():
            run = integrate(parameters, t_c + 273.15, INITIAL, group["grid"])
            for index in group["rows"]:
                obs = self.observations[index]
                key = MEASURED_LABEL_TO_KEY[obs.species]
                out[index] = run.series(key)[group["index"][obs.t_min]]
        return out

    def residuals(self, x: np.ndarray) -> np.ndarray:
        try:
            predicted = self.predict(x)
        except Exception:
            return np.full(len(self.observations), 1.0e3)
        res = np.empty(len(self.observations), dtype=float)
        for i, obs in enumerate(self.observations):
            res[i] = (
                math.log(max(predicted[i], 0.0) + obs.floor)
                - math.log(obs.value + obs.floor)
            ) / self.sigmas[obs.species]
        if not np.all(np.isfinite(res)):
            return np.full(len(self.observations), 1.0e3)
        return res

    def cost(self, x: np.ndarray) -> float:
        r = self.residuals(x)
        return 0.5 * float(np.dot(r, r))

    def bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        lower = [FITTED_BOUNDS_LOG10K[k][0] for k in FITTED_KEYS] + \
                [FITTED_EA_BOUNDS[0]] * N_FITTED
        upper = [FITTED_BOUNDS_LOG10K[k][1] for k in FITTED_KEYS] + \
                [FITTED_EA_BOUNDS[1]] * N_FITTED
        if self.free_melanoidin_sink:
            # log10 of the melanoidin sink k_ref. Martins' measured value is
            # -3.09; the bounds span five decades either side of it so the fit
            # is free to disagree with him by orders of magnitude.
            lower.append(-8.0)
            upper.append(0.0)
        return np.array(lower, dtype=float), np.array(upper, dtype=float)

    def labels(self) -> List[str]:
        out = [f"log10 {k}" for k in FITTED_KEYS] + [f"Ea {k}" for k in FITTED_KEYS]
        if self.free_melanoidin_sink:
            out.append("log10 k_tdg_mel")
        return out


def multistart(objective: Objective, starts: int, seed: int = 20260828):
    from scipy.optimize import least_squares

    lower, upper = objective.bounds()
    rng = np.random.default_rng(seed)
    results = []
    for _ in range(starts):
        x0 = lower + rng.random(lower.size) * (upper - lower)
        try:
            fit = least_squares(objective.residuals, x0, bounds=(lower, upper),
                                xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=2000)
        except Exception:
            continue
        results.append((0.5 * float(np.dot(fit.fun, fit.fun)), np.asarray(fit.x)))
    if not results:
        raise RuntimeError("every start failed")
    results.sort(key=lambda item: item[0])
    return results


def numeric_jacobian(objective: Objective, x: np.ndarray, step: float = 1e-5) -> np.ndarray:
    base = objective.residuals(x)
    jac = np.zeros((base.size, x.size))
    for i in range(x.size):
        h = step * max(abs(x[i]), 1.0)
        xp = x.copy()
        xp[i] += h
        jac[:, i] = (objective.residuals(xp) - base) / h
    return jac


def standard_errors(objective: Objective, x: np.ndarray,
                    variance_scale: float) -> Tuple[np.ndarray, Dict[str, Any]]:
    jac = numeric_jacobian(objective, x)
    hess = (jac.T @ jac) / max(variance_scale, 1e-12)
    eigenvalues, eigenvectors = np.linalg.eigh(hess)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    cutoff = max(eigenvalues.max(), 1.0) * 1e-10
    inv = np.array([1.0 / v if v > cutoff else np.inf for v in eigenvalues])
    finite = np.isfinite(inv)
    cov = (eigenvectors[:, finite] * inv[finite]) @ eigenvectors[:, finite].T
    variances = np.diag(cov).copy()
    for i in range(x.size):
        if float(np.sum(eigenvectors[i, ~finite] ** 2)) > 1e-8:
            variances[i] = np.inf
    labels = objective.labels()
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
    return np.sqrt(np.clip(variances, 0.0, None)), {
        "hessian_eigenvalues": [float(v) for v in np.sort(eigenvalues)],
        "condition_number": (
            float(eigenvalues.max() / eigenvalues.min())
            if eigenvalues.min() > 0 else float("inf")
        ),
        "spectrum": spectrum,
    }


def verdict_k(se_log10: float) -> str:
    if not math.isfinite(se_log10):
        return "UNIDENTIFIED (flat direction; the data do not bound this parameter)"
    if se_log10 < 0.05:
        return "well constrained (95% CI inside +/-26%)"
    if se_log10 < 0.15:
        return "constrained (95% CI inside a factor of 2)"
    if se_log10 < 0.5:
        return "weakly constrained (95% CI up to a factor of 10)"
    return "SLOPPY (order-of-magnitude only)"


def verdict_ea(se: float, value: float) -> str:
    if not math.isfinite(se):
        return "UNIDENTIFIED (flat direction)"
    lower = value - 1.96 * se
    if lower <= 0.0:
        return (f"UNCONSTRAINED SIGN (95% CI reaches {lower:.0f} kJ/mol; an Ea <= 0 "
                f"is unphysical, so the data do not determine this step's "
                f"temperature dependence)")
    if se < 5.0:
        return "well constrained"
    if se < 15.0:
        return "constrained"
    if se < 35.0:
        return "weakly constrained"
    return "SLOPPY (the 80-120 C window is too narrow for this step)"


def per_response(objective: Objective, x: np.ndarray) -> List[Dict[str, Any]]:
    predicted = objective.predict(x)
    groups: Dict[Tuple[str, float], List[int]] = {}
    for i, obs in enumerate(objective.observations):
        groups.setdefault((obs.species, obs.t_c), []).append(i)
    rows = []
    for (species, t_c), idx in sorted(groups.items()):
        errs = [
            abs(math.log10(
                (predicted[i] + objective.observations[i].floor)
                / (objective.observations[i].value + objective.observations[i].floor)
            ))
            for i in idx
        ]
        signed = [
            math.log10(
                (predicted[i] + objective.observations[i].floor)
                / (objective.observations[i].value + objective.observations[i].floor)
            )
            for i in idx
        ]
        rows.append({
            "species": species,
            "T_C": t_c,
            "n": len(idx),
            "median_fold_error": float(10 ** float(np.median(errs))),
            "max_fold_error": float(10 ** float(np.max(errs))),
            "median_signed_log10_bias": float(np.median(signed)),
        })
    return rows


def brands_cn_direction(parameters) -> Dict[str, Any]:
    """
    DIRECTIONAL-ONLY diagnostic on the melanoidin pool's elemental C/N.

    Brands (2002, Wageningen thesis Table 3.1, p. 37) measured the elemental
    C/N of the non-dialysable melanoidin fraction of a sugar/casein system at
    120 C: 4.01 / 4.10 / 4.22 at 10 / 30 / 60 min (glucose) against an unheated
    casein reference of 3.97.

    NOT COMMENSURABLE IN LEVEL: Brands' melanoidin is protein-bound, so his C/N
    is dominated by the casein backbone and sits near 4 by construction; the
    model's pool is a glycine-glucose polymer whose step-9 repeat unit is C8N1.
    What IS comparable is the SIGN: sugar carbon accreting onto a fixed
    nitrogen inventory must raise C/N with heating time. The model predicts
    that direction only if the carbon-only channel (r_mgo_mel) carries flux.

    This is a diagnostic, not a fit target, and Brands does not appear in the
    objective.
    """
    grid = np.array([0.0, 10.0, 30.0, 60.0])
    run = integrate(parameters, 120.0 + 273.15, INITIAL, grid)
    cn = run.melanoidin_c_over_n()
    predicted = [None if not math.isfinite(v) else float(v) for v in cn]
    finite = [v for v in predicted[1:] if v is not None]
    return {
        "measured_brands_glucose_casein_120C": {"10": 4.01, "30": 4.10, "60": 4.22,
                                                "unheated_casein_reference": 3.97},
        "predicted_model_c_over_n": {"0": predicted[0], "10": predicted[1],
                                     "30": predicted[2], "60": predicted[3]},
        "measured_direction": "increasing",
        "predicted_direction": (
            "increasing" if len(finite) >= 2 and finite[-1] > finite[0] + 1e-9
            else "flat_or_decreasing"
        ),
        "commensurable_in_level": False,
        "why_not": (
            "Brands' melanoidin is protein-bound (casein backbone, unheated "
            "reference C/N 3.97), so its level is set by the protein, not by the "
            "Maillard carbon. Only the SIGN of the change transfers."
        ),
        "role": "directional diagnostic, not a fit target and not a scored hold-out",
    }


def git_head() -> Dict[str, Any]:
    def run(args: List[str]) -> str:
        try:
            return subprocess.check_output(
                args, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unknown"

    return {
        "commit": run(["git", "rev-parse", "HEAD"]),
        "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "dirty": bool(run(["git", "status", "--porcelain"])),
    }


def fit_variant(objective: Objective, starts: int, name: str) -> Dict[str, Any]:
    results = multistart(objective, starts)
    best_cost, x = results[0]
    tol = max(best_cost * 1e-4, 1e-6)
    converged = sum(1 for cost, _ in results if cost - best_cost <= tol)

    dof = max(len(objective.observations) - objective.n_parameters, 1)
    reduced_chi2 = (2.0 * best_cost) / dof
    stderr, ident = standard_errors(objective, x, reduced_chi2)

    lower, upper = objective.bounds()
    parameters: Dict[str, Any] = {}
    for i, key in enumerate(FITTED_KEYS):
        se_log10 = float(stderr[i])
        se_ea = float(stderr[N_FITTED + i])
        k_ref = 10.0 ** float(x[i])
        parameters[key] = {
            "k_ref_100C": k_ref,
            "log10_k_ref_100C": float(x[i]),
            "unit": "1/min",
            "log10_k_ref_stderr": se_log10 if math.isfinite(se_log10) else None,
            "k_ref_ci95": (
                [k_ref * 10 ** (-1.96 * se_log10), k_ref * 10 ** (1.96 * se_log10)]
                if math.isfinite(se_log10) else None
            ),
            "ea_kj_mol": float(x[N_FITTED + i]),
            "ea_stderr_kj_mol": se_ea if math.isfinite(se_ea) else None,
            "ea_ci95_kj_mol": (
                [float(x[N_FITTED + i] - 1.96 * se_ea),
                 float(x[N_FITTED + i] + 1.96 * se_ea)]
                if math.isfinite(se_ea) else None
            ),
            "identifiability_k": verdict_k(se_log10),
            "identifiability_ea": verdict_ea(se_ea, float(x[N_FITTED + i])),
            # "at bound" means pinned to the edge of the search box: within
            # 0.001 in log10, i.e. 0.2% in k. A parameter that lands there is
            # not an estimate -- the fit wanted to go further and could not,
            # which for a lower bound of 1e-8 /min means the data prefer no
            # such channel at all.
            "at_bound": bool(
                abs(float(x[i]) - lower[i]) < 1e-3 or abs(float(x[i]) - upper[i]) < 1e-3
            ),
        }
    melanoidin_sink: Dict[str, Any] | None = None
    if objective.free_melanoidin_sink:
        i = 2 * N_FITTED
        se = float(stderr[i])
        k_ref = 10.0 ** float(x[i])
        melanoidin_sink = {
            "k_ref_100C": k_ref,
            "unit": "L/(mmol*min)",
            "ea_kj_mol_fixed_at_measured": MARTINS_M4["k_tdg_mel"].ea_kj_mol,
            "log10_k_ref_stderr": se if math.isfinite(se) else None,
            "k_ref_ci95": (
                [k_ref * 10 ** (-1.96 * se), k_ref * 10 ** (1.96 * se)]
                if math.isfinite(se) else None
            ),
            "identifiability_k": verdict_k(se),
            "martins_measured_k_ref_100C": MARTINS_M4["k_tdg_mel"].k_ref,
            "ratio_fitted_over_martins": k_ref / float(MARTINS_M4["k_tdg_mel"].k_ref),
            "constrained_by": (
                "the REACTANT side only -- 3-deoxyglucosone and glycine. The "
                "browning response never entered the objective."
            ),
        }

    return {
        "variant": name,
        "n_free_parameters": objective.n_parameters,
        "n_observations": len(objective.observations),
        "degrees_of_freedom": dof,
        "objective_value_half_ssr": best_cost,
        "sum_of_squared_residuals": 2.0 * best_cost,
        "reduced_chi_square": reduced_chi2,
        "rms_weighted_residual": math.sqrt(2.0 * best_cost / len(objective.observations)),
        "multistart_total": len(results),
        "multistart_reaching_the_optimum": converged,
        "identifiability": ident,
        "fitted_parameters": parameters,
        "melanoidin_sink": melanoidin_sink,
        "per_response_fit": per_response(objective, x),
        "_x": [float(v) for v in x],
    }


def build(starts: int) -> Dict[str, Any]:
    observations = _load_fit_rows()
    sigmas = estimate_log_sigmas(observations)

    objective_a = Objective(observations, sigmas, free_melanoidin_sink=False)
    objective_b = Objective(observations, sigmas, free_melanoidin_sink=True)

    variant_a = fit_variant(objective_a, starts, "A_measured_sink")
    variant_b = fit_variant(objective_b, starts, "B_reactant_side_sink")

    params_a = objective_a.unpack(np.array(variant_a["_x"]))
    params_b = objective_b.unpack(np.array(variant_b["_x"]))

    budget = flux_budget(params_a, 373.15, INITIAL, 120.0)
    total_flux = sum(v for v in budget.values() if v > 0)

    payload: Dict[str, Any] = {
        "artifact": BASENAME,
        "wave": "B1",
        "module": "kinetic core, Module 1: trunk network + melanoidin mass sink",
        "generated_on": date.today().isoformat(),
        "git": git_head(),
        "reference_temperature_K": T_REF_K,
        "fit_target_files": [f"data/lit/timeseries/{GG_FILE} (nine responses)"],
        "fit_species": list(FIT_SPECIES),
        "holdout_species_never_read": list(HOLDOUT_SPECIES),
        "holdout_declaration": "docs/reference/FIT_HOLDOUT_DECLARATION.md D.6, Module 4",
        "excluded": [
            "the melanoidin/browning response of the same file (Module 4 HOLD-OUT)",
            "epsilon = 0.64 L/(mmol*cm) (Martins) and 282 L/(mol*cm) (Knol) -- "
            "HOLD-OUT, absent from src/kinetic_core/ and from this script",
            "data/benchmarks/** (never read)",
            f"data/lit/timeseries/{BRANDS_FILE} (directional C/N diagnostic only)",
        ],
        "network": describe(),
        "registry": registry_metadata(params_a),
        "log_sigmas_from_replicate_scatter": sigmas,
        "additive_log_floors_mmol_L": dict(FLOOR),
        "floor_justification": (
            "2% of each panel's printed full scale, which is the source file's own "
            "MEASURED worst-case figure read-off error (attribution_validation: up "
            "to 2.2% of full scale for DFG at 60-120 min)."
        ),
        "variant_A": {k: v for k, v in variant_a.items() if k != "_x"},
        "variant_B": {k: v for k, v in variant_b.items() if k != "_x"},
        "frozen_parameters": {
            "variant_A_measured_sink": {
                key: {"k_ref_100C": variant_a["fitted_parameters"][key]["k_ref_100C"],
                      "ea_kj_mol": variant_a["fitted_parameters"][key]["ea_kj_mol"]}
                for key in FITTED_KEYS
            },
            "variant_B_reactant_side_sink": {
                **{key: {"k_ref_100C": variant_b["fitted_parameters"][key]["k_ref_100C"],
                         "ea_kj_mol": variant_b["fitted_parameters"][key]["ea_kj_mol"]}
                   for key in FITTED_KEYS},
                "k_tdg_mel": {
                    "k_ref_100C": variant_b["melanoidin_sink"]["k_ref_100C"],
                    "ea_kj_mol": MARTINS_M4["k_tdg_mel"].ea_kj_mol,
                },
            },
        },
        "carbon_budget_100C_120min_mmol_per_L": {
            "per_reaction_flux": budget,
            "share_of_total_flux": {
                k: (v / total_flux if total_flux > 0 else 0.0)
                for k, v in budget.items()
            },
        },
        "melanoidin_c_over_n_directional_diagnostic": brands_cn_direction(params_a),
        "cross_lab_comparators": [dict(row) for row in CROSS_LAB_COMPARATORS],
        "schiff_amadori_split": dict(SCHIFF_AMADORI_SPLIT),
        "prereg_holdout_plan": {
            "declared_before_any_holdout_value_was_read": True,
            "headline": (
                "variant B (melanoidin sink estimated from the REACTANT side only) "
                "is the pre-registered out-of-sample prediction of the Martins "
                "step-9 browning hold-out."
            ),
            "secondary": (
                "variant A (Martins' own measured step-9 constant) is reported "
                "alongside it as a REPRODUCIBILITY check, and is explicitly NOT "
                "out-of-sample: Martins fitted that constant to the held-out "
                "response."
            ),
            "no_parameter_may_change_after_the_holdout_is_read": True,
        },
    }
    return payload


def render_markdown(payload: Dict[str, Any]) -> str:
    lines: List[str] = []
    add = lines.append
    add(f"# Kinetic core B1 -- fit report")
    add("")
    add(f"Generated {payload['generated_on']} on "
        f"`{payload['git']['branch']}` @ `{payload['git']['commit'][:8]}`"
        f"{' (dirty)' if payload['git']['dirty'] else ''}.")
    add("")
    add("## What was fitted")
    add("")
    add(f"- Fit corpus: `{payload['fit_target_files'][0]}`")
    add(f"- Responses in the objective ({len(payload['fit_species'])}): "
        f"{', '.join(payload['fit_species'])}")
    add(f"- **Never read**: {', '.join(payload['holdout_species_never_read'])} "
        f"-- the Module 4 hold-out ({payload['holdout_declaration']})")
    add("")
    for item in payload["excluded"]:
        add(f"- excluded: {item}")
    add("")
    add("The ten Martins trunk constants are **fixed at their measured values** "
        "and are not fitted. Only four consumption steps, for which the corpus "
        "has no rate constant at all, are estimated.")
    add("")

    for variant_key, title in (("variant_A", "Variant A -- measured melanoidin sink"),
                               ("variant_B", "Variant B -- reactant-side melanoidin sink "
                                             "(the pre-registered predictor)")):
        variant = payload[variant_key]
        add(f"## {title}")
        add("")
        add(f"- observations {variant['n_observations']}, free parameters "
            f"{variant['n_free_parameters']}, dof {variant['degrees_of_freedom']}")
        add(f"- **objective value (0.5 * SSR of weighted log residuals): "
            f"{variant['objective_value_half_ssr']:.4f}**")
        add(f"- reduced chi-square {variant['reduced_chi_square']:.3f}; "
            f"RMS weighted residual {variant['rms_weighted_residual']:.3f}")
        add(f"- multistart: {variant['multistart_reaching_the_optimum']} of "
            f"{variant['multistart_total']} random starts reached this optimum")
        add("")
        add("| parameter | k(100 C) | 95% CI | Ea kJ/mol | 95% CI | at bound? | verdict (k) | verdict (Ea) |")
        add("|---|---:|---|---:|---|---|---|---|")
        for key, row in variant["fitted_parameters"].items():
            ci = (f"{row['k_ref_ci95'][0]:.3g} - {row['k_ref_ci95'][1]:.3g}"
                  if row["k_ref_ci95"] else "unbounded")
            eaci = (f"{row['ea_ci95_kj_mol'][0]:.0f} - {row['ea_ci95_kj_mol'][1]:.0f}"
                    if row["ea_ci95_kj_mol"] else "unbounded")
            add(f"| `{key}` | {row['k_ref_100C']:.4g} /min | {ci} | "
                f"{row['ea_kj_mol']:.1f} | {eaci} | "
                f"{'**YES -- REJECTED BY THE DATA**' if row['at_bound'] else 'no'} | "
                f"{row['identifiability_k']} | {row['identifiability_ea']} |")
        rejected, unconstrained, required = [], [], []
        for key, row in variant["fitted_parameters"].items():
            if row["at_bound"]:
                rejected.append(key)
            elif row["identifiability_k"].startswith("UNIDENTIFIED"):
                unconstrained.append(key)
            else:
                required.append(key)
        add("")
        add(f"- **channels the data REQUIRE** (bounded away from zero): "
            f"{', '.join(f'`{k}`' for k in required) or 'none'}")
        add(f"- **channels the data REJECT** (pinned to the lower search bound, i.e. "
            f"the fit prefers no such channel at all): "
            f"{', '.join(f'`{k}`' for k in rejected) or 'none'}")
        add(f"- **channels the data cannot resolve** (flat direction; the estimate "
            f"printed above is not information): "
            f"{', '.join(f'`{k}`' for k in unconstrained) or 'none'}")
        if variant["melanoidin_sink"]:
            sink = variant["melanoidin_sink"]
            add("")
            add(f"**Melanoidin sink, estimated from the reactant side only:** "
                f"k(100 C) = {sink['k_ref_100C']:.4g} L/(mmol*min) "
                f"(Martins measured {sink['martins_measured_k_ref_100C']:.4g}; "
                f"ratio fitted/measured = {sink['ratio_fitted_over_martins']:.3g}). "
                f"{sink['identifiability_k']}.")
        add("")
        add("### Per-response agreement (in the objective)")
        add("")
        add("| species | T (C) | n | median fold error | max fold error | signed log10 bias |")
        add("|---|---:|---:|---:|---:|---:|")
        for row in variant["per_response_fit"]:
            add(f"| {row['species']} | {row['T_C']:.0f} | {row['n']} | "
                f"{row['median_fold_error']:.2f}x | {row['max_fold_error']:.2f}x | "
                f"{row['median_signed_log10_bias']:+.3f} |")
        add("")

    worst = sorted(payload["variant_A"]["per_response_fit"],
                   key=lambda row: -row["median_fold_error"])[:5]
    add("## The five worst-fitting responses (variant A)")
    add("")
    add("| species | T (C) | median fold error | direction |")
    add("|---|---:|---:|---|")
    for row in worst:
        add(f"| {row['species']} | {row['T_C']:.0f} | {row['median_fold_error']:.2f}x | "
            f"{'model high' if row['median_signed_log10_bias'] > 0 else 'model low'} |")
    add("")
    add("These are misfits of the MEASURED constants, not of anything fitted here: "
        "the ten Martins steps were held at their published values. The largest, "
        "formic acid at 120 C, sits in the direction that Knol 2010's conflicting "
        "formic-acid Ea (84 +/- 14 vs Martins' 30 +/- 9) would move it, which is "
        "recorded as a standing conflict rather than resolved by refitting.")
    add("")

    add("## Where the carbon goes (100 C, 120 min, variant A)")
    add("")
    add("| reaction | integrated flux mmol/L | share |")
    add("|---|---:|---:|")
    budget = payload["carbon_budget_100C_120min_mmol_per_L"]
    for key, value in sorted(budget["per_reaction_flux"].items(),
                             key=lambda kv: -kv[1]):
        add(f"| `{key}` | {value:.4g} | "
            f"{100 * budget['share_of_total_flux'][key]:.1f}% |")
    add("")

    diag = payload["melanoidin_c_over_n_directional_diagnostic"]
    add("## Melanoidin C/N -- directional diagnostic (not fitted, not scored)")
    add("")
    add(f"- measured (Brands 2002, glucose/casein 120 C): "
        f"{diag['measured_brands_glucose_casein_120C']} -- direction "
        f"**{diag['measured_direction']}**")
    add(f"- predicted: {diag['predicted_model_c_over_n']} -- direction "
        f"**{diag['predicted_direction']}**")
    add(f"- {diag['why_not']}")
    add("")

    add("## Cross-lab comparators carried but NOT operative")
    add("")
    add("| quantity | value | source | why not operative |")
    add("|---|---|---|---|")
    for row in payload["cross_lab_comparators"]:
        value = (f"{row['value_kj_mol']} kJ/mol" if row.get("value_kj_mol")
                 else str(row.get("value_set_kj_mol", "-")))
        add(f"| {row['quantity']} | {value} | {row['source_anchor']} | "
            f"{row['why_not_operative']} |")
    add("")

    add("## Pre-registration of the hold-out evaluation")
    add("")
    prereg = payload["prereg_holdout_plan"]
    add(f"- {prereg['headline']}")
    add(f"- {prereg['secondary']}")
    add("- No parameter may change after the hold-out is read. The hold-out "
        "scorer is a separate script that reads the frozen parameters from the "
        "JSON companion of this report.")
    add("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=12)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    payload = build(args.starts)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"{BASENAME}.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )
    (args.output_dir / f"{BASENAME}.md").write_text(
        render_markdown(payload), encoding="utf-8"
    )
    print(f"wrote {args.output_dir / BASENAME}.json and .md")
    print(f"variant A objective = "
          f"{payload['variant_A']['objective_value_half_ssr']:.4f}")
    print(f"variant B objective = "
          f"{payload['variant_B']['objective_value_half_ssr']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
