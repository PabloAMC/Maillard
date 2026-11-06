#!/usr/bin/env python
"""
build_s1_profiles.py -- the profile likelihoods S1 asks for (plan sec. 7 rule S1; owner 2026-11-06).

WHY. The S1 audit (``results/validation/incumbent_evidence_tiers.json``) leaves 28 fitted coordinates at
``fitted_unidentified`` with ``identification.status = locally_pinned_profile_missing``: every local criterion passes
(no band edge, no flat slice, a CI inside the search box, the fit's own label), but S1's fitted-identified tier needs
"a converged fit whose MARGINAL profile is pinned with every other unidentified coordinate profiled over its band",
and no artifact on disk profiled them (slice profiles hold the others at the optimum; a Laplace sigma is local). This
generator computes that profile.

WHAT. For each of the 28 coordinates, on the fit that froze its SHIPPED value, with that fit's own residual function
(imported from the frozen wave generator, never edited), rows, sigmas, free set and search bounds:

1. re-fit every free coordinate from the frozen point with the engine's fit driver (``src.kinetic_core.fitting.fit``:
   the codebase's trf settings, its forward Jacobian with absolute per-kind steps, a central-difference polish,
   asserted convergence flagged), giving the optimum x_hat and S_min = r . r of today's engine on those rows;
2. walk the coordinate outward from x_hat on both sides (steps from its Laplace interval, growing), and at each grid
   value psi re-optimise EVERY OTHER free coordinate of that fit (warm start from the neighbouring point): the profile
   S_p(psi) = min over the others of r . r;
3. read the 95 % interval by likelihood ratio, Delta(psi) = S_p(psi) - S_min, three calibrations:
   * ``raw``     Delta <= chi2_1(0.95) = 3.841 (Wilks, the rows' sigmas as stated);
   * ``scaled``  Delta / max(1, chi2_red) <= 3.841 (the residual variance inflated by the fit's reduced chi-square
                 when it exceeds one, the convention of every fit's own Laplace covariance; THE VERDICT USES THIS);
   * ``F``       Delta / chi2_red <= F(0.95; 1, dof) (the variance estimated from the dof the fit leaves; reported as
                 a small-sample check when chi2_red > 1);
   edges by signed-root interpolation (exact for a quadratic profile). The calibration is ASYMPTOTIC (Wilks; regular,
   interior optimum, Gaussian residuals at the stated sigmas); it is not checked by simulation here;
4. judge: pinned on both sides by the likelihood (the scaled interval closes strictly inside the search bounds, which
   never count as a band, S4), the refit and every profile point inside the interval converged with no failed pot,
   and, where a measured or class band is on record, the interval inside it. None of the 28 has one: each row's band
   is the fit's own Laplace sigma (a fit output), the search bounds, a fixed value or nothing, so band containment is
   recorded as undefined and the closure is the likelihood's alone.

WHAT IT IS NOT. Profiling a fit on its own rows is IDENTIFICATION, never credit (item 2, spent data). It changes no
tier: ``docs/reference/evidence_tier_changes.yml`` is the register, and an upgrade needs this artifact AND a second
reader. Sulfur-lane coordinates (B60) and the pH-drift pKa run on the sulfur balance B0, falsified under S3 (PR #94,
M1b part (i), amendment A2): their profiles are CONDITIONAL on a falsified structure and propose no upgrade. Only the
fit's own free coordinates are profiled: a constant the fit held fixed (a measured value, an earlier fit's output)
stays fixed, and each fit's record lists what it held (``held_fixed``), for the second reader.

THE FED 3-DEOXY TRIANGLE. The S1 record attributes k_tdg_ddg, k_ddg_tdg, k_ddg_dgal, k_dgal_ddg and k_ddg_hmf to B39,
but the shipped values are B41's optimum (``parameters_dicarbonyl.SHIPPING_FIT_REPORT``: B39's rows and objective
re-run with the formic-acid exit's pH term on). B41 froze them, so B41 is profiled; the record's attribution is
reported as a finding.

JAX. No fit here has its objective ported to ``src/kinetic_core_jax/fit.py`` (a ``JaxProblem`` needs each fit's rows
as network observables; the ports so far are the lane networks and their own equivalence problems), so every profile
runs on the engine's fit driver.

Usage (one fit per call; each merges its coordinates into the artifact and re-renders the markdown):
    python scripts/generators/build_s1_profiles.py --fit B20 [--workers 4]
    python scripts/generators/build_s1_profiles.py --fit B60 --time-point          # time one profile point, project
    python scripts/generators/build_s1_profiles.py --render                        # markdown from the json
    python scripts/generators/build_s1_profiles.py --check                         # exit 1 if the markdown is stale
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import multiprocessing as mp
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
GENERATORS = ROOT / "scripts" / "generators"

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fitting  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT_JSON = V / "s1_profiles.json"
OUT_MD = V / "s1_profiles.md"
TIERS_JSON = V / "incumbent_evidence_tiers.json"
PROFILE_MISSING_TEXT = "the profile requirement is unmet"

CHI2_1_95 = 3.841458820694124
#: the conditional label for every coordinate whose fit runs on the sulfur balance B0
CONDITIONAL_B0 = ("conditional on B0 (falsified under S3: PR #94, M1b part (i), amendment A2); no upgrade proposed")
SPENT_DATA = ("identification only, never credit: the profile is on the fit's own rows (item 2, spent data)")
RULE_S1 = ("fitted-identified (a converged fit whose marginal profile is pinned with every other unidentified "
           "coordinate profiled over its band; any band-edge, flat-profile or band-spanning-CI criterion overrides a "
           "reported stderr; ...)")
RULE_CODE = ("src/evidence_tiers.py _fitted_tier: status locally_pinned_profile_missing = 'no artifact profiles this "
             "coordinate with the other unidentified coordinates profiled over their bands (slice profiles hold them "
             "at the optimum; Laplace is local)'")

#: The 28 coordinates: S1 row id -> (fit profiled, that fit's coordinate name). The row set is asserted equal to the
#: S1 record's rows whose basis says the profile requirement is unmet.
COORDINATES: Dict[str, Tuple[str, str]] = {
    "k_mgo_mel.log10_k_ref": ("B1", "log10_k_mgo_mel"),
    "k_aa_frag.log10_k_ref": ("B1", "log10_k_aa_frag"),
    "k_tdg_ddg.log10_k_ref": ("B41", "log10_k_tdg_ddg_100C"),
    "k_ddg_tdg.log10_k_ref": ("B41", "log10_k_ddg_tdg_100C"),
    "k_ddg_dgal.log10_k_ref": ("B41", "log10_k_ddg_dgal_100C"),
    "k_dgal_ddg.log10_k_ref": ("B41", "log10_k_dgal_ddg_100C"),
    "k_ddg_hmf.log10_k_ref": ("B41", "log10_k_ddg_hmf_100C"),
    "k_mgo_ak.log10_k_ref": ("B18", "log10_k_mgo_ak_100C"),
    "k_glyc.log10_k_ref": ("B20", "log10_k_glyc_100C"),
    "k_flp_cml.log10_k_ref": ("B20", "log10_k_flp_cml_100C"),
    "k_flp_cel.log10_k_ref": ("B20", "log10_k_flp_cel_100C"),
    "k_flp_decay.log10_k_ref": ("B20", "log10_k_flp_decay_100C"),
    "k_cml_loss.log10_k_ref": ("B20", "log10_k_cml_loss_100C"),
    "k_acr_dp.log10_k_ref": ("B3", "k_acr_dp"),
    "k_acr_dp.ea_kj_mol": ("B3", "Ea_acr_dp"),
    "b18.ph_slope_above_7_decades_per_unit": ("B18", "ph_slope_above_7_decades_per_unit"),
    "b18.ph_slope_below_7_decades_per_unit": ("B18", "ph_slope_below_7_decades_per_unit"),
    "k_arp_tdp.log10_k_ref": ("B60", "k_arp_tdp"),
    "k_tdp_fur.log10_k_ref": ("B60", "k_tdp_fur"),
    "k_nf_mft.log10_k_ref": ("B60", "k_nf_mft"),
    "k_nf_mp3p.log10_k_ref": ("B60", "k_nf_mp3p"),
    "k_mgo_mp.log10_k_ref": ("B60", "k_mgo_mp"),
    "k_fft_decay.log10_k_ref": ("B60", "k_fft_decay"),
    "k_osone_decay.log10_k_ref": ("B60", "k_osone_decay"),
    "k_pent_caramel.log10_k_ref": ("B60", "k_pent_caramel"),
    "k_cys_thermal.log10_k_ref": ("B60", "k_cys_thermal"),
    "k_thiolate_loss.log10_k_ref": ("B60", "k_thiolate_loss"),
    "b8.ph_drift.arp_secondary_ammonium_pKa": ("B60", "ph_arp_secondary_ammonium_pKa"),
}
#: the work order (task brief): trunk -> acrylamide -> b18 -> sulfur
GROUP_OF = {cid: ("sulfur" if fit == "B60" else "acrylamide" if fit == "B3"
                  else "b18" if cid.startswith("b18.") else "trunk") for cid, (fit, _) in COORDINATES.items()}
GROUP_ORDER = ("trunk", "acrylamide", "b18", "sulfur")


# =====================================================================================================================
# 1. The profile (pure: a residual function, bounds and kinds; no chemistry)
# =====================================================================================================================

@dataclass
class FitSpec:
    """One frozen fit: its free vector, bounds, kinds, frozen point and residual function r(x)."""

    name: str
    keys: Tuple[str, ...]
    kinds: Tuple[str, ...]
    lower: np.ndarray
    upper: np.ndarray
    x_frozen: np.ndarray
    residuals: Callable[[np.ndarray], np.ndarray]
    n_rows: int
    meta: Dict[str, Any] = field(default_factory=dict)
    #: fitting.fit settings for every (re)fit of this spec
    budget: int = 2000
    central_polish: int = 2
    polish_budget: int = fitting.POLISH_BUDGET
    #: the failed-evaluation residual of this fit's own residual function (rows at it are counted as failed pots)
    failed_values: Tuple[float, ...] = (25.0,)

    @property
    def dof(self) -> int:
        return max(1, self.n_rows - len(self.keys))


def _sub_problem(spec: FitSpec, fixed: Mapping[int, float]) -> Tuple[fitting.Problem, List[int]]:
    """The fit over every free coordinate not in ``fixed`` (index -> value), the others at their fixed values."""
    free = [j for j in range(len(spec.keys)) if j not in fixed]
    base = np.asarray(spec.x_frozen, dtype=float).copy()
    for j, v in fixed.items():
        base[j] = v

    def residuals(theta, rows):
        x = base.copy()
        x[free] = theta
        return spec.residuals(x)

    prob = fitting.Problem(keys=tuple(spec.keys[j] for j in free), kinds=tuple(spec.kinds[j] for j in free),
                           lower=np.asarray(spec.lower, dtype=float)[free], upper=np.asarray(spec.upper, dtype=float)[free],
                           rows=tuple(range(spec.n_rows)), residuals=residuals,
                           failed_residual=float(spec.failed_values[0]))
    return prob, free


def _n_failed(spec: FitSpec, r: np.ndarray) -> int:
    return int(sum(int(np.sum(np.isclose(np.abs(r), v))) for v in spec.failed_values))


def refit(spec: FitSpec, fixed: Mapping[int, float], start_full: np.ndarray) -> Dict[str, Any]:
    """Re-optimise every coordinate not in ``fixed`` from ``start_full``; returns the full point and its record."""
    prob, free = _sub_problem(spec, fixed)
    x = np.asarray(start_full, dtype=float).copy()
    for j, v in fixed.items():
        x[j] = v
    t0 = time.time()
    if free:
        f = fitting.fit(prob, x[free], budget=spec.budget, central_polish=spec.central_polish,
                        polish_budget=spec.polish_budget, require_converged="flag")
        x[free] = f.theta
        r = np.asarray(f.residuals, dtype=float)
        pol = f.polish or {}
        rec = {"status": f.status, "converged": bool(f.converged), "evaluations": f.evaluations,
               "polish_flagged": bool(pol.get("flagged", False)),
               "polish": {k: pol.get(k) for k in ("carried", "exhausted", "raised_cost", "forward_cost",
                                                  "polished_cost")} if pol else None}
        # the profile's convergence: the optimiser stopped on a tolerance (status > 0) and the central polish was not
        # exhausted. A polish that RAISED the cost is not a failure here: the forward endpoint, itself stopped on a
        # tolerance, is carried (fitting.fit's rule) and the polish found nothing lower; the driver's own flag
        # ("converged") is reported beside it.
        rec["converged_profile"] = bool(f.status > 0 and not pol.get("exhausted", False))
    else:
        r = np.asarray(spec.residuals(x), dtype=float)
        rec = {"status": 1, "converged": True, "converged_profile": True, "evaluations": 1, "polish_flagged": False}
    rec.update({"x": [float(v) for v in x], "S": float(r @ r), "failed_rows": _n_failed(spec, r),
                "seconds": round(time.time() - t0, 2)})
    return rec


def offsets_for(x0: float, lo: float, hi: float, h0: float, direction: int, max_points: int = 16) -> List[float]:
    """Grid values from x0 toward one bound: steps h0 x (1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, ...), the bound itself
    last when the walk reaches it."""
    mult, m = [], 1.0
    while len(mult) < max_points - 1:
        mult.append(m)
        m = m + 1.0 if m < 6 else math.ceil(m * 1.3)
    edge = hi if direction > 0 else lo
    out = []
    for k in mult:
        v = x0 + direction * k * h0
        if (direction > 0 and v >= edge - 1e-12) or (direction < 0 and v <= edge + 1e-12):
            break
        out.append(v)
    out.append(edge)
    return out


def walk(spec: FitSpec, index: int, x_hat: np.ndarray, s_min: float, direction: int, h0: float,
         stop_delta: float, beyond: int = 1, max_points: int = 16) -> List[Dict[str, Any]]:
    """One side of the profile: from x_hat outward, re-fitting the others at each value (warm start), until the
    profile rises ``stop_delta`` above ``s_min`` (plus ``beyond`` more points) or the search bound is reached."""
    grid = offsets_for(float(x_hat[index]), float(spec.lower[index]), float(spec.upper[index]), h0, direction,
                       max_points)
    points, start, extra = [], np.asarray(x_hat, dtype=float).copy(), None
    for v in grid:
        rec = refit(spec, {index: float(v)}, start)
        if not rec["converged_profile"] or rec["failed_rows"]:
            # one retry with twice the polish (calls and budget) from the same warm start; the lower S is kept
            from dataclasses import replace

            again = refit(replace(spec, central_polish=max(2, 2 * spec.central_polish),
                                  polish_budget=2 * spec.polish_budget), {index: float(v)}, start)
            first = {k: rec[k] for k in ("S", "status", "converged", "converged_profile", "polish")}
            if (again["converged_profile"] and not again["failed_rows"]) or again["S"] < rec["S"]:
                rec = again
            rec["retried"] = {"first_attempt": first}
        rec["value"] = float(v)
        rec["delta"] = rec["S"] - s_min
        points.append(rec)
        start = np.asarray(rec["x"], dtype=float)
        if extra is None and rec["delta"] > stop_delta:
            extra = beyond
        if extra is not None:
            if extra <= 0:
                break
            extra -= 1
    return points


def _crossing(v0: float, d0: float, side: Sequence[Tuple[float, float]], threshold: float) -> Tuple[float, bool]:
    """The first outward crossing of ``threshold`` by Delta, by signed-root interpolation (exact for a quadratic);
    (edge, True) or (last walked value, False) when the side never crosses."""
    pv, pd = v0, max(d0, 0.0)
    for v, d in side:
        if d > threshold:
            a, b = math.sqrt(pd), math.sqrt(d)
            t = (math.sqrt(threshold) - a) / (b - a) if b > a else 1.0
            return pv + t * (v - pv), True
        pv, pd = v, max(d, 0.0)
    return pv, False


def interval(v0: float, low_side: Sequence[Tuple[float, float]], high_side: Sequence[Tuple[float, float]],
             threshold: float, d0: float = 0.0) -> Dict[str, Any]:
    """The likelihood-ratio interval {psi : Delta(psi) <= threshold} connected to v0, from the two walked sides
    (each ordered outward)."""
    lo, lo_closed = _crossing(v0, d0, low_side, threshold)
    hi, hi_closed = _crossing(v0, d0, high_side, threshold)
    return {"threshold_delta": threshold, "lower": lo, "upper": hi, "closed_lower": lo_closed,
            "closed_upper": hi_closed, "half_width": 0.5 * (hi - lo)}


#: the spec a forked chain reads (a FitSpec's residual function is a closure and does not pickle; a forked worker
#: inherits this global, so only plain arguments travel)
_ACTIVE_SPEC: Optional[FitSpec] = None


def chain_task(args) -> Tuple[int, int, List[Dict[str, Any]]]:
    index, x_hat, s_min, direction, h0, stop_delta, beyond, max_points = args
    return index, direction, walk(_ACTIVE_SPEC, index, x_hat, s_min, direction, h0, stop_delta, beyond, max_points)


def _map(fn, items, workers: int):
    if workers <= 1 or len(items) <= 1:
        return [fn(it) for it in items]
    with mp.get_context("fork").Pool(processes=min(workers, len(items))) as pool:
        return pool.map(fn, items, chunksize=1)


def thresholds(chi2_red: float, dof: int) -> Dict[str, Optional[float]]:
    from scipy.stats import f as f_dist

    scale = max(1.0, chi2_red)
    out = {"raw": CHI2_1_95, "scaled": CHI2_1_95 * scale,
           "F": float(f_dist.ppf(0.95, 1, dof)) * chi2_red if chi2_red > 1.0 else None}
    return out


def profile_fit(spec: FitSpec, indices: Sequence[int], h0: Mapping[int, float], workers: int = 1,
                beyond: int = 1, max_points: int = 16) -> Dict[str, Any]:
    """The refit, then every requested coordinate's two-sided profile, with the three intervals."""
    t0 = time.time()
    full = refit(spec, {}, spec.x_frozen)
    x_hat, s_hat = np.asarray(full["x"]), full["S"]
    chi2_red = s_hat / spec.dof
    thr = thresholds(chi2_red, spec.dof)
    stop = max(v for v in thr.values() if v is not None)
    global _ACTIVE_SPEC
    _ACTIVE_SPEC = spec
    tasks = [(i, x_hat, s_hat, d, float(h0[i]), stop, beyond, max_points) for i in indices for d in (-1, +1)]
    chains = _map(chain_task, tasks, workers)
    sides = {(i, d): pts for i, d, pts in chains}
    s_min = min([s_hat] + [p["S"] for pts in sides.values() for p in pts])
    r_frozen = np.asarray(spec.residuals(np.asarray(spec.x_frozen, dtype=float)), dtype=float)
    out = {"refit": {**full, "x_frozen": [float(v) for v in spec.x_frozen], "S_at_frozen": float(r_frozen @ r_frozen)},
           "S_min": s_min, "S_refit": s_hat, "lower_optimum_found_by_profile": bool(s_min < s_hat - 1e-6 * max(1, s_hat)),
           "chi2_red": s_min / spec.dof, "dof": spec.dof, "coordinates": {}, "seconds": round(time.time() - t0, 1)}
    thr = thresholds(s_min / spec.dof, spec.dof)
    out["thresholds_delta"] = thr
    for i in indices:
        low = sides[(i, -1)]
        high = sides[(i, +1)]
        for p in low + high:
            p["delta"] = p["S"] - s_min
        d0 = s_hat - s_min
        lo_pairs = [(p["value"], p["delta"]) for p in low]
        hi_pairs = [(p["value"], p["delta"]) for p in high]
        ivs = {name: (None if t is None else interval(float(x_hat[i]), lo_pairs, hi_pairs, t, d0))
               for name, t in thr.items()}
        out["coordinates"][spec.keys[i]] = {"index": i, "x_hat": float(x_hat[i]), "h0": float(h0[i]),
                                            "low_side": low, "high_side": high, "intervals": ivs}
    return out


# =====================================================================================================================
# 2. The fits (each imports its frozen wave generator; nothing in a wave generator is edited)
# =====================================================================================================================

def _load(stem: str):
    spec = importlib.util.spec_from_file_location(f"s1p_{stem}", GENERATORS / f"{stem}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _report(stem: str) -> dict:
    return json.loads((V / f"{stem}.json").read_text(encoding="utf-8"))


def spec_b1() -> FitSpec:
    m = _load("generate_kinetic_core_b1_fit")
    from src.kinetic_core import FITTED_KEYS

    obs = m._load_fit_rows()
    obj = m.Objective(obs, m.estimate_log_sigmas(obs), free_melanoidin_sink=False)
    lower, upper = obj.bounds()
    fp = _report("kinetic_core_b1_fit_report")["variant_A"]["fitted_parameters"]
    x = np.array([fp[k]["log10_k_ref_100C"] for k in FITTED_KEYS] + [fp[k]["ea_kj_mol"] for k in FITTED_KEYS])
    keys = tuple(f"log10_{k}" for k in FITTED_KEYS) + tuple(f"ea_{k}" for k in FITTED_KEYS)
    return FitSpec("B1", keys, ("log10k",) * len(FITTED_KEYS) + ("Ea",) * len(FITTED_KEYS), lower, upper, x,
                   obj.residuals, len(obs), failed_values=(1.0e3,),
                   meta={"report": "results/validation/kinetic_core_b1_fit_report.json (variant_A)",
                         "generator": "scripts/generators/generate_kinetic_core_b1_fit.py (Objective, variant A)",
                         "objective": "B1 variant A: log-space residuals with additive floors over Martins 2005's nine "
                                      "fit responses at 80/100/120 C; sigmas from the replicate scatter",
                         "held_fixed": "Martins & van Boekel 2005's ten measured trunk steps (incl. the melanoidin "
                                       "step 9), fixed at their measured values"})


def spec_b41() -> FitSpec:
    from src.kinetic_core import config
    from src.kinetic_core import trunk_conditions as TC

    config.update({"trunk.three_deoxy_exit_ph_term": True,
                   "trunk.three_deoxy_exit_ph": {"k_tdg_fa": TC.THREE_DEOXY_EXIT_PH["k_tdg_fa"]}})
    m = _load("generate_kinetic_core_b39_fit")
    rep = _report("kinetic_core_b41_fit_report")
    assert rep["quick"] is False
    x = np.array([rep["frozen_parameters"]["fed_3deoxy"][k] for k in m.KEYS])
    return FitSpec("B41", tuple(m.KEYS), ("log10k",) * len(m.KEYS), np.asarray(m.LOWER), np.asarray(m.UPPER), x,
                   lambda v: m.residuals(np.asarray(v, dtype=float), False), len(m.DECISIVE),
                   meta={"report": "results/validation/kinetic_core_b41_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b39_fit.py residuals (quick=False), "
                                      "run as generate_kinetic_core_b41_fit.py runs it (the formic-acid exit's pH "
                                      "term on)",
                         "objective": "B39's twelve decisive rows (Mittelmaier 2011 fed maxima, times of maximum and "
                                      "shares; Zhang 2021 within-study ratios) on log10(pred/target)/sigma",
                         "attribution_note": "the S1 record names B39; the shipped values are B41's optimum "
                                             "(parameters_dicarbonyl.SHIPPING_FIT_REPORT)",
                         "held_fixed": "k_tdg_fa (a B1 fit row, Martins), k_tdg_mgo (Kocadagli), k_glc_tdg (the "
                                       "glucose entry); the declared barrier of the three new steps; every other trunk "
                                       "constant at its shipped value"})


def spec_b18() -> FitSpec:
    m = _load("generate_kinetic_core_b18_fit")
    fr = _report("kinetic_core_b18_fit_report")["frozen_parameters"]["pyrazine"]
    x = np.array([fr[k] for k in m.KEYS])
    kinds = ("log10k", "Ea", "log10k", "Ea", "ph_slope", "ph_slope")
    return FitSpec("B18", tuple(m.KEYS), kinds, np.asarray(m.LOWER), np.asarray(m.UPPER), x, m.residuals, len(m.ROWS),
                   meta={"report": "results/validation/kinetic_core_b18_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b18_fit.py residuals",
                         "objective": "Zhou 2024's six mean formation rates and Leahy 1989's four within-study pH "
                                      "ratios, log10(pred/target)/sigma",
                         "held_fixed": "the declared aminoketone condensation (K_COND); every trunk constant at its "
                                       "B1-operative value (incl. B1's fitted k_mgo_mel, which consumes MGO in the "
                                       "fed-MGO pots)"})


def spec_b20() -> FitSpec:
    m = _load("generate_kinetic_core_b20_fit")
    fr = _report("kinetic_core_b20_fit_report")["frozen_parameters"]["glycation"]
    x = np.array([fr[k] for k in m.KEYS])
    return FitSpec("B20", tuple(m.KEYS), ("log10k",) * len(m.KEYS), np.asarray(m.LOWER), np.asarray(m.UPPER), x,
                   m.residuals, len(m.ROWS),
                   meta={"report": "results/validation/kinetic_core_b20_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b20_fit.py residuals",
                         "objective": "Nguyen 2016 Table 1 (M1)'s ten printed rate constants at 120 / 130 C, "
                                      "log10(pred/target)/sigma; no integration (each row is one constant at its "
                                      "temperature through its declared barrier)",
                         "held_fixed": "the five declared barriers; no other constant enters these rows",
                         "rows_withdrawn_note": "every B20 row is flagged in the S1 record: fitted under a network "
                                                "the source contradicts (CF-072)"})


def spec_b3() -> FitSpec:
    m = _load("generate_kinetic_core_b3_fit")
    from src.kinetic_core.parameters_acrylamide import FITTED_ACRYLAMIDE_BOUNDS_LOG10K, FITTED_ACRYLAMIDE_EA_BOUNDS

    fr = _report("kinetic_core_b3_fit_report")["frozen_parameters"]
    keys = tuple(m.PARAM_ORDER) + tuple(m.EA_ORDER)
    x = np.array([fr["log10_k_ref_at_160C"][k] for k in m.PARAM_ORDER] + [fr["fitted_Ea_kJ_mol"][k] for k in m.EA_ORDER])
    lower = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][0] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[0]] * m.N_EA)
    upper = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][1] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[1]] * m.N_EA)
    m._progress = lambda r: None  # the wave's progress printer (stdout only)
    return FitSpec("B3", keys, ("log10k",) * m.N_K + ("Ea",) * m.N_EA, lower, upper, x,
                   lambda v: m.residuals(np.asarray(v, dtype=float), True), len(m.ACTIVE_FIT_ROWS),
                   meta={"report": "results/validation/kinetic_core_b3_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b3_fit.py residuals (quick=True, "
                                      "the setting the frozen fit ran at)",
                         "objective": "B3's 30 declared fit rows (Claeys 2005, De Vleeschouwer 2009, Knol 2005, the "
                                      "SPI extrusion benchmark): log10 rate / ppb rows, linear Ea rows, one-sided "
                                      "ceiling",
                         "held_fixed": "B1's four trunk fits (inherited frozen, incl. k_mgo_mel and k_aa_frag); the "
                                       "measured acrylamide constants; the scavenging channels' measured Ea_E2",
                         "rows_withdrawn_note": "the SPI extrusion row (spi_extrusion_130C_acrylamide_ppb) is "
                                                "withdrawn in the S1 record (CF-131) but is in B3's frozen objective; "
                                                "the profile is of the objective as frozen"})


def spec_b60() -> FitSpec:
    m = _load("generate_kinetic_core_b60_fit")
    from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402  (on sys.path via the B60 module)

    rep = _report("kinetic_core_b60_fit_report")
    lower_full, upper_full = m.full_bounds()
    x_full = np.clip(frozen_vector(rep), lower_full, upper_full)
    idx = list(m.FREE_INDEX)
    keys = tuple(m.ALL_KEYS[i] for i in idx)

    def kind(k: str) -> str:
        return ("Ea" if k.startswith("Ea_") else "acid_yield" if k == "ph_acid_yield_per_sink_event"
                else "pka" if k == "ph_arp_secondary_ammonium_pKa" else "log10k")

    def residuals(v):
        xf = x_full.copy()
        xf[idx] = v
        return m.residual_vector(xf, True)

    return FitSpec("B60", keys, tuple(kind(k) for k in keys), lower_full[idx], upper_full[idx], x_full[idx], residuals,
                   len(m.B23.ACTIVE_FIT_ROWS), central_polish=0, budget=600,
                   meta={"report": "results/validation/kinetic_core_b60_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b60_fit.py residual_vector (quick, "
                                      "B60's verified quick setting)",
                         "objective": "B51's objective on B47's 52 corrected rows, B60's pH-node grid",
                         "held_fixed": "B60's 25 frozen sulfur coordinates (frozen_keys), incl. Ea_lumped_formation",
                         "conditional": CONDITIONAL_B0})


SPECS: Dict[str, Callable[[], FitSpec]] = {"B1": spec_b1, "B41": spec_b41, "B18": spec_b18, "B20": spec_b20,
                                           "B3": spec_b3, "B60": spec_b60}
#: S = r . r at each fit's frozen point as its own report records it (B1, B3 and B60 report 0.5 r . r)
REPORTED_S: Dict[str, Callable[[], float]] = {
    "B1": lambda: 2.0 * _report("kinetic_core_b1_fit_report")["variant_A"]["objective_value_half_ssr"],
    "B41": lambda: min(m["cost"] for m in _report("kinetic_core_b41_fit_report")["members"]),
    "B18": lambda: _report("kinetic_core_b18_fit_report")["objective"]["final_cost"],
    "B20": lambda: min(m["cost"] for m in _report("kinetic_core_b20_fit_report")["members"]),
    "B3": lambda: 2.0 * _report("kinetic_core_b3_fit_report")["objective"]["final_cost"],
    "B60": lambda: 2.0 * _report("kinetic_core_b60_fit_report")["objective"]["final_cost"],
}
#: a frozen objective counts as reproduced when S at the frozen point is within this fraction of its report's (B3's
#: quick integrator setting alone moves it 0.1 %)
REPRODUCED_REL = 0.01
#: the smallest first step per kind (a coordinate's first grid step is its Laplace 95 % half-width / 2.5, at least this)
MIN_STEP = {"log10k": 0.02, "Ea": 1.0, "ph_slope": 0.01, "pka": 0.05, "acid_yield": 0.005}


# =====================================================================================================================
# 3. Verdicts
# =====================================================================================================================

def tier_rows() -> Dict[str, dict]:
    rows = json.loads(TIERS_JSON.read_text(encoding="utf-8"))["rows"]
    sel = {r["id"]: r for r in rows if PROFILE_MISSING_TEXT in (r.get("basis") or "")}
    assert set(sel) == set(COORDINATES), (sorted(set(sel) ^ set(COORDINATES)))
    return sel


def s4_band(row: dict) -> Tuple[Optional[Tuple[float, float]], str]:
    """A measured or class band for the coordinate, or None and why not (S4: a search bound is never a band)."""
    b = row.get("band")
    ci = (row.get("identification") or {}).get("ci95_halfwidth")
    if not b:
        return None, "no band on record"
    if row.get("band_is_search_bound"):
        return None, f"{b.get('source')}: samples the fit's search bounds (never a band, S4)"
    if b.get("distribution") == "fixed":
        return None, f"{b.get('source')}: held fixed by the envelope (no band)"
    if b.get("sigma") is not None and ci is not None and abs(1.96 * float(b["sigma"]) - float(ci)) < 1e-6 * max(1.0, ci):
        return None, (f"{b.get('source')}: sigma {float(b['sigma']):.4g} is the fit's own Laplace sigma (a fit output, "
                      "not a measured or class band)")
    if b.get("lo") is not None and b.get("hi") is not None:
        return (float(b["lo"]), float(b["hi"])), str(b.get("source"))
    return None, f"{b.get('source')}: no two-sided measured or class band"


def judge(cid: str, row: dict, fit_rec: dict, coord: dict, spec_bounds: Tuple[float, float]) -> Dict[str, Any]:
    iv = coord["intervals"]["scaled"]
    lo_b, hi_b = spec_bounds
    band, band_note = s4_band(row)
    inside = None if band is None else bool(band[0] <= iv["lower"] and iv["upper"] <= band[1])
    pinned = bool(iv["closed_lower"] and iv["closed_upper"])
    region = [p for p in coord["low_side"] + coord["high_side"]
              if iv["lower"] - 1e-9 <= p["value"] <= iv["upper"] + 1e-9 or p["delta"] <= iv["threshold_delta"]]
    walked = coord["low_side"] + coord["high_side"]
    non_conv = [p["value"] for p in walked if not p["converged_profile"]]
    driver_flagged = [p["value"] for p in walked if not p["converged"]]
    failed = [p["value"] for p in walked if p["failed_rows"]]
    conv_inside = all(p["converged_profile"] and not p["failed_rows"] for p in region)
    refit_ok = bool(fit_rec["refit"]["converged_profile"] and not fit_rec["refit"]["failed_rows"])
    shape = ("pinned" if pinned else "flat" if not (iv["closed_lower"] or iv["closed_upper"]) else
             "one-sided (open " + ("below" if not iv["closed_lower"] else "above") + ", to the search bound)")
    conditional = GROUP_OF[cid] == "sulfur"
    shipped = row.get("value")
    shipped_inside = None if shipped is None else bool(iv["lower"] - 1e-9 <= float(shipped) <= iv["upper"] + 1e-9)
    meets = bool(pinned and conv_inside and refit_ok and inside is not False and shipped_inside is not False)
    width = None if band is None else band[1] - band[0]
    spans_band = None if width is None else bool(2 * iv["half_width"] >= width)
    if spans_band:
        meets = False
    return {"shape": shape, "pinned_both_sides_by_likelihood": pinned, "inside_s4_band": inside, "s4_band": band,
            "s4_band_note": band_note, "ci_spans_s4_band": spans_band, "refit_converged": refit_ok,
            "shipped_value_inside_interval": shipped_inside,
            "converged_inside_interval": conv_inside, "non_converged_points": non_conv,
            "driver_flagged_points_polish_raised_cost": [v for v in driver_flagged if v not in non_conv],
            "points_with_failed_pots": failed, "meets_profile_requirement": meets,
            "conditional": CONDITIONAL_B0 if conditional else None,
            "fit_rows_flagged_in_s1_record": list(row.get("fit_rows_withdrawn_or_corrected") or []),
            "fit_objective_reproduced": bool(fit_rec["frozen_objective_reproduced"]),
            "proposed_upgrade": bool(meets and not conditional and not row.get("fit_rows_withdrawn_or_corrected")
                                     and fit_rec["frozen_objective_reproduced"]),
            "search_bounds": [lo_b, hi_b]}


# =====================================================================================================================
# 4. Running one fit, and the artifact
# =====================================================================================================================

def _h0(row: dict, kind: str, lo: float, hi: float) -> float:
    ci = (row.get("identification") or {}).get("ci95_halfwidth")
    h = (float(ci) / 2.5) if ci else 0.1 * (hi - lo)
    return float(min(max(h, MIN_STEP[kind]), (hi - lo) / 4.0))


def selected(fit: str, only: Optional[Sequence[str]]) -> List[str]:
    ids = [cid for cid, (f, _) in COORDINATES.items() if f == fit]
    return [c for c in ids if not only or c in only]


def run_fit(fit: str, only: Optional[Sequence[str]], workers: int) -> Dict[str, Any]:
    rows = tier_rows()
    spec = SPECS[fit]()
    ids = selected(fit, only)
    idx = {cid: spec.keys.index(COORDINATES[cid][1]) for cid in ids}
    h0 = {idx[c]: _h0(rows[c], spec.kinds[idx[c]], float(spec.lower[idx[c]]), float(spec.upper[idx[c]])) for c in ids}
    rec = profile_fit(spec, [idx[c] for c in ids], h0, workers=workers)
    fit_block = {"fit": fit, **spec.meta, "keys": list(spec.keys), "kinds": list(spec.kinds),
                 "lower": [float(v) for v in spec.lower], "upper": [float(v) for v in spec.upper],
                 "n_rows": spec.n_rows, "n_free": len(spec.keys), "dof": spec.dof,
                 "driver": {"fit": "src.kinetic_core.fitting.fit", "settings": fitting.LEAST_SQUARES_SETTINGS,
                            "jacobian": "forward differences, absolute per-kind steps (fitting.STEPS / EXTRA_STEPS)",
                            "central_polish_calls": spec.central_polish, "budget_per_fit": spec.budget,
                            "workers": workers},
                 "refit": rec["refit"], "S_reported_at_freeze": REPORTED_S[fit](),
                 "frozen_objective_reproduced": bool(abs(rec["refit"]["S_at_frozen"] - REPORTED_S[fit]())
                                                     <= REPRODUCED_REL * max(1.0, REPORTED_S[fit]())),
                 "S_refit": rec["S_refit"], "S_min": rec["S_min"],
                 "lower_optimum_found_by_profile": rec["lower_optimum_found_by_profile"],
                 "chi2_red": rec["chi2_red"], "thresholds_delta": rec["thresholds_delta"], "seconds": rec["seconds"]}
    coords = {}
    for cid in ids:
        key = COORDINATES[cid][1]
        c = rec["coordinates"][key]
        verdict = judge(cid, rows[cid], fit_block, c, (float(spec.lower[idx[cid]]), float(spec.upper[idx[cid]])))
        coords[cid] = {"id": cid, "lane": rows[cid]["lane"], "group": GROUP_OF[cid], "fit": fit, "fit_coordinate": key,
                       "shipped_value": rows[cid]["value"], "s1_fit_named": rows[cid]["identification"]["fit"],
                       "laplace_ci95_halfwidth": rows[cid]["identification"].get("ci95_halfwidth"),
                       "x_frozen": float(spec.x_frozen[idx[cid]]), **{k: c[k] for k in ("x_hat", "h0")},
                       "grid": sorted([p["value"] for p in c["low_side"] + c["high_side"]] + [c["x_hat"]]),
                       "profile": sorted([{"value": c["x_hat"], "S": rec["S_refit"], "delta": rec["S_refit"] - rec["S_min"],
                                           "converged": fit_block["refit"]["converged_profile"], "driver_converged": fit_block["refit"]["converged"], "status": fit_block["refit"]["status"],
                                           "evaluations": fit_block["refit"]["evaluations"],
                                           "failed_rows": fit_block["refit"]["failed_rows"], "others": None}]
                                         + [{"value": p["value"], "S": p["S"], "delta": p["delta"], "converged": p["converged_profile"], "driver_converged": p["converged"],
                                             "status": p["status"], "evaluations": p["evaluations"],
                                             "failed_rows": p["failed_rows"], "polish": p.get("polish"),
                                             "retried": p.get("retried"),
                                             "others": {k: round(v, 6) for j, (k, v) in enumerate(zip(spec.keys, p["x"]))
                                                        if j != idx[cid]}}
                                            for p in c["low_side"] + c["high_side"]], key=lambda q: q["value"]),
                       "intervals": c["intervals"], "verdict": verdict, "spent_data": SPENT_DATA}
    return {"fit": fit_block, "coordinates": coords}


def time_point(fit: str, only: Optional[Sequence[str]], workers: int) -> Dict[str, Any]:
    """Seconds for ONE profile point (the others re-optimised from the frozen point at psi = frozen + h0), and the
    projected wall-clock for the fit's coordinates at ~``points`` points each over ``workers`` chains."""
    rows = tier_rows()
    spec = SPECS[fit]()
    ids = selected(fit, only)
    cid = ids[0]
    i = spec.keys.index(COORDINATES[cid][1])
    h = _h0(rows[cid], spec.kinds[i], float(spec.lower[i]), float(spec.upper[i]))
    t0 = time.time()
    r = spec.residuals(spec.x_frozen)
    t_eval = time.time() - t0
    rec = refit(spec, {i: float(spec.x_frozen[i] + h)}, spec.x_frozen)
    points = 2 * 9
    serial = rec["seconds"] * (points * len(ids) + 1)
    chains = 2 * len(ids)
    wall = rec["seconds"] * (1 + 9 * math.ceil(chains / max(1, workers)))
    return {"fit": fit, "coordinate": cid, "seconds_one_evaluation": round(t_eval, 2),
            "S_at_frozen": float(np.asarray(r) @ np.asarray(r)), "seconds_one_profile_point": rec["seconds"],
            "evaluations_one_point": rec["evaluations"], "converged": rec["converged"], "status": rec["status"],
            "polish_flagged": rec["polish_flagged"], "S_point": rec["S"],
            "projected_hours_serial": round(serial / 3600.0, 2),
            "projected_hours_wall": round(wall / 3600.0, 2), "assumed_points_per_side": 9, "workers": workers}


def _f(x, nd=3):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        return f"{x:.{nd}f}" if (abs(x) >= 10 ** (-nd) or x == 0) else f"{x:.2e}"
    return str(x)


def summary(payload: Dict[str, Any]) -> Dict[str, Any]:
    coords = payload["coordinates"]
    by_group: Dict[str, Dict[str, Any]] = {}
    for g in GROUP_ORDER:
        ids = [c for c in COORDINATES if GROUP_OF[c] == g]
        done = [c for c in ids if c in coords]
        by_group[g] = {"coordinates": len(ids), "profiled": len(done),
                       "meet_requirement": [c for c in done if coords[c]["verdict"]["meets_profile_requirement"]],
                       "flat_or_one_sided": [c for c in done if not coords[c]["verdict"]["pinned_both_sides_by_likelihood"]],
                       "not_profiled": [c for c in ids if c not in coords]}
    by_lane: Dict[str, Dict[str, int]] = {}
    for c, rec in coords.items():
        d = by_lane.setdefault(rec["lane"], {"profiled": 0, "meet_requirement": 0})
        d["profiled"] += 1
        d["meet_requirement"] += int(rec["verdict"]["meets_profile_requirement"])
    return {"by_work_group": by_group, "by_s1_lane": by_lane,
            "proposed_upgrades": [c for c in COORDINATES if c in coords and coords[c]["verdict"]["proposed_upgrade"]],
            "conditional_no_upgrade": [c for c in COORDINATES if c in coords and coords[c]["verdict"]["conditional"]]}


def markdown(p: Dict[str, Any]) -> str:
    s = p["summary"]
    out = ["# S1 profile likelihoods", "",
           f"_Generated by `{p['provenance']['generated_by']}`; the rule: plan sec. 7 S1 and `src/evidence_tiers.py`. "
           "Identification only, never credit (the fits' own rows). Changes no tier: an upgrade needs this artifact and a "
           "second reader, entered in `docs/reference/evidence_tier_changes.yml`._", "",
           "## Method", "", p["method"]["what"], "", f"Calibration: {p['method']['calibration']}", "",
           f"Judgement: {p['method']['judgement']}", "", "## Summary", "",
           "| work group | coordinates | profiled | meet the requirement | flat or one-sided | not profiled |",
           "|---|---:|---:|---|---|---|"]
    for g in GROUP_ORDER:
        d = s["by_work_group"][g]
        out.append(f"| {g} | {d['coordinates']} | {d['profiled']} | {', '.join(d['meet_requirement']) or '—'} | "
                   f"{', '.join(d['flat_or_one_sided']) or '—'} | {', '.join(d['not_profiled']) or '—'} |")
    out += ["", "Proposed upgrades (each needs a second reader before it enters the register): "
            + (", ".join(s["proposed_upgrades"]) or "none"), ""]
    if s["conditional_no_upgrade"]:
        out += [f"Conditional, no upgrade proposed ({CONDITIONAL_B0}): " + ", ".join(s["conditional_no_upgrade"]), ""]
    for f in p.get("findings", []):
        out.append(f"- **Finding.** {f}")
    for name, nr in (p.get("not_run") or {}).items():
        out.append(f"- **Not run: {name}.** {nr['reason']} Measured: {nr['measured']} Projection: {nr['projection']}")
    out += ["", "## Fits", "", "| fit | rows | free | dof | S in its report | S at frozen now | reproduced | S refit | "
            "S min | chi2_red | refit converged | Delta thresholds (raw / scaled / F) | wall s |",
            "|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|---:|"]
    for name, fb in p["fits"].items():
        t = fb["thresholds_delta"]
        out.append(f"| {name} | {fb['n_rows']} | {fb['n_free']} | {fb['dof']} | {_f(fb['S_reported_at_freeze'])} | "
                   f"{_f(fb['refit']['S_at_frozen'])} | {_f(fb['frozen_objective_reproduced'])} | "
                   f"{_f(fb['S_refit'])} | {_f(fb['S_min'])} | {_f(fb['chi2_red'])} | {_f(fb['refit']['converged_profile'])} | "
                   f"{_f(t['raw'], 2)} / {_f(t['scaled'], 2)} / {_f(t['F'], 2)} | {_f(fb['seconds'], 0)} |")
    out += ["", "Held fixed by each fit (not profiled; for the second reader):", ""]
    for name, fb in p["fits"].items():
        out.append(f"- **{name}** ({fb['report']}): {fb['held_fixed']}."
                   + (f" {fb['attribution_note']}." if fb.get("attribution_note") else "")
                   + (f" {fb['rows_withdrawn_note']}." if fb.get("rows_withdrawn_note") else ""))
    out += ["", "## Coordinates", "",
            "| coordinate | fit | shipped | x_hat | Laplace ±95 | LR 95 % (scaled) | raw | F | shape | shipped inside | "
            "inside S4 band | converged | meets | note |", "|---|---|---:|---:|---:|---|---|---|---|---|---|---|---|---|"]
    for cid in COORDINATES:
        c = p["coordinates"].get(cid)
        if c is None:
            out.append(f"| {cid} | {COORDINATES[cid][0]} | | | | not profiled | | | | | | | | |")
            continue
        iv, v = c["intervals"], c["verdict"]

        def rng(x):
            if x is None:
                return "—"
            return (f"[{_f(x['lower'])}{'' if x['closed_lower'] else ' (open)'}, "
                    f"{_f(x['upper'])}{'' if x['closed_upper'] else ' (open)'}]")
        notes = []
        if v["conditional"]:
            notes.append("conditional on B0")
        if v.get("fit_rows_flagged_in_s1_record"):
            notes.append("fit rows flagged in the S1 record: " + "; ".join(v["fit_rows_flagged_in_s1_record"]))
        if not v.get("fit_objective_reproduced", True):
            notes.append("the fit's frozen objective is not reproduced by today's engine (re-derive first)")
        if v["inside_s4_band"] is None:
            notes.append("no S4 band: " + v["s4_band_note"])
        notes.append("PROPOSED UPGRADE" if v.get("proposed_upgrade") else "not proposed")
        note = " / ".join(notes)
        out.append(f"| {cid} | {c['fit']} | {_f(c['shipped_value'])} | {_f(c['x_hat'])} | "
                   f"{_f(c['laplace_ci95_halfwidth'])} | {rng(iv['scaled'])} | {rng(iv['raw'])} | {rng(iv['F'])} | "
                   f"{v['shape']} | {_f(v['shipped_value_inside_interval'])} | {_f(v['inside_s4_band'])} | {_f(v['converged_inside_interval'] and v['refit_converged'])} | "
                   f"**{_f(v['meets_profile_requirement'])}** | {note} |")
    out += ["", "## Profiles", "", "Delta = S_p(psi) - S_min (S = r . r on the fit's own rows); c = the others' fit "
            "converged; f = rows at the failed residual.", ""]
    for cid in COORDINATES:
        c = p["coordinates"].get(cid)
        if c is None:
            continue
        pts = "; ".join(f"{_f(q['value'])}: {_f(q['delta'], 2)}{'' if q['converged'] else ' (not c)'}"
                        f"{' f' + str(q['failed_rows']) if q['failed_rows'] else ''}" for q in c["profile"])
        out.append(f"- **{cid}** ({c['fit']} `{c['fit_coordinate']}`, search bounds "
                   f"[{_f(c['verdict']['search_bounds'][0])}, {_f(c['verdict']['search_bounds'][1])}]): {pts}")
    out.append("")
    return "\n".join(out)


#: fits whose profile was timed and not run (task brief: stop after the first three lanes when the sulfur lane would
#: take more than ~6 h wall). Run ``--fit B60 --workers 4`` to compute them; the record then replaces this entry.
NOT_RUN: Dict[str, Dict[str, str]] = {
    "B60": {"coordinates": "the 10 sulfur rate constants and b8.ph_drift.arp_secondary_ammonium_pKa",
            "reason": ("the sulfur-lane profiles (" + CONDITIONAL_B0 + ") would take far more than the ~6 h wall "
                       "budget, so the run stopped after the trunk, acrylamide and b18 lanes."),
            "measured": ("one profile point (22 coordinates re-optimised at one value of k_arp_tdp, B60's quick "
                         "setting, 2 pot workers) had not finished after 34 min of wall time and was stopped; B60's "
                         "own fit spent 7601 s on 246 evaluations (31 s each), so one forward Jacobian (23 columns) "
                         "costs about 12 min."),
            "projection": ("about 18 points per coordinate x 11 coordinates = ~200 points at >= 34 min each: >= 110 h "
                           "serial, >= 28 h wall on 4 chains, plus the full refit. Porting B60's objective to the JAX "
                           "fit path (the network carries sulfur since A2 + H3) is the route to a feasible run.")},
}

METHOD = {
    "what": ("For each coordinate, on the fit that froze its shipped value and that fit's own residual function, rows, "
             "sigmas, free set and search bounds: a full re-fit from the frozen point (the engine's fit driver, "
             "`src.kinetic_core.fitting.fit`, with a central-difference polish), then a walk outward on both sides "
             "(first step = the Laplace 95 % half-width / 2.5, growing) re-optimising every OTHER free coordinate of "
             "the fit at each grid value, warm-started from the neighbouring point, until the profile rises past the "
             "largest threshold (one point beyond) or reaches the search bound. S_min is the lowest S seen (the refit "
             "or any profile point). JAX was not used: none of these fits' objectives is a `JaxProblem`."),
    "calibration": ("likelihood ratio, Delta = S_p(psi) - S_min with S = r . r. raw: Delta <= chi2_1(0.95) = 3.841 "
                    "(Wilks, the rows' stated sigmas). scaled (the verdict): Delta / max(1, chi2_red) <= 3.841, the "
                    "residual variance inflated by the fit's reduced chi-square when above one, as each fit's own "
                    "Laplace covariance does. F: Delta / chi2_red <= F(0.95; 1, dof), reported when chi2_red > 1 as a "
                    "small-sample check. Edges by signed-root interpolation. Asymptotic (regular interior optimum, "
                    "Gaussian residuals); not checked by simulation."),
    "judgement": ("meets the profile requirement = the scaled interval closes on BOTH sides by the likelihood strictly "
                  "inside the search bounds (a search bound never closes it: S4), the refit and every profile point "
                  "inside the interval converged (the optimiser stopped on a tolerance and the central polish was not "
                  "exhausted; a polish that raised the cost carries the forward endpoint, as fitting.fit does, and is "
                  "listed; a non-converged point is retried once with twice the polish) with no failed pot, the "
                  "SHIPPED value inside the interval (where the "
                  "engine changed under the fit's rows since it froze, today's optimum can sit elsewhere: S1's "
                  "'re-derived for the coordinates whose rows a structural change moves'), and the interval inside a measured or class band "
                  "where one is on record (none of the 28 has one: their bands are the fit's own Laplace sigma, the "
                  "search bounds, a fixed value or nothing; containment is then undefined and the closure is the "
                  "likelihood's alone). Sulfur-lane and pH-drift coordinates: " + CONDITIONAL_B0 + ". A coordinate "
                  "that meets the requirement is PROPOSED as an upgrade only when it is not conditional, the S1 "
                  "record flags none of its fit's rows as withdrawn or corrected (a profile on contested rows pins the "
                  "coordinate to them, not to the chemistry), and today's engine reproduces the fit's frozen objective "
                  "(else the fit needs re-deriving before any tier moves)."),
}


def load_payload() -> Dict[str, Any]:
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text(encoding="utf-8"))
    return {"artifact": "s1_profiles", "fits": {}, "coordinates": {}}


def finalize(p: Dict[str, Any]) -> Dict[str, Any]:
    inputs = [TIERS_JSON] + sorted({ROOT / fb["report"].split(" ")[0] for fb in p["fits"].values()})
    p["provenance"] = provenance.provenance_block("s1_profiles", generated_by="scripts/generators/build_s1_profiles.py",
                                                  inputs=inputs)
    p["rule"] = {"plan_sec7_S1": RULE_S1, "evidence_tiers": RULE_CODE}
    p["method"] = METHOD
    p["caveats"] = {"spent_data": SPENT_DATA, "sulfur_and_ph_drift": CONDITIONAL_B0,
                    "held_fixed": "only the fit's own free coordinates are re-optimised; what each fit held fixed is "
                                  "listed per fit (held_fixed) for the second reader",
                    "tiers": "no tier is changed here; proposed upgrades need a second reader before they enter "
                             "docs/reference/evidence_tier_changes.yml"}
    findings = []
    if "B41" in p["fits"]:
        findings.append("The S1 record attributes the fed 3-deoxy triangle (k_tdg_ddg, k_ddg_tdg, k_ddg_dgal, "
                        "k_dgal_ddg, k_ddg_hmf) to B39 and reads B39's Laplace, but the shipped values are B41's "
                        "optimum (parameters_dicarbonyl.SHIPPING_FIT_REPORT); B41 is profiled here.")
    for name, fb in p["fits"].items():
        if not fb["frozen_objective_reproduced"]:
            findings.append(f"{name}: today's engine does NOT reproduce the fit's frozen objective (S = r . r at the "
                            f"frozen point {fb['refit']['S_at_frozen']:.3f} now, {fb['S_reported_at_freeze']:.3f} in "
                            "its report): the engine changed under its rows since it froze, so its coordinates need "
                            "re-deriving (S1) and the profile is of today's objective.")
        moved = [k for k, a, b in zip(fb["keys"], fb["refit"]["x_frozen"], fb["refit"]["x"]) if abs(a - b) > 1e-3]
        if fb["refit"]["S_at_frozen"] - fb["S_min"] > 0.05 * max(1.0, fb["S_min"]) or fb["lower_optimum_found_by_profile"]:
            findings.append(f"{name}: today's engine re-optimises the fit's objective below its frozen point (S at frozen "
                            f"{fb['refit']['S_at_frozen']:.3f}, S_min {fb['S_min']:.3f}; coordinates moved > 1e-3: "
                            f"{', '.join(moved) or 'none'}); the profile is around today's optimum.")
    p["findings"] = findings
    p["not_run"] = {k: v for k, v in NOT_RUN.items() if k not in p["fits"]}
    p["summary"] = summary(p)
    p["coordinates"] = {c: p["coordinates"][c] for c in COORDINATES if c in p["coordinates"]}
    return p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", choices=sorted(SPECS))
    ap.add_argument("--only", nargs="*", help="S1 row ids (default: every coordinate of the fit)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--time-point", action="store_true")
    ap.add_argument("--render", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        ok = OUT_JSON.exists() and OUT_MD.exists() and OUT_MD.read_text(encoding="utf-8") == markdown(load_payload())
        print("s1_profiles.md is current" if ok else "s1_profiles.md is STALE: run --render")
        return 0 if ok else 1
    if a.render:
        OUT_MD.write_text(markdown(load_payload()), encoding="utf-8")
        print(f"wrote {data_paths.rel(OUT_MD)}")
        return 0
    if not a.fit:
        ap.error("--fit, --render or --check")
    if a.time_point:
        print(json.dumps(time_point(a.fit, a.only, a.workers), indent=2))
        return 0
    res = run_fit(a.fit, a.only, a.workers)
    p = load_payload()
    p["fits"][a.fit] = res["fit"]
    p["coordinates"].update(res["coordinates"])
    p = finalize(p)
    OUT_JSON.write_text(json.dumps(p, indent=2, default=float) + "\n", encoding="utf-8")
    OUT_MD.write_text(markdown(p), encoding="utf-8")
    for cid, c in res["coordinates"].items():
        iv = c["intervals"]["scaled"]
        print(f"{cid}: [{iv['lower']:.4g}{'' if iv['closed_lower'] else ' open'}, {iv['upper']:.4g}"
              f"{'' if iv['closed_upper'] else ' open'}] {c['verdict']['shape']} meets={c['verdict']['meets_profile_requirement']}")
    print(f"wrote {data_paths.rel(OUT_JSON)} and {data_paths.rel(OUT_MD)} ({res['fit']['seconds']} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
