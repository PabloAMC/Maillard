#!/usr/bin/env python
"""
build_frozen_fit_reproducibility.py -- does today's engine still reproduce the frozen fits the shipped constants rest
on? (plan sec. 7 S1: "re-derived for the coordinates whose rows a structural change moves"; 2026-11-06)

WHY. The S1 profile run (draft PR #106, ``s1_profiles.md``) found that today's engine no longer reproduces some frozen
fits' objectives at their frozen points (B1: r . r 5026 against 3980 in its report; B18: 295 against 2.64). Since those
fits froze, the engine and the data changed under them: M0 data-correction batches, withdrawn rows (CF-131 and
others), the pH solve moved into the integrator (#93), the a_w window refusal (#77). A shipped constant whose fit no
longer reproduces rests on an objective that no longer exists.

WHAT. For every fit that froze at least one shipped coordinate (the S1 record's ``identification.fit`` of the
fitted-* rows, plus the analogue-tier rows whose value is a fit's output, B21), with that fit's OWN residual function
(imported from its frozen wave generator, never edited), rows, sigmas and bounds, at its FROZEN point:

1. ONE objective evaluation: S = r . r today against the S its report records; the row count then and now; each row's
   residual then (from the report, where it records rows) and now, and why it moved: the row was withdrawn or added,
   its target or sigma was corrected, or (same target and sigma) the model changed under it;
2. the forward-difference Jacobian at the frozen point (one evaluation per free coordinate, the fitting module's
   absolute per-kind steps), the gradient of S, the Gauss-Newton decrement (the decrease the local quadratic model
   predicts), and ONE bounded Gauss-Newton step evaluated (backtracked at most twice), which is a lower bound on how
   far below the frozen point today's optimum lies. No optimiser runs: this is not a refit;
3. the class:
   * ``REPRODUCES``: |S_today - S_recorded| <= tol x max(1, S_recorded), the tolerance the fit's own record states,
     else 1 %, and the same row count;
   * ``DRIFTED-BENIGN``: it does not reproduce, but the frozen point is still near-optimal on today's objective: the
     evaluated one-step decrease stays below the joint 95 % likelihood-ratio threshold chi2_p(0.95) x s2, AND every
     shipped coordinate the fit froze lies inside today's local 95 % region (its quadratic marginal Delta_j =
     delta_S' C_SS^-1 delta_S <= chi2_|S|(0.95) x s2, delta the Gauss-Newton step, C = (J'J)^+);
     and the local quadratic model those tests rest on holds at the evaluated step (trust ratio, evaluated over
     predicted decrease, >= 0.25);
   * ``DRIFTED-MATERIAL``: any of the three fails: the shipped constants it froze are not supported by any current
     objective (a failed trust ratio: near-optimality cannot be shown, which is no support).
   s2 = max(1, chi2_red AT THE FREEZE): the fit's own variance scaling when it froze. Today's chi2_red is not used,
   because a drift that inflates the misfit would otherwise widen its own tolerance.

WHAT IT IS NOT. Identification bookkeeping, never credit (the fits' own rows, item 2's spent data). It changes no
tier, register or constant; it reports what S1 / S4 imply for each drifted fit (the PR that adds it lists them).

Usage (each fit runs in its own process: several generators install module state at import):
    python scripts/generators/build_frozen_fit_reproducibility.py --fit B20 [--workers 3]
    python scripts/generators/build_frozen_fit_reproducibility.py --all [--workers 3]
    python scripts/generators/build_frozen_fit_reproducibility.py --time B60     # seconds for one evaluation
    python scripts/generators/build_frozen_fit_reproducibility.py --probe B18    # one evaluation + row sets, JSON (the drift gate)
    python scripts/generators/build_frozen_fit_reproducibility.py --render       # markdown from the json
    python scripts/generators/build_frozen_fit_reproducibility.py --check        # exit 1 if the markdown is stale
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import subprocess
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

from src import data_paths, evidence_tiers, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT_JSON = V / "frozen_fit_reproducibility.json"
OUT_MD = V / "frozen_fit_reproducibility.md"
TIERS_JSON = V / "incumbent_evidence_tiers.json"
GENERATED_BY = "scripts/generators/build_frozen_fit_reproducibility.py"

#: the default reproduction tolerance (relative, on S = r . r) when a fit's record states none
DEFAULT_REL_TOL = 0.01
#: forward-difference steps per coordinate kind: the fitting module's (B8's) absolute steps, plus a logit step for
#: B6's unconstrained simplex logits
STEPS: Dict[str, float] = {"Ea": 2.0, "log10k": 0.05, "acid_yield": 0.01, "pka": 0.1, "ph_slope": 0.05, "logit": 0.05}
CLASSES = ("REPRODUCES", "DRIFTED-BENIGN", "DRIFTED-MATERIAL")


# =====================================================================================================================
# 1. The audit (pure: a residual function, bounds and kinds; no chemistry). Unit-tested on a toy problem.
# =====================================================================================================================

def chi2_95(dof: int) -> float:
    from scipy.stats import chi2

    return float(chi2.ppf(0.95, max(1, int(dof))))


def forward_jacobian(residuals: Callable[[np.ndarray], np.ndarray], x: np.ndarray, r0: np.ndarray,
                     lower: np.ndarray, upper: np.ndarray, kinds: Sequence[str],
                     map_fn: Callable[[Callable, List[np.ndarray]], List[np.ndarray]] = None) -> Tuple[np.ndarray, List[float]]:
    """Forward differences with an absolute step per kind, stepping inward at an upper bound."""
    x = np.asarray(x, dtype=float)
    points, steps = [], []
    for j, kind in enumerate(kinds):
        h = STEPS[kind]
        if x[j] + h > upper[j]:
            h = -h
        p = x.copy()
        p[j] += h
        points.append(p)
        steps.append(h)
    mapper = map_fn or (lambda f, pts: [f(p) for p in pts])
    values = mapper(residuals, points)
    jac = np.zeros((len(r0), len(x)))
    for j, (h, r) in enumerate(zip(steps, values)):
        jac[:, j] = (np.asarray(r, dtype=float) - r0) / h
    return jac, steps


def gauss_newton_step(r0: np.ndarray, jac: np.ndarray, x: np.ndarray, lower: np.ndarray,
                      upper: np.ndarray) -> np.ndarray:
    """The bounded Gauss-Newton step: argmin ||r0 + J d|| subject to lower <= x + d <= upper."""
    from scipy.optimize import lsq_linear

    lo = np.asarray(lower, dtype=float) - x
    hi = np.asarray(upper, dtype=float) - x
    lo = np.minimum(lo, 0.0)
    hi = np.maximum(hi, 0.0)
    sol = lsq_linear(jac, -np.asarray(r0, dtype=float), bounds=(lo, hi), method="bvls", lsmr_tol="auto")
    return np.asarray(sol.x, dtype=float)


def near_optimality(residuals: Callable[[np.ndarray], np.ndarray], x0: np.ndarray, r0: np.ndarray,
                    jac: np.ndarray, lower: np.ndarray, upper: np.ndarray, s2: float,
                    coordinates: Mapping[str, Sequence[int]], backtracks: int = 2) -> Dict[str, Any]:
    """Gradient, Gauss-Newton decrement, one evaluated bounded step, and each shipped coordinate's quadratic
    marginal Delta at the frozen point (from today's local optimum estimate x0 + d)."""
    x0 = np.asarray(x0, dtype=float)
    r0 = np.asarray(r0, dtype=float)
    s0 = float(r0 @ r0)
    grad = 2.0 * jac.T @ r0
    jtj = jac.T @ jac
    cov = np.linalg.pinv(jtj, rcond=1e-12)
    # the unconstrained GN decrement r' J (J'J)^+ J' r: the decrease the quadratic model predicts with no bounds
    g_half = jac.T @ r0
    decrement = float(g_half @ cov @ g_half)
    d = gauss_newton_step(r0, jac, x0, lower, upper)
    predicted = float(s0 - np.sum((r0 + jac @ d) ** 2))
    tried = []
    best = {"S": s0, "scale": 0.0}
    scale = 1.0
    for _ in range(backtracks + 1):
        x1 = np.clip(x0 + scale * d, lower, upper)
        r1 = np.asarray(residuals(x1), dtype=float)
        s1 = float(r1 @ r1)
        tried.append({"scale": scale, "S": s1})
        if s1 < best["S"]:
            best = {"S": s1, "scale": scale}
        if s1 < s0:
            break
        scale *= 0.5
    actual = s0 - best["S"]
    marg = {}
    for cid, idx in coordinates.items():
        idx = list(idx)
        if not idx:
            continue
        dS = d[idx]
        c_ss = cov[np.ix_(idx, idx)]
        delta = float(dS @ np.linalg.pinv(c_ss, rcond=1e-12) @ dS)
        thr = chi2_95(len(idx)) * s2
        sd = [float(math.sqrt(max(cov[i, i], 0.0) * s2)) for i in idx]
        marg[cid] = {"indices": idx, "gn_displacement": [float(v) for v in dS],
                     "local_sd_scaled": sd, "delta_quadratic": delta, "threshold": thr,
                     "inside_local_95": bool(delta <= thr)}
    return {"S_frozen": s0, "gradient_norm": float(np.linalg.norm(grad)),
            "gradient": [float(v) for v in grad], "gn_decrement_unbounded": decrement,
            "gn_step": [float(v) for v in d], "gn_predicted_decrease_bounded": predicted,
            "step_evaluations": tried, "one_step_decrease": float(actual),
            "one_step_scale": best["scale"], "s2": float(s2), "coordinates": marg}


#: the trust ratio (evaluated / quadratic-predicted decrease at the evaluated step) below which the local quadratic
#: model is not evidence of near-optimality (the trust-region convention: a step is "acceptable" above 1/4)
MIN_TRUST_RATIO = 0.25


def trust_ratio(opt: Mapping[str, Any]) -> Optional[float]:
    """Evaluated over predicted decrease at the evaluated step. The quadratic model along the step d is
    q(t) = -(a t + b t^2), a = grad . d (= 2 r0' J d), b = ||J d||^2 = -q(1) - a; it is read off the stored record
    (gradient, step, predicted decrease at t = 1), so a re-render needs no evaluation."""
    g = np.asarray(opt["gradient"], dtype=float)
    d = np.asarray(opt["gn_step"], dtype=float)
    a = float(g @ d)
    b = -float(opt["gn_predicted_decrease_bounded"]) - a
    evals = opt.get("step_evaluations") or []
    if not evals:
        return None
    best = min(evals, key=lambda e: e["S"])
    t = float(best["scale"])
    predicted = -(a * t + b * t * t)
    actual = float(opt["S_frozen"]) - float(best["S"])
    if predicted <= 1e-12 * max(1.0, float(opt["S_frozen"])):
        return None  # nothing to predict: the frozen point is a stationary point of the model
    return actual / predicted


def classify(s_recorded: float, s_today: float, n_then: Optional[int], n_now: int, rel_tol: float,
             opt: Optional[Mapping[str, Any]], n_free: int) -> Dict[str, Any]:
    """REPRODUCES / DRIFTED-BENIGN / DRIFTED-MATERIAL (see the module docstring)."""
    diff = s_today - s_recorded
    reproduced = bool(abs(diff) <= rel_tol * max(1.0, abs(s_recorded)) and (n_then is None or n_then == n_now))
    joint_thr = None if opt is None else chi2_95(n_free) * float(opt["s2"])
    joint_ok = None if opt is None else bool(opt["one_step_decrease"] <= joint_thr)
    outside = [] if opt is None else [c for c, m in opt["coordinates"].items() if not m["inside_local_95"]]
    rho = None if opt is None else trust_ratio(opt)
    model_ok = rho is None or rho >= MIN_TRUST_RATIO
    reasons = []
    if reproduced:
        cls = "REPRODUCES"
    elif opt is None:
        cls = "DRIFTED-MATERIAL"
        reasons.append("no near-optimality record")
    else:
        if not joint_ok:
            reasons.append(f"one evaluated step lowers S by {opt['one_step_decrease']:.3g}, beyond the joint "
                           f"threshold {joint_thr:.3g}")
        if outside:
            reasons.append("shipped coordinates outside today's local 95 % region: " + ", ".join(outside))
        if not model_ok:
            reasons.append(f"the local quadratic model fails at the step (trust ratio {rho:.2f} < {MIN_TRUST_RATIO}), "
                           "so near-optimality cannot be shown")
        cls = "DRIFTED-BENIGN" if not reasons else "DRIFTED-MATERIAL"
    return {"class": cls, "reproduced": reproduced, "S_change": diff,
            "S_change_rel": diff / max(1.0, abs(s_recorded)), "rel_tol": rel_tol,
            "joint_threshold": joint_thr, "joint_near_optimal": joint_ok,
            "coordinates_outside_local_95": outside, "trust_ratio": rho, "local_model_valid": model_ok,
            "material_reasons": reasons}


def reclassify(rec: Dict[str, Any]) -> Dict[str, Any]:
    """The class from a stored fit record (so a rule change re-renders without evaluating anything)."""
    return classify(rec["S_recorded"], rec["S_today"], rec["n_rows_then"], rec["n_rows_now"], rec["rel_tol"],
                    rec["near_optimality"], rec["n_free"])


def row_changes(then: Mapping[str, Mapping[str, Any]], now: Mapping[str, Mapping[str, Any]],
                top: int = 8) -> Dict[str, Any]:
    """Per-row comparison. Each side maps row id -> {r, target?, sigma?}. Returns every row's record, sorted by
    |r_now^2 - r_then^2|, and the counts of withdrawn / added / corrected / model-moved rows."""
    out = []
    for rid in list(then) + [k for k in now if k not in then]:
        a, b = then.get(rid), now.get(rid)
        rec: Dict[str, Any] = {"id": rid}
        if a is not None and b is None:
            rec.update({"kind": "withdrawn", "r_then": a.get("r"), "r_now": None,
                        "delta_r2": -float(a.get("r") or 0.0) ** 2})
        elif a is None and b is not None:
            rec.update({"kind": "added", "r_then": None, "r_now": b["r"], "delta_r2": float(b["r"]) ** 2})
        else:
            changed = []
            for key in ("target", "sigma"):
                ta, tb = a.get(key), b.get(key)
                if ta is not None and tb is not None and not math.isclose(float(ta), float(tb), rel_tol=1e-6,
                                                                          abs_tol=1e-12):
                    changed.append(f"{key} {float(ta):.4g} -> {float(tb):.4g}")
            ra = a.get("r")
            rec.update({"r_then": ra, "r_now": b["r"],
                        "delta_r2": (float(b["r"]) ** 2 - float(ra) ** 2) if ra is not None else None})
            if changed:
                rec["kind"] = "corrected"
                rec["changed"] = changed
            elif ra is None:
                rec["kind"] = "not recorded then"
            elif abs(float(b["r"]) - float(ra)) > 1e-3 * max(1.0, abs(float(ra))):
                rec["kind"] = "model change"
            else:
                rec["kind"] = "unchanged"
        out.append(rec)
    out.sort(key=lambda q: -abs(q["delta_r2"] or 0.0))
    counts: Dict[str, int] = {}
    for q in out:
        counts[q["kind"]] = counts.get(q["kind"], 0) + 1
    return {"counts": counts, "top": out[:top], "all": out}


# =====================================================================================================================
# 2. The fits. Each spec imports its frozen wave generator (nothing in a wave generator is edited) and states: the
#    frozen point, bounds, kinds, residual function, row ids and rows now; the recorded S and rows then; the
#    tolerance its record states; which shipped coordinates (S1 ids) map to which fit coordinates.
# =====================================================================================================================

@dataclass
class FitSpec:
    name: str
    keys: Tuple[str, ...]
    kinds: Tuple[str, ...]
    lower: np.ndarray
    upper: np.ndarray
    x_frozen: np.ndarray
    residuals: Callable[[np.ndarray], np.ndarray]
    rows_now: Dict[str, Dict[str, Any]]           # id -> {target, sigma} (r filled at evaluation)
    rows_then: Optional[Dict[str, Dict[str, Any]]]  # id -> {r, target, sigma}, or None if the report keeps no rows
    s_recorded: float
    n_rows_then: Optional[int]
    chi2_red_then: float
    rel_tol: float = DEFAULT_REL_TOL
    tol_source: str = "none stated in the record: 1 % (task default)"
    coordinate_keys: Dict[str, List[str]] = field(default_factory=dict)   # S1 id -> fit keys
    meta: Dict[str, Any] = field(default_factory=dict)
    groups: Optional[Callable[[np.ndarray], Dict[str, Any]]] = None  # optional group-level comparison (B1)
    #: attribution probes: each {label, cite, apply} where apply() installs one post-freeze change's reversal and
    #: returns its undo; S is evaluated once under each (diagnostic, never a fit)
    attribution: List[Dict[str, Any]] = field(default_factory=list)


def _load(stem: str, alias: Optional[str] = None):
    if str(GENERATORS) not in sys.path:
        sys.path.insert(0, str(GENERATORS))
    spec = importlib.util.spec_from_file_location(alias or f"ffr_{stem}", GENERATORS / f"{stem}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _report(stem: str) -> dict:
    return json.loads((V / f"{stem}.json").read_text(encoding="utf-8"))


def _s1_rows() -> List[dict]:
    return json.loads(TIERS_JSON.read_text(encoding="utf-8"))["rows"]


#: The S1 record's ``identification.fit`` label -> the fit audited here (kept in src/evidence_tiers.py, which notes a
#: drifted fit on its coordinates). The fed 3-deoxy triangle is attributed to B39 by the record, but its shipped values
#: are B41's optimum, so B41 is the fit they rest on; B39 is audited too, as the record's attribution.
S1_FIT_LABEL = evidence_tiers.FIT_OF_IDENTIFICATION_LABEL
FIT_ORDER = ("B1", "B3", "B6", "B7", "B18", "B20", "B21", "B39", "B41", "B2.3", "B60")


def shipped_coordinates() -> Dict[str, List[dict]]:
    """Fit -> the S1 rows whose shipped value that fit froze."""
    out: Dict[str, List[dict]] = {}
    for r in _s1_rows():
        ident = r.get("identification") or {}
        fit = S1_FIT_LABEL.get(ident.get("fit"))
        if fit is None:
            continue
        out.setdefault(fit, []).append(r)
    return out


def _rows_dict(ids: Sequence[str], targets: Sequence[Any], sigmas: Sequence[Any]) -> Dict[str, Dict[str, Any]]:
    return {i: {"target": t, "sigma": s} for i, t, s in zip(ids, targets, sigmas)}


# --- attribution probes: reverse one post-freeze change at a time and evaluate S once ---------------------------------

#: The trunk reaction groups added after the first trunk fit, with the commit that added each (the history of
#: src/kinetic_core/network.py).
TRUNK_GROUPS: Tuple[Tuple[str, str], ...] = (
    ("DICARBONYL_REACTIONS", "B13, d6e89620 (2026-09-07): glucosone, glyoxal and diacetyl join the trunk"),
    ("PYRAZINE_REACTIONS", "B18, 34da5b51 (2026-09-08): the aminoketone pyrazine steps"),
    ("GLYCATION_REACTIONS", "B20, efb76fb6 (2026-09-09): the glycation arm"),
    ("AQUEOUS_GLYOXAL_REACTIONS", "B21, 5c5f2857 (2026-09-09): the aqueous glucosone route to glyoxal"),
    ("METHIONINE_REACTIONS", "B22, 489f9e5f (2026-09-09): the methionine chain"),
    ("PROLINE_REACTIONS", "B24, 62fe149a (2026-09-09): 2-acetyl-1-pyrroline from proline"),
    ("FED_3DEOXY_REACTIONS", "B39 / B41, fb0dbd30 (2026-09-11): the fed 3-deoxyglucosone triangle"),
)


def _zeroed(parameters: Mapping[str, Any], keys) -> Dict[str, Any]:
    from dataclasses import replace

    out = dict(parameters)
    for k in keys:
        if k in out and hasattr(out[k], "k_ref"):
            out[k] = replace(out[k], k_ref=0.0)
    return out


def _group_parameter_keys(group: str, exclude_reactions=()) -> List[str]:
    from src.kinetic_core import network

    skip = set(exclude_reactions)
    return sorted({r.parameter_key for r in getattr(network, group) if r.key not in skip})


def _knockout_via_module_function(m, attr: str, keys: Sequence[str]) -> Callable[[], Callable[[], None]]:
    """Wrap a generator's module-level parameter builder so that ``keys`` are zeroed (a reaction knock-out)."""
    def apply():
        original = getattr(m, attr)
        setattr(m, attr, lambda *a, **k: _zeroed(original(*a, **k), keys))
        return lambda: setattr(m, attr, original)
    return apply


def _knockout_in_dict(d: Dict[str, Any], keys: Sequence[str]) -> Callable[[], Callable[[], None]]:
    def apply():
        saved = {k: d[k] for k in keys if k in d}
        d.update(_zeroed(d, keys))
        return lambda: d.update(saved)
    return apply


def _config_probe(updates: Mapping[str, Any]) -> Callable[[], Callable[[], None]]:
    def apply():
        from src.kinetic_core import config

        snapshot = config._snapshot()
        config.update(dict(updates))
        return lambda: config._restore(snapshot)
    return apply


def _cross_pair_factor(value: float) -> Callable[[], Callable[[], None]]:
    """CF-077 (19c66e99, 2026-11-01) reversed: the AKG + AKM cross pair back to statistical factor ``value``."""
    def apply():
        from src.kinetic_core import network

        rx = [r for r in network.TRUNK_REACTIONS if r.key == "r_ak_mpz"]
        saved = [r.statistical_factor for r in rx]
        for r in rx:
            object.__setattr__(r, "statistical_factor", float(value))

        def undo():
            for r, v in zip(rx, saved):
                object.__setattr__(r, "statistical_factor", v)
        return undo
    return apply


def _all_of(*applies: Callable[[], Callable[[], None]]) -> Callable[[], Callable[[], None]]:
    def apply():
        undos = [a() for a in applies]
        return lambda: [u() for u in reversed(undos)]
    return apply


def run_attribution(spec: "FitSpec", x: np.ndarray, s_today: float) -> List[Dict[str, Any]]:
    out = []
    for probe in spec.attribution:
        undo = probe["apply"]()
        try:
            r = np.asarray(spec.residuals(x), dtype=float)
            s = float(r @ r)
        finally:
            undo()
        out.append({"label": probe["label"], "cite": probe["cite"], "S": s, "S_change_vs_today": s - s_today,
                    "recovers_recorded": bool(abs(s - spec.s_recorded) <= spec.rel_tol * max(1.0, spec.s_recorded))})
    return out


def spec_b1() -> FitSpec:
    m = _load("generate_kinetic_core_b1_fit")
    from src.kinetic_core import FITTED_KEYS

    obs = m._load_fit_rows()
    sig = m.estimate_log_sigmas(obs)
    obj = m.Objective(obs, sig, free_melanoidin_sink=False)
    lower, upper = obj.bounds()
    rep = _report("kinetic_core_b1_fit_report")
    a = rep["variant_A"]
    fp = a["fitted_parameters"]
    x = np.array([fp[k]["log10_k_ref_100C"] for k in FITTED_KEYS] + [fp[k]["ea_kj_mol"] for k in FITTED_KEYS])
    keys = tuple(f"log10_{k}" for k in FITTED_KEYS) + tuple(f"ea_{k}" for k in FITTED_KEYS)
    ids = [f"{o.species}@{o.t_c:g}C/{o.t_min:g}min#{i}" for i, o in enumerate(obs)]
    rows_now = _rows_dict(ids, [o.value for o in obs], [sig[o.species] for o in obs])
    then_sig = rep.get("log_sigmas_from_replicate_scatter") or {}

    def groups(r: np.ndarray) -> Dict[str, Any]:
        pred = obj.predict(x)
        now: Dict[Tuple[str, float], List[int]] = {}
        for i, o in enumerate(obs):
            now.setdefault((o.species, o.t_c), []).append(i)
        then = {(g["species"], float(g["T_C"])): g for g in a["per_response_fit"]}
        recs = []
        for key in sorted(set(now) | set(then)):
            idx = now.get(key, [])
            signed = [math.log10((pred[i] + obs[i].floor) / (obs[i].value + obs[i].floor)) for i in idx]
            g = then.get(key, {})
            recs.append({"group": f"{key[0]} @ {key[1]:g} C", "n_then": g.get("n"), "n_now": len(idx),
                         "median_bias_dex_then": g.get("median_signed_log10_bias"),
                         "median_bias_dex_now": float(np.median(signed)) if signed else None,
                         "r2_now": float(np.sum(np.asarray(r)[idx] ** 2)) if idx else 0.0})
        for q in recs:
            a_, b_ = q["median_bias_dex_then"], q["median_bias_dex_now"]
            q["bias_shift_dex"] = None if a_ is None or b_ is None else b_ - a_
        recs.sort(key=lambda q: -abs(q["bias_shift_dex"] or 0.0))
        sig_rows = [{"species": s, "sigma_then": then_sig.get(s), "sigma_now": sig.get(s)} for s in sorted(sig)]
        return {"groups": recs, "sigmas": sig_rows}

    coords = {}
    for k in FITTED_KEYS:
        coords[f"{k}.log10_k_ref"] = [f"log10_{k}"]
        coords[f"{k}.ea_kj_mol"] = [f"ea_{k}"]
    from src.kinetic_core import network

    recorded = {r["key"] for r in rep["network"]["reactions"]}
    probes, all_keys = [], []
    late_base = [r.key for r in network.REACTIONS if r.key not in recorded]
    if late_base:
        keys_ = sorted({r.parameter_key for r in network.REACTIONS if r.key in late_base})
        all_keys += keys_
        probes.append({"label": f"knock out the base-table steps added after the freeze ({len(late_base)} steps)",
                       "cite": "network.REACTIONS against the B1 report's network.reactions (B7 furanic channels, "
                               "92b746d6, 2026-08-29, and later)",
                       "apply": _knockout_via_module_function(m, "operative_parameters", keys_)})
    for group, cite in TRUNK_GROUPS:
        keys_ = _group_parameter_keys(group, exclude_reactions=recorded)
        if not keys_:
            continue
        all_keys += keys_
        probes.append({"label": f"knock out {group} ({len(keys_)} constants)", "cite": cite,
                       "apply": _knockout_via_module_function(m, "operative_parameters", keys_)})
    probes.append({"label": "knock out every trunk step added after B1 froze", "cite": "all of the above",
                   "apply": _knockout_via_module_function(m, "operative_parameters", sorted(set(all_keys)))})
    return FitSpec("B1", keys, ("log10k",) * len(FITTED_KEYS) + ("Ea",) * len(FITTED_KEYS), lower, upper, x,
                   obj.residuals, rows_now, None, float(a["sum_of_squared_residuals"]), int(a["n_observations"]),
                   float(a["reduced_chi_square"]), coordinate_keys=coords, groups=groups, attribution=probes,
                   meta={"report": "results/validation/kinetic_core_b1_fit_report.json (variant_A)",
                         "generator": "scripts/generators/generate_kinetic_core_b1_fit.py (Objective, variant A)",
                         "objective": "Martins 2005's nine responses at 80/100/120 C, log-space residuals with additive "
                                      "floors over the per-species replicate scatter (re-estimated from today's data)",
                         "rows_then_note": "the report keeps per-response summaries, not rows: the comparison is per "
                                           "(species, temperature) group and per-species sigma"})


def spec_b3() -> FitSpec:
    m = _load("generate_kinetic_core_b3_fit")
    from src.kinetic_core.parameters_acrylamide import FITTED_ACRYLAMIDE_BOUNDS_LOG10K, FITTED_ACRYLAMIDE_EA_BOUNDS

    m._progress = lambda r: None
    rep = _report("kinetic_core_b3_fit_report")
    fr = rep["frozen_parameters"]
    keys = tuple(m.PARAM_ORDER) + tuple(m.EA_ORDER)
    x = np.array([fr["log10_k_ref_at_160C"][k] for k in m.PARAM_ORDER] + [fr["fitted_Ea_kJ_mol"][k] for k in m.EA_ORDER])
    lower = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][0] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[0]] * m.N_EA)
    upper = np.array([FITTED_ACRYLAMIDE_BOUNDS_LOG10K[k][1] for k in m.PARAM_ORDER] + [FITTED_ACRYLAMIDE_EA_BOUNDS[1]] * m.N_EA)
    rows_now = {r["id"]: {"target": r["target"], "sigma": r["sigma"]} for r in m.ACTIVE_FIT_ROWS}
    rows_then = {r["id"]: {"r": r["weighted_residual"], "target": r["observed"], "sigma": r["sigma"]}
                 for r in rep["rows"]}
    coords = {f"{k}.log10_k_ref": [k] for k in m.PARAM_ORDER}
    coords.update({f"k_{k[3:]}.ea_kj_mol": [k] for k in m.EA_ORDER})
    o = rep["objective"]
    return FitSpec("B3", keys, ("log10k",) * m.N_K + ("Ea",) * m.N_EA, lower, upper, x,
                   lambda v: m.residuals(np.asarray(v, dtype=float), True), rows_now, rows_then,
                   2.0 * float(o["final_cost"]), int(o["n_rows"]), float(o["reduced_chi_square"]),
                   coordinate_keys=coords,
                   meta={"report": "results/validation/kinetic_core_b3_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b3_fit.py residuals (quick=True, the "
                                      "setting the frozen fit ran at)",
                         "objective": "B3's declared fit rows (Claeys 2005, De Vleeschouwer 2009, Knol 2005, the SPI "
                                      "extrusion benchmark): log10 rate / ppb rows, linear Ea rows, a one-sided ceiling",
                         "cost_convention": "the report's final_cost is 0.5 r . r"})


def _b6_vector(fr: Mapping[str, Any]) -> np.ndarray:
    from src.kinetic_core.lipid import POSITION_PRODUCTS

    bm = fr["branch_model"]["simplexes"]
    v: List[float] = []
    for cell in ("13_ct", "13_tt", "9_ct", "9_tt"):
        prods = POSITION_PRODUCTS[cell.split("_")[0]]
        s = [float(bm[cell][p]) for p in prods]
        v += [math.log(s[1] / s[0]), math.log(s[2] / s[0])]

    def logit(p: float) -> float:
        p = min(max(p, 1e-300), 1.0 - 1e-16)
        return max(-60.0, min(60.0, math.log(p / (1.0 - p))))

    comp = fr["system_compositions"]
    mixed = comp["mixed_ct_tt_9_13"]
    f13 = mixed["LOOH_13_ct"] + mixed["LOOH_13_tt"]
    g_tt = mixed["LOOH_13_tt"] / f13 if f13 > 0 else 0.0
    v += [logit(f13), logit(g_tt), logit(comp["pure_ct_13"]["LOOH_13_ct"]), logit(comp["tt_9_13"]["LOOH_13_tt"])]
    return np.array(v, dtype=float)


def spec_b6() -> FitSpec:
    from src.kinetic_core import lipid

    rep = _report("kinetic_core_b6_fit_report")
    x = _b6_vector(rep["frozen_parameters"])
    ids = [f"{s}:{p}" for s, obs in lipid.FRANKEL_ZERO_ADDITIVE.items() for p in lipid.FRANKEL_SLATE]
    rows_now = {f"{s}:{p}": {"target": float(obs[p]), "sigma": 1.0}
                for s, obs in lipid.FRANKEL_ZERO_ADDITIVE.items() for p in lipid.FRANKEL_SLATE}
    rows_then = {f"{t['system']}:{t['product']}": {"r": t["residual_pp"], "target": t["measured_relative_percent"],
                                                    "sigma": 1.0} for t in rep["fit_table"]}
    assert len(ids) == 18
    q = rep["fit_quality"]
    keys = tuple(f"logit_{c}_{i}" for c in ("13_ct", "13_tt", "9_ct", "9_tt") for i in (1, 2)) + tuple(
        f"logit_{k}" for k in ("mixed_f13", "mixed_g_tt", "pure_ct_13_f13", "tt_9_13_f13"))
    coords = {f"lipid.branch_simplex.{c}": [f"logit_{c}_1", f"logit_{c}_2"] for c in ("13_ct", "13_tt", "9_ct", "9_tt")}
    inf = np.full(12, np.inf)
    return FitSpec("B6", keys, ("logit",) * 12, -inf, inf, x,
                   lambda v: np.asarray(lipid.fit_residuals(v), dtype=float), rows_now, rows_then,
                   float(q["sum_squared_residuals"]), int(q["n_fit_values"]),
                   float(q["sum_squared_residuals"]) / float(q["degrees_of_freedom"]), coordinate_keys=coords,
                   meta={"report": "results/validation/kinetic_core_b6_fit_report.json",
                         "generator": "src/kinetic_core/lipid.py fit_residuals (B6's own residual function)",
                         "objective": "Frankel 1989's eighteen zero-additive relative-percent values, residuals in "
                                      "percentage points (no sigma: s2 is the fit's own chi2_red, 3 dof)",
                         "frozen_point_note": "the 12 logits are rebuilt from the report's simplex shares and "
                                              "compositions (exact up to the +/-60 logit clip)"})


def spec_b7() -> FitSpec:
    m = _load("generate_kinetic_core_b7_fit")
    from src.kinetic_core.furanic import blank1997_fit_cells

    rep = _report("kinetic_core_b7_fit_report")
    f = rep["fit"]
    cells = blank1997_fit_cells("fit")
    sigma = math.log10(1.10)
    ids = [f"{c.sugar}/glycine HDMF" for c in cells]

    def residuals(v: np.ndarray) -> np.ndarray:
        try:
            pred = m._blank_prediction_ug_per_mmol(10.0 ** float(v[0]))
        except SystemExit as exc:  # the predictor refuses the Blank condition (out of envelope)
            residuals.refusal = str(exc)
            return np.full(len(cells), 25.0)
        if pred <= 0.0:
            return np.full(len(cells), 25.0)
        return np.array([(math.log10(pred) - math.log10(c.ug_per_mmol_sugar)) / sigma for c in cells])

    residuals.refusal = None
    rows_now = {i: {"target": c.ug_per_mmol_sugar, "sigma": sigma} for i, c in zip(ids, cells)}
    rows_then = {k: {"r": v / float(f["sigma_log10"]), "target": f["targets_ug_per_mmol"][k], "sigma": f["sigma_log10"]}
                 for k, v in f["residual_log10_fold"].items()}
    return FitSpec("B7", ("log10_k_dpo_af",), ("log10k",), np.array([-np.inf]), np.array([np.inf]),
                   np.array([float(f["log10_value"])]), residuals, rows_now, rows_then, float(f["objective"]),
                   len(f["residual_log10_fold"]), float(f["objective"]) / max(1, len(f["residual_log10_fold"]) - 1),
                   coordinate_keys={"k_dpo_af.log10_k_ref": ["log10_k_dpo_af"]},
                   meta={"report": "results/validation/kinetic_core_b7_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b7_fit.py _blank_prediction_ug_per_mmol "
                                      "(Blank 1997's three HDMF cells, sigma log10(1.10))",
                         "objective": "Blank 1997's three pentose/glycine HDMF yields at 90 C, 1 h, pH 6: one "
                                      "prediction (one generic aldopentose) against three cells, log10 residuals",
                         "refusal": lambda: residuals.refusal})


def _recorded_bounds(rep: dict, keys: Sequence[str], lower, upper) -> Tuple[np.ndarray, np.ndarray]:
    """The bounds the report records (a {key: [lo, hi]} map), else the module's (which may have been re-centred on
    the shipped value since the freeze)."""
    b = rep.get("bounds")
    if isinstance(b, dict) and all(k in b for k in keys):
        return (np.array([float(b[k][0]) for k in keys]), np.array([float(b[k][1]) for k in keys]))
    return np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)


def _wave_rows_then(rep: dict, resid_key: str, sigma_key: str = "sigma_log") -> Dict[str, Dict[str, Any]]:
    return {r["id"]: {"r": float(r[resid_key]) / float(r[sigma_key]), "target": r["target"], "sigma": r[sigma_key]}
            for r in rep["rows"]}


def spec_b18() -> FitSpec:
    m = _load("generate_kinetic_core_b18_fit")
    rep = _report("kinetic_core_b18_fit_report")
    x = np.array([rep["frozen_parameters"]["pyrazine"][k] for k in m.KEYS])
    lower, upper = _recorded_bounds(rep, m.KEYS, m.LOWER, m.UPPER)
    o = rep["objective"]
    rows_now = {r["id"]: {"target": r["target"], "sigma": r["sigma_log"]} for r in m.ROWS}
    sig_then = {r["id"]: r for r in o["rows"]}
    rows_then = {rid: {"r": float(v) / float(sig_then[rid]["sigma_log"]), "target": sig_then[rid]["target"],
                       "sigma": sig_then[rid]["sigma_log"]} for rid, v in rep["residual_by_row_dex"].items()}
    coords = {"k_go_ak.log10_k_ref": ["log10_k_go_ak_100C"], "k_go_ak.ea_kj_mol": ["ea_go_ak_kj_mol"],
              "k_mgo_ak.log10_k_ref": ["log10_k_mgo_ak_100C"], "k_mgo_ak.ea_kj_mol": ["ea_mgo_ak_kj_mol"],
              "b18.ph_slope_above_7_decades_per_unit": ["ph_slope_above_7_decades_per_unit"],
              "b18.ph_slope_below_7_decades_per_unit": ["ph_slope_below_7_decades_per_unit"]}
    late = [(g, c) for g, c in TRUNK_GROUPS if g not in ("DICARBONYL_REACTIONS", "PYRAZINE_REACTIONS")]
    ko = {g: _group_parameter_keys(g) for g, _ in late}
    exit_off = _config_probe({"trunk.three_deoxy_exit_ph_term": False})
    cf077 = _cross_pair_factor(1.0)
    probes = [{"label": "the 3-DG exit's pH term off", "cite": "B41, fb0dbd30 (2026-09-11): default on since B41 ships",
               "apply": exit_off},
              {"label": "the AKG + AKM cross pair back to statistical factor 1", "cite": "CF-077, 19c66e99 (2026-11-01)",
               "apply": cf077}]
    probes += [{"label": f"knock out {g} ({len(ko[g])} constants)", "cite": c,
                "apply": _knockout_in_dict(m.BASE_PARAMETERS, ko[g])} for g, c in late]
    probes.append({"label": "all of the above together", "cite": "every post-freeze trunk change probed here",
                   "apply": _all_of(exit_off, cf077, _knockout_in_dict(m.BASE_PARAMETERS,
                                                                       sorted({k for v in ko.values() for k in v})))})
    return FitSpec("B18", tuple(m.KEYS), ("log10k", "Ea", "log10k", "Ea", "ph_slope", "ph_slope"), lower, upper, x,
                   m.residuals, rows_now, rows_then, float(o["final_cost"]), int(o["n_rows"]),
                   float(o["final_cost"]) / max(1, int(o["n_rows"]) - int(o["n_free_parameters"])),
                   coordinate_keys=coords, attribution=probes,
                   meta={"report": "results/validation/kinetic_core_b18_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b18_fit.py residuals",
                         "objective": "Zhou 2024's six mean formation rates and Leahy 1989's four within-study pH "
                                      "ratios, log10(pred/target)/sigma"})


def _simple_wave(name: str, stem: str, block: str, coord_fmt: Callable[[str], str], resid_key: str,
                 cost_of: Callable[[dict], float], objective: str, quick_arg: Optional[bool] = None) -> FitSpec:
    m = _load(stem)
    rep = _report(stem.replace("generate_", "") + "_report")
    x = np.array([rep["frozen_parameters"][block][k] for k in m.KEYS])
    lower, upper = _recorded_bounds(rep, m.KEYS, m.LOWER, m.UPPER)
    rows_now = {r["id"]: {"target": r["target"], "sigma": r["sigma_log"]} for r in m.ROWS}
    rows_then = _wave_rows_then(rep, resid_key)
    s_rec = cost_of(rep)
    n_then = len(rows_then)
    resid = (lambda v: m.residuals(np.asarray(v, dtype=float), quick_arg)) if quick_arg is not None else m.residuals
    coords: Dict[str, List[str]] = {}
    for k in m.KEYS:
        coords[coord_fmt(k)] = [k]
    return FitSpec(name, tuple(m.KEYS), ("log10k",) * len(m.KEYS), lower, upper, x, resid, rows_now, rows_then,
                   s_rec, n_then, s_rec / max(1, n_then - len(m.KEYS)), coordinate_keys=coords,
                   meta={"report": f"results/validation/{stem.replace('generate_', '')}_report.json",
                         "generator": f"scripts/generators/{stem}.py residuals"
                                      + ("" if quick_arg is None else f" (quick={quick_arg})"),
                         "objective": objective})


def _k_of(coord: str) -> str:
    return coord.replace("log10_", "").replace("_aqueous_100C", "").replace("_100C", "") + ".log10_k_ref"


def spec_b20() -> FitSpec:
    return _simple_wave("B20", "generate_kinetic_core_b20_fit", "glycation", _k_of, "residual_dex",
                        lambda rep: min(mm["cost"] for mm in rep["members"]),
                        "Nguyen 2016 Table 1 (M1)'s ten printed rate constants at 120 / 130 C, log10(pred/target)/"
                        "sigma; each row one constant through its declared barrier (no integration)")


def spec_b21() -> FitSpec:
    return _simple_wave("B21", "generate_kinetic_core_b21_fit", "aqueous_glyoxal", _k_of, "residual_dex",
                        lambda rep: float(rep["objective"]["final_cost"]),
                        "Hamzalioglu 2026's six printed first-order constants (110-140 C), log10 residuals; two "
                        "log10 constants at 100 C free, barriers declared")


def _b39_like(name: str) -> FitSpec:
    if name == "B41":
        from src.kinetic_core import config
        from src.kinetic_core import trunk_conditions as TC

        config.update({"trunk.three_deoxy_exit_ph_term": True,
                       "trunk.three_deoxy_exit_ph": {"k_tdg_fa": TC.THREE_DEOXY_EXIT_PH["k_tdg_fa"]}})
    rep_stem = "kinetic_core_b41_fit" if name == "B41" else "kinetic_core_b39_fit"
    m = _load("generate_kinetic_core_b39_fit")
    rep = _report(rep_stem + "_report")
    assert rep["quick"] is False
    x = np.array([rep["frozen_parameters"]["fed_3deoxy"][k] for k in m.KEYS])
    lower, upper = _recorded_bounds(rep, m.KEYS, m.LOWER, m.UPPER)
    rows_now = {r["id"]: {"target": r["target"], "sigma": r["sigma_log"]} for r in m.DECISIVE}
    rows_then = {r["id"]: {"r": float(r["after_dex"]) / float(r["sigma_log"]), "target": r["target"],
                           "sigma": r["sigma_log"]} for r in rep["rows"] if r.get("decisive", True)}
    s_rec = min(mm["cost"] for mm in rep["members"])
    coords = {_k_of(k): [k] for k in m.KEYS}
    meta = {"report": f"results/validation/{rep_stem}_report.json",
            "generator": "scripts/generators/generate_kinetic_core_b39_fit.py residuals (quick=False, as frozen)"
                         + (", run as generate_kinetic_core_b41_fit.py runs it (the formic-acid exit's pH term on)"
                            if name == "B41" else ""),
            "objective": "B39's twelve decisive rows (Mittelmaier 2011 fed maxima, times of maximum and shares; "
                         "Zhang 2021 within-study ratios), log10(pred/target)/sigma"}
    if name == "B39":
        meta["attribution_note"] = ("the S1 record names B39 for the fed 3-deoxy triangle, but the shipped values are "
                                    "B41's optimum (parameters_dicarbonyl.SHIPPING_FIT_REPORT): B39 froze no shipped "
                                    "coordinate and is audited as the record's attribution")
    probes = ([] if name == "B41" else
              [{"label": "the 3-DG exit's pH term off (B39's own physics)",
                "cite": "B41, fb0dbd30 (2026-09-11): trunk.three_deoxy_exit_ph_term defaults on since B41 ships",
                "apply": _config_probe({"trunk.three_deoxy_exit_ph_term": False})}])
    return FitSpec(name, tuple(m.KEYS), ("log10k",) * len(m.KEYS), lower, upper, x,
                   lambda v: m.residuals(np.asarray(v, dtype=float), False), rows_now, rows_then, float(s_rec),
                   len(rows_then), float(rep["laplace"].get("chi2_reduced", s_rec / max(1, len(rows_then) - len(m.KEYS)))),
                   coordinate_keys=coords if name == "B41" else {}, meta=meta, attribution=probes)


def spec_b39() -> FitSpec:
    return _b39_like("B39")


def spec_b41() -> FitSpec:
    return _b39_like("B41")


def _b23_then(rep: dict, floor: float, sigma_ph: Optional[float] = None) -> Dict[str, Dict[str, Any]]:
    out = {}
    for r in rep["rows"]:
        p, t, s = r["predicted"], float(r["observed"]), float(r["sigma_log"])
        if p is None or not math.isfinite(float(p)):
            v = 25.0
        elif r["kind"] == "ph_endpoint":
            v = (float(p) - t) / (sigma_ph or s)
        else:
            lr = math.log10((float(p) + floor) / (t + floor))
            v = (max(0.0, lr) if r["kind"] in ("ceiling", "peak_fraction_ceiling")
                 else min(0.0, lr) if r["kind"] == "floor" else lr) / s
        out[r["id"]] = {"r": max(-25.0, min(25.0, v)), "target": t, "sigma": s}
    return out


def _sulfur_kind(k: str) -> str:
    return ("Ea" if k.startswith("Ea_") else "acid_yield" if k == "ph_acid_yield_per_sink_event"
            else "pka" if k == "ph_arp_secondary_ammonium_pKa" else "log10k")


def _sulfur_coords(rows: Sequence[dict], keys: Sequence[str], decay_families: Mapping[str, Sequence[str]],
                   on_formation: Sequence[str]) -> Dict[str, List[str]]:
    family_of = {k: f for f, ks in decay_families.items() for k in ks}
    out: Dict[str, List[str]] = {}
    for r in rows:
        cid, const = r["id"], r["constant"]
        if cid == "b8.ph_drift.acid_yield_per_sink_event":
            fk = "ph_acid_yield_per_sink_event"
        elif cid == "b8.ph_drift.arp_secondary_ammonium_pKa":
            fk = "ph_arp_secondary_ammonium_pKa"
        elif cid.endswith(".log10_k_ref"):
            fk = const
        elif cid.endswith(".ea_kj_mol"):
            fam = family_of.get(const)
            fk = (f"Ea_decay_{fam}" if fam and const not in on_formation else "Ea_lumped_formation")
        else:
            continue
        if fk in keys:
            out[cid] = [fk]
    return out


def spec_b23() -> FitSpec:
    m = _load("generate_kinetic_core_b2_3_fit")
    m._progress = lambda r: None
    rep = _report("kinetic_core_b2_3_fit_report")
    fr = rep["frozen_parameters"]
    keys = tuple(m.PARAM_ORDER) + tuple(m.EXTRA_PARAM_NAMES)
    x = np.array([fr["log10_k_ref_at_145C"][k] for k in m.PARAM_ORDER] + [fr["lumped_formation_Ea_kJ_mol"]]
                 + [fr["decay_Ea_kJ_mol"][f] for f in m.DECAY_FAMILY_ORDER]
                 + [fr["ph_drift"]["acid_yield_per_sink_event"], fr["ph_drift"]["arp_secondary_ammonium_pKa"]])
    lower, upper = _b23_bounds(m)
    rows_now = {r["id"]: {"target": r["target"], "sigma": r["sigma_log"]} for r in m.ACTIVE_FIT_ROWS}
    rows_then = _b23_then(rep, m.FLOOR)
    o = rep["objective"]
    s1 = [r for r in shipped_coordinates().get("B2.3", [])]
    coords = _sulfur_coords(s1, keys, fr["decay_families"], fr["decay_keys_kept_on_formation_Ea"])
    from src.kinetic_core.ph_state import PhStateCorrections

    b51_off = _config_probe({"ph.b51": PhStateCorrections()})
    grid_old = _config_probe({"ph.node_grid": "uniform", "ph.fixed_point": "fixed_sweeps"})
    budget_old = _config_probe({"sulfur.flux_budget": "linear_rederived"})
    probes = [{"label": "B51's pH-state corrections off (B50 physics)",
               "cite": "B51 ships 2026-10-27 (ph.b51 defaults all on; kinetic_core_b51_prereg.md)", "apply": b51_off},
              {"label": "B60's pH-node grid and relaxed fixed point back to uniform / fixed sweeps",
               "cite": "B60 ships 2026-10-28 (ph.node_grid, ph.fixed_point: 'frozen waves set uniform')",
               "apply": grid_old},
              {"label": "B60's closed flux budget back to linear_rederived",
               "cite": "B60 (sulfur.flux_budget: 'frozen waves set linear_rederived')", "apply": budget_old},
              {"label": "all three: the freeze-era pH physics and numerics", "cite": "B51 + B60",
               "apply": _all_of(b51_off, grid_old, budget_old)}]
    return FitSpec("B2.3", keys, tuple(_sulfur_kind(k) for k in keys), lower, upper, np.clip(x, lower, upper),
                   lambda v: m.residuals(np.asarray(v, dtype=float), True), rows_now, rows_then,
                   2.0 * float(o["final_cost"]), int(o["n_rows"]), float(o["reduced_chi_square"]),
                   coordinate_keys=coords, attribution=probes,
                   meta={"report": "results/validation/kinetic_core_b2_3_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b2_3_fit.py residuals (quick=True, the "
                                      "setting the optimiser ran at), loaded alone (no later wave installed)",
                         "objective": "B2.2's 58 sulfur-lane fit rows after the charge-conservation fix, log10 "
                                      "residuals (pH endpoints in pH units, one-sided ceilings)",
                         "cost_convention": "the report's final_cost is 0.5 r . r",
                         "shipped_via": "B60 carries these 25 coordinates frozen at B2.3's values (B60 free_set."
                                        "frozen_keys); the shipped vector's other coordinates are B60's"})


def _b23_bounds(m) -> Tuple[np.ndarray, np.ndarray]:
    """B2.3's search bounds, exactly as its fit() builds them."""
    lo = ([m.FITTED_SULFUR_BOUNDS_LOG10K[k][0] for k in m.PARAM_ORDER]
          + [m.LUMPED_FORMATION_EA_BOUNDS[0], m.DECAY_EA_BOUNDS[0], m.DECAY_EA_BOUNDS[0], m.ACID_YIELD_BOUNDS[0],
             m.ARP_AMINE_PKA_BOUNDS[0]])
    hi = ([m.FITTED_SULFUR_BOUNDS_LOG10K[k][1] for k in m.PARAM_ORDER]
          + [m.LUMPED_FORMATION_EA_BOUNDS[1], m.DECAY_EA_BOUNDS[1], m.DECAY_EA_BOUNDS[1], m.ACID_YIELD_BOUNDS[1],
             m.ARP_AMINE_PKA_BOUNDS[1]])
    return np.array(lo, dtype=float), np.array(hi, dtype=float)


def spec_b60() -> FitSpec:
    m = _load("generate_kinetic_core_b60_fit")
    from generate_kinetic_core_b8_laplace import frozen_vector  # noqa: E402  (on sys.path via the B60 module)

    rep = _report("kinetic_core_b60_fit_report")
    lower_full, upper_full = m.full_bounds()
    x_full = np.clip(frozen_vector(rep), lower_full, upper_full)
    idx = list(m.FREE_INDEX)
    keys = tuple(m.ALL_KEYS[i] for i in idx)

    def residuals(v):
        xf = x_full.copy()
        xf[idx] = v
        return m.residual_vector(xf, True)

    best = min(rep["members"], key=lambda q: q["cost"])
    ids = [r["id"] for r in m.B23.ACTIVE_FIT_ROWS]
    rows_now = {r["id"]: {"target": r["target"], "sigma": r.get("sigma_log")} for r in m.B23.ACTIVE_FIT_ROWS}
    rows_then = {rid: {"r": float(best["residual_by_row"][rid])} for rid in best["row_ids"]}
    fr = rep["frozen_parameters"]
    s1 = shipped_coordinates().get("B60", [])
    coords = _sulfur_coords(s1, keys, fr["decay_families"], fr["decay_keys_kept_on_formation_Ea"])
    o = rep["objective"]
    s_rec = 2.0 * float(o["final_cost"])
    assert len(ids) == len(rows_now)
    return FitSpec("B60", keys, tuple(_sulfur_kind(k) for k in keys), lower_full[idx], upper_full[idx], x_full[idx],
                   residuals, rows_now, rows_then, s_rec, int(o["n_rows"]),
                   s_rec / max(1, int(o["n_rows"]) - int(o["n_free_parameters"])), coordinate_keys=coords,
                   rel_tol=DEFAULT_REL_TOL,
                   meta={"report": "results/validation/kinetic_core_b60_fit_report.json",
                         "generator": "scripts/generators/generate_kinetic_core_b60_fit.py residual_vector (quick, "
                                      "B60's verified quick setting)",
                         "objective": "B51's objective on B47's 52 corrected rows, B60's pH-node grid",
                         "cost_convention": "the report's final_cost is 0.5 r . r",
                         "rows_then_note": "the report's best member keeps residuals only (no targets or sigmas)",
                         "conditional": "the sulfur lane runs on the sulfur balance B0, falsified under S3 (PR #94, "
                                        "M1b part (i), amendment A2)"})


SPECS: Dict[str, Callable[[], FitSpec]] = {"B1": spec_b1, "B3": spec_b3, "B6": spec_b6, "B7": spec_b7, "B18": spec_b18,
                                           "B20": spec_b20, "B21": spec_b21, "B39": spec_b39, "B41": spec_b41,
                                           "B2.3": spec_b23, "B60": spec_b60}


# =====================================================================================================================
# 3. Running one fit
# =====================================================================================================================

def _map_parallel(fn, points: List[np.ndarray]) -> List[np.ndarray]:
    from src.kinetic_core import parallel

    return parallel.map_ordered(fn, points)


def audit_fit(spec: FitSpec, shipped_ids: Sequence[str], workers: int = 1) -> Dict[str, Any]:
    t0 = time.time()
    x0 = np.asarray(spec.x_frozen, dtype=float)
    r0 = np.asarray(spec.residuals(x0), dtype=float)
    t_eval = time.time() - t0
    s_today = float(r0 @ r0)
    ids = list(spec.rows_now)
    assert len(ids) == len(r0), (spec.name, len(ids), len(r0))
    now = {rid: {**spec.rows_now[rid], "r": float(v)} for rid, v in zip(ids, r0)}
    changes = row_changes(spec.rows_then, now) if spec.rows_then is not None else None
    groups = spec.groups(r0) if spec.groups is not None else None
    os.environ["MAILLARD_WORKERS"] = str(max(1, int(workers)))
    t1 = time.time()
    jac, steps = forward_jacobian(spec.residuals, x0, r0, spec.lower, spec.upper, spec.kinds,
                                  _map_parallel if workers > 1 else None)
    os.environ["MAILLARD_WORKERS"] = "1"
    t_jac = time.time() - t1
    s2 = max(1.0, float(spec.chi2_red_then))
    coord_idx = {cid: [spec.keys.index(k) for k in spec.coordinate_keys.get(cid, []) if k in spec.keys]
                 for cid in shipped_ids}
    t2 = time.time()
    opt = near_optimality(spec.residuals, x0, r0, jac, spec.lower, spec.upper, s2,
                          {c: i for c, i in coord_idx.items() if i})
    t_step = time.time() - t2
    verdict = classify(spec.s_recorded, s_today, spec.n_rows_then, len(r0), spec.rel_tol, opt, len(spec.keys))
    t3 = time.time()
    attribution = run_attribution(spec, x0, s_today) if not verdict["reproduced"] else []
    t_attr = time.time() - t3
    meta = {k: (v() if callable(v) else v) for k, v in spec.meta.items()}
    unmapped = [c for c, i in coord_idx.items() if not i]
    return {"fit": spec.name, **meta, "keys": list(spec.keys), "kinds": list(spec.kinds),
            "x_frozen": [float(v) for v in x0], "lower": [float(v) for v in spec.lower],
            "upper": [float(v) for v in spec.upper], "n_free": len(spec.keys),
            "S_recorded": spec.s_recorded, "S_today": s_today, "n_rows_then": spec.n_rows_then, "n_rows_now": len(r0),
            "chi2_red_then": spec.chi2_red_then, "chi2_red_today": s_today / max(1, len(r0) - len(spec.keys)),
            "rel_tol": spec.rel_tol, "tol_source": spec.tol_source,
            "shipped_coordinates": list(shipped_ids), "shipped_coordinates_unmapped": unmapped,
            "row_changes": changes, "groups": groups, "attribution": attribution, "jacobian_steps": steps,
            "near_optimality": opt, "verdict": verdict,
            "seconds": {"one_evaluation": round(t_eval, 2), "jacobian": round(t_jac, 1),
                        "step": round(t_step, 1), "attribution": round(t_attr, 1),
                        "total": round(time.time() - t0, 1)},
            "workers": workers}


def run_one(name: str, workers: int) -> Dict[str, Any]:
    shipped = shipped_coordinates()
    spec = SPECS[name]()
    ids = [r["id"] for r in shipped.get(name, [])]
    return audit_fit(spec, ids, workers=workers)


def probe(name: str) -> Dict[str, Any]:
    """ONE objective evaluation at the frozen point, with the row set then and now: what
    ``scripts/ci/frozen_fit_drift_gate.py`` judges (no Jacobian, no step, no attribution)."""
    t0 = time.time()
    spec = SPECS[name]()
    t1 = time.time()
    r = np.asarray(spec.residuals(np.asarray(spec.x_frozen, dtype=float)), dtype=float)
    t2 = time.time()
    return {"fit": name, "S_recorded": spec.s_recorded, "S_today": float(r @ r), "rel_tol": spec.rel_tol,
            "tol_source": spec.tol_source, "n_free": len(spec.keys), "n_rows_then": spec.n_rows_then,
            "n_rows_now": len(r), "rows_then": None if spec.rows_then is None else sorted(spec.rows_then),
            "rows_now": sorted(spec.rows_now),
            "seconds": {"setup": round(t1 - t0, 2), "one_evaluation": round(t2 - t1, 2)}}


#: the line prefix ``--probe`` prints its JSON under (wave generators may print at import)
PROBE_PREFIX = "FROZEN_FIT_PROBE "


def time_one(name: str) -> Dict[str, Any]:
    p = probe(name)
    return {"fit": name, "seconds_one_evaluation": p["seconds"]["one_evaluation"], "S_today": p["S_today"],
            "S_recorded": p["S_recorded"], "n_free": p["n_free"], "n_rows": p["n_rows_now"]}


# =====================================================================================================================
# 4. The artifact
# =====================================================================================================================

#: Why rows moved, by fit: the citations a reader needs (hand-kept; each entry names the evidence it rests on).
WHY: Dict[str, str] = {
    "B1": ("a model change, no data change: the 354 rows, their targets and the per-species sigmas are as frozen. The "
           "drift is the trunk steps added after B1 froze (2026-08-28), which B1's objective now integrates: knocking "
           "them all out recovers the recorded S exactly, and B21's aqueous glucosone step alone (r_ama_g, AMA -> "
           "glucosone + Gly, 5c5f2857, 2026-09-09) carries ~96 % of it: it drains DFG in Martins' glucose-glycine pots, "
           "shifting every response's median bias by -0.04 to -0.09 dex. B21's constant is Hamzalioglu 2026's "
           "lactulosyl-lysine-in-milk value carried to fructosyl-glycine in water (analogue tier), and B1 was never "
           "re-derived under it. The 11 base-table steps added since (the B7 furanic channels, 92b746d6, and later) carry "
           "123 more and the fed 3-deoxy triangle (B39 / B41) offsets 14 (the probes are not additive)."),
    "B18": ("a model change on the four Leahy 1989 pH-ratio rows (the six Zhou 2024 rows are unchanged): B21's aqueous "
            "glucosone step (5c5f2857, 2026-09-09, the day after B18 froze) moves them most (knocking it out alone "
            "takes S 295 -> 27), CF-077's cross-pair factor (19c66e99) the rest (-70 alone; CF-077 itself records "
            "2.64 -> 15.4 at the fit commit); with every probed change reversed S is 5.6, so ~3 of the 292 are not "
            "attributed by these probes. CF-077 already flags the frozen pH slopes conditional; this record adds that "
            "B21 moved them more."),
    "B39": ("a model change by design: the 3-DG exit's pH term (B41, fb0dbd30) now defaults on, and B39 is B41's "
            "term-off ablation; switching it off recovers B39's recorded S exactly. B39 froze no shipped coordinate."),
    "B3": ("+0.11 %: the Claeys rate and Ea rows move by <= 0.04 sigma (model change). The SPI extrusion row withdrawn "
           "under CF-131 (branch m0-data-corrections-batch1, not yet on main) is still in the generator's rows, so the "
           "row count is 30 then and now; when CF-131 lands, B3's objective loses a row with r = -3.98 (r^2 15.8, 15 % "
           "of S) and B3 must be re-classified."),
    "B2.3": ("a model change, no data change (56 of 58 rows moved, none withdrawn or re-targeted): B2.3's generator "
             "now integrates under the pH physics and numerics shipped after it froze, B51's pH-state corrections "
             "(2026-10-27) and B60's log-and-linear node grid with a relaxed fixed point (2026-10-28), which it never "
             "sets back; at B2.3's quick 9-node / 2-sweep setting the relaxed fixed point does not converge (residual "
             "0.05 pH, flagged). Reversing all three probes takes S 112.9 -> 45.1, not to 16.35: the rest is the "
             "measured sulfur barriers and channels shipped since (B8 onward), not probed here. The Kang MFT / FFT "
             "and Hofmann 1998 furfural + H2S rows move most (+2 to +5 sigma). The shipped vector is not B2.3's point: "
             "B60 re-fitted 23 of its 48 coordinates and carries the other 25 frozen at B2.3's values, so B60's "
             "objective (which reproduces) is the one these 25 now sit in, held, never free."),
    "B7": "+0.01 %: the one Blank 1997 prediction moved by 0.001 dex through today's sulfur-lane parameters (B60).",
}
#: What S1 / S4 imply for each class (the PR body lists the coordinates; nothing here moves a tier).
IMPLICATION = {
    "REPRODUCES": "none: the frozen objective still exists; the shipped coordinates keep their S1 tier.",
    "DRIFTED-BENIGN": "re-derive the fit's report on today's rows when next touched (S1 're-derived for the "
                      "coordinates whose rows a structural change moves'); the shipped point is still inside today's "
                      "local 95 % region, so no shipped output moves.",
    "DRIFTED-MATERIAL": "the shipped coordinates it froze rest on no current objective: under S1 they cannot be "
                        "fitted-identified until re-derived, so they stay fitted-unidentified and are carried under S4 "
                        "as bands (a measured or class band, else refused), or the fit is re-run through a case "
                        "(items 1 - 9) on today's rows.",
}


def _load_payload() -> Dict[str, Any]:
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text(encoding="utf-8"))
    return {"artifact": "frozen_fit_reproducibility", "fits": {}}


def _inputs(fits: Sequence[str]) -> List[Path]:
    stems = {"B1": "b1", "B3": "b3", "B6": "b6", "B7": "b7", "B18": "b18", "B20": "b20", "B21": "b21", "B39": "b39",
             "B41": "b41", "B2.3": "b2_3", "B60": "b60"}
    return [TIERS_JSON] + [V / f"kinetic_core_{stems[f]}_fit_report.json" for f in FIT_ORDER if f in fits]


def write_payload(payload: Dict[str, Any]) -> None:
    payload["provenance"] = provenance.provenance_block("frozen_fit_reproducibility", generated_by=GENERATED_BY,
                                                        inputs=_inputs(list(payload["fits"])))
    payload["method"] = METHOD
    for rec in payload["fits"].values():
        rec["verdict"] = reclassify(rec)
    payload["summary"] = summary(payload)
    OUT_JSON.write_text(json.dumps(payload, indent=1, default=str) + "\n", encoding="utf-8")
    OUT_MD.write_text(markdown(payload), encoding="utf-8")


METHOD = {
    "what": ("For every fit that froze a shipped coordinate (the S1 record's identification.fit; B41 for the fed "
             "3-deoxy triangle the record names B39; B21 for its two analogue-tier rows), the fit's own residual "
             "function imported from its frozen wave generator, at its frozen point: one objective evaluation "
             "(S = r . r, rows then and now), the forward-difference Jacobian (one evaluation per free coordinate, "
             "absolute per-kind steps: log10 k 0.05, Ea 2 kJ/mol, pH slope 0.05, pKa 0.1, acid yield 0.01, logit 0.05), "
             "and one bounded Gauss-Newton step evaluated (backtracked at most twice). No optimiser runs."),
    "classes": ("REPRODUCES: |S_today - S_recorded| <= tol x max(1, S_recorded) (the record's stated tolerance, else "
                "1 %) and the same row count. DRIFTED-BENIGN: not reproduced, but the evaluated one-step decrease is "
                "below chi2_p(0.95) x s2 (p = the fit's free coordinates) AND every shipped coordinate's quadratic "
                "marginal Delta_j = d_S' C_SS^-1 d_S <= chi2_|S|(0.95) x s2 (d the Gauss-Newton step, C = (J'J)^+), "
                "AND the quadratic model those tests use holds at the evaluated step (trust ratio = evaluated / "
                "predicted decrease >= 0.25). DRIFTED-MATERIAL: any of the three fails (a failed trust ratio means "
                "near-optimality cannot be shown, which S1 treats as no support). s2 = max(1, chi2_red at the freeze), never today's chi2_red (a drift "
                "must not widen its own tolerance)."),
    "limits": ("Local: the one-step decrease is a LOWER bound on S_frozen - S_min (a material verdict on it is "
               "certain), while a benign verdict rests on the quadratic model at the frozen point, which a strongly "
               "curved objective can understate; the S1 profile run's refits (PR #106) are the global check where "
               "they exist. Rows then come from each report's own row table; B1 keeps none (per-response summaries "
               "only) and B60 keeps residuals without targets. Identification bookkeeping, never credit; changes no "
               "tier, register or constant."),
}


def summary(p: Dict[str, Any]) -> Dict[str, Any]:
    by_class: Dict[str, List[str]] = {c: [] for c in CLASSES}
    material_coords: Dict[str, List[str]] = {}
    for name in FIT_ORDER:
        f = p["fits"].get(name)
        if f is None:
            continue
        by_class[f["verdict"]["class"]].append(name)
        if f["verdict"]["class"] == "DRIFTED-MATERIAL":
            material_coords[name] = list(f["shipped_coordinates"])
    return {"by_class": by_class, "not_run": [n for n in FIT_ORDER if n not in p["fits"]],
            "shipped_coordinates_on_material_fits": material_coords,
            "n_shipped_coordinates_on_material_fits": sum(len(v) for v in material_coords.values())}


def _f(x, nd=3):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, float):
        if x != 0 and (abs(x) < 10 ** (-nd) or abs(x) >= 1e6):
            return f"{x:.2e}"
        return f"{x:.{nd}f}"
    return str(x)


def markdown(p: Dict[str, Any]) -> str:
    s = p["summary"]
    out = ["# Frozen-fit reproducibility", "",
           f"_Generated by `{GENERATED_BY}`. Does today's engine, on today's data, still reproduce each frozen fit "
           "that a shipped constant rests on, and is the frozen point still near-optimal? One objective evaluation "
           "per fit plus one Jacobian and one Gauss-Newton step; no refits. Changes no tier, register or constant._",
           "", "## Method", "", p["method"]["what"], "", "**Classes.** " + p["method"]["classes"], "",
           "**Limits.** " + p["method"]["limits"], "", "## Summary", ""]
    for c in CLASSES:
        out.append(f"- **{c}**: {', '.join(s['by_class'][c]) or 'none'}")
    if s["not_run"]:
        out.append(f"- not run: {', '.join(s['not_run'])}")
    out += ["", f"Shipped coordinates resting on DRIFTED-MATERIAL fits: **{s['n_shipped_coordinates_on_material_fits']}**"
            + "".join(f"; {k}: {len(v)}" for k, v in s["shipped_coordinates_on_material_fits"].items()), "",
            "| fit | class | rows then / now | S recorded | S today | change | tol | one-step decrease | joint thr. "
            "| shipped coords | outside local 95 % | grad norm | eval s | wall s |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|"]
    for name in FIT_ORDER:
        f = p["fits"].get(name)
        if f is None:
            continue
        v, o = f["verdict"], f["near_optimality"]
        out.append(f"| {name} | **{v['class']}** | {_f(f['n_rows_then'])} / {f['n_rows_now']} | {_f(f['S_recorded'])} | "
                   f"{_f(f['S_today'])} | {100 * v['S_change_rel']:+.1f} % | {100 * v['rel_tol']:.0f} % | "
                   f"{_f(o['one_step_decrease'])} | {_f(v['joint_threshold'], 2)} | {len(f['shipped_coordinates'])} | "
                   f"{', '.join(v['coordinates_outside_local_95']) or '—'} | {_f(o['gradient_norm'], 2)} | "
                   f"{_f(f['seconds']['one_evaluation'], 1)} | {_f(f['seconds']['total'], 0)} |")
    out += ["", "What each class implies (plan sec. 7 S1 / S4; nothing here moves a tier):", ""]
    for c in CLASSES:
        out.append(f"- **{c}**: {IMPLICATION[c]}")
    out += ["", "## Fits", ""]
    for name in FIT_ORDER:
        f = p["fits"].get(name)
        if f is None:
            continue
        v, o = f["verdict"], f["near_optimality"]
        out += [f"### {name}: {v['class']}", "",
                f"- report: `{f.get('report')}`; residual function: {f.get('generator')}",
                f"- objective: {f.get('objective')}"]
        for k in ("cost_convention", "attribution_note", "shipped_via", "conditional", "rows_then_note",
                  "frozen_point_note"):
            if f.get(k):
                out.append(f"- {k.replace('_', ' ')}: {f[k]}")
        if f.get("refusal"):
            out.append(f"- the predictor refused at the frozen point: {f['refusal']}")
        out.append(f"- S = r . r: recorded {_f(f['S_recorded'])}, today {_f(f['S_today'])} "
                   f"({100 * v['S_change_rel']:+.2f} %; tolerance {100 * v['rel_tol']:.0f} %, {f['tol_source']}); rows "
                   f"{_f(f['n_rows_then'])} then, {f['n_rows_now']} now; chi2_red {_f(f['chi2_red_then'])} at the "
                   f"freeze, {_f(f['chi2_red_today'])} today (s2 = {_f(o['s2'])})")
        out.append(f"- near-optimality: |grad S| {_f(o['gradient_norm'])}; Gauss-Newton decrement {_f(o['gn_decrement_unbounded'])} "
                   f"(unbounded), {_f(o['gn_predicted_decrease_bounded'])} (bounded step, predicted); evaluated one-step "
                   f"decrease {_f(o['one_step_decrease'])} at scale {_f(o['one_step_scale'], 2)} against the joint threshold "
                   f"{_f(v['joint_threshold'], 2)} (chi2_{f['n_free']}(0.95) x s2); trust ratio {_f(v.get('trust_ratio'), 2)}")
        for reason in v.get("material_reasons") or []:
            out.append(f"- material: {reason}")
        if WHY.get(name) or f.get("why"):
            out.append(f"- why: {f.get('why') or WHY[name]}")
        if f.get("attribution"):
            out += ["", "Attribution probes (one post-freeze change reversed at a time, S evaluated once; diagnostic, "
                        "not a fit):", "", "| probe | source | S | vs today | recovers the record |",
                    "|---|---|---:|---:|---|"]
            for q in f["attribution"]:
                out.append(f"| {q['label']} | {q['cite']} | {_f(q['S'])} | {q['S_change_vs_today']:+.3f} | "
                           f"{_f(q['recovers_recorded'])} |")
            out.append("")
        if f["row_changes"] is not None:
            rc = f["row_changes"]
            out.append("- rows by kind: " + ", ".join(f"{k} {n}" for k, n in sorted(rc["counts"].items())))
            out += ["", "| row | kind | r then | r now | change in r^2 | note |", "|---|---|---:|---:|---:|---|"]
            for q in rc["top"]:
                out.append(f"| `{q['id']}` | {q['kind']} | {_f(q['r_then'])} | {_f(q['r_now'])} | {_f(q['delta_r2'])} | "
                           f"{'; '.join(q.get('changed', [])) or ''} |")
        if f.get("groups"):
            g = f["groups"]
            out += ["", "| group | n then / now | median bias then (dex) | now | shift | r^2 now |",
                    "|---|---|---:|---:|---:|---:|"]
            for q in g["groups"][:10]:
                out.append(f"| {q['group']} | {_f(q['n_then'])} / {q['n_now']} | {_f(q['median_bias_dex_then'])} | "
                           f"{_f(q['median_bias_dex_now'])} | {_f(q['bias_shift_dex'])} | {_f(q['r2_now'], 1)} |")
            moved = [q for q in g["sigmas"] if q["sigma_then"] is not None
                     and not math.isclose(float(q["sigma_then"]), float(q["sigma_now"]), rel_tol=1e-6)]
            out.append("")
            out.append("- per-species sigma (replicate scatter) then -> now: "
                       + (", ".join(f"{q['species']} {_f(q['sigma_then'])} -> {_f(q['sigma_now'])}" for q in moved)
                          or "unchanged"))
        if o["coordinates"]:
            out += ["", "| shipped coordinate | fit key | frozen | GN displacement | local sd (scaled) | Delta_j | "
                        "threshold | inside |", "|---|---|---:|---:|---:|---:|---:|---|"]
            for cid in f["shipped_coordinates"]:
                mrec = o["coordinates"].get(cid)
                if mrec is None:
                    continue
                keys = [f["keys"][i] for i in mrec["indices"]]
                xs = [f["x_frozen"][i] for i in mrec["indices"]]
                out.append(f"| {cid} | {', '.join(keys)} | {', '.join(_f(v_) for v_ in xs)} | "
                           f"{', '.join(_f(v_) for v_ in mrec['gn_displacement'])} | "
                           f"{', '.join(_f(v_) for v_ in mrec['local_sd_scaled'])} | {_f(mrec['delta_quadratic'], 2)} | "
                           f"{_f(mrec['threshold'], 2)} | {_f(mrec['inside_local_95'])} |")
        if f.get("shipped_coordinates_unmapped"):
            out.append(f"- shipped coordinates with no free key in this fit: {', '.join(f['shipped_coordinates_unmapped'])}")
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", choices=list(SPECS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--time", choices=list(SPECS))
    ap.add_argument("--probe", choices=list(SPECS), help="one evaluation at the frozen point, as JSON (the drift gate)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--render", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        p = _load_payload()
        ok = OUT_MD.exists() and OUT_MD.read_text(encoding="utf-8") == markdown(p)
        print("frozen_fit_reproducibility.md is " + ("current" if ok else "STALE: run --render"))
        return 0 if ok else 1
    if a.render:
        p = _load_payload()
        write_payload(p)
        print(f"wrote {data_paths.rel(OUT_MD)}")
        return 0
    if a.time:
        print(json.dumps(time_one(a.time), indent=1))
        return 0
    if a.probe:
        print(PROBE_PREFIX + json.dumps(probe(a.probe)), flush=True)
        return 0
    if a.all:
        for name in FIT_ORDER:
            print(f"[frozen-fit] {name} ...", flush=True)
            subprocess.run([sys.executable, __file__, "--fit", name, "--workers", str(a.workers)], check=True)
        return 0
    if a.fit:
        rec = run_one(a.fit, a.workers)
        if a.fit in WHY:
            rec["why"] = WHY[a.fit]
        p = _load_payload()
        p["fits"][a.fit] = rec
        write_payload(p)
        v = rec["verdict"]
        print(f"{a.fit}: {v['class']} | S {rec['S_recorded']:.3f} -> {rec['S_today']:.3f} | rows "
              f"{rec['n_rows_then']} -> {rec['n_rows_now']} | one-step -{rec['near_optimality']['one_step_decrease']:.3f} "
              f"| outside {v['coordinates_outside_local_95']} | {rec['seconds']}")
        return 0
    ap.error("one of --fit, --all, --time, --render, --check")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
