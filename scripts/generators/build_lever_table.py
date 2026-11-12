#!/usr/bin/env python
"""
The lever table from data (STATUS_AND_PLAN sec. 0, step 2c). Rules: docs/reference/lever_table/PREREG.md.

For each odorant x lever-class (x band) cell, the effect of one declared lever step, estimated from the resolved
double-keyed within-study contrasts in ``docs/reference/lever_table/contrasts/`` and nothing else. No engine output
enters. Writes ``results/validation/lever_table.{json,md}`` (all keyed data, the product) and
``lever_table_development_only.{json,md}`` (DEVELOPMENT laboratories of split v1 only: step 5's rival).

THE FIT. Every model is Gaussian with a few variance scales, so the posterior is computed exactly: the location
parameters (mu, or the per-level effects of a categorical cell, prior N(0, 3^2) dex) are integrated in closed form
and the scales (tau_lab, tau_study where PREREG sec. 5 allows it, sigma_resid) on a grid with half-normal priors.
Integrating the location with a wide prior and maximising over the scales is REML (the check PREREG asks for).
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, compound_keys, data_paths, provenance  # noqa: E402

OUT_ALL = data_paths.VALIDATION_DIR / "lever_table.json"
OUT_DEV = data_paths.VALIDATION_DIR / "lever_table_development_only.json"

PRIOR_SCALE = 0.5            # PREREG sec. 5: half-normal scale of every variance scale, dex
SENSITIVITY_SCALES = (0.25, 1.0)
LOCATION_SD = 3.0            # the wide normal prior on mu / level effects, dex
CV_DEFAULT, N_DEFAULT = 0.20, 2   # PREREG sec. 5 via step 0: 20 % per replicate, n = 2 where unstated
DIGITISED_SD = 0.05          # PREREG sec. 3
MIN_LABS, P_SAME_SIGN, LOLO_MIN_LABS, LOLO_FLAG = 3, 0.8, 4, 0.7
PRIOR_DOMINATED = 0.30
GRID = np.concatenate([[0.0], np.geomspace(0.01, 3.0, 34)])
#: log width of each grid cell. GRID is log-spaced, so a point's prior weight is the half-normal density times the
#: width of the cell it stands for; without it the grid applies ~ half-normal / s, which piles mass near zero and
#: makes sparse cells' new-laboratory intervals too narrow (fixed 2026-11-11; PREREG sec. 5 asks for half-normal)
GRID_LOG_WIDTH = np.log(np.gradient(GRID))

NUMERIC = {"pH": ("+1 pH unit", lambda x: x),
           "temperature": ("+10 C", lambda x: x / 10.0),
           "time": ("a doubling of time", math.log2),
           "precursor dose": ("a doubling of dose", math.log2),
           "water": ("+0.1 a_w", lambda x: x / 0.1),
           "buffer concentration": ("a doubling of concentration", math.log2)}
CATEGORICAL_REFERENCE = {"amino acid identity": "glycine",   # PREREG amendment A2
                         "sugar identity": "glucose", "sulfur source": "cysteine", "buffer": "phosphate",
                         "lipid": "none", "oxygen": "inert", "matrix": "buffer"}
#: PREREG sec. 4, "sugar identity: the named sugar": the levels that are sugars or sugar sources a formulator adds
#: (phosphorylated pentoses and IMP release ribose). A fed intermediate (furfural, a deoxyosone, norfuraneol) is not
#: a sugar level and is excluded from the cell, counted.
SUGAR_LEVELS = {"glucose", "fructose", "ribose", "xylose", "arabinose", "rhamnose", "galactose", "mannose", "maltose",
                "lactose", "sucrose", "ribose 5-phosphate", "glucose-6-phosphate", "glucose 6-phosphate",
                "inosine 5'-monophosphate", "imp", "guanosine 5'-monophosphate", "gmp", "none"}
PH_BANDS = ((-math.inf, 4.5), (4.5, 6.0), (6.0, 7.5), (7.5, math.inf))
TIME_BANDS = ((-math.inf, 30.0), (30.0, 120.0), (120.0, math.inf))       # for temperature contrasts
TEMP_BANDS = ((-math.inf, 100.0), (100.0, 130.0), (130.0, math.inf))     # for time contrasts


# ---------------------------------------------------------------------------------------------------- loading
def _num(x) -> Optional[float]:
    if x is None or isinstance(x, bool):
        return None
    try:
        return float(str(x).strip().split()[0].replace(",", ""))
    except (ValueError, IndexError):
        return None


def load_sources(directory: Path = data_paths.LEVER_TABLE_CONTRASTS_DIR) -> List[Dict[str, Any]]:
    """The resolved files only (``<source>.yml`` with ``reader: resolved``); reader files and the template are not
    data."""
    out = []
    for f in sorted(directory.glob("*.yml")):
        if f.name == "TEMPLATE.yml" or ".reader_" in f.name:
            continue
        doc = yaml.safe_load(f.read_text()) or {}
        if doc.get("reader") != "resolved":
            raise SystemExit(f"{f.name}: not a resolved file (reader: {doc.get('reader')!r})")
        for k in ("source", "contrasts"):
            if k not in doc:
                raise SystemExit(f"{f.name}: missing {k!r}")
        out.append(doc)
    return out


def lab_partition() -> Tuple[Dict[str, Tuple[str, str]], Dict[str, Any]]:
    """source key (lower case) -> (laboratory, partition) from split v1, plus new_labs.yml (TEST by PREREG sec. 2)."""
    split = yaml.safe_load(data_paths.CONTRAST_SPLIT.read_text())
    dev = set(split["split"]["development_labs"])
    out: Dict[str, Tuple[str, str]] = {}
    for c in split["contrasts"]:
        lab = c["lab"]
        out.setdefault(str(c["study"]).lower(), (lab, "development" if lab in dev else "test"))
    new = {}
    if data_paths.LEVER_TABLE_NEW_LABS.exists():
        new = yaml.safe_load(data_paths.LEVER_TABLE_NEW_LABS.read_text()) or {}
        for src, rec in new.items():
            lab = rec["lab"] if isinstance(rec, Mapping) else str(rec)
            out.setdefault(str(src).lower(), (lab, "development" if lab in dev else "test"))
    return out, {"split_sha256": split.get("sha256"), "new_labs": len(new)}


# ---------------------------------------------------------------------------------------------- observations
def arm_log_and_var(arm: Mapping[str, Any], digitised: bool) -> Tuple[Optional[float], Optional[float], str]:
    """(log10 value, its variance, censor state) of one arm; None where the arm enters no estimate."""
    censor = arm.get("censor") or "none"
    v = _num(arm.get("value"))
    if censor in ("nd_with_limit", "below_value"):
        v = _num(arm.get("limit"))
    if censor == "nd_no_limit" or v is None or v <= 0:
        return None, None, censor
    n = _num(arm.get("n")) or N_DEFAULT
    sd = _num(arm.get("sd"))
    if sd is not None and sd > 0 and censor == "none":
        var = (sd / v / math.log(10.0)) ** 2 / max(n, 1.0)
    else:
        var = math.log10(1.0 + CV_DEFAULT) ** 2 / N_DEFAULT
    if digitised:
        var += DIGITISED_SD ** 2
    return math.log10(v), var, censor


def _trace_fill(arms: List[Dict[str, Any]]) -> None:
    """PREREG sec. 3: 'tr' with no limit is censored at the lowest quantified value of the compound in the table
    (here: in the contrast, the narrowest table the keyed file carries)."""
    quantified = [_num(a.get("value")) for a in arms if (a.get("censor") or "none") == "none"]
    quantified = [q for q in quantified if q and q > 0]
    for a in arms:
        if a.get("censor") == "trace" and _num(a.get("limit")) is None:
            a["limit"] = min(quantified) if quantified else None
            a["censor"] = "below_value" if quantified else "nd_no_limit"


def _lever_x(lever_class: str, arm: Mapping[str, Any]) -> Optional[float]:
    raw = _num(arm.get("lever_value"))
    cond = arm.get("conditions") or {}
    if raw is None:
        raw = {"pH": _num(cond.get("ph")), "temperature": _num(cond.get("temp_C")),
               "time": _num(cond.get("time_min"))}.get(lever_class)
    if raw is None:
        return None
    try:
        return NUMERIC[lever_class][1](raw)
    except (ValueError, ZeroDivisionError):
        return None


def _band(bands, value) -> Optional[str]:
    if value is None:
        return None
    for lo, hi in bands:
        if lo <= value < hi:
            if lo == -math.inf:
                return f"<{hi:g}"
            if hi == math.inf:
                return f">{lo:g}"
            return f"{lo:g}-{hi:g}"
    return None


def _common(arms, key) -> Optional[float]:
    vals = {_num((a.get("conditions") or {}).get(key)) for a in arms}
    vals.discard(None)
    return vals.pop() if len(vals) == 1 else None


def numeric_observation(c: Mapping[str, Any], arms: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """One within-study slope per contrast (PREREG sec. 4): two arms give the normalised difference, a ladder the
    weighted least-squares slope over its longest monotone run."""
    lc = c["lever_class"]
    pts = []
    for a in arms:
        x = _lever_x(lc, a)
        ly, var, censor = arm_log_and_var(a, bool(c.get("digitised")))
        if x is not None and ly is not None:
            pts.append((x, ly, var, censor, _raw_lever(lc, a)))
    pts.sort()
    if len(pts) < 2 or len({p[0] for p in pts}) < 2:
        return None
    diffs = [np.sign(b[1] - a[1]) for a, b in zip(pts, pts[1:])]
    monotone = all(d >= 0 for d in diffs) or all(d <= 0 for d in diffs)
    run = pts
    if not monotone:          # the longest monotone run
        best = []
        for i in range(len(pts)):
            for j in range(i + 2, len(pts) + 1):
                seg = pts[i:j]
                ds = [np.sign(b[1] - a[1]) for a, b in zip(seg, seg[1:])]
                if (all(d >= 0 for d in ds) or all(d <= 0 for d in ds)) and len(seg) > len(best):
                    best = seg
        run = best
    x = np.array([p[0] for p in run])
    y = np.array([p[1] for p in run])
    w = 1.0 / np.array([p[2] for p in run])
    xb = float(np.sum(w * x) / np.sum(w))
    sxx = float(np.sum(w * (x - xb) ** 2))
    if sxx <= 0:
        return None
    slope = float(np.sum(w * (x - xb) * (y - np.sum(w * y) / np.sum(w))) / sxx)
    raw_x = [p[4] for p in run]
    mid = (min(raw_x) + max(raw_x)) / 2.0
    if lc == "pH":
        band = _band(PH_BANDS, mid)
        basis = {(a.get("conditions") or {}).get("ph_basis") or "unstated" for a in arms}
        sub = "held" if basis == {"held"} else ("final" if basis == {"final"} else "initial or label")
        band = f"pH {band}, {sub}"
    elif lc == "temperature":
        band = f"time {_band(TIME_BANDS, _common(arms, 'time_min')) or 'unstated'} min"
    elif lc == "time":
        band = f"T {_band(TEMP_BANDS, _common(arms, 'temp_C')) or 'unstated'} C"
    else:
        band = "all"
    return {"kind": "numeric", "y": slope, "v": 1.0 / sxx, "band": band, "non_monotone": not monotone,
            "range": [min(raw_x), max(raw_x)], "censored_arms": sum(1 for p in run if p[3] != "none"),
            "arms_used": len(run)}


def _raw_lever(lc: str, arm: Mapping[str, Any]) -> Optional[float]:
    raw = _num(arm.get("lever_value"))
    if raw is None:
        cond = arm.get("conditions") or {}
        raw = {"pH": _num(cond.get("ph")), "temperature": _num(cond.get("temp_C")),
               "time": _num(cond.get("time_min"))}.get(lc)
    return raw


def _inverse(lc, x):
    return {"pH": x, "temperature": x * 10.0, "time": 2.0 ** x, "precursor dose": 2.0 ** x,
            "water": x * 0.1, "buffer concentration": 2.0 ** x}[lc]


def _level(s: Any) -> str:
    t = str(s).strip().lower()
    for pre in ("l-", "d-"):
        if t.startswith(pre):
            t = t[len(pre):]
    return t


def categorical_observation(c: Mapping[str, Any], arms: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Pairs against the study's first-listed arm; the shared arm's variance is their covariance (PREREG sec. 3)."""
    vals = []
    for a in arms:
        ly, var, censor = arm_log_and_var(a, bool(c.get("digitised")))
        if ly is not None:
            vals.append((_level(a.get("lever_value") or a.get("label")), ly, var, censor))
    if len(vals) < 2:
        return None
    if c["lever_class"] == "sugar identity":
        dropped = [v for v in vals if v[0] not in SUGAR_LEVELS]
        vals = [v for v in vals if v[0] in SUGAR_LEVELS]
        if len(vals) < 2:
            return {"kind": "excluded", "reason": "levels are not sugars (fed intermediates)", "n": len(dropped)}
    ref = vals[0]
    pairs = [{"level": v[0], "against": ref[0], "y": v[1] - ref[1], "v": v[2] + ref[2]} for v in vals[1:]
             if v[0] != ref[0]]
    if not pairs:
        return None
    return {"kind": "categorical", "pairs": pairs, "shared_var": ref[2], "band": "all",
            "censored_arms": sum(1 for v in vals if v[3] != "none"), "arms_used": len(vals)}


def observations(sources, labs) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    obs, skipped = [], defaultdict(int)
    for doc in sources:
        src = str(doc["source"])
        lab, part = labs.get(src.lower(), (f"unmapped:{src}", "test"))
        if lab.startswith("unmapped:"):
            skipped["source without a laboratory (counted as its own TEST laboratory)"] += 1
        for c in doc["contrasts"] or []:
            lc = c.get("lever_class")
            if c.get("odorant_unlisted"):
                skipped["odorants outside the declared lists (PREREG amendment A3)"] += 1
                continue
            key = compound_keys.resolve(str(c.get("odorant", "")))
            if key is None:
                raise SystemExit(f"{src} {c.get('id')}: odorant {c.get('odorant')!r} does not resolve "
                                 "(src/compound_keys.py); fix the name or add the alias")
            q = str(c.get("quantity") or "amount")
            if q != "amount":
                skipped[{"share": "shares of a total (PREREG amendment A1)",
                         "fd_factor": "flavour-dilution factors, an odour potency not an amount (PREREG amendment A1)"}
                        .get(q, "rates, half-lives or storage losses (PREREG amendment A1)")] += 1
                continue
            if lc == "mixed":
                skipped["mixed-lever contrasts (no single-lever cell)"] += 1
                continue
            arms = [dict(a) for a in (c.get("arms") or [])]
            _trace_fill(arms)
            if lc in NUMERIC:
                o = numeric_observation(c, arms)
            elif lc in CATEGORICAL_REFERENCE:
                o = categorical_observation(c, arms)
            else:
                skipped[f"lever class {lc!r} not in PREREG sec. 4"] += 1
                continue
            if o is None:
                skipped["fewer than two usable arms"] += 1
                continue
            if o.get("kind") == "excluded":
                skipped[f"sugar-identity contrasts whose levels are not sugars (fed intermediates)"] += 1
                continue
            if c.get("system") == "fed_odorant":      # stability of an added odorant, not its formation
                o["band"] = f"{o['band']}, fed odorant (stability)"
            o.update({"source": src, "contrast": c.get("id"), "compound": key.display_name, "lever_class": lc,
                      "lab": lab, "partition": part, "focal": bool(c.get("focal")),
                      "peak_area": c.get("quant_type") != "absolute",
                      "matrix": sorted({str((a.get("conditions") or {}).get("matrix") or "unstated") for a in arms}),
                      "temps": [t for t in (_num((a.get("conditions") or {}).get("temp_C")) for a in arms) if t]})
            obs.append(o)
    return obs, dict(skipped)


# ------------------------------------------------------------------------------------------------- the model
def _halfnormal_logpdf(x, scale):
    return -0.5 * (x / scale) ** 2


def _grid_log_prior(scale: float) -> np.ndarray:
    """log prior weight of each GRID point under half-normal(scale): the density times the cell width."""
    return _halfnormal_logpdf(GRID, scale) + GRID_LOG_WIDTH


def fit_gaussian(y: np.ndarray, S: np.ndarray, X: np.ndarray, comps: Sequence[np.ndarray], scale: float
                 ) -> Dict[str, Any]:
    """y ~ N(X beta, S + sum_k s_k^2 K_k), beta ~ N(0, LOCATION_SD^2 I), s_k ~ half-normal(scale) on GRID.
    Returns the grid weights, each grid point's Gaussian posterior of beta, and the REML optimum."""
    n, p = X.shape
    P0 = np.eye(p) / LOCATION_SD ** 2
    pts, logw, means, covs, loglik = [], [], [], [], []
    log_prior = _grid_log_prior(scale)
    for idx in itertools.product(range(len(GRID)), repeat=len(comps)):
        s = tuple(float(GRID[i]) for i in idx)
        V = S.copy()
        for sk, K in zip(s, comps):
            V = V + sk ** 2 * K
        try:
            L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError:
            continue
        Vi = np.linalg.inv(V)
        A = X.T @ Vi @ X + P0
        Ai = np.linalg.inv(A)
        b = Ai @ (X.T @ Vi @ y)
        # log p(y | s) with beta integrated: N(0, V + X Sigma0 X^T)
        M = V + LOCATION_SD ** 2 * X @ X.T
        sign, logdet = np.linalg.slogdet(M)
        ll = -0.5 * (logdet + y @ np.linalg.solve(M, y) + n * math.log(2 * math.pi))
        lp = float(sum(log_prior[i] for i in idx))
        pts.append(s); logw.append(ll + lp); loglik.append(ll); means.append(b); covs.append(Ai)
    logw = np.array(logw)
    w = np.exp(logw - logw.max())
    w /= w.sum()
    reml = pts[int(np.argmax(loglik))]
    return {"points": np.array(pts), "w": w, "means": np.array(means), "covs": np.array(covs), "reml": reml}


def _mixture_quantiles(mu: np.ndarray, sd: np.ndarray, w: np.ndarray, qs=(0.05, 0.5, 0.95)) -> List[float]:
    lo, hi = float(np.min(mu - 8 * sd)), float(np.max(mu + 8 * sd))
    grid = np.linspace(lo, hi, 4001)
    cdf = (w[:, None] * stats.norm.cdf((grid[None, :] - mu[:, None]) / sd[:, None])).sum(0)
    return [float(np.interp(q, cdf, grid)) for q in qs]


def summarise(fit, j: int, extra_var_idx: Sequence[int]) -> Dict[str, Any]:
    """Coordinate j of beta: its 90 % interval, and the new-laboratory prediction (beta_j plus the scales in
    ``extra_var_idx``, i.e. a new laboratory's offset and a pot's residual)."""
    m = fit["means"][:, j]
    v = fit["covs"][:, j, j]
    q = _mixture_quantiles(m, np.sqrt(v), fit["w"])
    ev = v + sum(fit["points"][:, k] ** 2 for k in extra_var_idx)
    qn = _mixture_quantiles(m, np.sqrt(ev), fit["w"])
    sgn = 1.0 if q[1] >= 0 else -1.0
    p_same = float((fit["w"] * stats.norm.cdf(sgn * m / np.sqrt(ev))).sum())
    return {"estimate": q[1], "ci90": [q[0], q[2]], "new_lab_pi90": [qn[0], qn[2]], "p_same_sign_new_lab": p_same,
            "new_lab_pi90_width": qn[2] - qn[0]}


def _incidence(keys: Sequence[str]) -> np.ndarray:
    u = sorted(set(keys))
    Z = np.zeros((len(keys), len(u)))
    for i, k in enumerate(keys):
        Z[i, u.index(k)] = 1.0
    return Z @ Z.T


def numeric_fit(rows, scale):
    y = np.array([r["y"] for r in rows])
    S = np.diag([r["v"] for r in rows])
    X = np.ones((len(rows), 1))
    labs = [r["lab"] for r in rows]
    studies = [r["source"] for r in rows]
    per_lab = defaultdict(set)
    for r in rows:
        per_lab[r["lab"]].add(r["source"])
    with_study = sum(1 for s in per_lab.values() if len(s) >= 2) >= 2
    comps = [_incidence(labs)] + ([_incidence(studies)] if with_study else []) + [np.eye(len(rows))]
    fit = fit_gaussian(y, S, X, comps, scale)
    extra = list(range(len(comps)))      # a new lab's offset, a new study's, a pot's residual
    out = summarise(fit, 0, extra)
    out["tau_study_fitted"] = with_study
    out["reml_scales"] = [float(x) for x in fit["reml"]]
    return out


def categorical_fit(rows, reference, scale):
    """Network meta-analysis: fixed effect per level (reference = 0), one heterogeneity scale per pair with
    correlation 1/2 between pairs of one multi-arm study, measurement covariance from shared arms."""
    pairs = [(r, p) for r in rows for p in r["pairs"]]
    levels = sorted({p["level"] for _, p in pairs} | {p["against"] for _, p in pairs})
    # connectivity to the reference
    adj = defaultdict(set)
    for _, p in pairs:
        adj[p["level"]].add(p["against"]); adj[p["against"]].add(p["level"])
    comp, todo = {reference}, [reference]
    while todo:
        for nb in adj[todo.pop()]:
            if nb not in comp:
                comp.add(nb); todo.append(nb)
    if reference not in levels:
        # PREREG amendment A3: no study carries the declared reference; report against the most-connected level
        counts = {lv: sum(1 for _, p in pairs if lv in (p["level"], p["against"])) for lv in levels}
        fallback = sorted(levels, key=lambda lv: (-counts[lv], lv))[0]
        out = categorical_fit(rows, fallback, scale)
        out["reference_fallback"] = {"declared": reference, "used": fallback}
        return out
    pairs = [(r, p) for r, p in pairs if p["level"] in comp]
    cols = [lv for lv in levels if lv in comp and lv != reference]
    n = len(pairs)
    X = np.zeros((n, len(cols)))
    y = np.zeros(n)
    S = np.zeros((n, n))
    K = np.eye(n)
    for i, (r, p) in enumerate(pairs):
        y[i] = p["y"]
        if p["level"] in cols:
            X[i, cols.index(p["level"])] += 1.0
        if p["against"] in cols:
            X[i, cols.index(p["against"])] -= 1.0
        S[i, i] = p["v"]
        for k in range(i):
            if pairs[k][0] is r:
                S[i, k] = S[k, i] = r["shared_var"]
                K[i, k] = K[k, i] = 0.5
    fit = fit_gaussian(y, S, X, [K], scale)
    out = {"connected": True, "reference": reference, "levels": {}, "reml_scales": [float(x) for x in fit["reml"]]}
    for j, lv in enumerate(cols):
        s = summarise(fit, j, [0])
        s["labs"] = len({r["lab"] for r, p in pairs if lv in (p["level"], p["against"])})
        s["contrasts"] = sum(1 for _, p in pairs if lv in (p["level"], p["against"]))
        out["levels"][lv] = s
    for lv in levels:
        if lv not in comp:
            out["levels"][lv] = {"reading": "not connected"}
    return out


def reading(n_labs: int, s: Mapping[str, Any]) -> str:
    if n_labs < MIN_LABS:
        return "too few laboratories"
    ci = s["ci90"]
    if (ci[0] > 0 or ci[1] < 0) and s["p_same_sign_new_lab"] >= P_SAME_SIGN:
        return "answered by the literature"
    return "direction unclear"


def _width(s):
    return s["new_lab_pi90"][1] - s["new_lab_pi90"][0]


def lolo(rows) -> Optional[Dict[str, Any]]:
    labs = sorted({r["lab"] for r in rows})
    if len(labs) < LOLO_MIN_LABS:
        return None
    inside, errs, nulls = [], [], []
    for lab in labs:
        train = [r for r in rows if r["lab"] != lab]
        s = numeric_fit(train, PRIOR_SCALE)
        lo, hi = s["new_lab_pi90"]
        for r in rows:
            if r["lab"] == lab:
                half = 1.645 * math.sqrt(r["v"])
                inside.append(lo - half <= r["y"] <= hi + half)
                errs.append(abs(r["y"] - s["estimate"]))
                nulls.append(abs(r["y"]))
    cov = float(np.mean(inside))
    return {"labs": len(labs), "coverage_90": cov, "flag_too_narrow": cov < LOLO_FLAG,
            "median_abs_error": float(np.median(errs)), "median_abs_no_change": float(np.median(nulls))}


def cell_table(obs: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    cells = defaultdict(list)
    for o in obs:
        cells[(o["compound"], o["lever_class"], o["band"])].append(o)
    out = []
    for (compound, lc, band), rows in sorted(cells.items()):
        labs = sorted({r["lab"] for r in rows})
        base = {"compound": compound, "lever_class": lc, "band": band, "contrasts": len(rows),
                "studies": len({r["source"] for r in rows}), "labs": len(labs), "lab_names": labs,
                "peak_area_share": float(np.mean([r["peak_area"] for r in rows])),
                "focal_share": float(np.mean([r["focal"] for r in rows])),
                "matrices": sorted({m for r in rows for m in r["matrix"]}),
                "temperature_range_C": [min((t for r in rows for t in r["temps"]), default=None),
                                        max((t for r in rows for t in r["temps"]), default=None)],
                "rows": [{k: r.get(k) for k in ("source", "contrast", "lab", "y", "v", "range", "non_monotone",
                                                "pairs", "focal", "peak_area")} for r in rows]}
        if lc in NUMERIC:
            base["step"] = NUMERIC[lc][0]
            s = numeric_fit(rows, PRIOR_SCALE)
            sens = {str(sc): numeric_fit(rows, sc) for sc in SENSITIVITY_SCALES}
            s["prior_dominated"] = any(abs(_width(t) - _width(s)) / max(_width(s), 1e-9) > PRIOR_DOMINATED
                                       for t in sens.values())
            s["sensitivity"] = {k: {"estimate": v["estimate"], "new_lab_pi90": v["new_lab_pi90"]}
                                for k, v in sens.items()}
            for name, sub in (("without_peak_areas", [r for r in rows if not r["peak_area"]]),
                              ("incidental_only", [r for r in rows if not r["focal"]])):
                s[name] = ({"contrasts": len(sub), **{k: v for k, v in numeric_fit(sub, PRIOR_SCALE).items()
                                                      if k in ("estimate", "ci90")}} if sub else {"contrasts": 0})
            s["non_monotone_contrasts"] = sum(1 for r in rows if r["non_monotone"])
            base["fit"] = s
            base["reading"] = reading(len(labs), s)
            base["lolo"] = lolo(rows)
        else:
            ref = CATEGORICAL_REFERENCE[lc]
            f = categorical_fit(rows, ref, PRIOR_SCALE)
            base["reference"] = f.get("reference", ref)
            if f.get("reference_fallback"):
                base["reference_note"] = (f"declared reference {ref!r} absent from every study of this cell; reported "
                                          f"against {base['reference']!r} (PREREG amendment A3)")
            for lv, s in f["levels"].items():
                if "estimate" in s:
                    s["reading"] = reading(s["labs"], s)
            base["fit"] = f
            base["reading"] = "per level"
        out.append(base)
    return out


# ------------------------------------------------------------------------------------------------- rendering
def _f(x, nd=2):
    return "-" if x is None else f"{x:+.{nd}f}" if isinstance(x, float) else str(x)


def render(payload: Mapping[str, Any]) -> str:
    L = [f"# {payload['title']}", "",
         "*Generated by `scripts/generators/build_lever_table.py` under `docs/reference/lever_table/PREREG.md`. "
         "Every number is estimated from double-keyed within-study contrasts in the literature; no engine output "
         "enters. Read the new-laboratory interval: it is what a formulator should expect in their own pot.*", "",
         "**What this table cannot do.** Published contrasts are a selected sample (a lever that did nothing is "
         "under-reported), so a pooled effect is an upper estimate of the typical one. Most contrasts are model "
         "systems in buffer; a cell names its matrices and does not claim to transfer to a plant-protein product. "
         "A cell reports what happened, not why, and predicts no combination of levers that no study varied "
         "together.", "",
         f"- Sources {payload['counts']['sources']}, contrasts used {payload['counts']['observations']}, "
         f"laboratories {payload['counts']['labs']}; not used: "
         + ("; ".join(f"{k} {v}" for k, v in payload['skipped'].items()) or "none") + ".", ""]
    L += ["## Numeric levers", "",
          "| odorant | lever (step) | band | reading | estimate (90 % CI) | new laboratory 90 % | P(same sign) | "
          "labs / studies / contrasts | peak areas | matrices | LOLO coverage |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in payload["cells"]:
        if "step" not in c:
            continue
        s = c["fit"]
        flags = (" (prior-dominated)" if s.get("prior_dominated") else "") + \
                (f" ({s['non_monotone_contrasts']} non-monotone)" if s.get("non_monotone_contrasts") else "")
        lo = c.get("lolo")
        if c["reading"] == "too few laboratories":      # PREREG sec. 5: listed, nothing pooled
            raw = ", ".join(f"{r['y']:+.2f}" for r in c["rows"])
            L.append(f"| {c['compound']} | {c['lever_class']} ({c['step']}) | {c['band']} | {c['reading']} | "
                     f"not pooled; per contrast: {raw} | - | - | {c['labs']} / {c['studies']} / {c['contrasts']} | "
                     f"{c['peak_area_share']:.0%} | {', '.join(c['matrices'])} | - |")
            continue
        L.append(f"| {c['compound']} | {c['lever_class']} ({c['step']}) | {c['band']} | {c['reading']}{flags} | "
                 f"{_f(s['estimate'])} ({_f(s['ci90'][0])} to {_f(s['ci90'][1])}) | "
                 f"{_f(s['new_lab_pi90'][0])} to {_f(s['new_lab_pi90'][1])} | {s['p_same_sign_new_lab']:.2f} | "
                 f"{c['labs']} / {c['studies']} / {c['contrasts']} | {c['peak_area_share']:.0%} | "
                 f"{', '.join(c['matrices'])} | "
                 + ("-" if not lo else f"{lo['coverage_90']:.2f}") + " |")
    L += ["", "## Categorical levers (each level against the reference)", "",
          "| odorant | lever | level vs reference | reading | estimate (90 % CI) | new laboratory 90 % | "
          "P(same sign) | labs / contrasts |", "|---|---|---|---|---|---|---|---|"]
    for c in payload["cells"]:
        if "step" in c:
            continue
        for lv, s in c["fit"]["levels"].items():
            if "estimate" not in s:
                L.append(f"| {c['compound']} | {c['lever_class']} | {lv} vs {c['reference']} | not connected to "
                         f"the reference | - | - | - | - |")
                continue
            if s["reading"] == "too few laboratories":
                L.append(f"| {c['compound']} | {c['lever_class']} | {lv} vs {c['reference']} | too few laboratories | "
                         f"not pooled (network estimate {_f(s['estimate'])}, shown for orientation only) | - | - | "
                         f"{s['labs']} / {s['contrasts']} |")
                continue
            L.append(f"| {c['compound']} | {c['lever_class']} | {lv} vs {c['reference']} | {s['reading']} | "
                     f"{_f(s['estimate'])} ({_f(s['ci90'][0])} to {_f(s['ci90'][1])}) | "
                     f"{_f(s['new_lab_pi90'][0])} to {_f(s['new_lab_pi90'][1])} | {s['p_same_sign_new_lab']:.2f} | "
                     f"{s['labs']} / {s['contrasts']} |")
    L += ["", "Estimates are log10 ratios (dex) per lever step; +0.30 dex is about 2x. 'answered by the literature' "
          "needs at least three laboratories, a 90 % interval of the estimate that excludes zero, and P(a new "
          "laboratory sees the same sign) >= 0.8 (PREREG sec. 5).", ""]
    return "\n".join(L)


def build(obs, title, skipped, labs_meta) -> Dict[str, Any]:
    cells = cell_table(obs)
    return {"artifact": "lever_table", "title": title, "prereg": "docs/reference/lever_table/PREREG.md",
            "split": labs_meta, "skipped": skipped,
            "counts": {"sources": len({o["source"] for o in obs}), "observations": len(obs),
                       "labs": len({o["lab"] for o in obs})},
            "cells": cells}


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.parse_args(argv)
    sources = load_sources()
    labs, meta = lab_partition()
    obs, skipped = observations(sources, labs)
    inputs = sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml")) + [data_paths.CONTRAST_SPLIT]
    for out, sel, title in ((OUT_ALL, obs, "The lever table from data (all keyed contrasts)"),
                            (OUT_DEV, [o for o in obs if o["partition"] == "development"],
                             "The lever table from data, DEVELOPMENT laboratories only (step 5's rival)")):
        payload = build(sel, title, skipped, meta)
        payload["provenance"] = provenance.provenance_block(
            out.stem, generated_by="scripts/generators/build_lever_table.py",
            inputs=[p for p in inputs if p.name != "TEMPLATE.yml" and ".reader_" not in p.name])
        artifact_io.write_artifact(payload, out, render=render)
        print(f"wrote {out.relative_to(ROOT)} ({len(payload['cells'])} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
