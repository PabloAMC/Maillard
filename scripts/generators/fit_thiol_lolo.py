#!/usr/bin/env python
"""
fit_thiol_lolo.py -- leave-one-laboratory-out adequacy of the frozen thiol lane (plan step B8; roadmap rows 1 and 9).

Pre-registration: results/validation/kinetic_core_thiol_lolo_prereg.md (written before any fold was fitted). The
admitted structure R3En (EX-024) is refitted without each laboratory -- every one of its 16 coordinates, R3En's own
measurement model with every series' extra noise FIXED at s042_e_conv's, warm-started from the full optimum
(s042_r3en_fixed), converged (a budget-exhausted fit is continued; a fold that never converges gives no verdict) --
and the held-out laboratory's rows are predicted with a between-laboratory error per method class estimated on the
training laboratories only. Where the training rows do not pin a coordinate the full fit identified ("leak-prone"),
a reset arm refits with it at the arm's declared start. Scored per fold and pooled: held-out log10 error, coverage of
the 80 / 90 / 95 % intervals with their widths, and the paired held-out log density against two chemistry-free nulls
(typical yield; no change). The verdict follows the prereg's sec. 7 (the Acceptance rule).

A laboratory is the fit row's ``laboratory`` key (data/keys/laboratories.yml, one research group); a paper that
reprints another laboratory's data counts as that laboratory (``COPY_LAB``).

    python scripts/generators/fit_thiol_lolo.py --list                    # the folds, no engine run
    python scripts/generators/fit_thiol_lolo.py --fold jiangnan_zhang_x   # one fold (resumes from its checkpoint)
    python scripts/generators/fit_thiol_lolo.py --all                     # every unfinished fold, each in its own
                                                                          # process, then the summary
    python scripts/generators/fit_thiol_lolo.py --summary                 # the summary from the fold records
Writes results/validation/kinetic_core_thiol_lolo_<lab>[_<tag>].json and kinetic_core_thiol_lolo[_<tag>].{json,md}.
"""
from __future__ import annotations

import argparse
import contextlib
import json
import math
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows, fitting, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_deoxyosone_submodel as sm3  # noqa: E402
import fit_joint_deoxyosone_thiol as fj  # noqa: E402

ARM = "R3En"
STRECKER_SHARE = 0.042
FULL_FIT = "s042_r3en_fixed"        # EX-024: the admitted optimum
NOISE_FROM = "s042_e_conv"          # EX-024: the fixed series noise R3En was fitted under
BUDGET = 1200
FTOL = 1e-4
MAX_CONTINUATIONS = 3
LEVELS = (0.80, 0.90, 0.95)
GATE_LEVEL = 0.90
GATE_BAND = (0.80, 0.97)
#: a paper that reprints another laboratory's data (conflicts register, duplicate_confirmed) counts as the
#: ORIGINAL's laboratory: Kang 2026 (Zhejiang Gongshang) reprints Jiangnan's Zhai 2023 Food Chem data (CF-005)
COPY_LAB = {"kang2026": "jiangnan_zhang_x"}
NULLS = ("typical_yield", "no_change")
STEM = "kinetic_core_thiol_lolo"
PREREG = data_paths.VALIDATION_DIR / f"{STEM}_prereg.md"


# ----------------------------------------------------------------------------------------------------------------
# The fit table and the folds
# ----------------------------------------------------------------------------------------------------------------

def paper_of(row) -> str:
    return row.ledger_id.split("_")[0]


def lab_of(row) -> str:
    """The row's laboratory (one research group); a reprint counts as its original's laboratory."""
    return COPY_LAB.get(paper_of(row), row.laboratory)


def make_folds(rows: Sequence[Any], lab: Callable[[Any], str] = lab_of) -> Dict[str, Tuple[tuple, tuple]]:
    """laboratory -> (training rows, held-out rows): the held-out rows are exactly that laboratory's, the training
    rows every other laboratory's (row order kept)."""
    labs = sorted({lab(r) for r in rows})
    return {name: (tuple(r for r in rows if lab(r) != name), tuple(r for r in rows if lab(r) == name))
            for name in labs}


def _fit_record(suffix: str) -> Dict[str, Any]:
    return json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{suffix}.json")
                      .read_text(encoding="utf-8"))


def load_fit_table():
    """R3En's fit rows, selected exactly as fit_joint_deoxyosone_thiol.main does; refuses to run unless the count per
    data set x observable equals the full fit's record."""
    arms = fj.ARMS[ARM]
    predictable = set(sulfur.SULFUR_STATE_KEYS)
    jobs, rows_all, pots_all = {}, [], {}
    extra = ((fj.ACTZ_TABLE,) if "actz" in arms else ()) + ((fj.ENTRY_TABLE, fj.ENTRY_PH_TABLE) if "entry" in arms else ())

    def wanted(table, pot):
        if table == fj.ENTRY_TABLE:
            return pot.startswith(fj.ENTRY_FIT_PREFIX)
        if table == fj.ENTRY_PH_TABLE:
            return pot in fj.ENTRY_PH_HOLDOUT_POTS or pot.startswith(fj.ENTRY_PH_SCORE_PREFIX)
        return True

    for table in fj.TABLES + extra:
        pots, rows = fit_rows.load(table)
        for j in probe.plan(pots, rows):
            if "skip" in j or any(pots[j["pot"]]["initial_mmol_L"].get(h) for h in fj.HEXOSES) or not wanted(table, j["pot"]):
                continue
            jobs[f"{table}:{j['pot']}"] = j
            pots_all[f"{table}:{j['pot']}"] = pots[j["pot"]]
        rows_all += [r for r in rows if r.species in predictable or (table == fj.ENTRY_PH_TABLE and r.species == "pH")]
    rows_all = list(fj.dedupe(r for r in rows_all if fj.key_of(r) in jobs))

    def is_holdout(r):
        if r.table == fj.ENTRY_PH_TABLE:
            return r.pot in fj.ENTRY_PH_HOLDOUT_POTS
        return ((r.ledger_id.startswith("zhai2023") and float(pots_all[fj.key_of(r)]["t_c"]) >= 140.0)
                or r.ledger_id.startswith(fj.HOLDOUT_LEDGER_PREFIXES) or r.pot in fj.HOLDOUT_POTS)

    def is_score_only(r):
        if r.table == fj.ENTRY_PH_TABLE:
            return True
        if r.table == fj.ACTZ_TABLE:
            return float(pots_all[fj.key_of(r)]["label_ph"]) != fj.ACTZ_FIT_PH
        return r.ledger_id.startswith(fj.SCORE_ONLY_LEDGER_PREFIXES) or r.pot in fj.SCORE_ONLY_POTS

    data = tuple(r for r in rows_all if not is_holdout(r) and not is_score_only(r))
    counts: Dict[str, int] = defaultdict(int)
    for r in data:
        counts[f"{r.ledger_id} | {r.species}"] += 1
    expected = {k: v["n"] for k, v in _fit_record(FULL_FIT)["raw_at_optimum"].items()}
    if dict(counts) != expected:
        diff = sorted(set(counts.items()) ^ set(expected.items()))[:6]
        raise SystemExit(f"the rebuilt fit table differs from {FULL_FIT}'s: {diff}")
    return data, jobs, pots_all


@contextlib.contextmanager
def engine_state():
    """R3En's switches (as fit_joint_deoxyosone_thiol.main sets them for --arm R3En), restored on exit."""
    saved = (fj.STRECKER_SHARE, dict(fj.DDP_LAW_STATIC), dict(fj.EA_STATIC), dict(fj.TDP_EXIT_STATIC), sulfur.PH_COUPLING)
    arms = fj.ARMS[ARM]
    fj.STRECKER_SHARE = STRECKER_SHARE
    fj.DDP_LAW_STATIC = {"law": "none"} if "ddp_law_none" in arms else {}
    fj.EA_STATIC = {}
    fj.TDP_EXIT_STATIC = fj.tdp_exit_static(arms)
    sulfur.PH_COUPLING = "inline"
    try:
        yield
    finally:
        fj.STRECKER_SHARE, fj.DDP_LAW_STATIC, fj.EA_STATIC, fj.TDP_EXIT_STATIC, sulfur.PH_COUPLING = saved


# ----------------------------------------------------------------------------------------------------------------
# The measurement model and the predictive scoring
# ----------------------------------------------------------------------------------------------------------------

def level_of(row) -> str:
    """Absolute rows (sida, calibrated) share a between-laboratory level per paper x species (D4); every other row
    keeps its calibration series."""
    if float(sm3.SCALE_TAU[fj.method_of(row)]) == 0.0:
        return f"{paper_of(row)}|{row.species}"
    return fj.series_of(row)


def _model(extra_noise=None, between_tau=None, level=None) -> fitting.MeasurementModel:
    return fitting.MeasurementModel(series_of=fj.series_of, sigma_of=lambda r: r.sigma, method_of=fj.method_of,
                                    scale_tau=sm3.SCALE_TAU, extra_noise=dict(extra_noise or {}),
                                    level_of=level, between_tau=dict(between_tau or {}))


def training_terms(raw_train: np.ndarray, train: Sequence[Any], n_constants: int) -> Dict[str, Any]:
    """A model's pooled series noise per method class and its between-laboratory tau per absolute class, from its
    residuals on the training rows only (prereg sec. 4.5)."""
    ok = [i for i, v in enumerate(raw_train) if np.isfinite(v)]
    rows = tuple(train[i] for i in ok)
    r = np.asarray(raw_train, dtype=float)[ok]
    noise = fj.noise_with_floor(_model(), r, rows, n_constants)
    per_class: Dict[str, List[float]] = defaultdict(list)
    for series in sorted({fj.series_of(x) for x in rows}):
        method = fj.method_of(next(x for x in rows if fj.series_of(x) == series))
        per_class[method].append(float(noise.get(series, fj.NOISE_FLOOR_DEX)))
    class_noise = {m: float(statistics.median(v)) for m, v in per_class.items()}
    taus = fitting.estimate_between_tau(_model(noise, level=level_of), r, rows)
    return {"class_noise": class_noise, "between_tau": taus, "n_training_rows_scored": len(rows)}


def predictive_model(terms: Dict[str, Any], held: Sequence[Any], fixed_noise: Optional[Dict[str, float]] = None
                     ) -> Tuple[fitting.MeasurementModel, Dict[str, Any]]:
    """The held-out rows' predictive model: pooled class noise (or, for the in-sample sensitivity, the fixed noise),
    the training tau per absolute class, the largest training tau for an absolute class with none."""
    notes = {}
    noise = {}
    for row in held:
        s = fj.series_of(row)
        if fixed_noise is not None:
            noise[s] = float(fixed_noise.get(s, fj.NOISE_FLOOR_DEX))
        else:
            noise[s] = float(terms["class_noise"].get(fj.method_of(row), fj.NOISE_FLOOR_DEX))
    taus = dict(terms["between_tau"])
    for method in sorted({fj.method_of(r) for r in held}):
        if float(sm3.SCALE_TAU[method]) == 0.0 and method not in taus:
            taus[method] = max(taus.values()) if taus else 0.0
            notes[method] = f"no training tau; the fold's largest ({taus[method]:.3g}) used"
    return _model(noise, taus, level_of), {"between_tau_used": taus, "fallbacks": notes}


def row_intervals(pm: fitting.MeasurementModel, r: np.ndarray, rows: Sequence[Any]) -> List[Dict[str, Any]]:
    """Per row: the residual centred on what the model integrates out (free groups: their weighted mean) and its
    predictive variance -- one-level formulas of MeasurementModel.coverage; None for a row coverage does not score."""
    sigma = np.array([pm.sigma(x) for x in rows])
    groups: Dict[str, List[int]] = defaultdict(list)
    for i, x in enumerate(rows):
        groups[pm.level_key(x)].append(i)
    out: List[Dict[str, Any]] = [{"centred": None, "var": None} for _ in rows]
    for idx in groups.values():
        tau = pm.tau(pm.method_of(rows[idx[0]]))
        s = sigma[idx]
        if math.isinf(tau):
            if len(idx) < 2:
                continue
            w = 1.0 / s ** 2
            centred = r[idx] - float(np.sum(w * r[idx]) / np.sum(w))
            var = s ** 2 * (1.0 - w / np.sum(w))
        else:
            centred, var = r[idx], s ** 2 + tau * tau
        for j, i in enumerate(idx):
            out[i] = {"centred": float(centred[j]), "var": float(var[j])}
    return out


def score(pm: fitting.MeasurementModel, r: np.ndarray, rows: Sequence[Any]) -> Dict[str, Any]:
    """Log density, coverage with widths, and error statistics of one model on the paired held-out rows."""
    from scipy.stats import norm

    r = np.asarray(r, dtype=float)
    per_row = row_intervals(pm, r, rows)
    cov = {}
    for level in LEVELS:
        q = float(norm.ppf(0.5 + level / 2.0))
        flags = [abs(p["centred"]) <= q * math.sqrt(p["var"]) for p in per_row if p["var"] is not None]
        widths = [2.0 * q * math.sqrt(p["var"]) for p in per_row if p["var"] is not None]
        inside, scored = pm.coverage(r, rows, level)
        assert (inside, scored) == (sum(flags), len(flags)), "row intervals disagree with MeasurementModel.coverage"
        cov[f"{level:.2f}"] = {"inside": inside, "scored": scored,
                               "coverage": inside / scored if scored else None,
                               "median_width_dex": float(np.median(widths)) if widths else None}
    offsets = pm.offsets(r, rows)
    after = np.array([v - (offsets.get(pm.level_key(x), 0.0) if float(pm.scale_tau[pm.method_of(x)]) > 0.0 else 0.0)
                      for v, x in zip(r, rows)])
    return {"log_density": float(pm.log_density(r, rows)), "n": len(rows), "coverage": cov,
            "error": _error_stats(r, after, rows), "per_row": per_row}


def _stats(v: np.ndarray) -> Dict[str, float]:
    v = np.asarray(v, dtype=float)
    return {"n": int(v.size), "mean": float(np.mean(v)), "rms": float(np.sqrt(np.mean(v ** 2))),
            "median_abs": float(np.median(np.abs(v)))} if v.size else {"n": 0}


def _error_stats(raw: np.ndarray, after: np.ndarray, rows: Sequence[Any]) -> Dict[str, Any]:
    by_class: Dict[str, List[int]] = defaultdict(list)
    for i, x in enumerate(rows):
        by_class[fj.method_of(x)].append(i)
    return {"raw": _stats(raw), "after_level": _stats(after),
            "by_method": {m: {"raw": _stats(raw[idx]), "after_level": _stats(after[idx])} for m, idx in sorted(by_class.items())}}


# ----------------------------------------------------------------------------------------------------------------
# The null baselines (prereg sec. 5)
# ----------------------------------------------------------------------------------------------------------------

def limiting_precursor(row, pots_all) -> float:
    charge = pots_all[fj.key_of(row)]["initial_mmol_L"]
    return min(float(v) for v in charge.values() if v and float(v) > 0.0)


def null_residuals(kind: str, train: Sequence[Any], rows: Sequence[Any], pots_all) -> Tuple[np.ndarray, int]:
    """Predicted minus measured (dex) under a null fitted on ``train``; NaN where the species has no training row.
    Returns the residuals and the number of constants (one median per species)."""
    def shift(x):
        return math.log10(limiting_precursor(x, pots_all)) if kind == "typical_yield" else 0.0

    by_species: Dict[str, List[float]] = defaultdict(list)
    for x in train:
        if x.value > 0:
            by_species[x.species].append(math.log10(x.value) - shift(x))
    med = {s: float(statistics.median(v)) for s, v in by_species.items()}
    out = [med[x.species] + shift(x) - math.log10(x.value) if (x.species in med and x.value > 0) else float("nan")
           for x in rows]
    return np.array(out, dtype=float), len(med)


# ----------------------------------------------------------------------------------------------------------------
# One fold
# ----------------------------------------------------------------------------------------------------------------

def _out(lab: Optional[str], tag: str) -> Path:
    return data_paths.VALIDATION_DIR / (STEM + (f"_{lab}" if lab else "") + (f"_{tag}" if tag else "") + ".json")


def _write(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(provenance.dumps_compact(payload, default=float) + "\n", encoding="utf-8")


def _provenance(name: str) -> Dict[str, Any]:
    tables = fj.TABLES + (fj.ACTZ_TABLE, fj.ENTRY_TABLE, fj.ENTRY_PH_TABLE)
    return provenance.provenance_block(
        name, generated_by="scripts/generators/fit_thiol_lolo.py",
        inputs=[fit_rows.FIT_ROWS_DIR / f"{t}.yml" for t in tables]
               + [fit_rows.LABORATORIES, ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                  ROOT / "src/kinetic_core/fitting.py", ROOT / "scripts/generators/build_submodel_probe.py",
                  ROOT / "scripts/generators/fit_joint_deoxyosone_thiol.py", PREREG,
                  data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{FULL_FIT}.json",
                  data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{NOISE_FROM}.json"])


def converge(problem: fitting.Problem, start: np.ndarray, budget: int, max_continuations: int,
             history: List[Dict[str, Any]], keys) -> fitting.Fit:
    """fitting.fit, continued from its own end point while the budget runs out (prereg sec. 4.1)."""
    f_ = fitting.fit(problem, start, budget=budget, ftol=FTOL)
    history.append(f_.as_dict(keys))
    n = 0
    while f_.budget_exhausted and n < max_continuations:
        f_ = fitting.fit(problem, f_.theta, budget=budget, ftol=FTOL)
        history.append(f_.as_dict(keys))
        n += 1
    return f_


def run_fold(lab: str, budget: int, max_continuations: int, tag: str) -> Dict[str, Any]:
    t0 = time.time()
    data, jobs, pots_all = load_fit_table()
    folds = make_folds(data)
    if lab not in folds:
        raise SystemExit(f"unknown laboratory {lab!r}; folds: {sorted(folds)}")
    train, held = folds[lab]
    assert not ({lab_of(r) for r in train} & {lab}) and {lab_of(r) for r in held} == {lab}
    out = _out(lab, tag)
    state: Dict[str, Any] = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {}
    if state.get("complete"):
        print(f"{lab}: complete already ({out.name})")
        return state
    keys, kinds, lower, upper, start_declared = fj.coordinates(False, fj.ARMS[ARM])
    full = _fit_record(FULL_FIT)
    theta_full = np.array([full["fit"]["theta"][k] for k in keys], dtype=float)
    sigma_full = np.array([full["laplace"]["sigma"][k] for k in keys], dtype=float)
    ident_full = {k: bool(full["laplace"]["identified"][k]) and not bool(full["laplace"]["at_bound"][k]) for k in keys}
    fixed_noise = dict(_fit_record(NOISE_FROM)["series_extra_noise"])
    missing = {fj.series_of(r) for r in data} - set(fixed_noise)
    if missing:
        raise SystemExit(f"{NOISE_FROM}: no noise for series {sorted(missing)[:5]}")
    fit_model = _model(fixed_noise)
    raw_many = fj.make_raw_many(jobs, keys)

    def raw(theta, rows):
        return raw_many([theta], rows)[0]

    problem = fitting.Problem(keys=keys, kinds=kinds, lower=lower, upper=upper, rows=train,
                              residuals=lambda th, rs: fit_model.whiten(raw(th, rs), rs),
                              evaluate_many=lambda pts, rs: [fit_model.whiten(v, rs) for v in raw_many(pts, rs)],
                              group_of=sm3.paper_of, nuisance_dof=fit_model.nuisance_dof)
    state.update({"artifact": out.stem, "laboratory": lab, "arm": ARM, "full_fit": FULL_FIT, "noise_from": NOISE_FROM,
                  "budget": budget, "ftol": FTOL, "max_continuations": max_continuations, "tag": tag,
                  "training_rows": len(train), "held_out_rows": len(held),
                  "training_pots": len({fj.key_of(r) for r in train}), "held_out_pots": len({fj.key_of(r) for r in held}),
                  "training_laboratories": sorted({lab_of(r) for r in train}),
                  "held_out_papers": sorted({paper_of(r) for r in held})})

    def variant(name: str, start: np.ndarray) -> Dict[str, Any]:
        rec = state.get(name) or {}
        if rec.get("fit") and not rec["fit"]["budget_exhausted"]:
            return rec
        history = list(rec.get("history", []))
        resume = np.array([rec["fit"]["theta"][k] for k in keys]) if rec.get("fit") else start
        f_ = converge(problem, resume, budget, max_continuations, history, keys)
        rec = {"start": dict(zip(keys, map(float, start))), "history": history, "fit": f_.as_dict(keys),
               "converged": not f_.budget_exhausted}
        state[name] = rec
        _write(out, state)   # checkpoint: a restart skips a converged fit
        return rec

    def laplace_of(name: str) -> None:
        rec = state[name]
        if "laplace" in rec:
            return
        theta = np.array([rec["fit"]["theta"][k] for k in keys])
        lap = fitting.laplace(problem, theta)
        rec["laplace"] = lap
        rec["shift_dex"] = dict(zip(keys, map(float, theta - theta_full)))
        rec["shift_full_sigma"] = {k: float(d / s) if np.isfinite(s) and s > 0 else None
                                   for k, d, s in zip(keys, theta - theta_full, sigma_full)}
        _write(out, state)

    variant("warm", theta_full)
    laplace_of("warm")
    leak = [k for k in keys if ident_full[k]
            and (not state["warm"]["laplace"]["identified"][k] or state["warm"]["laplace"]["at_bound"][k])]
    state["leak_prone"] = leak
    if leak:
        reset_start = np.array([start_declared[i] if k in leak else theta_full[i] for i, k in enumerate(keys)])
        variant("reset", reset_start)
        laplace_of("reset")
    # ---- scoring (prereg secs. 4.4 - 6)
    thetas = {"warm": state["warm"]["fit"]["theta"], **({"reset": state["reset"]["fit"]["theta"]} if leak else {})}
    engine_held = {v: raw(np.array([t[k] for k in keys]), held) for v, t in thetas.items()}
    engine_train = {v: raw(np.array([t[k] for k in keys]), train) for v, t in thetas.items()}
    in_sample_held = raw(theta_full, held)
    nulls_held, nulls_train, n_null = {}, {}, {}
    for kind in NULLS:
        nulls_held[kind], n_null[kind] = null_residuals(kind, train, held, pots_all)
        nulls_train[kind], _ = null_residuals(kind, train, train, pots_all)
    all_held = list(engine_held.values()) + list(nulls_held.values())
    mask = np.all([np.isfinite(v) for v in all_held], axis=0)
    paired = tuple(r for r, m in zip(held, mask) if m)
    left_out = {"engine_non_finite": int(sum(not np.all([np.isfinite(engine_held[v][i]) for v in engine_held])
                                             for i in range(len(held)))),
                "null_no_answer": {k: int(np.sum(~np.isfinite(v))) for k, v in nulls_held.items()}}
    models: Dict[str, Dict[str, Any]] = {}
    for v in engine_held:
        terms = training_terms(engine_train[v], train, len(keys))
        pm, info = predictive_model(terms, paired)
        models[f"engine_{v}"] = {**score(pm, engine_held[v][mask], paired), "terms": terms, **info}
        # diagnostic, not a gate (prereg sec. 6 addendum): the engine on EVERY finite held-out row, the species no
        # null can answer (no training laboratory measured them) included
        ok = np.isfinite(engine_held[v])
        rows_ok = tuple(r for r, m in zip(held, ok) if m)
        pm_all, _ = predictive_model(terms, rows_ok)
        s_all = score(pm_all, engine_held[v][ok], rows_ok)
        models[f"engine_{v}_all_rows"] = {k: s_all[k] for k in ("log_density", "n", "coverage", "error")}
        if v == "warm":   # the in-sample-noise sensitivity (held-out series keep s042_e_conv's noise)
            pm_s, info_s = predictive_model(terms, paired, fixed_noise=fixed_noise)
            s = score(pm_s, engine_held[v][mask], paired)
            models["engine_warm_fixed_noise_sensitivity"] = {k: s[k] for k in ("log_density", "n", "coverage")}
    for kind in NULLS:
        terms = training_terms(nulls_train[kind], train, n_null[kind])
        pm, info = predictive_model(terms, paired)
        models[f"null_{kind}"] = {**score(pm, nulls_held[kind][mask], paired), "terms": terms, **info}
    in_sample = np.asarray(in_sample_held, dtype=float)
    state.update({
        "provenance": _provenance(out.stem),
        "paired_rows": len(paired), "left_out": left_out, "models": models,
        "paired_log_density_difference": {
            f"engine_{v}_minus_null_{kind}": models[f"engine_{v}"]["log_density"] - models[f"null_{kind}"]["log_density"]
            for v in engine_held for kind in NULLS},
        "in_sample_raw_at_full_optimum": _stats(in_sample[np.isfinite(in_sample)]),
        "rows": [{"ledger_id": r.ledger_id, "species": r.species, "t_min": float(r.t_min), "method": fj.method_of(r),
                  "engine_raw": {v: float(engine_held[v][i]) for v in engine_held},
                  "null_raw": {k: float(nulls_held[k][i]) for k in NULLS}, "paired": bool(mask[i])}
                 for i, r in enumerate(held)],
        "failed_pot_runs": len(fj.FAILURES), "failures_sample": fj.FAILURES[:10],
        "complete": all(state[v]["converged"] for v in thetas),
        "wall_seconds_this_run": round(time.time() - t0, 1),
    })
    _write(out, state)
    print(f"{lab}: warm {state['warm']['fit']['status']} ({sum(h['evaluations'] for h in state['warm']['history'])} ev), "
          f"leak-prone {leak}; paired rows {len(paired)}; complete {state['complete']}; {state['wall_seconds_this_run']} s")
    return state


# ----------------------------------------------------------------------------------------------------------------
# Summary and verdict (prereg secs. 6 - 7)
# ----------------------------------------------------------------------------------------------------------------

def _pooled(folds: Dict[str, Dict[str, Any]], pick: Callable[[Dict[str, Any]], str], model_of: Callable[[str], str]
            ) -> Dict[str, Any]:
    """Pooled coverage (summed counts) and the median width over the pooled rows, per level."""
    out = {}
    for level in LEVELS:
        from scipy.stats import norm
        q = float(norm.ppf(0.5 + level / 2.0))
        inside = scored = 0
        widths: List[float] = []
        for lab, rec in folds.items():
            m = rec["models"][model_of(pick(rec))]
            c = m["coverage"][f"{level:.2f}"]
            inside += c["inside"]
            scored += c["scored"]
            widths += [2.0 * q * math.sqrt(p["var"]) for p in m["per_row"] if p["var"] is not None]
        out[f"{level:.2f}"] = {"inside": inside, "scored": scored, "coverage": inside / scored if scored else None,
                               "median_width_dex": float(np.median(widths)) if widths else None}
    return out


def _difference(folds, pick, null: str) -> Dict[str, Any]:
    d = [rec["models"][f"engine_{pick(rec)}"]["log_density"] - rec["models"][f"null_{null}"]["log_density"]
         for rec in folds.values()]
    n_rows = sum(rec["paired_rows"] for rec in folds.values())
    mean = float(np.mean(d))
    se = float(np.std(d, ddof=1) / math.sqrt(len(d))) if len(d) > 1 else float("inf")
    verdict = "engine better" if mean > 2 * se else ("null better" if mean < -2 * se else "indistinguishable")
    return {"per_fold": dict(zip(folds, d)), "mean": mean, "se": se, "per_row_mean": float(np.sum(d) / n_rows),
            "verdict": verdict}


def summarise(tag: str) -> Dict[str, Any]:
    data, _, _ = load_fit_table()
    labs = sorted(make_folds(data))
    folds, missing = {}, []
    for lab in labs:
        p = _out(lab, tag)
        rec = json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
        if rec is None or "models" not in rec:
            missing.append(lab)
        else:
            folds[lab] = rec
    sets = {"warm": lambda rec: "warm", "reset_where_leak_prone": lambda rec: "reset" if rec.get("leak_prone") else "warm"}
    result: Dict[str, Any] = {"artifact": STEM + (f"_{tag}" if tag else ""), "provenance": _provenance(STEM),
                              "prereg": PREREG.name, "folds_expected": labs, "folds_missing": missing,
                              "folds_not_converged": [lab for lab, rec in folds.items() if not rec.get("complete")]}
    if folds:
        result["sets"] = {}
        for name, pick in sets.items():
            pooled = {m: _pooled(folds, pick, (lambda v, m=m: f"engine_{v}" if m == "engine" else f"null_{m}"))
                      for m in ("engine",) + NULLS}
            diffs = {null: _difference(folds, pick, null) for null in NULLS}
            c = pooled["engine"][f"{GATE_LEVEL:.2f}"]["coverage"]
            result["sets"][name] = {"pooled_coverage": pooled, "paired_log_density": diffs,
                                    "calibrated": c is not None and GATE_BAND[0] <= c <= GATE_BAND[1],
                                    "beats": {null: diffs[null]["verdict"] == "engine better" for null in NULLS}}
        s = result["sets"].values()
        calibrated = all(x["calibrated"] for x in s)
        beats_both = all(all(x["beats"].values()) for x in s)
        if missing or result["folds_not_converged"]:
            verdict = "INCOMPLETE"
        elif calibrated and beats_both:
            verdict = "ADEQUATE"
        elif calibrated:
            verdict = "CALIBRATED, NOT BETTER THAN A NULL"
        else:
            verdict = "NOT ADEQUATE"
        result["verdict"] = verdict
    out = _out(None, tag)
    _write(out, result)
    out.with_suffix(".md").write_text(render(result, folds) + "\n", encoding="utf-8")
    print(f"summary: {result.get('verdict', 'no folds')}; missing {missing}")
    return result


def _f(x, fmt="{:.2f}") -> str:
    return "-" if x is None else fmt.format(x)


def render(result: Dict[str, Any], folds: Dict[str, Dict[str, Any]]) -> str:
    lines = ["# Leave-one-laboratory-out adequacy of the frozen thiol lane (R3En)", "",
             f"*Generated by `scripts/generators/fit_thiol_lolo.py`; pre-registration `{PREREG.name}`. Residual = "
             "predicted minus measured (dex). Each fold refits R3En's 16 coordinates on the other laboratories under "
             f"{NOISE_FROM}'s fixed series noise, warm-started from {FULL_FIT}; held-out rows are scored with the "
             "training laboratories' pooled series noise and between-laboratory tau. 'reset' = leak-prone coordinates "
             "refitted from the arm's declared start.*", "",
             f"**Verdict (prereg sec. 7): {result.get('verdict', 'none')}**"
             + (f"; folds missing {result['folds_missing']}" if result["folds_missing"] else "")
             + (f"; not converged {result['folds_not_converged']}" if result["folds_not_converged"] else ""), ""]
    for name, s in (result.get("sets") or {}).items():
        lines += [f"## Pooled, set '{name}'", "",
                  "| model | 80 % cov (width dex) | 90 % cov (width dex) | 95 % cov (width dex) |", "|---|---|---|---|"]
        for m, pc in s["pooled_coverage"].items():
            lines.append(f"| {m} | " + " | ".join(
                f"{_f(pc[f'{lv:.2f}']['coverage'])} ({_f(pc[f'{lv:.2f}']['median_width_dex'])}; "
                f"{pc[f'{lv:.2f}']['inside']}/{pc[f'{lv:.2f}']['scored']})" for lv in LEVELS) + " |")
        lines += ["", "| engine minus null | mean over folds | se | per row | verdict |", "|---|---|---|---|---|"]
        lines += [f"| {null} | {d['mean']:+.1f} | {d['se']:.1f} | {d['per_row_mean']:+.3f} | {d['verdict']} |"
                  for null, d in s["paired_log_density"].items()]
        lines += ["", f"- calibrated (90 % coverage in {GATE_BAND}): {s['calibrated']}; beats: {s['beats']}", ""]
    lines += ["## Per fold", "",
              "| laboratory | train / held rows | paired | warm fit | leak-prone | variant | raw RMS held (in-sample) | "
              "after-level RMS | 90 % cov (width) | log density engine / typical / no change |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for lab, rec in folds.items():
        w = rec["warm"]["fit"]
        for v in ("warm", "reset"):
            m = rec["models"].get(f"engine_{v}")
            if m is None:
                continue
            c = m["coverage"]["0.90"]
            lines.append(
                f"| {lab} | {rec['training_rows']} / {rec['held_out_rows']} | {rec['paired_rows']} | "
                f"status {w['status']}, cost {w['cost']:.1f}, {sum(h['evaluations'] for h in rec['warm']['history'])} ev | "
                f"{', '.join(rec.get('leak_prone', [])) or '-'} | {v} | {m['error']['raw'].get('rms', float('nan')):.2f} "
                f"({rec['in_sample_raw_at_full_optimum'].get('rms', float('nan')):.2f}) | "
                f"{m['error']['after_level'].get('rms', float('nan')):.2f} | {_f(c['coverage'])} ({_f(c['median_width_dex'])}) | "
                f"{m['log_density']:.1f} / {rec['models']['null_typical_yield']['log_density']:.1f} / "
                f"{rec['models']['null_no_change']['log_density']:.1f} |")
    lines += ["", "## Engine on every finite held-out row (diagnostic, not a gate: includes species no null answers)", "",
              "| laboratory | variant | rows | 80 % cov (width) | 90 % cov (width) | 95 % cov (width) | raw RMS |",
              "|---|---|---|---|---|---|---|"]
    for lab, rec in folds.items():
        for v in ("warm", "reset"):
            m = rec["models"].get(f"engine_{v}_all_rows")
            if m is None:
                continue
            lines.append(f"| {lab} | {v} | {m['n']} | " + " | ".join(
                f"{_f(m['coverage'][f'{lv:.2f}']['coverage'])} ({_f(m['coverage'][f'{lv:.2f}']['median_width_dex'])})"
                for lv in LEVELS) + f" | {m['error']['raw'].get('rms', float('nan')):.2f} |")
    lines += ["", "## Per fold: error by method class (engine, warm; raw mean / RMS, after-level RMS)", ""]
    for lab, rec in folds.items():
        bm = rec["models"]["engine_warm"]["error"]["by_method"]
        lines.append(f"- **{lab}**: " + "; ".join(
            f"{m} n {e['raw']['n']}: {e['raw']['mean']:+.2f} / {e['raw']['rms']:.2f}, {e['after_level']['rms']:.2f}"
            for m, e in bm.items()))
    lines += ["", "## Per fold: between-laboratory tau and pooled series noise (engine, warm; from training rows)", ""]
    for lab, rec in folds.items():
        m = rec["models"]["engine_warm"]
        lines.append(f"- **{lab}**: tau {{{', '.join(f'{k}: {v:.2f}' for k, v in m['between_tau_used'].items())}}}; "
                     f"noise {{{', '.join(f'{k}: {v:.2f}' for k, v in m['terms']['class_noise'].items())}}}"
                     + (f"; fallbacks {m['fallbacks']}" if m["fallbacks"] else ""))
    lines += ["", "## Per fold: coordinate shifts from the full optimum (warm; dex)", ""]
    for lab, rec in folds.items():
        sh = rec["warm"]["shift_dex"]
        lines.append(f"- **{lab}**: " + ", ".join(f"{k} {v:+.2f}" for k, v in sh.items() if abs(v) >= 0.05))
    return "\n".join(lines)


# ----------------------------------------------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--list", action="store_true", help="print the folds and exit (no engine run)")
    g.add_argument("--fold", metavar="LAB", help="run one fold (resumes from its checkpoint)")
    g.add_argument("--all", action="store_true", help="every unfinished fold, each in its own process; then the summary")
    g.add_argument("--summary", action="store_true", help="the summary from the fold records")
    parser.add_argument("--budget", type=int, default=BUDGET)
    parser.add_argument("--max-continuations", type=int, default=MAX_CONTINUATIONS)
    parser.add_argument("--tag", default="", help="suffix for every output stem (e.g. smoke)")
    args = parser.parse_args(argv)
    if args.list:
        data, jobs, _ = load_fit_table()
        for lab, (train, held) in make_folds(data).items():
            print(f"{lab}: held {len(held)} rows / {len({fj.key_of(r) for r in held})} pots "
                  f"({', '.join(sorted({paper_of(r) for r in held}))}); train {len(train)} rows / "
                  f"{len({fj.key_of(r) for r in train})} pots")
        return 0
    if args.all:
        data, _, _ = load_fit_table()
        # smallest held-out laboratory first, so the largest (most expensive to predict) fold runs last
        order = sorted(make_folds(data).items(), key=lambda kv: len(kv[1][1]))
        failed = []
        for lab, _ in order:
            p = _out(lab, args.tag)
            if p.exists() and json.loads(p.read_text(encoding="utf-8")).get("complete"):
                print(f"{lab}: complete already")
                continue
            cmd = [sys.executable, str(Path(__file__).resolve()), "--fold", lab, "--budget", str(args.budget),
                   "--max-continuations", str(args.max_continuations)] + (["--tag", args.tag] if args.tag else [])
            print(f"== {lab} ({time.strftime('%H:%M:%S')})", flush=True)
            if subprocess.run(cmd, cwd=str(ROOT)).returncode != 0:
                failed.append(lab)
        with engine_state():
            summarise(args.tag)
        if failed:
            print(f"folds that exited with an error: {failed}")
        return 1 if failed else 0
    with engine_state():
        if args.summary:
            summarise(args.tag)
        else:
            run_fold(args.fold, args.budget, args.max_continuations, args.tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
