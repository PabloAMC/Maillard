"""
build_decision_relevance_step0 -- step 0 of the decision-relevance pre-registration (FROZEN v3.1).

What it computes (``results/validation/decision_relevance_prereg.md`` sec. 5.0, frozen in commit 9e02ce3a):

* **Stage 0c** (sec. 5.0.1): the engine's RATIO ERROR sigma_r between arms of one study, per lane and pooled, from
  residuals already on disk. A noise-model estimate on spent and seen data: it credits nothing, fits nothing, scores
  no rival and moves no shipped output.
* **Stage 0a** (sec. 5.0.2): each decision's ORACLE CEILING -- what a perfect design could show a formulator under
  that ratio error -- over a 512-world prior ensemble (the shipped envelope sampler, B1 / B13 at PR #108's points
  pinned at 5e61fed0), and the layer-1 reading of sec. 7.1.

Nothing here simulates a design, the early variant or the full run.

    python scripts/generators/build_decision_relevance_step0.py --stage 0c          # seconds, no engine
    python scripts/generators/build_decision_relevance_step0.py --stage 0a --workers 2   # hours, the engine
    python scripts/generators/build_decision_relevance_step0.py --render            # the .md from the .json
    python scripts/generators/build_decision_relevance_step0.py --check             # stage 0c and the .md current
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "decision_relevance_step0.json"
PREREG = "results/validation/decision_relevance_prereg.md"
FREEZE_COMMIT = "9e02ce3a"
CONTRAST_SCORES = data_paths.VALIDATION_DIR / "contrast_scores.json"
PANEL_SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
CONTRAST_BENCHMARK = data_paths.VALIDATION_DIR / "kinetic_core_contrast_benchmark_score_v1.json"
PR108_COMMIT = "5e61fed0"
PR108_POINTS = "results/validation/ph_balance_s3_rival_points.json"

OMEGA_DEX = 0.419          # H0 (EX-059); only for the derived rho diagnostic
TAU_ABS_DEX = 0.964
G_MIN = 5                  # sec. 5.0.1, item 3
N_UNSTATED = 2             # sec. 7 item 5: an unstated n is 2
CV_SIDA, CV_OTHER = 0.10, 0.20
Z90 = 1.6448536269514722
WORTH_DEX = 0.3
SIGN_DEX = 0.01
BAR_G, BAR_WRONG = 0.6, 0.1
N_WORLDS, SEED = 512, 20261107
DELTA_DRAWS = 1000
ARM_TIMEOUT_S = 300
EA_ISO_KJ_MOL = 76.0       # PR #108's ph_balance_refusal.EA_ISO_KJ_MOL
SHIPPED_SIGMA_R = math.sqrt(2.0) * OMEGA_DEX

LANE_OF_OUTPUT = {"MFT": "sulfur", "FFT": "sulfur", "FUR": "sulfur", "MFT/FFT": "sulfur", "HEX": "lipid",
                  "DMP": "trunk"}
TARGET = {"MFT": "2-methyl-3-furanthiol (MFT)", "FFT": "2-furfurylthiol (FFT)", "FUR": "furfural",
          "HEX": "hexanal", "DMP": "2,5-dimethylpyrazine"}


# ====================================================================================================== stage 0c
def _lstsq_arm_residuals(pairs: Sequence[Tuple[str, str, float]]) -> Dict[str, float]:
    """Arm residuals r_i (centred, sum 0) from pair residuals e_ab = r_a - r_b, by least squares on the pair graph.
    One connected component per call."""
    arms = sorted({a for a, _, _ in pairs} | {b for _, b, _ in pairs})
    idx = {a: i for i, a in enumerate(arms)}
    A = np.zeros((len(pairs) + 1, len(arms)))
    y = np.zeros(len(pairs) + 1)
    for k, (a, b, e) in enumerate(pairs):
        A[k, idx[a]], A[k, idx[b]], y[k] = 1.0, -1.0, e
    A[-1, :] = 1.0                       # the centring constraint
    r, *_ = np.linalg.lstsq(A, y, rcond=None)
    return {a: float(r[idx[a]]) for a in arms}


def _components(pairs: Sequence[Tuple[str, str, float]]) -> List[List[Tuple[str, str, float]]]:
    parent: Dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _ in pairs:
        parent[find(a)] = find(b)
    comp: Dict[str, List] = defaultdict(list)
    for p in pairs:
        comp[find(p[0])].append(p)
    return list(comp.values())


def _arm_measurement_var(cv: float, n: int = N_UNSTATED) -> float:
    return (math.log10(1.0 + cv)) ** 2 / n


def groups_from_pairs(pairs: Sequence[Mapping[str, Any]], arm_var: Mapping[str, float]) -> List[Dict[str, Any]]:
    """Within-study groups (study, compound, lane, out_of_sample) from pair records with keys study, compound, lane,
    a, b, e (= engine_dex - measured_dex), oos. Each connected component is one group."""
    bucket: Dict[Tuple, List[Tuple[str, str, float]]] = defaultdict(list)
    for p in pairs:
        bucket[(p["study"], p["compound"], p["lane"], bool(p["oos"]))].append((p["a"], p["b"], float(p["e"])))
    out = []
    for (study, compound, lane, oos), prs in sorted(bucket.items(), key=lambda kv: str(kv[0])):
        for comp in _components(prs):
            r = _lstsq_arm_residuals(comp)
            k = len(r)
            out.append({"study": study, "compound": compound, "lane": lane, "out_of_sample": oos, "k": k,
                        "ss": float(sum(v * v for v in r.values())), "df": k - 1,
                        "residuals": [float(v) for v in r.values()],
                        "measurement_var": [float(arm_var.get(a, _arm_measurement_var(CV_OTHER))) for a in r]})
    return out


def sigma_r_estimate(groups: Sequence[Mapping[str, Any]], *, robust: bool = False) -> Dict[str, Any]:
    """Sec. 5.0.1's estimator on a set of groups: pooled within-group variance, measurement part, the ratio error
    with its floor, the t-interval over study-level estimates, and leave-one-study-out."""
    from scipy import stats

    groups = [g for g in groups if g["df"] > 0]
    studies = sorted({g["study"] for g in groups})
    G = len(studies)
    if not groups:
        return {"studies": 0, "groups": 0, "status": "UNDECIDED (no groups)"}
    mvars = [v for g in groups for v in g["measurement_var"]]
    sigma_m2 = float(np.mean(mvars))
    floor = math.sqrt(2.0 * sigma_m2)

    def pooled(gs):
        df = sum(g["df"] for g in gs)
        return sum(g["ss"] for g in gs) / df if df else float("nan")

    def to_sigma_r(v):
        if not math.isfinite(v):
            return float("nan"), False
        raw2 = 2.0 * (v - sigma_m2)
        if raw2 <= floor ** 2:
            return floor, True
        return math.sqrt(raw2), False

    v_bar = pooled(groups)
    if robust:
        # sec. 5.0.1's robust form: 1.4826 x the MAD of the within-group residuals, each centred residual scaled by
        # sqrt(k / (k - 1)) so that it carries the arm-level variance
        adj = np.array([v * math.sqrt(g["k"] / g["df"]) for g in groups for v in g["residuals"]])
        v_bar = float((1.4826 * np.median(np.abs(adj))) ** 2)
    sigma_r, floored = to_sigma_r(v_bar)
    per_study = {}
    for s in studies:
        gs = [g for g in groups if g["study"] == s]
        per_study[s] = {"df": sum(g["df"] for g in gs), "v": pooled(gs), "groups": len(gs)}
    out: Dict[str, Any] = {
        "studies": G, "study_names": studies, "groups": len(groups), "df": sum(g["df"] for g in groups),
        "sigma_w2": v_bar, "sigma_m2": sigma_m2, "sigma_floor": floor, "sigma_r": sigma_r, "floored": floored,
        "rho_derived": 1.0 - sigma_r ** 2 / (2.0 * OMEGA_DEX ** 2), "per_study": per_study,
    }
    if G >= 2:
        w = np.array([per_study[s]["df"] for s in studies], float)
        w = w / w.sum()
        v = np.array([per_study[s]["v"] for s in studies], float)
        se = math.sqrt(G / (G - 1.0) * float(np.sum(w ** 2 * (v - v_bar) ** 2)))
        t = float(stats.t.ppf(0.95, G - 1))
        lo, hi = v_bar - t * se, v_bar + t * se
        out["t_interval_90"] = [to_sigma_r(lo)[0] if lo > 0 else floor, to_sigma_r(hi)[0]]
        loso = {}
        for s in studies:
            loso[s] = to_sigma_r(pooled([g for g in groups if g["study"] != s]))[0]
        out["loso"] = loso
        out["loso_min"], out["loso_max"] = min(loso.values()), max(loso.values())
    out["status"] = "ESTIMATED" if G >= G_MIN else f"UNDECIDED (G = {G} < G_min = {G_MIN})"
    return out


def _quant_cv(panel: Mapping[str, Any]) -> Dict[str, float]:
    """benchmark id -> the declared per-replicate CV from its quantification source (isotope dilution: SIDA)."""
    out = {}
    for b in panel.get("benchmarks", []):
        q = str(b.get("quantification_source") or "").lower()
        out[b["benchmark_id"]] = CV_SIDA if ("isotope" in q or "sida" in q or "stable" in q) else CV_OTHER
    return out


_SC = None


def _contrast_scorer():
    """scripts/generators/score_contrast_benchmark.py, for its pot specs (source 3); imported lazily."""
    global _SC
    if _SC is None:
        gen = str(ROOT / "scripts" / "generators")
        if gen not in sys.path:
            sys.path.insert(0, gen)
        import score_contrast_benchmark as sc   # noqa: E402

        _SC = sc
    return _SC


_TABLES: Dict[str, Any] = {}


def source3_spec(table: str, pid: str):
    sc = _contrast_scorer()
    if table not in _TABLES:
        _TABLES[table] = sc.load_table(table)
    return sc.pot_spec(table, pid, _TABLES[table]["pots"][pid])


def source3_lane(table: str, pid: str, species: str) -> str:
    """Amendment A1, fix (a): the engine lane that integrates `species` in this pot (declaration only, no run);
    the frozen generator hard-coded "sulfur"."""
    from src.kinetic_core import engine

    sc = _contrast_scorer()
    spec, tg = source3_spec(table, pid), sc.target_name(species)
    if spec is None or tg is None:
        return "unmapped"
    d = engine.declare_envelope(spec, [tg])
    lanes = list(d.lanes or ((d.lane or engine.TRUNK),))
    return next((ln for ln in lanes if ln in engine.MAILLARD_LANES), lanes[0] if lanes else "unmapped")


def collect_pairs(*, as_frozen: bool = False) -> Tuple[List[Dict[str, Any]], Dict[str, float], Dict[str, int]]:
    """The stage-0c pair units. ``as_frozen`` reproduces the generator as first run (2f1ef226), whose source 3
    was hard-coded to lane "sulfur" and out_of_sample True; the default applies amendment A1's fix (a): the engine's
    lane per pot and species, and out of sample when the row's exposure is "clean"."""
    cs = json.loads(CONTRAST_SCORES.read_text())
    panel = json.loads(PANEL_SCORES.read_text())
    cb = json.loads(CONTRAST_BENCHMARK.read_text())
    cv = _quant_cv(panel)
    arm_var = {bid: _arm_measurement_var(c) for bid, c in cv.items()}
    pairs, seen_arms = [], set()
    for p in cs["pairs"]:
        pairs.append({"study": p["study"], "compound": p["compound"], "lane": p["lane"], "a": p["a"], "b": p["b"],
                      "m": p["measured_dex"], "p": p["engine_dex"], "e": p["engine_dex"] - p["measured_dex"],
                      "oos": p["out_of_sample"], "source": 1})
        seen_arms |= {p["a"], p["b"]}
    n_src3, skipped3 = 0, 0
    for r in cb["primary"]["decision_set"]["rows"]:
        if r.get("m") is None or r.get("p") is None or r.get("kind") == "censored" or r["a"] in seen_arms:
            skipped3 += 1
            continue
        if r.get("sigma_m") is not None:
            for a in (r["a"], r["b"]):
                arm_var[a] = float(r["sigma_m"]) ** 2 / 2.0
        lane = "sulfur" if as_frozen else source3_lane(r["table"], r["a"], r["species"])
        oos = True if as_frozen else r.get("exposure") == "clean"
        pairs.append({"study": r["table"], "compound": r["species"], "lane": lane, "a": r["a"], "b": r["b"],
                      "m": float(r["m"]), "p": float(r["p"]), "e": float(r["p"]) - float(r["m"]), "oos": oos,
                      "source": 3})
        n_src3 += 1
    return pairs, arm_var, {"pairs_source1": len(cs["pairs"]), "pairs_source3": n_src3, "source3_skipped": skipped3}


def stage_0c(*, as_frozen: bool = False, pairs=None, arm_var=None, counts=None) -> Dict[str, Any]:
    if pairs is None:
        pairs, arm_var, counts = collect_pairs(as_frozen=as_frozen)
    n_src3, skipped3 = counts["pairs_source3"], counts["source3_skipped"]
    cs_pairs = counts["pairs_source1"]
    groups = groups_from_pairs(pairs, arm_var)
    oos = [g for g in groups if g["out_of_sample"]]
    lanes = sorted({g["lane"] for g in groups})
    result = {
        "units": {"pairs_source1": cs_pairs, "pairs_source3": n_src3, "source3_skipped": skipped3,
                  "source2": "no group of its own: contrast_scores.json is built from the same rows "
                             "(build_null_model_baselines._rows) with the same study key, so its groups are source 1's",
                  "groups_all": len(groups), "groups_out_of_sample": len(oos)},
        "measurement_rule": "no source here prints a per-arm SD: each arm takes the declared per-replicate CV "
                            "(SIDA 10 % where its benchmark's quantification source is isotope dilution, else 20 %) "
                            "with the unstated n = 2; a contrast-benchmark row that prints sigma_m (a ratio SD) gives "
                            "each arm sigma_m^2 / 2",
        "primary": {lane: sigma_r_estimate([g for g in oos if g["lane"] == lane]) for lane in lanes},
        "pooled": sigma_r_estimate(oos),
        "sensitivity_in_sample_included": {lane: sigma_r_estimate([g for g in groups if g["lane"] == lane])
                                           for lane in lanes},
        "sensitivity_robust": {lane: sigma_r_estimate([g for g in oos if g["lane"] == lane], robust=True)
                               for lane in lanes},
        "groups": groups,
    }
    return result


# ====================================================================================================== amendment A1
#: the two studies whose pH pairs the post-freeze check found scored under the falsified B0 balance
A1_FLAGGED_STUDIES = ("heldout_xu2007_tab1", "hofmann_schieberle_1998_free_mft_fft")
_BUNDLE_INDEX: Dict[str, Path] = {}


def _bundle_path(bid: str) -> Path:
    if not _BUNDLE_INDEX:
        for f in data_paths.BENCHMARKS_DIR.rglob("*.json"):
            _BUNDLE_INDEX.setdefault(f.stem, f)
    return _BUNDLE_INDEX[bid]


def unit_spec(source: int, study: str, arm: str, compound: str):
    """(FormulationSpec, engine target) for one arm of a stage-0c pair, built as its own scorer builds it."""
    from src import compound_keys, data_access
    from src.kinetic_core import panel

    if source == 3:
        return source3_spec(study, arm), _contrast_scorer().target_name(compound)
    bench = data_access.load_json(_bundle_path(arm))
    targets = list(panel.bundle_targets(bench))
    tg = next((t for t in targets if compound_keys.canonical_name(t) == compound), compound)
    return panel.core_spec(bench), tg


def _rescore_point(task: Mapping[str, Any]) -> Dict[str, Any]:
    from src.kinetic_core import config, engine

    changes, drift = ({}, None)
    if task["rival"] != "B0":
        points = _W["points"]
        point = points["rivals"][task["rival"]][task["point_index"]]
        changes, drift = _rival_changes(task["rival"], point,
                                        tuple(float(x) for x in points["davidek_breakpoints_in_situ"]))
    out: Dict[str, Optional[float]] = {}
    t0 = time.perf_counter()
    signal.signal(signal.SIGALRM, _alarm)
    with config.applied(config.Config.shipped().replace(changes)):
        for u in task["units"]:
            key = "|".join([str(u["source"]), u["study"], u["arm"], u["compound"]])
            try:
                spec, tg = unit_spec(u["source"], u["study"], u["arm"], u["compound"])
                if spec is None or tg is None:
                    out[key] = None
                    continue
                signal.alarm(ARM_TIMEOUT_S)
                kw = {} if drift is None else {"draw": engine.CoreDraw(ph_drift=drift)}
                run = engine.predict(spec, [tg], **kw)
                signal.alarm(0)
                mk = run.declaration.mapped_targets.get(tg) if run.answered else None
                out[key] = None if mk is None else float(run.species_mmol_per_l.get(mk, 0.0))
            except _Timeout:
                out[key] = None
            except Exception:  # noqa: BLE001  (recorded as a missing value and counted)
                signal.alarm(0)
                out[key] = None
    return {"point": task["name"], "rival": task["rival"], "values": out, "seconds": time.perf_counter() - t0}


def _init_rescore(points):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _W["points"] = points


def a1_rescore(pairs: Sequence[Mapping[str, Any]], workers: int, studies: Optional[Sequence[str]] = None
               ) -> Dict[str, Any]:
    """Re-score the out-of-sample sulfur-lane pairs (of ``studies``, or all) under B0 (a parity check against the
    on-disk values) and under each of PR #108's rival points at 5e61fed0."""
    import multiprocessing as mp

    points, sha = pr108_points()
    sel = [p for p in pairs if p["oos"] and p["lane"] == "sulfur" and (studies is None or p["study"] in studies)]
    units = sorted({(p["source"], p["study"], arm, p["compound"]) for p in sel for arm in (p["a"], p["b"])})
    unit_dicts = [{"source": s, "study": st, "arm": a, "compound": c} for s, st, a, c in units]
    tasks = [{"name": "B0", "rival": "B0", "point_index": None, "units": unit_dicts}]
    for rival, pts in points["rivals"].items():
        for k, pt in enumerate(pts):
            tasks.append({"name": pt["name"], "rival": rival, "point_index": k, "units": unit_dicts})
    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=_init_rescore, initargs=(points,)) as pool:
        results = list(pool.imap_unordered(_rescore_point, tasks))
    return {"points_sha256": sha, "pairs": len(sel), "units": len(units), "points": len(tasks) - 1,
            "wall_hours": (time.time() - t0) / 3600.0, "core_seconds": float(sum(r["seconds"] for r in results)),
            "results": sorted(results, key=lambda r: (r["rival"] != "B0", r["point"]))}


def _pair_p(p: Mapping[str, Any], values: Mapping[str, Optional[float]]) -> Optional[float]:
    def v(arm):
        return values.get("|".join([str(p["source"]), p["study"], arm, p["compound"]]))

    va, vb = v(p["a"]), v(p["b"])
    if not va or not vb or va <= 0 or vb <= 0:
        return None
    return math.log10(va / vb)


def a1_estimates(pairs, arm_var, counts, rescore, studies: Optional[Sequence[str]]) -> Dict[str, Any]:
    """sigma_r with the rescored pairs' engine log ratios from each rival point. The rival-scored reading pools,
    per group, the within-group sum of squares averaged over the 23 points (each stage-0a world is one point),
    with the other pairs unchanged."""
    by_point = [r for r in rescore["results"] if r["rival"] != "B0"]
    b0 = next(r for r in rescore["results"] if r["rival"] == "B0")
    parity = []
    target = [p for p in pairs if p["oos"] and p["lane"] == "sulfur" and (studies is None or p["study"] in studies)]
    for p in target:
        pb0 = _pair_p(p, b0["values"])
        if pb0 is not None:
            parity.append(abs(pb0 - p["p"]))
    per_point_groups, missing = [], 0
    for r in by_point:
        new = []
        for p in pairs:
            q = dict(p)
            if p in target:
                pr = _pair_p(p, r["values"])
                if pr is None:
                    missing += 1
                    continue
                q["p"], q["e"] = pr, pr - p["m"]
            new.append(q)
        per_point_groups.append(groups_from_pairs(new, arm_var))
    # average the sums of squares of the same group over the points where it is present
    acc: Dict[Tuple, Dict[str, Any]] = {}
    for gs in per_point_groups:
        for g in gs:
            key = (g["study"], g["compound"], g["lane"], g["out_of_sample"], g["k"])
            a = acc.setdefault(key, dict(g, ss_list=[], res_list=[]))
            a["ss_list"].append(g["ss"])
            a["res_list"].append(g["residuals"])
    groups = []
    for a in acc.values():
        g = {k: v for k, v in a.items() if k not in ("ss_list", "res_list")}
        g["ss"] = float(np.mean(a["ss_list"]))
        g["residuals"] = list(np.mean(np.array(a["res_list"]), axis=0)) if len({len(x) for x in a["res_list"]}) == 1 \
            else a["residuals"]
        groups.append(g)
    oos = [g for g in groups if g["out_of_sample"]]
    lanes = sorted({g["lane"] for g in groups})
    mean_p = {}
    for p in target:
        ps = [x for x in (_pair_p(p, r["values"]) for r in by_point) if x is not None]
        if ps:
            mean_p[(p["study"], p["compound"], p["a"], p["b"])] = float(np.mean(ps))
    thiol = [g for g in oos if g["lane"] == "sulfur" and is_thiol(g["compound"])]
    return {"parity_b0_max_abs_dex": max(parity) if parity else None, "parity_b0_n": len(parity),
            "pairs_rescored": len(target), "pair_point_values_missing": missing,
            "primary": {lane: sigma_r_estimate([g for g in oos if g["lane"] == lane]) for lane in lanes},
            "pooled": sigma_r_estimate(oos), "mean_rival_p": {"|".join(k): v for k, v in mean_p.items()},
            "thiol_only": {"primary": {"sulfur": sigma_r_estimate(thiol)},
                           "pooled": sigma_r_estimate(thiol + [g for g in oos if g["lane"] != "sulfur"])}}


def is_thiol(compound: str) -> bool:
    """A thiol or its disulfide (MFT, FFT, their dimers, methanethiol), by code or name."""
    c = str(compound)
    return c.upper() in {"MFT", "FFT", "MFTD", "FFTD", "MESH"} or any(
        k in c.lower() for k in ("thiol", "disulfide", "furanthiol"))


def calibration_diagnostics(pairs: Sequence[Mapping[str, Any]], drop: Sequence[str] = ()) -> Dict[str, Any]:
    """Amendment A1, item 3: the slope beta of observed on predicted log ratios (OLS with intercept; a 90 % interval
    from the study-clustered (CR1) standard error with t on G - 1 df) and the leave-one-study-out calibrated contrast
    (beta through the origin fitted on the other studies; median |m - beta p| against the no-change median |m|)."""
    from scipy import stats

    sel = [p for p in pairs if p["oos"] and p["lane"] == "sulfur" and p["study"] not in drop]
    if len(sel) < 3:
        return {"pairs": len(sel)}
    m = np.array([p["m"] for p in sel])
    x = np.array([p["p"] for p in sel])
    st = np.array([p["study"] for p in sel])
    X = np.column_stack([np.ones_like(x), x])
    beta = np.linalg.lstsq(X, m, rcond=None)[0]
    res = m - X @ beta
    bread = np.linalg.inv(X.T @ X)
    G = len(set(st))
    meat = sum(np.outer(X[st == s].T @ res[st == s], X[st == s].T @ res[st == s]) for s in set(st))
    n, k = len(m), 2
    cr1 = G / (G - 1) * (n - 1) / (n - k) if G > 1 else float("nan")
    se = math.sqrt(float((bread @ meat @ bread)[1, 1] * cr1)) if G > 1 else float("nan")
    t = float(stats.t.ppf(0.95, G - 1)) if G > 1 else float("nan")
    errs, nulls = [], []
    for s in sorted(set(st)):
        tr, te = st != s, st == s
        b = float(np.dot(x[tr], m[tr]) / np.dot(x[tr], x[tr])) if np.dot(x[tr], x[tr]) > 0 else 0.0
        errs += list(np.abs(m[te] - b * x[te]))
        nulls += list(np.abs(m[te]))
    return {"pairs": n, "studies": G, "beta": float(beta[1]), "beta_ci90": [float(beta[1] - t * se),
                                                                              float(beta[1] + t * se)],
            "loso_calibrated_median_abs_dex": float(np.median(errs)), "no_change_median_abs_dex": float(np.median(nulls))}


def sigma_r_for(lane: str, s0c: Mapping[str, Any]) -> Tuple[Optional[Dict[str, Any]], str]:
    """Sec. 5.0.1's rule: the lane's own sigma_r at >= G_min studies, else the pooled one, else UNDECIDED."""
    own = s0c["primary"].get(lane)
    if own and own.get("studies", 0) >= G_MIN:
        return own, f"{lane} (own, {own['studies']} studies)"
    pooled = s0c["pooled"]
    if pooled.get("studies", 0) >= G_MIN:
        return pooled, f"pooled ({pooled['studies']} studies; {lane}'s own: {own['studies'] if own else 0})"
    return None, "UNDECIDED (fewer than G_min studies, own and pooled)"


# ====================================================================================================== stage 0a
def oracle(L: np.ndarray, sigma_r: float, rng: np.random.Generator, draws: int = DELTA_DRAWS) -> Dict[str, Any]:
    """The oracle of sec. 5.0.2 on a vector of model levers L (one per world): shown when |L| > z sigma_r, correct
    when sign(L) = sign(T), T = L + delta, delta ~ N(0, sigma_r^2)."""
    L = np.asarray(L, float)
    L = L[np.isfinite(L)]
    if L.size == 0:
        return {"worlds": 0, "G": None, "P_right": None, "P_wrong": None, "eligible": False}
    T = L[:, None] + rng.normal(0.0, sigma_r, size=(L.size, draws))
    shown = (np.abs(L) > Z90 * sigma_r)[:, None] & np.ones_like(T, bool)
    correct = np.sign(T) == np.sign(L)[:, None]
    worth = np.abs(T) >= WORTH_DEX
    signed = np.abs(T) > SIGN_DEX
    n_worth = int(worth.sum())
    if n_worth == 0:
        return {"worlds": int(L.size), "G": None, "P_right": None, "P_wrong": None, "eligible": False,
                "note": "no truth worth acting on"}
    p_right = float((shown & correct & worth).sum() / n_worth)
    p_wrong_worth = float((shown & ~correct & worth).sum() / n_worth)
    p_wrong = float((shown & ~correct & signed).sum() / max(int(signed.sum()), 1))
    G = p_right - p_wrong_worth
    return {"worlds": int(L.size), "G": G, "P_right": p_right, "P_wrong_worth": p_wrong_worth, "P_wrong": p_wrong,
            "P_worth": float(worth.mean()), "eligible": bool(G >= BAR_G and p_wrong <= BAR_WRONG)}


def sigma_star(L: np.ndarray, seed: int) -> Optional[float]:
    """The largest sigma_r on a 400-point log grid in [0.01, 3] dex at which the oracle meets the bar."""
    best = None
    for s in np.geomspace(0.01, 3.0, 400):
        if oracle(L, float(s), np.random.default_rng(seed), draws=200)["eligible"]:
            best = float(s)
    return best


def _base_specs() -> Dict[str, Dict[str, Any]]:
    pea = {"matrix": "pea_isolate", "protein_g_per_l": 50.0, "protein_type": "pea_iso", "aw": 0.98}
    br_s = dict(name="BR_S", precursors={"L-Cysteine": 10.0, "D-Ribose": 10.0}, temp_C=120.0, time_min=60.0,
                ph=5.5, **pea)
    br_x = dict(name="BR_X", precursors={"L-Cysteine": 12.0, "D-Ribose": 12.0}, temp_C=125.0, time_min=45.0,
                ph=5.8, **pea)
    br_t = dict(name="BR_T", precursors={"D-Glucose": 50.0, "Glycine": 50.0}, temp_C=120.0, time_min=60.0,
                ph=6.5, **pea)

    def mod(base, name, pre=None, **kw):
        s = json.loads(json.dumps(base))
        s["name"] = name
        if pre:
            s["precursors"].update(pre)
        s.update(kw)
        return s

    return {
        "BR_S": br_s, "S_cys20": mod(br_s, "S_cys20", {"L-Cysteine": 20.0}),
        "BR_G": mod(br_s, "BR_G", None), "S_ph50": mod(br_s, "S_ph50", ph=5.0), "S_ph65": mod(br_s, "S_ph65", ph=6.5),
        "S_t120": mod(br_s, "S_t120", time_min=120.0), "S_T140": mod(br_s, "S_T140", temp_C=140.0),
        "S_ppi": mod(br_s, "S_ppi", {"pea protein isolate": 1.0}),
        "S_ppi_rooh": mod(br_s, "S_ppi_rooh", {"pea protein isolate": 1.0, "linoleate hydroperoxide": 1.0}),
        "S_gsh2": mod(br_s, "S_gsh2", {"glutathione": 2.0}), "S_thi2": mod(br_s, "S_thi2", {"thiamine": 2.0}),
        "BR_T": br_t, "T_ph75": mod(br_t, "T_ph75", ph=7.5),
        "BR_X": br_x, "X_cys18": mod(br_x, "X_cys18", {"L-Cysteine": 18.0}), "X_ph63": mod(br_x, "X_ph63", ph=6.3),
        "X_T135": mod(br_x, "X_T135", temp_C=135.0), "X_t90": mod(br_x, "X_t90", time_min=90.0),
        "X_cys5": mod(br_x, "X_cys5", {"L-Cysteine": 5.0}), "X_cys8": mod(br_x, "X_cys8", {"L-Cysteine": 8.0}),
        "X_ph70": mod(br_x, "X_ph70", ph=7.0), "X_100": mod(br_x, "X_100", temp_C=100.0, time_min=90.0),
        "X_100_gsh1": mod(br_x, "X_100_gsh1", {"glutathione": 1.0}, temp_C=100.0, time_min=90.0),
    }


def _fix_br_g(specs):
    g = specs["BR_G"]
    g["precursors"] = {"L-Cysteine": 10.0, "D-Glucose": 10.0}
    return specs


ARM_SPECS = _fix_br_g(_base_specs())
ARM_TARGETS = {a: (["DMP"] if a.startswith(("BR_T", "T_")) else ["HEX"] if a.startswith("S_ppi")
                   else ["MFT", "FFT"]) for a in ARM_SPECS}

#: (id, arm A, arm B, output)
DECISIONS: Tuple[Tuple[str, str, str, str], ...] = (
    ("D01", "BR_S", "S_cys20", "MFT"), ("D02", "BR_S", "S_cys20", "FFT"), ("D05", "BR_S", "BR_G", "FFT"),
    ("D06", "BR_S", "S_ph50", "FFT"), ("D07", "BR_S", "S_ph50", "MFT"), ("D08", "BR_S", "S_ph65", "MFT/FFT"),
    ("D09", "BR_S", "S_t120", "MFT"), ("D10", "BR_S", "S_T140", "FFT"), ("D14", "S_ppi", "S_ppi_rooh", "HEX"),
    ("D15", "BR_S", "S_gsh2", "MFT"), ("D16", "BR_S", "S_thi2", "MFT"), ("D18", "BR_T", "T_ph75", "DMP"),
    ("O01", "BR_X", "X_cys18", "MFT"), ("O02", "BR_X", "X_cys18", "FFT"), ("O03", "BR_X", "X_ph63", "FFT"),
    ("O04", "BR_X", "X_T135", "MFT"), ("O05", "BR_X", "X_t90", "MFT"), ("O08", "X_cys5", "X_cys8", "MFT"),
    ("O09", "BR_X", "X_ph70", "FFT"), ("O10", "X_100", "X_100_gsh1", "MFT"),
)
NOT_COMPUTABLE = {
    "D03": "the falsified hexose MFT route; Q-hex not built (A2b)",
    "D04": "the falsified hexose MFT route; Q-hex not built (A2b)",
    "D11": "Q-buf: R0 falsified at pH 5 - 7, the buffer switch not built",
    "D12": "Q-buf: R0 falsified at pH 5 - 7, the buffer switch not built",
    "D13": "Q-fat: R0 falsified, no rival built",
    "O06": "Q-buf: R0 falsified at pH 5 - 7, the buffer switch not built",
    "O07": "Q-fat: R0 falsified, no rival built",
    "D17": "structural refusal (no a_w term on the sulfur lane)",
    "D19": "structural refusal (no Strecker species)",
}


def pr108_points() -> Tuple[Dict[str, Any], str]:
    raw = subprocess.run(["git", "show", f"{PR108_COMMIT}:{PR108_POINTS}"], cwd=ROOT, check=True,
                         capture_output=True).stdout
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def _rival_changes(rival: str, point: Mapping[str, Any], breakpoints: Tuple[float, ...]):
    from src.kinetic_core import ph_balance_rivals as pbr
    from src.kinetic_core.ph_state import PhDrift

    v = point["coordinates"]
    changes: Dict[str, Any] = {"ph.coupling": "inline",
                               "ph.amine_survival": pbr.AmineSurvival(True, float(v["phi_cys"]), float(v["phi_arp"]))}
    if rival == "B13":
        changes["ph.acid_supply"] = pbr.AcidSupply(True, float(v["y5"]), float(v["y6"]), float(v["y_hi"]),
                                                   breakpoints, float(v["k_iso"]), EA_ISO_KJ_MOL, "hold_total")
        drift = PhDrift(acid_yield=0.0, arp_amine_pka=float(v["arp_amine_pka"]))
    else:
        drift = PhDrift(acid_yield=float(v["acid_yield"]), arp_amine_pka=float(v["arp_amine_pka"]))
    return changes, drift


def world_plan(n: int, points: Mapping[str, Any]) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(SEED + 1)
    plan = []
    for w in range(n):
        rival = "B1" if w % 2 == 0 else "B13"
        k = int(rng.integers(len(points["rivals"][rival])))
        plan.append({"world": w, "rival": rival, "point": points["rivals"][rival][k]["name"], "point_index": k})
    return plan


_W: Dict[str, Any] = {}


def _init_worker(points: Mapping[str, Any]):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from src.kinetic_core import uncertainty

    _W["draws"] = uncertainty.sample_draws(N_WORLDS, seed=SEED)
    _W["points"] = points


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def run_world(item: Mapping[str, Any]) -> Dict[str, Any]:
    import dataclasses

    from src.comparative_cli import spec_to_core, validate_spec
    from src.kinetic_core import config, engine

    points = _W["points"]
    draw = _W["draws"][item["world"]]
    point = points["rivals"][item["rival"]][item["point_index"]]
    changes, drift = _rival_changes(item["rival"], point, tuple(float(x) for x in points["davidek_breakpoints_in_situ"]))
    core = dataclasses.replace(draw.core, ph_drift=drift)
    values: Dict[str, Dict[str, Optional[float]]] = {}
    failed: Dict[str, str] = {}
    t0 = time.perf_counter()
    signal.signal(signal.SIGALRM, _alarm)
    with config.applied(config.Config.shipped().replace(changes)):
        for arm, spec in ARM_SPECS.items():
            targets = [TARGET[t] for t in ARM_TARGETS[arm]]
            try:
                signal.alarm(ARM_TIMEOUT_S)
                run = engine.predict(spec_to_core(validate_spec(spec, label=arm)), targets, draw=core,
                                     size_declared_bands=False)
                signal.alarm(0)
            except _Timeout:
                failed[arm] = "timeout"
                continue
            except Exception as exc:  # noqa: BLE001  (a failed arm is recorded, never silently read)
                signal.alarm(0)
                failed[arm] = f"{type(exc).__name__}: {exc}"[:200]
                continue
            if not run.answered:
                failed[arm] = "refused: " + "; ".join(run.declaration.reasons)[:200]
                continue
            vals = {}
            for t in ARM_TARGETS[arm]:
                key = run.declaration.mapped_targets.get(TARGET[t])
                vals[t] = None if key is None else float(run.species_mmol_per_l.get(key, 0.0))
            values[arm] = vals
    return {"world": item["world"], "rival": item["rival"], "point": item["point"], "values": values,
            "failed": failed, "seconds": time.perf_counter() - t0}


def lever(rec: Mapping[str, Any], a: str, b: str, out: str) -> Optional[float]:
    va, vb = rec["values"].get(a), rec["values"].get(b)
    if va is None or vb is None:
        return None

    def lg(v, k):
        x = v.get(k)
        return math.log10(x) if x is not None and x > 0 else None

    if out == "MFT/FFT":
        parts = [lg(vb, "MFT"), lg(vb, "FFT"), lg(va, "MFT"), lg(va, "FFT")]
        return None if None in parts else (parts[0] - parts[1]) - (parts[2] - parts[3])
    x, y = lg(vb, out), lg(va, out)
    return None if x is None or y is None else x - y


def stage_0a(s0c: Mapping[str, Any], workers: int, n: int = N_WORLDS, cache: Optional[Path] = None) -> Dict[str, Any]:
    import multiprocessing as mp

    points, points_sha = pr108_points()
    plan = world_plan(n, points)
    t0 = time.time()
    records: List[Dict[str, Any]] = []
    done = set()
    if cache and cache.exists():
        for line in cache.read_text().splitlines():
            r = json.loads(line)
            records.append(r)
            done.add(r["world"])
    todo = [p for p in plan if p["world"] not in done]
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init_worker, initargs=(points,)) as pool:
        for r in pool.imap_unordered(run_world, todo):
            records.append(r)
            if cache:
                with cache.open("a") as fh:
                    fh.write(json.dumps(r) + "\n")
            print(f"world {r['world']} {r['rival']} {r['seconds']:.0f} s, failed {len(r['failed'])}", flush=True)
    records.sort(key=lambda r: r["world"])
    wall = time.time() - t0
    return summarise_0a(records, s0c, points_sha, wall, n)


def summarise_0a(records, s0c, points_sha, wall, n) -> Dict[str, Any]:
    decisions = {}
    sigma_star_list = []
    levers: Dict[str, np.ndarray] = {}
    for did, a, b, out in DECISIONS:
        L = np.array([v for v in (lever(r, a, b, out) for r in records) if v is not None], float)
        levers[did] = L
        lane = LANE_OF_OUTPUT[out]
        est, rule = sigma_r_for(lane, s0c)
        row: Dict[str, Any] = {"arms": [a, b], "output": out, "lane": lane, "sigma_r_rule": rule,
                               "worlds_with_lever": int(L.size), "worlds_missing": int(len(records) - L.size)}
        if L.size:
            row["L_quantiles"] = {q: float(np.quantile(L, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
            row["P_abs_L_over_shipped_halfwidth"] = float(np.mean(np.abs(L) > Z90 * SHIPPED_SIGMA_R))
        seed = int(hashlib.sha256(did.encode()).hexdigest()[:8], 16)
        if est is None or not L.size:
            row["eligible"] = None
            row["status"] = "UNDECIDED" if est is None else "no lever computed"
        else:
            at = {"point": est["sigma_r"], "loso_min": est.get("loso_min", est["sigma_r"]),
                  "loso_max": est.get("loso_max", est["sigma_r"])}
            row["oracle"] = {k: dict(oracle(L, s, np.random.default_rng(seed)), sigma_r=s) for k, s in at.items()}
            elig = {k: v["eligible"] for k, v in row["oracle"].items()}
            row["eligible"] = elig["point"]
            row["eligible_agrees_across_loso"] = len(set(elig.values())) == 1
            row["status"] = "eligible" if elig["point"] else "omega-bound"
        row["oracle_as_shipped"] = dict(oracle(L, SHIPPED_SIGMA_R, np.random.default_rng(seed)),
                                        sigma_r=SHIPPED_SIGMA_R) if L.size else None
        row["sigma_star"] = sigma_star(L, seed) if L.size else None
        if est is not None and row["sigma_star"]:
            row["reduction_factor_needed"] = est["sigma_r"] / row["sigma_star"]
        decisions[did] = row
        if did.startswith("D"):
            sigma_star_list.append(row["sigma_star"] or 0.0)
    d_rows = {k: v for k, v in decisions.items() if k.startswith("D")}
    o_rows = {k: v for k, v in decisions.items() if k.startswith("O")}
    E = sum(1 for v in d_rows.values() if v.get("eligible"))
    undecided_d = [k for k, v in d_rows.items() if v.get("eligible") is None]
    loso_flip = [k for k, v in d_rows.items() if v.get("eligible_agrees_across_loso") is False]
    ss = sorted(sigma_star_list, reverse=True)
    majority_k = math.ceil(len(d_rows) / 2)
    e_at = {k: sum(1 for v in d_rows.values() if ((v.get("oracle") or {}).get(k) or {}).get("eligible"))
            for k in ("point", "loso_min", "loso_max")}
    # prereg sec. 7.1's sensitivities sigma_r in {0.593 (as shipped), 0.42, 0.19} (amendment A1, fix (b))
    e_fixed = {}
    for s_fix in (round(SHIPPED_SIGMA_R, 3), 0.42, 0.19):
        elig = [d for d in d_rows if levers[d].size and oracle(
            levers[d], s_fix, np.random.default_rng(int(hashlib.sha256(d.encode()).hexdigest()[:8], 16)))["eligible"]]
        e_fixed[f"{s_fix:.3f}"] = {"E": len(elig), "eligible": elig}
    if undecided_d:
        reading = (f"UNDECIDED: sigma_r is UNDECIDED (fewer than G_min = {G_MIN} studies) for "
                   f"{', '.join(undecided_d)}; step 0 reading (envelope incomplete)")
    elif len({n_e < 2 for n_e in e_at.values()}) > 1:
        reading = (f"UNDECIDED: the layer-1 reading differs across the leave-one-study-out bracket "
                   f"(eligible at point / LOSO_min / LOSO_max: {e_at['point']} / {e_at['loso_min']} / "
                   f"{e_at['loso_max']})")
    elif E < 2:
        reading = ("CONFIRMED (omega-bound): the lever table cannot show directions under the current discrepancy, "
                   "whatever the experiments.")
    else:
        reading = (f"NOT omega-bound at step 0: {E} of {len(d_rows)} computable decisions eligible "
                   f"(on {E} of 17); the design comparison (layer 2) is informative")
    failed_arms = sum(len(r["failed"]) for r in records)
    return {
        "worlds": len(records), "worlds_planned": n, "points_record": {"commit": PR108_COMMIT, "path": PR108_POINTS,
                                                                        "sha256": points_sha},
        "rivals": {k: sum(1 for r in records if r["rival"] == k) for k in ("B1", "B13")},
        "failed_arm_runs": failed_arms, "arm_runs": len(records) * len(ARM_SPECS),
        "failed_examples": sorted({v for r in records for v in r["failed"].values()})[:10],
        "wall_hours": wall / 3600.0, "core_seconds": float(sum(r["seconds"] for r in records)),
        "decisions": decisions, "not_computable": NOT_COMPUTABLE,
        "eligible": sorted(k for k, v in d_rows.items() if v.get("eligible")),
        "eligible_off_grid": sorted(k for k, v in o_rows.items() if v.get("eligible")),
        "E_step0": E, "E_at": e_at, "E_at_fixed_sigma_r": e_fixed, "computable_D": len(d_rows),
        "loso_flips": loso_flip,
        "sigma_star_for_2_eligible": ss[1] if len(ss) > 1 else None,
        "sigma_star_for_majority": ss[majority_k - 1] if len(ss) >= majority_k else None,
        "layer1_step0_reading": reading,
    }


# ====================================================================================================== record
def _f(x, p=3):
    return "-" if x is None else (f"{x:.{p}f}" if isinstance(x, float) else str(x))


def _rho_text(e: Mapping[str, Any]) -> str:
    r = e.get("rho_derived")
    if r is None:
        return "-"
    if r < -1.0:
        return "impossible: omega and sigma_r disagree (the engine over-responds to levers)"
    return f"{r:.2f}"


def _sigma_table(c: Mapping[str, Any]) -> List[str]:
    L = ["| lane | studies | groups | df | sigma_w^2 | sigma_m^2 | floor | **sigma_r** (dex) | 90 % t-interval | "
         "LOSO min - max | rho (derived diagnostic) | status |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for lane, e in list(c["primary"].items()) + [("POOLED", c["pooled"])]:
        if not e.get("groups"):
            L.append(f"| {lane} | 0 | 0 | | | | | | | | | {e['status']} |")
            continue
        ti = e.get("t_interval_90")
        L.append(f"| {lane} | {e['studies']} | {e['groups']} | {e['df']} | {_f(e['sigma_w2'])} | "
                 f"{_f(e['sigma_m2'], 4)} | {_f(e['sigma_floor'])} | **{_f(e['sigma_r'])}**"
                 f"{' (floor)' if e['floored'] else ''} | {('%.3f - %.3f' % tuple(ti)) if ti else '-'} | "
                 f"{('%.3f - %.3f' % (e['loso_min'], e['loso_max'])) if 'loso_min' in e else '-'} | "
                 f"{_rho_text(e)} | {e['status']} |")
    L += ["", "Studies per lane: " + "; ".join(
        f"{lane}: {', '.join(e.get('study_names', []))}" for lane, e in c["primary"].items()), ""]
    s = c["primary"].get("sulfur") or {}
    if s.get("per_study"):
        L += ["Per-study sigma_r, sulfur lane (sqrt(2 (v_s - sigma_m^2))): " + "; ".join(
            f"{st} {math.sqrt(max(2.0 * (v['v'] - s['sigma_m2']), 0.0)):.2f} (df {v['df']})"
            for st, v in s["per_study"].items()) + ".", ""]
    return L


def _q(q, key):
    return q.get(key) if key in q else q.get(str(key))


def _decision_table(a: Mapping[str, Any]) -> List[str]:
    L = ["| decision | arms | output | sigma_r rule | worlds | L 5 / 50 / 95 % (dex) | G* (point) | P_wrong | "
         "eligible | same across LOSO | G* as shipped | sigma* | reduction needed |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for did, r in a["decisions"].items():
        q = r.get("L_quantiles") or {}
        o = (r.get("oracle") or {}).get("point") or {}
        sh = r.get("oracle_as_shipped") or {}
        status = r.get("status")
        if did == "D14" and r.get("eligible"):
            status = "eligible (trivial: adding ROOH makes hexanal)"
        L.append(f"| {did} | {r['arms'][0]} -> {r['arms'][1]} | {r['output']} | {r['sigma_r_rule']} | "
                 f"{r['worlds_with_lever']} | {_f(_q(q, 0.05), 2)} / {_f(_q(q, 0.5), 2)} / {_f(_q(q, 0.95), 2)} | "
                 f"{_f(o.get('G'))} | {_f(o.get('P_wrong'))} | {status} | "
                 f"{r.get('eligible_agrees_across_loso', '-')} | {_f(sh.get('G'))} | {_f(r.get('sigma_star'))} | "
                 f"{_f(r.get('reduction_factor_needed'), 2)} |")
    ef = a.get("E_at_fixed_sigma_r") or {}
    trivial = " (D14's eligibility is trivial: adding ROOH makes hexanal)" if "D14" in a["eligible"] else ""
    L += ["", f"Eligible (of the {a['computable_D']} computable D-decisions): {len(a['eligible'])} of "
              f"{a['computable_D']}: {', '.join(a['eligible']) or 'none'}{trivial}; off-grid eligible: "
              f"{', '.join(a['eligible_off_grid']) or 'none'}. Eligible at point / LOSO_min / LOSO_max: "
              f"{a['E_at']['point']} / {a['E_at']['loso_min']} / {a['E_at']['loso_max']}.",
          "E at the pre-registered sigma_r sensitivities (sec. 7.1): " + "; ".join(
              f"sigma_r = {k}: E = {v['E']} ({', '.join(v['eligible']) or 'none'})" for k, v in ef.items()) + ".",
          f"The ratio error at which >= 2 D-decisions would be eligible: {_f(a['sigma_star_for_2_eligible'])} dex; "
          f"a majority of the computable ones: {_f(a['sigma_star_for_majority'])} dex.", "",
          f"> Layer 1 (sec. 7.1), step 0 reading (envelope incomplete; the full run's recomputation binds): "
          f"**{a['layer1_step0_reading']}**", ""]
    return L


def render(payload: Mapping[str, Any]) -> str:
    L = ["# Decision relevance, step 0: the ratio error and the oracle ceiling", "",
         f"*Generated by `scripts/generators/build_decision_relevance_step0.py` under the FROZEN v3.1 pre-registration "
         f"(`{PREREG}`, frozen in {FREEZE_COMMIT}), sec. 5.0, with its post-freeze amendment A1 (below). Step 0 only: "
         f"no design, early variant or full run is simulated. Stage 0c is a noise-model estimate on spent and seen "
         f"data; it credits nothing and moves no shipped output.*", ""]
    a1 = payload.get("amendment_A1")
    if a1:
        b = a1["readings"][a1["binding"]]
        L += ["## The binding result (amendment A1)", "",
              f"> **{b['stage_0a']['layer1_step0_reading']}** (binding reading: {a1['binding_label']}; sulfur "
              f"sigma_r = {_f(b['sigma']['primary']['sulfur']['sigma_r'])} dex, eligible "
              f"{len(b['stage_0a']['eligible'])} of {b['stage_0a']['computable_D']}: "
              f"{', '.join(b['stage_0a']['eligible']) or 'none'})", "",
              "| reading | sulfur sigma_r (dex) | LOSO | pooled | E (point / LOSO_min / LOSO_max) | eligible | layer 1 |",
              "|---|---|---|---|---|---|---|"]
        for key, r in a1["readings"].items():
            s, p = r["sigma"]["primary"]["sulfur"], r["sigma"]["pooled"]
            z = r["stage_0a"]
            L.append(f"| {r['label']}{' **(binding)**' if key == a1['binding'] else ''} | {_f(s['sigma_r'])} | "
                     f"{_f(s.get('loso_min'))} - {_f(s.get('loso_max'))} | {_f(p['sigma_r'])} | "
                     f"{z['E_at']['point']} / {z['E_at']['loso_min']} / {z['E_at']['loso_max']} | "
                     f"{', '.join(z['eligible']) or 'none'} | {z['layer1_step0_reading']} |")
        L += ["", "## Amendment A1 (post-freeze correction, 2026-11-07)", ""] + [f"{x}" for x in a1["declaration"]] + [""]
        rs = a1["rescore"]
        L += [f"Re-scoring: {rs['pairs']} out-of-sample sulfur-lane pairs ({rs['units']} arms) under B0 and PR #108's "
              f"{rs['points']} rival points (sha256 {rs['points_sha256'][:12]}); {rs['core_seconds'] / 3600:.2f} core-h, "
              f"{rs['wall_hours']:.2f} h wall at 2 processes. B0 parity with the on-disk engine values: max "
              f"{_f(a1['readings']['rival_all']['parity_b0_max_abs_dex'], 4)} dex over "
              f"{a1['readings']['rival_all']['parity_b0_n']} pairs.", ""]
        for key, r in a1["readings"].items():
            L += [f"### Reading: {r['label']}", ""] + _sigma_table(r["sigma"]) + _decision_table(r["stage_0a"])
        L += ["### Lever skill out of sample (A1 item 3; sulfur lane, out-of-sample pairs)", "",
              "| scoring | pairs | studies | slope beta of observed on predicted | 90 % CI (CR1, t on G - 1) | "
              "LOSO calibrated median abs error | no-change median abs |", "|---|---|---|---|---|---|---|"]
        for lab, d in a1["calibration"].items():
            ci = d.get("beta_ci90") or [None, None]
            L.append(f"| {lab} | {d.get('pairs')} | {d.get('studies', '-')} | {_f(d.get('beta'), 2)} | "
                     f"{_f(ci[0], 2)} - {_f(ci[1], 2)} | {_f(d.get('loso_calibrated_median_abs_dex'), 2)} | "
                     f"{_f(d.get('no_change_median_abs_dex'), 2)} |")
        L += ["", a1["calibration_reading"], ""]
    c0 = payload.get("stage_0c_as_first_run")
    if c0:
        L += ["## As first run (superseded by A1): stage 0c with the frozen generator's source-3 lanes", ""] + \
             _sigma_table(c0)
    a = payload.get("stage_0a")
    if a:
        L += ["## As first run (superseded by A1): stage 0a", "",
              f"{a['worlds']} of {a['worlds_planned']} worlds (B1 {a['rivals']['B1']}, B13 {a['rivals']['B13']}; "
              f"PR #108's points record at {a['points_record']['commit']}, sha256 {a['points_record']['sha256'][:12]}); "
              f"failed arm runs {a['failed_arm_runs']} of {a['arm_runs']}. Wall-clock {a['wall_hours']:.2f} h at 2 "
              f"processes ({a['core_seconds'] / 3600:.2f} core-h).", "",
              "Not computable at step 0: " + "; ".join(f"{k} ({v})" for k, v in a["not_computable"].items()), ""] + \
             _decision_table(a)
    c = payload.get("stage_0c")
    if c:
        L += ["## Stage 0c, B0-scored with A1's lane fix (the `--check` reference)", "",
              f"Units: {c['units']['pairs_source1']} pairs from `contrast_scores.json` (source 1), "
              f"{c['units']['pairs_source3']} answered within-study rows from the contrast benchmark (source 3; "
              f"{c['units']['source3_skipped']} censored or duplicate rows skipped); source 2: {c['units']['source2']}. "
              f"{c['units']['groups_out_of_sample']} of {c['units']['groups_all']} groups are out of sample.", "",
              f"Measurement: {c['measurement_rule']}.", ""] + _sigma_table(c) + [
              "Sensitivities (not deciding): with in-sample units, sigma_r = " + ", ".join(
                  f"{lane} {_f(e.get('sigma_r'))} ({e.get('studies', 0)} studies)"
                  for lane, e in c["sensitivity_in_sample_included"].items()) +
              "; robust (MAD) form: " + ", ".join(f"{lane} {_f(e.get('sigma_r'))}"
                                                  for lane, e in c["sensitivity_robust"].items()) + ".", "",
              f"As shipped today (iid omega): sigma_r = sqrt 2 x {OMEGA_DEX} = {SHIPPED_SIGMA_R:.3f} dex.", ""]
    return "\n".join(L) + "\n"


A1_DECLARATION = (
    "**Dated 2026-11-07, written before the re-scoring ran. Reason:** an independent check of step 0 found that "
    "stage 0c took sigma_r from residuals scored by the shipped engine, whose pH balance (B0) M1b part (i) falsified, "
    "while stage 0a's worlds use the carried balance rivals B1 / B13. The two pH-lever studies whose residuals set the "
    "sulfur sigma_r (xu2007, pH 5.6 / 7 / 8; Hofmann 1998, pH 3 / 7) were scored under B0. Their misses are therefore "
    "partly B0's pH failure, which the step-0 worlds have already removed. Counting them in sigma_r counts B0's failure "
    "twice. The frozen v3.1 text is not edited; this amendment is recorded here and in its revision log.",
    "",
    "**What A1 changes:**",
    "- **(a0) The check's concern, and what the prereg's lanes say.** The check asked that non-thiol species (Meynier "
    "NF; xu2007 furfural, 2-acetylthiazole) not set the thiol sigma_r. The prereg's lanes are ENGINE lanes "
    "(\"thiol = sulfur lane\", sec. 5.0.1), and the engine integrates every source-3 row on the sulfur lane, so (a) "
    "below moves no row. The concern is therefore met by a sensitivity, reported and not binding: \"rival-scored, "
    "thiol species only\". It takes the sulfur-lane groups of thiols and their disulfides only (MFT, FFT, their dimers, "
    "methanethiol). Changing the binding lane definition after the freeze would be a change of the frozen rule, which "
    "this amendment does not make.",
    "- **(a) A generator fix.** Source 3's rows are given their engine lane (from `declare_envelope` on each pot and "
    "species, no run) and their own exposure, where the generator as first run hard-coded \"sulfur\" and out of "
    "sample. The prereg's lanes are engine lanes (sec. 5.0.1).",
    "- **(b) Reporting fixes.** E is reported at sigma_r in {0.593, 0.42, 0.19}, which sec. 7.1 lists. A derived rho "
    "below -1 is printed as \"impossible: omega and sigma_r disagree (the engine over-responds to levers)\".",
    "- **(c) The rival-scored reading.** Every out-of-sample sulfur-lane pair is re-scored under each of PR #108's 23 "
    "rival points at 5e61fed0, with the shipped constants otherwise. Per group, the within-group sum of squares is "
    "averaged over the points, since each stage-0a world is one point. sigma_r then follows from sec. 5.0.1's "
    "estimator unchanged (G_min, floor, t-interval, LOSO).",
    "- **(d) Two rival-scored variants:**",
    "  - **\"rival-scored, all sulfur units\":** every out-of-sample sulfur pair re-scored;",
    "  - **\"rival-scored, xu2007 + Hofmann only\":** the check's variant, with the other studies kept B0-scored.",
    "- **(e) Which reading binds, chosen before any number:** \"rival-scored, all sulfur units\". A stage-0a world is a "
    "B1 or B13 engine, so every residual that sets its discrepancy must come from that engine. Re-scoring only the two "
    "flagged studies would leave the other studies' B0 residuals in. The other two readings are reported beside it.",
    "  - Stage 0a's eligibility is recomputed under each reading from the same 512 cached worlds, with no new world.",
    "  - Its point value and both LOSO ends bind (sec. 5.0.1): if they disagree, layer 1 is UNDECIDED.",
    "- **(f) A diagnostic** (A1 item 3): the slope of observed on predicted log ratios, and the leave-one-study-out "
    "calibrated contrast against no change. It is reported, and decides nothing.",
)


def amendment_a1(payload: Mapping[str, Any], workers: int, cache: Optional[Path]) -> Dict[str, Any]:
    """Amendment A1 (see A1_DECLARATION): the rival-scored sigma_r readings and stage 0a's eligibility under each,
    from the cached 512 worlds (no new world is run)."""
    if cache is None or not cache.exists():
        raise SystemExit("amendment A1 recomputes stage 0a from the cached worlds: pass --cache")
    records = sorted((json.loads(x) for x in cache.read_text().splitlines()), key=lambda r: r["world"])
    if len(records) != N_WORLDS:
        raise SystemExit(f"the cache holds {len(records)} worlds, not {N_WORLDS}")
    points_sha = payload["stage_0a"]["points_record"]["sha256"]
    pairs, arm_var, counts = collect_pairs()
    rescore = a1_rescore(pairs, workers)
    b0_fixed = stage_0c(pairs=pairs, arm_var=arm_var, counts=counts)
    two = a1_estimates(pairs, arm_var, counts, rescore, A1_FLAGGED_STUDIES)
    allr = a1_estimates(pairs, arm_var, counts, rescore, None)
    readings = {
        "b0_scored": {"label": "B0-scored (A1's lane fix only)", "sigma": b0_fixed},
        "rival_two_study": {"label": "rival-scored, xu2007 + Hofmann only", "sigma": two},
        "rival_all": {"label": "rival-scored, all sulfur units", "sigma": allr},
        "rival_all_thiol_only": {"label": "rival-scored, thiol species only (sensitivity)", "sigma": allr["thiol_only"]},
    }
    for r in readings.values():
        r["stage_0a"] = summarise_0a(records, r["sigma"], points_sha, 0.0, N_WORLDS)
    readings["rival_two_study"]["parity_b0_max_abs_dex"] = two["parity_b0_max_abs_dex"]
    readings["rival_all"]["parity_b0_max_abs_dex"] = allr["parity_b0_max_abs_dex"]
    readings["rival_all"]["parity_b0_n"] = allr["parity_b0_n"]
    # the lever-skill diagnostic: B0-scored, B0-scored without the two flagged studies, and rival-scored (mean p)
    rival_pairs = []
    for p in pairs:
        q = dict(p)
        k = "|".join([p["study"], p["compound"], p["a"], p["b"]])
        if k in allr["mean_rival_p"]:
            q["p"] = allr["mean_rival_p"][k]
            q["e"] = q["p"] - q["m"]
        rival_pairs.append(q)
    calib = {
        "B0-scored": calibration_diagnostics(pairs),
        "B0-scored, without xu2007 and Hofmann": calibration_diagnostics(pairs, drop=A1_FLAGGED_STUDIES),
        "rival-scored (mean over the 23 points)": calibration_diagnostics(rival_pairs),
    }
    return {"date": "2026-11-07", "declaration": list(A1_DECLARATION), "binding": "rival_all",
            "binding_label": readings["rival_all"]["label"],
            "rescore": {k: v for k, v in rescore.items() if k != "results"},
            "rescore_values": rescore["results"], "readings": readings, "calibration": calib,
            "calibration_reading": (
                "The check reported beta = 0.18 (90 % CI -0.06 to 0.42), 0.37 without the B0-pH studies, and a "
                "LOSO calibrated contrast of 0.63 vs 0.57 dex no-change (0.50 vs 0.49 without them). The numbers "
                "above are this generator's own, by the method stated in the table. Read together: out of sample, "
                "the engine shows no lever skill even after shrinkage. A diagnostic; it decides nothing.\n\n"
                "**What the re-scoring found.** Under B1 / B13 the two flagged studies still miss as much as under "
                "B0 (per-study sigma_r xu2007 2.78 -> 2.62, Hofmann 2.51 -> 2.58 dex). Their misses are not B0's "
                "pH-balance failure: xu2007 is a declared-pyrophosphate design, where the balance rivals move little. "
                "So the check's premise does not hold, and dropping the two studies (which would give ~0.5 dex) has "
                "no basis in the rule. The binding sulfur sigma_r is the rival-scored 1.23 dex.")}


def write(payload):
    artifact_io.write_artifact(payload, OUT, render=render)


def base_payload() -> Dict[str, Any]:
    return {"artifact": "decision_relevance_step0",
            "provenance": provenance.provenance_block(
                "decision_relevance_step0", generated_by="scripts/generators/build_decision_relevance_step0.py",
                inputs=[CONTRAST_SCORES, PANEL_SCORES, CONTRAST_BENCHMARK]),
            "prereg": PREREG, "freeze_commit": FREEZE_COMMIT}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--stage", choices=["0c", "0a", "all"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--worlds", type=int, default=N_WORLDS)
    ap.add_argument("--cache", type=Path, default=None, help="JSON-lines world cache (resumable)")
    ap.add_argument("--render", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--amendment", choices=["A1"], help="post-freeze amendment A1 (needs --cache with the 512 worlds)")
    args = ap.parse_args(argv)
    md = OUT.with_suffix(".md")
    if args.check:
        on_disk = json.loads(OUT.read_text())
        fresh = json.loads(json.dumps(stage_0c(), default=float))
        ok = json.loads(json.dumps(on_disk.get("stage_0c"), default=float)) == fresh
        ok0 = True
        if "stage_0c_as_first_run" in on_disk:
            ok0 = json.loads(json.dumps(on_disk["stage_0c_as_first_run"], default=float)) == json.loads(
                json.dumps(stage_0c(as_frozen=True), default=float))
        ok_md = md.read_text() == render(on_disk)
        print(f"stage 0c {'current' if ok else 'STALE'}; as-first-run 0c {'current' if ok0 else 'STALE'}; "
              f"md {'current' if ok_md else 'STALE'}")
        return 0 if ok and ok0 and ok_md else 1
    if args.render:
        md.write_text(render(json.loads(OUT.read_text())))
        return 0
    payload = json.loads(OUT.read_text()) if OUT.exists() else {}
    payload.update(base_payload())
    if args.stage in ("0c", "all"):
        payload["stage_0c"] = stage_0c()
        write(payload)
    if args.stage in ("0a", "all"):
        if "stage_0c" not in payload:
            raise SystemExit("stage 0a needs stage 0c first (sec. 5.6: 0c, then 0a)")
        payload["stage_0a"] = stage_0a(payload["stage_0c"], args.workers, args.worlds, args.cache)
        write(payload)
    if args.amendment == "A1":
        payload["amendment_A1"] = amendment_a1(payload, args.workers, args.cache)
        payload["stage_0c"] = stage_0c()
        payload["stage_0c_as_first_run"] = stage_0c(as_frozen=True)
        write(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
