#!/usr/bin/env python
"""
Step-0 addendum by lever class (STATUS_AND_PLAN sec. 0, step 1; 2026-11-10). A DIAGNOSTIC on spent data.

WHY. Step 0 of the decision-relevance pre-registration (``decision_relevance_step0.md``) reports ONE ratio error,
sigma_r, for the sulfur lane and reads it as "the lever table cannot show directions whatever the experiments".
An independent read (STATUS_AND_PLAN sec. 4.1, 2026-11-10; a read-only recomputation) found three things the record
does not report: the error scales with the predicted effect and is carried by two pH studies; sign accuracy; and
whether the shipped PARAMETER bands could account for the error ("structure, not parameters" was asserted, not
tested). This generator reproduces that read so nothing rests on an unreproduced recomputation.

WHAT IT READS. Exactly step 0's units: ``collect_pairs`` of ``build_decision_relevance_step0.py`` on the three
inputs whose sha256 the step-0 record pins (read from git history if the working copy has moved since). Nothing in
the frozen pre-registration or the step-0 record is edited; this is a new record beside them.

WHAT IT COMPUTES (all on the out-of-sample pairs; sulfur lane unless stated).
  A. Lever class of every pair, from its two arms' specs (pH, temperature, time, precursor identity, precursor
     dose, buffer, matrix; "mixed" when two or more differ). Per class: pairs, studies, sigma_r (step 0's
     estimator), RMS error of the engine, RMS of "no change", median |predicted|.
  B. Error against the predicted effect: the same statistics in bins of |predicted log ratio|.
  C. Sign accuracy (sign(predicted) == sign(measured); invariant to the arbitrary pair orientation), at two cuts
     declared here: |p| >= 0.3 dex; and |p| >= 0.05 with |m| >= 0.1. Exact binomial p against 0.5.
  D. Each of A - C with and without the two studies amendment A1 flagged (xu2007, Hofmann pH 3 / 7), and the
     slope of observed on predicted (step 0's own ``calibration_diagnostics``) both ways.
  E. The shipped parameter envelope on the residual pairs: every arm re-run under N_DRAWS draws of the shipped
     envelope sampler (``uncertainty.sample_draws``; B0 balance, as every scored pair was). Per draw, sigma_r over
     the same groups. Reported: the parameter spread of each pair's predicted ratio, its RMS against sigma_r, the
     central run's parity with the recorded values, and the distribution of sigma_r across draws.

DECLARED READING OF E (written before E was run). Let s_min be the 5th percentile of the per-draw sulfur sigma_r
and s_par the RMS over pairs of the across-draw SD of the predicted log ratio.
  * s_min <= 0.6 dex: worlds inside the shipped bands remove most of the error; "structure, not parameters" is
    REFUTED as stated (parameters, inside their bands, can account for it).
  * s_min >= 1.0 dex and s_par < 0.5 x sigma_r: no world inside the bands comes close; the error is not
    within-band parameter uncertainty. SUPPORTED, conditional on the bands being honest (S4).
  * otherwise: PARTLY (both contribute); the record says which pairs move.
A world's sigma_r is chosen post hoc from N_DRAWS, so s_min is optimistic; the 5th percentile, not the minimum,
is the reading. The bands themselves come from fits on overlapping data, so E bounds within-band variation only.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
if str(ROOT) not in sys.path:  # the generators' own convention (scripts/ is not a package)
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "decision_relevance_step0_addendum.json"
STEP0 = data_paths.VALIDATION_DIR / "decision_relevance_step0.json"
N_DRAWS, SEED = 64, 20261110
SIGN_CUT_P, SIGN_CUT_P_LOOSE, SIGN_CUT_M = 0.3, 0.05, 0.1
EFFECT_BINS = ((0.0, 0.05), (0.05, 0.3), (0.3, 1.0), (1.0, math.inf))
ARM_TIMEOUT_S = 300


def _step0_module():
    spec = importlib.util.spec_from_file_location(
        "build_decision_relevance_step0", ROOT / "scripts" / "generators" / "build_decision_relevance_step0.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S0 = _step0_module()


# ----------------------------------------------------------------------------------------------- pinned inputs
def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def pinned_bytes(rel: str, sha: str) -> bytes:
    """The bytes of ``rel`` whose sha256 is ``sha``: the working copy if it matches, else the newest git revision
    of the file that does."""
    p = ROOT / rel
    if p.exists() and _sha(p.read_bytes()) == sha:
        return p.read_bytes()
    revs = subprocess.run(["git", "log", "--all", "--format=%H", "--", rel], cwd=ROOT, capture_output=True,
                          text=True, check=True).stdout.split()
    for rev in revs:
        b = subprocess.run(["git", "show", f"{rev}:{rel}"], cwd=ROOT, capture_output=True).stdout
        if b and _sha(b) == sha:
            return b
    raise SystemExit(f"no revision of {rel} has sha256 {sha}")


def step0_pairs() -> Tuple[List[Dict[str, Any]], Dict[str, float], Dict[str, int], List[Dict[str, str]]]:
    """Step 0's units, rebuilt from the exact inputs its record pins."""
    rec = json.loads(STEP0.read_text())
    pins = rec["provenance"]["inputs"]
    attr = {"results/validation/contrast_scores.json": "CONTRAST_SCORES",
            "results/validation/core_panel_scores.json": "PANEL_SCORES",
            "results/validation/kinetic_core_contrast_benchmark_score_v1.json": "CONTRAST_BENCHMARK"}
    import tempfile
    tmp = Path(tempfile.mkdtemp(prefix="step0_addendum_inputs_"))
    try:
        for pin in pins:
            b = pinned_bytes(pin["path"], pin["sha256"])
            f = tmp / Path(pin["path"]).name
            f.write_bytes(b)
            setattr(S0, attr[pin["path"]], f)
        pairs, arm_var, counts = S0.collect_pairs()
    finally:
        for f in tmp.glob("*"):
            f.unlink()
        tmp.rmdir()
    return pairs, arm_var, counts, pins


# ------------------------------------------------------------------------------------------------- lever class
def _thermal(spec) -> Tuple[float, float]:
    segs = spec.process.thermal.segments
    return float(max(t for _, t in segs)), float(sum(d for d, _ in segs))


def lever_class(spec_a, spec_b) -> Tuple[str, List[str]]:
    diff = []
    if abs(float(spec_a.process.ph) - float(spec_b.process.ph)) > 1e-6:
        diff.append("pH")
    ta, tb = _thermal(spec_a), _thermal(spec_b)
    if abs(ta[0] - tb[0]) > 1e-6:
        diff.append("temperature")
    if abs(ta[1] - tb[1]) > 1e-6:
        diff.append("time")
    ka = {k for k, v in spec_a.precursors.items() if float(v) > 0}
    kb = {k for k, v in spec_b.precursors.items() if float(v) > 0}
    if ka != kb:
        diff.append("precursor identity")
    elif any(abs(float(spec_a.precursors[k]) - float(spec_b.precursors[k])) > 1e-9 for k in ka):
        diff.append("precursor dose")
    if repr(spec_a.process.buffer) != repr(spec_b.process.buffer):
        diff.append("buffer")
    if spec_a.process.matrix != spec_b.process.matrix:
        diff.append("matrix")
    if not diff:
        return "none detected", diff
    return (diff[0] if len(diff) == 1 else "mixed"), diff


def classify(pairs: Sequence[Dict[str, Any]]) -> None:
    for p in pairs:
        try:
            sa, _ = S0.unit_spec(p["source"], p["study"], p["a"], p["compound"])
            sb, _ = S0.unit_spec(p["source"], p["study"], p["b"], p["compound"])
            p["lever"], p["lever_fields"] = lever_class(sa, sb) if sa is not None and sb is not None \
                else ("spec unavailable", [])
        except Exception as exc:  # noqa: BLE001  (recorded, counted)
            p["lever"], p["lever_fields"] = "spec unavailable", [type(exc).__name__]


# -------------------------------------------------------------------------------------------------- statistics
def _rms(x) -> Optional[float]:
    x = np.asarray(list(x), float)
    return float(math.sqrt(np.mean(x ** 2))) if x.size else None


def summary(sel: Sequence[Mapping[str, Any]], arm_var: Mapping[str, float]) -> Dict[str, Any]:
    if not sel:
        return {"pairs": 0}
    est = S0.sigma_r_estimate(S0.groups_from_pairs(sel, arm_var))
    return {"pairs": len(sel), "studies": len({p["study"] for p in sel}),
            "study_names": sorted({p["study"] for p in sel}),
            "sigma_r": est.get("sigma_r"), "sigma_r_df": est.get("df"),
            "rms_engine_error": _rms(p["e"] for p in sel), "rms_no_change": _rms(p["m"] for p in sel),
            "median_abs_predicted": float(np.median([abs(p["p"]) for p in sel])),
            "median_abs_measured": float(np.median([abs(p["m"]) for p in sel]))}


def sign_accuracy(sel: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    from scipy import stats

    def one(rows):
        n = len(rows)
        k = sum(1 for p in rows if np.sign(p["p"]) == np.sign(p["m"]))
        return {"n": n, "right": k, "share": (k / n) if n else None,
                "binom_p_two_sided": float(stats.binomtest(k, n, 0.5).pvalue) if n else None}
    return {f"|p| >= {SIGN_CUT_P}": one([p for p in sel if abs(p["p"]) >= SIGN_CUT_P]),
            f"|p| >= {SIGN_CUT_P_LOOSE} and |m| >= {SIGN_CUT_M}":
                one([p for p in sel if abs(p["p"]) >= SIGN_CUT_P_LOOSE and abs(p["m"]) >= SIGN_CUT_M])}


def by_key(sel, key, arm_var) -> Dict[str, Any]:
    buckets: Dict[str, List] = defaultdict(list)
    for p in sel:
        buckets[key(p)].append(p)
    return {k: summary(v, arm_var) for k, v in sorted(buckets.items())}


def effect_bin(p) -> str:
    a = abs(p["p"])
    for lo, hi in EFFECT_BINS:
        if lo <= a < hi:
            return f"{lo:g} - {hi:g} dex"
    return "?"


def section_a_to_d(pairs, arm_var) -> Dict[str, Any]:
    sul = [p for p in pairs if p["oos"] and p["lane"] == "sulfur"]
    flagged = set(S0.A1_FLAGGED_STUDIES)
    out = {}
    for name, sel in (("all", sul), ("without the two A1-flagged pH studies", [p for p in sul
                                                                               if p["study"] not in flagged])):
        out[name] = {"overall": summary(sel, arm_var),
                     "by_lever_class": by_key(sel, lambda p: p["lever"], arm_var),
                     "by_predicted_effect": by_key(sel, effect_bin, arm_var),
                     "sign_accuracy": sign_accuracy(sel)}
    out["slope"] = {"all": S0.calibration_diagnostics(pairs),
                    "without the two A1-flagged pH studies": S0.calibration_diagnostics(pairs, drop=sorted(flagged))}
    out["other_lanes_out_of_sample"] = by_key([p for p in pairs if p["oos"] and p["lane"] != "sulfur"],
                                              lambda p: p["lane"], arm_var)
    return out


# ------------------------------------------------------------------------------------------- E. the envelope
_W: Dict[str, Any] = {}


def _init(draws):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _W["draws"] = draws


def _run_world(task: Mapping[str, Any]) -> Dict[str, Any]:
    from src.kinetic_core import engine

    draw = None if task["index"] is None else _W["draws"][task["index"]].core
    out: Dict[str, Optional[float]] = {}
    t0 = time.perf_counter()
    signal.signal(signal.SIGALRM, S0._alarm)
    for u in task["units"]:
        key = "|".join([str(u["source"]), u["study"], u["arm"], u["compound"]])
        try:
            spec, tg = S0.unit_spec(u["source"], u["study"], u["arm"], u["compound"])
            if spec is None or tg is None:
                out[key] = None
                continue
            signal.alarm(ARM_TIMEOUT_S)
            kw = {} if draw is None else {"draw": draw, "size_declared_bands": False}
            run = engine.predict(spec, [tg], **kw)
            signal.alarm(0)
            mk = run.declaration.mapped_targets.get(tg) if run.answered else None
            out[key] = None if mk is None else float(run.species_mmol_per_l.get(mk, 0.0))
        except S0._Timeout:
            out[key] = None
        except Exception:  # noqa: BLE001  (a missing value, counted)
            signal.alarm(0)
            out[key] = None
    return {"index": task["index"], "values": out, "seconds": time.perf_counter() - t0}


def section_e(pairs, arm_var, workers: int, n_draws: int) -> Dict[str, Any]:
    import multiprocessing as mp
    from src.kinetic_core import uncertainty

    sel = [p for p in pairs if p["oos"] and p["lane"] == "sulfur"]
    units = sorted({(p["source"], p["study"], arm, p["compound"]) for p in sel for arm in (p["a"], p["b"])})
    unit_dicts = [{"source": s, "study": st, "arm": a, "compound": c} for s, st, a, c in units]
    draws = uncertainty.sample_draws(n_draws, seed=SEED)
    tasks = [{"index": None, "units": unit_dicts}] + [{"index": i, "units": unit_dicts} for i in range(n_draws)]
    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(draws,)) as pool:
        results = list(pool.imap_unordered(_run_world, tasks))
    central = next(r for r in results if r["index"] is None)
    worlds = sorted((r for r in results if r["index"] is not None), key=lambda r: r["index"])

    parity = [abs(x - p["p"]) for p in sel if (x := S0._pair_p(p, central["values"])) is not None]
    per_pair, per_world_sigma, missing = [], [], 0
    pw = {id(p): [] for p in sel}
    for w in worlds:
        new = []
        for p in pairs:
            q = dict(p)
            if p in sel:
                x = S0._pair_p(p, w["values"])
                if x is None:
                    missing += 1
                    continue
                pw[id(p)].append(x)
                q["p"], q["e"] = x, x - p["m"]
            new.append(q)
        est = S0.sigma_r_estimate([g for g in S0.groups_from_pairs(new, arm_var)
                                   if g["out_of_sample"] and g["lane"] == "sulfur"])
        per_world_sigma.append(est.get("sigma_r"))
    for p in sel:
        xs = np.array(pw[id(p)])
        if xs.size < 2:
            continue
        q5, q95 = np.percentile(xs, [5, 95])
        per_pair.append({"study": p["study"], "compound": p["compound"], "a": p["a"], "b": p["b"],
                         "lever": p["lever"], "measured": p["m"], "recorded_predicted": p["p"],
                         "draw_mean": float(xs.mean()), "draw_sd": float(xs.std(ddof=1)),
                         "draw_q05": float(q5), "draw_q95": float(q95),
                         "measured_inside_90": bool(q5 <= p["m"] <= q95)})
    sig = np.array([s for s in per_world_sigma if s is not None and math.isfinite(s)])
    point = S0.sigma_r_estimate([g for g in S0.groups_from_pairs(sel, arm_var)])["sigma_r"]
    s_par = _rms(r["draw_sd"] for r in per_pair)
    s_min = float(np.percentile(sig, 5)) if sig.size else None
    if s_min is None:
        reading = "NOT COMPUTED"
    elif s_min <= 0.6:
        reading = "REFUTED as stated: worlds inside the shipped bands remove most of the error"
    elif s_min >= 1.0 and s_par is not None and s_par < 0.5 * point:
        reading = ("SUPPORTED, conditional on honest bands: no world inside the shipped bands comes close, and "
                   "within-band parameter spread is small next to sigma_r")
    else:
        reading = "PARTLY: within-band parameters move the error but do not remove it"
    return {"n_draws": n_draws, "seed": SEED, "sampler": "src.kinetic_core.uncertainty.sample_draws (shipped priors)",
            "units": len(units), "pairs": len(sel), "wall_hours": (time.time() - t0) / 3600.0,
            "core_seconds": float(sum(r["seconds"] for r in results)),
            "central_parity_max_abs_dex": max(parity) if parity else None, "central_parity_n": len(parity),
            "pair_world_values_missing": missing,
            "sigma_r_point": point,
            "sigma_r_by_world": {"q05": s_min, "median": float(np.median(sig)) if sig.size else None,
                                 "q95": float(np.percentile(sig, 95)) if sig.size else None,
                                 "min": float(sig.min()) if sig.size else None, "worlds": int(sig.size)},
            "rms_parameter_sd_of_predicted_ratio": s_par,
            "coverage_of_measured_by_parameter_90": (float(np.mean([r["measured_inside_90"] for r in per_pair]))
                                                     if per_pair else None),
            "declared_reading": reading, "per_pair": per_pair}


# ------------------------------------------------------------------------------------------------------ record
def _f(x, nd=2):
    return "-" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}"


def _table(rows: Mapping[str, Mapping[str, Any]], head: str) -> List[str]:
    L = [f"| {head} | pairs | studies | sigma_r | RMS engine error | RMS no change | median abs predicted | "
         "median abs measured |", "|---|---|---|---|---|---|---|---|"]
    for k, s in rows.items():
        if not s.get("pairs"):
            continue
        L.append(f"| {k} | {s['pairs']} | {s['studies']} | {_f(s['sigma_r'])} | {_f(s['rms_engine_error'])} | "
                 f"{_f(s['rms_no_change'])} | {_f(s['median_abs_predicted'])} | {_f(s['median_abs_measured'])} |")
    return L


def render(payload: Mapping[str, Any]) -> str:
    ad = payload["a_to_d"]
    L = ["# Decision relevance, step-0 addendum: the ratio error by lever class",
         "",
         "*Generated by `scripts/generators/build_decision_relevance_step0_addendum.py` (STATUS_AND_PLAN sec. 0, "
         "step 1). A diagnostic on spent data, beside the step-0 record "
         "(`decision_relevance_step0.md`), which it does not edit. Its units are exactly step 0's, rebuilt from the "
         "three inputs that record pins by sha256. It credits nothing and moves no shipped output.*", ""]
    if payload.get("e"):
        e = payload["e"]
        L += ["## The headline", "",
              f"> **Parameters within the shipped bands (E): {e['declared_reading']}.** Sulfur sigma_r "
              f"{_f(e['sigma_r_point'])} dex at the shipped constants; across {e['sigma_r_by_world']['worlds']} "
              f"envelope worlds 5 / 50 / 95 %: {_f(e['sigma_r_by_world']['q05'])} / "
              f"{_f(e['sigma_r_by_world']['median'])} / {_f(e['sigma_r_by_world']['q95'])} dex; RMS within-band SD "
              f"of a predicted ratio {_f(e['rms_parameter_sd_of_predicted_ratio'])} dex.", ""]
    for name, block in ((k, v) for k, v in ad.items() if k not in ("slope", "other_lanes_out_of_sample")):
        L += [f"## Sulfur lane, out of sample: {name}", "",
              *(_table({"all pairs": block["overall"]}, "set")), "",
              "### By lever class (A)", "", *_table(block["by_lever_class"], "lever class"), "",
              "### By predicted effect (B)", "", *_table(block["by_predicted_effect"], "abs predicted log ratio"), "",
              "### Sign accuracy (C)", "", "| cut | right / n | share | binomial p vs 0.5 |", "|---|---|---|---|"]
        for cut, s in block["sign_accuracy"].items():
            L.append(f"| {cut} | {s['right']} / {s['n']} | {_f(s['share'])} | {_f(s['binom_p_two_sided'], 3)} |")
        L.append("")
    L += ["## Slope of observed on predicted (D; step 0's own estimator)", "",
          "| set | pairs | studies | beta | 90 % CI | LOSO calibrated median abs | no-change median abs |",
          "|---|---|---|---|---|---|---|"]
    for k, s in ad["slope"].items():
        ci = s.get("beta_ci90") or [None, None]
        L.append(f"| {k} | {s.get('pairs')} | {s.get('studies')} | {_f(s.get('beta'))} | {_f(ci[0])} - {_f(ci[1])} "
                 f"| {_f(s.get('loso_calibrated_median_abs_dex'))} | {_f(s.get('no_change_median_abs_dex'))} |")
    L += ["", "## Other lanes, out of sample", "", *_table(ad["other_lanes_out_of_sample"], "lane"), ""]
    if payload.get("e"):
        e = payload["e"]
        L += ["## E. The shipped parameter envelope on the residual pairs", "",
              f"- {e['n_draws']} draws of the shipped priors (seed {e['seed']}), {e['units']} arms, {e['pairs']} "
              f"pairs; {_f(e['wall_hours'])} h wall, {_f(e['core_seconds'] / 3600.0)} core-h.",
              f"- Parity of the central run with the recorded predictions: max {_f(e['central_parity_max_abs_dex'], 4)} "
              f"dex over {e['central_parity_n']} pairs; pair-world values missing: {e['pair_world_values_missing']}.",
              f"- Share of measured ratios inside the parameter-only 90 % band: "
              f"{_f(e['coverage_of_measured_by_parameter_90'])} (measurement noise not added).",
              "- Share inside the parameter-only 90 % band, by lever class: " + "; ".join(
                  f"{lv} {sum(r['measured_inside_90'] for r in e['per_pair'] if r['lever'] == lv)} / "
                  f"{sum(1 for r in e['per_pair'] if r['lever'] == lv)}"
                  for lv in sorted({r['lever'] for r in e['per_pair']})) + ".",
              "- Declared reading (written before the run, in the generator's docstring): 5th-percentile world "
              "sigma_r <= 0.6 dex refutes \"structure, not parameters\" as stated; >= 1.0 dex with within-band "
              "spread below half of sigma_r supports it, conditional on honest bands; otherwise PARTLY.", "",
              "| study | compound | lever | measured | recorded predicted | draw 5 / 50 / 95 % | inside |",
              "|---|---|---|---|---|---|---|"]
        for r in sorted(e["per_pair"], key=lambda r: (r["study"], r["compound"])):
            L.append(f"| {r['study']} | {r['compound']} | {r['lever']} | {_f(r['measured'])} | "
                     f"{_f(r['recorded_predicted'])} | {_f(r['draw_q05'])} / {_f(r['draw_mean'])} / "
                     f"{_f(r['draw_q95'])} | {'yes' if r['measured_inside_90'] else 'no'} |")
        L.append("")
    L += ["## What this does not do", "",
          "- It does not re-run step 0's oracle or change its verdict under the frozen rule.",
          "- E varies parameters inside the SHIPPED bands only; those bands come from fits on overlapping data.",
          "- Lever classes are read from the specs; a pair whose arms differ in two fields is \"mixed\".", ""]
    return "\n".join(L)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--draws", type=int, default=N_DRAWS)
    ap.add_argument("--skip-envelope", action="store_true", help="sections A - D only (fast)")
    ap.add_argument("--render", action="store_true", help="re-render the markdown from the JSON on disk")
    args = ap.parse_args(argv)
    if args.render:
        OUT.with_suffix(".md").write_text(render(json.loads(OUT.read_text())))
        return 0

    pairs, arm_var, counts, pins = step0_pairs()
    classify(pairs)
    payload: Dict[str, Any] = {"artifact": "decision_relevance_step0_addendum", "step0_record": str(STEP0.relative_to(ROOT)),
                               "units": counts, "a_to_d": section_a_to_d(pairs, arm_var),
                               "pair_levers": [{k: p[k] for k in ("source", "study", "compound", "lane", "a", "b",
                                                                   "p", "m", "oos", "lever", "lever_fields")}
                                               for p in pairs]}
    if not args.skip_envelope:
        payload["e"] = section_e(pairs, arm_var, args.workers, args.draws)
    payload["provenance"] = provenance.provenance_block(
        "decision_relevance_step0_addendum",
        generated_by="scripts/generators/build_decision_relevance_step0_addendum.py",
        inputs=[STEP0, ROOT / "scripts" / "generators" / "build_decision_relevance_step0.py"])
    payload["step0_pinned_inputs"] = pins
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {OUT.relative_to(ROOT)} and .md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
