#!/usr/bin/env python
"""
score_deoxyosone_signatures.py -- B7a's scale-free signatures of the deoxyosone partition (2026-10-31).

results/validation/kinetic_core_dpo_partition_prereg.md secs. 3 and 7. A calibration scale shared by the 3- and
1-deoxypentosone of one figure (one method, OPD quinoxalines) cancels in their RATIO, so the ratio's shape is a test a
fitted scale cannot absorb. Declared before any B7a arm was read:
  (i)   the time slope of log10(TDP / DPO) against log10(t), per pot, |predicted - measured| <= 0.3;
  (ii)  the pH slope of log10(TDP / DPO), Zhai 2021's two ladders (TTCA; TTCA + xylose; 100 C), by segment
        (5.5 - 6, 6 - 7, 7 - 8), each within +-0.15 dex per pH unit;
  (iii) cysteine's double ratio, Zhai 2019: (TDP / DPO) with cysteine over without -- REPORTED, not decisive (CF-071);
  (vi)  TDP's decline after its measured peak, Zhai 2023 JAFC (TTCA, 120 / 140 C): the apparent first-order rate,
        log10(predicted / measured) within +-0.3. Decisive on those two pots; the same rate on every other pot with a
        measured TDP series (JAFC 100 C, the + Cys pots, Zhai 2023 FC 140 C) and DPO's own post-peak rate beside it
        are REPORTED, not decisive (added after B7a's design review, 2026-10-31, before any arm was read).
Also reported: the engine's in-situ pH at 60 min in each Zhai 2021 pot (the pH the split law reads).
Each is computed from the two-reader table (measured) and the engine at a fit's optimum (predicted), the fit's arm
and fixed barriers read from its record.

    python scripts/generators/score_deoxyosone_signatures.py --fit s042_r3en_fixed
Writes results/validation/kinetic_core_deoxyosone_signatures_<fit>.{json,md}.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402
import score_mechanism_signatures as sc  # noqa: E402

TABLE = "deoxyosone_submodel"
TIME_SLOPE_POTS = ("zhai2019_ttca_120c", "zhai2019_ttca_cys_120c", "zhai2021_ttca_ph55_100c_dx",
                   "zhai2021_ttca_ph70_100c_dx", "zhai2021_ttca_ph80_100c_dx", "zhai2023fc_dx_ttca_xyl_140c",
                   "zhai2023jafc_ttca_120C_dx", "zhai2023jafc_ttca_140C_dx")
PH_LADDERS = {"TTCA": ("zhai2021_ttca_ph55_100c_dx", "zhai2021_ttca_ph60_100c_dx", "zhai2021_ttca_ph70_100c_dx",
                       "zhai2021_ttca_ph80_100c_dx"),
              "TTCA + xylose": ("zhai2021_ttca_xyl_ph55_100c_dx", "zhai2021_ttca_xyl_ph60_100c_dx",
                                "zhai2021_ttca_xyl_ph70_100c_dx", "zhai2021_ttca_xyl_ph80_100c_dx")}
PH_OF = {"ph55": 5.5, "ph60": 6.0, "ph70": 7.0, "ph80": 8.0}
DOUBLE_RATIO = ("zhai2019_ttca_cys_120c", "zhai2019_ttca_120c")
DECLINE_POTS = ("zhai2023jafc_ttca_120C_dx", "zhai2023jafc_ttca_140C_dx")
#: reported beside the decisive pair (design review, 2026-10-31): not counted in the (vi) verdict
REPORTED_DECLINE_POTS = ("zhai2023jafc_ttca_100C_dx", "zhai2023jafc_ttca_cys_100C_dx", "zhai2023jafc_ttca_cys_120C_dx",
                         "zhai2023jafc_ttca_cys_140C_dx", "zhai2023fc_dx_ttca_xyl_140c")
PH_REPORT_MIN = 60.0
TOL_TIME, TOL_PH, TOL_DECLINE = 0.3, 0.15, 0.3


def _slope(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.polyfit(x, y, 1)[0]) if len(x) >= 2 else float("nan")


def score(trial: Dict[str, Any], strecker_share: float, plans, rows, strict: bool = True):
    """Every B7a signature at one trial (the arm already configured): (results, pass counts); with strict=False a
    pot the engine cannot run makes the draw (None, None) instead of stopping."""
    ph_at_60: Dict[str, float] = {}
    meas: Dict[Tuple[str, str, float], float] = {}
    pred: Dict[Tuple[str, str, float], float] = {}
    for pot, job in plans.items():
        j = dict(job, **trial, strecker_share=float(strecker_share))
        res = joint.run_task(j)
        if "failed" in res:
            if strict:
                raise SystemExit(f"{pot}: {res['failed']}")
            return None, None
        if pot.startswith("zhai2021") and res.get("ph_in_situ") and max(j["times"]) >= PH_REPORT_MIN:
            ph_at_60[pot] = float(np.interp(PH_REPORT_MIN, j["times"], res["ph_in_situ"]))
        for r in rows:
            if r.pot == pot and r.species in ("TDP", "DPO"):
                meas[(pot, r.species, r.t_min)] = r.value
                pred[(pot, r.species, r.t_min)] = probe.predict(r, j, res)

    def ratio(src, pot):
        times = sorted(t for (p, s, t) in src if p == pot and s == "TDP" and (p, "DPO", t) in src)
        return times, [math.log10(src[(pot, "TDP", t)] / src[(pot, "DPO", t)])
                       if src[(pot, "TDP", t)] > 0 and src[(pot, "DPO", t)] > 0 else float("nan") for t in times]

    out: Dict[str, Any] = {"time_slope": [], "ph_segments": [], "double_ratio": None, "tdp_decline": []}
    for pot in TIME_SLOPE_POTS:
        tm, ym = ratio(meas, pot)
        tp, yp = ratio(pred, pot)
        sm = _slope(np.log10(tm), ym)
        sp = _slope(np.log10(tp), yp)
        out["time_slope"].append({"pot": pot, "measured": sm, "predicted": sp,
                                  "pass": bool(math.isfinite(sp) and abs(sp - sm) <= TOL_TIME)})
    for ladder, ladder_pots in PH_LADDERS.items():
        def level(src, pot):
            _, y = ratio(src, pot)
            return float(np.nanmean(y))
        phs = [PH_OF[p.split("_")[-3]] for p in ladder_pots]
        lm = [level(meas, p) for p in ladder_pots]
        lp = [level(pred, p) for p in ladder_pots]
        for a in range(len(phs) - 1):
            dm = (lm[a + 1] - lm[a]) / (phs[a + 1] - phs[a])
            dp = (lp[a + 1] - lp[a]) / (phs[a + 1] - phs[a])
            out["ph_segments"].append({"ladder": ladder, "segment": f"{phs[a]:g} - {phs[a + 1]:g}", "measured": dm,
                                       "predicted": dp, "pass": bool(math.isfinite(dp) and abs(dp - dm) <= TOL_PH)})
    with_cys, without = DOUBLE_RATIO
    dm = float(np.nanmean(ratio(meas, with_cys)[1]) - np.nanmean(ratio(meas, without)[1]))
    dp = float(np.nanmean(ratio(pred, with_cys)[1]) - np.nanmean(ratio(pred, without)[1]))
    out["double_ratio"] = {"measured": 10 ** dm, "predicted": 10 ** dp, "decisive": False}
    def decline(pot, species):
        """The apparent first-order rate after the species' MEASURED peak, measured and predicted at the same times."""
        series = sorted((t, v) for (p, sp, t), v in meas.items() if p == pot and sp == species)
        if not series:
            return {"note": f"no measured {species}"}
        peak = max(range(len(series)), key=lambda i: series[i][1])
        after = series[peak:]
        if len(after) < 2:
            return {"note": "no measured decline after the peak"}
        km = -_slope([t for t, _ in after], [math.log(v) for _, v in after])
        kp = -_slope([t for t, _ in after], [math.log(max(pred[(pot, species, t)], 1e-300)) for t, _ in after])
        return {"peak_min": after[0][0], "measured_per_min": km, "predicted_per_min": kp,
                "within_tolerance": bool(km > 0 and kp > 0 and abs(math.log10(kp / km)) <= TOL_DECLINE)}

    for pot in DECLINE_POTS + REPORTED_DECLINE_POTS:
        tdp = decline(pot, "TDP")
        entry = {"pot": pot, "decisive": pot in DECLINE_POTS, **tdp, "dpo": decline(pot, "DPO")}
        if "within_tolerance" in tdp:
            entry["pass"] = tdp["within_tolerance"]
        out["tdp_decline"].append(entry)
    out["ph_in_situ_at_60_min"] = ph_at_60

    passed = {"time_slope": sum(s["pass"] for s in out["time_slope"]),
              "ph_segments": sum(s["pass"] for s in out["ph_segments"]),
              "tdp_decline": sum(s.get("pass", False) for s in out["tdp_decline"] if s["decisive"])}
    return out, passed


def chance(arm: str, n: int, seed: int, strecker_share: float = 0.042) -> int:
    """The chance baseline of these signatures (STATUS_AND_PLAN rule d; B7a review S2): pass counts when every
    coordinate of the arm is drawn uniformly within its declared box -- what a constant set chosen without the data
    passes. Written before any B7a arm is credited."""
    from src.kinetic_core import parallel
    sc._configure(arm, strecker_share)
    keys, _kinds, lower, upper, _start = joint.coordinates(False, joint.ARMS[arm])
    rng = np.random.default_rng(seed)
    draws = [lower + (upper - lower) * rng.random(len(keys)) for _ in range(n)]
    pots, rows = fit_rows.load(TABLE)
    plans = {j["pot"]: j for j in probe.plan(pots, rows) if j["pot"] in _wanted()}

    def one(x):
        try:
            _, passed = score(joint.trial(list(x), keys), strecker_share, plans, rows, strict=False)
            return passed if passed is not None else {"error": "a pot failed"}
        except Exception as exc:  # noqa: BLE001 - counted as a failed draw
            return {"error": repr(exc)[:200]}

    results = parallel.map_ordered(one, draws)
    ok = [r for r in results if "error" not in r]
    dist = {fam: {str(k): sum(1 for r in ok if r[fam] == k) for k in sorted({r[fam] for r in ok})}
            for fam in ("time_slope", "ph_segments", "tdp_decline")}
    thresholds = {"time_slope": 5, "ph_segments": 4, "tdp_decline": 2}
    rate = {fam: sum(1 for r in ok if r[fam] >= t) / max(1, len(ok)) for fam, t in thresholds.items()}
    stem = data_paths.VALIDATION_DIR / f"kinetic_core_deoxyosone_signatures_chance_{arm.lower()}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/score_deoxyosone_signatures.py",
                   inputs=[ROOT / "src/kinetic_core/sulfur.py", fit_rows.FIT_ROWS_DIR / f"{TABLE}.yml"]),
               "arm": arm, "draws": n, "seed": seed, "runnable_draws": len(ok), "pass_count_distribution": dist,
               "chance_rate_of_the_decisive_threshold": rate, "thresholds": thresholds,
               "errors_sample": [r["error"] for r in results if "error" in r][:5]}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    lines = [f"# Deoxyosone partition signatures: chance baseline, arm {arm}", "",
             f"*Generated by `scripts/generators/score_deoxyosone_signatures.py --chance {n} --arm {arm}` (seed {seed}); "
             f"every coordinate uniform within the arm's box. {len(ok)} of {n} draws ran.*", "",
             "| family | decisive threshold | chance rate of reaching it | passes per draw |", "|---|---|---|---|"]
    lines += [f"| {fam} | >= {t} | {rate[fam]:.2f} | {dist[fam]} |" for fam, t in thresholds.items()]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: {rate}")
    return 0


def _wanted():
    return (set(TIME_SLOPE_POTS) | {p for ps in PH_LADDERS.values() for p in ps} | set(DOUBLE_RATIO) | set(DECLINE_POTS)
            | set(REPORTED_DECLINE_POTS))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", help="the joint fit's stem suffix (s042_r3en_fixed, ...)")
    parser.add_argument("--chance", type=int, metavar="N", help="the chance baseline from N uniform draws (with --arm)")
    parser.add_argument("--arm", choices=tuple(joint.ARMS), help="the arm whose box the chance baseline draws from")
    parser.add_argument("--seed", type=int, default=20261031)
    args = parser.parse_args(argv)
    if args.chance:
        if not args.arm:
            raise SystemExit("--chance needs --arm")
        return chance(args.arm, args.chance, args.seed)
    if not args.fit:
        raise SystemExit("--fit or --chance")
    fit_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json"
    fit = json.loads(fit_path.read_text())
    theta, arm = fit["fit"]["theta"], fit.get("arm") or sc._arm_of(args.fit)
    sc._configure(arm, float(fit["strecker_share"]))
    if fit.get("ea_static") is not None:
        joint.EA_STATIC = dict(fit["ea_static"])
    if fit.get("tdp_exit") is not None:    # B7a R-T's constants, from the fit's record (before B7a: from the arm)
        joint.TDP_EXIT_STATIC = dict(fit["tdp_exit"])
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    pots, rows = fit_rows.load(TABLE)
    wanted = (set(TIME_SLOPE_POTS) | {p for ps in PH_LADDERS.values() for p in ps} | set(DOUBLE_RATIO) | set(DECLINE_POTS)
              | set(REPORTED_DECLINE_POTS))
    plans = {j["pot"]: j for j in probe.plan(pots, rows) if j["pot"] in wanted}
    out, passed = score(trial, float(fit["strecker_share"]), plans, rows)
    stem = data_paths.VALIDATION_DIR / f"kinetic_core_deoxyosone_signatures_{args.fit}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/score_deoxyosone_signatures.py",
                   inputs=[fit_path, ROOT / "src/kinetic_core/sulfur.py", fit_rows.FIT_ROWS_DIR / f"{TABLE}.yml"]),
               "fit": data_paths.rel(fit_path), "arm": arm, "within_tau": fit.get("within_tau"),
               "fit_converged": not fit["fit"].get("budget_exhausted", False), "passed": passed, **out}
    stem.with_suffix(".json").write_text(provenance.dumps_compact(payload) + "\n")
    lines = [f"# Deoxyosone partition signatures (B7a): {args.fit} (arm {arm}, within_tau {fit.get('within_tau')})", "",
             "*Generated by `scripts/generators/score_deoxyosone_signatures.py`; declared in "
             "`results/validation/kinetic_core_dpo_partition_prereg.md` secs. 3 and 7 before any B7a arm was read.*", "",
             f"**Time slope {passed['time_slope']} / {len(out['time_slope'])} pots; pH segments "
             f"{passed['ph_segments']} / {len(out['ph_segments'])}; TDP decline {passed['tdp_decline']} / "
             f"{len(DECLINE_POTS)}.**" + ("" if payload["fit_converged"] else
                                                 " The fit stopped on its budget: no verdict (rule a)."), "",
             f"## (i) d log10(TDP / DPO) / d log10(t), tolerance {TOL_TIME}", "", "| pot | measured | predicted | pass |",
             "|---|---|---|---|"]
    lines += [f"| {s['pot']} | {s['measured']:+.2f} | {s['predicted']:+.2f} | {'yes' if s['pass'] else 'NO'} |"
              for s in out["time_slope"]]
    lines += ["", f"## (ii) d log10(TDP / DPO) / d pH by segment (Zhai 2021, 100 C), tolerance {TOL_PH}", "",
              "| ladder | segment | measured | predicted | pass |", "|---|---|---|---|---|"]
    lines += [f"| {s['ladder']} | {s['segment']} | {s['measured']:+.2f} | {s['predicted']:+.2f} | "
              f"{'yes' if s['pass'] else 'NO'} |" for s in out["ph_segments"]]
    dr = out["double_ratio"]
    lines += ["", "## (iii) cysteine's double ratio (Zhai 2019; reported, not decisive -- CF-071)", "",
              f"measured {dr['measured']:.2f}, predicted {dr['predicted']:.2f}", "",
              f"## (vi) TDP's apparent first-order decline after its measured peak, tolerance {TOL_DECLINE} dex "
              "(decisive: Zhai 2023 JAFC TTCA 120 / 140 C; the rest reported), with DPO's own beside it", "",
              "| pot | role | TDP from (min) | TDP measured /min | TDP predicted /min | pass | DPO from (min) | "
              "DPO measured /min | DPO predicted /min | DPO within |", "|---|---|---|---|---|---|---|---|---|---|"]

    def _cells(d):
        if "note" in d:
            return [d["note"], "-", "-", "-"]
        return [f"{d['peak_min']:g}", f"{d['measured_per_min']:.3g}", f"{d['predicted_per_min']:.3g}",
                "yes" if d["within_tolerance"] else "no"]

    for s in out["tdp_decline"]:
        tdp = _cells(s)
        if s["decisive"]:
            tdp[3] = "yes" if s.get("pass") else "NO"
        lines.append("| " + " | ".join([s["pot"], "decisive" if s["decisive"] else "reported", *tdp,
                                         *_cells(s["dpo"])]) + " |")
    lines += ["", f"## The engine's in-situ pH at {PH_REPORT_MIN:g} min (Zhai 2021, 100 C; the pH the split reads)", "",
              "| pot | in-situ pH |", "|---|---|"]
    lines += [f"| {k} | {v:.2f} |" for k, v in sorted(out["ph_in_situ_at_60_min"].items())]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: {passed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
