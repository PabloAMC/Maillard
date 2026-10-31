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
        log10(predicted / measured) within +-0.3.
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
TOL_TIME, TOL_PH, TOL_DECLINE = 0.3, 0.15, 0.3


def _slope(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.polyfit(x, y, 1)[0]) if len(x) >= 2 else float("nan")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", required=True, help="the joint fit's stem suffix (s042_r3en_fixed, ...)")
    args = parser.parse_args(argv)
    fit_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json"
    fit = json.loads(fit_path.read_text())
    theta, arm = fit["fit"]["theta"], fit.get("arm") or sc._arm_of(args.fit)
    sc._configure(arm, float(fit["strecker_share"]))
    if fit.get("ea_static") is not None:
        joint.EA_STATIC = dict(fit["ea_static"])
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    pots, rows = fit_rows.load(TABLE)
    wanted = set(TIME_SLOPE_POTS) | {p for ps in PH_LADDERS.values() for p in ps} | set(DOUBLE_RATIO) | set(DECLINE_POTS)
    plans = {j["pot"]: j for j in probe.plan(pots, rows) if j["pot"] in wanted}
    meas: Dict[Tuple[str, str, float], float] = {}
    pred: Dict[Tuple[str, str, float], float] = {}
    for pot, job in plans.items():
        j = dict(job, **trial, strecker_share=float(fit["strecker_share"]))
        res = joint.run_task(j)
        if "failed" in res:
            raise SystemExit(f"{pot}: {res['failed']}")
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
    for pot in DECLINE_POTS:
        tdp = sorted((t, v) for (p, s, t), v in meas.items() if p == pot and s == "TDP")
        peak = max(range(len(tdp)), key=lambda i: tdp[i][1])
        after = tdp[peak:]
        if len(after) < 2:
            out["tdp_decline"].append({"pot": pot, "note": "no measured decline after the peak"})
            continue
        km = -_slope([t for t, _ in after], [math.log(v) for _, v in after])
        kp = -_slope([t for t, _ in after], [math.log(max(pred[(pot, "TDP", t)], 1e-300)) for t, _ in after])
        ok = km > 0 and kp > 0 and abs(math.log10(kp / km)) <= TOL_DECLINE
        out["tdp_decline"].append({"pot": pot, "peak_min": after[0][0], "measured_per_min": km,
                                   "predicted_per_min": kp, "pass": bool(ok)})

    passed = {"time_slope": sum(s["pass"] for s in out["time_slope"]),
              "ph_segments": sum(s["pass"] for s in out["ph_segments"]),
              "tdp_decline": sum(s.get("pass", False) for s in out["tdp_decline"])}
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
             f"{len(out['tdp_decline'])}.**" + ("" if payload["fit_converged"] else
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
              f"## (vi) TDP's apparent first-order decline after its measured peak (Zhai 2023 JAFC), tolerance "
              f"{TOL_DECLINE} dex", "", "| pot | from (min) | measured /min | predicted /min | pass |", "|---|---|---|---|---|"]
    lines += [f"| {s['pot']} | {s.get('peak_min', '-')} | {s.get('measured_per_min', float('nan')):.3g} | "
              f"{s.get('predicted_per_min', float('nan')):.3g} | {'yes' if s.get('pass') else 'NO'} |"
              for s in out["tdp_decline"]]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}: {passed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
