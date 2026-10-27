#!/usr/bin/env python
"""
Wave B48 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b48_prereg.md` sec. 4. The union panel and the directional panel are
scored twice, each in a child process: with the 2,4-decadienal breakage step OFF and ON
(``DECADIENAL_BREAKAGE.enabled``). Writes `results/validation/kinetic_core_b48_ship_rule.{json,md}`.

  T1  hold-out: pure 2,4-decadienal, 60 min at 180 C, share to hexanal within [6 %, 21 %]
      (Zamora 2007's 9-14 %, widened 1.5x each side for the medium)
  T2  scored lipid rows: within-3x not lower; median |log10 error| not worse by more than 0.05 dex
  T3  directional panel, independent claims agreeing: not fewer
  T4  recorded: every scored row that moves by more than 0.01 dex
Ship rule: SHIP if T1, T2 and T3 hold.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_b48_ship_rule.json"
HOLDOUT_BAND = (0.06, 0.21)


def _child(enabled: bool) -> Dict[str, Any]:
    from src.kinetic_core import lipid, parameters_lipid

    step = dataclasses.replace(parameters_lipid.DECADIENAL_BREAKAGE, enabled=enabled)
    parameters_lipid.DECADIENAL_BREAKAGE = step
    lipid.DECADIENAL_BREAKAGE = step
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    rows = []
    for b in score_panel()["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"),
                         "signed_log10": r.get("signed_log10_error", r.get("log10_error")),
                         "within": bool(r.get("within_band"))})
    ind = score_directional()["summary"]["independent"]["total"]
    return {"enabled": enabled, "rows": rows, "directional": {"agree": ind["agree"], "evaluable": ind["evaluable"]}}


def _score(enabled: bool) -> Dict[str, Any]:
    done = subprocess.run([sys.executable, __file__, "--child", "on" if enabled else "off"], cwd=ROOT,
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"child ({enabled}) failed:\n{done.stderr[-3000:]}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def _summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    folds = [r["fold"] for r in rows if r["fold"] is not None and math.isfinite(r["fold"]) and r["fold"] > 0]
    return {"rows": len(rows), "within_3x": sum(r["within"] for r in rows),
            "median_abs_log10": statistics.median(abs(math.log10(f)) for f in folds) if folds else None}


def _t1() -> Dict[str, Any]:
    from src.kinetic_core.parameters_lipid import DECADIENAL_BREAKAGE as D

    share = D.share_to_hexanal(180.0, 60.0)
    return {"share_to_hexanal_180C_60min": share, "zamora2007_printed": [0.09, 0.14], "band": list(HOLDOUT_BAND),
            "pass": HOLDOUT_BAND[0] <= share <= HOLDOUT_BAND[1]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", choices=("on", "off"), default=None)
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child == "on")))
        return 0
    off, on = _score(False), _score(True)
    lipid_off = [r for r in off["rows"] if r["lane"] == "lipid"]
    lipid_on = [r for r in on["rows"] if r["lane"] == "lipid"]
    t1 = _t1()
    t2 = {"off": _summary(lipid_off), "on": _summary(lipid_on)}
    t2["pass"] = (t2["on"]["within_3x"] >= t2["off"]["within_3x"]
                  and t2["on"]["median_abs_log10"] <= t2["off"]["median_abs_log10"] + 0.05)
    t3 = {"off": off["directional"], "on": on["directional"], "pass": on["directional"]["agree"] >= off["directional"]["agree"]}
    by_key = {(r["benchmark"], r["compound"]): r for r in on["rows"]}
    moved = []
    for r in off["rows"]:
        o = by_key.get((r["benchmark"], r["compound"]))
        if o and r["fold"] and o["fold"] and math.isfinite(r["fold"]) and math.isfinite(o["fold"]):
            delta = abs(math.log10(o["fold"])) - abs(math.log10(r["fold"]))
            if abs(delta) > 0.01:
                moved.append({"benchmark": r["benchmark"], "compound": r["compound"], "fold_off": r["fold"],
                              "fold_on": o["fold"], "delta_abs_log10": delta})
    verdict = "SHIP" if (t1["pass"] and t2["pass"] and t3["pass"]) else "DO NOT SHIP"
    payload = {
        "artifact": "kinetic_core_b48_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b48_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b48_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "parameters_lipid.py", ROOT / "src" / "kinetic_core" / "lipid.py"]),
        "prereg": data_paths.rel(V / "kinetic_core_b48_prereg.md"),
        "T1_holdout_zamora2007": t1, "T2_lipid_rows": t2, "T3_directional": t3, "T4_rows_moved": moved,
        "union_panel": {"off": _summary(off["rows"]), "on": _summary(on["rows"])},
        "verdict": verdict, "rule": "SHIP if T1 (Zamora 2007 hold-out), T2 (lipid rows) and T3 (directional) hold",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: T1 {100 * t1['share_to_hexanal_180C_60min']:.1f}% in [6, 21]; T2 lipid {t2['off']['within_3x']}/{t2['off']['rows']} "
          f"-> {t2['on']['within_3x']}/{t2['on']['rows']} (median {t2['off']['median_abs_log10']:.3f} -> {t2['on']['median_abs_log10']:.3f}); "
          f"T3 {t3['off']['agree']} -> {t3['on']['agree']}; rows moved > 0.01 dex: {len(moved)}")
    return 0


def render(payload: Dict[str, Any]) -> str:
    t1, t2, t3, moved = (payload[k] for k in ("T1_holdout_zamora2007", "T2_lipid_rows", "T3_directional", "T4_rows_moved"))
    u = payload["union_panel"]
    lines = [
        "# Wave B48 ship rule — the 2,4-decadienal breakage step", "",
        f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
        f"## Verdict: **{payload['verdict']}**", "",
        "| test | step off | step on | pass |", "| --- | --- | --- | --- |",
        f"| T1 Zamora 2007 hold-out (180 C, 60 min, share to hexanal; printed 9-14 %, band 6-21 %) | — | "
        f"{100 * t1['share_to_hexanal_180C_60min']:.1f} % | {t1['pass']} |",
        f"| T2 lipid rows within 3x | {t2['off']['within_3x']}/{t2['off']['rows']} (median {t2['off']['median_abs_log10']:.3f} dex) "
        f"| {t2['on']['within_3x']}/{t2['on']['rows']} (median {t2['on']['median_abs_log10']:.3f} dex) | {t2['pass']} |",
        f"| T3 directional, independent agree | {t3['off']['agree']}/{t3['off']['evaluable']} | {t3['on']['agree']}/{t3['on']['evaluable']} | {t3['pass']} |",
        f"| union panel within 3x (recorded) | {u['off']['within_3x']}/{u['off']['rows']} | {u['on']['within_3x']}/{u['on']['rows']} | — |",
        "", "### T4 rows that moved by more than 0.01 dex", "",
    ]
    if moved:
        lines += ["| bundle | compound | fold off | fold on | change in abs log10 |", "| --- | --- | ---: | ---: | ---: |"]
        lines += [f"| {m['benchmark']} | {m['compound']} | {m['fold_off']:.3g} | {m['fold_on']:.3g} | {m['delta_abs_log10']:+.3f} |" for m in moved]
    else:
        lines.append("None.")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
