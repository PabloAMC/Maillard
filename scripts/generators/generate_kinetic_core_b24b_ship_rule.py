#!/usr/bin/env python
"""
Wave B24b -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-10).

`results/validation/kinetic_core_b24b_prereg.md` section 4, computed from the frozen B24b report.

  T1  the switch, DECISIVE: all three AP : ATHP ratios within 0.3 dex AND strictly increasing with
      the methylglyoxal charge. The ordering is the finding; getting it wrong fails outright
  T2  the pyrroline excess within 0.5 dex
  T3  B24's two fed rows within 0.3 dex of where B24 left them
  T4  identification, reported
  T4b the pH transfer, reported: B18's slopes on a ladder they were not fitted on
  T5  nothing else moves
Ship rule: SHIP if T1, T2, T3 and T5 hold; T4 and T4b reported.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
REPORT = V / "kinetic_core_b24b_fit_report.json"
OUT = V / "kinetic_core_b24b_ship_rule.json"


def main() -> int:
    rep = json.loads(REPORT.read_text(encoding="utf-8"))
    resid = rep["residual_by_row_dex"]
    pred = rep["predicted_by_row"]
    switch_ids = sorted((k for k in pred if k.startswith("switch_mgo")),
                        key=lambda k: float(k.replace("switch_mgo", "")))
    ordered = all(pred[a] < pred[b] for a, b in zip(switch_ids, switch_ids[1:]))
    worst_switch = max((abs(resid[k]), k) for k in switch_ids)
    model_span = pred[switch_ids[-1]] / pred[switch_ids[0]]
    printed_span = 12.8 / 0.16
    t1 = {"ordering_increasing": bool(ordered), "residuals_dex": {k: resid[k] for k in switch_ids},
          "worst": worst_switch[1], "worst_dex": worst_switch[0],
          "model_span_x": model_span, "printed_span_x": printed_span,
          "pass": bool(ordered and worst_switch[0] <= 0.3)}
    t2 = {"dex": resid["pyrl_excess_ap_molpct"],
          "predicted_mol_pct": pred["pyrl_excess_ap_molpct"],
          "pass": bool(abs(resid["pyrl_excess_ap_molpct"]) <= 0.5)}
    fed = rep["checks"]["b24_fed_rows"]
    t3 = {"rows": fed, "worst_dex": max(abs(r["dex"]) for r in fed.values()),
          "pass": bool(all(abs(r["dex"]) <= 0.3 for r in fed.values()))}
    lap = rep["laplace"]
    t4 = {"sigma": lap["sigma"], "identified": lap["identified"], "on_bound": lap["on_bound"],
          "reported_only": True}
    t4b = {"ladder": rep["checks"]["ph_ladder_vs_b18_transfer"],
           "worst_dex": max(abs(r["dex"]) for r in rep["checks"]["ph_ladder_vs_b18_transfer"].values()),
           "reported_only": True}
    try:
        tracked = json.loads(subprocess.check_output(
            ["git", "show", "HEAD:" + data_paths.rel(V / "core_panel_scores.json")], cwd=ROOT, text=True))
        from src.kinetic_core import scoring
        live = scoring.score_panel()
        same = (json.dumps([b.get("compounds") for b in tracked["benchmarks"]], sort_keys=True, default=str)
                == json.dumps([b.get("compounds") for b in live["benchmarks"]], sort_keys=True, default=str))
    except Exception as exc:  # pragma: no cover
        same, exc_msg = False, str(exc)
    t5 = {"panel_identical": bool(same), "pass": bool(same)}

    ships = bool(t1["pass"] and t2["pass"] and t3["pass"] and t5["pass"])
    payload = {
        "artifact": "kinetic_core_b24b_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b24b_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b24b_ship_rule.py", inputs=[REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b24b_prereg.md"),
        "rule": "SHIP if T1 (the switch), T2 (the pyrroline excess), T3 (B24's fed rows) and T5 hold",
        "T1": t1, "T2": t2, "T3": t3, "T4": t4, "T4b": t4b, "T5": t5,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=_render)
    print(payload["verdict"], "| T1", t1["pass"], "ordered", t1["ordering_increasing"],
          "| T2", t2["pass"], round(t2["dex"], 2), "| T3", t3["pass"], "| T5", t5["pass"])
    return 0


def _render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4, T4b, T5 = p["T1"], p["T2"], p["T3"], p["T4"], p["T4b"], p["T5"]
    L = [f"# Wave B24b ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration "
         f"`{p['prereg']}`.*", "", "| test | result | pass |", "|---|---|---|",
         f"| T1 the switch | ordering increasing {T1['ordering_increasing']}; worst {T1['worst']} "
         f"{T1['worst_dex']:+.2f} dex; model spans {T1['model_span_x']:.1f}x against a printed "
         f"{T1['printed_span_x']:.0f}x | {T1['pass']} |",
         f"| T2 pyrroline excess | {T2['predicted_mol_pct']:.3g} mol % against 0.33, "
         f"{T2['dex']:+.2f} dex | {T2['pass']} |",
         f"| T3 B24's fed rows | worst {T3['worst_dex']:+.3f} dex | {T3['pass']} |",
         f"| T4 identification | sigma {json.dumps({k: round(v, 2) for k, v in T4['sigma'].items()})}; "
         f"on bound {T4['on_bound']} | reported |",
         f"| T4b the pH transfer | worst {T4b['worst_dex']:+.2f} dex | reported |",
         f"| T5 nothing else moves | panel identical {T5['panel_identical']} | {T5['pass']} |", ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
