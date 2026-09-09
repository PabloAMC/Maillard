#!/usr/bin/env python
"""
Wave B24 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-09).

`results/validation/kinetic_core_b24_prereg.md` sec. 4, from the frozen B24 report and the live panel.

  T1  the two fed-pyrroline rows within 0.3 dex; the three proline rows within 0.5 dex and in the printed order
  T2  no currently scored panel row moves by more than 0.05 dex
  T3  the apparent barrier of the whole cascade against Chan & Reineccius 1994's 60.2 kJ/mol (reported)
  T4  Hofmann's excess-pyrroline experiment (reported, expected to fail)
  T5  identification: Laplace sigma below one decade on both coordinates, neither on its bound
Ship rule: SHIP if T1, T2 and T5 hold; T3 and T4 are reported.
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
REPORT = V / "kinetic_core_b24_fit_report.json"
SCORECARD = V / "core_panel_scores.json"
OUT = V / "kinetic_core_b24_ship_rule.json"


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def t1(rep):
    rows = {r["id"]: r for r in rep["rows"]}
    fed = {k: v["residual_dex"] for k, v in rows.items() if v["basis"] == "PYRL"}
    pro = {k: v for k, v in rows.items() if v["basis"] == "PRO"}
    order_ok = pro["hof_t9_pro400_mgo4"]["predicted_mol_pct"] < pro["hof_t9_pro400_mgo40"]["predicted_mol_pct"] < pro["hof_t9_pro400_mgo400"]["predicted_mol_pct"]
    pro_dex = {k: v["residual_dex"] for k, v in pro.items()}
    return {"fed_pyrroline_dex": fed, "proline_dex": pro_dex, "proline_order_ok": bool(order_ok),
            "pass": bool(all(abs(v) <= 0.3 for v in fed.values()) and all(abs(v) <= 0.5 for v in pro_dex.values()) and order_ok)}


def _predicted_leaves(o, path=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(_predicted_leaves(v, f"{path}/{k}"))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(_predicted_leaves(v, f"{path}[{i}]"))
    elif isinstance(o, (int, float)) and not isinstance(o, bool) and "pred" in path.lower():
        out[path] = float(o)
    return out


def t2():
    from src.kinetic_core import scoring

    try:
        tracked = json.loads(subprocess.check_output(["git", "show", "HEAD:" + data_paths.rel(SCORECARD)], cwd=ROOT, text=True))
    except Exception as exc:  # pragma: no cover
        return {"status": f"tracked scorecard unavailable: {exc}", "pass": False}
    live = scoring.score_panel()
    old = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in tracked["benchmarks"]}
    new = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in live["benchmarks"]}
    deltas = [abs(math.log10(new[bid][path] / v)) for bid, leaves in old.items() for path, v in leaves.items() if new.get(bid, {}).get(path) and v > 0]
    worst = max(deltas, default=0.0)
    return {"n_compared": len(deltas), "worst_dex": worst, "pass": bool(worst < 0.05)}


def render(p):
    T1, T2, T3, T4, T5 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"]
    return "\n".join([f"# Wave B24 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
                      f"Cost {p['cost']:.2f} on {p['n_rows']} rows; optimum {json.dumps({k: round(v, 3) for k, v in p['frozen'].items()})}", "",
                      "| test | result | pass |", "|---|---|---|",
                      f"| T1 rows | fed 1-pyrroline {json.dumps({k: round(v, 2) for k, v in T1['fed_pyrroline_dex'].items()})}; proline {json.dumps({k: round(v, 2) for k, v in T1['proline_dex'].items()})}; order {T1['proline_order_ok']} | {T1['pass']} |",
                      f"| T2 panel | {T2.get('n_compared')} numbers, worst change {T2.get('worst_dex', float('nan')):.2e} dex | {T2['pass']} |",
                      f"| T3 apparent barrier | model {T3['model_kj_mol']:.0f} kJ/mol vs Chan & Reineccius 1994's {T3['chan1994b_kj_mol']} | reported |",
                      f"| T4 excess pyrroline | model {T4['model_mol_pct_of_mgo']:.3g} mol % of methylglyoxal vs printed {T4['printed_mol_pct_of_mgo']} ({T4['dex']:+.2f} dex) | reported |",
                      f"| T5 identification | sigma {json.dumps({k: round(v, 2) for k, v in T5['sigma'].items()})}; on bound {any(T5['on_bound'].values())} | {T5['pass']} |", ""])


def main() -> int:
    rep = _read(REPORT)
    T1, T2 = t1(rep), t2()
    T3 = rep["diagnostics"]["apparent_ea_glucose_proline_pot"]
    T4 = rep["diagnostics"]["hofmann_expt3_pyrl10_mgo2"]
    lap = rep["laplace"]
    T5 = {"sigma": lap["sigma"], "on_bound": lap["on_bound"], "pass": bool(all(s < 1.0 for s in lap["sigma"].values()) and not any(lap["on_bound"].values()))}
    ships = bool(T1["pass"] and T2["pass"] and T5["pass"])
    payload = {"artifact": "kinetic_core_b24_ship_rule",
               "provenance": provenance.provenance_block("kinetic_core_b24_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b24_ship_rule.py", inputs=[REPORT]),
               "prereg": data_paths.rel(V / "kinetic_core_b24_prereg.md"),
               "rule": "SHIP if T1 (fed rows within 0.3 dex, proline rows within 0.5 and in order), T2 (no panel row moves 0.05 dex) and T5 (both identified, off their bounds) hold; T3 and T4 reported",
               "cost": rep["objective"]["final_cost"], "n_rows": rep["objective"]["n_rows"], "frozen": rep["frozen_parameters"]["proline"],
               "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5, "verdict": "SHIP" if ships else "DO NOT SHIP"}
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], "| T2", T2["pass"], "| T5", T5["pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
