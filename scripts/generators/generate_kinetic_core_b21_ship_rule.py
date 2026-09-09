#!/usr/bin/env python
"""
Wave B21 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-09).

`results/validation/kinetic_core_b21_prereg.md` sec. 4, computed from the frozen B21 report (which carries
the before / after diagnostics) and the live panel, written to `kinetic_core_b21_ship_rule.{json,md}`.

  T1  the decisive rows (glucosone formation at 110, 120, 140 C; glyoxal formation at 120 C) within 0.3 dex
  T2  nothing shipped breaks: (a) the B1 browning hold-out median fold below 2 and within-3x fraction 1.0;
      (b) no currently scored panel row moves by more than 0.3 dex (the changes listed)
  T3  Quan 2020's glyoxal level at 21 min within the printed range widened by 0.5 dex, at 100 and 130 C
  T4  Xia 2022: glyoxal above methylglyoxal at 130 C / 80 min (reported)
  T5  Leahy 1989's total pyrazine miss in decades against B18's 2.9 (reported)
  T6  identification: Laplace sigma below one decade on both coordinates, neither on its bound
Ship rule: SHIP if T1, T2, T3 and T6 hold; T4 and T5 are reported.
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
REPORT = V / "kinetic_core_b21_fit_report.json"
SCORECARD = V / "core_panel_scores.json"
OUT = V / "kinetic_core_b21_ship_rule.json"
DECISIVE = ("ham_k_ama_g_110C", "ham_k_ama_g_120C", "ham_k_ama_g_140C", "ham_k_g_go_120C")


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def t1(rep):
    rows = {r["id"]: r["residual_dex"] for r in rep["rows"]}
    dec = {k: rows[k] for k in DECISIVE}
    worst = max(dec.items(), key=lambda kv: abs(kv[1]))
    return {"decisive_rows_dex": dec, "reported_rows_dex": {k: v for k, v in rows.items() if k not in DECISIVE},
            "worst_row": worst[0], "worst_dex": worst[1], "pass": bool(all(abs(v) <= 0.3 for v in dec.values()))}


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


def t2(rep):
    b = rep["diagnostics"]["before_glass"]["b1_browning_holdout"]
    a = rep["diagnostics"]["after_candidate"]["b1_browning_holdout"]
    browning = {"before": b, "after": a, "pass": bool(a["median_fold_error"] < 2.0 and a["fraction_within_3x"] >= 1.0)}
    from src.kinetic_core import scoring

    try:
        tracked = json.loads(subprocess.check_output(["git", "show", "HEAD:" + data_paths.rel(SCORECARD)], cwd=ROOT, text=True))
    except Exception as exc:  # pragma: no cover
        return {"browning": browning, "panel": {"status": f"tracked scorecard unavailable: {exc}", "pass": False}, "pass": False}
    live = scoring.score_panel()
    old = {bb["benchmark_id"]: _predicted_leaves(bb.get("compounds")) for bb in tracked["benchmarks"]}
    new = {bb["benchmark_id"]: _predicted_leaves(bb.get("compounds")) for bb in live["benchmarks"]}
    deltas = []
    for bid, leaves in old.items():
        for path, v in leaves.items():
            w = new.get(bid, {}).get(path)
            if w is None or v <= 0 or w <= 0:
                continue
            deltas.append((math.log10(w / v), bid, path))
    moved = sorted([d for d in deltas if abs(d[0]) > 0.01], key=lambda d: -abs(d[0]))
    worst = max((abs(d[0]) for d in deltas), default=0.0)
    panel = {"n_compared": len(deltas), "worst_abs_dex": worst, "moved_over_0.01_dex": [{"dex": d[0], "benchmark": d[1], "path": d[2]} for d in moved[:20]],
             "pass": bool(worst <= 0.3)}
    return {"browning": browning, "panel": panel, "pass": bool(browning["pass"] and panel["pass"])}


def t3(rep):
    q = rep["diagnostics"]["after_candidate"]["quan2020_glyoxal"]
    return {"quan": q, "pass": bool(all(v["within_range_widened_0.5_dex"] for v in q.values()))}


def t4(rep):
    return {"xia": rep["diagnostics"]["after_candidate"]["xia2022_ordering_130C_80min"], "reported_only": True}


def t5(rep):
    return {"leahy": rep["diagnostics"]["after_candidate"]["leahy1989_total_pyrazine_95C_2h"], "reported_only": True}


def t6(rep):
    lap = rep["laplace"]
    return {"sigma_dex": lap["sigma"], "on_bound": lap["on_bound"],
            "pass": bool(all(s < 1.0 for s in lap["sigma"].values()) and not any(lap["on_bound"].values()))}


def render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4, T5, T6 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"], p["T6"]
    bb, ba = T2["browning"]["before"], T2["browning"]["after"]
    L = [f"# Wave B21 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         f"Cost {p['cost']:.2f} on {p['n_rows']} rows; fitted (log10 at 100 C): {json.dumps({k: round(v, 3) for k, v in p['frozen'].items()})}", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 decisive rows | worst {T1['worst_row']} {T1['worst_dex']:+.2f} dex; {json.dumps({k: round(v, 2) for k, v in T1['decisive_rows_dex'].items()})}; reported {json.dumps({k: round(v, 2) for k, v in T1['reported_rows_dex'].items()})} | {T1['pass']} |",
         f"| T2a browning hold-out | median fold {bb['median_fold_error']:.2f} -> {ba['median_fold_error']:.2f}; within 3x {bb['fraction_within_3x']:.2f} -> {ba['fraction_within_3x']:.2f} | {T2['browning']['pass']} |",
         f"| T2b panel | {T2['panel'].get('n_compared')} predicted numbers; worst change {T2['panel'].get('worst_abs_dex', float('nan')):.3f} dex; moved over 0.01 dex: {len(T2['panel'].get('moved_over_0.01_dex', []))} | {T2['panel']['pass']} |",
         f"| T3 Quan 2020 glyoxal | {json.dumps({k: {'model_mmol_l': round(v['model_go_mmol_l'], 4), 'printed': v['printed_range'], 'dex_from_range': round(v['dex_from_range'], 2)} for k, v in T3['quan'].items()})} | {T3['pass']} |",
         f"| T4 Xia 2022 ordering | glyoxal {T4['xia']['go_mmol_l']:.2f} vs methylglyoxal {T4['xia']['mgo_mmol_l']:.2f} mmol/L at 130 C / 80 min; glyoxal above: {T4['xia']['go_above_mgo']} | reported |",
         f"| T5 Leahy total pyrazine | model {T5['leahy']['model_ug_per_l']:.3g} ug/L vs 13100 ({T5['leahy']['dex']:+.2f} dex; B18 {T5['leahy']['b18_dex']:+.2f}) | reported |",
         f"| T6 identification | sigma {json.dumps({k: round(v, 3) for k, v in T6['sigma_dex'].items()})}; on bound {any(T6['on_bound'].values())} | {T6['pass']} |", ""]
    if T2["panel"].get("moved_over_0.01_dex"):
        L += ["## Panel rows that moved (log10, live minus tracked)", ""] + [f"- {m['benchmark']} {m['path']}: {m['dex']:+.3f}" for m in T2["panel"]["moved_over_0.01_dex"]] + [""]
    return "\n".join(L)


def main() -> int:
    rep = _read(REPORT)
    T1, T2, T3, T4, T5, T6 = t1(rep), t2(rep), t3(rep), t4(rep), t5(rep), t6(rep)
    ships = bool(T1["pass"] and T2["pass"] and T3["pass"] and T6["pass"])
    payload = {"artifact": "kinetic_core_b21_ship_rule",
               "provenance": provenance.provenance_block("kinetic_core_b21_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b21_ship_rule.py", inputs=[REPORT]),
               "prereg": data_paths.rel(V / "kinetic_core_b21_prereg.md"),
               "rule": "SHIP if T1 (decisive rows within 0.3 dex), T2 (browning hold-out median below 2 and within 3x, no panel row moves 0.3 dex), T3 (Quan's glyoxal within the widened range) and T6 (both coordinates identified) hold; T4 and T5 reported",
               "cost": rep["objective"]["final_cost"], "n_rows": rep["objective"]["n_rows"], "frozen": rep["frozen_parameters"]["aqueous_glyoxal"],
               "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5, "T6": T6, "verdict": "SHIP" if ships else "DO NOT SHIP"}
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], "| T2", T2["pass"], T2["panel"].get("worst_abs_dex"), "| T3", T3["pass"], "| T6", T6["pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
