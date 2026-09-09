#!/usr/bin/env python
"""
Wave B22 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-09).

`results/validation/kinetic_core_b22_prereg.md` sec. 4, computed from the frozen B22 report and the live
panel, written to `kinetic_core_b22_ship_rule.{json,md}`.

  T1  the six methional and methanethiol rows within 0.3 dex, the three disulfide rows within 0.5, the
      release barrier off its bounds
  T2  no currently scored panel row moves by more than 0.05 dex
  T3  Deng 2022's methional at 120 min within one decade and rising from 30 to 120 min (reported)
  T4  identification: Laplace sigma below one decade on every coordinate; the identity ratio inside its band
  T5  Chin & Lindsay 1994's methanethiol half-life at 30 C (reported)
Ship rule: SHIP if T1, T2 and T4 hold; T3 and T5 are reported.
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
REPORT = V / "kinetic_core_b22_fit_report.json"
SCORECARD = V / "core_panel_scores.json"
OUT = V / "kinetic_core_b22_ship_rule.json"


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def t1(rep):
    rows = {r["id"]: r for r in rep["rows"]}
    tight = {k: v["residual_dex"] for k, v in rows.items() if v["species"] in ("MTAL", "MSH")}
    loose = {k: v["residual_dex"] for k, v in rows.items() if v["species"] == "DMDS"}
    ea_on_bound = rep["laplace"]["on_bound"]["ea_mtal_msh_kj_mol"]
    worst = max(tight.items(), key=lambda kv: abs(kv[1]))
    return {"methional_and_methanethiol_dex": tight, "disulfide_dex": loose, "worst_row": worst[0], "worst_dex": worst[1],
            "release_barrier_on_bound": ea_on_bound,
            "pass": bool(all(abs(v) <= 0.3 for v in tight.values()) and all(abs(v) <= 0.5 for v in loose.values()) and not ea_on_bound)}


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
    deltas = [abs(math.log10(new[bid][path] / v)) for bid, leaves in old.items() for path, v in leaves.items()
              if new.get(bid, {}).get(path) and v > 0]
    worst = max(deltas, default=0.0)
    return {"n_compared": len(deltas), "worst_dex": worst, "pass": bool(worst < 0.05)}


def t3(rep):
    d = rep["diagnostics"]["deng2022_methional_120C"]
    return {"deng_120min_dex": d["120min"]["dex"], "rising_30_to_120": d["rising_30_to_120"], "all": {k: v for k, v in d.items() if isinstance(v, dict)},
            "pass": bool(abs(d["120min"]["dex"]) <= 1.0 and d["rising_30_to_120"]), "reported_only": True}


def t4(rep):
    lap = rep["laplace"]
    thresh = {"log10_identity_ratio_met_over_gly": 1.0, "log10_k_mtal_msh_100C": 1.0, "ea_mtal_msh_kj_mol": 60.0, "log10_k_msh_dmds_100C": 1.0}
    ratio = rep["frozen_parameters"]["methionine"]["log10_identity_ratio_met_over_gly"]
    inside = bool(-2.0 + 1e-3 < ratio < 2.0 - 1e-3)
    return {"sigma": lap["sigma"], "on_bound": lap["on_bound"], "ratio_inside_band": inside,
            "pass": bool(all(lap["sigma"][k] < thresh[k] for k in thresh) and inside and not any(lap["on_bound"].values()))}


def t5(rep):
    return {**rep["diagnostics"]["chin1994_methanethiol_half_life_30C"], "reported_only": True}


def render(p):
    T1, T2, T3, T4, T5 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"]
    L = [f"# Wave B22 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         f"Cost {p['cost']:.0f} on {p['n_rows']} rows; optimum {json.dumps({k: round(v, 3) for k, v in p['frozen'].items()})}", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 rows | methional and methanethiol {json.dumps({k: round(v, 2) for k, v in T1['methional_and_methanethiol_dex'].items()})}; disulfide {json.dumps({k: round(v, 2) for k, v in T1['disulfide_dex'].items()})}; release barrier on bound {T1['release_barrier_on_bound']} | {T1['pass']} |",
         f"| T2 panel | {T2.get('n_compared')} numbers, worst change {T2.get('worst_dex', float('nan')):.2e} dex | {T2['pass']} |",
         f"| T3 Deng 2022 | methional at 120 min {T3['deng_120min_dex']:+.2f} dex; rising 30 -> 120 min {T3['rising_30_to_120']}; all {json.dumps({k: round(v['dex'], 2) for k, v in T3['all'].items()})} | reported ({T3['pass']}) |",
         f"| T4 identification | sigma {json.dumps({k: round(v, 2) for k, v in T4['sigma'].items()})}; on bound {json.dumps(T4['on_bound'])}; ratio inside band {T4['ratio_inside_band']} | {T4['pass']} |",
         f"| T5 Chin & Lindsay 1994 | methanethiol half-life at 30 C: model {T5['model_min']:.3g} min vs {T5['chin_min_with_cu']} with copper | reported |", ""]
    return "\n".join(L)


def main() -> int:
    rep = _read(REPORT)
    T1, T2, T3, T4, T5 = t1(rep), t2(), t3(rep), t4(rep), t5(rep)
    ships = bool(T1["pass"] and T2["pass"] and T4["pass"])
    payload = {"artifact": "kinetic_core_b22_ship_rule",
               "provenance": provenance.provenance_block("kinetic_core_b22_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b22_ship_rule.py", inputs=[REPORT]),
               "prereg": data_paths.rel(V / "kinetic_core_b22_prereg.md"),
               "rule": "SHIP if T1 (methional and methanethiol rows within 0.3 dex, disulfide within 0.5, release barrier off its bounds), T2 (no panel row moves 0.05 dex) and T4 (every coordinate identified, the ratio inside its band) hold; T3 and T5 reported",
               "cost": rep["objective"]["final_cost"], "n_rows": rep["objective"]["n_rows"], "frozen": rep["frozen_parameters"]["methionine"],
               "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5, "verdict": "SHIP" if ships else "DO NOT SHIP"}
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], T1["worst_dex"], "| T2", T2["pass"], "| T4", T4["pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
