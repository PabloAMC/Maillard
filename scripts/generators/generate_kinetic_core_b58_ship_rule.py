#!/usr/bin/env python
"""
Wave B58 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b58_prereg.md` sec. 3. Three arms, each in its own child process, with the
B51 corrections on (what ships) and the B58 TTCA formation step switched per arm:

  baseline  B51's fit, step OFF  (what ships)
  (A)       the B58 refit, step ON
  (B)       B51's fit, step ON   (no refit)

  T1  mechanism: (a) Zhang 2025b's pot (Cys + Xyl 40 + 40 mM, 0.2 M phosphate pH 6.5, 120 C, 50 min): TTCA
      within a factor 2 of 45 % of the charged cysteine; (b) Zhai 2019's pot (Xyl 165 + Cys 82.6 mM, water,
      pH 7.4, 90 C, 160 min): TTCA within a factor 2 of 58 mM
  T2  hold-out sulfur rows: within-3x not lower than the baseline's by more than 2; median |log10 error| not
      worse by more than max(0.1 dex, one bootstrap SE of the baseline median)
  T3  Hofmann ribose + xylose pH 5 rows: median |log10 error| not worse by more than the same tolerance
  T4  directional panel, independent agree: not fewer than the baseline's minus 3
Verdict: SHIP (A) if (A) passes T1-T4; else SHIP (B) if (B) does; else DO NOT SHIP.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import random
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
B51_REPORT = V / "kinetic_core_b51_fit_report.json"
B58_REPORT = V / "kinetic_core_b58_fit_report.json"
OUT = V / "kinetic_core_b58_ship_rule.json"
T3_BUNDLES = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5")
ARMS = {"baseline": ("off", B51_REPORT), "A_refit": ("on", B58_REPORT), "B_frozen_b51": ("on", B51_REPORT)}
PROBES = {
    "zhang2025b": dict(precursors={"L-Cysteine": 40.0, "D-Xylose": 40.0}, temp_c=120.0, minutes=50.0, ph=6.5,
                       buffer=("phosphate", 0.2), cys0=40.0, measured_ttca=18.0, measured_label="45 % of Cys"),
    "zhai2019": dict(precursors={"L-Cysteine": 82.6, "D-Xylose": 165.2}, temp_c=90.0, minutes=160.0, ph=7.4,
                     buffer=("none", 0.0), cys0=82.6, measured_ttca=58.0, measured_label="58 mM plateau"),
}


def _child(mode: str, report: Path) -> Dict[str, Any]:
    from src.kinetic_core import engine, parameters_sulfur as ps
    from src.kinetic_core.ph_state import BufferSpec

    ps.TTCA_FORMATION = dataclasses.replace(ps.TTCA_FORMATION, enabled=(mode == "on"))
    engine._B2_FIT_REPORT = Path(report)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    rows = []
    for b in score_panel()["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"), "within": bool(r.get("within_band"))})
    ind = score_directional()["summary"]["independent"]["total"]
    probes = {}
    for name, p in PROBES.items():
        kind, molar = p["buffer"]
        buf = (BufferSpec(kind="phosphate", phosphate_mol_l=molar, declared=True, source=name) if kind == "phosphate"
               else BufferSpec(kind="none", declared=True, source=name))
        spec = engine.FormulationSpec(
            name=name, precursors=dict(p["precursors"]),
            process=engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(p["temp_c"], p["minutes"]),
                                       ph=p["ph"], buffer=buf))
        decl = engine.declare_envelope(spec, ["2-methyl-3-furanthiol"])
        state, _ = engine._integrate_program(engine.SULFUR, engine.core_parameters(engine.SULFUR),
                                             dict(decl.mapped_precursors), spec.process)
        probes[name] = {"ttca_mM": float(state.get("TTCA", 0.0)), "cys_free_mM": float(state.get("Cys", 0.0)),
                        "ttca_share_of_cys": float(state.get("TTCA", 0.0)) / p["cys0"]}
    return {"mode": mode, "report": str(report), "rows": rows, "probes": probes,
            "directional_independent": {"agree": ind["agree"], "evaluable": ind["evaluable"]}}


def _run(*args: str) -> Dict[str, Any]:
    done = subprocess.run([sys.executable, __file__, *args], cwd=ROOT, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"child {args} failed:\n{done.stderr[-3000:]}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def _abs_logs(rows: List[Dict[str, Any]]) -> List[float]:
    return [abs(math.log10(r["fold"])) for r in rows if r["fold"] and math.isfinite(r["fold"]) and r["fold"] > 0]


def _summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    logs = _abs_logs(rows)
    return {"rows": len(rows), "within_3x": sum(r["within"] for r in rows),
            "median_abs_log10": statistics.median(logs) if logs else None}


def _bootstrap_se_median(values: List[float], n: int = 2000, seed: int = 58) -> float:
    if len(values) < 2:
        return 0.0
    rng = random.Random(seed)
    meds = [statistics.median(rng.choices(values, k=len(values))) for _ in range(n)]
    return statistics.pstdev(meds)


def _tests(base: Dict[str, Any], arm: Dict[str, Any]) -> Dict[str, Any]:
    t1 = {}
    for name, p in PROBES.items():
        pred = arm["probes"][name]["ttca_mM"]
        lo, hi = p["measured_ttca"] / 2.0, min(p["measured_ttca"] * 2.0, p["cys0"])
        t1[name] = {"predicted_ttca_mM": pred, "measured_ttca_mM": p["measured_ttca"], "band_mM": [lo, hi],
                    "baseline_predicted_ttca_mM": base["probes"][name]["ttca_mM"], "pass": bool(lo <= pred <= hi)}
    t1["pass"] = all(v["pass"] for k, v in t1.items() if k != "pass")
    hold = lambda s: [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]
    hof = lambda s: [r for r in s["rows"] if r["benchmark"] in T3_BUNDLES]
    tol2 = max(0.1, _bootstrap_se_median(_abs_logs(hold(base))))
    tol3 = max(0.1, _bootstrap_se_median(_abs_logs(hof(base))))
    t2 = {"baseline": _summary(hold(base)), "arm": _summary(hold(arm)), "tolerance_dex": tol2}
    t2["pass"] = bool(t2["arm"]["within_3x"] >= t2["baseline"]["within_3x"] - 2
                      and t2["arm"]["median_abs_log10"] <= t2["baseline"]["median_abs_log10"] + tol2)
    t3 = {"baseline": _summary(hof(base)), "arm": _summary(hof(arm)), "tolerance_dex": tol3}
    t3["pass"] = bool(t3["arm"]["median_abs_log10"] <= t3["baseline"]["median_abs_log10"] + tol3)
    t4 = {"baseline": base["directional_independent"], "arm": arm["directional_independent"]}
    t4["pass"] = bool(arm["directional_independent"]["agree"] >= base["directional_independent"]["agree"] - 3)
    return {"T1_mechanism": t1, "T2_holdout_sulfur": t2, "T3_hofmann_ph5": t3, "T4_directional": t4,
            "union_panel": {"baseline": _summary(base["rows"]), "arm": _summary(arm["rows"])},
            "pass": bool(t1["pass"] and t2["pass"] and t3["pass"] and t4["pass"])}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", nargs=2, metavar=("MODE", "REPORT"))
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child[0], Path(args.child[1]))))
        return 0
    scored = {arm: _run("--child", mode, str(report)) for arm, (mode, report) in ARMS.items()}
    arms = {arm: _tests(scored["baseline"], scored[arm]) for arm in ("A_refit", "B_frozen_b51")}
    if arms["A_refit"]["pass"]:
        verdict = "SHIP (A): the B58 refit"
    elif arms["B_frozen_b51"]["pass"]:
        verdict = "SHIP (B): B51's fit with the TTCA formation step"
    else:
        verdict = "DO NOT SHIP"
    payload = {
        "artifact": "kinetic_core_b58_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b58_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b58_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "sulfur.py", ROOT / "src" / "kinetic_core" / "parameters_sulfur.py",
                    B51_REPORT, B58_REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b58_prereg.md"),
        "arms": arms, "probes": {a: scored[a]["probes"] for a in scored}, "verdict": verdict,
        "rule": "SHIP (A) if (A) passes T1-T4; else SHIP (B) if (B) does; else DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    a, b = arms["A_refit"], arms["B_frozen_b51"]
    print(f"{verdict} | A: T1 {a['T1_mechanism']['pass']} T2 {a['T2_holdout_sulfur']['pass']} T3 {a['T3_hofmann_ph5']['pass']} "
          f"T4 {a['T4_directional']['pass']} | B: T1 {b['T1_mechanism']['pass']} T2 {b['T2_holdout_sulfur']['pass']} "
          f"T3 {b['T3_hofmann_ph5']['pass']} T4 {b['T4_directional']['pass']}")
    return 0


def _fmt(s: Dict[str, Any]) -> str:
    m = s["median_abs_log10"]
    return f"{s['within_3x']}/{s['rows']} ({m:.2f} dex)" if m is not None else f"{s['within_3x']}/{s['rows']}"


def render(payload: Dict[str, Any]) -> str:
    a, b = payload["arms"]["A_refit"], payload["arms"]["B_frozen_b51"]
    lines = ["# Wave B58 ship rule — TTCA formation as a measured equilibrium", "",
             f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
             f"## Verdict: **{payload['verdict']}**", "",
             "| test | baseline (B51, step off) | (A) B58 refit | pass | (B) B51, step on | pass |",
             "| --- | --- | --- | --- | --- | --- |"]
    for name in PROBES:
        ta, tb = a["T1_mechanism"][name], b["T1_mechanism"][name]
        lines.append(f"| T1 {name}: TTCA mM (measured {ta['measured_ttca_mM']:g}, band {ta['band_mM'][0]:g}-{ta['band_mM'][1]:g}) | "
                     f"{ta['baseline_predicted_ttca_mM']:.3g} | {ta['predicted_ttca_mM']:.3g} | {ta['pass']} | {tb['predicted_ttca_mM']:.3g} | {tb['pass']} |")
    lines.append(f"| T2 hold-out sulfur rows (tolerance {a['T2_holdout_sulfur']['tolerance_dex']:.2f} dex) | {_fmt(a['T2_holdout_sulfur']['baseline'])} | "
                 f"{_fmt(a['T2_holdout_sulfur']['arm'])} | {a['T2_holdout_sulfur']['pass']} | {_fmt(b['T2_holdout_sulfur']['arm'])} | {b['T2_holdout_sulfur']['pass']} |")
    lines.append(f"| T3 Hofmann ribose + xylose pH 5 (tolerance {a['T3_hofmann_ph5']['tolerance_dex']:.2f} dex) | {_fmt(a['T3_hofmann_ph5']['baseline'])} | "
                 f"{_fmt(a['T3_hofmann_ph5']['arm'])} | {a['T3_hofmann_ph5']['pass']} | {_fmt(b['T3_hofmann_ph5']['arm'])} | {b['T3_hofmann_ph5']['pass']} |")
    d0, da, db = a["T4_directional"]["baseline"], a["T4_directional"]["arm"], b["T4_directional"]["arm"]
    lines.append(f"| T4 directional, independent agree | {d0['agree']}/{d0['evaluable']} | {da['agree']}/{da['evaluable']} | "
                 f"{a['T4_directional']['pass']} | {db['agree']}/{db['evaluable']} | {b['T4_directional']['pass']} |")
    lines.append(f"| union panel within 3x (recorded) | {_fmt(a['union_panel']['baseline'])} | {_fmt(a['union_panel']['arm'])} | — | "
                 f"{_fmt(b['union_panel']['arm'])} | — |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
