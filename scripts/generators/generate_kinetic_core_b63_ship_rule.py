#!/usr/bin/env python
"""
Wave B63 -- THE PRE-REGISTERED RULE, evaluated (2026-10-28).

`results/validation/kinetic_core_b63_prereg.md`. Under WAVES.md "Evidence decides, scores diagnose" the
verdict is the evidence audit; its one mechanical gate is V1 below. T1-T4 are recorded DIAGNOSIS. Three arms,
each in its own child process on what ships otherwise, ``parameters_sulfur.SINGLE_TEMPERATURE_BARRIERS``
switched per arm:

  baseline  the shipped fit, switch OFF
  (A)       the B63 refit, switch ON
  (B)       the shipped fit, switch ON (no refit)

  V1  every route-sharing single-temperature constant carries its class median in the engine, k(145 C) unchanged

  T1  Yiltirak's four FFT hold-out rows: median |log10 error| falls by >= 0.3 dex against the baseline
  T2  hold-out sulfur rows: within-3x not lower than the baseline's by more than 2; median |log10 error| not worse
      by more than max(0.1 dex, one bootstrap SE of the baseline median)
  T3  Hofmann ribose + xylose pH 5 rows: median |log10 error| not worse by more than the same tolerance
  T4  directional panel, independent agree: not fewer than the baseline's minus 3
Verdict: SHIP (A) if (A) passes T1-T4; else SHIP (B) if (B) does; else DO NOT SHIP.
Recorded, not gating: Kang's furfural 120/100 C ratio, the Zhou 2023c and Resconi furfural rows, Yiltirak's FFT
slope 130/100 C, each arm's fit cost.
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
B63_REPORT = V / "kinetic_core_b63_fit_report.json"
OUT = V / "kinetic_core_b63_ship_rule.json"
FFT = "2-Furfurylthiol (FFT)"
YILTIRAK = ("mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026", "mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026",
            "mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026", "mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026")
T3_BUNDLES = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5")
KANG_FUR = ("kang_100C_FUR", "kang_120C_FUR")


def _shipped_report() -> Path:
    from src.kinetic_core import engine

    return Path(engine._B2_FIT_REPORT)


def _child(mode: str, report: Path) -> Dict[str, Any]:
    from src.kinetic_core import engine, parameters_sulfur as ps

    ps.SINGLE_TEMPERATURE_BARRIERS = dataclasses.replace(ps.SINGLE_TEMPERATURE_BARRIERS, enabled=(mode == "on"))
    engine._B2_FIT_REPORT = Path(report)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    rows = []
    for b in score_panel()["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"), "within": bool(r.get("within_band")),
                         "predicted": r.get("predicted"), "measured": r.get("measured")})
    ind = score_directional()["summary"]["independent"]["total"]
    return {"mode": mode, "report": str(report), "rows": rows,
            "directional_independent": {"agree": ind["agree"], "evaluable": ind["evaluable"]}}


def _cost_child(mode: str, report: Path) -> Dict[str, Any]:
    """The 52-row objective at a report's vector with the switch set; also Kang's furfural predictions."""
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    from src.kinetic_core import parameters_sulfur as ps

    ps.SINGLE_TEMPERATURE_BARRIERS = dataclasses.replace(ps.SINGLE_TEMPERATURE_BARRIERS, enabled=(mode == "on"))
    import generate_kinetic_core_b63_fit as FIT  # the shipped objective; its own switch is overridden next line
    from generate_kinetic_core_b8_laplace import frozen_vector

    ps.SINGLE_TEMPERATURE_BARRIERS = dataclasses.replace(ps.SINGLE_TEMPERATURE_BARRIERS, enabled=(mode == "on"))
    r = FIT.residual_vector(frozen_vector(json.loads(Path(report).read_text(encoding="utf-8"))), True)
    ids = [row["id"] for row in FIT.B23.ACTIVE_FIT_ROWS]
    sig = {row["id"]: float(row["sigma_log"]) for row in FIT.B23.ACTIVE_FIT_ROWS}
    tgt = {row["id"]: float(row["target"]) for row in FIT.B23.ACTIVE_FIT_ROWS}
    kang = {k: tgt[k] * 10 ** (float(r[ids.index(k)]) * sig[k]) for k in KANG_FUR if k in ids}
    return {"mode": mode, "report": str(report), "cost": 0.5 * float(r @ r), "kang_fur_pred": kang,
            "kang_fur_target": {k: tgt[k] for k in KANG_FUR if k in tgt}}


def _v1_child() -> Dict[str, Any]:
    """V1: with the switch on, every route-sharing target carries its class median and k(145 C) is unchanged."""
    from src.kinetic_core import ea_class_prior, engine, parameters_sulfur as ps

    off = engine.core_parameters(engine.SULFUR)
    ps.SINGLE_TEMPERATURE_BARRIERS = dataclasses.replace(ps.SINGLE_TEMPERATURE_BARRIERS, enabled=True)
    on = engine.core_parameters(engine.SULFUR)
    bad = []
    for key, _klass, shipped in ea_class_prior.single_temperature_targets():
        if shipped is None:
            continue
        want = ps.ROUTE_CLASS_MEDIAN_EA_KJ_MOL[ps.FORMATION_ROUTE_OF[key]]
        if abs(on[key].ea_kj_mol - want) > 1e-9 or on[key].k_ref != off[key].k_ref:
            bad.append(key)
    moved = sorted(k for k in on if getattr(on[k], "ea_kj_mol", None) != getattr(off[k], "ea_kj_mol", None))
    return {"pass": not bad, "bad": bad, "moved": moved}


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


def _bootstrap_se_median(values: List[float], n: int = 2000, seed: int = 62) -> float:
    if len(values) < 2:
        return 0.0
    rng = random.Random(seed)
    return statistics.pstdev([statistics.median(rng.choices(values, k=len(values))) for _ in range(n)])


def _yil(s):
    return [r for r in s["rows"] if r["benchmark"] in YILTIRAK and r["compound"] == FFT]


def _hold(s):
    return [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]


def _hof(s):
    return [r for r in s["rows"] if r["benchmark"] in T3_BUNDLES]


def _tests(base: Dict[str, Any], arm: Dict[str, Any]) -> Dict[str, Any]:
    b1, a1 = _summary(_yil(base)), _summary(_yil(arm))
    t1 = {"baseline": b1, "arm": a1, "drop_dex": b1["median_abs_log10"] - a1["median_abs_log10"],
          "rows": {r["benchmark"]: r["fold"] for r in _yil(arm)}}
    t1["pass"] = bool(a1["rows"] == 4 and t1["drop_dex"] >= 0.3)
    tol2 = max(0.1, _bootstrap_se_median(_abs_logs(_hold(base))))
    tol3 = max(0.1, _bootstrap_se_median(_abs_logs(_hof(base))))
    t2 = {"baseline": _summary(_hold(base)), "arm": _summary(_hold(arm)), "tolerance_dex": tol2}
    t2["pass"] = bool(t2["arm"]["within_3x"] >= t2["baseline"]["within_3x"] - 2
                      and t2["arm"]["median_abs_log10"] <= t2["baseline"]["median_abs_log10"] + tol2)
    t3 = {"baseline": _summary(_hof(base)), "arm": _summary(_hof(arm)), "tolerance_dex": tol3}
    t3["pass"] = bool(t3["arm"]["median_abs_log10"] <= t3["baseline"]["median_abs_log10"] + tol3)
    t4 = {"baseline": base["directional_independent"], "arm": arm["directional_independent"]}
    t4["pass"] = bool(arm["directional_independent"]["agree"] >= base["directional_independent"]["agree"] - 3)
    return {"T1_yiltirak_fft": t1, "T2_holdout_sulfur": t2, "T3_hofmann_ph5": t3, "T4_directional": t4,
            "union_panel": {"baseline": _summary(base["rows"]), "arm": _summary(arm["rows"])},
            "pass": bool(t1["pass"] and t2["pass"] and t3["pass"] and t4["pass"])}


def _recorded(scored: Dict[str, Dict[str, Any]], costs: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"fit_cost": {a: c["cost"] for a, c in costs.items()}, "kang_fur_120_over_100": {},
                           "furfural_rows": {}, "yiltirak_fft_slope_130_over_100": {}}
    for arm, c in costs.items():
        p, t = c["kang_fur_pred"], c["kang_fur_target"]
        if len(p) == 2:
            out["kang_fur_120_over_100"][arm] = {"model": p["kang_120C_FUR"] / p["kang_100C_FUR"],
                                                 "measured": t["kang_120C_FUR"] / t["kang_100C_FUR"]}
    for arm, s in scored.items():
        out["furfural_rows"][arm] = {r["benchmark"]: r["fold"] for r in s["rows"]
                                     if str(r["compound"]).lower() == "furfural"}
        y = {r["benchmark"]: r for r in _yil(s)}
        lo, hi = y.get(YILTIRAK[0]), y.get(YILTIRAK[3])
        if lo and hi and lo["predicted"] and hi["predicted"]:
            out["yiltirak_fft_slope_130_over_100"][arm] = {"model": hi["predicted"] / lo["predicted"],
                                                           "measured": hi["measured"] / lo["measured"]}
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", nargs=2, metavar=("MODE", "REPORT"))
    parser.add_argument("--cost", nargs=2, metavar=("MODE", "REPORT"))
    parser.add_argument("--v1", action="store_true")
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child[0], Path(args.child[1]))))
        return 0
    if args.cost:
        print(json.dumps(_cost_child(args.cost[0], Path(args.cost[1]))))
        return 0
    if args.v1:
        print(json.dumps(_v1_child()))
        return 0
    shipped = _shipped_report()
    arms_def = {"baseline": ("off", shipped), "A_refit": ("on", B63_REPORT), "B_frozen_shipped": ("on", shipped)}
    scored = {arm: _run("--child", mode, str(report)) for arm, (mode, report) in arms_def.items()}
    arms = {arm: _tests(scored["baseline"], scored[arm]) for arm in ("A_refit", "B_frozen_shipped")}
    v1 = _run("--v1")
    verdict = ("ACCEPT (A) on the evidence audit: V1 passes; T1-T4 recorded as diagnosis" if v1["pass"]
               else "HOLD: V1 fails (the engine does not carry the declared barriers)")
    costs = {arm: _run("--cost", mode, str(report)) for arm, (mode, report) in arms_def.items()}
    payload = {
        "artifact": "kinetic_core_b63_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b63_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b63_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "sulfur.py", ROOT / "src" / "kinetic_core" / "parameters_sulfur.py",
                    shipped, B63_REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b63_prereg.md"), "baseline_report": data_paths.rel(shipped),
        "V1": v1, "arms": arms, "recorded": _recorded(scored, costs), "verdict": verdict,
        "rule": "WAVES.md 'Evidence decides, scores diagnose': accept on the audit; V1 is the mechanical gate; T1-T4 are diagnosis",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    names = ("T1_yiltirak_fft", "T2_holdout_sulfur", "T3_hofmann_ph5", "T4_directional")
    a, b = arms["A_refit"], arms["B_frozen_shipped"]
    print(f"{verdict} | A: " + " ".join(f"{n[:2]} {a[n]['pass']}" for n in names)
          + " | B: " + " ".join(f"{n[:2]} {b[n]['pass']}" for n in names))
    return 0


def _fmt(s: Dict[str, Any]) -> str:
    m = s["median_abs_log10"]
    return f"{s['within_3x']}/{s['rows']} ({m:.2f} dex)" if m is not None else f"{s['within_3x']}/{s['rows']}"


def render(payload: Dict[str, Any]) -> str:
    a, b = payload["arms"]["A_refit"], payload["arms"]["B_frozen_shipped"]
    rec = payload["recorded"]
    lines = ["# Wave B63 — single-temperature barriers at their measured class's median", "",
             f"V1: {'PASS' if payload['V1']['pass'] else 'FAIL ' + str(payload['V1']['bad'])}; barriers moved on "
             f"{len(payload['V1']['moved'])} constants.", "",
             f"*Pre-registration: `{payload['prereg']}`. Baseline: `{payload['baseline_report']}`. Rule: {payload['rule']}.*", "",
             f"## Verdict: **{payload['verdict']}**", "",
             "| diagnostic | baseline (switch off) | (A) B63 refit | would pass | (B) shipped fit, switch on | would pass |",
             "| --- | --- | --- | --- | --- | --- |"]
    ya, yb = a["T1_yiltirak_fft"], b["T1_yiltirak_fft"]
    lines.append(f"| T1 Yiltirak FFT, median dex (must fall >= 0.3) | {ya['baseline']['median_abs_log10']:.2f} | "
                 f"{ya['arm']['median_abs_log10']:.2f} (-{ya['drop_dex']:.2f}) | {ya['pass']} | "
                 f"{yb['arm']['median_abs_log10']:.2f} (-{yb['drop_dex']:.2f}) | {yb['pass']} |")
    for key, label in (("T2_holdout_sulfur", "T2 hold-out sulfur rows"), ("T3_hofmann_ph5", "T3 Hofmann ribose + xylose pH 5")):
        lines.append(f"| {label} (tolerance {a[key]['tolerance_dex']:.2f} dex) | {_fmt(a[key]['baseline'])} | "
                     f"{_fmt(a[key]['arm'])} | {a[key]['pass']} | {_fmt(b[key]['arm'])} | {b[key]['pass']} |")
    d0, da, db = a["T4_directional"]["baseline"], a["T4_directional"]["arm"], b["T4_directional"]["arm"]
    lines.append(f"| T4 directional, independent agree | {d0['agree']}/{d0['evaluable']} | {da['agree']}/{da['evaluable']} | "
                 f"{a['T4_directional']['pass']} | {db['agree']}/{db['evaluable']} | {b['T4_directional']['pass']} |")
    lines.append(f"| union panel within 3x (recorded) | {_fmt(a['union_panel']['baseline'])} | {_fmt(a['union_panel']['arm'])} | — | "
                 f"{_fmt(b['union_panel']['arm'])} | — |")
    lines += ["", "## Recorded, not gating", "", "| arm | fit cost | Kang furfural 120/100 (measured 1.7) | Yiltirak FFT 130/100 |",
              "| --- | --- | --- | --- |"]
    for arm, cost in rec["fit_cost"].items():
        k = rec["kang_fur_120_over_100"].get(arm, {})
        y = rec["yiltirak_fft_slope_130_over_100"].get(arm, {})
        lines.append(f"| {arm} | {cost:.3f} | {k.get('model', float('nan')):.2f} (measured {k.get('measured', float('nan')):.2f}) | "
                     f"{y.get('model', float('nan')):.3g} (measured {y.get('measured', float('nan')):.3g}) |")
    lines += ["", "Furfural rows, fold error:", ""]
    for arm, rows in rec["furfural_rows"].items():
        lines.append(f"- {arm}: " + ", ".join(f"{k} {v:.3g}" for k, v in rows.items()))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
