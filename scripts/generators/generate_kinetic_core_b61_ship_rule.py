#!/usr/bin/env python
"""
Wave B61 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-28).

`results/validation/kinetic_core_b61_prereg.md` sec. 3 (and Amendment 1). Three arms, each in its own child
process on what ships otherwise (B51's corrections, B60's numerics), with the B61 hydroxide-catalysed norfuraneol
decay switched per arm:

  baseline  B60's fit, step OFF  (what ships)
  (A)       the B61 refit, step ON
  (B)       B60's fit, step ON   (no refit)

  T1  the measurement is what the engine does: NF 1.75 mM alone, 0.2 M phosphate, 100 C, 60 min -- label pH 8.0:
      effective first-order loss within 0.1 dex of 0.032 1/min; label pH 7.0: >= 1.8e-3 1/min
  T2  the four Zhou 2023c MFT rows: median |log10 error| falls by >= 0.3 dex against the baseline
  T3  hold-out sulfur rows: within-3x not lower than the baseline's by more than 2; median |log10 error| not worse
      by more than max(0.1 dex, one bootstrap SE of the baseline median)
  T4  Hofmann ribose + xylose pH 5 rows: median |log10 error| not worse by more than the same tolerance
  T5  directional panel, independent agree: not fewer than the baseline's minus 3
Verdict: SHIP (A) if (A) passes T1-T5; else SHIP (B) if (B) does; else DO NOT SHIP.
Recorded, not gating: the Hofmann pH 3 / 5 / 7 sweep, Zhang 2025b and Yiltirak MFT, each arm's fit cost against
B60's, and T2 at the carbohydrate class's 10th / 90th percentile barrier.
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
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
B60_REPORT = V / "kinetic_core_b60_fit_report.json"
B61_REPORT = V / "kinetic_core_b61_fit_report.json"
OUT = V / "kinetic_core_b61_ship_rule.json"
ZHOU_PREFIX = "mp_holdout_arp_glutathione_"
MFT = "2-Methyl-3-furanthiol (MFT)"
T4_BUNDLES = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5")
SWEEP = {3.0: "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3", 5.0: "hofmann1998_ribose_cysteine_145C_20min_pH5",
         7.0: "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7"}
ARMS = {"baseline": ("off", B60_REPORT), "A_refit": ("on", B61_REPORT), "B_frozen_b60": ("on", B60_REPORT)}
MIKAMI_K = 0.032
MIKAMI_PH7_FLOOR = 1.8e-3
COST_B60 = 27.13196188533709
ADMISSIBLE_RISE = 3.59


def _mikami_k(label_ph: float) -> float:
    """The engine's effective first-order NF loss over 0-60 min in Mikami's NF-alone pot."""
    from src.kinetic_core import engine
    from src.kinetic_core.ph_state import BufferSpec

    captured = {}
    original = engine.integrate_sulfur

    def _capture(parameters, temperature_k, state, times, **kw):
        captured["run"] = original(parameters, temperature_k, state, times, **kw)
        return captured["run"]

    engine.integrate_sulfur = _capture
    try:
        spec = engine.FormulationSpec(name="mikami2017", precursors={"Norfuraneol": 1.75}, process=engine.ProcessSpec(
            thermal=engine.ThermalProgram.isothermal(100.0, 60.0), ph=label_ph,
            buffer=BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True, source="Mikami 2017")))
        decl = engine.declare_envelope(spec, ["2-methyl-3-furanthiol"])
        engine._integrate_program(engine.SULFUR, engine.core_parameters(engine.SULFUR),
                                  dict(decl.mapped_precursors), spec.process)
    finally:
        engine.integrate_sulfur = original
    run = captured["run"]
    return math.log(run.series("NF")[0] / run.final("NF")) / 60.0


def _child(mode: str, report: Path, ea: Optional[float]) -> Dict[str, Any]:
    from src.kinetic_core import engine, parameters_sulfur as ps

    step = dataclasses.replace(ps.NF_HYDROXIDE_DECAY, enabled=(mode == "on"))
    if ea is not None:
        step = dataclasses.replace(step, ea_kj_mol=float(ea))
    ps.NF_HYDROXIDE_DECAY = step
    engine._B2_FIT_REPORT = Path(report)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    rows = []
    for b in score_panel()["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"), "within": bool(r.get("within_band")),
                         "predicted": r.get("predicted"), "measured": r.get("measured")})
    out = {"mode": mode, "report": str(report), "ea_kj_mol": step.ea_kj_mol, "rows": rows,
           "mikami_k": {"8.0": _mikami_k(8.0), "7.0": _mikami_k(7.0)}}
    if ea is None:
        ind = score_directional()["summary"]["independent"]["total"]
        out["directional_independent"] = {"agree": ind["agree"], "evaluable": ind["evaluable"]}
    return out


def _cost_child(mode: str, report: Path) -> Dict[str, Any]:
    """The 52-row objective (B60's quick setting) at a report's vector, with the step switched."""
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))
    from src.kinetic_core import parameters_sulfur as ps

    import generate_kinetic_core_b60_fit as B60
    from generate_kinetic_core_b8_laplace import frozen_vector

    ps.NF_HYDROXIDE_DECAY = dataclasses.replace(ps.NF_HYDROXIDE_DECAY, enabled=(mode == "on"))
    r = B60.residual_vector(frozen_vector(json.loads(Path(report).read_text(encoding="utf-8"))), True)
    return {"mode": mode, "report": str(report), "cost": 0.5 * float(r @ r)}


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


def _bootstrap_se_median(values: List[float], n: int = 2000, seed: int = 61) -> float:
    if len(values) < 2:
        return 0.0
    rng = random.Random(seed)
    meds = [statistics.median(rng.choices(values, k=len(values))) for _ in range(n)]
    return statistics.pstdev(meds)


def _zhou(s: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [r for r in s["rows"] if r["benchmark"].startswith(ZHOU_PREFIX) and r["compound"] == MFT]


def _hold(s: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]


def _hof(s: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [r for r in s["rows"] if r["benchmark"] in T4_BUNDLES]


def _t2(base: Dict[str, Any], arm: Dict[str, Any]) -> Dict[str, Any]:
    b, a = _summary(_zhou(base)), _summary(_zhou(arm))
    drop = b["median_abs_log10"] - a["median_abs_log10"]
    return {"baseline": b, "arm": a, "drop_dex": drop, "pass": bool(len(_zhou(arm)) == 4 and drop >= 0.3),
            "rows": {r["benchmark"]: r["fold"] for r in _zhou(arm)}}


def _tests(base: Dict[str, Any], arm: Dict[str, Any]) -> Dict[str, Any]:
    k8, k7 = arm["mikami_k"]["8.0"], arm["mikami_k"]["7.0"]
    t1 = {"k_ph8": k8, "k_ph7": k7, "dex_off_ph8": math.log10(k8 / MIKAMI_K),
          "pass": bool(abs(math.log10(k8 / MIKAMI_K)) <= 0.1 and k7 >= MIKAMI_PH7_FLOOR)}
    tol3 = max(0.1, _bootstrap_se_median(_abs_logs(_hold(base))))
    tol4 = max(0.1, _bootstrap_se_median(_abs_logs(_hof(base))))
    t3 = {"baseline": _summary(_hold(base)), "arm": _summary(_hold(arm)), "tolerance_dex": tol3}
    t3["pass"] = bool(t3["arm"]["within_3x"] >= t3["baseline"]["within_3x"] - 2
                      and t3["arm"]["median_abs_log10"] <= t3["baseline"]["median_abs_log10"] + tol3)
    t4 = {"baseline": _summary(_hof(base)), "arm": _summary(_hof(arm)), "tolerance_dex": tol4}
    t4["pass"] = bool(t4["arm"]["median_abs_log10"] <= t4["baseline"]["median_abs_log10"] + tol4)
    t5 = {"baseline": base["directional_independent"], "arm": arm["directional_independent"]}
    t5["pass"] = bool(arm["directional_independent"]["agree"] >= base["directional_independent"]["agree"] - 3)
    t2 = _t2(base, arm)
    return {"T1_measurement": t1, "T2_zhou_mft": t2, "T3_holdout_sulfur": t3, "T4_hofmann_ph5": t4,
            "T5_directional": t5, "union_panel": {"baseline": _summary(base["rows"]), "arm": _summary(arm["rows"])},
            "pass": bool(t1["pass"] and t2["pass"] and t3["pass"] and t4["pass"] and t5["pass"])}


def _recorded(scored: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    def fold(s, bench, compound):
        return next((r["fold"] for r in s["rows"] if r["benchmark"] == bench and r["compound"] == compound), None)

    out: Dict[str, Any] = {"hofmann_sweep": {}, "zhang_yiltirak_mft": {}}
    for arm, s in scored.items():
        out["hofmann_sweep"][arm] = {f"pH {ph:g}": {c: fold(s, b, c) for c in (MFT, "2-Furfurylthiol (FFT)")}
                                     for ph, b in SWEEP.items()}
        out["zhang_yiltirak_mft"][arm] = {r["benchmark"]: r["fold"] for r in s["rows"] if r["compound"] == MFT and (
            "Zhang2025b" in r["benchmark"] or "Yiltirak" in r["benchmark"])}
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", nargs=2, metavar=("MODE", "REPORT"))
    parser.add_argument("--ea", type=float, default=None)
    parser.add_argument("--cost", nargs=2, metavar=("MODE", "REPORT"))
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child[0], Path(args.child[1]), args.ea)))
        return 0
    if args.cost:
        print(json.dumps(_cost_child(args.cost[0], Path(args.cost[1]))))
        return 0
    from src.kinetic_core import ea_class_prior

    scored = {arm: _run("--child", mode, str(report)) for arm, (mode, report) in ARMS.items()}
    arms = {arm: _tests(scored["baseline"], scored[arm]) for arm in ("A_refit", "B_frozen_b60")}
    if arms["A_refit"]["pass"]:
        verdict = "SHIP (A): the B61 refit"
    elif arms["B_frozen_b60"]["pass"]:
        verdict = "SHIP (B): B60's fit with the hydroxide-catalysed norfuraneol decay"
    else:
        verdict = "DO NOT SHIP"
    carb = ea_class_prior.class_values()[ea_class_prior.CARBOHYDRATE]
    sensitivity = {}
    for q in (0.1, 0.9):
        ea = ea_class_prior.empirical_quantile(carb, q)
        s = _run("--child", "on", str(B60_REPORT), "--ea", f"{ea:.6f}")
        sensitivity[f"q{int(q * 100)}"] = {"ea_kj_mol": ea, "T2": _t2(scored["baseline"], s)}
    costs = {arm: _run("--cost", mode, str(report))["cost"] for arm, (mode, report) in ARMS.items()}
    recorded = _recorded(scored)
    recorded["fit_cost"] = {arm: {"cost": c, "rise_over_b60": c - COST_B60, "within_admissible_rise": c - COST_B60 <= ADMISSIBLE_RISE}
                            for arm, c in costs.items()}
    recorded["ea_sensitivity_b60_vector_step_on"] = sensitivity
    payload = {
        "artifact": "kinetic_core_b61_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b61_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b61_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "sulfur.py", ROOT / "src" / "kinetic_core" / "parameters_sulfur.py",
                    B60_REPORT, B61_REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b61_prereg.md"),
        "arms": arms, "recorded": recorded, "verdict": verdict,
        "rule": "SHIP (A) if (A) passes T1-T5; else SHIP (B) if (B) does; else DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    a, b = arms["A_refit"], arms["B_frozen_b60"]
    names = ("T1_measurement", "T2_zhou_mft", "T3_holdout_sulfur", "T4_hofmann_ph5", "T5_directional")
    print(f"{verdict} | A: " + " ".join(f"{n[:2]} {a[n]['pass']}" for n in names)
          + " | B: " + " ".join(f"{n[:2]} {b[n]['pass']}" for n in names))
    return 0


def _fmt(s: Dict[str, Any]) -> str:
    m = s["median_abs_log10"]
    return f"{s['within_3x']}/{s['rows']} ({m:.2f} dex)" if m is not None else f"{s['within_3x']}/{s['rows']}"


def render(payload: Dict[str, Any]) -> str:
    a, b = payload["arms"]["A_refit"], payload["arms"]["B_frozen_b60"]
    rec = payload["recorded"]
    lines = ["# Wave B61 ship rule — norfuraneol's hydroxide-catalysed decay (Mikami 2017)", "",
             f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
             f"## Verdict: **{payload['verdict']}**", "",
             "| test | baseline (B60, step off) | (A) B61 refit | pass | (B) B60, step on | pass |",
             "| --- | --- | --- | --- | --- | --- |"]
    ta, tb = a["T1_measurement"], b["T1_measurement"]
    lines.append(f"| T1 Mikami pot: k at label pH 8 (0.032) / pH 7 (>= 1.8e-3), 1/min | — | {ta['k_ph8']:.4f} / {ta['k_ph7']:.2e} | "
                 f"{ta['pass']} | {tb['k_ph8']:.4f} / {tb['k_ph7']:.2e} | {tb['pass']} |")
    za, zb = a["T2_zhou_mft"], b["T2_zhou_mft"]
    lines.append(f"| T2 Zhou 2023c MFT, median dex (must fall >= 0.3) | {za['baseline']['median_abs_log10']:.2f} | "
                 f"{za['arm']['median_abs_log10']:.2f} (-{za['drop_dex']:.2f}) | {za['pass']} | "
                 f"{zb['arm']['median_abs_log10']:.2f} (-{zb['drop_dex']:.2f}) | {zb['pass']} |")
    lines.append(f"| T3 hold-out sulfur rows (tolerance {a['T3_holdout_sulfur']['tolerance_dex']:.2f} dex) | {_fmt(a['T3_holdout_sulfur']['baseline'])} | "
                 f"{_fmt(a['T3_holdout_sulfur']['arm'])} | {a['T3_holdout_sulfur']['pass']} | {_fmt(b['T3_holdout_sulfur']['arm'])} | {b['T3_holdout_sulfur']['pass']} |")
    lines.append(f"| T4 Hofmann ribose + xylose pH 5 (tolerance {a['T4_hofmann_ph5']['tolerance_dex']:.2f} dex) | {_fmt(a['T4_hofmann_ph5']['baseline'])} | "
                 f"{_fmt(a['T4_hofmann_ph5']['arm'])} | {a['T4_hofmann_ph5']['pass']} | {_fmt(b['T4_hofmann_ph5']['arm'])} | {b['T4_hofmann_ph5']['pass']} |")
    d0, da, db = a["T5_directional"]["baseline"], a["T5_directional"]["arm"], b["T5_directional"]["arm"]
    lines.append(f"| T5 directional, independent agree | {d0['agree']}/{d0['evaluable']} | {da['agree']}/{da['evaluable']} | "
                 f"{a['T5_directional']['pass']} | {db['agree']}/{db['evaluable']} | {b['T5_directional']['pass']} |")
    lines.append(f"| union panel within 3x (recorded) | {_fmt(a['union_panel']['baseline'])} | {_fmt(a['union_panel']['arm'])} | — | "
                 f"{_fmt(b['union_panel']['arm'])} | — |")
    lines += ["", "## Recorded, not gating", "", "| arm | fit cost | rise over B60 (admissible <= 3.59) |", "| --- | --- | --- |"]
    for arm, c in rec["fit_cost"].items():
        lines.append(f"| {arm} | {c['cost']:.3f} | {c['rise_over_b60']:+.3f} ({'inside' if c['within_admissible_rise'] else 'OUTSIDE'}) |")
    lines += ["", "Hofmann ribose + Cys sweep, fold error MFT / FFT:", "", "| arm | pH 3 | pH 5 | pH 7 |", "| --- | --- | --- | --- |"]
    for arm, sw in rec["hofmann_sweep"].items():
        cells = [" / ".join("-" if v is None else f"{v:.3g}" for v in sw[p].values()) for p in ("pH 3", "pH 5", "pH 7")]
        lines.append(f"| {arm} | " + " | ".join(cells) + " |")
    lines += ["", "T2 at the carbohydrate class's 10th / 90th percentile barrier (B60's vector, step on):", ""]
    for q, v in rec["ea_sensitivity_b60_vector_step_on"].items():
        lines.append(f"- {q}: Ea {v['ea_kj_mol']:.1f} kJ/mol, Zhou MFT median {v['T2']['arm']['median_abs_log10']:.2f} dex "
                     f"(drop {v['T2']['drop_dex']:.2f}, pass {v['T2']['pass']})")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
