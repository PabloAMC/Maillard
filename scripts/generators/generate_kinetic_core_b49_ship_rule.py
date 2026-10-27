#!/usr/bin/env python
"""
Wave B49 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b49_prereg.md` sec. 4. B9 and B49 are each scored in a child
process whose engine reads that wave's fit report (``engine._B2_FIT_REPORT``), so every number below
comes from the code path a user calls; the engine on disk is not switched. Writes
`results/validation/kinetic_core_b49_ship_rule.{json,md}`.

  T1  the mechanism, held out: Whitfield 2001's pot (NF + cysteine 50/50 mM, 0.5 M phosphate pH 6.5,
      140 C, 60 min): predicted MFT <= 0.003 mol %, and MFT falls >= 10x from the pH 4.5 twin
  T2  Maillard-path hold-out sulfur rows: within-3x not lower than B9's by more than one row, and the
      median |log10 error| not worse by more than 0.1 dex
  T3  the four Hofmann pH-5 bundles (out of sample since B9): within-3x not lower than B9's, median
      |log10 error| not worse by more than 0.1 dex
  T4  the directional panel, independent claims agreeing: not fewer than B9's
  T5  recorded: Laplace identified coordinates; the MFT sink at 140 C; B49's cost
Ship rule: SHIP if T1, T2, T3 and T4 hold.
"""
from __future__ import annotations

import argparse
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
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
B9_REPORT = V / "kinetic_core_b9_fit_report.json"
B49_REPORT = V / "kinetic_core_b49_fit_report.json"
B49_LAPLACE = V / "kinetic_core_b49_laplace_covariance.json"
B9_LAPLACE = V / "kinetic_core_b9_laplace_covariance.json"
OUT = V / "kinetic_core_b49_ship_rule.json"
HOFMANN_PH5 = tuple(
    f"hofmann1998_{sugar}_cysteine_145C_20min_pH5" for sugar in ("ribose", "xylose", "glucose", "fructose")
)


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def _child(report: Path) -> Dict[str, Any]:
    """Score the union panel and the directional panel with the engine reading ``report``."""
    from src.kinetic_core import engine

    engine._B2_FIT_REPORT = Path(report)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    panel = score_panel()
    rows = []
    for b in panel["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"),
                         "within": bool(r.get("within_band"))})
    directional = score_directional()
    ind = directional["summary"]["independent"]["total"]
    return {"report": str(report), "rows": rows,
            "directional_independent": {"agree": ind["agree"], "evaluable": ind["evaluable"]}}


def _score(report: Path) -> Dict[str, Any]:
    done = subprocess.run([sys.executable, __file__, "--child", str(report)], cwd=ROOT,
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"child for {report} failed:\n{done.stderr[-3000:]}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def _summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    folds = [r["fold"] for r in rows if r["fold"] is not None and math.isfinite(r["fold"]) and r["fold"] > 0]
    return {"rows": len(rows), "within_3x": sum(r["within"] for r in rows),
            "median_abs_log10": statistics.median(abs(math.log10(f)) for f in folds) if folds else None,
            "detail": rows}


def _t1() -> Dict[str, Any]:
    """The pH-6.5 hold-out, integrated by the fit harness under B9's and B49's vectors."""
    import generate_kinetic_core_b49_fit as B49
    from generate_kinetic_core_b8_laplace import frozen_vector
    from src.kinetic_core.ph_state import BufferSpec

    B23 = B49.B23

    def mft_molpct(x, build, ph):
        params = build(x)
        _f, _e, _d, drift = B23.unpack(x)
        buf = BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="Whitfield & Mottram 1999 / 2001")
        run = B23.integrate_sulfur(params, 140.0 + B23.CELSIUS, {"NF": 50.0, "Cys": 50.0}, np.array([0.0, 60.0]),
                                   ph=ph, buffer_spec=buf, ph_drift=drift, ph_nodes=B23.PH_NODES_CAREFUL,
                                   ph_iterations=B23.PH_ITERS_CAREFUL, rtol=1e-8, atol=1e-14)
        return 100.0 * run.final("MFT") / 50.0

    import numpy as np  # noqa: F811

    x49 = frozen_vector(_read(B49_REPORT))
    x9 = frozen_vector(_read(B9_REPORT))
    out = {}
    for tag, x, build in (("b9", x9, B49._B23_BUILD_ORIGINAL), ("b49", x49, B49.build_parameters)):
        lo, hi = mft_molpct(x, build, 4.5), mft_molpct(x, build, 6.5)
        out[tag] = {"mft_molpct_pH4_5": lo, "mft_molpct_pH6_5": hi, "fall_factor": lo / hi if hi > 0 else float("inf")}
    out["measured"] = {"mft_molpct_pH4_5": 0.150, "mft_molpct_pH6_5": "< 0.001 (not detected)", "fall_factor": ">= 150"}
    r = out["b49"]
    out["pass"] = bool(r["mft_molpct_pH6_5"] <= 0.003 and r["fall_factor"] >= 10.0)
    import numpy as np  # noqa: F811
    r49 = B49.residual_vector(x49, True)
    out["b49_cost"] = 0.5 * float(np.dot(r49, r49))
    return out


def _t5() -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for tag, path in (("b9", B9_LAPLACE), ("b49", B49_LAPLACE)):
        if path.exists():
            lap = _read(path)
            out[tag] = {"identified": int(sum(bool(v) for v in lap["identified"])), "free": len(lap["identified"]),
                        "chi2_red": lap.get("chi2_red")}
        else:
            out[tag] = "no laplace artifact"
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", default=None)
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(Path(args.child))))
        return 0
    for p in (B9_REPORT, B49_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    s9, s49 = _score(B9_REPORT), _score(B49_REPORT)

    def holdout_sulfur(s):
        return [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]

    def hofmann(s):
        return [r for r in s["rows"] if r["benchmark"] in HOFMANN_PH5]

    t2 = {"b9": _summary(holdout_sulfur(s9)), "b49": _summary(holdout_sulfur(s49))}
    t2["pass"] = (t2["b49"]["within_3x"] >= t2["b9"]["within_3x"] - 1
                  and t2["b49"]["median_abs_log10"] <= t2["b9"]["median_abs_log10"] + 0.1)
    t3 = {"b9": _summary(hofmann(s9)), "b49": _summary(hofmann(s49))}
    t3["pass"] = (t3["b49"]["within_3x"] >= t3["b9"]["within_3x"]
                  and t3["b49"]["median_abs_log10"] <= t3["b9"]["median_abs_log10"] + 0.1)
    t4 = {"b9": s9["directional_independent"], "b49": s49["directional_independent"]}
    t4["pass"] = t4["b49"]["agree"] >= t4["b9"]["agree"]
    all9, all49 = _summary(s9["rows"]), _summary(s49["rows"])
    t1 = _t1()
    verdict = "SHIP" if (t1["pass"] and t2["pass"] and t3["pass"] and t4["pass"]) else "DO NOT SHIP"
    payload = {
        "artifact": "kinetic_core_b49_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b49_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b49_ship_rule.py",
            inputs=[p for p in (B9_REPORT, B49_REPORT, B9_LAPLACE, B49_LAPLACE) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b49_prereg.md"),
        "T1_mechanism_pH6_5": t1, "T2_holdout_sulfur": t2, "T3_hofmann_pH5": t3, "T4_directional": t4, "T5_laplace": _t5(),
        "union_panel": {"b9": {k: v for k, v in all9.items() if k != "detail"},
                        "b49": {k: v for k, v in all49.items() if k != "detail"}},
        "verdict": verdict,
        "rule": "SHIP if T1 (pH 6.5 mechanism, held out), T2 (hold-out sulfur rows), T3 (Hofmann pH-5 bundles) and T4 (directional) hold",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: T1 {t1['pass']} (B9 fall {t1['b9']['fall_factor']:.3g}x, B49 fall {t1['b49']['fall_factor']:.3g}x, "
          f"B49 MFT pH 6.5 {t1['b49']['mft_molpct_pH6_5']:.3g} mol%); T2 {t2['b9']['within_3x']}/{t2['b9']['rows']} -> {t2['b49']['within_3x']}/{t2['b49']['rows']} "
          f"(median {t2['b9']['median_abs_log10']:.2f} -> {t2['b49']['median_abs_log10']:.2f} dex); "
          f"T3 {t3['b9']['within_3x']}/{t3['b9']['rows']} -> {t3['b49']['within_3x']}/{t3['b49']['rows']} "
          f"(median {t3['b9']['median_abs_log10']:.2f} -> {t3['b49']['median_abs_log10']:.2f}); "
          f"T4 {t4['b9']['agree']}/{t4['b9']['evaluable']} -> {t4['b49']['agree']}/{t4['b49']['evaluable']}; "
          f"union {all9['within_3x']}/{all9['rows']} -> {all49['within_3x']}/{all49['rows']}")
    return 0


def _rows_md(title: str, t: Dict[str, Any]) -> List[str]:
    out = [f"### {title}", "", "| bundle | compound | B9 fold | B49 fold |", "| --- | --- | ---: | ---: |"]
    b49 = {(r["benchmark"], r["compound"]): r for r in t["b49"]["detail"]}
    for r in t["b9"]["detail"]:
        o = b49.get((r["benchmark"], r["compound"]), {})
        out.append(f"| {r['benchmark']} | {r['compound']} | {r['fold']:.3g} | {o.get('fold', float('nan')):.3g} |")
    return out + [""]


def render(payload: Dict[str, Any]) -> str:
    t1, t2, t3, t4, t5 = (payload[k] for k in ("T1_mechanism_pH6_5", "T2_holdout_sulfur", "T3_hofmann_pH5", "T4_directional", "T5_laplace"))
    u = payload["union_panel"]
    lines = [
        "# Wave B49 ship rule — norfuraneol's other sulfide products", "",
        f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
        f"## Verdict: **{payload['verdict']}**", "",
        "| test | B9 | B49 | pass |", "| --- | --- | --- | --- |",
        f"| T1 NF + Cys, 140 C: MFT fall pH 4.5 -> 6.5 (measured >= 150x; MFT at 6.5 < 0.001 mol %) "
        f"| {t1['b9']['fall_factor']:.3g}x ({t1['b9']['mft_molpct_pH6_5']:.3g} mol %) "
        f"| {t1['b49']['fall_factor']:.3g}x ({t1['b49']['mft_molpct_pH6_5']:.3g} mol %) | {t1['pass']} |",
        f"| T2 hold-out sulfur rows within 3x | {t2['b9']['within_3x']}/{t2['b9']['rows']} (median {t2['b9']['median_abs_log10']:.2f} dex) "
        f"| {t2['b49']['within_3x']}/{t2['b49']['rows']} (median {t2['b49']['median_abs_log10']:.2f} dex) | {t2['pass']} |",
        f"| T3 Hofmann pH-5 rows within 3x | {t3['b9']['within_3x']}/{t3['b9']['rows']} (median {t3['b9']['median_abs_log10']:.2f} dex) "
        f"| {t3['b49']['within_3x']}/{t3['b49']['rows']} (median {t3['b49']['median_abs_log10']:.2f} dex) | {t3['pass']} |",
        f"| T4 directional, independent agree | {t4['b9']['agree']}/{t4['b9']['evaluable']} | {t4['b49']['agree']}/{t4['b49']['evaluable']} | {t4['pass']} |",
        f"| T5 Laplace identified (recorded) | {t5['b9'] if isinstance(t5['b9'], str) else str(t5['b9']['identified']) + '/' + str(t5['b9']['free'])} "
        f"| {t5['b49'] if isinstance(t5['b49'], str) else str(t5['b49']['identified']) + '/' + str(t5['b49']['free'])} | — |",
        f"| union panel within 3x (recorded) | {u['b9']['within_3x']}/{u['b9']['rows']} | {u['b49']['within_3x']}/{u['b49']['rows']} | — |",
        "",
    ]
    return "\n".join(lines + _rows_md("T2 rows", t2) + _rows_md("T3 rows", t3))


if __name__ == "__main__":
    raise SystemExit(main())
