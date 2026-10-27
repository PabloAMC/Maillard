#!/usr/bin/env python
"""
Wave B47 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b47_prereg.md` sec. 5. B9 and B47 are each scored in a child
process whose engine reads that wave's fit report (``engine._B2_FIT_REPORT``), so every number below
comes from the code path a user calls; the engine on disk is not switched. Writes
`results/validation/kinetic_core_b47_ship_rule.{json,md}`.

  T1  recorded: B47's cost on its 52 rows, and B9's vector on the same 52 rows
  T2  Maillard-path hold-out sulfur rows: within-3x not lower than B9's by more than one row, and the
      median |log10 error| not worse by more than 0.1 dex
  T3  the four Hofmann pH-5 bundles (out of sample since B9): within-3x not lower than B9's, median
      |log10 error| not worse by more than 0.1 dex
  T4  the directional panel, independent claims agreeing: not fewer than B9's
  T5  recorded: Laplace identified coordinates at the B47 optimum
Ship rule: SHIP if T2, T3 and T4 hold.
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
B47_REPORT = V / "kinetic_core_b47_fit_report.json"
B47_LAPLACE = V / "kinetic_core_b47_laplace_covariance.json"
B9_LAPLACE = V / "kinetic_core_b9_laplace_covariance.json"
OUT = V / "kinetic_core_b47_ship_rule.json"
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
    import numpy as np

    import generate_kinetic_core_b47_fit as B47
    from generate_kinetic_core_b8_laplace import frozen_vector

    r47 = B47.residual_vector(frozen_vector(_read(B47_REPORT)), True)
    r9 = B47.residual_vector(frozen_vector(_read(B9_REPORT)), True)
    return {"rows": int(len(r47)), "b47_cost": 0.5 * float(np.dot(r47, r47)),
            "b9_vector_on_b47_rows_cost": 0.5 * float(np.dot(r9, r9))}


def _t5() -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for tag, path in (("b9", B9_LAPLACE), ("b47", B47_LAPLACE)):
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
    for p in (B9_REPORT, B47_REPORT):
        if not p.exists():
            raise SystemExit(f"{p} missing")
    s9, s47 = _score(B9_REPORT), _score(B47_REPORT)

    def holdout_sulfur(s):
        return [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]

    def hofmann(s):
        return [r for r in s["rows"] if r["benchmark"] in HOFMANN_PH5]

    t2 = {"b9": _summary(holdout_sulfur(s9)), "b47": _summary(holdout_sulfur(s47))}
    t2["pass"] = (t2["b47"]["within_3x"] >= t2["b9"]["within_3x"] - 1
                  and t2["b47"]["median_abs_log10"] <= t2["b9"]["median_abs_log10"] + 0.1)
    t3 = {"b9": _summary(hofmann(s9)), "b47": _summary(hofmann(s47))}
    t3["pass"] = (t3["b47"]["within_3x"] >= t3["b9"]["within_3x"]
                  and t3["b47"]["median_abs_log10"] <= t3["b9"]["median_abs_log10"] + 0.1)
    t4 = {"b9": s9["directional_independent"], "b47": s47["directional_independent"]}
    t4["pass"] = t4["b47"]["agree"] >= t4["b9"]["agree"]
    all9, all47 = _summary(s9["rows"]), _summary(s47["rows"])
    verdict = "SHIP" if (t2["pass"] and t3["pass"] and t4["pass"]) else "DO NOT SHIP"
    payload = {
        "artifact": "kinetic_core_b47_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b47_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b47_ship_rule.py",
            inputs=[p for p in (B9_REPORT, B47_REPORT, B9_LAPLACE, B47_LAPLACE) if p.exists()]),
        "prereg": data_paths.rel(V / "kinetic_core_b47_prereg.md"),
        "T1_cost": _t1(), "T2_holdout_sulfur": t2, "T3_hofmann_pH5": t3, "T4_directional": t4, "T5_laplace": _t5(),
        "union_panel": {"b9": {k: v for k, v in all9.items() if k != "detail"},
                        "b47": {k: v for k, v in all47.items() if k != "detail"}},
        "verdict": verdict,
        "rule": "SHIP if T2 (hold-out sulfur rows), T3 (Hofmann pH-5 bundles) and T4 (directional) hold",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: T2 {t2['b9']['within_3x']}/{t2['b9']['rows']} -> {t2['b47']['within_3x']}/{t2['b47']['rows']} "
          f"(median {t2['b9']['median_abs_log10']:.2f} -> {t2['b47']['median_abs_log10']:.2f} dex); "
          f"T3 {t3['b9']['within_3x']}/{t3['b9']['rows']} -> {t3['b47']['within_3x']}/{t3['b47']['rows']} "
          f"(median {t3['b9']['median_abs_log10']:.2f} -> {t3['b47']['median_abs_log10']:.2f}); "
          f"T4 {t4['b9']['agree']}/{t4['b9']['evaluable']} -> {t4['b47']['agree']}/{t4['b47']['evaluable']}; "
          f"union {all9['within_3x']}/{all9['rows']} -> {all47['within_3x']}/{all47['rows']}")
    return 0


def _rows_md(title: str, t: Dict[str, Any]) -> List[str]:
    out = [f"### {title}", "", "| bundle | compound | B9 fold | B47 fold |", "| --- | --- | ---: | ---: |"]
    b47 = {(r["benchmark"], r["compound"]): r for r in t["b47"]["detail"]}
    for r in t["b9"]["detail"]:
        o = b47.get((r["benchmark"], r["compound"]), {})
        out.append(f"| {r['benchmark']} | {r['compound']} | {r['fold']:.3g} | {o.get('fold', float('nan')):.3g} |")
    return out + [""]


def render(payload: Dict[str, Any]) -> str:
    t1, t2, t3, t4, t5 = (payload[k] for k in ("T1_cost", "T2_holdout_sulfur", "T3_hofmann_pH5", "T4_directional", "T5_laplace"))
    u = payload["union_panel"]
    lines = [
        "# Wave B47 ship rule — the sulfur objective, corrected against its PDFs", "",
        f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
        f"## Verdict: **{payload['verdict']}**", "",
        "| test | B9 | B47 | pass |", "| --- | --- | --- | --- |",
        f"| T1 cost on B47's {t1['rows']} rows (recorded) | {t1['b9_vector_on_b47_rows_cost']:.2f} | {t1['b47_cost']:.2f} | — |",
        f"| T2 hold-out sulfur rows within 3x | {t2['b9']['within_3x']}/{t2['b9']['rows']} (median {t2['b9']['median_abs_log10']:.2f} dex) "
        f"| {t2['b47']['within_3x']}/{t2['b47']['rows']} (median {t2['b47']['median_abs_log10']:.2f} dex) | {t2['pass']} |",
        f"| T3 Hofmann pH-5 rows within 3x | {t3['b9']['within_3x']}/{t3['b9']['rows']} (median {t3['b9']['median_abs_log10']:.2f} dex) "
        f"| {t3['b47']['within_3x']}/{t3['b47']['rows']} (median {t3['b47']['median_abs_log10']:.2f} dex) | {t3['pass']} |",
        f"| T4 directional, independent agree | {t4['b9']['agree']}/{t4['b9']['evaluable']} | {t4['b47']['agree']}/{t4['b47']['evaluable']} | {t4['pass']} |",
        f"| T5 Laplace identified (recorded) | {t5['b9'] if isinstance(t5['b9'], str) else str(t5['b9']['identified']) + '/' + str(t5['b9']['free'])} "
        f"| {t5['b47'] if isinstance(t5['b47'], str) else str(t5['b47']['identified']) + '/' + str(t5['b47']['free'])} | — |",
        f"| union panel within 3x (recorded) | {u['b9']['within_3x']}/{u['b9']['rows']} | {u['b47']['within_3x']}/{u['b47']['rows']} | — |",
        "",
    ]
    return "\n".join(lines + _rows_md("T2 rows", t2) + _rows_md("T3 rows", t3))


if __name__ == "__main__":
    raise SystemExit(main())
