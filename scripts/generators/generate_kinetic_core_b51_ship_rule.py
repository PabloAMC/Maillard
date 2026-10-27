#!/usr/bin/env python
"""
Wave B51 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b51_prereg.md` sec. 5. Three arms, each scored in its own child process
whose engine reads that arm's ``ph_state.B51`` switches and fit report, so every number comes from the
code path a user calls:

  baseline  B9's frozen vector, switches OFF  (what ships, plus the B51 data corrections)
  (A)       the B51 refit, switches ON
  (B)       B9's frozen vector, switches ON   (corrected physics, no refit)

  T1  recorded: cost of every arm on B47's 52 rows; per-switch attribution of B9's vector on its 54 rows
  T2  Maillard-path hold-out sulfur rows: within-3x not lower than the baseline's by more than one row,
      median |log10 error| not worse by more than 0.1 dex
  T3  Hofmann ribose pH 5 (trust loop) + xylose pH 5 (hold-out): within-3x not lower, median |log10 error|
      not worse by more than 0.1 dex
  T4  directional panel, independent claims agreeing: not fewer than the baseline's
  T5  recorded: the pH-axis directional subset; Laplace identified coordinates for (A) when computed
Verdict: SHIP (A) if (A) passes T2-T4; else SHIP (B) if (B) does; else DO NOT SHIP.
Writes `results/validation/kinetic_core_b51_ship_rule.{json,md}`.
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
B51_REPORT = V / "kinetic_core_b51_fit_report.json"
B51_LAPLACE = V / "kinetic_core_b51_laplace_covariance.json"
OUT = V / "kinetic_core_b51_ship_rule.json"
T3_BUNDLES = ("hofmann1998_ribose_cysteine_145C_20min_pH5", "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5")
SWITCHES = ("heat_capacity", "ionic_strength", "cysteine_gkl", "h2s_suleimenov", "label_scale_tables",
            "sid_across_segments")
ARMS = {
    "baseline": ("off", B9_REPORT),
    "A_refit": ("on", B51_REPORT),
    "B_frozen_b9": ("on", B9_REPORT),
}


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def _set_switches(mode: str) -> None:
    from src.kinetic_core import ph_state

    ph_state.B51 = ph_state.B51_ALL_ON if mode == "on" else ph_state.PhStateCorrections()


def _child(mode: str, report: Path) -> Dict[str, Any]:
    _set_switches(mode)
    from src.kinetic_core import engine

    engine._B2_FIT_REPORT = Path(report)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.scoring import score_panel

    rows = []
    for b in score_panel()["benchmarks"]:
        for r in b["compounds"]:
            rows.append({"benchmark": b["benchmark_id"], "panel": b["panel"], "lane": r.get("lane"),
                         "compound": r.get("compound"), "fold": r.get("fold_error"),
                         "within": bool(r.get("within_band"))})
    summary = score_directional()["summary"]["independent"]
    ph = summary["by_category"].get("ph", {})
    return {"mode": mode, "report": str(report), "rows": rows,
            "directional_independent": {"agree": summary["total"]["agree"],
                                        "evaluable": summary["total"]["evaluable"]},
            "directional_ph": {"agree": ph.get("agree"), "evaluable": ph.get("evaluable"),
                               "misses": ph.get("misses")}}


def _attribution() -> Dict[str, Any]:
    """B9's frozen vector on B9's own 54 rows, with each switch on alone (B47 is NOT imported here)."""
    import numpy as np

    import generate_kinetic_core_b9_fit as B9
    from generate_kinetic_core_b8_laplace import frozen_vector
    from src.kinetic_core import ph_state

    x = frozen_vector(_read(B9_REPORT))
    out = {}
    for name in ("none",) + SWITCHES + ("all",):
        if name == "none":
            ph_state.B51 = ph_state.PhStateCorrections()
        elif name == "all":
            ph_state.B51 = ph_state.B51_ALL_ON
        else:
            ph_state.B51 = ph_state.PhStateCorrections(**{name: True})
        r = B9.residual_vector(x, True)
        out[name] = 0.5 * float(np.dot(r, r))
    return out


def _costs_on_b47_rows() -> Dict[str, Any]:
    import numpy as np

    import generate_kinetic_core_b47_fit as B47
    from generate_kinetic_core_b8_laplace import frozen_vector
    from src.kinetic_core import ph_state

    x9 = frozen_vector(_read(B9_REPORT))
    xa = frozen_vector(_read(B51_REPORT))
    out = {}
    for arm, mode, x in (("baseline", "off", x9), ("A_refit", "on", xa), ("B_frozen_b9", "on", x9)):
        ph_state.B51 = ph_state.B51_ALL_ON if mode == "on" else ph_state.PhStateCorrections()
        r = B47.residual_vector(x, True)
        out[arm] = 0.5 * float(np.dot(r, r))
    return out


def _run(*args: str) -> Dict[str, Any]:
    done = subprocess.run([sys.executable, __file__, *args], cwd=ROOT, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"child {args} failed:\n{done.stderr[-3000:]}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def _summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    folds = [r["fold"] for r in rows if r["fold"] is not None and math.isfinite(r["fold"]) and r["fold"] > 0]
    return {"rows": len(rows), "within_3x": sum(r["within"] for r in rows),
            "median_abs_log10": statistics.median(abs(math.log10(f)) for f in folds) if folds else None}


def _tests(base: Dict[str, Any], arm: Dict[str, Any]) -> Dict[str, Any]:
    def t2_rows(s):
        return [r for r in s["rows"] if r["panel"] == "maillard_path_holdout" and r["lane"] == "sulfur"]

    def t3_rows(s):
        return [r for r in s["rows"] if r["benchmark"] in T3_BUNDLES]

    t2 = {"baseline": _summary(t2_rows(base)), "arm": _summary(t2_rows(arm))}
    t2["pass"] = bool(t2["arm"]["within_3x"] >= t2["baseline"]["within_3x"] - 1
                      and t2["arm"]["median_abs_log10"] <= t2["baseline"]["median_abs_log10"] + 0.1)
    t3 = {"baseline": _summary(t3_rows(base)), "arm": _summary(t3_rows(arm))}
    t3["pass"] = bool(t3["arm"]["within_3x"] >= t3["baseline"]["within_3x"]
                      and t3["arm"]["median_abs_log10"] <= t3["baseline"]["median_abs_log10"] + 0.1)
    t4 = {"baseline": base["directional_independent"], "arm": arm["directional_independent"]}
    t4["pass"] = bool(arm["directional_independent"]["agree"] >= base["directional_independent"]["agree"])
    return {"T2_holdout_sulfur": t2, "T3_hofmann_ph5": t3, "T4_directional": t4,
            "T5_directional_ph": {"baseline": base["directional_ph"], "arm": arm["directional_ph"]},
            "union_panel": {"baseline": _summary(base["rows"]), "arm": _summary(arm["rows"])},
            "pass": bool(t2["pass"] and t3["pass"] and t4["pass"])}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", nargs=2, metavar=("MODE", "REPORT"))
    parser.add_argument("--attribution", action="store_true")
    parser.add_argument("--costs", action="store_true")
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child[0], Path(args.child[1]))))
        return 0
    if args.attribution:
        print(json.dumps(_attribution()))
        return 0
    if args.costs:
        print(json.dumps(_costs_on_b47_rows()))
        return 0

    scored = {arm: _run("--child", mode, str(report)) for arm, (mode, report) in ARMS.items()}
    arms = {arm: _tests(scored["baseline"], scored[arm]) for arm in ("A_refit", "B_frozen_b9")}
    t1 = {"cost_on_b47_52_rows": _run("--costs"), "attribution_b9_vector_on_b9_54_rows": _run("--attribution")}
    laplace = None
    if B51_LAPLACE.exists():
        lap = _read(B51_LAPLACE)
        laplace = {"identified": lap.get("identified_count", lap.get("n_identified")), "free": lap.get("n_free")}
    if arms["A_refit"]["pass"]:
        verdict = "SHIP (A): the B51 refit"
    elif arms["B_frozen_b9"]["pass"]:
        verdict = "SHIP (B): B9's vector under the corrected pH state"
    else:
        verdict = "DO NOT SHIP"
    moved = []
    base_by = {(r["benchmark"], r["compound"]): r for r in scored["baseline"]["rows"]}
    for arm in ("A_refit", "B_frozen_b9"):
        for r in scored[arm]["rows"]:
            b = base_by.get((r["benchmark"], r["compound"]))
            if b and b["fold"] and r["fold"] and math.isfinite(b["fold"]) and math.isfinite(r["fold"]):
                d = abs(math.log10(r["fold"])) - abs(math.log10(b["fold"]))
                if abs(d) > 0.05:
                    moved.append({"arm": arm, "benchmark": r["benchmark"], "compound": r["compound"],
                                  "fold_baseline": b["fold"], "fold_arm": r["fold"], "delta_abs_log10": d})
    payload = {
        "artifact": "kinetic_core_b51_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b51_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b51_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "ph_state.py", ROOT / "src" / "kinetic_core" / "parameters_sulfur.py",
                    ROOT / "src" / "kinetic_core" / "sulfur.py", B9_REPORT, B51_REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b51_prereg.md"),
        "T1": t1, "arms": arms, "laplace_A": laplace, "rows_moved_over_0_05_dex": moved,
        "verdict": verdict,
        "rule": "SHIP (A) if (A) passes T2, T3, T4; else SHIP (B) if (B) passes them; else DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    a, b = arms["A_refit"], arms["B_frozen_b9"]
    print(f"{verdict} | A: T2 {a['T2_holdout_sulfur']['pass']} T3 {a['T3_hofmann_ph5']['pass']} T4 {a['T4_directional']['pass']}"
          f" | B: T2 {b['T2_holdout_sulfur']['pass']} T3 {b['T3_hofmann_ph5']['pass']} T4 {b['T4_directional']['pass']}"
          f" | costs {t1['cost_on_b47_52_rows']}")
    return 0


def _fmt(s: Dict[str, Any]) -> str:
    med = s["median_abs_log10"]
    return f"{s['within_3x']}/{s['rows']} (median {med:.2f} dex)" if med is not None else f"{s['within_3x']}/{s['rows']}"


def render(payload: Dict[str, Any]) -> str:
    t1, arms = payload["T1"], payload["arms"]
    lines = [
        "# Wave B51 ship rule — the pH state, corrected against its thermodynamic sources", "",
        f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
        f"## Verdict: **{payload['verdict']}**", "",
        "| test | baseline (B9, switches off) | (A) B51 refit | pass | (B) B9 vector, switches on | pass |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    a, b = arms["A_refit"], arms["B_frozen_b9"]
    for key, label in (("T2_holdout_sulfur", "T2 hold-out sulfur rows within 3x"),
                       ("T3_hofmann_ph5", "T3 Hofmann ribose + xylose pH 5")):
        lines.append(f"| {label} | {_fmt(a[key]['baseline'])} | {_fmt(a[key]['arm'])} | {a[key]['pass']} | "
                     f"{_fmt(b[key]['arm'])} | {b[key]['pass']} |")
    d0, da, db = a["T4_directional"]["baseline"], a["T4_directional"]["arm"], b["T4_directional"]["arm"]
    lines.append(f"| T4 directional, independent agree | {d0['agree']}/{d0['evaluable']} | {da['agree']}/{da['evaluable']} | "
                 f"{a['T4_directional']['pass']} | {db['agree']}/{db['evaluable']} | {b['T4_directional']['pass']} |")
    p0, pa, pb = a["T5_directional_ph"]["baseline"], a["T5_directional_ph"]["arm"], b["T5_directional_ph"]["arm"]
    lines.append(f"| T5 pH-axis claims (recorded) | {p0['agree']}/{p0['evaluable']} | {pa['agree']}/{pa['evaluable']} | — | "
                 f"{pb['agree']}/{pb['evaluable']} | — |")
    lines.append(f"| union panel within 3x (recorded) | {_fmt(a['union_panel']['baseline'])} | {_fmt(a['union_panel']['arm'])} | — | "
                 f"{_fmt(b['union_panel']['arm'])} | — |")
    c = t1["cost_on_b47_52_rows"]
    lines += ["", f"**T1 (recorded), cost on B47's 52 rows:** baseline {c['baseline']:.2f}, (A) {c['A_refit']:.2f}, "
              f"(B) {c['B_frozen_b9']:.2f}. Laplace for (A): {payload['laplace_A'] or 'not computed'}.", "",
              "**Per-switch attribution** — B9's frozen vector on its own 54 rows, one switch on at a time:", "",
              "| switch | cost |", "| --- | ---: |"]
    lines += [f"| {k} | {v:.3f} |" for k, v in t1["attribution_b9_vector_on_b9_54_rows"].items()]
    lines += ["", "### Rows that moved by more than 0.05 dex", ""]
    if payload["rows_moved_over_0_05_dex"]:
        lines += ["| arm | bundle | compound | fold, baseline | fold, arm | change in abs log10 |",
                  "| --- | --- | --- | ---: | ---: | ---: |"]
        lines += [f"| {m['arm']} | {m['benchmark']} | {m['compound']} | {m['fold_baseline']:.3g} | {m['fold_arm']:.3g} | "
                  f"{m['delta_abs_log10']:+.3f} |" for m in payload["rows_moved_over_0_05_dex"]]
    else:
        lines.append("None.")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
