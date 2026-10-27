#!/usr/bin/env python
"""
Wave B51.1 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b51_1_prereg.md` sec. 3. B9's frozen vector scored twice, each in a child
process: all B51 switches off (baseline, what ships) and `label_scale_tables` + `sid_across_segments` on
(candidate). T2-T4 are B51's tests (`generate_kinetic_core_b51_ship_rule._tests`, imported, not copied).
Writes `results/validation/kinetic_core_b51_1_ship_rule.{json,md}`.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

import generate_kinetic_core_b51_ship_rule as B51R  # noqa: E402
from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_b51_1_ship_rule.json"
CANDIDATE = {"label_scale_tables": True, "sid_across_segments": True}


def _config(mode: str):
    from src.kinetic_core import ph_state

    return ph_state.PhStateCorrections(**CANDIDATE) if mode == "on" else ph_state.PhStateCorrections()


def _child(mode: str) -> Dict[str, Any]:
    from src.kinetic_core import engine, ph_state

    ph_state.B51 = _config(mode)
    engine._B2_FIT_REPORT = B51R.B9_REPORT
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
    return {"mode": mode, "rows": rows,
            "directional_independent": {"agree": summary["total"]["agree"], "evaluable": summary["total"]["evaluable"]},
            "directional_ph": {"agree": ph.get("agree"), "evaluable": ph.get("evaluable"), "misses": ph.get("misses")}}


def _cost(mode: str) -> float:
    import numpy as np

    import generate_kinetic_core_b9_fit as B9
    from generate_kinetic_core_b8_laplace import frozen_vector
    from src.kinetic_core import ph_state

    ph_state.B51 = _config(mode)
    r = B9.residual_vector(frozen_vector(B51R._read(B51R.B9_REPORT)), True)
    return 0.5 * float(np.dot(r, r))


def _run(*args: str) -> Any:
    done = subprocess.run([sys.executable, __file__, *args], cwd=ROOT, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"child {args} failed:\n{done.stderr[-3000:]}")
    return json.loads(done.stdout.strip().splitlines()[-1])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", choices=("on", "off"))
    parser.add_argument("--cost", choices=("on", "off"))
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child)))
        return 0
    if args.cost:
        print(json.dumps(_cost(args.cost)))
        return 0
    base, cand = _run("--child", "off"), _run("--child", "on")
    tests = B51R._tests(base, cand)
    t1 = {"baseline": _run("--cost", "off"), "candidate": _run("--cost", "on")}
    verdict = "SHIP" if tests["pass"] else "DO NOT SHIP"
    by = {(r["benchmark"], r["compound"]): r for r in base["rows"]}
    moved = []
    for r in cand["rows"]:
        b = by.get((r["benchmark"], r["compound"]))
        if b and b["fold"] and r["fold"] and b["fold"] > 0 and r["fold"] > 0:
            import math

            d = abs(math.log10(r["fold"])) - abs(math.log10(b["fold"]))
            if abs(d) > 0.01:
                moved.append({"benchmark": r["benchmark"], "compound": r["compound"], "fold_baseline": b["fold"],
                              "fold_candidate": r["fold"], "delta_abs_log10": d})
    payload = {
        "artifact": "kinetic_core_b51_1_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b51_1_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b51_1_ship_rule.py",
            inputs=[ROOT / "src" / "kinetic_core" / "sulfur.py", B51R.B9_REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b51_1_prereg.md"),
        "candidate_switches": CANDIDATE, "T1_cost_b9_rows": t1, "tests": tests, "rows_moved_over_0_01_dex": moved,
        "verdict": verdict, "rule": "SHIP if T2, T3 and T4 hold",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict} | T2 {tests['T2_holdout_sulfur']['pass']} T3 {tests['T3_hofmann_ph5']['pass']} "
          f"T4 {tests['T4_directional']['pass']} | cost {t1['baseline']:.3f} -> {t1['candidate']:.3f} | moved {len(moved)}")
    return 0


def render(payload: Dict[str, Any]) -> str:
    t, c = payload["tests"], payload["T1_cost_b9_rows"]
    f = B51R._fmt
    lines = [
        "# Wave B51.1 ship rule — the 25 °C scale for the measured pH tables, and the SID across segments", "",
        f"*Pre-registration: `{payload['prereg']}`. Rule: {payload['rule']}.*", "",
        f"## Verdict: **{payload['verdict']}**", "",
        "| test | baseline (B9, all off) | candidate (two switches on) | pass |", "| --- | --- | --- | --- |",
        f"| T1 B9 cost on its 54 rows (recorded) | {c['baseline']:.3f} | {c['candidate']:.3f} | — |",
        f"| T2 hold-out sulfur rows within 3x | {f(t['T2_holdout_sulfur']['baseline'])} | {f(t['T2_holdout_sulfur']['arm'])} | {t['T2_holdout_sulfur']['pass']} |",
        f"| T3 Hofmann ribose + xylose pH 5 | {f(t['T3_hofmann_ph5']['baseline'])} | {f(t['T3_hofmann_ph5']['arm'])} | {t['T3_hofmann_ph5']['pass']} |",
        f"| T4 directional, independent agree | {t['T4_directional']['baseline']['agree']}/{t['T4_directional']['baseline']['evaluable']} | "
        f"{t['T4_directional']['arm']['agree']}/{t['T4_directional']['arm']['evaluable']} | {t['T4_directional']['pass']} |",
        f"| pH-axis claims (recorded) | {t['T5_directional_ph']['baseline']['agree']}/{t['T5_directional_ph']['baseline']['evaluable']} | "
        f"{t['T5_directional_ph']['arm']['agree']}/{t['T5_directional_ph']['arm']['evaluable']} | — |",
        f"| union panel within 3x (recorded) | {f(t['union_panel']['baseline'])} | {f(t['union_panel']['arm'])} | — |",
        "", "### Rows that moved by more than 0.01 dex", "",
    ]
    if payload["rows_moved_over_0_01_dex"]:
        lines += ["| bundle | compound | fold, baseline | fold, candidate | change in abs log10 |", "| --- | --- | ---: | ---: | ---: |"]
        lines += [f"| {m['benchmark']} | {m['compound']} | {m['fold_baseline']:.3g} | {m['fold_candidate']:.3g} | {m['delta_abs_log10']:+.3f} |"
                  for m in payload["rows_moved_over_0_01_dex"]]
    else:
        lines.append("None.")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
