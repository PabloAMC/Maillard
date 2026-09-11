#!/usr/bin/env python
"""
Wave B22b -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-10).

`results/validation/kinetic_core_b22b_prereg.md` section 4.

  T1  the shape, DECISIVE: all five rows within 0.4 dex AND the model peaks between 60 and 180 min
  T2  the two-arm direction check -- reported as VACUOUS, see below
  T3  identification
  T4  nothing else moves
Ship rule: SHIP if T1, T2, T3 and T4 hold.

T2 CANNOT BE EVALUATED AS WRITTEN, and that is this rule's own finding. The pre-registration said
the binary Met + Glc arm "rests on charging methionine as glycine for the Amadori chemistry". It
does not rest on anything: wave B22 did not ship, so its four steps are inert at zero, and there is
no route from methionine and glucose to methional in the model at all. The binary arm is
STRUCTURALLY zero, so the ratio is a division by nothing rather than a wrong number. Recorded as
vacuous instead of failed, because those are different and only one of them is about this wave.
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
REPORT = V / "kinetic_core_b22b_fit_report.json"
OUT = V / "kinetic_core_b22b_ship_rule.json"


def main() -> int:
    rep = json.loads(REPORT.read_text(encoding="utf-8"))
    resid = rep["residual_by_row_dex"]
    peak = float(rep["peak_time_min"])
    worst = max((abs(v), k) for k, v in resid.items())
    t1 = {"residuals_dex": resid, "worst_row_min": worst[1], "worst_dex": worst[0],
          "model_peak_min": peak, "printed_peak_min": 120.0,
          "peak_in_window": bool(60.0 <= peak <= 180.0),
          "pass": bool(worst[0] <= 0.4 and 60.0 <= peak <= 180.0)}
    arm = rep["two_arm_check"]
    model_binary = [r["model_binary_umol_l"] for r in arm["rows"].values()]
    t2 = {"state": "VACUOUS",
          "why": ("wave B22 did not ship, so its four steps are inert at zero and the model has NO "
                  "route from methionine and glucose to methional. The binary arm is structurally "
                  "zero, so the ratio divides by nothing. The pre-registration assumed the arm "
                  "rested on a declared substitution; it rests on nothing."),
          "model_binary_umol_l": model_binary,
          "printed_arp_over_binary": [r["printed_arp_over_binary"] for r in arm["rows"].values()],
          "pass": None}
    lap = rep["laplace"]
    t3 = {"sigma": lap["sigma"], "identified": lap["identified"], "on_bound": lap["on_bound"],
          "pass": bool(all(lap["identified"].values()) and not any(lap["on_bound"].values()))}
    try:
        tracked = json.loads(subprocess.check_output(
            ["git", "show", "HEAD:" + data_paths.rel(V / "core_panel_scores.json")], cwd=ROOT, text=True))
        from src.kinetic_core import scoring
        live = scoring.score_panel()
        # 2026-09-11 (review of PR #16): compared as strings, byte for byte, which provenance.py
        # documents as wrong for a scorecard -- the same panel differs by ~7e-8 relative between
        # arm64 and x86 runners. The tolerance-aware comparison the freshness gate uses.
        same = not provenance.payload_differences(
            [b.get("compounds") for b in tracked["benchmarks"]],
            [b.get("compounds") for b in live["benchmarks"]],
        )
    except Exception:
        same = False
    t4 = {"panel_identical": bool(same), "pass": bool(same)}

    ships = bool(t1["pass"] and t3["pass"] and t4["pass"] and t2["pass"])
    payload = {
        "artifact": "kinetic_core_b22b_ship_rule",
        "provenance": provenance.provenance_block(
            "kinetic_core_b22b_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b22b_ship_rule.py", inputs=[REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b22b_prereg.md"),
        "rule": "SHIP if T1 (the shape and its peak), T2, T3 and T4 hold",
        "T1": t1, "T2": t2, "T3": t3, "T4": t4,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=_render)
    print(payload["verdict"], "| T1", t1["pass"], "peak", peak, "| T2", t2["state"],
          "| T3", t3["pass"], "| T4", t4["pass"])
    return 0


def _render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4 = p["T1"], p["T2"], p["T3"], p["T4"]
    L = [f"# Wave B22b ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}.*", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 the shape | worst {T1['worst_dex']:+.2f} dex at {T1['worst_row_min']} min; model peaks "
         f"at {T1['model_peak_min']:.0f} min against a printed 120 | {T1['pass']} |",
         f"| T2 the two-arm check | **{T2['state']}** | -- |",
         f"| T3 identification | sigma {json.dumps({k: round(v, 3) for k, v in T3['sigma'].items()})} | {T3['pass']} |",
         f"| T4 nothing else moves | panel identical {T4['panel_identical']} | {T4['pass']} |", "",
         f"> **T2 is vacuous, not failed.** {T2['why']}", ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
