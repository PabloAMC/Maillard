#!/usr/bin/env python
"""
Wave B50 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-10-27).

`results/validation/kinetic_core_b50_prereg.md` sec. 3. The union and directional panels are scored in
child processes with thiamine's H2S branch OFF (share forced to zero) and ON. Writes
`results/validation/kinetic_core_b50_ship_rule.{json,md}`.
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
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
OUT = V / "kinetic_core_b50_ship_rule.json"


def _switch(on: bool) -> None:
    from src.kinetic_core import sulfur

    if not on:
        sulfur.thiamine_h2s_share = lambda ph: 0.0


def _child(on: bool) -> Dict[str, Any]:
    _switch(on)
    from src.kinetic_core.directional import score_panel as score_directional
    from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram, predict
    from src.kinetic_core.ph_state import BufferSpec
    from src.kinetic_core.scoring import score_panel

    rows = [{"benchmark": b["benchmark_id"], "lane": r.get("lane"), "compound": r.get("compound"),
             "fold": r.get("fold_error"), "within": bool(r.get("within_band"))}
            for b in score_panel()["benchmarks"] for r in b["compounds"]]
    ind = score_directional()["summary"]["independent"]["total"]
    buf = BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="B50 probe")
    probe = {}
    for ph in (5.0, 6.5, 7.0):
        spec = FormulationSpec("ribose_thiamine", {"D-Ribose": 50.0, "Thiamine": 10.0},
                               ProcessSpec(ThermalProgram.isothermal(150.0, 10.0), ph=ph, buffer=buf))
        run = predict(spec, ["FFT", "MFT"])
        probe[str(ph)] = dict(run.concentrations_ug_per_l) if run.answered else None
    return {"rows": rows, "directional": {"agree": ind["agree"], "evaluable": ind["evaluable"]}, "probe": probe}


def _score(on: bool) -> Dict[str, Any]:
    done = subprocess.run([sys.executable, __file__, "--child", "on" if on else "off"], cwd=ROOT,
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(done.stderr[-3000:])
    return json.loads(done.stdout.strip().splitlines()[-1])


def _summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    folds = [r["fold"] for r in rows if r["fold"] and math.isfinite(r["fold"]) and r["fold"] > 0]
    return {"rows": len(rows), "within_3x": sum(r["within"] for r in rows),
            "median_abs_log10": statistics.median(abs(math.log10(f)) for f in folds) if folds else None}


def _t1() -> Dict[str, Any]:
    """B9's cost on its 54 rows with the branch ON and OFF, in the same code (see prereg sec. 4 amendment)."""
    import numpy as np

    import generate_kinetic_core_b9_fit as B9
    from generate_kinetic_core_b8_laplace import frozen_vector
    from src.kinetic_core import sulfur

    x = frozen_vector(json.loads((V / "kinetic_core_b9_fit_report.json").read_text(encoding="utf-8")))
    r = B9.residual_vector(x, True)
    on = 0.5 * float(np.dot(r, r))
    saved = sulfur.thiamine_h2s_share
    try:
        sulfur.thiamine_h2s_share = lambda ph: 0.0
        r = B9.residual_vector(x, True)
        off = 0.5 * float(np.dot(r, r))
    finally:
        sulfur.thiamine_h2s_share = saved
    return {"b9_cost": on, "reference": off, "pass": abs(on - off) <= 1e-9 * off}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--child", choices=("on", "off"))
    args = parser.parse_args(argv)
    if args.child:
        print(json.dumps(_child(args.child == "on")))
        return 0
    off, on = _score(False), _score(True)
    t1 = _t1()
    t2 = {"off": _summary(off["rows"]), "on": _summary(on["rows"]),
          "sulfur_off": _summary([r for r in off["rows"] if r["lane"] == "sulfur"]),
          "sulfur_on": _summary([r for r in on["rows"] if r["lane"] == "sulfur"])}
    t2["pass"] = (t2["on"]["within_3x"] >= t2["off"]["within_3x"]
                  and t2["sulfur_on"]["median_abs_log10"] <= t2["sulfur_off"]["median_abs_log10"] + 0.05)
    t3 = {"off": off["directional"], "on": on["directional"], "pass": on["directional"]["agree"] >= off["directional"]["agree"]}
    verdict = "SHIP" if (t1["pass"] and t2["pass"] and t3["pass"]) else "DO NOT SHIP"
    payload = {"artifact": "kinetic_core_b50_ship_rule",
               "provenance": provenance.provenance_block("kinetic_core_b50_ship_rule",
                                                         generated_by="scripts/generators/generate_kinetic_core_b50_ship_rule.py",
                                                         inputs=[ROOT / "src" / "kinetic_core" / "parameters_sulfur.py"]),
               "prereg": data_paths.rel(V / "kinetic_core_b50_prereg.md"),
               "T1_b9_reproduced": t1, "T2_panel": t2, "T3_directional": t3,
               "T4_ribose_thiamine_probe": {"off": off["probe"], "on": on["probe"]}, "verdict": verdict}
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"{verdict}: T1 {t1['pass']} ({t1['b9_cost']!r}); T2 union {t2['off']['within_3x']}/{t2['off']['rows']} -> "
          f"{t2['on']['within_3x']}/{t2['on']['rows']}, sulfur median {t2['sulfur_off']['median_abs_log10']:.3f} -> "
          f"{t2['sulfur_on']['median_abs_log10']:.3f}; T3 {t3['off']['agree']} -> {t3['on']['agree']}; probe FFT on {on['probe']}")
    return 0


def render(payload: Dict[str, Any]) -> str:
    t1, t2, t3, t4 = (payload[k] for k in ("T1_b9_reproduced", "T2_panel", "T3_directional", "T4_ribose_thiamine_probe"))
    lines = ["# Wave B50 ship rule — thiamine releases H2S", "", f"*Pre-registration: `{payload['prereg']}`.*", "",
             f"## Verdict: **{payload['verdict']}**", "", "| test | branch off | branch on | pass |", "| --- | --- | --- | --- |",
             f"| T1 B9 cost, branch off vs on (same code) | {t1['reference']:.12g} | {t1['b9_cost']:.12g} | {t1['pass']} |",
             f"| T2 union within 3x | {t2['off']['within_3x']}/{t2['off']['rows']} | {t2['on']['within_3x']}/{t2['on']['rows']} | {t2['pass']} |",
             f"| T2 sulfur median abs log10 | {t2['sulfur_off']['median_abs_log10']:.3f} | {t2['sulfur_on']['median_abs_log10']:.3f} | |",
             f"| T3 directional, independent agree | {t3['off']['agree']}/{t3['off']['evaluable']} | {t3['on']['agree']}/{t3['on']['evaluable']} | {t3['pass']} |",
             "", "### T4 ribose 50 mM + thiamine 10 mM, 150 C, 10 min, 0.5 M phosphate (ug/L)", "",
             "| pH | FFT off | FFT on | MFT off | MFT on |", "| --- | ---: | ---: | ---: | ---: |"]
    for ph in ("5.0", "6.5", "7.0"):
        o, n = t4["off"].get(ph) or {}, t4["on"].get(ph) or {}
        lines.append(f"| {ph} | {o.get('FFT', float('nan')):.3g} | {n.get('FFT', float('nan')):.3g} | {o.get('MFT', float('nan')):.3g} | {n.get('MFT', float('nan')):.3g} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
