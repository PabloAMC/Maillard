#!/usr/bin/env python
"""
Wave B20 -- THE PRE-REGISTERED SHIP RULE, evaluated (2026-09-09).

`results/validation/kinetic_core_b20_prereg.md` sec. 4, computed from the frozen B20 report and the
live engine, written to `results/validation/kinetic_core_b20_ship_rule.{json,md}`.

  T1  the decisive rows (the prereg's list: k3 at 120 and 130 C, k7 at 130, k8 at 120 and 130, k9 at 130,
      k11 at 130) within 0.3 dex of the printed value; the rest reported
  T2  the panel is untouched: no currently scored row moves by more than 0.05 dex
  T3  Nguyen's pot at 120 C / 30 min: CML within the printed 0.025-0.135 mmol/L and CEL below CML
  T4  another laboratory: Berk 2021 (180 C, dry sesame) and Hamzalioglu 2026 (milk, 110-140 C), in decades; reported
  T5  identification: Laplace sigma below one decade on every coordinate, none on its bound
Ship rule: SHIP if T1, T2 and T5 hold; T3 and T4 are reported.
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
REPORT = V / "kinetic_core_b20_fit_report.json"
SCORECARD = V / "core_panel_scores.json"
OUT = V / "kinetic_core_b20_ship_rule.json"
#: The prereg's decisive rows (sec. 4, T1), by row id.
DECISIVE = ("nguyen_k_glyc_120C", "nguyen_k_glyc_130C", "nguyen_k_flp_cml_130C", "nguyen_k_flp_decay_120C",
            "nguyen_k_flp_decay_130C", "nguyen_k_flp_cel_130C", "nguyen_k_cml_loss_130C")


def _read(p: Path) -> Dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def t1(report):
    rows = {r["id"]: r for r in report["rows"]}
    decisive = {rid: rows[rid]["residual_dex"] for rid in DECISIVE}
    reported = {rid: r["residual_dex"] for rid, r in rows.items() if rid not in DECISIVE}
    worst = max(decisive.items(), key=lambda kv: abs(kv[1]))
    return {"decisive_rows_dex": decisive, "reported_rows_dex": reported, "worst_row": worst[0], "worst_dex": worst[1],
            "pass": bool(all(abs(v) <= 0.3 for v in decisive.values()))}


def _predicted_leaves(o, path=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(_predicted_leaves(v, f"{path}/{k}"))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(_predicted_leaves(v, f"{path}[{i}]"))
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        if "pred" in path.lower():
            out[path] = float(o)
    return out


def t2():
    """The tracked scorecard (git HEAD) against the live one: every predicted number's log10 change."""
    from src.kinetic_core import scoring

    try:
        tracked = json.loads(subprocess.check_output(["git", "show", "HEAD:" + data_paths.rel(SCORECARD)], cwd=ROOT, text=True))
    except Exception as exc:  # pragma: no cover
        return {"status": f"tracked scorecard unavailable: {exc}", "pass": False}
    live = scoring.score_panel()
    old = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in tracked["benchmarks"]}
    new = {b["benchmark_id"]: _predicted_leaves(b.get("compounds")) for b in live["benchmarks"]}
    deltas = []
    for bid, leaves in old.items():
        for path, v in leaves.items():
            w = new.get(bid, {}).get(path)
            if w is None or v <= 0 or w <= 0:
                continue
            deltas.append((abs(math.log10(w / v)), bid, path))
    deltas.sort(reverse=True)
    worst = deltas[0] if deltas else (0.0, None, None)
    refused_old = sum(len(b.get("refused_compounds", [])) for b in tracked["benchmarks"])
    refused_new = sum(len(b.get("refused_compounds", [])) for b in live["benchmarks"])
    return {"n_compared": len(deltas), "worst_dex": worst[0], "worst_benchmark": worst[1], "worst_path": worst[2],
            "n_rows_moved_over_0.01_dex": sum(1 for d in deltas if d[0] > 0.01), "refused_rows_before_after": [refused_old, refused_new],
            "pass": bool(worst[0] < 0.05)}


def t3(report):
    pot = report["diagnostics"]["nguyen_pot_30min"]["120C"]
    lo, hi = report["diagnostics"]["nguyen_printed_cml_range_mmol_l"]
    return {"cml_mmol_l": pot["CML"], "cel_mmol_l": pot["CEL"], "printed_range": [lo, hi], "lysine_lost_fraction": pot["lysine_lost_fraction"],
            "pass": bool(lo <= pot["CML"] <= hi and pot["CEL"] < pot["CML"])}


def t4(report):
    d = report["diagnostics"]
    return {"berk2021_180C_dex": d["berk2021_flp_to_cml_180C"]["dex"],
            "hamzalioglu2026_dex": {t: v["dex"] for t, v in d["hamzalioglu2026_laclys_to_cml"].items()},
            "troise2015": d["troise2015_direction"], "reported_only": True}


def t5(report):
    lap = report["laplace"]
    return {"sigma_dex": lap["sigma"], "on_bound": lap["on_bound"],
            "pass": bool(all(s < 1.0 for s in lap["sigma"].values()) and not any(lap["on_bound"].values()))}


def render(p: Dict[str, Any]) -> str:
    T1, T2, T3, T4, T5 = p["T1"], p["T2"], p["T3"], p["T4"], p["T5"]
    L = [f"# Wave B20 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
         f"Cost {p['cost']:.2f} on {p['n_rows']} rows; fitted constants (log10 at 100 C): {json.dumps({k: round(v, 3) for k, v in p['frozen'].items()})}", "",
         "| test | result | pass |", "|---|---|---|",
         f"| T1 decisive rows | worst {T1['worst_row']} {T1['worst_dex']:+.2f} dex; all decisive {json.dumps({k: round(v, 2) for k, v in T1['decisive_rows_dex'].items()})}; reported {json.dumps({k: round(v, 2) for k, v in T1['reported_rows_dex'].items()})} | {T1['pass']} |",
         f"| T2 panel untouched | {T2.get('n_compared')} predicted numbers compared; worst change {T2.get('worst_dex', float('nan')):.2e} dex; refused rows {T2.get('refused_rows_before_after')} | {T2['pass']} |",
         f"| T3 Nguyen's pot 120 C / 30 min | CML {T3['cml_mmol_l']:.4f} mmol/L (printed {T3['printed_range']}), CEL {T3['cel_mmol_l']:.4f}, lysine lost {100 * T3['lysine_lost_fraction']:.1f} % | {T3['pass']} |",
         f"| T4 other laboratories | Berk 2021 180 C (dry sesame) {T4['berk2021_180C_dex']:+.2f} dex; Hamzalioglu 2026 milk {json.dumps({k: round(v, 2) for k, v in T4['hamzalioglu2026_dex'].items()})} dex; {T4['troise2015']} | reported |",
         f"| T5 identification | sigma {json.dumps({k: round(v, 3) for k, v in T5['sigma_dex'].items()})}; on bound {any(T5['on_bound'].values())} | {T5['pass']} |",
         ""]
    return "\n".join(L)


def main() -> int:
    if not REPORT.exists():
        raise SystemExit(f"{REPORT} missing")
    rep = _read(REPORT)
    T1, T2, T3, T4, T5 = t1(rep), t2(), t3(rep), t4(rep), t5(rep)
    ships = bool(T1["pass"] and T2["pass"] and T5["pass"])
    payload = {
        "artifact": "kinetic_core_b20_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b20_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b20_ship_rule.py",
                                                  inputs=[REPORT]),
        "prereg": data_paths.rel(V / "kinetic_core_b20_prereg.md"),
        "rule": "SHIP if T1 (the prereg's decisive rows within 0.3 dex), T2 (no scored panel row moves 0.05 dex) and T5 (every coordinate identified, off its bound) hold; T3 and T4 reported",
        "cost": rep["objective"]["final_cost"], "n_rows": rep["objective"]["n_rows"], "frozen": rep["frozen_parameters"]["glycation"],
        "T1": T1, "T2": T2, "T3": T3, "T4": T4, "T5": T5,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    artifact_io.write_artifact(payload, OUT, render=render)
    print(payload["verdict"], "| T1", T1["pass"], T1["worst_dex"], "| T2", T2["pass"], T2.get("worst_dex"), "| T3", T3["pass"], "| T5", T5["pass"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
