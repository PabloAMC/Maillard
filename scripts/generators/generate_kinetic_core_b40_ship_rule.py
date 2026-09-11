#!/usr/bin/env python
"""
B40 ship rule (2026-09-11), judged on FROZEN artifacts: the B40 fit report, and the panel
scorecard BEFORE (results/validation/_b40_baseline/core_panel_scores_before_b40.json, the tracked
panel at commit c365533, before the fit) and AFTER (.../core_panel_scores_after_b40.json, the panel
regenerated with the fitted literals installed). Pre-registration sec. 6:

  SHIP if P1 (the fed rows fit), P4 (the Leitzen hold-out improves without breaking) and P6 (no row
  now within 3x leaves the band) hold. P2, P3, P5 are reported.

Run inside docker:  PYTHONPATH=/workspace python scripts/generators/generate_kinetic_core_b40_ship_rule.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402

VAL = data_paths.VALIDATION_DIR
PREREG = VAL / "kinetic_core_b40_prereg.md"
FIT = VAL / "kinetic_core_b40_fit_report.json"
BEFORE = VAL / "_b40_baseline" / "core_panel_scores_before_b40.json"
AFTER = VAL / "_b40_baseline" / "core_panel_scores_after_b40.json"
OUT_JSON = VAL / "kinetic_core_b40_ship_rule.json"
OUT_MD = VAL / "kinetic_core_b40_ship_rule.md"
LEITZEN = "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021"


def _rows(panel: Dict[str, Any]) -> Dict[Tuple[str, str], Dict[str, Any]]:
    out = {}
    for b in panel["benchmarks"]:
        for c in b.get("compounds", []):
            out[(b["benchmark_id"], c["compound"])] = c
    return out


def main() -> int:
    for f in (FIT, BEFORE, AFTER):
        assert f.exists(), f"missing frozen input {f}"
    fit = json.loads(FIT.read_text()); before = json.loads(BEFORE.read_text()); after = json.loads(AFTER.read_text())
    rb, ra = _rows(before), _rows(after)
    # P4 -- the hold-out the fit never read
    def fe(rows, comp):
        r = rows.get((LEITZEN, comp)); return None if r is None else float(r["fold_error"])
    p4 = {"3,4-dideoxyglucosone": (fe(rb, "3,4-dideoxyglucosone"), fe(ra, "3,4-dideoxyglucosone")),
          "3-deoxyglucosone": (fe(rb, "3-deoxyglucosone"), fe(ra, "3-deoxyglucosone")),
          "5-Hydroxymethylfurfural (HMF)": (fe(rb, "5-Hydroxymethylfurfural (HMF)"), fe(ra, "5-Hydroxymethylfurfural (HMF)"))}
    hmf_b, hmf_a = p4["5-Hydroxymethylfurfural (HMF)"]
    p4_pass = (p4["3,4-dideoxyglucosone"][1] is not None and p4["3,4-dideoxyglucosone"][1] < 10.0
               and hmf_a is not None and hmf_a <= max(hmf_b or 1e9, 3.0) * 1.03)
    tdg_ok = p4["3-deoxyglucosone"][1] is not None and p4["3-deoxyglucosone"][1] <= 3.0
    # P6 -- no row now within 3x leaves the band
    within_before = {k for k, r in rb.items() if r.get("within_band")}
    within_after = {k for k, r in ra.items() if r.get("within_band")}
    left = sorted(within_before - within_after); entered = sorted(within_after - within_before)
    p6_pass = not left
    # everything else that moved, reported
    moved = []
    for k in sorted(set(rb) & set(ra)):
        b, a = float(rb[k]["fold_error"]), float(ra[k]["fold_error"])
        if abs(a - b) / max(b, 1e-12) > 0.01:
            moved.append({"benchmark": k[0], "compound": k[1], "fold_before": b, "fold_after": a})
    pr = fit["predictions"]
    p1 = bool(pr["P1_rows_fit_with_tmax"]["held"])
    ships = p1 and p4_pass and p6_pass and tdg_ok
    counts = {"before": {"within": len(within_before), "rows": len(rb)}, "after": {"within": len(within_after), "rows": len(ra)}}
    payload = {
        "artifact": "kinetic_core_b40_ship_rule",
        "provenance": provenance.provenance_block("kinetic_core_b40_ship_rule", generated_by="scripts/generators/generate_kinetic_core_b40_ship_rule.py",
                                                  wave="B40", inputs=[PREREG, FIT, BEFORE, AFTER]),
        "rule": "SHIP if P1 and P4 and P6 and Leitzen 3-DG <= 3x (P3 decisive half); P2, P5 reported (prereg sec. 5)",
        "P1": pr["P1_rows_fit_with_tmax"], "P2": pr["P2_k_tdg_ddg_moves_less_than_b39"], "P3_leitzen_3dg_within_3x": {"fold_before_after": p4["3-deoxyglucosone"], "predicted": "between 2x and 6x", "decisive_rule": "<= 3x", "pass": tdg_ok},
        "P4": {"leitzen_fold_before_after": p4, "rule": "3,4-DGE after < 10x; HMF after within 3 % of before (or <= 3x)", "pass": p4_pass},
        "P5": pr["P5_pinned"],
        "P6": {"left_band": left, "entered_band": entered, "counts": counts, "pass": p6_pass},
        "other_rows_moved": moved,
        "verdict": "SHIP" if ships else "DO NOT SHIP",
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    L = [f"# Wave B40 ship rule: {payload['verdict']}", "", f"_{payload['rule']}. Judged on frozen artifacts: the fit report and the panel before/after._", "",
         f"- **P1** rows fit: **{'HELD' if p1 else 'REFUTED'}** — maxima within 0.3 dex {pr['P1_rows_fit_with_tmax']['maxima_within_0p3_dex']}, t_max in brackets {pr['P1_rows_fit_with_tmax']['tmax_in_brackets']}, χ²_red {pr['P1_rows_fit_with_tmax']['chi2_reduced']:.2f}",
         f"- **P2** k_tdg_ddg moves less than in B39: {'HELD' if payload['P2']['held'] else 'REFUTED'} ({payload['P2']}) — reported",
         f"- **P3** Leitzen 3-DG at or under 3× (decisive): **{'HELD' if tdg_ok else 'REFUTED'}** — {p4['3-deoxyglucosone']}",
         f"- **P4** the Leitzen hold-out, never read by the fit: **{'HELD' if p4_pass else 'REFUTED'}**", "",
         "| Leitzen 2021 row | fold error before | after |", "|---|---:|---:|"]
    for k, (b, a) in p4.items():
        L.append(f"| {k} | {b if b is None else f'{b:.2f}x'} | {a if a is None else f'{a:.2f}x'} |")
    L += ["", f"- **P5** ≥ 3 pinned: {'HELD' if payload['P5']['held'] else 'REFUTED'} ({payload['P5']['n_pinned']} pinned; collinear {payload['P5']['collinear_pairs']}) — reported",
          f"- **P6** no row within 3× leaves the band: **{'HELD' if p6_pass else 'REFUTED'}** — left {left}; entered {entered}; within {counts['before']['within']}/{counts['before']['rows']} → {counts['after']['within']}/{counts['after']['rows']}", ""]
    if moved:
        L += ["## Every other row that moved by more than 1 %", "", "| benchmark | compound | before | after |", "|---|---|---:|---:|"]
        L += [f"| {m['benchmark']} | {m['compound']} | {m['fold_before']:.2f}x | {m['fold_after']:.2f}x |" for m in moved]
    OUT_MD.write_text("\n".join(L) + "\n")
    print(payload["verdict"], "| P1", p1, "| 3-DG<=3x", tdg_ok, "| P4", p4_pass, {k: (None if b is None else round(b, 2), None if a is None else round(a, 2)) for k, (b, a) in p4.items()},
          "| P6", p6_pass, counts, "| moved", len(moved))
    return 0


if __name__ == "__main__":
    sys.exit(main())
