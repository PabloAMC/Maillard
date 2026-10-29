#!/usr/bin/env python
"""
B55 / B56 -- THE PRE-REGISTERED SHIP RULES (2026-10-28).

``results/validation/kinetic_core_b55_b56_prereg.md`` sec. 5, written before any coverage number under
either change existed. Four envelope arms of the SAME code (n = 200, seed 0, the union panel), written by
``generate_core_prediction_uncertainty.py --b55 {on,off} --b56 {on,off}`` into
``results/validation/_b55_b56_arms/``:

  A baseline (off, off)   B b55 (on, off)   C b56 (off, on)   D b55_b56 (on, on)

B55 (B vs A) SHIPS iff
  T1  every row of an isothermal 145 C pot, and every row off the sulfur lane, has p5 / p50 / p95 equal
      to A within a relative 1e-9;
  T2  over the sulfur rows away from 145 C the median change in 90 % width is >= +0.10 dex and no such
      row narrows by more than 0.05 dex;
  T3  out-of-sample coverage (fixed denominator) in B >= A - 1 row.
B56 (D vs B if B55 ships, else C vs A) SHIPS iff
  T4  out-of-sample 90 % coverage (fixed denominator) lies in [70 %, 95 %]; the median out-of-sample
      width and the per-lane coverage and width are reported, with no threshold.

FIXED DENOMINATOR: every out-of-sample literature row with an envelope formed in either compared arm; a
not-evaluable (zero-width) interval, or a row missing from an arm, counts as a miss.

Writes ``results/validation/kinetic_core_b55_b56_ship_rule.{json,md}``.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
ARMS_DIR = V / "_b55_b56_arms"
ARMS = {
    "A": ARMS_DIR / "core_prediction_uncertainty_baseline.json",
    "B": ARMS_DIR / "core_prediction_uncertainty_b55.json",
    "C": ARMS_DIR / "core_prediction_uncertainty_b56.json",
    "D": ARMS_DIR / "core_prediction_uncertainty_b55_b56.json",
}
OUT = V / "kinetic_core_b55_b56_ship_rule.json"
PREREG = "results/validation/kinetic_core_b55_b56_prereg.md"
SULFUR_T_REF_C = 145.0
REL_TOL = 1e-9
T2_MEDIAN_WIDENING_DEX = 0.10
T2_MAX_NARROWING_DEX = 0.05
T4_BAND = (0.70, 0.95)
MIN_EVALUABLE = 0.01   # uncertainty.MIN_EVALUABLE_CI_WIDTH_LOG10


Key = Tuple[str, str]


def _rows(path: Path) -> Dict[Key, Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    out: Dict[Key, Dict[str, Any]] = {}
    for b in payload["benchmarks"]:
        for c in b["compounds"]:
            out[(b["benchmark_id"], c["compound"])] = {
                **c,
                "bench_file": b["bench_file"],
                "signal_origin": b["signal_origin"],
                "evidence_role": b["evidence_role"],
            }
    return out


def _isothermal_at(bench_file: str) -> Optional[float]:
    """The pot's temperature if it is held at one temperature, else None."""
    from src.kinetic_core.panel import core_spec, load_bundle

    spec = core_spec(load_bundle(ROOT / bench_file), use_buffer=True)
    temps = {float(t) for _minutes, t in spec.process.thermal.segments}
    return temps.pop() if len(temps) == 1 else None


def _out_of_sample(row: Dict[str, Any]) -> bool:
    return (row["signal_origin"] == "external_literature" and row["evidence_role"] != "fit_recovery"
            and not row.get("in_core_fit"))


def _hit(row: Optional[Dict[str, Any]]) -> bool:
    if row is None:
        return False
    w = row.get("ci_width_log10")
    return bool(row["inside_ci"]) and w is not None and w > MIN_EVALUABLE


def _coverage(arm: Dict[Key, Dict[str, Any]], keys: List[Key], rows_meta: Dict[Key, Dict[str, Any]]) -> Dict[str, Any]:
    hits = sum(1 for k in keys if _hit(arm.get(k)))
    widths = [arm[k]["ci_width_log10"] for k in keys if k in arm and arm[k].get("ci_width_log10") is not None]
    by_lane: Dict[str, Dict[str, Any]] = {}
    for lane in sorted({rows_meta[k]["lane"] for k in keys}):
        lk = [k for k in keys if rows_meta[k]["lane"] == lane]
        lw = [arm[k]["ci_width_log10"] for k in lk if k in arm and arm[k].get("ci_width_log10") is not None]
        by_lane[lane] = {
            "hits": sum(1 for k in lk if _hit(arm.get(k))), "total": len(lk),
            "median_width_log10": statistics.median(lw) if lw else None,
        }
    return {
        "hits": hits, "total": len(keys), "rate": hits / len(keys) if keys else None,
        "median_width_log10": statistics.median(widths) if widths else None,
        "by_lane": by_lane,
    }


def _same(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    for f in ("predicted_p5", "predicted_p50", "predicted_p95"):
        x, y = float(a[f]), float(b[f])
        if abs(x - y) > REL_TOL * max(abs(x), abs(y), 1e-300):
            return False
    return True


def build() -> Dict[str, Any]:
    arms = {name: _rows(path) for name, path in ARMS.items()}
    universe = sorted(set().union(*[set(a) for a in arms.values()]))
    meta: Dict[Key, Dict[str, Any]] = {}
    for k in universe:
        row = next(a[k] for a in arms.values() if k in a)
        meta[k] = {"lane": row["lane"], "oos": _out_of_sample(row), "bench_file": row["bench_file"],
                   "T_C": _isothermal_at(row["bench_file"])}
    oos = [k for k in universe if meta[k]["oos"]]
    A, B, C, D = arms["A"], arms["B"], arms["C"], arms["D"]

    # ---- B55 ---------------------------------------------------------------
    t1_rows = [k for k in universe if meta[k]["lane"] != "sulfur" or meta[k]["T_C"] == SULFUR_T_REF_C]
    t1_fail = [f"{k[0]} / {k[1]}" for k in t1_rows if (k in A) != (k in B) or (k in A and not _same(A[k], B[k]))]
    t2_rows = [k for k in universe if meta[k]["lane"] == "sulfur" and meta[k]["T_C"] != SULFUR_T_REF_C
               and k in A and k in B]
    t2_table = []
    for k in t2_rows:
        wa, wb = A[k]["ci_width_log10"], B[k]["ci_width_log10"]
        t2_table.append({
            "benchmark_id": k[0], "compound": k[1], "T_C": meta[k]["T_C"],
            "width_A": wa, "width_B": wb, "delta_width": (wb - wa) if (wa is not None and wb is not None) else None,
            "p50_A": A[k]["predicted_p50"], "p50_B": B[k]["predicted_p50"], "measured": A[k]["measured"],
            "inside_A": _hit(A[k]), "inside_B": _hit(B[k]),
        })
    deltas = [r["delta_width"] for r in t2_table if r["delta_width"] is not None]
    t2_median = statistics.median(deltas) if deltas else None
    t2_min = min(deltas) if deltas else None
    cov = {name: _coverage(arms[name], oos, meta) for name in arms}
    t1 = not t1_fail
    t2 = bool(deltas) and t2_median >= T2_MEDIAN_WIDENING_DEX and t2_min >= -T2_MAX_NARROWING_DEX
    t3 = cov["B"]["hits"] >= cov["A"]["hits"] - 1
    b55_ship = t1 and t2 and t3

    # ---- B56 ---------------------------------------------------------------
    b56_arm, b56_ref = ("D", "B") if b55_ship else ("C", "A")
    rate = cov[b56_arm]["rate"]
    t4 = rate is not None and T4_BAND[0] <= rate <= T4_BAND[1]
    b56_ship = t4

    sigma = json.loads((V / "kinetic_core_b56_structural_discrepancy.json").read_text(encoding="utf-8"))
    own = {name: json.loads(path.read_text(encoding="utf-8"))["summary"] for name, path in ARMS.items()}
    return {
        "artifact": "kinetic_core_b55_b56_ship_rule",
        "prereg": PREREG,
        "provenance": provenance.provenance_block(
            "kinetic_core_b55_b56_ship_rule",
            generated_by="scripts/generators/generate_kinetic_core_b55_b56_ship_rule.py",
            inputs=(*ARMS.values(), V / "kinetic_core_b56_structural_discrepancy.json"),
        ),
        "arms": {name: data_paths.rel(path) for name, path in ARMS.items()},
        "fixed_denominator_out_of_sample_rows": len(oos),
        "coverage_out_of_sample": cov,
        "artifact_own_definitions": {
            name: {
                "honest_literature": [s["honest_literature_coverage"][x] for x in ("hits", "total", "not_evaluable")],
                "honest_literature_median_width": s["honest_literature_coverage"]["median_ci_width_log10"],
                "out_of_sample": [s["out_of_sample_literature_coverage"][x] for x in ("hits", "total", "not_evaluable")],
            }
            for name, s in own.items()
        },
        "b55": {
            "T1_invariance_at_145C_and_off_sulfur": {"pass": t1, "rows_checked": len(t1_rows), "failures": t1_fail},
            "T2_widening_away_from_145C": {
                "pass": t2, "median_delta_width": t2_median, "min_delta_width": t2_min,
                "thresholds": {"median_at_least": T2_MEDIAN_WIDENING_DEX, "no_row_narrows_by_more_than": T2_MAX_NARROWING_DEX},
                "rows": t2_table,
            },
            "T3_no_coverage_loss": {"pass": t3, "hits_A": cov["A"]["hits"], "hits_B": cov["B"]["hits"],
                                     "total": len(oos)},
            "verdict": "SHIP" if b55_ship else "DO NOT SHIP",
        },
        "b56": {
            "sigma_log10_by_lane": sigma["sigma_log10_by_lane"],
            "judged_arm": b56_arm, "reference_arm": b56_ref,
            "T4_out_of_sample_coverage_in_band": {
                "pass": t4, "hits": cov[b56_arm]["hits"], "total": cov[b56_arm]["total"], "rate": rate,
                "band": list(T4_BAND), "median_width_log10": cov[b56_arm]["median_width_log10"],
                "reference_rate": cov[b56_ref]["rate"], "reference_median_width_log10": cov[b56_ref]["median_width_log10"],
            },
            "verdict": "SHIP" if b56_ship else "DO NOT SHIP",
        },
    }


def _f(x: Any, nd: int = 3) -> str:
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def render(p: Dict[str, Any]) -> str:
    b55, b56, cov = p["b55"], p["b56"], p["coverage_out_of_sample"]
    out = [
        "# B55 / B56 ship rules",
        "",
        f"Pre-registration: `{p['prereg']}`. Four arms of the same code (n = 200, seed 0). Out-of-sample "
        f"coverage on a FIXED denominator of {p['fixed_denominator_out_of_sample_rows']} rows "
        "(a not-evaluable interval counts as a miss).",
        "",
        f"**B55: {b55['verdict']}.** T1 {'pass' if b55['T1_invariance_at_145C_and_off_sulfur']['pass'] else 'FAIL'} "
        f"({b55['T1_invariance_at_145C_and_off_sulfur']['rows_checked']} rows checked, "
        f"{len(b55['T1_invariance_at_145C_and_off_sulfur']['failures'])} moved); "
        f"T2 {'pass' if b55['T2_widening_away_from_145C']['pass'] else 'FAIL'} (median Δwidth "
        f"{_f(b55['T2_widening_away_from_145C']['median_delta_width'])} dex, min "
        f"{_f(b55['T2_widening_away_from_145C']['min_delta_width'])}); "
        f"T3 {'pass' if b55['T3_no_coverage_loss']['pass'] else 'FAIL'} "
        f"({b55['T3_no_coverage_loss']['hits_A']} → {b55['T3_no_coverage_loss']['hits_B']} of {b55['T3_no_coverage_loss']['total']}).",
        "",
        f"**B56: {b56['verdict']}** (arm {b56['judged_arm']} vs {b56['reference_arm']}). T4 "
        f"{'pass' if b56['T4_out_of_sample_coverage_in_band']['pass'] else 'FAIL'}: "
        f"{b56['T4_out_of_sample_coverage_in_band']['hits']}/{b56['T4_out_of_sample_coverage_in_band']['total']} = "
        f"{_f(b56['T4_out_of_sample_coverage_in_band']['rate'])} against [0.70, 0.95]; median width "
        f"{_f(b56['T4_out_of_sample_coverage_in_band']['median_width_log10'])} dex (reference "
        f"{_f(b56['T4_out_of_sample_coverage_in_band']['reference_rate'])}, "
        f"{_f(b56['T4_out_of_sample_coverage_in_band']['reference_median_width_log10'])} dex). σ_lane: "
        + ", ".join(f"{k} {v:.3f}" for k, v in b56["sigma_log10_by_lane"].items()) + " dex.",
        "",
        "## Out-of-sample coverage by arm and lane (fixed denominator)",
        "",
        "| arm | hits / total | rate | median width | "
        + " | ".join(sorted(cov["A"]["by_lane"])) + " |",
        "|---|---|---|---|" + "---|" * len(cov["A"]["by_lane"]),
    ]
    for name in ("A", "B", "C", "D"):
        c = cov[name]
        lanes = " | ".join(
            f"{c['by_lane'][l]['hits']}/{c['by_lane'][l]['total']} ({_f(c['by_lane'][l]['median_width_log10'], 2)})"
            for l in sorted(c["by_lane"])
        )
        out.append(f"| {name} | {c['hits']}/{c['total']} | {_f(c['rate'])} | {_f(c['median_width_log10'])} | {lanes} |")
    out += [
        "",
        "Artifact's own definitions (hits, total, not evaluable):",
        "",
    ]
    for name, d in p["artifact_own_definitions"].items():
        out.append(f"- {name}: honest literature {d['honest_literature']} (median width "
                   f"{_f(d['honest_literature_median_width'])}), out of sample {d['out_of_sample']}")
    out += [
        "",
        "## T2 rows: sulfur away from 145 °C (A → B)",
        "",
        "| benchmark | compound | T (°C) | width A | width B | Δ | p50 A | p50 B | measured | in A | in B |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in b55["T2_widening_away_from_145C"]["rows"]:
        out.append(
            f"| {r['benchmark_id']} | {r['compound']} | {_f(r['T_C'], 0)} | {_f(r['width_A'])} | {_f(r['width_B'])} | "
            f"{_f(r['delta_width'])} | {_f(r['p50_A'])} | {_f(r['p50_B'])} | {_f(r['measured'])} | "
            f"{'yes' if r['inside_A'] else 'no'} | {'yes' if r['inside_B'] else 'no'} |"
        )
    if b55["T1_invariance_at_145C_and_off_sulfur"]["failures"]:
        out += ["", "## T1 failures", ""] + [f"- {x}" for x in b55["T1_invariance_at_145C_and_off_sulfur"]["failures"]]
    out.append("")
    return "\n".join(out)


def main() -> int:
    payload = build()
    json_path, md_path = artifact_io.write_artifact(payload, OUT, render=render, sort_keys=True)
    print(f"wrote {json_path}\nwrote {md_path}")
    print(f"B55: {payload['b55']['verdict']}   B56: {payload['b56']['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
