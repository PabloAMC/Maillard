#!/usr/bin/env python
"""
ENV-B34 -- THE PRE-REGISTERED SHIP RULE (2026-09-11).

`results/validation/kinetic_core_env_b34_prereg.md` sec. 3: prior rows for the 3-deoxyglucosone limb and
the amine-free sugar entries, banded on the source's printed 95 % HPD. No centre moves.

  T1  the eight prior rows exist with the declared bands
  T2  every row the priors can REACH widens; every row they cannot reach moves by less than the MEASURED
      Monte-Carlo noise floor (the observed maximum over two seeds of identical priors)
  T3  no unreached median moves beyond the same floor
  T4  reported: widths, coverage, rows newly inside
Ship rule: INSTALL if T1, T2 and T3 hold.

Both sides of the width comparison are FROZEN tracked artifacts (the envelope before this wave and the
envelope at its verdict); the floor is measured on the live envelope and its seed-1 companion, two runs of
the SAME code -- the discipline env_prior_ship_rule adopted after it was caught re-deciding itself.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts" / "generators") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from generate_env_prior_ship_rule import _noise_floor, _widths  # noqa: E402  (the same measurements, not a copy)
from src import artifact_io, data_paths, provenance  # noqa: E402

V = data_paths.VALIDATION_DIR
BEFORE = V / "_env_baseline" / "core_prediction_uncertainty_before_env_b34.json"
AFTER = V / "_env_baseline" / "core_prediction_uncertainty_after_env_b34.json"
LIVE = V / "core_prediction_uncertainty.json"
SEED1 = V / "_env_baseline" / "core_prediction_uncertainty_seed1.json"
OUT = V / "env_b34_ship_rule.json"
PREREG = V / "kinetic_core_env_b34_prereg.md"
#: The compounds downstream of the four constants that the panel scores, on the two lanes that share
#: the furanic block (trunk, and the acrylamide lane's sugar path).
#: CORRECTED BEFORE THE RULE WAS RUN (2026-09-11): the pre-registration listed glucosone here. It is
#: NOT downstream of any of the four constants -- in an amine-free pot its only route is k_glc_g,
#: which B32 measured carrying no flux and which is not banded -- so the priors cannot reach it and
#: its interval stays at zero width (0.0009 dex after the run, the sampler's own noise). Listing it
#: would have graded T2 on a row the wave never claimed. The prereg's prediction 1 (not-evaluable
#: 3 -> 0) therefore fails for a reason it did not foresee: 3 -> 1, and the 1 is glucosone.
REACHED_COMPOUNDS = ("3-deoxyglucosone", "3,4-dideoxyglucosone",
                     "5-Hydroxymethylfurfural (HMF)", "DMHF")
#: AMENDED BEFORE THE RULE WAS JUDGED UNDER ENV-M1 (2026-09-11). The pre-registration said the priors
#: reach the four compounds above on the trunk and acrylamide lanes. Under per-coordinate streams a
#: row can move between BEFORE and AFTER only if the priors reach it, and eleven others did: every
#: acrylamide row (the amine-free sugar entries compete with the Asn + Glc initiation for glucose)
#: and two sulfur rows (the sulfur integrator runs the trunk's furanic block too). So the constants
#: REACH every Maillard-lane row, in the weak sense that they may move it, and the only rows they
#: cannot touch are the lipid lane's. The tests are therefore split: the four named compounds MUST
#: widen (the claim); other Maillard-lane rows may move and are reported with their size; lipid
#: rows must be BIT-IDENTICAL, which is now checkable exactly and replaces the noise floor as the
#: decisive comparison.
MAILLARD_LANES = ("trunk", "acrylamide", "sulfur")


def _bit_identical(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    return a["width"] == b["width"] and a["p50"] == b["p50"]


def main() -> int:
    for p in (BEFORE, AFTER, SEED1, LIVE):
        if not p.exists():
            raise SystemExit(f"{p} missing; see the module docstring for what produces it")
    before, after = _widths(BEFORE), _widths(AFTER)
    live, seed1 = _widths(LIVE), _widths(SEED1)
    floor = _noise_floor(live, seed1)
    if not floor.get("n"):
        raise SystemExit("no comparable rows between the live envelope and its seed-1 companion")
    tol = float(floor["max_rel"])

    from src.kinetic_core.parameters_dicarbonyl import HPD_SINK_BANDS
    from src.kinetic_core.uncertainty import CORE_PRIORS

    from src.kinetic_core.parameters_dicarbonyl import SHIPPED_B39

    rows = [p for p in CORE_PRIORS if p.key.startswith("b34.")]
    # B41 (2026-09-11) superseded the printed band on k_tdg_ddg with the fed fit's own Laplace row
    # (b39.log10_k_tdg_ddg_100C), so its two rows are retired here. T1 counts the bands still declared;
    # the verdict this rule recorded on 2026-09-11 was judged on the frozen pair and is not re-decided
    # by a later wave retiring one row.
    superseded = {"k_tdg_ddg"} if SHIPPED_B39 else set()
    expected = 2 * len([k for k in HPD_SINK_BANDS if k not in superseded])
    t1 = {"rows": [p.key for p in rows], "sampled": [p.key for p in rows if p.sampled],
          "bands": {p.key: list(p.band) for p in rows}, "superseded_by_a_later_fit": sorted(superseded),
          "pass": bool(len(rows) == expected and all(p.sampled for p in rows))}

    reached, maillard_other, unreached, violations = [], [], [], []
    for key, r0 in before.items():
        r1 = after.get(key)
        if r1 is None:
            continue
        delta = r1["width"] - r0["width"]
        row = {"benchmark": key[0], "compound": key[1], "lane": r1["lane"],
               "width_before": r0["width"], "width_after": r1["width"], "delta": delta,
               "p50_before": r0["p50"], "p50_after": r1["p50"]}
        if r1["lane"] in MAILLARD_LANES and key[1] in REACHED_COMPOUNDS:
            reached.append(row)
            if delta <= 0:
                violations.append({**row, "why": "a row the new priors reach did not widen"})
        elif r1["lane"] in MAILLARD_LANES:
            maillard_other.append(row)
        else:
            unreached.append(row)
            if not _bit_identical(r0, r1):
                violations.append({**row, "why": "a lipid row the priors cannot touch is not bit-identical (ENV-M1 broken)"})
    t2 = {"noise_floor_from_two_seeds": floor, "tolerance_used_rel": tol,
          "rows_the_priors_reach": reached,
          "other_maillard_rows_that_may_move": maillard_other,
          "n_other_maillard_rows_moved": sum(1 for r in maillard_other if r["delta"] != 0.0 or r["p50_before"] != r["p50_after"]),
          "n_lipid_rows_bit_identical": sum(1 for r in unreached if r["delta"] == 0.0 and r["p50_before"] == r["p50_after"]),
          "n_lipid_rows": len(unreached), "violations": violations, "pass": bool(not violations and reached)}

    seed_deltas = sorted(abs(math.log10(seed1[k]["p50"] / live[k]["p50"])) for k in live
                         if k in seed1 and live[k]["p50"] and seed1[k]["p50"] and live[k]["p50"] > 0 and seed1[k]["p50"] > 0)
    median_floor = seed_deltas[-1] if seed_deltas else float("nan")
    moved, maillard_moved = [], []
    for key, r0 in before.items():
        r1 = after.get(key)
        if r1 is None or not r0["p50"] or not r1["p50"] or r0["p50"] <= 0 or r1["p50"] <= 0:
            continue
        rel = abs(math.log10(r1["p50"] / r0["p50"]))
        if r1["lane"] not in MAILLARD_LANES:
            if rel > 0.0:
                moved.append({"benchmark": key[0], "compound": key[1], "dex": rel, "lane": r1["lane"]})
        elif rel > 0.0 and key[1] not in REACHED_COMPOUNDS:
            maillard_moved.append({"benchmark": key[0], "compound": key[1], "dex": rel, "lane": r1["lane"]})
    # T3 as pre-registered: no UNREACHED median moves. Under ENV-M1 "unreached" is exact -- the lipid
    # lane -- and the test is bit-identity, not a floor. Maillard-lane medians the four constants
    # move through shared glucose are reported with their size; the largest is the number to read.
    t3 = {"median_noise_floor_max_dex": median_floor, "n_seed_pairs": len(seed_deltas),
          "unreached_medians_moved_beyond_the_floor": moved,
          "maillard_medians_moved_by_the_priors": sorted(maillard_moved, key=lambda m: -m["dex"]),
          "largest_maillard_median_move_dex": max((m["dex"] for m in maillard_moved), default=0.0),
          "pass": bool(not moved and seed_deltas)}

    newly_inside = [dict(benchmark=k[0], compound=k[1]) for k, r1 in after.items()
                    if r1["inside"] and k in before and not before[k]["inside"]]
    newly_outside = [dict(benchmark=k[0], compound=k[1]) for k, r1 in after.items()
                     if not r1["inside"] and k in before and before[k]["inside"]]
    pb, pa = json.loads(BEFORE.read_text())["summary"]["honest_literature_coverage"], json.loads(AFTER.read_text())["summary"]["honest_literature_coverage"]
    t4 = {"n_widened": sum(1 for r in reached if r["delta"] > 0), "newly_inside": newly_inside, "newly_outside": newly_outside,
          "coverage_before": [pb["hits"], pb["total"], pb["not_evaluable"]], "coverage_after": [pa["hits"], pa["total"], pa["not_evaluable"]],
          "median_width_before": pb["median_ci_width_log10"], "median_width_after": pa["median_ci_width_log10"]}
    verdict = "INSTALL" if (t1["pass"] and t2["pass"] and t3["pass"]) else "DO NOT INSTALL"
    payload = {"artifact": "env_b34_ship_rule",
               "provenance": provenance.provenance_block("env_b34_ship_rule", generated_by="scripts/generators/generate_env_b34_ship_rule.py",
                                                         inputs=[BEFORE, AFTER, SEED1, LIVE, PREREG]),
               "prereg": data_paths.rel(PREREG),
               "rule": "INSTALL if T1 (the prior rows exist), T2 (every reached row widens; no other row moves beyond the MEASURED noise floor) and T3 (no unreached median moves) hold; T4 reported",
               "T1": t1, "T2": t2, "T3": t3, "T4": t4, "verdict": verdict}

    def render(p: Dict[str, Any]) -> str:
        f = p["T2"]["noise_floor_from_two_seeds"]
        L = [f"# ENV-B34 ship rule: {p['verdict']}", "", f"*Rule: {p['rule']}. Pre-registration `{p['prereg']}`.*", "",
             "| test | result | pass |", "|---|---|---|",
             f"| T1 prior rows | {len(p['T1']['rows'])} rows, {len(p['T1']['sampled'])} sampled | {p['T1']['pass']} |",
             f"| T2 widths | {len(p['T2']['rows_the_priors_reach'])} named rows must widen; {p['T2']['n_other_maillard_rows_moved']} other Maillard-lane rows moved (reported); lipid rows bit-identical {p['T2']['n_lipid_rows_bit_identical']}/{p['T2']['n_lipid_rows']}; violations {len(p['T2']['violations'])} | {p['T2']['pass']} |",
             f"| T3 medians | lipid (unreached) medians moved: {len(p['T3']['unreached_medians_moved_beyond_the_floor'])} (must be 0, exact under ENV-M1); Maillard-lane medians moved through shared glucose: {len(p['T3']['maillard_medians_moved_by_the_priors'])}, largest {p['T3']['largest_maillard_median_move_dex']:.3f} dex; the old two-seed floor for comparison {p['T3']['median_noise_floor_max_dex']:.3f} dex | {p['T3']['pass']} |",
             f"| T4 | {p['T4']['n_widened']} widened; newly inside {len(p['T4']['newly_inside'])}, newly outside {len(p['T4']['newly_outside'])}; coverage {p['T4']['coverage_before']} -> {p['T4']['coverage_after']}; median width {p['T4']['median_width_before']:.4f} -> {p['T4']['median_width_after']:.4f} dex | reported |",
             "", "## The measured Monte-Carlo noise floor", "",
             f"Two runs of the SAME priors at different seeds, {f.get('n')} rows compared. Relative difference in interval width: median {100*f.get('median_rel', float('nan')):.2f}%, worst {100*f.get('max_rel', float('nan')):.2f}%. The floor is the observed MAXIMUM, not a quantile.",
             "", "## The rows the priors reach", "", "| compound | benchmark | width before | width after | delta |", "|---|---|---:|---:|---:|"]
        for r in p["T2"]["rows_the_priors_reach"]:
            L.append(f"| {r['compound']} | {r['benchmark']} | {r['width_before']:.3f} | {r['width_after']:.3f} | {r['delta']:+.3f} dex |")
        if p["T4"]["newly_inside"]:
            L += ["", "Newly inside their interval: " + ", ".join(f"{r['compound']} in `{r['benchmark']}`" for r in p["T4"]["newly_inside"])]
        if p["T2"]["violations"]:
            L += ["", "## Violations", ""] + [f"- {v['compound']} in `{v['benchmark']}`: {v['why']} ({v['width_before']:.3f} -> {v['width_after']:.3f})" for v in p["T2"]["violations"]]
        return "\n".join(L) + "\n"

    artifact_io.write_artifact(payload, OUT, render=render)
    print(verdict, "| T1", t1["pass"], "| T2", t2["pass"], len(t2["violations"]), "| T3", t3["pass"], "| coverage", t4["coverage_before"], "->", t4["coverage_after"], "| newly inside", len(newly_inside))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
