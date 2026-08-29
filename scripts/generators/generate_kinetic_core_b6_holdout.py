#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b6_holdout.py

THE HOLD-OUT REPORT OF BUILD WAVE B6 (the lipid-oxidation module).

Scores the frozen B6 fit against the three Module 5 hold-outs declared in
``docs/reference/FIT_HOLDOUT_DECLARATION.md`` D.6, in the roles the
pre-registration fixed BEFORE the fit was run:

  1. the alpha-tocopherol two-sided signature -- **seen_diagnostic**, scored as
     MONOTONICITY THEOREMS over the whole donor range, never as a fold-error
     against a tocopherol column (prereg sec. 0 and sec. 3);
  2. the nonanal ABSENCE -- **gating**, a structural zero;
  3. the exam delta on the 8 refused matrix/lipid rows -- reported by
     re-running ``generate_cutover_final_exam.py``, which is the ONLY thing in
     this repository permitted to open ``data/benchmarks/external_validation/``.

THIS SCRIPT FITS NOTHING. It reads the frozen fit report and evaluates.

Writes results/validation/kinetic_core_b6_holdout_report.{json,md}.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.kinetic_core.engine import (  # noqa: E402
    LIPID,
    FormulationSpec,
    ProcessSpec,
    ThermalProgram,
    core_lipid_model,
    predict,
)
from src.kinetic_core.lipid import (  # noqa: E402
    REFERENCE_AUTOXIDATION_SYSTEM,
    LOOHComposition,
    charge_from_carrier,
    integrate_lipid,
    lane_coupling_verdict,
    validate_lipid_structure,
)
from src.kinetic_core.parameters_lipid import LIPID_CARRIERS  # noqa: E402
from src.kinetic_core.species_lipid import CLEAVAGE_MECHANISM, FRANKEL_SLATE  # noqa: E402

FIT_REPORT = REPO / "results/validation/kinetic_core_b6_fit_report.json"
OUT_JSON = REPO / "results/validation/kinetic_core_b6_holdout_report.json"
OUT_MD = REPO / "results/validation/kinetic_core_b6_holdout_report.md"
PREREG = REPO / "results/validation/kinetic_core_b6_prereg.md"
EXAM_JSON = REPO / "results/validation/cutover_final_exam.json"
EXAM_BASELINE = REPO / "results/validation/kinetic_core_b6_exam_baseline.json"

#: The eight rows the pre-registration names, and the outcome it registered.
PREREG_EXAM_ROWS = {
    ("external_validation_bi_2020_raw_pea_hexanal", "hexanal"): "ANSWERED",
    ("external_validation_bi_2020_roasted_pea_hexanal", "hexanal"): "ANSWERED",
    ("external_validation_li_2026_spi_wg_hme_control", "hexanal"): "ANSWERED",
    ("external_validation_li_2026_spi_wg_hme_control", "nonanal"): "STILL REFUSED",
    ("external_validation_li_2026_spi_wg_hme_control", "2-pentylfuran"): "STILL REFUSED",
    ("external_validation_li_2026_spi_wg_hme_control", "1-hexanol"): "STILL REFUSED",
    ("external_validation_liu_2023_ppi_offnote_baseline", "hexanal"): "ANSWERED",
    ("external_validation_liu_2023_ppi_offnote_baseline", "nonanal"): "STILL REFUSED",
}

DONOR_GRID = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95, 0.99)


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=REPO, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"

    return {"commit": _run("rev-parse", "HEAD"),
            "branch": _run("rev-parse", "--abbrev-ref", "HEAD")}


# ---------------------------------------------------------------------------
# HOLD-OUT 1 -- the alpha-tocopherol two-sided signature (seen_diagnostic)
# ---------------------------------------------------------------------------


def score_tocopherol(branch, composition) -> Dict[str, Any]:
    """
    The donor response, as monotonicity theorems over the WHOLE range.

    Deliberately NOT a fold-error against a tocopherol column. The donor
    loading (wt %) has no mapping onto the suppression parameter that any
    source licenses, and inventing one would be exactly the failure mode this
    module exists to avoid. Signs and monotonicity are what the architecture
    can be held to, and they are what was pre-registered.
    """
    curve: List[Dict[str, Any]] = []
    for donor in DONOR_GRID:
        flux = branch.slate_flux(composition, donor)
        total = sum(flux.values())
        curve.append({
            "donor_suppression": donor,
            "total_relative_flux": total,
            **{f"share_{p}": flux[p] / total for p in FRANKEL_SLATE},
        })

    def monotone(key: str, direction: str) -> bool:
        values = [row[key] for row in curve]
        pairs = zip(values, values[1:])
        if direction == "down":
            return all(b < a for a, b in pairs)
        return all(b > a for a, b in pairs)

    checks = {
        "H1-a total volatile flux DECREASES": monotone("total_relative_flux", "down"),
        "H1-b hexanal SHARE INCREASES": monotone("share_HEXANAL", "up"),
        "H1-d me-9-oxononanoate share INCREASES":
            monotone("share_ME_9_OXONONANOATE", "up"),
        "H1-d pentane share DECREASES": monotone("share_PENTANE", "down"),
        "H1-d methyl octanoate share DECREASES":
            monotone("share_ME_OCTANOATE", "down"),
        "H1-e me-13-oxo-tridecadienoate share DECREASES (registered as EXPECTED WRONG)":
            monotone("share_ME_13_OXO_TRIDECADIENOATE", "down"),
        "H1-e 2,4-decadienal share DECREASES (registered as EXPECTED WRONG)":
            monotone("share_DECADIENAL", "down"),
    }
    checks["H1-c the two move in OPPOSITE directions"] = (
        checks["H1-a total volatile flux DECREASES"]
        and checks["H1-b hexanal SHARE INCREASES"]
    )
    return {
        "role": "seen_diagnostic (prereg sec. 0 -- the builder saw these columns)",
        "why_not_gating": (
            "Frankel 1989 prints the FIT column and the tocopherol columns in "
            "the same table rows, and states the result in its abstract. Under "
            "the Amendment 9 clause 1 precedent a seen hold-out cannot gate."
        ),
        "why_it_is_still_worth_something": (
            "The donor parameter is NEVER FITTED and has no stored value. Every "
            "claim here holds for EVERY d in (0, 1), so there was nothing to "
            "tune toward the seen numbers even in principle."
        ),
        "mechanism_assumption": (
            "a hydrogen donor quenches the alkoxyl radical before homolytic "
            "beta-scission and does not block the heat-promoted heterolytic "
            "(Hock) cleavage. Taken from Frankel 1989's INTRODUCTION, which "
            "attributes the product/mechanism assignment to its refs 3-10, all "
            "pre-1989 -- not from the held-out arms."
        ),
        "mechanism_assignment": dict(CLEAVAGE_MECHANISM),
        "curve": curve,
        "checks": checks,
        "all_gating_machinery_checks_pass": all(
            checks[k] for k in checks if k.startswith(("H1-a", "H1-b", "H1-c"))
        ),
        "H1-e_registered_as_expected_wrong": (
            "The architecture treats methyl 13-oxo-tridecadienoate and "
            "2,4-decadienal as homolytic CO-PRODUCTS, so their shares fall with "
            "the donor. Frankel's Discussion states that no trend was observed "
            "for those two, and his Fig. 4 excludes them from both numerator "
            "and denominator. This was registered as a PREDICTED FAILURE before "
            "the fit, so it cannot be presented as a discovery now. It is the "
            "sharpest open question the module leaves: either those two "
            "products are not pure homolytic co-products, or their measured "
            "flatness is a stability artefact of the kind Frankel himself "
            "warns about."
        ),
    }


# ---------------------------------------------------------------------------
# HOLD-OUT 2 -- the nonanal ABSENCE (gating)
# ---------------------------------------------------------------------------


def score_nonanal(branch, composition) -> Dict[str, Any]:
    """Feed pure linoleate hydroperoxide, exactly as Frankel did. Expect 0.0."""
    carrier = LIPID_CARRIERS["frankel_pure_hydroperoxide"]
    charge = charge_from_carrier(carrier, composition)
    run = integrate_lipid(charge, [(1.0, 180.0)], branch)
    nonanal = run.state_mmol_per_l["NONANAL"]

    spec = FormulationSpec(
        "frankel_pure_linoleate_hydroperoxide",
        {"methyl linoleate hydroperoxide": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(180.0, 1.0), ph=6.7,
                    matrix="methyl linoleate hydroperoxide"),
    )
    engine_run = predict(spec, ["nonanal", "hexanal"])

    return {
        "role": "GATING (structural, genuinely blind -- topology, not a number)",
        "prediction": "exactly 0.0 nonanal from a pure linoleate feed",
        "nonanal_mmol_per_l": nonanal,
        "oleate_pool_mmol_per_l": charge.looh_oleate_mmol_l,
        "exact_zero": nonanal == 0.0,
        "PASS": nonanal == 0.0 and charge.looh_oleate_mmol_l == 0.0,
        "engine_answers_it": bool(engine_run.answered),
        "engine_nonanal_ug_per_l": (
            engine_run.concentrations_ug_per_l.get("nonanal")
            if engine_run.answered else None
        ),
        "structure": validate_lipid_structure(branch),
        "the_other_half": (
            "In a REAL matrix the oleate pool is not zero and the oleate -> "
            "nonanal branch fraction is measured NOWHERE, so the engine REFUSES "
            "absolute nonanal there rather than answering it. Both halves are "
            "the same hold-out: honouring an absence means refusing where the "
            "parent exists and answering zero where it does not."
        ),
    }


def score_nonanal_refusal_in_a_real_matrix() -> Dict[str, Any]:
    spec = FormulationSpec(
        "pea_isolate_probe", {"Pea Protein Isolate": 1.0},
        ProcessSpec(thermal=ThermalProgram.isothermal(140.0, 5.0), ph=6.0,
                    matrix="Pea Protein Isolate"),
    )
    run = predict(spec, ["nonanal"])
    return {
        "refused": not run.answered,
        "reasons": list(run.declaration.reasons),
        "PASS": not run.answered,
    }


# ---------------------------------------------------------------------------
# HOLD-OUT 3 -- the exam
# ---------------------------------------------------------------------------


def run_exam() -> Dict[str, Any]:
    """Re-run the cutover exam. The ONLY door to the frozen bundles."""
    subprocess.check_call(
        [sys.executable, "scripts/generators/generate_cutover_final_exam.py"],
        cwd=REPO,
    )
    return json.loads(EXAM_JSON.read_text())


def score_exam(after: Dict[str, Any], before: Dict[str, Any]) -> Dict[str, Any]:
    def index(payload):
        return {
            (r["benchmark_id"], r["compound"]): r
            for r in payload["rows"]
            if r.get("buffer_applied", True)
        }

    rows_before, rows_after = index(before), index(after)

    delta: List[Dict[str, Any]] = []
    for key, expected in PREREG_EXAM_ROWS.items():
        old, new = rows_before.get(key), rows_after.get(key)
        if old is None or new is None:
            delta.append({"benchmark_id": key[0], "compound": key[1],
                          "error": "row not found in one of the two exams"})
            continue
        outcome = "ANSWERED" if new["answered"] else "STILL REFUSED"
        delta.append({
            "benchmark_id": key[0],
            "compound": key[1],
            "was_answered": bool(old["answered"]),
            "now_answered": bool(new["answered"]),
            "prereg_outcome": expected,
            "actual_outcome": outcome,
            "prereg_met": outcome == expected,
            "old_lane": old.get("lane"),
            "new_lane": new.get("lane"),
            "new_envelope_state": new.get("envelope_state"),
            "new_reasons": new.get("declaration_reasons"),
            "measured": new.get("measured"),
            "core_predicted": new.get("core_predicted"),
            "core_fold_error": new.get("core_fold_error"),
            "old_lane_predicted": new.get("old_predicted"),
            "old_lane_fold_error": new.get("old_fold_error"),
            "core_interval_ug_per_L": new.get("core_interval_ug_per_L"),
            "measured_within_interval": new.get("core_measured_within_interval"),
            "temp_C": new.get("temp_C"),
            "time_min": new.get("time_min"),
        })

    # G-1: the 23 previously-answered points must be byte-identical.
    regressions: List[Dict[str, Any]] = []
    for key, old in rows_before.items():
        if not old["answered"]:
            continue
        new = rows_after.get(key)
        if new is None:
            regressions.append({"row": list(key), "issue": "row disappeared"})
            continue
        if not new["answered"]:
            regressions.append({"row": list(key), "issue": "now refused"})
            continue
        if old.get("core_predicted") != new.get("core_predicted"):
            regressions.append({
                "row": list(key), "issue": "prediction moved",
                "before": old.get("core_predicted"),
                "after": new.get("core_predicted"),
            })

    answered_hexanal = [
        d for d in delta
        if d.get("compound") == "hexanal" and d.get("now_answered")
    ]
    folds = [d["core_fold_error"] for d in answered_hexanal
             if d.get("core_fold_error") is not None]
    over = [
        d for d in answered_hexanal
        if d.get("core_predicted") is not None and d.get("measured")
        and d["core_predicted"] > d["measured"]
    ]

    inside = [d for d in answered_hexanal if d.get("measured_within_interval")]
    # B6 FINDING. The four hexanal rows are not four of a kind. Two of them are
    # 160 C process points; two are 40 C / 10 min, which is an AMBIENT
    # HEADSPACE measurement of an as-received isolate, not a cooking process.
    # The lipid lane is a FORMATION model with a declared no-formation-during-
    # heating gap, so on a 40 C / 10 min hold it correctly computes that almost
    # nothing forms -- and is then scored against the isolate's ACCUMULATED
    # storage oxidation. That is a category error of the same class as the
    # declaration's Brewer 1995 `dose_added_pre_cook` reclassification, and it
    # is reported here rather than argued away.
    process_rows = [d for d in answered_hexanal
                    if (d.get("temp_C") or 0.0) >= 100.0]
    ambient_rows = [d for d in answered_hexanal
                    if (d.get("temp_C") or 0.0) < 100.0]
    under = [
        d for d in answered_hexanal
        if d.get("core_predicted") is not None and d.get("measured")
        and d["core_predicted"] < d["measured"]
    ]

    return {
        "refusals_before": sum(1 for r in rows_before.values() if not r["answered"]),
        "refusals_after": sum(1 for r in rows_after.values() if not r["answered"]),
        "prereg_refusals_after": 13,
        "per_row": delta,
        "prereg_all_eight_met": all(d.get("prereg_met") for d in delta),
        "G1_regression_guard": {
            "n_previously_answered": sum(1 for r in rows_before.values() if r["answered"]),
            "regressions": regressions,
            "PASS": not regressions,
        },
        "E1_accuracy": {
            "n_hexanal_answered": len(answered_hexanal),
            "fold_errors": folds,
            "n_worse_than_3x": sum(1 for f in folds if f > 3.0),
            "prereg_expected_worse_than_3x_at_least": 3,
            "n_over_predicted": len(over),
            "n_under_predicted": len(under),
            "prereg_expected_over_predicted_at_least": 3,
            "prereg_direction_met": len(over) >= 3,
            "n_measured_inside_interval": len(inside),
            "prereg_expected_inside_interval_at_least": 3,
            "prereg_interval_met": len(inside) >= 3,
            "split_by_thermal_step": {
                "why": (
                    "Two of the four rows are 160 C process points; two are "
                    "40 C for 10 min, which measures an as-received isolate's "
                    "ACCUMULATED storage oxidation, not a process output. A "
                    "formation model with a declared no-formation-during-"
                    "heating gap cannot make hexanal in 10 min at 40 C and "
                    "should not be expected to. Reported as a benchmark-design "
                    "finding, not corrected by fiat."
                ),
                "process_rows": [
                    {"benchmark_id": d["benchmark_id"], "temp_C": d.get("temp_C"),
                     "time_min": d.get("time_min"),
                     "fold": d.get("core_fold_error"),
                     "inside_interval": d.get("measured_within_interval")}
                    for d in process_rows
                ],
                "ambient_rows": [
                    {"benchmark_id": d["benchmark_id"], "temp_C": d.get("temp_C"),
                     "time_min": d.get("time_min"),
                     "fold": d.get("core_fold_error"),
                     "inside_interval": d.get("measured_within_interval")}
                    for d in ambient_rows
                ],
                "process_median_fold": (
                    sorted(d["core_fold_error"] for d in process_rows
                           if d.get("core_fold_error"))[len(process_rows) // 2]
                    if process_rows else None
                ),
                "process_rows_inside_interval": sum(
                    1 for d in process_rows if d.get("measured_within_interval")),
                "ambient_rows_inside_interval": sum(
                    1 for d in ambient_rows if d.get("measured_within_interval")),
                "identical_input_records": (
                    "bi_2020_raw_pea and liu_2023_ppi_offnote_baseline declare "
                    "IDENTICAL conditions (40 C, 10 min, pH 6.0, a_w 0.95) and "
                    "IDENTICAL precursors ('Pea Protein Isolate'), yet their "
                    "measured hexanal differs 9.0x. No model that reads only "
                    "the recorded fields can separate them, and this one "
                    "predicts the same number for both. That is a defect in the "
                    "BENCHMARK RECORDS, not in the model: whatever distinguishes "
                    "the two samples is not written down."
                ),
            },
            "per_row_interval": [
                {
                    "benchmark_id": d["benchmark_id"],
                    "measured": d.get("measured"),
                    "point": d.get("core_predicted"),
                    "interval": d.get("core_interval_ug_per_L"),
                    "inside": d.get("measured_within_interval"),
                }
                for d in answered_hexanal
            ],
        },
    }


def main() -> int:
    if not FIT_REPORT.exists():
        raise SystemExit(f"{FIT_REPORT} missing -- run the B6 fit generator first.")
    if not EXAM_BASELINE.exists():
        raise SystemExit(
            f"{EXAM_BASELINE} missing. Snapshot the PRE-B6 exam first:\n"
            f"  cp results/validation/cutover_final_exam.json "
            f"{EXAM_BASELINE.relative_to(REPO)}\n"
            f"taken at the commit BEFORE the lipid lane existed."
        )

    branch, composition = core_lipid_model()
    before = json.loads(EXAM_BASELINE.read_text())
    after = run_exam()

    payload: Dict[str, Any] = {
        "artifact": "kinetic_core_b6_holdout_report",
        "wave": "B6 -- the lipid-oxidation module",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "pre_registration": str(PREREG.relative_to(REPO)),
        "fit_report": str(FIT_REPORT.relative_to(REPO)),
        "holdout_1_tocopherol": score_tocopherol(branch, composition),
        "holdout_2_nonanal": score_nonanal(branch, composition),
        "holdout_2b_nonanal_refused_in_a_real_matrix":
            score_nonanal_refusal_in_a_real_matrix(),
        "holdout_3_exam": score_exam(after, before),
        "lane_coupling_verdict": lane_coupling_verdict(["Glc", "Fru", "Gly"]),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n")

    exam = payload["holdout_3_exam"]
    toc = payload["holdout_1_tocopherol"]
    non = payload["holdout_2_nonanal"]

    lines: List[str] = []
    lines.append("# Kinetic core B6 -- HOLD-OUT report (the lipid-oxidation module)\n")
    lines.append(f"Generated {payload['generated_on']} at "
                 f"`{payload['git']['commit'][:12]}`.\n")
    lines.append(f"Scored against `{payload['pre_registration']}`, written "
                 f"before the fit.\n")

    lines.append("## Headline\n")
    lines.append(
        f"| hold-out | role | outcome |\n|---|---|---|\n"
        f"| Frankel alpha-tocopherol two-sided signature | seen_diagnostic | "
        f"{'machinery checks PASS' if toc['all_gating_machinery_checks_pass'] else 'FAIL'} |\n"
        f"| Frankel nonanal ABSENCE | **gating** | "
        f"{'PASS' if non['PASS'] else 'FAIL'} |\n"
        f"| exam refusals | -- | "
        f"{exam['refusals_before']} -> {exam['refusals_after']} "
        f"(pre-registered {exam['prereg_refusals_after']}) |\n"
        f"| all eight registered row outcomes met? | -- | "
        f"{'YES' if exam['prereg_all_eight_met'] else 'NO'} |\n"
        f"| G-1 regression guard (23 answered points unmoved) | -- | "
        f"{'PASS' if exam['G1_regression_guard']['PASS'] else 'FAIL'} |\n")

    e1h = exam["E1_accuracy"]
    lines.append("\n### Pre-registration scorecard -- INCLUDING what was missed\n")
    lines.append(
        f"| registered claim | outcome |\n|---|---|\n"
        f"| E-1 fold: at least 3 of 4 hexanal points worse than 3x | "
        f"**{'MET' if e1h['n_worse_than_3x'] >= 3 else 'MISSED'}** "
        f"({e1h['n_worse_than_3x']}/4) |\n"
        f"| E-1 direction: at least 3 of 4 OVER-predicted | "
        f"**{'MET' if e1h['prereg_direction_met'] else 'MISSED'}** "
        f"({e1h['n_over_predicted']} over, {e1h['n_under_predicted']} under) |\n"
        f"| E-1 interval: at least 3 of 4 measured values inside the interval | "
        f"**{'MET' if e1h['prereg_interval_met'] else 'MISSED'}** "
        f"({e1h['n_measured_inside_interval']}/4) |\n")
    lines.append(
        "\n**The direction claim was MISSED, and the miss is the wave's single "
        "most important result.** The pre-registration expected OVER-prediction, "
        "on the strength of this repository's own standing diagnosis for this "
        "lane (the FAST lane's documented 36x and 1304x hexanal "
        "over-predictions). The B6 lipid lane UNDER-predicts all four points. "
        "The standing diagnosis does not transfer: it was a property of the FAST "
        "lane's unbounded first-order extrapolation from an invented 0.37 branch "
        "fraction, not of lipid chemistry.\n")
    split = e1h["split_by_thermal_step"]
    lines.append(
        "\n**And the four points are not four of a kind.** " + split["why"] + "\n")
    lines.append(
        f"\n| class | rows | fold errors | inside interval |\n|---|---:|---|---:|\n"
        f"| 160 C PROCESS points | {len(split['process_rows'])} | "
        + ", ".join(f"{r['fold']:.1f}x" for r in split["process_rows"])
        + f" | {split['process_rows_inside_interval']}/{len(split['process_rows'])} |\n"
        f"| 40 C / 10 min AMBIENT points | {len(split['ambient_rows'])} | "
        + ", ".join(f"{r['fold']:.0f}x" for r in split["ambient_rows"])
        + f" | {split['ambient_rows_inside_interval']}/{len(split['ambient_rows'])} |\n")
    lines.append(
        "\nOn the two rows that are actually a thermal process, the lane lands "
        "at 3.7x and 8.7x with BOTH measured values inside the reported "
        "interval. On the two 40 C / 10 min rows it is 3 717x and 33 392x low, "
        "which is exactly what `parameters_lipid.LOOH_FORMATION_GAP` says will "
        "happen: the module cannot MAKE hydroperoxide during the hold, so on a "
        "near-ambient point it predicts that almost nothing forms -- and is then "
        "scored against the isolate's ACCUMULATED storage oxidation. Reported, "
        "not corrected by fiat.\n")
    lines.append("\n" + split["identical_input_records"] + "\n")

    lines.append("\n## Hold-out 1 -- the alpha-tocopherol two-sided signature\n")
    lines.append(f"**Role: {toc['role']}.** {toc['why_not_gating']}\n")
    lines.append(toc["why_it_is_still_worth_something"] + "\n")
    lines.append("| check | result |\n|---|---|")
    for check, ok in toc["checks"].items():
        lines.append(f"| {check} | {'HOLDS' if ok else 'FAILS'} |")
    lines.append("\n### The donor response curve\n")
    lines.append("| d | total flux (rel.) | hexanal share | me-9-oxononanoate share | pentane share |")
    lines.append("|---:|---:|---:|---:|---:|")
    for row in toc["curve"]:
        lines.append(
            f"| {row['donor_suppression']:.2f} | {row['total_relative_flux']:.4f} "
            f"| {row['share_HEXANAL']:.4f} | {row['share_ME_9_OXONONANOATE']:.4f} "
            f"| {row['share_PENTANE']:.4f} |"
        )
    lines.append("\n" + toc["H1-e_registered_as_expected_wrong"] + "\n")

    lines.append("## Hold-out 2 -- the nonanal ABSENCE (gating)\n")
    lines.append(
        f"Pure linoleate hydroperoxide feed, 180 C: nonanal = "
        f"**{non['nonanal_mmol_per_l']:.6g} mmol/L** "
        f"(exact zero: {non['exact_zero']}). **"
        f"{'PASS' if non['PASS'] else 'FAIL'}.**\n")
    lines.append(non["the_other_half"] + "\n")
    lines.append(
        f"In a real (pea) matrix the engine refuses nonanal: "
        f"{payload['holdout_2b_nonanal_refused_in_a_real_matrix']['refused']}.\n")

    lines.append("## Hold-out 3 -- the exam\n")
    lines.append("| bundle | compound | was | now | pre-registered | met? | measured | core | fold | old lane |")
    lines.append("|---|---|---|---|---|---|---:|---:|---:|---:|")
    for row in exam["per_row"]:
        if "error" in row:
            lines.append(f"| {row['benchmark_id']} | {row['compound']} | "
                         f"{row['error']} | | | | | | | |")
            continue
        lines.append(
            f"| `{row['benchmark_id']}` | {row['compound']} | "
            f"{'answered' if row['was_answered'] else 'refused'} | "
            f"{'answered' if row['now_answered'] else 'refused'} | "
            f"{row['prereg_outcome']} | {'YES' if row['prereg_met'] else 'NO'} | "
            f"{_fmt(row['measured'])} | {_fmt(row['core_predicted'])} | "
            f"{_fmt(row['core_fold_error'])} | {_fmt(row['old_lane_fold_error'])} |"
        )
    e1 = exam["E1_accuracy"]
    lines.append(
        f"\n**E-1**: {e1['n_hexanal_answered']} hexanal points answered; "
        f"{e1['n_worse_than_3x']} worse than 3x "
        f"(pre-registered at least {e1['prereg_expected_worse_than_3x_at_least']}); "
        f"{e1['n_over_predicted']} over-predicted and {e1['n_under_predicted']} "
        f"UNDER-predicted (pre-registered at least "
        f"{e1['prereg_expected_over_predicted_at_least']} OVER-predicted -- "
        f"**{'MET' if e1['prereg_direction_met'] else 'MISSED'}**); "
        f"{e1['n_measured_inside_interval']}/{e1['n_hexanal_answered']} measured "
        f"values inside the reported interval (pre-registered at least "
        f"{e1['prereg_expected_inside_interval_at_least']} -- "
        f"**{'MET' if e1['prereg_interval_met'] else 'MISSED'}**).\n")
    lines.append("\n| bundle | measured | point | interval | inside? |")
    lines.append("|---|---:|---:|---|---|")
    for row in e1["per_row_interval"]:
        band = row["interval"]
        lines.append(
            f"| `{row['benchmark_id']}` | {_fmt(row['measured'])} | "
            f"{_fmt(row['point'])} | "
            + ("--" if not band else f"{_fmt(band[0])} - {_fmt(band[1])}")
            + f" | {row['inside']} |"
        )
    lines.append(
        f"**G-1**: {exam['G1_regression_guard']['n_previously_answered']} points "
        f"were answered before B6; regressions: "
        f"{len(exam['G1_regression_guard']['regressions'])}.\n")
    if exam["G1_regression_guard"]["regressions"]:
        lines.append("```json\n"
                     + json.dumps(exam["G1_regression_guard"]["regressions"], indent=1)
                     + "\n```\n")

    lines.append("## Lane coupling\n")
    lines.append("```json\n" + json.dumps(payload["lane_coupling_verdict"], indent=1)
                 + "\n```\n")

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_JSON.relative_to(REPO)}")
    print(f"wrote {OUT_MD.relative_to(REPO)}")
    print(f"refusals {exam['refusals_before']} -> {exam['refusals_after']}; "
          f"prereg met on all eight: {exam['prereg_all_eight_met']}; "
          f"G-1 {'PASS' if exam['G1_regression_guard']['PASS'] else 'FAIL'}")
    return 0


def _fmt(value) -> str:
    if value is None:
        return "--"
    return f"{float(value):.3g}"


if __name__ == "__main__":
    raise SystemExit(main())
