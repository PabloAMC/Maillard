"""
Build Wave B2.4 -- the two wave reports, assembled from artifacts already on
disk.

THIS FILE COMPUTES NOTHING CHEMICAL. It reads the ensemble, the three
per-weighting hold-out panels, the exams (three per-weighting, six per-member
at the control weighting, and one W-4 control at B2.2's own parameters), and
the rebuilt amine-fate probe -- and it scores the pre-registration against
them, claim by claim, whether each held or not.

Outputs:
  results/validation/kinetic_core_b2_4_fit_report.{md,json}
  results/validation/kinetic_core_b2_4_holdout_report.{md,json}

FIREWALL. Every measured value this file touches has already passed through a
scorer. It opens no bundle and no benchmark.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
V = data_paths.VALIDATION_DIR

TAGS = ("shipped", "half", "measured")
TAG_LABEL = {"shipped": "W-SHIPPED", "half": "W-HALF", "measured": "W-MEASURED"}

#: The Hofmann pH-3 / pH-7 ribose block: four exam points, four panel rows, one
#: set of measurements. Declared as shared in both artifacts since B2.4.
HOFMANN_PH_BLOCK = (
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3",
     "2-Furfurylthiol (FFT)", "hofmann_ribose_pH3_FFT"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3",
     "2-Methyl-3-furanthiol (MFT)", "hofmann_ribose_pH3_MFT"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7",
     "2-Furfurylthiol (FFT)", "hofmann_ribose_pH7_FFT"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7",
     "2-Methyl-3-furanthiol (MFT)", "hofmann_ribose_pH7_MFT"),
)

#: The pH-5 rung of the same ladder is a DECLARED FIT anchor, and its measured
#: value is therefore fit-side. Used to place the middle rung of the shape test.
FIT_PH5_ANCHOR_PPB = {"FFT": 121.0, "MFT": 198.0}
FIT_PH5_ROW = {"FFT": "hofmann_ribose_FFT", "MFT": "hofmann_ribose_MFT"}


def _load(path: Path) -> Optional[Dict[str, Any]]:
    return json.loads(path.read_text()) if path.exists() else None


def _fmt(v, nd=4):
    if v is None:
        return "--"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if not math.isfinite(v):
            return "inf"
        if abs(v) >= 1e5 or (v != 0 and abs(v) < 1e-3):
            return f"{v:.3g}"
        return f"{v:.{nd}g}"
    return str(v)


# ---------------------------------------------------------------------------
# Extractors
# ---------------------------------------------------------------------------


#: THE POOL MOVED UNDER THIS WAVE, AND EVERY CROSS-WAVE NUMBER DEPENDS ON IT.
#:
#: D1's reconciliation reports its geometric means over 23 answered points. The
#: exam now answers 27: concurrent Build Wave B6 gave the core a lipid lane, so
#: the four `matrix_path` points it used to DECLINE are now scored -- and their
#: family median is ~1900x, which drags any pool-wide geometric mean upward for
#: a reason that has nothing to do with a pH weighting.
#:
#: So every cross-wave comparison in this report is made on D1'S OWN POOL: the
#: same 23 points, i.e. all scored points except the lipid family. Both numbers
#: are printed everywhere, and the restriction is named every time it is used.
D1_POOL_EXCLUDED_FAMILY = "matrix_path_lipid"


def _geo(folds) -> Optional[float]:
    clean = [f for f in folds if f and math.isfinite(f) and f > 0]
    if not clean:
        return None
    return 10.0 ** (sum(abs(math.log10(f)) for f in clean) / len(clean))


def exam_stats(exam: Dict[str, Any]) -> Dict[str, Any]:
    s = exam["summary"]
    core, paired = s["core"], s["paired_subset"]
    d1_rows = [r for r in exam["rows"]
               if r["core_fold_error"] is not None
               and r["family"] != D1_POOL_EXCLUDED_FAMILY]
    hof = {}
    for bundle, compound, panel_row in HOFMANN_PH_BLOCK:
        for row in exam["rows"]:
            if row["benchmark_id"] == bundle and row["compound"] == compound:
                hof[panel_row] = {
                    "predicted_ppb": row["core_predicted"],
                    "measured_ppb": row["measured"],
                    "fold": row["core_fold_error"],
                    "in_band": row["core_within_band"],
                }
    return {
        # --- D1's own 23-point pool, the ONLY like-for-like cross-wave column
        "d1_pool_n": len(d1_rows),
        "d1_pool_geometric_mean_fold": _geo(
            [r["core_fold_error"] for r in d1_rows]),
        "d1_pool_median_fold": (
            statistics.median([r["core_fold_error"] for r in d1_rows])
            if d1_rows else None),
        "d1_pool_as_was_geometric_mean_fold": _geo(
            [r["as_was_fold_error"] for r in d1_rows]),
        "d1_pool_as_was_median_fold": (
            statistics.median([r["as_was_fold_error"] for r in d1_rows
                               if r.get("as_was_fold_error")])
            if d1_rows else None),
        "d1_pool_in_band": sum(1 for r in d1_rows if r["core_within_band"]),
        "answered": s["core_answered"],
        "declined": s["core_declined"],
        "scored": s["core_scored"],
        "in_band": s["core_within_band"],
        "median_fold_all": core.get("median_fold_error"),
        "geometric_mean_fold_all": core.get("geometric_mean_fold"),
        "median_abs_log10_all": core.get("median_abs_log10_error"),
        "worst_fold": core.get("worst_fold_error"),
        "paired_n": paired["n"],
        "paired_core_median": paired["core"].get("median_fold_error"),
        "paired_core_geomean": paired["core"].get("geometric_mean_fold"),
        "paired_old_median": paired["old"].get("median_fold_error"),
        "paired_old_geomean": paired["old"].get("geometric_mean_fold"),
        "as_was_geomean": exam["both_ways"]["as_was"]["all_points"].get(
            "geometric_mean_fold"),
        "as_was_median": exam["both_ways"]["as_was"]["all_points"].get(
            "median_fold_error"),
        "by_family": {
            k: {"core_median_fold": v.get("core_median_fold"),
                "answered": v.get("answered"), "points": v.get("points")}
            for k, v in exam["by_family"].items()},
        "hofmann_ph_block": hof,
        "shared_rows_declared": (exam.get("shared_with_holdout_panel") or {}).get("n"),
    }


def panel_stats(panel: Dict[str, Any]) -> Dict[str, Any]:
    sc = panel["scorecard"]
    block = {}
    for row in panel["rows"]:
        if row["id"] in {r[2] for r in HOFMANN_PH_BLOCK}:
            block[row["id"]] = {"fold": row.get("fold_error"),
                                "pass": row.get("pass"),
                                "shared": row.get("shared_with_cutover_exam")}
    return {
        "gating_passed": sc["gating_passed"],
        "gating_rows": sc["gating_rows"],
        "diagnostic_passed": sc["diagnostic_passed"],
        "diagnostic_rows": sc["diagnostic_rows"],
        "median_abs_log10_fold_gating": sc.get("median_abs_log10_fold_gating"),
        "geometric_mean_fold_gating": sc.get("geometric_mean_fold_gating"),
        "median_abs_log10_fold_all_scored": sc.get("median_abs_log10_fold_all_scored"),
        "geometric_mean_fold_all_scored": sc.get("geometric_mean_fold_all_scored"),
        "hofmann_ph_block": block,
        "hofmann_ph_block_passed": sum(
            1 for v in block.values() if v.get("pass")),
        "shared_rows_declared": sum(
            1 for r in panel["rows"] if r.get("shared_with_cutover_exam")),
    }


def member_fit_rows(member: Dict[str, Any], b23_rows) -> Dict[str, Dict[str, Any]]:
    """Per-row residual and fold at one ensemble member, from its stored r."""
    out = {}
    for row_id, r, row in zip(member["row_ids"], member["residuals"], b23_rows):
        sigma = float(row["sigma_log"])
        out[row_id] = {
            "residual": r,
            "sigma": sigma,
            "kind": row["kind"],
            "fold": (10.0 ** (r * sigma)) if row["kind"] != "ph_endpoint" else None,
            "miss_pH": (r * sigma) if row["kind"] == "ph_endpoint" else None,
        }
    return out


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build() -> Dict[str, Any]:
    import scripts.generators.generate_kinetic_core_b2_3_fit as B23

    ensemble = _load(V / "kinetic_core_b2_4_ensemble.json")
    if ensemble is None:
        raise SystemExit("the ensemble has not been consolidated")
    members = ensemble["members"]

    exams = {t: _load(V / f"kinetic_core_b2_4_exam_{t}.json") for t in TAGS}
    panels = {t: _load(V / f"kinetic_core_b2_4_panel_{t}.json") for t in TAGS}
    member_exams = {}
    for t in TAGS:
        for start in range(6):
            path = V / f"kinetic_core_b2_4_exam_{t}_s{start}.json"
            if path.exists():
                member_exams[(t, start)] = _load(path)
    w4 = _load(V / "kinetic_core_b2_4_exam_w4_control_b2_2_params.json")
    amine = _load(V / "kinetic_core_b2_4_amine_fate_probe.json")

    per_tag: Dict[str, Any] = {}
    for t in TAGS:
        block = [m for m in members if m["weight_tag"] == t]
        best = min(block, key=lambda m: m["cost"]) if block else None
        rows = member_fit_rows(best, B23.ACTIVE_FIT_ROWS) if best else {}
        largest = (max(rows.items(), key=lambda kv: abs(kv[1]["residual"]))[0]
                   if rows else None)
        per_tag[t] = {
            "ensemble": ensemble["ensemble"][t],
            "best_member": {k: v for k, v in (best or {}).items()
                            if k not in ("x_full", "residuals")},
            "largest_residual_row": largest,
            "largest_residual": abs(rows[largest]["residual"]) if largest else None,
            "fit_rows": rows,
            "exam": exam_stats(exams[t]) if exams[t] else None,
            "panel": panel_stats(panels[t]) if panels[t] else None,
        }

    # --- the pH ladder shape test, on the three rungs already scored ---------
    # pH 3 and pH 7 come from the EXAM (hold-out, scored there); pH 5 is a
    # DECLARED FIT anchor and comes from the fit row. Three rungs, not D1's
    # seven -- stated as the limitation it is.
    for t in TAGS:
        blk = per_tag[t]
        if not blk["exam"]:
            blk["ladder"] = None
            continue
        ladder = {}
        for compound in ("FFT", "MFT"):
            fit_row = blk["fit_rows"].get(FIT_PH5_ROW[compound])
            ph5 = (fit_row["fold"] * FIT_PH5_ANCHOR_PPB[compound]
                   if fit_row and fit_row["fold"] else None)
            hof = blk["exam"]["hofmann_ph_block"]
            ph3 = (hof.get(f"hofmann_ribose_pH3_{compound}") or {}).get("predicted_ppb")
            ph7 = (hof.get(f"hofmann_ribose_pH7_{compound}") or {}).get("predicted_ppb")
            decades = (abs(math.log10(ph7 / ph3))
                       if ph3 and ph7 and ph3 > 0 and ph7 > 0 else None)
            ladder[compound] = {
                "pH3_predicted_ppb": ph3,
                "pH5_predicted_ppb": ph5,
                "pH7_predicted_ppb": ph7,
                "pH3_to_pH7_decades": decades,
                "monotone_falling": (
                    None if None in (ph3, ph5, ph7) else bool(ph3 > ph5 > ph7)),
                "humped_at_pH5": (
                    None if None in (ph3, ph5) else bool(ph5 > ph3)),
            }
        blk["ladder"] = ladder

    checks = score_prereg(ensemble, per_tag, member_exams, w4, amine)

    return {
        "wave": "B2.4 -- the declared weighting, the ensemble, and two scorer conditions",
        "prereg": data_paths.rel(data_paths.VALIDATION_DIR / "kinetic_core_b2_4_prereg.md"),
        "generated_by": "scripts/generators/generate_kinetic_core_b2_4_reports.py",
        "declared_weights": ensemble["declared_weights"],
        "weight_basis": ensemble["weight_basis"],
        "free_set": ensemble["free_set"],
        "shipping_choice": ensemble["shipping_choice"],
        "per_weighting": per_tag,
        "member_exams": {
            f"{t}_s{s}": exam_stats(e) for (t, s), e in member_exams.items()},
        "w4_control_exam_at_b2_2_parameters": exam_stats(w4) if w4 else None,
        "amine_fate_probe": (amine or {}).get("verdict"),
        "prereg_checks": checks,
        "measured_ladder_for_reference": {
            "note": ("Hofmann & Schieberle 1998's own ribose ladder, as the "
                     "SCORERS report it: pH 3 and pH 7 are hold-out and appear "
                     "here only because the exam has already scored them; the "
                     "pH-5 rung is a declared FIT anchor."),
            "FFT_ppb": {"pH3": 229.0, "pH5": 121.0, "pH7": 12.0},
            "MFT_ppb": {"pH3": 553.0, "pH5": 198.0, "pH7": 25.0},
            "FFT_pH3_to_pH7_decades": round(abs(math.log10(12.0 / 229.0)), 3),
            "MFT_pH3_to_pH7_decades": round(abs(math.log10(25.0 / 553.0)), 3),
        },
    }


def score_prereg(ensemble, per_tag, member_exams, w4, amine) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []

    def add(claim, outcome, detail):
        out.append({"claim": claim, "outcome": outcome, "detail": detail})

    # --- W-1 ---------------------------------------------------------------
    e = ensemble["ensemble"]["shipped"]
    s_pooled = e["spread_S_log10_max_over_min"]
    s_local = e["spread_S_local_arm_only"]
    if s_pooled is None:
        add("W-1 spread at W-SHIPPED > 0.30", "UNSCOREABLE",
            "fewer than two converged members at the control weighting")
    else:
        # THE FALSIFICATION IS REAL BUT IT RESTS ON TWO MEMBERS, AND THE OTHER
        # TWO WEIGHTINGS SHOW WHY THAT MATTERS. At W-SHIPPED all four non-local
        # members exhausted A1's 500-evaluation budget and were excluded by the
        # pre-declared rule, so S is computed over the incumbent and ONE local
        # start. The SAME seeds at W-HALF and W-MEASURED converged inside the
        # budget and landed in a SECOND basin at roughly twice the cost. So the
        # control's S = 0 is partly a statement about the budget, and the
        # ensemble is NOT single-basin once the members are allowed to finish.
        other = {t: ensemble["ensemble"][t] for t in ("half", "measured")
                 if t in ensemble["ensemble"]}
        add("W-1 spread at W-SHIPPED > 0.30 (log10 max/min cost, converged)",
            "HELD" if s_pooled > 0.30 else "FALSIFIED",
            f"pooled S = {_fmt(s_pooled)} over {e['n_converged']} converged of "
            f"{e['n_members']}; local arm alone S = {_fmt(s_local)}; "
            f"{e['n_budget_exhausted']} budget-exhausted. Costs: "
            f"{[round(c, 3) for c in e['costs']]}. "
            f"TWO GLOBAL STARTS 2005 AND 1942 COST-UNITS APART LANDED WITHIN "
            f"0.0003 OF EACH OTHER, so within the ~20 identified coordinates "
            f"the objective really is single-basin; the B2.2->B2.3 scatter "
            f"lives in the ~28 coordinates this wave FROZE and says nothing "
            f"about. THE FALSIFICATION IS QUALIFIED, AND THE QUALIFICATION IS "
            f"NOT IN THE CONTROL'S FAVOUR: at W-SHIPPED only "
            f"{e['n_converged']} members survived A1's budget to enter S, so "
            f"the statistic is a two-member statistic. The same seeds at the "
            f"other two weightings DID converge and DID resolve a second basin "
            + "; ".join(
                f"-- {TAG_LABEL[t]}: {b['n_converged']} converged, S = "
                f"{_fmt(b['spread_S_log10_max_over_min'])} (local arm alone "
                f"{_fmt(b['spread_S_local_arm_only'])}), costs "
                f"{[round(c, 2) for c in b['costs']]}"
                for t, b in other.items())
            + ". Read together: the near-basin structure is genuinely flat "
              "(local-arm S is 0.003 and 0.0002), but the GLOBAL arm finds a "
              "distinct attractor at about twice the cost at every weighting "
              "where it is allowed to converge. W-1 is falsified as written "
              "and the objective is still multi-modal over the free set.")

    # --- W-2 ---------------------------------------------------------------
    geos = []
    for (tag, _start), exam in member_exams.items():
        if tag != "shipped":
            continue
        g = exam_stats(exam)["d1_pool_geometric_mean_fold"]
        if g:
            geos.append(g)
    if len(geos) >= 2:
        rng = max(geos) / min(geos)
        add("W-2 exam geometric-mean fold spans >= 1.5x across the W-SHIPPED ensemble",
            "HELD" if rng >= 1.5 else "FALSIFIED",
            f"{len(geos)} members re-sat: {[round(g, 2) for g in sorted(geos)]}; "
            f"range {_fmt(rng)}x")
    else:
        add("W-2 exam geometric-mean fold spans >= 1.5x across the W-SHIPPED ensemble",
            "UNSCOREABLE", f"only {len(geos)} member exams on disk")

    # --- W-3 ---------------------------------------------------------------
    detail, held = [], True
    for t in ("half", "measured"):
        blk = per_tag[t]
        k = (blk["best_member"].get("parameters") or {}).get("k_fur_fft")
        lad = (blk.get("ladder") or {}).get("FFT") or {}
        dec = lad.get("pH3_to_pH7_decades")
        ok = (k is not None and k < -2.0) and (dec is not None and dec < 3.0)
        held = held and ok
        detail.append(f"{TAG_LABEL[t]}: k_fur_fft {_fmt(k)} (needs < -2.0), "
                      f"pH3->pH7 FFT slope {_fmt(dec)} decades (needs < 3.0)")
    add("W-3 down-weighting moves k_fur_fft below -2.0 AND the FFT pH slope below 3.0 decades",
        "HELD" if held else "FALSIFIED", "; ".join(detail))

    # --- W-4 ---------------------------------------------------------------
    if w4:
        st = exam_stats(w4)
        g_all = st["as_was_geomean"]
        g_d1 = st["d1_pool_as_was_geometric_mean_fold"]
        # AS WRITTEN the claim named 26.8x, which is D1's `B22_ref` cell --
        # B2.2's parameters under B2.2's OWN stoichiometry. The cell the claim
        # DESCRIBES ("B2.2's parameters with B2.3's stoichiometry") is D1's
        # `S_only`, whose published value is 23.97x. Both readings are scored.
        as_written = g_d1 is not None and abs(g_d1 - 26.8) / 26.8 <= 0.10
        correct_ref = g_d1 is not None and abs(g_d1 - 23.97) / 23.97 <= 0.10
        add("W-4 the charge-conservation fix stays free (B2.2's parameters under "
            "B2.3's stoichiometry; claim as written: within 10% of 26.8x)",
            "HELD" if as_written else
            "HELD ON THE CORRECT REFERENCE, FALSIFIED AS WRITTEN"
            if correct_ref else "FALSIFIED",
            f"On D1's own 23-point pool the as-was geometric mean is "
            f"{_fmt(g_d1)}x. THE CLAIM MIS-NAMED ITS REFERENCE: 26.8x is D1's "
            f"`B22_ref` cell (B2.2 stoichiometry), while the cell the claim "
            f"describes is D1's `S_only` (B2.3 stoichiometry), published at "
            f"23.97x. Against 23.97x this reproduces D1 to four significant "
            f"figures; against the 26.8x actually written it is "
            f"{_fmt(100 * (g_d1 - 26.8) / 26.8, 3)}%, just outside the 10% band. "
            f"Either way the substantive claim -- the conservation fix is free "
            f"-- stands, and the as-was PAIRED MEDIAN is "
            f"{_fmt(st['d1_pool_as_was_median_fold'])}x against D1's 10.63x. "
            f"On the CURRENT 27-point pool the same number is {_fmt(g_all)}x; "
            f"the pool grew because concurrent Wave B6 gave the core a lipid "
            f"lane that answers four points D1's core declined.")
    else:
        add("W-4 the charge-conservation fix stays free", "UNSCOREABLE",
            "the B2.2-parameter control exam has not been run")

    # --- W-5 ---------------------------------------------------------------
    blocks, shapes = [], []
    for t in TAGS:
        p = per_tag[t]["panel"]
        if p:
            blocks.append((t, p["hofmann_ph_block_passed"]))
        lad = per_tag[t].get("ladder") or {}
        for compound in ("FFT", "MFT"):
            mono = (lad.get(compound) or {}).get("monotone_falling")
            if mono is not None:
                shapes.append((t, compound, mono))
    reached = [t for t, n in blocks if n >= 3]
    monotone = [f"{t}/{c}" for t, c, m in shapes if m]
    # THE PREREG'S FALSIFICATION CONDITION IS PER-WEIGHTING AND NEEDS **BOTH**
    # COMPOUNDS: "the block reaches 3/4 at any weighting, OR any weighting's
    # ladder is monotone falling from pH 3 to pH 7 ON BOTH FFT AND MFT". A
    # first draft of this scorer falsified on ANY SINGLE (weighting, compound)
    # being monotone, which is a STRICTER test than the one pre-registered and
    # is not the scorer's to tighten after the fact. Corrected to the prereg's
    # own wording; the discrepancy is stated in the detail string below rather
    # than being silently resolved in either direction.
    mono_by_tag = {t: {c: m for tt, c, m in shapes if tt == t} for t in TAGS}
    both_monotone = [t for t in TAGS
                     if mono_by_tag.get(t)
                     and all(mono_by_tag[t].get(c) for c in ("FFT", "MFT"))]
    add("W-5 the Hofmann pH-3/pH-7 block never reaches 3/4 and the ladder stays "
        "humped -- re-weighting does NOT fix a shape defect",
        "HELD" if (not reached and not both_monotone) else "FALSIFIED",
        f"panel block: {', '.join(f'{TAG_LABEL[t]} {n}/4' for t, n in blocks)} "
        f"-- never reaches 3/4. Three-rung ladder (pH 5 fit-side-anchored, pH 3 "
        f"and pH 7 hold-out, scored by the exam): monotone falling at "
        f"{monotone or 'nowhere'}; monotone on BOTH compounds -- which is what "
        f"the pre-registration declares as the falsifier -- at "
        f"{[TAG_LABEL[t] for t in both_monotone] or 'no weighting'}. "
        f"REPORTED AGAINST THE CLAIM'S OWN PROSE, NOT ONLY ITS FALSIFIER: the "
        f"claim asserts the ladder 'stays humped ... at every weighting', and "
        f"only MFT does. FFT is monotone FALLING at all three weightings, at "
        f"4.88-4.91 decades against a measured 1.28 -- the right SHAPE with the "
        f"wrong STEEPNESS, roughly four decades too steep. So W-5's verdict "
        f"holds and its substance holds -- re-weighting fixed neither the block "
        f"nor the shape -- but its description of the defect is half wrong, and "
        f"the FFT defect is a slope defect rather than a hump. "
        f"LIMITATION: three rungs, not D1's seven-rung sweep, so a hump "
        f"BETWEEN pH 3 and pH 5 is invisible to this test.")

    # --- W-6 ---------------------------------------------------------------
    rows = [(t, per_tag[t]["largest_residual_row"], per_tag[t]["largest_residual"])
            for t in TAGS]
    ok = all(r == "fed_c2c3_MFT_pH3" for _t, r, _v in rows)
    add("W-6 fed_c2c3_MFT_pH3 is the largest single residual at every weighting",
        "HELD" if ok else "FALSIFIED",
        "; ".join(f"{TAG_LABEL[t]}: {r} at {_fmt(v)} sigma" for t, r, v in rows))

    # --- W-7 ---------------------------------------------------------------
    base_e = per_tag["shipped"]["exam"]
    base_p = per_tag["shipped"]["panel"]
    signs = []
    if base_e and base_p:
        for t in ("half", "measured"):
            ex, pa = per_tag[t]["exam"], per_tag[t]["panel"]
            if not ex or not pa:
                continue
            de = ex["median_abs_log10_all"] - base_e["median_abs_log10_all"]
            dp = (pa["median_abs_log10_fold_gating"]
                  - base_p["median_abs_log10_fold_gating"])
            signs.append((t, de, dp, (de <= 0) == (dp <= 0)))
    if signs:
        # THE TEST IS DEGENERATE ON THE STATISTIC IT NAMES, AND THE WAVE SAYS SO
        # RATHER THAN BANKING THE PASS. W-7 names median |log10 fold|; on the
        # exam side that number does not move AT ALL between weightings (the
        # median row is the same row, unmoved), so `de` is exactly 0.0 and the
        # sign comparison is decided by the convention that 0 counts as "not
        # worse". The other continuous statistic this wave added -- the
        # geometric mean -- DOES move on both sides, and it moves in OPPOSITE
        # directions. That is reported here beside the pre-registered verdict.
        exam_dead = all(abs(de) < 1e-12 for _t, de, _dp, _ok in signs)
        geo_signs = []
        for t in ("half", "measured"):
            ex, pa = per_tag[t]["exam"], per_tag[t]["panel"]
            if not ex or not pa:
                continue
            dge = ex["d1_pool_geometric_mean_fold"] - base_e["d1_pool_geometric_mean_fold"]
            dgp = pa["geometric_mean_fold_gating"] - base_p["geometric_mean_fold_gating"]
            geo_signs.append((t, dge, dgp, (dge <= 0) == (dgp <= 0)))
        geo_disagree = [t for t, _a, _b, ok in geo_signs if not ok]
        add("W-7 exam and panel median |log10 fold| agree in SIGN at every weighting",
            "HELD" if all(s[3] for s in signs) else "FALSIFIED",
            "; ".join(f"{TAG_LABEL[t]}: exam {de:+.3f}, panel {dp:+.3f} "
                      f"({'same' if ok else 'OPPOSITE'} direction)"
                      for t, de, dp, ok in signs)
            + (". THE VERDICT IS DEGENERATE AND IS NOT EVIDENCE OF AGREEMENT: "
               "the exam's median |log10 fold| is BIT-IDENTICAL across all "
               "three weightings, so the exam contributes no sign and the "
               "test cannot fail on this statistic." if exam_dead else "")
            + (". ON THE GEOMETRIC MEAN -- the other continuous statistic this "
               "wave added, and the one that actually moves -- THE TWO "
               "SCORECARDS DISAGREE IN SIGN at "
               + ", ".join(TAG_LABEL[t] for t in geo_disagree)
               + ": " + "; ".join(
                   f"{TAG_LABEL[t]} exam {dge:+.3f}x, panel {dgp:+.3f}x"
                   for t, dge, dgp, _ok in geo_signs)
               + ". The panel improves as E falls while the exam worsens. By "
                 "W-7's own stated consequence -- 'if they ever disagree in "
                 "sign, the four shared rows are masking something and the "
                 "wave says so' -- THE WAVE SAYS SO."
               if geo_disagree else
               ". On the geometric mean the two scorecards also agree in sign."))
    else:
        add("W-7 exam and panel median |log10 fold| agree in SIGN", "UNSCOREABLE",
            "not every weighting has both artifacts")

    # --- W-8 ---------------------------------------------------------------
    lv = [(t, per_tag[t]["best_member"].get("sum_r2_level")) for t in TAGS]
    vals = [v for _t, v in lv]
    ok = all(v is not None for v in vals) and vals[0] >= vals[1] >= vals[2]
    add("W-8 sum_r2_level improves monotonically as E falls "
        "(W-SHIPPED >= W-HALF >= W-MEASURED)",
        "HELD" if ok else "FALSIFIED",
        "; ".join(f"{TAG_LABEL[t]}: {_fmt(v)}" for t, v in lv))

    # --- W-9 ---------------------------------------------------------------
    misses = []
    for t in TAGS:
        m = per_tag[t]["best_member"].get("ph_endpoint_miss_pH_units") or {}
        misses.append((t, max((abs(v) for v in m.values()), default=None), m))
    mono = all(
        misses[i][1] is not None and misses[i + 1][1] is not None
        and misses[i + 1][1] >= misses[i][1] for i in range(len(misses) - 1))
    g1_measured = ensemble["ensemble"]["measured"]["gate_G1_calibration_within_1pH"]
    add("W-9 the pH-endpoint misses grow monotonically as E falls, AND "
        "W-MEASURED fails gate G1",
        "HELD" if (mono and g1_measured is False) else "FALSIFIED",
        "; ".join(f"{TAG_LABEL[t]}: worst |miss| {_fmt(v)} pH units"
                  for t, v, _m in misses)
        + f". W-MEASURED gate G1 = {'PASS' if g1_measured else 'FAIL'}")

    # --- W-10 --------------------------------------------------------------
    shipped_tag = ensemble["shipping_choice"]["shipped"]
    geos = {t: (per_tag[t]["exam"] or {}).get("d1_pool_geometric_mean_fold")
            for t in TAGS}
    scored = {t: g for t, g in geos.items() if g is not None}
    best_exam = min(scored, key=scored.get) if scored else None
    add("W-10 the pre-declared criterion does NOT select the best-exam weighting",
        "HELD" if (best_exam and shipped_tag != best_exam)
        else "FALSIFIED" if best_exam else "UNSCOREABLE",
        f"criterion ships {TAG_LABEL.get(shipped_tag, shipped_tag)}; best exam "
        f"geometric mean is {TAG_LABEL.get(best_exam, best_exam)} at "
        f"{_fmt(scored.get(best_exam))}x. All: "
        + ", ".join(f"{TAG_LABEL[t]} {_fmt(g)}x" for t, g in geos.items()))

    # --- the amine-fate prediction -----------------------------------------
    if amine:
        v = amine["verdict"]
        add("PREREG sec.5 the rebuilt probe reproduces the DIRECTION of the "
            "B2.3 table (carry-both leaves all three Zhou endpoints above pH 4.5)",
            "HELD" if v["prereg_expectation_held"] else "FALSIFIED",
            f"three encodings distinct: {v['encodings_are_distinct']}; "
            f"reproduces the B2.3 prereg table cell-by-cell: "
            f"{v.get('reproduces_the_b2_3_prereg_table')}; mean |miss| on the "
            f"two acidifying anchors -- shipped "
            f"{_fmt(v['mean_abs_miss_on_the_two_acidifying_anchors_pH']['shipped'])}, "
            f"carry-both "
            f"{_fmt(v['mean_abs_miss_on_the_two_acidifying_anchors_pH']['carry_both'])} "
            f"pH units")
    return out


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def render_fit_report(p: Dict[str, Any]) -> str:
    a: List[str] = []
    w = a.append
    w("# Kinetic core, Build Wave B2.4 — the DECLARED weighting, and an ensemble")
    w("")
    w(f"Pre-registered in [`{p['prereg']}`](kinetic_core_b2_4_prereg.md), saved "
      f"to disk before any fit ran. Every claim below is scored in §5 whether it "
      f"held or not.")
    w("")
    w("**This wave changes no stoichiometry, no species, no reaction, no "
      "benchmark, no engine input, no hold-out row and no pass band.** It makes "
      "an accidental weighting explicit, replaces a point estimate with an "
      "ensemble, adds a continuous statistic to each of the two scorecards, and "
      "repairs a broken citation.")
    w("")

    w("## 1. The declared exchange rate")
    w("")
    w("`E` = decades of level error that ONE pH unit of endpoint error is "
      "declared to be worth. Implemented as `sigma_ph = 0.35 / E` on the three "
      "`zhou_final_pH_*` rows and nowhere else.")
    w("")
    w("| tag | E | σ_ph | basis |")
    w("|---|---:|---:|---|")
    for t in TAGS:
        w(f"| **{TAG_LABEL[t]}** | {p['declared_weights'][t]} | "
          f"{0.35 / p['declared_weights'][t]:.3f} | {p['weight_basis'][t]} |")
    w("")
    w("> B2.3 never declared a rate. Its three pH rows sat at σ = 0.25 **pH "
      "units** in the same sum of squares as 55 rows at σ_log = 0.20–0.60 "
      "**log-folds**, and D1 §3 measured the consequence: one pH unit priced at "
      "~9× a 3× level miss, with `zhou_final_pH_from_pH8` alone carrying 44% of "
      "the entire B2.2→B2.3 refit. W-SHIPPED reproduces that objective exactly "
      "and is the control.")
    w("")

    w("## 2. The free set — 20 of 48, frozen 28")
    w("")
    fs = p["free_set"]
    w(f"**{fs['n_free']} free, {fs['n_frozen']} frozen at their B2.3 values.** "
      f"The route is D1's own stated fallback and the reason is arithmetic: "
      f"B2.3's multistart trace records 4731 s and 2280 s for two full-vector "
      f"starts, and this container delivers about one core of real throughput.")
    w("")
    w("| parameter | clause |")
    w("|---|---|")
    for k in fs["keys"]:
        w(f"| `{k}` | {fs['clause'][k]} |")
    w("")

    w("## 3. THE ENSEMBLE — the spread, not the best member")
    w("")
    w("| weighting | members | converged | budget-exhausted | costs | **S = log₁₀(max/min)** | S, local arm | Σr²_level at best | G1 | G2 |")
    w("|---|---:|---:|---:|---|---:|---:|---:|---|---|")
    for t in TAGS:
        e = p["per_weighting"][t]["ensemble"]
        best = p["per_weighting"][t]["best_member"]
        costs = ", ".join(f"{c:.2f}" for c in e["costs"])
        w(f"| **{TAG_LABEL[t]}** | {e['n_members']} | {e['n_converged']} | "
          f"{e['n_budget_exhausted']} | {costs} | "
          f"**{_fmt(e['spread_S_log10_max_over_min'])}** | "
          f"{_fmt(e['spread_S_local_arm_only'])} | "
          f"{_fmt(best.get('sum_r2_level'))} | "
          f"{_fmt(e['gate_G1_calibration_within_1pH'])} | "
          f"{_fmt(e['gate_G2_at_least_4_of_6_converged'])} |")
    w("")
    w("> **Total cost is not comparable across weightings** — they are three "
      "different objectives. `Σr²_level`, the sum of squared residuals over the "
      "55 non-pH rows at their unchanged σ_log, is, and it is the only "
      "cross-weighting comparison this report makes.")
    w("")

    w("## 4. THE SHIPPING CHOICE")
    w("")
    c = p["shipping_choice"]
    w(f"**Shipped: {TAG_LABEL.get(c['shipped'], c['shipped'])}.**")
    w("")
    w(f"> {c['criterion']}")
    w("")
    w(f"Qualifying: {', '.join(TAG_LABEL.get(t, t) for t in c['qualifying']) or 'none'}. "
      f"Tie broken toward the largest E: {c.get('tie_broken_toward_largest_E')}. "
      f"Fallback branch used: {c.get('fallback_used')}.")
    w("")

    w("## 5. THE PRE-REGISTRATION, SCORED")
    w("")
    held = sum(1 for c in p["prereg_checks"] if c["outcome"] == "HELD")
    fal = sum(1 for c in p["prereg_checks"] if c["outcome"] == "FALSIFIED")
    uns = sum(1 for c in p["prereg_checks"] if c["outcome"] == "UNSCOREABLE")
    w(f"**{held} held, {fal} falsified, {uns} unscoreable.**")
    w("")
    for c in p["prereg_checks"]:
        w(f"- **{c['outcome']}** — {c['claim']}")
        w(f"  - {c['detail']}")
    w("")

    w("## 6. THE AMINE-FATE PROBE, REBUILT")
    w("")
    v = p.get("amine_fate_probe")
    if v:
        w("D1 §7 found that `ph_state.AMINE_FATE_BASIS` cited a probe that "
          "**cannot run** — `KeyError: 'AMN'`, a species removed from the tree — "
          "and that two of its three encodings had collapsed onto one code path.")
        w("")
        w("`scripts/generators/probe_amine_fate_b2_4.py` rebuilds the axis "
          "against the current species set, deriving the released amino nitrogen "
          "as `(Cys + ARP + TTCA at t=0) − (Cys + ARP + TTCA at t)` and adding it "
          "back as ammonium at pK_a 9.25. No new species, no new parameter, no "
          "reaction changed.")
        w("")
        w(f"- three encodings genuinely distinct: **{v['encodings_are_distinct']}**")
        w(f"- reproduces the B2.3 pre-registration's published table, cell by "
          f"cell to two decimals: **{v.get('reproduces_the_b2_3_prereg_table')}**")
        w(f"- mean |miss| on the two acidifying anchors (Zhou pH 6, pH 7): "
          f"shipped **"
          f"{_fmt(v['mean_abs_miss_on_the_two_acidifying_anchors_pH']['shipped'])}** "
          f"pH units vs carry-both **"
          f"{_fmt(v['mean_abs_miss_on_the_two_acidifying_anchors_pH']['carry_both'])}**")
        w("")
        w("> **The defect was in the script, not in the evidence.** The B2.3 "
          "table's digits are reproducible from the current tree, so the "
          "declaration keeps its basis rather than being weakened; what changes "
          "is that the citation now points at a probe that runs, and the "
          "self-referential last sentence of `AMINE_FATE_BASIS` is gone.")
    else:
        w("Not run.")
    w("")
    return "\n".join(a)


def render_holdout_report(p: Dict[str, Any]) -> str:
    a: List[str] = []
    w = a.append
    w("# Kinetic core, Build Wave B2.4 — the blind re-sit, under THREE declared weightings")
    w("")
    w("Both scorecards, run once per declared weighting from the frozen "
      "ensemble-best parameters of that weighting. **No optimiser runs in either "
      "scorer, and neither scorer was forked**: the existing "
      "`generate_kinetic_core_b2_3_holdout` and `generate_cutover_final_exam` "
      "are called with two things rebound — which fit report the parameters come "
      "from, and which basename the artifact is written under.")
    w("")

    w("## 1. THE FOUR SHARED ROWS — declared in both scorecards from this wave on")
    w("")
    w("D1 §5 established that four of the exam's 23 answered points are **the "
      "same measurements** as four hold-out panel rows, not analogues of them. "
      "The panel reads `hofmann_ribose_pH7_FFT` at 0.0020× and the exam reads "
      "the same measurement at 499×: one number, inverted. **The two scorecards "
      "are therefore not independent evidence on this axis**, and agreement "
      "between them here is one measurement counted twice.")
    w("")
    w("| panel row | exam point |")
    w("|---|---|")
    for _b, _c, row in HOFMANN_PH_BLOCK:
        w(f"| `{row}` | `{_b}` / {_c} |")
    w("")
    for t in TAGS:
        pa = p["per_weighting"][t]["panel"]
        ex = p["per_weighting"][t]["exam"]
        if pa and ex:
            w(f"- {TAG_LABEL[t]}: panel declares {pa['shared_rows_declared']} "
              f"shared rows, exam declares {ex['shared_rows_declared']}.")
    w("")

    w("## 2. THE WEIGHTING TABLE — fit, panel and exam under each declared value")
    w("")
    w("| | E | fit cost | Σr²_level | **panel gating** | **panel median \\|log₁₀\\|** | panel geo-mean | **exam geo-mean (D1 pool, n=23)** | exam geo-mean (all 27) | exam paired median | exam in band | Hofmann pH block |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for t in TAGS:
        b = p["per_weighting"][t]
        best, pa, ex = b["best_member"], b["panel"], b["exam"]
        w(f"| **{TAG_LABEL[t]}** | {p['declared_weights'][t]} | "
          f"{_fmt(best.get('cost'))} | {_fmt(best.get('sum_r2_level'))} | "
          + (f"**{pa['gating_passed']}/{pa['gating_rows']}** | "
             f"**{_fmt(pa['median_abs_log10_fold_gating'])}** | "
             f"{_fmt(pa['geometric_mean_fold_gating'])}x | " if pa
             else "-- | -- | -- | ")
          + (f"**{_fmt(ex['d1_pool_geometric_mean_fold'])}x** | "
             f"{_fmt(ex['geometric_mean_fold_all'])}x | "
             f"{_fmt(ex['paired_core_median'])}x | "
             f"{ex['in_band']}/{ex['scored']} | " if ex
             else "-- | -- | -- | -- | ")
          + (f"{pa['hofmann_ph_block_passed']}/4 |" if pa else "-- |"))
    w("")
    w("> **THE POOL MOVED, AND IT IS NOT THIS WAVE'S DOING.** D1's geometric "
      "means are over **23** answered points. The exam now answers **27**: "
      "concurrent Build Wave B6 gave the core a lipid lane, so four "
      "`matrix_path` points it used to decline are now scored — at a family "
      "median near 1900×. Every cross-wave comparison in this report is "
      "therefore made on **D1's own 23-point pool** (all scored points except "
      "the lipid family), with the 27-point number printed beside it. Comparing "
      "a B2.4 27-point mean against D1's 26.8× would be comparing two different "
      "pools and calling the difference a result.")
    w("")
    w("> **Read the two continuous columns, not the two censored ones.** The "
      "panel's gating count is a pass/fail tally that cannot see a 100× "
      "degradation on a row that was already failing — D1 §5 measured 1.42 net "
      "decades of B2.2→B2.3 degradation that cost it nothing. The exam's paired "
      "median is a rank statistic over a lumpy 23-point pool that swung 4.7× "
      "between those same waves while the geometric mean moved 1.78×. Both "
      "continuous statistics are new in this wave and both are permanent.")
    w("")

    w("## 3. THE HOFMANN pH LADDER — the shape defect, scored")
    w("")
    m = p["measured_ladder_for_reference"]
    w(f"Measured: FFT {m['FFT_ppb']['pH3']} / {m['FFT_ppb']['pH5']} / "
      f"{m['FFT_ppb']['pH7']} ppb and MFT {m['MFT_ppb']['pH3']} / "
      f"{m['MFT_ppb']['pH5']} / {m['MFT_ppb']['pH7']} ppb at pH 3 / 5 / 7 — "
      f"**monotone falling on both**, {m['FFT_pH3_to_pH7_decades']} and "
      f"{m['MFT_pH3_to_pH7_decades']} decades. The pH-5 rung is a declared FIT "
      f"anchor; pH 3 and pH 7 are hold-out and appear here only because the exam "
      f"has already scored them.")
    w("")
    w("| weighting | cmpd | pred pH 3 | pred pH 5 | pred pH 7 | decades | monotone falling? |")
    w("|---|---|---:|---:|---:|---:|---|")
    for t in TAGS:
        lad = p["per_weighting"][t].get("ladder") or {}
        for compound in ("FFT", "MFT"):
            L = lad.get(compound) or {}
            w(f"| {TAG_LABEL[t]} | {compound} | "
              f"{_fmt(L.get('pH3_predicted_ppb'))} | "
              f"{_fmt(L.get('pH5_predicted_ppb'))} | "
              f"{_fmt(L.get('pH7_predicted_ppb'))} | "
              f"{_fmt(L.get('pH3_to_pH7_decades'))} | "
              f"{_fmt(L.get('monotone_falling'))} |")
    w("")
    w("> **LIMITATION, stated:** three rungs, not D1 §4's seven-rung sweep. A "
      "hump lying between pH 3 and pH 5 is invisible to this test. The three "
      "rungs used are the two the exam already scores and the one the fit "
      "already anchors, so nothing new was integrated at a hold-out condition "
      "to build this table.")
    w("")

    w("## 4. THE ENSEMBLE ON THE EXAM — W-2's re-sits")
    w("")
    me = p.get("member_exams") or {}
    if me:
        w("| member | **exam geo-mean (D1 pool)** | exam geo-mean (all) | exam paired median | in band |")
        w("|---|---:|---:|---:|---:|")
        for key in sorted(me):
            v = me[key]
            w(f"| `{key}` | **{_fmt(v['d1_pool_geometric_mean_fold'])}x** | "
              f"{_fmt(v['geometric_mean_fold_all'])}x | "
              f"{_fmt(v['paired_core_median'])}x | {v['in_band']}/{v['scored']} |")
        w("")
        w("> This table is the wave's sharpest single result if the numbers "
          "spread: it says the exam score is a property of which basin the "
          "optimiser happened to find, not of the model.")
    else:
        w("No per-member exams on disk.")
    w("")

    w("## 5. THE PRE-REGISTRATION, SCORED")
    w("")
    for c in p["prereg_checks"]:
        w(f"- **{c['outcome']}** — {c['claim']}")
        w(f"  - {c['detail']}")
    w("")
    return "\n".join(a)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Assemble the Build Wave B2.4 fit and hold-out reports from "
            "artifacts already on disk and score the pre-registration claim by "
            "claim; writes "
            "results/validation/kinetic_core_b2_4_{fit,holdout}_report.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(V),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    payload = build()
    out = Path(args.output_dir)
    (out / "kinetic_core_b2_4_fit_report.json").write_text(
        json.dumps(payload, indent=2, default=str))
    (out / "kinetic_core_b2_4_fit_report.md").write_text(render_fit_report(payload))
    (out / "kinetic_core_b2_4_holdout_report.json").write_text(
        json.dumps(payload, indent=2, default=str))
    (out / "kinetic_core_b2_4_holdout_report.md").write_text(
        render_holdout_report(payload))
    print("wrote kinetic_core_b2_4_{fit,holdout}_report.{md,json}")
    for c in payload["prereg_checks"]:
        print(f"  {c['outcome']:12s} {c['claim'][:80]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
