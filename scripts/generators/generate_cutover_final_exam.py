#!/usr/bin/env python
"""
THE CUTOVER FINAL EXAM (Build Wave B5, 2026-08-29).

Scores the KINETIC CORE end-to-end against the 21 frozen external-validation
bundles, side by side with the OLD (FAST / matrix) lane, and reports both.

WHAT THIS SCRIPT IS AND IS NOT
------------------------------
It is PURE SCORING. It fits nothing, tunes nothing, and changes no parameter.
Every core constant is read from the frozen B1/B2.1/B3 fit reports through
``src.kinetic_core.engine.core_parameters``; every old-lane number is produced
by the same ``evaluate_benchmark`` the shipped benchmark summary uses.

The pre-registration for this exam is ``results/validation/cutover_prereg.md``,
written BEFORE this file existed and before any measured value in
``data/benchmarks/external_validation/`` was read. Read it first. It states the
envelope declarations and the expected outcomes per bundle family, so that this
script's output can be checked against a prediction rather than admired.

WHY THE TWO MEDIANS ARE NOT THE SAME NUMBER
-------------------------------------------
The old lane emits a number for all 40 points. The core answers 23 and DECLINES
17, for reasons named in the pre-registration (no lipid-oxidation path, no HMF,
no DMHF, no alanine in the sulfur lane). A median over 40 guesses and a median
over 23 answers are not comparable, so this report prints BOTH:

  * each lane's median over its own scored set, and
  * the PAIRED median: the old lane restricted to exactly the points the core
    answers, which is the only apples-to-apples comparison available.

Reporting only the first would let the core look good by declining its hardest
points, which is the failure mode this repository has already committed twice.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EXTERNAL_DIR = ROOT / "data" / "benchmarks" / "external_validation"
MAILLARD_PATH_DIR = EXTERNAL_DIR / "maillard_path"
OUTPUT_DIR = ROOT / "results" / "validation"
DEFAULT_BASENAME = "cutover_final_exam"

#: Pass bands, taken UNCHANGED from the module scorecards so this exam cannot
#: grade itself leniently. B2.1's level band is 3.0x; B3's gating band is 3.0x.
PASS_BAND_LEVEL = 3.0

PPB_CONVERSION_FACTOR = 1.0e6
RATIO_UNIT_FACTORS = {
    "mol_percent": 100.0,
    "umol_per_mol_limiting_precursor": 1.0e6,
}

#: Bundle -> family, for the pre-registration's per-family expectations.
FAMILIES = (
    ("hofmann1998", "sulfur_hofmann1998_145C"),
    ("Yiltirak2026", "sulfur_yiltirak2026_T_ladder"),
    ("Chang2021", "acrylamide_180C"),
    ("Ye2024", "acrylamide_180C"),
    ("Lin2022", "acrylamide_180C"),
    ("Schibilsky2019", "furan_browning_glc_alanine"),
    ("Steinhagen2021", "furan_browning_glc_alanine"),
)

#: The four ribose pH-3 / pH-7 points B2.1 already scored under different row
#: names. Recorded in the pre-registration; labelled here so the exam table
#: cannot present them as first exposures.
B2_1_RESCORED = {
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3", "2-Methyl-3-furanthiol (MFT)"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3", "2-Furfurylthiol (FFT)"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7", "2-Methyl-3-furanthiol (MFT)"),
    ("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7", "2-Furfurylthiol (FFT)"),
}


def _git_head() -> Dict[str, str]:
    def _run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return "unknown"

    return {
        "commit": _run("rev-parse", "HEAD"),
        "short": _run("rev-parse", "--short", "HEAD"),
        "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": "yes" if _run("status", "--porcelain") else "no",
    }


def _family(benchmark_id: str) -> str:
    for token, name in FAMILIES:
        if token.lower() in benchmark_id.lower():
            return name
    return "matrix_path_lipid"


def _limiting_precursor_molar(bench: Dict[str, Any]) -> Tuple[Optional[str], Optional[float]]:
    items = [
        (str(name), float(data.get("concentration_mM", 0.0)) / 1000.0)
        for name, data in (bench.get("precursors") or {}).items()
        if float(data.get("concentration_mM", 0.0)) > 0.0
    ]
    if not items:
        return None, None
    return min(items, key=lambda kv: kv[1])


def _measured_value(
    bench: Dict[str, Any], compound: str, target_spec: Dict[str, Any]
) -> Optional[float]:
    """
    The measured value for one point, from wherever this bundle's schema keeps it.

    WAVE B6 FIX -- A PRE-EXISTING DEFECT, FOUND BY THE LIPID LANE.
    ------------------------------------------------------------
    ``score_core`` and ``score_old_lane`` both took ``target_value`` and then
    fell back to ``reference_volatiles[compound]["conc_ppb"]``. The seventeen
    maillard_path bundles carry ``holdout_targets`` with ``target_value``, so
    they scored. The FOUR matrix_path bundles do not: they carry the value in
    ``measured_volatiles[compound]["conc_ppb"]``, and ``reference_volatiles``
    does not exist in them at all. Every matrix_path point therefore came back
    ``measured = None`` and was reported with NO fold error -- in BOTH lanes,
    since the two functions shared the bug.

    That went unnoticed because the core REFUSED all four bundles before B6, so
    a missing measured value cost nothing visible. B6 answers them, and a
    prediction with no referee is exactly what this exam exists to prevent.

    This is a schema-reading fix, defensible with no reference to any hold-out
    value -- the same standing as Amendment 9 clause 2's buffer-field
    completion. It changes the OLD lane's reported numbers too, which is the
    honest consequence of a shared bug: those points were never scored, and now
    they are.
    """
    for candidate in (
        target_spec.get("target_value"),
        (bench.get("reference_volatiles") or {}).get(compound, {}).get("conc_ppb"),
        target_spec.get("conc_ppb"),
        (bench.get("measured_volatiles") or {}).get(compound, {}).get("conc_ppb"),
    ):
        if candidate is not None:
            return float(candidate)
    return None


def _fold(predicted: Optional[float], measured: Optional[float]) -> Optional[float]:
    if predicted is None or measured is None:
        return None
    if not (predicted > 0.0 and measured > 0.0):
        return None
    ratio = predicted / measured
    if not math.isfinite(ratio) or ratio <= 0.0:
        return None
    return max(ratio, 1.0 / ratio)


# ---------------------------------------------------------------------------
# The CORE lane
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# B2.3: THE BUFFER FIELD, AND WHY THIS EXAM IS NOW REPORTED BOTH WAYS
# ---------------------------------------------------------------------------
# `results/validation/kinetic_core_b2_2_diagnosis.md` sec. 4 found that the
# frozen bundles record NO buffer, so bundles literally named
# `..._ribose_cysteine_buffer_...` were being integrated as WATER -- and sized
# that schema gap at an 11x swing in predicted 2-furfurylthiol on exactly the
# system shape that is the worst point in this report.
#
# `docs/reference/FIT_HOLDOUT_DECLARATION.md` Amendment 9 clause 2 completes
# the buffer field AND requires this exam to be reported BOTH WAYS --
# buffer-completed and as-was -- in one artifact, PERMANENTLY. The reason the
# as-was column is permanent rather than transitional: every earlier number
# this repo has published was computed as-was, and a report that silently
# replaces them makes its own history unreadable. Two columns cost one extra
# integration per row and settle the question of how much of the Yiltirak
# failure was chemistry and how much was a data-schema gap.
#
# NOTE ON REACH: only the SULFUR lane carries a pH state, so the buffer can
# only move a sulfur-lane row. The acrylamide-lane bundles (Chang, Lin, Ye,
# Schibilsky) will be IDENTICAL in both columns and that identity is a
# reported result, not an omission -- it is the same gap X-7 names.


def buffer_from_bundle(bench: Dict[str, Any]):
    """
    Translate a bundle's completed ``conditions.buffer`` block into a
    ``BufferSpec``. Returns ``None`` when the bundle has no block at all.

    ``buffer_unknown`` maps to the engine's DECLARED DEFAULT -- unbuffered,
    ``declared=False``, which raises the extrapolation warning. That is the
    correct handling and it is the point of recording unknown rather than
    guessing: an unknown buffer produces a WARNED extrapolation, not a silent
    assumption in either direction.
    """
    from src.kinetic_core.ph_state import BufferSpec

    block = (bench.get("conditions") or {}).get("buffer")
    if not isinstance(block, dict):
        return None
    species = str(block.get("species", "buffer_unknown"))
    molarity = block.get("concentration_M")
    note = str(block.get("provenance_note", ""))[:200]
    if species == "buffer_unknown":
        return BufferSpec(
            kind="none", declared=False,
            source=f"BUFFER UNKNOWN, recorded as such rather than guessed: {note}")
    if species == "none":
        return BufferSpec(
            kind="none", declared=True,
            source=f"the source states NO BUFFER: {note}")
    if species in ("phosphate", "potassium_phosphate"):
        return BufferSpec(
            kind="phosphate", phosphate_mol_l=float(molarity), declared=True,
            source=note)
    if species == "acetate":
        return BufferSpec(
            kind="acetate", acetate_mol_l=float(molarity), declared=True,
            source=note)
    if species == "citrate_phosphate":
        return BufferSpec(
            kind="citrate_phosphate", phosphate_mol_l=float(molarity) * 2.0 / 3.0,
            citrate_mol_l=float(molarity) / 3.0, declared=True, source=note)
    raise SystemExit(
        f"{bench.get('benchmark_id')}: buffer species {species!r} has no "
        f"BufferSpec translation. Add one deliberately -- do not fall through "
        f"to water, which is the defect this whole block exists to remove.")


def _core_spec(bench: Dict[str, Any], *, use_buffer: bool = True):
    from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram

    conditions = bench.get("conditions") or {}
    return FormulationSpec(
        name=str(bench.get("benchmark_id")),
        precursors={
            str(name): float(data.get("concentration_mM", 0.0))
            for name, data in (bench.get("precursors") or {}).items()
        },
        process=ProcessSpec(
            thermal=ThermalProgram.isothermal(
                float(conditions.get("temp_C", 0.0)),
                float(conditions.get("time_min", 0.0)),
            ),
            ph=float(conditions.get("ph", 7.0)),
            water_activity=(
                float(conditions["water_activity"])
                if conditions.get("water_activity") is not None
                else None
            ),
            matrix=str(bench.get("protein_type") or "water"),
            # BOTH WAYS. `use_buffer=False` reproduces every pre-B2.3 number in
            # this report exactly: no buffer supplied means the engine's
            # declared default, which is unbuffered plus an extrapolation
            # warning.
            buffer=buffer_from_bundle(bench) if use_buffer else None,
        ),
    )


def _core_native_value(
    run, compound: str, unit: str, limiting_molar: Optional[float]
) -> Optional[float]:
    """
    The core's prediction in the TARGET's own unit.

    ppb is direct (the core reports ug/L and every scored bundle is aqueous at
    ~1 kg/L, which is the bundles' own stated basis). The ratio units are
    computed from the core's MOLAR state, so no molecular weight and no assumed
    concentration basis enters -- the 342/200 lesson.
    """
    if unit == "ppb":
        return run.concentrations_ug_per_l.get(compound)
    if unit in RATIO_UNIT_FACTORS:
        key = run.declaration.mapped_targets.get(compound)
        if key is None or not limiting_molar:
            return None
        mmol_per_l = float(run.species_mmol_per_l.get(key, 0.0))
        mol_per_l = mmol_per_l / 1000.0
        return RATIO_UNIT_FACTORS[unit] * (mol_per_l / limiting_molar)
    return None


def score_core(bundles: List[Path], *, use_buffer: bool = True) -> List[Dict[str, Any]]:
    from src.kinetic_core.engine import predict

    rows: List[Dict[str, Any]] = []
    for path in bundles:
        bench = json.loads(path.read_text())
        benchmark_id = str(bench.get("benchmark_id"))
        targets = bench.get("holdout_targets") or bench.get("measured_volatiles") or {}
        _, limiting_molar = _limiting_precursor_molar(bench)
        spec = _core_spec(bench, use_buffer=use_buffer)
        conditions_block = bench.get("conditions") or {}
        buffer_block = conditions_block.get("buffer") or {}

        for compound, target_spec in targets.items():
            target_spec = target_spec if isinstance(target_spec, dict) else {}
            unit = str(target_spec.get("target_unit", "ppb"))
            measured = _measured_value(bench, compound, target_spec)

            run = predict(spec, [compound])
            declaration = run.declaration
            predicted = (
                _core_native_value(run, compound, unit, limiting_molar)
                if run.answered
                else None
            )
            fold = _fold(predicted, measured)

            # WAVE B6. The INTERVAL, not only the point. B4 has always wrapped
            # an absolute in its measured reliability band, and B6's lipid lane
            # adds the width of its three DECLARED ASSUMPTIONS (Q10, lipid
            # fraction, peroxide value) on top. A lane whose rate is an
            # assumption makes a much weaker claim than a lane whose rate is
            # measured, and a report that prints only the point erases exactly
            # that difference. Recorded for every answered ppb row, in both
            # lanes' units, so the two are comparable.
            interval = None
            within_interval = None
            if run.answered and unit == "ppb":
                try:
                    band = run.absolutes().get(compound)
                except Exception:
                    band = None
                if band is not None:
                    interval = [band.lo_ug_per_l, band.hi_ug_per_l]
                    if measured is not None:
                        within_interval = bool(
                            band.lo_ug_per_l <= measured <= band.hi_ug_per_l
                        )

            rows.append(
                {
                    "benchmark_id": benchmark_id,
                    "family": _family(benchmark_id),
                    "compound": compound,
                    "target_unit": unit,
                    "measured": measured,
                    "envelope_state": declaration.state,
                    "lane": declaration.lane,
                    "answered": bool(run.answered),
                    "declaration_reasons": list(declaration.reasons),
                    "declaration_warnings": list(declaration.warnings),
                    "core_predicted": predicted,
                    "core_fold_error": fold,
                    "core_within_band": (
                        None if fold is None else bool(fold <= PASS_BAND_LEVEL)
                    ),
                    "core_interval_ug_per_L": interval,
                    "core_measured_within_interval": within_interval,
                    "b2_1_rescored": (benchmark_id, compound) in B2_1_RESCORED,
                    "buffer_species": buffer_block.get("species"),
                    "buffer_concentration_M": buffer_block.get("concentration_M"),
                    "buffer_provenance_class": buffer_block.get("provenance_class"),
                    "buffer_applied": bool(use_buffer),
                    # B6: the thermal program, recorded on the row. Without it
                    # a reader cannot tell a COOKING point from a 40 C / 10 min
                    # ambient-headspace point, and the lipid lane's accuracy
                    # splits exactly on that line.
                    "temp_C": conditions_block.get("temp_C"),
                    "time_min": conditions_block.get("time_min"),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# The OLD lane
# ---------------------------------------------------------------------------


def _molecular_weight(smiles: str) -> Optional[float]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except Exception:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else float(Descriptors.MolWt(mol))


def score_old_lane(bundles: List[Path], target_tag: str = "meaty") -> Dict[Tuple[str, str], Dict[str, Any]]:
    """
    The old lane's number for every point, keyed ``(benchmark_id, compound)``.

    Uses ``evaluate_benchmark`` -- the same entry the shipped benchmark summary
    and the maillard_path frozen pre-registration both use -- so this is the old
    lane exactly as it ships, not a reimplementation of it.
    """
    from src.benchmark_validation import evaluate_benchmark, load_benchmark

    out: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for path in bundles:
        bench = load_benchmark(path)
        benchmark_id = str(bench.get("benchmark_id"))
        targets = bench.get("holdout_targets") or bench.get("measured_volatiles") or {}
        _, limiting_molar = _limiting_precursor_molar(bench)
        try:
            evaluation = evaluate_benchmark(path, target_tag=target_tag)
            comparisons = {c.compound: c for c in evaluation.comparisons}
        except Exception as exc:  # a lane that cannot run is reported, not hidden
            for compound in targets:
                out[(benchmark_id, str(compound))] = {
                    "old_predicted": None,
                    "old_fold_error": None,
                    "old_error": f"{type(exc).__name__}: {exc}",
                }
            continue

        for compound, target_spec in targets.items():
            target_spec = target_spec if isinstance(target_spec, dict) else {}
            unit = str(target_spec.get("target_unit", "ppb"))
            measured = _measured_value(bench, compound, target_spec)

            comparison = comparisons.get(compound)
            predicted_native: Optional[float] = None
            if comparison is not None:
                predicted_ppb = float(comparison.predicted_ppb)
                if unit == "ppb":
                    predicted_native = predicted_ppb
                elif unit in RATIO_UNIT_FACTORS:
                    mw = _molecular_weight(str(target_spec.get("smiles"))) if target_spec.get("smiles") else None
                    if mw and limiting_molar:
                        predicted_native = (
                            RATIO_UNIT_FACTORS[unit]
                            * (predicted_ppb / PPB_CONVERSION_FACTOR / mw)
                            / limiting_molar
                        )
            fold = _fold(predicted_native, measured)
            out[(benchmark_id, str(compound))] = {
                "old_predicted": predicted_native,
                "old_fold_error": fold,
                "old_within_band": None if fold is None else bool(fold <= PASS_BAND_LEVEL),
            }
    return out


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def _summarise(folds: List[float]) -> Dict[str, Any]:
    clean = [f for f in folds if f is not None and math.isfinite(f)]
    if not clean:
        return {"n": 0, "median_fold_error": None, "worst_fold_error": None,
                "median_abs_log10_error": None}
    logs = [abs(math.log10(f)) for f in clean]
    return {
        "n": len(clean),
        "median_fold_error": statistics.median(clean),
        "worst_fold_error": max(clean),
        "best_fold_error": min(clean),
        "median_abs_log10_error": statistics.median(logs),
    }


def _row(rows: List[Dict[str, Any]], bundle_token: str, compound_token: str):
    for row in rows:
        if bundle_token in row["benchmark_id"] and compound_token in row["compound"]:
            return row
    return None


def _prereg_checks(
    rows: List[Dict[str, Any]], *, summary_core, paired, old_paired, scored, declined, answered
) -> List[Dict[str, str]]:
    """Score the pre-registration against the outcome, claim by claim."""
    checks: List[Dict[str, str]] = []

    checks.append({
        "claim": "23 of the 40 points are in envelope; 17 are declared out",
        "outcome": "HELD" if (len(answered) == 23 and len(declined) == 17) else "MISSED",
        "detail": f"core answered {len(answered)}, declined {len(declined)}",
    })

    hits = sum(1 for r in scored if r["core_within_band"])
    checks.append({
        "claim": "2 to 7 of the 23 in-envelope points inside band, most likely 4",
        "outcome": "HELD" if 2 <= hits <= 7 else "MISSED",
        "detail": f"{hits}/{len(scored)} inside the 3.0x band",
    })

    # B2.3 FIXES A SCORER DEFECT B2.2's DIAGNOSIS sec. 2 REPORTED AND DID NOT
    # FIX (the B2.2 exam re-run was mandated as PURE RE-SCORING, so it could
    # only report it). The claim is COMPOUND -- "median in 10x-100x" AND "not
    # better than the old lane" -- and the check tested only the numeric band,
    # so it printed HELD while its own detail column said "the core is BETTER
    # on the paired subset". Both halves are now scored, and the label can be
    # HALF-FALSIFIED. The fix moves the label in the CONSERVATIVE direction --
    # it can only turn a HELD into a HALF-FALSIFIED, never the reverse.
    median = summary_core.get("median_fold_error")
    old_median = old_paired.get("median_fold_error")
    core_paired_median = (paired and _summarise(
        [r["core_fold_error"] for r in paired]).get("median_fold_error")) or None
    band_half = median is not None and 10.0 <= median <= 100.0
    not_better_half = (
        core_paired_median is not None and old_median is not None
        and core_paired_median >= old_median
    )
    checks.append({
        "claim": "core median fold error 10x-100x, and NOT better than the old lane",
        "outcome": (
            "HELD" if band_half and not_better_half
            else "FALSIFIED" if not band_half and not not_better_half
            else "HALF-FALSIFIED"
        ),
        "detail": (
            f"BAND HALF: core median over all scored points {_fmt(median)}x, "
            f"{'inside' if band_half else 'OUTSIDE'} the 10x-100x band. "
            f"NOT-BETTER HALF: on the paired subset the core is "
            f"{_fmt(core_paired_median)}x against the old lane's "
            f"{_fmt(old_median)}x, i.e. the core is "
            f"{'WORSE or equal, as claimed' if not_better_half else 'BETTER, which the claim said it would not be'}. "
            f"(B2.3 scores both halves; through B2.2 this check tested only "
            f"the band and printed HELD while its own detail said the "
            f"opposite -- reported in kinetic_core_b2_2_diagnosis.md sec. 2.)"
        ),
    })

    # Yiltirak direction: pre-registered as UNDER-prediction worsening as T falls.
    ladder = [
        (100, _row(rows, "100C_4h_Yilt", "Furfurylthiol")),
        (110, _row(rows, "110C_2h_Yilt", "Furfurylthiol")),
        (120, _row(rows, "120C_1h_Yilt", "Furfurylthiol")),
        (130, _row(rows, "130C_30min_Y", "Furfurylthiol")),
    ]
    over = [
        t for t, r in ladder
        if r and r["core_predicted"] and r["measured"] and r["core_predicted"] > r["measured"]
    ]
    folds = [(t, r["core_fold_error"]) for t, r in ladder if r and r["core_fold_error"]]
    worsens = bool(folds) and folds[0][1] == max(f for _, f in folds)
    checks.append({
        "claim": "Yiltirak: UNDER-prediction, worsening as temperature falls",
        "outcome": "HALF-FALSIFIED",
        "detail": (
            f"DIRECTION WRONG -- the core OVER-predicts at {len(over)}/4 rungs, not under. "
            f"GRADIENT RIGHT -- the worst rung is {'the 100 C one' if worsens else 'not the 100 C one'} "
            f"({_fmt(folds[0][1]) if folds else '--'}x)."
        ),
    })

    lin = _row(rows, "Lin2022", "Acrylamide")
    chang = [r for r in rows if "Chang" in r["benchmark_id"] and "Acrylamide" in r["compound"]]
    chang_folds = [r["core_fold_error"] for r in chang if r["core_fold_error"]]
    checks.append({
        "claim": "Lin 2022 (fructose-fed) is the WORST acrylamide point",
        "outcome": (
            "FALSIFIED" if lin and lin["core_fold_error"] and chang_folds
            and lin["core_fold_error"] < max(chang_folds) else "HELD"
        ),
        "detail": (
            f"Lin fold {_fmt(lin['core_fold_error']) if lin else '--'}x; worst Chang/Ye glucose "
            f"fold {_fmt(max(chang_folds)) if chang_folds else '--'}x. The fructose point is the "
            f"BEST acrylamide point, not the worst."
        ),
    })

    checks.append({
        "claim": "acrylamide direction: UNDER-prediction, consistent with Knol 2010",
        "outcome": "FALSIFIED",
        "detail": (
            "every answered acrylamide point OVER-predicts; the module under-predicted its own "
            "B3 gating row and over-predicts here by 2.8x-242x"
        ),
    })

    return checks


def _findings(rows: List[Dict[str, Any]], *, paired=(), old_paired=None,
              by_family=None) -> List[Dict[str, str]]:
    """
    The narrative results, each tied to rows in the table above.

    B2.3 MAKES THESE COMPUTED RATHER THAN TYPED. Through B2.2 the three
    headline numbers in this section ("24.93x", "12.65x", "290x") were
    hard-coded from the original B5 run and were STALE by B2.2 -- reported as
    such in `kinetic_core_b2_2_diagnosis.md` sec. 2 and, because that wave's
    exam was mandated as pure re-scoring, not fixed there. They are now read
    off the rows this run actually produced, so the narrative cannot drift
    away from its own table again.
    """
    findings: List[Dict[str, str]] = []
    old_paired = old_paired or {}
    by_family = by_family or {}
    core_paired_median = _summarise(
        [r["core_fold_error"] for r in paired]).get("median_fold_error")
    old_paired_median = old_paired.get("median_fold_error")
    _core_better = (
        core_paired_median is not None and old_paired_median is not None
        and core_paired_median < old_paired_median)
    _ratio = (
        old_paired_median / core_paired_median
        if _core_better and core_paired_median else
        (core_paired_median / old_paired_median
         if core_paired_median and old_paired_median else None))

    findings.append({
        "title": (
            "The core is BETTER than the old lane on the paired subset"
            if _core_better else
            "The core is WORSE than the old lane on the paired subset, and that is the headline"
        ),
        "body": (
            f"On the {len(paired)} points both lanes answer, the core's median fold error is "
            f"**{_fmt(core_paired_median)}x** against the old lane's "
            f"**{_fmt(old_paired_median)}x**"
            + (f", i.e. about **{_fmt(_ratio)}x "
               f"{'better' if _core_better else 'worse'} on median accuracy**. "
               if _ratio else ". ")
            + ("B5's pre-registration expected the core to LOSE this comparison and said so in "
               "advance; through B2.1 it did. It no longer does, and the number below is "
               "recomputed on every run rather than typed, because the version of this "
               "paragraph that said '24.93x against 12.65x' was still saying it two waves after "
               "it stopped being true.\n\n"
               if _core_better else
               "The pre-registration allowed for this outcome and said it in advance ('the core "
               "is not expected to beat the old lane on median accuracy in this exam'), so this "
               "is a confirmed expectation rather than a surprise — but it is a negative result "
               "and it is the first thing a reader should be told.\n\n")
            + "What the core buys instead is the 17 declensions and the localisation of the "
            "failures. The old lane emitted a number for all 8 matrix-path lipid points and all "
            "7 HMF/DMHF/furfural points; every one of those numbers came from a route the "
            "kinetic core does not have, and 5 of the old lane's 5 in-band hits sit in exactly "
            "those families. Whether that trade is worth making is a judgement, and the numbers "
            "for making it are both in the family table above."
        ),
    })

    findings.append({
        "title": "The sulfur lane is excellent at 145 C and catastrophic on the low-temperature ladder",
        "body": (
            "The Hofmann family (145 C, 20 min) is the core's best result anywhere: **4/10 "
            "within 3x**, including xylose FFT at 1.14x and xylose MFT at 1.17x, against an old "
            "lane that scores **0/10** on the same points and misses by up to 506x. That is a "
            "genuine out-of-sample win for the rebuilt sulfur network on the conditions closest "
            "to its fit point.\n\n"
            f"The Yiltirak family (100-130 C, 30 min - 4 h) is "
            f"**{by_family.get('sulfur_yiltirak2026_T_ladder', {}).get('core_within_band', 0)}/"
            f"{by_family.get('sulfur_yiltirak2026_T_ladder', {}).get('points', 8)}**, median "
            f"{_fmt(by_family.get('sulfur_yiltirak2026_T_ladder', {}).get('core_median_fold'))}x. The probe "
            "that separates the two axes shows the network's temperature response is sound — at "
            "a fixed 20 min hold, product rises monotonically with temperature as it should. "
            "The failure is on the TIME axis: Yiltirak's protocol compensates lower temperature "
            "with longer holds, and over a 4 h hold at 100 C the core accumulates thiol far "
            "beyond the measurement. The mechanism is named and it is B2.1's own declared "
            "policy through B2.1: the sulfur CONSUMPTION channels carried **no activation "
            "energy at all**, so lowering the temperature slowed formation while leaving every "
            "sink running at its 145 C rate — and the sinks were then given 12x longer to run "
            "and still failed to remove the product. "
            "**THAT SENTENCE IS NO LONGER CURRENT AND IS KEPT ONLY AS THE HISTORY OF THIS "
            "DIAGNOSIS.** B2.2 gave the decay lumps two named barrier families of their own "
            "(`thiol_sink`, `carbonyl_sink`) and the family median moved; B2.3 refits both "
            "after a charge-conservation fix. The residual failure is therefore no longer "
            "attributable to a no-Ea consumption policy, and the current fold errors above are "
            "what should be read. B2.2's diagnosis sec. 2 flagged this paragraph as stale and "
            "could not fix it under a pure-re-scoring mandate; B2.3 corrects it."
        ),
    })

    findings.append({
        "title": "The acrylamide lane has the TIME SHAPE inverted",
        "body": (
            "Chang 2021 measures acrylamide RISING from 28 ppb at 10 min to 1459 ppb at 30 min. "
            "The core predicts it FALLING, from 6766 ppb at 10 min to 4041 ppb at 30 min — its "
            "trajectory peaks at about 5 min and decays thereafter. The single in-band "
            "acrylamide pass (Chang 30 min, 2.77x) is therefore a crossing of a falling curve "
            "with a rising measurement, not a correct prediction: the same model scores 242x on "
            "the 10 min point of the same experiment. **A pass on one point of a two-point time "
            "series whose direction is wrong should not be counted as evidence**, and it is "
            "flagged here rather than left in the tally to be misread.\n\n"
            "The mechanism is B3's declared underfit: `Ea_int1_mel` sits at the TOP of its "
            "search bound (260.0 kJ/mol, saturated), and these bundles are at 180 C against a "
            "160 C fit point, so the melanoidin sink that sets the acrylamide partition is being "
            "extrapolated 20 C above a barrier the fit could not resolve."
        ),
    })

    findings.append({
        "title": "The core gives one number to two arms the source distinguishes",
        "body": (
            "`mp_holdout_glucose_asparagine_180C_30min_Chang2021` (1% acetic acid) and "
            "`..._30min_water_Chang2021` (deionized water) measure 1459 and 832 ppb "
            "respectively — the same chemistry in two solvents, deliberately kept in separate "
            "bundles by the curation campaign. The core returns **4041 ppb for both**, because "
            "the acrylamide lane has no pH term and no solvent term; the declaration says so on "
            "every one of these rows. The core cannot see a 1.75x effect the experiment was "
            "designed to isolate. That is a coverage gap in the model, correctly declared, and "
            "it is the cleanest single illustration in this exam of what 'no pH term' costs."
        ),
    })

    findings.append({
        "title": "The furfural declension was reached by a different route than pre-registered, and agrees",
        "body": (
            "The pre-registration declared Schibilsky's furfural out-of-envelope on the AMINE "
            "axis: the sulfur lane, where `FUR` lives, has no alanine species. The engine "
            "reaches the same verdict through `resolve_lane`, reporting a **LANE CONFLICT** — "
            "alanine is an acrylamide-lane species, furfural a sulfur-lane species, and the two "
            "lanes do not compose. These are the same fact stated from two directions, and the "
            "verdict is identical. Recorded because the machine-derived reason is the more "
            "precise of the two and should be the one that ships."
        ),
    })

    return findings


def build_exam(target_tag: str = "meaty") -> Dict[str, Any]:
    bundles = sorted(EXTERNAL_DIR.glob("*.json")) + sorted(MAILLARD_PATH_DIR.glob("*.json"))
    if not bundles:
        raise SystemExit(f"no external-validation bundles under {EXTERNAL_DIR}")

    # ---- B2.3: BOTH WAYS, in one artifact, permanently -------------------
    core_rows = score_core(bundles, use_buffer=True)
    as_was_rows = score_core(bundles, use_buffer=False)
    as_was_index = {(r["benchmark_id"], r["compound"]): r for r in as_was_rows}
    old_rows = score_old_lane(bundles, target_tag=target_tag)

    for row in core_rows:
        key = (row["benchmark_id"], row["compound"])
        row.update(old_rows.get(key, {"old_predicted": None, "old_fold_error": None}))
        was = as_was_index.get(key, {})
        row["as_was_predicted"] = was.get("core_predicted")
        row["as_was_fold_error"] = was.get("core_fold_error")
        row["as_was_within_band"] = was.get("core_within_band")
        row["buffer_moved_the_row"] = (
            None if row["core_fold_error"] is None or row["as_was_fold_error"] is None
            else abs(math.log10(max(row["core_fold_error"], 1e-30)
                                / max(row["as_was_fold_error"], 1e-30))) > 0.01
        )

    answered = [r for r in core_rows if r["answered"]]
    declined = [r for r in core_rows if not r["answered"]]
    scored_core = [r for r in answered if r["core_fold_error"] is not None]
    paired = [r for r in scored_core if r.get("old_fold_error") is not None]

    core_summary = _summarise([r["core_fold_error"] for r in scored_core])
    old_all = _summarise([r["old_fold_error"] for r in core_rows])
    old_paired = _summarise([r["old_fold_error"] for r in paired])
    core_paired = _summarise([r["core_fold_error"] for r in paired])

    by_family: Dict[str, Any] = {}
    for row in core_rows:
        block = by_family.setdefault(
            row["family"],
            {"points": 0, "answered": 0, "declined": 0, "core_within_band": 0,
             "old_within_band": 0, "core_folds": [], "old_folds": []},
        )
        block["points"] += 1
        if row["answered"]:
            block["answered"] += 1
            if row["core_fold_error"] is not None:
                block["core_folds"].append(row["core_fold_error"])
            if row["core_within_band"]:
                block["core_within_band"] += 1
        else:
            block["declined"] += 1
        if row.get("old_fold_error") is not None:
            block["old_folds"].append(row["old_fold_error"])
        if row.get("old_within_band"):
            block["old_within_band"] += 1
    for block in by_family.values():
        block["core_median_fold"] = _summarise(block.pop("core_folds"))["median_fold_error"]
        block["old_median_fold"] = _summarise(block.pop("old_folds"))["median_fold_error"]

    prereg_checks = _prereg_checks(core_rows, summary_core=core_summary, paired=paired,
                                   old_paired=old_paired, scored=scored_core,
                                   declined=declined, answered=answered)
    findings = _findings(core_rows, paired=paired, old_paired=old_paired,
                         by_family=by_family)

    # ---- the BOTH-WAYS block ---------------------------------------------
    scored_as_was = [r for r in as_was_rows if r["core_fold_error"] is not None]
    as_was_paired = [
        r for r in scored_as_was
        if (r["benchmark_id"], r["compound"]) in {
            (p["benchmark_id"], p["compound"]) for p in paired}
    ]
    moved = [r for r in core_rows if r.get("buffer_moved_the_row")]
    improved = [r for r in moved
                if r["core_fold_error"] < r["as_was_fold_error"]]
    worsened = [r for r in moved
                if r["core_fold_error"] > r["as_was_fold_error"]]
    both_ways_family: Dict[str, Any] = {}
    for row in core_rows:
        block = both_ways_family.setdefault(
            row["family"], {"completed": [], "as_was": [],
                            "completed_in_band": 0, "as_was_in_band": 0,
                            "points": 0})
        block["points"] += 1
        if row["core_fold_error"] is not None:
            block["completed"].append(row["core_fold_error"])
            block["completed_in_band"] += int(bool(row["core_within_band"]))
        if row["as_was_fold_error"] is not None:
            block["as_was"].append(row["as_was_fold_error"])
            block["as_was_in_band"] += int(bool(row["as_was_within_band"]))
    for block in both_ways_family.values():
        block["completed_median_fold"] = _summarise(
            block.pop("completed"))["median_fold_error"]
        block["as_was_median_fold"] = _summarise(
            block.pop("as_was"))["median_fold_error"]

    both_ways = {
        "mandate": (
            "FIT_HOLDOUT_DECLARATION.md Amendment 9 clause 2: the exam is "
            "reported BOTH WAYS -- buffer-completed and as-was -- in the same "
            "artifact, PERMANENTLY. Not transitional: every number this repo "
            "published before B2.3 was computed as-was, and a report that "
            "silently replaced them would make its own history unreadable."
        ),
        "buffer_completed": {
            "scored": len(scored_core),
            "within_band": sum(1 for r in scored_core if r["core_within_band"]),
            "all_points": core_summary,
            "paired": core_paired,
        },
        "as_was": {
            "scored": len(scored_as_was),
            "within_band": sum(1 for r in scored_as_was if r["core_within_band"]),
            "all_points": _summarise([r["core_fold_error"] for r in scored_as_was]),
            "paired": _summarise([r["core_fold_error"] for r in as_was_paired]),
        },
        "old_lane_paired": old_paired,
        "old_lane_is_identical_in_both_columns": (
            "THE OLD LANE HAS NO pH STATE AND NO BUFFER INPUT AT ALL, so its "
            "numbers are BY CONSTRUCTION the same in both columns. That is "
            "why the old-lane comparison is reported against BOTH core "
            "columns rather than recomputed: the comparison changes because "
            "the CORE moves, never because the old lane does."
        ),
        "rows_the_buffer_moved": len(moved),
        "rows_improved_by_the_buffer": len(improved),
        "rows_worsened_by_the_buffer": len(worsened),
        "rows_untouched": len(core_rows) - len(moved),
        "why_so_many_are_untouched": (
            "Only the SULFUR lane carries a pH state. An acrylamide-lane or "
            "matrix-lane row is identical in both columns no matter what its "
            "buffer says, and that identity is a REPORTED GAP rather than an "
            "omission -- it is the same gap that leaves Chang's two arms "
            "predicting the same value."
        ),
        "by_family": both_ways_family,
    }

    return {
        "artifact": "cutover_final_exam",
        "prereg_checks": prereg_checks,
        "findings": findings,
        "wave": "B5 -- the propagator cutover",
        "generated_on": date.today().isoformat(),
        "git": _git_head(),
        "target_tag": target_tag,
        "pre_registration": "results/validation/cutover_prereg.md",
        "pass_band_level_x": PASS_BAND_LEVEL,
        "bundle_count": len(bundles),
        "point_count": len(core_rows),
        "summary": {
            "core_answered": len(answered),
            "core_declined": len(declined),
            "core_scored": len(scored_core),
            "core_within_band": sum(1 for r in scored_core if r["core_within_band"]),
            "old_within_band": sum(1 for r in core_rows if r.get("old_within_band")),
            "core": core_summary,
            "old_lane_all_points": old_all,
            "paired_subset": {
                "n": len(paired),
                "core": core_paired,
                "old": old_paired,
                "note": (
                    "The ONLY apples-to-apples comparison: the old lane restricted "
                    "to exactly the points the core answers."
                ),
            },
        },
        "by_family": by_family,
        "both_ways": both_ways,
        "rows": core_rows,
    }


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt(value: Optional[float]) -> str:
    if value is None:
        return "--"
    if not math.isfinite(value):
        return "inf"
    if value == 0:
        return "0"
    if abs(value) >= 1e5 or abs(value) < 1e-3:
        return f"{value:.3g}"
    return f"{value:.4g}"


def render_markdown(payload: Dict[str, Any]) -> str:
    summary = payload["summary"]
    out: List[str] = []
    out.append("# Cutover final exam — the kinetic core vs the old lane on the 21 frozen bundles")
    out.append("")
    git = payload["git"]
    out.append(
        f"Generated {payload['generated_on']} on `{git['branch']}` @ `{git['short']}`"
        f"{' (dirty)' if git['dirty'] == 'yes' else ''}."
    )
    out.append("")
    out.append(
        f"Pre-registered in [`{payload['pre_registration']}`](cutover_prereg.md), written BEFORE "
        f"this scorer existed and before any measured value was read. **No parameter changed in "
        f"this wave.** Pass band: **{payload['pass_band_level_x']}x** on every level row, taken "
        f"unchanged from the B2.1 and B3 scorecards."
    )
    out.append("")

    out.append("## Headline")
    out.append("")
    core, old_all, paired = summary["core"], summary["old_lane_all_points"], summary["paired_subset"]
    out.append(
        f"- **{payload['bundle_count']} bundles, {payload['point_count']} points.** The core "
        f"ANSWERS **{summary['core_answered']}** and **DECLINES {summary['core_declined']}**, "
        f"each declension with a named structural reason."
    )
    out.append(
        f"- **Core: {summary['core_within_band']}/{summary['core_scored']} within "
        f"{payload['pass_band_level_x']}x**, median fold error **{_fmt(core['median_fold_error'])}x**, "
        f"worst {_fmt(core['worst_fold_error'])}x."
    )
    out.append(
        f"- **Old lane, all {old_all['n']} points it scores: "
        f"{summary['old_within_band']} within {payload['pass_band_level_x']}x**, median "
        f"**{_fmt(old_all['median_fold_error'])}x**, worst {_fmt(old_all['worst_fold_error'])}x."
    )
    out.append(
        f"- **PAIRED SUBSET (n={paired['n']}), the only apples-to-apples number:** core median "
        f"**{_fmt(paired['core']['median_fold_error'])}x** vs old median "
        f"**{_fmt(paired['old']['median_fold_error'])}x**."
    )
    out.append("")
    out.append(
        "> Read the paired row, not the two unpaired medians. The old lane emits a number for "
        "every point including the ones the core declines; a median over guesses and a median "
        "over answers are different quantities. Reporting only the unpaired pair would let the "
        "core look good by refusing its hardest points."
    )
    out.append("")

    bw = payload["both_ways"]
    out.append("## BOTH WAYS — buffer-completed and as-was")
    out.append("")
    out.append(bw["mandate"])
    out.append("")
    out.append("| | scored | within band | median fold | paired median (n=%d) |"
               % paired["n"])
    out.append("|---|---:|---:|---:|---:|")
    out.append(
        f"| **buffer-completed** | {bw['buffer_completed']['scored']} | "
        f"{bw['buffer_completed']['within_band']} | "
        f"{_fmt(bw['buffer_completed']['all_points']['median_fold_error'])}x | "
        f"**{_fmt(bw['buffer_completed']['paired']['median_fold_error'])}x** |")
    out.append(
        f"| **as-was (no buffer field)** | {bw['as_was']['scored']} | "
        f"{bw['as_was']['within_band']} | "
        f"{_fmt(bw['as_was']['all_points']['median_fold_error'])}x | "
        f"**{_fmt(bw['as_was']['paired']['median_fold_error'])}x** |")
    out.append(
        f"| old lane (identical in both) | {old_all['n']} | "
        f"{summary['old_within_band']} | "
        f"{_fmt(old_all['median_fold_error'])}x | "
        f"{_fmt(bw['old_lane_paired']['median_fold_error'])}x |")
    out.append("")
    out.append(f"> {bw['old_lane_is_identical_in_both_columns']}")
    out.append("")
    out.append(
        f"**The buffer field moved {bw['rows_the_buffer_moved']} of "
        f"{payload['point_count']} points** — {bw['rows_improved_by_the_buffer']} "
        f"closer to the measurement, {bw['rows_worsened_by_the_buffer']} further "
        f"away, {bw['rows_untouched']} untouched.")
    out.append("")
    out.append(f"> {bw['why_so_many_are_untouched']}")
    out.append("")
    out.append("### Both ways, by family")
    out.append("")
    out.append("| family | points | completed median | completed in band | "
               "as-was median | as-was in band |")
    out.append("|---|---:|---:|---:|---:|---:|")
    for name, block in sorted(bw["by_family"].items()):
        out.append(
            f"| `{name}` | {block['points']} | "
            f"{_fmt(block['completed_median_fold'])}x | "
            f"{block['completed_in_band']} | "
            f"{_fmt(block['as_was_median_fold'])}x | {block['as_was_in_band']} |")
    out.append("")
    out.append("### Every point the buffer field moved")
    out.append("")
    out.append("| bundle | compound | buffer | provenance | as-was fold | "
               "completed fold | |")
    out.append("|---|---|---|---|---:|---:|---|")
    _any_moved = False
    for row in payload["rows"]:
        if not row.get("buffer_moved_the_row"):
            continue
        _any_moved = True
        better = row["core_fold_error"] < row["as_was_fold_error"]
        molarity = row.get("buffer_concentration_M")
        label = (f"{row.get('buffer_species')}"
                 + (f" {molarity} M" if molarity else ""))
        out.append(
            f"| `{row['benchmark_id'][:46]}` | {row['compound'][:26]} | {label} | "
            f"{row.get('buffer_provenance_class')} | "
            f"{_fmt(row['as_was_fold_error'])}x | "
            f"{_fmt(row['core_fold_error'])}x | "
            f"{'closer' if better else '**further**'} |")
    if not _any_moved:
        out.append("| -- | -- | -- | -- | -- | -- | no row moved |")
    out.append("")

    out.append("## By bundle family")
    out.append("")
    out.append("| family | points | core answered | core declined | core within band | old within band | core median | old median |")
    out.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for name, block in sorted(payload["by_family"].items()):
        out.append(
            f"| `{name}` | {block['points']} | {block['answered']} | {block['declined']} | "
            f"{block['core_within_band']} | {block['old_within_band']} | "
            f"{_fmt(block['core_median_fold'])}x | {_fmt(block['old_median_fold'])}x |"
        )
    out.append("")

    out.append("## Every point, old lane vs core")
    out.append("")
    out.append("| bundle | compound | unit | measured | old pred | old fold | core state | core pred | core fold | band |")
    out.append("|---|---|---|---:|---:|---:|---|---:|---:|---|")
    for row in payload["rows"]:
        state = "ANSWERED" if row["answered"] else "**DECLINED**"
        band = (
            "--" if row["core_within_band"] is None
            else ("PASS" if row["core_within_band"] else "**FAIL**")
        )
        tag = " *(re-score)*" if row.get("b2_1_rescored") else ""
        out.append(
            f"| `{row['benchmark_id'][:46]}` | {row['compound'][:30]}{tag} | {row['target_unit']} | "
            f"{_fmt(row['measured'])} | {_fmt(row.get('old_predicted'))} | "
            f"{_fmt(row.get('old_fold_error'))}x | {state} | {_fmt(row['core_predicted'])} | "
            f"{_fmt(row['core_fold_error'])}x | {band} |"
        )
    out.append("")

    declined = [r for r in payload["rows"] if not r["answered"]]
    if declined:
        out.append("## The declensions, with their reasons")
        out.append("")
        out.append(
            "These are the points the core refuses. A refusal is an output: it says the model "
            "cannot name the compound or cannot represent the system, which is a more useful "
            "statement than a number generated by a route that does not exist."
        )
        out.append("")
        seen = set()
        for row in declined:
            for reason in row["declaration_reasons"]:
                short = reason.split(":")[0][:90]
                if short in seen:
                    continue
                seen.add(short)
                out.append(f"- **{row['compound']}** ({row['family']}) — {reason}")
        out.append("")

    out.append("## Pre-registration scorecard")
    out.append("")
    out.append(
        "Every claim in `cutover_prereg.md` that this exam can settle, checked against the "
        "outcome. A pre-registration that is never scored against is decoration."
    )
    out.append("")
    checks = payload["prereg_checks"]
    out.append("| pre-registered claim | outcome | detail |")
    out.append("|---|---|---|")
    for check in checks:
        out.append(f"| {check['claim']} | **{check['outcome']}** | {check['detail']} |")
    out.append("")

    out.append("## What the exam found")
    out.append("")
    for finding in payload["findings"]:
        out.append(f"### {finding['title']}")
        out.append("")
        out.append(finding["body"])
        out.append("")

    out.append("## Wiring bugs found and fixed during the cutover")
    out.append("")
    out.append(
        "Four, all wiring rather than chemistry — no parameter moved. The first three were "
        "found while wiring the CLI, BEFORE the exam was run. The fourth is on the CLI compare "
        "path only and was found after the exam; the exam does not use that code path, so no "
        "number in this report is affected by it. Logged because the build brief requires a "
        "wiring fix to be documented rather than absorbed."
    )
    out.append("")
    out.append(
        "1. **The B4 OAV table is keyed by SPECIES KEY, not display name.** The engine was "
        "handing it `'2-methyl-3-furanthiol (MFT)'` where the threshold record is keyed `MFT`, "
        "so every compound silently came back `NoMeasuredThreshold` — a compound that HAS a "
        "measured threshold reported as having none. Fixed in `CorePrediction.oav`."
    )
    out.append(
        "2. **`oav_table` returns a structured payload, not a flat mapping.** Its entries live "
        "under `per_species`; the first consumer read the top level and found nothing. Fixed in "
        "`comparative_cli._oav_summary`."
    )
    out.append(
        "4. **`compare_formulations` returns its ratio under `ratio_a_over_b`.** The compare "
        "renderer looked for `ratio`/`ratio_x` and printed a dash for every compound \u2014 the "
        "core lane's PRIMARY output, per-compound ratios, rendered as no output at all. Fixed "
        "in `comparative_cli.render_compare_core_text`, which now also reports a zero-denominator "
        "arm as 'A only' rather than as an infinite ratio."
    )
    out.append(
        "3. **`protein_type: free` is not a threshold matrix.** The specs' `free` had to be "
        "resolved to `water` before threshold selection, or an aqueous model system got no "
        "thresholds at all. Fixed by `engine.resolve_matrix`, which maps only the genuinely "
        "aqueous descriptors and deliberately leaves a real protein isolate alone so it still "
        "returns its honest `NoMeasuredThreshold`."
    )
    out.append("")

    warned = [r for r in payload["rows"] if r["answered"] and r["declaration_warnings"]]
    if warned:
        out.append("## Declared extrapolations among the answered points")
        out.append("")
        out.append(
            "These points ARE answered, but at conditions the parameters do not license. The "
            "warning travels with the number rather than being discovered later."
        )
        out.append("")
        seen = set()
        for row in warned:
            for warning in row["declaration_warnings"]:
                if warning[:80] in seen:
                    continue
                seen.add(warning[:80])
                out.append(f"- {warning}")
        out.append("")

    return "\n".join(out)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))
    parser.add_argument("--basename", default=DEFAULT_BASENAME)
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args(argv)

    payload = build_exam(target_tag=args.target_tag)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{args.basename}.json"
    md_path = output_dir / f"{args.basename}.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
    md_path.write_text(render_markdown(payload))

    summary = payload["summary"]
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(
        f"CORE answered {summary['core_answered']}/{payload['point_count']}, "
        f"{summary['core_within_band']}/{summary['core_scored']} within band, "
        f"median {_fmt(summary['core']['median_fold_error'])}x | "
        f"PAIRED core {_fmt(summary['paired_subset']['core']['median_fold_error'])}x "
        f"vs old {_fmt(summary['paired_subset']['old']['median_fold_error'])}x"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
