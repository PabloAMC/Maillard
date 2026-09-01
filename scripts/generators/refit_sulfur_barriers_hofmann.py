#!/usr/bin/env python3
"""Refit the sulfur-branch FAST barriers against Hofmann1998 — and only Hofmann1998.

Context (2026-08-27, Wave H)
---------------------------
Wave G1 replaced the fabricated one-step MFT shortcut
(3-deoxyosone + 2 H2S -> 2-methyl-3-furanthiol) with the accepted three-step route

    Amadori --(2,3-enolisation)--> 1-deoxyosone
            --(cyclodehydration)--> norfuraneol
            --(+ H2S, reductone-mediated)--> 2-methyl-3-furanthiol

(van den Ouweland & Peer 1975, 10.1021/jf60199a045).  Absolute sulfur yields fell
5-40x as a result — Hofmann1998 MFT went from ratio 1.02 to 7.83x UNDER.

WITHDRAWN 2026-08-27 (Wave S2c) -- READ THIS BEFORE THE PARAGRAPH BELOW.
`cys_ribose_140C_Hofmann1998` was NEVER a literature constraint. Wave S2b traced its
MFT 342 / FFT 200 ppb to data/benchmarks/maillard_validation_benchmarks.md section 1.3,
an abstract-reconstructed range table committed in c7efbbc -- the SAME commit that
created the benchmark JSON -- whose row gives MFT `~0.02-0.05` mol % and FFT
`~0.01-0.03` mol %. On the benchmark's declared (unattested) 10 mM basis with MW 114.17:
0.0300 mol % -> 342.5 -> 342 ppb, and the FFT band's geometric mean 0.017321 mol % ->
197.8 -> 200 ppb. Both are interior points of two invented, overlapping bands (~90%
confidence, arithmetic exact). THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.
DO NOT RE-RUN THIS SCRIPT against that benchmark: it needs a real target first (the ILL
pack in tasks/audit_remediation.md "## Wave S2b" section (f), then a rebuild in native
mol %). This script's own output record, results/validation/sulfur_barrier_refit_hofmann.
{json,md}, is annotated rather than retracted, because the constant it moved
(`thiol_addition_norfuraneol`) sits on a family no step emits; its companion
sulfur_barrier_refit_pentodiulose IS retracted and its constant reverted.

The forensics record (tasks/audit_remediation.md, "Re-anchor the WHOLE sulfur
branch") established that after the Mottram1994 / Farmer1999 quarantine,
`cys_ribose_140C_Hofmann1998` is the ONLY surviving literature constraint on the
sulfur branch.  Everything else that touches it is either synthetic
(Internal2026 / ProtocolPilot2026), a hold-out (data/benchmarks/external_validation/)
or has its own per-lane observability history (the xylose HVP lanes).  This script
therefore fits against that one benchmark and asserts the exclusions rather than
merely intending them.

Prerequisite fixed in the same wave
-----------------------------------
Before 2026-08-27 the two knobs this script exists to fit had EXACTLY ZERO
derivative, because `Enolisation_2_3_Amadori` — the step that opens the whole new
route — was accidentally collecting the amine-nucleophile ionisation correction
and the Amadori water-activity correction through substring matching on its own
name (see `src/conditions.py::_releases_rather_than_attacks_with_the_amine`).  At
Hofmann's pH 5.0 / aw 0.98 that was +6.06 kcal/mol, which kept the real route
~1600x below the demoted shortcut.  With that fixed the norfuraneol route becomes
the selected MFT path (at unchanged predictions) and the knobs below become
identifiable.

Knobs and their defensible ranges
---------------------------------
Every range below is bounded by values ALREADY IN `FAST_BARRIERS` for the same
mechanistic class, so no knob can leave the envelope the table itself spans.  All
four are ESTIMATED tier (see their comments in src/barrier_constants.py).

    thiol_addition_norfuraneol  [23.30, 29.65]  H2S/thiol addition to a carbonyl:
                                low = thiohemiacetal_formation (23.30),
                                high = thiol_addition_hexose (29.65)
    furanone_cyclisation        [26.00, 30.00]  cyclodehydration class: dehydration /
                                2,3-enolisation / furanone_formation are all 28.0, +/- 2
    thiohemiacetal_formation    [21.00, 26.80]  bounded above by thiol_dehydration (26.80)
    thiol_dehydration           [23.30, 28.00]  bounded below by thiohemiacetal_formation
                                (23.30), above by dehydration (28.00)

`thiol_addition` is REPORTED but NOT FITTED: after Wave G1 it labels only the
DEMOTED one-step shortcut (`Thiol_Addition_Legacy_Shortcut`) and the lumped
`Thiol_Addition_H2` step.  Fitting it would tune the fabricated route that the
wave exists to retire.

Decision rules (stated up front, applied mechanically)
------------------------------------------------------
1. IDENTIFIABILITY.  A knob whose objective profile is flat to
   `FLAT_PROFILE_DEX` over its whole defensible range is reported as
   NOT IDENTIFIABLE by this benchmark and keeps its incumbent value.
2. MATERIALITY.  A knob moves only if the best achievable objective gain is at
   least `MIN_OBJECTIVE_GAIN_DEX`.  Chasing less than that on a two-row
   benchmark is fitting noise.
3. CONSERVATIVE EDGE.  Among the values within `INDIFFERENCE_DEX` of the profile
   minimum, the LARGEST barrier is adopted — i.e. the least sulfur-favourable
   point the Hofmann data cannot distinguish from the optimum.  This is what
   keeps a saturating profile from dragging a constant to the bottom of its
   range and calling that a fit.

Usage
-----
    python scripts/generators/refit_sulfur_barriers_hofmann.py
    python scripts/generators/refit_sulfur_barriers_hofmann.py --apply   # rewrite the
                                                                        # constants
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src import barrier_constants as barrier_constants_module
from src.barrier_constants import FAST_BARRIERS
from src.benchmark_validation import evaluate_benchmark
from src.uncertainty_propagation import _benchmark_signal_origin

TARGET_BENCHMARK = data_paths.benchmark_path("cys_ribose_140C_Hofmann1998")
BARRIER_SOURCE = ROOT / "src" / "barrier_constants.py"

UNSCORED_ROW_PENALTY = 3.0        # dex — same convention as refit_projection_constants.py
GRID_STEP_KCAL = 0.05
FLAT_PROFILE_DEX = 1.0e-6         # below this a profile is numerically flat
MIN_OBJECTIVE_GAIN_DEX = 0.02     # materiality threshold for moving a constant
INDIFFERENCE_DEX = 0.01           # width of the "cannot distinguish" band

# knob -> (low, high, basis for the range)
KNOBS: Dict[str, Tuple[float, float, str]] = {
    "thiol_addition_norfuraneol": (
        23.30,
        29.65,
        "H2S/thiol addition to a carbonyl; low = thiohemiacetal_formation (23.30), "
        "high = thiol_addition_hexose (29.65), both already in FAST_BARRIERS",
    ),
    "furanone_cyclisation": (
        26.00,
        30.00,
        "cyclodehydration class; dehydration / 2,3-enolisation / furanone_formation "
        "are all 28.0 in FAST_BARRIERS, taken +/- 2 kcal/mol",
    ),
    "thiohemiacetal_formation": (
        21.00,
        26.80,
        "thiohemiacetal formation must stay below its own dehydration step, "
        "thiol_dehydration (26.80)",
    ),
    "thiol_dehydration": (
        23.30,
        28.00,
        "bounded below by thiohemiacetal_formation (23.30) and above by the generic "
        "dehydration barrier (28.00)",
    ),
}

REPORTED_NOT_FITTED = {
    "thiol_addition": (
        "After Wave G1 this key labels only `Thiol_Addition_Legacy_Shortcut` (the "
        "demoted fabricated one-step MFT route) and the lumped `Thiol_Addition_H2` "
        "step. Fitting it would tune the route this wave exists to retire; its "
        "derivative is reported instead."
    ),
}


# ---------------------------------------------------------------------------
# fit-target guards
# ---------------------------------------------------------------------------
def _guarded_target() -> Path:
    path = TARGET_BENCHMARK
    assert path.exists(), f"fit target missing: {path}"
    assert "external_validation" not in path.parts, (
        f"hold-out benchmark reached the fit-target selector: {path}"
    )
    assert "Internal2026" not in path.name and "ProtocolPilot2026" not in path.name, (
        f"synthetic benchmark reached the fit-target selector: {path}"
    )
    origin = _benchmark_signal_origin(path)
    assert origin == "external_literature", (
        f"fit target must be literature-sourced, got signal origin {origin!r}"
    )
    return path


# ---------------------------------------------------------------------------
# objective
# ---------------------------------------------------------------------------
def _score(path: Path) -> Tuple[float, List[Dict[str, Any]]]:
    evaluation = evaluate_benchmark(path)
    rows: List[Dict[str, Any]] = []
    if not evaluation.supported:
        return UNSCORED_ROW_PENALTY, rows
    for comparison in evaluation.comparisons:
        scored = (
            comparison.matched_name is not None
            and comparison.measured_ppb > 0.0
            and comparison.predicted_ppb > 0.0
        )
        error = (
            abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb))
            if scored
            else UNSCORED_ROW_PENALTY
        )
        rows.append(
            {
                "compound": comparison.compound,
                "measured_ppb": comparison.measured_ppb,
                "predicted_ppb": comparison.predicted_ppb,
                "ratio": (
                    comparison.predicted_ppb / comparison.measured_ppb
                    if comparison.measured_ppb
                    else None
                ),
                "scored": scored,
                "abs_log10_error": error,
            }
        )
    if not rows:
        return UNSCORED_ROW_PENALTY, rows
    return sum(row["abs_log10_error"] for row in rows) / len(rows), rows


def _evaluate_at(path: Path, values: Dict[str, float]) -> Tuple[float, List[Dict[str, Any]]]:
    saved = {key: FAST_BARRIERS[key] for key in values}
    try:
        for key, value in values.items():
            FAST_BARRIERS[key] = (float(value), saved[key][1])
        return _score(path)
    finally:
        for key, entry in saved.items():
            FAST_BARRIERS[key] = entry


def _frange(low: float, high: float, step: float) -> List[float]:
    count = int(round((high - low) / step))
    return [round(low + index * step, 6) for index in range(count + 1)]


# ---------------------------------------------------------------------------
# source rewriting
# ---------------------------------------------------------------------------
def _rewrite_constant(text: str, key: str, new_value: float, note: str) -> str:
    pattern = re.compile(
        r'("%s"\s*:\s*\(\s*)(-?\d+(?:\.\d+)?)(\s*,\s*")' % re.escape(key)
    )
    match = pattern.search(text)
    if match is None:
        raise SystemExit(f"could not locate FAST_BARRIERS entry for {key!r}")
    replacement = f"{match.group(1)}{new_value:.2f}{match.group(3)}"
    return text[: match.start()] + replacement + text[match.end():]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument(
        "--apply",
        action="store_true",
        help="rewrite the adopted values into src/barrier_constants.py",
    )
    args = parser.parse_args()

    path = _guarded_target()
    incumbent = {key: float(FAST_BARRIERS[key][0]) for key in KNOBS}

    baseline_objective, baseline_rows = _score(path)
    print(f"Fit target: {path.name}")
    print(f"Incumbent objective J = {baseline_objective:.4f} dex")
    for row in baseline_rows:
        print(
            f"   {row['compound']:26s} measured {row['measured_ppb']:.4g} ppb  "
            f"predicted {row['predicted_ppb']:.4g} ppb  ratio {row['ratio']:.4g}"
        )

    record: Dict[str, Any] = {
        "generated_by": "scripts/generators/refit_sulfur_barriers_hofmann.py",
        "fit_target": path.name,
        "fit_target_doi": "10.1021/jf9705983",
        "objective": (
            "mean |log10(predicted_ppb / measured_ppb)| over the matched rows of "
            "cys_ribose_140C_Hofmann1998"
        ),
        "unscored_row_penalty_dex": UNSCORED_ROW_PENALTY,
        "forbidden_as_fit_targets": [
            "data/benchmarks/external_validation/** (hold-out; excluded by assertion)",
            "*Internal2026*  (synthetic reproducibility lane; excluded by assertion)",
            "*ProtocolPilot2026*  (synthetic reproducibility lane; excluded by assertion)",
            "spi_hvp_xylose_120C_PMC9905368 / wheat_gluten_hvp_xylose_120C_PMC9905368 "
            "(QUARANTINED 2026-08-27 as fabricated; they were already excluded here before "
            "that, because they are constrained by their own per-lane "
            "upstream_observability_factor)",
            "thiamine_cys_xylose_145C_Cerny2008 (VALUES_NEED_RE_EXAMINATION)",
            "thiamine_cys_glucose_120C_Bolton1994 (kept as an honest failure)",
        ],
        # 2026-08-27 (Wave I). Contamination review of this record, forced by the
        # quarantine of the two PMC9905368 benchmarks as fabricated.
        "contamination_review": {
            "date": "2026-08-27",
            "verdict": "THIS REFIT IS UNCONTAMINATED — its result STANDS.",
            "basis": (
                "The fit target is the single benchmark cys_ribose_140C_Hofmann1998 "
                "(10.1021/jf9705983, Hofmann & Schieberle, JAFC 46:235-241), a real, verified "
                "paper. The two quarantined PMC9905368 files were already on this script's "
                "forbidden list before the quarantine and contributed exactly zero rows to the "
                "objective, so the adopted value thiol_addition_norfuraneol 28.60 -> 26.85 and "
                "every 'keep the incumbent' decision here are unaffected."
            ),
            "what_does_NOT_stand": (
                "The separate hydrolysate-observability re-derivation "
                "(scripts/generators/rederive_hydrolysate_observability.py) IS retracted: its "
                "only targets were the two quarantined files. Its one applied value, the "
                "Methional base_factor, is reverted 0.05623 -> 0.0045. That layer sits "
                "DOWNSTREAM of these barriers and does not feed back into this objective."
            ),
            "SUPERSEDED_BY_A_LATER_WAVE_I_FIX": (
                "READ THIS BEFORE CITING ANY CONCLUSION FROM AN EARLIER VERSION OF THIS "
                "RECORD. Wave H's headline finding -- 'the profiles saturate; no barrier "
                "value in any defensible range reproduces the measured absolute yields; "
                "the residual is an ALLOCATION deficit, not a barrier deficit' -- was "
                "measured against a network in which the MFT route was throttled by a "
                "defect. The norfuraneol + H2S step consumed a pool reducing equivalent "
                "whose only source in a ribose/cysteine system was pyrazine chemistry, so "
                "MFT was structurally starved regardless of barrier value -- which is "
                "exactly what a saturating profile looks like. Wave I sourced that "
                "reductant from the thiol redox couple instead, and MFT went from 5.58x "
                "under to 1.45x under on this very benchmark WITH NO BARRIER CHANGED. "
                "The saturation was an artifact of the defect. `thiol_addition_norfuraneol` "
                "still ships at Wave H's 26.85, which is now a value fitted against a "
                "network that no longer exists. Re-running this refit is an OPEN OWNER "
                "ITEM and was deliberately NOT done inside Wave I: a refit here would be "
                "a recalibration event on top of a chemistry change, and the two must not "
                "be entangled in one pass."
            ),
            "standing_caveat_unrelated_to_the_quarantine": (
                "Hofmann1998's own validation contract (max_ratio 1.45 / mean_abs_log10_error "
                "0.09) is the 3rd-tightest in the collection -- the same fitting-tell pattern "
                "flagged for the quarantined files. The tolerance has NOT been widened, and "
                "the panel currently fails it, so nothing is being hidden; but a single "
                "benchmark with a suspiciously tight contract is a thin anchor for the whole "
                "sulfur branch, and that remains an open owner item."
            ),
        },
        "decision_rules": {
            "grid_step_kcal_per_mol": GRID_STEP_KCAL,
            "flat_profile_dex": FLAT_PROFILE_DEX,
            "min_objective_gain_dex": MIN_OBJECTIVE_GAIN_DEX,
            "indifference_dex": INDIFFERENCE_DEX,
            "conservative_edge": (
                "among values within indifference_dex of the profile minimum, adopt the "
                "LARGEST barrier (the least sulfur-favourable point the data cannot "
                "distinguish from the optimum)"
            ),
        },
        # What this refit ACTUALLY changed when it was first run, 2026-08-27. Kept here
        # because the script is idempotent: re-running it after the fit has been applied
        # reports "no further move" (which is the right convergence check, and is what the
        # committed artifact shows), and that would otherwise erase the record of the move.
        "applied_history": [
            {
                "date": "2026-08-27",
                "knob": "thiol_addition_norfuraneol",
                "from": 28.60,
                "to": 26.85,
                "objective_before": 0.6987,
                "objective_after": 0.6298,
                "hofmann_mft_fold_error_before": 7.83,
                "hofmann_mft_fold_error_after": 5.58,
                "hofmann_fft_fold_error_before": 3.19,
                "hofmann_fft_fold_error_after": 3.26,
                "basis": (
                    "conservative edge of the 23.30-26.85 indifference band over the "
                    "defensible range [23.30, 29.65]; the 28.60 it replaced had been "
                    "inherited from `thiol_addition` and from that key's PRE-Wave-G1 "
                    "Hofmann window, which did not survive the MFT route change"
                ),
            }
        ],
        "incumbent_values": incumbent,
        "incumbent_objective": baseline_objective,
        "incumbent_rows": baseline_rows,
        "reported_not_fitted": {},
        "profiles": {},
        "adopted_values": {},
    }

    # ---- reported-but-not-fitted derivatives -----------------------------
    for key, reason in REPORTED_NOT_FITTED.items():
        value = float(FAST_BARRIERS[key][0])
        probe: Dict[str, Any] = {"reason_not_fitted": reason, "incumbent": value, "probe": []}
        for delta in (-2.0, -1.0, 1.0, 2.0):
            objective, _rows = _evaluate_at(path, {key: value + delta})
            probe["probe"].append(
                {"value": round(value + delta, 4), "objective": objective}
            )
        span = max(item["objective"] for item in probe["probe"]) - min(
            item["objective"] for item in probe["probe"]
        )
        probe["objective_span_over_probe_dex"] = span
        probe["identifiable"] = span > FLAT_PROFILE_DEX
        record["reported_not_fitted"][key] = probe
        print(
            f"\n[reported, not fitted] {key} = {value:.2f}  "
            f"objective span over +/-2 kcal = {span:.4f} dex"
        )

    # ---- coordinate-wise profiles ---------------------------------------
    adopted: Dict[str, float] = dict(incumbent)
    for key, (low, high, basis) in KNOBS.items():
        grid = _frange(low, high, GRID_STEP_KCAL)
        profile = []
        for value in grid:
            objective, _rows = _evaluate_at(path, {key: value})
            profile.append({"value": value, "objective": objective})
        objectives = [point["objective"] for point in profile]
        best = min(objectives)
        worst = max(objectives)
        span = worst - best
        identifiable = span > FLAT_PROFILE_DEX
        gain = baseline_objective - best

        entry: Dict[str, Any] = {
            "range_kcal_per_mol": [low, high],
            "range_basis": basis,
            "incumbent": incumbent[key],
            "profile_min_objective": best,
            "profile_max_objective": worst,
            "profile_span_dex": span,
            "identifiable": identifiable,
            "achievable_gain_dex": gain,
            "profile": profile,
        }

        if not identifiable:
            entry["decision"] = "NOT IDENTIFIABLE — flat profile; incumbent kept"
        elif gain < MIN_OBJECTIVE_GAIN_DEX:
            entry["decision"] = (
                f"IMMATERIAL — best achievable gain {gain:.4f} dex < "
                f"{MIN_OBJECTIVE_GAIN_DEX} dex; incumbent kept"
            )
        else:
            indifferent = [
                point["value"]
                for point in profile
                if point["objective"] <= best + INDIFFERENCE_DEX
            ]
            chosen = max(indifferent)
            adopted[key] = chosen
            entry["indifference_band_kcal_per_mol"] = [min(indifferent), max(indifferent)]
            entry["adopted"] = chosen
            entry["decision"] = (
                f"MOVED {incumbent[key]:.2f} -> {chosen:.2f} "
                f"(conservative edge of the {min(indifferent):.2f}-{max(indifferent):.2f} "
                f"indifference band; gain {gain:.4f} dex)"
            )
        record["profiles"][key] = entry
        print(f"\n{key}: {entry['decision']}")
        print(
            f"   range [{low:.2f}, {high:.2f}]  profile span {span:.4f} dex  "
            f"min J {best:.4f}"
        )

    # ---- joint evaluation at the adopted point ---------------------------
    joint_objective, joint_rows = _evaluate_at(path, adopted)
    record["adopted_values"] = adopted
    record["objective_at_adopted"] = joint_objective
    record["rows_at_adopted"] = joint_rows
    record["objective_gain_dex"] = baseline_objective - joint_objective

    # 2026-08-27 (Wave I). Fit leverage — see src/fit_target_index.py and
    # scripts/ci/fit_target_gate.py. This one is the uncomfortable case and it is stated
    # rather than smoothed: four free barriers against a two-row benchmark.
    _fitted_rows = sum(1 for row in joint_rows if row.get("scored", True))
    _free = len(KNOBS)
    record["fit_leverage"] = {
        "free_parameters": _free,
        "fitted_rows": _fitted_rows,
        "parameters_per_row": (_free / _fitted_rows) if _fitted_rows else None,
        "class": "per_row_recovery",
        "interpretation": (
            f"{_free} free barriers against {_fitted_rows} rows. The fit has more freedom "
            "than the data constrains, so cys_ribose_140C_Hofmann1998 cannot be scored as "
            "out-of-sample evidence: its rows are excluded from the literature-coverage "
            "numerator AND denominator and reported in the fitted-row bucket instead."
        ),
        "why_this_exclusion_does_not_flatter_the_model": (
            "Excluding a fitted row removes its MISSES as well as its hits, and this row is "
            "a miss: even with four free barriers the fit cannot reproduce the measurement "
            "(MFT stays ~5.6x under, the profiles saturate inside their defensible ranges). "
            "The fitted-row bucket therefore has to be read, not skipped -- 'the model "
            "fails a benchmark it was fitted to' is a stronger negative result than any "
            "coverage number. See `irreducible_residual` below."
        ),
        "what_would_fix_the_leverage": (
            "A second independent literature benchmark of the sulfur branch. After the "
            "Mottram1994 / Farmer1999 quarantine there is exactly one, and it is this one."
        ),
    }

    # ---- what the fit cannot fix -----------------------------------------
    floor_point = {key: KNOBS[key][0] for key in KNOBS}
    floor_objective, floor_rows = _evaluate_at(path, floor_point)
    record["irreducible_residual"] = {
        "definition": (
            "objective with EVERY fitted knob pinned at the bottom of its defensible "
            "range — the best the sulfur branch can do on this benchmark by barrier "
            "moves alone"
        ),
        "values": floor_point,
        "objective": floor_objective,
        "rows": floor_rows,
        "reading": (
            "The profiles saturate well inside the defensible ranges: the MFT ceiling is "
            "set by the shared upstream span (Schiff/Amadori) and by the volatile-budget "
            "allocation, in which furfural — unmeasured in this benchmark — takes ~78% of "
            "a total budget that is itself the right order of magnitude (~1050 ppb against "
            "a measured MFT+FFT of 542 ppb). No barrier value in any defensible range "
            "reproduces the measured absolute sulfur yields; the remaining factor of ~5 on "
            "MFT and ~3 on FFT is an allocation deficit, not a barrier deficit. That is the "
            "finding this refit exists to produce, and it is reported rather than tuned away."
        ),
    }

    print()
    print(f"Objective at adopted point: {joint_objective:.4f} dex "
          f"(incumbent {baseline_objective:.4f})")
    for row in joint_rows:
        print(
            f"   {row['compound']:26s} predicted {row['predicted_ppb']:.4g} ppb  "
            f"ratio {row['ratio']:.4g}  ({1.0 / row['ratio']:.2f}x under)"
        )
    print(f"Floor of the defensible ranges: {floor_objective:.4f} dex")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "sulfur_barrier_refit_hofmann.json"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    lines = [
        "# Sulfur-branch barrier refit against Hofmann1998 (Wave H)",
        "",
        f"Fit target: `{path.name}` (DOI 10.1021/jf9705983) — the ONLY surviving "
        "literature constraint on the sulfur branch after the Mottram1994 / Farmer1999 "
        "quarantine.",
        "",
        "> **CORRECTED 2026-08-27 (Wave S2c): the line above is false in the word "
        "\"literature\".** Wave S2b traced this benchmark's MFT 342 / FFT 200 ppb to "
        "`data/benchmarks/maillard_validation_benchmarks.md` §1.3, an abstract-reconstructed "
        "range table committed in the same commit as the benchmark JSON; both values are "
        "interior points of two invented, overlapping mol % bands. **The sulfur branch has "
        "zero absolute literature anchors.** Do not re-run this script against this "
        "benchmark.",
        "",
        f"Objective: `{record['objective']}`",
        "",
        "Synthetic (Internal2026 / ProtocolPilot2026) lanes and the "
        "`external_validation/` hold-out are excluded **by assertion**, not convention. "
        "The xylose HVP lanes are excluded because they are constrained by their own "
        "per-lane `upstream_observability_factor`, re-derived separately.",
        "",
        "## Decision rules",
        "",
        f"1. **Identifiability** — a profile flat to {FLAT_PROFILE_DEX} dex over the whole "
        "defensible range keeps its incumbent.",
        f"2. **Materiality** — a knob moves only for a gain of at least "
        f"{MIN_OBJECTIVE_GAIN_DEX} dex.",
        f"3. **Conservative edge** — among values within {INDIFFERENCE_DEX} dex of the "
        "profile minimum, the LARGEST barrier is adopted.",
        "",
        "## What this refit changed (2026-08-27)",
        "",
        "| Knob | From | To | Objective | Hofmann MFT | Hofmann FFT |",
        "| --- | --- | --- | --- | --- | --- |",
    ] + [
        f"| `{item['knob']}` | {item['from']:.2f} | {item['to']:.2f} | "
        f"{item['objective_before']:.4f} -> {item['objective_after']:.4f} dex | "
        f"{item['hofmann_mft_fold_error_before']:.2f}x -> "
        f"{item['hofmann_mft_fold_error_after']:.2f}x under | "
        f"{item['hofmann_fft_fold_error_before']:.2f}x -> "
        f"{item['hofmann_fft_fold_error_after']:.2f}x under |"
        for item in record["applied_history"]
    ] + [
        "",
        "This script is idempotent: re-run after the fit is applied and it reports no "
        "further move, which is the convergence check. The section below is that re-run.",
        "",
        "## Result",
        "",
        "| Knob | Range | Basis | Incumbent | Profile span (dex) | Decision |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for key, entry in record["profiles"].items():
        low, high = entry["range_kcal_per_mol"]
        lines.append(
            f"| `{key}` | [{low:.2f}, {high:.2f}] | {entry['range_basis']} | "
            f"{entry['incumbent']:.2f} | {entry['profile_span_dex']:.4f} | "
            f"{entry['decision']} |"
        )
    for key, entry in record["reported_not_fitted"].items():
        lines.append(
            f"| `{key}` | reported, not fitted | {entry['reason_not_fitted']} | "
            f"{entry['incumbent']:.2f} | {entry['objective_span_over_probe_dex']:.4f} | "
            "NOT FITTED |"
        )
    lines += [
        "",
        f"Objective: **{baseline_objective:.4f} dex** (incumbent) -> "
        f"**{joint_objective:.4f} dex** (adopted).",
        "",
        "| Compound | Measured ppb | Predicted ppb (incumbent) | Predicted ppb (adopted) | Fold error (adopted) |",
        "| --- | --- | --- | --- | --- |",
    ]
    by_compound = {row["compound"]: row for row in baseline_rows}
    for row in joint_rows:
        before = by_compound.get(row["compound"], {})
        fold = 1.0 / row["ratio"] if row["ratio"] else float("nan")
        lines.append(
            f"| {row['compound']} | {row['measured_ppb']:.4g} | "
            f"{before.get('predicted_ppb', float('nan')):.4g} | "
            f"{row['predicted_ppb']:.4g} | {fold:.2f}x under |"
        )
    lines += [
        "",
        "## What this fit cannot fix",
        "",
        record["irreducible_residual"]["reading"],
        "",
        f"With every fitted knob pinned at the bottom of its defensible range the "
        f"objective is still **{floor_objective:.4f} dex**.",
        "",
    ]
    # 2026-08-27 (Wave I): record whether the adopted point is actually SHIPPED. Without
    # this the record reads as if `adopted_values` were live, which they are not when the
    # script is run in report-only mode -- and the Wave I re-run WAS report-only, on
    # purpose (see `applied_to_runtime` below).
    shipped_now = {key: float(FAST_BARRIERS[key][0]) for key in KNOBS}
    would_move = {
        key: {"shipped": shipped_now[key], "adopted": adopted[key]}
        for key in KNOBS
        if abs(adopted[key] - shipped_now[key]) > 1e-9
    }
    record["applied_to_runtime"] = bool(args.apply and would_move)
    record["shipped_constants"] = shipped_now
    if would_move and not args.apply:
        record["adopted_but_NOT_APPLIED"] = {
            "constants": would_move,
            "objective_shipped_dex": baseline_objective,
            "objective_if_applied_dex": joint_objective,
            "why_not_applied": (
                "Wave I (2026-08-27) re-ran this record because the network changed "
                "underneath it: decoupling the MFT step from the pyrazine-supplied hydrogen "
                "pool took cys_ribose_140C_Hofmann1998 from 5.58x under to 1.45x under with "
                "NO barrier changed. Applying a barrier refit on top of a chemistry change, "
                "in the same pass, would entangle the two and make it impossible to say "
                "afterwards which one produced the agreement. The record is written; the "
                "constants are not moved. Applying it is an OPEN OWNER ITEM."
            ),
        }
    json_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    if would_move and not args.apply:
        lines.extend([
            "",
            "## NOT APPLIED",
            "",
            "**The adopted point above is NOT what ships.** This run was report-only. Shipped "
            "vs adopted:",
            "",
            "| Constant | Shipped | Adopted here |",
            "| --- | ---: | ---: |",
        ])
        for key, pair in sorted(would_move.items()):
            lines.append(f"| `{key}` | {pair['shipped']:.2f} | {pair['adopted']:.2f} |")
        lines.extend([
            "",
            record["adopted_but_NOT_APPLIED"]["why_not_applied"],
        ])

    md_path = output_dir / "sulfur_barrier_refit_hofmann.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")

    if args.apply:
        moved = {
            key: value
            for key, value in adopted.items()
            if abs(value - incumbent[key]) > 1e-9
        }
        if not moved:
            print("No constant moved; src/barrier_constants.py left untouched.")
            return 0
        text = BARRIER_SOURCE.read_text(encoding="utf-8")
        for key, value in moved.items():
            text = _rewrite_constant(text, key, value, FAST_BARRIERS[key][1])
        BARRIER_SOURCE.write_text(text, encoding="utf-8")
        print(f"Applied to {BARRIER_SOURCE}: " + ", ".join(
            f"{key} {incumbent[key]:.2f} -> {value:.2f}" for key, value in moved.items()
        ))
        print("Update each moved constant's inline comment by hand — the value and its "
              "stated provenance must stay in sync.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
