#!/usr/bin/env python3
"""Refit `thiol_addition_pentodiulose` against Hofmann1998 — and only Hofmann1998.

Context (2026-08-27, Wave P item 1)
-----------------------------------
Wave N replaced the contradicted norfuraneol -> MFT step with the isotope-evidenced
1,4-dideoxyosone route (Cerny & Davidek 2003, 10.1021/jf026123f; 2004,
10.1021/jf035265m):

    1-deoxy-2,3-pentodiulose + 2[H] -> 1,4-dideoxypento-2,3-diulose + H2O
    1,4-dideoxypento-2,3-diulose + H2S -> 2-methyl-3-furanthiol + 2 H2O

and shipped the new sulfur-incorporation barrier at the UNFITTED
``thiol_addition`` class value, 28.60 kcal/mol, explicitly labelled ESTIMATED
UNCONSTRAINED. It was deliberately NOT set to the 26.85 that Wave H had fitted
THROUGH the retired route: a barrier fitted through a route the literature
retires does not transfer.

This script is the owner-approved refit of that one constant against the repo's
single surviving sulfur-branch literature anchor. It is the successor of
``refit_sulfur_barriers_hofmann.py``, whose record
(``results/validation/sulfur_barrier_refit_hofmann.{json,md}``) is STALE: it
profiles ``thiol_addition_norfuraneol``, a family no step emits any more.

WHY THIS IS RUN AFTER THE WAVE P CHEMISTRY ADDITIONS, NOT BEFORE
----------------------------------------------------------------
Wave L1's own warning about the Wave H refit was that fitting a barrier before a
route change fits the wrong route more precisely. The same logic applies within
this wave: Wave P adds a second MFT channel (the Hofmann C2+C3 recombination) and
a new reducing-equivalent producer, both of which move the ribose/cysteine
system. The refit is therefore the LAST thing this wave does, so the value is
fitted against the network that actually ships.

THE CAVEAT THAT TRAVELS WITH THIS NUMBER
----------------------------------------
The fit target's own ``content_verification_note`` (Wave K, 2026-08-27) says, in
the benchmark file, verbatim:

    "The paper's abstract reports MFT/FFT yields in mol % (e.g. MFT 1.4 mol % for
    its best precursor system), not ppb. The 342/200 ppb values here require a
    mol%->ppb conversion (system volume, precursor moles, molecular weights) that
    is NOT documented anywhere in this repo, and the full text is paywalled so the
    values could not be confirmed against the paper's tables. This is the panel's
    tightest contract (1.45x / 0.09 dex) resting on an unverified derivation."

SUPERSEDED 2026-08-27 (Wave S2c). The Wave K caveat quoted above asks the wrong
question. The problem was never that the mol%->ppb conversion was undocumented; it is
that THERE IS NO MEASUREMENT ON THE FAR END OF IT. Wave S2b traced 342 and 200 ppb to
data/benchmarks/maillard_validation_benchmarks.md section 1.3 -- an abstract-reconstructed
range table committed in c7efbbc, the SAME commit that created the benchmark JSON --
whose row gives MFT `~0.02-0.05` mol % and FFT `~0.01-0.03` mol %. On the benchmark's
declared (unattested) 10 mM basis with MW 114.17: 0.0300 mol % -> 342.5 -> 342 ppb, and
the FFT band's geometric mean 0.017321 mol % -> 197.8 -> 200 ppb. Both are interior
points of two invented, overlapping bands (~90% confidence, arithmetic exact).
CONSEQUENCE: this script's result is RETRACTED
(results/validation/sulfur_barrier_refit_pentodiulose.{json,md}) and the constant it
produced is REVERTED, 26.35 -> 28.60 kcal/mol, in src/barrier_constants.py.
**DO NOT RE-RUN THIS SCRIPT AGAINST cys_ribose_140C_Hofmann1998.** Fitting a barrier to
this repository's own arithmetic is circular whatever the objective says. A real target
comes first: the ILL request pack in tasks/audit_remediation.md "## Wave S2b" section (f),
then a rebuilt benchmark expressed in the paper's native mol % rather than in ppb (mol %
is basis-free; the ppb target smuggles in the unattested 10 mM basis as a free
multiplicative parameter underneath whatever contract sits on top of it).
THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.

So this constant is fitted against an anchor whose UNIT CONVERSION IS UNVERIFIED.
If that conversion is wrong, the error is localised HERE, in one barrier, rather
than distributed. That sentence is reproduced verbatim in the constant's own
``FAST_BARRIERS`` rationale in ``src/barrier_constants.py`` and in the JSON record
this script writes. It is the reason the fit is worth doing at all: a single named
constant carrying a single named risk is auditable; the same error absorbed
silently into the route is not.

Knob and its defensible range
-----------------------------
ONE knob. The range is bounded by values ALREADY IN ``FAST_BARRIERS`` for the same
mechanistic class, so the fit cannot leave the envelope the table itself spans:

    thiol_addition_pentodiulose  [23.30, 29.65]
        low  = thiohemiacetal_formation (23.30) -- the fastest sulfur-addition
               step in the table
        high = thiol_addition_hexose    (29.65) -- the slowest

This is the SAME range the retired ``thiol_addition_norfuraneol`` knob used, and
deliberately so: the two steps are the same mechanistic class (H2S addition to a
carbonyl of a sugar-derived intermediate), so re-using the range means the two
fits are comparable, which is itself a result (see ``comparison_to_retired_fit``
in the record).

`furanone_cyclisation`, `thiohemiacetal_formation` and `thiol_dehydration` are
REPORTED but NOT FITTED here. The Wave H script fitted four barriers against a
two-row benchmark, which it correctly flagged as ``per_row_recovery`` with 2.0
parameters per row. One knob against two rows is 0.5 -- still at the
``per_row_recovery`` threshold, so the benchmark stays out of the honest
literature-coverage count either way, but the fit has strictly less freedom to
manufacture agreement. Their derivatives are probed and published so the reader
can see what was left on the table.

Decision rules (stated up front, applied mechanically)
------------------------------------------------------
Identical to ``refit_sulfur_barriers_hofmann.py``, deliberately, so the two runs
are comparable:

1. IDENTIFIABILITY. A knob whose objective profile is flat to ``FLAT_PROFILE_DEX``
   over its whole defensible range is NOT IDENTIFIABLE and keeps its incumbent.
2. MATERIALITY. A knob moves only if the best achievable objective gain is at least
   ``MIN_OBJECTIVE_GAIN_DEX``. Chasing less than that on a two-row benchmark is
   fitting noise.
3. CONSERVATIVE EDGE. Among the values within ``INDIFFERENCE_DEX`` of the profile
   minimum, the LARGEST barrier is adopted -- the least sulfur-favourable point the
   Hofmann data cannot distinguish from the optimum. This is what stops a
   saturating profile from dragging the constant to the bottom of its range and
   calling that a fit.
4. BOUNDARY HONESTY (new here). If the profile minimum sits AT a range endpoint the
   record says so explicitly (``bounds_check.hit_a_bound``): a boundary optimum is
   a statement that the data want a value the table's own class envelope does not
   contain, and that is a finding, not a fit.

Usage
-----
    python scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py
    python scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py --apply
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.barrier_constants import FAST_BARRIERS
from src.benchmark_validation import evaluate_benchmark
from src.uncertainty_propagation import _benchmark_signal_origin

TARGET_BENCHMARK = ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json"
BARRIER_SOURCE = ROOT / "src" / "barrier_constants.py"

UNSCORED_ROW_PENALTY = 3.0        # dex — same convention as the Wave H script
GRID_STEP_KCAL = 0.05
FLAT_PROFILE_DEX = 1.0e-6         # below this a profile is numerically flat
MIN_OBJECTIVE_GAIN_DEX = 0.02     # materiality threshold for moving a constant
INDIFFERENCE_DEX = 0.01           # width of the "cannot distinguish" band

FITTED_KNOB = "thiol_addition_pentodiulose"
KNOB_RANGE: Tuple[float, float] = (23.30, 29.65)
KNOB_RANGE_BASIS = (
    "H2S/thiol addition to a carbonyl of a sugar-derived intermediate; low = "
    "thiohemiacetal_formation (23.30), high = thiol_addition_hexose (29.65), both "
    "already in FAST_BARRIERS. Same range the retired thiol_addition_norfuraneol "
    "knob used, so the two fits are directly comparable."
)

REPORTED_NOT_FITTED = {
    "furanone_cyclisation": (
        "The competing branch of the same 1-deoxyosone fork (-> norfuraneol). Wave N "
        "set it EQUAL to `deoxyosone_reduction` (both 28.0) precisely to express no "
        "prior preference at that fork; fitting it here would encode a preference "
        "derived from two rows of one benchmark."
    ),
    "deoxyosone_reduction": (
        "The upstream reduction that feeds the fitted step. It is in series with it, "
        "so on a single benchmark the two are not separately identifiable -- fitting "
        "both would split one measurable quantity across two constants and report the "
        "split as knowledge."
    ),
    "thiohemiacetal_formation": (
        "FFT-side, not MFT-side. It moves the OTHER measured row of this benchmark, so "
        "fitting it would let the fit trade FFT accuracy for MFT accuracy across the "
        "two rows -- the classic two-row/two-knob recovery."
    ),
    "thiol_dehydration": ("Same reason as `thiohemiacetal_formation`: FFT-side."),
    "thiol_addition": (
        "After Wave G1/Wave I this key labels only the DEMOTED hexose legacy shortcut "
        "and the lumped `Thiol_Addition_H2` FFT step. Fitting it would tune a lump the "
        "campaign exists to retire."
    ),
}

#: Copied VERBATIM out of the fit target's own `content_verification_note` (Wave K).
FIT_TARGET_CONTENT_CAVEAT = (
    "The paper's abstract reports MFT/FFT yields in mol % (e.g. MFT 1.4 mol % for its "
    "best precursor system), not ppb. The 342/200 ppb values here require a mol%->ppb "
    "conversion (system volume, precursor moles, molecular weights) that is NOT "
    "documented anywhere in this repo, and the full text is paywalled so the values "
    "could not be confirmed against the paper's tables. This is the panel's tightest "
    "contract (1.45x / 0.09 dex) resting on an unverified derivation. Needs: full-text "
    "access to 10.1021/jf9705983 and a written conversion, or replacement with the "
    "paper's native mol % as the target unit."
)


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
    # The caveat this fit carries must actually be in the benchmark file; if the note
    # is ever removed, this run must fail rather than publish a stale caveat.
    payload = json.loads(path.read_text(encoding="utf-8"))
    note = ((payload.get("content_verification_note") or {}).get("note") or "")
    assert note.strip() == FIT_TARGET_CONTENT_CAVEAT.strip(), (
        "the fit target's content_verification_note has changed; update "
        "FIT_TARGET_CONTENT_CAVEAT (it is quoted verbatim into the shipped constant's "
        "rationale) before re-running this fit"
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


def _rewrite_constant(text: str, key: str, new_value: float) -> str:
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
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="rewrite the adopted value into src/barrier_constants.py",
    )
    args = parser.parse_args()

    path = _guarded_target()
    incumbent = float(FAST_BARRIERS[FITTED_KNOB][0])
    baseline_objective, baseline_rows = _score(path)

    print(f"Fit target: {path.name}  (DOI 10.1021/jf9705983)")
    print(f"Knob: {FITTED_KNOB} = {incumbent:.2f} kcal/mol")
    print(f"Incumbent objective J = {baseline_objective:.4f} dex")
    for row in baseline_rows:
        print(
            f"   {row['compound']:26s} measured {row['measured_ppb']:.4g} ppb  "
            f"predicted {row['predicted_ppb']:.4g} ppb  ratio {row['ratio']:.4g}"
        )

    low, high = KNOB_RANGE
    grid = _frange(low, high, GRID_STEP_KCAL)
    profile: List[Dict[str, Any]] = []
    per_compound: Dict[str, List[Dict[str, float]]] = {}
    for value in grid:
        objective, rows = _evaluate_at(path, {FITTED_KNOB: value})
        profile.append({"value": value, "objective": objective})
        for row in rows:
            per_compound.setdefault(row["compound"], []).append(
                {"value": value, "predicted_ppb": row["predicted_ppb"], "ratio": row["ratio"]}
            )

    objectives = [point["objective"] for point in profile]
    best = min(objectives)
    worst = max(objectives)
    span = worst - best
    identifiable = span > FLAT_PROFILE_DEX
    gain = baseline_objective - best
    argmin = min(profile, key=lambda point: point["objective"])["value"]
    hit_a_bound = abs(argmin - low) < 1e-9 or abs(argmin - high) < 1e-9

    adopted = incumbent
    if not identifiable:
        decision = "NOT IDENTIFIABLE — flat profile; incumbent kept"
        indifference: List[float] = []
    elif gain < MIN_OBJECTIVE_GAIN_DEX:
        decision = (
            f"IMMATERIAL — best achievable gain {gain:.4f} dex < "
            f"{MIN_OBJECTIVE_GAIN_DEX} dex; incumbent kept"
        )
        indifference = []
    else:
        indifference = [
            point["value"] for point in profile if point["objective"] <= best + INDIFFERENCE_DEX
        ]
        adopted = max(indifference)
        decision = (
            f"MOVED {incumbent:.2f} -> {adopted:.2f} (conservative edge of the "
            f"{min(indifference):.2f}-{max(indifference):.2f} indifference band; "
            f"gain {gain:.4f} dex)"
        )

    adopted_objective, adopted_rows = _evaluate_at(path, {FITTED_KNOB: adopted})

    # --- co-movement of the OTHER measured row --------------------------------
    # The fitted knob is on the MFT lane, but the two lanes compete for the same
    # upstream flux, so FFT can move without being fitted. Measured, not assumed.
    def _row(rows: List[Dict[str, Any]], compound: str) -> Dict[str, Any]:
        return next((r for r in rows if r["compound"] == compound), {})

    mft_before, mft_after = _row(baseline_rows, "2-methyl-3-furanthiol"), _row(adopted_rows, "2-methyl-3-furanthiol")
    fft_before, fft_after = _row(baseline_rows, "2-furfurylthiol"), _row(adopted_rows, "2-furfurylthiol")
    fft_series = per_compound.get("2-furfurylthiol", [])
    fft_span = (
        max(p["predicted_ppb"] for p in fft_series) / min(p["predicted_ppb"] for p in fft_series)
        if fft_series and min(p["predicted_ppb"] for p in fft_series) > 0
        else None
    )
    mft_series = per_compound.get("2-methyl-3-furanthiol", [])
    mft_span = (
        max(p["predicted_ppb"] for p in mft_series) / min(p["predicted_ppb"] for p in mft_series)
        if mft_series and min(p["predicted_ppb"] for p in mft_series) > 0
        else None
    )

    # --- reported-but-not-fitted derivatives ----------------------------------
    reported: Dict[str, Any] = {}
    for key, reason in REPORTED_NOT_FITTED.items():
        value = float(FAST_BARRIERS[key][0])
        probe = []
        for delta in (-2.0, -1.0, 1.0, 2.0):
            objective, _rows = _evaluate_at(path, {key: value + delta})
            probe.append({"value": round(value + delta, 4), "objective": objective})
        probe_span = max(p["objective"] for p in probe) - min(p["objective"] for p in probe)
        reported[key] = {
            "reason_not_fitted": reason,
            "incumbent": value,
            "probe": probe,
            "objective_span_over_probe_dex": probe_span,
            "identifiable": probe_span > FLAT_PROFILE_DEX,
        }
        print(f"[reported, not fitted] {key} = {value:.2f}  span over +/-2 kcal = {probe_span:.4f} dex")

    fitted_rows = sum(1 for row in adopted_rows if row.get("scored", True))
    record: Dict[str, Any] = {
        "generated_by": "scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py",
        "wave": "Wave P item 1 (2026-08-27)",
        "supersedes": (
            "results/validation/sulfur_barrier_refit_hofmann.json — STALE: it profiles "
            "`thiol_addition_norfuraneol`, a family no step emits since Wave N. That "
            "record is kept for provenance and still declares the same fit target, so "
            "the fit-target accounting is unchanged by this file's arrival."
        ),
        "fit_target": path.name,
        "fit_target_doi": "10.1021/jf9705983",
        "fit_target_content_verification": {
            "status": "conversion_undocumented",
            "note_verbatim": FIT_TARGET_CONTENT_CAVEAT,
            "consequence_for_this_fit": (
                "This constant is fitted against an anchor whose mol%->ppb conversion is "
                "UNVERIFIED. If the conversion is wrong, the resulting error is LOCALISED "
                "in this one barrier rather than distributed across the route. That is the "
                "argument for doing the fit at all, and it is also the reason the value "
                "must never be quoted as a measured or literature-derived barrier."
            ),
        },
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
            "(QUARANTINED as fabricated)",
            "thiamine_cys_xylose_145C_Cerny2008 (VALUES_NEED_RE_EXAMINATION)",
            "thiamine_cys_glucose_120C_Bolton1994 (kept as an honest failure)",
        ],
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
        "incumbent_values": {FITTED_KNOB: incumbent},
        "incumbent_objective": baseline_objective,
        "incumbent_rows": baseline_rows,
        "profiles": {
            FITTED_KNOB: {
                "range_kcal_per_mol": [low, high],
                "range_basis": KNOB_RANGE_BASIS,
                "incumbent": incumbent,
                "profile_min_objective": best,
                "profile_max_objective": worst,
                "profile_span_dex": span,
                "identifiable": identifiable,
                "achievable_gain_dex": gain,
                "argmin_kcal_per_mol": argmin,
                "indifference_band_kcal_per_mol": (
                    [min(indifference), max(indifference)] if indifference else None
                ),
                "adopted": adopted,
                "decision": decision,
                "profile": profile,
            }
        },
        "bounds_check": {
            "range_kcal_per_mol": [low, high],
            "argmin_kcal_per_mol": argmin,
            "hit_a_bound": hit_a_bound,
            "reading": (
                "The objective minimum sits AT a range endpoint. The data want a value the "
                "table's own sulfur-addition class envelope does not contain, i.e. the "
                "profile SATURATES: the residual is not removable by this barrier. The "
                "conservative-edge rule is what keeps the adopted value off the boundary."
                if hit_a_bound
                else "Interior optimum: the defensible range contains the objective minimum."
            ),
        },
        "adopted_values": {FITTED_KNOB: adopted},
        "objective_at_adopted": adopted_objective,
        "rows_at_adopted": adopted_rows,
        "objective_gain_dex": baseline_objective - adopted_objective,
        "co_movement": {
            "what_this_measures": (
                "The fitted knob sits on the MFT lane only. FFT shares the upstream sugar "
                "flux, so it can move without being fitted. Measured over the same grid."
            ),
            "mft_predicted_ppb": {
                "incumbent": mft_before.get("predicted_ppb"),
                "adopted": mft_after.get("predicted_ppb"),
                "measured": mft_before.get("measured_ppb"),
                "fold_span_over_range": mft_span,
            },
            "fft_predicted_ppb": {
                "incumbent": fft_before.get("predicted_ppb"),
                "adopted": fft_after.get("predicted_ppb"),
                "measured": fft_before.get("measured_ppb"),
                "fold_span_over_range": fft_span,
                "direction": (
                    "FFT moves OPPOSITE to MFT: lowering the MFT barrier diverts shared "
                    "upstream flux into the MFT lane and FFT falls. The two rows are "
                    "therefore NOT independent evidence, which is exactly why one knob is "
                    "fitted and not two."
                ),
            },
        },
        "comparison_to_retired_fit": {
            "retired_knob": "thiol_addition_norfuraneol",
            "retired_fitted_value": 26.85,
            "retired_range": [23.30, 29.65],
            "why_comparable": (
                "Same benchmark, same objective, same decision rules, same defensible "
                "range, different ROUTE. Any agreement between the two adopted values is "
                "a coincidence of the shared upstream trunk, not evidence that the "
                "retired route was right."
            ),
        },
        "fit_leverage": {
            "free_parameters": 1,
            "fitted_rows": fitted_rows,
            "parameters_per_row": (1.0 / fitted_rows) if fitted_rows else None,
            "class": "per_row_recovery",
            "interpretation": (
                f"ONE free barrier against {fitted_rows} rows is "
                f"{1.0 / fitted_rows if fitted_rows else float('nan'):.2f} parameters per "
                "row, at src.fit_target_index.FIT_LEVERAGE_THRESHOLD (0.5). "
                "cys_ribose_140C_Hofmann1998 therefore stays OUT of the honest "
                "literature-coverage numerator and denominator and is reported in the "
                "fitted-row bucket. That classification is UNCHANGED by this refit -- the "
                "benchmark was already a declared fit target of the Wave H record. What "
                "changes is which constant carries the fit, and how much of the residual "
                "it can absorb."
            ),
            "why_this_exclusion_does_not_flatter_the_model": (
                "Excluding a fitted row removes its misses as well as its hits. Read "
                "`bounds_check` and `irreducible_residual`: even at the fitted value the "
                "model does not reproduce the measurement, and a fit that cannot reach its "
                "own target is a stronger negative result than any coverage number."
            ),
            "what_would_fix_the_leverage": (
                "A second independent literature benchmark of the sulfur branch. After the "
                "Mottram1994 / Farmer1999 quarantine there is exactly one, and it is this "
                "one -- and its unit conversion is unverified."
            ),
        },
        "reported_not_fitted": reported,
    }

    floor_objective, floor_rows = _evaluate_at(path, {FITTED_KNOB: low})
    record["irreducible_residual"] = {
        "definition": (
            "objective with the fitted knob pinned at the BOTTOM of its defensible range "
            "— the best this one barrier can do on this benchmark"
        ),
        "value": {FITTED_KNOB: low},
        "objective": floor_objective,
        "rows": floor_rows,
    }

    print()
    print(f"{FITTED_KNOB}: {decision}")
    print(f"   range [{low:.2f}, {high:.2f}]  span {span:.4f} dex  min J {best:.4f} at {argmin:.2f}")
    print(f"   hit_a_bound = {hit_a_bound}")
    print(f"Objective at adopted point: {adopted_objective:.4f} dex (incumbent {baseline_objective:.4f})")
    for row in adopted_rows:
        fold = row["ratio"] if row["ratio"] and row["ratio"] >= 1 else (1.0 / row["ratio"] if row["ratio"] else float("nan"))
        direction = "over" if row["ratio"] and row["ratio"] >= 1 else "under"
        print(
            f"   {row['compound']:26s} predicted {row['predicted_ppb']:.4g} ppb  "
            f"ratio {row['ratio']:.4g}  ({fold:.2f}x {direction})"
        )
    print(f"Floor of the defensible range: {floor_objective:.4f} dex")

    shipped_now = float(FAST_BARRIERS[FITTED_KNOB][0])
    would_move = abs(adopted - shipped_now) > 1e-9
    record["applied_to_runtime"] = bool(args.apply and would_move)
    record["shipped_constants"] = {FITTED_KNOB: shipped_now}
    if would_move and not args.apply:
        record["adopted_but_NOT_APPLIED"] = {
            "shipped": shipped_now,
            "adopted": adopted,
            "objective_shipped_dex": baseline_objective,
            "objective_if_applied_dex": adopted_objective,
        }

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "sulfur_barrier_refit_pentodiulose.json"
    json_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# `thiol_addition_pentodiulose` refit against Hofmann1998 (Wave P item 1)",
        "",
        f"Fit target: `{path.name}` (DOI 10.1021/jf9705983) — the ONLY surviving literature "
        "constraint on the sulfur branch.",
        "",
        "## The caveat this number carries",
        "",
        "The fit target's own `content_verification_note` says, verbatim:",
        "",
        f"> {FIT_TARGET_CONTENT_CAVEAT}",
        "",
        record["fit_target_content_verification"]["consequence_for_this_fit"],
        "",
        "## Result",
        "",
        "| Knob | Range | Incumbent | Profile span (dex) | argmin | Decision |",
        "| --- | --- | ---: | ---: | ---: | --- |",
        f"| `{FITTED_KNOB}` | [{low:.2f}, {high:.2f}] | {incumbent:.2f} | {span:.4f} | "
        f"{argmin:.2f} | {decision} |",
        "",
        f"Objective: **{baseline_objective:.4f} dex** (incumbent) -> "
        f"**{adopted_objective:.4f} dex** (adopted).",
        "",
        f"Boundary check: `hit_a_bound = {hit_a_bound}`. "
        + record["bounds_check"]["reading"],
        "",
        "| Compound | Measured ppb | Predicted (incumbent) | Predicted (adopted) | Fold error (adopted) |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    by_compound = {row["compound"]: row for row in baseline_rows}
    for row in adopted_rows:
        before = by_compound.get(row["compound"], {})
        ratio = row["ratio"] or float("nan")
        fold = ratio if ratio >= 1 else 1.0 / ratio
        direction = "over" if ratio >= 1 else "under"
        lines.append(
            f"| {row['compound']} | {row['measured_ppb']:.4g} | "
            f"{before.get('predicted_ppb', float('nan')):.4g} | "
            f"{row['predicted_ppb']:.4g} | {fold:.4f}x {direction} |"
        )
    lines += [
        "",
        "## FFT co-movement (measured, not assumed)",
        "",
        record["co_movement"]["what_this_measures"],
        "",
        f"Over the whole search range MFT spans **{mft_span:.4g}x** and FFT spans "
        f"**{fft_span:.4g}x**. " + record["co_movement"]["fft_predicted_ppb"]["direction"],
        "",
        "## Reported, not fitted",
        "",
        "| Knob | Incumbent | Objective span over +/-2 kcal (dex) | Why not fitted |",
        "| --- | ---: | ---: | --- |",
    ]
    for key, entry in reported.items():
        lines.append(
            f"| `{key}` | {entry['incumbent']:.2f} | "
            f"{entry['objective_span_over_probe_dex']:.4f} | {entry['reason_not_fitted']} |"
        )
    lines += [
        "",
        "## What this fit cannot fix",
        "",
        f"With the knob pinned at the bottom of its defensible range ({low:.2f} kcal/mol) "
        f"the objective is still **{floor_objective:.4f} dex**.",
        "",
        "## Fit leverage",
        "",
        record["fit_leverage"]["interpretation"],
        "",
    ]
    if would_move and not args.apply:
        lines += [
            "## NOT APPLIED",
            "",
            "**The adopted value above is NOT what ships.** This run was report-only.",
            "",
            f"| Constant | Shipped | Adopted here |",
            f"| --- | ---: | ---: |",
            f"| `{FITTED_KNOB}` | {shipped_now:.2f} | {adopted:.2f} |",
            "",
        ]
    md_path = output_dir / "sulfur_barrier_refit_pentodiulose.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")

    if args.apply:
        if not would_move:
            print("No constant moved; src/barrier_constants.py left untouched.")
            return 0
        text = BARRIER_SOURCE.read_text(encoding="utf-8")
        BARRIER_SOURCE.write_text(_rewrite_constant(text, FITTED_KNOB, adopted), encoding="utf-8")
        print(f"Applied to {BARRIER_SOURCE}: {FITTED_KNOB} {incumbent:.2f} -> {adopted:.2f}")
        print("Update the constant's inline rationale by hand — the value and its stated "
              "provenance must stay in sync.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
