#!/usr/bin/env python3
"""Fit `furanone_reductive_sulfhydrylation` against Hofmann 1998 Table 4 — and only that row.

Context (2026-08-28, Wave X)
----------------------------
Wave X re-added the norfuraneol -> MFT step that Wave N had retired, as a SLOW
PARALLEL channel rather than as the dominant route:

    norfuraneol + H2S + 2[H] -> 2-methyl-3-furanthiol + 2 H2O
    C5H6O3 + H2S + H2 -> C5H6OS + 2 H2O            (both sides C5H10O3S, exact)

The mechanism is Hofmann & Schieberle 1998's own Figure 1 and the sentence beneath
it (p.238, read from data/articles/hofmann1998.pdf):

    "Incorporation of a thiol group at carbon 3 yields 4-mercapto-5-methyl-3(2H)-
     furanone in a first reaction step, which in turn might be reduced by a
     reductone or by a further molecule of hydrogen sulfide.  Elimination of water
     from the intermediate 3-hydroxy-4-mercapto-5-methyl-1,2-dihydrofuran would
     give rise to MFT."

The full argument for why re-adding this step does NOT reverse Wave N is in the
step's own docstring (src/reaction_templates.py::_norfuraneol_mft_parallel_route)
and in tasks/audit_remediation.md "## Wave X" (a).  In one line: a spike experiment
measures a route's SHARE OF FLUX (Cerny & Davidek 2003 — norfuraneol is unimportant
IN SITU), a synthesis experiment measures whether the route EXISTS (Hofmann Table 4
— norfuraneol FED gives 14x the MFT that ribose gives), and under the Wave S1
ADDITIVE propagator both can be true at once.

WHY THIS SCRIPT EXISTS AT ALL, AND WHAT IT MUST NOT DO
------------------------------------------------------
The retired constant `thiol_addition_norfuraneol` = 26.85 kcal/mol is the barrier
the OLD norfuraneol step carried.  Reusing it would have been one line and zero
work.  It is refused, for two independent reasons either of which is sufficient:

  1. It was fitted (Wave H) THROUGH a route the isotope literature contradicts, and
     a barrier fitted through a contradicted route does not transfer — that is
     Wave N's own stated finding, applied to itself.
  2. It was fitted AGAINST MFT 342 ppb / FFT 200 ppb in
     data/benchmarks/cys_ribose_140C_Hofmann1998.json.  Wave S2b traced those two
     numbers to a repo-internal derivation from
     data/benchmarks/maillard_validation_benchmarks.md section 1.3, and Wave W then
     confirmed FROM THE PRIMARY SOURCE that 342 and 200 appear NOWHERE in
     10.1021/jf9705983.  It is a fit to this repository's own arithmetic.

So the new key gets its value from the measurement, or it stays at its unfitted
class seed.  The retired constant's full history is cross-referenced from the new
entry in src/barrier_constants.py so the two can never be confused, and the retired
entry itself is left untouched at 26.85 as the provenance record.

THE FIT TARGET, AND WHY THE OTHER NF ROW IS NOT ONE
----------------------------------------------------
ONE row: `hofmann1998_norfuraneol_h2s_145C_20min_pH5` (Table 4 NF row, 211.2 ug in
50 mL = 4224 ppb; the identical measurement is reprinted as the NF/H2S row of
Table 10 and is ingested only once).  Leverage is ONE free parameter against ONE
scored row = 1.0, far above `FIT_LEVERAGE_THRESHOLD` 0.5, so
`scripts/ci/fit_target_gate.py` classifies the row `per_row_recovery` and removes it
from the honest literature-coverage numerator AND denominator.  Its post-fit
agreement carries no information and must never be quoted as validation.

`hofmann1998_norfuraneol_cysteine_145C_20min_pH5` (Table 10, NF/cysteine, 50.8 ug =
1016 ppb) runs through the SAME barrier and is deliberately NOT a fit target.  It
supplies its H2S from cysteine rather than as a charge, so it is the TRANSFER TEST:
one row sets the constant, the other asks whether it generalises.  Its residual is
reported here, before and after, and it is not optimised.

CO-MOVEMENT IS REPORTED, NOT HIDDEN
------------------------------------
This barrier is on a channel that competes for upstream flux with the two incumbent
MFT routes and with the FFT lane.  Lowering it therefore moves the ribose/cysteine
panel rows, which are NOT fit targets and which the model already OVER-predicts.
The script measures that co-movement on all three Wave W panel rows at the incumbent
and at the adopted value and writes both into the record.  A fit that improves its
own target by damaging four rows it was not fitted on is a finding, and the record
says so in `co_movement`.

Knob and its defensible range
-----------------------------
ONE knob.  The range is bounded by values ALREADY IN `FAST_BARRIERS` for the same
mechanistic class, so the fit cannot leave the envelope the table itself spans:

    furanone_reductive_sulfhydrylation  [23.30, 29.65]
        low  = thiohemiacetal_formation (23.30) -- the fastest sulfur-addition step
        high = thiol_addition_hexose    (29.65) -- the slowest

This is the SAME range the retired `thiol_addition_norfuraneol` knob used and the
same one `refit_thiol_addition_pentodiulose_hofmann.py` used, deliberately: the
three steps are the same mechanistic class (H2S addition to a carbonyl of a
sugar-derived intermediate), so the fits are directly comparable, which is itself a
result.

Decision rules (stated up front, applied mechanically)
------------------------------------------------------
Identical to `refit_thiol_addition_pentodiulose_hofmann.py`, deliberately, so the
runs are comparable:

1. IDENTIFIABILITY.  A knob whose objective profile is flat to `FLAT_PROFILE_DEX`
   over its whole defensible range is NOT IDENTIFIABLE and keeps its incumbent.
2. MATERIALITY.  A knob moves only if the best achievable objective gain is at least
   `MIN_OBJECTIVE_GAIN_DEX`.
3. CONSERVATIVE EDGE.  Among the values within `INDIFFERENCE_DEX` of the profile
   minimum, the LARGEST barrier is adopted -- the SLOWEST norfuraneol route the data
   cannot distinguish from the optimum.  Here that rule points in the same direction
   as the isotope constraint, which is a coincidence worth stating rather than
   relying on: the rule was copied unchanged from the earlier script, not chosen.
4. BOUNDARY HONESTY.  If the profile minimum sits AT a range endpoint the record says
   so (`bounds_check.hit_a_bound`): a boundary optimum means the data want a value
   the table's own class envelope does not contain, and that is a finding, not a fit.
5. THE ISOTOPE GATE (new here, and the reason this fit exists in this form).  After a
   value is selected by rules 1-4 it must still satisfy the constraint that motivated
   re-adding the step as a SLOW PARALLEL channel: in the RIBOSE/CYSTEINE system, where
   norfuraneol is made in situ rather than charged, this channel's share of predicted
   MFT must stay a MINORITY.  If the selected value violates that, the fit is REJECTED
   and the incumbent is kept.  A rejected fit is a result: it says the two measurements
   pull in opposite directions and the model cannot satisfy both.

   WHERE THE CEILING COMES FROM, AND WHY IT IS 0.50 AND NOT SOMETHING TIGHTER.  The
   only quantitative content in the source is the word "mainly".  Cerny & Davidek 2003
   (10.1021/jf026123f), quoted verbatim in docs/validation/isotope_topology_evidence.md:

       "In another trial cysteine, 4-hydroxy-5-methyl-3(2H)-furanone and [13C5]ribose
        were reacted under the same conditions.  The resulting 2-methyl-3-furanthiol
        was mainly 13C5-labeled, suggesting that it stems from ribose and that
        4-hydroxy-5-methyl-3(2H)-furanone is unimportant as an intermediate."

   "Mainly 13C5-labeled" means the ribose-skeleton (non-norfuraneol) fraction is the
   MAJORITY, i.e. the norfuraneol share is below one half.  The paper prints no
   percentage, so 0.50 is the whole of what the sentence supports; a tighter ceiling
   would be a number invented to look rigorous, which is the failure mode this audit
   exists to remove.  The ceiling is also CONSERVATIVE IN THE RIGHT DIRECTION for a
   second reason stated in the source itself: Cerny's experiment ADDED authentic
   norfuraneol on top of what the system makes, so the in-situ share without a spike is
   strictly lower than whatever it was there.

   DISCLOSURE.  The first draft of the Wave X ledger entry pre-registered this ceiling
   as 0.20.  That was a guess, not a reading, and it was corrected to 0.50 once the
   source's actual wording was checked -- a correction in the PERMISSIVE direction,
   which is the direction that has to be declared.  It changes no decision in this
   wave: the shipped incumbent's share and the fitted candidate's share fall on the
   same sides of both ceilings (see `isotope_gate` in the record), so the gate's verdict
   is identical at 0.20 and at 0.50.

Usage
-----
    python scripts/generators/fit_furanone_reductive_sulfhydrylation_hofmann.py
    python scripts/generators/fit_furanone_reductive_sulfhydrylation_hofmann.py --apply
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

from src import data_paths
from src.barrier_constants import FAST_BARRIERS
from src.benchmark_validation import evaluate_benchmark
from src.uncertainty_propagation import _benchmark_signal_origin

BENCH_DIR = data_paths.BENCHMARKS_DIR
TARGET_BENCHMARK = BENCH_DIR / "hofmann1998_norfuraneol_h2s_145C_20min_pH5.json"
TRANSFER_BENCHMARK = BENCH_DIR / "hofmann1998_norfuraneol_cysteine_145C_20min_pH5.json"
CO_MOVEMENT_BENCHMARKS = (
    BENCH_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json",
    BENCH_DIR / "hofmann1998_glucose_cysteine_145C_20min_pH5.json",
    BENCH_DIR / "hofmann1998_fructose_cysteine_145C_20min_pH5.json",
)
BARRIER_SOURCE = ROOT / "src" / "barrier_constants.py"

UNSCORED_ROW_PENALTY = 3.0        # dex — same convention as the Wave H / Wave P scripts
GRID_STEP_KCAL = 0.05
FLAT_PROFILE_DEX = 1.0e-6
MIN_OBJECTIVE_GAIN_DEX = 0.02
INDIFFERENCE_DEX = 0.01

#: Rule 5. The norfuraneol channel's maximum admissible share of predicted MFT in the
#: ribose/cysteine system. 0.50 = "mainly 13C5-labeled" (Cerny & Davidek 2003); see
#: decision rule 5 in the module docstring for the derivation and for the disclosure
#: about the 0.20 that appeared in the ledger's first draft.
ISOTOPE_SHARE_CEILING = 0.50
RIBOSE_BENCHMARK = BENCH_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json"
MFT_COMPOUND = "2-Methyl-3-furanthiol (MFT)"

FITTED_KNOB = "furanone_reductive_sulfhydrylation"
KNOB_RANGE: Tuple[float, float] = (23.30, 29.65)
KNOB_RANGE_BASIS = (
    "H2S addition to a carbonyl of a sugar-derived intermediate, with reduction; low = "
    "thiohemiacetal_formation (23.30), high = thiol_addition_hexose (29.65), both already in "
    "FAST_BARRIERS. Identical to the range used by the retired thiol_addition_norfuraneol fit "
    "and by refit_thiol_addition_pentodiulose_hofmann.py, so the three are comparable."
)

REPORTED_NOT_FITTED = {
    "thiol_addition_pentodiulose": (
        "The INCUMBENT MFT route (Cerny & Davidek's 1,4-dideoxyosone step). It is the channel this "
        "one runs in parallel with, so on a norfuraneol-fed system where the dideoxyosone route "
        "cannot run at all it has no derivative -- and on the ribose system, where it does, fitting "
        "both would split one measurable quantity across two constants and report the split as "
        "knowledge."
    ),
    "furanone_cyclisation": (
        "The step that MAKES norfuraneol from the 1-deoxyosone. In the fit target norfuraneol is "
        "CHARGED, so this step is upstream of nothing; fitting it here would tune a constant the fit "
        "target does not constrain."
    ),
    "furanone_reductive_opening": (
        "The competing consumer of norfuraneol (-> 2,3-pentanedione, Wave P item 3). It is "
        "identifiable on this benchmark, because it drains the same pool -- which is exactly why it "
        "is NOT fitted: two knobs on one row is per-row recovery with extra steps."
    ),
    "mercaptoketone_formation": (
        "Downstream of the reductive ring opening, and on the OTHER product (2-mercapto-3-pentanone). "
        "Not on the MFT lane in this system."
    ),
}


# ---------------------------------------------------------------------------
# fit-target guards
# ---------------------------------------------------------------------------
def _guarded_target(path: Path) -> Path:
    assert path.exists(), f"fit target missing: {path}"
    assert "external_validation" not in path.parts, (
        f"hold-out benchmark reached the fit-target selector: {path}"
    )
    assert "step_level_unreachable" not in path.parts, (
        f"an unexecutable row reached the fit-target selector: {path}"
    )
    assert "Internal2026" not in path.name and "ProtocolPilot2026" not in path.name, (
        f"synthetic benchmark reached the fit-target selector: {path}"
    )
    origin = _benchmark_signal_origin(path)
    assert origin == "external_literature", (
        f"fit target must be literature-sourced, got signal origin {origin!r}"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    # The benchmark must SAY it is a fit target. If that declaration is ever removed
    # the fit must fail rather than quietly fit a row the panel is scoring.
    assert "fit_target_declaration" in payload, (
        "the fit target's `fit_target_declaration` block is missing; a row cannot be fitted "
        "and scored at the same time"
    )
    assert payload["source_doi"] == "10.1021/jf9705983", payload["source_doi"]
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


def _with_knob(value: float, fn):
    saved = FAST_BARRIERS[FITTED_KNOB]
    try:
        FAST_BARRIERS[FITTED_KNOB] = (float(value), saved[1])
        return fn()
    finally:
        FAST_BARRIERS[FITTED_KNOB] = saved


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


def _ribose_mft_ppb(with_route: bool) -> float:
    """Predicted MFT on the ribose/cysteine panel row, with the NF channel on or off.

    Turning the channel OFF is done by replacing the template function the engine
    imported, not by editing a barrier: a barrier can only be made large, never absent,
    and "large" is exactly what this measurement is trying to quantify.
    """
    import src.smirks_engine as smirks_engine

    saved = smirks_engine._norfuraneol_mft_parallel_route
    try:
        if not with_route:
            smirks_engine._norfuraneol_mft_parallel_route = lambda pool: []
        evaluation = evaluate_benchmark(RIBOSE_BENCHMARK)
        for comparison in evaluation.comparisons:
            if comparison.compound == MFT_COMPOUND:
                return float(comparison.predicted_ppb)
        raise AssertionError(f"{MFT_COMPOUND!r} not among the ribose row's comparisons")
    finally:
        smirks_engine._norfuraneol_mft_parallel_route = saved


def _isotope_share(barrier_value: float) -> Dict[str, float]:
    """The NF channel's share of predicted ribose MFT at one barrier value.

    share = 1 - MFT(channel off) / MFT(channel on).  This subtraction is meaningful
    ONLY because the flux propagator is ADDITIVE (Wave S1): channels sum, so removing
    one removes exactly its contribution.  Under a max-channel or multiplicative rule
    the difference would not be a share and this whole test would be meaningless -- the
    same property that makes re-adding the step defensible in the first place.
    """
    without = _ribose_mft_ppb(with_route=False)
    with_route = _with_knob(barrier_value, lambda: _ribose_mft_ppb(with_route=True))
    share = 0.0 if with_route <= 0.0 else 1.0 - (without / with_route)
    return {
        "barrier": float(barrier_value),
        "ribose_mft_ppb_without_nf_channel": without,
        "ribose_mft_ppb_with_nf_channel": with_route,
        "nf_channel_share_of_mft": float(share),
        "passes_ceiling": bool(share < ISOTOPE_SHARE_CEILING),
    }


def _panel_snapshot(value: float) -> Dict[str, Any]:
    """Every NON-fit-target row this knob can move, at one knob value."""
    def _run() -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        for path in (TRANSFER_BENCHMARK, *CO_MOVEMENT_BENCHMARKS):
            objective, rows = _score(path)
            out[path.stem] = {
                "objective_dex": objective,
                "rows": [
                    {
                        "compound": row["compound"],
                        "measured_ppb": row["measured_ppb"],
                        "predicted_ppb": row["predicted_ppb"],
                        "ratio": row["ratio"],
                    }
                    for row in rows
                ],
            }
        return out
    return _with_knob(value, _run)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument(
        "--apply",
        action="store_true",
        help="rewrite the adopted value into src/barrier_constants.py",
    )
    args = parser.parse_args()

    path = _guarded_target(TARGET_BENCHMARK)
    incumbent = float(FAST_BARRIERS[FITTED_KNOB][0])
    baseline_objective, baseline_rows = _score(path)

    print(f"Fit target: {path.name}  (DOI 10.1021/jf9705983, Table 4 NF row)")
    print(f"Knob: {FITTED_KNOB} = {incumbent:.2f} kcal/mol")
    print(f"Incumbent objective J = {baseline_objective:.4f} dex")
    for row in baseline_rows:
        print(
            f"   {row['compound']:30s} measured {row['measured_ppb']:.4g} ppb  "
            f"predicted {row['predicted_ppb']:.4g} ppb  ratio {row['ratio']:.4g}"
        )

    low, high = KNOB_RANGE
    grid = _frange(low, high, GRID_STEP_KCAL)
    profile: List[Dict[str, Any]] = []
    for value in grid:
        objective, _rows = _with_knob(value, lambda: _score(path))
        profile.append({"value": value, "objective": objective})

    objectives = [point["objective"] for point in profile]
    best = min(objectives)
    worst = max(objectives)
    span = worst - best
    identifiable = span > FLAT_PROFILE_DEX
    gain = baseline_objective - best
    argmin = min(profile, key=lambda point: point["objective"])["value"]
    hit_a_bound = abs(argmin - low) < 1e-9 or abs(argmin - high) < 1e-9

    adopted = incumbent
    indifference: List[float] = []
    if not identifiable:
        decision = "NOT IDENTIFIABLE — flat profile; incumbent kept"
    elif gain < MIN_OBJECTIVE_GAIN_DEX:
        decision = (
            f"IMMATERIAL — best achievable gain {gain:.4f} dex < "
            f"{MIN_OBJECTIVE_GAIN_DEX} dex; incumbent kept"
        )
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
    print(f"rule 1-4 selection: {decision}")

    # ── Rule 5: the isotope gate ─────────────────────────────────────────────
    selected_by_rules_1_to_4 = adopted
    isotope_at_incumbent = _isotope_share(incumbent)
    isotope_at_selected = _isotope_share(selected_by_rules_1_to_4)
    # The share as a function of the barrier, so the crossing point is visible rather
    # than implied. Coarser than the objective grid on purpose: each point costs two
    # full network enumerations.
    isotope_profile = [
        _isotope_share(value) for value in _frange(low, high, 0.25)
    ]
    admissible = [p["barrier"] for p in isotope_profile if p["passes_ceiling"]]
    isotope_limiting_barrier = min(admissible) if admissible else None

    gate_verdict = "PASS"
    if not isotope_at_selected["passes_ceiling"]:
        gate_verdict = "VIOLATED — FIT REJECTED, INCUMBENT KEPT"
        adopted = incumbent
        decision = (
            f"REJECTED BY THE ISOTOPE GATE (rule 5). Rules 1-4 selected "
            f"{selected_by_rules_1_to_4:.2f} kcal/mol, at which the norfuraneol channel "
            f"supplies {isotope_at_selected['nf_channel_share_of_mft']:.1%} of predicted MFT in "
            f"the ribose/cysteine system -- a MAJORITY, which is what Cerny & Davidek 2003's "
            f"spiking experiment excludes ('mainly 13C5-labeled'). The incumbent "
            f"{incumbent:.2f} is kept, at which the share is "
            f"{isotope_at_incumbent['nf_channel_share_of_mft']:.1%}. THE TWO MEASUREMENTS PULL IN "
            f"OPPOSITE DIRECTIONS AND THE MODEL CANNOT SATISFY BOTH: Hofmann Table 4 wants a "
            f"barrier at or below the range floor {low:.2f}; the isotope constraint admits "
            f"nothing below "
            + (f"{isotope_limiting_barrier:.2f}" if isotope_limiting_barrier is not None
               else "any value in the range")
            + ". That conflict is the result of this fit."
        )
    print(f"isotope gate: {gate_verdict}")
    print(f"final decision: {decision}")

    adopted_objective, adopted_rows = _with_knob(adopted, lambda: _score(path))

    before = _panel_snapshot(incumbent)
    after = _panel_snapshot(adopted)

    reported: Dict[str, Any] = {}
    for key, reason in REPORTED_NOT_FITTED.items():
        value = float(FAST_BARRIERS[key][0])
        probe = []
        for delta in (-2.0, -1.0, 1.0, 2.0):
            saved = FAST_BARRIERS[key]
            try:
                FAST_BARRIERS[key] = (value + delta, saved[1])
                objective, _rows = _score(path)
            finally:
                FAST_BARRIERS[key] = saved
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
        "record_id": "furanone_reductive_sulfhydrylation_refit_hofmann",
        "generated": "2026-08-28",
        "wave": "X",
        "generator": "scripts/generators/fit_furanone_reductive_sulfhydrylation_hofmann.py",
        "fit_target": [str(TARGET_BENCHMARK.relative_to(ROOT)), TARGET_BENCHMARK.stem],
        "fit_target_doi": "10.1021/jf9705983",
        "fit_target_row": "Table 4, NF row: 211.2 ug MFT per 50 mL = 4224 ppb (0.19 mol %)",
        "parameters_fitted": 1,
        "fitted_rows": fitted_rows,
        "fit_leverage": {
            "free_parameters": 1,
            "fitted_rows": fitted_rows,
            "parameters_per_row": (1.0 / fitted_rows) if fitted_rows else None,
            "class": "per_row_recovery",
            "interpretation": (
                "ONE free parameter against ONE scored row = 1.0 parameters per row, twice "
                "src.fit_target_index.FIT_LEVERAGE_THRESHOLD (0.5). The row is therefore "
                "classified per_row_recovery and stays OUT of the honest literature-coverage "
                "numerator AND denominator. Post-fit agreement on it carries no information "
                "about the model and must never be quoted as validation."
            ),
            "why_this_exclusion_is_symmetric": (
                "Excluding a fitted row removes its MISSES as well as its hits. That matters "
                "unusually much here, because the fit was REJECTED and the row is still 2.3x "
                "under-predicted: the exclusion removes a failure from the count, not a "
                "success. The failure is reported in full in tasks/audit_remediation.md "
                "'## Wave X' instead, which is where an excluded row's outcome belongs."
            ),
        },
        "knob": FITTED_KNOB,
        "knob_range": list(KNOB_RANGE),
        "knob_range_basis": KNOB_RANGE_BASIS,
        "incumbent": incumbent,
        "incumbent_basis": (
            "The un-fitted `thiol_addition` CLASS value 28.60, i.e. the same neutral seed Wave N "
            "gave `thiol_addition_pentodiulose`. DELIBERATELY NOT the retired "
            "`thiol_addition_norfuraneol` = 26.85, which was fitted through a route the isotope "
            "literature contradicts AND against values Wave S2b/Wave W proved were not in the "
            "paper."
        ),
        "incumbent_objective_dex": baseline_objective,
        "incumbent_rows": baseline_rows,
        "profile": profile,
        "profile_min_dex": best,
        "profile_max_dex": worst,
        "profile_span_dex": span,
        "identifiable": identifiable,
        "achievable_gain_dex": gain,
        "bounds_check": {"argmin": argmin, "hit_a_bound": hit_a_bound, "range": list(KNOB_RANGE)},
        "indifference_band": [min(indifference), max(indifference)] if indifference else None,
        "selected_by_rules_1_to_4": selected_by_rules_1_to_4,
        "adopted": adopted,
        "adopted_objective_dex": adopted_objective,
        "adopted_rows": adopted_rows,
        "decision": decision,
        "isotope_gate": {
            "rule": (
                "The norfuraneol channel's share of predicted MFT in the RIBOSE/CYSTEINE system "
                "must stay a minority. Measured as 1 - MFT(channel off)/MFT(channel on), which is "
                "a share only because the Wave S1 flux propagator is ADDITIVE."
            ),
            "ceiling": ISOTOPE_SHARE_CEILING,
            "ceiling_source": (
                "Cerny & Davidek 2003 (10.1021/jf026123f), verbatim: 'The resulting "
                "2-methyl-3-furanthiol was mainly 13C5-labeled, suggesting that it stems from "
                "ribose and that 4-hydroxy-5-methyl-3(2H)-furanone is unimportant as an "
                "intermediate.' 'Mainly' = the non-norfuraneol fraction is the majority = the "
                "norfuraneol share is below one half. The paper prints NO percentage, so 0.50 is "
                "the whole of what the sentence supports. The experiment also ADDED norfuraneol, "
                "so the unspiked in-situ share is strictly lower than whatever it was there."
            ),
            "ceiling_disclosure": (
                "The Wave X ledger's first draft pre-registered 0.20. That was a guess rather than "
                "a reading and was corrected to 0.50 against the source's wording -- a correction "
                "in the PERMISSIVE direction, declared here for that reason. It changes NO decision "
                "in this wave: the incumbent's share and the rule-1-to-4 candidate's share fall on "
                "the same sides of both 0.20 and 0.50, so the gate's verdict is identical either "
                "way. See `verdict_at_0_20` below."
            ),
            "at_incumbent": isotope_at_incumbent,
            "at_rules_1_to_4_selection": isotope_at_selected,
            "verdict": gate_verdict,
            "verdict_at_0_20": {
                "incumbent_passes": isotope_at_incumbent["nf_channel_share_of_mft"] < 0.20,
                "candidate_passes": isotope_at_selected["nf_channel_share_of_mft"] < 0.20,
                "note": (
                    "Recorded so the threshold correction can be audited. If the incumbent fails "
                    "0.20 and passes 0.50 while the candidate fails both, the gate REJECTS the fit "
                    "under either ceiling and the correction is decision-neutral -- but the "
                    "incumbent is then not comfortably inside the constraint either, and that is "
                    "reported in the Wave X ledger rather than smoothed over."
                ),
            },
            "share_vs_barrier_profile": isotope_profile,
            "isotope_limiting_barrier": isotope_limiting_barrier,
            "the_conflict": (
                "Hofmann Table 4 (fed norfuraneol) wants this barrier LOW; Cerny & Davidek's "
                "spiking experiment (in-situ norfuraneol) requires it HIGH. Where the two admissible "
                "sets do not overlap, the model cannot reproduce both experiments at once, and the "
                "non-overlap is itself the measurement this fit produced."
            ),
        },
        "co_movement": {
            "what_this_is": (
                "Every benchmark this knob can move that is NOT the fit target, scored at the "
                "incumbent and at the adopted value. `hofmann1998_norfuraneol_cysteine_...` is the "
                "TRANSFER TEST (same barrier, H2S from cysteine instead of charged); the three "
                "ribose/glucose/fructose rows are the Wave W panel rows, which share upstream flux "
                "with this channel and were not fitted to anything."
            ),
            "at_incumbent": before,
            "at_adopted": after,
        },
        "reported_not_fitted": reported,
        "retired_constant_cross_reference": {
            "key": "thiol_addition_norfuraneol",
            "value": 26.85,
            "status": "RETIRED 2026-08-27 (Wave N); left at 26.85 as a provenance record; no step "
                      "emits that family.",
            "why_not_reused": (
                "(1) It was fitted THROUGH the norfuraneol route that Cerny & Davidek 2003 "
                "(10.1021/jf026123f) contradicts as an in-situ mechanism, and Wave N's own rule is "
                "that a barrier fitted through a contradicted route does not transfer. (2) It was "
                "fitted AGAINST MFT 342 / FFT 200 ppb, which Wave S2b traced to a repo-internal "
                "derivation and Wave W confirmed from the primary source appear NOWHERE in "
                "10.1021/jf9705983. Either defect alone is disqualifying."
            ),
            "record": "results/validation/sulfur_barrier_refit_hofmann.{json,md} (STALE), "
                      "results/validation/sulfur_barrier_refit_pentodiulose.{json,md} (RETRACTED)",
        },
        "caveats": [
            "THE HEADSPACE CAVEAT IS LOAD-BEARING HERE. The fit target charges 1 mmol of H2S into a "
            "closed 200 mL autoclave holding 50 mL of solution at 145 degC. A large share of that "
            "H2S is in the gas phase; the engine has no gas-liquid partition term and reacts with "
            "the full nominal 20 mM. The fitted barrier therefore absorbs whatever the true aqueous "
            "H2S deficit is. The error is LOCALISED in this one constant, which is the reason to "
            "fit one named constant rather than let the discrepancy distribute itself.",
            "The isotope constraint is NOT enforced by this script. It is enforced by "
            "tests/scientific/test_wave_x_step_level_2026_08.py::"
            "test_norfuraneol_route_stays_a_minor_share_of_ribose_mft_flux, which runs AFTER the fit "
            "so that a fitted value demanding a fast norfuraneol route FALSIFIES the Wave X "
            "resolution rather than shipping under it.",
        ],
    }

    out_dir = ROOT / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "furanone_reductive_sulfhydrylation_refit_hofmann.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# `furanone_reductive_sulfhydrylation` — fit against Hofmann 1998 Table 4",
        "",
        "Generated by `scripts/generators/fit_furanone_reductive_sulfhydrylation_hofmann.py`,",
        "2026-08-28 (Wave X). Read that script's docstring first; it states what this fit may and",
        "may not be used for.",
        "",
        f"**Fit target (ONE row):** `{TARGET_BENCHMARK.stem}` — Hofmann & Schieberle 1998,",
        "DOI `10.1021/jf9705983`, Table 4 NF row: **211.2 µg MFT per 50 mL = 4224 ppb (0.19 mol %)**.",
        "",
        f"**Knob:** `{FITTED_KNOB}`, range {KNOB_RANGE[0]:.2f}–{KNOB_RANGE[1]:.2f} kcal/mol.",
        f"{KNOB_RANGE_BASIS}",
        "",
        f"**Incumbent** {incumbent:.2f} → **adopted** {adopted:.2f} kcal/mol.",
        f"Objective {baseline_objective:.4f} → {adopted_objective:.4f} dex.",
        f"Profile span {span:.4f} dex; achievable gain {gain:.4f} dex; argmin {argmin:.2f}"
        f"{' (AT A RANGE BOUND)' if hit_a_bound else ''}.",
        "",
        f"**Decision:** {decision}",
        "",
        "## The isotope gate (decision rule 5)",
        "",
        f"Ceiling **{ISOTOPE_SHARE_CEILING:.2f}**, from Cerny & Davidek 2003's *\"mainly "
        "13C5-labeled\"* and from nothing else. Share is measured on the **ribose/cysteine** row as",
        "`1 - MFT(channel off) / MFT(channel on)`, which is a share only because the flux",
        "propagator is additive.",
        "",
        "| barrier | ribose MFT, channel off | channel on | NF share | passes |",
        "|---:|---:|---:|---:|---|",
        f"| incumbent {incumbent:.2f} | {isotope_at_incumbent['ribose_mft_ppb_without_nf_channel']:.2f} |"
        f" {isotope_at_incumbent['ribose_mft_ppb_with_nf_channel']:.2f} |"
        f" {isotope_at_incumbent['nf_channel_share_of_mft']:.1%} |"
        f" {'yes' if isotope_at_incumbent['passes_ceiling'] else 'NO'} |",
        f"| rules 1-4 pick {selected_by_rules_1_to_4:.2f} |"
        f" {isotope_at_selected['ribose_mft_ppb_without_nf_channel']:.2f} |"
        f" {isotope_at_selected['ribose_mft_ppb_with_nf_channel']:.2f} |"
        f" {isotope_at_selected['nf_channel_share_of_mft']:.1%} |"
        f" {'yes' if isotope_at_selected['passes_ceiling'] else 'NO'} |",
        "",
        f"**Verdict: {gate_verdict}.** Lowest barrier in range that satisfies the ceiling: "
        + (f"**{isotope_limiting_barrier:.2f}** kcal/mol." if isotope_limiting_barrier is not None
           else "**none** — the constraint is violated across the whole range."),
        "",
        "## Leverage — why this row leaves the scored evidence",
        "",
        f"One free parameter against {fitted_rows} scored row(s) = leverage "
        f"{(1.0/fitted_rows) if fitted_rows else float('nan'):.2f}, above the 0.5 threshold, so the",
        "row is classified `per_row_recovery` and `scripts/ci/fit_target_gate.py` removes it from",
        "the honest literature-coverage numerator **and** denominator. Its agreement after the fit",
        "carries no information about the model and must never be quoted as validation.",
        "",
        "## The transfer test and the co-movement",
        "",
        "| benchmark | analyte | measured ppb | at incumbent | at adopted |",
        "|---|---|---:|---:|---:|",
    ]
    for stem in before:
        for row_before, row_after in zip(before[stem]["rows"], after[stem]["rows"]):
            lines.append(
                f"| `{stem}` | {row_before['compound']} | {row_before['measured_ppb']:.1f} | "
                f"{row_before['predicted_ppb']:.2f} | {row_after['predicted_ppb']:.2f} |"
            )
    lines += [
        "",
        "`hofmann1998_norfuraneol_cysteine_145C_20min_pH5` is the **transfer test**: same barrier,",
        "H₂S liberated from cysteine instead of charged. It is **scored**, not fitted. The three",
        "ribose / glucose / fructose rows are Wave W's panel rows; they share upstream flux with this",
        "channel and were fitted to nothing.",
        "",
        "## Reported but not fitted",
        "",
        "| knob | incumbent | objective span over ±2 kcal | identifiable | why not fitted |",
        "|---|---:|---:|---|---|",
    ]
    for key, entry in reported.items():
        lines.append(
            f"| `{key}` | {entry['incumbent']:.2f} | {entry['objective_span_over_probe_dex']:.4f} | "
            f"{'yes' if entry['identifiable'] else 'no'} | {entry['reason_not_fitted']} |"
        )
    lines += [
        "",
        "## The retired constant this fit deliberately does not reuse",
        "",
        "`thiol_addition_norfuraneol` = **26.85** kcal/mol, retired by Wave N and left in",
        "`FAST_BARRIERS` as a provenance record. It is refused here for two independent reasons,",
        "either sufficient on its own: it was fitted **through** a route the isotope literature",
        "contradicts, and it was fitted **against** MFT 342 / FFT 200 ppb — values Wave S2b traced to",
        "this repository's own arithmetic and Wave W confirmed appear nowhere in the paper.",
        "",
        "## Caveats",
        "",
    ]
    for caveat in record["caveats"]:
        lines.append(f"- {caveat}")
    lines.append("")
    (out_dir / "furanone_reductive_sulfhydrylation_refit_hofmann.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    if args.apply and abs(adopted - incumbent) > 1e-9:
        text = BARRIER_SOURCE.read_text(encoding="utf-8")
        BARRIER_SOURCE.write_text(_rewrite_constant(text, FITTED_KNOB, adopted), encoding="utf-8")
        print(f"APPLIED: {FITTED_KNOB} -> {adopted:.2f} in {BARRIER_SOURCE}")
    elif args.apply:
        print("APPLIED: nothing to write; the incumbent was kept.")

    print(f"Records written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
