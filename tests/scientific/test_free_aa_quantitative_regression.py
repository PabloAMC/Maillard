import sys
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark, evaluate_benchmark_payload


# cys_ribose_140C_Hofmann1998 is the only free-amino-acid panel benchmark with a verified
# source, so it is the only entry here.
#
# Two further entries used to sit below it: cys_ribose_150C_Mottram1994 and
# cys_glucose_150C_Farmer1999. They were quarantined on 2026-08-26 for provenance and then
# DELETED on the same date, once source recovery established that no literature source exists
# for either and that their values were invented (findings kept in
# data/benchmarks/quarantined/README.md). Ratio bounds against invented reference values are
# not drift guards -- they are guards against drifting away from a fabrication -- so they were
# removed rather than retargeted. The coverage regression they also carried now lives in
# tests/scientific/test_pentose_hexose_sulfur_ordering.py on synthetic formulations.
# RE-DERIVED 2026-08-27 (Wave H). This table used to hold one-sided upper bounds (MFT
# <= 2.0x, FFT <= 1.55x) and the test was named "...stay_locally_calibrated". Neither is
# true any more, and the honest statement is a two-sided pin on a known, quantified
# under-prediction rather than a widened ceiling:
#
#   * Wave G1 replaced the fabricated one-step 3-deoxyosone + 2 H2S MFT shortcut with the
#     accepted 1-deoxyosone -> norfuraneol -> MFT route (van den Ouweland & Peer 1975)
#     and deleted the fabricated lipid radical chain. Sulfur absolute yields fell 5-40x;
#     this benchmark went to 7.83x under on MFT and 3.19x on FFT.
#   * Wave H refit the sulfur branch against THIS benchmark and nothing else (it is the
#     only surviving literature constraint on that branch). The refit moved
#     `thiol_addition_norfuraneol` 28.60 -> 26.85 and recovered MFT to 5.58x under.
#     It also established, by scanning every candidate knob over its full defensible
#     range, that NO barrier value can do better: `furanone_cyclisation` and
#     `thiohemiacetal_formation` have exactly zero derivative here, and the objective
#     saturates at 0.61 dex with every knob pinned at the bottom of its range. The
#     remaining ~5x is a volatile-budget ALLOCATION deficit -- furfural, which this
#     experiment did not measure, takes ~78% of a total budget that is itself the right
#     order of magnitude (~1050 ppb against a measured MFT+FFT of 542 ppb).
#     Full profile: results/validation/sulfur_barrier_refit_hofmann.md.
#
# The benchmark's own contract (1.45x / 0.09 dex) is UNTOUCHED, so it fails visibly in
# every regenerated artifact. Closing this gap means fixing the allocation layer, not
# widening these numbers.
#
# CORRECTED 2026-08-27 (Wave S2c): the sentence above calling this benchmark "the only
# surviving literature constraint on that branch" is FALSE and is kept because it records
# what the repo believed. It was never a literature constraint -- see the Wave S2c block on
# the MFT band below. The contract described as UNTOUCHED is now RETIRED (not widened), and
# the tier is demoted PRIMARY -> REFERENCE.
# THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.
#
# RE-DERIVED AGAIN 2026-08-27 (Wave I fix 8 + fix 12). The Wave H numbers above are now
# WRONG, and the reason matters more than the new numbers do, so it is recorded here
# rather than folded away:
#
#   * Wave I fix 8 (red-team finding H4) found that the norfuraneol + H2S -> MFT step
#     was pool-gated on the `[HH]` 2[H] reducing-equivalent token whose ONLY producer
#     reachable from a ribose/cysteine system was `Aminoketone_Condensation` -- the
#     pyrazine aromatisation. The entire sulfur branch was therefore running behind a
#     pyrazine bottleneck, and disabling the pyrazine lane drove MFT to exactly 0.0 ppb.
#     Giving the token its own source (`2 cysteine -> cystine + 2[H]`,
#     src/reaction_templates.py::_thiol_reductant_pool) moved this benchmark to
#     MFT 345.04 ppb / FFT 187.49 ppb.
#   * Wave I fix 12 then restricted the demoted one-step MFT lump
#     (`Thiol_Addition_Legacy_Shortcut`) to hexoses, as its own docstring had claimed
#     all along. It had been double-producing pentose MFT alongside the accepted
#     norfuraneol route. That took MFT to 235.32 ppb and, with less competition for the
#     sulfur budget, took FFT UP to 219.96 ppb.
#
# NO BARRIER WAS TOUCHED in Wave I. `thiol_addition_norfuraneol` is still the 26.85 that
# Wave H fitted. THE WAVE H DIAGNOSIS QUOTED ABOVE IS THEREFORE SUPERSEDED: "no barrier
# value can do better", "the objective saturates at 0.61 dex" and "the remaining ~5x is a
# volatile-budget ALLOCATION deficit" were all measured on a network whose sulfur branch
# was throttled by the H4 coupling. results/validation/sulfur_barrier_refit_hofmann.md
# must be re-run before it is cited again, and 26.85 is now an UNREFIT constant that
# happens to sit near a benchmark it was fitted to under different mechanics. The
# agreement below is a coincidence of two independent choices and MUST NOT be reported as
# a calibration.
#
# The direction also changed, which is why this table now carries one: MFT still
# under-predicts, FFT now OVER-predicts slightly. The bands stay deliberately wide (they
# are drift guards, not accuracy claims) and the benchmark's own contract still FAILS --
# as of the 2026-08-27 Wave S1b pH/aw routing repair, max_ratio 2.2086 against a 1.45
# threshold (OUTSIDE) and MALE 0.2352 against 0.09 (OUTSIDE), i.e. it fails BOTH criteria
# by more than before. (History: 1.4864 / 0.1267 after the Wave S1 additive propagator;
# 1.4110 / 0.0935 after the Wave P refit -- one criterion;
# 2.2519 / 0.2192 after the Wave N route correction; 1.4533 / 0.1019 under Wave I.)
# Nothing here turns a failing benchmark into a passing one, and where an earlier wave's
# improvement was FIT RECOVERY on a declared fit target it is still labelled as such --
# see the per-compound note below.
BENCHMARK_EXPECTED_FOLD_ERRORS = {
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json": {
        # compound -> (low, pinned, high, direction). `ratio` is the symmetric fold
        # error, so the direction has to be pinned separately.
        #
        # RE-DERIVED 2026-08-27 (Wave N -- MFT ROUTE CORRECTION). Both bands moved because
        # the norfuraneol -> MFT step was RETIRED on isotope evidence: Cerny & Davidek 2003
        # (10.1021/jf026123f) spiked authentic norfuraneol into a [13C5]ribose/cysteine
        # system and the MFT came out mainly 13C5-labelled, i.e. "4-hydroxy-5-methyl-3(2H)-
        # furanone is unimportant as an intermediate"; Cerny & Davidek 2004
        # (10.1021/jf035265m) positionally confirmed 1,4-dideoxypento-2,3-diulose with
        # [1-13C]ribose. The route is now Deoxyosone_Reduction + Thiol_Addition_Pentodiulose.
        #
        # THE NUMBERS GOT WORSE AND WERE NOT CLAWED BACK. The 1.453x had been bought with
        # `thiol_addition_norfuraneol` = 26.85 kcal/mol, a barrier Wave H fitted THROUGH the
        # contradicted route; its replacement ships at the un-fitted sulfur-addition class
        # value 28.60 (ESTIMATED, UNCONSTRAINED). Nothing here was tuned to recover the lost
        # agreement, and the benchmark's own contract (1.45x / 0.09 dex) is still UNTOUCHED,
        # so it fails visibly in every regenerated artifact -- now at max_ratio 2.2519 and
        # MALE 0.2192 rather than 1.4533 / 0.1019.
        #
        # RE-DERIVED 2026-08-27 (Wave P item 1 -- THE OWNER-APPROVED REFIT OF THE
        # CORRECTED ROUTE). Wave N deliberately shipped `thiol_addition_pentodiulose` at
        # the un-fitted class value 28.60 and left the refit as an owner decision; Wave P
        # took it, and ran it LAST in the wave so the fit sees the network that ships.
        # ONE free parameter over [23.30, 29.65], decision rules identical to Wave H's
        # script. Record: results/validation/sulfur_barrier_refit_pentodiulose.{json,md};
        # generator: scripts/generators/refit_thiol_addition_pentodiulose_hofmann.py.
        #
        # BOTH ROWS IMPROVED, AND THAT IS NOT A VALIDATION RESULT. This benchmark is a
        # declared FIT TARGET at `per_row_recovery` leverage; a fitted row getting closer
        # to its own target is arithmetic, not evidence. Two things make the improvement
        # readable anyway, and both are pinned here:
        #   * the profile MINIMUM sits at the range FLOOR (argmin 23.30, hit_a_bound =
        #     true), so 26.35 is the conservative edge of the indifference band and the
        #     residual is NOT removable by this barrier -- with the knob pinned at the
        #     floor the objective is still 0.0836 dex;
        #   * FFT was NOT fitted and moved anyway, in the OPPOSITE direction, because the
        #     two lanes draw on the same upstream sugar flux. That coupling is exactly why
        #     ONE knob was fitted against these two rows rather than two.
        #
        # AND THE CONTRACT STILL FAILS, by a hair and on the second criterion only:
        # max_ratio 1.4110 is now INSIDE the 1.45 threshold, but MALE 0.0935 is outside
        # 0.09. The contract is UNTOUCHED and the benchmark still reports `scale-gap`.
        #
        # RE-DERIVED 2026-08-27 (Wave S1 -- THE ADDITIVE FLUX PROPAGATOR). NO CONSTANT
        # MOVED. `src/recommend.py::predict_from_steps` used to keep the lowest-span route
        # per product and discard every parallel route's flux; it now SUMS kinetically
        # distinct routes (deduplicated on the route's full step-set) before the fixed
        # volatile budget is allocated. Both compounds here are reached by two enumerated
        # routes, so BOTH rose:
        #     MFT 242.38 -> 283.59 ppb vs 342   1.4110x under -> 1.2060x under  (better)
        #     FFT 217.99 -> 297.28 ppb vs 200   1.0900x over  -> 1.4864x over   (WORSE)
        # THE BENCHMARK GOT WORSE OVERALL AND IS PINNED WORSE: MALE 0.0935 -> 0.1267 dex
        # and max_ratio 1.4110 -> 1.4864, so the untouched contract (1.45x / 0.09 dex) now
        # fails on BOTH criteria again instead of one. Nothing was clawed back -- the two
        # lanes share their upstream trunk, so any barrier that pushed FFT down would take
        # MFT down with it, and re-fitting a barrier to absorb a propagator change is
        # exactly the move this campaign exists to remove.
        #
        # RE-DERIVED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR).
        # NO CONSTANT MOVED. `get_ph_multiplier`, the enolisation route-selection term, had
        # never been called on the prediction path; it now is. At this benchmark's pH 5.0
        # it boosts `Enolisation_1_2` (the 3-deoxyosone -> furfural/FFT arm) 4.5x and leaves
        # `Enolisation_2_3_Amadori` (the 1-deoxyosone -> MFT arm) at ~1.0, so the fixed
        # volatile budget moves from MFT to FFT:
        #     MFT 283.59 -> 154.85 ppb vs 342   1.2060x under -> 2.2086x under  (WORSE)
        #     FFT 297.28 -> 267.50 ppb vs 200   1.4864x over  -> 1.3375x over   (better)
        # DO NOT READ THE MFT ROW AS A CONFLICT WITH A MEASUREMENT. Wave S1b's first
        # draft of this comment did, and Wave S2b (same day) showed that was wrong: the
        # 342 / 200 ppb targets were derived INSIDE THIS REPOSITORY from
        # data/benchmarks/maillard_validation_benchmarks.md section 1.3, an
        # abstract-reconstructed range table (MFT ~0.02-0.05 mol %, FFT ~0.01-0.03 mol %)
        # committed in the SAME commit as the benchmark file. 342 and 200 are interior
        # points of two INVENTED and OVERLAPPING bands (MFT 228-571 ppb, FFT 114-342 ppb),
        # so the MFT > FFT ordering is an artifact of midpoint selection, and the
        # 1.45x / 0.09 dex contract is ~1.7x tighter than its own source band. The
        # mechanism above is real and the degradation is real; what is NOT available is a
        # measurement to say the model is wrong about it. See
        # `content_verification_note.wave_s2_followup` in the benchmark file and the
        # '## Wave S2b' section of tasks/audit_remediation.md.
        #
        # RE-DERIVED 2026-08-27 (Wave S2c -- THE ANCHOR RETIREMENT AND THE BARRIER REVERT).
        # THIS IS NOT A CHEMISTRY CHANGE. Wave S2b settled where 342 / 200 came from and the
        # answer is: from this repository. The two values are interior points of two invented,
        # OVERLAPPING mol % bands in data/benchmarks/maillard_validation_benchmarks.md
        # section 1.3, an abstract-reconstructed table committed in c7efbbc -- the SAME commit
        # that created the benchmark JSON. On the file's declared (unattested) 10 mM basis with
        # MW 114.17: 0.0300 mol % -> 342.5 -> 342 ppb, and the FFT band's geometric mean
        # 0.017321 mol % -> 197.8 -> 200 ppb. ~90% confidence, arithmetic exact.
        # CONSEQUENCE EXECUTED HERE: `thiol_addition_pentodiulose` REVERTED 26.35 -> 28.60 (the
        # un-fitted Wave N class value) and results/validation/sulfur_barrier_refit_pentodiulose
        # RETRACTED, because that refit's sole target was this benchmark. The benchmark's own
        # 1.45x / 0.09 dex contract is RETIRED -- not widened -- and its tier demoted
        # PRIMARY -> REFERENCE; it now inherits the global free-precursor default 1.5x / 0.10 dex
        # and fails it by more: max_ratio 4.3797, MALE 0.4041 dex.
        # WHAT THESE BANDS MEAN NOW, stated so nobody over-reads them: they are drift guards on
        # the model's own output, measured against a number this repository invented. They can
        # detect a silent change in the sulfur lane. They cannot say anything about accuracy,
        # in either direction, and no wave may cite this benchmark's fold error as agreement or
        # disagreement with literature. THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.
        # 2026-08-27 (Wave S2c): measured 78.09 vs 342 ppb (was 154.85 under Wave S1b, 283.59
        # under Wave S1, 242.38 under Wave P, 151.87 under Wave N, 235.32 under Wave I). The
        # band carries the SAME RELATIVE SPAN as the Wave S1b pin (x1.244 up, clipped at 1.00
        # down), so it is a re-centring on a much WORSE number, not a loosening.
        "2-methyl-3-furanthiol": (1.00, 4.380, 5.45, "under"),
        # 2026-08-27 (Wave S1b): measured 267.50 vs 200 ppb (was 297.28 under Wave S1,
        # 217.99 under Wave P, 243.72 under Wave N, 219.96 under Wave I). Still an
        # OVER-prediction; the direction did not change, and this is the half of the pair
        # that IMPROVED. The band is TIGHTENED to the Wave S1 pin's relative span applied
        # to the smaller number (x1.312 -> upper 1.76), so improving here does not buy
        # slack for the next wave.
        # 2026-08-27 (Wave S2c): measured 293.67 vs 200 ppb (was 267.50 under Wave S1b, 297.28
        # under Wave S1, 217.99 under Wave P, 243.72 under Wave N, 219.96 under Wave I). This
        # half got WORSE again under the barrier revert -- the two lanes share their upstream
        # trunk, so lifting the MFT-lane barrier hands budget share to the FFT lane. Same
        # relative span as the Wave S1b pin (x1.3154 up), applied to the larger number: a
        # re-centring on a worse value, not a loosening.
        "2-furfurylthiol": (1.00, 1.468, 1.93, "over"),
    },
}


def test_primary_free_amino_acid_benchmarks_land_in_the_pinned_fold_error_bands():
    for bench_path, limits in BENCHMARK_EXPECTED_FOLD_ERRORS.items():
        evaluation = evaluate_benchmark(bench_path)

        assert evaluation.supported, evaluation.reason
        assert evaluation.coverage == 1.0

        comparisons = {comparison.compound: comparison for comparison in evaluation.comparisons}
        for comparison in evaluation.comparisons:
            print(f"DEBUG: {evaluation.benchmark_id} | {comparison.compound} | Measured: {comparison.measured_ppb:.2f} | Predicted: {comparison.predicted_ppb:.2f} | Ratio: {comparison.ratio:.3f}")
        for compound, (low, pinned, high, direction) in limits.items():
            comparison = comparisons[compound]
            observed = "under" if comparison.predicted_ppb < comparison.measured_ppb else "over"
            assert observed == direction, (
                f"{evaluation.benchmark_id} {compound} was pinned as an {direction.upper()}-"
                f"prediction ({pinned:.3f}x) on 2026-08-27 (Wave N) and now {observed}-predicts "
                f"({comparison.predicted_ppb:.2f} vs {comparison.measured_ppb:.2f} ppb) -- "
                "re-derive this expectation with a dated cause, do not relax it."
            )
            assert low <= comparison.ratio <= high, (
                f"{evaluation.benchmark_id} {compound}: fold error {comparison.ratio:.3f}x "
                f"is outside the pinned band [{low}, {high}] around {pinned:.3f}x. "
                "Something beyond the 2026-08-27 approved chemistry changes (Wave N MFT "
                "route correction included) has moved the "
                "sulfur branch -- investigate before re-pinning. If the gap has genuinely "
                "closed, this test should become a calibration assertion again."
            )


def test_ribose_beats_glucose_for_equal_condition_mft_prediction():
    base_payload = {
        "benchmark_id": "donor_identity_equal_condition_mft",
        "source_doi": "synthetic_equal_condition_regression",
        "precursors": {
            "L-Cysteine": {"concentration_mM": 10.0},
            "D-Ribose": {"concentration_mM": 10.0},
        },
        "conditions": {
            "temp_C": 150.0,
            "ph": 5.5,
            "water_activity": 0.95,
            "time_min": 60.0,
        },
        "metadata": {
            "tier": "PRIMARY",
            "family": "free_aa_sugar_discrimination",
            "execution_path": "free_precursor",
        },
        "validation_contract": {
            "scale_thresholds": {
                "max_ratio": 8.0,
                "mean_abs_log10_error": 1.25,
            }
        },
        "measured_volatiles": {
            "2-methyl-3-furanthiol": {"conc_ppb": 60.0, "uncertainty_pct": 30},
        },
    }

    ribose_eval = evaluate_benchmark_payload(base_payload)
    glucose_payload = deepcopy(base_payload)
    glucose_payload["benchmark_id"] = "donor_identity_equal_condition_mft_glucose"
    glucose_payload["precursors"] = {
        "L-Cysteine": {"concentration_mM": 10.0},
        "D-Glucose": {"concentration_mM": 10.0},
    }
    glucose_eval = evaluate_benchmark_payload(glucose_payload)

    assert ribose_eval.supported, ribose_eval.reason
    assert glucose_eval.supported, glucose_eval.reason
    assert ribose_eval.comparisons[0].compound == "2-methyl-3-furanthiol"
    assert ribose_eval.comparisons[0].predicted_ppb > glucose_eval.comparisons[0].predicted_ppb