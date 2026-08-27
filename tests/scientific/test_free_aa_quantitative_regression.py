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
BENCHMARK_EXPECTED_FOLD_ERRORS = {
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json": {
        # compound -> (low, pinned, high); all UNDER-predictions
        "2-methyl-3-furanthiol": (4.0, 5.584, 7.5),
        "2-furfurylthiol": (2.4, 3.255, 4.4),
    },
}


def test_primary_free_amino_acid_benchmarks_under_predict_by_the_pinned_amount():
    for bench_path, limits in BENCHMARK_EXPECTED_FOLD_ERRORS.items():
        evaluation = evaluate_benchmark(bench_path)

        assert evaluation.supported, evaluation.reason
        assert evaluation.coverage == 1.0

        comparisons = {comparison.compound: comparison for comparison in evaluation.comparisons}
        for comparison in evaluation.comparisons:
            print(f"DEBUG: {evaluation.benchmark_id} | {comparison.compound} | Measured: {comparison.measured_ppb:.2f} | Predicted: {comparison.predicted_ppb:.2f} | Ratio: {comparison.ratio:.3f}")
        for compound, (low, pinned, high) in limits.items():
            comparison = comparisons[compound]
            assert comparison.predicted_ppb < comparison.measured_ppb, (
                f"{evaluation.benchmark_id} {compound} is a known UNDER-prediction "
                f"({pinned:.3f}x). It now over-predicts "
                f"({comparison.predicted_ppb:.2f} vs {comparison.measured_ppb:.2f} ppb) -- "
                "re-derive this expectation, do not relax it."
            )
            assert low <= comparison.ratio <= high, (
                f"{evaluation.benchmark_id} {compound}: fold error {comparison.ratio:.3f}x "
                f"is outside the pinned band [{low}, {high}] around {pinned:.3f}x. "
                "Something beyond the 2026-08-27 approved chemistry changes has moved the "
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