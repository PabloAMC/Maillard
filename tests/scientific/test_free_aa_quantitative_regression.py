import sys
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark, evaluate_benchmark_payload


BENCHMARK_LIMITS = {
    ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json": {
        "2-methyl-3-furanthiol": 2.0,
        "2-furfurylthiol": 1.55,
    },
    ROOT / "data" / "benchmarks" / "cys_ribose_150C_Mottram1994.json": {
        "2-methyl-3-furanthiol": 1.35,
        "bis(2-methyl-3-furyl) disulfide": 1.05,
        "furfural": 1.85,
    },
    ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json": {
        "2-methyl-3-furanthiol": 1.5,
        "furfural": 1.25,
        "pyrazine": 1.25,
    },
}


def test_primary_free_amino_acid_benchmarks_stay_locally_calibrated():
    for bench_path, limits in BENCHMARK_LIMITS.items():
        evaluation = evaluate_benchmark(bench_path)

        assert evaluation.supported, evaluation.reason
        assert evaluation.coverage == 1.0

        ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}
        for comparison in evaluation.comparisons:
            print(f"DEBUG: {evaluation.benchmark_id} | {comparison.compound} | Measured: {comparison.measured_ppb:.2f} | Predicted: {comparison.predicted_ppb:.2f} | Ratio: {comparison.ratio:.3f}")
        for compound, max_ratio in limits.items():
            assert ratios[compound] <= max_ratio, (
                f"{evaluation.benchmark_id} {compound}: ratio {ratios[compound]:.3f} exceeds {max_ratio:.3f}"
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