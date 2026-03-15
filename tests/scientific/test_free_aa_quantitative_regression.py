import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark


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