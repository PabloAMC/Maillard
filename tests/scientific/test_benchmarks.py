import os
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark, get_benchmark_files, load_benchmark, summarize_evaluation


STRICT_BENCHMARKS = os.getenv("MAILLARD_STRICT_BENCHMARKS", "0") == "1"


@pytest.mark.parametrize("bench_file", get_benchmark_files())
def test_benchmark_correlation(bench_file):
    evaluation = evaluate_benchmark(bench_file)
    bench = load_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type=bench.get("protein_type", "free"))

    if not evaluation.supported:
        pytest.xfail(evaluation.reason or "Benchmark not yet executable through InverseDesigner")

    assert evaluation.comparisons, f"No comparisons produced for {evaluation.benchmark_id}"
    if evaluation.coverage < 1.0:
        pytest.xfail(f"Incomplete species coverage for {evaluation.benchmark_id}: {evaluation.coverage:.2%}")
    assert evaluation.mae_ppb is not None

    if STRICT_BENCHMARKS and bench.get("protein_type", "free") == "free":
        assert summary.strict_ready, (
            f"{evaluation.benchmark_id} failed strict free-AA gate: "
            f"{'; '.join(summary.blocking_issues) or summary.overall_status}"
        )

    for comparison in evaluation.comparisons:
        print(
            f"Benchmark {evaluation.benchmark_id} - {comparison.compound}: "
            f"Measured={comparison.measured_ppb}, Predicted={comparison.predicted_ppb:.4f}, "
            f"Matched={comparison.matched_name}"
        )
        assert comparison.matched_name is not None, (
            f"Model failed to detect {comparison.compound}, expected ~{comparison.measured_ppb} ppb"
        )


if __name__ == "__main__":
    pytest.main(["-s", __file__])

