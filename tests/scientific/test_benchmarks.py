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
def test_benchmark_is_executable_and_every_measured_compound_is_detected(bench_file):
    """Executability + species-detection smoke check. NOT an accuracy or correlation test.

    RENAMED 2026-08-27 (Wave J2, red-team finding: deceptive test names). The previous name
    was ``test_benchmark_correlation``, which is what a reader of a green CI badge would
    take as evidence that predicted ppb tracks measured ppb across the panel. It never
    asserted anything of the kind -- there is no correlation coefficient, no ratio bound and
    no error bound anywhere below. What it actually establishes, per benchmark file, is:

      * the file loads and the benchmark runs end-to-end through the evaluator;
      * it produces at least one comparison row and a finite MAE;
      * every compound the source measured is DETECTED by the network (``matched_name`` is
        not None) -- i.e. the model can reach the species, at any magnitude whatsoever.

    Detection is a coverage property, not an agreement property. The panel's actual
    agreement is far worse than the green tick here suggests and is pinned honestly in
    ``tests/scientific/test_honest_headline_guards.py``: 0/6 predictive benchmarks are free
    of blocking gaps, 0/14 are strict-ready, and the worst quantitative ratio in this same
    panel is ~1204x. Read that file, not this one, for how well the model predicts.
    """
    evaluation = evaluate_benchmark(bench_file)
    bench = load_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type=bench.get("protein_type", "free"))

    # TIGHTENED 2026-08-27 (Wave J2): these two branches were `pytest.xfail(...)` calls.
    # `pytest.xfail()` is imperative -- it aborts the test at that line -- so a benchmark
    # that became unexecutable, or that lost the ability to detect some of its own measured
    # compounds, silently converted this test from a pass into an xfail rather than into a
    # failure. That is a self-excusing escape hatch on precisely the property under test.
    # Measured on the current tree (2026-08-27, results/validation/benchmark_summary.json):
    # all 14 panel benchmarks are `supported` and all sit at coverage == 1.0, so both
    # branches were dead code and hardening them costs nothing and lands green.
    assert evaluation.supported, (
        f"{evaluation.benchmark_id} is no longer executable through MaillardPipeline: "
        f"{evaluation.reason or 'no reason given'}"
    )

    assert evaluation.comparisons, f"No comparisons produced for {evaluation.benchmark_id}"
    assert evaluation.coverage == 1.0, (
        f"Incomplete species coverage for {evaluation.benchmark_id}: {evaluation.coverage:.2%}. "
        f"Every measured compound in a panel benchmark must be reachable by the network; "
        f"a coverage regression must go red, not xfail."
    )
    assert evaluation.mae_ppb is not None

    if STRICT_BENCHMARKS and summary.tier == "PRIMARY":
        assert summary.strict_ready, (
            f"{evaluation.benchmark_id} (Tier: {summary.tier}, Path: {summary.execution_path}) failed "
            f"strict validation gate: {'; '.join(summary.blocking_issues) or summary.overall_status}"
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

