import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    evaluate_benchmark,
    snapshot_all_benchmark_targets,
    snapshot_benchmark_targets,
    summarize_evaluation,
)


def test_pea_isolate_matrix_only_benchmark_is_executable_with_full_coverage():
    bench_file = ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json"

    evaluation = evaluate_benchmark(bench_file)

    assert evaluation.supported is True
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 3


def test_pea_isolate_matrix_only_benchmark_preserves_measured_ordering_without_entering_strict_gate():
    bench_file = ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json"

    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type="pea_iso")
    ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}
    predicted = {comparison.compound: comparison.predicted_ppb for comparison in evaluation.comparisons}

    assert predicted["2-pentylfuran"] > predicted["hexanal"] > predicted["hexanol"]
    assert ratios["hexanal"] <= 1.25
    assert ratios["2-pentylfuran"] <= 1.2
    assert ratios["hexanol"] <= 1.35
    assert summary.strict_ready is False


def test_matrix_only_benchmark_is_deliberately_excluded_from_target_snapshots():
    bench_file = ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json"

    assert snapshot_benchmark_targets(bench_file) == []

    rows = snapshot_all_benchmark_targets([
        bench_file,
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
    ])

    assert rows
    assert {row.benchmark_id for row in rows} == {"cys_glucose_150C_Farmer1999"}