import sys
from pathlib import Path
from dataclasses import replace


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    BenchmarkEvaluation,
    CompoundComparison,
    render_benchmark_summary_markdown,
    summarize_benchmarks,
    summarize_evaluation,
)


def test_benchmark_summary_separates_supported_and_unsupported_cases():
    summaries = summarize_benchmarks([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    ])

    by_id = {summary.benchmark_id: summary for summary in summaries}

    assert by_id["cys_glucose_150C_Farmer1999"].supported is True
    assert by_id["cys_glucose_150C_Farmer1999"].overall_status in {"pass", "scale-gap", "ranking-gap", "partial-pass"}

    assert by_id["pea_isolate_40C_PratapSingh2021"].supported is True
    assert by_id["pea_isolate_40C_PratapSingh2021"].strict_ready is False
    assert by_id["pea_isolate_40C_PratapSingh2021"].overall_status in {"pass", "partial-pass", "scale-gap", "ranking-gap"}


def test_benchmark_summary_markdown_includes_gap_labels():
    summaries = summarize_benchmarks([
        ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    ])

    markdown = render_benchmark_summary_markdown(summaries)

    assert "Benchmark Summary" in markdown
    assert "cys_glucose_150C_Farmer1999" in markdown
    assert "pea_isolate_40C_PratapSingh2021" in markdown
    assert "matrix-only intake path is executable" in markdown
    assert "Strict Ready" in markdown


def test_strict_gate_summary_reflects_threshold_failures():
    evaluation = BenchmarkEvaluation(
        benchmark_id="synthetic_benchmark",
        bench_file=ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=[
            CompoundComparison("cmp1", 10.0, 10.0, "cmp1", None, 1.0),
            CompoundComparison("cmp2", 20.0, 10.0, "cmp2", None, 1.0),
            CompoundComparison("cmp3", 30.0, 10.0, "cmp3", None, 1.0),
        ],
        pearson_r=0.80,
        mae_ppb=10.0,
    )

    summary = summarize_evaluation(evaluation, protein_type="free")

    assert summary.strict_ready is False
    assert summary.overall_status == "ranking-gap"
    assert any("ranking" in issue for issue in summary.blocking_issues)

    strict_pass = summarize_evaluation(
        replace(evaluation, pearson_r=0.95, comparisons=[
            CompoundComparison("cmp1", 10.0, 10.0, "cmp1", None, 1.0),
            CompoundComparison("cmp2", 20.0, 18.0, "cmp2", None, 1.0),
            CompoundComparison("cmp3", 30.0, 28.0, "cmp3", None, 1.0),
        ]),
        protein_type="free",
    )

    assert strict_pass.strict_ready is True
    assert strict_pass.blocking_issues == []