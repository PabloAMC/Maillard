import sys
from pathlib import Path
from dataclasses import replace


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    build_family_lane_validation_artifact,
    render_family_lane_validation_markdown,
    summarize_benchmarks,
    summarize_evaluation,
)
from src.benchmark_types import BenchmarkEvaluation, CompoundComparison
from src.presentation import render_benchmark_summary_markdown


def test_benchmark_summary_separates_supported_and_unsupported_cases(benchmark_summary_matrix_scope):
    summaries = benchmark_summary_matrix_scope

    by_id = {summary.benchmark_id: summary for summary in summaries}

    assert by_id["cys_glucose_150C_Farmer1999"].supported is True
    assert by_id["cys_glucose_150C_Farmer1999"].benchmark_engine == "fast_observable"
    assert by_id["cys_glucose_150C_Farmer1999"].cantera_role == "diagnostic_reference_only"
    assert by_id["cys_glucose_150C_Farmer1999"].thermodynamic_gating_policy == "diagnostic_only"
    assert by_id["cys_glucose_150C_Farmer1999"].overall_status in {"pass", "scale-gap", "ranking-gap", "pass-no-ranking"}

    assert by_id["pea_isolate_40C_PratapSingh2021"].supported is True
    assert by_id["pea_isolate_40C_PratapSingh2021"].benchmark_engine == "matrix_intake_headspace"
    assert by_id["pea_isolate_40C_PratapSingh2021"].cantera_role == "not_authoritative"
    assert by_id["pea_isolate_40C_PratapSingh2021"].thermodynamic_gating_policy == "not_applicable"
    assert by_id["pea_isolate_40C_PratapSingh2021"].process_state == "ambient_slurry"
    assert by_id["pea_isolate_40C_PratapSingh2021"].ranking_contract_status == "pass"
    assert by_id["pea_isolate_40C_PratapSingh2021"].strict_ready is False
    assert by_id["pea_isolate_40C_PratapSingh2021"].overall_status in {"pass", "pass-no-ranking", "scale-gap", "ranking-gap"}

    assert by_id["soy_isolate_40C_PratapSingh2021"].supported is True
    assert by_id["soy_isolate_40C_PratapSingh2021"].benchmark_engine == "matrix_intake_headspace"
    assert by_id["soy_isolate_40C_PratapSingh2021"].cantera_role == "not_authoritative"
    assert by_id["soy_isolate_40C_PratapSingh2021"].process_state == "ambient_slurry"
    assert by_id["soy_isolate_40C_PratapSingh2021"].ranking_contract_status == "pass"
    assert by_id["soy_isolate_40C_PratapSingh2021"].strict_ready is False
    assert by_id["soy_isolate_40C_PratapSingh2021"].overall_status in {"pass", "pass-no-ranking", "scale-gap", "ranking-gap"}


def test_benchmark_summary_markdown_includes_gap_labels(benchmark_summary_matrix_scope):
    summaries = benchmark_summary_matrix_scope

    markdown = render_benchmark_summary_markdown(summaries)

    assert "Benchmark Summary" in markdown
    assert "cys_glucose_150C_Farmer1999" in markdown
    assert "pea_isolate_40C_PratapSingh2021" in markdown
    assert "soy_isolate_40C_PratapSingh2021" in markdown
    assert "matrix-only intake path is executable" in markdown
    assert "ambient_slurry" in markdown
    assert "Ranking Contract" in markdown
    assert "Strict Ready" in markdown
    assert "Cantera Role" in markdown
    assert "Thermo Policy" in markdown
    assert "Mean |log10 ratio|" in markdown
    assert "diagnostic_reference_only" in markdown
    assert "Chemistry Families" in markdown
    assert "Payload Roles" in markdown


def test_family_lane_validation_artifact_groups_benchmarks_by_family_and_lane(family_lane_validation_payload):
    payload = family_lane_validation_payload

    assert payload["summary"]["benchmark_count"] == 3
    assert payload["summary"]["lane_count"] >= 2
    assert any(row["execution_path"] == "free_precursor" for row in payload["lanes"])
    assert any(row["execution_path"] == "matrix_only" for row in payload["lanes"])
    assert any("benchmark_payload" in row["payload_roles"] for row in payload["families"])

    markdown = render_family_lane_validation_markdown(payload)
    assert "Family Lane Validation" in markdown
    assert "Lane Summary" in markdown
    assert "benchmark_payload" in markdown


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


def test_two_point_benchmark_surfaces_pass_without_ranking_when_scale_holds():
    evaluation = BenchmarkEvaluation(
        benchmark_id="two_point_benchmark",
        bench_file=ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=[
            CompoundComparison("cmp1", 10.0, 9.5, "cmp1", None, 1.0),
            CompoundComparison("cmp2", 20.0, 19.0, "cmp2", None, 1.0),
        ],
        pearson_r=None,
        mae_ppb=0.75,
    )

    summary = summarize_evaluation(evaluation, protein_type="free")

    assert summary.ranking_status == "n/a"
    assert summary.scale_status == "pass"
    assert summary.overall_status == "pass-no-ranking"
    assert summary.strict_ready is True