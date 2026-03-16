import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import BenchmarkEvaluation, CompoundComparison, summarize_evaluation
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT


def test_all_benchmark_json_files_carry_explicit_validation_metadata():
    benchmark_dir = ROOT / "data" / "benchmarks"
    benchmark_files = sorted(path for path in benchmark_dir.glob("*.json"))

    assert benchmark_files
    for bench_file in benchmark_files:
        with open(bench_file, "r", encoding="utf-8") as handle:
            bench = json.load(handle)
        metadata = bench.get("metadata")
        assert metadata is not None, f"{bench_file.name} is missing metadata"
        assert metadata.get("tier"), f"{bench_file.name} is missing metadata.tier"
        assert metadata.get("family"), f"{bench_file.name} is missing metadata.family"
        assert metadata.get("execution_path"), f"{bench_file.name} is missing metadata.execution_path"


def test_validation_contract_centralizes_replication_axes_and_thresholds():
    contract = DEFAULT_VALIDATION_CONTRACT

    assert "coverage, ranking, and scale" in contract.replication_meaning
    assert "ordering or trend" in contract.directional_validity
    assert "Pearson and ratio thresholds" in contract.quantitative_replication
    assert "ranking recipes or interventions" in contract.formulation_utility
    assert contract.thresholds.min_matched_for_ranking == 3
    assert contract.thresholds.ranking_threshold == 0.85
    assert contract.thresholds.free_aa_ratio_threshold == 1.5
    assert contract.thresholds.matrix_ratio_threshold == 2.0


def test_strict_gate_eligibility_is_centralized_in_validation_contract():
    comparison_rows = [
        CompoundComparison("cmp1", 10.0, 10.0, "cmp1", None, 1.0),
        CompoundComparison("cmp2", 20.0, 19.0, "cmp2", None, 1.0),
        CompoundComparison("cmp3", 30.0, 29.0, "cmp3", None, 1.0),
    ]

    free_eval = BenchmarkEvaluation(
        benchmark_id="free_eval",
        bench_file=ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=comparison_rows,
        pearson_r=0.95,
        mae_ppb=1.0,
    )
    matrix_eval = BenchmarkEvaluation(
        benchmark_id="matrix_eval",
        bench_file=ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=comparison_rows,
        pearson_r=0.95,
        mae_ppb=1.0,
    )

    assert summarize_evaluation(free_eval, protein_type="free").strict_ready is True
    matrix_summary = summarize_evaluation(matrix_eval, protein_type="pea_iso")
    assert matrix_summary.strict_ready is False
    assert any("strict release gate" in issue for issue in matrix_summary.blocking_issues)