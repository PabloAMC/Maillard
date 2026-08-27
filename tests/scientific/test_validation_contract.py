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
    assert "Pearson, ratio, and log-scale error thresholds" in contract.quantitative_replication
    assert "ranking recipes or interventions" in contract.formulation_utility
    assert "FAST observable projection" in contract.benchmark_policy
    assert "Cantera remains a diagnostic reference lane" in contract.benchmark_policy
    assert contract.thresholds.min_matched_for_ranking == 3
    assert contract.thresholds.ranking_threshold == 0.85
    assert contract.thresholds.free_aa_ratio_threshold == 1.5
    assert contract.thresholds.matrix_ratio_threshold == 2.0
    assert contract.thresholds.free_aa_mean_abs_log10_error_threshold == 0.10
    assert contract.thresholds.matrix_mean_abs_log10_error_threshold == 0.12


def test_validation_contract_registers_execution_policy_for_fast_and_matrix_paths():
    contract = DEFAULT_VALIDATION_CONTRACT

    free_policy = contract.policy_for_execution_path("free_precursor")
    assert free_policy.benchmark_engine == "fast_observable"
    assert free_policy.comparator_signal == "predicted_ppb"
    assert free_policy.cantera_role == "diagnostic_reference_only"
    assert free_policy.target_snapshot_policy == "included"
    assert free_policy.thermodynamic_gating_policy == "diagnostic_only"

    matrix_policy = contract.policy_for_execution_path("matrix_only")
    assert matrix_policy.benchmark_engine == "matrix_intake_headspace"
    assert matrix_policy.comparator_signal == "predicted_ppb"
    assert matrix_policy.cantera_role == "not_authoritative"
    assert matrix_policy.target_snapshot_policy == "excluded"
    assert matrix_policy.thermodynamic_gating_policy == "not_applicable"


def test_strict_gate_eligibility_is_centralized_in_validation_contract():
    comparison_rows = [
        CompoundComparison("cmp1", 10.0, 10.0, "cmp1", None, 1.0),
        CompoundComparison("cmp2", 20.0, 19.0, "cmp2", None, 1.0),
        CompoundComparison("cmp3", 30.0, 29.0, "cmp3", None, 1.0),
    ]

    free_eval = BenchmarkEvaluation(
        benchmark_id="free_eval",
        bench_file=ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
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
    assert matrix_summary.cantera_role == "not_authoritative"
    assert any("strict release gate" in issue for issue in matrix_summary.blocking_issues)


def test_benchmark_specific_scale_thresholds_override_global_defaults():
    # Retargeted twice: originally cys_glucose_150C_Farmer1999, then Parker 2012 (both
    # quarantined 2026-08-26 for unlocatable sources; Farmer since deleted). It now runs on
    # cys_ribose_140C_Hofmann1998, the surviving verified free-precursor benchmark.
    #
    # This is a sharper exercise of the override than either predecessor. The global default
    # for protein_type="free" is max_ratio 1.5 (BenchmarkThresholds.free_aa_ratio_threshold);
    # this benchmark's file pins 1.45. The observed ratio below (1.48) sits BETWEEN the two,
    # so it passes under the global default and fails only if the benchmark-specific value is
    # actually being read -- the override is the whole difference between pass and fail.
    #
    # The mean |log10 ratio| is held out of the way on purpose: mean(log10(1.48), 0) = 0.085,
    # under both the file's 0.09 and the global 0.10, so the failure below isolates the ratio.
    comparison_rows = [
        CompoundComparison("cmp1", 100.0, 148.0, "cmp1", None, 1.0),
        CompoundComparison("cmp2", 100.0, 100.0, "cmp2", None, 1.0),
    ]

    evaluation = BenchmarkEvaluation(
        benchmark_id="override_eval",
        bench_file=ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=comparison_rows,
        pearson_r=None,
        mae_ppb=15.0,
    )

    summary = summarize_evaluation(evaluation, protein_type="free")

    assert summary.max_ratio == 1.48
    assert summary.mean_abs_log10_error is not None
    assert summary.mean_abs_log10_error < DEFAULT_VALIDATION_CONTRACT.thresholds.free_aa_mean_abs_log10_error_threshold
    # Would pass under the global default...
    assert summary.max_ratio < DEFAULT_VALIDATION_CONTRACT.thresholds.ratio_threshold_for("free")
    # ...but the benchmark-specific 1.45 is tighter, and it is what gets enforced.
    assert summary.scale_status == "fail"
    assert any("max ratio 1.480 > 1.45" in issue for issue in summary.blocking_issues)