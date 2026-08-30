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
    """A PRIMARY free_precursor file is strict-gate eligible; other tiers and paths are not.

    RETARGETED 2026-08-27 (Wave S2c). The free-precursor fixture was
    `cys_ribose_140C_Hofmann1998.json`; it is now
    `acrylamide_spi_extrusion_130C_ACSRef3.json`. CAUSE: Hofmann's `metadata.tier` was demoted
    PRIMARY -> REFERENCE after Wave S2b showed its MFT 342 / FFT 200 ppb are a repo-internal
    derivation (interior points of two invented mol % bands in
    data/benchmarks/maillard_validation_benchmarks.md section 1.3) rather than a measurement
    from 10.1021/jf9705983. `strict_gate_tiers` is ("PRIMARY",), so the file stopped being
    strict-gate eligible and could no longer serve as the POSITIVE fixture here.
    THE LOSS IS TURNED INTO A GUARD rather than deleted: the demotion is now asserted
    explicitly at the bottom of this test, so a future edit that quietly restores
    `tier: PRIMARY` on that benchmark fails here instead of silently re-admitting a
    fabricated anchor to the strict gate.
    """
    comparison_rows = [
        CompoundComparison("cmp1", 10.0, 10.0, "cmp1", None, 1.0),
        CompoundComparison("cmp2", 20.0, 19.0, "cmp2", None, 1.0),
        CompoundComparison("cmp3", 30.0, 29.0, "cmp3", None, 1.0),
    ]

    free_eval = BenchmarkEvaluation(
        benchmark_id="free_eval",
        bench_file=ROOT / "data" / "benchmarks" / "acrylamide_spi_extrusion_130C_ACSRef3.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=comparison_rows,
        pearson_r=0.95,
        mae_ppb=1.0,
    )
    demoted_eval = BenchmarkEvaluation(
        benchmark_id="demoted_eval",
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

    # Wave S2c: the TIER half of the gate, asserted on the file the demotion was made for.
    # Everything else about this benchmark is identical to `free_eval` -- same
    # execution_path (free_precursor), same synthetic rows, same protein_type -- so the ONLY
    # thing standing between it and strict_ready is `metadata.tier`.
    demoted_summary = summarize_evaluation(demoted_eval, protein_type="free")
    assert demoted_summary.tier == "REFERENCE"
    assert demoted_summary.strict_ready is False, (
        "cys_ribose_140C_Hofmann1998 is strict-gate eligible again. Its two 'measured' "
        "values are a repo-internal derivation with no verifiable source (Wave S2b); it was "
        "demoted PRIMARY -> REFERENCE by Wave S2c and must not be re-promoted until it is "
        "rebuilt from the paper's own table in native mol %. See its metadata.tier_history "
        "and tasks/audit_remediation.md '## Wave S2c'."
    )
    assert not DEFAULT_VALIDATION_CONTRACT.is_strict_gate_eligible(
        tier="REFERENCE", execution_path="free_precursor"
    )


def test_benchmark_specific_scale_thresholds_override_global_defaults():
    # Retargeted three times: originally cys_glucose_150C_Farmer1999, then Parker 2012 (both
    # quarantined 2026-08-26 for unlocatable sources; Farmer since deleted), then
    # cys_ribose_140C_Hofmann1998.
    #
    # RETARGETED AGAIN 2026-08-27 (Wave S2c), TO furosine_extrusion_crossover_140C_
    # RamirezJimenez2000, AND THE TEST'S DIRECTION IS FLIPPED. Read why, because the reason is
    # a finding rather than bookkeeping. Wave S2c RETIRED the 1.45x / 0.09 dex contract on
    # cys_ribose_140C_Hofmann1998: Wave S2b had shown the benchmark's MFT 342 / FFT 200 ppb to
    # be a repo-internal derivation (interior points of two invented mol % bands in
    # data/benchmarks/maillard_validation_benchmarks.md section 1.3), and the contract was
    # ~1.7x TIGHTER than the 2.5x spread of the band its own target came from. Nothing looser
    # was invented to replace it; the file now inherits the global default.
    #
    # THE CONSEQUENCE FOR THIS TEST: after that retirement, NO free-precursor benchmark in
    # data/benchmarks/ carries a contract TIGHTER than the global default any more. The
    # survivors are acrylamide 1.5/0.20 (equal on ratio), cml_cel 1.8/0.25, furosine 2.0/0.30
    # and Bolton 3.0/0.48 -- every one of them looser. So the override can no longer be
    # exercised in the tighten-then-fail direction against a shipped file, and manufacturing a
    # fixture file to preserve the old shape would test a fixture rather than the panel.
    # The override is therefore exercised in the LOOSEN-then-pass direction instead, which
    # isolates exactly the same mechanism: the observed ratio 1.70 sits BETWEEN the global
    # default (1.5, would FAIL) and this file's own 2.0 (PASSES), so the assertion below can
    # only hold if the benchmark-specific value is actually being read. The MALE does the same
    # work independently -- mean(log10(1.70), 0) = 0.1152 is OVER the global 0.10 and UNDER
    # the file's 0.30 -- so both criteria are on the correct side of the override at once.
    # The fail direction is kept too, at ratio 2.50 against the file's own 2.0, so this test
    # still proves the override can reject and has not simply become permissive.
    comparison_rows = [
        CompoundComparison("cmp1", 100.0, 170.0, "cmp1", None, 1.0),
        CompoundComparison("cmp2", 100.0, 100.0, "cmp2", None, 1.0),
    ]

    evaluation = BenchmarkEvaluation(
        benchmark_id="override_eval",
        bench_file=ROOT
        / "data"
        / "benchmarks"
        / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=comparison_rows,
        pearson_r=None,
        mae_ppb=15.0,
    )

    summary = summarize_evaluation(evaluation, protein_type="free")

    assert summary.max_ratio == 1.70
    assert summary.mean_abs_log10_error is not None
    # Would FAIL on both criteria under the global default...
    assert summary.max_ratio > DEFAULT_VALIDATION_CONTRACT.thresholds.ratio_threshold_for("free")
    assert (
        summary.mean_abs_log10_error
        > DEFAULT_VALIDATION_CONTRACT.thresholds.free_aa_mean_abs_log10_error_threshold
    )
    # ...but the benchmark-specific 2.0 / 0.30 is looser, and it is what gets enforced.
    assert summary.scale_status == "pass"
    assert summary.blocking_issues == []

    # And the override still REJECTS: 2.50 is outside this file's own 2.0.
    rejecting = BenchmarkEvaluation(
        benchmark_id="override_eval_reject",
        bench_file=ROOT
        / "data"
        / "benchmarks"
        / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json",
        supported=True,
        reason=None,
        predicted_ppb={},
        comparisons=[
            CompoundComparison("cmp1", 100.0, 250.0, "cmp1", None, 1.0),
            CompoundComparison("cmp2", 100.0, 100.0, "cmp2", None, 1.0),
        ],
        pearson_r=None,
        mae_ppb=15.0,
    )
    rejected = summarize_evaluation(rejecting, protein_type="free")
    assert rejected.scale_status == "fail"
    assert any("max ratio 2.500 > 2.00" in issue for issue in rejected.blocking_issues)

    # The retired Hofmann contract must stay retired: no scale_thresholds block, and no
    # replacement quietly installed. If someone re-adds one it has to come from a real
    # measurement's reported precision, not from a number this repository picked.
    hofmann = json.loads(
        (ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json").read_text(
            encoding="utf-8"
        )
    )
    contract = hofmann.get("validation_contract", {})
    assert "scale_thresholds" not in contract, (
        "cys_ribose_140C_Hofmann1998 has a scale contract again. It was RETIRED by Wave S2c "
        "because it was ~1.7x tighter than the spread of the invented band its own target "
        "value was interpolated from. Do not re-add one until the benchmark is rebuilt from "
        "the paper's table in native mol %."
    )
    assert contract["RETIRED"]["retired_thresholds"] == {
        "max_ratio": 1.45,
        "mean_abs_log10_error": 0.09,
    }