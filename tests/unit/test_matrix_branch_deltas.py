from pathlib import Path

from src.benchmark_validation import (
    MatrixBenchmarkDelta,
    MatrixBenchmarkEvidence,
    compare_matrix_benchmark_delta_sets,
    render_matrix_branch_deltas_markdown,
)


def test_compare_matrix_benchmark_delta_sets_reports_added_and_removed_rows():
    current_rows = [
        MatrixBenchmarkDelta(
            benchmark_id="pea_meaty",
            bench_file=Path("pea_meaty.json"),
            protein_type="pea_iso",
            execution_path="matrix_precursor_augmented",
            process_state="aqueous_pre_extrusion_model",
            reference_signal_origin="reference_volatiles",
            ranking_contract_status="pass",
            compound="2-furfurylthiol",
            role="desirable_marker",
            reference_ppb=0.3,
            predicted_ppb=0.3,
            abs_delta_ppb=0.0,
            pct_delta=0.0,
            ratio=1.0,
            calibration_source="class_fallback",
            calibration_evidence_strength="heuristic",
            calibration_fallback_mode="class_level",
        )
    ]
    baseline_rows = [
        MatrixBenchmarkDelta(
            benchmark_id="pea_off",
            bench_file=Path("pea_off.json"),
            protein_type="pea_iso",
            execution_path="matrix_only",
            process_state="ambient_slurry",
            reference_signal_origin="measured_volatiles",
            ranking_contract_status="pass",
            compound="hexanal",
            role="adverse_marker",
            reference_ppb=260.0,
            predicted_ppb=260.0,
            abs_delta_ppb=0.0,
            pct_delta=0.0,
            ratio=1.0,
            calibration_source="literature",
            calibration_evidence_strength="literature_anchored",
            calibration_fallback_mode="compound_specific",
        )
    ]
    current_evidence = [
        MatrixBenchmarkEvidence(
            benchmark_id="pea_meaty",
            bench_file=Path("pea_meaty.json"),
            protein_type="pea_iso",
            execution_path="matrix_precursor_augmented",
            process_state="aqueous_pre_extrusion_model",
            reference_signal_origin="reference_volatiles",
            source_origin="internal_reproducibility_candidate",
            source_reference="internal_reproducibility_candidate:docker_frozen_model_candidate",
            target_profile="meaty_positive",
            external_data_status="internal_reference_only",
            promotable=False,
            promotion_blocker="missing external quantitative matrix evidence for meaty-positive targets",
        )
    ]
    baseline_evidence = [
        MatrixBenchmarkEvidence(
            benchmark_id="pea_off",
            bench_file=Path("pea_off.json"),
            protein_type="pea_iso",
            execution_path="matrix_only",
            process_state="ambient_slurry",
            reference_signal_origin="measured_volatiles",
            source_origin="external_literature",
            source_reference="10.3390/molecules26134104",
            target_profile="adverse_only",
            external_data_status="external_quantitative",
            promotable=False,
            promotion_blocker="benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present",
        )
    ]

    rows = compare_matrix_benchmark_delta_sets(
        current_rows,
        baseline_rows,
        current_evidence=current_evidence,
        baseline_evidence=baseline_evidence,
    )
    markdown = render_matrix_branch_deltas_markdown(rows, base_ref="main")

    assert {row.change_type for row in rows} == {"added", "removed"}
    assert any(row.benchmark_id == "pea_meaty" and row.change_type == "added" for row in rows)
    assert any(row.benchmark_id == "pea_off" and row.change_type == "removed" for row in rows)
    assert "Matrix Benchmark Branch Comparison vs main" in markdown