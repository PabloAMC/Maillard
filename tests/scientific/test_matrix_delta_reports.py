import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (  # noqa: E402
    build_matrix_benchmark_deltas,
    evaluate_benchmark,
    snapshot_benchmark_targets,
    summarize_evaluation,
)
from src.presentation import render_matrix_benchmark_deltas_markdown


PEA_MEATY = ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json"
SOY_MEATY = ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_Internal2026.json"


@pytest.mark.slow
def test_matrix_precursor_augmented_benchmark_is_reproducible_but_not_strict_ready():
    evaluation = evaluate_benchmark(PEA_MEATY)
    summary = summarize_evaluation(evaluation, protein_type="pea_iso")

    assert evaluation.supported is True
    assert evaluation.reference_signal_origin == "reference_volatiles"
    assert summary.execution_path == "matrix_precursor_augmented"
    assert summary.ranking_contract_status == "pass"
    assert summary.strict_ready is False
    assert summary.reference_signal_origin == "reference_volatiles"


@pytest.mark.slow
def test_matrix_precursor_augmented_benchmark_exposes_target_snapshots():
    rows = snapshot_benchmark_targets(SOY_MEATY)
    names = {row.target_name for row in rows}

    assert rows
    assert "2-Furfurylthiol (FFT)" in names
    assert "2-Methyl-3-furanthiol (MFT)" in names


def test_matrix_benchmark_delta_report_covers_matrix_only_and_meaty_candidates(matrix_benchmark_delta_rows):
    rows = matrix_benchmark_delta_rows
    markdown = render_matrix_benchmark_deltas_markdown(rows)

    benchmark_ids = {row.benchmark_id for row in rows}
    assert "pea_isolate_40C_PratapSingh2021" in benchmark_ids
    assert "soy_isolate_40C_PratapSingh2021" in benchmark_ids
    assert "pea_isolate_ribose_cysteine_100C_45min_Internal2026" in benchmark_ids
    assert "soy_isolate_ribose_cysteine_100C_45min_Internal2026" in benchmark_ids
    assert "Matrix Benchmark Deltas" in markdown
    assert "reference_volatiles" in markdown
    assert "matrix_precursor_augmented" in markdown
