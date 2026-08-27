import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    evaluate_benchmark,
    get_matrix_ranking_contract,
    snapshot_all_benchmark_targets,
    snapshot_benchmark_targets,
    summarize_evaluation,
)


MATRIX_ONLY_BENCHMARKS = {
    "pea_iso": ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    "soy_iso": ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
}
TRIKUSUMA_BENCHMARK = ROOT / "data" / "benchmarks" / "pea_isolate_uht_140C_Trikusuma2019.json"


@pytest.mark.parametrize("bench_file", MATRIX_ONLY_BENCHMARKS.values())
def test_matrix_only_benchmark_is_executable_with_full_coverage(bench_file):

    evaluation = evaluate_benchmark(bench_file)

    assert evaluation.supported is True
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 3


@pytest.mark.parametrize(
    ("protein_type", "bench_file", "expected_limits"),
    [
        (
            "pea_iso",
            MATRIX_ONLY_BENCHMARKS["pea_iso"],
            {"hexanal": 1.25, "2-pentylfuran": 1.2, "hexanol": 1.35},
        ),
        (
            "soy_iso",
            MATRIX_ONLY_BENCHMARKS["soy_iso"],
            {"hexanal": 1.1, "2-pentylfuran": 1.05, "hexanol": 1.05},
        ),
    ],
)
def test_matrix_only_benchmark_preserves_measured_ordering_without_entering_strict_gate(
    protein_type,
    bench_file,
    expected_limits,
):

    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type=protein_type)
    ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}
    predicted = {comparison.compound: comparison.predicted_ppb for comparison in evaluation.comparisons}

    assert predicted["2-pentylfuran"] > predicted["hexanal"] > predicted["hexanol"]
    assert ratios["hexanal"] <= expected_limits["hexanal"]
    assert ratios["2-pentylfuran"] <= expected_limits["2-pentylfuran"]
    assert ratios["hexanol"] <= expected_limits["hexanol"]
    assert summary.strict_ready is False


def test_matrix_only_benchmark_is_deliberately_excluded_from_target_snapshots():
    bench_file = MATRIX_ONLY_BENCHMARKS["pea_iso"]
    soy_bench_file = MATRIX_ONLY_BENCHMARKS["soy_iso"]

    assert snapshot_benchmark_targets(bench_file) == []
    assert snapshot_benchmark_targets(soy_bench_file) == []

    rows = snapshot_all_benchmark_targets([
        bench_file,
        soy_bench_file,
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ])

    assert rows
    assert {row.benchmark_id for row in rows} == {"cys_ribose_140C_Hofmann1998"}


def test_matrix_only_benchmark_exposes_ranking_contract_and_calibration_metadata():
    bench_file = MATRIX_ONLY_BENCHMARKS["soy_iso"]

    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type="soy_iso")
    contract = get_matrix_ranking_contract(bench_file)

    assert contract["process_state"] == "ambient_slurry"
    assert contract["calibration_mode"] == "compound_specific_headspace"
    assert [item["name"] for item in contract["observable_targets"]] == [
        "2-pentylfuran",
        "hexanal",
        "hexanol",
    ]

    hexanal_meta = evaluation.projection_metadata["Hexanal"]
    assert hexanal_meta["calibration_evidence_strength"] == "literature_anchored"
    assert hexanal_meta["calibration_fallback_mode"] == "compound_specific"
    assert hexanal_meta["process_state"] == "ambient_slurry"
    assert hexanal_meta["evidence_state"] == "externally_benchmarked"
    assert hexanal_meta["target_class"] == "adverse_lipid_markers"
    assert summary.ranking_contract_status == "pass"
    assert summary.adverse_markers == ["2-pentylfuran", "hexanal", "hexanol"]


def test_trikusuma_heated_pea_matrix_benchmark_is_quantitatively_supported():
    evaluation = evaluate_benchmark(TRIKUSUMA_BENCHMARK)
    summary = summarize_evaluation(evaluation, protein_type="pea_iso")
    predicted = {comparison.compound: comparison.predicted_ppb for comparison in evaluation.comparisons}
    ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}

    assert predicted["hexanal"] > predicted["2-pentylfuran"] > predicted["nonanal"]
    assert ratios["hexanal"] <= 1.05
    assert ratios["2-pentylfuran"] <= 1.05
    assert ratios["nonanal"] <= 1.05
    assert summary.process_state == "heated_matrix"
    assert summary.ranking_contract_status == "pass"
    assert summary.overall_status == "pass"
    assert summary.strict_ready is False