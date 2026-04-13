from pathlib import Path

from src.benchmark_validation import evaluate_benchmark
from src.safety import predict_furosine


ROOT = Path(__file__).resolve().parents[2]


def _comparison_map(evaluation):
    return {comparison.compound: comparison for comparison in evaluation.comparisons}


def test_hvp_xylose_benchmarks_are_executable_and_rank_wheat_above_spi():
    spi_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "spi_hvp_xylose_120C_PMC9905368.json")
    wheat_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "wheat_gluten_hvp_xylose_120C_PMC9905368.json")

    assert spi_eval.supported, spi_eval.reason
    assert wheat_eval.supported, wheat_eval.reason
    assert spi_eval.coverage == 1.0
    assert wheat_eval.coverage == 1.0

    spi_mft = _comparison_map(spi_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    wheat_mft = _comparison_map(wheat_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    assert wheat_mft > spi_mft
    assert max(comparison.ratio for comparison in spi_eval.comparisons) <= 1.5
    assert max(comparison.ratio for comparison in wheat_eval.comparisons) <= 1.5


def test_fast_acrylamide_extrusion_benchmark_is_locally_calibrated():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "acrylamide_spi_extrusion_130C_ACSRef3.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    comparison = evaluation.comparisons[0]
    assert comparison.compound == "acrylamide"
    assert comparison.ratio <= 1.35


def test_cml_cel_and_furosine_benchmarks_execute_with_expected_directionality():
    age_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "cml_cel_commercial_pbma_Foods2023.json")
    furosine_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json")

    assert age_eval.supported, age_eval.reason
    assert age_eval.coverage == 1.0
    assert furosine_eval.supported, furosine_eval.reason
    assert furosine_eval.coverage == 1.0

    age_by_compound = _comparison_map(age_eval)
    assert age_by_compound["Nε-(Carboxyethyl)lysine (CEL)"].predicted_ppb > age_by_compound["Nε-(Carboxymethyl)lysine (CML)"].predicted_ppb

    peak = predict_furosine(140.0, 20.0, lysine_mM=35.0, reducing_sugar_mM=35.0, protein_type="pea_iso", water_activity=0.55)
    severe = predict_furosine(165.0, 20.0, lysine_mM=35.0, reducing_sugar_mM=35.0, protein_type="pea_iso", water_activity=0.55)
    assert peak > severe
    assert max(comparison.ratio for comparison in furosine_eval.comparisons) <= 2.0


def test_resconi_identity_gap_subset_benchmark_executes_on_supported_marker_subset():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "resconi_2023_pbma_beef_identity_benchmark.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "furfural"
    assert comparison.matched_name is not None
    assert comparison.ratio <= 1.5