import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.inverse_design import InverseDesigner


def test_free_amino_acid_benchmark_does_not_inject_lipid_oxidation_products():
    bench = load_benchmark(ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json")
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)

    result = InverseDesigner(target_tag="meaty").evaluate_single(formulation, conditions)

    assert result.predicted_ppb.get("Hexanal") is None
    assert result.predicted_ppb.get("Nonanal") is None
    assert result.predicted_ppb.get("2-Pentylfuran") is None


def test_free_amino_acid_benchmark_predicted_ppb_excludes_input_precursors():
    bench = load_benchmark(ROOT / "data" / "benchmarks" / "cys_ribose_150C_Mottram1994.json")
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)

    result = InverseDesigner(target_tag="meaty").evaluate_single(formulation, conditions)

    assert result.predicted_ppb.get("D-Ribose") is None
    assert result.predicted_ppb.get("L-Cysteine") is None