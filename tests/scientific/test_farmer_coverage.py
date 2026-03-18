import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark


def test_farmer_benchmark_has_full_species_coverage():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "cys_glucose_150C_Farmer1999.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0

    matched = {comparison.compound: comparison.matched_name for comparison in evaluation.comparisons}
    assert matched["2-methyl-3-furanthiol"] is not None
    assert matched["furfural"] is not None
    assert matched["pyrazine"] is not None