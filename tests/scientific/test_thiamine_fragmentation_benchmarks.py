import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark  # noqa: E402


def test_hofmann_1996_thiamine_fragmentation_benchmark_is_executable_and_inside_strict_band():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "thiamine_cys_ribose_100C_Hofmann1996.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "2-Methyl-3-furanthiol (MFT)"
    assert comparison.ratio <= 1.5


def test_cerny_2008_thiamine_fragmentation_benchmark_is_executable_and_inside_strict_band():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "thiamine_cys_xylose_145C_Cerny2008.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert evaluation.reference_signal_origin == "reference_volatiles"
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "2-Methyl-3-furanthiol (MFT)"
    assert comparison.ratio <= 1.5