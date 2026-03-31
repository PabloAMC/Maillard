import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_mottram_benchmark_has_full_species_coverage(mottram_evaluation):
    evaluation = mottram_evaluation

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0

    matched = {comparison.compound: comparison.matched_name for comparison in evaluation.comparisons}
    assert matched["2-methyl-3-furanthiol"] is not None
    assert matched["bis(2-methyl-3-furyl) disulfide"] is not None
    assert matched["furfural"] is not None