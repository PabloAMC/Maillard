from copy import deepcopy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.inverse_design import InverseDesigner


def test_formulation_time_minutes_changes_predicted_ppb():
    bench = load_benchmark(ROOT / "data" / "benchmarks" / "cys_ribose_150C_Mottram1994.json")
    formulation_short = benchmark_to_formulation(bench)
    formulation_long = deepcopy(formulation_short)
    conditions = benchmark_to_conditions(bench)

    formulation_short["time_minutes"] = 5.0
    formulation_long["time_minutes"] = 60.0

    designer = InverseDesigner(target_tag="meaty")
    result_short = designer.evaluate_single(formulation_short, conditions)
    result_long = designer.evaluate_single(formulation_long, conditions)

    assert result_short.predicted_ppb["furfural"] != result_long.predicted_ppb["furfural"]
    assert result_short.predicted_ppb["2-methyl-3-furanthiol"] != result_long.predicted_ppb["2-methyl-3-furanthiol"]