from copy import deepcopy
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.recommend import _canon


MOTTRAM_BENCH = ROOT / "data" / "benchmarks" / "cys_ribose_150C_Mottram1994.json"
FURFURAL_CANON = _canon("O=Cc1ccco1")
MFT_CANON = _canon("Cc1occc1S")


def _evaluate_benchmark_variant(
    designer,
    *,
    temp_c: float | None = None,
    time_minutes: float | None = None,
    scale: float = 1.0,
):
    bench = load_benchmark(MOTTRAM_BENCH)
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)

    if temp_c is not None:
        conditions.temperature_celsius = temp_c
        formulation["temp"] = temp_c
    if time_minutes is not None:
        formulation["time_minutes"] = time_minutes
    if scale != 1.0:
        formulation["molar_ratios"] = {
            key: float(value) * scale for key, value in formulation["molar_ratios"].items()
        }

    return designer.evaluate_single(formulation, conditions)


@pytest.mark.slow
def test_formulation_time_minutes_changes_predicted_ppb(pipeline_meaty):
    result_short = _evaluate_benchmark_variant(pipeline_meaty, temp_c=150.0, time_minutes=5.0)
    result_long = _evaluate_benchmark_variant(pipeline_meaty, temp_c=150.0, time_minutes=60.0)

    assert result_long.predicted_ppb["furfural"] > result_short.predicted_ppb["furfural"]
    assert result_long.predicted_ppb["2-methyl-3-furanthiol"] > result_short.predicted_ppb["2-methyl-3-furanthiol"]
    assert (
        result_long.projection_metadata[FURFURAL_CANON]["projection_time_factor"]
        > result_short.projection_metadata[FURFURAL_CANON]["projection_time_factor"]
    )
    assert (
        result_long.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
        > result_short.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
    )


@pytest.mark.slow
def test_benchmark_temperature_increase_raises_proxy_and_observable_outputs(pipeline_meaty):
    result_cool = _evaluate_benchmark_variant(pipeline_meaty, temp_c=130.0, time_minutes=60.0)
    result_hot = _evaluate_benchmark_variant(pipeline_meaty, temp_c=170.0, time_minutes=60.0)

    assert result_hot.predicted_proxy_ppb["furfural"] > result_cool.predicted_proxy_ppb["furfural"]
    assert result_hot.predicted_proxy_ppb["2-methyl-3-furanthiol"] > result_cool.predicted_proxy_ppb["2-methyl-3-furanthiol"]
    assert result_hot.predicted_ppb["furfural"] > result_cool.predicted_ppb["furfural"]
    assert result_hot.predicted_ppb["2-methyl-3-furanthiol"] > result_cool.predicted_ppb["2-methyl-3-furanthiol"]
    assert (
        result_hot.projection_metadata[FURFURAL_CANON]["projection_temperature_factor"]
        > result_cool.projection_metadata[FURFURAL_CANON]["projection_temperature_factor"]
    )
    assert (
        result_hot.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
        > result_cool.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
    )


@pytest.mark.slow
def test_benchmark_precursor_loading_scales_proxy_and_observable_outputs(pipeline_meaty):
    result_low = _evaluate_benchmark_variant(pipeline_meaty, temp_c=150.0, time_minutes=60.0, scale=0.5)
    result_high = _evaluate_benchmark_variant(pipeline_meaty, temp_c=150.0, time_minutes=60.0, scale=4.0)

    assert result_high.predicted_proxy_ppb["furfural"] > result_low.predicted_proxy_ppb["furfural"]
    assert result_high.predicted_proxy_ppb["2-methyl-3-furanthiol"] > result_low.predicted_proxy_ppb["2-methyl-3-furanthiol"]
    assert result_high.predicted_ppb["furfural"] > result_low.predicted_ppb["furfural"]
    assert result_high.predicted_ppb["2-methyl-3-furanthiol"] > result_low.predicted_ppb["2-methyl-3-furanthiol"]
    assert (
        result_high.projection_metadata[FURFURAL_CANON]["limiting_precursor_molar"]
        > result_low.projection_metadata[FURFURAL_CANON]["limiting_precursor_molar"]
    )
    assert (
        result_high.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
        > result_low.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
    )