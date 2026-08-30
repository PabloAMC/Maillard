from copy import deepcopy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.pipeline import MaillardPipeline
from src.recommend import _canon


# Retargeted from cys_ribose_150C_Mottram1994 (quarantined 2026-08-26, see
# data/benchmarks/quarantined/README.md) to the surviving free_precursor cys+ribose panel
# benchmark. These tests assert monotonicity in temperature, time and precursor loading;
# they use the benchmark only as a well-formed free-precursor fixture and never compare
# against its measured values.
FREE_PRECURSOR_BENCH = ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json"
FURFURAL_CANON = _canon("O=Cc1ccco1")
MFT_CANON = _canon("Cc1occc1S")


def _evaluate_benchmark_variant(*, temp_c: float | None = None, time_minutes: float | None = None, scale: float = 1.0):
    bench = load_benchmark(FREE_PRECURSOR_BENCH)
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

    designer = MaillardPipeline(target_tag="meaty")
    return designer.evaluate_single(formulation, conditions)


def test_formulation_time_minutes_changes_predicted_ppb():
    result_short = _evaluate_benchmark_variant(temp_c=150.0, time_minutes=5.0)
    result_long = _evaluate_benchmark_variant(temp_c=150.0, time_minutes=60.0)

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


def test_benchmark_temperature_increase_raises_volatile_budget():
    """The total volatile budget and the temperature factor are monotonic in temperature.

    This is the claim that actually holds model-wide. It is asserted over a wide bracket
    (130 -> 170 C) precisely because it is not regime-dependent.
    """
    result_cool = _evaluate_benchmark_variant(temp_c=130.0, time_minutes=60.0)
    result_hot = _evaluate_benchmark_variant(temp_c=170.0, time_minutes=60.0)

    assert (
        result_hot.projection_metadata[FURFURAL_CANON]["projection_temperature_factor"]
        > result_cool.projection_metadata[FURFURAL_CANON]["projection_temperature_factor"]
    )
    assert (
        result_hot.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
        > result_cool.projection_metadata[FURFURAL_CANON]["total_volatile_budget_molar"]
    )


def test_benchmark_temperature_increase_raises_proxy_and_observable_outputs():
    """Per-species output rises monotonically with temperature.

    HISTORY (resolved 2026-08-27). This test used to be named
    `..._below_turnover` and was bracketed to 110 -> 140 C because individual predicted
    ppb values were NOT monotonic: for cysteine + pentose systems the model's furfural
    output peaked around 140-150 C and then declined (493 ppb at 130 C, 508 at 140, 506
    at 150, 489 at 170). The docstring flagged the turnover as an open science question --
    real chemistry (competing degradation and melanoidin routes above ~150 C) or an
    artefact of the budget-share arithmetic?

    It was the artefact. Two defects in the projection layer produced it, both fixed in
    the 2026-08-27 retune (audit remediation Part 1):

      1. The volatile budget was scaled by a severity sigmoid centred on 110 C with width
         18 C, which saturates by construction -- 0.966 at 170 C, 0.988 at 190 C. The
         budget therefore grew only ~1.11x from 150 C to 190 C while the model's own
         Arrhenius drive grew ~20x. A nearly-constant numerator divided by a still-growing
         softmax denominator is a peak.
      2. Boltzmann selectivity on `path_span` was applied twice -- once explicitly and
         once inside the pathway flux it multiplied by -- so the denominator grew at an
         effective temperature of T/2.19, sharpening the competition artificially fast.

    The budget's thermal dependence is now a first-order conversion extent under an
    apparent Arrhenius rate, and selectivity is applied once at the physical temperature.
    The furfural curve is monotone across the whole range (62 ppb at 110 C, 352 at 130,
    1688 at 150, 7154 at 170, 25770 at 190), so the bracket is restored to the full
    130 -> 170 C span the original test wanted to assert.
    """
    result_cool = _evaluate_benchmark_variant(temp_c=130.0, time_minutes=60.0)
    result_hot = _evaluate_benchmark_variant(temp_c=170.0, time_minutes=60.0)

    assert result_hot.predicted_proxy_ppb["furfural"] > result_cool.predicted_proxy_ppb["furfural"]
    assert result_hot.predicted_proxy_ppb["2-methyl-3-furanthiol"] > result_cool.predicted_proxy_ppb["2-methyl-3-furanthiol"]
    assert result_hot.predicted_ppb["furfural"] > result_cool.predicted_ppb["furfural"]
    assert result_hot.predicted_ppb["2-methyl-3-furanthiol"] > result_cool.predicted_ppb["2-methyl-3-furanthiol"]


def test_furfural_output_is_monotonic_in_temperature_across_the_full_range():
    """Regression guard for the resolved furfural turnover (2026-08-27).

    The old artefact peaked near 145-150 C. Sweep the full 110 -> 190 C range and require
    strict monotonicity, so a future change to the budget scale or the allocation cannot
    quietly reintroduce it.
    """
    temperatures = [110.0, 130.0, 150.0, 170.0, 190.0]
    furfural = [
        _evaluate_benchmark_variant(temp_c=temp, time_minutes=60.0).predicted_ppb["furfural"]
        for temp in temperatures
    ]

    for cooler, hotter, low_t, high_t in zip(furfural, furfural[1:], temperatures, temperatures[1:]):
        assert hotter > cooler, (
            f"Furfural fell from {cooler:.1f} ppb at {low_t:.0f} C to {hotter:.1f} ppb at "
            f"{high_t:.0f} C. The projection-layer turnover artefact has returned -- check "
            "the budget's temperature dependence and the Boltzmann de-duplication in "
            "src/projection.py / src/recommend.py before touching this test."
        )


def test_benchmark_precursor_loading_scales_proxy_and_observable_outputs():
    result_low = _evaluate_benchmark_variant(temp_c=150.0, time_minutes=60.0, scale=0.5)
    result_high = _evaluate_benchmark_variant(temp_c=150.0, time_minutes=60.0, scale=4.0)

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