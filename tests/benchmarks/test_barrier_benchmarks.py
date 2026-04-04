import pytest
from scipy.stats import spearmanr

from src.authority_benchmark_data import (
    load_phase33_barrier_benchmarks,
    load_phase35_double_hybrid_benchmarks,
)

from tests.benchmarks._lane_policy import (
    HAS_PHASE33_DATASET,
    HAS_PHASE35_DATASET,
    PHASE33_SKIP_REASON,
    PHASE35_SKIP_REASON,
)

# These benchmarks represent the research frontier (Phase 3). 
# Most tests are currently skipped as they depend on ongoing DFT batch runs.

REACTION_FAMILIES = [
    ("amadori_rearrangement", 20.0),
    ("enolisation_1_2", 20.0),
    ("enolisation_2_3", 25.0),
    ("strecker_decarboxylation", 24.0),
    ("fft_formation", 25.0),
    ("retro_aldol", 23.0),
    ("dha_elimination", 27.0),
    ("trapping_hexanal", 18.0),
    ("pyrazine_aminoketone", 28.0),
]

PHASE33_BENCHMARKS = load_phase33_barrier_benchmarks() if HAS_PHASE33_DATASET else {}
PHASE35_BENCHMARKS = load_phase35_double_hybrid_benchmarks() if HAS_PHASE35_DATASET else {}

@pytest.mark.slow
@pytest.mark.optional_dft_authority_lane
@pytest.mark.parametrize("family, expected_best", REACTION_FAMILIES)
@pytest.mark.skipif(not HAS_PHASE33_DATASET, reason=PHASE33_SKIP_REASON)
def test_barrier_benchmark_wb97mv(family, expected_best, literature_barriers_dict):
    """
    Benchmark mounted wB97M-V fixture values against literature targets.
    """
    benchmark = PHASE33_BENCHMARKS[family]
    literature = benchmark["literature"]
    barrier = benchmark["wb97mv_kcal_mol"]

    assert literature == literature_barriers_dict[family]
    assert literature["low"] < barrier < literature["high"]
    assert abs(barrier - expected_best) <= 1.5

@pytest.mark.slow
@pytest.mark.optional_dft_authority_lane
@pytest.mark.parametrize("family, _", REACTION_FAMILIES)
@pytest.mark.skipif(not HAS_PHASE35_DATASET, reason=PHASE35_SKIP_REASON)
def test_barrier_benchmark_revdsd_pbep86(family, _, literature_barriers_dict):
    """
    Validate mounted wB97M-V values against double-hybrid comparison fixtures.
    """
    phase33 = PHASE33_BENCHMARKS[family]
    phase35 = PHASE35_BENCHMARKS[family]
    literature = literature_barriers_dict[family]

    assert phase35["wb97mv_kcal_mol"] == phase33["wb97mv_kcal_mol"]
    assert abs(phase35["wb97mv_kcal_mol"] - phase35["revdsd_pbep86_d4_kcal_mol"]) < 2.0
    assert literature["low"] < phase35["revdsd_pbep86_d4_kcal_mol"] < literature["high"]

@pytest.mark.slow
@pytest.mark.optional_dft_authority_lane
@pytest.mark.skipif(not HAS_PHASE33_DATASET, reason=PHASE33_SKIP_REASON)
def test_barrier_correlation_dft_vs_xtb(literature_barriers_dict, xtb_barriers_dict):
    """
    Verify that rank-order correlation between DFT and xTB is preserved.
    """
    families = [family for family, _ in REACTION_FAMILIES]
    dft_values = [PHASE33_BENCHMARKS[family]["wb97mv_kcal_mol"] for family in families]
    xtb_values = [xtb_barriers_dict[family] for family in families]

    for family in families:
        assert PHASE33_BENCHMARKS[family]["xtb_reference_kcal_mol"] == xtb_barriers_dict[family]
        assert PHASE33_BENCHMARKS[family]["literature"] == literature_barriers_dict[family]

    correlation = spearmanr(dft_values, xtb_values).statistic
    assert correlation is not None
    assert correlation > 0.95
