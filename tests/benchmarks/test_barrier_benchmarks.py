import pytest

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

@pytest.mark.slow
@pytest.mark.parametrize("family, expected_best", REACTION_FAMILIES)
def test_barrier_benchmark_wb97mv(family, expected_best, literature_barriers_dict):
    """
    Benchmark wB97M-V barriers against literature targets.
    Currently acts as a placeholder for Phase 3.3 batch execution.
    """
    pytest.skip(f"Benchmark run for {family} is pending Phase 3.3 data.")
    # barrier = compute_barrier(..., functional='wB97M-V')
    # lit = literature_barriers_dict[family]
    # assert lit['low'] < barrier < lit['high']

@pytest.mark.slow
@pytest.mark.parametrize("family, _", REACTION_FAMILIES)
def test_barrier_benchmark_revdsd_pbep86(family, _, literature_barriers_dict):
    """
    Validate wB97M-V against Double-Hybrid revDSD-PBEP86-D4 (Phase 3.5).
    Ensures single-hybrid barriers are not drifting significantly.
    """
    pytest.skip(f"Double-hybrid verification for {family} is pending Phase 3.5.")
    # barrier_wb97mv = compute_barrier(..., functional='wB97M-V')
    # barrier_dsd = compute_barrier(..., functional='revDSD-PBEP86-D4')
    # assert abs(barrier_wb97mv - barrier_dsd) < 2.0

@pytest.mark.slow
def test_barrier_correlation_dft_vs_xtb(literature_barriers_dict, xtb_barriers_dict):
    """
    Verify that rank-order correlation between DFT and xTB is preserved.
    """
    pytest.skip("Full network correlation check pending Phase 3.3 completion.")
