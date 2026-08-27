import pytest

from src.barrier_constants import (  # noqa: E402
    DEFAULT_BARRIER,
    DEFAULT_REFERENCE_PREEXPONENTIAL,
    arrhenius_rate_constant,
    effective_barrier_from_rate_constant,
    get_arrhenius_params,
    get_barrier,
    get_reference_pre_exponential,
)

def test_get_barrier_exact_match():
    # Test valid exact matches and normalized variants
    assert get_barrier("schiff_condensation")[0] == 15.0
    assert get_barrier("Schiff Condensation")[0] == 15.0 # Normalizes to schiff_condensation
    assert get_barrier("1,2-enolisation")[0] == 28.0
    assert get_barrier("1_2-enolisation")[0] == 28.0
    assert get_barrier("cysteine_thermolysis")[0] == 30.0
    
    # Test invalid or unknown matches fall back
    assert get_barrier("unknown_reaction")[0] == DEFAULT_BARRIER
    assert get_barrier("")[0] == DEFAULT_BARRIER

def test_get_arrhenius_params_valid():
    # Test that valid families from the YAML load correctly
    # Note: tests depend on the actual contents of arrhenius_params.yml
    
    # Schiff base
    params = get_arrhenius_params("schiff_condensation")
    assert params is not None
    A, Ea_kcal, _, _ = params
    assert A == 1.5e11
    # Audit 2026-08-26: corrected to the thesis value (Table 5.2, Glu+Gly->DFG).
    assert abs(Ea_kcal - 23.18) < 0.1 # 97.0 kJ/mol / 4.184 = 23.18 kcal/mol
    
    # Amadori
    params = get_arrhenius_params("amadori_rearrangement")
    # A_value is 1.0e11 in YAML
    assert params is not None
    A, Ea_kcal, quality, _ = params
    assert A == 1.0e11
    assert quality == "literature_estimated"
    
    # Dehydration (fructose)
    params = get_arrhenius_params("dehydration")
    assert params is not None
    A, Ea_kcal, _, _ = params
    assert A == 3.09e11
    assert abs(Ea_kcal - 29.6) < 0.1 # 123.88 kJ/mol / 4.184

def test_get_arrhenius_params_invalid():
    # Unknown families should return None
    assert get_arrhenius_params("magic_reaction") is None
    assert get_arrhenius_params("") is None

def test_get_arrhenius_params_normalization():
    # Should normalize "schiff" to "schiff_condensation"
    assert get_arrhenius_params("schiff") == get_arrhenius_params("schiff_condensation")


def test_enolisation_uses_enolisation_arrhenius_contract_not_fructose_dehydration():
    params = get_arrhenius_params("1,2-enolisation")
    assert params is not None
    A, ea_kcal, quality, _ = params

    assert A == 1e12
    assert quality == "literature_estimated"
    assert abs(ea_kcal - (120.0 / 4.184)) < 0.1


def test_fast_enolisation_barrier_stays_close_to_arrhenius_reference():
    fast_barrier, _ = get_barrier("1,2-enolisation")
    params = get_arrhenius_params("1,2-enolisation")
    assert params is not None
    _, arrhenius_barrier, _, _ = params

    assert abs(fast_barrier - arrhenius_barrier) < 1.0


def test_reference_pre_exponential_uses_literature_then_fallback():
    assert get_reference_pre_exponential("schiff_condensation") == 1.5e11
    assert get_reference_pre_exponential("magic_reaction") == DEFAULT_REFERENCE_PREEXPONENTIAL


def test_arrhenius_round_trip_preserves_family_barrier():
    barrier = 13.62
    temperature_kelvin = 423.15
    family = "schiff_condensation"

    rate_constant = arrhenius_rate_constant(barrier, temperature_kelvin, family=family)
    recovered = effective_barrier_from_rate_constant(rate_constant, temperature_kelvin, family=family)

    assert recovered == pytest.approx(barrier, rel=1e-9)
