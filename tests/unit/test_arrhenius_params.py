from src.barrier_constants import get_barrier, get_arrhenius_params, DEFAULT_BARRIER  # noqa: E402

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
    assert abs(Ea_kcal - 13.62) < 0.1 # 57.0 kJ/mol / 4.184 = 13.62 kcal/mol
    
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
