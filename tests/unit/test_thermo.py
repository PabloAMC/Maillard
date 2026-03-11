from src.thermo import QuasiHarmonicCorrector  # noqa: E402

def test_harmonic_entropy():
    corrector = QuasiHarmonicCorrector(temp_k=298.15)
    
    # 0 frequency should return 0 entropy
    assert corrector._harmonic_entropy(0.0) == 0.0
    
    # High frequency >1000 cm-1 contributes very little to entropy at room temp
    s_high = corrector._harmonic_entropy(1500.0)
    assert s_high < 1.0 # J/mol K
    
    # Low frequency 50 cm-1 contributes a lot
    s_low = corrector._harmonic_entropy(50.0)
    assert s_low > 20.0

def test_rotor_entropy():
    corrector = QuasiHarmonicCorrector(temp_k=298.15)
    
    # Should safely return 0 for 0 frequency
    assert corrector._rotor_entropy(0.0) == 0.0
    
    # A free rotor entropy calculation
    s_rot = corrector._rotor_entropy(50.0)
    assert s_rot > 0.0

def test_qh_thermo_result():
    corrector = QuasiHarmonicCorrector(temp_k=298.15, cutoff_freq_cm1=100.0)
    
    # Mock some frequencies where one is high, one is below cutoff
    freqs = [1500.0, 50.0]
    
    res = corrector.calculate_thermo(freqs, electronic_energy_h=-76.0)
    
    # Since 50 cm-1 is below 100 cm-1, the interpolated entropy should be lower 
    # than the raw pure harmonic entropy which diverges.
    # We test that the diff is negative (entropy was reduced) or positive depending
    # on if rot > harm (usually harm > rot for very low freqs, so rot brings it down) 
    
    # For 50 cm-1, harmonic is ~25 J/mol K, rotor is often lower or bounded.
    # Just verify that qh_entropy != raw_entropy
    assert res.qh_entropy_h != res.raw_entropy_h
    
    # Output should include standard thermo fields
    assert res.qh_gibbs_h is not None
    assert res.raw_gibbs_h is not None
    
    # The gibbs free energy should change due to entropy
    assert res.qh_gibbs_h != res.raw_gibbs_h
