import pytest
from src.skala_refiner import SkalaRefiner

# Methane XYZ (simplest test case)
METHANE_XYZ = """5

C 0.000000 0.000000 0.000000
H 0.627600 0.627600 0.627600
H -0.627600 -0.627600 0.627600
H -0.627600 0.627600 -0.627600
H 0.627600 -0.627600 -0.627600
"""

@pytest.mark.slow
def test_skala_single_point():
    """Run a fast, extremely low-basis single point calculation to verify PySCF plumbing."""
    # Use STO-3G for speed in testing; turn Skala off if not installed (fallback to basic DFT/HF)
    refiner = SkalaRefiner(basis='sto-3g', solvent_name=None, use_skala=False)
    
    # We use LDA (slater) for the fastest possible mock run
    refiner.fallback_xc = 'lda'
    
    res = refiner.single_point(METHANE_XYZ)
    
    # Methane STO-3G LDA energy is roughly -39.7 to -39.9 Hartree
    assert res.converged is True
    assert -41.0 < res.energy_hartree < -38.0
    
@pytest.mark.slow
def test_skala_optimization():
    """Test geometric optimization interface with a fast basis set."""
    refiner = SkalaRefiner(basis='sto-3g', solvent_name=None, use_skala=False)
    refiner.fallback_xc = 'lda'
    
    res = refiner.optimize_geometry(METHANE_XYZ, max_steps=5) # cap steps for test speed
    
    assert res.converged is True
    assert res.energy_hartree < 0.0
    assert "C" in res.optimized_xyz
    assert "H" in res.optimized_xyz
