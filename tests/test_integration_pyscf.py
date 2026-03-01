import pytest
from src.skala_refiner import SkalaRefiner

def pyscf_installed():
    try:
        import pyscf
        return True
    except ImportError:
        return False

@pytest.mark.slow
@pytest.mark.skipif(not pyscf_installed(), reason="PySCF not installed")
def test_integration_pyscf_water():
    """
    Validates PySCF / Skala backend plumbing with a live single point calculation.
    """
    xyz = "3\n\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n"
    
    # Use lowest possible basis set and LDA just to test plumbing without Skala dependency
    refiner = SkalaRefiner(basis='sto-3g', solvent_name=None, use_skala=False)
    refiner.fallback_xc = 'lda' 
    
    res = refiner.single_point(xyz)
    
    assert res.converged
    # Approximate Hartree-Fock/LDA energy of water in STO-3G is ~ -75.0 Hartree
    assert -80.0 < res.energy_hartree < -70.0
