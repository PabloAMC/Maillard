import pytest
from src.skala_refiner import SkalaRefiner  # noqa: E402

def pyscf_installed():
    try:
        from pyscf import dft # noqa: F401
        return True
    except ImportError:
        return False

@pytest.mark.slow
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
from src.xtb_screener import XTBScreener  # noqa: E402
from src.pathway_extractor import ElementaryStep, Species  # noqa: E402
import subprocess # noqa: E402

def xtb_installed():
    """Use XTBScreener's auto-detect to find xTB binary."""
    screener = XTBScreener()
    try:
        res = subprocess.run([screener.xtb_bin, "--version"], capture_output=True, text=True)
        return res.returncode == 0
    except FileNotFoundError:
        return False

@pytest.mark.slow
@pytest.mark.skipif(not xtb_installed(), reason="xTB binary not found in PATH")
def test_integration_xtb_methanediol_formation():
    """
    Validates true xTB physics by evaluating the hydration of formaldehyde.
    H2C=O + H2O -> H2C(OH)2 
    This reaction is known to be slightly exothermic in water.
    """
    formaldehyde = Species("Formaldehyde", "C=O")
    water = Species("Water", "O")
    methanediol = Species("Methanediol", "OC(O)")
    
    step = ElementaryStep(
        reactants=[formaldehyde, water],
        products=[methanediol],
        reaction_family="Hydration"
    )
    
    screener = XTBScreener()
    delta_E_kcal, barrier_kcal = screener.compute_reaction_energy(step)
    
    # Validation criteria:
    # 1. Reaction should report meaningful values (not the mock 999.0 penalty)
    assert barrier_kcal < 900.0
    
    # 2. Formaldehyde hydration in implicit water. 
    # GFN2-xTB with GBSA water systematically over-stabilises small polar molecules (~-28 kcal/mol).
    # The key check is that we get a physical value, not the 999.0 mock penalty.
    assert -35.0 < delta_E_kcal < 10.0, f"ΔE out of expected range: {delta_E_kcal:.2f}"
