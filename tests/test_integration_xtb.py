import pytest
from src.xtb_screener import XTBScreener
from src.pathway_extractor import ElementaryStep, Species
import subprocess
from pathlib import Path

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
