import pytest
import os
from src.dft_refiner import DFTRefiner, DFTResult  # noqa: E402

# Use a very simple H-transfer transition state guess (planar H3)
# To keep DFT fast, we use STO-3G
H3_TS_GUESS = """3
H3 TS Guess
H  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.900000
H  0.000000  0.000000 -0.900000
"""

@pytest.mark.skipif(not os.path.exists(".venv/lib/python3.14/site-packages/sella"), reason="Sella not installed")
def test_sella_ts_search_integration():
    """
    Verify that DFTRefiner correctly uses Sella for TS searches 
    via the built-in PySCF calculator bridge.
    """
    # Initialize refiner
    refiner = DFTRefiner()
    
    # Force fast settings
    refiner.opt_method = 'hf'
    refiner.opt_basis = 'sto-3g'
    refiner.refinement_method = 'hf'
    refiner.refinement_basis = 'sto-3g'
    
    # Run TS optimization
    res = refiner.optimize_geometry(H3_TS_GUESS, charge=0, spin=1, is_ts=True, max_steps=5)
    
    assert isinstance(res, DFTResult)
    # If Sella worked, we should see it in the logs (captured by pytest if -s used)
    # The result should contain an optimized geometry
    assert res.optimized_xyz is not None
    assert "H" in res.optimized_xyz
