import pytest
import os
from src.dft_refiner import DFTRefiner, DFTResult  # noqa: E402

# Simple water molecule for fast testing
WATER_XYZ = """3
water
O 0.000000 0.000000 0.000000
H 0.000000 0.000000 0.958400
H 0.000000 0.924000 -0.252900
"""

@pytest.mark.skipif(not os.path.exists("conda_env/bin/crest"), reason="CREST binary not found")
def test_dft_refiner_explicit_solvation_integration():
    """
    Verify that DFTRefiner correctly triggers SolvationEngine.
    We use --fast mode (HF/def2-svp) to keep it quick.
    """
    refiner = DFTRefiner(use_explicit_solvent=True, n_water=1)
    
    # We override the methods to be as fast as possible for a smoke test
    refiner.opt_method = 'hf'
    refiner.opt_basis = 'sto-3g'
    refiner.refinement_method = 'hf'
    refiner.refinement_basis = 'sto-3g'
    
    # Run optimization
    # Note: we use is_ts=False for simplicity in this integration check
    res = refiner.optimize_geometry(WATER_XYZ, is_ts=False, max_steps=2)
    
    assert isinstance(res, DFTResult)
    # The optimized XYZ should have 6 atoms (3 solute + 3 water)
    lines = res.optimized_xyz.strip().split('\n')
    assert int(lines[0]) == 6
    assert res.converged is True or res.converged is False # Just check it returns a bool
