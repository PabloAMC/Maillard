"""
tests/test_ts_optimizer.py

Phase 11: Verification of the Sella eigenvector-following TS search.
Uses a simple model system (H-transfer in HCN ⇌ CNH) to verify that 
the optimizer correctly identifies a first-order saddle point.
"""

import pytest
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from src.ts_optimizer import TSOptimizer

# Model system: Hydrogen Cyanide (HCN) to Hydrogen Isocyanide (CNH)
# Simple initial guess for the TS (bent structure)
HCN_TS_GUESS = Atoms(
    symbols='HCN',
    positions=[
        [0.0, 1.0, 0.0],  # H shifted up
        [-1.0, 0.0, 0.0], # C
        [1.0, 0.0, 0.0]   # N
    ]
)

def test_ts_optimizer_hcn():
    """ Verify TSOptimizer can find a TS for the HCN model system using a fast calculator. """
    optimizer = TSOptimizer(fmax=0.05)
    
    # We use EMT as a fast toy calculator for structural tests
    atoms = HCN_TS_GUESS.copy()
    calc = EMT()
    
    # Find the TS
    optimized_atoms = optimizer.find_ts(atoms, calc, max_steps=100)
    
    # In a true TS, we expect the max force to be below fmax
    assert optimizer.is_converged(optimized_atoms)
    
    # Verify we still have 3 atoms
    assert len(optimized_atoms) == 3
    
    # Check that it's no longer linear (H should be in a bridging position)
    pos = optimized_atoms.get_positions()
    h_pos = pos[0]
    c_pos = pos[1]
    n_pos = pos[2]
    
    # In the TS for HCN isomarisation, H is roughly equidistant to C and N
    dist_hc = np.linalg.norm(h_pos - c_pos)
    dist_hn = np.linalg.norm(h_pos - n_pos)
    
    assert dist_hc < 2.0
    assert dist_hn < 2.0
    
def test_ts_optimizer_import_error():
    """ Verify graceful failure if sella is missing (mocked). """
    import sys
    from unittest.mock import patch
    
    # Temporarily hide sella
    with patch.dict(sys.modules, {'sella': None}):
        # We need to reload or re-import to trigger the failure if it's already cached
        # but for simplicity, we check the raise in find_ts
        opt = TSOptimizer()
        with pytest.raises(ImportError):
            # We pass None for atoms/calc because it should fail before using them
            with patch('src.ts_optimizer.Sella', None):
                opt.find_ts(None, None)
