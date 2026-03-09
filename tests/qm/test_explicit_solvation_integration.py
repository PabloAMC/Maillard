import pytest
import os
from src.dft_refiner import DFTRefiner

# Simple formaldehyde model
FORMALDEHYDE_XYZ = """4
Formaldehyde
O       0.000000    0.000000    1.200000
C       0.000000    0.000000   -0.000000
H       0.000000    0.940000   -0.580000
H       0.000000   -0.940000   -0.580000
"""

def test_dft_refiner_explicit_solvation():
    """
    Smoke test: Run DFTRefiner with 1 explicit water molecule.
    Uses 'Fast Mode' (HF/STO-3G) to keep the test quick.
    """
    # Initialize refiner with fast-mode equivalents
    refiner = DFTRefiner(use_explicit_solvent=True, n_water=1)
    refiner.opt_method = 'hf'
    refiner.opt_basis = 'sto-3g'
    refiner.refinement_method = 'hf'
    refiner.refinement_basis = 'sto-3g'
    
    # We only run a few steps to verify the pipeline doesn't crash
    # and that the solvent cluster is handled.
    print("\nRunning fast-mode explicit solvation smoke test...")
    
    # We'll use is_ts=True to trigger the constraint logic
    result = refiner.optimize_geometry(FORMALDEHYDE_XYZ, is_ts=True, max_steps=2)
    
    # Assertions
    assert result.converged is False or result.converged is True # We don't care about convergence in 2 steps
    assert "Solvated Cluster Guess" in result.optimized_xyz or len(result.optimized_xyz.split('\n')) > 6
    
    # Verify atom count (Formaldehyde 4 + Water 3 = 7 atoms)
    lines = result.optimized_xyz.strip().split('\n')
    n_atoms = int(lines[0])
    assert n_atoms == 7, f"Expected 7 atoms (4+3), got {n_atoms}"
    
    print("Smoke test passed: Explicit solvation cluster generated and passed to geomeTRIC.")

if __name__ == "__main__":
    test_dft_refiner_explicit_solvation()
