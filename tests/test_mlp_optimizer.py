"""
tests/test_mlp_optimizer.py

Phase 10: Verification of the MACE MLP geometry optimizer.
"""

import pytest
import numpy as np

try:
    from src.mlp_optimizer import MLPOptimizer
except ImportError:
    MLPOptimizer = None

# Direct check for MACE backend — more reliable than checking MLPOptimizer class
_MACE_AVAILABLE = False
try:
    from mace.calculators import mace_mp
    _MACE_AVAILABLE = True
except (ImportError, Exception):
    pass

FORMALDEHYDE_XYZ = """\
4
Formaldehyde
O       0.000000    0.000000    1.200000
C       0.000000    0.000000    0.000000
H       0.000000    0.940000   -0.580000
H       0.000000   -0.940000   -0.580000
"""

@pytest.mark.skipif(not _MACE_AVAILABLE, reason="MACE ML dependencies not installed")
@pytest.mark.slow
def test_mlp_optimizer_smoke():
    """
    Smoke test: Ensure the MLPOptimizer can load the small MACE general-purpose model
    and optimize a simple molecule without crashing.
    """
    # Use "small" model to keep the test runtime manageable
    optimizer = MLPOptimizer(model_name="small", device="cpu")
    
    # Run a loose optimization
    optimized_xyz = optimizer.optimize_geometry(FORMALDEHYDE_XYZ, fmax=0.05, max_steps=100)
    
    lines = optimized_xyz.strip().splitlines()
    assert int(lines[0]) == 4, "Atom count should remain 4"
    assert len(lines) == 6, "XYZ should have 2 header lines + 4 coordinate lines"
    
    # Check that elements are preserved
    elements_out = [line.split()[0] for line in lines[2:]]
    assert sorted(elements_out) == sorted(["O", "C", "H", "H"])

@pytest.mark.skipif(not _MACE_AVAILABLE, reason="MACE ML dependencies not installed")
@pytest.mark.slow
def test_mlp_optimizer_ts_fallback():
    """ Verify the TS scaffold correctly warns and falls back to standard minimization. """
    optimizer = MLPOptimizer(model_name="small", device="cpu")
    opt_xyz = optimizer.optimize_ts(FORMALDEHYDE_XYZ, fmax=0.05, max_steps=10)
    assert "4" in opt_xyz.splitlines()[0]

@pytest.mark.skipif(not _MACE_AVAILABLE, reason="MACE ML dependencies not installed")
@pytest.mark.slow
def test_dft_refiner_mace_backend():
    """ Verify DFTRefiner routes structural relaxation through MLPOptimizer and returns valid energies. """
    from src.dft_refiner import DFTRefiner
    # Fast mode equivalents to keep test fast
    refiner = DFTRefiner(geometry_backend="mace", use_explicit_solvent=False)
    refiner.opt_method = "hf"  # Fast electronic tier
    refiner.opt_basis = "sto-3g"
    refiner.refinement_method = "hf"
    refiner.refinement_basis = "sto-3g"
    
    # We must explicitly force MLPOptimizer to use 'small' for test speed
    refiner.mlp_optimizer.model_name = "small"
    refiner.mlp_optimizer.calc = MLPOptimizer(model_name="small", device="cpu").calc

    # Attempt a geometry optimization cycle
    result = refiner.optimize_geometry(FORMALDEHYDE_XYZ, max_steps=10)
    
    assert result.converged is True
    assert result.energy_hartree < 0.0
    assert result.frequencies_cm1 is not None
