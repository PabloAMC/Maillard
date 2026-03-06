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

@pytest.mark.slow
def test_skala_ts_optimization_mock(monkeypatch):
    """Verify TS optimization interface."""
    # We mock geometric_solver to avoid actual expensive TS search in unit test
    class MockMol:
        def tostring(self, format): return METHANE_XYZ
    
    def mock_optimize(mf, **kwargs):
        assert kwargs.get('custom_engine_params', {}).get('transition') is True
        return MockMol()
        
    import pyscf.geomopt.geometric_solver
    monkeypatch.setattr(pyscf.geomopt.geometric_solver, "optimize", mock_optimize)
    
    refiner = SkalaRefiner(basis='sto-3g', solvent_name=None, use_skala=False)
    refiner.fallback_xc = 'lda'
    
    # This will call our mock_optimize
    res = refiner.optimize_ts(METHANE_XYZ)
    assert "C" in res.optimized_xyz
    assert res.converged is True

def test_skala_import_fallback(capsys):
    """Test that it falls back to B3LYP if Skala is missing."""
    refiner = SkalaRefiner(use_skala=True)
    # Force 'import skala' to fail by monkeypatching sys.modules if needed, 
    # but the code already has a try-except. 
    # Let's just verify what happens when Skala is NOT installed (default state usually).
    
    # We can check the print output
    from pyscf import gto
    mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
    mf = refiner._build_mf(mol)
    
    if mf.xc != 'SKALA':
        captured = capsys.readouterr()
        # The print happens in _build_mf if import fails
        # Actually mf.xc will be 'b3lyp' (fallback_xc)
        assert mf.xc == 'b3lyp'
