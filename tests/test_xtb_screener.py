import pytest
from src.xtb_screener import XTBScreener
from src.pathway_extractor import ElementaryStep, Species
import tempfile
from pathlib import Path

def test_generate_3d_xyz():
    screener = XTBScreener()
    xyz = screener.generate_3d_xyz("O") # Water
    
    assert xyz is not None
    assert "O" in xyz
    assert "H" in xyz
    lines = xyz.strip().split("\n")
    assert lines[0] == "3" # 3 atoms in water

def test_xtb_screener_execution(monkeypatch):
    """Mock the subprocess call to xTB to test parsing logic without binary dependency."""
    
    class MockProcess:
        returncode = 0
        stdout = "          | TOTAL ENERGY               -13.064547948281 Eh   |\n"
        stderr = ""

    def mock_run(*args, **kwargs):
        # Create a fake xtbopt.xyz
        run_dir = kwargs.get("cwd")
        if run_dir:
            (run_dir / "xtbopt.xyz").write_text("3\n\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n")
        return MockProcess()

    import subprocess
    monkeypatch.setattr(subprocess, "run", mock_run)
    
    screener = XTBScreener()
    res = screener.optimize_species("O")
    
    assert res.energy_hartree == -13.064547948281
    assert "O 0.0 0.0 0.0" in res.optimized_xyz

def test_run_xtb_path(monkeypatch):
    """Test NEB path parsing logic from mocked file output."""
    class MockProcess:
        returncode = 0
        stdout = "NEB run finished"
        stderr = ""

    def mock_run(*args, **kwargs):
        run_dir = kwargs.get("cwd")
        if run_dir:
            # Fake xtbpath.xyz with 3 images, energy in comments
            # image 1: E = -10.0
            # image 2: E = -9.0 (TS)
            # image 3: E = -11.0
            content = (
                "3\n energy: -10.0\nA 0 0 0\nB 1 1 1\nC 2 2 2\n"
                "3\n energy: -9.0\nA 0 0 0\nB 1 1 1\nC 2 2 2\n"
                "3\n energy: -11.0\nA 0 0 0\nB 1 1 1\nC 2 2 2\n"
            )
            (run_dir / "xtbpath.xyz").write_text(content)
        return MockProcess()

    import subprocess
    monkeypatch.setattr(subprocess, "run", mock_run)
    
    screener = XTBScreener()
    with tempfile.TemporaryDirectory() as td:
        barrier_kcal = screener._run_xtb_path("fake_r", "fake_p", Path(td))
        
    # max_E = -9.0, start_E = -10.0 -> barrier = 1.0 Hartree
    # 1.0 * 627.509 = 627.509 kcal/mol
    assert abs(barrier_kcal - 627.509) < 0.1

def test_compute_reaction_energy(monkeypatch):
    """Test thermodynamics pipeline with mocked xTB energies."""
    
    fake_energies = {
        "O": -13.0,
        "C(=O)O": -20.0,
        "C": -4.0 # Methane
    }
    
    def mock_optimize_species(self, smiles):
        from src.xtb_screener import XTBResult
        E = fake_energies.get(smiles, -10.0) 
        return XTBResult(energy_hartree=E, optimized_xyz="mock")

    def mock_run_xtb_path(self, r_xyz, p_xyz, run_dir):
        return 12.5  # Mock NEB barrier directly

    monkeypatch.setattr(XTBScreener, "optimize_species", mock_optimize_species)
    monkeypatch.setattr(XTBScreener, "_run_xtb_path", mock_run_xtb_path)
    
    step = ElementaryStep(
        reactants=[Species("Water", "O"), Species("Methane", "C")], 
        products=[Species("Acid", "C(=O)O")]                        
    )
    
    screener = XTBScreener()
    delta_E, barrier = screener.compute_reaction_energy(step)
    
    # E = -20 - (-13 - 4) = -3 Hartree = -1882.527 kcal/mol
    assert abs(delta_E - (-1882.527)) < 0.1
    assert barrier == 12.5
