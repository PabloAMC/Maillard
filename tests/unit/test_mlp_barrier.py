"""
tests/unit/test_mlp_barrier.py

Verifies MACE-OFF activation energy estimation.
"""

import pytest
from src.mlp_barrier import MLPBarrier  # noqa: E402
from src.results_db import ResultsDB  # noqa: E402
from src.dft_refiner import DFTRefiner  # noqa: E402

# Simple balanced stoichiometry: H2 + H -> H + H2 (as a dummy check)
H_MOL = """2
H2 molecule
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""

H_MOL_STRETCHED = """2
H2 stretched
H 0.0 0.0 0.0
H 0.0 0.0 1.5
"""

# Unbalanced stoichiometry
H3_MOL = """3
H3 molecule
H 0.0 0.0 0.0
H 0.0 0.0 1.0
H 1.0 0.0 0.0
"""

def test_mlp_energy_calculation():
    """Verify that MACE can compute a sanity energy for H2."""
    mlp = MLPBarrier(model="medium", device="cpu")
    energy = mlp.get_energy(H_MOL)
    assert isinstance(energy, float)
    assert energy < 0  # Bound system should have negative energy

def test_mlp_barrier_stoichiometry_check():
    """Verify that mismatched atom counts return None instead of huge numbers."""
    mlp = MLPBarrier(model="medium", device="cpu")
    
    # Matching stoichiometry
    barrier = mlp.estimate_barrier(H_MOL, H_MOL_STRETCHED)
    assert barrier is not None
    assert barrier > 0
    
    # Mismatched stoichiometry
    barrier_mismatch = mlp.estimate_barrier(H_MOL, H3_MOL)
    assert barrier_mismatch is None

def test_dft_refiner_mlp_integration(tmp_path):
    """Verify that DFTRefiner correctly uses MLP and logs to ResultsDB."""
    db_file = tmp_path / "test_maillard.db"
    refiner = DFTRefiner(db_path=str(db_file))
    
    reaction_meta = {
        "reactants": ["[H][H]"],
        "products": ["[H][H]"],
        "family": "H2_stretch"
    }
    
    barrier_kcal = refiner.calculate_mlp_barrier(H_MOL, H_MOL_STRETCHED, reaction_meta=reaction_meta)
    
    assert barrier_kcal is not None
    assert barrier_kcal > 0
    
    # Check if it was saved to DB
    db = ResultsDB(db_path=str(db_file))
    res = db.find_barrier(["[H][H]"], ["[H][H]"])
    assert res is not None
    assert res["method"] == "mace-off"
    assert abs(res["delta_g_kcal"] - barrier_kcal) < 1e-5

if __name__ == "__main__":
    pytest.main([__file__])
