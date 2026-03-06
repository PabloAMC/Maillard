import pytest
import os
import yaml
import tempfile
from pathlib import Path

from src.inverse_design import InverseDesigner
from src.conditions import ReactionConditions

@pytest.fixture
def mock_grid_path():
    """Create a temporary formulation_grid.yml with known concentration ratios to test Boltzmann scoring."""
    grid_data = {
        "formulations": [
            {
                "name": "Low Cysteine Base",
                "sugars": ["ribose"],
                "amino_acids": ["cysteine"],
                "molar_ratios": {
                    "cysteine": 0.5,
                    "ribose": 1.0
                }
            },
            {
                "name": "High Cysteine Base",
                "sugars": ["ribose"],
                "amino_acids": ["cysteine"],
                "molar_ratios": {
                    "cysteine": 2.0,
                    "ribose": 1.0
                }
            }
        ]
    }
    
    fd, path = tempfile.mkstemp(suffix=".yml")
    with os.fdopen(fd, 'w') as f:
        yaml.dump(grid_data, f)
        
    yield Path(path)
    os.remove(path)

def test_concentration_boltzmann_scoring(mock_grid_path, monkeypatch):
    """
    Formally verifies Phase 8.D Sub-task B:
    Two identical formulations (except for concentration ratio) must receive 
    different scores, with the higher-concentration precursor yielding a higher score.
    """
    import src.inverse_design as inv
    monkeypatch.setattr(inv, "GRID_FILE", mock_grid_path)
    
    designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
    
    cond = ReactionConditions(pH=5.5, temperature_celsius=150.0)
    results = designer.evaluate_all(cond)
    
    assert len(results) == 2
    
    # Extract the scores for both runs
    high_cys = next(r for r in results if r.name == "High Cysteine Base")
    low_cys = next(r for r in results if r.name == "Low Cysteine Base")
    
    # High cysteine formulation should have a strictly greater target score 
    # than the low cysteine formulation because the limiting reagent is 4x higher.
    assert high_cys.target_score > low_cys.target_score, (
        f"High Cys score ({high_cys.target_score}) must be > Low Cys score ({low_cys.target_score})"
    )
