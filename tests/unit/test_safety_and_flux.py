
import math
import pytest
from src.inverse_design import InverseDesigner  # noqa: E402
from src.smirks_engine import ReactionConditions  # noqa: E402

def test_safety_scoring_acrylamide():
    """Verify that asparagine-rich formulations trigger safety penalties."""
    designer = InverseDesigner(target_tag="roasted", minimize_tag="beany")
    
    # Mock some grid entries if they don't exist, but we can use the real ones if they serve
    # Or just run evaluate_all and check results
    conds = ReactionConditions(pH=7.0, temperature_celsius=160)
    results = designer.evaluate_all(conds)
    
    # Check if any formulation has Asparagine in the grid and verify penalty
    # We'll create a manual test case to be sure
    class MockFormulation:
        def __init__(self, name, sugars, amino_acids):
            self.name = name
            self.sugars = sugars
            self.amino_acids = amino_acids
            self.ph = 7.0
            self.temp = 160
            
        def get(self, key, default=None):
            return getattr(self, key, default)

    # Injecting a custom formulation into the grid for testing
    asparagine_form = {
        "name": "AcrylamideRisk",
        "sugars": ["D-Glucose"],
        "amino_acids": ["L-Asparagine"],
        "ph": 7.0,
        "temp": 160
    }
    glycine_form = {
        "name": "SafetyBaseline",
        "sugars": ["D-Glucose"],
        "amino_acids": ["Glycine"],
        "ph": 7.0,
        "temp": 160
    }
    
    designer.grid = [asparagine_form, glycine_form]
    results = designer.evaluate_all(conds)
    
    # Find results
    acry_res = next(r for r in results if r.name == "AcrylamideRisk")
    safe_res = next(r for r in results if r.name == "SafetyBaseline")
    
    # AcrylamideRisk should have a safety score and Acrylamide flagged
    assert acry_res.safety_score > 0
    assert "Acrylamide" in acry_res.flagged_toxics
    
    # SafetyBaseline might have HMF, but should NOT have Acrylamide
    assert "Acrylamide" not in safe_res.flagged_toxics
    # AcrylamideRisk should be way more dangerous than baseline due to acrylamide priority
    assert acry_res.safety_score > safe_res.safety_score

def test_concentration_aware_ranking():
    """Verify that doubling precursor concentration doubles the flux/score."""
    designer = InverseDesigner(target_tag="meaty")
    
    # Formulation 1: 0.1 M Cysteine
    cys_low = {
        "name": "LowCys",
        "sugars": ["D-Ribose"],
        "amino_acids": ["L-Cysteine"],
        "molar_ratios": {"D-Ribose": 0.1, "L-Cysteine": 0.1},
        "ph": 6.0,
        "temp": 140
    }
    # Formulation 2: 1.0 M Cysteine
    cys_high = {
        "name": "HighCys",
        "sugars": ["D-Ribose"],
        "amino_acids": ["L-Cysteine"],
        "molar_ratios": {"D-Ribose": 1.0, "L-Cysteine": 1.0},
        "ph": 6.0,
        "temp": 140
    }
    
    designer.grid = [cys_low, cys_high]
    conds = ReactionConditions(pH=6.0, temperature_celsius=140)
    results = designer.evaluate_all(conds)
    
    res_low = next(r for r in results if r.name == "LowCys")
    res_high = next(r for r in results if r.name == "HighCys")
    
    # target_score now uses Sensory radar intensity (Stevens' Law)
    # Intensity ~ Conc^0.5
    # High (1.0M) is 10x conc of Low (0.1M)
    # Target score should increase by sqrt(10) ~ 3.16x per Stevens' Law
    assert res_high.target_score == pytest.approx(math.sqrt(10.0) * res_low.target_score, rel=0.1)

if __name__ == "__main__":
    pytest.main([__file__])
