"""
tests/unit/test_sensory_full.py

Verifies the unified SensoryDatabase and psychophysical mixing.
"""

import pytest
from src.sensory import SensoryDatabase, SensoryPredictor  # noqa: E402
from src.headspace import HeadspaceModel  # noqa: E402

def test_sensory_database_loading():
    """Verify that we load compounds from multiple YAML files."""
    db = SensoryDatabase()
    
    # Check FFT (from desirable_targets.yml)
    fft = db.find_entry("2-Furfurylthiol (FFT)")
    assert fft is not None
    assert fft["threshold_ppm"] == 0.00001 # 0.01 ppb / 1000
    
    # Check Hexanal (from off_flavour_targets.yml)
    hexanal = db.find_entry("Hexanal")
    assert hexanal is not None
    assert hexanal["threshold_ppm"] == 0.0045 # 4.5 ppb / 1000
    
    # Check CML (from toxic_markers.yml - might not have threshold)
    cml = db.find_entry("Nε-(Carboxymethyl)lysine (CML)")
    assert cml is not None

def test_psychophysical_scaling():
    """Verify Stevens' Law scaling: Intensity = (C/Threshold)^0.5"""
    predictor = SensoryPredictor()
    
    # threshold = 1.0 ppb, conc = 100.0 ppb -> OAV = 100 -> Intensity = 10
    mock_entry = {
        "name": "TestCompound",
        "threshold_ppb": 1.0,
        "descriptors": ["test"],
        "smiles": "C_TEST"
    }
    predictor.db.compounds["TestCompound"] = mock_entry
    predictor.db.smiles_map["C_TEST"] = mock_entry
    predictor.db.chemical_to_descriptor["C_TEST"] = {"odt": 1.0, "descriptor": "test"}
    
    # predict_profile expects SMILES 
    profile = predictor.predict_profile({"C_TEST": 100.0})
    assert profile["TestCompound"][0] == pytest.approx(10.0)

def test_radar_aggregation():
    """Verify that intensities are correctly grouped into categories."""
    predictor = SensoryPredictor()
    
    # Mock concentrations in ppb
    # 2-Furfurylthiol (FFT) threshold ~0.01 ppb
    # Hexanal threshold ~4.5 ppb
    mock_conc = {
        "2-Furfurylthiol (FFT)": 0.1,  # OAV = 10, Intensity = sqrt(10) ~ 3.16
    }
    
    radar = predictor.get_radar_data(mock_conc)
    
    # radar[tag] is now (score, uncertainty)
    assert radar["roasted"][0] == pytest.approx(3.162, rel=1e-3)

def test_sensory_headspace_integration():
    """Verify end-to-end headspace aware sensory scoring."""
    hs = HeadspaceModel()
    predictor = SensoryPredictor(headspace=hs)
    
    # 1000 ppb (1 ppm) Hexanal in 10% fat matrix
    # Expected intensity ~0.27 (see calculations in original test)
    profile = predictor.predict_profile({"Hexanal": 1000.0}, temp_c=25.0, fat_fraction=0.1)
    # Using .get because sub-threshold compounds might be omitted from profile
    assert profile.get("Hexanal", (0.0, 0.0))[0] < 0.4
    
    # Same concentration, but 0% fat -> much higher intensity (~1.82)
    profile_pure = predictor.predict_profile({"Hexanal": 1000.0}, temp_c=25.0, fat_fraction=0.0)
    assert profile_pure.get("Hexanal", (0.0, 0.0))[0] > 1.5

if __name__ == "__main__":
    pytest.main([__file__])
