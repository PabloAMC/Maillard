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
    
    # threshold = 1.0, conc = 100.0 -> OAV = 100 -> Intensity = 10
    # Custom mock DB entry for predictable numbers
    mock_entry = {
        "name": "TestCompound",
        "threshold_ppm": 1.0,
        "descriptors": ["test"],
        "smiles": "C"
    }
    predictor.db.compounds["TestCompound"] = mock_entry
    
    profile = predictor.predict_profile({"TestCompound": 100.0})
    # profile["TestCompound"] is now (intensity, uncertainty)
    assert profile["TestCompound"][0] == pytest.approx(10.0)

def test_radar_aggregation():
    """Verify that intensities are correctly grouped into categories."""
    predictor = SensoryPredictor()
    
    # Mock data
    # 2-Furfurylthiol (FFT) is in 'roasted' tag in sensory_tags.yml
    # Hexanal is in 'beany' tag
    
    mock_conc = {
        "2-Furfurylthiol (FFT)": 0.0001, # OAV = 10, Intensity = sqrt(10) ~ 3.16
        "Hexanal": 0.045                # OAV = 10, Intensity = sqrt(10) ~ 3.16
    }
    
    radar = predictor.get_radar_data(mock_conc)
    
    # radar[tag] is now (score, uncertainty)
    assert radar["roasted"][0] == pytest.approx(3.162, rel=1e-3)
    assert radar["beany"][0] == pytest.approx(3.162, rel=1e-3)

def test_sensory_headspace_integration():
    """Verify end-to-end headspace aware sensory scoring."""
    hs = HeadspaceModel()
    predictor = SensoryPredictor(headspace=hs)
    
    # 1.0 ppm Hexanal in 10% fat matrix
    # Kaw(25C) = 0.015, Kfat = 450
    # Kaw_eff = 0.015 / (1 + 450*0.1) = 0.015 / 46 = 0.000326
    # conc_air = 1.0 * 0.000326 = 3.26e-4
    # Hexanal ODT = 4.5 ppb = 0.0045 ppm
    # OAV = 3.26e-4 / 4.5e-3 = 0.072
    # Intensity = sqrt(0.072) ~ 0.27
    
    profile = predictor.predict_profile({"Hexanal": 1.0}, temp_c=25.0, fat_fraction=0.1)
    assert profile["Hexanal"][0] < 0.3
    
    # Same concentration, but 0% fat -> much higher intensity
    profile_pure = predictor.predict_profile({"Hexanal": 1.0}, temp_c=25.0, fat_fraction=0.0)
    # Kaw = 0.015, conc_air = 0.015, OAV = 0.015 / 0.0045 = 3.33, Intensity = sqrt(3.33) ~ 1.82
    assert profile_pure["Hexanal"][0] > 1.5

if __name__ == "__main__":
    pytest.main([__file__])
