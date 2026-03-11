"""
tests/unit/test_headspace.py

Verifies headspace partitioning and matrix effect logic.
"""

import pytest
from src.headspace import HeadspaceModel  # noqa: E402

def test_headspace_temperature_scaling():
    """Verify that volatility increases with temperature."""
    model = HeadspaceModel()
    # Hexanal at 25C vs 100C
    c_25 = model.predict_headspace({"Hexanal": 1.0}, 25.0)["Hexanal"]
    c_100 = model.predict_headspace({"Hexanal": 1.0}, 100.0)["Hexanal"]
    
    assert c_100 > c_25 * 5.0  # Volatility should increase significantly

def test_lipid_suppression():
    """Verify that hydrophobic compounds are suppressed by fat."""
    model = HeadspaceModel()
    matrix = {"Hexanal": 1.0, "Furfural": 1.0}
    
    # 0% fat
    air_0 = model.predict_headspace(matrix, 25.0, fat_fraction=0.0)
    # 10% fat
    air_10 = model.predict_headspace(matrix, 25.0, fat_fraction=0.1)
    
    # Hexanal is hydrophobic (Kfat ~ 450)
    # 1 / (1 + 450*0.1) = 1/46th
    assert air_10["Hexanal"] < air_0["Hexanal"] / 40.0
    
    # Furfural is polar (Kfat ~ 5)
    # 1 / (1 + 5*0.1) = 1/1.5
    assert air_10["Furfural"] > air_0["Furfural"] / 2.0
    
    # Hexanal should be MUCH more suppressed than Furfural
    ratio_hex = air_0["Hexanal"] / air_10["Hexanal"]
    ratio_fur = air_0["Furfural"] / air_10["Furfural"]
    assert ratio_hex > ratio_fur * 10.0

def test_protein_sequestration():
    """Verify that protein fraction reduces headspace concentrations."""
    model = HeadspaceModel()
    matrix = {"Methional": 1.0}
    
    # 0% protein
    air_0 = model.predict_headspace(matrix, 25.0, protein_fraction=0.0)
    # 20% protein (typical for meat analog)
    air_20 = model.predict_headspace(matrix, 25.0, protein_fraction=0.2)
    
    # Since Methional Kprot=2.0 -> 1 / (1 + 2.0*0.2) = 1/1.4
    assert air_20["Methional"] == pytest.approx(air_0["Methional"] / 1.4)

if __name__ == "__main__":
    pytest.main([__file__])
