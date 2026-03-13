import pytest
from src.conditions import ReactionConditions  # noqa: E402

def test_ph_multipliers():
    # Acidic conditions
    cond_acid = ReactionConditions(pH=5.0)
    assert cond_acid.get_ph_multiplier("1,2-enolisation") > 1.0
    assert cond_acid.get_ph_multiplier("2,3-enolisation") < 1.0
    
    # Alkaline conditions
    cond_alkaline = ReactionConditions(pH=8.0)
    # With smooth sigmoids, we use approx for the tails
    assert cond_alkaline.get_ph_multiplier("1,2-enolisation") == pytest.approx(1.0, abs=0.1)
    assert cond_alkaline.get_ph_multiplier("2,3-enolisation") > 4.0 # Heavily favored

def test_water_activity():
    assert ReactionConditions(water_activity=0.7).get_water_activity_multiplier() == 1.0
    assert ReactionConditions(water_activity=0.2).get_water_activity_multiplier() < 0.5
    assert ReactionConditions(water_activity=0.95).get_water_activity_multiplier() < 0.5

def test_arrhenius():
    cond_hot = ReactionConditions(temperature_celsius=200.0) # 473.15 K
    cond_cold = ReactionConditions(temperature_celsius=25.0) # 298.15 K
    
    # For a 15 kcal/mol barrier
    rate_hot = cond_hot.get_arrhenius_multiplier(15.0)
    rate_cold = cond_cold.get_arrhenius_multiplier(15.0)
    
    assert rate_hot > rate_cold

def test_heme_catalysis():
    """Verify that heme catalyst increases multipliers for specific families."""
    cond_none = ReactionConditions(metal_catalyst=None)
    cond_heme = ReactionConditions(metal_catalyst="heme")
    
    # 1. Lipid Oxidation promotion
    mult_none = cond_none.get_ph_multiplier("Lipid_Oxidation")
    mult_heme = cond_heme.get_ph_multiplier("Lipid_Oxidation")
    assert mult_heme == pytest.approx(mult_none * 5.0)
    
    # 2. Pyrazine promotion
    mult_none_p = cond_none.get_ph_multiplier("Pyrazine_Formation")
    mult_heme_p = cond_heme.get_ph_multiplier("Pyrazine_Formation")
    assert mult_heme_p == pytest.approx(mult_none_p * 1.5)
    
    # 3. Non-catalyzed family (e.g. Schiff base) should NOT be boosted by heme specifically beyond pH effect
    mult_none_s = cond_none.get_ph_multiplier("Schiff_Base")
    mult_heme_s = cond_heme.get_ph_multiplier("Schiff_Base")
    assert mult_heme_s == mult_none_s
