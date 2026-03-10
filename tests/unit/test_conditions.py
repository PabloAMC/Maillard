import pytest
import math
from src.conditions import ReactionConditions

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
