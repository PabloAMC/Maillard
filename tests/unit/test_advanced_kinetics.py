import pytest
from src.kinetics import KineticsEngine  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402

def test_ph_scaling_sigmoid():
    """Verify that pH scaling is smooth and monotonic for enolisation families."""
    KineticsEngine(temperature_k=423.15)
    
    # 1. 1,2-enolisation (Acidic favored)
    # Check monotonicity: pH 5.0 > pH 6.0 > pH 7.0
    ph_range = [5.0, 5.5, 6.0, 6.5, 7.0]
    multipliers = []
    for ph in ph_range:
        cond = ReactionConditions(pH=ph)
        multipliers.append(cond.get_ph_multiplier("Enolisation_1_2"))
    
    # Should be strictly decreasing
    assert all(x > y for x, y in zip(multipliers, multipliers[1:]))
    assert multipliers[0] > 3.0  # High at acidic
    assert multipliers[-1] < 2.0 # Low at neutral
    
    # 2. 2,3-enolisation (Alkaline favored)
    # Check monotonicity: pH 5.0 < pH 7.0 < pH 9.0
    ph_range_alk = [5.0, 6.0, 7.0, 8.0, 9.0]
    multipliers_alk = []
    for ph in ph_range_alk:
        cond = ReactionConditions(pH=ph)
        multipliers_alk.append(cond.get_ph_multiplier("Enolisation_2_3"))
        
    # Should be strictly increasing
    assert all(x < y for x, y in zip(multipliers_alk, multipliers_alk[1:]))
    assert multipliers_alk[0] < 1.0 # Suppressed at acidic
    assert multipliers_alk[-1] > 5.0 # High at alkaline

def test_ph_scaling_gaussian():
    """Verify that Schiff base formation has a peak response around pH 5.5."""
    cond_low = ReactionConditions(pH=4.0)
    cond_peak = ReactionConditions(pH=5.5)
    cond_high = ReactionConditions(pH=7.5)
    
    m_low = cond_low.get_ph_multiplier("Schiff_Base")
    m_peak = cond_peak.get_ph_multiplier("Schiff_Base")
    m_high = cond_high.get_ph_multiplier("Schiff_Base")
    
    # Peak should be highest
    assert m_peak > m_low
    assert m_peak > m_high
    assert m_peak == pytest.approx(3.0) # 1.0 + 2.0 * Gaussian(0)

def test_solvent_scaling():
    """Verify Kirkwood-Onsager solvent scaling."""
    engine = KineticsEngine(temperature_k=423.15)
    
    # Water (high epsilon) vs Lipid (low epsilon)
    cond_water = ReactionConditions(solvent_name="water")
    cond_lipid = ReactionConditions(solvent_name="lipid")
    
    # Polar transitions are faster in polar solvents (lower barrier)
    k_water = engine.get_rate_constant(25.0, conditions=cond_water)
    k_lipid = engine.get_rate_constant(25.0, conditions=cond_lipid)
    
    # f(78.4) approx 0.49
    # f(2.0) approx 0.20
    # Shift = 5.0 * (0.20 - 0.49) = -1.45 kcal/mol (barrier increases by 1.45 for lipid)
    assert k_water > k_lipid

def test_thermo_gating():
    """Verify Joback-based thermo-gating."""
    engine = KineticsEngine()
    
    # Clearly unphysical: Highly complex molecule from simple ones without energy input
    # (Using simple SMILES for Joback stability)
    reactants = ["C"] # Methane
    products = ["CCCCCCCCCC"] # Decane
    
    is_feasible, dg = engine.is_reaction_feasible(reactants, products, threshold_kcal_mol=30.0)
    
    assert not is_feasible
    assert dg > 30.0

def test_feasible_reaction():
    """Verify that feasible reactions are not gated."""
    engine = KineticsEngine()
    
    # Simple isomerisation: Acetaldehyde to Vinyl Alcohol (likely near zero or slightly endo)
    reactants = ["CC=O"]
    products = ["C=CO"]
    
    is_feasible, dg = engine.is_reaction_feasible(reactants, products, threshold_kcal_mol=30.0)
    
    assert is_feasible
    assert dg < 30.0
