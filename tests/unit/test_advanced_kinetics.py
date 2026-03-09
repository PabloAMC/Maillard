import pytest
import numpy as np
from src.kinetics import KineticsEngine
from src.conditions import ReactionConditions

def test_ph_scaling():
    """Verify that pH scaling affects rate constants for specific families."""
    engine = KineticsEngine(temperature_k=423.15)
    
    # 1,2-enolisation (Favored at pH < 6)
    cond_acid = ReactionConditions(pH=5.0)
    cond_neutral = ReactionConditions(pH=7.0)
    
    k_acid = engine.get_rate_constant(25.0, conditions=cond_acid, reaction_family="Enolisation_1_2")
    k_neutral = engine.get_rate_constant(25.0, conditions=cond_neutral, reaction_family="Enolisation_1_2")
    
    # Should be 5x faster in acidic conditions per get_ph_multiplier
    assert k_acid / k_neutral == pytest.approx(5.0, rel=0.1)

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
