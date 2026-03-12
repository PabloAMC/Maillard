import pytest
import sys
from pathlib import Path
from src.inverse_design import InverseDesigner
from src.conditions import ReactionConditions

def test_peptide_accessibility_blind_spot():
    """
    BLIND SPOT: The engine assumes 100% amino acid availability.
    In raw pea protein, amino acids are buried. 
    EXPECTATION: High hydrolysis should increase flavor score.
    CURRENT: Score remains constant regardless of 'hydrolysis' state.
    """
    designer = InverseDesigner(target_tag="meaty")
    cond = ReactionConditions(temperature_celsius=100)
    
    # Formulation with a hypothetical 'degree_of_hydrolysis'
    low_hydrolysis = {
        "name": "IntactProtein",
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"lysine": 1.0, "cysteine": 0.5},
        "degree_of_hydrolysis": 0.1 # 10%
    }
    
    high_hydrolysis = {
        "name": "HydrolyzedProtein",
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"lysine": 1.0, "cysteine": 0.5},
        "degree_of_hydrolysis": 0.9 # 90%
    }
    
    res_low = designer.evaluate_single(low_hydrolysis, cond)
    res_high = designer.evaluate_single(high_hydrolysis, cond)
    
    # This currently fails because the engine doesn't look at 'degree_of_hydrolysis'
    # We use xfail to document this gap.
    pytest.xfail("Engine does not yet support peptide accessibility scaling.")
    assert res_high.target_score > res_low.target_score

def test_matrix_inhibition_blind_spot():
    """
    BLIND SPOT: Volatiles are 'trapped' by fiber/starch.
    EXPECTATION: High fiber content should decrease sensory radar scores.
    CURRENT: Radar scores depend only on chemical concentration.
    """
    designer = InverseDesigner(target_tag="meaty")
    cond_clear = ReactionConditions(protein_fraction=1.0) # Pure solution
    cond_matrix = ReactionConditions(protein_fraction=0.1, matrix_fiber=0.5) # High bread/pea matrix
    
    form = {
        "name": "BaseForm",
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"lysine": 1.0, "cysteine": 0.5}
    }
    
    res_clear = designer.evaluate_single(form, cond_clear)
    res_matrix = designer.evaluate_single(form, cond_matrix)
    
    pytest.xfail("Engine does not yet support matrix inhibition (volatile partitioning).")
    assert res_matrix.radar["meaty"][0] < res_clear.radar["meaty"][0]

def test_metal_catalysis_blind_spot():
    """
    BLIND SPOT: Non-heme iron (common in pea) accelerates pyrazines.
    EXPECTATION: Presence of Iron should lower pyrazine barriers.
    CURRENT: Only temperature and pH are taken into account.
    """
    designer = InverseDesigner(target_tag="roasted") # Pyrazines
    cond_no_iron = ReactionConditions(temperature_celsius=120)
    cond_iron = ReactionConditions(temperature_celsius=120, metal_catalyst="Fe2+")
    
    form = {
        "name": "NuttyForm",
        "amino_acids": ["glycine", "alanine"],
        "sugars": ["glucose"]
    }
    
    res1 = designer.evaluate_single(form, cond_no_iron)
    res2 = designer.evaluate_single(form, cond_iron)
    
    pytest.xfail("Engine does not yet support non-heme metal catalysis (Iron).")
    assert res2.target_score > res1.target_score
