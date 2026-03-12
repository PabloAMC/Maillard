import pytest
from src.bayesian_optimizer import FormulationOptimizer
from src.inverse_design import InverseDesigner, ReactionConditions

def test_intervention_barrier_shift():
    designer = InverseDesigner(target_tag="meaty")
    cond = ReactionConditions(temperature_celsius=180)
    
    # Asparagine + Glucose (Acrylamide precursors)
    form_base = {
        "name": "Base",
        "sugars": ["glucose"],
        "amino_acids": ["asparagine"],
        "molar_ratios": {"glucose": 0.1, "asparagine": 0.1},
        "interventions": []
    }
    
    form_calcium = {
        "name": "Calcium",
        "sugars": ["glucose"],
        "amino_acids": ["asparagine"],
        "molar_ratios": {"glucose": 0.1, "asparagine": 0.1},
        "interventions": ["calcium_carbonate"]
    }
    
    res_base = designer.evaluate_single(form_base, cond)
    res_calcium = designer.evaluate_single(form_calcium, cond)
    
    # Calcium intervention should result in a lower safety penalty for acrylamide
    # (Because the barrier was increased)
    assert res_calcium.safety_score < res_base.safety_score

def test_bo_suggests_intervention():
    # Test that BO can run with interventions enabled
    optimizer = FormulationOptimizer(target_tag="meaty", risk_aversion=5.0) # High risk aversion
    
    # Mocking it to run only a few trials
    study = optimizer.optimize(
        fixed_sugars=["glucose"],
        fixed_amino_acids=["asparagine"],
        n_trials=5
    )
    
    # At least check that the trials contain interventions
    trials = study.get_trials()
    interventions = [t.params.get("intervention") for t in trials]
    assert any(i in ["calcium_carbonate", "rosemary_extract"] for i in interventions)
