import pytest

from src.bayesian_optimizer import FormulationOptimizer  # noqa: E402
from src.inverse_design import FormulationResult  # noqa: E402

def test_formulation_optimizer_initialization():
    """Verify the optimizer initializes with the correct targets."""
    optimizer = FormulationOptimizer(target_tag="meaty", minimize_tag="beany", risk_aversion=2.0)
    assert optimizer.target_tag == "meaty"
    assert optimizer.minimize_tag == "beany"
    assert optimizer.risk_aversion == 2.0
    assert optimizer.target_tag == "meaty"

def test_optimization_execution(monkeypatch):
    """
    Verify that the objective function generates correct bounds and the study 
    successfully completes. We mock the evaluate_all function to make it run instantly.
    """
    optimizer = FormulationOptimizer(target_tag="meaty")
    
    # Mock InverseDesigner.evaluate_single to return a predictable result based on temp
    # R.8: designer is now created per-trial, so we monkeypatch the class method
    def mock_evaluate_single(self, formulation, cond):
        temp = cond.temperature_celsius
        score = 100.0 - abs(temp - 150.0)
        
        return FormulationResult(
            name="MockResult",
            target_score=score,
            off_flavour_risk=0.0,
            lysine_budget=0.0,
            trapping_efficiency=0.0,
            detected_targets=["FFT"],
            detected_minimize=[],
            radar={"meaty": score},
            safety_score=0.0,
            flagged_toxics=[]
        )

    from src.inverse_design import InverseDesigner
    monkeypatch.setattr(InverseDesigner, "evaluate_single", mock_evaluate_single)
    
    # Run a short study
    study = optimizer.optimize(["ribose"], ["cysteine"], n_trials=5)
    
    assert study is not None
    assert len(study.trials) == 5
    
    # The best temperature should be close to 150 due to our mock objective
    # Note: with 5 iterations it might not perfectly hit 150, but it shouldn't crash
    best_temp = study.best_trial.params["temp"]
    assert 100.0 <= best_temp <= 200.0
    
    # Verify metadata tracking
    assert "target_score" in study.best_trial.user_attrs
    assert "safety_score" in study.best_trial.user_attrs

if __name__ == "__main__":
    pytest.main([__file__])
