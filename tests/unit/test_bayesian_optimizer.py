import pytest

from src.bayesian_optimizer import FormulationOptimizer  # noqa: E402
from src.pipeline import FormulationResult  # noqa: E402

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
    
    # Mock MaillardPipeline.evaluate_single to return a predictable result based on temp
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
            flagged_toxics=[],
            texture_risk=0.0
        )

    from src.pipeline import MaillardPipeline
    monkeypatch.setattr(MaillardPipeline, "evaluate_single", mock_evaluate_single)
    
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


def test_meaty_quality_penalty_reduces_objective(monkeypatch):
    optimizer = FormulationOptimizer(target_tag="meaty")

    class DummyTrial:
        number = 0

        def __init__(self):
            self.user_attrs = {}
            # Mirror optuna.Trial: sampled values are readable off `.params`,
            # which is how the objective now rebuilds the formulation
            # (src.bayesian_optimizer.formulation_from_params).
            self.params = {}

        def suggest_float(self, name, low, high, log=False):
            values = {
                "sugar_conc": 0.1,
                "aa_conc_sulfur": 0.1,
                "aa_conc_branched": 0.1,
                "aa_conc_basic": 0.1,
                "aa_conc_other": 0.1,
                "ph": 5.5,
                "temp": 150.0,
                "aw": 0.5,
                "time_minutes": 30.0,
                "intervention_dose": 0.0,
            }
            self.params[name] = values[name]
            return values[name]

        def suggest_categorical(self, name, options):
            values = {
                "intervention_agent": "none",
                "pre_processing": "none",
            }
            self.params[name] = values[name]
            return values[name]

        def set_user_attr(self, key, value):
            self.user_attrs[key] = value

    def mock_evaluate_single(self, formulation, cond):
        return FormulationResult(
            name="QualityPenaltyProbe",
            target_score=10.0,
            off_flavour_risk=0.0,
            safety_score=0.0,
            meaty_quality_penalty=1.75,
            mft_to_furfural_ratio=0.001,
            avg_uncertainty=0.0,
        )

    from src.pipeline import MaillardPipeline
    monkeypatch.setattr(MaillardPipeline, "evaluate_single", mock_evaluate_single)

    trial = DummyTrial()
    value = optimizer.objective(trial, ["ribose"], ["cysteine"], None)

    assert value == pytest.approx(8.25)
    assert trial.user_attrs["meaty_quality_penalty"] == pytest.approx(1.75)
    assert trial.user_attrs["mft_to_furfural_ratio"] == pytest.approx(0.001)


def test_furanone_penalty_reduces_objective(monkeypatch):
    optimizer = FormulationOptimizer(target_tag="meaty")

    class DummyTrial:
        number = 0

        def __init__(self):
            self.user_attrs = {}
            # Mirror optuna.Trial: sampled values are readable off `.params`,
            # which is how the objective now rebuilds the formulation
            # (src.bayesian_optimizer.formulation_from_params).
            self.params = {}

        def suggest_float(self, name, low, high, log=False):
            values = {
                "sugar_conc": 0.1,
                "aa_conc_sulfur": 0.1,
                "aa_conc_branched": 0.1,
                "aa_conc_basic": 0.1,
                "aa_conc_other": 0.1,
                "ph": 5.5,
                "temp": 150.0,
                "aw": 0.5,
                "time_minutes": 30.0,
                "intervention_dose": 0.0,
            }
            self.params[name] = values[name]
            return values[name]

        def suggest_categorical(self, name, options):
            values = {
                "intervention_agent": "none",
                "pre_processing": "none",
            }
            self.params[name] = values[name]
            return values[name]

        def set_user_attr(self, key, value):
            self.user_attrs[key] = value

    def mock_evaluate_single(self, formulation, cond):
        return FormulationResult(
            name="FuranonePenaltyProbe",
            target_score=10.0,
            off_flavour_risk=0.0,
            safety_score=0.0,
            furanone_penalty=0.35,
            avg_uncertainty=0.0,
        )

    from src.pipeline import MaillardPipeline
    monkeypatch.setattr(MaillardPipeline, "evaluate_single", mock_evaluate_single)

    trial = DummyTrial()
    value = optimizer.objective(trial, ["ribose"], ["alanine"], None)

    assert value == pytest.approx(9.65)
    assert trial.user_attrs["furanone_penalty"] == pytest.approx(0.35)

if __name__ == "__main__":
    pytest.main([__file__])
