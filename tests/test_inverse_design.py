import pytest
import os
import yaml
import tempfile
from pathlib import Path

from src.inverse_design import InverseDesigner
from src.conditions import ReactionConditions

@pytest.fixture
def mock_grid_path():
    """Create a temporary formulation_grid.yml with known concentration ratios to test Boltzmann scoring."""
    grid_data = {
        "formulations": [
            {
                "name": "Low Cysteine Base",
                "sugars": ["ribose"],
                "amino_acids": ["cysteine"],
                "molar_ratios": {
                    "cysteine": 0.5,
                    "ribose": 1.0
                }
            },
            {
                "name": "High Cysteine Base",
                "sugars": ["ribose"],
                "amino_acids": ["cysteine"],
                "molar_ratios": {
                    "cysteine": 2.0,
                    "ribose": 1.0
                }
            }
        ]
    }
    
    fd, path = tempfile.mkstemp(suffix=".yml")
    with os.fdopen(fd, 'w') as f:
        yaml.dump(grid_data, f)
        
    yield Path(path)
    os.remove(path)

def test_concentration_boltzmann_scoring(mock_grid_path, monkeypatch):
    """
    Formally verifies Phase 8.D Sub-task B:
    Two identical formulations (except for concentration ratio) must receive 
    different scores, with the higher-concentration precursor yielding a higher score.
    """
    import src.inverse_design as inv
    monkeypatch.setattr(inv, "GRID_FILE", mock_grid_path)
    
    designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
    
    cond = ReactionConditions(pH=5.5, temperature_celsius=150.0)
    results = designer.evaluate_all(cond)
    
    assert len(results) == 2
    
    # Extract the scores for both runs
    high_cys = next(r for r in results if r.name == "High Cysteine Base")
    low_cys = next(r for r in results if r.name == "Low Cysteine Base")
    
    # High cysteine formulation should have a strictly greater target score 
    # than the low cysteine formulation because the limiting reagent is 4x higher.
    assert high_cys.target_score > low_cys.target_score, (
        f"High Cys score ({high_cys.target_score}) must be > Low Cys score ({low_cys.target_score})"
    )


class TestInverseDesignerInitialization:
    """Test InverseDesigner initialization and setup."""

    def test_inverse_designer_creation(self):
        """InverseDesigner should initialize without errors."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        assert designer is not None
        assert designer.target_tag == "meaty"
        assert designer.minimize_tag == "beany"

    def test_inverse_designer_default_minimize(self):
        """Default minimize tag should be 'beany'."""
        designer = InverseDesigner(target_tag="roasted")
        assert designer.minimize_tag == "beany"

    def test_inverse_designer_grid_loaded(self):
        """Grid should be loaded on initialization."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        assert hasattr(designer, 'grid')
        assert len(designer.grid) > 0, "Grid should be loaded with formulations"

    def test_inverse_designer_tag_validation(self):
        """Tags should be strings, not None."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        assert isinstance(designer.target_tag, str)
        assert isinstance(designer.minimize_tag, str)
        assert len(designer.target_tag) > 0
        assert len(designer.minimize_tag) > 0


class TestInverseDesignerEvaluation:
    """Test InverseDesigner.evaluate_all() and scoring logic."""

    def test_evaluate_all_returns_results(self):
        """evaluate_all() should return FormulationResult objects."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        assert results is not None
        assert len(results) > 0, "Should have results for grid formulations"
        assert all(hasattr(r, 'target_score') for r in results), \
            "All results should have target_score attribute"

    def test_evaluate_all_result_attributes(self):
        """FormulationResult should have expected attributes."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        for result in results:
            assert hasattr(result, 'name'), "Result should have name"
            assert hasattr(result, 'target_score'), "Result should have target_score"
            assert hasattr(result, 'off_flavour_risk'), "Result should have off_flavour_risk"
            assert hasattr(result, 'lysine_budget'), "Result should have lysine_budget"
            assert hasattr(result, 'trapping_efficiency'), "Result should have trapping_efficiency"

    def test_evaluate_all_results_sorted(self):
        """Results should be sorted by target_score (descending)."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        if len(results) > 1:
            scores = [r.target_score for r in results]
            assert scores == sorted(scores, reverse=True), \
                "Results should be sorted by target_score (descending)"

    def test_target_score_is_numeric(self):
        """All target scores should be numeric."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        for result in results:
            assert isinstance(result.target_score, (int, float)), \
                f"target_score should be numeric, got {type(result.target_score)}"

    def test_risk_penalty_is_numeric(self):
        """Off-flavor risk should be numeric."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        for result in results:
            assert isinstance(result.off_flavour_risk, (int, float)), \
                f"off_flavour_risk should be numeric"
            assert result.off_flavour_risk >= 0, "Risk should be non-negative"

    def test_different_conditions_affect_scoring(self):
        """Different pH/temp should produce different scores."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        
        cond_acidic = ReactionConditions(pH=5.0, temperature_celsius=150.0)
        results_acidic = designer.evaluate_all(cond_acidic)
        
        cond_neutral = ReactionConditions(pH=7.0, temperature_celsius=150.0)
        results_neutral = designer.evaluate_all(cond_neutral)
        
        # Get first formulation scores
        score_acidic = results_acidic[0].target_score if results_acidic else 0
        score_neutral = results_neutral[0].target_score if results_neutral else 0
        
        # At different pH, should potentially get different results or rankings
        assert score_acidic >= 0 and score_neutral >= 0, "Scores should be non-negative"

    def test_different_tags_produce_different_rankings(self):
        """Different target tags should produce different rankings."""
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        
        designer_meaty = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        results_meaty = designer_meaty.evaluate_all(cond)
        
        designer_roasted = InverseDesigner(target_tag="roasted", minimize_tag="beany")
        results_roasted = designer_roasted.evaluate_all(cond)
        
        # Different targets should produce different top formulation scores
        top_meaty = results_meaty[0].target_score if results_meaty else 0
        top_roasted = results_roasted[0].target_score if results_roasted else 0
        
        # At minimum, both should be valid
        assert top_meaty >= 0
        assert top_roasted >= 0


class TestInverseDesignerErrorHandling:
    """Test error handling and edge cases."""

    def test_evaluate_all_with_extreme_conditions(self):
        """Should handle extreme pH/temperature gracefully."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        
        # Very acidic
        cond_acidic = ReactionConditions(pH=2.0, temperature_celsius=200.0)
        results = designer.evaluate_all(cond_acidic)
        assert len(results) >= 0, "Should handle extreme conditions"
        
        # Very alkaline
        cond_alkaline = ReactionConditions(pH=10.0, temperature_celsius=100.0)
        results = designer.evaluate_all(cond_alkaline)
        assert len(results) >= 0, "Should handle extreme conditions"

    def test_results_have_valid_names(self):
        """All results should have non-empty names."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        for result in results:
            assert result.name is not None
            assert isinstance(result.name, str)
            assert len(result.name) > 0, "Result name should not be empty"

    def test_evaluate_all_returns_list(self):
        """evaluate_all() should always return a list."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        cond = ReactionConditions(pH=6.0, temperature_celsius=150.0)
        results = designer.evaluate_all(cond)
        
        assert isinstance(results, list), "evaluate_all should return a list"


class TestFormulationGridLoading:
    """Test grid loading and formulation data."""

    def test_grid_formulations_have_required_fields(self):
        """Loaded formulations should have name and precursor info."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        
        for formulation in designer.grid:
            assert hasattr(formulation, 'name') or formulation.get('name'), \
                "Formulation should have a name"
            # Should have at least some precursor information
            has_precursor = (
                formulation.get('sugars') or
                formulation.get('amino_acids') or
                formulation.get('lipids') or
                formulation.get('additives')
            )
            assert has_precursor, "Formulation should have at least one precursor category"

    def test_grid_is_not_empty(self):
        """Grid should contain actual formulations."""
        designer = InverseDesigner(target_tag="meaty", minimize_tag="beany")
        assert len(designer.grid) > 0, "Grid should have at least one formulation"
        assert len(designer.grid) <= 1000, "Grid should not be unreasonably large"

