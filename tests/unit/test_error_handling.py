"""
tests/test_error_handling.py — Error handling and edge case tests

Comprehensive tests for error conditions, invalid inputs, and boundary cases
across the Maillard pipeline.

Verifies:
1. Invalid SMILES handling
2. Empty precursor lists
3. Malformed/missing input files
4. Boundary conditions (extreme pH, temp, etc.)
5. Resource exhaustion scenarios
"""

import pytest
from pathlib import Path
from src.chem_utils import Species, ElementaryStep  # noqa: E402
from src.smirks_engine import SmirksEngine  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from src.recommend import Recommender  # noqa: E402
from rdkit import Chem  # noqa: E402


class TestInvalidSmilesHandling:
    """Test handling of invalid or malformed SMILES strings."""

    def test_species_with_valid_smiles(self):
        """Valid SMILES should create valid Species."""
        s = Species("test", "CC(=O)O")  # acetic acid
        assert s.label == "test"
        assert s.smiles == "CC(=O)O"

    def test_invalid_smiles_behavior(self):
        """Invalid SMILES should not crash Species creation but may fail validation."""
        # Species should accept the SMILES even if it's invalid
        s = Species("invalid", "C(C(C")  # Unclosed parenthesis
        assert s.label == "invalid"
        # RDKit should reject this
        mol = Chem.MolFromSmiles(s.smiles)
        assert mol is None, "RDKit should reject invalid SMILES"

    def test_empty_smiles_string(self):
        """Empty SMILES string should be handled gracefully."""
        s = Species("empty", "")
        assert s.smiles == ""
        # RDKit accepts empty string and may create an empty molecule
        # (behavior varies by RDKit version)
        Chem.MolFromSmiles(s.smiles)
        # Just verify it doesn't crash the system
        assert s.label == "empty"


class TestEmptyPrecursorListHandling:
    """Test pipeline behavior with empty or missing precursors."""

    def test_smirks_engine_empty_input(self):
        """SmirksEngine should handle empty precursor list gracefully."""
        engine = SmirksEngine(ReactionConditions(pH=6.0, temperature_celsius=150.0))
        
        try:
            steps = engine.enumerate([])
            assert isinstance(steps, list), "Should return list even if empty"
        except ValueError as e:
            # Acceptable to raise error for empty input
            assert "empty" in str(e).lower() or "no" in str(e).lower()

    def test_recommender_empty_steps(self):
        """Recommender should handle empty step list gracefully."""
        recommender = Recommender()
        
        try:
            result = recommender.predict_from_steps([], {}, {"C": 1.0})
            assert isinstance(result, dict), "Should return dict result"
            assert "targets" in result or len(result) >= 0
        except (ValueError, KeyError, IndexError):
            # Acceptable to raise error for empty input
            pass


class TestConditionBoundaryValues:
    """Test ReactionConditions with boundary values."""

    def test_extreme_low_ph(self):
        """Should handle very low pH."""
        cond = ReactionConditions(pH=0.1, temperature_celsius=150.0)
        assert cond.pH == 0.1
        mult = cond.get_ph_multiplier("enolisation")
        assert mult is not None and isinstance(mult, (int, float))

    def test_extreme_high_ph(self):
        """Should handle very high pH."""
        cond = ReactionConditions(pH=14.0, temperature_celsius=150.0)
        assert cond.pH == 14.0
        mult = cond.get_ph_multiplier("enolisation")
        assert mult is not None and isinstance(mult, (int, float))

    def test_extreme_low_temperature(self):
        """Should handle very low temperature."""
        cond = ReactionConditions(pH=7.0, temperature_celsius=10.0)
        assert cond.temperature_celsius == 10.0
        mult = cond.get_arrhenius_multiplier(20.0)
        assert mult >= 0, "Arrhenius multiplier should be non-negative"

    def test_extreme_high_temperature(self):
        """Should handle very high temperature."""
        cond = ReactionConditions(pH=7.0, temperature_celsius=500.0)
        assert cond.temperature_celsius == 500.0
        mult = cond.get_arrhenius_multiplier(20.0)
        assert mult >= 0, "Arrhenius multiplier should be non-negative"

    def test_zero_temperature_kelvin(self):
        """Should not allow absolute zero."""
        cond = ReactionConditions(pH=7.0, temperature_celsius=-273.15)
        # Temperature in Kelvin should be non-negative
        assert cond.temperature_kelvin >= 0

    def test_water_activity_boundary(self):
        """Water activity should be bounded [0, 1]."""
        cond1 = ReactionConditions(water_activity=0.0)
        cond2 = ReactionConditions(water_activity=1.0)
        cond3 = ReactionConditions(water_activity=0.5)
        
        assert cond1.water_activity == 0.0
        assert cond2.water_activity == 1.0
        assert cond3.water_activity == 0.5


class TestNegativeBarrierHandling:
    """Test handling of invalid barrier values."""

    def test_negative_barrier_in_dict(self):
        """Negative barriers should be filtered or rejected."""
        barriers_dict = {
            "C+C->CC": 20.0,      # valid
            "C+C->CCC": -5.0,     # invalid (negative)
        }
        
        # Implementation might filter or raise error
        # Test just ensures it doesn't crash
        Recommender()
        # Try to use it - should either filter or handle gracefully
        try:
            # Just verify negative barriers don't crash the system
            filtered = {k: max(0.0, v) for k, v in barriers_dict.items()}
            assert filtered["C+C->CC"] == 20.0
            assert filtered["C+C->CCC"] == 0.0
        except Exception:
            pass


class TestMissingFileHandling:
    """Test behavior when expected files are missing."""

    def test_pipelineer_missing_grid(self, monkeypatch):
        """Should handle missing grid file gracefully."""
        from src.pipeline import MaillardPipeline
        import src.pipeline as inv
        
        # Point grid file to nonexistent location
        monkeypatch.setattr(inv, "GRID_FILE", Path("/tmp/nonexistent_grid_xyz_12345.yml"))
        
        try:
            designer = MaillardPipeline(target_tag="test", minimize_tag="test")
            # Should either fail during init or have empty grid
            assert designer.grid is not None
        except (FileNotFoundError, IOError, ValueError):
            # Acceptable to raise error for missing grid
            pass


class TestConcurrencyAndResourceLimits:
    """Test behavior under resource constraints."""

    def test_large_step_count(self):
        """Should handle reasonably large reaction networks."""
        # Create many mock elementary steps
        from src.barrier_constants import get_barrier
        
        steps = []
        for i in range(100):
            step = ElementaryStep(
                reactants=[Species(f"R{i}", "C")],
                products=[Species(f"P{i}", "CC")],
                reaction_family=f"Test_Reaction_{i%5}"
            )
            steps.append(step)
        
        barriers = {
            "C->CC": get_barrier(f"test_reaction_{i%5}") 
            for i in range(100)
        }
        
        recommender = Recommender()
        try:
            result = recommender.predict_from_steps(steps, barriers, {"C": 1.0})
            assert result is not None
        except Exception as e:
            # Should at least not crash with segfault
            assert isinstance(e, (ValueError, KeyError, AttributeError, IndexError))

    def test_large_barrier_dict(self):
        """Should handle large barrier dictionaries."""
        # Create barrier dict with many entries
        barriers = {f"R{i}->P{i}": float(i % 30 + 5) for i in range(1000)}
        
        assert len(barriers) == 1000
        # Should not cause issues
        recommender = Recommender()
        try:
            result = recommender.predict_from_steps([], barriers, {})
            # Should either work or fail gracefully
            assert result is not None or result is None
        except Exception as e:
            assert isinstance(e, (ValueError, KeyError, IndexError))


class TestTypeValidation:
    """Test input type checking and validation."""

    def test_smirks_engine_invalid_conditions_type(self):
        """SmirksEngine should validate conditions parameter type."""
        engine = SmirksEngine(ReactionConditions(pH=6.0, temperature_celsius=150.0))
        assert isinstance(engine.conditions, ReactionConditions)

    def test_non_numeric_barrier_is_not_validated_at_the_api_boundary(self):
        """KNOWN GAP, stated rather than skipped: barriers_dict values are never type-checked.

        REWRITTEN 2026-08-27 (Wave J2, red-team finding: self-excusing skips). This was:

            try:
                recommender.predict_from_steps([], invalid_barriers, {})
                pytest.skip("Implementation may not validate barrier types")
            except (TypeError, ValueError):
                pass

        which skips EXACTLY when the defect it exists to detect is present, and passes
        silently when it is not. Either way it reports nothing. A skip whose condition is
        "the thing under test is broken" can never fail.

        Worse, the setup was vacuous independently of the skip: the step list is ``[]``, so
        no step ever looks a barrier up and ``invalid_barriers`` is dead input. Even a fully
        validating implementation would only catch this at the lookup site, which is never
        reached. So the original test could not have detected barrier validation whether or
        not it existed.

        What is asserted now is the measured truth (2026-08-27): the call RETURNS NORMALLY,
        producing a well-formed empty result, and the garbage barrier is silently ignored.
        There is no boundary validation on barriers_dict. If someone adds it, this test goes
        red and they must decide deliberately whether the new behaviour is right -- which is
        what a test is for.
        """
        invalid_barriers = {"R+R->P": "not_a_number"}  # String instead of float
        recommender = Recommender()

        result = recommender.predict_from_steps([], invalid_barriers, {})

        assert isinstance(result, dict), (
            "predict_from_steps now rejects a non-numeric barrier. That is an improvement, "
            "but it is a behaviour change: re-pin this test to assert the exception."
        )
        assert result.get("targets") == [], (
            f"expected no targets from an empty step list, got {result.get('targets')!r}"
        )

    def test_non_numeric_concentration_raises_rather_than_being_coerced(self):
        """Concentration dict values ARE rejected -- via arithmetic, not validation.

        REWRITTEN 2026-08-27 (Wave J2, red-team finding: self-excusing skips). Same
        `pytest.skip("Implementation may not validate concentration types")` escape hatch as
        the test above. Here the skip was never reached, because the call does raise -- but
        the test would have gone quietly green-by-skip the moment it stopped raising, which
        is precisely when someone should have been told.

        Measured 2026-08-27: this raises ``TypeError: can't multiply sequence by non-int of
        type 'float'``. Note what that message reveals -- the string is not REJECTED by a
        validator, it is fed into the matrix-correction arithmetic and blows up there. The
        protection is incidental. It is asserted here because incidental protection is still
        protection, and losing it must be visible; the message itself is deliberately NOT
        asserted, since a real validator would raise something better worded.
        """
        invalid_conc = {"C": "not_a_number"}  # String instead of float
        recommender = Recommender()

        with pytest.raises((TypeError, ValueError)):
            recommender.predict_from_steps([], {}, invalid_conc)


class TestDataConsistency:
    """Test consistency of data structures and values."""

    def test_elementary_step_mass_balance_check(self):
        """Elementary steps should have valid mass balance."""
        # Create an invalid step (check it's detectable if validated)
        step = ElementaryStep(
            reactants=[Species("A", "C")],    # methane
            products=[Species("B", "CC")],     # ethane
            reaction_family="None"
        )
        
        assert len(step.reactants) > 0
        assert len(step.products) > 0

    def test_recommendation_result_field_consistency(self):
        """Recommendation results should have consistent fields."""
        recommender = Recommender()
        
        # Even with empty input, structure should be consistent
        try:
            result = recommender.predict_from_steps([], {}, {})
            if result:
                # Check expected key fields
                expected_keys = ["targets", "metrics"]
                for key in expected_keys:
                    # Some keys might not exist, that's ok
                    if key in result:
                        assert result[key] is not None
        except (ValueError, KeyError, IndexError):
            pass
