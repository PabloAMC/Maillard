import pytest

from src.conditions import ReactionConditions


def test_reaction_conditions_interpolate_dynamic_profiles():
    conditions = ReactionConditions(
        pH=6.0,
        temperature_celsius=80.0,
        water_activity=0.8,
        temperature_profile=((0.0, 80.0), (10.0, 180.0)),
        water_activity_profile=((0.0, 0.8), (10.0, 0.4)),
    )

    assert conditions.temperature_celsius_at(5.0) == pytest.approx(130.0)
    assert conditions.water_activity_at(5.0) == pytest.approx(0.6)


def test_reaction_conditions_reject_unsorted_profiles():
    with pytest.raises(ValueError, match="strictly increasing"):
        ReactionConditions(
            temperature_profile=((0.0, 80.0), (5.0, 120.0), (5.0, 140.0)),
        )


def test_reaction_conditions_require_initial_and_terminal_points():
    with pytest.raises(ValueError, match="initial point and a terminal point"):
        ReactionConditions(
            temperature_profile=((0.0, 80.0),),
        )