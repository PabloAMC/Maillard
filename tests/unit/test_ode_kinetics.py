import pytest

from src.conditions import ReactionConditions
from src.ode_kinetics import simulate_kinetic_trace
from src.pathway_extractor import ElementaryStep, Species


def test_simulate_kinetic_trace_conserves_simple_mass_balance():
    reactant = Species("A", "C")
    product = Species("B", "CO")
    step = ElementaryStep(
        reactants=[reactant],
        products=[product],
        reaction_family="test_family",
        rate_constant_k=0.2,
    )

    result = simulate_kinetic_trace(
        [step],
        {},
        {reactant.smiles: 1.0},
        ReactionConditions(),
        time_minutes=1.0,
    )

    assert result.successful is True
    assert result.final_concentrations[reactant.smiles] < 1.0
    assert result.final_concentrations[product.smiles] > 0.0
    total_final = result.final_concentrations[reactant.smiles] + result.final_concentrations[product.smiles]
    assert total_final == pytest.approx(1.0, rel=1e-4)
    assert result.time_minutes[0] == 0.0
    assert result.time_minutes[-1] == pytest.approx(1.0)
    assert result.summary.reaction_count == 1
    assert result.summary.fallback_to_projection is False


def test_simulate_kinetic_trace_parallel_channels_favor_lower_barrier_path():
    precursor = Species("A", "CC=O")
    fast_product = Species("B", "C=CO")
    slow_product = Species("C", "C1CO1")
    fast_step = ElementaryStep(reactants=[precursor], products=[fast_product], reaction_family="fast")
    slow_step = ElementaryStep(reactants=[precursor], products=[slow_product], reaction_family="slow")

    result = simulate_kinetic_trace(
        [fast_step, slow_step],
        {
            f"{precursor.smiles}->{fast_product.smiles}": 12.0,
            f"{precursor.smiles}->{slow_product.smiles}": 24.0,
        },
        {precursor.smiles: 1.0},
        ReactionConditions(temperature_celsius=150.0),
        time_minutes=5.0,
    )

    assert result.final_concentrations[fast_product.smiles] > result.final_concentrations[slow_product.smiles]


def test_simulate_kinetic_trace_keeps_zero_initialized_species_at_zero_when_disconnected():
    precursor = Species("A", "CC")
    product = Species("B", "C=C")
    disconnected = Species("D", "O")
    step = ElementaryStep(reactants=[precursor], products=[product], reaction_family="test_family", rate_constant_k=0.1)

    result = simulate_kinetic_trace(
        [step],
        {},
        {precursor.smiles: 1.0, disconnected.smiles: 0.0},
        ReactionConditions(),
        time_minutes=1.0,
    )

    assert result.final_concentrations[disconnected.smiles] == pytest.approx(0.0)


def test_simulate_kinetic_trace_uses_dynamic_temperature_profile():
    precursor = Species("A", "CC=O")
    product = Species("B", "C=CO")
    step = ElementaryStep(reactants=[precursor], products=[product], reaction_family="test_family")
    barriers = {f"{precursor.smiles}->{product.smiles}": 20.0}

    cool_result = simulate_kinetic_trace(
        [step],
        barriers,
        {precursor.smiles: 1.0},
        ReactionConditions(temperature_celsius=80.0),
        time_minutes=10.0,
    )
    ramp_result = simulate_kinetic_trace(
        [step],
        barriers,
        {precursor.smiles: 1.0},
        ReactionConditions(
            temperature_celsius=80.0,
            temperature_profile=((0.0, 80.0), (10.0, 160.0)),
        ),
        time_minutes=10.0,
    )

    assert ramp_result.final_concentrations[product.smiles] > cool_result.final_concentrations[product.smiles]
    assert ramp_result.summary.used_dynamic_profiles is True