from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from src.conditions import ReactionConditions
from src.pathway_extractor import ElementaryStep
from src.projection import _canon


@dataclass(frozen=True)
class SpeciesState:
    index: int
    canonical_smiles: str
    label: str
    initial_concentration: float


@dataclass(frozen=True)
class ReactionChannel:
    reactant_indices: Tuple[int, ...]
    product_indices: Tuple[int, ...]
    reactant_stoichiometry: Tuple[Tuple[int, int], ...]
    product_stoichiometry: Tuple[Tuple[int, int], ...]
    reaction_family: str
    barrier_kcal: Optional[float]
    barrier_uncertainty_kcal: float
    source_quality: str
    base_rate_constant: float


@dataclass(frozen=True)
class KineticTrace:
    time_minutes: List[float]
    concentrations_by_species: Dict[str, List[float]]
    final_concentrations: Dict[str, float]
    peak_concentrations: Dict[str, float]
    integrated_concentrations: Dict[str, float]


@dataclass(frozen=True)
class KineticRunSummary:
    solver_name: str
    success: bool
    time_horizon_minutes: float
    species_count: int
    reaction_count: int
    valid_channel_count: int
    tracked_species: int
    used_dynamic_profiles: bool
    used_pruning: bool
    concentration_floor: float
    solver_fallback_used: bool = False
    fallback_to_projection: bool = False
    fallback_reason: Optional[str] = None


@dataclass(frozen=True)
class KineticSimulationResult:
    species_states: Tuple[SpeciesState, ...]
    reaction_channels: Tuple[ReactionChannel, ...]
    trace: KineticTrace
    summary: KineticRunSummary

    @property
    def time_minutes(self) -> List[float]:
        return self.trace.time_minutes

    @property
    def final_concentrations(self) -> Dict[str, float]:
        return self.trace.final_concentrations

    @property
    def peak_concentrations(self) -> Dict[str, float]:
        return self.trace.peak_concentrations

    @property
    def integrated_concentrations(self) -> Dict[str, float]:
        return self.trace.integrated_concentrations

    @property
    def trajectories(self) -> Dict[str, List[float]]:
        return self.trace.concentrations_by_species

    @property
    def successful(self) -> bool:
        return self.summary.success

    @property
    def solver(self) -> str:
        return self.summary.solver_name


def _species_key(smiles: str) -> str:
    return _canon(smiles) or smiles


def _reaction_key(step: ElementaryStep) -> str:
    reactants = "+".join(sorted(_species_key(species.smiles) for species in step.reactants))
    products = "+".join(sorted(_species_key(species.smiles) for species in step.products))
    return f"{reactants}->{products}"


def _stoichiometry(species_list, index_by_species: Dict[str, int]) -> Dict[int, int]:
    counts: Dict[int, int] = {}
    for species in species_list:
        canon = _species_key(species.smiles)
        idx = index_by_species[canon]
        counts[idx] = counts.get(idx, 0) + 1
    return counts


def _resolve_barrier(step: ElementaryStep, barriers_dict: Dict[str, float]) -> Tuple[Optional[float], float]:
    barrier_data = barriers_dict.get(_reaction_key(step))
    if barrier_data is None:
        return None, float(step.barrier_uncertainty_kcal)
    if isinstance(barrier_data, tuple):
        return float(barrier_data[0]), float(barrier_data[1])
    return float(barrier_data), float(step.barrier_uncertainty_kcal)


def _resolve_rate_constant(
    step: ElementaryStep,
    barriers_dict: Dict[str, float],
    conditions: ReactionConditions,
) -> float:
    if step.rate_constant_k is not None and step.rate_constant_k > 0.0:
        return float(step.rate_constant_k)

    barrier, _ = _resolve_barrier(step, barriers_dict)
    if barrier is None:
        return float(conditions.get_rate_constant(step.reaction_family or "unknown"))
    return float(
        conditions.get_rate_constant(
            step.reaction_family or "unknown",
            ea_override_kcal=barrier,
        )
    )


def build_kinetic_system(
    steps: List[ElementaryStep],
    barriers_dict: Dict[str, float],
    initial_concentrations: Dict[str, float],
    conditions: ReactionConditions,
) -> Tuple[Tuple[SpeciesState, ...], Tuple[ReactionChannel, ...], np.ndarray]:
    species_catalog: Dict[str, str] = {}
    for smiles in initial_concentrations:
        canon = _species_key(smiles)
        species_catalog.setdefault(canon, canon)
    for step in steps:
        for species in [*step.reactants, *step.products]:
            canon = _species_key(species.smiles)
            species_catalog.setdefault(canon, species.label or canon)

    ordered_species = sorted(species_catalog)
    index_by_species = {canon: idx for idx, canon in enumerate(ordered_species)}
    initial_state = np.zeros(len(ordered_species), dtype=float)
    species_states: List[SpeciesState] = []
    for canon in ordered_species:
        initial_concentration = 0.0
        for smiles, concentration in initial_concentrations.items():
            if _species_key(smiles) == canon:
                initial_concentration = max(float(concentration), 0.0)
                break
        initial_state[index_by_species[canon]] = initial_concentration
        species_states.append(
            SpeciesState(
                index=index_by_species[canon],
                canonical_smiles=canon,
                label=species_catalog.get(canon, canon),
                initial_concentration=initial_concentration,
            )
        )

    reaction_channels: List[ReactionChannel] = []
    for step in steps:
        reactants = _stoichiometry(step.reactants, index_by_species)
        products = _stoichiometry(step.products, index_by_species)
        if not reactants and not products:
            continue
        barrier_kcal, barrier_uncertainty_kcal = _resolve_barrier(step, barriers_dict)
        base_rate_constant = _resolve_rate_constant(step, barriers_dict, conditions)
        if base_rate_constant <= 0.0:
            continue
        reaction_channels.append(
            ReactionChannel(
                reactant_indices=tuple(sorted(reactants)),
                product_indices=tuple(sorted(products)),
                reactant_stoichiometry=tuple(sorted((idx, int(coeff)) for idx, coeff in reactants.items())),
                product_stoichiometry=tuple(sorted((idx, int(coeff)) for idx, coeff in products.items())),
                reaction_family=step.reaction_family or "unknown",
                barrier_kcal=barrier_kcal,
                barrier_uncertainty_kcal=barrier_uncertainty_kcal,
                source_quality=step.source_quality,
                base_rate_constant=float(base_rate_constant),
            )
        )
    return tuple(species_states), tuple(reaction_channels), initial_state


def simulate_kinetic_trace(
    steps: List[ElementaryStep],
    barriers_dict: Dict[str, float],
    initial_concentrations: Dict[str, float],
    conditions: ReactionConditions,
    time_minutes: float,
) -> KineticSimulationResult:
    species_states, reaction_channels, initial_state = build_kinetic_system(
        steps,
        barriers_dict,
        initial_concentrations,
        conditions,
    )
    if not species_states:
        empty_trace = KineticTrace([], {}, {}, {}, {})
        empty_summary = KineticRunSummary(
            solver_name="solve_ivp_bdf",
            success=True,
            time_horizon_minutes=max(float(time_minutes or 0.0), 0.0),
            species_count=0,
            reaction_count=0,
            valid_channel_count=0,
            tracked_species=0,
            used_dynamic_profiles=conditions.has_dynamic_profile,
            used_pruning=False,
            concentration_floor=1.0e-18,
            fallback_to_projection=True,
            fallback_reason="empty_species_network",
        )
        return KineticSimulationResult(tuple(), tuple(), empty_trace, empty_summary)

    total_seconds = max(float(time_minutes or 0.0), 0.0) * 60.0
    if total_seconds <= 0.0 or not reaction_channels:
        time_axis = np.array([0.0, max(float(time_minutes or 0.0), 0.0)], dtype=float)
        concentrations = np.repeat(initial_state[:, None], len(time_axis), axis=1)
        fallback_reason = "no_valid_reaction_channels" if not reaction_channels else None
        return _build_result(
            species_states,
            reaction_channels,
            time_axis,
            concentrations,
            solver_name="solve_ivp_bdf",
            success=True,
            used_dynamic_profiles=conditions.has_dynamic_profile,
            concentration_floor=1.0e-18,
            fallback_to_projection=not reaction_channels,
            fallback_reason=fallback_reason,
        )

    point_count = max(2, min(25, int(max(time_minutes, 1.0)) + 1))
    evaluation_times = np.linspace(0.0, total_seconds, point_count)
    concentration_floor = 1.0e-18

    def rhs(_time_seconds: float, state: np.ndarray) -> np.ndarray:
        non_negative_state = np.maximum(state, 0.0)
        delta = np.zeros_like(non_negative_state)
        current_time_minutes = _time_seconds / 60.0
        for channel in reaction_channels:
            if conditions.has_dynamic_profile:
                rate = conditions.get_rate_constant(
                    channel.reaction_family,
                    ea_override_kcal=channel.barrier_kcal,
                    time_minutes=current_time_minutes,
                )
            else:
                rate = channel.base_rate_constant
            for idx, coefficient in channel.reactant_stoichiometry:
                concentration = non_negative_state[idx]
                if concentration <= concentration_floor:
                    rate = 0.0
                    break
                rate *= concentration ** coefficient
            if rate <= 0.0:
                continue
            for idx, coefficient in channel.reactant_stoichiometry:
                delta[idx] -= coefficient * rate
            for idx, coefficient in channel.product_stoichiometry:
                delta[idx] += coefficient * rate
        return delta

    solution = solve_ivp(
        rhs,
        (0.0, total_seconds),
        initial_state,
        method="BDF",
        t_eval=evaluation_times,
        atol=1e-9,
        rtol=1e-6,
    )
    solver_name = "solve_ivp_bdf"
    solver_fallback_used = False
    if not solution.success or solution.y.size == 0:
        solution = solve_ivp(
            rhs,
            (0.0, total_seconds),
            initial_state,
            method="RK45",
            t_eval=evaluation_times,
            atol=1e-9,
            rtol=1e-6,
        )
        solver_name = "solve_ivp_rk45"
        solver_fallback_used = True
    if not solution.success or solution.y.size == 0:
        concentrations = np.repeat(initial_state[:, None], len(evaluation_times), axis=1)
        return _build_result(
            species_states,
            reaction_channels,
            evaluation_times / 60.0,
            concentrations,
            solver_name=solver_name,
            success=False,
            used_dynamic_profiles=conditions.has_dynamic_profile,
            concentration_floor=concentration_floor,
            solver_fallback_used=solver_fallback_used,
            fallback_to_projection=True,
            fallback_reason="solver_failure",
        )

    concentrations = np.maximum(solution.y, 0.0)
    return _build_result(
        species_states,
        reaction_channels,
        solution.t / 60.0,
        concentrations,
        solver_name=solver_name,
        success=True,
        used_dynamic_profiles=conditions.has_dynamic_profile,
        concentration_floor=concentration_floor,
        solver_fallback_used=solver_fallback_used,
    )


def _build_result(
    species_states: Tuple[SpeciesState, ...],
    reaction_channels: Tuple[ReactionChannel, ...],
    time_axis_minutes: np.ndarray,
    concentrations: np.ndarray,
    *,
    solver_name: str,
    success: bool,
    used_dynamic_profiles: bool,
    concentration_floor: float,
    solver_fallback_used: bool = False,
    fallback_to_projection: bool = False,
    fallback_reason: Optional[str] = None,
) -> KineticSimulationResult:
    time_grid = [float(value) for value in time_axis_minutes]
    concentration_by_species: Dict[str, List[float]] = {}
    final_concentrations: Dict[str, float] = {}
    peak_concentrations: Dict[str, float] = {}
    integrated_concentrations: Dict[str, float] = {}
    tracked_species = 0
    for species in species_states:
        trace = [float(value) for value in concentrations[species.index]]
        concentration_by_species[species.canonical_smiles] = trace
        final_concentrations[species.canonical_smiles] = trace[-1]
        peak_concentrations[species.canonical_smiles] = max(trace)
        integrated_concentrations[species.canonical_smiles] = float(np.trapezoid(trace, time_grid)) if len(time_grid) >= 2 else float(trace[-1])
        if peak_concentrations[species.canonical_smiles] > concentration_floor:
            tracked_species += 1
    trace = KineticTrace(
        time_minutes=time_grid,
        concentrations_by_species=concentration_by_species,
        final_concentrations=final_concentrations,
        peak_concentrations=peak_concentrations,
        integrated_concentrations=integrated_concentrations,
    )
    summary = KineticRunSummary(
        solver_name=solver_name,
        success=success,
        time_horizon_minutes=float(time_grid[-1]) if time_grid else 0.0,
        species_count=len(species_states),
        reaction_count=len(reaction_channels),
        valid_channel_count=len(reaction_channels),
        tracked_species=tracked_species,
        used_dynamic_profiles=used_dynamic_profiles,
        used_pruning=tracked_species < len(species_states),
        concentration_floor=concentration_floor,
        solver_fallback_used=solver_fallback_used,
        fallback_to_projection=fallback_to_projection,
        fallback_reason=fallback_reason,
    )
    return KineticSimulationResult(species_states, reaction_channels, trace, summary)