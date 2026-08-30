from dataclasses import dataclass
from typing import Dict, Any, Optional

from src.chem_utils import canonicalize_smiles

@dataclass(frozen=True)
class ProjectionBudget:
    limiting_precursor_molar: float
    load_factor: float
    temperature_factor: float
    time_factor: float
    severity: float
    volatile_yield_fraction: float
    total_volatile_budget_molar: float
    limiting_precursor_name: str
    kinetic_drive: float = 0.0
    conversion_extent: float = 0.0


@dataclass(frozen=True)
class ProjectionStrategy:
    name: str
    precursor_concentration_unit: str
    limiting_pool_to_molar_factor: float
    baseline_volatile_yield_fraction: float
    conversion_ceiling_fraction: float
    apparent_activation_energy_kj_mol: float
    reference_temperature_kelvin: float
    reference_conversion_time_min: float
    ppb_conversion_factor: float
    ppb_basis: str
    notes: str


# Residence time assumed when a formulation does not pin one. Mirrors the
# `float(time_minutes or 60.0)` fallback already used for matrix process-state
# classification in recommend.py, so the two layers agree on "unspecified hold".
DEFAULT_NOMINAL_TIME_MINUTES = 60.0

_GAS_CONSTANT_J_PER_MOL_K = 8.314462618


DEFAULT_PROJECTION_STRATEGY = ProjectionStrategy(
    name="precursor_limited_observable_v2_arrhenius",
    precursor_concentration_unit="mM",
    limiting_pool_to_molar_factor=1.0e-3,
    # Fitted 2026-08-27 against literature-sourced benchmark rows only; see
    # scripts/generators/refit_projection_constants.py for the objective and the
    # reproducible fit record (results/validation/projection_constant_refit.json).
    baseline_volatile_yield_fraction=1.0e-6,
    # Mass-conservation ceiling: the volatile budget can never exceed the limiting
    # precursor pool it is allocated from.
    conversion_ceiling_fraction=1.0,
    # Apparent activation energy for Maillard volatile formation. Anchored to the
    # `enolisation` entry in data/lit/arrhenius_params.yml (120 kJ/mol) -- the
    # Amadori-enolisation gateway that every downstream volatile route passes
    # through -- and inside the 97-160 kJ/mol envelope spanned by the other
    # volatile-forming steps in that file (schiff 97, dehydration 124, pyrazine 138,
    # retro-aldol 160). Held FIXED during the refit: with ~10 literature rows there
    # is not enough temperature leverage to fit it without overfitting, and pinning
    # it to literature is what makes the temperature dependence honest.
    apparent_activation_energy_kj_mol=120.0,
    reference_temperature_kelvin=423.15,  # 150 C, the panel's high-temperature anchor
    # RETAINED 2026-08-27 (Wave H) — deliberately NOT re-fitted. After the Wave G1
    # chemistry rebuild the refit optimum moved to 2512 min (objective 0.8777 -> 0.7543
    # dex over the literature rows), but tau_ref is a single GLOBAL SCALE on the volatile
    # budget: at that optimum the Hofmann1998 MFT residual collapses to 0.048 dex while
    # Resconi furfural grows to 19x OVER. The budget is already the right order of
    # magnitude (~1050 ppb at the Hofmann conditions against 542 ppb of measured
    # MFT+FFT); the sulfur deficit is an ALLOCATION deficit, and closing it with a global
    # multiplier would hide that finding rather than fix it. Full reasoning and both
    # objectives: results/validation/projection_constant_refit.md.
    reference_conversion_time_min=1.2589e4,  # fitted (see above): 10**4.10 min
    ppb_conversion_factor=1.0e6,
    ppb_basis="aqueous_mass_equivalent_ppb",
    notes=(
        "Allocates a conservative volatile budget from the limiting precursor pool, then converts "
        "M to ppb via MW assuming dilute aqueous density (~1 kg/L) before matrix/headspace projection. "
        "The budget's thermal dependence is a first-order conversion extent driven by an apparent "
        "Arrhenius rate (v2, 2026-08-27); it replaces the v1 severity sigmoid, which saturated at "
        "~1.11x its 150 C value by 190 C and so capped all high-temperature chemistry."
    ),
)


def _canon(smi: str) -> str:
    return canonicalize_smiles(smi, fallback_to_original=True, strip_salts=True) or smi

def _thermal_severity(temperature_kelvin: float, time_minutes: Optional[float]) -> float:
    return _projection_temperature_factor(temperature_kelvin) * _projection_time_factor(time_minutes)


def _projection_temperature_factor(temperature_kelvin: float) -> float:
    """Bounded process-state severity index in (0, 1).

    NOTE (2026-08-27 retune): this sigmoid used to set the SIZE of the volatile
    budget, which was wrong -- it saturates by construction (0.966 at 170 C, 0.988
    at 190 C), so the budget grew only 1.11x from 150 C to 190 C while the model's
    own Arrhenius drive grows ~20x over the same span. The visible symptom was a
    spurious furfural maximum near 145-150 C: a nearly-constant budget divided by a
    still-growing softmax denominator.

    It survives only as what it always actually was -- a dimensionless "how hot a
    process is this" descriptor, consumed by melanoidin trapping, the low-severity
    depth bias and the direct-sulfur bonus in recommend.py. Those are state
    heuristics that SHOULD saturate. The budget's thermal dependence now lives in
    `_projection_kinetic_drive`.
    """
    import math

    temp_c = temperature_kelvin - 273.15
    return 1.0 / (1.0 + math.exp(-(temp_c - 110.0) / 18.0))


def _projection_kinetic_drive(
    temperature_kelvin: float,
    time_minutes: Optional[float],
    strategy: "ProjectionStrategy" = None,  # type: ignore[assignment]
) -> float:
    """Dimensionless first-order conversion drive k(T) * t for the volatile budget.

    k(T) is an apparent Arrhenius rate referenced to `reference_temperature_kelvin`
    with the literature apparent activation energy carried on the strategy, so the
    budget tracks the same exponential temperature dependence as the reaction
    network it is allocating over instead of an independent saturating sigmoid.
    """
    import math

    if strategy is None:
        strategy = DEFAULT_PROJECTION_STRATEGY
    if temperature_kelvin <= 0.0:
        return 0.0
    effective_time = (
        DEFAULT_NOMINAL_TIME_MINUTES if time_minutes is None else max(float(time_minutes), 0.0)
    )
    if effective_time <= 0.0:
        return 0.0
    tau_ref = max(float(strategy.reference_conversion_time_min), 1.0e-12)
    ea_j = float(strategy.apparent_activation_energy_kj_mol) * 1.0e3
    t_ref = float(strategy.reference_temperature_kelvin)
    exponent = -(ea_j / _GAS_CONSTANT_J_PER_MOL_K) * (1.0 / temperature_kelvin - 1.0 / t_ref)
    # Guard the extremes: the network itself is only meaningful over ~room
    # temperature to ~300 C, and exp() overflow here would be silent nonsense.
    exponent = max(min(exponent, 60.0), -60.0)
    return (effective_time / tau_ref) * math.exp(exponent)


def _projection_conversion_extent(
    temperature_kelvin: float,
    time_minutes: Optional[float],
    strategy: "ProjectionStrategy" = None,  # type: ignore[assignment]
) -> float:
    """First-order extent of conversion, 1 - exp(-k t), bounded by mass conservation."""
    import math

    if strategy is None:
        strategy = DEFAULT_PROJECTION_STRATEGY
    drive = _projection_kinetic_drive(temperature_kelvin, time_minutes, strategy)
    return 1.0 - math.exp(-min(drive, 700.0))


def _projection_time_factor(time_minutes: Optional[float]) -> float:
    import math

    if time_minutes is None:
        return 1.0
    return 1.0 - math.exp(-max(time_minutes, 0.0) / 25.0)


def _estimate_projection_budget(
    corrected_initial: Dict[str, float],
    temperature_kelvin: float,
    time_minutes: Optional[float],
    strategy: ProjectionStrategy = DEFAULT_PROJECTION_STRATEGY,
) -> ProjectionBudget:
    if not corrected_initial:
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none", 0.0, 0.0)

    # Find limiting precursor (minimum positive molarity)
    positive_items = [(k, float(v)) for k, v in corrected_initial.items() if float(v) > 0.0]
    if not positive_items:
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none", 0.0, 0.0)

    limiting_name, limiting_val = min(positive_items, key=lambda x: x[1])
    limiting_precursor_molar = limiting_val * strategy.limiting_pool_to_molar_factor

    load_factor = _relative_precursor_load_factor(corrected_initial)
    # `severity` is retained as the bounded process-state index (melanoidin trapping,
    # depth bias, direct-sulfur bonus). It no longer scales the budget -- see the
    # docstring on _projection_temperature_factor.
    temperature_factor = _projection_temperature_factor(temperature_kelvin)
    time_factor = _projection_time_factor(time_minutes)
    severity = temperature_factor * time_factor

    # Budget scale: a first-order conversion extent under an apparent Arrhenius
    # drive. In the low-conversion regime the panel actually occupies (k*t ~ 1e-3)
    # this is essentially linear in k(T)*t, i.e. the budget inherits the true
    # exponential temperature dependence; it only rolls over once the limiting
    # precursor pool is genuinely being consumed (k*t ~ 1, above ~240 C at 60 min),
    # which is honest saturation rather than a parameterisation artifact.
    kinetic_drive = _projection_kinetic_drive(temperature_kelvin, time_minutes, strategy)
    conversion_extent = _projection_conversion_extent(temperature_kelvin, time_minutes, strategy)
    volatile_yield_fraction = (
        strategy.baseline_volatile_yield_fraction
        + (
            max(strategy.conversion_ceiling_fraction - strategy.baseline_volatile_yield_fraction, 0.0)
            * conversion_extent
        )
    )
    total_volatile_budget_molar = limiting_precursor_molar * volatile_yield_fraction * max(load_factor, 0.0)

    return ProjectionBudget(
        limiting_precursor_molar=float(limiting_precursor_molar),
        load_factor=float(max(load_factor, 0.0)),
        temperature_factor=float(temperature_factor),
        time_factor=float(time_factor),
        severity=float(severity),
        volatile_yield_fraction=float(volatile_yield_fraction),
        total_volatile_budget_molar=float(max(total_volatile_budget_molar, 0.0)),
        limiting_precursor_name=limiting_name,
        kinetic_drive=float(kinetic_drive),
        conversion_extent=float(conversion_extent),
    )


def _temporal_accessibility(total_tau_minutes: float, time_minutes: Optional[float]) -> float:
    import math

    if time_minutes is None:
        return 1.0
    if time_minutes <= 0.0:
        return 0.0
    if total_tau_minutes <= 0.0:
        return 1.0
    return 1.0 - math.exp(-time_minutes / total_tau_minutes)


def _relative_precursor_load_factor(corrected_initial: Dict[str, float]) -> float:
    import math

    positive_values = [max(float(value), 0.0) for value in corrected_initial.values() if float(value) > 0.0]
    if not positive_values:
        return 0.0

    limiting_value = min(positive_values)
    if limiting_value <= 0.0:
        return 0.0

    normalized = [max(value / limiting_value, 1.0e-12) for value in positive_values]
    return math.exp(sum(math.log(value) for value in normalized) / len(normalized))


def _projection_strategy_metadata(strategy: ProjectionStrategy) -> Dict[str, Any]:
    return {
        "projection_strategy_name": strategy.name,
        "projection_precursor_unit": strategy.precursor_concentration_unit,
        "projection_ppb_basis": strategy.ppb_basis,
        "projection_limiting_pool_to_molar_factor": float(strategy.limiting_pool_to_molar_factor),
        "projection_baseline_volatile_yield_fraction": float(strategy.baseline_volatile_yield_fraction),
        "projection_conversion_ceiling_fraction": float(strategy.conversion_ceiling_fraction),
        "projection_apparent_activation_energy_kj_mol": float(strategy.apparent_activation_energy_kj_mol),
        "projection_reference_temperature_kelvin": float(strategy.reference_temperature_kelvin),
        "projection_reference_conversion_time_min": float(strategy.reference_conversion_time_min),
        "projection_ppb_conversion_factor": float(strategy.ppb_conversion_factor),
        "projection_strategy_notes": strategy.notes,
    }


