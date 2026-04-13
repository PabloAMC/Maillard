from dataclasses import dataclass
from typing import Dict, Any, Optional

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


@dataclass(frozen=True)
class ProjectionStrategy:
    name: str
    precursor_concentration_unit: str
    limiting_pool_to_molar_factor: float
    baseline_volatile_yield_fraction: float
    severity_volatile_yield_slope: float
    ppb_conversion_factor: float
    ppb_basis: str
    notes: str


DEFAULT_PROJECTION_STRATEGY = ProjectionStrategy(
    name="precursor_limited_observable_v1",
    precursor_concentration_unit="mM",
    limiting_pool_to_molar_factor=1.0e-3,
    baseline_volatile_yield_fraction=1.0e-6,
    severity_volatile_yield_slope=1.5e-3,
    ppb_conversion_factor=1.0e6,
    ppb_basis="aqueous_mass_equivalent_ppb",
    notes=(
        "Allocates a conservative volatile budget from the limiting precursor pool, then converts "
        "M to ppb via MW assuming dilute aqueous density (~1 kg/L) before matrix/headspace projection."
    ),
)

def _thermal_severity(temperature_kelvin: float, time_minutes: Optional[float]) -> float:
    return _projection_temperature_factor(temperature_kelvin) * _projection_time_factor(time_minutes)


def _projection_temperature_factor(temperature_kelvin: float) -> float:
    import math

    temp_c = temperature_kelvin - 273.15
    return 1.0 / (1.0 + math.exp(-(temp_c - 110.0) / 18.0))


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
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none")

    # Find limiting precursor (minimum positive molarity)
    positive_items = [(k, float(v)) for k, v in corrected_initial.items() if float(v) > 0.0]
    if not positive_items:
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none")

    limiting_name, limiting_val = min(positive_items, key=lambda x: x[1])
    limiting_precursor_molar = limiting_val * strategy.limiting_pool_to_molar_factor
    
    load_factor = _relative_precursor_load_factor(corrected_initial)
    temperature_factor = _projection_temperature_factor(temperature_kelvin)
    time_factor = _projection_time_factor(time_minutes)
    severity = temperature_factor * time_factor
    volatile_yield_fraction = (
        strategy.baseline_volatile_yield_fraction
        + strategy.severity_volatile_yield_slope * severity
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
        limiting_precursor_name=limiting_name
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
        "projection_severity_volatile_yield_slope": float(strategy.severity_volatile_yield_slope),
        "projection_ppb_conversion_factor": float(strategy.ppb_conversion_factor),
        "projection_strategy_notes": strategy.notes,
    }


