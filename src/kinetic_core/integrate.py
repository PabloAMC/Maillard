"""
src/kinetic_core/integrate.py

STIFF-SAFE INTEGRATION OF THE KINETIC CORE.
===========================================

The network is stiff by construction: at 100 C the fastest constant
(1-deoxyglucosone consumption, 1.45 /min) is ~10^5 times the slowest
(fructose fragmentation, 4.4e-5 /min), and the pinned Amadori rearrangement is
another 45x above the condensation. The default method is therefore LSODA
(automatic stiff/non-stiff switching); BDF is offered for callers who want to
force the implicit branch.

NON-NEGATIVITY is guarded in two places, deliberately, because either alone is
insufficient:

  * inside the right-hand side, where reactant concentrations are clipped at
    zero before the rate law is evaluated, so a solver probe into slightly
    negative territory cannot produce a negative rate that then drives the
    state further negative (the runaway mode);
  * on the OUTPUT, where a residual negative of order the absolute tolerance
    is clipped for reporting.

Clipping the output alone would hide a runaway; clipping the RHS alone leaves
-1e-13 values in the returned array that then break a log-space objective.

TEMPERATURE AND TIME RANGE
--------------------------
Validated 100-200 C and seconds to hours (see
tests/unit/test_kinetic_core_b1_network.py). Note the difference between what
the integrator can DO and what the parameters LICENSE: every operative rate
constant was measured over 80-120 C, so integrating at 200 C is a numerically
sound extrapolation of an experimentally unsupported barrier. The run metadata
says so, on every run, via ``extrapolation_warnings``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np

from .network import REACTIONS, derivatives, rate_constants_at, reaction_rates
from .parameters import KineticParameter, registry_metadata
from .species import (
    INDEX,
    N_SPECIES,
    SPECIES_KEYS,
    initial_state,
    melanoidin_c_over_n,
    melanoidin_repeat_units,
    total_carbon,
    total_nitrogen,
)


@dataclass(frozen=True)
class CoreRun:
    """One integrated experiment, with its conservation audit attached."""

    times_min: np.ndarray
    #: shape (n_times, n_species), mmol/L (or mmol element/L for the pools)
    concentrations: np.ndarray
    temperature_k: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def series(self, species_key: str) -> np.ndarray:
        return self.concentrations[:, INDEX[species_key]]

    def carbon(self) -> np.ndarray:
        return np.array([total_carbon(row) for row in self.concentrations])

    def nitrogen(self) -> np.ndarray:
        return np.array([total_nitrogen(row) for row in self.concentrations])

    def melanoidin_mmol_L(self) -> np.ndarray:
        """Melanoidin in mmol of step-9 repeat units per litre."""
        return np.array(
            [melanoidin_repeat_units(row) for row in self.concentrations]
        )

    def melanoidin_c_over_n(self) -> np.ndarray:
        return np.array([melanoidin_c_over_n(row) for row in self.concentrations])

    def carbon_closure(self) -> Dict[str, float]:
        """Relative drift of the total-carbon invariant over the run."""
        carbon = self.carbon()
        c0 = float(carbon[0])
        drift = float(np.max(np.abs(carbon - c0)))
        return {
            "initial_carbon_mmol_L": c0,
            "max_abs_drift_mmol_L": drift,
            "max_relative_drift": drift / c0 if c0 > 0 else float("nan"),
        }

    def nitrogen_closure(self) -> Dict[str, float]:
        nitrogen = self.nitrogen()
        n0 = float(nitrogen[0])
        drift = float(np.max(np.abs(nitrogen - n0)))
        return {
            "initial_nitrogen_mmol_L": n0,
            "max_abs_drift_mmol_L": drift,
            "max_relative_drift": drift / n0 if n0 > 0 else float("nan"),
        }


def _extrapolation_warnings(
    parameters: Mapping[str, KineticParameter], temperature_c: float
) -> list:
    out = []
    for key, parameter in sorted(parameters.items()):
        low, high = parameter.temperature_range_c
        if temperature_c < low - 1e-9 or temperature_c > high + 1e-9:
            out.append(
                f"{key}: evaluated at {temperature_c:.1f} C, measured over "
                f"{low:.0f}-{high:.0f} C ({parameter.source_anchor})"
            )
    return out


def integrate(
    parameters: Mapping[str, KineticParameter],
    temperature_k: float,
    initial: Mapping[str, float],
    times_min: Sequence[float],
    *,
    method: str = "LSODA",
    rtol: float = 1e-9,
    atol: float = 1e-12,
    max_step: Optional[float] = None,
) -> CoreRun:
    """
    Integrate the kinetic core at one temperature.

    ``initial`` maps species keys to mmol/L; anything omitted starts at zero.
    ``times_min`` is the output grid in minutes, non-decreasing, starting at or
    after zero.
    """
    from scipy.integrate import solve_ivp

    if method not in ("LSODA", "BDF", "Radau"):
        raise ValueError(
            f"method {method!r} is not a stiff-capable solver; use LSODA, BDF or Radau"
        )
    grid = np.asarray(times_min, dtype=float)
    if grid.ndim != 1 or grid.size == 0:
        raise ValueError("times_min must be a non-empty 1-D sequence")
    if np.any(np.diff(grid) < 0):
        raise ValueError("times_min must be non-decreasing")
    if grid[0] < 0:
        raise ValueError("times_min must start at or after 0")

    k = rate_constants_at(parameters, float(temperature_k))
    y0 = initial_state(initial)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        # Reactant clipping lives in reaction_rates(); this keeps a solver probe
        # into negative territory from generating a negative rate.
        return derivatives(y, k)

    kwargs: Dict[str, Any] = {}
    if max_step is not None:
        kwargs["max_step"] = float(max_step)

    solution = solve_ivp(
        rhs,
        (0.0, float(max(grid[-1], 1e-12))),
        y0,
        t_eval=grid,
        method=method,
        rtol=rtol,
        atol=atol,
        **kwargs,
    )
    if not solution.success:
        raise RuntimeError(f"kinetic-core integration failed: {solution.message}")

    raw = solution.y.T
    # Guard 2: clip residual negatives of order atol on the OUTPUT. Anything
    # larger than that is a real failure and is raised rather than hidden.
    worst_negative = float(np.min(raw)) if raw.size else 0.0
    if worst_negative < -max(1e3 * atol, 1e-8):
        raise RuntimeError(
            f"kinetic-core integration produced a state of {worst_negative:.3e}, "
            f"far below the absolute tolerance {atol:.1e}: this is a genuine "
            f"non-negativity failure, not solver noise."
        )
    concentrations = np.clip(raw, 0.0, None)

    temperature_c = float(temperature_k) - 273.15
    metadata: Dict[str, Any] = {
        "method": method,
        "rtol": rtol,
        "atol": atol,
        "temperature_C": temperature_c,
        "temperature_K": float(temperature_k),
        "n_species": N_SPECIES,
        "n_reactions": len(REACTIONS),
        "species_order": list(SPECIES_KEYS),
        "worst_raw_negative_before_clip": worst_negative,
        "extrapolation_warnings": _extrapolation_warnings(parameters, temperature_c),
        "rate_constants_at_T": {key: float(value) for key, value in sorted(k.items())},
    }
    metadata.update(registry_metadata(parameters))

    run = CoreRun(
        times_min=grid,
        concentrations=concentrations,
        temperature_k=float(temperature_k),
        metadata=metadata,
    )
    metadata["carbon_closure"] = run.carbon_closure()
    metadata["nitrogen_closure"] = run.nitrogen_closure()
    return run


def flux_budget(
    parameters: Mapping[str, KineticParameter],
    temperature_k: float,
    initial: Mapping[str, float],
    duration_min: float,
    *,
    n_points: int = 401,
) -> Dict[str, float]:
    """
    Time-integrated flux through every reaction, mmol/L, over ``duration_min``.

    Used by the reports to say where the carbon actually went, rather than only
    what the end state looks like.
    """
    grid = np.linspace(0.0, float(duration_min), int(n_points))
    run = integrate(parameters, temperature_k, initial, grid)
    k = rate_constants_at(parameters, float(temperature_k))
    rates = np.array([reaction_rates(row, k) for row in run.concentrations])
    integral = np.trapezoid(rates, grid, axis=0)
    return {
        reaction.key: float(integral[j]) for j, reaction in enumerate(REACTIONS)
    }
