"""
src/trunk_kinetics.py

A MINIMAL, DEDICATED ODE INTEGRATOR FOR THE SUGAR-TRUNK SUBSYSTEM.
=================================================================

WHY THIS MODULE EXISTS (Wave S3, 2026-08-27)
--------------------------------------------
This repository has two kinetics lanes, and *neither* of them can be honestly
fitted to concentration-vs-time data:

  * The **FAST screening lane** (`src/barrier_constants.FAST_BARRIERS` ->
    `src/recommend.py` -> `src/benchmark_validation.py`) is not an ODE
    integrator at all. It propagates relative fluxes through an enumerated
    network; there is no time axis and no absolute rate constant in it.

  * The **Cantera export lane** (`data/lit/arrhenius_params.yml` ->
    `src/cantera_export.py` -> `scripts/run_cantera_kinetics.py`) *does*
    integrate ODEs, but it cannot represent the experiment. Measured at
    HEAD b2a1b20 on the glucose/glycine system, four blockers, each verified
    by running the lane rather than by reading it:

      1. **THE INITIAL MOLARITY IS DISCARDED.** `src/kinetics.py:373` does
         `phase.X = initial_concentrations`; Cantera normalises `X` to sum to
         one, so only the *ratio* of the precursors survives and the absolute
         concentration is then set by the phase molar density implied by
         Girolami molar-volume estimates. Measured: feeding 0.02 M, 0.2 M and
         2.0 M glucose/glycine gives BIT-IDENTICAL trajectories (Amadori at
         9000 s = 1.93585 kmol/m3 in all three cases), against a t=0 glucose
         concentration of 2.546 mol/L that nobody asked for. A bimolecular
         rate constant cannot be fitted through a lane that does not know the
         concentration.
      2. **FOUR OF THE TEN MEASURED SPECIES DO NOT EXIST IN THE NETWORK.**
         The SMIRKS enumeration for glucose + glycine emits no fructose (no
         Lobry de Bruyn-Alberda van Ekenstein isomerisation step), no formic
         acid, no acetic acid and no melanoidins.
      3. **THE FAMILIES ARE LUMPED IN THE YAML.** `_arrhenius_yaml_key` maps
         `1,2-enolisation`, `2,3-enolisation` and `enolisation_intermediate`
         onto a single `enolisation` A/Ea pair, so three physically distinct
         and separately-measured steps cannot take three different rates.
      4. **THE YAML BARRIER IS NOT EVEN THE OPERATIVE NUMBER.**
         `src/cantera_export.py:161` does
         `barrier_kcal = max(barrier_kcal, float(barrier_lit))` -- a max
         against the FAST barrier, not an assignment -- and every step is
         then marked reversible with its equilibrium set by Joback-estimated
         NASA polynomials. A fitted Ea would be entangled with all of that.

Distorting either shipped lane to make it fittable would have been worse than
the disease. So the fitted subsystem lives here, in the smallest, most
inspectable form that can reproduce the measurement: explicit mass-action
ODEs, `scipy.solve_ivp`, no thermo estimation, no gating, no unit games.

WHAT THIS MODULE IS NOT
-----------------------
**Nothing in the shipped prediction path imports this module.** It is the
calibration lane and only the calibration lane. It does not change any
prediction, any benchmark score or any FAST barrier. See
`results/validation/trunk_rate_calibration.md` for what did and did not
propagate, and the "two-lane split" note in
`src/barrier_constants.py` / `data/lit/arrhenius_params.yml`.

TEMPERATURE PARAMETERISATION
----------------------------
Rate constants are carried as (k_ref, Ea) at a reference temperature rather
than as (A, Ea). This is Martins' own reparameterisation (thesis eqs 5.3-5.6,
`10.1016/j.foodchem.2004.04.006`): over an 80-120 C window the Arrhenius A and
Ea are almost perfectly correlated, and fitting A directly produces a fit whose
parameters are individually meaningless. With

    k(T) = k_ref * exp( -(Ea/R) * (1/T - 1/T_ref) )

k_ref and Ea are close to orthogonal, which is what makes the per-parameter
confidence intervals below worth printing. T_ref is 373.15 K (100 C), the
midpoint of the 80/100/120 C data and the same T_av Martins used.

UNITS ARE NATIVE THROUGHOUT
---------------------------
Concentrations mmol/L, time minutes, temperature Kelvin, Ea kJ/mol.
Unimolecular k in 1/min; the one bimolecular k in L/(mmol*min). These are the
units the source data are printed in; nothing is converted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np

# Gas constant in kJ/(mol K) -- matches the kJ/mol unit that Ea is carried in.
R_KJ = 8.314462618e-3

# Reference temperature for the reparameterised Arrhenius form: 100 C.
# This is the midpoint of the 80/100/120 C fit corpus and equals the T_av of
# Martins' own five-temperature reparameterisation.
T_REF_K = 373.15


# ---------------------------------------------------------------------------
# The trunk
# ---------------------------------------------------------------------------
# Eight mass-action steps. The `family` column records which repository
# reaction family each step corresponds to, and -- just as important -- which
# steps correspond to NO repository family. Those are declared as structural
# gaps rather than quietly folded into a neighbouring barrier.
#
#   step key      transformation                     repo family
#   ------------  ---------------------------------  --------------------------
#   k_schiff      Glc + Gly -> SB                     schiff_condensation
#   k_amadori     SB -> DFG                           amadori_rearrangement
#   k_dfg_3dg     DFG -> 3DG + Gly                    enolisation_intermediate
#   k_dfg_1dg     DFG -> 1DG + Gly                    2,3-enolisation
#   k_dfg_other   DFG -> X (no Gly released)          NONE (structural gap)
#   k_3dg_out     3DG -> products                     1,2-enolisation (lump)
#   k_1dg_out     1DG -> products                     furanone_amino_acid_reduction (lump)
#   k_glc_other   Glc -> Fru / acids                  NONE (structural gap)
#
# The two "NONE" rows are the honest part. `k_glc_other` lumps the Lobry de
# Bruyn-Alberda van Ekenstein isomerisation and the amine-independent sugar
# fragmentation to formic/acetic acid: at 100 C those account for more than
# half of the measured glucose loss (glucose falls 70 mmol/L in 150 min while
# glycine falls only 29), and the repository's network contains neither. If it
# were omitted the condensation step would be forced to absorb all of the
# glucose loss and every downstream constant would inherit the error.
# `k_dfg_other` is the DFG channel that does NOT regenerate glycine; it is
# required by two independent observations (see the fit report).

STEP_KEYS: Tuple[str, ...] = (
    "k_schiff",
    "k_amadori",
    "k_dfg_3dg",
    "k_dfg_1dg",
    "k_dfg_other",
    "k_3dg_out",
    "k_1dg_out",
    "k_glc_other",
)

STEP_FAMILY: Mapping[str, str | None] = {
    "k_schiff": "schiff_condensation",
    "k_amadori": "amadori_rearrangement",
    "k_dfg_3dg": "enolisation_intermediate",
    "k_dfg_1dg": "2,3-enolisation",
    "k_dfg_other": None,
    "k_3dg_out": "1,2-enolisation",
    "k_1dg_out": "furanone_amino_acid_reduction",
    "k_glc_other": None,
}

STEP_UNIT: Mapping[str, str] = {
    "k_schiff": "L/(mmol*min)",
    "k_amadori": "1/min",
    "k_dfg_3dg": "1/min",
    "k_dfg_1dg": "1/min",
    "k_dfg_other": "1/min",
    "k_3dg_out": "1/min",
    "k_1dg_out": "1/min",
    "k_glc_other": "1/min",
}

# State vector order. `SB` (the Schiff base / E1 intermediate) is NOT measured
# by any experiment in the fit corpus and is carried only so that the
# condensation can be split into two steps at all -- which is the whole point
# of the Schiff/Amadori authority-inversion test.
SPECIES: Tuple[str, ...] = ("Glc", "Gly", "SB", "DFG", "TDG", "ODG")
_IDX = {name: i for i, name in enumerate(SPECIES)}


def arrhenius_k(k_ref: float, ea_kj_mol: float, temperature_k: float) -> float:
    """
    Reparameterised Arrhenius: k(T) = k_ref * exp(-(Ea/R)(1/T - 1/T_ref)).

    `k_ref` is the rate constant at `T_REF_K` (373.15 K), in whatever unit the
    step carries. Returns the same unit.
    """
    return float(k_ref) * float(
        np.exp(-(ea_kj_mol / R_KJ) * (1.0 / float(temperature_k) - 1.0 / T_REF_K))
    )


def rate_constants_at(
    params: Mapping[str, Tuple[float, float]], temperature_k: float
) -> Dict[str, float]:
    """Evaluate every step's rate constant at `temperature_k`.

    `params` maps a step key to `(k_ref, Ea_kJ_per_mol)`.
    """
    return {
        key: arrhenius_k(params[key][0], params[key][1], temperature_k)
        for key in STEP_KEYS
    }


def derivatives(state: np.ndarray, k: Mapping[str, float]) -> np.ndarray:
    """Mass-action right-hand side for the eight-step trunk.

    Written out explicitly rather than assembled from a stoichiometric matrix:
    at this size the explicit form is the documentation.
    """
    glc, gly, sb, dfg, tdg, odg = state

    r_schiff = k["k_schiff"] * glc * gly          # Glc + Gly -> SB
    r_amadori = k["k_amadori"] * sb               # SB -> DFG
    r_3dg = k["k_dfg_3dg"] * dfg                  # DFG -> 3DG + Gly
    r_1dg = k["k_dfg_1dg"] * dfg                  # DFG -> 1DG + Gly
    r_other = k["k_dfg_other"] * dfg              # DFG -> X (no Gly)
    r_3dg_out = k["k_3dg_out"] * tdg              # 3DG -> products
    r_1dg_out = k["k_1dg_out"] * odg              # 1DG -> products
    r_glc_other = k["k_glc_other"] * glc          # Glc -> Fru / acids

    return np.array(
        [
            -r_schiff - r_glc_other,                       # d[Glc]
            -r_schiff + r_3dg + r_1dg,                     # d[Gly]
            r_schiff - r_amadori,                          # d[SB]
            r_amadori - r_3dg - r_1dg - r_other,           # d[DFG]
            r_3dg - r_3dg_out,                             # d[3DG]
            r_1dg - r_1dg_out,                             # d[1DG]
        ],
        dtype=float,
    )


@dataclass(frozen=True)
class TrunkRun:
    """One integrated experiment."""

    times_min: np.ndarray
    #: shape (len(times_min), len(SPECIES))
    concentrations: np.ndarray
    #: cumulative glycine RELEASED by DFG degradation (the DFG-degradation
    #: experiment measures this directly, and it is not the same quantity as
    #: the `Gly` state when glycine is also being consumed).
    glycine_released: np.ndarray

    def series(self, species: str) -> np.ndarray:
        return self.concentrations[:, _IDX[species]]


def integrate(
    params: Mapping[str, Tuple[float, float]],
    temperature_k: float,
    initial: Mapping[str, float],
    times_min: Sequence[float],
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> TrunkRun:
    """
    Integrate the trunk at one temperature.

    `initial` maps species names in `SPECIES` to mmol/L; anything omitted
    starts at zero. `times_min` is the output grid in minutes and must be
    non-decreasing and start at or after 0.

    The integrator is stiff-capable (`LSODA`) because the fitted Amadori step
    is fast relative to the condensation, which is exactly the stiffness the
    Schiff/Amadori split question is about.
    """
    from scipy.integrate import solve_ivp

    k = rate_constants_at(params, temperature_k)

    y0 = np.zeros(len(SPECIES), dtype=float)
    for name, value in initial.items():
        if name not in _IDX:
            raise KeyError(f"unknown species {name!r}; expected one of {SPECIES}")
        y0[_IDX[name]] = float(value)

    grid = np.asarray(times_min, dtype=float)
    if grid.ndim != 1 or grid.size == 0:
        raise ValueError("times_min must be a non-empty 1-D sequence")
    if np.any(np.diff(grid) < 0):
        raise ValueError("times_min must be non-decreasing")
    if grid[0] < 0:
        raise ValueError("times_min must start at or after 0")

    # Augment the state with the cumulative glycine-release integral so that
    # it is produced by the same solver, on the same error control, rather
    # than by post-hoc quadrature on the output grid.
    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        d = derivatives(y[: len(SPECIES)], k)
        released = (k["k_dfg_3dg"] + k["k_dfg_1dg"]) * y[_IDX["DFG"]]
        return np.concatenate([d, [released]])

    t_span = (0.0, float(max(grid[-1], 1e-12)))
    sol = solve_ivp(
        rhs,
        t_span,
        np.concatenate([y0, [0.0]]),
        t_eval=grid,
        method="LSODA",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"trunk integration failed: {sol.message}")

    y = sol.y.T
    return TrunkRun(
        times_min=grid,
        concentrations=np.clip(y[:, : len(SPECIES)], 0.0, None),
        glycine_released=np.clip(y[:, len(SPECIES)], 0.0, None),
    )


def carbon_like_balance(run: TrunkRun) -> np.ndarray:
    """
    Sum of the sugar-derived pools that the trunk conserves.

    Every step except `k_3dg_out`, `k_1dg_out`, `k_dfg_other` and
    `k_glc_other` maps one sugar-derived molecule onto one sugar-derived
    molecule, so `Glc + SB + DFG + 3DG + 1DG` is non-increasing and decreases
    only through those four declared sinks. Used by the unit tests to catch a
    stoichiometry typo in `derivatives`.
    """
    idx = [_IDX[n] for n in ("Glc", "SB", "DFG", "TDG", "ODG")]
    return run.concentrations[:, idx].sum(axis=1)
