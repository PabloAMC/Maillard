"""
Observation operators: what an experiment measured, computed from the predicted state (sub-model 1, 2026-10-29).

A fit compares a measured row with the model's prediction OF THAT MEASUREMENT, not with the nearest state variable.
A pH read on a cooled sample is the charge balance at 25 C, not the in-situ pH at 120 C; a pH-stat's base uptake is
the strong base that holds the pH while acid is made; a titratable acidity is the base that takes a sample to an end
point; a diluted sample's pH is the charge balance of the diluted composition. Each operator is written once, here,
against the same charge balance the pH state uses (``ph_state``), and tested against a hand calculation
(tests/unit/test_observations.py). Nothing here is fitted.
"""
from __future__ import annotations

import dataclasses
from typing import Mapping

from .ph_state import (
    LABEL_TEMPERATURE_K,
    BufferSpec,
    PhDrift,
    ph_of_state,
    strong_ion_difference,
)

#: the end point of a titratable acidity (Brands 2001: titrated to pH 8.3, phenolphthalein's)
TITRATION_END_PH = 8.3


def ph_cooled(state: Mapping[str, float], sid_mmol_l: float, drift: PhDrift, buffer_spec: BufferSpec,
              temperature_k: float = LABEL_TEMPERATURE_K) -> float:
    """The pH of a sample cooled to ``temperature_k`` (25 C by default): the charge balance at the pot's conserved
    strong-ion difference. Ice-cooled samples read at room temperature are this, not the in-situ pH."""
    return float(ph_of_state(state, temperature_k, sid_mmol_l, drift, buffer_spec))


def diluted(buffer_spec: BufferSpec, factor: float) -> BufferSpec:
    """The buffer of a sample diluted ``factor``-fold in water: every molarity divided by the factor."""
    scaled = {f.name: getattr(buffer_spec, f.name) / factor
              for f in dataclasses.fields(buffer_spec) if f.name.endswith("_mol_l")}
    return dataclasses.replace(buffer_spec, **scaled)


def ph_diluted(state: Mapping[str, float], sid_mmol_l: float, drift: PhDrift, buffer_spec: BufferSpec,
               factor: float, temperature_k: float = LABEL_TEMPERATURE_K) -> float:
    """The pH of a sample diluted ``factor``-fold in water and read at ``temperature_k`` (Zhao 2017b: 0.5 mL to 3 mL,
    read at 25 C). Every solute, the buffer and the strong ions are diluted alike; the activity coefficients follow the
    lower ionic strength."""
    if factor < 1.0:
        raise ValueError(f"a dilution factor is at least 1, got {factor}")
    thinned = {k: float(v) / factor for k, v in state.items()}
    return float(ph_of_state(thinned, temperature_k, sid_mmol_l / factor, drift, diluted(buffer_spec, factor)))


def base_to_hold(target_ph: float, state: Mapping[str, float], reference_state: Mapping[str, float],
                 drift: PhDrift, buffer_spec: BufferSpec, temperature_k: float) -> float:
    """Strong base (mmol/L) that brings ``state`` to ``target_ph`` at ``temperature_k``, minus what the same end point
    takes for ``reference_state``. Two measurements are this quantity:

    * a pH-stat's NaOH uptake (Davidek 2002): the base added since t = 0 to hold the set pH at the reaction
      temperature -- ``target_ph`` the set pH, the reference the t = 0 state;
    * a titratable acidity relative to the unheated sample (Brands 2001): base to pH 8.3 on the cooled sample minus
      the same for the unheated one -- the pot's own strong ions cancel in the difference.
    """
    after = strong_ion_difference(target_ph, state, drift, buffer_spec, label_temperature_k=temperature_k)
    before = strong_ion_difference(target_ph, reference_state, drift, buffer_spec, label_temperature_k=temperature_k)
    return float(after - before)


def formic_plus_acetic(state: Mapping[str, float], drift: PhDrift) -> float:
    """Formic + acetic acid, mmol/L: the trunk's measured pools (FA, AA) plus the lumped organic-acid pool at the
    calibrated yield (``ACID``, the pentose sink's acid, which the charge balance counts at the same yield)."""
    return float(state.get("FA", 0.0) + state.get("AA", 0.0) + float(drift.acid_yield) * state.get("ACID", 0.0))
