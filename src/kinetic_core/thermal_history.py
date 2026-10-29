"""
src/kinetic_core/thermal_history.py -- THE POT'S TEMPERATURE AS ITS SOURCE RECORDS IT (wave B64, 2026-10-28).

WHY THIS EXISTS
===============
Every bundle carries a nominal ``temp_C`` and ``time_min``, and until B64 every pot was integrated as
ISOTHERMAL at ``temp_C`` from t = 0. The sources say otherwise. The sulfur lane's anchor laboratory runs
its reference pot "by raising the temperature within 20 min from 20 to 145 °C" in a 200 mL laboratory
autoclave (Hofmann & Schieberle 1995, 1997, verbatim), and the 1998 paper that supplies the fit's
145 °C rows writes the same pot as "20 min at 145 °C"; van Seeventer 2001 heats "from room temperature
to 130 °C in 10 min" and then holds 20 min. A 20 -> 145 °C ramp exposes a pot to about 2-4 min at
145 °C for barriers of 64-100 kJ/mol, not 20.

WHAT IT DOES
============
It reads a bundle's ``conditions.heating`` block and returns the ``ThermalProgram`` it describes:
a piecewise-linear profile (``ThermalProgram.from_profile``) from the stated start temperature, through
the stated ramp, hold and cool-down. A bundle with no block, or whose block says ``isothermal``, gets
``ThermalProgram.isothermal(temp_C, time_min)``, exactly what it had before. Nothing here is fitted:
every number is the source's, and a missing one is refused rather than guessed (``HeatingUnstated``).

THE BLOCK
=========
``time_convention`` says what the source's stated time MEANS:

``isothermal``           the source states, or the bundle declares, a hold from t = 0 (the pre-B64 reading)
``ramp_only``            the stated time IS the heat-up: start_C -> temp_C over ``ramp_min`` (= time_min)
``ramp_plus_hold``       a heat-up of ``ramp_min``, then ``hold_min`` at temp_C (both stated)
``hold_after_reaching``  ``time_min`` counts from when the pot reached temp_C, after a ``ramp_min`` come-up
``includes_come_up``     ``time_min`` counts from insertion; the first ``ramp_min`` of it is the come-up

and optionally ``cool_min`` / ``cool_to_C``, a linear cool-down after the hold (a quench is omitted).
``start_C`` defaults to 20 °C, the temperature every source here starts its ramps at.
"""
from __future__ import annotations

from typing import Any, Mapping, Optional

TIME_CONVENTIONS = (
    "isothermal",
    "ramp_only",
    "ramp_plus_hold",
    "hold_after_reaching",
    "includes_come_up",
)
#: the ramp start when a source states a ramp but not where it starts ("room temperature")
DEFAULT_START_C = 20.0


class HeatingUnstated(ValueError):
    """A heating block names a convention whose numbers it does not carry."""


def _number(block: Mapping[str, Any], key: str, convention: str) -> float:
    value = block.get(key)
    if value is None:
        raise HeatingUnstated(f"heating convention {convention!r} needs {key!r}, which the block does not state")
    return float(value)


def profile_points(heating: Optional[Mapping[str, Any]], temp_c: float, time_min: float):
    """The breakpoints ``((t_min, T_C), ...)`` of the recorded history, or None for an isothermal pot."""
    if not heating:
        return None
    convention = str(heating.get("time_convention", "isothermal"))
    if convention not in TIME_CONVENTIONS:
        raise ValueError(f"unknown heating time_convention {convention!r}; one of {TIME_CONVENTIONS}")
    if convention == "isothermal":
        return None
    start = float(heating.get("start_C", DEFAULT_START_C))
    if convention == "ramp_only":
        ramp, hold = float(heating.get("ramp_min", time_min)), 0.0
    elif convention == "ramp_plus_hold":
        ramp, hold = _number(heating, "ramp_min", convention), _number(heating, "hold_min", convention)
    elif convention == "hold_after_reaching":
        ramp, hold = _number(heating, "ramp_min", convention), float(time_min)
    else:  # includes_come_up
        ramp = _number(heating, "ramp_min", convention)
        if ramp > time_min:
            raise ValueError(f"a {ramp:g} min come-up inside a {time_min:g} min run")
        hold = float(time_min) - ramp
    points = [(0.0, start), (ramp, float(temp_c))]
    if hold > 0.0:
        points.append((ramp + hold, float(temp_c)))
    cool = heating.get("cool_min")
    if cool:
        points.append((points[-1][0] + float(cool), _number(heating, "cool_to_C", convention)))
    return tuple(points)


def thermal_program(heating: Optional[Mapping[str, Any]], temp_c: float, time_min: float):
    """The ``ThermalProgram`` a bundle's heating block describes (isothermal when it describes none)."""
    from .engine import ThermalProgram

    points = profile_points(heating, temp_c, time_min)
    if points is None:
        return ThermalProgram.isothermal(float(temp_c), float(time_min))
    return ThermalProgram.from_profile(points)
