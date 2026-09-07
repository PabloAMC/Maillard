"""
src/kinetic_core/acrylamide_conditions.py -- WATER ACTIVITY ON THE ACRYLAMIDE LANE (wave B14, 2026-09-07).

WHY THIS EXISTS
===============
Until B14 the acrylamide lane carried NO water-activity term: a_w was metadata, and the engine
refused every comparison that moved it. Its own source measures the axis. De Vleeschouwer, Van
der Plancken, Van Loey & Hendrickx 2008 (JAFC 56:6460; `devleeschouwer2008_extraction.md`)
equilibrated equimolar asparagine-glucose powders to a_w 0.88 / 0.92 / 0.96 / 0.99 at 4 C, heated
them at 120-200 C and fitted the same multiresponse scheme at each a_w. The lane's shipped
constants ARE that paper's a_w 0.92 column (k_Fref 3.57e-3 +/- 1.38e-3 /min, Ea_F 159.2 +/- 29.5;
k_INT 1.70 +/- 1.05, Ea_INT 117.5 +/- 25.2 -- `parameters_acrylamide.py` cites them through the
2009 Part I paper, which reprints them). The authors' finding, at the 95 % HPD level: the
formation and elimination parameters "did not change significantly within the range of water
activities tested", nor with a potato matrix. Point estimates still move (Table 2, k_Fref at
160 C: 2.29 / 3.57 / 3.45 / 1.45 x 1e-3 /min), and that spread is what the band carries.

WHAT THIS MODULE DOES, AND DOES NOT DO
======================================
Inside the MEASURED WINDOW a_w 0.88-0.99 the lane carries a DECLARED FLAT term: the multiplier
on the acrylamide-forming step `k_int1_acr` is 1.0 (the source's finding), with an envelope
band (0.41, 1.39) = the union of the four point estimates relative to the shipped 0.92 column
(0.41-1.0) and that column's own 95 % HPD (0.61-1.39). A comparison that moves a_w within the
window is therefore ANSWERED, and answered flat within the band. Outside the window nothing is
measured -- the corpus's dry-side claims (extrusion at a_w 0.3-0.6) sit exactly there -- so the
engine keeps refusing a_w moves that leave the window, and a single run outside it is flagged.
At a_w None the term is inert (the run is at the constants' own a_w by assumption).
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

#: The measured window (De Vleeschouwer 2008 Table 1: KCl / Sr(NO3)2 / KNO3 / K2SO4 at 4 C).
AW_WINDOW: Tuple[float, float] = (0.88, 0.99)
#: The a_w the shipped constants were measured at (the Table 2 column the lane carries).
REFERENCE_AW = 0.92
#: The step the term scales: the acrylamide-forming step (k_Fref of the source).
SCALED_KEYS: Tuple[str, ...] = ("k_int1_acr",)
#: Table 2, k_Fref at 160 C (1e-3 /min) per a_w, and its 95 % HPD half-widths.
KF_TABLE: Dict[float, Tuple[float, float]] = {0.88: (2.29, 0.43), 0.92: (3.57, 1.38), 0.96: (3.45, 1.19), 0.99: (1.45, 0.42)}
#: Table 2, k_Eref at 160 C (/min) per a_w -- recorded, not applied (the lane's elimination
#: constant is a declared band, not this source's).
KE_TABLE: Dict[float, Tuple[float, float]] = {0.88: (0.11, 0.02), 0.92: (0.10, 0.04), 0.96: (0.09, 0.04), 0.99: (0.05, 0.03)}
#: The declared multiplier: FLAT (the source's 95 % HPD finding).
AW_MULTIPLIER = 1.0
#: The envelope band on the multiplier: point estimates relative to the 0.92 column (0.41-1.0)
#: united with the 0.92 column's own HPD (0.61-1.39).
AW_SCALE_BAND: Tuple[float, float] = (
    round(min(v[0] for v in KF_TABLE.values()) / KF_TABLE[REFERENCE_AW][0], 2),
    round((KF_TABLE[REFERENCE_AW][0] + KF_TABLE[REFERENCE_AW][1]) / KF_TABLE[REFERENCE_AW][0], 2),
)
AW_SOURCE = (
    "De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2008, J. Agric. Food Chem. 56:6460, "
    "Table 2 (equimolar Asn-Glc, a_w 0.88-0.99 at 4 C, 120-200 C, T_ref 160 C)"
)


def in_window(aw: Optional[float]) -> Optional[bool]:
    """True inside the measured window, False outside, None when no a_w was given."""
    if aw is None:
        return None
    return AW_WINDOW[0] - 1e-9 <= float(aw) <= AW_WINDOW[1] + 1e-9


def aw_multiplier(aw: Optional[float], scale: Optional[float] = None) -> float:
    """The multiplier on the acrylamide-forming step: 1.0 (flat) inside the window, drawn across
    the band by the envelope; exactly 1.0 at a_w None and outside the window (no term there)."""
    if in_window(aw) is not True:
        return 1.0
    return float(AW_MULTIPLIER if scale is None else scale)


def declarations(process) -> List[str]:
    """What a run must print about its water activity on this lane."""
    aw = getattr(process, "water_activity", None)
    state = in_window(aw)
    if state is None:
        return []
    if state:
        return [
            f"WATER ACTIVITY (B14): a_w {float(aw):.2f} is inside the measured window "
            f"{AW_WINDOW[0]:.2f}-{AW_WINDOW[1]:.2f}; the acrylamide lane carries a DECLARED FLAT "
            f"a_w term on {', '.join(SCALED_KEYS)} (multiplier {AW_MULTIPLIER:g}; envelope band "
            f"{AW_SCALE_BAND[0]:g}-{AW_SCALE_BAND[1]:g}, the source's own spread). Source: {AW_SOURCE}."
        ]
    return [
        f"WATER ACTIVITY (B14): a_w {float(aw):.2f} is OUTSIDE the measured window "
        f"{AW_WINDOW[0]:.2f}-{AW_WINDOW[1]:.2f}; the acrylamide lane's constants were measured at a_w "
        f"{REFERENCE_AW:.2f} and found flat only across that window ({AW_SOURCE}). No a_w term "
        "applies here: the value is recorded and changes no rate, and a comparison that moves a_w "
        "across this boundary is refused."
    ]


def apply(parameters: Mapping[str, Any], process, *, aw_scale: Optional[float] = None) -> Tuple[Dict[str, Any], List[str]]:
    """A COPY of the acrylamide parameter dict with the term applied, plus its declarations."""
    aw = getattr(process, "water_activity", None)
    factor = aw_multiplier(aw, aw_scale)
    out: Dict[str, Any] = dict(parameters)
    if factor != 1.0:
        for key in SCALED_KEYS:
            p = out.get(key)
            if p is not None and getattr(p, "k_ref", None) is not None:
                out[key] = replace(p, k_ref=float(p.k_ref) * factor)
    return out, declarations(process)


__all__ = ["AW_WINDOW", "REFERENCE_AW", "SCALED_KEYS", "KF_TABLE", "KE_TABLE", "AW_MULTIPLIER", "AW_SCALE_BAND",
           "AW_SOURCE", "in_window", "aw_multiplier", "declarations", "apply"]
