"""
src/kinetic_core/acrylamide_conditions.py -- WATER ACTIVITY AND pH ON THE ACRYLAMIDE LANE
(wave B14, 2026-09-07: a flat a_w term inside 0.88-0.99; wave B15, 2026-09-07: the window extended
to 0.34, an elimination a_w shape, and a pH factor -- all DECLARED from De Vleeschouwer's own series).

WHY THIS EXISTS
===============
Until B14 the acrylamide lane carried NO water-activity term and NO pH term, and the engine refused
every comparison that moved either axis. The lane's own laboratory measured both axes on the same
equimolar asparagine-glucose system with the same multiresponse fit:

  * a_w, high side -- De Vleeschouwer 2008 (JAFC 56:6460; `devleeschouwer2008_extraction.md`):
    powders at a_w 0.88 / 0.92 / 0.96 / 0.99, 120-200 C. The lane's shipped constants ARE that
    paper's 0.92 column. Formation and elimination "did not change significantly" (95 % HPD).
  * a_w, dry side -- De Vleeschouwer 2007 (Biotechnol. Prog. 23:722; `devleeschouwer2007_extraction.md`):
    powders at a_w 0.34 / 0.59 / 0.73 / 0.82 / 0.88 / 0.92. Formation k_F "varies only slightly"
    (0.71-1.09 of the 0.92 value); elimination k_E has a MINIMUM at a_w 0.82 (0.33 of the 0.92 value)
    and the Maillard competition a maximum there.
  * pH -- De Vleeschouwer 2006 (JAFC 54:7847; `devleeschouwer2006_extraction.md`): 0.1 M Asn + Glc in
    0.05 M phosphate at initial pH 4 / 6 / 8, 120-200 C. ln k is linear in pH: formation slope
    0.5414 +/- 0.106, elimination 0.3442 +/- 0.013 per pH unit in NATURAL-log units (the paper's
    "log-linear" fit; the two-point check ln(37.5/4.30)/4 = 0.54 confirms the base), i.e. 0.235 and
    0.149 decades per pH unit (potato matrix: 0.187 / 0.148).

WHAT THIS MODULE DOES, AND DOES NOT DO
======================================
`apply(parameters, process, ...)` returns a COPY of the acrylamide parameter dict with three declared
factors on named steps, plus the declarations the run must print:

  * FORMATION a_w multiplier on `k_int1_acr`: exactly 1.0 (flat) inside the measured window
    0.34-0.99; envelope band (0.41, 1.39) = the two papers' point-estimate spreads and the 0.92
    column's HPD. No term outside the window (recorded, refused across the boundary).
  * ELIMINATION a_w multiplier on `k_acr_dp`: piecewise-linear through the 2007 k_E ratios relative
    to a_w 0.92 for the dry side (0.76 / 0.60 / 0.33 / 0.37 at 0.34 / 0.59 / 0.73 / 0.82), joining
    1.0 at a_w 0.88 -- the floor of the 2008 window where the same laboratory measured the constants
    FLAT (2007's own 0.88 point, 0.66 +/- 0.33, is within one SE of 1) -- and 1.0 at and above 0.88;
    the envelope scales the deficit (1 - multiplier) by a factor in (0, 1.2), so the band reaches
    from "no effect" to 1.2x the shape.
  * pH factors 10^(s (pH - 6.8)) on `k_asn_glc` (s = 0.235 decades per unit, band 0.11-0.28) and on
    `k_acr_dp` (s = 0.149, band 0.12-0.16) inside the measured window pH 4-8, HELD at the window edge outside it
    with a warning. The reference is the lane's declared network pH (6.8), where the factor is 1.

At a_w None and pH 6.8 every factor is exactly 1.0, so every earlier wave reproduces bit for bit.
All factors are within-study ratios installed as constants with bands -- the standing of a measured
barrier override, never a fitted coordinate. The 2006 constants are indexed to INITIAL pH while the
pot drifts by up to two units during heating; the lane has no pH trajectory, so this is an
initial-pH factor by construction and says so.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Water activity
# ---------------------------------------------------------------------------

#: The measured window: De Vleeschouwer 2007 (0.34-0.92) and 2008 (0.88-0.99), joined.
AW_WINDOW: Tuple[float, float] = (0.34, 0.99)
#: The a_w the shipped constants were measured at (the 2008 Table 2 column the lane carries).
REFERENCE_AW = 0.92
#: The step the FORMATION term scales.
SCALED_KEYS: Tuple[str, ...] = ("k_int1_acr",)
#: The step the ELIMINATION term scales (the lane's first-order acrylamide degradation).
ELIMINATION_KEYS: Tuple[str, ...] = ("k_acr_dp",)
#: 2008 Table 2, k_Fref at 160 C (1e-3 /min) per a_w, with 95 % HPD half-widths.
KF_TABLE: Dict[float, Tuple[float, float]] = {0.88: (2.29, 0.43), 0.92: (3.57, 1.38), 0.96: (3.45, 1.19), 0.99: (1.45, 0.42)}
#: 2008 Table 2, k_Eref at 160 C (/min) -- recorded, not applied (the 2007 series carries the shape).
KE_TABLE: Dict[float, Tuple[float, float]] = {0.88: (0.11, 0.02), 0.92: (0.10, 0.04), 0.96: (0.09, 0.04), 0.99: (0.05, 0.03)}
#: 2007 Table 2, k_Fref (1e-3 /M/min) and k_Eref (1e-3 /min) per a_w, with SEs.
KF_2007: Dict[float, Tuple[float, float]] = {0.34: (0.647, 0.168), 0.59: (0.736, 0.110), 0.73: (0.675, 0.138),
                                             0.82: (0.996, 0.259), 0.88: (0.703, 0.058), 0.92: (0.913, 0.082)}
KE_2007: Dict[float, Tuple[float, float]] = {0.34: (371.0, 333.0), 0.59: (291.0, 96.4), 0.73: (162.0, 36.5),
                                             0.82: (178.0, 41.4), 0.88: (321.0, 161.0), 0.92: (485.0, 186.0)}
#: The declared FORMATION multiplier: FLAT (both papers' finding).
AW_MULTIPLIER = 1.0
#: The envelope band on the formation multiplier: 2008 point estimates relative to the 0.92 column
#: (0.41-1.0) united with that column's HPD (0.61-1.39); the 2007 spread (0.71-1.09) lies inside it.
AW_SCALE_BAND: Tuple[float, float] = (
    round(min(v[0] for v in KF_TABLE.values()) / KF_TABLE[REFERENCE_AW][0], 2),
    round((KF_TABLE[REFERENCE_AW][0] + KF_TABLE[REFERENCE_AW][1]) / KF_TABLE[REFERENCE_AW][0], 2),
)
#: The ELIMINATION shape: 2007 k_E relative to a_w 0.92 on the dry side, joining the 2008 flat
#: window (1.0) at a_w 0.88. 2007's own 0.88 point (0.66 +/- 0.33) is within one SE of 1.
AW_ELIMINATION_TABLE: Tuple[Tuple[float, float], ...] = tuple(
    (aw, round(KE_2007[aw][0] / KE_2007[0.92][0], 2)) for aw in sorted(KE_2007) if aw <= 0.82
) + ((0.88, 1.0),)
#: The envelope scales the elimination DEFICIT (1 - multiplier): 0 = no effect, 1.2 = 1.2x the shape.
AW_ELIMINATION_SCALE_BAND: Tuple[float, float] = (0.0, 1.2)
AW_SOURCE = (
    "De Vleeschouwer et al. 2008 JAFC 56:6460 Table 2 (a_w 0.88-0.99) and 2007 Biotechnol. Prog. 23:722 "
    "Table 2 (a_w 0.34-0.92); equimolar Asn-Glc powders, 120-200 C, T_ref 160 C"
)

# ---------------------------------------------------------------------------
# pH
# ---------------------------------------------------------------------------

#: The lane's declared network pH: every factor is exactly 1 there.
REFERENCE_PH = 6.8
#: De Vleeschouwer 2006 Table 1, phosphate: the measured window of initial pH.
PH_WINDOW: Tuple[float, float] = (4.0, 8.0)
#: The paper's ln-linear slopes per pH unit (phosphate; potato matrix 0.4312 +/- 0.168 / 0.3397 +/- 0.073).
LN_SLOPE_FORMATION, LN_SLOPE_FORMATION_SE = 0.5414, 0.106
LN_SLOPE_ELIMINATION, LN_SLOPE_ELIMINATION_SE = 0.3442, 0.013
LN_SLOPE_FORMATION_POTATO, LN_SLOPE_FORMATION_POTATO_SE = 0.4312, 0.168
LN_SLOPE_ELIMINATION_POTATO, LN_SLOPE_ELIMINATION_POTATO_SE = 0.3397, 0.073
_LN10 = 2.302585092994046
#: ... in decades per pH unit, the exponent the factor uses.
PH_EXPONENT_FORMATION = round(LN_SLOPE_FORMATION / _LN10, 3)          # 0.235
PH_EXPONENT_ELIMINATION = round(LN_SLOPE_ELIMINATION / _LN10, 3)      # 0.149
#: Bands: the potato slope minus its SE to the phosphate slope plus its SE, in decades.
PH_EXPONENT_FORMATION_BAND: Tuple[float, float] = (
    round((LN_SLOPE_FORMATION_POTATO - LN_SLOPE_FORMATION_POTATO_SE) / _LN10, 3),
    round((LN_SLOPE_FORMATION + LN_SLOPE_FORMATION_SE) / _LN10, 3),
)
PH_EXPONENT_ELIMINATION_BAND: Tuple[float, float] = (
    round((LN_SLOPE_ELIMINATION_POTATO - LN_SLOPE_ELIMINATION_POTATO_SE) / _LN10, 3),
    round((LN_SLOPE_ELIMINATION + LN_SLOPE_ELIMINATION_SE) / _LN10, 3),
)
#: The steps the pH factors scale: initiation (formation) and the first-order elimination.
PH_FORMATION_KEYS: Tuple[str, ...] = ("k_asn_glc",)
PH_ELIMINATION_KEYS: Tuple[str, ...] = ("k_acr_dp",)
#: 2006 Table 1, phosphate: k_Fref (1e-3 /M/min), k_Eref (1e-3 /min) per initial pH.
KF_PH_2006: Dict[float, Tuple[float, float]] = {8.0: (37.5, 4.21), 6.0: (8.78, 1.34), 4.0: (4.30, 0.636)}
KE_PH_2006: Dict[float, Tuple[float, float]] = {8.0: (333.6, 41.4), 6.0: (175.3, 36.1), 4.0: (84.2, 18.8)}
PH_SOURCE = (
    "De Vleeschouwer et al. 2006 JAFC 54:7847 Table 1 (0.1 M Asn + Glc, 0.05 M phosphate, initial pH "
    "4 / 6 / 8, 120-200 C): ln k linear in pH, slopes 0.5414 +/- 0.106 (formation) and 0.3442 +/- 0.013 "
    "(elimination) per pH unit = 0.235 / 0.149 decades"
)


def in_window(aw: Optional[float]) -> Optional[bool]:
    """True inside the measured a_w window, False outside, None when no a_w was given."""
    if aw is None:
        return None
    return AW_WINDOW[0] - 1e-9 <= float(aw) <= AW_WINDOW[1] + 1e-9


def aw_multiplier(aw: Optional[float], scale: Optional[float] = None) -> float:
    """The FORMATION multiplier: 1.0 (flat) inside the window, drawn across the band by the envelope;
    exactly 1.0 at a_w None and outside the window (no term there)."""
    if in_window(aw) is not True:
        return 1.0
    return float(AW_MULTIPLIER if scale is None else scale)


def aw_elimination_multiplier(aw: Optional[float], deficit_scale: float = 1.0) -> float:
    """The ELIMINATION multiplier: the 2007 shape relative to a_w 0.92, held below the table, 1.0 at
    and above 0.92 and at a_w None; the deficit (1 - m) scaled by the envelope's draw."""
    if aw is None or in_window(aw) is not True:
        return 1.0
    a = float(aw)
    table = AW_ELIMINATION_TABLE
    if a <= table[0][0]:
        m = table[0][1]          # the window floor IS the first table point
    elif a >= table[-1][0]:
        m = 1.0
    else:
        m = 1.0
        for (a0, m0), (a1, m1) in zip(table, table[1:]):
            if a0 <= a <= a1:
                m = m0 + (m1 - m0) * (a - a0) / (a1 - a0)
                break
    return 1.0 - float(deficit_scale) * (1.0 - m)


def ph_in_window(ph: Optional[float]) -> Optional[bool]:
    if ph is None:
        return None
    return PH_WINDOW[0] - 1e-9 <= float(ph) <= PH_WINDOW[1] + 1e-9


def ph_factor(ph: Optional[float], exponent: float) -> float:
    """10^(exponent (pH - 6.8)), with the pH HELD at the window edge outside the measured window."""
    if ph is None:
        return 1.0
    p = min(max(float(ph), PH_WINDOW[0]), PH_WINDOW[1])
    return 10.0 ** (float(exponent) * (p - REFERENCE_PH))


def declarations(process) -> List[str]:
    """What a run must print about its water activity and pH on this lane."""
    out: List[str] = []
    aw = getattr(process, "water_activity", None)
    state = in_window(aw)
    if state is True:
        out.append(
            f"WATER ACTIVITY (B14/B15): a_w {float(aw):.2f} is inside the measured window "
            f"{AW_WINDOW[0]:.2f}-{AW_WINDOW[1]:.2f}; the acrylamide lane carries a DECLARED FLAT formation "
            f"term on {', '.join(SCALED_KEYS)} (multiplier {AW_MULTIPLIER:g}; envelope band "
            f"{AW_SCALE_BAND[0]:g}-{AW_SCALE_BAND[1]:g}) and a declared elimination shape on "
            f"{', '.join(ELIMINATION_KEYS)} (x{aw_elimination_multiplier(aw):.2f} here; minimum at a_w 0.82). "
            f"Source: {AW_SOURCE}."
        )
    elif state is False:
        out.append(
            f"WATER ACTIVITY (B14/B15): a_w {float(aw):.2f} is OUTSIDE the measured window "
            f"{AW_WINDOW[0]:.2f}-{AW_WINDOW[1]:.2f} ({AW_SOURCE}). No a_w term applies here: the value is "
            "recorded and changes no rate, and a comparison that moves a_w across this boundary is refused."
        )
    ph = getattr(process, "ph", None)
    if ph is not None and abs(float(ph) - REFERENCE_PH) > 1e-9:
        held = ph_in_window(ph) is False
        out.append(
            f"pH (B15): the acrylamide lane carries a DECLARED INITIAL-pH factor: x{ph_factor(ph, PH_EXPONENT_FORMATION):.2f} "
            f"on {', '.join(PH_FORMATION_KEYS)} (10^({PH_EXPONENT_FORMATION:g} (pH - {REFERENCE_PH:g})), band "
            f"{PH_EXPONENT_FORMATION_BAND[0]:g}-{PH_EXPONENT_FORMATION_BAND[1]:g}) and "
            f"x{ph_factor(ph, PH_EXPONENT_ELIMINATION):.2f} on {', '.join(PH_ELIMINATION_KEYS)} "
            f"(10^({PH_EXPONENT_ELIMINATION:g} (pH - {REFERENCE_PH:g})), band "
            f"{PH_EXPONENT_ELIMINATION_BAND[0]:g}-{PH_EXPONENT_ELIMINATION_BAND[1]:g}), measured window pH "
            f"{PH_WINDOW[0]:g}-{PH_WINDOW[1]:g}"
            + (f"; pH {float(ph):g} is OUTSIDE it and the factor is HELD at the window edge (extrapolation)" if held else "")
            + f". The source indexes its constants to the INITIAL pH of a pot that drifts by up to two units "
            f"while heating; this lane has no pH trajectory. Source: {PH_SOURCE}."
        )
    return out


def apply(
    parameters: Mapping[str, Any],
    process,
    *,
    aw_scale: Optional[float] = None,
    aw_elimination_deficit_scale: float = 1.0,
    ph_exponent_formation: float = PH_EXPONENT_FORMATION,
    ph_exponent_elimination: float = PH_EXPONENT_ELIMINATION,
) -> Tuple[Dict[str, Any], List[str]]:
    """A COPY of the acrylamide parameter dict with the declared terms applied, plus its declarations."""
    aw = getattr(process, "water_activity", None)
    ph = getattr(process, "ph", None)
    factors: Dict[str, float] = {}
    f_aw = aw_multiplier(aw, aw_scale)
    for key in SCALED_KEYS:
        factors[key] = factors.get(key, 1.0) * f_aw
    e_aw = aw_elimination_multiplier(aw, aw_elimination_deficit_scale)
    for key in ELIMINATION_KEYS:
        factors[key] = factors.get(key, 1.0) * e_aw
    for key in PH_FORMATION_KEYS:
        factors[key] = factors.get(key, 1.0) * ph_factor(ph, ph_exponent_formation)
    for key in PH_ELIMINATION_KEYS:
        factors[key] = factors.get(key, 1.0) * ph_factor(ph, ph_exponent_elimination)
    out: Dict[str, Any] = dict(parameters)
    for key, factor in factors.items():
        if factor != 1.0:
            p = out.get(key)
            if p is not None and getattr(p, "k_ref", None) is not None:
                out[key] = replace(p, k_ref=float(p.k_ref) * factor)
    return out, declarations(process)


__all__ = [
    "AW_WINDOW", "REFERENCE_AW", "SCALED_KEYS", "ELIMINATION_KEYS", "KF_TABLE", "KE_TABLE", "KF_2007", "KE_2007",
    "AW_MULTIPLIER", "AW_SCALE_BAND", "AW_ELIMINATION_TABLE", "AW_ELIMINATION_SCALE_BAND", "AW_SOURCE",
    "REFERENCE_PH", "PH_WINDOW", "PH_EXPONENT_FORMATION", "PH_EXPONENT_ELIMINATION",
    "PH_EXPONENT_FORMATION_BAND", "PH_EXPONENT_ELIMINATION_BAND", "PH_FORMATION_KEYS", "PH_ELIMINATION_KEYS",
    "KF_PH_2006", "KE_PH_2006", "PH_SOURCE",
    "in_window", "aw_multiplier", "aw_elimination_multiplier", "ph_in_window", "ph_factor", "declarations", "apply",
]
