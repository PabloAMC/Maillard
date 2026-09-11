"""
src/kinetic_core/trunk_conditions.py -- WATER ACTIVITY AND pH ON THE TRUNK LANE (wave B12, 2026-09-07).

WHY THIS EXISTS
===============
Until B12 the trunk lane (glucose / fructose / glycine -> Amadori, deoxyosones, HMF, DMHF)
carried NO water-activity term and NO pH term: the engine refused every comparison that moved
either axis (four moisture claims and every trunk pH claim on the directional panel were "not
evaluable"). The corpus holds two measured shapes for exactly these axes, both read in full and
neither used by any fit until now:

  * WATER ACTIVITY -- Pereyra Gonzales, Naranjo, Leiva & Malec 2010 (Int. Dairy J. 20:40;
    `pereyragonzales2010_extraction.md` sec. 3): first-order loss of available lysine in skim
    milk powder at initial a_w 0.33 / 0.43 / 0.52 / 0.69 / 0.85 / 0.98 and 37 / 50 / 60 C, with
    95 % CI. The Maillard rate at FIXED dry-basis composition is ~3.5x faster at a_w 0.4-0.7
    than in solution (a_w 0.98) and falls toward the glass. Bell 1995 (`bell1995_extraction.md`
    sec. 4a) deconfounds the same axis in a PVP/glucose/glycine model at FIXED MOLALITY: the
    rate is flat from a_w 0.54 to 0.96 and ~5x lower in the glass at 0.33. The two agree on the
    glass collapse and disagree on the plateau because they hold different things fixed --
    dry-basis composition (the practitioner's variable) versus aqueous molality (the mass-action
    variable). The multiplier below is Pereyra Gonzales's NET shape; its band reaches down to
    Bell's no-effect plateau, so the interval spans the disagreement instead of hiding it.

  * pH -- Martins & van Boekel 2003, Part II (Carbohydr. Res. 338:1665; `martins2003_extraction.md`
    sec. 5): the degradation of the glucose/glycine Amadori compound (DFG) fitted at 100 and
    120 C x pH 5.5 and 6.8. Every DFG-consuming step is faster at pH 6.8: k1 (1,2-enolisation)
    x3.0 / x8.0, k2 (2,3-enolisation) x15.6 / x7.3, k3 (direct fragmentation) x8.6 / x9.8 at
    100 / 120 C, i.e. 0.37-0.92 decades per pH unit, mean 0.69. The three Amadori-decay steps of
    the trunk (`r_ama_tdg`, `r_ama_odg`, `r_ama_mgo`) take that exponent about the reference
    pH 6.8 of every Martins 2005 constant. The amine-free caramelisation entries are NOT given
    a pH term: no source in the corpus contrasts them across pH (Agcam 2022 is one pH, 3.5).

WHAT THIS MODULE DOES, AND DOES NOT DO
======================================
`apply(parameters, process)` returns a COPY of the trunk parameter dict with `k_ref` scaled on
the named steps, plus the declarations the run must print. At a_w None or 0.98 and pH 6.8 every
factor is exactly 1.0, so every wave up to B10 is reproduced bit for bit (pinned by
`tests/unit/test_kinetic_core_b12.py`). The factors are DECLARED CONSTANTS with bands, not
fitted: they come from measured within-study ratios, the same standing as a measured Ea
override. They apply to the trunk lane run on its own; the copy of the trunk steps inside the
sulfur network keeps its own (pH-state) machinery and no a_w term -- adopting these there is a
sulfur wave, not this one. The acrylamide lane carries neither term (De Vleeschouwer 2009, on
disk and unread, is its a_w source).
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

#: The trunk's reference conditions: every Martins 2005 constant was measured at pH 6.8 in a
#: 0.1 M phosphate solution (a_w ~0.98).
REFERENCE_PH = 6.8
REFERENCE_AW = 0.98

# ---------------------------------------------------------------------------
# Water activity
# ---------------------------------------------------------------------------

#: Pereyra Gonzales 2010 Table 1, k(a_w) / k(0.98) at 50 and 60 C, averaged; the 37 C row is
#: excluded at a_w 0.33 (glassy: 0.018x) and its ratios elsewhere agree with 50/60 C within 20 %.
#:   50 C: 7.14/1.89 = 3.78 (0.33), 3.83 (0.43), 3.36 (0.52), 3.75 (0.69), 2.41 (0.85), 1 (0.98)
#:   60 C: 17.46/7.54 = 2.32 (0.33), 3.63 (0.43), 3.45 (0.52), 3.61 (0.69), 2.76 (0.85), 1 (0.98)
AW_MULTIPLIER_TABLE: Tuple[Tuple[float, float], ...] = (
    (0.33, 3.05),
    (0.43, 3.73),
    (0.52, 3.40),
    (0.69, 3.68),
    (0.85, 2.58),
    (0.98, 1.00),
)
#: The band, expressed as a SCALE on the excess (m - 1): 0 is Bell 1995's fixed-molality plateau
#: (no effect), 1 the Pereyra Gonzales centre, 1.2 the source's own 95 % CI. The envelope draws the
#: scale uniformly over AW_SCALE_BAND (`CoreDraw.trunk_aw_scale`); the multiplier is 1 + (m - 1) s.
AW_SCALE_BAND = (0.0, 1.2)
AW_BAND_LOW_MULTIPLIER = 1.0
AW_BAND_HIGH_FACTOR = 1.2
#: Below the lowest measured a_w the multiplier is HELD at the 0.33 value with a glass warning:
#: Pereyra Gonzales's 37 C point (0.018x) and Bell's glass (0.2x) say the collapse is real and
#: temperature-dependent, and nothing here sizes it.
AW_TABLE_FLOOR = 0.33
#: The steps the multiplier scales: the amine-sugar condensation is what lysine loss measures.
AW_STEPS: Tuple[str, ...] = ("k_schiff",)
AW_SOURCE = (
    "Pereyra Gonzales, Naranjo, Leiva & Malec 2010, Int. Dairy J. 20:40-45, Table 1: first-order "
    "available-lysine loss in skim milk powder at a_w 0.33-0.98, 37-60 C, 95 % CI "
    "(pereyragonzales2010_extraction.md sec. 3); plateau alternative from Bell 1995, J. Food Sci. "
    "(bell1995_extraction.md sec. 4a)"
)


def aw_multiplier(water_activity: Optional[float]) -> float:
    """Piecewise-linear in a_w through the table; 1.0 at None or >= 0.98; held below 0.33."""
    if water_activity is None:
        return 1.0
    aw = float(water_activity)
    if aw >= REFERENCE_AW:
        return 1.0
    if aw <= AW_TABLE_FLOOR:
        return AW_MULTIPLIER_TABLE[0][1]
    for (a0, m0), (a1, m1) in zip(AW_MULTIPLIER_TABLE[:-1], AW_MULTIPLIER_TABLE[1:]):
        if a0 <= aw <= a1:
            return m0 + (m1 - m0) * (aw - a0) / (a1 - a0)
    raise AssertionError(aw)


def aw_band(water_activity: Optional[float]) -> Tuple[float, float]:
    m = aw_multiplier(water_activity)
    if m == 1.0:
        return (1.0, 1.0)
    return (1.0 + (m - 1.0) * AW_SCALE_BAND[0], 1.0 + (m - 1.0) * AW_SCALE_BAND[1])


# ---------------------------------------------------------------------------
# pH
# ---------------------------------------------------------------------------

#: Martins & van Boekel 2003 Part II Table 3, k(pH 6.8) / k(pH 5.5) per DFG-consuming step:
#:   100 C: k1 0.57/0.19 = 3.0, k2 1.56/0.10 = 15.6, k3 1.55/0.18 = 8.6
#:   120 C: k1 8.89/1.11 = 8.0, k2 6.29/0.86 = 7.3, k3 8.62/0.88 = 9.8
#: -> decades per pH unit over 1.3 units: 0.37, 0.92, 0.72 (100 C); 0.69, 0.66, 0.76 (120 C).
PH_EXPONENT_DECADES_PER_UNIT = 0.69
PH_EXPONENT_BAND = (0.37, 0.92)
#: The measured window; outside it the factor is an extrapolation and says so.
PH_MEASURED_WINDOW = (5.5, 6.8)
#: The three Amadori-decay steps of the trunk.
PH_STEPS: Tuple[str, ...] = ("k_ama_tdg", "k_ama_odg", "k_ama_mgo")
PH_SOURCE = (
    "Martins & van Boekel 2003, Carbohydr. Res. 338:1665 (Part II), Table 3: DFG degradation "
    "constants at 100 / 120 C x pH 5.5 / 6.8, 95 % HPD (martins2003_extraction.md sec. 5)"
)

# ---------------------------------------------------------------------------
# B40 (2026-09-11): THE 3-DEOXYGLUCOSONE EXITS' pH TERM. B12 scaled the three Amadori-decay steps
# from Martins 2003 Table 3 and left the 3-DG exits at their pH-6.8 values everywhere. The same
# table prints those exits at pH 5.5: k6 (3-DG -> formic acid) is 14x slower at 100 C and 7x at
# 120 C; k5 (3-DG -> fragments) 6.6x at 100 C. Declared slopes in decades per pH unit, reference
# pH 6.8, the same measured window. INERT until B40 ships (THREE_DEOXY_EXIT_PH_TERM False).
# Pre-registration: results/validation/kinetic_core_b40_prereg.md.
# ---------------------------------------------------------------------------
#: SHIPPED by wave B41 (2026-09-11, kinetic_core_b41_ship_rule.json: SHIP). B40 tried the term on BOTH
#: exits and the Leitzen hold-out rejected the one on k_tdg_mgo (its methylglyoxal row went 1.28x ->
#: 33x) while keeping the one on k_tdg_fa; B41 ships the formic-acid exit's term alone.
THREE_DEOXY_EXIT_PH_TERM: bool = True
#: key -> (exponent, band, source rows)
THREE_DEOXY_EXIT_PH: Mapping[str, Tuple[float, Tuple[float, float], str]] = {
    "k_tdg_fa": (0.77, (0.65, 0.89), "Martins 2003 Table 3 k6: 1.9e-3 vs 2.74e-2 (100 C), 4.30e-2 vs 3.04e-1 (120 C), pH 5.5 vs 6.8"),
}
#: The declaration B40 made and the hold-out rejected, kept as the record and NOT applied: Martins'
#: k5 is a lumped "3-DG -> fragments" step, and transferring its pH slope to Kocadagli's amine-free
#: methylglyoxal route moved Leitzen 2021's methylglyoxal from 1.28x to 33x.
THREE_DEOXY_EXIT_PH_REJECTED_B40: Mapping[str, Tuple[float, Tuple[float, float], str]] = {
    "k_tdg_mgo": (0.63, (0.37, 0.92), "Martins 2003 Table 3 k5 at 100 C only: 1.38e-2 vs 9.07e-2; REJECTED by the Leitzen methylglyoxal row in B40"),
}


def ph_factor(ph: Optional[float], exponent: float = PH_EXPONENT_DECADES_PER_UNIT) -> float:
    if ph is None:
        return 1.0
    return 10.0 ** (float(exponent) * (float(ph) - REFERENCE_PH))


def ph_band(ph: Optional[float]) -> Tuple[float, float]:
    if ph is None or abs(float(ph) - REFERENCE_PH) < 1e-9:
        return (1.0, 1.0)
    lo, hi = (ph_factor(ph, e) for e in PH_EXPONENT_BAND)
    return (min(lo, hi), max(lo, hi))


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------


def _scaled(parameters: Mapping[str, Any], keys: Sequence[str], factor: float) -> Dict[str, Any]:
    out = dict(parameters)
    for key in keys:
        p = out.get(key)
        if p is None:
            raise KeyError(f"{key!r}: not in the trunk parameter dict")
        if p.k_ref is None:
            raise ValueError(f"{key!r}: unpopulated; cannot scale")
        out[key] = replace(p, k_ref=float(p.k_ref) * float(factor))
    return out


def pyrazine_factor(process, slopes=None) -> Tuple[float, List[str]]:
    """B18: the pyrazine step's own pH factor (two slopes, knot at 7, reference pH 8) and its declaration."""
    from .parameters_pyrazine import (
        PYRAZINE_PH_SLOPES, PYRAZINE_PH_SOURCE, PYRAZINE_REFERENCE_PH, pyrazine_ph_factor,
    )

    ph = getattr(process, "ph", None)
    s = tuple(PYRAZINE_PH_SLOPES if slopes is None else slopes)
    f = pyrazine_ph_factor(ph, s)
    notes: List[str] = []
    if abs(f - 1.0) > 1e-12:
        notes.append(
            f"PYRAZINE pH TERM (B18): pH {float(ph):g} scales the three pyrazine steps by x{f:.3g} relative to "
            f"pH {PYRAZINE_REFERENCE_PH:g} ({s[0]:.2f} decades per unit above 7, {s[1]:.2f} below; fitted on "
            f"{PYRAZINE_PH_SOURCE}). Below pH 5 and above 9 the term is an extrapolation of a three-point ladder."
        )
    return f, notes


def factors(
    process, *, aw_scale: float = 1.0, ph_exponent: float = PH_EXPONENT_DECADES_PER_UNIT,
) -> Tuple[float, float, List[str]]:
    """(a_w multiplier, pH factor, declarations) for one process; 1.0, 1.0, [] at the references."""
    warnings: List[str] = []
    aw = getattr(process, "water_activity", None)
    m = aw_multiplier(aw)
    m_eff = 1.0
    if m != 1.0:
        m_eff = 1.0 + (m - 1.0) * float(aw_scale)
        lo, hi = aw_band(aw)
        warnings.append(
            f"WATER ACTIVITY TERM (B12): a_w {float(aw):.2f} scales the amine-sugar condensation "
            f"by x{m_eff:.2f} (declared band {lo:.2f}-{hi:.2f}; {AW_SOURCE}). The shape is the net "
            f"effect at FIXED dry-basis composition measured at 37-60 C in one matrix; the band's "
            f"floor is the fixed-molality no-effect alternative. Dehydration steps carry no a_w term."
        )
        if float(aw) <= AW_TABLE_FLOOR:
            warnings.append(
                f"a_w {float(aw):.2f} is at or below the lowest measured point (0.33): the multiplier is "
                f"HELD there. Both sources measure a further collapse in the glass (0.02-0.2x) that "
                f"depends on temperature and is not sized here."
            )
    ph = getattr(process, "ph", None)
    f = ph_factor(ph, ph_exponent)
    if abs(f - 1.0) > 1e-12:
        lo, hi = ph_band(ph)
        warnings.append(
            f"pH TERM (B12): pH {float(ph):g} scales the three Amadori-decay steps by x{f:.2f} "
            f"(10^({ph_exponent:.2f} per pH unit); declared band {lo:.2f}-{hi:.2f}; {PH_SOURCE}). "
            f"The amine-free caramelisation entries and every downstream step carry no pH term: no "
            f"source contrasts them across pH."
        )
        if not (PH_MEASURED_WINDOW[0] <= float(ph) <= PH_MEASURED_WINDOW[1]):
            warnings.append(
                f"pH {float(ph):g} is outside the measured window {PH_MEASURED_WINDOW[0]}-"
                f"{PH_MEASURED_WINDOW[1]}: the factor is an extrapolation of a two-point contrast."
            )
    return m_eff, f, warnings


def declarations(process) -> List[str]:
    """The declarations alone, for the envelope declaration a run prints."""
    return factors(process)[2]


# ===========================================================================
# WAVE B29 (2026-09-10): THE FIRST OXYGEN AXIS ON THE TRUNK
# ===========================================================================
# Until this wave every trunk prediction was at whatever oxygen the fits happened to have, and
# nothing in the output said so. Two constants were already flagged as oxidative -- `k_glc_g`
# carries `oxidative_entry_in_air` and `r_ama_g` is the Amadori route to the same glucosone -- so
# the branch existed and had no lever.
#
# WHY THE LEVER GOES ON EXACTLY THESE TWO STEPS, AND IT IS THE SOURCE'S OWN CONTROL THAT SAYS SO.
# Hofmann & Schieberle's companion paper feeds the DICARBONYLS directly and finds the Strecker
# ALDEHYDE oxygen-INDEPENDENT (while the Strecker ACID is not). So the aldehyde's oxygen dependence
# is entirely UPSTREAM of the dicarbonyl -- in getting from the sugar or the Amadori compound to it
# -- which is where these two entries sit. Putting the lever anywhere else would contradict a
# measured control.
#: The two oxidative entries, and nothing else.
OXIDATIVE_ENTRY_STEPS: Tuple[str, ...] = ("k_ama_g", "k_glc_g")
#: AIR IS THE REFERENCE AND ITS FACTOR IS EXACTLY 1, by definition rather than by fit: every fit row
#: in this model was run in a closed vial in air, so air is the atmosphere the constants already
#: describe. A pot that declares no atmosphere is air, and the answer is bit-for-bit what it was
#: before this wave existed.
ATMOSPHERE_REFERENCE = "air"
ATMOSPHERE_VALUES: Tuple[str, ...] = ("argon", "air", "air_cu")
#: The fitted multipliers. Empty until a B29 report supplies them; at the default every atmosphere
#: other than air RAISES a refusal rather than silently returning the air answer.
ATMOSPHERE_FACTORS: Dict[str, float] = {"air": 1.0}
ATMOSPHERE_SOURCE = (
    "Hofmann & Schieberle 2000b Table 2 (ARP-Phe and glucose + Phe, 1 mmol each in 10 mL of 0.5 M "
    "phosphate pH 7.0, 100 C, 120 min, closed vial under argon / air / air + 5 mmol/L Cu(II)): the "
    "Strecker aldehyde is 9.2x higher under air than argon from the Amadori compound and 3.5x from "
    "the sugar pot, and 2.5x / 1.9x higher again with copper. hofmann2000b_extraction.md Table 2"
)
ATMOSPHERE_TRANSFER = (
    "Hofmann's amino acid is PHENYLALANINE and this model's is glycine, so what transfers is not a "
    "yield but the ratio of one pot TO ITSELF under two atmospheres. The argument that such a ratio "
    "transfers is the companion paper's fed-dicarbonyl control -- the aldehyde is oxygen-independent "
    "once the dicarbonyl is supplied -- which places the sensitivity in the sugar chemistry, upstream "
    "of the amino acid. It is an argument from a measurement, not an assumption, and it is testable: "
    "one multiplier has to explain BOTH 9.2 and 3.5."
)


def atmosphere_factor(process, factors_table=None) -> Tuple[float, List[str]]:
    """
    The multiplier on the two oxidative entries for this pot's declared atmosphere.

    A pot with no declared atmosphere is AIR and gets exactly 1.0 with no warning: that is what
    every constant in this model already describes. Any other atmosphere with no fitted factor
    RAISES, because silently returning the air answer for an argon pot would be inventing a
    number -- the failure mode this repository refuses everywhere else.
    """
    table = dict(ATMOSPHERE_FACTORS if factors_table is None else factors_table)
    name = getattr(process, "atmosphere", None)
    if name is None or str(name) == ATMOSPHERE_REFERENCE:
        return 1.0, []
    key = str(name)
    if key not in ATMOSPHERE_VALUES:
        raise ValueError(
            f"atmosphere {key!r} is not one of {ATMOSPHERE_VALUES}. The axis has three settings "
            f"because one paper measured three; it is not a continuous oxygen partial pressure."
        )
    if key not in table:
        raise ValueError(
            f"atmosphere {key!r} has no fitted multiplier. Wave B29 measures it "
            f"(kinetic_core_b29_prereg.md); until that report ships, only {ATMOSPHERE_REFERENCE!r} "
            f"can be answered. Returning the air answer for an argon pot would invent a number."
        )
    return float(table[key]), [
        f"ATMOSPHERE {key}: the two oxidative entries to glucosone are scaled by "
        f"{table[key]:.3g}. {ATMOSPHERE_TRANSFER}"
    ]


def apply(
    parameters: Mapping[str, Any],
    process,
    *,
    aw_scale: float = 1.0,
    ph_exponent: float = PH_EXPONENT_DECADES_PER_UNIT,
    pyrazine_slopes=None,
    atmosphere_factors=None,
) -> Tuple[Dict[str, Any], List[str]]:
    """
    The trunk parameter dict with the condition terms applied, and the declarations.

    ``aw_scale`` and ``ph_exponent`` are the envelope's hooks: a draw moves the multiplier
    within its band and the exponent within its band. At the defaults the factors are the
    declared centres; at a_w None / >= 0.98 and pH 6.8 they are exactly 1.0. ``pyrazine_slopes``
    (B18) is the pyrazine step's own two-slope pH term, the fitted slopes by default; the fit
    generator passes candidates. It scales the three pyrazine steps only, when they are present.
    """
    m_eff, f, warnings = factors(process, aw_scale=aw_scale, ph_exponent=ph_exponent)
    out: Dict[str, Any] = dict(parameters)
    if m_eff != 1.0:
        out = _scaled(out, AW_STEPS, m_eff)
    if abs(f - 1.0) > 1e-12:
        out = _scaled(out, PH_STEPS, f)
        if THREE_DEOXY_EXIT_PH_TERM:
            ph = getattr(process, "ph", None)
            for key, (exponent, band, source) in THREE_DEOXY_EXIT_PH.items():
                if key in out and ph is not None:
                    fx = ph_factor(ph, exponent)
                    out = _scaled(out, (key,), fx)
                    warnings.append(
                        f"pH TERM ON THE 3-DG EXIT (B40): pH {ph:g} scales {key} by x{fx:.3g} "
                        f"(10^({exponent} per pH unit); declared band {band[0]}-{band[1]}; {source})."
                    )
    from .parameters_pyrazine import PYRAZINE_PH_STEPS

    if all(key in out for key in PYRAZINE_PH_STEPS):
        fp, notes = pyrazine_factor(process, pyrazine_slopes)
        if abs(fp - 1.0) > 1e-12:
            out = _scaled(out, PYRAZINE_PH_STEPS, fp)
            # B22: methionine's Strecker steps are glycine's times an identity ratio; the same pH term.
            from .parameters_methionine import METHIONINE_PH_STEPS

            if all(key in out for key in METHIONINE_PH_STEPS):
                out = _scaled(out, METHIONINE_PH_STEPS, fp)
            # B24: the proline Strecker step, likewise.
            from .parameters_proline import PROLINE_PH_STEPS

            if all(key in out for key in PROLINE_PH_STEPS):
                out = _scaled(out, PROLINE_PH_STEPS, fp)
        warnings = list(warnings) + notes
    # B29: the oxygen axis. Exactly 1.0 for air, which is every fit row in this model.
    fo, onotes = atmosphere_factor(process, atmosphere_factors)
    if abs(fo - 1.0) > 1e-12 and all(key in out for key in OXIDATIVE_ENTRY_STEPS):
        out = _scaled(out, OXIDATIVE_ENTRY_STEPS, fo)
    warnings = list(warnings) + onotes
    return out, warnings


__all__ = [
    "AW_BAND_HIGH_FACTOR", "AW_BAND_LOW_MULTIPLIER", "AW_SCALE_BAND", "AW_MULTIPLIER_TABLE", "AW_SOURCE", "AW_STEPS",
    "AW_TABLE_FLOOR", "PH_EXPONENT_BAND", "PH_EXPONENT_DECADES_PER_UNIT", "PH_MEASURED_WINDOW",
    "PH_SOURCE", "PH_STEPS", "REFERENCE_AW", "REFERENCE_PH", "apply", "aw_band", "aw_multiplier", "pyrazine_factor",
    "ATMOSPHERE_FACTORS", "ATMOSPHERE_REFERENCE", "ATMOSPHERE_SOURCE", "ATMOSPHERE_TRANSFER",
    "ATMOSPHERE_VALUES", "OXIDATIVE_ENTRY_STEPS", "atmosphere_factor",
    "declarations", "factors", "ph_band", "ph_factor",
]
