"""
src/kinetic_core/parameters_pyrazine.py -- THE PYRAZINE STEP (Build Wave B18, 2026-09-08).

Five steps on the trunk lane (network.PYRAZINE_REACTIONS): two Strecker deaminations (glyoxal or
methylglyoxal + glycine -> aminoketone + CO2 + formaldehyde; second order, dicarbonyl x glycine,
L/(mmol*min) at the trunk's 100 C reference) and three aminoketone condensations (2 AKG -> pyrazine,
2 AKM -> 2,5-dimethylpyrazine, AKG + AKM -> 2-methylpyrazine) on ONE shared constant declared fast.
With the condensation fast the aminoketones sit at steady state and the pyrazine formation rate is
the Strecker rate over two, which is the "zero-order in product" rate Zhou 2024 printed; the rate
law dicarbonyl x amino acid is Yu 2018's (yu2018_extraction.md sec. 4).

WHAT IS FITTED (wave B18, pre-registered in results/validation/kinetic_core_b18_prereg.md): the two
Strecker constants' log10 k at 100 C and their barriers, against Zhou 2024's six formation rates
(fed 20 mmol/L glyoxal or methylglyoxal + 20 mmol/L alanine, water, initial pH 8, 100 / 110 / 120 C;
zhou2024_extraction.md Table 2), and the two slopes of the pH term against Leahy & Reineccius 1989's
four within-study pH ratios (leahy1989a_extraction.md). The fitted values below are FROZEN LITERALS
asserted equal to the B18 fit report by tests/unit/test_kinetic_core_b18.py (the B7 pattern), so
every lane's parameter set carries them without reading a report at import.

WHAT IS DECLARED. (i) Alanine -> glycine: Zhou fed alanine, the trunk holds glycine; the products
are the same (the ring carbons are the dicarbonyl's, the amino acid gives the nitrogen) and the rate
is not measured for the pair; a +/- 0.5 dex band travels on every answer that uses the step
(`PYRAZINE_TRANSFER_BAND_DECADES`). (ii) The condensation is FAST: Jousse 2002 write the aminoketone
condensation "I + I -> pyrazines" as "fast" (jousse2002_extraction.md, Table of R10) and Zhou's
products rise linearly from t = 0; `k_cond` is set where the aminoketone lifetime is minutes at the
fitted fluxes, its value is a declared assumption, and the ship rule reports the sensitivity. The
mixed pyrazine then follows the two pools statistically: rate_MPZ = 2 sqrt(rate_PZ rate_DMP), which
is the pre-registration's "geometric mean" made mechanistic (with its factor two). (iii) The pH term
is piecewise-linear in log10 k against pH with a knot at 7, exactly 1 at the trunk's reference pH
6.8 where the stored constants apply, on the two Strecker steps (trunk_conditions.pyrazine_factor);
Zhou's pH-8 pots carry its factor in the fit like any other pot.
"""
from __future__ import annotations

import math
from typing import Dict, Mapping, Tuple

from .parameters import AW_OF_MEASUREMENT, KineticParameter

#: The pH at which the stored constants apply: the trunk's own reference (trunk_conditions.REFERENCE_PH),
#: so every condition term is exactly 1 at the references. Zhou 2024's pots (initial pH 8.0) carry the
#: term's factor from 6.8 to 8 in the fit, as any other pot does.
PYRAZINE_REFERENCE_PH = 6.8
#: The pH of the pots the Strecker constants were fitted on (Zhou 2024, unbuffered, initial 8.0).
PYRAZINE_FIT_PH = 8.0
#: The knot of the two-slope pH term (Leahy's arms are 9, 7 and 5).
PYRAZINE_PH_KNOT = 7.0
#: The declared alanine -> glycine transfer band, in decades, on every pyrazine answer.
PYRAZINE_TRANSFER_BAND_DECADES = 0.5
#: The declared fast condensation, L/(mmol*min), no barrier: at the fitted Strecker fluxes the
#: aminoketone steady-state lifetime is 1-3 min at 100-120 C, so the Strecker step is rate-determining.
K_COND_DECLARED_L_PER_MMOL_MIN = 1.0e4

# ---------------------------------------------------------------------------
# THE FROZEN B18 VALUES. Replaced by the fit; asserted against the report.
# Start values (prior centres): Zhou 2024's printed rates re-expressed as second-order constants
# at 20 mM x 20 mM (6.98e-8 and 8.75e-9 L/(mmol*min) at 100 C, zhou2024_extraction.md sec. 4),
# times two because the measured product rate is the Strecker rate over two, brought from Zhou's
# pH 8 to the trunk's reference pH 6.8 with the prior slopes; the printed barriers 100.59 and
# 111.66 kJ/mol; and Leahy's ratios (k(9)/k(7) ~ 2.7, k(7)/k(5) ~ 14) as slopes of about 0.22 and
# 0.58 decades per pH unit.
# ---------------------------------------------------------------------------
#: The prior centres the fit started from (Zhou at pH 8 brought to 6.8 with slopes 0.22 / 0.58, times two):
#: log10 k_go_ak -6.52, log10 k_mgo_ak -7.42. The B18 fit report's optimum (2026-09-08; cost 2.64 on 10
#: rows, reduced chi2 0.66; both starts agree to 2e-8):
FROZEN_B18: Mapping[str, float] = {
    "log10_k_go_ak_100C": -6.5415452891122525,
    "ea_go_ak_kj_mol": 103.0999999999812,
    "log10_k_mgo_ak_100C": -7.529714282148093,
    "ea_mgo_ak_kj_mol": 114.89999999999999,
    "ph_slope_above_7_decades_per_unit": 0.19744840125387295,
    "ph_slope_below_7_decades_per_unit": 0.5795142879227316,
}
#: What the ship rule found and every pyrazine answer must carry (kinetic_core_b18_ship_rule.md):
#: (a) the B13 dry-glass glyoxal sink empties a fed 20 mM glyoxal pot in two hours at 100 C, so
#: log10 k_go_ak is 0.59 dex higher than it would be over a constant pool (the nosink variant);
#: (b) from a sugar + amine pot the pyrazine yield follows the trunk's dicarbonyl SUPPLY, which misses
#: Leahy 1989's 95 C total by 2.9 decades, makes almost no pyrazine (glyoxal), and runs 300-460 kJ/mol
#: apparent barriers against 100-180 measured: the step is measured, the supply is not.
PYRAZINE_SUPPLY_CAVEAT = (
    "PYRAZINES (B18): the two Strecker constants are MEASURED on fed glyoxal / methylglyoxal + alanine at "
    "100-120 C (Zhou 2024; alanine -> glycine declared, +/- 0.5 dex). From a sugar + amine pot the yield "
    "follows the trunk's dicarbonyl SUPPLY. Since B21 (2026-09-09) the glyoxal supply in water is fitted "
    "(the Amadori compound's route to glucosone, Hamzalioglu 2026; Quan 2020's glyoxal level reproduced "
    "within 0.35 dex) and pyrazine itself is no longer absent; the methylglyoxal supply is Martins' aqueous "
    "step. Against Leahy 1989 (lysine + glucose, 95 C, 2 h, pH 9) the model's total is still 2.8 decades "
    "low with apparent barriers of 300-460 kJ/mol against 150-180 measured, so what remains is the "
    "Strecker step itself at 95 C and pH 9, or lysine against glycine, not the glyoxal. Trust a pyrazine "
    "number from this model only as a fed-dicarbonyl statement."
)
PYRAZINE_SINK_CAVEAT = (
    "PYRAZINES (B18): the glyoxal Strecker constant was fitted with the B13 glyoxal sink in force (a dry-glass "
    "rate at 180 C, barrier fixed to zero), which removes 98 % of a fed 20 mM glyoxal in two hours at 100 C; "
    "with that sink zeroed the constant is 0.59 dex lower (the wave's information-only variant). "
    "A pyrazine answer inherits that conditionality until a wave measures the dicarbonyl sinks in water."
)

_ZHOU = ("Zhou, Hu, Cui, Hussain et al. 2024, J. Agric. Food Chem. 72:18630 (doi 10.1021/acs.jafc.4c03706), "
         "Table 2: pyrazine 0.0279 / 0.0791 / 0.1507 and 2,5-dimethylpyrazine 0.0035 / 0.0100 / 0.0230 umol L-1 min-1 "
         "at 100 / 110 / 120 C from 20 mM alanine + 20 mM glyoxal or methylglyoxal, water, initial pH 8; "
         "Ea 100.59 / 111.66 kJ/mol")
_ZHOU_DOSSIER = "zhou2024_extraction.md sec. 4 (the second-order re-expression is the dossier's, marked derived_assumption)"
_LEAHY = ("Leahy & Reineccius 1989, ACS Symp. Ser. 409 ch. 18 (doi 10.1021/bk-1989-0409.ch018), Table I: "
          "lysine + glucose pyrazine 3.596 / 1.346 / 0.0938 ppm/h at pH 9 / 7 / 5, 95 C")
_JOUSSE = ("Jousse, Jongen, Agterof, Russell & Braat 2002, J. Food Sci. 67:2534, Table 1 / reaction R10: the aminoketone "
           "condensation I + I -> pyrazines, second order, 'fast' (jousse2002_extraction.md sec. 3)")


def _strecker(key: str, transformation: str, log10_k: float, ea: float, note: str) -> KineticParameter:
    return KineticParameter(
        key=key, transformation=transformation, k_ref=10.0 ** float(log10_k), ea_kj_mol=float(ea),
        unit="L/(mmol*min)", order=2, evidence_class="derived_from_fit_data",
        source_anchor=_ZHOU, dossier_anchor=_ZHOU_DOSSIER,
        conditions="water, initial pH 8 (unbuffered), 100-120 C, 20 + 20 mmol/L, conversion < 0.2 %",
        ph_of_measurement=PYRAZINE_FIT_PH, temperature_range_c=(100.0, 120.0),
        rate_transfer="licensed_at_measurement_ph_only", aw_of_measurement=AW_OF_MEASUREMENT,
        flags=("b18_pyrazine", "alanine_to_glycine_declared", "fitted_wave_b18", "conditional_on_b13_dicarbonyl_sinks"), note=note,
    )


K_COND: KineticParameter = KineticParameter(
    key="k_cond", transformation="2 aminoketones -> pyrazine (condensation, dehydration, oxidation; net; shared by the three pairings)",
    k_ref=K_COND_DECLARED_L_PER_MMOL_MIN, ea_kj_mol=0.0, unit="L/(mmol*min)", order=2,
    evidence_class="bounded_from_a_timescale_bracket",
    source_anchor="DECLARED FAST: " + _JOUSSE + "; Zhou 2024's pyrazines rise linearly from t = 0 (Figure 2, figure-only)",
    dossier_anchor="results/validation/kinetic_core_b18_prereg.md sec. 2; jousse2002_extraction.md",
    conditions="not rate-determining by construction; the ship rule reports the sensitivity to a tenfold change",
    ph_of_measurement=None, temperature_range_c=(75.0, 150.0), rate_transfer="not_licensed",
    aw_of_measurement=AW_OF_MEASUREMENT,
    flags=("b18_pyrazine", "declared_fast", "no_measured_rate", "shared_by_three_pairings"),
    note="A declared assumption, not a measurement: any value above about 1e3 gives the same pyrazine yields; the "
         "mixed pyrazine follows the two aminoketone pools statistically because the constant is shared.",
)


def with_fitted_pyrazine(log10_k_go_ak: float, ea_go_ak: float, log10_k_mgo_ak: float, ea_mgo_ak: float) -> Dict[str, KineticParameter]:
    """The pyrazine block at arbitrary Strecker values (the fit generator's hook and the report reader's)."""
    return {
        "k_go_ak": _strecker("k_go_ak", "glyoxal + glycine -> aminoacetaldehyde + CO2 + HCHO (Strecker, net)", log10_k_go_ak, ea_go_ak,
                             "B18 fit row source: the three glyoxal-route rates. Alanine -> glycine declared (+/- 0.5 dex)."),
        "k_mgo_ak": _strecker("k_mgo_ak", "methylglyoxal + glycine -> aminoacetone + CO2 + HCHO (Strecker, net)", log10_k_mgo_ak, ea_mgo_ak,
                              "B18 fit row source: the three methylglyoxal-route rates. Alanine -> glycine declared (+/- 0.5 dex)."),
        "k_cond": K_COND,
    }


PYRAZINE_PARAMETERS: Mapping[str, KineticParameter] = with_fitted_pyrazine(
    FROZEN_B18["log10_k_go_ak_100C"], FROZEN_B18["ea_go_ak_kj_mol"],
    FROZEN_B18["log10_k_mgo_ak_100C"], FROZEN_B18["ea_mgo_ak_kj_mol"],
)
PYRAZINE_KEYS: Tuple[str, ...] = tuple(PYRAZINE_PARAMETERS)
#: The steps the pyrazine pH term scales: the two Strecker steps (the amine-dependent ones).
PYRAZINE_PH_STEPS: Tuple[str, ...] = ("k_go_ak", "k_mgo_ak")
PYRAZINE_PH_SLOPES: Tuple[float, float] = (
    FROZEN_B18["ph_slope_above_7_decades_per_unit"], FROZEN_B18["ph_slope_below_7_decades_per_unit"],
)
PYRAZINE_PH_SOURCE = _LEAHY + " (four within-study ratios, k(7)/k(9) and k(5)/k(9) for pyrazine and methylpyrazine)"


def _log10_shape(ph: float, slopes: Tuple[float, float]) -> float:
    """The piecewise-linear log10 k against pH, knot at 7, zero at the knot."""
    s_hi, s_lo = float(slopes[0]), float(slopes[1])
    x = float(ph)
    return s_hi * (x - PYRAZINE_PH_KNOT) if x >= PYRAZINE_PH_KNOT else s_lo * (x - PYRAZINE_PH_KNOT)


def pyrazine_ph_factor(ph, slopes: Tuple[float, float] = PYRAZINE_PH_SLOPES) -> float:
    """10^(log10 k(pH) - log10 k(6.8)): two slopes with a knot at 7; exactly 1 at the trunk's reference pH."""
    if ph is None:
        return 1.0
    return 10.0 ** (_log10_shape(ph, slopes) - _log10_shape(PYRAZINE_REFERENCE_PH, slopes))


#: What would replace each declared decision (read by the wishlist through the flags).
PYRAZINE_WISHLIST: Mapping[str, str] = {
    "k_cond": "an aminoketone time course, or methylpyrazine from fed glyoxal + methylglyoxal with an amino acid (the mixed condensation)",
    "k_go_ak": "the same fed-glyoxal ladder with glycine (Zhou 2024 used alanine): the transfer band would close",
    "k_mgo_ak": "the same fed-methylglyoxal ladder with glycine; and any of the three at pH 5-6 in a buffer",
}

__all__ = ["FROZEN_B18", "K_COND", "PYRAZINE_FIT_PH", "PYRAZINE_SINK_CAVEAT", "PYRAZINE_SUPPLY_CAVEAT", "K_COND_DECLARED_L_PER_MMOL_MIN", "PYRAZINE_KEYS", "PYRAZINE_PARAMETERS",
           "PYRAZINE_PH_KNOT", "PYRAZINE_PH_SLOPES", "PYRAZINE_PH_SOURCE", "PYRAZINE_PH_STEPS",
           "PYRAZINE_REFERENCE_PH", "PYRAZINE_TRANSFER_BAND_DECADES", "PYRAZINE_WISHLIST",
           "pyrazine_ph_factor", "with_fitted_pyrazine"]
