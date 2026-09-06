"""
src/kinetic_core/parameters_dicarbonyl.py -- THE DICARBONYL TRIO (Build Wave B13, 2026-09-07).

Glucosone, glyoxal and diacetyl on the trunk lane, with the constants Kocadagli & Gokmen
2016 (J. Agric. Food Chem. 64:6446, doi 10.1021/acs.jafc.6b01862) fitted to their amine-free
glucose glass at 160 / 180 / 200 C, re-referenced from the paper's T_b = 180 C to the core's
100 C through the same helper Build Wave B7 used for the furanic channel
(`parameters_furanic._kocadagli`). Dossier: `kocadagli2016jafc_extraction.md` sec. 4 (Table 2,
glucose system, k_b in min^-1 x 10^3 with 95 % HPD).

WHY THESE THREE. Glyoxal is the CML precursor and methylglyoxal (already on the trunk) the
CEL precursor: the panel's two AGE rows are refused today for want of the species. Diacetyl
is the buttery odorant and the 3-mercapto-2-butanone precursor Yiltirak 2026 quantifies.
Glucosone is the oxidative entry that makes glyoxal at all.

WHAT IS NOT CLAIMED. No fit row of this repository measures any of the three; the constants
are one laboratory's, one amine-free matrix's, 160-200 C, extrapolated downward. Two sinks
carry the authors' own decisions (glyoxal barrier fixed to zero; diacetyl sink rate zero) and
are flagged. The wishlist lists the measurement that would replace each. The sulfur network
does NOT carry these steps (B9's topology is frozen); a sulfur wave that wants diacetyl for
3-mercapto-2-butanone adopts them then.
"""
from __future__ import annotations

from typing import Mapping, Tuple

from .parameters import KineticParameter
from .parameters_furanic import _kocadagli

#: Kocadagli 2016 JAFC Table 2, GLUCOSE system: (k_b x 1000 at 180 C, HPD), (Ea, HPD).
DICARBONYL_PARAMETERS: Mapping[str, KineticParameter] = {
    "k_glc_g": _kocadagli(
        "k_glc_g", "glucose -> glucosone", 0.069, 125.9, 4.9, 9,
        flags=("b13_dicarbonyl", "oxidative_entry_in_air"),
        note="k_b 0.069 +/- 0.005 x 1e-3 /min at 180 C; Ea 125.9 +/- 4.9. The smallest "
             "entry on the trunk by two decades; at 100-145 C it is 1e-6 to 1e-5 /min.",
    ),
    "k_g_go": _kocadagli(
        "k_g_go", "glucosone -> glyoxal + C4", 737.0, 93.8, 6.4, 10,
        flags=("b13_dicarbonyl",),
        note="k_b 737 +/- 58.9 x 1e-3 /min at 180 C; Ea 93.8 +/- 6.4. Glucosone is transient: "
             "this step is four decades faster than its formation.",
    ),
    "k_odg_da": _kocadagli(
        "k_odg_da", "1-deoxyglucosone -> diacetyl + C2", 12.2, 150.8, 8.8, 12,
        flags=("b13_dicarbonyl",),
        note="k_b 12.2 +/- 1.12 x 1e-3 /min at 180 C; Ea 150.8 +/- 8.8 -- the steepest barrier on "
             "the trunk, so diacetyl is a high-temperature product.",
    ),
    "k_go_sink": _kocadagli(
        "k_go_sink", "glyoxal -> unassigned (P3)", 32.6, 0.0, None, 15,
        flags=("b13_dicarbonyl", "ea_fixed_to_zero_by_authors"),
        note="k_b 32.6 +/- 8.83 x 1e-3 /min at 180 C; the authors FIXED the barrier to zero "
             "during estimation, so the sink runs at its 180 C rate at every temperature. "
             "Declared, flagged; the wishlist asks for a glyoxal loss rate at two temperatures.",
    ),
    "k_da_sink": _kocadagli(
        "k_da_sink", "diacetyl -> unassigned (P5)", 0.0, 0.0, None, 17,
        flags=("b13_dicarbonyl", "rate_zero_in_source"),
        note="Kocadagli step 17: 0 +/- 0 (blank Ea). Diacetyl accumulates in the source's "
             "glass; carried at zero as a PREDICTION the data may reject.",
    ),
}

DICARBONYL_KEYS: Tuple[str, ...] = tuple(DICARBONYL_PARAMETERS)

#: What would replace each declared decision (read by the wishlist through the flags).
DICARBONYL_WISHLIST: Mapping[str, str] = {
    "k_go_sink": "glyoxal loss from a glucose/glycine pot at two temperatures (a barrier for the sink)",
    "k_da_sink": "diacetyl loss from a glucose/glycine pot heated alone (the source measured none)",
    "k_glc_g": "glucosone in a glucose/glycine solution at 100-145 C (the entry is extrapolated from a 160-200 C glass)",
}

__all__ = ["DICARBONYL_KEYS", "DICARBONYL_PARAMETERS", "DICARBONYL_WISHLIST"]
