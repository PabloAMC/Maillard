"""
src/kinetic_core/species_lipid.py

THE LIPID-OXIDATION STATE VECTOR (Build Wave B6, 2026-08-29).
=============================================================================

Module 5 of the kinetic-core rebuild. This module declares the species; the
branch topology lives in ``lipid.py`` and every parameter in
``parameters_lipid.py``.

WHY THE HYDROPEROXIDE POOL IS SPLIT FOUR WAYS
---------------------------------------------
Frankel 1989 -- the ONLY measured branch distribution in the corpus, and this
module's declared FIT source -- does not measure "the hexanal share". It
measures three DIFFERENT distributions from three DIFFERENT hydroperoxide
preparations, and the hexanal share is 11 / 13 / 20 % across them at zero
additive. A single scalar branch fraction is falsified by the fit column itself.

The state therefore carries the hydroperoxide pool resolved along the two axes
Frankel resolved it along:

  * POSITION -- 9-OOH or 13-OOH. This is what decides WHICH products are even
    reachable: pentane and hexanal and methyl 13-oxo-9,11-tridecadienoate come
    from the 13-hydroperoxide; methyl octanoate, 2,4-decadienal and methyl
    9-oxononanoate come from the 9-hydroperoxide. These are STRUCTURAL ZEROS,
    not small numbers.
  * GEOMETRY -- cis,trans or trans,trans. Frankel's central structural finding
    is that the two geometries behave differently (his pure cis,trans-13 arm
    shows a distribution that his trans,trans arms do not), and the corpus
    contains a separate zero-additive column for each.

WHAT IS *NOT* HERE, AND WHY THAT IS THE POINT
---------------------------------------------
  * **NONANAL is present as a species and gets NO edge from any linoleate
    hydroperoxide.** Nonanal is the C9 fragment of the OLEATE double bond. It
    appears in no table, figure or sentence of Frankel 1989, which fed pure
    methyl linoleate hydroperoxides. That absence is a declared HOLD-OUT
    (negative test) and it is honoured STRUCTURALLY: ``LOOH_OL`` is the only
    pool with a nonanal edge, and its branch fraction is ``None`` -- unmeasured
    anywhere in the fit corpus -- so a request for absolute nonanal in a real
    (oleate-bearing) matrix is REFUSED rather than answered with the shipped
    FAST-lane 0.15.
  * **2-pentylfuran and 1-hexanol are NOT species.** No measured branch
    fraction exists for the alkylfuran route, and no aldehyde-reduction step is
    measured anywhere in the corpus. They stay on the engine's
    unrepresented-compound list with sharper reasons.
  * **2-nonenal and methyl 12-oxo-10-dodecenoate are carried as NAMED
    CO-PRODUCTS with no quantitation.** Frankel's introduction names them as the
    Hock partners of hexanal and methyl 9-oxononanoate, but his six-product
    slate does not include them, so no share can be fitted for them. They are
    recorded in ``NAMED_UNQUANTIFIED_COPRODUCTS`` so their absence from the
    slate is legible rather than silent.

UNITS: mmol/L, minutes, Kelvin -- B1's units, unchanged.
"""

from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .species import Species

# ---------------------------------------------------------------------------
# 1. The hydroperoxide pools
# ---------------------------------------------------------------------------
# Carbon counts are for the METHYL ESTER as Frankel prepared it (methyl
# linoleate hydroperoxide, C19H34O4 -> 19 carbons), because every branch
# fraction in the FIT column is measured on the methyl ester and two of the six
# measured products (methyl octanoate, methyl 9-oxononanoate) still carry the
# ester carbon. Using the free-acid count would break the carbon closure by
# exactly one carbon per molecule, silently.

LIPID_PRECURSORS: Tuple[Species, ...] = (
    Species("LOOH_13_ct", "methyl cis,trans-13-hydroperoxy-octadecadienoate",
            19, 0, "reactant", True,
            "Frankel 1989 prepared this one PURE, by soy lipoxygenase oxidation "
            "of linoleic acid followed by silicic acid chromatography and "
            "esterification. Its zero-additive column is one of the three FIT "
            "columns."),
    Species("LOOH_13_tt", "methyl trans,trans-13-hydroperoxy-octadecadienoate",
            19, 0, "reactant", True,
            "Separated by semi-preparative reversed-phase HPLC from the "
            "autoxidation mixture; enters the FIT column only as part of the "
            "trans,trans 9+13 mixture."),
    Species("LOOH_9_ct", "methyl cis,trans-9-hydroperoxy-octadecadienoate",
            19, 0, "reactant", True),
    Species("LOOH_9_tt", "methyl trans,trans-9-hydroperoxy-octadecadienoate",
            19, 0, "reactant", True),
    Species("LOOH_OL", "oleate hydroperoxide pool (8-/9-/10-/11-OOH, lumped)",
            19, 0, "reactant", False,
            "MEASURED BY NOTHING IN THE FIT CORPUS. Carried so that nonanal has "
            "a structurally correct parent rather than an invented one, and so "
            "that a real (oleate-bearing) matrix can be told apart from "
            "Frankel's pure-linoleate feed. Its branch fractions are None and "
            "any request that needs them is refused."),
)


# ---------------------------------------------------------------------------
# 2. The MEASURED product slate -- Frankel 1989's six columns, and nothing else
# ---------------------------------------------------------------------------

LIPID_PRODUCTS: Tuple[Species, ...] = (
    Species("PENTANE", "pentane", 5, 0, "product", True,
            "Homolytic beta-scission product of the 13-alkoxyl radical. A C5 "
            "ALKANE measured by headspace GC after cryo-trapping at -65 C: its "
            "recovery is not commensurable with a C14 dienoate's, which is why "
            "this module does NOT impose the 1:1 stoichiometric pairing with "
            "methyl 13-oxo-9,11-tridecadienoate that the mechanism implies."),
    Species("HEXANAL", "hexanal", 6, 0, "product", True,
            "The module's headline product. Reachable from the 13-hydroperoxide "
            "by BOTH routes (heterolytic Hock cleavage, and homolytic pathway "
            "B), which is why the hydrogen-donor term raises its SHARE while "
            "lowering the total."),
    Species("ME_OCTANOATE", "methyl octanoate", 9, 0, "product", True,
            "Homolytic product of the 9-alkoxyl radical."),
    Species("DECADIENAL", "trans,trans-2,4-decadienal", 10, 0, "product", True,
            "The alpha,beta-UNSATURATED member of the slate, and the test case "
            "for B4's unsaturation penalty. Doubly conjugated; k2 sec. A.8 "
            "records it as the extreme outlier in all three media."),
    Species("ME_9_OXONONANOATE", "methyl 9-oxononanoate", 10, 0, "product", True,
            "Frankel: 'the 2nd most abundant, up to 47 %'. Reachable by both "
            "routes from the 9-hydroperoxide, the mirror of hexanal."),
    Species("ME_13_OXO_TRIDECADIENOATE", "methyl 13-oxo-9,11-tridecadienoate",
            14, 0, "product", True,
            "The C14 half of the 13-OOH homolytic scission. Frankel's own "
            "Discussion warns that the conjugated dienals 'may be more reactive "
            "and less stable than their saturated counterparts under the "
            "conditions used'."),
)


# ---------------------------------------------------------------------------
# 3. Species that exist so that a REFUSAL can be precise
# ---------------------------------------------------------------------------

LIPID_STRUCTURAL: Tuple[Species, ...] = (
    Species("NONANAL", "nonanal", 9, 0, "product", False,
            "STRUCTURAL ZERO from linoleate. Frankel 1989 fed pure methyl "
            "linoleate hydroperoxides and nonanal appears in no table, figure "
            "or sentence -- the declared HOLD-OUT negative test. It has exactly "
            "one incoming edge in this network, from LOOH_OL, whose branch "
            "fraction is unmeasured. So: exactly 0.0 from a linoleate feed, and "
            "REFUSED in an oleate-bearing matrix."),
    Species("LIPID_FRAG_C", "unassigned lipid fragment carbon", 1, 0, "pool", False,
            "B1's FRAG_C discipline, applied here. Every scission has two "
            "halves and Frankel's slate quantifies only some of them (the Hock "
            "partners 2-nonenal and methyl 12-oxo-10-dodecenoate are named in "
            "his introduction and measured in none of his tables). The "
            "unreported carbon is routed here so the carbon balance closes as "
            "an EQUALITY and the size of the unquantified pool is visible."),
)


#: Named in Frankel 1989's introduction as the Hock partners of the two
#: heterolytic products, and quantified in none of his tables. Carried as text
#: so that "the slate has six members" is a legible statement about the
#: MEASUREMENT rather than about the chemistry.
NAMED_UNQUANTIFIED_COPRODUCTS: Mapping[str, str] = {
    "trans-2-nonenal": (
        "the Hock partner of methyl 9-oxononanoate from the 9-hydroperoxide "
        "(Frankel 1989 introduction, attributing refs 4-7). Not in the "
        "six-product slate; no share is fitted for it, and none is invented."
    ),
    "methyl 12-oxo-10-dodecenoate": (
        "the Hock partner of HEXANAL from the 13-hydroperoxide (same source). "
        "Its absence from the slate is why the hexanal share cannot be read as "
        "a mass yield: the slate's denominator is the sum of six MEASURED "
        "peaks, not the scission's full carbon."
    ),
}


LIPID_SPECIES: Tuple[Species, ...] = (
    LIPID_PRECURSORS + LIPID_PRODUCTS + LIPID_STRUCTURAL
)

LIPID_KEYS: Tuple[str, ...] = tuple(s.key for s in LIPID_SPECIES)
LIPID_INDEX: Mapping[str, int] = {s.key: i for i, s in enumerate(LIPID_SPECIES)}
N_LIPID_SPECIES = len(LIPID_SPECIES)

#: The six MEASURED columns, in Frankel's own table order. The fitter's input
#: array is asserted to be exactly (3 systems x these 6), which is the firewall
#: that keeps the tocopherol columns out of the objective.
FRANKEL_SLATE: Tuple[str, ...] = tuple(s.key for s in LIPID_PRODUCTS)

#: position -> the products only that position can make. STRUCTURAL ZEROS.
POSITION_PRODUCTS: Mapping[str, Tuple[str, ...]] = {
    "13": ("PENTANE", "HEXANAL", "ME_13_OXO_TRIDECADIENOATE"),
    "9": ("ME_OCTANOATE", "DECADIENAL", "ME_9_OXONONANOATE"),
}

#: The four linoleate pools, decomposed as (position, geometry).
LOOH_POOLS: Mapping[str, Tuple[str, str]] = {
    "LOOH_13_ct": ("13", "ct"),
    "LOOH_13_tt": ("13", "tt"),
    "LOOH_9_ct": ("9", "ct"),
    "LOOH_9_tt": ("9", "tt"),
}

#: Frankel's Scheme 1 / introduction assignment of each measured product to a
#: cleavage MECHANISM. Taken from the INTRODUCTION, which attributes it to refs
#: 3-10 (all pre-1989) -- NOT from the tocopherol arms, which are hold-out.
#:
#:   "homolytic" -- reachable ONLY by homolytic beta-scission of the alkoxyl
#:                  radical, so a hydrogen donor suppresses it;
#:   "both"      -- reachable by heterolytic (Hock) cleavage as well, so a
#:                  hydrogen donor does not remove it.
CLEAVAGE_MECHANISM: Mapping[str, str] = {
    "PENTANE": "homolytic",
    "ME_13_OXO_TRIDECADIENOATE": "homolytic",
    "ME_OCTANOATE": "homolytic",
    "DECADIENAL": "homolytic",
    "HEXANAL": "both",
    "ME_9_OXONONANOATE": "both",
}


# ---------------------------------------------------------------------------
# 4. Molecular weights -- arithmetic, not parameters, but a place a 1000x hides
# ---------------------------------------------------------------------------

MOLECULAR_WEIGHT_G_PER_MOL: Mapping[str, float] = {
    "LOOH_13_ct": 310.47,   # C19H34O4, methyl linoleate hydroperoxide
    "LOOH_13_tt": 310.47,
    "LOOH_9_ct": 310.47,
    "LOOH_9_tt": 310.47,
    "LOOH_OL": 312.49,      # C19H36O4, methyl oleate hydroperoxide
    "PENTANE": 72.15,       # C5H12
    "HEXANAL": 100.16,      # C6H12O
    "ME_OCTANOATE": 158.24,  # C9H18O2
    "DECADIENAL": 152.23,   # C10H16O
    "ME_9_OXONONANOATE": 186.25,   # C10H18O3
    "ME_13_OXO_TRIDECADIENOATE": 238.32,  # C14H22O3
    "NONANAL": 142.24,      # C9H18O
}


def mmol_per_litre_to_ug_per_litre(key: str, mmol_per_litre: float) -> float:
    """mmol/L -> ug/L. Unit-tested; never inlined at a call site."""
    return float(mmol_per_litre) * MOLECULAR_WEIGHT_G_PER_MOL[key] * 1.0e3


def ug_per_litre_to_mmol_per_litre(key: str, ug_per_litre: float) -> float:
    return float(ug_per_litre) / MOLECULAR_WEIGHT_G_PER_MOL[key] * 1.0e-3


# ---------------------------------------------------------------------------
# 5. Bridge to the B4 output layer
# ---------------------------------------------------------------------------
# B4 keys its structural registry, its thresholds and its unsaturation penalty
# on its OWN compound keys. This is the ONLY mapping between the two, and the
# unit tests assert that every key on the right-hand side exists in B4's
# COMPOUND_STRUCTURE -- so a rename on either side fails loudly instead of
# silently returning NoMeasuredThreshold for a compound that has one (the
# wiring bug the B5 cutover already found once).

B4_COMPOUND_KEY: Mapping[str, str] = {
    "HEXANAL": "hexanal",
    "NONANAL": "nonanal",
    "DECADIENAL": "tt_2_4_decadienal",
}

#: Products with NO B4 structural record, and therefore no threshold, no
#: binding class and no unsaturation gate. Recorded rather than defaulted.
NO_B4_RECORD: Mapping[str, str] = {
    "PENTANE": "an alkane; the B4 registry has no hydrocarbon class at all",
    "ME_OCTANOATE": "a fatty-acid methyl ester; B4's 'ester' class is measured "
                    "on short-chain flavour esters (isoamyl acetate and "
                    "friends) and transferring it across eleven carbons is not "
                    "licensed by any source",
    "ME_9_OXONONANOATE": "an oxo-ester; no class in B4",
    "ME_13_OXO_TRIDECADIENOATE": "an oxo-dienoate ester; no class in B4",
}


def initial_lipid_state(charge: Mapping[str, float]) -> Dict[str, float]:
    """A full zeroed lipid state with ``charge`` (mmol/L) applied."""
    state = {key: 0.0 for key in LIPID_KEYS}
    for key, value in charge.items():
        if key not in state:
            raise KeyError(f"{key!r} is not a lipid species")
        state[key] = float(value)
    return state


def total_lipid_carbon(state: Mapping[str, float]) -> float:
    """Total carbon in mmol/L across the lipid state, INCLUDING the frag pool."""
    return sum(
        float(state.get(s.key, 0.0)) * s.carbon for s in LIPID_SPECIES
    )


__all__ = [
    "B4_COMPOUND_KEY",
    "CLEAVAGE_MECHANISM",
    "FRANKEL_SLATE",
    "LIPID_INDEX",
    "LIPID_KEYS",
    "LIPID_PRODUCTS",
    "LIPID_SPECIES",
    "LOOH_POOLS",
    "MOLECULAR_WEIGHT_G_PER_MOL",
    "NAMED_UNQUANTIFIED_COPRODUCTS",
    "NO_B4_RECORD",
    "N_LIPID_SPECIES",
    "POSITION_PRODUCTS",
    "initial_lipid_state",
    "mmol_per_litre_to_ug_per_litre",
    "total_lipid_carbon",
    "ug_per_litre_to_mmol_per_litre",
]
