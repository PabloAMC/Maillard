"""
src/kinetic_core/species_acrylamide.py

THE ACRYLAMIDE / SAFETY EXTENSION OF THE STATE VECTOR (Build Wave B3).
=============================================================================

Build Wave B1 declared thirteen sulfur-free trunk species in ``species.py``;
Build Wave B2 appended twenty-nine sulfur species in ``species_sulfur.py``.
This module APPENDS the acrylamide block and publishes the CONCATENATED table,
in exactly the same way and for exactly the same reason: the first
``len(SULFUR_STATE)`` entries of ``ACRYLAMIDE_STATE`` are B1's followed by
B2's, in their order, so a trunk state vector is a prefix of a sulfur state
vector is a prefix of an acrylamide state vector, and every index either wave
published is still valid.

WHY THE ACRYLAMIDE BLOCK SITS ON TOP OF THE SULFUR BLOCK AND NOT ON THE TRUNK
----------------------------------------------------------------------------
Because of ONE species: ``Cys``. De Vleeschouwer 2009 Part II measures the
bimolecular acrylamide + cysteine adduct constant ``k_E2`` -- at 2.4 % relative
standard error the tightest parameter in the whole corpus -- and that reaction
consumes the SAME cysteine pool that B2's thermolysis, Nedvidek and
2-acetylthiazole channels consume. Declaring a second, private ``Cys`` here
would let the two modules each spend the same cysteine, which is precisely the
double-counting the rebuild exists to remove. So the acrylamide block shares
B2's cysteine species, by index, and the adduct it forms carries a sulfur atom
that the three-element balance checker polices.

WHAT IS NEW HERE
----------------
Nine species: the asparagine lane (Asn, its Schiff base, the Amadori-analog
intermediate, aspartic acid), acrylamide itself, the measured cysteine adduct,
and the three COMPETING AMINO ACIDS of Claeys 2005's competitor systems
(glutamine, lysine, alanine -- cysteine, the fourth, already exists).

The competitors are carried as REAL SPECIES with their real atom counts rather
than as one generic "competitor" pseudo-species, because they do not have the
same formula and because the whole point of the competition term is that it is
a mass-action channel over a SHARED PRECURSOR POOL -- glucose -- and not a
per-amino-acid multiplier on an acrylamide yield. A lumped competitor would
have had to be given an invented atom count, and the balance checker would then
have been policing an invention.

WHAT IS *NOT* HERE
------------------
No sugar other than glucose and fructose (both are B1's). There is no sucrose
species, which is why the declared sucrose HOLD-OUT is reported UNSCOREABLE
rather than given a number: sucrose hydrolysis has a measured constant in the
corpus but it sits in the HOLD-OUT column, so the module may not carry it.

UNITS: mmol/L, minutes, Kelvin, kJ/mol -- B1's and B2's units, unchanged.
"""

from __future__ import annotations

from typing import Dict, Mapping, Tuple

from .species import Species
from .species_sulfur import SULFUR_STATE

# ---------------------------------------------------------------------------
# The acrylamide block
# ---------------------------------------------------------------------------

ACRYLAMIDE_SPECIES: Tuple[Species, ...] = (
    # ---- the asparagine lane ----------------------------------------------
    Species("Asn", "L-asparagine", 4, 2, "reactant", True,
            "THE acrylamide precursor. Its amide nitrogen is the one that ends "
            "up in the acrylamide; its alpha-nitrogen leaves as ammonia on "
            "deamidation to aspartic acid, which is the competing fate De "
            "Vleeschouwer 2009 I measures separately (k_Asp, 26.43e-3 min^-1 at "
            "160 C, Ea 105.4 kJ/mol -- 'remarkably sugar-independent', 105-109 "
            "across all three sugars)."),
    Species("SBA", "Schiff base of asparagine + glucose", 10, 2, "intermediate", False,
            "NOT MEASURED by any experiment in the fit corpus, and the source "
            "REFUSES the split: De Vleeschouwer's Scheme 4 lumps everything "
            "between the reactants and the acrylamide-forming step into one "
            "intermediate, 'Int1'. Carried as a state for the same reason B1 "
            "carries the glycine Schiff base -- so a later module has something "
            "to attach to -- and parameterised as a composite whose "
            "rate-determining constant is the measured bimolecular one. See "
            "ASN_SCHIFF_AMADORI_SPLIT in parameters_acrylamide.py."),
    Species("INT1", "Amadori-analog Maillard intermediate (De Vleeschouwer's Int1)",
            10, 2, "intermediate", False,
            "The branch point of the whole module: it partitions between "
            "acrylamide (measured, k_F) and everything else (NOT measured -- De "
            "Vleeschouwer's k_M is one of the parameters its own authors mark "
            "'NO PHYSICAL MEANING', so the competing sink is FITTED here). That "
            "partition is what sets the acrylamide yield, and it is the single "
            "most consequential fitted number in this wave."),
    Species("Asp", "L-aspartic acid", 4, 1, "intermediate", True,
            "The deamidation product. It has a consumption term (k_X, "
            "De Vleeschouwer II's cysteine column, 0.04 min^-1, Ea 97.2) -- the "
            "only k_X in the corpus that is neither 'Indeterminate' nor "
            "physically absurd (Part I's glucose value is indeterminate and its "
            "Ea_X is 668.9 kJ/mol)."),
    Species("ACR", "acrylamide", 3, 1, "product", True,
            "odour-irrelevant and toxicologically central. Formation AND "
            "elimination are both explicit steps; see the module docstring of "
            "acrylamide.py for why an acrylamide pool without a real "
            "elimination term is the specific defect this wave removes."),
    # ---- the measured scavenging adduct ------------------------------------
    Species("ACRCYS", "S-(2-carbamoylethyl)cysteine (the acrylamide-cysteine adduct)",
            6, 2, "pool", False,
            "TERMINAL. The Michael adduct of acrylamide onto the cysteine "
            "thiol. It is the only acrylamide sink in the corpus with a "
            "measured bimolecular constant (k_E2 = 49.36 +/- 1.18 M^-1 min^-1 "
            "at 160 C, Ea 51.3 +/- 1.5 kJ/mol) and the only one whose product "
            "is identified, which is why it is a named species while the fitted "
            "amine adducts are routed to the fragment pools. Nothing in the "
            "corpus measures its further reaction, so it is terminal by "
            "declaration rather than by chemistry, and the flag says so.",
            sulfur=1),
    # ---- the Claeys 2005 competitor amino acids ----------------------------
    # Cysteine, the fourth competitor, is B2's `Cys` and is NOT redeclared.
    Species("Gln", "L-glutamine (Claeys competitor)", 5, 2, "reactant", True,
            "The competitor that PROMOTES acrylamide (Claeys 2005 T2: k_F "
            "1.640e-3 against the control's 0.451e-3). This module has NO "
            "promotion mechanism and does not acquire one -- see the "
            "'deliberate under-fit' block in parameters_acrylamide.py. The "
            "reason is that glutamine's promotion carries the B5.5 "
            "sign-crossing whose other half is a declared HOLD-OUT."),
    Species("Lys", "L-lysine (Claeys competitor)", 6, 2, "reactant", True,
            "The competitor with the largest measured effect on acrylamide "
            "ELIMINATION (Claeys k_E 280.2e-3 against 111.1e-3). Its "
            "epsilon-amino group is a Michael nucleophile, i.e. the same "
            "family the k_E2 licence names."),
    Species("Ala", "L-alanine (Claeys competitor)", 3, 1, "reactant", True,
            "The NULL competitor: Claeys' alanine row is statistically "
            "indistinguishable from the control on all four constants (same "
            "significance letter on k_F and k_E). It is in the panel precisely "
            "so that the fitted competition channels have a row that must come "
            "out at ~zero, and the bounds allow zero."),
)

#: The concatenated table. B1's entries, then B2's, then B3's -- so every index
#: either earlier wave published is preserved exactly.
ACRYLAMIDE_STATE: Tuple[Species, ...] = SULFUR_STATE + ACRYLAMIDE_SPECIES

ACRYLAMIDE_STATE_KEYS: Tuple[str, ...] = tuple(s.key for s in ACRYLAMIDE_STATE)
ACRYLAMIDE_INDEX: Mapping[str, int] = {
    s.key: i for i, s in enumerate(ACRYLAMIDE_STATE)
}
ACRYLAMIDE_BY_KEY: Mapping[str, Species] = {s.key: s for s in ACRYLAMIDE_STATE}
N_ACRYLAMIDE_STATE = len(ACRYLAMIDE_STATE)

#: The competitor amino acids of the Claeys 2005 competitor systems, in the
#: order the source's Table 2 prints them. `Cys` is B2's species.
COMPETITOR_KEYS: Tuple[str, ...] = ("Gln", "Cys", "Lys", "Ala")

#: Terminal in THIS block: nothing may consume them.
ACRYLAMIDE_TERMINAL_POOLS: Tuple[str, ...] = ("ACRCYS",)


# ---------------------------------------------------------------------------
# Molecular weights -- arithmetic, not parameters
# ---------------------------------------------------------------------------
# Used ONLY to convert a printed ug/kg (ppb) anchor into the module's mmol/L and
# back. Acrylamide's is the one that matters: every safety number in the
# repository is printed in ppb and every rate in this module runs in mmol/L, so
# this is exactly the place a factor of 1000 can hide. It is unit-tested by a
# round trip and by an independent hand-checked value.
ACRYLAMIDE_MW_G_PER_MOL: Mapping[str, float] = {
    "Asn": 132.12,    # C4H8N2O3
    "Asp": 133.10,    # C4H7NO4
    "SBA": 276.24,    # C10H16N2O8, Asn + Glc - H2O
    "INT1": 276.24,
    "ACR": 71.08,     # C3H5NO
    "ACRCYS": 192.24,  # C6H12N2O3S, acrylamide + cysteine
    "Gln": 146.14,    # C5H10N2O3
    "Lys": 146.19,    # C6H14N2O2
    "Ala": 89.09,     # C3H7NO2
}


def acrylamide_ppb(mmol_per_litre: float) -> float:
    """
    mmol/L of acrylamide -> ppb (ug/kg, taken as ug/L at unit density).

    The density assumption is stated rather than hidden: every acrylamide
    reference value in this repository is a ug/kg of a food, and every model
    system in the fit corpus is either an aqueous solution or a powder whose
    printed molarity already embeds its own basis. Converting between them is
    an arithmetic identity at unit density and an approximation otherwise; the
    fit report carries the flag.
    """
    return float(mmol_per_litre) * ACRYLAMIDE_MW_G_PER_MOL["ACR"] * 1.0e3


def ppb_to_mmol_per_litre(ppb: float) -> float:
    """The inverse of ``acrylamide_ppb``."""
    return float(ppb) / ACRYLAMIDE_MW_G_PER_MOL["ACR"] * 1.0e-3


# ---------------------------------------------------------------------------
# Element vectors and invariants over the EXTENDED state
# ---------------------------------------------------------------------------


def acrylamide_element_vector(element: str) -> Tuple[int, ...]:
    return tuple(getattr(s, element) for s in ACRYLAMIDE_STATE)


def total_element_acrylamide(state, element: str) -> float:
    return float(
        sum(
            n * float(state[i])
            for i, n in enumerate(acrylamide_element_vector(element))
        )
    )


def initial_acrylamide_state(concentrations: Mapping[str, float]):
    """Build an extended state vector from a {species_key: mmol/L} mapping."""
    import numpy as np

    y0 = np.zeros(N_ACRYLAMIDE_STATE, dtype=float)
    for key, value in concentrations.items():
        if key not in ACRYLAMIDE_INDEX:
            raise KeyError(
                f"unknown species {key!r}; expected one of {ACRYLAMIDE_STATE_KEYS}"
            )
        if float(value) < 0.0:
            raise ValueError(f"negative initial concentration for {key!r}")
        y0[ACRYLAMIDE_INDEX[key]] = float(value)
    return y0


def acrylamide_state_as_dict(state) -> Dict[str, float]:
    return {
        key: float(state[ACRYLAMIDE_INDEX[key]]) for key in ACRYLAMIDE_STATE_KEYS
    }
