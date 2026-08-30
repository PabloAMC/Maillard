"""
src/kinetic_core/species.py

THE STATE VECTOR OF THE MASS-ACTION KINETIC CORE, WITH ITS CARBON AND NITROGEN
BOOKKEEPING (Build Wave B1, 2026-08-28; trunk block extended by B7).
=============================================================================

WAVE: B1 for the thirteen trunk state variables; B7 added the furanic block
(HMF, DMHF and their partners). The sulfur, acrylamide and lipid lanes keep
their own state in ``species_sulfur``, ``species_acrylamide`` and
``species_lipid`` -- one file per lane, so a lane's atom counts and its
molecular weights cannot drift apart.
EXAM: none of its own. Nothing here is fitted and nothing here is a rate, so
there is nothing to score; what this file is checked against is the CONSERVATION
INVARIANT, enforced at import by ``network.validate_balance()``.
DECLARED GAPS: two of the thirteen trunk variables are not molecular
concentrations (see the melanoidin note below), so they have no molecular weight
and cannot be converted to ug/L -- the engine reports them in their own unit
rather than inventing one. A species with no entry in
``MOLECULAR_WEIGHT_G_PER_MOL`` is in that category by construction, not by
omission.

Every entry below carries the atom counts that make the conservation invariant
computable. Nothing in this file is fitted, and nothing in it is a rate.

Two of the thirteen state variables are NOT molecular concentrations:

  * ``MEL_C`` / ``MEL_N``  -- the MELANOIDIN MASS SINK, carried in
    mmol of ELEMENT per litre rather than mmol of "melanoidin molecule" per
    litre. Melanoidins are a polydisperse polymer class with no molecular
    weight, so a molar concentration of them is only meaningful relative to a
    declared repeat unit. The pool is therefore carried elementally, and the
    repeat-unit molarity that Martins' browning readout uses is DERIVED from it
    (see ``melanoidin_repeat_units``), not stored.

  * ``FRAG_C`` -- unassigned fragment carbon. Several MEASURED steps in the
    trunk report only one of their products (Martins measures the formic acid
    from 3-deoxyglucosone but not the C5 residue that leaves with it). The
    unreported carbon is not thrown away; it is routed here so that the total
    carbon balance closes exactly and the size of the unassigned pool is
    visible rather than hidden. FRAG_C is an accounting pool, NOT a chemical
    species, and nothing consumes it.

UNITS: concentrations mmol/L, time minutes, temperature Kelvin, Ea kJ/mol.
These are the units the source data are printed in; nothing is converted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Tuple


@dataclass(frozen=True)
class Species:
    """One state variable, with the atom counts the invariants need."""

    key: str
    label: str
    carbon: int
    nitrogen: int
    #: "reactant" | "intermediate" | "product" | "pool"
    role: str
    #: is this species MEASURED in the fit corpus?
    measured: bool
    note: str = ""
    #: SULFUR atoms per unit. Added by Build Wave B2 (the sulfur module) as a
    #: keyword field with a zero default, so that every B1 trunk entry above --
    #: all of which are sulfur-free -- is unchanged and still constructs
    #: positionally. The sulfur-bearing species live in ``species_sulfur.py``.
    sulfur: int = 0


SPECIES: Tuple[Species, ...] = (
    Species("Glc", "D-glucose", 6, 0, "reactant", True),
    Species("Fru", "D-fructose", 6, 0, "intermediate", True,
            "Lobry de Bruyn-Alberda van Ekenstein partner of glucose; measured"),
    Species("Gly", "glycine (the amine)", 2, 1, "reactant", True,
            "the amine pool. Glycine in the Martins system; any alpha-amino acid "
            "with the same carbon/nitrogen count substitutes without changing the "
            "bookkeeping, but NOT without changing the rates (see the epsilon "
            "amine-specificity note in parameters.py)"),
    Species("SB", "Schiff base / condensation intermediate", 8, 1, "intermediate", False,
            "NOT MEASURED by any experiment in the fit corpus. Carried structurally "
            "so the condensation can be written as two elementary steps; the source "
            "(Martins 2005 T1) REFUSES the split, so the composite is what is "
            "parameterised -- see SCHIFF_AMADORI_SPLIT in parameters.py"),
    Species("AMA", "Amadori compound (DFG, N-(1-deoxy-D-fructos-1-yl)glycine)", 8, 1,
            "intermediate", True, "measured as 'DFG' in the Martins figures"),
    Species("TDG", "3-deoxyglucosone (3-DG)", 6, 0, "intermediate", True),
    Species("ODG", "1-deoxyglucosone (1-DG)", 6, 0, "intermediate", True),
    Species("MGO", "methylglyoxal", 3, 0, "intermediate", True),
    Species("FA", "formic acid", 1, 0, "product", True),
    Species("AA", "acetic acid", 2, 0, "product", True),
    Species("MEL_C", "melanoidin pool, CARBON", 1, 0, "pool", False,
            "mmol of carbon per litre held in the terminal melanoidin polymer"),
    Species("MEL_N", "melanoidin pool, NITROGEN", 0, 1, "pool", False,
            "mmol of nitrogen per litre held in the terminal melanoidin polymer"),
    Species("FRAG_C", "unassigned fragment carbon", 1, 0, "pool", False,
            "carbon leaving a measured step in an unmeasured co-product. An "
            "accounting pool, not a species; nothing consumes it"),
    # -----------------------------------------------------------------------
    # BUILD WAVE B7 -- THE FURANIC CHANNEL. Appended, never interleaved, so
    # that every B1/B2/B3 index above is unchanged and a pre-B7 state vector is
    # still a prefix of this one.
    #
    # These five live on the TRUNK rather than in a lane of their own because
    # their parents already do: HMF's two sources are fructose and
    # 3-deoxyglucosone, and DMHF's are 1-deoxyglucosone and methylglyoxal.
    # All four parents are B1 trunk species, so the furanic block is reachable
    # from the trunk, the sulfur and the acrylamide lanes without any lane
    # composing with any other. See ``furanic.py`` for the topology and
    # ``parameters_furanic.py`` for every constant and its source.
    # -----------------------------------------------------------------------
    Species("INT", "undetermined cyclic intermediate of fructose dehydration",
            6, 0, "intermediate", False,
            "Kocadagli & Gokmen 2016 (JAFC 10.1021/acs.jafc.6b01862) call this "
            "'Int' and say in as many words that it is UNDETERMINED and "
            "UNMEASURED. Its concentration scale is therefore NOT identified: "
            "only the product k7*[Int] is constrained by their data, so "
            "k_fru_int and k_int_hmf are carried as a PAIR and neither is "
            "transferable on its own. This is K5a constraint C2 and it is why "
            "no rate constant on this limb may ever be compared in magnitude "
            "with one on the measured 3-DG limb (K5a MUST-NOT #3)."),
    Species("DDG", "3,4-dideoxyglucosone (3,4-DG)", 6, 0, "intermediate", True,
            "SEMI-QUANTITATED against the 3-DG response factor in both "
            "Kocadagli papers (author-declared, K5a C22), so its absolute "
            "scale carries an unknown multiplicative error and both edges that "
            "touch it inherit it. Carried as a state variable rather than "
            "lumped away because 3-DG -> 3,4-DG is the RATE-DETERMINING STEP "
            "of the 3-DG limb in two independent matrices (K5a C3)."),
    Species("HMF", "5-hydroxymethylfurfural (5-HMF)", 6, 0, "product", True,
            "THE compound of the K5a cluster. NOT norfuraneol: two papers the "
            "repo already holds (whitfield1999, whitfield2001) and two in the "
            "K5b cluster (blank1996 'HMF (3)', apriyantono1993 'HMFone') use "
            "an HMF-shaped token to mean 4-hydroxy-5-methyl-3(2H)-furanone, "
            "which is species ``NF`` on the sulfur lane and a different "
            "molecule. See ``furanic.NAMING_TRAPS``."),
    Species("AF", "acetylformoin (4,5-dihydroxy-2,6-heptanedione, cyclised)",
            6, 0, "intermediate", False,
            "The DMHF progenitor on the INTACT-SKELETON edge, and the species "
            "that STRUCTURALLY SEPARATES the two DMHF routes: Wang & Ho 2008 "
            "fed [13C6]glucose + [12C3]methylglyoxal and observed NO "
            "[12C6]acetylformoin, so the methylglyoxal route does not pass "
            "through it (K5b B7). Unmeasured in magnitude anywhere."),
    Species("DMHF", "2,5-dimethyl-4-hydroxy-3(2H)-furanone (furaneol, HDMF)",
            6, 0, "product", True,
            "Written HDMF by Blank 1996/1997 and Poisson 2019 and DMHF by "
            "Wang & Ho 2008 and Shu & Ho 1988; the repo uses DMHF. A DIFFERENT "
            "COMPOUND from norfuraneol (``NF``, C5): the pre-B7 engine refused "
            "every DMHF request on exactly that ground and the refusal was "
            "correct."),
)

SPECIES_KEYS: Tuple[str, ...] = tuple(s.key for s in SPECIES)
INDEX: Mapping[str, int] = {s.key: i for i, s in enumerate(SPECIES)}
BY_KEY: Mapping[str, Species] = {s.key: s for s in SPECIES}

N_SPECIES = len(SPECIES)

#: species whose atom counts are per-MOLECULE (everything except the elemental pools)
MOLECULAR_KEYS: Tuple[str, ...] = tuple(
    s.key for s in SPECIES if s.role != "pool"
)

#: mmol/L of measured species <-> state key, for the Martins figure labels
MEASURED_LABEL_TO_KEY: Mapping[str, str] = {
    "glucose": "Glc",
    "fructose": "Fru",
    "glycine": "Gly",
    "DFG": "AMA",
    "3-DG": "TDG",
    "1-DG": "ODG",
    "methylglyoxal": "MGO",
    "formic_acid": "FA",
    "acetic_acid": "AA",
    # 'melanoidins' is DELIBERATELY ABSENT. It is the Module 4 hold-out
    # (docs/reference/FIT_HOLDOUT_DECLARATION.md D.6, "Martins 2005 browning,
    # step 9, epsilon 0.64"). The hold-out scorer maps it explicitly, at
    # scoring time, and the fit objective cannot reach it through this table.
}


def carbon_vector() -> Tuple[int, ...]:
    """Carbon atoms per unit of each state variable."""
    return tuple(s.carbon for s in SPECIES)


def nitrogen_vector() -> Tuple[int, ...]:
    """Nitrogen atoms per unit of each state variable."""
    return tuple(s.nitrogen for s in SPECIES)


def sulfur_vector() -> Tuple[int, ...]:
    """
    Sulfur atoms per unit of each state variable.

    Every B1 trunk species is sulfur-free, so this is all zeros for the trunk.
    It exists so that the balance checker can be run over the same three
    elements on the trunk and on the sulfur extension without branching.
    """
    return tuple(s.sulfur for s in SPECIES)


def total_carbon(state) -> float:
    """Total carbon, mmol C/L, summed over every pool including the sinks."""
    return float(sum(c * float(state[i]) for i, c in enumerate(carbon_vector())))


def total_nitrogen(state) -> float:
    """Total nitrogen, mmol N/L, summed over every pool including the sinks."""
    return float(sum(n * float(state[i]) for i, n in enumerate(nitrogen_vector())))


def total_sulfur(state) -> float:
    """Total sulfur, mmol S/L, summed over every pool including the sinks."""
    return float(sum(s * float(state[i]) for i, s in enumerate(sulfur_vector())))


#: Carbon atoms in one melanoidin REPEAT UNIT as Martins' step 9 writes it:
#: 3-deoxyglucosone (C6) + glycine (C2) -> one browning-active unit.
#: Source anchor: Martins & van Boekel 2005 Table 2 step 9, "3-DG + Gly ->
#: melanoidins" (data/lit/extraction_dossiers/k3_final_parameter_inventory.md
#: line 119). This is the ONLY basis on which the elemental pool can be turned
#: back into the molar concentration the browning readout is expressed in.
MELANOIDIN_REPEAT_UNIT_CARBON = 8
MELANOIDIN_REPEAT_UNIT_NITROGEN = 1


def melanoidin_repeat_units(state) -> float:
    """
    Melanoidin concentration in mmol of REPEAT UNITS per litre.

    Martins' browning response is a concentration of melanoidin "molecules",
    obtained as A470 / epsilon. The only definition of a melanoidin molecule
    his scheme supplies is the product of step 9, i.e. one 3-DG plus one
    glycine. The repeat-unit count is therefore the NITROGEN pool: every step-9
    event contributes exactly one nitrogen, and every carbon-only addition to
    the polymer (e.g. a trapped methylglyoxal) grows an existing unit rather
    than creating a new one.

    Returning MEL_N rather than MEL_C/8 is a deliberate choice and it matters:
    with carbon-only additions present the two differ, and the nitrogen count
    is the one that tracks "how many step-9 units exist".
    """
    return float(state[INDEX["MEL_N"]])


def melanoidin_c_over_n(state) -> float:
    """Predicted elemental C/N of the melanoidin pool (NaN before any forms)."""
    n = float(state[INDEX["MEL_N"]])
    if n <= 0.0:
        return float("nan")
    return float(state[INDEX["MEL_C"]]) / n


def initial_state(concentrations: Mapping[str, float]):
    """Build a state vector from a {species_key: mmol/L} mapping."""
    import numpy as np

    y0 = np.zeros(N_SPECIES, dtype=float)
    for key, value in concentrations.items():
        if key not in INDEX:
            raise KeyError(f"unknown species {key!r}; expected one of {SPECIES_KEYS}")
        if float(value) < 0.0:
            raise ValueError(f"negative initial concentration for {key!r}")
        y0[INDEX[key]] = float(value)
    return y0


def state_as_dict(state) -> Dict[str, float]:
    return {key: float(state[INDEX[key]]) for key in SPECIES_KEYS}
