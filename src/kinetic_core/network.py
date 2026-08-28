"""
src/kinetic_core/network.py

THE MASS-ACTION REACTION NETWORK OF THE TRUNK, WITH THE MELANOIDIN MASS SINK.
=============================================================================

FIFTEEN steps over the thirteen state variables of ``species.py``. Every step
is written as an explicit reactant->product stoichiometry and the module
REFUSES TO IMPORT if any of them fails to balance carbon or nitrogen. The
conservation invariant is therefore a property of the network's construction,
not something the tests hope to observe.

WHAT IS DIFFERENT FROM THE SEED (``src/trunk_kinetics.py``)
-----------------------------------------------------------
The seed carries six species, eight steps, two undeclared lumped sinks and no
product pools: fructose, the acids, methylglyoxal and the melanoidins are
absent, and its deoxyosone sinks discard their carbon. This network:

  * carries all ten MEASURED Martins responses as states;
  * gives EVERY intermediate a formation term and a consumption term --
    including methylglyoxal, which Martins' own scheme lets accumulate forever;
  * routes every atom of carbon that leaves a measured step in an unmeasured
    co-product into an explicit ``FRAG_C`` pool instead of deleting it;
  * terminates in an explicit MELANOIDIN MASS SINK carried elementally.

THE MELANOIDIN MASS SINK
------------------------
Two channels feed it, and only two:

  1. ``r_tdg_mel``  3-DG + Gly -> melanoidin. MEASURED: Martins step 9,
     X = 8.12e-4 L/(mmol*min) at 100 C, Ea 95.2 +/- 2.3 kJ/mol. Contributes
     8 carbon and 1 nitrogen per event -- the C6 of the deoxyosone plus the C2
     and the N of the amine, which is the stoichiometry the step's own reaction
     equation states.
  2. ``r_mgo_mel``  methylglyoxal -> melanoidin carbon. FITTED HERE (no
     literature value); contributes 3 carbon and no nitrogen.

Because (2) adds carbon without nitrogen, the pool's predicted C/N RISES with
heating time. That direction is measurable and is checked against Brands' 2002
elemental C/N series as a directional-only diagnostic (see the fit report);
it is not fitted to.

NOTHING LEAVES THE SINK. It is terminal by construction, which is what makes
the total-carbon invariant an equality rather than an inequality.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from .parameters import (
    MARTINS_M4,
    SCHIFF_AMADORI_SPLIT,
    KineticParameter,
)
from .species import (
    BY_KEY,
    INDEX,
    N_SPECIES,
    SPECIES_KEYS,
)


@dataclass(frozen=True)
class Reaction:
    """One elementary step, with an explicit balanced stoichiometry."""

    key: str
    #: {species_key: stoichiometric coefficient}, positive integers
    reactants: Mapping[str, int]
    products: Mapping[str, int]
    #: key into the parameter registry, or None for a derived rate
    parameter_key: Optional[str]
    note: str = ""

    @property
    def order(self) -> int:
        return int(sum(self.reactants.values()))

    def atom_balance(
        self, element: str, lookup: Optional[Mapping[str, object]] = None
    ) -> Tuple[int, int]:
        """
        (left, right) atom counts for 'carbon', 'nitrogen' or 'sulfur'.

        ``lookup`` defaults to the trunk's species table. Build Wave B2 passes
        the EXTENDED table (trunk + sulfur species) so that the same Reaction
        type and the same balance arithmetic serve both networks; nothing about
        the trunk's own behaviour changes when the argument is omitted.
        """
        table = BY_KEY if lookup is None else lookup

        def side(mapping: Mapping[str, int]) -> int:
            return sum(
                coefficient * getattr(table[key], element)
                for key, coefficient in mapping.items()
            )

        return side(self.reactants), side(self.products)


# ---------------------------------------------------------------------------
# The network
# ---------------------------------------------------------------------------
# Carbon accounting note, once, because it is the whole point of the FRAG_C
# pool: Martins measures ONE product of several of these steps and not the
# rest. Step 5 reports the formic acid from a C6 deoxyosone; the other five
# carbons exist and are unmeasured. Writing "3-DG -> formic acid" as if five
# carbons vanished is precisely the defect this rebuild removes, so the residue
# is routed to FRAG_C and its size is reported.

REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_schiff", {"Glc": 1, "Gly": 1}, {"SB": 1}, "k_schiff",
        "Martins step 1, first half. The source refuses the Schiff/Amadori "
        "split; the pair is a composite (see parameters.SCHIFF_AMADORI_SPLIT).",
    ),
    Reaction(
        "r_amadori", {"SB": 1}, {"AMA": 1}, None,
        "Martins step 1, second half. Rate DERIVED from r_schiff by the pinned "
        "44.9x ratio, not independently parameterised. Irreversible and "
        "sink-free, so every condensed molecule reaches the Amadori pool and "
        "the pair reduces exactly to Martins' one-step step 1.",
    ),
    Reaction("r_glc_fru", {"Glc": 1}, {"Fru": 1}, "k_glc_fru", "Martins step 2"),
    Reaction("r_fru_glc", {"Fru": 1}, {"Glc": 1}, "k_fru_glc", "Martins step 3"),
    Reaction(
        "r_ama_tdg", {"AMA": 1}, {"TDG": 1, "Gly": 1}, "k_ama_tdg",
        "Martins step 4. Regenerates the amine.",
    ),
    Reaction(
        "r_tdg_fa", {"TDG": 1}, {"FA": 1, "FRAG_C": 5}, "k_tdg_fa",
        "Martins step 5. The C5 residue is unmeasured and is routed to FRAG_C.",
    ),
    Reaction(
        "r_ama_mgo", {"AMA": 1}, {"MGO": 1, "Gly": 1, "FRAG_C": 3}, "k_ama_mgo",
        "Martins step 6. Releases the amine in Martins' scheme; the C3 residue "
        "of the sugar skeleton is unmeasured and is routed to FRAG_C.",
    ),
    Reaction("r_ama_odg", {"AMA": 1}, {"ODG": 1, "Gly": 1}, "k_ama_odg", "Martins step 7"),
    Reaction(
        "r_odg_aa", {"ODG": 1}, {"AA": 1, "FRAG_C": 4}, "k_odg_aa",
        "Martins step 8. C4 residue unmeasured -> FRAG_C.",
    ),
    Reaction(
        "r_tdg_mel", {"TDG": 1, "Gly": 1}, {"MEL_C": 8, "MEL_N": 1}, "k_tdg_mel",
        "Martins step 9 -- THE MEASURED MELANOIDIN MASS SINK. One deoxyosone "
        "(C6) plus one glycine (C2, N1) enter the terminal polymer.",
    ),
    Reaction(
        "r_fru_acids", {"Fru": 1}, {"FA": 1, "AA": 1, "FRAG_C": 3}, "k_fru_acids",
        "Martins step 10. C3 residue unmeasured -> FRAG_C.",
    ),
    # ---- the fitted extension: consumption terms the corpus has no rate for --
    Reaction(
        "r_glc_frag", {"Glc": 1}, {"FRAG_C": 6}, "k_glc_frag",
        "FITTED. The amine-independent sugar lane; its existence is measured "
        "(Martins' formic:acetic inversion without amine) but its rate is not.",
    ),
    Reaction(
        "r_mgo_mel", {"MGO": 1}, {"MEL_C": 3}, "k_mgo_mel",
        "FITTED. Methylglyoxal into the melanoidin pool -- the consumption term "
        "Martins' scheme omits entirely. Carbon-only, so it RAISES the pool's "
        "predicted C/N.",
    ),
    Reaction("r_fa_frag", {"FA": 1}, {"FRAG_C": 1}, "k_fa_frag", "FITTED."),
    Reaction("r_aa_frag", {"AA": 1}, {"FRAG_C": 2}, "k_aa_frag", "FITTED."),
)

REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in REACTIONS)


# ---------------------------------------------------------------------------
# Construction-time validation
# ---------------------------------------------------------------------------


#: The elements the balance checker enforces. SULFUR was added by Build Wave B2
#: (the sulfur module). Every B1 trunk species carries ``sulfur = 0``, so the
#: extra element is a no-op on the trunk and a real constraint on the sulfur
#: network in ``sulfur.py``, which calls this same function.
BALANCED_ELEMENTS: Tuple[str, ...] = ("carbon", "nitrogen", "sulfur")


def validate_balance(reactions: Sequence[Reaction] = REACTIONS) -> None:
    """Raise unless every reaction balances carbon, nitrogen AND sulfur."""
    for reaction in reactions:
        for element in BALANCED_ELEMENTS:
            left, right = reaction.atom_balance(element)
            if left != right:
                raise ValueError(
                    f"{reaction.key}: {element} does not balance "
                    f"({left} -> {right}). Every step must conserve atoms; the "
                    f"unmeasured residue belongs in FRAG_C, not in nowhere."
                )
        for key in list(reaction.reactants) + list(reaction.products):
            if key not in INDEX:
                raise ValueError(f"{reaction.key}: unknown species {key!r}")
        if BY_KEY.get("FRAG_C") and "FRAG_C" in reaction.reactants:
            raise ValueError(
                f"{reaction.key}: FRAG_C is an accounting pool and must never be "
                f"a reactant."
            )
        if "MEL_C" in reaction.reactants or "MEL_N" in reaction.reactants:
            raise ValueError(
                f"{reaction.key}: the melanoidin sink is TERMINAL; nothing may "
                f"consume it."
            )


validate_balance()


def stoichiometric_matrix(
    reactions: Sequence[Reaction] = REACTIONS,
) -> np.ndarray:
    """(n_species, n_reactions) net stoichiometry."""
    matrix = np.zeros((N_SPECIES, len(reactions)), dtype=float)
    for j, reaction in enumerate(reactions):
        for key, coefficient in reaction.reactants.items():
            matrix[INDEX[key], j] -= float(coefficient)
        for key, coefficient in reaction.products.items():
            matrix[INDEX[key], j] += float(coefficient)
    return matrix


STOICHIOMETRY: np.ndarray = stoichiometric_matrix()


# ---------------------------------------------------------------------------
# Rate evaluation
# ---------------------------------------------------------------------------


def rate_constants_at(
    parameters: Mapping[str, KineticParameter], temperature_k: float
) -> Dict[str, float]:
    """
    Evaluate every reaction's rate constant at ``temperature_k``.

    ``parameters`` must contain every ``parameter_key`` the network references.
    The Amadori rearrangement has no parameter of its own: its constant is
    DERIVED from the condensation by the pinned split ratio, and this is the
    only derived rate in the module.
    """
    out: Dict[str, float] = {}
    for reaction in REACTIONS:
        if reaction.parameter_key is None:
            continue
        parameter = parameters.get(reaction.parameter_key)
        if parameter is None:
            raise KeyError(
                f"{reaction.key}: no parameter {reaction.parameter_key!r} supplied"
            )
        if parameter.k_ref is None or parameter.ea_kj_mol is None:
            raise ValueError(
                f"{reaction.key}: parameter {parameter.key!r} is unpopulated "
                f"(evidence_class={parameter.evidence_class}). The fitted steps "
                f"must be given values before the network can be integrated; "
                f"there is no silent default."
            )
        out[reaction.key] = parameter.k_at(temperature_k)

    schiff = out["r_schiff"]  # L/(mmol*min)
    out["r_amadori"] = (
        float(SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
        * schiff
        * float(SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"])
    )
    return out


def reaction_rates(state: np.ndarray, k: Mapping[str, float]) -> np.ndarray:
    """Mass-action rate of every reaction, in mmol/(L*min)."""
    y = np.clip(np.asarray(state, dtype=float), 0.0, None)
    rates = np.empty(len(REACTIONS), dtype=float)
    for j, reaction in enumerate(REACTIONS):
        value = k[reaction.key]
        for key, coefficient in reaction.reactants.items():
            value *= y[INDEX[key]] ** coefficient
        rates[j] = value
    return rates


def derivatives(state: np.ndarray, k: Mapping[str, float]) -> np.ndarray:
    """d(state)/dt, mmol/(L*min)."""
    return STOICHIOMETRY @ reaction_rates(state, k)


def describe() -> Dict[str, object]:
    """A machine-readable description of the network, for the reports."""
    return {
        "species": [
            {
                "key": s.key,
                "label": BY_KEY[s.key].label,
                "carbon": BY_KEY[s.key].carbon,
                "nitrogen": BY_KEY[s.key].nitrogen,
                "role": BY_KEY[s.key].role,
                "measured_in_fit_corpus": BY_KEY[s.key].measured,
            }
            for s in (BY_KEY[key] for key in SPECIES_KEYS)
        ],
        "reactions": [
            {
                "key": r.key,
                "equation": " + ".join(
                    f"{c if c > 1 else ''}{s}" for s, c in r.reactants.items()
                )
                + " -> "
                + " + ".join(f"{c if c > 1 else ''}{s}" for s, c in r.products.items()),
                "order": r.order,
                "parameter_key": r.parameter_key,
                "carbon_balance": list(r.atom_balance("carbon")),
                "nitrogen_balance": list(r.atom_balance("nitrogen")),
                "note": r.note,
            }
            for r in REACTIONS
        ],
        "measured_parameter_keys": sorted(MARTINS_M4),
        "melanoidin_sink_channels": ["r_tdg_mel (measured)", "r_mgo_mel (fitted here)"],
    }
