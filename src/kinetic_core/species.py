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
    # ---- Build Wave B13 (2026-09-07): the dicarbonyl trio, TRUNK-ONLY ----------
    # Kocadagli & Gokmen 2016 JAFC (amine-free glucose glass, 160-200 C) measure
    # Glc -> glucosone -> glyoxal and 1-DG -> diacetyl with barriers. Appended at
    # the END so every existing index is unchanged; the sulfur state vector
    # (SPECIES + SULFUR_SPECIES) gains three inert entries and the sulfur network
    # keeps B9's topology (network.DICARBONYL_REACTIONS runs on the trunk only).
    Species("G", "glucosone (D-arabino-hexos-2-ulose)", 6, 0, "intermediate", False,
            "B13. The oxidative entry: glucose -> glucosone -> glyoxal (Kocadagli "
            "2016 JAFC steps 9-10). Not measured in any fit row of this repository."),
    Species("GO", "glyoxal", 2, 0, "product", False,
            "B13. The CML precursor (glyoxal + lysine). Sink: Kocadagli step 15, "
            "barrier FIXED to zero by the authors."),
    Species("DA", "2,3-butanedione (diacetyl)", 4, 0, "product", False,
            "B13. From 1-deoxyglucosone (Kocadagli step 12, Ea 150.8); the "
            "3-mercapto-2-butanone precursor once a sulfur wave adopts it. Sink: "
            "Kocadagli step 17, rate 0 in the source."),
    # ---- Build Wave B18 (2026-09-08): the pyrazine step, TRUNK-ONLY -----------------
    # Zhou 2024 (JAFC 72:18630) measured pyrazine and 2,5-dimethylpyrazine formation from
    # fed glyoxal / methylglyoxal + alanine at 100-120 C; Leahy & Reineccius 1989 the pH
    # ladder. Two aminoketone intermediates (the Strecker products) and three pyrazines,
    # appended at the END so every existing index is unchanged; the sulfur and acrylamide
    # state vectors leave them out (TRUNK_ONLY_KEYS). Pre-registered in
    # results/validation/kinetic_core_b18_prereg.md; constants in parameters_pyrazine.py.
    Species("PZ", "pyrazine", 4, 2, "product", False,
            "B18. Two aminoacetaldehydes condense (rule R28). The ring carbons are the "
            "dicarbonyl's; the nitrogens the amino acid's (Zhou 2024's isotope labelling, "
            "zhou2024_extraction.md Table 1)."),
    Species("DMP", "2,5-dimethylpyrazine", 6, 2, "product", False,
            "B18. Two aminoacetones condense; the panel's roasted marker."),
    Species("MPZ", "2-methylpyrazine", 5, 2, "product", False,
            "B18. Aminoacetaldehyde + aminoacetone, the mixed condensation. Key MPZ because "
            "MP is the sulfur lane's 1-mercapto-2-propanone."),
    Species("AKG", "aminoacetaldehyde (glyoxal's Strecker aminoketone)", 2, 1, "intermediate", False,
            "B18. Glyoxal + glycine -> aminoacetaldehyde + CO2 + formaldehyde (Strecker, rule "
            "R07); the rate-determining step. Never measured; a steady-state intermediate."),
    Species("AKM", "aminoacetone (methylglyoxal's Strecker aminoketone)", 3, 1, "intermediate", False,
            "B18. Methylglyoxal + glycine -> aminoacetone + CO2 + formaldehyde (Strecker, rule "
            "R07). As AKG."),
    # ---- B20 (2026-09-09): THE GLYCATION ARM, trunk-only. Protein-bound lysine as a reactant:
    # the matrix layer's amine pool becomes a species, the sugar glycates it to the bound Amadori
    # compound, which oxidises to CML, goes to CEL, or decays back to the sugar path (3-DG) and
    # returns the lysine. Rates from Nguyen 2016 (casein + glucose in water, 120 / 130 C), barriers
    # from Berk 2021 and the trunk's own Amadori steps. Pre-registered in
    # results/validation/kinetic_core_b20_prereg.md; constants in parameters_glycation.py.
    Species("LYSP", "protein-bound lysine residue (epsilon-amine site; counted as lysine)", 6, 2, "reactant", True,
            "B20. Charged from the spec's protein loading and the matrix's amine density times the "
            "declared available fraction; zero without a loading, so every earlier pot is unchanged."),
    Species("FLP", "fructosyl-lysine, protein-bound (the bound Amadori compound; furosine's parent)", 12, 2,
            "intermediate", True, "B20. Nguyen 2016's AP; measured as furosine after acid hydrolysis."),
    Species("CML", "N-epsilon-(carboxymethyl)lysine (CML), protein-bound", 8, 2, "product", True,
            "B20. From the bound Amadori compound's oxidative cleavage (Nguyen 2016 k7, Berk 2021 k8); "
            "the glyoxal route fits to zero in both laboratories and is not written."),
    Species("CEL", "N-epsilon-(carboxyethyl)lysine (CEL), protein-bound", 9, 2, "product", True,
            "B20. From the bound Amadori compound via methylglyoxal, lumped as Nguyen 2016 fitted it (k9)."),

    # ---- B22 (2026-09-09): THE METHIONINE CHAIN, trunk-only. Methionine as the Strecker substrate on
    # glyoxal and methylglyoxal (the aminoketones are glycine's, the aldehyde is methionine's), the
    # aldehyde's retro-Michael release of methanethiol, and the disulfide on an apparent constant.
    # Pre-registered in results/validation/kinetic_core_b22_prereg.md; constants in parameters_methionine.py.
    Species("MET", "L-methionine (the Strecker substrate)", 5, 1, "reactant", True,
            "B22. Charged from the spec; also charged, declared, as glycine at the same molarity for the "
            "Amadori chemistry that makes the dicarbonyls.", sulfur=1),
    Species("MTAL", "methional (3-(methylthio)propanal)", 4, 0, "product", True,
            "B22. The Strecker aldehyde of methionine; Pan 2025's zero-order rates at 100-140 C.", sulfur=1),
    Species("MSH", "methanethiol made from methional (the sugar path's pool; the sulfur lane's MESH is thiamine's)", 1, 0,
            "product", True, "B22. The retro-Michael release from methional; acrolein to the fragment pool.", sulfur=1),
    Species("DMDS", "dimethyl disulfide", 2, 0, "product", True,
            "B22. Two methanethiols on an APPARENT constant: the sugar path tracks no oxidant.", sulfur=2),
    # ---- B24 (2026-09-09): 2-ACETYL-1-PYRROLINE FROM PROLINE, trunk-only. The Strecker of proline on
    # methylglyoxal gives 1-pyrroline (the ring nitrogen stays), and methylglyoxal acylates it.
    # Pre-registered in results/validation/kinetic_core_b24_prereg.md; constants in parameters_proline.py.
    Species("PRO", "L-proline (the pyrroline source)", 5, 1, "reactant", True,
            "B24. Charged from the spec; also charged, declared, as glycine for the Amadori chemistry."),
    Species("PYRL", "1-pyrroline", 4, 1, "intermediate", True,
            "B24. Proline's Strecker product; fed by Hofmann & Schieberle 1998b."),
    Species("AP", "2-acetyl-1-pyrroline", 6, 1, "product", True,
            "B24. The acylation of 1-pyrroline by methylglyoxal (Hofmann & Schieberle 1998b, Table 7)."),
)
SPECIES_KEYS: Tuple[str, ...] = tuple(s.key for s in SPECIES)
#: B13: species whose steps exist on the trunk integrator only. The sulfur and acrylamide
#: state vectors leave them out, so those lanes keep the shape their fits were run on.
TRUNK_ONLY_KEYS: Tuple[str, ...] = ("G", "GO", "DA", "PZ", "DMP", "MPZ", "AKG", "AKM", "LYSP", "FLP", "CML", "CEL",
                                     "MET", "MTAL", "MSH", "DMDS", "PRO", "PYRL", "AP")
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

#: THE STRUCTURE ABOVE IS FALSIFIED AND THE ANSWER IT PRODUCES IS NOT (2026-09-09,
#: mundt2004_extraction.md, read in the reading audit).
#:
#: Six carbons from 3-deoxyglucosone plus two from an INTACT glycine set a structural FLOOR
#: of C/N = 8.0 on this pool. Mundt & Wedzicha 2004 measure 7.64 +/- 0.21 on a dialysed
#: glucose-glycine polymer with no protein in it (Table 1, MW > 12 500, n = 4), by two
#: independent methods that agree: CHN microanalysis, and a 14C reconstruction giving
#: whole glycine : DECARBOXYLATED glycine : glucose = 0.289 : 0.662 : 1, i.e. C/N 7.61.
#: The floor is missed by about 1.7 analytical standard deviations, and the radiochemistry
#: says exactly why: about two thirds of the incorporated glycine arrives decarboxylated
#: and contributes ONE carbon per nitrogen, not two.
#:
#: What is falsified is the repeat unit's STRUCTURE, not the number the trunk reports. The
#: measurement is at 70 C and pH 5.5, one point with no series of any kind, and its authors
#: say (citing others, not measuring it) that amino-acid incorporation falls -- and C/N
#: rises -- as temperature rises. So 7.64 is a LOWER BOUND for a 120 C polymer and the
#: trunk's 8.42 to 9.94 clears it. This is a bracket, not a match.
#:
#: Nothing is changed here. A repeat unit that mixed decarboxylated and intact glycine
#: would need a second nitrogen-bearing pool and a branching ratio, neither of which any
#: source on disk measures at cooking temperature; that is a wave, and it is in the backlog.
MELANOIDIN_REPEAT_UNIT_FALSIFYING_MEASUREMENT = (
    "Mundt & Wedzicha 2004 (J. Agric. Food Chem.), glucose 0.25 M + glycine 0.25 M, 0.2 M "
    "acetate pH 5.5, 70.0 C, dialysed MW > 12 500: C/N 7.64 +/- 0.21 by microanalysis and "
    "7.61 by 14C reconstruction, against this unit's structural floor of 8.0."
)

#: CORRECTED 2026-09-09, THE SAME DAY, AND THE CORRECTION MATTERS MORE THAN THE ORIGINAL.
#:
#: The note above was written reading the 7.64 as a LOWER BOUND for a 120 C polymer, on its
#: authors' stated temperature direction, and concluding that the trunk's 8.42 to 9.94 clears
#: it. That conclusion was reached without checking the corpus for a measurement at cooking
#: temperature in the same system. There is one, it was already on this disk and already
#: dossiered, and it says the opposite: the trunk is LOW, not comfortably above a floor.
#:
#: Martins & van Boekel 2003 (Food Chem. 83:135, doi 10.1016/S0308-8146(03)00219-X;
#: martins2003c_extraction.md sections 7 and 8), glucose + glycine, MEASURED microanalysis:
#:
#:      120 C, pH 6.8   C/N = 11, flat over 15 to 60 min
#:      100 C, pH 6.8   C/N = 15, 14, 11, 11 over 30 to 180 min
#:      100 C, pH 5.5   C/N = 19 at 60 min, 16 at 180 min
#:
#: and its Table 2 compiles nine more literature values from 7 to 13, including Cammerer &
#: Kroh 1995's own pair for this system: C/N 7 at 60 C and 9 at 100 C (cammerer1994_extraction.md
#: is that paper, also on disk). The authors' verdict on their own compilation is that the
#: literature values "are not consistent, either with pH or temperature".
#:
#: SO THE HONEST STATEMENT IS THIS. At 70 C one laboratory measures below the structural floor
#: and falsifies the repeat unit; at 100 to 120 C, in the same sugar and the same amine, the
#: nearest measurements sit at 11 to 19 while the trunk predicts 8.42 to 9.94. The model is not
#: bracketed above a floor -- it is between two measurements that disagree with each other by
#: more than it disagrees with either. Nothing is changed here on that basis: a repeat unit
#: that mixed decarboxylated and intact glycine needs a second nitrogen pool and a branching
#: ratio, and the spread across these sources is wider than any one of them justifies fitting
#: to. It is a wave, it is in the backlog, and the C/N diagnostic should be read as a spread
#: rather than as a bound until it runs.
MELANOIDIN_REPEAT_UNIT_SAME_SYSTEM_AT_COOKING_TEMPERATURE = (
    "Martins & van Boekel 2003 (Food Chem. 83:135), glucose + glycine, measured C/N: 11 at "
    "120 C pH 6.8; 15 falling to 11 at 100 C pH 6.8; 19 falling to 16 at 100 C pH 5.5. Its "
    "Table 2 compiles nine further literature values from 7 to 13 and its authors call them "
    "inconsistent with both pH and temperature. The trunk predicts 8.42 to 9.94."
)


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
