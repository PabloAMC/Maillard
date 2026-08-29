"""
src/kinetic_core/furanic.py

THE FURANIC CHANNEL'S ANALYSIS LAYER (Build Wave B7, 2026-08-29).
=============================================================================

The species are in ``species.py`` and ``species_sulfur.py``, the topology is in
``network.py`` and ``sulfur.py``, and every constant is in
``parameters_furanic.py``. What lives HERE is everything that is a STATEMENT
ABOUT the channel rather than a part of it:

  * the declared-FIT tables, transcribed once, in the units their papers print;
  * the structural constraints, as executable predicates rather than prose;
  * the naming traps, because two of them have already caused corpus-level
    errors and a third is one grep away;
  * the limb-share diagnostic, which exists to make it VISIBLE that the HMF
    branch is a consequence of the pools and not a parameter.

WHY THERE IS A "LIMB SHARE" FUNCTION AND NOT A "BRANCH FRACTION" CONSTANT
-------------------------------------------------------------------------
K5a's governing finding is that the fructose-vs-3-DG question is TWO
QUANTITIES that were being read as one. The fructose limb wins on FLUX --
established by model discrimination in three independent matrices, where
deleting it under-predicts HMF badly. The 3-DG limb wins on PUBLISHED RATE
CONSTANT in four of six comparisons -- because its species are the MEASURED
ones and the fructose limb's intermediate is unmeasured in all five networks,
so a constant on it is identified only up to the pool scale. Those are not in
conflict, and only the flux statement is a chemistry claim.

``hmf_limb_shares`` therefore reports a FLUX share, computed from the state at
a given time, and there is no constant anywhere in this module that could be
edited to change it. The unit tests assert exactly that.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from .network import FURANIC_REACTION_KEYS, REACTIONS as TRUNK_REACTIONS
from .parameters_furanic import (
    BLANK_PH,
    BLANK_STATED_MAX_SD_FRACTION,
    BLANK_SUGAR_LOADING_MMOL_PER_L,
    BLANK_TEMPERATURE_C,
    BLANK_TIME_MIN,
    EDGE_C_ZERO_BY_DECLARATION,
    FURANIC_PARAMETERS,
    FURANONE_EA_ASSUMPTION,
    MU_M_HAZARD,
    PROHIBITED_DERIVATIONS,
)

# ---------------------------------------------------------------------------
# 1. THE NAMING TRAPS -- three of them, two already realised as corpus errors
# ---------------------------------------------------------------------------

NAMING_TRAPS: Mapping[str, str] = {
    "HMF means NORFURANEOL in four papers the repo holds": (
        "whitfield1999.pdf (48 hits) and whitfield2001.pdf (65 hits) use "
        "'HMF' for 4-hydroxy-5-methyl-3(2H)-furanone -- NORFURANEOL, species "
        "``NF`` on the sulfur lane, a C5 compound. Their own titles say so: "
        "'Reaction of 4-hydroxy-5-methyl-3(2H)-furanone (HMF) with cysteine or "
        "hydrogen sulfide at pH 4.5'. blank1996 writes 'HMF (3)' for the same "
        "compound and apriyantono1993 writes 'HMFone', defining it four pages "
        "AFTER the compound first appears in its Table 1. NONE of the four "
        "means 5-hydroxymethylfurfural. A grep-driven ingestion of 'HMF' "
        "across this corpus silently merges the furanone lane into the furanic "
        "lane, and k3 sec. D's 'likely a typo' note about Whitfield 2001 is "
        "not a typo -- it is this collision."
    ),
    "apriyantono1993 names three confusable molecules and measures one": (
        "It writes HMFone (= norfuraneol, which it measures), separately "
        "mentions 5-hydroxymethyl-2-furfural (the real 5-HMF, cited to Tressl "
        "1989, NOT measured) and separately mentions "
        "2,5-dimethyl-4-hydroxy-3(2H)-furanone (furaneol, cited to Shu & Ho "
        "1988, NOT measured). Three confusably named molecules in one paper."
    ),
    "HDMF and DMHF are the same molecule with the letters reordered": (
        "Furaneol is written HDMF by Blank 1996, Blank 1997 and Poisson 2019 "
        "and DMHF by Wang & Ho 2008 and Shu & Ho 1988. The repo uses DMHF. "
        "Homofuraneol is HEMF in Blank 1996 and EHMF in Blank 1997 -- same "
        "lab, same first author, one year apart, same molecule."
    ),
    "the two Kocadagli PDFs are swapped relative to the naive guess": (
        "Kocada2016.pdf (SHORTER stem, 810 kB) is JAFC 10.1021/acs.jafc."
        "6b01862, glucose +/- NaCl, and is the HMF FIT source. "
        "Kocadagli2016.pdf (LONGER stem, 613 kB) is Food Chem 10.1016/"
        "j.foodchem.2016.05.150, glucose/wheat flour -- and its text layer is "
        "CIPHER-GARBLED, so a grep for 'HMF', 'Ea' or any English word returns "
        "zero hits regardless of content. Any negative grep result on that "
        "file is void."
    ),
    "ug per mmol and mg per mol are the same unit": MU_M_HAZARD,
}


# ---------------------------------------------------------------------------
# 2. THE HMF NODE'S TWO PARALLEL SOURCES -- named, so the topology is legible
# ---------------------------------------------------------------------------

#: reaction key -> which limb it belongs to. The ONLY place the network's
#: topology is partitioned into limbs, and it partitions REACTIONS, never
#: assigns a share.
HMF_LIMB_SOURCES: Mapping[str, str] = {
    "r_ddg_hmf": "three_deoxyglucosone_limb",
    "r_int_hmf": "fructose_limb",
}

#: The four independent networks that agree on this topology. Recorded because
#: adopting an architecture on four papers' agreement is a different act from
#: inventing one, and a later wave should be able to tell them apart.
HMF_SOURCE_TOPOLOGY_CORROBORATION: Tuple[str, ...] = (
    "Kocadagli & Gokmen 2016 JAFC (glucose melt): "
    "d[HMF]/dt = k5*[3,4-DG] + k7*[Int] - k18*[HMF]",
    "Kocadagli & Gokmen 2016 Food Chem (wheat flour): same source topology",
    "Goncuoglu Tas & Gokmen 2017 (hazelnut): same source topology",
    "Gursul Aktag & Gokmen 2020 (fruit juice): "
    "d[HMF]/dt = k8*[3-DG] + k11*[FFC], NO sink",
    "Sen & Gokmen 2022 (nuts/seeds): two parallel sources, bimolecular sink",
    "Han 2025: two parallel sources, first-order sink to melanoidin",
)


def hmf_limb_shares(
    state: Mapping[str, float], rate_constants: Mapping[str, float]
) -> Dict[str, Any]:
    """
    The INSTANTANEOUS flux share of each HMF-forming limb, from the state.

    THIS IS A DIAGNOSTIC, NOT A PARAMETER. It reads the two terminal
    concentrations and the two terminal rate constants and divides. Change the
    sugar, the amine, the temperature or the time and it changes, because the
    pools change -- which is precisely K5a sec. 3.1 Rule 1: every paper that
    explains why ITS limb wins appeals to precursor SUPPLY, and in four of six
    published comparisons the terminal rate constant points the other way.

    Returns zero shares and ``defined = False`` before any HMF has formed,
    rather than a NaN that would propagate into a mean.
    """
    dg = float(rate_constants.get("r_ddg_hmf", 0.0)) * float(state.get("DDG", 0.0))
    fru = float(rate_constants.get("r_int_hmf", 0.0)) * float(state.get("INT", 0.0))
    total = dg + fru
    if total <= 0.0:
        return {
            "defined": False,
            "three_deoxyglucosone_limb": 0.0,
            "fructose_limb": 0.0,
            "flux_mmol_per_l_min": 0.0,
            "note": "no HMF flux at this state; the share is undefined, not zero",
        }
    return {
        "defined": True,
        "three_deoxyglucosone_limb": dg / total,
        "fructose_limb": fru / total,
        "flux_mmol_per_l_min": total,
        "note": (
            "an INSTANTANEOUS FLUX share computed from the pools, not a stored "
            "branch fraction. K5a MUST-NOT #1."
        ),
    }


# ---------------------------------------------------------------------------
# 3. THE DECLARED-FIT TABLE -- Blank, Fay, Lakner & Schlosser 1997
# ---------------------------------------------------------------------------
# JAFC 45:2642-2648, article ID JF960997I, CCC line S0021-8561(96)00997-1, all
# verified from the primary document by the K5b wave (which also cleared
# blank1996's 85 %-confidence DOI flag).
#
# ROLE: FIT. Amendment 12, verbatim: "blank1997's 39 SIDA cells = FIT (pentose
# formation edge -- the channel's calibration)."
#
# UNITS: ug per mmol of SUGAR. See MU_M_HAZARD before comparing any of these
# numbers with Wang & Ho's mg/mol of METHYLGLYOXAL.
#
# LEAKAGE WARNING, from K5b sec. 8.4 and honoured by the B7 hold-out scorer:
# the conditions `xylose/glycine, pH 6, phosphate, 1:1` and `xylose/L-alanine,
# pH 6, phosphate, 1:1` EACH APPEAR IN FOUR OF THE FIVE TABLES. The pH-6 column
# of Table 4 IS Table 1. So the pH-ladder hold-out (H2) is scored on the SLOPE
# CONTRAST only -- never on a level -- and the sugar-ordering hold-out is
# scored ordinally.

@dataclass(frozen=True)
class Blank1997Cell:
    """One isotope-dilution cell, in the units the paper prints."""

    table: str
    sugar: str
    amine: str
    compound: str          # "HDMF" (= DMHF) or "EHMF" (= HEMF)
    ug_per_mmol_sugar: float
    role: str              # "fit" | "fit_unrepresentable" | "structural"
    reason: str = ""


#: Table 1 -- six systems x two furanones, plus footnote c's zero-amine
#: control. Table 2 -- the 4-aminobutyric-acid (GABA) fragmentation baseline,
#: which is what makes the 73/27 Strecker split a NON-ISOTOPIC measurement.
_BLANK1997_CELLS: Tuple[Blank1997Cell, ...] = (
    # ---- Table 1, glycine systems: THE FIT ROWS --------------------------
    Blank1997Cell("T1", "arabinose", "glycine", "HDMF", 5.1, "fit"),
    Blank1997Cell("T1", "xylose", "glycine", "HDMF", 2.6, "fit"),
    Blank1997Cell("T1", "ribose", "glycine", "HDMF", 3.6, "fit"),
    # ---- Table 1, glycine systems, the homofuranone ----------------------
    Blank1997Cell(
        "T1", "arabinose", "glycine", "EHMF", 1.3, "fit_unrepresentable",
        "HEMF needs a C2 Strecker donor. The core has no route to it: the "
        "sulfur lane (which owns the pentose) carries no alanine, and the "
        "acrylamide lane (which owns alanine) carries no pentose. Reported as "
        "an unrepresentable FIT row rather than dropped."),
    Blank1997Cell(
        "T1", "xylose", "glycine", "EHMF", 0.3, "fit_unrepresentable",
        "as above. NOTE that it is 0.3 and NOT ZERO -- this single cell is why "
        "Amendment 8's 'HEMF requires alanine' truth table had to be demoted "
        "to a preference by Amendment 12."),
    Blank1997Cell(
        "T1", "ribose", "glycine", "EHMF", 0.7, "fit_unrepresentable", "as above"),
    # ---- Table 1, alanine systems ----------------------------------------
    Blank1997Cell(
        "T1", "arabinose", "alanine", "HDMF", 1.2, "fit_unrepresentable",
        "no lane carries pentose + alanine; see the lane-conflict note"),
    Blank1997Cell(
        "T1", "xylose", "alanine", "HDMF", 0.9, "fit_unrepresentable", "as above"),
    Blank1997Cell(
        "T1", "ribose", "alanine", "HDMF", 1.6, "fit_unrepresentable", "as above"),
    Blank1997Cell(
        "T1", "arabinose", "alanine", "EHMF", 6.8, "fit_unrepresentable", "as above"),
    Blank1997Cell(
        "T1", "xylose", "alanine", "EHMF", 7.5, "fit_unrepresentable", "as above"),
    Blank1997Cell(
        "T1", "ribose", "alanine", "EHMF", 10.0, "fit_unrepresentable", "as above"),
    # ---- Table 1 footnote c: THE CEILING, not a positive control ---------
    Blank1997Cell(
        "T1c", "xylose", "none", "HDMF", 0.01, "structural",
        "UPPER BOUND, reported as '<0.01'. Amendment 12 SIGN-REVERSES "
        "blank1996's proposed hold-out #8: pentose alone is a NEGATIVE control "
        "with a ceiling, >=260x below the xylose/glycine cell, not a positive "
        "one."),
    Blank1997Cell(
        "T1c", "xylose", "none", "EHMF", 0.01, "structural", "as above"),
    # ---- Table 2: the GABA fragmentation baseline ------------------------
    Blank1997Cell(
        "T2", "xylose", "4-aminobutyric acid", "HDMF", 0.4, "structural",
        "GABA is a gamma-amino acid and 'does not decompose by Strecker "
        "deamination', so this is the FRAGMENTATION-ONLY baseline."),
    Blank1997Cell(
        "T2", "xylose", "4-aminobutyric acid", "EHMF", 0.1, "structural", "as above"),
    Blank1997Cell(
        "T2", "xylose", "4-aminobutyric acid + glycine", "HDMF", 1.5, "structural",
        "the increment over the GABA baseline IS the Strecker channel: "
        "(1.5 - 0.4)/1.5 = 73 %."),
    Blank1997Cell(
        "T2", "xylose", "4-aminobutyric acid + alanine", "EHMF", 3.2, "structural",
        "(3.2 - 0.1)/3.2 = 97 %."),
)


def blank1997_fit_cells(role: Optional[str] = None) -> Tuple[Blank1997Cell, ...]:
    """The declared-FIT table, optionally filtered to one role."""
    if role is None:
        return _BLANK1997_CELLS
    return tuple(c for c in _BLANK1997_CELLS if c.role == role)


def blank1997_conditions() -> Dict[str, Any]:
    """Blank's single condition set, in the core's own units."""
    return {
        "sugar_mmol_per_l": BLANK_SUGAR_LOADING_MMOL_PER_L,
        "amine_mmol_per_l": BLANK_SUGAR_LOADING_MMOL_PER_L,
        "temperature_c": BLANK_TEMPERATURE_C,
        "time_min": BLANK_TIME_MIN,
        "ph": BLANK_PH,
        "buffer": "0.2 M phosphate, set-point only; hold not stated and no "
                  "final pH reported",
        "stated_max_sd_fraction": BLANK_STATED_MAX_SD_FRACTION,
        "replication": "n >= 2 assays x 2 injections, maximum SD <= 10 %",
        "unit": "ug per mmol of SUGAR (= mg per mol of sugar, exactly)",
    }


def blank1997_structural_summary() -> Dict[str, Any]:
    """
    The DERIVED structural quantities, computed from the table above.

    Computed rather than transcribed so that the numbers in the B7 reports
    cannot drift from the cells they come from.
    """
    def cell(table: str, sugar: str, amine: str, compound: str) -> float:
        for c in _BLANK1997_CELLS:
            if (c.table, c.sugar, c.amine, c.compound) == (table, sugar, amine, compound):
                return c.ug_per_mmol_sugar
        raise KeyError((table, sugar, amine, compound))

    gaba_h = cell("T2", "xylose", "4-aminobutyric acid", "HDMF")
    gaba_gly_h = cell("T2", "xylose", "4-aminobutyric acid + glycine", "HDMF")
    gaba_e = cell("T2", "xylose", "4-aminobutyric acid", "EHMF")
    gaba_ala_e = cell("T2", "xylose", "4-aminobutyric acid + alanine", "EHMF")

    ehmf_pref = {
        sugar: cell("T1", sugar, "alanine", "EHMF") / cell("T1", sugar, "glycine", "EHMF")
        for sugar in ("arabinose", "xylose", "ribose")
    }
    return {
        "strecker_share_HDMF": (gaba_gly_h - gaba_h) / gaba_gly_h,
        "fragmentation_share_HDMF": gaba_h / gaba_gly_h,
        "strecker_share_EHMF": (gaba_ala_e - gaba_e) / gaba_ala_e,
        "corroborating_isotopomer_split": {
            "source": "Blank & Fay 1996, 13C labelling, xylose/glycine, 90 C",
            "strecker": 0.70,
            "fragmentation": 0.30,
            "note": (
                "73/27 by amino-acid substitution against 70/30 by isotopomer "
                "distribution: two completely different methods, three "
                "percentage points apart, on the same system. The strongest "
                "corroboration in the K5b cluster."
            ),
        },
        "hemf_alanine_preference_by_sugar": ehmf_pref,
        "hemf_alanine_preference_range": (
            min(ehmf_pref.values()), max(ehmf_pref.values())
        ),
        "hdmf_sugar_ordering": ("arabinose", "ribose", "xylose"),
        "ehmf_sugar_ordering": ("ribose", "xylose", "arabinose"),
        "hdmf_ehmf_selectivity_swing_xylose": (
            cell("T1", "xylose", "glycine", "HDMF")
            / cell("T1", "xylose", "glycine", "EHMF")
        ) / (
            cell("T1", "xylose", "alanine", "HDMF")
            / cell("T1", "xylose", "alanine", "EHMF")
        ),
        "pentose_alone_ceiling_ug_per_mmol": 0.01,
    }


# ---------------------------------------------------------------------------
# 4. THE STRUCTURAL CONSTRAINTS, AS PREDICATES
# ---------------------------------------------------------------------------

FURANONE_STRUCTURAL_CONSTRAINTS: Mapping[str, str] = {
    "edge_B_is_acetylformoin_free": (
        "No reaction that has methylglyoxal among its REACTANTS may have "
        "acetylformoin among its PRODUCTS. Wang & Ho 2008 sec. 4.2, verbatim: "
        "'no [12C6]acetylformoin was observed suggesting that acetylformoin "
        "was not a precursor during the DMHF formation from MG'. This is a "
        "measured NEGATIVE, it resolves an ambiguity Poisson 2019 explicitly "
        "leaves open (its pathway b vs c), and it is enforced by "
        "``validate_furanic_structure`` at import."
    ),
    "the_edges_do_not_mix": (
        "In both of Wang & Ho's CAMOLA experiments only [13C6] and [12C6] "
        "isotopomers were ever observed -- never [13C3]. No DMHF is assembled "
        "from one hexose-derived C3 plus one MG-derived C3. The network "
        "honours it: Edge B consumes TWO methylglyoxals and nothing else."
    ),
    "hexose_DMHF_keeps_the_intact_C6_skeleton": (
        "Measured twice: Wang & Ho ([13C6]:[12C6] = 1:1 with [13C1]-[13C5] "
        "absent, and acetylformoin itself at m/z 144:150 = 1:1) and Poisson "
        "2019 in a real bean (intact share 87.3-100 % at all nine roast times, "
        "with 2,3-butanedione's intact share collapsing 25.4 -> 0.4 % in the "
        "SAME runs as an internal positive control). The hexose edge is "
        "therefore FIRST ORDER IN THE DEOXYOSONE ALONE and takes no carbon "
        "from any amine."
    ),
    "pentose_DMHF_needs_a_Strecker_carbon": (
        "C5 + C1 -> C6. The pentose edge is bimolecular in the amine because "
        "the sixth carbon comes from the amino acid's Strecker aldehyde, and "
        "73 % of it does (Blank 1997's GABA control), corroborated at 70 % by "
        "13C labelling."
    ),
    "norfuraneol_dominates_DMHF_at_the_deoxypentosone_fork": (
        "Blank 1997 p. 2646: 'in all samples analyzed, 4-hydroxy-5-methyl-"
        "3(2H)-furanone was the main reaction product (data not shown)'. "
        "Corroborated across two papers and six systems and QUANTIFIED IN "
        "NEITHER, so this is an ORDINAL ceiling and can never be a ratio. "
        "Apriyantono's near-null norfuraneol cell must NOT be scored against "
        "it (K5b sec. 8.5): both terms are at the detection floor, the amine "
        "is lysine, and the isolation is SDE."
    ),
    "DMHF_from_pentose_alone_is_a_CEILING": (
        "< 0.01 ug/mmol for BOTH furanones, in phosphate AND in water -- "
        ">=260x below the xylose/glycine cell and reported only as an upper "
        "bound. Amendment 12 SIGN-REVERSES blank1996's proposal, which had it "
        "as a positive control."
    ),
    "HEMF_alanine_preference_is_a_RATIO_not_a_SWITCH": (
        "Amendment 12 demotes Amendment 8's truth table. Blank 1997 measures "
        "EHMF at 0.3-1.3 ug/mmol in ALL THREE glycine systems, so a model that "
        "emits HEMF from pentose + glycine is RIGHT."
    ),
    "a_3DG_only_HMF_node_is_falsified": (
        "Three independent matrices, one lab, one method, three "
        "model-discrimination tests, one answer: deleting the fructose/cation "
        "edge under-predicts HMF badly ('by no means', 'remarkably "
        "underestimated', 'far below the experimental values'). K5a C1."
    ),
    "the_3DG_limb_and_the_fructose_limb_have_MIRROR_IMAGE_internal_structure": (
        "3-DG -> 3,4-DG is the RDS of the 3-DG limb and 3,4-DG -> HMF is fast "
        "(2-5x, two matrices). Fructose -> cation is FAST and cation -> HMF is "
        "the RDS. K5a C3/C4. A model that lumps either limb into one edge "
        "loses this, so neither is lumped."
    ),
}


def validate_furanic_structure(
    reactions: Sequence[Any] = TRUNK_REACTIONS,
    sulfur_reactions: Optional[Sequence[Any]] = None,
) -> None:
    """
    Enforce the structural constraints that are properties of the TOPOLOGY.

    Called at import so that a later edit which, say, routes methylglyoxal
    through acetylformoin fails loudly instead of quietly contradicting a
    measured null.
    """
    from .parameters_furanic import MEASURED_FURANIC

    all_reactions = list(reactions) + list(sulfur_reactions or ())

    # 1. Edge B is acetylformoin-free -- Wang & Ho's measured null.
    for reaction in all_reactions:
        if "MGO" in reaction.reactants and "AF" in reaction.products:
            raise ValueError(
                f"{reaction.key}: routes methylglyoxal through acetylformoin. "
                f"Wang & Ho 2008 sec. 4.2 measured that this does NOT happen "
                f"(no [12C6]acetylformoin in the MG-spiked run). "
                f"{FURANONE_STRUCTURAL_CONSTRAINTS['edge_B_is_acetylformoin_free']}"
            )

    # 2. The two HMF sources exist, are PARALLEL, and are exactly two.
    sources = {
        r.key for r in all_reactions
        if "HMF" in r.products and r.key in HMF_LIMB_SOURCES
    }
    if sources != set(HMF_LIMB_SOURCES):
        raise ValueError(
            f"the HMF node must carry EXACTLY the two parallel sources four "
            f"independent networks agree on (K5a sec. 8.1); found {sorted(sources)}"
        )

    # 3. The hexose furanone edge takes no amine carbon (intact C6, measured).
    hexose = next(r for r in all_reactions if r.key == "r_odg_af")
    if set(hexose.reactants) != {"ODG"}:
        raise ValueError(
            "r_odg_af must be first order in the 1-deoxyosone ALONE: the "
            "hexose DMHF route keeps the intact C6 skeleton, measured by two "
            "independent CAMOLA experiments with an internal positive control."
        )

    # 4. Edge C is present and its rate is exactly zero.
    edge_c = MEASURED_FURANIC["k_dmhf_h2s"]
    if edge_c.k_ref != 0.0:
        raise ValueError(
            "Edge C must ship at EXACTLY zero. Shu & Ho 1988 reports a GC area "
            "percent and no consumption of any kind; "
            f"{EDGE_C_ZERO_BY_DECLARATION['why']}"
        )

    # 5. No constant anywhere in the module encodes an HMF branch fraction.
    #    A dimensionless number in (0,1) attached to an HMF-forming edge would
    #    be exactly the thing K5a MUST-NOT #1 forbids.
    for key in ("k_ddg_hmf", "k_int_hmf"):
        parameter = MEASURED_FURANIC[key]
        if parameter.unit != "1/min":
            raise ValueError(
                f"{key} must be a first-order RATE, not a dimensionless share. "
                "The HMF branch is set by the pools; there is no fraction."
            )


# ---------------------------------------------------------------------------
# 5. The declared-assumption band, and its bookkeeping
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FuranoneYieldModel:
    """
    The frozen furanone calibration: one number, and what it was fitted to.

    Kept as an object rather than a bare float so that the B7 fit report, the
    engine and the tests all read the SAME provenance, and so that a future
    wave that adds a second fitted furanone constant has somewhere to put it.
    """

    k_dpo_af_l_per_mmol_min: float
    fitted_on: Tuple[str, ...]
    residual_log10_fold: Mapping[str, float]
    sigma_log10: float

    def as_dict(self) -> Dict[str, Any]:
        return {
            "k_dpo_af_l_per_mmol_min": self.k_dpo_af_l_per_mmol_min,
            "fitted_on": list(self.fitted_on),
            "residual_log10_fold": dict(self.residual_log10_fold),
            "sigma_log10": self.sigma_log10,
        }


def furanone_yield_model_from_dict(payload: Mapping[str, Any]) -> FuranoneYieldModel:
    return FuranoneYieldModel(
        k_dpo_af_l_per_mmol_min=float(payload["k_dpo_af_l_per_mmol_min"]),
        fitted_on=tuple(payload.get("fitted_on", ())),
        residual_log10_fold=dict(payload.get("residual_log10_fold", {})),
        sigma_log10=float(payload.get("sigma_log10", 0.0)),
    )


def furanone_declared_assumption_corners() -> Tuple[float, float, float]:
    """The three partition-barrier values the interval is priced at."""
    ea = float(FURANONE_EA_ASSUMPTION["value_kj_mol"])
    band = float(FURANONE_EA_ASSUMPTION["band_kj_mol"])
    return (ea - band, ea, ea + band)


def ug_per_mmol_from_state(
    concentration_mmol_per_l: float,
    sugar_charge_mmol_per_l: float,
    molecular_weight_g_per_mol: float,
) -> float:
    """
    mmol/L of a furanone -> ug per mmol of SUGAR, Blank 1997's unit.

    Written out because it is the conversion the whole DMHF calibration turns
    on and because the target unit is the one MU_M_HAZARD is about.
    """
    if sugar_charge_mmol_per_l <= 0.0:
        raise ValueError("a per-mmol-of-sugar yield needs a sugar charge")
    ug_per_l = (
        float(concentration_mmol_per_l) * float(molecular_weight_g_per_mol) * 1.0e3
    )
    return ug_per_l / float(sugar_charge_mmol_per_l)


def describe_furanic() -> Dict[str, Any]:
    """A machine-readable description of the channel, for the reports."""
    return {
        "module": "kinetic_core.furanic",
        "wave": "B7",
        "reactions_on_the_trunk": list(FURANIC_REACTION_KEYS),
        "reactions_on_the_sulfur_lane": ["r_dpo_af", "r_hmf_cys", "r_dmhf_h2s"],
        "hmf_limb_sources": dict(HMF_LIMB_SOURCES),
        "hmf_source_topology_corroboration": list(
            HMF_SOURCE_TOPOLOGY_CORROBORATION
        ),
        "structural_constraints": dict(FURANONE_STRUCTURAL_CONSTRAINTS),
        "naming_traps": dict(NAMING_TRAPS),
        "prohibited_derivations": dict(PROHIBITED_DERIVATIONS),
        "edge_c": dict(EDGE_C_ZERO_BY_DECLARATION),
        "furanone_ea_assumption": dict(FURANONE_EA_ASSUMPTION),
        "blank1997_conditions": blank1997_conditions(),
        "blank1997_structural_summary": blank1997_structural_summary(),
        "n_parameters": len(FURANIC_PARAMETERS),
        "declared_gaps": list(DECLARED_GAPS),
    }


DECLARED_GAPS: Tuple[str, ...] = (
    "NO usable activation energy exists for either HMF-forming edge in any "
    "REAL FOOD MATRIX. The only reproducible one (Int -> HMF, 145-152 kJ/mol) "
    "is from an amine-free freeze-dried glucose melt over an UNMEASURED "
    "intermediate pool, and all three real-matrix systems in the corpus "
    "destroy it (wheat flour Ea = -97.6 with R^2 = 0.322; hazelnut "
    "non-monotonic; five nuts/seeds with k = 0 in six of ten cells). K5a G1.",
    "NO HMF-sink constant exists above 50 C that survives audit. Hamzalioglu "
    "covers 5-50 C; the hazelnut first-order sink covers 150-160 C but is a "
    "lumped decay with no amine dependence. THE 50-150 C WINDOW IS EMPTY, so "
    "this module has no effective HMF sink at cooking temperature and must be "
    "expected to OVER-PREDICT HMF. K5a G2.",
    "The fructose-limb intermediate has NEVER been measured in any of the five "
    "published networks. Until [Int] or [FFC] is quantified, no absolute rate "
    "constant on that limb is recoverable from any of them. K5a G3.",
    "NO rate constant, NO activation energy, NO time course and NO temperature "
    "series exists anywhere in the DMHF cluster. All five papers are "
    "single-temperature. K5b sec. 10.",
    "There is NO absolute HEXOSE DMHF yield in any of the five K5b papers. The "
    "intact-C6 STRUCTURE is settled twice over and the MAGNITUDE is measured "
    "nowhere; this module transfers the pentose level by a declared rule and "
    "flags every hexose answer.",
    "There is NO measurement of DMHF CONSUMPTION. Edge C ships at exactly "
    "zero. Haleva-Toledo et al. 1999 (JAFC 47:4140-4145) is the identified "
    "source that would close it, and it would serve the HMF node at the same "
    "time.",
    "NO pH ladder exists anywhere in the HMF cluster. Six distinct pH values "
    "appear across the seven papers and NO SINGLE PAPER VARIES pH, so any pH "
    "term on the HMF node would be fitted across labs and matrices at once -- "
    "which k3 sec. B.2 already forbids at family level. This module therefore "
    "carries NO pH term on any HMF edge. K5a G8.",
    "3-deoxyglucosone and 3-deoxyGALACTOSONE were not chromatographically "
    "resolved in ANY of the five Gokmen multiresponse papers (author-declared "
    "once, same method throughout). Every 3-DG number ingested here carries "
    "that ambiguity. K5a G6.",
    "3,4-dideoxyglucosone is SEMI-QUANTITATED against the 3-DG response "
    "factor, so both edges that touch it inherit an unknown multiplicative "
    "scale. K5a C22.",
)


validate_furanic_structure()


__all__ = [
    "Blank1997Cell",
    "DECLARED_GAPS",
    "FURANONE_STRUCTURAL_CONSTRAINTS",
    "FuranoneYieldModel",
    "HMF_LIMB_SOURCES",
    "HMF_SOURCE_TOPOLOGY_CORROBORATION",
    "NAMING_TRAPS",
    "blank1997_conditions",
    "blank1997_fit_cells",
    "blank1997_structural_summary",
    "describe_furanic",
    "furanone_declared_assumption_corners",
    "furanone_yield_model_from_dict",
    "hmf_limb_shares",
    "ug_per_mmol_from_state",
    "validate_furanic_structure",
]
