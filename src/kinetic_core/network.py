"""
src/kinetic_core/network.py

THE MASS-ACTION REACTION NETWORK OF THE TRUNK, WITH THE MELANOIDIN MASS SINK
(Build Wave B1, 2026-08-28; extended by B7).
=============================================================================

WAVE: B1 for the fifteen trunk steps; B7 hung the FURANIC CHANNEL here (the HMF
node and DMHF Edges A and B, including ``r_mgo_dmhf``), which is why the sulfur
lane's residuals moved in B7 without any sulfur constant being touched -- the
sulfur lane runs this trunk.
EXAM: the cutover final exam (results/validation/cutover_final_exam.md) scores
this network's predictions; B1's own fit and hold-out records are
results/validation/kinetic_core_b1_{fit,holdout}_report.md.
DECLARED GAPS: the melanoidin sink is an ELEMENTAL pool, not a molecule, so it
has no molecular weight and is reported in its own unit; the furanic edges
carry NO activation energy from any source (all five papers of the cluster are
single-temperature), so their partition barrier is a declared assumption priced
by re-integration rather than a measurement; and ``k_mgo_dmhf``'s level is a
digitised bar-chart value at one temperature. The balance invariant, by
contrast, is not a gap: this module REFUSES TO IMPORT if any step fails to
balance carbon or nitrogen.

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
from .parameters_furanic import FURANIC_PARAMETERS
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
    # =======================================================================
    # BUILD WAVE B7 -- THE FURANIC CHANNEL, eleven steps.
    # =======================================================================
    # It hangs HERE, on the trunk, rather than in a lane of its own, because
    # all four of its parents are trunk species: fructose and 3-deoxyglucosone
    # for HMF, 1-deoxyglucosone and methylglyoxal for DMHF. Living on the trunk
    # means the sulfur and acrylamide lanes inherit it without any lane
    # composing with any other -- there is no new lane conflict to resolve.
    #
    # THE HMF NODE'S ARCHITECTURE IS NOT THIS MODULE'S INVENTION. Four
    # independent groups fitted four multiresponse networks in four matrices
    # (Kocadagli x2, Goncuoglu Tas, Gursul Aktag, Sen, Han) and ALL FOUR write
    # the same source topology: EXACTLY TWO PARALLEL FIRST-ORDER INPUTS, one
    # from the 3-DG/3,4-DG chain and one from the fructose/cation chain. K5a
    # sec. 8.1 calls that the strongest architectural agreement in the cluster
    # and says it should be adopted without modification. It is.
    #
    # THERE IS NO BRANCH FRACTION HERE AND THERE CANNOT BE ONE. The share each
    # limb takes is whatever the dynamic Fru and TDG pools make it, which is
    # what every paper that explains its own verdict actually appeals to:
    # pool size (Gursul Aktag's fructose-rich juices), a starved 3-DG source
    # (Kocadagli's k3, the smallest constant in his table), or a drained
    # cation pool (Sen's k3/k20 at 10-300x the dehydration step). K5a sec. 3.1
    # Rule 1, and MUST-NOT #1.
    Reaction(
        "r_glc_tdg", {"Glc": 1}, {"TDG": 1}, "k_glc_tdg",
        "B7. Kocadagli JAFC step 3, the AMINE-FREE entry to the 3-DG limb. "
        "B1's trunk reaches 3-deoxyglucosone only through the Amadori "
        "compound, so before B7 a sugar-only pot had no 3-DG limb at all and "
        "a glucose/alanine pot had no deoxyosone of any kind -- the trunk's "
        "only amine is glycine.",
    ),
    Reaction(
        "r_tdg_ddg", {"TDG": 1}, {"DDG": 1}, "k_tdg_ddg",
        "B7. Kocadagli JAFC step 4. THE RATE-DETERMINING STEP OF THE 3-DG "
        "LIMB (K5a C3), corroborated in two independent matrices. Its product "
        "is semi-quantitated against the 3-DG response factor, so both edges "
        "that touch DDG inherit an unknown multiplicative scale (C22).",
    ),
    Reaction(
        "r_ddg_hmf", {"DDG": 1}, {"HMF": 1}, "k_ddg_hmf",
        "B7. Kocadagli JAFC step 5. FAST relative to its parent (2-5x, two "
        "matrices) and carrying the authors' OWN Ea = 0: their k runs "
        "160 -> 110 -> 137 across 160 -> 180 -> 200 C and they fixed the "
        "barrier to zero with the footnote 'does not follow Arrhenius "
        "equation'. No usable Ea exists for this edge in any paper of the "
        "cluster.",
    ),
    Reaction(
        "r_fru_int", {"Fru": 1}, {"INT": 1}, "k_fru_int",
        "B7. Kocadagli JAFC step 6. The fructose limb's FAST entry step -- the "
        "mirror image of the 3-DG limb, whose fast step is the second one "
        "(K5a C3/C4). Deleting this edge is what THREE independent model-"
        "discrimination tests in THREE matrices reject: 'did not fit ... by no "
        "means', 'remarkably underestimated', 'far below the experimental "
        "values'. A 3-DG-only HMF node is falsified three times over (C1).",
    ),
    Reaction(
        "r_int_hmf", {"INT": 1}, {"HMF": 1}, "k_int_hmf",
        "B7. Kocadagli JAFC step 7. The fructose limb's RATE-DETERMINING step. "
        "[Int] is UNMEASURED, so this constant and k_fru_int are identified "
        "only up to a common pool scale and neither may be compared in "
        "magnitude with a constant on the measured 3-DG limb (C2).",
    ),
    Reaction(
        "r_fru_odg", {"Fru": 1}, {"ODG": 1}, "k_fru_odg",
        "B7. Kocadagli JAFC step 8, the AMINE-FREE entry to the furanone limb "
        "(2,3-enolisation). K5a sec. 5 row 4: the parent of 1-DG has four "
        "different answers in four matrices from one lab. B1 carries the "
        "Amadori parent; this adds the melt's, so the split between them is "
        "set by the pools rather than by a constant.",
    ),
    Reaction(
        "r_tdg_mgo", {"TDG": 1}, {"MGO": 1, "FRAG_C": 3}, "k_tdg_mgo",
        "B7. Kocadagli JAFC step 11, the AMINE-FREE methylglyoxal source and "
        "therefore DMHF Edge B's feed in a sugar-only pot. The C3 residue is "
        "unmeasured and goes to FRAG_C, B1's discipline unchanged. The parent "
        "of MGO switches between matrices in the same lab (K5a sec. 5 row 3); "
        "both parents now exist and neither is hard-coded as dominant.",
    ),
    Reaction(
        "r_hmf_self", {"HMF": 1}, {"FRAG_C": 6}, "k_hmf_self",
        "B7. HMF self-degradation, from Hamzalioglu's model-free control "
        "(0.9 % in 7 days at 5 C, pH 3.5). ONE TEMPERATURE, so Ea = 0 by "
        "declaration -- which makes this sink negligible at cooking "
        "temperature and leaves the model with NO validated HMF sink there. "
        "K5a declared gap G2: the 50-150 C window is empty. Pre-registered "
        "consequence: HMF is expected to be OVER-predicted.",
    ),
    Reaction(
        "r_odg_af", {"ODG": 1}, {"AF": 1}, "k_odg_af",
        "B7. DMHF Edge A, HEXOSE arm: FIRST ORDER IN THE DEOXYOSONE ALONE, "
        "because the C6 skeleton stays intact and needs no Strecker carbon. "
        "That is measured twice over -- Wang & Ho's CAMOLA ([13C1]-[13C5] "
        "ABSENT) and Poisson's in-bean CAMOLA (intact share 87-100 % at all "
        "nine roast times, against 2,3-butanedione's collapsing 25.4 -> 0.4 % "
        "in the same runs). The STRUCTURE is measured; the LEVEL is a declared "
        "transfer from the pentose fit, because no hexose DMHF magnitude "
        "exists anywhere in the cluster.",
    ),
    Reaction(
        "r_af_dmhf", {"AF": 1}, {"DMHF": 1}, "k_af_dmhf",
        "B7. Acetylformoin reduction to DMHF. Its RATE has no source; the "
        "constant is a declared, unconstrained 'not rate-limiting' assumption "
        "in the register of thiol_addition_pentodiulose, swept over three "
        "decades in the fit report.",
    ),
    Reaction(
        "r_mgo_dmhf", {"MGO": 2}, {"DMHF": 1}, "k_mgo_dmhf",
        "B7. DMHF Edge B, the C3 + C3 recombination. IT DOES NOT PASS THROUGH "
        "ACETYLFORMOIN and the network makes that structurally impossible: no "
        "reaction with MGO as a reactant has AF among its products. That is "
        "Wang & Ho's measured null (no [12C6]acetylformoin in the MG-spiked "
        "run) and it resolves an ambiguity Poisson 2019 leaves open. The LEVEL "
        "is a digitised bar-chart prior and is flagged as such everywhere.",
    ),
)

REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in REACTIONS)

#: Build Wave B13 (2026-09-07): the dicarbonyl trio, TRUNK-ONLY. These five steps run
#: when the trunk lane integrates on its own (`TRUNK_REACTIONS`); the sulfur network
#: imports `REACTIONS` and keeps exactly the topology wave B9 was fitted on, and the
#: acrylamide network builds its own. Constants: `parameters_dicarbonyl.py`
#: (Kocadagli & Gokmen 2016 JAFC, glucose glass, T_b 180 C, re-referenced to 100 C).
DICARBONYL_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_glc_g", {"Glc": 1}, {"G": 1}, "k_glc_g",
        "Kocadagli step 9, Glc -> glucosone (oxidative; the source's amine-free glass "
        "holds air). Ea 125.9 +/- 4.9 kJ/mol.",
    ),
    Reaction(
        "r_g_go", {"G": 1}, {"GO": 1, "FRAG_C": 4}, "k_g_go",
        "Kocadagli step 10, glucosone -> glyoxal + C4 residue (unmeasured -> FRAG_C). "
        "Ea 93.8 +/- 6.4.",
    ),
    Reaction(
        "r_odg_da", {"ODG": 1}, {"DA": 1, "FRAG_C": 2}, "k_odg_da",
        "Kocadagli step 12, 1-deoxyglucosone -> diacetyl + C2 residue. Ea 150.8 +/- 8.8.",
    ),
    Reaction(
        "r_go_sink", {"GO": 1}, {"FRAG_C": 2}, "k_go_sink",
        "Kocadagli step 15, GO -> P3. Rate measured at 180 C; barrier FIXED TO ZERO by the "
        "authors, carried as such and flagged.",
    ),
    Reaction(
        "r_da_sink", {"DA": 1}, {"FRAG_C": 4}, "k_da_sink",
        "Kocadagli step 17, DA -> P5: rate 0 +/- 0 in the source (diacetyl accumulates). "
        "Carried at ZERO as a prediction, not left undefined.",
    ),
)
#: Build Wave B18 (2026-09-08): the pyrazine step, TRUNK-ONLY, five steps. Two Strecker
#: deaminations (dicarbonyl + glycine -> aminoketone + CO2 + formaldehyde; hypothesis-layer rule
#: R07), second order and RATE-DETERMINING, then three aminoketone condensations (rule R28) on one
#: shared constant DECLARED FAST (Jousse 2002's "I + I -> pyrazines: fast", jousse2002_extraction.md
#: Table R10), so the measured pyrazine rate is the Strecker rate over two and the mixed pyrazine
#: follows the two aminoketone pools statistically. Each glycine leaves its two carbons as carbon
#: dioxide and formaldehyde, booked to the unassigned fragment pool. Constants:
#: `parameters_pyrazine.py`; pre-registration `results/validation/kinetic_core_b18_prereg.md`.
PYRAZINE_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_go_ak", {"GO": 1, "Gly": 1}, {"AKG": 1, "FRAG_C": 2}, "k_go_ak",
        "B18. glyoxal + glycine -> aminoacetaldehyde + CO2 + HCHO (Strecker, net). FITTED to Zhou "
        "2024's three-temperature pyrazine formation rates on fed glyoxal + alanine (alanine -> "
        "glycine declared): the rate-determining step of the pyrazine route.",
    ),
    Reaction(
        "r_mgo_ak", {"MGO": 1, "Gly": 1}, {"AKM": 1, "FRAG_C": 2}, "k_mgo_ak",
        "B18. methylglyoxal + glycine -> aminoacetone + CO2 + HCHO (Strecker, net). FITTED to Zhou "
        "2024's three-temperature 2,5-dimethylpyrazine rates on fed methylglyoxal + alanine.",
    ),
    Reaction(
        "r_akg_pz", {"AKG": 2}, {"PZ": 1}, "k_cond",
        "B18. 2 aminoacetaldehyde -> pyrazine (condensation, dehydration, oxidation; net). DECLARED "
        "FAST (shared k_cond): not rate-determining, sensitivity reported in the ship rule.",
    ),
    Reaction(
        "r_akm_dmp", {"AKM": 2}, {"DMP": 1}, "k_cond",
        "B18. 2 aminoacetone -> 2,5-dimethylpyrazine (net). Shared declared k_cond.",
    ),
    Reaction(
        "r_ak_mpz", {"AKG": 1, "AKM": 1}, {"MPZ": 1}, "k_cond",
        "B18. aminoacetaldehyde + aminoacetone -> 2-methylpyrazine (net). Shared declared k_cond: the "
        "mixed pyrazine follows the two pools statistically (2 sqrt of the two homo rates); no "
        "source measures the mixed condensation.",
    ),
)
#: Build Wave B20 (2026-09-09): THE GLYCATION ARM, TRUNK-ONLY, five steps on protein-bound lysine.
#: The sugar glycates the bound lysine to the bound Amadori compound (second order, Nguyen 2016 k3);
#: the Amadori compound oxidises to CML (k7; Berk 2021's barrier), goes to CEL via methylglyoxal
#: lumped (k9), or decays to 3-deoxyglucosone and gives the lysine back (k8, the dominant loss, the
#: trunk's own Amadori-decay form); CML is lost into the melanoidin pools (k11). Every constant is
#: fitted on Nguyen 2016's printed rates; the barriers are declared from measured ones. With no
#: protein loading LYSP is zero and the five steps carry no flux, so every earlier pot reproduces.
GLYCATION_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_glc_lysp", {"Glc": 1, "LYSP": 1}, {"FLP": 1}, "k_glyc",
        "B20. glucose + bound lysine -> bound fructosyl-lysine (Schiff base and Amadori rearrangement "
        "lumped, as Nguyen 2016 fitted k3; second order, L/(mmol*min)).",
    ),
    Reaction(
        "r_flp_cml", {"FLP": 1}, {"CML": 1, "FRAG_C": 4}, "k_flp_cml",
        "B20. bound fructosyl-lysine -> CML + C4 fragments (oxidative cleavage; Nguyen 2016 k7, Berk 2021 k8).",
    ),
    Reaction(
        "r_flp_cel", {"FLP": 1}, {"CEL": 1, "FRAG_C": 3}, "k_flp_cel",
        "B20. bound fructosyl-lysine -> CEL + C3 fragments (via methylglyoxal, lumped; Nguyen 2016 k9).",
    ),
    Reaction(
        "r_flp_decay", {"FLP": 1}, {"TDG": 1, "LYSP": 1}, "k_flp_decay",
        "B20. bound fructosyl-lysine -> 3-deoxyglucosone + bound lysine (the Amadori decay that returns "
        "the amine, the trunk's r_ama_tdg form; Nguyen 2016 k8, 'AP -> MRPs', the dominant loss).",
    ),
    Reaction(
        "r_cml_loss", {"CML": 1}, {"MEL_C": 8, "MEL_N": 2}, "k_cml_loss",
        "B20. CML -> melanoidin pools (Nguyen 2016 k11; sets the CML plateau with k7).",
    ),
)
#: Build Wave B21 (2026-09-09): the aqueous glucosone route. In water the glucosone comes from the
#: Amadori compound (Hamzalioglu 2026), not from the sugar as in the B13 glass; the step returns the
#: amine like the trunk's other Amadori decays. Pre-registered in kinetic_core_b21_prereg.md.
AQUEOUS_GLYOXAL_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_ama_g", {"AMA": 1}, {"G": 1, "Gly": 1}, "k_ama_g",
        "B21. Amadori (DFG) -> glucosone + glycine (oxidative cleavage; the aqueous glucosone entry, "
        "Hamzalioglu 2026 step 4). FITTED on four first-order constants at 110-140 C.",
    ),
)
#: Build Wave B22 (2026-09-09): THE METHIONINE CHAIN, TRUNK-ONLY. The Strecker step of B18 with
#: methionine as the amino acid: the dicarbonyl keeps its carbons in the aminoketone (glycine's
#: AKG / AKM by construction), methionine leaves as methional and CO2 (to FRAG_C). Then the
#: retro-Michael release of methanethiol (acrolein to FRAG_C) and the disulfide on an apparent
#: constant. Pre-registered in kinetic_core_b22_prereg.md; constants in parameters_methionine.py.
METHIONINE_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_go_met", {"GO": 1, "MET": 1}, {"AKG": 1, "MTAL": 1, "FRAG_C": 1}, "k_go_met",
        "B22. glyoxal + methionine -> aminoacetaldehyde + methional + CO2 (Strecker, net). The identity "
        "ratio to glycine's k_go_ak is FITTED on Pan 2025's methional rates; barrier and pH term are B18's.",
    ),
    Reaction(
        "r_mgo_met", {"MGO": 1, "MET": 1}, {"AKM": 1, "MTAL": 1, "FRAG_C": 1}, "k_mgo_met",
        "B22. methylglyoxal + methionine -> aminoacetone + methional + CO2 (Strecker, net). Same ratio.",
    ),
    Reaction(
        "r_mtal_msh", {"MTAL": 1}, {"MSH": 1, "FRAG_C": 3}, "k_mtal_msh",
        "B22. methional -> methanethiol + acrolein (retro-Michael; acrolein to the fragment pool). FITTED "
        "on Pan 2025's methanethiol rates, barrier free within 20-150 kJ/mol.",
    ),
    Reaction(
        "r_msh_dmds", {"MSH": 2}, {"DMDS": 1}, "k_msh_dmds",
        "B22. 2 methanethiol -> dimethyl disulfide, an APPARENT second-order constant (the pot's internal "
        "oxidant is not tracked on this lane; Xu 2010). FITTED on Pan 2025's disulfide rates.",
    ),
    Reaction(
        "r_marp_mtal", {"MARP": 1}, {"MTAL": 1, "MEL_C": 7, "MEL_N": 1}, "k_marp_mtal",
        "B22b. The methionine Amadori compound decomposes to methional, first order. THE ROUTE DENG'S OWN "
        "EXPERIMENT NAMES: the fed Amadori compound gives 1.4 to 2.6 times more methional than methionine "
        "plus glucose, so the dicarbonyl arrives inside the molecule rather than as a free pool. Carbon "
        "closes as 11 = 4 + 7 and the NITROGEN goes to the melanoidin pool, not to the fragment pool, "
        "which holds carbon only: the Strecker aldehyde takes no nitrogen with it and the residue is a "
        "nitrogen-bearing sugar fragment the source does not measure. FITTED on Deng 2022 Table 1's "
        "five-point time course at 120 C.",
    ),
    Reaction(
        "r_marp_loss", {"MARP": 1}, {"MEL_C": 11, "MEL_N": 1, "MEL_S": 1}, "k_marp_loss",
        "B22b. The Amadori compound's own competing loss, first order -- what every other Amadori compound "
        "in this model has. It is here because Deng's series RISES to 120 minutes and then FALLS, and a "
        "single first-order decomposition of a fed pool saturates rather than falling. The sink keeps both "
        "the carbon and the nitrogen; it is accounting, and the source names no product.",
    ),
)
#: Build Wave B24 (2026-09-09): 2-ACETYL-1-PYRROLINE FROM PROLINE, TRUNK-ONLY. Pre-registered in
#: kinetic_core_b24_prereg.md; constants in parameters_proline.py.
PROLINE_REACTIONS: Tuple[Reaction, ...] = (
    Reaction(
        "r_mgo_pro", {"MGO": 1, "PRO": 1}, {"PYRL": 1, "ACETOL": 1, "FRAG_C": 1}, "k_mgo_pro",
        "B24, AMENDED BY B24b (2026-09-10). methylglyoxal + proline -> 1-pyrroline + HYDROXYACETONE + CO2 "
        "(Strecker of a secondary amine; the ring nitrogen stays in the pyrroline). B24 routed the "
        "hydroxyacetone to the fragment pool, which is why its arm had no competing branch; it is a species "
        "now. Carbon closes as 3 + 5 = 4 + 3 + 1. FITTED on Hofmann & Schieberle 1998b Table 9.",
    ),
    Reaction(
        "r_pyrl_ha_athp", {"PYRL": 1, "ACETOL": 1}, {"ATHP": 1}, "k_ha_athp",
        "B24b. 1-pyrroline + hydroxyacetone -> 2-acetyltetrahydropyridine. THE BRANCH IS EXCLUSIVE AND THE "
        "SOURCE SAYS SO IN WORDS: Schieberle & Hofmann 2005 state that hydroxyacetone gives only this product "
        "and methylglyoxal only 2-acetyl-1-pyrroline, so the two are not competing rates on one substrate -- "
        "they are competing claims on the METHYLGLYOXAL. pH-gated on Schieberle & Hofmann 2005 Table 2 "
        "(<0.1 / 0.9 / 10.8 / 38.4 ug at pH 3 / 5 / 7 / 9). Carbon closes as 4 + 3 = 7.",
    ),
    Reaction(
        "r_pyrl_loss", {"PYRL": 1}, {"MEL_C": 4, "MEL_N": 1}, "k_pyrl_loss",
        "B24b. 1-pyrroline's own loss, first order. B24 had none, and with 1-pyrroline fed in fivefold "
        "excess it made 41 mol % of the methylglyoxal into the product against a printed 0.33 -- 2.1 decades "
        "out. Sized on that experiment (Hofmann & Schieberle 1998b Table 7, experiment 3). The lost "
        "1-pyrroline goes to the MELANOIDIN pools and not to the fragment pool, because it carries a "
        "NITROGEN and the fragment pool holds carbon only; a nitrogen-bearing residue that browning "
        "does not account for would leave the nitrogen balance open. What this does NOT claim is that "
        "the loss is browning: the source measures a disappearance and names no product, so this is "
        "the accounting sink that keeps both atoms, not a mechanism.",
    ),
    Reaction(
        "r_pyrl_ap", {"PYRL": 1, "MGO": 1}, {"AP": 1, "FRAG_C": 1}, "k_pyrl_ap",
        "B24. 1-pyrroline + methylglyoxal -> 2-acetyl-1-pyrroline + CO2 (acylation at C-2, oxidation in air; "
        "net). FITTED on Hofmann & Schieberle 1998b Table 7's fed-pyrroline yields.",
    ),
)
TRUNK_REACTIONS: Tuple[Reaction, ...] = (REACTIONS + DICARBONYL_REACTIONS + PYRAZINE_REACTIONS + GLYCATION_REACTIONS
                                         + AQUEOUS_GLYOXAL_REACTIONS + METHIONINE_REACTIONS + PROLINE_REACTIONS)
TRUNK_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in TRUNK_REACTIONS)
DICARBONYL_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in DICARBONYL_REACTIONS)
PYRAZINE_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in PYRAZINE_REACTIONS)
GLYCATION_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in GLYCATION_REACTIONS)
AQUEOUS_GLYOXAL_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in AQUEOUS_GLYOXAL_REACTIONS)
METHIONINE_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in METHIONINE_REACTIONS)
PROLINE_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in PROLINE_REACTIONS)

#: Build Wave B7's eleven steps, named so a report can say which part of the
#: trunk is B1's and which is B7's without counting.
FURANIC_REACTION_KEYS: Tuple[str, ...] = (
    "r_glc_tdg", "r_tdg_ddg", "r_ddg_hmf", "r_fru_int", "r_int_hmf",
    "r_fru_odg", "r_tdg_mgo", "r_hmf_self", "r_odg_af", "r_af_dmhf",
    "r_mgo_dmhf",
)

#: B1's fifteen, i.e. everything the B1 fit report was fitted against.
B1_REACTION_KEYS: Tuple[str, ...] = tuple(
    key for key in REACTION_KEYS if key not in FURANIC_REACTION_KEYS
)


# ---------------------------------------------------------------------------
# Construction-time validation
# ---------------------------------------------------------------------------


#: The elements the balance checker enforces. SULFUR was added by Build Wave B2
#: (the sulfur module). Every B1 trunk species carries ``sulfur = 0``, so the
#: extra element is a no-op on the trunk and a real constraint on the sulfur
#: network in ``sulfur.py``, which calls this same function.
BALANCED_ELEMENTS: Tuple[str, ...] = ("carbon", "nitrogen", "sulfur")


def validate_balance(reactions: Sequence[Reaction] = None) -> None:
    reactions = TRUNK_REACTIONS if reactions is None else reactions
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


# ---------------------------------------------------------------------------
# B2.3: the trunk's half of the CENTRE LEDGER
# ---------------------------------------------------------------------------
# See ph_state.validate_charge_closure. The trunk owns five steps that move a
# titratable centre and every one of them is a Martins-measured acid step; the
# trunk's own amine (glycine, and the Schiff base and Amadori compound that
# carry its carboxyl) is invisible to the charge balance, which is a DECLARED
# GAP in ph_state.UNTRACKED_TITRATABLE rather than a licence -- the trunk lane
# has no pH state, so the gap costs no prediction today and becomes a defect
# the day it gets one.
#
# Both extended networks REUSE this table, so a trunk step can only be declared
# in one place: sulfur.CENTRE_LEDGER and acrylamide.ACRYLAMIDE_CENTRE_LEDGER
# both start from it.

TRUNK_CENTRE_LEDGER: Mapping[str, Mapping[str, object]] = {
    "r_tdg_fa": {"carboxyl": +1, "basis": (
        "Martins 2005 step 5: 3-deoxyglucosone -> FORMIC ACID + C5 residue. A "
        "neutral acid formed from a neutral deoxyosone; nothing is titrated "
        "into existence, the molecule simply now has a dissociable proton.")},
    "r_odg_aa": {"carboxyl": +1, "basis": (
        "Martins 2005 step 8: 1-deoxyglucosone -> ACETIC ACID + C4 residue. "
        "As r_tdg_fa.")},
    "r_fru_acids": {"carboxyl": +2, "basis": (
        "Martins 2005 step 10: fructose -> formic + acetic + C3 residue. TWO "
        "neutral acids from one neutral sugar.")},
    "r_fa_frag": {"carboxyl": -1, "basis": (
        "FITTED decomposition of formic acid to unassigned fragment carbon "
        "(decarbonylation / decarboxylation). The acid group is genuinely "
        "destroyed and the step genuinely consumes a proton equivalent. "
        "TRUNK-LANE ONLY: the dynamic pH state is sulfur-lane only, so no "
        "scored pH observable depends on this today.")},
    "r_aa_frag": {"carboxyl": -1, "basis": (
        "FITTED decomposition of acetic acid to unassigned fragment carbon, "
        "exactly as r_fa_frag. The acid group is genuinely destroyed, the step "
        "genuinely consumes a proton equivalent, and the same trunk-lane "
        "caveat applies: no scored pH observable depends on it today.")},
}


def _validate_trunk_charge_closure() -> None:
    """Deferred so that ``network`` need not import ``ph_state`` at module top."""
    from .ph_state import validate_charge_closure

    validate_charge_closure(TRUNK_REACTIONS, TRUNK_CENTRE_LEDGER)


_validate_trunk_charge_closure()


def stoichiometric_matrix(
    reactions: Sequence[Reaction] = None,
) -> np.ndarray:
    """(n_species, n_reactions) net stoichiometry."""
    reactions = TRUNK_REACTIONS if reactions is None else reactions
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
    for reaction in TRUNK_REACTIONS:
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
    rates = np.empty(len(TRUNK_REACTIONS), dtype=float)
    for j, reaction in enumerate(TRUNK_REACTIONS):
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
