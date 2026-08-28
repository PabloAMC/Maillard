"""
src/kinetic_core/sulfur.py

THE SULFUR REACTION NETWORK: FORMATION *AND* CONSUMPTION OF THE THIOL AROMA
COMPOUNDS, AS A MASS-ACTION EXTENSION OF BUILD WAVE B1's CORE.
==========================================================================

Modules 1 and 2 of ``docs/reference/FIT_HOLDOUT_DECLARATION.md``. The network
is B1's fifteen trunk steps PLUS the thirty-four steps below, over B1's
thirteen species PLUS the twenty-eight of ``species_sulfur.py``. Carbon,
nitrogen AND SULFUR balance is enforced at import: the module refuses to load
if any step fails.

THE FOUR ARCHITECTURAL COMMITMENTS, each a measured finding
-----------------------------------------------------------
1. **NO FIXED BRANCH FRACTIONS, ANYWHERE.** Cerny 2007 Table 5 measures a 2x
   change in precursor loading moving the xylose share of MFT from 15% to 46% --
   one lab, one method, one pH, one temperature. Every split this module can be
   asked about (thiamine vs sugar, dimer vs monomer, norfuraneol vs intact
   skeleton, FFT vs the rest of the furfural flux) is computed by
   ``branch_shares`` as a RATIO OF TIME-INTEGRATED MASS-ACTION FLUXES. There is
   no fraction constant to find, because the routes have different reaction
   ORDERS -- the thiamine route is first order in thiamine while the sugar
   route is second order in (deoxyosone x sulfide) -- so a concentration change
   moves the split as a matter of arithmetic. A model with fixed fractions
   fails that hold-out by construction, which is the point of holding it out.

2. **THIOL CONSUMPTION IS FOUR NAMED CHANNELS, NOT ONE ARRHENIUS.** Inventory
   sec. B.7: four papers, four temperatures, four dominant sinks, each
   excluding the others'. The channels are declared in
   ``parameters_sulfur.THIOL_CHANNELS`` with their own temperature windows and
   ``ea_kj_mol = None``. Deriving an Ea by pairing the 30 C and 115 C rates is
   a PROHIBITED DERIVATION (sec. C.1) and no code path reaches it.

3. **DIMERISATION IS NOT AROMA LOSS.** The MFT dimer is 15.6x more potent than
   its own monomer, so the ~7-10% of MFT-equivalents that sit in it carry an
   OAV that matches the monomer's. It is a tracked SPECIES with its own
   threshold (``species_sulfur.odour_activity_values``), never a loss term.

4. **pH IS STRUCTURAL, NOT FITTED.** Acid-catalysed 1,2-enolisation, base-
   catalysed 2,3-enolisation, the measured neutral-H2S speciation and Zheng &
   Ho's measured four-pH cysteine thermolysis set. Zero fitted pH parameters.
   See ``parameters_sulfur.PH_TERM_PROVENANCE``.

WHERE THE TRUNK FEEDS THIS MODULE, AND WHERE IT DOES NOT
--------------------------------------------------------
B1's fifteen steps run unchanged inside this network and its methylglyoxal,
glucose and fructose pools are shared. But B1's trunk is GLUCOSE/GLYCINE
specific: its Amadori route needs glycine, and the sulfur systems' amine is
CYSTEINE. In a ribose/cysteine pot B1's steps are therefore near-silent, and
the pentose/cysteine trunk steps this module adds (``r_pent_dpo`` through
``r_dpo_c2c3``) are what actually run. Both are present; neither is duplicated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from .network import BALANCED_ELEMENTS, REACTIONS as TRUNK_REACTIONS, Reaction
from .parameters import KineticParameter
from .parameters_sulfur import (
    MEASURED_SULFUR,
    THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL,
    SulfurParameter,
    T_REF_S_K,
    cysteine_thermolysis_k,
    oligomerisation_rate,
    ph_factor,
    sulfur_registry_metadata,
)
from .species_sulfur import (
    N_SULFUR_STATE,
    SITE_POOLS,
    SULFUR_BY_KEY,
    SULFUR_INDEX,
    SULFUR_STATE_KEYS,
    TERMINAL_POOLS,
    initial_sulfur_state,
    odour_activity_values,
    total_element,
)

# ---------------------------------------------------------------------------
# OUT OF SCOPE for Build Wave B2 -- declared, with what each strands
# ---------------------------------------------------------------------------
OUT_OF_SCOPE: Tuple[Dict[str, str], ...] = (
    {
        "lane": "pyrazines",
        "what_is_stranded": (
            "Amrani-Hemaimi 1995 Table 2's 40 isotope-fraction cells (a FIT "
            "row) and the alanine-vs-glycine ethyl-pyrazine ON/OFF switch (a "
            "star HOLD-OUT), plus Zhou 2023 Table 1's pyrazine and "
            "methylpyrazine columns."
        ),
        "why": (
            "Pyrazines are a nitrogen-heterocycle lane with its own dicarbonyl "
            "and amino-acid bookkeeping; nothing in it is a thiol. Adding it "
            "would double the network for no gain on Modules 1-2. Deferred to a "
            "later wave, EXPLICITLY -- the on/off switch hold-out is NOT scored "
            "below and is not claimed."
        ),
    },
    {
        "lane": "Sotolon and 2-acetyl-2-thiazoline",
        "what_is_stranded": (
            "Hofmann 1996b's Sotolon anchors (13.5 / 764.7 / 273.1 ug, a FIT "
            "row) and the oxidant series (a star HOLD-OUT)."
        ),
        "why": (
            "Sotolon's precursors (butane-2,3-dione + hydroxyacetaldehyde) and "
            "2-acetyl-2-thiazoline's (HDT/ATD) are a separate sub-network, and "
            "the AT numbers are figure-derived LOWER BOUNDS whose Ea inequality "
            "(sec. B10.4) is explicitly licensed as a MODEL CHANGE and not a "
            "parameter transfer. Deferred."
        ),
    },
    {
        "lane": "butane-2,3-dione and the Cerny 2004 in-situ ratios",
        "what_is_stranded": (
            "Cerny 2004's 54:46 butanedione split, 65:35 thiazole and 87:13 "
            "methylthio splits (FIT rows)."
        ),
        "why": (
            "Butane-2,3-dione is not in the state vector. The one Cerny 2004 "
            "constraint this module DOES carry is structural rather than "
            "numeric: 'the C2+C3 route was not relevant at 95 C', which is why "
            "that route competes rather than being assigned a share."
        ),
    },
    {
        "lane": "the Bornhorst norfuraneol Ea and the whole alkaline block",
        "what_is_stranded": "Bornhorst 2017/2017b Ea, matrix pair and structural zero.",
        "why": (
            "pH 8.4-9.5 against this module's pH 4.5-7 -- the alkaline-pH wall, "
            "sec. C.13, 'neither the rates nor the M-2 values transfer to a "
            "pH-5 sulfur benchmark'. Carried as an UNVALIDATED PRIOR in "
            "parameters_sulfur.ALKALINE_PRIORS and operative nowhere."
        ),
    },
    {
        "lane": "ribose vs xylose",
        "what_is_stranded": "the 1.38x MFT / 1.26x FFT gap between them at pH 5.",
        "why": (
            "One generic aldopentose. No step in the corpus is measured "
            "separately for the two sugars, so a per-sugar constant would be "
            "fitted to exactly the two numbers it then reproduces. The gap is "
            "reported as a fit residual, not absorbed."
        ),
    },
)


# ---------------------------------------------------------------------------
# THE NETWORK
# ---------------------------------------------------------------------------

SULFUR_REACTIONS: Tuple[Reaction, ...] = (
    # ================= SULFUR SUPPLY =====================================
    Reaction(
        "r_cys_h2s", {"Cys": 1}, {"H2S": 1, "FRAG_C": 3, "FRAG_N": 1}, "k_cys_h2s",
        "MEASURED, and the only pH-resolved rate in the module. Zheng & Ho 1994, "
        "first-order H2S release from aqueous cysteine at pH 3/5/7/9, 80-110 C, "
        "with the CORRECTED prefactors (9.8e11 / 1.9e12 / 2.4e13 / 1.0e12 1/s), "
        "not the repo's refuted 1.0e14. The C3 and N1 that leave with the "
        "sulfide are unmeasured and are routed to FRAG_C / FRAG_N.",
    ),
    Reaction(
        "r_h2s_loss", {"H2S": 1}, {"FRAG_S": 1}, "k_h2s_loss",
        "FITTED consumption of the sulfide pool -- volatilisation and every "
        "non-thiol sink. Without it H2S accumulates forever, which is the "
        "'nothing is ever consumed' defect this rebuild exists to remove.",
    ),
    # ================= THE PENTOSE / CYSTEINE TRUNK ======================
    Reaction(
        "r_pent_dpo", {"PENT": 1, "Cys": 1}, {"DPO": 1, "Cys": 1}, "k_pent_dpo",
        "2,3-enolisation to the 1-deoxypentosone. Amine-CATALYSED: the amine is "
        "regenerated, and the free-amine fraction supplies the BASE-catalysis "
        "pH factor.",
    ),
    Reaction(
        "r_pent_tdp", {"PENT": 1, "Cys": 1}, {"TDP": 1, "Cys": 1}, "k_pent_tdp",
        "1,2-enolisation to the 3-deoxypentosone. ACID-catalysed, which is what "
        "makes furfural and hence FFT fall monotonically with pH in all three "
        "papers that measure it.",
    ),
    Reaction(
        "r_pent_caramel", {"PENT": 1}, {"TDP": 1}, "k_pent_caramel",
        "THE AMINE-INDEPENDENT (caramelisation) 1,2-enolisation, and it is not "
        "optional. Hofmann 1998 T3/T6 FEED ribose + H2S with NO AMINE AT ALL "
        "and still measure FFT (0.008 mol%) and MFT (0.01 mol%). A network in "
        "which every enolisation is amine-catalysed predicts exactly zero for "
        "those two rows -- which is what the first fit attempt did, and it is "
        "how this step was found. Acid-catalysed, so it carries the same pH "
        "factor as its amine-catalysed twin. Ajandouz 2008 sec. 3.4 sizes the "
        "lane independently at 25-80% of A294 and 7-55% of A420, with the share "
        "RISING with temperature.",
    ),
    Reaction(
        "r_pent_thermal_dpo", {"PENT": 1}, {"DPO": 1}, "k_pent_thermal",
        "The amine-independent 2,3-enolisation, the partner of r_pent_caramel "
        "and required by the same two fed-precursor rows. Bounded to allow "
        "~zero, so the data may reject it.",
    ),
    Reaction(
        "r_arp_dpo", {"ARP": 1}, {"DPO": 1, "FRAG_C": 3, "FRAG_N": 1}, "k_arp_dpo",
        "The fed Amadori's 2,3-enolisation. Base-catalysed. Zhou 2023's feed "
        "partitions between this and r_arp_tdp, which is why MFT (this branch) "
        "and FFT (the acid branch) move in OPPOSITE pH directions in one pot.",
    ),
    Reaction("r_arp_tdp", {"ARP": 1}, {"TDP": 1, "FRAG_C": 3, "FRAG_N": 1},
             "k_arp_tdp", "The fed Amadori's 1,2-enolisation. Acid-catalysed."),
    Reaction("r_dpo_nf", {"DPO": 1}, {"NF": 1}, "k_dpo_nf", ""),
    Reaction(
        "r_dpo_ptr", {"DPO": 1}, {"PTR": 1}, "k_dpo_ptr",
        "Nedvidek 1992 Scheme 2's OTHER partition product, verified by two "
        "negative controls. It exists so that funnelling the whole "
        "1-deoxyosone pool into norfuraneol -- which over-supplies everything "
        "downstream -- is impossible by construction.",
    ),
    Reaction(
        "r_dpo_ddp", {"DPO": 1, "Cys": 1},
        {"DDP": 1, "FRAG_C": 3, "FRAG_N": 1, "FRAG_S": 1}, "k_dpo_ddp",
        "Nedvidek 1992 Scheme 3: 1-deoxyosone + alpha-amino acid -> "
        "1,4-dideoxyosone + RCHO + CO2 + NH3, balance verified C7H13NO6 both "
        "sides for glycine. THE AMINO ACID IS A REAGENT: this step CONSUMES "
        "cysteine, putting the sulfur donor and the carbon skeleton in "
        "competition for one pool. The cysteine's own sulfur leaves in the "
        "Strecker aldehyde (mercaptoacetaldehyde) and is routed to FRAG_S "
        "rather than being silently promoted to H2S.",
    ),
    Reaction(
        "r_dpo_c2c3", {"DPO": 1}, {"HA": 1, "MGO": 1}, "k_dpo_c2c3",
        "The C2+C3 retro-aldol split, as ONE step producing BOTH fragments so "
        "they arrive in the 1:1 stoichiometry the route needs.",
    ),
    Reaction("r_tdp_fur", {"TDP": 1}, {"FUR": 1}, "k_tdp_fur", ""),
    Reaction(
        "r_glc_c2c3", {"Glc": 1}, {"HA": 1, "MGO": 1, "FRAG_C": 1}, "k_glc_ha",
        "THE HEXOSE ENTRY, and the structural reason pentose beats hexose. "
        "Hexoses reach MFT ONLY through fragmentation to C2+C3; pentoses reach "
        "it through fragmentation AND through the intact C5 skeleton. The "
        "measured 10.4x MFT advantage of ribose over glucose (Hofmann 1998 T1) "
        "therefore emerges from having one extra route, with no fitted "
        "sugar-reactivity factor anywhere in the module.",
    ),
    Reaction(
        "r_glc_fur", {"Glc": 1}, {"FUR": 1, "FRAG_C": 1}, "k_glc_fur",
        "THE HEXOSE ROUTE TO FURFURAL, and it is also not optional. Hofmann "
        "1998 T1 measures FFT at 28 ppb from glucose and 32 ppb from fructose, "
        "and FFT's only precursor in this network is furfural, which is a C5. A "
        "hexose therefore has to shed a carbon to reach it. Without this step "
        "the model predicts EXACTLY ZERO FFT from either hexose -- which is "
        "what the first fit attempt did, and it is how this step was found. "
        "Bounded to allow ~zero, so the data may still reject it.",
    ),
    Reaction("r_arp_decay", {"ARP": 1}, {"FRAG_C": 8, "FRAG_N": 1}, "k_osone_decay", ""),
    Reaction("r_osone_decay_dpo", {"DPO": 1}, {"FRAG_C": 5}, "k_osone_decay", ""),
    Reaction("r_osone_decay_tdp", {"TDP": 1}, {"FRAG_C": 5}, "k_osone_decay", ""),
    Reaction("r_osone_decay_ddp", {"DDP": 1}, {"FRAG_C": 5}, "k_osone_decay", ""),
    Reaction("r_ptr_decay", {"PTR": 1}, {"FRAG_C": 5}, "k_osone_decay", ""),
    # ================= THIOL FORMATION ====================================
    Reaction(
        "r_ddp_mft", {"DDP": 1, "H2S": 1}, {"MFT": 1}, "k_ddp_mft",
        "THE INTACT-SKELETON ROUTE. Cerny 2003 T2: 49% unlabelled / 46% 13C5 "
        "with no fragment pattern, 'pathways via ribose fragmentation were not "
        "relevant' => ~93% of MFT carries the intact pentose skeleton at 95 C.",
    ),
    Reaction("r_nf_mft", {"NF": 1, "H2S": 1}, {"MFT": 1}, "k_nf_mft", ""),
    Reaction(
        "r_nf_mp3p", {"NF": 1, "H2S": 1}, {"MP3P": 1}, "k_nf_mp3p",
        "Norfuraneol's DOMINANT fate: mercaptoketones : MFT = 16.3 : 1 from fed "
        "NF, MFT only 2.6% of everything fed NF produces (Whitfield 1999). "
        "2-mercapto-3-pentanone is 96% unlabelled from Cerny's NF spike -- the "
        "sharpest single-species NF-route marker in the corpus.",
    ),
    Reaction(
        "r_ha_mp_mft", {"HA": 1, "MP": 1}, {"MFT": 1}, "k_ha_mp_mft",
        "THE C2+C3 ROUTE. Hofmann 1998 T10's single most effective measured MFT "
        "route (268.1 ug, 0.24 mol% aqueous; 1553.9 ug, 1.39 mol% dry).",
    ),
    Reaction(
        "r_mgo_mp", {"MGO": 1, "H2S": 1}, {"MP": 1}, "k_mgo_mp",
        "Hofmann 1998 T7. 1.8 mol% at 1:1 H2S and 4.0 mol% at 1:2 -- doubling "
        "H2S gives 2.2x product, which is SUPER-linear and which a first-order-"
        "in-H2S mass action law under-predicts. The residual is reported rather "
        "than absorbed into an invented reaction order.",
    ),
    Reaction(
        "r_fur_fft", {"FUR": 1, "H2S": 1}, {"FFT": 1}, "k_fur_fft",
        "Furfural is 60x ribose and 7x the 3-deoxyosone as an FFT source "
        "(Hofmann 1998 T3).",
    ),
    Reaction(
        "r_fur_decay", {"FUR": 1}, {"FRAG_C": 5}, "k_fur_decay",
        "THE LARGE UNIDENTIFIED FURFURAL SINK. Yaghmur 2005: the FFT branch is "
        "<=1.2% of the furfural flux, so nearly all the furfural that "
        "disappears goes somewhere the corpus never identifies. That somewhere "
        "is this step, and the 1.2% ceiling is a declared FIT row that "
        "constrains the ratio between the two.",
    ),
    Reaction("r_nf_decay", {"NF": 1}, {"FRAG_C": 5}, "k_nf_decay", ""),
    Reaction("r_mp_decay", {"MP": 1}, {"FRAG_C": 3, "FRAG_S": 1}, "k_thiol_decay", ""),
    Reaction("r_mp3p_decay", {"MP3P": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_thiol_decay", ""),
    Reaction("r_mp2p_decay", {"MP2P": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_thiol_decay", ""),
    Reaction("r_hmp_decay", {"HMP": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_thiol_decay", ""),
    # ================= THE THIAMINE ROUTE =================================
    Reaction(
        "r_thi_hmp", {"THI": 1}, {"HMP": 1, "FRAG_C": 7, "FRAG_N": 4}, "k_thi_hmp",
        "Thiamine's cleavage to 5-hydroxy-3-mercapto-2-pentanone, the "
        "thiamine -> MFT intermediate. FIRST ORDER IN THIAMINE, against the "
        "sugar route's second order in (deoxyosone x sulfide) -- which is "
        "exactly why the branch fraction MOVES when the loading changes.",
    ),
    Reaction("r_thi_mesh", {"THI": 1}, {"MESH": 1, "FRAG_C": 11, "FRAG_N": 4},
             "k_thi_mesh",
             "The methanethiol source in Zhang's Fig. 1 system, which contains "
             "no methionine."),
    Reaction("r_hmp_mft", {"HMP": 1}, {"MFT": 1}, "k_hmp_mft", ""),
    Reaction(
        "r_hmp_mp2p", {"HMP": 1}, {"MP2P": 1}, "k_hmp_mp2p",
        "3-mercapto-2-pentanone, the THIAMINE-diagnostic isomer: Cerny 2007 T4 "
        "has it 77-90% thiamine-derived while its isomer 2-mercapto-3-pentanone "
        "is 94->95% xylose-derived, in the same pot. The pair is an isomer-split "
        "diagnostic no lumped 'mercaptopentanone' node can express.",
    ),
    # ================= 2-ACETYLTHIAZOLE ===================================
    Reaction(
        "r_cys_actz", {"Cys": 1, "MGO": 1}, {"ACTZ": 1, "FRAG_C": 1}, "k_cys_actz",
        "Zhou 2023 measures 2-acetylthiazole in BOTH an ARP-fed and an MGO-fed "
        "system at pH 8 (582 vs 665 ug/L). That pair is a step-level, "
        "two-system constraint on the ARP -> MGO flux with no free downstream "
        "parameter, and it is a declared HOLD-OUT.",
    ),
    # ================= THIOL CONSUMPTION: FOUR NAMED CHANNELS =============
    # Channel 1 -- covalent addition to matrix electrophiles. 25-30 C. MEASURED.
    Reaction(
        "ch_thioether_mft", {"MFT": 1, "MELE": 1}, {"BND": 1}, "k_thioether",
        "CHANNEL 1, dominant at 25-30 C. Bimolecular on a DEPLETABLE site pool, "
        "which is what lets the model express Hofmann's 80 C brew being SLOWER "
        "than his own 30 C models. Disulfide branch measured at <1.5% of the "
        "flux at 30 C, which fixes the stoichiometry as 1:1 and the order as "
        "first in thiol.",
    ),
    Reaction("ch_thioether_fft", {"FFT": 1, "MELE": 1}, {"BND": 1}, "k_thioether", ""),
    Reaction(
        "ch_thioether_release", {"BND": 1}, {"MFT": 1, "MELE": 1}, None,
        "THE CONJUGATION IS REVERSIBLE, and leaving it out would have been a "
        "real error. Stack 2018 measures thiol-quinone conjugation as an "
        "EQUILIBRIUM (K ~ 1e2-1e3 M^-1) -- a SEQUESTRATION equilibrium, not an "
        "irreversible sink -- which is the mechanism behind Sun 2019's finding "
        "that ADDING cysteine RELEASES bound FFT. An irreversible-only sink "
        "over-predicts permanent aroma loss on long timescales. "
        "ITS RATE IS DERIVED, NOT FITTED: k_reverse = k_thioether / K, so the "
        "step adds no free parameter and cannot be tuned. "
        "TWO LIMITATIONS, STATED: the release is written to MFT because the "
        "bound adduct is a single lumped C5 species that serves both thiols, so "
        "a pot containing both would slowly convert FFT into MFT through it; "
        "the module is therefore only valid for release where one thiol "
        "dominates, and BND is reported alongside so the size of that pool is "
        "visible. And the equilibrium constant is measured on quinone "
        "conjugates, i.e. on one of the two channels MELE lumps.",
    ),
    # Channel 2 -- acid-catalysed C-5 oligomerisation. 50 C. NO MEASURED RATE.
    Reaction(
        "ch_oligomer_mft", {"MFT": 1}, {"OLG": 1}, "k_oligomer",
        "CHANNEL 2, dominant at 50 C: acid-catalysed electrophilic "
        "oligomerisation at C-5 (85% H-D exchange at C-5 against 10% at C-4; "
        "air ~ argon, so NOT oxidative; the MFT mass balance fails to close as "
        "thiol + disulfide while MB, FFT and DMFT all close). ITS RATE IS ZERO "
        "HERE: van Seeventer 2001 Table 1 is the only measurement and it is a "
        "declared HOLD-OUT, so fitting it would be reading the hold-out. The "
        "consequence -- the model predicts no 50 C oligomerisation and will "
        "fail that row -- is pre-registered, not discovered afterwards.",
    ),
    Reaction("ch_oligomer_fft", {"FFT": 1}, {"OLG": 1}, "k_oligomer", ""),
    # Channel 3 -- oxidative dimerisation. 115-120 C.
    Reaction(
        "ch_dimer_mft", {"MFT": 2, "OX": 1}, {"MFTD": 1}, "k_dimer_mft",
        "CHANNEL 3, dominant at 115-120 C -- the channel the 30 C system rules "
        "out at <1.5%. SECOND order in thiol and FIRST in oxidant equivalents, "
        "because Zhang 2024 Fig. 1 shows the branch responds to the additive's "
        "REDOX STATE and not to its concentration: at near-matched molar sulfur "
        "(123.8 vs 124.8 mM S) cystine sends 54.2% of the MFT pool to the dimer "
        "and cysteine only 8.6%. NOT AROMA LOSS -- the dimer is 15.6x more "
        "potent than the monomer.",
    ),
    Reaction("ch_dimer_fft", {"FFT": 2, "OX": 1}, {"FFTD": 1}, "k_dimer_fft", ""),
    # Channel 4 -- radical coupling to methanethiol. 115 C.
    Reaction(
        "ch_mmft", {"MFT": 1, "MESH": 1}, {"MMFT": 1}, "k_mmft",
        "CHANNEL 4: MFT + MeSH radical coupling. Requires a methanethiol "
        "partner, which the 30 C and 50 C systems do not have -- which is why "
        "it is a separate channel and not a term in a lumped sink.",
    ),
    # Terminal decay of the aroma species, bounded to allow ~zero.
    Reaction("ch_mft_decay", {"MFT": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_mft_decay", ""),
    Reaction("ch_fft_decay", {"FFT": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_fft_decay", ""),
    Reaction("ch_mftd_decay", {"MFTD": 1}, {"FRAG_C": 10, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_fftd_decay", {"FFTD": 1}, {"FRAG_C": 10, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_mmft_decay", {"MMFT": 1}, {"FRAG_C": 6, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_mesh_decay", {"MESH": 1}, {"FRAG_C": 1, "FRAG_S": 1}, "k_thiol_decay", ""),
    Reaction("ch_actz_decay", {"ACTZ": 1}, {"FRAG_C": 5, "FRAG_N": 1, "FRAG_S": 1},
             "k_thiol_decay", ""),
)

#: The full network: B1's trunk first, then the sulfur block.
FULL_REACTIONS: Tuple[Reaction, ...] = TRUNK_REACTIONS + SULFUR_REACTIONS
FULL_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in FULL_REACTIONS)

#: Which pH factor multiplies which reaction. Assigning a step to a factor is a
#: MECHANISM claim (acid- or base-catalysed enolisation, sulfide nucleophile)
#: and carries no fitted number; see parameters_sulfur.ph_factor.
REACTION_PH_FACTOR: Mapping[str, str] = {
    "r_arp_dpo": "base",
    "r_arp_tdp": "acid",
    "r_pent_dpo": "base",
    "r_dpo_ddp": "base",
    "r_pent_tdp": "acid",
    "r_pent_caramel": "acid",
    "r_glc_fur": "acid",
    "r_ddp_mft": "neutral_h2s",
    "r_nf_mft": "neutral_h2s",
    "r_nf_mp3p": "neutral_h2s",
    "r_mgo_mp": "neutral_h2s",
    "r_fur_fft": "neutral_h2s",
}


# ---------------------------------------------------------------------------
# Construction-time validation, over THREE elements
# ---------------------------------------------------------------------------


def validate_sulfur_balance(
    reactions: Sequence[Reaction] = FULL_REACTIONS,
) -> None:
    """
    Raise unless every step conserves carbon, nitrogen AND sulfur.

    The same arithmetic B1 uses, over the extended species table. The sulfur
    element is the new constraint and it is a real one: it is what forces the
    Nedvidek step to account for the cysteine sulfur it consumes, and what
    stops a thiol being made out of nothing.
    """
    for reaction in reactions:
        for element in BALANCED_ELEMENTS:
            left, right = reaction.atom_balance(element, SULFUR_BY_KEY)
            if left != right:
                raise ValueError(
                    f"{reaction.key}: {element} does not balance "
                    f"({left} -> {right}). Every step must conserve atoms; the "
                    f"unmeasured residue belongs in FRAG_C / FRAG_N / FRAG_S, "
                    f"not in nowhere."
                )
        for key in list(reaction.reactants) + list(reaction.products):
            if key not in SULFUR_INDEX:
                raise ValueError(f"{reaction.key}: unknown species {key!r}")
        for pool in ("FRAG_C",) + TERMINAL_POOLS:
            if pool in reaction.reactants:
                raise ValueError(
                    f"{reaction.key}: {pool} is a terminal accounting pool and "
                    f"must never be a reactant."
                )
        if "MEL_C" in reaction.reactants or "MEL_N" in reaction.reactants:
            raise ValueError(f"{reaction.key}: the melanoidin sink is TERMINAL.")


validate_sulfur_balance()


def sulfur_stoichiometric_matrix(
    reactions: Sequence[Reaction] = FULL_REACTIONS,
) -> np.ndarray:
    matrix = np.zeros((N_SULFUR_STATE, len(reactions)), dtype=float)
    for j, reaction in enumerate(reactions):
        for key, coefficient in reaction.reactants.items():
            matrix[SULFUR_INDEX[key], j] -= float(coefficient)
        for key, coefficient in reaction.products.items():
            matrix[SULFUR_INDEX[key], j] += float(coefficient)
    return matrix


SULFUR_STOICHIOMETRY: np.ndarray = sulfur_stoichiometric_matrix()


# ---------------------------------------------------------------------------
# Rate evaluation
# ---------------------------------------------------------------------------


def sulfur_rate_constants_at(
    parameters: Mapping[str, Any],
    temperature_k: float,
    ph: float,
) -> Dict[str, float]:
    """
    Every reaction's rate constant at ``temperature_k`` and ``ph``.

    THREE THINGS HAPPEN HERE AND NOWHERE ELSE:

      * the CONSUMPTION channels are NOT Arrhenius-scaled, because they have no
        activation energy (policy 2). ``SulfurParameter.k_at`` returns k_ref
        unchanged when ``ea_kj_mol`` is None;
      * the pH factor of a step multiplies its constant. The factor is
        normalised to 1 at pH 5, so a fitted k_ref is always "k at pH 5";
      * cysteine thermolysis is recomputed from Zheng & Ho's MEASURED four-pH
        Arrhenius set rather than being scaled by a pH factor, because for that
        one step the pH dependence is measured directly.
    """
    out: Dict[str, float] = {}
    for reaction in FULL_REACTIONS:
        key = reaction.parameter_key
        if key is None:
            continue
        if key == "k_oligomer":
            out[reaction.key] = oligomerisation_rate()
            continue
        if key == "k_cys_h2s":
            out[reaction.key] = cysteine_thermolysis_k(ph, float(temperature_k))
            continue
        parameter = parameters.get(key)
        if parameter is None:
            raise KeyError(f"{reaction.key}: no parameter {key!r} supplied")
        if parameter.k_ref is None:
            raise ValueError(
                f"{reaction.key}: parameter {key!r} is unpopulated "
                f"(evidence_class={parameter.evidence_class}). The fitted steps "
                f"must be given values before the network can be integrated; "
                f"there is no silent default."
            )
        if isinstance(parameter, KineticParameter) and parameter.ea_kj_mol is None:
            raise ValueError(f"{reaction.key}: trunk parameter {key!r} has no Ea")
        value = parameter.k_at(float(temperature_k))
        value *= ph_factor(REACTION_PH_FACTOR.get(reaction.key, ""), ph)
        out[reaction.key] = value

    # DERIVED RATE 1: the reverse of the thioether conjugation. Stack 2018
    # measures the equilibrium, not the reverse rate, so the reverse constant is
    # k_forward / K and is not a free parameter. Units check: k_thioether is
    # L/(mmol*min) and K is L/mmol, so the quotient is 1/min, which is the
    # first-order unit the release step needs.
    out["ch_thioether_release"] = (
        out["ch_thioether_mft"] / THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL
    )

    # DERIVED RATE 2: B1's Amadori rearrangement, pinned to the condensation.
    # Reproduced here so the trunk behaves identically inside the extended
    # network.
    from .parameters import SCHIFF_AMADORI_SPLIT

    out["r_amadori"] = (
        float(SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
        * out["r_schiff"]
        * float(SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"])
    )
    return out


# Precomputed reactant layout. The right-hand side is called tens of thousands
# of times per integration and hundreds of integrations per fit iteration, so
# the per-call dict iteration and attribute lookup of the naive form dominate
# the wall clock. This flattens the same arithmetic to index/exponent tuples
# built once at import; the RESULT is identical and the unit tests pin that.
_REACTANT_LAYOUT: Tuple[Tuple[Tuple[int, int], ...], ...] = tuple(
    tuple((SULFUR_INDEX[key], int(c)) for key, c in r.reactants.items())
    for r in FULL_REACTIONS
)
_RATE_KEY_ORDER: Tuple[str, ...] = tuple(r.key for r in FULL_REACTIONS)


def rate_vector(k: Mapping[str, float]) -> np.ndarray:
    """The rate constants as a positional array in FULL_REACTIONS order."""
    return np.array([float(k[key]) for key in _RATE_KEY_ORDER], dtype=float)


def sulfur_reaction_rates(state: np.ndarray, k) -> np.ndarray:
    """
    Mass-action rate of every reaction, mmol/(L*min).

    ``k`` may be the {reaction_key: constant} mapping or the positional array
    from ``rate_vector``; the array form skips a dict lookup per reaction per
    call and is what the integrator uses.
    """
    y = np.asarray(state, dtype=float)
    if y.min() < 0.0:
        y = np.clip(y, 0.0, None)
    kv = k if isinstance(k, np.ndarray) else rate_vector(k)
    rates = np.empty(len(_REACTANT_LAYOUT), dtype=float)
    for j, layout in enumerate(_REACTANT_LAYOUT):
        value = kv[j]
        for index, coefficient in layout:
            value = value * (y[index] if coefficient == 1 else y[index] ** coefficient)
        rates[j] = value
    return rates


def sulfur_derivatives(state: np.ndarray, k) -> np.ndarray:
    return SULFUR_STOICHIOMETRY @ sulfur_reaction_rates(state, k)


# ---------------------------------------------------------------------------
# Integration, with an optional MEASURED pH trajectory
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SulfurRun:
    """One integrated sulfur experiment, with its conservation audit attached."""

    times_min: np.ndarray
    concentrations: np.ndarray
    temperature_k: float
    ph_initial: float
    ph_final: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def series(self, species_key: str) -> np.ndarray:
        return self.concentrations[:, SULFUR_INDEX[species_key]]

    def final(self, species_key: str) -> float:
        return float(self.concentrations[-1, SULFUR_INDEX[species_key]])

    def element_closure(self, element: str) -> Dict[str, float]:
        totals = np.array(
            [total_element(row, element) for row in self.concentrations]
        )
        t0 = float(totals[0])
        drift = float(np.max(np.abs(totals - t0)))
        return {
            f"initial_{element}_mmol_L": t0,
            "max_abs_drift_mmol_L": drift,
            "max_relative_drift": drift / t0 if t0 > 0 else float("nan"),
        }

    def odour_activity(self) -> Dict[str, object]:
        return odour_activity_values(self.concentrations[-1])


def _sulfur_extrapolation_warnings(
    parameters: Mapping[str, Any], temperature_c: float, ph: float
) -> list:
    out = []
    for key, parameter in sorted(parameters.items()):
        low, high = parameter.temperature_range_c
        if temperature_c < low - 1e-9 or temperature_c > high + 1e-9:
            no_ea = getattr(parameter, "ea_kj_mol", 0.0) is None
            out.append(
                f"{key}: evaluated at {temperature_c:.1f} C, measured over "
                f"{low:.0f}-{high:.0f} C"
                + (
                    " -- AND IT HAS NO ACTIVATION ENERGY, so the constant is "
                    "HELD FIXED rather than extrapolated. This is a declared "
                    "gap (inventory sec. C.1), not a modelling choice."
                    if no_ea
                    else f" ({parameter.source_anchor[:60]}...)"
                )
            )
        measured_ph = getattr(parameter, "ph_of_measurement", None)
        if measured_ph is not None and abs(float(measured_ph) - float(ph)) > 1.0:
            out.append(
                f"{key}: evaluated at pH {ph:.2f}, measured at pH "
                f"{float(measured_ph):.2f}"
            )
    return out


def integrate_sulfur(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    times_min: Sequence[float],
    *,
    ph: float = 5.0,
    ph_final: Optional[float] = None,
    method: str = "LSODA",
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> SulfurRun:
    """
    Integrate the sulfur network at one temperature and one pH TRAJECTORY.

    ``ph`` is the initial pH. ``ph_final``, when given, is the MEASURED final
    pH: the pH then moves linearly in time between the two and every
    pH-dependent constant is re-evaluated inside the right-hand side.

    That is not a convenience. Zhou 2023's pH labels are INITIAL pH of an
    UNBUFFERED system whose pH-7 run ENDS at 3.42 and whose pH-6 run ends at
    3.22 -- within 0.2 units of each other (inventory sec. C.11, "the most
    important methodological number in the paper"). Both endpoints are measured.
    FIT_HOLDOUT_DECLARATION.md sec.5 decision 1 scores the pH-6 hold-out column
    as DIAGNOSTIC ONLY "until the model carries a pH-trajectory state"; this is
    that state. A buffered system (every Hofmann and Cerny row, 0.5 M
    phosphate) is run with ``ph_final=None``, i.e. clamped.
    """
    from scipy.integrate import solve_ivp

    if method not in ("LSODA", "BDF", "Radau"):
        raise ValueError(f"method {method!r} is not a stiff-capable solver")
    grid = np.asarray(times_min, dtype=float)
    if grid.ndim != 1 or grid.size == 0:
        raise ValueError("times_min must be a non-empty 1-D sequence")
    if np.any(np.diff(grid) < 0):
        raise ValueError("times_min must be non-decreasing")
    if grid[0] < 0:
        raise ValueError("times_min must start at or after 0")

    duration = float(max(grid[-1], 1e-12))
    ph0 = float(ph)
    ph1 = ph0 if ph_final is None else float(ph_final)
    y0 = initial_sulfur_state(initial)

    if ph1 == ph0:
        k_fixed = rate_vector(
            sulfur_rate_constants_at(parameters, float(temperature_k), ph0)
        )

        def rhs(_t: float, y: np.ndarray) -> np.ndarray:
            return sulfur_derivatives(y, k_fixed)
    else:
        # The pH trajectory is smooth and the rate constants depend on it only
        # through three closed-form factors, so the constants are pre-tabulated
        # on a fine grid in pH and interpolated inside the right-hand side.
        # Re-evaluating the whole registry at every solver step is what makes a
        # trajectory run 30x slower than a clamped one; this makes it 1.1x.
        n_grid = 41
        ph_grid = np.linspace(ph0, ph1, n_grid)
        k_grid = np.array([
            rate_vector(sulfur_rate_constants_at(parameters, float(temperature_k), p))
            for p in ph_grid
        ])

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            frac = min(max(float(t) / duration, 0.0), 1.0)
            position = frac * (n_grid - 1)
            lo = min(int(position), n_grid - 2)
            w = position - lo
            k_t = k_grid[lo] * (1.0 - w) + k_grid[lo + 1] * w
            return sulfur_derivatives(y, k_t)

    solution = solve_ivp(
        rhs, (0.0, duration), y0, t_eval=grid, method=method, rtol=rtol, atol=atol
    )
    if not solution.success:
        raise RuntimeError(f"sulfur-network integration failed: {solution.message}")

    raw = solution.y.T
    worst_negative = float(np.min(raw)) if raw.size else 0.0
    if worst_negative < -max(1e3 * atol, 1e-8):
        raise RuntimeError(
            f"sulfur-network integration produced a state of {worst_negative:.3e}, "
            f"far below the absolute tolerance {atol:.1e}: a genuine "
            f"non-negativity failure, not solver noise."
        )
    concentrations = np.clip(raw, 0.0, None)

    temperature_c = float(temperature_k) - 273.15
    metadata: Dict[str, Any] = {
        "method": method,
        "temperature_C": temperature_c,
        "temperature_K": float(temperature_k),
        "ph_initial": ph0,
        "ph_final": ph1,
        "ph_is_clamped": ph1 == ph0,
        "n_species": N_SULFUR_STATE,
        "n_reactions": len(FULL_REACTIONS),
        "species_order": list(SULFUR_STATE_KEYS),
        "worst_raw_negative_before_clip": worst_negative,
        "extrapolation_warnings": _sulfur_extrapolation_warnings(
            parameters, temperature_c, ph0
        ),
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
    }
    metadata.update(sulfur_registry_metadata(
        {k: v for k, v in parameters.items() if isinstance(v, SulfurParameter)}
    ))

    run = SulfurRun(
        times_min=grid,
        concentrations=concentrations,
        temperature_k=float(temperature_k),
        ph_initial=ph0,
        ph_final=ph1,
        metadata=metadata,
    )
    for element in BALANCED_ELEMENTS:
        metadata[f"{element}_closure"] = run.element_closure(element)
    return run


def sulfur_flux_budget(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    duration_min: float,
    *,
    ph: float = 5.0,
    ph_final: Optional[float] = None,
    n_points: int = 201,
) -> Dict[str, float]:
    """Time-integrated flux through every reaction, mmol/L."""
    grid = np.linspace(0.0, float(duration_min), int(n_points))
    run = integrate_sulfur(
        parameters, temperature_k, initial, grid, ph=ph, ph_final=ph_final
    )
    ph0, ph1 = run.ph_initial, run.ph_final
    rows = []
    for i, t in enumerate(grid):
        frac = min(max(t / max(grid[-1], 1e-12), 0.0), 1.0)
        k_t = sulfur_rate_constants_at(
            parameters, float(temperature_k), ph0 + frac * (ph1 - ph0)
        )
        rows.append(sulfur_reaction_rates(run.concentrations[i], k_t))
    integral = np.trapezoid(np.array(rows), grid, axis=0)
    return {r.key: float(integral[j]) for j, r in enumerate(FULL_REACTIONS)}


# ---------------------------------------------------------------------------
# BRANCH SHARES -- computed, never stored
# ---------------------------------------------------------------------------


def branch_shares(flux: Mapping[str, float]) -> Dict[str, float]:
    """
    Every route split this module can be asked about, as a RATIO OF FLUXES.

    Not one of these is a parameter. Each is recomputed from the time-integrated
    mass-action fluxes of a particular run, so each MOVES when the loading,
    the pH or the temperature moves. That is the whole architectural point:
    Cerny 2007 Table 5 measures a 2x precursor change moving the xylose share of
    MFT from 15% to 46%, and a model with a stored fraction cannot reproduce it.
    """
    def share(numerator: float, denominator: float) -> float:
        return float(numerator) / float(denominator) if denominator > 0 else float("nan")

    mft_from_thiamine = flux.get("r_hmp_mft", 0.0)
    mft_from_intact = flux.get("r_ddp_mft", 0.0)
    mft_from_nf = flux.get("r_nf_mft", 0.0)
    mft_from_c2c3 = flux.get("r_ha_mp_mft", 0.0)
    mft_total = mft_from_thiamine + mft_from_intact + mft_from_nf + mft_from_c2c3

    fur_to_fft = flux.get("r_fur_fft", 0.0)
    fur_total = fur_to_fft + flux.get("r_fur_decay", 0.0)

    mft_made = mft_total
    mft_to_dimer = 2.0 * flux.get("ch_dimer_mft", 0.0)
    mft_to_mmft = flux.get("ch_mmft", 0.0)
    mft_to_thioether = flux.get("ch_thioether_mft", 0.0)
    mft_to_oligomer = flux.get("ch_oligomer_mft", 0.0)
    mft_to_decay = flux.get("ch_mft_decay", 0.0)

    return {
        "MFT_share_thiamine_route": share(mft_from_thiamine, mft_total),
        "MFT_share_sugar_routes": share(
            mft_from_intact + mft_from_nf + mft_from_c2c3, mft_total
        ),
        "MFT_share_intact_skeleton_route": share(mft_from_intact, mft_total),
        "MFT_share_norfuraneol_route": share(mft_from_nf, mft_total),
        "MFT_share_C2_plus_C3_route": share(mft_from_c2c3, mft_total),
        "FFT_share_of_furfural_flux": share(fur_to_fft, fur_total),
        "MFT_consumed_share_dimerisation": share(mft_to_dimer, mft_made),
        "MFT_consumed_share_MMFT": share(mft_to_mmft, mft_made),
        "MFT_consumed_share_thioether": share(mft_to_thioether, mft_made),
        "MFT_consumed_share_oligomerisation": share(mft_to_oligomer, mft_made),
        "MFT_consumed_share_unassigned_decay": share(mft_to_decay, mft_made),
        "note": (
            "Every value above is a ratio of time-integrated mass-action fluxes "
            "computed for THIS run. None is stored anywhere. Perturb a "
            "precursor and they all move -- which is the pre-registered "
            "response to Cerny 2007 Table 5."
        ),
    }


def describe_sulfur() -> Dict[str, object]:
    """Machine-readable description of the sulfur network, for the reports."""
    return {
        "n_species": N_SULFUR_STATE,
        "n_reactions": len(FULL_REACTIONS),
        "trunk_reactions": len(TRUNK_REACTIONS),
        "sulfur_reactions": len(SULFUR_REACTIONS),
        "species": [
            {
                "key": s.key,
                "label": s.label,
                "carbon": s.carbon,
                "nitrogen": s.nitrogen,
                "sulfur": s.sulfur,
                "role": s.role,
                "measured_in_fit_corpus": s.measured,
            }
            for s in (SULFUR_BY_KEY[key] for key in SULFUR_STATE_KEYS)
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
                "ph_factor": REACTION_PH_FACTOR.get(r.key, ""),
                "carbon_balance": list(r.atom_balance("carbon", SULFUR_BY_KEY)),
                "nitrogen_balance": list(r.atom_balance("nitrogen", SULFUR_BY_KEY)),
                "sulfur_balance": list(r.atom_balance("sulfur", SULFUR_BY_KEY)),
                "note": r.note,
            }
            for r in FULL_REACTIONS
        ],
        "site_pools_zero_atom": list(SITE_POOLS),
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
        "reference_temperature_K": T_REF_S_K,
        "measured_sulfur_parameters": sorted(MEASURED_SULFUR),
    }
