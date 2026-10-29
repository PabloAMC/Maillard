"""
src/kinetic_core/sulfur.py

THE SULFUR REACTION NETWORK: FORMATION *AND* CONSUMPTION OF THE THIOL AROMA
COMPOUNDS, AS A MASS-ACTION EXTENSION OF BUILD WAVE B1's CORE
(Build Wave B2, 2026-08-28; refitted through B2.1-B2.4).
==========================================================================

WAVE: B2 built this lane; B2.1-B2.4 refitted it (B2.2 made pH a computed state
rather than a prescribed one, B2.3 fixed charge conservation, B2.4 declared the
objective weighting and shipped the 18-member ensemble). B7 did NOT touch it and
its residuals moved anyway, because four of B7's furanic steps drain trunk
species this lane reads -- disclosed at
tests/unit/test_kinetic_core_b2_4.py rather than absorbed by a refit.
EXAM: results/validation/kinetic_core_b2_4_{fit,holdout}_report.md is the
current record; the cutover final exam scores it out of sample.
DECLARED GAPS: the lane has ZERO ABSOLUTE LITERATURE ANCHORS -- Wave S2c
retracted the last one when its MFT 342 / FFT 200 ppb targets were traced to a
repo-internal reconstruction rather than to the paper. Its RATIOS are therefore
first-class and its ABSOLUTES are not. See ``parameters_sulfur``'s own
``PROHIBITED_DERIVATIONS`` for the named fits that must never be attempted, and
the ``furanone_reductive_sulfhydrylation`` entry in ``src/barrier_constants.py``
for a fit that was RUN and REJECTED by its own pre-registered isotope gate.

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

2. **THIOL CONSUMPTION IS NAMED CHANNELS, NOT ONE ARRHENIUS.** Inventory
   sec. B.7: four papers, four temperatures, four dominant sinks, each
   excluding the others'. The channels are declared in
   ``parameters_sulfur.THIOL_CHANNELS`` with their own temperature windows and
   ``ea_kj_mol = None``. Deriving an Ea by pairing the 30 C and 115 C rates is
   a PROHIBITED DERIVATION (sec. C.1) and no code path reaches it.
   *B2.1 narrows this to the NAMED channels it was always about.* The residual
   ``*_decay`` lumps have never been measured at any temperature, so there is
   no cross-lab pairing to prohibit -- and holding them fixed at 145 C is the
   strong claim, not the weak one. See ``NAMED_CHANNEL_KEYS`` versus
   ``UNASSIGNED_SINK_KEYS``.

3. **DIMERISATION IS NOT AROMA LOSS.** The MFT dimer is 15.6x more potent than
   its own monomer, so the ~7-10% of MFT-equivalents that sit in it carry an
   OAV that matches the monomer's. It is a tracked SPECIES with its own
   threshold (``species_sulfur.odour_activity_values``), never a loss term.

4. **pH IS STRUCTURAL, NOT FITTED.** Acid- and base-catalysed enolisation, the
   measured sulfide and amine speciation, and Zheng & Ho's measured four-pH
   cysteine thermolysis set. Zero fitted pH parameters -- still zero after
   B2.1. See ``parameters_sulfur.PH_TERM_PROVENANCE``.

WHAT BUILD WAVE B2.1 CHANGED, AND WHY
-------------------------------------
B2's pre-registered scorecard was 7/27 gating with a NAMED dominant defect: its
structural pH slope was about one decade per pH unit, which gave the right
SHAPES and far too much SLOPE, so every level three pH units from the fitted pH
collapsed or exploded. Five changes, each traceable to a source:

* **pH, part 1 -- speciation at reaction temperature.** Every pKa is now
  evaluated by van 't Hoff at the temperature the reaction runs at. B2 used
  25 C values at 115-145 C and called it a declared approximation; the cysteine
  alpha-ammonium moves 10.28 -> 8.09 and the sulfide 7.05 -> 5.94, i.e. from
  outside the operating window to inside it.
* **pH, part 2 -- no lane is 100% catalysed.** Each pH-sensitive step now runs
  in parallel with an uncatalysed or oppositely-charged twin
  (``r_arp_tdp_thermal``, ``r_arp_dpo_thermal``, ``r_ddp_mft_hs``,
  ``r_fur_fft_hs``, and ``r_pent_caramel`` with its acid tag removed). The
  observable slope is now a ratio of two ordinary rate constants fitted on FIT
  rows -- not a pH exponent.
* **pH, part 3 -- consumption carries pH too.** Kumazawa 2003's seven-point
  121 C grid (declared FIT) measures 2-furfurylthiol survival collapsing from
  99.5% to under 0.5% between pH 3 and 7 with no formation involved at all.
  B2 had to attribute all of that to formation. ``ch_thiolate_loss_*`` and the
  thiolate tag on the dimerisation channels carry it where it belongs.
* **Temperature.** Kang 2026 SI Table S4's 100 and 120 C rungs are FIT, giving
  the module its first in-system temperature axis (B2's single lumped Ea landed
  on its bound and was reported as a failure to determine). Kang's measured
  Ea(free-cysteine depletion) = 55.1 kJ/mol is a FIXED barrier on the new
  ``r_cys_thermal``. No new global sulfur Ea was added.
* **Stack 2018, corrected.** Both the forward and the reverse constant are
  measured, so K = 5.64 M^-1 at 19.4 C and dH = -28.5 kJ/mol. B2 used 316 M^-1
  with no temperature dependence, on this paper's stated authority; the paper
  does not support it. And BND is split into BND and BND_F, which retires B2's
  own recorded limitation that the shared reservoir converted FFT into MFT.
* **A new channel.** ``ch_protein_ss_*``: FFT is consumed by protein
  disulfides (Amendment 5). Rate-bounded rather than fitted, and identically
  zero in every scored row because no system in either panel contains protein.

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

from .network import (
    BALANCED_ELEMENTS,
    REACTIONS as TRUNK_REACTIONS,
    TRUNK_CENTRE_LEDGER,
    Reaction,
)
from .parameters import KineticParameter
from .parameters_furanic import TEMPERATURE_CAP_K as _FURANIC_TEMPERATURE_CAP_K
from .parameters_sulfur import (
    MEASURED_SULFUR,
    SulfurParameter,
    T_REF_S_K,
    cysteine_thermolysis_k,
    glutathione_thermolysis_k,
    thiamine_h2s_share,
    oligomerisation_rate,
    ph_factor,
    sulfur_registry_metadata,
    thiol_adduct_equilibrium_l_per_mmol,
)
from .ph_state import (
    AMINE_FATE_BASIS,
    BUFFER_ABSENT_WARNING,
    DEFAULT_BUFFER,
    LABEL_TEMPERATURE_K,
    BufferSpec,
    PhDrift,
    buffer_capacity_mmol_per_ph,
    centre_delta,
    ph_and_ionic_strength_of_state,
    ph_of_state,
    strong_ion_difference,
    validate_charge_closure,
)
from . import ph_state as _ph
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
    # 2026-10-27: glutathione releases H2S too (Zheng & Ho 1994 Table III, measured, pH-resolved).
    Reaction(
        "r_gsh_h2s", {"GSH": 1}, {"H2S": 1, "FRAG_C": 10, "FRAG_N": 3}, "k_gsh_h2s",
        "Glutathione thermolysis to H2S, as measured; the residue is unassigned. Only this channel is "
        "modelled.",
    ),
    Reaction(
        "r_cys_h2s", {"Cys": 1},
        {"H2S": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1}, "k_cys_h2s",
        "B2.3 CHARGE FIX: the products were {H2S, FRAG_C 3, FRAG_N 1}, which "
        "deleted cysteine's alpha-carboxylate AND its ammonium into "
        "untitratable pools. Cysteine's thermal beta-elimination route is "
        "textbook -- H2S + NH3 + pyruvic acid -- so both centres have a named "
        "destination and neither is invented: the carboxyl goes to CBX (as "
        "pyruvate; CBX's declared pKa 4.25 approximates pyruvic acid's 2.49) "
        "and the nitrogen leaves as NH3 -- whose FATE, not whose release, "
        "is what CENTRE_LEDGER declares (ph_state.AMINE_FATE_BASIS). "
        "MEASURED, and the only pH-resolved rate in the module. Zheng & Ho 1994, "
        "first-order H2S release from aqueous cysteine at pH 3/5/7/9, 80-110 C, "
        "with the CORRECTED prefactors (9.8e11 / 1.9e12 / 2.4e13 / 1.0e12 1/s), "
        "not the repo's refuted 1.0e14. The C3 and N1 that leave with the "
        "sulfide are unmeasured and are routed to FRAG_C / FRAG_N.",
    ),
    Reaction(
        "r_cys_thermal", {"Cys": 1},
        {"CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "FRAG_S": 1},
        "k_cys_thermal",
        "B2.3 CHARGE FIX: the products were {FRAG_C 3, FRAG_N 1, FRAG_S 1}. "
        "This step is a LUMP over cystine formation, non-sulfide degradation "
        "and self-condensation, and the largest named member -- cystine -- "
        "RETAINS both of cysteine's titratable centres, so carrying them is "
        "the minimally-perturbing encoding as well as the conservative one. "
        "At pH 5-8 a zwitterion's two centres nearly cancel, so this fix is "
        "close to charge-neutral by construction; that it is close to neutral "
        "is precisely why the pair had to move TOGETHER. "
        "B2.1, NEW, AND IT CARRIES A MEASURED ACTIVATION ENERGY. Kang 2026 SI "
        "Fig. S4 tracks FREE cysteine, 10 mmol/L, pH 7, sealed, and measures "
        "16.2% and 38.7% conversion by 120 min at 100 and 120 C -- the "
        "corpus's first activation energy for cysteine CONSUMPTION under "
        "Maillard conditions, Ea = 55.1 kJ/mol, R^2 0.994. "
        "IT IS NOT THE SAME REACTION AS r_cys_h2s AND MUST NOT BE CONFUSED "
        "WITH IT. Zheng & Ho measure H2S RELEASE (Ea 133 kJ/mol) and Zheng's "
        "own constants account for only ~2.7% conversion at 100 C against "
        "Kang's 16.2%, so most of what consumes cysteine in a sealed aqueous "
        "system releases no sulfide at all. This step is that remainder: "
        "cystine formation, non-sulfide degradation and self-condensation, "
        "lumped. Its RATE is fitted; its BARRIER is Kang's measurement, so the "
        "temperature dependence of the cysteine pool is not a free shape. "
        "THE MEASUREMENT'S OWN QUALIFICATIONS TRAVEL WITH IT: the curves are "
        "not cleanly first order (120 and 140 C decelerate, 100 C mildly "
        "accelerates), the figure is DIGITISED rather than printed, n is not "
        "stated, and the identification of the system as free cysteine rather "
        "than the TTCA-bound moiety is 85% confident on the dossier's own "
        "reading.",
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
        "r_arp_dpo", {"ARP": 1},
        {"DPO": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1}, "k_arp_dpo",
        "B2.3 CHARGE FIX, AND THIS IS THE SITE B2.2's DIAGNOSIS DID NOT NAME. "
        "sec. 3 listed the four cysteine sinks; it missed that all FOUR "
        "Amadori enolisations delete the Amadori's own amino-acid moiety "
        "(C3 N1 -- the alanine) into FRAG_C / FRAG_N, taking a carboxylate and "
        "an ammonium with it. On Zhou's pot, whose entire charge inventory is "
        "20 mmol/L of fed ARP, that is the DOMINANT leak, not a minor one. It "
        "is also the mechanical reason the B2.2 fit drove the Amadori ammonium "
        "pKa to 10.77 (diagnosis sec. 5 contradiction 2): a pKa that high "
        "makes ARP electrically invisible, which is the cheapest way for an "
        "optimiser to stop a leak it cannot repair. The liberated amino acid "
        "acid's CARBOXYL goes to CBX, its amine centre is declared destroyed "
        "in CENTRE_LEDGER, and the sugar residue stays in FRAG_C. "
        "The fed Amadori's 2,3-enolisation. Base-catalysed. Zhou 2023's feed "
        "partitions between this and r_arp_tdp, which is why MFT (this branch) "
        "and FFT (the acid branch) move in OPPOSITE pH directions in one pot.",
    ),
    Reaction("r_arp_tdp", {"ARP": 1},
             {"TDP": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1},
             "k_arp_tdp", "The fed Amadori's 1,2-enolisation. Acid-catalysed. "
                          "B2.3 charge fix, as r_arp_dpo."),
    Reaction(
        "r_arp_tdp_thermal", {"ARP": 1},
        {"TDP": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1},
        "k_arp_tdp_th",
        "B2.3 charge fix, as r_arp_dpo. "
        "B2.1: THE UNCATALYSED 1,2-ENOLISATION OF THE FED AMADORI, and it is "
        "the fix for B2's dominant defect rather than a new chemical claim. If "
        "a product has exactly one route and that route is fully rate-limited "
        "by [H+], then the product is one decade per pH unit BY CONSTRUCTION, "
        "with no chemistry left over -- which is precisely what B2's hold-out "
        "report diagnosed. The uncatalysed route is the same reaction without "
        "the proton, it is what keeps caramelisation running in neutral and "
        "alkaline solution, and its free-sugar analogue (r_pent_caramel) was "
        "already in the network and already forced there by Hofmann's "
        "amine-free fed rows. With both present the OBSERVABLE pH slope is "
        "sub-decade and its crossover is a ratio of two ordinary rate "
        "constants fitted on FIT rows -- not a pH parameter.",
    ),
    Reaction(
        "r_arp_dpo_thermal", {"ARP": 1},
        {"DPO": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1},
        "k_arp_dpo_th",
        "The base-lane partner of r_arp_tdp_thermal, for the same reason. "
        "B2.3 charge fix, as r_arp_dpo.",
    ),
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
        {"DDP": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "FRAG_S": 1}, "k_dpo_ddp",
        "B2.3 CHARGE FIX, AND IT IS THE ONE THE SOURCE WRITES OUT FOR US: "
        "Nedvidek's own equation ends '+ RCHO + CO2 + NH3'. The CO2 is the "
        "cysteine carboxyl leaving as an acid and the NH3 is its amino "
        "nitrogen leaving. The CO2 half is transcription and goes to CBX; "
        "the NH3 half is where the module must make a claim, and it makes "
        "it in CENTRE_LEDGER on ph_state.AMINE_FATE_BASIS rather than here. "
        "(CBX carries the mixture's pKa 4.25 against "
        "carbonic acid's 6.35 -- the declared approximation on the CBX "
        "species; a carbonate pool would need a CO2 partition model that no "
        "bundle in this corpus supports.) "
        "Nedvidek 1992 Scheme 3: 1-deoxyosone + alpha-amino acid -> "
        "1,4-dideoxyosone + RCHO + CO2 + NH3, balance verified C7H13NO6 both "
        "sides for glycine. THE AMINO ACID IS A REAGENT: this step CONSUMES "
        "cysteine, putting the sulfur donor and the carbon skeleton in "
        "competition for one pool. The cysteine's own sulfur leaves in the "
        "Strecker aldehyde (mercaptoacetaldehyde) and is routed to FRAG_S "
        "rather than being silently promoted to H2S.",
    ),
    Reaction(
        "r_dpo_ddp_h2s", {"DPO": 1, "Cys": 1},
        {"DPO": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "H2S": 1}, "k_dpo_ddp_h2s",
        "B59a (revised 2026-10-29 after the independent review). The Strecker step's other fate: its aminoreductone "
        "(the dicarbonyl having taken cysteine's amine and 2 H) reduces the mercaptoacetaldehyde to H2S + "
        "acetaldehyde (Kobayasi & Fujimaki 1965, 'reductive cleavage of the C-S bond' by the enaminol) and, "
        "oxidised back to the imine, hydrolyses to THE SAME 1-deoxyosone + NH3. Net: Cys + H2O -> H2S + CH3CHO + "
        "CO2 + NH3, the deoxyosone a catalyst; H, O and electrons balance (the first version emitted the 2e-reduced "
        "DDP with the H2S and counted the reduction twice). Reductant-limited by construction: one H2S per Strecker "
        "event at most. DERIVED: STRECKER_SULFUR.h2s_share x the Strecker constant, zero when the switch is off.",
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
    # B2.2: every deoxyosone-sink event now deposits ONE acid EQUIVALENT in the
    # terminal ACID pool, on the direct analogy of Martins 2005's own measured
    # hexose steps (3-deoxyglucosone -> formic acid, step 5; 1-deoxyglucosone ->
    # acetic acid, step 8). NOTHING CONSUMES `ACID` and no rate changed: the
    # only observable it can move is the pH, through
    # ph_state.titratable_inventory, and it enters there scaled by the single
    # calibrated acid-yield fraction. Carbon balances step by step (the pool is
    # counted on a formic basis, 1 C, and the residue stays in FRAG_C).
    Reaction("r_arp_decay", {"ARP": 1},
             {"ACID": 1, "CBX": 1, "FRAG_C": 6, "FRAG_N": 1},
             "k_osone_decay",
             "B2.3 charge fix. The step keeps its ONE lumped-sink acid "
             "equivalent (ACID, scaled by the fitted yield) and additionally "
             "carries the Amadori's own amino-acid CARBOXYL into CBX at full "
             "strength -- it is not a fraction of anything -- while its amine "
             "centre is declared destroyed in CENTRE_LEDGER."),
    Reaction("r_osone_decay_dpo", {"DPO": 1}, {"ACID": 1, "FRAG_C": 4},
             "k_osone_decay", ""),
    Reaction("r_osone_decay_tdp", {"TDP": 1}, {"ACID": 1, "FRAG_C": 4},
             "k_osone_decay", ""),
    Reaction("r_osone_decay_ddp", {"DDP": 1}, {"ACID": 1, "FRAG_C": 4},
             "k_osone_decay", ""),
    Reaction("r_ptr_decay", {"PTR": 1}, {"ACID": 1, "FRAG_C": 4},
             "k_osone_decay", ""),
    # ================= THIOL FORMATION ====================================
    Reaction(
        "r_ddp_mft", {"DDP": 1, "H2S": 1}, {"MFT": 1}, "k_ddp_mft",
        "THE INTACT-SKELETON ROUTE. Cerny 2003 T2: 49% unlabelled / 46% 13C5 "
        "with no fragment pattern, 'pathways via ribose fragmentation were not "
        "relevant' => ~93% of MFT carries the intact pentose skeleton at 95 C.",
    ),
    Reaction(
        "r_ddp_mft_hs", {"DDP": 1, "H2S": 1}, {"MFT": 1}, "k_ddp_mft_hs",
        "B2.1: THE SAME ADDITION THROUGH THE HYDROSULFIDE ANION. A nucleophile "
        "with two protonation states adds through both; B2 carried only the "
        "neutral molecule, which is a mechanism claim that runs against the "
        "elementary chemistry (HS- is by far the better nucleophile and H2S "
        "dominates only by abundance below the pKa) and which is half of why "
        "B2's thiols collapsed at high pH. Two branches, two constants, both "
        "fitted on FIT rows: the crossover is a rate ratio, not a pH parameter.",
    ),
    Reaction("r_nf_mft", {"NF": 1, "H2S": 1}, {"MFT": 1}, "k_nf_mft", ""),
    Reaction(
        "r_nf_mp3p", {"NF": 1, "H2S": 1}, {"MP3P": 1}, "k_nf_mp3p",
        "Norfuraneol's DOMINANT fate: mercaptoketones : MFT = 16.3 : 1 from fed "
        "NF, MFT only 2.6% of everything fed NF produces (Whitfield 1999). "
        "2-mercapto-3-pentanone is 96% unlabelled from Cerny's NF spike -- the "
        "sharpest single-species NF-route marker in the corpus.",
    ),
    # ---- B49 (2026-10-27): norfuraneol's other sulfide products -------------------------
    # Whitfield & Mottram 1999 / 2001 Table 2: fed NF + H2S gives mostly dithiolanones,
    # thiophenones and thiophenes, and the split moves away from thiols between pH 4.5 and
    # 6.5 even when H2S is fed directly. One lumped step through the HYDROSULFIDE ANION, so its
    # pH response is the declared mechanism and not a fitted slope. INERT (zero) unless a B49
    # report supplies `k_nf_shet`, so every earlier wave reproduces bit for bit.
    Reaction(
        "r_nf_shet", {"NF": 1, "H2S": 1}, {"SHET": 1}, "k_nf_shet",
        "B49. Norfuraneol + HS- -> sulfur heterocycles (lumped). Zero until B49 supplies it.",
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
        "r_fur_fft_hs", {"FUR": 1, "H2S": 1}, {"FFT": 1}, "k_fur_fft_hs",
        "The FFT partner of r_ddp_mft_hs.",
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
    # WAVE B50 (2026-10-27): thiamine also releases H2S at pH >= 6.5 (Dwivedi & Arnold 1971). A
    # DERIVED rate, not a constant: share(pH) / (1 - share) x the two thiamine steps' total, so the
    # measured share of degraded thiamine goes to sulfide. Zero at pH < 6.25 -- every fit pot.
    Reaction("r_thi_h2s", {"THI": 1}, {"H2S": 1, "FRAG_C": 12, "FRAG_N": 4}, "k_thi_h2s",
             "B50. Thiamine -> H2S, as the measured share of thiamine's thermal loss (Dwivedi & Arnold 1971)."),
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
        "r_cys_actz", {"Cys": 1, "MGO": 1}, {"ACTZ": 1, "CBX": 1}, "k_cys_actz",
        "B2.3 CHARGE FIX, AND IT IS THE ONE PLACE THE FIX IS ASYMMETRIC. The "
        "carboxyl is CARRIED (thiazole formation from cysteine and a "
        "dicarbonyl is a decarboxylation, so the group leaves as CO2 and "
        "lands in CBX -- it replaces the FRAG_C 1 the step used to emit). The "
        "AMINE IS NOT CARRIED AND MUST NOT BE: the cysteine nitrogen ends up "
        "in the THIAZOLE RING of the product, already accounted by the "
        "nitrogen balance, and a thiazolium's pKa near 2.5 is outside this "
        "module's window. That centre is genuinely destroyed and the step "
        "genuinely releases a proton. CENTRE_LEDGER declares the -1 rather "
        "than depositing a titratable amine equivalent that would "
        "manufacture base. "
        "Zhou 2023 measures 2-acetylthiazole in BOTH an ARP-fed and an MGO-fed "
        "system at pH 8 (582 vs 665 ug/L). That pair is a step-level, "
        "two-system constraint on the ARP -> MGO flux with no free downstream "
        "parameter, and it is a declared HOLD-OUT.",
    ),
    # ================= B2.1: THE FED TTCA INTERMEDIATE =====================
    Reaction(
        "r_ttca_cys", {"TTCA": 1}, {"Cys": 1, "PENT": 1}, "k_ttca_cys",
        "B2.1: TTCA's ring-opening back to free cysteine and the pentose. Kang "
        "2026 charges purified TTCA at 10 mmol/L and measures free cysteine "
        "RISING THEN FALLING at all three temperatures, with the maximum moving "
        "earlier as T rises -- consecutive A -> B -> C kinetics. This is the "
        "release half; the module's existing cysteine sinks are the consumption "
        "half, so the transient emerges rather than being imposed.",
    ),
    Reaction(
        "r_ttca_deg", {"TTCA": 1},
        {"DPO": 1, "CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "FRAG_S": 1},
        "k_ttca_deg",
        "B2.3 CHANGES B2.2's FIX IN TWO WAYS AND BOTH ARE DELIBERATE. (i) The "
        "acid equivalent moves from ACID to CBX: TTCA's carboxyl is a carried "
        "certainty, not a fraction of the lumped osone sink, so it must not "
        "be multiplied by the fitted acid yield -- B2.2 had it scaled by "
        "0.968 for no stated reason. (ii) The ring nitrogen's centre is now "
        "DECLARED DESTROYED in CENTRE_LEDGER instead of being dropped into "
        "FRAG_N in silence. This step leans on that declaration LESS than any "
        "other, which is worth stating because Kang's whole ladder runs on "
        "it: TTCA's ring nitrogen is a secondary amine of pKa 6.24, already "
        "only ~15% protonated at Kang's initial pH 7, where an alpha-ammonium "
        "would be essentially fully protonated. "
        "B2.1: TTCA's OTHER fate, and it is what the 16.3 mol% ceiling is made "
        "of. Peak free cysteine is 1.63 mmol/L against 10 mmol/L TTCA charged, "
        "so >=84% of the cysteine moiety leaves by a route that never passes "
        "through free cysteine. A network with only the release route would "
        "have to violate that bound; this step is how it can satisfy it. "
        "B2.2 ADDS ONE ACID EQUIVALENT AND IT IS A CHARGE-CONSERVATION FIX, "
        "not a chemistry claim: TTCA is a CARBOXYLIC ACID and at Kang's "
        "initial pH 7 its carboxylate carries ~1 equivalent of the NaOH that "
        "was used to set that pH. A degradation route that deletes the "
        "carboxylate while leaving the sodium behind manufactures strong base "
        "out of nothing -- measured doing exactly that in this wave's probe, "
        "which drove Kang's pot to pH 11 against a measured 4.9. The carboxyl "
        "survives its molecule and is carried here. (The release route "
        "r_ttca_cys needs no such term: the carboxyl leaves inside cysteine, "
        "which the charge balance already tracks.)",
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
    Reaction("ch_thioether_fft", {"FFT": 1, "MELE": 1}, {"BND_F": 1},
             "k_thioether",
             "B2.1: FFT now has its OWN adduct reservoir -- see "
             "ch_thioether_release_fft."),
    Reaction(
        "ch_thioether_release", {"BND": 1}, {"MFT": 1, "MELE": 1}, None,
        "THE CONJUGATION IS REVERSIBLE, and B2 was right to say so and wrong "
        "about the number. Stack 2018 measures BOTH the forward and the reverse "
        "constant at four temperatures, so K = k1/k-1 is MEASURED, not a range: "
        "5.64 M^-1 at 19.4 C rising to 10.57 M^-1 at 4.6 C. B2 used 316 M^-1 on "
        "this paper's authority; the paper does not support it, and the "
        "correction is 56x. "
        "AND K FALLS WITH TEMPERATURE -- the addition is exothermic, dH = -28.5 "
        "kJ/mol, confirmed to the last digit by two independent Arrhenius "
        "refits -- so HEATING RELEASES BOUND THIOL. That is the measured "
        "mechanism behind Hofmann's 80 C brew being slower than his own 30 C "
        "models, and B2 could not express it because its K had no temperature "
        "dependence at all. "
        "ITS RATE IS STILL DERIVED, NOT FITTED: k_reverse = k_thioether(T) / "
        "K(T). "
        "B2.1 ALSO RETIRES B2's STATED LIMITATION that the shared adduct "
        "silently converted FFT into MFT: BND and BND_F are separate species "
        "and each releases the thiol that was bound.",
    ),
    Reaction(
        "ch_thioether_release_fft", {"BND_F": 1}, {"FFT": 1, "MELE": 1}, None,
        "The FFT half of the split. Same derived rate, same measured K(T).",
    ),
    # Channel 1b -- B2.1: thiol-disulfide exchange with matrix PROTEIN.
    Reaction(
        "ch_protein_ss_fft", {"FFT": 1, "PROT_SS": 1}, {"PRB": 1},
        "k_protein_ss",
        "B2.1, NEW CHANNEL, required by FIT_HOLDOUT_DECLARATION.md Amendment 5. "
        "anantharamkrishnan2020b Fig. 7c: 2-furfurylthiol cleaves "
        "beta-lactoglobulin's disulfide linkages to give BOTH 1:1 (+114 Da) and "
        "2:1 (+228 Da) adducts within 24 h at ambient temperature. The sulfur "
        "branch had no protein sink at all. "
        "IT IS STOICHIOMETRIC, NOT CATALYTIC -- capacity is the protein's own "
        "disulfide inventory -- which is why the partner is a titrated SITE "
        "POOL. PROT_SS IS ZERO IN EVERY SYSTEM IN BOTH SCORED PANELS, so this "
        "channel carries EXACTLY ZERO FLUX in every fit and hold-out row, by "
        "mass action and not by a flag. Its rate is BOUNDED from the paper's "
        "single timescale bracket, never fitted: no FIT row contains protein.",
    ),
    Reaction(
        "ch_protein_ss_mft", {"MFT": 1, "PROT_SS": 1}, {"PRB": 1},
        "k_protein_ss",
        "The MFT arm, and it is an INFERENCE rather than a measurement: "
        "anantharamkrishnan measures FFT and n-propanethiol, and the dossier "
        "flags the extension to the furan-3-thiol class as its own inference. "
        "Carried at the same bounded rate with that flag attached.",
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
        "ch_dimer_mft", {"MFT": 2, "OX": 1}, {"MFTD": 1, "OXV": 1}, "k_dimer_mft",
        "CHANNEL 3, dominant at 115-120 C -- the channel the 30 C system rules "
        "out at <1.5%. SECOND order in thiol and FIRST in oxidant equivalents, "
        "because Zhang 2024 Fig. 1 shows the branch responds to the additive's "
        "REDOX STATE and not to its concentration: at near-matched molar sulfur "
        "(123.8 vs 124.8 mM S) cystine sends 54.2% of the MFT pool to the dimer "
        "and cysteine only 8.6%. NOT AROMA LOSS -- the dimer is 15.6x more "
        "potent than the monomer.",
    ),
    Reaction("ch_dimer_fft", {"FFT": 2, "OX": 1}, {"FFTD": 1, "OXV": 1}, "k_dimer_fft", ""),
    # Channel 4 -- radical coupling to methanethiol. 115 C.
    Reaction(
        "ch_mmft", {"MFT": 1, "MESH": 1}, {"MMFT": 1}, "k_mmft",
        "CHANNEL 4: MFT + MeSH radical coupling. Requires a methanethiol "
        "partner, which the 30 C and 50 C systems do not have -- which is why "
        "it is a separate channel and not a term in a lumped sink.",
    ),
    # Channel 5 -- B2.1: the pH-dependent, thiolate-mediated oxidative sink.
    Reaction(
        "ch_thiolate_loss_fft", {"FFT": 1}, {"FRAG_C": 5, "FRAG_S": 1},
        "k_thiolate_loss",
        "B2.1, NEW, AND IT IS WHERE THE FFT-VERSUS-pH SLOPE ACTUALLY LIVES. "
        "Kumazawa 2003 heats 1 ppm 2-furfurylthiol at 121 C / 10 min across "
        "seven pH levels and measures survival falling 99.5 -> 96.2 -> 89.1 -> "
        "79.5 -> 45.1 -> 11 -> <0.5 % from pH 3.0 to 7.0, with the difurfuryl "
        "disulfide growing monotonically over the same grid AND a mass-balance "
        "shortfall the authors attribute to a Fenton-type route. B2 had NO pH "
        "dependence on any consumption channel, so it was forced to load the "
        "whole observed pH response of FFT onto FORMATION at a decade per pH "
        "unit. A measurable part of that response is CONSUMPTION. "
        "Thiol oxidation proceeds through the thiolate, so the factor is the "
        "thiolate fraction at reaction temperature: a mechanism claim with no "
        "fitted pH number. This step is the NON-disulfide half (the shortfall); "
        "the disulfide half is ch_dimer_*, which now carries the same factor.",
    ),
    Reaction("ch_thiolate_loss_mft", {"MFT": 1}, {"FRAG_C": 5, "FRAG_S": 1},
             "k_thiolate_loss", ""),
    # Terminal decay of the aroma species, bounded to allow ~zero.
    Reaction("ch_mft_decay", {"MFT": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_mft_decay", ""),
    Reaction("ch_fft_decay", {"FFT": 1}, {"FRAG_C": 5, "FRAG_S": 1}, "k_fft_decay", ""),
    Reaction("ch_mftd_decay", {"MFTD": 1}, {"FRAG_C": 10, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_fftd_decay", {"FFTD": 1}, {"FRAG_C": 10, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_mmft_decay", {"MMFT": 1}, {"FRAG_C": 6, "FRAG_S": 2}, "k_dimer_decay", ""),
    Reaction("ch_mesh_decay", {"MESH": 1}, {"FRAG_C": 1, "FRAG_S": 1}, "k_thiol_decay", ""),
    Reaction("ch_actz_decay", {"ACTZ": 1}, {"FRAG_C": 5, "FRAG_N": 1, "FRAG_S": 1},
             "k_thiol_decay", ""),
    # =======================================================================
    # BUILD WAVE B7 -- the three furanic edges that need a SULFUR-LANE partner
    # =======================================================================
    # The rest of the furanic channel lives on the trunk (network.py). These
    # three are here and only here because their partners -- cysteine, sulfide
    # and the 1-deoxyPENTOSONE -- exist only on this lane.
    Reaction(
        "r_dpo_af", {"DPO": 1, "Gly": 1}, {"AF": 1, "FRAG_C": 1, "FRAG_N": 1},
        "k_dpo_af",
        "B7. DMHF Edge A, PENTOSE arm, and THE ONLY FITTED CONSTANT IN THE "
        "FURANIC CHANNEL. C5 + C1 -> C6: the pentose cannot make a C6 furanone "
        "from its own skeleton, so the sixth carbon is the STRECKER ALDEHYDE "
        "of the amino acid. Blank & Fay 1996 measured that split by 13C "
        "labelling (~70 % Strecker / ~30 % fragmentation) and Blank 1997 "
        "reproduced it NON-ISOTOPICALLY at 73/27 using 4-aminobutyric acid, "
        "which does not undergo Strecker deamination, as a fragmentation-only "
        "baseline. Two methods, three percentage points apart, same system. "
        "DECLARED DEFECT: this step is FIRST ORDER in the amine and Blank's "
        "Table 5 measures the apparent order at n = 0.35-0.71. That is "
        "hold-out H8 and it is declared in advance, not discovered.",
    ),
    Reaction(
        "r_hmf_cys", {"HMF": 1, "Cys": 1},
        {"HMFAD": 1, "CBX": 1, "FRAG_N": 1}, "k_hmf_cys",
        "B7. THE ONLY MEASURED HMF SINK IN THE CORPUS THAT SURVIVES AUDIT "
        "(Hamzalioglu & Gokmen 2018; K5a sec. 4.1 refuses the other seven, one "
        "of them because its own authors tested and rejected it). Written "
        "BIMOLECULAR AND COMPETITIVE rather than as a fixed first-order decay, "
        "which is K5a MUST-NOT #4: at 160-200 C the edge is out-competed for "
        "the amine pool by glyoxal, 1-DG and MGO, and writing it first order "
        "would make it uncompetable. "
        "Its constant is CLAMPED at 323.15 K -- the data span 5-50 C and "
        "extrapolating them to cooking temperature is a named prohibited "
        "derivation. The carboxyl is CARRIED into CBX and the amine centre is "
        "declared destroyed in CENTRE_LEDGER, the same treatment every other "
        "cysteine-consuming step in this network gets.",
    ),
    Reaction(
        "r_dmhf_h2s", {"DMHF": 1, "H2S": 1}, {"DMHFS": 1}, "k_dmhf_h2s",
        "B7. DMHF Edge C, THE SINK, AND ITS RATE IS EXACTLY ZERO. The "
        "stoichiometry balances exactly (C6H8O3 + H2S -> C6H8O2S + H2O) and it "
        "is the exact C6 counterpart of this module's existing C5 edge "
        "furanone_reductive_sulfhydrylation. What does not exist is a "
        "magnitude: Shu & Ho 1988 report a GC area percent and no residual "
        "DMHF, no conversion, no mass balance and no molar yield of anything. "
        "k5b_dmhf_synthesis.md sec. 8.6 names the failure mode a constant "
        "fitted to its 6.0 % would repeat. The edge ships structural.",
    ),
    # =======================================================================
    # WAVE B11 (2026-09-07) -- OXYGEN AS AN INPUT (kinetic_core_b11_prereg.md)
    # =======================================================================
    # Dissolved oxygen (OX, ambient units) is resupplied from the headspace
    # reservoir (OXR) into the vacancy the consumers create (OXV), and consumed by
    # cysteine autoxidation and by the reductone pool. With k_cys_ox = k_red_ox = 0
    # and OXR = 0 -- the defaults every wave before B11 runs at -- the four steps
    # carry no flux and every earlier artefact reproduces exactly.
    Reaction(
        "ox_supply", {"OXR": 1, "OXV": 1}, {"OX": 1}, "k_ox_supply",
        "B11. Gas-liquid resupply: fast (a declared constant), proportional to the "
        "vacancy, so OX sits at saturation while the reservoir lasts and falls with "
        "consumption once it is spent.",
    ),
    Reaction(
        "ch_cys_ox", {"Cys": 1, "OX": 1},
        {"CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "FRAG_S": 1, "OXV": 1}, "k_cys_ox",
        "B11. Cysteine autoxidation (to cystine and beyond), first order in cysteine "
        "and in dissolved oxygen, through the thiolate. The products mirror "
        "r_cys_thermal's routing (carboxyl carried as CBX). NO measurement in the "
        "corpus pins its rate (Bagiyan 2004 prints initial rates only); B11 fits it "
        "and expects the objective not to identify it.",
    ),
    Reaction(
        "ch_red_ox_dpo", {"DPO": 1, "OX": 1}, {"FRAG_C": 5, "OXV": 1}, "k_red_ox",
        "B11. The reductone pool consumes oxygen: the pentose deoxyosone arm. Shares "
        "k_red_ox with the norfuraneol arm; unmeasured, fitted, expected unidentified.",
    ),
    Reaction(
        "ch_red_ox_nf", {"NF": 1, "OX": 1}, {"FRAG_C": 5, "OXV": 1}, "k_red_ox",
        "B11. The norfuraneol arm of the reductone oxygen sink.",
    ),
    # ---- B17 (2026-09-08): the disulfide gives the thiol back --------------------------
    # One shared constant, `k_dimer_release`, INERT (zero) unless a B17 report supplies it,
    # so every wave before B17 reproduces bit for bit (the B11 discipline). Its barrier is
    # the dimerisation's own measured one (Zhang 2026 k17, MEASURED_EA_OVERRIDES), so the
    # equilibrium constant, not the two rates, carries the temperature dependence
    # (kinetic_core_b17_prereg.md sec. 2, variant b).
    Reaction(
        "ch_dimer_release_mft", {"MFTD": 1}, {"MFT": 2}, "k_dimer_release",
        "B17. bis(2-methyl-3-furyl) disulfide -> 2 MFT (reduction back to the thiol; the "
        "reducing partner is the pot's reductone pool, not tracked). Zero until B17 supplies it.",
    ),
    Reaction(
        "ch_dimer_release_fft", {"FFTD": 1}, {"FFT": 2}, "k_dimer_release",
        "B17. bis(2-furfuryl) disulfide -> 2 FFT; shares k_dimer_release.",
    ),
    # ---- B17 variant (a) (2026-09-09): the pot makes its own electrophile sites --------------
    # A parallel branch of the deoxyosone decay in which the browning carbon carries ONE matrix
    # electrophile site (MELE, zero atoms) into the pool the measured thioether channel drains.
    # Its rate is the site yield times k_osone_decay (same barrier, the carbonyl-sink family);
    # at the inert default yield of zero every wave before B17a reproduces bit for bit. The
    # carbon goes to FRAG_C in full (5 per osone) and no ACID is made here, so the pH drift is
    # untouched. Pre-registration: kinetic_core_b17_prereg.md sec. 2, variant (a).
    Reaction(
        "ch_mele_from_dpo", {"DPO": 1}, {"FRAG_C": 5, "MELE": 1}, "k_mele_site",
        "B17a. 1-deoxypentosone -> browning fragments + one electrophile site.",
    ),
    Reaction(
        "ch_mele_from_tdp", {"TDP": 1}, {"FRAG_C": 5, "MELE": 1}, "k_mele_site",
        "B17a. 3-deoxypentosone -> browning fragments + one electrophile site.",
    ),
    Reaction(
        "ch_mele_from_ddp", {"DDP": 1}, {"FRAG_C": 5, "MELE": 1}, "k_mele_site",
        "B17a. 1,4-dideoxypentosone -> browning fragments + one electrophile site.",
    ),
    # ---- B25 (2026-09-09): the thiol sink, variant (c) -----------------------------------------
    # The meaty thiols add IRREVERSIBLY to the pot's own unsaturated carbonyls, the two deoxypentosones
    # the lane carries while the sugar lasts; the adduct joins the terminal oligomer pool (OLG, C5 S1;
    # the osone's carbons to FRAG_C). One shared second-order constant with its own fitted barrier;
    # at the inert default of zero every wave before B25 reproduces bit for bit. Pre-registration:
    # kinetic_core_b25_prereg.md.
    Reaction("ch_add_mft_dpo", {"MFT": 1, "DPO": 1}, {"OLG": 1, "FRAG_C": 5}, "k_add",
             "B25. MFT + 1-deoxypentosone -> irreversible adduct (thiol-Michael on the enol; net)."),
    Reaction("ch_add_mft_tdp", {"MFT": 1, "TDP": 1}, {"OLG": 1, "FRAG_C": 5}, "k_add",
             "B25. MFT + 3-deoxypentosone -> irreversible adduct."),
    Reaction("ch_add_fft_dpo", {"FFT": 1, "DPO": 1}, {"OLG": 1, "FRAG_C": 5}, "k_add",
             "B25. FFT + 1-deoxypentosone -> irreversible adduct."),
    Reaction("ch_add_fft_tdp", {"FFT": 1, "TDP": 1}, {"OLG": 1, "FRAG_C": 5}, "k_add",
             "B25. FFT + 3-deoxypentosone -> irreversible adduct."),
    # ---- B27 (2026-09-11): the oxidant the three refused sink waves were never given ------------
    # Whitfield & Mottram 1999 p. 1631 and Figure 6: the alpha-dicarbonyl reduced to a hydroxyalkanone
    # on the way to a mercaptoketone is the redox system that makes the disulfides, because aerial
    # oxidation is unfavourable at these H2S levels. The lane already carries that flux as r_nf_mp3p,
    # so the oxidant is delivered as a CO-PRODUCT of a parallel branch of it: this step runs at
    # phi x k_nf_mp3p and r_nf_mp3p at (1 - phi) x k_nf_mp3p (parameters_sulfur.apply_dicarbonyl_redox),
    # same barrier and pH factor, so the total mercaptoketone flux is unchanged for every phi and one
    # OX equivalent (one disulfide) is made per mercaptoketone formed through this branch. phi is the
    # ONE fitted coordinate; at the inert default of zero every wave before B27 reproduces bit for bit.
    # OX here is chemistry-derived and NOT bounded by the O2 saturation that ox_supply respects; it
    # shares the pool because the consumers (ch_dimer_*, ch_cys_ox, ch_red_ox_*) cannot tell them apart.
    # Pre-registration: kinetic_core_b27_prereg.md sec. 3 and sec. 10.
    Reaction("ch_redox_mp3p", {"NF": 1, "H2S": 1}, {"MP3P": 1, "OX": 1}, "k_redox_mp3p",
             "B27. NF + H2S -> 2-mercapto-3-pentanone + one oxidant equivalent (the dicarbonyl's "
             "reduction on the way to the mercaptoketone oxidises a thiol; Whitfield 1999 Figure 6)."),
    # B58 (2026-10-27): the condensation r_ttca_cys reverses. Rate DERIVED, k = K x k_ttca_cys, with K
    # measured by Cao 2017 (parameters_sulfur.TTCA_FORMATION); zero while the wave's switch is off.
    Reaction("r_ttca_form", {"Cys": 1, "PENT": 1}, {"TTCA": 1}, "k_ttca_form",
             "B58. Cysteine + pentose -> TTCA, the forward half of the measured TTCA equilibrium "
             "(Cao 2017, K ~ 60 L/mol); detailed balance against r_ttca_cys."),
    # B61 (2026-10-28): norfuraneol's hydroxide-catalysed decay, measured at pH 8 / 100 C by Mikami 2017
    # (parameters_sulfur.NF_HYDROXIDE_DECAY); in parallel with the pH-independent r_nf_decay, zero while
    # the wave's switch is off.
    Reaction("r_nf_decay_oh", {"NF": 1}, {"FRAG_C": 5}, "k_nf_decay_oh",
             "B61. Norfuraneol -> fragments through its enolate, specific-base catalysed (Mikami 2017)."),
    # B59b (2026-10-29): DE ROOS 2005's ROUTE. TTCA is a reservoir in fast equilibrium with the open Schiff base;
    # the Schiff base rearranges slowly to the cysteine Amadori compound, whose enolisation gives the
    # deoxyosones and RELEASES the cysteine. All three are DERIVED rates (CYS_AMADORI_ROUTE), zero when off.
    Reaction("r_cama_form", {"TTCA": 1}, {"CAMA": 1}, "k_cama_form",
             "B59b. TTCA -> (Schiff base) -> cysteine Amadori compound. First order in the TTCA pool, which holds "
             "~89 % of the cysteine within 0.5 h at 90 C (Cao 2017). Constant anchored on Cao's own blank pot."),
    Reaction("r_cama_dpo", {"CAMA": 1}, {"DPO": 1, "Cys": 1}, "k_cama_dpo",
             "B59b. The Amadori compound's 2,3-enolisation to the 1-deoxyosone, the amino acid released "
             "(de Roos 2005). Share of the measured decay: Zhai 2019's fed-TTCA 3-DX / 1-DX = 2.1-2.6."),
    Reaction("r_cama_tdp", {"CAMA": 1}, {"TDP": 1, "Cys": 1}, "k_cama_tdp",
             "B59b. The 1,2-enolisation to the 3-deoxyosone, the amino acid released (de Roos 2005)."),
)

def inactive_reactions() -> frozenset:
    """Reactions in the network whose wave switch is OFF: their rate is structurally zero, so nothing reading
    the network (the hypothesis layer) may count them as modelled. B58's TTCA formation and B61's
    hydroxide-catalysed norfuraneol decay sit in SULFUR_REACTIONS either way."""
    from . import parameters_sulfur as _ps

    off = set()
    if not (_ps.TTCA_FORMATION.enabled or _ps.CYS_AMADORI_ROUTE.enabled):
        off.add("r_ttca_form")
    if not _ps.CYS_AMADORI_ROUTE.enabled:
        off.update({"r_cama_form", "r_cama_dpo", "r_cama_tdp"})
    else:
        off.update({"r_pent_dpo", "r_pent_tdp"})
    if not _ps.NF_HYDROXIDE_DECAY.enabled:
        off.add("r_nf_decay_oh")
    return frozenset(off)


#: The full network: B1's trunk first, then the sulfur block.
#: B39-B41 (2026-09-11): the fed 3-deoxy triangle (the reverse hydration and the epimer) is trunk chemistry
#: on the furanic channel this lane already carries, and its two existing constants (k_tdg_ddg, k_ddg_hmf)
#: are shared through the operative table -- so the lane composes the three steps too, or its 3-DG limb
#: would run the fitted constants without the return paths they were fitted with.
from .network import FED_3DEOXY_REACTIONS as _FED_3DEOXY_REACTIONS

FULL_REACTIONS: Tuple[Reaction, ...] = TRUNK_REACTIONS + _FED_3DEOXY_REACTIONS + SULFUR_REACTIONS
FULL_REACTION_KEYS: Tuple[str, ...] = tuple(r.key for r in FULL_REACTIONS)


# ---------------------------------------------------------------------------
# B2.3: THE CENTRE LEDGER -- every step that moves a titratable centre, and why
# ---------------------------------------------------------------------------
# `ph_state.validate_charge_closure` refuses the import unless this table
# matches the stoichiometry above EXACTLY: every step whose net carboxyl or
# amine count is non-zero appears here with that delta and a written basis, and
# no entry here describes a step that moves nothing. After the B2.3 fix every
# entry falls into one of three classes, and the fact that the FOURTH class is
# EMPTY is the whole result:
#
#   * CREATION of a neutral acid from neutral precursors -- the trunk's
#     measured Martins acids, and the lumped osone sink. Charge-balanced as
#     written; the dissociation belongs to the equilibrium, not to the
#     stoichiometry.
#   * DECLARED DESTRUCTION of a liberated amino nitrogen, on the shared
#     ph_state.AMINE_FATE_BASIS -- a chemistry claim with a measured basis and
#     a named error direction, not an exemption.
#   * DECLARED DESTRUCTION of one nitrogen into a thiazole ring (r_cys_actz),
#     whose basis is the product's own pKa.
#   * ZERO undeclared deletions of a carboxylate or an amine. That class is
#     what B2.2's diagnosis sec. 3 named four members of and this wave's audit
#     found fourteen of; it is now empty, and it cannot be refilled without
#     this table refusing the import.

CENTRE_LEDGER: Mapping[str, Mapping[str, Any]] = {
    # ---- B1 trunk: Martins' MEASURED acid-forming steps -------------------
    # Reused verbatim from network.TRUNK_CENTRE_LEDGER so a trunk step is
    # declared in exactly one place.
    **TRUNK_CENTRE_LEDGER,
    # ---- the LUMPED deoxyosone sink ---------------------------------------
    # These five are why `ACID` exists at all and why its yield is the pH
    # model's one calibrated constant: Martins measures that PART of the
    # deoxyosone flux terminates as formic/acetic acid and part as browning
    # polymer, and nothing in the corpus measures the pentose analogue's split.
    # FOUR of the five are declared here; the fifth, `r_arp_decay`, is declared
    # once below with the B2.3 amine movement as well (2026-09-04: a pre-B2.3
    # copy of it stood here, silently shadowed by the live entry and asserting
    # the opposite -- that the Amadori amine is CARRIED, not destroyed. Python
    # collapses a repeated key before `validate_charge_closure` can see it, so
    # the validator could not catch it; `test_no_dict_literal_repeats_a_key`
    # now does).
    "r_osone_decay_dpo": {"carboxyl": +1, "basis": (
        "the lumped deoxyosone sink's acid equivalent; see r_arp_decay.")},
    "r_osone_decay_tdp": {"carboxyl": +1, "basis": (
        "the lumped deoxyosone sink's acid equivalent; see r_arp_decay.")},
    "r_osone_decay_ddp": {"carboxyl": +1, "basis": (
        "the lumped deoxyosone sink's acid equivalent; see r_arp_decay.")},
    "r_ptr_decay": {"carboxyl": +1, "basis": (
        "the lumped deoxyosone sink's acid equivalent; see r_arp_decay.")},
    # ---- THE DECLARED AMINE DESTRUCTIONS ----------------------------------
    # Eight steps liberate an amino nitrogen, and every one of them declares
    # the centre destroyed on the SAME shared basis (ph_state.AMINE_FATE_BASIS)
    # rather than on eight separate hand-waves. Their CARBOXYL halves are all
    # carried into CBX, which is the half B2.2 was measured deleting.
    #
    # ASYMMETRY IS THE POINT. B2.2's diagnosis sec. 3 said "carry the
    # carboxylate"; taken alone that would MANUFACTURE ACID, because an
    # alpha-amino acid at pH 5-8 is a zwitterion whose two centres nearly
    # cancel. B2.3 does not take it alone: the amine's fate is stated, its
    # evidence is named, and the direction of the error if the statement is
    # wrong is named too.
    "r_cys_h2s": {"amine": -1, "basis": (
        "cysteine thermolysis, Zheng & Ho's route: H2S + NH3 + pyruvic acid. "
        "The carboxyl is CARRIED (CBX, as pyruvate). " + AMINE_FATE_BASIS)},
    "ch_cys_ox": {"amine": -1, "basis": (
        "B11 cysteine autoxidation, routed exactly as r_cys_thermal: the carboxyl is "
        "CARRIED (CBX), the amine leaves with the fragment nitrogen. " + AMINE_FATE_BASIS)},
    "r_cys_thermal": {"amine": -1, "basis": (
        "Kang's lumped non-sulfide cysteine consumption (cystine formation, "
        "self-condensation, degradation). Carboxyl CARRIED. " + AMINE_FATE_BASIS)},
    "r_dpo_ddp": {"amine": -1, "basis": (
        "Nedvidek 1992 Scheme 3 writes the leaving groups out: '+ RCHO + CO2 "
        "+ NH3'. The CO2 is the carboxyl and it is CARRIED (CBX). The NH3 is "
        "released and is exactly the nitrogen the shared basis is about -- "
        "released, then consumed by the Maillard chemistry this step is part "
        "of. " + AMINE_FATE_BASIS)},
    "r_dpo_ddp_h2s": {"amine": -1, "basis": (
        "as r_dpo_ddp (B59a, revised): the Strecker event whose aminoreductone reduces the mercaptoacetaldehyde to "
        "H2S and returns as the deoxyosone; the amine leaves as NH3 either way. Carboxyl CARRIED. "
        + AMINE_FATE_BASIS)},
    "r_arp_dpo": {"amine": -1, "basis": (
        "the Amadori's amino-acid moiety leaves in the 2,3-enolisation. "
        "Carboxyl CARRIED. " + AMINE_FATE_BASIS)},
    "r_arp_tdp": {"amine": -1, "basis": (
        "as r_arp_dpo, 1,2-enolisation. " + AMINE_FATE_BASIS)},
    "r_arp_dpo_thermal": {"amine": -1, "basis": (
        "as r_arp_dpo, uncatalysed. " + AMINE_FATE_BASIS)},
    "r_arp_tdp_thermal": {"amine": -1, "basis": (
        "as r_arp_tdp, uncatalysed. " + AMINE_FATE_BASIS)},
    "r_ttca_deg": {"amine": -1, "basis": (
        "TTCA's thiazolidine ring opens and its nitrogen leaves with the C2S1 "
        "residue. Carboxyl CARRIED (CBX, at FULL strength -- B2.2 put it in "
        "ACID, where the fitted acid yield scaled a stoichiometric "
        "certainty). NOTE, because this step is the one Kang's whole ladder "
        "runs on: the ring nitrogen is a SECONDARY AMINE of pKa 6.24, i.e. "
        "already only ~15% protonated at Kang's initial pH 7, so this "
        "particular declaration moves less charge than the alpha-ammonium "
        "ones do. " + AMINE_FATE_BASIS)},
    "r_arp_decay": {"carboxyl": +1, "amine": -1, "basis": (
        "TWO movements in one step and they are unrelated. The +1 carboxyl is "
        "the LUMPED SINK's own new acid (ACID, scaled by the fitted yield), on "
        "the analogy of Martins' measured steps 5 and 8. The -1 amine is the "
        "Amadori's amino-acid nitrogen leaving; its carboxyl is separately "
        "CARRIED into CBX at full strength. " + AMINE_FATE_BASIS)},
    "r_cys_actz": {"amine": -1, "basis": (
        "THE ONE DESTRUCTION THAT NEEDS NO SHARED BASIS, because the product "
        "itself is the evidence: the cysteine nitrogen is in the THIAZOLE RING "
        "of 2-acetylthiazole (the atom balance carries it there) and a "
        "thiazolium's conjugate acid sits near pKa 2.5, below anything these "
        "pots reach. The centre is destroyed at the moment of ring closure and "
        "the step releases a proton equivalent. Its carboxyl IS carried (CBX): "
        "thiazole formation from cysteine and a dicarbonyl is a "
        "decarboxylation.")},
    # ---- B7 ---------------------------------------------------------------
    "r_hmf_cys": {"amine": -1, "basis": (
        "THE HMF SINK. Cysteine's alpha-carboxyl is CARRIED into CBX at full "
        "strength; its alpha-amine is declared destroyed on the shared basis. "
        "A CAVEAT THIS DECLARATION CANNOT RESOLVE, stated instead of hidden: "
        "Hamzalioglu confirms FOUR adducts at delta < 1 ppm -- thiol-Michael, "
        "amine-Michael, Schiff base and a 2 HMF : 1 Cys adduct -- and only two "
        "of the four actually consume the amine. The thiol-Michael adduct "
        "keeps both centres, exactly as ACRCYS does on the acrylamide lane. "
        "With no measured split among the four, the conservative direction is "
        "the one that cannot manufacture base, which is this one; if the thiol "
        "route dominates, this declaration destroys an amine that survived and "
        "the pot is predicted slightly more ACIDIC than it should be. "
        + AMINE_FATE_BASIS)},
}

#: Which pH factor multiplies which reaction. Assigning a step to a factor is a
#: MECHANISM claim (acid- or base-catalysed enolisation, sulfide nucleophile)
#: and carries no fitted number; see parameters_sulfur.ph_factor.
REACTION_PH_FACTOR: Mapping[str, str] = {
    "r_arp_dpo": "base",
    "r_arp_tdp": "acid",
    "r_pent_dpo": "base",
    "r_dpo_ddp": "base",
    "r_pent_tdp": "acid",
    # B2.1: r_pent_caramel LOSES its acid factor. It is the module's declared
    # AMINE-INDEPENDENT thermal enolisation -- the step Hofmann's amine-free fed
    # rows forced into existence -- and caramelisation is not specific-acid
    # catalysed: it runs at neutral and alkaline pH, faster if anything.
    # Labelling it "acid" made BOTH routes to the 3-deoxyosone proportional to
    # [H+], which left the free-sugar lane with no uncatalysed route at all and
    # is a second instance of the defect B2's report named once.
    "r_pent_caramel": "",
    "r_pent_thermal_dpo": "",
    "r_arp_tdp_thermal": "",
    "r_arp_dpo_thermal": "",
    "r_glc_fur": "acid",
    "r_ddp_mft": "neutral_h2s",
    "r_ddp_mft_hs": "hs_anion",
    "r_nf_mft": "neutral_h2s",
    "r_nf_mp3p": "neutral_h2s",
    # B49: the heterocycle branch adds through the hydrosulfide anion; its rise with pH is what
    # takes norfuraneol away from the thiols (Whitfield 2001, pH 6.5).
    "r_nf_shet": "hs_anion",
    # B27: the redox branch of r_nf_mp3p carries the SAME factor, so the split leaves the total flux exact.
    "ch_redox_mp3p": "neutral_h2s",
    "r_mgo_mp": "neutral_h2s",
    "r_fur_fft": "neutral_h2s",
    "r_fur_fft_hs": "hs_anion",
    # B2.1: thiol oxidation to the disulfide proceeds through the THIOLATE.
    # This is a mechanism claim and it costs no parameter. Kumazawa's grid
    # measures the disulfide growing monotonically with pH while the thiol
    # collapses, in the same tubes, which is exactly the signature.
    "ch_dimer_mft": "thiolate",
    "ch_dimer_fft": "thiolate",
    "ch_thiolate_loss_mft": "thiolate",
    "ch_thiolate_loss_fft": "thiolate",
    # B11: cysteine autoxidation also proceeds through the thiolate.
    "ch_cys_ox": "thiolate",
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

#: B2.3: THE FOURTH IMPORT-TIME INVARIANT, alongside carbon, nitrogen and
#: sulfur. It walks the WHOLE network -- B1's trunk and the sulfur block -- and
#: refuses the import unless every step's titratable-centre movement is
#: declared in CENTRE_LEDGER. Making this an import-time invariant rather than
#: a test is deliberate: the B2.2 defect survived a full unit suite because
#: nothing in the suite knew to look for it, and a module that will not load is
#: a stronger guarantee than a test somebody has to remember to run.
validate_charge_closure(FULL_REACTIONS, CENTRE_LEDGER)


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


#: Sub-model 1 (2026-10-29), pre-registered with the switch OFF: the trunk's MEASURED pH terms (Martins & van Boekel
#: 2003 Part II Table 3, ``trunk_conditions``: the three Amadori-decay steps and the 3-DG -> formic acid exit) applied
#: inside the sulfur lane at the pH the lane computes, on the label scale the table was measured on. The trunk lane
#: applies them at the declared pH; the sulfur lane, the one lane whose pH moves, never did -- so its acid supply was
#: the same at pH 5 and pH 8 (Davidek 2002's pH-stat measures 3 vs 79 mmol/L). The no-refit probe that found it:
#: results/validation/kinetic_core_ph_acid_submodel_probe*.md. The evidence audit decides whether it ships.
TRUNK_PH_TERMS: bool = False


def _trunk_ph_exponents() -> Dict[str, float]:
    from . import trunk_conditions as _tc

    out = {key: float(_tc.PH_EXPONENT_DECADES_PER_UNIT) for key in _tc.PH_STEPS}
    if _tc.THREE_DEOXY_EXIT_PH_TERM:
        out.update({key: float(exponent) for key, (exponent, _band, _source) in _tc.THREE_DEOXY_EXIT_PH.items()})
    return out


def sulfur_rate_constants_at(
    parameters: Mapping[str, Any],
    temperature_k: float,
    ph: float,
    *,
    ph_label: Optional[float] = None,
    ionic_strength_m: float = 0.0,
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
    # B51: the MEASURED pH-indexed tables were measured on 25 C labels; under B51 they read the
    # cooled pH of the state (``ph_label``), every speciation factor reads the in-situ pH.
    ph_table = float(ph_label) if (ph_label is not None and _ph.B51.label_scale_tables) else ph
    trunk_ph = _trunk_ph_exponents() if TRUNK_PH_TERMS else {}
    for reaction in FULL_REACTIONS:
        key = reaction.parameter_key
        if key is None:
            continue
        if key == "k_oligomer":
            out[reaction.key] = oligomerisation_rate()
            continue
        if key == "k_cys_h2s":
            out[reaction.key] = cysteine_thermolysis_k(ph_table, float(temperature_k))
            continue
        if key == "k_thi_h2s":
            continue  # B50: derived below from the two thiamine steps
        if key == "k_ttca_form":
            continue  # B58: derived below from r_ttca_cys by detailed balance
        if key == "k_nf_decay_oh":
            continue  # B61: a measured constant, below
        if key == "k_dpo_ddp_h2s":
            continue  # B59a: a measured share of r_dpo_ddp, below
        if key in ("k_cama_form", "k_cama_dpo", "k_cama_tdp"):
            continue  # B59b: measured / source-anchored, below
        if key == "k_gsh_h2s":
            out[reaction.key] = glutathione_thermolysis_k(ph_table, float(temperature_k))
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
        # B7: A HARD TEMPERATURE CEILING, for the one constant that has one.
        # Hamzalioglu's HMF + cysteine Ea is fitted to three points spanning
        # 5-50 C with R^2 = 0.87 and a slope resting on the single point at
        # which the authors declare the pseudo-first-order assumption
        # compromised; cooking is 120-200 C. K5a sec. 7.3 registers
        # extrapolating it as a PROHIBITED DERIVATION. The constant is
        # therefore HELD at its 50 C value rather than extrapolated -- which
        # UNDER-states the sink, and that direction is stated in the B7 report
        # because it makes the module's HMF over-prediction worse rather than
        # flattering it. A mapping rather than an ``if`` so the next capped
        # constant cannot be added by editing a branch.
        _cap = _FURANIC_TEMPERATURE_CAP_K.get(key)
        value = parameter.k_at(
            float(temperature_k) if _cap is None
            else min(float(temperature_k), float(_cap))
        )
        value *= ph_factor(
            REACTION_PH_FACTOR.get(reaction.key, ""), ph, float(temperature_k), ionic_strength_m
        )
        if key in trunk_ph:
            from .trunk_conditions import REFERENCE_PH as _TRUNK_REFERENCE_PH

            value *= 10.0 ** (trunk_ph[key] * (ph_table - _TRUNK_REFERENCE_PH))
        out[reaction.key] = value

    # DERIVED RATE 1: the reverse of the thioether conjugation. Stack 2018
    # measures BOTH constants, so K = k1/k-1 is measured and so is its
    # temperature dependence (dH = -28.5 kJ/mol, exothermic). The reverse
    # constant is k_forward(T) / K(T) and is not a free parameter. Units check:
    # k_thioether is L/(mmol*min) and K is L/mmol, so the quotient is 1/min,
    # which is the first-order unit the release step needs.
    #
    # B2.1 CHANGED TWO THINGS HERE. B2 divided by a fixed 0.3162 L/mmol taken
    # from a "K ~ 1e2-1e3 M^-1" range that Stack does not print; the measured
    # value is 5.64e-3 L/mmol at 19.4 C, 56x smaller. And B2's K was
    # temperature-independent, so it could not express the one food-relevant
    # thing Stack measures: heating pushes the equilibrium back toward FREE
    # thiol.
    # DERIVED RATE 0 (B50): thiamine's sulfide branch, the measured share of its total thermal loss.
    _share = thiamine_h2s_share(ph_table)
    out["r_thi_h2s"] = (_share / (1.0 - _share)) * (out["r_thi_hmp"] + out["r_thi_mesh"]) if _share > 0.0 else 0.0

    # DERIVED RATE 0b (B58): TTCA formation by detailed balance with its measured equilibrium constant.
    from . import parameters_sulfur as _ps

    _ttca = _ps.TTCA_FORMATION
    out["r_ttca_form"] = float(_ttca.k_l_per_mmol) * out["r_ttca_cys"] if _ttca.enabled else 0.0

    # DERIVED RATE 0e (B59b): de Roos 2005's route through TTCA and the cysteine Amadori compound, and the
    # direct cysteine-catalysed pentose enolisations it replaces.
    _route = _ps.CYS_AMADORI_ROUTE
    if _route.enabled:
        out["r_ttca_form"] = float(_ttca.k_l_per_mmol) * out["r_ttca_cys"]
        out["r_cama_form"] = _route.k_form(float(temperature_k))
        _decay = _route.k_decay(float(temperature_k), float(ph))
        out["r_cama_dpo"] = _route.dpo_share * _decay
        out["r_cama_tdp"] = (1.0 - _route.dpo_share) * _decay
        out["r_pent_dpo"] = 0.0
        out["r_pent_tdp"] = 0.0
    else:
        out["r_cama_form"] = out["r_cama_dpo"] = out["r_cama_tdp"] = 0.0

    # DERIVED RATE 0d (B59a): the Strecker sulfur released as H2S, a share of the Strecker step (Kobayasi 1965).
    _strecker = _ps.STRECKER_SULFUR
    if _strecker.enabled:
        out["r_dpo_ddp_h2s"] = float(_strecker.h2s_share) * out["r_dpo_ddp"]
        out["r_dpo_ddp"] = (1.0 - float(_strecker.h2s_share)) * out["r_dpo_ddp"]
    else:
        out["r_dpo_ddp_h2s"] = 0.0

    # DERIVED RATE 0c (B61): norfuraneol's hydroxide-catalysed decay, Mikami 2017's measured constant.
    _nf_oh = _ps.NF_HYDROXIDE_DECAY
    out["r_nf_decay_oh"] = (
        _nf_oh.k_ph5(float(temperature_k)) * ph_factor("hydroxide", ph, float(temperature_k), ionic_strength_m)
        if _nf_oh.enabled else 0.0
    )

    _k_adduct = thiol_adduct_equilibrium_l_per_mmol(float(temperature_k))
    out["ch_thioether_release"] = out["ch_thioether_mft"] / _k_adduct
    out["ch_thioether_release_fft"] = out["ch_thioether_fft"] / _k_adduct

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

#: C5 (2026-10-28): the same layout as arrays, for the inline-coupled integration's vectorised rates and
#: analytic mass-action Jacobian. Padded positions point at an appended 1.0 with exponent 0.
_MAX_REACTANTS = max(len(layout) for layout in _REACTANT_LAYOUT)
_R_IDX = np.full((len(_REACTANT_LAYOUT), _MAX_REACTANTS), N_SULFUR_STATE, dtype=int)
_R_POW = np.zeros((len(_REACTANT_LAYOUT), _MAX_REACTANTS), dtype=float)
for _j, _layout in enumerate(_REACTANT_LAYOUT):
    for _p, (_index, _coefficient) in enumerate(_layout):
        _R_IDX[_j, _p], _R_POW[_j, _p] = _index, float(_coefficient)
_R_ROWS = np.repeat(np.arange(len(_REACTANT_LAYOUT)), _MAX_REACTANTS).reshape(_R_IDX.shape)


def _rates_vectorised(y: np.ndarray, kv: np.ndarray) -> np.ndarray:
    """``sulfur_reaction_rates`` in array form (same values to rounding)."""
    terms = np.append(np.clip(y, 0.0, None), 1.0)[_R_IDX] ** _R_POW
    return kv * np.prod(terms, axis=1)


def _rate_jacobian(y: np.ndarray, kv: np.ndarray) -> np.ndarray:
    """d(rate_j)/d(y_i) of mass action at fixed constants, shape (n_reactions, n_species)."""
    base = np.append(np.clip(y, 0.0, None), 1.0)[_R_IDX]
    terms = base ** _R_POW
    out = np.zeros((len(_REACTANT_LAYOUT), N_SULFUR_STATE + 1))
    for p in range(_MAX_REACTANTS):
        others = np.prod(np.delete(terms, p, axis=1), axis=1)
        power = _R_POW[:, p]
        d = kv * power * np.where(power > 0.0, base[:, p] ** np.maximum(power - 1.0, 0.0), 0.0) * others
        np.add.at(out, (_R_ROWS[:, p], _R_IDX[:, p]), d)
    return out[:, :N_SULFUR_STATE]


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
    #: B2.2: the COMPUTED pH at every point of ``times_min``. It is ``None``
    #: for a clamped or prescribed-trajectory run, so that a caller cannot
    #: mistake an interpolation for a prediction.
    ph_series: Optional[np.ndarray] = None
    #: B60: (node_times, k_nodes) -- the rate vectors the solver interpolated between, so a flux
    #: audit (``reaction_flux_integrals``) uses exactly the constants the trajectory saw.
    rate_nodes: Optional[Tuple[np.ndarray, np.ndarray]] = None

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


#: B60 (2026-10-28): WHERE THE pH NODES SIT. "uniform" is the grid every wave through B58 ran on:
#: ``ph_nodes`` points equally spaced in time. It is not converged. In a cysteine + pentose pot the
#: in-situ pH falls most of its way in the first ~0.1 % of the run (Hofmann ribose pH 5, 145 C: 5.27
#: -> 4.47 within 0.07 of 20 min), while the first uniform node sits at duration/40, so every rate
#: constant through the transient is evaluated at the wrong pH; doubling the nodes moved that pot's MFT
#: 39x and its final cysteine 220x. "log_and_linear" adds ``ph_nodes`` points geometrically spaced from
#: 1e-5 x duration, which resolves the transient (41 vs 81 per family agree within 2 % on MFT and FFT
#: in the three probe pots; results/validation/kinetic_core_b60_node_convergence.json). Frozen waves
#: set "uniform" to reproduce their numbers. Pre-registration: kinetic_core_b60_prereg.md.
PH_NODE_GRID: str = "log_and_linear"


#: B60 amendment 1: HOW THE pH FIXED POINT IS ITERATED. "fixed_sweeps" is every wave through B58: plain
#: Picard sweeps, stopped after ``ph_iterations`` whatever the residual. In the Kang TTCA pots it
#: oscillates between two branches (final pH 5.52 / 6.55 / 5.56 / 6.37 after sweeps 2-5, residual ~1 pH)
#: and settles only after ~20-40 sweeps. "relaxed_converged" under-relaxes each sweep by
#: PH_FIXED_POINT_RELAXATION and iterates until the residual is below PH_FIXED_POINT_TOLERANCE (at least
#: ``ph_iterations`` sweeps, at most PH_FIXED_POINT_MAX_SWEEPS). omega = 0.7 reached < 1e-3 pH within ~8
#: sweeps on every probe system, where omega = 1 needed 40; the converged answer is the same for every omega
#: tried (1.0 / 0.7 / 0.5 / 0.3), so the damping selects no branch. The tolerance is 1e-3 pH, not 1e-4: the
#: residual floors at 3e-5 - 6e-4 pH from the ODE tolerances. A run left above PH_FIXED_POINT_FLAG is marked
#: ``ph_fixed_point_converged = False`` in its metadata.
PH_FIXED_POINT: str = "relaxed_converged"
PH_FIXED_POINT_RELAXATION: float = 0.7
PH_FIXED_POINT_TOLERANCE: float = 1e-3
PH_FIXED_POINT_MAX_SWEEPS: int = 40
PH_FIXED_POINT_FLAG: float = 1e-2


class PhFixedPointNotConverged(UserWarning):
    """The in-situ pH fixed point ended above ``PH_FIXED_POINT_FLAG``: the trajectory, and every pH-dependent
    rate integrated on it, is not the converged answer. Never silent (2026-10-28)."""


#: B60: HOW THE SHARE ROWS' FLUX BUDGET IS INTEGRATED. "closed" is ``reaction_flux_integrals`` on a
#: linear-and-log output grid with the run's own node constants (closes to < 1 %; see
#: ``sulfur_flux_budget``). "linear_rederived" is every wave through B58, which did not close; frozen
#: waves set it, with PH_NODE_GRID = "uniform" and PH_FIXED_POINT = "fixed_sweeps", to reproduce their numbers.
FLUX_BUDGET: str = "closed"


def ph_node_times(duration: float, n_nodes: int) -> np.ndarray:
    """The times at which the dynamic pH (and with it every rate constant) is evaluated."""
    linear = np.linspace(0.0, duration, n_nodes)
    if PH_NODE_GRID == "uniform":
        return linear
    if PH_NODE_GRID != "log_and_linear":
        raise ValueError(f"unknown PH_NODE_GRID {PH_NODE_GRID!r}")
    return np.unique(np.concatenate([linear, np.geomspace(duration * 1e-5, duration, n_nodes)]))


#: B64: the largest temperature step between two rate nodes along a recorded profile, in kelvin. The
#: rate vector is linear in time between nodes; across 2 K at 100 kJ/mol a rate changes by ~12 %, which
#: the straight line follows to ~0.2 %.
PROFILE_NODE_STEP_K: float = 2.0


def _checked_profile(
    temperature_profile: Optional[Tuple[Sequence[float], Sequence[float]]], duration: float
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """A ``(times_min, temperatures_K)`` profile as arrays, checked; None passes through."""
    if temperature_profile is None:
        return None
    times = np.asarray(temperature_profile[0], dtype=float)
    temps = np.asarray(temperature_profile[1], dtype=float)
    if times.ndim != 1 or times.shape != temps.shape or times.size < 2:
        raise ValueError("a temperature profile is two 1-D arrays of equal length >= 2")
    if times[0] != 0.0 or np.any(np.diff(times) <= 0.0):
        raise ValueError(
            "a temperature profile starts at t = 0 with strictly increasing times (no jumps: "
            "the rate vector is linear in time between nodes, so a jump would be smeared)")
    if abs(times[-1] - duration) > 1e-9 * max(1.0, duration):
        raise ValueError(f"the profile ends at {times[-1]:g} min, the run at {duration:g} min")
    if np.any(temps <= 0.0):
        raise ValueError("profile temperatures are in kelvin")
    return times, temps


def profile_node_times(times: np.ndarray, temps: np.ndarray, duration: float) -> np.ndarray:
    """Node times that cut every linear piece of a profile into steps of at most PROFILE_NODE_STEP_K."""
    out = [np.array([0.0, duration])]
    for t0, t1, c0, c1 in zip(times[:-1], times[1:], temps[:-1], temps[1:]):
        n = max(1, int(np.ceil(abs(c1 - c0) / PROFILE_NODE_STEP_K - 1e-12)))
        out.append(np.linspace(t0, t1, n + 1))
    return np.unique(np.clip(np.concatenate(out), 0.0, duration))


def rate_vector_at(node_times: np.ndarray, k_nodes: np.ndarray, t: float) -> np.ndarray:
    """The rate vector at time ``t``: linear in time between the two nodes that bracket it."""
    n_nodes = len(node_times)
    lo = min(max(int(np.searchsorted(node_times, t, side="right")) - 1, 0), n_nodes - 2)
    span = float(node_times[lo + 1] - node_times[lo])
    w = min(max((float(t) - float(node_times[lo])) / span, 0.0), 1.0) if span > 0.0 else 0.0
    return k_nodes[lo] * (1.0 - w) + k_nodes[lo + 1] * w


def _solve_on_time_nodes(
    k_nodes: np.ndarray,
    node_times: np.ndarray,
    y0: np.ndarray,
    t_eval: np.ndarray,
    duration: float,
    method: str,
    rtol: float,
    atol: float,
):
    """
    Integrate with a rate vector that is PIECEWISE LINEAR IN TIME between the
    supplied nodes. This is B2.1's pre-tabulation trick generalised from "the
    pH moves linearly between two measured endpoints" to "the pH follows an
    arbitrary computed profile", and it costs the same: one registry evaluation
    per node, not one per solver step.
    """
    from scipy.integrate import solve_ivp

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return sulfur_derivatives(y, rate_vector_at(node_times, k_nodes, t))

    return solve_ivp(
        rhs, (0.0, duration), y0, t_eval=t_eval, method=method, rtol=rtol, atol=atol
    )


#: C5 (2026-10-28): HOW THE DYNAMIC pH COUPLES TO THE CHEMISTRY. "fixed_point_nodes" is every wave
#: through B63: the pH is held at time nodes, the network integrated, the pH recomputed from the states,
#: and the sweep repeated to a fixed point (B60). "inline" solves the charge balance at every
#: right-hand-side evaluation from the current state, so the rate constants see the pH the state
#: implies at that instant: no node grid, no sweeps, and no way for either to be unconverged (the two
#: defects B60 found). Same physics; the answer is the limit the node fixed point converges to.
#: Frozen waves keep "fixed_point_nodes"; B64 runs on "inline" (kinetic_core_b64_prereg.md).
PH_COUPLING: str = "fixed_point_nodes"

#: The state keys the charge balance reads (``ph_state.titratable_inventory``).
_TITRATABLE_KEYS = ("ACID", "FA", "AA", "Cys", "ARP", "TTCA", "CBX", "CAMA")


def _integrate_inline_ph(parameters, y0, eval_times, duration, temperature_of, sid, drift, spec,
                         method, rtol, atol):
    """
    The inline-coupled integration. Returns (solution, ph_hot, ionic, ph_cool) at ``eval_times`` and the
    rate vectors there (for the flux audit). The pH solve is warm-started from the previous call;
    the Jacobian is mass action at the current constants (the pH's own dependence on the state is
    left out of it, which an implicit solver tolerates: it changes the step, not the answer).
    """
    from scipy.integrate import solve_ivp

    from . import ph_kernel
    from .ph_state import (solve_ph_and_ionic_strength_compiled, solve_ph_and_ionic_strength_warm,
                           titratable_inventory)

    warm_solve = solve_ph_and_ionic_strength_compiled if ph_kernel.COMPILED else solve_ph_and_ionic_strength_warm

    index = [SULFUR_INDEX[k] for k in _TITRATABLE_KEYS]
    use_label = bool(_ph.B51.label_scale_tables)
    memo = {"ph": None, "i": 0.0, "cool": None}

    def state_ph(y, temperature):
        inventory = titratable_inventory(
            {k: float(max(y[i], 0.0)) for k, i in zip(_TITRATABLE_KEYS, index)}, drift, spec)
        if memo["ph"] is None:
            from .ph_state import solve_ph_and_ionic_strength
            ph, ionic = solve_ph_and_ionic_strength(temperature, sid, inventory)
        else:
            ph, ionic = warm_solve(temperature, sid, inventory, memo["ph"], memo["i"])
        cool = None
        if use_label:
            guess = memo["cool"] if memo["cool"] is not None else ph
            cool = warm_solve(LABEL_TEMPERATURE_K, sid, inventory, guess, ionic)[0]
        memo.update(ph=ph, i=ionic, cool=cool)
        return ph, ionic, cool

    def constants(t, y):
        temperature = temperature_of(t)
        ph, ionic, cool = state_ph(y, temperature)
        return rate_vector(sulfur_rate_constants_at(
            parameters, temperature, float(ph), ph_label=cool, ionic_strength_m=float(ionic)))

    def rhs(t, y):
        return SULFUR_STOICHIOMETRY @ _rates_vectorised(y, constants(t, y))

    def jac(t, y):
        return SULFUR_STOICHIOMETRY @ _rate_jacobian(y, constants(t, y))

    kwargs = {"jac": jac} if method in ("BDF", "Radau", "LSODA") else {}
    solution = solve_ivp(rhs, (0.0, duration), y0, t_eval=eval_times, method=method,
                         rtol=rtol, atol=atol, **kwargs)
    if not solution.success:
        raise RuntimeError(f"sulfur-network integration failed: {solution.message}")
    states = np.clip(solution.y.T, 0.0, None)
    memo.update(ph=None, i=0.0, cool=None)
    ph_hot, ionic, ph_cool, k_rows = [], [], [], []
    for t, y in zip(eval_times, states):
        temperature = temperature_of(float(t))
        p, i, c = state_ph(y, temperature)
        ph_hot.append(p), ionic.append(i), ph_cool.append(c)
        k_rows.append(rate_vector(sulfur_rate_constants_at(
            parameters, temperature, float(p), ph_label=c, ionic_strength_m=float(i))))
    return solution, np.array(ph_hot), np.array(ionic), ph_cool, np.array(k_rows)


def integrate_sulfur(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    times_min: Sequence[float],
    *,
    ph: float = 5.0,
    ph_final: Optional[float] = None,
    buffer_spec: Optional[BufferSpec] = None,
    ph_drift: Optional[PhDrift] = None,
    method: str = "LSODA",
    rtol: float = 1e-9,
    atol: float = 1e-12,
    ph_nodes: int = 41,
    ph_iterations: int = 3,
    sid_mmol_l: Optional[float] = None,
    temperature_profile: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
) -> SulfurRun:
    """
    Integrate the sulfur network at one temperature and one pH TRAJECTORY.

    B64 (2026-10-28): or along a recorded TEMPERATURE PROFILE. ``temperature_profile`` is
    ``(times_min, temperatures_K)``, piecewise linear from t = 0; the rate constants (and the
    in-situ pH they depend on) are then evaluated at T(t) on the same time nodes the pH uses,
    with nodes added so that no node interval spans more than ``PROFILE_NODE_STEP_K``.
    ``temperature_k`` is kept for the run's metadata and its extrapolation warnings, and
    should be the profile's peak. Without a profile every number is exactly what it was.

    THREE pH MODES, and which one is in force is always recorded in the run's
    metadata so that a number can never be mistaken for a different kind of
    number:

    * **CLAMPED** (``ph_drift is None`` and ``ph_final is None``, or an explicit
      ``BufferSpec(kind="clamped")``). The pH is held at ``ph``. This is B2's
      behaviour and it remains the behaviour for any caller that does not ask
      for more.

    * **PRESCRIBED** (``ph_final`` given, ``ph_drift`` not given). B2.1's path:
      the pH moves linearly in time between a measured initial and a measured
      final value. It is INTERPOLATION BETWEEN TWO MEASUREMENTS, not a model --
      it cannot be run where nobody published the endpoint -- and it is kept
      only so that B2.1's reports remain reproducible.

    * **DYNAMIC** (``ph_drift`` given). B2.2: the pH is solved from a CHARGE
      BALANCE against the state at every point of the trajectory
      (``ph_state.ph_of_state``), driven by the organic-acid pool the network
      itself produces and resisted by the declared buffer capacity. This is the
      state ``FIT_HOLDOUT_DECLARATION.md`` sec.5 decision 1 says the Zhou pH-6
      column stays DIAGNOSTIC "until the model carries a pH-trajectory state".
      If ``ph_final`` is also supplied it is IGNORED and the fact is recorded:
      a measured endpoint is a thing to be PREDICTED once the state exists.

    The dynamic mode is a fixed point. The pH feeds the rate constants, the
    rate constants make acid, the acid feeds the pH. Because the feedback is
    weak in most pots it usually converges in a few sweeps under B60's relaxed
    iteration (up to 10 on the fit systems), but not always: a run whose residual
    stays above ``PH_FIXED_POINT_FLAG`` emits ``PhFixedPointNotConverged`` (and records
    ``ph_fixed_point_converged = False``), and the headline guards fail on any scored pot
    that does (independent review, 2026-10-28: a diagnostic nothing reads is decoration).
    """
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
    y0 = initial_sulfur_state(initial)
    t_k = float(temperature_k)
    profile = _checked_profile(temperature_profile, duration)

    def _with_profile_nodes(nodes: np.ndarray) -> np.ndarray:
        if profile is None:
            return nodes
        return np.unique(np.concatenate([nodes, profile_node_times(*profile, duration)]))

    def _node_temperatures(nodes: np.ndarray) -> np.ndarray:
        if profile is None:
            return np.full(len(nodes), t_k)
        return np.interp(nodes, profile[0], profile[1])

    spec = buffer_spec if buffer_spec is not None else DEFAULT_BUFFER
    dynamic = ph_drift is not None and not spec.is_clamped
    ph_notes: list = []
    if buffer_spec is None and ph_drift is not None:
        ph_notes.append(BUFFER_ABSENT_WARNING)
    if dynamic and ph_final is not None:
        ph_notes.append(
            "a measured final pH was supplied AND the dynamic pH state is "
            "active; the measured value was IGNORED and the endpoint is a "
            "PREDICTION."
        )

    ph_series: Optional[np.ndarray] = None
    ph_convergence: Optional[float] = None
    sid: Optional[float] = None

    if dynamic:
        assert ph_drift is not None
        node_times = _with_profile_nodes(ph_node_times(duration, max(int(ph_nodes), 5)))
        node_tk = _node_temperatures(node_times)
        n_nodes = len(node_times)
        # every requested output time, plus the internal nodes, in one sorted
        # evaluation grid so that no extra integration is needed
        eval_times = np.unique(np.concatenate([grid, node_times]))
        initial_map = {
            key: float(y0[SULFUR_INDEX[key]]) for key in SULFUR_STATE_KEYS
        }
        # ``ph`` is the DECLARED LABEL, read on a cooled sample. It is used
        # exactly once -- here -- and never again; from this point on the pH is
        # solved from the conserved SID at whatever temperature is being asked
        # about.
        # B51: a later thermal segment of the SAME pot passes the SID its first segment computed,
        # so acid made earlier is not re-titrated away at the segment boundary.
        if sid_mmol_l is not None and _ph.B51.sid_across_segments:
            sid = float(sid_mmol_l)
        else:
            sid = strong_ion_difference(ph0, initial_map, ph_drift, spec)

        def _node_state(state_map, temperature):
            p_hot, i_hot = ph_and_ionic_strength_of_state(state_map, temperature, sid, ph_drift, spec)
            p_cool = (ph_of_state(state_map, LABEL_TEMPERATURE_K, sid, ph_drift, spec)
                      if _ph.B51.label_scale_tables else None)
            return p_hot, i_hot, p_cool

        if profile is None or PH_COUPLING == "inline":
            p0, i0, c0 = _node_state(initial_map, float(node_tk[0]))
            node_ph = np.full(n_nodes, p0)
            node_i = np.full(n_nodes, i0)
            node_cool = [c0] * n_nodes
        else:
            # B64: the starting guess at each node is the initial state's pH at that node's T.
            starts = [_node_state(initial_map, float(tk)) for tk in node_tk]
            node_ph = np.array([s[0] for s in starts])
            node_i = np.array([s[1] for s in starts])
            node_cool = [s[2] for s in starts]
        solution = None
        relaxed = PH_FIXED_POINT == "relaxed_converged"
        if not relaxed and PH_FIXED_POINT != "fixed_sweeps":
            raise ValueError(f"unknown PH_FIXED_POINT {PH_FIXED_POINT!r}")
        if PH_COUPLING not in ("fixed_point_nodes", "inline"):
            raise ValueError(f"unknown PH_COUPLING {PH_COUPLING!r}")
        min_sweeps = max(int(ph_iterations), 1)
        max_sweeps = max(min_sweeps, PH_FIXED_POINT_MAX_SWEEPS) if relaxed else min_sweeps
        sweeps_run = 0
        if PH_COUPLING == "inline":
            # The audit grid is the node grid's times; nothing is held there, it only records.
            if profile is None:
                def temperature_of(t):
                    return t_k
            else:
                def temperature_of(t):
                    return float(np.interp(t, profile[0], profile[1]))
            solution, node_ph, node_i, node_cool, k_nodes = _integrate_inline_ph(
                parameters, y0, eval_times, duration, temperature_of, sid, ph_drift, spec,
                method, rtol, atol)
            node_times = eval_times
            n_nodes = len(node_times)
            ph_convergence = 0.0
            max_sweeps = 0
        for sweep in range(max_sweeps):
            sweeps_run = sweep + 1
            k_nodes = np.array([
                rate_vector(sulfur_rate_constants_at(
                    parameters, float(tk), float(p), ph_label=c, ionic_strength_m=float(i)))
                for p, i, c, tk in zip(node_ph, node_i, node_cool, node_tk)
            ])
            solution = _solve_on_time_nodes(
                k_nodes, node_times, y0, eval_times, duration, method, rtol, atol
            )
            if not solution.success:
                raise RuntimeError(
                    f"sulfur-network integration failed: {solution.message}"
                )
            states = np.clip(solution.y.T, 0.0, None)
            fresh = np.empty(n_nodes)
            for j, t_node in enumerate(node_times):
                index = int(np.searchsorted(eval_times, t_node))
                index = min(index, len(eval_times) - 1)
                row = states[index]
                fresh[j], node_i[j], node_cool[j] = _node_state(
                    {k: float(row[SULFUR_INDEX[k]]) for k in SULFUR_STATE_KEYS}, float(node_tk[j]))
            ph_convergence = float(np.max(np.abs(fresh - node_ph)))
            if relaxed:
                node_ph = node_ph + PH_FIXED_POINT_RELAXATION * (fresh - node_ph)
                if ph_convergence < PH_FIXED_POINT_TOLERANCE and sweeps_run >= min_sweeps:
                    break
            else:
                node_ph = fresh
                if ph_convergence < 1e-4:
                    break
        assert solution is not None
        keep = np.searchsorted(eval_times, grid)
        keep = np.clip(keep, 0, len(eval_times) - 1)
        raw = np.clip(solution.y.T, 0.0, None)[keep]
        worst_negative = float(np.min(solution.y)) if solution.y.size else 0.0
        ph_series = np.interp(grid, node_times, node_ph)
        ph1 = float(ph_series[-1])
    else:
        ph1 = ph0 if ph_final is None else float(ph_final)
        if profile is not None:
            # B64: clamped or prescribed pH along a recorded profile; the pH, where it moves, is
            # B2.1's straight line in time between the two measured endpoints.
            n_nodes = max(int(ph_nodes), 5)
            node_times = _with_profile_nodes(np.linspace(0.0, duration, n_nodes))
            node_ph_line = ph0 + (ph1 - ph0) * node_times / duration
            k_nodes = np.array([
                rate_vector(sulfur_rate_constants_at(parameters, float(tk), float(p)))
                for tk, p in zip(_node_temperatures(node_times), node_ph_line)
            ])
        elif ph1 == ph0:
            k_fixed = rate_vector(sulfur_rate_constants_at(parameters, t_k, ph0))
            node_times = np.array([0.0, duration])
            k_nodes = np.array([k_fixed, k_fixed])
        else:
            # B2.1's path: the pH is smooth and enters only through closed-form
            # factors, so the constants are pre-tabulated and interpolated.
            n_nodes = max(int(ph_nodes), 5)
            node_times = np.linspace(0.0, duration, n_nodes)
            k_nodes = np.array([
                rate_vector(sulfur_rate_constants_at(parameters, t_k, float(p)))
                for p in np.linspace(ph0, ph1, n_nodes)
            ])
        solution = _solve_on_time_nodes(
            k_nodes, node_times, y0, grid, duration, method, rtol, atol
        )
        if not solution.success:
            raise RuntimeError(
                f"sulfur-network integration failed: {solution.message}"
            )
        raw = solution.y.T
        worst_negative = float(np.min(raw)) if raw.size else 0.0

    if worst_negative < -max(1e3 * atol, 1e-8):
        raise RuntimeError(
            f"sulfur-network integration produced a state of {worst_negative:.3e}, "
            f"far below the absolute tolerance {atol:.1e}: a genuine "
            f"non-negativity failure, not solver noise."
        )
    concentrations = np.clip(raw, 0.0, None)

    temperature_c = t_k - 273.15
    metadata: Dict[str, Any] = {
        "method": method,
        "temperature_C": temperature_c,
        "temperature_K": t_k,
        "temperature_profile_K": (
            None if profile is None
            else [[float(a), float(b)] for a, b in zip(profile[0], profile[1])]
        ),
        "ph_initial": ph0,
        "ph_final": ph1,
        "ph_is_clamped": (not dynamic) and ph1 == ph0,
        "ph_mode": (
            "dynamic_charge_balance" if dynamic
            else ("prescribed_measured_endpoints" if ph1 != ph0 else "clamped")
        ),
        "ph_notes": ph_notes,
        "buffer": spec.as_dict(),
        "n_species": N_SULFUR_STATE,
        "n_reactions": len(FULL_REACTIONS),
        "species_order": list(SULFUR_STATE_KEYS),
        "worst_raw_negative_before_clip": worst_negative,
        "extrapolation_warnings": _sulfur_extrapolation_warnings(
            parameters, temperature_c, ph0
        ),
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
    }
    if dynamic:
        assert ph_drift is not None
        final_map = {
            key: float(concentrations[-1, SULFUR_INDEX[key]])
            for key in SULFUR_STATE_KEYS
        }
        ph_final_cooled = ph_of_state(
            final_map, LABEL_TEMPERATURE_K, float(sid), ph_drift, spec
        )
        metadata.update({
            "ph_drift_constants": ph_drift.as_dict(),
            "strong_ion_difference_mmol_L": float(sid) if sid is not None else None,
            "ph_label_declared": ph0,
            "ph_label_temperature_K": LABEL_TEMPERATURE_K,
            "ph_initial_in_situ": float(ph_series[0]),  # type: ignore[index]
            "ph_final_in_situ": ph1,
            #: THE ONLY QUANTITY COMPARABLE WITH A PUBLISHED ENDPOINT.
            "ph_final_cooled": float(ph_final_cooled),
            "ph_series": [float(v) for v in ph_series],  # type: ignore[union-attr]
            "ph_fixed_point_residual": ph_convergence,
            "ph_fixed_point_sweeps": sweeps_run,
            "ph_fixed_point_converged": bool(ph_convergence <= PH_FIXED_POINT_FLAG),
            "ph_coupling": PH_COUPLING,
            "ph_numerics_b60": {"node_grid": PH_NODE_GRID, "fixed_point": PH_FIXED_POINT,
                                "relaxation": PH_FIXED_POINT_RELAXATION, "n_nodes": int(n_nodes)},
            "ionic_strength_M_final": float(node_i[-1]),
            "ph_state_corrections_b51": _ph.B51.as_dict(),
            "ph_drop_units_cooled": float(ph0 - ph_final_cooled),
            "buffer_capacity_initial_mmol_per_pH": buffer_capacity_mmol_per_ph(
                {key: float(y0[SULFUR_INDEX[key]]) for key in SULFUR_STATE_KEYS},
                t_k, ph0, ph_drift, spec),
            "buffer_capacity_final_mmol_per_pH": buffer_capacity_mmol_per_ph(
                final_map, t_k, ph1, ph_drift, spec),
            "acid_pool_final_mmol_L": final_map.get("ACID", 0.0),
        })
        if ph_convergence > PH_FIXED_POINT_FLAG:
            import warnings

            warnings.warn(
                f"pH fixed point not converged: residual {ph_convergence:.3g} pH after {sweeps_run} sweeps "
                f"(flag {PH_FIXED_POINT_FLAG:g})", PhFixedPointNotConverged, stacklevel=2)
    metadata.update(sulfur_registry_metadata(
        {k: v for k, v in parameters.items() if isinstance(v, SulfurParameter)}
    ))

    run = SulfurRun(
        times_min=grid,
        concentrations=concentrations,
        temperature_k=t_k,
        ph_initial=ph0,
        ph_final=ph1,
        metadata=metadata,
        ph_series=ph_series,
        rate_nodes=(np.asarray(node_times, dtype=float), np.asarray(k_nodes, dtype=float)),
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
    buffer_spec: Optional[BufferSpec] = None,
    ph_drift: Optional[PhDrift] = None,
    n_points: int = 201,
) -> Dict[str, float]:
    """Time-integrated flux through every reaction, mmol/L.

    B60: the rates are integrated on the union of a linear and a log-spaced output grid (a
    cysteine + pentose transient lasts seconds; a linear grid of 101 points over 20 min does not see
    it) and evaluated with the run's own node constants (``SulfurRun.rate_nodes``), so they are the
    rates the trajectory saw -- B51's ionic strength and label-scale pH included. Before B60 the
    budget re-derived the constants from the pH series alone and did not close (in Cao 2017's pot it
    booked 2 575 mM of 3-deoxypentosone made from 400 mM of xylose).
    """
    duration = float(duration_min)
    if FLUX_BUDGET == "linear_rederived":
        return _legacy_flux_budget(parameters, temperature_k, initial, duration, ph, ph_final,
                                   buffer_spec, ph_drift, n_points)
    if FLUX_BUDGET != "closed":
        raise ValueError(f"unknown FLUX_BUDGET {FLUX_BUDGET!r}")
    grid = np.unique(np.concatenate([
        np.linspace(0.0, duration, int(n_points)),
        np.geomspace(duration * 1e-7, duration, int(n_points)),
    ]))
    run = integrate_sulfur(
        parameters, temperature_k, initial, grid, ph=ph, ph_final=ph_final,
        buffer_spec=buffer_spec, ph_drift=ph_drift,
    )
    return reaction_flux_integrals(run)


def _legacy_flux_budget(parameters, temperature_k, initial, duration, ph, ph_final, buffer_spec,
                        ph_drift, n_points) -> Dict[str, float]:
    """The budget every wave through B58 scored its share rows on (``FLUX_BUDGET = "linear_rederived"``):
    a linear output grid, constants re-derived from the pH series alone. Kept only so frozen waves
    reproduce their numbers; it does not close."""
    grid = np.linspace(0.0, duration, int(n_points))
    run = integrate_sulfur(
        parameters, temperature_k, initial, grid, ph=ph, ph_final=ph_final,
        buffer_spec=buffer_spec, ph_drift=ph_drift,
    )
    ph0, ph1 = run.ph_initial, run.ph_final
    rows = []
    for i, t in enumerate(grid):
        if run.ph_series is not None:
            ph_here = float(run.ph_series[i])
        else:
            frac = min(max(t / max(grid[-1], 1e-12), 0.0), 1.0)
            ph_here = ph0 + frac * (ph1 - ph0)
        k_t = sulfur_rate_constants_at(parameters, float(temperature_k), ph_here)
        rows.append(sulfur_reaction_rates(run.concentrations[i], k_t))
    integral = np.trapezoid(np.array(rows), grid, axis=0)
    return {r.key: float(integral[j]) for j, r in enumerate(FULL_REACTIONS)}


def reaction_flux_integrals(run: "SulfurRun") -> Dict[str, float]:
    """Integral over the run of every reaction's rate, mmol/L, by the trapezoid rule on the run's
    own output grid with the constants the solver used. Only as accurate as that grid is dense:
    ``flux_closure`` says how well it closes."""
    if run.rate_nodes is None:
        raise ValueError("the run carries no rate nodes")
    node_times, k_nodes = run.rate_nodes
    rows = [
        sulfur_reaction_rates(run.concentrations[i], rate_vector_at(node_times, k_nodes, float(t)))
        for i, t in enumerate(run.times_min)
    ]
    integral = np.trapezoid(np.array(rows), run.times_min, axis=0)
    return {r.key: float(integral[j]) for j, r in enumerate(FULL_REACTIONS)}


def flux_closure(run: "SulfurRun", flux: Mapping[str, float], species: Sequence[str]) -> Dict[str, float]:
    """For each species, |(final - initial) - net integrated production| over the larger of the
    species' initial amount and its gross turnover: the relative error of the flux budget."""
    out: Dict[str, float] = {}
    for key in species:
        net, gross = 0.0, 0.0
        for rxn in FULL_REACTIONS:
            nu = rxn.products.get(key, 0) - rxn.reactants.get(key, 0)
            if nu:
                net += nu * flux[rxn.key]
                gross += abs(nu) * abs(flux[rxn.key])
        change = run.final(key) - float(run.concentrations[0, SULFUR_INDEX[key]])
        scale = max(float(run.concentrations[0, SULFUR_INDEX[key]]), gross, 1e-30)
        out[key] = abs(change - net) / scale
    return out


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
    # B2.1: THE INTACT-SKELETON ROUTE IS NOW TWO PARALLEL STEPS, one through
    # neutral H2S and one through HS-. They are the same ROUTE -- same carbon
    # skeleton, same precursor, same product -- differing only in which
    # protonation state of the sulfide does the attack, so a route share that
    # counted one and not the other would be wrong twice over: it would
    # under-report the intact-skeleton share, and it would make that share move
    # with pH for a reason that has nothing to do with the route mix.
    mft_from_intact = flux.get("r_ddp_mft", 0.0) + flux.get("r_ddp_mft_hs", 0.0)
    mft_from_nf = flux.get("r_nf_mft", 0.0)
    mft_from_c2c3 = flux.get("r_ha_mp_mft", 0.0)
    mft_total = mft_from_thiamine + mft_from_intact + mft_from_nf + mft_from_c2c3

    fur_to_fft = flux.get("r_fur_fft", 0.0) + flux.get("r_fur_fft_hs", 0.0)
    fur_total = fur_to_fft + flux.get("r_fur_decay", 0.0)

    mft_made = mft_total
    mft_to_dimer = 2.0 * flux.get("ch_dimer_mft", 0.0)
    mft_to_mmft = flux.get("ch_mmft", 0.0)
    mft_to_thioether = flux.get("ch_thioether_mft", 0.0)
    mft_to_oligomer = flux.get("ch_oligomer_mft", 0.0)
    mft_to_decay = flux.get("ch_mft_decay", 0.0)
    mft_to_thiolate = flux.get("ch_thiolate_loss_mft", 0.0)
    mft_to_protein = flux.get("ch_protein_ss_mft", 0.0)

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
        "MFT_consumed_share_thiolate_oxidation": share(mft_to_thiolate, mft_made),
        "MFT_consumed_share_protein_disulfide": share(mft_to_protein, mft_made),
        "MFT_share_intact_skeleton_via_neutral_H2S": share(
            flux.get("r_ddp_mft", 0.0), mft_from_intact
        ),
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
                # B2.3: the FOURTH invariant, per step. Zero for both centres
                # means the step moves no titratable capacity; anything else
                # must match CENTRE_LEDGER or the module will not import.
                "titratable_centre_delta": centre_delta(r.reactants, r.products),
                "centre_ledger_basis": (
                    CENTRE_LEDGER.get(r.key, {}).get("basis")),
                "note": r.note,
            }
            for r in FULL_REACTIONS
        ],
        "site_pools_zero_atom": list(SITE_POOLS),
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
        "reference_temperature_K": T_REF_S_K,
        "measured_sulfur_parameters": sorted(MEASURED_SULFUR),
    }
