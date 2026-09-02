#!/usr/bin/env python3
"""
scripts/generators/generate_kinetic_core_b2_3_fit.py

THE FIT STAGE OF BUILD WAVE B2.3 (kinetic core, Modules 1-2: sulfur formation
and thiol consumption).

WHAT B2.3 CHANGES, AND WHY -- ONE INVARIANT, NO NEW PARAMETER
===========================================================================
B2.3 IS A REFIT AT UNCHANGED STRUCTURE-MINUS-ONE-BUG, ON THE UNCHANGED
DECLARED FIT ROWS, WITH THE SAME OPTIMISER SETTINGS. It adds NO free
parameter, NO rate constant, NO species that any rate law can see, and it
touches no fit row.

What it fixes is a CONSERVATION DEFECT that B2.2's own post-freeze diagnosis
(`results/validation/kinetic_core_b2_2_diagnosis.md` sec. 3) found and
deliberately did NOT apply, because applying a stoichiometry change after
reading a hold-out scorecard is a refit in disguise. The defect: sink steps
delete a carboxylate- or amine-bearing solute without saying where the
titratable centre went, so the strong ion that had titrated it balances
nothing and the pot drifts alkaline for a bookkeeping reason. Kang's 120 C
pot ran to a predicted pH 11.4 against a measured 4.9.

FIT_HOLDOUT_DECLARATION.md Amendment 9 pre-declares this wave, including its
TRIGGER DISCLOSURE: the defect was noticed via hold-out SCORING, so the Kang
140 C MFT row is demoted to seen-diagnostic PERMANENTLY and its B2.2 nominal
pass is ruled not-counted. The fix itself is a conservation law and is
defensible with no reference to any hold-out value.

ONE new TERMINAL accounting pool carries the half of the defect that is a
certainty: `CBX`, carried carboxylate equivalents. The OTHER half -- the
liberated amino nitrogen -- is DECLARED NON-TITRATABLE at the point of
release rather than carried, because carrying it is arithmetically
incompatible with Amendment 7's own FIT anchors (40 mmol/L of amino nitrogen
finishing at a measured pH 3.22 / 3.42). That asymmetry is the wave's main
scientific claim and it is stated, evidenced and error-bounded in
`ph_state.AMINE_FATE_BASIS`. Nothing consumes CBX, no rate law reads it, and
the only observable it can move is the pH. A fourth import-time invariant -- `ph_state.validate_charge_closure`,
alongside carbon, nitrogen and sulfur -- now refuses the import unless every
step's titratable-centre movement is DECLARED in the network's
`CENTRE_LEDGER`, which is what makes the defect unreachable rather than
merely repaired.

WHAT B2.2 CHANGED, AND WHY -- BOTH CHANGES WERE B2.1's OWN DIAGNOSIS, ACTED ON
===========================================================================
1. **THE DECAY LUMPS GET THEIR OWN BARRIERS.** B2.1 gave every `*_decay` lump
   the LUMPED FORMATION Ea. That fixed eleven low-temperature rows and broke
   the top of the ladder: `kinetic_core_b2_1_diagnosis.md` sec. 2a measured the
   140 C pot going precursor-exhausted while sinks that accelerate exactly as
   fast as their source keep running, and its closing sentence is "the
   `*_decay` lumps should not share the FORMATION barrier. They need their
   own". They now have two, in two NAMED FAMILIES chosen so that the declared
   FIT rows can identify them (parameters_sulfur.DECAY_FAMILIES):
     * `thiol_sink` -- k_mft_decay, k_fft_decay, k_thiol_decay, k_thiolate_loss.
       Identified by Kumazawa 2003's 121 C survival grid (a FORMATION-FREE read
       of the sink), Hofmann's 145 C pH-5 anchors and Kang's declared 100/120 C
       rungs.
     * `carbonyl_sink` -- k_fur_decay, k_osone_decay, k_nf_decay. Identified by
       Kang's furfural ladder at 100/120 C and Hofmann T5's 145 C in-situ
       furfural and norfuraneol.
   `k_dimer_decay` and `k_h2s_loss` deliberately KEEP the formation Ea: no fit
   row measures a disulfide or a sulfide sink at two temperatures, so a barrier
   for either would be identified by nothing.

2. **THE pH IS A COMPUTED STATE.** B2.1 moved the pH linearly between a
   MEASURED initial and a MEASURED final value -- interpolation between two
   answers, unusable anywhere nobody published the endpoint. B2.2 solves the pH
   from a CHARGE BALANCE against the state (`src/kinetic_core/ph_state.py`):
   the organic-acid pool the network itself produces, against the declared
   buffer capacity, at textbook pKa corrected to reaction temperature. TWO
   constants are calibrated and they are calibrated on the THREE anchors
   FIT_HOLDOUT_DECLARATION.md Amendment 7 declares FIT for exactly this
   purpose -- Zhou 2023 Fig. 2's ARP+Cys endpoints, which are PROCESS-STATE
   observables and not volatile levels.

NET NEW FREE PARAMETERS: FOUR (two decay barriers, two drift constants).

Estimates the sulfur network's rate constants -- NONE of which has a literature
value -- by least squares on the DECLARED FIT ROWS ONLY, and writes
`results/validation/kinetic_core_b2_3_fit_report.{json,md}` together with the
FROZEN parameter set that the hold-out scorer reads.

===========================================================================
THE FIT CORPUS -- and nothing else
===========================================================================
Every row below is marked FIT in docs/reference/FIT_HOLDOUT_DECLARATION.md
D.3 (Module 1), D.4 (Module 2) or D.6 (the Module 9 rows the declaration
assigns to the sulfur lane). Each carries its source anchor in the row itself.

  * Hofmann & Schieberle 1998 Table 1/2, **pH 5.0 rows only** -- the SIDA
    absolute anchors (ribose, xylose, glucose, fructose)
  * Hofmann 1998 Tables 3/4/6/7/8/10 -- the step-level fed-precursor yields
  * Hofmann 1998 Table 5 -- the in-situ norfuraneol and furfural levels
  * Zhou 2023 Table 1, **pH 7 column only** (ARP and ARP+Cys)
  * Whitfield 1999 -- fed norfuraneol at pH 4.5
  * van Seeventer 2001 -- the precursor-conversion pair (55% / 75%)
  * Zhang 2024 Figure 1 -- the four-additive redox series, as RATIOS
  * Cerny 2007b Table 3 -- the full-ternary 54:46 MFT isotope split
  * Cerny 2007 Table 4 -- the per-species isotope fractions, EXCEPT its MFT
    column (see THE DISJOINTNESS PROBLEM below)
  * Yaghmur 2005 -- the FFT-share ceiling, as a one-sided constraint

WHAT BUILD WAVE B2.1 ADDS TO THE FIT CORPUS
-------------------------------------------
Three blocks, all declared FIT before this wave ran:

  * **Kang 2026 SI Table S4, the 100 and 120 C rungs** (Amendment 5) -- MFT,
    FFT and furfural in a purified-TTCA pot, Tier A with printed calibration
    curves. This is the FIRST TEMPERATURE AXIS INSIDE A SINGLE SYSTEM the
    sulfur module has ever had: B2's fit panel was one temperature per system
    with the feedstock changing alongside it, which is why its single lumped Ea
    landed on its bound and was reported as a failure to determine.
  * **Kang 2026 SI Fig. S4, the free-cysteine ladder at 100 and 120 C**
    (Amendment 5) -- plus the 16.3 mol% TTCA -> free-cysteine yield CEILING as
    a one-sided bound, and Kang's measured Ea of 55.1 kJ/mol carried as a FIXED
    barrier on r_cys_thermal rather than as a fitted shape.
  * **Kumazawa 2003 Figure 3, six of seven pH levels at 121 C** (Amendment 4)
    -- 2-furfurylthiol survival with NO FORMATION PRESENT. It is the only
    measurement in the corpus of a thiol's pH response that cannot be confounded
    with the formation lane, and it is what pays for the thiolate-mediated
    consumption channel B2 did not have.

WHAT IS NEVER READ
------------------
No file under `data/benchmarks/` and no frozen `mp_holdout_*` bundle is opened
by this script. The hold-out VALUES do not appear in it: Zhou's pH-6 and pH-8
columns, Whitfield 2001's pH-6.5 collapse, Cerny 2003's 95 C isotope shares,
Hofmann's dry-180 rows, Hofmann's pH-3/pH-7 aqueous rows, Meynier's directional
rows, van Seeventer's 50 C zero-order rates, Zhang's Figure 2 fractions, Zhou's
120 C dimer fractions, Cerny 2007's 1x/2x pair, **Kang's 140 C column and its
140 C cysteine rung, and Sun 2019's pH-9 column** live ONLY in the hold-out
scorer. `tests/unit/test_kinetic_core_b2_3.py` greps this file and the
`src/kinetic_core/` sulfur modules for those literals, and separately walks this
script's own SYSTEMS dictionary asserting that no INTEGRATED CONDITION is a
hold-out condition -- which is the stronger check, because a value can leak
through a system definition without its literal ever being typed.

FIREWALL DISCLOSURE, because a hold-out whose exposure is undisclosed is worth
less than one whose exposure is stated: `kang2026_SI_extraction.md` PRINTS the
140 C column, and the build brief directed this wave to read that dossier --
it is also where the 100 and 120 C FIT rows live. Those hold-out values were
therefore SEEN before this fit ran. The same is true of Sun 2019's pH-9 column,
which the dossier prints beside the pH-3 column.

===========================================================================
THE DISJOINTNESS PROBLEM THIS SCRIPT FOUND, AND HOW IT RESOLVED IT
===========================================================================
The declaration puts **Cerny 2007 Table 4 isotope splits across pH** in the FIT
column and **Cerny 2007 Table 5's 1x/2x concentration pair** in the HOLD-OUT
column. But Table 4's MFT entry at pH 5 is "85% thiamine-derived", and Table
5's 1x arm IS "85 : 15 thiamine : xylose". The two columns therefore share a
number.

RESOLUTION, applied here and reported: **Table 4's MFT column is EXCLUDED from
the objective.** The fit uses only Cerny 2007b's unambiguous full-ternary
54 : 46 (a different paper, a different system composition, declared FIT) plus
Table 4's OTHER species -- FFT, 2-mercapto-3-pentanone, 3-mercapto-2-pentanone
-- which are not duplicated anywhere in the hold-out column. This is a
CONTRADICTION IN THE DECLARATION, not a modelling choice, and it is reported as
one rather than resolved silently in either direction.

===========================================================================
WHAT IS FITTED, AND WHAT IS NOT
===========================================================================
NOT FITTED (measured, fixed):
    k_cys_h2s     Zheng & Ho 1994's pH-resolved cysteine thermolysis set, as
                  the SOURCE-CONSISTENT MATCHED (Ea, A) PAIRS of
                  `data/lit/extraction_dossiers/zheng1994_extraction.md` sec. 5b
                  -- (131.2, 9.79e11), (133.0, 1.93e12), (135.5, 2.36e13),
                  (123.3, 1.04e12) in kJ/mol and 1/s at pH 3/5/7/9. These
                  reproduce the paper's own Table I everywhere within 3x, which
                  is unit-tested. The Cantera lane's row (which
                  shipped 130.4 kJ/mol with A = 1.0e14 1/s until Wave B8
                  repaired it to the matched pH-5 pair on 2026-08-30) is NOT
                  used here and never was: that pair ran ~15x faster than its
                  own source at pH 7 / 100 C. pH 3.0 is treated as a
                  SEPARATE MECHANISM and is not interpolated through, on the
                  authors' own statement and on their own data, which invert
                  the pH ordering there at 100 C.
    k_thioether   Hofmann 2002 / Charles-Bernard 2005, recast bimolecular on
                  the titrated site density. NO activation energy, by policy.
    k_oligomer    ZERO. Its only measurement is a declared HOLD-OUT.
    every pH shape (acid/base catalysis, H2S speciation) -- structural, not
                  fitted; there are ZERO fitted pH parameters.

FITTED: the log10 of every remaining rate constant at 145 C, plus ONE lumped
activation energy shared by the whole FORMATION lane. The CONSUMPTION lane gets
no Ea at all -- inventory sec. C.1 / B.7, and it is a prohibited derivation.

Each is started from RANDOM points inside deliberately wide bounds, so any
agreement with a literature number is a result and not an initialisation.

Usage:
    python scripts/generators/generate_kinetic_core_b2_3_fit.py [--starts N] [--quick]
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.kinetic_core import operative_parameters  # noqa: E402
from src.kinetic_core.ph_state import (  # noqa: E402
    ACID_YIELD_BOUNDS,
    ARP_AMINE_PKA_BOUNDS,
    PH_STATE_PROVENANCE,
    BufferSpec,
    PhDrift,
)
from src.kinetic_core.parameters_sulfur import (  # noqa: E402
    CONSUMPTION_KEYS,
    DECAY_EA_BOUNDS,
    DECAY_FAMILIES,
    DECAY_FAMILY_IDENTIFYING_ROWS,
    DECAY_KEYS_ON_FORMATION_EA,
    FITTED_SULFUR_BOUNDS_LOG10K,
    FITTED_SULFUR_KEYS,
    LUMPED_FORMATION_EA_BOUNDS,
    KANG_TTCA_CEILING_ANCHOR,
    KANG_TTCA_FREE_CYS_YIELD_CEILING_MOL_PCT,
    MEASURED_SULFUR,
    NO_EA_KEYS,
    T_REF_S_K,
    YAGHMUR_CEILING_ANCHOR,
    YAGHMUR_FFT_SHARE_CEILING,
    sulfur_registry_metadata,
    with_fitted_sulfur,
)
from src.kinetic_core.species_sulfur import (  # noqa: E402
    mmol_per_litre_to_ug_per_litre,
    ug_per_litre_to_mmol_per_litre,
)
from src.kinetic_core.sulfur import (  # noqa: E402
    OUT_OF_SCOPE,
    branch_shares,
    describe_sulfur,
    integrate_sulfur,
    sulfur_flux_budget,
)

# The B1 trunk fit, frozen, from results/validation/kinetic_core_b1_fit_report.json.
# The sulfur module inherits it unchanged: B2 does not re-fit any trunk constant.
B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

CELSIUS = 273.15

#: An "ambient aerobic oxidant charge", mmol/L of oxidant equivalents. It is a
#: MODELLING CONVENTION, not a measurement: the oxidant pool is an abstract
#: equivalent (zero atoms) and only the PRODUCT k_dimer * [OX] is identifiable,
#: so this constant is absorbed into the fitted k_dimer and moves nothing. It is
#: stated rather than hidden because the ratio between the ambient charge and a
#: cystine charge IS identifiable, and that ratio is what Zhang Fig. 1 measures.
OX_AMBIENT_MMOL_L = 1.0

#: Charles-Bernard's titrated site density x a nominal matrix loading. Only the
#: model-system rows that CONTAIN a matrix carry any MELE at all; the Hofmann,
#: Cerny, Whitfield and Zhou model solutions contain none, which is why the
#: thioether channel is silent in them and the dimerisation channel is not.
MELE_COFFEE_MMOL_L = 9.0 * 12.5  # 9 mmol/g x 12.5 g/L melanoidin = 112.5 mmol/L


# ===========================================================================
# THE DECLARED BUFFERS -- B2.2, and every one is a quotation
# ===========================================================================
# The buffer is now an INPUT to the integrator, not an assumption inside it.
# Each spec below is what the SOURCE says, and nothing else. Where a source
# says nothing, the system gets BUFFER_NONE and the run carries the declared
# extrapolation warning -- it is never silently clamped.
BUFFER_HOFMANN = BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=True,
    source=("Hofmann & Schieberle 1998 Methods: 0.5 mol/L phosphate buffer, "
            "quoted in zhou2023_extraction.md sec. 3 ('Hofmann 1998's pH "
            "series... used 0.5 mol/L phosphate buffer and therefore held pH')"))

#: ⚠️ A DECLARED ASSUMPTION, NOT A QUOTATION. Cerny 2007b, Whitfield 1999,
#: Zhang 2024 and van Seeventer 2001 are buffered aqueous model solutions of
#: the same class as Hofmann's, but NO dossier this wave read prints their
#: buffer or its molarity. They are given Hofmann's spec because that is the
#: only buffered composition the corpus states, and the assumption is carried
#: in the source string so that it travels into the fit report rather than
#: living in a reviewer's memory. It is NUMERICALLY IDENTICAL to
#: BUFFER_HOFMANN; only the provenance differs.
BUFFER_PHOSPHATE_ASSUMED = BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=False,
    source=("ASSUMED, NOT QUOTED: no dossier read by this wave states the "
            "buffer of this system. Given Hofmann's 0.5 M phosphate because "
            "it is the only buffered composition the corpus prints. If the "
            "real buffer is stronger the modelled pH drift is too large here, "
            "and if it is weaker the drift is too small."))
BUFFER_MCILVAINE = BufferSpec(
    kind="citrate_phosphate", phosphate_mol_l=0.2, citrate_mol_l=0.1,
    declared=True,
    source=("Kumazawa 2003 Methods, kumazawa2003_extraction.md sec. 2.3: "
            "'mixed 1/10 M citric acid and 1/5 M Na2HPO4' = 0.1 M citric acid "
            "+ 0.2 M disodium hydrogen phosphate (McIlvaine)"))
BUFFER_NONE = BufferSpec(
    kind="none", declared=True,
    source=("the source states DEIONISED WATER and no buffer "
            "(zhou2023_extraction.md sec. 1.2; kang2026_SI_extraction.md sec. 4)"))


# ===========================================================================
# THE SYSTEMS
# ===========================================================================
# Each is one pot: initial concentrations in mmol/L, a temperature, a hold time,
# a pH (and, for the unbuffered Zhou systems, a MEASURED final pH).

SYSTEMS: Dict[str, Dict[str, Any]] = {
    # ---- Hofmann & Schieberle 1998, the SIDA absolute anchors, pH 5.0 only --
    # cysteine 3.3 mmol + carbohydrate 10.0 mmol in 100 mL = 33 / 100 mmol/L,
    # 0.5 mol/L phosphate (BUFFERED, so pH is clamped), pH 5.0, 145 C / 20 min.
    "hofmann_pentose_pH5": dict(
        initial={"PENT": 100.0, "Cys": 33.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann & Schieberle 1998 Tables 1-2, ribose/xylose pH 5.0 row",
    ),
    "hofmann_glucose_pH5": dict(
        initial={"Glc": 100.0, "Cys": 33.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann & Schieberle 1998 Tables 1-2, glucose pH 5.0 row",
    ),
    "hofmann_fructose_pH5": dict(
        initial={"Fru": 100.0, "Cys": 33.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann & Schieberle 1998 Tables 1-2, fructose pH 5.0 row",
    ),
    # ---- Hofmann 1998 step-level fed-precursor systems ---------------------
    # 50 mL, 1 mmol precursor + 1 mmol co-reactant = 20 mmol/L each.
    "fed_ribose_h2s": dict(
        initial={"PENT": 20.0, "H2S": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T3/T6, ribose + H2S"),
    "fed_tdp_h2s": dict(
        initial={"TDP": 20.0, "H2S": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T3, 3-deoxyribosulose + H2S"),
    "fed_furfural_h2s": dict(
        initial={"FUR": 20.0, "H2S": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T3, furan-2-aldehyde + H2S"),
    "fed_nf_h2s": dict(
        initial={"NF": 20.0, "H2S": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T4, norfuraneol + H2S"),
    "fed_nf_cys": dict(
        initial={"NF": 20.0, "Cys": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T4, norfuraneol + cysteine"),
    "fed_c2c3": dict(
        initial={"HA": 20.0, "MP": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T10, hydroxyacetaldehyde + mercapto-2-propanone"),
    "fed_c2c3_pH3": dict(
        initial={"HA": 20.0, "MP": 20.0}, t_c=145.0, minutes=20.0, ph=3.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T10 pH ladder, pH 3.0"),
    "fed_c2c3_pH7": dict(
        initial={"HA": 20.0, "MP": 20.0}, t_c=145.0, minutes=20.0, ph=7.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T10 pH ladder, pH 7.0"),
    "fed_thiamine": dict(
        initial={"THI": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T8, thiamin"),
    "fed_mgo_h2s_1to1": dict(
        initial={"MGO": 20.0, "H2S": 20.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T7, 2-oxopropanal + H2S 1:1"),
    "fed_mgo_h2s_1to2": dict(
        initial={"MGO": 20.0, "H2S": 40.0}, t_c=145.0, minutes=20.0, ph=5.0, buffer=BUFFER_HOFMANN,
        anchor="Hofmann 1998 T7, 2-oxopropanal + H2S 1:2"),
    # ---- Zhou 2023, the pH 7 column ONLY -----------------------------------
    # ARP 20 mmol/L +/- cysteine 20 mmol/L, DEIONISED WATER (UNBUFFERED),
    # 120 C / 60 min. The pH LABEL is the INITIAL pH; the MEASURED final pH is
    # 4.13 (ARP alone) and 3.42 (ARP + Cys). Both endpoints are measured and
    # both are used -- see sulfur.integrate_sulfur's docstring and inventory
    # sec. C.11.
    "zhou_arp_cys_pH7": dict(
        initial={"ARP": 20.0, "Cys": 20.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=120.0, minutes=60.0, ph=7.0, buffer=BUFFER_NONE,
        anchor="Zhou 2023 Table 1, ARP + Cys, pH 7 column"),
    "zhou_arp_pH7": dict(
        initial={"ARP": 20.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=120.0, minutes=60.0, ph=7.0, buffer=BUFFER_NONE,
        anchor="Zhou 2023 Table 1, ARP alone, pH 7 column"),
    # ---- Whitfield 1999, fed norfuraneol at pH 4.5, 140 C / 60 min ---------
    "whitfield_nf_cys": dict(
        initial={"NF": 20.0, "Cys": 20.0}, t_c=140.0, minutes=60.0, ph=4.5, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="Whitfield & Mottram 1999 Table 1, NF + cysteine, pH 4.5"),
    "whitfield_nf_h2s": dict(
        initial={"NF": 20.0, "H2S": 40.0}, t_c=140.0, minutes=60.0, ph=4.5, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="Whitfield & Mottram 1999 Table 1, NF + H2S 1:2, pH 4.5"),
    # ---- van Seeventer 2001, the precursor-conversion pair ------------------
    "vanseeventer_130C": dict(
        initial={"PENT": 100.0, "Cys": 33.0}, t_c=130.0, minutes=20.0, ph=5.0, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="van Seeventer 2001, Hofmann process flavour, 130 C / 20 min"),
    # ---- Zhang 2024 Figure 1, the redox series, 115 C / 60 min, pH 4.9 -----
    # VB1 15 mg/mL = 44.5 mM thiamine, xylose at equal weight = 99.9 mM.
    # Cys 123.8 mM (REDUCED) vs cystine 62.4 mM (OXIDISED) at near-matched
    # molar sulfur -- the axis the whole figure is about.
    "zhang_fig1_cys": dict(
        initial={"THI": 44.5, "PENT": 99.9, "Cys": 123.8, "OX": OX_AMBIENT_MMOL_L},
        t_c=115.0, minutes=60.0, ph=4.9, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="Zhang 2024 Fig. 1, VB1-Xyl-Cys"),
    "zhang_fig1_gcys": dict(
        # cystine is BOTH a sulfur source (2 S per molecule, released as it is
        # reduced) and an OXIDANT (one disulfide = one oxidant equivalent).
        initial={"THI": 44.5, "PENT": 99.9, "Cys": 124.8, "OX": 62.4},
        t_c=115.0, minutes=60.0, ph=4.9, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="Zhang 2024 Fig. 1, VB1-Xyl-GCys (cystine)"),
    # ---- Cerny 2007b, the full ternary, 145 C / 20 min, pH 5.00 ------------
    "cerny_ternary": dict(
        initial={"Cys": 99.9, "THI": 99.9, "PENT": 299.7}, t_c=145.0,
        minutes=20.0, ph=5.0, buffer=BUFFER_PHOSPHATE_ASSUMED,
        anchor="Cerny & Briffod 2007b system A, cys:thiamine:xylose 1:1:3"),
    # =====================================================================
    # B2.1 NEW SYSTEMS
    # =====================================================================
    # ---- Kang 2026 SI Table S4: the sulfur branch's FIRST temperature -----
    #      ladder. PURIFIED TTCA 10 mmol/L in water, sealed, initial pH 7.0,
    #      120 min, at 100 / 120 / 140 C. FIT_HOLDOUT_DECLARATION.md
    #      Amendment 5 puts 100 and 120 C in the FIT column and 140 C in the
    #      GATING HOLD-OUT column, so ONLY the two lower rungs appear here.
    #      The pH label is INITIAL: the main paper measures the pH-7 run
    #      finishing at 4.9, so the trajectory state carries it, exactly as it
    #      does for Zhou.
    "kang_ttca_100": dict(
        initial={"TTCA": 10.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=100.0, minutes=120.0, ph=7.0, buffer=BUFFER_NONE, grid_points=25,
        anchor="Kang 2026 SI Table S4, 100 C column (TTCA 10 mM, pH 7, 120 min)"),
    "kang_ttca_120": dict(
        initial={"TTCA": 10.0, "OX": OX_AMBIENT_MMOL_L},
        t_c=120.0, minutes=120.0, ph=7.0, buffer=BUFFER_NONE, grid_points=25,
        anchor="Kang 2026 SI Table S4, 120 C column"),
    # ---- Kang 2026 SI Fig. S4: the FREE-CYSTEINE control system ------------
    #      A SEPARATE POT from the TTCA ladder: free cysteine 10 mmol/L, pH 7,
    #      sealed, 0-120 min, at 100 / 120 / 140 C. The dossier reads it as a
    #      free-Cys control at 85% confidence and the main text uses it that
    #      way ("thermal degradation of TTCA AND free Cys"); the alternative
    #      reading -- that it tracks the TTCA-bound moiety -- cannot be excluded
    #      from the published record and is flagged. Only 100 and 120 C are
    #      FIT; 140 C belongs to the hold-out column.
    "kang_freecys_100": dict(
        initial={"Cys": 10.0}, t_c=100.0, minutes=120.0, ph=7.0,
        buffer=BUFFER_NONE,
        anchor="Kang 2026 SI Fig. S4, free cysteine 10 mM, pH 7, 100 C"),
    "kang_freecys_120": dict(
        initial={"Cys": 10.0}, t_c=120.0, minutes=120.0, ph=7.0,
        buffer=BUFFER_NONE,
        anchor="Kang 2026 SI Fig. S4, free cysteine 10 mM, pH 7, 120 C"),
    # ---- Kumazawa 2003 Figure 3: 2-furfurylthiol SURVIVAL across seven pH --
    #      levels at 121 C / 10 min. 1 ppm FFT = 8.76e-3 mmol/L in mixed
    #      citrate / phosphate buffer (so the pH is CLAMPED), canned WITHOUT
    #      deoxidisation, i.e. air present, so the oxidant pool is the ambient
    #      charge. Six of the seven levels are used; pH 7.0 is left-censored at
    #      <0.5% residual and is excluded rather than fitted as a level.
    #      Declared FIT by FIT_HOLDOUT_DECLARATION.md Amendment 4.
    **{
        f"kumazawa_pH{ph_label}": dict(
            initial={"FFT": 8.7597e-3, "OX": OX_AMBIENT_MMOL_L},
            t_c=121.0, minutes=10.0, ph=ph_value,
            buffer=BUFFER_MCILVAINE,
            anchor=f"Kumazawa 2003 Fig. 3, pH {ph_value}, 121 C / 10 min")
        for ph_label, ph_value in (
            ("3_0", 3.0), ("4_0", 4.0), ("5_0", 5.0), ("5_4", 5.4),
            ("6_0", 6.0), ("6_4", 6.4),
        )
    },
    # ---- B2.2: Amendment 7's THREE DECLARED pH-DRIFT ANCHORS ---------------
    #      FIT_HOLDOUT_DECLARATION.md Amendment 7 (2026-08-29): "Zhou 2023's
    #      three measured final-pH endpoints (initial 6/7/8 -> final
    #      3.22/3.42/5.07, unbuffered, 120 C) are declared FIT as calibration
    #      anchors for the drift model... They are process-state observables,
    #      not volatile levels -- no volatile hold-out is affected."
    #
    #      ⚠️ THESE SYSTEMS SIT ON THE HELD-OUT pH-6 AND pH-8 COLUMNS. That is
    #      licensed for the ENDPOINT ONLY. The firewall test asserts, per row,
    #      that EVERY fit row referencing a `zhou_drift_*` system has
    #      kind == "ph_endpoint" -- so no volatile level from those columns can
    #      enter the objective even by accident.
    **{
        f"zhou_drift_pH{label}": dict(
            initial={"ARP": 20.0, "Cys": 20.0, "OX": OX_AMBIENT_MMOL_L},
            t_c=120.0, minutes=60.0, ph=ph_value, buffer=BUFFER_NONE,
            drift_only=True,
            anchor=f"Zhou 2023 Fig. 2, ARP + Cys, initial pH {ph_value}, "
                   f"120 C / 60 min, unbuffered -- FINAL pH ONLY")
        for label, ph_value in (("6", 6.0), ("7", 7.0), ("8", 8.0))
    },
}


# ===========================================================================
# THE FIT ROWS
# ===========================================================================
# kind:
#   "conc"        target is a concentration, mmol/L, of `species` at t_end
#   "molpct"      target is mol% of `basis_mmol_L` charged
#   "ratio"       target is species_a / species_b at t_end (mass = molar here)
#   "conversion"  target is 1 - final/initial of `species`
#   "share"       target is a branch share from branch_shares()
#   "ceiling"     one-sided: a branch share must not EXCEED the target


def _ppb(species: str, ug_per_l: float) -> float:
    return ug_per_litre_to_mmol_per_litre(species, ug_per_l)


FIT_ROWS: Tuple[Dict[str, Any], ...] = (
    # --- Hofmann 1998 SIDA absolute anchors, pH 5.0 -------------------------
    dict(id="hofmann_ribose_FFT", system="hofmann_pentose_pH5", kind="conc",
         benchmark_id="hofmann1998_ribose_cysteine_145C_20min_pH5", benchmark_compound="2-Furfurylthiol (FFT)",
         species="FFT", target=_ppb("FFT", 121.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, ribose pH 5.0: 12.1 ug/100 mL = 121 ppb",
         note="the sulfur branch's primary absolute anchor"),
    dict(id="hofmann_ribose_MFT", system="hofmann_pentose_pH5", kind="conc",
         benchmark_id="hofmann1998_ribose_cysteine_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", target=_ppb("MFT", 198.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, ribose pH 5.0: 19.8 ug/100 mL = 198 ppb"),
    dict(id="hofmann_xylose_FFT", system="hofmann_pentose_pH5", kind="conc",
         benchmark_id="hofmann1998_xylose_cysteine_145C_20min_pH5", benchmark_compound="2-Furfurylthiol (FFT)",
         species="FFT", target=_ppb("FFT", 96.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, xylose pH 5.0: 9.6 ug/100 mL = 96 ppb",
         note="SAME prediction as the ribose row: the module carries ONE "
              "generic aldopentose. The 1.26x gap is a declared limitation and "
              "shows up here as an irreducible residual."),
    dict(id="hofmann_xylose_MFT", system="hofmann_pentose_pH5", kind="conc",
         benchmark_id="hofmann1998_xylose_cysteine_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", target=_ppb("MFT", 143.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, xylose pH 5.0: 14.3 ug/100 mL = 143 ppb"),
    dict(id="hofmann_glucose_FFT", system="hofmann_glucose_pH5", kind="conc",
         benchmark_id="hofmann1998_glucose_cysteine_145C_20min_pH5", benchmark_compound="2-Furfurylthiol (FFT)",
         species="FFT", target=_ppb("FFT", 28.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, glucose pH 5.0: 2.8 ug/100 mL = 28 ppb"),
    dict(id="hofmann_glucose_MFT", system="hofmann_glucose_pH5", kind="conc",
         benchmark_id="hofmann1998_glucose_cysteine_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", target=_ppb("MFT", 19.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, glucose pH 5.0: 1.9 ug/100 mL = 19 ppb"),
    dict(id="hofmann_fructose_FFT", system="hofmann_fructose_pH5", kind="conc",
         benchmark_id="hofmann1998_fructose_cysteine_145C_20min_pH5", benchmark_compound="2-Furfurylthiol (FFT)",
         species="FFT", target=_ppb("FFT", 32.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, fructose pH 5.0: 3.2 ug/100 mL = 32 ppb"),
    dict(id="hofmann_fructose_MFT", system="hofmann_fructose_pH5", kind="conc",
         benchmark_id="hofmann1998_fructose_cysteine_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", target=_ppb("MFT", 25.0), sigma_log=0.35,
         anchor="Hofmann 1998 T1, fructose pH 5.0: 2.5 ug/100 mL = 25 ppb"),
    # --- Hofmann 1998 T5, the in-situ levels in the same pot ----------------
    dict(id="hofmann_ribose_NF_insitu", system="hofmann_pentose_pH5", kind="conc",
         species="NF", target=54530.0 / 114.10 * 1e-2, sigma_log=0.35,
         anchor="Hofmann 1998 T5, ribose pH 5.0: 54 530 ug/100 mL = 4.78 mmol/L",
         note="ribose makes ~2750x more norfuraneol than MFT in one pot, which "
              "is the authors' own reason for doubting NF as THE key MFT "
              "intermediate. Getting this level right while keeping MFT small "
              "is what forces NF to be a NODE and not a funnel."),
    dict(id="hofmann_ribose_FUR_insitu", system="hofmann_pentose_pH5", kind="conc",
         species="FUR", target=67.5 / 96.08 * 1e-2, sigma_log=0.5,
         anchor="Hofmann 1998 T5, ribose pH 5.0: 67.5 ug/100 mL furan-2-aldehyde"),
    # --- Hofmann 1998 step-level fed-precursor mol% -------------------------
    dict(id="fed_ribose_h2s_FFT", system="fed_ribose_h2s", kind="molpct",
         species="FFT", basis=20.0, target=0.008, sigma_log=0.45,
         anchor="Hofmann 1998 T3: ribose + H2S -> FFT, 9.2 ug, 0.008 mol%"),
    dict(id="fed_ribose_h2s_MFT", system="fed_ribose_h2s", kind="molpct",
         species="MFT", basis=20.0, target=0.01, sigma_log=0.45,
         anchor="Hofmann 1998 T6: ribose + H2S -> MFT, 15.1 ug, 0.01 mol%"),
    dict(id="fed_tdp_h2s_FFT", system="fed_tdp_h2s", kind="molpct",
         species="FFT", basis=20.0, target=0.08, sigma_log=0.45,
         anchor="Hofmann 1998 T3: 3-deoxyribosulose + H2S -> FFT, 78.6 ug, 0.08 mol%"),
    dict(id="fed_furfural_h2s_FFT", system="fed_furfural_h2s", kind="molpct",
         benchmark_id="hofmann1998_furan2aldehyde_h2s_145C_20min_pH5", benchmark_compound="2-Furfurylthiol (FFT)",
         species="FFT", basis=20.0, target=0.48, sigma_log=0.45,
         anchor="Hofmann 1998 T3: furan-2-aldehyde + H2S -> FFT, 550.8 ug, 0.48 mol%",
         note="furfural is 60x ribose and 7x the 3-deoxyosone as an FFT source"),
    dict(id="fed_nf_h2s_MFT", system="fed_nf_h2s", kind="molpct",
         benchmark_id="hofmann1998_norfuraneol_h2s_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", basis=20.0, target=0.19, sigma_log=0.45,
         anchor="Hofmann 1998 T4: norfuraneol + H2S -> MFT, 211.2 ug, 0.19 mol%"),
    dict(id="fed_nf_cys_MFT", system="fed_nf_cys", kind="molpct",
         benchmark_id="hofmann1998_norfuraneol_cysteine_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", basis=20.0, target=0.05, sigma_log=0.45,
         anchor="Hofmann 1998 T4: norfuraneol + cysteine -> MFT, 50.8 ug, 0.05 mol%"),
    dict(id="fed_c2c3_MFT", system="fed_c2c3", kind="molpct",
         benchmark_id="hofmann1998_c2c3_recombination_145C_20min_pH5", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", basis=20.0, target=0.24, sigma_log=0.45,
         anchor="Hofmann 1998 T10: C2 + C3 -> MFT, 268.1 ug, 0.24 mol%",
         note="the single most effective MFT route measured anywhere in the corpus"),
    dict(id="fed_c2c3_MFT_pH3", system="fed_c2c3_pH3", kind="molpct",
         benchmark_id="hofmann1998_c2c3_recombination_145C_20min_pH3", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", basis=20.0, target=0.01, sigma_log=0.6,
         anchor="Hofmann 1998 T10 pH ladder: 15.5 ug at pH 3.0"),
    dict(id="fed_c2c3_MFT_pH7", system="fed_c2c3_pH7", kind="molpct",
         benchmark_id="hofmann1998_c2c3_recombination_145C_20min_pH7", benchmark_compound="2-Methyl-3-furanthiol (MFT)",
         species="MFT", basis=20.0, target=0.27, sigma_log=0.6,
         anchor="Hofmann 1998 T10 pH ladder: 311.5 ug at pH 7.0",
         note="THE C2+C3 ROUTE HAS NO pH FACTOR IN THIS MODULE. Assigning one "
              "would be inventing a mechanism to fit two numbers, so these two "
              "rows are carried in the objective and the residual is REPORTED. "
              "They are the module's honest pH-shape cost on this lane."),
    dict(id="fed_thiamine_MFT", system="fed_thiamine", kind="molpct",
         species="MFT", basis=20.0, target=0.01, sigma_log=0.45,
         anchor="Hofmann 1998 T8: thiamin -> MFT, 8.2 ug, 0.01 mol%",
         note="thiamin is 30x weaker than the C2+C3 route in a FED experiment; "
              "Cerny's ternary makes it the MAJORITY route in situ. Both are "
              "true and the difference is concentration -- which is exactly why "
              "no fixed branch fraction can be right."),
    dict(id="fed_mgo_h2s_1to1_MP", system="fed_mgo_h2s_1to1", kind="molpct",
         species="MP", basis=20.0, target=1.8, sigma_log=0.4,
         anchor="Hofmann 1998 T7: 2-oxopropanal + H2S 1:1 -> 1650 ug, 1.8 mol%"),
    dict(id="fed_mgo_h2s_1to2_MP", system="fed_mgo_h2s_1to2", kind="molpct",
         species="MP", basis=20.0, target=4.0, sigma_log=0.4,
         anchor="Hofmann 1998 T7: 2-oxopropanal + H2S 1:2 -> 3600 ug, 4.0 mol%",
         note="doubling H2S gives 2.2x product. A first-order-in-H2S mass "
              "action law gives less than 2x once the sulfide is consumed, so "
              "this row is where the SUPER-LINEARITY shows up as a residual "
              "rather than being absorbed into an invented reaction order."),
    # --- Zhou 2023, pH 7 column only ---------------------------------------
    dict(id="zhou_pH7_MFT", system="zhou_arp_cys_pH7", kind="conc",
         species="MFT", target=_ppb("MFT", 1588.57), sigma_log=0.5,
         anchor="Zhou 2023 T1 pH 7: MFT 1588.57 +/- 21.24 ug/L",
         note="HS-SPME with external calibration, NOT SIDA. Flagged "
              "absolute_concentration: false (inventory sec. B9.1); carried at "
              "a wide sigma for that reason."),
    dict(id="zhou_pH7_FFT", system="zhou_arp_cys_pH7", kind="conc",
         species="FFT", target=_ppb("FFT", 757.965), sigma_log=0.5,
         anchor="Zhou 2023 T1 pH 7: FFT 757.965 +/- 13.03 ug/L"),
    dict(id="zhou_pH7_MFT_over_FFT", system="zhou_arp_cys_pH7", kind="ratio",
         species="MFT", species_b="FFT", target=2.096, sigma_log=0.25,
         anchor="Zhou 2023 T1 pH 7, MFT/FFT = 2.096 [D]",
         note="RATIO-ONLY, and the strongest row in the Zhou block: it is "
              "immune to the HS-SPME absolute-calibration defect and it "
              "cross-validates against Hofmann 1998 T1 pH 5 (1.636) across two "
              "labs, two methods and two feedstocks."),
    dict(id="zhou_pH7_dimer_over_MFT", system="zhou_arp_cys_pH7", kind="ratio",
         species="MFTD", species_b="MFT", target=0.0323, sigma_log=0.4,
         anchor="Zhou 2023 T1 pH 7: 102.59 ug/L dimer against 1588.57 MFT "
                "(molar ratio 0.0323; the inventory's 6.5% is in THIOL "
                "EQUIVALENTS, i.e. 2x this, because one dimer carries two "
                "thiols)"),
    dict(id="zhou_pH7_ACTZ", system="zhou_arp_cys_pH7", kind="conc",
         species="ACTZ", target=_ppb("ACTZ", 11.70), sigma_log=0.6,
         anchor="Zhou 2023 T1 pH 7: 2-acetylthiazole 11.70 +/- 2.14 ug/L"),
    dict(id="zhou_pH7_FUR_arp_alone", system="zhou_arp_pH7", kind="conc",
         species="FUR", target=_ppb("FUR", 1339.37), sigma_log=0.5,
         anchor="Zhou 2023 T1 pH 7, ARP ALONE: 2-furfural 1339.37 +/- 83.04 ug/L"),
    # --- Whitfield 1999, fed NF at pH 4.5 -----------------------------------
    dict(id="whitfield_nf_cys_MFT", system="whitfield_nf_cys", kind="molpct",
         species="MFT", basis=20.0, target=0.150, sigma_log=0.5,
         anchor="Whitfield & Mottram 1999 T1: NF + cysteine, pH 4.5, 140 C, "
                "0.150 mol% (DHS, response factors ASSUMED 1)"),
    dict(id="whitfield_nf_h2s_MFT", system="whitfield_nf_h2s", kind="molpct",
         species="MFT", basis=20.0, target=0.120, sigma_log=0.5,
         anchor="Whitfield & Mottram 1999 T1: NF + H2S 1:2, pH 4.5, 0.120 mol%",
         note="two labs, two methods, agreement within ~1.6x on the H2S channel "
              "(against Hofmann's 0.19 mol% at 145 C)"),
    dict(id="whitfield_mercaptoketone_over_MFT", system="whitfield_nf_cys",
         kind="ratio", species="MP3P", species_b="MFT", target=16.3, sigma_log=0.4,
         anchor="Whitfield 1999: mercaptoketones : MFT = 16.3 : 1 from fed NF",
         note="norfuraneol's DOMINANT fate is not MFT. A ratio, so immune to "
              "the response-factor caveat."),
    # --- van Seeventer 2001, the reactant side ------------------------------
    dict(id="vanseeventer_cys_conversion", system="vanseeventer_130C",
         kind="conversion", species="Cys", target=0.55, sigma_log=0.30,
         anchor="van Seeventer 2001: cysteine 33 -> 15 mM (55%), 130 C / 20 min",
         note="a REACTANT-side constraint: it cannot be traded against any "
              "product-side residual, so it adds information without competing."),
    dict(id="vanseeventer_ribose_conversion", system="vanseeventer_130C",
         kind="conversion", species="PENT", target=0.75, sigma_log=0.30,
         anchor="van Seeventer 2001: ribose 100 -> 25 mM (75%), 130 C / 20 min"),
    # --- Zhang 2024 Figure 1, the redox axis, as RATIOS ---------------------
    dict(id="zhang_fig1_cys_dimer_over_MFT", system="zhang_fig1_cys", kind="ratio",
         species="MFTD", species_b="MFT", target=0.0429, sigma_log=0.4,
         anchor="Zhang 2024 Fig. 1: 0.115 ng/mL dimer against 1.34 MFT = 8.6% "
                "BY MASS = 0.0429 molar (the dimer's MW is 2x the monomer's)"),
    dict(id="zhang_fig1_gcys_dimer_over_MFT", system="zhang_fig1_gcys", kind="ratio",
         species="MFTD", species_b="MFT", target=0.2711, sigma_log=0.4,
         anchor="Zhang 2024 Fig. 1: 1.09 against 2.01 = 54.2% by mass = 0.2711 molar",
         note="THE REDOX ROW. At near-matched molar sulfur (123.8 vs 124.8 mM "
              "S) the oxidised additive sends 6.3x more of the MFT pool to the "
              "dimer. This is the only thing in the corpus that identifies the "
              "oxidant term, which is why it is a FIT row and Fig. 2 is not."),
    dict(id="zhang_fig1_cys_MMFT_over_MFT", system="zhang_fig1_cys", kind="ratio",
         species="MMFT", species_b="MFT", target=0.0213, sigma_log=0.5,
         anchor="Zhang 2024 Fig. 1: MMFT 0.04 against MFT 1.34 ng/mL = 0.0213 molar"),
    # --- Cerny 2007b, the full-ternary branch split -------------------------
    dict(id="cerny_ternary_thiamine_share", system="cerny_ternary", kind="share",
         share_key="MFT_share_thiamine_route", target=0.54, sigma_log=0.20,
         anchor="Cerny 2007b Table 3, full ternary: 54% unlabelled (thiamine) : "
                "46% 13C5-labelled (xylose)",
         note="INGEST THE TABLE, NOT THE PROSE: both Cerny papers' running text "
              "says '54% labelled' where the tables say 46% labelled / 54% "
              "unlabelled (inventory sec. F row 15). The prose error repeats in "
              "BOTH papers."),
    dict(id="cerny_FFT_xylose_share", system="cerny_ternary", kind="share",
         share_key="MFT_share_thiamine_route", target=None, sigma_log=None,
         skip=True,
         anchor="placeholder -- see cerny_fft_from_sugar below"),
    dict(id="cerny_MP3P_sugar_share", system="cerny_ternary", kind="ratio",
         species="MP3P", species_b="MP2P", target=None, sigma_log=None, skip=True,
         anchor="placeholder -- superseded by the isomer-split rows below"),
    dict(id="cerny_isomer_split", system="cerny_ternary", kind="ratio",
         species="MP3P", species_b="MP2P", target=1.0, sigma_log=0.6,
         anchor="Cerny 2007 T4: 2-mercapto-3-pentanone is 94->95% XYLOSE-derived "
                "while 3-mercapto-2-pentanone is 77-90% THIAMINE-derived; "
                "Whitfield 1999 T1 measures the two isomers at 74.5 : 77.5, "
                "i.e. ~1:1 in amount",
         note="THE ISOMER-SPLIT DIAGNOSTIC. It constrains the two routes to "
              "deliver comparable amounts of their OWN diagnostic isomer, which "
              "is a much sharper statement than either yield alone. Cerny 2007 "
              "Table 4's MFT COLUMN is deliberately EXCLUDED -- see the module "
              "docstring's disjointness note."),
    # =======================================================================
    # B2.1 NEW FIT ROWS
    # =======================================================================
    # --- Kang 2026 SI Table S4, the 100 and 120 C rungs ONLY ---------------
    # Tier A, external calibration curves printed (MFT R^2 0.9989, FFT 0.9992,
    # furfural 0.9923), raster-verified mu-g/L, arithmetic closure to +/-0.003
    # on all nine subtotals, and the 120 C column independently anchored to the
    # main paper's own Table 1 control row on THREE quantities to three
    # decimals. The SDs are NOT used: Table S5's pH-7 column reproduces Table
    # S4's 100 C column mean-for-mean with DIFFERENT SDs, so one of the two is
    # mis-pasted and neither set is trustworthy. sigma_log below is the
    # dossier's replacement uncertainty (sec. 7d: +/-15% relative for Tier A)
    # widened for the fact that this module carries ONE generic aldopentose and
    # a lumped TTCA topology.
    dict(id="kang_100C_MFT", system="kang_ttca_100", kind="conc",
         species="MFT", target=_ppb("MFT", 1.237), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 100 C: MFT 1.237 ug/L (Tier A, R^2 0.9989)",
         note="THE FIRST TEMPERATURE LADDER THE SULFUR BRANCH HAS EVER HAD. "
              "B2 had no temperature axis inside any single system, which is "
              "why its one lumped formation Ea landed on its bound and was "
              "reported as a failure to determine."),
    dict(id="kang_120C_MFT", system="kang_ttca_120", kind="conc",
         species="MFT", target=_ppb("MFT", 1.388), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 120 C: MFT 1.388 ug/L; the same value "
                "appears in the main paper's Table 1 control row, which "
                "independently anchors this column"),
    dict(id="kang_100C_FFT", system="kang_ttca_100", kind="conc",
         species="FFT", target=_ppb("FFT", 3.734), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 100 C: FFT 3.734 ug/L (Tier A, R^2 0.9992)"),
    dict(id="kang_120C_FFT", system="kang_ttca_120", kind="conc",
         species="FFT", target=_ppb("FFT", 4.107), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 120 C: FFT 4.107 ug/L"),
    dict(id="kang_100C_FUR", system="kang_ttca_100", kind="conc",
         species="FUR", target=_ppb("FUR", 3.381), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 100 C: furfural 3.381 ug/L",
         note="the best-behaved Arrhenius series in the table, and the one the "
              "main paper's own '3.3 times that at 100 C' sentence "
              "independently confirms"),
    dict(id="kang_120C_FUR", system="kang_ttca_120", kind="conc",
         species="FUR", target=_ppb("FUR", 5.793), sigma_log=0.4,
         anchor="Kang 2026 SI T-S4, 120 C: furfural 5.793 ug/L"),
    dict(id="kang_100C_MFT_over_FFT", system="kang_ttca_100", kind="ratio",
         species="MFT", species_b="FFT", target=1.0 / 3.02, sigma_log=0.3,
         anchor="Kang 2026 SI sec. 7b: FFT/MFT = 3.02 at 100 C",
         note="RATIO, immune to the calibration axis. The dossier characterises "
              "this ratio on three axes and finds it INVARIANT to pH (2.73-3.02 "
              "over 2.5 units) and only weakly sensitive to temperature -- but "
              "spanning 62x across nitrogen co-substrates. It is a shape "
              "constraint, not a level."),
    # The cysteine ladder: MEASURED conversions at the two fitted temperatures,
    # carrying Kang's measured Ea(free-Cys depletion) = 55.1 kJ/mol implicitly.
    # They are entered as CONVERSIONS rather than as an Ea because this module
    # computes cysteine depletion from its own competing channels; adding a
    # separate lumped Cys sink with Kang's own k would DOUBLE-COUNT the
    # depletion the network already predicts.
    dict(id="kang_100C_cys_conversion", system="kang_freecys_100",
         kind="conversion", species="Cys",
         target=0.162, sigma_log=0.3,
         anchor="Kang 2026 SI Fig. S4 (digitised): free-cysteine conversion "
                "16.2% at 100 C / 120 min; the three-point Arrhenius gives "
                "Ea = 55.1 kJ/mol, R^2 0.994, robust to the fitting window "
                "at +/-10%"),
    dict(id="kang_120C_cys_conversion", system="kang_freecys_120",
         kind="conversion", species="Cys",
         target=0.387, sigma_log=0.3,
         anchor="Kang 2026 SI Fig. S4 (digitised): 38.7% at 120 C / 120 min",
         note="THE 140 C RUNG (62.6%) IS NOT HERE. It belongs to the gating "
              "hold-out column."),
    dict(id="kang_free_cys_yield_ceiling", system="kang_ttca_120",
         kind="peak_fraction_ceiling", species="Cys", basis=10.0,
         target=KANG_TTCA_FREE_CYS_YIELD_CEILING_MOL_PCT / 100.0,
         sigma_log=0.25,
         anchor=KANG_TTCA_CEILING_ANCHOR,
         note="ONE-SIDED MASS-BALANCE CEILING, and it is the kind of constraint "
              "a network violates silently. Being below it costs nothing."),
    # --- Kumazawa 2003 Fig. 3, the seven-point pH grid at 121 C ------------
    # RESIDUAL RATIOS (heated / unheated half of the SAME solution), so they are
    # ratio-scale and immune to the absolute-calibration defects that qualify
    # most level rows in this corpus. THESE ARE THE ROWS THAT PAY FOR THE
    # THIOLATE CHANNEL: they are the only measurement in the corpus of a thiol's
    # pH response with formation held out of the picture entirely.
    *[
        dict(id=f"kumazawa_FFT_retention_pH{label}",
             system=f"kumazawa_pH{label}", kind="retention",
             species="FFT", basis=8.7597e-3, target=value, sigma_log=0.20,
             anchor=f"Kumazawa 2003 Fig. 3: {value * 100:.1f}% of 1 ppm "
                    f"2-furfurylthiol survives 121 C / 10 min at pH {ph_v}")
        for label, ph_v, value in (
            ("3_0", 3.0, 0.995), ("4_0", 4.0, 0.962), ("5_0", 5.0, 0.891),
            ("5_4", 5.4, 0.795), ("6_0", 6.0, 0.451), ("6_4", 6.4, 0.11),
        )
    ],
    # =======================================================================
    # B2.2 NEW FIT ROWS -- Amendment 7's pH-drift anchors, and ONLY these
    # =======================================================================
    # THREE rows, TWO calibrated constants. They are scored in pH UNITS, not in
    # log-fold like every other row, because a pH IS a logarithm already;
    # sigma_ph is 0.25 units -- four times the dossier's +/-0.06 digitisation
    # error on Fig. 2, widened for the fact that the model carries one generic
    # aldopentose Amadori and a lumped acid pool.
    #
    # THE PREDICTED QUANTITY IS THE **COOLED** pH, because that is what a pH
    # meter reads. See ph_state's "TWO pH SCALES" block.
    *[
        dict(id=f"zhou_final_pH_from_pH{label}", system=f"zhou_drift_pH{label}",
             kind="ph_endpoint", target=value, sigma_log=0.25,
             anchor=f"Zhou 2023 Fig. 2 [D, +/-0.06]: ARP + Cys, initial pH "
                    f"{label}.0 -> final pH {value} (120 C / 60 min, "
                    f"unbuffered). DECLARED FIT by Amendment 7.",
             note="a PROCESS-STATE observable. No volatile level from the "
                  "pH-6 or pH-8 columns enters this objective.")
        for label, value in (("6", 3.22), ("7", 3.42), ("8", 5.07))
    ],
    # --- Yaghmur 2005, the one-sided ceiling --------------------------------
    dict(id="yaghmur_fft_share_ceiling", system="hofmann_pentose_pH5",
         kind="ceiling", share_key="FFT_share_of_furfural_flux",
         target=YAGHMUR_FFT_SHARE_CEILING, sigma_log=0.3,
         anchor=YAGHMUR_CEILING_ANCHOR,
         note="ONE-SIDED. Being below the ceiling costs nothing; exceeding it "
              "is penalised. A one-sided constraint is the only honest way to "
              "use a bound."),
)

ACTIVE_FIT_ROWS: Tuple[Dict[str, Any], ...] = tuple(
    row for row in FIT_ROWS if not row.get("skip")
)


# ===========================================================================
# THE OBJECTIVE
# ===========================================================================

PARAM_ORDER: Tuple[str, ...] = FITTED_SULFUR_KEYS
N_K = len(PARAM_ORDER)

#: The five NON-RATE parameters, in vector order after the N_K rate constants.
#: FOUR of them are new in B2.2 and each is named, bounded and reported.
EXTRA_PARAM_NAMES: Tuple[str, ...] = (
    "Ea_lumped_formation",
    "Ea_decay_thiol_sink",
    "Ea_decay_carbonyl_sink",
    "ph_acid_yield_per_sink_event",
    "ph_arp_secondary_ammonium_pKa",
)
N_EXTRA = len(EXTRA_PARAM_NAMES)
DECAY_FAMILY_ORDER: Tuple[str, ...] = ("thiol_sink", "carbonyl_sink")


def unpack(x: np.ndarray) -> Tuple[Dict[str, float], float, Dict[str, float], PhDrift]:
    fitted = {key: float(x[i]) for i, key in enumerate(PARAM_ORDER)}
    formation_ea = float(x[N_K])
    decay_ea = {
        family: float(x[N_K + 1 + i])
        for i, family in enumerate(DECAY_FAMILY_ORDER)
    }
    drift = PhDrift(
        acid_yield=float(np.clip(x[N_K + 3], *ACID_YIELD_BOUNDS)),
        arp_amine_pka=float(np.clip(x[N_K + 4], *ARP_AMINE_PKA_BOUNDS)),
    )
    return fitted, formation_ea, decay_ea, drift


def build_parameters(x: np.ndarray) -> Dict[str, Any]:
    """Assemble the full operative parameter set from the fit vector."""
    fitted, formation_ea, decay_ea, _drift = unpack(x)
    parameters: Dict[str, Any] = dict(operative_parameters(B1_FITTED))
    parameters.update(MEASURED_SULFUR)
    parameters.update(with_fitted_sulfur(fitted, formation_ea, decay_ea))
    return parameters


#: Fixed-point settings for the dynamic pH. The quick setting is what the
#: optimiser runs on; the careful setting is what the frozen report is
#: regenerated with, and the difference between the two is REPORTED rather than
#: assumed away.
PH_NODES_QUICK, PH_ITERS_QUICK = 9, 2
PH_NODES_CAREFUL, PH_ITERS_CAREFUL = 41, 4


def run_systems(parameters: Mapping[str, Any], quick: bool,
                drift: PhDrift = None) -> Dict[str, Any]:
    """Integrate every system once and cache what the rows need."""
    out: Dict[str, Any] = {}
    for name, spec in SYSTEMS.items():
        minutes = float(spec["minutes"])
        # B2.1: a system carrying a TRANSIENT-INTERMEDIATE row needs its whole
        # trajectory, not just the endpoint, because the quantity measured is a
        # PEAK. Kang's free cysteine rises then falls at every temperature, and
        # the 16.3 mol% ceiling is on the maximum of that curve.
        points = int(spec.get("grid_points", 0))
        grid = (
            np.linspace(0.0, minutes, points) if points >= 3
            else np.array([0.0, minutes])
        )
        run = integrate_sulfur(
            parameters,
            float(spec["t_c"]) + CELSIUS,
            spec["initial"],
            grid,
            ph=float(spec["ph"]),
            buffer_spec=spec.get("buffer"),
            ph_drift=drift,
            ph_nodes=PH_NODES_QUICK if quick else PH_NODES_CAREFUL,
            ph_iterations=PH_ITERS_QUICK if quick else PH_ITERS_CAREFUL,
            rtol=1e-5 if quick else 1e-8,
            atol=1e-14,
        )
        out[name] = {"run": run, "spec": spec}
    return out


def _needs_flux(rows) -> set:
    return {r["system"] for r in rows if r["kind"] in ("share", "ceiling")}


def observables(
    parameters: Mapping[str, Any], quick: bool, drift: PhDrift = None
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """Predicted value of every active fit row."""
    cache = run_systems(parameters, quick, drift)
    fluxes: Dict[str, Dict[str, float]] = {}
    for name in _needs_flux(ACTIVE_FIT_ROWS):
        spec = SYSTEMS[name]
        fluxes[name] = sulfur_flux_budget(
            parameters,
            float(spec["t_c"]) + CELSIUS,
            spec["initial"],
            float(spec["minutes"]),
            ph=float(spec["ph"]),
            buffer_spec=spec.get("buffer"),
            ph_drift=drift,
            n_points=21 if quick else 101,
        )

    predicted: Dict[str, float] = {}
    for row in ACTIVE_FIT_ROWS:
        run = cache[row["system"]]["run"]
        kind = row["kind"]
        if kind == "conc":
            predicted[row["id"]] = run.final(row["species"])
        elif kind == "molpct":
            predicted[row["id"]] = 100.0 * run.final(row["species"]) / float(row["basis"])
        elif kind == "ratio":
            denominator = run.final(row["species_b"])
            predicted[row["id"]] = (
                run.final(row["species"]) / denominator if denominator > 0 else np.nan
            )
        elif kind == "cross_system_ratio":
            # ADDED BY WAVE B8, ADDITIVE ONLY: no row of B2.1/B2.2/B2.3/B2.4
            # uses this kind, so every earlier artefact stays bit-reproducible.
            # It exists because Amendment 17 re-bands the Kang/Zhai ladder to
            # SEMI-QUANT (f' = 1, n = 1), and a single-IS semi-quant source
            # licenses SHAPE and RATIO but not absolute magnitude. The
            # response-factor-free form of "MFT rises 1.12x from 100 to 120 C"
            # is a ratio of the SAME species between TWO systems -- which is
            # exactly what this kind scores and what no existing kind could.
            denominator = cache[row["system_b"]]["run"].final(row["species"])
            predicted[row["id"]] = (
                run.final(row["species"]) / denominator if denominator > 0
                else np.nan
            )
        elif kind == "conversion":
            initial = float(SYSTEMS[row["system"]]["initial"].get(row["species"], 0.0))
            predicted[row["id"]] = (
                (initial - run.final(row["species"])) / initial if initial > 0 else np.nan
            )
        elif kind in ("share", "ceiling"):
            shares = branch_shares(fluxes[row["system"]])
            predicted[row["id"]] = float(shares[row["share_key"]])
        elif kind == "retention":
            # fraction of the CHARGED species surviving to t_end
            predicted[row["id"]] = run.final(row["species"]) / float(row["basis"])
        elif kind == "ph_endpoint":
            # THE COOLED pH -- what a pH meter reads. See ph_state.
            predicted[row["id"]] = float(run.metadata["ph_final_cooled"])
        elif kind == "peak_fraction_ceiling":
            peak = float(np.max(run.series(row["species"])))
            predicted[row["id"]] = peak / float(row["basis"])
        else:
            raise ValueError(f"unknown row kind {kind!r}")
    return predicted, {"cache": cache, "fluxes": fluxes}


FLOOR = 1e-12


_PROGRESS = {"n": 0, "best": float("inf"), "t0": time.time()}


def residuals(x: np.ndarray, quick: bool = True) -> np.ndarray:
    try:
        parameters = build_parameters(x)
        _f, _e, _d, drift = unpack(x)
        predicted, _ = observables(parameters, quick, drift)
    except Exception:
        _progress(np.full(len(ACTIVE_FIT_ROWS), 25.0))
        return np.full(len(ACTIVE_FIT_ROWS), 25.0)
    out = np.empty(len(ACTIVE_FIT_ROWS))
    for i, row in enumerate(ACTIVE_FIT_ROWS):
        p = predicted[row["id"]]
        t = float(row["target"])
        sigma = float(row["sigma_log"])
        if not np.isfinite(p):
            out[i] = 25.0
            continue
        if row["kind"] == "ph_endpoint":
            # a pH IS a logarithm; the residual is in pH UNITS, not in folds.
            out[i] = (p - t) / sigma
        elif row["kind"] in ("ceiling", "peak_fraction_ceiling"):
            # one-sided: no penalty below the ceiling
            out[i] = max(0.0, math.log10((p + FLOOR) / (t + FLOOR))) / sigma
        else:
            out[i] = math.log10((p + FLOOR) / (t + FLOOR)) / sigma
    out = np.clip(out, -25.0, 25.0)
    _progress(out)
    return out


def _progress(r: np.ndarray) -> None:
    """Print a progress line every 100 evaluations, so a long fit is visible."""
    _PROGRESS["n"] += 1
    if _PROGRESS["n"] % 50 == 0:
        # The first attempt at this fit was lost to the container's memory
        # ceiling after ~2 h: repeated stiff solve_ivp calls accumulate, and at
        # ~300 000 integrations the process reached 4.3 GB of an 8 GB box.
        # Collecting explicitly is cheap next to an integration and it is the
        # difference between a fit that finishes and one that does not.
        gc.collect()
    cost = 0.5 * float(np.dot(r, r))
    _PROGRESS["best"] = min(_PROGRESS["best"], cost)
    if _PROGRESS["n"] % 100 == 0:
        print(f"    [{_PROGRESS['n']:5d} evals, "
              f"{time.time() - _PROGRESS['t0']:6.0f}s] best cost "
              f"{_PROGRESS['best']:.3f}", flush=True)


# ===========================================================================
# THE FIT
# ===========================================================================


#: Coordinate-descent sweeps before the Gauss-Newton polish. ZERO is the
#: shipped setting and the reason is measured, not aesthetic: once the three
#: STRUCTURALLY MISSING steps that the first fit attempt exposed were added
#: (`r_pent_caramel` and `r_pent_thermal_dpo`, without which a fed ribose + H2S
#: pot with no amine can react at all, and `r_glc_fur`, without which no hexose
#: can reach 2-furfurylthiol), no fit row sits further than ~4.5 log units from
#: its target at the uniform start, and the residual surface has a usable
#: gradient everywhere. The descent is retained, and can be switched on, because
#: it is what DIAGNOSED that flat surface -- it made no progress at all, and a
#: line search that cannot move a single coordinate is unambiguous evidence
#: that the model, not the optimiser, is what is stuck.
DESCENT_SWEEPS = 0


def _cost(x: np.ndarray, quick: bool) -> float:
    r = residuals(x, quick)
    return 0.5 * float(np.dot(r, r))


def coordinate_descent(
    x: np.ndarray, lower: np.ndarray, upper: np.ndarray, quick: bool,
    sweeps: int, rng,
) -> Tuple[np.ndarray, float]:
    """
    A robust greedy line search, one parameter at a time.

    WHY THIS AND NOT A GRADIENT SEARCH FROM A RANDOM POINT. The residual
    surface of a 34-parameter mass-action network over 22 pots is flat in most
    directions at a random start: a step that has not yet been switched on
    produces a prediction of ~1e-13 against a target of ~1e-3, so its
    FINITE-DIFFERENCE gradient is numerically indistinguishable from zero and a
    trust-region method terminates after a handful of evaluations without
    having moved. A coordinate line search does not need a gradient; it just
    evaluates. It is slower per unit of progress and it cannot exploit
    curvature, which is why the Gauss-Newton polish still runs afterwards -- but
    it reliably gets OFF the flat region, which the polish alone does not.

    NOTHING ABOUT THE STARTING POINT ENCODES A LITERATURE VALUE, because there
    are no literature values: not one of these 34 constants is printed anywhere
    in the corpus. The start is a single uniform log10 k, identical for every
    step, plus a mid-range lumped Ea.
    """
    x = x.copy()
    best = _cost(x, quick)
    if sweeps <= 0:
        return x, best
    for sweep in range(sweeps):
        order = rng.permutation(len(x))
        for i in order:
            span = (upper[i] - lower[i])
            steps = [0.75, 0.25, -0.25, -0.75] if sweep == 0 else [0.3, 0.1, -0.1, -0.3]
            scale = span / 14.0 if i < N_K else span / 6.0
            improved = True
            while improved:
                improved = False
                for step in steps:
                    trial = x.copy()
                    trial[i] = float(np.clip(x[i] + step * scale * 4.0,
                                             lower[i], upper[i]))
                    if trial[i] == x[i]:
                        continue
                    c = _cost(trial, quick)
                    if c < best - 1e-9:
                        x, best, improved = trial, c, True
                        break
        print(f"  sweep {sweep}: cost {best:.4f}", flush=True)
    return x, best


def fit(starts: int, quick: bool, max_nfev: int, seed: int = 20260828) -> Dict[str, Any]:
    from scipy.optimize import least_squares

    lower = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][0] for k in PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[0], DECAY_EA_BOUNDS[0], DECAY_EA_BOUNDS[0],
           ACID_YIELD_BOUNDS[0], ARP_AMINE_PKA_BOUNDS[0]]
    )
    upper = np.array(
        [FITTED_SULFUR_BOUNDS_LOG10K[k][1] for k in PARAM_ORDER]
        + [LUMPED_FORMATION_EA_BOUNDS[1], DECAY_EA_BOUNDS[1], DECAY_EA_BOUNDS[1],
           ACID_YIELD_BOUNDS[1], ARP_AMINE_PKA_BOUNDS[1]]
    )
    rng = np.random.default_rng(seed)
    best = None
    trace: List[Dict[str, Any]] = []
    for start in range(starts):
        # A single uniform log10 k for every step, jittered per start. See
        # coordinate_descent's docstring for why this is not an informed start:
        # there is no literature value for any of these 34 constants.
        base = -3.0 + rng.normal(0.0, 0.5)
        # The two DECAY barriers are started from INDEPENDENT random draws
        # spanning their whole bound, deliberately NOT from the formation
        # barrier's start: if they end up near it, that is a result. The two
        # drift constants are started uniformly inside their own bounds.
        x0 = np.concatenate([
            np.full(N_K, base),
            [120.0 + rng.normal(0.0, 20.0)],
            rng.uniform(DECAY_EA_BOUNDS[0], DECAY_EA_BOUNDS[1], size=2),
            [rng.uniform(*ACID_YIELD_BOUNDS)],
            [rng.uniform(*ARP_AMINE_PKA_BOUNDS)],
        ])
        x0 = np.clip(x0, lower, upper)
        t0 = time.time()
        x0, c0 = coordinate_descent(x0, lower, upper, quick, DESCENT_SWEEPS, rng)
        print(f"  start {start}: after descent cost {c0:.4f}  "
              f"{round(time.time() - t0, 1)}s", flush=True)
        result = least_squares(
            residuals, x0, bounds=(lower, upper), args=(quick,),
            method="trf", x_scale="jac", max_nfev=max_nfev,
            diff_step=3e-2, ftol=1e-6, xtol=1e-8, verbose=0,
        )
        cost = float(result.cost)
        trace.append({
            "start": start,
            "cost_after_coordinate_descent": c0,
            "cost": cost,
            "nfev": int(result.nfev),
            "seconds": round(time.time() - t0, 1),
            "status": int(result.status),
        })
        print(f"  start {start}: cost {cost:.4f}  nfev {result.nfev}  "
              f"{trace[-1]['seconds']}s", flush=True)
        if best is None or cost < best.cost:
            best = result
            # CHECKPOINT. The first two attempts at this fit were lost to the
            # container's memory ceiling before any report was written; the
            # best vector so far is dumped after every start so that a loss
            # costs one start rather than all of them. It is a recovery
            # artefact, not an input: nothing reads it back automatically.
            try:
                (ROOT / "scratch").mkdir(exist_ok=True)
                (ROOT / "scratch/b23_best_x.json").write_text(json.dumps({
                    "start": start, "cost": cost, "x": best.x.tolist(),
                }))
            except Exception:  # a checkpoint must never break a fit
                pass
    assert best is not None
    return {"result": best, "trace": trace, "bounds": (lower, upper)}


class _FrozenResult:
    """
    A least-squares result reconstructed at an ALREADY-FROZEN parameter vector,
    with NO optimisation.

    Used by ``--reuse-frozen``, which exists for one specific and legitimate
    situation: the network gained a step whose rate is DERIVED rather than
    fitted (the reversible thiol release, whose constant is k_forward / K from
    Stack 2018's measured equilibrium), and that step carries ZERO FLUX on
    every system in the fit panel because none of them contains a matrix
    electrophile pool. The frozen parameters are therefore still the parameters
    of the new network, and the only thing that needs regenerating is the
    report's description of it.

    THIS IS NOT A REFIT AND CANNOT BECOME ONE. No optimiser runs; the parameter
    vector is read from the existing report and returned unchanged. The mode
    ASSERTS that every fit-row prediction is unchanged and refuses to write if
    any moved, so it cannot silently launder a network change into the freeze.
    """

    def __init__(self, x: np.ndarray, quick: bool):
        self.x = np.asarray(x, dtype=float)
        r0 = residuals(self.x, quick)
        self.fun = r0
        self.cost = 0.5 * float(np.dot(r0, r0))
        self.nfev = 1
        self.status = 0
        step = 3e-2
        jac = np.empty((len(r0), len(self.x)))
        for i in range(len(self.x)):
            xp = self.x.copy()
            xp[i] += step
            jac[:, i] = (residuals(xp, quick) - r0) / step
        self.jac = jac


def intervals(result, n_residuals: int) -> Dict[str, Any]:
    """Gauss-Newton 95% intervals, with flat directions reported as such."""
    jac = np.asarray(result.jac)
    dof = max(n_residuals - jac.shape[1], 1)
    chi2_red = 2.0 * float(result.cost) / dof
    out: Dict[str, Any] = {}
    try:
        jtj = jac.T @ jac
        covariance = np.linalg.pinv(jtj) * chi2_red
        sigmas = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    except Exception:
        sigmas = np.full(jac.shape[1], np.inf)
    names = list(PARAM_ORDER) + list(EXTRA_PARAM_NAMES)
    thresholds = {
        "Ea_lumped_formation": 60.0,
        "Ea_decay_thiol_sink": 60.0,
        "Ea_decay_carbonyl_sink": 60.0,
        "ph_acid_yield_per_sink_event": 0.5,
        "ph_arp_secondary_ammonium_pKa": 1.5,
    }
    for i, name in enumerate(names):
        s = float(sigmas[i])
        unidentified = (not np.isfinite(s)) or s > thresholds.get(name, 3.0)
        out[name] = {
            "value": float(result.x[i]),
            "ci95_halfwidth": None if unidentified else 1.96 * s,
            "identified": not unidentified,
        }
    out["_reduced_chi_square"] = chi2_red
    out["_dof"] = dof
    return out


# ===========================================================================
# REPORTING
# ===========================================================================


def main() -> int:
    parser = argparse.ArgumentParser()
    # TWO starts, seeded, ftol 1e-6 -- IDENTICAL to what B2.2 actually shipped
    # (its report says two, not three; see its diagnosis sec. 7). The container
    # memory ceiling that killed B2.2's first two attempts is real and a third
    # start was budgeted and not run there either. Changing the optimiser
    # budget between a fit and its refit would confound the comparison this
    # wave exists to make.
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--max-nfev", dest="max_nfev", type=int, default=250)
    parser.add_argument(
        "--reuse-frozen", action="store_true",
        help=("Regenerate the report at the ALREADY-FROZEN parameters without "
              "running any optimiser. Refuses to write if any fit-row "
              "prediction has moved. See _FrozenResult."),
    )
    parser.add_argument("--quick", action="store_true", default=True)
    parser.add_argument("--careful", dest="quick", action="store_false")
    args = parser.parse_args()

    print(f"B2.3 sulfur fit: {len(ACTIVE_FIT_ROWS)} declared FIT rows, "
          f"{N_K + N_EXTRA} free parameters", flush=True)
    out_json = data_paths.VALIDATION_DIR / "kinetic_core_b2_3_fit_report.json"
    previous = None
    if args.reuse_frozen:
        if not out_json.exists():
            raise SystemExit("--reuse-frozen needs an existing fit report")
        previous = json.loads(out_json.read_text())
        frozen_in = previous["frozen_parameters"]
        x_frozen = np.array(
            [frozen_in["log10_k_ref_at_145C"][k] for k in PARAM_ORDER]
            + [frozen_in["lumped_formation_Ea_kJ_mol"]]
            + [frozen_in["decay_Ea_kJ_mol"][f] for f in DECAY_FAMILY_ORDER]
            + [frozen_in["ph_drift"]["acid_yield_per_sink_event"],
               frozen_in["ph_drift"]["arp_secondary_ammonium_pKa"]]
        )
        print("  REUSING the frozen parameters; no optimiser will run.", flush=True)
        result = _FrozenResult(x_frozen, args.quick)
        fitted = {"result": result,
                  "trace": previous["objective"].get("multistart_trace", []),
                  "bounds": None}
    else:
        fitted = fit(args.starts, args.quick, args.max_nfev)
        result = fitted["result"]
    parameters = build_parameters(result.x)
    _fk, _fea, _dea, drift = unpack(result.x)
    predicted, extra = observables(parameters, args.quick, drift)
    ci = intervals(result, len(ACTIVE_FIT_ROWS))

    rows_out = []
    for row in ACTIVE_FIT_ROWS:
        p = float(predicted[row["id"]])
        t = float(row["target"])
        fold = (p / t) if (t > 0 and p > 0) else float("nan")
        rows_out.append({
            "id": row["id"],
            "system": row["system"],
            "kind": row["kind"],
            "observed": t,
            "predicted": p,
            "fold_error": fold,
            "log10_residual": (
                math.log10((p + FLOOR) / (t + FLOOR)) if np.isfinite(p) else None
            ),
            "sigma_log": row["sigma_log"],
            "source_anchor": row["anchor"],
            "note": row.get("note", ""),
        })

    frozen = {
        "log10_k_ref_at_145C": {k: float(result.x[i]) for i, k in enumerate(PARAM_ORDER)},
        "lumped_formation_Ea_kJ_mol": float(result.x[N_K]),
        "decay_Ea_kJ_mol": {
            family: float(result.x[N_K + 1 + i])
            for i, family in enumerate(DECAY_FAMILY_ORDER)
        },
        "decay_families": {k: list(v) for k, v in DECAY_FAMILIES.items()},
        "decay_keys_kept_on_formation_Ea": list(DECAY_KEYS_ON_FORMATION_EA),
        "ph_drift": drift.as_dict(),
        "reference_temperature_K": T_REF_S_K,
    }

    ph_state_rows = []
    for name, entry in extra["cache"].items():
        meta = entry["run"].metadata
        if meta.get("ph_mode") != "dynamic_charge_balance":
            continue
        ph_state_rows.append({
            "system": name,
            "buffer": meta["buffer"]["kind"],
            "buffer_declared_by_source": meta["buffer"]["declared_by_source"],
            "buffer_source": meta["buffer"]["source"],
            "label_pH": meta["ph_label_declared"],
            "in_situ_start": meta["ph_initial_in_situ"],
            "in_situ_end": meta["ph_final_in_situ"],
            "cooled_end": meta["ph_final_cooled"],
            "acid_pool_mmol_L": meta["acid_pool_final_mmol_L"],
            "SID_mmol_L": meta["strong_ion_difference_mmol_L"],
            "buffer_capacity_initial_mmol_per_pH":
                meta["buffer_capacity_initial_mmol_per_pH"],
            "fixed_point_residual_pH": meta["ph_fixed_point_residual"],
        })

    shares = {
        name: branch_shares(flux) for name, flux in extra["fluxes"].items()
    }

    payload: Dict[str, Any] = {
        "wave": "B2.3 -- charge-conservation fix (no new parameter), refit on the unchanged B2.2 FIT rows",
        "generated_by": "scripts/generators/generate_kinetic_core_b2_3_fit.py",
        "declaration": ("docs/reference/FIT_HOLDOUT_DECLARATION.md D.3, D.4, D.6 "
                        "and Amendment 7"),
        "network": describe_sulfur(),
        "objective": {
            "form": (
                "weighted least squares on log10 of the observable: "
                "r = log10((pred + 1e-12)/(obs + 1e-12)) / sigma_log, with "
                "one-sided treatment for the single CEILING row"
            ),
            "n_rows": len(ACTIVE_FIT_ROWS),
            "n_free_parameters": N_K + N_EXTRA,
            "new_free_parameters_vs_b2_1": 4,
            "ph_endpoint_rows_scored_in_pH_units_not_folds": True,
            "final_cost": float(result.cost),
            "reduced_chi_square": ci["_reduced_chi_square"],
            "dof": ci["_dof"],
            "multistart_trace": fitted["trace"],
            "starts_are_random_inside_wide_bounds": True,
        },
        "frozen_parameters": frozen,
        "parameter_intervals": {
            k: v for k, v in ci.items() if not k.startswith("_")
        },
        "rows": rows_out,
        "branch_shares_at_the_fit": shares,
        "registry_metadata": sulfur_registry_metadata(
            {k: v for k, v in parameters.items()
             if hasattr(v, "channel")}
        ),
        "out_of_scope": [dict(r) for r in OUT_OF_SCOPE],
        "decay_barrier_families": {
            family: {
                "keys": list(keys),
                "fitted_Ea_kJ_mol": float(
                    result.x[N_K + 1 + DECAY_FAMILY_ORDER.index(family)]),
                "identified_by": DECAY_FAMILY_IDENTIFYING_ROWS[family],
            }
            for family, keys in DECAY_FAMILIES.items()
        },
        "ph_state": {
            "provenance": PH_STATE_PROVENANCE,
            "per_system": ph_state_rows,
        },
        "declaration_contradiction_found": {
            "what": (
                "Cerny 2007 Table 4's MFT column (FIT) and Cerny 2007 Table 5's "
                "1x arm (HOLD-OUT) are the SAME NUMBER: 85% thiamine-derived / "
                "85:15. The declaration's disjointness rule 1 is violated by "
                "that pair."
            ),
            "resolved_how": (
                "Table 4's MFT column is EXCLUDED from this objective. The fit "
                "uses Cerny 2007b's full-ternary 54:46 (a different paper and a "
                "different system composition) and Table 4's non-MFT species. "
                "Reported, not resolved silently."
            ),
        },
    }

    if previous is not None:
        # THE TOLERANCE IS THE INTEGRATOR'S OWN, NOT ZERO. Adding a reaction
        # changes the solver's internal step sequence, so two runs of the same
        # chemistry differ at the level of the ODE tolerance (rtol 1e-5 in the
        # quick mode this is run in). The question the guard must answer is
        # "did any FLUX change", not "is the arithmetic bit-identical", so the
        # threshold is 1e-4 relative -- ten times the solver's own resolution,
        # and four orders below the smallest chemically meaningful change.
        tolerance = 1e-4
        moved, worst = [], 0.0
        for row, before in zip(rows_out, previous["rows"]):
            a, b = float(row["predicted"]), float(before["predicted"])
            scale = max(abs(a), abs(b), 1e-30)
            relative = abs(a - b) / scale
            worst = max(worst, relative)
            if relative > tolerance:
                moved.append((row["id"], b, a))
        payload["reuse_frozen_verification"] = {
            "rows_checked": len(rows_out),
            "tolerance_relative": tolerance,
            "worst_relative_move": worst,
            "rows_that_moved": [
                {"id": i, "before": b, "after": a} for i, b, a in moved
            ],
            "what_this_proves": (
                "The network step added after the freeze -- the reversible "
                "thiol release, whose rate is DERIVED from Stack 2018's "
                "measured equilibrium and not fitted -- carries zero flux on "
                "every system in the fit panel, because none of them contains "
                "a matrix electrophile pool. The frozen parameters are "
                "therefore unchanged AND still optimal for the new network."
            ),
        }
        if moved:
            raise SystemExit(
                f"--reuse-frozen: {len(moved)} fit-row predictions MOVED "
                f"({moved[:3]}...). The network change is NOT inert on the fit "
                f"panel, so the freeze does not carry over and a real refit is "
                f"required. Refusing to write."
            )
        print(f"  verified: all {len(rows_out)} fit rows unchanged; worst "
              f"relative move {worst:.2e} against a solver rtol of 1e-5.",
              flush=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, default=str))

    lines: List[str] = []
    a = lines.append
    a("# Kinetic core, Build Wave B2.3 -- the charge-conservation fix, refitted")
    a("")
    a("Modules 1 (sulfur formation) and 2 (thiol consumption) of "
      "`docs/reference/FIT_HOLDOUT_DECLARATION.md`.")
    a("")
    a(f"- network: **{payload['network']['n_species']} species, "
      f"{payload['network']['n_reactions']} reactions** "
      f"({payload['network']['trunk_reactions']} inherited from B1's trunk, "
      f"{payload['network']['sulfur_reactions']} new), carbon, nitrogen AND "
      f"sulfur balance enforced at import")
    a(f"- objective: **{len(ACTIVE_FIT_ROWS)} declared FIT rows**, "
      f"**{N_K + N_EXTRA} free parameters** (FOUR more than B2.1: two decay "
      f"barriers and two pH-drift constants), final cost "
      f"**{float(result.cost):.3f}**, reduced chi-square "
      f"**{ci['_reduced_chi_square']:.2f}**")
    a(f"- fitted pH parameters: **0**. The pH shape is structural "
      f"(acid/base catalysis + measured H2S speciation + Zheng & Ho's measured "
      f"four-pH thermolysis set).")
    a(f"- branch-fraction constants: **0**. Every split is a ratio of "
      f"time-integrated mass-action fluxes.")
    a(f"- activation energies in the NAMED CONSUMPTION channels: **0**, by "
      f"policy (inventory sec. C.1 / B.7) -- unchanged.")
    a(f"- activation energies on the DECAY LUMPS: **two**, one per named "
      f"family, both fitted. `k_dimer_decay` and `k_h2s_loss` keep B2.1's "
      f"lumped formation Ea because nothing measures either at two "
      f"temperatures.")
    a("")
    a("## The two decay barriers")
    a("")
    a("| family | keys | fitted Ea (kJ/mol) | identified by |")
    a("|---|---|---:|---|")
    for family, keys in DECAY_FAMILIES.items():
        value = float(result.x[N_K + 1 + DECAY_FAMILY_ORDER.index(family)])
        a(f"| `{family}` | {', '.join('`' + k + '`' for k in keys)} | "
          f"{value:.1f} | {DECAY_FAMILY_IDENTIFYING_ROWS[family][:110]} |")
    a("")
    a(f"The lumped FORMATION barrier is **{float(result.x[N_K]):.1f} kJ/mol** "
      f"for comparison. A decay barrier BELOW it means the sinks lose ground "
      f"to formation as temperature rises (which is the direction B2.1's own "
      f"140 C diagnosis asked for); ABOVE it means they gain. Neither bound "
      f"was tightened toward either answer and both barriers were started "
      f"from independent uniform draws over 20-250 kJ/mol.")
    a("")
    a("## The pH-trajectory state")
    a("")
    a(f"Two calibrated constants: acid yield per sink event "
      f"**{drift.acid_yield:.3f}**, Amadori secondary-ammonium pKa "
      f"**{drift.arp_amine_pka:.2f}**. Both are calibrated on Amendment 7's "
      f"THREE declared drift anchors and on nothing else; every other "
      f"acid-base constant in `ph_state.py` is textbook and immovable.")
    a("")
    a("| system | buffer | label pH | in-situ start | in-situ end | cooled end "
      "| acid, mM | SID, mM | beta0, mM/pH |")
    a("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in ph_state_rows:
        mark = "" if r["buffer_declared_by_source"] else " ⚠assumed"
        a(f"| `{r['system']}` | {r['buffer']}{mark} | {r['label_pH']:.2f} | "
          f"{r['in_situ_start']:.2f} | {r['in_situ_end']:.2f} | "
          f"{r['cooled_end']:.2f} | {r['acid_pool_mmol_L']:.1f} | "
          f"{r['SID_mmol_L']:.1f} | "
          f"{r['buffer_capacity_initial_mmol_per_pH']:.1f} |")
    a("")
    n_identified = sum(
        1 for k, v in ci.items() if not k.startswith("_") and v["identified"]
    )
    ea_value = ci["Ea_lumped_formation"]["value"]
    ea_at_bound = (
        abs(ea_value - LUMPED_FORMATION_EA_BOUNDS[0]) < 0.5
        or abs(ea_value - LUMPED_FORMATION_EA_BOUNDS[1]) < 0.5
    )
    a("## Read this before reading anything else")
    a("")
    a(f"**{n_identified} of {N_K + N_EXTRA} parameters are individually "
      f"identified.** With {len(ACTIVE_FIT_ROWS)} rows against "
      f"{N_K + N_EXTRA} free "
      f"parameters the fit has {ci['_dof']} degrees of freedom, so the "
      f"row-by-row agreement below is NOT evidence that the model is right. "
      f"What the corpus determines is a set of RATIOS and SHAPES -- branch "
      f"shares, MFT/FFT, dimer/monomer, the pH ordering -- and those are what "
      f"the hold-out report scores. No individual rate constant in this table "
      f"should be quoted as a measured quantity, cited elsewhere, or carried "
      f"into another module as if it were one.")
    a("")
    if ea_at_bound:
        a(f"**The single lumped formation activation energy landed ON ITS "
          f"BOUND at {ea_value:.1f} kJ/mol.** That is not a barrier estimate, "
          f"it is the fit saying it wants LESS temperature dependence than the "
          f"bound allows. The reason is visible in the corpus: the fit panel "
          f"spans 115-145 C but each temperature comes with a DIFFERENT "
          f"feedstock, buffer and lab, so temperature is confounded with system "
          f"and a cross-system Ea is not identifiable even in principle. "
          f"Reported as a failure to determine, not as a value.")
        a("")
    a("## Row-by-row")
    a("")
    a("| row | observed | predicted | fold | source |")
    a("|---|---:|---:|---:|---|")
    for r in rows_out:
        fold = r["fold_error"]
        fold_s = f"{fold:.2f}x" if np.isfinite(fold) else "n/a"
        a(f"| `{r['id']}` | {r['observed']:.4g} | {r['predicted']:.4g} | "
          f"{fold_s} | {r['source_anchor'][:70]} |")
    a("")
    a("## Parameters")
    a("")
    a("| parameter | log10 k(145 C) | 95% half-width | identified? |")
    a("|---|---:|---:|---|")
    for name in list(PARAM_ORDER) + list(EXTRA_PARAM_NAMES):
        entry = ci[name]
        hw = entry["ci95_halfwidth"]
        a(f"| `{name}` | {entry['value']:.3f} | "
          f"{'--' if hw is None else f'{hw:.3f}'} | "
          f"{'yes' if entry['identified'] else '**UNIDENTIFIED**'} |")
    a("")
    a("## Branch shares at the fitted point (computed, never stored)")
    a("")
    for name, s in shares.items():
        a(f"**{name}**: " + ", ".join(
            f"{k} = {v:.4g}" for k, v in s.items() if isinstance(v, float)
        ))
        a("")
    a("## Contradiction found in the declaration")
    a("")
    a(payload["declaration_contradiction_found"]["what"])
    a("")
    a(payload["declaration_contradiction_found"]["resolved_how"])
    a("")
    a("## Out of scope for this wave")
    a("")
    for row in OUT_OF_SCOPE:
        a(f"- **{row['lane']}** -- strands: {row['what_is_stranded']} "
          f"Reason: {row['why']}")
    a("")

    out_md = data_paths.VALIDATION_DIR / "kinetic_core_b2_3_fit_report.md"
    out_md.write_text("\n".join(lines))
    print(f"wrote {out_json}")
    print(f"wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
