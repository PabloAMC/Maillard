"""
src/kinetic_core/engine.py -- THE PROPAGATOR ENTRY POINT (Build Wave B5, 2026-08-29).

THE CUTOVER. This module is the single door between a user-facing formulation +
process specification and the mass-action kinetic core. Before B5 the core was a
calibration lane that nothing in the shipped path imported; from B5 it is the
shipped prediction path, and this file is where that binding happens.

WHAT THIS MODULE DOES
---------------------
Maps a ``FormulationSpec`` (precursors in mM) plus a ``ProcessSpec``
(temperature program, time, pH, matrix descriptor) onto ONE of the core's three
networks, integrates it, and emits the B4 output layer's objects: absolute
concentrations with reliability intervals, OAV tables, per-compound ratios
between formulations, rankings, and residual decompositions.

THE THING THIS MODULE EXISTS TO PREVENT
---------------------------------------
Silent numbers. The core is three DISJOINT networks over a NAMED species list;
a great many perfectly reasonable requests fall outside all three. Every such
request produces an explicit ``EnvelopeDeclaration`` with a named reason, and
NO NUMBER. Asking a declared-out prediction for an absolute raises
``OutOfEnvelope`` rather than returning a plausible-looking float.

THE FOUR LANES, AND HOW THEY COMPOSE
------------------------------------
Step counts below are as of B7 (2026-08-29). They are stated for orientation
only -- ``engine_metadata()`` COUNTS them at call time rather than quoting
these, because the literals that used to live there went stale across B6 and B7
and shipped wrong counts into every artifact for two waves (Q1).

  * ``TRUNK``      -- ``REACTIONS`` (26 steps). Glc/Fru/Gly -> melanoidins,
                      plus B7's furanic channel (HMF, DMHF).
                      No pH term, no a_w term.
  * ``SULFUR``     -- ``FULL_REACTIONS`` (93 steps) = trunk + sulfur. Adds the
                      pentoses, cysteine, thiamine, MFT/FFT/furfural. Carries a
                      pH trajectory.
  * ``ACRYLAMIDE`` -- ``FULL_ACRYLAMIDE_REACTIONS`` (42 steps) = trunk +
                      acrylamide. Adds asparagine and the acrylamide block.
  * ``LIPID``      -- B6's hydroperoxide pool and Frankel 1989's six-product
                      slate. This is the one lane that DOES compose: it
                      co-integrates with any ONE Maillard lane as a direct sum,
                      on the asserted-disjoint-species condition ``predict()``
                      checks at runtime.

The sulfur STEPS are deliberately absent from the acrylamide lane
(``acrylamide.OUT_OF_SCOPE``): composing them would spend the same cysteine
twice. A request whose targets span both lanes is therefore not a hard case, it
is an UNANSWERABLE case, and ``resolve_lane`` declares it rather than picking a
lane silently.

NO PARAMETERS LIVE HERE. Every constant is read from the frozen B1/B2.1/B3 fit
reports. This module fits nothing, tunes nothing, and contains no numeric
chemistry of its own.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Callable

import numpy as np

from src import data_paths

from . import operative_parameters
from .acrylamide import integrate_acrylamide
from .integrate import integrate
from .matrix_oav import (
    absolute_concentration,
    compare_formulations,
    decompose_residual,
    oav_table,
    predict_matrix_shift,
)
from .parameters import NETWORK_PH
from .parameters_acrylamide import MEASURED_ACRYLAMIDE, with_fitted_acrylamide
from .parameters_sulfur import (
    MEASURED_SULFUR, OX_AMBIENT_MMOL_L, OX_RESERVOIR_DEFAULT_UNITS, OX_SAT_MMOL_L,
    oxygen_parameters, with_fitted_sulfur,
)
from . import acrylamide_conditions, trunk_conditions
from .species import SPECIES_KEYS
from .species_acrylamide import ACRYLAMIDE_INDEX, acrylamide_ppb
from .species_sulfur import (
    MOLECULAR_WEIGHT_G_PER_MOL,
    SULFUR_INDEX,
    mmol_per_litre_to_ug_per_litre,
)
from .ph_state import BUFFER_ABSENT_WARNING, DEFAULT_BUFFER, BufferSpec, PhDrift
from .sulfur import integrate_sulfur

CELSIUS = 273.15

_B1_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b1_fit_report.json"
#: THE SULFUR LANE'S FROZEN PARAMETERS -- and this is a CUTOVER, stated here
#: rather than buried in a wave report. Build Wave B2.3 refits B2.2's own 48
#: parameters on B2.2's own 58 FIT rows after a CONSERVATION FIX (see
#: `ph_state.validate_charge_closure` and `sulfur.CENTRE_LEDGER`), so where a
#: B2.3 report exists it SUPERSEDES B2.2's: the B2.2 numbers were fitted
#: against a network that manufactured strong base out of bookkeeping, and
#: preferring them would be preferring a known defect. B2.2's report is kept
#: as the fallback so that every earlier artefact stays regenerable on a
#: checkout that has not run B2.3.
#: WAVE B8 (FIT_HOLDOUT_DECLARATION.md Amendments 16-18) supersedes B2.3 for
#: the same reason B2.3 superseded B2.2: preferring the predecessor would be
#: preferring a known defect. B2.3's vector carries
#: `Ea_decay_thiol_sink` = 216.1 kJ/mol, a value REFUTED by measurement --
#: Gigl 2021 measures the covalent-capture channel at k(333)/k(279) = 67.2 and
#: 216 kJ/mol predicts 4.6e6 for that ratio -- and it carries no barrier at all
#: on the two disulfide channels that Zhang 2026 k17 measures at 122.2.
#:
#: THE PROMOTION WAS DECLARED BLIND, in `kinetic_core_b8_prereg.md` sec. 2,
#: BEFORE any B8 score existed, and it is NOT contingent on a scorecard: what
#: B8 carries is four measured barriers and the removal of a refuted one. It
#: ships even where it scores worse, and where it scores worse the B8 hold-out
#: report says so. (It does score worse on the hold-out panel: 12/32 -> 8/30.)
_B2_FIT_REPORT_CANDIDATES = (
    data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json",  # 2026-09-03: fit/validate split
    data_paths.VALIDATION_DIR / "kinetic_core_b8_fit_report.json",
    data_paths.VALIDATION_DIR / "kinetic_core_b2_3_fit_report.json",
    data_paths.VALIDATION_DIR / "kinetic_core_b2_2_fit_report.json",
)
_B2_FIT_REPORT = next(
    (p for p in _B2_FIT_REPORT_CANDIDATES if p.exists()),
    _B2_FIT_REPORT_CANDIDATES[-1],
)
_B3_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b3_fit_report.json"


# ---------------------------------------------------------------------------
# Lanes
# ---------------------------------------------------------------------------

TRUNK = "trunk"
SULFUR = "sulfur"
ACRYLAMIDE = "acrylamide"
#: Build Wave B6. The lipid lane is the FOURTH lane and the FIRST one that
#: CO-INTEGRATES with a Maillard lane rather than conflicting with it. The
#: ruling and its condition live in ``lipid.lane_coupling_verdict`` and are
#: pre-registered in ``results/validation/kinetic_core_b6_prereg.md`` sec. 6.
LIPID = "lipid"

LANES: Tuple[str, ...] = (TRUNK, SULFUR, ACRYLAMIDE, LIPID)

#: The lanes that consume the same cysteine and therefore cannot compose.
MAILLARD_LANES: Tuple[str, ...] = (TRUNK, SULFUR, ACRYLAMIDE)


class OutOfEnvelope(RuntimeError):
    """
    Raised when a caller asks a DECLARED-OUT prediction for a number.

    Carrying this as an exception rather than a NaN is deliberate: a NaN
    propagates into a mean and disappears, and every one of this repository's
    documented accuracy defects began as a number that should not have existed.
    """

    def __init__(self, message: str, declaration: "EnvelopeDeclaration") -> None:
        super().__init__(message)
        self.declaration = declaration


# ---------------------------------------------------------------------------
# The species vocabulary -- the ONLY place a user-facing name becomes a species
# ---------------------------------------------------------------------------

#: Precursor synonyms -> core species key. Everything not in this table is an
#: unmapped precursor and produces a declaration, never a guess.
PRECURSOR_ALIASES: Mapping[str, str] = {
    "glucose": "Glc",
    "d-glucose": "Glc",
    "dextrose": "Glc",
    "fructose": "Fru",
    "d-fructose": "Fru",
    "glycine": "Gly",
    "gly": "Gly",
    "ribose": "PENT",
    "d-ribose": "PENT",
    "xylose": "PENT",
    "d-xylose": "PENT",
    "pentose": "PENT",
    "arabinose": "PENT",
    "cysteine": "Cys",
    "l-cysteine": "Cys",
    "cys": "Cys",
    "thiamine": "THI",
    "vitamin b1": "THI",
    "asparagine": "Asn",
    "l-asparagine": "Asn",
    "asn": "Asn",
    "glutamine": "Gln",
    "l-glutamine": "Gln",
    "lysine": "Lys",
    "l-lysine": "Lys",
    # B22 (2026-09-09): methionine, the Strecker substrate of the methionine chain, trunk lane
    "methionine": "MET",
    "l-methionine": "MET",
    # B24 (2026-09-09): proline, the 1-pyrroline source, trunk lane
    "proline": "PRO",
    "l-proline": "PRO",
    "1-pyrroline": "PYRL",
    "alanine": "Ala",
    "l-alanine": "Ala",
    "methylglyoxal": "MGO",
    # B13 (2026-09-07): the dicarbonyl trio
    "glyoxal": "GO",
    "glucosone": "G",
    # B37 (2026-09-11): the 3-deoxy series becomes CHARGEABLE, so that a fed-dicarbonyl pot can be
    # expressed at all. B13 made glyoxal, glucosone, diacetyl and methylglyoxal chargeable and left
    # these three as targets only; Mittelmaier et al. 2011 charge pure 3-DG at 120 C and follow
    # 3,4-DGE, which is the experiment docs/guides/EXPERIMENTS.md asks for by name and which no
    # spec could express until now. INERT: no bundle, benchmark, directional claim or fit row
    # charges any of these three (tests/unit/test_kinetic_core_b37.py holds that).
    "3-deoxyglucosone": "TDG",
    "3-dg": "TDG",
    "3,4-dideoxyglucosone": "DDG",
    "3,4-dideoxyglucosone-3-ene": "DDG",
    "3,4-dge": "DDG",
    "1-deoxyglucosone": "ODG",
    "1-dg": "ODG",
    # B39 (2026-09-11): the epimer Mittelmaier 2011 feeds
    "3-deoxygalactosone": "DGAL",
    "3-dgal": "DGAL",
    "diacetyl": "DA",
    "2,3-butanedione": "DA",
    "butane-2,3-dione": "DA",
    "norfuraneol": "NF",
    "amadori": "AMA",
    "arp": "ARP",
    # W6 (2026-09-07): the xylose-cysteine thiazolidine the sulfur lane carries; Zhai 2020 shows the
    # group's "Cys-Amadori" intermediate is ~94 % TTCA, so those names charge TTCA.
    "ttca": "TTCA",
    "2-threityl-thiazolidine-4-carboxylic acid": "TTCA",
    "2-(tetrahydroxybutyl)thiazolidine-4-carboxylic acid": "TTCA",
    "cys-amadori": "TTCA",
    "cysteine amadori": "TTCA",
    "cysteine-xylose amadori": "TTCA",
    "xylose-cysteine amadori": "TTCA",
}

#: Target-compound synonyms -> core species key.
TARGET_ALIASES: Mapping[str, str] = {
    # B13 (2026-09-07): the trunk's dicarbonyls, answerable on the trunk lane only
    "glyoxal": "GO",
    "glucosone": "G",
    "diacetyl": "DA",
    # B20 (2026-09-09): the glycation arm, trunk lane only, on protein-bound lysine
    "cml": "CML",
    "carboxymethyllysine": "CML",
    "n-epsilon-(carboxymethyl)lysine": "CML",
    "nε-(carboxymethyl)lysine (cml)": "CML",
    "nε-(carboxymethyl)lysine": "CML",
    "cel": "CEL",
    "carboxyethyllysine": "CEL",
    "n-epsilon-(carboxyethyl)lysine": "CEL",
    "nε-(carboxyethyl)lysine (cel)": "CEL",
    "nε-(carboxyethyl)lysine": "CEL",
    "fructosyl-lysine": "FLP",
    "fructosyllysine": "FLP",
    "fructoselysine": "FLP",
    "bound lysine": "LYSP",
    "protein-bound lysine": "LYSP",
    # B24 (2026-09-09): 2-acetyl-1-pyrroline, trunk lane only
    "2-acetyl-1-pyrroline": "AP",
    "2-acetyl-1-pyrroline (2-ap)": "AP",
    "2-ap": "AP",
    "acetylpyrroline": "AP",
    "1-pyrroline": "PYRL",
    # B22 (2026-09-09): the methionine chain, trunk lane only
    "methional": "MTAL",
    "3-(methylthio)propanal": "MTAL",
    "3-(methylthio)propionaldehyde": "MTAL",
    "methional (3-(methylthio)propanal)": "MTAL",
    "methanethiol from methional": "MSH",
    "dimethyl disulfide": "DMDS",
    "dimethyl disulfide (dmds)": "DMDS",
    "dmds": "DMDS",
    # B18 (2026-09-08): the pyrazine step, trunk lane only
    "pyrazine": "PZ",
    "2,5-dimethylpyrazine": "DMP",
    "2,5-dimethyl pyrazine": "DMP",
    "dimethylpyrazine": "DMP",
    "methylpyrazine": "MPZ",
    "2-methylpyrazine": "MPZ",
    "2,3-butanedione": "DA",
    "butane-2,3-dione": "DA",
    "methylglyoxal": "MGO",
    "acrylamide": "ACR",
    "2-furfurylthiol": "FFT",
    "2-furfurylthiol (fft)": "FFT",
    "furfurylthiol": "FFT",
    "fft": "FFT",
    "2-methyl-3-furanthiol": "MFT",
    "2-methyl-3-furanthiol (mft)": "MFT",
    "mft": "MFT",
    "bis(2-methyl-3-furyl) disulfide": "MFTD",
    "mft dimer": "MFTD",
    "furfural": "FUR",
    "2-furaldehyde": "FUR",
    # -- B6, the lipid lane. Frankel 1989's six-product slate, plus nonanal. --
    "hexanal": "HEXANAL",
    "n-hexanal": "HEXANAL",
    "nonanal": "NONANAL",
    "pentane": "PENTANE",
    "2,4-decadienal": "DECADIENAL",
    "trans,trans-2,4-decadienal": "DECADIENAL",
    "(e,e)-2,4-decadienal": "DECADIENAL",
    "methyl octanoate": "ME_OCTANOATE",
    "methyl 9-oxononanoate": "ME_9_OXONONANOATE",
    "methyl 13-oxo-9,11-tridecadienoate": "ME_13_OXO_TRIDECADIENOATE",
    # -- B28, 2026-09-09: the alkylfuran the lane refused until Frankel 1981 --
    "2-pentylfuran": "PENTYLFURAN",
    "2-pentyl furan": "PENTYLFURAN",
    "2-pentylfuran (pentyl furan)": "PENTYLFURAN",
    "pentylfuran": "PENTYLFURAN",
    "methanethiol": "MESH",
    "2-acetylthiazole": "ACTZ",
    "norfuraneol": "NF",
    "hydrogen sulfide": "H2S",
    "melanoidins": "MEL_N",
    # -- B7, the furanic channel. Both compounds left UNREPRESENTED_COMPOUNDS
    # in the same wave that gave them a route; the pre-B7 refusals were correct
    # and are quoted in the B7 report so the change of verdict is legible.
    "5-hydroxymethylfurfural": "HMF",
    "5-hydroxymethylfurfural (hmf)": "HMF",
    "5-hmf": "HMF",
    "hmf": "HMF",
    "dmhf": "DMHF",
    "hdmf": "DMHF",
    "furaneol": "DMHF",
    "2,5-dimethyl-4-hydroxy-3(2h)-furanone": "DMHF",
    "3,4-dideoxyglucosone": "DDG",
    "3-deoxyglucosone": "TDG",
    "3-dg": "TDG",
    "1-deoxyglucosone": "ODG",
    "1-dg": "ODG",
    "3-deoxygalactosone": "DGAL",
    "3-dgal": "DGAL",
    "acetylformoin": "AF",
}

#: Compounds a user may plausibly ask for that the core CANNOT NAME, each with
#: the reason. Being on this list is what turns a request into a declaration
#: instead of a KeyError, and the reason is what makes the declaration useful.
UNREPRESENTED_COMPOUNDS: Mapping[str, str] = {
    # -- B7: 5-HMF, DMHF and furaneol LEFT this list. Both pre-B7 refusals
    # were CORRECT at the time and are quoted verbatim in the B7 report:
    # "the hexose-dehydration route that forms it was never parameterised"
    # and "reporting NF as DMHF would be a species substitution the corpus
    # does not license". The first is now false (Kocadagli's amine-free
    # glucose system is ingested whole); the second is still true and is
    # honoured by DMHF being its OWN species with its own route, never an
    # alias of NF.
    #
    # HEMF / homofuraneol did NOT leave the list, and the reason is sharper
    # and different from the pre-B7 one -- exactly the discipline B6 used for
    # 1-hexanol and 2-pentylfuran.
    "hemf": (
        "2-ethyl-4-hydroxy-5-methyl-3(2H)-furanone (homofuraneol) needs a C2 "
        "Strecker donor -- alanine -- and the core cannot put alanine and a "
        "pentose in the same lane: the pentose lives on the sulfur lane and "
        "alanine only on the acrylamide lane, which do not compose. Blank 1997 "
        "measures HEMF at 6.8-10.0 ug/mmol in pentose/alanine systems and at "
        "0.3-1.3 in pentose/glycine ones -- a 5.2-25x PREFERENCE, not a switch "
        "(docs/reference/FIT_HOLDOUT_DECLARATION.md, amendment 12, which corrected amendment 8 on "
        "exactly this) -- so the compound is real, the route is understood, "
        "and the lane algebra is what refuses. Refused rather than answered "
        "with a DMHF number wearing a different name."
    ),
    "homofuraneol": (
        "see HEMF: the core cannot put alanine and a pentose in the same lane."
    ),
    "2-ethyl-4-hydroxy-5-methyl-3(2h)-furanone": (
        "see HEMF: the core cannot put alanine and a pentose in the same lane."
    ),
    "2,5-dimethyl-4-hydroxy-3(2h)-thiophenone": (
        "DMHF's ring-oxygen-to-sulfur swap IS a species (DMHFS) and its edge "
        "IS in the network, balanced -- but its RATE IS EXACTLY ZERO. Shu & Ho "
        "1988 is the only fed-precursor DMHF + cysteine experiment in the "
        "corpus and it reports a GC AREA PERCENT with no internal standard, no "
        "residual DMHF, no conversion and no molar yield of anything. Fitting "
        "a constant to its 6.0 % is a named prohibited derivation "
        "(k5b_dmhf_synthesis.md sec. 8.6, the thiol_addition_pentodiulose "
        "failure class). Haleva-Toledo 1999 would close it."
    ),
    # -- B6: hexanal and nonanal LEFT this list. The lipid lane forms both. --
    # 1-hexanol and 2-pentylfuran did NOT: the lane exists now, and the reason
    # they are still refused is sharper and different. A wave that un-refused
    # them would have invented two branch fractions.
    "1-hexanol": (
        "The lipid lane exists and forms the SIX products Frankel 1989 "
        "measured, but 1-hexanol is not one of them and NO aldehyde-reduction "
        "step is measured anywhere in the corpus -- in a thermally processed "
        "extrudate the reductant pool is not even identified. The retired screening "
        "lane emitted a number for it; this lane refuses. See "
        "parameters_lipid.PROHIBITED_DERIVATIONS."
    ),
    # 2-pentylfuran and "2-pentyl furan" LEFT this table in WAVE B28 (2026-09-09), came BACK the
    # same evening on a diagnosis that was wrong, and left again on 2026-09-10. The full record,
    # because the wrong step is the instructive one:
    #
    #   The old reason was "no branch fraction for the linoleate -> alkylfuran route is measured
    #   anywhere in the fit corpus". Frankel, Neff & Selke 1981 Table III measures it, in the same
    #   laboratory and by the same injector-port method as the slate this lane is fitted on, and it
    #   ships as PENTYLFURAN_PER_HEXANAL = 0.16 -- a ratio, so none of 1981's own denominator
    #   travels with it. Un-refusing it produced predictions about 1e5 below measurement. I read
    #   that as the lipid lane not being what produces these rows' hexanal, wrote that reasoning
    #   into five places, and restored the refusal.
    #
    #   It was a missing entry in _TARGET_LANE. Without one the concentration loop reports the
    #   species in mmol/L rather than ug/L. The lane makes 3.65e-5 mmol/L of it, which is 5.0 ug/L,
    #   against 163 measured -- a 32x miss, in line with this panel's median, not a degenerate
    #   answer at all. The lane was right, the ratio was right, and the dictionary was short one
    #   line.
    #
    #   WHAT THE EPISODE IS WORTH KEEPING FOR: a plausible mechanistic story explained a unit bug
    #   for a day. The ship rule's new degeneracy test caught that something was wrong and was
    #   right to; the diagnosis of WHY was mine and it was wrong.
    "propanal": (
        "The lipid lane forms no propanal. Propanal is an alpha-LINOLENATE "
        "scission product; Frankel 1989 fed linoleate only, so the FIT column "
        "contains no propanal share, and Schroen's 7 % is a property of "
        "RAPESEED OIL's fatty-acid profile rather than a transferable branch "
        "fraction."
    ),
    "2-nonenal": (
        "Named in Frankel 1989's introduction as the Hock partner of methyl "
        "9-oxononanoate, and quantified in none of his tables. No share can be "
        "fitted for it; see species_lipid.NAMED_UNQUANTIFIED_COPRODUCTS."
    ),
}

#: Which lane each target species is reachable in.
#: B34 (2026-09-11). THE SPECIES THAT ARE NOT MOLECULES, named rather than caught by an `else`.
#: These are elemental or lumped accounting pools -- a mole of "melanoidin nitrogen" is a mole of N
#: atoms, not of any compound -- so a molar mass would have to be invented. They are reported in
#: mmol/L on purpose. Everything else that reaches the reporting loop MUST have a molar mass; see
#: `_concentrations`, where a missing one now raises instead of silently changing the unit.
_REPORTED_IN_MMOL_PER_L: frozenset = frozenset({
    "MEL_C", "MEL_N", "MEL_S", "FRAG_C", "FRAG_N", "FRAG_S",
    "OX", "OXR", "OXV", "OLG", "MELE", "PROT_SS", "ACID", "CBX", "SB", "LYS_SITES",
})

_TARGET_LANE: Mapping[str, str] = {
    "ACR": ACRYLAMIDE,
    "FFT": SULFUR,
    "MFT": SULFUR,
    "MFTD": SULFUR,
    "FUR": SULFUR,
    "MESH": SULFUR,
    "ACTZ": SULFUR,
    "NF": SULFUR,
    "H2S": SULFUR,
    "MEL_N": TRUNK,
    # -- B7, the furanic channel. TRUNK, and that is load-bearing rather than
    # arbitrary: a TRUNK target adds NO lane requirement in ``resolve_lanes``,
    # so asking for HMF alongside acrylamide or alongside a thiol does not
    # create a lane conflict. The channel's parents are all trunk species and
    # the trunk network runs inside every lane, so HMF and DMHF are answerable
    # wherever their precursors are.
    "HMF": TRUNK,
    "DMHF": TRUNK,
    "DDG": TRUNK,
    "AF": TRUNK,
    # -- B13, the dicarbonyl trio: TRUNK species whose STEPS run only when the trunk
    # integrates on its own (network.TRUNK_REACTIONS); declare_envelope refuses them on
    # any other lane rather than returning the inert zero the sulfur state would carry.
    "G": TRUNK,
    "GO": TRUNK,
    "DA": TRUNK,
    "MGO": TRUNK,
    # -- B18, the pyrazine step: trunk-only as the dicarbonyls are
    "PZ": TRUNK,
    "DMP": TRUNK,
    "CML": TRUNK,
    "CEL": TRUNK,
    "FLP": TRUNK,
    "LYSP": TRUNK,
    "MTAL": TRUNK,
    "MSH": TRUNK,
    "DMDS": TRUNK,
    "AP": TRUNK,
    "PYRL": TRUNK,
    "MPZ": TRUNK,
    # -- B6, the lipid lane ------------------------------------------------
    "HEXANAL": LIPID,
    "NONANAL": LIPID,
    # WAVE B28, ADDED 2026-09-10 AND THIS OMISSION COST A WHOLE DIAGNOSIS. Without a lane here the
    # concentration loop falls through to its last branch and reports the species in mmol/L instead
    # of ug/L -- a factor of about 1.4e5 for this compound. Read as a prediction it looked like the
    # lipid lane making almost none of it, and a refusal was restored on that reading. It was a
    # missing dictionary entry.
    "PENTYLFURAN": LIPID,
    "PENTANE": LIPID,
    "DECADIENAL": LIPID,
    "ME_OCTANOATE": LIPID,
    "ME_9_OXONONANOATE": LIPID,
    "ME_13_OXO_TRIDECADIENOATE": LIPID,
}

#: B13: the species whose steps exist on the trunk integrator only.
# The trunk's optional arms and their target keys are one table now (trunk_arms.py); the names
# below are re-exported so nothing that imported them from here has to change.
from .trunk_arms import (  # noqa: E402
    DICARBONYL_TARGET_KEYS, GLYCATION_TARGET_KEYS, METHIONINE_TARGET_KEYS, PROLINE_TARGET_KEYS,
    PYRAZINE_TARGET_KEYS, TRUNK_ARMS, TRUNK_ONLY_TARGET_KEYS, named_targets,
)

#: Which lane each precursor species REQUIRES (absent = available in all lanes).
_PRECURSOR_LANE: Mapping[str, str] = {
    "PENT": SULFUR,
    "Cys": SULFUR,
    "THI": SULFUR,
    "ARP": SULFUR,
    "NF": SULFUR,
    "Asn": ACRYLAMIDE,
    "Gln": ACRYLAMIDE,
    "Lys": ACRYLAMIDE,
    "Ala": ACRYLAMIDE,
}

#: B31 (2026-09-10). THE LINE BETWEEN A COOK AND A HEADSPACE INCUBATION, expressed
#: as the fraction of the hydroperoxide pool that decomposes over the whole thermal
#: program. It is NOT a fitted quantity and NOT a tuned one.
#:
#: What it separates, computed on the panel's own conditions with the lane's own
#: anchored decomposition constant:
#:
#:     40 C, 10 min  (the four HS-SPME incubations)      3.826e-3
#:     140 C, 6 s    (Trikusuma UHT, the mildest cook)   0.2578
#:     160 C, 25 s   (Li 2026 extrusion)                 0.9994
#:     160 C, 30 min (Bi 2020 roasted pea)               1.000
#:
#: The gap between the first row and the second is a factor of 67 and NOTHING IN THE
#: PANEL LIES INSIDE IT. Thresholds of 0.01, 0.05 and 0.10 all give the identical
#: verdict on every row, which is what a threshold that is not doing any fitting looks
#: like. The verdict also survives the Q10 band end to end: at q10 = 2.0, the corner
#: that slows the hot pots most, the two sides are 2.824e-3 and 2.855e-2 -- still on
#: opposite sides of the line.
#:
#: The 1 % figure is chosen as the round number at the bottom of that empty gap, and
#: it is used ONLY to refuse, never to scale anything.
UNCOOKED_LOOH_CONVERSION_LIMIT = 0.01

#: B6. A LIPID CARRIER is not a precursor species: it is a matrix declaration
#: that resolves to a hydroperoxide pool through
#: ``parameters_lipid.LIPID_CARRIERS``, whose lipid fraction and peroxide value
#: are DECLARED ASSUMPTIONS with bands, not measurements. They are kept out of
#: ``mapped_precursors`` deliberately -- nothing may charge a Maillard network
#: with a protein isolate, which was the correct half of the pre-B6 refusal.
LIPID_CARRIER_ALIASES: Mapping[str, str] = {
    "pea protein isolate": "pea_protein_isolate",
    "pea protein": "pea_protein_isolate",
    "ppi": "pea_protein_isolate",
    "soy protein isolate": "soy_protein_isolate",
    "soy protein": "soy_protein_isolate",
    "spi": "soy_protein_isolate",
    "soy protein concentrate": "soy_protein_isolate",
    "methyl linoleate hydroperoxide": "frankel_pure_hydroperoxide",
    "methyl linoleate hydroperoxides": "frankel_pure_hydroperoxide",
    "linoleate hydroperoxide": "frankel_pure_hydroperoxide",
    "lipid hydroperoxide": "frankel_pure_hydroperoxide",
}


#: The compounds each lane can be asked to REPORT, in the display names the
#: engine's vocabulary maps. Used when a caller does not name targets.
LANE_DEFAULT_TARGETS: Mapping[str, Tuple[str, ...]] = {
    TRUNK: ("melanoidins", "5-HMF", "DMHF"),
    SULFUR: (
        "2-methyl-3-furanthiol (MFT)",
        "2-furfurylthiol (FFT)",
        "bis(2-methyl-3-furyl) disulfide",
        "furfural",
        "2-acetylthiazole",
        "methanethiol",
    ),
    ACRYLAMIDE: ("acrylamide",),
    LIPID: (
        "hexanal",
        "pentane",
        "2,4-decadienal",
        "methyl octanoate",
        "methyl 9-oxononanoate",
        "methyl 13-oxo-9,11-tridecadienoate",
        "nonanal",                # wave B28: answered, on a declared anchor
        "2-pentylfuran",          # wave B28, restored 2026-09-10 once the unit bug was found
    ),
}


#: Matrix descriptors that ARE the aqueous reference, under other names. A
#: free-amino-acid model system in buffer is water as far as an odour threshold
#: is concerned; a protein isolate is NOT, and is left alone so that
#: ``select_threshold`` returns its honest ``NoMeasuredThreshold``.
_AQUEOUS_MATRIX_ALIASES: Tuple[str, ...] = (
    "free", "water", "aqueous", "buffer", "none", "",
)


def resolve_matrix(descriptor: Optional[str]) -> str:
    """Map a spec's matrix descriptor onto a B4 threshold matrix, or pass through."""
    normalised = _norm(descriptor or "")
    if normalised in _AQUEOUS_MATRIX_ALIASES:
        return "water"
    return normalised


def _norm(name: str) -> str:
    return " ".join(str(name).strip().lower().split())


def default_targets_for(precursors: Mapping[str, float]) -> Tuple[str, ...]:
    """
    The compounds the core can report for this charge, when the caller names none.

    Resolves the lane from the precursors alone and returns that lane's
    reportable products. A charge that maps to no core species at all returns an
    empty tuple, which then produces an out-of-envelope declaration rather than
    an empty success.
    """
    keys = []
    carriers = []
    for name in precursors:
        key = PRECURSOR_ALIASES.get(_norm(name))
        if key is not None:
            keys.append(key)
            continue
        carrier = LIPID_CARRIER_ALIASES.get(_norm(name))
        if carrier is not None:
            carriers.append(carrier)
    if not keys and not carriers:
        return ()
    lanes, reasons = resolve_lanes(keys, [], carriers)
    if reasons or not lanes:
        return ()
    out: Tuple[str, ...] = ()
    for lane in lanes:
        out = out + LANE_DEFAULT_TARGETS.get(lane, ())
    # 2026-09-11 (review of PR #16): the docstring says "the compounds the core can report for this
    # charge", and the table is per LANE, so a ribose + cysteine pot used to be asked for
    # methanethiol (a methionine product) and answered it as 0.0. Filter each lane's list by what
    # the charge can actually reach; the lipid lane's defaults are the carrier's business.
    charged = {k: float(v) for k, v in zip(keys, (precursors[n] for n in precursors if PRECURSOR_ALIASES.get(_norm(n))))}
    reachable_by_lane = {
        lane: _reachable_species(lane, {k for k, v in charged.items() if v > 0.0}, None)
        for lane in lanes if lane in MAILLARD_LANES
    }
    kept = []
    for name in out:
        key = TARGET_ALIASES.get(_norm(name))
        lane = _TARGET_LANE.get(key)
        if lane == LIPID or key is None:
            kept.append(name)
            continue
        if any(key in reach for reach in reachable_by_lane.values()):
            kept.append(name)
    return tuple(kept)


# ---------------------------------------------------------------------------
# Specs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ThermalProgram:
    """
    T(t) as a piecewise-constant program: ``((duration_min, temperature_C), ...)``.

    Piecewise-constant rather than an arbitrary callable because that is what
    the integrators support without re-deriving their rate constants inside the
    right-hand side, and because every process specification this repository
    ingests (isothermal holds, extrusion zones, an autoclave ramp) is expressed
    as segments. A finer ramp is expressed as more segments; nothing is
    interpolated behind the caller's back.
    """

    segments: Tuple[Tuple[float, float], ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError("a thermal program needs at least one segment")
        for duration, _ in self.segments:
            if float(duration) < 0.0:
                raise ValueError("segment durations must be non-negative")

    @classmethod
    def isothermal(cls, temperature_c: float, time_min: float) -> "ThermalProgram":
        return cls(((float(time_min), float(temperature_c)),))

    @property
    def total_minutes(self) -> float:
        return float(sum(d for d, _ in self.segments))

    @property
    def peak_temperature_c(self) -> float:
        return float(max(t for _, t in self.segments))

    @property
    def min_temperature_c(self) -> float:
        return float(min(t for _, t in self.segments))

    def describe(self) -> str:
        if len(self.segments) == 1:
            d, t = self.segments[0]
            return f"isothermal {t:.1f} C for {d:g} min"
        return " -> ".join(f"{t:.1f} C x {d:g} min" for d, t in self.segments)


@dataclass(frozen=True)
class ProcessSpec:
    """Everything about the process that is not a precursor charge."""

    thermal: ThermalProgram
    ph: float = NETWORK_PH
    #: Measured FINAL pH, when the system is unbuffered and the source reports
    #: it. Only the sulfur lane can use it; declared as ignored elsewhere.
    ph_final: Optional[float] = None
    water_activity: Optional[float] = None
    #: B2.2: THE BUFFER IS NOW AN INPUT. ``None`` means "no buffer was
    #: declared", which resolves to ``ph_state.DEFAULT_BUFFER`` (unbuffered)
    #: and raises an extrapolation warning -- a pot whose buffer nobody
    #: recorded is a pot whose pH trajectory is being extrapolated. Supply
    #: ``BufferSpec(kind="clamped")`` to get B2's fixed-pH behaviour back
    #: explicitly rather than by accident.
    buffer: Optional[BufferSpec] = None
    #: The two calibrated pH-drift constants. ``None`` disables the dynamic pH
    #: state entirely, which is what keeps every B2/B2.1 artefact reproducible.
    ph_drift: Optional[PhDrift] = None
    #: Free-text matrix descriptor, matched against the B4 threshold matrices.
    matrix: str = "water"
    #: B11 (2026-09-07): the pot's physical state (`vessel.VesselSpec`, from the bundle's
    #: conditions.vessel). ``None`` = no vessel recorded: the sulfur lane charges the declared
    #: default reservoir and says so.
    vessel: Optional[Any] = None
    #: 2026-09-08 (the matrix layer): protein loading in g/L and, optionally, the isolate's own site
    #: densities in mmol per gram ({free_thiol_mmol_per_g, disulfide_mmol_per_g, amine_mmol_per_g}).
    #: With a named matrix on file, or with protein_sites, the loading charges the reactive-site
    #: pools (matrix_sites.resolve); without either, nothing is charged and the answer says so.
    protein_g_per_l: Optional[float] = None
    protein_sites: Optional[Mapping[str, float]] = None
    #: B31 (2026-09-10): what the pot STARTS with, in ug/L, for compounds the raw material carries
    #: in rather than the cook making. Added to the integrated concentration BEFORE the
    #: matrix-binding factor, because the protein cannot tell a carried molecule from a made one.
    #: ``None`` or absent means zero, so every pot that declares nothing is bit-for-bit unchanged.
    #: Only a level the source PRINTS as an unheated control of the same pot may be put here.
    carried_volatiles: Optional[Mapping[str, float]] = None
    #: B29 (2026-09-10): the pot's atmosphere -- "argon", "air" or "air_cu". ``None`` means air,
    #: which is what every fit row in this model was run in, so a spec that says nothing gets
    #: exactly the answer it got before the axis existed. Anything else with no fitted multiplier
    #: RAISES rather than quietly returning the air answer.
    atmosphere: Optional[str] = None

    @property
    def time_min(self) -> float:
        return self.thermal.total_minutes


@dataclass(frozen=True)
class FormulationSpec:
    """A named precursor charge, in mM."""

    name: str
    precursors: Mapping[str, float]
    process: ProcessSpec

    def __post_init__(self) -> None:
        for key, value in self.precursors.items():
            if float(value) < 0.0:
                raise ValueError(f"negative charge for {key!r}")


# ---------------------------------------------------------------------------
# The envelope
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EnvelopeDeclaration:
    """
    The engine's verdict on whether a request is answerable, and why.

    ``state`` is one of:
      * ``in_envelope``            -- every precursor and target is a species in
                                      one lane, and conditions are inside the
                                      parameters' measured range;
      * ``in_envelope_extrapolated`` -- answerable, but at conditions outside
                                      what the parameters license. A number is
                                      emitted AND the warnings are attached.
      * ``out_of_envelope``        -- not answerable. NO number is emitted.
    """

    state: str
    lane: Optional[str]
    reasons: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    unmapped_precursors: Tuple[str, ...] = ()
    unrepresented_targets: Tuple[Tuple[str, str], ...] = ()
    mapped_precursors: Mapping[str, float] = field(default_factory=dict)
    mapped_targets: Mapping[str, str] = field(default_factory=dict)
    #: 2026-09-11: requested names whose species no chain of the lane's reactions can reach from
    #: the charge, when OTHER requested targets can be. They are dropped from the answer and listed
    #: in `run_metadata["refused_targets"]` by name; a request where NONE is reachable is refused
    #: whole, in `reasons`, like any other out-of-envelope pot.
    unreachable_targets: Tuple[str, ...] = ()
    #: B6. Every lane this request needs. ``lane`` stays the PRIMARY (Maillard)
    #: lane so that every pre-B6 caller is unchanged; ``lanes`` is the tuple the
    #: propagator actually runs, and it has more than one member only for a
    #: CO-INTEGRATED lipid + Maillard request.
    lanes: Tuple[str, ...] = ()
    #: B6. The lipid carriers this charge declares, and where each came from.
    lipid_carriers: Tuple[str, ...] = ()

    @property
    def is_answerable(self) -> bool:
        return self.state in ("in_envelope", "in_envelope_extrapolated")

    @property
    def summary(self) -> str:
        """
        The verdict in ONE LINE, for a header, a log line or a card title. Q1.

        Every consumer of this object was writing its own version of this
        sentence, and they had drifted: the CLI said "REFUSED", the HTML report
        said "out of envelope", and the explain subcommand said neither. The
        wording is fixed here so that a refusal reads the same wherever it is
        printed, and so that an EXTRAPOLATED answer never renders as a plain
        answer just because a caller forgot to check ``warnings``.

        It is DERIVED, never stored: there is no state in which the summary and
        the fields it summarises can disagree.
        """
        lanes = ", ".join(self.lanes or ((self.lane,) if self.lane else ())) or "no lane"
        if self.state == "out_of_envelope":
            reason = self.reasons[0] if self.reasons else "no reason recorded"
            more = (
                f" (+{len(self.reasons) - 1} more)" if len(self.reasons) > 1 else ""
            )
            return f"REFUSED, no number emitted -- {reason}{more}"
        n_targets = len(self.mapped_targets)
        if self.state == "in_envelope_extrapolated":
            first = self.warnings[0] if self.warnings else "conditions outside the fit range"
            more = (
                f" (+{len(self.warnings) - 1} more)" if len(self.warnings) > 1 else ""
            )
            return (
                f"ANSWERED but EXTRAPOLATED on the {lanes} lane, {n_targets} "
                f"target(s) -- {first}{more}"
            )
        return f"answered in envelope on the {lanes} lane, {n_targets} target(s)"

    def as_dict(self) -> Dict[str, Any]:
        return {
            "state": self.state,
            "summary": self.summary,
            "lane": self.lane,
            "lanes": list(self.lanes or ((self.lane,) if self.lane else ())),
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "unmapped_precursors": list(self.unmapped_precursors),
            "unrepresented_targets": [
                {"compound": c, "reason": r} for c, r in self.unrepresented_targets
            ],
            "mapped_precursors": dict(self.mapped_precursors),
            "mapped_targets": dict(self.mapped_targets),
            "lipid_carriers": list(self.lipid_carriers),
        }


def resolve_lanes(
    precursor_keys: Sequence[str],
    target_keys: Sequence[str],
    lipid_carriers: Sequence[str] = (),
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """
    Every lane this request needs, or the reason no combination can carry it.

    B6 CHANGES THE RULE FOR EXACTLY ONE PAIR. The three Maillard lanes still
    refuse to compose with each other -- the acrylamide network deliberately
    omits every sulfur step, and summing them spends the same cysteine twice.
    The LIPID lane composes with any ONE of them, as a DIRECT SUM, because the
    species sets are disjoint and the only candidate coupling (the
    aldehyde-lysine covalent channel) is inert by ruling. The verdict is not
    hard-coded here: it is asked of ``lipid.lane_coupling_verdict`` on every
    call, so that enabling the covalent sink makes co-integration stop working
    rather than silently start double-counting.
    """
    required = set()
    for key in target_keys:
        lane = _TARGET_LANE.get(key)
        if lane is not None and lane != TRUNK:
            required.add(lane)
    for key in precursor_keys:
        lane = _PRECURSOR_LANE.get(key)
        if lane is not None:
            required.add(lane)
    if lipid_carriers:
        required.add(LIPID)

    maillard = sorted(required & set(MAILLARD_LANES))
    if len(maillard) > 1:
        return (), (
            "LANE CONFLICT: this request needs both the "
            + " and ".join(maillard)
            + " lanes at once. They do not compose -- the acrylamide network "
            "deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), "
            "because composing them would spend the same cysteine twice. No "
            "single integration can answer it.",
        )

    if LIPID in required and maillard:
        from .lipid import lane_coupling_verdict
        from .species import SPECIES_KEYS

        verdict = lane_coupling_verdict(list(SPECIES_KEYS))
        if not verdict["may_cointegrate"]:
            return (), (
                "LANE CONFLICT (lipid + " + maillard[0] + "): " + verdict["reason"],
            )
        return (maillard[0], LIPID), ()

    if LIPID in required:
        return (LIPID,), ()
    if not maillard:
        return (TRUNK,), ()
    return (maillard[0],), ()


def resolve_lane(
    precursor_keys: Sequence[str],
    target_keys: Sequence[str],
    lipid_carriers: Sequence[str] = (),
) -> Tuple[Optional[str], Tuple[str, ...]]:
    """
    The PRIMARY lane, or the reason none can carry this request.

    Unchanged for every pre-B6 request: with no lipid target and no lipid
    carrier this returns exactly what it returned before. When a lipid request
    is co-integrated with a Maillard lane, the Maillard lane is the primary
    (it owns the pH state, the buffer and the thermal warnings); a lipid-only
    request returns ``"lipid"``.
    """
    lanes, reasons = resolve_lanes(precursor_keys, target_keys, lipid_carriers)
    if reasons or not lanes:
        return None, reasons
    return lanes[0], ()


#: Warning tag the scorers look for (see :func:`unidentified_routes`).
HEXOSE_ENTRY_UNIDENTIFIED = "HEXOSE ENTRY UNIDENTIFIED"
#: Species keys of the sugars that reach the thiols only through the unidentified entry.
_HEXOSE_KEYS = ("Glc", "Fru")
#: Thiols whose only hexose route is that entry, and the products made from them (the two
#: disulfides and the methanethiol coupling product), which inherit the floor artefact: a ratio
#: of 1e27 for the dimer on a glucose arm is the same non-number as the thiol's (2026-09-09).
_HEXOSE_ENTRY_TARGETS = ("MFT", "FFT", "MFTD", "FFTD", "MMFT")


def unidentified_routes(
    mapped_precursors: Mapping[str, float], mapped_targets: Mapping[str, str]
) -> Tuple[str, ...]:
    """Target KEYS (``MFT``/``FFT`` and their disulfides and coupling product) whose formation
    from this charge runs only through the unidentified hexose entry: a hexose is charged, no
    pentose and no thiamine are, and the target is a thiol or made from one. Empty for every
    other request."""
    charged = {k for k, v in mapped_precursors.items() if float(v) > 0.0}
    if not any(k in charged for k in _HEXOSE_KEYS) or "PENT" in charged or "THI" in charged:
        return ()
    return tuple(sorted({key for key in mapped_targets.values() if key in _HEXOSE_ENTRY_TARGETS}))


def declared_unidentified(declaration: "EnvelopeDeclaration", compound: str) -> bool:
    """Whether ``compound`` (a bundle target name) is one the declaration flagged as running
    through the unidentified hexose entry."""
    if not any(str(w).startswith(HEXOSE_ENTRY_UNIDENTIFIED) for w in declaration.warnings):
        return False
    return declaration.mapped_targets.get(str(compound)) in _HEXOSE_ENTRY_TARGETS


def _lane_reactions(lane: str):
    """The reaction tuple the integrator will run for this lane."""
    if lane == SULFUR:
        from .sulfur import FULL_REACTIONS
        return FULL_REACTIONS
    if lane == ACRYLAMIDE:
        from .acrylamide import FULL_ACRYLAMIDE_REACTIONS
        return FULL_ACRYLAMIDE_REACTIONS
    from .network import TRUNK_REACTIONS
    return TRUNK_REACTIONS


def _ambient_seeds(lane: str, process) -> set:
    """Species `_integrate_program` charges on its own, without a precursor: the oxidant pool
    and its reservoir on the sulfur lane, and the protein pools when a loading is stated."""
    seeds = set()
    if lane == SULFUR:
        seeds |= {"OX", "OXR", "OXV"}
    try:
        from .matrix_sites import resolve as _resolve_sites
        charged, _ = _resolve_sites(process)
    except Exception:  # noqa: BLE001 - a malformed loading is reported by the matrix layer itself
        charged = None
    if charged is not None:
        if lane == TRUNK and charged.amine > 0:
            seeds.add("LYSP")
        if lane == SULFUR and charged.disulfide > 0:
            seeds.add("PROT_SS")
    return seeds


def _reachable_species(lane: str, charged: set, process) -> set:
    """Forward closure: every species some chain of the lane's reactions can make from `charged`."""
    reachable = set(charged) | _ambient_seeds(lane, process)
    reactions = _lane_reactions(lane)
    grew = True
    while grew:
        grew = False
        for r in reactions:
            if r.products and set(r.reactants) <= reachable and not set(r.products) <= reachable:
                reachable |= set(r.products)
                grew = True
    return reachable


def _unreachable_targets(lane: str, mapped_precursors, mapped_targets, process):
    """The requested names (caller's spelling) whose species no chain of the lane's reactions can
    reach from the charge. Lipid targets are the lipid lane's business and are never listed."""
    charged = {k for k, v in mapped_precursors.items() if float(v) > 0.0}
    reachable = _reachable_species(lane, charged, process)
    return sorted(
        name for name, key in mapped_targets.items()
        if _TARGET_LANE.get(key) != LIPID and key not in reachable
    )


def _carried_by_species(carried: Mapping[str, Any]) -> Dict[str, float]:
    """{species key: ug/L} from a `carried_volatiles` mapping, resolved through TARGET_ALIASES
    (the same table a target request comes in on). A name the table does not know is dropped
    HERE, and only here, so that every consumer sees the same declaration; a negative amount is
    a data error and is dropped too. Zero is kept: it is a declared 'not detected'."""
    out: Dict[str, float] = {}
    for name, amount in (carried or {}).items():
        key = TARGET_ALIASES.get(_norm(str(name)))
        if key is None or float(amount) < 0.0:
            continue
        out[key] = out.get(key, 0.0) + float(amount)
    return out


def declare_envelope(
    spec: FormulationSpec, targets: Sequence[str]
) -> EnvelopeDeclaration:
    """
    Decide, BEFORE integrating, whether this request is answerable.

    Everything that can make a request unanswerable is checked here and nowhere
    else, so that there is exactly one place to read to know what the core will
    refuse.
    """
    reasons: list = []
    warnings: list = []

    # --- precursors ------------------------------------------------------
    mapped_precursors: Dict[str, float] = {}
    unmapped: list = []
    carriers: list = []
    for name, value in spec.precursors.items():
        key = PRECURSOR_ALIASES.get(_norm(name))
        if key is None:
            carrier = LIPID_CARRIER_ALIASES.get(_norm(name))
            if carrier is not None:
                if carrier not in carriers:
                    carriers.append(carrier)
                warnings.append(
                    f"{name!r} is a LIPID CARRIER, not a precursor species. Its "
                    f"declared charge ({float(value):g}) is IGNORED -- 'mM of a "
                    f"protein isolate' has no defensible molar basis -- and the "
                    f"hydroperoxide pool comes instead from the carrier "
                    f"registry's declared lipid fraction and peroxide value, "
                    f"both of which are DECLARED ASSUMPTIONS with bands. It "
                    f"charges NO Maillard network: an isolate is still not a "
                    f"small-molecule precursor."
                )
                continue
            unmapped.append(str(name))
            continue
        mapped_precursors[key] = mapped_precursors.get(key, 0.0) + float(value)
    if unmapped:
        reasons.append(
            "UNMAPPED PRECURSORS "
            + ", ".join(repr(u) for u in sorted(unmapped))
            + ": not a species in any core lane. The core is a named "
            "small-molecule network; an intact protein, an isolate or a flour "
            "is not a precursor it can charge."
        )

    # --- targets ---------------------------------------------------------
    mapped_targets: Dict[str, str] = {}
    unrepresented: list = []
    for compound in targets:
        norm = _norm(compound)
        if norm in UNREPRESENTED_COMPOUNDS:
            unrepresented.append((str(compound), UNREPRESENTED_COMPOUNDS[norm]))
            continue
        key = TARGET_ALIASES.get(norm)
        if key is None:
            unrepresented.append(
                (
                    str(compound),
                    "not a species in any core lane, and not on the named "
                    "unrepresented-compound list either: the engine has no "
                    "vocabulary entry for it.",
                )
            )
            continue
        mapped_targets[str(compound)] = key
    if unrepresented:
        reasons.append(
            "UNREPRESENTED TARGETS: "
            + "; ".join(f"{c} -- {r}" for c, r in unrepresented)
        )

    # --- the lipid charge -------------------------------------------------
    # A lipid target with no carrier in the precursor list falls back to the
    # MATRIX descriptor, because that is where a real bundle records "this is a
    # pea protein isolate". The fallback is only consulted when a lipid target
    # was actually asked for, so it cannot switch a Maillard-only request onto
    # the lipid lane behind the caller's back.
    lipid_targets = [
        key for key in mapped_targets.values() if _TARGET_LANE.get(key) == LIPID
    ]
    if lipid_targets and not carriers:
        from_matrix = LIPID_CARRIER_ALIASES.get(_norm(spec.process.matrix or ""))
        if from_matrix is not None:
            carriers.append(from_matrix)
            warnings.append(
                f"no lipid carrier was named among the precursors; the MATRIX "
                f"descriptor {spec.process.matrix!r} was used instead. Its "
                f"lipid fraction and peroxide value are declared assumptions."
            )

    # --- routes the primary evidence does not identify --------------------
    # 2026-09-04 (after wave B9). Hexoses reach MFT and FFT only through the
    # fragmentation entry r_glc_c2c3 / r_glc_fur, whose rate constants no
    # step-level measurement in the corpus constrains; B9 (primary evidence
    # only) put them at the floor of their declared bands. A number computed
    # from a coordinate sitting on an arbitrary floor is not a prediction, so
    # a hexose-only charge asked for a thiol gets the number AND a declaration
    # that the scorers treat as NOT EVALUABLE. Thiamine has its own MFT route
    # (Bolton 1994), so a charge that carries thiamine is not affected; a
    # pentose charge is not affected because the intact-C5 route is fitted.
    unidentified = unidentified_routes(mapped_precursors, mapped_targets)
    if unidentified:
        warnings.append(
            f"{HEXOSE_ENTRY_UNIDENTIFIED} ({', '.join(sorted(unidentified))}): the only "
            "route from a hexose to these thiols is the C2+C3 fragmentation entry, "
            "whose rate constants no primary measurement identifies (the primary-evidence refit left them on "
            "their band floor). The number below is a floor artefact, not a fit; the "
            "scorecard and the envelope list this row as not evaluable, and the ordering "
            "'pentose above hexose' is the structural claim the model does support."
        )

    # --- lane ------------------------------------------------------------
    lanes, lane_reasons = resolve_lanes(
        list(mapped_precursors), list(mapped_targets.values()), carriers
    )
    lane = lanes[0] if lanes else None
    reasons.extend(lane_reasons)

    # B13: the dicarbonyl steps are trunk-only (the sulfur and acrylamide networks keep the
    # topology their fits were run on), so a dicarbonyl target on another lane is refused
    # by name instead of answered with the inert zero those state vectors carry.
    # The trunk's optional arms, one table (trunk_arms.py): not-shipped refusals, missing-precursor
    # refusals and the lane-conflict clause, in the orders the hand-written blocks emitted them.
    arm_targets = {arm.label: named_targets(arm, mapped_targets) for arm in TRUNK_ARMS}
    for arm in sorted(TRUNK_ARMS, key=lambda a: a.refusal_order):
        found = arm_targets[arm.label]
        if not found or lane != TRUNK:
            continue
        if arm.shipped is not None:
            shipped, why = arm.shipped()
            if not shipped:
                reasons.append(arm.label + " " + ", ".join(repr(c) for c in found) + ": " + why)
        if arm.required_precursors is not None and all(
                mapped_precursors.get(k, 0.0) <= 0.0 for k in arm.required_precursors):
            reasons.append(arm.label + " " + ", ".join(repr(c) for c in found) + arm.missing_precursor_message)
    # Two irregular checks stay explicit: this compound has no species key, so it is matched on the
    # raw target string ...
    if any(str(c).strip().lower() in ("dimethyl trisulfide", "dmts") for c in targets):
        from .parameters_methionine import METHIONINE_NO_DMTS_REASON

        reasons.append("UNREPRESENTED TARGET 'dimethyl trisulfide' (wave B22): " + METHIONINE_NO_DMTS_REASON)
    # ... and the glycation arm refuses on the MATRIX LAYER's charged amine sites, not on a precursor.
    glycation = arm_targets["GLYCATION TARGETS"]
    if glycation and lane == TRUNK:
        from .matrix_sites import resolve as _resolve_sites_for_glycation
        from .parameters_glycation import GLYCATION_NO_PROTEIN_REASON

        try:
            _charged, _ = _resolve_sites_for_glycation(spec.process)
        except Exception:  # noqa: BLE001 - a malformed loading is reported by the matrix layer itself
            _charged = None
        if _charged is None or _charged.amine <= 0:
            reasons.append(GLYCATION_NO_PROTEIN_REASON + " Targets: " + ", ".join(repr(c) for c in glycation) + ".")
    if lane == TRUNK and not lane_reasons:
        # 2026-09-11 (review of PR #16). Since B20 a stated protein loading charges the bound-lysine
        # pool on EVERY trunk run, and the glycation arm recycles glucose through fructosyl-lysine
        # back to 3-deoxyglucosone. That moves answers that never asked about glycation -- 5-HMF in
        # a glucose/glycine pot rises by about half at 30 g/L of pea isolate -- and the glycation
        # caveat was attached only when a glycation target was requested. Say it on every loaded
        # trunk answer instead.
        from .matrix_sites import resolve as _resolve_sites_for_loading

        try:
            _loaded, _ = _resolve_sites_for_loading(spec.process)
        except Exception:  # noqa: BLE001 - a malformed loading is reported by the matrix layer itself
            _loaded = None
        if _loaded is not None and _loaded.amine > 0 and not any(
            k in GLYCATION_TARGET_KEYS for k in mapped_targets.values()
        ):
            warnings.append(
                "A PROTEIN LOADING IS STATED, SO THE BOUND-LYSINE POOL IS CHARGED "
                f"({float(_loaded.amine):.3g} mmol/L of amine sites) and the glycation steps run on this "
                "trunk answer even though no glycation product was asked for: glucose is recycled "
                "through fructosyl-lysine to 3-deoxyglucosone, which raises the sugar-path products "
                "downstream of it. The same pot with no loading gives the unloaded number. Read the "
                "loading as an input that moved this answer, not as decoration."
            )
    if any(arm_targets.values()) and lane is not None and lane != TRUNK:
        named = [arm.label + " " + ", ".join(repr(c) for c in arm_targets[arm.label]) + f" (wave {arm.wave})"
                 for arm in sorted(TRUNK_ARMS, key=lambda a: a.conflict_order) if arm_targets[arm.label]]
        reasons.append(
            " and ".join(named)
            + f" run on the trunk lane only: the {lane} lane's network keeps the topology its fit was run "
            "on and carries these species inert. Ask for them in a sugar + amine pot that resolves to the "
            "trunk, or bring a measurement."
        )

    # --- the lipid lane's own refusals ------------------------------------
    if LIPID in lanes:
        from .parameters_lipid import LIPID_CARRIERS, oleate_fraction
        from .parameters_lipid_b28 import (
            OLEATE_MOLAR_ANCHOR_BAND, OLEATE_MOLAR_ANCHOR_CENTRE)

        if not carriers:
            reasons.append(
                "The lipid lane was selected but the charge declares NO LIPID "
                "CARRIER. Every product in this lane comes from a hydroperoxide "
                "pool, and the pool's size is an INPUT (an oxidation-state "
                "proxy): there is no route that makes a lipid aldehyde from a "
                "sugar or an amino acid. Name a carrier "
                f"({', '.join(sorted(set(LIPID_CARRIER_ALIASES.values())))}) or "
                "supply a peroxide value."
            )
        elif "NONANAL" in lipid_targets:
            oleate = max(
                oleate_fraction(LIPID_CARRIERS[c]) for c in carriers
                if c in LIPID_CARRIERS
            )
            if oleate > 0.0:
                # WAVE B28 (2026-09-09). This branch used to REFUSE, on the
                # ground that "the oleate -> nonanal branch fraction is measured
                # NOWHERE in the fit corpus". It is measured now (Frankel 1981
                # Table II), so the refusal has become a WARNING -- and the
                # warning is not decoration. The measured quantity is a SHARE of
                # a peak-area slate; turning it into an absolute needs a molar
                # anchor that no source supplies for oleate, so the answer rests
                # on a declared assumption with a wide band. Anyone reading a
                # nonanal number out of this lane must see that in the same
                # breath as the number.
                warnings.append(
                    "NONANAL RESTS ON A DECLARED ANCHOR, NOT A MEASURED YIELD. "
                    f"This matrix is {100.0 * oleate:.0f} % oleate by fatty-acid "
                    "share. The oleate -> nonanal SHARE is measured -- 15 % of "
                    "the slate (Selke 1978, republished by Frankel 1981) and "
                    "10 % (Frankel 1981's own photosensitized column, the only "
                    "independent determination). But Frankel 1981 prints PEAK "
                    "AREAS: no internal standard, no response factors, no "
                    "replicates. No absolute yield from an oleate hydroperoxide "
                    "exists anywhere in the corpus, so the absolute here assumes "
                    "the named-product molar yield per oleate hydroperoxide is "
                    f"{OLEATE_MOLAR_ANCHOR_CENTRE:g} times the measured one per "
                    "LINOLEATE hydroperoxide, banded "
                    f"{OLEATE_MOLAR_ANCHOR_BAND[0]:g} to "
                    f"{OLEATE_MOLAR_ANCHOR_BAND[1]:g}. Read the interval, not "
                    "the point. Frankel 1989's silence on nonanal remains a "
                    "declared hold-out and is still honoured: nonanal from a "
                    "LINOLEATE feed is exactly zero, by construction."
                )
        # -- B31 T3 (2026-09-10): A POT THAT WAS NEVER COOKED ------------------
        # Four panel pots hold at 40 C for ten minutes and are not cooks at all: the
        # 40 C / 10 min block is the HS-SPME headspace incubation, and each bundle's own
        # vessel provenance says so in as many words ("never heated", "an UNHEATED protein
        # powder", "No cook"). What they measure is what the raw material ARRIVED WITH.
        # A formation model asked to make 1260 ug/kg of hexanal out of a flour nobody
        # heated is not being tested on its chemistry, and scoring the miss as a chemistry
        # failure -- 3357x, 6078x, 3717x, 33392x -- misreports what is wrong.
        #
        # TWO INDEPENDENT DECLARATIONS HAVE TO AGREE before a row is refused, because
        # either alone is unsafe:
        #   1. the bundle's vessel says `closure: "no cook"`. This is a datum the bundles
        #      recorded months before this wave and it owes nothing to any prediction --
        #      but the string is OVERLOADED. Three hot bundles carry it meaning "no vessel
        #      to record" (a synthetic snapshot; two commercial products whose conditions
        #      block is a proxy operating point). It cannot be the whole rule.
        #   2. the thermal load cannot form what was measured: the fraction of the
        #      hydroperoxide pool that decomposes over the program is below 1 %. The four
        #      unheated pots sit at 3.8e-3; the coldest real cook in the panel (140 C for
        #      6 s) sits at 0.258, sixty-seven times higher. NOTHING IN THE PANEL LIES
        #      BETWEEN THEM, which is why the threshold is not a tuned knob: 0.01, 0.05 and
        #      0.10 all refuse the same seven rows. The verdict also survives the Q10 band
        #      -- at the worst corner (q10 = 2.0, which slows the hot pots most) the two
        #      sides are still 2.8e-3 and 2.9e-2, on opposite sides of the line.
        # The three hot "no cook" bundles fail clause 2 and are untouched.
        #
        # THE REFUSAL IS CONDITIONAL AND NAMES ITS OWN CURE. Declare what the pot started
        # with (`conditions.carried_volatiles`, this wave's other half) and the row is
        # answered. Trikusuma does exactly that and is scored. What is refused is the pot
        # for which NO source on this disk prints a starting state, and the honest report
        # of that is "cannot be asked", not a four-decade miss.
        #
        # This SHRINKS the panel: seven rows leave, all of them misses, so within-3x goes
        # 7/46 -> 7/39 and out-of-sample 6/45 -> 6/42 on arithmetic alone. That is exactly
        # the shape of a self-serving rule and is flagged here rather than buried. What
        # makes it not one: the criterion is condition-side and was fixed before any error
        # was looked at, and it leaves every lipid miss in a pot that WAS cooked standing
        # -- including the panel's largest, 2-pentylfuran at 366x in li 2026.
        if carriers and lipid_targets:
            vessel = getattr(spec.process, "vessel", None)
            closure = _norm(str(getattr(vessel, "closure", "") or ""))
            if closure == "no cook":
                from .parameters_lipid import k_looh_decomp_per_min

                exponent = sum(
                    k_looh_decomp_per_min(float(temperature_c)) * float(duration)
                    for duration, temperature_c in spec.process.thermal.segments
                )
                extent = 1.0 - math.exp(-exponent)
                if extent < UNCOOKED_LOOH_CONVERSION_LIMIT:
                    # A compound counts as declared when its SPECIES was declared, under any of
                    # the names the alias table accepts, and a declared 0.0 counts (2026-09-11:
                    # this compared raw strings and dropped zeros, so a level declared under one
                    # spelling and requested under another read as undeclared).
                    carried_declared = set(_carried_by_species(
                        getattr(spec.process, "carried_volatiles", None) or {}
                    ))
                    undeclared = sorted(
                        {
                            name for name, key in mapped_targets.items()
                            if _TARGET_LANE.get(key) == LIPID and key not in carried_declared
                        }
                    )
                    if undeclared:
                        reasons.append(
                            "THIS POT WAS NEVER COOKED, so what it measures is what the raw "
                            "material ARRIVED WITH, and this lane models FORMATION. The "
                            "bundle's own vessel says so (closure = 'no cook'), and the "
                            "physics agrees: over this thermal program only "
                            f"{100.0 * extent:.2f} % of the hydroperoxide pool decomposes, "
                            f"against {100.0 * UNCOOKED_LOOH_CONVERSION_LIMIT:.0f} % taken as "
                            "the floor for a cook and 25.8 % for the mildest real cook in the "
                            "panel. There is no thermal step here to model, so "
                            + ", ".join(undeclared)
                            + " is refused rather than answered with a formation from zero. "
                            "THE CURE IS A DECLARED STARTING STATE: put the level the source "
                            "prints for the unheated material in conditions.carried_volatiles "
                            "and the row is answered. Only a level the source PRINTS may go "
                            "there; nothing may be inferred from another paper's isolate."
                        )

    if lipid_targets and LIPID not in lanes and not lane_reasons:
        reasons.append(
            "a lipid product was requested but the lipid lane was not selected"
        )
    # -- B35 (2026-09-11): A POT WITH NO MAILLARD PRECURSOR ANSWERS NO MAILLARD TARGET ---------
    # A matrix-only charge declares a protein isolate and nothing else. The isolate is a LIPID
    # CARRIER and is deliberately kept out of `mapped_precursors` (see LIPID_CARRIER_ALIASES), so
    # the trunk, sulfur and acrylamide networks are integrated from an all-zero state and every
    # species in them stays zero BY CONSTRUCTION. Until this clause they were reported as 0.0 with
    # no refusal and no warning -- the audit found furaneol and furfural scored that way against
    # measurements of 2780 and 327 ug/kg, and 5-HMF would have been the same.
    #
    # THIS IS THE THIRD TIME THIS FAMILY OF BUG HAS COST A WAVE. B28 spent a day on 2-pentylfuran
    # reported in the wrong unit; B34 found the same silent unit fallback still catching two more
    # species; this is the same idea one level up -- an absence of a prediction dressed as one. The
    # repo's own words for it, from B28's record: "A near-zero is the absence of a prediction
    # dressed as one, so the refusal was restored with a sharper reason naming what would lift it."
    #
    # The lipid lane is exempt because its charge IS the carrier: it needs no free precursor.
    non_lipid_targets = sorted(
        name for name, key in mapped_targets.items() if _TARGET_LANE.get(key) != LIPID
    )
    if non_lipid_targets and not mapped_precursors and not lane_reasons:
        reasons.append(
            "THIS POT CHARGES NO PRECURSOR THAT COULD MAKE "
            + ", ".join(repr(c) for c in non_lipid_targets)
            + ". The charge declares only a matrix/lipid carrier, which is not a precursor: it "
            "resolves to a hydroperoxide pool for the lipid lane and charges NOTHING on the trunk, "
            "sulfur or acrylamide networks, so every species there is zero by construction rather "
            "than by prediction. Refused rather than answered with that zero. THE CURE IS A CHARGE: "
            "declare the sugar and amino acid this matrix brings to the cook, and the question "
            "becomes answerable."
        )

    # A target whose lane needs a precursor species this charge cannot supply.
    if lane is not None and not unmapped:
        if lane == SULFUR and not (
            # 2026-09-11: PENT was in this set. A pentose carries no sulfur, so a ribose-only
            # charge asked for a thiol was answered 0.0 instead of refused here.
            {"Cys", "THI", "ARP", "H2S", "TTCA"} & {k for k, v in mapped_precursors.items() if v > 0.0}   # W6: TTCA carries its cysteine sulfur
        ):
            if set(mapped_targets.values()) & {"MFT", "FFT", "MFTD", "MESH", "ACTZ"}:
                reasons.append(
                    "The sulfur lane was selected but the charge contains NO "
                    "sulfur source (no cysteine, thiamine or sulfide). Every "
                    "thiol in the core is built from a charged sulfur atom; "
                    "there is no route that makes one from a sugar alone."
                )
        if lane == ACRYLAMIDE and "Asn" not in mapped_precursors:
            if "ACR" in set(mapped_targets.values()):
                reasons.append(
                    "The acrylamide lane was selected but the charge contains "
                    "NO asparagine. Acrylamide in this network comes only from "
                    "the Asn + Glc initiation."
                )

    # --- 2026-09-11 (review of PR #16): THE ONE RULE THE THREE GUARDS WERE PROJECTIONS OF ------
    # B28 (a target reported in the wrong unit), B34 (a silent unit fallback) and B35 (a pot with
    # no precursor answering 0.0) were each patched where they bit. The invariant underneath all
    # three is that a requested target must be REACHABLE from what is charged, in the network of
    # the lane that will run: if no chain of the lane's reactions leads from the charged species
    # (plus the lane's ambient seeds) to the target's species, the integrator returns exactly
    # zero for it BY CONSTRUCTION, and that zero is not a prediction. The review found the same
    # thing waiting on every lane -- a cysteine-only pot answering the thiols, thiamine alone
    # answering furfurylthiol, asparagine alone answering acrylamide, glycine alone answering
    # HMF, and a zero-amount charge slipping past the B35 clause because its key was present.
    # One forward closure over the lane's reaction tuple catches all of them, and the next one.
    unreachable_targets: Tuple[str, ...] = ()
    if lane in MAILLARD_LANES and not unmapped and not lane_reasons and not any(
        "CHARGES NO PRECURSOR" in r for r in reasons
    ):
        unreachable = _unreachable_targets(lane, mapped_precursors, mapped_targets, spec.process)
        answerable = [n for n, k in mapped_targets.items() if _TARGET_LANE.get(k) != LIPID and n not in unreachable]
        if unreachable and answerable:
            # A mixed request (the CLI's default target list, say) answers what it can and refuses
            # the rest BY NAME rather than failing the whole pot for one compound it cannot make.
            unreachable_targets = tuple(unreachable)
            warnings.append(
                "NOT ANSWERED, BY NAME: " + ", ".join(repr(c) for c in unreachable)
                + f" -- in the {lane} lane's network no chain of reactions leads from what is charged "
                "to it, so its integrated value would be exactly zero by construction. It is left out "
                "of the answer and listed under refused_targets; the other targets are answered."
            )
        elif unreachable:
            charged = sorted(k for k, v in mapped_precursors.items() if v > 0.0)
            reasons.append(
                "THIS POT CHARGES NO PRECURSOR THAT COULD MAKE "
                + ", ".join(repr(c) for c in unreachable)
                + f": in the {lane} lane's network no chain of reactions leads from what is charged "
                + (f"({', '.join(charged)})" if charged else "(nothing above zero)")
                + " to it, so the integrator would return exactly zero by construction rather than "
                "by prediction. Refused rather than answered with that zero. THE CURE IS A CHARGE: "
                "declare the precursor this compound is made from."
            )

    # --- conditions ------------------------------------------------------
    peak = spec.process.thermal.peak_temperature_c
    low = spec.process.thermal.min_temperature_c
    if peak > 200.0 or low < 80.0:
        warnings.append(
            f"temperature program spans {low:.1f}-{peak:.1f} C; the integrator "
            f"is validated over 100-200 C and every operative rate constant was "
            f"measured over 80-120 C. This is a numerically sound extrapolation "
            f"of an experimentally unsupported barrier."
        )
    # B15 (2026-09-07): the acrylamide lane's declared initial-pH factor is printed by
    # acrylamide_conditions.declarations below, together with its a_w terms.
    if lane == TRUNK:
        # B12: the trunk carries a declared a_w term and a declared pH term (Amadori decay).
        warnings.extend(trunk_conditions.declarations(spec.process))
    if lane == ACRYLAMIDE:
        # B14/B15: the declared a_w terms (window 0.34-0.99) and the declared initial-pH factor.
        warnings.extend(acrylamide_conditions.declarations(spec.process))
    if spec.process.ph_final is not None and lane != SULFUR:
        warnings.append(
            f"a final pH was supplied but the {lane} lane has no pH trajectory; "
            f"it is ignored."
        )
    # B2.2: the buffer is an input with a declared default, and its ABSENCE is
    # an extrapolation rather than a silent assumption.
    if lane == SULFUR:
        # B11: the oxygen reservoir this run will be charged with, and where it came from.
        # The vessel's ABSENCE is an extrapolation only once it changes a rate, i.e. once the
        # shipped report consumes oxygen; with inert consumers the reservoir is inert too.
        reservoir, basis = oxygen_reservoir_units(spec.process)
        consumers = shipped_oxygen_consumers()
        if any(v > 0.0 for v in consumers.values()) and getattr(spec.process, "vessel", None) is None:
            warnings.append(
                f"OXYGEN RESERVOIR (B11): {reservoir:.0f} ambient units per litre of liquid "
                f"({basis}; 1 unit = {OX_SAT_MMOL_L:g} mmol/L dissolved at saturation). The "
                f"shipped report consumes oxygen (k_cys_ox {consumers['k_cys_ox']:.2e}, "
                f"k_red_ox {consumers['k_red_ox']:.2e} per unit per min), so the vessel is an "
                "input this run did not receive: declare it (conditions.vessel) to leave the flag."
            )
        if spec.process.buffer is None:
            warnings.append(BUFFER_ABSENT_WARNING)
        elif spec.process.buffer.is_clamped:
            warnings.append(
                "the buffer spec CLAMPS the pH. The dynamic pH state is "
                "switched off for this run and the declared pH is held for the "
                "whole hold, which is an assumption about the experiment, not "
                "a prediction about it."
            )


    # B6: the lipid lane is ALWAYS an extrapolation, and says so first.
    if LIPID in lanes:
        from .parameters_lipid import K_LOOH_DECOMP_ANCHOR, Q10_ASSUMPTION

        warnings.insert(0, Q10_ASSUMPTION.warning)
        peak = spec.process.thermal.peak_temperature_c
        warnings.insert(
            1,
            f"the lipid lane's rate anchor was measured at "
            f"{K_LOOH_DECOMP_ANCHOR.temperature_of_measurement_c:g} C and this "
            f"program peaks at {peak:.1f} C: "
            f"{Q10_ASSUMPTION.decades_of_extrapolation(peak):.1f} decades of "
            f"10 C, a factor of "
            f"{Q10_ASSUMPTION.factor(peak, Q10_ASSUMPTION.lo):.3g}-"
            f"{Q10_ASSUMPTION.factor(peak, Q10_ASSUMPTION.hi):.3g} on the rate.",
        )
        if abs(float(spec.process.ph) - float(K_LOOH_DECOMP_ANCHOR.ph_of_measurement)) > 1e-9:
            warnings.append(
                f"pH {spec.process.ph:g} was supplied; the lipid lane carries NO "
                f"pH term (its anchor is a single pH-6.7 emulsion). The pH is "
                f"recorded and IGNORED."
            )

    # The arms' own declarations, from the same table: each answer that names an arm's targets
    # carries that arm's caveats, and an amine the trunk charges as glycine says so.
    target_keys_here = set(mapped_targets.values())
    for arm in sorted(TRUNK_ARMS, key=lambda a: a.warning_order):
        if arm.charged_as_glycine is not None:
            pkey, on_trunk, on_other = arm.charged_as_glycine
            amount = mapped_precursors.get(pkey, 0.0)
            if amount > 0.0 and lane == TRUNK:
                warnings.append(on_trunk.format(amount=amount, lane=lane))
            if amount > 0.0 and on_other is not None and lane is not None and lane != TRUNK:
                warnings.append(on_other.format(amount=amount, lane=lane))
        if arm.target_caveats is not None and target_keys_here & arm.target_keys:
            warnings.extend(arm.target_caveats())
        if arm.label == "PYRAZINE TARGETS" and target_keys_here & (arm.target_keys | {"GO", "G"}):
            # B21's aqueous glyoxal supply is declared on glyoxal, glucosone AND the pyrazines, which
            # is a wider set than the pyrazine arm's own targets; it sits here so it keeps its place.
            from .parameters_dicarbonyl import AQUEOUS_GLYOXAL_CAVEAT

            warnings.append(AQUEOUS_GLYOXAL_CAVEAT)

    # --- B7: the furanic channel's own declarations -----------------------
    # Every one of these is an EXTRAPOLATION WARNING, not a refusal, and each
    # names the source that limits it. A caller who reads them knows exactly
    # what the number does and does not rest on.
    furanic_keys = set(mapped_targets.values()) & {"HMF", "DMHF", "AF", "DDG"}
    if furanic_keys:
        from .parameters_furanic import (
            FURANONE_EA_ASSUMPTION,
            HMF_SINK_NO_EXTRAPOLATION_ABOVE_K,
        )

        peak_k = spec.process.thermal.peak_temperature_c + CELSIUS
        if "HMF" in furanic_keys:
            warnings.append(
                "5-HMF: the two formation limbs are ingested WHOLE from "
                "Kocadagli & Gokmen 2016's AMINE-FREE amorphous glucose melt "
                "at 160-200 C. This program runs at "
                f"{spec.process.thermal.peak_temperature_c:.0f} C in an "
                "aqueous or matrix system, so both the temperature and the "
                "physical state are extrapolations. The furanic extraction dossier, sec. 6.2: that limb's "
                "activation energy reproduces four independent ways in the "
                "melt and COLLAPSES in all three real-matrix systems in the "
                "corpus."
            )
            warnings.append(
                "5-HMF: THE MODEL HAS NO VALIDATED SINK AT COOKING "
                "TEMPERATURE. The only audit-surviving HMF sink in the corpus "
                "(Hamzalioglu 2018, HMF + cysteine) is measured over 5-50 C "
                "and is CLAMPED at "
                f"{HMF_SINK_NO_EXTRAPOLATION_ABOVE_K - CELSIUS:.0f} C rather "
                "than extrapolated, and HMF self-degradation is a "
                "single-temperature 0.9 %-per-7-days control carried with no "
                "activation energy. The furanic extraction dossier's declared gap G2: the 50-150 C window "
                "is empty. EXPECT HMF TO BE OVER-PREDICTED."
            )
            if peak_k > HMF_SINK_NO_EXTRAPOLATION_ABOVE_K and (
                "Cys" in mapped_precursors
            ):
                warnings.append(
                    "5-HMF + cysteine: the sink constant is HELD at its 50 C "
                    "value for this whole program. Holding it UNDER-states the "
                    "sink; extrapolating it is a named prohibited derivation "
                    "(furanic extraction dossier, sec. 7.3), and the direction is stated rather than "
                    "chosen for convenience."
                )
        if "DMHF" in furanic_keys or "AF" in furanic_keys:
            warnings.append(str(FURANONE_EA_ASSUMPTION["warning"]))
            warnings.append(
                "DMHF: the LEVEL of the hexose route is a DECLARED TRANSFER "
                "from the pentose calibration. There is no absolute hexose "
                "DMHF yield in any of the five papers of the cluster -- the "
                "intact-C6 structure is settled twice over by CAMOLA and the "
                "magnitude is measured nowhere. Blank 1997's 39 cells are all "
                "pentose; Wang & Ho's nine are all per mole of methylglyoxal."
            )
            warnings.append(
                "DMHF: the Edge B (methylglyoxal, C3+C3) level is DIGITISED "
                "FROM A BAR CHART with no text layer, by external-standard "
                "HPLC with no recovery correction and an unstated pH hold -- "
                "three transmission defects deep, carried as a PRIOR ONLY. Its "
                "bracket (below detection in situ; 8-13 % in a real bean; 20 % "
                "at a 1.4 M methylglyoxal spike) is a hold-out, not a fit."
            )
            warnings.append(
                "DMHF: the CYSTEINE SINK (Edge C) is present, balanced, and "
                "runs at EXACTLY ZERO. No measurement of DMHF consumption "
                "exists anywhere; fitting one to Shu & Ho's 6.0 % GC area is a "
                "named prohibited derivation. Any DMHF number here is a "
                "FORMATION number with no sink."
            )

    if reasons:
        return EnvelopeDeclaration(
            state="out_of_envelope",
            lane=lane,
            lanes=tuple(lanes),
            reasons=tuple(reasons),
            warnings=tuple(warnings),
            unmapped_precursors=tuple(sorted(unmapped)),
            unrepresented_targets=tuple(unrepresented),
            mapped_precursors=mapped_precursors,
            mapped_targets=mapped_targets,
            lipid_carriers=tuple(carriers),
        )
    return EnvelopeDeclaration(
        state="in_envelope_extrapolated" if warnings else "in_envelope",
        lane=lane,
        lanes=tuple(lanes),
        warnings=tuple(warnings),
        mapped_precursors=mapped_precursors,
        mapped_targets=mapped_targets,
        lipid_carriers=tuple(carriers),
        unreachable_targets=unreachable_targets,
    )


# ---------------------------------------------------------------------------
# Frozen parameters
# ---------------------------------------------------------------------------


#: Parsed fit reports keyed by (path, mtime_ns, size): one stat per call instead of one
#: read + parse. 2026-09-03 (envelope cost): every predict re-read five reports from disk;
#: inside the Docker bind mount six envelope workers serialised on the file-sharing layer
#: and the pool gave no speed-up at all. A regenerated report (new mtime) is re-read.
_REPORT_CACHE: Dict[Tuple[str, int, int], Dict[str, Any]] = {}


def _read(path: Path) -> Dict[str, Any]:
    try:
        stat = path.stat()
    except FileNotFoundError:
        raise SystemExit(
            f"{path} not found. The engine never fits anything; it reads the "
            f"frozen fit reports. Regenerate them first."
        ) from None
    key = (str(path), stat.st_mtime_ns, stat.st_size)
    cached = _REPORT_CACHE.get(key)
    if cached is None:
        cached = json.loads(path.read_text())
        for stale in [k for k in _REPORT_CACHE if k[0] == key[0]]:
            del _REPORT_CACHE[stale]  # one live version per path
        _REPORT_CACHE[key] = cached
    return cached


def b1_fitted(variant: str = "variant_A_measured_sink") -> Dict[str, Tuple[float, float]]:
    """
    B1's four fitted trunk constants, as ``{key: (k_ref_100C, Ea)}``.

    ``variant_A_measured_sink`` is the default because it is what B2.1 and B3
    both inherited (their fit reports pin the identical four pairs); variant B
    is B1's out-of-sample browning variant and is offered for the browning lane.
    """
    frozen = _read(_B1_FIT_REPORT)["frozen_parameters"][variant]
    return {
        key: (float(v["k_ref_100C"]), float(v["ea_kj_mol"]))
        for key, v in frozen.items()
    }


def core_ph_drift() -> PhDrift:
    """
    B2.2's two FROZEN pH-drift constants, read from the fit report.

    THE ENGINE NEVER CONSTRUCTS ITS OWN. A caller may override the drift on a
    ProcessSpec (for a sensitivity study), but the shipped default is the
    frozen calibration and nothing else.
    """
    frozen = _read(_B2_FIT_REPORT)["frozen_parameters"]["ph_drift"]
    return PhDrift(
        acid_yield=float(frozen["acid_yield_per_sink_event"]),
        arp_amine_pka=float(frozen["arp_secondary_ammonium_pKa"]),
    )


_B7_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b7_fit_report.json"

#: The B1 variant every lane inherits (B2.1 and B3 pin the identical four pairs).
B1_VARIANT = "variant_A_measured_sink"


def frozen_parameters(lane: str) -> Dict[str, Any]:
    """
    The FIT-REPORT-SPACE parameter vector one lane reads, as ONE dict.

    Retirement step B2. The dict is shaped like the fit reports' own
    ``frozen_parameters`` blocks, merged: the B1 variant block under its own
    name (``variant_A_measured_sink``: ``{key: {k_ref_100C, ea_kj_mol}}``),
    then the lane's own block (B8's ``log10_k_ref_at_145C`` /
    ``lumped_formation_Ea_kJ_mol`` / ``decay_Ea_kJ_mol``, or B3's
    ``log10_k_ref_at_160C`` / ``fitted_Ea_kJ_mol``), and B7's ``k_dpo_af``.
    The keys are disjoint across reports, so one dict serves every lane.

    A Monte-Carlo draw perturbs THIS dict and hands it to
    :func:`core_parameters` -- so a draw moves the fit's own coordinates
    (one shared lumped Ea stays one number; ``k_odg_af`` keeps following
    ``k_dpo_af``; ``MEASURED_EA_OVERRIDES`` and ``NO_EA_KEYS`` are honoured
    by ``with_fitted_sulfur`` exactly as in the fit) rather than editing
    operative constants one by one.
    """
    if lane not in MAILLARD_LANES:
        raise ValueError(
            f"{lane!r} has no fit-report parameter vector; the lipid lane's "
            "frozen state is a branch model (core_lipid_model)."
        )
    out: Dict[str, Any] = {
        B1_VARIANT: {
            key: {"k_ref_100C": float(v["k_ref_100C"]), "ea_kj_mol": float(v["ea_kj_mol"])}
            for key, v in _read(_B1_FIT_REPORT)["frozen_parameters"][B1_VARIANT].items()
        }
    }
    if lane == SULFUR:
        frozen = _read(_B2_FIT_REPORT)["frozen_parameters"]
        out["log10_k_ref_at_145C"] = {
            k: float(v) for k, v in frozen["log10_k_ref_at_145C"].items()
        }
        out["lumped_formation_Ea_kJ_mol"] = float(frozen["lumped_formation_Ea_kJ_mol"])
        out["decay_Ea_kJ_mol"] = {
            k: float(v) for k, v in (frozen.get("decay_Ea_kJ_mol") or {}).items()
        }
        if frozen.get("formation_Ea_by_route_kJ_mol"):
            # B10: two barriers by route. The lumped value above is kept for every
            # reader that predates the split (it equals the sugar-trunk route).
            out["formation_Ea_by_route_kJ_mol"] = {
                k: float(v) for k, v in frozen["formation_Ea_by_route_kJ_mol"].items()
            }
        if frozen.get("oxygen"):
            out["oxygen"] = {k: float(v) for k, v in frozen["oxygen"].items()}
        if frozen.get("oxygen_log10_k"):
            out["oxygen_log10_k"] = {k: float(v) for k, v in frozen["oxygen_log10_k"].items()}
        if frozen.get("dimer_release_log10_k"):
            # B17: the disulfide-release constant, log10 (a shipped B17 report or a draw)
            out["dimer_release_log10_k"] = {k: float(v) for k, v in frozen["dimer_release_log10_k"].items()}
        if frozen.get("mele_site_log10_yield"):
            # B17a: log10 of the electrophile-site yield per osone decayed (a shipped B17a report or a draw)
            out["mele_site_log10_yield"] = {k: float(v) for k, v in frozen["mele_site_log10_yield"].items()}
        if frozen.get("thiol_addition"):
            # B25: the thiols' addition to the deoxypentosones, log10 k at 145 C and its barrier
            out["thiol_addition"] = {k: float(v) for k, v in frozen["thiol_addition"].items()}
        if frozen.get("dicarbonyl_redox"):
            # B27: log10 of the oxidant yield per mercaptoketone formed (a shipped B27 report or a draw)
            out["dicarbonyl_redox"] = {k: float(v) for k, v in frozen["dicarbonyl_redox"].items()}
    if lane == ACRYLAMIDE:
        frozen = _read(_B3_FIT_REPORT)["frozen_parameters"]
        out["log10_k_ref_at_160C"] = {
            k: float(v) for k, v in frozen["log10_k_ref_at_160C"].items()
        }
        out["fitted_Ea_kJ_mol"] = {
            k: float(v) for k, v in frozen["fitted_Ea_kJ_mol"].items()
        }
    if _B7_FIT_REPORT.exists():
        out["k_dpo_af"] = float(_read(_B7_FIT_REPORT)["frozen_parameters"]["k_dpo_af"])
    return out


def core_parameters(
    lane: str, *, frozen: Optional[Mapping[str, Any]] = None
) -> Dict[str, Any]:
    """
    The full operative parameter set for one lane, from the frozen reports.

    ``frozen`` (retirement step B2) is a dict shaped like
    :func:`frozen_parameters`; any block it carries REPLACES the report's, any
    block it omits is read from the report. With ``frozen=None`` the result is
    byte-identical to what this function returned before B2: the B1 pairs are
    passed straight from the report to ``operative_parameters`` and no
    furanic block is touched.
    """
    if lane == LIPID:
        raise ValueError(
            "the lipid lane has no mass-action parameter dictionary: its "
            "frozen state is a BRANCH MODEL plus a rate ASSUMPTION. Call "
            "core_lipid_model() instead -- the distinction is the module's "
            "whole point."
        )
    if lane not in MAILLARD_LANES:
        raise ValueError(f"unknown lane {lane!r}")

    override = dict(frozen or {})

    if B1_VARIANT in override:
        b1 = {
            key: (float(v["k_ref_100C"]), float(v["ea_kj_mol"]))
            for key, v in override[B1_VARIANT].items()
        }
    else:
        b1 = b1_fitted()
    parameters = dict(operative_parameters(b1))

    if lane == SULFUR:
        report = None
        if not {"log10_k_ref_at_145C", "lumped_formation_Ea_kJ_mol",
                "decay_Ea_kJ_mol", "formation_Ea_by_route_kJ_mol", "oxygen",
                "oxygen_log10_k", "dimer_release_log10_k", "mele_site_log10_yield", "thiol_addition",
                "dicarbonyl_redox"} <= set(override):
            report = _read(_B2_FIT_REPORT)["frozen_parameters"]
        pick = lambda key: override[key] if key in override else report[key]  # noqa: E731
        # B10: a report (or a draw) that carries the two route barriers uses them;
        # a wave before B10 carries only the lumped value and every route gets it.
        routes: Dict[str, float] = {}
        if report is not None and report.get("formation_Ea_by_route_kJ_mol"):
            routes.update(report["formation_Ea_by_route_kJ_mol"])
        routes.update(override.get("formation_Ea_by_route_kJ_mol") or {})
        formation = routes if routes else pick("lumped_formation_Ea_kJ_mol")
        parameters.update(MEASURED_SULFUR)
        parameters.update(
            with_fitted_sulfur(
                pick("log10_k_ref_at_145C"),
                formation,
                pick("decay_Ea_kJ_mol"),
            )
        )
        # B11: the oxygen consumers from the report's "oxygen" block (or a draw's), else the
        # inert zero defaults MEASURED_SULFUR already carries.
        oxygen: Dict[str, float] = {}
        if report is not None and report.get("oxygen"):
            oxygen.update(report["oxygen"])
        oxygen.update(override.get("oxygen") or {})
        # the Laplace draw moves the fit's own coordinates, which are log10 (block oxygen_log10_k)
        for key, value in (override.get("oxygen_log10_k") or {}).items():
            oxygen[key] = 10.0 ** float(value)
        if oxygen:
            parameters.update(oxygen_parameters(
                k_cys_ox=float(oxygen.get("k_cys_ox", 0.0)),
                k_red_ox=float(oxygen.get("k_red_ox", 0.0)),
            ))
        # B17: the disulfide-release constant from the report's (or a draw's) log10 block; the inert
        # zero MEASURED_SULFUR carries otherwise.
        release: Dict[str, float] = {}
        if report is not None and report.get("dimer_release_log10_k"):
            release.update(report["dimer_release_log10_k"])
        release.update(override.get("dimer_release_log10_k") or {})
        if release:
            from .parameters_sulfur import dimer_release_parameters

            parameters.update(dimer_release_parameters(k_dimer_release=10.0 ** float(release["k_dimer_release"])))
        # B17a: the site yield from the report's (or a draw's) log10 block; k_mele_site = yield x k_osone_decay
        # at 145 C with the carbonyl-sink family's barrier. The inert zero MEASURED_SULFUR carries otherwise.
        site: Dict[str, float] = {}
        if report is not None and report.get("mele_site_log10_yield"):
            site.update(report["mele_site_log10_yield"])
        site.update(override.get("mele_site_log10_yield") or {})
        if site:
            from .parameters_sulfur import mele_site_parameters

            k_osone = 10.0 ** float(pick("log10_k_ref_at_145C")["k_osone_decay"])
            ea_family = float(pick("decay_Ea_kJ_mol")["carbonyl_sink"])
            parameters.update(mele_site_parameters(k_mele_site=(10.0 ** float(site["mele_site_yield"])) * k_osone,
                                                   ea_kj_mol=ea_family))
        # B25: the addition constant and its barrier from the report's (or a draw's) block; inert zero otherwise.
        add: Dict[str, float] = {}
        if report is not None and report.get("thiol_addition"):
            add.update(report["thiol_addition"])
        add.update(override.get("thiol_addition") or {})
        if add:
            from .parameters_sulfur import thiol_addition_parameters

            parameters.update(thiol_addition_parameters(k_add=10.0 ** float(add["log10_k_add_145C"]), ea_kj_mol=float(add["ea_add_kj_mol"])))
        # B27: the dicarbonyl redox yield from the report's (or a draw's) block. Applied LAST because it
        # rescales k_nf_mp3p, which the fitted block above has just set; inert (phi = 0) otherwise.
        redox: Dict[str, float] = {}
        if report is not None and report.get("dicarbonyl_redox"):
            redox.update(report["dicarbonyl_redox"])
        redox.update(override.get("dicarbonyl_redox") or {})
        if redox:
            from .parameters_sulfur import apply_dicarbonyl_redox

            apply_dicarbonyl_redox(parameters, float(redox["log10_ox_yield_per_mercaptoketone"]))
    if lane == ACRYLAMIDE:
        report = None
        if not {"log10_k_ref_at_160C", "fitted_Ea_kJ_mol"} <= set(override):
            report = _read(_B3_FIT_REPORT)["frozen_parameters"]
        pick = lambda key: override[key] if key in override else report[key]  # noqa: E731
        parameters.update(MEASURED_ACRYLAMIDE)
        parameters.update(
            with_fitted_acrylamide(
                pick("log10_k_ref_at_160C"), pick("fitted_Ea_kJ_mol")
            )
        )
    if "k_dpo_af" in override:
        # Only an EXPLICIT override touches the furanic block: the frozen
        # literal in parameters_furanic is asserted equal to the B7 report by
        # a unit test, so the default path leaves it exactly as
        # operative_parameters installed it.
        from .parameters_furanic import with_fitted_furanic

        parameters.update(with_fitted_furanic(float(override["k_dpo_af"])))
    if "disputed_sinks" in override:
        # ENV-B13 (2026-09-10): the four trunk sinks a second laboratory disputes, at drawn values.
        # The envelope had no prior row for any of them, so every published interval asserted them
        # with certainty -- including two whose own authors flagged them as decisions and which the
        # second laboratory refutes. No centre moves; the envelope integrates across the disagreement.
        from .parameters_dicarbonyl import with_disputed_sinks

        parameters.update(with_disputed_sinks(override["disputed_sinks"]))
    if "aqueous_glyoxal" in override:
        # B21: the two aqueous glyoxal-supply constants, log10 at 100 C (the fit generator's candidates).
        from .parameters_dicarbonyl import AQUEOUS_GLYOXAL_COORDINATES, with_aqueous_glyoxal

        b = override["aqueous_glyoxal"]
        parameters.update(with_aqueous_glyoxal(*[float(b[k]) for k in AQUEOUS_GLYOXAL_COORDINATES]))
    # B39 (2026-09-11): the fed 3-deoxy triangle. The block is ALWAYS applied -- with the three new
    # steps at k = 0 until the fit ships, so nothing moves -- and an override carries the fit
    # generator's candidates and the envelope's draws. Applied after the disputed-sink draws so a
    # fitted centre wins over a printed band on the same constant.
    from .parameters_dicarbonyl import FED_3DEOXY_PARAMETERS, with_fed_3deoxy

    parameters.update(FED_3DEOXY_PARAMETERS)
    if "fed_3deoxy" in override:
        parameters.update(with_fed_3deoxy({k: (None if v is None else float(v)) for k, v in override["fed_3deoxy"].items()}))
    if "proline" in override:
        from .parameters_proline import PROLINE_COORDINATES, with_fitted_proline

        b = override["proline"]
        parameters.update(with_fitted_proline(*[float(b[k]) for k in PROLINE_COORDINATES]))
    if "methionine" in override:
        # B22: the four methionine coordinates (the fit generator's candidates, a later draw).
        from .parameters_methionine import METHIONINE_COORDINATES, with_fitted_methionine

        b = override["methionine"]
        parameters.update(with_fitted_methionine(*[float(b[k]) for k in METHIONINE_COORDINATES]))
    if "glycation" in override:
        # B20, the same discipline: the frozen literals in parameters_glycation are the default; an
        # explicit block of the five log10 constants at 100 C replaces them.
        from .parameters_glycation import GLYCATION_COORDINATES, with_fitted_glycation

        b = override["glycation"]
        parameters.update(with_fitted_glycation(*[float(b[k]) for k in GLYCATION_COORDINATES]))
    if "pyrazine" in override:
        # B18, the same discipline: the frozen literals in parameters_pyrazine are the default;
        # an explicit block {log10_k_go_ak_100C, ea_go_ak_kj_mol, log10_k_mgo_ak_100C,
        # ea_mgo_ak_kj_mol} (the fit generator's candidates, a later draw) replaces them.
        from .parameters_pyrazine import with_fitted_pyrazine

        b = override["pyrazine"]
        parameters.update(with_fitted_pyrazine(
            float(b["log10_k_go_ak_100C"]), float(b["ea_go_ak_kj_mol"]),
            float(b["log10_k_mgo_ak_100C"]), float(b["ea_mgo_ak_kj_mol"]),
        ))
    return parameters


_B6_FIT_REPORT = data_paths.VALIDATION_DIR / "kinetic_core_b6_fit_report.json"


def core_lipid_model():
    """
    B6's FROZEN branch model plus the default hydroperoxide-pool composition.

    Returns ``(BranchModel, LOOHComposition)``. The composition default is
    Frankel's AUTOXIDATION column as fitted -- the closest thing the corpus has
    to "what an oxidising food lipid's hydroperoxide pool looks like". It is a
    FIT quantity, not an assumption.
    """
    from .lipid import LOOHComposition, branch_model_from_dict

    frozen = _read(_B6_FIT_REPORT)["frozen_parameters"]
    branch = branch_model_from_dict(frozen["branch_model"])
    cells = frozen["default_pool_composition"]
    composition = LOOHComposition(
        f13_ct=float(cells["LOOH_13_ct"]),
        f13_tt=float(cells["LOOH_13_tt"]),
        f9_ct=float(cells["LOOH_9_ct"]),
        f9_tt=float(cells["LOOH_9_tt"]),
    )
    return branch, composition


# ---------------------------------------------------------------------------
# The prediction
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CoreDraw:
    """
    ONE Monte-Carlo draw of everything the core lets a sampler move.

    Retirement step B2. Every field is ``None`` by default, and a draw whose
    fields are all ``None`` reproduces the deterministic prediction exactly;
    that is asserted by a unit test, not assumed.

    * ``maillard`` -- a FIT-REPORT-SPACE override, shaped like
      :func:`frozen_parameters`, handed to :func:`core_parameters`. ``None``
      means the frozen reports.
    * ``q10`` -- the lipid lane's Q10 literal (declared band [2, 3]).
    * ``lipid_fraction_scale`` / ``peroxide_scale`` -- multiplicative scales on
      each carrier's declared lipid mass fraction and peroxide value (1.0 is
      the declared centre; the result is clipped to the carrier's own band).
    * ``furanone_partition_ea_kj_mol`` -- the offset on the furanone
      PARTITION barrier (declared band +/-50 kJ/mol), applied through the same
      helper the corner re-integration uses.
    * ``ph_drift`` -- a ``PhDrift`` for the sulfur lane, or ``None`` for the
      frozen calibration. A spec's own ``ph_drift`` wins over the draw's.
    """

    maillard: Optional[Mapping[str, Any]] = None
    #: B12 (2026-09-07): the trunk's declared condition bands. ``trunk_aw_scale`` scales the
    #: water-activity multiplier's excess over 1 (band trunk_conditions.AW_SCALE_BAND);
    #: ``trunk_ph_exponent`` is the Amadori-decay pH exponent (band PH_EXPONENT_BAND).
    #: ``None`` = the declared centres.
    trunk_aw_scale: Optional[float] = None
    trunk_ph_exponent: Optional[float] = None
    #: B14 (2026-09-07): the acrylamide lane's declared flat a_w multiplier inside the measured
    #: window (band acrylamide_conditions.AW_SCALE_BAND). ``None`` = the declared centre, 1.0.
    acrylamide_aw_scale: Optional[float] = None
    #: B15: the elimination a_w deficit scale (band AW_ELIMINATION_SCALE_BAND) and the two pH exponents.
    acrylamide_aw_elimination_scale: Optional[float] = None
    acrylamide_ph_exponent_formation: Optional[float] = None
    acrylamide_ph_exponent_elimination: Optional[float] = None
    #: B11: a multiplicative scale on the sulfur lane's oxygen reservoir (declared band
    #: OX_RESERVOIR_SCALE_BAND, spanning the saturation band and the unrecorded-vessel band).
    oxygen_reservoir_scale: Optional[float] = None
    q10: Optional[float] = None
    lipid_fraction_scale: Optional[float] = None
    peroxide_scale: Optional[float] = None
    furanone_partition_ea_kj_mol: Optional[float] = None
    ph_drift: Optional[PhDrift] = None

    @property
    def is_centre(self) -> bool:
        """True when every field is ``None`` -- the deterministic point."""
        return all(
            getattr(self, name) is None for name in self.__dataclass_fields__
        )


@dataclass(frozen=True)
class CorePrediction:
    """
    One integrated formulation, with its envelope declaration attached.

    ``concentrations_ug_per_l`` is EMPTY when the declaration is
    ``out_of_envelope``. That is the point of the type: there is no state in
    which a refused request carries a number.
    """

    spec: FormulationSpec
    declaration: EnvelopeDeclaration
    concentrations_ug_per_l: Mapping[str, float] = field(default_factory=dict)
    species_mmol_per_l: Mapping[str, float] = field(default_factory=dict)
    run_metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def answered(self) -> bool:
        return self.declaration.is_answerable

    def require(self, compound: str) -> float:
        """The absolute for ``compound``, or raise. Never a silent fallback."""
        if not self.answered:
            raise OutOfEnvelope(
                f"{self.spec.name}: refused -- "
                + " | ".join(self.declaration.reasons),
                self.declaration,
            )
        if compound not in self.concentrations_ug_per_l:
            raise OutOfEnvelope(
                f"{self.spec.name}: {compound!r} was not among the requested "
                f"targets, or is not a core species.",
                self.declaration,
            )
        return float(self.concentrations_ug_per_l[compound])

    def absolutes(self) -> Dict[str, Any]:
        """
        Every answered concentration wrapped in its B4 reliability band.

        B6: a LIPID compound also carries the width of its three DECLARED
        ASSUMPTIONS (Q10, lipid fraction, peroxide value), computed by
        re-integration at both corners and added in quadrature with B4's
        measured reliability band. A lipid absolute therefore reports a much
        wider interval than a sulfur one, which is the honest difference
        between a lane whose rate is measured and a lane whose rate is not.
        """
        widths = dict(self.run_metadata.get("lipid_extra_decades") or {})
        # B7: the furanone edges carry NO activation energy from any source
        # (all five papers of the cluster are single-temperature), so their
        # partition barrier is a DECLARED ASSUMPTION and is priced the same
        # way B6 prices its Q10 -- by re-integrating at both corners.
        furanic = dict(self.run_metadata.get("furanic_extra_decades") or {})
        widths.update(furanic)
        # 2026-09-08: a per-laboratory calibration adds its response factor's uncertainty (half-width
        # in decades) to the compounds it scaled; see calibration.Calibration.apply_factors.
        calibrated = dict(self.run_metadata.get("calibration_extra_decades") or {})
        for compound, extra in calibrated.items():
            widths[compound] = math.hypot(float(widths.get(compound, 0.0)), float(extra))
        # 2026-09-08: the matrix layer's declared binding brackets, priced at their corners
        for compound, extra in dict(self.run_metadata.get("matrix_extra_decades") or {}).items():
            widths[compound] = math.hypot(float(widths.get(compound, 0.0)), float(extra))
        return {
            compound: absolute_concentration(
                value,
                via_partition=True,
                extra_decades=float(widths.get(compound, 0.0)),
                provenance=(
                    f"kinetic core {self.declaration.lane} lane"
                    + (
                        "; +declared-assumption band (furanone partition "
                        "barrier, +/-50 kJ/mol) sized by re-integration"
                        if compound in furanic else
                        "; +declared-assumption band (Q10, lipid fraction, "
                        "peroxide value) sized by re-integration"
                        if compound in widths else ""
                    )
                ),
            )
            for compound, value in self.concentrations_ug_per_l.items()
        }

    def oav(self, matrix: Optional[str] = None) -> Dict[str, object]:
        """
        The B4 OAV table, with intervals, in the spec's matrix.

        Keyed by SPECIES KEY, because that is how the B4 threshold tables are
        keyed (``MFT``, ``FFT``, ``FUR``, ``ACTZ``, ``MFTD``). Handing them a
        display name silently returns ``NoMeasuredThreshold`` for a compound
        that has one -- a wiring bug found and fixed during the B5 cutover.
        """
        from .keyspaces import keys_for

        # B6: feed the ALREADY-WIDENED AbsoluteConcentration, not the bare
        # float. odour_activity auto-wraps a float in B4's measured band alone,
        # which would silently drop the lipid lane's declared-assumption width
        # from the OAV interval -- the one place the honesty could leak out.
        wrapped = self.absolutes()
        by_key: Dict[str, Any] = {}
        for compound, value in wrapped.items():
            keys = keys_for(compound, self.declaration.mapped_targets)
            if keys.b4 is None:
                # B6: no structural record, so no threshold, no binding class
                # and no unsaturation gate. Dropped from the OAV table rather
                # than defaulted -- ``NO_B4_RECORD`` says why for each, and
                # ``interval_rows()`` still carries the compound's interval.
                continue
            by_key[keys.b4] = value
        return oav_table(
            by_key,
            matrix=matrix or resolve_matrix(self.spec.process.matrix),
            temperature_c=self.spec.process.thermal.peak_temperature_c,
        )

    def interval_rows(self, matrix: Optional[str] = None) -> Tuple[Dict[str, Any], ...]:
        """
        One row per answered compound, CARRYING ITS OWN INTERVAL. Q1.

        The report layer used to reconstruct a row's interval by looking the
        compound up in the OAV table. That silently loses the interval of every
        compound the OAV table drops -- which is the whole ``NO_B4_RECORD`` set,
        i.e. four of the lipid lane's seven products. A compound with no
        measured odour threshold still has a perfectly well-defined
        concentration interval, and refusing to print it was an accident of
        where the number was stored, not a statement about the evidence.

        So the interval is attached HERE, next to the point it belongs to, and
        the OAV table is consulted only for the OAV. ``oav`` is ``None`` only
        when the compound is not in the table; ``no_b4_reason`` then says why,
        in the words ``NO_B4_RECORD`` records.

        Ordered by descending concentration, like :meth:`ranking`.
        """
        from .keyspaces import keys_for

        table = self.oav(matrix) if self.answered else {}
        per_species = dict(table.get("per_species") or {})
        wrapped = self.absolutes()
        rows: list = []
        for compound, value in self.ranking():
            keys = keys_for(compound, self.declaration.mapped_targets)
            absolute = wrapped.get(compound)
            entry = per_species.get(keys.b4) if keys.b4 else None
            rows.append(
                {
                    "compound": compound,
                    "species_key": keys.species,
                    "b4_key": keys.b4,
                    "lane": _TARGET_LANE.get(keys.species) or self.declaration.lane,
                    "predicted_ug_per_l": float(value),
                    "interval_ug_per_l": (
                        [absolute.lo_ug_per_l, absolute.hi_ug_per_l]
                        if absolute is not None else [None, None]
                    ),
                    "band_x": absolute.band_x if absolute is not None else None,
                    "interval_provenance": (
                        absolute.provenance if absolute is not None else None
                    ),
                    "oav": dict(entry) if isinstance(entry, Mapping) else None,
                    "no_b4_reason": keys.no_b4_reason,
                }
            )
        return tuple(rows)

    def ranking(self) -> Tuple[Tuple[str, float], ...]:
        """Compounds ordered by descending concentration."""
        return tuple(
            sorted(
                self.concentrations_ug_per_l.items(),
                key=lambda kv: kv[1],
                reverse=True,
            )
        )

    def as_dict(self) -> Dict[str, Any]:
        return {
            "formulation": self.spec.name,
            "declaration": self.declaration.as_dict(),
            "concentrations_ug_per_l": dict(self.concentrations_ug_per_l),
            "species_mmol_per_l": dict(self.species_mmol_per_l),
            "run_metadata": dict(self.run_metadata),
        }


def _run_lipid_lane(
    spec: FormulationSpec,
    declaration: EnvelopeDeclaration,
    *,
    q10: Optional[float] = None,
    lipid_scale: Optional[float] = None,
    pv_scale: Optional[float] = None,
    corners: bool = True,
) -> Tuple[Dict[str, float], Dict[str, Any], Dict[str, float]]:
    """
    Run the B6 lipid lane and size its interval BY RE-INTEGRATION.

    B2: the three declared assumptions are ARGUMENTS. ``q10`` is the literal
    (``None`` = the declared default), ``lipid_scale`` and ``pv_scale`` scale
    each carrier's declared centre (``None`` = 1.0, and the scaled value is
    clipped to the carrier's own declared band). ``corners=False`` skips the
    two corner re-integrations -- a Monte-Carlo draw prices the bands by
    sampling them and must not ALSO price them by re-integration.

    THE INTERVAL IS NOT A NOMINAL WIDTH. The lipid lane's absolute scale rests
    on three declared assumptions -- the Q10, the carrier's lipid fraction and
    its peroxide value -- and the honest way to price them is to run the model
    at both corners of all three and report the span. That also exposes a real
    property of the kinetics: at process temperature the hydroperoxide pool is
    EXHAUSTED within the hold, so the Q10 band largely cancels and what is left
    is the pool band. A nominal width could not have shown that.
    """
    from .lipid import charge_from_carrier, integrate_lipid
    from .parameters_lipid import LIPID_CARRIERS, Q10_ASSUMPTION

    branch, composition = core_lipid_model()
    segments = list(spec.process.thermal.segments)
    carrier_keys = [c for c in declaration.lipid_carriers if c in LIPID_CARRIERS]
    if not carrier_keys:
        raise OutOfEnvelope(
            f"{spec.name}: the lipid lane ran with no carrier", declaration
        )

    def _run(q10_value, lipid_of, pv_of):
        state: Dict[str, float] = {}
        runs = []
        for key in carrier_keys:
            carrier = LIPID_CARRIERS[key]
            charge = charge_from_carrier(
                carrier, composition,
                lipid_fraction=lipid_of(carrier),
                peroxide_value_meq_per_kg=pv_of(carrier),
            )
            run = integrate_lipid(charge, segments, branch, q10=q10_value)
            runs.append(run)
            for species_key, value in run.state_mmol_per_l.items():
                state[species_key] = state.get(species_key, 0.0) + value
        return state, runs

    def _scaled(centre, lo, hi, scale):
        if scale is None:
            return centre
        return min(max(float(centre) * float(scale), float(lo)), float(hi))

    point, point_runs = _run(
        q10,
        lambda c: _scaled(c.lipid_mass_fraction, c.lipid_lo, c.lipid_hi, lipid_scale),
        lambda c: _scaled(c.peroxide_value_meq_per_kg, c.pv_lo, c.pv_hi, pv_scale),
    )

    extra_decades: Dict[str, float] = {}
    low: Dict[str, float] = {}
    high: Dict[str, float] = {}
    if corners:
        low, _ = _run(Q10_ASSUMPTION.lo, lambda c: c.lipid_lo, lambda c: c.pv_lo)
        high, _ = _run(Q10_ASSUMPTION.hi, lambda c: c.lipid_hi, lambda c: c.pv_hi)
        for key, value in point.items():
            lo, hi = low.get(key, 0.0), high.get(key, 0.0)
            if value > 0.0 and lo > 0.0 and hi > 0.0:
                extra_decades[key] = 0.5 * abs(math.log10(hi / lo))

    metadata = {
        "carriers": carrier_keys,
        "branch_model": branch.as_dict(),
        "pool_composition": composition.as_dict(),
        "q10_default": Q10_ASSUMPTION.default,
        "q10_band": [Q10_ASSUMPTION.lo, Q10_ASSUMPTION.hi],
        "interval_method": (
            "RE-INTEGRATION at both corners of the three declared assumptions "
            "(Q10, lipid fraction, peroxide value). Not a nominal width."
        ),
        "declared_assumption_decades": dict(extra_decades),
        "runs": [dict(r.metadata) for r in point_runs],
        "warnings": sorted({w for r in point_runs for w in r.warnings}),
        "refusals": {k: v for r in point_runs for k, v in r.refusals.items()},
        "lower_corner_mmol_per_l": low,
        "upper_corner_mmol_per_l": high,
    }
    if not corners:
        metadata["interval_method"] = (
            "corner re-integration SKIPPED (size_declared_bands=False): the "
            "declared assumptions are being sampled by the caller."
        )
    if q10 is not None or lipid_scale is not None or pv_scale is not None:
        metadata["draw"] = {
            "q10": q10, "lipid_fraction_scale": lipid_scale, "peroxide_scale": pv_scale,
        }
    return point, metadata, extra_decades


#: B7. The compounds whose absolute level rests on the furanone-partition
#: barrier -- a DECLARED ASSUMPTION, because no activation energy for any
#: furanone family exists in the accessible literature on any edge.
FURANONE_BANDED_KEYS: Tuple[str, ...] = ("DMHF", "AF")

#: The edges the assumption sits on. ``k_af_dmhf`` is NOT here: it inherits a
#: measured Ea (Martins' 1-DG -> acetic acid, corroborated by Knol 2010) and is
#: swept separately, over three decades, in the B7 fit report.
_FURANONE_PARTITION_EDGES: Tuple[str, ...] = (
    "k_dpo_af", "k_odg_af", "k_mgo_dmhf",
)


def _furanone_corner_parameters(
    parameters: Mapping[str, Any], ea_offset_kj_mol: float
) -> Dict[str, Any]:
    """The operative set with the furanone PARTITION barrier shifted."""
    from dataclasses import replace

    out = dict(parameters)
    for key in _FURANONE_PARTITION_EDGES:
        parameter = out.get(key)
        if parameter is None or getattr(parameter, "ea_kj_mol", None) is None:
            continue
        out[key] = replace(
            parameter,
            ea_kj_mol=float(parameter.ea_kj_mol) + float(ea_offset_kj_mol),
        )
    return out


def shipped_oxygen_consumers() -> Dict[str, float]:
    """The two B11 oxygen consumers as the shipped report carries them (zero = inert)."""
    frozen = _read(_B2_FIT_REPORT)["frozen_parameters"]
    block = frozen.get("oxygen") or {}
    return {"k_cys_ox": float(block.get("k_cys_ox", 0.0)), "k_red_ox": float(block.get("k_red_ox", 0.0))}


def oxygen_reservoir_units(process) -> Tuple[float, str]:
    """
    B11: the headspace oxygen reservoir in ambient units per litre of liquid, and its basis.

    From the process's vessel block when the fill and vessel volumes are stated and the
    atmosphere is air; the declared default (Hofmann 1998's 100 mL pot) with a basis that says
    so otherwise. An open vessel or a continuous process gets a LARGE reservoir (oxygen is not
    limited by a headspace): ten times the default, said so.
    """
    from . import vessel as _vessel

    spec = getattr(process, "vessel", None)
    if spec is None:
        return OX_RESERVOIR_DEFAULT_UNITS, "no vessel recorded: the declared default reservoir"
    if spec.atmosphere in ("open", "continuous_process"):
        return 10.0 * OX_RESERVOIR_DEFAULT_UNITS, f"{spec.atmosphere}: oxygen not limited by a headspace (10x the default)"
    if spec.atmosphere == "inert_gas":
        return 0.0, "inert-gas atmosphere stated: no reservoir"
    o2 = spec.o2_mmol()
    if o2 is None or spec.fill_mL is None or spec.fill_mL <= 0:
        return OX_RESERVOIR_DEFAULT_UNITS, f"vessel {spec.atmosphere} with volumes unstated: the declared default reservoir"
    mmol_per_l = o2 / (float(spec.fill_mL) / 1000.0)
    return mmol_per_l / OX_SAT_MMOL_L, f"{spec.headspace_mL:.0f} mL headspace over {spec.fill_mL:g} mL: {mmol_per_l:.1f} mmol O2 per litre"


def _integrate_program(
    lane: str,
    parameters: Mapping[str, Any],
    initial: Mapping[str, float],
    process: ProcessSpec,
    *,
    ph_drift: Optional[PhDrift] = None,
    trunk_draw: Optional[Tuple[Optional[float], Optional[float]]] = None,
    acrylamide_draw: Optional[Dict[str, Optional[float]]] = None,
    reservoir_scale: Optional[float] = None,
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """
    Integrate a piecewise-constant thermal program, chaining the state across
    segments, and return the FINAL state as ``{species_key: mmol/L}``.

    ``trunk_draw`` (B12) = ``(aw_scale, ph_exponent)`` from a ``CoreDraw``; ``None`` entries
    mean the declared centres.

    ``ph_drift`` (B2) is consulted only when the process declares none: the
    order is spec, then draw, then the frozen calibration.
    """
    state: Dict[str, float] = dict(initial)
    metadata: Dict[str, Any] = {"segments": [], "lane": lane}
    if lane == SULFUR:
        # B10 (2026-09-06): every fit system since B2.3 was integrated with the
        # ambient oxidant pool charged at OX_AMBIENT_MMOL_L; the engine charged
        # nothing, so the two oxidant channels carried flux in the fit and none
        # in use (B11 prereg sec. 2.1). Charged here so fit and deployment agree;
        # a caller that passes its own "OX" (a fit generator, wave B11's vessel
        # charge) is left alone. Effect on every panel row: below 1 % at trace
        # thiol (the 2026-09-06 probe).
        state.setdefault("OX", OX_AMBIENT_MMOL_L)
        # 2026-09-08 (the matrix layer): the protein disulfide pool from the spec's protein loading
        # and the matrix's site densities; zero, as before, when no loading is stated.
        if "PROT_SS" not in state:
            from .matrix_sites import resolve as _resolve_sites

            charged, _note = _resolve_sites(process)
            if charged is not None and charged.disulfide > 0:
                state["PROT_SS"] = float(charged.disulfide)
        # B11 (2026-09-07): the headspace reservoir, in ambient units per litre of
        # liquid, from the process's vessel block; the declared default otherwise.
        # Inert while the report's consumers are zero (every wave before B11).
        if "OXR" not in state:
            reservoir, _basis = oxygen_reservoir_units(process)
            state["OXR"] = reservoir * (1.0 if reservoir_scale is None else float(reservoir_scale))
        state.setdefault("OXV", 0.0)
    if lane != TRUNK:
        # 2026-09-11 (review of PR #16). Methionine, proline and 1-pyrroline are TRUNK species. On
        # the sulfur or acrylamide lane the declaration promises they are "recorded and not
        # charged" (trunk_arms.py) -- and then this function handed them to an integrator whose
        # state vector does not contain them, which raised KeyError("unknown species 'MET'") on any
        # cysteine + ribose + methionine pot asked for a thiol. Recorded means dropped here.
        for key in ("MET", "PRO", "PYRL"):
            if key in state:
                metadata.setdefault("recorded_not_charged", {})[key] = float(state.pop(key))
    if lane == TRUNK and state.get("PRO", 0.0) > 0.0:
        # B24: proline as the Amadori amine too, declared (kinetic_core_b24_prereg.md sec. 2).
        state["Gly"] = float(state.get("Gly", 0.0)) + float(state["PRO"])
    if lane == TRUNK and state.get("MET", 0.0) > 0.0:
        # B22 (2026-09-09): methionine is the Strecker substrate (MET) and, DECLARED, the amine of the
        # Amadori chemistry that makes the dicarbonyls, charged as glycine at the same molarity (the
        # same amine plays both roles in turn; kinetic_core_b22_prereg.md sec. 2, declaration v).
        state["Gly"] = float(state.get("Gly", 0.0)) + float(state["MET"])
    if lane == TRUNK and "LYSP" not in state:
        # B20 (2026-09-09): the bound-lysine pool from the spec's protein loading and the matrix's
        # amine density, times the declared available fraction (the band centre); zero, as before,
        # when no loading is stated, so the glycation steps carry no flux.
        from .matrix_sites import resolve as _resolve_sites
        from .parameters_glycation import available_fraction

        try:
            charged, _note = _resolve_sites(process)
        except Exception:  # noqa: BLE001 - a malformed loading was already refused by the declaration
            charged = None
        if charged is not None and charged.amine > 0:
            state["LYSP"] = float(charged.amine) * available_fraction(charged.amine_band)

    for index, (duration, temperature_c) in enumerate(process.thermal.segments):
        grid = np.array([0.0, float(duration)])
        if lane == SULFUR:
            run = integrate_sulfur(
                parameters,
                float(temperature_c) + CELSIUS,
                state,
                grid,
                ph=float(process.ph),
                ph_final=(
                    float(process.ph_final)
                    if process.ph_final is not None
                    else None
                ),
                # B2.2: the sulfur lane runs the DYNAMIC pH state by default,
                # on the FROZEN calibration. A caller gets the old clamped
                # behaviour only by asking for it explicitly with
                # BufferSpec(kind="clamped") -- never by omission.
                buffer_spec=(
                    process.buffer if process.buffer is not None
                    else DEFAULT_BUFFER
                ),
                ph_drift=(
                    process.ph_drift if process.ph_drift is not None
                    else (ph_drift if ph_drift is not None else core_ph_drift())
                ),
                rtol=1e-8,
                atol=1e-14,
            )
            keys = list(SULFUR_INDEX)
        elif lane == ACRYLAMIDE:
            # B14 (2026-09-07): the declared flat a_w term inside the measured window scales the
            # acrylamide-forming step before integration; exactly 1.0 at the centre and outside.
            d = acrylamide_draw or {}
            acr_parameters, condition_terms = acrylamide_conditions.apply(
                parameters, process,
                aw_scale=d.get("aw_scale"),
                aw_elimination_deficit_scale=1.0 if d.get("aw_elimination_scale") is None else float(d["aw_elimination_scale"]),
                ph_exponent_formation=(acrylamide_conditions.PH_EXPONENT_FORMATION
                                       if d.get("ph_exponent_formation") is None else float(d["ph_exponent_formation"])),
                ph_exponent_elimination=(acrylamide_conditions.PH_EXPONENT_ELIMINATION
                                         if d.get("ph_exponent_elimination") is None else float(d["ph_exponent_elimination"])),
            )
            if condition_terms:
                metadata.setdefault("condition_terms", list(condition_terms))
            run = integrate_acrylamide(
                acr_parameters,
                float(temperature_c) + CELSIUS,
                state,
                grid,
                water_activity=process.water_activity,
                rtol=1e-8,
                atol=1e-14,
            )
            keys = list(ACRYLAMIDE_INDEX)
        else:
            # B12 (2026-09-07): the trunk's declared water-activity and pH terms scale the
            # named steps' k_ref before integration; exactly 1.0 at the references.
            aw_scale, ph_exponent = (trunk_draw or (None, None))
            trunk_parameters, condition_terms = trunk_conditions.apply(
                parameters, process,
                aw_scale=1.0 if aw_scale is None else float(aw_scale),
                ph_exponent=(trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT
                             if ph_exponent is None else float(ph_exponent)),
            )
            metadata.setdefault("condition_terms", list(condition_terms))
            run = integrate(
                trunk_parameters,
                float(temperature_c) + CELSIUS,
                state,
                grid,
                rtol=1e-8,
                atol=1e-14,
            )
            keys = list(SPECIES_KEYS)

        state = {key: float(run.concentrations[-1, i]) for i, key in enumerate(keys)}
        metadata["segments"].append(
            {
                "index": index,
                "duration_min": float(duration),
                "temperature_C": float(temperature_c),
                "extrapolation_warnings": list(
                    run.metadata.get("extrapolation_warnings", [])
                ),
                # B2.2: the pH is now an OUTPUT of the sulfur lane, not only an
                # input, so it travels with the segment that produced it.
                "ph_mode": run.metadata.get("ph_mode"),
                "ph_in_situ_start": run.metadata.get("ph_initial_in_situ"),
                "ph_in_situ_end": run.metadata.get("ph_final_in_situ"),
                "ph_cooled_end": run.metadata.get("ph_final_cooled"),
                "ph_notes": list(run.metadata.get("ph_notes", [])),
            }
        )
    return state, metadata


def predict(
    spec: FormulationSpec,
    targets: Sequence[str],
    *,
    parameters: Optional[Mapping[str, Any]] = None,
    draw: Optional[CoreDraw] = None,
    size_declared_bands: bool = True,
) -> CorePrediction:
    """
    THE ENTRY POINT. Map ``spec`` onto a lane, integrate, emit B4 objects.

    An out-of-envelope request returns a prediction carrying the declaration
    and NO concentrations. It does not raise here -- a caller scoring a panel
    needs to record the refusal alongside the answers -- but every accessor
    that would hand back a number raises instead.

    B2. ``draw`` moves the sampled quantities (see :class:`CoreDraw`);
    ``size_declared_bands=False`` skips the furanone-corner and lipid lo/hi
    re-integrations, so a Monte-Carlo caller that samples those bands does
    not ALSO price them by re-integration. ``parameters`` is still the raw
    operative override and cannot be combined with ``draw.maillard``. With
    the defaults the output is byte-identical to the pre-B2 engine.
    """
    if parameters is not None and draw is not None and draw.maillard is not None:
        raise ValueError(
            "pass either an operative `parameters` override or a fit-report-"
            "space `draw.maillard`, not both: they would silently shadow each "
            "other."
        )
    declaration = declare_envelope(spec, targets)
    if not declaration.is_answerable:
        return CorePrediction(spec=spec, declaration=declaration)

    lanes = declaration.lanes or ((declaration.lane or TRUNK),)
    maillard_lane = next((l for l in lanes if l in MAILLARD_LANES), None)

    final_state: Dict[str, float] = {}
    metadata: Dict[str, Any] = {"segments": [], "lane": declaration.lane,
                                "lanes": list(lanes)}
    if maillard_lane is not None:
        operative = (
            dict(parameters) if parameters is not None
            else core_parameters(
                maillard_lane,
                frozen=draw.maillard if draw is not None else None,
            )
        )
        if draw is not None and draw.furanone_partition_ea_kj_mol is not None:
            operative = _furanone_corner_parameters(
                operative, float(draw.furanone_partition_ea_kj_mol)
            )
        trunk_draw = (
            (draw.trunk_aw_scale, draw.trunk_ph_exponent) if draw is not None else None
        )
        acrylamide_draw = (
            {"aw_scale": draw.acrylamide_aw_scale, "aw_elimination_scale": draw.acrylamide_aw_elimination_scale,
             "ph_exponent_formation": draw.acrylamide_ph_exponent_formation,
             "ph_exponent_elimination": draw.acrylamide_ph_exponent_elimination}
            if draw is not None else None
        )
        final_state, metadata = _integrate_program(
            maillard_lane, operative, dict(declaration.mapped_precursors), spec.process,
            ph_drift=draw.ph_drift if draw is not None else None,
            trunk_draw=trunk_draw,
            acrylamide_draw=acrylamide_draw,
            reservoir_scale=draw.oxygen_reservoir_scale if draw is not None else None,
        )
        metadata["lanes"] = list(lanes)

    # B7. PRICE THE FURANONE ASSUMPTION BY RE-INTEGRATION, not by nominating a
    # width. Two extra integrations at the corners of the declared +/-50 kJ/mol
    # partition barrier. The result is often much narrower than the barrier
    # alone would suggest, because the deoxyosone POOL that feeds the edge is
    # itself depleting -- which a nominal width could not have shown.
    furanic_decades: Dict[str, float] = {}
    if size_declared_bands and maillard_lane is not None and (
        set(declaration.mapped_targets.values()) & set(FURANONE_BANDED_KEYS)
    ):
        from .parameters_furanic import FURANONE_PARTITION_EA_BAND_KJ_MOL

        band = float(FURANONE_PARTITION_EA_BAND_KJ_MOL)
        corners = []
        for offset in (-band, +band):
            corner_state, _ = _integrate_program(
                maillard_lane,
                _furanone_corner_parameters(operative, offset),
                dict(declaration.mapped_precursors),
                spec.process,
                trunk_draw=trunk_draw,
                acrylamide_draw=acrylamide_draw,
            )
            corners.append(corner_state)
        for key in FURANONE_BANDED_KEYS:
            lo = min(float(c.get(key, 0.0)) for c in corners)
            hi = max(float(c.get(key, 0.0)) for c in corners)
            if lo > 0.0 and hi > 0.0 and float(final_state.get(key, 0.0)) > 0.0:
                furanic_decades[key] = 0.5 * abs(math.log10(hi / lo))
        metadata["furanone_partition_band_kj_mol"] = band
        metadata["furanone_partition_corner_mmol_per_l"] = [
            {k: float(c.get(k, 0.0)) for k in FURANONE_BANDED_KEYS} for c in corners
        ]

    extra_decades: Dict[str, float] = {}
    if LIPID in lanes:
        lipid_state, lipid_metadata, extra_decades = _run_lipid_lane(
            spec, declaration,
            q10=draw.q10 if draw is not None else None,
            lipid_scale=draw.lipid_fraction_scale if draw is not None else None,
            pv_scale=draw.peroxide_scale if draw is not None else None,
            corners=size_declared_bands,
        )
        overlap = set(lipid_state) & set(final_state)
        if overlap:
            raise AssertionError(
                "lipid and Maillard states overlap on "
                f"{sorted(overlap)} -- the direct-sum co-integration ruling "
                "assumed disjoint species sets and that assumption has broken."
            )
        final_state.update(lipid_state)
        metadata["lipid"] = lipid_metadata

    concentrations: Dict[str, float] = {}
    for compound, key in declaration.mapped_targets.items():
        mmol = float(final_state.get(key, 0.0))
        if _TARGET_LANE.get(key) == LIPID:
            from .species_lipid import (
                mmol_per_litre_to_ug_per_litre as _lipid_ug,
            )

            concentrations[compound] = _lipid_ug(key, mmol)
            continue
        if key == "ACR":
            concentrations[compound] = acrylamide_ppb(mmol)
        elif key in MOLECULAR_WEIGHT_G_PER_MOL:
            concentrations[compound] = mmol_per_litre_to_ug_per_litre(key, mmol)
        elif key in _REPORTED_IN_MMOL_PER_L:
            # The elemental and lumped pools have no molecular weight because they are not
            # molecules; they are reported in their own unit rather than given an invented one.
            concentrations[compound] = mmol
        else:
            # B34 (2026-09-11). THIS USED TO BE THE `else`, AND THAT COST TWO WAVES. A species with
            # no molar mass was silently reported in mmol/L, which reads as a prediction between
            # one and five ORDERS too small: 2-pentylfuran in B28 (diagnosed for a day as a routing
            # problem), then 3-deoxyglucosone and methylglyoxal the moment B34 asked for them.
            # A missing weight is now a bug report, not a quiet change of unit.
            raise KeyError(
                f"{compound!r} resolves to species {key!r}, which has no molecular weight and is "
                "not a declared unit-less pool. Reporting it in mmol/L would look like a prediction "
                f"{'a factor of its molar mass'} too small. Add it to "
                "species_sulfur.MOLECULAR_WEIGHT_G_PER_MOL, or to engine._REPORTED_IN_MMOL_PER_L "
                "if it is genuinely not a molecule."
            )

    # 2026-09-08 (the matrix layer): declared binding of aldehydes and HMF to the charged protein
    # sites, applied after integration as a pseudo-first-order factor over the thermal programme,
    # its bracket priced as an interval width. Nothing happens without a stated protein loading.
    from .matrix_sites import bound_fraction as _bound_fraction
    from .matrix_sites import resolve as _resolve_sites

    # B31: what the pot started with. Added BEFORE the binding factor below, because the protein
    # binds a carried molecule and a made one alike. Zero unless the bundle declares it, and the
    # declaration may only quote a printed unheated control of the same pot.
    carried = dict(getattr(spec.process, "carried_volatiles", None) or {})
    carried_applied: Dict[str, float] = {}
    if carried:
        # 2026-09-11 (review of PR #16): matched THROUGH THE ALIAS TABLE, not by raw string. A
        # level declared as '2-pentylfuran' and a target requested as '2-pentyl furan' are the
        # same species and used to miss each other silently. And a declared 0.0 is a declaration
        # ("not detected" in the unheated control), not an absence of one.
        by_species = _carried_by_species(carried)
        for compound, key in declaration.mapped_targets.items():
            if compound not in concentrations or key not in by_species:
                continue
            amount = by_species[key]
            concentrations[compound] = float(concentrations[compound]) + amount
            carried_applied[compound] = amount
    charged_sites, sites_note = _resolve_sites(spec.process)
    binding: Dict[str, Any] = {}
    if charged_sites is not None:
        for compound, key in declaration.mapped_targets.items():
            result = _bound_fraction(key, charged_sites, spec.process.thermal.segments)
            if result is None or compound not in concentrations:
                continue
            concentrations[compound] = concentrations[compound] * result["remaining_fraction"]
            binding[compound] = result
    metadata["matrix_sites"] = charged_sites.as_dict() if charged_sites is not None else None
    if sites_note:
        metadata["matrix_sites_note"] = sites_note
    metadata["matrix_binding"] = binding
    metadata["carried_volatiles"] = carried_applied
    metadata["matrix_extra_decades"] = {c: r["extra_decades"] for c, r in binding.items()}
    metadata["ph"] = float(spec.process.ph)
    metadata["ph_final"] = spec.process.ph_final
    metadata["buffer"] = (
        spec.process.buffer.as_dict() if spec.process.buffer is not None
        else DEFAULT_BUFFER.as_dict()
    )
    metadata["ph_drift_constants"] = (
        spec.process.ph_drift.as_dict() if spec.process.ph_drift is not None
        else (core_ph_drift().as_dict() if maillard_lane == SULFUR else None)
    )
    metadata["matrix"] = spec.process.matrix
    metadata["thermal_program"] = spec.process.thermal.describe()
    if draw is not None and not draw.is_centre:
        metadata["draw"] = {
            "maillard_override_blocks": sorted(draw.maillard or {}),
            "q10": draw.q10,
            "lipid_fraction_scale": draw.lipid_fraction_scale,
            "peroxide_scale": draw.peroxide_scale,
            "furanone_partition_ea_kj_mol": draw.furanone_partition_ea_kj_mol,
            "ph_drift": draw.ph_drift.as_dict() if draw.ph_drift is not None else None,
        }
    if not size_declared_bands:
        metadata["declared_bands_sized"] = False
    # B6: the declared-assumption band, re-keyed from species key to the
    # caller's own compound name so ``absolutes()`` can find it.
    if extra_decades:
        metadata["lipid_extra_decades"] = {
            compound: extra_decades[key]
            for compound, key in declaration.mapped_targets.items()
            if key in extra_decades
        }
    if furanic_decades:
        metadata["furanic_extra_decades"] = {
            compound: furanic_decades[key]
            for compound, key in declaration.mapped_targets.items()
            if key in furanic_decades
        }

    if declaration.unreachable_targets:
        metadata["refused_targets"] = {
            name: "no chain of the lane's reactions leads from the charge to this species; "
                  "its integrated value is exactly zero by construction and is not reported"
            for name in declaration.unreachable_targets
        }
        for name in declaration.unreachable_targets:
            concentrations.pop(name, None)

    return CorePrediction(
        spec=spec,
        declaration=declaration,
        concentrations_ug_per_l=concentrations,
        species_mmol_per_l=final_state,
        run_metadata=metadata,
    )


# ---------------------------------------------------------------------------
# The comparative surface -- the layer's PRIMARY output
# ---------------------------------------------------------------------------


#: Lanes whose parameters carry NO pH term (declared in their parameter modules).
#: B12 (2026-09-07): the trunk gained a declared pH term on its Amadori-decay steps.
NO_PH_TERM_LANES = frozenset({LIPID})
#: B15: lanes whose pH term is measured only inside a window (outside it the factor is held and a
#: comparison that leaves the window is refused).
PH_TERM_WINDOWS = {ACRYLAMIDE: acrylamide_conditions.PH_WINDOW}
#: Lanes that carry a water-activity term (B12: the trunk's declared multiplier).
AW_TERM_LANES = frozenset({TRUNK, ACRYLAMIDE})
#: B14: lanes whose a_w term exists only inside a measured window (outside it the axis is refused).
AW_TERM_WINDOWS = {ACRYLAMIDE: acrylamide_conditions.AW_WINDOW}


def _lanes_of(declaration) -> Tuple[str, ...]:
    lanes = getattr(declaration, "lanes", None)
    if lanes:
        return tuple(str(x) for x in lanes)
    return (str(declaration.lane),) if declaration.lane else ()


def axis_refusal(spec_a, spec_b, declaration_a, declaration_b) -> Optional[str]:
    """
    2026-09-03 (owner decision, step 5): a comparison that moves an axis the resolved
    lane carries no term for is REFUSED, not answered with two identical numbers.

    * water activity differs between the arms: NO lane carries an a_w term.
    * pH differs and every resolved lane is trunk / acrylamide / lipid: those lanes are
      homogeneous in pH by declaration; only the sulfur lane carries a pH trajectory.

    Returns the refusal reason, or None when the comparison is answerable.
    """
    pa, pb = spec_a.process, spec_b.process
    aw_a, aw_b = pa.water_activity, pb.water_activity
    lanes = set(_lanes_of(declaration_a)) | set(_lanes_of(declaration_b))
    if aw_a is not None and aw_b is not None and abs(float(aw_a) - float(aw_b)) > 1e-9:
        if not (lanes & AW_TERM_LANES):
            return (
                f"REFUSED -- the two arms differ in WATER ACTIVITY and the resolved lane(s) "
                f"({', '.join(sorted(lanes)) or 'none'}) carry no a_w term; the model would return "
                "identical arms and call it a comparison. The trunk lane carries a declared a_w term "
                "(B12) and the acrylamide lane a declared flat one inside a_w 0.88-0.99 (B14). Hold "
                "a_w fixed, or bring a measurement."
            )
        for lane_name, (lo, hi) in AW_TERM_WINDOWS.items():
            if lane_name in lanes and not all(lo - 1e-9 <= float(a) <= hi + 1e-9 for a in (aw_a, aw_b)):
                return (
                    f"REFUSED -- the two arms differ in WATER ACTIVITY ({float(aw_a):.2f} vs {float(aw_b):.2f}) "
                    f"and the {lane_name} lane's a_w term is measured only inside {lo:.2f}-{hi:.2f} "
                    f"(De Vleeschouwer 2008, B14): outside it the lane carries no term and would return "
                    "identical arms. Keep both arms inside the window, or bring a measurement."
                )
    if abs(float(pa.ph) - float(pb.ph)) > 1e-9:
        if lanes and lanes <= NO_PH_TERM_LANES:
            return (
                f"REFUSED -- the two arms differ in pH and the resolved lane(s) "
                f"({', '.join(sorted(lanes))}) carry NO pH term by declaration; the model would "
                "return identical arms. The sulfur lane carries a pH trajectory, the trunk a "
                "declared Amadori-decay pH term (B12) and the acrylamide lane a declared initial-pH "
                "factor inside pH 4-8 (B15)."
            )
        for lane_name, (lo, hi) in PH_TERM_WINDOWS.items():
            if lane_name in lanes and not all(lo - 1e-9 <= float(v) <= hi + 1e-9 for v in (pa.ph, pb.ph)):
                return (
                    f"REFUSED -- the two arms differ in pH ({float(pa.ph):g} vs {float(pb.ph):g}) and the "
                    f"{lane_name} lane's pH factor is measured only inside pH {lo:g}-{hi:g} (De Vleeschouwer "
                    "2006, B15): outside it the factor is held at the window edge and the arms would not be a "
                    "comparison. Keep both arms inside the window, or bring a measurement."
                )
    return None


def compare(
    spec_a: FormulationSpec,
    spec_b: FormulationSpec,
    targets: Sequence[str],
    *,
    predict_fn: Optional[Callable[..., "CorePrediction"]] = None,
) -> Dict[str, Any]:
    """
    Per-compound RATIOS between two formulations, via the B4 layer.

    A ratio is the layer's primary unit because the two dominant error sources
    -- the HS-SPME calibration offset and the air/water partition constant --
    are shared between the arms and cancel exactly in a within-run ratio. If
    EITHER arm is out of envelope, no ratio is emitted for the affected
    compounds; a ratio against a refusal is not a ratio.

    Q1: each arm now also carries its OWN OAV table and its own interval rows
    (``oav_table_a``/``oav_table_b``, ``rows_a``/``rows_b``). Before this, a
    compare returned the arms as plain ``as_dict()`` payloads, which drop the
    object and therefore drop ``.oav()`` and ``.absolutes()`` -- so the report
    layer rebuilt the OAV table by hand from the run dict. That copy had
    ALREADY drifted: it was written in B6 and never taught about B7's furanone
    declared-assumption band, so a compare drew narrower intervals than a
    predict of the identical arm. The tables are emitted here, from the live
    objects, so there is exactly one implementation to keep correct.
    """
    # 2026-09-08: a caller may supply the predictor (a per-laboratory calibration wraps `predict`);
    # the default is byte-identical to the plain engine.
    _predict = predict_fn or predict
    run_a = _predict(spec_a, targets)
    run_b = _predict(spec_b, targets)

    if not (run_a.answered and run_b.answered):
        return {
            "comparable": False,
            "declaration_a": run_a.declaration.as_dict(),
            "declaration_b": run_b.declaration.as_dict(),
            "reason": (
                "at least one arm is out of envelope; a ratio against a "
                "refusal is not a ratio."
            ),
        }
    refusal = axis_refusal(spec_a, spec_b, run_a.declaration, run_b.declaration)
    if refusal is not None:
        return {
            "comparable": False,
            "declaration_a": run_a.declaration.as_dict(),
            "declaration_b": run_b.declaration.as_dict(),
            "reason": refusal,
            "axis_refusal": True,
        }

    shared = sorted(
        set(run_a.concentrations_ug_per_l) & set(run_b.concentrations_ug_per_l)
    )
    payload = compare_formulations(
        {c: run_a.concentrations_ug_per_l[c] for c in shared},
        {c: run_b.concentrations_ug_per_l[c] for c in shared},
        label_a=spec_a.name,
        label_b=spec_b.name,
    )
    # 2026-09-04: a row whose formation in either arm runs through an UNIDENTIFIED route
    # (declared on the arm, see unidentified_routes) is not a ratio between two predictions
    # but between a prediction and a band-floor artefact. It is reported as undefined, with the
    # arm named, so that "1e13x higher in the pentose arm" never reaches a table.
    for row in payload.get("rows", []):
        arms = [label for label, decl in (("a", run_a.declaration), ("b", run_b.declaration))
                if declared_unidentified(decl, str(row["compound"]))]
        if arms:
            row["ratio_a_over_b"] = None
            row["direction"] = "undefined"
            row["within_reliability_band"] = False
            row["unidentified_arm"] = arms[0] if len(arms) == 1 else "both"
            row["note"] = (
                f"arm {row['unidentified_arm'].upper()}: {HEXOSE_ENTRY_UNIDENTIFIED} -- its number is a "
                "band-floor artefact, so no ratio is claimed; the model supports the ordering "
                "'pentose above hexose' structurally, not a magnitude"
            )
    rows_ = payload.get("rows", [])
    payload["n_undefined"] = sum(1 for r in rows_ if r.get("direction") == "undefined")
    payload["n_resolved"] = sum(1 for r in rows_ if r.get("direction") != "undefined" and not r.get("within_reliability_band"))
    return {
        "comparable": True,
        "ratios": payload,
        "declaration_a": run_a.declaration.as_dict(),
        "declaration_b": run_b.declaration.as_dict(),
        "run_a": run_a.as_dict(),
        "run_b": run_b.as_dict(),
        # Q1: the arms' OWN B4 output, from the live objects. See the docstring.
        "oav_table_a": dict(run_a.oav()),
        "oav_table_b": dict(run_b.oav()),
        "rows_a": [dict(r) for r in run_a.interval_rows()],
        "rows_b": [dict(r) for r in run_b.interval_rows()],
    }


def residual_report(
    compound: str, matrix: str, measured_ratio: float, ph: Optional[float] = None
) -> Dict[str, Any]:
    """
    The B4 residual decomposition for one compound in one matrix, surfaced.

    "Measured shift Nx, the model's named terms explain Mx, the rest is
    unexplained residual" -- which is the layer's honest output on a matrix it
    has no constant for.
    """
    prediction = predict_matrix_shift(compound, matrix, ph=ph)
    decomposition = decompose_residual(prediction, float(measured_ratio))
    return {
        "compound": compound,
        "matrix": matrix,
        "measured_ratio": float(measured_ratio),
        "predicted_shift": getattr(prediction, "ratio", None),
        "model_state": getattr(prediction, "state", None),
        "decomposition": decomposition,
    }


def fit_report_paths() -> Tuple[Path, ...]:
    """The frozen fit reports this engine reads, in lane order (B1, sulfur wave, B3, B6, B7).

    The scorecard, the envelope and the directional artifact list exactly these as their
    provenance inputs; a report on disk the engine does not read is not a parameter source.
    """
    return (_B1_FIT_REPORT, _B2_FIT_REPORT, _B3_FIT_REPORT, _B6_FIT_REPORT, _B7_FIT_REPORT)


def engine_metadata() -> Dict[str, Any]:
    """
    What this engine is, for embedding in every artifact it produces.

    Q1: THE STEP COUNTS ARE NOW COUNTED, NOT TRANSCRIBED. They were written out
    as literals at B5 -- "15 steps", "79 steps", "31 steps" -- and B6 and B7
    then added edges without updating them, so every artifact this engine has
    produced since carried step counts that were wrong by 11, 14 and 11
    respectively, and a ``"wave": "B5"`` stamp two waves out of date. A
    provenance field that silently describes an older model than the one that
    produced the number is worse than no field, because it is quoted in
    good faith. Counting them at call time makes the class of error impossible.
    """
    from .acrylamide import FULL_ACRYLAMIDE_REACTIONS
    from .network import REACTIONS
    from .sulfur import FULL_REACTIONS

    return {
        "module": "src/kinetic_core/engine.py",
        "wave": "furanic channels (HMF, DMHF), fit wave B7",
        "lanes": list(LANES),
        "lane_networks": {
            TRUNK: f"REACTIONS ({len(REACTIONS)} steps), no pH term, no a_w term",
            SULFUR: (
                f"FULL_REACTIONS ({len(FULL_REACTIONS)} steps) = trunk + sulfur, "
                "pH trajectory"
            ),
            ACRYLAMIDE: (
                f"FULL_ACRYLAMIDE_REACTIONS ({len(FULL_ACRYLAMIDE_REACTIONS)} "
                "steps) = trunk + acrylamide; sulfur STEPS deliberately absent"
            ),
            LIPID: (
                "a hydroperoxide pool resolved by position (9-/13-) and "
                "geometry (cis,trans / trans,trans), decomposing first-order "
                "into Frankel 1989's six-product measured slate. The "
                "DISTRIBUTION is fitted and frozen; the RATE is a declared, "
                "bounded ASSUMPTION and every prediction says so."
            ),
        },
        "lanes_compose": False,
        "lipid_lane_cointegrates": {
            "rule": "direct sum with any ONE Maillard lane",
            "why": (
                "disjoint species sets, and the only candidate coupling (the "
                "aldehyde-lysine covalent channel) is INERT BY RULING "
                "(docs/reference/FIT_HOLDOUT_DECLARATION.md, amendment 6, ruling 2). Checked at "
                "every call by lipid.lane_coupling_verdict, not hard-coded."
            ),
            "condition": (
                "revisit the moment the aldehyde-lysine Ea on food proteins is "
                "measured -- the amine pool then becomes genuinely shared"
            ),
        },
        "lipid_rate_is_an_assumption": True,
        "parameters_from": [
            data_paths.rel(_B1_FIT_REPORT),
            data_paths.rel(_B2_FIT_REPORT),
            data_paths.rel(_B3_FIT_REPORT),
        ],
        "fits_anything": False,
        "network_ph": NETWORK_PH,
        "unrepresented_compounds": sorted(set(UNREPRESENTED_COMPOUNDS)),
    }


__all__ = [
    "ACRYLAMIDE",
    "LIPID",
    "LIPID_CARRIER_ALIASES",
    "MAILLARD_LANES",
    "core_lipid_model",
    "resolve_lanes",
    "CoreDraw",
    "CorePrediction",
    "EnvelopeDeclaration",
    "FormulationSpec",
    "LANES",
    "LANE_DEFAULT_TARGETS",
    "default_targets_for",
    "OutOfEnvelope",
    "PRECURSOR_ALIASES",
    "ProcessSpec",
    "SULFUR",
    "TARGET_ALIASES",
    "axis_refusal",
    "AW_TERM_LANES",
    "AW_TERM_WINDOWS",
    "PH_TERM_WINDOWS",
    "NO_PH_TERM_LANES",
    "TRUNK",
    "ThermalProgram",
    "UNREPRESENTED_COMPOUNDS",
    "b1_fitted",
    "compare",
    "core_parameters",
    "frozen_parameters",
    "declare_envelope",
    "engine_metadata",
    "predict",
    "residual_report",
    "resolve_lane",
    "resolve_matrix",
]
