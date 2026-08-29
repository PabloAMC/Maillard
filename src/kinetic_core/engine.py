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
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

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
from .parameters_sulfur import MEASURED_SULFUR, with_fitted_sulfur
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

_ROOT = Path(__file__).resolve().parents[2]
_B1_FIT_REPORT = _ROOT / "results/validation/kinetic_core_b1_fit_report.json"
#: THE SULFUR LANE'S FROZEN PARAMETERS -- and this is a CUTOVER, stated here
#: rather than buried in a wave report. Build Wave B2.3 refits B2.2's own 48
#: parameters on B2.2's own 58 FIT rows after a CONSERVATION FIX (see
#: `ph_state.validate_charge_closure` and `sulfur.CENTRE_LEDGER`), so where a
#: B2.3 report exists it SUPERSEDES B2.2's: the B2.2 numbers were fitted
#: against a network that manufactured strong base out of bookkeeping, and
#: preferring them would be preferring a known defect. B2.2's report is kept
#: as the fallback so that every earlier artefact stays regenerable on a
#: checkout that has not run B2.3.
_B2_FIT_REPORT_CANDIDATES = (
    _ROOT / "results/validation/kinetic_core_b2_3_fit_report.json",
    _ROOT / "results/validation/kinetic_core_b2_2_fit_report.json",
)
_B2_FIT_REPORT = next(
    (p for p in _B2_FIT_REPORT_CANDIDATES if p.exists()),
    _B2_FIT_REPORT_CANDIDATES[-1],
)
_B3_FIT_REPORT = _ROOT / "results/validation/kinetic_core_b3_fit_report.json"


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
    "alanine": "Ala",
    "l-alanine": "Ala",
    "methylglyoxal": "MGO",
    "norfuraneol": "NF",
    "amadori": "AMA",
    "arp": "ARP",
}

#: Target-compound synonyms -> core species key.
TARGET_ALIASES: Mapping[str, str] = {
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
        "(FIT_HOLDOUT_DECLARATION Amendment 12 corrected Amendment 8 on "
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
        "The B6 lipid lane exists and forms the SIX products Frankel 1989 "
        "measured, but 1-hexanol is not one of them and NO aldehyde-reduction "
        "step is measured anywhere in the corpus -- in a thermally processed "
        "extrudate the reductant pool is not even identified. The FAST lane "
        "emits a number for it; this lane refuses. See "
        "parameters_lipid.PROHIBITED_DERIVATIONS."
    ),
    "2-pentylfuran": (
        "The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's "
        "six-product slate and no branch fraction for the linoleate -> "
        "alkylfuran route is measured anywhere in the fit corpus. The FAST "
        "lane's shipped 0.08 has no source. Refused rather than invented."
    ),
    "2-pentyl furan": (
        "The B6 lipid lane exists, but 2-pentylfuran is NOT in Frankel 1989's "
        "six-product slate and no branch fraction for the linoleate -> "
        "alkylfuran route is measured anywhere in the fit corpus. The FAST "
        "lane's shipped 0.08 has no source. Refused rather than invented."
    ),
    "propanal": (
        "The B6 lipid lane forms no propanal. Propanal is an alpha-LINOLENATE "
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
    # -- B6, the lipid lane ------------------------------------------------
    "HEXANAL": LIPID,
    "NONANAL": LIPID,
    "PENTANE": LIPID,
    "DECADIENAL": LIPID,
    "ME_OCTANOATE": LIPID,
    "ME_9_OXONONANOATE": LIPID,
    "ME_13_OXO_TRIDECADIENOATE": LIPID,
}

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
    return out


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

    # --- lane ------------------------------------------------------------
    lanes, lane_reasons = resolve_lanes(
        list(mapped_precursors), list(mapped_targets.values()), carriers
    )
    lane = lanes[0] if lanes else None
    reasons.extend(lane_reasons)

    # --- the lipid lane's own refusals ------------------------------------
    if LIPID in lanes:
        from .parameters_lipid import LIPID_CARRIERS, oleate_fraction

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
                reasons.append(
                    "UNREPRESENTED TARGETS: nonanal -- the lipid lane exists "
                    "and nonanal is a species in it, but its ONLY parent is the "
                    "OLEATE hydroperoxide pool and the oleate -> nonanal branch "
                    f"fraction is measured NOWHERE in the fit corpus. This "
                    f"matrix is {100.0 * oleate:.0f} % oleate by fatty-acid "
                    "share, so the pool is not zero. Frankel 1989 fed linoleate "
                    "only and nonanal appears in no table, figure or sentence "
                    "of it -- that ABSENCE is a declared hold-out, and honouring "
                    "it means refusing here rather than carrying the FAST "
                    "lane's unsourced 'nonanal 0.15' forward."
                )
    if lipid_targets and LIPID not in lanes and not lane_reasons:
        reasons.append(
            "a lipid product was requested but the lipid lane was not selected"
        )

    # A target whose lane needs a precursor species this charge cannot supply.
    if lane is not None and not unmapped:
        if lane == SULFUR and not (
            {"Cys", "THI", "PENT", "ARP", "H2S"} & set(mapped_precursors)
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
    if lane in (TRUNK, ACRYLAMIDE) and abs(float(spec.process.ph) - NETWORK_PH) > 1e-9:
        warnings.append(
            f"pH {spec.process.ph:g} was supplied, but the {lane} lane carries "
            f"NO pH term at all -- its parameters are homogeneous at pH "
            f"{NETWORK_PH:g}. The pH is recorded and IGNORED; it changes no "
            f"rate. Any pH sensitivity in the measurement is unmodelled."
        )
    if spec.process.water_activity is not None and lane == ACRYLAMIDE:
        warnings.append(
            "water activity is METADATA ONLY on the acrylamide lane: it changes "
            "no rate. The corpus spans a_w 0.35-1.0 without measuring the axis."
        )
    if spec.process.ph_final is not None and lane != SULFUR:
        warnings.append(
            f"a final pH was supplied but the {lane} lane has no pH trajectory; "
            f"it is ignored."
        )
    # B2.2: the buffer is an input with a declared default, and its ABSENCE is
    # an extrapolation rather than a silent assumption.
    if lane == SULFUR:
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
                "physical state are extrapolations. K5a sec. 6.2: that limb's "
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
                "activation energy. K5a declared gap G2: the 50-150 C window "
                "is empty. EXPECT HMF TO BE OVER-PREDICTED."
            )
            if peak_k > HMF_SINK_NO_EXTRAPOLATION_ABOVE_K and (
                "Cys" in mapped_precursors
            ):
                warnings.append(
                    "5-HMF + cysteine: the sink constant is HELD at its 50 C "
                    "value for this whole program. Holding it UNDER-states the "
                    "sink; extrapolating it is a named prohibited derivation "
                    "(K5a sec. 7.3), and the direction is stated rather than "
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
    )


# ---------------------------------------------------------------------------
# Frozen parameters
# ---------------------------------------------------------------------------


def _read(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise SystemExit(
            f"{path} not found. The engine never fits anything; it reads the "
            f"frozen fit reports. Regenerate them first."
        )
    return json.loads(path.read_text())


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


def core_parameters(lane: str) -> Dict[str, Any]:
    """The full operative parameter set for one lane, from the frozen reports."""
    if lane == TRUNK:
        return dict(operative_parameters(b1_fitted()))
    if lane == SULFUR:
        frozen = _read(_B2_FIT_REPORT)["frozen_parameters"]
        parameters = dict(operative_parameters(b1_fitted()))
        parameters.update(MEASURED_SULFUR)
        parameters.update(
            with_fitted_sulfur(
                frozen["log10_k_ref_at_145C"],
                frozen["lumped_formation_Ea_kJ_mol"],
                frozen["decay_Ea_kJ_mol"],
            )
        )
        return parameters
    if lane == ACRYLAMIDE:
        frozen = _read(_B3_FIT_REPORT)["frozen_parameters"]
        parameters = dict(operative_parameters(b1_fitted()))
        parameters.update(MEASURED_ACRYLAMIDE)
        parameters.update(
            with_fitted_acrylamide(
                frozen["log10_k_ref_at_160C"], frozen["fitted_Ea_kJ_mol"]
            )
        )
        return parameters
    if lane == LIPID:
        raise ValueError(
            "the lipid lane has no mass-action parameter dictionary: its "
            "frozen state is a BRANCH MODEL plus a rate ASSUMPTION. Call "
            "core_lipid_model() instead -- the distinction is the module's "
            "whole point."
        )
    raise ValueError(f"unknown lane {lane!r}")


_B6_FIT_REPORT = _ROOT / "results/validation/kinetic_core_b6_fit_report.json"


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
    spec: FormulationSpec, declaration: EnvelopeDeclaration
) -> Tuple[Dict[str, float], Dict[str, Any], Dict[str, float]]:
    """
    Run the B6 lipid lane and size its interval BY RE-INTEGRATION.

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

    def _run(q10, lipid_scale, pv_scale):
        state: Dict[str, float] = {}
        runs = []
        for key in carrier_keys:
            carrier = LIPID_CARRIERS[key]
            charge = charge_from_carrier(
                carrier, composition,
                lipid_fraction=lipid_scale(carrier),
                peroxide_value_meq_per_kg=pv_scale(carrier),
            )
            run = integrate_lipid(charge, segments, branch, q10=q10)
            runs.append(run)
            for species_key, value in run.state_mmol_per_l.items():
                state[species_key] = state.get(species_key, 0.0) + value
        return state, runs

    point, point_runs = _run(
        None, lambda c: c.lipid_mass_fraction, lambda c: c.peroxide_value_meq_per_kg
    )
    low, _ = _run(Q10_ASSUMPTION.lo, lambda c: c.lipid_lo, lambda c: c.pv_lo)
    high, _ = _run(Q10_ASSUMPTION.hi, lambda c: c.lipid_hi, lambda c: c.pv_hi)

    extra_decades: Dict[str, float] = {}
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


def _integrate_program(
    lane: str,
    parameters: Mapping[str, Any],
    initial: Mapping[str, float],
    process: ProcessSpec,
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """
    Integrate a piecewise-constant thermal program, chaining the state across
    segments, and return the FINAL state as ``{species_key: mmol/L}``.
    """
    state: Dict[str, float] = dict(initial)
    metadata: Dict[str, Any] = {"segments": [], "lane": lane}

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
                    else core_ph_drift()
                ),
                rtol=1e-8,
                atol=1e-14,
            )
            keys = list(SULFUR_INDEX)
        elif lane == ACRYLAMIDE:
            run = integrate_acrylamide(
                parameters,
                float(temperature_c) + CELSIUS,
                state,
                grid,
                water_activity=process.water_activity,
                rtol=1e-8,
                atol=1e-14,
            )
            keys = list(ACRYLAMIDE_INDEX)
        else:
            run = integrate(
                parameters,
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
) -> CorePrediction:
    """
    THE ENTRY POINT. Map ``spec`` onto a lane, integrate, emit B4 objects.

    An out-of-envelope request returns a prediction carrying the declaration
    and NO concentrations. It does not raise here -- a caller scoring a panel
    needs to record the refusal alongside the answers -- but every accessor
    that would hand back a number raises instead.
    """
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
            else core_parameters(maillard_lane)
        )
        final_state, metadata = _integrate_program(
            maillard_lane, operative, dict(declaration.mapped_precursors), spec.process
        )
        metadata["lanes"] = list(lanes)

    # B7. PRICE THE FURANONE ASSUMPTION BY RE-INTEGRATION, not by nominating a
    # width. Two extra integrations at the corners of the declared +/-50 kJ/mol
    # partition barrier. The result is often much narrower than the barrier
    # alone would suggest, because the deoxyosone POOL that feeds the edge is
    # itself depleting -- which a nominal width could not have shown.
    furanic_decades: Dict[str, float] = {}
    if maillard_lane is not None and (
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
        lipid_state, lipid_metadata, extra_decades = _run_lipid_lane(spec, declaration)
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
        else:
            # No molecular weight is defined for the elemental melanoidin pool;
            # it is reported in its own unit rather than given an invented one.
            concentrations[compound] = mmol

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


def compare(
    spec_a: FormulationSpec,
    spec_b: FormulationSpec,
    targets: Sequence[str],
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
    run_a = predict(spec_a, targets)
    run_b = predict(spec_b, targets)

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

    shared = sorted(
        set(run_a.concentrations_ug_per_l) & set(run_b.concentrations_ug_per_l)
    )
    payload = compare_formulations(
        {c: run_a.concentrations_ug_per_l[c] for c in shared},
        {c: run_b.concentrations_ug_per_l[c] for c in shared},
        label_a=spec_a.name,
        label_b=spec_b.name,
    )
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
        "wave": "B7 -- furanic channels (HMF, DMHF); propagator cutover at B5",
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
                "B6: a hydroperoxide pool resolved by position (9-/13-) and "
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
                "(FIT_HOLDOUT_DECLARATION Amendment 6 ruling 2). Checked at "
                "every call by lipid.lane_coupling_verdict, not hard-coded."
            ),
            "condition": (
                "revisit the moment the aldehyde-lysine Ea on food proteins is "
                "measured -- the amine pool then becomes genuinely shared"
            ),
        },
        "lipid_rate_is_an_assumption": True,
        "parameters_from": [
            str(_B1_FIT_REPORT.relative_to(_ROOT)),
            str(_B2_FIT_REPORT.relative_to(_ROOT)),
            str(_B3_FIT_REPORT.relative_to(_ROOT)),
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
    "TRUNK",
    "ThermalProgram",
    "UNREPRESENTED_COMPOUNDS",
    "b1_fitted",
    "compare",
    "core_parameters",
    "declare_envelope",
    "engine_metadata",
    "predict",
    "residual_report",
    "resolve_lane",
    "resolve_matrix",
]
