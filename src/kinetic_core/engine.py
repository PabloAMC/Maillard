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

THE THREE LANES, AND WHY THEY DO NOT COMPOSE
--------------------------------------------
  * ``TRUNK``      -- ``REACTIONS`` (15 steps). Glc/Fru/Gly -> melanoidins.
                      No pH term, no a_w term.
  * ``SULFUR``     -- ``FULL_REACTIONS`` (79 steps) = trunk + sulfur. Adds the
                      pentoses, cysteine, thiamine, MFT/FFT/furfural. Carries a
                      pH trajectory.
  * ``ACRYLAMIDE`` -- ``FULL_ACRYLAMIDE_REACTIONS`` (31 steps) = trunk +
                      acrylamide. Adds asparagine and the acrylamide block.

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

LANES: Tuple[str, ...] = (TRUNK, SULFUR, ACRYLAMIDE)


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
    "methanethiol": "MESH",
    "2-acetylthiazole": "ACTZ",
    "norfuraneol": "NF",
    "hydrogen sulfide": "H2S",
    "melanoidins": "MEL_N",
}

#: Compounds a user may plausibly ask for that the core CANNOT NAME, each with
#: the reason. Being on this list is what turns a request into a declaration
#: instead of a KeyError, and the reason is what makes the declaration useful.
UNREPRESENTED_COMPOUNDS: Mapping[str, str] = {
    "5-hydroxymethylfurfural": (
        "5-HMF is not a species in any core lane. The hexose-dehydration route "
        "that forms it was never parameterised: no dataset in the fit corpus "
        "measures it."
    ),
    "5-hydroxymethylfurfural (hmf)": (
        "5-HMF is not a species in any core lane. The hexose-dehydration route "
        "that forms it was never parameterised: no dataset in the fit corpus "
        "measures it."
    ),
    "hmf": (
        "5-HMF is not a species in any core lane. The hexose-dehydration route "
        "that forms it was never parameterised: no dataset in the fit corpus "
        "measures it."
    ),
    "dmhf": (
        "2,5-dimethyl-4-hydroxy-3(2H)-furanone is not a core species. The core "
        "carries NORFURANEOL (4-hydroxy-5-methyl-3(2H)-furanone, key NF), which "
        "is a different compound on a different (pentose) route. Reporting NF "
        "as DMHF would be a species substitution the corpus does not license."
    ),
    "furaneol": (
        "2,5-dimethyl-4-hydroxy-3(2H)-furanone is not a core species; see DMHF."
    ),
    "hexanal": (
        "The kinetic core has NO lipid-oxidation path. Hexanal is a lipid "
        "hydroperoxide beta-scission product and no core lane forms it."
    ),
    "nonanal": (
        "The kinetic core has NO lipid-oxidation path. Nonanal is the C9 "
        "fragment of the OLEATE double bond and no core lane forms it."
    ),
    "1-hexanol": (
        "The kinetic core has NO lipid-oxidation path, and no alcohol-reduction "
        "step anywhere."
    ),
    "2-pentylfuran": (
        "The kinetic core has NO lipid-oxidation path. 2-pentylfuran is a "
        "linoleate-derived alkylfuran and no core lane forms it."
    ),
    "2-pentyl furan": (
        "The kinetic core has NO lipid-oxidation path. 2-pentylfuran is a "
        "linoleate-derived alkylfuran and no core lane forms it."
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


#: The compounds each lane can be asked to REPORT, in the display names the
#: engine's vocabulary maps. Used when a caller does not name targets.
LANE_DEFAULT_TARGETS: Mapping[str, Tuple[str, ...]] = {
    TRUNK: ("melanoidins",),
    SULFUR: (
        "2-methyl-3-furanthiol (MFT)",
        "2-furfurylthiol (FFT)",
        "bis(2-methyl-3-furyl) disulfide",
        "furfural",
        "2-acetylthiazole",
        "methanethiol",
    ),
    ACRYLAMIDE: ("acrylamide",),
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
    for name in precursors:
        key = PRECURSOR_ALIASES.get(_norm(name))
        if key is not None:
            keys.append(key)
    if not keys:
        return ()
    lane, reasons = resolve_lane(keys, [])
    if lane is None or reasons:
        return ()
    return LANE_DEFAULT_TARGETS.get(lane, ())


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

    @property
    def is_answerable(self) -> bool:
        return self.state in ("in_envelope", "in_envelope_extrapolated")

    def as_dict(self) -> Dict[str, Any]:
        return {
            "state": self.state,
            "lane": self.lane,
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "unmapped_precursors": list(self.unmapped_precursors),
            "unrepresented_targets": [
                {"compound": c, "reason": r} for c, r in self.unrepresented_targets
            ],
            "mapped_precursors": dict(self.mapped_precursors),
            "mapped_targets": dict(self.mapped_targets),
        }


def resolve_lane(
    precursor_keys: Sequence[str], target_keys: Sequence[str]
) -> Tuple[Optional[str], Tuple[str, ...]]:
    """
    Pick the one lane that can carry this request, or explain why none can.

    Returns ``(lane_or_None, reasons)``.
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

    if len(required) > 1:
        return None, (
            "LANE CONFLICT: this request needs both the "
            + " and ".join(sorted(required))
            + " lanes at once. They do not compose -- the acrylamide network "
            "deliberately omits every sulfur step (acrylamide.OUT_OF_SCOPE), "
            "because composing them would spend the same cysteine twice. No "
            "single integration can answer it.",
        )
    if not required:
        return TRUNK, ()
    return required.pop(), ()


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
    for name, value in spec.precursors.items():
        key = PRECURSOR_ALIASES.get(_norm(name))
        if key is None:
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

    # --- lane ------------------------------------------------------------
    lane, lane_reasons = resolve_lane(
        list(mapped_precursors), list(mapped_targets.values())
    )
    reasons.extend(lane_reasons)

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


    if reasons:
        return EnvelopeDeclaration(
            state="out_of_envelope",
            lane=lane,
            reasons=tuple(reasons),
            warnings=tuple(warnings),
            unmapped_precursors=tuple(sorted(unmapped)),
            unrepresented_targets=tuple(unrepresented),
            mapped_precursors=mapped_precursors,
            mapped_targets=mapped_targets,
        )
    return EnvelopeDeclaration(
        state="in_envelope_extrapolated" if warnings else "in_envelope",
        lane=lane,
        warnings=tuple(warnings),
        mapped_precursors=mapped_precursors,
        mapped_targets=mapped_targets,
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
    raise ValueError(f"unknown lane {lane!r}")


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
        """Every answered concentration wrapped in its B4 reliability band."""
        return {
            compound: absolute_concentration(
                value,
                via_partition=True,
                provenance=f"kinetic core {self.declaration.lane} lane",
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
        by_key = {
            self.declaration.mapped_targets.get(compound, compound): value
            for compound, value in self.concentrations_ug_per_l.items()
        }
        return oav_table(
            by_key,
            matrix=matrix or resolve_matrix(self.spec.process.matrix),
            temperature_c=self.spec.process.thermal.peak_temperature_c,
        )

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

    lane = declaration.lane or TRUNK
    operative = dict(parameters) if parameters is not None else core_parameters(lane)
    final_state, metadata = _integrate_program(
        lane, operative, dict(declaration.mapped_precursors), spec.process
    )

    concentrations: Dict[str, float] = {}
    for compound, key in declaration.mapped_targets.items():
        mmol = float(final_state.get(key, 0.0))
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
        else (core_ph_drift().as_dict() if lane == SULFUR else None)
    )
    metadata["matrix"] = spec.process.matrix
    metadata["thermal_program"] = spec.process.thermal.describe()

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
    """What this engine is, for embedding in every artifact it produces."""
    return {
        "module": "src/kinetic_core/engine.py",
        "wave": "B5 -- the propagator cutover",
        "lanes": list(LANES),
        "lane_networks": {
            TRUNK: "REACTIONS (15 steps), no pH term, no a_w term",
            SULFUR: "FULL_REACTIONS (79 steps) = trunk + sulfur, pH trajectory",
            ACRYLAMIDE: (
                "FULL_ACRYLAMIDE_REACTIONS (31 steps) = trunk + acrylamide; "
                "sulfur STEPS deliberately absent"
            ),
        },
        "lanes_compose": False,
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
