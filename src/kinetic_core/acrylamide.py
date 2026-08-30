"""
src/kinetic_core/acrylamide.py

THE ACRYLAMIDE / SAFETY NETWORK (Build Wave B3), Module 3 of the kinetic-core
rebuild.
=============================================================================

SIXTEEN steps bolted onto B1's fifteen-step trunk, over the extended state
vector of ``species_acrylamide.py``. Same contract as B1 and B2: explicit
reactant->product stoichiometry, carbon AND nitrogen AND sulfur balance
enforced at import, a formation term and a consumption term for every
intermediate, no DFT, no invented pH or a_w dependence.

THE ONE DEFECT THIS MODULE EXISTS TO REMOVE
-------------------------------------------
The shipped FAST lane's acrylamide had NO ELIMINATION TERM. An accumulating
pool with a formation rate and no sink is linear in time for as long as its
precursors last, so it responds to a longer bake by growing without limit --
which is why the audit measured it ~40x under-responsive to time relative to
data that has a real maximum. Acrylamide's actual kinetics are two consecutive
reactions with a MAXIMUM: Claeys 2005 measures the peak at 40-50 min at 160 C
and the elimination constant is 246x the formation constant. Here acrylamide
has three sinks and the peak is a consequence of the network rather than a
correction applied to it.

WHY THE SULFUR REACTIONS ARE *NOT* COMPOSED IN
----------------------------------------------
The state vector is B2's, so cysteine here IS B2's cysteine. The REACTIONS are
not: this network is B1's trunk plus B3's block, and B2's twenty-odd sulfur
steps are absent. That is deliberate and it is a chemistry decision, not a
scoping convenience.

De Vleeschouwer 2009 Part II measures a cysteine sink IN THE SYSTEM THIS MODULE
FITS (k_Y = 0.35 min^-1 at 160 C, Ea 110.5) and this module runs it. B2's four
thiol-consumption channels are measured in coffee brew at 30 C, in storage at
50 C and in thiamine/xylose pots at 115-120 C, and B2's formation steps are
fitted at pH 5 / 145 C. Running both blocks at once would spend the same
cysteine twice, through two channel sets calibrated on different systems, and
the acrylamide scavenging rate -- which is FIRST ORDER IN CYSTEINE -- would
inherit the error directly. The composed network is a later wave's problem and
it needs a decision about which cysteine sink is which; ``OUT_OF_SCOPE`` names
what that strands.

WHAT THIS MODULE CANNOT DO, STATED UP FRONT
-------------------------------------------
  * No sucrose, so the declared sucrose HOLD-OUT is UNSCOREABLE rather than
    wrong.
  * No glutamine PROMOTION mechanism. See DELIBERATE_UNDERFITS in
    parameters_acrylamide.py: the shape belongs to a hold-out.
  * No water-activity dependence, on a fit panel that spans a_w 0.35-1.0.
    Every run reports the gap between its own a_w and each parameter's.
  * No pH dependence, and specifically not the shipped ~2000x pH factor.

UNITS: mmol/L, minutes, Kelvin, kJ/mol.
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
from .parameters import SCHIFF_AMADORI_SPLIT
from .parameters_acrylamide import (
    ASN_SCHIFF_AMADORI_SPLIT,
    AcrylamideParameter,
    MEASURED_ACRYLAMIDE,
    T_REF_A_K,
    acrylamide_registry_metadata,
)
from .species_acrylamide import (
    ACRYLAMIDE_BY_KEY,
    ACRYLAMIDE_INDEX,
    ACRYLAMIDE_STATE_KEYS,
    ACRYLAMIDE_TERMINAL_POOLS,
    N_ACRYLAMIDE_STATE,
    acrylamide_ppb,
    initial_acrylamide_state,
    total_element_acrylamide,
)
from .ph_state import AMINE_FATE_BASIS as _AMINE_FATE_BASIS
from .species_sulfur import TERMINAL_POOLS as SULFUR_TERMINAL_POOLS


# ---------------------------------------------------------------------------
# What this wave leaves on the table
# ---------------------------------------------------------------------------
OUT_OF_SCOPE: Tuple[Dict[str, str], ...] = (
    {
        "lane": "sucrose (and any non-reducing sugar)",
        "what_is_stranded": (
            "De Vleeschouwer 2009 I's sucrose system, whose hydrolysis and "
            "isomerisation constants are measured. It is a declared HOLD-OUT, "
            "so the module may not carry it -- but that means the hold-out is "
            "reported UNSCOREABLE rather than predicted."
        ),
        "why": (
            "Adding a sucrose species would require its hydrolysis rate, and "
            "the only measurement of that rate sits in the hold-out column. "
            "Inventing one to be able to score the hold-out would defeat the "
            "hold-out."
        ),
    },
    {
        "lane": "the composed sulfur + acrylamide network",
        "what_is_stranded": (
            "the interaction between B2's four thiol sinks and B3's cysteine "
            "channels; any prediction of MFT/FFT and acrylamide in the same "
            "pot; the cysteine-rich extrusion systems where both matter."
        ),
        "why": (
            "Two independently calibrated cysteine sink sets would consume the "
            "same pool. Resolving which sink applies at which temperature and "
            "a_w is a decision about the corpus, not a coding step."
        ),
    },
    {
        "lane": "glutamine promotion of acrylamide",
        "what_is_stranded": (
            "Claeys' glutamine competitor row, which the model under-predicts "
            "by roughly its measured promotion factor."
        ),
        "why": (
            "The promotion's temperature SHAPE is a declared hold-out "
            "(inventory sec. B5.5). A term fitted to the FIT half would be a "
            "term built toward the hold-out."
        ),
    },
    {
        "lane": "asparagine's OTHER Maillard products",
        "what_is_stranded": (
            "everything between Int1 and the melanoidin pool: Strecker "
            "aldehydes, pyrazines, the 3-aminopropionamide route."
        ),
        "why": (
            "The Int1 partition is measured only as an aggregate (and its "
            "literature constant is one the authors mark 'no physical "
            "meaning'), so resolving the non-acrylamide branch would be fitting "
            "structure to one number."
        ),
    },
)


# ---------------------------------------------------------------------------
# The network
# ---------------------------------------------------------------------------
# Carbon accounting, once: several steps in this block report ONE product and
# not the rest -- De Vleeschouwer measures the acrylamide leaving Int1 but not
# the C7/N1 residue that leaves with it, and none of the three labs identifies a
# single acrylamide degradation product. Following B1, the unreported atoms are
# routed to the FRAG_C / FRAG_N / FRAG_S pools rather than deleted, so the
# balances close as equalities and the size of the unassigned pool is visible.

ACRYLAMIDE_REACTIONS: Tuple[Reaction, ...] = (
    # ---- the asparagine trunk branch --------------------------------------
    Reaction(
        "a_asn_glc_sb", {"Asn": 1, "Glc": 1}, {"SBA": 1}, "k_asn_glc",
        "THE bimolecular initiation. De Vleeschouwer 2009 I k_INTg. It is what "
        "makes the acrylamide yield respond to precursor CONCENTRATION instead "
        "of being a fixed fraction.",
    ),
    Reaction(
        "a_sb_int1", {"SBA": 1}, {"INT1": 1}, None,
        "Rate DERIVED from a_asn_glc_sb by the pinned ratio of "
        "ASN_SCHIFF_AMADORI_SPLIT, not independently parameterised. "
        "Irreversible and sink-free, so every condensed molecule reaches Int1 "
        "and the pair reduces exactly to the measured one-step composite.",
    ),
    Reaction(
        "a_int1_acr", {"INT1": 1}, {"ACR": 1, "FRAG_C": 7, "FRAG_N": 1},
        "k_int1_acr",
        "ACRYLAMIDE FORMATION. The C7 and the second nitrogen leave in "
        "co-products the source does not measure and are routed to the "
        "fragment pools.",
    ),
    Reaction(
        "a_int1_mel", {"INT1": 1}, {"MEL_C": 10, "MEL_N": 2}, "k_int1_mel",
        "FITTED. The non-acrylamide fate of the intermediate, and the partition "
        "that sets the yield. Its literature value exists and is REFUSED (the "
        "authors mark k_M 'no physical meaning'). Without this step every "
        "condensed molecule becomes acrylamide.",
    ),
    Reaction(
        "a_asn_asp", {"Asn": 1}, {"Asp": 1, "FRAG_N": 1}, "k_asn_asp",
        "Deamidation, the competing fate of asparagine. MEASURED.",
    ),
    Reaction(
        "a_asp_sink", {"Asp": 1}, {"FRAG_C": 4, "FRAG_N": 1}, "k_asp_sink",
        "Aspartic acid consumption. MEASURED (the one usable k_X in the "
        "corpus); products unidentified, so they go to the fragment pools.",
    ),
    # ---- ACRYLAMIDE ELIMINATION -- the term the old lane did not have ------
    Reaction(
        "a_acr_dp", {"ACR": 1}, {"FRAG_C": 3, "FRAG_N": 1}, "k_acr_dp",
        "First-order degradation. FITTED to three independent measurements of "
        "the same constant (Claeys k_E, De Vleeschouwer k_E, Knol k6 at five "
        "temperatures). No lab identifies the products.",
    ),
    Reaction(
        "a_acr_cys", {"ACR": 1, "Cys": 1}, {"ACRCYS": 1}, "k_acr_cys",
        "MEASURED bimolecular scavenging -- De Vleeschouwer II k_E2, the "
        "tightest parameter in the corpus. This is the step that couples the "
        "acrylamide lane to B2's cysteine pool.",
    ),
    Reaction(
        "a_cys_sink", {"Cys": 1},
        {"CBX": 1, "FRAG_C": 2, "FRAG_N": 1, "FRAG_S": 1},
        "k_cys_sink",
        "MEASURED cysteine sink of the SAME system (De Vleeschouwer II k_Y). "
        "The product is not identified by the source, and the note says so. "
        "B2.3 CHARGE FIX, identical in form and in reasoning to the sulfur "
        "lane's r_cys_thermal: an unidentified sink is not a licence to delete "
        "cysteine's carboxylate and ammonium, so both centres are CARRIED "
        "(CBX) while its amine centre is DECLARED destroyed in the ledger "
        "rather than dropped into an untitratable pool in silence. The "
        "acrylamide lane has no pH state today, so this changes no prediction "
        "here -- it exists so that the day it gets one, the defect B2.2 found "
        "in the sulfur lane is not waiting in this one.",
    ),
    # ---- competition: consumption of the SHARED sugar pool ------------------
    Reaction(
        "a_cys_glc", {"Cys": 1, "Glc": 1}, {"MEL_C": 9, "MEL_N": 1, "FRAG_S": 1},
        "k_cys_glc",
        "MEASURED competition (De Vleeschouwer II k_INT2). The template for "
        "the three fitted competitor channels below: a competing amino acid "
        "acts by removing glucose from the asparagine lane.",
    ),
    Reaction(
        "a_gln_glc", {"Gln": 1, "Glc": 1}, {"MEL_C": 11, "MEL_N": 2}, "k_gln_glc",
        "FITTED. Bounds allow ~zero.",
    ),
    Reaction(
        "a_lys_glc", {"Lys": 1, "Glc": 1}, {"MEL_C": 12, "MEL_N": 2}, "k_lys_glc",
        "FITTED. Bounds allow ~zero.",
    ),
    Reaction(
        "a_ala_glc", {"Ala": 1, "Glc": 1}, {"MEL_C": 9, "MEL_N": 1}, "k_ala_glc",
        "FITTED. Claeys' alanine row is indistinguishable from the control, so "
        "this one is expected at the bottom of its range.",
    ),
    # ---- competition: Michael-acceptor scavenging of acrylamide -------------
    Reaction(
        "a_acr_gln", {"ACR": 1, "Gln": 1}, {"FRAG_C": 8, "FRAG_N": 3}, "k_acr_gln",
        "FITTED pre-exponential; BARRIER NOT FITTED -- held at the measured "
        "Ea_E2 = 51.3 kJ/mol under the inventory's family licence.",
    ),
    Reaction(
        "a_acr_lys", {"ACR": 1, "Lys": 1}, {"FRAG_C": 9, "FRAG_N": 3}, "k_acr_lys",
        "FITTED pre-exponential, measured Ea_E2.",
    ),
    Reaction(
        "a_acr_ala", {"ACR": 1, "Ala": 1}, {"FRAG_C": 6, "FRAG_N": 2}, "k_acr_ala",
        "FITTED pre-exponential, measured Ea_E2.",
    ),
)

#: The full network: B1's trunk first, then the acrylamide block. B2's sulfur
#: STEPS are absent by decision (see the module docstring); B2's sulfur SPECIES
#: are present, because the state vector is shared.
FULL_ACRYLAMIDE_REACTIONS: Tuple[Reaction, ...] = (
    TRUNK_REACTIONS + ACRYLAMIDE_REACTIONS
)
FULL_ACRYLAMIDE_REACTION_KEYS: Tuple[str, ...] = tuple(
    r.key for r in FULL_ACRYLAMIDE_REACTIONS
)

#: Every step that CONSUMES acrylamide. The fit and hold-out reports read this
#: rather than a hand-maintained list, so a channel added later cannot be
#: silently left out of the apparent elimination constant.
ACRYLAMIDE_SINK_REACTIONS: Tuple[str, ...] = tuple(
    r.key for r in ACRYLAMIDE_REACTIONS if "ACR" in r.reactants
)

#: Every step that FORMS acrylamide.
ACRYLAMIDE_SOURCE_REACTIONS: Tuple[str, ...] = tuple(
    r.key for r in ACRYLAMIDE_REACTIONS if "ACR" in r.products
)


# ---------------------------------------------------------------------------
# Construction-time validation, over THREE elements
# ---------------------------------------------------------------------------

_TERMINAL = tuple(
    sorted(set(("FRAG_C", "MEL_C", "MEL_N")
               + SULFUR_TERMINAL_POOLS + ACRYLAMIDE_TERMINAL_POOLS))
)


def validate_acrylamide_balance(
    reactions: Sequence[Reaction] = FULL_ACRYLAMIDE_REACTIONS,
) -> None:
    """Raise unless every step conserves carbon, nitrogen AND sulfur."""
    for reaction in reactions:
        for element in BALANCED_ELEMENTS:
            left, right = reaction.atom_balance(element, ACRYLAMIDE_BY_KEY)
            if left != right:
                raise ValueError(
                    f"{reaction.key}: {element} does not balance "
                    f"({left} -> {right}). Every step must conserve atoms; the "
                    f"unmeasured residue belongs in FRAG_C / FRAG_N / FRAG_S, "
                    f"not in nowhere."
                )
        for key in list(reaction.reactants) + list(reaction.products):
            if key not in ACRYLAMIDE_INDEX:
                raise ValueError(f"{reaction.key}: unknown species {key!r}")
        for pool in _TERMINAL:
            if pool in reaction.reactants:
                raise ValueError(
                    f"{reaction.key}: {pool} is a terminal pool and must never "
                    f"be a reactant."
                )


validate_acrylamide_balance()


# ---------------------------------------------------------------------------
# B2.3: the acrylamide lane's half of the CENTRE LEDGER
# ---------------------------------------------------------------------------
# The audit is run here for the same reason it is run in `sulfur`: Amendment 9
# says ALL sites, and this lane consumes B2's cysteine through three steps. Two
# of them needed no change (a_acr_cys hands both centres to a real amino acid,
# S-(2-carbamoylethyl)cysteine, which is now in the centre table; a_cys_sink is
# fixed above), and the third is a DECLARED destruction whose basis is the same
# one the trunk's melanoidin sink already relies on.
ACRYLAMIDE_CENTRE_LEDGER: Mapping[str, Mapping[str, object]] = {
    **TRUNK_CENTRE_LEDGER,
    "a_cys_sink": {"amine": -1, "basis": (
        "De Vleeschouwer II's k_Y cysteine sink, whose product the source does "
        "not identify. Its CARBOXYL is carried (CBX) -- an unidentified "
        "product is not a licence to delete a carboxylate -- and its amine is "
        "declared destroyed on the shared basis. " + _AMINE_FATE_BASIS)},
    "a_cys_glc": {"carboxyl": -1, "amine": -1, "basis": (
        "MELANOIDIN INCORPORATION. Cysteine's carboxylate and ammonium enter "
        "the terminal browning polymer, which really does retain them -- "
        "melanoidins are polyanionic -- but the pool is carried ELEMENTALLY "
        "(mmol of atom, not of repeat unit), so no per-unit centre count is "
        "definable and the charge balance cannot see them. This is the SAME "
        "declared gap as the trunk's r_tdg_mel + glycine, recorded in "
        "ph_state.UNTRACKED_TITRATABLE['MEL_N']. It is a gap, not a licence: "
        "the acrylamide lane carries no pH state, so it costs no prediction "
        "today and becomes a defect the day it gets one.")},
}


def _validate_acrylamide_charge_closure() -> None:
    from .ph_state import validate_charge_closure

    validate_charge_closure(FULL_ACRYLAMIDE_REACTIONS, ACRYLAMIDE_CENTRE_LEDGER)


_validate_acrylamide_charge_closure()


def acrylamide_stoichiometric_matrix(
    reactions: Sequence[Reaction] = FULL_ACRYLAMIDE_REACTIONS,
) -> np.ndarray:
    matrix = np.zeros((N_ACRYLAMIDE_STATE, len(reactions)), dtype=float)
    for j, reaction in enumerate(reactions):
        for key, coefficient in reaction.reactants.items():
            matrix[ACRYLAMIDE_INDEX[key], j] -= float(coefficient)
        for key, coefficient in reaction.products.items():
            matrix[ACRYLAMIDE_INDEX[key], j] += float(coefficient)
    return matrix


ACRYLAMIDE_STOICHIOMETRY: np.ndarray = acrylamide_stoichiometric_matrix()


# ---------------------------------------------------------------------------
# Rate evaluation
# ---------------------------------------------------------------------------


def acrylamide_rate_constants_at(
    parameters: Mapping[str, Any], temperature_k: float
) -> Dict[str, float]:
    """
    Every reaction's rate constant at ``temperature_k``.

    There is NO pH argument, deliberately (see the parameter registry's policy
    3). Two constants are DERIVED rather than looked up: B1's Amadori
    rearrangement and this module's asparagine Schiff-base rearrangement, both
    pinned to their own condensation step.
    """
    out: Dict[str, float] = {}
    for reaction in FULL_ACRYLAMIDE_REACTIONS:
        key = reaction.parameter_key
        if key is None:
            continue
        parameter = parameters.get(key)
        if parameter is None:
            raise KeyError(f"{reaction.key}: no parameter {key!r} supplied")
        if parameter.k_ref is None or parameter.ea_kj_mol is None:
            raise ValueError(
                f"{reaction.key}: parameter {key!r} is unpopulated "
                f"(evidence_class={parameter.evidence_class}). The fitted steps "
                f"must be given values before the network can be integrated; "
                f"there is no silent default."
            )
        out[reaction.key] = parameter.k_at(float(temperature_k))

    # DERIVED 1 -- B1's Amadori rearrangement, pinned to its condensation.
    out["r_amadori"] = (
        float(SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
        * out["r_schiff"]
        * float(SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"])
    )
    # DERIVED 2 -- the asparagine Schiff-base rearrangement, same device.
    out["a_sb_int1"] = (
        float(ASN_SCHIFF_AMADORI_SPLIT["ratio_amadori_over_schiff_pseudo_first_order"])
        * out["a_asn_glc_sb"]
        * float(ASN_SCHIFF_AMADORI_SPLIT["amine_loading_mmol_L_for_the_ratio"])
    )
    return out


_A_REACTANT_LAYOUT: Tuple[Tuple[Tuple[int, int], ...], ...] = tuple(
    tuple((ACRYLAMIDE_INDEX[key], int(c)) for key, c in r.reactants.items())
    for r in FULL_ACRYLAMIDE_REACTIONS
)
_A_RATE_KEY_ORDER: Tuple[str, ...] = tuple(
    r.key for r in FULL_ACRYLAMIDE_REACTIONS
)


def acrylamide_rate_vector(k: Mapping[str, float]) -> np.ndarray:
    """The rate constants as a positional array in network order."""
    return np.array([float(k[key]) for key in _A_RATE_KEY_ORDER], dtype=float)


def acrylamide_reaction_rates(state: np.ndarray, k) -> np.ndarray:
    """Mass-action rate of every reaction, mmol/(L*min)."""
    y = np.asarray(state, dtype=float)
    if y.min() < 0.0:
        y = np.clip(y, 0.0, None)
    kv = k if isinstance(k, np.ndarray) else acrylamide_rate_vector(k)
    rates = np.empty(len(_A_REACTANT_LAYOUT), dtype=float)
    for j, layout in enumerate(_A_REACTANT_LAYOUT):
        value = kv[j]
        for index, coefficient in layout:
            value = value * (y[index] if coefficient == 1 else y[index] ** coefficient)
        rates[j] = value
    return rates


def acrylamide_derivatives(state: np.ndarray, k) -> np.ndarray:
    return ACRYLAMIDE_STOICHIOMETRY @ acrylamide_reaction_rates(state, k)


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AcrylamideRun:
    """One integrated acrylamide experiment, with its conservation audit."""

    times_min: np.ndarray
    concentrations: np.ndarray
    temperature_k: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def series(self, species_key: str) -> np.ndarray:
        return self.concentrations[:, ACRYLAMIDE_INDEX[species_key]]

    def final(self, species_key: str) -> float:
        return float(self.concentrations[-1, ACRYLAMIDE_INDEX[species_key]])

    def acrylamide_ppb(self) -> np.ndarray:
        """The acrylamide trajectory in ppb -- the unit every source prints."""
        return np.array([acrylamide_ppb(v) for v in self.series("ACR")])

    def peak_acrylamide_ppb(self) -> float:
        return float(np.max(self.acrylamide_ppb()))

    def time_of_peak_min(self) -> float:
        return float(self.times_min[int(np.argmax(self.series("ACR")))])

    def molar_yield_on_asparagine(self) -> float:
        """Peak acrylamide as a fraction of the asparagine initially charged."""
        asn0 = float(self.metadata.get("initial_asparagine_mmol_L", 0.0))
        if asn0 <= 0.0:
            return float("nan")
        return float(np.max(self.series("ACR"))) / asn0

    def element_closure(self, element: str) -> Dict[str, float]:
        totals = np.array(
            [total_element_acrylamide(row, element) for row in self.concentrations]
        )
        t0 = float(totals[0])
        drift = float(np.max(np.abs(totals - t0)))
        return {
            f"initial_{element}_mmol_L": t0,
            "max_abs_drift_mmol_L": drift,
            "max_relative_drift": drift / t0 if t0 > 0 else float("nan"),
        }


def _acrylamide_warnings(
    parameters: Mapping[str, Any], temperature_c: float, water_activity: Optional[float]
) -> list:
    """
    Temperature AND water-activity extrapolation warnings.

    The a_w half is what makes this module's warnings different from B1's. Its
    fit panel spans a_w 0.35 (extrusion) to 1.0 (Claeys), the De Vleeschouwer
    constants that carry most of the chemistry were measured at 0.92, and there
    is no a_w term. Every prediction therefore reports how far it is from the
    water activity its constants were measured at, on every run.
    """
    out = []
    for key, parameter in sorted(parameters.items()):
        low, high = parameter.temperature_range_c
        if temperature_c < low - 1e-9 or temperature_c > high + 1e-9:
            out.append(
                f"{key}: evaluated at {temperature_c:.1f} C, measured over "
                f"{low:.0f}-{high:.0f} C ({parameter.source_anchor[:60]}...)"
            )
        measured_aw = getattr(parameter, "aw_of_measurement", None)
        if (
            water_activity is not None
            and measured_aw is not None
            and abs(float(measured_aw) - float(water_activity)) > 0.1
        ):
            out.append(
                f"{key}: evaluated at a_w {float(water_activity):.2f}, measured "
                f"at a_w {float(measured_aw):.2f} -- and this module has NO a_w "
                f"term, because nothing in its fit corpus measures one step at "
                f"two water activities."
            )
    return out


def integrate_acrylamide(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    times_min: Sequence[float],
    *,
    water_activity: Optional[float] = None,
    method: str = "LSODA",
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> AcrylamideRun:
    """
    Integrate the acrylamide network at one temperature.

    ``water_activity`` is METADATA ONLY: it changes no rate. It is accepted so
    that the a_w gap between a run and its parameters is reported rather than
    silently ignored, which is the honest encoding of a corpus that spans
    a_w 0.35-1.0 without measuring the axis.
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
    y0 = initial_acrylamide_state(initial)
    k_fixed = acrylamide_rate_vector(
        acrylamide_rate_constants_at(parameters, float(temperature_k))
    )

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        return acrylamide_derivatives(y, k_fixed)

    solution = solve_ivp(
        rhs, (0.0, duration), y0, t_eval=grid, method=method, rtol=rtol, atol=atol
    )
    if not solution.success:
        raise RuntimeError(f"acrylamide-network integration failed: {solution.message}")

    raw = solution.y.T
    worst_negative = float(np.min(raw)) if raw.size else 0.0
    if worst_negative < -max(1e3 * atol, 1e-8):
        raise RuntimeError(
            f"acrylamide-network integration produced a state of "
            f"{worst_negative:.3e}, far below the absolute tolerance "
            f"{atol:.1e}: a genuine non-negativity failure, not solver noise."
        )
    concentrations = np.clip(raw, 0.0, None)

    temperature_c = float(temperature_k) - 273.15
    metadata: Dict[str, Any] = {
        "method": method,
        "temperature_C": temperature_c,
        "temperature_K": float(temperature_k),
        "water_activity_of_the_run": water_activity,
        "n_species": N_ACRYLAMIDE_STATE,
        "n_reactions": len(FULL_ACRYLAMIDE_REACTIONS),
        "species_order": list(ACRYLAMIDE_STATE_KEYS),
        "initial_asparagine_mmol_L": float(initial.get("Asn", 0.0)),
        "worst_raw_negative_before_clip": worst_negative,
        "extrapolation_warnings": _acrylamide_warnings(
            parameters, temperature_c, water_activity
        ),
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
    }
    metadata.update(acrylamide_registry_metadata(
        {k: v for k, v in parameters.items() if isinstance(v, AcrylamideParameter)}
    ))

    run = AcrylamideRun(
        times_min=grid,
        concentrations=concentrations,
        temperature_k=float(temperature_k),
        metadata=metadata,
    )
    for element in BALANCED_ELEMENTS:
        metadata[f"{element}_closure"] = run.element_closure(element)
    return run


def acrylamide_flux_budget(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    duration_min: float,
    *,
    n_points: int = 201,
    rtol: float = 1e-8,
    atol: float = 1e-14,
) -> Dict[str, float]:
    """Time-integrated flux through every reaction, mmol/L."""
    grid = np.linspace(0.0, float(duration_min), int(n_points))
    run = integrate_acrylamide(
        parameters, temperature_k, initial, grid, rtol=rtol, atol=atol
    )
    k = acrylamide_rate_vector(
        acrylamide_rate_constants_at(parameters, float(temperature_k))
    )
    rates = np.array([acrylamide_reaction_rates(row, k) for row in run.concentrations])
    integral = np.trapezoid(rates, grid, axis=0)
    return {
        reaction.key: float(integral[j])
        for j, reaction in enumerate(FULL_ACRYLAMIDE_REACTIONS)
    }


# ---------------------------------------------------------------------------
# THE APPARENT LUMPED CONSTANTS -- how a mass-action model is compared to a
# two-consecutive-first-order fit
# ---------------------------------------------------------------------------
# Claeys 2005 reports ONE pair of first-order constants per system, from the
# model  dC_R/dt = -k_F*C_R ;  dC_AA/dt = k_F*C_R - k_E*C_AA.
# This network does not have those constants; it has elementary steps. The two
# are compared by asking what k_F and k_E a first-order description of THIS
# model's own trajectory would need, which is a flux-weighted average:
#
#     k_F_app = (total acrylamide FORMED over the window) / (integral of [Asn])
#     k_E_app = (total acrylamide REMOVED over the window) / (integral of [ACR])
#
# For a system that really is first order these are EXACT, at any window --
# which is the property the unit test pins. For this network they are not
# exact, and the two places they differ are stated rather than hidden:
#
#   * the formation constant is biased LOW at low temperature, because the
#     intermediate needs time to build up and a first-order law has no
#     induction period. The bias is largest at 140 C and smallest at 200 C, so
#     it inflates the apparent formation ACTIVATION ENERGY;
#   * the elimination constant in a competitor system is a genuine average over
#     a falling competitor concentration, which is exactly what Claeys' own
#     single fitted k_E is.
#
# `window_min` is fixed at the same value for every temperature, so no
# literature value enters the choice of window.


def apparent_lumped_constants(
    parameters: Mapping[str, Any],
    temperature_k: float,
    initial: Mapping[str, float],
    window_min: float,
    *,
    n_points: int = 61,
    water_activity: Optional[float] = None,
    rtol: float = 1e-8,
    atol: float = 1e-14,
) -> Dict[str, float]:
    """
    The apparent first-order formation and elimination constants of a run.

    Returns k_F_app and k_E_app in 1/min, plus the peak acrylamide, the time of
    the peak and the molar yield, so that one integration serves every
    observable the fit and hold-out reports need.
    """
    grid = np.linspace(0.0, float(window_min), int(n_points))
    run = integrate_acrylamide(
        parameters, temperature_k, initial, grid,
        water_activity=water_activity, rtol=rtol, atol=atol,
    )
    k = acrylamide_rate_vector(
        acrylamide_rate_constants_at(parameters, float(temperature_k))
    )
    rates = np.array([acrylamide_reaction_rates(row, k) for row in run.concentrations])
    index = {key: j for j, key in enumerate(_A_RATE_KEY_ORDER)}

    formation = np.zeros(len(grid))
    for key in ACRYLAMIDE_SOURCE_REACTIONS:
        formation = formation + rates[:, index[key]]
    elimination = np.zeros(len(grid))
    per_channel: Dict[str, float] = {}
    for key in ACRYLAMIDE_SINK_REACTIONS:
        channel = rates[:, index[key]]
        elimination = elimination + channel
        per_channel[key] = float(np.trapezoid(channel, grid))

    asn = run.series("Asn")
    acr = run.series("ACR")
    int_formation = float(np.trapezoid(formation, grid))
    int_elimination = float(np.trapezoid(elimination, grid))
    int_asn = float(np.trapezoid(asn, grid))
    int_acr = float(np.trapezoid(acr, grid))

    return {
        "k_F_app_per_min": int_formation / int_asn if int_asn > 0 else float("nan"),
        "k_E_app_per_min": int_elimination / int_acr if int_acr > 0 else float("nan"),
        "peak_acrylamide_ppb": run.peak_acrylamide_ppb(),
        "time_of_peak_min": run.time_of_peak_min(),
        "final_acrylamide_ppb": acrylamide_ppb(float(acr[-1])),
        "molar_yield_on_asparagine": run.molar_yield_on_asparagine(),
        "elimination_flux_by_channel_mmol_L": per_channel,
        "window_min": float(window_min),
        "n_extrapolation_warnings": len(run.metadata["extrapolation_warnings"]),
    }


def apparent_activation_energy(
    k_low: float, temperature_low_k: float, k_high: float, temperature_high_k: float
) -> float:
    """
    The two-point Arrhenius barrier implied by an apparent constant, kJ/mol.

    TWO POINTS IS A DELIBERATE CHOICE HERE AND IT IS NOT THE PROHIBITED
    DERIVATION OF B2's sec. C.1. There, two rates belonging to two DIFFERENT
    reactions were being paired. Here the two rates are the SAME model
    observable at two temperatures, computed from the same parameter set, and
    the quantity being formed is a PREDICTION to be compared with a published
    Ea -- not a parameter that enters the model.
    """
    from math import log

    from .parameters import R_KJ

    if k_low <= 0.0 or k_high <= 0.0:
        return float("nan")
    return float(
        R_KJ
        * log(float(k_high) / float(k_low))
        / (1.0 / float(temperature_low_k) - 1.0 / float(temperature_high_k))
    )


def describe_acrylamide() -> Dict[str, object]:
    """A machine-readable description of the network, for the reports."""
    return {
        "n_species": N_ACRYLAMIDE_STATE,
        "n_reactions": len(FULL_ACRYLAMIDE_REACTIONS),
        "trunk_reactions": len(TRUNK_REACTIONS),
        "acrylamide_reactions": len(ACRYLAMIDE_REACTIONS),
        "sulfur_reactions_composed_in": 0,
        "acrylamide_source_reactions": list(ACRYLAMIDE_SOURCE_REACTIONS),
        "acrylamide_sink_reactions": list(ACRYLAMIDE_SINK_REACTIONS),
        "reference_temperature_K": T_REF_A_K,
        "measured_parameter_keys": sorted(MEASURED_ACRYLAMIDE),
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
                "carbon_balance": list(r.atom_balance("carbon", ACRYLAMIDE_BY_KEY)),
                "nitrogen_balance": list(r.atom_balance("nitrogen", ACRYLAMIDE_BY_KEY)),
                "sulfur_balance": list(r.atom_balance("sulfur", ACRYLAMIDE_BY_KEY)),
                "note": r.note,
            }
            for r in ACRYLAMIDE_REACTIONS
        ],
        "out_of_scope": [dict(row) for row in OUT_OF_SCOPE],
    }
