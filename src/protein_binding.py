"""Measured protein-flavour binding physics for the matrix observability lane.

2026-08-27 (Wave S4).  See tasks/audit_remediation.md, "Wave S4".

WHAT THIS MODULE IS FOR
-----------------------
The matrix-only lane predicts a total volatile concentration in a protein slurry and then
multiplies it by an *observability factor* to get the number a paper reports.  Until this
wave every one of those factors was FITTED: back-solved from the benchmark it is then
scored against, and in two cases (the 1-hexanol pair `0.143 / 0.063`) back-solved from
values that appear in no publication at all (Wave T3).

This module is the alternative, and it contains NO FITTED PARAMETER.  It computes the
observability factor from measured protein-binding data:

    f_free = 1 / (1 + K_eff * c_p)

a single-site Langmuir isotherm in the dilute-ligand limit, where `c_p` is the protein
concentration in the SAME mass units `K_eff` was measured in.  Every constant comes from
``data/lit/binding_constants.yml``, where each carries a verbatim quote from a retrieved
full text.

THE TWO HALVES OF THE PHYSICS, AND WHY THE SECOND ONE MATTERS AS MUCH AS THE FIRST
---------------------------------------------------------------------------------
`f_free` alone does not tell you the observability factor.  You also have to know WHICH
concentration the source you are being scored against actually reported:

  * A **matrix-matched** quantification (standards spiked into the slurry, or isotope
    dilution) measures the TOTAL analyte.  Binding is already inside the response factor,
    so the correct observability factor is **1.0** and `f_free` must NOT be applied.
    Pratap-Singh 2021 is this case ("The pea protein sample was spiked with 1, 2, 3, 4,
    5 uL of stock standard and 5 uL of internal standard hexanal-d12 was used to generate
    the standard curve").
  * A **water-calibrated** quantification (external curve built in water, applied to a
    protein slurry) reads back the FREE concentration.  The correct factor is `f_free`.
    Liu 2021 is this case ("pea protein solutions were replaced with 5 mL DI water and
    were spiked with known compounds at five different levels").
  * **unknown** -> no correction is applied, and the row says so.  Guessing here would be
    the same class of error this campaign exists to remove.

WHAT IS NOT MODELLED, AND WHY
-----------------------------
1. **No denaturation / heat term.**  The mission allowed one only if a source quantifies
   it.  Two sources address it and they DISAGREE IN SIGN:
     - Barallat-Perez 2023 (commercial isolates, hexanal): "increasing the temperature
       also led to a slight increase in the level of hexanal binding" -- qualitative, no
       number.
     - Heng 2005 (purified pea vicilin, octanal): 96% bound non-heated vs 32% bound after
       90 C / 30 min -- a number, but for a purified fraction, and pointing the other way.
   A term whose sign the literature disputes is not a term.  Recorded in the YAML under
   `denaturation_effect_evidence`; NOT applied.
2. **No low-moisture / extrudate lane.**  This is an aqueous-phase partition model.  At
   a_w 0.35 there is no bulk aqueous phase for it to describe, and volatile retention in a
   dry extrudate is sorption, not aqueous binding.  Those lanes get
   `in_domain = False` and factor 1.0.  The out-of-domain value is still computed and
   reported by the comparison generator so the size of what is being declined is visible.
3. **No covalent term.**  Aldehydes form Schiff bases with lysine; that channel is real
   (it is why heptanal/octanal in Barallat-Perez 2023 span 13.7-52.8% over one CH2, which
   no Pow-only model reproduces) and nothing retrieved quantifies it per gram of pea or
   soy protein.  The a_p route was fitted on KETONES, so for an aldehyde it is a LOWER
   BOUND on binding and therefore an UPPER BOUND on `f_free`.

BOUNDARY WITH THE REST OF THE PIPELINE (stated because the mission required it)
------------------------------------------------------------------------------
This module touches the OBSERVABILITY layer only.  It does not read, write or influence
the volatile-budget / allocation machinery (`src/projection.py`,
`MATRIX_BENCHMARK_BASE_MARKER_YIELDS`, `src/lipid_oxidation.py`), which is a separate
[P] workstream.  A consequence worth stating plainly: because the base marker yields were
themselves built from the (partly fabricated) 260 / 380 / 80 / 120 ppb values, the
absolute scale of the matrix lane still lives in those yields, and no binding model can
repair it from the observability side.
"""

from __future__ import annotations

import math
import os
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Dict, Iterator, Mapping, Optional

import yaml

from src import data_paths

BINDING_CONSTANTS_PATH = data_paths.BINDING_CONSTANTS

# --- Observability modes -------------------------------------------------------------
#
# `fitted_factors` is the INCUMBENT and remains the default.  Nothing about the shipped
# prediction path changes unless a caller explicitly selects another mode.
MODE_FITTED = "fitted_factors"
MODE_BINDING = "binding_physics"
MODE_BINDING_ALL = "binding_physics_out_of_domain"
# `unit_observability` is the NULL MODEL: observability factor 1.0 everywhere, i.e. "the
# model predicts the total and the paper reports the total". It exists because without it
# a win for `binding_physics` over `fitted_factors` is unattributable: on most lanes the
# binding mode has no binding data to apply and reduces to exactly this, so any
# improvement there measures the fitted factors being worse than nothing, not the binding
# physics being right. Scoring all three separates the two claims.
MODE_UNIT = "unit_observability"
OBSERVABILITY_MODES = (MODE_FITTED, MODE_BINDING, MODE_BINDING_ALL, MODE_UNIT)

_MODE_ENV = "MAILLARD_MATRIX_OBSERVABILITY_MODE"
_MODE_OVERRIDE: Optional[str] = None


def observability_mode() -> str:
    """The active matrix observability mode.

    Resolution order: explicit process-level override (``use_observability_mode``), then
    the ``MAILLARD_MATRIX_OBSERVABILITY_MODE`` environment variable, then the shipped
    default ``fitted_factors``.
    """
    if _MODE_OVERRIDE is not None:
        return _MODE_OVERRIDE
    raw = str(os.environ.get(_MODE_ENV, "") or "").strip()
    if raw in OBSERVABILITY_MODES:
        return raw
    return MODE_FITTED


@contextmanager
def use_observability_mode(mode: str) -> Iterator[str]:
    """Temporarily select an observability mode (generators and tests only)."""
    global _MODE_OVERRIDE
    if mode not in OBSERVABILITY_MODES:
        raise ValueError(f"unknown observability mode {mode!r}; expected one of {OBSERVABILITY_MODES}")
    previous = _MODE_OVERRIDE
    _MODE_OVERRIDE = mode
    try:
        yield mode
    finally:
        _MODE_OVERRIDE = previous


def binding_mode_active() -> bool:
    """True when the observability factor should come from this module, not the registry.

    Includes the null model: `unit_observability` also bypasses the fitted registry (it
    returns 1.0), so the same no-double-count assertion covers it.
    """
    return observability_mode() in {MODE_BINDING, MODE_BINDING_ALL, MODE_UNIT}


# --- Data loading --------------------------------------------------------------------

_CACHE: Dict[str, Any] = {}


def load_binding_constants() -> Dict[str, Any]:
    """Parsed ``data/lit/binding_constants.yml`` (cached)."""
    if "payload" not in _CACHE:
        if not BINDING_CONSTANTS_PATH.exists():
            _CACHE["payload"] = {}
        else:
            with open(BINDING_CONSTANTS_PATH, "r", encoding="utf-8") as handle:
                _CACHE["payload"] = yaml.safe_load(handle) or {}
    return _CACHE["payload"]


def _records() -> list[Mapping[str, Any]]:
    return [r for r in load_binding_constants().get("records", []) if isinstance(r, Mapping)]


def octanol_water_partition(compound: str) -> Optional[Mapping[str, Any]]:
    """The Pow entry for a compound, or None when the file has no verified value."""
    normalized = str(compound).strip().lower()
    for entry in load_binding_constants().get("octanol_water_partition", {}).get("entries", []):
        if str(entry.get("compound", "")).strip().lower() == normalized:
            return entry
    return None


def matrix_protein_loading(lane_id: Optional[str]) -> Optional[Mapping[str, Any]]:
    """The declared protein loading / quantification basis for a benchmark lane."""
    if not lane_id:
        return None
    for entry in load_binding_constants().get("matrix_protein_loadings", []):
        if str(entry.get("lane_id", "")) == str(lane_id):
            return entry
    return None


def hydrophobic_interaction_parameter(protein_source: str, compound_class: str = "ketone") -> Optional[Mapping[str, Any]]:
    """The Harrison-Hills a_p record (L per g protein) for a protein source."""
    for record in _records():
        if record.get("quantity") != "hydrophobic_interaction_parameter_ap":
            continue
        if str(record.get("protein_source")) != str(protein_source):
            continue
        if str(record.get("compound_class")) != str(compound_class):
            continue
        return record
    return None


# --- The model -----------------------------------------------------------------------

# Isolate protein purity is needed because `a_p` is per gram of PROTEIN while several
# sources state their loading per gram of ISOLATE POWDER.  Rather than invent a purity,
# the model works on whatever basis the loading record declares and reports the basis it
# used; where the basis is powder the resulting c_p is an UPPER bound on the protein
# concentration, so f_free is a LOWER bound.  The size of that bound is <=1.3x (commercial
# plant isolates run 77-93 wt% protein in the sources this file cites).
POWDER_TO_PROTEIN_UPPER_BOUND = 1.0


@dataclass(frozen=True)
class FreeFraction:
    """Result of a binding calculation, with everything needed to audit it."""

    f_free: float
    in_domain: bool
    mechanism: str
    reasons: tuple[str, ...] = ()
    k_eff_l_per_g: Optional[float] = None
    protein_concentration_g_per_l: Optional[float] = None
    protein_basis: Optional[str] = None
    pow_value: Optional[float] = None
    pow_basis: Optional[str] = None
    record_ids: tuple[str, ...] = ()
    quantification_basis: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "binding_f_free": float(self.f_free),
            "binding_in_domain": bool(self.in_domain),
            "binding_mechanism": self.mechanism,
            "binding_reasons": list(self.reasons),
            "binding_k_eff_l_per_g": self.k_eff_l_per_g,
            "binding_protein_concentration_g_per_l": self.protein_concentration_g_per_l,
            "binding_protein_basis": self.protein_basis,
            "binding_pow": self.pow_value,
            "binding_pow_basis": self.pow_basis,
            "binding_record_ids": list(self.record_ids),
            "binding_quantification_basis": self.quantification_basis,
        }


def free_fraction_from_ap(
    compound: str,
    *,
    protein_source: str,
    protein_concentration_g_per_l: float,
    protein_basis: str = "g_protein",
    compound_class: str = "ketone",
) -> FreeFraction:
    """f_free from the Harrison-Hills hydrophobic partition parameter.

        K_eff = a_p * Pow          [L per g protein]
        f_free = 1 / (1 + K_eff * c_p)

    `compound_class` selects WHICH a_p record is used and is the class the parameter was
    MEASURED on -- not the class of `compound`.  Passing "ketone" for an aldehyde is a
    declared transfer and makes the answer an upper bound on f_free (aldehydes bind at
    least as strongly as the ketone of the same Pow, never less).
    """
    reasons: list[str] = []
    ap_record = hydrophobic_interaction_parameter(protein_source, compound_class)
    if ap_record is None:
        return FreeFraction(
            f_free=1.0,
            in_domain=False,
            mechanism="none",
            reasons=(f"no a_p record for protein_source={protein_source!r} class={compound_class!r}",),
        )
    pow_entry = octanol_water_partition(compound)
    if pow_entry is not None and pow_entry.get("pow") is None:
        pow_entry = None
    if pow_entry is None:
        return FreeFraction(
            f_free=1.0,
            in_domain=False,
            mechanism="none",
            reasons=(f"no content-verified octanol-water partition coefficient for {compound!r}",),
            record_ids=(str(ap_record.get("record_id")),),
        )

    a_p = float(ap_record["value"])
    pow_value = float(pow_entry["pow"])
    k_eff = a_p * pow_value
    c_p = float(protein_concentration_g_per_l)
    f_free = 1.0 / (1.0 + k_eff * c_p)

    if str(pow_entry.get("basis")) != "experimental":
        reasons.append(f"Pow for {compound} is {pow_entry.get('basis')}, not an experimental measurement")
    if protein_basis == "g_isolate_powder":
        reasons.append(
            "c_p is on an ISOLATE-POWDER basis while a_p is per gram of PROTEIN; c_p is "
            "therefore an upper bound (<=1.3x) and f_free a lower bound"
        )
    ap_range = (ap_record.get("conditions") or {}).get("protein_range_g_isolate_per_kg")
    if ap_range and c_p > float(ap_range[1]):
        reasons.append(
            f"c_p {c_p:.1f} g/L extrapolates {c_p / float(ap_range[1]):.1f}x above the "
            f"highest protein level a_p was fitted at ({ap_range[1]} g/kg); the model is "
            "linear in c_p so this is an extrapolation of the fit range, not of the form"
        )

    return FreeFraction(
        f_free=f_free,
        in_domain=True,
        mechanism="harrison_hills_ap_pow",
        reasons=tuple(reasons),
        k_eff_l_per_g=k_eff,
        protein_concentration_g_per_l=c_p,
        protein_basis=protein_basis,
        pow_value=pow_value,
        pow_basis=str(pow_entry.get("basis")),
        record_ids=(str(ap_record.get("record_id")), f"pow:{compound}"),
    )


def k_eff_from_percent_bound(percent_bound: float, protein_concentration_g_per_l: float) -> float:
    """Invert a percent-bound-at-conditions record into K_eff (L/g).

    Exact arithmetic under the same single-site form:  p/(1-p) = K_eff * c_p.
    This is NOT an extrapolation into an affinity constant -- it is the same equation
    read backwards at the conditions the percentage was measured at, and it is only
    valid to re-apply at a different c_p if the single-site form holds, which is the
    model's stated assumption.
    """
    p = float(percent_bound) / 100.0
    if not (0.0 < p < 1.0):
        raise ValueError(f"percent_bound must be strictly between 0 and 100, got {percent_bound}")
    return (p / (1.0 - p)) / float(protein_concentration_g_per_l)


def resolve_binding_context(bench: Mapping[str, Any]) -> Dict[str, Any]:
    """Build the per-lane binding context from a benchmark bundle + the YAML.

    Returns a dict with `lane_id`, `protein_source`, `protein_concentration_g_per_l`,
    `protein_basis`, `aqueous`, `quantification_basis` and `reasons`.  Everything in it
    comes from the YAML's `matrix_protein_loadings`, i.e. from each lane's OWN source --
    never from this repository's guess.
    """
    lane_id = str(bench.get("benchmark_id", "") or "")
    entry = matrix_protein_loading(lane_id)
    if entry is None:
        return {
            "lane_id": lane_id,
            "protein_source": str(bench.get("protein_type", "") or ""),
            "protein_concentration_g_per_l": None,
            "protein_basis": None,
            "aqueous": None,
            "quantification_basis": "unknown",
            "reasons": [f"no declared protein loading for lane {lane_id!r} in data/lit/binding_constants.yml"],
        }
    return {
        "lane_id": lane_id,
        "protein_source": str(entry.get("protein_source") or bench.get("protein_type") or ""),
        "protein_concentration_g_per_l": entry.get("protein_concentration_g_per_L"),
        "protein_basis": entry.get("protein_basis"),
        "aqueous": entry.get("aqueous"),
        "quantification_basis": str(entry.get("quantification_basis", "unknown")),
        "reasons": [],
    }


def observability_factor(
    compound: str,
    *,
    context: Mapping[str, Any],
    mode: Optional[str] = None,
) -> FreeFraction:
    """The binding-physics observability factor for one compound on one lane.

    Returns 1.0 (with `in_domain=False` and a stated reason) whenever the physics does not
    apply or the data to apply it does not exist.  It NEVER falls back to a fitted factor:
    a hybrid of the two would be unattributable, and the no-double-count assertion in
    `src/benchmark_validation.py` forbids it.
    """
    active = mode or observability_mode()
    reasons = list(context.get("reasons") or [])

    if active == MODE_UNIT:
        return FreeFraction(
            f_free=1.0,
            in_domain=True,
            mechanism="unit_null_model",
            reasons=("null model: no observability correction of any kind",),
            quantification_basis=str(context.get("quantification_basis", "unknown")),
        )

    basis = str(context.get("quantification_basis", "unknown"))
    if basis == "matrix_matched":
        return FreeFraction(
            f_free=1.0,
            in_domain=True,
            mechanism="matrix_matched_reports_total",
            reasons=tuple(reasons + [
                "the reference quantification is matrix-matched, so it reports the TOTAL "
                "analyte and no binding correction applies"
            ]),
            quantification_basis=basis,
        )
    if basis == "unknown":
        return FreeFraction(
            f_free=1.0,
            in_domain=False,
            mechanism="none",
            reasons=tuple(reasons + [
                "the reference's quantification basis is not established, so it is not "
                "known whether the reported number is a total or a free concentration"
            ]),
            quantification_basis=basis,
        )

    aqueous = context.get("aqueous")
    if aqueous is False and active != MODE_BINDING_ALL:
        return FreeFraction(
            f_free=1.0,
            in_domain=False,
            mechanism="none",
            reasons=tuple(reasons + [
                "lane is not an aqueous dispersion (low water activity / structured "
                "solid); an aqueous-partition Langmuir does not describe it"
            ]),
            quantification_basis=basis,
        )

    c_p = context.get("protein_concentration_g_per_l")
    if c_p is None:
        return FreeFraction(
            f_free=1.0,
            in_domain=False,
            mechanism="none",
            reasons=tuple(reasons + ["no declared protein loading for this lane"]),
            quantification_basis=basis,
        )

    result = free_fraction_from_ap(
        compound,
        protein_source=str(context.get("protein_source") or ""),
        protein_concentration_g_per_l=float(c_p),
        protein_basis=str(context.get("protein_basis") or "g_protein"),
        compound_class="ketone",
    )
    return FreeFraction(
        f_free=result.f_free,
        in_domain=result.in_domain,
        mechanism=result.mechanism,
        reasons=tuple(reasons + list(result.reasons) + [
            "a_p was measured on 2-alkanones; applying it to a non-ketone is a declared "
            "class transfer and makes f_free an UPPER bound (aldehydes also bind covalently)"
        ]),
        k_eff_l_per_g=result.k_eff_l_per_g,
        protein_concentration_g_per_l=result.protein_concentration_g_per_l,
        protein_basis=result.protein_basis,
        pow_value=result.pow_value,
        pow_basis=result.pow_basis,
        record_ids=result.record_ids,
        quantification_basis=basis,
    )


# --- Independent cross-checks (no fitting, reported not tuned) ------------------------

def cross_check_percent_bound() -> list[Dict[str, Any]]:
    """Predict every `percent_bound_at_conditions` record from the a_p route and compare.

    This is the wave's zero-parameter validation of the binding model: the a_p values were
    fitted by Snel et al. on 2-alkanones in demineralised water with APCI-MS; the
    percent-bound records come from other laboratories, other methods (static and dynamic
    headspace depletion) and in three cases other CHEMICAL CLASSES.  Nothing here is
    tuned; the residuals are whatever they are.
    """
    out: list[Dict[str, Any]] = []
    for record in _records():
        if record.get("quantity") != "percent_bound_at_conditions":
            continue
        conditions = record.get("conditions") or {}
        c_p_powder = conditions.get("protein_concentration_g_per_L")
        if c_p_powder is None:
            continue
        compound = str(record.get("compound"))
        measured = float(record["value"])
        purity = conditions.get("protein_purity_fraction")
        c_p = float(c_p_powder) * float(purity) if purity else float(c_p_powder)
        row: Dict[str, Any] = {
            "record_id": record.get("record_id"),
            "compound": compound,
            "compound_class": record.get("compound_class"),
            "protein_source_declared": record.get("protein_source"),
            "measured_percent_bound": measured,
            "protein_concentration_g_per_l_used": c_p,
            "protein_basis": conditions.get("protein_basis"),
        }
        for source in ("pea_iso", "soy_iso"):
            pred = free_fraction_from_ap(
                compound,
                protein_source=source,
                protein_concentration_g_per_l=c_p,
                protein_basis=str(conditions.get("protein_basis") or "g_protein"),
            )
            if pred.in_domain:
                row[f"predicted_percent_bound_{source}"] = 100.0 * (1.0 - pred.f_free)
                row[f"residual_points_{source}"] = 100.0 * (1.0 - pred.f_free) - measured
            else:
                row[f"predicted_percent_bound_{source}"] = None
                row[f"unavailable_because_{source}"] = "; ".join(pred.reasons)
        out.append(row)
    return out


def describe_model() -> Dict[str, Any]:
    """Machine-readable description of the model, for artifacts and gates."""
    payload = load_binding_constants()
    return {
        "module": "src/protein_binding.py",
        "form": "f_free = 1 / (1 + K_eff * c_p)  [single-site Langmuir, dilute-ligand limit]",
        "k_eff_route": "K_eff = a_p * Pow  (Harrison & Hills partition model, a_p in L per g protein)",
        "fitted_parameters": 0,
        "fitted_parameters_note": (
            "Every constant is a literature transcription with a verbatim quote in "
            "data/lit/binding_constants.yml. Nothing in this module was chosen to make any "
            "prediction match any benchmark."
        ),
        "denaturation_modelled": False,
        "denaturation_not_modelled_because": (
            "the two retrieved sources disagree in sign (Barallat-Perez 2023: heat raises "
            "hexanal binding, qualitative; Heng 2005: heat lowers octanal binding to "
            "purified vicilin, 96% -> 32%)"
        ),
        "covalent_binding_modelled": False,
        "record_count": len(_records()),
        "source_count": len(payload.get("sources", []) or []),
        "modes": list(OBSERVABILITY_MODES),
        "default_mode": MODE_FITTED,
    }


__all__ = [
    "MODE_FITTED",
    "MODE_BINDING",
    "MODE_BINDING_ALL",
    "MODE_UNIT",
    "OBSERVABILITY_MODES",
    "FreeFraction",
    "binding_mode_active",
    "cross_check_percent_bound",
    "describe_model",
    "free_fraction_from_ap",
    "hydrophobic_interaction_parameter",
    "k_eff_from_percent_bound",
    "load_binding_constants",
    "matrix_protein_loading",
    "observability_factor",
    "observability_mode",
    "octanol_water_partition",
    "resolve_binding_context",
    "use_observability_mode",
]
