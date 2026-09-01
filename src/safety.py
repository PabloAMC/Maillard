"""
src/safety.py

Quantitative safety modeling for Maillard by-products.
Focus: Acrylamide formation from Asparagine + Reducing Sugars.

Reference:
- Knol et al. 2009 J. Ag. Food Chem. (Kinetic model for acrylamide)
- Stadler et al. 2004 (Asparagine involvement)

UNITS CONTRACT (declared once here, 2026-08-27; see `UNITS` below)
------------------------------------------------------------------
Every quantity crossing this module's boundary has exactly one unit, and it is
declared in the module-level `UNITS` dict. The short version:

* INPUTS  `*_mM`            -> millimolar (mM) free-precursor concentration.
* INPUTS  `temp_C`          -> degrees Celsius; `time_min` -> minutes.
* OUTPUTS `predict_acrylamide(...).acrylamide_ppb`
          `predict_cml(...)`, `predict_cel(...)`, `predict_furosine(...)`
          -> ppb of the FOOD MATRIX, i.e. ug/kg. ppb == ug/kg == 1e-3 mg/kg.
* OUTPUT  `evaluate_formulation_safety(...)[0]` -> dimensionless risk in [0, 1].

HONESTY NOTE (2026-08-27 unit-collision repair). Before this date the CML/CEL/
furosine predictors were consumed as ppb by the benchmark lane and as mg/kg (or
mg per 100 g protein) by src/literature_runtime.py. The two benchmark JSONs had
mg/kg source values written into `conc_ppb` fields, which manufactured ~1.0-1.2x
agreement out of a 1000x unit error. The unit is now declared as ppb on this
side, the benchmark files carry true ppb, and the runtime anchors are converted
to ppb before any ratio is taken. The consequence is that the CML/CEL/furosine
predictors are now visibly ~1e3x/2e2x BELOW the literature bands they claim to
match: their empirical scale constants were tuned against the mis-unit-ed files
and have never been calibrated against a real ug/kg measurement. That failure is
left visible in the benchmark panel rather than papered over; do not "fix" it by
rescaling the constants without a real calibration set.
"""

import json
import math
import re
import warnings
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple, Optional

from src.extrusion import normalize_moisture_regime
from src import data_paths
from src import data_access


SAFETY_REFERENCE_PAYLOAD_PATH = data_paths.SAFETY_REFERENCE_PAYLOADS


#: Single source of truth for the units of everything this module exchanges.
#: Consumers (benchmark lane, literature runtime, reporting) must read this
#: rather than inferring a unit from a field name.
UNITS: Dict[str, str] = {
    "precursor_concentration": "mM",
    "temperature": "degC",
    "time": "min",
    "predict_acrylamide.acrylamide_ppb": "ppb (ug/kg food)",
    "predict_acrylamide.uncertainty_ppb": "ppb (ug/kg food)",
    "predict_cml": "ppb (ug/kg food)",
    "predict_cel": "ppb (ug/kg food)",
    "predict_furosine": "ppb (ug/kg food)",
    "evaluate_formulation_safety.risk": "dimensionless [0, 1]",
    "extrusion total_damage_load.furosine_mg_per_kg": "mg/kg food",
    "extrusion total_damage_load.lal_mg_per_kg": "mg/kg food",
}

#: mg/kg -> ppb. Kept explicit so conversions are greppable.
MG_PER_KG_TO_PPB = 1000.0

#: Declared protein basis used to convert `mg per 100 g protein` anchors (the
#: Ramirez-Jimenez furosine crossover point) into ppb of food. 20 g protein per
#: 100 g product is the mid-point of the commercial PBMA burger/extrudate range
#: (typically 15-25 g/100 g). It is a DECLARED basis, not a measured one, and it
#: is the only place the conversion is defined.
FUROSINE_ANCHOR_PROTEIN_FRACTION = 0.20


def mg_per_100g_protein_to_ppb(value: float, protein_fraction: float = FUROSINE_ANCHOR_PROTEIN_FRACTION) -> float:
    """Convert `mg analyte / 100 g protein` to ppb (ug/kg) of food.

    mg/100 g protein * (protein_fraction g protein / g food) = mg/100 g food
    mg/100 g food * 10 = mg/kg food; mg/kg * 1000 = ppb.
    e.g. 8.7 mg/100 g protein at 20% protein -> 8.7 * 0.20 * 10 * 1000 = 17400 ppb.
    """
    return float(value) * float(protein_fraction) * 10.0 * MG_PER_KG_TO_PPB


def _load_safety_reference_payloads() -> dict:
    return data_access.load_json(SAFETY_REFERENCE_PAYLOAD_PATH)


SAFETY_REFERENCE_PAYLOADS = _load_safety_reference_payloads()


_DEFAULT_SAFETY_VISIBILITY_BY_KIND = {
    "industrial_endpoint_reference": "default",
    "finished_product_reference": "default",
    "precursor_correlation_reference": "default",
    "kinetic_model_reference": "extended",
    "contextual_process_modifier": "extended",
}


def _normalize_safety_visibility(entry: Dict[str, Any]) -> str:
    explicit = str(entry.get("report_visibility", "")).strip().lower()
    if explicit in {"default", "extended"}:
        return explicit
    kind = str(entry.get("kind", "")).strip()
    return _DEFAULT_SAFETY_VISIBILITY_BY_KIND.get(kind, "extended")


#: Analyte strings in the payload use four different separators for composite
#: analytes: "CML/CEL", "furosine_cml", "cml_cel_acrylamide", "CML_CEL".
#: Splitting on "/" and "+" alone (the pre-2026-08-27 behaviour) made every
#: underscore-joined composite invisible to the reference-context builder.
_ANALYTE_SEPARATORS = re.compile(r"[+/_,;\s]+")


def _entry_analyte_field(entry: Dict[str, Any]) -> str:
    """Resolve the analyte string of a payload entry.

    Seven entries in data/lit/safety_reference_payloads.json carry no `analyte`
    field at all (they use the older `category` schema). They were previously
    invisible to every analyte-filtered query. `category` is used as the
    fallback; entries with neither are reported as excluded rather than silently
    dropped (see `build_safety_reference_context`).
    """
    for key in ("analyte", "category"):
        value = str(entry.get(key, "") or "").strip()
        if value:
            return value
    return ""


def _analyte_tokens(analyte_field: str) -> set:
    normalized = str(analyte_field).strip().lower()
    if not normalized:
        return set()
    tokens = {token.strip() for token in _ANALYTE_SEPARATORS.split(normalized) if token.strip()}
    tokens.add(normalized)
    return tokens


def _analyte_matches(entry_analyte: str, requested_analyte: Optional[str]) -> bool:
    if not requested_analyte:
        return True
    normalized_requested = str(requested_analyte).strip().lower()
    return normalized_requested in _analyte_tokens(entry_analyte)


def _reference_band_ppb(*, analyte: str) -> Dict[str, float]:
    """Collect the ug/kg (== ppb) reference band for an analyte from the payload.

    Returns {"attention_ppb", "action_ppb", "source_ids"}: the lowest lower bound
    and the highest upper bound over every `matrix_reference_ranges` row whose
    units are ug_per_kg. These are the bands that replaced the dead 100-ppb
    literal that used to gate `SafetyResult.flagged`.
    """
    lows: List[float] = []
    highs: List[float] = []
    source_ids: List[str] = []
    for entry in SAFETY_REFERENCE_PAYLOADS.get("entries", []):
        if not _analyte_matches(_entry_analyte_field(entry), analyte):
            continue
        for row in entry.get("matrix_reference_ranges", []) or []:
            if str(row.get("units", "")).strip().lower() != "ug_per_kg":
                continue
            values = [
                float(row[key])
                for key in ("min", "max", "mean", "point")
                if row.get(key) is not None
            ]
            if not values:
                continue
            lows.append(min(values))
            highs.append(max(values))
            entry_id = str(entry.get("id", ""))
            if entry_id and entry_id not in source_ids:
                source_ids.append(entry_id)
    if not lows or not highs:
        return {"attention_ppb": 0.0, "action_ppb": 0.0, "source_ids": []}
    return {
        "attention_ppb": float(min(lows)),
        "action_ppb": float(max(highs)),
        "source_ids": source_ids,
    }


_ACRYLAMIDE_BAND = _reference_band_ppb(analyte="acrylamide")

#: Lowest acrylamide level actually reported in a plant-protein product by the
#: payload band (Fu et al. 2023 mixed PBMA lower bound, 31.81 ug/kg). A
#: prediction whose upper uncertainty edge reaches this level is flagged: it is
#: at the level of real finished products.
#: NOTE ON REGULATION: Commission Regulation (EU) 2017/2158 Annex IV has NO meat
#: analogue / plant-protein category. The closest REAL Annex IV categories for an
#: extruded cereal-legume matrix are "breakfast cereals (maize/oat/spelt/barley/
#: rice based)" at 150 ug/kg and "crackers other than potato based" at 400 ug/kg.
#: The band used here sits below both; any Annex IV comparison is ANALOGICAL, not
#: regulatory, and no Annex IV number is used as a threshold in this module.
ACRYLAMIDE_ATTENTION_PPB = float(_ACRYLAMIDE_BAND["attention_ppb"] or 31.81)

#: Top of the observed industrial band (Squeo et al. 2023 soy wet-extraction
#: upper bound, 748 ug/kg). The risk sub-score saturates here.
ACRYLAMIDE_ACTION_PPB = float(_ACRYLAMIDE_BAND["action_ppb"] or 748.0)

#: Reporting floor for "the analyte was assessed and is present". Declared, not
#: payload-derived: it is a routine LC-MS/MS acrylamide reporting floor and sits
#: three decades below every band in the payload.
ACRYLAMIDE_REPORTING_FLOOR_PPB = 1.0

#: Furosine attention level, converted from the only furosine anchor in the
#: payload (Ramirez-Jimenez 2000 crossover peak, 8.7 mg/100 g protein) on the
#: declared 20% protein basis -> 17400 ppb == 17.4 mg/kg.
FUROSINE_ATTENTION_MG_PER_KG = mg_per_100g_protein_to_ppb(8.7) / MG_PER_KG_TO_PPB
#: Saturation at 5x the mild-extrusion peak. DECLARED severity ceiling, not
#: source-derived: the payload has no severe-extrusion furosine band.
FUROSINE_ACTION_MG_PER_KG = 5.0 * FUROSINE_ATTENTION_MG_PER_KG

#: PROVISIONAL: data/lit/safety_reference_payloads.json contains NO lysinoalanine
#: reference band, so these two levels are declared placeholders chosen around the
#: ingredient-inherited LAL baselines carried by src/extrusion.py (26-88 mg/kg).
#: They are NOT band-derived and must be replaced once a LAL payload lands.
LAL_PROVISIONAL_ATTENTION_MG_PER_KG = 50.0
LAL_PROVISIONAL_ACTION_MG_PER_KG = 250.0

#: Damage markers (furosine, LAL) are only flagged once the formulation has seen
#: real thermal exposure. DECLARED onset, not source-derived: below ~100 C an
#: aqueous/high-moisture matrix does not accumulate Amadori-derived furosine on
#: process timescales. Without this gate a 25 C / 0 min formulation was flagged
#: for Furosine and LAL purely from the ingredient-inherited baseline that
#: src/extrusion.py attaches to every profile. This is deliberately stricter than
#: src.extrusion.AMBIENT_THERMAL_THRESHOLD_C (50 C, "is this a process at all"):
#: this gate asks the narrower question "is this a process that forms Amadori-
#: derived damage markers". Where an extrusion profile supplies
#: `damage_load_attribution`, the ingredient-inherited share of the total load is
#: available there and is NOT process-attributable.
DAMAGE_MARKER_THERMAL_ONSET_C = 100.0


def get_safety_reference_payload(reference_id: str = "squeo_2023_pbpi_acrylamide") -> Optional[dict]:
    for entry in SAFETY_REFERENCE_PAYLOADS.get("entries", []):
        if str(entry.get("id", "")) == reference_id:
            return entry
    return None


def get_safety_reference_entries(
    *,
    analyte: Optional[str] = None,
    visibility: str = "default",
) -> List[dict]:
    requested_visibility = str(visibility).strip().lower()
    entries: List[dict] = []
    for entry in SAFETY_REFERENCE_PAYLOADS.get("entries", []):
        entry_analyte = _entry_analyte_field(entry).strip().lower()
        if not _analyte_matches(entry_analyte, analyte):
            continue
        entry_visibility = _normalize_safety_visibility(entry)
        if requested_visibility == "default" and entry_visibility != "default":
            continue
        entries.append(dict(entry, report_visibility=entry_visibility))
    return entries


def build_safety_reference_context(*, analyte: str = "acrylamide") -> Dict[str, Any]:
    """Reference context for one analyte, plus an explicit exclusion disclosure.

    2026-08-27: composite analyte strings are now split on `_` as well as `+`
    and `/`, and entries that carry no `analyte` field fall back to `category`.
    Entries that resolve to no analyte at all cannot be filtered and are listed
    under `excluded_entries` so a report can state that evidence was excluded
    instead of silently dropping it.
    """
    def _summarize(entry: Dict[str, Any]) -> Dict[str, Any]:
        supports = entry.get("what_it_supports", []) or []
        summary = str(supports[0]) if supports else str(entry.get("comment", ""))
        return {
            "id": str(entry.get("id", "unknown")),
            "kind": str(entry.get("kind", "unknown")),
            "report_visibility": str(entry.get("report_visibility", _normalize_safety_visibility(entry))),
            "source_citation": str(entry.get("source_citation", "unknown")),
            "summary": summary,
        }

    default_entries = [_summarize(entry) for entry in get_safety_reference_entries(analyte=analyte, visibility="default")]
    extended_entries = [_summarize(entry) for entry in get_safety_reference_entries(analyte=analyte, visibility="all") if _normalize_safety_visibility(entry) == "extended"]
    # Exclusion disclosure. Two classes, both previously invisible:
    #  (a) entries with neither an `analyte` nor a `category` field - unmatchable;
    #  (b) legacy `category`-schema entries whose coarse category ("ages",
    #      "protein_damage_markers") does not resolve to the requested analyte.
    included_ids = {row["id"] for row in default_entries} | {row["id"] for row in extended_entries}
    excluded_entries: List[Dict[str, Any]] = []
    for entry in SAFETY_REFERENCE_PAYLOADS.get("entries", []):
        entry_id = str(entry.get("id", "unknown"))
        if entry_id in included_ids:
            continue
        if entry.get("analyte"):
            continue  # a real analyte that simply is not the requested one
        analyte_field = _entry_analyte_field(entry)
        excluded_entries.append(
            {
                "id": entry_id,
                "reason": "no_analyte_or_category_field" if not analyte_field else "category_only_entry_not_matched",
                "category": analyte_field,
                "source_citation": str(entry.get("source_citation", "unknown")),
            }
        )
    exclusion_note = ""
    if excluded_entries:
        exclusion_note = (
            f"{len(excluded_entries)} safety reference entries were EXCLUDED from this "
            f"{analyte} context: they carry no `analyte` field, and their legacy `category` "
            "value does not resolve to the requested analyte. Their evidence is not "
            "represented above; see excluded_entries."
        )
    out = {
        "default_entries": default_entries,
        "extended_entries": extended_entries,
        "excluded_entries": excluded_entries,
        "exclusion_note": exclusion_note,
    }
    # Wave B8 (Amendment 16 clause 2): the unsourced-Arrhenius flag TRAVELS on
    # the payload, not only to stderr -- the `literature_runtime` shape. It is
    # attached for the analytes those four pairs actually reach.
    if str(analyte).strip().lower() in {"cml", "cel", "furosine", "ages",
                                        "protein_damage_markers"}:
        out["unsourced_arrhenius_pairs"] = dict(SAFETY_ARRHENIUS_PROVENANCE)
    return out


def get_safety_reference_range(matrix_family: str, reference_id: str = "squeo_2023_pbpi_acrylamide") -> Optional[dict]:
    payload = get_safety_reference_payload(reference_id)
    if payload is None:
        return None
    for item in payload.get("matrix_reference_ranges", []):
        if str(item.get("matrix_family", "")) == matrix_family:
            return item
    return None


_NAME_TOKEN_SPLIT = re.compile(r"[^a-z0-9]+")


def _name_tokens(name: Any) -> List[str]:
    return [token for token in _NAME_TOKEN_SPLIT.split(str(name).strip().lower()) if token]


def _has_word(name: Any, word: str) -> bool:
    """Word-boundary ingredient matching (plural-tolerant).

    Substring matching used to parse "asnase" (asparaginase, an acrylamide
    MITIGATION enzyme) as `asn` == asparagine, i.e. as an acrylamide PRECURSOR
    at the enzyme's own concentration. Whole-token matching removes that class
    of error while still matching "L-Asparagine", "D-Glucose" and "reducing
    sugars".
    """
    target = str(word).strip().lower()
    if not target:
        return False
    tokens = _name_tokens(name)
    if " " in target:
        return target in " ".join(tokens)
    return any(token == target or token == target + "s" for token in tokens)


def _has_any_word(name: Any, words: List[str]) -> bool:
    return any(_has_word(name, word) for word in words)


def _normalize_runtime_tokens(values: Any) -> List[str]:
    if not isinstance(values, list):
        return []
    return [str(value).strip().lower() for value in values if str(value).strip()]


#: Half-saturation constant for amino-acid capping mitigation (mM). DECLARED,
#: not source-derived: PMC12648097 reports a suppression endpoint but not a dose
#: response, so a Michaelis-shaped saturating curve is assumed with a 10 mM half
#: point, normalised so a 50 mM capping dose reaches ~100% of the ceiling.
MITIGATION_HALF_SATURATION_MM = 10.0
MITIGATION_REFERENCE_DOSE_MM = 50.0
#: Ceilings retained from the presence-only implementation (0.65 for cysteine +
#: glycine together, 0.20 for one alone). They are deliberately well below the
#: 0.978 suppression fraction stored in the payload, which was measured in a
#: capped model system at an unreported (high) dose.
MITIGATION_CEILING_BOTH = 0.65
MITIGATION_CEILING_SINGLE = 0.20


def _mitigation_dose_response(total_mM: float) -> float:
    """Saturating dose factor in [0, 1) for amino-acid capping mitigation.

    VALID RANGE NOTE: the underlying source (PMC12648097 / 10.1016/j.fochx.2025.103252)
    characterises capping at a single, high, unreported dose in a model system. Only
    the SHAPE here is defensible - monotone, saturating, and near-zero at trace
    additions. Absolute values below ~1 mM and above ~50 mM are extrapolation.
    """
    dose = max(0.0, float(total_mM))
    if dose <= 0.0:
        return 0.0
    saturating = dose / (MITIGATION_HALF_SATURATION_MM + dose)
    reference = MITIGATION_REFERENCE_DOSE_MM / (MITIGATION_HALF_SATURATION_MM + MITIGATION_REFERENCE_DOSE_MM)
    return float(min(1.0, saturating / reference))


def _resolve_acrylamide_runtime_adjustment(
    precursors: Dict[str, float],
    modifiers: Dict[str, Any],
) -> Dict[str, Any]:
    runtime_context = modifiers.get("__runtime_context__", {}) if isinstance(modifiers.get("__runtime_context__", {}), dict) else {}
    additives = _normalize_runtime_tokens(runtime_context.get("additives", []))
    interventions = _normalize_runtime_tokens(runtime_context.get("interventions", []))

    # Dose-dependent since 2026-08-27: a trace of cysteine no longer earns the
    # same 65% discount as a 50 mM capping dose.
    cys_mM = 0.0
    gly_mM = 0.0
    for name, concentration in precursors.items():
        try:
            value = max(0.0, float(concentration))
        except (TypeError, ValueError):
            continue
        if _has_word(name, "cysteine") or _has_word(name, "cys"):
            cys_mM += value
        if _has_word(name, "glycine") or _has_word(name, "gly"):
            gly_mM += value
    cys_present = cys_mM > 0.0
    gly_present = gly_mM > 0.0
    saponin_active = any("saponin" in token or "quillaja" in token for token in additives + interventions)

    mitigation_fraction = 0.0
    mitigation_ceiling = 0.0
    reference_ids: List[str] = []
    if cys_present and gly_present:
        mitigation_ceiling = MITIGATION_CEILING_BOTH
        reference_ids.append("pmc_12648097_acrylamide_mitigation")
    elif cys_present or gly_present:
        mitigation_ceiling = MITIGATION_CEILING_SINGLE
        reference_ids.append("pmc_12648097_acrylamide_mitigation")
    dose_factor = _mitigation_dose_response(cys_mM + gly_mM)
    mitigation_fraction = mitigation_ceiling * dose_factor

    process_modifier = 1.0
    if saponin_active:
        process_modifier = 0.88
        reference_ids.append("kocadagli_2016_saponin_acrylamide_modifier")

    return {
        "mitigation_fraction": float(mitigation_fraction),
        "mitigation_ceiling": float(mitigation_ceiling),
        "mitigation_dose_factor": float(dose_factor),
        "mitigator_dose_mM": float(cys_mM + gly_mM),
        "process_modifier": float(process_modifier),
        "reference_ids": reference_ids,
    }

@dataclass
class SafetyResult:
    """Acrylamide prediction. Units: see `UNITS` (ppb == ug/kg of food).

    `flagged` is True when the upper uncertainty edge of the prediction reaches
    the finished-product reference band (`ACRYLAMIDE_ATTENTION_PPB`).
    `assessed` is False when the analyte could not be evaluated at all (no
    asparagine and/or no reducing sugar); a 0.0 in that case means "not
    assessed", NOT "measured clean".
    """

    acrylamide_ppb: float
    uncertainty_ppb: float
    flagged: bool
    description: str
    assessed: bool = True


def _acrylamide_moisture_factor(
    water_activity: Optional[float],
    moisture_regime: Optional[str],
) -> float:
    if water_activity is None:
        return 1.0

    aw = max(0.05, min(0.98, float(water_activity)))
    regime = normalize_moisture_regime(moisture_regime, aw)
    if regime == "lme":
        progress = max(0.0, min(1.0, (aw - 0.20) / 0.20))
        return 0.70 + 0.70 * progress

    progress = max(0.0, min(1.0, (aw - 0.50) / 0.45))
    return 1.20 - 0.50 * progress


def _append_unique(flagged: List[str], marker: str) -> None:
    if marker not in flagged:
        flagged.append(marker)


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


def _sigmoid(value: float, center: float, width: float) -> float:
    if width <= 0.0:
        return 1.0 if value >= center else 0.0
    exponent = _clamp((float(value) - center) / width, -60.0, 60.0)
    return 1.0 / (1.0 + math.exp(-exponent))


def _resolve_effective_temp_c(temp_C: float, effective_temp_c: Optional[float] = None) -> float:
    return float(temp_C if effective_temp_c is None else effective_temp_c)


def _arrhenius_rate(*, temp_c: float, pre_exponential: float, activation_energy_kj_mol: float) -> float:
    temperature_kelvin = float(temp_c) + 273.15
    return float(pre_exponential) * math.exp(-(float(activation_energy_kj_mol) * 1000.0) / (8.314 * temperature_kelvin))


# ===========================================================================
# THE FOUR UNCITED (Ea, A) PAIRS -- LABELLED, NOT MOVED (Wave B8)
# ===========================================================================
# FIT_HOLDOUT_DECLARATION.md Amendment 16 clause 2, ratified 2026-08-30:
# "safety.py's 4 uncited (Ea, A) pairs: to be LABELLED unsourced per the Wave T3
# convention in the next wave; conduct disclosed here."
#
# `results/validation/prefactor_audit.md` enumerates 53 shipped (Ea, A) pairs
# and classes exactly four of them `hardcoded_in_code_uncited` -- the strictest
# class in the audit, meaning there is not even a comment naming a paper:
#
#   safety_dicarbonyl_pool_formation      Ea 52,  A 4.0e6   (_estimate_dicarbonyl_pools)
#   safety_dicarbonyl_pool_elimination    Ea 74,  A 2.5e9   (_estimate_dicarbonyl_pools)
#   safety_furosine_formation             Ea 50,  A 8.0e5   (predict_furosine)
#   safety_furosine_elimination           Ea 73,  A 1.8e10  (predict_furosine)
#
# The Wave T3 convention (src/matrix_correction.py::_warn_if_registry_unsourced,
# src/literature_runtime.py::PROTEIN_SOURCE_PROVENANCE) is: state the defect,
# name what depends on it downstream, declare that NO VALUE IS SUBSTITUTED OR
# RESCALED, emit a RuntimeWarning, and carry a provenance record the flag can
# travel on. All five parts are honoured below. FOUR PAIRS ARE LABELLED AND
# ZERO NUMBERS MOVE.
#
# WHY NOTHING IS SUBSTITUTED: there is no source to substitute FROM. Inventing a
# citation-shaped value would be the defect again, one layer deeper. And these
# constants are not free-standing measurements -- they are the shape parameters
# of a formation/elimination signal whose SCALE is separately fitted (`scale=`),
# so "correcting" one member of a pair without the others would move predictions
# without improving provenance.
#
# ⚠ DO NOT REFLOW THE FOUR CALL SITES. `generate_prefactor_audit.py` re-parses
# those literals out of this file POSITIONALLY, with a regex that requires the
# four keywords on four consecutive lines in the printed order, and treats the
# two matches as (dicarbonyl, furosine) in file order. Hoisting them into named
# constants or reordering them silently drops rows from the audit.
SAFETY_UNCITED_ARRHENIUS_PAIRS: Tuple[str, ...] = (
    "safety_dicarbonyl_pool_formation",
    "safety_dicarbonyl_pool_elimination",
    "safety_furosine_formation",
    "safety_furosine_elimination",
)

#: Wave T3 flag value, verbatim, so the repo-wide census idiom matches.
SAFETY_ARRHENIUS_SOURCE_STATUS = "no_verifiable_source"

#: The provenance record. Emitted on the safety payload so the status travels
#: with the numbers it contaminates instead of living only in a warning nobody
#: sees (the `literature_runtime.PROTEIN_SOURCE_PROVENANCE` shape).
SAFETY_ARRHENIUS_PROVENANCE: Dict[str, object] = {
    "parameters": list(SAFETY_UNCITED_ARRHENIUS_PAIRS),
    "source_status": SAFETY_ARRHENIUS_SOURCE_STATUS,
    "where": "src/safety.py::_estimate_dicarbonyl_pools, ::predict_furosine",
    "value_basis": (
        "hardcoded_in_code_uncited -- no paper, no table and not even a comment "
        "names a source for any of the four (Ea, A) pairs"
    ),
    "declared_upstream": "results/validation/prefactor_audit.md (the enumeration)",
    "affects": [
        "predict_cml", "predict_cel", "predict_furosine",
        "the Amadori/dicarbonyl pool shape every one of them divides by",
    ],
    "warning": (
        "UNSOURCED SHAPE CONSTANTS. Differences these produce between two "
        "process conditions are not evidence; they reproduce a curve someone "
        "wrote down. Labelled by Wave B8 per FIT_HOLDOUT_DECLARATION.md "
        "Amendment 16 clause 2. VALUES ARE NOT SUBSTITUTED OR RESCALED."
    ),
    "known_miscalibration": (
        "predict_furosine sits ~2.0e2x below its own Ramirez-Jimenez anchor "
        "(see its docstring). That is REPORTED, not corrected, and it is an "
        "independent defect from the missing provenance."
    ),
}

#: Emitted at most once per process, and it names only the PARAMETER KEYS, not
#: their values -- the `literature_runtime._assess_concentration_unit` design
#: rule, and it matters here because both call sites run inside optimiser
#: sweeps, where a per-value message would emit thousands of times.
_SAFETY_ARRHENIUS_WARNED = False


def _warn_safety_arrhenius_unsourced(consumer: str) -> bool:
    """Surface the four uncited (Ea, A) pairs at first use. Wave T3 shape."""
    global _SAFETY_ARRHENIUS_WARNED
    if _SAFETY_ARRHENIUS_WARNED:
        return True
    _SAFETY_ARRHENIUS_WARNED = True
    warnings.warn(
        f"{consumer}: the Arrhenius pairs "
        f"{', '.join(SAFETY_UNCITED_ARRHENIUS_PAIRS)} are declared "
        f"source_status='{SAFETY_ARRHENIUS_SOURCE_STATUS}' -- no paper, no "
        "table and no comment names a source for any of them "
        "(results/validation/prefactor_audit.md). They are nonetheless LIVE: "
        "they set the Amadori/dicarbonyl pool shape that predict_cml, "
        "predict_cel and predict_furosine all divide by, so any CML, CEL or "
        "furosine difference between two process conditions is unanchored on "
        "that axis. Values are NOT substituted or rescaled.",
        RuntimeWarning,
        stacklevel=3,
    )
    return True


def _formation_elimination_signal(
    precursor_drive: float,
    *,
    temp_c: float,
    time_min: float,
    formation_pre_exponential: float,
    formation_ea_kj_mol: float,
    elimination_pre_exponential: float,
    elimination_ea_kj_mol: float,
    scale: float = 1.0,
    provenance_consumer: Optional[str] = None,
) -> float:
    if provenance_consumer is not None:
        _warn_safety_arrhenius_unsourced(provenance_consumer)
    if precursor_drive <= 0.0 or time_min <= 0.0:
        return 0.0

    time_seconds = float(time_min) * 60.0
    k_form = _arrhenius_rate(
        temp_c=temp_c,
        pre_exponential=formation_pre_exponential,
        activation_energy_kj_mol=formation_ea_kj_mol,
    )
    k_elim = _arrhenius_rate(
        temp_c=temp_c,
        pre_exponential=elimination_pre_exponential,
        activation_energy_kj_mol=elimination_ea_kj_mol,
    )
    if k_elim < 1e-12:
        signal = float(precursor_drive) * k_form * time_seconds
    else:
        signal = float(precursor_drive) * (k_form / k_elim) * (1.0 - math.exp(-k_elim * time_seconds))
    return max(0.0, float(scale) * signal)


def _estimate_dicarbonyl_pools(
    *,
    lysine_mM: float,
    reducing_sugar_mM: float,
    temp_C: float,
    time_min: float,
    water_activity: Optional[float] = None,
    effective_temp_c: Optional[float] = None,
    polyphenol_factor: float = 0.0,
) -> Tuple[float, float, float]:
    if lysine_mM <= 0.0 or reducing_sugar_mM <= 0.0 or time_min <= 0.0:
        return 0.0, 0.0, 0.0

    thermal_temp_c = _resolve_effective_temp_c(temp_C, effective_temp_c)
    aw = 0.55 if water_activity is None else _clamp(water_activity, 0.05, 0.98)
    aw_norm = _clamp((aw - 0.25) / 0.65, 0.0, 1.0)
    precursor_drive = math.sqrt(float(lysine_mM) * float(reducing_sugar_mM)) / 45.0
    amadori_pool = _formation_elimination_signal(
        precursor_drive,
        temp_c=thermal_temp_c,
        time_min=time_min,
        formation_pre_exponential=4.0e6,
        formation_ea_kj_mol=52.0,
        elimination_pre_exponential=2.5e9,
        elimination_ea_kj_mol=74.0,
        scale=260.0 * (0.92 + 0.22 * aw_norm),
        # Wave B8: labels the two dicarbonyl-pool pairs above as
        # source_status='no_verifiable_source'. No value is changed.
        provenance_consumer="src.safety._estimate_dicarbonyl_pools",
    )

    high_temp_drive = _sigmoid(thermal_temp_c, 138.0, 8.0)
    go_share = _clamp(0.62 + 0.20 * aw_norm - 0.34 * high_temp_drive, 0.12, 0.82)
    mgo_share = _clamp(0.28 + 0.44 * high_temp_drive - 0.20 * aw_norm, 0.12, 0.82)
    total_share = go_share + mgo_share
    if total_share > 0.95:
        scale = 0.95 / total_share
        go_share *= scale
        mgo_share *= scale

    polyphenol_window = _clamp(polyphenol_factor, 0.0, 1.0)
    go_pool = amadori_pool * go_share * (1.0 - 0.10 * polyphenol_window)
    mgo_pool = amadori_pool * mgo_share * (1.0 - 0.247 * polyphenol_window)
    return amadori_pool, go_pool, mgo_pool


def predict_cml(
    lysine_mM: float,
    reducing_sugar_mM: float,
    temp_C: float,
    time_min: float,
    *,
    water_activity: Optional[float] = None,
    effective_temp_c: Optional[float] = None,
    polyphenol_factor: float = 0.0,
    soy_isoflavone_factor: float = 0.0,
) -> float:
    """Predict N-epsilon-(carboxymethyl)lysine.

    Units: inputs in mM / degC / min; RETURN VALUE IS ppb (ug/kg of food), per
    the module `UNITS` contract. Known miscalibration (2026-08-27): this
    predictor's empirical constants were tuned against benchmark files that held
    mg/kg values in ppb fields, so it currently sits ~1.2e3x below the commercial
    PBMA CML band (16.5-47.6 mg/kg = 1.65e4-4.76e4 ppb). Reported, not rescaled.
    """
    thermal_temp_c = _resolve_effective_temp_c(temp_C, effective_temp_c)
    aw = 0.55 if water_activity is None else _clamp(water_activity, 0.05, 0.98)
    aw_norm = _clamp((aw - 0.25) / 0.65, 0.0, 1.0)
    _amadori_pool, go_pool, _mgo_pool = _estimate_dicarbonyl_pools(
        lysine_mM=lysine_mM,
        reducing_sugar_mM=reducing_sugar_mM,
        temp_C=temp_C,
        time_min=time_min,
        water_activity=water_activity,
        effective_temp_c=effective_temp_c,
        polyphenol_factor=polyphenol_factor,
    )
    formation_window = _sigmoid(thermal_temp_c, 122.0, 10.0)
    degradation_window = _sigmoid(thermal_temp_c, 162.0, 7.0)
    moisture_gain = 0.80 + 0.35 * aw_norm
    isoflavone_guardrail = 1.0 - 0.55 * _clamp(soy_isoflavone_factor, 0.0, 1.0)
    return max(0.0, go_pool * 0.40 * formation_window * moisture_gain * (1.0 - 0.45 * degradation_window) * isoflavone_guardrail)


def predict_cel(
    lysine_mM: float,
    reducing_sugar_mM: float,
    temp_C: float,
    time_min: float,
    *,
    water_activity: Optional[float] = None,
    effective_temp_c: Optional[float] = None,
    polyphenol_factor: float = 0.0,
    soy_isoflavone_factor: float = 0.0,
) -> float:
    """Predict N-epsilon-(carboxyethyl)lysine.

    Units: inputs in mM / degC / min; RETURN VALUE IS ppb (ug/kg of food), per
    the module `UNITS` contract. Known miscalibration (2026-08-27): ~1.1e3x below
    the commercial PBMA CEL band (25.2-86.2 mg/kg = 2.52e4-8.62e4 ppb) for the
    same reason as `predict_cml`. Reported, not rescaled.
    """
    thermal_temp_c = _resolve_effective_temp_c(temp_C, effective_temp_c)
    aw = 0.45 if water_activity is None else _clamp(water_activity, 0.05, 0.98)
    aw_norm = _clamp((aw - 0.25) / 0.65, 0.0, 1.0)
    _amadori_pool, _go_pool, mgo_pool = _estimate_dicarbonyl_pools(
        lysine_mM=lysine_mM,
        reducing_sugar_mM=reducing_sugar_mM,
        temp_C=temp_C,
        time_min=time_min,
        water_activity=water_activity,
        effective_temp_c=effective_temp_c,
        polyphenol_factor=polyphenol_factor,
    )
    formation_window = _sigmoid(thermal_temp_c, 132.0, 9.0)
    high_temp_gain = 0.90 + 0.30 * _sigmoid(thermal_temp_c, 142.0, 8.0)
    moisture_penalty = 1.12 - 0.40 * aw_norm
    isoflavone_guardrail = 1.0 - 0.65 * _clamp(soy_isoflavone_factor, 0.0, 1.0)
    return max(0.0, mgo_pool * 0.44 * formation_window * high_temp_gain * moisture_penalty * isoflavone_guardrail)


def predict_furosine(
    temp_C: float,
    time_min: float,
    *,
    lysine_mM: Optional[float] = None,
    reducing_sugar_mM: Optional[float] = None,
    protein_type: str = "free",
    water_activity: Optional[float] = None,
    effective_temp_c: Optional[float] = None,
) -> float:
    """Predict furosine.

    Units: inputs in mM / degC / min; RETURN VALUE IS ppb (ug/kg of food), per
    the module `UNITS` contract. Known miscalibration (2026-08-27): the
    Ramirez-Jimenez mild-extrusion peak is 8.7 mg per 100 g protein = 17400 ppb
    on the declared 20% protein basis (`mg_per_100g_protein_to_ppb`), so this
    predictor sits ~2.0e2x below its own anchor. Reported, not rescaled.
    """
    if time_min <= 0.0:
        return 0.0

    thermal_temp_c = _resolve_effective_temp_c(temp_C, effective_temp_c)
    aw = 0.55 if water_activity is None else _clamp(water_activity, 0.05, 0.98)
    aw_norm = _clamp((aw - 0.25) / 0.65, 0.0, 1.0)
    protein_key = str(protein_type).strip().lower()
    protein_factor = {
        "free": 0.80,
        "pea_conc": 0.95,
        "pea_iso": 1.00,
        "soy_conc": 1.05,
        "soy_iso": 1.12,
        "myco": 0.75,
    }.get(protein_key, 1.0)
    if lysine_mM is not None and reducing_sugar_mM is not None and lysine_mM > 0.0 and reducing_sugar_mM > 0.0:
        amadori_pool, _go_pool, _mgo_pool = _estimate_dicarbonyl_pools(
            lysine_mM=lysine_mM,
            reducing_sugar_mM=reducing_sugar_mM,
            temp_C=temp_C,
            time_min=time_min,
            water_activity=water_activity,
            effective_temp_c=effective_temp_c,
        )
        moisture_gain = 0.94 + 0.12 * aw_norm
        return max(0.0, amadori_pool * 0.435 * protein_factor * moisture_gain)

    precursor_drive = protein_factor * (0.95 + 0.18 * aw_norm)
    return _formation_elimination_signal(
        precursor_drive,
        temp_c=thermal_temp_c,
        time_min=time_min,
        formation_pre_exponential=8.0e5,
        formation_ea_kj_mol=50.0,
        elimination_pre_exponential=1.8e10,
        elimination_ea_kj_mol=73.0,
        scale=120.0,
        # Wave B8: labels the two furosine pairs above as
        # source_status='no_verifiable_source'. No value is changed.
        provenance_consumer="src.safety.predict_furosine",
    )

#: pH response shape constants for acrylamide formation (declared; see the
#: comment in `predict_acrylamide`). The decline term is normalised at pH 6.0 so
#: the previously calibrated pH-6 behaviour is bit-for-bit preserved and only the
#: alkaline tail changes.
ACRYLAMIDE_PH_DECLINE_MIDPOINT = 7.6
ACRYLAMIDE_PH_DECLINE_SLOPE = 1.0
ACRYLAMIDE_PH_NORMALIZATION_POINT = 6.0
#: Provenance for the pH response: this payload entry was previously loaded but
#: never consumed anywhere in the module.
ACRYLAMIDE_PH_REFERENCE_ID = (
    "de_vleeschouwer_2006_acrylamide_aqueous"
    if get_safety_reference_payload("de_vleeschouwer_2006_acrylamide_aqueous")
    else ""
)


def _acrylamide_ph_decline(pH: float) -> float:
    exponent = _clamp(ACRYLAMIDE_PH_DECLINE_SLOPE * (float(pH) - ACRYLAMIDE_PH_DECLINE_MIDPOINT), -60.0, 60.0)
    return 1.0 / (1.0 + 10.0 ** exponent)


def _acrylamide_ph_factor(pH: float) -> float:
    """Bounded pH response: rises with amine deprotonation, peaks near pH 8, falls above.

    Below the normalisation point (pH 6.0) the decline term is capped at 1.0 so
    the acid side is governed by deprotonation alone, exactly as before.
    """
    nucleophilicity = 1.0 / (1.0 + 10.0 ** _clamp(8.8 - float(pH), -60.0, 60.0))
    decline = _acrylamide_ph_decline(pH) / _acrylamide_ph_decline(ACRYLAMIDE_PH_NORMALIZATION_POINT)
    return float(nucleophilicity * min(1.0, decline))


def predict_acrylamide(
    asparagine_mM: float,
    reducing_sugar_mM: float,
    temp_C: float,
    time_min: float,
    pH: float = 6.0,
    ea_modifier_kcal: float = 0.0,
    water_activity: Optional[float] = None,
    moisture_regime: Optional[str] = None,
    effective_temp_c: Optional[float] = None,
) -> SafetyResult:
    """
    Implements formation-elimination kinetics for acrylamide (Knol/Parker model).
    
    Kinetics:
    dA/dt = kf * [Asn] * [Sugar] - ke * [A]
    
    Analytical Solution:
    A(t) = (kf * [Asn] * [Sugar] / ke) * (1 - exp(-ke * t))
    
    Units: asparagine_mM / reducing_sugar_mM in mM, temp_C in degC, time_min in
    minutes; the returned acrylamide_ppb and uncertainty_ppb are ppb (ug/kg of
    food). See the module-level `UNITS` contract.

    References:
    - Knol et al. 2009 (Formation kinetics)
    - Parker et al. 2012 (Elimination/Degradation kinetics)
    - De Vleeschouwer et al. 2006, JAFC 54:7847 (pH-resolved formation/elimination)
    """
    if asparagine_mM <= 0 or reducing_sugar_mM <= 0:
        # NOT a clean zero: the analyte could not be assessed at all.
        return SafetyResult(0.0, 0.0, False, "Not assessed (no asparagine and/or reducing sugar)", assessed=False)

    # Constants
    R = 8.314 # J/mol/K
    thermal_temp_c = temp_C if effective_temp_c is None else float(effective_temp_c)
    T_K = thermal_temp_c + 273.15
    time_sec = time_min * 60.0
    MW_AA = 71.08
    
    # 1. Formation (kf)
    # Base Ea for formation from Knol 2009 (~129 kJ/mol)
    Ea_f = 129000.0 + (ea_modifier_kcal * 4184.0)
    A_f = 1.6e13 # L/mol/s 
    
    # pH effect on formation: Asparagine amine nucleophilicity (pKa ~8.8)
    # Most reactive in slightly alkaline, sharply drops below pH 5.
    # 2026-08-27: the bare deprotonation sigmoid rose monotonically with pH, so
    # the model kept gaining acrylamide at pH 9-10 where real systems lose it.
    # A high-pH decline term is applied so the response PEAKS near pH 8 and falls
    # away above it. Provenance: de_vleeschouwer_2006_acrylamide_aqueous in
    # data/lit/safety_reference_payloads.json establishes that both formation and
    # elimination rate constants are pH-resolved bands in this region; that entry
    # pins the EXISTENCE and direction of the effect, not the optimum, so the
    # decline midpoint (7.6) and slope are declared shape constants.
    f_ph = _acrylamide_ph_factor(pH)
    moisture_factor = _acrylamide_moisture_factor(water_activity, moisture_regime)
    kf = A_f * math.exp(-Ea_f / (R * T_K)) * f_ph * moisture_factor
    
    # 2. Elimination (ke)
    # Acrylamide degrades at high T. Ea_e typically ~90-110 kJ/mol
    Ea_e = 105000.0 
    A_e = 5.0e8 # s^-1
    ke = A_e * math.exp(-Ea_e / (R * T_K))
    
    # 3. Resolve Concentrations (mM to mol/L)
    asn_molar = asparagine_mM / 1000.0
    sugar_molar = reducing_sugar_mM / 1000.0
    
    # Avoid division by zero if ke is extremely small
    if ke < 1e-12:
        # Fallback to simple first-order formation
        aa_molar = kf * asn_molar * sugar_molar * time_sec
    else:
        # Analytical solution
        aa_molar = (kf * asn_molar * sugar_molar / ke) * (1 - math.exp(-ke * time_sec))
        
    # Convert mol/L -> ppm (mg/kg assuming rho=1) -> ppb
    aa_ppm = aa_molar * MW_AA * 1000.0
    aa_ppb = aa_ppm * 1000.0
    
    # Heuristic uncertainty (25% based on literature variability)
    unc = aa_ppb * 0.25

    # Band-derived flag (2026-08-27). The old 100-ppb literal was described as an
    # "EU benchmark for meat analogues"; no such category exists in Reg. (EU)
    # 2017/2158 Annex IV, and the literal was never traceable to anything. The
    # threshold is now the lowest acrylamide level actually reported in a
    # plant-protein product by the reference payload
    # (ACRYLAMIDE_ATTENTION_PPB), and the prediction's own uncertainty edge is
    # what is compared against it.
    exceeds_band = (aa_ppb + unc) >= ACRYLAMIDE_ATTENTION_PPB

    return SafetyResult(
        acrylamide_ppb=aa_ppb,
        uncertainty_ppb=unc,
        flagged=exceeds_band,
        description=(
            f"At or above the finished-product reference band ({ACRYLAMIDE_ATTENTION_PPB:.1f} ppb)"
            if exceeds_band
            else f"Below the finished-product reference band ({ACRYLAMIDE_ATTENTION_PPB:.1f} ppb)"
        ),
        assessed=True,
    )

def _band_severity(value: float, *, attention: float, action: float) -> float:
    """Map a level onto [0, 1] against a reference band.

    Shape (declared): linear ramp 0 -> 0.35 from zero up to the attention level
    (the lowest level seen in real products), then a logarithmic climb 0.35 -> 1.0
    from the attention level to the action level (the top of the observed
    industrial band), saturating at 1.0 above it. Monotone everywhere, exactly
    0.0 only at 0, and bounded for any input - a 1e9-fold span in the input still
    lands inside [0, 1].
    """
    level = max(0.0, float(value))
    if level <= 0.0:
        return 0.0
    attention = max(float(attention), 1e-12)
    action = max(float(action), attention * 1.0000001)
    if level < attention:
        return 0.35 * (level / attention)
    span = math.log10(action / attention)
    progress = _clamp(math.log10(level / attention) / span, 0.0, 1.0)
    return 0.35 + 0.65 * progress


def _combine_risk(sub_scores: List[float]) -> float:
    """Severity-dominant aggregation, guaranteed to stay in [0, 1].

    The worst analyte sets the floor; co-occurring analytes add a bounded
    crowding term that can only consume a quarter of the remaining headroom.
    """
    scores = [_clamp(score, 0.0, 1.0) for score in sub_scores if score > 0.0]
    if not scores:
        return 0.0
    worst = max(scores)
    if len(scores) > 1:
        remainder = sorted(scores, reverse=True)[1:]
        crowding = sum(remainder) / len(remainder)
        return _clamp(worst + (1.0 - worst) * 0.25 * crowding, 0.0, 1.0)
    return _clamp(worst, 0.0, 1.0)


def _thermal_exposure_active(
    temp_C: float,
    time_min: float,
    extrusion_process: Dict[str, Any],
) -> bool:
    """True when the formulation actually saw a damage-relevant thermal history."""
    if float(time_min) <= 0.0:
        return False
    temperatures = [float(temp_C)]
    effective = extrusion_process.get("effective_temperature_celsius")
    if effective is not None:
        temperatures.append(float(effective))
    sterilization = extrusion_process.get("sterilization", {})
    if isinstance(sterilization, dict) and sterilization.get("enabled"):
        sterilization_temp = sterilization.get("temperature_celsius")
        if sterilization_temp is not None:
            temperatures.append(float(sterilization_temp))
    return max(temperatures) >= DAMAGE_MARKER_THERMAL_ONSET_C


def evaluate_formulation_safety(
    precursors: Dict[str, float],
    temp_C: float,
    time_min: float,
    pH: float,
    modifiers: Optional[Dict[str, Any]] = None
) -> Tuple[float, List[str]]:
    """Aggregated safety score and flagged toxins.

    Units: `precursors` values in mM, `temp_C` in degC, `time_min` in minutes.
    RETURNS a genuine [0, 1] risk (1.0 = at or beyond the top of every reference
    band, 0.0 = nothing assessed above zero) plus the list of flagged analytes.

    A 0.0 with an empty flag list means NOTHING WAS ASSESSED (no asparagine, no
    reducing sugar, no thermal exposure) - it is not a clean bill of health. The
    per-analyte assessment status is surfaced by `SafetyResult.assessed` and by
    the family-12 runtime lane's `acrylamide_assessment_status`.

    2026-08-27: the old aggregate was `log10(ppb / 1e-20) / 10` summed across
    analytes - unbounded, routinely outside [0, 1] despite its docstring, and so
    log-compressed that a 1e9-fold acrylamide difference moved it by 0.9. Each
    analyte is now scored against the reference bands loaded from
    data/lit/safety_reference_payloads.json (see `_band_severity`).
    """
    sub_scores: List[float] = []
    flagged: List[str] = []
    mods = modifiers or {}
    extrusion_process = mods.get("__extrusion_process__", {}) if isinstance(mods.get("__extrusion_process__", {}), dict) else {}

    # 1. Acrylamide Check
    asn_conc = 0.0
    sugar_conc = 0.0
    # Sucrose is deliberately NOT in the reducing-sugar list: it is a
    # non-reducing disaccharide, and substituting sucrose for glucose/fructose is
    # itself a standard acrylamide mitigation. Counting it as a reducing sugar
    # penalised the mitigation. Lactose and maltose stay: both are reducing.
    reducing_sugars = ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar", "lactose"]
    for name, conc in precursors.items():
        try:
            value = float(conc)
        except (TypeError, ValueError):
            continue
        # Word-boundary matching: "asn" used to match "asnase" (asparaginase).
        if _has_word(name, "asparagine") or _has_word(name, "asn"):
            asn_conc = value
        if _has_any_word(name, reducing_sugars):
            sugar_conc += value

    if asn_conc > 0 and sugar_conc > 0:
        # Resolve Ea modifier for Acrylamide. Modifier values are not guaranteed
        # numeric (runtime context dicts live in the same mapping under dunder
        # keys), so a non-numeric hit is skipped instead of crashing with a
        # TypeError deep inside the Arrhenius term.
        ea_mod = 0.0
        for k, v in mods.items():
            key = str(k)
            if key.startswith("__") or "acrylamide" not in key.lower():
                continue
            try:
                ea_mod = float(v)
            except (TypeError, ValueError):
                warnings.warn(
                    f"safety: ignoring non-numeric acrylamide Ea modifier {key!r}={v!r}; "
                    "expected a float in kcal/mol.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                ea_mod = 0.0
            break

        aa_res = predict_acrylamide(
            asn_conc,
            sugar_conc,
            temp_C,
            time_min,
            pH,
            ea_mod,
            water_activity=extrusion_process.get("water_activity"),
            moisture_regime=extrusion_process.get("moisture_regime"),
            effective_temp_c=extrusion_process.get("effective_temperature_celsius"),
        )
        runtime_adjustment = _resolve_acrylamide_runtime_adjustment(precursors, mods)
        mitigation = float(runtime_adjustment.get("mitigation_fraction", 0.0))
        effective_acrylamide_ppb = (
            aa_res.acrylamide_ppb
            * (1.0 - mitigation)
            * float(runtime_adjustment.get("process_modifier", 1.0))
        )
        effective_uncertainty_ppb = aa_res.uncertainty_ppb * (1.0 - mitigation)
        if aa_res.assessed and effective_acrylamide_ppb >= ACRYLAMIDE_REPORTING_FLOOR_PPB:
            _append_unique(flagged, "Acrylamide")
            # SafetyResult.flagged is now consumed (it used to be dead code):
            # re-evaluated on the mitigated level, it marks band exceedance.
            if (effective_acrylamide_ppb + effective_uncertainty_ppb) >= ACRYLAMIDE_ATTENTION_PPB:
                _append_unique(flagged, "Acrylamide above reference band")
        if aa_res.assessed:
            sub_scores.append(
                _band_severity(
                    effective_acrylamide_ppb,
                    attention=ACRYLAMIDE_ATTENTION_PPB,
                    action=ACRYLAMIDE_ACTION_PPB,
                )
            )

    total_damage = extrusion_process.get("total_damage_load", {}) if isinstance(extrusion_process, dict) else {}
    furosine = float(total_damage.get("furosine_mg_per_kg", 0.0) or 0.0)
    lal = float(total_damage.get("lal_mg_per_kg", 0.0) or 0.0)

    # Flag hygiene (2026-08-27): the damage load carried by an extrusion profile
    # includes an ingredient-inherited baseline that src/extrusion.py attaches
    # regardless of temperature or time, so `> 0.0` flagged Furosine and LAL for a
    # 25 C / 0 min formulation. Flags now require actual thermal exposure. The
    # ROOT CAUSE (unconditional `pre_extrusion_damage_load` + `active: True` for
    # any lme/hme regime) lives in src/extrusion.py and is reported separately.
    thermal_exposure = _thermal_exposure_active(temp_C, time_min, extrusion_process)

    if furosine > 0.0 and thermal_exposure:
        _append_unique(flagged, "Furosine")
        sub_scores.append(
            _band_severity(
                furosine,
                attention=FUROSINE_ATTENTION_MG_PER_KG,
                action=FUROSINE_ACTION_MG_PER_KG,
            )
        )

    if lal > 0.0 and thermal_exposure:
        _append_unique(flagged, "LAL")
        sub_scores.append(
            _band_severity(
                lal,
                attention=LAL_PROVISIONAL_ATTENTION_MG_PER_KG,
                action=LAL_PROVISIONAL_ACTION_MG_PER_KG,
            )
        )

    return _combine_risk(sub_scores), flagged
