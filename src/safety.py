"""
src/safety.py

Quantitative safety modeling for Maillard by-products.
Focus: Acrylamide formation from Asparagine + Reducing Sugars.

Reference:
- Knol et al. 2009 J. Ag. Food Chem. (Kinetic model for acrylamide)
- Stadler et al. 2004 (Asparagine involvement)
"""

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional

from src.extrusion import normalize_moisture_regime


ROOT = Path(__file__).resolve().parents[1]
SAFETY_REFERENCE_PAYLOAD_PATH = ROOT / "data" / "lit" / "safety_reference_payloads.json"


def _load_safety_reference_payloads() -> dict:
    with open(SAFETY_REFERENCE_PAYLOAD_PATH, "r", encoding="utf-8") as handle:
        return json.load(handle)


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


def _analyte_matches(entry_analyte: str, requested_analyte: Optional[str]) -> bool:
    if not requested_analyte:
        return True
    normalized_entry = str(entry_analyte).strip().lower()
    normalized_requested = str(requested_analyte).strip().lower()
    if normalized_entry == normalized_requested:
        return True
    tokens = {
        token.strip()
        for token in normalized_entry.replace("+", "/").split("/")
        if token.strip()
    }
    return normalized_requested in tokens


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
        entry_analyte = str(entry.get("analyte", "")).strip().lower()
        if not _analyte_matches(entry_analyte, analyte):
            continue
        entry_visibility = _normalize_safety_visibility(entry)
        if requested_visibility == "default" and entry_visibility != "default":
            continue
        entries.append(dict(entry, report_visibility=entry_visibility))
    return entries


def build_safety_reference_context(*, analyte: str = "acrylamide") -> Dict[str, List[Dict[str, Any]]]:
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
    return {
        "default_entries": default_entries,
        "extended_entries": extended_entries,
    }


def get_safety_reference_range(matrix_family: str, reference_id: str = "squeo_2023_pbpi_acrylamide") -> Optional[dict]:
    payload = get_safety_reference_payload(reference_id)
    if payload is None:
        return None
    for item in payload.get("matrix_reference_ranges", []):
        if str(item.get("matrix_family", "")) == matrix_family:
            return item
    return None


def _normalize_runtime_tokens(values: Any) -> List[str]:
    if not isinstance(values, list):
        return []
    return [str(value).strip().lower() for value in values if str(value).strip()]


def _resolve_acrylamide_runtime_adjustment(
    precursors: Dict[str, float],
    modifiers: Dict[str, Any],
) -> Dict[str, Any]:
    runtime_context = modifiers.get("__runtime_context__", {}) if isinstance(modifiers.get("__runtime_context__", {}), dict) else {}
    additives = _normalize_runtime_tokens(runtime_context.get("additives", []))
    interventions = _normalize_runtime_tokens(runtime_context.get("interventions", []))

    precursor_names = [str(name).strip().lower() for name in precursors]
    cys_present = any("cysteine" in name or name == "cys" for name in precursor_names)
    gly_present = any("glycine" in name or name == "gly" for name in precursor_names)
    saponin_active = any("saponin" in token or "quillaja" in token for token in additives + interventions)

    mitigation_fraction = 0.0
    reference_ids: List[str] = []
    if cys_present and gly_present:
        mitigation_fraction = 0.65
        reference_ids.append("pmc_12648097_acrylamide_mitigation")
    elif cys_present or gly_present:
        mitigation_fraction = 0.20
        reference_ids.append("pmc_12648097_acrylamide_mitigation")

    process_modifier = 1.0
    if saponin_active:
        process_modifier = 0.88
        reference_ids.append("kocadagli_2016_saponin_acrylamide_modifier")

    return {
        "mitigation_fraction": float(mitigation_fraction),
        "process_modifier": float(process_modifier),
        "reference_ids": reference_ids,
    }

@dataclass
class SafetyResult:
    acrylamide_ppb: float
    uncertainty_ppb: float
    flagged: bool
    description: str


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
) -> float:
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
    )

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
    
    References:
    - Knol et al. 2009 (Formation kinetics)
    - Parker et al. 2012 (Elimination/Degradation kinetics)
    """
    if asparagine_mM <= 0 or reducing_sugar_mM <= 0:
        return SafetyResult(0.0, 0.0, False, "No precursors")

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
    # Most reactive in slightly alkaline, sharply drops below pH 5
    f_ph = 1.0 / (1.0 + 10**(8.8 - pH))
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
    
    # EU benchmark for meat analogues/cereals is often ~300-750 ppb
    # We use a conservative 100 ppb threshold for flagging
    is_safe = aa_ppb < 100.0
    
    return SafetyResult(
        acrylamide_ppb=aa_ppb,
        uncertainty_ppb=unc,
        flagged=not is_safe,
        description="High acrylamide risk" if not is_safe else "Normal levels"
    )

def evaluate_formulation_safety(
    precursors: Dict[str, float],
    temp_C: float,
    time_min: float,
    pH: float,
    modifiers: Optional[Dict[str, Any]] = None
) -> Tuple[float, List[str]]:
    """
    Aggregated safety score and flagged toxins.
    1.0 = Max danger, 0.0 = Safe (though we don't cap in scientific mode)
    """
    total_risk = 0.0
    flagged = []
    mods = modifiers or {}
    extrusion_process = mods.get("__extrusion_process__", {}) if isinstance(mods.get("__extrusion_process__", {}), dict) else {}
    
    # 1. Acrylamide Check
    asn_conc = 0.0
    sugar_conc = 0.0
    for name, conc in precursors.items():
        n_low = name.lower()
        if "asparagine" in n_low or "asn" in n_low:
            asn_conc = conc
        if any(s in n_low for s in ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar", "sucrose", "lactose"]):
            sugar_conc += conc
            
    if asn_conc > 0 and sugar_conc > 0:
        # Resolve Ea modifier for Acrylamide
        ea_mod = 0.0
        for k, v in mods.items():
            if "acrylamide" in k.lower():
                ea_mod = v
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
        effective_acrylamide_ppb = (
            aa_res.acrylamide_ppb
            * (1.0 - float(runtime_adjustment.get("mitigation_fraction", 0.0)))
            * float(runtime_adjustment.get("process_modifier", 1.0))
        )
        # Threshold for detection - ensure it's high enough to be seen but low enough to catch precursors
        if effective_acrylamide_ppb > 1e-25:
            _append_unique(flagged, "Acrylamide")
            # Logarithmic risk scaling: ensure we don't saturate for small differences
            # log10(1e-15 / 1e-20) / 10 = 0.5
            # log10(1.2e-16 / 1e-20) / 10 = 0.4
            risk_raw = math.log10(effective_acrylamide_ppb / 1e-20) / 10.0
            total_risk += max(0.01, risk_raw)

    total_damage = extrusion_process.get("total_damage_load", {}) if isinstance(extrusion_process, dict) else {}
    furosine = float(total_damage.get("furosine_mg_per_kg", 0.0) or 0.0)
    lal = float(total_damage.get("lal_mg_per_kg", 0.0) or 0.0)

    if furosine > 0.0:
        _append_unique(flagged, "Furosine")
        total_risk += min(1.0, furosine / 40.0)

    if lal > 0.0:
        _append_unique(flagged, "LAL")
        total_risk += min(1.25, lal / 120.0)
            
    return total_risk, flagged
