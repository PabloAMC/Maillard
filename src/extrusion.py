from __future__ import annotations

import math
from typing import Dict, List, Mapping, Optional, Sequence


_MOISTURE_REGIME_ALIASES = {
    "lme": "lme",
    "low": "lme",
    "low_moisture": "lme",
    "low-moisture": "lme",
    "low moisture": "lme",
    "hme": "hme",
    "high": "hme",
    "high_moisture": "hme",
    "high-moisture": "hme",
    "high moisture": "hme",
}


_PRE_EXTRUSION_DAMAGE_BASELINES = {
    "free": {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0},
    "pea_conc": {"furosine_mg_per_kg": 12.0, "lal_mg_per_kg": 42.0},
    "pea_iso": {"furosine_mg_per_kg": 18.0, "lal_mg_per_kg": 68.0},
    "soy_conc": {"furosine_mg_per_kg": 15.0, "lal_mg_per_kg": 54.0},
    "soy_iso": {"furosine_mg_per_kg": 24.0, "lal_mg_per_kg": 88.0},
    "myco": {"furosine_mg_per_kg": 8.0, "lal_mg_per_kg": 26.0},
}


_DEFAULT_TRANSPORT_PANEL = [
    "2-methyl-3-furanthiol",
    "2-furfurylthiol",
    "Hexanal",
    "2-Pentylfuran",
    "2,5-Dimethylpyrazine",
]


_VOLATILE_CLASS_PARAMETERS = {
    "sulfur": {
        "shear_release_gain": 0.22,
        "aggregation_penalty": 0.18,
        "aggregation_center_c": 152.0,
        "max_die_stripping_fraction": 0.24,
        "die_flash_center_c": 78.0,
        "relative_vapor_pressure": 0.92,
    },
    "aldehyde": {
        "shear_release_gain": 0.14,
        "aggregation_penalty": 0.10,
        "aggregation_center_c": 145.0,
        "max_die_stripping_fraction": 0.40,
        "die_flash_center_c": 72.0,
        "relative_vapor_pressure": 1.00,
    },
    "furan": {
        "shear_release_gain": 0.11,
        "aggregation_penalty": 0.08,
        "aggregation_center_c": 148.0,
        "max_die_stripping_fraction": 0.26,
        "die_flash_center_c": 78.0,
        "relative_vapor_pressure": 0.78,
    },
    "pyrazine": {
        "shear_release_gain": 0.18,
        "aggregation_penalty": 0.11,
        "aggregation_center_c": 158.0,
        "max_die_stripping_fraction": 0.10,
        "die_flash_center_c": 92.0,
        "relative_vapor_pressure": 0.36,
    },
    "alcohol": {
        "shear_release_gain": 0.10,
        "aggregation_penalty": 0.07,
        "aggregation_center_c": 145.0,
        "max_die_stripping_fraction": 0.20,
        "die_flash_center_c": 76.0,
        "relative_vapor_pressure": 0.72,
    },
    "unknown": {
        "shear_release_gain": 0.09,
        "aggregation_penalty": 0.07,
        "aggregation_center_c": 150.0,
        "max_die_stripping_fraction": 0.16,
        "die_flash_center_c": 82.0,
        "relative_vapor_pressure": 0.55,
    },
}


def _sigmoid(x: float, center: float, k: float) -> float:
    try:
        return 1.0 / (1.0 + math.exp(-k * (float(x) - float(center))))
    except OverflowError:
        return 0.0 if x < center else 1.0


def classify_volatile_compound(name: str) -> str:
    compound = str(name).strip().lower()
    if any(token in compound for token in ["thiol", "sulfide", "mercapt", "mft", "fft"]):
        return "sulfur"
    if any(token in compound for token in ["anal", "enal"]):
        return "aldehyde"
    if "pyrazine" in compound or "thiazole" in compound:
        return "pyrazine"
    if "furan" in compound:
        return "furan"
    if any(token in compound for token in ["hexanol", "octen-3-ol", " alcohol", "-ol"]):
        return "alcohol"
    return "unknown"


def estimate_die_exit_temperature_celsius(
    effective_temperature_celsius: float,
    moisture_regime: Optional[str],
    *,
    die_exit_temperature_celsius: Optional[float] = None,
) -> float:
    if die_exit_temperature_celsius is not None:
        return float(die_exit_temperature_celsius)

    regime = normalize_moisture_regime(moisture_regime)
    effective_temp = float(effective_temperature_celsius)
    cooling_drop = 35.0 if regime == "hme" else 12.0
    lower_bound = 60.0 if regime == "hme" else 90.0
    return max(lower_bound, min(effective_temp, effective_temp - cooling_drop))


def estimate_mean_residence_time_seconds(
    *,
    sme_kj_per_kg: float,
    moisture_regime: Optional[str],
    screw_speed_rpm: Optional[float] = None,
    feed_rate_kg_per_h: Optional[float] = None,
) -> float:
    regime = normalize_moisture_regime(moisture_regime)
    base_seconds = 42.0 if regime == "hme" else 28.0
    rpm = None if screw_speed_rpm is None else max(50.0, float(screw_speed_rpm))
    feed = None if feed_rate_kg_per_h is None else max(0.5, float(feed_rate_kg_per_h))
    speed_factor = 1.0 if rpm is None else max(0.70, min(1.35, 280.0 / rpm))
    feed_factor = 1.0 if feed is None else max(0.75, min(1.35, 4.6 / feed))
    sme_factor = max(0.85, min(1.10, 1.02 - 0.0005 * max(0.0, float(sme_kj_per_kg) - 120.0)))
    return max(12.0, min(90.0, base_seconds * speed_factor * feed_factor * sme_factor))


def estimate_relative_rtd_spread(
    *,
    sme_kj_per_kg: float,
    moisture_regime: Optional[str],
    screw_speed_rpm: Optional[float] = None,
) -> float:
    regime = normalize_moisture_regime(moisture_regime)
    base_spread = 0.24 if regime == "hme" else 0.30
    sme_broadening = (0.03 if regime == "hme" else 0.06) * _sigmoid(float(sme_kj_per_kg), 180.0, 0.03)
    rpm_compaction = 0.0
    if screw_speed_rpm is not None:
        rpm_compaction = 0.03 * _sigmoid(float(screw_speed_rpm), 260.0, 0.04)
    return max(0.15, min(0.45, base_spread + sme_broadening - rpm_compaction))


def evaluate_extrusion_rtd_need(extrusion_process: Mapping[str, object]) -> Dict[str, object]:
    regime = normalize_moisture_regime(str(extrusion_process.get("moisture_regime", "hme")))
    sme = float(extrusion_process.get("sme_kj_per_kg", 0.0) or 0.0)
    spread = float(extrusion_process.get("relative_rtd_spread", 0.0) or 0.0)
    mean_residence = float(extrusion_process.get("mean_residence_time_seconds", 0.0) or 0.0)

    if regime == "lme" and (sme >= 180.0 or spread >= 0.33):
        decision = "simple_rtd_recommended"
        reason = "Low-moisture or broad-distribution cases can distort time-at-temperature enough that the sequential-zone surrogate is no longer the best scientist-facing summary."
    else:
        decision = "sequential_zone_sufficient_for_current_use_case"
        reason = "For the current HME scientist workflow, the sequential-zone model already captures the dominant thermal history while screw-geometry-free RTD would add parameters the repo cannot yet benchmark."

    return {
        "decision": decision,
        "reason": reason,
        "mean_residence_time_seconds": mean_residence,
        "relative_spread": spread,
    }


def compute_extrusion_headspace_adjustment(
    compound_name: str,
    extrusion_process: Mapping[str, object],
) -> Dict[str, float | str]:
    compound_class = classify_volatile_compound(compound_name)
    params = _VOLATILE_CLASS_PARAMETERS.get(compound_class, _VOLATILE_CLASS_PARAMETERS["unknown"])
    regime = normalize_moisture_regime(str(extrusion_process.get("moisture_regime", "hme")), float(extrusion_process.get("water_activity", 0.8) or 0.8))
    sme = float(extrusion_process.get("sme_kj_per_kg", 0.0) or 0.0)
    water_activity = float(extrusion_process.get("water_activity", 0.8) or 0.8)
    effective_temp = float(extrusion_process.get("effective_temperature_celsius", extrusion_process.get("jacket_temperature_celsius", 120.0)) or 120.0)
    die_temp = estimate_die_exit_temperature_celsius(
        effective_temp,
        regime,
        die_exit_temperature_celsius=extrusion_process.get("die_exit_temperature_celsius"),
    )

    sme_term = max(0.0, min(1.0, sme / 180.0))
    release_gain = float(params["shear_release_gain"]) * sme_term * (1.0 if regime == "hme" else 0.82)
    aggregation_penalty = float(params["aggregation_penalty"]) * _sigmoid(effective_temp, float(params["aggregation_center_c"]), 0.09)
    if regime == "lme":
        aggregation_penalty *= 1.15
    if compound_class == "sulfur":
        aggregation_penalty *= 1.10
    shear_release_factor = max(0.75, min(1.45, 1.0 + release_gain - aggregation_penalty))

    moisture_flash_factor = 0.80 + 0.45 * max(0.0, min(1.0, (water_activity - 0.45) / 0.35))
    die_flash_activation = _sigmoid(die_temp, float(params["die_flash_center_c"]), 0.09)
    die_stripping_fraction = min(
        float(params["max_die_stripping_fraction"]),
        float(params["max_die_stripping_fraction"]) * die_flash_activation * moisture_flash_factor,
    )
    combined_factor = max(0.30, min(1.35, shear_release_factor * (1.0 - die_stripping_fraction)))

    return {
        "compound": str(compound_name),
        "compound_class": compound_class,
        "relative_vapor_pressure_proxy": float(params["relative_vapor_pressure"]),
        "shear_release_factor": float(shear_release_factor),
        "die_stripping_fraction": float(die_stripping_fraction),
        "combined_headspace_factor": float(combined_factor),
    }


def _build_transport_summary(extrusion_process: Mapping[str, object]) -> Dict[str, object]:
    panel = {
        compound: compute_extrusion_headspace_adjustment(compound, extrusion_process)
        for compound in _DEFAULT_TRANSPORT_PANEL
    }
    return {
        "panel": panel,
        "dominant_tradeoff": "die_loss_vs_precursor_release",
    }


def normalize_moisture_regime(
    moisture_regime: Optional[str],
    water_activity: Optional[float] = None,
) -> str:
    if moisture_regime is not None:
        normalized = _MOISTURE_REGIME_ALIASES.get(str(moisture_regime).strip().lower())
        if normalized is not None:
            return normalized

    aw = 0.8 if water_activity is None else float(water_activity)
    return "hme" if aw >= 0.65 else "lme"


def compute_sme_temperature_offset(
    sme_kj_per_kg: float,
    moisture_regime: Optional[str] = None,
    water_activity: Optional[float] = None,
) -> float:
    sme = max(0.0, float(sme_kj_per_kg))
    if sme <= 0.0:
        return 0.0

    regime = normalize_moisture_regime(moisture_regime, water_activity)
    slope = 0.16 if regime == "lme" else 0.10
    return min(40.0, max(5.0, sme * slope))


def pre_extrusion_damage_load(protein_type: str) -> Dict[str, float]:
    base = _PRE_EXTRUSION_DAMAGE_BASELINES.get(str(protein_type).strip().lower())
    if base is None:
        base = _PRE_EXTRUSION_DAMAGE_BASELINES["free"]
    return dict(base)


def autoclave_damage_increment(
    temperature_celsius: Optional[float],
    time_minutes: float,
    *,
    protein_type: str,
    moisture_regime: Optional[str] = None,
    water_activity: Optional[float] = None,
) -> Dict[str, float]:
    if temperature_celsius is None or time_minutes <= 0.0:
        return {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0}

    temp_c = float(temperature_celsius)
    if temp_c < 121.0:
        return {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0}

    regime = normalize_moisture_regime(moisture_regime, water_activity)
    severity = max(0.0, temp_c - 121.0) / 5.0 + math.log1p(float(time_minutes) / 10.0)
    moisture_factor = 1.15 if regime == "hme" else 0.90
    protein_factor = 1.10 if "iso" in str(protein_type).lower() else 1.0

    return {
        "furosine_mg_per_kg": 5.0 * severity * moisture_factor * protein_factor,
        "lal_mg_per_kg": 9.0 * severity * protein_factor,
    }


def build_sequential_isothermal_zones(
    base_temperature_celsius: float,
    effective_temperature_celsius: float,
    *,
    zone_temperatures: Optional[Sequence[float]] = None,
    zone_time_fractions: Optional[Sequence[float]] = None,
    zone_count: int = 3,
) -> List[Dict[str, float]]:
    if zone_temperatures:
        temperatures = [float(value) for value in zone_temperatures]
    else:
        count = max(2, int(zone_count))
        inlet = max(60.0, min(float(base_temperature_celsius), float(effective_temperature_celsius)) - 20.0)
        outlet = float(effective_temperature_celsius)
        if count == 2:
            temperatures = [inlet, outlet]
        else:
            step = (outlet - inlet) / float(count - 1)
            temperatures = [inlet + index * step for index in range(count)]

    if zone_time_fractions:
        fractions = [max(0.0, float(value)) for value in zone_time_fractions]
    else:
        if len(temperatures) == 2:
            fractions = [0.45, 0.55]
        elif len(temperatures) == 3:
            fractions = [0.20, 0.35, 0.45]
        else:
            fractions = [1.0 / float(len(temperatures))] * len(temperatures)

    total = sum(fractions) or 1.0
    normalized = [value / total for value in fractions]
    return [
        {
            "zone_index": float(index + 1),
            "temperature_celsius": float(temp),
            "time_fraction": float(fraction),
        }
        for index, (temp, fraction) in enumerate(zip(temperatures, normalized))
    ]


def build_extrusion_process_profile(
    *,
    base_temperature_celsius: float,
    water_activity: float,
    protein_type: str,
    sme_kj_per_kg: float = 0.0,
    moisture_regime: Optional[str] = None,
    screw_speed_rpm: Optional[float] = None,
    feed_rate_kg_per_h: Optional[float] = None,
    die_exit_temperature_celsius: Optional[float] = None,
    sterilization_temperature_celsius: Optional[float] = None,
    sterilization_time_minutes: float = 0.0,
    zone_temperatures: Optional[Sequence[float]] = None,
    zone_time_fractions: Optional[Sequence[float]] = None,
) -> Dict[str, object]:
    regime = normalize_moisture_regime(moisture_regime, water_activity)
    offset_c = compute_sme_temperature_offset(
        sme_kj_per_kg,
        moisture_regime=regime,
        water_activity=water_activity,
    )
    effective_temp_c = float(base_temperature_celsius) + offset_c
    base_load = pre_extrusion_damage_load(protein_type)
    sterilization_increment = autoclave_damage_increment(
        sterilization_temperature_celsius,
        sterilization_time_minutes,
        protein_type=protein_type,
        moisture_regime=regime,
        water_activity=water_activity,
    )
    total_damage_load = {
        key: float(base_load.get(key, 0.0)) + float(sterilization_increment.get(key, 0.0))
        for key in ("furosine_mg_per_kg", "lal_mg_per_kg")
    }
    zones = build_sequential_isothermal_zones(
        base_temperature_celsius,
        effective_temp_c,
        zone_temperatures=zone_temperatures,
        zone_time_fractions=zone_time_fractions,
    )
    resolved_die_temp = estimate_die_exit_temperature_celsius(
        effective_temp_c,
        regime,
        die_exit_temperature_celsius=die_exit_temperature_celsius,
    )
    mean_residence = estimate_mean_residence_time_seconds(
        sme_kj_per_kg=sme_kj_per_kg,
        moisture_regime=regime,
        screw_speed_rpm=screw_speed_rpm,
        feed_rate_kg_per_h=feed_rate_kg_per_h,
    )
    relative_rtd_spread = estimate_relative_rtd_spread(
        sme_kj_per_kg=sme_kj_per_kg,
        moisture_regime=regime,
        screw_speed_rpm=screw_speed_rpm,
    )

    profile = {
        "active": offset_c > 0.0 or bool(sterilization_temperature_celsius) or regime in {"lme", "hme"},
        "model": "sequential_isothermal_zones",
        "moisture_regime": regime,
        "water_activity": float(water_activity),
        "sme_kj_per_kg": float(max(0.0, sme_kj_per_kg)),
        "screw_speed_rpm": None if screw_speed_rpm is None else float(screw_speed_rpm),
        "feed_rate_kg_per_h": None if feed_rate_kg_per_h is None else float(feed_rate_kg_per_h),
        "jacket_temperature_celsius": float(base_temperature_celsius),
        "sme_temperature_offset_celsius": float(offset_c),
        "effective_temperature_celsius": float(effective_temp_c),
        "die_exit_temperature_celsius": float(resolved_die_temp),
        "mean_residence_time_seconds": float(mean_residence),
        "relative_rtd_spread": float(relative_rtd_spread),
        "pre_extrusion_damage": base_load,
        "sterilization": {
            "enabled": sterilization_temperature_celsius is not None and sterilization_time_minutes > 0.0,
            "temperature_celsius": None if sterilization_temperature_celsius is None else float(sterilization_temperature_celsius),
            "time_minutes": float(max(0.0, sterilization_time_minutes)),
            "damage_increment": sterilization_increment,
        },
        "total_damage_load": total_damage_load,
        "zone_profile": zones,
    }
    profile["rtd_assessment"] = evaluate_extrusion_rtd_need(profile)
    profile["volatile_transport"] = _build_transport_summary(profile)
    return profile