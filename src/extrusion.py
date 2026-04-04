from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence


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

    return {
        "active": offset_c > 0.0 or bool(sterilization_temperature_celsius) or regime in {"lme", "hme"},
        "model": "sequential_isothermal_zones",
        "moisture_regime": regime,
        "water_activity": float(water_activity),
        "sme_kj_per_kg": float(max(0.0, sme_kj_per_kg)),
        "jacket_temperature_celsius": float(base_temperature_celsius),
        "sme_temperature_offset_celsius": float(offset_c),
        "effective_temperature_celsius": float(effective_temp_c),
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