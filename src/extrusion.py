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


# ---------------------------------------------------------------------------
# SME response calibration (recalibrated 2026-08-27, audit item 3.1)
#
# The previous forms saturated across the whole real twin-screw window: the
# temperature offset hit its 40 C cap at 250 (LME) / 400 (HME) kJ/kg, the
# headspace shear term hit min(1, SME/180) at 180 kJ/kg, the residence-time
# factor hit its 0.85 floor at 460 kJ/kg and the RTD-spread sigmoid was
# effectively 1.0 from 280 kJ/kg upwards.  Industrial twin-screw operation runs
# at roughly 300-800 kJ/kg, so every SME quantity in this module was constant
# over the entire window of interest.
#
# The temperature offset is now an explicit melt energy balance,
#     dT = retained_fraction * SME / cp,
# closed by a tanh soft ceiling so nothing is ever exactly saturated.
#
# Anchors and tiers:
#  * cp of a protein melt: 3.1 kJ/(kg K) at HME moisture (~57 wt%, the Li et al.
#    2026 control in data/lit/benchmark_intake_registry.json
#    `pmc_2026_hme_hexanal_baseline`), 2.2 kJ/(kg K) at LME moisture
#    (~25-30 wt%).  Water-weighted estimate from the moisture levels the repo
#    carries; ESTIMATED TIER (no calorimetric anchor in the repo).
#  * retained fraction: the share of mechanical dissipation that shows up as a
#    melt-temperature rise rather than leaving through barrel cooling and
#    evaporation.  Set so a 145-160 C barrel lands inside the 140-180 C
#    extrusion window that `de_leyn_2019_thiamine_retention`
#    (data/lit/process_state_calibrations.json) reports for structured soy, and
#    so the drier LME regime heats more per unit SME, per the moisture-as-
#    lubricant behaviour in `ilo_1996_maize_sme_lysine_damage`
#    (Ilo & Berghofer 2003, data/lit/computational_priors.json).
#    ESTIMATED TIER: directionally anchored, magnitude not fitted to data.
#  * SME_WINDOW_REFERENCE_KJ_PER_KG = 300 and the 450 kJ/kg half-saturation
#    constants below are placed at the low end and the middle of the real
#    300-800 kJ/kg window so the responses discriminate inside it.
#    ESTIMATED TIER; the only SME-resolved repo anchor
#    (`raman_sds_extrusion_disulfide_severity`) is directional only
#    ("exact_values_extracted": false).
_SME_MELT_HEAT_CAPACITY_KJ_PER_KG_K = {"hme": 3.1, "lme": 2.2}
_SME_MELT_RETAINED_FRACTION = {"hme": 0.10, "lme": 0.16}
_SME_OFFSET_SOFT_CEILING_C = 120.0
SME_WINDOW_REFERENCE_KJ_PER_KG = 300.0
SME_WINDOW_UPPER_KJ_PER_KG = 800.0
_SME_HALF_SATURATION_KJ_PER_KG = 450.0

# Below this barrel temperature there is no thermal processing to speak of, so
# the extrusion profile is reported as inactive and contributes no damage beyond
# the ingredient baseline.  ESTIMATED TIER: furosine/LAL formation is negligible
# below ~50 C on process time scales, and every extrusion anchor in the repo
# sits at 140-180 C.
AMBIENT_THERMAL_THRESHOLD_C = 50.0


def _sigmoid(x: float, center: float, k: float) -> float:
    try:
        return 1.0 / (1.0 + math.exp(-k * (float(x) - float(center))))
    except OverflowError:
        return 0.0 if x < center else 1.0


def classify_volatile_compound(name: str) -> str:
    compound = str(name).strip().lower()
    if any(token in compound for token in ["thiol", "sulfide", "mercapt", "mft", "fft"]):
        return "sulfur"
    # Furanic names first: furfural / 2-furaldehyde are furan-ring carbonyls and
    # partition like the furan class, and neither spelling contains "furan".
    if any(token in compound for token in ["furan", "furfur", "fural"]):
        return "furan"
    # "-aldehyde" spellings (benzaldehyde, phenylacetaldehyde, acetaldehyde)
    # contain neither "anal" nor "enal" and used to fall through to "unknown".
    if any(token in compound for token in ["anal", "enal", "aldehyde"]):
        return "aldehyde"
    if "pyrazine" in compound or "thiazole" in compound:
        return "pyrazine"
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
    # 280 rpm / 4.6 kg/h are the Li et al. (2026) HME control settings
    # (`pmc_2026_hme_hexanal_baseline`), so both factors are 1.0 at the anchor.
    speed_factor = 1.0 if rpm is None else max(0.70, min(1.35, 280.0 / rpm))
    feed_factor = 1.0 if feed is None else max(0.75, min(1.35, 4.6 / feed))
    # Power-law shortening with SME instead of the old linear ramp, which hit
    # its 0.85 floor at 460 kJ/kg and was therefore flat over most of the real
    # window.  Normalised to 1.0 at SME_WINDOW_REFERENCE_KJ_PER_KG; the 60 kJ/kg
    # softening constant keeps it finite at SME -> 0.  ESTIMATED TIER.
    sme = max(0.0, float(sme_kj_per_kg))
    # Upper clamp 1.05 keeps the no-shear default close to the previous
    # behaviour (42.0 -> 44.1 s at SME 0); the clamps are outside the real
    # 300-800 window, where the factor runs 1.00 -> 0.83 unclamped.
    sme_factor = max(0.60, min(1.05, ((SME_WINDOW_REFERENCE_KJ_PER_KG + 60.0) / (sme + 60.0)) ** 0.22))
    return max(12.0, min(90.0, base_seconds * speed_factor * feed_factor * sme_factor))


def estimate_relative_rtd_spread(
    *,
    sme_kj_per_kg: float,
    moisture_regime: Optional[str],
    screw_speed_rpm: Optional[float] = None,
) -> float:
    regime = normalize_moisture_regime(moisture_regime)
    base_spread = 0.24 if regime == "hme" else 0.30
    # Sigmoid re-centred on the middle of the real 300-800 kJ/kg window (was
    # centred at 180 with k=0.03, i.e. >0.97 across the entire window).
    sme_broadening = (0.03 if regime == "hme" else 0.06) * _sigmoid(
        float(sme_kj_per_kg), _SME_HALF_SATURATION_KJ_PER_KG, 0.006
    )
    rpm_compaction = 0.0
    if screw_speed_rpm is not None:
        rpm_compaction = 0.03 * _sigmoid(float(screw_speed_rpm), 260.0, 0.04)
    return max(0.15, min(0.45, base_spread + sme_broadening - rpm_compaction))


def evaluate_extrusion_rtd_need(extrusion_process: Mapping[str, object]) -> Dict[str, object]:
    regime = normalize_moisture_regime(str(extrusion_process.get("moisture_regime", "hme")))
    sme = float(extrusion_process.get("sme_kj_per_kg", 0.0) or 0.0)
    spread = float(extrusion_process.get("relative_rtd_spread", 0.0) or 0.0)
    mean_residence = float(extrusion_process.get("mean_residence_time_seconds", 0.0) or 0.0)

    # Threshold moved from 180 to the low edge of the real twin-screw window
    # (2026-08-27): 180 kJ/kg sat below every industrial LME setting, so the
    # "simple RTD recommended" branch fired for essentially all LME runs.
    if regime == "lme" and (sme >= SME_WINDOW_REFERENCE_KJ_PER_KG or spread >= 0.33):
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

    # Saturating (Michaelis-Menten style) shear-release term with a
    # half-saturation constant in the middle of the real 300-800 kJ/kg window.
    # The old min(1, SME/180) was pinned at 1.0 for every industrial setting.
    # ESTIMATED TIER.
    sme_term = max(0.0, sme) / (max(0.0, sme) + _SME_HALF_SATURATION_KJ_PER_KG)
    release_gain = float(params["shear_release_gain"]) * sme_term * (1.0 if regime == "hme" else 0.82)
    aggregation_penalty = float(params["aggregation_penalty"]) * _sigmoid(effective_temp, float(params["aggregation_center_c"]), 0.09)
    if regime == "lme":
        aggregation_penalty *= 1.15
    if compound_class == "sulfur":
        aggregation_penalty *= 1.10
    shear_release_factor = max(0.75, min(1.45, 1.0 + release_gain - aggregation_penalty))

    # Moisture modulation is now a bounded (0, 1] scaling of the class ceiling
    # instead of a >1 multiplier clipped by min(): previously the product
    # activation * moisture_factor exceeded 1 for any aw above ~0.60 once the
    # die was hotter than ~140 C, so the min() pinned the fraction at the class
    # maximum and water activity had no effect at all in exactly the regime
    # where flash-off is strongest.  Range 0.55-1.00 over aw 0.30-0.85.
    # ESTIMATED TIER (steam flash-off at the die is directional in the repo's
    # HME anchors; no quantitative aw-resolved stripping dataset is carried).
    moisture_flash_factor = 0.55 + 0.45 * max(0.0, min(1.0, (water_activity - 0.30) / 0.55))
    die_flash_activation = _sigmoid(die_temp, float(params["die_flash_center_c"]), 0.09)
    die_stripping_fraction = float(params["max_die_stripping_fraction"]) * die_flash_activation * moisture_flash_factor
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
    """Melt-temperature rise above the barrel setpoint caused by shear heating.

    Recalibrated 2026-08-27 (audit item 3.1).  Melt energy balance
    ``dT = retained_fraction * SME / cp`` with a tanh soft ceiling, replacing
    ``min(40, max(5, slope * SME))`` which was saturated from 250 kJ/kg (LME) /
    400 kJ/kg (HME) upwards — i.e. across nearly all of the real twin-screw
    window — and jumped discontinuously to 5 C for any SME above zero.

    Representative values (C): HME 300/500/800 kJ/kg -> 9.7 / 16.0 / 25.4;
    LME 300/500/800 -> 21.5 / 35.1 / 54.0.  With a 140-160 C barrel this keeps
    the effective temperature inside the 140-180 C structured-extrusion window
    reported by `de_leyn_2019_thiamine_retention`.
    """
    sme = max(0.0, float(sme_kj_per_kg))
    if sme <= 0.0:
        return 0.0

    regime = normalize_moisture_regime(moisture_regime, water_activity)
    cp = _SME_MELT_HEAT_CAPACITY_KJ_PER_KG_K[regime]
    retained = _SME_MELT_RETAINED_FRACTION[regime]
    adiabatic_rise_c = retained * sme / cp
    # Soft ceiling: strictly increasing everywhere, so the response never goes
    # flat, while keeping physically absurd offsets bounded at extreme SME.
    return _SME_OFFSET_SOFT_CEILING_C * math.tanh(adiabatic_rise_c / _SME_OFFSET_SOFT_CEILING_C)


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

    # Thermal activation ramp instead of the old `if temp_c < 121: return 0`
    # cliff, which jumped from zero to the full time term at 121.00 C (a 25 min
    # hold went 0 -> ~7 mg/kg furosine across 0.01 C).  Centred on the 121 C
    # retort reference with a ~3 C width; the sigmoid is the smoothing device,
    # the severity term below still carries the temperature dependence.
    # ESTIMATED TIER: the width is a smoothing choice, not a fitted constant.
    activation = _sigmoid(temp_c, 121.0, 1.2)
    if activation < 1e-6:
        return {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0}

    regime = normalize_moisture_regime(moisture_regime, water_activity)
    severity = activation * (max(0.0, temp_c - 121.0) / 5.0 + math.log1p(float(time_minutes) / 10.0))
    moisture_factor = 1.15 if regime == "hme" else 0.90
    protein_factor = 1.10 if "iso" in str(protein_type).lower() else 1.0

    return {
        "furosine_mg_per_kg": 5.0 * severity * moisture_factor * protein_factor,
        # Moisture factor now applied to LAL as well (2026-08-27): LAL forms by
        # base-catalysed beta-elimination followed by Michael addition of lysine,
        # a hydrolytic route that literature report 12
        # (data/lit/extrusion_damage_reference_payloads.json, `extrusion_r12_lal`)
        # associates with wet high-temperature extrusion, and
        # `ilo_1996_maize_sme_lysine_damage` reports moisture modulating lysine
        # damage.  Using the same factor as furosine is ESTIMATED TIER: the
        # direction is anchored, the magnitude is not separately fitted.
        "lal_mg_per_kg": 9.0 * severity * moisture_factor * protein_factor,
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
        # zip() below silently truncated to the shorter list, so a caller who
        # supplied 5 barrel temperatures and 3 time fractions got a 3-zone
        # profile with the last two zones dropped and no warning.
        if len(fractions) != len(temperatures):
            raise ValueError(
                "zone_time_fractions and zone_temperatures must have the same length "
                f"(got {len(fractions)} fractions for {len(temperatures)} zones)"
            )
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
    process_time_minutes: Optional[float] = None,
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
    # --- Thermal-exposure gate (2026-08-27, G2 finding) --------------------
    # `active` used to be `... or regime in {"lme", "hme"}`, and the regime is
    # ALWAYS one of those two, so the profile was unconditionally active: a
    # 25 C / 0 min formulation still attached an extrusion process state and
    # therefore still raised Furosine/LAL safety flags — from the ingredient
    # baseline alone, attributed to a process that never ran.
    # The gate is now real thermal exposure; with no thermal exposure the
    # profile adds ZERO damage beyond the ingredient baseline, and the two
    # contributions are reported separately either way.
    sterilization_enabled = sterilization_temperature_celsius is not None and sterilization_time_minutes > 0.0
    has_process_time = process_time_minutes is None or float(process_time_minutes) > 0.0
    thermal_exposure_reasons: List[str] = []
    if float(base_temperature_celsius) >= AMBIENT_THERMAL_THRESHOLD_C:
        thermal_exposure_reasons.append("barrel_temperature_above_ambient_threshold")
    if float(max(0.0, sme_kj_per_kg)) > 0.0:
        thermal_exposure_reasons.append("mechanical_energy_input")
    if sterilization_enabled:
        thermal_exposure_reasons.append("sterilization_hold")
    thermal_exposure = bool(thermal_exposure_reasons) and has_process_time
    if not has_process_time:
        thermal_exposure_reasons = []

    process_damage_increment = (
        dict(sterilization_increment)
        if thermal_exposure
        else {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0}
    )
    total_damage_load = {
        key: float(base_load.get(key, 0.0)) + float(process_damage_increment.get(key, 0.0))
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
        "active": thermal_exposure,
        "thermal_exposure": thermal_exposure,
        "thermal_exposure_basis": thermal_exposure_reasons,
        "ambient_thermal_threshold_celsius": AMBIENT_THERMAL_THRESHOLD_C,
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
        # Damage attribution is explicit: what the ingredient already carried
        # before any processing, versus what this process state added.
        "process_damage_increment": process_damage_increment,
        "damage_load_attribution": {
            "ingredient_baseline": dict(base_load),
            "process_increment": dict(process_damage_increment),
            "note": (
                "total_damage_load = ingredient_baseline + process_increment. "
                "With no thermal exposure the process increment is exactly zero, "
                "so any non-zero total is ingredient provenance, not processing."
            ),
        },
        "total_damage_load": total_damage_load,
        "zone_profile": zones,
    }
    profile["rtd_assessment"] = evaluate_extrusion_rtd_need(profile)
    profile["volatile_transport"] = _build_transport_summary(profile)
    return profile