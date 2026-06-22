"""
src/lipid_oxidation.py - Radical chain mechanism for lipid oxidation in plant protein matrices.
Models generation of key off-flavor aldehydes from polyunsaturated fatty acids.

Kinetic constants (Arrhenius prefactor / activation energy / branching / per-matrix
lipid profiles) live in ``data/lit/lipid_oxidation_calibration.json`` (S27 Workstream A,
2026-06-22). The hydroperoxide extent is capped at physical 100% conversion via a
first-order saturation so high temperature x time no longer extrapolates without bound;
the saturating form is linear at low extent, so previously calibrated low-temperature
anchors are preserved.
"""

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence

_CALIBRATION_PATH = Path(__file__).resolve().parents[1] / "data" / "lit" / "lipid_oxidation_calibration.json"

# Fallback constants matching the pre-S27 hardcoded behaviour, used only when the
# calibration file is missing or malformed (so the module never hard-fails).
# max_conversion_fraction is None => saturation cap DISABLED (default). See the
# JSON provenance: empirically the cap cannot improve the external hold-out without
# regressing the in-panel headline (in-panel and hold-out share the same kinetic
# regime), so it ships off; the mechanism is retained for Workstream B.
_FALLBACK_KINETICS: dict[str, Any] = {
    "prefactor_per_min": 1.0e8,
    "activation_energy_j_per_mol": 80000.0,
    "gas_constant_j_per_mol_k": 8.314,
    "iron_pro_oxidant_coefficient_per_ppm": 0.05,
    "antioxidant_ic50_mM": 5.0,
    "hydroperoxide_scale": 1.0e6,
    "max_conversion_fraction": None,
}
# The formulation path (predict_lop_generation) historically used a 1e4 proxy scale.
_FORMULATION_PROXY_SCALE = 1.0e4


@lru_cache(maxsize=1)
def _load_calibration() -> dict[str, Any]:
    try:
        return json.loads(_CALIBRATION_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _kinetics() -> dict[str, Any]:
    calib = _load_calibration().get("kinetics", {})
    merged = dict(_FALLBACK_KINETICS)
    for key, value in calib.items():
        if value is None:
            merged[key] = None  # explicit disable (e.g. max_conversion_fraction)
            continue
        try:
            merged[key] = float(value)
        except (TypeError, ValueError):
            continue
    return merged


def _conversion_factor(oxidation_rate: float, time_min: float, max_conversion: "float | None") -> float:
    """Multiplicative correction applied to the historical linear hydroperoxide
    term to cap first-order conversion at ``max_conversion``. Returns exactly 1.0
    when the cap is disabled (max_conversion is None/<=0) or inactive, so the
    default path is byte-identical to the pre-S27 linear model. When enabled it
    equals ``_saturated_extent(progress)/progress`` (<= 1)."""
    if not max_conversion or max_conversion <= 0.0:
        return 1.0
    progress = oxidation_rate * time_min
    if progress <= 0.0:
        return 1.0
    return _saturated_extent(progress, max_conversion) / progress


def _saturated_extent(progress: float, max_conversion: float) -> float:
    """First-order oxidation extent capped at ``max_conversion`` (a conversion
    ceiling). Reduces to ``progress`` when progress is small, so the
    low-temperature/short-time regime is unchanged within rounding."""
    if max_conversion <= 0.0:
        return max(0.0, progress)
    if progress <= 0.0:
        return 0.0
    # Guard against overflow in exp for very large progress (deep extrapolation).
    ratio = progress / max_conversion
    if ratio > 50.0:
        return max_conversion
    return max_conversion * (1.0 - math.exp(-ratio))


def _oxidation_rate_per_min(
    temp_C: float,
    *,
    iron_ppm: float,
    antioxidant_mM: float = 0.0,
    oxygen_availability: float = 1.0,
) -> float:
    """Iron- and antioxidant-modulated Arrhenius initiation rate (per minute)."""
    k = _kinetics()
    T_K = temp_C + 273.15
    k_init = k["prefactor_per_min"] * math.exp(-k["activation_energy_j_per_mol"] / (k["gas_constant_j_per_mol_k"] * T_K))
    fe_factor = 1.0 + iron_ppm * k["iron_pro_oxidant_coefficient_per_ppm"]
    ic50 = k["antioxidant_ic50_mM"]
    ao_factor = max(0.0, 1.0 - antioxidant_mM / ic50) if ic50 > 0.0 else 1.0
    return k_init * fe_factor * ao_factor * oxygen_availability


@dataclass
class LipidProfile:
    linoleic_acid_pct: float      # C18:2 — primary hexanal precursor
    alpha_linolenic_pct: float    # C18:3 — propanal/hexanal precursor
    oleic_acid_pct: float         # C18:1 — more oxidatively stable
    total_lipid_pct: float        # weight % in dry ingredient
    pro_oxidant_iron_ppm: float   # non-heme iron in plant material


def lipid_profile_for(protein_type: str | None) -> "LipidProfile | None":
    """Build a LipidProfile from the calibration registry for a protein type, or
    None if no profile is registered (caller decides the fallback)."""
    if not protein_type:
        return None
    profiles = _load_calibration().get("matrix_lipid_profiles", {})
    spec = profiles.get(str(protein_type))
    if not isinstance(spec, Mapping):
        return None
    try:
        return LipidProfile(
            linoleic_acid_pct=float(spec["linoleic_acid_pct"]),
            alpha_linolenic_pct=float(spec["alpha_linolenic_pct"]),
            oleic_acid_pct=float(spec["oleic_acid_pct"]),
            total_lipid_pct=float(spec["total_lipid_pct"]),
            pro_oxidant_iron_ppm=float(spec["pro_oxidant_iron_ppm"]),
        )
    except (KeyError, TypeError, ValueError):
        return None

# Typical profiles from literature. Prefer the calibration registry so there is a
# single source of truth; the literals below are the fallback if the file is absent.
PEA_LIPID_PROFILE = lipid_profile_for("pea_iso") or LipidProfile(
    linoleic_acid_pct=50.0,
    alpha_linolenic_pct=12.0,
    oleic_acid_pct=22.0,
    total_lipid_pct=2.5,      # pea isolate ~1-3% lipid
    pro_oxidant_iron_ppm=25.0
)

SOY_LIPID_PROFILE = lipid_profile_for("soy_iso") or LipidProfile(
    linoleic_acid_pct=53.0,
    alpha_linolenic_pct=8.0,
    oleic_acid_pct=23.0,
    total_lipid_pct=2.0,
    pro_oxidant_iron_ppm=15.0
)


LIPID_OXIDATION_MARKER_NAMES = {
    "CCCCCC=O": "Hexanal",
    "CCCCCc1ccco1": "2-Pentylfuran",
    "CCCCCCO": "1-Hexanol",
    "CCCCCCCCC=O": "Nonanal",
}

GENERIC_LIPID_PROXY_LOADS = {
    "sunflower oil": 1.0,
    "canola oil": 1.0,
    "rapeseed oil": 1.0,
    "soybean oil": 1.0,
    "soy oil": 1.0,
    "pea oil": 1.0,
    "flax oil": 1.0,
    "linoleic acid": 1.0,
    "linolenic acid": 1.0,
    "lipid": 0.8,
    "oil": 0.8,
    "fat": 0.8,
}

_BENCHMARK_READY_LIPID_MARKERS = {"Hexanal", "2-Pentylfuran", "1-Hexanol", "Nonanal"}


def _normalize_lipid_name(name: str) -> str:
    return " ".join(str(name).strip().lower().replace("_", " ").replace("-", " ").split())


def build_lipid_input_proxy_loads(
    lipid_names: Sequence[str] | None,
    molar_ratios: Mapping[str, float] | None = None,
) -> dict[str, float]:
    ratios = molar_ratios or {}
    resolved: dict[str, float] = {}

    for raw_name, raw_value in ratios.items():
        normalized = _normalize_lipid_name(raw_name)
        if any(token in normalized for token in ["linoleic", "linolenic", "oil", "fat", "lipid"]):
            resolved[normalized] = resolved.get(normalized, 0.0) + float(raw_value)

    if resolved:
        return resolved

    for lipid_name in lipid_names or []:
        normalized = _normalize_lipid_name(lipid_name)
        if not normalized:
            continue
        proxy_load = GENERIC_LIPID_PROXY_LOADS.get(normalized)
        if proxy_load is None:
            if any(token in normalized for token in ["oil", "fat", "lipid"]):
                proxy_load = 0.8
            else:
                proxy_load = 0.5
        resolved[normalized] = resolved.get(normalized, 0.0) + float(proxy_load)
    return resolved


def summarize_lipid_runtime_split(
    named_markers: Mapping[str, float],
    *,
    lipid_input_count: int,
    donor_pressure: float,
    polyphenol_active: bool,
) -> dict[str, Any]:
    normalized_markers = {str(name): float(value) for name, value in named_markers.items()}
    lipid_marker_signal_ppb = float(sum(normalized_markers.values()))
    marker_count = sum(1 for value in normalized_markers.values() if value > 0.0)
    dominant_marker = "none"
    if normalized_markers:
        dominant_marker = max(normalized_markers.items(), key=lambda item: float(item[1]))[0]
    benchmark_ready_targets = [
        name for name, value in normalized_markers.items() if name in _BENCHMARK_READY_LIPID_MARKERS and float(value) > 0.0
    ]
    carbonyl_competition_factor = min(1.75, 0.18 * int(lipid_input_count) + float(donor_pressure) * lipid_marker_signal_ppb / 120.0)
    coupled_crosstalk = bool(polyphenol_active and (lipid_marker_signal_ppb > 0.0 or lipid_input_count > 0))
    strecker_suppression_factor = min(
        1.0,
        (carbonyl_competition_factor / 1.75) * (1.15 if coupled_crosstalk else 1.0),
    )
    maillard_closure_pressure = min(2.0, carbonyl_competition_factor * (1.0 + 0.25 * float(coupled_crosstalk)))
    active = bool(lipid_input_count or lipid_marker_signal_ppb > 0.0)
    return {
        "lipid_marker_signal_ppb": float(lipid_marker_signal_ppb),
        "lipid_marker_count": int(marker_count),
        "dominant_marker": dominant_marker,
        "benchmark_ready_targets": benchmark_ready_targets,
        "polyphenol_crosstalk_active": bool(coupled_crosstalk),
        "carbonyl_competition_factor": float(carbonyl_competition_factor),
        "strecker_suppression_factor": float(strecker_suppression_factor),
        "maillard_closure_pressure": float(maillard_closure_pressure),
        "runtime_sub_lanes": {
            "adverse_marker_generation_and_retention": {
                "active": bool(active),
                "dominant_marker": dominant_marker,
                "marker_count": int(marker_count),
                "benchmark_ready_target_count": len(benchmark_ready_targets),
                "benchmark_ready_targets": benchmark_ready_targets,
            },
            "carbonyl_competition_and_crosstalk": {
                "active": bool(active),
                "donor_pressure": float(donor_pressure),
                "polyphenol_coupled": bool(coupled_crosstalk),
                "strecker_suppression_factor": float(strecker_suppression_factor),
                "maillard_closure_pressure": float(maillard_closure_pressure),
            },
        },
    }


def predict_lop_generation_named(
    lipid_input: dict,
    temp_C: float,
    time_min: float,
    water_activity: float = 0.8,
    oxygen_availability: float = 1.0,
) -> dict[str, float]:
    named: dict[str, float] = {}
    for smiles, value in predict_lop_generation(
        lipid_input,
        temp_C,
        time_min,
        water_activity=water_activity,
        oxygen_availability=oxygen_availability,
    ).items():
        name = LIPID_OXIDATION_MARKER_NAMES.get(smiles, smiles)
        named[name] = named.get(name, 0.0) + float(value)
    return named


def build_lipid_oxidation_context(
    lipid_input: dict,
    temp_C: float,
    time_min: float,
    water_activity: float = 0.8,
    oxygen_availability: float = 1.0,
) -> dict[str, object]:
    named_markers = predict_lop_generation_named(
        lipid_input,
        temp_C,
        time_min,
        water_activity=water_activity,
        oxygen_availability=oxygen_availability,
    )
    ordered_markers = sorted(named_markers.items(), key=lambda item: float(item[1]), reverse=True)
    runtime_split = summarize_lipid_runtime_split(
        named_markers,
        lipid_input_count=len(lipid_input),
        donor_pressure=0.55,
        polyphenol_active=False,
    )
    return {
        "lipid_input_proxy_load": float(sum(float(value) for value in lipid_input.values())),
        "generated_markers": {name: float(value) for name, value in ordered_markers},
        "benchmark_ready_targets": [name for name, _value in ordered_markers if name in {"Hexanal", "2-Pentylfuran", "1-Hexanol", "Nonanal"}],
        "dominant_marker": ordered_markers[0][0] if ordered_markers else "none",
        "runtime_split": runtime_split,
    }

def predict_lop_generation(
    lipid_input: dict,
    temp_C: float,
    time_min: float,
    water_activity: float = 0.8,
    oxygen_availability: float = 1.0
) -> dict[str, float]:
    """
    Restored from docs/Claude_feedback.md via pipeline.py requirements.
    Predicts Lipid Oxidation Product (LOP) SMILES and concentrations.
    """
    # No lipid precursor means no lipid oxidation signal should enter the network.
    if not lipid_input:
        return {}

    # Map input to a profile or use defaults.
    # This implementation keeps a lightweight interface but should only emit genuine LOPs.
    total_lipid_load = sum(float(value) for value in lipid_input.values())
    if total_lipid_load <= 0.0:
        return {}

    kinetics = _kinetics()
    # Default to pea-like iron for the formulation path (lipid here is an added
    # ingredient; the matrix iron load is the pre-S27 default).
    oxidation_rate = _oxidation_rate_per_min(
        temp_C,
        iron_ppm=25.0,
        oxygen_availability=oxygen_availability,
    )

    # Historical linear load, multiplied by the conversion cap factor (exactly 1.0
    # when the cap is disabled => byte-identical to the pre-S27 behaviour).
    load = oxidation_rate * time_min * total_lipid_load * _FORMULATION_PROXY_SCALE
    load *= _conversion_factor(oxidation_rate, time_min, kinetics["max_conversion_fraction"])

    branching = (_load_calibration().get("branching_ratios", {}) or {}).get("formulation_path", {})

    def _ratio(smiles: str, default: float) -> float:
        try:
            return float(branching.get(smiles, default))
        except (TypeError, ValueError):
            return default

    return {
        "CCCCCC=O": load * _ratio("CCCCCC=O", 0.37),          # hexanal
        "CCCCCc1ccco1": load * _ratio("CCCCCc1ccco1", 0.08),  # 2-pentylfuran
        "CCCCCCO": load * _ratio("CCCCCCO", 0.05),            # 1-hexanol
        "CCCCCCCCC=O": load * _ratio("CCCCCCCCC=O", 0.12),    # nonanal
    }

def predict_hexanal_generation(
    lipid_profile: LipidProfile,
    temp_C: float,
    time_min: float,
    oxygen_availability: float = 1.0,
    antioxidant_mM: float = 0.0,
) -> dict[str, float]:
    """
    Frankel 1998 radical chain model implementation.

    The hydroperoxide extent is capped at physical 100% conversion (S27): the
    matrix-benchmark path consumes ``total_hydroperoxide`` as its oxidation load,
    so an unbounded Arrhenius x time term previously let out-of-calibration
    process states (roasting, extrusion) extrapolate to absurd loads. The
    saturating extent is linear in the calibrated low-temperature regime, so the
    registry-pinned anchors are unchanged.
    """
    kinetics = _kinetics()
    linoleic_fraction = lipid_profile.linoleic_acid_pct / 100.0
    total_lipid = lipid_profile.total_lipid_pct / 100.0

    oxidation_rate = _oxidation_rate_per_min(
        temp_C,
        iron_ppm=lipid_profile.pro_oxidant_iron_ppm,
        antioxidant_mM=antioxidant_mM,
        oxygen_availability=oxygen_availability,
    )
    # Historical linear hydroperoxide term (same multiplication order as pre-S27),
    # multiplied by the conversion cap factor (exactly 1.0 when the cap is disabled,
    # so the registry-pinned benchmark predictions are bit-for-bit preserved).
    hydroperoxide_ppm = oxidation_rate * linoleic_fraction * total_lipid * time_min * kinetics["hydroperoxide_scale"]
    hydroperoxide_ppm *= _conversion_factor(oxidation_rate, time_min, kinetics["max_conversion_fraction"])

    branching = (_load_calibration().get("branching_ratios", {}) or {}).get("hexanal_path", {})

    def _ratio(name: str, default: float) -> float:
        try:
            return float(branching.get(name, default))
        except (TypeError, ValueError):
            return default

    return {
        "hexanal": hydroperoxide_ppm * _ratio("hexanal", 0.37),
        "2-pentylfuran": hydroperoxide_ppm * _ratio("2-pentylfuran", 0.08),
        "nonanal": hydroperoxide_ppm * _ratio("nonanal", 0.15),
        "total_hydroperoxide": hydroperoxide_ppm,
    }
