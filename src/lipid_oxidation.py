"""
src/lipid_oxidation.py - Radical chain mechanism for lipid oxidation in plant protein matrices.
Models generation of key off-flavor aldehydes from polyunsaturated fatty acids.
"""

import numpy as np
from dataclasses import dataclass
from typing import Mapping, Sequence

@dataclass
class LipidProfile:
    linoleic_acid_pct: float      # C18:2 — primary hexanal precursor
    alpha_linolenic_pct: float    # C18:3 — propanal/hexanal precursor
    oleic_acid_pct: float         # C18:1 — more oxidatively stable
    total_lipid_pct: float        # weight % in dry ingredient
    pro_oxidant_iron_ppm: float   # non-heme iron in plant material

# Typical profiles from literature
PEA_LIPID_PROFILE = LipidProfile(
    linoleic_acid_pct=50.0,
    alpha_linolenic_pct=12.0,
    oleic_acid_pct=22.0,
    total_lipid_pct=2.5,      # pea isolate ~1-3% lipid
    pro_oxidant_iron_ppm=25.0
)

SOY_LIPID_PROFILE = LipidProfile(
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
    return {
        "lipid_input_proxy_load": float(sum(float(value) for value in lipid_input.values())),
        "generated_markers": {name: float(value) for name, value in ordered_markers},
        "benchmark_ready_targets": [name for name, _value in ordered_markers if name in {"Hexanal", "2-Pentylfuran", "1-Hexanol", "Nonanal"}],
        "dominant_marker": ordered_markers[0][0] if ordered_markers else "none",
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
    
    T_K = temp_C + 273.15
    Ea_init = 80000.0  # J/mol
    R = 8.314
    
    k_init = 1e8 * np.exp(-Ea_init / (R * T_K))
    fe_factor = 1.0 + 25.0 * 0.05 # Defaulting to pea-like iron
    
    # Simple rate-based generation
    oxidation_rate = k_init * fe_factor * oxygen_availability
    
    # MFT / Hexanal SMILES for crosstalk
    # hexanal = "CCCCCC=O"
    # 2-pentylfuran = "CCCCCc1ccco1"
    
    # Scale based on time and temp
    load = oxidation_rate * time_min * total_lipid_load * 1e4
    
    return {
        "CCCCCC=O": load * 0.37,          # hexanal
        "CCCCCc1ccco1": load * 0.08,      # 2-pentylfuran
        "CCCCCCO": load * 0.05,           # 1-hexanol
        "CCCCCCCCC=O": load * 0.12        # nonanal
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
    """
    T_K = temp_C + 273.15
    Ea_init = 80000  # J/mol
    R = 8.314

    k_init = 1e8 * np.exp(-Ea_init / (R * T_K))
    fe_factor = 1.0 + lipid_profile.pro_oxidant_iron_ppm * 0.05
    ao_factor = max(0.0, 1.0 - antioxidant_mM / 5.0)

    linoleic_fraction = lipid_profile.linoleic_acid_pct / 100.0
    total_lipid = lipid_profile.total_lipid_pct / 100.0

    oxidation_rate = k_init * fe_factor * ao_factor * oxygen_availability
    hydroperoxide_ppm = oxidation_rate * linoleic_fraction * total_lipid * time_min * 1e6

    return {
        "hexanal": hydroperoxide_ppm * 0.37,
        "2-pentylfuran": hydroperoxide_ppm * 0.08,
        "nonanal": hydroperoxide_ppm * 0.15,
        "total_hydroperoxide": hydroperoxide_ppm,
    }
