#!/usr/bin/env python3
"""
src/recommend.py — Maillard Reaction Pathway Recommender

Tier 4 Pipeline Integration:
Loads the Tier 1 xTB screening results and matches them against user-defined
or canonical precursors to recommend actionable formulation adjustments.
"""

import json
import yaml
import math
import numpy as np
import pandas as pd
from pathlib import Path

from typing import List, Dict, Set, Optional, Any, Tuple
from typing import Mapping

from src import data_paths
from src import data_access

from src.logger import get_logger
logger = get_logger(__name__)

from src.projection import (
    ProjectionBudget, ProjectionStrategy, DEFAULT_PROJECTION_STRATEGY,
    _thermal_severity, _projection_temperature_factor, _projection_time_factor,
    _estimate_projection_budget, _temporal_accessibility, _relative_precursor_load_factor,
    _projection_strategy_metadata
)
from src.matrix_targets import get_compound_panel_entry
from src.projection_metadata import ProjectionMetadataMap, make_projection_metadata_row

from src.barrier_constants import arrhenius_rate_constant, get_reference_pre_exponential
from src.headspace import HeadspaceModel
from src.matrix_calibration_registry import describe_matrix_calibration, determine_matrix_process_state
from src.matrix_correction import (
    ProteinType,
    apply_matrix_correction,
    classify_accessibility_state,
    classify_volatile_matrix_family,
    describe_compound_matrix_retention,
    get_protein_source_profile,
    get_volatile_class_retention_factor,
    resolve_compound_matrix_retention,
    resolve_matrix_correction,
)
from src.chem_utils import Species
try:
    from rdkit import Chem
except ImportError:
    Chem = None


if Chem is not None:
    _CARBOXYLIC_ACID_SMARTS = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    _PRIMARY_AMINE_SMARTS = Chem.MolFromSmarts("[NX3;H1,H2;!$(NC=O)]")
    _IMINE_SMARTS = Chem.MolFromSmarts("[CX3]=[NX2;!R]")
else:
    _CARBOXYLIC_ACID_SMARTS = None
    _PRIMARY_AMINE_SMARTS = None
    _IMINE_SMARTS = None

_HENRY_CONSTANTS_PATH = data_paths.HENRY_CONSTANTS
_NON_OBSERVABLE_KAW_THRESHOLD = 1.0e-8
_LOW_HEADSPACE_REFERENCE_KAW = 1.0e-5


def _normalize_chemical_name(name: str) -> str:
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())


def _load_henry_lookup() -> Dict[str, Dict[str, Any]]:
    raw = data_access.load_yaml(_HENRY_CONSTANTS_PATH) or {}
    constants = raw.get("constants", [])
    lookup: Dict[str, Dict[str, Any]] = {}
    for entry in constants:
        if not entry.get("name"):
            continue
        lookup[_normalize_chemical_name(entry["name"])] = entry
    return lookup


_HENRY_LOOKUP = _load_henry_lookup()
_HEADSPACE_MODEL = HeadspaceModel(str(_HENRY_CONSTANTS_PATH))

def _trunc(s: str, max_len: int) -> str:
    """Pad or truncate string for fixed-width columns."""
    if s is None:
        s = "-"
    # Handle invisible unicode characters (like emojis) which mess up standard len()
    # Simple approximation: standard str.ljust, but we'll try to keep it simple.
    s = str(s)
    # Emojis count as 1 char but render wider in some terminals, 
    # but ljust treats them as 1 char. 
    if len(s) > max_len:
        return s[:max_len-3] + "..."
    return s.ljust(max_len)

# Imports moved to top





# Build canonical SMILES lookup for targets
from src.chem_utils import canonicalize_smiles
def _canon(smi: str) -> str:
    return canonicalize_smiles(smi, fallback_to_original=True, strip_salts=True)

def _weight(barrier_kcal, temp_kelvin=423.15): # Default 150C
    import math
    if barrier_kcal >= 99.0: 
        return 0.0
    R = 0.001987
    return math.exp(-barrier_kcal / (R * temp_kelvin))

def _integrate_arrhenius(barrier_kcal: float, 
                         ramp_data: List[Tuple[float, float]], 
                         pre_exponential: Optional[float] = None) -> float:
    """
    SOTA: Numerically integrate the Arrhenius propensity over a temperature ramp.
    Returns the integrated yield approximation: 1 - exp(-integral(k(t) dt)).
    
    ramp_data: List of (time_minutes, temp_kelvin)
    """
    if not ramp_data or barrier_kcal >= 99.0:
        return 0.0
    
    # Convert minutes to seconds for S-1 pre-exponential
    times_min, temps_k = zip(*ramp_data)
    times_sec = np.array(times_min) * 60.0
    temps_k = np.array(temps_k)
    
    # Interpolate to 100 points for fidelity
    t_fine = np.linspace(times_sec[0], times_sec[-1], 100)
    T_fine = np.interp(t_fine, times_sec, temps_k)
    
    reference_pre_exponential = pre_exponential
    if reference_pre_exponential is None:
        reference_pre_exponential = get_reference_pre_exponential()

    k_fine = np.array(
        [arrhenius_rate_constant(barrier_kcal, temp_k, multiplier=1.0) for temp_k in T_fine]
    )
    if reference_pre_exponential != get_reference_pre_exponential():
        k_fine *= float(reference_pre_exponential) / get_reference_pre_exponential()
    
    # Integrate using Trapezoidal rule (standard for discrete steps)
    integral_k_dt = np.trapezoid(k_fine, t_fine)
    
    # Yield approximation (saturation handling)
    return 1.0 - np.exp(-integral_k_dt)

def _load_ramp(ramp_path: str) -> List[Tuple[float, float]]:
    """Load a temperature ramp from CSV (columns ``time`` in minutes, ``temp`` in Celsius).

    Raises instead of degrading: until 2026-09-01 a missing or malformed ramp file was
    logged as a warning and the run silently continued isothermally.
    """
    path = Path(ramp_path)
    if not path.exists():
        raise FileNotFoundError(f"temperature ramp file not found: {ramp_path}")
    df = pd.read_csv(path)
    if 'time' not in df.columns or 'temp' not in df.columns:
        raise ValueError(
            f"temperature ramp {ramp_path} must have 'time' and 'temp' columns; found {list(df.columns)}"
        )
    return list(zip(df['time'], df['temp'] + 273.15))


def _mw_from_smiles(smiles: str) -> float:
    if Chem is None:
        return 100.0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 100.0
    try:
        from rdkit.Chem import Descriptors
        return float(Descriptors.MolWt(mol))
    except Exception:
        return 100.0


_BUDGET_EXCLUDED_CANONICAL = {
    "O",
    "O=C=O",
    "[HH]",
    "[S]",
    "S",
    "N",
    "C=O",
}


def _carbon_count(smiles: str) -> int:
    if Chem is None:
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)


def _has_reactive_nonvolatile_functionality(smiles: str) -> bool:
    if Chem is None:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    if any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms()):
        return True
    if _CARBOXYLIC_ACID_SMARTS is not None and mol.HasSubstructMatch(_CARBOXYLIC_ACID_SMARTS):
        return True
    if _PRIMARY_AMINE_SMARTS is not None and mol.HasSubstructMatch(_PRIMARY_AMINE_SMARTS):
        return True
    if _IMINE_SMARTS is not None and mol.HasSubstructMatch(_IMINE_SMARTS):
        return True
    return False


def _henry_entry_for_species(species: Species, target_lookup: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    candidate_names: List[str] = []
    canon = _canon(species.smiles)
    target_info = target_lookup.get(canon)
    if target_info:
        candidate_names.append(target_info["name"])
    if species.label:
        candidate_names.append(species.label)

    alias_map = {
        "furfural": "Furfural",
        "hmf": "5-Hydroxymethylfurfural (HMF)",
        "2-methyl-3-furanthiol": "2-Methyl-3-furanthiol (MFT)",
        "2-furfurylthiol": "2-Furfurylthiol (FFT)",
        "2,5-dimethylpyrazine": "2,5-Dimethylpyrazine",
        "2,3-dimethylpyrazine": "2,3-Dimethylpyrazine",
        "hydrogen sulfide": "Hydrogen Sulfide",
        "dimethyl disulfide": "Dimethyl disulfide",
        "dimethyl trisulfide": "Dimethyl trisulfide",
        "acrylamide": "Acrylamide",
    }
    for name in list(candidate_names):
        normalized = _normalize_chemical_name(name)
        if normalized in alias_map:
            candidate_names.append(alias_map[normalized])

    for name in candidate_names:
        entry = _HENRY_LOOKUP.get(_normalize_chemical_name(name))
        if entry is not None:
            return entry
    return None


def _is_observable_target_species(species: Species, target_lookup: Dict[str, Dict[str, Any]]) -> bool:
    entry = _henry_entry_for_species(species, target_lookup)
    if entry is None:
        return True
    return float(entry.get("Kaw_25c", 0.01)) >= _NON_OBSERVABLE_KAW_THRESHOLD


def _headspace_observability_metadata(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    entry = _henry_entry_for_species(species, target_lookup)
    if entry is None:
        return {
            "headspace_observable": True,
            "headspace_class": "assumed_observable",
            "henry_kaw_25c": None,
            "henry_source_name": None,
        }

    kaw_25c = float(entry.get("Kaw_25c", 0.01))
    observable = kaw_25c >= _NON_OBSERVABLE_KAW_THRESHOLD
    return {
        "headspace_observable": observable,
        "headspace_class": "observable" if observable else "low_headspace",
        "henry_kaw_25c": kaw_25c,
        "henry_source_name": entry.get("name"),
    }


def _headspace_observability_factor(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
    temperature_kelvin: float,
    fat_fraction: float = 0.0,
    protein_fraction: float = 0.0,
) -> float:
    entry = _henry_entry_for_species(species, target_lookup)
    if entry is None:
        return 1.0

    kaw_25c = float(entry.get("Kaw_25c", 0.01))
    if kaw_25c >= _NON_OBSERVABLE_KAW_THRESHOLD:
        return 1.0

    compound_name = str(entry.get("name", species.label or ""))
    try:
        kaw_at_temp = float(_HEADSPACE_MODEL.get_kaw_at_temp(compound_name, temperature_kelvin))
    except Exception:
        kaw_at_temp = kaw_25c

    try:
        reference_kaw = float(_HEADSPACE_MODEL.get_kaw_at_temp("Furfural", temperature_kelvin))
    except Exception:
        reference_kaw = _LOW_HEADSPACE_REFERENCE_KAW

    reference_kaw = max(reference_kaw, _LOW_HEADSPACE_REFERENCE_KAW)
    intrinsic_factor = max(1.0e-6, min(1.0, kaw_at_temp / reference_kaw))

    matrix_factor = 1.0
    temp_c = temperature_kelvin - 273.15
    if fat_fraction > 0.0 or protein_fraction > 0.0:
        try:
            baseline_air = float(_HEADSPACE_MODEL.predict_headspace({compound_name: 1.0}, temp_c, fat_fraction=0.0, protein_fraction=0.0).get(compound_name, 1.0))
            matrix_air = float(_HEADSPACE_MODEL.predict_headspace({compound_name: 1.0}, temp_c, fat_fraction=fat_fraction, protein_fraction=protein_fraction).get(compound_name, baseline_air))
            if baseline_air > 0.0:
                matrix_factor = max(1.0e-3, min(1.0, matrix_air / baseline_air))
        except Exception:
            matrix_factor = 1.0

    # Keep the penalty conservative so we do not destabilize the validated
    # free-amino-acid benchmarks while still reflecting that near-nonvolatile
    # species should not consume the same observable headspace budget.
    return intrinsic_factor * matrix_factor


def _resolve_output_matrix_context(
    protein_type: ProteinType,
    fat_fraction: float,
    protein_fraction: float,
) -> Tuple[float, float, bool]:
    fat = max(float(fat_fraction), 0.0)
    protein = max(float(protein_fraction), 0.0)

    if protein_type == ProteinType.FREE_AMINO_ACID:
        return fat, 0.0, fat > 0.0

    # ReactionConditions defaults protein_fraction to 1.0, which in this codebase
    # often means "unspecified" rather than a real volumetric matrix fraction.
    fractions_explicit = fat > 0.0 or (0.0 < protein < 0.999)
    if not fractions_explicit:
        return 0.0, 0.0, False
    return fat, protein, True


_MELANOIDIN_TRAPPING_PROFILES = {
    _normalize_chemical_name("2-methyl-3-furanthiol"): {"slope": 0.85, "floor": 0.20},
    _normalize_chemical_name("2-furfurylthiol"): {"slope": 1.10, "floor": 0.08},
    _normalize_chemical_name("bis(2-methyl-3-furyl) disulfide"): {"slope": 0.55, "floor": 0.35},
}

# Fraction of the modelled volatile taken to survive to the headspace in a
# hydrolysed-vegetable-protein matrix. These constrain the PRODUCT of the lane, not
# any barrier.
#
# WAVE I REVERT, 2026-08-27 — READ THIS BEFORE TOUCHING THESE NUMBERS.
# Wave H (2026-08-27) re-derived the Methional factor 0.0045 -> 0.05623 against the two
# xylose HVP benchmarks `spi_hvp_xylose_120C_PMC9905368` and
# `wheat_gluten_hvp_xylose_120C_PMC9905368`. The cold-start red team then established
# that BOTH of those benchmarks are fabricated: the cited paper
# 10.1007/s10068-022-01194-w uses glucose/fructose at pH 7.5 for 90 min, reports only
# RELATIVE PEAK AREAS, and never mentions FFT or MFT — so it cannot be the source of any
# absolute ppb value in those files. The Wave H fit therefore had no measurement behind
# it, and the two "methional rows inside the CI" it produced were a self-fit scored
# against its own fit target. Both files are now quarantined
# (data/benchmarks/quarantined/) and this constant is REVERTED to its pre-Wave-H value.
# `results/validation/hydrolysate_observability_rederivation.{json,md}` is marked
# RETRACTED and must not be cited as a warrant for any value here.
#   * Methional  0.05623 -> 0.0045 (revert; the 0.05623 was fitted to fabricated data).
#   * 2-Furfurylthiol / 2-Methyl-3-furanthiol / bis(2-methyl-3-furyl) disulfide: Wave H
#     did NOT move these (their unconstrained optima were 8.65 and 3.49, above the
#     physical maximum of 1.0 for a surviving fraction; the disulfide had no comparator
#     at all), so there is nothing to revert. They remain UNCONSTRAINED LEGACY
#     ESTIMATES, never fitted values.
# NET STATUS after this revert: with the xylose HVP lanes quarantined, EVERY entry in
# this table is an unconstrained legacy estimate. There is no surviving literature
# constraint on any of them.
_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES = {
    _normalize_chemical_name("Methional"): {"base_factor": 0.0045, "source_sensitive": False},
    _normalize_chemical_name("2-Furfurylthiol"): {"base_factor": 0.13, "source_sensitive": True},
    _normalize_chemical_name("2-Methyl-3-furanthiol"): {"base_factor": 0.13, "source_sensitive": True},
    _normalize_chemical_name("bis(2-methyl-3-furyl) disulfide"): {"base_factor": 0.18, "source_sensitive": True},
}


def _resolve_melanoidin_trapping_factor(
    compound_name: str,
    *,
    protein_type: ProteinType,
    process_state: str,
    projection_severity: float,
    family_upstream_contract: Optional[Mapping[str, Any]] = None,
) -> float:
    profile = _MELANOIDIN_TRAPPING_PROFILES.get(_normalize_chemical_name(compound_name))
    if profile is None:
        return 1.0

    family_lanes = (family_upstream_contract or {}).get("family_lanes", {}) or {}
    thiamine_lane = family_lanes.get("03", {}) or {}
    sulfur_lane = family_lanes.get("05", {}) or {}
    melanoidin_lane = family_lanes.get("16", {}) or {}
    free_supported_context = bool(thiamine_lane.get("active", False) or sulfur_lane.get("active", False))
    if protein_type == ProteinType.FREE_AMINO_ACID and not (bool(melanoidin_lane.get("active", False)) and free_supported_context):
        return 1.0

    severity = max(0.0, min(1.0, float(projection_severity)))
    factor = 1.0 - float(profile["slope"]) * severity
    if process_state == "heated_matrix":
        factor *= 0.92
    return max(float(profile["floor"]), min(1.0, factor))


def _resolve_upstream_observability_factor(
    compound_name: str,
    *,
    protein_source: Optional[str],
    family_upstream_contract: Optional[Mapping[str, Any]],
) -> float:
    if not family_upstream_contract:
        return 1.0

    family_lanes = family_upstream_contract.get("family_lanes", {}) or {}
    thiamine_lane = family_lanes.get("03", {}) or {}
    sulfur_lane = family_lanes.get("05", {}) or {}
    pretreatment_lane = family_lanes.get("10", {}) or {}
    thiamine_active = bool(thiamine_lane.get("active", False))
    sulfur_active = bool(sulfur_lane.get("active", False)) or bool(pretreatment_lane.get("precursor_release_active", False))
    if not (sulfur_active or thiamine_active):
        return 1.0

    normalized_compound = _normalize_chemical_name(compound_name)
    if thiamine_active and normalized_compound in {
        _normalize_chemical_name("2-Methyl-3-furanthiol"),
        _normalize_chemical_name("2-Furfurylthiol"),
        _normalize_chemical_name("bis(2-methyl-3-furyl) disulfide"),
    }:
        thiamine_calibration = family_upstream_contract.get("thiamine_calibration", {}) or {}
        thiamine_factor = float(thiamine_calibration.get("observable_efficiency_factor", 1.0) or 1.0)
    else:
        thiamine_factor = 1.0

    if not sulfur_active:
        return max(1.0e-4, min(1.0, thiamine_factor))

    profile = _HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES.get(_normalize_chemical_name(compound_name))
    if profile is None:
        return max(1.0e-4, min(1.0, thiamine_factor))

    factor = float(profile.get("base_factor", 1.0))
    sulfur_lane = family_lanes.get("05", {}) or {}
    if profile.get("source_sensitive") and str(sulfur_lane.get("peptide_mode", "")) == "hydrolysate_supported":
        selected_peptide_ratio = float(sulfur_lane.get("selected_peptide_ratio", 1.0) or 1.0)
        peptide_accessibility_factor = float(sulfur_lane.get("peptide_accessibility_factor", 1.0) or 1.0)
        hydrolysate_release_uplift = 1.0 + max(0.0, selected_peptide_ratio - 1.0) * peptide_accessibility_factor
        factor *= max(1.0, min(selected_peptide_ratio, hydrolysate_release_uplift))
    if profile.get("source_sensitive"):
        source_profile = get_protein_source_profile(protein_source)
        if source_profile is not None:
            factor *= float(source_profile.hydrolysate_observability_bias)
    return max(1.0e-4, min(1.0, factor * thiamine_factor))


def _registry_compound_name(
    canon: str,
    species_label: Optional[str],
    target_lookup: Dict[str, Dict[str, Any]],
) -> str:
    """Resolve the NAME under which a tracked species is looked up in the name-keyed registries.

    WAVE S1 FIX (2026-08-27) — closes the Wave O finding "the matrix observability
    registry is unreachable from the recommend path".

    CAUSE. Species reach `predict_from_steps` from two places. Network species carry a
    human label ("2-methyl-3-furanthiol"). Species INJECTED into `corrected_initial`
    — the lipid-oxidation markers on the `matrix_precursor_augmented` lane — are
    materialised at the top of the relaxation as `Species(species_name_lookup.get(canon,
    canon), canon)`, and when no step ever named them the label falls back to the
    CANONICAL SMILES. `_apply_output_projection` then called
    `describe_matrix_calibration("CCCCCC=O", ...)`, which normalises to the literal
    string "ccccccc=o", matches no record and no class anchor, and returns
    `calibration_observable_factor=None` -> the caller's `or 1.0` silently applied a
    factor of ONE. Measured by Wave O on
    `soy_isolate_ribose_cysteine_100C_45min_Internal2026`: the snapshot was BIT-IDENTICAL
    after a 4.32x change to the soy heated-matrix hexanal factor, i.e. the internal
    snapshots could not detect drift in `src/matrix_calibration_registry.py` at all.

    FIX. Names, not SMILES, are the registry's key, so resolve to a name at the lookup
    boundary: keep the species' own label when it IS a name, and fall back to the
    compound-database name for that canonical SMILES when the label is merely the SMILES
    string again. `target_lookup` is keyed by canonical SMILES and its `name` is the
    entry the desirable/off-flavour/toxic YAML databases use, which is the same spelling
    `MATRIX_BENCHMARK_BASE_MARKER_YIELDS` (and therefore the registry) uses on the
    `matrix_only` lane. The two lanes now agree on the key.

    LIMITATION, stated. A species that is neither labelled nor present in the compound
    databases still reaches the registry as a SMILES and still misses. That is not
    reachable today for any scored row (only target-database species are projected to
    ppb), but the fallback is left as the SMILES rather than raising, so such a species
    degrades to the pre-fix behaviour rather than crashing.
    """
    label = str(species_label or "").strip()
    if label:
        try:
            label_is_smiles = _canon(label) == canon
        except Exception:  # pragma: no cover - canonicaliser is total, belt and braces
            label_is_smiles = False
        if not label_is_smiles:
            return label
    entry = target_lookup.get(canon)
    if entry:
        registered = str(entry.get("name") or "").strip()
        if registered:
            return registered
    return label or canon


def _apply_output_projection(
    raw_concentrations: Dict[str, float],
    species_catalog: Dict[str, Species],
    target_lookup: Dict[str, Dict[str, Any]],
    temperature_kelvin: float,
    protein_type: str,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    denaturation_state: float = 0.5,
    fat_fraction: float = 0.0,
    protein_fraction: float = 0.0,
    projection_budget: Optional[ProjectionBudget] = None,
    projection_strategy: ProjectionStrategy = DEFAULT_PROJECTION_STRATEGY,
    species_name_lookup: Optional[Dict[str, str]] = None,
    protein_source: Optional[str] = None,
    family_upstream_contract: Optional[Mapping[str, Any]] = None,
) -> Tuple[Dict[str, float], ProjectionMetadataMap]:
    p_type = ProteinType(protein_type)
    fat_eff, protein_eff, explicit_matrix_fractions = _resolve_output_matrix_context(
        p_type,
        fat_fraction,
        protein_fraction,
    )

    if explicit_matrix_fractions:
        fallback_matrix_factor = 1.0
    else:
        fallback_matrix_factor = float(resolve_matrix_correction(p_type, denaturation_state).volatile_retention)

    observable_ppb: Dict[str, float] = {}
    projection_metadata: ProjectionMetadataMap = {}
    process_state = determine_matrix_process_state(
        temperature_celsius=temperature_kelvin - 273.15,
        time_minutes=float(time_minutes or 60.0),
        water_activity=water_activity,
    )
    accessibility_state = classify_accessibility_state(
        protein_type,
        denaturation_state,
        dominant_source="denaturation_state_arg",
    )

    budget_metadata = {}
    if projection_budget is not None:
        budget_metadata = {
            "limiting_precursor_molar": float(projection_budget.limiting_precursor_molar),
            "projection_load_factor": float(projection_budget.load_factor),
            "projection_temperature_factor": float(projection_budget.temperature_factor),
            "projection_time_factor": float(projection_budget.time_factor),
            "projection_severity": float(projection_budget.severity),
            "projection_kinetic_drive": float(projection_budget.kinetic_drive),
            "projection_conversion_extent": float(projection_budget.conversion_extent),
            "volatile_yield_fraction": float(projection_budget.volatile_yield_fraction),
            "total_volatile_budget_molar": float(projection_budget.total_volatile_budget_molar),
        }
    budget_metadata.update(_projection_strategy_metadata(projection_strategy))
    projection_severity = float(
        budget_metadata.get("projection_severity", _thermal_severity(temperature_kelvin, time_minutes))
    )
    budget_metadata.setdefault("projection_severity", projection_severity)

    for canon, raw_value in raw_concentrations.items():
        species = species_catalog.get(canon)
        if species is None:
            observable_ppb[canon] = raw_value * fallback_matrix_factor
            # WAVE S1 (2026-08-27): resolve SMILES-labelled species to their registry name.
            compound_name = _registry_compound_name(
                canon,
                species_name_lookup.get(canon) if species_name_lookup else None,
                target_lookup,
            )
            calibration = describe_matrix_calibration(
                compound_name,
                protein_type=protein_type,
                process_state=process_state,
            )
            panel_entry = get_compound_panel_entry(compound_name) or {}
            calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
            melanoidin_factor = _resolve_melanoidin_trapping_factor(
                compound_name,
                protein_type=p_type,
                process_state=process_state,
                projection_severity=projection_severity,
                family_upstream_contract=family_upstream_contract,
            )
            observable_ppb[canon] *= calibration_factor * melanoidin_factor
            projection_metadata[canon] = make_projection_metadata_row(
                compound=compound_name,
                proxy_ppb=float(raw_value),
                observable_ppb=float(observable_ppb[canon]),
                extras={
                    "matrix_factor": float(fallback_matrix_factor),
                    "base_matrix_factor": float(fallback_matrix_factor),
                    "class_matrix_factor": 1.0,
                    "dynamic_retention_factor": 1.0,
                    "reversible_release_factor": 1.0,
                    "temporal_attenuation_factor": 1.0,
                    "extrusion_moisture_factor": 1.0,
                    "extrusion_structure_factor": 1.0,
                    "headspace_factor": 1.0,
                    "calibration_factor": calibration_factor,
                    "melanoidin_trapping_factor": float(melanoidin_factor),
                    "volatile_class": "other",
                    "process_state": process_state,
                    "accessibility_profile": accessibility_state.profile,
                    "accessibility_warning": accessibility_state.accessibility_warning,
                    "accessibility_dominant_source": accessibility_state.dominant_source,
                    **panel_entry,
                    **calibration,
                    **budget_metadata,
                },
            )
            continue

        if explicit_matrix_fractions:
            effective_matrix_factor = 1.0
            class_matrix_factor = 1.0
        else:
            retention_description = describe_compound_matrix_retention(
                species.label or canon,
                protein_type=p_type,
                denaturation_state=denaturation_state,
                smiles=species.smiles,
                temperature_celsius=temperature_kelvin - 273.15,
                time_minutes=time_minutes,
                water_activity=water_activity,
                process_state=process_state,
                protein_source=protein_source,
            )
            class_matrix_factor = float(retention_description.get("class_matrix_factor", 1.0))
            effective_matrix_factor = float(retention_description.get("matrix_factor", 1.0))

        headspace_factor = _headspace_observability_factor(
            species,
            target_lookup,
            temperature_kelvin,
            fat_fraction=fat_eff,
            protein_fraction=protein_eff,
        )
        # WAVE S1 (2026-08-27): resolve SMILES-labelled species to their registry name.
        # See `_registry_compound_name` — this is the boundary at which the Wave O
        # "registry unreachable on matrix_precursor_augmented" defect was located.
        compound_name = _registry_compound_name(canon, species.label, target_lookup)
        calibration = describe_matrix_calibration(
            compound_name,
            protein_type=protein_type,
            process_state=process_state,
        )
        panel_entry = get_compound_panel_entry(compound_name) or {}
        calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
        melanoidin_factor = _resolve_melanoidin_trapping_factor(
            compound_name,
            protein_type=p_type,
            process_state=process_state,
            projection_severity=projection_severity,
            family_upstream_contract=family_upstream_contract,
        )
        upstream_observability_factor = _resolve_upstream_observability_factor(
            compound_name,
            protein_source=protein_source,
            family_upstream_contract=family_upstream_contract,
        )
        observable_value = (
            raw_value
            * effective_matrix_factor
            * headspace_factor
            * calibration_factor
            * melanoidin_factor
            * upstream_observability_factor
        )
        observable_ppb[canon] = observable_value
        projection_metadata[canon] = make_projection_metadata_row(
            compound=compound_name,
            proxy_ppb=float(raw_value),
            observable_ppb=float(observable_value),
            extras={
                "matrix_factor": float(effective_matrix_factor),
                "base_matrix_factor": float(retention_description.get("base_matrix_factor", fallback_matrix_factor)),
                "class_matrix_factor": float(class_matrix_factor),
                "dynamic_retention_factor": float(retention_description.get("dynamic_retention_factor", 1.0)),
                "retention_runtime_mode": retention_description.get("retention_runtime_mode", "static_class_profile"),
                "retention_reference_sources": retention_description.get("retention_reference_sources", []),
                "reversible_release_factor": float(retention_description.get("reversible_release_factor", 1.0)),
                "temporal_attenuation_factor": float(retention_description.get("temporal_attenuation_factor", 1.0)),
                "extrusion_moisture_factor": float(retention_description.get("extrusion_moisture_factor", 1.0)),
                "extrusion_structure_factor": float(retention_description.get("extrusion_structure_factor", 1.0)),
                "headspace_factor": float(headspace_factor),
                "calibration_factor": float(calibration_factor),
                "melanoidin_trapping_factor": float(melanoidin_factor),
                "upstream_observability_factor": float(upstream_observability_factor),
                "browning_index": float(projection_severity),
                "browning_narrative": "melanoidin-linked sulfur trapping surrogate" if melanoidin_factor < 1.0 else "no explicit browning-linked sulfur penalty",
                "volatile_class": classify_volatile_matrix_family(compound_name, smiles=species.smiles),
                "process_state": process_state,
                "accessibility_profile": accessibility_state.profile,
                "accessibility_warning": accessibility_state.accessibility_warning,
                "accessibility_dominant_source": accessibility_state.dominant_source,
                **panel_entry,
                **calibration,
                **budget_metadata,
            },
        )

    return observable_ppb, projection_metadata


def _is_budget_relevant_species(species: Species, target_lookup: Dict[str, Dict[str, Any]]) -> bool:
    canon = _canon(species.smiles)
    if not canon:
        return False
    if canon in _BUDGET_EXCLUDED_CANONICAL:
        return False
    if canon in target_lookup:
        return True
    if not species.is_volatile:
        return False
    if _carbon_count(canon) < 2:
        return False
    if _has_reactive_nonvolatile_functionality(canon):
        return False
    return True


def _is_ppb_output_species(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
    exogenous_reactants: Set[str],
) -> bool:
    canon = _canon(species.smiles)
    if not canon or canon in exogenous_reactants:
        return False
    if canon not in target_lookup:
        return False
    return _is_budget_relevant_species(species, target_lookup)


def _select_accumulating_projection_species(
    tracked_species: Dict[str, Tuple[float, float, int, float, float]],
    species_catalog: Dict[str, Species],
    target_lookup: Dict[str, Dict[str, Any]],
    exogenous_reactants: Set[str],
) -> Set[str]:
    """Which tracked species the budget projection lets accumulate.

    2026-08-27 (Wave T4): two parameters REMOVED, AST-verified unused in the body
    — `steps` and `downstream_margin_kcal: float = 0.25`. Both are leftovers of a
    pre-Wave-S1 selection heuristic that compared each candidate against its
    downstream competitors by a barrier margin; Wave S1's additive propagator
    replaced that comparison entirely and the parameters were never unwired. The
    magic default is the part worth removing: a reader tuning `0.25` would have
    been tuning nothing, which is the same false-affordance class as the Wave I
    false zeros. Behaviour is unchanged by construction — the body never read
    either name.
    """
    candidate_canons: Set[str] = set()
    for canon, (span, _conc, depth, _weight, _unc) in tracked_species.items():
        species = species_catalog.get(canon)
        if species is None:
            continue
        if depth <= 0 or span == float("inf"):
            continue
        if _is_ppb_output_species(species, target_lookup, exogenous_reactants):
            candidate_canons.add(canon)

    if not candidate_canons:
        return set()
    return candidate_canons


def _route_channel_id(path: List[Dict[str, Any]]) -> str:
    """Channel identity of a route: its FULL ORDERED STEP-SET.

    2026-08-27 (Wave S1). Two enumerated routes to the same product are the same kinetic
    channel iff they are the same sequence of elementary steps. Anything else is a
    genuinely distinct branch and its conductance adds.

    THE ALTERNATIVE RULE WAS IMPLEMENTED, MEASURED AND REJECTED, and the reason is the
    load-bearing argument of this change, so it is written here rather than in a ledger
    nobody reads next to the code.

      Rejected rule: "two routes that share their RATE-LIMITING STEP are one channel;
      take the max, not the sum."  Measured on cys_ribose_140C_Hofmann1998: BOTH MFT
      routes -- the pentodiulose lane and the Hofmann & Schieberle C2+C3 lane -- have the
      SAME highest-barrier step, the `Amadori_Rearrangement` at 29.06 kcal/mol, which sits
      on the shared cysteine/ribose trunk and is the slowest step of essentially every
      route in the network. Under that rule MFT keeps exactly its old value (242.38 ppb)
      and the whole panel moves 3 rows by <1.15x, i.e. the propagator stays
      winner-takes-all in all but name.

      And the rule is WRONG, not merely inert. Take X --(slow, R_c)--> Y, then
      Y --(fast, R_1)--> P, Y --(fast, R_2)--> P, Y --(fast, R_3)--> Q. At steady state
      the trunk fixes the total flux and the branches PARTITION it by conductance:
      P's share is (1/R_1 + 1/R_2) / (1/R_1 + 1/R_2 + 1/R_3) = 2/3 for equal branches.
      This propagator's per-route weight is pool * exp(-span/RT) with
      exp(span/RT) = sum_i exp(barrier_i/RT), so when R_c dominates every route's weight
      collapses to the SAME trunk value ~ pool/exp(R_c/RT). SUMMING P's two branches then
      reproduces the 2/3 partition; taking the max returns 1/2. The shared bottleneck is
      the reason the sum is right, not a reason to suppress it -- and the mass the sum
      appears to create is not created at all, because the allocation layer normalises by
      the total activity against a fixed volatile budget.

    LIMITATION OF THE RULE THAT SHIPS, stated honestly. Every distinct producing step is
    treated as an independent branch. The model carries no pool depletion, so two branches
    drawing on the same scarce intermediate are still summed as if that intermediate were
    unlimited; the only thing standing between that and minted mass is the global budget
    cap, which bounds the TOTAL but not the split. Second, because `best_paths` keeps one
    upstream route per reactant, the step-set of a candidate is determined by its terminal
    step, so in the current enumeration this rule reduces to "one contribution per distinct
    producing step". Third, additivity is applied where flux is CONSUMED (at the product)
    and is not propagated, so an intermediate with two routes still hands its own products
    a single span: parallelism does not compound along a chain.
    """
    if not path:
        return ""
    return "|".join(str(trace.get("step_key", "")) for trace in path)


def _project_weighted_flux_to_ppb(
    steps: List[Any],
    tracked_species: Dict[str, Tuple[float, float, int, float, float]],
    best_paths: Dict[str, List[Dict[str, Any]]],
    species_catalog: Dict[str, Species],
    corrected_initial: Dict[str, float],
    target_lookup: Dict[str, Dict[str, Any]],
    exogenous_reactants: Set[str],
    temperature_kelvin: float,
    time_minutes: Optional[float],
    projection_strategy: ProjectionStrategy = DEFAULT_PROJECTION_STRATEGY,
    projection_budget: Optional[ProjectionBudget] = None,
    channel_flux_totals: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    if projection_budget is None:
        projection_budget = _estimate_projection_budget(
            corrected_initial,
            temperature_kelvin,
            time_minutes,
            strategy=projection_strategy,
        )
    if projection_budget.total_volatile_budget_molar <= 0.0:
        return {}

    severity = projection_budget.severity
    total_volatile_budget_molar = projection_budget.total_volatile_budget_molar

    projected_species = _select_accumulating_projection_species(
        tracked_species,
        species_catalog,
        target_lookup,
        exogenous_reactants,
    )
    if not projected_species:
        return {}

    # ADDITIVE FLUX (Wave S1, 2026-08-27). `channel_flux_totals[canon]` is the SUM over
    # kinetically distinct routes to `canon` of each route's Boltzmann flux, deduplicated on
    # route step-set (see `predict_from_steps`). When it is absent -- the only callers
    # that omit it are direct unit-test calls -- this falls back to the pre-Wave-S1
    # winner-takes-all flux held in `tracked_species`, so the old behaviour is still
    # reachable and testable, but it is NOT what ships.
    def _flux_for(canon: str, winner_weight: float) -> float:
        if channel_flux_totals is None:
            return winner_weight
        return float(channel_flux_totals.get(canon, winner_weight))

    candidate_entries = {
        canon: (span, depth, _flux_for(canon, weight), best_paths.get(canon, []))
        for canon, (span, conc, depth, weight, unc) in tracked_species.items()
        if canon in projected_species and depth > 0 and span < float("inf")
    }
    if not candidate_entries:
        return {}

    min_depth = min(depth for _span, depth, _weight, _path in candidate_entries.values())
    max_weight = max(max(weight, 0.0) for _canon, (_span, _depth, weight, _path) in candidate_entries.items())

    # BOLTZMANN DE-DUPLICATION (2026-08-27).
    # `weight` (tracking[canon][3]) is already a flux: reactant_flux_pool *
    # co_reactant_factor * exp(-path_span / RT) * temporal_accessibility. The previous
    # code ALSO multiplied in span_activity = exp(-(span - best_span) / (0.65 * RT))
    # and then softened the flux with a ^0.65 exponent, so the net span discrimination
    # was exp(-(1/0.65 + 0.65) * dspan / RT) = exp(-2.19 * dspan / RT) -- selectivity
    # evaluated at an effective temperature of T/2.19 (~193 K at a 150 C bake), with no
    # physical justification for either the second factor or the exponent.
    # Selectivity is now applied exactly ONCE, at the physical temperature, by using the
    # relative pathway flux directly. Everything else in the allocation (depth bias,
    # direct-sulfur bonus) is an explicit non-Boltzmann heuristic keyed on severity.
    depth_bias_strength = max(0.0, 0.85 - severity) * 1.0
    activities = {}
    for canon, (span, depth, weight, best_path) in candidate_entries.items():
        if max_weight > 0.0:
            flux_activity = max(max(weight, 0.0) / max_weight, 1.0e-12)
        else:
            flux_activity = 1.0
        depth_activity = math.exp(-depth_bias_strength * max(depth - min_depth, 0))
        terminal_family = ""
        if best_path:
            terminal_family = str(best_path[-1].get("family", "")).lower().replace("-", "_").replace(" ", "_")
        direct_sulfur_bonus = 1.0
        # AUDIT 2026-08-27 (Wave H): this equality test is now DEAD. Wave G1 renamed
        # every terminal sulfur family, and none of the survivors is spelled
        # "thiol_addition". So an unconstrained legacy heuristic silently
        # switched itself off in the same commit that changed the chemistry. Measured
        # size: at most 1.68x, and only 1.007x at the Hofmann1998 conditions where the
        # sulfur residual is 5.6x, so it is NOT the cause of the sulfur collapse.
        # Left dead ON PURPOSE rather than re-pointed: `direct_sulfur_bonus` is an
        # unconstrained fit from the quarantined-target era and is the exact knob that
        # would absorb the MFT residuals this panel exists to expose (see
        # tasks/audit_remediation.md, "NOT refit, flagged for a dedicated workstream").
        # Re-pointing it belongs to that workstream, with a refit, not to this one.
        #
        # 2026-08-27 (Wave T4) — THE FAMILY NAMES ABOVE WERE UPDATED, THE KNOB WAS
        # NOT TOUCHED. The Wave H text named `Thiol_Addition_Norfuraneol` and
        # `Thiol_Addition_Legacy_Shortcut` as "the routes that reach MFT today".
        # Neither is emitted any more: Wave N RETIRED the norfuraneol route on the
        # Cerny & Davidek isotope evidence, and `Thiol_Addition_Legacy_Shortcut` has
        # zero source literals and zero runtime emissions. The families that reach
        # MFT now are `Thiol_Addition_Pentodiulose` (Wave N) and
        # `Mercaptoketone_Cyclodehydration` (Wave P) — still neither spelled
        # "thiol_addition", so the CONCLUSION is unchanged and in fact stronger: the
        # branch is more thoroughly dead than when Wave H measured it. The knob stays
        # dead deliberately, per the paragraph above; only the two stale names here
        # were wrong.
        if terminal_family == "thiol_addition":
            direct_sulfur_bonus += 0.8 * max(0.0, 0.85 - severity)
        activities[canon] = flux_activity * depth_activity * direct_sulfur_bonus

    total_activity = sum(activities.values())
    if total_activity <= 0.0:
        return {}

    projected_ppb: Dict[str, float] = {}
    for canon, activity in activities.items():
        mol_fraction = activity / total_activity
        molar_concentration = total_volatile_budget_molar * mol_fraction
        projected_ppb[canon] = molar_concentration * _mw_from_smiles(canon) * projection_strategy.ppb_conversion_factor
    return projected_ppb


class Recommender:
    def __init__(self, results_path: Optional[Path] = None):
        self.results_path = results_path
        self.screening_data = self._load_results() if results_path else {}
        self.toxic_markers = self._load_toxic_markers()

    def predict(self, precursor_names: List[str]) -> List[str]:
        """Backward-compatible lightweight recommendation entrypoint."""
        if not precursor_names:
            return []

        from src.conditions import ReactionConditions
        from src.precursor_resolver import resolve
        from src.smirks_engine import SmirksEngine

        precursors = [resolve(name) for name in precursor_names]
        conditions = ReactionConditions()
        engine = SmirksEngine(conditions)
        steps = engine.enumerate(precursors, max_generations=4)

        heuristic_barriers: Dict[str, Tuple[float, float]] = {}
        for step in steps:
            rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
            heuristic_barriers[rxn_key] = (40.0, 5.0)

        initial_concentrations = {_canon(species.smiles): 1.0 for species in precursors}
        result = self.predict_from_steps(
            steps,
            heuristic_barriers,
            initial_concentrations,
            temperature_kelvin=conditions.temperature_kelvin,
        )

        summaries = [f"Precursors: {', '.join(precursor_names)}"]
        summaries.extend(
            f"{target['name']} ({target['type']}, span={target['span']:.2f} kcal/mol)"
            for target in result.get("targets", [])
        )
        return summaries
        
        
    def _load_yaml_db(self, filename: str) -> dict:
        data = data_access.load_yaml(data_paths.SPECIES_DIR / filename) or {}
        return {item["name"]: item for item in data.get("compounds", [])}

    def _load_toxic_markers(self) -> dict:
        return self._load_yaml_db("toxic_markers.yml")
        
    def _load_desirable(self) -> dict:
        return self._load_yaml_db("desirable_targets.yml")
        
    def _load_off_flavours(self) -> dict:
        return self._load_yaml_db("off_flavour_targets.yml")
        
    def _load_results(self) -> dict:
        if self.results_path is None or not self.results_path.exists():
            raise FileNotFoundError(
                f"Screening results not found at {self.results_path}. "
                "Please run `python scripts/run_curated_screening.py` first."
            )
            
        with open(self.results_path, "r") as f:
            data = json.load(f)
            
        # Map pathway name to span
        return {item["pathway"]: item["energetic_span_kcal"] for item in data}



    def predict_from_steps(self, 
                           steps: List[Any], 
                           barriers_dict: Dict[str, float], 
                           initial_concentrations: Dict[str, float], 
                           temperature_kelvin: float = 423.15, 
                           time_minutes: Optional[float] = None,
                           water_activity: Optional[float] = None,
                           protein_type: str = "free",
                           denaturation_state: float = 0.5,
                           fat_fraction: float = 0.0,
                           protein_fraction: float = 0.0,
                           temp_ramp_csv: Optional[str] = None,
                           protein_source: Optional[str] = None,
                           family_upstream_contract: Optional[Mapping[str, Any]] = None):
        """
        Dynamically predict active pathways given a list of generated ElementarySteps
        and their computed barriers from xTB or Hammond fallback.
        
        If temp_ramp_csv provided (SOTA), integrates propensity over the ramp.
        """
        p_type = ProteinType(protein_type)
        projection_strategy = DEFAULT_PROJECTION_STRATEGY
        ramp_data = _load_ramp(temp_ramp_csv) if temp_ramp_csv else None
        
        # Phase 7: Apply Matrix Accessibility Corrections to precursors
        # We need to map labels to concentrations for apply_matrix_correction
        # Note: apply_matrix_correction handles the scaling.
        _, corrected_initial = apply_matrix_correction(
            predicted_concentrations={}, 
            reactive_amino_acids={k: v for k, v in initial_concentrations.items()},
            protein_type=p_type,
            denaturation_state=denaturation_state,
            protein_source=protein_source
        )
        
        desirable = self._load_desirable()
        off_flavours = self._load_off_flavours()
        toxic = self._load_toxic_markers()
        
        target_lookup = {}
        for db, t_type in [(desirable, "desirable"), (off_flavours, "competing"), (toxic, "toxic")]:
            for name, data in db.items():
                if data.get("smiles"):
                    can = _canon(data["smiles"])
                    existing = target_lookup.get(can)
                    if existing is None:
                        target_lookup[can] = {
                            "name": name,
                            "type": t_type,
                            "roles": [t_type],
                            "data": data,
                        }
                        continue
                    existing_roles = set(existing.get("roles", [existing.get("type")]))
                    existing_roles.add(t_type)
                    existing["roles"] = sorted(existing_roles)
                    if existing.get("type") == "toxic" and t_type != "toxic":
                        existing["name"] = name
                        existing["type"] = t_type
                        existing["data"] = data

        species_name_lookup: Dict[str, str] = {}
        species_catalog: Dict[str, Species] = {}
        reactant_species: Set[str] = set()
        product_species: Set[str] = set()
        for step in steps:
            for species in [*step.reactants, *step.products]:
                can = _canon(species.smiles)
                if can:
                    species_catalog.setdefault(can, species)
            for species in [*step.reactants, *step.products]:
                can = _canon(species.smiles)
                if can and species.label:
                    species_name_lookup.setdefault(can, species.label)
            for species in step.reactants:
                can = _canon(species.smiles)
                if can:
                    reactant_species.add(can)
            for species in step.products:
                can = _canon(species.smiles)
                if can:
                    product_species.add(can)

        # tracking dict: canon_smiles -> (span, concentration, depth, weight, uncertainty)
        tracking = {}
        best_paths: Dict[str, List[Dict[str, Any]]] = {}
        # Pre-calculate exp(0) for initial precursors
        import math
        for s, conc in corrected_initial.items():
            canon = _canon(s)
            tracking[canon] = (0.0, conc, 0, conc * 1.0, 0.0)
            best_paths[canon] = []
            species_catalog[canon] = Species(species_name_lookup.get(canon, canon), canon)

        exogenous_reactants = set(corrected_initial)

        def _evaluate_step(step) -> Optional[Dict[str, Any]]:
            """Span / flux / trace this step confers on its products from the CURRENT `tracking`.

            Extracted verbatim from the relaxation loop body on 2026-08-27 (Wave S1) so the
            SAME arithmetic can be replayed once more after the fixpoint is reached, for the
            additive-channel pass. Returns None when a reactant is not yet reachable.
            """
            r_smiles = [r.smiles for r in step.reactants]
            p_smiles = [p.smiles for p in step.products]
            step_key = f"{'+'.join(sorted(r_smiles))}->{'+'.join(sorted(p_smiles))}"
            barrier_data = barriers_dict.get(step_key, (99.0, 5.0))
            barrier, step_unc = barrier_data if isinstance(barrier_data, tuple) else (barrier_data, 5.0)

            r_canons = [_canon(r.smiles) for r in step.reactants]
            p_canons = [_canon(p.smiles) for p in step.products]

            # The distance to fire this step is the MAX distance of all its reactants.
            # We propagate two separate notions:
            # - `conc`: user-specified precursor abundance proxy.
            # - `weight`: pathway-available reactive flux proxy.
            # Products must inherit from available flux, not directly from the original precursor pool.
            max_r_dist = 0.0
            max_r_unc = 0.0
            min_r_conc = float('inf')
            max_r_depth = 0
            dominant_reactant = None
            dominant_reactant_span = -1.0

            for r in r_canons:
                if r not in tracking:
                    return None
                r_span, r_conc, r_depth, r_weight, r_unc = tracking[r]
                max_r_dist = max(max_r_dist, r_span)
                max_r_unc = max(max_r_unc, r_unc)
                min_r_conc = min(min_r_conc, r_conc)
                max_r_depth = max(max_r_depth, r_depth)
                if r_span >= dominant_reactant_span:
                    dominant_reactant = r
                    dominant_reactant_span = r_span

            # Path properties to products via this step
            # Expert Refinement (R.7): Use sequential bottleneck (microkinetics)
            # Instead of max(barriers), we use the cumulative resistance:
            # exp(G_eff/RT) = sum(exp(G_i/RT))
            RT = 0.001987 * temperature_kelvin

            # To avoid overflow, we use the log-sum-exp trick:
            # ln(sum(exp(x_i))) = x_max + ln(sum(exp(x_i - x_max)))
            # x_max = max(max_r_dist, barrier)
            # For sequential bottleneck, we propagate uncertainty as the max of reactant/step uncertainties
            # (since they are typically dominated by the rate-limiting step's error)
            x_max = max(max_r_dist, barrier)
            path_span = x_max + RT * math.log(math.exp((max_r_dist - x_max)/RT) + math.exp((barrier - x_max)/RT))
            path_unc = max(max_r_unc, step_unc)

            path_depth = max_r_depth + 1

            # Phase G: Concentration-Aware Weighting
            # Flux = (product of reactant concs) * exp(-barrier/RT)
            # But for the cumulative pathway, we use the bottleneck span

            # Use the least-available upstream precursor/intermediate pool once.
            # The cumulative span already captures pathway resistance, so reusing an
            # already-discounted upstream weight here would double-penalize deep routes.
            reactant_flux_pool = min_r_conc
            if not math.isfinite(reactant_flux_pool):
                reactant_flux_pool = 0.0

            # Additional co-reactant availability factor keeps multi-reactant steps sensitive
            # to precursor abundance without turning units into concentration^n.
            reference_concentration = 10.0
            co_reactant_factor = 1.0
            for r in r_canons:
                if r not in exogenous_reactants:
                    continue
                normalized_conc = tracking[r][1] / (tracking[r][1] + reference_concentration)
                co_reactant_factor *= normalized_conc

            if ramp_data:
                # SOTA: Integrated Propensity
                integrated_propensity = _integrate_arrhenius(path_span, ramp_data)
                path_weight = reactant_flux_pool * co_reactant_factor * integrated_propensity
            else:
                # Path weight (Flux approximation)
                path_weight = reactant_flux_pool * co_reactant_factor * math.exp(-path_span / RT)

                # Phase Q.1: Temporal FAST Mode (Fallback)
                if time_minutes is not None:
                    # Characteristic time approx (seconds)
                    tau_sec = math.exp(path_span / RT) / get_reference_pre_exponential()
                    tau_min = tau_sec / 60.0

                    # Number of steps increases characteristic time roughly linearly
                    total_tau = tau_min * path_depth

                    # Finite-duration accessibility: slow routes scale roughly linearly
                    # when t << tau and saturate smoothly when t >> tau.
                    path_weight *= _temporal_accessibility(total_tau, time_minutes)

            # Product concentration proxy should inherit the precursor-limited pool,
            # not the exponentially discounted pathway score.
            path_conc = reactant_flux_pool
            base_path = list(best_paths.get(dominant_reactant, [])) if dominant_reactant else []
            step_trace = {
                "family": step.reaction_family or "unknown",
                "barrier": barrier,
                "path_span": path_span,
                # 2026-08-27 (Wave S1): the step's own identity, so a path can be reduced to
                # its rate-limiting step. This is the key the additive propagator dedupes on.
                "step_key": step_key,
                "reactants": [species_name_lookup.get(r, r) for r in r_canons],
                "reactant_canons": list(r_canons),
                "products": [species_name_lookup.get(p, p) for p in p_canons],
            }
            return {
                "p_canons": p_canons,
                "path_span": path_span,
                "path_conc": path_conc,
                "path_depth": path_depth,
                "path_weight": path_weight,
                "path_unc": path_unc,
                "candidate_path": base_path + [step_trace],
            }

        changed = True
        iterations = 0
        max_iterations = len(steps) + 1  # Longest possible path

        while changed and iterations < max_iterations:
            changed = False
            iterations += 1

            for step in steps:
                evaluated = _evaluate_step(step)
                if evaluated is None:
                    continue

                path_span = evaluated["path_span"]
                path_conc = evaluated["path_conc"]
                path_depth = evaluated["path_depth"]
                path_weight = evaluated["path_weight"]
                path_unc = evaluated["path_unc"]
                candidate_path = evaluated["candidate_path"]

                # Relaxation: we primarily want the lowest span path.
                for p in evaluated["p_canons"]:
                    # Update if new span is lower OR if span is same but weight is higher
                    p_key = p # Assuming p is the canonical SMILES string
                    if p_key not in tracking:
                        tracking[p_key] = (float('inf'), 0.0, 0, 0.0, 0.0) # Initialize if not present

                    current_span, current_conc, current_depth, current_weight, current_unc = tracking[p_key]

                    if path_span < current_span:
                        tracking[p_key] = (path_span, path_conc, path_depth, path_weight, path_unc)
                        best_paths[p_key] = candidate_path
                        changed = True
                    elif path_span == current_span and path_weight > current_weight:
                        tracking[p_key] = (path_span, path_conc, path_depth, path_weight, path_unc)
                        best_paths[p_key] = candidate_path
                        changed = True

        # ── ADDITIVE FLUX OVER PARALLEL CHANNELS (Wave S1, 2026-08-27) ────────────────
        # THE DEFECT THIS REPLACES. The relaxation above is winner-takes-all: `tracking[p]`
        # ends up holding the LOWEST-SPAN route to p and that route's flux alone, so the
        # allocation layer saw one channel per product. Wave P measured the consequence on
        # cys_ribose_140C_Hofmann1998: the pentodiulose MFT lane alone predicted 242.38 ppb,
        # the Hofmann & Schieberle C2+C3 lane alone 71.02 ppb, and BOTH together 242.38 ppb
        # -- adding a real, literature-evidenced second route to the flagship compound moved
        # the flagship number by EXACTLY ZERO. Real kinetics sums parallel channels.
        #
        # WHAT IS SUMMED. One more sweep over the steps, from the CONVERGED tracking state
        # (not accumulated during relaxation -- that would add the same channel once per
        # iteration). Each step contributes the Boltzmann-weighted flux it confers on its
        # products, exactly the `path_weight` the relaxation already computes: the
        # rate-limiting-span Boltzmann factor times the available upstream pool.
        #
        # THE DEDUPE RULE. Channel identity is the route's FULL ORDERED STEP-SET, and within a
        # channel the largest flux wins rather than the sum, so no route can be counted twice.
        # `_route_channel_id` carries the full argument, including the rate-limiting-step rule
        # that was implemented, MEASURED and rejected: both MFT routes here share the trunk
        # `Amadori_Rearrangement` as their slowest step, so that rule left the flagship number
        # unmoved AND under-counts the conductance partition it was meant to protect.
        #
        # WHAT IS *NOT* CHANGED. Span, depth, concentration proxy and `best_paths` still come
        # from the lowest-span representative route, so every span-driven behaviour
        # downstream (ranking by span, path traces, `depth_activity`) is bit-identical. Only
        # the flux proxy the budget allocation reads becomes additive. Consequence, and it is
        # a real limitation: additivity is applied where the flux is CONSUMED (at the product,
        # by the allocation layer), not propagated -- an INTERMEDIATE with two routes still
        # hands its downstream products a single span, so parallelism does not compound along
        # a chain.
        channel_flux: Dict[str, Dict[str, float]] = {}
        for step in steps:
            evaluated = _evaluate_step(step)
            if evaluated is None:
                continue
            channel_id = _route_channel_id(evaluated["candidate_path"])
            weight = float(evaluated["path_weight"])
            for p in evaluated["p_canons"]:
                by_channel = channel_flux.setdefault(p, {})
                if weight > by_channel.get(channel_id, float("-inf")):
                    by_channel[channel_id] = weight
        summed_channel_flux: Dict[str, float] = {
            p: float(sum(by_channel.values())) for p, by_channel in channel_flux.items()
        }

        # Project ranked FAST outputs onto a bounded volatile ppb budget.
        projection_budget = _estimate_projection_budget(
            corrected_initial,
            temperature_kelvin,
            time_minutes,
            strategy=projection_strategy,
        )
        raw_concentrations = _project_weighted_flux_to_ppb(
            steps,
            tracking,
            best_paths,
            species_catalog,
            corrected_initial,
            target_lookup,
            exogenous_reactants,
            temperature_kelvin,
            time_minutes,
            projection_strategy=projection_strategy,
            projection_budget=projection_budget,
            channel_flux_totals=summed_channel_flux,
        )

        # Ensure injected targets (e.g. Hexanal from lipid oxidation) are included in raw_concentrations
        for canon, conc in corrected_initial.items():
            if canon in target_lookup and canon not in raw_concentrations:
                raw_concentrations[canon] = conc

        observable_volatiles, projection_metadata = _apply_output_projection(
            raw_concentrations,
            species_catalog,
            target_lookup,
            temperature_kelvin,
            protein_type=protein_type,
            time_minutes=time_minutes,
            water_activity=water_activity,
            denaturation_state=denaturation_state,
            fat_fraction=fat_fraction,
            protein_fraction=protein_fraction,
            projection_budget=projection_budget,
            projection_strategy=projection_strategy,
            species_name_lookup=species_name_lookup,
            protein_source=protein_source,
            family_upstream_contract=family_upstream_contract,
        )
        proxy_volatiles = dict(raw_concentrations)
        
        # Add names to the dictionary for downstream sensory and benchmark matching
        final_volatiles = {}
        final_proxy_volatiles = {}
        for p_canon, conc in observable_volatiles.items():
            final_volatiles[p_canon] = conc
            final_proxy_volatiles[p_canon] = proxy_volatiles.get(p_canon, conc)
            species_name = species_name_lookup.get(p_canon)
            if species_name:
                final_volatiles[species_name] = conc
                final_proxy_volatiles[species_name] = proxy_volatiles.get(p_canon, conc)
            t_info = target_lookup.get(p_canon)
            if t_info:
                final_volatiles[t_info["name"]] = conc
                final_proxy_volatiles[t_info["name"]] = proxy_volatiles.get(p_canon, conc)

        # Identify which targets were produced
        active_pathways = []
        for p_canon, (span, conc, depth, weight, unc) in tracking.items():
            t_info = target_lookup.get(p_canon)
            if t_info and span < float('inf') and span >= 0.0:
                species = species_catalog.get(p_canon, Species(t_info["name"], p_canon))
                observability = _headspace_observability_metadata(species, target_lookup)
                

        # Update active_pathways with final ppb
        active_pathways = []
        active_targets = {
            p_canon: t_info
            for p_canon, t_info in target_lookup.items()
            if p_canon in tracking and tracking[p_canon][0] < float('inf') and tracking[p_canon][0] >= 0.0
        }

        class MockTarget:
            def __init__(self, label):
                self.label = label

        for p_canon, t_info in active_targets.items():
            _, _, depth, _, _ = tracking.get(p_canon, (0, 0, 0, 0, 0))
            span = tracking[p_canon][0]
            
            species = species_catalog.get(p_canon, Species(t_info["name"], p_canon))
            observability = _headspace_observability_metadata(species, target_lookup)
            
            p_dict = {
                "name": t_info["name"],
                "span": span,
                "concentration": observable_volatiles.get(p_canon, 0.0), # Use corrected
                "proxy_concentration": final_proxy_volatiles.get(p_canon, 0.0),
                "weighted_flux": observable_volatiles.get(p_canon, 0.0), # Observable output budget proxy for ranking tables
                "span_uncertainty": tracking[p_canon][4],
                "depth": depth,
                "target": MockTarget(t_info["name"]),
                "type": t_info["type"],
                "penalty": "LOW",
                "toxicity": None,
                "sensory": t_info["data"].get("sensory_desc", "-"),
                "threshold": t_info["data"].get("odour_threshold_ug_per_kg", None),
                "roles": t_info.get("roles", [t_info["type"]]),
                "projection": projection_metadata.get(p_canon, {}),
                **observability,
            }
            
            if t_info["type"] == "toxic":
                p_dict["toxicity"] = {
                    "name": t_info["name"],
                    "risk": t_info["data"].get("health_risk", "Unknown"),
                    "priority": t_info["data"].get("priority", "high").upper()
                }
            
            active_pathways.append(p_dict)
            
        # Sort targets: lower span first, then higher concentration
        active_pathways.sort(key=lambda x: (x["span"], -x["concentration"]))

        # --- Advanced Diagnostics: Precursor Attribution ---
        precursor_attribution = {} # {precursor_name: contribution_score}
        # We estimate contribution based on how many targets a precursor "feeds" 
        # weighted by the target's observable concentration.
        for p_canon, t_info in active_targets.items():
            path = best_paths.get(p_canon, [])
            target_conc = observable_volatiles.get(p_canon, 0.0)
            if target_conc <= 0: continue
            
            # Find all precursors in this path
            path_precursors = set()
            for step_tr in path:
                for r_canon in step_tr.get("reactant_canons", []):
                    # Initial precursors have span 0.0 in tracking
                    if r_canon in tracking and tracking[r_canon][0] == 0.0:
                        path_precursors.add(r_canon)
            
            if path_precursors:
                split_conc = target_conc / len(path_precursors)
                for prec_can in path_precursors:
                    p_name = species_name_lookup.get(prec_can, prec_can)
                    precursor_attribution[p_name] = precursor_attribution.get(p_name, 0.0) + split_conc

        # --- Advanced Diagnostics: Suppressed Compounds ---
        suppressed_compounds = []
        for canon, meta in projection_metadata.items():
            proxy = meta.get("proxy_ppb", 0.0)
            obs = meta.get("observable_ppb", 0.0)
            if proxy > 5.0 and (obs / proxy) < 0.3: # Significant suppression
                suppressed_compounds.append({
                    "name": species_name_lookup.get(canon, canon),
                    "proxy_ppb": proxy,
                    "observable_ppb": obs,
                    "reduction_factor": 1.0 - (obs / proxy),
                    "primary_cause": "headspace" if meta.get("headspace_factor", 1.0) < meta.get("matrix_factor", 1.0) else "matrix"
                })
        suppressed_compounds.sort(key=lambda x: x["reduction_factor"], reverse=True)
            
        # ── PBMA Metrics: Lipid Trapping Efficiency ──
        # Find which initial pool members are lipids
        lipid_pool_canons = []
        lysine_can = _canon("NCCCCC(N)C(=O)O")
        has_lysine = lysine_can in tracking

        for can in initial_concentrations.keys():
            if can in target_lookup and target_lookup[can]["name"] in off_flavours:
                lipid_pool_canons.append(can)

        trapping_results = {}
        for lipid_can in lipid_pool_canons:
            lipid_name = target_lookup[lipid_can]["name"]
            
            # Find all Schiff bases derived from this lipid
            sb_weights = 0.0
            for step in steps:
                if step.reaction_family == "Lipid_Schiff_Base":
                    step_r_canons = [_canon(r.smiles) for r in step.reactants]
                    if lipid_can in step_r_canons:
                        # Path barrier for this step
                        max_r_dist = 0.0
                        reachable = True
                        for r_smi in [r.smiles for r in step.reactants]:
                            rc = _canon(r_smi)
                            if rc not in tracking:
                                reachable = False
                                break
                            max_r_dist = max(max_r_dist, tracking[rc][0])
                        
                        if reachable:
                            step_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                            barrier_data = barriers_dict.get(step_key, (99.0, 5.0))
                            barrier = barrier_data[0] if isinstance(barrier_data, tuple) else barrier_data
                            path_barrier = max(max_r_dist, barrier)
                            
                            if ramp_data:
                                sb_weights += _integrate_arrhenius(path_barrier, ramp_data)
                            else:
                                sb_weights += _weight(path_barrier, temperature_kelvin)
            
            if ramp_data:
                persistence_w = _integrate_arrhenius(30.0, ramp_data)
            else:
                persistence_w = _weight(30.0, temperature_kelvin)
            if sb_weights + persistence_w > 0:
                eff = 100.0 * sb_weights / (persistence_w + sb_weights)
            else:
                eff = 0.0
            trapping_results[lipid_name] = eff

        # ── PBMA Metrics: Lysine Budget (DHA Competition) ──
        lysine_budget = 0.0
        if has_lysine:
            w_maillard = 0.0
            w_dha = 0.0
            
            for step in steps:
                step_r_canons = [_canon(r.smiles) for r in step.reactants]
                if lysine_can in step_r_canons:
                    # Path barrier — must re-initialise per step
                    max_r_dist = 0.0
                    reachable = True
                    for r_smi in [r.smiles for r in step.reactants]:
                        rc = _canon(r_smi)
                        if rc not in tracking:
                            reachable = False
                            break
                        max_r_dist = max(max_r_dist, tracking[rc][0])
                    
                    if not reachable: 
                        continue
                    
                    step_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                    barrier_data = barriers_dict.get(step_key, (99.0, 5.0))
                    barrier = barrier_data[0] if isinstance(barrier_data, tuple) else barrier_data
                    path_barrier = max(max_r_dist, barrier)
                    
                    if ramp_data:
                        weight = _integrate_arrhenius(path_barrier, ramp_data)
                    else:
                        weight = _weight(path_barrier, temperature_kelvin)
                    
                    if step.reaction_family in ["Schiff_Base_Formation", "Lipid_Schiff_Base"]:
                        w_maillard += weight
                    elif step.reaction_family == "DHA_Crosslinking":
                        w_dha += weight
            
            if w_maillard + w_dha > 0:
                lysine_budget = 100.0 * w_dha / (w_maillard + w_dha)
            else:
                lysine_budget = 0.0

        metrics = {
            "trapping_efficiency": trapping_results,
            "lysine_budget_dha": lysine_budget,
            "bottleneck": {
                "precursor": species_name_lookup.get(projection_budget.limiting_precursor_name, projection_budget.limiting_precursor_name),
                "severity": float(projection_budget.severity),
                "load_factor": float(projection_budget.load_factor)
            },
            "precursor_attribution": precursor_attribution,
            "suppressed_compounds": suppressed_compounds,
            "thiamine_pathway_active": any("thiamine" in str(name).lower() for name in species_name_lookup.values()),
        }

        return {
            "targets": active_pathways,
            "metrics": metrics,
            "predicted_ppb": final_volatiles,
            "predicted_proxy_ppb": final_proxy_volatiles,
            "projection_metadata": projection_metadata,
            "projection_context": {
                "limiting_precursor_molar": float(projection_budget.limiting_precursor_molar),
                "projection_load_factor": float(projection_budget.load_factor),
                "projection_temperature_factor": float(projection_budget.temperature_factor),
                "projection_time_factor": float(projection_budget.time_factor),
                "projection_severity": float(projection_budget.severity),
                "projection_kinetic_drive": float(projection_budget.kinetic_drive),
                "projection_conversion_extent": float(projection_budget.conversion_extent),
                "water_activity": None if water_activity is None else float(water_activity),
                "volatile_yield_fraction": float(projection_budget.volatile_yield_fraction),
                "total_volatile_budget_molar": float(projection_budget.total_volatile_budget_molar),
                **_projection_strategy_metadata(projection_strategy),
            },
            "debug_paths": best_paths,
            # 2026-08-27 (Wave S1): per-product flux BROKEN OUT BY CHANNEL, so "the model has
            # N routes to X" is auditable against what the allocation layer actually summed.
            # Keys are canonical SMILES. Inner keys are CHANNEL IDS: the full ordered
            # step-set of the route, joined on "|" (`_route_channel_id`), which is the
            # key the additive propagator deduplicates on.
            #
            # 2026-08-27 (Wave T4) CORRECTION. This line used to read "inner keys are
            # rate-limiting `step_key`s". That is the description of the REJECTED rule
            # — "routes sharing a rate-limiting step are one channel, take the max" —
            # which `_route_channel_id`'s own docstring records as implemented,
            # MEASURED and REJECTED (both MFT routes share the trunk Amadori step as
            # their slowest, so it left the flagship number unmoved and under-counts
            # the conductance partition a shared bottleneck actually produces). A
            # reader auditing "the model has N routes to X" against this payload would
            # have mis-read every key. No behaviour change; the payload was always the
            # full step-set.
            "debug_channel_flux": channel_flux,
            "species_names": species_name_lookup,
        }



