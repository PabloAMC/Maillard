"""
src/projection.py — Projection budget, strategy helpers, and output-projection logic.

This module centralises:
  - ProjectionBudget / ProjectionStrategy dataclasses and budget estimation helpers
    (originally here)
  - All output-projection functions extracted from ``src/recommend.py`` as part
    of the Priority-2 monolith decomposition (P2.2).

``src/recommend.py`` re-exports every public symbol from this module so that
downstream ``from src.recommend import X`` calls continue to work unchanged.
New code should import directly from this module.
"""
from __future__ import annotations

import math
import yaml
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]

# ── Optional RDKit ────────────────────────────────────────────────────────────
try:
    from rdkit import Chem
except ImportError:
    Chem = None  # type: ignore[assignment]

if Chem is not None:
    _CARBOXYLIC_ACID_SMARTS = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    _PRIMARY_AMINE_SMARTS = Chem.MolFromSmarts("[NX3;H1,H2;!$(NC=O)]")
    _IMINE_SMARTS = Chem.MolFromSmarts("[CX3]=[NX2;!R]")
else:
    _CARBOXYLIC_ACID_SMARTS = None
    _PRIMARY_AMINE_SMARTS = None
    _IMINE_SMARTS = None

# ── Domain imports (none of these import from recommend/projection → no cycles) ──
from src.chem_utils import canonicalize_smiles
from src.headspace import HeadspaceModel
from src.matrix_calibration_registry import (
    describe_matrix_calibration,
    determine_matrix_process_state,
)
from src.matrix_correction import (
    ProteinType,
    classify_accessibility_state,
    classify_volatile_matrix_family,
    describe_compound_matrix_retention,
    resolve_matrix_correction,
)
from src.matrix_targets import get_compound_panel_entry
from src.pathway_extractor import Species
from src.projection_metadata import ProjectionMetadataMap, make_projection_metadata_row

@dataclass(frozen=True)
class ProjectionBudget:
    limiting_precursor_molar: float
    load_factor: float
    temperature_factor: float
    time_factor: float
    severity: float
    volatile_yield_fraction: float
    total_volatile_budget_molar: float
    limiting_precursor_name: str


@dataclass(frozen=True)
class ProjectionStrategy:
    name: str
    precursor_concentration_unit: str
    limiting_pool_to_molar_factor: float
    baseline_volatile_yield_fraction: float
    severity_volatile_yield_slope: float
    ppb_conversion_factor: float
    ppb_basis: str
    notes: str


DEFAULT_PROJECTION_STRATEGY = ProjectionStrategy(
    name="precursor_limited_observable_v1",
    precursor_concentration_unit="mM",
    limiting_pool_to_molar_factor=1.0e-3,
    baseline_volatile_yield_fraction=1.0e-6,
    severity_volatile_yield_slope=1.5e-3,
    ppb_conversion_factor=1.0e6,
    ppb_basis="aqueous_mass_equivalent_ppb",
    notes=(
        "Allocates a conservative volatile budget from the limiting precursor pool, then converts "
        "M to ppb via MW assuming dilute aqueous density (~1 kg/L) before matrix/headspace projection."
    ),
)

def _thermal_severity(temperature_kelvin: float, time_minutes: Optional[float]) -> float:
    return _projection_temperature_factor(temperature_kelvin) * _projection_time_factor(time_minutes)


def _projection_temperature_factor(temperature_kelvin: float) -> float:
    import math

    temp_c = temperature_kelvin - 273.15
    return 1.0 / (1.0 + math.exp(-(temp_c - 110.0) / 18.0))


def _projection_time_factor(time_minutes: Optional[float]) -> float:
    import math

    if time_minutes is None:
        return 1.0
    return 1.0 - math.exp(-max(time_minutes, 0.0) / 25.0)


def _estimate_projection_budget(
    corrected_initial: Dict[str, float],
    temperature_kelvin: float,
    time_minutes: Optional[float],
    strategy: ProjectionStrategy = DEFAULT_PROJECTION_STRATEGY,
) -> ProjectionBudget:
    if not corrected_initial:
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none")

    # Find limiting precursor (minimum positive molarity)
    positive_items = [(k, float(v)) for k, v in corrected_initial.items() if float(v) > 0.0]
    if not positive_items:
        return ProjectionBudget(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "none")

    limiting_name, limiting_val = min(positive_items, key=lambda x: x[1])
    limiting_precursor_molar = limiting_val * strategy.limiting_pool_to_molar_factor
    
    load_factor = _relative_precursor_load_factor(corrected_initial)
    temperature_factor = _projection_temperature_factor(temperature_kelvin)
    time_factor = _projection_time_factor(time_minutes)
    severity = temperature_factor * time_factor
    volatile_yield_fraction = (
        strategy.baseline_volatile_yield_fraction
        + strategy.severity_volatile_yield_slope * severity
    )
    total_volatile_budget_molar = limiting_precursor_molar * volatile_yield_fraction * max(load_factor, 0.0)

    return ProjectionBudget(
        limiting_precursor_molar=float(limiting_precursor_molar),
        load_factor=float(max(load_factor, 0.0)),
        temperature_factor=float(temperature_factor),
        time_factor=float(time_factor),
        severity=float(severity),
        volatile_yield_fraction=float(volatile_yield_fraction),
        total_volatile_budget_molar=float(max(total_volatile_budget_molar, 0.0)),
        limiting_precursor_name=limiting_name
    )


def _temporal_accessibility(total_tau_minutes: float, time_minutes: Optional[float]) -> float:
    import math

    if time_minutes is None:
        return 1.0
    if time_minutes <= 0.0:
        return 0.0
    if total_tau_minutes <= 0.0:
        return 1.0
    return 1.0 - math.exp(-time_minutes / total_tau_minutes)


def _relative_precursor_load_factor(corrected_initial: Dict[str, float]) -> float:
    import math

    positive_values = [max(float(value), 0.0) for value in corrected_initial.values() if float(value) > 0.0]
    if not positive_values:
        return 0.0

    limiting_value = min(positive_values)
    if limiting_value <= 0.0:
        return 0.0

    normalized = [max(value / limiting_value, 1.0e-12) for value in positive_values]
    return math.exp(sum(math.log(value) for value in normalized) / len(normalized))


def _projection_strategy_metadata(strategy: ProjectionStrategy) -> Dict[str, Any]:
    return {
        "projection_strategy_name": strategy.name,
        "projection_precursor_unit": strategy.precursor_concentration_unit,
        "projection_ppb_basis": strategy.ppb_basis,
        "projection_limiting_pool_to_molar_factor": float(strategy.limiting_pool_to_molar_factor),
        "projection_baseline_volatile_yield_fraction": float(strategy.baseline_volatile_yield_fraction),
        "projection_severity_volatile_yield_slope": float(strategy.severity_volatile_yield_slope),
        "projection_ppb_conversion_factor": float(strategy.ppb_conversion_factor),
        "projection_strategy_notes": strategy.notes,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Functions below were extracted from src/recommend.py (P2.2 decomposition).
# src/recommend.py re-exports every symbol so existing imports continue to work.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Chemical name normalisation ───────────────────────────────────────────────

def _normalize_chemical_name(name: str) -> str:
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())


def _canon(smi: str) -> str:
    return canonicalize_smiles(smi, fallback_to_original=True, strip_salts=True)


# ── Henry's Law / headspace lookups ──────────────────────────────────────────

_HENRY_CONSTANTS_PATH = ROOT / "data" / "lit" / "henry_constants.yml"
_NON_OBSERVABLE_KAW_THRESHOLD = 1.0e-8
_LOW_HEADSPACE_REFERENCE_KAW = 1.0e-5


def _load_henry_lookup() -> Dict[str, Dict[str, Any]]:
    if not _HENRY_CONSTANTS_PATH.exists():
        return {}
    with open(_HENRY_CONSTANTS_PATH, "r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    constants = raw.get("constants", [])
    lookup: Dict[str, Dict[str, Any]] = {}
    for entry in constants:
        if not entry.get("name"):
            continue
        lookup[_normalize_chemical_name(entry["name"])] = entry
    return lookup


_HENRY_LOOKUP: Dict[str, Dict[str, Any]] = _load_henry_lookup()
_HEADSPACE_MODEL = HeadspaceModel(str(_HENRY_CONSTANTS_PATH))

# ── RDKit-backed chemical helpers ─────────────────────────────────────────────

_BUDGET_EXCLUDED_CANONICAL: Set[str] = {
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


def _mw_from_smiles(smiles: str) -> float:
    if Chem is None:
        return 100.0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 100.0
    try:
        from rdkit.Chem import Descriptors  # type: ignore[import]
        return float(Descriptors.MolWt(mol))
    except Exception:
        return 100.0


# ── Henry / headspace helpers ─────────────────────────────────────────────────

def _henry_entry_for_species(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
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


def _is_observable_target_species(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
) -> bool:
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
            baseline_air = float(
                _HEADSPACE_MODEL.predict_headspace(
                    {compound_name: 1.0}, temp_c, fat_fraction=0.0, protein_fraction=0.0
                ).get(compound_name, 1.0)
            )
            matrix_air = float(
                _HEADSPACE_MODEL.predict_headspace(
                    {compound_name: 1.0}, temp_c, fat_fraction=fat_fraction, protein_fraction=protein_fraction
                ).get(compound_name, baseline_air)
            )
            if baseline_air > 0.0:
                matrix_factor = max(1.0e-3, min(1.0, matrix_air / baseline_air))
        except Exception:
            matrix_factor = 1.0

    # Keep the penalty conservative so we do not destabilize the validated
    # free-amino-acid benchmarks while still reflecting that near-nonvolatile
    # species should not consume the same observable headspace budget.
    return intrinsic_factor * matrix_factor


# ── Melanoidin trapping ───────────────────────────────────────────────────────

_MELANOIDIN_TRAPPING_PROFILES: Dict[str, Dict[str, float]] = {
    _normalize_chemical_name("2-methyl-3-furanthiol"): {"slope": 0.85, "floor": 0.20},
    _normalize_chemical_name("2-furfurylthiol"): {"slope": 1.10, "floor": 0.08},
    _normalize_chemical_name("bis(2-methyl-3-furyl) disulfide"): {"slope": 0.55, "floor": 0.35},
}


def _resolve_melanoidin_trapping_factor(
    compound_name: str,
    *,
    protein_type: ProteinType,
    process_state: str,
    projection_severity: float,
) -> float:
    if protein_type == ProteinType.FREE_AMINO_ACID:
        return 1.0

    profile = _MELANOIDIN_TRAPPING_PROFILES.get(_normalize_chemical_name(compound_name))
    if profile is None:
        return 1.0

    severity = max(0.0, min(1.0, float(projection_severity)))
    factor = 1.0 - float(profile["slope"]) * severity
    if process_state == "heated_matrix":
        factor *= 0.92
    return max(float(profile["floor"]), min(1.0, factor))


# ── Matrix output context ─────────────────────────────────────────────────────

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


# ── Budget selection helpers ──────────────────────────────────────────────────

def _is_budget_relevant_species(
    species: Species,
    target_lookup: Dict[str, Dict[str, Any]],
) -> bool:
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
    steps: List[Any],
    tracked_species: Dict[str, Tuple[float, float, int, float, float]],
    species_catalog: Dict[str, Species],
    target_lookup: Dict[str, Dict[str, Any]],
    exogenous_reactants: Set[str],
    downstream_margin_kcal: float = 0.25,
) -> Set[str]:
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
        steps,
        tracked_species,
        species_catalog,
        target_lookup,
        exogenous_reactants,
    )
    if not projected_species:
        return {}

    candidate_entries = {
        canon: (span, depth, weight, best_paths.get(canon, []))
        for canon, (span, conc, depth, weight, unc) in tracked_species.items()
        if canon in projected_species and depth > 0 and span < float("inf")
    }
    if not candidate_entries:
        return {}

    best_span = min(span for span, _depth, _weight, _path in candidate_entries.values())
    min_depth = min(depth for _span, depth, _weight, _path in candidate_entries.values())
    span_window_kcal = max(0.35, 0.65 * 0.001987 * temperature_kelvin)
    max_weight = max(max(weight, 0.0) for _canon, (_span, _depth, weight, _path) in candidate_entries.items())

    # At lower thermal severity, short terminal routes should retain a mild advantage over
    # deeper ones when the final ppb budget is allocated across competing outputs.
    depth_bias_strength = max(0.0, 0.85 - severity) * 1.0
    activities: Dict[str, float] = {}
    for canon, (span, depth, weight, best_path) in candidate_entries.items():
        span_activity = math.exp(-(span - best_span) / span_window_kcal)
        if max_weight > 0.0:
            relative_weight = max(weight, 0.0) / max_weight
            flux_activity = max(relative_weight, 1.0e-6) ** 0.65
        else:
            flux_activity = 1.0
        depth_activity = math.exp(-depth_bias_strength * max(depth - min_depth, 0))
        terminal_family = ""
        if best_path:
            terminal_family = str(best_path[-1].get("family", "")).lower().replace("-", "_").replace(" ", "_")
        direct_sulfur_bonus = 1.0
        if terminal_family == "thiol_addition":
            direct_sulfur_bonus += 0.8 * max(0.0, 0.85 - severity)
        activities[canon] = span_activity * flux_activity * depth_activity * direct_sulfur_bonus

    total_activity = sum(activities.values())
    if total_activity <= 0.0:
        return {}

    projected_ppb: Dict[str, float] = {}
    for canon, activity in activities.items():
        mol_fraction = activity / total_activity
        molar_concentration = total_volatile_budget_molar * mol_fraction
        projected_ppb[canon] = molar_concentration * _mw_from_smiles(canon) * projection_strategy.ppb_conversion_factor
    return projected_ppb


# ── Full observable output projection ─────────────────────────────────────────

def _apply_output_projection(
    raw_concentrations: Dict[str, float],
    species_catalog: Dict[str, Species],
    target_lookup: Dict[str, Dict[str, Any]],
    temperature_kelvin: float,
    protein_type: str,
    process_state: Optional[str] = None,
    time_minutes: Optional[float] = None,
    water_activity: Optional[float] = None,
    denaturation_state: float = 0.5,
    fat_fraction: float = 0.0,
    protein_fraction: float = 0.0,
    projection_budget: Optional[ProjectionBudget] = None,
    projection_strategy: ProjectionStrategy = DEFAULT_PROJECTION_STRATEGY,
    species_name_lookup: Optional[Dict[str, str]] = None,
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
    resolved_process_state = process_state or determine_matrix_process_state(
        temperature_celsius=temperature_kelvin - 273.15,
        time_minutes=float(time_minutes or 60.0),
        water_activity=water_activity,
    )
    accessibility_state = classify_accessibility_state(
        protein_type,
        denaturation_state,
        dominant_source="denaturation_state_arg",
    )

    budget_metadata: Dict[str, Any] = {}
    if projection_budget is not None:
        budget_metadata = {
            "limiting_precursor_molar": float(projection_budget.limiting_precursor_molar),
            "projection_load_factor": float(projection_budget.load_factor),
            "projection_temperature_factor": float(projection_budget.temperature_factor),
            "projection_time_factor": float(projection_budget.time_factor),
            "projection_severity": float(projection_budget.severity),
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
            compound_name = species_name_lookup.get(canon, canon) if species_name_lookup else canon
            calibration = describe_matrix_calibration(
                compound_name,
                protein_type=protein_type,
                process_state=resolved_process_state,
            )
            panel_entry = get_compound_panel_entry(compound_name) or {}
            calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
            melanoidin_factor = _resolve_melanoidin_trapping_factor(
                compound_name,
                protein_type=p_type,
                process_state=process_state,
                projection_severity=projection_severity,
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
                    "process_state": resolved_process_state,
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
            retention_description: Dict[str, Any] = {}
        else:
            retention_description = describe_compound_matrix_retention(
                species.label or canon,
                protein_type=p_type,
                denaturation_state=denaturation_state,
                smiles=species.smiles,
                temperature_celsius=temperature_kelvin - 273.15,
                time_minutes=time_minutes,
                water_activity=water_activity,
                process_state=resolved_process_state,
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
        compound_name = species.label or canon
        calibration = describe_matrix_calibration(
            compound_name,
            protein_type=protein_type,
            process_state=resolved_process_state,
        )
        panel_entry = get_compound_panel_entry(compound_name) or {}
        calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
        melanoidin_factor = _resolve_melanoidin_trapping_factor(
            compound_name,
            protein_type=p_type,
            process_state=process_state,
            projection_severity=projection_severity,
        )
        observable_value = raw_value * effective_matrix_factor * headspace_factor * calibration_factor * melanoidin_factor
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
                "browning_index": float(projection_severity),
                "browning_narrative": (
                    "melanoidin-linked sulfur trapping surrogate"
                    if melanoidin_factor < 1.0
                    else "no explicit browning-linked sulfur penalty"
                ),
                "volatile_class": classify_volatile_matrix_family(compound_name, smiles=species.smiles),
                "process_state": resolved_process_state,
                "accessibility_profile": accessibility_state.profile,
                "accessibility_warning": accessibility_state.accessibility_warning,
                "accessibility_dominant_source": accessibility_state.dominant_source,
                **panel_entry,
                **calibration,
                **budget_metadata,
            },
        )

    return observable_ppb, projection_metadata

