#!/usr/bin/env python3
"""
src/recommend.py — Maillard Reaction Pathway Recommender

Tier 4 Pipeline Integration:
Loads the Tier 1 xTB screening results and matches them against user-defined
or canonical precursors to recommend actionable formulation adjustments.
"""

import sys
import json
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Set, Optional, Any, Tuple

from data.reactions.curated_pathways import PATHWAYS, PATHWAY_METADATA
from src.matrix_correction import ProteinType, apply_matrix_correction
from src.pathway_extractor import Species
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

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_HENRY_CONSTANTS_PATH = ROOT / "data" / "lit" / "henry_constants.yml"
_NON_OBSERVABLE_KAW_THRESHOLD = 1.0e-8


def _normalize_chemical_name(name: str) -> str:
    return " ".join(str(name).lower().replace("_", " ").replace("-", " ").split())


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


_HENRY_LOOKUP = _load_henry_lookup()

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


@dataclass
class PrecursorSystem:
    name: str
    precursors: List[str]
    notes: str


SYSTEMS = [
    PrecursorSystem(
        "Ribose + Cysteine (Savory Base)",
        ["D-ribose", "L-cysteine"],
        "Classic model system for meaty flavors."
    ),
    PrecursorSystem(
        "Glucose + Glycine (Baked Base)",
        ["D-glucose", "glycine"],
        "Classic model for baked/roasted notes."
    ),
    PrecursorSystem(
        "Ribose + Cysteine + Leucine",
        ["D-ribose", "L-cysteine", "L-leucine"],
        "Complex system targeting Strecker aldehydes."
    ),
    PrecursorSystem(
        "Plant-Based Deficient (Glucose + Lysine + Hexanal)",
        ["D-glucose", "L-lysine", "hexanal"],
        "Mimics a legume base undergoing lipid oxidation."
    ),
    PrecursorSystem(
        "Ribose + Cysteine + Lysine (DHA Penaly Test)",
        ["D-ribose", "L-cysteine", "L-lysine"],
        "Tests if the DHA cross-linking pathway penalises the FFT pathway."
    )
]


# Build canonical SMILES lookup for targets
def _canon(smi):
    if Chem is None:
        return smi
    try:
        can = set(Chem.MolToSmiles(Chem.MolFromSmiles(smi)).split("."))
        # just return the largest fragment if disconnected
        return max(can, key=len)
    except ImportError:
        # RDKit not available, return original SMILES
        return smi
    except Exception as e:
        print(f"Warning: RDKit conformer generation failed: {e}")
        return smi

def _weight(barrier_kcal, temp_kelvin=423.15): # Default 150C
    import math
    if barrier_kcal >= 99.0: 
        return 0.0
    R = 0.001987
    return math.exp(-barrier_kcal / (R * temp_kelvin))

def _integrate_arrhenius(barrier_kcal: float, 
                         ramp_data: List[Tuple[float, float]], 
                         pre_exponential: float = 1e11) -> float:
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
    
    R = 0.001987 # kcal/mol/K
    k_fine = pre_exponential * np.exp(-barrier_kcal / (R * T_fine))
    
    # Integrate using Trapezoidal rule (standard for discrete steps)
    integral_k_dt = np.trapezoid(k_fine, t_fine)
    
    # Yield approximation (saturation handling)
    return 1.0 - np.exp(-integral_k_dt)

def _load_ramp(ramp_path: str) -> List[Tuple[float, float]]:
    """Load a temperature ramp from CSV (time_min, temp_c)."""
    try:
        df = pd.read_csv(ramp_path)
        # Standard schema: 'time' in minutes, 'temp' in Celsius
        if 'time' not in df.columns or 'temp' not in df.columns:
            return []
        return list(zip(df['time'], df['temp'] + 273.15))
    except Exception as e:
        print(f"Warning: Failed to load ramp {ramp_path}: {e}")
        return []

def _weight(barrier_kcal, temp_kelvin=423.15): # Default 150C
    import math
    if barrier_kcal >= 99.0: 
        return 0.0
    R = 0.001987
    return math.exp(-barrier_kcal / (R * temp_kelvin))


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


def _thermal_severity(temperature_kelvin: float, time_minutes: Optional[float]) -> float:
    import math

    temp_c = temperature_kelvin - 273.15
    temp_factor = 1.0 / (1.0 + math.exp(-(temp_c - 110.0) / 18.0))
    if time_minutes is None:
        time_factor = 1.0
    else:
        time_factor = 1.0 - math.exp(-max(time_minutes, 0.0) / 25.0)
    return temp_factor * time_factor


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
) -> Dict[str, float]:
    import math

    if not corrected_initial:
        return {}

    limiting_precursor_molar = max(min(corrected_initial.values()), 0.0) / 1000.0
    if limiting_precursor_molar <= 0.0:
        return {}

    load_factor = _relative_precursor_load_factor(corrected_initial)
    if load_factor <= 0.0:
        return {}

    severity = _thermal_severity(temperature_kelvin, time_minutes)
    volatile_yield_fraction = 1.0e-6 + 1.5e-3 * severity
    total_volatile_budget_molar = limiting_precursor_molar * volatile_yield_fraction * load_factor

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
    activities = {}
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
        projected_ppb[canon] = molar_concentration * _mw_from_smiles(canon) * 1e6
    return projected_ppb


class Recommender:
    def __init__(self, results_path: Optional[Path] = None):
        self.results_path = results_path
        self.screening_data = self._load_results() if results_path else {}
        self.toxic_markers = self._load_toxic_markers()
        
        
    def _load_yaml_db(self, filename: str) -> dict:
        path = ROOT / "data" / "species" / filename
        if not path.exists():
            return {}
        with open(path, "r") as f:
            data = yaml.safe_load(f)
        return {item["name"]: item for item in data.get("compounds", [])}

    def _load_toxic_markers(self) -> dict:
        return self._load_yaml_db("toxic_markers.yml")
        
    def _load_desirable(self) -> dict:
        return self._load_yaml_db("desirable_targets.yml")
        
    def _load_off_flavours(self) -> dict:
        return self._load_yaml_db("off_flavour_targets.yml")
        
    def _load_results(self) -> dict:
        if self.results_path is None or not self.results_path.exists():
            print(f"ERROR: Screening results not found at {self.results_path}")
            print("Please run `python scripts/run_curated_screening.py` first.")
            sys.exit(1)
            
        with open(self.results_path, "r") as f:
            data = json.load(f)
            
        # Map pathway name to span
        return {item["pathway"]: item["energetic_span_kcal"] for item in data}

    def _get_pathway_requirements(self, pathway_name: str) -> Set[str]:
        """Extract the exogenous reactants required for a pathway."""
        steps = PATHWAYS[pathway_name]
        produced_intermediates = set()
        required_exogenous = set()
        
        for step in steps:
            for reactant in step.reactants:
                if reactant.label not in produced_intermediates:
                    required_exogenous.add(reactant.label)
            for product in step.products:
                produced_intermediates.add(product.label)
                
        return required_exogenous

    def predict_from_steps(self, 
                           steps: List[Any], 
                           barriers_dict: Dict[str, float], 
                           initial_concentrations: Dict[str, float], 
                           temperature_kelvin: float = 423.15, 
                           time_minutes: Optional[float] = None,
                           protein_type: str = "free",
                           denaturation_state: float = 0.5,
                           temp_ramp_csv: Optional[str] = None):
        """
        Dynamically predict active pathways given a list of generated ElementarySteps
        and their computed barriers from xTB or Hammond fallback.
        
        If temp_ramp_csv provided (SOTA), integrates propensity over the ramp.
        """
        p_type = ProteinType(protein_type)
        ramp_data = _load_ramp(temp_ramp_csv) if temp_ramp_csv else None
        
        # Phase 7: Apply Matrix Accessibility Corrections to precursors
        # We need to map labels to concentrations for apply_matrix_correction
        # Note: apply_matrix_correction handles the scaling.
        _, corrected_initial = apply_matrix_correction(
            predicted_concentrations={}, 
            reactive_amino_acids={k: v for k, v in initial_concentrations.items()},
            protein_type=p_type,
            denaturation_state=denaturation_state
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

        changed = True
        iterations = 0
        max_iterations = len(steps) + 1  # Longest possible path
        
        while changed and iterations < max_iterations:
            changed = False
            iterations += 1
            
            for step in steps:
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
                reachable = True
                dominant_reactant = None
                dominant_reactant_span = -1.0
                
                for r in r_canons:
                    if r not in tracking:
                        reachable = False
                        break
                    r_span, r_conc, r_depth, r_weight, r_unc = tracking[r]
                    max_r_dist = max(max_r_dist, r_span)
                    max_r_unc = max(max_r_unc, r_unc)
                    min_r_conc = min(min_r_conc, r_conc)
                    max_r_depth = max(max_r_depth, r_depth)
                    if r_span >= dominant_reactant_span:
                        dominant_reactant = r
                        dominant_reactant_span = r_span
                    
                if not reachable:
                    continue
                    
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
                
                path_conc = min_r_conc
                path_depth = max_r_depth + 1
                
                # Phase G: Concentration-Aware Weighting
                # Flux = (product of reactant concs) * exp(-barrier/RT)
                # But for the cumulative pathway, we use the bottleneck span
                import math
                
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
                    RT = 0.001987 * temperature_kelvin
                    # Path weight (Flux approximation)
                    path_weight = reactant_flux_pool * co_reactant_factor * math.exp(-path_span / RT)
                    
                    # Phase Q.1: Temporal FAST Mode (Fallback)
                    if time_minutes is not None:
                        # Characteristic time approx (seconds)
                        # Using A ~ 1e11 (from new Arrhenius data average)
                        tau_sec = math.exp(path_span / RT) / 1e11
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
                    "reactants": [species_name_lookup.get(r, r) for r in r_canons],
                    "products": [species_name_lookup.get(p, p) for p in p_canons],
                }
                candidate_path = base_path + [step_trace]
                
                # Relaxation: we primarily want the lowest span path. 
                for p in p_canons:
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

        # Project ranked FAST outputs onto a bounded volatile ppb budget.
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
        )

        # Phase 7: Apply Matrix Retention Corrections to outputs
        corrected_volatiles, _ = apply_matrix_correction(
            predicted_concentrations=raw_concentrations,
            reactive_amino_acids={},
            protein_type=p_type,
            denaturation_state=denaturation_state
        )
        
        # Add names to the dictionary for downstream sensory and benchmark matching
        final_volatiles = {}
        for p_canon, conc in corrected_volatiles.items():
            final_volatiles[p_canon] = conc
            species_name = species_name_lookup.get(p_canon)
            if species_name:
                final_volatiles[species_name] = conc
            t_info = target_lookup.get(p_canon)
            if t_info:
                final_volatiles[t_info["name"]] = conc

        # Identify which targets were produced
        active_pathways = []
        for p_canon, (span, conc, depth, weight, unc) in tracking.items():
            t_info = target_lookup.get(p_canon)
            if t_info and span < float('inf') and span >= 0.0:
                
                # We mock a Species object for the old table formatter
                class MockTarget:
                    def __init__(self, label):
                        self.label = label
                
                p_dict = {
                    "name": t_info["name"],
                    "span": span,
                    "concentration": final_volatiles.get(p_canon, 0.0), # Use corrected
                    "weighted_flux": final_volatiles.get(p_canon, 0.0), # Use corrected
                    "span_uncertainty": tracking[p_canon][4],
                    "depth": depth,
                    "target": MockTarget(t_info["name"]),
                    "type": t_info["type"],
                    "penalty": "LOW",
                    "toxicity": None,
                    "sensory": t_info["data"].get("sensory_desc", "-"),
                    "threshold": t_info["data"].get("odour_threshold_ug_per_kg", None)
                }
                
                if t_info["type"] == "toxic":
                    p_dict["toxicity"] = {
                        "name": t_info["name"],
                        "risk": t_info["data"].get("health_risk", "Unknown"),
                        "priority": t_info["data"].get("priority", "high").upper()
                    }
                
                active_pathways.append(p_dict)
                
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
            "lysine_budget_dha": lysine_budget
        }

        return {
            "targets": active_pathways,
            "metrics": metrics,
            "predicted_ppb": final_volatiles,
            "debug_paths": best_paths,
            "species_names": species_name_lookup,
        }

    def predict(self, pool: List[str]):
        """Predict the outcome for a given pool of precursors (static curated)."""
        available_species = set(pool)
        
        # Ubiquitous molecules present in Maillard reaction environments:
        # water is the solvent, H2 and NH3 are common by-products that accumulate,
        # CO2 is released in decarboxylation steps. These should not block pathway
        # activation since they are always available in any food-chemistry system.
        IMPLICIT_AMBIENT = {"water", "hydrogen", "ammonia", "CO2"}
        available_species |= IMPLICIT_AMBIENT
        
        active_pathways = []
        
        # Iteratively activate pathways (since one pathway can feed another)
        added_new = True
        while added_new:
            added_new = False
            for p_name, steps in PATHWAYS.items():
                if p_name in [p["name"] for p in active_pathways]:
                    continue
                    
                reqs = self._get_pathway_requirements(p_name)
                
                if reqs.issubset(available_species):
                    # Activate!
                    span = self.screening_data.get(p_name, float('inf'))
                    meta = PATHWAY_METADATA.get(p_name, {})
                    
                    active_pathways.append({
                        "name": p_name,
                        "span": span,
                        "target": meta.get("target", None),
                        "type": meta.get("type", "unknown")
                    })
                    
                    # Add its products to the available pool (so Strecker can fire)
                    for step in steps:
                        for prod in step.products:
                            if prod.label not in available_species:
                                available_species.add(prod.label)
                                added_new = True
                                
        # Sort active pathways by kinetic probability (energetic span)
        active_pathways.sort(key=lambda x: x["span"])
        
        # Calculate penalties and extract toxicity
        for p in active_pathways:
            p["penalty"] = "LOW"
            p["toxicity"] = None
            
            # Toxicity check
            tox_flag = PATHWAY_METADATA.get(p["name"], {}).get("toxicity_flag")
            if tox_flag:
                marker = self.toxic_markers.get(tox_flag, {})
                p["toxicity"] = {
                    "name": tox_flag,
                    "risk": marker.get("health_risk", "Unknown risk"),
                    "priority": marker.get("priority", "medium").upper()
                }
            
            # Penalty check for desirable pathways
            if p["type"] == "desirable":
                desirable_span = p["span"]
                desirable_consumes = set(PATHWAY_METADATA.get(p["name"], {}).get("consumes", []))
                
                penalty_score = 0.0
                for comp_p in active_pathways:
                    if comp_p["type"] in ["competing", "masking"]:
                        comp_consumes = set(PATHWAY_METADATA.get(comp_p["name"], {}).get("consumes", []))
                        shared = desirable_consumes.intersection(comp_consumes)
                        for _ in shared:
                            # Faster competing pathway (lower span) = higher penalty
                            penalty_score += desirable_span / max(0.1, comp_p["span"])
                            
                if penalty_score < 0.5:
                    p["penalty"] = "LOW"
                elif penalty_score <= 1.0:
                    p["penalty"] = "MEDIUM"
                else:
                    p["penalty"] = "HIGH"
        
        return active_pathways


def main():
    print("======================================================")
    print("      Maillard Formulation Recommender Engine")
    print("======================================================\n")
    
    results_path = ROOT / "results" / "curated_screening_results.json"
    recommender = Recommender(results_path)
    
    for system in SYSTEMS:
        print("-" * 60)
        print(f"System: {system.name}")
        print(f"Input:  {', '.join(system.precursors)}")
        print(f"Notes:  {system.notes}")
        print("-" * 60)
        
        active = recommender.predict(system.precursors)
        
        if not active:
            print("  [!] No pathways active. The precursors do not react under these rules.")
            continue
            
        print("  Active Pathways:")
        
        # Table Header
        print("    ┌" + "─"*24 + "┬" + "─"*18 + "┬" + "─"*15 + "┬" + "─"*15 + "┬" + "─"*22 + "┐")
        print("    │ PREDICTED COMPOUND     │ PATHWAY TYPE     │ BARRIER (ΔE‡) │ PENALTY RISK  │ TOXICITY ALERT       │")
        print("    ├" + "─"*24 + "┼" + "─"*18 + "┼" + "─"*15 + "┼" + "─"*15 + "┼" + "─"*22 + "┤")
        
        for p in active:
            target_str = p['target'].label if p['target'] else "Unknown"
            
            # Formatting tags based on type
            tag = ""
            if p['type'] == 'desirable':
                tag = "[✅ AROMA]"
            elif p['type'] == 'competing':
                tag = "[⚠️ COMPETING]"
            elif p['type'] == 'masking':
                tag = "[🛡️ MASKING]"
            
            barrier_str = f"{p['span']:.1f} kcal"
            penalty_str = p['penalty']
            
            tox_str = "-"
            if p.get('toxicity'):
                meta = p['toxicity']
                tox_str = f"[{meta['priority']}] {meta['name']}"
                
            # Note: emojis can throw off terminal alignment slightly due to double-width rendering,
            # but we use a loose truncation to handle it fine in most modern terminals.
            # Emojis take up 1 char in len() but 2 visual slots.
            # We will pad manually accounting for the emojis in the tags.
            # The exact visual alignment might be slightly off by 1 space per emoji.
            
            # Truncate and pad
            col1 = _trunc(target_str, 22)
            col2 = _trunc(tag, 16)
            col3 = _trunc(barrier_str, 13)
            col4 = _trunc(penalty_str, 13)
            col5 = _trunc(tox_str, 20)
            
            print(f"    │ {col1} │ {col2} │ {col3} │ {col4} │ {col5} │")
            
        print("    └" + "─"*24 + "┴" + "─"*18 + "┴" + "─"*15 + "┴" + "─"*15 + "┴" + "─"*22 + "┘")

    print("\n" + "═"*85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - Confidence values (xTB ΔE‡ barriers) reflect relative kinetic rankings only.")
    print("    - Absolute yield predictions require higher-level Tier 2 DFT (Skala) and Cantera")
    print("      microkinetic modeling to account for temporal concentration profiles.")
    print("═"*85 + "\n")


if __name__ == "__main__":
    main()
