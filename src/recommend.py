#!/usr/bin/env python3
"""
src/recommend.py — Maillard Reaction Pathway Recommender

Tier 4 Pipeline Integration:
Loads the Tier 1 xTB screening results and matches them against user-defined
or canonical precursors to recommend actionable formulation adjustments.

P2.2 refactor: projection helpers and output-projection logic now live in
``src/projection.py``.  YAML target loading now uses the cached registry in
``src/target_registry.py``.  All symbols are re-exported from this module so
that existing imports continue to work unchanged.
"""

import sys
import json
import math
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Set, Optional, Any, Tuple

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.logger import get_logger
logger = get_logger(__name__)

# ── Re-export everything from projection.py (backward compat + direct access) ─
from src.projection import (
    # Original projection dataclasses / budget helpers
    ProjectionBudget, ProjectionStrategy, DEFAULT_PROJECTION_STRATEGY,
    _thermal_severity, _projection_temperature_factor, _projection_time_factor,
    _estimate_projection_budget, _temporal_accessibility, _relative_precursor_load_factor,
    _projection_strategy_metadata,
    # P2.2 extracted symbols – keep re-exported so downstream imports work
    _canon,
    _normalize_chemical_name,
    _HENRY_CONSTANTS_PATH, _NON_OBSERVABLE_KAW_THRESHOLD, _LOW_HEADSPACE_REFERENCE_KAW,
    _HENRY_LOOKUP, _HEADSPACE_MODEL,
    _BUDGET_EXCLUDED_CANONICAL,
    _carbon_count, _has_reactive_nonvolatile_functionality, _mw_from_smiles,
    _henry_entry_for_species,
    _is_observable_target_species,
    _headspace_observability_metadata, _headspace_observability_factor,
    _MELANOIDIN_TRAPPING_PROFILES, _resolve_melanoidin_trapping_factor,
    _resolve_output_matrix_context,
    _is_budget_relevant_species, _is_ppb_output_species,
    _select_accumulating_projection_species,
    _project_weighted_flux_to_ppb,
    _apply_output_projection,
)

# ── Cached YAML target registry ───────────────────────────────────────────────
from src.target_registry import get_toxic_markers, get_desirable_targets, get_off_flavour_targets

from data.reactions.curated_pathways import PATHWAYS, PATHWAY_METADATA
from src.barrier_constants import arrhenius_rate_constant, get_reference_pre_exponential
from src.conditions import ReactionConditions
from src.matrix_correction import (
    ProteinType,
    apply_matrix_correction,
)
from src.ode_kinetics import simulate_kinetic_trace
from src.pathway_extractor import Species


def _kinetic_observable_for_target(
    canon: str,
    kinetic_trace: Any,
) -> Tuple[float, float, float, str]:
    final_concentration = float(kinetic_trace.final_concentrations.get(canon, 0.0))
    peak_concentration = float(kinetic_trace.peak_concentrations.get(canon, final_concentration))
    integrated_exposure = float(kinetic_trace.integrated_concentrations.get(canon, final_concentration))
    return final_concentration, peak_concentration, integrated_exposure, "final_concentration"


def _raw_concentrations_from_kinetic_trace(
    kinetic_trace: Any,
) -> Dict[str, float]:
    return {canon: float(value) for canon, value in kinetic_trace.final_concentrations.items()}


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
    """Load a temperature ramp from CSV (time_min, temp_c)."""
    try:
        df = pd.read_csv(ramp_path)
        # Standard schema: 'time' in minutes, 'temp' in Celsius
        if 'time' not in df.columns or 'temp' not in df.columns:
            return []
        return list(zip(df['time'], df['temp'] + 273.15))
    except Exception as e:
        logger.warning(f"Failed to load ramp {ramp_path}: {e}")
        return []


# ── Functions above moved to src/projection.py (P2.2). Re-exported via import at top. ──


class Recommender:
    def __init__(self, results_path: Optional[Path] = None):
        self.results_path = results_path
        self.screening_data = self._load_results() if results_path else {}
        self.toxic_markers = get_toxic_markers()

    def _load_yaml_db(self, filename: str) -> dict:
        """Deprecated: use target_registry functions instead."""
        from src.target_registry import get_toxic_markers as _gt, get_desirable_targets as _gd, get_off_flavour_targets as _go
        mapping = {"toxic_markers.yml": _gt, "desirable_targets.yml": _gd, "off_flavour_targets.yml": _go}
        fn = mapping.get(filename)
        return fn() if fn is not None else {}

    def _load_toxic_markers(self) -> dict:
        return get_toxic_markers()

    def _load_desirable(self) -> dict:
        return get_desirable_targets()

    def _load_off_flavours(self) -> dict:
        return get_off_flavour_targets()

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
                           water_activity: Optional[float] = None,
                           protein_type: str = "free",
                           denaturation_state: float = 0.5,
                           fat_fraction: float = 0.0,
                           protein_fraction: float = 0.0,
                           process_state: Optional[str] = None,
                           temp_ramp_csv: Optional[str] = None,
                           prediction_mode: str = "projection",
                           reaction_conditions: Optional[ReactionConditions] = None,
                           kinetic_trace: Optional[Any] = None,
                           kinetic_summary: Optional[Any] = None):
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
            denaturation_state=denaturation_state
        )
        
        desirable = get_desirable_targets()
        off_flavours = get_off_flavour_targets()
        toxic = get_toxic_markers()
        
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

        kinetic_result = None
        if kinetic_trace is not None:
            kinetic_result = kinetic_trace
        elif prediction_mode == "kinetic":
            ramp_profile = _load_ramp(temp_ramp_csv) if temp_ramp_csv else None
            kinetics_conditions = reaction_conditions or ReactionConditions(
                temperature_celsius=temperature_kelvin - 273.15,
                water_activity=0.8 if water_activity is None else water_activity,
                fat_fraction=fat_fraction,
                protein_fraction=protein_fraction,
                protein_type=protein_type,
                prediction_mode="kinetic",
                temperature_profile=ramp_profile,
            )
            kinetic_result = simulate_kinetic_trace(
                steps,
                barriers_dict,
                corrected_initial,
                kinetics_conditions,
                float(time_minutes or 0.0),
            )
            kinetic_summary = kinetic_result.summary

        # tracking dict: canon_smiles -> (span, concentration, depth, weight, uncertainty)
        tracking = {}
        best_paths: Dict[str, List[Dict[str, Any]]] = {}
        # Pre-calculate exp(0) for initial precursors
        import math
        for s, conc in corrected_initial.items():
            canon = _canon(s)
            tracking[canon] = (0.0, conc, 0, conc * 1.0, 0.0)
            best_paths[canon] = []
            target_name = target_lookup.get(canon, {}).get("name")
            species_catalog[canon] = Species(species_name_lookup.get(canon) or target_name or canon, canon)

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
                for r in set(r_canons):
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
                    "reactants": [species_name_lookup.get(r, r) for r in r_canons],
                    "reactant_canons": list(r_canons),
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
        projection_budget = _estimate_projection_budget(
            corrected_initial,
            temperature_kelvin,
            time_minutes,
            strategy=projection_strategy,
        )
        if kinetic_result is None:
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
            )
        else:
            raw_concentrations = _raw_concentrations_from_kinetic_trace(kinetic_result)

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
            process_state=process_state,
            time_minutes=time_minutes,
            water_activity=water_activity,
            denaturation_state=denaturation_state,
            fat_fraction=fat_fraction,
            protein_fraction=protein_fraction,
            projection_budget=projection_budget,
            projection_strategy=projection_strategy,
            species_name_lookup=species_name_lookup,
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
            final_kinetic_conc = 0.0
            peak_kinetic_conc = 0.0
            integrated_kinetic_exposure = 0.0
            selected_observable = "projection"
            if kinetic_result is not None:
                (
                    final_kinetic_conc,
                    peak_kinetic_conc,
                    integrated_kinetic_exposure,
                    selected_observable,
                ) = _kinetic_observable_for_target(p_canon, kinetic_result)
            
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
                "peak_concentration": peak_kinetic_conc,
                "integrated_exposure": integrated_kinetic_exposure,
                "selected_observable": selected_observable,
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
                "water_activity": None if water_activity is None else float(water_activity),
                "volatile_yield_fraction": float(projection_budget.volatile_yield_fraction),
                "total_volatile_budget_molar": float(projection_budget.total_volatile_budget_molar),
                "prediction_engine": "time_resolved_microkinetic" if kinetic_result is not None else "path_span_projection",
                **_projection_strategy_metadata(projection_strategy),
            },
            "kinetic_metadata": {
                "prediction_engine": "time_resolved_microkinetic" if kinetic_result is not None else "path_span_projection",
                "time_grid_minutes": [] if kinetic_result is None else kinetic_result.time_minutes,
                "solver": None if kinetic_result is None else kinetic_result.solver,
                "successful": True if kinetic_result is None else kinetic_result.successful,
                "tracked_species": 0 if kinetic_result is None else len(kinetic_result.trajectories),
                "reaction_count": 0 if kinetic_summary is None else int(getattr(kinetic_summary, "reaction_count", 0)),
                "species_count": 0 if kinetic_summary is None else int(getattr(kinetic_summary, "species_count", 0)),
                "time_horizon_minutes": 0.0 if kinetic_summary is None else float(getattr(kinetic_summary, "time_horizon_minutes", 0.0)),
                "used_dynamic_profiles": False if kinetic_summary is None else bool(getattr(kinetic_summary, "used_dynamic_profiles", False)),
                "used_pruning": False if kinetic_summary is None else bool(getattr(kinetic_summary, "used_pruning", False)),
                "concentration_floor": None if kinetic_summary is None else float(getattr(kinetic_summary, "concentration_floor", 0.0)),
                "solver_fallback_used": False if kinetic_summary is None else bool(getattr(kinetic_summary, "solver_fallback_used", False)),
                "fallback_to_projection": False if kinetic_summary is None else bool(getattr(kinetic_summary, "fallback_to_projection", False)),
                "fallback_reason": None if kinetic_summary is None else getattr(kinetic_summary, "fallback_reason", None),
                "observable_surface": "end_state_default",
            },
            "debug_paths": best_paths,
            "species_names": species_name_lookup,
        }

    def predict(self, pool: List[str]):
        """
        [DEPRECATED] Static pathway estimation logic from Phase 1.
        Superseded by predict_from_steps() returning ElementaryStep flows.
        """
        logger.warning("Recommender.predict() is deprecated. Do not use for new implementations.")
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

