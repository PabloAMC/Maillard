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
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Set, Optional

from data.reactions.curated_pathways import PATHWAYS, PATHWAY_METADATA
try:
    from rdkit import Chem
except ImportError:
    Chem = None

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

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
        if not self.results_path.exists():
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

    def predict_from_steps(self, steps: List[any], barriers_dict: Dict[str, float], initial_concentrations: Dict[str, float], temperature_kelvin: float = 423.15, time_minutes: Optional[float] = None):
        """
        Dynamically predict active pathways given a list of generated ElementarySteps
        and their computed barriers from xTB or Hammond fallback.
        Works without hardcoded pathways!
        """
        desirable = self._load_desirable()
        off_flavours = self._load_off_flavours()
        toxic = self._load_toxic_markers()
        
        target_lookup = {}
        for db, t_type in [(desirable, "desirable"), (off_flavours, "competing"), (toxic, "toxic")]:
            for name, data in db.items():
                if data.get("smiles"):
                    can = _canon(data["smiles"])
                    target_lookup[can] = {"name": name, "type": t_type, "data": data}

        # tracking dict: canon_smiles -> (span, concentration, depth, weight, uncertainty)
        tracking = {}
        # Pre-calculate exp(0) for initial precursors
        import math
        for s, conc in initial_concentrations.items():
            tracking[_canon(s)] = (0.0, conc, 0, conc * 1.0, 0.0)

        changed = True
        iterations = 0
        max_iterations = len(steps) + 1  # Longest possible path
        
        while changed and iterations < max_iterations:
            changed = False
            iterations += 1
            
            for step in steps:
                step_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                barrier_data = barriers_dict.get(step_key, (99.0, 5.0))
                barrier, step_unc = barrier_data if isinstance(barrier_data, tuple) else (barrier_data, 5.0)
                
                r_canons = [_canon(r.smiles) for r in step.reactants]
                p_canons = [_canon(p.smiles) for p in step.products]
                
                # The distance to fire this step is the MAX distance of all its reactants
                # The limiting concentration is the MIN of all reactants
                # The depth is the MAX depth of all reactants + 1
                max_r_dist = 0.0
                max_r_unc = 0.0
                min_r_conc = float('inf')
                max_r_depth = 0
                reachable = True
                
                for r in r_canons:
                    if r not in tracking:
                        reachable = False
                        break
                    r_span, r_conc, r_depth, r_weight, r_unc = tracking[r]
                    max_r_dist = max(max_r_dist, r_span)
                    max_r_unc = max(max_r_unc, r_unc)
                    min_r_conc = min(min_r_conc, r_conc)
                    max_r_depth = max(max_r_depth, r_depth)
                    
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
                RT = 0.001987 * temperature_kelvin
                
                # Product of reactant concentrations
                reactant_conc_product = 1.0
                for r in r_canons:
                    reactant_conc_product *= tracking[r][1]
                
                # Path weight (Flux approximation)
                path_weight = reactant_conc_product * math.exp(-path_span / RT)
                
                # Phase Q.1: Temporal FAST Mode.
                # If a time is provided, penalise pathways that require many steps.
                # Approximate 1st order characteristic timescale: tau ~ exp(Ea/RT) / A
                # We assume a generic pre-exponential for this rough heuristic.
                if time_minutes is not None:
                    # Characteristic time approx (seconds)
                    # Using A ~ 1e11 (from new Arrhenius data average)
                    tau_sec = math.exp(path_span / RT) / 1e11
                    tau_min = tau_sec / 60.0
                    
                    # Number of steps increases characteristic time roughly linearly
                    total_tau = tau_min * path_depth
                    
                    # Weight decay: if total_tau >> time_minutes, weight drops exponentially
                    time_penalty = math.exp(-total_tau / time_minutes)
                    path_weight *= time_penalty
                
                # Relaxation: we primarily want the lowest span path. 
                for p in p_canons:
                    # Update if new span is lower OR if span is same but weight is higher
                    p_key = p # Assuming p is the canonical SMILES string
                    if p_key not in tracking:
                        tracking[p_key] = (float('inf'), 0.0, 0, 0.0, 0.0) # Initialize if not present
                    
                    current_span, current_conc, current_depth, current_weight, current_unc = tracking[p_key]

                    if path_span < current_span:
                        tracking[p_key] = (path_span, path_conc, path_depth, path_weight, path_unc)
                        changed = True
                    elif path_span == current_span and path_weight > current_weight:
                        tracking[p_key] = (path_span, path_conc, path_depth, path_weight, path_unc)
                        changed = True

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
                    "concentration": conc,
                    "weighted_flux": weight,
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
                            sb_weights += _weight(path_barrier)
            
            persistence_w = _weight(30.0)
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
                    weight = _weight(path_barrier)
                    
                    if step.reaction_family in ["Schiff_Base_Formation", "Lipid_Schiff_Base"]:
                        w_maillard += weight
                    elif step.reaction_family == "DHA_Crosslinking":
                        w_dha += weight
            
            if w_maillard + w_dha > 0:
                lysine_budget = 100.0 * w_dha / (w_maillard + w_dha)

        return {
            "targets": active_pathways,
            "metrics": {
                "trapping_efficiency": trapping_results,
                "lysine_budget_dha": lysine_budget
            }
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
