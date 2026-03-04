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
import networkx as nx
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Set, Optional

try:
    from rdkit import Chem
except ImportError:
    pass

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

from data.reactions.curated_pathways import PATHWAYS, PATHWAY_METADATA


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

    def predict_from_steps(self, steps: List[any], barriers_dict: Dict[str, float], initial_pool_smiles: List[str]):
        """
        Dynamically predict active pathways given a list of generated ElementarySteps
        and their computed barriers from xTB or Hammond fallback.
        Works without hardcoded pathways!
        """
        desirable = self._load_desirable()
        off_flavours = self._load_off_flavours()
        toxic = self._load_toxic_markers()
        
        # Build canonical SMILES lookup for targets
        def _canon(smi):
            try:
                can = set(Chem.MolToSmiles(Chem.MolFromSmiles(smi)).split("."))
                # just return the largest fragment if disconnected
                return max(can, key=len)
            except:
                return smi

        target_lookup = {}
        for db, t_type in [(desirable, "desirable"), (off_flavours, "competing"), (toxic, "toxic")]:
            for name, data in db.items():
                if data.get("smiles"):
                    can = _canon(data["smiles"])
                    target_lookup[can] = {"name": name, "type": t_type, "data": data}

        # Build hypergraph relaxation to find min-max barrier
        distances = {}
        for s in initial_pool_smiles:
            distances[_canon(s)] = 0.0

        changed = True
        iterations = 0
        max_iterations = len(steps) + 1  # Longest possible path
        
        while changed and iterations < max_iterations:
            changed = False
            iterations += 1
            
            for step in steps:
                step_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                barrier = barriers_dict.get(step_key, 99.0)
                
                r_canons = [_canon(r.smiles) for r in step.reactants]
                p_canons = [_canon(p.smiles) for p in step.products]
                
                # The distance to fire this step is the MAX distance of all its reactants
                max_r_dist = 0.0
                reachable = True
                for r in r_canons:
                    if r not in distances:
                        reachable = False
                        break
                    max_r_dist = max(max_r_dist, distances[r])
                    
                if not reachable:
                    continue
                    
                # Path barrier to products via this step
                path_barrier = max(max_r_dist, barrier)
                
                # Relaxation
                for p in p_canons:
                    if p not in distances or path_barrier < distances[p]:
                        distances[p] = path_barrier
                        changed = True

        # Identify which targets were produced
        active_pathways = []
        for p_canon, span in distances.items():
            if p_canon in target_lookup and span < float('inf') and span >= 0.0:
                t_info = target_lookup[p_canon]
                
                # We mock a Species object for the old table formatter
                class MockTarget:
                    def __init__(self, label):
                        self.label = label
                
                p_dict = {
                    "name": t_info["name"],
                    "span": span,
                    "target": MockTarget(t_info["name"]),
                    "type": t_info["type"],
                    "penalty": "LOW",
                    "toxicity": None
                }
                
                if t_info["type"] == "toxic":
                    p_dict["toxicity"] = {
                        "name": t_info["name"],
                        "risk": t_info["data"].get("health_risk", "Unknown"),
                        "priority": t_info["data"].get("priority", "high").upper()
                    }
                
                active_pathways.append(p_dict)
                
        # Sort active pathways by kinetic probability (energetic span)
        active_pathways.sort(key=lambda x: x["span"])
        
        # Penalties: very rough heuristic for dynamic pathways
        # If there are competing off-flavours with lower barriers, increase penalty
        for p in active_pathways:
            if p["type"] == "desirable":
                comp_spans = [cp["span"] for cp in active_pathways if cp["type"] == "competing"]
                if not comp_spans:
                    p["penalty"] = "LOW"
                else:
                    min_comp = min(comp_spans)
                    if min_comp < p["span"]:
                        p["penalty"] = "HIGH"
                    elif min_comp < p["span"] + 5.0:
                        p["penalty"] = "MEDIUM"
                        
        return active_pathways

    def predict(self, pool: List[str]):
        """Predict the outcome for a given pool of precursors (static curated)."""
        available_species = set(pool)
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
            if p['type'] == 'desirable': tag = "[✅ AROMA]"
            elif p['type'] == 'competing': tag = "[⚠️ COMPETING]"
            elif p['type'] == 'masking': tag = "[🛡️ MASKING]"
            
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
