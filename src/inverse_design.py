import yaml
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass

from src.smirks_engine import SmirksEngine, ReactionConditions
from src.recommend import Recommender
from src.precursor_resolver import resolve_many
from src.barrier_constants import get_barrier, HEME_CATALYST_REDUCTION, HEME_CATALYST_FAMILIES

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
GRID_FILE = ROOT / "data" / "formulation_grid.yml"
TAGS_FILE = ROOT / "data" / "species" / "sensory_tags.yml"

@dataclass
class FormulationResult:
    name: str
    target_score: float
    off_flavour_risk: float
    lysine_budget: float
    trapping_efficiency: float
    detected_targets: List[str]
    detected_minimize: List[str]

class InverseDesigner:
    def __init__(self, target_tag: str, minimize_tag: str = "beany"):
        self.target_tag = target_tag.lower()
        self.minimize_tag = minimize_tag.lower()
        self.grid = self._load_grid()
        self.tags = self._load_tags()
        
        if self.target_tag not in self.tags:
            raise ValueError(f"Unknown target tag: '{self.target_tag}'. Available: {list(self.tags.keys())}")
        if self.minimize_tag and self.minimize_tag not in self.tags:
             raise ValueError(f"Unknown minimize tag: '{self.minimize_tag}'. Available: {list(self.tags.keys())}")

    def _load_grid(self) -> List[Dict]:
        if not GRID_FILE.exists():
            return []
        with open(GRID_FILE, "r") as f:
            data = yaml.safe_load(f)
            return data.get("formulations", [])
            
    def _load_tags(self) -> Dict[str, List[str]]:
        if not TAGS_FILE.exists():
            return {}
        with open(TAGS_FILE, "r") as f:
            data = yaml.safe_load(f)
            return data.get("tags", {})
            
    def _score_targets(self, targets_found: List[dict], target_list: List[str], conditions: ReactionConditions) -> Tuple[float, List[str]]:
        """
        Score a list of found targets against a desired list of compound names.
        Score uses Boltzmann weighting: [conc] * exp(-barrier / RT) * (0.8^depth)
        """
        import math
        score = 0.0
        detected = []
        
        # T in Kelvin, R in kcal/(mol*K)
        temp_k = conditions.temperature_celsius + 273.15
        RT = 0.001987 * temp_k
        
        for t in targets_found:
            name = t["name"]
            if name in target_list:
                barrier = t.get("span", 99.0)
                conc = t.get("concentration", 1.0)
                depth = t.get("depth", 1)
                
                # Baseline 0 score for unachievable pathways
                if barrier >= 99.0:
                    continue
                    
                boltzmann_weight = math.exp(-barrier / RT)
                depth_penalty = 0.8 ** depth # 20% penalty per reaction step
                
                # Scale up by 1e6 for readability (otherwise scores are tiny)
                impact = conc * boltzmann_weight * depth_penalty * 1e6
                score += impact
                detected.append(name)
        return score, detected

    def evaluate_all(self, global_conditions: ReactionConditions) -> List[FormulationResult]:
        """
        Run the generative pipeline for every formulation in the grid.
        Returns a ranked list of results.
        """
        results = []
        target_compounds = self.tags.get(self.target_tag, [])
        minimize_compounds = self.tags.get(self.minimize_tag, []) if self.minimize_tag else []

        for form in self.grid:
            name = form.get("name", "Unknown")
            sugars = form.get("sugars", [])
            amino_acids = form.get("amino_acids", [])
            additives = form.get("additives", [])
            lipids = form.get("lipids", [])
            catalyst = form.get("catalyst", None)
            
            # Combine all names
            all_names = sugars + amino_acids + additives + lipids
            try:
                precursors = resolve_many(all_names)
            except ValueError as e:
                print(f"Skipping '{name}': {e}")
                continue
                
            # Create a localized conditions object (e.g. to apply catalyst override if specified)
            cond = ReactionConditions(
                pH=form.get("ph", global_conditions.pH),
                temperature_celsius=form.get("temp", global_conditions.temperature_celsius),
                water_activity=form.get("aw", global_conditions.water_activity)
            )
            
            # FAST mode heuristic barrier overrides
            heuristic_barriers = {}
            # Apply heme catalyst heuristic to Strecker/Pyrazines if formulation has it
            apply_heme = (catalyst == "heme" or global_conditions.pH > 7)  # Note: just relying on simple logic for demo
            
            engine = SmirksEngine(cond)
            steps = engine.enumerate(precursors, max_generations=4)
            
            # Calculate barriers from centralised constants
            for step in steps:
                bar = get_barrier(step.reaction_family)
                    
                ph_mult = cond.get_ph_multiplier(step.reaction_family or "")
                if ph_mult > 1.0:
                    bar /= ph_mult
                elif ph_mult < 1.0:
                    bar += (1.0 - ph_mult) * 10.0
                    
                # Heme heuristic application
                if apply_heme and step.reaction_family in HEME_CATALYST_FAMILIES:
                    bar -= HEME_CATALYST_REDUCTION
                    
                rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                heuristic_barriers[rxn_key] = max(0.0, bar)
                
            # Build canonical concentrations map
            from src.recommend import _canon
            initial_concentrations = {}
            ratios = form.get("molar_ratios", {})
            for p in precursors:
                # Default ratio is 1.0 if not specified
                qty = 1.0
                for k, v in ratios.items():
                    if k.lower() in p.label.lower() or p.label.lower() in k.lower():
                        qty = float(v)
                        break
                initial_concentrations[_canon(p.smiles)] = qty
                
            recommender = Recommender()
            rec_result = recommender.predict_from_steps(steps, heuristic_barriers, initial_concentrations)
            
            # Score against tags
            t_score, t_detected = self._score_targets(rec_result["targets"], target_compounds, cond)
            m_score, m_detected = self._score_targets(rec_result["targets"], minimize_compounds, cond)
            
            trap_dict = rec_result["metrics"].get("trapping_efficiency", {})
            trap_avg = sum(trap_dict.values()) / len(trap_dict) if trap_dict else 0.0

            results.append(FormulationResult(
                name=name,
                target_score=t_score,
                off_flavour_risk=m_score,
                lysine_budget=rec_result["metrics"].get("lysine_budget_dha", 0.0),
                trapping_efficiency=trap_avg,
                detected_targets=t_detected,
                detected_minimize=m_detected
            ))
            
        # Rank by target_score DESC, off_flavour_risk ASC
        results.sort(key=lambda x: (x.target_score, -x.off_flavour_risk), reverse=True)
        return results
