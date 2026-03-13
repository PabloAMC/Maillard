import yaml
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, field

from src.smirks_engine import SmirksEngine, ReactionConditions, Species  # noqa: E402
from src.recommend import Recommender  # noqa: E402
from src.precursor_resolver import resolve_many  # noqa: E402
from src.barrier_constants import HEME_CATALYST_REDUCTION, HEME_CATALYST_FAMILIES  # noqa: E402
from src.results_db import ResultsDB  # noqa: E402
from src.sensory import SensoryPredictor  # noqa: E402
from src.safety import evaluate_formulation_safety  # noqa: E402
from src.lipid_oxidation import predict_lop_generation  # noqa: E402

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
GRID_FILE = ROOT / "data" / "formulation_grid.yml"
TAGS_FILE = ROOT / "data" / "species" / "sensory_tags.yml"

@dataclass
class FormulationResult:
    name: str
    target_score: float
    off_flavour_risk: float
    lysine_budget: float = 0.0
    trapping_efficiency: float = 0.0
    detected_targets: List[str] = field(default_factory=list)
    detected_minimize: List[str] = field(default_factory=list)
    radar: Dict[str, tuple] = field(default_factory=dict)
    safety_score: float = 0.0
    flagged_toxics: List[str] = field(default_factory=list)
    texture_risk: float = 0.0
    predicted_ppb: Dict[str, float] = field(default_factory=dict)
    avg_uncertainty: float = 5.0


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
        
        # Initialize Database connection for unified barrier lookup
        self.db = ResultsDB()
        self.sensory = SensoryPredictor()

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
        score = 0.0
        detected = []
        
        # T in Kelvin, R in kcal/(mol*K)
        
        for t in targets_found:
            name = t["name"]
            if name in target_list:
                barrier = t.get("span", 99.0)
                flux = t.get("weighted_flux", 0.0)
                depth = t.get("depth", 1)
                
                # Baseline 0 score for unachievable pathways
                if barrier >= 99.0:
                    continue
                    
                depth_penalty = 0.8 ** depth # 20% penalty per reaction step
                
                # Scale up by 1e6 for readability (otherwise scores are tiny)
                impact = flux * depth_penalty * 1e6
                score += impact
                detected.append(name)
        return score, detected

    def _score_texture_risk(self, precursors: List[Species], r_sugars: List[str]) -> float:
        """
        Heuristic for DHA-Maillard competition.
        High risk if [Sugar] / ([Serine] + [Cysteine]) is low.
        Returns a risk score 0-100.
        """
        sugar_total = 0.0
        dha_precursor_total = 0.0
        
        # Check species labels or smiles
        for p in precursors:
            label = p.label.lower()
            # Crude molar estimation from tags/labels
            if any(s in label for s in ["serine", "cysteine"]):
                dha_precursor_total += 1.0 # Base stoichiometry
            if any(s in label for s in ["glucose", "ribose", "xylose", "fructose", "maltose", "sucrose"]):
                sugar_total += 1.0
                
        # If no Ser/Cys, no DHA risk
        if dha_precursor_total == 0:
            return 0.0
            
        ratio = sugar_total / dha_precursor_total
        # If ratio < 0.5, risk is high (DHA dominates)
        # Risk = 100 * (1 - ratio/0.5) clamped
        risk = max(0.0, 100.0 * (1.0 - ratio / 0.5))
        return risk


    def evaluate_single(self, formulation: Dict, global_conditions: ReactionConditions) -> FormulationResult:
        """
        R.8.2: Evaluate a single formulation without touching self.grid.
        """
        old_grid = self.grid
        try:
            self.grid = [formulation]
            results = self.evaluate_all(global_conditions)
            if not results:
                raise ValueError("Evaluation failed for formulation")
            return results[0]
        finally:
            self.grid = old_grid

    def evaluate_all(self, global_conditions: ReactionConditions, grid_override: Optional[List[Dict]] = None) -> List[FormulationResult]:
        """
        Run the generative pipeline for every formulation in the grid.
        Returns a ranked list of results.
        """
        results = []
        target_compounds = self.tags.get(self.target_tag, [])
        minimize_compounds = self.tags.get(self.minimize_tag, []) if self.minimize_tag else []

        eval_grid = grid_override if grid_override is not None else self.grid
        print(f"DEBUG: evaluate_all starting for {len(eval_grid)} formulations. First 2 names: {[f.get('name') for f in eval_grid[:2]]}")
        for form in eval_grid:
            name = form.get("name", "Unknown")
            protein_type = form.get("protein_type", global_conditions.protein_type)
            denaturation_state = form.get("denaturation_state", 0.5)
            
            sugars = form.get("sugars", [])
            amino_acids = form.get("amino_acids", [])
            additives = form.get("additives", [])
            lipids = form.get("lipids", [])
            catalyst = form.get("catalyst", None)
            interventions = form.get("interventions", [])
            
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
                water_activity=form.get("aw", global_conditions.water_activity),
                fat_fraction=global_conditions.fat_fraction,
                protein_fraction=global_conditions.protein_fraction,
                protein_type=protein_type
            )
            
            # FAST mode heuristic barrier overrides
            heuristic_barriers = {}
            # Apply heme catalyst heuristic to Strecker/Pyrazines if formulation has it
            apply_heme = (catalyst == "heme" or global_conditions.pH > 7)  # Note: just relying on simple logic for demo
            
            # Phase 20: Pre-calculate intervention modifiers from library
            modifiers = {}
            if interventions:
                import yaml
                LIBRARY_PATH = ROOT / "data" / "interventions.yml"
                if LIBRARY_PATH.exists():
                    with open(LIBRARY_PATH, "r") as f:
                        lib_data = yaml.safe_load(f)
                        intervention_lib = {a["name"]: a for a in lib_data.get("interventions", [])}
                    
                    print(f"DEBUG: {name} raw interventions: {interventions}")
                    for inter in interventions:
                        agent_name = inter.get("name")
                        dose = inter.get("dose", 1.0)
                        print(f"DEBUG: loading intervention {agent_name} with dose {dose}")
                        agent_data = intervention_lib.get(agent_name)
                        if agent_data:
                            for mech in agent_data.get("mechanisms", []):
                                target = mech.get("target_family", "")
                                delta = mech.get("delta_barrier_per_unit", 0.0)
                                modifiers[target] = delta * dose
                print(f"DEBUG_INV: {name} final modifiers: {modifiers}")
            
            engine = SmirksEngine(cond)
            steps = engine.enumerate(precursors, max_generations=4)
            
            # Calculate barriers from unified source (DB first, then heuristic)
            for step in steps:
                reactants = [s.smiles for s in step.reactants]
                products = [s.smiles for s in step.products]
                
                bar, source, unc = self.db.get_best_barrier(reactants, products, step.reaction_family or "unknown")
                    
                # Phase 7: Use new Arrhenius microkinetics with ionization corrections
                # k = get_rate_constant returns a rate constant. 
                # We need to map this back to an effective barrier for the span logic:
                # k = A * exp(-E_eff / RT) => E_eff = -RT * ln(k / A)
                # But even better: we can just use the rate constant directly in future iterations.
                # For now, we adjust the barrier.
                
                k = cond.get_rate_constant(step.reaction_family or "unknown", ea_override_kcal=bar)
                
                # Effective barrier including pH and aw effects:
                # Ea_eff = -RT * ln(k / A)
                import math
                RT = 0.001987 * cond.temperature_kelvin
                A = 1e11 # Must match conditions.py
                if k > 0:
                    bar_eff = -RT * math.log(k / A)
                else:
                    bar_eff = 99.0
                    
                if apply_heme and step.reaction_family in HEME_CATALYST_FAMILIES:
                    bar_eff -= HEME_CATALYST_REDUCTION
                
                # Phase 20: Apply intervention modifiers
                family = step.reaction_family or ""
                for mod_fam, delta in modifiers.items():
                    if mod_fam.lower() in family.lower():
                        bar_eff += delta
                    
                rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
                heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)
                
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
                
            # Create a name-to-concentration map for the safety module
            name_ratios = {p.label: initial_concentrations[_canon(p.smiles)] for p in precursors}
                
            # PHASE L: Lipid Oxidation Crosstalk (Fix 7)
            lipids_input = {k: v for k, v in ratios.items() if any(l in k.lower() for l in ["linoleic", "linolenic", "oil", "fat"])}
            generated_lops = predict_lop_generation(
                lipids_input, 
                cond.temperature_celsius, 
                form.get("time_minutes", 60.0),
                cond.water_activity
            )
            
            # Merge LOPs into initial_concentrations
            for lop_smi, lop_conc in generated_lops.items():
                initial_concentrations[_canon(lop_smi)] = initial_concentrations.get(_canon(lop_smi), 0.0) + lop_conc
                
            recommender = Recommender()
            rec_result = recommender.predict_from_steps(
                steps, 
                heuristic_barriers, 
                initial_concentrations, 
                temperature_kelvin=cond.temperature_kelvin,
                protein_type=protein_type,
                denaturation_state=denaturation_state
            )
            
            # Score against tags
            t_score_legacy, t_detected = self._score_targets(rec_result["targets"], target_compounds, cond)
            m_score_legacy, m_detected = self._score_targets(rec_result["targets"], minimize_compounds, cond)
            
            # PHASE F: Safety evaluation
            # Replace old heuristic with new quantitative model
            safety_val, flagged = evaluate_formulation_safety(
                name_ratios, 
                cond.temperature_celsius, 
                form.get("time_minutes", 60.0), 
                cond.pH,
                modifiers=modifiers
            )
            s_penalty = safety_val * 2.0 # Weight for optimization
            
            trap_dict = rec_result["metrics"].get("trapping_efficiency", {})
            trap_avg = sum(trap_dict.values()) / len(trap_dict) if trap_dict else 0.0

            # Calculate perceived sensory profile (Phase C/D)
            # Kinetic-aware weighting: Conc_eff = Conc_limiting * exp(-(PathSpan - 20) / RT)
            # We use 20 kcal/mol as a baseline to make results visible in ppm units.
            RT = 0.001987 * cond.temperature_kelvin
            # Populate conc_map for legacy compatibility (though radar scores now use ppb)
            conc_map = rec_result["predicted_ppb"]

            # PHASE G: Sensory Scoring
            radar_scores = self.sensory.get_radar_data(
                rec_result["predicted_ppb"], 
                protein_type=protein_type,
                temp_c=cond.temperature_celsius,
                fat_fraction=cond.fat_fraction,
                protein_fraction=cond.protein_fraction
            )
            
            # Calculate average uncertainty for the optimizer
            valid_spans = [t["span_uncertainty"] for t in rec_result["targets"] if t["span_uncertainty"] > 0]
            avg_unc = sum(valid_spans) / len(valid_spans) if valid_spans else 5.0
            
            # Use radar score for the target category as the official t_score
            t_score = radar_scores.get(self.target_tag, (0.0, 0))[0]
            # Use radar score for the minimize category as the official m_score
            m_score = radar_scores.get(self.minimize_tag, (0.0, 0))[0] if self.minimize_tag else 0.0

            results.append(FormulationResult(
                name=name,
                target_score=t_score,
                off_flavour_risk=m_score,
                lysine_budget=rec_result["metrics"].get("lysine_budget_dha", 0.0),
                trapping_efficiency=trap_avg,
                detected_targets=t_detected,
                detected_minimize=m_detected,
                radar=radar_scores,
                safety_score=s_penalty,
                flagged_toxics=flagged,
                texture_risk=self._score_texture_risk(precursors, sugars),
                predicted_ppb=conc_map,
                avg_uncertainty=avg_unc
            ))

            
        # Rank by (target_score - risk_aversion * safety_score - texture_penalty) DESC, off_flavour_risk ASC
        risk_aversion = 1.0 # Default
        texture_aversion = 0.01 # 100 risk = 1.0 score penalty
        results.sort(
            key=lambda x: (x.target_score - risk_aversion * x.safety_score - texture_aversion * x.texture_risk, -x.off_flavour_risk), 
            reverse=True
        )
        return results
