import yaml
import logging
import math
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field

from src.smirks_engine import SmirksEngine, ReactionConditions, Species  # noqa: E402
from src.recommend import Recommender  # noqa: E402
from src.precursor_resolver import resolve_many  # noqa: E402
from src.barrier_constants import HEME_CATALYST_REDUCTION, HEME_CATALYST_FAMILIES, effective_barrier_from_rate_constant  # noqa: E402
from src.results_db import ResultsDB  # noqa: E402
from src.sensory import SensoryPredictor  # noqa: E402
from src.safety import evaluate_formulation_safety  # noqa: E402
from src.lipid_oxidation import build_lipid_input_proxy_loads, predict_lop_generation  # noqa: E402
from src.matrix_correction import build_matrix_explainability, resolve_effective_denaturation_state  # noqa: E402
from src.literature_runtime import build_family_upstream_contract, build_flavor_axis_summary
from src.pre_processor import PreProcessor
from src.projection_metadata import ProjectionMetadataMap
from src.formulation import Formulation

# Locate data files
ROOT = Path(__file__).resolve().parents[1]
GRID_FILE = ROOT / "data" / "formulation_grid.yml"
TAGS_FILE = ROOT / "data" / "species" / "sensory_tags.yml"

_MFT_ALIASES = (
    "2-Methyl-3-furanthiol (MFT)",
    "2-methyl-3-furanthiol",
    "mft",
)
_FURFURAL_ALIASES = (
    "Furfural",
    "furfural",
)


def _max_predicted_signal(predicted_ppb: Dict[str, float], aliases: Tuple[str, ...]) -> float:
    signal = 0.0
    for alias in aliases:
        value = predicted_ppb.get(alias)
        if value is not None:
            signal = max(signal, float(value))
    return signal


def _compute_meaty_quality_constraint(predicted_ppb: Dict[str, float]) -> Tuple[float, float]:
    """
    Penalize sulfur-poor runs where furfural dominates the meaty branch.

    The threshold is intentionally broad and only activates when the observed
    MFT/furfural balance drops by more than an order of magnitude below the
    expected qualitative literature window for a meat-like outcome.
    """
    mft_ppb = _max_predicted_signal(predicted_ppb, _MFT_ALIASES)
    furfural_ppb = _max_predicted_signal(predicted_ppb, _FURFURAL_ALIASES)

    if furfural_ppb <= 0.0:
        return (1.0e6 if mft_ppb > 0.0 else 0.0, 0.0)
    if mft_ppb <= 0.0:
        return 0.0, 2.5

    ratio = mft_ppb / furfural_ppb
    target_ratio = 1.0e-2
    penalty = max(0.0, math.log10(target_ratio / max(ratio, 1.0e-9))) * 1.25
    return ratio, penalty

@dataclass
class FormulationResult:
    name: str
    target_score: float
    off_flavour_risk: float
    lysine_budget: float = 0.0
    trapping_efficiency: float = 0.0
    detected_targets: List[str] = field(default_factory=list)
    detected_minimize: List[str] = field(default_factory=list)
    radar: Dict[str, Tuple[float, float, float]] = field(default_factory=dict)
    safety_score: float = 0.0
    flagged_toxics: List[str] = field(default_factory=list)
    texture_risk: float = 0.0
    safety_metrics: Dict[str, float] = field(default_factory=dict)
    predicted_ppb: Dict[str, float] = field(default_factory=dict)
    predicted_proxy_ppb: Dict[str, float] = field(default_factory=dict)
    projection_metadata: ProjectionMetadataMap = field(default_factory=dict)
    avg_uncertainty: float = 5.0
    effective_denaturation_state: float = 0.5
    matrix_explainability: Dict[str, Any] = field(default_factory=dict)
    confidence_metadata: Dict[str, Any] = field(default_factory=dict)
    targets: List[Dict[str, Any]] = field(default_factory=list)
    bottleneck_precursor: str = "none"
    bottleneck_severity: float = 0.0
    precursor_contributions: Dict[str, float] = field(default_factory=dict)
    suppressed_compounds: List[Dict[str, Any]] = field(default_factory=list)
    mft_to_furfural_ratio: float = 0.0
    meaty_quality_penalty: float = 0.0
    strecker_balance_score: float = 0.0
    strecker_gap_penalty: float = 0.0
    pyrazine_propensity: float = 0.0
    pyrazine_burden: float = 0.0
    pyrazine_penalty: float = 0.0
    furanone_penalty: float = 0.0
    ranking_score: float = 0.0
    flavor_axis_summary: Dict[str, Any] = field(default_factory=dict)


def compute_ranking_score(
    result: FormulationResult,
    *,
    risk_aversion: float = 1.0,
    texture_aversion: float = 0.01,
) -> float:
    return (
        float(result.target_score)
        - risk_aversion * float(result.safety_score)
        - texture_aversion * float(result.texture_risk)
        - float(result.meaty_quality_penalty)
        - float(result.strecker_gap_penalty)
        - float(result.pyrazine_penalty)
        - float(result.furanone_penalty)
    )


class MaillardPipeline:
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
        self.pre_processor = PreProcessor()

    def _load_grid(self) -> List[Formulation]:
        if not GRID_FILE.exists():
            return []
        with open(GRID_FILE, "r") as f:
            data = yaml.safe_load(f)
            raw_forms = data.get("formulations", [])
            return [Formulation.from_dict(f) for f in raw_forms]
            
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


    def evaluate_single(self, formulation: Formulation, global_conditions: ReactionConditions) -> FormulationResult:
        """
        R.8.2: Evaluate a single formulation without touching self.grid.
        """
        if isinstance(formulation, dict):
            formulation = Formulation.from_dict(formulation)
        old_grid = self.grid
        try:
            self.grid = [formulation]
            results = self.evaluate_all(global_conditions)
            if not results:
                raise ValueError(f"Evaluation failed for formulation {formulation.name}")
            return results[0]
        finally:
            self.grid = old_grid

    def evaluate_all(self, global_conditions: ReactionConditions, grid_override: Optional[List[Any]] = None) -> List[FormulationResult]:
        """
        Run the generative pipeline for every formulation in the grid.
        Returns a ranked list of results.
        """
        results = []
        target_compounds = self.tags.get(self.target_tag, [])
        minimize_compounds = self.tags.get(self.minimize_tag, []) if self.minimize_tag else []

        eval_grid = grid_override if grid_override is not None else self.grid
        preview_names = [
            preview_form.get("name", "Unknown") if isinstance(preview_form, dict) else getattr(preview_form, "name", "Unknown")
            for preview_form in eval_grid[:2]
        ]
        logging.getLogger(__name__).debug(
            "evaluate_all starting for %d formulations. First 2 names: %s",
            len(eval_grid),
            preview_names,
        )
        for form in eval_grid:
            if isinstance(form, dict):
                form = Formulation.from_dict(form)
            name = form.name
            protein_type = form.protein_type or global_conditions.protein_type
            # Use global_conditions.pH as fallback when formulation does not specify pH
            base_pH = form.ph if form.ph is not None else global_conditions.pH
            denaturation_state = resolve_effective_denaturation_state(
                protein_type=protein_type,
                temperature_celsius=form.temperature,
                time_minutes=getattr(form, "time_minutes", 60.0), # Adding a fallback for non-schema fields if any
                pH=base_pH,
                explicit_denaturation_state=getattr(form, "denaturation_state", None),
            )
            
            sugars = form.sugars
            amino_acids = form.amino_acids
            additives = form.additives
            lipids = form.lipids
            catalyst = form.catalyst
            interventions = form.interventions
            raw_ratios = form.molar_ratios

            family_upstream_contract = build_family_upstream_contract(
                sugars=sugars,
                amino_acids=amino_acids,
                additives=additives,
                interventions=interventions,
                protein_type=protein_type,
                pH=base_pH,
                thiamine_availability=form.thiamine_availability,
                molar_ratios=raw_ratios,
            )
            effective_ratios = dict(family_upstream_contract.get("effective_molar_ratios", {})) or raw_ratios.copy()
            if family_upstream_contract.get("pretreatment_active"):
                effective_ratios = self.pre_processor.apply(effective_ratios, interventions)

            added_precursors = [
                precursor
                for precursor in family_upstream_contract.get("added_precursors", [])
                if precursor not in sugars + amino_acids + additives + lipids
            ]
            for precursor_name, ratio_value in (family_upstream_contract.get("added_precursor_ratios", {}) or {}).items():
                effective_ratios.setdefault(str(precursor_name), float(ratio_value))
            
            # Combine all names
            all_names = sugars + amino_acids + additives + lipids + added_precursors
            try:
                precursors = resolve_many(all_names)
            except ValueError as e:
                logging.getLogger(__name__).warning("Skipping '%s': %s", name, e)
                continue
                
            # Create a localized conditions object (e.g. to apply catalyst override if specified)
            effective_pH = family_upstream_contract.get("effective_pH", base_pH)
            cond = ReactionConditions(
                pH=effective_pH,
                temperature_celsius=form.temperature,
                water_activity=form.water_activity,
                fat_fraction=global_conditions.fat_fraction,
                protein_fraction=global_conditions.protein_fraction,
                protein_type=protein_type
            )
            
            # FAST mode heuristic barrier overrides
            heuristic_barriers = {}
            # Apply heme catalyst heuristic to Strecker/Pyrazines if formulation has it
            apply_heme = (catalyst == "heme" or cond.pH > 7)  # Note: just relying on simple logic for demo
            
            # Phase 20: Pre-calculate intervention modifiers from library
            modifiers = {}
            if interventions:
                import yaml
                LIBRARY_PATH = ROOT / "data" / "interventions.yml"
                if LIBRARY_PATH.exists():
                    with open(LIBRARY_PATH, "r") as f:
                        lib_data = yaml.safe_load(f)
                        intervention_lib = {a["name"]: a for a in lib_data.get("interventions", [])}
                    logging.getLogger(__name__).debug("%s raw interventions: %s", name, interventions)
                    for inter in interventions:
                        if isinstance(inter, dict):
                            agent_name = inter.get("name")
                            dose = inter.get("dose", 1.0)
                        else:
                            agent_name = str(inter)
                            dose = 1.0
                        logging.getLogger(__name__).debug("loading intervention %s with dose %s", agent_name, dose)
                        agent_data = intervention_lib.get(agent_name)
                        if agent_data:
                            for mech in agent_data.get("mechanisms", []):
                                target = mech.get("target_family", "")
                                delta = mech.get("delta_barrier_per_unit", 0.0)
                                modifiers[target] = delta * dose
                logging.getLogger(__name__).debug("%s final modifiers: %s", name, modifiers)
            
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
                
                bar_eff = effective_barrier_from_rate_constant(
                    k,
                    cond.temperature_kelvin,
                    step.reaction_family or "unknown",
                )
                    
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
            ratios = effective_ratios
            
            for p in precursors:
                # Default ratio is 1.0 if not specified
                qty = 1.0
                for k, v in ratios.items():
                    if k.lower() in p.label.lower() or p.label.lower() in k.lower():
                        qty = float(v)
                        break
                # Keep the historical raw precursor pool here. recommend.predict_from_steps()
                # already owns matrix accessibility scaling, and pre-scaling in the pipeline
                # double-applies corrections and breaks legacy benchmark contracts.
                initial_concentrations[_canon(p.smiles)] = qty
                
            # Create a name-to-concentration map for the safety module
            name_ratios = {p.label: initial_concentrations[_canon(p.smiles)] for p in precursors}
                
            # PHASE L: Lipid Oxidation Crosstalk (Fix 7)
            lipids_input = build_lipid_input_proxy_loads(lipids, ratios)
            generated_lops = predict_lop_generation(
                lipids_input, 
                cond.temperature_celsius, 
                form.time_minutes,
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
                time_minutes=form.time_minutes,
                water_activity=cond.water_activity,
                protein_type=protein_type,
                denaturation_state=denaturation_state,
                fat_fraction=cond.fat_fraction,
                protein_fraction=cond.protein_fraction,
            )
            
            # Score against tags
            t_score_legacy, t_detected = self._score_targets(rec_result["targets"], target_compounds, cond)
            m_score_legacy, m_detected = self._score_targets(rec_result["targets"], minimize_compounds, cond)
            
            # PHASE F: Safety evaluation
            # Replace old heuristic with new quantitative model
            safety_val, flagged, safety_result = evaluate_formulation_safety(
                name_ratios, 
                cond.temperature_celsius, 
                form.time_minutes, 
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
            
            mft_to_furfural_ratio, meaty_quality_penalty = _compute_meaty_quality_constraint(
                conc_map
            )
            flavor_axis_summary = build_flavor_axis_summary(
                projection_metadata=rec_result.get("projection_metadata", {}),
                sugars=sugars,
                amino_acids=amino_acids,
                additives=additives,
                lipids=lipids,
                interventions=interventions,
                protein_type=protein_type,
                pH=cond.pH,
                thiamine_availability=form.thiamine_availability,
                molar_ratios=effective_ratios,
                family_upstream_contract=family_upstream_contract,
            )
            strecker_gap_penalty = float(flavor_axis_summary.get("strecker_gap_penalty", 0.0))
            pyrazine_penalty = float(flavor_axis_summary.get("pyrazine_penalty", 0.0))
            furanone_penalty = float(flavor_axis_summary.get("furanone_penalty", 0.0))
            family_adjustments = flavor_axis_summary.get("family_lane_adjustments", {}) or {}
            family_target_delta = float(family_adjustments.get("target_score_delta", 0.0) or 0.0)
            family_maillard_closure_delta = float(family_adjustments.get("maillard_closure_delta", 0.0) or 0.0)
            family_offnote_delta = float(family_adjustments.get("off_flavour_risk_delta", 0.0) or 0.0)

            # Use radar score for the target category as the official t_score
            t_score = max(0.0, radar_scores.get(self.target_tag, (0.0, 0))[0] + family_target_delta + family_maillard_closure_delta)
            # Use radar score for the minimize category as the official m_score
            m_score = max(0.0, (radar_scores.get(self.minimize_tag, (0.0, 0))[0] if self.minimize_tag else 0.0) + family_offnote_delta)

            result = FormulationResult(
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
                safety_metrics={
                    "acrylamide_ppb": safety_result.acrylamide_ppb,
                    "furosine_mg_100g": safety_result.furosine_mg_100g,
                    "cml_mg_kg": safety_result.cml_mg_kg,
                    "cel_mg_kg": safety_result.cel_mg_kg,
                    "lal_mg_100g": safety_result.lal_mg_100g
                },
                predicted_ppb=conc_map,
                predicted_proxy_ppb=rec_result.get("predicted_proxy_ppb", {}),
                projection_metadata=rec_result.get("projection_metadata", {}),
                avg_uncertainty=avg_unc,
                effective_denaturation_state=denaturation_state,
                matrix_explainability=build_matrix_explainability(
                    protein_type=protein_type,
                    effective_denaturation_state=denaturation_state,
                    temperature_celsius=cond.temperature_celsius,
                    time_minutes=form.time_minutes,
                    pH=cond.pH,
                    dominant_source=(
                        "denaturation_state_arg"
                        if form.denaturation_state is not None
                        else "estimated_from_conditions"
                    ),
                ),
                targets=rec_result.get("targets", []),
                bottleneck_precursor=rec_result["metrics"].get("bottleneck", {}).get("precursor", "none"),
                bottleneck_severity=rec_result["metrics"].get("bottleneck", {}).get("severity", 0.0),
                precursor_contributions=rec_result["metrics"].get("precursor_attribution", {}),
                suppressed_compounds=rec_result["metrics"].get("suppressed_compounds", []),
                mft_to_furfural_ratio=mft_to_furfural_ratio,
                meaty_quality_penalty=meaty_quality_penalty,
                strecker_balance_score=float(flavor_axis_summary.get("strecker_balance_score", 0.0)),
                strecker_gap_penalty=strecker_gap_penalty,
                pyrazine_propensity=float(flavor_axis_summary.get("pyrazine_propensity", 0.0)),
                pyrazine_burden=float(flavor_axis_summary.get("pyrazine_burden", 0.0)),
                pyrazine_penalty=pyrazine_penalty,
                furanone_penalty=furanone_penalty,
                flavor_axis_summary=flavor_axis_summary,
            )
            results.append(result)

            
        # Rank by the full scientist-facing recommendation objective, then by lower off-flavour risk.
        risk_aversion = 1.0 # Default
        texture_aversion = 0.01 # 100 risk = 1.0 score penalty
        for result in results:
            result.ranking_score = compute_ranking_score(
                result,
                risk_aversion=risk_aversion,
                texture_aversion=texture_aversion,
            )
        results.sort(
            key=lambda x: (
                x.ranking_score,
                -x.off_flavour_risk,
            ), 
            reverse=True
        )
        return results
