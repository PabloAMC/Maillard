import optuna
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[1]

from src.inverse_design import InverseDesigner  # noqa: E402
from src.smirks_engine import ReactionConditions  # noqa: E402
from src.pre_processor import PreProcessor  # noqa: E402

class FormulationOptimizer:
    """
    Bayesian Optimizer for Maillard reaction formulations.
    Searches the continuous parameter space (concentrations, pH, temp) 
    to maximize the Pareto-ranked sensory outcome minus safety penalties.
    """
    def __init__(self, target_tag: str, minimize_tag: str = "beany", risk_aversion: float = 1.0, protein_type: str = "free", denaturation_state: float = 0.5):
        self.target_tag = target_tag
        self.minimize_tag = minimize_tag
        self.risk_aversion = risk_aversion
        self.protein_type = protein_type
        self.denaturation_state = denaturation_state
        self.study = None

    def objective(self, trial: optuna.Trial, fixed_sugars: List[str], fixed_amino_acids: List[str], fixed_lipids: Optional[List[str]] = None) -> float:
        fixed_lipids = fixed_lipids or []
        
        # 1. Sample continuous parameters
        sugar_conc = trial.suggest_float("sugar_conc", 0.01, 1.0, log=True)
        
        # Phase N: Sample independent concentrations based on amino acid class
        aa_conc_sulfur = trial.suggest_float("aa_conc_sulfur", 0.01, 1.0, log=True)
        aa_conc_branched = trial.suggest_float("aa_conc_branched", 0.01, 1.0, log=True)
        aa_conc_basic = trial.suggest_float("aa_conc_basic", 0.01, 1.0, log=True)
        aa_conc_other = trial.suggest_float("aa_conc_other", 0.01, 1.0, log=True)

        ph = trial.suggest_float("ph", 3.0, 9.0)
        temp = trial.suggest_float("temp", 100.0, 200.0)
        aw = trial.suggest_float("aw", 0.3, 0.95)
        # Time is logged for future kinetic-model expansion, currently not used in FAST bounds
        time_mins = trial.suggest_float("time_minutes", 10.0, 120.0)
        
        # Phase 20: Suggest interventions from library
        import yaml
        LIBRARY_PATH = ROOT / "data" / "interventions.yml"
        if LIBRARY_PATH.exists():
            with open(LIBRARY_PATH, "r") as f:
                lib_data = yaml.safe_load(f)
                agents = [a["name"] for a in lib_data.get("interventions", [])] + ["none"]
        else:
            agents = ["none"]
            
        agent = trial.suggest_categorical("intervention_agent", agents)
        agent_dose = trial.suggest_float("intervention_dose", 0.0, 1.0) if agent != "none" else 0.0
        
        interventions = [{"name": agent, "dose": agent_dose}] if agent != "none" else []
        
        # Phase 21: Pre-processing options (could also be part of the trial)
        pre_processing = trial.suggest_categorical("pre_processing", ["none", "yeast_fermentation", "protease_hydrolysis", "both"])
        pre_steps = []
        if pre_processing == "both":
            pre_steps = ["yeast_fermentation", "protease_hydrolysis"]
        elif pre_processing != "none":
            pre_steps = [pre_processing]

        # 2. Setup the single evaluation condition
        cond = ReactionConditions(
            pH=ph, 
            temperature_celsius=temp, 
            water_activity=aw,
            protein_type=self.protein_type
        )
        
        # 3. Create a custom grid override
        molar_ratios = {}
        for s in fixed_sugars:
            molar_ratios[s] = sugar_conc
            
        for a in fixed_amino_acids:
            name_lower = a.lower()
            if name_lower in ["cysteine", "methionine"]:
                molar_ratios[a] = aa_conc_sulfur
            elif name_lower in ["leucine", "isoleucine", "valine"]:
                molar_ratios[a] = aa_conc_branched
            elif name_lower in ["lysine", "arginine", "histidine"]:
                molar_ratios[a] = aa_conc_basic
            else:
                molar_ratios[a] = aa_conc_other
                
        for L in fixed_lipids:
            molar_ratios[L] = 0.1 # Example fixed trace lipid
            
        # Phase 21: Apply Pre-Processing
        processor = PreProcessor()
        molar_ratios = processor.apply(molar_ratios, pre_steps)

        formulation = {
            "name": f"Trial_{trial.number}",
            "sugars": fixed_sugars,
            "amino_acids": fixed_amino_acids,
            "lipids": fixed_lipids,
            "molar_ratios": molar_ratios,
            "ph": ph,
            "temp": temp,
            "aw": aw,
            "time_minutes": time_mins,
            "interventions": interventions,
            "protein_type": self.protein_type,
            "denaturation_state": self.denaturation_state
        }
        
        # 4. Evaluate using the robust pipeline without mutating global state (R.8 fix)
        designer = InverseDesigner(self.target_tag, self.minimize_tag)
        res = designer.evaluate_single(formulation, cond)
        
        # 5. Objective Calculation
        # Maximize: target_score - risk_aversion * safety_score
        # But we also want to penalize extreme off-flavours if they spike too high
        off_flavour_penalty = res.off_flavour_risk * 0.5 
        
        # The following lines are from the user's provided change.
        # It seems there was an intention to use 'best_res' from a list of 'results',
        # but in this 'objective' method, we have a single 'res' object.
        # Assuming 'best_res' should refer to 'res' in this context for the subsequent calculations.
        # best_res = sorted(results, key=lambda x: x.target_score, reverse=True)[0] # This line is commented out as 'results' is not defined.
            
        # Uncertainty-aware scoring (Fix 5)
        # High span uncertainty (from missing/low-quality barriers) penalizes the score
        # to prevent the optimizer from exploiting model blind spots.
        
        # Average uncertainty across detected targets
        # total_unc = sum(t["span_uncertainty"] for t in res.predicted_ppb.values() if isinstance(t, dict) and "span_uncertainty" in t)
        # wait, best_res.predicted_ppb is the conc_map from evaluate_all. 
        # It doesn't have the uncertainty. 
        # I need to expose average uncertainty in FormulationResult.
        
        target_val = res.target_score
        safety_penalty = self.risk_aversion * res.safety_score
        off_flavor_penalty = 0.5 * res.off_flavour_risk
        
        # Heuristic uncertainty penalty: -0.1 per kcal of span uncertainty
        unc_penalty = res.avg_uncertainty * 0.1
        
        final_objective = target_val - safety_penalty - off_flavor_penalty - unc_penalty
        
        trial.set_user_attr("target_score", target_val)
        trial.set_user_attr("safety_score", res.safety_score)
        trial.set_user_attr("off_flavour_risk", res.off_flavour_risk)
        trial.set_user_attr("avg_uncertainty", res.avg_uncertainty)
        trial.set_user_attr("flagged_toxics", res.flagged_toxics)
        
        return final_objective

    def optimize(self, 
                 fixed_sugars: List[str], 
                 fixed_amino_acids: List[str], 
                 fixed_lipids: Optional[List[str]] = None,
                 n_trials: int = 50) -> optuna.Study:
        """Runs the Expected Improvement (EI) Bayesian optimization."""
        # Suppress noisy standard INFO level logging from optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        
        # Create a study maximizing the combined score
        self.study = optuna.create_study(direction="maximize")
        
        self.study.optimize(
            lambda trial: self.objective(trial, fixed_sugars, fixed_amino_acids, fixed_lipids), 
            n_trials=n_trials,
            catch=(optuna.exceptions.TrialPruned,)
        )
        return self.study
