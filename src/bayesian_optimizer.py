import optuna
from typing import List, Optional

from src.inverse_design import InverseDesigner  # noqa: E402
from src.smirks_engine import ReactionConditions  # noqa: E402
from src.pre_processor import PreProcessor  # noqa: E402

class FormulationOptimizer:
    """
    Bayesian Optimizer for Maillard reaction formulations.
    Searches the continuous parameter space (concentrations, pH, temp) 
    to maximize the Pareto-ranked sensory outcome minus safety penalties.
    """
    def __init__(self, target_tag: str, minimize_tag: str = "beany", risk_aversion: float = 1.0):
        self.target_tag = target_tag
        self.minimize_tag = minimize_tag
        self.risk_aversion = risk_aversion
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
        
        # Phase 20: Suggest interventions
        intervention = trial.suggest_categorical("intervention", ["none", "calcium_carbonate", "rosemary_extract"])
        interventions = [intervention] if intervention != "none" else []
        
        # Phase 21: Pre-processing options (could also be part of the trial)
        pre_processing = trial.suggest_categorical("pre_processing", ["none", "yeast_fermentation", "protease_hydrolysis", "both"])
        pre_steps = []
        if pre_processing == "both":
            pre_steps = ["yeast_fermentation", "protease_hydrolysis"]
        elif pre_processing != "none":
            pre_steps = [pre_processing]

        # 2. Setup the single evaluation condition
        cond = ReactionConditions(pH=ph, temperature_celsius=temp, water_activity=aw)
        
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
            "interventions": interventions
        }
        
        # 4. Evaluate using the robust pipeline without mutating global state (R.8 fix)
        designer = InverseDesigner(self.target_tag, self.minimize_tag)
        res = designer.evaluate_single(formulation, cond)
        
        # 5. Objective Calculation
        # Maximize: target_score - risk_aversion * safety_score
        # But we also want to penalize extreme off-flavours if they spike too high
        off_flavour_penalty = res.off_flavour_risk * 0.5 
        
        score = res.target_score - (self.risk_aversion * res.safety_score) - off_flavour_penalty
        
        # Store useful metadata for the CLI output
        trial.set_user_attr("target_score", res.target_score)
        trial.set_user_attr("safety_score", res.safety_score)
        trial.set_user_attr("off_flavour_risk", res.off_flavour_risk)
        trial.set_user_attr("flagged_toxics", list(res.flagged_toxics))
        
        return score

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
