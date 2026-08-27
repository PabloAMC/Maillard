import optuna
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

ROOT = Path(__file__).resolve().parents[1]

from src.pipeline import MaillardPipeline  # noqa: E402
from src.smirks_engine import ReactionConditions  # noqa: E402
from src.pre_processor import PreProcessor  # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# Trial-parameter contract
#
# The Optuna search space is deliberately CLASS-based: one concentration knob
# per amino-acid class, not one per precursor. That means `trial.params` keys
# ("aa_conc_sulfur", "ph", …) are NOT precursor names, and feeding them to the
# pipeline as `molar_ratios` is wrong twice over — the concentrations never
# match a precursor label (so everything falls back to a flat 1.0) and "ph"
# substring-matches "L-Phenylalanine", leaking the trial pH in as that amino
# acid's concentration. The mapping below is the single source of truth used by
# both `objective()` and `scripts/optimize_formulation.py`.
# ─────────────────────────────────────────────────────────────────────────────

#: Trial params that carry a precursor concentration, in millimolar (mM) — the
#: unit `molar_ratios` values are consumed as downstream (see src/safety.py).
CONCENTRATION_PARAM_KEYS = (
    "sugar_conc",
    "aa_conc_sulfur",
    "aa_conc_branched",
    "aa_conc_basic",
    "aa_conc_other",
)
#: Trial params that are process conditions, not concentrations.
CONDITION_PARAM_KEYS = ("ph", "temp", "aw", "time_minutes")
#: Trial params whose value is a string; float() on these is what made the CLI
#: crash unconditionally at the end of every optimization run.
CATEGORICAL_PARAM_KEYS = ("intervention_agent", "pre_processing")
#: Remaining numeric knobs that are neither concentrations nor conditions.
AUXILIARY_PARAM_KEYS = ("intervention_dose",)

_SULFUR_AMINO_ACIDS = ("cysteine", "methionine")
_BRANCHED_AMINO_ACIDS = ("leucine", "isoleucine", "valine")
_BASIC_AMINO_ACIDS = ("lysine", "arginine", "histidine")

#: Fixed trace loading (mM) applied to every lipid off-flavour precursor.
FIXED_LIPID_CONCENTRATION_MM = 0.1


def amino_acid_concentration_param(name: str) -> str:
    """Return the trial-parameter key that governs this amino acid's loading."""
    lowered = str(name).strip().lower()
    if lowered in _SULFUR_AMINO_ACIDS:
        return "aa_conc_sulfur"
    if lowered in _BRANCHED_AMINO_ACIDS:
        return "aa_conc_branched"
    if lowered in _BASIC_AMINO_ACIDS:
        return "aa_conc_basic"
    return "aa_conc_other"


def build_molar_ratios(
    params: Mapping[str, Any],
    fixed_sugars: List[str],
    fixed_amino_acids: List[str],
    fixed_lipids: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Map class-level trial concentrations onto real precursor names (mM)."""
    molar_ratios: Dict[str, float] = {}
    sugar_conc = float(params.get("sugar_conc", 1.0))
    for sugar in fixed_sugars:
        molar_ratios[sugar] = sugar_conc
    for amino_acid in fixed_amino_acids:
        key = amino_acid_concentration_param(amino_acid)
        molar_ratios[amino_acid] = float(params.get(key, 1.0))
    for lipid in (fixed_lipids or []):
        molar_ratios[lipid] = FIXED_LIPID_CONCENTRATION_MM
    return molar_ratios


def interventions_from_params(params: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Rebuild the intervention list a trial actually used."""
    agent = params.get("intervention_agent", "none")
    if not agent or str(agent) == "none":
        return []
    return [{"name": str(agent), "dose": float(params.get("intervention_dose", 0.0) or 0.0)}]


def pre_processing_steps_from_params(params: Mapping[str, Any]) -> List[str]:
    """Rebuild the pre-processing steps a trial actually used."""
    choice = str(params.get("pre_processing", "none"))
    if choice == "both":
        return ["yeast_fermentation", "protease_hydrolysis"]
    if choice in {"", "none"}:
        return []
    return [choice]


def formulation_from_params(
    params: Mapping[str, Any],
    *,
    name: str,
    fixed_sugars: List[str],
    fixed_amino_acids: List[str],
    fixed_lipids: Optional[List[str]] = None,
    protein_type: str = "free",
    denaturation_state: Optional[float] = None,
    pre_processor: Optional[PreProcessor] = None,
) -> Dict[str, Any]:
    """Rebuild the exact formulation dict a set of trial params describes.

    Used by `objective()` to evaluate a trial and by the CLI to re-evaluate the
    winning trial. Because both go through here, the "Optimized Parameter Set"
    the CLI prints is the formulation that was actually scored.
    """
    fixed_lipids = fixed_lipids or []
    molar_ratios = build_molar_ratios(params, fixed_sugars, fixed_amino_acids, fixed_lipids)
    processor = pre_processor or PreProcessor()
    molar_ratios = processor.apply(molar_ratios, pre_processing_steps_from_params(params))

    return {
        "name": name,
        "sugars": list(fixed_sugars),
        "amino_acids": list(fixed_amino_acids),
        "lipids": list(fixed_lipids),
        "molar_ratios": molar_ratios,
        "ph": float(params.get("ph", 6.0)),
        "temp": float(params.get("temp", 150.0)),
        "aw": float(params.get("aw", 0.8)),
        "time_minutes": float(params.get("time_minutes", 60.0)),
        "interventions": interventions_from_params(params),
        "protein_type": protein_type,
        "denaturation_state": denaturation_state,
    }


class FormulationOptimizer:
    """
    Bayesian Optimizer for Maillard reaction formulations.
    Searches the continuous parameter space (concentrations, pH, temp) 
    to maximize the Pareto-ranked sensory outcome minus safety penalties.
    """
    def __init__(self, target_tag: str, minimize_tag: str = "beany", risk_aversion: float = 1.0, protein_type: str = "free", denaturation_state: Optional[float] = None, seed: Optional[int] = None):
        self.target_tag = target_tag
        self.minimize_tag = minimize_tag
        self.risk_aversion = risk_aversion
        self.protein_type = protein_type
        self.denaturation_state = denaturation_state
        self.seed = seed
        self.study = None

    def objective(self, trial: optuna.Trial, fixed_sugars: List[str], fixed_amino_acids: List[str], fixed_lipids: Optional[List[str]] = None) -> float:
        fixed_lipids = fixed_lipids or []
        
        # 1. Sample continuous parameters (concentrations are in mM).
        # Values are read back off `trial.params` by `formulation_from_params`
        # below, so the sampled space and the evaluated formulation cannot drift.
        trial.suggest_float("sugar_conc", 0.01, 1.0, log=True)

        # Phase N: Sample independent concentrations based on amino acid class
        trial.suggest_float("aa_conc_sulfur", 0.01, 1.0, log=True)
        trial.suggest_float("aa_conc_branched", 0.01, 1.0, log=True)
        trial.suggest_float("aa_conc_basic", 0.01, 1.0, log=True)
        trial.suggest_float("aa_conc_other", 0.01, 1.0, log=True)

        ph = trial.suggest_float("ph", 3.0, 9.0)
        temp = trial.suggest_float("temp", 100.0, 200.0)
        aw = trial.suggest_float("aw", 0.3, 0.95)
        trial.suggest_float("time_minutes", 10.0, 120.0)


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
        if agent != "none":
            trial.suggest_float("intervention_dose", 0.0, 1.0)

        # Phase 21: Pre-processing options (could also be part of the trial)
        trial.suggest_categorical("pre_processing", ["none", "yeast_fermentation", "protease_hydrolysis", "both"])

        # 2. Setup the single evaluation condition
        cond = ReactionConditions(
            pH=ph,
            temperature_celsius=temp,
            water_activity=aw,
            protein_type=self.protein_type
        )

        # 3. Create a custom grid override.
        # Built through the shared mapper so the CLI can rebuild this exact
        # formulation from `trial.params` when it re-scores the winner.
        formulation = formulation_from_params(
            trial.params,
            name=f"Trial_{trial.number}",
            fixed_sugars=fixed_sugars,
            fixed_amino_acids=fixed_amino_acids,
            fixed_lipids=fixed_lipids,
            protein_type=self.protein_type,
            denaturation_state=self.denaturation_state,
        )


        # 4. Evaluate using the robust pipeline without mutating global state (R.8 fix)
        designer = MaillardPipeline(self.target_tag, self.minimize_tag)
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
        quality_penalty = res.meaty_quality_penalty
        strecker_penalty = res.strecker_gap_penalty
        pyrazine_penalty = res.pyrazine_penalty
        furanone_penalty = res.furanone_penalty
        
        # Heuristic uncertainty penalty: -0.1 per kcal of span uncertainty
        unc_penalty = res.avg_uncertainty * 0.1
        
        final_objective = target_val - safety_penalty - off_flavor_penalty - quality_penalty - strecker_penalty - pyrazine_penalty - furanone_penalty - unc_penalty
        
        trial.set_user_attr("target_score", target_val)
        trial.set_user_attr("safety_score", res.safety_score)
        trial.set_user_attr("off_flavour_risk", res.off_flavour_risk)
        trial.set_user_attr("meaty_quality_penalty", res.meaty_quality_penalty)
        trial.set_user_attr("mft_to_furfural_ratio", res.mft_to_furfural_ratio)
        trial.set_user_attr("strecker_gap_penalty", res.strecker_gap_penalty)
        trial.set_user_attr("pyrazine_penalty", res.pyrazine_penalty)
        trial.set_user_attr("pyrazine_burden", res.pyrazine_burden)
        trial.set_user_attr("furanone_penalty", res.furanone_penalty)
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
        sampler = optuna.samplers.TPESampler(seed=self.seed) if self.seed is not None else None
        self.study = optuna.create_study(direction="maximize", sampler=sampler)
        
        self.study.optimize(
            lambda trial: self.objective(trial, fixed_sugars, fixed_amino_acids, fixed_lipids), 
            n_trials=n_trials,
            catch=(optuna.exceptions.TrialPruned,)
        )
        return self.study
