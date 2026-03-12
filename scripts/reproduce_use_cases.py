import sys
from pathlib import Path

# Setup environment
ROOT = Path("/Users/pabloantoniomorenocasares/Developer/Maillard")
sys.path.insert(0, str(ROOT))

from src.pre_processor import PreProcessor
from src.bayesian_optimizer import FormulationOptimizer
from src.inverse_design import InverseDesigner
from src.conditions import ReactionConditions

def run_pea_protein_reproduction():
    print("\n[REPRODUCING] Premium Roast Pea Protein")
    pp = PreProcessor()
    raw_ppi = {"Hexanal": 1.2, "Nonanal": 0.8}
    cleaned = pp.apply(raw_ppi, [{"yeast_fermentation": {"time_hours": 4}}])
    
    designer = InverseDesigner(target_tag="meaty")
    cond = ReactionConditions(pH=5.6, temperature_celsius=105.0)
    best_form = {
        "name": "PeaBest",
        "sugars": ["glucose", "ribose"],
        "amino_acids": ["lysine", "cysteine"],
        "molar_ratios": {"glucose": 0.5, "ribose": 0.2, "lysine": 0.5, "cysteine": 0.2},
    }
    res = designer.evaluate_single(best_form, cond)
    print(f"Meaty Score: {res.target_score:.2f}")
    print(f"Detected: {res.detected_targets}")

def run_nutty_reproduction():
    print("\n[REPRODUCING] Roasted Nutty Profile")
    designer = InverseDesigner(target_tag="roasted")
    # Increase Temp to 160C to get over the pyrazine barrier
    cond = ReactionConditions(pH=8.5, temperature_celsius=160.0)
    form = {
        "name": "RoastedForm",
        "sugars": ["glucose"],
        "amino_acids": ["glycine", "alanine"],
        "molar_ratios": {"glucose": 1.0, "glycine": 0.5, "alanine": 0.5},
    }
    res = designer.evaluate_single(form, cond)
    print(f"Roasted Score: {res.target_score:.2f}")
    print(f"Detected: {res.detected_targets}")

def run_toxicity_decoupling():
    print("\n[REPRODUCING] Toxicity-Flavor Decoupling in High-Heat Searing")
    designer = InverseDesigner(target_tag="roasted") # Searing
    
    # 1. Baseline: High pH (Risky)
    baseline_cond = ReactionConditions(pH=8.0, temperature_celsius=175.0)
    # 2. Optimized: Low pH (Safe)
    safe_cond = ReactionConditions(pH=5.8, temperature_celsius=175.0)
    
    form = {
        "name": "SearingPatty",
        "sugars": ["glucose"],
        "amino_acids": ["glycine", "asparagine"],
        "molar_ratios": {"glucose": 1.0, "glycine": 0.5, "asparagine": 0.5}
    }
    
    res_risky = designer.evaluate_single(form, baseline_cond)
    res_safe = designer.evaluate_single(form, safe_cond)
    
    print(f"--- Baseline (pH 8.0) ---")
    print(f"Roasted Score: {res_risky.target_score:.2f}")
    print(f"Safety Penalty: {res_risky.safety_score:.2f}")
    print(f"Flagged Toxics: {res_risky.flagged_toxics}")
    
    print(f"--- Optimized (pH 5.8) ---")
    print(f"Roasted Score: {res_safe.target_score:.2f}")
    print(f"Safety Penalty: {res_safe.safety_score:.2f}")
    print(f"Flagged Toxics: {res_safe.flagged_toxics}")

if __name__ == "__main__":
    run_pea_protein_reproduction()
    run_nutty_reproduction()
    run_toxicity_decoupling()
