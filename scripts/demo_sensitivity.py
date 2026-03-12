import sys
from pathlib import Path
import math

# Setup environment
ROOT = Path("/Users/pabloantoniomorenocasares/Developer/Maillard")
sys.path.insert(0, str(ROOT))

from src.inverse_design import InverseDesigner
from src.conditions import ReactionConditions

def calculate_sensitivity():
    print("--- Sensitivity Analysis: Meaty Target (Pea Protein) ---")
    designer = InverseDesigner(target_tag="meaty")
    cond = ReactionConditions(pH=5.6, temperature_celsius=105.0)
    
    base_form = {
        "name": "Base",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine"],
        "molar_ratios": {"ribose": 0.1, "cysteine": 0.1},
    }
    
    # 1. Base Score
    res_base = designer.evaluate_single(base_form, cond)
    s_base = res_base.target_score
    
    # 2. Add 10% more Cysteine
    form_cys = base_form.copy()
    form_cys["molar_ratios"] = {"ribose": 0.1, "cysteine": 0.11}
    res_cys = designer.evaluate_single(form_cys, cond)
    s_cys = res_cys.target_score
    
    # 3. Add 10% more Ribose
    form_rib = base_form.copy()
    form_rib["molar_ratios"] = {"ribose": 0.11, "cysteine": 0.1}
    res_rib = designer.evaluate_single(form_rib, cond)
    s_rib = res_rib.target_score
    
    # 4. Sensitivity Coefficient S = (%Delta Score) / (%Delta Input)
    # %Delta Input is 0.1 (10% increase)
    sens_cys = ((s_cys - s_base) / s_base) / 0.1
    sens_rib = ((s_rib - s_base) / s_base) / 0.1
    
    print(f"Base Meaty Score: {s_base:.2f}")
    print(f"Cysteine Sensitivity (S_cys): {sens_cys:.3f}")
    print(f"Ribose Sensitivity (S_rib):   {sens_rib:.3f}")
    
    print("\nINTERPRETATION:")
    if sens_rib < 0.001:
        print("SCIENTIFIC FINDING: The system is perfectly Ribose-saturated or Cysteine-limited.")
        print("ACTION: Do NOT add more Ribose. It will not increase flavor at this point.")
    elif sens_cys > sens_rib:
        print(f"Cysteine is {sens_cys/sens_rib:.1f}x more impactful than Ribose at this point.")
        print("ACTION: Prioritize increasing Cysteine dosage for maximum flavor lift.")
    else:
        print(f"Ribose is {sens_rib/sens_cys:.1f}x more impactful than Cysteine at this point.")
        print("ACTION: Prioritize increasing Ribose dosage for maximum flavor lift.")

if __name__ == "__main__":
    calculate_sensitivity()
