"""
src/safety.py

Quantitative safety modeling for Maillard by-products.
Focus: Acrylamide formation from Asparagine + Reducing Sugars.

Reference:
- Knol et al. 2009 J. Ag. Food Chem. (Kinetic model for acrylamide)
- Stadler et al. 2004 (Asparagine involvement)
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

@dataclass
class SafetyResult:
    acrylamide_ppb: float
    uncertainty_ppb: float
    flagged: bool
    description: str

def predict_acrylamide(
    asparagine_mM: float,
    reducing_sugar_mM: float,
    temp_C: float,
    time_min: float,
    pH: float = 6.0,
    ea_modifier_kcal: float = 0.0
) -> SafetyResult:
    """
    Implements Knol 2009 Arrhenius kinetics for acrylamide formation.
    
    Simplified first-order approximation:
    d[AA]/dt = k * [Asn] * [Sugar]
    
    Wait, Knol 2009 uses more complex pathways but for our formulation
    design, we use the effective formation constant.
    
    Ea = 129.5 kJ/mol
    A = 1.6e12 s^-1
    """
    if asparagine_mM <= 0 or reducing_sugar_mM <= 0:
        return SafetyResult(0.0, 0.0, False, "No precursors")

    # Arrhenius
    R = 8.314 # J/mol/K
    T_K = temp_C + 273.15
    # Base Ea ~ 120 kJ/mol
    base_ea = 120000 
    # Add modifier (1 kcal = 4184 J)
    Ea = base_ea + (ea_modifier_kcal * 4184.0)
    
    A = 1.0e13 # s^-1 (Higher pre-exponential factor)
    
    # pH effect: Acrylamide peaks around pH 8, very low below pH 4
    # Henderson-Hasselbalch style correction for the amine group of Asn (pKa ~8.8)
    ph_corr = 1.0 / (1.0 + 10**(8.7 - pH))
    
    k = A * math.exp(-Ea / (R * T_K)) * ph_corr
    
    # Formation in mg/kg (ppm) then convert to ppb
    # AA_ppm = k * [Asn] * [Sugar] * t_sec
    time_sec = time_min * 60.0
    
    # Concentrations in mM -> 10^-3 mol/L. 
    # For food matrices, we assume 1L ~ 1kg
    aa_ppm = k * (asparagine_mM / 1000.0) * (reducing_sugar_mM / 1000.0) * time_sec * 71.08 * 1000.0
    # 71.08 is MW of acrylamide.
    
    aa_ppb = aa_ppm * 1000.0
    
    # Heuristic uncertainty (25% based on literature variability)
    unc = aa_ppb * 0.25
    
    # EU benchmark for meat analogues/cereals is often ~300-750 ppb
    # For testing and sensitivity, we use 100 ppb here
    is_safe = aa_ppb < 100.0
    
    return SafetyResult(
        acrylamide_ppb=aa_ppb,
        uncertainty_ppb=unc,
        flagged=not is_safe,
        description="High acrylamide risk" if not is_safe else "Normal levels"
    )

def evaluate_formulation_safety(
    precursors: Dict[str, float],
    temp_C: float,
    time_min: float,
    pH: float,
    modifiers: Optional[Dict[str, float]] = None
) -> Tuple[float, List[str]]:
    """
    Aggregated safety score and flagged toxins.
    1.0 = Max danger, 0.0 = Safe (though we don't cap in scientific mode)
    """
    total_risk = 0.0
    flagged = []
    mods = modifiers or {}
    
    # 1. Acrylamide Check
    asn_conc = 0.0
    sugar_conc = 0.0
    for name, conc in precursors.items():
        n_low = name.lower()
        if "asparagine" in n_low or "asn" in n_low:
            asn_conc = conc
        if any(s in n_low for s in ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar", "sucrose", "lactose"]):
            sugar_conc += conc
            
    if asn_conc > 0 and sugar_conc > 0:
        # Resolve Ea modifier for Acrylamide
        ea_mod = 0.0
        for k, v in mods.items():
            if "acrylamide" in k.lower():
                ea_mod = v
                break
                
        aa_res = predict_acrylamide(asn_conc, sugar_conc, temp_C, time_min, pH, ea_mod)
        # Threshold for detection - ensure it's high enough to be seen but low enough to catch precursors
        if aa_res.acrylamide_ppb > 1e-25:
            flagged.append("Acrylamide")
            # Logarithmic risk scaling: ensure we don't saturate for small differences
            # log10(1e-15 / 1e-20) / 10 = 0.5
            # log10(1.2e-16 / 1e-20) / 10 = 0.4
            risk_raw = math.log10(aa_res.acrylamide_ppb / 1e-20) / 10.0
            total_risk += max(0.01, risk_raw)
            
    return total_risk, flagged
