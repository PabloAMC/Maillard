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
    Implements formation-elimination kinetics for acrylamide (Knol/Parker model).
    
    Kinetics:
    dA/dt = kf * [Asn] * [Sugar] - ke * [A]
    
    Analytical Solution:
    A(t) = (kf * [Asn] * [Sugar] / ke) * (1 - exp(-ke * t))
    
    References:
    - Knol et al. 2009 (Formation kinetics)
    - Parker et al. 2012 (Elimination/Degradation kinetics)
    """
    if asparagine_mM <= 0 or reducing_sugar_mM <= 0:
        return SafetyResult(0.0, 0.0, False, "No precursors")

    # Constants
    R = 8.314 # J/mol/K
    T_K = temp_C + 273.15
    time_sec = time_min * 60.0
    MW_AA = 71.08
    
    # 1. Formation (kf)
    # Base Ea for formation from Knol 2009 (~129 kJ/mol)
    Ea_f = 129000.0 + (ea_modifier_kcal * 4184.0)
    A_f = 1.6e13 # L/mol/s 
    
    # pH effect on formation: Asparagine amine nucleophilicity (pKa ~8.8)
    # Most reactive in slightly alkaline, sharply drops below pH 5
    f_ph = 1.0 / (1.0 + 10**(8.8 - pH))
    kf = A_f * math.exp(-Ea_f / (R * T_K)) * f_ph
    
    # 2. Elimination (ke)
    # Acrylamide degrades at high T. Ea_e typically ~90-110 kJ/mol
    Ea_e = 105000.0 
    A_e = 5.0e8 # s^-1
    ke = A_e * math.exp(-Ea_e / (R * T_K))
    
    # 3. Resolve Concentrations (mM to mol/L)
    asn_molar = asparagine_mM / 1000.0
    sugar_molar = reducing_sugar_mM / 1000.0
    
    # Avoid division by zero if ke is extremely small
    if ke < 1e-12:
        # Fallback to simple first-order formation
        aa_molar = kf * asn_molar * sugar_molar * time_sec
    else:
        # Analytical solution
        aa_molar = (kf * asn_molar * sugar_molar / ke) * (1 - math.exp(-ke * time_sec))
        
    # Convert mol/L -> ppm (mg/kg assuming rho=1) -> ppb
    aa_ppm = aa_molar * MW_AA * 1000.0
    aa_ppb = aa_ppm * 1000.0
    
    # Heuristic uncertainty (25% based on literature variability)
    unc = aa_ppb * 0.25
    
    # EU benchmark for meat analogues/cereals is often ~300-750 ppb
    # We use a conservative 100 ppb threshold for flagging
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
