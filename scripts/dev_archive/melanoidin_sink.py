from typing import Dict
import math

# Reference constants for basic polymerization kinetics
# Melanoidin formation is generally a high-barrier condensation reaction
# Favored at high temperatures and longer cooking times.
MELANOIDIN_KINETICS = {
    # Alpha-dicarbonyl polymerization (e.g. Diacetyl, Methylglyoxal)
    "alpha_dicarbonyl": {
        "Ea_kJ": 115.0,  # High barrier for non-enzymatic browning polymers
        "ln_A": 22.5     # Pre-exponential factor derived for late-stage browning
    },
    # Furan/Furfural polymerization (often condenses with amines to form pigment)
    "furan": {
        "Ea_kJ": 105.0,
        "ln_A": 20.0
    },
    # Strecker aldehydes (often evaporate, but can condense if trapped)
    "strecker_aldehyde": {
        "Ea_kJ": 125.0,
        "ln_A": 24.0
    },
    # Generic polymer sink for other intermediates
    "generic_intermediate": {
        "Ea_kJ": 120.0,
        "ln_A": 21.0
    }
}

def classify_intermediate_for_sink(compound_name: str) -> str:
    """Classify the chemical name into a functional category matching our kinetic sink model."""
    name = compound_name.lower()
    
    if any(x in name for x in ["furan", "furfural", "hydroxymethylfurfural"]):
        return "furan"
    elif any(x in name for x in ["diacetyl", "methylglyoxal", "glyoxal", "butanedione"]):
        return "alpha_dicarbonyl"
    elif any(x in name for x in ["aldehyde", "anal", "propanal", "butanal", "pentanal", "hexanal"]):
        # Typical strecker or lipid oxidation aldehydes
        return "strecker_aldehyde"
    else:
        return "generic_intermediate"

def calculate_melanoidin_trapping_fraction(
    compound_name: str,
    temperature_celsius: float,
    time_minutes: float,
    pH: float
) -> float:
    """
    Calculates the fraction of an intermediate compound that *remains* 
    after continuous polymerization into melanoidin sinks.
    
    Returns:
        Fraction remaining (0.0 to 1.0), where 1.0 means no polymerization.
    """
    if time_minutes <= 0.0:
        return 1.0
        
    intermediate_class = classify_intermediate_for_sink(compound_name)
    kinetics = MELANOIDIN_KINETICS[intermediate_class]
    
    T_K = temperature_celsius + 273.15
    R_kJ = 0.008314  # Gas constant in kJ/(mol*K)
    
    # Calculate Arrhenius rate constant k (min^-1)
    # k = A * exp(-Ea / RT) = exp(ln A - Ea / RT)
    ln_k = kinetics["ln_A"] - (kinetics["Ea_kJ"] / (R_kJ * T_K))
    try:
        k_poly = math.exp(ln_k)
    except OverflowError:
        k_poly = float('inf')  # Reaction is essentially instantaneous
        
    # pH effect on polymerization:
    # Browning (melanoidin formation) is dramatically faster at alkaline pH.
    # We apply an empirical exponential multiplier anchored at pH 6.0.
    ph_multiplier = math.exp(0.5 * (pH - 6.0))
    k_poly *= ph_multiplier
    
    # First-order decay: [A] = [A]0 * exp(-k * t)
    # The fraction remaining is exp(-k * t)
    exponent = -k_poly * time_minutes
    
    # Cap to avoid math domain errors
    if exponent < -50:
        return 0.0  # Fully consumed by the sink
        
    return math.exp(exponent)
