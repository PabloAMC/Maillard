"""
src/headspace.py

Phase D: Headspace & Volatility Model.
Converts matrix concentrations to predicted air-phase (headspace) concentrations.
"""

import math
import yaml
from pathlib import Path
from typing import Dict, Any, Optional

class HeadspaceModel:
    """
    Models the partitioning of volatiles between the food matrix and air.
    Accounts for temperature (Van't Hoff) and matrix suppression (lipids/proteins).
    """
    
    def __init__(self, constants_path: str = "data/lit/henry_constants.yml"):
        self.constants_path = Path(constants_path)
        self.data = self._load_constants()
        self.R = 0.008314  # kJ/(mol*K)

    def _load_constants(self) -> Dict[str, Dict]:
        if not self.constants_path.exists():
            return {}
        with open(self.constants_path, "r") as f:
            raw = yaml.safe_load(f)
            return {c["name"]: c for c in raw.get("constants", [])}

    def get_kaw_at_temp(self, name: str, temp_k: float) -> float:
        """
        Calculates the dimensionless air-water partition coefficient at temp_k.
        Uses Van't Hoff / Clausius-Clapeyron extrapolation.
        Kaw(T) = Kaw(Tr) * exp(-dH/R * (1/T - 1/Tr))
        """
        entry = self.data.get(name)
        if not entry:
            return 0.01  # Default fallback volatility
            
        kaw_298 = entry["Kaw_25c"]
        dh = entry["delta_H_sol_kj_mol"]
        
        # Extrapolate: Kaw(T) = Kaw(Tr) * exp(dH_sol/R * (1/temp_k - 1/Tr))
        # Since dH_sol is negative, Kaw increases as temp_k increases.
        exponent = (dh / self.R) * (1.0 / temp_k - 1.0 / 298.15)
        return kaw_298 * math.exp(exponent)

    def predict_headspace(self, 
                          matrix_concentrations: Dict[str, float], 
                          temp_c: float, 
                          fat_fraction: float = 0.0,
                          protein_fraction: float = 0.0) -> Dict[str, float]:
        """
        Predicts air-phase concentrations (ppm).
        
        Equation: C_air = C_total * Kaw_eff
        Kaw_eff = Kaw(T) / (1 + Kfat * phi_fat + Kprot * phi_prot)
        """
        temp_k = temp_c + 273.15
        air_concs = {}
        
        for name, c_total in matrix_concentrations.items():
            entry = self.data.get(name)
            kaw_base = self.get_kaw_at_temp(name, temp_k)
            
            if entry:
                k_fat = entry.get("Kfat", 1.0)
                # Protein binding is simplified here; polar compounds bind more to protein.
                # In this version, we use a generic Kprot = 5.0 for non-H2S compounds.
                k_prot = 5.0 if name != "Hydrogen Sulfide" else 0.0
                
                # Effective Kaw accounting for matrix sequestration
                denom = 1.0 + (k_fat * fat_fraction) + (k_prot * protein_fraction)
                kaw_eff = kaw_base / denom
                
                air_concs[name] = c_total * kaw_eff
            else:
                # Basic fallback
                air_concs[name] = c_total * kaw_base
                
        return air_concs

if __name__ == "__main__":
    model = HeadspaceModel()
    # 100 ppm total hexanal (highly hydrophobic)
    matrix = {"Hexanal": 100.0}
    
    print("Hexanal Headspace Projection at 25°C:")
    # No fat
    c_air_pure = model.predict_headspace(matrix, 25.0, fat_fraction=0.0)["Hexanal"]
    print(f"  Water matrix: {c_air_pure:.4f} ppm")
    
    # 10% fat
    c_air_fat = model.predict_headspace(matrix, 25.0, fat_fraction=0.1)["Hexanal"]
    print(f"  10% Fat matrix: {c_air_fat:.4f} ppm (Suppression: {c_air_pure/c_air_fat:.1f}x)")
