"""
src/headspace.py

Phase D: Headspace & Volatility Model.
Converts matrix concentrations to predicted air-phase (headspace) concentrations.
"""

import math
import yaml
from pathlib import Path
from typing import Dict, Optional, List

from src.matrix_correction import ProteinType, resolve_matrix_correction

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

    def _extract_properties(self, smiles: str) -> Dict[str, float]:
        """Estimated logP and MW using RDKit for binding characterization."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"logP": 1.0, "MW": 100.0}
        return {
            "logP": Crippen.MolLogP(mol),
            "MW": Descriptors.MolWt(mol)
        }

    def _get_kprot_for_compound(self, name: str, smiles: Optional[str] = None) -> float:
        """
        Calculates empirical protein binding constant (Kprot) for a volatile.
        Phase 1: Partitioning Correction.
        Uses a hydrophobicity-driven model (logP) + molecular weight scaling.
        Kprot = alpha * logP + beta * MW
        """
        if smiles:
            props = self._extract_properties(smiles)
            logp = props["logP"]
            mw = props["MW"]
            
            # Hydrophobic binding (binding affinity increases with logP)
            # Scaling factors estimated from Mcclements et al. (flavor-protein interactions)
            kprot_hydrophobic = max(0, 15.0 * logp)
            # Entropy/Size factor
            kprot_size = 0.05 * mw
            
            return kprot_hydrophobic + kprot_size
            
        n = name.lower()
        # Fallback for named compounds without SMILES lookup
        if "sulfide" in n or "h2s" == n:
            return 0.0      
        elif "thiol" in n:
            return 8.0      
        elif "methional" in n:
            return 2.0      
        elif "pyrazine" in n or "thiazole" in n:
            return 12.0     
        elif "aldehyde" in n or "anal" in n or "fural" in n:
            return 45.0     # Carbonyl-amine covalent binding (Schiff)
        else:
            return 10.0

    def _matrix_retention_fallback(
        self,
        protein_type: Optional[str],
        fat_fraction: float,
        protein_fraction: float,
        denaturation_state: float,
    ) -> float:
        if fat_fraction > 0.0 or protein_fraction > 0.0 or not protein_type:
            return 1.0

        try:
            p_type = ProteinType(protein_type)
        except ValueError:
            return 1.0

        if p_type == ProteinType.FREE_AMINO_ACID:
            return 1.0
        return float(resolve_matrix_correction(p_type, denaturation_state).volatile_retention)

    def get_matrix_ph_release_factor(
        self,
        name: str,
        *,
        protein_type: Optional[str],
        pH: Optional[float],
    ) -> float:
        """
        Empirical pH-dependent headspace release factor for plant matrices.

        The current calibration is intentionally narrow:
        - anchored to pH 6.0 so the existing Pratap-Singh matrix-only baselines stay stable
        - only applied to pea/soy matrix types
        - only applied to acid-sensitive lipid-derived off-flavour markers
        - tuned so pH 4.5 vs 6.5 yields roughly the ~1.6x enhancement reported by
          the Pouvreau pea-isolate benchmark family for aldehydes/furans
        """
        if pH is None or not protein_type:
            return 1.0

        try:
            p_type = ProteinType(protein_type)
        except ValueError:
            return 1.0

        if p_type not in {
            ProteinType.PEA_ISOLATE,
            ProteinType.PEA_CONCENTRATE,
            ProteinType.SOY_ISOLATE,
            ProteinType.SOY_CONCENTRATE,
        }:
            return 1.0

        compound = name.lower()
        acid_sensitive = any(token in compound for token in ["anal", "enal", "furan"])
        if not acid_sensitive:
            return 1.0

        # Reference pH is the ambient plant-isolate baseline already used in the
        # current executable matrix-only benchmarks.
        centered_delta = 6.0 - float(pH)
        factor = math.exp(0.235 * centered_delta)
        return max(0.75, min(1.6, factor))

    def predict_headspace(self, 
                          matrix_concentrations: Dict[str, float], 
                          temp_c: float, 
                          fat_fraction: float = 0.0,
                          protein_fraction: float = 0.0,
                          protein_type: Optional[str] = None,
                          denaturation_state: float = 0.5,
                          pH: Optional[float] = None) -> Dict[str, float]:
        """
        Predicts air-phase concentrations (ppm).
        
        Equation: C_air = C_total * Kaw_eff
        Kaw_eff = Kaw(T) / (1 + Kfat * phi_fat + Kprot * phi_prot)

        If no explicit matrix fractions are provided, `protein_type` can supply
        a conservative fallback retention calibrated from the matrix-correction
        literature estimates already used elsewhere in the project.

        `pH` is optional and currently only applies an empirical plant-matrix
        release correction for acid-sensitive aldehydes/furans in pea/soy systems.
        """
        temp_k = temp_c + 273.15
        air_concs = {}
        matrix_retention = self._matrix_retention_fallback(
            protein_type,
            fat_fraction,
            protein_fraction,
            denaturation_state,
        )
        
        for name, c_total in matrix_concentrations.items():
            entry = self.data.get(name)
            kaw_base = self.get_kaw_at_temp(name, temp_k)
            ph_release_factor = self.get_matrix_ph_release_factor(
                name,
                protein_type=protein_type,
                pH=pH,
            )
            
            if entry:
                k_fat = entry.get("Kfat", 1.0)
                k_prot = self._get_kprot_for_compound(name)
                
                # Effective Kaw accounting for matrix sequestration
                denom = 1.0 + (k_fat * fat_fraction) + (k_prot * protein_fraction)
                kaw_eff = kaw_base / denom
                
                air_concs[name] = c_total * kaw_eff * matrix_retention * ph_release_factor
            else:
                # Basic fallback
                air_concs[name] = c_total * kaw_base * matrix_retention * ph_release_factor
                
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
