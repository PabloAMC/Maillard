"""
src/pre_processor.py — Enzymatic and Biological Pre-Processing Simulation

This module simulates the "matrix cleanup" phase (e.g., fermentation, 
enzymatic hydrolysis) that occurs before thermal processing.
"""

from typing import Dict

class PreProcessor:
    """
    Simulates changes to the precursor pool based on pre-processing interventions.
    """
    
    @staticmethod
    def simulate_yeast_cleaning(molar_ratios: Dict[str, float], efficiency: float = 0.8) -> Dict[str, float]:
        """
        Simulates yeast action (e.g., Saccharomyces cerevisiae) which can 
        metabolize beany/rancid aldehydes into less impactful alcohols.
        
        Library expansion (Phase 2):
        - Hexanal -> Hexanol (Beany -> Mild)
        - Nonanal -> Nonanol (Green/Fatty -> Mild)
        - 2,4-Decadienal -> 2,4-Decadienol (Deep Fatty/Rancid -> Mild)
        """
        new_ratios = molar_ratios.copy()
        target_map = {
            "hexanal": "hexanol",
            "nonanal": "nonanol",
            "decadienal": "decadienol"
        }
        
        for k in list(new_ratios.keys()):
            k_lower = k.lower()
            for aldehyde, alcohol in target_map.items():
                if aldehyde in k_lower:
                    val = new_ratios[k]
                    # Apply biotransformation efficiency
                    new_ratios[k] = val * (1.0 - efficiency)
                    new_ratios[alcohol] = new_ratios.get(alcohol, 0.0) + val * efficiency
        return new_ratios

    @staticmethod
    def simulate_protease_hydrolysis(molar_ratios: Dict[str, float]) -> Dict[str, float]:
        """
        Simulates protease treatment (e.g., Papain, Alcalase) which increases 
        free amino acid concentrations by hydrolyzing proteins.
        """
        new_ratios = molar_ratios.copy()
        # Heuristic: 2x increase in all free amino acids as they are liberated from the matrix
        # (In a real model, this would subtract from a 'protein_bound' pool)
        for k in list(new_ratios.keys()):
            # Detect common amino acids
            if any(aa in k.lower() for aa in ["lysine", "leucine", "isoleucine", "valine", "phenylalanine", "methionine", "cysteine"]):
                new_ratios[k] *= 2.0
        return new_ratios

    def apply(self, molar_ratios: Dict[str, float], interventions: list) -> Dict[str, float]:
        """Apply a set of pre-processing interventions."""
        ratios = molar_ratios
        
        # Check for fermentation specific intervention with time/efficiency
        for inter in interventions:
            if "yeast_fermentation" in str(inter):
                eff = 0.8
                # Handle both {"yeast_fermentation": {"time_hours": 5}} and {"yeast_fermentation": True}
                if isinstance(inter, dict):
                    params = inter.get("yeast_fermentation")
                    if isinstance(params, dict):
                        eff = params.get("efficiency", 0.8)
                        if "time_hours" in params:
                            import math
                            t = params["time_hours"]
                            eff = 1.0 - math.exp(-0.4 * t)
                    else:
                        eff = inter.get("efficiency", 0.8)
                
                ratios = self.simulate_yeast_cleaning(ratios, efficiency=eff)
                break
                
        if "protease_hydrolysis" in interventions:
            ratios = self.simulate_protease_hydrolysis(ratios)
        return ratios
