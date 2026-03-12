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
    def simulate_yeast_cleaning(molar_ratios: Dict[str, float]) -> Dict[str, float]:
        """
        Simulates yeast action (e.g., Saccharomyces cerevisiae) which can 
        metabolize beany aldehydes (hexanal) into less impactful alcohols.
        """
        new_ratios = molar_ratios.copy()
        # Yeast ADH (alcohol dehydrogenase) converts hexanal -> hexanol
        for k in list(new_ratios.keys()):
            if "hexanal" in k.lower():
                val = new_ratios[k]
                # Assume 80% conversion in a standard fermentation step
                new_ratios[k] = val * 0.2
                new_ratios["hexanol"] = new_ratios.get("hexanol", 0.0) + val * 0.8
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
        if "yeast_fermentation" in interventions:
            ratios = self.simulate_yeast_cleaning(ratios)
        if "protease_hydrolysis" in interventions:
            ratios = self.simulate_protease_hydrolysis(ratios)
        return ratios
