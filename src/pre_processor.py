"""
src/pre_processor.py — Enzymatic and Biological Pre-Processing Simulation

This module simulates the "matrix cleanup" phase (e.g., fermentation, 
enzymatic hydrolysis) that occurs before thermal processing.
"""

from typing import Any, Dict, Optional, Tuple

#: The pre-processing interventions this module knows how to simulate.
#: Matching is exact against these tokens — never by substring — so that
#: ``"no_yeast_fermentation"`` (a deliberate opt-out spelling) and
#: ``{"yeast_fermentation": False}`` do NOT switch the intervention on.
KNOWN_INTERVENTIONS = ("yeast_fermentation", "protease_hydrolysis")


def _falsy_flag(value: Any) -> bool:
    """True when ``value`` explicitly disables an intervention."""
    if value is None:
        return True
    if isinstance(value, bool):
        return not value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value) == 0.0
    if isinstance(value, str):
        return value.strip().lower() in {"", "false", "no", "off", "none", "0"}
    return False


def resolve_intervention(entry: Any, name: str) -> Tuple[bool, Optional[Dict[str, Any]]]:
    """Decide whether ``entry`` activates the intervention called ``name``.

    Accepted spellings (all exact-token, never substring):

    * ``"yeast_fermentation"`` — bare string
    * ``{"yeast_fermentation": True}`` / ``{"yeast_fermentation": {...}}``
    * ``{"name": "yeast_fermentation", ...}`` — the library-style record used
      by ``data/interventions.yml``

    Returns ``(active, params)`` where ``params`` is the parameter mapping when
    one was supplied (``{"time_hours": 5}``, ``{"efficiency": 0.6}``, …).
    """
    target = str(name).strip().lower()

    if isinstance(entry, str):
        return (entry.strip().lower() == target, None)

    if isinstance(entry, dict):
        # library-style record: {"name": ..., "dose": ...}
        entry_name = entry.get("name")
        if isinstance(entry_name, str):
            if entry_name.strip().lower() != target:
                return (False, None)
            return (not _falsy_flag(entry.get("enabled", True)), dict(entry))

        # flag-style record: {"yeast_fermentation": True | {...} | False}
        for key, value in entry.items():
            if str(key).strip().lower() != target:
                continue
            if _falsy_flag(value):
                return (False, None)
            if isinstance(value, dict):
                return (True, dict(value))
            return (True, dict(entry))
        return (False, None)

    return (False, None)


def intervention_is_active(interventions: Any, name: str) -> bool:
    """True when any entry in ``interventions`` activates ``name``."""
    for entry in interventions or []:
        active, _ = resolve_intervention(entry, name)
        if active:
            return True
    return False


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
        """Apply a set of pre-processing interventions.

        Matching is exact-token (see :func:`resolve_intervention`). Before
        2026-08-27 this used ``"yeast_fermentation" in str(inter)``, which
        turned the intervention ON for ``"no_yeast_fermentation"`` and for
        ``{"yeast_fermentation": False}`` — i.e. explicitly disabling it
        enabled it.
        """
        ratios = molar_ratios

        # Check for fermentation specific intervention with time/efficiency
        for inter in interventions or []:
            active, params = resolve_intervention(inter, "yeast_fermentation")
            if not active:
                continue
            eff = 0.8
            if params:
                if "efficiency" in params:
                    eff = float(params["efficiency"])
                if "time_hours" in params:
                    import math
                    t = float(params["time_hours"])
                    eff = 1.0 - math.exp(-0.4 * t)
            ratios = self.simulate_yeast_cleaning(ratios, efficiency=eff)
            break

        if intervention_is_active(interventions, "protease_hydrolysis"):
            ratios = self.simulate_protease_hydrolysis(ratios)
        return ratios
