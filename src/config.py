"""
src/config.py

Centralized default configurations for the Maillard framework.
Prevents magic constants from leaking across scripts and modules.
"""

from dataclasses import dataclass

@dataclass(frozen=True)
class PipelineDefaults:
    # Physical Conditions
    default_ph: float = 6.0
    default_temp_c: float = 150.0
    default_aw: float = 1.0
    default_time_minutes: float = 60.0
    default_protein_type: str = "free"
    default_denaturation_state: float = 0.5

    # Inverse Design / Recommendations
    default_target_tag: str = "meaty"
    default_minimize_tag: str = "beany"
    
    # Microkinetics
    default_cantera_time_sec: float = 600.0

DEFAULTS = PipelineDefaults()
