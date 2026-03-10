from dataclasses import dataclass

@dataclass
class ReactionConditions:
    """
    Environmental parameters governing the Maillard cascade.
    """
    pH: float = 7.0
    temperature_celsius: float = 150.0 
    water_activity: float = 0.8
    fat_fraction: float = 0.0      # Added phase D
    protein_fraction: float = 0.15  # Added phase D
    dielectric_constant: float = 78.4 # Default: Water
    solvent_name: str = "water"
    
    def __post_init__(self):
        """Set dielectric constant based on solvent name if provided."""
        presets = {
            "water": 78.4,
            "ethanol": 24.5,
            "methanol": 32.7,
            "lipid": 2.0,
            "benzene": 2.3,
            "dimethyl_sulfoxide": 46.7
        }
        if self.solvent_name.lower() in presets:
            self.dielectric_constant = presets[self.solvent_name.lower()]

    @property
    def temperature_kelvin(self) -> float:
        return self.temperature_celsius + 273.15
        
    def _sigmoid(self, x: float, center: float, k: float) -> float:
        """Helper for sigmoid transitions."""
        import math
        try:
            return 1.0 / (1.0 + math.exp(-k * (x - center)))
        except OverflowError:
            return 0.0 if x < center else 1.0

    def _gaussian(self, x: float, center: float, sigma: float) -> float:
        """Helper for peaked responses (e.g. Schiff base)."""
        import math
        return math.exp(-0.5 * ((x - center) / sigma) ** 2)

    def get_ph_multiplier(self, reaction_family: str) -> float:
        """
        Calculates a kinetic multiplier based on pH using smooth sigmoid/Gaussian curves.
        
        Physics:
        - 1,2-enolisation (furans): Favored at acidic pH (pH < 6).
        - 2,3-enolisation (pyrazines): Favored at alkaline pH (pH > 6).
        - Schiff base: Optimal at pH 5.5 (amine nucleophilicity vs protonation balance).
        """
        if not reaction_family:
            return 1.0
            
        fam = reaction_family.lower()
        
        # 1. Schiff Base: Gaussian peak at pH 5.5
        if any(x in fam for x in ["schiff", "condensation"]):
            return 1.0 + 2.0 * self._gaussian(self.pH, 5.5, 1.0)

        # 2. 2,3-enolisation / Pyrazine / Strecker / Amadori / Heyns (Alkaline favored)
        elif any(x in fam for x in ["2,3", "2_3", "pyrazine", "strecker", "amadori", "heyns"]):
            base_mult = 0.2 + 8.0 * self._sigmoid(self.pH, 6.5, 1.5)
            return max(0.01, base_mult)

        # 3. 1,2-enolisation / Furan / Thiol / Thio / Cysteine / [Generic Enolisation] (Acidic favored)
        elif any(x in fam for x in ["1,2", "1_2", "furan", "thiol", "thio", "cysteine", "enolisation"]):
            return 1.0 + 4.0 * (1.0 - self._sigmoid(self.pH, 6.0, 2.0))
            
        return 1.0
        
    def get_arrhenius_multiplier(self, activation_barrier_kcal: float) -> float:
        """
        Returns the relative rate multiplier e^(-Ea / RT).
        Uses R = 1.987 cal/(mol*K).
        """
        if activation_barrier_kcal <= 0:
            return 1.0
            
        R_kcal = 0.001987
        exponent = - (activation_barrier_kcal) / (R_kcal * self.temperature_kelvin)
        
        import math
        try:
            return math.exp(exponent)
        except OverflowError:
            return 0.0
            
    def get_water_activity_multiplier(self) -> float:
        """
        Maillard reaction peaks at aw 0.6-0.8.
        Falls off rapidly due to diffusion limits (glass state) at low aw,
        and dilution at high aw.
        """
        if 0.6 <= self.water_activity <= 0.8:
            return 1.0
        elif self.water_activity < 0.6:
            # Drop off linearly to 0 at aw 0
            return max(0.01, self.water_activity / 0.6)
        else:
            # Greater than 0.8, dilute
            return max(0.1, 1.0 - (self.water_activity - 0.8) * 4)
