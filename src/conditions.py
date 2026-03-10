from dataclasses import dataclass

@dataclass
class ReactionConditions:
    """
    Environmental parameters governing the Maillard cascade.
    """
    pH: float = 7.0
    temperature_celsius: float = 150.0 
    water_activity: float = 0.8
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
        
    def get_ph_multiplier(self, reaction_family: str) -> float:
        """
        Calculates a kinetic multiplier based on pH according to known Maillard chemistry.
        - Acidic pH (< 6) favors 1,2-enolisation (furan pathways).
        - Alkaline/Neutral pH (>= 6) favors 2,3-enolisation (pyrazine/Strecker pathways).
        """
        if not reaction_family:
            return 1.0
            
        fam = reaction_family.lower()
        
        # 1,2-enolisation to furans (favored at acidic pH)
        if "1,2" in fam or "1_2" in fam or "furan" in fam or "thiol" in fam:
            if self.pH < 6.0:
                return 5.0 # Accelerated
            return 1.0
            
        # 2,3-enolisation to pyrazines/strecker (favored at alkaline pH)
        if "2,3" in fam or "2_3" in fam or "pyrazine" in fam or "strecker" in fam:
            if self.pH >= 6.0:
                base = 1.0
                # Exponential increase with alkalinity
                base += (self.pH - 6.0) * 2.0 
                return max(1.0, base)
            return 0.2 # Suppressed at low pH
            
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
