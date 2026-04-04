from dataclasses import dataclass
from typing import Optional, Sequence

from src.barrier_constants import arrhenius_rate_constant, get_donor_reactivity_multiplier
from src.extrusion import build_extrusion_process_profile, normalize_moisture_regime

@dataclass
class ReactionConditions:
    """
    Environmental parameters governing the Maillard cascade.
    """
    def __init__(self,
                 pH: float = 6.0,
                 temperature_celsius: float = 120.0,
                 water_activity: float = 0.8,
                 fat_fraction: float = 0.0,
                 protein_fraction: float = 1.0,
                 dielectric_constant: float = 78.4, # Default: Water
                 solvent_name: str = "water",
                 matrix_fiber: float = 0.0, # Placeholder for blind spot
                 metal_catalyst: Optional[str] = None, # Placeholder for blind spot
                 protein_type: str = "free",
                 sme_kj_per_kg: float = 0.0,
                 moisture_regime: Optional[str] = None,
                 sterilization_temperature_celsius: Optional[float] = None,
                 sterilization_time_minutes: float = 0.0,
                 barrel_zone_temperatures: Optional[Sequence[float]] = None,
                 barrel_zone_time_fractions: Optional[Sequence[float]] = None,
                 ):
        self.pH = pH
        self.temperature_celsius = temperature_celsius
        self.water_activity = water_activity
        self.fat_fraction = fat_fraction
        self.protein_fraction = protein_fraction
        self.dielectric_constant = dielectric_constant
        self.solvent_name = solvent_name
        self.matrix_fiber = matrix_fiber
        self.metal_catalyst = metal_catalyst
        self.protein_type = protein_type
        self.sme_kj_per_kg = float(sme_kj_per_kg)
        self.moisture_regime = moisture_regime
        self.sterilization_temperature_celsius = sterilization_temperature_celsius
        self.sterilization_time_minutes = float(sterilization_time_minutes)
        self.barrel_zone_temperatures = None if barrel_zone_temperatures is None else [float(value) for value in barrel_zone_temperatures]
        self.barrel_zone_time_fractions = None if barrel_zone_time_fractions is None else [float(value) for value in barrel_zone_time_fractions]
        self.__post_init__()

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
        self.moisture_regime = normalize_moisture_regime(self.moisture_regime, self.water_activity)

    @property
    def extrusion_profile(self) -> dict[str, object]:
        return build_extrusion_process_profile(
            base_temperature_celsius=self.temperature_celsius,
            water_activity=self.water_activity,
            protein_type=self.protein_type,
            sme_kj_per_kg=self.sme_kj_per_kg,
            moisture_regime=self.moisture_regime,
            sterilization_temperature_celsius=self.sterilization_temperature_celsius,
            sterilization_time_minutes=self.sterilization_time_minutes,
            zone_temperatures=self.barrel_zone_temperatures,
            zone_time_fractions=self.barrel_zone_time_fractions,
        )

    @property
    def effective_temperature_celsius(self) -> float:
        return float(self.extrusion_profile.get("effective_temperature_celsius", self.temperature_celsius))

    @property
    def temperature_kelvin(self) -> float:
        return self.effective_temperature_celsius + 273.15

    @property
    def nominal_temperature_kelvin(self) -> float:
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
        mult = 1.0
        if any(x in fam for x in ["schiff", "condensation"]):
            mult *= (1.0 + 2.0 * self._gaussian(self.pH, 5.5, 1.0))

        # 2. 2,3-enolisation / Pyrazine / Strecker / Amadori / Heyns (Alkaline favored)
        elif any(x in fam for x in ["2,3", "2_3", "pyrazine", "strecker", "amadori", "heyns", "nitrogen_heterocycle"]):
            base_mult = 0.2 + 8.0 * self._sigmoid(self.pH, 6.5, 1.5)
            mult *= max(0.01, base_mult)

        # 3. 1,2-enolisation / Furan / Thiol / Thio / Cysteine / [Generic Enolisation] (Acidic favored)
        elif any(x in fam for x in ["1,2", "1_2", "furan", "thiol", "thio", "cysteine", "enolisation", "oxygen_heterocycle"]):
            mult *= (1.0 + 4.0 * (1.0 - self._sigmoid(self.pH, 6.0, 2.0)))
            
        # 4. Metal Catalysis (Heme/Iron Enrichment)
        # Heme is a potent promoter of lipid oxidation and certain Maillard steps.
        if self.metal_catalyst and self.metal_catalyst.lower() == "heme":
            if "oxidation" in fam or "radical" in fam:
                mult *= 5.0  # Heme significantly accelerates radical generation
            elif "pyrazine" in fam:
                mult *= 1.5  # Iron catalysis of nitrogen cyclization
                
        return mult

    def get_rate_constant(
        self,
        pathway_type: str,
        ea_override_kcal: float = None,
        *,
        reactant_labels: Optional[Sequence[str]] = None,
    ) -> float:
        """
        Arrhenius rate constant: k = A * exp(-Ea / RT)
        ea_override_kcal: barrier in kcal/mol from results_db.py
        
        Default Activation Energies (kJ/mol):
        - amadori_formation: 75.0
        - strecker_degradation: 85.0
        - thiol_formation: 95.0
        - pyrazine_formation: 110.0
        - furfural_formation: 70.0
        - acrylamide_formation: 130.0
        """
        T_K = self.temperature_kelvin
        
        # Ea in kJ/mol
        ACTIVATION_ENERGIES_KJ = {
            "amadori_formation": 75.0,
            "strecker_degradation": 85.0,
            "thiol_formation": 95.0,       
            "pyrazine_formation": 110.0,
            "furfural_formation": 70.0,
            "acrylamide_formation": 130.0,
        }
        
        if ea_override_kcal is not None:
            Ea_J = ea_override_kcal * 4184.0 # kcal/mol -> J/mol
        else:
            # Match family to default kJ/mol map
            fam = pathway_type.lower()
            ea_kj = 90.0 # Default fallback
            for key, val in ACTIVATION_ENERGIES_KJ.items():
                if key in fam or fam in key:
                    ea_kj = val
                    break
            Ea_J = ea_kj * 1000.0
        
        barrier_kcal = Ea_J / 4184.0
        
        # Apply pH ionization correction
        ph_factor = self._ionization_correction(pathway_type)
        
        # Apply water activity correction (Labuza)
        aw_factor = self._water_activity_correction(pathway_type)
        donor_factor = get_donor_reactivity_multiplier(pathway_type, reactant_labels=reactant_labels)
        
        return arrhenius_rate_constant(
            barrier_kcal,
            T_K,
            family=pathway_type,
            multiplier=ph_factor * aw_factor * donor_factor,
        )

    def _ionization_correction(self, pathway_type: str) -> float:
        """
        Henderson-Hasselbalch based correction for reactive species ionization.
        Replaces arbitrary sigmoids.
        """
        fam = pathway_type.lower()
        # Amine reactions require the deprotonated form (R-NH2)
        if any(x in fam for x in ("amadori", "strecker", "schiff")):
            pKa = 8.0  # Common alpha-amino pKa
            return 1.0 / (1.0 + 10**(pKa - self.pH))
        
        # Pyrazines often peak at slightly alkaline pH
        elif "pyrazine" in fam:
            return 1.0 / (1.0 + 10**(6.5 - self.pH))
            
        return 1.0

    def _water_activity_correction(self, pathway_type: str) -> float:
        """
        Correction for water activity effect on Maillard reaction rate.
        Based on Labuza & Saltmarch (1981).
        Rate peaks around aw=0.65, drops at both extremes.
        """
        aw = self.water_activity
        fam = pathway_type.lower()
        
        if any(x in fam for x in ("amadori", "strecker", "pyrazine")):
            if aw < 0.2:
                return 0.1
            elif aw < 0.65:
                # Linear ramp
                return 0.1 + (aw - 0.2) / (0.65 - 0.2) * 0.9
            else:
                # Product inhibition/dilution beyond peak
                return 1.0 - (aw - 0.65) / (1.0 - 0.65) * 0.4
                
        elif "furfural" in fam:
            # Furfural formation (acid-catalyzed dehydration) favored at lower aw
            return max(0.1, 1.3 - aw)
            
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
        Calculates the kinetic multiplier for water activity (aw).
        
        Physics:
        Maillard reactions peak at aw ≈ 0.6-0.8 (typically around 0.7).
        - At low aw (glassy state), diffusion becomes limiting.
        - At high aw, water acts as a diluent and product inhibitor (Le Chatelier).
        
        Replaced the piecewise linear model with a continuous Gaussian curve:
        mult = exp(-0.5 * ((aw - 0.7) / 0.15)^2)
        """
        # Center peak around aw = 0.7, with a standard deviation of 0.15
        val = self._gaussian(self.water_activity, center=0.7, sigma=0.15)
        # Prevent it from dropping exactly to zero to avoid numerical singularities in ODEs
        return max(0.01, val)
