"""
src/kinetics.py

Transition State Theory (TST) and Microkinetic modeling utilities.
Converts DFT barriers (Delta G‡) into rate constants and temporal fluxes.
"""

import numpy as np
from scipy.constants import kilo, calorie_th, Planck, Boltzmann, gas_constant
from typing import Dict, List, Optional, Tuple, Set
from src.conditions import ReactionConditions
from src.thermo import JobackEstimator

class KineticsEngine:
    def __init__(self, temperature_k: float = 423.15):
        self.T = temperature_k
        # Convert kcal/mol to Joules per molecule
        self.kcal_to_j_per_mol = calorie_th * kilo
        self.kb = Boltzmann
        self.h = Planck
        self.R = gas_constant

    def get_rate_constant(self, delta_g_kcal_mol: float, 
                          temperature_k: Optional[float] = None,
                          conditions: Optional[ReactionConditions] = None,
                          reaction_family: Optional[str] = None) -> float:
        """
        Calculate the first-order rate constant k using Eyring-Polanyi equation:
        k = (kB * T / h) * exp(-Delta G‡ / RT)
        
        Includes Phase 12 enhancements:
        - pH scaling from ReactionConditions
        - Kirkwood-Onsager Solvent Scaling
        """
        T = temperature_k or self.T
        if conditions:
            T = conditions.temperature_kelvin
            
        # 1. Apply pH/Environmental Multipliers
        multiplier = 1.0
        if conditions and reaction_family:
            multiplier *= conditions.get_ph_multiplier(reaction_family)
            multiplier *= conditions.get_water_activity_multiplier()
            
        # 2. Kirkwood-Onsager Solvent Scaling (Phase 12.1)
        # Simplified: Polar transitions (Maillard) are accelerated in polar solvents.
        # Barriers are adjusted by f(epsilon) = (eps-1)/(2eps+1).
        # We assume the input barrier is for water (eps=78.4, f=0.49).
        # Adjusted Barrier = Base_Barrier - Shift * (f(eps) - f(water))
        if conditions:
            eps = conditions.dielectric_constant
            f_eps = (eps - 1) / (2 * eps + 1)
            f_water = (78.4 - 1) / (2 * 78.4 + 1)
            # Empirical shift: 5 kcal/mol sensitivity to solvent polarity
            barrier_shift = 5.0 * (f_eps - f_water)
            delta_g_kcal_mol -= barrier_shift

        delta_g_j = delta_g_kcal_mol * self.kcal_to_j_per_mol
        
        pre_exponential = (self.kb * T) / self.h
        exponential = np.exp(-delta_g_j / (self.R * T))
        
        return pre_exponential * exponential * multiplier

    def simulate_flux(self, initial_conc: float, rate_constant: float, time_steps_sec: np.ndarray) -> np.ndarray:
        """
        Simple first-order decay simulation: [A](t) = [A]0 * exp(-kt)
        """
        return initial_conc * np.exp(-rate_constant * time_steps_sec)

    def simulate_competition(self, initial_conc: float, rates_dict: dict, time_steps_sec: np.ndarray) -> dict:
        """
        Simulate a branching reaction: A -> B (k1), A -> C (k2), etc.
        Returns a dictionary of concentrations over time.
        """
        k_total = sum(rates_dict.values())
        conc_a = initial_conc * np.exp(-k_total * time_steps_sec)
        
        results = {"A": conc_a}
        for product, k in rates_dict.items():
            # [Product](t) = [A]0 * (k / k_total) * (1 - exp(-k_total * t))
            conc_p = initial_conc * (k / k_total) * (1 - np.exp(-k_total * time_steps_sec))
            results[product] = conc_p
            
        return results

    def is_reaction_feasible(self, reactants: List[str], products: List[str], 
                             threshold_kcal_mol: float = 30.0) -> Tuple[bool, float]:
        """
        Phase 12.3: Dynamic Thermo-Gating.
        Uses Joback Group Additivity to estimate Delta G for a reaction.
        If Delta G > threshold, the reaction is considered unphysical (non-spontaneous).
        
        Returns:
            (is_feasible, delta_g_estimated_kcal_mol)
        """
        try:
            r_energy = sum(JobackEstimator.estimate(s)["G298"] for s in reactants)
            p_energy = sum(JobackEstimator.estimate(s)["G298"] for s in products)
            
            # Convert J/mol to kcal/mol
            delta_g = (p_energy - r_energy) / 4184.0
            
            return delta_g <= threshold_kcal_mol, delta_g
        except Exception as e:
            # If Joback fails (unsupported groups), default to feasible for safety
            return True, 0.0

    def simulate_network_cantera(self, mechanism_yaml: str, initial_concentrations: Dict[str, float], 
                                 time_span: tuple, temperature_k: Optional[float] = None,
                                 temperature_profile: Optional[List[Tuple[float, float]]] = None) -> Dict[str, np.ndarray]:
        """
        Phase 12/15: Rigorous ODE integration of a reaction network using Cantera.
        Supports optional non-isothermal temperature profiles.
        
        Args:
            mechanism_yaml: Path to the Cantera YAML mechanism file.
            initial_concentrations: Dict mapping species names to initial molarities.
            time_span: Tuple (t_start, t_end) in seconds.
            temperature_k: Constant temperature (ignored if temperature_profile is provided).
            temperature_profile: List of (time_sec, temp_k) points for interpolation.
            
        Returns:
            Dict mapping species names to concentration arrays.
        """
        try:
            import cantera as ct
        except ImportError:
            raise ImportError("Cantera is not installed. Please run 'pip install cantera'.")

        # 1. Load the mechanism
        gas = ct.Solution(mechanism_yaml)
        
        # Initial T
        if temperature_profile:
            t_pts, t_vals = zip(*temperature_profile)
            T_init = float(np.interp(time_span[0], t_pts, t_vals))
        else:
            T_init = temperature_k or self.T
            
        gas.TP = T_init, 101325  # Standard pressure
        
        # 2. Set initial state
        gas.X = initial_concentrations
        
        # 3. Create a batch reactor
        # We use an IdealGasReactor but we will manually update its temperature
        # at each step to simulate a defined T(t) profile.
        # clone=False suppresses Cantera 3.2+ deprecation warnings
        r = ct.IdealGasConstPressureReactor(gas, clone=False)
        r.energy_enabled = False
        sim = ct.ReactorNet([r])
        
        # 4. Integrate
        n_points = 100
        time_points = np.linspace(time_span[0], time_span[1], n_points)
        results = {name: [] for name in gas.species_names}
        # Also store mole fractions for T/P invariant conversion math
        for name in gas.species_names:
            results[f"{name}_X"] = []
            
        results["time"] = []
        results["temperature"] = []
        
        for t in time_points:
            # Update temperature if profile provided
            if temperature_profile:
                T_curr = float(np.interp(t, t_pts, t_vals))
                gas.TP = T_curr, 101325
                r.syncState() 
            
            sim.advance(t)
            
            results["time"].append(t)
            results["temperature"].append(gas.T)
            for i, name in enumerate(gas.species_names):
                results[name].append(gas.concentrations[i])
                results[f"{name}_X"].append(gas.X[i])
                
        # Convert to numpy arrays
        return {k: np.array(v) for k, v in results.items()}

if __name__ == "__main__":
    # Quick sanity check for a 20 kcal/mol barrier at 150 C (423 K)
    engine = KineticsEngine(temperature_k=423.15)
    k = engine.get_rate_constant(20.0)
    print(f"Rate constant for 20 kcal/mol at 150C: {k:.2e} s^-1")
    
    # 50% conversion time (t1/2 = ln(2)/k)
    t_half_min = (np.log(2) / k) / 60
    print(f"Half-life: {t_half_min:.2f} minutes")

if __name__ == "__main__":
    # Quick sanity check for a 20 kcal/mol barrier at 150 C (423 K)
    engine = KineticsEngine(temperature_k=423.15)
    k = engine.get_rate_constant(20.0)
    print(f"Rate constant for 20 kcal/mol at 150C: {k:.2e} s^-1")
    
    # 50% conversion time (t1/2 = ln(2)/k)
    t_half_min = (np.log(2) / k) / 60
    print(f"Half-life: {t_half_min:.2f} minutes")
