"""
src/kinetics.py

Transition State Theory (TST) and Microkinetic modeling utilities.
Converts DFT barriers (Delta G‡) into rate constants and temporal fluxes.
"""

import numpy as np
from scipy.constants import kilo, calorie_th, Planck, Boltzmann, gas_constant

class KineticsEngine:
    def __init__(self, temperature_k: float = 423.15):
        self.T = temperature_k
        # Convert kcal/mol to Joules per molecule
        self.kcal_to_j_per_mol = calorie_th * kilo
        self.kb = Boltzmann
        self.h = Planck
        self.R = gas_constant

    def get_rate_constant(self, delta_g_kcal_mol: float) -> float:
        """
        Calculate the first-order rate constant k using Eyring-Polanyi equation:
        k = (kB * T / h) * exp(-Delta G‡ / RT)
        """
        delta_g_j = delta_g_kcal_mol * self.kcal_to_j_per_mol
        
        pre_exponential = (self.kb * self.T) / self.h
        exponential = np.exp(-delta_g_j / (self.R * self.T))
        
        return pre_exponential * exponential

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

    def simulate_network_cantera(self, mechanism_yaml: str, initial_concentrations: Dict[str, float], 
                                 time_span: tuple, temperature_k: Optional[float] = None) -> Dict[str, np.ndarray]:
        """
        Phase 12: Rigorous ODE integration of a reaction network using Cantera.
        
        Args:
            mechanism_yaml: Path to the Cantera YAML mechanism file.
            initial_concentrations: Dict mapping species names to initial molarities.
            time_span: Tuple (t_start, t_end) in seconds.
            temperature_k: Temperature for the simulation (defaults to self.T).
            
        Returns:
            Dict mapping species names to concentration arrays.
        """
        try:
            import cantera as ct
        except ImportError:
            raise ImportError("Cantera is not installed. Please run 'pip install cantera'.")

        T = temperature_k or self.T
        
        # 1. Load the mechanism
        gas = ct.Solution(mechanism_yaml)
        gas.TP = T, 101325  # Standard pressure
        
        # 2. Set initial state
        # Cantera typically uses mass/mole fractions. For simplicity in Maillard 
        # liquid phases, we approximate molarity as relative mole fractions if sum is small.
        gas.X = initial_concentrations
        
        # 3. Create a batch reactor
        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])
        
        # 4. Integrate
        time_points = np.linspace(time_span[0], time_span[1], 100)
        results = {name: [] for name in gas.species_names}
        results["time"] = []
        
        for t in time_points:
            sim.advance(t)
            results["time"].append(t)
            for i, name in enumerate(gas.species_names):
                results[name].append(gas.concentrations[i])
                
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
