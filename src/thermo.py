import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
import math

# Universal gas constant in J/(mol K)
R_GAS = 8.314462618

@dataclass
class ThermoResult:
    """Contains raw and quasi-harmonic corrected thermodynamic properties."""
    temperature_k: float
    zpe_kcal_mol: float
    
    # Raw harmonic properties (Hartree)
    raw_enthalpy_h: float
    raw_entropy_h: float
    raw_gibbs_h: float
    
    # Quasi-harmonic corrected properties (Hartree)
    qh_enthalpy_h: float
    qh_entropy_h: float
    qh_gibbs_h: float
    
    # Differences (kcal/mol)
    entropy_diff_kcal_mol: float


class QuasiHarmonicCorrector:
    """
    Applies the Grimme quasi-harmonic approximation to vibrational frequencies.
    
    Standard harmonic oscillator approximations break down for low-frequency modes 
    (typically <100 cm⁻¹), overestimating the vibrational entropy contribution.
    This class corrects modes below a cutoff frequency by interpolating them towards 
    a free-rotor approximation, preventing artifactual ΔS explosions in Gibbs free energy.
    
    Reference: Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964.
    """
    
    def __init__(self, cutoff_freq_cm1: float = 100.0, temp_k: float = 298.15):
        """
        Args:
            cutoff_freq_cm1: The frequency threshold (cm⁻¹) below which corrections are applied.
                             Grimme's standard is 100 cm⁻¹.
            temp_k: Temperature in Kelvin.
        """
        self.cutoff_freq_cm1 = cutoff_freq_cm1
        self.temp_k = temp_k
        
        # Physical constants
        self.h = 6.62607015e-34      # Planck constant (J*s)
        self.k_b = 1.380649e-23      # Boltzmann constant (J/K)
        self.c = 29979245800.0       # Speed of light (cm/s)
        self.N_A = 6.02214076e23     # Avogadro constant (mol⁻¹)
        self.hartree_to_joules = 4.3597447222071e-18
        self.hartree_to_kcal = 627.509
        self.joules_mol_to_hartree = 1.0 / (self.hartree_to_joules * self.N_A)

    def _harmonic_entropy(self, freq_cm1: float) -> float:
        """Calculate standard harmonic oscillator entropy for a single mode (J / mol K)."""
        if freq_cm1 <= 0.0:
            return 0.0
            
        freq_hz = freq_cm1 * self.c
        x = (self.h * freq_hz) / (self.k_b * self.temp_k)
        
        try:
            exp_term = math.exp(-x)
            sv = self.k_b * ( x / (math.exp(x) - 1.0) - math.log(1.0 - exp_term) )
            return sv * self.N_A
        except OverflowError:
            return 0.0

    def _rotor_entropy(self, freq_cm1: float, moment_of_inertia: float = 1.0e-44) -> float:
        """
        Calculate free-rotor entropy for a single mode (J / mol K).
        Standard approach assumes an average moment of inertia (Bav) for all low modes.
        Grimme uses Bav = 1e-44 kg*m^2.
        """
        # Bav = h / (8 * pi^2 * c * v_0)
        # We simplify to the established damping function by Grimme
        if freq_cm1 <= 0.0:
            return 0.0
            
        freq_hz = freq_cm1 * self.c
        mu = (self.h * freq_hz) / (self.k_b * self.temp_k)
        
        # Free rotor approximation
        sr = self.k_b * (0.5 + 0.5 * math.log( (8.0 * math.pi**3 * moment_of_inertia * self.k_b * self.temp_k) / (self.h**2) ))
        
        return sr * self.N_A

    def calculate_thermo(self, freqs_cm1: List[float], electronic_energy_h: float = 0.0) -> ThermoResult:
        """
        Calculate thermodynamic properties given a list of vibrational frequencies.
        
        Args:
            freqs_cm1: List of vibrational frequencies in cm⁻¹. Should exclude imaginary/negative freqs.
            electronic_energy_h: The electronic energy in Hartrees (used as the base for Enthalpy/Gibbs).
        
        Returns:
            ThermoResult object containing all thermodynamic properties.
        """
        # Constants
        R = R_GAS  # J/(mol K)
        
        real_freqs = [float(np.real(f)) for f in freqs_cm1 if np.isreal(f) and float(np.real(f)) > 0.01]
        
        # Zero Point Energy (J/mol)
        zpe_j_mol = sum(0.5 * self.h * (f * self.c) * self.N_A for f in real_freqs)
        zpe_hartree = zpe_j_mol * self.joules_mol_to_hartree
        zpe_kcal = zpe_hartree * self.hartree_to_kcal
        
        # Vibrational Enthalpy limit: H_vib = ZPE + \int_0^T C_v dT
        # Standard harmonic approximation for enthalpy is usually fine regardless of low freqs
        h_vib_term = sum( ( (self.h * f * self.c * self.N_A) / (math.exp((self.h * f * self.c)/(self.k_b * self.temp_k)) - 1.0) ) for f in real_freqs )
        
        thermal_energy_j_mol = zpe_j_mol + h_vib_term
        
        # H = E_elec + E_thm + RT (assuming ideal gas behavior)
        enthalpy_h = electronic_energy_h + (thermal_energy_j_mol * self.joules_mol_to_hartree) + (R * self.temp_k * self.joules_mol_to_hartree)
        
        # Calculate Entropies
        raw_S_vib_j_mol_k = 0.0
        qh_S_vib_j_mol_k = 0.0
        
        for f in real_freqs:
            s_harm = self._harmonic_entropy(f)
            raw_S_vib_j_mol_k += s_harm
            
            # Grimme damping: alpha = 1 / (1 + (cutoff / f)^4 )
            alpha = 1.0 / (1.0 + (self.cutoff_freq_cm1 / f)**4)
            s_rot = self._rotor_entropy(f)
            
            # Interpolated entropy
            s_qh = (alpha * s_harm) + ((1.0 - alpha) * s_rot)
            qh_S_vib_j_mol_k += s_qh
            
        # Add translational and rotational entropy (assumed constant between standard and QH for simplicity here, 
        # normally retrieved from full thermo outputs, but we isolate vibrational correction here)
        # Since we just want the *difference* or the corrected Gibbs, the correction is solely in the vibrational entropy term.
        
        # We compute the *delta* S from harmonic to QH, and apply it.
        delta_S_j_mol_k = qh_S_vib_j_mol_k - raw_S_vib_j_mol_k
        
        # Return partial properties just based on vibration/ZPE (relative Gibbs, not absolute, unless translational/rotational added)
        # Note: A full implementation integrates seamlessly with PySCF's thermo module by replacing its S_vib.
        # Here we provide the correction delta which can be applied directly to a raw `Thermodynamics` output.
        
        raw_entropy_h = raw_S_vib_j_mol_k * self.joules_mol_to_hartree
        qh_entropy_h = qh_S_vib_j_mol_k * self.joules_mol_to_hartree
        
        raw_gibbs_h = enthalpy_h - (self.temp_k * raw_entropy_h)
        qh_gibbs_h = enthalpy_h - (self.temp_k * qh_entropy_h)
        
        entropy_diff_kcal = (qh_entropy_h - raw_entropy_h) * self.temp_k * self.hartree_to_kcal
        
        return ThermoResult(
            temperature_k=self.temp_k,
            zpe_kcal_mol=zpe_kcal,
            raw_enthalpy_h=enthalpy_h,
            raw_entropy_h=raw_entropy_h,
            raw_gibbs_h=raw_gibbs_h,
            qh_enthalpy_h=enthalpy_h,
            qh_entropy_h=qh_entropy_h,
            qh_gibbs_h=qh_gibbs_h,
            entropy_diff_kcal_mol=entropy_diff_kcal
        )
