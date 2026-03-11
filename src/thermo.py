import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any
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
from rdkit import Chem
from scipy.constants import gas_constant

# Joback Group Contributions
# Format: { 'smarts': (dfH, dfG, cp_a, cp_b, cp_c, cp_d) }
# Units: dfH, dfG in kJ/mol. Cp = a + bT + cT^2 + dT^3 in J/(mol K)
JOBACK_GROUPS = {
    "carboxyl": ("[CX3](=O)[OX2H1]", -426.72, -387.87, 2.41e1, 4.27e-1, -2.88e-4, 7.40e-8),
    "aldehyde": ("[CX3H1]=O", -162.03, -143.48, 3.09e1, -3.36e-2, 1.60e-4, -8.67e-8),
    "ketone": ("[CX3](=O)[#6]", -132.18, -120.48, 1.35e1, 2.20e-1, -1.14e-4, 2.50e-8),
    "hydroxyl": ("O[H]", -208.04, -189.20, 2.57e1, -6.91e-2, 1.77e-4, -9.88e-8),
    "ether": ("[OX2]([#6])[#6]", -132.22, -105.00, 2.55e1, -6.32e-2, 1.11e-4, -5.48e-8),
    "primary_amine": ("[NH2;X1]", -22.02, 14.07, 2.69e1, -4.12e-2, 1.64e-4, -9.76e-8),
    "secondary_amine": ("[NH;X2]", 53.47, 89.04, -1.21, 2.33e-1, -1.74e-4, 4.66e-8),
    "thiol": ("[SH;X1]", -11.33, 8.44, 3.53e1, -7.58e-2, 1.85e-4, -1.03e-7),
    "sulfide": ("[SX2]", 68.07, 80.24, 3.43e1, -1.28e-2, 1.86e-4, -1.01e-7),
    "quaternary_c": ("[C;X4]", 8.25, 20.97, -3.74e1, 1.30e0, -1.02e-3, 2.69e-7),
    "methine": ("[CH;X3]", -6.12, 7.93, -2.30e1, 1.10e0, -7.26e-4, 1.76e-7),
    "methylene": ("[CH2;X2]", -20.64, -8.42, -9.09, 9.50e-1, -5.44e-4, 1.19e-7),
    "methyl": ("[CH3;X1]", -45.83, -43.85, 1.95e1, 8.08e-1, -4.60e-4, 9.67e-8),
    "alkene_methylene": ("[CH2;X2]=C", -1.37, 15.05, 2.36e1, -3.81e-2, 1.72e-4, -1.03e-7),
    "alkene_methine": ("[CH;X2]=C", 8.64, 33.25, -1.61, 2.35e-1, -9.85e-5, 1.53e-8),
    "ring_methylene": ("[CH2;X2;R]", -4.82, 11.37, -6.03, 8.54e-1, -4.80e-4, 1.05e-7),
    "ring_methine": ("[CH;X3;R]", 7.26, 26.15, -2.05e1, 9.71e-1, -5.98e-4, 1.45e-7),
}

class JobackEstimator:
    @staticmethod
    def estimate(smiles: str) -> Dict[str, any]:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        
        dfH = 68.29 # Base kJ/mol
        dfG = 53.88 # Base kJ/mol
        cp_coeffs = np.zeros(4)
        
        matched_indices = set()
        priority_keys = list(JOBACK_GROUPS.keys())
        
        for k in priority_keys:
            smarts, h, g, a, b, c, d = JOBACK_GROUPS[k]
            patt = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(patt)
            
            for m in matches:
                if m[0] not in matched_indices:
                    matched_indices.add(m[0])
                    dfH += h
                    dfG += g
                    cp_coeffs[0] += a
                    cp_coeffs[1] += b
                    cp_coeffs[2] += c
                    cp_coeffs[3] += d
        
        # Apply Joback polynomial offsets
        cp_coeffs[0] -= 37.93
        cp_coeffs[1] += 0.210
        cp_coeffs[2] -= 3.91e-4
        cp_coeffs[3] += 2.06e-7
        
        return {
            "H298": dfH * 1000.0, # J/mol
            "G298": dfG * 1000.0, # J/mol
            "cp_coeffs": cp_coeffs
        }

# Literature NASA-7 polynomials for common small molecules (low T range: 200-1000K)
# Source: Burcat & Ruscic (2005) or NIST
SMALL_MOLECULE_THERMO = {
    "O": [4.19864056, -2.03643410e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12, -3.02937267e4, -0.84903220], # H2O
    "O=C=O": [2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699544e-13, -4.83719697e4, 9.90105222], # CO2
    "N": [2.73147413, 5.92211918e-3, -8.62541814e-7, -1.50343719e-9, 8.01021461e-13, -6.69033346e3, 4.61741644],  # NH3
    "[HH]": [2.34433112, 7.98052075e-3, -1.94781510e-5, 2.01572094e-8, -7.37611761e-12, -8.54730620e2, -3.94046088], # H2
    "S": [3.32730594, 3.25049333e-3, 4.35401667e-7, -2.13203333e-9, 8.44026667e-13, -3.42000000e3, 7.32000000],  # H2S
}

def get_nasa_coefficients(smiles: str) -> List[float]:
    """
    Returns 14 NASA polynomial coefficients for the given SMILES (two ranges).
    Low range: 300-1000 K. High range: 1000-3000 K.
    Returns: [low_a1..a7, high_a1..a7]
    """
    # 1. Specialized Literature Overrides
    if smiles in SMALL_MOLECULE_THERMO:
        # Standard NASA-7 for small molecules usually has 14 or 15 coeffs. 
        # We'll duplicate the single range we have into both for these.
        c = SMALL_MOLECULE_THERMO[smiles]
        return c + c
        
    if smiles == "[H][H]": smiles = "[HH]" # Normalization
    
    # 2. Joback Estimation (Low Range: 300-1000K)
    res = JobackEstimator.estimate(smiles)
    R = gas_constant
    
    T = np.linspace(300, 1000, 20)
    a, b, c, d = res["cp_coeffs"]
    Cp = a + b*T + c*T**2 + d*T**3 
    
    # Cp/R = a1 + a2T + a3T^2 + a4T^3 + a5T^4
    low_coeffs = np.polyfit(T, Cp/R, 4)[::-1]
    a1, a2, a3, a4, a5 = low_coeffs
    
    # a6 (Enthalpy offset)
    T0 = 298.15
    H0 = res["H298"]
    a6 = (H0 / R) - (a1*T0 + a2*T0**2/2 + a3*T0**3/3 + a4*T0**4/4 + a5*T0**5/5)
    
    # a7 (Entropy offset)
    S0 = (res["H298"] - res["G298"]) / T0
    a7 = (S0 / R) - (a1*np.log(T0) + a2*T0 + a3*T0**2/2 + a4*T0**3/3 + a5*T0**4/4)
    low_range = [float(x) for x in [a1, a2, a3, a4, a5, a6, a7]]
    
    # 3. Continuous High Range (1000-3000K)
    # We maintain the 1000K values but zero out the slopes to avoid wild extrapolation.
    # To ensure continuity in H and S at Tmid (1000K), we calculate high-range a6 and a7.
    T_mid = 1000.0
    
    # Low-range Cp/R, H/RT, S/R at T_mid
    a1, a2, a3, a4, a5, a6_low, a7_low = low_range
    cp_r_mid = a1 + a2*T_mid + a3*T_mid**2 + a4*T_mid**3 + a5*T_mid**4
    h_rt_mid = a1 + a2*T_mid/2 + a3*T_mid**2/3 + a4*T_mid**3/4 + a5*T_mid**4/5 + a6_low/T_mid
    s_r_mid = a1*np.log(T_mid) + a2*T_mid + a3*T_mid**2/2 + a4*T_mid**3/3 + a5*T_mid**4/4 + a7_low
    
    # High range a1 set to match Cp/R at Tmid
    high_a1 = cp_r_mid
    high_6 = (h_rt_mid * T_mid) - (high_a1 * T_mid)
    high_7 = s_r_mid - (high_a1 * np.log(T_mid))
    
    high_range = [float(high_a1), 0.0, 0.0, 0.0, 0.0, float(high_6), float(high_7)]
    
    return low_range + high_range
