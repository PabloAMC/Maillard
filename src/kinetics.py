"""
src/kinetics.py

Transition State Theory (TST) and Microkinetic modeling utilities.
Converts DFT barriers (Delta G‡) into rate constants and temporal fluxes.
"""

import numpy as np
from scipy.constants import kilo, calorie_th, Planck, Boltzmann, gas_constant
from typing import Dict, List, Optional, Tuple, Sequence
from src.logger import get_logger
logger = get_logger(__name__)

from src.conditions import ReactionConditions  # noqa: E402
from src.barrier_constants import get_donor_reactivity_multiplier  # noqa: E402
from src.thermo import JobackEstimator  # noqa: E402

# ---------------------------------------------------------------------------
# Thermodynamic basis bookkeeping (audit 2026-08-26, item 3.3)
# ---------------------------------------------------------------------------
# Joback returns *formation* properties at 298.15 K: dfH298 and dfG298, i.e.
# both are referenced to the constituent elements in their standard states.
# Therefore
#     s298 := (dfH298 - dfG298) / 298.15  ==  dfS298  ==  S(compound) - sum_e n_e S(element_e)
# is likewise a formation-basis entropy, NOT the compound's absolute entropy.
#
# `get_reaction_thermo` builds per-species
#     h_i(T) = dfH298_i + int_298^T Cp_i dT          (= H_i(T) - sum_e n_ei H_e(298))
#     s_i(T) = dfS298_i + int_298^T Cp_i/T' dT'      (= S_i(T) - sum_e n_ei S_e(298))
# so each species carries a *constant* element offset frozen at 298.15 K.
# Summing over a reaction with stoichiometric coefficients nu_i:
#     sum_i nu_i h_i(T) = dH_rxn(T) - sum_e (delta n_e) H_e(298)
#     sum_i nu_i s_i(T) = dS_rxn(T) - sum_e (delta n_e) S_e(298)
# The element offsets cancel EXACTLY when the reaction is atom balanced
# (delta n_e = 0 for every element), so dG_rxn(T) = dH_rxn(T) - T dS_rxn(T) is
# correct on a single, consistent basis at every temperature. The audit's
# "mismatched bases" concern is therefore a false positive for balanced steps
# (verified numerically in tests/unit/test_thermo_basis_and_guards.py).
#
# What is NOT safe is an atom-UNBALANCED step: there the leftover element terms
# survive, and they are large (a silently dropped H2O leaves roughly
# +58 kcal/mol of dfH residual plus a temperature-proportional element-entropy
# residual). Such a dG is not a property of the chemistry at all, so it is
# labelled explicitly in the returned payload.
#
# NEUTRALIZE_UNBALANCED_THERMO_GATE controls what happens to that number.
#   False (default, current shipped behaviour): the raw element-contaminated
#     value is returned unchanged. Downstream gates (benchmark_validation,
#     cantera_export) then usually kill unbalanced steps,
#     which is accidentally protective — tests/unit/test_advanced_kinetics.py
#     ::test_thermo_gating relies on exactly this to reject mass-creating
#     reactions (CH4 -> decane).
#   True: unbalanced steps return dG = 0 (gate-neutral) and rely on an explicit
#     mass-balance check upstream instead. Flipping this is a science-owner
#     decision, because it removes the accidental mass-creation guard.
# Measured exposure on the current panel: 197/198 enumerated SMIRKS steps are
# atom balanced; the single unbalanced step
# (thiamine Additive_Thermal_Degradation, dG = 17.4 kcal/mol) is below the
# 30 kcal/mol gate, and thermodynamic gating currently resolves to "off" for
# every benchmark in data/benchmarks. Impact of flipping today: zero.
NEUTRALIZE_UNBALANCED_THERMO_GATE = False

# Standard-state absolute entropies of the elements at 298.15 K, J/(mol K).
# Source: CODATA / NIST-JANAF. Used only to expose the element reference that
# the formation basis carries, and by the unit tests that assert the basis.
ELEMENT_STANDARD_ENTROPY_J_MOL_K = {
    "C": 5.74,      # graphite
    "H": 130.68 / 2.0,   # 1/2 S(H2, g)
    "O": 205.15 / 2.0,   # 1/2 S(O2, g)
    "N": 191.61 / 2.0,   # 1/2 S(N2, g)
    "S": 32.05,     # rhombic
    "P": 41.09,     # white
    "Cl": 223.08 / 2.0,
}


def reaction_element_balance(
    reactants: Sequence[str], products: Sequence[str]
) -> Dict[str, object]:
    """
    Element/charge balance of a reaction written as SMILES lists.

    Returns a payload with `balanced` (bool), `element_imbalance`
    (products minus reactants, per element symbol, non-zero entries only) and
    `charge_imbalance` (net formal charge of products minus reactants).
    `balanced` is None when a SMILES cannot be parsed.
    """
    from collections import Counter
    from rdkit import Chem

    def _tally(smiles_list: Sequence[str]):
        atoms: "Counter[str]" = Counter()
        charge = 0
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None, None
            charge += Chem.GetFormalCharge(mol)
            mol = Chem.AddHs(mol)
            for atom in mol.GetAtoms():
                atoms[atom.GetSymbol()] += 1
        return atoms, charge

    r_atoms, r_charge = _tally(reactants)
    p_atoms, p_charge = _tally(products)
    if r_atoms is None or p_atoms is None:
        return {"balanced": None, "element_imbalance": {}, "charge_imbalance": None}

    imbalance = {
        symbol: p_atoms.get(symbol, 0) - r_atoms.get(symbol, 0)
        for symbol in set(r_atoms) | set(p_atoms)
        if p_atoms.get(symbol, 0) != r_atoms.get(symbol, 0)
    }
    charge_imbalance = p_charge - r_charge
    return {
        "balanced": not imbalance and charge_imbalance == 0,
        "element_imbalance": imbalance,
        "charge_imbalance": charge_imbalance,
    }


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
                          reaction_family: Optional[str] = None,
                          reactant_labels: Optional[Sequence[str]] = None) -> float:
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
        if reaction_family:
            multiplier *= get_donor_reactivity_multiplier(reaction_family, reactant_labels=reactant_labels)
            
        # 2. Kirkwood-Onsager Solvent Scaling (Phase 12.1 + O.1)
        # Simplified: Polar transitions (Maillard) are accelerated in polar solvents.
        # Barriers are adjusted by f(epsilon) = (eps-1)/(2eps+1).
        # We assume the input barrier is for water (eps=78.4, f=0.49).
        # Adjusted Barrier = Base_Barrier - Shift * (f(eps) - f(water))
        if conditions:
            eps = conditions.dielectric_constant
            f_eps = (eps - 1) / (2 * eps + 1)
            f_water = (78.4 - 1) / (2 * 78.4 + 1)
            
            # Phase O.1: Gate sensitivity by reaction family
            sensitivity = 0.0
            if reaction_family:
                rf = reaction_family.lower()
                if "condensation" in rf or "addition" in rf or "schiff" in rf:
                    sensitivity = 5.0  # Highly dependent on solvent stabilization
                elif "elimination" in rf or "dehydration" in rf or "thermolysis" in rf:
                    sensitivity = 3.0  # Polar transition states, but mostly unimolecular
                elif "cleavage" in rf or "strecker" in rf:
                    sensitivity = 1.0  # Mostly entropic, less charge separation
                else:
                    sensitivity = 2.0  # Generic default
            else:
                sensitivity = 5.0 # Fallback for legacy
                
            barrier_shift = sensitivity * (f_eps - f_water)
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

    def get_reaction_thermo(self, reactants: List[str], products: List[str], 
                            temperature_k: Optional[float] = None) -> Dict[str, float]:
        """
        Phase 1: Thermodynamic Governance.
        Calculates Delta G, Delta H, and Delta S for a reaction using Joback increments.

        Basis (audit 2026-08-26, item 3.3 — see the module header for the full
        derivation): Joback supplies formation properties at 298.15 K, so both
        the enthalpy and the entropy carried per species are element-referenced
        at 298.15 K. Enthalpy and entropy are therefore on the *same* basis, and
        the element reference cancels exactly when the reaction is atom and
        charge balanced. The reaction quantities returned here are then the true
            dH_rxn(T) = dH_rxn(298) + int_298^T dCp_rxn dT
            dS_rxn(T) = dS_rxn(298) + int_298^T dCp_rxn/T dT
            dG_rxn(T) = dH_rxn(T) - T * dS_rxn(T)
        For an atom-unbalanced step the element reference does not cancel; the
        payload flags this via `atom_balanced` / `thermo_basis` and
        NEUTRALIZE_UNBALANCED_THERMO_GATE decides whether the contaminated
        number is still handed to the caller's gate.
        """
        T = temperature_k or self.T

        balance = reaction_element_balance(reactants, products)
        atom_balanced = balance.get("balanced")

        try:
            r_res = [JobackEstimator.estimate(s) for s in reactants]
            p_res = [JobackEstimator.estimate(s) for s in products]

            def _sensible_h(cp_coeffs, temp: float) -> float:
                """int_298^T Cp dT for Cp = a + bT + cT^2 + dT^3 (J/mol)."""
                def int_cp(t):
                    return (cp_coeffs[0]*t + 0.5*cp_coeffs[1]*t**2
                            + (1/3)*cp_coeffs[2]*t**3 + 0.25*cp_coeffs[3]*t**4)
                return int_cp(temp) - int_cp(298.15)

            def _sensible_s(cp_coeffs, temp: float) -> float:
                """int_298^T Cp/T dT for the same polynomial (J/(mol K))."""
                def int_cp_over_t(t):
                    return (cp_coeffs[0]*np.log(t) + cp_coeffs[1]*t
                            + 0.5*cp_coeffs[2]*t**2 + (1/3)*cp_coeffs[3]*t**3)
                return int_cp_over_t(temp) - int_cp_over_t(298.15)

            def side_totals(res_list, temp: float) -> Tuple[float, float]:
                """Element-referenced H(T) and S(T) summed over one side."""
                total_h = 0.0
                total_s = 0.0
                for r in res_list:
                    cp_coeffs = r["cp_coeffs"]
                    # H on the formation basis: dfH(298) plus sensible heat.
                    total_h += r["H298"] + _sensible_h(cp_coeffs, temp)
                    # S on the SAME formation basis: dfS(298) = (dfH - dfG)/298.15
                    # plus the sensible entropy of the same species.
                    s298 = (r["H298"] - r["G298"]) / 298.15
                    total_s += s298 + _sensible_s(cp_coeffs, temp)
                return total_h, total_s

            rh, rs = side_totals(r_res, T)
            ph, ps = side_totals(p_res, T)

            delta_h = ph - rh
            delta_s = ps - rs
            delta_g = delta_h - T * delta_s

            payload = {
                "delta_g_j_mol": delta_g,
                "delta_g_kcal_mol": delta_g / 4184.0,
                "delta_h_j_mol": delta_h,
                "delta_h_kcal_mol": delta_h / 4184.0,
                "delta_s_j_mol_k": delta_s,
                "is_spontaneous": delta_g < 0,
                "atom_balanced": atom_balanced,
                "element_imbalance": balance.get("element_imbalance", {}),
                "charge_imbalance": balance.get("charge_imbalance"),
                "thermo_basis": (
                    "formation_298K_element_referenced"
                    if atom_balanced
                    else "element_reference_uncancelled"
                ),
            }

            if atom_balanced is not True:
                # The element reference survives: dG here also contains the
                # formation free energy of the atoms that appear or vanish.
                logger.debug(
                    "Unbalanced reaction thermo %s -> %s: imbalance=%s charge=%s dG=%.2f kcal/mol",
                    reactants, products, balance.get("element_imbalance"),
                    balance.get("charge_imbalance"), payload["delta_g_kcal_mol"],
                )
                payload["delta_g_kcal_mol_raw"] = payload["delta_g_kcal_mol"]
                if NEUTRALIZE_UNBALANCED_THERMO_GATE:
                    payload["delta_g_j_mol"] = 0.0
                    payload["delta_g_kcal_mol"] = 0.0
                    payload["is_spontaneous"] = True
                    payload["thermo_basis"] = "unavailable_unbalanced"

            return payload
        except Exception as exc:
            logger.debug("Joback thermo failed for %s -> %s: %s", reactants, products, exc)
            # Gate-neutral fallback: an unknown dG must not silently kill a step.
            return {
                "delta_g_j_mol": 0.0,
                "delta_g_kcal_mol": 0.0,
                "delta_h_j_mol": 0.0,
                "delta_h_kcal_mol": 0.0,
                "delta_s_j_mol_k": 0.0,
                "is_spontaneous": True,
                "atom_balanced": atom_balanced,
                "element_imbalance": balance.get("element_imbalance", {}),
                "charge_imbalance": balance.get("charge_imbalance"),
                "thermo_basis": "unavailable_estimator_error",
            }

    def is_reaction_feasible(self, reactants: List[str], products: List[str], 
                             threshold_kcal_mol: float = 30.0,
                             temperature_k: Optional[float] = None) -> Tuple[bool, float]:
        """
        Phase 12.3: Dynamic Thermo-Gating.
        If Delta G > threshold, the reaction is considered unphysical (non-spontaneous).
        """
        thermo = self.get_reaction_thermo(reactants, products, temperature_k)
        dg = thermo.get("delta_g_kcal_mol", 0.0)
        return dg <= threshold_kcal_mol, dg

    def simulate_network_cantera(self, mechanism_yaml: str, initial_concentrations: Dict[str, float], 
                                 time_span: Tuple[float, float], temperature_k: Optional[float] = None,
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
            # Shielded top-level-ish import for linter visibility while preserving env-safety
            ct = None

        # 1. Load the mechanism
        phase = ct.Solution(mechanism_yaml)
        
        # Initial T
        if temperature_profile:
            t_pts, t_vals = zip(*temperature_profile)
            T_init = float(np.interp(time_span[0], t_pts, t_vals))
        else:
            T_init = temperature_k or self.T
            
        phase.TP = T_init, 101325  # Standard pressure
        
        # 2. Set initial state
        phase.X = initial_concentrations
        
        # 3. Create a batch reactor
        # R.1: Use ConstPressureMoleReactor — compatible with incompressible
        # condensed-phase models (ideal-condensed). Unlike MoleReactor (constant
        # volume), this reactor varies volume at constant pressure, avoiding
        # the "setDensity not available" error on incompressible phases.
        r = ct.ConstPressureMoleReactor(phase, energy='off', clone=False)
        sim = ct.ReactorNet([r])

        
        # 4. Integrate
        n_points = 100
        time_points = np.linspace(time_span[0], time_span[1], n_points)
        results = {name: [] for name in phase.species_names}
        # Also store mole fractions for T/P invariant conversion math
        for name in phase.species_names:
            results[f"{name}_X"] = []
            
        results["time"] = []
        results["temperature"] = []
        
        for t in time_points:
            # Update temperature if profile provided
            if temperature_profile:
                T_curr = float(np.interp(t, t_pts, t_vals))
                phase.TP = T_curr, 101325
                r.syncState() 
            
            sim.advance(t)
            r.syncState()
            
            results["time"].append(t)
            results["temperature"].append(phase.T)
            for i, name in enumerate(phase.species_names):
                results[name].append(phase.concentrations[i])
                results[f"{name}_X"].append(phase.X[i])

                
        # Convert to numpy arrays
        return {k: np.array(v) for k, v in results.items()}

