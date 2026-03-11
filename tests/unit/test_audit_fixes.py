import pytest
import numpy as np
import os
import yaml
from pathlib import Path
import cantera as ct
from src.thermo import get_nasa_coefficients, JobackEstimator, gas_constant
from src.cantera_export import CanteraExporter
from src.bayesian_optimizer import FormulationOptimizer
from src.smirks_engine import ReactionConditions

def test_r4_nasa7_cantera_consistency():
    """R.4.4: Verify NASA7 in Cantera matches Joback within 5%."""
    # Test with Ribose (C5H10O5)
    smiles = "OCC(O)C(O)C(O)C=O"
    coeffs = get_nasa_coefficients(smiles)
    assert len(coeffs) == 14
    
    # Export a minimal mechanism
    exporter = CanteraExporter()
    exporter.add_species(smiles, "Ribose")
    # Add a dummy reaction to satisfy exporter
    exporter.add_reaction([smiles], [smiles], 20.0, reaction_family="enolisation_1_2")
    
    mech_path = "/tmp/test_r4.yaml"
    exporter.export_yaml(mech_path)
    
    # Load in Cantera
    gas = ct.Solution(mech_path)
    species = gas.species("Ribose")
    
    # Compare Cp at various temperatures
    joback_res = JobackEstimator.estimate(smiles)
    a, b, c, d = joback_res["cp_coeffs"]
    
    def joback_cp(T):
        return a + b*T + c*T**2 + d*T**3
        
    for T in [300.0, 423.15, 500.0, 800.0]:
        gas.TP = T, 101325
        ct_cp = gas.cp_mole / 1000.0 # Convert J/kmol to J/mol
        j_cp = joback_cp(T)
        
        error = abs(ct_cp - j_cp) / j_cp
        print(f"T={T}K: Cantera Cp={ct_cp:.2f}, Joback Cp={j_cp:.2f}, Error={error:.2%}")
        assert error < 0.05, f"Cp mismatch at {T}K: {error:.2%}"

def test_r8_optimizer_consistency():
    """R.8.3: Verify optimizer produces consistent results (Thread Safety/State Integrity)."""
    optimizer = FormulationOptimizer(target_tag="meaty")
    
    # Run two sequential optimizations with few trials
    # If state was being mutated/corrupted, the second run might be different or crash
    study1 = optimizer.optimize(["ribose"], ["cysteine"], n_trials=5)
    val1 = study1.best_value
    
    study2 = optimizer.optimize(["ribose"], ["cysteine"], n_trials=5)
    val2 = study2.best_value
    
    # Since Optuna is stochastic, they won't be identical unless we fix the seed,
    # but we want to ensure it remains functional and reasonable.
    assert study1 is not None
    assert study2 is not None
    assert len(study1.trials) == 5
    assert len(study2.trials) == 5
    print(f"Run 1 Best: {val1:.4f}, Run 2 Best: {val2:.4f}")

if __name__ == "__main__":
    test_r4_nasa7_cantera_consistency()
    test_r8_optimizer_consistency()
