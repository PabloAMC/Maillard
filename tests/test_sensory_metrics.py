import pytest
import sys
from pathlib import Path

# Add project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.recommend import Recommender
from src.smirks_engine import SmirksEngine, ReactionConditions
from src.precursor_resolver import resolve

def test_trapping_efficiency_calculation():
    """Verify that Lysine traps Hexanal more efficiently than Glycine."""
    recommender = Recommender()
    engine = SmirksEngine(ReactionConditions(pH=6.0, temperature_celsius=150.0))
    
    # System A: Hexanal + Glycine
    hexanal = resolve("hexanal")
    glycine = resolve("glycine")
    steps_a = engine.enumerate([hexanal, glycine])
    
    # Mock barriers (simple heuristic matching FAST mode)
    barriers_a = {
        f"{'+'.join(sorted(r.smiles for r in s.reactants))}->{'+'.join(sorted(p.smiles for p in s.products))}": 15.0 # Schiff base
        for s in steps_a if s.reaction_family == "Lipid_Schiff_Base"
    }
    
    results_a = recommender.predict_from_steps(steps_a, barriers_a, [hexanal.smiles, glycine.smiles])
    eff_a = results_a["metrics"]["trapping_efficiency"]["Hexanal"]
    
    # System B: Hexanal + Lysine
    lysine = resolve("lysine")
    steps_b = engine.enumerate([hexanal, lysine])
    barriers_b = {
        f"{'+'.join(sorted(r.smiles for r in s.reactants))}->{'+'.join(sorted(p.smiles for p in s.products))}": 15.0 # Schiff base
        for s in steps_b if s.reaction_family == "Lipid_Schiff_Base"
    }
    
    results_b = recommender.predict_from_steps(steps_b, barriers_b, [hexanal.smiles, lysine.smiles])
    eff_b = results_b["metrics"]["trapping_efficiency"]["Hexanal"]
    
    print(f"\nTrapping Efficiency - Glycine: {eff_a:.2f}%")
    print(f"Trapping Efficiency - Lysine: {eff_b:.2f}%")
    
    # Lysine has 2 amino groups and should generate more Schiff base variants/steps
    # resulting in a higher Boltzmann sum.
    assert eff_b > eff_a

def test_sensory_metadata_presence():
    """Verify that predictions include sensory descriptors and thresholds."""
    recommender = Recommender()
    engine = SmirksEngine(ReactionConditions(pH=5.0, temperature_celsius=150.0))
    
    ribose = resolve("ribose")
    cysteine = resolve("cysteine")
    steps = engine.enumerate([ribose, cysteine])
    
    # Just need 1 target to check
    barriers = {}
    for s in steps:
        key = f"{'+'.join(sorted(r.smiles for r in s.reactants))}->{'+'.join(sorted(p.smiles for p in s.products))}"
        barriers[key] = 40.0
        
    results = recommender.predict_from_steps(steps, barriers, [ribose.smiles, cysteine.smiles])
    targets = results["targets"]
    
    fft_found = False
    for t in targets:
        if "Furfurylthiol" in t["name"]:
            fft_found = True
            assert "sensory" in t
            assert "threshold" in t
            assert t["threshold"] is not None
            assert "Roasted" in t["sensory"]
            
    assert fft_found, "FFT not found in Ribose+Cysteine products"

if __name__ == "__main__":
    # For manual debugging
    test_trapping_efficiency_calculation()
    test_sensory_metadata_presence()
    print("All sensory metric tests passed!")
