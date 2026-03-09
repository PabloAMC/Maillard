import pytest
import sys
from pathlib import Path

# Add project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.recommend import Recommender
from src.smirks_engine import SmirksEngine, ReactionConditions
from src.precursor_resolver import resolve

def test_lysine_budget_competition():
    """Verify that Lysine Budget increases when Serine (DHA precursor) is added."""
    recommender = Recommender()
    engine = SmirksEngine(ReactionConditions(pH=6.0, temperature_celsius=150.0))
    
    # System A: Glucose + Lysine (Maillard only)
    glucose = resolve("glucose")
    lysine = resolve("lysine")
    steps_a = engine.enumerate([glucose, lysine])
    
    # Mock barriers (simple)
    barriers_a = {
        f"{'+'.join(sorted(r.smiles for r in s.reactants))}->{'+'.join(sorted(p.smiles for p in s.products))}": 20.0
        for s in steps_a
    }
    
    initial_conc_a = {glucose.smiles: 1.0, lysine.smiles: 1.0}
    results_a = recommender.predict_from_steps(steps_a, barriers_a, initial_conc_a)
    budget_a = results_a["metrics"]["lysine_budget_dha"]
    
    # System B: Glucose + Lysine + Serine (Maillard + DHA)
    serine = resolve("serine")
    steps_b = engine.enumerate([glucose, lysine, serine])
    barriers_b = {
        f"{'+'.join(sorted(r.smiles for r in s.reactants))}->{'+'.join(sorted(p.smiles for p in s.products))}": 20.0
        for s in steps_b
    }
    
    initial_conc_b = {glucose.smiles: 1.0, lysine.smiles: 1.0, serine.smiles: 1.0}
    results_b = recommender.predict_from_steps(steps_b, barriers_b, initial_conc_b)
    budget_b = results_b["metrics"]["lysine_budget_dha"]
    
    # With FAST mode heuristics, DHA (18 kcal) competes with Schiff Base (15 kcal)
    # No Serine -> No DHA steps -> Budget = 0
    # With Serine -> DHA steps exist -> Budget > 0
    assert budget_a == 0.0
    assert budget_b > 0.0

if __name__ == "__main__":
    test_lysine_budget_competition()
    print("Lysine budget test passed!")
