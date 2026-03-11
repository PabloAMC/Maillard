import pytest
from src.barrier_constants import FAST_BARRIERS, get_barrier
from src.recommend import Recommender
from src.smirks_engine import SmirksEngine
from src.conditions import ReactionConditions
from src.precursor_resolver import resolve_many

def test_barrier_constants_coverage():
    """Ensure we have constants for the major families we generate."""
    expected_families = [
        "amadori", "schiff", "strecker", "cysteine", 
        "retro", "enolisation", "beta", "thiazole"
    ]
    for fam in expected_families:
        bar, _ = get_barrier(fam)
        assert bar < 40.0, f"Missing barrier calibration for {fam}"
        
def run_model_system(sugars, aminos, ph=6.0, temp=150.0):
    precursors = resolve_many(sugars + aminos)
    cond = ReactionConditions(pH=ph, temperature_celsius=temp)
    engine = SmirksEngine(cond)
    steps = engine.enumerate(precursors, max_generations=4)
    
    barriers_dict = {}
    for step in steps:
        bar, unc = get_barrier(step.reaction_family)
        ph_mult = cond.get_ph_multiplier(step.reaction_family or "")
        bar /= ph_mult
        
        rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        barriers_dict[rxn_key] = (max(0.0, bar), unc)
        
    recommender = Recommender()
    # Convert list of precursor SMILES to dict with default concentration 1.0
    initial_conc = {p.smiles: 1.0 for p in precursors}
    res = recommender.predict_from_steps(steps, barriers_dict, initial_conc)
    
    # Sort targets by barrier ascending (lowest barrier = most favoured output)
    targets = sorted(res["targets"], key=lambda t: t["span"])
    return [t["name"] for t in targets]

class TestLiteratureValidationGate:
    """
    Phase 8.C.5: The Literature Validation Gate.
    Verifies the pipeline correctly ranks volatiles for 3 known model systems
    using the calibrated FAST-mode heuristic barriers.
    """
    
    def test_ribose_cysteine_system(self):
        """
        Model System 1: Ribose + Cysteine at pH 5.0, 150C
        Literature expectation: FFT (2-furfurylthiol) is dominant.
        """
        top_volatiles = run_model_system(["ribose"], ["cysteine"], ph=5.0)
        assert len(top_volatiles) > 0, "No volatiles generated"
        
        # FFT should be in the top 5 (it should ideally be #1 or #2, but ties can push it down)
        top_5 = top_volatiles[:5]
        assert any("FFT" in t for t in top_5), f"FFT missing from top 5. Got: {top_5}"

    def test_glucose_glycine_system(self):
        """
        Model System 2: Glucose + Glycine at pH 7.0, 150C
        Literature expectation: Pyrazines/HMF/Pyruvaldehyde are favoured. FFT should NOT be present.
        """
        top_volatiles = run_model_system(["glucose"], ["glycine"], ph=7.0)
        assert len(top_volatiles) > 0, "No volatiles generated"
        
        # FFT impossible without sulfur
        assert not any("FFT" in t for t in top_volatiles)
        
        # HMF or pyruvaldehyde-derived products/pyrazines should be highly ranked
        top_5 = top_volatiles[:5]
        assert any(t in ["HMF", "Pyruvaldehyde", "2,5-Dimethylpyrazine"] for t in top_5), f"Missing typical hexose/neutral-AA outputs. Got: {top_5}"
        
    def test_ribose_cysteine_leucine_system(self):
        """
        Model System 3: Ribose + Cysteine + Leucine at pH 5.0, 150C
        Literature expectation: FFT and Strecker aldehyde (3-methylbutanal) co-dominant.
        """
        top_volatiles = run_model_system(["ribose"], ["cysteine", "leucine"], ph=5.0)
        
        top_5 = top_volatiles[:5]
        
        # Should see both the sulfur route and the Strecker route
        assert any("FFT" in t for t in top_5)
        assert any("3-Methylbutanal" in t for t in top_5)
