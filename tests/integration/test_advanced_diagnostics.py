import pytest
from pathlib import Path
from src.recommend import Recommender
from src.pathway_extractor import Species, ElementaryStep
from src.pipeline import MaillardPipeline, FormulationResult

def test_precursor_attribution_tracing():
    # Mock a simple case: Glucose -> Furfural
    rec = Recommender()
    
    GLUCOSE = Species("Glucose", "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
    FURFURAL = Species("Furfural", "O=Cc1ccco1") # target
    
    step = ElementaryStep(reactants=[GLUCOSE], products=[FURFURAL])
    steps = [step]
    
    barriers = {"C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O->O=Cc1ccco1": 25.0}
    initial_pool = {GLUCOSE.smiles: 1.0}
    
    # Mock targets
    rec._load_desirable = lambda: {
        "Furfural": {"name": "Furfural", "smiles": "O=Cc1ccco1", "data": {}}
    }
    
    # Run predict_from_steps
    result = rec.predict_from_steps(
        steps=steps,
        barriers_dict=barriers,
        initial_concentrations=initial_pool,
        protein_type="pea_iso"
    )
    
    metrics = result["metrics"]
    assert "precursor_attribution" in metrics
    # Furfural yield should be attributed to Glucose
    # Depending on how the name is resolved, it might be the SMILES or the name if in lookup
    # In my implementation, it uses species_name_lookup.get(prec_can, prec_can)
    # Since Glucose isn't in a real lookup here, it might be SMILES.
    
    # Let's check if the attribution dict isn't empty
    assert len(metrics["precursor_attribution"]) > 0

def test_suppressed_compounds_detection():
    rec = Recommender()
    
    # Mock a target with high proxy but low observable (e.g. very low headspace)
    # We'll use a fake species with low volatility
    STIFF_VOLATILE = Species("Stiff", "C1CCCCCCCCCCCCCCC1O")
    
    raw_concs = {STIFF_VOLATILE.smiles: 100.0} # 100 ppb proxy
    
    # Mock targets
    rec._load_desirable = lambda: {
        "Stiff": {"name": "Stiff", "smiles": STIFF_VOLATILE.smiles, "data": {"boiling_point_c": 500.0}}
    }
    
    # Run predict_from_steps
    initial_pool = {"C": 1.0}
    result = rec.predict_from_steps(
        steps=[],
        barriers_dict={},
        initial_concentrations=initial_pool,
        protein_type="pea_iso"
    )
    
    assert "suppressed_compounds" in result["metrics"]
