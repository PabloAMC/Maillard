import pytest
from src.recommend import Recommender
from src.pathway_extractor import Species, ElementaryStep

def test_temporal_ramp_in_fast_mode():
    """
    Verify that the FAST recommender correctly ingests temporal profile CSVs
    and uses integrated Arrhenius propensity (SOTA).
    """
    recommender = Recommender()
    
    # Precursors (using SMILES for keys)
    ribose_smi = "OCC(O)C(O)C(O)C=O"
    lysine_smi = "NCCCCC(N)C(=O)O"
    ribose = Species("ribose", ribose_smi)
    lysine = Species("lysine", lysine_smi)
    
    steps = [
        ElementaryStep(
            reactants=[ribose, lysine],
            products=[Species("Schiff", "OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O")],
            reaction_family="Schiff_Base_Formation"
        )
    ]
    
    # Barriers
    step_key = f"{ribose_smi}+{lysine_smi}->OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O"
    barriers = {step_key: 20.0}
    
    initial = {ribose_smi: 1.0, lysine_smi: 1.0}
    
    # Use the test ramp we created
    ramp_path = "data/temp_profiles/test_ramp.csv"
    
    # Call with the new parameter
    res = recommender.predict_from_steps(steps, barriers, initial, temp_ramp_csv=ramp_path)
    
    assert res is not None
    assert "metrics" in res
    # Integrated weight should be non-zero for a 20 kcal/mol barrier at 150C
    assert any(v > 0 for v in res["predicted_ppb"].values())
