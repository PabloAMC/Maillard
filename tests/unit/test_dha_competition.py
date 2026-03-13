import pytest
from src.recommend import Recommender
from src.pathway_extractor import Species, ElementaryStep
from pathlib import Path

def test_lysine_budget_scaling():
    """
    Verify that lysine_budget_dha increases when DHA pathway is active.
    This tests the competition model in recommend.py.
    """
    recommender = Recommender()
    
    # Precursors
    ribose = Species("ribose", "OCC(O)C(O)C(O)C=O")
    lysine = Species("lysine", "NCCCCC(N)C(=O)O")
    cysteine = Species("cysteine", "NC(CS)C(=O)O")
    
    # Steps
    # 1. Maillard (Schiff Base)
    step_maillard = ElementaryStep(
        reactants=[ribose, lysine],
        products=[Species("Schiff", "OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O")],
        reaction_family="Schiff_Base_Formation"
    )
    
    # 2. DHA formation from Cysteine
    dha = Species("dha", "C=C(N)C(=O)O")
    step_dha_form = ElementaryStep(
        reactants=[cysteine],
        products=[dha, Species("H2S", "S")],
        reaction_family="Beta_Elimination"
    )
    
    # 3. DHA consuming Lysine
    step_dha_cons = ElementaryStep(
        reactants=[dha, lysine],
        products=[Species("LAL", "NC(CCCCNCC(N)C(=O)O)C(=O)O")],
        reaction_family="DHA_Crosslinking"
    )
    
    steps = [step_maillard, step_dha_form, step_dha_cons]
    
    # Barriers (make DHA slightly slower so it competes)
    # KEYS MUST BE SORTED BY SMILES TO MATCH recommend.py logic
    barriers = {
        "NCCCCC(N)C(=O)O+OCC(O)C(O)C(O)C=O->OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O": 20.0,
        "NC(CS)C(=O)O->C=C(N)C(=O)O+S": 22.0,
        "C=C(N)C(=O)O+NCCCCC(N)C(=O)O->NC(CCCCNCC(N)C(=O)O)C(=O)O": 21.0
    }
    
    initial = {
        "OCC(O)C(O)C(O)C=O": 1.0, 
        "NCCCCC(N)C(=O)O": 1.0, 
        "NC(CS)C(=O)O": 1.0
    }
    
    res = recommender.predict_from_steps(steps, barriers, initial)
    
    budget = res["metrics"]["lysine_budget_dha"]
    # We expect some budget to be consumed by DHA
    assert budget > 0.0
    assert budget < 100.0
    print(f"Lysine Budget (DHA %): {budget:.2f}%")
