import pytest
from src.recommend import Recommender
from src.pathway_extractor import Species, ElementaryStep

def test_temporal_ramp_gap_in_fast_mode():
    """
    DOCUMENTED GAP: The FAST recommender does not yet ingest temporal profile CSVs. 
    It evaluates at a single temperature point.
    """
    recommender = Recommender()
    
    # Precursors
    ribose = Species("ribose", "OCC(O)C(O)C(O)C=O")
    lysine = Species("lysine", "NCCCCC(N)C(=O)O")
    
    steps = [
        ElementaryStep(
            reactants=[ribose, lysine],
            products=[Species("Schiff", "OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O")],
            reaction_family="Schiff_Base_Formation"
        )
    ]
    
    barriers = {
        "OCC(O)C(O)C(O)C=O+NCCCCC(N)C(=O)O->OCC(O)C(O)C(O)C=NCCCCC(N)C(=O)O": 20.0
    }
    
    initial = {"ribose": 1.0, "lysine": 1.0}
    
    # We want to pass a CSV ramp like the CLI supports
    ramp_path = "data/temp_profiles/isothermal_150.csv"
    
    # Current predict_from_steps only takes float temperature_kelvin
    # Passing a path should either work (if we implement it) or fail clearly.
    # For now, we xfail to mark this feature request.
    pytest.xfail("Recommender.predict_from_steps does not yet support temperature ramp CSVs.")
    
    # This will fail at runtime if we try to pass a string where it expects a float
    res = recommender.predict_from_steps(steps, barriers, initial, temperature_kelvin=ramp_path)
    assert res is not None
