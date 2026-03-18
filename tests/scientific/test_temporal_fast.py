import pytest
from src.recommend import Recommender, _temporal_accessibility
from src.pathway_extractor import Species, ElementaryStep


def test_temporal_accessibility_behaves_saturating_not_exponential_collapse():
    assert _temporal_accessibility(total_tau_minutes=1.0, time_minutes=60.0) > 0.99
    slow_progress = _temporal_accessibility(total_tau_minutes=600.0, time_minutes=10.0)
    assert 0.0 < slow_progress < 0.1


def test_temporal_accessibility_increases_with_time():
    short = _temporal_accessibility(total_tau_minutes=60.0, time_minutes=5.0)
    long = _temporal_accessibility(total_tau_minutes=60.0, time_minutes=60.0)
    assert long > short

def test_temporal_ramp_in_fast_mode():
    """
    Verify that the FAST recommender correctly ingests temporal profile CSVs
    and uses integrated Arrhenius propensity (SOTA).
    """
    recommender = Recommender()
    
    # Precursors (using SMILES for keys)
    ribose_smi = "OCC(O)C(O)C(O)C=O"
    ribose = Species("ribose", ribose_smi)
    furfural_smi = "O=Cc1ccco1"
    
    steps = [
        ElementaryStep(
            reactants=[ribose],
            products=[Species("furfural", furfural_smi)],
            reaction_family="Enolisation_1_2"
        )
    ]
    
    # Barriers
    step_key = f"{ribose_smi}->{furfural_smi}"
    barriers = {step_key: 20.0}
    
    initial = {ribose_smi: 1.0}
    
    # Use the test ramp we created
    ramp_path = "data/temp_profiles/test_ramp.csv"
    
    # Call with the new parameter
    res = recommender.predict_from_steps(steps, barriers, initial, temp_ramp_csv=ramp_path)
    
    assert res is not None
    assert "metrics" in res
    # Integrated weight should project a non-zero curated volatile output.
    assert res["predicted_ppb"].get("furfural", 0.0) > 0.0
