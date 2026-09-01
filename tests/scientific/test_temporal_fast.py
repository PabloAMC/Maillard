import pytest
from src.recommend import Recommender, _temporal_accessibility
from src.chem_utils import Species, ElementaryStep


def test_temporal_accessibility_behaves_saturating_not_exponential_collapse():
    assert _temporal_accessibility(total_tau_minutes=1.0, time_minutes=60.0) > 0.99
    slow_progress = _temporal_accessibility(total_tau_minutes=600.0, time_minutes=10.0)
    assert 0.0 < slow_progress < 0.1


def test_temporal_accessibility_increases_with_time():
    short = _temporal_accessibility(total_tau_minutes=60.0, time_minutes=5.0)
    long = _temporal_accessibility(total_tau_minutes=60.0, time_minutes=60.0)
    assert long > short

def test_temporal_ramp_in_fast_mode(tmp_path, caplog):
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
    
    # A real ramp file. Until 2026-09-01 this test pointed at a path that did not exist
    # (data/temp_profiles/test_ramp.csv); _load_ramp logged a warning, returned [], and the
    # test silently exercised the NON-ramp path while claiming to test ramps.
    ramp_path = tmp_path / "test_ramp.csv"
    ramp_path.write_text("time,temp\n0,25\n10,100\n30,140\n45,140\n", encoding="utf-8")

    with caplog.at_level("WARNING"):
        res = recommender.predict_from_steps(steps, barriers, initial, temp_ramp_csv=str(ramp_path))
    assert "Failed to load ramp" not in caplog.text
    
    assert res is not None
    assert "metrics" in res
    # Integrated weight should project a non-zero curated volatile output.
    assert res["predicted_ppb"].get("furfural", 0.0) > 0.0
