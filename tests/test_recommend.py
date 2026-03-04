import pytest
from pathlib import Path
from src.recommend import Recommender

def test_recommender_canonical_systems():
    """
    Test that the Recommender properly identifies active pathways 
    for the canonical model systems.
    """
    results_path = Path(__file__).resolve().parents[1] / "results" / "curated_screening_results.json"
    
    # Needs the json file to exist
    if not results_path.exists():
        pytest.skip(f"Screening results not found at {results_path}. Cannot test recommender.")
        
    recommender = Recommender(results_path)
    
    # 1. Savory/Meat System (Ribose + Cysteine)
    savory_active = recommender.predict(["D-ribose", "L-cysteine"])
    savory_names = [p["name"] for p in savory_active]
    assert "C_S_Maillard_FFT" in savory_names
    # Shouldn't trigger Glycine pathways
    assert "A_Core_Maillard_Ribose_Gly" not in savory_names
    
    # 2. Baked System (Glucose + Glycine)
    baked_active = recommender.predict(["D-glucose", "glycine"])
    baked_names = [p["name"] for p in baked_active]
    assert "A_Core_Maillard_Glucose_Gly" in baked_names
    
    # Toxicity check: Glucose+Glycine produces HMF
    hmf_pathway = next(p for p in baked_active if p["name"] == "A_Core_Maillard_Glucose_Gly")
    assert hmf_pathway["toxicity"] is not None
    assert hmf_pathway["toxicity"]["name"] == "5-Hydroxymethylfurfural (HMF)"
    
    # 3. Strecker System (Ribose + Cys + Leu)
    cplx_active = recommender.predict(["D-ribose", "L-cysteine", "L-leucine"])
    cplx_names = [p["name"] for p in cplx_active]
    assert "C_S_Maillard_FFT" in cplx_names
    assert "B_Strecker_Leu" in cplx_names
    
    # 4. Off-flavor Trapping (Hexanal + Lysine)
    trap_active = recommender.predict(["hexanal", "L-lysine"])
    trap_names = [p["name"] for p in trap_active]
    assert "D_Offflavour_Trapping_Lys" in trap_names

def test_recommender_penalties():
    """Test that competing pathways correctly apply penalties to desirable ones."""
    results_path = Path(__file__).resolve().parents[1] / "results" / "curated_screening_results.json"
    if not results_path.exists():
        pytest.skip(f"Screening results not found")
        
    recommender = Recommender(results_path)
    
    # System with Cysteine and Lysine -> Triggers E_DHA_Competition
    # AND Ribose -> Triggers C_S_Maillard_FFT
    # Since DHA competes for Cysteine, FFT pathway should get a penalty.
    active = recommender.predict(["D-ribose", "L-cysteine", "L-lysine"])
    
    active_names = [p["name"] for p in active]
    assert "E_DHA_Competition" in active_names
    assert "C_S_Maillard_FFT" in active_names
    
    fft_pathway = next(p for p in active if p["name"] == "C_S_Maillard_FFT")
    dha_pathway = next(p for p in active if p["name"] == "E_DHA_Competition")
    
    # Verify toxicity flag on DHA
    assert dha_pathway.get("toxicity") is not None
    assert "Lysinoalanine" in dha_pathway["toxicity"]["name"]
    
    # Verify penalty on FFT (since they both share L-cysteine)
    assert fft_pathway["penalty"] in ["MEDIUM", "HIGH"]
