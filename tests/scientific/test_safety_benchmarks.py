import pytest
from src.safety import predict_acrylamide, evaluate_formulation_safety

def test_acrylamide_non_monotonic_behavior():
    """
    Verify that at extremely high temperatures or long times,
    acrylamide levels off due to elimination kinetics.
    """
    asparagine = 10.0
    glucose = 10.0
    pH = 7.0
    
    # 180C for 20 mins
    res_20min = predict_acrylamide(asparagine, glucose, 180, 20, pH)
    # 180C for 200 mins
    res_200min = predict_acrylamide(asparagine, glucose, 180, 200, pH)
    
    # Without elimination, 200min would be exactly 10x 20min.
    # With elimination, it should be significantly less than 10x.
    ratio = res_200min.acrylamide_ppb / res_20min.acrylamide_ppb if res_20min.acrylamide_ppb > 0 else 0
    assert ratio < 10.0
    print(f"Non-monotonic Ratio (200min/20min): {ratio:.2f}")

def test_acrylamide_ph_dependency():
    """Verify that acrylamide formation increases with pH."""
    asparagine = 10.0
    glucose = 10.0
    temp = 150
    time = 20
    
    res_ph5 = predict_acrylamide(asparagine, glucose, temp, time, 5.0)
    res_ph8 = predict_acrylamide(asparagine, glucose, temp, time, 8.0)
    
    assert res_ph8.acrylamide_ppb > res_ph5.acrylamide_ppb

def test_evaluate_formulation_safety_integration():
    """Verify integration of acrylamide in the formulation safety scoring."""
    precursors = {"Asparagine": 1.0, "Glucose": 5.0}
    risk, flagged = evaluate_formulation_safety(precursors, 180, 15, 6.5)
    
    assert risk > 0
    assert "Acrylamide" in flagged
