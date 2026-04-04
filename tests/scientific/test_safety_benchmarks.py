import pytest
from src.safety import (
    evaluate_formulation_safety,
    get_safety_reference_entries,
    predict_acrylamide,
    predict_cel,
    predict_cml,
    predict_furosine,
)

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


def test_acrylamide_moisture_dependency_flips_between_lme_and_hme():
    lme_low = predict_acrylamide(10.0, 10.0, 150, 20, 6.5, water_activity=0.25, moisture_regime="lme")
    lme_high = predict_acrylamide(10.0, 10.0, 150, 20, 6.5, water_activity=0.38, moisture_regime="lme")
    hme_low = predict_acrylamide(10.0, 10.0, 150, 20, 6.5, water_activity=0.55, moisture_regime="hme")
    hme_high = predict_acrylamide(10.0, 10.0, 150, 20, 6.5, water_activity=0.90, moisture_regime="hme")

    assert lme_high.acrylamide_ppb > lme_low.acrylamide_ppb
    assert hme_high.acrylamide_ppb < hme_low.acrylamide_ppb


def test_extrusion_damage_load_is_folded_into_safety_result():
    precursors = {"Asparagine": 1.0, "Glucose": 5.0}
    risk, flagged = evaluate_formulation_safety(
        precursors,
        150,
        15,
        6.5,
        modifiers={
            "__extrusion_process__": {
                "moisture_regime": "hme",
                "water_activity": 0.55,
                "effective_temperature_celsius": 165.0,
                "total_damage_load": {
                    "furosine_mg_per_kg": 20.0,
                    "lal_mg_per_kg": 50.0,
                },
            }
        },
    )

    assert risk > 0
    assert "Furosine" in flagged
    assert "LAL" in flagged


def test_cml_and_cel_follow_report12_directionality():
    cml_dry = predict_cml(45.0, 45.0, 145.0, 20.0, water_activity=0.35)
    cml_wet = predict_cml(45.0, 45.0, 145.0, 20.0, water_activity=0.75)
    cel_low_temp = predict_cel(45.0, 45.0, 135.0, 20.0, water_activity=0.45)
    cel_high_temp = predict_cel(45.0, 45.0, 165.0, 20.0, water_activity=0.45)
    cel_dry = predict_cel(45.0, 45.0, 150.0, 20.0, water_activity=0.30)
    cel_wet = predict_cel(45.0, 45.0, 150.0, 20.0, water_activity=0.85)

    assert cml_wet > cml_dry
    assert cel_high_temp > cel_low_temp
    assert cel_dry > cel_wet


def test_furosine_exhibits_crossover_behavior():
    peak = predict_furosine(
        140.0,
        20.0,
        lysine_mM=35.0,
        reducing_sugar_mM=35.0,
        protein_type="pea_iso",
        water_activity=0.55,
    )
    severe = predict_furosine(
        165.0,
        20.0,
        lysine_mM=35.0,
        reducing_sugar_mM=35.0,
        protein_type="pea_iso",
        water_activity=0.55,
    )

    assert 60.0 < peak < 110.0
    assert peak > severe


def test_composite_age_reference_entries_are_runtime_visible():
    cml_entries = get_safety_reference_entries(analyte="cml", visibility="all")
    cel_entries = get_safety_reference_entries(analyte="cel", visibility="all")
    furosine_entries = get_safety_reference_entries(analyte="furosine", visibility="all")

    assert any(entry["id"] == "foods_2023_pba_cml_cel_benchmark" for entry in cml_entries)
    assert any(entry["id"] == "pmc_2024_pba_cml_cel_ranges" for entry in cel_entries)
    assert any(entry["id"] == "ramirez_jimenez_2000_furosine_crossover" for entry in furosine_entries)
