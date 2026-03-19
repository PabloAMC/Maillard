"""
tests/unit/test_headspace.py

Verifies headspace partitioning and matrix effect logic.
"""

import pytest
from src.headspace import HeadspaceModel  # noqa: E402
from src.matrix_calibration_registry import get_matrix_runtime_composition_policy

def test_headspace_temperature_scaling():
    """Verify that volatility increases with temperature."""
    model = HeadspaceModel()
    # Hexanal at 25C vs 100C
    c_25 = model.predict_headspace({"Hexanal": 1.0}, 25.0)["Hexanal"]
    c_100 = model.predict_headspace({"Hexanal": 1.0}, 100.0)["Hexanal"]
    
    assert c_100 > c_25 * 5.0  # Volatility should increase significantly

def test_lipid_suppression():
    """Verify that hydrophobic compounds are suppressed by fat."""
    model = HeadspaceModel()
    matrix = {"Hexanal": 1.0, "Furfural": 1.0}
    
    # 0% fat
    air_0 = model.predict_headspace(matrix, 25.0, fat_fraction=0.0)
    # 10% fat
    air_10 = model.predict_headspace(matrix, 25.0, fat_fraction=0.1)
    
    # Hexanal is hydrophobic (Kfat ~ 450)
    # 1 / (1 + 450*0.1) = 1/46th
    assert air_10["Hexanal"] < air_0["Hexanal"] / 40.0
    
    # Furfural is polar (Kfat ~ 5)
    # 1 / (1 + 5*0.1) = 1/1.5
    assert air_10["Furfural"] > air_0["Furfural"] / 2.0
    
    # Hexanal should be MUCH more suppressed than Furfural
    ratio_hex = air_0["Hexanal"] / air_10["Hexanal"]
    ratio_fur = air_0["Furfural"] / air_10["Furfural"]
    assert ratio_hex > ratio_fur * 10.0

def test_protein_sequestration():
    """Verify that protein fraction reduces headspace concentrations."""
    model = HeadspaceModel()
    matrix = {"Methional": 1.0}
    
    # 0% protein
    air_0 = model.predict_headspace(matrix, 25.0, protein_fraction=0.0)
    # 20% protein (typical for meat analog)
    air_20 = model.predict_headspace(matrix, 25.0, protein_fraction=0.2)
    
    # Since Methional Kprot=2.0 -> 1 / (1 + 2.0*0.2) = 1/1.4
    assert air_20["Methional"] == pytest.approx(air_0["Methional"] / 1.4)


def test_matrix_retention_fallback_uses_pea_and_soy_profiles_when_fractions_are_unspecified():
    model = HeadspaceModel()
    matrix = {"Furfural": 1.0}

    air_free = model.predict_headspace(matrix, 25.0, protein_type="free")
    air_pea = model.predict_headspace(matrix, 25.0, protein_type="pea_iso")
    air_soy = model.predict_headspace(matrix, 25.0, protein_type="soy_iso")

    # Furfural is a furan. For pea_iso at denaturation 0.5: 
    # base=0.50, class_factor=0.945 -> effective=0.4725
    assert air_pea["Furfural"] == pytest.approx(air_free["Furfural"] * 0.4725)
    # For soy_iso at denaturation 0.5:
    # base=0.55, class_factor=0.975 -> effective=0.53625
    assert air_soy["Furfural"] == pytest.approx(air_free["Furfural"] * 0.53625)
    assert air_soy["Furfural"] > air_pea["Furfural"]


def test_denaturation_state_relaxes_headspace_fallback_when_fractions_are_unspecified():
    model = HeadspaceModel()
    matrix = {"Furfural": 1.0}

    air_native = model.predict_headspace(matrix, 25.0, protein_type="pea_iso", denaturation_state=0.0)
    air_mid = model.predict_headspace(matrix, 25.0, protein_type="pea_iso", denaturation_state=0.5)
    air_denatured = model.predict_headspace(matrix, 25.0, protein_type="pea_iso", denaturation_state=1.0)

    assert air_native["Furfural"] < air_mid["Furfural"] < air_denatured["Furfural"]


def test_acidic_ph_increases_plant_matrix_release_for_acid_sensitive_off_flavour_markers():
    model = HeadspaceModel()
    matrix = {
        "Hexanal": 1.0,
        "Nonanal": 1.0,
        "2-Pentylfuran": 1.0,
        "2,5-Dimethylpyrazine": 1.0,
    }

    air_acid = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=4.5)
    air_neutral = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=6.5)

    assert air_acid["Hexanal"] / air_neutral["Hexanal"] == pytest.approx(1.6, rel=0.08)
    assert air_acid["Nonanal"] / air_neutral["Nonanal"] == pytest.approx(1.6, rel=0.08)
    assert air_acid["2-Pentylfuran"] / air_neutral["2-Pentylfuran"] == pytest.approx(1.6, rel=0.08)
    assert air_acid["2,5-Dimethylpyrazine"] == pytest.approx(air_neutral["2,5-Dimethylpyrazine"])


def test_pratap_singh_headspace_calibration_carries_soy_release_gap_in_headspace_layer():
    model = HeadspaceModel()

    assert model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="pea_iso",
        pH=6.0,
    ) == pytest.approx(1.0)
    assert model.get_matrix_benchmark_headspace_factor(
        "2-Pentylfuran",
        protein_type="pea_iso",
        pH=6.0,
    ) == pytest.approx(1.0)
    assert model.get_matrix_benchmark_headspace_factor(
        "1-Hexanol",
        protein_type="pea_iso",
        pH=6.0,
    ) == pytest.approx(1.0)

    assert model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="soy_iso",
        pH=6.0,
    ) == pytest.approx(0.453 / 0.205)
    assert model.get_matrix_benchmark_headspace_factor(
        "2-Pentylfuran",
        protein_type="soy_iso",
        pH=6.0,
    ) == pytest.approx(2.972 / 0.502)
    assert model.get_matrix_benchmark_headspace_factor(
        "1-Hexanol",
        protein_type="soy_iso",
        pH=6.0,
    ) == pytest.approx(0.143 / 0.063)


def test_shu_heated_soy_calibration_suppresses_2_pentylfuran_to_detection_floor_surrogate():
    model = HeadspaceModel()

    ambient = model.get_matrix_benchmark_headspace_factor(
        "2-Pentylfuran",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=40.0,
        time_minutes=10.0,
    )
    heated = model.get_matrix_benchmark_headspace_factor(
        "2-Pentylfuran",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=120.0,
        time_minutes=20.0,
    )

    assert ambient == pytest.approx(2.972 / 0.502)
    assert heated == pytest.approx((2.972 / 0.502) * 0.03)
    assert heated < ambient * 0.05


def test_dynamic_soy_hexanal_release_increases_with_temperature_before_heavy_time_attenuation():
    model = HeadspaceModel()

    cool = model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=25.0,
        time_minutes=2.0,
    )
    hot = model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=95.0,
        time_minutes=2.0,
    )

    assert hot > cool


def test_dynamic_soy_hexanal_release_at_95c_is_temporally_attenuated():
    model = HeadspaceModel()

    short = model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=95.0,
        time_minutes=2.0,
    )
    long = model.get_matrix_benchmark_headspace_factor(
        "Hexanal",
        protein_type="soy_iso",
        pH=6.0,
        temperature_celsius=95.0,
        time_minutes=30.0,
    )

    assert long < short


def test_runtime_composition_policy_preserves_ambient_hexanal_baseline_but_activates_thermal_states():
    ambient = get_matrix_runtime_composition_policy(
        "Hexanal",
        protein_type="soy_iso",
        process_state="ambient_slurry",
    )
    intermediate = get_matrix_runtime_composition_policy(
        "Hexanal",
        protein_type="soy_iso",
        process_state="intermediate_matrix",
    )
    heated = get_matrix_runtime_composition_policy(
        "Hexanal",
        protein_type="soy_iso",
        process_state="heated_matrix",
    )

    assert ambient["mode"] == "static_observable_calibration"
    assert intermediate["mode"] == "compose_dynamic_retention"
    assert heated["mode"] == "compose_dynamic_retention"


def test_explicit_matrix_fractions_override_retention_fallback():
    model = HeadspaceModel()
    matrix = {"Methional": 1.0}

    air_with_fraction_only = model.predict_headspace(
        matrix,
        25.0,
        protein_fraction=0.2,
    )
    air_with_fraction_and_type = model.predict_headspace(
        matrix,
        25.0,
        protein_fraction=0.2,
        protein_type="pea_iso",
    )

    assert air_with_fraction_and_type["Methional"] == pytest.approx(
        air_with_fraction_only["Methional"]
    )

if __name__ == "__main__":
    pytest.main([__file__])
