import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pipeline import MaillardPipeline
from src.smirks_engine import ReactionConditions


def _matrix_accessibility_probe(
    protein_type: str,
    denaturation_state: float | None,
    *,
    temp_c: float = 140.0,
    time_minutes: float = 30.0,
    pH: float = 5.5,
):
    designer = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    formulation = {
        "name": f"matrix_probe_{protein_type}_{denaturation_state}",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine", "lysine"],
        "molar_ratios": {
            "cysteine": 1200.0,
            "lysine": 800.0,
            "ribose": 1000.0,
        },
        "time_minutes": time_minutes,
        "temp": temp_c,
        "ph": pH,
        "aw": 0.95,
        "protein_type": protein_type,
    }
    if denaturation_state is not None:
        formulation["denaturation_state"] = denaturation_state
    conditions = ReactionConditions(pH=pH, temperature_celsius=temp_c, water_activity=0.95)
    return designer.evaluate_single(formulation, conditions)


def test_matrix_accessibility_penalizes_sulfurous_meaty_prediction_in_recommendation_path():
    result_free = _matrix_accessibility_probe("free", 0.0)
    result_soy = _matrix_accessibility_probe("soy_iso", 0.0)
    result_pea = _matrix_accessibility_probe("pea_iso", 0.0)

    assert result_free.target_score > result_soy.target_score > result_pea.target_score
    assert result_free.predicted_ppb["2-Methyl-3-furanthiol (MFT)"] > result_soy.predicted_ppb["2-Methyl-3-furanthiol (MFT)"] > result_pea.predicted_ppb["2-Methyl-3-furanthiol (MFT)"]
    assert result_free.predicted_ppb["2-Furfurylthiol (FFT)"] > result_soy.predicted_ppb["2-Furfurylthiol (FFT)"] > result_pea.predicted_ppb["2-Furfurylthiol (FFT)"]


def test_pea_denaturation_recovers_signal_without_collapsing_back_to_free_amino_acid_behavior():
    result_free = _matrix_accessibility_probe("free", 0.0)
    result_pea_native = _matrix_accessibility_probe("pea_iso", 0.0)
    result_pea_denatured = _matrix_accessibility_probe("pea_iso", 1.0)

    assert result_pea_denatured.target_score > result_pea_native.target_score
    assert result_pea_denatured.predicted_ppb["2-Methyl-3-furanthiol (MFT)"] > result_pea_native.predicted_ppb["2-Methyl-3-furanthiol (MFT)"]
    assert result_pea_denatured.predicted_ppb["2-Furfurylthiol (FFT)"] > result_pea_native.predicted_ppb["2-Furfurylthiol (FFT)"]

    assert result_pea_denatured.target_score < result_free.target_score
    assert result_pea_denatured.predicted_ppb["2-Methyl-3-furanthiol (MFT)"] < result_free.predicted_ppb["2-Methyl-3-furanthiol (MFT)"]
    assert result_pea_denatured.predicted_ppb["2-Furfurylthiol (FFT)"] < result_free.predicted_ppb["2-Furfurylthiol (FFT)"]


def test_automatic_denaturation_inference_recovers_signal_without_manual_override():
    pea_native = _matrix_accessibility_probe("pea_iso", 0.0, temp_c=90.0, time_minutes=30.0)
    pea_auto = _matrix_accessibility_probe("pea_iso", None, temp_c=90.0, time_minutes=30.0)
    pea_denatured = _matrix_accessibility_probe("pea_iso", 1.0, temp_c=90.0, time_minutes=30.0)

    assert 0.4 < pea_auto.effective_denaturation_state < 0.9
    assert pea_native.effective_denaturation_state == 0.0
    assert pea_denatured.effective_denaturation_state == 1.0
    assert pea_native.target_score < pea_auto.target_score < pea_denatured.target_score
    assert (
        pea_native.predicted_ppb["2-Furfurylthiol (FFT)"]
        < pea_auto.predicted_ppb["2-Furfurylthiol (FFT)"]
        < pea_denatured.predicted_ppb["2-Furfurylthiol (FFT)"]
    )


def test_automatic_denaturation_stays_low_for_cool_pea_matrix_conditions():
    pea_cool = _matrix_accessibility_probe("pea_iso", None, temp_c=40.0, time_minutes=10.0, pH=6.0)
    pea_hot = _matrix_accessibility_probe("pea_iso", None, temp_c=140.0, time_minutes=30.0, pH=5.5)

    assert pea_cool.effective_denaturation_state < 0.2
    assert pea_hot.effective_denaturation_state > 0.9
    assert pea_hot.target_score > pea_cool.target_score