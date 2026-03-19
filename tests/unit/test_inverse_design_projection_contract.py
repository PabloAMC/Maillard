from unittest.mock import MagicMock

import pytest

from src.inverse_design import InverseDesigner
from src.smirks_engine import ReactionConditions, Species
from src.pathway_extractor import Species as OutputSpecies
from src.recommend import _apply_output_projection, _canon


def test_inverse_design_preserves_proxy_and_projection_metadata(monkeypatch):
    designer = InverseDesigner(target_tag="meaty")
    formulation = {
        "name": "Projection contract",
        "protein_type": "free",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine"],
        "molar_ratios": {"ribose": 1.0, "cysteine": 1.0},
    }
    conditions = ReactionConditions(temperature_celsius=150.0, pH=5.0)

    monkeypatch.setattr("src.inverse_design.resolve_many", lambda names: [
        Species("ribose", "O=CC(O)C(O)C(O)CO"),
        Species("cysteine", "NC(CS)C(=O)O"),
    ])

    fake_engine = MagicMock()
    fake_engine.enumerate.return_value = []
    monkeypatch.setattr("src.inverse_design.SmirksEngine", lambda cond: fake_engine)
    monkeypatch.setattr("src.inverse_design.predict_lop_generation", lambda *args, **kwargs: {})
    monkeypatch.setattr(designer.db, "get_best_barrier", lambda *args, **kwargs: (20.0, "mock", 2.0))
    monkeypatch.setattr(designer.sensory, "get_radar_data", lambda *args, **kwargs: {"meaty": (1.0, 1)})
    monkeypatch.setattr(
        "src.inverse_design.evaluate_formulation_safety",
        lambda *args, **kwargs: (0.0, []),
    )

    projection_metadata = {
        "furfural": {
            "proxy_ppb": 120.0,
            "matrix_factor": 0.8,
            "headspace_factor": 0.5,
            "observable_ppb": 48.0,
        }
    }

    fake_recommender = MagicMock()
    fake_recommender.predict_from_steps.return_value = {
        "targets": [
            {
                "name": "furfural",
                "span_uncertainty": 1.5,
                "projection": projection_metadata["furfural"],
            }
        ],
        "metrics": {"trapping_efficiency": {}, "lysine_budget_dha": 0.0},
        "predicted_ppb": {"furfural": 48.0},
        "predicted_proxy_ppb": {"furfural": 120.0},
        "projection_metadata": projection_metadata,
        "projection_context": {"total_volatile_budget_molar": 1.2e-6},
    }
    monkeypatch.setattr("src.inverse_design.Recommender", lambda: fake_recommender)

    result = designer.evaluate_single(formulation, conditions)

    assert result.predicted_ppb["furfural"] == 48.0
    assert result.predicted_proxy_ppb["furfural"] == 120.0
    assert result.projection_metadata["furfural"]["observable_ppb"] == 48.0
    assert result.projection_metadata["furfural"]["proxy_ppb"] == 120.0


def test_output_projection_applies_compound_specific_matrix_calibration():
    species = OutputSpecies("Hexanal", "CCCCCC=O")
    canon = _canon(species.smiles)

    observable, metadata = _apply_output_projection(
        {canon: 100.0},
        {canon: species},
        {},
        120.0 + 273.15,
        protein_type="soy_iso",
        time_minutes=20.0,
        denaturation_state=0.5,
    )

    row = metadata[canon]

    assert row["calibration_observable_factor"] == pytest.approx((0.453 / 0.205) * (1.0 - 0.7060))
    assert row["calibration_factor"] == pytest.approx(row["calibration_observable_factor"])
    assert observable[canon] == pytest.approx(
        100.0 * row["matrix_factor"] * row["headspace_factor"] * row["calibration_factor"]
    )