from unittest.mock import MagicMock

import pytest

from src.pipeline import FormulationResult, MaillardPipeline
from src.projection_utils import build_projection_rows
from src.smirks_engine import ReactionConditions, Species
from src.pathway_extractor import Species as OutputSpecies
from src.recommend import _apply_output_projection, _canon
from src.safety import SafetyResult


def test_pipeline_preserves_proxy_and_projection_metadata(monkeypatch):
    designer = MaillardPipeline(target_tag="meaty")
    formulation = {
        "name": "Projection contract",
        "protein_type": "free",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine"],
        "molar_ratios": {"ribose": 1.0, "cysteine": 1.0},
    }
    conditions = ReactionConditions(temperature_celsius=150.0, pH=5.0)

    monkeypatch.setattr("src.pipeline.resolve_many", lambda names: [
        Species("ribose", "O=CC(O)C(O)C(O)CO"),
        Species("cysteine", "NC(CS)C(=O)O"),
    ])

    fake_engine = MagicMock()
    fake_engine.enumerate.return_value = []
    monkeypatch.setattr("src.pipeline.SmirksEngine", lambda cond: fake_engine)
    monkeypatch.setattr("src.pipeline.predict_lop_generation", lambda *args, **kwargs: {})
    monkeypatch.setattr(designer.db, "get_best_barrier", lambda *args, **kwargs: (20.0, "mock", 2.0))
    monkeypatch.setattr(designer.sensory, "get_radar_data", lambda *args, **kwargs: {"meaty": (1.0, 1)})
    monkeypatch.setattr(
        "src.pipeline.evaluate_formulation_safety",
        lambda *args, **kwargs: (0.0, [], SafetyResult(acrylamide_ppb=0.0, flagged=False, description="Mock")),
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
    monkeypatch.setattr("src.pipeline.Recommender", lambda: fake_recommender)

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
    assert row["accessibility_profile"] == "peptide_bound"
    assert row["accessibility_warning"] is True
    assert row["accessibility_dominant_source"] == "denaturation_state_arg"
    assert observable[canon] == pytest.approx(
        100.0 * row["matrix_factor"] * row["headspace_factor"] * row["calibration_factor"]
    )


def test_pipeline_uses_explicit_thiamine_availability_metadata(monkeypatch):
    designer = MaillardPipeline(target_tag="meaty")
    formulation = {
        "name": "Thiamine provenance contract",
        "protein_type": "soy_iso",
        "sugars": ["ribose"],
        "amino_acids": ["cysteine"],
        "molar_ratios": {"ribose": 1.0, "cysteine": 1.0},
        "thiamine_availability": {"available": True, "source": "pbma_fortified"},
    }
    conditions = ReactionConditions(temperature_celsius=150.0, pH=5.0, protein_type="soy_iso")

    monkeypatch.setattr("src.pipeline.resolve_many", lambda names: [
        Species("ribose", "O=CC(O)C(O)C(O)CO"),
        Species("cysteine", "NC(CS)C(=O)O"),
    ])

    fake_engine = MagicMock()
    fake_engine.enumerate.return_value = []
    monkeypatch.setattr("src.pipeline.SmirksEngine", lambda cond: fake_engine)
    monkeypatch.setattr("src.pipeline.predict_lop_generation", lambda *args, **kwargs: {})
    monkeypatch.setattr(designer.db, "get_best_barrier", lambda *args, **kwargs: (20.0, "mock", 2.0))
    monkeypatch.setattr(designer.sensory, "get_radar_data", lambda *args, **kwargs: {"meaty": (1.0, 1), "beany": (0.0, 1)})
    monkeypatch.setattr(
        "src.pipeline.evaluate_formulation_safety",
        lambda *args, **kwargs: (0.0, [], SafetyResult(acrylamide_ppb=0.0, flagged=False, description="Mock")),
    )

    fake_recommender = MagicMock()
    fake_recommender.predict_from_steps.return_value = {
        "targets": [],
        "metrics": {"trapping_efficiency": {}, "lysine_budget_dha": 0.0},
        "predicted_ppb": {},
        "predicted_proxy_ppb": {},
        "projection_metadata": {},
        "projection_context": {"total_volatile_budget_molar": 0.0},
    }
    monkeypatch.setattr("src.pipeline.Recommender", lambda: fake_recommender)

    result = designer.evaluate_single(formulation, conditions)

    assert result.flavor_axis_summary["thiamine_pathway_active"] is True
    assert result.flavor_axis_summary["thiamine_provenance_mode"] == "mixed_thiamine_plus_pentose"
    assert result.flavor_axis_summary["thiamine_availability_source"] == "pbma_fortified"


def test_pipeline_applies_family_upstream_contract_before_prediction(monkeypatch):
    designer = MaillardPipeline(target_tag="meaty")
    formulation = {
        "name": "Family upstream contract",
        "protein_type": "soy_iso",
        "sugars": ["ribose", "glucose"],
        "amino_acids": ["cysteine"],
        "molar_ratios": {"ribose": 1.0, "glucose": 1.0, "cysteine": 0.5},
        "interventions": ["yeast_fermentation", "protease_hydrolysis"],
        "thiamine_availability": {"available": True, "source": "pbma_fortified"},
    }
    conditions = ReactionConditions(temperature_celsius=150.0, pH=5.5, protein_type="soy_iso")

    captured_names = []
    captured_engine_ph = {}

    def fake_resolve_many(names):
        captured_names[:] = list(names)
        mapping = {
            "ribose": Species("ribose", "O=CC(O)C(O)C(O)CO"),
            "glucose": Species("glucose", "O=CC(O)C(O)C(O)C(O)CO"),
            "cysteine": Species("cysteine", "NC(CS)C(=O)O"),
            "thiamine": Species("thiamine", "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"),
        }
        return [mapping[name] for name in names]

    monkeypatch.setattr("src.pipeline.resolve_many", fake_resolve_many)

    class FakeEngine:
        def __init__(self, cond):
            captured_engine_ph["value"] = cond.pH

        def enumerate(self, precursors, max_generations=4):
            return []

    monkeypatch.setattr("src.pipeline.SmirksEngine", FakeEngine)
    monkeypatch.setattr("src.pipeline.predict_lop_generation", lambda *args, **kwargs: {})
    monkeypatch.setattr(designer.db, "get_best_barrier", lambda *args, **kwargs: (20.0, "mock", 2.0))
    monkeypatch.setattr(designer.sensory, "get_radar_data", lambda *args, **kwargs: {"meaty": (1.0, 1), "beany": (0.0, 1)})
    monkeypatch.setattr(
        "src.pipeline.evaluate_formulation_safety",
        lambda *args, **kwargs: (0.0, [], SafetyResult(acrylamide_ppb=0.0, flagged=False, description="Mock")),
    )

    fake_recommender = MagicMock()
    fake_recommender.predict_from_steps.return_value = {
        "targets": [],
        "metrics": {"trapping_efficiency": {}, "lysine_budget_dha": 0.0},
        "predicted_ppb": {},
        "predicted_proxy_ppb": {},
        "projection_metadata": {},
        "projection_context": {"total_volatile_budget_molar": 0.0},
    }
    monkeypatch.setattr("src.pipeline.Recommender", lambda: fake_recommender)

    result = designer.evaluate_single(formulation, conditions)

    assert "thiamine" in captured_names
    assert captured_engine_ph["value"] == pytest.approx(5.15)

    initial_concentrations = fake_recommender.predict_from_steps.call_args.args[2]
    ribose_smiles = "O=CC(O)C(O)C(O)CO"
    glucose_smiles = "O=CC(O)C(O)C(O)C(O)CO"
    cysteine_smiles = "NC(CS)C(=O)O"
    thiamine_smiles = "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"

    assert initial_concentrations[ribose_smiles] > initial_concentrations[glucose_smiles]
    assert initial_concentrations[cysteine_smiles] == pytest.approx(1.0)
    assert initial_concentrations[thiamine_smiles] > 0.0
    assert result.flavor_axis_summary["family_upstream_contract"]["effective_pH"] == pytest.approx(5.15)


def test_projection_rows_preserve_panel_role_kind_and_modeling_regimes():
    result = FormulationResult(
        name="panel metadata probe",
        target_score=1.0,
        off_flavour_risk=0.0,
        targets=[{"name": "Hexanal", "concentration": 12.0}],
        projection_metadata={
            "hex": {
                "compound": "Hexanal",
                "proxy_ppb": 18.0,
                "observable_ppb": 12.0,
                "target_class": "adverse_lipid_markers",
                "panel_role": "constrained",
                "observable_kind": "volatile",
                "modeling_regimes": ["matrix_hydrated", "extrusion_like"],
            }
        },
    )

    rows = build_projection_rows(result)

    assert rows[0]["panel_role"] == "constrained"
    assert rows[0]["observable_kind"] == "volatile"
    assert rows[0]["modeling_regimes"] == ["matrix_hydrated", "extrusion_like"]