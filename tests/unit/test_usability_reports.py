from pathlib import Path

from src.conditions import ReactionConditions
from src.pipeline import FormulationResult
from src.projection_metadata import normalize_projection_metadata_row
from src.reporting import generate_comparison_report, generate_report
from src.usability_reports import (
    DomainWarning,
    assess_formulation_confidence,
    build_confidence_package,
    build_formulation_explainability_payload,
    build_validated_envelope_report,
)
from src.presentation import (
    render_formulation_explainability_markdown,
    render_validated_envelope_markdown,
    render_projection_rows_markdown,
    render_provenance_markdown,
)
from src.projection_utils import build_projection_rows


def test_formulation_explainability_payload_surfaces_matrix_and_projection_context():
    formulation = {
        "name": "Pea explainability probe",
        "protein_type": "pea_iso",
    }
    result = FormulationResult(
        name="Pea explainability probe",
        target_score=2.5,
        off_flavour_risk=0.4,
        safety_score=0.1,
        texture_risk=5.0,
        lysine_budget=12.0,
        trapping_efficiency=7.0,
        predicted_ppb={"furfural": 48.0},
        predicted_proxy_ppb={"furfural": 120.0},
        projection_metadata={
            "furfural": {
                "proxy_ppb": 120.0,
                "observable_ppb": 48.0,
                "proxy_to_observable_ratio": 0.4,
                "matrix_factor": 0.8,
                "dynamic_retention_factor": 1.1,
                "headspace_factor": 0.5,
                "volatile_class": "furan",
                "process_state": "heated_matrix",
                "retention_runtime_mode": "temporal_attenuation",
                "calibration_source": "class_fallback",
                "calibration_evidence_strength": "heuristic",
                "calibration_fallback_mode": "class_level",
            }
        },
        effective_denaturation_state=0.62,
        matrix_explainability={
            "effective_denaturation_state": 0.62,
            "lysine_accessibility": 0.40,
            "cysteine_accessibility": 0.05,
            "bulk_volatile_retention": 0.51,
            "denaturation_source": "test source",
        },
        strecker_balance_score=0.3,
        pyrazine_burden=12.0,
        flavor_axis_summary={
            "strecker_balance_score": 0.3,
            "pyrazine_burden": 12.0,
            "thiamine_pathway_active": False,
        },
        targets=[{"name": "furfural", "concentration": 48.0}],
    )

    payload = build_formulation_explainability_payload(
        formulation,
        result,
        target_tag="meaty",
        minimize_tag="beany",
    )
    markdown = render_formulation_explainability_markdown(payload)

    assert payload["matrix_explainability"]["effective_denaturation_state"] == 0.62
    assert payload["top_projection_rows"][0]["volatile_class"] == "furan"
    assert payload["top_projection_rows"][0]["observable_ratio"] == 0.4
    assert payload["top_projection_rows"][0]["process_state"] == "heated_matrix"
    assert payload["top_projection_rows"][0]["retention_runtime_mode"] == "temporal_attenuation"
    assert "Formulation Explainability" in markdown
    assert "Effective denaturation state" in markdown
    assert "Obs/Proxy" in markdown
    assert "Calibration" in markdown
    assert "Flavor Axis" in markdown


def test_projection_metadata_normalizer_fills_schema_defaults_for_sparse_rows():
    row = normalize_projection_metadata_row(
        {
            "proxy_ppb": 120.0,
            "observable_ppb": 48.0,
            "process_state": "heated_matrix",
        },
        compound_fallback="furfural",
    )

    assert row["compound"] == "furfural"
    assert row["matrix_factor"] == 1.0
    assert row["headspace_factor"] == 1.0
    assert row["retention_runtime_mode"] == "static_class_profile"
    assert row["calibration_source"] == "class_fallback"
    assert row["proxy_to_observable_ratio"] == 0.4


def test_build_projection_rows_resolves_target_projection_when_metadata_key_is_canonical():
    result = FormulationResult(
        name="canonical-key probe",
        target_score=1.0,
        off_flavour_risk=0.0,
        targets=[
            {
                "name": "2-Methyl-3-furanthiol (MFT)",
                "concentration": 18.0,
                "projection": {
                    "compound": "2-Methyl-3-furanthiol (MFT)",
                    "proxy_ppb": 30.0,
                    "observable_ppb": 18.0,
                    "proxy_to_observable_ratio": 0.6,
                    "matrix_factor": 0.75,
                    "headspace_factor": 0.8,
                    "process_state": "heated_matrix",
                    "retention_runtime_mode": "dynamic_release",
                    "calibration_source": "compound_specific",
                    "calibration_evidence_strength": "literature_anchored",
                    "calibration_fallback_mode": "compound_specific",
                },
            }
        ],
        projection_metadata={
            "c1oc(cc1)CS": {
                "compound": "2-Methyl-3-furanthiol (MFT)",
                "proxy_ppb": 30.0,
                "observable_ppb": 18.0,
            }
        },
    )

    rows = build_projection_rows(result)

    assert len(rows) == 1
    assert rows[0]["compound"] == "2-Methyl-3-furanthiol (MFT)"
    assert rows[0]["observable_ppb"] == 18.0
    assert rows[0]["observable_ratio"] == 0.6
    assert rows[0]["retention_runtime_mode"] == "dynamic_release"


def test_validated_envelope_report_mentions_strict_ready_and_matrix_scope():
    report = build_validated_envelope_report(target_tag="meaty")
    markdown = render_validated_envelope_markdown(report)

    assert report.supported_benchmarks >= 1
    assert any("Matrix benchmarks" in warning for warning in report.warnings)
    assert "Validated Envelope" in markdown
    assert "Strict-ready benchmarks" in markdown
    assert "Matrix-only executable benchmarks" in markdown


def test_assess_formulation_confidence_marks_primary_free_precursor_as_high():
    result = FormulationResult(
        name="trusted free run",
        target_score=10.0,
        off_flavour_risk=0.0,
        avg_uncertainty=2.2,
        matrix_explainability={"protein_type": "free"},
    )

    assessment = assess_formulation_confidence(
        result,
        [],
        precursor_names=["glucose", "cysteine"],
        protein_type="free",
    )

    assert assessment.tier == "high"
    assert assessment.benchmark_neighborhood == "primary_free_precursor"
    assert assessment.score >= 85.0


def test_assess_formulation_confidence_marks_matrix_sparse_case_as_exploratory():
    result = FormulationResult(
        name="matrix exploratory run",
        target_score=4.0,
        off_flavour_risk=1.0,
        avg_uncertainty=8.7,
        matrix_explainability={"protein_type": "pea_iso"},
    )
    warnings = [
        DomainWarning(level="WARNING", category="PRECURSORS", message="Sparse analogies."),
        DomainWarning(level="CAUTION", category="MATRIX", message="Matrix-only support."),
    ]

    assessment = assess_formulation_confidence(
        result,
        warnings,
        precursor_names=["thiamine", "pea peptide hydrolysate"],
        protein_type="pea_iso",
    )

    assert assessment.tier in {"low", "exploratory"}
    assert assessment.benchmark_neighborhood == "matrix_intake_only"
    assert assessment.score < 65.0


def test_generate_report_includes_confidence_metadata(tmp_path: Path):
    result = FormulationResult(
        name="report confidence probe",
        target_score=12.0,
        off_flavour_risk=0.0,
        matrix_explainability={"protein_type": "free", "effective_denaturation_state": 1.0},
        confidence_metadata={
            "tier": "medium",
            "score": 72.0,
            "benchmark_neighborhood": "free_precursor_partial_analogy",
            "prediction_mode": "ranking_supported",
            "recommended_posture": "Reliable for ranking.",
            "dominant_factors": ["Only part of the precursor set is benchmark-anchored."],
            "calibration_diagnostics": {
                "supported_envelope": False,
                "summary": "Recommendation extrapolates beyond the strongest support.",
                "extrapolation_axes": ["benchmark_neighborhood"],
            },
            "compound_confidence": [
                {
                    "compound": "furfural",
                    "observable_ppb": 12.0,
                    "tier": "medium",
                    "score": 68.0,
                    "prediction_mode": "ranking_supported",
                }
            ],
            "aggregate_confidence": {
                "meaty": {
                    "score": 4.2,
                    "support_count": 2,
                    "tier": "medium",
                    "prediction_mode": "ranking_supported",
                }
            },
            "sensitivity_summary": {
                "mode": "local_oat",
                "evaluated_perturbations": 4,
                "ranking_drivers": [
                    {
                        "input": "cysteine",
                        "perturbation": "remove_amino_acids",
                        "decision_delta": -2.0,
                        "safety_delta": 0.0,
                    }
                ],
                "safety_drivers": [
                    {
                        "input": "temp",
                        "perturbation": "increase_temp",
                        "safety_delta": 1.5,
                    }
                ],
            },
        },
        projection_metadata={
            "furfural": {
                "compound": "furfural",
                "proxy_ppb": 18.0,
                "observable_ppb": 12.0,
                "proxy_to_observable_ratio": 0.67,
                "matrix_factor": 0.85,
                "headspace_factor": 0.7,
                "melanoidin_trapping_factor": 1.0,
                "process_state": "ambient_slurry",
                "retention_runtime_mode": "static_class_profile",
                "calibration_source": "Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
                "calibration_evidence_strength": "literature_anchored",
                "calibration_fallback_mode": "compound_specific",
            }
        },
        strecker_balance_score=0.25,
        strecker_gap_penalty=0.4,
        pyrazine_propensity=0.7,
        pyrazine_burden=14.0,
        pyrazine_penalty=0.2,
        flavor_axis_summary={
            "strecker_balance_score": 0.25,
            "strecker_gap_penalty": 0.4,
            "pyrazine_signal_ppb": 7.0,
            "pyrazine_propensity": 0.7,
            "pyrazine_burden": 14.0,
            "pyrazine_penalty": 0.2,
            "thiamine_pathway_active": False,
            "thiamine_availability_source": "native_matrix_default_inactive",
            "thiamine_provenance_mode": "inactive",
            "lincoln_crosstalk_prior": {"summary": "inactive"},
        },
    )

    out_dir = generate_report(result, [], {"sugars": "glucose"}, output_dir=tmp_path / "report")
    json_text = (out_dir / "report.json").read_text()
    markdown_text = (out_dir / "report.md").read_text()

    assert "confidence_metadata" in json_text
    assert "Confidence & Support" in markdown_text
    assert "free_precursor_partial_analogy" in markdown_text
    assert "Calibration Diagnostics" in markdown_text
    assert "Compound Confidence" in markdown_text
    assert "Aggregate Sensory Confidence" in markdown_text
    assert "Sensitivity Summary" in markdown_text
    assert "Benchmark Neighborhood" in markdown_text
    assert "Calibration Summary" in markdown_text
    assert "Compound Evidence Ladder" in markdown_text
    assert "Missing Data" in markdown_text
    assert "Safety Reference Context" in markdown_text
    assert "Flavor Reference Policy" in markdown_text
    assert "Projection Calibration" in markdown_text
    assert "Flavor Axis Diagnostics" in markdown_text
    assert "projection_metadata" in json_text
    assert "compound_evidence_ladder" in json_text
    assert "calibration_summary" in json_text
    assert "missing_data_summary" in json_text
    assert "benchmark_neighborhood_summary" in json_text
    assert "safety_reference_summary" in json_text
    assert "flavor_reference_policy" in json_text
    assert '"provenance"' in json_text
    assert "## 5. Provenance" in markdown_text


def test_generate_comparison_report_includes_provenance(tmp_path: Path):
    first = FormulationResult(
        name="A",
        target_score=5.0,
        off_flavour_risk=1.0,
        safety_score=0.3,
        mft_to_furfural_ratio=0.02,
        meaty_quality_penalty=0.4,
        strecker_balance_score=0.2,
        strecker_gap_penalty=0.8,
        pyrazine_burden=35.0,
        pyrazine_penalty=0.7,
        projection_metadata={
            "mft": {
                "compound": "2-Methyl-3-furanthiol (MFT)",
                "observable_ppb": 8.0,
                "melanoidin_trapping_factor": 0.52,
                "volatile_class": "sulfur",
                "calibration_source": "Pratap-Singh 2021 soy-vs-pea ambient slurry release ratio",
                "calibration_evidence_strength": "literature_anchored",
                "calibration_fallback_mode": "compound_specific",
            }
        },
        confidence_metadata={"benchmark_neighborhood": "matrix_intake_only", "prediction_mode": "directional_only"},
        flavor_axis_summary={"thiamine_pathway_active": False, "thiamine_availability_source": "native_matrix_default_inactive", "furanone_expected": ["HEMF"]},
    )
    second = FormulationResult(
        name="B",
        target_score=3.0,
        off_flavour_risk=0.5,
        safety_score=0.1,
        mft_to_furfural_ratio=0.10,
        meaty_quality_penalty=0.1,
        strecker_balance_score=0.7,
        strecker_gap_penalty=0.0,
        pyrazine_burden=8.0,
        pyrazine_penalty=0.0,
        projection_metadata={
            "fft": {
                "compound": "2-Furfurylthiol (FFT)",
                "observable_ppb": 11.0,
                "melanoidin_trapping_factor": 0.91,
                "volatile_class": "sulfur",
                "calibration_source": "class_fallback",
                "calibration_evidence_strength": "heuristic",
                "calibration_fallback_mode": "class_level",
            }
        },
        confidence_metadata={"benchmark_neighborhood": "free_precursor_partial_analogy", "prediction_mode": "ranking_supported"},
        flavor_axis_summary={"thiamine_pathway_active": True, "thiamine_availability_source": "pbma_fortified", "furanone_expected": ["HEMF", "DMHF"], "furanone_missing": ["DMHF"]},
    )

    out_dir = generate_comparison_report(
        [first, second],
        [{"name": "A", "ph": 5.5, "temp": 105.0}, {"name": "B", "ph": 5.8, "temp": 120.0}],
        output_dir=tmp_path / "comparison",
        campaign_metadata={"name": "shareable-screen"},
    )
    json_text = (out_dir / "comparison.json").read_text()
    markdown_text = (out_dir / "comparison.md").read_text()

    assert '"provenance"' in json_text
    assert '"campaign"' in json_text
    assert '"mft_to_furfural_ratio"' in json_text
    assert '"strecker_balance_score"' in json_text
    assert '"pyrazine_burden"' in json_text
    assert '"sulfur_trapping_summary"' in json_text
    assert '"benchmark_neighborhood_summary"' in json_text
    assert '"safety_reference_summary"' in json_text
    assert '"flavor_reference_policy"' in json_text
    assert "Cross-Marker Context" in markdown_text
    assert "Calibration Contrast" in markdown_text
    assert "MFT/Furfural Ratio" in markdown_text
    assert "Sulfur Trapping" in markdown_text
    assert "Strecker Support" in markdown_text
    assert "Benchmark Neighborhood" in markdown_text
    assert "## 5. Provenance" in markdown_text


def test_build_confidence_package_adds_compound_aggregate_and_sensitivity_sections():
    class StubDesigner:
        def evaluate_single(self, formulation, conditions):
            target_shift = 4.0 if "cysteine" in formulation.get("amino_acids", []) else 0.8
            safety_shift = 0.0
            if float(formulation.get("temp", conditions.temperature_celsius)) > 150.0:
                safety_shift = 1.5
            return FormulationResult(
                name="variant",
                target_score=target_shift,
                off_flavour_risk=0.5,
                safety_score=safety_shift,
                texture_risk=2.0,
            )

    result = FormulationResult(
        name="baseline",
        target_score=4.0,
        off_flavour_risk=0.5,
        safety_score=0.2,
        texture_risk=2.0,
        avg_uncertainty=4.0,
        radar={"meaty": (4.0, 2), "beany": (0.6, 1)},
        targets=[
            {
                "name": "furfural",
                "span_uncertainty": 4.5,
            }
        ],
        projection_metadata={
            "furfural": {
                "compound": "furfural",
                "proxy_ppb": 40.0,
                "observable_ppb": 20.0,
                "proxy_to_observable_ratio": 0.5,
                "matrix_factor": 0.8,
                "headspace_factor": 0.6,
            }
        },
        precursor_contributions={"cysteine": 20.0, "glucose": 10.0},
    )
    formulation = {
        "name": "baseline",
        "sugars": ["glucose"],
        "amino_acids": ["cysteine"],
        "additives": [],
        "lipids": [],
        "protein_type": "free",
        "ph": 6.0,
        "temp": 150.0,
        "aw": 0.8,
        "time_minutes": 60.0,
    }

    payload = build_confidence_package(
        result,
        [],
        precursor_names=["glucose", "cysteine"],
        protein_type="free",
        formulation=formulation,
        baseline_conditions=ReactionConditions(pH=6.0, temperature_celsius=150.0, water_activity=0.8),
        designer=StubDesigner(),
    )

    assert payload["prediction_mode"] == "benchmark_supported_quantitative"
    assert payload["compound_confidence"][0]["compound"] == "furfural"
    assert payload["aggregate_confidence"]["meaty"]["tier"] in {"medium", "high"}
    assert payload["sensitivity_summary"]["mode"] == "local_oat"
    assert payload["sensitivity_summary"]["evaluated_perturbations"] >= 4
    assert payload["sensitivity_summary"]["ranking_drivers"][0]["input"] in {"cysteine", "temp", "ph", "time_minutes"}