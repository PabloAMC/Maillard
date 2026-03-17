from src.inverse_design import FormulationResult
from src.usability_reports import (
    build_formulation_explainability_payload,
    build_validated_envelope_report,
    render_formulation_explainability_markdown,
    render_validated_envelope_markdown,
)


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
                "headspace_factor": 0.5,
                "volatile_class": "furan",
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
    assert "Formulation Explainability" in markdown
    assert "Effective denaturation state" in markdown
    assert "Obs/Proxy" in markdown


def test_validated_envelope_report_mentions_strict_ready_and_matrix_scope():
    report = build_validated_envelope_report(target_tag="meaty")
    markdown = render_validated_envelope_markdown(report)

    assert report.supported_benchmarks >= 1
    assert any("Matrix benchmarks" in warning for warning in report.warnings)
    assert "Validated Envelope" in markdown
    assert "Strict-ready benchmarks" in markdown
    assert "Matrix-only executable benchmarks" in markdown