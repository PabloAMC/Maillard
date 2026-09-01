from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.data_ingest import build_ingest_preview, build_intake_payload, write_ingest_artifacts
from src.matrix_experiment_intake import build_matrix_experiment_support_delta_artifact, load_matrix_experiment_intake


ROOT = Path(__file__).resolve().parents[2]


def test_build_intake_payload_from_csv_matches_expected_contract(tmp_path: Path):
    csv_path = tmp_path / "scientist_export.csv"
    frame = pd.DataFrame(
        [
            {
                "Analyte Name": "2-furfurylthiol",
                "Observed ppb": 0.0056,
                "RSD %": 12.0,
                "Temperature C": 100.0,
                "pH": 5.6,
                "Aw": 0.95,
                "Time Minutes": 45.0,
                "Reference": "Example pea isolate external matrix package",
                "DOI": "10.0000/example-pea-matrix-package",
                "Measurement Date": "2026-04-01",
                "Notes": "Imported from GC-MS export",
                "Target Benchmark ID": "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
            },
            {
                "Analyte Name": "2-methyl-3-furanthiol",
                "Observed ppb": 0.0037,
                "RSD %": 12.0,
                "Temperature C": 100.0,
                "pH": 5.6,
                "Aw": 0.95,
                "Time Minutes": 45.0,
                "Reference": "Example pea isolate external matrix package",
                "DOI": "10.0000/example-pea-matrix-package",
                "Measurement Date": "2026-04-01",
                "Notes": "Imported from GC-MS export",
                "Target Benchmark ID": "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
            },
        ]
    )
    frame.to_csv(csv_path, index=False)

    payload = build_intake_payload(
        csv_path,
        protein_type="pea_iso",
        process_state="aqueous_pre_extrusion_model",
        experiment_id="pea_iso_csv_roundtrip",
        source_kind="external_literature",
        evidence_class="calibration_candidate",
        matrix_format="pea isolate slurry with exogenous ribose and cysteine",
        precursor_entries=[
            "Pea Protein Isolate=1000.0",
            "D-Ribose=1.0",
            "L-Cysteine=1.0",
        ],
    )

    assert payload["experiment_id"] == "pea_iso_csv_roundtrip"
    assert payload["protein_type"] == "pea_iso"
    assert payload["process_state"] == "aqueous_pre_extrusion_model"
    assert payload["conditions"] == {
        "temp_C": 100.0,
        "ph": 5.6,
        "water_activity": 0.95,
        "time_min": 45.0,
    }
    assert payload["formulation"]["precursors"] == {
        "Pea Protein Isolate": {"concentration_mM": 1000.0},
        "D-Ribose": {"concentration_mM": 1.0},
        "L-Cysteine": {"concentration_mM": 1.0},
    }
    assert payload["measured_volatiles"] == {
        "2-furfurylthiol": {"conc_ppb": 0.0056, "uncertainty_pct": 12.0},
        "2-methyl-3-furanthiol": {"conc_ppb": 0.0037, "uncertainty_pct": 12.0},
    }
    assert payload["provenance"]["source_reference"] == "Example pea isolate external matrix package"
    assert payload["benchmark_alignment"]["target_benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_Internal2026"


def test_build_ingest_preview_surfaces_coverage_and_shift_deltas(monkeypatch):
    payload = load_matrix_experiment_intake(ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml")
    support_delta = build_matrix_experiment_support_delta_artifact(payload)
    aligned_benchmark_id = support_delta["aligned_benchmark"]["benchmark_id"]
    candidate_benchmark_id = payload["experiment_id"]

    baseline_uncertainty = {
        "summary": {
            "ci_coverage_hits": 39,
            "matched_compound_count": 48,
            "ci_coverage_rate": 39 / 48,
        },
        "benchmarks": [
            {
                "benchmark_id": aligned_benchmark_id,
                "compounds": [
                    {
                        "compound": "2-furfurylthiol",
                        "predicted_p5": 1.0,
                        "predicted_p50": 10.0,
                        "predicted_p95": 100.0,
                    }
                ],
            }
        ],
    }
    candidate_uncertainty = {
        "summary": {
            "ci_coverage_hits": 40,
            "matched_compound_count": 49,
            "ci_coverage_rate": 40 / 49,
        },
        "benchmarks": [
            {
                "benchmark_id": aligned_benchmark_id,
                "compounds": baseline_uncertainty["benchmarks"][0]["compounds"],
            },
            {
                "benchmark_id": candidate_benchmark_id,
                "compounds": [
                    {
                        "compound": "2-furfurylthiol",
                        "predicted_p5": 2.0,
                        "predicted_p50": 30.0,
                        "predicted_p95": 60.0,
                    }
                ],
            },
        ],
    }

    calls = {"count": 0}

    def _fake_propagate(*args, **kwargs):
        calls["count"] += 1
        return baseline_uncertainty if calls["count"] == 1 else candidate_uncertainty

    monkeypatch.setattr("src.data_ingest.propagate_benchmarks", _fake_propagate)

    preview = build_ingest_preview(payload, preview_samples=8, preview_seed=3)

    assert preview["coverage_preview"]["participates_in_trust_loop"] is True
    assert preview["coverage_preview"]["delta"]["inside_ci_count"] == 1
    assert preview["coverage_preview"]["delta"]["matched_compound_count"] == 1
    assert preview["compound_deltas"]["top_envelope_tightening_vs_aligned"][0]["compound"] == "2-furfurylthiol"
    assert preview["compound_deltas"]["top_median_shift_vs_aligned"][0]["compound"] == "2-furfurylthiol"


def test_write_ingest_artifacts_is_idempotent(tmp_path: Path):
    payload = load_matrix_experiment_intake(ROOT / "data" / "protocols" / "example_matrix_experiment_intake.yaml")
    preview = {
        "experiment_id": payload["experiment_id"],
        "protein_type": payload["protein_type"],
        "process_state": payload["process_state"],
        "benchmark_preview": {
            "candidate_benchmark_id": payload["experiment_id"],
            "benchmarks_added": [payload["experiment_id"]],
            "aligned_benchmark_id": "pea_isolate_ribose_cysteine_100C_45min_Internal2026",
        },
        "coverage_preview": {"participates_in_trust_loop": False, "reason": "synthetic test"},
        "promotion_preview": {"readiness_change": "evidence_strengthened_not_yet_promotable"},
        "compound_deltas": {},
        "support_delta": {},
        "audit": {},
    }

    first = write_ingest_artifacts(payload, preview, output_dir=tmp_path, persist_intake=True)
    second = write_ingest_artifacts(payload, preview, output_dir=tmp_path, persist_intake=True)

    assert first.intake_yaml == second.intake_yaml
    assert first.preview_json.read_text(encoding="utf-8") == second.preview_json.read_text(encoding="utf-8")
    assert first.preview_md.read_text(encoding="utf-8") == second.preview_md.read_text(encoding="utf-8")
    assert first.support_delta_json.read_text(encoding="utf-8") == second.support_delta_json.read_text(encoding="utf-8")
    assert first.support_delta_md.read_text(encoding="utf-8") == second.support_delta_md.read_text(encoding="utf-8")
    yaml_payload = first.intake_yaml
    assert yaml_payload is not None
    normalized = load_matrix_experiment_intake(yaml_payload)
    assert normalized["experiment_id"] == payload["experiment_id"]


def test_lane_h_csv_templates_round_trip_through_build_intake_payload():
    """Lane H (sprint 2026-05-10b): the bundled CSV templates must parse
    cleanly so a scientist can run the 30-minute walkthrough end-to-end."""
    template_dir = ROOT / "docs" / "templates"
    for filename, expected_compounds in (
        ("hs_spme_gc_ms_template.csv", {"hexanal", "nonanal", "2-pentylfuran"}),
        ("sida_template.csv", {"2-methyl-3-furanthiol"}),
    ):
        template_path = template_dir / filename
        assert template_path.exists(), f"missing template: {template_path}"
        payload = build_intake_payload(
            template_path,
            protein_type="pea_iso",
            process_state="heated_matrix",
            experiment_id=f"template_round_trip_{template_path.stem}",
            source_kind="external_literature",
            evidence_class="calibration_candidate",
            matrix_format="pea isolate template ingest round trip",
            precursor_entries=[
                "Pea Protein Isolate=1000.0",
                "D-Ribose=1.0",
                "L-Cysteine=1.0",
            ],
        )
        assert payload["protein_type"] == "pea_iso"
        ingested = {name.lower() for name in payload["measured_volatiles"].keys()}
        # Template compounds should all be picked up via fuzzy header matching.
        assert expected_compounds.issubset(ingested), (
            f"{filename}: expected {expected_compounds} ingested, got {ingested}"
        )
