import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_wave1 import (  # noqa: E402
    build_wave1_dft_ingestion_artifact,
    build_wave1_dft_promotion_payload,
    build_wave1_computational_gap_plan_artifact,
    build_wave1_dft_job_manifest,
    build_wave1_xtb_job_manifest,
    render_wave1_dft_ingestion_markdown,
)


def test_wave1_plan_surfaces_targets_ceilings_and_write_back_artifacts():
    payload = build_wave1_computational_gap_plan_artifact()

    assert payload["summary"]["target_count"] == 4
    by_id = {row["id"]: row for row in payload["targets"]}

    assert by_id["hexanal_radical_quench"]["promotion_ceiling"] == "ranking_only"
    assert "data/lit/computational_priors.json" in by_id["hexanal_radical_quench"]["write_back_artifacts"]
    assert by_id["hexanal_radical_quench"]["dft"]["spin"] == 1
    assert by_id["asparagine_sugar_explicit_water_cluster"]["promotion_ceiling"] == "uncertainty_narrowing_only"
    assert by_id["asparagine_sugar_explicit_water_cluster"]["dft"]["use_explicit_solvent"] is True


def test_wave1_manifests_capture_xtb_and_dft_readiness_without_wet_lab_overclaim():
    xtb_manifest = build_wave1_xtb_job_manifest()
    dft_manifest = build_wave1_dft_job_manifest()

    assert xtb_manifest["summary"]["job_count"] == 4
    assert xtb_manifest["summary"]["ready_count"] >= 3
    assert xtb_manifest["summary"]["seed_required_count"] == 1

    assert dft_manifest["summary"]["job_count"] == 4
    assert dft_manifest["summary"]["seed_required_count"] == 1
    assert any(job["status"] == "blocked_missing_xtb_ts" for job in dft_manifest["jobs"])


def test_dft_ingestion_report_handles_missing_results_file():
    manifest_payload = build_wave1_dft_job_manifest()
    payload = build_wave1_dft_ingestion_artifact(manifest_payload=manifest_payload)
    markdown = render_wave1_dft_ingestion_markdown(payload)

    assert payload["summary"]["available_result_count"] == 0
    assert "No DFT results available yet." in markdown


def test_wave1_dft_promotion_payload_only_promotes_completed_non_fast_results():
    priors_payload = {
        "dft_kinetic_priors": {
            "entries": [
                {
                    "reaction_key": "hexanal_radical_quench",
                    "current_tier": "xtb_derived_gfn2",
                    "promotion_ceiling": "ranking_only",
                },
                {
                    "reaction_key": "lysinoalanine_crosslink",
                    "current_tier": "xtb_derived_gfn2",
                    "promotion_ceiling": "bounded_calibration",
                },
            ]
        }
    }
    execution_payload = {
        "jobs": [
            {
                "reaction_key": "hexanal_radical_quench",
                "status": "completed",
                "fast_mode": False,
                "promotion_ready": True,
                "barrier_kcal_mol": 8.5,
            },
            {
                "reaction_key": "lysinoalanine_crosslink",
                "status": "completed",
                "fast_mode": True,
                "promotion_ready": False,
                "barrier_kcal_mol": 11.2,
            },
        ]
    }

    payload = build_wave1_dft_promotion_payload(
        priors_payload=priors_payload,
        execution_payload=execution_payload,
    )

    assert payload["summary"]["promoted_count"] == 1
    assert payload["promotions"][0]["reaction_key"] == "hexanal_radical_quench"