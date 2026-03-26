import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chemistry_benchmark_validator import (  # noqa: E402
    build_adoption_decisions_from_assessment,
    build_mlp_assessment_artifact,
    evaluate_candidate_against_reaction_benchmark,
    render_mlp_assessment_markdown,
)
from src.mlp_adoption_contract import (  # noqa: E402
    MLPModelCandidate,
    build_adoption_note_payload,
    render_adoption_note_markdown,
)
from src.reaction_benchmark import load_reaction_benchmark_entries  # noqa: E402


def test_validator_quarantines_nonphysical_barrier_predictions(tmp_path):
    benchmark_file = tmp_path / "wild_results.json"
    benchmark_file.write_text(
        json.dumps(
            {
                "strecker": {"estimated_barrier_kcal": 10000.0},
                "enolisation": {"estimated_barrier_kcal": 5000.0},
                "retro_aldol": {"estimated_barrier_kcal": 2500.0},
            }
        ),
        encoding="utf-8",
    )
    candidate = MLPModelCandidate(
        candidate_id="toy_candidate",
        model_family="toy",
        model_name="bad",
        chemistry_domain="general_organic",
        materials_first=False,
        proposed_role="local_barrier_surrogate",
        target_motif_families=["carbonyl_fragmentation", "proton_transfer_dehydration"],
        benchmark_visible_gap="toy gap",
        expected_speedup=10.0,
        likely_failure_modes=["nonphysical_absolute_energies"],
        fallback_comparator="xTB_plus_selective_DFT",
        benchmark_results_path=str(benchmark_file),
        geometry_benchmark_path=None,
        status="results_available",
    )

    row = evaluate_candidate_against_reaction_benchmark(candidate, load_reaction_benchmark_entries())

    assert row["decision"] == "quarantine"
    assert "nonphysical_barrier_predictions" in row["stop_reasons"]
    assert row["nonphysical_reaction_count"] >= 1


def test_validator_defers_geometry_candidate_without_geometry_evidence():
    candidate = MLPModelCandidate(
        candidate_id="geom_only",
        model_family="mace_mp",
        model_name="small",
        chemistry_domain="general_organic",
        materials_first=False,
        proposed_role="geom_preopt",
        target_motif_families=["carbonyl_fragmentation"],
        benchmark_visible_gap="geometry acceleration",
        expected_speedup=5.0,
        likely_failure_modes=["geometry_drift_on_sulfur_species"],
        fallback_comparator="r2SCAN_3c_geometry_plus_wB97MV_single_point",
        benchmark_results_path=None,
        geometry_benchmark_path=None,
        external_evidence_id="aimnet2_external_2026q1",
        status="shortlist_candidate",
    )

    row = evaluate_candidate_against_reaction_benchmark(candidate, load_reaction_benchmark_entries())

    assert row["decision"] == "defer"
    assert "missing_geometry_benchmark_evidence" in row["stop_reasons"]
    assert row["external_prior_strength"] == "strong"
    assert row["external_selection_priority"] == "high"


def test_validator_reads_candidate_indexed_geometry_assessment(tmp_path):
    geometry_file = tmp_path / "geometry_assessment.json"
    geometry_file.write_text(
        json.dumps(
            {
                "candidate_assessments": [
                    {
                        "candidate_id": "geom_only",
                        "available": True,
                        "backend_available": True,
                        "reason": "evaluated",
                        "max_rmsd_angstrom": 0.22,
                        "mean_rmsd_angstrom": 0.12,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    candidate = MLPModelCandidate(
        candidate_id="geom_only",
        model_family="mace_mp",
        model_name="small",
        chemistry_domain="general_organic",
        materials_first=False,
        proposed_role="geom_preopt",
        target_motif_families=["carbonyl_fragmentation"],
        benchmark_visible_gap="geometry acceleration",
        expected_speedup=5.0,
        likely_failure_modes=["geometry_drift_on_sulfur_species"],
        fallback_comparator="r2SCAN_3c_geometry_plus_wB97MV_single_point",
        benchmark_results_path=None,
        geometry_benchmark_path=str(geometry_file),
        external_evidence_id="mace_omol_external_2026q1",
        status="shortlist_candidate",
    )

    row = evaluate_candidate_against_reaction_benchmark(candidate, load_reaction_benchmark_entries())

    assert row["decision"] == "adopt_offline"
    assert row["stop_reasons"] == []


def test_validator_can_build_default_geometry_assessment_when_file_is_not_materialized(monkeypatch, tmp_path):
    from src import chemistry_benchmark_validator as validator_module

    monkeypatch.setattr(
        validator_module,
        "_repo_root",
        lambda: ROOT,
    )

    def _fake_build_geometry_assessment_artifact():
        return {
            "candidate_assessments": [
                {
                    "candidate_id": "mace_mp_small",
                    "available": True,
                    "backend_available": True,
                    "reason": "evaluated",
                    "max_rmsd_angstrom": 0.21,
                    "mean_rmsd_angstrom": 0.12,
                }
            ]
        }

    import src.geometry_benchmark_validator as geometry_validator

    monkeypatch.setattr(
        geometry_validator,
        "build_geometry_assessment_artifact",
        _fake_build_geometry_assessment_artifact,
    )

    candidate = MLPModelCandidate(
        candidate_id="mace_mp_small",
        model_family="mace_mp",
        model_name="small",
        chemistry_domain="general_organic",
        materials_first=False,
        proposed_role="geom_preopt",
        target_motif_families=["carbonyl_fragmentation"],
        benchmark_visible_gap="geometry acceleration",
        expected_speedup=5.0,
        likely_failure_modes=["geometry_drift_on_sulfur_species"],
        fallback_comparator="r2SCAN_3c_geometry_plus_wB97MV_single_point",
        benchmark_results_path=None,
        geometry_benchmark_path=str(tmp_path / "missing" / "mlp_geometry_assessment.json"),
        external_evidence_id="mace_mp_external_2026q1",
        status="shortlist_candidate",
    )

    row = evaluate_candidate_against_reaction_benchmark(candidate, load_reaction_benchmark_entries())

    assert row["decision"] == "adopt_offline"
    assert "missing_geometry_benchmark_evidence" not in row["stop_reasons"]


def test_validator_distinguishes_backend_unavailable_from_missing_ts_seed_evidence(tmp_path):
    ts_seed_file = tmp_path / "ts_seed_assessment.json"
    ts_seed_file.write_text(
        json.dumps(
            {
                "candidate_assessments": [
                    {
                        "candidate_id": "aimnet2_shortlist",
                        "available": False,
                        "backend_available": False,
                        "reason": "AIMNet2 backend not installed",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    candidate = MLPModelCandidate(
        candidate_id="aimnet2_shortlist",
        model_family="aimnet2",
        model_name="aimnet2_wb97m",
        chemistry_domain="organic_qm_surrogate",
        materials_first=False,
        proposed_role="ts_initialization",
        target_motif_families=["proton_transfer_rearrangement"],
        benchmark_visible_gap="ts initialization",
        expected_speedup=5.0,
        likely_failure_modes=["missing_transition_state_support"],
        fallback_comparator="sella_plus_pyscf",
        benchmark_results_path=None,
        geometry_benchmark_path=None,
        ts_seed_benchmark_path=str(ts_seed_file),
        external_evidence_id="aimnet2_external_2026q1",
        status="shortlist_candidate",
    )

    row = evaluate_candidate_against_reaction_benchmark(candidate, load_reaction_benchmark_entries())

    assert row["decision"] == "defer"
    assert "ts_seed_backend_unavailable" in row["stop_reasons"]
    assert "missing_ts_seed_benchmark_evidence" not in row["stop_reasons"]


def test_mlp_assessment_and_adoption_notes_render_current_registry():
    payload = build_mlp_assessment_artifact()
    markdown = render_mlp_assessment_markdown(payload)
    decisions = build_adoption_decisions_from_assessment(payload)
    note_payload = build_adoption_note_payload(
        decisions,
        benchmark_set_id=payload["summary"]["benchmark_set_id"],
    )
    note_markdown = render_adoption_note_markdown(note_payload)

    assert payload["summary"]["candidate_count"] >= 3
    assert any(row["decision"] == "quarantine" for row in payload["candidates"])
    assert "Offline ML Accelerator Assessment" in markdown
    assert "External Prior" in markdown
    assert note_payload["summary"]["default_policy"] == "no_default_mlp_until_reaction_benchmark_passes"
    assert "Offline ML Accelerator Adoption Notes" in note_markdown