import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.computational_gap_refinement as refinement  # noqa: E402

from src.computational_gap_refinement import (  # noqa: E402
    build_computational_gap_dft_ingestion_artifact,
    build_computational_gap_dft_job_manifest,
    build_computational_gap_dft_promotion_payload,
    build_computational_gap_refinement_plan_artifact,
    build_computational_gap_xtb_job_manifest,
    render_computational_gap_dft_ingestion_markdown,
)


@pytest.fixture
def mocked_geometry_inputs(monkeypatch: pytest.MonkeyPatch) -> None:
    available_paths = {
        "data/geometries/xtb_inputs/hexanal_radical_quench/reactant.xyz",
        "data/geometries/xtb_inputs/hexanal_radical_quench/product.xyz",
        "data/geometries/xtb_inputs/hexanal_radical_quench/run_xtb.sh",
        "data/geometries/xtb_inputs/hexanal_radical_quench/xtbpath_ts.xyz",
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/reactant.xyz",
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/product.xyz",
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/run_xtb.sh",
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/xtbpath_ts.xyz",
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/reactant.xyz",
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/product.xyz",
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/run_xtb.sh",
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/xtbpath_ts.xyz",
        "data/geometries/xtb_inputs/quinone_cys_michael/reactant.xyz",
        "data/geometries/xtb_inputs/quinone_cys_michael/product.xyz",
        "data/geometries/xtb_inputs/quinone_cys_michael/run_xtb.sh",
        "data/geometries/xtb_inputs/quinone_cys_michael/xtbpath_ts.xyz",
        "data/geometries/xtb_inputs/pe_schiff_base/reactant.xyz",
        "data/geometries/xtb_inputs/pe_schiff_base/product.xyz",
        "data/geometries/xtb_inputs/pe_schiff_base/run_xtb.sh",
        "data/geometries/xtb_inputs/pe_schiff_base/xtbpath_ts.xyz",
        "data/geometries/xtb_inputs/pe_amadori/reactant.xyz",
        "data/geometries/xtb_inputs/pe_amadori/product.xyz",
        "data/geometries/xtb_inputs/pe_amadori/run_xtb.sh",
        "data/geometries/xtb_inputs/pe_amadori/xtbpath_ts.xyz",
    }
    atom_counts = {
        "data/geometries/xtb_inputs/hexanal_radical_quench/reactant.xyz": 21,
        "data/geometries/xtb_inputs/hexanal_radical_quench/product.xyz": 21,
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/reactant.xyz": 35,
        "data/geometries/xtb_inputs/lysinoalanine_crosslink/product.xyz": 35,
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/reactant.xyz": 21,
        "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl/product.xyz": 21,
        "data/geometries/xtb_inputs/quinone_cys_michael/reactant.xyz": 26,
        "data/geometries/xtb_inputs/quinone_cys_michael/product.xyz": 26,
        "data/geometries/xtb_inputs/pe_schiff_base/reactant.xyz": 36,
        "data/geometries/xtb_inputs/pe_schiff_base/product.xyz": 36,
        "data/geometries/xtb_inputs/pe_amadori/reactant.xyz": 33,
        "data/geometries/xtb_inputs/pe_amadori/product.xyz": 33,
    }

    monkeypatch.setattr(refinement, "_path_exists", lambda relative_path: relative_path in available_paths)

    def fake_xyz_atom_count(path: Path) -> int | None:
        relative = path.resolve(strict=False).relative_to(refinement.ROOT).as_posix()
        return atom_counts.get(relative)

    monkeypatch.setattr(refinement, "_xyz_atom_count", fake_xyz_atom_count)
    monkeypatch.setattr(
        refinement,
        "infer_forming_bond_metadata",
        lambda source_row, *, reactant_relative_path, product_relative_path: {
            "available": reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz"),
            "method": "geometry_topology_diff",
            "atom_indices_zero_based": [0, 6] if reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz") else [],
            "atom_indices_one_based": [1, 7] if reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz") else [],
            "atom_symbols": ["C", "N"] if reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz") else [],
            "family_recipe_source": "data/rmg_extensions/families/DHA_Crosslinking/groups.py" if reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz") else None,
            "family_recipe_form_bonds": [{"left_label": "*1", "order": 1, "right_label": "*4"}] if reactant_relative_path.endswith("lysinoalanine_crosslink/reactant.xyz") else [],
            "all_new_bonds": [],
        },
    )


def test_computational_gap_plan_surfaces_targets_ceilings_and_write_back_artifacts(mocked_geometry_inputs: None):
    payload = build_computational_gap_refinement_plan_artifact()

    assert payload["summary"]["target_count"] == 7
    by_id = {row["id"]: row for row in payload["targets"]}

    assert by_id["hexanal_radical_quench"]["promotion_ceiling"] == "ranking_only"
    assert "data/lit/computational_priors.json" in by_id["hexanal_radical_quench"]["write_back_artifacts"]
    assert by_id["hexanal_radical_quench"]["dft"]["spin"] == 1
    assert by_id["hexanal_radical_quench"]["surrogate"]["available"] is False
    assert by_id["lysinoalanine_crosslink"]["xtb"]["status"] == "ready"
    assert by_id["lysinoalanine_crosslink"]["xtb"]["atom_count_match"] is True
    assert by_id["lysinoalanine_crosslink"]["surrogate"]["available"] is True
    assert by_id["lysinoalanine_crosslink"]["surrogate"]["kind"] == "dha_michael_addition"
    assert by_id["lysinoalanine_crosslink"]["surrogate"]["barrier_kcal_mol"] == 16.0
    assert by_id["lysinoalanine_crosslink"]["forming_bond"]["available"] is True
    assert by_id["lysinoalanine_crosslink"]["forming_bond"]["atom_indices_zero_based"] == [0, 6]
    assert by_id["aa_ring_open_dicarbonyl"]["xtb"]["status"] == "ready"
    assert by_id["aa_ring_open_dicarbonyl"]["xtb"]["atom_count_match"] is True
    assert by_id["aa_ring_open_dicarbonyl"]["surrogate"]["available"] is True
    assert by_id["aa_ring_open_dicarbonyl"]["surrogate"]["short_label"] == "Family 14 HCW surrogate"
    assert by_id["quinone_cys_michael"]["xtb"]["status"] == "ready"
    assert by_id["quinone_cys_michael"]["xtb"]["atom_count_match"] is True
    assert by_id["quinone_cys_michael"]["surrogate"]["available"] is False
    assert by_id["pe_schiff_base"]["xtb"]["status"] == "ready"
    assert by_id["pe_schiff_base"]["dft"]["status"] == "ready_for_dft"
    assert by_id["pe_schiff_base"]["promotion_ceiling"] == "bounded_calibration"
    assert by_id["pe_amadori"]["xtb"]["status"] == "ready"
    assert by_id["pe_amadori"]["dft"]["status"] == "ready_for_dft"
    assert by_id["pe_amadori"]["promotion_ceiling"] == "bounded_calibration"
    assert by_id["asparagine_sugar_explicit_water_cluster"]["promotion_ceiling"] == "uncertainty_narrowing_only"
    assert by_id["asparagine_sugar_explicit_water_cluster"]["dft"]["use_explicit_solvent"] is True


def test_computational_gap_manifests_capture_xtb_and_dft_readiness_without_wet_lab_overclaim(mocked_geometry_inputs: None):
    xtb_manifest = build_computational_gap_xtb_job_manifest()
    dft_manifest = build_computational_gap_dft_job_manifest()

    xtb_by_id = {job["target_id"]: job for job in xtb_manifest["jobs"]}
    dft_by_id = {job["target_id"]: job for job in dft_manifest["jobs"]}

    assert xtb_manifest["summary"]["job_count"] == 7
    assert xtb_manifest["summary"]["ready_count"] == 6
    assert xtb_manifest["summary"]["seed_required_count"] == 1
    assert xtb_manifest["summary"]["blocked_atom_count_mismatch_count"] == 0
    assert xtb_by_id["hexanal_radical_quench"]["status"] == "ready"
    assert xtb_by_id["hexanal_radical_quench"]["surrogate_available"] is False
    assert xtb_by_id["lysinoalanine_crosslink"]["status"] == "ready"
    assert xtb_by_id["lysinoalanine_crosslink"]["surrogate_available"] is True
    assert xtb_by_id["lysinoalanine_crosslink"]["surrogate_kind"] == "dha_michael_addition"
    assert xtb_by_id["aa_ring_open_dicarbonyl"]["status"] == "ready"
    assert xtb_by_id["aa_ring_open_dicarbonyl"]["surrogate_short_label"] == "Family 14 HCW surrogate"
    assert xtb_by_id["quinone_cys_michael"]["status"] == "ready"
    assert xtb_by_id["quinone_cys_michael"]["surrogate_available"] is False
    assert xtb_by_id["pe_schiff_base"]["status"] == "ready"
    assert xtb_by_id["pe_amadori"]["status"] == "ready"

    assert dft_manifest["summary"]["job_count"] == 7
    assert dft_manifest["summary"]["seed_required_count"] == 1
    assert dft_manifest["summary"]["blocked_atom_count_mismatch_count"] == 0
    assert dft_manifest["summary"]["ready_for_dft_count"] == 6
    assert dft_by_id["hexanal_radical_quench"]["status"] == "ready_for_dft"
    assert dft_by_id["lysinoalanine_crosslink"]["status"] == "ready_for_dft"
    assert dft_by_id["lysinoalanine_crosslink"]["surrogate_barrier_kcal_mol"] == 16.0
    assert dft_by_id["lysinoalanine_crosslink"]["forming_bond_atoms"] == [0, 6]
    assert dft_by_id["aa_ring_open_dicarbonyl"]["status"] == "ready_for_dft"
    assert dft_by_id["aa_ring_open_dicarbonyl"]["surrogate_kind"] == "ascorbic_upstream_dicarbonyl_source"
    assert dft_by_id["quinone_cys_michael"]["status"] == "ready_for_dft"
    assert dft_by_id["pe_schiff_base"]["status"] == "ready_for_dft"
    assert dft_by_id["pe_amadori"]["status"] == "ready_for_dft"


def test_dft_ingestion_report_handles_missing_results_file():
    manifest_payload = build_computational_gap_dft_job_manifest()
    payload = build_computational_gap_dft_ingestion_artifact(manifest_payload=manifest_payload)
    markdown = render_computational_gap_dft_ingestion_markdown(payload)

    assert payload["summary"]["available_result_count"] == 0
    assert "No DFT results available yet." in markdown


def test_computational_gap_dft_promotion_payload_only_promotes_completed_non_fast_results():
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

    payload = build_computational_gap_dft_promotion_payload(
        priors_payload=priors_payload,
        execution_payload=execution_payload,
    )

    assert payload["summary"]["promoted_count"] == 1
    assert payload["promotions"][0]["reaction_key"] == "hexanal_radical_quench"