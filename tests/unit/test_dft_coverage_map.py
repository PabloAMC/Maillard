import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dft_coverage_map import build_dft_coverage_map_artifact, render_dft_coverage_map_markdown  # noqa: E402


def test_dft_coverage_map_tracks_family_governance_by_visibility_and_role():
    payload = build_dft_coverage_map_artifact()

    assert payload["summary"]["family_count"] == 10
    by_family = {row["family_id"]: row for row in payload["families"]}

    assert by_family["amino_acid_sugar_core"]["literature_status"] == "literature_arrhenius_available"
    assert by_family["amino_acid_sugar_core"]["benchmark_visibility"] == "benchmark_visible_reaction_family"
    assert "beta_elimination" in by_family["amino_acid_sugar_core"]["benchmark_visible_lanes"]
    assert "thiamine_degradation" in by_family["thiamine_fragmentation_support"]["benchmark_visible_lanes"]
    assert "glutathione_cleavage" in by_family["glutathione_and_peptide_support"]["benchmark_visible_lanes"]
    assert "geom_preopt" in by_family["amino_acid_sugar_core"]["mlp_roles"]
    assert by_family["lipid_oxidation_and_carbonylic_crosstalk"]["benchmark_visibility"] == "benchmark_visible_reaction_family"
    assert by_family["lipid_oxidation_and_carbonylic_crosstalk"]["xtb_status"] == "screen_missing_or_unanchored_barriers_first"
    assert by_family["lipid_oxidation_and_carbonylic_crosstalk"]["next_compute_action"] == (
        "close_literature_and_observable_gaps_then_xTB_screen_lipid_crosstalk_steps_before_selective_DFT"
    )
    assert by_family["alternative_protein_matrix_scope"]["dft_status"] == "not_primary_use_case_for_family_wide_dft"
    assert by_family["fermentation_pretreatment"]["next_compute_action"] == (
        "defer_compute_escalation_until_the_family_is_benchmark_visible"
    )


def test_dft_coverage_map_markdown_surfaces_policy_and_roles():
    markdown = render_dft_coverage_map_markdown(build_dft_coverage_map_artifact())

    assert "DFT Coverage Map" in markdown
    assert "MLP Roles" in markdown
    assert "Lipid oxidation and carbonylic crosstalk" in markdown
    assert "Thiamine degradation and sulfur support" in markdown
    assert "Glutathione and low-molecular-weight peptide support" in markdown
    assert "close_literature_and_observable_gaps_then_xTB_screen_lipid_crosstalk_steps_before_selective_DFT" in markdown
    assert "stay_literature_and_calibration_first_no_family_wide_qm" in markdown