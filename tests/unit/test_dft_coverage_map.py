import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dft_coverage_map import build_dft_coverage_map_artifact, build_c4_c5_dft_status, render_dft_coverage_map_markdown  # noqa: E402


def test_dft_coverage_map_tracks_family_governance_by_visibility_and_role():
    payload = build_dft_coverage_map_artifact()

    assert payload["summary"]["family_count"] >= 10
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
    assert "Glutathione and peptide support" in markdown
    assert "close_literature_and_observable_gaps_then_xTB_screen_lipid_crosstalk_steps_before_selective_DFT" in markdown
    assert "stay_literature_and_calibration_first_no_family_wide_qm" in markdown


# ── C4/C5 reaction-level DFT tier tests ───────────────────────────────────────

C4_C5_EXPECTED_REACTIONS = {
    "hexanal_radical_quench",
    "quinone_cys_michael",
    "aa_ring_open_dicarbonyl",
    "pe_schiff_base",
    "pe_amadori",
    "lysinoalanine_crosslink",
}


def test_c4_c5_dft_status_returns_all_six_reactions():
    """build_c4_c5_dft_status() must return exactly the 6 C4/C5 targets."""
    rows = build_c4_c5_dft_status()
    keys = {r["reaction_key"] for r in rows}
    assert C4_C5_EXPECTED_REACTIONS == keys, (
        f"Missing reactions: {C4_C5_EXPECTED_REACTIONS - keys}; "
        f"unexpected: {keys - C4_C5_EXPECTED_REACTIONS}"
    )


def test_c4_c5_dft_status_sentinel_logic():
    """Before real DFT runs, dft_filled must be False (sentinel 99.0 in FAST_BARRIERS[*_dft])."""
    rows = build_c4_c5_dft_status()
    for row in rows:
        assert not row["dft_filled"], (
            f"Reaction {row['reaction_key']} reports dft_filled=True before DFT completes"
        )
        assert row["dft_barrier_kcal"] is None
        # active_barrier must be a real number from xTB or literature anchors
        assert row["active_barrier_kcal"] is not None
        assert isinstance(row["active_barrier_kcal"], float)
        assert 0 < row["active_barrier_kcal"] < 50


def test_c4_c5_dft_status_tiers_are_consistent():
    """Evidence tiers must be recognised values and target_tier must be selective_dft_anchor."""
    valid_tiers = {
        "xtb_derived_gfn2",
        "literature_derived_transfer",
        "mlp_screen_mace",
        "selective_dft_anchor",
    }
    rows = build_c4_c5_dft_status()
    for row in rows:
        assert row["current_tier"] in valid_tiers, (
            f"{row['reaction_key']}: unexpected current_tier={row['current_tier']!r}"
        )
        assert row["target_tier"] == "selective_dft_anchor", (
            f"{row['reaction_key']}: target_tier should be selective_dft_anchor"
        )
        assert row["uncertainty_kj"] > 0


def test_c4_c5_summary_fields_in_artifact():
    """Top-level artifact must expose c4_c5_reaction_count and fill counters."""
    payload = build_dft_coverage_map_artifact()
    summary = payload["summary"]
    assert "c4_c5_reaction_count" in summary
    assert summary["c4_c5_reaction_count"] == 6
    assert "c4_c5_dft_filled_count" in summary
    assert "c4_c5_dft_pending_count" in summary
    # Before any DFT jobs complete, pending == 6
    assert summary["c4_c5_dft_pending_count"] == 6
    assert "c4_c5_reactions" in payload
    assert len(payload["c4_c5_reactions"]) == 6