import pytest

from src.literature_runtime import (
    build_family_upstream_contract,
    build_flavor_axis_summary,
    describe_retention_runtime,
    query_family_runtime_priors,
    query_flavor_reference_entries,
    query_retention_reference_entries,
)


def test_build_flavor_axis_summary_penalizes_sulfur_positive_but_strecker_poor_profile():
    projection_metadata = {
        "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 12.0},
        "furfural": {"compound": "Furfural", "observable_ppb": 200.0},
        "pyraz": {"compound": "2,5-Dimethylpyrazine", "observable_ppb": 40.0},
    }

    summary = build_flavor_axis_summary(
        projection_metadata=projection_metadata,
        sugars=["fructose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=8.5,
    )

    assert summary["strecker_balance_score"] < 0.4
    assert summary["strecker_gap_penalty"] > 0.0
    assert summary["pyrazine_burden"] > summary["pyrazine_signal_ppb"]
    assert summary["pyrazine_penalty"] > 0.0


def test_build_flavor_axis_summary_penalizes_expected_but_missing_furanones():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["alanine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
    )

    assert summary["furanone_expected"] == ["HEMF", "DMHF"]
    assert summary["furanone_observed"] == []
    assert summary["furanone_support_score"] == 0.0
    assert summary["furanone_penalty"] > 0.0
    assert summary["furanone_missing"] == ["HEMF", "DMHF"]


def test_build_flavor_axis_summary_keeps_native_soy_thiamine_inactive_by_default():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
    )

    assert summary["thiamine_pathway_active"] is False
    assert summary["thiamine_provenance_mode"] == "inactive"
    assert summary["thiamine_availability_source"] == "native_matrix_default_inactive"
    assert summary["thiamine_availability_explicit"] is False


def test_build_flavor_axis_summary_uses_explicit_thiamine_metadata_for_pbma_like_context():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
    )

    assert summary["thiamine_pathway_active"] is True
    assert summary["thiamine_provenance_mode"] == "mixed_thiamine_plus_pentose"
    assert summary["thiamine_availability_source"] == "pbma_fortified"
    assert summary["thiamine_availability_explicit"] is True
    assert summary["thiamine_mft_fraction_estimate"] > 0.0


def test_describe_retention_runtime_for_soy_hexanal_is_temperature_and_time_sensitive():
    cool = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=25.0,
        time_minutes=2.0,
        process_state="ambient_slurry",
    )
    hot_short = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=95.0,
        time_minutes=2.0,
        process_state="heated_matrix",
    )
    hot_long = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=95.0,
        time_minutes=30.0,
        process_state="heated_matrix",
    )

    assert hot_short["dynamic_retention_factor"] > cool["dynamic_retention_factor"]
    assert hot_long["dynamic_retention_factor"] < hot_short["dynamic_retention_factor"]
    assert "Ince et al. (2024)" in " ".join(hot_short["retention_reference_sources"])


def test_describe_retention_runtime_for_pea_markers_carries_structured_ph_release_provenance():
    hexanal = describe_retention_runtime(
        "Hexanal",
        protein_type="pea_iso",
        temperature_celsius=25.0,
        time_minutes=10.0,
        process_state="ambient_slurry",
    )
    pentylfuran = describe_retention_runtime(
        "2-Pentylfuran",
        protein_type="pea_iso",
        temperature_celsius=25.0,
        time_minutes=10.0,
        process_state="ambient_slurry",
    )

    assert hexanal["retention_runtime_mode"] == "direct_binding_plus_ph_release_reference"
    assert "Karolkowski et al. (2021)" in " ".join(hexanal["retention_reference_sources"])
    assert pentylfuran["retention_runtime_mode"] == "ph_release_reference"
    assert "Karolkowski et al. (2021)" in " ".join(pentylfuran["retention_reference_sources"])


def test_describe_retention_runtime_for_extrusion_states_is_aw_sensitive():
    hydrated = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=150.0,
        time_minutes=3.0,
        water_activity=0.60,
        process_state="aqueous_pre_extrusion_model",
    )
    dry = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=165.0,
        time_minutes=3.0,
        water_activity=0.35,
        process_state="extrusion_structured",
    )

    assert "extrusion_" in hydrated["retention_runtime_mode"]
    assert "extrusion_" in dry["retention_runtime_mode"]
    assert dry["dynamic_retention_factor"] < hydrated["dynamic_retention_factor"]
    assert dry["extrusion_moisture_factor"] < hydrated["extrusion_moisture_factor"]


def test_build_flavor_axis_summary_surfaces_secondary_strecker_reference_markers():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "three_mb": {"compound": "3-Methylbutanal", "observable_ppb": 11.0},
            "benz": {"compound": "Benzaldehyde", "observable_ppb": 30.0},
        },
        sugars=["glucose"],
        amino_acids=["leucine", "phenylalanine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=6.0,
    )

    assert summary["strecker_signals_ppb"]["3_methylbutanal"] == 11.0
    assert summary["strecker_signals_ppb"]["benzaldehyde"] == 30.0
    assert summary["strecker_reference_targets_ppb"]["3_methylbutanal"] > 0.0
    assert summary["strecker_reference_targets_ppb"]["benzaldehyde"] > 0.0


def test_build_flavor_axis_summary_surfaces_priority_family_lanes():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "hex": {"compound": "Hexanal", "observable_ppb": 65.0},
            "furan": {"compound": "2-Pentylfuran", "observable_ppb": 18.0},
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 7.0},
        },
        sugars=["ribose", "glucose"],
        amino_acids=["cysteine"],
        additives=["yeast extract", "grape seed extract"],
        lipids=["sunflower oil"],
        interventions=["yeast_fermentation", "protease_hydrolysis"],
        protein_type="soy_iso",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
    )

    assert {"02", "07", "08", "10"}.issubset(set(summary["active_family_lanes"]))
    assert summary["family_lane_summary"]["02"]["lipid_marker_signal_ppb"] > 0.0
    assert "Hexanal" in summary["family_lane_summary"]["02"]["benchmark_ready_targets"]
    assert "lincoln_2025_polyphenol_crosstalk_v1" in summary["family_lane_summary"]["02"]["competition_prior_ids"]
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["adverse_marker_generation_and_retention"]["retention_reference_count"] > 0
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["carbonyl_competition_and_crosstalk"]["active"] is True
    assert summary["family_lane_summary"]["02"]["maillard_closure_pressure"] > 0.0
    assert summary["family_lane_summary"]["07"]["dominant_donor_class"] == "pentose"
    assert summary["family_lane_summary"]["08"]["suppression_pressure_active"] is True
    assert summary["family_lane_summary"]["10"]["precursor_release_active"] is True
    assert summary["family_lane_summary"]["10"]["off_note_cleanup_active"] is True
    assert summary["family_lane_adjustments"]["maillard_closure_delta"] < 0.0


def test_build_flavor_axis_summary_surfaces_all_family_lanes_and_adjustments():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "hex": {"compound": "Hexanal", "observable_ppb": 65.0},
            "furan": {"compound": "2-Pentylfuran", "observable_ppb": 18.0},
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 7.0},
            "hmf": {"compound": "5-Hydroxymethylfurfural (HMF)", "observable_ppb": 22.0},
            "hemf": {"compound": "HEMF", "observable_ppb": 3.5},
        },
        sugars=["ribose", "glucose"],
        amino_acids=["cysteine"],
        additives=["yeast extract", "grape seed extract", "glutathione", "IMP", "GMP"],
        lipids=["sunflower oil"],
        interventions=["yeast_fermentation", "protease_hydrolysis"],
        protein_type="myco",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
    )

    assert summary["active_family_lanes"] == ["01", "02", "07", "10", "08", "03", "04", "05", "09", "06"]
    assert summary["family_lane_summary"]["03"]["thiamine_support_score"] > 0.0
    assert summary["family_lane_summary"]["04"]["nucleotide_support_active"] is True
    assert summary["family_lane_summary"]["05"]["glutathione_active"] is True
    assert summary["family_lane_summary"]["06"]["matrix_scope_active"] is True
    assert summary["family_lane_summary"]["09"]["severity_signal_ppb"] > 0.0
    assert "per_lane" in summary["family_lane_adjustments"]
    assert "thiamine_fragmentation_support" in summary["family_prior_bundle"]
    assert any(marker["marker_id"] == "thiamineavailability" for marker in summary["family_state_markers"])
    assert any(marker["marker_id"] == "caramelizationseverity" for marker in summary["family_state_markers"])
    assert summary["family_upstream_contract"]["dominant_donor_class"] == "pentose"
    assert summary["family_upstream_contract"]["effective_pH"] == pytest.approx(5.15)
    assert "thiamine" in summary["family_upstream_contract"]["added_precursors"]


def test_build_family_upstream_contract_reweights_donors_and_adds_bounded_thiamine():
    contract = build_family_upstream_contract(
        sugars=["ribose", "glucose"],
        amino_acids=["cysteine"],
        additives=[],
        interventions=["yeast_fermentation", "protease_hydrolysis"],
        protein_type="soy_iso",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
        molar_ratios={"ribose": 1.0, "glucose": 1.0, "cysteine": 0.5},
    )

    assert contract["dominant_donor_class"] == "pentose"
    assert contract["donor_pool_factors"]["ribose"] > contract["donor_pool_factors"]["glucose"]
    assert contract["effective_molar_ratios"]["ribose"] > contract["effective_molar_ratios"]["glucose"]
    assert contract["effective_pH"] == pytest.approx(5.15)
    assert contract["pretreatment_active"] is True
    assert contract["thiamine_mode"] == "mixed_thiamine_plus_pentose"
    assert contract["added_precursors"] == ["thiamine"]
    assert contract["added_precursor_ratios"]["thiamine"] > 0.0


def test_query_family_runtime_priors_returns_family_aware_strecker_prior():
    rows = query_family_runtime_priors(
        runtime_lane="strecker_crosstalk",
        family="02",
        compound_name="catechin",
    )

    assert rows
    assert rows[0]["runtime_lane"] == "strecker_crosstalk"
    assert rows[0]["prior_section"] == "strecker_crosstalk_priors"
    assert rows[0]["family"]["slr_family"] == "02"


def test_query_flavor_reference_entries_separates_scoring_targets_from_reference_only_entries():
    scoring_rows = query_flavor_reference_entries(family="01", scoring_only=True)
    methional_rows = query_flavor_reference_entries(entry_id="hernandez_2023_methional_panel")

    assert any(row["id"] == "hernandez_2023_methional_panel" for row in scoring_rows)
    assert methional_rows[0]["pipeline_role"] in {"primary_target", "secondary_marker", "diagnostic_marker", "optimization_constraint"}


def test_query_retention_reference_entries_returns_family_and_matrix_filtered_rows():
    soy_hexanal_rows = query_retention_reference_entries(
        family="02",
        compound_name="Hexanal",
        protein_type="soy_iso",
    )

    assert soy_hexanal_rows
    assert any(row["id"] == "xu_2023_spi_hexanal_temporal_profile" for row in soy_hexanal_rows)
    assert all(row["matrix_family"].startswith("soy") for row in soy_hexanal_rows)