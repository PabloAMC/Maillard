import pytest

from src.literature_runtime import (
    build_family_upstream_contract,
    build_flavor_axis_summary,
    describe_retention_runtime,
    query_benchmark_intake_entries,
    query_dft_kinetic_priors,
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
    assert summary["family_lane_summary"]["02"]["dominant_marker"] == "Hexanal"
    assert "Hexanal" in summary["family_lane_summary"]["02"]["benchmark_ready_targets"]
    assert summary["family_lane_summary"]["02"]["primary_benchmark_id"] == "trikusuma_2019"
    assert summary["family_lane_summary"]["02"]["benchmark_marker_targets_ug_per_l"]["hexanal"] == pytest.approx(782.0)
    assert "lincoln_2025_polyphenol_crosstalk_v1" in summary["family_lane_summary"]["02"]["competition_prior_ids"]
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["adverse_marker_generation_and_retention"]["retention_reference_count"] > 0
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["adverse_marker_generation_and_retention"]["benchmark_anchor_ids"] == ["trikusuma_2019"]
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["adverse_marker_generation_and_retention"]["dominant_marker"] == "Hexanal"
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["carbonyl_competition_and_crosstalk"]["active"] is True
    assert summary["family_lane_summary"]["02"]["runtime_sub_lanes"]["carbonyl_competition_and_crosstalk"]["donor_pressure"] == pytest.approx(1.0)
    assert summary["family_lane_summary"]["02"]["maillard_closure_pressure"] > 0.0
    assert summary["family_lane_summary"]["11"]["competition_window_active"] is True
    assert summary["family_lane_summary"]["11"]["kinetic_prior_ids"] == ["hexanal_radical_quench"]
    assert summary["family_lane_summary"]["11"]["kinetic_prior"]["active_arrhenius_key"] == "hexanal_radical_quench_xtb_derived"
    assert set(summary["family_lane_summary"]["11"]["benchmark_anchor_ids"]) == {"pmc_2026_hme_hexanal_baseline", "acs_2020_raw_pea_hexanal_baseline"}
    assert summary["family_lane_summary"]["11"]["selected_benchmark_anchor_id"] == "pmc_2026_hme_hexanal_baseline"
    assert summary["family_lane_summary"]["11"]["hexanal_suppression_fraction"] > 0.0
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

    assert summary["active_family_lanes"] == ["01", "02", "07", "10", "08", "11", "13", "03", "04", "05", "09", "06"]
    assert summary["family_lane_summary"]["11"]["kinetic_prior"]["barrier_kj_mol"] == pytest.approx(31.72)
    assert summary["family_lane_summary"]["11"]["hexanal_baseline_anchors_ug_per_kg"]["acs_2020_raw_pea_hexanal_baseline"] == pytest.approx(1260.0)
    assert summary["family_lane_summary"]["03"]["thiamine_support_score"] > 0.0
    assert summary["family_lane_summary"]["04"]["nucleotide_support_active"] is True
    assert summary["family_lane_summary"]["05"]["glutathione_active"] is True
    assert summary["family_lane_summary"]["05"]["benchmark_anchor_ids"] == ["nishimura_abe_2024"]
    assert summary["family_lane_summary"]["05"]["sulfur_peptide_support_score"] > 0.75
    assert summary["family_lane_summary"]["05"]["pyrazine_tradeoff_ratio_vs_free_cysteine"] == pytest.approx(0.75)
    assert summary["family_lane_summary"]["06"]["matrix_scope_active"] is True
    assert summary["family_lane_summary"]["06"]["source_id"] == "mycoprotein"
    assert summary["family_lane_summary"]["06"]["process_state_anchor_ids"] == ["asen_2022", "li_2025"]
    assert summary["family_lane_summary"]["06"]["structural_gap_ids"] == ["ellman_opa_dsc_same_experiment"]
    assert summary["family_lane_summary"]["06"]["matrix_uncertainty_factor"] < 0.75
    assert summary["family_lane_summary"]["13"]["polyphenol_active"] is True
    assert summary["family_lane_summary"]["13"]["kinetic_prior_ids"] == ["quinone_cys_michael"]
    assert summary["family_lane_summary"]["13"]["cysteine_depletion_factor"] > 0.0
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
    assert contract["effective_molar_ratios"]["ribose"] == pytest.approx(1.0)
    assert contract["effective_molar_ratios"]["glucose"] == pytest.approx(1.0)
    assert contract["effective_pH"] == pytest.approx(5.15)
    assert contract["pretreatment_active"] is True
    assert contract["thiamine_mode"] == "mixed_thiamine_plus_pentose"
    assert contract["added_precursors"] == ["thiamine"]
    assert contract["added_precursor_ratios"]["thiamine"] > 0.0


def test_build_family_upstream_contract_calibrates_thiamine_by_ph_and_extrusion_survival():
    optimal = build_family_upstream_contract(
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=["thiamine"],
        protein_type="soy_iso",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
        molar_ratios={"ribose": 1.0, "cysteine": 1.0, "thiamine": 0.1},
    )
    alkaline = build_family_upstream_contract(
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=["thiamine"],
        protein_type="soy_iso",
        pH=8.0,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
        molar_ratios={"ribose": 1.0, "cysteine": 1.0, "thiamine": 0.1},
    )
    extruded = build_family_upstream_contract(
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=["thiamine"],
        protein_type="soy_iso",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
        process_state="extrusion_structured",
        temperature_celsius=180.0,
        time_minutes=2.0,
        water_activity=0.75,
        molar_ratios={"ribose": 1.0, "cysteine": 1.0, "thiamine": 0.1},
    )

    assert optimal["thiamine_fraction_baseline"] == pytest.approx(0.5)
    assert optimal["thiamine_fraction_estimate"] > alkaline["thiamine_fraction_estimate"]
    assert extruded["thiamine_fraction_estimate"] < optimal["thiamine_fraction_estimate"]
    assert optimal["effective_molar_ratios"]["thiamine"] > 0.1
    assert extruded["effective_molar_ratios"]["thiamine"] < optimal["effective_molar_ratios"]["thiamine"]
    assert extruded["thiamine_calibration"]["extrusion_survival_factor"] == pytest.approx(0.04)
    assert "cerny_guntz_dubini_2008" in extruded["thiamine_calibration"]["benchmark_anchor_ids"]
    assert optimal["added_precursor_ratios"] == {}
    assert extruded["added_precursor_ratios"] == {}


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


def test_query_family_runtime_priors_returns_family_07_pyrazine_control_prior():
    rows = query_family_runtime_priors(
        runtime_lane="pyrazine_control",
        family="07",
    )

    assert rows
    assert rows[0]["id"] == "slr_pyrazine_control_surface_v1"
    assert rows[0]["prior_section"] == "pyrazine_control_priors"
    assert rows[0]["family"]["slr_family"] == "01"


def test_query_family_runtime_priors_returns_family_05_sulfur_peptide_prior():
    rows = query_family_runtime_priors(
        runtime_lane="sulfur_peptide_support",
        family="05",
    )

    assert rows
    assert rows[0]["id"] == "wang_xu_glutathione_peptide_support_v1"
    assert rows[0]["prior_section"] == "sulfur_peptide_priors"
    assert rows[0]["gsh_mft_ratio_vs_free_cysteine"] == pytest.approx(2.25)


def test_build_flavor_axis_summary_calibrates_family_05_support_lane_for_soy_hydrolysate_context():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 9.0},
            "fft": {"compound": "2-Furfurylthiol (FFT)", "observable_ppb": 3.0},
        },
        sugars=["xylose"],
        amino_acids=["cysteine"],
        additives=["glutathione", "soy hydrolysate"],
        lipids=[],
        interventions=["protease_hydrolysis"],
        protein_type="soy_iso",
        pH=6.0,
    )

    family_05 = summary["family_lane_summary"]["05"]

    assert family_05["glutathione_active"] is True
    assert family_05["hydrolysate_active"] is True
    assert family_05["peptide_mode"] == "hydrolysate_supported"
    assert family_05["prior_ids"] == ["wang_xu_glutathione_peptide_support_v1"]
    assert family_05["retention_reference_ids"]


def test_build_flavor_axis_summary_calibrates_family_06_matrix_scope_lane_for_mycoprotein():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 9.0},
            "fft": {"compound": "2-Furfurylthiol (FFT)", "observable_ppb": 3.0},
        },
        sugars=["xylose"],
        amino_acids=["cysteine"],
        additives=["glutathione", "soy hydrolysate"],
        lipids=[],
        interventions=["protease_hydrolysis"],
        protein_type="myco",
        pH=6.0,
    )

    family_06 = summary["family_lane_summary"]["06"]

    assert family_06["matrix_scope_active"] is True
    assert family_06["benchmark_anchor_ids"] == [
        "pmc9905368_spi_hvp_xylose_benchmark",
        "pmc9905368_wheat_gluten_hvp_xylose_benchmark",
    ]
    assert family_06["selected_benchmark_anchor_id"] == "pmc9905368_spi_hvp_xylose_benchmark"
    assert family_06["benchmark_transfer_mode"] == "nearest_source_transfer"
    assert family_06["process_state_transfer_confidence"] > 0.0


def test_build_flavor_axis_summary_calibrates_family_07_donor_lane_with_pyrazine_controls():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "pyraz": {"compound": "2,5-Dimethylpyrazine", "observable_ppb": 18.0},
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 6.0},
        },
        sugars=["ribose", "glucose"],
        amino_acids=["cysteine"],
        additives=["yeast extract"],
        lipids=[],
        interventions=["protease_hydrolysis"],
        protein_type="soy_iso",
        pH=8.7,
    )

    family_07 = summary["family_lane_summary"]["07"]

    assert family_07["dominant_donor_class"] == "pentose"
    assert family_07["prior_ids"] == ["slr_pyrazine_control_surface_v1"]
    assert "laemont_2023_pyrazine_ph_direction" in family_07["benchmark_anchor_ids"]
    assert family_07["peptide_intensification_active"] is True
    assert family_07["pyrazine_pressure_score"] > 0.5


def test_build_flavor_axis_summary_calibrates_family_08_guardrail_lane_with_safety_reference():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "hex": {"compound": "Hexanal", "observable_ppb": 65.0},
        },
        sugars=["glucose"],
        amino_acids=["cysteine"],
        additives=["grape seed extract", "calcium carbonate"],
        lipids=["sunflower oil"],
        protein_type="soy_iso",
        pH=6.2,
    )

    family_08 = summary["family_lane_summary"]["08"]

    assert family_08["suppression_pressure_active"] is True
    assert family_08["benchmark_anchor_ids"] == ["squeo_2023"]
    assert family_08["crosstalk_prior_ids"] == ["lincoln_2025_polyphenol_crosstalk_v1"]
    assert family_08["safety_reference_ids"] == ["squeo_2023_pbpi_acrylamide"]
    assert family_08["acrylamide_reference_mean_ug_per_kg"] == pytest.approx(451.0)
    assert family_08["suppression_pressure_score"] > 0.3


def test_build_flavor_axis_summary_calibrates_family_09_with_furanone_and_carbonyl_anchors():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "furf": {"compound": "Furfural", "observable_ppb": 28.0},
            "hmf": {"compound": "5-Hydroxymethylfurfural (HMF)", "observable_ppb": 12.0},
            "hemf": {"compound": "HEMF", "observable_ppb": 2.4},
        },
        sugars=["ribose"],
        amino_acids=["alanine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
    )

    family_09 = summary["family_lane_summary"]["09"]

    assert family_09["prior_ids"] == ["blank_fay_1996_furanone_expectation_v1"]
    assert "hernandez_2023_furfural_ratio_anchor" in family_09["carbonyl_anchor_ids"]
    assert "blank_fay_1996_hemf_mechanistic_anchor" in family_09["furanone_anchor_ids"]
    assert family_09["furanone_support_active"] is True
    assert family_09["furfural_reference_ratio"] > 1.0


def test_build_flavor_axis_summary_calibrates_family_10_with_pretreatment_anchors():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=["yeast extract", "soy hydrolysate"],
        lipids=[],
        interventions=["protease_hydrolysis"],
        protein_type="soy_iso",
        pH=5.8,
    )

    family_10 = summary["family_lane_summary"]["10"]

    assert family_10["benchmark_anchor_ids"] == ["nishimura_abe_2024", "matoba_1988_nucleotide_hydrolysis"]
    assert family_10["selected_benchmark_anchor_id"] == "nishimura_abe_2024"
    assert "wang_xu_glutathione_peptide_support_v1" in family_10["prior_ids"]
    assert "hofmann_1997_beef_mft_band" in family_10["flavor_anchor_ids"]
    assert family_10["free_amino_acid_enrichment_score"] > 0.7


def test_build_flavor_axis_summary_calibrates_family_12_damage_lane_with_safety_anchors():
    summary = build_flavor_axis_summary(
        projection_metadata={},
        sugars=["glucose"],
        amino_acids=["lysine", "asparagine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=6.1,
        process_state="extrusion_structured",
        temperature_celsius=140.0,
        time_minutes=25.0,
        water_activity=0.55,
        molar_ratios={"glucose": 1.0, "lysine": 0.8, "asparagine": 0.2},
    )

    family_12 = summary["family_lane_summary"]["12"]

    assert family_12["damage_guardrail_active"] is True
    assert set(family_12["benchmark_anchor_ids"]) == {
        "acs_2022_pba_lysine_loss_benchmark",
        "acs_ref3_spi_acrylamide_fast_kinetics",
        "foods_2023_cml_cel_proxy_benchmark",
        "ramirez_jimenez_2000_furosine_crossover_benchmark",
    }
    assert family_12["selected_benchmark_anchor_id"] in set(family_12["benchmark_anchor_ids"])
    assert family_12["predicted_acrylamide_ppb"] > 0.0
    assert family_12["predicted_cml_proxy"] > 0.0
    assert family_12["predicted_cel_proxy"] > 0.0
    assert family_12["predicted_furosine_proxy"] > 0.0
    assert family_12["damage_burden_score"] > 0.0
    assert any(marker["family_lane"]["slr_family"] == "12" for marker in summary["family_state_markers"])


def test_build_family_upstream_contract_surfaces_family_13_precursor_sink_on_polyphenol_context():
    contract = build_family_upstream_contract(
        sugars=["ribose"],
        amino_acids=["cysteine", "lysine"],
        additives=["grape seed extract"],
        protein_type="soy_iso",
        pH=6.8,
        temperature_celsius=120.0,
        molar_ratios={"ribose": 1.0, "cysteine": 0.5, "lysine": 0.5},
    )

    family_13 = contract["family_lanes"]["13"]

    assert family_13["active"] is True
    assert family_13["kinetic_prior_ids"] == ["quinone_cys_michael"]
    assert family_13["cysteine_depletion_factor"] > family_13["lysine_depletion_factor"]
    assert contract["effective_molar_ratios"]["cysteine"] < 0.5
    assert contract["summary"]["polyphenol_precursor_sink_active"] is True


def test_build_flavor_axis_summary_surfaces_family_14_bounded_ascorbic_dicarbonyl_source():
    summary = build_flavor_axis_summary(
        projection_metadata={},
        sugars=["glucose"],
        amino_acids=["lysine"],
        additives=["ascorbic acid"],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
        process_state="extrusion_structured",
        temperature_celsius=145.0,
        time_minutes=12.0,
        water_activity=0.64,
    )

    family_14 = summary["family_lane_summary"]["14"]

    assert family_14["active"] is True
    assert family_14["kinetic_prior_ids"] == ["aa_ring_open_dicarbonyl"]
    assert family_14["effective_ea_kj_mol"] < 31.70
    assert family_14["water_activity_modulation_factor"] > 0.9
    assert family_14["dicarbonyl_source_pressure"] > 0.0
    assert family_14["pentosidine_load"] > 0.0
    assert any(marker["family_lane"]["slr_family"] == "14" for marker in summary["family_state_markers"])


def test_build_family_upstream_contract_surfaces_family_15_and_reweights_available_sugar_pool():
    contract = build_family_upstream_contract(
        sugars=["ribose", "glucose"],
        amino_acids=["cysteine"],
        additives=["soy lecithin"],
        lipids=["sunflower oil"],
        protein_type="soy_iso",
        pH=6.2,
        process_state="heated_matrix",
        temperature_celsius=145.0,
        time_minutes=20.0,
        water_activity=0.78,
        molar_ratios={"ribose": 1.0, "glucose": 1.0, "cysteine": 1.0},
    )

    family_15 = contract["family_lanes"]["15"]

    assert family_15["active"] is True
    assert family_15["kinetic_prior_ids"] == ["pe_schiff_base", "pe_amadori"]
    assert family_15["sugar_sink_fraction"] > 0.0
    assert family_15["available_sugar_retention_factor"] < 1.0
    assert contract["effective_molar_ratios"]["ribose"] < 1.0
    assert contract["summary"]["phospholipid_sugar_sink_active"] is True


def test_build_flavor_axis_summary_surfaces_family_16_bounded_melanoidin_trapping_lane():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "fft": {"compound": "2-Furfurylthiol (FFT)", "observable_ppb": 5.0},
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 4.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.7,
        process_state="extrusion_structured",
        temperature_celsius=150.0,
        time_minutes=6.0,
        water_activity=0.62,
    )

    family_16 = summary["family_lane_summary"]["16"]

    assert family_16["active"] is True
    assert family_16["fft_fold_reduction_anchor"] == pytest.approx(16.0)
    assert family_16["melanoidin_mass"] > 0.0
    assert family_16["thiol_scavenging_factor"] > 0.0
    assert any(marker["family_lane"]["slr_family"] == "16" for marker in summary["family_state_markers"])


def test_query_benchmark_intake_and_dft_kinetic_priors_surface_family_02_and_11_contracts():
    family_02_rows = query_benchmark_intake_entries(
        family="02",
        primary_only=True,
        matrix_family="pea_uht_beverage",
    )
    family_03_rows = query_benchmark_intake_entries(
        family="03",
        primary_only=True,
    )
    family_04_rows = query_benchmark_intake_entries(
        family="04",
        primary_only=True,
    )
    family_11_kinetics = query_dft_kinetic_priors(family="11", reaction_key="hexanal_radical_quench")

    assert family_02_rows
    assert family_02_rows[0]["id"] == "trikusuma_2019"
    assert family_02_rows[0]["key_values"]["tracked_uht_markers_ug_per_l"]["hexanal"] == pytest.approx(782.0)
    assert {row["id"] for row in family_03_rows} == {
        "cerny_guntz_dubini_2008",
        "de_leyn_2019",
        "hofmann_schieberle_grosch_1996",
    }
    assert {row["id"] for row in family_04_rows} == {
        "matoba_1988_nucleotide_hydrolysis",
        "mouritsen_2024_umami_thresholds",
        "nakamura_1988_imp_ribose_release",
        "yamaguchi_ninomiya_2000_euc_anchor",
    }
    assert family_11_kinetics
    assert family_11_kinetics[0]["id"] == "hexanal_radical_quench"
    assert family_11_kinetics[0]["active_arrhenius_key"] == "hexanal_radical_quench_xtb_derived"


def test_build_family_upstream_contract_calibrates_family_04_tradeoff_and_adds_bounded_ribose():
    mild = build_family_upstream_contract(
        sugars=[],
        amino_acids=["cysteine"],
        additives=["IMP"],
        protein_type="soy_iso",
        pH=7.0,
        process_state="heated_matrix",
        temperature_celsius=100.0,
        time_minutes=30.0,
    )
    severe = build_family_upstream_contract(
        sugars=[],
        amino_acids=["cysteine"],
        additives=["IMP", "GMP"],
        protein_type="soy_iso",
        pH=7.0,
        process_state="extrusion_structured",
        temperature_celsius=150.0,
        time_minutes=3.0,
        water_activity=0.45,
    )

    mild_family_04 = mild["family_lanes"]["04"]
    severe_family_04 = severe["family_lanes"]["04"]

    assert mild_family_04["nucleotide_survival_factor"] > severe_family_04["nucleotide_survival_factor"]
    assert mild_family_04["umami_support_factor"] > severe_family_04["umami_support_factor"]
    assert severe_family_04["ribose_delivery_factor"] > mild_family_04["ribose_delivery_factor"]
    assert severe_family_04["ribose_shift_active"] is True
    assert severe["nucleotide_calibration"]["benchmark_anchor_ids"] == [
        "matoba_1988_nucleotide_hydrolysis",
        "mouritsen_2024_umami_thresholds",
        "nakamura_1988_imp_ribose_release",
        "yamaguchi_ninomiya_2000_euc_anchor",
    ]
    assert severe["added_precursors"] == ["ribose"]
    assert severe["added_precursor_ratios"]["ribose"] > 0.1


def test_build_flavor_axis_summary_surfaces_family_04_benchmark_context_and_ribose_shift():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 6.0},
        },
        sugars=[],
        amino_acids=["cysteine"],
        additives=["IMP", "GMP"],
        lipids=[],
        protein_type="soy_iso",
        pH=7.0,
        process_state="extrusion_structured",
        temperature_celsius=150.0,
        time_minutes=3.0,
        water_activity=0.45,
    )

    family_04 = summary["family_lane_summary"]["04"]
    assert family_04["nucleotide_support_active"] is True
    assert family_04["nucleotide_survival_factor"] < 0.6
    assert family_04["ribose_delivery_factor"] > 0.3
    assert family_04["ribose_shift_active"] is True
    assert family_04["umami_reference_mode"] == "hydrolyzing_nucleotide_pool"
    assert "Matoba, Terao & Fujimaki (1988), JAFC 36:1033" in family_04["benchmark_anchor_citations"]


def test_build_flavor_axis_summary_surfaces_family_03_benchmark_context_and_extrusion_penalty():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.5,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
        process_state="extrusion_structured",
        temperature_celsius=180.0,
        time_minutes=2.0,
        water_activity=0.75,
    )

    family_03 = summary["family_lane_summary"]["03"]
    assert family_03["extrusion_reference_id"] == "de_leyn_2019_thiamine_retention"
    assert "Cerny & Guntz-Dubini (2008), JAFC 56:5138" in family_03["benchmark_anchor_citations"]
    assert family_03["extrusion_survival_factor"] == pytest.approx(0.04)
    assert family_03["thiamine_reference_yield_mode"] == "mixed_system_optimal_window"
    assert family_03["thiamine_support_score"] < 0.1


def test_query_flavor_reference_entries_separates_scoring_targets_from_reference_only_entries():
    scoring_rows = query_flavor_reference_entries(family="01", scoring_only=True)
    methional_rows = query_flavor_reference_entries(entry_id="hernandez_2023_methional_panel")

    assert any(row["id"] == "hernandez_2023_methional_panel" for row in scoring_rows)
    assert methional_rows[0]["pipeline_role"] in {"primary_target", "secondary_marker", "diagnostic_marker", "optimization_constraint"}


def test_query_flavor_reference_entries_exposes_raw_pea_hexanal_as_reference_only_family_11_anchor():
    rows = query_flavor_reference_entries(entry_id="bi_2020_raw_pea_hexanal_point")

    assert rows
    assert rows[0]["slr_family_source"] == "11"
    assert rows[0]["section"] == "off_note_reference_anchors"
    assert rows[0]["pipeline_role"] == "reference_only"
    assert rows[0]["matrix_context"] == "raw_pea_flour"


def test_query_flavor_reference_entries_exposes_hme_hexanal_control_as_reference_only_family_11_anchor():
    rows = query_flavor_reference_entries(entry_id="li_2026_spi_wg_hme_hexanal_control_point")

    assert rows
    assert rows[0]["slr_family_source"] == "11"
    assert rows[0]["section"] == "off_note_reference_anchors"
    assert rows[0]["pipeline_role"] == "reference_only"
    assert rows[0]["matrix_context"] == "spi_wheat_gluten_hme_control_57pct_moisture"


def test_query_retention_reference_entries_returns_family_and_matrix_filtered_rows():
    soy_hexanal_rows = query_retention_reference_entries(
        family="02",
        compound_name="Hexanal",
        protein_type="soy_iso",
    )

    assert soy_hexanal_rows
    assert any(row["id"] == "xu_2023_spi_hexanal_temporal_profile" for row in soy_hexanal_rows)
    assert all(row["matrix_family"].startswith("soy") for row in soy_hexanal_rows)