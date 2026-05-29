import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.reporting import _build_scientific_surface


def _load(relative_path: str) -> dict:
    with open(ROOT / relative_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def test_track0_and_track1_artifacts_exist_and_are_exposed_in_reporting():
    surface = _build_scientific_surface(ROOT)

    assert surface["slr_incorporation_matrix"] == "data/lit/slr_incorporation_matrix.json"
    assert surface["chemistry_family_scope_registry"] == "data/lit/chemistry_family_scope_registry.json"
    assert surface["family_ingestion_plan_registry"] == "data/lit/family_ingestion_plan.json"
    assert surface["family_identifier_contract"] == "results/validation/family_identifier_contract.md"
    assert surface["family_strategy_policy"] == "results/validation/family_strategy_policy.md"
    assert surface["flavor_reference_payloads"] == "data/lit/flavor_reference_payloads.json"
    assert surface["matrix_family_coverage_registry"] == "data/lit/matrix_family_coverage_registry.json"
    assert surface["retention_reference_payloads"] == "data/lit/retention_reference_payloads.json"
    assert surface["process_gap_registry"] == "data/lit/process_gap_registry.json"
    assert surface["chemistry_family_scope"] == "results/validation/chemistry_family_scope.md"
    assert surface["family_ingestion_plan"] == "results/validation/family_ingestion_plan.md"
    assert surface["family_identifier_contract_json"] == "results/validation/family_identifier_contract.json"
    assert surface["family_strategy_policy_json"] == "results/validation/family_strategy_policy.json"
    assert surface["family_payload_coverage"] == "results/validation/family_payload_coverage.md"
    assert surface["matrix_family_coverage"] == "results/validation/matrix_family_coverage.md"
    assert surface["reaction_benchmark_set"] == "data/lit/reaction_benchmark_set.json"
    assert surface["mlp_candidate_registry"] == "data/lit/mlp_candidate_registry.json"
    assert surface["mlp_external_benchmark_evidence"] == "data/lit/mlp_external_benchmark_evidence.json"
    assert surface["geometry_benchmark_set"] == "data/lit/geometry_benchmark_set.json"
    assert surface["literature_backlog"] == "results/validation/literature_backlog.md"
    assert surface["literature_backlog_json"] == "results/validation/literature_backlog.json"
    assert surface["deep_research_runtime_queue"] == "results/validation/deep_research_runtime_queue.md"
    assert surface["deep_research_runtime_queue_json"] == "results/validation/deep_research_runtime_queue.json"


def test_family_ingestion_plan_registry_prioritizes_first_wave_extension_lanes():
    payload = _load("data/lit/family_ingestion_plan.json")
    by_slr = {entry["slr_family"]: entry for entry in payload["families"]}

    assert by_slr["02"]["strategic_posture"] == "immediate_expansion_lane"
    assert by_slr["07"]["runtime_concept"] == "carbonyl_donor_hierarchy"
    assert by_slr["10"]["runtime_concept"] == "fermentation_pretreatment_node"
    assert by_slr["08"]["preferred_payload_types"][2] == "safety_payload"
    assert by_slr["06"]["strategic_posture"] == "matrix_scope_lane"


def test_family_metadata_is_present_on_literature_payload_surfaces():
    flavor = _load("data/lit/flavor_reference_payloads.json")
    retention = _load("data/lit/retention_reference_payloads.json")
    priors = _load("data/lit/computational_priors.json")
    intake = _load("data/lit/benchmark_intake_registry.json")
    gaps = _load("data/lit/process_gap_registry.json")
    panel = _load("data/lit/matrix_decision_panel.json")

    assert flavor["section_family_metadata"]["sulfur_reference_anchors"]["chemistry_family"] == "amino_acid_sugar_core"
    assert retention["section_family_metadata"]["aldehydes"]["chemistry_family"] == "lipid_oxidation_and_carbonylic_crosstalk"
    assert priors["section_family_metadata"]["thiamine_pathway_priors"]["slr_family_source"] == "03"
    assert intake["eligible_references"][0]["payload_role"] == "benchmark_intake"
    assert gaps["entries"][0]["payload_role"] == "structural_gap_entry"
    assert panel["target_class_family_metadata"]["umami_support_markers"]["chemistry_family"] == "nucleotide_and_ribose_support"


def test_slr_incorporation_matrix_covers_new_track1_sources():
    payload = _load("data/lit/slr_incorporation_matrix.json")
    by_id = {entry["paper_id"]: entry for entry in payload["entries"]}

    required_ids = {
        "foods_2022_spi_free_sh",
        "karolkowski_2021_ppi_ph_release",
        "xu_2023_spi_temporal_release",
        "jafc_3c05991_hexanal_binding",
        "ince_2024_hexanal_glycinin",
        "zhang_2026_spi_vsc_retention",
        "foods_2023_pbma_acrylamide_ages",
        "choi_2024_garden_pea_acrylamide",
        "knol_2009_acrylamide_kinetics",
        "trikusuma_2019_pea_uht_aroma_panel",
        "lincoln_2025_polyphenol_crosstalk",
        "troise_2018_soy_amadori",
        "laemont_barringer_2023_pyrazine_ph",
        "wang_2021_peptide_pyrazines",
        "blank_fay_1996_furanones",
    }

    assert required_ids.issubset(by_id)
    assert by_id["foods_2022_spi_free_sh"]["incorporation_status"] == "encoded_modeled_not_shown"
    assert by_id["karolkowski_2021_ppi_ph_release"]["incorporation_status"] == "encoded_modeled_shown"
    assert by_id["xu_2023_spi_temporal_release"]["incorporation_status"] == "encoded_modeled_shown"
    assert by_id["jafc_3c05991_hexanal_binding"]["incorporation_status"] == "encoded_modeled_shown"
    assert by_id["trikusuma_2019_pea_uht_aroma_panel"]["incorporation_status"] == "encoded_modeled_shown"
    assert by_id["lincoln_2025_polyphenol_crosstalk"]["incorporation_status"] == "encoded_modeled_shown"
    assert by_id["choi_2024_garden_pea_acrylamide"]["incorporation_status"] == "encoded_not_yet_shown"


def test_benchmark_intake_registry_encodes_trikusuma_and_lincoln_artifacts():
    payload = _load("data/lit/benchmark_intake_registry.json")
    by_id = {entry["id"]: entry for entry in payload["eligible_references"]}

    assert by_id["cerny_guntz_dubini_2008"]["runtime_artifacts"][0]["artifact_id"] == "thiamine_cys_xylose_145C_Cerny2008"
    assert by_id["trikusuma_2019"]["status"] == "encoded"
    assert by_id["trikusuma_2019"]["key_values"]["tracked_uht_markers_ug_per_l"]["hexanal"] == 782.0
    assert by_id["trikusuma_2019"]["runtime_artifacts"][0]["artifact_id"] == "pea_isolate_uht_140C_Trikusuma2019"
    for entry_id, artifact_id in {
        "pmc_2026_hme_hexanal_baseline": "li_2026_spi_wg_hme_hexanal_control_point",
        "acs_2022_pba_lysine_loss_benchmark": "acs_2022_pba_lysine_loss",
        "pmc11049305_spirulina_offnote_anchor": "pmc11049305_spirulina_beta_ionone_oav_floor",
        "pmc12155365_sunflower_roasted_anchor": "pmc12155365_sunflower_4_vinylguaiacol_fd_point",
        "pmc_2024_pba_cml_cel_ranges_anchor": "pmc_2024_pba_cml_cel_ranges",
        "pmc_12648097_acrylamide_mitigation_anchor": "pmc_12648097_acrylamide_mitigation",
        "pmid_1904866_pentosidine_equivalence_anchor": "pmid_1904866_aa_pentosidine_equivalence_v1",
        "pmc5992167_amadori_pe_burden_anchor": "pmc5992167_amadori_pe_food_matrix_burden",
        "wang_2012_gsh_xylose_sulfur_uplift": "wang_xu_glutathione_peptide_support_v1",
        "ohsu_2025_kokumi_casr_anchor": "ohsu_2025_kokumi_casr_support_v1",
        "soladoye_2020_sous_vide_euc_anchor": "soladoye_2020_low_temp_euc_window_v1",
        "ahlberg_2021_yeast_extract_grade_anchor": "ahlberg_2021_yeast_extract_nucleotide_grade_window_v1",
        "cui_2022_mushroom_nucleotide_anchor": "cui_2022_mushroom_gmp_euc_window_v1",
        "voelker_2021_thiamine_kinetics": "voelker_2021_thiamine_arrhenius_v1",
        "arabshahi_1988_aw_thiamine_kinetics": "arabshahi_1988_aw_dependent_thiamine_ea_v1",
        "huang_2022_thiamine_metal_catalysis": "huang_2022_thiamine_metal_catalysis_v1",
        "blank_grosch_1991_hdmf_anchor": "blank_grosch_1991_beef_hdmf_band",
        "liu_2023_ppi_offnote_baseline": "liu_2023_ppi_ibmp_band",
        "marquez_ruiz_2014_oleic_oav_anchor": "marquez_ruiz_2014_oleic_nonanal_oav_band",
        "messina_2022_pbma_oil_oav_anchor": "messina_2022_pbma_oil_oav_panel",
        "ref41_ppi_sulfur_binding": "ref41_ppi_sulfur_volatile_binding_v1",
        "acs_jafc_3c08432_crosstalk_cleanup_link": "rizzello_2024_lactic_fermentation_cleanup",
        "maillard_van_boekel_1992_sugar_reactivity_hierarchy": "maillard_van_boekel_1992_sugar_reactivity_hierarchy_v1",
        "blank_1997_rhamnose_proline_hdmf_anchor": "blank_1997_rhamnose_proline_hdmf_uplift_v1",
        "brands_2002_mgo_hdmf_anchor": "brands_2002_mgo_hdmf_c3_route_v1",
        "wang_2022_lab_hexanal_cleanup_anchor": "wang_2022_lab_hexanal_cleanup_oav_target",
        "bhandari_1998_beta_cd_aldehyde_binding_anchor": "bhandari_1998_beta_cd_aldehyde_binding_v1",
        "zhang_2022_unsaturated_aldehyde_potency_anchor": "zhang_2022_unsaturated_aldehyde_offnote_potency_v1",
        "frontiers_2022_hcw_aa_arrhenius_anchor": "frontiers_2022_hcw_aa_arrhenius_v1",
        "scielo_brasil_aa_crosslink_hierarchy_anchor": "scielo_brasil_aa_crosslink_hierarchy_v1",
        "pmc9351765_crosspy_trapping_anchor": "pmc9351765_crosspy_mft_scavenging_v1",
        "uspto_ptacts_2023_yeast_extract_anchor": "uspto_ptacts_2023_yeast_extract_mft_oav_band",
        "wageningen_ref9_hme_rework_hydration_anchor": "wageningen_ref9_hme_rework_hydration_collapse",
        "acs_foodscitech_2024_hme_firmness_anchor": "acs_foodscitech_2024_hme_firmness_window",
        "jafc_2019_ref21_pea_gum_arabic_architecture_anchor": "jafc_2019_ref21_pea_gum_arabic_architecture_state",
    }.items():
        assert any(item["artifact_id"] == artifact_id for item in by_id[entry_id]["runtime_artifacts"])
    assert any(item["artifact_id"] == "mottram_2001_lipid_aldehyde_mft_quench_v1" for item in by_id["mottram_2001_mft_quench_buffering_anchor"]["runtime_artifacts"])
    assert any(item["artifact_id"] == "mottram_2001_carnosine_buffered_mft_uplift" for item in by_id["mottram_2001_mft_quench_buffering_anchor"]["runtime_artifacts"])
    assert by_id["yeo_mottram_2023_lecithin_crosstalk_anchor"]["runtime_artifacts"][0]["artifact_id"] == "yeo_mottram_2023_soy_lecithin_thiophene_uplift"
    assert by_id["comunian_2021_thiamine_encapsulation"]["runtime_artifacts"][0]["artifact_id"] == "comunian_2021_thiamine_encapsulation"
    assert by_id["lincoln_2025"]["status"] == "ready_for_directional_prior_encoding"
    assert by_id["lincoln_2025"]["runtime_artifacts"][0]["artifact_id"] == "lincoln_2025_polyphenol_crosstalk_v1"
    assert by_id["blank_devaud_grosch_2003_g6p_hdmf_prior"]["citation_aliases"] == ["JAFC DOI:10.1021/jf034037p (2003)"]
    assert "polyphenol" in json.dumps(by_id["lincoln_2025"]).lower()


def test_hme_family11_intake_registry_captures_pdf_method_but_keeps_closure_blockers_explicit():
    payload = _load("data/lit/benchmark_intake_registry.json")
    by_id = {entry["id"]: entry for entry in payload["eligible_references"]}

    hme = by_id["pmc_2026_hme_hexanal_baseline"]

    assert hme["doi"] == "10.3390/foods15050912"
    assert hme["key_values"]["extrusion_moisture_wt_pct"] == 57.0
    assert hme["key_values"]["barrel_temp_profile_C"][-1] == 160.0
    assert hme["key_values"]["wg_neutral_protease_u_per_g"] == 250.0
    assert hme["key_values"]["flavour_optimum_pretreatment_time_min"] == 40.0
    assert hme["key_values"]["2-pentylfuran_control_ug_per_kg"] == 221.51
    assert hme["key_values"]["1-hexanol_control_ug_per_kg"] == 20.04
    assert any("final-blend pH" in item for item in hme["what_it_does_not_support"])
    assert any("water-activity" in item for item in hme["what_it_does_not_support"])


def test_operational_benchmark_intake_registry_excludes_markdown_backlog_candidates():
    payload = _load("data/lit/benchmark_intake_registry.json")

    for entry in payload["eligible_references"]:
        assert "extracted_from_markdown" not in (entry.get("observable_panel_tags") or [])
        assert entry.get("status") != "pending_json_payload"
        assert entry.get("matrix_family") != "unknown_see_details"


def test_computational_priors_include_lincoln_crosstalk_prior():
    payload = _load("data/lit/computational_priors.json")
    prior = next(entry for entry in payload["strecker_crosstalk_priors"] if entry["id"] == "lincoln_2025_polyphenol_crosstalk_v1")
    thiamine_prior_ids = {entry["id"] for entry in payload["thiamine_pathway_priors"]}
    donor_prior_ids = {entry["id"] for entry in payload["carbonyl_donor_priors"]}
    lipid_offnote_prior_ids = {entry["id"] for entry in payload["lipid_offnote_priors"]}
    ascorbic_prior_ids = {entry["id"] for entry in payload["ascorbic_pathway_priors"]}
    melanoidin_prior_ids = {entry["id"] for entry in payload["melanoidin_trapping_priors"]}
    nucleotide_prior_ids = {entry["id"] for entry in payload["nucleotide_pathway_priors"]}

    assert prior["effect_direction"] == "suppress_strecker_and_moderate_oxidative_crosstalk"
    assert "glucose" in prior["required_sugars"]
    assert "catechin" in prior["polyphenol_examples"]
    assert {"voelker_2021_thiamine_arrhenius_v1", "arabshahi_1988_aw_dependent_thiamine_ea_v1", "huang_2022_thiamine_metal_catalysis_v1"}.issubset(thiamine_prior_ids)
    assert {"soladoye_2020_low_temp_euc_window_v1", "ahlberg_2021_yeast_extract_nucleotide_grade_window_v1", "cui_2022_mushroom_gmp_euc_window_v1"}.issubset(nucleotide_prior_ids)
    assert "maillard_van_boekel_1992_sugar_reactivity_hierarchy_v1" in donor_prior_ids
    assert "blank_1997_rhamnose_proline_hdmf_uplift_v1" in donor_prior_ids
    assert any(entry["id"] == "brands_2002_mgo_hdmf_c3_route_v1" for entry in payload["furanone_priors"])
    assert {"mottram_2001_lipid_aldehyde_mft_quench_v1", "zhang_2022_unsaturated_aldehyde_offnote_potency_v1"}.issubset(lipid_offnote_prior_ids)
    assert {"frontiers_2022_hcw_aa_arrhenius_v1", "scielo_brasil_aa_crosslink_hierarchy_v1", "pmid_1904866_aa_pentosidine_equivalence_v1"}.issubset(ascorbic_prior_ids)
    assert {"pmc9351765_crosspy_mft_scavenging_v1", "jafc_2019_ref21_pea_gum_arabic_architecture_v1"}.issubset(melanoidin_prior_ids)
    assert any(entry["id"] == "bhandari_1998_beta_cd_aldehyde_binding_v1" for entry in payload["retention_binding_priors"])


def test_computational_priors_promote_mycoprotein_to_a_first_class_matrix_family():
    payload = _load("data/lit/computational_priors.json")

    accessibility = next(entry for entry in payload["accessibility_windows"] if entry["protein_type"] == "myco")
    correction = next(entry for entry in payload["matrix_corrections"] if entry["protein_type"] == "myco")

    assert accessibility["provenance_tier"] == "literature_bounded_provisional"
    assert accessibility["uncertainty_posture"] == "directional_only"
    assert "extrusion_structured" in accessibility["process_state_applicability"]
    assert correction["provenance_tier"] == "literature_bounded_provisional"
    assert correction["uncertainty_posture"] == "directional_only"


def test_flavor_reference_payloads_cover_sulfur_strecker_pyrazine_and_furanones():
    payload = _load("data/lit/flavor_reference_payloads.json")

    assert any(entry["compound"] == "2-methyl-3-furanthiol" for entry in payload["sulfur_reference_anchors"])
    assert any(entry["id"] == "uspto_ptacts_2023_yeast_extract_mft_oav_band" for entry in payload["sulfur_reference_anchors"])
    assert {"2-methylbutanal", "3-methylbutanal", "benzaldehyde"}.issubset(
        {entry["compound"] for entry in payload["strecker_reference_anchors"]}
    )
    assert any(entry["compound"] == "methylpyrazine" for entry in payload["pyrazine_reference_anchors"])
    assert any(entry["compound"] == "furfural" for entry in payload["carbonyl_reference_anchors"])
    assert {"HEMF", "HDMF"}.issubset({entry["compound"] for entry in payload["furanone_reference_anchors"]})
    assert any(entry["compound"] == "3-isobutyl-2-methoxypyrazine" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "liu_2023_ppi_hexanal_band" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "wang_2022_lab_hexanal_cleanup_oav_target" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "marquez_ruiz_2014_oleic_nonanal_oav_band" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "messina_2022_pbma_oil_oav_panel" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "pmc11049305_spirulina_beta_ionone_oav_floor" for entry in payload["off_note_reference_anchors"])
    assert any(entry["id"] == "pmc12155365_sunflower_2_methylbutanal_fd_point" for entry in payload["strecker_reference_anchors"])
    assert any(entry["id"] == "pmc12155365_sunflower_4_vinylguaiacol_fd_point" for entry in payload["off_note_reference_anchors"])
    assert all(entry.get("pipeline_role") for section in payload.values() if isinstance(section, list) for entry in section)


def test_retention_reference_payloads_capture_reversibility_and_censoring():
    payload = _load("data/lit/retention_reference_payloads.json")

    pea_aldehydes = payload["aldehydes"]["pea_protein"]
    karolkowski_hexanal = next(entry for entry in pea_aldehydes if entry["id"] == "karolkowski_2021_ppi_hexanal_ph_release")
    assert karolkowski_hexanal["runtime_surrogate"]["type"] == "ph_release_modifier"
    assert karolkowski_hexanal["provenance_tier"] == "literature_derived_transfer"

    soy_aldehydes = payload["aldehydes"]["soy_protein"]
    ince_entry = next(entry for entry in soy_aldehydes if entry["id"] == "ince_2024_glycinin_hexanal_binding")
    assert ince_entry["reversibility_assumption"] == "explicitly_reversible_non_covalent"
    xu_hexanal = next(entry for entry in soy_aldehydes if entry["id"] == "xu_2023_spi_hexanal_temporal_profile")
    assert xu_hexanal["runtime_surrogate"]["type"] == "temporal_attenuation"

    pea_furan = payload["furans"]["pea_protein"][0]
    assert pea_furan["id"] == "karolkowski_2021_ppi_2_pentylfuran_native_panel"

    soy_furans = payload["furans"]["soy_protein"]
    xu_furan = next(entry for entry in soy_furans if entry["id"] == "xu_2023_spi_2_pentylfuran_temporal_profile")
    assert xu_furan["runtime_surrogate"]["floor"] == 0.55

    heated_furan = next(entry for entry in soy_furans if entry["id"] == "shu_2024_heated_spi_2_pentylfuran_censored")
    assert heated_furan["numeric_reference"]["value"] == "not_detected"


def test_process_gap_registry_separates_structural_gaps_from_codable_payloads():
    payload = _load("data/lit/process_gap_registry.json")
    gap_ids = {entry["gap_id"] for entry in payload["entries"]}

    assert "intact_spi_ppi_quantified_mft_fft" in gap_ids
    assert "pipi_spi_pyrazine_absolute_quantitation" in gap_ids
    assert "aqueous_pipi_spi_acrylamide_kinetics" in gap_ids


def test_safety_reference_payloads_distinguish_endpoints_from_kinetic_provenance():
    payload = _load("data/lit/safety_reference_payloads.json")
    by_id = {entry["id"]: entry for entry in payload["entries"]}

    assert by_id["foods_2023_pbma_acrylamide_ages"]["kind"] == "finished_product_reference"
    assert by_id["choi_2024_pea_acrylamide_asparagine"]["kind"] == "precursor_correlation_reference"
    assert by_id["knol_2009_acrylamide_kinetics"]["kind"] == "kinetic_model_reference"
    assert by_id["squeo_2023_pbpi_acrylamide"]["report_visibility"] == "default"
    assert by_id["knol_2009_acrylamide_kinetics"]["report_visibility"] == "extended"