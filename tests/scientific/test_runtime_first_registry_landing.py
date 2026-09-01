import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_prior_registry import query_family_prior_entries
from src.safety import get_safety_reference_payload


def _load_json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8"))


def test_runtime_first_batch_is_landed_in_operational_registries():
    process_payload = _load_json("data/lit/process_state_calibrations.json")
    computational_payload = _load_json("data/lit/computational_priors.json")
    process_ids = {entry["id"] for entry in process_payload["entries"]}
    assert {
        "acs_jafc_3c05991_ppi_spi_partitioning",
        "malia_2025_pea_free_sh_crosscheck",
        "raman_sds_extrusion_disulfide_severity",
        "rizzello_2024_lactic_fermentation_cleanup",
        "zhao_2022_moromi_precursor_release",
        "ordoudi_2014_hmf_peak_window",
        "hrncirik_2014_coconut_oil_thermal_profile",
        "comunian_2021_thiamine_encapsulation",
        "mottram_2001_carnosine_buffered_mft_uplift",
        "yeo_mottram_2023_soy_lecithin_thiophene_uplift",
        "wageningen_ref9_hme_rework_hydration_collapse",
        "acs_foodscitech_2024_hme_firmness_window",
        "jafc_2019_ref21_pea_gum_arabic_architecture_state",
        "pmc5992167_amadori_pe_food_matrix_burden",
    }.issubset(process_ids)

    family11_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="lipid_maillard_crosstalk")}
    assert "acs_jafc_3c02618_mft_disulfide_trapping_v1" in family11_priors
    assert "acs_jafc_0c01925_protein_binding_hierarchy_v1" in family11_priors
    assert "ref41_ppi_sulfur_volatile_binding_v1" in family11_priors

    donor_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="carbonyl_donor_hierarchy")}
    assert "blank_1997_rhamnose_proline_hdmf_uplift_v1" in donor_priors

    support_priors = {row["id"] for row in query_family_prior_entries(supporting_family="08")}
    assert "nakagawa_2004_isoflavone_dicarbonyl_sink_v1" in support_priors
    assert "blank_2001_epoxydecenal_guardrail_v1" in support_priors
    assert "hidalgo_zamora_2004_4hne_pentylpyrrole_v1" in support_priors
    assert "mottram_2001_lipid_aldehyde_mft_quench_v1" in support_priors
    assert "zhang_2022_unsaturated_aldehyde_offnote_potency_v1" in support_priors
    assert "bhandari_1998_beta_cd_aldehyde_binding_v1" in support_priors

    thiamine_prior_ids = {entry["id"] for entry in computational_payload["thiamine_pathway_priors"]}
    assert {"cerny_2007_thiamine_split_v1", "voelker_2021_thiamine_arrhenius_v1", "arabshahi_1988_aw_dependent_thiamine_ea_v1", "huang_2022_thiamine_metal_catalysis_v1"}.issubset(thiamine_prior_ids)

    nucleotide_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="nucleotide_and_ribose_support")}
    assert "liardon_1991_r5p_donor_potency_v1" in nucleotide_priors
    assert "aliani_2005_donor_potency_nucleotide_context_v1" in nucleotide_priors
    assert "blank_devaud_grosch_2003_g6p_hdmf_uplift_v1" in nucleotide_priors
    assert "soladoye_2020_low_temp_euc_window_v1" in nucleotide_priors
    assert "ahlberg_2021_yeast_extract_nucleotide_grade_window_v1" in nucleotide_priors
    assert "cui_2022_mushroom_gmp_euc_window_v1" in nucleotide_priors

    sulfur_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="glutathione_and_peptide_support")}
    assert "huang_2021_sulfur_oav_support_v1" in sulfur_priors
    assert "ohsu_2025_kokumi_casr_support_v1" in sulfur_priors

    capping_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="polyphenol_amino_capping")}
    assert "jafc_2019_ref24_polyphenol_thiol_capping_v1" in capping_priors

    caramelization_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="carbohydrate_pyrolysis_and_caramelization")}
    assert "glomb_1995_3dg_fragmentation_stoichiometry_v1" in caramelization_priors

    ascorbic_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="ascorbic_acid_maillard")}
    assert {"frontiers_2022_hcw_aa_arrhenius_v1", "scielo_brasil_aa_crosslink_hierarchy_v1", "pmid_1904866_aa_pentosidine_equivalence_v1"}.issubset(ascorbic_priors)

    melanoidin_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="melanoidin_polymerization")}
    assert "pmc9351765_crosspy_mft_scavenging_v1" in melanoidin_priors
    assert "jafc_2019_ref21_pea_gum_arabic_architecture_v1" in melanoidin_priors

    furanone_prior_ids = {entry["id"] for entry in computational_payload["furanone_priors"]}
    assert "mottram_nobrega_2002_furanone_sulfur_bridge_v1" in furanone_prior_ids
    assert "brands_2002_mgo_hdmf_c3_route_v1" in furanone_prior_ids

    dft_prior_ids = {entry["reaction_key"] for entry in computational_payload["dft_kinetic_priors"]["entries"]}
    assert {"pyrraline_from_3dg", "furosine_from_3dg", "pe_schiff_base", "pe_amadori"}.issubset(dft_prior_ids)

    backlog_payload = _load_json("results/literature/deep_research_backlog.json")
    backlog_rows = {row["citation"]: row for row in backlog_payload["items"]}
    for citation in {
        "Rizzello et al. (2024)",
        "Zhao et al. (2022)",
        "J. Agric. Food Chem. 2019 (Ref. 24)",
        "Liardon, de Weck-Gaudard & Philippossian (1991)",
        "Huang et al. (2021)",
        "JAFC DOI:10.1021/jf034037p (2003)",
        "Blank et al. (2001)",
        "Ordoudi et al. (2014 / PMC12484514)",
        "Hrncirik & Zeelenberg (2014)",
        "Aliani & Farmer (2005)",
        "Blank, Devaud & Grosch (2003)",
        "Glomb & Monnier (1995)",
        "Hidalgo & Zamora (2004)",
        "ACS APTS (Ref. 24)",
        "Mottram & Nobrega (2002 / Chapter 9 review)",
        "PMC PMCID:PMC4419266 (Ref. 18)",
        "Comunian et al. (2021)",
        "Voelker, Taylor & Mauer (2021)",
        "Arabshahi & Lund (1988)",
        "Huang (2022)",
        "Wang, Z. et al. (2012)",
        "Ohsu et al. (2025)",
        "Soladoye et al. (2020)",
        "Ahlberg & Mohammadi (2021)",
        "Cui et al. (2022)",
        "Blank & Grosch (1991)",
        "Liu, Y. (2023 thesis)",
        "Marquez-Ruiz et al. (2014)",
        "DOI ref. 41 in raw/11_maillard_lipid_crosstalk.md",
        "ACS JAFC 3c08432",
        "van Boekel (2001)",
        "Mottram et al. (2001)",
        "Wang et al. (2022)",
        "Yeo & Mottram (2023)",
        "Zhang et al. (2022)",
        "Frontiers Nutr. 2022 (Ref. 22)",
        "SciELO Brasil (Ref. 4)",
        "PMC PMCID:PMC9351765 (Ref. 12)",
        "USPTO PTACTS (2023 compiled)",
        "PMC 2024 (PMC12451096)",
        "PMC PMCID:PMC12648097 (Ref. 5)",
        "Wageningen Ref. 9",
        "ACS Food Sci. Technol. 2024 (Ref. 1)",
        "J. Agric. Food Chem. 2019 (Ref. 21)",
        "Blank et al. (1997)",
        "Bhandari et al. (1998)",
        "Brands & van Boekel (2002)",
        "PMC11049305 (2024)",
        "PMC12155365 (2025)",
        "PubMed PMID:1904866 (Ref. 5)",
        "PMC PMCID:PMC5992167 (Refs. 16/17)",
    }:
        assert backlog_rows[citation]["status"] == "RUNTIME_BOUND"

    mitigation_reference = get_safety_reference_payload("pmc_12648097_acrylamide_mitigation")
    saponin_reference = get_safety_reference_payload("kocadagli_2016_saponin_acrylamide_modifier")
    assert mitigation_reference is not None
    assert saponin_reference is not None
    assert mitigation_reference["report_visibility"] == "extended"
    assert saponin_reference["kind"] == "contextual_process_modifier"