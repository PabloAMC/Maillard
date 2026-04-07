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
    process_ids = {entry["id"] for entry in process_payload["entries"]}
    assert "acs_jafc_3c05991_ppi_spi_partitioning" in process_ids
    assert "malia_2025_pea_free_sh_crosscheck" in process_ids
    assert "raman_sds_extrusion_disulfide_severity" in process_ids
    assert "rizzello_2024_lactic_fermentation_cleanup" in process_ids
    assert "zhao_2022_moromi_precursor_release" in process_ids
    assert "ordoudi_2014_hmf_peak_window" in process_ids
    assert "hrncirik_2014_coconut_oil_thermal_profile" in process_ids

    family11_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="lipid_maillard_crosstalk")}
    assert "acs_jafc_3c02618_mft_disulfide_trapping_v1" in family11_priors
    assert "acs_jafc_0c01925_protein_binding_hierarchy_v1" in family11_priors

    support_priors = {row["id"] for row in query_family_prior_entries(supporting_family="08")}
    assert "nakagawa_2004_isoflavone_dicarbonyl_sink_v1" in support_priors
    assert "blank_2001_epoxydecenal_guardrail_v1" in support_priors
    assert "hidalgo_zamora_2004_4hne_pentylpyrrole_v1" in support_priors

    nucleotide_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="nucleotide_and_ribose_support")}
    assert "liardon_1991_r5p_donor_potency_v1" in nucleotide_priors
    assert "aliani_2005_donor_potency_nucleotide_context_v1" in nucleotide_priors
    assert "blank_devaud_grosch_2003_g6p_hdmf_uplift_v1" in nucleotide_priors

    sulfur_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="glutathione_and_peptide_support")}
    assert "huang_2021_sulfur_oav_support_v1" in sulfur_priors

    capping_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="polyphenol_amino_capping")}
    assert "jafc_2019_ref24_polyphenol_thiol_capping_v1" in capping_priors

    caramelization_priors = {row["id"] for row in query_family_prior_entries(chemistry_family="carbohydrate_pyrolysis_and_caramelization")}
    assert "glomb_1995_3dg_fragmentation_stoichiometry_v1" in caramelization_priors

    backlog_payload = _load_json("data/lit/deep_research_backlog.json")
    backlog_rows = {row["citation"]: row for row in backlog_payload["items"]}
    assert backlog_rows["Rizzello et al. (2024)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Zhao et al. (2022)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["J. Agric. Food Chem. 2019 (Ref. 24)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Liardon, de Weck-Gaudard & Philippossian (1991)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Huang et al. (2021)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Blank et al. (2001)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Ordoudi et al. (2014 / PMC12484514)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Hrncirik & Zeelenberg (2014)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Aliani & Farmer (2005)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Blank, Devaud & Grosch (2003)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Glomb & Monnier (1995)"]["status"] == "RUNTIME_BOUND"
    assert backlog_rows["Hidalgo & Zamora (2004)"]["status"] == "RUNTIME_BOUND"

    mitigation_reference = get_safety_reference_payload("pmc_12648097_acrylamide_mitigation")
    saponin_reference = get_safety_reference_payload("kocadagli_2016_saponin_acrylamide_modifier")
    assert mitigation_reference is not None
    assert saponin_reference is not None
    assert mitigation_reference["report_visibility"] == "extended"
    assert saponin_reference["kind"] == "contextual_process_modifier"