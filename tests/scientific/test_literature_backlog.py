import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_intake_registry import build_literature_backlog_artifact, render_literature_backlog_markdown


def test_literature_backlog_queues_are_exclusive_and_surface_minimum_primary_experiment():
    payload = build_literature_backlog_artifact(ROOT)
    summary = payload["summary"]
    encoded_rows = {row["id"]: row for row in payload["encoded_reference_rows"]}

    assert summary["encoded_reference_count"] >= 5
    assert summary["queue_conflict_count"] == 0
    assert payload["conflicts"]["encoded_and_ready_ids"] == []

    queue_ids = {row["id"] for row in payload["ready_runtime"] + payload["ready_benchmark"]}
    encoded_ids = set(encoded_rows)
    assert queue_ids.isdisjoint(encoded_ids)
    for entry_id in {
        "rizzello_2024_fermentation_cleanup",
        "zhao_2022_moromi_precursor_release",
        "ordoudi_2014_hmf_peak_window",
    }:
        assert encoded_rows[entry_id]["triage_status"] == "ready_calibration"
    for entry_id in {
        "acs_jafc_3c05991_ppi_spi_partitioning",
        "acs_jafc_3c02618_binding_prior",
        "acs_jafc_0c01925_protein_binding_hierarchy",
        "resconi_2023_pbma_beef_identity_benchmark",
        "acs_apts_ref24_3dg_arrhenius_anchor",
        "mottram_nobrega_2002_furanone_bridge",
        "pmc_4419266_pe_interfacial_maillard_kinetics",
        "comunian_2021_thiamine_encapsulation",
        "voelker_2021_thiamine_kinetics",
        "huang_2022_thiamine_metal_catalysis",
        "wang_2012_gsh_xylose_sulfur_uplift",
        "blank_grosch_1991_hdmf_anchor",
        "liu_2023_ppi_offnote_baseline",
        "pmc11049305_spirulina_offnote_anchor",
        "pmc12155365_sunflower_roasted_anchor",
        "pmc_2024_pba_cml_cel_ranges_anchor",
        "pmc_12648097_acrylamide_mitigation_anchor",
        "frontiers_2022_hcw_aa_arrhenius_anchor",
        "pmid_1904866_pentosidine_equivalence_anchor",
        "scielo_brasil_aa_crosslink_hierarchy_anchor",
        "pmc5992167_amadori_pe_burden_anchor",
        "pmc9351765_crosspy_trapping_anchor",
        "uspto_ptacts_2023_yeast_extract_anchor",
        "wageningen_ref9_hme_rework_hydration_anchor",
        "acs_foodscitech_2024_hme_firmness_anchor",
    }:
        assert encoded_rows[entry_id]["encoding_status"] == "encoded_runtime_artifact"

    wet_lab_gaps = {row["gap_id"] for row in payload["wet_lab_blocked"]}
    assert "ppi_meaty_positive_matrix_benchmark" in wet_lab_gaps
    assert "spi_meaty_positive_matrix_benchmark" in wet_lab_gaps
    assert "meaty_off_flavour_safety_tradeoff_panel" in wet_lab_gaps

    minimum_primary_experiment = payload["minimum_primary_experiment"]
    assert minimum_primary_experiment["matrices"] == ["PPI 5% buffered slurry", "SPI 5% buffered slurry"]
    assert "HS-SPME-GC-MS" in minimum_primary_experiment["instrumentation"]

    markdown = render_literature_backlog_markdown(payload)
    assert "Literature Backlog" in markdown
    assert "Ready Runtime" in markdown
    assert "Ready Benchmark" in markdown
    assert "no benchmark-ready backlog remains in the curated intake registry" in markdown
    assert "Wet-Lab Blocked" in markdown
    assert "Minimum Primary Experiment" in markdown