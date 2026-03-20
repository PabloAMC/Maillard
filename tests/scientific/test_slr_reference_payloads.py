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
    assert surface["flavor_reference_payloads"] == "data/lit/flavor_reference_payloads.json"
    assert surface["retention_reference_payloads"] == "data/lit/retention_reference_payloads.json"
    assert surface["process_gap_registry"] == "data/lit/process_gap_registry.json"


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

    assert by_id["trikusuma_2019"]["status"] == "ready_for_intake_encoding"
    assert by_id["trikusuma_2019"]["key_values"]["tracked_uht_markers_ug_per_l"]["hexanal"] == 782.0
    assert by_id["trikusuma_2019"]["runtime_artifacts"][0]["artifact_id"] == "pea_isolate_uht_140C_Trikusuma2019"
    assert by_id["lincoln_2025"]["status"] == "ready_for_directional_prior_encoding"
    assert by_id["lincoln_2025"]["runtime_artifacts"][0]["artifact_id"] == "lincoln_2025_polyphenol_crosstalk_v1"
    assert "polyphenol" in json.dumps(by_id["lincoln_2025"]).lower()


def test_computational_priors_include_lincoln_crosstalk_prior():
    payload = _load("data/lit/computational_priors.json")
    prior = next(entry for entry in payload["strecker_crosstalk_priors"] if entry["id"] == "lincoln_2025_polyphenol_crosstalk_v1")

    assert prior["effect_direction"] == "suppress_strecker_and_moderate_oxidative_crosstalk"
    assert "glucose" in prior["required_sugars"]
    assert "catechin" in prior["polyphenol_examples"]


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
    assert any(entry["compound"] == "2-methylbutanal" for entry in payload["strecker_reference_anchors"])
    assert any(entry["compound"] == "3-methylbutanal" for entry in payload["strecker_reference_anchors"])
    assert any(entry["compound"] == "benzaldehyde" for entry in payload["strecker_reference_anchors"])
    assert any(entry["compound"] == "methylpyrazine" for entry in payload["pyrazine_reference_anchors"])
    assert any(entry["compound"] == "furfural" for entry in payload["carbonyl_reference_anchors"])
    assert any(entry["compound"] == "HEMF" for entry in payload["furanone_reference_anchors"])
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