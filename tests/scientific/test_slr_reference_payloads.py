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
        "jafc_3c05991_hexanal_binding",
        "ince_2024_hexanal_glycinin",
        "zhang_2026_spi_vsc_retention",
        "foods_2023_pbma_acrylamide_ages",
        "choi_2024_garden_pea_acrylamide",
        "knol_2009_acrylamide_kinetics",
        "troise_2018_soy_amadori",
        "laemont_barringer_2023_pyrazine_ph",
        "wang_2021_peptide_pyrazines",
        "blank_fay_1996_furanones",
    }

    assert required_ids.issubset(by_id)
    assert by_id["foods_2022_spi_free_sh"]["incorporation_status"] == "encoded_modeled_not_shown"
    assert by_id["jafc_3c05991_hexanal_binding"]["incorporation_status"] == "encoded_not_yet_modeled"
    assert by_id["choi_2024_garden_pea_acrylamide"]["incorporation_status"] == "encoded_not_yet_shown"


def test_flavor_reference_payloads_cover_sulfur_strecker_pyrazine_and_furanones():
    payload = _load("data/lit/flavor_reference_payloads.json")

    assert any(entry["compound"] == "2-methyl-3-furanthiol" for entry in payload["sulfur_reference_anchors"])
    assert any(entry["compound"] == "2-methylbutanal" for entry in payload["strecker_reference_anchors"])
    assert any(entry["compound"] == "methylpyrazine" for entry in payload["pyrazine_reference_anchors"])
    assert any(entry["compound"] == "furfural" for entry in payload["carbonyl_reference_anchors"])
    assert any(entry["compound"] == "HEMF" for entry in payload["furanone_reference_anchors"])


def test_retention_reference_payloads_capture_reversibility_and_censoring():
    payload = _load("data/lit/retention_reference_payloads.json")

    soy_aldehydes = payload["aldehydes"]["soy_protein"]
    ince_entry = next(entry for entry in soy_aldehydes if entry["id"] == "ince_2024_glycinin_hexanal_binding")
    assert ince_entry["reversibility_assumption"] == "explicitly_reversible_non_covalent"

    heated_furan = payload["furans"]["soy_protein"][0]
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