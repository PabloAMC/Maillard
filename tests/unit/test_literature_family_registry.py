import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_family_registry import (  # noqa: E402
    build_family_payload_coverage_artifact,
    get_family_prior_entries,
    iter_flavor_reference_entries,
    iter_matrix_decision_panel_entries,
)


def test_family_payload_coverage_tracks_all_slr_families():
    payload = build_family_payload_coverage_artifact()

    assert payload["summary"]["family_count"] == 10
    by_slr = {row["slr_family"]: row for row in payload["families"]}
    assert by_slr["01"]["total_primary_payload_count"] > 0
    assert by_slr["02"]["total_primary_payload_count"] > 0
    assert by_slr["03"]["total_primary_payload_count"] > 0
    assert by_slr["04"]["total_primary_payload_count"] > 0
    assert by_slr["10"]["total_supporting_payload_count"] > 0


def test_family_registry_returns_family_specific_prior_and_reference_entries():
    thiamine_priors = get_family_prior_entries(family="03")
    sulfur_refs = list(iter_flavor_reference_entries(family="03"))
    nucleotide_panel = list(iter_matrix_decision_panel_entries(family="04"))

    assert any(row["id"] == "cerny_2007_thiamine_split_v1" for row in thiamine_priors)
    assert any(row["id"] == "hofmann_1997_beef_mft_band" for row in sulfur_refs)
    assert any(row["canonical_name"] == "imp" for row in nucleotide_panel)