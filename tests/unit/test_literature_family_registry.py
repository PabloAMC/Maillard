import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.literature_family_registry as literature_family_registry  # noqa: E402
from src.literature_family_registry import (  # noqa: E402
    build_family_payload_coverage_artifact,
    get_family_prior_entries,
    iter_flavor_reference_entries,
    iter_matrix_decision_panel_entries,
)
from src.family_ingestion_plan import load_family_ingestion_plan  # noqa: E402


def test_family_payload_coverage_tracks_all_slr_families():
    payload = build_family_payload_coverage_artifact()
    plan = load_family_ingestion_plan()

    assert payload["summary"]["family_count"] == len(plan["families"])
    by_slr = {row["slr_family"]: row for row in payload["families"]}
    assert by_slr["01"]["total_primary_payload_count"] > 0
    assert by_slr["02"]["total_primary_payload_count"] > 0
    assert by_slr["03"]["total_primary_payload_count"] > 0
    assert by_slr["04"]["total_primary_payload_count"] > 0
    assert by_slr["07"]["has_runtime_support"] is True
    assert by_slr["10"]["total_supporting_payload_count"] > 0
    assert by_slr["11"]["total_primary_payload_count"] > 0
    assert by_slr["12"]["total_primary_payload_count"] > 0
    assert by_slr["13"]["total_primary_payload_count"] > 0
    assert by_slr["14"]["total_primary_payload_count"] > 0


def test_family_registry_returns_family_specific_prior_and_reference_entries():
    donor_priors = get_family_prior_entries(family="07")
    guardrail_priors = get_family_prior_entries(family="08")
    thiamine_priors = get_family_prior_entries(family="03")
    furanone_priors = get_family_prior_entries(family="09")
    nucleotide_priors = get_family_prior_entries(family="04")
    sulfur_peptide_priors = get_family_prior_entries(family="05")
    sulfur_refs = list(iter_flavor_reference_entries(family="03"))
    nucleotide_panel = list(iter_matrix_decision_panel_entries(family="04"))

    assert any(row["id"] == "maillard_van_boekel_1992_sugar_reactivity_hierarchy_v1" for row in donor_priors)
    assert any(row["id"] == "blank_1997_rhamnose_proline_hdmf_uplift_v1" for row in donor_priors)
    assert any(row["id"] == "bhandari_1998_beta_cd_aldehyde_binding_v1" for row in guardrail_priors)
    assert any(row["id"] == "cerny_2007_thiamine_split_v1" for row in thiamine_priors)
    assert any(row["id"] == "arabshahi_1988_aw_dependent_thiamine_ea_v1" for row in thiamine_priors)
    assert any(row["id"] == "brands_2002_mgo_hdmf_c3_route_v1" for row in furanone_priors)
    assert any(row["id"] == "matoba_1988_nucleotide_hydrolysis_v1" for row in nucleotide_priors)
    assert any(row["id"] == "soladoye_2020_low_temp_euc_window_v1" for row in nucleotide_priors)
    assert any(row["id"] == "ahlberg_2021_yeast_extract_nucleotide_grade_window_v1" for row in nucleotide_priors)
    assert any(row["id"] == "cui_2022_mushroom_gmp_euc_window_v1" for row in nucleotide_priors)
    assert any(row["id"] == "wang_xu_glutathione_peptide_support_v1" for row in sulfur_peptide_priors)
    assert any(row["id"] == "ohsu_2025_kokumi_casr_support_v1" for row in sulfur_peptide_priors)
    assert any(row["id"] == "hofmann_1997_beef_mft_band" for row in sulfur_refs)
    assert any(row["canonical_name"] == "imp" for row in nucleotide_panel)


def test_iter_benchmark_intake_entries_skips_pending_rows_by_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
        registry_path = tmp_path / "benchmark_intake_registry.json"
        registry_path.write_text(
                """
                {
                    "eligible_references": [
                        {
                            "id": "pending_entry",
                            "citation": "Pending Entry",
                            "chemistry_family": "amino_acid_sugar_core",
                            "slr_family_source": "01",
                            "status": "pending_json_payload"
                        },
                        {
                            "id": "ready_entry",
                            "citation": "Ready Entry",
                            "chemistry_family": "amino_acid_sugar_core",
                            "slr_family_source": "01",
                            "status": "ready_for_intake_encoding"
                        }
                    ]
                }
                """,
                encoding="utf-8",
        )

        monkeypatch.setattr(literature_family_registry, "BENCHMARK_INTAKE_REGISTRY_PATH", registry_path)

        assert [row["id"] for row in literature_family_registry.iter_benchmark_intake_entries()] == ["ready_entry"]
        assert [row["id"] for row in literature_family_registry.iter_benchmark_intake_entries(include_pending=True)] == [
                "pending_entry",
                "ready_entry",
        ]