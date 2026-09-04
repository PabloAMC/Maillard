"""The chemistry-family axis is a closed vocabulary (16 ids in the scope registry)."""
from __future__ import annotations

from src import data_access, data_paths


def _canonical_ids() -> set[str]:
    scope = data_access.load_json(data_paths.CHEMISTRY_FAMILY_SCOPE_REGISTRY)
    return {str(f["family_id"]) for f in scope["families"]}


def test_scope_registry_and_ingestion_plan_agree_on_the_sixteen_families():
    plan = data_access.load_json(data_paths.FAMILY_INGESTION_PLAN)
    assert {str(f["family_id"]) for f in plan["families"]} == _canonical_ids()
    assert len(_canonical_ids()) == 16


def test_every_intake_reference_uses_a_canonical_chemistry_family():
    """2026-09-02: five aliases (11 records) were rewritten; an alias map in
    literature_family_registry used to hide them."""
    intake = data_access.load_json(data_paths.BENCHMARK_INTAKE_REGISTRY)
    canonical = _canonical_ids()
    bad = sorted(
        {r.get("chemistry_family") for r in intake["eligible_references"]} - canonical - {None}
    )
    assert not bad, f"non-canonical chemistry_family ids in the intake registry: {bad}"


def test_process_gap_registry_uses_canonical_chemistry_families():
    gaps = data_access.load_json(data_paths.PROCESS_GAP_REGISTRY)
    canonical = _canonical_ids()
    assert all(str(e.get("chemistry_family")) in canonical for e in gaps["entries"])
