"""Tests for :mod:`src.external_validation` (S20.5 Lane A.1)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.external_validation import (
    BENCHMARK_DIR,
    FLAVOR_REFERENCE_PATH,
    build_inventory,
    build_holdout_bundles,
    write_holdout_bundles,
    EXTERNAL_VALIDATION_EVIDENCE_CLASS,
    render_markdown,
    write_artifact,
)
from src.benchmark_validation import assess_matrix_benchmark_evidence, evaluate_benchmark_payload


@pytest.fixture(scope="module")
def inventory():
    return build_inventory()


def test_inventory_has_at_least_one_executable_candidate(inventory):
    """Lane A.1 acceptance: there must be at least one externally
    scoreable candidate today; the rest are honestly classified as
    blocked. If this drops to zero, either the alias map regressed or
    the flavor reference payload changed shape."""

    groups = inventory.by_eligibility()
    assert len(groups.get("executable_candidate", [])) >= 1, (
        "expected at least one executable_candidate; "
        "check the Lane A.2 hold-out bundle specs and compound aliases "
        "against data/lit/flavor_reference_payloads.json"
    )


def test_inventory_marks_already_benchmarked_pbma_furfural_as_redundant(inventory):
    row = next(c for c in inventory.candidates if c.anchor_id == "hernandez_2023_furfural_ratio_anchor")
    assert row.eligibility == "redundant_with_panel"


def test_inventory_promotes_registry_backed_holdouts_to_executable(inventory):
    anchor_ids = {c.anchor_id for c in inventory.by_eligibility().get("executable_candidate", [])}
    assert "bi_2020_raw_pea_hexanal_point" in anchor_ids
    assert "liu_2023_ppi_hexanal_band" in anchor_ids
    assert "li_2026_spi_wg_hme_2_pentylfuran_control_point" in anchor_ids


def test_every_candidate_has_provenance(inventory):
    for c in inventory.candidates:
        assert c.anchor_id, "anchor_id is required"
        assert c.section.endswith("_anchors"), c.section
        assert c.eligibility in {
            "executable_candidate",
            "narrative_only",
            "redundant_with_panel",
        }, c.eligibility
        # Compound + matrix should always be carried through even when
        # numeric extraction fails — they drive the eligibility reason.
        assert c.compound, c.anchor_id
        assert c.matrix_context, c.anchor_id


def test_executable_candidates_have_numeric_value(inventory):
    """An executable_candidate must carry a finite point value (mid of
    band when applicable) — otherwise Lane A.2 cannot synthesise a
    benchmark JSON for it."""

    for c in inventory.by_eligibility().get("executable_candidate", []):
        assert c.point_value is not None, c.anchor_id
        assert c.point_value > 0, c.anchor_id


def test_panel_compounds_loaded(inventory):
    # The current calibration panel should expose at least the canonical
    # MFT / FFT / hexanal / furfural axes; if this list goes empty we are
    # almost certainly looking at an empty data/benchmarks directory.
    panel = set(inventory.panel_compounds)
    assert any("furfural" in c for c in panel)
    assert any("hexanal" in c for c in panel)


def test_render_markdown_starts_with_banner(inventory):
    text = render_markdown(inventory)
    first_line = text.splitlines()[0]
    assert first_line.startswith("<!-- Auto-regenerated"), first_line
    assert "generate_external_validation_inventory.py" in first_line


def test_jsonable_round_trip(inventory):
    payload = inventory.to_jsonable()
    # Make sure the artifact serialises cleanly.
    raw = json.dumps(payload, sort_keys=True)
    restored = json.loads(raw)
    assert restored["summary"]["total_candidates"] == len(inventory.candidates)


def test_write_artifact_creates_paired_files(inventory, tmp_path: Path):
    paths = write_artifact(inventory, output_dir=tmp_path, basename="ext_val_test")
    assert paths["markdown"].exists()
    assert paths["json"].exists()
    md = paths["markdown"].read_text(encoding="utf-8")
    assert md.startswith("<!-- Auto-regenerated")
    js = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert js["summary"]["total_candidates"] == len(inventory.candidates)


def test_source_paths_resolve():
    assert FLAVOR_REFERENCE_PATH.exists(), FLAVOR_REFERENCE_PATH
    assert BENCHMARK_DIR.exists(), BENCHMARK_DIR


def test_holdout_bundles_materialize_as_external_validation_only():
    bundles = build_holdout_bundles()
    assert len(bundles) >= 4
    assert sum(bundle.matched_compound_count() for bundle in bundles) >= 7

    for bundle in bundles:
        assert bundle.intake_payload["evidence_class"] == EXTERNAL_VALIDATION_EVIDENCE_CLASS
        assert bundle.benchmark_payload["evidence_class"] == EXTERNAL_VALIDATION_EVIDENCE_CLASS
        assert bundle.benchmark_payload["source_metadata"]["evidence_class"] == EXTERNAL_VALIDATION_EVIDENCE_CLASS
        assert bundle.benchmark_payload["metadata"]["evidence_class"] == EXTERNAL_VALIDATION_EVIDENCE_CLASS

        evidence = assess_matrix_benchmark_evidence(bundle.benchmark_payload)
        assert evidence.external_data_status == EXTERNAL_VALIDATION_EVIDENCE_CLASS
        assert evidence.promotable is False


def test_holdout_bundles_evaluate_against_runtime_surface():
    bundles = build_holdout_bundles()
    for bundle in bundles:
        evaluation = evaluate_benchmark_payload(bundle.benchmark_payload, benchmark_id=bundle.bundle_id)
        assert evaluation.supported is True, bundle.bundle_id
        assert evaluation.comparisons, bundle.bundle_id
        assert any(row.matched_name is not None for row in evaluation.comparisons), bundle.bundle_id


def test_regenerating_holdout_bundles_over_the_committed_files_is_a_noop(tmp_path: Path):
    """The committed hold-out bundles are frozen evidence (2026-09-01, cleaning Phase 0).

    They carry primary-source corrections (Wave W values/uncertainties, `conditions.buffer`)
    that `_HOLDOUT_BUNDLE_SPECS` does not know about, so a regeneration must leave them
    untouched. Seed scratch copies of the committed JSON + YAML, regenerate on top, and
    require byte-identical files.
    """
    from src.external_validation import (
        EXTERNAL_VALIDATION_BENCHMARK_DIR,
        EXTERNAL_VALIDATION_PROTOCOL_DIR,
    )

    benchmark_dir = tmp_path / "benchmarks"
    protocol_dir = tmp_path / "protocols"
    benchmark_dir.mkdir()
    protocol_dir.mkdir()
    committed = {}
    for src_dir, dst_dir, pattern in (
        (EXTERNAL_VALIDATION_BENCHMARK_DIR, benchmark_dir, "*.json"),
        (EXTERNAL_VALIDATION_PROTOCOL_DIR, protocol_dir, "*.yaml"),
    ):
        for path in sorted(src_dir.glob(pattern)):
            text = path.read_text(encoding="utf-8")
            (dst_dir / path.name).write_text(text, encoding="utf-8")
            committed[dst_dir / path.name] = text
    assert committed, "no committed hold-out bundles found"

    written = write_holdout_bundles(
        build_holdout_bundles(), protocol_dir=protocol_dir, benchmark_dir=benchmark_dir
    )
    assert written["protocols"] == [] and written["benchmarks"] == []
    assert len(written["skipped"]) == len(committed)
    for path, before in committed.items():
        assert path.read_text(encoding="utf-8") == before, f"regeneration changed {path.name}"


def test_overwriting_holdout_bundles_carries_curated_keys_forward(tmp_path: Path):
    """With overwrite=True the spec wins on shared keys, but curated-only keys survive."""
    from src.external_validation import EXTERNAL_VALIDATION_BENCHMARK_DIR

    benchmark_dir = tmp_path / "benchmarks"
    benchmark_dir.mkdir()
    for path in EXTERNAL_VALIDATION_BENCHMARK_DIR.glob("*.json"):
        (benchmark_dir / path.name).write_text(path.read_text(encoding="utf-8"), encoding="utf-8")

    write_holdout_bundles(
        build_holdout_bundles(),
        protocol_dir=tmp_path / "protocols",
        benchmark_dir=benchmark_dir,
        overwrite=True,
    )
    for path in EXTERNAL_VALIDATION_BENCHMARK_DIR.glob("*.json"):
        before = json.loads(path.read_text(encoding="utf-8"))
        after = json.loads((benchmark_dir / path.name).read_text(encoding="utf-8"))
        assert after["conditions"]["buffer"] == before["conditions"]["buffer"], path.name
        for compound, block in before["measured_volatiles"].items():
            for key in block:
                assert key in after["measured_volatiles"][compound], f"{path.name}: {compound}.{key} lost"
