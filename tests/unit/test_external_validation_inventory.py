"""Tests for :mod:`src.external_validation` (S20.5 Lane A.1)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.external_validation import (
    BENCHMARK_DIR,
    FLAVOR_REFERENCE_PATH,
    build_inventory,
    render_markdown,
    write_artifact,
)


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
        "check src.external_validation._EXECUTABLE_MATRIX_TAGS and "
        "_COMPOUND_ALIASES against data/lit/flavor_reference_payloads.json"
    )


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
