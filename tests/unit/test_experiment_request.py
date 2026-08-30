"""Unit tests for src.experiment_request."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import pytest
import yaml

from src.experiment_request import (
    _build_intake_payload,
    _slug,
    build_requests,
    render_index_markdown,
    write_index,
)
from src.experiment_value import ExperimentCandidate


def _candidate(**overrides: Any) -> ExperimentCandidate:
    base: Dict[str, Any] = dict(
        rank=1,
        benchmark_id="cys_ribose_140C_Hofmann1998",
        compound="2-methyl-3-furanthiol",
        measured_ppb=10.0,
        predicted_p5=0.001,
        predicted_p50=0.5,
        predicted_p95=100.0,
        inside_ci=False,
        envelope_miss_log10=0.0,
        ci_width_log10=4.0,
        odour_threshold_ug_per_kg=1e-4,
        decision_relevance=5.0,
        voi_score=10.0,
        suggested_doe_template="blocking_benchmark_gap",
        rationale="critical meaty odorant",
    )
    base.update(overrides)
    return ExperimentCandidate(**base)


def test_slug_normalises():
    assert _slug("Cys+Ribose 150°C / Mottram 1994") == "cys_ribose_150_c_mottram_1994"
    assert _slug("") == "request"


def test_build_intake_payload_carries_voi_metadata():
    payload = _build_intake_payload(
        request_id="requested_test_123",
        candidate=_candidate(),
        bench={"protein_type": "soy_iso", "conditions": {"temp_C": 100.0}},
        goal="meaty aroma",
        budget_label="3 lab days",
    )
    assert payload["status"] == "pending_lab"
    assert payload["protein_type"] == "soy_iso"
    assert payload["conditions"] == {"temp_C": 100.0}
    assert payload["measured_volatiles"]["2-methyl-3-furanthiol"]["conc_ppb"] is None
    meta = payload["request_metadata"]
    assert meta["originating_voi_rank"] == 1
    assert meta["goal"] == "meaty aroma"
    assert meta["budget_label"] == "3 lab days"


def test_build_requests_writes_files_with_protein_filter(tmp_path: Path):
    candidates = [
        _candidate(rank=1, benchmark_id="bench_soy", compound="2-methyl-3-furanthiol"),
        _candidate(rank=2, benchmark_id="bench_pea", compound="2-furfurylthiol"),
    ]
    bench_index = {
        "bench_soy": {"protein_type": "soy_iso", "conditions": {"temp_C": 110.0}},
        "bench_pea": {"protein_type": "pea_iso", "conditions": {"temp_C": 110.0}},
    }
    protocols_dir = tmp_path / "protocols"
    requests_dir = tmp_path / "requests"
    requests = build_requests(
        candidates,
        top_n=5,
        protein_type="soy",
        protocols_dir=protocols_dir,
        requests_dir=requests_dir,
        benchmark_index=bench_index,
    )
    assert len(requests) == 1
    assert requests[0].benchmark_id == "bench_soy"
    intake = yaml.safe_load(requests[0].intake_yaml_path.read_text(encoding="utf-8"))
    assert intake["protein_type"] == "soy_iso"
    assert requests[0].protocol_md_path.exists()
    assert "VoI rank" in requests[0].protocol_md_path.read_text(encoding="utf-8")

    index_path = write_index(requests, requests_dir=requests_dir)
    assert index_path.exists()
    assert "soy" in index_path.read_text(encoding="utf-8").lower()


def test_render_index_markdown_empty_message():
    md = render_index_markdown([])
    assert "No requests" in md


def test_intake_payload_includes_analytical_context_defaults():
    payload = _build_intake_payload(
        request_id="requested_test_ac",
        candidate=_candidate(compound="2-methyl-3-furanthiol"),
        bench={"protein_type": "soy_iso", "conditions": {"temp_C": 100.0}},
        goal=None,
        budget_label=None,
    )
    ac = payload["analytical_context"]
    assert ac["headspace_method"] == "HS-SPME-GC-MS"
    assert ac["quantification_mode"] == "internal_standard_calibrated"
    assert ac["replicates"] == 3
    assert ac["non_detect_policy"] == "report_lod_and_do_not_backfill"
    # Compound-specific IS hint resolves; lipid anchor included for
    # blocking_benchmark_gap multi-band template.
    assert any("13C-2-methyl-3-furanthiol" in s for s in ac["internal_standards"])
    assert "hexanal-d12" in ac["internal_standards"]


def test_intake_payload_internal_standards_fallback_for_unknown_compound():
    payload = _build_intake_payload(
        request_id="requested_test_unknown",
        candidate=_candidate(
            compound="some_uncatalogued_volatile",
            suggested_doe_template="missing_kinetic_dataset",
        ),
        bench=None,
        goal=None,
        budget_label=None,
    )
    ac = payload["analytical_context"]
    assert ac["headspace_method"] == "LC-MS_MS_timecourse_quench"
    # Falls back to a slug-keyed placeholder rather than silently omitting.
    assert any(
        s.startswith("compound_specific_internal_standard_for_")
        for s in ac["internal_standards"]
    )


def test_protocol_markdown_emits_cro_checklist_and_analytical_block(tmp_path: Path):
    candidates = [
        _candidate(rank=1, benchmark_id="bench_soy", compound="2-furfurylthiol"),
    ]
    bench_index = {
        "bench_soy": {"protein_type": "soy_iso", "conditions": {"temp_C": 110.0}},
    }
    requests = build_requests(
        candidates,
        top_n=1,
        protocols_dir=tmp_path / "protocols",
        requests_dir=tmp_path / "requests",
        benchmark_index=bench_index,
    )
    md = requests[0].protocol_md_path.read_text(encoding="utf-8")
    # Both new sections are present and ordered before "## Data return".
    assert "## Analytical context" in md
    assert "## CRO send-to-lab checklist" in md
    assert md.index("## Analytical context") < md.index("## CRO send-to-lab checklist")
    assert md.index("## CRO send-to-lab checklist") < md.index("## Data return")
    # Checklist mirrors intake schema field names (so a returned YAML lands
    # cleanly via scripts/ingest_results.py).
    assert "internal_standard_calibrated" in md
    assert "report_lod_and_do_not_backfill" in md
    assert "measured_volatiles.2-furfurylthiol" in md
    # Compound-specific IS hint shows up in the analytical block.
    assert "13C-2-furfurylthiol" in md
