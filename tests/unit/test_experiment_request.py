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
        benchmark_id="cys_ribose_150C_Mottram1994",
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
