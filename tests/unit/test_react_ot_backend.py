"""Tests for the React-OT backend plugin and adapter scaffolding.

These tests must pass *without* the secondary `react_ot` conda env or any
React-OT install, so they exercise only the safety-net behaviour:
graceful failure of probe_backend, NotImplementedError from the
single-seed contract, registry wiring, and adapter construction.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mlp_adoption_contract import (  # noqa: E402
    DEFAULT_MLP_CANDIDATE_REGISTRY,
    MLPModelCandidate,
    load_mlp_candidates,
)
from src.mlp_backend_adapters import (  # noqa: E402
    BackendAvailability,
    ReactOTBackendAdapter,
    build_candidate_adapter,
)
from src import react_ot_backend  # noqa: E402


REACT_OT_CANDIDATE_ID = "react_ot_ts_v1"


def _candidate(**overrides) -> MLPModelCandidate:
    base = {
        "candidate_id": REACT_OT_CANDIDATE_ID,
        "model_family": "react_ot",
        "model_name": "deepprinciple_react_ot_sb_pretrained",
        "chemistry_domain": "general_organic",
        "materials_first": False,
        "proposed_role": "ts_initialization",
        "target_motif_families": ["carbonyl_amine_condensation"],
        "benchmark_visible_gap": "test gap",
        "expected_speedup": 5.0,
        "likely_failure_modes": ["benchmark_gap"],
        "fallback_comparator": "sella_plus_pyscf",
        "backend_locator": "python:src.react_ot_backend",
    }
    base.update(overrides)
    return MLPModelCandidate(**base)


def test_registry_includes_react_ot_candidate():
    candidates = load_mlp_candidates()
    matching = [c for c in candidates if c.candidate_id == REACT_OT_CANDIDATE_ID]
    assert len(matching) == 1, "registry must declare react_ot_ts_v1 exactly once"
    candidate = matching[0]
    assert candidate.model_family == "react_ot"
    assert candidate.proposed_role == "ts_initialization"
    assert candidate.backend_locator == "python:src.react_ot_backend"


def test_registry_payload_has_status_marker():
    payload = json.loads(Path(DEFAULT_MLP_CANDIDATE_REGISTRY).read_text())
    rows = [row for row in payload["candidates"] if row["candidate_id"] == REACT_OT_CANDIDATE_ID]
    assert rows, "registry json missing react_ot row"
    assert rows[0]["status"] == "secondary_env_required"


def test_build_candidate_adapter_returns_react_ot_adapter():
    adapter = build_candidate_adapter(_candidate())
    assert isinstance(adapter, ReactOTBackendAdapter)
    assert adapter.model_family == "react_ot"


def test_probe_availability_reports_missing_backend_locator():
    adapter = build_candidate_adapter(_candidate(backend_locator=None))
    probe = adapter.probe_availability()
    assert isinstance(probe, BackendAvailability)
    assert probe.available is False
    assert probe.backend_available is False
    assert "python:src.react_ot_backend" in probe.reason


def test_probe_availability_reports_secondary_env_missing(monkeypatch):
    """When the secondary env is not installed, probe must fail gracefully
    with a reason that points the operator at the setup script."""
    def _fake_probe(model_name, locator):
        return False, "react-ot not importable in secondary env 'react_ot'; run scripts/setup_react_ot_env.sh"

    monkeypatch.setattr(react_ot_backend, "probe_backend", _fake_probe)
    adapter = build_candidate_adapter(_candidate())
    probe = adapter.probe_availability()
    assert probe.available is False
    assert probe.backend_available is False
    assert "scripts/setup_react_ot_env.sh" in probe.reason


def test_probe_availability_reports_contract_gap_when_backend_present(monkeypatch):
    """When the secondary env *is* installed, the adapter still reports
    available=False because the TS-seed-benchmark contract feeds a single
    seed, not a reactant+product pair. backend_available should be True
    so the validator records that React-OT is present-but-skipped."""
    def _fake_probe(model_name, locator):
        return True, "react-ot importable in secondary env 'react_ot'"

    monkeypatch.setattr(react_ot_backend, "probe_backend", _fake_probe)
    adapter = build_candidate_adapter(_candidate())
    probe = adapter.probe_availability()
    assert probe.available is False
    assert probe.backend_available is True
    assert "scripts/recover_ts_react_ot_seed.py" in probe.reason


def test_prepare_ts_seed_raises_with_actionable_message():
    adapter = build_candidate_adapter(_candidate())
    with pytest.raises(NotImplementedError) as exc_info:
        adapter.prepare_ts_seed("3\nHCN\nH 0 0 0\nC 0 0 1\nN 0 0 2\n")
    assert "scripts/recover_ts_react_ot_seed.py" in str(exc_info.value)


def test_predict_ts_returns_missing_checkpoint_without_env_var(monkeypatch):
    monkeypatch.delenv(react_ot_backend.DEFAULT_CHECKPOINT_ENV, raising=False)
    result = react_ot_backend.predict_ts_from_reactant_product(
        reactant_xyz="dummy", product_xyz="dummy", checkpoint=None
    )
    assert result["status"] == "missing_checkpoint"


def test_predict_ts_returns_checkpoint_not_found(tmp_path, monkeypatch):
    fake_checkpoint = tmp_path / "missing.ckpt"
    monkeypatch.delenv(react_ot_backend.DEFAULT_CHECKPOINT_ENV, raising=False)
    result = react_ot_backend.predict_ts_from_reactant_product(
        reactant_xyz="dummy",
        product_xyz="dummy",
        checkpoint=str(fake_checkpoint),
    )
    assert result["status"] == "checkpoint_not_found"
    assert result["checkpoint"] == str(fake_checkpoint)


def test_provenance_file_present_and_well_formed():
    provenance_path = ROOT / "models" / "external" / "react_ot" / "provenance.json"
    assert provenance_path.exists(), "models/external/react_ot/provenance.json must exist"
    payload = json.loads(provenance_path.read_text())
    assert payload["candidate_id"] == REACT_OT_CANDIDATE_ID
    assert payload["trust_posture"]["energy_use_allowed"] is False
    assert payload["trust_posture"]["is_runtime_authority"] is False


def test_react_ot_inference_helper_script_exists():
    script_path = ROOT / "scripts" / "_react_ot_inference.py"
    assert script_path.exists(), "scripts/_react_ot_inference.py must exist"
    pilot_path = ROOT / "scripts" / "recover_ts_react_ot_seed.py"
    colab_bundle_path = ROOT / "scripts" / "prepare_react_ot_colab_bundle.py"
    colab_import_path = ROOT / "scripts" / "import_react_ot_colab_artifacts.py"
    smoke_path = ROOT / "scripts" / "run_react_ot_smoke.py"
    setup_path = ROOT / "scripts" / "setup_react_ot_env.sh"
    notebook_path = ROOT / "notebooks" / "react_ot_colab_gpu.ipynb"
    notebook_readme_path = ROOT / "notebooks" / "README.md"
    assert pilot_path.exists()
    assert colab_bundle_path.exists()
    assert colab_import_path.exists()
    assert smoke_path.exists()
    assert setup_path.exists()
    assert notebook_path.exists()
    assert notebook_readme_path.exists()
