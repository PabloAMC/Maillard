from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mlp_adoption_contract import MLPModelCandidate  # noqa: E402
from src.mlp_backend_adapters import build_candidate_adapter  # noqa: E402


def _candidate(model_family: str, proposed_role: str) -> MLPModelCandidate:
    backend_locator = None
    if model_family == "mace_omol":
        backend_locator = "builtin:large"
    elif model_family == "mace_mp":
        backend_locator = "builtin:small"
    return MLPModelCandidate(
        candidate_id=f"{model_family}_{proposed_role}",
        model_family=model_family,
        model_name="small",
        chemistry_domain="general_organic",
        materials_first=False,
        proposed_role=proposed_role,
        target_motif_families=["carbonyl_fragmentation"],
        benchmark_visible_gap="test gap",
        expected_speedup=5.0,
        likely_failure_modes=["benchmark_gap"],
        fallback_comparator="fallback",
        backend_locator=backend_locator,
    )


def test_adapter_registry_returns_explicit_adapters_for_chemistry_first_shortlist():
    mace_omol = build_candidate_adapter(_candidate("mace_omol", "geom_preopt"))
    aimnet2 = build_candidate_adapter(_candidate("aimnet2", "ts_initialization"))
    orbmol = build_candidate_adapter(_candidate("orbmol", "ts_initialization"))

    assert mace_omol.model_family == "mace_omol"
    assert aimnet2.model_family == "aimnet2"
    assert orbmol.model_family == "orbmol"


def test_unsupported_adapter_reports_reason():
    adapter = build_candidate_adapter(_candidate("unknown_model", "geom_preopt"))
    probe = adapter.probe_availability()

    assert probe.available is False
    assert "No explicit backend adapter" in probe.reason


def test_ts_backend_plugin_locator_enables_fast_mock_probe(monkeypatch):
    monkeypatch.setenv("MAILLARD_AIMNET2_BACKEND", "python:tests.helpers.mock_ts_backend")
    candidate = MLPModelCandidate(
        candidate_id="aimnet2_shortlist",
        model_family="aimnet2",
        model_name="aimnet2_wb97m",
        chemistry_domain="organic_qm_surrogate",
        materials_first=False,
        proposed_role="ts_initialization",
        target_motif_families=["proton_transfer_rearrangement"],
        benchmark_visible_gap="test gap",
        expected_speedup=5.0,
        likely_failure_modes=["benchmark_gap"],
        fallback_comparator="fallback",
        backend_locator="env:MAILLARD_AIMNET2_BACKEND",
    )

    adapter = build_candidate_adapter(candidate)
    probe = adapter.probe_availability()
    prepared = adapter.prepare_ts_seed("3\nHCN\nH 0 1 0\nC -1 0 0\nN 1 0 0\n")

    assert probe.available is True
    assert "mock backend ready" in probe.reason
    assert "HCN" in prepared