from __future__ import annotations

import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import assess_matrix_benchmark_evidence, build_benchmark_index, load_benchmark
from src.matrix_experiment_intake import load_matrix_experiment_intake, materialize_matrix_experiment_benchmark


PILOT_CASES = [
    (
        ROOT / "data" / "protocols" / "pea_iso_protocol_pilot_intake.yaml",
        ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
        "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026",
    ),
    (
        ROOT / "data" / "protocols" / "soy_iso_protocol_pilot_intake.yaml",
        ROOT / "data" / "benchmarks" / "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026.json",
        "soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026",
    ),
]


@pytest.mark.parametrize(("pilot_intake", "pilot_benchmark", "benchmark_id"), PILOT_CASES)
def test_protocol_pilot_benchmark_matches_materialized_intake_payload(pilot_intake, pilot_benchmark, benchmark_id):
    intake = load_matrix_experiment_intake(pilot_intake)
    expected = materialize_matrix_experiment_benchmark(intake)
    actual = load_benchmark(pilot_benchmark)

    assert actual == expected
    assert actual["benchmark_id"] == benchmark_id
    assert actual["source_metadata"]["origin"] == "internal_experiment"
    assert "measured_volatiles" in actual


@pytest.mark.parametrize(("_pilot_intake", "pilot_benchmark", "benchmark_id"), PILOT_CASES)
def test_protocol_pilot_benchmark_is_usable_but_not_external_promotion_evidence(_pilot_intake, pilot_benchmark, benchmark_id):
    entries = build_benchmark_index([pilot_benchmark])
    assert len(entries) == 1
    assert entries[0].benchmark_id == benchmark_id
    assert entries[0].execution_path == "matrix_precursor_augmented"

    evidence = assess_matrix_benchmark_evidence(pilot_benchmark)
    assert evidence.reference_signal_origin == "measured_volatiles"
    assert evidence.external_data_status == "quantitative_unspecified_origin"
    assert evidence.promotable is False