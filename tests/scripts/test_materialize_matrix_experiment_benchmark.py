import json
import sys
from pathlib import Path

from scripts.generators import materialize_matrix_experiment_benchmark as materializer


ROOT = Path(__file__).resolve().parents[2]


def test_materialize_script_writes_benchmark_json(tmp_path, monkeypatch):
    input_payload = ROOT / "data/protocols/pea_iso_protocol_pilot_intake.yaml"
    output_path = tmp_path / "pilot_benchmark.json"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "materialize_matrix_experiment_benchmark.py",
            str(input_payload),
            str(output_path),
        ],
    )

    assert materializer.main() == 0
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["benchmark_id"] == "pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026"
    assert payload["evidence_class"] == "diagnostic_only"
    assert payload["source_metadata"]["evidence_class"] == "diagnostic_only"