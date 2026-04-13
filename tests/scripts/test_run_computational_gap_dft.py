import json
import sys
import types
from pathlib import Path

import pytest

from scripts import run_computational_gap_dft as runner


def test_runner_overwrites_stale_execution_payload(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "hexanal_radical_quench",
                        "reaction_key": "hexanal_radical_quench",
                        "status": "ready_for_dft",
                        "reactant_path": "missing/reactant.xyz",
                        "ts_guess_path": "missing/ts.xyz",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    output_path.write_text(
        json.dumps({"jobs": [{"target_id": "stale_failed_job", "status": "failed"}]}),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_dft.py",
            "--manifest",
            str(manifest_path),
            "--output",
            str(output_path),
        ],
    )

    assert runner.main() == 0

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["jobs"] == [
        {
            "target_id": "hexanal_radical_quench",
            "reaction_key": "hexanal_radical_quench",
            "status": "seed_required",
            "execute_mode": False,
            "fast_mode": False,
            "method_chain": runner.DEFAULT_METHOD_CHAIN,
            "promotion_ready": False,
        }
    ]


def test_runner_persists_interrupted_status(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    reactant_path = tmp_path / "reactant.xyz"
    ts_path = tmp_path / "ts.xyz"
    reactant_path.write_text("1\nreactant\nH 0.0 0.0 0.0\n", encoding="utf-8")
    ts_path.write_text("1\nts\nH 0.0 0.0 0.1\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "hexanal_radical_quench",
                        "reaction_key": "hexanal_radical_quench",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 1,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    class FakeRefiner:
        def __init__(self, **_: object) -> None:
            pass

        def calculate_barrier(self, *_: object, **__: object) -> float:
            raise KeyboardInterrupt()

    monkeypatch.setitem(sys.modules, "src.dft_refiner", types.SimpleNamespace(DFTRefiner=FakeRefiner))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_dft.py",
            "--manifest",
            str(manifest_path),
            "--output",
            str(output_path),
            "--target",
            "hexanal_radical_quench",
            "--execute",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        runner.main()

    assert exc_info.value.code == 130
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["summary"]["interrupted_count"] == 1
    assert payload["jobs"][0]["status"] == "interrupted"
    assert payload["jobs"][0]["error"] == "Interrupted by user"


def test_runner_marks_job_running_before_heavy_compute(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    reactant_path = tmp_path / "reactant.xyz"
    ts_path = tmp_path / "ts.xyz"
    reactant_path.write_text("1\nreactant\nH 0.0 0.0 0.0\n", encoding="utf-8")
    ts_path.write_text("1\nts\nH 0.0 0.0 0.1\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "hexanal_radical_quench",
                        "reaction_key": "hexanal_radical_quench",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 1,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    class FakeRefiner:
        def __init__(self, **_: object) -> None:
            pass

        def calculate_barrier(self, *_: object, **__: object) -> float:
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            assert payload["jobs"][0]["status"] == "running"
            assert payload["jobs"][0]["started_at"]
            assert payload["jobs"][0]["updated_at"]
            assert payload["jobs"][0]["elapsed_seconds"] == 0
            assert payload["summary"]["running_count"] == 1
            assert payload["summary"]["active_targets"] == ["hexanal_radical_quench"]
            raise KeyboardInterrupt()

    monkeypatch.setitem(sys.modules, "src.dft_refiner", types.SimpleNamespace(DFTRefiner=FakeRefiner))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_dft.py",
            "--manifest",
            str(manifest_path),
            "--output",
            str(output_path),
            "--target",
            "hexanal_radical_quench",
            "--execute",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        runner.main()

    assert exc_info.value.code == 130


def test_runner_passes_string_family_metadata(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    reactant_path = tmp_path / "reactant.xyz"
    ts_path = tmp_path / "ts.xyz"
    reactant_path.write_text("1\nreactant\nH 0.0 0.0 0.0\n", encoding="utf-8")
    ts_path.write_text("1\nts\nH 0.0 0.0 0.1\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "hexanal_radical_quench",
                        "reaction_key": "hexanal_radical_quench",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 1,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    class FakeRefiner:
        def __init__(self, **_: object) -> None:
            pass

        def calculate_barrier(self, *_: object, **kwargs: object) -> float:
            reaction_meta = kwargs.get("reaction_meta", {})
            assert isinstance(reaction_meta.get("family"), str)
            assert reaction_meta.get("family") == "hexanal_radical_quench"
            return 17.88

    monkeypatch.setitem(sys.modules, "src.dft_refiner", types.SimpleNamespace(DFTRefiner=FakeRefiner))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_computational_gap_dft.py",
            "--manifest",
            str(manifest_path),
            "--output",
            str(output_path),
            "--target",
            "hexanal_radical_quench",
            "--execute",
        ],
    )

    assert runner.main() == 0
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["summary"]["completed_count"] == 1
    assert payload["jobs"][0]["status"] == "completed"