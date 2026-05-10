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


def test_runner_blocks_bad_ts_guess_before_heavy_compute(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    reactant_path = tmp_path / "reactant.xyz"
    ts_path = tmp_path / "ts.xyz"
    reactant_path.write_text("1\nreactant\nH 0.0 0.0 0.0\n", encoding="utf-8")
    ts_path.write_text("1\nts\nO 0.0 0.0 0.1\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "aa_ring_open_dicarbonyl",
                        "reaction_key": "aa_ring_open_dicarbonyl",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 0,
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
            raise AssertionError("DFT should not start when preflight blocks the TS guess")

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
            "aa_ring_open_dicarbonyl",
            "--execute",
        ],
    )

    assert runner.main() == 0
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["summary"]["blocked_count"] == 1
    assert payload["jobs"][0]["status"] == "blocked_bad_ts_guess"
    assert payload["jobs"][0]["restart_recommended"] is True
    assert "element_sequence_mismatch" in payload["jobs"][0]["preflight"]["blockers"]


def test_runner_records_progress_phases(tmp_path, monkeypatch):
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
            progress_callback = kwargs["progress_callback"]
            progress_callback("reactant_geometry_optimization", {"label": "Reactant geometry optimization"})
            payload = json.loads(output_path.read_text(encoding="utf-8"))
            assert payload["jobs"][0]["phase"] == "reactant_geometry_optimization"
            progress_callback("ts_single_point", {"label": "TS high-level single-point"})
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
    assert payload["jobs"][0]["status"] == "completed"
    assert payload["jobs"][0]["phase"] == "completed"
    phases = [entry["phase"] for entry in payload["jobs"][0]["phase_history"]]
    assert "preflight_passed" in phases
    assert "reactant_geometry_optimization" in phases
    assert "ts_single_point" in phases


def test_runner_repairs_bad_ts_guess_from_xtb_frames(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    runner_dir = tmp_path / "data/geometries/xtb_inputs/aa_ring_open_dicarbonyl"
    runner_dir.mkdir(parents=True)
    reactant_path = runner_dir / "reactant.xyz"
    ts_path = runner_dir / "xtbpath_ts.xyz"
    frame_path = runner_dir / "xtbpath_10.xyz"
    reactant_path.write_text("2\nreactant\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7\n", encoding="utf-8")
    ts_path.write_text("2\nts\nH 0.0 0.0 0.0\nO 0.0 0.0 0.7\n", encoding="utf-8")
    frame_path.write_text("2\nframe\nH 0.0 0.0 0.0\nH 0.0 0.0 0.8\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "aa_ring_open_dicarbonyl",
                        "reaction_key": "aa_ring_open_dicarbonyl",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 0,
                    }
                ]
            }
        ),
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
            "--target",
            "aa_ring_open_dicarbonyl",
            "--preflight-only",
        ],
    )

    assert runner.main() == 0
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["jobs"][0]["status"] == "ready_for_dft"
    assert payload["jobs"][0]["ts_guess_path_used"].endswith("xtbpath_10.xyz")
    assert payload["jobs"][0]["preflight_repair"]["original_blockers"] == ["element_sequence_mismatch"]


def test_runner_uses_forming_bond_recovery_before_blocking(tmp_path, monkeypatch):
    manifest_path = tmp_path / "manifest.json"
    output_path = tmp_path / "execution.json"
    reactant_path = tmp_path / "reactant.xyz"
    ts_path = tmp_path / "ts.xyz"
    reactant_path.write_text("2\nreactant\nC 0.0 0.0 0.0\nN 0.0 0.0 3.0\n", encoding="utf-8")
    ts_path.write_text("2\nts\nC 0.0 0.0 0.0\nO 0.0 0.0 0.7\n", encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "jobs": [
                    {
                        "target_id": "lysinoalanine_crosslink",
                        "reaction_key": "lysinoalanine_crosslink",
                        "status": "ready_for_dft",
                        "reactant_path": str(reactant_path),
                        "ts_guess_path": str(ts_path),
                        "charge": 0,
                        "spin": 0,
                        "forming_bond_atoms": [0, 1],
                        "forming_bond": {"available": True, "atom_indices_zero_based": [0, 1]},
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    class FakeRefiner:
        def __init__(self, **_: object) -> None:
            pass

        def recover_ts_guess_from_forming_bond(self, *_: object, **__: object) -> dict[str, object]:
            return {
                "used": True,
                "strategy": "forming_bond_relaxed_scan",
                "ts_xyz": "2\nscan\nC 0.0 0.0 0.0\nN 0.0 0.0 1.4\n",
                "imaginary_mode_count": 1,
                "initial_imaginary_mode_count": 3,
                "scan": {"scan_point_count": 5},
                "steric_repair": None,
            }

        def calculate_barrier(self, *_: object, **__: object) -> float:
            return 12.34

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
            "lysinoalanine_crosslink",
            "--execute",
        ],
    )

    assert runner.main() == 0
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["jobs"][0]["status"] == "completed"
    assert payload["jobs"][0]["ts_guess_path_used"] == "forming_bond_relaxed_scan"
    assert payload["jobs"][0]["forming_bond_recovery"]["used"] is True
