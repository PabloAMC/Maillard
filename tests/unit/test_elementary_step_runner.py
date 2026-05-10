"""Unit tests for ElementaryStepRunner with mocked refiner / GSM / xTB.

We don't exercise the real DFT/xTB stack here — we patch the heavyweight
dependencies and verify the runner's orchestration logic:

  - aggregation (max_step / sum / unknown → fallback)
  - status propagation (completed / partial / failed)
  - early-exit on TS-validation failure (n_imag, mode mismatch)
  - dry-endpoint missing → ts_guess_failed
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

from src.elementary_step_runner import (
    ElementaryStepRunner,
    StepResult,
    TargetResult,
    find_target,
    load_multistep_targets,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DUMMY_XYZ = """3
water
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0
"""


def _write_xyz(path: Path, text: str = DUMMY_XYZ) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _mk_step_spec(tmp_path: Path, step_id: str, *, with_microsolv: bool = False) -> dict:
    r = tmp_path / f"{step_id}_R.xyz"
    p = tmp_path / f"{step_id}_P.xyz"
    _write_xyz(r)
    _write_xyz(p)
    spec = {
        "step_id": step_id,
        "name": f"step {step_id}",
        # Absolute paths: REPO_ROOT / abs_path returns abs_path, so the runner
        # resolves them correctly without requiring tmp_path under REPO_ROOT.
        "dry_reactant_xyz": str(r),
        "dry_product_xyz": str(p),
        "charge": 0,
        "spin": 0,
    }
    if with_microsolv:
        spec["microsolvation"] = {
            "topology": "proton_shuttle",
            "shuttles": [{"donor_atom": 0, "acceptor_atom": 1, "h_atom": 2}],
        }
        spec["ts_mode_expected_atoms"] = [0, 1, 2]
    else:
        spec["microsolvation"] = {"topology": "none"}
    return spec


class _GSMOk:
    def __init__(self, ts_xyz: str = DUMMY_XYZ):
        self.converged = True
        self.ts_xyz = ts_xyz
        self.reason = "ok"
        self.elapsed_s = 0.1
        self.peak_index = 5
        self.n_iters = 3


class _GSMFail:
    converged = False
    ts_xyz = None
    reason = "did not converge"
    elapsed_s = 0.05
    peak_index = None
    n_iters = 1


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_load_and_find_target():
    payload = load_multistep_targets()
    target = find_target(payload, "pe_amadori_water_shuttle")
    assert target["target_id"] == "pe_amadori_water_shuttle"
    assert len(target["elementary_steps"]) >= 1


def test_find_target_missing_raises():
    payload = {"targets": [{"target_id": "x"}]}
    with pytest.raises(KeyError):
        find_target(payload, "nope")


def test_run_step_dry_endpoint_missing(tmp_path):
    refiner = MagicMock()
    runner = ElementaryStepRunner(refiner)
    spec = {
        "step_id": "s1", "name": "s1",
        "dry_reactant_xyz": "does/not/exist_R.xyz",
        "dry_product_xyz": "does/not/exist_P.xyz",
        "microsolvation": {"topology": "none"},
    }
    res = runner.run_step(spec, tmp_path / "s1")
    assert res.status == "ts_guess_failed"
    assert "missing dry endpoint" in res.reason
    assert res.barrier_kcal_mol is None


def test_run_step_gsm_failure(tmp_path):
    refiner = MagicMock()
    runner = ElementaryStepRunner(refiner)
    spec = _mk_step_spec(tmp_path, "s1")

    with patch("src.elementary_step_runner.GSMRunner") as MockGSM:
        MockGSM.return_value.run_de_gsm.return_value = _GSMFail()
        res = runner.run_step(spec, tmp_path / "s1_run")
    assert res.status == "ts_guess_failed"
    assert res.gsm["converged"] is False
    refiner.calculate_robust_barrier.assert_not_called()


def test_run_step_ts_validation_n_imag_zero(tmp_path):
    refiner = MagicMock()
    runner = ElementaryStepRunner(refiner)
    spec = _mk_step_spec(tmp_path, "s1")

    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb") as mock_probe:
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        mock_probe.return_value = {"n_imag": 0, "lowest_freq_cm": 25.0}
        res = runner.run_step(spec, tmp_path / "s1_run")
    assert res.status == "ts_validation_failed"
    assert res.n_imag == 0
    refiner.calculate_robust_barrier.assert_not_called()


def test_run_step_ts_mode_mismatch(tmp_path):
    refiner = MagicMock()
    runner = ElementaryStepRunner(refiner, ts_mode_threshold=0.3)
    spec = _mk_step_spec(tmp_path, "s1", with_microsolv=False)
    spec["ts_mode_expected_atoms"] = [0]  # but mode lives entirely on atom 2

    # 3-atom system: vector has 9 dofs; concentrate on atom 2
    v0 = np.zeros(9)
    v0[6:9] = 1.0  # atom 2

    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb") as mock_probe, \
         patch("src.elementary_step_runner.compute_xtb_ts_mode") as mock_mode:
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        mock_probe.return_value = {"n_imag": 1, "lowest_freq_cm": -800.0}
        mock_mode.return_value = v0
        res = runner.run_step(spec, tmp_path / "s1_run")

    assert res.status == "ts_validation_failed"
    assert res.mode_pass is False
    assert res.mode_concentration < 0.3
    refiner.calculate_robust_barrier.assert_not_called()


def test_run_step_completed(tmp_path):
    refiner = MagicMock()
    refiner.calculate_robust_barrier.return_value = 18.5
    runner = ElementaryStepRunner(refiner)
    spec = _mk_step_spec(tmp_path, "s1", with_microsolv=False)
    # Provide expected_atoms so the mode validation runs
    spec["ts_mode_expected_atoms"] = [0, 1, 2]
    v0 = np.ones(9) / 3.0  # uniform → concentration on [0,1,2] = 1.0

    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb") as mock_probe, \
         patch("src.elementary_step_runner.compute_xtb_ts_mode") as mock_mode:
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        mock_probe.return_value = {"n_imag": 1, "lowest_freq_cm": -1200.0}
        mock_mode.return_value = v0
        res = runner.run_step(spec, tmp_path / "s1_run")

    assert res.status == "completed"
    assert res.barrier_kcal_mol == pytest.approx(18.5)
    assert res.mode_pass is True
    assert res.mode_concentration == pytest.approx(1.0)
    refiner.calculate_robust_barrier.assert_called_once()


def test_run_target_aggregation_max(tmp_path):
    refiner = MagicMock()
    refiner.calculate_robust_barrier.side_effect = [12.0, 22.0]
    runner = ElementaryStepRunner(refiner)
    target = {
        "target_id": "demo",
        "aggregation": "max_step",
        "literature": {"barrier_kcal_mol": 20.0},
        "elementary_steps": [
            _mk_step_spec(tmp_path, "s1"),
            _mk_step_spec(tmp_path, "s2"),
        ],
    }

    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb") as mock_probe:
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        mock_probe.return_value = {"n_imag": 1, "lowest_freq_cm": -1000.0}
        result = runner.run_target(target, tmp_path / "out")

    assert result.status == "completed"
    assert result.aggregated_barrier_kcal_mol == pytest.approx(22.0)
    assert result.gap_kcal_mol == pytest.approx(2.0)
    assert len(result.steps) == 2
    assert all(s.status == "completed" for s in result.steps)


def test_run_target_partial_when_one_step_fails(tmp_path):
    refiner = MagicMock()
    refiner.calculate_robust_barrier.return_value = 15.0
    runner = ElementaryStepRunner(refiner)
    target = {
        "target_id": "demo",
        "aggregation": "max_step",
        "literature": {"barrier_kcal_mol": 20.0},
        "elementary_steps": [
            _mk_step_spec(tmp_path, "s1"),
            _mk_step_spec(tmp_path, "s2"),
        ],
    }

    # First step passes probe, second fails (n_imag=0)
    probe_seq = [
        {"n_imag": 1, "lowest_freq_cm": -900.0},
        {"n_imag": 0, "lowest_freq_cm": 30.0},
    ]
    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb", side_effect=probe_seq):
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        result = runner.run_target(target, tmp_path / "out")

    assert result.status == "partial"
    assert result.aggregated_barrier_kcal_mol is None
    assert result.steps[0].status == "completed"
    assert result.steps[1].status == "ts_validation_failed"


def test_run_target_failed_when_all_fail(tmp_path):
    refiner = MagicMock()
    runner = ElementaryStepRunner(refiner)
    target = {
        "target_id": "demo",
        "aggregation": "max_step",
        "literature": {"barrier_kcal_mol": 20.0},
        "elementary_steps": [_mk_step_spec(tmp_path, "s1")],
    }
    with patch("src.elementary_step_runner.GSMRunner") as MockGSM:
        MockGSM.return_value.run_de_gsm.return_value = _GSMFail()
        result = runner.run_target(target, tmp_path / "out")

    assert result.status == "failed"
    assert result.aggregated_barrier_kcal_mol is None
    assert result.gap_kcal_mol is None


def test_step_result_to_dict_serialisable():
    s = StepResult(step_id="x", name="x", status="completed", barrier_kcal_mol=10.0, reason="ok")
    d = s.to_dict()
    assert d["step_id"] == "x"
    assert d["barrier_kcal_mol"] == 10.0


def test_run_target_unknown_aggregation_falls_back_with_note(tmp_path):
    refiner = MagicMock()
    refiner.calculate_robust_barrier.side_effect = [10.0, 20.0]
    runner = ElementaryStepRunner(refiner)
    target = {
        "target_id": "demo",
        "aggregation": "weird",
        "literature": {"barrier_kcal_mol": 15.0},
        "elementary_steps": [
            _mk_step_spec(tmp_path, "s1"),
            _mk_step_spec(tmp_path, "s2"),
        ],
    }
    with patch("src.elementary_step_runner.GSMRunner") as MockGSM, \
         patch("src.elementary_step_runner.probe_ts_guess_xtb") as mock_probe:
        MockGSM.return_value.run_de_gsm.return_value = _GSMOk()
        mock_probe.return_value = {"n_imag": 1, "lowest_freq_cm": -900.0}
        result = runner.run_target(target, tmp_path / "out")
    assert result.status == "completed"
    assert result.aggregated_barrier_kcal_mol == pytest.approx(20.0)
    assert "unknown aggregation" in result.notes
