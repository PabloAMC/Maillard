"""Unit tests for ``src.gsm_backend``.

These tests cover the lightweight surface (XYZ helpers, GSMResult schema,
GSMRunner construction).  A full pyGSM execution is exercised in the
smoke layer (``tests/qm`` or manual scripts) because it requires a heavy
xTB stack and several seconds of compute.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.gsm_backend import (
    GSMResult,
    GSMRunner,
    _geom_to_xyz,
    _import_pygsm,
    _parse_xyz,
    _xyz_to_endpoints_file,
)


WATER_XYZ = (
    "3\n"
    "water\n"
    "O   0.0000000000   0.0000000000   0.0000000000\n"
    "H   0.7570000000   0.0000000000   0.0000000000\n"
    "H   0.0000000000   0.7570000000   0.0000000000\n"
)


def test_parse_xyz_returns_pygsm_geom_layout():
    geom = _parse_xyz(WATER_XYZ)
    assert geom == [
        ["O", 0.0, 0.0, 0.0],
        ["H", 0.757, 0.0, 0.0],
        ["H", 0.0, 0.757, 0.0],
    ]


def test_parse_xyz_rejects_truncated_block():
    with pytest.raises(ValueError):
        _parse_xyz("3\nbroken\nO 0 0 0\n")


def test_geom_to_xyz_round_trip_preserves_atoms_and_coords():
    geom = _parse_xyz(WATER_XYZ)
    xyz = _geom_to_xyz(geom, comment="round-trip")
    second = _parse_xyz(xyz)
    assert second[0][0] == "O"
    assert second[1][0] == "H"
    assert second[2][0] == "H"
    for original, replayed in zip(geom, second):
        for axis in range(1, 4):
            assert original[axis] == pytest.approx(replayed[axis], abs=1e-9)


def test_xyz_to_endpoints_file_concatenates_two_frames(tmp_path: Path):
    path = tmp_path / "endpoints.xyz"
    _xyz_to_endpoints_file(WATER_XYZ, WATER_XYZ, path)
    contents = path.read_text()
    assert contents.count("water\n") == 2  # two frames written


def test_gsmrunner_construction_does_not_import_pygsm(tmp_path: Path):
    runner = GSMRunner(work_dir=tmp_path, charge=1, spin=1, n_nodes=7,
                       max_iters=4, timeout_s=30.0)
    assert runner.charge == 1
    assert runner.spin == 1
    assert runner.n_nodes == 7
    assert runner.max_iters == 4
    assert runner.work_dir == tmp_path


def test_gsmresult_audit_dict_omits_bulky_xyz_payload(tmp_path: Path):
    result = GSMResult(
        converged=True,
        ts_xyz=WATER_XYZ,
        ts_energy_eh_xtb=-76.4,
        peak_index=5,
        n_iters=12,
        elapsed_s=1.23,
        audit_dir=str(tmp_path),
        reason="ok",
    )
    payload = result.to_audit_dict()
    assert "ts_xyz" not in payload
    assert payload["converged"] is True
    assert payload["peak_index"] == 5
    assert payload["audit_dir"] == str(tmp_path)


def test_import_pygsm_returns_required_symbols():
    """Smoke check that pyGSM is installed and importable in the test env."""
    pytest.importorskip("pyGSM")
    symbols = _import_pygsm()
    for required in ("DE_GSM", "ASELoT", "Molecule", "PES",
                     "eigenvector_follow", "Topology"):
        assert required in symbols, f"missing pyGSM symbol: {required}"


def test_gsmrunner_coerces_relative_work_dir_to_absolute(tmp_path: Path, monkeypatch):
    """Construction with a relative work_dir must yield an absolute path so
    later os.chdir into the scratch directory does not lose track of it.
    Regression test for the gsm_retry/endpoints.xyz FileNotFoundError."""
    monkeypatch.chdir(tmp_path)
    runner = GSMRunner(work_dir=Path("rel/sub/gsm_retry"))
    assert runner.work_dir.is_absolute(), runner.work_dir
    assert runner.work_dir == tmp_path / "rel" / "sub" / "gsm_retry"


def test_run_de_gsm_creates_workdir_and_endpoints_before_pygsm_call(
    tmp_path: Path, monkeypatch
):
    """`run_de_gsm` must mkdir the work_dir and write endpoints.xyz before any
    pyGSM import is attempted, even when given a non-existent relative path.
    Regression for ``FileNotFoundError: gsm_retry/endpoints.xyz``."""
    monkeypatch.chdir(tmp_path)

    captured: dict = {}

    def fake_build_gsm(self, endpoints_path):
        captured["endpoints_path"] = Path(endpoints_path)
        captured["exists"] = Path(endpoints_path).is_file()
        captured["work_dir_cwd"] = Path.cwd()
        # Short-circuit: simulate a pyGSM failure after endpoints exist.
        raise RuntimeError("stub: pyGSM not invoked in unit test")

    monkeypatch.setattr(GSMRunner, "_build_gsm", fake_build_gsm, raising=True)
    runner = GSMRunner(work_dir=Path("never_existed/gsm_retry"))
    result = runner.run_de_gsm(WATER_XYZ, WATER_XYZ)

    assert captured["exists"] is True, "endpoints.xyz must exist before pyGSM call"
    assert captured["endpoints_path"].is_absolute()
    assert captured["endpoints_path"].name == "endpoints.xyz"
    assert captured["work_dir_cwd"] == runner.work_dir, (
        "pyGSM must be invoked with cwd inside work_dir"
    )
    assert result.converged is False
    assert "stub: pyGSM not invoked" in result.reason
    # Audit dir must be the absolute work_dir, never a relative remnant.
    assert Path(result.audit_dir).is_absolute()
    assert Path(result.audit_dir) == runner.work_dir


@pytest.mark.parametrize("bad_input", [("", "x"), ("x", ""), (None, "x"), ("x", None)])
def test_run_de_gsm_rejects_empty_endpoints(tmp_path: Path, bad_input):
    runner = GSMRunner(work_dir=tmp_path / "wd")
    with pytest.raises(ValueError):
        runner.run_de_gsm(bad_input[0], bad_input[1])
