from __future__ import annotations

from src.dft_refiner import compute_irc, DFTRefiner


def test_compute_irc_exposes_structured_api(monkeypatch):
    def _fake_generate_irc(self, ts_xyz, charge=0, spin=0, step_size=0.05):
        assert charge == 1
        assert spin == 2
        assert step_size == 0.1
        return "backward_xyz", "forward_xyz"

    def _fake_single_point(self, xyz_content, xc_method='wB97M-V', basis='def2-tzvp', charge=0, spin=0):
        return {
            "ts_xyz": -100.0,
            "backward_xyz": -101.5,
            "forward_xyz": -101.0,
        }[xyz_content]

    monkeypatch.setattr(DFTRefiner, "generate_irc", _fake_generate_irc)
    monkeypatch.setattr(DFTRefiner, "single_point", _fake_single_point)

    payload = compute_irc(
        "ts_xyz",
        charge=1,
        spin=2,
        step_size=0.1,
    )

    assert payload["path_type"] == "double_optimization_proxy"
    assert payload["backward_endpoint"] == "backward_xyz"
    assert payload["forward_endpoint"] == "forward_xyz"
    assert payload["energies"] == [-101.5, -100.0, -101.0]
    assert payload["max_energy"] == -100.0
    assert payload["max_point"] == "ts_xyz"


def test_compute_irc_requires_at_least_one_direction():
    try:
        compute_irc("ts_xyz", forward=False, backward=False)
    except ValueError as exc:
        assert "At least one IRC direction" in str(exc)
    else:
        raise AssertionError("compute_irc should reject requests with no directions enabled")