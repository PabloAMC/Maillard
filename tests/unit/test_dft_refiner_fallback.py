from types import SimpleNamespace

from src.dft_refiner import DFTRefiner


MINIMAL_XYZ = """1

H 0.0 0.0 0.0
"""


def test_optimize_geometry_retries_after_scf_gradient_failure(monkeypatch):
    refiner = DFTRefiner(solvent_name=None)
    calls = []

    def fake_opt_backend(
        xyz_content,
        charge,
        spin,
        is_ts,
        max_steps,
        eff_use_explicit,
        n_atoms_solute,
        *,
        harden_scf=False,
        use_solvent_optimization=False,
    ):
        calls.append((xyz_content, harden_scf, use_solvent_optimization))
        if len(calls) == 1:
            raise RuntimeError("Nuclear gradients of <pyscf.scf.hf.RKS_Scanner object at 0x0> not converged")
        return True, object(), MINIMAL_XYZ

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_preoptimize_seed_xyz", lambda *args, **kwargs: MINIMAL_XYZ)
    monkeypatch.setattr(refiner, "_run_solvent_scf_with_retry", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert result.converged is True
    assert len(calls) == 2
    assert calls[0][1:] == (False, False)
    assert calls[1][1:] == (True, False)


def test_detects_scf_gradient_failure_pattern():
    refiner = DFTRefiner(solvent_name=None)

    assert refiner._is_scf_gradient_failure(RuntimeError("Nuclear gradients of <scanner> not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("SCF not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("some unrelated error")) is False


def test_optimize_geometry_uses_solvent_retry_after_second_scf_gradient_failure(monkeypatch):
    refiner = DFTRefiner(solvent_name=None)
    calls = []

    def fake_opt_backend(
        xyz_content,
        charge,
        spin,
        is_ts,
        max_steps,
        eff_use_explicit,
        n_atoms_solute,
        *,
        harden_scf=False,
        use_solvent_optimization=False,
    ):
        calls.append((harden_scf, use_solvent_optimization))
        if len(calls) < 3:
            raise RuntimeError("Nuclear gradients of <pyscf.scf.hf.RKS_Scanner object at 0x0> not converged")
        return True, object(), MINIMAL_XYZ

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_preoptimize_seed_xyz", lambda *args, **kwargs: MINIMAL_XYZ)
    monkeypatch.setattr(refiner, "_run_solvent_scf_with_retry", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert result.converged is True
    assert calls == [
        (False, False),
        (True, False),
        (True, True),
    ]