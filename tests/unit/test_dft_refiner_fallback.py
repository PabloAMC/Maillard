from types import SimpleNamespace

from src.dft_refiner import DFTRefiner


MINIMAL_XYZ = """1

H 0.0 0.0 0.0
"""


def test_optimize_geometry_retries_with_soscf_after_scf_failure(monkeypatch):
    refiner = DFTRefiner(solvent_name=None)
    # Disable MLP prerelax so failure path is exercised
    refiner._prerelax_mlp = None
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
        use_soscf=False,
    ):
        calls.append((use_soscf, use_solvent_optimization))
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
    # First call: normal SCF, no SOSCF, no solvent
    assert calls[0] == (False, False)
    # Second call: SOSCF + solvent for non-TS retry
    assert calls[1] == (True, True)


def test_detects_scf_gradient_failure_pattern():
    refiner = DFTRefiner(solvent_name=None)

    assert refiner._is_scf_gradient_failure(RuntimeError("Nuclear gradients of <scanner> not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("SCF not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("some unrelated error")) is False


def test_mlp_prerelax_called_for_non_ts(monkeypatch):
    """MLP pre-relaxation runs for reactants and feeds into the PySCF backend."""
    refiner = DFTRefiner(solvent_name=None)
    prerelax_calls = []

    def fake_mlp_prerelax(xyz_content):
        prerelax_calls.append(xyz_content)
        return "2\nrelaxed\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *, harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
    ):
        return True, object(), xyz_content

    monkeypatch.setattr(refiner, "_mlp_prerelax", fake_mlp_prerelax)
    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_run_solvent_scf_with_retry", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert len(prerelax_calls) == 1


def test_mlp_prerelax_skipped_for_ts(monkeypatch):
    """MLP pre-relaxation is skipped for TS optimizations."""
    refiner = DFTRefiner(solvent_name=None)
    prerelax_calls = []

    def fake_mlp_prerelax(xyz_content):
        prerelax_calls.append(xyz_content)
        return xyz_content

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *, harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
    ):
        return True, object(), xyz_content

    monkeypatch.setattr(refiner, "_mlp_prerelax", fake_mlp_prerelax)
    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_run_solvent_scf_with_retry", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    refiner.optimize_geometry(MINIMAL_XYZ, is_ts=True, max_steps=5)

    assert len(prerelax_calls) == 0