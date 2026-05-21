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
        label="default",
        harden_scf=False,
        use_solvent_optimization=False,
        use_soscf=False,
        smearing_sigma=None,
        level_shift=0.0, **kwargs,
    ):
        calls.append((use_soscf, harden_scf, smearing_sigma))
        if len(calls) == 1:
            raise RuntimeError("Nuclear gradients of <pyscf.scf.hf.RKS_Scanner object at 0x0> not converged")
        return True, object(), MINIMAL_XYZ

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert result.converged is True
    assert len(calls) == 2
    # First call: standard-DIIS (no SOSCF, no hardening, no smearing)
    assert calls[0] == (False, False, None)
    # Second call: hardened-DIIS+smearing
    assert calls[1] == (False, True, 0.005)


def test_detects_scf_gradient_failure_pattern():
    refiner = DFTRefiner(solvent_name=None)

    assert refiner._is_scf_gradient_failure(RuntimeError("Nuclear gradients of <scanner> not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("SCF not converged")) is True
    assert refiner._is_scf_gradient_failure(RuntimeError("some unrelated error")) is False


def test_non_ts_runtime_exception_uses_best_available_geometry(monkeypatch):
    refiner = DFTRefiner(solvent_name=None)
    refiner._prerelax_mlp = None
    opt_calls = []
    fallback_calls = []

    def fake_opt_backend(
        xyz_content,
        charge,
        spin,
        is_ts,
        max_steps,
        eff_use_explicit,
        n_atoms_solute,
        *,
        label="default",
        harden_scf=False,
        use_solvent_optimization=False,
        use_soscf=False,
        smearing_sigma=None,
        level_shift=0.0,
        **kwargs,
    ):
        opt_calls.append(label)
        raise RuntimeError("geomeTRIC callback exploded")

    def fake_fallback(xyz, *, charge, spin, basis, label):
        fallback_calls.append((xyz, label))
        return False, SimpleNamespace(natm=1), xyz

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, charge, spin: xyz)
    monkeypatch.setattr(refiner, "_best_available_geometry_fallback", fake_fallback)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert opt_calls == ["default"]
    assert fallback_calls == [(MINIMAL_XYZ, "Reactant refinement")]
    assert result.converged is False
    assert result.optimized_xyz == MINIMAL_XYZ


def test_mlp_prerelax_called_for_non_ts(monkeypatch):
    """MLP pre-relaxation runs for reactants and feeds into the PySCF backend."""
    refiner = DFTRefiner(solvent_name=None)
    prerelax_calls = []

    def fake_mlp_prerelax(xyz_content):
        prerelax_calls.append(xyz_content)
        return "2\nrelaxed\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *, label="default",
        harden_scf=False, use_solvent_optimization=False,
        use_soscf=False, smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        return True, object(), xyz_content

    monkeypatch.setattr(refiner, "_mlp_prerelax", fake_mlp_prerelax)
    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    assert len(prerelax_calls) == 1


def test_mlp_prerelax_runs_for_ts(monkeypatch):
    """MLP pre-relaxation also runs for TS (geometry cleanup)."""
    refiner = DFTRefiner(solvent_name=None)
    prerelax_calls = []

    def fake_mlp_prerelax(xyz_content):
        prerelax_calls.append(xyz_content)
        return xyz_content

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *, label="default",
        harden_scf=False, use_solvent_optimization=False,
        use_soscf=False, smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        return True, object(), xyz_content

    monkeypatch.setattr(refiner, "_mlp_prerelax", fake_mlp_prerelax)
    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)

    refiner.optimize_geometry(MINIMAL_XYZ, is_ts=True, max_steps=5)

    assert len(prerelax_calls) == 1

    assert len(prerelax_calls) == 1


# ---------------------------------------------------------------------------
# Coverage gap closure: items 4 from the DFT Ladder Bulletproofing backlog —
# non-TS best-geometry SP fallback, authoritative resume, SP exception fallback.
# ---------------------------------------------------------------------------

_FALLBACK_XYZ = """1

H 0.5 0.5 0.5
"""


def test_non_ts_sp_uses_best_available_fallback_geometry(monkeypatch):
    """When the DFT ladder fails on a non-TS, the high-level single-point and
    Hessian must be evaluated on the *fallback* geometry produced by
    `_best_available_geometry_fallback`, not the original input geometry.

    Regression for: maintenance backlog item "Add focused unit coverage for
    non-TS best-geometry SP fallback" (tasks/todo.md → DFT Ladder
    Bulletproofing). Previous coverage in
    `test_non_ts_runtime_exception_uses_best_available_geometry` only checked
    that the fallback was called — not that downstream stages consumed its
    output. A regression that fed the original xyz into the SP would have
    silently passed.
    """
    refiner = DFTRefiner(solvent_name=None)
    refiner._prerelax_mlp = None

    fallback_mol = SimpleNamespace(natm=1)

    def fake_opt_backend(*args, **kwargs):
        raise RuntimeError("geomeTRIC blew up; not an SCF gradient failure")

    def fake_fallback(xyz, *, charge, spin, basis, label):
        # Simulate: best-available geometry differs from the original input.
        return False, fallback_mol, _FALLBACK_XYZ

    sp_calls = []
    build_mf_calls = []
    hessian_calls = []

    def fake_single_point(xyz, **kwargs):
        sp_calls.append(xyz)
        return -0.7

    def fake_build_mf(mol, **kwargs):
        build_mf_calls.append(mol)
        return SimpleNamespace(e_tot=-0.7)

    def fake_hessian(mf):
        hessian_calls.append(mf)
        return ([100.0], -0.6, -0.5)

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, c, s: xyz)
    monkeypatch.setattr(refiner, "_best_available_geometry_fallback", fake_fallback)
    monkeypatch.setattr(refiner, "single_point", fake_single_point)
    monkeypatch.setattr(refiner, "_build_mf", fake_build_mf)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", fake_hessian)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    # The single-point and Hessian/Thermo paths must consume the *fallback*
    # xyz, not the original `MINIMAL_XYZ`. The Hessian receives the same `mf`
    # object that was built from the fallback xyz.
    assert sp_calls == [_FALLBACK_XYZ]
    assert len(build_mf_calls) == 1
    assert len(hessian_calls) == 1
    assert result.optimized_xyz == _FALLBACK_XYZ
    assert result.converged is False
    assert result.energy_hartree == -0.7


def test_geometry_ready_resume_skips_mlp_xtb_and_dft_ladder(monkeypatch, tmp_path):
    """When a checkpoint marks the geometry as ready for post-processing, the
    refiner must skip MLP pre-relaxation, xTB optimization, and the entire
    DFT SCF strategy ladder — and go straight to single-point + Hessian.

    Authoritative-resume coverage for the DFT Ladder Bulletproofing backlog.
    """
    refiner = DFTRefiner(solvent_name=None)

    mlp_calls: list = []
    xtb_calls: list = []
    backend_calls: list = []

    monkeypatch.setattr(
        refiner,
        "_mlp_prerelax",
        lambda xyz: (mlp_calls.append(xyz), xyz)[1],
    )
    monkeypatch.setattr(
        refiner,
        "_run_xtb_optimization",
        lambda xyz, c, s: (xtb_calls.append(xyz), xyz)[1],
    )

    def fake_opt_backend(*args, **kwargs):
        backend_calls.append(args)
        raise AssertionError("DFT ladder must not be entered on geometry-ready resume.")

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_setup_mol", lambda *a, **kw: SimpleNamespace(natm=1))
    monkeypatch.setattr(refiner, "_build_mf", lambda *a, **kw: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda xyz, **kw: -0.42)
    monkeypatch.setattr(
        refiner,
        "_run_hessian_and_thermo",
        lambda mf: ([200.0], -0.41, -0.40),
    )

    refiner._set_optimization_checkpoint_context(
        checkpoint_dir=tmp_path,
        phase_prefix="reactant",
        resume_state={
            "geometry_ready_for_postprocessing": True,
            "mlp_completed": True,
            "xtb_completed": True,
            "converged": True,
        },
    )

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    # No re-relaxation work: MLP, xTB, and the DFT ladder must all be skipped.
    assert mlp_calls == [], "MLP pre-relax must not run on geometry-ready resume"
    assert xtb_calls == [], "xTB optimization must not run on geometry-ready resume"
    assert backend_calls == [], "DFT SCF ladder must not be entered on geometry-ready resume"
    # Post-processing still produces a result with the resumed geometry.
    assert result.optimized_xyz == MINIMAL_XYZ
    assert result.converged is True
    assert result.energy_hartree == -0.42


def test_single_point_recovers_from_solvent_scf_exception(monkeypatch):
    """`single_point` must catch a runtime exception from the primary SCF
    (with implicit solvent) and retry without solvent + hardened SCF.

    Coverage for "Harden the high-level single-point stage against runtime
    exceptions, not only SCF non-convergence" in the DFT Ladder
    Bulletproofing backlog. A regression that re-raised the primary
    exception would have lost a valid optimized geometry.
    """
    refiner = DFTRefiner(solvent_name="water")

    monkeypatch.setattr(refiner, "_setup_mol", lambda *a, **kw: SimpleNamespace(natm=1))
    monkeypatch.setattr(refiner, "_log_spin_integrity", lambda mf: None)

    build_calls: list = []

    class _FailingSCF:
        converged = False

        def scf(self):
            raise RuntimeError("ddCOSMO crashed on weird geometry")

    class _RecoveredSCF:
        converged = True

        def scf(self):
            return -76.123

    def fake_build_mf(mol, *, xc_method, use_solvent, harden_scf=False, **kwargs):
        build_calls.append(
            {"use_solvent": use_solvent, "harden_scf": harden_scf}
        )
        if use_solvent and not harden_scf:
            return _FailingSCF()
        # Recovery path: vacuum + hardened SCF.
        return _RecoveredSCF()

    monkeypatch.setattr(refiner, "_build_mf", fake_build_mf)

    energy = refiner.single_point(MINIMAL_XYZ)

    assert energy == -76.123
    # Two `_build_mf` calls: primary (solvent, no hardening), then recovery
    # (vacuum, hardened).
    assert build_calls == [
        {"use_solvent": True, "harden_scf": False},
        {"use_solvent": False, "harden_scf": True},
    ]