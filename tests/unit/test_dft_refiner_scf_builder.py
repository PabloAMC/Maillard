import sys
from types import SimpleNamespace

from src.dft_refiner import DFTRefiner


def test_build_mf_drops_smearing_when_soscf_is_requested(monkeypatch):
    """SOTA SCF safety contract: smearing and SOSCF/Newton must never coexist.

    PySCF's analytical gradient/Hessian paths (used by Sella and geomeTRIC)
    assume integer occupations; fractional `mo_occ` from Fermi smearing
    produces shape-mismatched boolean masks in newton_ah / density-fit code.
    `_build_mf` must drop smearing with a warning and run Newton cleanly.
    """
    refiner = DFTRefiner(solvent_name="water")
    warnings: list[str] = []

    class FakeMF:
        def __init__(self):
            self.events = []
            self.conv_tol = None
            self.max_cycle = None
            self.level_shift = 0.0
            self.damp = 0.0
            self.diis_space = 0
            self.direct_scf_tol = 0.0
            self.init_guess = None
            self.with_solvent = None

        def ddCOSMO(self):
            self.events.append("ddCOSMO")
            self.with_solvent = SimpleNamespace(eps=None)
            return self

        def newton(self):
            self.events.append("newton")
            wrapped = FakeNewtonMF(self)
            wrapped.events = self.events
            wrapped.with_solvent = self.with_solvent
            return wrapped

    class FakeNewtonMF:
        def __init__(self, base):
            self._scf = base
            self.events = []
            self.with_solvent = None
            self.max_cycle = None
            self.level_shift = 0.0

    def fake_smearing_(*args, **kwargs):
        raise AssertionError("smearing must not be applied when SOSCF is requested")

    fake_scf = SimpleNamespace(
        RHF=lambda mol: FakeMF(),
        UHF=lambda mol: FakeMF(),
        addons=SimpleNamespace(smearing_=fake_smearing_),
    )
    fake_dft = SimpleNamespace(
        RKS=lambda mol: FakeMF(),
        UKS=lambda mol: FakeMF(),
    )
    fake_pyscf = SimpleNamespace(scf=fake_scf, dft=fake_dft)

    monkeypatch.setitem(sys.modules, "pyscf", fake_pyscf)
    monkeypatch.setitem(sys.modules, "pyscf.scf", fake_scf)
    monkeypatch.setitem(sys.modules, "pyscf.dft", fake_dft)
    monkeypatch.setattr("src.dft_refiner.logger.warning", lambda message: warnings.append(message))

    mf = refiner._build_mf(
        SimpleNamespace(spin=0),
        xc_method="r2SCAN",
        use_solvent=True,
        use_soscf=True,
        smearing_sigma=0.01,
    )

    assert isinstance(mf, FakeNewtonMF)
    assert mf.events == ["ddCOSMO", "newton"]
    assert any("Refusing to combine Fermi smearing with SOSCF" in w for w in warnings)


def test_build_mf_drops_smearing_when_soscf_is_requested_open_shell(monkeypatch):
    """Same contract for UKS: smearing is dropped, Newton still runs."""
    refiner = DFTRefiner(solvent_name="water")
    warnings: list[str] = []

    class FakeMF:
        def __init__(self):
            self.events = []
            self.conv_tol = None
            self.max_cycle = None
            self.level_shift = 0.0
            self.damp = 0.0
            self.diis_space = 0
            self.direct_scf_tol = 0.0
            self.init_guess = None
            self.with_solvent = None

        def ddCOSMO(self):
            self.events.append("ddCOSMO")
            self.with_solvent = SimpleNamespace(eps=None)
            return self

        def newton(self):
            self.events.append("newton")
            wrapped = FakeNewtonMF(self)
            wrapped.events = self.events
            wrapped.with_solvent = self.with_solvent
            return wrapped

    class FakeNewtonMF:
        def __init__(self, base):
            self._scf = base
            self.events = []
            self.with_solvent = None
            self.max_cycle = None
            self.level_shift = 0.0

    def fake_smearing_(*args, **kwargs):
        raise AssertionError("smearing must not be applied when SOSCF is requested")

    fake_scf = SimpleNamespace(
        RHF=lambda mol: FakeMF(),
        UHF=lambda mol: FakeMF(),
        addons=SimpleNamespace(smearing_=fake_smearing_),
    )
    fake_dft = SimpleNamespace(
        RKS=lambda mol: FakeMF(),
        UKS=lambda mol: FakeMF(),
    )
    fake_pyscf = SimpleNamespace(scf=fake_scf, dft=fake_dft)

    monkeypatch.setitem(sys.modules, "pyscf", fake_pyscf)
    monkeypatch.setitem(sys.modules, "pyscf.scf", fake_scf)
    monkeypatch.setitem(sys.modules, "pyscf.dft", fake_dft)
    monkeypatch.setattr("src.dft_refiner.logger.warning", lambda message: warnings.append(message))

    mf = refiner._build_mf(
        SimpleNamespace(spin=1),
        xc_method="r2SCAN",
        use_solvent=True,
        use_soscf=True,
        smearing_sigma=0.01,
    )

    assert isinstance(mf, FakeNewtonMF)
    assert mf.events == ["ddCOSMO", "newton"]
    assert any("Refusing to combine Fermi smearing with SOSCF" in w for w in warnings)


def test_build_mf_keeps_newton_for_non_smearing_refinement(monkeypatch):
    refiner = DFTRefiner(solvent_name="water")

    class FakeMF:
        def __init__(self):
            self.events = []
            self.conv_tol = None
            self.max_cycle = None
            self.level_shift = 0.0
            self.damp = 0.0
            self.diis_space = 0
            self.direct_scf_tol = 0.0
            self.init_guess = None
            self.with_solvent = None

        def ddCOSMO(self):
            self.events.append("ddCOSMO")
            self.with_solvent = SimpleNamespace(eps=None)
            return self

        def newton(self):
            self.events.append("newton")
            wrapped = FakeNewtonMF(self)
            wrapped.events = self.events
            wrapped.with_solvent = self.with_solvent
            return wrapped

    class FakeNewtonMF:
        def __init__(self, base):
            self._scf = base
            self.events = []
            self.with_solvent = None
            self.max_cycle = None
            self.level_shift = 0.0

    def fake_smearing_(mf, sigma, method):
        mf.events.append(f"smearing:{sigma}:{method}")
        return mf

    fake_scf = SimpleNamespace(
        RHF=lambda mol: FakeMF(),
        UHF=lambda mol: FakeMF(),
        addons=SimpleNamespace(smearing_=fake_smearing_),
    )
    fake_dft = SimpleNamespace(
        RKS=lambda mol: FakeMF(),
        UKS=lambda mol: FakeMF(),
    )
    fake_pyscf = SimpleNamespace(scf=fake_scf, dft=fake_dft)

    monkeypatch.setitem(sys.modules, "pyscf", fake_pyscf)
    monkeypatch.setitem(sys.modules, "pyscf.scf", fake_scf)
    monkeypatch.setitem(sys.modules, "pyscf.dft", fake_dft)

    mf = refiner._build_mf(
        SimpleNamespace(spin=0),
        xc_method="r2SCAN",
        use_solvent=True,
        use_soscf=True,
        smearing_sigma=None,
    )

    assert isinstance(mf, FakeNewtonMF)
    assert mf.events == ["ddCOSMO", "newton"]


def test_fermi_smearing_guard_skips_newton_wrapped_object(monkeypatch):
    warnings = []
    refiner = DFTRefiner(solvent_name=None)

    class FakeBaseMF:
        pass

    class FakeNewtonMF:
        def __init__(self):
            self._scf = FakeBaseMF()

    def fake_smearing_(*args, **kwargs):
        raise AssertionError("smearing_ should not be called for Newton-wrapped objects")

    monkeypatch.setattr("src.dft_refiner.logger.warning", lambda message: warnings.append(message))

    mf = refiner._apply_fermi_smearing(
        FakeNewtonMF(),
        SimpleNamespace(addons=SimpleNamespace(smearing_=fake_smearing_)),
        0.01,
    )

    assert isinstance(mf, FakeNewtonMF)
    assert any("Refusing to apply Fermi smearing" in message for message in warnings)


def test_repeated_force_plateau_detects_static_geometry_with_identical_forces():
    positions = __import__("numpy").array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])
    forces = __import__("numpy").array([[0.10, 0.00, 0.00], [-0.10, 0.00, 0.00]])

    assert DFTRefiner._is_repeated_force_plateau(
        positions,
        positions.copy(),
        forces,
        forces.copy(),
        -1.2345,
        -1.2345,
    ) is True


def test_repeated_force_plateau_rejects_meaningful_force_change():
    positions = __import__("numpy").array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]])
    forces_a = __import__("numpy").array([[0.10, 0.00, 0.00], [-0.10, 0.00, 0.00]])
    forces_b = __import__("numpy").array([[0.14, 0.00, 0.00], [-0.14, 0.00, 0.00]])

    assert DFTRefiner._is_repeated_force_plateau(
        positions,
        positions.copy(),
        forces_a,
        forces_b,
        -1.2345,
        -1.2344,
    ) is False


def test_orbital_gap_instability_detects_inverted_closed_shell_gap():
    mf = SimpleNamespace(
        mol=SimpleNamespace(spin=0),
        mo_energy=__import__('numpy').array([-0.40, -0.10, -0.12, 0.20]),
        mo_occ=__import__('numpy').array([2.0, 2.0, 0.0, 0.0]),
    )

    result = DFTRefiner._orbital_gap_instability(mf)

    assert result is not None
    assert result['inverted'] is True
    assert result['near_degenerate'] is True
    assert float(result['gap_ev']) < 0.0


def test_orbital_gap_instability_ignores_open_shell():
    mf = SimpleNamespace(
        mol=SimpleNamespace(spin=1),
        mo_energy=__import__('numpy').array([-0.40, -0.10, 0.05]),
        mo_occ=__import__('numpy').array([1.0, 1.0, 0.0]),
    )

    assert DFTRefiner._orbital_gap_instability(mf) is None


def test_optimize_geometry_strategy_escalation_first_wins(monkeypatch):
    """First strategy succeeds → only one backend call."""
    refiner = DFTRefiner(solvent_name=None)
    backend_calls = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, charge, spin: xyz)

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        backend_calls.append({
            "harden_scf": harden_scf,
            "use_solvent_optimization": use_solvent_optimization,
            "use_soscf": use_soscf,
            "smearing_sigma": smearing_sigma,
        })
        return True, SimpleNamespace(), xyz_content

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    result = refiner.optimize_geometry("1\n\nH 0.0 0.0 0.0\n", is_ts=False, max_steps=5)

    assert result.converged is True
    # Only the first strategy ("standard-DIIS") was needed.
    assert len(backend_calls) == 1
    assert backend_calls[0] == {
        "harden_scf": False,
        "use_solvent_optimization": True,
        "use_soscf": False,
        "smearing_sigma": None,
    }


def test_optimize_geometry_strategy_escalation_falls_through(monkeypatch):
    """First strategy crashes with SCF gradient failure → escalates to second."""
    refiner = DFTRefiner(solvent_name=None)
    backend_calls = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, charge, spin: xyz)

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        backend_calls.append({
            "harden_scf": harden_scf,
            "use_soscf": use_soscf,
            "smearing_sigma": smearing_sigma,
        })
        if len(backend_calls) == 1:
            raise RuntimeError("Nuclear gradients of <Scanner> not converged")
        return True, SimpleNamespace(), xyz_content

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    result = refiner.optimize_geometry("1\n\nH 0.0 0.0 0.0\n", is_ts=False, max_steps=5)

    assert result.converged is True
    assert len(backend_calls) == 2
    # Second strategy is "hardened-DIIS+smearing"
    assert backend_calls[1]["harden_scf"] is True
    assert backend_calls[1]["smearing_sigma"] == 0.005


def test_optimize_geometry_partial_progress_feeds_next_strategy(monkeypatch):
    """Non-converged but non-crashing strategy feeds its best geometry forward."""
    refiner = DFTRefiner(solvent_name=None)
    backend_inputs = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, charge, spin: xyz)

    improved_xyz = "2\nimproved\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        backend_inputs.append(xyz_content)
        if len(backend_inputs) == 1:
            # First strategy: not converged, but returns partial geometry
            return False, SimpleNamespace(), improved_xyz
        # Second strategy: converges using the improved geometry
        return True, SimpleNamespace(), xyz_content

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    result = refiner.optimize_geometry("1\n\nH 0.0 0.0 0.0\n", is_ts=False, max_steps=5)

    assert result.converged is True
    assert len(backend_inputs) == 2
    # Second call received the improved geometry from the first strategy
    assert backend_inputs[1] == improved_xyz


def test_optimize_geometry_ts_uses_single_point_fallback_after_strategy_exhaustion(monkeypatch):
    refiner = DFTRefiner(solvent_name=None)
    backend_calls = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)
    # The TS-safe xTB pre-relax now runs when MLP reverts on a TS guess. For this
    # ladder-exhaustion test we keep the geometry unchanged so the assertion below
    # can compare against the raw input XYZ.
    monkeypatch.setattr(
        refiner, "_run_xtb_ts_safe_prerelax", lambda xyz, charge, spin, max_cycles=20: xyz
    )
    monkeypatch.setattr(refiner, "_setup_mol", lambda *args, **kwargs: SimpleNamespace(natm=2))

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        backend_calls.append(
            (use_soscf, smearing_sigma, use_solvent_optimization, harden_scf, level_shift)
        )
        raise RuntimeError("SCF not converged")

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))

    result = refiner.optimize_geometry("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", is_ts=True, max_steps=5)

    assert result.converged is False
    assert result.optimized_xyz == "2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
    # SOTA TS ladder: smearing+Newton is forbidden (PySCF analytical derivatives
    # break under fractional occupations). Plateau-breaking uses level_shift
    # instead, then Newton runs clean on the rotated state.
    assert backend_calls == [
        (False, None, False, False, 0.0),       # standard-DIIS
        (False, 0.005, False, True, 0.0),       # hardened-DIIS+smearing (DIIS-only)
        (False, None, False, False, 0.5),       # level-shifted-DIIS
        (True, None, False, False, 0.3),        # SOSCF/Newton + level_shift
    ]


def test_optimize_geometry_ts_invokes_xtb_safe_when_mlp_reverts(monkeypatch):
    """When all MLPs revert (return input unchanged) for a TS guess, the ladder
    must invoke the new TS-safe xTB pre-relax to relieve clashes before DFT.
    """
    refiner = DFTRefiner(solvent_name=None)
    ts_safe_calls = []

    # MLP returns the input unchanged (drift rejection scenario).
    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)

    relaxed_xyz = "2\nrelaxed\nH 0.0 0.0 0.0\nH 0.0 0.0 0.80\n"

    def fake_ts_safe(xyz, charge, spin, max_cycles=20):
        ts_safe_calls.append((xyz, charge, spin, max_cycles))
        return relaxed_xyz

    monkeypatch.setattr(refiner, "_run_xtb_ts_safe_prerelax", fake_ts_safe)

    seen_inputs = []

    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        seen_inputs.append(xyz_content)
        return True, SimpleNamespace(), xyz_content

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    refiner.optimize_geometry("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", is_ts=True, max_steps=5)

    assert len(ts_safe_calls) == 1, "TS-safe xTB pre-relax must be invoked when MLP reverts on a TS guess"
    # The DFT ladder must receive the xTB-relaxed geometry, not the input.
    assert seen_inputs[0] == relaxed_xyz


def test_optimize_geometry_ts_skips_xtb_safe_when_mlp_succeeds(monkeypatch):
    """When the MLP produces a usable update, the TS-safe xTB pre-relax is NOT
    invoked (avoid double-pre-relaxation that could destroy saddle character).
    """
    refiner = DFTRefiner(solvent_name=None)
    ts_safe_calls = []

    updated_xyz = "2\nupdated\nH 0.0 0.0 0.0\nH 0.0 0.0 0.80\n"
    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: updated_xyz)

    def fake_ts_safe(xyz, charge, spin, max_cycles=20):
        ts_safe_calls.append((xyz, charge, spin, max_cycles))
        return xyz

    monkeypatch.setattr(refiner, "_run_xtb_ts_safe_prerelax", fake_ts_safe)
    monkeypatch.setattr(
        refiner,
        "_optimize_with_pyscf_backend",
        lambda *a, **kw: (True, SimpleNamespace(), updated_xyz),
    )
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    refiner.optimize_geometry("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", is_ts=True, max_steps=5)

    assert ts_safe_calls == [], "TS-safe xTB must NOT run when MLP already produced an updated geometry"


def test_optimize_geometry_reactant_does_not_invoke_ts_safe_xtb(monkeypatch):
    """The TS-safe xTB pre-relax is gated on is_ts=True only. Reactants use the
    standard ALPB-water xTB optimization."""
    refiner = DFTRefiner(solvent_name=None)
    ts_safe_calls = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)  # MLP reverts
    monkeypatch.setattr(refiner, "_run_xtb_optimization", lambda xyz, c, s: xyz)

    def fake_ts_safe(xyz, charge, spin, max_cycles=20):
        ts_safe_calls.append((xyz, charge, spin, max_cycles))
        return xyz

    monkeypatch.setattr(refiner, "_run_xtb_ts_safe_prerelax", fake_ts_safe)
    monkeypatch.setattr(
        refiner,
        "_optimize_with_pyscf_backend",
        lambda *a, **kw: (True, SimpleNamespace(), "1\n\nH 0.0 0.0 0.0\n"),
    )
    monkeypatch.setattr(refiner, "_build_mf", lambda *args, **kwargs: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *args, **kwargs: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *args, **kwargs: ([100.0], -0.9, -0.8))

    refiner.optimize_geometry("1\n\nH 0.0 0.0 0.0\n", is_ts=False, max_steps=5)

    assert ts_safe_calls == [], "TS-safe xTB must never run for reactants (is_ts=False)"


def test_secondary_prerelax_backend_logs_when_unavailable(monkeypatch):
    """When the secondary MLP (UBio) cannot be imported, the failure must be
    surfaced (logged + cached) so silent fallbacks are visible in the run log."""
    refiner = DFTRefiner(solvent_name=None)
    refiner._secondary_prerelax_mlp = None
    refiner._secondary_prerelax_mlp_error = None

    # Force the import inside _secondary_prerelax_backend to fail.
    import builtins
    real_import = builtins.__import__

    def failing_import(name, *args, **kwargs):
        if name == "src.ubio_optimizer" or name.endswith("ubio_optimizer"):
            raise ImportError("No module named 'torch_cluster'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", failing_import)

    # Capture warnings emitted by the module logger directly.
    from src import dft_refiner as dft_refiner_module
    captured = []
    monkeypatch.setattr(
        dft_refiner_module.logger,
        "warning",
        lambda msg, *a, **kw: captured.append(str(msg) % a if a else str(msg)),
    )

    result = refiner._secondary_prerelax_backend()

    assert result is None
    assert refiner._secondary_prerelax_mlp_error is not None
    assert "torch_cluster" in refiner._secondary_prerelax_mlp_error
    assert any("MLP FALLBACK" in m for m in captured), (
        f"Expected 'MLP FALLBACK' warning, captured: {captured}"
    )
    assert any("torch_cluster" in m for m in captured)

    # Subsequent calls must short-circuit without re-emitting the warning.
    captured.clear()
    second = refiner._secondary_prerelax_backend()
    assert second is None
    assert captured == []
