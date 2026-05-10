"""Tests for TSOptimizer watchdog (timeout, stall detection) and
BuiltinPySCFCalc factory-based MF construction."""

import pytest
import numpy as np
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from src.ts_optimizer import TSOptimizer, TSOptimizationStall


# ---------------------------------------------------------------------------
# TSOptimizer: timeout & stall
# ---------------------------------------------------------------------------

class _FakeAtoms:
    """Minimal stand-in for ASE Atoms with a calculator."""
    def __init__(self, forces):
        self._forces = forces
        self.calc = None

    def get_forces(self):
        return self._forces

    @property
    def positions(self):
        return np.zeros((2, 3))


class _FakeSella:
    """Fake Sella dynamics that calls observers on each step."""
    def __init__(self, atoms, logfile=None, trajectory=None, **kwargs):
        self.atoms = atoms
        self._observers = []

    def attach(self, fn):
        self._observers.append(fn)

    def run(self, fmax, steps):
        for _ in range(steps):
            for obs in self._observers:
                obs()


def test_ts_optimizer_timeout(monkeypatch, tmp_path):
    """Watchdog raises TSOptimizationStall when wall-time is exceeded."""
    monkeypatch.setattr("src.ts_optimizer.Sella", _FakeSella)

    # Force time.monotonic() to jump past the timeout
    call_count = 0
    def fake_monotonic():
        nonlocal call_count
        call_count += 1
        # First call is the start time (0), second call triggers timeout
        return 0.0 if call_count <= 1 else 99999.0

    monkeypatch.setattr("src.ts_optimizer.time.monotonic", fake_monotonic)

    forces = np.ones((2, 3))
    atoms = _FakeAtoms(forces)
    atoms.calc = SimpleNamespace(results={"forces": forces})

    ts_opt = TSOptimizer(fmax=0.01, max_steps=100)
    with pytest.raises(TSOptimizationStall, match="exceeded"):
        ts_opt.find_ts(atoms, atoms.calc, timeout_seconds=10)


def test_ts_optimizer_stall_detection(monkeypatch, tmp_path):
    """Watchdog raises TSOptimizationStall when forces plateau."""
    monkeypatch.setattr("src.ts_optimizer.Sella", _FakeSella)

    # Fixed time that never hits timeout
    monkeypatch.setattr("src.ts_optimizer.time.monotonic", lambda: 0.0)

    # Constant forces → stall after stall_window steps
    flat_forces = np.full((2, 3), 1.0)
    atoms = _FakeAtoms(flat_forces)
    atoms.calc = SimpleNamespace(results={"forces": flat_forces})

    ts_opt = TSOptimizer(fmax=0.001, max_steps=200)
    with pytest.raises(TSOptimizationStall, match="stalled"):
        ts_opt.find_ts(atoms, atoms.calc, stall_window=5)


def test_ts_optimizer_converges_normally(monkeypatch):
    """No error when optimization runs within limits."""
    monkeypatch.setattr("src.ts_optimizer.Sella", _FakeSella)
    monkeypatch.setattr("src.ts_optimizer.time.monotonic", lambda: 0.0)

    # Decreasing forces → no stall
    step = [0]
    base_forces = np.full((2, 3), 1.0)

    class ProgressAtoms(_FakeAtoms):
        pass

    atoms = ProgressAtoms(base_forces)
    # Make cached forces decrease each step
    decreasing = [np.full((2, 3), 1.0 / (i + 1)) for i in range(10)]

    class FakeDynWithProgress(_FakeSella):
        def run(self, fmax, steps):
            for i in range(min(steps, 10)):
                atoms.calc.results = {"forces": decreasing[i]}
                for obs in self._observers:
                    obs()

    monkeypatch.setattr("src.ts_optimizer.Sella", FakeDynWithProgress)
    atoms.calc = SimpleNamespace(results={"forces": base_forces})

    ts_opt = TSOptimizer(fmax=0.01, max_steps=10)
    result = ts_opt.find_ts(atoms, atoms.calc, stall_window=8)
    assert result is atoms  # returns the same atoms object


# ---------------------------------------------------------------------------
# DFTRefiner: Sella stall → geomeTRIC fallback
# ---------------------------------------------------------------------------

def test_sella_stall_triggers_geometric_fallback(monkeypatch):
    """When Sella raises TSOptimizationStall, _optimize_with_pyscf_backend
    falls through to geomeTRIC and returns its result."""
    from src.dft_refiner import DFTRefiner

    refiner = DFTRefiner(solvent_name=None)
    geometric_calls = []

    # Make ts_optimizer raise TSOptimizationStall
    class StallOptimizer:
        fmax = 0.05
        def find_ts(self, atoms, calc, **kwargs):
            raise TSOptimizationStall("test stall")
        def is_converged(self, atoms):
            return False

    refiner.ts_optimizer = StallOptimizer()

    def fake_geometric(mf, *, label, max_steps, is_ts_fallback, eff_use_explicit, n_atoms_solute, **kw):
        geometric_calls.append({"label": label, "is_ts_fallback": is_ts_fallback})
        return True, SimpleNamespace(), "2\nopt\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"

    monkeypatch.setattr(refiner, "_optimize_with_geometric_kernel", fake_geometric)
    monkeypatch.setattr(refiner, "_setup_mol", lambda *a, **kw: SimpleNamespace(
        spin=0, atom_coords=lambda: np.zeros((2, 3)),
        atom_symbol=lambda i: "H", natm=2,
        set_geom_=lambda *a, **kw: None,
    ))
    monkeypatch.setattr(refiner, "_build_mf", lambda *a, **kw: SimpleNamespace(
        mol=SimpleNamespace(set_geom_=lambda *a, **kw: None),
        max_cycle=50,
    ))

    conv, mol_opt, opt_xyz = refiner._optimize_with_pyscf_backend(
        "2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        charge=0, spin=0, is_ts=True, max_steps=10,
        eff_use_explicit=False, n_atoms_solute=0, label="test",
    )

    assert conv is True
    assert len(geometric_calls) == 1
    assert geometric_calls[0]["is_ts_fallback"] is True


def test_sella_stall_propagates_to_strategy_escalation(monkeypatch):
    """Full integration: Sella stalls on strategy 1, strategy 2 succeeds via geomeTRIC."""
    from src.dft_refiner import DFTRefiner

    refiner = DFTRefiner(solvent_name=None)
    backend_calls = []

    monkeypatch.setattr(refiner, "_mlp_prerelax", lambda xyz: xyz)
    # MLP returns input unchanged, which would normally trigger the TS-safe xTB
    # pre-relax. Stub it out so this test stays focused on ladder escalation.
    monkeypatch.setattr(
        refiner, "_run_xtb_ts_safe_prerelax", lambda xyz, charge, spin, max_cycles=20: xyz
    )

    # First call: Sella stall + geomeTRIC fails.  Second call: succeeds.
    def fake_opt_backend(
        xyz_content, charge, spin, is_ts, max_steps,
        eff_use_explicit, n_atoms_solute, *,
        label="default", harden_scf=False,
        use_solvent_optimization=False, use_soscf=False,
        smearing_sigma=None, level_shift=0.0, **kwargs,
    ):
        backend_calls.append({"smearing_sigma": smearing_sigma, "harden_scf": harden_scf})
        if len(backend_calls) == 1:
            # Simulates: Sella stalled, geomeTRIC also didn't converge,
            # but produced partial geometry.
            return False, SimpleNamespace(), "2\npartial\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
        return True, SimpleNamespace(), xyz_content

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_build_mf", lambda *a, **kw: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "single_point", lambda *a, **kw: -0.7)
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *a, **kw: ([100.0], -0.9, -0.8))

    result = refiner.optimize_geometry("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", is_ts=True, max_steps=5)

    assert result.converged is True
    assert len(backend_calls) == 2
    # Second strategy used smearing
    assert backend_calls[1]["smearing_sigma"] == 0.005
    assert backend_calls[1]["harden_scf"] is True


# ---------------------------------------------------------------------------
# Watchdog callback: energy plateau detection
# ---------------------------------------------------------------------------

def test_watchdog_callback_detects_energy_plateau():
    """The _watchdog_callback inside _optimize_with_geometric_kernel should
    detect a plateau when PySCF passes envs with 'energy' key (not 'e_tot').
    After _PLATEAU_LOOKBACK (5) cycles with spread < 1e-4 Ha it should raise."""
    from src.dft_refiner import OptimizationTimeoutError

    # Simulate the closure state that the real callback captures
    _plateau_energies: list = []
    _PLATEAU_LOOKBACK = 5
    _PLATEAU_ABS_TOL = 1e-4

    def watchdog(envs):
        e_tot = None
        if isinstance(envs, dict):
            e_tot = envs.get("energy", envs.get("e_tot", None))
        if e_tot is not None:
            _plateau_energies.append(float(e_tot))
            if len(_plateau_energies) >= _PLATEAU_LOOKBACK:
                recent = _plateau_energies[-_PLATEAU_LOOKBACK:]
                spread = max(recent) - min(recent)
                if spread < _PLATEAU_ABS_TOL:
                    raise OptimizationTimeoutError(
                        f"Energy plateau (spread={spread:.2e} Ha over {_PLATEAU_LOOKBACK} cycles)"
                    )

    # First 4 calls — should NOT raise
    for i in range(4):
        watchdog({"energy": -817.409 + i * 1e-5})

    # 5th call with energy within 1e-4 spread → should raise
    with pytest.raises(OptimizationTimeoutError, match="plateau"):
        watchdog({"energy": -817.409 + 2e-5})

    # Verify it would NOT have triggered with 'e_tot' key missing and no 'energy'
    _plateau_energies.clear()
    for i in range(6):
        watchdog({"other_key": -817.0})  # no 'energy' key → no tracking → no raise
