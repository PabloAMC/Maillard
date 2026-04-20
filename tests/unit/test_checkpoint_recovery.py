"""Tests for the 3-tier checkpoint / recovery architecture in calculate_barrier."""

import json
import pathlib
from types import SimpleNamespace

import pytest

from src.dft_refiner import DFTRefiner, DFTResult


# ───────── helpers ──────────────────────────────────────────────────────────

MINIMAL_XYZ = "1\n\nH 0.0 0.0 0.0\n"

FAKE_RESULT = DFTResult(
    method="wB97M-V/def2-tzvp",
    energy_hartree=-76.123456,
    gibbs_free_energy_hartree=-76.0,
    quasi_harmonic_gibbs_hartree=-76.01,
    optimized_xyz=MINIMAL_XYZ,
    converged=True,
    frequencies_cm1=[100.0, 200.0, -350.0],
)


def _make_refiner_with_fakes(monkeypatch):
    """Return a DFTRefiner with all heavy-compute methods stubbed out."""
    refiner = DFTRefiner(solvent_name=None)
    refiner._prerelax_mlp = None

    opt_calls: list = []

    def fake_optimize(xyz, charge=0, spin=0, is_ts=False, max_steps=100,
                      label="default", use_explicit_solvent=None, n_water=None,
                      progress_callback=None, ts_refine_ctx=None):
        opt_calls.append({"is_ts": is_ts, "max_steps": max_steps, "xyz": xyz})
        return DFTResult(
            method="wB97M-V/def2-tzvp",
            energy_hartree=-76.0,
            gibbs_free_energy_hartree=-75.9,
            quasi_harmonic_gibbs_hartree=-75.91,
            optimized_xyz=xyz,
            converged=True,
            frequencies_cm1=[-350.0] if is_ts else [100.0],
        )

    monkeypatch.setattr(refiner, "optimize_geometry", fake_optimize)
    return refiner, opt_calls


# ───────── DFTResult serialization ──────────────────────────────────────────

def test_dft_result_round_trip():
    """to_dict → from_dict preserves all fields."""
    reconstructed = DFTResult.from_dict(FAKE_RESULT.to_dict())
    assert reconstructed.method == FAKE_RESULT.method
    assert reconstructed.energy_hartree == FAKE_RESULT.energy_hartree
    assert reconstructed.gibbs_free_energy_hartree == FAKE_RESULT.gibbs_free_energy_hartree
    assert reconstructed.quasi_harmonic_gibbs_hartree == FAKE_RESULT.quasi_harmonic_gibbs_hartree
    assert reconstructed.optimized_xyz == FAKE_RESULT.optimized_xyz
    assert reconstructed.converged == FAKE_RESULT.converged
    assert reconstructed.frequencies_cm1 == FAKE_RESULT.frequencies_cm1


def test_dft_result_from_dict_missing_optional():
    """from_dict tolerates absent optional fields."""
    minimal = {"method": "hf/sto-3g", "energy_hartree": -1.0}
    r = DFTResult.from_dict(minimal)
    assert r.gibbs_free_energy_hartree is None
    assert r.frequencies_cm1 is None
    assert r.converged is False


def test_dft_result_json_serializable():
    """to_dict output can survive json.dumps → json.loads."""
    raw = json.dumps(FAKE_RESULT.to_dict())
    loaded = json.loads(raw)
    r = DFTResult.from_dict(loaded)
    assert r.energy_hartree == FAKE_RESULT.energy_hartree


# ───────── Tier 1: full result checkpoint ───────────────────────────────────

def test_tier1_reactant_result_skips_all_compute(monkeypatch, tmp_path):
    """If reactant_result.json exists, optimize_geometry is never called for reactant."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    # Seed tier-1 checkpoint for reactant only
    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    (ckpt / "reactant_result.json").write_text(json.dumps(FAKE_RESULT.to_dict()))

    barrier = refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    # Only TS should have been computed (1 call), not reactant
    assert len(opt_calls) == 1
    assert opt_calls[0]["is_ts"] is True
    assert isinstance(barrier, float)


def test_tier1_both_results_skip_all_compute(monkeypatch, tmp_path):
    """If both result JSONs exist, zero optimize_geometry calls are made."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    (ckpt / "reactant_result.json").write_text(json.dumps(FAKE_RESULT.to_dict()))

    ts_result = DFTResult(
        method="wB97M-V/def2-tzvp",
        energy_hartree=-76.0,
        gibbs_free_energy_hartree=-75.8,
        quasi_harmonic_gibbs_hartree=-75.81,
        optimized_xyz=MINIMAL_XYZ,
        converged=True,
        frequencies_cm1=[-350.0],
    )
    (ckpt / "ts_result.json").write_text(json.dumps(ts_result.to_dict()))

    barrier = refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    assert len(opt_calls) == 0
    assert isinstance(barrier, float)


# ───────── Tier 2: geometry-only checkpoint ─────────────────────────────────

def test_tier2_geometry_triggers_max_steps_zero(monkeypatch, tmp_path):
    """If reactant_optimized.xyz exists but no JSON, optimize_geometry is called with max_steps=0."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    (ckpt / "reactant_optimized.xyz").write_text(MINIMAL_XYZ)

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    # Two calls: reactant (max_steps=0) + TS (full)
    assert len(opt_calls) == 2
    reactant_call = [c for c in opt_calls if not c["is_ts"]][0]
    assert reactant_call["max_steps"] == 0


# ───────── Tier 3: no checkpoint → full pipeline ───────────────────────────

def test_tier3_no_checkpoint_runs_full_pipeline(monkeypatch, tmp_path):
    """With an empty checkpoint dir, both phases run fully."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    assert len(opt_calls) == 2
    reactant_call = [c for c in opt_calls if not c["is_ts"]][0]
    ts_call = [c for c in opt_calls if c["is_ts"]][0]
    assert reactant_call["max_steps"] == 100
    assert ts_call["max_steps"] == 100


# ───────── Checkpoint creation: both XYZ + JSON saved ──────────────────────

def test_fresh_run_creates_both_checkpoints(monkeypatch, tmp_path):
    """A fresh run should save both XYZ and JSON for each phase."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    assert (ckpt / "reactant_optimized.xyz").exists()
    assert (ckpt / "reactant_result.json").exists()
    assert (ckpt / "ts_optimized.xyz").exists()
    assert (ckpt / "ts_result.json").exists()

    # JSON can be deserialized back
    r = DFTResult.from_dict(json.loads((ckpt / "reactant_result.json").read_text()))
    assert r.energy_hartree == -76.0


# ───────── Corrupt tier-1 falls through to tier-2 ─────────────────────────

def test_corrupt_json_falls_through_to_tier2(monkeypatch, tmp_path):
    """If the JSON checkpoint is corrupt, the system falls through to geometry checkpoint."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    (ckpt / "reactant_result.json").write_text("NOT JSON")
    (ckpt / "reactant_optimized.xyz").write_text(MINIMAL_XYZ)

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    reactant_call = [c for c in opt_calls if not c["is_ts"]][0]
    assert reactant_call["max_steps"] == 0  # tier-2: geometry-only resume


def test_partial_progress_resume_uses_inflight_geometry(monkeypatch, tmp_path):
    """If a partial inflight checkpoint exists, resume from that geometry instead of seed XYZ."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    inflight_xyz = "1\n\nH 1.0 0.0 0.0\n"
    (ckpt / "reactant_inflight.xyz").write_text(inflight_xyz)
    (ckpt / "reactant_progress.json").write_text(json.dumps({
        "phase_prefix": "reactant",
        "mlp_completed": True,
        "xtb_completed": True,
        "next_strategy_index": 2,
        "geometry_ready_for_postprocessing": False,
    }))

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    reactant_call = [c for c in opt_calls if not c["is_ts"]][0]
    assert reactant_call["max_steps"] == 100
    assert reactant_call["xyz"] == inflight_xyz


def test_geometry_ready_progress_skips_to_postprocessing(monkeypatch, tmp_path):
    """If partial progress says geometry is ready, resume with max_steps=0."""
    refiner, opt_calls = _make_refiner_with_fakes(monkeypatch)

    ckpt = tmp_path / "ckpt"
    ckpt.mkdir()
    inflight_xyz = "1\n\nH 2.0 0.0 0.0\n"
    (ckpt / "reactant_inflight.xyz").write_text(inflight_xyz)
    (ckpt / "reactant_progress.json").write_text(json.dumps({
        "phase_prefix": "reactant",
        "geometry_ready_for_postprocessing": True,
        "current_stage": "hessian",
    }))

    refiner.calculate_barrier(
        MINIMAL_XYZ, MINIMAL_XYZ,
        checkpoint_dir=str(ckpt),
        reaction_meta={"target_id": "test"},
    )

    reactant_call = [c for c in opt_calls if not c["is_ts"]][0]
    assert reactant_call["max_steps"] == 0
    assert reactant_call["xyz"] == inflight_xyz


# ───────── Non-TS strategy exhaustion → fallback, not raise ────────────────

def test_non_ts_strategy_exhaustion_uses_fallback(monkeypatch):
    """Non-TS should fall back to best-available geometry instead of raising."""
    refiner = DFTRefiner(solvent_name=None)
    refiner._prerelax_mlp = None

    call_count = 0
    def fake_opt_backend(xyz_content, charge, spin, is_ts, max_steps,
                         eff_use_explicit, n_atoms_solute, *, label="default",
                         harden_scf=False, use_solvent_optimization=False,
                         use_soscf=False, smearing_sigma=None, level_shift=0.0, **kwargs):
        nonlocal call_count
        call_count += 1
        raise RuntimeError("Nuclear gradients of <scanner> not converged")

    def fake_fallback(xyz, *, charge, spin, basis, label):
        fake_mol = SimpleNamespace(natm=1)
        return False, fake_mol, xyz

    monkeypatch.setattr(refiner, "_optimize_with_pyscf_backend", fake_opt_backend)
    monkeypatch.setattr(refiner, "_best_available_geometry_fallback", fake_fallback)
    monkeypatch.setattr(refiner, "_build_mf", lambda *a, **kw: SimpleNamespace(e_tot=-1.0))
    monkeypatch.setattr(refiner, "_run_hessian_and_thermo", lambda *a, **kw: ([100.0], -0.9, -0.8))
    monkeypatch.setattr(refiner, "single_point", lambda *a, **kw: -0.7)

    result = refiner.optimize_geometry(MINIMAL_XYZ, is_ts=False, max_steps=5)

    # All 4 non-TS strategies should have been tried
    assert call_count == 4
    # Should NOT have raised — got a DFTResult back
    assert isinstance(result, DFTResult)
