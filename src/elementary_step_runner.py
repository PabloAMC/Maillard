"""Elementary-step runner for multi-step transition-state campaigns.

For reactions whose mechanism is a sequence of elementary steps (e.g.
Schiff base → enol → Amadori), this runner:

1. Reads the per-step dry reactant/product geometries.
2. Optionally adds explicit-water "shuttles" (microsolvation) so each
   step has a 5-membered TS topology instead of the highly strained
   3-membered one.
3. Runs DE-GSM (pyGSM) to obtain a TS guess on the hydrated R/P pair.
4. Validates the guess at xTB level: ``n_imag >= 1`` and the imaginary
   mode must be concentrated on the donor / acceptor / migrating-H atoms
   ("mode-coordinate sanity"); see :func:`src.xtb_backend.validate_ts_mode`.
5. Refines TS + computes the barrier via
   :meth:`src.dft_refiner.DFTRefiner.calculate_robust_barrier` (xTB geometry
   + DFT single-point at r2SCAN-3c).
6. Aggregates per-step barriers into a target-level barrier (default:
   ``max_step``).

The schema is defined in
``data/lit/computational_gap_multistep_targets.json``.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .gsm_backend import GSMRunner
from .microsolvation import HydratedEndpoints, build_hydrated_endpoints
from .xtb_backend import compute_xtb_ts_mode, probe_ts_guess_xtb, validate_ts_mode

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class StepResult:
    step_id: str
    name: str
    status: str  # "completed" | "ts_guess_failed" | "ts_validation_failed" | "barrier_failed"
    barrier_kcal_mol: Optional[float]
    reason: str
    n_imag: Optional[int] = None
    lowest_freq_cm: Optional[float] = None
    mode_concentration: Optional[float] = None
    mode_pass: Optional[bool] = None
    mode_top_atoms: List[List[float]] = field(default_factory=list)
    hydrated_endpoints: Optional[Dict[str, Any]] = None
    gsm: Optional[Dict[str, Any]] = None
    work_dir: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class TargetResult:
    target_id: str
    aggregation: str
    aggregated_barrier_kcal_mol: Optional[float]
    literature_barrier_kcal_mol: Optional[float]
    gap_kcal_mol: Optional[float]
    steps: List[StepResult]
    status: str
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target_id": self.target_id,
            "aggregation": self.aggregation,
            "aggregated_barrier_kcal_mol": self.aggregated_barrier_kcal_mol,
            "literature_barrier_kcal_mol": self.literature_barrier_kcal_mol,
            "gap_kcal_mol": self.gap_kcal_mol,
            "status": self.status,
            "notes": self.notes,
            "steps": [s.to_dict() for s in self.steps],
        }


# ---------------------------------------------------------------------------
# Config loader
# ---------------------------------------------------------------------------

DEFAULT_MULTISTEP_REGISTRY = REPO_ROOT / "data" / "lit" / "computational_gap_multistep_targets.json"


def load_multistep_targets(path: Optional[Path] = None) -> Dict[str, Any]:
    p = Path(path) if path else DEFAULT_MULTISTEP_REGISTRY
    return json.loads(p.read_text(encoding="utf-8"))


def find_target(payload: Dict[str, Any], target_id: str) -> Dict[str, Any]:
    for t in payload.get("targets", []):
        if str(t.get("target_id")) == target_id:
            return t
    raise KeyError(f"target_id={target_id!r} not found in multistep registry")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

class ElementaryStepRunner:
    """Run a multi-step DFT campaign decomposed into elementary steps."""

    def __init__(
        self,
        refiner: Any,  # src.dft_refiner.DFTRefiner; left untyped to avoid circular import
        *,
        gsm_nodes: int = 11,
        gsm_max_iters: int = 80,
        gsm_timeout_s: float = 1800.0,
        ts_mode_threshold: float = 0.3,
    ) -> None:
        self.refiner = refiner
        self.gsm_nodes = int(gsm_nodes)
        self.gsm_max_iters = int(gsm_max_iters)
        self.gsm_timeout_s = float(gsm_timeout_s)
        self.ts_mode_threshold = float(ts_mode_threshold)

    # ──────────────────────────────────────────────────────────────────
    # Single elementary step
    # ──────────────────────────────────────────────────────────────────
    def run_step(self, step_spec: Dict[str, Any], work_dir: Path) -> StepResult:
        step_id = str(step_spec["step_id"])
        name = str(step_spec.get("name", step_id))
        charge = int(step_spec.get("charge", 0))
        spin = int(step_spec.get("spin", 0))
        work_dir.mkdir(parents=True, exist_ok=True)

        dry_r_path = REPO_ROOT / step_spec["dry_reactant_xyz"]
        dry_p_path = REPO_ROOT / step_spec["dry_product_xyz"]
        if not dry_r_path.exists() or not dry_p_path.exists():
            return StepResult(
                step_id=step_id, name=name,
                status="ts_guess_failed",
                barrier_kcal_mol=None,
                reason=f"missing dry endpoint: R={dry_r_path.exists()} P={dry_p_path.exists()}",
                work_dir=str(work_dir),
            )

        # Phase 1: microsolvation → hydrated R/P
        microsolv = step_spec.get("microsolvation") or {}
        topology = str(microsolv.get("topology", "none"))
        if topology == "proton_shuttle":
            shuttles = microsolv.get("shuttles", [])
            shuttle_pairs = [(int(s["donor_atom"]), int(s["acceptor_atom"])) for s in shuttles]
            logger.info("[%s] building hydrated endpoints with %d proton-shuttle waters",
                        step_id, len(shuttle_pairs))
            hydrated = build_hydrated_endpoints(
                dry_r_path, dry_p_path,
                work_dir=work_dir / "hydrated",
                proton_shuttles=shuttle_pairs,
                charge=charge,
            )
            r_path = hydrated.reactant_xyz_path
            p_path = hydrated.product_xyz_path
            hydrated_summary = {
                "n_solute": hydrated.n_solute,
                "n_total": hydrated.n_total,
                "heavy_rmsd_after_align": hydrated.heavy_rmsd_after_align,
                "reactant_min_nonbonded_a": hydrated.clash_report_reactant.min_distance_a,
                "product_min_nonbonded_a": hydrated.clash_report_product.min_distance_a,
                "reactant_clash": hydrated.clash_report_reactant.has_clash,
                "product_clash": hydrated.clash_report_product.has_clash,
            }
        elif topology in ("none", ""):
            r_path, p_path = dry_r_path, dry_p_path
            hydrated_summary = {"topology": "none"}
        else:
            return StepResult(
                step_id=step_id, name=name,
                status="ts_guess_failed",
                barrier_kcal_mol=None,
                reason=f"unsupported microsolvation.topology={topology!r}",
                work_dir=str(work_dir),
            )

        # Phase 2: DE-GSM TS guess
        gsm_dir = work_dir / "gsm"
        runner = GSMRunner(
            work_dir=gsm_dir,
            charge=charge, spin=spin, solvent="water",
            n_nodes=self.gsm_nodes, max_iters=self.gsm_max_iters,
            timeout_s=self.gsm_timeout_s,
        )
        reactant_xyz = r_path.read_text(encoding="utf-8")
        product_xyz = p_path.read_text(encoding="utf-8")
        logger.info("[%s] running DE-GSM (n_nodes=%d, max_iters=%d)",
                    step_id, self.gsm_nodes, self.gsm_max_iters)
        gsm_result = runner.run_de_gsm(reactant_xyz, product_xyz)
        gsm_summary = {
            "converged": gsm_result.converged,
            "reason": gsm_result.reason,
            "elapsed_s": gsm_result.elapsed_s,
            "peak_index": gsm_result.peak_index,
            "n_iters": gsm_result.n_iters,
        }
        if not gsm_result.converged or not gsm_result.ts_xyz:
            return StepResult(
                step_id=step_id, name=name,
                status="ts_guess_failed",
                barrier_kcal_mol=None,
                reason=f"DE-GSM did not produce a usable TS guess: {gsm_result.reason}",
                hydrated_endpoints=hydrated_summary,
                gsm=gsm_summary,
                work_dir=str(work_dir),
            )
        ts_xyz = gsm_result.ts_xyz
        (work_dir / "ts_guess.xyz").write_text(ts_xyz, encoding="utf-8")

        # Phase 3: xTB validation of TS guess (n_imag + mode-coord sanity)
        logger.info("[%s] validating TS guess at xTB level", step_id)
        probe = probe_ts_guess_xtb(ts_xyz, charge=charge, spin=spin)
        n_imag = probe.get("n_imag")
        lowest = probe.get("lowest_freq_cm")
        mode_concentration = None
        mode_pass: Optional[bool] = None
        mode_top: List[List[float]] = []
        expected = step_spec.get("ts_mode_expected_atoms")
        if expected:
            v0 = compute_xtb_ts_mode(ts_xyz, charge=charge, spin=spin)
            if v0 is not None:
                mode_check = validate_ts_mode(
                    np.asarray(v0, dtype=float),
                    expected_atoms=list(expected),
                    threshold=self.ts_mode_threshold,
                )
                mode_concentration = mode_check["concentration"]
                mode_pass = mode_check["pass"]
                mode_top = [[float(i), float(f)] for i, f in mode_check["top_atoms"]]
                np.save(str(work_dir / "ts_v0_xtb.npy"), v0)

        if n_imag is None or n_imag < 1:
            return StepResult(
                step_id=step_id, name=name,
                status="ts_validation_failed",
                barrier_kcal_mol=None,
                reason=f"xTB n_imag={n_imag} (need ≥1); lowest_freq={lowest}",
                n_imag=n_imag, lowest_freq_cm=lowest,
                mode_concentration=mode_concentration, mode_pass=mode_pass,
                mode_top_atoms=mode_top,
                hydrated_endpoints=hydrated_summary, gsm=gsm_summary,
                work_dir=str(work_dir),
            )
        if expected and mode_pass is False:
            return StepResult(
                step_id=step_id, name=name,
                status="ts_validation_failed",
                barrier_kcal_mol=None,
                reason=(
                    f"TS mode does not match expected atoms "
                    f"(concentration={mode_concentration:.3f} < {self.ts_mode_threshold})"
                ),
                n_imag=n_imag, lowest_freq_cm=lowest,
                mode_concentration=mode_concentration, mode_pass=mode_pass,
                mode_top_atoms=mode_top,
                hydrated_endpoints=hydrated_summary, gsm=gsm_summary,
                work_dir=str(work_dir),
            )

        # Phase 4: barrier (xTB-Sella TS refine + DFT SP)
        logger.info("[%s] computing barrier via calculate_robust_barrier", step_id)
        ckpt = work_dir / "barrier_ckpt"
        ckpt.mkdir(parents=True, exist_ok=True)
        # Persist v0 in the checkpoint dir so calculate_robust_barrier can pick it up.
        if (work_dir / "ts_v0_xtb.npy").exists():
            (ckpt / "ts_v0_xtb.npy").write_bytes((work_dir / "ts_v0_xtb.npy").read_bytes())
        try:
            barrier = self.refiner.calculate_robust_barrier(
                reactant_xyz, ts_xyz,
                charge=charge, spin=spin,
                checkpoint_dir=str(ckpt),
                product_xyz=product_xyz,
                reaction_meta={
                    "target_id": step_id,
                    "reaction_key": step_id,
                    "family": "multistep_elementary",
                    "wave": "computational_gap_refinement_multistep",
                },
            )
        except Exception as exc:
            logger.exception("[%s] calculate_robust_barrier failed", step_id)
            return StepResult(
                step_id=step_id, name=name,
                status="barrier_failed",
                barrier_kcal_mol=None,
                reason=f"calculate_robust_barrier raised: {exc}",
                n_imag=n_imag, lowest_freq_cm=lowest,
                mode_concentration=mode_concentration, mode_pass=mode_pass,
                mode_top_atoms=mode_top,
                hydrated_endpoints=hydrated_summary, gsm=gsm_summary,
                work_dir=str(work_dir),
            )

        return StepResult(
            step_id=step_id, name=name,
            status="completed",
            barrier_kcal_mol=float(barrier),
            reason="ok",
            n_imag=n_imag, lowest_freq_cm=lowest,
            mode_concentration=mode_concentration, mode_pass=mode_pass,
            mode_top_atoms=mode_top,
            hydrated_endpoints=hydrated_summary, gsm=gsm_summary,
            work_dir=str(work_dir),
        )

    # ──────────────────────────────────────────────────────────────────
    # Full target (sequence of steps)
    # ──────────────────────────────────────────────────────────────────
    def run_target(
        self,
        target_spec: Dict[str, Any],
        output_dir: Path,
    ) -> TargetResult:
        target_id = str(target_spec["target_id"])
        aggregation = str(target_spec.get("aggregation", "max_step"))
        lit = target_spec.get("literature") or {}
        lit_barrier = lit.get("barrier_kcal_mol")
        steps_specs = list(target_spec.get("elementary_steps", []))
        output_dir.mkdir(parents=True, exist_ok=True)

        results: List[StepResult] = []
        for spec in steps_specs:
            step_id = str(spec["step_id"])
            wd = output_dir / step_id
            logger.info("=" * 72)
            logger.info("[%s] starting elementary step", step_id)
            res = self.run_step(spec, wd)
            results.append(res)
            (wd / "step_result.json").write_text(
                json.dumps(res.to_dict(), indent=2), encoding="utf-8",
            )
            logger.info("[%s] status=%s barrier=%s reason=%s",
                        step_id, res.status, res.barrier_kcal_mol, res.reason)

        completed = [r for r in results if r.status == "completed" and r.barrier_kcal_mol is not None]
        notes = ""
        if not completed:
            agg = None
            status = "failed"
            notes = "no elementary step completed successfully"
        elif len(completed) != len(results):
            agg = None
            status = "partial"
            notes = f"only {len(completed)}/{len(results)} steps completed"
        else:
            barriers = [r.barrier_kcal_mol for r in completed]
            if aggregation == "max_step":
                agg = max(barriers)
            elif aggregation == "sum":
                agg = sum(barriers)
            else:
                agg = max(barriers)
                notes = f"unknown aggregation={aggregation!r}; defaulted to max_step"
            status = "completed"

        gap = None
        if agg is not None and lit_barrier is not None:
            gap = float(agg) - float(lit_barrier)

        result = TargetResult(
            target_id=target_id,
            aggregation=aggregation,
            aggregated_barrier_kcal_mol=agg,
            literature_barrier_kcal_mol=lit_barrier,
            gap_kcal_mol=gap,
            steps=results,
            status=status,
            notes=notes,
        )
        (output_dir / "target_result.json").write_text(
            json.dumps(result.to_dict(), indent=2), encoding="utf-8",
        )
        return result


__all__ = [
    "ElementaryStepRunner",
    "StepResult",
    "TargetResult",
    "load_multistep_targets",
    "find_target",
]
