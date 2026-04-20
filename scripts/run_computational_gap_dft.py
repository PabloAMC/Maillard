#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import sys
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_refinement import build_computational_gap_dft_job_manifest  # noqa: E402
from src.dft_geometry_preflight import build_dft_geometry_preflight, preflight_candidate_score, repair_steric_clash  # noqa: E402
from src.xtb_path_quality import path_frame_index  # noqa: E402
from src.xyz_common import extract_xyz_last_frame  # noqa: E402

DEFAULT_METHOD_CHAIN = "wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp + ddCOSMO(water)"
FAST_METHOD_CHAIN = "HF/sto-3g exploratory fast mode"

_write_lock = threading.Lock()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _elapsed_seconds(started_at: str | None) -> int | None:
    if not started_at:
        return None
    started = datetime.fromisoformat(started_at)
    return max(0, int((datetime.now(timezone.utc) - started).total_seconds()))


def _build_payload(results: List[Dict[str, Any]], *, execute_mode: bool) -> Dict[str, Any]:
    running_jobs = [row for row in results if row.get("status") == "running"]
    completed_count = sum(1 for row in results if row.get("status") == "completed")
    return {
        "summary": {
            "job_count": len(results),
            "completed_count": completed_count,
            "promotion_ready_count": sum(1 for row in results if row.get("promotion_ready")),
            "invalid_saddle_count": sum(
                1 for row in results if row.get("barrier_quality") == "invalid_saddle"
            ),
            "gsm_attempted_count": sum(
                1 for row in results if isinstance(row.get("gsm_attempt"), dict)
            ),
            "gsm_recovered_true_saddle_count": sum(
                1
                for row in results
                if isinstance(row.get("gsm_attempt"), dict)
                and row["gsm_attempt"].get("recovered_true_saddle")
            ),
            "ready_for_dft_count": sum(1 for row in results if row.get("status") == "ready_for_dft"),
            "blocked_count": sum(
                1
                for row in results
                if row.get("status") in {"blocked_missing_xtb_ts", "seed_required", "blocked_bad_ts_guess"}
            ),
            "failed_count": sum(1 for row in results if row.get("status") == "failed"),
            "interrupted_count": sum(1 for row in results if row.get("status") == "interrupted"),
            "running_count": len(running_jobs),
            "remaining_count": max(0, len(results) - completed_count),
            "active_targets": [row.get("target_id") for row in running_jobs],
            "execute_mode": execute_mode,
            "last_updated_at": _utc_now(),
        },
        "jobs": results,
    }


def _write_output(output_path: Path, results: List[Dict[str, Any]], *, execute_mode: bool) -> Dict[str, Any]:
    payload = _build_payload(results, execute_mode=execute_mode)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with _write_lock:
        output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def _load_manifest(path: Path) -> Dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return build_computational_gap_dft_job_manifest()


def _resolve(relative_path: str | None) -> Path | None:
    if not relative_path:
        return None
    return ROOT / relative_path


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _maybe_repair_ts_guess(reactant_xyz: str, ts_path: Path, initial_report: Dict[str, Any]) -> Dict[str, Any] | None:
    candidate_rows: List[Dict[str, Any]] = []
    for candidate_path in sorted(ts_path.parent.glob("xtbpath_*.xyz"), key=path_frame_index):
        if path_frame_index(candidate_path) < 0 or candidate_path == ts_path:
            continue
        raw_content = candidate_path.read_text(encoding="utf-8")
        last_frame = extract_xyz_last_frame(raw_content)
        candidate_report = build_dft_geometry_preflight(reactant_xyz, last_frame)
        candidate_rows.append(
            {
                "path": candidate_path,
                "report": candidate_report,
                "last_frame": last_frame,
            }
        )

    passing_rows = [row for row in candidate_rows if row["report"].get("quality_gate_passed")]
    if not passing_rows:
        return None

    selected = max(
        passing_rows,
        key=lambda row: (
            *preflight_candidate_score(row["report"]),
            path_frame_index(row["path"]),
        ),
    )
    return {
        "selected_path": selected["path"],
        "selected_report": selected["report"],
        "selected_last_frame": selected["last_frame"],
        "original_report": initial_report,
        "candidate_count": len(candidate_rows),
    }


def _append_phase(row: Dict[str, Any], phase: str, *, detail: Dict[str, Any] | None = None) -> None:
    history = row.setdefault("phase_history", [])
    if history and history[-1].get("phase") == phase:
        return
    history.append({"phase": phase, "at": _utc_now(), "detail": detail or {}})


def _set_phase(row: Dict[str, Any], phase: str, *, detail: Dict[str, Any] | None = None) -> None:
    row["phase"] = phase
    row["phase_updated_at"] = _utc_now()
    if detail:
        row["phase_detail"] = dict(detail)
    _append_phase(row, phase, detail=detail)


def _build_reaction_meta(job: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "target_id": str(job.get("target_id", "unknown")),
        "reaction_key": str(job.get("reaction_key", "unknown")),
        "family": str(job.get("reaction_key", "unknown")),
        "wave": "computational_gap_refinement",
        "forming_bond_atoms": list(job.get("forming_bond_atoms", [])),
        "forming_bond": dict(job.get("forming_bond", {}) or {}),
    }


def _read_product_xyz(job: Dict[str, Any]) -> str | None:
    """Read product XYZ from the job manifest, if available."""
    product_path = _resolve(job.get("product_path"))
    if product_path and product_path.exists():
        return product_path.read_text(encoding="utf-8")
    return None


def _resolve_path_frames_dir(job: Dict[str, Any]) -> str | None:
    """Return the directory containing xtbpath_*.xyz frames, if it exists."""
    ts_path = _resolve(job.get("ts_guess_path"))
    if ts_path and ts_path.exists():
        parent = ts_path.parent
        # Check if parent has xtbpath_1.xyz (sign of path frames)
        if (parent / "xtbpath_1.xyz").exists():
            return str(parent)
    return None


def _apply_steric_repair(xyz: str, row: Dict[str, Any], label: str) -> str:
    """Run steric-clash repair and record metadata in *row* under *label*."""
    repair = repair_steric_clash(xyz)
    if repair.get("repaired"):
        row[f"{label}_preflight_repair"] = {
            "repair_count": len(repair.get("repairs", [])),
            "repair_complete": bool(repair.get("repair_complete", True)),
            "repairs": list(repair.get("repairs", [])),
        }
        return repair["xyz"]
    return xyz


def _start_heartbeat(
    output_path: Path,
    results: List[Dict[str, Any]],
    row: Dict[str, Any],
    *,
    execute_mode: bool,
    interval_seconds: int,
) -> tuple[threading.Event, threading.Thread]:
    stop_event = threading.Event()

    def _heartbeat() -> None:
        while not stop_event.wait(interval_seconds):
            row["updated_at"] = _utc_now()
            row["elapsed_seconds"] = _elapsed_seconds(row.get("started_at"))
            _write_output(output_path, results, execute_mode=execute_mode)

    thread = threading.Thread(target=_heartbeat, name=f"dft-heartbeat-{row.get('target_id', 'job')}", daemon=True)
    thread.start()
    return stop_event, thread


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="results/computational_gap_refinement/computational_gap_dft_job_manifest.json")
    parser.add_argument("--output", default="results/computational_gap_refinement/computational_gap_dft_execution.json")
    parser.add_argument("--target", default="all")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--irc", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--ignore-preflight-failures", action="store_true")
    parser.add_argument("--heartbeat-seconds", type=int, default=30)
    parser.add_argument("--ts-guess-refine", choices=["auto", "always", "off"], default="auto",
                        help="P0 TS-guess refinement policy (default: auto = refine if fmax > threshold)")
    parser.add_argument("--ts-guess-fmax-threshold", type=float, default=1.0,
                        help="fmax threshold (eV/Å) for triggering TS-guess refinement (default: 1.0)")
    parser.add_argument("--enable-gsm-fallback", dest="enable_gsm_fallback",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Run pyGSM DE-GSM as a hard fallback when the TS fails the saddle gate. "
                             "Default: enabled. Use --no-enable-gsm-fallback to disable.")
    parser.add_argument("--gsm-nodes", type=int, default=11,
                        help="Number of pyGSM nodes (default: 11).")
    parser.add_argument("--gsm-max-iters", type=int, default=80,
                        help="Maximum DE-GSM optimization passes (default: 80).")
    parser.add_argument("--gsm-timeout-s", type=float, default=1800.0,
                        help="Wall-clock cap (seconds) for a single pyGSM run (default: 1800).")
    parser.add_argument("--strategy", choices=["full_dft", "robust_sp"], default="full_dft",
                        help="Barrier estimation strategy. 'robust_sp' skips DFT geometry "
                             "optimisation: refines TS at xTB-Sella then runs DFT single-points. "
                             "Much faster (~30 min vs hours) at the cost of less precise geometry.")
    args = parser.parse_args()

    manifest = _load_manifest(Path(args.manifest))
    jobs = list(manifest.get("jobs", []))
    if args.target != "all":
        jobs = [job for job in jobs if str(job.get("target_id", "")) == args.target]
    output_path = Path(args.output)

    refiner = None
    if args.execute:
        from src.dft_refiner import DFTRefiner  # noqa: E402

        refiner = DFTRefiner(solvent_name="water", temp_k=423.15, db_path="results/maillard_results.db")
        if args.fast:
            refiner.opt_basis = "sto-3g"
            refiner.opt_method = "hf"
            refiner.refinement_method = "hf"
            refiner.refinement_basis = "sto-3g"

    results: List[Dict[str, Any]] = []
    _write_output(output_path, results, execute_mode=bool(args.execute))
    for job in jobs:
        reactant_path = _resolve(job.get("reactant_path"))
        ts_path = _resolve(job.get("ts_guess_path"))
        row: Dict[str, Any] = {
            "target_id": job.get("target_id"),
            "reaction_key": job.get("reaction_key"),
            "status": str(job.get("status", "unknown")),
            "execute_mode": bool(args.execute),
            "fast_mode": bool(args.fast),
            "method_chain": FAST_METHOD_CHAIN if args.fast else DEFAULT_METHOD_CHAIN,
            "promotion_ready": False,
        }
        results.append(row)
        heartbeat_stop: threading.Event | None = None
        heartbeat_thread: threading.Thread | None = None
        try:
            if not reactant_path or not reactant_path.exists():
                row["status"] = "seed_required"
            elif not ts_path or not ts_path.exists():
                row["status"] = "blocked_missing_xtb_ts"
            else:
                reactant_xyz = reactant_path.read_text(encoding="utf-8")
                ts_xyz = ts_path.read_text(encoding="utf-8")
                reactant_xyz = _apply_steric_repair(reactant_xyz, row, "reactant")
                ts_xyz = _apply_steric_repair(ts_xyz, row, "ts")
                row["preflight"] = build_dft_geometry_preflight(reactant_xyz, ts_xyz)

                repaired_ts_xyz = ts_xyz
                repaired_ts_path = ts_path
                if not row["preflight"].get("quality_gate_passed"):
                    repair = _maybe_repair_ts_guess(reactant_xyz, ts_path, row["preflight"])
                    if repair is not None:
                        repaired_ts_path = repair["selected_path"]
                        repaired_ts_xyz = repair["selected_last_frame"]
                        row["preflight_repair"] = {
                            "selected_path": _display_path(repaired_ts_path),
                            "candidate_count": int(repair["candidate_count"]),
                            "original_blockers": list(repair["original_report"].get("blockers", [])),
                        }
                        row["preflight"] = dict(repair["selected_report"])
                        row["ts_guess_path_used"] = _display_path(repaired_ts_path)

                if not row["preflight"].get("quality_gate_passed") and not args.ignore_preflight_failures:
                    relaxed_scan_recovery = None
                    if args.execute and refiner is not None and job.get("forming_bond_atoms"):
                        relaxed_scan_recovery = refiner.recover_ts_guess_from_forming_bond(
                            reactant_xyz,
                            repaired_ts_xyz,
                            charge=int(job.get("charge", 0) or 0),
                            spin=int(job.get("spin", 0) or 0),
                            reaction_meta=_build_reaction_meta(job),
                            force_scan=True,
                        )
                        relaxed_scan_xyz = str(relaxed_scan_recovery.get("ts_xyz", repaired_ts_xyz))
                        relaxed_scan_report = build_dft_geometry_preflight(reactant_xyz, relaxed_scan_xyz)
                        row["forming_bond_recovery"] = {
                            "strategy": relaxed_scan_recovery.get("strategy"),
                            "used": bool(relaxed_scan_recovery.get("used")),
                            "imaginary_mode_count": relaxed_scan_recovery.get("imaginary_mode_count"),
                            "initial_imaginary_mode_count": relaxed_scan_recovery.get("initial_imaginary_mode_count"),
                            "scan": relaxed_scan_recovery.get("scan"),
                            "steric_repair": relaxed_scan_recovery.get("steric_repair"),
                        }
                        if relaxed_scan_report.get("quality_gate_passed"):
                            repaired_ts_xyz = relaxed_scan_xyz
                            row["preflight"] = relaxed_scan_report
                            row["ts_guess_path_used"] = "forming_bond_relaxed_scan"

                if not row["preflight"].get("quality_gate_passed") and not args.ignore_preflight_failures:
                    row["status"] = "blocked_bad_ts_guess"
                    row["blocker_note"] = ", ".join(row["preflight"].get("blockers", []))
                    row["restart_recommended"] = True
                elif not args.execute or args.preflight_only:
                    row["status"] = "ready_for_dft"
                    row["command_preview"] = f"python scripts/run_computational_gap_dft.py --target {job.get('target_id')} --execute"
                else:
                    # Gate 3: Pre-execution RMSD sanity check
                    from src.dft_geometry_preflight import _parse_xyz  # noqa: E402
                    from math import sqrt as msqrt
                    _r_atoms, _r_coords = _parse_xyz(reactant_xyz)
                    _t_atoms, _t_coords = _parse_xyz(repaired_ts_xyz)
                    if len(_r_coords) == len(_t_coords) and len(_r_coords) >= 3:
                        _rmsd = msqrt(sum(
                            (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
                            for a, b in zip(_r_coords, _t_coords)
                        ) / len(_r_coords))
                        row["ts_reactant_rmsd_angstrom"] = round(_rmsd, 4)
                        if _rmsd < 0.5:
                            row["status"] = "blocked_ts_is_reactant"
                            row["blocker_note"] = (
                                f"TS guess RMSD vs reactant is only {_rmsd:.4f} Å. "
                                "The TS is essentially the reactant; DFT would produce a meaningless barrier. "
                                "Regenerate the TS guess or use a different xTB path frame."
                            )
                            _write_output(output_path, results, execute_mode=bool(args.execute))
                            continue

                    row["status"] = "running"
                    row["started_at"] = _utc_now()
                    row["updated_at"] = row["started_at"]
                    row["elapsed_seconds"] = 0
                    _set_phase(row, "preflight_passed", detail={"label": "DFT preflight passed"})
                    _write_output(output_path, results, execute_mode=bool(args.execute))
                    heartbeat_stop, heartbeat_thread = _start_heartbeat(
                        output_path,
                        results,
                        row,
                        execute_mode=bool(args.execute),
                        interval_seconds=max(1, int(args.heartbeat_seconds)),
                    )

                    def _progress_callback(phase: str, detail: Dict[str, Any]) -> None:
                        _set_phase(row, phase, detail=detail)
                        row["updated_at"] = _utc_now()
                        row["elapsed_seconds"] = _elapsed_seconds(row.get("started_at"))
                        _write_output(output_path, results, execute_mode=bool(args.execute))

                    # Fragment spins from manifest (for bimolecular reactions)
                    _frag_spins_raw = job.get("fragment_spins")
                    _frag_spins = [int(s) for s in _frag_spins_raw] if _frag_spins_raw else None

                    if args.strategy == "robust_sp":
                        barrier = refiner.calculate_robust_barrier(
                            reactant_xyz,
                            repaired_ts_xyz,
                            charge=int(job.get("charge", 0) or 0),
                            spin=int(job.get("spin", 0) or 0),
                            checkpoint_dir=f"data/geometries/dft_checkpoints/{job.get('target_id', 'unknown')}",
                            progress_callback=_progress_callback,
                            product_xyz=_read_product_xyz(job),
                            reaction_meta=_build_reaction_meta(job),
                            enable_gsm_fallback=bool(args.enable_gsm_fallback),
                            gsm_nodes=int(args.gsm_nodes),
                            gsm_max_iters=int(args.gsm_max_iters),
                            gsm_timeout_s=float(args.gsm_timeout_s),
                        )
                    else:
                        barrier = refiner.calculate_barrier(
                            reactant_xyz,
                            repaired_ts_xyz,
                            charge=int(job.get("charge", 0) or 0),
                            spin=int(job.get("spin", 0) or 0),
                            run_irc=bool(args.irc),
                            reaction_meta=_build_reaction_meta(job),
                            checkpoint_dir=f"data/geometries/dft_checkpoints/{job.get('target_id', 'unknown')}",
                            progress_callback=_progress_callback,
                            product_xyz=_read_product_xyz(job),
                            ts_refine_policy=args.ts_guess_refine,
                            ts_refine_fmax_threshold=args.ts_guess_fmax_threshold,
                            path_frames_dir=_resolve_path_frames_dir(job),
                            enable_gsm_fallback=bool(args.enable_gsm_fallback),
                            gsm_nodes=int(args.gsm_nodes),
                            gsm_max_iters=int(args.gsm_max_iters),
                            gsm_timeout_s=float(args.gsm_timeout_s),
                            fragment_spins=_frag_spins,
                        )
                    row["status"] = "completed"
                    row["barrier_kcal_mol"] = float(barrier)
                    row["barrier_kj_mol"] = round(float(barrier) * 4.184, 2)
                    row["barrier_strategy"] = args.strategy

                    # P5: barrier quality classifier.  Read the per-phase
                    # progress files written by DFTRefiner to decide whether
                    # the result is publication-grade or qualitative-only.
                    quality = "converged"
                    quality_notes: List[str] = []
                    ckpt_dir = ROOT / "data" / "geometries" / "dft_checkpoints" / str(job.get("target_id", "unknown"))
                    for phase in ("reactant", "ts"):
                        p = ckpt_dir / f"{phase}_progress.json"
                        if not p.exists():
                            continue
                        try:
                            ps = json.loads(p.read_text(encoding="utf-8"))
                        except Exception:
                            continue
                        if not ps.get("converged", False):
                            quality = "qualitative_only"
                            status = ps.get("last_strategy_status", "unknown")
                            quality_notes.append(f"{phase}:{status}")
                    # Phase-2 hardening: hard gate on imaginary-frequency
                    # validation. A "barrier" computed from a minimum or a
                    # higher-order saddle is physically meaningless and must
                    # NOT be promoted, regardless of geometry/SCF convergence.
                    ts_validation_path = ckpt_dir / "ts_validation.json"
                    ts_validation: Optional[Dict[str, Any]] = None
                    if ts_validation_path.exists():
                        try:
                            ts_validation = json.loads(ts_validation_path.read_text(encoding="utf-8"))
                        except Exception as exc:
                            quality_notes.append(f"ts_validation:read_error:{exc}")
                    if ts_validation is not None:
                        row["ts_validation"] = ts_validation
                        if not ts_validation.get("is_true_saddle", False):
                            quality = "invalid_saddle"
                            verdict = ts_validation.get("verdict", "unknown")
                            n_sig = ts_validation.get("n_significant_imaginary", 0)
                            n_imag = ts_validation.get("n_imaginary", 0)
                            quality_notes.append(
                                f"ts:{verdict} n_imag={n_imag} n_sig_imag={n_sig}"
                            )
                    else:
                        # Missing validation artifact = cannot prove it's a TS.
                        # Be conservative.
                        quality = "invalid_saddle"
                        quality_notes.append("ts:no_validation_artifact")
                    # Surface pyGSM fallback bookkeeping (when invoked).
                    gsm_path = ckpt_dir / "gsm_attempt.json"
                    if gsm_path.exists():
                        try:
                            row["gsm_attempt"] = json.loads(gsm_path.read_text(encoding="utf-8"))
                            quality_notes.append(
                                "gsm:"
                                + ("recovered" if row["gsm_attempt"].get("recovered_true_saddle") else "no_recovery")
                            )
                        except Exception as exc:
                            quality_notes.append(f"gsm_attempt:read_error:{exc}")
                    row["barrier_quality"] = quality
                    if quality_notes:
                        row["barrier_quality_notes"] = quality_notes
                    # Only promote when fully converged AND TS validation passed
                    # AND not in fast mode.
                    row["promotion_ready"] = (
                        not bool(args.fast) and quality == "converged"
                    )
                    _set_phase(row, "completed", detail={
                        "label": "DFT barrier completed",
                        "barrier_kcal_mol": float(barrier),
                        "barrier_quality": quality,
                    })

                    # P0: TS-guess refinement audit trail
                    _refine_ckpt = ckpt_dir / "ts_xtb_refinement.json"
                    if _refine_ckpt.exists():
                        try:
                            row["ts_guess_refinement"] = json.loads(
                                _refine_ckpt.read_text(encoding="utf-8")
                            )
                        except Exception:
                            pass
        except KeyboardInterrupt:  # pragma: no cover - interactive safeguard
            row["status"] = "interrupted"
            row["error"] = "Interrupted by user"
            _set_phase(row, "interrupted", detail={"label": "Interrupted by user"})
            row["updated_at"] = _utc_now()
            row["elapsed_seconds"] = _elapsed_seconds(row.get("started_at"))
            if heartbeat_stop is not None:
                heartbeat_stop.set()
            if heartbeat_thread is not None:
                heartbeat_thread.join(timeout=1)
            payload = _write_output(output_path, results, execute_mode=bool(args.execute))
            print(json.dumps(payload["summary"], indent=2))
            print(f"Wrote {output_path}")
            raise SystemExit(130)
        except Exception as exc:  # pragma: no cover - exercised only in real QM runs
            from src.dft_refiner import BadTSGuessRejected
            if isinstance(exc, BadTSGuessRejected):
                row["status"] = "blocked_bad_ts_guess"
                row["blocker_note"] = f"P0 TS refinement rejected: {exc}"
                _set_phase(row, "blocked_bad_ts_guess", detail={"label": "TS guess rejected by P0", "error": str(exc)})
            else:
                row["status"] = "failed"
                row["error"] = str(exc)
                _set_phase(row, "failed", detail={"label": "DFT calculation failed", "error": str(exc)})
        finally:
            if heartbeat_stop is not None:
                heartbeat_stop.set()
            if heartbeat_thread is not None:
                heartbeat_thread.join(timeout=1)
            if row.get("started_at"):
                row["updated_at"] = _utc_now()
                row["elapsed_seconds"] = _elapsed_seconds(row.get("started_at"))
        _write_output(output_path, results, execute_mode=bool(args.execute))

    payload = _write_output(output_path, results, execute_mode=bool(args.execute))
    print(json.dumps(payload["summary"], indent=2))
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
