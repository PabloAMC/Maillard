#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_refinement import build_computational_gap_dft_job_manifest  # noqa: E402


DEFAULT_METHOD_CHAIN = "wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp + ddCOSMO(water)"
FAST_METHOD_CHAIN = "HF/sto-3g exploratory fast mode"


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
            "ready_for_dft_count": sum(1 for row in results if row.get("status") == "ready_for_dft"),
            "blocked_count": sum(1 for row in results if row.get("status") in {"blocked_missing_xtb_ts", "seed_required"}),
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
    parser.add_argument("--heartbeat-seconds", type=int, default=30)
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
        if args.execute and reactant_path and reactant_path.exists() and ts_path and ts_path.exists():
            row["status"] = "running"
            row["started_at"] = _utc_now()
            row["updated_at"] = row["started_at"]
            row["elapsed_seconds"] = 0
            _write_output(output_path, results, execute_mode=bool(args.execute))
            heartbeat_stop, heartbeat_thread = _start_heartbeat(
                output_path,
                results,
                row,
                execute_mode=bool(args.execute),
                interval_seconds=max(1, int(args.heartbeat_seconds)),
            )
        try:
            if not reactant_path or not reactant_path.exists():
                row["status"] = "seed_required"
            elif not ts_path or not ts_path.exists():
                row["status"] = "blocked_missing_xtb_ts"
            elif not args.execute:
                row["status"] = "ready_for_dft"
                row["command_preview"] = f"python scripts/run_computational_gap_dft.py --target {job.get('target_id')} --execute"
            else:
                reactant_xyz = reactant_path.read_text(encoding="utf-8")
                ts_xyz = ts_path.read_text(encoding="utf-8")
                barrier = refiner.calculate_barrier(
                    reactant_xyz,
                    ts_xyz,
                    charge=int(job.get("charge", 0) or 0),
                    spin=int(job.get("spin", 0) or 0),
                    run_irc=bool(args.irc),
                    reaction_meta={
                        "family": str(job.get("reaction_key", "unknown")),
                        "wave": "computational_gap_refinement",
                    },
                )
                row["status"] = "completed"
                row["barrier_kcal_mol"] = float(barrier)
                row["barrier_kj_mol"] = round(float(barrier) * 4.184, 2)
                row["promotion_ready"] = not bool(args.fast)
        except KeyboardInterrupt:  # pragma: no cover - interactive safeguard
            row["status"] = "interrupted"
            row["error"] = "Interrupted by user"
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
            row["status"] = "failed"
            row["error"] = str(exc)
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