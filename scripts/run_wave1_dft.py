#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_wave1 import build_wave1_dft_job_manifest  # noqa: E402


DEFAULT_METHOD_CHAIN = "wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp + ddCOSMO(water)"
FAST_METHOD_CHAIN = "HF/sto-3g exploratory fast mode"


def _load_manifest(path: Path) -> Dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return build_wave1_dft_job_manifest()


def _resolve(relative_path: str | None) -> Path | None:
    if not relative_path:
        return None
    return ROOT / relative_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="results/computational_gap_wave1/wave1_dft_job_manifest.json")
    parser.add_argument("--output", default="results/computational_gap_wave1/wave1_dft_execution.json")
    parser.add_argument("--target", default="all")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--irc", action="store_true")
    args = parser.parse_args()

    manifest = _load_manifest(Path(args.manifest))
    jobs = list(manifest.get("jobs", []))
    if args.target != "all":
        jobs = [job for job in jobs if str(job.get("target_id", "")) == args.target]

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
        if not reactant_path or not reactant_path.exists():
            row["status"] = "seed_required"
            results.append(row)
            continue
        if not ts_path or not ts_path.exists():
            row["status"] = "blocked_missing_xtb_ts"
            results.append(row)
            continue
        if not args.execute:
            row["status"] = "ready_for_dft"
            row["command_preview"] = f"python scripts/run_wave1_dft.py --target {job.get('target_id')} --execute"
            results.append(row)
            continue

        try:
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
                    "wave": "computational_gap_wave1",
                },
            )
            row["status"] = "completed"
            row["barrier_kcal_mol"] = float(barrier)
            row["barrier_kj_mol"] = round(float(barrier) * 4.184, 2)
            row["promotion_ready"] = not bool(args.fast)
        except Exception as exc:  # pragma: no cover - exercised only in real QM runs
            row["status"] = "failed"
            row["error"] = str(exc)
        results.append(row)

    payload = {
        "summary": {
            "job_count": len(results),
            "completed_count": sum(1 for row in results if row.get("status") == "completed"),
            "promotion_ready_count": sum(1 for row in results if row.get("promotion_ready")),
            "ready_for_dft_count": sum(1 for row in results if row.get("status") == "ready_for_dft"),
            "blocked_count": sum(1 for row in results if row.get("status") in {"blocked_missing_xtb_ts", "seed_required"}),
            "failed_count": sum(1 for row in results if row.get("status") == "failed"),
            "execute_mode": bool(args.execute),
        },
        "jobs": results,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(json.dumps(payload["summary"], indent=2))
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())