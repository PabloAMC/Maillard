#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_wave1 import build_wave1_xtb_job_manifest  # noqa: E402


def _load_manifest(path: Path) -> Dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return build_wave1_xtb_job_manifest()


def _resolve(relative_path: str | None) -> Path | None:
    if not relative_path:
        return None
    return ROOT / relative_path


def _path_frame_index(path: Path) -> int:
    match = re.search(r"xtbpath_(\d+)\.xyz$", path.name)
    return int(match.group(1)) if match else -1


def _materialize_xtb_outputs(runner_dir: Path) -> List[Path]:
    frame_paths = sorted(runner_dir.glob("xtbpath_*.xyz"), key=_path_frame_index)
    if not frame_paths:
        return []

    path_bundle = runner_dir / "xtbpath.xyz"
    ts_guess = runner_dir / "xtbpath_ts.xyz"

    if not path_bundle.exists():
        frame_chunks = [frame_path.read_text(encoding="utf-8").strip() for frame_path in frame_paths]
        path_bundle.write_text("\n".join(chunk for chunk in frame_chunks if chunk) + "\n", encoding="utf-8")

    if not ts_guess.exists():
        midpoint_frame = frame_paths[len(frame_paths) // 2]
        ts_guess.write_text(midpoint_frame.read_text(encoding="utf-8"), encoding="utf-8")

    return [path_bundle, ts_guess]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="results/computational_gap_wave1/wave1_xtb_job_manifest.json")
    parser.add_argument("--output", default="results/computational_gap_wave1/wave1_xtb_execution.json")
    parser.add_argument("--target", default="all")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    manifest = _load_manifest(Path(args.manifest))
    jobs = list(manifest.get("jobs", []))
    if args.target != "all":
        jobs = [job for job in jobs if str(job.get("target_id", "")) == args.target]

    results: List[Dict[str, Any]] = []
    for job in jobs:
        runner_path = _resolve(job.get("runner_script"))
        required_inputs = [_resolve(path) for path in job.get("required_inputs", [])]
        expected_outputs = [_resolve(path) for path in job.get("expected_outputs", [])]
        row: Dict[str, Any] = {
            "target_id": job.get("target_id"),
            "reaction_key": job.get("reaction_key"),
            "status": str(job.get("status", "unknown")),
            "execute_mode": bool(args.execute),
            "expected_outputs_present": [str(path.relative_to(ROOT)) for path in expected_outputs if path and path.exists()],
        }
        if not runner_path or not runner_path.exists() or any(path is None or not path.exists() for path in required_inputs):
            row["status"] = "seed_required"
            results.append(row)
            continue
        if not args.execute:
            row["status"] = str(job.get("status", "ready"))
            row["command_preview"] = f"cd {runner_path.parent.relative_to(ROOT)} && bash {runner_path.name}"
            results.append(row)
            continue

        try:
            completed = subprocess.run(
                ["bash", runner_path.name],
                cwd=runner_path.parent,
                capture_output=True,
                text=True,
                timeout=max(1, int(args.timeout)),
                check=False,
            )
            row["returncode"] = int(completed.returncode)
            row["stdout_tail"] = "\n".join(completed.stdout.splitlines()[-10:])
            row["stderr_tail"] = "\n".join(completed.stderr.splitlines()[-10:])
            synthesized_outputs = _materialize_xtb_outputs(runner_path.parent)
            if synthesized_outputs:
                expected_outputs = list(dict.fromkeys([*expected_outputs, *synthesized_outputs]))
            if completed.returncode == 0 and any(path is not None and path.exists() for path in expected_outputs):
                row["status"] = "completed"
            elif completed.returncode == 0:
                row["status"] = "completed_missing_expected_outputs"
            else:
                row["status"] = "failed"
        except subprocess.TimeoutExpired:
            row["status"] = "timeout"
        results.append(row)

    payload = {
        "summary": {
            "job_count": len(results),
            "completed_count": sum(1 for row in results if row.get("status") == "completed"),
            "seed_required_count": sum(1 for row in results if row.get("status") == "seed_required"),
            "failed_count": sum(1 for row in results if row.get("status") in {"failed", "timeout"}),
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