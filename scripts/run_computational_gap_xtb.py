#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_refinement import build_computational_gap_xtb_job_manifest  # noqa: E402
from src.computational_gap_proxy_readiness import PROXY_LANES  # noqa: E402
from src.xtb_path_quality import assess_xtb_path_quality  # noqa: E402


DEFAULT_OUTPUT_PATH = "results/computational_gap_refinement/computational_gap_xtb_execution.json"


def _load_manifest(path: Path) -> Dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return build_computational_gap_xtb_job_manifest()


def _resolve(relative_path: str | None) -> Path | None:
    if not relative_path:
        return None
    return ROOT / relative_path


def _proxy_job(target_id: str) -> Dict[str, Any] | None:
    for lane in PROXY_LANES:
        if lane.get("target_id") != target_id:
            continue
        geometry_dir = str(lane["geometry_dir"])
        return {
            "target_id": target_id,
            "reaction_key": target_id,
            "status": "proxy_ready",
            "runner_script": f"{geometry_dir}/run_xtb.sh",
            "required_inputs": [
                f"{geometry_dir}/reactant.xyz",
                f"{geometry_dir}/product.xyz",
            ],
            "expected_outputs": [
                f"{geometry_dir}/xtbpath.xyz",
                f"{geometry_dir}/xtbpath_ts.xyz",
            ],
        }
    return None


def _output_path_for_target(output_arg: str, target_id: str) -> Path:
    if target_id != "all" and output_arg == DEFAULT_OUTPUT_PATH:
        return Path(f"results/computational_gap_refinement/{target_id}_xtb_execution.json")
    return Path(output_arg)


def _xyz_atom_count(path: Path | None) -> int | None:
    if path is None or not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                return int(stripped)
            except ValueError:
                return None
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="results/computational_gap_refinement/computational_gap_xtb_job_manifest.json")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--target", default="all")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--timeout", type=int, default=10800)
    args = parser.parse_args()

    manifest = _load_manifest(Path(args.manifest))
    jobs = list(manifest.get("jobs", []))
    if args.target != "all":
        jobs = [job for job in jobs if str(job.get("target_id", "")) == args.target]
        if not jobs:
            proxy_job = _proxy_job(str(args.target))
            if proxy_job is not None:
                jobs = [proxy_job]

    results: List[Dict[str, Any]] = []
    for job in jobs:
        runner_path = _resolve(job.get("runner_script"))
        required_inputs = [_resolve(path) for path in job.get("required_inputs", [])]
        expected_outputs = [_resolve(path) for path in job.get("expected_outputs", [])]
        reactant_atom_count = _xyz_atom_count(required_inputs[0] if len(required_inputs) >= 1 else None)
        product_atom_count = _xyz_atom_count(required_inputs[1] if len(required_inputs) >= 2 else None)
        atom_count_match = (
            reactant_atom_count is not None
            and product_atom_count is not None
            and reactant_atom_count == product_atom_count
        )
        row: Dict[str, Any] = {
            "target_id": job.get("target_id"),
            "reaction_key": job.get("reaction_key"),
            "status": str(job.get("status", "unknown")),
            "execute_mode": bool(args.execute),
            "reactant_atom_count": reactant_atom_count,
            "product_atom_count": product_atom_count,
            "atom_count_match": atom_count_match,
            "expected_outputs_present": [str(path.relative_to(ROOT)) for path in expected_outputs if path and path.exists()],
        }
        if not runner_path or not runner_path.exists() or any(path is None or not path.exists() for path in required_inputs):
            row["status"] = "seed_required"
            results.append(row)
            continue
        if not atom_count_match:
            row["status"] = "blocked_atom_count_mismatch"
            row["blocker_note"] = (
                f"Reactant/product atom counts differ ({reactant_atom_count} vs {product_atom_count}); "
                "xTB path search requires a balanced mapped pair."
            )
            results.append(row)
            continue
        cached_outputs_complete = bool(expected_outputs) and all(
            path is not None and path.exists() for path in expected_outputs
        )
        if args.execute and cached_outputs_complete:
            xtb_quality = assess_xtb_path_quality(runner_path.parent, materialize_missing=True)
            row["status"] = "completed_cached"
            row["cached_outputs_used"] = True
            row["quality_gate_passed"] = bool(xtb_quality["quality_gate_passed"])
            row["failure_markers"] = list(xtb_quality["failure_markers"])
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
            xtb_quality = assess_xtb_path_quality(runner_path.parent, materialize_missing=True)
            if xtb_quality["synthesized_outputs"]:
                expected_outputs = list(dict.fromkeys([*expected_outputs, *xtb_quality["synthesized_outputs"]]))
                row["synthesized_outputs"] = [
                    str(path.relative_to(ROOT)) for path in xtb_quality["synthesized_outputs"]
                ]
            row["quality_gate_passed"] = bool(xtb_quality["quality_gate_passed"])
            row["failure_markers"] = list(xtb_quality["failure_markers"])
            if completed.returncode == 0 and any(path is not None and path.exists() for path in expected_outputs):
                row["status"] = "completed" if xtb_quality["quality_gate_passed"] else "completed_with_quality_warnings"
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
            "completed_count": sum(1 for row in results if row.get("status") in {"completed", "completed_cached"}),
            "completed_cached_count": sum(1 for row in results if row.get("status") == "completed_cached"),
            "completed_with_quality_warnings_count": sum(1 for row in results if row.get("status") == "completed_with_quality_warnings"),
            "seed_required_count": sum(1 for row in results if row.get("status") == "seed_required"),
            "blocked_atom_count_mismatch_count": sum(1 for row in results if row.get("status") == "blocked_atom_count_mismatch"),
            "failed_count": sum(1 for row in results if row.get("status") in {"failed", "timeout"}),
            "execute_mode": bool(args.execute),
        },
        "jobs": results,
    }
    output_path = _output_path_for_target(str(args.output), str(args.target))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(json.dumps(payload["summary"], indent=2))
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())