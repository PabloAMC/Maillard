from __future__ import annotations

import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src.artifact_io import repo_root


SKIP_LINE_PATTERN = re.compile(r"^SKIPPED \[(?P<count>\d+)\] (?P<path>.+?):(?P<line>\d+): (?P<reason>.+)$")


def _classify_dependency(reason: str, path: str) -> str:
    normalized = reason.lower()
    if "not installed" in normalized or "binary not found" in normalized or "not in path" in normalized:
        return "missing_optional_backend"
    if "dataset not mounted" in normalized or "phase 3.3" in normalized:
        return "missing_external_dataset"
    if "fixtures not mounted" in normalized or "fixture triplets not mounted" in normalized:
        return "missing_external_dataset"
    if "phase 3.5" in normalized:
        return "long_running_campaign"
    if "implementation pending" in normalized or "implementation missing" in normalized:
        return "not_implemented_module"
    if "screening results not found" in normalized or "benchmark file not found" in normalized:
        return "missing_external_dataset"
    return "unknown"


def _classify_lane(path: str, reason: str) -> str:
    normalized_path = path.replace("\\", "/")
    normalized_reason = reason.lower()
    if normalized_path.endswith("test_quasi_harmonic_correction.py"):
        return "deterministic_helper_lane"
    if "/tests/qm/" in normalized_path and "mace" in normalized_reason:
        return "optional_mlp_acceleration_lane"
    if "/tests/qm/" in normalized_path and any(token in normalized_reason for token in ["xtb", "crest", "pyscf", "sella", "solvation"]):
        return "optional_dft_authority_lane"
    if normalized_path.endswith("test_barrier_benchmarks.py") or normalized_path.endswith("test_irc_validation.py"):
        return "optional_dft_authority_lane"
    return "default_deterministic_lane"


def _owner_for_class(dependency_class: str) -> str:
    return {
        "missing_optional_backend": "environment_backend_owner",
        "missing_external_dataset": "benchmark_dataset_owner",
        "not_implemented_module": "qm_infrastructure_owner",
        "long_running_campaign": "offline_research_campaign_owner",
    }.get(dependency_class, "unassigned")


def _unblock_criteria(dependency_class: str, lane: str) -> str:
    if dependency_class == "missing_optional_backend":
        return "install the backend in the selected environment and rerun the capability probe"
    if dependency_class == "missing_external_dataset":
        return "mount the declared dataset or generate the required benchmark fixture before rerunning the lane"
    if dependency_class == "not_implemented_module":
        return "implement the exported helper or authority-lane entry point and convert the placeholder to executable assertions"
    if dependency_class == "long_running_campaign":
        return "complete the offline campaign and write the resulting artifact back into the benchmark lane"
    if lane == "deterministic_helper_lane":
        return "promote the helper to executable deterministic coverage"
    return "declare the missing contract explicitly and rerun the lane"


def run_skip_scan(test_targets: Optional[Iterable[str]] = None) -> Dict[str, Any]:
    root = repo_root()
    targets = list(test_targets or ["tests/benchmarks", "tests/qm"])
    command = [sys.executable, "-m", "pytest", *targets, "-rs", "-q"]
    result = subprocess.run(
        command,
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    output = (result.stdout or "") + ("\n" + result.stderr if result.stderr else "")
    if result.returncode not in {0, 5}:
        raise RuntimeError(f"pytest skip scan failed with code {result.returncode}\n{output}")
    return {
        "command": " ".join(command),
        "returncode": result.returncode,
        "output": output,
    }


def build_skip_policy_registry(test_targets: Optional[Iterable[str]] = None) -> Dict[str, Any]:
    scan = run_skip_scan(test_targets)
    rows: List[Dict[str, Any]] = []
    lane_counts: Dict[str, int] = defaultdict(int)
    class_counts: Dict[str, int] = defaultdict(int)

    for line in scan["output"].splitlines():
        match = SKIP_LINE_PATTERN.match(line.strip())
        if not match:
            continue
        count = int(match.group("count"))
        path = match.group("path")
        reason = match.group("reason")
        dependency_class = _classify_dependency(reason, path)
        lane = _classify_lane(path, reason)
        rows.append(
            {
                "path": path,
                "line": int(match.group("line")),
                "reason": reason,
                "count": count,
                "dependency_class": dependency_class,
                "lane": lane,
                "owner": _owner_for_class(dependency_class),
                "unblock_criteria": _unblock_criteria(dependency_class, lane),
            }
        )
        lane_counts[lane] += count
        class_counts[dependency_class] += count

    rows.sort(key=lambda row: (-int(row["count"]), str(row["path"]), str(row["reason"])))
    return {
        "summary": {
            "scan_command": scan["command"],
            "skip_count": sum(int(row["count"]) for row in rows),
            "skip_cluster_count": len(rows),
            "dependency_class_counts": dict(sorted(class_counts.items())),
            "lane_counts": dict(sorted(lane_counts.items())),
            "quasi_harmonic_helper_active": True,
            "default_ci_contract": "default lanes stay deterministic; optional QM and MLP lanes remain explicitly gated",
        },
        "rows": rows,
    }


def render_skip_policy_registry_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# S9 Skip Registry",
        "",
        f"Scan command: {summary.get('scan_command', 'unknown')}",
        f"Total skips: {int(summary.get('skip_count', 0))}",
        f"Skip clusters: {int(summary.get('skip_cluster_count', 0))}",
        f"Quasi-harmonic helper active: {summary.get('quasi_harmonic_helper_active', False)}",
        f"Default CI contract: {summary.get('default_ci_contract', 'unknown')}",
        "",
        "## Skip Clusters",
        "",
        "| File | Line | Count | Class | Lane | Owner | Reason | Unblock Criteria |",
        "| --- | ---: | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("rows", []):
        lines.append(
            f"| {row.get('path', 'unknown')} | {int(row.get('line', 0))} | {int(row.get('count', 0))} | "
            f"{row.get('dependency_class', 'unknown')} | {row.get('lane', 'unknown')} | {row.get('owner', 'unknown')} | "
            f"{row.get('reason', 'unknown')} | {row.get('unblock_criteria', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Lane Policy",
        "",
        f"Dependency classes: {', '.join(f'{key}={value}' for key, value in summary.get('dependency_class_counts', {}).items()) or 'none'}",
        f"Lane counts: {', '.join(f'{key}={value}' for key, value in summary.get('lane_counts', {}).items()) or 'none'}",
    ])
    return "\n".join(lines) + "\n"