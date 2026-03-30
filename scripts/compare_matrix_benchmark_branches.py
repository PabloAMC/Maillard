#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    build_matrix_benchmark_deltas,
    build_matrix_benchmark_evidence_audit,
    compare_matrix_benchmark_delta_sets,
    render_matrix_branch_deltas_markdown,
)


def _git_output(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True)


def _benchmark_paths_for_ref(ref: str) -> list[str]:
    rows = _git_output("ls-tree", "-r", "--name-only", ref, "data/benchmarks").splitlines()
    return [row for row in rows if row.endswith(".json")]


def _materialize_ref_benchmarks(ref: str, target_root: Path) -> list[Path]:
    paths: list[Path] = []
    for relative_path in _benchmark_paths_for_ref(ref):
        blob = _git_output("show", f"{ref}:{relative_path}")
        destination = target_root / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(blob, encoding="utf-8")
        paths.append(destination)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-ref", default="main")
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    current_rows = build_matrix_benchmark_deltas(target_tag=args.target_tag)
    current_evidence = build_matrix_benchmark_evidence_audit()

    with tempfile.TemporaryDirectory(prefix="matrix-branch-compare-") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        baseline_files = _materialize_ref_benchmarks(args.base_ref, tmp_dir)
        baseline_rows = build_matrix_benchmark_deltas(baseline_files, target_tag=args.target_tag)
        baseline_evidence = build_matrix_benchmark_evidence_audit(baseline_files)

    rows = compare_matrix_benchmark_delta_sets(
        current_rows,
        baseline_rows,
        current_evidence=current_evidence,
        baseline_evidence=baseline_evidence,
    )
    markdown = render_matrix_branch_deltas_markdown(rows, base_ref=args.base_ref)

    payload = {
        "base_ref": args.base_ref,
        "current_row_count": len(current_rows),
        "baseline_row_count": len(baseline_rows),
        "changes": [asdict(row) for row in rows],
    }

    markdown_path = output_dir / "matrix_branch_delta_report.md"
    json_path = output_dir / "matrix_branch_delta_report.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())