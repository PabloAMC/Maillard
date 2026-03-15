#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import render_benchmark_targets_markdown, snapshot_all_benchmark_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    snapshots = snapshot_all_benchmark_targets(target_tag=args.target_tag)
    markdown = render_benchmark_targets_markdown(snapshots)
    payload = []
    for snapshot in snapshots:
        row = asdict(snapshot)
        row["bench_file"] = str(snapshot.bench_file)
        payload.append(row)

    markdown_path = output_dir / "benchmark_targets.md"
    json_path = output_dir / "benchmark_targets.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    print(f"Low-headspace rows: {sum(1 for row in snapshots if row.headspace_class == 'low_headspace')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())