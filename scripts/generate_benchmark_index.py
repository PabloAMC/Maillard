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

from src.benchmark_validation import build_benchmark_index, render_benchmark_index_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = build_benchmark_index(target_tag=args.target_tag)
    markdown = render_benchmark_index_markdown(entries)
    payload = []
    for entry in entries:
        row = asdict(entry)
        row["bench_file"] = str(entry.bench_file)
        payload.append(row)

    markdown_path = output_dir / "benchmark_index.md"
    json_path = output_dir / "benchmark_index.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())