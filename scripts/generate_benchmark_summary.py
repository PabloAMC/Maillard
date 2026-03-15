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

from src.benchmark_validation import render_benchmark_summary_markdown, summarize_benchmarks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = summarize_benchmarks(target_tag=args.target_tag)
    markdown = render_benchmark_summary_markdown(summaries)
    json_payload = []
    for summary in summaries:
        row = asdict(summary)
        row["bench_file"] = str(summary.bench_file)
        json_payload.append(row)

    markdown_path = output_dir / "benchmark_summary.md"
    json_path = output_dir / "benchmark_summary.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    strict_ready = [summary.benchmark_id for summary in summaries if summary.strict_ready]
    if strict_ready:
        print("Strict-ready benchmarks: " + ", ".join(strict_ready))
    else:
        print("Strict-ready benchmarks: none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())