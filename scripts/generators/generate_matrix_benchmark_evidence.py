#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import build_matrix_benchmark_evidence_audit
from src.presentation import render_matrix_benchmark_evidence_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = build_matrix_benchmark_evidence_audit()
    markdown = render_matrix_benchmark_evidence_markdown(rows)
    payload = []
    for row in rows:
        item = asdict(row)
        item["bench_file"] = str(row.bench_file)
        payload.append(item)

    markdown_path = output_dir / "matrix_benchmark_evidence.md"
    json_path = output_dir / "matrix_benchmark_evidence.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())