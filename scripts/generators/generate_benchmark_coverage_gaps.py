#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths
from src.benchmark_coverage_gaps import build_benchmark_coverage_gap_rows, render_benchmark_coverage_gap_markdown


def _build_rows() -> list[dict[str, object]]:
    return build_benchmark_coverage_gap_rows()


def _render_markdown(rows: list[dict[str, object]]) -> str:
    return render_benchmark_coverage_gap_markdown(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    args = parser.parse_args()

    rows = _build_rows()
    markdown = _render_markdown(rows)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    markdown_path = output_dir / "benchmark_coverage_gaps.md"
    json_path = output_dir / "benchmark_coverage_gaps.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
