#!/usr/bin/env python3
"""Generate matrix target support artifacts for P0.4."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import json

from src import data_paths
from src.benchmark_validation import build_matrix_target_status_artifact
from src.presentation import render_matrix_target_status_markdown


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the P0.4 matrix target-support status artifact; writes "
            "results/validation/matrix_target_status.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(data_paths.VALIDATION_DIR),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    artifact = build_matrix_target_status_artifact()
    markdown = render_matrix_target_status_markdown(artifact)

    json_path = output_dir / "matrix_target_status.json"
    markdown_path = output_dir / "matrix_target_status.md"
    json_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
