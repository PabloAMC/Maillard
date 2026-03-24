#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_experiment_intake import build_matrix_experiment_support_delta_artifact
from src.presentation import render_matrix_experiment_support_delta_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    experiment_path = Path(args.experiment)
    payload = build_matrix_experiment_support_delta_artifact(experiment_path)
    markdown = render_matrix_experiment_support_delta_markdown(payload)

    stem = experiment_path.stem
    markdown_path = output_dir / f"{stem}_support_delta.md"
    json_path = output_dir / f"{stem}_support_delta.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())