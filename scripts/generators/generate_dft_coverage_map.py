#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dft_coverage_map import build_dft_coverage_map_artifact, render_dft_coverage_map_markdown


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the DFT barrier-coverage map; writes "
            "data/lit/dft_coverage_map.json and "
            "results/validation/dft_coverage_map.{json,md}."
        )
    )
    parser.add_argument(
        "--data-dir",
        default=str(ROOT / "data" / "lit"),
        help="directory the lit-side copy is written to",
    )
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "results" / "validation"),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_dft_coverage_map_artifact()
    markdown = render_dft_coverage_map_markdown(payload)

    (data_dir / "dft_coverage_map.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (output_dir / "dft_coverage_map.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (output_dir / "dft_coverage_map.md").write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())