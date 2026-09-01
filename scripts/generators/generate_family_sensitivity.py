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
from src.family_sensitivity import build_family_sensitivity_artifact, render_family_sensitivity_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    parser.add_argument("--target-tag", default="meaty")
    parser.add_argument("--delta-kcal", type=float, default=3.0)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_family_sensitivity_artifact(target_tag=args.target_tag, delta_kcal=args.delta_kcal)
    markdown = render_family_sensitivity_markdown(payload)

    markdown_path = output_dir / "family_sensitivity.md"
    json_path = output_dir / "family_sensitivity.json"

    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())