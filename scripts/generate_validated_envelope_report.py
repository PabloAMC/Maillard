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

from src.usability_reports import build_validated_envelope_report, render_validated_envelope_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report = build_validated_envelope_report(target_tag=args.target_tag)
    markdown = render_validated_envelope_markdown(report)

    markdown_path = output_dir / "validated_envelope.md"
    json_path = output_dir / "validated_envelope.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(asdict(report), indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())