#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.inverse_design import InverseDesigner
from src.smirks_engine import ReactionConditions
from src.usability_reports import build_formulation_explainability_payload, render_formulation_explainability_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True)
    parser.add_argument("--target-tag", default="meaty")
    parser.add_argument("--minimize-tag", default="beany")
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    designer = InverseDesigner(target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    formulation = next((row for row in designer.grid if row.get("name") == args.name), None)
    if formulation is None:
        raise SystemExit(f"Unknown formulation name: {args.name}")

    conditions = ReactionConditions(
        pH=formulation.get("ph", 6.0),
        temperature_celsius=formulation.get("temp", 150.0),
        water_activity=formulation.get("aw", 0.95),
        protein_type=formulation.get("protein_type", "free"),
    )
    result = designer.evaluate_single(formulation, conditions)
    payload = build_formulation_explainability_payload(
        formulation,
        result,
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
    )
    markdown = render_formulation_explainability_markdown(payload)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^a-z0-9_]+", "_", formulation.get("name", "formulation").lower()).strip("_")
    markdown_path = output_dir / f"{stem}_explainability.md"
    json_path = output_dir / f"{stem}_explainability.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())