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

from src.benchmark_validation import audit_all_thermodynamic_gating
from src.presentation import render_thermodynamic_gating_audit_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = audit_all_thermodynamic_gating(target_tag=args.target_tag)
    markdown = render_thermodynamic_gating_audit_markdown(rows)
    payload = []
    for row in rows:
        item = asdict(row)
        item["bench_file"] = str(row.bench_file)
        payload.append(item)

    markdown_path = output_dir / "thermodynamic_gating_audit.md"
    json_path = output_dir / "thermodynamic_gating_audit.json"
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    print(f"Material improvements: {sum(1 for row in rows if row.material_improvement)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())