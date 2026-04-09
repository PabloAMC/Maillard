#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_wave1 import (  # noqa: E402
    apply_wave1_dft_promotions,
    render_wave1_dft_promotion_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execution", default="results/computational_gap_wave1/wave1_dft_execution.json")
    parser.add_argument("--priors", default="data/lit/computational_priors.json")
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    payload = apply_wave1_dft_promotions(
        execution_path=Path(args.execution),
        priors_path=Path(args.priors),
        write_changes=not bool(args.dry_run),
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    markdown_path = output_dir / "wave1_dft_promotion_report.md"
    json_path = output_dir / "wave1_dft_promotion_report.json"
    markdown_path.write_text(render_wave1_dft_promotion_markdown(payload), encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown_path.read_text(encoding="utf-8"))
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())