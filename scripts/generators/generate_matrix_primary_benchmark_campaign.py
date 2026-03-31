#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_primary_benchmark_campaign import (  # noqa: E402
    build_matrix_primary_benchmark_campaign_artifact,
    render_matrix_primary_benchmark_campaign_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_matrix_primary_benchmark_campaign_artifact()
    markdown = render_matrix_primary_benchmark_campaign_markdown(payload)

    (output_dir / "matrix_primary_benchmark_campaign.md").write_text(markdown, encoding="utf-8")
    (output_dir / "matrix_primary_benchmark_campaign.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())