#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dft_coverage_map import build_dft_coverage_map_artifact, render_dft_coverage_map_markdown


def main() -> int:
    data_dir = ROOT / "data" / "lit"
    output_dir = ROOT / "results" / "validation"
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