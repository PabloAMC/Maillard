#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.hexanal_nonanal_calibration import (  # noqa: E402
    build_hexanal_nonanal_calibration_artifact,
    render_hexanal_nonanal_calibration_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=data_paths.rel(data_paths.VALIDATION_DIR))
    args = parser.parse_args()

    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_hexanal_nonanal_calibration_artifact()
    markdown = render_hexanal_nonanal_calibration_markdown(payload)

    (output_dir / "hexanal_nonanal_calibration_closure.md").write_text(markdown, encoding="utf-8")
    (output_dir / "hexanal_nonanal_calibration_closure.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())