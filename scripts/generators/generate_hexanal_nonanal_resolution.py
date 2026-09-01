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
from src.hexanal_nonanal_resolution import (  # noqa: E402
    build_hexanal_nonanal_resolution_artifact,
    render_hexanal_nonanal_resolution_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the hexanal / nonanal attribution resolution artifact; "
            "writes results/validation/hexanal_nonanal_resolution.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(data_paths.VALIDATION_DIR),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = build_hexanal_nonanal_resolution_artifact()
    (output_dir / "hexanal_nonanal_resolution.md").write_text(
        render_hexanal_nonanal_resolution_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "hexanal_nonanal_resolution.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())