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
from src.chemistry_family_scope import (  # noqa: E402
    build_chemistry_family_scope_artifact,
    render_chemistry_family_scope_markdown,
)


OUTPUT_DIR = data_paths.VALIDATION_DIR


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the chemistry-family scope artifact; writes "
            "results/validation/chemistry_family_scope.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(OUTPUT_DIR),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    payload = build_chemistry_family_scope_artifact()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "chemistry_family_scope.json"
    markdown_path = output_dir / "chemistry_family_scope.md"

    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    with open(markdown_path, "w", encoding="utf-8") as handle:
        handle.write(render_chemistry_family_scope_markdown(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())