#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_family_next_action import (  # noqa: E402
    build_matrix_family_next_action_artifact,
    render_matrix_family_next_action_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the matrix-family next-action artifact; writes "
            "results/validation/matrix_family_next_action.{json,md}."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(ROOT / "results" / "validation"),
        help="directory the artifacts are written to",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_matrix_family_next_action_artifact()
    markdown = render_matrix_family_next_action_markdown(payload)

    (output_dir / "matrix_family_next_action.md").write_text(markdown, encoding="utf-8")
    (output_dir / "matrix_family_next_action.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())