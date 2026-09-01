#!/usr/bin/env python3
"""Generate matrix-family priority ranking artifacts."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.matrix_family_priority_ranking import (  # noqa: E402
    build_matrix_family_priority_ranking_artifact,
    render_matrix_family_priority_ranking_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the matrix-family priority ranking; writes "
            "results/validation/matrix_family_priority_ranking.{json,md}."
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

    artifact = build_matrix_family_priority_ranking_artifact()
    markdown = render_matrix_family_priority_ranking_markdown(artifact)

    json_path = output_dir / "matrix_family_priority_ranking.json"
    markdown_path = output_dir / "matrix_family_priority_ranking.md"
    json_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())