#!/usr/bin/env python3
"""Generate matrix-family coverage artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.matrix_family_coverage import (  # noqa: E402
    build_matrix_family_coverage_artifact,
    render_matrix_family_coverage_markdown,
)


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)

    artifact = build_matrix_family_coverage_artifact()
    markdown = render_matrix_family_coverage_markdown(artifact)

    json_path = output_dir / "matrix_family_coverage.json"
    markdown_path = output_dir / "matrix_family_coverage.md"
    json_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    markdown_path.write_text(markdown, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())