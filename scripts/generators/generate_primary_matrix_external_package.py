#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.primary_matrix_external_package import (  # noqa: E402
    build_primary_matrix_external_package_artifact,
    render_primary_matrix_external_package_markdown,
)


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = build_primary_matrix_external_package_artifact()
    (output_dir / "primary_matrix_external_package.md").write_text(
        render_primary_matrix_external_package_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "primary_matrix_external_package.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())