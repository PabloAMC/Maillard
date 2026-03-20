#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.chemistry_family_scope import (  # noqa: E402
    build_chemistry_family_scope_artifact,
    render_chemistry_family_scope_markdown,
)


OUTPUT_DIR = ROOT / "results" / "validation"


def main() -> int:
    payload = build_chemistry_family_scope_artifact()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    json_path = OUTPUT_DIR / "chemistry_family_scope.json"
    markdown_path = OUTPUT_DIR / "chemistry_family_scope.md"

    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    with open(markdown_path, "w", encoding="utf-8") as handle:
        handle.write(render_chemistry_family_scope_markdown(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())