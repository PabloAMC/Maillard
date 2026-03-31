#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.scope_gap_guard import (  # noqa: E402
    build_scope_gap_guard_artifact,
    render_scope_gap_guard_markdown,
)


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = build_scope_gap_guard_artifact()
    (output_dir / "scope_gap_guard.md").write_text(
        render_scope_gap_guard_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "scope_gap_guard.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())