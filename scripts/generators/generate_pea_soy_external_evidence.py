#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pea_soy_external_evidence import (  # noqa: E402
    build_pea_soy_external_evidence_artifact,
    render_pea_soy_external_evidence_markdown,
)


def main() -> int:
    output_dir = ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = build_pea_soy_external_evidence_artifact()
    (output_dir / "pea_soy_external_evidence.md").write_text(
        render_pea_soy_external_evidence_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "pea_soy_external_evidence.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())