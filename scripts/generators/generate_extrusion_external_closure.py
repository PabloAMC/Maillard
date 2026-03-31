#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.extrusion_external_closure import (  # noqa: E402
    build_extrusion_external_closure_artifact,
    render_extrusion_external_closure_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_extrusion_external_closure_artifact()
    markdown = render_extrusion_external_closure_markdown(payload)

    (output_dir / "extrusion_external_closure.md").write_text(markdown, encoding="utf-8")
    (output_dir / "extrusion_external_closure.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())