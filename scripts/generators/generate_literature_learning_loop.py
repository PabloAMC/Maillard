#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.literature_learning_loop import (
    build_literature_learning_loop_payload,
    render_literature_learning_loop_markdown,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_literature_learning_loop_payload(ROOT)
    markdown = render_literature_learning_loop_markdown(payload)

    markdown_path = output_dir / "literature_learning_loop.md"
    json_path = output_dir / "literature_learning_loop.json"
    template_path = output_dir / "literature_runtime_templates.json"

    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    template_path.write_text(json.dumps(payload.get("runtime_templates", []), indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {template_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())