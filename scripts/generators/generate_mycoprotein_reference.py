#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mycoprotein_reference import (  # noqa: E402
    build_mycoprotein_reference_artifact,
    render_mycoprotein_reference_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the mycoprotein reference artifact; writes "
            "results/validation/mycoprotein_reference.{json,md}."
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
    payload = build_mycoprotein_reference_artifact()
    (output_dir / "mycoprotein_reference.md").write_text(
        render_mycoprotein_reference_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "mycoprotein_reference.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())