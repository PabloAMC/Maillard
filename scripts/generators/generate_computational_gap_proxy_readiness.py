#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.computational_gap_proxy_readiness import (  # noqa: E402
    build_computational_gap_proxy_readiness_artifact,
    render_computational_gap_proxy_readiness_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Assess whether the computational-gap proxies are ready to stand in "
            "for the missing measurements; writes "
            "results/validation/computational_gap_proxy_readiness.{json,md}."
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

    payload = build_computational_gap_proxy_readiness_artifact()
    markdown = render_computational_gap_proxy_readiness_markdown(payload)

    (output_dir / "computational_gap_proxy_readiness.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (output_dir / "computational_gap_proxy_readiness.md").write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())