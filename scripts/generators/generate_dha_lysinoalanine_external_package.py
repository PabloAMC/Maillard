#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.dha_lysinoalanine_external_package import (  # noqa: E402
    build_dha_lysinoalanine_external_package_artifact,
    render_dha_lysinoalanine_external_package_markdown,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Assemble the DHA / lysinoalanine external-evidence package; writes "
            "results/validation/dha_lysinoalanine_external_package.{json,md}."
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
    payload = build_dha_lysinoalanine_external_package_artifact()
    (output_dir / "dha_lysinoalanine_external_package.md").write_text(
        render_dha_lysinoalanine_external_package_markdown(payload),
        encoding="utf-8",
    )
    (output_dir / "dha_lysinoalanine_external_package.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())