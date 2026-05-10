#!/usr/bin/env python3
"""S20.5 Lane A.1 — generate the external validation candidate inventory.

Wraps :mod:`src.external_validation` as a CLI. Run inside the Docker
container per ``agents.md``::

    ./scripts/docker_maillard.sh run \\
        "python scripts/generators/generate_external_validation_inventory.py"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.external_validation import build_inventory, write_artifact


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="results/validation",
        help="Directory to write the .md and .json artifacts to.",
    )
    parser.add_argument(
        "--basename",
        default="external_validation_candidates",
        help="Filename stem for the artifact pair.",
    )
    args = parser.parse_args()

    inventory = build_inventory()
    paths = write_artifact(
        inventory,
        output_dir=Path(args.output_dir),
        basename=args.basename,
    )
    print(json.dumps({k: str(v) for k, v in paths.items()}, indent=2))
    summary = inventory.to_jsonable()["summary"]
    print(
        f"Total {summary['total_candidates']} | "
        f"executable {summary['executable_candidate']} | "
        f"narrative {summary['narrative_only']} | "
        f"redundant {summary['redundant_with_panel']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
