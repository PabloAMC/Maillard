#!/usr/bin/env python
"""
results/validation/network_hypotheses.{json,md}: what the cited reaction rules propose from each lane's
reference charge, placed against the engine's reactions (modelled / mechanism known / proposed).
Steps and products only; no rates. Pre-registered in results/validation/network_hypotheses_prereg.md.

    python scripts/generators/generate_network_hypotheses.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io  # noqa: E402
from src.network_hypotheses.report import OUTPUT_JSON, build, render_markdown  # noqa: E402


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(OUTPUT_JSON), help="JSON path; the markdown twin is written beside it")
    args = parser.parse_args(argv)
    payload = build()
    json_path, md_path = artifact_io.write_artifact(payload, Path(args.output), render=render_markdown)
    s = payload["summary"]
    print(f"wrote {json_path} and {md_path}: {s['rules']} rules, {s['steps']} steps {s['steps_by_placement']}, products {s['products_by_kind']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
