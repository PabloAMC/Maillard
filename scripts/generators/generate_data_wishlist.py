#!/usr/bin/env python3
"""Write results/validation/data_wishlist.{json,md}: what to measure next and what it unlocks.

Reads the tracked scorecard, slice profile, Laplace covariance, directional scorecard and
value-of-information ranking (src/data_wishlist.py). Cheap (< 1 s); the artifact-freshness
gate regenerates it and compares.

    python scripts/generators/generate_data_wishlist.py            # write the pair
    python scripts/generators/generate_data_wishlist.py --output /tmp/w.json
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_paths  # noqa: E402
from src.data_wishlist import OUTPUT_JSON, build, render_markdown  # noqa: E402


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default=str(OUTPUT_JSON), help="JSON path; the markdown twin is written beside it")
    args = parser.parse_args(argv)
    payload = build()
    json_path, md_path = artifact_io.write_artifact(payload, Path(args.output), render=render_markdown)
    s = payload["summary"]
    print(
        f"wrote {data_paths.rel(json_path)} and {data_paths.rel(md_path)}: "
        f"{s['unidentified_coordinates']} unidentified coordinates | {s['not_evaluable_rows']} rows not evaluable | "
        f"{s['refused_row_count']} rows refused | {s['axes_below_trust']} axes below trust"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
