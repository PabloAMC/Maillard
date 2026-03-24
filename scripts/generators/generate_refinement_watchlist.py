#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.refinement_watchlist import build_refinement_watchlist, render_refinement_watchlist_markdown


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/validation")
    parser.add_argument("--target-tag", default="meaty")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = build_refinement_watchlist(target_tag=args.target_tag)
    markdown = render_refinement_watchlist_markdown(payload)

    markdown_path = output_dir / "refinement_watchlist.md"
    json_path = output_dir / "refinement_watchlist.json"
    jobs_path = output_dir / "offline_dft_jobs.json"

    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    jobs_path.write_text(json.dumps(payload.get("offline_jobs", []), indent=2), encoding="utf-8")

    print(markdown)
    print(f"Wrote {markdown_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {jobs_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())