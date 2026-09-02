#!/usr/bin/env python
"""
THE CORE ON THE DIRECTIONAL CLAIMS PANEL (step B7).

Writes ``results/validation/core_directional_scores.{json,md}``: every claim of
``docs/validation/directional_claims_panel.yml`` scored on the kinetic core through
the same front door a user calls. See ``src/kinetic_core/directional.py`` for the rules.

Usage:
    python scripts/generators/generate_core_directional_scores.py [--output PATH.json]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(data_paths.CORE_DIRECTIONAL_SCORES))
    parser.add_argument("--panel", default=str(data_paths.DIRECTIONAL_CLAIMS_PANEL))
    args = parser.parse_args(argv)

    from src.kinetic_core.directional import score_panel, write_artifact

    started = time.perf_counter()
    payload = score_panel(args.panel)
    payload["summary"]["wall_seconds"] = round(time.perf_counter() - started, 1)
    json_path, md_path = write_artifact(payload, Path(args.output))
    s = payload["summary"]
    ind, allc = s["independent"]["total"], s["all_claims"]["total"]
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(
        f"headline (independent) {ind['agree']}/{ind['evaluable']} ({ind['not_evaluable']} n.e.) | "
        f"all claims {allc['agree']}/{allc['evaluable']} ({allc['not_evaluable']} n.e.) | "
        f"wall={s['wall_seconds']}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
