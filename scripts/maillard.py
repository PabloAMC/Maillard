#!/usr/bin/env python
"""`python scripts/maillard.py <verb>`: the container's entry to the front door, src/cli.py.
Outside the container, `pip install -e .` gives the same thing as the `maillard` command."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
