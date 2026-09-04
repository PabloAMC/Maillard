"""data/README.md is generated and must be current (every tracked data file described, no ghosts)."""
from __future__ import annotations

import subprocess
import sys

from src import data_paths


def test_data_readme_is_current():
    completed = subprocess.run(
        [sys.executable, "scripts/generators/build_data_readme.py", "--check"],
        cwd=data_paths.REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
