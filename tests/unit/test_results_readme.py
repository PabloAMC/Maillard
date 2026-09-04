"""results/README.md is generated and current (2026-09-03, the data/ treatment for results/)."""
from __future__ import annotations

import subprocess
import sys

from src import data_paths


def test_results_readme_is_current():
    proc = subprocess.run(
        [sys.executable, "scripts/generators/build_results_readme.py", "--check"],
        cwd=str(data_paths.REPO_ROOT), capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_every_tracked_results_file_is_described_and_every_entry_matches():
    sys.path.insert(0, str(data_paths.REPO_ROOT / "scripts" / "generators"))
    import build_results_readme as brr

    files = brr.tracked_files()
    assert files, "git ls-files results returned nothing"
    assert all(brr.covering_key(rel) is not None for rel in files)
    # build() raises SystemExit naming every undescribed file or dangling entry
    text = brr.build()
    assert "Tracked files:" in text and "core_panel_scores" in text
