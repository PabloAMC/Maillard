"""
tests/unit/test_cli_scripts.py

Smoke tests to ensure that the command-line orchestrators parse arguments, 
load configs, and exit gracefully under expected conditions (like --help or dry runs).
"""

import subprocess
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

@pytest.mark.slow
def test_run_pipeline_help():
    """Verify run_pipeline.py parses --help successfully."""
    cmd = ["python", "scripts/run_pipeline.py", "--help"]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0
    assert "Maillard formulation screening pipeline" in result.stdout
    assert "--prediction-mode" in result.stdout

@pytest.mark.slow
def test_run_pipeline_dry_run():
    """Verify run_pipeline.py validates inputs correctly without running integration."""
    cmd = [
        "python", "scripts/run_pipeline.py",
        "--sugars", "ribose",
        "--amino-acids", "cysteine",
        "--prediction-mode", "kinetic",
        "--dry-run"
    ]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0
    assert "Dry-run complete" in result.stdout

@pytest.mark.slow
def test_run_cantera_help():
    """Verify run_cantera_kinetics.py parses --help successfully."""
    cmd = ["python", "scripts/run_cantera_kinetics.py", "--help"]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0
    assert "Run Cantera microkinetic simulation" in result.stdout

@pytest.mark.slow
def test_calibrate_barriers_importable():
    """Verify calibrate_barriers.py can be loaded without syntax or config errors."""
    cmd = ["python", "-c", "import scripts.calibrate_barriers"]
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0
