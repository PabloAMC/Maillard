"""
tests/test_pbma_integration.py — Integration tests for Phase 7.3 PBMA enhancements.

Verifies:
1. Heme catalyst heuristic reduces barriers for Strecker and Pyrazine pathways.
2. Precursor reporting in the final recommendation table.
3. End-to-end pipeline execution for complex PBMA formulations.
"""

import pytest
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT = ROOT / "scripts" / "run_pipeline.py"

def test_heme_catalyst_heuristic():
    """
    Verify that the --catalyst heme flag reduces barriers for specific pathways.
    We compare a control run vs a heme run for Glucose + Glycine.
    """
    cmd_control = [sys.executable, str(PIPELINE_SCRIPT), "--sugars", "glucose", "--amino-acids", "glycine", "--ph", "7.0"]
    cmd_heme = [sys.executable, str(PIPELINE_SCRIPT), "--sugars", "glucose", "--amino-acids", "glycine", "--ph", "7.0", "--catalyst", "heme"]
    
    result_control = subprocess.run(cmd_control, capture_output=True, text=True, check=True)
    result_heme = subprocess.run(cmd_heme, capture_output=True, text=True, check=True)
    
    # We look for the 2,5-Dimethylpyrazine barrier in the output table.
    # Logic in run_pipeline.py: adjusted_bar -= 7.0 if catalyst == "heme"
    
    def extract_barrier(output, compound_name):
        for line in output.splitlines():
            if compound_name in line:
                # Extract e.g. "33.0 kcal"
                parts = line.split("│")
                if len(parts) > 3:
                    barrier_str = parts[3].strip().split()[0]
                    return float(barrier_str)
        return None

    barrier_control = extract_barrier(result_control.stdout, "2,5-Dimethylpyrazine")
    barrier_heme = extract_barrier(result_heme.stdout, "2,5-Dimethylpyrazine")
    
    assert barrier_control is not None
    assert barrier_heme is not None
    assert barrier_heme < barrier_control
    assert barrier_control - barrier_heme == pytest.approx(7.0)

def test_lipid_precursor_reporting():
    """
    Verify that input lipids (precursors) now appear in the results table as "COMPETING".
    """
    cmd = [sys.executable, str(PIPELINE_SCRIPT), "--sugars", "glucose", "--amino-acids", "glycine", "--lipids", "hexanal"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    assert "Hexanal" in result.stdout
    assert "[⚠️ COMPETING]" in result.stdout

def test_advanced_formulation_cli():
    """
    Verify that the CLI accepts all new arguments without crashing.
    """
    cmd = [
        sys.executable, str(PIPELINE_SCRIPT), 
        "--sugars", "ribose", 
        "--amino-acids", "cysteine", 
        "--additives", "thiamine,glutathione", 
        "--lipids", "hexanal", 
        "--catalyst", "heme",
        "--ph", "5.5"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert "Generated" in result.stdout
    assert "Predicted Targets" in result.stdout
