import pytest
import sys
from src.recommend import main
from unittest.mock import patch
import io

def test_recommend_script_execution():
    """Ensure the prototype integration script executes without crashing."""
    
    # Capture standard output to prevent cluttering test logs
    captured_output = io.StringIO()
    
    with patch('sys.stdout', new=captured_output):
        # Using patch on SkalaRefiner.refine_barrier just in case PySCF isn't fully set up on the test system
        with patch('src.skala_refiner.SkalaRefiner.refine_barrier', return_value=15.0):
            main()
            
    output = captured_output.getvalue()
    
    # Assertions on expected stages being reached
    assert "Maillard Reaction Pathway Recommender" in output
    assert "Setting Conditions: pH 6.5" in output
    assert "Extracting Pathways" in output
    assert "xTB Screening" in output
    assert "Savory_FFT_Pathway" in output
    assert "DFT/Skala Refinement of Top Candidates" in output
    assert "Final Recommendation: Savory_FFT_Pathway is highly favoured" in output
