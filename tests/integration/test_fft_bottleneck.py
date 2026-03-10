"""
tests/integration/test_fft_bottleneck.py
"""
import pytest
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import tempfile
import pandas as pd
from scripts.run_cantera_kinetics import run_simulation

@pytest.mark.slow
def test_fft_bottleneck_resolution():
    """
    Test Phase 16: Verify FFT bottleneck is resolved in Ribose+Cys+Leu system.
    """
    precursors = {"ribose": 0.1, "cysteine": 0.05, "leucine": 0.05}
    
    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "fft_test")
        
        # We need the barrier DB
        db_path = "results/maillard_results.db"
        if not os.path.exists(db_path):
            pytest.skip("Barrier database not found. Skipping full simulation.")
            
        run_simulation(
            barriers_json=db_path,
            precursors=precursors,
            temp_c=150.0,
            time_sec=3600.0,
            from_smirks=True,
            track_species=["2-furfurylthiol", "H2S"],
            output_prefix=out_prefix,
            verbose_reactions=False,
            no_gating=True,
            pH=6.0,
            solvent="water"
        )
        
        csv_path = f"{out_prefix}_results.csv"
        assert os.path.exists(csv_path)
        
        df = pd.read_csv(csv_path)
        
        assert "2-furfurylthiol" in df.columns
        final_fft = df["2-furfurylthiol"].iloc[-1]
        
        # Minimum threshold enforced by the 2-step H2S pathway
        assert final_fft > 1e-35, f"FFT yield {final_fft} is too low, bottleneck might have regressed."
        
        assert "H2S" in df.columns
        final_h2s = df["H2S"].iloc[-1]
        assert final_h2s > 0, "H2S should be produced from cysteine."
