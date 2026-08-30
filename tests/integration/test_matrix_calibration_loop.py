import json
from pathlib import Path
from tempfile import TemporaryDirectory

from src.matrix_correction import MATRIX_CORRECTIONS, VOLATILE_CLASS_RETENTION_PROFILES, ProteinType
from src.matrix_experiment_intake import calibrate_from_intake

def test_calibration_loop_improves_mae():
    """
    Tests that the scipy L-BFGS-B optimizer successfully lowers the 
    Log-MAE error between pipeline predictions and experimental intake payloads.
    Uses a synthetic r11 SPI baseline payload with an artificially high target
    to force the optimizer to increase volatile retention / class factors.
    """
    payload = {
        "experiment_id": "test_r11_synthetic_hexanal",
        "source_kind": "synthetic_diagnostic",
        "protein_type": "soy_iso",
        "process_state": "extrusion",
        "conditions": {
            "temp_C": 140.0,
            "ph": 7.0,
            "water_activity": 0.5,
            "time_min": 1.0
        },
        "formulation": {
            "precursors": {
                "L-Cysteine": {"concentration_mM": 10.0},
                "D-Ribose": {"concentration_mM": 10.0}
            }
        },
        "measured_volatiles": {
            # Let's say MFT is measured at 500 ppb (very high for SPI)
            "2-methyl-3-furanthiol": {"conc_ppb": 500.0, "uncertainty_pct": 10},
            "2-furfurylthiol": {"conc_ppb": 200.0, "uncertainty_pct": 10}
        },
        "provenance": {
            "origin": "test_suite",
            "source_reference": "Internal Test"
        }
    }

    # Backup original constants
    orig_corr = MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE]
    orig_prof = VOLATILE_CLASS_RETENTION_PROFILES[ProteinType.SOY_ISOLATE]

    # Mock the guardrails so we don't run pytest inside pytest
    import src.matrix_calibration_optimizer
    original_guardrail = src.matrix_calibration_optimizer._run_guardrail_tests
    # 2026-08-27 (Wave I): `_run_guardrail_tests` now takes the candidate
    # constants (and optional test paths) so that the guardrail subprocess
    # validates the CANDIDATE rather than the on-disk baseline -- see FIX 17 in
    # src/matrix_calibration_optimizer.py. The stub has to accept them.
    src.matrix_calibration_optimizer._run_guardrail_tests = lambda *args, **kwargs: True

    try:
        # Run calibration
        result = calibrate_from_intake(payload, dry_run=False)

        assert result is not None, "Calibration failed or returned None"
        assert "baseline_log_mae" in result
        assert "calibrated_log_mae" in result
        assert "improvement_pct" in result
        
        # Assert MAE dropped
        assert result["calibrated_log_mae"] < result["baseline_log_mae"]
        assert result["improvement_pct"] > 0.0

    finally:
        # Restore state
        MATRIX_CORRECTIONS[ProteinType.SOY_ISOLATE] = orig_corr
        VOLATILE_CLASS_RETENTION_PROFILES[ProteinType.SOY_ISOLATE] = orig_prof
        src.matrix_calibration_optimizer._run_guardrail_tests = original_guardrail

