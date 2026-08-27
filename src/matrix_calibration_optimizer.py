import json
import logging
import math
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize

from src.matrix_correction import (
    MATRIX_CORRECTIONS,
    VOLATILE_CLASS_RETENTION_PROFILES,
    ProteinType,
    MatrixCorrection,
    VolatileClassRetentionProfile,
    ACCESSIBILITY_LITERATURE_WINDOWS
)
from src.pipeline import MaillardPipeline
from src.smirks_engine import ReactionConditions

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
CALIBRATION_HISTORY_DIR = ROOT / "data" / "calibration_history"

# Measurements below ~0.1 ppb are at or under the quantitation limit of the GC-O
# methods the intake payloads come from, so the measurement side of the log error is
# floored there. Predictions are floored only to keep log10 defined (see the dated
# note in `_compute_prediction_error`): flooring them at a physically meaningful
# concentration flattens the objective and silently disables the calibration loop.
MEASUREMENT_QUANTITATION_FLOOR_PPB = 0.1
PREDICTION_LOG_GUARD_PPB = 1.0e-12


def _monkey_patch_matrix_constants(
    protein_type: ProteinType,
    volatile_retention_denatured: float,
    aldehyde_factor: float,
    sulfur_factor: float,
    pyrazine_factor: float,
):
    """
    Temporarily overrides the global matrix constants in src.matrix_correction.
    """
    # 1. Update MATRIX_CORRECTIONS
    if protein_type in MATRIX_CORRECTIONS:
        orig = MATRIX_CORRECTIONS[protein_type]
        MATRIX_CORRECTIONS[protein_type] = MatrixCorrection(
            protein_type=orig.protein_type,
            lysine_accessibility=orig.lysine_accessibility,
            cysteine_accessibility=orig.cysteine_accessibility,
            lysine_accessibility_native=orig.lysine_accessibility_native,
            lysine_accessibility_denatured=orig.lysine_accessibility_denatured,
            cysteine_accessibility_native=orig.cysteine_accessibility_native,
            cysteine_accessibility_denatured=orig.cysteine_accessibility_denatured,
            volatile_retention=orig.volatile_retention,
            volatile_retention_native=orig.volatile_retention_native,
            volatile_retention_denatured=volatile_retention_denatured,
            source=orig.source
        )

    # 2. Update VOLATILE_CLASS_RETENTION_PROFILES
    if protein_type in VOLATILE_CLASS_RETENTION_PROFILES:
        orig_prof = VOLATILE_CLASS_RETENTION_PROFILES[protein_type]
        denatured_factors = dict(orig_prof.denatured_factors)
        denatured_factors["aldehyde"] = aldehyde_factor
        denatured_factors["sulfur"] = sulfur_factor
        denatured_factors["pyrazine"] = pyrazine_factor

        VOLATILE_CLASS_RETENTION_PROFILES[protein_type] = VolatileClassRetentionProfile(
            protein_type=orig_prof.protein_type,
            native_factors=orig_prof.native_factors,
            denatured_factors=denatured_factors,
            source=orig_prof.source
        )


def _compute_prediction_error(
    experiments: List[Dict[str, Any]],
    protein_type: ProteinType,
    x: np.ndarray
) -> float:
    """
    x[0] = volatile_retention_denatured
    x[1] = aldehyde_factor
    x[2] = sulfur_factor
    x[3] = pyrazine_factor
    """
    _monkey_patch_matrix_constants(
        protein_type,
        volatile_retention_denatured=x[0],
        aldehyde_factor=x[1],
        sulfur_factor=x[2],
        pyrazine_factor=x[3]
    )

    pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    
    total_error = 0.0
    count = 0

    for exp in experiments:
        cond = ReactionConditions(
            pH=float(exp["conditions"]["ph"]),
            temperature_celsius=float(exp["conditions"]["temp_C"]),
            water_activity=float(exp["conditions"]["water_activity"]),
            protein_type=exp["protein_type"]
        )
        
        flat_molar_ratios = {}
        for k, v in exp["precursors"].items():
            flat_molar_ratios[k] = float(v.get("concentration_mM", 1.0))
            
        sugars = []
        amino_acids = []
        additives = []
        for k in flat_molar_ratios.keys():
            kl = k.lower()
            if "ose" in kl or "sugar" in kl:
                sugars.append(k)
            elif "ine" in kl or "acid" in kl:
                amino_acids.append(k)
            else:
                additives.append(k)

        formulation = {
            "name": exp.get("benchmark_id", "calib_exp"),
            "sugars": sugars,
            "amino_acids": amino_acids,
            "lipids": [],
            "additives": additives,
            "molar_ratios": flat_molar_ratios,
            "ph": cond.pH,
            "temp": cond.temperature_celsius,
            "aw": cond.water_activity,
            "time_minutes": float(exp["conditions"]["time_min"]),
            "interventions": [],
            "protein_type": cond.protein_type,
            "denaturation_state": float(exp.get("denaturation_state", 1.0))
        }

        # Predict
        res = pipeline.evaluate_single(formulation, cond)
        predicted_ppb_map = getattr(res, "predicted_ppb", {})
        if not isinstance(predicted_ppb_map, dict):
            predicted_ppb_map = {}

        measured = exp.get("measured_volatiles", {})
        for compound, data in measured.items():
            meas_val = float(data.get("conc_ppb", 0.0))
            if meas_val <= 0:
                continue
                
            pred_val = float(predicted_ppb_map.get(compound, 0.0))
            if isinstance(pred_val, dict):
                pred_val = float(pred_val.get("value", 0.0))
                
            # Log-scaled Mean Absolute Error (Log-MAE)
            # We use log10 to handle the huge dynamic range of ppb (from 0.1 to 10000)
            #
            # AUDIT 2026-08-27 (Wave H). The prediction side used to be clamped at the
            # same 0.1 ppb floor as the measurement. That floor is right for a
            # MEASUREMENT (0.1 ppb is around the quantitation limit of the GC-O methods
            # these payloads come from) and wrong for a PREDICTION: it makes the
            # objective exactly constant for every prediction below 0.1 ppb, so the
            # L-BFGS-B gradient is identically zero and `calibrate_matrix_constants`
            # reports "did not improve the MAE" and reverts — silently, for a reason that
            # has nothing to do with the calibration. After the Wave G1 chemistry rebuild
            # dropped sulfur yields 5-40x, that is the normal case, not an edge case:
            # tests/integration/test_matrix_calibration_loop.py went red with a perfectly
            # flat objective at MAE 3.5 = (|log10(0.1/500)| + |log10(0.1/200)|)/2, i.e. the
            # value of the clamp itself rather than of any prediction.
            # The prediction floor is now a pure log(0) guard placed far below any
            # physically meaningful concentration, so under-prediction stays visible to
            # the optimiser as a gradient instead of being flattened into a constant.
            log_meas = math.log10(max(meas_val, MEASUREMENT_QUANTITATION_FLOOR_PPB))
            log_pred = math.log10(max(pred_val, PREDICTION_LOG_GUARD_PPB))
            
            total_error += abs(log_pred - log_meas)
            count += 1

    if count == 0:
        return 1e6
        
    return total_error / count


def _run_guardrail_tests() -> bool:
    """Runs the scientific benchmarks to ensure we haven't broken free-precursor kinetics."""
    logger.info("Running guardrail tests: pytest tests/scientific/test_benchmarks.py")
    try:
        env = {"MAILLARD_STRICT_BENCHMARKS": "1"}
        import os
        full_env = os.environ.copy()
        full_env.update(env)
        
        result = subprocess.run(
            ["pytest", "tests/scientific/test_benchmarks.py", "-q"],
            env=full_env,
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            return True
        else:
            logger.warning(f"Guardrail tests failed:\n{result.stdout}\n{result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error running guardrail tests: {e}")
        return False


def calibrate_matrix_constants(experiments: List[Dict[str, Any]], protein_type_str: str) -> Optional[Dict[str, Any]]:
    """
    Minimizes the error between predicted and measured ppb across all experiments
    for a specific protein type, updating matrix retention constants.
    """
    try:
        protein_type = ProteinType(protein_type_str)
    except ValueError:
        logger.error(f"Unsupported protein type for calibration: {protein_type_str}")
        return None

    orig_corr = MATRIX_CORRECTIONS.get(protein_type)
    orig_prof = VOLATILE_CLASS_RETENTION_PROFILES.get(protein_type)
    
    if not orig_corr or not orig_prof:
        logger.error(f"Missing matrix prior entries for {protein_type_str}")
        return None

    # Initial guess
    x0 = np.array([
        orig_corr.volatile_retention_denatured,
        orig_prof.denatured_factors.get("aldehyde", 1.0),
        orig_prof.denatured_factors.get("sulfur", 1.0),
        orig_prof.denatured_factors.get("pyrazine", 1.0)
    ])

    # Bounds: fractions must be between 0.01 and 1.0 (or >1 for generation, but retention is <= 1.0 generally)
    bounds = [
        (0.01, 1.0),   # volatile_retention_denatured
        (0.01, 2.0),   # aldehyde_factor (can be slightly > 1 if lipid oxidation contributes)
        (0.01, 1.0),   # sulfur_factor (strongly trapped)
        (0.01, 1.0)    # pyrazine_factor
    ]

    # Calculate pre-calibration baseline error
    baseline_mae = _compute_prediction_error(experiments, protein_type, x0)
    logger.info(f"Baseline Log-MAE: {baseline_mae:.4f}")

    logger.info("Starting L-BFGS-B optimization...")
    start_time = time.time()
    
    res = minimize(
        fun=lambda x: _compute_prediction_error(experiments, protein_type, x),
        x0=x0,
        bounds=bounds,
        method="L-BFGS-B",
        options={"ftol": 1e-4, "maxiter": 50}
    )

    elapsed = time.time() - start_time
    logger.info(f"Optimization finished in {elapsed:.1f}s. Success: {res.success}. Final Log-MAE: {res.fun:.4f}")

    if res.fun >= baseline_mae:
        logger.warning("Optimization did not improve the MAE. Reverting to baseline.")
        # Revert
        _compute_prediction_error(experiments, protein_type, x0)
        return None

    # Test guardrails with new parameters
    guardrails_passed = _run_guardrail_tests()
    if not guardrails_passed:
        logger.warning("Guardrail tests failed with new parameters. Reverting to baseline.")
        # Revert
        _compute_prediction_error(experiments, protein_type, x0)
        return None

    # Success! Write diff to data/calibration_history
    mae_drop_pct = (baseline_mae - res.fun) / baseline_mae * 100
    logger.info(f"Calibration successful! MAE dropped by {mae_drop_pct:.1f}%. Guardrails passed.")

    CALIBRATION_HISTORY_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    diff = {
        "timestamp": datetime.now().isoformat(),
        "protein_type": protein_type_str,
        "baseline_log_mae": float(baseline_mae),
        "calibrated_log_mae": float(res.fun),
        "improvement_pct": float(mae_drop_pct),
        "parameters": {
            "volatile_retention_denatured": {"old": float(x0[0]), "new": float(res.x[0])},
            "aldehyde_factor": {"old": float(x0[1]), "new": float(res.x[1])},
            "sulfur_factor": {"old": float(x0[2]), "new": float(res.x[2])},
            "pyrazine_factor": {"old": float(x0[3]), "new": float(res.x[3])}
        },
        "experiments_used": [exp.get("benchmark_id", "unknown") for exp in experiments]
    }
    
    out_path = CALIBRATION_HISTORY_DIR / f"calib_{protein_type_str}_{timestamp}.json"
    out_path.write_text(json.dumps(diff, indent=2))
    logger.info(f"Wrote calibration history to {out_path}")

    return diff
