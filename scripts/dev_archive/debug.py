from src.matrix_calibration_optimizer import _compute_prediction_error
from src.matrix_correction import ProteinType
import numpy as np
from src.pipeline import MaillardPipeline
from src.smirks_engine import ReactionConditions

pl = MaillardPipeline("meaty", "beany")
cond = ReactionConditions(pH=7.0, temperature_celsius=140.0, water_activity=0.5, protein_type="soy_iso")
formulation = {
    "name": "dbg",
    "sugars": ["D-Glucose"],
    "amino_acids": ["L-Methionine"],
    "lipids": [],
    "molar_ratios": {"D-Glucose": 5.0, "L-Methionine": 5.0, "Hexanal": 5.0},
    "ph": 7.0,
    "temp": 140.0,
    "aw": 0.5,
    "time_minutes": 1.0,
    "interventions": [],
    "protein_type": "soy_iso",
    "denaturation_state": 1.0
}
try:
    res = pl.evaluate_single(formulation, cond)
    print("--- PREDICTED PPB MAP ---")
    for k,v in res.predicted_ppb.items():
        print(f"{k}: {v['value'] if isinstance(v, dict) else v}")
except Exception as e:
    print("EVAL FAILED:", e)

print("\n--- OPTIMIZATION ERROR EVAL ---")
payload = {
    "experiment_id": "test_r11_synthetic_hexanal",
    "source_kind": "synthetic_diagnostic",
    "protein_type": "soy_iso",
    "process_state": "extrusion",
    "conditions": {
        "temp_C": 140.0, "ph": 7.0, "water_activity": 0.5, "time_min": 1.0
    },
    "precursors": {
        "L-Methionine": {"concentration_mM": 5.0},
        "D-Glucose": {"concentration_mM": 5.0},
        "Hexanal": {"concentration_mM": 5.0}
    },
    "measured_volatiles": {
        "methional": {"conc_ppb": 500.0, "uncertainty_pct": 10},
        "hexanal": {"conc_ppb": 605.6, "uncertainty_pct": 10}
    }
}
print("Base MAE:", _compute_prediction_error([payload], ProteinType.SOY_ISOLATE, np.array([0.5, 1.0, 1.0, 1.0])))
print("High MAE:", _compute_prediction_error([payload], ProteinType.SOY_ISOLATE, np.array([1.0, 2.0, 1.0, 1.0])))
