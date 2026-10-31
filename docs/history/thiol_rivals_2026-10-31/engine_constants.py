"""Read-only: the engine's constants for the thiol-loss and cysteine-fate rivals at the converged R3En (s042_r3en_fixed)."""
import json, math
from src import data_paths
from src.kinetic_core import engine
rec = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en_fixed.json").read_text())
theta = rec.get("theta") or rec.get("optimum") or {}
print("record keys:", list(rec)[:30])
base = engine.core_parameters(engine.SULFUR)
for k in ("k_mft_decay", "k_fft_decay", "k_dimer_mft", "k_ttca_cys", "k_ttca_deg", "k_cama_form", "k_cama_dpo", "k_cama_tdp", "k_cys_h2s", "k_h2s_loss"):
    if k in base:
        p = base[k]
        print(k, "k_ref", getattr(p, "k_ref", None), "ea", getattr(p, "ea_kj_mol", None), "T_ref", getattr(p, "t_ref_k", getattr(p, "temperature_ref_k", None)))
print([k for k in base if "cama" in k or "ttca" in k or k.startswith("k_cys")])
