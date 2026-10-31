import json, dataclasses
from src import data_paths
from src.kinetic_core import engine, parameters_sulfur as ps
rec = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_joint_deoxyosone_thiol_fit_s042_r3en_fixed.json").read_text())
fit = rec["fit"]; print("fit keys", list(fit)[:20])
for k in ("keys", "theta", "x", "optimum"):
    if k in fit: print(k, fit[k] if not isinstance(fit[k], dict) else {a: b for a, b in fit[k].items() if "mft" in a or "ttca" in a or "cama" in a or "cys" in a})
base = engine.core_parameters(engine.SULFUR)
for k in ("k_cys_ox", "k_cys_thermal"):
    p = base[k]; print(k, p.k_ref, p.ea_kj_mol)
print(dataclasses.asdict(ps.CysAmadoriRoute()))
print(dataclasses.asdict(ps.TtcaFormation()))
