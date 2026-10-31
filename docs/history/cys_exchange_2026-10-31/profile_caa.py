"""POST HOC diagnostic (labelled; not part of the verdict rule): how large must the amino-acid catalysis c_aa be for
E1 / E4 to pass the CF-076 core (T4 + T5 + T6 + T7)? c_aa pinned at each value, everything else searched in its band.
Single worker (the pinned band must live in this process)."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cys_exchange as m  # noqa: E402

OBS = m.observed()
CORE = ["T4", "T5", "T6", "T7"]


def best(rival, maxiter=20):
    dim = len({**m.COMMON, **m.CONFIG_KEYS[rival]})

    def f(u):
        r = m.predict((m.from_unit(rival, u), rival))
        return 1e3 if "error" in r else float(max(m.zmax(np.asarray(r[t])[None], OBS[t], t)[0] for t in CORE))

    res = differential_evolution(f, [(0, 1)] * dim, maxiter=maxiter, popsize=10, seed=5, polish=False, tol=1e-4)
    return float(res.fun), {k: float(f"{v:.3g}") for k, v in m.from_unit(rival, res.x).items()}


out = {}
for rival in sys.argv[1].split(","):
    for c in (1.0, 3.0, 10.0, 20.0):
        band = m.CONFIG_KEYS[rival]["c_aa"]
        m.CONFIG_KEYS[rival]["c_aa"] = (c, c, "pinned")
        z, p = best(rival)
        m.CONFIG_KEYS[rival]["c_aa"] = band
        out[f"{rival}_c_aa_{c:g}"] = {"best_core_max_z": z, "params": p}
        print(rival, c, round(z, 2), {k: p[k] for k in p if k in ("k0_100", "k_r45_100", "Ea", "c_P", "a")}, flush=True)
(m.OUT / "profile_caa.json").write_text(json.dumps(out, indent=1))
