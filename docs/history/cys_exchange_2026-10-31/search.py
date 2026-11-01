"""Feasibility search (prereg 'scoring' item 2; not a fit): min over each rival's bands of max |z| on a test set.
Primary = all scored tests (prereg). Sensitivity (POST HOC, labelled): the CF-076 core T4 + T5 + T6 + T7, and all
tests without T3 (T3's pH 7 point fails every rival, the same non-exchange xylose sink the prereg excluded at pH 8).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cys_exchange as m  # noqa: E402

OBS = m.observed()


class Objective:
    def __init__(self, rival, tests):
        self.rival, self.tests = rival, tests

    def __call__(self, u):
        r = m.predict((m.from_unit(self.rival, u), self.rival))
        if "error" in r:
            return 1e3
        return float(max(m.zmax(np.asarray(r[t])[None], OBS[t], t)[0] for t in self.tests))


def search(rival, tests, maxiter, seed=3):
    dim = len({**m.COMMON, **m.CONFIG_KEYS[rival]})
    res = differential_evolution(Objective(rival, tests), [(0, 1)] * dim, maxiter=maxiter, popsize=12, seed=seed,
                                 polish=False, workers=2, updating="deferred", tol=1e-4)
    p = m.from_unit(rival, res.x)
    r = m.predict((p, rival))
    return {"tests": tests, "best_max_z": float(res.fun), "evaluations": int(res.nfev),
            "z_by_test": {t: round(float(m.zmax(np.asarray(r[t])[None], OBS[t], t)[0]), 2) for t in m.TESTS},
            "params": {k: float(f"{v:.4g}") for k, v in p.items()},
            "pred": {t: np.asarray(r[t]).round(3).tolist() for t in m.TESTS}}


if __name__ == "__main__":
    rivals = sys.argv[1].split(",")
    maxiter = int(sys.argv[2])
    sets = {"primary_all": list(m.TESTS),
            "posthoc_core_T4567": ["T4", "T5", "T6", "T7"],
            "posthoc_all_but_T3": [t for t in m.TESTS if t != "T3"]}
    out_path = m.OUT / "search_results.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    for rv in rivals:
        for name, tests in sets.items():
            out.setdefault(rv, {})[name] = search(rv, tests, maxiter)
            print(rv, name, out[rv][name]["best_max_z"], out[rv][name]["z_by_test"], flush=True)
            out_path.write_text(json.dumps(out, indent=1))
