#!/usr/bin/env python
"""
zhou2023_arp_acid_sources.py -- where the acid comes from in Zhou 2023's unbuffered ARP pots (2026-10-30; experiment
register EX-009). Read-only diagnosis at a joint fit's optimum: for each ARP (+ Cys) pot at initial pH 6 / 7 / 8, the
measured and predicted final pH and the steps that make the lane's acids (ACID, FA, AA, CBX). Finding recorded in
kinetic_core_thiol_submodel_prereg.md sec. 7: the lane makes ~3.5 mM net acid at every initial pH and has no pentose
route to formic or acetic acid.

    PYTHONPATH=. python scripts/diagnostics/zhou2023_arp_acid_sources.py [--fit s042]
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)


def main(argv=None) -> int:
    import numpy as np

    import build_submodel_probe as probe
    import fit_joint_deoxyosone_thiol as joint
    from build_sulfur_flux_attribution import _ranked
    from src import data_paths
    from src.kinetic_core import fit_rows, sulfur

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", default="s042")
    args = parser.parse_args(argv)
    fit = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json").read_text())
    theta = fit["fit"]["theta"]
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    sulfur.PH_COUPLING = "inline"
    pots, rows = fit_rows.load("ph_acid_submodel")
    for j in probe.plan(pots, rows):
        if not j["pot"].startswith("zhou2023") or "skip" in j:
            continue
        job = dict(j, **trial, strecker_share=fit["strecker_share"])
        job["times"] = sorted({0.0, *job["times"], *np.geomspace(1e-5, max(job["times"]), 200).tolist()})
        joint.install(job)
        run = probe.integrate_pot(job)
        flux = sulfur.reaction_flux_integrals(run)
        res = {"states": [{k: float(run.series(k)[i]) for k in probe.SULFUR_STATE_KEYS} for i in range(len(job["times"]))],
               "sid": run.metadata.get("strong_ion_difference_mmol_L")}
        row = next(r for r in rows if r.pot == j["pot"] and r.species == "pH")
        print(f"{j['pot']}: measured final pH {row.value}, predicted {probe.predict(row, job, res):.2f}")
        for acid in ("ACID", "FA", "AA", "CBX"):
            made = _ranked(flux, acid)["made"]
            print(f"  {acid} made {made['total_mmol_L']:.3g}: " + ", ".join(f"{t['reaction']} {t['share']:.0%}" for t in made["top"][:4]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
