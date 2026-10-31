#!/usr/bin/env python
"""
deoxypentosone_partition.py -- where 1-deoxypentosone goes, at a joint fit's optimum (2026-10-31; experiment
register EX-018, plan step B7). Read-only: integrated flux of the pentose, its Amadori compounds, the two
deoxypentosones and norfuraneol in Davidek 2008's and Desclaux 2006's glycine pots and Schieberle 2000's cysteine pot.
Finding recorded in kinetic_core_pentose_amine_entry_prereg.md sec. 9: in the glycine pots norfuraneol takes 38 - 53 %
of the 1-deoxypentosone (measured ~3 %), for want of the other exits Desclaux's network has.

    PYTHONPATH=. python scripts/diagnostics/deoxypentosone_partition.py [--fit s042_e]
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

POTS = (("ph_acid_davidek_acids", "davidek2008_xyl_gly_120c_ph60"), ("ph_acid_davidek_acids", "davidek2008_xyl_gly_120c_ph80"),
        ("deoxyosone_submodel", "desclaux2006_xyl_gly_ph6_100c"), ("thiol_submodel", "schieberle2000_rib_cys_ph5_100C"))


def main(argv=None) -> int:
    import numpy as np

    import build_submodel_probe as probe
    import fit_joint_deoxyosone_thiol as joint
    from build_sulfur_flux_attribution import _ranked
    from src import data_paths
    from src.kinetic_core import fit_rows, sulfur

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", default="s042_e")
    args = parser.parse_args(argv)
    fit = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json").read_text())
    theta = fit["fit"]["theta"]
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    sulfur.PH_COUPLING = "inline"
    for table, pot in POTS:
        pots, rows = fit_rows.load(table)
        job = next(j for j in probe.plan(pots, rows) if j["pot"] == pot)
        job = dict(job, **trial, strecker_share=fit["strecker_share"])
        job["times"] = sorted({0.0, *job["times"], *np.geomspace(1e-4, max(job["times"]), 200).tolist()})
        joint.install(job)
        run = probe.integrate_pot(job)
        flux = sulfur.reaction_flux_integrals(run)
        print(f"\n== {pot} (in-situ pH {run.ph_series[0]:.2f} -> {run.ph_series[-1]:.2f})")
        for sp in ("PENT", "PAMA", "DPO", "TDP", "NF"):
            r = _ranked(flux, sp)
            for side in ("made", "used"):
                t = r[side]
                print(f"  {sp} {side} {t['total_mmol_L']:.3g}: " + ", ".join(f"{x['reaction']} {x['share']:.0%}" for x in t["top"][:5]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
