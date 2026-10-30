#!/usr/bin/env python
"""
diagnose_thiol_flux_at_joint_optimum.py -- where the thiol layer's misses come from (2026-10-30).

The joint sub-model 3 + 4 fit (fit_joint_deoxyosone_thiol.py) closes the dicarbonyls but leaves the absolute
thiol rows off in opposite directions in two ribose + cysteine laboratories at pH 5 - 5.5: Schieberle 2000
(Garching, SIDA) under by 0.85 dex in MFT and 2.0 in FFT, Yiltirak 2026 (Reading, SIDA) over by 2.35 and 1.9;
2-acetylthiazole in Zhou 2023's pot is about 1000x high. This generator integrates the diagnostic pots at the
fit's optimum on a log-dense grid and ranks, per species, the steps that make and use it (the integrated rate
over the run, ``sulfur.reaction_flux_integrals``; a budget that does not close to 1 % is refused), with the
precursors' half-times and every row's prediction against its measurement. It decides nothing: it says which
step a structural arm would have to change.

    python scripts/generators/diagnose_thiol_flux_at_joint_optimum.py [--fit s042]
Writes results/validation/kinetic_core_thiol_flux_diagnosis_<fit>.{json,md}.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for entry in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if entry not in sys.path:
        sys.path.insert(0, entry)

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import fit_rows, parameters_sulfur, sulfur  # noqa: E402

import build_submodel_probe as probe  # noqa: E402
import fit_joint_deoxyosone_thiol as joint  # noqa: E402
from build_sulfur_flux_attribution import _half_time, _ranked  # noqa: E402

#: (table, pot, why it is here)
POTS = (
    ("thiol_submodel", "schieberle2000_rib_cys_ph5_100C", "absolute, under-predicted (MFT -0.85, FFT -2.0)"),
    ("thiol_submodel", "yiltirak2026_rib_cys_ph55_buffer_100C", "absolute hold-out, over-predicted (MFT +2.35, FFT +1.9)"),
    ("thiol_submodel", "yiltirak2026_rib_cys_ph55_buffer_130C", "the same, short hot end of the equal-cook series"),
    ("thiol_submodel", "hofmann1998_cys_xylose_ph5", "absolute, Garching 145 C 20 min: the control that fits"),
    ("thiol_submodel", "hofmann1998_cys_ribose_ph3", "absolute hold-out, pH ladder low end"),
    ("thiol_submodel", "hofmann1998_cys_ribose_ph7", "absolute hold-out, pH ladder high end"),
    ("thiol_submodel", "zhai2023jafc_ttcacys_ph7_120C", "semi-quantitative, shape fits (sd 0.2)"),
    ("thiol_submodel", "zhou2022_axarp_cys_ph7_120C", "semi-quantitative, Zhou ARP + Cys"),
)
TRACK = ("Cys", "PENT", "ARP", "TTCA", "H2S", "NF", "FUR", "DPO", "MGO", "MFT", "FFT", "MFTD", "FFTD", "ACTZ", "OX")
CLOSURE_SPECIES = ("Cys", "PENT", "H2S", "MFT", "FFT")
CLOSURE_TOL = 0.01


def _grid(times):
    end = max(times)
    return sorted({0.0, *times, *np.geomspace(end * 1e-7, end, 301).tolist(), *np.linspace(0.0, end, 121).tolist()})


def _install(job) -> None:
    """The fit's engine state for this trial: route switches, NF -> MGO share, the Strecker H2S share."""
    joint.sm3.install(job)
    parameters_sulfur.STRECKER_SULFUR = dataclasses.replace(parameters_sulfur.STRECKER_SULFUR,
                                                           h2s_share=float(job["strecker_share"]))


def _diagnose(job, rows) -> Dict[str, Any]:
    measured_times = list(job["times"])
    job = dict(job, times=_grid(measured_times))
    run = probe.integrate_pot(job)
    flux = sulfur.reaction_flux_integrals(run)
    closure = sulfur.flux_closure(run, flux, CLOSURE_SPECIES)
    if max(closure.values()) > CLOSURE_TOL:
        raise SystemExit(f"{job['pot']}: flux budget does not close: {closure}")
    states = [{k: float(run.series(k)[i]) for k in probe.SULFUR_STATE_KEYS} for i in range(len(job["times"]))]
    result = {"states": states, "sid": run.metadata.get("strong_ion_difference_mmol_L")}
    predictions = []
    for r in rows:
        p = probe.predict(r, job, result)
        predictions.append({"species": r.species, "t_min": r.t_min, "unit": r.unit, "method": r.method,
                            "measured": r.value, "predicted": p,
                            "log10_ratio": math.log10(p / r.value) if p > 0 and r.value > 0 else None})
    ph = None if run.ph_series is None else [float(v) for v in run.ph_series]
    return {
        "t_c": job["t_c"], "label_ph": job["label_ph"], "charge_mmol_L": job["charge"], "minutes": max(measured_times),
        "half_time_min": {k: _half_time(run, k) for k in ("Cys", "PENT", "ARP", "TTCA") if k in job["charge"]},
        "ph_in_situ": None if ph is None else {"start": ph[0], "end": ph[-1]},
        "final_mmol_L": {k: run.final(k) for k in TRACK},
        "peak_mmol_L": {k: float(np.max(run.series(k))) for k in TRACK},
        "closure_relative_error": closure,
        "species": {k: _ranked(flux, k) for k in TRACK},
        "predictions": predictions,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--fit", default="s042", help="the joint fit's stem suffix (s042, s042_capped, s080)")
    args = parser.parse_args(argv)
    fit_path = data_paths.VALIDATION_DIR / f"kinetic_core_joint_deoxyosone_thiol_fit_{args.fit}.json"
    fit = json.loads(fit_path.read_text())
    theta = fit["fit"]["theta"]
    keys = tuple(theta)
    trial = joint.trial([theta[k] for k in keys], keys)
    share = float(fit["strecker_share"])
    sulfur.PH_COUPLING = "inline"

    out: Dict[str, Any] = {}
    loaded: Dict[str, Any] = {}
    for table, pot, why in POTS:
        if table not in loaded:
            loaded[table] = fit_rows.load(table)
        pots, rows = loaded[table]
        job = next(j for j in probe.plan(pots, rows) if j["pot"] == pot)
        job = dict(job, **trial, strecker_share=share)
        _install(job)
        out[pot] = dict(_diagnose(job, [r for r in rows if r.pot == pot and r.species in probe.SULFUR_STATE_KEYS]),
                        why=why, ledger_id=job["ledger_id"])
        print(f"  {pot}: closure {max(out[pot]['closure_relative_error'].values()):.1e}", flush=True)
    actz_job = dict(joint.ACTZ_FALSIFIER, **trial, strecker_share=share)
    _install(actz_job)
    out[actz_job["pot"]] = dict(_diagnose(actz_job, []), why="2-acetylthiazole falsifier (11.70 ug/L measured)",
                                ledger_id=actz_job["ledger_id"])

    stem = data_paths.VALIDATION_DIR / f"kinetic_core_thiol_flux_diagnosis_{args.fit}"
    payload = {"artifact": stem.name,
               "provenance": provenance.provenance_block(
                   stem.name, generated_by="scripts/generators/diagnose_thiol_flux_at_joint_optimum.py",
                   inputs=[fit_path, ROOT / "src/kinetic_core/sulfur.py"]),
               "fit": data_paths.rel(fit_path), "theta": theta, "strecker_share": share, "pots": out}
    stem.with_suffix(".json").write_text(json.dumps(payload, indent=2))

    lines = [f"# Thiol-layer flux diagnosis at the joint optimum ({args.fit})", "",
             f"*Generated by `scripts/generators/diagnose_thiol_flux_at_joint_optimum.py` from "
             f"`{data_paths.rel(fit_path)}`. Integrated reaction flux over each run (mmol/L); every budget closes to "
             f"{CLOSURE_TOL:.0%}. Ratios are predicted / measured in dex, before any semi-quantitative scale.*", ""]
    for pot, d in out.items():
        lines += [f"## {pot}", "", f"{d['why']}. {d['t_c']:g} C, label pH {d['label_ph']:g}, {d['minutes']:g} min, "
                  f"charge {d['charge_mmol_L']}; half-times (min) {d['half_time_min']}; in-situ pH {d['ph_in_situ']}.", ""]
        if d["predictions"]:
            lines += ["| species | t (min) | measured | predicted | log10 ratio |", "|---|---|---|---|---|"]
            lines += [f"| {p['species']} | {p['t_min']:g} | {p['measured']:.3g} {p['unit']} | {p['predicted']:.3g} | "
                      f"{'' if p['log10_ratio'] is None else format(p['log10_ratio'], '+.2f')} |" for p in d["predictions"]]
            lines.append("")
        for sp in ("H2S", "NF", "FUR", "MFT", "FFT", "ACTZ", "Cys"):
            s = d["species"][sp]
            for label in ("made", "used"):
                top = ", ".join(f"{t['reaction']} {t['share']:.0%}" for t in s[label]["top"][:4])
                lines.append(f"- {sp} {label} {s[label]['total_mmol_L']:.3g} mmol/L: {top}")
        lines += [f"- final (mmol/L): " + ", ".join(f"{k} {v:.3g}" for k, v in d["final_mmol_L"].items()), ""]
    stem.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {data_paths.rel(stem.with_suffix('.md'))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
