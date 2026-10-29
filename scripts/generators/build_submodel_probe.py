#!/usr/bin/env python
"""
build_ph_acid_submodel_probe.py -- sub-model 1's no-refit probe: the shipped engine against every verified pH and
acid row, before anything is fitted (docs/guides/OPERATING_PROCEDURE.md step 3; 2026-10-29).

Each pot of data/lit/fit_rows/ph_acid_submodel.yml runs once on the sulfur lane (the lane that carries the pH state
and the trunk's measured acid steps) at the shipped constants and pH coupling, on a grid of the pot's measured times.
Each row is then predicted by its observation operator (src/kinetic_core/observations.py): pH of the cooled sample,
pH after a 6x dilution, the pH-stat's NaOH uptake, titratable acid to pH 8.3, formic / acetic acid. Residuals are in
the fit's units (pH units; log10 ratio for amounts) and are reported by laboratory, data set, observable, temperature,
starting pH and time, raw and divided by each row's error.

Pots whose reactant the sulfur lane does not carry are skipped and listed: those flagged `engine_species_missing` and
those whose only reactants are outside the lane's state (Deng 2021's methionine ARP). Reactants the lane does not
carry beside one it does (casein lysine, hexanal, decadienal) are dropped and the pot is marked, not hidden. t = 0 pH
rows are the input the charge balance is built from; they are reported and kept out of the summaries.

    python scripts/generators/build_ph_acid_submodel_probe.py [--trunk-ph-dynamic] [--without KEY ...]
Writes results/validation/kinetic_core_ph_acid_submodel_probe.{json,md}. ``--trunk-ph-dynamic`` is the first candidate
arm, still with nothing refitted: the trunk lane's measured pH terms (trunk_conditions: Martins 2003 Table 3, the
Amadori-decay steps and the 3-DG -> formic exit), which the sulfur lane does not apply, evaluated at the pH the lane
computes (``sulfur.TRUNK_PH_TERMS``, switched on for this run only): ..._probe_trunk_ph_dynamic.{json,md}.
``--without KEY`` (repeatable) switches a step off (its k_ref set to zero), for a step whose basis the probe questions;
the arm's name records it (..._probe_trunk_ph_without_k_aa_frag.{json,md}).
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, fit_rows, fitting, observations, parallel, sulfur  # noqa: E402
from src.kinetic_core.ph_state import LABEL_TEMPERATURE_K, BufferSpec  # noqa: E402
from src.kinetic_core.sulfur import SULFUR_STATE_KEYS, integrate_sulfur  # noqa: E402

TABLE = "ph_acid_submodel"
OUT_STEM = data_paths.VALIDATION_DIR / "kinetic_core_ph_acid_submodel_probe"
CELSIUS = 273.15
#: Zhao 2017b: 0.5 mL of the reaction mixture made up to 3 mL before the pH was read
DILUTION = {"pH_6x_dilution_25C": 6.0}


def buffer_of(pot: Dict[str, Any]) -> BufferSpec:
    """The pot's BufferSpec from its recorded buffer: a pH-stat is clamped; a stated phosphate molarity is phosphate;
    'none' is a declared unbuffered pot. Protein (caseinate) buffering is not represented and is marked by the caller."""
    text = str(pot.get("buffer") or "")
    if "pH-stat" in text:
        return BufferSpec(kind="clamped", declared=True, source=text)
    m = re.search(r"phosphate\s+(\d+(?:\.\d+)?)\s*M", text)
    if m:
        return BufferSpec(kind="phosphate", phosphate_mol_l=float(m.group(1)), declared=True, source=text)
    if text.startswith("none"):
        return BufferSpec(kind="none", declared=True, source=text)
    raise ValueError(f"buffer {text!r} is not one the probe can map")


def plan(pots: Dict[str, dict], rows) -> List[Dict[str, Any]]:
    """One job per pot the lane can run; a skip record for the others."""
    jobs = []
    for pot_id, pot in pots.items():
        charge = {k: float(v) for k, v in (pot.get("initial_mmol_L") or {}).items() if v}
        carried = {k: v for k, v in charge.items() if k in SULFUR_STATE_KEYS}
        dropped = sorted(set(charge) - set(carried))
        times = sorted({0.0, *(r.t_min for r in rows if r.pot == pot_id)})
        job = {"pot": pot_id, "ledger_id": pot["ledger_id"], "t_c": float(pot["t_c"]), "label_ph": float(pot["label_ph"]),
               "charge": carried, "dropped": dropped, "times": times, "buffer": pot.get("buffer")}
        if pot.get("engine_species_missing"):
            job["skip"] = "engine_species_missing: " + ", ".join(dropped)
        elif not any(k in carried for k in ("Glc", "Fru", "AMA", "PENT", "ARP")):
            job["skip"] = "no sugar or Amadori compound the sulfur lane carries (" + ", ".join(dropped) + ")"
        notes = []
        if dropped and "skip" not in job:
            notes.append("not carried by the lane: " + ", ".join(dropped))
        if "caseinate" in str(pot.get("buffer")):
            notes.append("caseinate's own titratable groups not represented")
        job["notes"] = notes
        jobs.append(job)
    return jobs


def run_pot(job: Dict[str, Any]) -> Dict[str, Any]:
    """Integrate one pot on its time grid; return the state at every grid time and the pot's SID."""
    spec = buffer_of({"buffer": job["buffer"]})
    state0 = dict(job["charge"])
    state0.setdefault("OX", engine.OX_AMBIENT_MMOL_L)
    state0.setdefault("OXR", engine.OX_RESERVOIR_DEFAULT_UNITS)
    start = time.time()
    parameters = engine.core_parameters(engine.SULFUR)
    for key in job.get("without", ()):
        parameters[key] = dataclasses.replace(parameters[key], k_ref=0.0)
    run = integrate_sulfur(parameters, job["t_c"] + CELSIUS, state0, np.array(job["times"]),
                           ph=job["label_ph"], buffer_spec=spec, ph_drift=engine.core_ph_drift(), rtol=1e-8, atol=1e-14)
    states = [{k: float(run.series(k)[i]) for k in SULFUR_STATE_KEYS} for i in range(len(job["times"]))]
    return {"states": states, "sid": run.metadata.get("strong_ion_difference_mmol_L"),
            "converged": run.metadata.get("ph_fixed_point_converged"), "seconds": time.time() - start,
            "ph_in_situ": None if run.ph_series is None else [float(v) for v in run.ph_series]}


def predict(row, job: Dict[str, Any], result: Dict[str, Any]) -> float:
    """The row's observable from the predicted state, by its observation operator."""
    drift, spec = engine.core_ph_drift(), buffer_of({"buffer": job["buffer"]})
    i = job["times"].index(row.t_min)
    state, state0, sid = result["states"][i], result["states"][0], result["sid"]
    if row.species == "pH":
        return observations.ph_cooled(state, sid, drift, spec)
    if row.species in DILUTION:
        return observations.ph_diluted(state, sid, drift, spec, DILUTION[row.species])
    if row.species == "NaOH_uptake":
        # the pH-stat's electrode read and held the set pH at the reaction temperature (Davidek 2002)
        return observations.base_to_hold(job["label_ph"], state, state0, drift, _held(spec), job["t_c"] + CELSIUS)
    if row.species == "TA_to_pH8.3":
        return observations.base_to_hold(observations.TITRATION_END_PH, state, state0, drift, spec, LABEL_TEMPERATURE_K)
    if row.species == "FA+AA":
        return observations.formic_plus_acetic(state, drift)
    if row.species in ("FA", "AA"):
        return state[row.species]
    raise KeyError(f"no observation operator for {row.species!r}")


def _held(spec: BufferSpec) -> BufferSpec:
    """A pH-stat's buffer for the charge balance: the stated phosphate, or water; never 'clamped' (which has none)."""
    m = re.search(r"phosphate\s+(\d+(?:\.\d+)?)\s*M", spec.source)
    return BufferSpec(kind="phosphate", phosphate_mol_l=float(m.group(1)), declared=True) if m else BufferSpec(kind="none", declared=True)


class _Residuals:
    """fitting.diagnostics reads a problem's rows and evaluate(theta, rows); here the residuals are fixed."""

    def __init__(self, rows, values):
        self.rows, self._by = tuple(rows), {id(r): v for r, v in zip(rows, values)}

    def evaluate(self, theta, rows):
        return np.array([self._by[id(r)] for r in rows])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--trunk-ph-dynamic", action="store_true")
    parser.add_argument("--without", action="append", default=[])
    args = parser.parse_args(argv)
    sulfur.TRUNK_PH_TERMS = bool(args.trunk_ph_dynamic)   # set here, never at import; the workers fork after it
    arm = "_".join(["trunk_ph_dynamic"] * args.trunk_ph_dynamic + [f"without_{k}" for k in args.without]) or None
    out_json = Path(f"{OUT_STEM}{'_' + arm if arm else ''}.json")
    out_md = out_json.with_suffix(".md")
    start = time.time()
    pots, rows = fit_rows.load(TABLE)
    jobs = plan(pots, rows)
    for job in jobs:
        job["without"] = tuple(args.without)
    runnable = [j for j in jobs if "skip" not in j]
    results = dict(zip((j["pot"] for j in runnable), parallel.map_ordered(run_pot, runnable)))
    by_pot = {j["pot"]: j for j in jobs}
    records, scored = [], []
    for row in rows:
        job = by_pot[row.pot]
        if "skip" in job:
            continue
        predicted = predict(row, job, results[row.pot])
        residual = fit_rows.residual(row, predicted)
        record = {"pot": row.pot, "laboratory": row.laboratory, "data_set": job["ledger_id"], "observable": row.species,
                  "t_c": job["t_c"], "label_ph": job["label_ph"], "t_min": row.t_min, "measured": row.value,
                  "predicted": round(predicted, 6), "unit": row.unit, "sigma": round(row.sigma, 4),
                  "residual": round(residual, 4), "z": round(residual / row.sigma, 3),
                  "input_row": row.species == "pH" and row.t_min == 0.0}
        records.append(record)
        if not record["input_row"]:
            scored.append((row, record))
    attributes = {"laboratory": lambda r: r.laboratory, "data_set": lambda r: by_pot[r.pot]["ledger_id"],
                  "observable": lambda r: r.species, "temperature_C": lambda r: by_pot[r.pot]["t_c"],
                  "label_pH": lambda r: by_pot[r.pot]["label_ph"], "time_min": lambda r: r.t_min}
    scored_rows = [r for r, _ in scored]
    report = {
        "artifact": out_json.stem,
        "arm": arm or "shipped",
        "provenance": provenance.provenance_block(
            out_json.stem, generated_by="scripts/generators/build_ph_acid_submodel_probe.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{TABLE}.yml", ROOT / "src/kinetic_core/observations.py",
                    ROOT / "src/kinetic_core/ph_state.py", ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/parameters.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                    ROOT / "src/kinetic_core/trunk_conditions.py"]),
        "skipped": {j["pot"]: j["skip"] for j in jobs if "skip" in j},
        "notes": {j["pot"]: j["notes"] for j in jobs if j["notes"]},
        "runs": {p: {"seconds": round(r["seconds"], 1), "ph_fixed_point_converged": r["converged"],
                     "sid_mmol_L": r["sid"]} for p, r in results.items()},
        "raw": fitting.diagnostics(_Residuals(scored_rows, [rec["residual"] for _, rec in scored]), None, attributes),
        "whitened": fitting.diagnostics(_Residuals(scored_rows, [rec["z"] for _, rec in scored]), None, attributes),
        "rows": records,
        "wall_seconds": round(time.time() - start, 1),
    }
    out_json.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    out_md.write_text(markdown(report), encoding="utf-8")
    print(f"{len(scored)} rows scored on {len(results)} pots ({len(report['skipped'])} skipped) in {report['wall_seconds']} s "
          f"-> {data_paths.rel(out_md)}")
    return 0


def markdown(report: Dict[str, Any]) -> str:
    rows = [r for r in report["rows"] if not r["input_row"]]
    arm = "the shipped engine" if report["arm"] == "shipped" else "the shipped engine with " + " and ".join(
        "the trunk's measured pH terms applied at the lane's computed pH" if part == "trunk_ph_dynamic"
        else f"`{part[8:]}` switched off" for part in re.findall(r"trunk_ph_dynamic|without_k_[a-z_]+?(?=_without_|$)", report["arm"]))
    lines = [f"# Sub-model 1 no-refit probe: {arm} against the verified pH and acid rows", "",
             "*Generated by `scripts/generators/build_ph_acid_submodel_probe.py`. Shipped constants and pH coupling, "
             "nothing fitted. Residual = measured minus predicted, in pH units for pH and as log10(measured / predicted) "
             "for amounts; whitened = residual / the row's error. t = 0 pH rows (the charge balance's input) are listed "
             "in the JSON and kept out of these tables.*", "",
             "## Per data set (raw residuals)", "",
             "| data set | observable | n | mean | RMS | worst |", "|---|---|---|---|---|---|"]
    groups: Dict[tuple, list] = {}
    for r in rows:
        groups.setdefault((r["data_set"], r["observable"]), []).append(r["residual"])
    for (ds, ob), v in groups.items():
        worst = max(v, key=abs)
        lines.append(f"| {ds} | {ob} | {len(v)} | {np.mean(v):+.2f} | {math.sqrt(np.mean(np.square(v))):.2f} | {worst:+.2f} |")
    lines += ["", "## Skipped pots", ""] + [f"- `{p}`: {why}" for p, why in report["skipped"].items()]
    lines += ["", "## Pots run with a part of the system missing", ""] + \
             [f"- `{p}`: {'; '.join(n)}" for p, n in report["notes"].items()]
    lines += ["", fitting.diagnostics_markdown(report["raw"], "Raw residual diagnostics (pH units / dex)"),
              fitting.diagnostics_markdown(report["whitened"], "Whitened residual diagnostics (per row error)")]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
