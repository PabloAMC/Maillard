#!/usr/bin/env python
"""
build_submodel_probe.py -- a sub-model's no-refit probe: the shipped engine against every verified row of its keyed
fit table, before anything is fitted (docs/guides/OPERATING_PROCEDURE.md step 3; 2026-10-29). Written for sub-model 1
(pH and acid); any READY table whose species the sulfur lane carries runs the same way (sub-model 2, the route).

Each pot of data/lit/fit_rows/<table>.yml runs once on the sulfur lane (the lane that carries the pH state
and the trunk's measured acid steps) at the shipped constants and pH coupling, on a grid of the pot's measured times.
Each row is then predicted by its observation operator (src/kinetic_core/observations.py): pH of the cooled sample,
pH after a 6x dilution, the pH-stat's NaOH uptake, titratable acid to pH 8.3, formic / acetic acid. Residuals are in
the fit's units (pH units; log10 ratio for amounts) and are reported by laboratory, data set, observable, temperature,
starting pH and time, raw and divided by each row's error. A row whose species is a state variable of the lane (TTCA,
cysteine, pentose ...) is compared with that variable directly.

Pots whose reactant the sulfur lane does not carry are skipped and listed: those flagged `engine_species_missing` and
those with no reactant in the lane's state (Deng 2021's methionine ARP). Reactants the lane does not
carry beside one it does (casein lysine, hexanal, decadienal) are dropped and the pot is marked, not hidden. t = 0 pH
rows are the input the charge balance is built from; they are reported and kept out of the summaries.

    python scripts/generators/build_submodel_probe.py [--table NAME] [--trunk-ph-dynamic] [--without KEY ...]
Writes results/validation/kinetic_core_<table>_probe.{json,md} (default table ph_acid_submodel). ``--trunk-ph-dynamic`` is the first candidate
arm, still with nothing refitted: the trunk lane's measured pH terms (trunk_conditions: Martins 2003 Table 3, the
Amadori-decay steps and the 3-DG -> formic exit), which the sulfur lane does not apply, evaluated at the pH the lane
computes (``sulfur.TRUNK_PH_TERMS``, switched on for this run only): ..._probe_trunk_ph_dynamic.{json,md}.
``--switch-on NAME`` (repeatable) enables a pre-registered step of ``parameters_sulfur`` that ships off (B58
TTCA_FORMATION, B59b CYS_AMADORI_ROUTE, B59a STRECKER_SULFUR, B61 NF_HYDROXIDE_DECAY), for this run only.
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
from src.kinetic_core import engine, fit_rows, fitting, observations, parallel, parameters_sulfur, sulfur  # noqa: E402
from src.kinetic_core import ph_state  # noqa: E402
from src.kinetic_core.ph_state import LABEL_TEMPERATURE_K, BufferSpec  # noqa: E402
from src.kinetic_core.sulfur import SULFUR_STATE_KEYS, integrate_sulfur  # noqa: E402

CELSIUS = 273.15
#: Zhao 2017b: 0.5 mL of the reaction mixture made up to 3 mL before the pH was read
DILUTION = {"pH_6x_dilution_25C": 6.0}
#: observables computed from the state rather than read off one state variable
OPERATORS = ("pH", "NaOH_uptake", "TA_to_pH8.3", "FA+AA", *DILUTION)
#: pre-registered steps that ship off, by switch name -> the short name an arm carries
SWITCHES = {"TTCA_FORMATION": "b58", "CYS_AMADORI_ROUTE": "b59b", "STRECKER_SULFUR": "b59a", "NF_HYDROXIDE_DECAY": "b61"}
#: a series whose absolute level is unknown (fitting.SCALE_TAU = inf) is compared in shape only
SHAPE_ONLY = ("semi_quantitative", "relative")


#: McIlvaine (1921) citrate-phosphate: mL of 0.2 M Na2HPO4 per 20 mL of mix (the rest 0.1 M citric acid) by pH, the
#: standard recipe; used where a paper names McIlvaine buffer and its pH but not its mixing ratio.
MCILVAINE_ML_PHOSPHATE = ((3.0, 4.11), (4.0, 7.71), (5.0, 10.30), (6.0, 12.63), (7.0, 16.47), (8.0, 19.45))
_MOLARITY = r"(\d+(?:\.\d+)?)\s*M\b"


class UnmappedBuffer(ValueError):
    """A buffer whose composition the record does not state well enough to compute its capacity."""


def buffer_of(pot: Dict[str, Any], label_ph: float = None) -> BufferSpec:
    """The pot's BufferSpec from its recorded buffer: a pH-stat is clamped; a phosphate with a stated molarity is
    phosphate; McIlvaine at its stated pH is citrate-phosphate on the standard recipe; 'none' is a declared unbuffered
    pot. A phosphate of unprinted strength raises UnmappedBuffer (the pot is skipped, not guessed). Protein (caseinate)
    buffering is not represented and is marked by the caller."""
    text = str(pot.get("buffer") or "")
    if "pH-stat" in text:
        return BufferSpec(kind="clamped", declared=True, source=text)
    if "McIlvaine" in text:
        if label_ph is None:
            raise UnmappedBuffer(f"McIlvaine buffer without a pH: {text!r}")
        ml = float(np.interp(label_ph, *zip(*MCILVAINE_ML_PHOSPHATE)))
        return BufferSpec(kind="citrate_phosphate", phosphate_mol_l=0.2 * ml / 20.0, citrate_mol_l=0.1 * (20.0 - ml) / 20.0,
                          declared=True, source=text + " [mixing ratio from McIlvaine's table at the pot's pH]")
    if "phosphate" in text or "PO4" in text:
        m = re.search(_MOLARITY, text)
        if not m:
            raise UnmappedBuffer(f"phosphate of unprinted strength: {text!r}")
        return BufferSpec(kind="phosphate", phosphate_mol_l=float(m.group(1)), declared=True, source=text)
    if text.startswith("none"):
        return BufferSpec(kind="none", declared=True, source=text)
    raise UnmappedBuffer(f"buffer {text!r} is not one the probe can map")


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
        elif str(pot.get("engine_precursor_mapping", "")).startswith("unresolved"):
            job["skip"] = "engine_precursor_mapping " + pot["engine_precursor_mapping"]
        elif not carried:
            job["skip"] = "no reactant the sulfur lane carries (" + ", ".join(dropped) + ")"
        else:
            try:
                buffer_of(pot, job["label_ph"])
            except UnmappedBuffer as exc:
                job["skip"] = str(exc)
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
    start = time.time()
    run = integrate_pot(job)
    states = [{k: float(run.series(k)[i]) for k in SULFUR_STATE_KEYS} for i in range(len(job["times"]))]
    return {"states": states, "sid": run.metadata.get("strong_ion_difference_mmol_L"),
            "converged": run.metadata.get("ph_fixed_point_converged"), "seconds": time.time() - start,
            "ph_in_situ": None if run.ph_series is None else [float(v) for v in run.ph_series]}


def integrate_pot(job: Dict[str, Any]):
    """The sulfur lane's run for one pot on its time grid, with the job's trial constants and switches."""
    spec = buffer_of({"buffer": job["buffer"]}, job["label_ph"])
    state0 = dict(job["charge"])
    state0.setdefault("OX", engine.OX_AMBIENT_MMOL_L)
    state0.setdefault("OXR", engine.OX_RESERVOIR_DEFAULT_UNITS)
    parameters = engine.core_parameters(engine.SULFUR)
    for key in job.get("without", ()):
        parameters[key] = dataclasses.replace(parameters[key], k_ref=0.0)
    for key, k_ref in (job.get("k_ref") or {}).items():   # a fit's trial value (fit_ph_acid_submodel.py)
        parameters[key] = dataclasses.replace(parameters[key], k_ref=float(k_ref))
    for key, ea in (job.get("ea") or {}).items():
        parameters[key] = dataclasses.replace(parameters[key], ea_kj_mol=float(ea))
    clamp_label = job.get("clamp_label_ph")   # a pH-stat set on the label scale and held at its hot equivalent
    if spec.is_clamped and clamp_label is None:
        # a pH-stat holds the set pH at the reaction temperature (Davidek 2002: "measured at 90 C"); the measured
        # pH-indexed tables read the same state's pH on the 25 C label scale (design review B3)
        held, t_hot = _held(spec), job["t_c"] + CELSIUS
        sid = ph_state.strong_ion_difference(job["label_ph"], state0, _drift(job), held, label_temperature_k=t_hot)
        clamp_label = observations.ph_cooled(state0, sid, _drift(job), held)
    return integrate_sulfur(parameters, job["t_c"] + CELSIUS, state0, np.array(job["times"]),
                            ph=job["label_ph"], buffer_spec=spec, ph_drift=_drift(job), rtol=1e-8, atol=1e-14,
                            clamp_ph_label=clamp_label)


def _drift(job: Dict[str, Any]):
    """The pH-drift calibration: a fit's trial values when the job carries them, the frozen one otherwise."""
    return job.get("drift") or engine.core_ph_drift()


class PotTimeout(RuntimeError):
    """A pot integration over its wall-clock limit (a stiff corner of parameter space)."""


def run_pot_limited(job: Dict[str, Any], seconds: int = 120) -> Dict[str, Any]:
    """``run_pot`` under a wall-clock limit (2026-10-30: one stiff pot stalled a fit's pool for an hour). Unix
    alarm, so it works in forked workers; on time-out it raises PotTimeout, which a fit counts as a failed run."""
    import signal

    def _alarm(signum, frame):
        raise PotTimeout(f"pot {job.get('pot')} exceeded {seconds} s")

    previous = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(int(seconds))
    try:
        return run_pot(job)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def predict(row, job: Dict[str, Any], result: Dict[str, Any]) -> float:
    """The row's observable from the predicted state, by its observation operator."""
    drift, spec = _drift(job), buffer_of({"buffer": job["buffer"]}, job["label_ph"])
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
    if row.species in state:
        return observations.amount_in(row.unit, row.species, state[row.species], state0[row.species])
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
    parser.add_argument("--table", default="ph_acid_submodel")
    parser.add_argument("--trunk-ph-dynamic", action="store_true")
    parser.add_argument("--trunk-solutes-in-balance", action="store_true",
                        help="glycine, the Amadori compound and the Schiff base in the charge balance (design review B4)")
    parser.add_argument("--without", action="append", default=[])
    parser.add_argument("--switch-on", action="append", default=[], choices=SWITCHES)
    parser.add_argument("--scale-group", choices=("series", "figure"), default="series")
    parser.add_argument("--k-multiplier", action="append", default=[], metavar="KEY=LOG10",
                        help="a fitted arm's constant applied unfitted (k_ref x 10^LOG10), e.g. k_dpo_ddp=-4")
    parser.add_argument("--nf-mgo-share", type=float, default=None,
                        help="NF_HYDROXIDE_DECAY.mgo_share (sub-model 3), applied unfitted")
    parser.add_argument("--label", default=None, help="a name for this arm in the artifact's file name")
    parser.add_argument("--scale-tau", type=float, default=float("inf"))
    parser.add_argument("--ph-coupling", choices=("fixed_point_nodes", "inline"), default=sulfur.PH_COUPLING,
                        help="the shipped coupling by default; 'inline' to check it against the shipped one on this table")
    args = parser.parse_args(argv)
    sulfur.PH_COUPLING = args.ph_coupling
    ph_state.TRUNK_SOLUTES_IN_BALANCE = bool(args.trunk_solutes_in_balance)
    for name in args.switch_on:   # set here, never at import; the workers fork after it
        setattr(parameters_sulfur, name, dataclasses.replace(getattr(parameters_sulfur, name), enabled=True))
    if args.nf_mgo_share is not None:
        parameters_sulfur.NF_HYDROXIDE_DECAY = dataclasses.replace(parameters_sulfur.NF_HYDROXIDE_DECAY,
                                                                   mgo_share=float(args.nf_mgo_share))
    base = engine.core_parameters(engine.SULFUR)
    k_ref_override = {k: base[k].k_ref * 10.0 ** float(v) for k, v in (kv.split("=") for kv in args.k_multiplier)}
    sulfur.TRUNK_PH_TERMS = bool(args.trunk_ph_dynamic)   # set here, never at import; the workers fork after it
    arm = "_".join(([args.label] if args.label else []) + ["trunk_ph_dynamic"] * args.trunk_ph_dynamic + ["trunk_solutes"] * args.trunk_solutes_in_balance
                   + ([f"figure_scale_tau{args.scale_tau:g}"] if args.scale_group == "figure" else [])
                   + [f"with_{SWITCHES[n]}" for n in args.switch_on]
                   + [f"without_{k}" for k in args.without] + (["inline"] if args.ph_coupling == "inline" else [])) or None
    out_json = data_paths.VALIDATION_DIR / f"kinetic_core_{args.table}_probe{'_' + arm if arm else ''}.json"
    out_md = out_json.with_suffix(".md")
    start = time.time()
    pots, rows = fit_rows.load(args.table)
    jobs = plan(pots, rows)
    for job in jobs:
        job["without"] = tuple(args.without)
        job["k_ref"] = dict(k_ref_override)
    runnable = [j for j in jobs if "skip" not in j]
    results = dict(zip((j["pot"] for j in runnable), parallel.map_ordered(run_pot, runnable)))
    by_pot = {j["pot"]: j for j in jobs}
    records, scored, unpredictable = [], [], {}
    for row in rows:
        job = by_pot[row.pot]
        if "skip" in job:
            continue
        if row.species not in OPERATORS and row.species not in SULFUR_STATE_KEYS:
            unpredictable[row.species] = unpredictable.get(row.species, 0) + 1
            continue
        predicted = predict(row, job, results[row.pot])
        if not predicted > 0.0 and row.log_scale:
            unpredictable[f"{row.species} (predicted zero)"] = unpredictable.get(f"{row.species} (predicted zero)", 0) + 1
            continue
        residual = fit_rows.residual(row, predicted)
        record = {"pot": row.pot, "laboratory": row.laboratory, "data_set": job["ledger_id"], "observable": row.species,
                  "t_c": job["t_c"], "label_ph": job["label_ph"], "t_min": row.t_min, "measured": row.value,
                  "predicted": round(predicted, 6), "unit": row.unit, "sigma": float(f"{row.sigma:.4g}"),
                  "residual": round(residual, 4), "z": round(residual / row.sigma, 3),
                  "input_row": row.species == "pH" and row.t_min == 0.0}
        record["method"] = row.method
        records.append(record)
        if not record["input_row"]:
            scored.append((row, record))
    # semi-quantitative series: their unknown scale estimated and removed. Per series (ledger|pot|species, tau = inf)
    # by default: the series' mean log residual. With --scale-group figure, one scale per (ledger id, species) across
    # the figure's pots -- keeping the cross-pot contrasts (pH, temperature, +/- cysteine) -- shrunk toward 0 by a
    # prior of width --scale-tau dex: the posterior mean of the offset, tau^2 S / (1 + tau^2 W), S = sum r / s^2,
    # W = sum 1 / s^2 (sub-model 3 design review, B2).
    def group_of(row):
        return row.series if args.scale_group == "series" else f"{row.ledger_id}|{row.species}"

    by_series: Dict[str, list] = {}
    for row, record in scored:
        if row.method in SHAPE_ONLY:
            by_series.setdefault(group_of(row), []).append(record)
    for series in by_series.values():
        if math.isinf(args.scale_tau):
            offset = float(np.mean([r["residual"] for r in series]))
        else:
            w = np.array([1.0 / r["sigma"] ** 2 for r in series])
            ssum = float(np.sum(w * np.array([r["residual"] for r in series])))
            offset = args.scale_tau ** 2 * ssum / (1.0 + args.scale_tau ** 2 * float(np.sum(w)))
        for r in series:
            r["scale_offset_dex"] = round(offset, 4)
            r["residual"] = round(r["residual"] - offset, 4)
            r["z"] = round(r["residual"] / r["sigma"], 3)
    scored = [(row, rec) for row, rec in scored if not (row.method in SHAPE_ONLY and len(by_series[group_of(row)]) < 2)]
    attributes = {"laboratory": lambda r: r.laboratory, "data_set": lambda r: by_pot[r.pot]["ledger_id"],
                  "observable": lambda r: r.species, "temperature_C": lambda r: by_pot[r.pot]["t_c"],
                  "label_pH": lambda r: by_pot[r.pot]["label_ph"], "time_min": lambda r: r.t_min}
    scored_rows = [r for r, _ in scored]
    report = {
        "artifact": out_json.stem,
        "table": args.table,
        "arm": arm or "shipped",
        "ph_coupling": args.ph_coupling,
        "provenance": provenance.provenance_block(
            out_json.stem, generated_by="scripts/generators/build_submodel_probe.py",
            inputs=[fit_rows.FIT_ROWS_DIR / f"{args.table}.yml", ROOT / "src/kinetic_core/observations.py",
                    ROOT / "src/kinetic_core/ph_state.py", ROOT / "src/kinetic_core/sulfur.py",
                    ROOT / "src/kinetic_core/parameters.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
                    ROOT / "src/kinetic_core/trunk_conditions.py"]),
        "skipped": {j["pot"]: j["skip"] for j in jobs if "skip" in j},
        "rows_without_an_engine_prediction": unpredictable,
        "shape_only_series": {s: v[0]["scale_offset_dex"] for s, v in by_series.items()},
        "notes": {j["pot"]: j["notes"] for j in jobs if j["notes"]},
        "runs": {p: {"seconds": round(r["seconds"], 1), "ph_fixed_point_converged": r["converged"],
                     "sid_mmol_L": r["sid"]} for p, r in results.items()},
        "raw": fitting.diagnostics(_Residuals(scored_rows, [rec["residual"] for _, rec in scored]), None, attributes),
        "whitened": fitting.diagnostics(_Residuals(scored_rows, [rec["z"] for _, rec in scored]), None, attributes),
        "rows": records,
        "wall_seconds": round(time.time() - start, 1),
    }
    out_json.write_text(provenance.dumps_compact(report) + "\n", encoding="utf-8")
    out_md.write_text(markdown(report), encoding="utf-8")
    print(f"{len(scored)} rows scored on {len(results)} pots ({len(report['skipped'])} skipped) in {report['wall_seconds']} s "
          f"-> {data_paths.rel(out_md)}")
    return 0


def markdown(report: Dict[str, Any]) -> str:
    rows = [r for r in report["rows"] if not r["input_row"]]
    arm = "the shipped engine" if report["arm"] == "shipped" else "the shipped engine with " + " and ".join(
        "the trunk's measured pH terms applied at the lane's computed pH" if part == "trunk_ph_dynamic"
        else "glycine, the Amadori compound and the Schiff base in the charge balance" if part == "trunk_solutes"
        else f"semi-quantitative scales shared per figure (tau {part.split('tau')[1]} dex)" if part.startswith("figure_scale")
        else "the inline pH coupling" if part == "inline"
        else f"{part[5:].upper()} switched on" if part.startswith("with_")
        else f"`{part[8:]}` switched off"
        for part in re.findall(r"trunk_ph_dynamic|trunk_solutes|figure_scale_tau[0-9.inf]+|inline|with_b[0-9]+[a-z]?|without_k_[a-z_]+?(?=_without_|_inline|$)", report["arm"]))
    lines = [f"# No-refit probe of `{report['table']}`: {arm} against its verified rows", "",
             f"*Generated by `scripts/generators/build_submodel_probe.py --table {report['table']}`. Shipped constants, "
             f"pH coupling `{report.get('ph_coupling', 'fixed_point_nodes')}`, "
             "nothing fitted. Residual = predicted minus measured, in pH units for pH and as log10(predicted / measured) "
             "for amounts; whitened = residual / the row's error. t = 0 pH rows (the charge balance's input) are listed "
             "in the JSON and kept out of these tables. A semi-quantitative series is compared in shape only: its mean log "
             "residual (the unknown scale, in the JSON as `scale_offset_dex`) is removed, and a one-row series drops out.*", "",
             "## Per data set (raw residuals)", "",
             "| data set | observable | n | mean | RMS | worst |", "|---|---|---|---|---|---|"]
    groups: Dict[tuple, list] = {}
    for r in rows:
        groups.setdefault((r["data_set"], r["observable"]), []).append(r["residual"])
    for (ds, ob), v in groups.items():
        worst = max(v, key=abs)
        lines.append(f"| {ds} | {ob} | {len(v)} | {np.mean(v):+.2f} | {math.sqrt(np.mean(np.square(v))):.2f} | {worst:+.2f} |")
    lines += ["", "## Skipped pots", ""] + [f"- `{p}`: {why}" for p, why in report["skipped"].items()]
    if report.get("rows_without_an_engine_prediction"):
        lines += ["", "## Rows the lane cannot predict", ""] + [
            f"- {species}: {n}" for species, n in report["rows_without_an_engine_prediction"].items()]
    lines += ["", "## Pots run with a part of the system missing", ""] + \
             [f"- `{p}`: {'; '.join(n)}" for p, n in report["notes"].items()]
    lines += ["", fitting.diagnostics_markdown(report["raw"], "Raw residual diagnostics (pH units / dex)"),
              fitting.diagnostics_markdown(report["whitened"], "Whitened residual diagnostics (per row error)")]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
