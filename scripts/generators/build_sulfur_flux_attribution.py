#!/usr/bin/env python
"""
Where the sulfur lane's carbon and sulfur actually go (2026-10-28).

For eight reference pots (the engine path, the shipped calibration), integrate every reaction's rate over the
run and rank the producers and consumers of the trunk and thiol species. Each pot is re-integrated on a
log-dense output grid with the exact arguments the engine passed, and the rates use the run's own node
constants (``sulfur.reaction_flux_integrals``), so the budget is the one the trajectory saw; its closure
(initial - final against net integrated production) is checked and recorded, and the generator refuses to
write a budget that does not close to 1 %.

Also recorded per pot: the time by which half the charged cysteine and pentose are gone, and the in-situ pH
at 0.1 % and at the end of the run -- the two numbers that showed the trunk runs in seconds.

EVIDENCE GROUNDING (2026-10-28, owner: "scientifically grounded"; WAVES.md "Evidence decides, scores
diagnose"): the share of the pot's integrated reaction flux that runs through steps whose rate constant is a
measurement, a fit, or neither. Each step is graded by its parameter's ``evidence_class``. A step with no
parameter (a rate derived from others, e.g. a detailed-balance reverse) is graded "derived". This measures how
much of the answer rests on evidence, not whether the answer is right.

    python scripts/generators/build_sulfur_flux_attribution.py
    # writes results/validation/sulfur_flux_attribution.{json,md}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths, provenance  # noqa: E402
from src.kinetic_core import engine, sulfur  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402
from src.kinetic_core.species_sulfur import SULFUR_INDEX  # noqa: E402

OUT = data_paths.VALIDATION_DIR / "sulfur_flux_attribution.json"
OUT_MD = data_paths.VALIDATION_DIR / "sulfur_flux_attribution.md"
TRACK = ("Cys", "PENT", "TTCA", "DPO", "TDP", "DDP", "NF", "FUR", "H2S", "MFT", "FFT")
CLOSURE_SPECIES = ("Cys", "PENT", "H2S", "MFT", "FFT")
CLOSURE_TOL = 0.01
TOP = 5


def _phosphate(molar: float, source: str) -> BufferSpec:
    return BufferSpec(kind="phosphate", phosphate_mol_l=molar, declared=True, source=source)


#: name -> (precursors mmol/L, T C, minutes, label pH, buffer, what the corpus measured in that pot)
POTS = {
    "hofmann1998_ribose_cys_ph5_145C_20min": (
        {"L-Cysteine": 33.0, "D-Ribose": 100.0}, 145.0, 20.0, 5.0, _phosphate(0.5, "Hofmann 1998"),
        "thiol yields only (hofmann1998_extraction.md); no cysteine or sugar time course"),
    "hofmann1998_ribose_cys_ph3_145C_20min": (
        {"L-Cysteine": 33.0, "D-Ribose": 100.0}, 145.0, 20.0, 3.0, _phosphate(0.5, "Hofmann 1998"),
        "MFT 553 / FFT 229 ppb (the pH 3 arm of the sweep; MFT 553 -> 198 -> 25 and FFT 229 -> 121 -> 12 ppb at pH 3 / 5 / 7)"),
    "hofmann1998_ribose_cys_ph7_145C_20min": (
        {"L-Cysteine": 33.0, "D-Ribose": 100.0}, 145.0, 20.0, 7.0, _phosphate(0.5, "Hofmann 1998"),
        "MFT 25 / FFT 12 ppb (the pH 7 arm of the sweep)"),
    "cao2017_xylose_cys_90C_300min": (
        {"L-Cysteine": 400.0, "D-Xylose": 400.0}, 90.0, 300.0, 5.5, _phosphate(0.2, "Cao 2017"),
        "free Cys + TTCA + Cys-Amadori = 0.425 M of 0.40 at 5 h: no measurable cysteine destroyed; pH 5.49 -> 4.67 (0.5 h) -> 4.15 (cao2017_extraction.md Figs. 4-5)"),
    "zhang2025b_xylose_cys_ph65_120C_50min": (
        {"L-Cysteine": 40.0, "D-Xylose": 40.0}, 120.0, 50.0, 6.5, _phosphate(0.2, "Zhang 2025b"),
        "free Cys 1.0 + TTCA 18.2 + Cys-Amadori 1.6 = 20.8 of 40 mM at 50 min: ~48 % destroyed, xylose likewise (zhang2025b_extraction.md sec. 3.4)"),
    "zhai2019_xylose_cys_unbuffered_90C_160min": (
        {"L-Cysteine": 82.6, "D-Xylose": 165.2}, 90.0, 160.0, 7.4, BufferSpec(kind="none", declared=True, source="Zhai 2019, water + NaOH"),
        "TTCA plateau 54-63 mM (65-77 % of the cysteine) at 70-90 C (zhai2019_extraction.md)"),
    "yiltirak_ribose_cys_120C_60min": (
        {"L-Cysteine": 25.0, "D-Ribose": 25.0}, 120.0, 60.0, 5.5, _phosphate(0.5, "reference pot"),
        "thiol yields only"),
    "food_like_unbuffered_ribose_cys_140C_30min": (
        {"L-Cysteine": 10.0, "D-Ribose": 10.0}, 140.0, 30.0, 6.0, BufferSpec(kind="none", declared=True, source="reference pot"),
        "a constructed reference, no measurement"),
}


def _engine_arguments(name: str, pot) -> Tuple:
    """The (parameters, temperature, initial state, keyword arguments) the engine hands the sulfur lane."""
    prec, t_c, minutes, ph, buf, _ = pot
    captured: Dict = {}
    original = engine.integrate_sulfur

    def _capture(parameters, temperature_k, state, times, **kw):
        captured["args"] = (parameters, temperature_k, dict(state), dict(kw))
        return original(parameters, temperature_k, state, times, **kw)

    engine.integrate_sulfur = _capture
    try:
        spec = engine.FormulationSpec(name=name, precursors=prec, process=engine.ProcessSpec(
            thermal=engine.ThermalProgram.isothermal(t_c, minutes), ph=ph, buffer=buf))
        decl = engine.declare_envelope(spec, ["2-methyl-3-furanthiol"])
        engine._integrate_program(engine.SULFUR, engine.core_parameters(engine.SULFUR),
                                  dict(decl.mapped_precursors), spec.process)
    finally:
        engine.integrate_sulfur = original
    return captured["args"]


#: evidence_class -> the grounding bucket it counts toward
GROUNDING = {"measured_rate": "measured", "measured_activation_energy": "measured",
             "derived_from_fit_data": "fitted", "bounded_from_a_timescale_bracket": "bounded",
             "structural_constant": "measured"}


def _grade(rxn, parameters) -> str:
    if rxn.parameter_key is None or rxn.parameter_key not in parameters:
        return "derived"
    return GROUNDING.get(getattr(parameters[rxn.parameter_key], "evidence_class", ""), "unclassified")


def _grounding(flux: Dict[str, float], parameters) -> Dict[str, float]:
    """Share of the pot's total integrated reaction flux per evidence grade."""
    totals: Dict[str, float] = {}
    for rxn in sulfur.FULL_REACTIONS:
        f = flux.get(rxn.key, 0.0)
        if f > 0:
            g = _grade(rxn, parameters)
            totals[g] = totals.get(g, 0.0) + f
    whole = sum(totals.values())
    return {g: v / whole for g, v in sorted(totals.items())} if whole > 0 else {}


def _half_time(run, key: str):
    series = run.series(key)
    if series[0] <= 0:
        return None
    below = np.nonzero(series <= 0.5 * series[0])[0]
    return float(run.times_min[below[0]]) if below.size else None


def _ranked(flux: Dict[str, float], key: str) -> Dict[str, List]:
    made, used = [], []
    for rxn in sulfur.FULL_REACTIONS:
        nu = rxn.products.get(key, 0) - rxn.reactants.get(key, 0)
        f = flux[rxn.key]
        if nu > 0 and f > 0:
            made.append((rxn.key, nu * f))
        elif nu < 0 and f > 0:
            used.append((rxn.key, -nu * f))
    out = {}
    for label, rows in (("made", made), ("used", used)):
        total = sum(v for _, v in rows)
        rows.sort(key=lambda r: -r[1])
        out[label] = {"total_mmol_L": total,
                      "top": [{"reaction": k, "mmol_L": v, "share": v / total if total > 0 else 0.0} for k, v in rows[:TOP]]}
    return out


def main() -> int:
    pots = {}
    for name, pot in POTS.items():
        parameters, temperature_k, state, kw = _engine_arguments(name, pot)
        minutes = pot[2]
        grid = np.unique(np.concatenate([np.linspace(0.0, minutes, 201), np.geomspace(minutes * 1e-7, minutes, 401)]))
        run = sulfur.integrate_sulfur(parameters, temperature_k, state, grid, **kw)
        flux = sulfur.reaction_flux_integrals(run)
        closure = sulfur.flux_closure(run, flux, CLOSURE_SPECIES)
        if max(closure.values()) > CLOSURE_TOL:
            raise SystemExit(f"{name}: flux budget does not close: {closure}")
        ph_series = run.metadata.get("ph_series")
        early = int(np.searchsorted(run.times_min, 1e-3 * minutes))
        pots[name] = {
            "measured_in_the_corpus": pot[5],
            "temperature_C": pot[1], "minutes": minutes, "label_ph": pot[3],
            "charged_mmol_L": {k: float(run.concentrations[0, SULFUR_INDEX[k]]) for k in ("Cys", "PENT")},
            "final_mmol_L": {k: run.final(k) for k in TRACK},
            "half_time_min": {k: _half_time(run, k) for k in ("Cys", "PENT")},
            "ph_in_situ": {"start": float(ph_series[0]), "at_0.1_percent_of_run": float(ph_series[early]),
                           "end": float(ph_series[-1])} if ph_series else None,
            "closure_relative_error": closure,
            "fixed_point_residual_pH": run.metadata.get("ph_fixed_point_residual"),
            "species": {k: _ranked(flux, k) for k in TRACK},
            "evidence_grounding": _grounding(flux, parameters),
            "grade_of_reaction": {rxn.key: _grade(rxn, parameters) for rxn in sulfur.FULL_REACTIONS if flux.get(rxn.key, 0.0) > 0},
        }
        print(f"  {name}: Cys t1/2 {pots[name]['half_time_min']['Cys']}, closure {max(closure.values()):.1e}", flush=True)

    payload = {
        "artifact": "sulfur_flux_attribution",
        "provenance": provenance.provenance_block(
            "sulfur_flux_attribution", generated_by="scripts/generators/build_sulfur_flux_attribution.py",
            inputs=[engine._B2_FIT_REPORT, ROOT / "src/kinetic_core/sulfur.py"]),
        "calibration": data_paths.rel(engine._B2_FIT_REPORT),
        "ph_numerics": {"node_grid": sulfur.PH_NODE_GRID, "fixed_point": sulfur.PH_FIXED_POINT},
        "pots": pots,
    }
    OUT.write_text(json.dumps(payload, indent=2))

    lines = ["# Sulfur-lane flux attribution", "",
             f"*Generated by `scripts/generators/build_sulfur_flux_attribution.py` on the shipped calibration "
             f"(`{data_paths.rel(engine._B2_FIT_REPORT)}`). Every budget closes to "
             f"{CLOSURE_TOL:.0%} (checked). Shares are of each species' total production or consumption.*", ""]
    for name, p in pots.items():
        ph = p["ph_in_situ"]
        lines += [f"## {name}", "",
                  f"Measured: {p['measured_in_the_corpus']}.", "",
                  f"Half-life of the charged cysteine / pentose: {p['half_time_min']['Cys']} / "
                  f"{p['half_time_min']['PENT']} min of {p['minutes']:g}; in-situ pH "
                  + (f"{ph['start']:.2f} -> {ph['at_0.1_percent_of_run']:.2f} (0.1 % of the run) -> {ph['end']:.2f}." if ph else "clamped."),
                  "", "Evidence grounding of the flux: " + ", ".join(f"{g} {v:.0%}" for g, v in p["evidence_grounding"].items())
                  + " (step grades in brackets below: m measured, f fitted, b bounded, d derived).",
                  "", "| species | final mM | made by (share) | used by (share) |", "| --- | --- | --- | --- |"]
        for key in TRACK:
            s = p["species"][key]
            if s["made"]["total_mmol_L"] + s["used"]["total_mmol_L"] <= 1e-12:
                continue
            grade = p["grade_of_reaction"]
            fmt = lambda side: ", ".join(f"{r['reaction']} [{grade.get(r['reaction'], '?')[0]}] {r['share']:.0%}"  # noqa: E731
                                         for r in s[side]["top"][:3]) or "—"
            lines.append(f"| {key} | {p['final_mmol_L'][key]:.3g} | {fmt('made')} | {fmt('used')} |")
        lines.append("")
    OUT_MD.write_text("\n".join(lines))
    print(f"wrote {OUT} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
