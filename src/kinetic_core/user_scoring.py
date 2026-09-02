"""
src/kinetic_core/user_scoring.py -- BRING YOUR OWN MEASUREMENT (step 6, 2026-09-03).

A scientist with a GC-MS run has one question the panel cannot answer: "how does the
model do on MY system?" This module scores a user's measured concentrations against the
kinetic core exactly the way the panel scorecard scores a bundle (``scoring.py``: fold
error, the 3x band, the B4 reliability interval, refusals with their reason) and writes
the measurement as a BUNDLE-SHAPED record under ``results/user/`` -- never under
``data/`` -- so the next fit wave can pick it up as a declared row with its provenance
attached. Scoring is the loop's first half; refitting is a new wave by the rule in
``scripts/generators/WAVES.md``, and this module deliberately does not fit anything.

Input document (YAML or JSON), one or more systems::

    systems:
      - name: my_pea_isolate_run
        precursors: {L-Cysteine: 10.0, D-Ribose: 10.0}     # mM
        temp_C: 140.0
        time_min: 30.0
        ph: 5.0
        aw: 0.98
        matrix: water                                       # or pea_iso / soy_iso ...
        measured:
          2-methyl-3-furanthiol: {value: 198.0, unit: ppb, uncertainty_pct: 10}
          2-furfurylthiol: {value: 121.0, unit: ppb}
        quantification_class: stable_isotope_dilution_gcms  # headspace / extraction family
        source: {lab: "...", date: "2026-09-01", instrument: "HS-SPME GC-MS", notes: "..."}

Units: ``ppb`` (= ug/L in water) only, for now. A compound the core cannot represent is
returned as a refusal with the engine's reason; it is not silently dropped.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import yaml

from src import data_paths

from .engine import engine_metadata, predict
from .scoring import PASS_BAND_LEVEL, fold_error, summarise_folds

USER_RESULTS_DIR: Path = data_paths.RESULTS_ROOT / "user"

TEMPLATE = """\
# Maillard `score` spec: your own measured system(s), scored on the kinetic core.
# Precursors in mM; measured values in ppb (= ug/L in water). Nothing here is written
# under data/; the scored record lands under results/user/<name>.json.
systems:
  - name: my_ribose_cysteine_run
    precursors:
      L-Cysteine: 10.0
      D-Ribose: 10.0
    temp_C: 145.0
    time_min: 20.0
    ph: 5.0
    aw: 0.98
    matrix: water
    measured:
      2-methyl-3-furanthiol: {value: 198.0, unit: ppb, uncertainty_pct: 10}
      2-furfurylthiol: {value: 121.0, unit: ppb, uncertainty_pct: 10}
    quantification_class: stable_isotope_dilution_gcms
    source:
      lab: your lab
      date: 2026-09-01
      instrument: SIDA GC-MS after solvent extraction
      notes: phosphate buffer 0.5 M; sealed vessel
"""


class MeasurementSpecError(ValueError):
    """A user measurement document is unusable. Never repaired silently."""


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------


def load_document(path: Path | str) -> Dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    payload = json.loads(text) if str(path).endswith(".json") else yaml.safe_load(text)
    if not isinstance(payload, Mapping):
        raise MeasurementSpecError(f"{path}: must be a mapping with a `systems:` list")
    return dict(payload)


def validate_system(entry: Mapping[str, Any], *, label: str) -> Dict[str, Any]:
    from src.comparative_cli import validate_spec

    if not isinstance(entry, Mapping):
        raise MeasurementSpecError(f"{label}: must be a mapping")
    spec = validate_spec({k: v for k, v in entry.items() if k != "measured"}, label=label)
    measured = entry.get("measured")
    if not isinstance(measured, Mapping) or not measured:
        raise MeasurementSpecError(f"{label}: `measured` must be a non-empty mapping of compound -> {{value, unit}}")
    rows: Dict[str, Dict[str, Any]] = {}
    for compound, rec in measured.items():
        if isinstance(rec, (int, float)):
            rec = {"value": float(rec), "unit": "ppb"}
        if not isinstance(rec, Mapping) or rec.get("value") is None:
            raise MeasurementSpecError(f"{label}: measured[{compound!r}] needs a `value`")
        unit = str(rec.get("unit", "ppb"))
        if unit != "ppb":
            raise MeasurementSpecError(f"{label}: measured[{compound!r}] unit {unit!r}: only ppb is supported")
        value = float(rec["value"])
        if value < 0:
            raise MeasurementSpecError(f"{label}: measured[{compound!r}] is negative")
        rows[str(compound)] = {
            "value": value, "unit": unit,
            "uncertainty_pct": None if rec.get("uncertainty_pct") is None else float(rec["uncertainty_pct"]),
        }
    spec["measured"] = rows
    spec["quantification_class"] = entry.get("quantification_class")
    spec["source"] = {str(k): (v if isinstance(v, (int, float, str, bool)) or v is None else str(v)) for k, v in dict(entry.get("source") or {}).items()}
    if entry.get("matrix"):
        spec["matrix"] = str(entry["matrix"])
    return spec


def systems_of(document: Mapping[str, Any]) -> List[Dict[str, Any]]:
    systems = document.get("systems")
    if systems is None and "measured" in document:
        systems = [document]
    if not isinstance(systems, list) or not systems:
        raise MeasurementSpecError("the document needs a non-empty `systems:` list (or a single system with `measured`)")
    return [validate_system(s, label=f"systems[{i}]") for i, s in enumerate(systems)]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def score_system(spec: Mapping[str, Any], *, pass_band: float = PASS_BAND_LEVEL) -> Dict[str, Any]:
    from src.comparative_cli import spec_to_core

    core_spec = spec_to_core(spec)
    rows: List[Dict[str, Any]] = []
    refused: List[Dict[str, Any]] = []
    for compound, rec in spec["measured"].items():
        run = predict(core_spec, [compound])
        declaration = run.declaration
        if not run.answered or compound not in run.concentrations_ug_per_l:
            refused.append(
                {"compound": compound, "measured_ppb": rec["value"],
                 "reason": " | ".join(declaration.reasons) or "out of envelope",
                 "lane": declaration.lane, "envelope_state": declaration.state}
            )
            continue
        predicted = float(run.concentrations_ug_per_l[compound])
        fold = fold_error(predicted, rec["value"])
        interval = None
        inside = None
        try:
            band = run.absolutes().get(compound)
        except Exception:  # noqa: BLE001
            band = None
        if band is not None:
            interval = [float(band.lo_ug_per_l), float(band.hi_ug_per_l)]
            inside = bool(band.lo_ug_per_l <= rec["value"] <= band.hi_ug_per_l)
        rows.append(
            {
                "compound": compound,
                "measured_ppb": rec["value"],
                "uncertainty_pct": rec["uncertainty_pct"],
                "predicted_ppb": predicted,
                "fold_error": fold,
                "abs_log10_error": None if fold is None else math.log10(fold),
                "within_band": None if fold is None else bool(fold <= pass_band),
                "interval_ug_per_L": interval,
                "measured_within_interval": inside,
                "lane": declaration.lane,
                "envelope_state": declaration.state,
                "declaration_warnings": list(declaration.warnings),
            }
        )
    folds = [r["fold_error"] for r in rows]
    return {
        "name": spec["name"],
        "spec": {k: spec[k] for k in ("precursors", "temp_C", "time_min", "ph", "aw") if k in spec}
        | ({"matrix": spec["matrix"]} if spec.get("matrix") else {}),
        "quantification_class": spec.get("quantification_class"),
        "source": spec.get("source") or {},
        "rows": rows,
        "refused": refused,
        "within_band": sum(1 for r in rows if r["within_band"]),
        "scored": len(rows),
        "fold_summary": summarise_folds(folds),
        "pass_band_level": float(pass_band),
    }


def score_document(document: Mapping[str, Any], *, pass_band: float = PASS_BAND_LEVEL) -> Dict[str, Any]:
    systems = [score_system(s, pass_band=pass_band) for s in systems_of(document)]
    all_rows = [r for s in systems for r in s["rows"]]
    return {
        "artifact": "user_measurements_scored",
        "generated_on": time.strftime("%Y-%m-%d"),
        "engine": engine_metadata(),
        "pass_band_level": float(pass_band),
        "systems": systems,
        "summary": {
            "systems": len(systems),
            "scored_rows": len(all_rows),
            "refused_rows": sum(len(s["refused"]) for s in systems),
            "within_band": sum(1 for r in all_rows if r["within_band"]),
            "inside_interval": sum(1 for r in all_rows if r["measured_within_interval"]),
            "fold_summary": summarise_folds([r["fold_error"] for r in all_rows]),
        },
        "what_this_is_not": (
            "Scoring, not calibration. Nothing was refitted. To let a future fit wave read these "
            "rows, the bundle-shaped records under results/user/ carry the provenance a wave needs; "
            "the wave itself is a new pre-registered run (scripts/generators/WAVES.md)."
        ),
    }


# ---------------------------------------------------------------------------
# The bundle-shaped record (a candidate row for the next wave, never under data/)
# ---------------------------------------------------------------------------


def bundle_record(system: Mapping[str, Any]) -> Dict[str, Any]:
    spec = system["spec"]
    return {
        "benchmark_id": f"user_{system['name']}",
        "evidence_class": "user_measurement",
        "protein_type": "free" if str(spec.get("matrix", "water")) in ("water", "free") else str(spec["matrix"]),
        "precursors": {k: {"concentration_mM": float(v)} for k, v in spec["precursors"].items()},
        "conditions": {"temp_C": spec["temp_C"], "time_min": spec["time_min"], "ph": spec["ph"], "water_activity": spec.get("aw")},
        "measured_volatiles": {
            r["compound"]: {"conc_ppb": r["measured_ppb"], "uncertainty_pct": r["uncertainty_pct"], "value_status": "user_reported"}
            for r in list(system["rows"]) + [
                {"compound": x["compound"], "measured_ppb": x["measured_ppb"], "uncertainty_pct": None} for x in system["refused"]
            ]
        },
        "content_verification": {"quantification_class": system.get("quantification_class")},
        "source_metadata": {"origin": "user_measurement", **(system.get("source") or {})},
        "scored_on": {
            "engine": None,  # filled by write_records
            "within_band": system["within_band"], "scored": system["scored"],
            "refused": [x["compound"] for x in system["refused"]],
        },
        "next_wave_note": (
            "A fit wave that reads this record must list its benchmark_id in that wave's "
            "fit_target_ids (scripts/generators/generate_core_fit_targets.py) so the scorecard "
            "flags it in_core_fit."
        ),
    }


def write_records(payload: Mapping[str, Any], out_dir: Path | str = USER_RESULTS_DIR) -> List[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for system in payload["systems"]:
        record = bundle_record(system)
        record["scored_on"]["engine"] = payload["engine"]
        path = out_dir / f"{record['benchmark_id']}.json"
        path.write_text(json.dumps(record, indent=2, default=str) + "\n", encoding="utf-8")
        written.append(path)
    scored = out_dir / "scored.json"
    scored.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    written.append(scored)
    return written


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt(v: Optional[float]) -> str:
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "yes" if v else "no"
    return f"{v:.3g}"


def render_text(payload: Mapping[str, Any]) -> str:
    out: List[str] = []
    s = payload["summary"]
    out.append("=" * 96)
    out.append("  SCORE [KINETIC CORE]   your measurements against the model")
    out.append("=" * 96)
    out.append("")
    out.append(
        f"  {s['scored_rows']} row(s) scored, {s['refused_rows']} refused; within {payload['pass_band_level']:.0f}x: "
        f"{s['within_band']}/{s['scored_rows']}; inside the reliability interval: {s['inside_interval']}/{s['scored_rows']}; "
        f"median fold {_fmt(s['fold_summary']['median_fold_error'])}"
    )
    out.append("")
    for system in payload["systems"]:
        out.append(f"  {system['name']}  ({system['spec']['temp_C']:g} C, {system['spec']['time_min']:g} min, pH {system['spec']['ph']:g})")
        out.append(f"  {'compound':<34} {'measured':>10} {'predicted':>10} {'fold':>8}  within 3x  in interval  lane")
        out.append("  " + "-" * 92)
        for r in system["rows"]:
            out.append(
                f"  {r['compound']:<34} {_fmt(r['measured_ppb']):>10} {_fmt(r['predicted_ppb']):>10} "
                f"{_fmt(r['fold_error']):>8}  {_fmt(r['within_band']):<9}  {_fmt(r['measured_within_interval']):<11}  {r['lane']}"
            )
        for x in system["refused"]:
            out.append(f"  {x['compound']:<34} {_fmt(x['measured_ppb']):>10} {'REFUSED':>10}  {x['reason'][:70]}")
        out.append("")
    out.append("  " + payload["what_this_is_not"])
    out.append("")
    return "\n".join(out)


__all__ = [
    "MeasurementSpecError", "TEMPLATE", "USER_RESULTS_DIR", "bundle_record", "load_document",
    "render_text", "score_document", "score_system", "systems_of", "validate_system", "write_records",
]
