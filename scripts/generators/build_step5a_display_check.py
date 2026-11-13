#!/usr/bin/env python
"""
The step-5 display rule as a diagnostic on the engine as shipped (STATUS_AND_PLAN sec. 0, step 5a).
Rules: docs/reference/lever_table/STEP5A_PREREG.md; a rule here that differs from that file is a bug here.

For every lever x odorant cell with TEST contrasts from two or more laboratories (after the sources the engine has
already seen are removed), three predictions of each TEST contrast are compared on the log-ratio per lever step:
the engine, "no change" (zero), and the DEVELOPMENT-only lever table's estimate. The engine would be displayed in a
cell where its error is below both, the 90 % interval of each difference (bootstrap over laboratories) excluding zero.

Two modes, run in this order:

* ``--items``: the cell list, the exclusions, and for every TEST contrast whether it maps to an engine run, from the
  keyed conditions and ``engine.declare_envelope`` alone. No integration. Writes ``step5a_display_check_items``.
* ``--score``: integrates the mapped arms (``engine.compare``, the shipped path of a ratio, refusals included) and
  writes ``step5a_display_check``.

No constant or structure of the engine is changed; nothing ships on the result. This generator never opens the
all-data lever table.

    ./scripts/docker_maillard.sh run "python scripts/generators/build_step5a_display_check.py --items"
    ./scripts/docker_maillard.sh run "python scripts/generators/build_step5a_display_check.py --score"
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "scripts" / "generators"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import build_lever_table as LT  # noqa: E402
from src import artifact_io, compound_keys, data_paths, ingredients, provenance  # noqa: E402
from src.kinetic_core import engine  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

ITEMS_OUT = data_paths.VALIDATION_DIR / "step5a_display_check_items.json"
SCORE_OUT = data_paths.VALIDATION_DIR / "step5a_display_check.json"
PREREG = data_paths.LEVER_TABLE_DIR / "STEP5A_PREREG.md"
RIVAL = LT.OUT_DEV      # the DEVELOPMENT-only table; the all-data table is never opened here

#: Sources removed from the TEST evidence beyond those flagged ``engine_seen`` in new_labs.yml (prereg sec. 2).
ENGINE_FITTED_OR_SCORED = {
    "kumazawa2003": "fitted: its pH grid is a fit row of the thiol sub-model (data/lit/fit_rows/thiol_submodel.yml; "
                    "parameters_sulfur 'kumazawa2003_ph_grid')",
}
#: Prereg sec. 6: the fit estimates of the rival table as used, hashed (cell key, estimates), so that additive
#: fields in the table's file (PR 130) cannot move the rival unnoticed.
RIVAL_FIT_DIGEST = "e5a8dfe97195ec2d595d8c69e7aec5abef09938716e3c9018d86d4b5b1a13b3f"

MIN_LABS = 2
NOISE_DEX = math.log10(1.4)     # the directional panel's measurement-noise fold (build_contrast_scores.NOISE_DEX)
FLOOR_DEX = 1e-3                # an engine change below this over the whole contrast is numerical, read as zero
SEED = 20261112
N_BOOT = 2000
READINGS = ("front_door", "added_precursors", "added_precursors_ph_clamped")
PRIMARY = "front_door"

#: g/mol of the named reagent (free compound; thiamine as its hydrochloride, the reagent form).
MOLAR_MASS = {"glucose": 180.16, "fructose": 180.16, "ribose": 150.13, "xylose": 150.13, "arabinose": 150.13,
              "cysteine": 121.16, "methionine": 149.21, "glycine": 75.07, "proline": 115.13, "glutamine": 146.15,
              "lysine": 146.19, "alanine": 89.09, "asparagine": 132.12, "thiamine": 337.27, "glutathione": 307.32}
#: Prereg sec. 4.2: a keyed precursor that is no engine precursor is the MATRIX MATERIAL when its name holds one of
#: these words; any other unmapped name is an unmapped pure compound.
MATRIX_WORDS = ("hydrolysate", "protein", "extract", "powder", "meat", "trimmings", "carcass", "solids",
                "free amino acids")


#: Where the compound registry's spelling match would mislead: 'Pyrazines' is a class total, not the engine's
#: unsubstituted pyrazine (which the registry files under that spelling).
TARGET_OVERRIDES: Dict[str, Optional[str]] = {"Pyrazines": None, "Unsubstituted pyrazine": "pyrazine"}
#: The keyed matrices that are a solution of pure compounds (template: buffer | water | meat | plant protein | ...).
AQUEOUS_MATRICES = ("buffer", "water")


# ----------------------------------------------------------------------------------------------- names and amounts
def norm_name(key: Any) -> str:
    """A keyed precursor name as the engine's vocabulary spells it: lower case, no parenthetical, no L-/D- prefix."""
    t = re.sub(r"\([^)]*\)", "", str(key)).strip().lower()
    for pre in ("l-", "d-"):
        if t.startswith(pre):
            t = t[len(pre):]
    return t.strip()


def engine_precursor(key: Any) -> Optional[str]:
    """The front-door precursor name of a keyed precursor, or None."""
    name = norm_name(key)
    return name if name in engine.PRECURSOR_ALIASES else None


def target_alias(display_name: str) -> Optional[str]:
    """The engine target name of a lever-table odorant: the first target alias (dict order) that the compound
    registry resolves to the same compound."""
    if display_name in TARGET_OVERRIDES:
        return TARGET_OVERRIDES[display_name]
    for alias in engine.TARGET_ALIASES:
        if compound_keys.same_compound(alias, str(display_name)):
            return alias
    return None


def is_matrix_material(key: Any) -> bool:
    name = str(key).lower()
    return engine_precursor(key) is None and any(w in name for w in MATRIX_WORDS)


_NUM = r"(\d+(?:\.\d+)?)"


def parse_amount(text: Any) -> Optional[Tuple[str, float]]:
    """A printed amount as (kind, number): ('molar', mM), ('g_per_l', g/L), ('mass', g), ('per_protein', ratio).
    None where the text gives no amount. First rule that matches (prereg sec. 4.3)."""
    s = str(text).strip()
    m = re.match(rf"{_NUM}\s*mmol in {_NUM}\s*mL", s)
    if m:
        return "molar", 1000.0 * float(m.group(1)) / float(m.group(2))
    m = re.match(rf"{_NUM}\s*g\s*/\s*100\s*(?:mL|g)\b", s) or re.match(rf"{_NUM}\s*%", s)
    if m:
        return "g_per_l", 10.0 * float(m.group(1))
    m = re.match(rf"{_NUM}\s*(mM|mmol/L|M|mol/L)\b", s)
    if m:
        return "molar", float(m.group(1)) * (1.0 if m.group(2) in ("mM", "mmol/L") else 1000.0)
    m = re.match(rf"{_NUM}\s*g\b", s)
    if m:
        return "mass", float(m.group(1))
    m = re.match(rf"protein-to-sugar 1:{_NUM}", s)
    if m:
        return "per_protein", float(m.group(1))
    return None


def pot_total_g(precursors: Mapping[str, Any], buffer_text: str) -> Optional[float]:
    """The pot's mass in g (1 g/mL) behind plain masses, first rule that applies (prereg sec. 4.3): a printed total;
    the masses over a printed solids share; a printed water volume plus the masses; else the sum of the masses."""
    masses = [p[1] for k, v in precursors.items()
              for p in [parse_amount(v)] if p and p[0] == "mass" and norm_name(k) != "water"]
    for v in precursors.values():
        m = re.search(rf"total {_NUM}\s*g", str(v))
        if m:
            return float(m.group(1))
    share = None
    for k, v in precursors.items():
        if "total solids" in str(k).lower():
            m = re.match(rf"{_NUM}\s*%", str(v).strip())
            share = float(m.group(1)) / 100.0 if m else None
    m = re.search(rf"{_NUM}\s*%\s*w/v", buffer_text)
    if share is None and m:
        share = float(m.group(1)) / 100.0
    if share:
        return sum(masses) / share if masses else None
    m = re.search(rf"{_NUM}\s*mL", buffer_text)
    if m:
        return float(m.group(1)) + sum(masses)
    water = [p[1] for k, v in precursors.items() for p in [parse_amount(v)]
             if p and p[0] == "mass" and norm_name(k) == "water"]
    return (sum(masses) + sum(water)) or None


def protein_g_per_l(precursors: Mapping[str, Any]) -> Optional[float]:
    for k, v in precursors.items():
        m = re.match(rf"{_NUM}\s*%\s*w/w protein", str(v).strip())
        if m and is_matrix_material(k):
            return 10.0 * float(m.group(1))
    return None


# ------------------------------------------------------------------------------------- one keyed arm -> one engine run
def _buffer(text: str, clamped: bool) -> Optional[BufferSpec]:
    if clamped:
        return BufferSpec(kind="clamped", declared=True, source="step 5a sensitivity: pH held at the keyed value")
    t = text.strip().lower()
    m = re.match(rf"{_NUM}\s*m phosphate", t)
    if m:
        return BufferSpec(kind="phosphate", phosphate_mol_l=float(m.group(1)), declared=True, source=text)
    if t.startswith("none"):
        return BufferSpec(kind="none", declared=True, source=f"keyed: {text}")
    return None     # the engine's own default (unbuffered, with its warning)


def arm_spec(arm: Mapping[str, Any], lever_class: str, reading: str, name: str = "arm"
             ) -> Tuple[Optional[engine.FormulationSpec], str]:
    """(the engine run of one keyed arm, '') or (None, why not). Prereg sec. 4."""
    cond = arm.get("conditions") or {}
    raw = LT._num(arm.get("lever_value"))
    ph, t_c, minutes = (LT._num(cond.get(k)) for k in ("ph", "temp_C", "time_min"))
    if lever_class == "pH" and ph is None:
        ph = raw
    if lever_class == "temperature" and t_c is None:
        t_c = raw
    if lever_class == "time" and minutes is None:
        minutes = raw
    basis = str(cond.get("ph_basis") or "unstated")
    if ph is None or basis == "final":
        return None, "pH not printed" if ph is None else "only the final pH is printed"
    if t_c is None or minutes is None:
        return None, "temperature or time not printed"
    keyed = cond.get("precursors") or {}
    if not isinstance(keyed, Mapping) or not keyed:
        return None, "no precursor charge keyed"
    matrix = str(cond.get("matrix") or "").strip().lower()
    if reading == "front_door" and not matrix.startswith(AQUEOUS_MATRICES):
        return None, f"matrix material not on file: {matrix or 'unstated'}"
    buffer_text = str(cond.get("buffer") or "")
    total = pot_total_g(keyed, buffer_text)
    charges: Dict[str, float] = {}
    for key, value in keyed.items():
        name_n = norm_name(key)
        if name_n == "water":
            continue
        amount = parse_amount(value)
        ingredient = ingredients._key(name_n) in ingredients.table()
        if engine_precursor(key) is None and not ingredient:
            if is_matrix_material(key):
                if reading == "front_door":
                    return None, f"matrix material not on file: {norm_name(key)}"
                continue        # 'added precursors' readings: the material is left out, and the output says so
            return None, f"unmapped precursor: {name_n}"
        if amount is None:
            return None, f"charge not derivable: {name_n} ({str(value)[:40]})"
        kind, x = amount
        if kind == "mass":
            if not total:
                return None, f"charge not derivable: {name_n} (a mass with no pot size)"
            g_per_l = 1000.0 * x / total
        elif kind == "per_protein":
            protein = protein_g_per_l(keyed)
            if protein is None:
                return None, f"charge not derivable: {name_n} (a ratio to protein with no protein content)"
            g_per_l = x * protein
        elif kind == "g_per_l":
            g_per_l = x
        else:
            g_per_l = None
        if ingredient:
            if g_per_l is None:
                return None, f"charge not derivable: {name_n} (an ingredient needs a mass share)"
            if g_per_l > 0:
                expanded, _ = ingredients.expand({"ingredients": {name_n: g_per_l / 10.0}})
                for k2, v2 in expanded["precursors"].items():
                    charges[k2] = charges.get(k2, 0.0) + float(v2)
            continue
        if g_per_l is None:
            mm = x
        else:
            if name_n not in MOLAR_MASS:
                return None, f"charge not derivable: {name_n} (no molar mass declared)"
            mm = 1000.0 * g_per_l / MOLAR_MASS[name_n]
        if mm > 0:
            charges[name_n] = charges.get(name_n, 0.0) + mm
    if not charges:
        return None, "no engine precursor charged"
    process = engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(t_c, minutes), ph=ph,
                                 buffer=_buffer(buffer_text, reading.endswith("ph_clamped") or basis == "held"),
                                 vessel=None)
    return engine.FormulationSpec(name=name, precursors=charges, process=process), ""


def spec_key(spec: engine.FormulationSpec) -> str:
    p = spec.process
    species: Dict[str, float] = defaultdict(float)      # the engine's species, so that two names of one species agree
    for k, v in spec.precursors.items():
        species[engine.PRECURSOR_ALIASES.get(k, k)] += float(v)
    return json.dumps([sorted((k, round(v, 9)) for k, v in species.items()),
                       [list(s) for s in p.thermal.segments], p.ph, p.buffer.as_dict() if p.buffer else None],
                      sort_keys=True, default=str)


# ------------------------------------------------------------------------------------------------- the TEST items
def excluded_sources() -> Dict[str, str]:
    new = yaml.safe_load(data_paths.LEVER_TABLE_NEW_LABS.read_text()) or {}
    out = {str(k).lower(): "flagged engine_seen in new_labs.yml: " + str(v.get("note", "")).strip()
           for k, v in new.items() if isinstance(v, Mapping) and v.get("engine_seen")}
    out.update(ENGINE_FITTED_OR_SCORED)
    return out


def run_points(c: Mapping[str, Any], arms: Sequence[Mapping[str, Any]]) -> List[Tuple[float, float, float, int]]:
    """The arms ``LT.numeric_observation`` fits its slope on: (x, log10 y, variance, arm index), sorted by x, the
    longest monotone run where the ladder is not monotone. A copy of that function's selection; ``build_items``
    asserts that the slope it gives equals the lever table's."""
    lc = c["lever_class"]
    pts = []
    for i, a in enumerate(arms):
        x = LT._lever_x(lc, a)
        ly, var, _ = LT.arm_log_and_var(a, bool(c.get("digitised")))
        if x is not None and ly is not None:
            pts.append((x, ly, var, LT._raw_lever(lc, a), i))
    pts.sort()
    diffs = [np.sign(b[1] - a[1]) for a, b in zip(pts, pts[1:])]
    if not (all(d >= 0 for d in diffs) or all(d <= 0 for d in diffs)):
        best: List[Any] = []
        for i in range(len(pts)):
            for j in range(i + 2, len(pts) + 1):
                seg = pts[i:j]
                ds = [np.sign(b[1] - a[1]) for a, b in zip(seg, seg[1:])]
                if (all(d >= 0 for d in ds) or all(d <= 0 for d in ds)) and len(seg) > len(best):
                    best = seg
        pts = best
    return [(p[0], p[1], p[2], p[4]) for p in pts]


def wls_slope(x: Sequence[float], y: Sequence[float], var: Sequence[float]) -> float:
    x, y, w = np.asarray(x, float), np.asarray(y, float), 1.0 / np.asarray(var, float)
    xb = float(np.sum(w * x) / np.sum(w))
    return float(np.sum(w * (x - xb) * (y - np.sum(w * y) / np.sum(w))) / np.sum(w * (x - xb) ** 2))


def build_items(sources: Sequence[Mapping[str, Any]], labs: Mapping[str, Tuple[str, str]],
                excluded: Mapping[str, str]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """One item per TEST contrast of a numeric cell and per TEST pair of a categorical cell, for the cells with two
    or more TEST laboratories after the exclusions. Each item carries its keyed arms; no engine output."""
    obs = []
    for doc in sources:         # one contrast at a time, so that each observation keeps its keyed contrast
        for c in doc["contrasts"] or []:     # (ids repeat across a source's files)
            for o in LT.observations([{"source": doc["source"], "contrasts": [c]}], labs)[0]:
                obs.append({**o, "keyed": c})
    test = [o for o in obs if o["partition"] == "test"]
    kept = [o for o in test if o["source"].lower() not in excluded]
    cells = defaultdict(list)
    for o in kept:
        cells[(o["compound"], o["lever_class"], o["band"])].append(o)
    qualifying = {k: v for k, v in cells.items() if len({o["lab"] for o in v}) >= MIN_LABS}
    items = []
    for (compound, lc, band), rows in sorted(qualifying.items()):
        for o in rows:
            c = o["keyed"]
            arms = [dict(a) for a in (c.get("arms") or [])]
            LT._trace_fill(arms)
            base = {"cell": [compound, lc, band], "compound": compound, "lever_class": lc, "band": band,
                    "lab": o["lab"], "source": o["source"], "contrast": o["contrast"], "kind": o["kind"],
                    "censored_arms": int(o.get("censored_arms") or 0), "peak_area": bool(o["peak_area"])}
            if o["kind"] == "numeric":
                pts = run_points(c, arms)
                xs, ys, vs = [p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts]
                if abs(wls_slope(xs, ys, vs) - o["y"]) > 1e-9:
                    raise SystemExit(f"{o['source']} {o['contrast']}: run_points departs from the lever table's slope")
                items.append({**base, "m": o["y"], "span_steps": max(xs) - min(xs), "x": xs, "var": vs,
                              "arms": [arms[p[3]] for p in pts]})
            else:
                usable = [(LT._level(a.get("lever_value") or a.get("label")), a) for a in arms
                          if LT.arm_log_and_var(a, bool(c.get("digitised")))[0] is not None]
                if lc == "sugar identity":
                    usable = [u for u in usable if u[0] in LT.SUGAR_LEVELS]
                others = [u for u in usable[1:] if u[0] != usable[0][0]]      # the order of the table's pairs
                assert [u[0] for u in others] == [p["level"] for p in o["pairs"]], (o["source"], o["contrast"])
                for p, (_, arm) in zip(o["pairs"], others):
                    items.append({**base, "m": p["y"], "span_steps": 1.0, "level": p["level"], "against": p["against"],
                                  "arms": [usable[0][1], arm]})
    meta = {"test_observations": len(test), "test_observations_after_exclusions": len(kept),
            "test_cells": len(cells), "qualifying_cells": len(qualifying),
            "excluded_observations": dict(Counter(o["source"] for o in test if o["source"].lower() in excluded)),
            "test_laboratories": sorted({o["lab"] for o in kept}), "test_sources": sorted({o["source"] for o in kept})}
    return items, meta


def classify(item: Mapping[str, Any], reading: str) -> Tuple[str, Optional[List[engine.FormulationSpec]]]:
    """('simulable', the arms' engine runs) or (why not, None): the target, each arm's mapping, the envelope, and
    whether the lever is expressed at all in the engine's inputs. No integration."""
    target = target_alias(item["compound"])
    if target is None:
        return "no engine target for this odorant", None
    specs = []
    for i, arm in enumerate(item["arms"]):
        spec, why = arm_spec(arm, item["lever_class"], reading, name=f"{item['contrast']}#{i}")
        if spec is None:
            return why, None
        decl = engine.declare_envelope(spec, [target])
        if not decl.is_answerable:
            return "the engine refuses the arm: " + (decl.reasons[0] if decl.reasons else decl.state)[:160], None
        specs.append(spec)
    if len({spec_key(s) for s in specs}) < 2:
        return "the engine has no term for this lever (the arms map to identical inputs)", None
    return "simulable", specs


# ------------------------------------------------------------------------------------------------------ the rival
def rival_cells() -> Tuple[Dict[Tuple[str, str, str], Mapping[str, Any]], str, str]:
    """The DEVELOPMENT-only table's cells, the file's sha256, and the digest of the fit estimates used here."""
    raw = RIVAL.read_bytes()
    cells = {(c["compound"], c["lever_class"], c["band"]): c for c in json.loads(raw)["cells"]}
    used = []
    for key, c in sorted(cells.items()):
        f = c["fit"]
        if "step" in c:
            used.append([list(key), round(f["estimate"], 9)])
        else:
            used.append([list(key), f.get("reference"),
                         sorted((lv, round(s["estimate"], 9)) for lv, s in f["levels"].items() if "estimate" in s)])
    digest = hashlib.sha256(json.dumps(used, sort_keys=True).encode()).hexdigest()
    return cells, hashlib.sha256(raw).hexdigest(), digest


def rival_prediction(cells: Mapping[Tuple[str, str, str], Mapping[str, Any]], item: Mapping[str, Any]
                     ) -> Optional[float]:
    """The table's prediction of the item's log ratio per step, or None where it has no estimate (prereg sec. 5)."""
    c = cells.get(tuple(item["cell"]))
    if c is None:
        return None
    f = c["fit"]
    if item["kind"] == "numeric":
        return float(f["estimate"])

    def effect(level: str) -> Optional[float]:
        if level == f.get("reference"):
            return 0.0
        s = (f.get("levels") or {}).get(level) or {}
        return float(s["estimate"]) if "estimate" in s else None

    a, b = effect(item["level"]), effect(item["against"])
    return None if a is None or b is None else a - b


# ----------------------------------------------------------------------------------------------------- the engine
class Runner:
    """``engine.compare`` on two arms with the arms' integrations cached (one integration per distinct run)."""

    def __init__(self) -> None:
        self.cache: Dict[Tuple[str, str], Any] = {}
        self.integrations = 0

    def _predict(self, spec, targets):
        key = (spec_key(spec), ",".join(targets))
        if key not in self.cache:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.cache[key] = engine.predict(spec, targets)
            self.integrations += 1
        return self.cache[key]

    def log_ratio(self, spec_a, spec_b, target: str) -> Tuple[Optional[float], str]:
        """log10(a / b) of ``target`` as the shipped engine answers it, or (None, the refusal)."""
        out = engine.compare(spec_a, spec_b, [target], predict_fn=self._predict)
        if not out.get("comparable"):
            return None, str(out.get("reason"))[:200]
        rows = out["ratios"].get("rows") or []
        row = rows[0] if len(rows) == 1 else next((r for r in rows if str(r["compound"]).lower() == target), None)
        ratio = row.get("ratio_a_over_b") if row else None
        if not isinstance(ratio, (int, float)) or not math.isfinite(ratio) or ratio <= 0:
            return None, str((row or {}).get("note") or "no ratio (a zero or undefined arm)")[:200]
        return math.log10(ratio), ""


def engine_prediction(item: Mapping[str, Any], specs: Sequence[engine.FormulationSpec], runner: Runner
                      ) -> Tuple[Optional[float], str]:
    """The engine's log ratio per lever step for the item, by the item's own estimator (prereg sec. 5)."""
    target = target_alias(item["compound"])
    if item["kind"] == "categorical":
        p, why = runner.log_ratio(specs[1], specs[0], target)
    else:
        logs = [0.0]
        for s in specs[1:]:
            r, why = runner.log_ratio(s, specs[0], target)
            if r is None:
                return None, why
            logs.append(r)
        p, why = wls_slope(item["x"], logs, item["var"]), ""
    if p is None:
        return None, why
    return (0.0 if abs(p) * item["span_steps"] < FLOOR_DEX else p), ""


# ---------------------------------------------------------------------------------------------------- the scores
def _ci(boot: Sequence[float]) -> List[float]:
    lo, hi = np.percentile(np.asarray(boot, float), [5, 95], method="inverted_cdf")
    return [float(lo), float(hi)]


def compare_errors(rows: Sequence[Mapping[str, Any]], rng: np.random.Generator) -> Dict[str, Any]:
    """Median absolute error of the engine, 'no change' and the table on the same rows; each difference (engine minus
    rival, negative = the engine is closer) with its 90 % interval from a paired bootstrap over laboratories."""
    labs = sorted({r["lab"] for r in rows})
    by_lab = {lab: [r for r in rows if r["lab"] == lab] for lab in labs}

    def stats(pool):
        e, z, t = (np.median([r[k] for r in pool]) for k in ("err_engine", "err_no_change", "err_table"))
        return float(e), float(z), float(t)

    e, z, t = stats(rows)
    out = {"contrasts": len(rows), "laboratories": len(labs), "lab_names": labs,
           "median_abs_error": {"engine": e, "no_change": z, "table": t},
           "table_estimates": sum(1 for r in rows if r["table"] is not None),
           "engine_minus_no_change": {"difference": e - z}, "engine_minus_table": {"difference": e - t}}
    if len(labs) >= MIN_LABS:
        boot = []
        for _ in range(N_BOOT):
            pool = [r for j in rng.integers(0, len(labs), size=len(labs)) for r in by_lab[labs[j]]]
            be, bz, bt = stats(pool)
            boot.append((be - bz, be - bt))
        out["engine_minus_no_change"]["ci90"] = _ci([b[0] for b in boot])
        out["engine_minus_table"]["ci90"] = _ci([b[1] for b in boot])
    for key in ("engine_minus_no_change", "engine_minus_table"):
        out[key]["reading"] = verdict(out[key].get("ci90"))
    out["beats_both"] = all(out[k]["reading"] == "beats" for k in ("engine_minus_no_change", "engine_minus_table"))
    return out


def verdict(ci: Optional[Sequence[float]]) -> str:
    if ci is None:
        return "not judged (one laboratory)"
    return "beats" if ci[1] < 0 else "loses" if ci[0] > 0 else "ties"


def sign_accuracy(rows: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """On the rows whose measured change over the whole contrast clears the measurement noise: the share where each
    prediction has the measured sign. A prediction of exactly zero is a miss; 'no change' has no sign."""
    clear = [r for r in rows if abs(r["m"]) * r["span_steps"] > NOISE_DEX]

    def share(key):
        return (sum(1 for r in clear if r[key] and math.copysign(1, r[key]) == math.copysign(1, r["m"])) / len(clear)
                if clear else None)

    return {"contrasts": len(clear), "engine": share("engine"), "table": share("table")}


def score_rows(rows: Sequence[Dict[str, Any]], rng: np.random.Generator) -> Dict[str, Any]:
    """Errors on the rows with no censored arm (a censored arm makes the measured ratio a bound: prereg sec. 3);
    sign accuracy on all of them."""
    exact = [r for r in rows if not r["censored_arms"]]
    out = compare_errors(exact, rng) if exact else {"contrasts": 0, "laboratories": 0, "beats_both": False}
    out["sign_accuracy"] = sign_accuracy(rows)
    out["contrasts_with_a_censored_arm"] = len(rows) - len(exact)
    return out


def score_reading(items: Sequence[Mapping[str, Any]], reading: str, rival: Mapping, runner: Runner) -> Dict[str, Any]:
    rng = np.random.default_rng(SEED)
    rows, not_simulable = [], []
    for it in items:
        state, specs = classify(it, reading)
        ident = {k: it.get(k) for k in ("cell", "lab", "source", "contrast", "level", "against", "kind")}
        if specs is None:
            not_simulable.append({**ident, "why": state})
            continue
        p, why = engine_prediction(it, specs, runner)
        if p is None:
            not_simulable.append({**ident, "why": "the engine refuses the ratio: " + why})
            continue
        t = rival_prediction(rival, it)
        rows.append({**ident, "lever_class": it["lever_class"], "m": it["m"], "span_steps": it["span_steps"],
                     "censored_arms": it["censored_arms"], "engine": p, "table": t,
                     "err_engine": abs(p - it["m"]), "err_no_change": abs(it["m"]),
                     "err_table": abs((t or 0.0) - it["m"])})
    cells = []
    for cell in sorted({tuple(it["cell"]) for it in items}):
        mine = [r for r in rows if tuple(r["cell"]) == cell]
        total = sum(1 for it in items if tuple(it["cell"]) == cell)
        rec: Dict[str, Any] = {"cell": list(cell), "test_contrasts": total, "simulated": len(mine)}
        if mine:
            rec.update(score_rows(mine, rng))
        else:
            rec.update({"beats_both": False, "not_simulable": dict(Counter(
                n["why"] for n in not_simulable if tuple(n["cell"]) == cell).most_common())})
        # prereg sec. 5: only the front-door reading displays a cell
        rec["displayed"] = bool(reading == PRIMARY and rec["beats_both"])
        cells.append(rec)
    classes = {}
    for lc in sorted({it["lever_class"] for it in items}):
        mine = [r for r in rows if r["lever_class"] == lc]
        classes[lc] = {"test_contrasts": sum(1 for it in items if it["lever_class"] == lc), "simulated": len(mine),
                       **(score_rows(mine, rng) if mine else {})}
    return {"reading": reading, "test_contrasts": len(items), "simulated": len(rows),
            "not_simulable": dict(Counter(n["why"] for n in not_simulable).most_common()),
            "overall": score_rows(rows, rng) if rows else None, "lever_classes": classes, "cells": cells,
            "displayed_cells": [c["cell"] for c in cells if c.get("displayed")],
            "rows": rows, "not_simulable_rows": not_simulable}


# ------------------------------------------------------------------------------------------------------ payloads
def _why_class(why: str) -> str:
    """The reason without its instance (the compound or arm it names), for the counts."""
    return why.split(":")[0]


def _by_class(reasons: Mapping[str, int]) -> Dict[str, int]:
    out: Counter = Counter()
    for why, n in reasons.items():
        out[_why_class(why)] += n
    return dict(out.most_common())


def items_payload() -> Dict[str, Any]:
    sources = LT.load_sources()
    labs, _ = LT.lab_partition()
    excluded = excluded_sources()
    items, meta = build_items(sources, labs, excluded)
    _, rival_sha, rival_digest = rival_cells()
    cells = []
    for cell in sorted({tuple(it["cell"]) for it in items}):
        mine = [it for it in items if tuple(it["cell"]) == cell]
        rec = {"cell": list(cell), "engine_target": target_alias(cell[0]), "test_contrasts": len(mine),
               "laboratories": sorted({it["lab"] for it in mine}), "sources": sorted({it["source"] for it in mine})}
        for reading in READINGS[:2]:
            states = [classify(it, reading)[0] for it in mine]
            sim = [it for it, s in zip(mine, states) if s == "simulable"]
            rec[reading] = {"simulable": len(sim), "simulable_laboratories": len({it["lab"] for it in sim}),
                            "not_simulable": dict(Counter(s for s in states if s != "simulable").most_common())}
        cells.append(rec)
    odorants = {}
    for cell in cells:
        o = odorants.setdefault(cell["cell"][0], {"engine_target": cell["engine_target"], "cells": 0})
        o["cells"] += 1
    coverage = {}
    for reading in READINGS[:2]:
        states = [classify(it, reading)[0] for it in items]
        coverage[reading] = {
            "simulable_contrasts": sum(1 for s in states if s == "simulable"),
            "cells_with_two_simulable_laboratories": sum(1 for c in cells
                                                         if c[reading]["simulable_laboratories"] >= MIN_LABS),
            "not_simulable_by_reason": dict(Counter(_why_class(s) for s in states if s != "simulable").most_common())}
    return {"artifact": "step5a_display_check_items", "prereg": str(PREREG.relative_to(ROOT)),
            "excluded_sources": excluded, "counts": {**meta, "items": len(items)},
            "rival": {"file": str(RIVAL.relative_to(ROOT)), "sha256": rival_sha, "fit_digest": rival_digest},
            "coverage": coverage, "odorants": odorants, "cells": cells}


def score_payload() -> Dict[str, Any]:
    frozen = json.loads(ITEMS_OUT.read_text())
    live = items_payload()
    for key in ("excluded_sources", "counts", "coverage", "cells"):
        if frozen[key] != live[key]:
            raise SystemExit(f"the items record ({ITEMS_OUT.name}) no longer matches the keyed data at {key!r}; "
                             "the pre-registered item set would change. Stop.")
    rival, rival_sha, digest = rival_cells()
    if digest != RIVAL_FIT_DIGEST:
        raise SystemExit(f"the DEVELOPMENT-only table's fit estimates changed (digest {digest}); the rival is frozen "
                         "by the pre-registration. Stop.")
    labs, _ = LT.lab_partition()
    items, _ = build_items(LT.load_sources(), labs, excluded_sources())
    runner = Runner()
    readings = {r: score_reading(items, r, rival, runner) for r in READINGS}
    return {"artifact": "step5a_display_check", "prereg": str(PREREG.relative_to(ROOT)),
            "primary_reading": PRIMARY, "rival": {**live["rival"], "sha256_at_scoring": rival_sha},
            "excluded_sources": live["excluded_sources"], "counts": live["counts"],
            "engine": engine.engine_metadata(), "integrations": runner.integrations,
            "displayed_cells": readings[PRIMARY]["displayed_cells"], "readings": readings}


# ----------------------------------------------------------------------------------------------------- rendering
def _cell(c: Sequence[str]) -> str:
    return f"{c[0]} x {c[1]}" + ("" if c[2] == "all" else f" ({c[2]})")


def _reasons(d: Mapping[str, int], n: int = 3) -> str:
    return "; ".join(f"{k} ({v})" for k, v in list(d.items())[:n]) or "-"


def render_items(p: Mapping[str, Any]) -> str:
    c = p["counts"]
    L = ["# Step 5a, the item set: which TEST contrasts the engine can be run on", "",
         "*Generated by `scripts/generators/build_step5a_display_check.py --items` under "
         f"`{p['prereg']}`. Built from the keyed contrasts and the engine's envelope alone: no integration, no engine "
         "output. The numbers here go into the pre-registration before any score run.*", "",
         f"- TEST contrasts in the lever table: {c['test_observations']}; after the exclusions "
         f"{c['test_observations_after_exclusions']}, in {c['test_cells']} cells, {c['qualifying_cells']} of them "
         f"with two or more TEST laboratories ({c['items']} items: one per numeric contrast, one per categorical "
         "pair).",
         "- Excluded sources: " + "; ".join(f"`{k}` ({v}; {c['excluded_observations'].get(k, 0)} contrasts)"
                                            for k, v in p["excluded_sources"].items()) + ".",
         f"- Rival: `{p['rival']['file']}`, sha256 `{p['rival']['sha256']}`, fit digest `{p['rival']['fit_digest']}`.",
         ""]
    for reading, title in (("front_door", "Front-door reading (decides display)"),
                           ("added_precursors", "Added-precursors reading (secondary; the matrix material left out)")):
        cov = p["coverage"][reading]
        L += [f"## {title}", "",
              f"{cov['simulable_contrasts']} of {c['items']} items map to an engine run; "
              f"{cov['cells_with_two_simulable_laboratories']} of {c['qualifying_cells']} cells keep two laboratories.",
              "", "| why an item does not map | items |", "|---|---|"]
        L += [f"| {k} | {v} |" for k, v in cov["not_simulable_by_reason"].items()] + [""]
    L += ["## Odorants", "", "| odorant | engine target | cells |", "|---|---|---|"]
    L += [f"| {k} | {v['engine_target'] or 'none: not simulable'} | {v['cells']} |"
          for k, v in sorted(p["odorants"].items())]
    L += ["", "## Cells", "",
          "| cell | TEST contrasts | laboratories | front door: simulable (labs) | added precursors: simulable (labs) "
          "| main reasons (added precursors) |", "|---|---|---|---|---|---|"]
    for cell in p["cells"]:
        fd, ap = cell["front_door"], cell["added_precursors"]
        L.append(f"| {_cell(cell['cell'])} | {cell['test_contrasts']} | {len(cell['laboratories'])} | "
                 f"{fd['simulable']} ({fd['simulable_laboratories']}) | {ap['simulable']} "
                 f"({ap['simulable_laboratories']}) | {_reasons(ap['not_simulable'], 2)} |")
    return "\n".join(L) + "\n"


def _d(x: Optional[float]) -> str:
    return "-" if x is None else f"{x:+.2f}"


def _diff(d: Optional[Mapping[str, Any]]) -> str:
    if not d or "difference" not in d:
        return "-"
    ci = d.get("ci90")
    return f"{_d(d['difference'])}" + (f" ({_d(ci[0])} to {_d(ci[1])})" if ci else "") + f": {d['reading']}"


def _score_line(label: str, total: int, s: Optional[Mapping[str, Any]], simulated: int) -> str:
    if not s or not s.get("contrasts"):
        return f"| {label} | {total} | {simulated} | - | - | - | - | - | - | - |"
    e, sa = s["median_abs_error"], s["sign_accuracy"]
    sign = ("-" if sa["engine"] is None else f"{sa['engine']:.0%} / "
            + (f"{sa['table']:.0%}" if s["table_estimates"] else "no estimate") + f" (n = {sa['contrasts']})")
    return (f"| {label} | {total} | {simulated} | {s['laboratories']} | {e['engine']:.2f} | {e['no_change']:.2f} | "
            f"{e['table']:.2f} | {_diff(s['engine_minus_no_change'])} | {_diff(s['engine_minus_table'])} | {sign} |")


HEAD = ("| {what} | TEST contrasts | simulated | labs | engine | no change | table | engine - no change (90 %) | "
        "engine - table (90 %) | sign: engine / table |\n|---|---|---|---|---|---|---|---|---|---|")


def render_score(p: Mapping[str, Any]) -> str:
    c = p["counts"]
    primary = p["readings"][PRIMARY]
    L = ["# Step 5a: the step-5 display rule on the engine as shipped (a diagnostic)", "",
         "*Generated by `scripts/generators/build_step5a_display_check.py --score` under "
         f"`{p['prereg']}`. No constant or structure of the engine is changed and nothing ships on this result.*", "",
         f"**Result.** Of {c['qualifying_cells']} lever x odorant cells with TEST contrasts from two or more "
         f"laboratories ({c['items']} contrasts), the engine as a formulator would call it can be run on "
         f"**{primary['simulated']}** contrasts, and **{len(p['displayed_cells'])}** cells would be displayed"
         + (": " + "; ".join(_cell(x) for x in p["displayed_cells"]) if p["displayed_cells"] else "") + ".", "",
         "Errors are median absolute errors of the log10 ratio per lever step (dex); a difference is the engine's "
         "error minus the rival's, so negative means the engine is closer; 'beats' / 'loses' need the 90 % interval "
         "(paired bootstrap over laboratories) to exclude zero. With two laboratories that interval excludes zero "
         "only when each laboratory alone and both together agree. Where the DEVELOPMENT-only table has no estimate "
         "for a contrast its prediction is zero, the same as 'no change'.", "",
         "- Excluded sources: " + "; ".join(f"`{k}` ({v})" for k, v in p["excluded_sources"].items()) + ".",
         f"- Rival: `{p['rival']['file']}` (fit digest `{p['rival']['fit_digest']}`).", ""]
    titles = {"front_door": "Front-door reading (decides display): every keyed component must be an engine "
                            "precursor or an ingredient on file",
              "added_precursors": "Added-precursors reading (secondary, never displayed): the matrix material is "
                                  "left out, the added pure precursors are charged, buffer as keyed",
              "added_precursors_ph_clamped": "Added-precursors reading with the pH held at the keyed value "
                                             "(sensitivity)"}
    for reading in READINGS:
        r = p["readings"][reading]
        L += [f"## {titles[reading]}", "",
              f"{r['simulated']} of {r['test_contrasts']} contrasts simulated. Not simulated: "
              + _reasons(_by_class(r["not_simulable"]), 12) + ".", ""]
        if not r["simulated"]:
            continue
        L += [HEAD.format(what="lever class")]
        L += [_score_line(lc, s["test_contrasts"], s if s.get("contrasts") else None, s["simulated"])
              for lc, s in r["lever_classes"].items()]
        L += [_score_line("**all**", r["test_contrasts"], r["overall"], r["simulated"]), "",
              HEAD.format(what="cell")]
        L += [_score_line(_cell(x["cell"]) + (" **(displayed)**" if x.get("displayed") else
                                               " (beats both; secondary reading, not a display)"
                                               if x.get("beats_both") else ""),
                          x["test_contrasts"], x if x.get("contrasts") else None, x["simulated"])
              for x in r["cells"] if x["simulated"]]
        L += ["", "| contrast | cell | lab | measured | engine | table |", "|---|---|---|---|---|---|"]
        L += [f"| {x['contrast']}" + (f" ({x['level']} vs {x['against']})" if x.get("level") else "")
              + f" | {_cell(x['cell'])} | {x['lab']} | {_d(x['m'])}" + (" (bound)" if x["censored_arms"] else "")
              + f" | {_d(x['engine'])} | " + ("no estimate" if x["table"] is None else _d(x["table"])) + " |"
              for x in r["rows"]]
        L.append("")
    L += ["## Cells the engine cannot simulate (front-door reading)", "",
          "| cell | TEST contrasts | why |", "|---|---|---|"]
    L += [f"| {_cell(x['cell'])} | {x['test_contrasts']} | {_reasons(x.get('not_simulable') or {}, 2)} |"
          for x in primary["cells"] if not x["simulated"]]
    return "\n".join(L) + "\n"


def _inputs() -> List[Path]:
    files = sorted(data_paths.LEVER_TABLE_CONTRASTS_DIR.glob("*.yml"))
    return [q for q in files if q.name != "TEMPLATE.yml" and ".reader_" not in q.name] + [
        data_paths.CONTRAST_SPLIT, data_paths.LEVER_TABLE_NEW_LABS]


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--items", action="store_true")
    mode.add_argument("--score", action="store_true")
    args = ap.parse_args(argv)
    payload, out, render = ((items_payload(), ITEMS_OUT, render_items) if args.items
                            else (score_payload(), SCORE_OUT, render_score))
    # the rival's file is not a hashed input: its estimates are frozen by RIVAL_FIT_DIGEST instead (prereg sec. 6)
    payload["provenance"] = provenance.provenance_block(
        out.stem, generated_by="scripts/generators/build_step5a_display_check.py",
        inputs=_inputs() + ([PREREG, ITEMS_OUT] if args.score else []))
    artifact_io.write_artifact(payload, out, render=render)
    print(f"wrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
