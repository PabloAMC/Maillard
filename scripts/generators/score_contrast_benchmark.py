"""
The contrast benchmark, scoring step 1: coverage of the held-out contrasts (Q1), and the shipped engine on the clean
held-out pairs it can answer (Q2). Pre-registration: results/validation/kinetic_core_contrast_benchmark_score_prereg.md
(v6); every rule below is that document's, and a rule that differs from it is a bug here.

Two modes, run in this order and never the other way round:

* ``--items`` builds the item sets from the YAML and the envelope alone (``engine.declare_envelope``, no integration):
  the frozen name map, the Q1 coverage of both families (the 65 v2 TEST atlas cells; the one-factor pairs of the 15
  keyed held-out tables), the Q2 items with their exposure tags and ladders, and the power line. Writes
  ``kinetic_core_contrast_benchmark_items.{json,md}``. Its numbers go into the prereg BEFORE any score run.
* ``--score`` asserts the items record, the hand table's hash and the shipped configuration, then integrates one
  ``engine.predict`` per (pot, species) and writes ``kinetic_core_contrast_benchmark_score_v1.{json,md}``.

A diagnostic (STRATEGY_2026-10-29.md sec. 9): nothing ships, no constant moves, no split changes on its result.

    ./scripts/docker_maillard.sh run "python scripts/generators/score_contrast_benchmark.py --items"
    ./scripts/docker_maillard.sh run "python scripts/generators/score_contrast_benchmark.py --score"
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import statistics
import subprocess
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import artifact_io, data_access, data_paths, provenance  # noqa: E402
from src.kinetic_core import engine  # noqa: E402
from src.kinetic_core.config import Config, applied  # noqa: E402
from src.kinetic_core.fit_rows import FIT_ROWS_DIR  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402
from src.kinetic_core.sulfur import PhFixedPointNotConverged  # noqa: E402

SPLIT_V2 = ROOT / "docs" / "reference" / "contrast_benchmark" / "split_v2_per_change.yml"
ITEMS_OUT = data_paths.VALIDATION_DIR / "kinetic_core_contrast_benchmark_items.json"
SCORE_OUT = data_paths.VALIDATION_DIR / "kinetic_core_contrast_benchmark_score_v1.json"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_contrast_benchmark_score_prereg.md"

KEYED_TABLES = (
    "heldout_deroos2005_fig4", "heldout_du2023_tab2", "heldout_grosch1992_tab2", "heldout_larralde1995",
    "heldout_lee2024_fig2", "heldout_liu2023b_tab3", "heldout_meynier1995_tab1_2_4", "heldout_sedehi2013_tab1",
    "heldout_thornalley1999_figs1b_2e_2f_3c", "heldout_wang2022_tab2", "heldout_xu2007_tab1",
    "heldout_zamora2015_tab1", "heldout_zhang2021_table3", "heldout_zhao2017b_tab1", "heldout_zhao2019b_tab1_2",
)
SCORED_TABLES = ("heldout_grosch1992_tab2", "heldout_meynier1995_tab1_2_4", "heldout_xu2007_tab1")

#: Prereg v6 sec. 2C: the atlas cells each keyed table reproduces (README sec. 11). A tag, not a de-duplication: the
#: two Q1 families are reported apart and never pooled.
REPRODUCES_ATLAS = {
    "heldout_thornalley1999_figs1b_2e_2f_3c": ["F05#23@thornalley1999", "F06#11@thornalley1999", "F08#2@thornalley1999"],
    "heldout_zamora2015_tab1": ["F30#15@zamora2015"],
    "heldout_liu2023b_tab3": ["F23#6@liu2023b", "F24#6@liu2023b", "F28#26@liu2023b"],
    "heldout_zhao2019b_tab1_2": ["F23#7@zhao2019b", "F24#4@zhao2019b", "F25#23@zhao2019b"],
    "heldout_du2023_tab2": ["F23#8@du2023", "F24#5@du2023", "F26#24@du2023"],
    "heldout_wang2022_tab2": ["F28#25@wang2022"],
    "heldout_grosch1992_tab2": ["F27#15@grosch1992 (Table II only)", "F27#16@grosch1992"],
}

# --- Q1, atlas family: the hand table (prereg v6 sec. 2A) ------------------------------------------------------------

#: Classes 1-5 by hand, read from the atlas cell (contrast, quantity, tags), first that applies. `inputs_missing` only
#: where the cell says a charge / time / temperature / pH is not printed (incl. `ph_basis: not_measured`), or the system
#: is a food whose precursor charges the cell does not give. A cell is answerable if ANY of its legs is: such cells
#: carry a spec in ATLAS_SPECS naming the leg, and the envelope assigns class 6 or 7. No hand-classed cell has an
#: answerable leg.
ATLAS_COVERAGE: Dict[str, Tuple[str, str]] = {
    "F02#1@potman1989": ("operator_missing", "glucose-loss rate gain with phosphate (rate)"),
    "F03#8@zyzak1995": ("operator_missing", "product shares at a similar extent of fFL loss (percent of loss)"),
    "F03#9@zyzak1995": ("operator_missing", "fFL half-life air vs N2 + chelators (rate)"),
    "F04#2@kallen1971b": ("operator_missing", "thiazolidine formation rate vs phosphate (rate)"),
    "F04#26@mizutani1960": ("precursor_unmapped", "aldehyde thiazolidines of cysteine (formaldehyde ... furfural) have no engine name"),
    "F05#23@thornalley1999": ("operator_missing", "3-DG initial rate, uM/day (rate)"),
    "F06#11@thornalley1999": ("operator_missing", "3-DG formation rate vs phosphate / chelator (rate)"),
    "F08#2@thornalley1999": ("operator_missing", "glyoxal initial rate vs chelator / phosphate (rate)"),
    "F16#11@monforte2021": ("operator_missing", "headspace release rate (Weibull k)"),
    "F16#14@huang2017": ("inputs_missing", "absolute concentrations not stated"),
    "F16#15@balagiannis2009": ("species_absent", "2- / 3-methylbutanal: no engine target"),
    "F16#16@martin2001": ("inputs_missing", "re-infused potato slices: the glucose / amino-acid charges are not in the cell"),
    "F18#16@mandin1999": ("precursor_unmapped", "linoleic acid has no engine name"),
    "F19#0@amranihemaimi1995": ("operator_missing", "13C label share of pyrazines"),
    "F19#1@low2007": ("inputs_missing", "potato flake cakes: sugar charges not given"),
    "F19#13@vanlancker2012": ("precursor_unmapped", "dipeptides (GlyGly, SerGly, ...) have no engine name"),
    "F19#17@shibamoto1977": ("inputs_missing", "pH not printed (8 M NH4OH); also total pyrazines, no engine target"),
    "F20#1@amranihemaimi1995": ("operator_missing", "13C label share of 2,6-DMP"),
    "F20#2@low2007": ("inputs_missing", "potato cakes: sugar charges not given"),
    "F20#3@low2007": ("inputs_missing", "potato cakes: sugar charges not given"),
    "F20#13@shibamoto1977": ("inputs_missing", "pH not printed (8 M NH4OH); also a share of total pyrazines, no engine target"),
    "F21#4@hwang1995b": ("species_absent", "pyridines: no engine target"),
    "F21#6@kim1996b": ("operator_missing", "15N share of 2-pentylpyridine N"),
    "F21#8@zamora2020": ("species_absent", "alkylpyridines: no engine target"),
    "F21#9@zhang1989": ("species_absent", "2-pentylpyridine: no engine target"),
    "F21#10@elmore1997": ("operator_missing", "percent of chromatogram area"),
    "F22#5@mizutani1960": ("inputs_missing", "pH not measured (cell tag); also aldehyde thiazolidines of cysteine, no engine name"),
    "F22#6@obata1960": ("inputs_missing", "pH not measured (cell tag); also cystine, no engine name"),
    "F22#12@shu1985": ("species_absent", "trithiolanes / thiazoles / MAA-derived rings: no engine target"),
    "F22#14@zhang1991b": ("species_absent", "cyclic polysulfides, thiophenes: no engine target"),
    "F22#15@tai1997": ("species_absent", "total sulfur volatiles: no engine target"),
    "F23#5@wang2012": ("operator_missing", "isotopomer share of MFT"),
    "F23#6@liu2023b": ("inputs_missing", "reactant amounts not printed (README sec. 11)"),
    "F23#7@zhao2019b": ("operator_missing", "isotopomer share of MFT"),
    "F23#8@du2023": ("operator_missing", "isotopomer share of MFT"),
    "F24#4@zhao2019b": ("operator_missing", "isotopomer share of MFT"),
    "F24#5@du2023": ("operator_missing", "isotopomer share of MFT"),
    "F24#6@liu2023b": ("inputs_missing", "reactant amounts not printed (README sec. 11)"),
    "F25#17@us8920862": ("factor_not_translated", "structured lipid phase vs aqueous buffer: no engine input"),
    "F25#23@zhao2019b": ("operator_missing", "isotopomer share of FFT"),
    "F26#22@farmer1990": ("species_absent", "mercaptoketones: no engine target"),
    "F26#24@du2023": ("species_absent", "2-pentylthiophene: no engine target"),
    "F27#5@wang2012": ("operator_missing", "isotopomer share of MFT"),
    "F27#12@hincelin1992": ("operator_missing", "percent of the SDE isolate (open stripping)"),
    "F27#16@grosch1992": ("precursor_unmapped", "thiamine pyrophosphate has no engine name"),
    "F27#17@grosch1992": ("inputs_missing", "H2S dosed as gas flushes: the dissolved charge is not printed"),
    "F27#18@thomas2014": ("inputs_missing", "cooked ham: brine doses, not pot concentrations"),
    "F27#20@hartman1984b": ("precursor_unmapped", "cystine, ascorbic acid, MSG have no engine name"),
    "F28#8@ehrenberg1989": ("operator_missing", "cysteine autoxidation rate vs metals"),
    "F28#20@xu2011": ("precursor_unmapped", "lard (fresh / oxidised) has no engine name"),
    "F28#25@wang2022": ("precursor_unmapped", "the T x t grid pot carries 1 % methyl linoleate, no engine name (CF-028)"),
    "F28#26@liu2023b": ("inputs_missing", "reactant amounts not printed (README sec. 11)"),
    "F29#22@low2006": ("inputs_missing", "potato cake: sugar / Asn charges not given"),
    "F29#23@kocadagli2021": ("inputs_missing", "potato-flake dough: charges not given"),
    "F29#24@goncuoglu2026": ("inputs_missing", "extrusion; sugar not stated"),
    "F30#4@whitfield1988": ("species_absent", "2-pentylpyridine, alkylthiophenes, thiapyran: no engine target"),
    "F30#9@du2023": ("species_absent", "2-pentylthiophene: no engine target"),
    "F30#11@kim1996b": ("operator_missing", "15N share of pentylpyridine N"),
    "F30#15@zamora2015": ("species_absent", "phenylacetaldehyde / phenylethylamine: no engine target"),
    "F30#16@hidalgo1993": ("operator_missing", "browning / fluorescence rate constants"),
}

#: The atlas cells whose class the envelope decides: both arms as written in the cell (names are front-door aliases).
#: (precursors, t_c, minutes, label_ph, buffer kwargs, target, note). One arm per entry; the contrast's other arm(s) in
#: `arms`.
ATLAS_SPECS: Dict[str, Dict[str, Any]] = {
    "F16#17@mandin1999": {"target": "methional", "t_c": 102.5, "ph": 6.0, "buffer": {"kind": "clamped"},
                          "arms": [({"methionine": 140, "glucose": 72}, 120), ({"methionine": 140}, 120)],
                          "note": "Met 0.14 M +/- glucose 0.072 M, 2 h reflux at 100 - 105 C, pH-stat 6"},
    "F22#7@li2022": {"target": "hydrogen sulfide", "t_c": 90, "ph": 7.0,
                     "buffer": {"kind": "phosphate", "phosphate_mol_l": 0.01},
                     "arms": [({"cysteine": 5.7}, 10), ({"cysteine": 5.7}, 120)],
                     "note": "Cys 5.7 mM, 10 mM phosphate pH 7, 90 C, 10 vs 120 min"},
    "F26#11@meynier1995": {"target": "2-acetylthiazole", "t_c": 140, "ph": None,
                           "buffer": {"kind": "pyrophosphate", "pyrophosphate_mol_l": 0.2},
                           "arms": [({"cysteine": 50, "ribose": 35}, 60, 4.5), ({"cysteine": 50, "ribose": 35}, 60, 6.5)],
                           "note": "Meynier Table 3 (cysteine; NOT in the keyed file), the 2-acetylthiazole leg; "
                                   "3-mercapto-2-pentanone and thiazole have no engine target"},
    "F27#15@grosch1992": {"target": "mft", "t_c": 100, "ph": 5.7,
                          "buffer": {"kind": "pyrophosphate", "pyrophosphate_mol_l": 0.2},
                          "arms": [({"thiamine": 0.06}, 120), ({"ribose": 1.3, "cysteine": 0.65}, 120),
                                   ({"thiamine": 0.06, "cysteine": 0.65}, 120)],
                          "note": "Tables II, III, VIII; the keyed file holds Table II only"},
    "F29#11@quan2020": {"target": "acrylamide", "t_c": None, "ph": 7.0,
                        "buffer": {"kind": "phosphate", "phosphate_mol_l": 0.1},
                        "arms": [({"asparagine": 200, "glucose": 100}, 21, None, 100),
                                 ({"asparagine": 200, "glucose": 100}, 21, None, 130)],
                        "note": "the temperature leg (Asn 200 + Glc 100 mM, 0.1 M phosphate pH 7, 100 vs 130 C, 21 min); "
                                "the mixture leg charges Trp, which has no engine name"},
}

# --- the frozen name map (prereg v6 sec. 3) --------------------------------------------------------------------------


def _alias_for(key: str, aliases: Mapping[str, str]) -> Optional[str]:
    """Rule: a YAML key K maps to the FIRST alias (dict order) whose value equals K, case-insensitively; none -> unmapped."""
    for alias, value in aliases.items():
        if value.lower() == key.lower():
            return alias
    return None


def precursor_name(key: str) -> Optional[str]:
    return _alias_for(key, engine.PRECURSOR_ALIASES)


def target_name(key: str) -> Optional[str]:
    return _alias_for(key, engine.TARGET_ALIASES)


# --- loading -----------------------------------------------------------------------------------------------------------

#: Prereg v6 sec. 2B: the factor fields of the broad (Q1) pair rule. Composition keys count one factor each.
SCALAR_FACTORS = ("label_ph", "t_c", "minutes", "buffer", "phosphate_mM", "lipid", "phospholipid_pct", "aldehyde",
                  "post_heat", "atmosphere", "protein_mg_per_mL", "matrix")
#: Descriptive fields that restate a composition key: they count only when no composition key differs.
RESTATING = ("aldehyde", "lipid", "phospholipid_pct")
#: Class 5 (prereg v6): the factors this scorer does not translate into the engine's inputs, and why. Three have no
#: engine input at all; the rest exist on `ProcessSpec` but the keyed tables give them as free text or in a form the
#: engine does not take, and translating them is not this step's job.
NOT_TRANSLATED = {
    "lipid": "lipid: no input on the sulfur lane (a lipid carrier there is a lane conflict)",
    "phospholipid_pct": "phospholipid: no input on the sulfur lane (a lipid carrier there is a lane conflict)",
    "aldehyde": "lipid aldehyde: no engine name",
    "post_heat": "reheating: no engine input",
    "buffer": "buffer species / concentration: BufferSpec exists; the table's free text is not translated here",
    "phosphate_mM": "phosphate concentration: BufferSpec exists; not translated here",
    "atmosphere": "atmosphere: the engine takes argon / air / air_cu; the table's N2 / free text is not translated here",
    "protein_mg_per_mL": "protein: ProcessSpec.protein_g_per_l exists; not translated here",
    "matrix": "matrix: ProcessSpec.matrix exists; not translated here",
}
#: Units no engine observation operator produces (class 2).
OPERATOR_UNITS = {"percent_remaining": "percent remaining", "percent_of_isotopomers": "isotopomer %",
                  "percent_of_isotopologues": "isotopologue %", "per_s": "rate", "per_min": "rate",
                  "per_M_per_min": "rate", "per_M_per_s": "rate", "log10_per_M_per_s": "rate",
                  "kJ_per_mol": "activation energy", "per_pH": "pH slope", "ratio": "rate ratio", "uM_per_day": "rate"}
#: Prereg v6 sec. 2D: detection bounds for censored arms, in the table's own unit.
ND_BOUND = {"heldout_meynier1995_tab1_2_4": 0.01, "heldout_xu2007_tab1": 0.07}
TR_BOUND = {"heldout_meynier1995_tab1_2_4": ("trace_below_0p1", 0.1), "heldout_xu2007_tab1": ("trace_below_0p3", 0.3)}
EXPOSED = {("heldout_meynier1995_tab1_2_4", "FUR"), ("heldout_xu2007_tab1", "MFT"), ("heldout_xu2007_tab1", "FFT")}
#: Prereg v6 sec. 2D: each paper's own statement of the scatter of a printed value (relative SD of the reported
#: mean). Grosch prints SD and N per row (used per row). Meynier: "triplicate analyses, CV < 20 %" -> 0.20 / sqrt(3).
#: Xu: "mean CV < 20 %", n not stated -> 0.20. The row's read error (dex) is added in quadrature.
STATED_REL_SD = {"heldout_meynier1995_tab1_2_4": 0.20 / math.sqrt(3), "heldout_xu2007_tab1": 0.20}
#: Prereg v6 sec. 4: the numerical floor of a predicted log-ratio is measured per item, as |p(shipped) - p(tight)|
#: with the sulfur lane's pH fixed point tightened from its shipped 1e-3 pH (the lane's numerical floor; config.py
#: "ph.fixed_point.tolerance") to TIGHT_PH_TOLERANCE, near the residual floor of 3e-5 - 6e-4 pH. A |p| at or below
#: max(that, FLOAT_REL_TOL / ln 10) is numerical, not a prediction, and counts as p = 0.
TIGHT_PH_TOLERANCE = 1e-4
FLOAT_FLOOR_DEX = provenance.FLOAT_REL_TOL / math.log(10)
#: Descriptive only (a column, never a decision).
BAND_DEX = math.log10(3.0)
SEED = 20261101
N_BOOT = 2000
N_SIM = 2000
#: The per-step agreements at which the detection curve is reported (prereg v6 sec. 5); no point on it is a gate.
AGREEMENT_GRID = (0.5, 0.6, 0.7, 0.8, 0.9)


def load_table(name: str) -> Dict[str, Any]:
    return data_access.load_yaml(FIT_ROWS_DIR / f"{name}.yml")


def _num(x: Any) -> Optional[float]:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def composition(pot: Mapping[str, Any]) -> Dict[str, Any]:
    out = dict(pot.get("initial_mmol_L") or {})
    out.update({f"umol:{k}": v for k, v in (pot.get("amounts_umol") or {}).items()})
    return out


def factors_between(a: Mapping[str, Any], b: Mapping[str, Any], *, broad: bool) -> List[str]:
    ca, cb = composition(a), composition(b)
    if not broad:
        ca = dict(a.get("initial_mmol_L") or {})
        cb = dict(b.get("initial_mmol_L") or {})
    composition_diff = [k for k in sorted(set(ca) | set(cb)) if ca.get(k) != cb.get(k)]
    fields = SCALAR_FACTORS if broad else ("label_ph",)
    return composition_diff + [f for f in fields
                               if a.get(f) != b.get(f) and not (f in RESTATING and composition_diff)]


def species_arms(table: Mapping[str, Any]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """pot -> species -> its arm, from the item sources of prereg v6 sec. 2B in this order of precedence: a printed row;
    an upper limit (nd with the row's own bound); nd; tr; present-not-quantified or not-quantified (a source that is
    never scored, sec. 2D)."""
    out: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
    for r in table["rows"]:
        out[r["pot"]][r["species"]] = {"value": _num(r.get("value")), "unit": r.get("unit"), "sd": _num(r.get("sd")),
                                       "n": _num(r.get("n")), "read_error": _num(r.get("read_error")) or 0.0}
    for pid, pot in table["pots"].items():
        for sp, lim in (pot.get("upper_limits_ug") or {}).items():
            out[pid].setdefault(sp, {"nd": True, "bound": float(lim)})
        for sp in pot.get("not_detected") or []:
            out[pid].setdefault(sp, {"nd": True})
        for key, val in pot.items():
            if key.startswith("trace_below"):
                for sp in val or []:
                    out[pid].setdefault(sp, {"tr": True})
        for sp in list(pot.get("present_not_quantified") or []) + list(pot.get("not_quantified") or {}):
            out[pid].setdefault(sp, {"pnq": True})
    return out


def item_units(table: Mapping[str, Any], a: str, b: str, species: str) -> set:
    """Prereg v6 sec. 2 class 2: the units of the item's own rows in arms a and b; with no row in either arm, the units of
    that species' rows in the pots' ledger entry."""
    units = {r.get("unit") for r in table["rows"] if r["species"] == species and r["pot"] in (a, b)}
    if not units:
        ledger = table["pots"][a].get("ledger_id")
        units = {r.get("unit") for r in table["rows"]
                 if r["species"] == species and table["pots"][r["pot"]].get("ledger_id") == ledger}
    return units


# --- specs and the envelope --------------------------------------------------------------------------------------------


def scored_buffer(table: str) -> BufferSpec:
    assert table in SCORED_TABLES
    return BufferSpec(kind="pyrophosphate", pyrophosphate_mol_l=0.2, declared=True,
                      source=f"{table}: 0.2 M pyrophosphate (prereg v6 sec. 3)")


def pot_spec(table: str, pid: str, pot: Mapping[str, Any], *, clamped: bool = False) -> Optional[engine.FormulationSpec]:
    precursors = {}
    for key, amount in (pot.get("initial_mmol_L") or {}).items():
        name, value = precursor_name(key), _num(amount)
        if name is None or value is None:
            return None
        precursors[name] = value
    t_c, minutes, ph = _num(pot.get("t_c")), _num(pot.get("minutes")), _num(pot.get("label_ph"))
    if not precursors or t_c is None or minutes is None or ph is None:
        return None
    if table in SCORED_TABLES:
        buffer = BufferSpec(kind="clamped", declared=True, source="sensitivity: pH clamped at label_ph") if clamped \
            else scored_buffer(table)
    else:
        buffer = None     # Q1 only: the envelope's verdict; a keyed table's buffer is not translated here
    return engine.FormulationSpec(name=pid, precursors=precursors, process=engine.ProcessSpec(
        thermal=engine.ThermalProgram.isothermal(t_c, minutes), ph=ph, buffer=buffer, vessel=None))


def envelope(spec: engine.FormulationSpec, target: str) -> Tuple[bool, str]:
    d = engine.declare_envelope(spec, [target])
    return d.is_answerable, (d.reasons[0] if d.reasons else d.state)


def atlas_coverage() -> List[Dict[str, Any]]:
    split = data_access.load_yaml(SPLIT_V2)["split"]["errors"]
    errors_of: Dict[str, List[str]] = defaultdict(list)
    for e, v in split.items():
        for c in v["test_contrasts"]:
            errors_of[c].append(e)
    hand, specs = set(ATLAS_COVERAGE), set(ATLAS_SPECS)
    if hand & specs or (hand | specs) != set(errors_of):
        raise SystemExit(f"atlas family != the 65 v2 TEST contrasts: {sorted((hand | specs) ^ set(errors_of))}")
    items = []
    for cid in sorted(errors_of):
        if cid in ATLAS_COVERAGE:
            cls, reason = ATLAS_COVERAGE[cid]
        else:
            s = ATLAS_SPECS[cid]
            verdicts = []
            for arm in s["arms"]:
                prec, minutes = arm[0], arm[1]
                ph = arm[2] if len(arm) > 2 and arm[2] is not None else s["ph"]
                t_c = arm[3] if len(arm) > 3 else s["t_c"]
                buf = BufferSpec(declared=True, source=s["note"], **s["buffer"])
                spec = engine.FormulationSpec(name=cid, precursors=prec, process=engine.ProcessSpec(
                    thermal=engine.ThermalProgram.isothermal(t_c, minutes), ph=ph, buffer=buf))
                verdicts.append(envelope(spec, s["target"]))
            refused = [r for ok, r in verdicts if not ok]
            cls, reason = ("refused_other", refused[0]) if refused else ("answerable", s["note"])
        items.append({"contrast": cid, "errors": errors_of[cid], "class": cls, "reason": reason})
    return items


def atlas_table_sha256() -> str:
    blob = json.dumps({"hand": ATLAS_COVERAGE, "specs": ATLAS_SPECS}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


# --- Q1 keyed family and Q2 items --------------------------------------------------------------------------------------


def keyed_pairs(name: str, table: Mapping[str, Any], *, broad: bool):
    pots = table["pots"]
    for a, b in itertools.combinations(sorted(pots), 2):
        if pots[a].get("ledger_id") != pots[b].get("ledger_id"):
            continue
        f = factors_between(pots[a], pots[b], broad=broad)
        if len(f) == 1:
            yield a, b, f[0]


def classify_pair(name, table, a, b, factor, species, verdicts) -> Tuple[str, str]:
    pots = table["pots"]
    for pid in (a, b):
        p = pots[pid]
        charge = p.get("initial_mmol_L") or {}
        if not charge or any(_num(v) is None for v in charge.values()) or _num(p.get("t_c")) is None or _num(p.get("minutes")) is None \
                or _num(p.get("label_ph")) is None:
            return "inputs_missing", f"{pid}: no printed charge (initial_mmol_L), numeric time, temperature or label pH"
    operator = sorted(u for u in item_units(table, a, b, species) if u in OPERATOR_UNITS)
    if operator:
        return "operator_missing", f"{species} in {operator[0]} ({OPERATOR_UNITS[operator[0]]})"
    if target_name(species) is None:
        return "species_absent", f"{species}: no engine target"
    unmapped = sorted({k for pid in (a, b) for k in (pots[pid].get("initial_mmol_L") or {}) if precursor_name(k) is None})
    if unmapped:
        return "precursor_unmapped", f"no engine name for {', '.join(unmapped)}"
    if factor in NOT_TRANSLATED:
        return "factor_not_translated", NOT_TRANSLATED[factor]
    refused = [verdicts[(pid, species)][1] for pid in (a, b) if not verdicts[(pid, species)][0]]
    if refused:
        return "refused_other", refused[0]
    return "answerable", ""


def build_items() -> Dict[str, Any]:
    tables = {n: load_table(n) for n in KEYED_TABLES}
    name_map = {"precursors": {}, "targets": {}}
    for t in tables.values():
        for pot in t["pots"].values():
            for k in pot.get("initial_mmol_L") or {}:
                name_map["precursors"][k] = precursor_name(k)
        for r in t["rows"]:
            name_map["targets"][r["species"]] = target_name(r["species"])
    verdicts: Dict[Tuple[str, str], Tuple[bool, str]] = {}
    keyed, q2 = [], []
    for name, table in tables.items():
        arms = species_arms(table)
        specs = {pid: pot_spec(name, pid, pot) for pid, pot in table["pots"].items()}
        species = {sp for sps in arms.values() for sp in sps}
        for pid, spec in specs.items():
            for sp in species:
                tg = target_name(sp)
                verdicts[(pid, sp)] = envelope(spec, tg) if spec is not None and tg is not None \
                    else (False, "no spec or no target")
        for a, b, factor in keyed_pairs(name, table, broad=True):
            for sp in sorted(set(arms[a]) | set(arms[b])):
                cls, reason = classify_pair(name, table, a, b, factor, sp, verdicts)
                keyed.append({"table": name, "a": a, "b": b, "factor": factor, "species": sp, "class": cls,
                              "reason": reason, "reproduces_atlas": REPRODUCES_ATLAS.get(name, [])})
        if name in SCORED_TABLES:
            q2.extend(q2_items(name, table, arms, verdicts))
    return {"tables": tables, "name_map": name_map, "keyed": keyed, "q2": q2}


def q2_items(name, table, arms, verdicts) -> List[Dict[str, Any]]:
    """Prereg v6 sec. 2D: the narrow pair rule, items answerable in at least one arm, detected / censored / excluded."""
    pots = table["pots"]
    out = []
    for a, b, factor in keyed_pairs(name, table, broad=False):
        unmapped = sorted({k for pid in (a, b) for k in pots[pid]["initial_mmol_L"] if precursor_name(k) is None})
        for sp in sorted(set(arms[a]) | set(arms[b])):
            if target_name(sp) is None or unmapped:
                continue
            ok_a, ok_b = verdicts[(a, sp)][0], verdicts[(b, sp)][0]
            if not (ok_a or ok_b):
                continue
            item = {"table": name, "a": a, "b": b, "factor": factor, "species": sp,
                    "ladder": f"{pots[a]['ledger_id']}|{sp}|{factor}",
                    "exposure": "dev_exposed" if (name, sp) in EXPOSED else "clean",
                    "methionine": "MET" in pots[a]["initial_mmol_L"],
                    "refused_arm": None if ok_a and ok_b else (b if ok_a else a),
                    "factor_value": [_factor_value(pots[a], factor), _factor_value(pots[b], factor)]}
            item.update(_measured(name, arms[a].get(sp, {}), arms[b].get(sp, {})))
            out.append(item)
    return out


def _factor_value(pot, factor):
    return pot.get("label_ph") if factor == "label_ph" else (pot["initial_mmol_L"].get(factor) or 0.0)


def _bound(name: str, arm: Mapping[str, Any]) -> Optional[float]:
    if arm.get("nd"):
        return arm.get("bound", ND_BOUND.get(name))
    if arm.get("tr"):
        return TR_BOUND[name][1] if name in TR_BOUND else None
    return None


def sigma_dex(name: str, arm: Mapping[str, Any]) -> float:
    """The log10 scatter of one printed value: the paper's own SD of the reported mean, plus the read error."""
    v = arm["value"]
    if arm.get("sd") is not None:
        rel = arm["sd"] / v / math.sqrt(arm.get("n") or 1.0)
    else:
        rel = STATED_REL_SD[name]
    return math.hypot(rel / math.log(10), arm.get("read_error") or 0.0)


def dof(arm: Mapping[str, Any]) -> float:
    """Degrees of freedom of a printed value's scatter: N - 1 where the row prints SD and N (Grosch), else infinite (a
    stated CV)."""
    return arm["n"] - 1 if arm.get("sd") is not None and arm.get("n") and arm["n"] > 1 else math.inf


def welch_dof(sa: float, da: float, sb: float, db: float) -> float:
    terms = [s ** 4 / d for s, d in ((sa, da), (sb, db)) if math.isfinite(d)]
    return math.inf if not terms else (sa ** 2 + sb ** 2) ** 2 / sum(terms)


def _cdf(z: float, df: float) -> float:
    """Normal, or Student-t where the scatter comes from 2 - 4 trials (v5 review item 4)."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))) if math.isinf(df) else float(student_t.cdf(z, df))


def _measured(name: str, ea: Mapping[str, Any], eb: Mapping[str, Any]) -> Dict[str, Any]:
    """Prereg v6 sec. 2D. A pair's S weight w = 2 F(|m| / sigma_m) - 1, F the normal or (Welch) Student-t CDF: how
    surely the measured sign is real given the papers' own scatter (0 = a coin, 1 = certain). No cut-off: every
    detected pair enters S with its weight."""
    va, vb = ea.get("value"), eb.get("value")
    if va and vb and va > 0 and vb > 0:
        m = math.log10(va / vb)
        sa, sb = sigma_dex(name, ea), sigma_dex(name, eb)
        sig, df = math.hypot(sa, sb), welch_dof(sa, dof(ea), sb, dof(eb))
        return {"kind": "detected", "m": m, "sigma_m": sig, "dof": df if math.isfinite(df) else None,
                "weight": 2 * _cdf(abs(m) / sig, df) - 1}
    if ea.get("pnq") or eb.get("pnq"):
        return {"kind": "excluded", "why": "present, not quantified"}
    for printed, bound_arm, sign in ((ea, eb, 1.0), (eb, ea, -1.0)):
        v = printed.get("value")
        if v and v > 0 and _bound(name, bound_arm):
            m_low = math.log10(v / _bound(name, bound_arm))
            w = max(0.0, 2 * _cdf(m_low / sigma_dex(name, printed), dof(printed)) - 1)
            if w > 0:
                return {"kind": "censored", "m_sign": sign, "bound_ratio": v / _bound(name, bound_arm), "weight": w}
            return {"kind": "excluded", "why": f"printed value below the other arm's limit (ratio {v / _bound(name, bound_arm):.2f})"}
    if (va and not eb) or (vb and not ea):
        return {"kind": "excluded", "why": "one arm has no row and no censoring entry"}
    return {"kind": "excluded", "why": "both arms censored, or n.d. with no limit printed"}


def in_decision_set(it: Mapping[str, Any]) -> bool:
    return it["exposure"] == "clean" and not it["methionine"]


def s_items(q2) -> List[Dict[str, Any]]:
    return [i for i in q2 if in_decision_set(i) and i["kind"] in ("detected", "censored")]


def e_items(q2) -> List[Dict[str, Any]]:
    return [i for i in q2 if in_decision_set(i) and i["kind"] == "detected"]


def measured_sign(it) -> float:
    return math.copysign(1.0, it["m"]) if it["kind"] == "detected" else it["m_sign"]


# --- statistics ----------------------------------------------------------------------------------------------------------


def ladder_boot_counts(ladders: Sequence[str], rng: np.random.Generator) -> Tuple[List[str], np.ndarray]:
    names = sorted(set(ladders))
    draws = rng.integers(0, len(names), size=(N_BOOT, len(names)))
    counts = np.zeros((N_BOOT, len(names)))
    for j in range(len(names)):
        counts[:, j] = (draws == j).sum(axis=1)
    return names, counts


def interval_above_chance(items, hits: np.ndarray, counts: np.ndarray, names: List[str]) -> Tuple[bool, Dict[str, Any]]:
    """The weighted S and its 90 % ladder-bootstrap interval (inverted_cdf); True when the interval lies above 0.5."""
    idx = {n: j for j, n in enumerate(names)}
    L = len(names)
    hit_l, w_l = np.zeros(L), np.zeros(L)
    for it, h in zip(items, hits):
        hit_l[idx[it["ladder"]]] += it["weight"] * h
        w_l[idx[it["ladder"]]] += it["weight"]
    denom = counts @ w_l
    boot = np.where(denom > 0, (counts @ hit_l) / np.where(denom > 0, denom, 1.0), np.nan)
    lo, hi = np.nanpercentile(boot, [5, 95], method="inverted_cdf")
    s_w = float(hit_l.sum() / w_l.sum()) if w_l.sum() > 0 else float("nan")
    return bool(lo > 0.5), {"S": s_w, "ci90": [float(lo), float(hi)]}


def simulate_hits(items, rng: np.random.Generator, measured_values: Mapping[Tuple[str, str], float],
                  agreement: float, one_draw: bool = False) -> np.ndarray:
    """Prereg v6 sec. 5 power model. A presence / absence ladder: each pair agrees with prob `agreement` (with
    `one_draw`, one draw per ladder decides all its pairs: the bound for pairs that share one engine law). An ordered
    (pH) ladder: the engine's latent values walk over the pots in factor order; each adjacent step agrees with the
    measured adjacent step with prob `agreement`, magnitude Exp(1); pair signs follow. A refused arm is always a miss."""
    hits = np.zeros(len(items))
    by_ladder: Dict[str, List[int]] = defaultdict(list)
    for k, it in enumerate(items):
        by_ladder[it["ladder"]].append(k)
    for lad, ks in by_ladder.items():
        first = items[ks[0]]
        if first["factor"] != "label_ph":
            shared = rng.random() < agreement
            for k in ks:
                agree = shared if one_draw else rng.random() < agreement
                hits[k] = 0.0 if items[k]["refused_arm"] else float(agree)
            continue
        pots = sorted({(items[k]["factor_value"][0], items[k]["a"]) for k in ks} |
                      {(items[k]["factor_value"][1], items[k]["b"]) for k in ks})
        latent = {pots[0][1]: 0.0}
        for (_, p0), (_, p1) in zip(pots, pots[1:]):
            step = math.copysign(1.0, measured_values.get((p1, first["species"]), 0.0)
                                 - measured_values.get((p0, first["species"]), 0.0) or 1.0)
            agree = rng.random() < agreement
            latent[p1] = latent[p0] + (step if agree else -step) * rng.exponential()
        for k in ks:
            it = items[k]
            p = latent[it["a"]] - latent[it["b"]]
            hits[k] = 0.0 if it["refused_arm"] else float(math.copysign(1.0, p) == measured_sign(it))
    return hits


def power_line(q2, tables) -> Dict[str, Any]:
    """P(the 90 % interval of the weighted S lies above 0.5) over the agreement grid, both presence-ladder variants; at
    agreement 0.5 it is the false-pass rate. Context for reading the result, not a gate (prereg v6 sec. 5)."""
    items = s_items(q2)
    if not items:
        return {"items": 0, "ladders": 0, "power": 0.0}
    measured: Dict[Tuple[str, str], float] = {}
    for name in SCORED_TABLES:
        for pid, sps in species_arms(tables[name]).items():
            for sp, arm in sps.items():
                v = arm.get("value")
                measured[(pid, sp)] = math.log10(v) if v and v > 0 else -9.0
    names, counts = ladder_boot_counts([i["ladder"] for i in items], np.random.default_rng(SEED))
    curve: Dict[str, Dict[str, float]] = {}
    for variant, one_draw in (("per_pair", False), ("one_draw_per_presence_ladder", True)):
        curve[variant] = {}
        for agreement in AGREEMENT_GRID:
            rng = np.random.default_rng(SEED)
            passes = sum(interval_above_chance(items, simulate_hits(items, rng, measured, agreement, one_draw),
                                               counts, names)[0] for _ in range(N_SIM))
            curve[variant][f"{agreement:.1f}"] = passes / N_SIM
    return {"items": len(items), "ladders": len(names), "simulations": N_SIM, "bootstrap": N_BOOT, "seed": SEED,
            "false_pass_rate": curve["per_pair"]["0.5"], "curve": curve}


# --- score mode (prereg v6 secs. 3 - 6) -----------------------------------------------------------------------------


def engine_value(table: str, pid: str, pot: Mapping[str, Any], species: str, *, clamped: bool) -> Tuple[Optional[float], bool]:
    """(ug/L of `species` in `pid`, pH fixed point unconverged?). The value is None ONLY when the envelope refuses
    (prereg v6 sec. 4: a zero). Any other exception aborts the run: catching it would turn a bug into a 'zero' (the v4
    review's B8)."""
    spec = pot_spec(table, pid, pot, clamped=clamped)
    tg = target_name(species)
    if spec is None or tg is None:
        return None, False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pred = engine.predict(spec, [tg])
    unconverged = any(issubclass(w.category, PhFixedPointNotConverged) for w in caught)
    if not pred.answered:
        return None, unconverged
    return pred.require(tg), unconverged     # keyed by the requested alias, not the engine key


def log_ratio(va: Optional[float], vb: Optional[float]) -> Optional[float]:
    return math.log10(va / vb) if va and vb and va > 0 and vb > 0 else None


def p_of(va: Optional[float], vb: Optional[float], floor: float = 0.0) -> Optional[float]:
    """p, with the zero / refusal rule: None = a one-armed zero or refusal (a miss, error +inf); both arms = 0; a |p| at
    or below the item's numerical floor (sec. 4) = 0."""
    za, zb = not va or va <= 0, not vb or vb <= 0
    if za and zb:
        return 0.0
    if za or zb:
        return None
    p = math.log10(va / vb)
    return 0.0 if abs(p) <= max(floor, FLOAT_FLOOR_DEX) else p


def s_hit(it, p: Optional[float]) -> float:
    return 0.0 if p is None or p == 0.0 else float(math.copysign(1.0, p) == measured_sign(it))


def e_err(it, p: Optional[float]) -> float:
    return math.inf if p is None else abs(p - it["m"])


def _ci(boot: Sequence[float]) -> Tuple[float, float]:
    lo, hi = np.nanpercentile(np.array(boot, dtype=float), [5, 95], method="inverted_cdf")
    return float(lo), float(hi)


def ladder_ci(items, values: Sequence[float], stat, rng) -> Tuple[float, float]:
    """90 % ladder-cluster bootstrap of stat(items, values) over (item, value) pairs."""
    names, counts = ladder_boot_counts([i["ladder"] for i in items], rng)
    groups = defaultdict(list)
    for it, v in zip(items, values):
        groups[it["ladder"]].append((it, v))
    boot = []
    for row in counts:
        pool = [iv for n, c in zip(names, row) for _ in range(int(c)) for iv in groups[n]]
        boot.append(stat(pool) if pool else math.nan)
    return _ci(boot)


def lab_ci(items, values, stat, rng) -> Tuple[float, float]:
    """The degenerate lab-cluster bootstrap (3 labs = 3 tables), shown, never read."""
    tables = sorted({i["table"] for i in items})
    groups = {t: [(i, v) for i, v in zip(items, values) if i["table"] == t] for t in tables}
    boot = []
    for _ in range(N_BOOT):
        pick = rng.integers(0, len(tables), size=len(tables))
        pool = [iv for j in pick for iv in groups[tables[j]]]
        boot.append(stat(pool) if pool else math.nan)
    return _ci(boot)


def _median(xs):
    """The lower median: rank-based, so +/-inf errors never interpolate to NaN (prereg v6 sec. 5)."""
    xs = list(xs)
    return statistics.median_low(xs) if xs else math.nan


def _weighted_share(pool) -> float:
    w = sum(i["weight"] for i, _ in pool)
    return sum(i["weight"] * h for i, h in pool) / w if w > 0 else math.nan


def _median_of(pool) -> float:
    return _median(v for _, v in pool)


def _share(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else math.nan


def item_floor(it, preds, tight) -> float:
    """|p(shipped) - p(tight pH fixed point)|, the item's measured numerical floor (0 when either is not a ratio)."""
    key = lambda pid: (it["table"], pid, it["species"])
    ps, pt = (log_ratio(src[key(it["a"])], src[key(it["b"])]) for src in (preds, tight))
    return abs(ps - pt) if ps is not None and pt is not None else 0.0


def score_block(items: Sequence[Mapping[str, Any]], preds: Mapping[Tuple[str, str, str], Optional[float]],
                tight: Mapping[Tuple[str, str, str], Optional[float]], rng) -> Dict[str, Any]:
    s_set = [i for i in items if i["kind"] in ("detected", "censored")]
    e_set = [i for i in items if i["kind"] == "detected"]
    floors = {id(i): item_floor(i, preds, tight) for i in items}
    ps = {id(i): p_of(preds[(i["table"], i["a"], i["species"])], preds[(i["table"], i["b"], i["species"])],
                      floors[id(i)]) for i in items}
    hits = [s_hit(i, ps[id(i)]) for i in s_set]
    errs = [e_err(i, ps[id(i)]) for i in e_set]
    null_errs = [abs(i["m"]) for i in e_set]
    diffs = [n - e for n, e in zip(null_errs, errs)]
    s_pairs = list(zip(s_set, hits))
    by = lambda pairs, key: {k: [iv for iv in pairs if iv[0][key] == k] for k in sorted({iv[0][key] for iv in pairs})}
    return {
        "S": {"items": len(s_set), "weight_sum": sum(i["weight"] for i in s_set), "share_weighted": _weighted_share(s_pairs),
              "share_unweighted": _share(hits),
              "ci90_ladder": ladder_ci(s_set, hits, _weighted_share, rng) if s_set else None,
              "ci90_lab": lab_ci(s_set, hits, _weighted_share, rng) if s_set else None,
              "median_of_ladders": _median(_weighted_share(v) for v in by(s_pairs, "ladder").values()),
              "per_table": {t: _weighted_share(v) for t, v in by(s_pairs, "table").items()},
              "zero_or_refused_one_arm": sum(ps[id(i)] is None for i in s_set),
              "p_zero": sum(ps[id(i)] == 0.0 for i in s_set)},
        "E": {"items": len(e_set), "engine_median_dex": _median(errs), "no_change_median_dex": _median(null_errs),
              "engine_within_3x_share": _share(e < BAND_DEX for e in errs),
              "no_change_within_3x_share": _share(e < BAND_DEX for e in null_errs),
              "median_difference_dex": _median(diffs),
              "ci90_ladder_difference": ladder_ci(e_set, diffs, _median_of, rng) if e_set else None,
              "ci90_lab_difference": lab_ci(e_set, diffs, _median_of, rng) if e_set else None,
              "per_table_difference_dex": {t: _median_of(v) for t, v in by(list(zip(e_set, diffs)), "table").items()},
              "difference_median_of_ladders": _median(_median_of(v) for v in by(list(zip(e_set, diffs)), "ladder").values()),
              "per_table_engine_median_dex": {t: _median_of(v) for t, v in by(list(zip(e_set, errs)), "table").items()},
              "per_table_no_change_median_dex": {t: _median_of(v) for t, v in by(list(zip(e_set, null_errs)), "table").items()},
              "engine_median_of_ladders": _median(_median_of(v) for v in by(list(zip(e_set, errs)), "ladder").values()),
              "infinite_errors": sum(math.isinf(e) for e in errs)},
        "rows": [{**{k: i[k] for k in ("table", "a", "b", "factor", "species", "ladder", "exposure", "methionine", "kind")},
                  "weight": i.get("weight"), "sigma_m": i.get("sigma_m"), "numerical_floor_dex": floors[id(i)],
                  "m": i.get("m"), "m_sign": measured_sign(i) if i["kind"] in ("detected", "censored") else None,
                  "p": ps[id(i)]} for i in items],
    }


def _reading(ci: Optional[Sequence[float]], null: float, above: str, inside: str, below: str) -> str:
    if not ci:
        return "no items"
    lo, hi = ci
    return above if lo > null else (below if hi < null else inside)


def readings(block: Mapping[str, Any], power: Mapping[str, Any]) -> Dict[str, Any]:
    """Prereg v6 sec. 6: comparisons with their intervals, read by the sentences written there; no pass mark."""
    s, e = block["S"], block["E"]
    return {
        "direction": _reading(s["ci90_ladder"], 0.5,
                              "better than a coin on these held-out directions",
                              "not distinguishable from a coin on these data (see the detection curve for what they could detect)",
                              "worse than a coin on these held-out directions"),
        "log_ratio_error": _reading(e["ci90_ladder_difference"], 0.0,
                                    "smaller log-ratio error than predicting no change",
                                    "not distinguishable from predicting no change on these data",
                                    "larger log-ratio error than predicting no change"),
        "power_curve": power["curve"], "false_pass_rate": power["false_pass_rate"],
        "status": "PROVISIONAL (3 laboratories): a diagnostic; nothing ships, no constant moves, no split changes",
    }


def score_payload(items_record: Mapping[str, Any]) -> Dict[str, Any]:
    tables = {n: load_table(n) for n in SCORED_TABLES}
    q2 = items_record["q2_items"]
    allowed_none = {(i["table"], i["refused_arm"], i["species"]) for i in q2 if i["refused_arm"]}
    keys = sorted({(it["table"], pid, it["species"]) for it in q2 for pid in (it["a"], it["b"])})
    tight_config = Config.from_globals().replace({"ph.fixed_point.tolerance": TIGHT_PH_TOLERANCE})
    out: Dict[str, Any] = {}
    for label, clamped in (("primary", False), ("sensitivity_clamped_ph", True)):
        run = {k: engine_value(k[0], k[1], tables[k[0]]["pots"][k[1]], k[2], clamped=clamped) for k in keys}
        with applied(tight_config):
            tight = {k: engine_value(k[0], k[1], tables[k[0]]["pots"][k[1]], k[2], clamped=clamped)[0] for k in keys}
        preds = {k: v for k, (v, _) in run.items()}
        unexpected = sorted("|".join(k) for k, v in preds.items() if v is None and k not in allowed_none)
        if unexpected and not clamped:
            raise SystemExit(f"refusals the frozen envelope verdicts did not predict: {unexpected}")
        rng = np.random.default_rng(SEED)
        decision = [i for i in q2 if in_decision_set(i)]
        exposed = [i for i in q2 if i["exposure"] == "dev_exposed"]
        methionine = [i for i in q2 if i["methionine"]]
        out[label] = {"decision_set": score_block(decision, preds, tight, rng),
                      "dev_exposed": score_block(exposed, preds, tight, rng),
                      "methionine_all": score_block(methionine, preds, tight, rng),
                      "new_refusals": unexpected,
                      "ph_fixed_point_unconverged": sorted("|".join(k) for k, (_, u) in run.items() if u),
                      "predictions_ug_per_l": {"|".join(k): v for k, v in preds.items()},
                      "predictions_tight_ph_ug_per_l": {"|".join(k): v for k, v in tight.items()}}
    out["readings"] = readings(out["primary"]["decision_set"], items_record["power"])
    return out


# --- modes ---------------------------------------------------------------------------------------------------------------


def items_payload() -> Dict[str, Any]:
    built = build_items()
    q2 = built["q2"]
    tables = built["tables"]
    count_table = []
    for t in SCORED_TABLES:
        mine = [i for i in q2 if i["table"] == t]
        e = [i for i in e_items(q2) if i["table"] == t]
        s = [i for i in s_items(q2) if i["table"] == t]
        count_table.append({
            "table": t, "q2_items": len(mine),
            "exposed": sum(i["exposure"] == "dev_exposed" for i in mine), "methionine": sum(i["methionine"] for i in mine),
            "E_detected": len(e), "E_ladders": len({i["ladder"] for i in e}),
            "S_detected": sum(i["kind"] == "detected" for i in s), "S_censored": sum(i["kind"] == "censored" for i in s),
            "S_weight": round(sum(i["weight"] for i in s), 3),
            "S_ladders": len({i["ladder"] for i in s}), "S_refused_arm": sum(bool(i["refused_arm"]) for i in s)})
    atlas = atlas_coverage()
    return {
        "artifact": "kinetic_core_contrast_benchmark_items",
        "provenance": provenance.provenance_block(
            "kinetic_core_contrast_benchmark_items", generated_by="scripts/generators/score_contrast_benchmark.py",
            inputs=[SPLIT_V2] + [FIT_ROWS_DIR / f"{n}.yml" for n in KEYED_TABLES]),
        "prereg": "results/validation/kinetic_core_contrast_benchmark_score_prereg.md (v6)",
        "config_digest": Config.from_globals().digest(),
        "atlas_table_sha256": atlas_table_sha256(),
        "name_map": built["name_map"],
        "q1_atlas": {"classes": dict(Counter(i["class"] for i in atlas)), "items": atlas},
        "q1_keyed": {"classes": dict(Counter(i["class"] for i in built["keyed"])),
                     "by_table": {t: dict(Counter(i["class"] for i in built["keyed"] if i["table"] == t)) for t in KEYED_TABLES},
                     "pot_pairs": len({(i["table"], i["a"], i["b"]) for i in built["keyed"]}),
                     "items": built["keyed"]},
        "q2_counts": count_table,
        "q2_items": q2,
        "power": power_line(q2, tables),
    }


def render_items(p: Mapping[str, Any]) -> str:
    lines = ["# Contrast benchmark: the frozen item sets (no engine integration)", "",
             f"*Generated by `scripts/generators/score_contrast_benchmark.py --items` for {p['prereg']}. Config digest "
             f"`{p['config_digest'][:16]}`; hand table sha256 `{p['atlas_table_sha256'][:16]}`.*", "",
             "## Q1, atlas family (the 65 v2 TEST contrasts)", "", "| class | cells |", "| --- | --- |"]
    lines += [f"| {k} | {v} |" for k, v in sorted(p["q1_atlas"]["classes"].items(), key=lambda kv: -kv[1])]
    lines += ["", "| cell | errors | class | reason |", "| --- | --- | --- | --- |"]
    lines += [f"| {i['contrast']} | {', '.join(i['errors'])} | {i['class']} | {i['reason'][:160]} |" for i in p["q1_atlas"]["items"]]
    k = p["q1_keyed"]
    lines += ["", f"## Q1, keyed family ({k['pot_pairs']} one-factor pot pairs; items = pot pair x species)", "",
              "| table | " + " | ".join(sorted(k["classes"])) + " |", "| --- |" + " --- |" * len(k["classes"])]
    for t, c in k["by_table"].items():
        lines.append(f"| {t} | " + " | ".join(str(c.get(cl, 0)) for cl in sorted(k["classes"])) + " |")
    lines.append("| **all** | " + " | ".join(f"**{k['classes'][cl]}**" for cl in sorted(k["classes"])) + " |")
    lines += ["", "## Q2 items", "", "All Q2 items, then the decision set (clean, methionine excluded). "
              "Exposed and methionine overlap (Meynier Met FUR). S items are weighted by how surely the measured sign "
              "is real given the papers' own scatter (sum of weights shown).", "",
              "| table | all Q2 items | of which exposed | of which methionine | decision set: E detected / ladders | "
              "S detected | S censored | S weight | S ladders | S with a refused arm |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for c in p["q2_counts"]:
        lines.append(f"| {c['table']} | {c['q2_items']} | {c['exposed']} | {c['methionine']} | {c['E_detected']} / "
                     f"{c['E_ladders']} | {c['S_detected']} | {c['S_censored']} | {c['S_weight']:.2f} | {c['S_ladders']} | "
                     f"{c['S_refused_arm']} |")
    pw = p["power"]
    lines += ["", f"**What these data can detect:** P(the 90 % interval of the weighted S lies above 0.5) on {pw['items']} "
              f"S items in {pw['ladders']} ladders, by per-step agreement q ({pw.get('simulations')} simulations x "
              f"{pw.get('bootstrap')} ladder bootstraps, seed {pw.get('seed')}). At q = 0.5 it is the false-pass rate.", "",
              "| per-step agreement | " + " | ".join(pw["curve"]["per_pair"]) + " |",
              "| --- |" + " --- |" * len(pw["curve"]["per_pair"])]
    for variant, row in pw["curve"].items():
        lines.append(f"| {variant.replace('_', ' ')} | " + " | ".join(f"{v:.3f}" for v in row.values()) + " |")
    lines += ["",
              "## Q2 item list", "", "| table | a | b | factor | species | ladder | exposure | Met | kind | refused arm |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for i in p["q2_items"]:
        kind = i["kind"] + (f" (w {i['weight']:.2f})" if i.get("weight") is not None else "") + (f" ({i['why']})" if i.get("why") else "")
        lines.append(f"| {i['table'].replace('heldout_', '')} | {i['a']} | {i['b']} | {i['factor']} | {i['species']} | "
                     f"{i['ladder']} | {i['exposure']} | {'yes' if i['methionine'] else ''} | {kind} | {i['refused_arm'] or ''} |")
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--items", action="store_true", help="build the item sets and the power line; no integration")
    mode.add_argument("--score", action="store_true", help="integrate and score (after the prereg is frozen)")
    args = ap.parse_args(argv)
    if Config.from_globals().digest() != Config.shipped().digest():
        raise SystemExit("the engine globals are not the shipped Config")
    if args.items:
        payload = items_payload()
        artifact_io.write_artifact(payload, ITEMS_OUT, render=render_items)
        print(f"atlas {payload['q1_atlas']['classes']}; keyed {payload['q1_keyed']['classes']}; "
              f"detection curve {payload['power']['curve']['per_pair']} on {payload['power']['items']} S items")
        return 0
    record = data_access.load_json(ITEMS_OUT)
    rederived = json.loads(artifact_io.dump_json(items_payload()))
    if provenance.strip_volatile(rederived) != provenance.strip_volatile(record):
        raise SystemExit("the items re-derived now differ from the committed items record: re-freeze first")
    if record["atlas_table_sha256"] != atlas_table_sha256():
        raise SystemExit("ATLAS_COVERAGE / ATLAS_SPECS changed since the items were frozen")
    frozen_at = _committed_at(Path(__file__), ITEMS_OUT, ITEMS_OUT.with_suffix(".md"), PREREG)
    payload = {
        "artifact": "kinetic_core_contrast_benchmark_score_v1",
        "provenance": provenance.provenance_block(
            "kinetic_core_contrast_benchmark_score_v1", generated_by="scripts/generators/score_contrast_benchmark.py",
            inputs=[ITEMS_OUT, PREREG] + [FIT_ROWS_DIR / f"{n}.yml" for n in SCORED_TABLES]),
        "prereg": record["prereg"], "config_digest": Config.from_globals().digest(),
        "atlas_table_sha256": record["atlas_table_sha256"], "frozen_at_commit": frozen_at,
        "coverage": {"atlas": record["q1_atlas"]["classes"], "keyed": record["q1_keyed"]["classes"]},
        "power": record["power"],
        **score_payload(record),
    }
    artifact_io.write_artifact(payload, SCORE_OUT, render=render_score)
    print(json.dumps({k: v for k, v in payload["readings"].items() if k != "power_curve"}, indent=1))
    return 0


def _committed_at(*paths: Path) -> str:
    """The commit that froze the generator, the items record and the prereg; refuses uncommitted changes to any."""
    dirty = subprocess.run(["git", "status", "--porcelain", "--", *map(str, paths)], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.strip()
    if dirty:
        raise SystemExit(f"commit the frozen files before scoring:\n{dirty}")
    return subprocess.run(["git", "log", "-1", "--format=%H", "--", *map(str, paths)], cwd=ROOT,
                          capture_output=True, text=True, check=True).stdout.strip()


def render_score(p: Mapping[str, Any]) -> str:
    r = p["readings"]
    ci = lambda x: f"{x[0]:.2f} - {x[1]:.2f}" if x else ""
    lines = ["# Contrast benchmark, scoring step 1: the shipped engine on the clean held-out contrasts", "",
             f"*Generated by `scripts/generators/score_contrast_benchmark.py --score` under {p['prereg']}; items, hand "
             f"table and prereg frozen at `{p['frozen_at_commit'][:10]}`; Config digest `{p['config_digest'][:16]}`. "
             "A diagnostic: comparisons with their intervals, no pass mark (STATUS_AND_PLAN sec. 7).*", "",
             "## Readings (sentences frozen in the prereg, sec. 6)", "",
             f"- **Direction:** {r['direction']}.", f"- **Log-ratio error:** {r['log_ratio_error']}.",
             f"- **Detection curve** (P(interval above chance) by per-step agreement q): "
             + ", ".join(f"q {q}: {v:.2f}" for q, v in r["power_curve"]["per_pair"].items()) + ".",
             f"- **Status:** {r['status']}.", "",
             "## Coverage (Q1)", "", f"- Atlas family (65 v2 TEST cells): {p['coverage']['atlas']}",
             f"- Keyed family (pot pair x species): {p['coverage']['keyed']}", ""]
    for label in ("primary", "sensitivity_clamped_ph"):
        lines += [f"## {label.replace('_', ' ')}", ""]
        if p[label].get("new_refusals"):
            lines += [f"New refusals under this run (not in the frozen verdicts): {len(p[label]['new_refusals'])}.", ""]
        lines += [f"pH fixed point unconverged (warned) in {len(p[label]['ph_fixed_point_unconverged'])} of "
                  f"{len(p[label]['predictions_ug_per_l'])} predictions.", ""]
        lines += ["| set | S items (weight) | S weighted | S 90 % (ladders) | S 90 % (labs, degenerate) | S unweighted | "
                  "per table | E items | engine median | no-change median | difference 90 % (ladders) | "
                  "difference 90 % (labs) | within 3x engine / null | one-arm zero or refusal | p = 0 |", "|" + " --- |" * 15]
        for setname in ("decision_set", "dev_exposed", "methionine_all"):
            S, E = p[label][setname]["S"], p[label][setname]["E"]
            per = ", ".join(f"{t.replace('heldout_', '').split('_')[0]} {v:.2f}" for t, v in S["per_table"].items() if v == v)
            fx = lambda x: f"{x:.2f}" if x == x else ""
            lines.append(f"| {setname.replace('_', ' ')} | {S['items']} ({S['weight_sum']:.2f}) | {fx(S['share_weighted'])} | "
                         f"{ci(S['ci90_ladder'])} | {ci(S['ci90_lab'])} | {fx(S['share_unweighted'])} | {per} | {E['items']} | "
                         f"{fx(E['engine_median_dex'])} dex | {fx(E['no_change_median_dex'])} dex | "
                         f"{ci(E['ci90_ladder_difference'])} | {ci(E['ci90_lab_difference'])} | "
                         f"{fx(E['engine_within_3x_share'])} / {fx(E['no_change_within_3x_share'])} | "
                         f"{S['zero_or_refused_one_arm']} | {S['p_zero']} |")
        lines.append("")
    lines += ["## Items (primary)", "", "| set | table | a | b | factor | species | kind | weight | m | measured sign | p |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    fmt = lambda x: "" if x is None else (f"{x:+.2f}" if isinstance(x, float) else str(x))
    for setname in ("decision_set", "dev_exposed", "methionine_all"):
        for row in p["primary"][setname]["rows"]:
            w = row.get("weight")
            lines.append(f"| {setname.replace('_', ' ')} | {row['table'].replace('heldout_', '')} | {row['a']} | {row['b']} | "
                         f"{row['factor']} | {row['species']} | {row['kind']} | {'' if w is None else f'{w:.2f}'} | "
                         f"{fmt(row['m'])} | {fmt(row['m_sign'])} | "
                         f"{'none (zero or refusal in one arm)' if row['p'] is None else fmt(row['p'])} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
