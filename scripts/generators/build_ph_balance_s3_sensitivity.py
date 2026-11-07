"""
The S3 sensitivity run that fixes which shipped outputs depend on the falsified pH balance (2026-11-07).

WHY. M1b part (i) (results/validation/kinetic_core_m1b_ph_balance_pentose_steps_prereg.md secs. 20 - 21, PR #94)
falsified the shipped charge balance B0 + the shipped centre ledger under S3 on CF-084's controls; B1 (amine survival)
and B13 (B1 + acid supply) are carried, not credited (sec. 21.8.6). STATUS_AND_PLAN sec. 7, S3: "the outputs that
depend on it (those whose sensitivity to the falsified structure's switch exceeds the lane's within-lab term) are
refused, never shown at the incumbent's point". This generator measures that sensitivity on every output path the
engine can produce with the sulfur lane's pH solved from the balance, and writes the list the refusal implements.

WHAT IS RUN. Every recipe below is integrated through ``engine.predict`` (the path ``predict``, ``compare``, the
scorecard, the envelope and the directional panel use) at:

* ``B0``: the shipped engine (B60's constants and balance, ``fixed_point_nodes`` coupling): the incumbent's point;
* ``B0_inline``: the same balance on the inline coupling, a numerics diagnostic only (the rivals run inline);
* every rival point of ``results/validation/ph_balance_s3_rival_points.json`` (prereg sec. 21.10, R3; ``--points``):
  for B1 and for B13, the stored screen draws that pass all four signatures in at least one (A1 form x A2 choice in D)
  cell of sec. 21.1, represented by the passing draw at the minimum and at the maximum of each free coordinate, plus
  the rival's best draw (falsification record). These are ``ph_balance_refusal.rival_points()``, the engine's own.

``predict`` is called with ``balance_refusal=False``: the dependence of an output on the falsified balance can only be
measured at the falsified point. No shipped surface passes that flag.

THE RECIPES: (a) every scored panel bundle on the sulfur lane (the scorecard's and the envelope's rows); (b) every arm
of every directional claim that resolves to the sulfur lane (``core_directional_scores``; its pH levers among them);
(c) a stratified grid of the recipe space the CLI accepts: 11 precursor sets (cysteine- and ARP-containing, TTCA,
thiamine, glutathione, glycine, hexose) x 6 declared buffers (none; phosphate 0.05 / 0.2 / 0.5 M; pyrophosphate 0.2 M;
acetate 0.2 M) x 4 label pH (4.5, 5.5, 6.5, 7.5) x 2 thermal programmes (120 C / 60 min; 145 C / 20 min), plus the
same precursor x pH x thermal cells with a CLAMPED pH (the measured-trajectory input: the balance is not solved, so
the rivals do not apply and the clamped row is a control), and two lipid co-integrations (the lipid lane carries no pH
term). The whole space is run; no subset is needed (see the record's timing).

THE OUTPUTS AND THE TERMS (prereg sec. 21.10; record ``thresholds``). Each output is judged by the uncertainty it is
shipped with (R1, R2):

* an absolute concentration (ug/L): balance-dependent when its spread over B0 and every rival point, log10(max / min),
  exceeds tau = H0's 0.964 dex, the between-laboratory term a level is shipped with (sec. 21.8.6, confirmed); omega
  (0.419 dex) is reported as a sensitivity;
* a ratio between two recipes (``compare``, the contrast scores, the directional claims): the spread of the log ratio
  over the points against omega, the within-laboratory term a relative output is shipped with, and the direction
  refused when the points disagree in sign (S3);
* the cooled and in-situ pH: A1's run offset and per-point shape in quadrature, sqrt(0.45^2 + 0.16^2) = 0.478 pH (sec.
  20.2), in every pot; the five-pot bound 0.265 pH (sec. 21.8.1, a drift term) is reported as a sensitivity.

``--extremes-check`` (part of the full run) samples further passing draws on a set of recorded pots and reports whether
any widens the spread the extremes give (R3's monotonicity premise).

    python scripts/generators/build_ph_balance_s3_sensitivity.py --points     # the rivals' points record
    python scripts/generators/build_ph_balance_s3_sensitivity.py --time-one   # one recipe, all four points, timed
    python scripts/generators/build_ph_balance_s3_sensitivity.py --workers 4  # the full run -> the record
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
import statistics
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
for _p in (str(ROOT), str(ROOT / "scripts" / "generators")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from src import artifact_io, data_access, data_paths, provenance  # noqa: E402

GENERATED_BY = "scripts/generators/build_ph_balance_s3_sensitivity.py"
OUT = data_paths.VALIDATION_DIR / "ph_balance_s3_sensitivity.json"
#: every recipe row and every per-output / per-ratio entry, gzipped (the summary record stays small enough to read)
DETAIL = data_paths.VALIDATION_DIR / "ph_balance_s3_sensitivity_detail.json.gz"
PREREG = data_paths.VALIDATION_DIR / "kinetic_core_m1b_ph_balance_pentose_steps_prereg.md"
FALSIFICATION = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_falsification.json"
SCREEN = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_screen.json"
A2 = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_a2.json"
S3_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_a2_s3.json"

from src.kinetic_core.ph_balance_refusal import SIGN_FLOOR_DEX  # noqa: E402  (a log ratio this small has no sign)

#: a point that runs longer than this is recorded as failed (one directional arm ground for > 20 min on 2026-11-07)
POINT_TIMEOUT_S = 300
#: below this a concentration is zero, not a number (the directional panel's floor, directional.ZERO_FLOOR_UG_PER_L)
ZERO_FLOOR_UG_PER_L = 1e-6

TARGETS: Tuple[str, ...] = (
    "2-methyl-3-furanthiol (MFT)", "2-furfurylthiol (FFT)", "bis(2-methyl-3-furyl) disulfide", "furfural",
    "2-acetylthiazole", "methanethiol", "norfuraneol", "hydrogen sulfide", "5-HMF", "melanoidins",
)

#: the grid (c): precursor sets in mM
PRECURSOR_SETS: Mapping[str, Mapping[str, float]] = {
    "cys50_rib35": {"cysteine": 50.0, "ribose": 35.0},
    "cys10_rib35": {"cysteine": 10.0, "ribose": 35.0},
    "cys50_xyl35": {"cysteine": 50.0, "xylose": 35.0},
    "cys50_glc35": {"cysteine": 50.0, "glucose": 35.0},
    "cys20_rib20_gly20": {"cysteine": 20.0, "ribose": 20.0, "glycine": 20.0},
    "arp5": {"arp": 5.0},
    "arp20_cys20": {"arp": 20.0, "cysteine": 20.0},
    "ttca20": {"ttca": 20.0},
    "thi10_rib35": {"thiamine": 10.0, "ribose": 35.0},
    "thi10": {"thiamine": 10.0},
    "gsh20_rib35": {"glutathione": 20.0, "ribose": 35.0},
}
#: precursor class of each set (the summary's stratum)
PRECURSOR_CLASS: Mapping[str, str] = {
    "cys50_rib35": "cysteine + pentose", "cys10_rib35": "cysteine + pentose", "cys50_xyl35": "cysteine + pentose",
    "cys50_glc35": "cysteine + hexose", "cys20_rib20_gly20": "cysteine + pentose + amine", "arp5": "ARP",
    "arp20_cys20": "ARP + cysteine", "ttca20": "TTCA", "thi10_rib35": "thiamine + pentose", "thi10": "thiamine",
    "gsh20_rib35": "glutathione + pentose",
}
BUFFERS: Mapping[str, Optional[Mapping[str, Any]]] = {
    "none": {"kind": "none"},
    "phosphate_0.05M": {"kind": "phosphate", "phosphate_mol_l": 0.05},
    "phosphate_0.2M": {"kind": "phosphate", "phosphate_mol_l": 0.2},
    "phosphate_0.5M": {"kind": "phosphate", "phosphate_mol_l": 0.5},
    "pyrophosphate_0.2M": {"kind": "pyrophosphate", "pyrophosphate_mol_l": 0.2},
    "acetate_0.2M": {"kind": "acetate", "acetate_mol_l": 0.2},
}
CLAMPED = {"kind": "clamped"}
LABEL_PH: Tuple[float, ...] = (4.5, 5.5, 6.5, 7.5)
THERMAL: Mapping[str, Tuple[float, float]] = {"120C_60min": (120.0, 60.0), "145C_20min": (145.0, 20.0)}
#: the grid's levers, in the order a one-factor pair is looked for
LEVERS: Tuple[str, ...] = ("ph", "buffer", "thermal", "dose", "sugar")
DOSE_PAIRS = (("cys10_rib35", "cys50_rib35"),)
SUGAR_PAIRS = (("cys50_rib35", "cys50_xyl35"), ("cys50_rib35", "cys50_glc35"))


# ----------------------------------------------------------------------------------------------------------------
# Thresholds, from the records (never typed)
# ----------------------------------------------------------------------------------------------------------------
def thresholds() -> Dict[str, Any]:
    """The engine's own terms (``ph_balance_refusal.thresholds``; prereg sec. 21.10, R1 - R2), with their sources."""
    from src.kinetic_core import ph_balance_refusal

    th = dict(ph_balance_refusal.thresholds())
    fals = data_access.load_json(FALSIFICATION)
    th.update({
        "tau_source": ("R1: H0's tau_ABS(all), the between-laboratory term an absolute level is shipped with (prereg "
                       "sec. 21.8.6, confirmed): judges absolute concentration outputs"),
        "omega_source": ("R1: H0's omega(all), the within-laboratory term a relative output is shipped with (STATUS_AND_PLAN "
                         "sec. 7, S1 - S4 preamble): judges ratios, contrasts and directions, with the sign; reported as a "
                         "sensitivity on absolute levels"),
        "ph_level_source": (f"R2: prereg sec. 20.2 (A1), the run-level offset s_off = {float(fals['s_off']):g} and the "
                            f"per-point shape s_shape = {float(fals['s_shape']):g} pH of one laboratory's repeat, in "
                            "quadrature: a pH level's within-laboratory term, in every pot"),
        "ph_five_pot_source": ("sensitivity: sigma_max,5 = 0.2 / Phi^-1((1 + 0.05^(1/5)) / 2), the five-cysteine-pot bound "
                               "of prereg sec. 21.8.1 (a drift term, assumption-laden: sec. 21.9); never judges a level"),
    })
    return th


# ----------------------------------------------------------------------------------------------------------------
# The points (R3): B0, B0 inline (a numerics diagnostic), and each rival's screen-passing draws' extremes + best draw
# ----------------------------------------------------------------------------------------------------------------
POINTS_OUT = data_paths.VALIDATION_DIR / "ph_balance_s3_rival_points.json"
CARRIED = ("B1", "B13")


def passing_masks() -> Dict[str, Any]:
    """Per carried rival, which stored screen draws pass all four signatures in at least one (A1 form x A2 choice in
    D) cell of prereg sec. 21.1 (the union: a draw some defensible noise choice admits), from the stored draws."""
    import numpy as np

    import fit_m1b_ph_balance as FM

    a2 = data_access.load_json(A2)
    choices = [FM.A2Choice(c["name"], tuple(float(x) for x in c["sd"]), float(c["edge"]), bool(c["in_d"]),
                           str(c["basis"])) for c in a2["choices"] if c["in_d"]]
    ens = FM.load_screen_ensembles()
    pinned = ens["B0_pinned"]
    pinned_rise = float(pinned.f["C-Zhou-ARP"][0, -1] - pinned.f["C-Zhou-ARP"][0, 0])
    data = FM.a1_measured()
    out: Dict[str, Any] = {}
    for rival in CARRIED:
        e = ens[rival]
        mask = np.zeros(e.size, dtype=bool)
        for form in FM.A1_FORMS:
            for c in choices:
                mask |= FM.a2_signatures(e, data, form, c, pinned_rise)["all"]
        out[rival] = {"draws": e.draws, "mask": mask,
                      "coordinates": [c.name for c in FM.coordinates(rival, "primary")]}
    return out


def build_points() -> Dict[str, Any]:
    """The points record: per carried rival, the passing draw at the min and at the max of each free coordinate, and
    the best draw of the falsification record (sec. 20.4), de-duplicated by draw index."""
    import numpy as np

    masks = passing_masks()
    best = data_access.load_json(FALSIFICATION)["best_draws"]
    rivals: Dict[str, Any] = {}
    counts: Dict[str, Any] = {}
    for rival, m in masks.items():
        idx = np.flatnonzero(m["mask"])
        chosen: Dict[int, List[str]] = {}
        for coord in m["coordinates"]:
            vals = np.array([m["draws"][i][coord] for i in idx])
            for tag, j in (("min", int(idx[int(np.argmin(vals))])), ("max", int(idx[int(np.argmax(vals))]))):
                chosen.setdefault(j, []).append(f"{tag} {coord}")
        b = int(best[rival]["index"])
        chosen.setdefault(b, []).append("best draw (falsification record)")
        rivals[rival] = [{"name": f"{rival}_d{j}", "index": j, "why": why, "passes": bool(m["mask"][j]),
                          "coordinates": {k: float(v) for k, v in m["draws"][j].items()}}
                         for j, why in sorted(chosen.items())]
        counts[rival] = {"stored_draws": int(len(m["draws"])), "passing": int(idx.size), "points": len(rivals[rival]),
                         "coordinates": m["coordinates"],
                         "passing_box": {c: [float(min(m["draws"][i][c] for i in idx)),
                                             float(max(m["draws"][i][c] for i in idx))] for c in m["coordinates"]}}
    return {"artifact": "ph_balance_s3_rival_points",
            "provenance": provenance.provenance_block("ph_balance_s3_rival_points", generated_by=GENERATED_BY + " --points",
                                                      wave="S3 application: the pH balance (M1b part (i))",
                                                      # the points come from the screen records alone; the prereg
                                                      # reports them (sec. 21.10), so hashing it would be a cycle
                                                      inputs=[FALSIFICATION, SCREEN, A2,
                                                              data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_screen_draws.json.gz"]),
            "definition": ("prereg sec. 21.10, R3: each carried rival is represented by the stored screen draws that pass "
                           "all four signatures in at least one (A1 form x A2 choice in D) cell (sec. 21.1), through the "
                           "passing draw at the minimum and at the maximum of each free coordinate, plus its best draw"),
            "davidek_breakpoints_in_situ": data_access.load_json(SCREEN)["davidek_breakpoints_in_situ"],
            "counts": counts, "rivals": rivals}


def point_configs() -> Dict[str, Tuple[Dict[str, Any], Any]]:
    """{point: (Config changes, PhDrift or None)}; None drift = the shipped balance constants. The rivals' points are
    the engine's own (``ph_balance_refusal.rival_points``), so the record audits the rule the engine applies."""
    from src.kinetic_core import ph_balance_refusal

    out: Dict[str, Tuple[Dict[str, Any], Any]] = {"B0": ({}, None), "B0_inline": ({"ph.coupling": "inline"}, None)}
    out.update(ph_balance_refusal.rival_points())
    return out


def deciding(row: Optional[Mapping[str, Any]] = None) -> Tuple[str, ...]:
    """The points whose spread decides: the row's own (B0 and every rival point it ran; B0_inline is a numerics
    diagnostic), else the engine's."""
    if row is not None and row.get("points"):
        return tuple(p for p in row["points"] if p != "B0_inline")
    from src.kinetic_core import ph_balance_refusal

    return ph_balance_refusal.all_points()


def point_coordinates() -> Dict[str, Any]:
    rec = data_access.load_json(POINTS_OUT)
    return {"B0": "shipped (B60 balance constants, fixed_point_nodes)", "B0_inline": "shipped constants, inline",
            "rivals": {r: [{k: p[k] for k in ("name", "index", "why", "passes")} for p in pts]
                       for r, pts in rec["rivals"].items()}, "counts": rec["counts"]}


# ----------------------------------------------------------------------------------------------------------------
# Recipes
# ----------------------------------------------------------------------------------------------------------------
@dataclass
class Recipe:
    rid: str
    source: str                        # "panel" | "directional" | "grid" | "toy"
    spec: Any                          # engine.FormulationSpec
    targets: Tuple[str, ...]
    meta: Dict[str, Any] = field(default_factory=dict)


def _cli_spec(name: str, precursors: Mapping[str, float], temp_c: float, time_min: float, ph: float,
              buffer: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    spec: Dict[str, Any] = {"name": name, "precursors": dict(precursors), "temp_C": temp_c, "time_min": time_min,
                            "ph": ph}
    if buffer is not None:
        spec["buffer"] = dict(buffer, source="S3 sensitivity grid")
    return spec


def grid_recipes() -> List[Recipe]:
    from src.comparative_cli import spec_to_core

    out: List[Recipe] = []
    buffers = dict(BUFFERS, clamped=CLAMPED)
    for pset, buf, ph, th in itertools.product(PRECURSOR_SETS, buffers, LABEL_PH, THERMAL):
        temp_c, time_min = THERMAL[th]
        rid = f"grid:{pset}|{buf}|pH{ph:g}|{th}"
        spec = spec_to_core(_cli_spec(rid, PRECURSOR_SETS[pset], temp_c, time_min, ph, buffers[buf]))
        out.append(Recipe(rid, "grid", spec, TARGETS, {"pset": pset, "precursor_class": PRECURSOR_CLASS[pset],
                                                      "buffer": buf, "ph": ph, "thermal": th}))
    for name, extra in (("lipid:cys50_rib35+pea", {"pea protein isolate": 1.0}),):
        for buf in ("none", "phosphate_0.2M"):
            rid = f"{name}|{buf}|pH6.5|145C_20min"
            spec = spec_to_core(_cli_spec(rid, dict(PRECURSOR_SETS["cys50_rib35"], **extra), 145.0, 20.0, 6.5,
                                          BUFFERS[buf]))
            out.append(Recipe(rid, "grid", spec, TARGETS + ("hexanal", "pentane"),
                              {"pset": "cys50_rib35+lipid", "precursor_class": "lipid co-integration",
                               "buffer": buf, "ph": 6.5, "thermal": "145C_20min"}))
    return out


def panel_recipes() -> List[Recipe]:
    from src.kinetic_core import engine, panel

    out: List[Recipe] = []
    scored, _ = panel.panel_bundles()
    for path, tag in scored:
        bench = panel.load_bundle(path)
        spec = panel.core_spec(bench)
        targets = tuple(panel.bundle_targets(bench))
        decl = engine.declare_envelope(spec, list(targets))
        if engine.SULFUR not in (decl.lanes or (decl.lane,)):
            continue
        buffer = spec.process.buffer
        out.append(Recipe(f"panel:{bench.get('benchmark_id') or path.stem}", "panel", spec, targets,
                          {"bench_file": data_paths.rel(path), "panel": tag,
                           "study": panel.study_key(bench, str(path)),
                           "buffer": "clamped" if (buffer is not None and buffer.is_clamped)
                           else (buffer.kind if buffer is not None else "undeclared"),
                           "ph": float(spec.process.ph)}))
    return out


def directional_recipes() -> Tuple[List[Recipe], List[Dict[str, Any]]]:
    """Every arm of every runnable directional claim that resolves to the sulfur lane, and the claims."""
    from src.comparative_cli import spec_to_core
    from src.kinetic_core import directional, engine

    panel_doc = directional.load_panel()
    recipes: List[Recipe] = []
    claims: List[Dict[str, Any]] = []
    for claim in panel_doc["panel"]:
        arms = directional._arms(claim)
        observables = [str(o) for o in claim.get("observables") or []]
        if not arms or not observables:
            continue
        sweeps = [None]
        unstated = claim.get("unstated_inputs") or {}
        if unstated:
            (key, values), = list(unstated.items())
            sweeps = [(key, float(v)) for v in values]
        for sweep in sweeps:
            c = directional._with_input(claim, sweep[0], sweep[1]) if sweep else claim
            ids = []
            sulfur = False
            for i, (label, system) in enumerate(directional._arms(c)):
                spec = spec_to_core(directional._spec_from_system(label, system))
                decl = engine.declare_envelope(spec, observables)
                sulfur = sulfur or engine.SULFUR in (decl.lanes or (decl.lane,))
                suffix = f"|{sweep[0]}={sweep[1]:g}" if sweep else ""
                rid = f"directional:{claim['claim_id']}{suffix}|arm{i}:{label}"
                ids.append(rid)
                recipes.append(Recipe(rid, "directional", spec, tuple(observables),
                                      {"claim_id": str(claim["claim_id"]), "arm": i,
                                       "category": str(claim.get("category"))}))
            if not sulfur:
                recipes = recipes[: len(recipes) - len(ids)]
                continue
            claims.append({"claim_id": str(claim["claim_id"]), "category": str(claim.get("category")),
                           "claim_type": str(claim.get("claim_type")), "expected": str(claim.get("expected_relation")),
                           "observables": observables, "sweep": list(sweep) if sweep else None, "arms": ids})
    return recipes, claims


# ----------------------------------------------------------------------------------------------------------------
# One recipe at the four points
# ----------------------------------------------------------------------------------------------------------------
def _answerable_targets(spec, targets: Sequence[str]) -> List[str]:
    """The requested targets the engine answers, alone or together (as ``predict_core`` splits them)."""
    from src.kinetic_core import engine

    if engine.declare_envelope(spec, list(targets)).is_answerable:
        return list(targets)
    return [t for t in targets if engine.declare_envelope(spec, [t]).is_answerable]


def _beta_t0(spec, declaration) -> Optional[float]:
    """The pot's buffer capacity at its label pH, 25 C, t = 0, I = 0 (mmol/L per pH): a stratification variable."""
    from src.kinetic_core import engine, ph_state

    buffer = spec.process.buffer if spec.process.buffer is not None else ph_state.DEFAULT_BUFFER
    if buffer.is_clamped:
        return None
    try:
        return float(ph_state.buffer_capacity_mmol_per_ph(dict(declaration.mapped_precursors), 298.15,
                                                          float(spec.process.ph), engine.core_ph_drift(), buffer))
    except Exception:  # noqa: BLE001 - a stratification variable; recorded as missing
        return None


def evaluate_recipe(recipe: Recipe, points: Optional[Mapping[str, Tuple[Dict[str, Any], Any]]] = None,
                    predict_fn: Optional[Callable[..., Any]] = None) -> Dict[str, Any]:
    """Integrate one recipe at every point; returns per point the concentrations and the pH, or the failure.

    ``predict_fn(spec, targets, changes, drift)`` replaces the engine (the unit test's toy evaluator)."""
    from src.kinetic_core import config, engine

    points = points if points is not None else point_configs()
    t0 = time.time()
    targets = _answerable_targets(recipe.spec, recipe.targets) if predict_fn is None else list(recipe.targets)
    decl = engine.declare_envelope(recipe.spec, targets or list(recipe.targets))
    buffer = recipe.spec.process.buffer
    clamped = bool(buffer is not None and buffer.is_clamped)
    lanes = list(decl.lanes or ((decl.lane,) if decl.lane else ()))
    row: Dict[str, Any] = {"rid": recipe.rid, "source": recipe.source, "meta": dict(recipe.meta), "lanes": lanes,
                           "clamped": clamped, "declared_buffer": bool(buffer is not None and buffer.kind != "none"),
                           "targets": targets, "absolute_refused_h0": list(decl.absolute_refused_targets),
                           "unidentified_hexose": list(engine.unidentified_routes(decl.mapped_precursors,
                                                                                  decl.mapped_targets)),
                           "beta_t0_mmol_per_ph": _beta_t0(recipe.spec, decl), "points": {}}
    if not targets:
        row["refused_before_s3"] = list(decl.reasons)
        row["seconds"] = time.time() - t0
        return row
    rival_failed = False
    for name, (changes, drift) in points.items():
        if clamped and name != "B0":
            # the balance is not solved on a clamped pot, and the rivals run only on the inline dynamic pH
            row["points"][name] = {"not_applicable": "clamped pH: the balance is not solved"}
            continue
        if rival_failed and name not in ("B0", "B0_inline"):
            # the engine refuses the pot at its first failed rival point; the rest would only repeat the budget
            row["points"][name] = {"failed": "skipped: an earlier rival point failed (the engine refuses at the first)"}
            continue
        t = time.time()
        import signal

        def _timeout(signum, frame):
            raise TimeoutError(f"point ran longer than {POINT_TIMEOUT_S} s")

        previous = signal.signal(signal.SIGALRM, _timeout) if predict_fn is None else None
        if predict_fn is None:
            signal.alarm(POINT_TIMEOUT_S)
        try:
            if predict_fn is not None:
                result = predict_fn(recipe.spec, targets, changes, drift)
            else:
                with config.applied(config.Config.shipped().replace(changes)):
                    run = engine.predict(recipe.spec, targets,
                                         draw=engine.CoreDraw(ph_drift=drift) if drift is not None else None,
                                         balance_refusal=False)
                seg = (run.run_metadata.get("segments") or [{}])[-1]
                result = {"answered": bool(run.answered),
                          "concentrations_ug_per_l": {k: float(v) for k, v in run.concentrations_ug_per_l.items()},
                          "ph_cooled_end": seg.get("ph_cooled_end"), "ph_in_situ_end": seg.get("ph_in_situ_end"),
                          "ph_mode": seg.get("ph_mode")}
        except Exception as exc:  # noqa: BLE001 - a failed pot run is recorded, never dropped silently
            result = {"failed": f"{type(exc).__name__}: {exc}"[:300]}
        finally:
            if predict_fn is None:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, previous)
        result["seconds"] = round(time.time() - t, 3)
        row["points"][name] = result
        if result.get("failed") and name not in ("B0", "B0_inline"):
            rival_failed = True
    row["seconds"] = round(time.time() - t0, 3)
    return row


# ----------------------------------------------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------------------------------------------
def _values(row: Mapping[str, Any], output: str, points: Optional[Sequence[str]] = None) -> Optional[List[float]]:
    points = deciding(row) if points is None else points
    vals = []
    for p in points:
        r = (row.get("points") or {}).get(p) or {}
        if r.get("failed") or r.get("not_applicable") or not r.get("answered", False):
            return None
        if output in ("ph_cooled_end", "ph_in_situ_end"):
            v = r.get(output)
        else:
            v = (r.get("concentrations_ug_per_l") or {}).get(output)
        if v is None:
            return None
        vals.append(float(v))
    return vals


def log_spread(vals: Sequence[float]) -> float:
    """log10(max / min) over positive values; inf when one is zero and another is not; 0 when all are zero."""
    pos = [v for v in vals if v > ZERO_FLOOR_UG_PER_L]
    if not pos:
        return 0.0
    if len(pos) < len(vals):
        return math.inf
    return math.log10(max(pos) / min(pos))


def ph_tolerance(row: Mapping[str, Any], th: Mapping[str, Any]) -> float:
    """R2: the pH-level term, buffered or not."""
    return float(th["ph_level"])


def output_entries(rows: Sequence[Mapping[str, Any]], th: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """One entry per (recipe, output): its values at the points, the spread, the tolerance and the verdict."""
    out: List[Dict[str, Any]] = []
    for row in rows:
        if not row.get("targets"):
            continue
        if row.get("clamped"):
            for t in row["targets"]:
                out.append({"rid": row["rid"], "source": row["source"], "output": t, "kind": "absolute",
                            "dependent": False, "spread": 0.0, "tolerance": th["tau_dex"],
                            "why": "clamped pH: the balance is not solved (the measured-trajectory path)"})
            continue
        b0_inline = {}
        for t in list(row["targets"]) + ["ph_cooled_end", "ph_in_situ_end"]:
            vals = _values(row, t)
            kind = "ph" if t.startswith("ph_") else "absolute"
            entry: Dict[str, Any] = {"rid": row["rid"], "source": row["source"], "output": t, "kind": kind,
                                     "meta": row.get("meta"), "lanes": row.get("lanes")}
            if vals is None:
                entry.update(dependent=None, why="not computed at every point (failed, refused or not answered)")
                out.append(entry)
                continue
            if kind == "ph":
                spread = max(vals) - min(vals)
                tol = ph_tolerance(row, th)
                entry["dependent_sensitivity"] = bool(spread > th["ph_five_pot_sensitivity"])
            else:
                spread = log_spread(vals)
                tol = th["tau_dex"]
                entry["dependent_sensitivity"] = bool(spread > th["omega_dex"])
                entry["h0_absolute_refused"] = t in (row.get("absolute_refused_h0") or [])
                entry["hexose_floor"] = bool(row.get("unidentified_hexose"))
            inline = _values(row, t, ("B0", "B0_inline"))
            if inline is not None:
                b0_inline = (abs(inline[0] - inline[1]) if kind == "ph" else log_spread(inline))
                entry["b0_vs_inline"] = b0_inline
            entry.update(values=dict(zip(deciding(row), vals)), spread=spread, tolerance=tol,
                         dependent=bool(spread > tol))
            out.append(entry)
    return out


def _grid_index(rows: Sequence[Mapping[str, Any]]) -> Dict[Tuple[str, str, float, str], Mapping[str, Any]]:
    return {(r["meta"]["pset"], r["meta"]["buffer"], float(r["meta"]["ph"]), r["meta"]["thermal"]): r
            for r in rows if r["source"] == "grid" and r["meta"].get("pset") in PRECURSOR_SETS}


def lever_pairs(rows: Sequence[Mapping[str, Any]]) -> List[Tuple[str, Mapping[str, Any], Mapping[str, Any]]]:
    """(lever, row_a, row_b) for grid recipes that differ in exactly one lever (adjacent levels)."""
    idx = _grid_index(rows)
    pairs = []
    buffers = list(BUFFERS) + ["clamped"]
    for (pset, buf, ph, th), a in idx.items():
        i = LABEL_PH.index(ph)
        if i + 1 < len(LABEL_PH) and (pset, buf, LABEL_PH[i + 1], th) in idx:
            pairs.append(("ph", a, idx[(pset, buf, LABEL_PH[i + 1], th)]))
        if buf == "none":
            for other in buffers[1:]:
                if (pset, other, ph, th) in idx:
                    pairs.append(("buffer", a, idx[(pset, other, ph, th)]))
        if th == "120C_60min" and (pset, buf, ph, "145C_20min") in idx:
            pairs.append(("thermal", a, idx[(pset, buf, ph, "145C_20min")]))
        for lever, table in (("dose", DOSE_PAIRS), ("sugar", SUGAR_PAIRS)):
            for p, q in table:
                if pset == p and (q, buf, ph, th) in idx:
                    pairs.append((lever, a, idx[(q, buf, ph, th)]))
    return pairs


def _log_ratio(a: float, b: float) -> Optional[float]:
    if a <= ZERO_FLOOR_UG_PER_L or b <= ZERO_FLOOR_UG_PER_L:
        return None
    return math.log10(a / b)


def ratio_entry(lever: str, a: Mapping[str, Any], b: Mapping[str, Any], output: str,
                th: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    """The log ratio a / b of one output at every deciding point, its spread and its sign verdict."""
    clamped = a.get("clamped") and b.get("clamped")
    pts = ("B0",) if clamped else deciding(a)
    va, vb = _values(a, output, pts), _values(b, output, pts)
    if va is None or vb is None:
        return None
    lr = [_log_ratio(x, y) for x, y in zip(va, vb)]
    entry: Dict[str, Any] = {"lever": lever, "a": a["rid"], "b": b["rid"], "output": output,
                             "clamped_both": bool(clamped)}
    if any(v is None for v in lr):
        zero_pattern = [x <= ZERO_FLOOR_UG_PER_L for x in va] != [x <= ZERO_FLOOR_UG_PER_L for x in vb]
        entry.update(log_ratio=None, dependent=None, sign_flip=None,
                     why="a zero arm at some point: no ratio" + (" (the zero moves with the balance)"
                                                                 if zero_pattern else ""))
        return entry
    spread = max(lr) - min(lr)
    nonzero = [v for v in lr if abs(v) > SIGN_FLOOR_DEX]
    flip = bool(nonzero) and (min(nonzero) < 0.0 < max(nonzero))
    entry.update(log_ratio=dict(zip(pts, lr)), spread=spread, tolerance=th["omega_dex"],
                 sign_flip=flip, dependent=bool(spread > th["omega_dex"] or flip))
    if clamped:
        entry.update(dependent=False, sign_flip=False, spread=0.0, why="both arms clamped: the balance is not solved")
    return entry


def ratio_entries(rows: Sequence[Mapping[str, Any]], th: Mapping[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for lever, a, b in lever_pairs(rows):
        for t in sorted(set(a.get("targets") or ()) & set(b.get("targets") or ())):
            e = ratio_entry(lever, a, b, t, th)
            if e is not None:
                e["precursor_class"] = a["meta"].get("precursor_class")
                e["buffer_a"], e["buffer_b"] = a["meta"].get("buffer"), b["meta"].get("buffer")
                out.append(e)
    # the panel's within-study pairs (what contrast_scores and compare score on panel pots)
    by_study: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
    for r in rows:
        if r["source"] == "panel":
            by_study[str(r["meta"].get("study"))].append(r)
    for study, members in sorted(by_study.items()):
        for a, b in itertools.combinations(sorted(members, key=lambda r: r["rid"]), 2):
            for t in sorted(set(a.get("targets") or ()) & set(b.get("targets") or ())):
                e = ratio_entry("within_study", a, b, t, th)
                if e is not None:
                    e["study"] = study
                    out.append(e)
    return out


def claim_entries(rows: Sequence[Mapping[str, Any]], claims: Sequence[Mapping[str, Any]],
                  th: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Per directional claim on the sulfur lane: the direction of each step at every point, and whether it flips."""
    by_id = {r["rid"]: r for r in rows}
    out = []
    for c in claims:
        arm_rows = [by_id.get(i) for i in c["arms"]]
        entry: Dict[str, Any] = {k: c[k] for k in ("claim_id", "category", "claim_type", "expected", "sweep")}
        entry["observables"] = c["observables"]
        if any(r is None for r in arm_rows):
            entry.update(dependent=None, why="an arm was not run")
            out.append(entry)
            continue
        steps_by_point: Dict[str, List[Optional[float]]] = {}
        if c["claim_type"] == "ranking":
            row = arm_rows[0]
            for p in deciding(arm_rows[0]):
                vals = [(_values(row, o, (p,)) or [None])[0] for o in c["observables"]]
                steps_by_point[p] = ([None] * (len(vals) - 1) if any(v is None for v in vals) else
                                     [_log_ratio(vals[i + 1], vals[i]) for i in range(len(vals) - 1)])
        else:
            obs = c["observables"][0]
            for p in deciding(arm_rows[0]):
                vals = [(_values(r, obs, (p,)) or [None])[0] for r in arm_rows]
                steps_by_point[p] = ([None] * (len(vals) - 1) if any(v is None for v in vals) else
                                     [_log_ratio(vals[i + 1], vals[i]) for i in range(len(vals) - 1)])
        entry["log_steps"] = steps_by_point
        clamped = all(r.get("clamped") for r in arm_rows)
        if clamped:
            entry.update(dependent=False, sign_flip=False, why="every arm clamped: the balance is not solved")
            out.append(entry)
            continue
        if any(s is None for p in deciding(arm_rows[0]) for s in steps_by_point[p]):
            entry.update(dependent=None, sign_flip=None, why="a step is undefined at some point (zero, refused or failed)")
            out.append(entry)
            continue
        n = len(steps_by_point["B0"])
        flips = []
        spreads = []
        for i in range(n):
            vals = [steps_by_point[p][i] for p in deciding(arm_rows[0])]
            nonzero = [v for v in vals if abs(v) > SIGN_FLOOR_DEX]
            flips.append(bool(nonzero) and min(nonzero) < 0.0 < max(nonzero))
            spreads.append(max(vals) - min(vals))
        entry.update(sign_flip=any(flips), max_step_spread=max(spreads),
                     dependent=bool(any(flips) or max(spreads) > th["omega_dex"]))
        out.append(entry)
    return out


def _summarise(entries: Sequence[Mapping[str, Any]], key: Callable[[Mapping[str, Any]], Tuple]) -> List[Dict[str, Any]]:
    groups: Dict[Tuple, List[Mapping[str, Any]]] = defaultdict(list)
    for e in entries:
        groups[key(e)].append(e)
    out = []
    for k, es in sorted(groups.items(), key=lambda kv: tuple(str(x) for x in kv[0])):
        decided = [e for e in es if e.get("dependent") is not None]
        finite = [float(e["spread"]) for e in decided if e.get("spread") is not None and math.isfinite(e["spread"])]
        out.append({"class": list(k), "n": len(es), "n_decided": len(decided),
                    "n_dependent": sum(1 for e in decided if e["dependent"]),
                    "n_dependent_sensitivity": sum(1 for e in decided if e.get("dependent_sensitivity")),
                    "n_sign_flip": sum(1 for e in decided if e.get("sign_flip")),
                    "n_infinite_spread": sum(1 for e in decided if e.get("spread") == math.inf),
                    "median_spread": statistics.median(finite) if finite else None,
                    "min_spread": min(finite) if finite else None,
                    "max_spread": max(finite) if finite else None})
    return out


def analyse(rows: Sequence[Mapping[str, Any]], claims: Sequence[Mapping[str, Any]],
            th: Mapping[str, Any]) -> Dict[str, Any]:
    outputs = output_entries(rows, th)
    ratios = ratio_entries(rows, th)
    directional = claim_entries(rows, claims, th)
    dyn = [e for e in outputs if not str(e.get("why", "")).startswith("clamped")]

    def stratum(e):
        meta = e.get("meta") or {}
        return (e["kind"], e["output"], meta.get("buffer", "?"))

    summary = {
        "outputs_by_kind": _summarise(dyn, lambda e: (e["kind"],)),
        "outputs_by_output": _summarise(dyn, lambda e: (e["kind"], e["output"])),
        "outputs_by_output_and_buffer": _summarise(dyn, stratum),
        "outputs_by_source": _summarise(dyn, lambda e: (e["source"], e["kind"])),
        "ratios_by_lever": _summarise(ratios, lambda e: (e["lever"],)),
        "ratios_by_lever_and_output": _summarise(ratios, lambda e: (e["lever"], e["output"])),
        "directional_by_category": _summarise(directional, lambda e: (e["category"],)),
        "clamped_controls": {"n_outputs": sum(1 for e in outputs if str(e.get("why", "")).startswith("clamped")),
                             "all_independent": all(not e["dependent"] for e in outputs
                                                    if str(e.get("why", "")).startswith("clamped"))},
    }
    nondep = [e for e in dyn if e.get("dependent") is False]
    summary["non_dependent_dynamic_outputs"] = [
        {k: e.get(k) for k in ("rid", "output", "kind", "spread", "tolerance")} for e in nondep]
    summary["undecided_dynamic_outputs"] = sum(1 for e in dyn if e.get("dependent") is None)
    return {"outputs": outputs, "ratios": ratios, "directional": directional, "summary": summary}


# ----------------------------------------------------------------------------------------------------------------
# Record
# ----------------------------------------------------------------------------------------------------------------
def _fmt(x: Any, places: int = 3) -> str:
    if x is None:
        return "-"
    if isinstance(x, float):
        return "inf" if x == math.inf else f"{x:.{places}f}"
    return str(x)


def render(payload: Mapping[str, Any]) -> str:
    th = payload["thresholds"]
    s = payload["summary"]
    lines = [
        "# The S3 sensitivity of shipped outputs to the falsified pH balance",
        "",
        f"Generated by `{GENERATED_BY}`. Applies STATUS_AND_PLAN sec. 7, S3's refusal clause, to M1b part (i)'s "
        "falsification of B0 + the shipped centre ledger (prereg secs. 20 - 21; CF-084, CF-134).",
        "",
        "## Thresholds",
        "",
        f"* absolute concentrations: tau = {th['tau_dex']:.3f} dex ({th['tau_source']})",
        f"* ratios, contrasts, directions: omega = {th['omega_dex']:.3f} dex and the sign ({th['omega_source']})",
        f"* pH levels: {th['ph_level']:.3f} pH ({th['ph_level_source']})",
        f"* sensitivities: omega on absolute levels; {th['ph_five_pot_sensitivity']:.3f} pH on the pH "
        f"({th['ph_five_pot_source']})",
        "",
        "## Points",
        "",
    ]
    for k, v in payload["points"].items():
        lines.append(f"* `{k}`: {v}")
    t = payload["timing"]
    lines += ["", "## Cost", "",
              f"* recipes run: {t['recipes']}; workers {t['workers']}; wall {t['wall_seconds'] / 3600.0:.2f} h; "
              f"projected before the run: {_fmt(t.get('projected_hours'), 2)} h",
              f"* seconds per point, median: " + ", ".join(f"{k} {_fmt(v, 2)}" for k, v in t["median_seconds"].items()),
              f"* failed point runs: {payload['failed_point_runs']}", ""]

    def table(title, rows, head):
        lines.extend([f"## {title}", "", "| " + " | ".join(head) + " | n | decided | dependent | dependent at the "
                      "sensitivity term | sign flips | median spread | min | max |", "|" + "---|" * (len(head) + 8)])
        for r in rows:
            lines.append("| " + " | ".join(str(c) for c in r["class"]) + f" | {r['n']} | {r['n_decided']} | "
                         f"{r['n_dependent']} | {r.get('n_dependent_sensitivity', '-')} | {r['n_sign_flip']} | "
                         f"{_fmt(r['median_spread'])} | "
                         f"{_fmt(r['min_spread'])} | {_fmt(r['max_spread'])} |")
        lines.append("")

    table("Outputs by kind (dynamic-pH sulfur-lane recipes)", s["outputs_by_kind"], ["kind"])
    table("Outputs by output", s["outputs_by_output"], ["kind", "output"])
    table("Outputs by output and buffer", s["outputs_by_output_and_buffer"], ["kind", "output", "buffer"])
    table("Outputs by recipe source", s["outputs_by_source"], ["source", "kind"])
    table("Ratios by lever (compare / contrasts)", s["ratios_by_lever"], ["lever"])
    table("Ratios by lever and output", s["ratios_by_lever_and_output"], ["lever", "output"])
    table("Directional claims on the sulfur lane, by category", s["directional_by_category"], ["category"])
    lines += ["## Clamped controls", "",
              f"{s['clamped_controls']['n_outputs']} outputs on clamped-pH recipes; all independent of the balance by "
              f"construction: {s['clamped_controls']['all_independent']}.", "",
              "## Dynamic-pH outputs below the threshold", "",
              f"{len(s['non_dependent_dynamic_outputs'])} of the decided dynamic-pH outputs; "
              f"{s['undecided_dynamic_outputs']} not decided (a point failed or refused).", ""]
    x = payload.get("extremes_check") or {}
    if x:
        lines += ["## R3: do the extremes bound the passing draws?", "", x.get("summary", ""), ""]
    dn = payload.get("draws_needed") or {}
    if dn:
        lines += ["## R3: how many passing draws bound the band", "", dn.get("summary", ""), "",
                  "| draws per rival | share within the capture | share whose verdict agrees | (output, subset) pairs |",
                  "|---|---|---|---|"]
        lines += [f"| {r['draws_per_rival']} | {r['share_within_capture']:.3f} | {r['share_verdict_agrees']:.3f} | {r['n']} |"
                  for r in dn.get("table", [])]
        lines.append("")
    lines += ["## The list (the classes the engine's live rule refuses, read off the run)", ""] + [f"* {x}" for x in payload["the_list"]] + [""]
    return "\n".join(lines)


def the_list(analysis: Mapping[str, Any]) -> List[str]:
    """The classes the refusal implements, stated from the summary (no number typed)."""
    s = analysis["summary"]
    out = []
    for r in s["outputs_by_output"]:
        kind, output = r["class"]
        out.append(f"{kind} {output}: {r['n_dependent']} of {r['n_decided']} decided dynamic-pH recipes "
                   f"balance-dependent (spread median {_fmt(r['median_spread'])}, max {_fmt(r['max_spread'])})")
    for r in s["ratios_by_lever"]:
        out.append(f"ratio, lever {r['class'][0]}: {r['n_dependent']} of {r['n_decided']} dependent, "
                   f"{r['n_sign_flip']} direction flips")
    for r in s["directional_by_category"]:
        out.append(f"directional claims, {r['class'][0]}: {r['n_dependent']} of {r['n_decided']} dependent, "
                   f"{r['n_sign_flip']} direction flips")
    return out


#: the recorded pots the extremes check samples (one per precursor class, buffered and unbuffered)
EXTREMES_RIDS: Tuple[str, ...] = (
    "grid:cys50_rib35|none|pH5.5|120C_60min", "grid:cys50_rib35|pyrophosphate_0.2M|pH4.5|145C_20min",
    "grid:arp5|none|pH6.5|120C_60min", "grid:arp20_cys20|phosphate_0.2M|pH5.5|145C_20min",
    "grid:thi10_rib35|phosphate_0.05M|pH6.5|120C_60min", "grid:cys20_rib20_gly20|acetate_0.2M|pH7.5|145C_20min",
    "grid:gsh20_rib35|none|pH4.5|120C_60min", "grid:ttca20|phosphate_0.5M|pH5.5|120C_60min")
EXTREMES_SAMPLES = 24
EXTREMES_SEED = 20261107


def _extremes_task(item) -> Dict[str, Any]:
    """One recorded pot at a set of extra passing draws (the item carries their configs)."""
    from src.kinetic_core import config, engine

    recipe, extra = item
    targets = _answerable_targets(recipe.spec, recipe.targets)
    base = config.Config.from_globals()
    out: Dict[str, Any] = {}
    for name, (changes, drift) in extra.items():
        try:
            with config.applied(base.replace(changes)):
                run = engine.predict(recipe.spec, targets, draw=engine.CoreDraw(ph_drift=drift), balance_refusal=False)
            seg = (run.run_metadata.get("segments") or [{}])[-1]
            out[name] = {"concentrations_ug_per_l": dict(run.concentrations_ug_per_l),
                         "ph_cooled_end": seg.get("ph_cooled_end")}
        except Exception as exc:  # noqa: BLE001
            out[name] = {"failed": f"{type(exc).__name__}: {exc}"[:200]}
    return {"rid": recipe.rid, "targets": targets, "extra": out}


def extremes_check(rows: Sequence[Mapping[str, Any]], recipes: Sequence[Recipe], workers: int,
                   th: Mapping[str, Any]) -> Dict[str, Any]:
    """R3's premise, checked: on EXTREMES_RIDS, EXTREMES_SAMPLES further passing draws per rival (seeded, not among the
    points); does any widen an output's spread beyond what B0 + the extremes give, and does any verdict change?"""
    import numpy as np

    import fit_m1b_ph_balance as FM
    from src.kinetic_core import parallel

    masks = passing_masks()
    rec = data_access.load_json(POINTS_OUT)
    bp = rec["davidek_breakpoints_in_situ"]
    rng = np.random.default_rng(EXTREMES_SEED)
    extra: Dict[str, Tuple[Dict[str, Any], Any]] = {}
    for rival in CARRIED:
        taken = {int(p["index"]) for p in rec["rivals"][rival]}
        idx = [int(i) for i in np.flatnonzero(masks[rival]["mask"]) if int(i) not in taken]
        for j in rng.choice(idx, size=min(EXTREMES_SAMPLES, len(idx)), replace=False):
            extra[f"{rival}_s{int(j)}"] = FM.rival_changes(rival, "primary", masks[rival]["draws"][int(j)], False, bp)
    by_rid = {r["rid"]: r for r in rows}
    items = [(r, extra) for r in recipes if r.rid in EXTREMES_RIDS]
    t0 = time.time()
    results = parallel.map_ordered(_extremes_task, items, n_workers=workers)
    per: List[Dict[str, Any]] = []
    widened, flips, n = 0, 0, 0
    worst = 0.0
    for res in results:
        row = by_rid[res["rid"]]
        for t in res["targets"] + ["ph_cooled_end"]:
            base_vals = _values(row, t)
            if base_vals is None:
                continue
            ex = []
            for v in res["extra"].values():
                if v.get("failed"):
                    continue
                ex.append(v.get(t) if t == "ph_cooled_end" else (v.get("concentrations_ug_per_l") or {}).get(t))
            ex = [float(x) for x in ex if x is not None]
            if t == "ph_cooled_end":
                sp0 = max(base_vals) - min(base_vals)
                sp1 = max(base_vals + ex) - min(base_vals + ex)
                tol = th["ph_level"]
            else:
                sp0, sp1, tol = log_spread(base_vals), log_spread(base_vals + ex), th["tau_dex"]
            n += 1
            if sp1 > sp0 + 1e-9:
                widened += 1
                worst = max(worst, (sp1 - sp0) if math.isfinite(sp1 - sp0) else math.inf)
            if (sp0 > tol) != (sp1 > tol):
                flips += 1
            per.append({"rid": res["rid"], "output": t, "spread_extremes": sp0, "spread_with_samples": sp1,
                        "tolerance": tol, "verdict_changes": (sp0 > tol) != (sp1 > tol)})
    summary = (f"{len(items)} recorded pots x {len(extra)} further passing draws ({EXTREMES_SAMPLES} per rival, seed "
               f"{EXTREMES_SEED}): {widened} of {n} outputs' spreads widen beyond B0 + the extremes (largest widening "
               f"{_fmt(worst)}); {flips} verdicts change. {time.time() - t0:.0f} s.")
    return {"summary": summary, "n_outputs": n, "n_widened": widened, "largest_widening": worst,
            "n_verdict_changes": flips, "samples_per_rival": EXTREMES_SAMPLES, "seed": EXTREMES_SEED,
            "pots": list(EXTREMES_RIDS), "per_output": per}


DRAWS_SIZES: Tuple[int, ...] = (5, 10, 20, 40, 80)
DRAWS_REPEATS = 200
#: an output's spread is "captured" when a subset's spread is within this of the all-passing-draws spread (dex or pH)
DRAWS_CAPTURE = 0.05


def draws_needed(rows: Sequence[Mapping[str, Any]], recipes: Sequence[Recipe], workers: int,
                 th: Mapping[str, Any]) -> Dict[str, Any]:
    """When the extremes do not bound the band: on EXTREMES_RIDS, EVERY passing draw of both rivals is run, and random
    subsets of n draws per rival (with B0) are compared with the all-draws spread: the share of (output, subset) within
    DRAWS_CAPTURE of it, and the share whose verdict agrees, per n; the extremes (+ best draw) beside them."""
    import numpy as np

    import fit_m1b_ph_balance as FM
    from src.kinetic_core import parallel

    masks = passing_masks()
    rec = data_access.load_json(POINTS_OUT)
    bp = rec["davidek_breakpoints_in_situ"]
    configs: Dict[str, Tuple[Dict[str, Any], Any]] = {}
    groups: Dict[str, List[str]] = {}
    for rival in CARRIED:
        for j in np.flatnonzero(masks[rival]["mask"]):
            name = f"{rival}_a{int(j)}"
            configs[name] = FM.rival_changes(rival, "primary", masks[rival]["draws"][int(j)], False, bp)
            groups.setdefault(rival, []).append(name)
    ext = {r: [f"{r}_a{int(p['index'])}" for p in rec["rivals"][r]] for r in CARRIED}
    by_rid = {r["rid"]: r for r in rows}
    items = [(r, configs) for r in recipes if r.rid in EXTREMES_RIDS]
    t0 = time.time()
    results = parallel.map_ordered(_extremes_task, items, n_workers=workers)
    rng = np.random.default_rng(EXTREMES_SEED + 1)
    outputs = []
    for res in results:
        row = by_rid[res["rid"]]
        b0 = (row["points"].get("B0") or {})
        for t in res["targets"] + ["ph_cooled_end"]:
            def val(v):
                if v.get("failed"):
                    return None
                return v.get(t) if t == "ph_cooled_end" else (v.get("concentrations_ug_per_l") or {}).get(t)
            b0v = b0.get(t) if t == "ph_cooled_end" else (b0.get("concentrations_ug_per_l") or {}).get(t)
            vals = {n: val(v) for n, v in res["extra"].items()}
            if b0v is None or any(x is None for x in vals.values()):
                continue
            outputs.append((res["rid"], t, float(b0v), {n: float(x) for n, x in vals.items()}))

    def spread(t, xs):
        return (max(xs) - min(xs)) if t == "ph_cooled_end" else log_spread(xs)

    def tol(t):
        return th["ph_level"] if t == "ph_cooled_end" else th["tau_dex"]

    full = {(rid, t): spread(t, [b0v] + list(v.values())) for rid, t, b0v, v in outputs}
    table = []
    for n in DRAWS_SIZES + ("extremes",):
        captured = agree = total = 0
        for rid, t, b0v, v in outputs:
            reps = 1 if n == "extremes" else DRAWS_REPEATS
            for _ in range(reps):
                if n == "extremes":
                    names = ext["B1"] + ext["B13"]
                else:
                    names = [x for r in CARRIED for x in rng.choice(groups[r], size=min(n, len(groups[r])), replace=False)]
                sp = spread(t, [b0v] + [v[x] for x in names if x in v])
                f = full[(rid, t)]
                total += 1
                captured += int((sp >= f - DRAWS_CAPTURE) or (not math.isfinite(f) and not math.isfinite(sp)))
                agree += int((sp > tol(t)) == (f > tol(t)))
        table.append({"draws_per_rival": n, "share_within_capture": captured / total, "share_verdict_agrees": agree / total,
                      "n": total})
    need = next((r["draws_per_rival"] for r in table if r["draws_per_rival"] != "extremes"
                 and r["share_within_capture"] >= 0.95), None)
    summary = (f"every passing draw (B1 {len(groups['B1'])}, B13 {len(groups['B13'])}) on {len(items)} recorded pots, "
               f"{len(outputs)} outputs: the extremes (+ best draw) reach within {DRAWS_CAPTURE} of the all-draws spread "
               f"in {table[-1]['share_within_capture']:.2f} of outputs (verdict agrees in "
               f"{table[-1]['share_verdict_agrees']:.2f}); random subsets reach it in 95 % of (output, subset) at "
               f"{need if need is not None else 'more than ' + str(DRAWS_SIZES[-1])} draws per rival. "
               f"{time.time() - t0:.0f} s.")
    return {"summary": summary, "table": table, "draws_needed_for_95pct": need, "capture": DRAWS_CAPTURE,
            "repeats": DRAWS_REPEATS, "pots": list(EXTREMES_RIDS),
            "per_output": [{"rid": rid, "output": t, "spread_all_draws": full[(rid, t)]} for rid, t, _, _ in outputs]}


def write_detail(rows: Sequence[Mapping[str, Any]], analysis: Mapping[str, Any]) -> None:
    """The recipe rows and every per-output / per-ratio entry, gzipped (deterministic: no timestamp in the header)."""
    import gzip
    import json

    body = json.dumps({"recipes": list(rows), "outputs": analysis["outputs"], "ratios": analysis["ratios"]},
                      separators=(",", ":")).encode("utf-8")
    with open(DETAIL, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(body)


def read_detail() -> Dict[str, Any]:
    import gzip
    import json

    with gzip.open(DETAIL, "rb") as fh:
        return json.loads(fh.read().decode("utf-8"))


def inputs() -> List[Path]:
    return [PREREG, FALSIFICATION, SCREEN, A2, S3_RECORD, POINTS_OUT,
            ROOT / "src/kinetic_core/ph_balance_rivals.py", ROOT / "src/kinetic_core/ph_state.py",
            ROOT / "src/kinetic_core/sulfur.py", ROOT / "src/kinetic_core/parameters_sulfur.py",
            data_paths.VALIDATION_DIR / "kinetic_core_b60_fit_report.json",
            ROOT / "docs/validation/directional_claims_panel.yml"]


def _memory_free_percent() -> Optional[float]:
    import subprocess

    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, timeout=30).stdout
    except Exception:  # noqa: BLE001 - not on macOS
        return None
    for line in out.splitlines()[::-1]:
        if "free percentage" in line:
            return float(line.rsplit(":", 1)[1].strip().rstrip("%"))
    return None


def run(recipes: Sequence[Recipe], workers: int, chunk: int = 64, min_free: float = 35.0,
        log: Callable[[str], None] = print, checkpoint: Optional[Path] = None) -> List[Dict[str, Any]]:
    """Every recipe, in chunks behind the memory gate; with ``checkpoint`` (a JSON-lines file) the rows already
    written there are reused and each chunk is appended, so an interrupted run resumes."""
    import json

    from src.kinetic_core import parallel

    done: Dict[str, Dict[str, Any]] = {}
    if checkpoint is not None and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["rid"]] = r
        log(f"checkpoint: {len(done)} recipes reused from {checkpoint}")
    todo = [r for r in recipes if r.rid not in done]
    rows: List[Dict[str, Any]] = []
    for start in range(0, len(todo), chunk):
        while True:
            free = _memory_free_percent()
            if free is None or free >= min_free:
                break
            log(f"memory free {free:.0f} % < {min_free:.0f} %: waiting 60 s")
            time.sleep(60)
        batch = list(todo[start:start + chunk])
        new = parallel.map_ordered(evaluate_recipe, batch, n_workers=workers)
        rows.extend(new)
        if checkpoint is not None:
            with checkpoint.open("a") as fh:
                for r in new:
                    fh.write(json.dumps(r) + "\n")
        log(f"{len(rows) + len(done)} / {len(recipes)} recipes")
    by_rid = {r["rid"]: r for r in rows}
    by_rid.update({k: v for k, v in done.items() if k not in by_rid})
    return [by_rid[r.rid] for r in recipes]


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--workers", type=int, default=int(os.environ.get("MAILLARD_WORKERS", "1")))
    ap.add_argument("--time-one", action="store_true", help="time one recipe at every point and project")
    ap.add_argument("--points", action="store_true", help="write the rivals' points record (R3) and stop")
    ap.add_argument("--limit", type=int, default=None, help="run only the first N recipes (a smoke run; no record)")
    ap.add_argument("--checkpoint", type=Path, default=None, help="JSON-lines file of finished recipes (resume)")
    ap.add_argument("--draws-needed", action="store_true",
                    help="run every passing draw on the extremes-check pots and add the convergence table to the record")
    ap.add_argument("--reanalyse", action="store_true",
                    help="re-read the record's own recipe rows and rewrite its analysis (no integration)")
    ap.add_argument("--projected-hours", type=float, default=None, help="the --time-one projection, for the record")
    args = ap.parse_args(argv)
    if args.points:
        payload = build_points()
        artifact_io.write_artifact(payload, POINTS_OUT)
        print(f"wrote {data_paths.rel(POINTS_OUT)}: " + "; ".join(f"{r}: {c['passing']} of {c['stored_draws']} draws pass, "
                                                                f"{c['points']} points" for r, c in payload["counts"].items()))
        return 0
    th = thresholds()
    grid = grid_recipes()
    pan = panel_recipes()
    dirr, claims = directional_recipes()
    recipes = pan + dirr + grid
    print(f"recipes: panel {len(pan)}, directional arms {len(dirr)} ({len(claims)} claims), grid {len(grid)}",
          flush=True)
    if args.draws_needed:
        record = data_access.load_json(OUT)
        rows = read_detail()["recipes"]
        record["draws_needed"] = draws_needed(rows, recipes, args.workers, th)
        artifact_io.write_artifact(record, OUT, render=render)
        print(record["draws_needed"]["summary"])
        return 0
    if args.reanalyse:
        record = data_access.load_json(OUT)
        rows = read_detail()["recipes"] if DETAIL.exists() else record.pop("recipes")
        analysis = analyse(rows, claims, th)
        write_detail(rows, analysis)
        record.pop("analysis", None)
        record.update(summary=analysis["summary"], directional=analysis["directional"], the_list=the_list(analysis),
                      thresholds=th, detail_file=data_paths.rel(DETAIL),
                      provenance=provenance.provenance_block("ph_balance_s3_sensitivity", generated_by=GENERATED_BY,
                                                             wave="S3 application: the pH balance (M1b part (i))",
                                                             inputs=inputs()))
        record["timing"]["reanalysed"] = True
        artifact_io.write_artifact(record, OUT, render=render)
        print(f"re-analysed {data_paths.rel(OUT)} from its {len(rows)} recipe rows")
        return 0
    if args.time_one:
        probe = next(r for r in grid if r.meta["pset"] == "cys50_rib35" and r.meta["buffer"] == "pyrophosphate_0.2M")
        row = evaluate_recipe(probe)
        secs = {k: v.get("seconds") for k, v in row["points"].items()}
        per_recipe = sum(v for v in secs.values() if v)
        print(f"one recipe: {secs}; {per_recipe:.1f} s per recipe; projected {len(recipes)} recipes at "
              f"{max(args.workers, 1)} workers: {len(recipes) * per_recipe / max(args.workers, 1) / 3600.0:.2f} h")
        return 0
    if args.limit is not None:
        recipes = recipes[: args.limit]
    t0 = time.time()
    rows = run(recipes, args.workers, log=lambda m: print(m, flush=True), checkpoint=args.checkpoint)
    wall = time.time() - t0
    analysis = analyse(rows, claims, th)
    secs: Dict[str, List[float]] = defaultdict(list)
    failed = 0
    for r in rows:
        for p, v in (r.get("points") or {}).items():
            if v.get("failed"):
                failed += 1
            if v.get("seconds") is not None and not v.get("not_applicable"):
                secs[p].append(float(v["seconds"]))
    payload = {
        "artifact": "ph_balance_s3_sensitivity",
        "provenance": provenance.provenance_block("ph_balance_s3_sensitivity", generated_by=GENERATED_BY,
                                                  wave="S3 application: the pH balance (M1b part (i))",
                                                  inputs=inputs()),
        "definition": __doc__.split("    python")[0].strip(),
        "thresholds": th,
        "points": point_coordinates(),
        "timing": {"recipes": len(rows), "workers": args.workers, "wall_seconds": wall,
                   "projected_hours": args.projected_hours,
                   "median_seconds": {p: statistics.median(v) for p, v in secs.items() if v}},
        "failed_point_runs": failed,
        "complete": args.limit is None,
        "detail_file": data_paths.rel(DETAIL),
        "summary": analysis["summary"],
        "directional": analysis["directional"],
    }
    rival_names = set(point_configs()) - {"B0", "B0_inline"}
    per_recipe = [sum(float(v["seconds"]) for p, v in (r.get("points") or {}).items()
                      if p in rival_names and v.get("seconds") is not None and not v.get("not_applicable"))
                  for r in rows if not r.get("clamped") and r.get("targets")]
    payload["timing"]["rivals_seconds_per_pot"] = ({"median": statistics.median(per_recipe), "max": max(per_recipe),
                                                    "n_points": len(rival_names)} if per_recipe else None)
    payload["the_list"] = the_list(analysis)
    if args.limit is None:
        payload["extremes_check"] = extremes_check(rows, recipes, args.workers, th)
    if args.limit is not None:
        print("smoke run (--limit): no record written")
        print(render(payload)[:6000])
        return 0
    write_detail(rows, analysis)
    artifact_io.write_artifact(payload, OUT, render=render)
    print(f"wrote {data_paths.rel(OUT)} ({wall / 3600.0:.2f} h, {failed} failed point runs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
