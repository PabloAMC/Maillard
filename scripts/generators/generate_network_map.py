#!/usr/bin/env python
"""
scripts/generators/generate_network_map.py -- THE REACTION NETWORK, FROM THE
LIVE CODE (Build Wave V1).

WHAT IT MAKES
-------------
``docs/assets/network_map.html``: one self-contained page showing all FOUR
lanes of the kinetic core -- trunk, sulfur, acrylamide, lipid -- with

  * one node per SPECIES, from the live species tables;
  * one edge per REACTION, from the live reaction tuples;
  * edges COLOURED BY EVIDENCE CLASS, read off each step's parameter in the
    frozen registries: measured / fitted / derived / pinned;
  * a native tooltip on every node and edge carrying the SOURCE ANCHOR and the
    VALIDITY WINDOW (temperature range, pH of measurement, transfer licence).

Nothing here is a drawing of what the network is *supposed* to be. Every node,
every edge and every colour is read at generation time from
``src/kinetic_core/``, so a step added, retired or re-classified moves this
picture on the next run. That is the point: a hand-drawn network diagram is a
claim that decays silently, and this repository has already been bitten once by
prose that stopped matching its own artifacts.

FLUX MODE -- WHICH CHEMISTRY ACTUALLY DOMINATES *YOUR* PROCESS
--------------------------------------------------------------
Given a formulation + process spec, the page additionally sets EDGE WIDTH by
the TIME-INTEGRATED ABSOLUTE FLUX through each step (mmol/L over the whole
thermal program), computed by integrating the real network and trapezoid-
integrating the mass-action rates. A hairline is a step that carried almost no
material in your process; a thick edge is where your chemistry actually went.

The same network looks completely different at 145 C for 20 min in water and in
a three-zone extruder, and that difference is the most useful thing this
picture can show a formulator.

DETERMINISM
-----------
Byte-identical on repeated runs: no timestamp, no wall-clock, sorted iteration
throughout, fixed float formatting, and the short commit hash (not
``--dirty``) as the version stamp. ``tests/unit/test_v1_reports.py`` runs it
twice and compares bytes.

USAGE
-----
    python scripts/generators/generate_network_map.py
    python scripts/generators/generate_network_map.py --spec docs/examples/flux_ribose_cysteine_145C.json
    python scripts/generators/generate_network_map.py --all-examples
    python scripts/generators/generate_network_map.py --out /tmp/map.html
"""

from __future__ import annotations

import argparse
import html as _html
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402
from src.explain_compound import (  # noqa: E402
    DERIVED,
    EVIDENCE_MEANING,
    EVIDENCE_ORDER,
    FITTED,
    MEASURED,
    PINNED,
    evidence_class_of,
)

DEFAULT_OUT = data_paths.DOCS_ASSETS_DIR / "network_map.html"
EXAMPLES_DIR = data_paths.EXAMPLES_DIR
DEFAULT_EXAMPLES = (
    "flux_ribose_cysteine_145C.json",
    "flux_extrusion_pea_three_zone.json",
    "flux_acrylamide_fry_180C.json",
)

EVIDENCE_COLOUR: Mapping[str, str] = {
    MEASURED: "#2e6b4f",
    FITTED: "#2f5f8f",
    DERIVED: "#8a6d1f",
    PINNED: "#8f2b1e",
}

#: Species that are bookkeeping pools rather than chemistry. Drawn, but muted:
#: hiding them would make the carbon balance look like it does not close, and
#: this repository routes every unmeasured co-product into them ON PURPOSE.
POOL_ROLES = ("pool",)


def esc(value: Any) -> str:
    return _html.escape("" if value is None else str(value), quote=True)


def num(value: Optional[float], places: int = 4) -> str:
    """Fixed-precision formatting. Determinism depends on this being total."""
    if value is None:
        return "-"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != number or number in (float("inf"), float("-inf")):
        return str(number)
    if number == 0.0:
        return "0"
    if abs(number) >= 1e4 or abs(number) < 1e-3:
        return f"{number:.{places}e}"
    return f"{number:.{places}g}"


def _commit() -> str:
    """Short HEAD hash. Deliberately NOT --dirty: this file must be reproducible."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(ROOT), capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError):  # pragma: no cover - defensive
        return "unknown"
    return out.stdout.strip() if out.returncode == 0 else "unknown"


# ===========================================================================
# 1. READING THE LIVE NETWORK
# ===========================================================================


def lane_species(lane: str) -> Tuple[Any, ...]:
    if lane == "trunk":
        from src.kinetic_core.species import SPECIES

        return tuple(SPECIES)
    if lane == "sulfur":
        from src.kinetic_core.species_sulfur import SULFUR_STATE

        return tuple(SULFUR_STATE)
    if lane == "acrylamide":
        from src.kinetic_core.species_acrylamide import ACRYLAMIDE_STATE

        return tuple(ACRYLAMIDE_STATE)
    if lane == "lipid":
        from src.kinetic_core.species_lipid import LIPID_SPECIES

        return tuple(LIPID_SPECIES)
    raise ValueError(f"unknown lane {lane!r}")


def lane_reactions(lane: str) -> Tuple[Any, ...]:
    if lane == "trunk":
        from src.kinetic_core.network import REACTIONS

        return tuple(REACTIONS)
    if lane == "sulfur":
        from src.kinetic_core.sulfur import FULL_REACTIONS

        return tuple(FULL_REACTIONS)
    if lane == "acrylamide":
        from src.kinetic_core.acrylamide import FULL_ACRYLAMIDE_REACTIONS

        return tuple(FULL_ACRYLAMIDE_REACTIONS)
    return ()


def _lipid_pseudo_reactions() -> List[Dict[str, Any]]:
    """
    The lipid lane has no mass-action tuple; it has a BRANCH MODEL.

    Rendered as one pseudo-reaction per (pool, product) with a non-zero branch
    share, because that is what the lane actually does: a first-order pool
    decomposition distributed over Frankel's measured slate. The RATE edge is
    drawn separately and is the one edge on this whole page whose class is
    PINNED because the constant is an assumption rather than a measurement.
    """
    from src.kinetic_core.engine import core_lipid_model
    from src.kinetic_core.parameters_lipid import K_LOOH_DECOMP_ANCHOR, Q10_ASSUMPTION
    from src.kinetic_core.species_lipid import CLEAVAGE_MECHANISM, LOOH_POOLS

    branch, _composition = core_lipid_model()
    rows: List[Dict[str, Any]] = []
    for pool in sorted(LOOH_POOLS):
        position, geometry = LOOH_POOLS[pool]
        shares = branch.simplexes[(position, geometry)]
        for product in sorted(shares):
            share = float(shares[product])
            if share <= 0.0:
                continue
            rows.append(
                {
                    "key": f"lipid_{pool}_{product}",
                    "source": pool,
                    "target": product,
                    "evidence": FITTED,
                    "raw_class": "derived_from_fit_data",
                    "parameter_key": "branch_model",
                    "anchor": branch.provenance,
                    "window": (
                        "branch DISTRIBUTION is fitted to Frankel 1989's "
                        "zero-additive columns; structural zeros enforced"
                    ),
                    "detail": (
                        f"branch share of this pool {share:.4f}; cleavage "
                        f"mechanism {CLEAVAGE_MECHANISM.get(product, 'unassigned')}"
                    ),
                    "share": share,
                }
            )
    for pool in sorted(LOOH_POOLS):
        rows.append(
            {
                "key": f"lipid_rate_{pool}",
                "source": pool,
                "target": "LIPID_FRAG_C",
                "evidence": PINNED,
                "raw_class": "declared_assumption",
                "parameter_key": "k_LOOH_decomp",
                "anchor": K_LOOH_DECOMP_ANCHOR.source_anchor,
                "window": (
                    f"anchored at "
                    f"{K_LOOH_DECOMP_ANCHOR.temperature_of_measurement_c:g} C, "
                    f"pH {K_LOOH_DECOMP_ANCHOR.ph_of_measurement:g}; Q10 "
                    f"{Q10_ASSUMPTION.lo:g}-{Q10_ASSUMPTION.hi:g} licensed only "
                    f"over {Q10_ASSUMPTION.licensed_span_c[0]:g}-"
                    f"{Q10_ASSUMPTION.licensed_span_c[1]:g} C"
                ),
                "detail": (
                    "THE RATE IS AN ASSUMPTION, NOT A MEASUREMENT. Carbon not "
                    "named by the measured slate is routed here so the balance "
                    "closes as an equality."
                ),
                "share": None,
            }
        )
    return rows


def edge_records(lane: str) -> List[Dict[str, Any]]:
    """Every drawable edge of one lane, with its evidence and provenance."""
    if lane == "lipid":
        return _lipid_pseudo_reactions()

    from src.explain_compound import _parameter_record

    parameters = {}
    try:
        from src.kinetic_core.engine import core_parameters

        parameters = dict(core_parameters(lane))
    except Exception as exc:  # pragma: no cover - a missing fit report
        raise SystemExit(
            f"cannot read the frozen parameters for the {lane} lane: {exc}"
        ) from exc

    rows: List[Dict[str, Any]] = []
    for reaction in lane_reactions(lane):
        record = _parameter_record(reaction, parameters, lane)
        window_bits: List[str] = []
        if record.get("temperature_range_c"):
            low, high = record["temperature_range_c"]
            window_bits.append(f"valid {low:g}-{high:g} C")
        if record.get("ph_of_measurement") is not None:
            window_bits.append(f"measured at pH {record['ph_of_measurement']:g}")
        if record.get("ea_kj_mol") is not None:
            window_bits.append(f"Ea {float(record['ea_kj_mol']):g} kJ/mol")
        else:
            window_bits.append("NO activation energy: the constant is held fixed")
        if record.get("rate_transfer"):
            window_bits.append(f"transfer {record['rate_transfer']}")
        equation = (
            " + ".join(
                f"{c} {k}" if c != 1 else k
                for k, c in sorted(reaction.reactants.items())
            )
            + " -> "
            + " + ".join(
                f"{c} {k}" if c != 1 else k
                for k, c in sorted(reaction.products.items())
            )
        )
        for source in sorted(reaction.reactants):
            for target in sorted(reaction.products):
                rows.append(
                    {
                        "key": reaction.key,
                        "source": source,
                        "target": target,
                        "evidence": record["evidence"],
                        "raw_class": record.get("raw_evidence_class"),
                        "parameter_key": record.get("parameter_key"),
                        "anchor": record.get("source_anchor", ""),
                        "window": "; ".join(window_bits),
                        "detail": equation
                        + (f" | {record['note']}" if record.get("note") else "")
                        + (
                            " | flags: " + ", ".join(record["flags"])
                            if record.get("flags")
                            else ""
                        ),
                        "share": None,
                    }
                )
    return rows


# ===========================================================================
# 2. FLUX MODE
# ===========================================================================


def load_flux_spec(path: Path) -> Dict[str, Any]:
    """
    Read a flux spec. Accepts the CLI's own spec plus one extension.

    EXTENSION: ``thermal_segments: [[duration_min, temperature_C], ...]``, so
    that an extrusion profile can be written as the multi-zone program it
    physically is rather than collapsed to a single isothermal hold. The engine
    has always supported a piecewise-constant program; nothing but the front
    door was missing.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{path}: a flux spec must be a JSON object")
    if "thermal_segments" not in data:
        if "temp_C" not in data or "time_min" not in data:
            raise SystemExit(
                f"{path}: needs either 'thermal_segments' or both 'temp_C' and "
                f"'time_min'"
            )
        data["thermal_segments"] = [[float(data["time_min"]), float(data["temp_C"])]]
    data["thermal_segments"] = [
        [float(duration), float(temperature)]
        for duration, temperature in data["thermal_segments"]
    ]
    return data


def _spec_objects(spec: Mapping[str, Any]):
    from src.kinetic_core.engine import FormulationSpec, ProcessSpec, ThermalProgram

    process = ProcessSpec(
        thermal=ThermalProgram(
            tuple((float(d), float(t)) for d, t in spec["thermal_segments"])
        ),
        ph=float(spec.get("ph", 6.8)),
        water_activity=(
            float(spec["aw"]) if spec.get("aw") is not None else None
        ),
        matrix=str(spec.get("matrix") or spec.get("protein_type") or "water"),
    )
    return FormulationSpec(
        name=str(spec.get("name") or Path(str(spec.get("_source", "spec"))).stem),
        precursors={str(k): float(v) for k, v in dict(spec["precursors"]).items()},
        process=process,
    )


def flux_budget(spec: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Time-integrated absolute flux through every step, mmol/L, for one spec.

    Segment-chained: each thermal segment starts from the previous segment's
    final state, exactly as ``engine._integrate_program`` does, so a multi-zone
    program is integrated as a program and not as a set of independent holds.
    The rates are re-evaluated on a fine grid inside each segment and trapezoid-
    integrated; nothing is approximated by an end-point rate.
    """
    import numpy as np

    from src.kinetic_core.engine import (
        ACRYLAMIDE,
        CELSIUS,
        LIPID,
        SULFUR,
        TRUNK,
        core_parameters,
        declare_envelope,
        default_targets_for,
    )

    core_spec = _spec_objects(spec)
    targets = list(spec.get("targets") or default_targets_for(core_spec.precursors))
    declaration = declare_envelope(core_spec, targets)
    result: Dict[str, Any] = {
        "name": core_spec.name,
        "state": declaration.state,
        "lanes": list(declaration.lanes),
        "reasons": list(declaration.reasons),
        "warnings": list(declaration.warnings),
        "thermal": core_spec.process.thermal.describe(),
        "ph": core_spec.process.ph,
        "matrix": core_spec.process.matrix,
        "targets": targets,
        "flux": {},
        "n_points": 121,
    }
    if not declaration.is_answerable:
        return result

    lanes = list(declaration.lanes)
    maillard = next((lane for lane in lanes if lane in (TRUNK, SULFUR, ACRYLAMIDE)), None)
    flux: Dict[str, float] = {}
    grid_points = result["n_points"]

    if maillard is not None:
        parameters = core_parameters(maillard)
        state = dict(declaration.mapped_precursors)
        if maillard == SULFUR:
            from src.kinetic_core.ph_state import DEFAULT_BUFFER
            from src.kinetic_core.engine import core_ph_drift
            from src.kinetic_core.species_sulfur import SULFUR_INDEX
            from src.kinetic_core.sulfur import (
                FULL_REACTIONS,
                integrate_sulfur,
                sulfur_rate_constants_at,
                sulfur_reaction_rates,
            )

            keys = list(SULFUR_INDEX)
            reactions = FULL_REACTIONS
            for duration, temperature_c in core_spec.process.thermal.segments:
                grid = np.linspace(0.0, float(duration), grid_points)
                run = integrate_sulfur(
                    parameters, float(temperature_c) + CELSIUS, state, grid,
                    ph=float(core_spec.process.ph),
                    buffer_spec=DEFAULT_BUFFER,
                    ph_drift=core_ph_drift(),
                    rtol=1e-8, atol=1e-14,
                )
                rows = []
                for index in range(len(grid)):
                    if run.ph_series is not None:
                        ph_here = float(run.ph_series[index])
                    else:  # pragma: no cover - clamped buffers only
                        ph_here = float(core_spec.process.ph)
                    k_t = sulfur_rate_constants_at(
                        parameters, float(temperature_c) + CELSIUS, ph_here
                    )
                    rows.append(sulfur_reaction_rates(run.concentrations[index], k_t))
                integral = np.trapezoid(np.array(rows), grid, axis=0)
                for position, reaction in enumerate(reactions):
                    flux[reaction.key] = flux.get(reaction.key, 0.0) + float(
                        integral[position]
                    )
                state = {
                    key: float(run.concentrations[-1, i]) for i, key in enumerate(keys)
                }
        elif maillard == ACRYLAMIDE:
            from src.kinetic_core.acrylamide import (
                FULL_ACRYLAMIDE_REACTIONS,
                acrylamide_rate_constants_at,
                acrylamide_rate_vector,
                acrylamide_reaction_rates,
                integrate_acrylamide,
            )
            from src.kinetic_core.species_acrylamide import ACRYLAMIDE_INDEX

            keys = list(ACRYLAMIDE_INDEX)
            for duration, temperature_c in core_spec.process.thermal.segments:
                grid = np.linspace(0.0, float(duration), grid_points)
                run = integrate_acrylamide(
                    parameters, float(temperature_c) + CELSIUS, state, grid,
                    water_activity=core_spec.process.water_activity,
                    rtol=1e-8, atol=1e-14,
                )
                k = acrylamide_rate_vector(
                    acrylamide_rate_constants_at(
                        parameters, float(temperature_c) + CELSIUS
                    )
                )
                rates = np.array(
                    [acrylamide_reaction_rates(row, k) for row in run.concentrations]
                )
                integral = np.trapezoid(rates, grid, axis=0)
                for position, reaction in enumerate(FULL_ACRYLAMIDE_REACTIONS):
                    flux[reaction.key] = flux.get(reaction.key, 0.0) + float(
                        integral[position]
                    )
                state = {
                    key: float(run.concentrations[-1, i]) for i, key in enumerate(keys)
                }
        else:
            from src.kinetic_core.network import (
                REACTIONS,
                rate_constants_at,
                reaction_rates,
            )
            from src.kinetic_core.species import SPECIES_KEYS
            from src.kinetic_core.integrate import integrate

            keys = list(SPECIES_KEYS)
            for duration, temperature_c in core_spec.process.thermal.segments:
                grid = np.linspace(0.0, float(duration), grid_points)
                run = integrate(
                    parameters, float(temperature_c) + CELSIUS, state, grid,
                    rtol=1e-8, atol=1e-14,
                )
                k = rate_constants_at(parameters, float(temperature_c) + CELSIUS)
                rates = np.array(
                    [reaction_rates(row, k) for row in run.concentrations]
                )
                integral = np.trapezoid(rates, grid, axis=0)
                for position, reaction in enumerate(REACTIONS):
                    flux[reaction.key] = flux.get(reaction.key, 0.0) + float(
                        integral[position]
                    )
                state = {
                    key: float(run.concentrations[-1, i]) for i, key in enumerate(keys)
                }

    if LIPID in lanes:
        flux.update(_lipid_flux(core_spec, declaration))

    result["flux"] = {key: float(value) for key, value in sorted(flux.items())}
    return result


def _lipid_flux(core_spec, declaration) -> Dict[str, float]:
    """
    Flux through the lipid lane's branch edges, mmol/L.

    Each pool's DECOMPOSED amount over the whole program times that pool's
    branch share for the product, which is exactly the material that travelled
    down that edge. The RATE edges carry the pool's total decomposition.
    """
    from src.kinetic_core.engine import core_lipid_model
    from src.kinetic_core.lipid import charge_from_carrier, integrate_lipid
    from src.kinetic_core.parameters_lipid import LIPID_CARRIERS
    from src.kinetic_core.species_lipid import LOOH_POOLS

    branch, composition = core_lipid_model()
    carriers = [c for c in declaration.lipid_carriers if c in LIPID_CARRIERS]
    if not carriers:
        return {}
    segments = list(core_spec.process.thermal.segments)

    decomposed: Dict[str, float] = {pool: 0.0 for pool in LOOH_POOLS}
    for carrier_key in sorted(carriers):
        charge = charge_from_carrier(LIPID_CARRIERS[carrier_key], composition)
        run = integrate_lipid(charge, segments, branch)
        start = {
            "LOOH_13_ct": charge.looh_linoleate_mmol_l * composition.f13_ct,
            "LOOH_13_tt": charge.looh_linoleate_mmol_l * composition.f13_tt,
            "LOOH_9_ct": charge.looh_linoleate_mmol_l * composition.f9_ct,
            "LOOH_9_tt": charge.looh_linoleate_mmol_l * composition.f9_tt,
        }
        for pool in LOOH_POOLS:
            decomposed[pool] += start[pool] - float(run.state_mmol_per_l.get(pool, 0.0))

    out: Dict[str, float] = {}
    for pool in sorted(LOOH_POOLS):
        position, geometry = LOOH_POOLS[pool]
        out[f"lipid_rate_{pool}"] = decomposed[pool]
        for product, share in sorted(branch.simplexes[(position, geometry)].items()):
            out[f"lipid_{pool}_{product}"] = decomposed[pool] * float(share)
    return out


# ===========================================================================
# 3. LAYOUT -- deterministic layered placement, no libraries
# ===========================================================================


#: Nodes per column before a column is split into two. Keeps the sulfur lane's
#: fragment/product bulge from producing a 40-row column nothing can read.
MAX_PER_COLUMN = 12


def layer_assignment(
    node_keys: Sequence[str], edges: Sequence[Mapping[str, Any]]
) -> Dict[str, int]:
    """
    SHORTEST-PATH layering from the network's sources. Cycle-safe by design.

    Longest-path layering is the textbook choice and is wrong here: these
    networks contain real cycles (glucose <-> fructose, the thioether
    conjugation and its measured reverse), and relaxing a cycle unrolls it into
    as many columns as there are relaxation passes -- which produced a
    12 000-pixel-wide sulfur lane on the first attempt. A shortest-path layering
    assigns each node its distance from the nearest charged species, visits
    every node exactly once, and therefore cannot be inflated by a cycle.

    Deterministic: the frontier is sorted at every level.
    """
    keys = list(node_keys)
    in_edges: Dict[str, List[str]] = {key: [] for key in keys}
    out_edges: Dict[str, List[str]] = {key: [] for key in keys}
    for edge in edges:
        source, target = str(edge["source"]), str(edge["target"])
        if source not in in_edges or target not in in_edges or source == target:
            continue
        out_edges[source].append(target)
        in_edges[target].append(source)

    sources = sorted(key for key in keys if not in_edges[key])
    if not sources:  # every node is in a cycle; seed on out-degree instead
        sources = sorted(keys, key=lambda k: (-len(out_edges[k]), k))[:1]

    depth: Dict[str, int] = {key: 0 for key in sources}
    frontier = list(sources)
    while frontier:
        nxt: List[str] = []
        for key in sorted(frontier):
            for target in sorted(set(out_edges[key])):
                if target not in depth:
                    depth[target] = depth[key] + 1
                    nxt.append(target)
        frontier = nxt

    # Unreached nodes (only reachable against the arrow) go one column past the
    # deepest reached node rather than being silently dropped.
    tail = (max(depth.values()) + 1) if depth else 0
    for key in sorted(keys):
        depth.setdefault(key, tail)
    return depth


def layout(
    species: Sequence[Any], edges: Sequence[Mapping[str, Any]]
) -> Tuple[Dict[str, Tuple[float, float]], int, int]:
    keys = [s.key for s in species]
    depth = layer_assignment(keys, edges)
    grouped: Dict[int, List[str]] = {}
    for key in sorted(keys):
        grouped.setdefault(depth[key], []).append(key)

    # Split any over-full column into consecutive sub-columns, preserving order.
    columns: List[List[str]] = []
    for level in sorted(grouped):
        members = grouped[level]
        for start in range(0, len(members), MAX_PER_COLUMN):
            columns.append(members[start:start + MAX_PER_COLUMN])

    col_w, row_h = 168.0, 34.0
    left, top = 96.0, 44.0
    tallest = max((len(c) for c in columns), default=1)
    positions: Dict[str, Tuple[float, float]] = {}
    for index, members in enumerate(columns):
        block = len(members) * row_h
        offset = top + (tallest * row_h - block) / 2.0
        for row, key in enumerate(members):
            positions[key] = (left + index * col_w, offset + row * row_h)
    width = int(left + max(len(columns) - 1, 0) * col_w + 150)
    height = int(top + tallest * row_h + 60)
    return positions, width, height


# ===========================================================================
# 4. RENDERING
# ===========================================================================


NODE_W, NODE_H = 96.0, 21.0


def _node_title(species: Any, lane: str) -> str:
    bits = [
        f"{species.key} -- {species.label}",
        f"role: {species.role}",
        f"C {species.carbon}, N {species.nitrogen}, S {getattr(species, 'sulfur', 0)}",
        "MEASURED in the fit corpus" if species.measured
        else "not measured in the fit corpus",
    ]
    if getattr(species, "note", ""):
        bits.append(species.note)
    return "\n".join(bits)


def _edge_title(edge: Mapping[str, Any], flux: Optional[float], total: Optional[float]) -> str:
    bits = [
        f"{edge['key']}  [{edge['evidence'].upper()}]",
        edge.get("detail", ""),
        f"parameter: {edge.get('parameter_key') or '(derived, no parameter of its own)'}",
        f"evidence class in the registry: {edge.get('raw_class') or 'n/a'}",
        f"validity: {edge.get('window', '')}",
        f"source: {edge.get('anchor', '')}",
    ]
    if flux is not None:
        share = (100.0 * flux / total) if total else 0.0
        bits.append(
            f"FLUX in this process: {num(flux)} mmol/L ({share:.2f} % of the "
            f"lane's total step flux)"
        )
    return "\n".join(b for b in bits if b)


def _svg_lane(
    lane: str,
    edges: Sequence[Mapping[str, Any]],
    flux: Optional[Mapping[str, float]],
) -> str:
    species = lane_species(lane)
    by_key = {s.key: s for s in species}
    drawable = [e for e in edges if e["source"] in by_key and e["target"] in by_key]
    positions, width, height = layout(species, drawable)

    fluxes = {}
    if flux:
        for edge in drawable:
            value = abs(float(flux.get(edge["key"], 0.0)))
            fluxes[edge["key"]] = value
    total = sum(fluxes.values()) if fluxes else 0.0
    peak = max(fluxes.values()) if fluxes else 0.0

    def stroke_width(edge: Mapping[str, Any]) -> float:
        if not fluxes or peak <= 0.0:
            return 1.1
        value = fluxes.get(edge["key"], 0.0)
        if value <= 0.0:
            return 0.35
        # log scaling: five decades of flux compressed into 0.5-7 px, so a
        # dominant channel is unmistakable without erasing a minor one.
        decades = math.log10(peak / value) if value > 0 else 6.0
        return max(0.5, 7.0 - 1.3 * min(decades, 5.0))

    out: List[str] = [
        f'<svg class="net" viewBox="0 0 {width} {height}" width="100%" '
        f'role="img" aria-label="{esc(lane)} lane reaction network">'
    ]

    for edge in sorted(drawable, key=lambda e: (e["key"], e["source"], e["target"])):
        x0, y0 = positions[edge["source"]]
        x1, y1 = positions[edge["target"]]
        sx, sy = x0 + NODE_W / 2.0, y0
        tx, ty = x1 - NODE_W / 2.0, y1
        curve = max(abs(tx - sx) * 0.42, 26.0)
        if tx < sx:  # a back edge; bow it so it does not vanish under the nodes
            sy -= NODE_H / 2.0 + 3.0
            ty -= NODE_H / 2.0 + 3.0
        colour = EVIDENCE_COLOUR[edge["evidence"]]
        value = fluxes.get(edge["key"]) if fluxes else None
        out.append(
            f'<g class="e ev-{edge["evidence"]}" data-ev="{edge["evidence"]}">'
            f'<title>{esc(_edge_title(edge, value, total))}</title>'
            f'<path d="M {sx:.1f} {sy:.1f} C {sx + curve:.1f} {sy:.1f}, '
            f'{tx - curve:.1f} {ty:.1f}, {tx:.1f} {ty:.1f}" '
            f'fill="none" stroke="{colour}" stroke-width="{stroke_width(edge):.2f}" '
            f'stroke-opacity="0.62"/></g>'
        )

    for key in sorted(positions):
        species_row = by_key[key]
        x, y = positions[key]
        pool = species_row.role in POOL_ROLES
        klass = "n pool" if pool else ("n measured" if species_row.measured else "n")
        out.append(
            f'<g class="{klass}"><title>{esc(_node_title(species_row, lane))}</title>'
            f'<rect x="{x - NODE_W / 2:.1f}" y="{y - NODE_H / 2:.1f}" '
            f'width="{NODE_W}" height="{NODE_H}" rx="4"/>'
            f'<text x="{x:.1f}" y="{y + 4:.1f}" text-anchor="middle">{esc(key)}</text>'
            "</g>"
        )
    out.append("</svg>")
    return "".join(out)


CSS = """
:root{--ink:#1b1a17;--muted:#5f5b53;--rule:#ddd8cf;--paper:#fbfaf7;--card:#fff;--accent:#7a3e12}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
 font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.wrap{max-width:1180px;margin:0 auto;padding:26px 20px 60px}
h1{font-size:25px;margin:0 0 4px}
h2{font-size:18px;margin:32px 0 8px;padding-bottom:6px;border-bottom:2px solid var(--rule)}
h3{font-size:13px;margin:16px 0 6px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
p{margin:8px 0}
.sub{color:var(--muted);margin:0 0 14px}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.9em}
.card{background:var(--card);border:1px solid var(--rule);border-radius:7px;padding:12px 15px;margin:12px 0}
.rule{color:var(--muted);font-size:13.5px;border-left:3px solid var(--accent);padding-left:10px;margin:10px 0}
.note{color:var(--muted);font-size:13px}
.scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:7px;background:#fff;padding:6px}
svg.net{display:block;min-width:640px}
svg.net .n rect{fill:#f3f1ec;stroke:#b9b3a8;stroke-width:1}
svg.net .n.measured rect{fill:#e7efe9;stroke:#7fa792}
svg.net .n.pool rect{fill:#faf6ec;stroke:#d9cba6;stroke-dasharray:3 2}
svg.net .n text{font:11.5px ui-monospace,Menlo,monospace;fill:#26241f}
svg.net .e.dim path{stroke-opacity:.06}
.legend{display:flex;flex-wrap:wrap;gap:8px;margin:10px 0}
.legend label{display:inline-flex;align-items:center;gap:6px;font-size:13px;
 border:1px solid var(--rule);border-radius:16px;padding:3px 11px;background:#fff;cursor:pointer}
.sw{width:22px;height:3px;border-radius:2px;display:inline-block}
table.data{width:100%;border-collapse:collapse;font-size:13.5px}
table.data th{text-align:left;font-size:11.5px;text-transform:uppercase;letter-spacing:.05em;
 color:var(--muted);border-bottom:2px solid var(--rule);padding:5px 8px 5px 0}
table.data td{padding:5px 8px 5px 0;border-bottom:1px solid var(--rule);vertical-align:top}
table.data td.num,table.data th.num{text-align:right;font-family:ui-monospace,Menlo,monospace}
.bar{display:inline-block;height:9px;background:#8fa8bd;border-radius:2px;vertical-align:0}
footer{margin-top:40px;border-top:2px solid var(--rule);padding-top:14px;color:var(--muted);font-size:13px}
@media print{.legend{display:none}.scroll{overflow:visible}body{background:#fff}
 .card,svg.net,table.data{break-inside:avoid}}
"""

JS = """
(function(){
  var boxes=document.querySelectorAll('input[data-ev]');
  function apply(){
    var on={};
    boxes.forEach(function(b){on[b.getAttribute('data-ev')]=b.checked;});
    document.querySelectorAll('g.e').forEach(function(g){
      var ev=g.getAttribute('data-ev');
      if(on[ev]===false){g.classList.add('dim');}else{g.classList.remove('dim');}
    });
  }
  boxes.forEach(function(b){b.addEventListener('change',apply);});
  apply();
})();
"""


def _legend() -> str:
    items = "".join(
        f'<label><input type="checkbox" data-ev="{cls}" checked>'
        f'<span class="sw" style="background:{EVIDENCE_COLOUR[cls]}"></span>'
        f"<strong>{cls}</strong> &mdash; {esc(EVIDENCE_MEANING[cls])}</label>"
        for cls in EVIDENCE_ORDER
    )
    return f'<div class="legend">{items}</div>'


def _evidence_census(edges: Sequence[Mapping[str, Any]]) -> Dict[str, int]:
    steps: Dict[str, str] = {}
    for edge in edges:
        steps[edge["key"]] = edge["evidence"]
    census = {cls: 0 for cls in EVIDENCE_ORDER}
    for evidence in steps.values():
        census[evidence] += 1
    return census


def _flux_table(
    edges: Sequence[Mapping[str, Any]], flux: Mapping[str, float], top_n: int = 12
) -> str:
    by_step: Dict[str, Dict[str, Any]] = {}
    for edge in edges:
        by_step.setdefault(
            edge["key"],
            {"key": edge["key"], "evidence": edge["evidence"],
             "detail": edge.get("detail", ""), "anchor": edge.get("anchor", "")},
        )
    rows = []
    for key, entry in by_step.items():
        rows.append({**entry, "flux": abs(float(flux.get(key, 0.0)))})
    rows.sort(key=lambda r: (-r["flux"], r["key"]))
    total = sum(r["flux"] for r in rows)
    peak = rows[0]["flux"] if rows else 0.0
    body = "".join(
        '<tr><td><code>{key}</code></td><td>{ev}</td>'
        '<td class="num">{flux}</td><td class="num">{share}</td>'
        '<td><span class="bar" style="width:{px}px;background:{colour}"></span></td></tr>'.format(
            key=esc(row["key"]),
            ev=esc(row["evidence"]),
            flux=num(row["flux"]),
            share=f"{(100.0 * row['flux'] / total):.2f} %" if total else "-",
            px=int(round(150.0 * (row["flux"] / peak))) if peak else 0,
            colour=EVIDENCE_COLOUR[row["evidence"]],
        )
        for row in rows[:top_n]
    )
    return (
        '<table class="data"><thead><tr><th>step</th><th>evidence</th>'
        '<th class="num">flux, mmol/L</th><th class="num">share</th><th></th>'
        f"</tr></thead><tbody>{body}</tbody></table>"
    )


def build_page(
    flux_runs: Sequence[Mapping[str, Any]], *, commit: str
) -> str:
    from src.kinetic_core.engine import LANES, engine_metadata

    metadata = engine_metadata()
    parts: List[str] = [
        "<h1>Maillard kinetic core &mdash; reaction network map</h1>",
        '<p class="sub">Every node, every edge and every colour on this page is '
        "read from <code>src/kinetic_core/</code> at generation time. Nothing "
        "here is drawn by hand, so nothing here can silently stop being true.</p>",
        '<p class="rule">Hover any edge for its source anchor, its validity '
        "window and its evidence class; hover any node for its role and atom "
        "counts. Use the legend to dim a class of evidence &mdash; dimming "
        "everything but <strong>pinned</strong> shows you exactly how much of "
        "this network rests on constants nothing measured.</p>",
        _legend(),
    ]

    lane_edges: Dict[str, List[Dict[str, Any]]] = {
        lane: edge_records(lane) for lane in LANES
    }

    # ---- flux mode --------------------------------------------------------
    if flux_runs:
        parts.append("<h2>Flux mode &mdash; which chemistry dominates each process</h2>")
        parts.append(
            '<p class="note">Edge width below is the TIME-INTEGRATED ABSOLUTE '
            "FLUX through that step over the whole thermal program (mmol/L), on "
            "a log scale spanning five decades. A hairline carried almost "
            "nothing; a thick edge is where the material actually went. The "
            "same network looks completely different in each of these "
            "processes, which is the point.</p>"
        )
        for run in flux_runs:
            parts.append(f"<h3>{esc(run['name'])}</h3>")
            parts.append(
                '<div class="card"><table class="data"><tbody>'
                f"<tr><td>program</td><td>{esc(run['thermal'])}</td></tr>"
                f"<tr><td>pH / matrix</td><td>{num(run['ph'])} / "
                f"<code>{esc(run['matrix'])}</code></td></tr>"
                f"<tr><td>lanes</td><td>{esc(', '.join(run['lanes']) or 'none')}</td></tr>"
                f"<tr><td>envelope</td><td>{esc(run['state'])}</td></tr>"
                "</tbody></table>"
                + "".join(
                    f'<p class="note">refused: {esc(reason)}</p>'
                    for reason in run.get("reasons", [])
                )
                + "".join(
                    f'<p class="note">extrapolation: {esc(warning)}</p>'
                    for warning in run.get("warnings", [])
                )
                + "</div>"
            )
            if not run.get("flux"):
                parts.append(
                    '<p class="note">No flux is reported: the request is out of '
                    "envelope, and a refused request carries no numbers.</p>"
                )
                continue
            for lane in run["lanes"]:
                edges = lane_edges.get(lane) or []
                if not edges:
                    continue
                parts.append(f"<h3>{esc(lane)} lane, weighted by flux</h3>")
                parts.append(
                    '<div class="scroll">'
                    + _svg_lane(lane, edges, run["flux"])
                    + "</div>"
                )
                parts.append("<h3>Dominant steps</h3>")
                parts.append(_flux_table(edges, run["flux"]))

    # ---- the unweighted structural map ------------------------------------
    parts.append("<h2>The network itself, all four lanes</h2>")
    parts.append(
        '<p class="note">Unweighted: every step drawn at the same width, so this '
        "is the model's STRUCTURE rather than any one process's chemistry. The "
        "lanes do not compose &mdash; the sulfur steps are deliberately absent "
        "from the acrylamide lane, because composing them would spend the same "
        "cysteine twice.</p>"
        '<p class="note">A node with NO edges in a lane is a state variable that '
        "lane carries but never touches. The acrylamide lane is the clear case: "
        "it shares the sulfur lane's whole state vector while omitting every "
        "sulfur STEP, so every sulfur species appears there as an isolated node. "
        "That is the refusal to compose, drawn.</p>"
    )
    for lane in LANES:
        edges = lane_edges[lane]
        census = _evidence_census(edges)
        steps = len({e["key"] for e in edges})
        parts.append(f"<h2>{esc(lane)} lane</h2>")
        parts.append(
            f'<p class="note">{esc(metadata["lane_networks"].get(lane, ""))}</p>'
        )
        parts.append(
            f'<p class="note"><strong>{steps} steps</strong>, '
            + ", ".join(
                f"{count} {cls}" for cls, count in census.items() if count
            )
            + f'. {len(lane_species(lane))} species.</p>'
        )
        parts.append('<div class="scroll">' + _svg_lane(lane, edges, None) + "</div>")

    parts.append(
        "<footer>"
        f"<p>Generated by <code>scripts/generators/generate_network_map.py</code> "
        f"at commit <code>{esc(commit)}</code>. Deterministic: no timestamp, "
        "sorted iteration, fixed float formatting &mdash; two runs on the same "
        "tree produce byte-identical files.</p>"
        "<p>Parameters read from "
        + ", ".join(f"<code>{esc(p)}</code>" for p in metadata["parameters_from"])
        + ". This generator fits nothing and tunes nothing.</p>"
        "<p>Self-contained: no network requests, no external assets. The only "
        "script on the page is the inline legend filter, and the page is fully "
        "readable with JavaScript disabled.</p>"
        "</footer>"
    )

    return (
        "<!doctype html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        "<title>Maillard kinetic core - reaction network map</title>"
        f"<style>{CSS}</style></head><body><div class=\"wrap\">"
        + "".join(parts)
        + f'</div><script>{JS}</script></body></html>\n'
    )


# ===========================================================================
# 5. CLI
# ===========================================================================


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Render the kinetic core's reaction network from the live code, "
            "optionally weighted by the time-integrated flux of one or more "
            "process specs."
        )
    )
    parser.add_argument(
        "--spec",
        action="append",
        default=None,
        metavar="PATH.json",
        help="a formulation + process spec to run in FLUX MODE; repeatable",
    )
    parser.add_argument(
        "--all-examples",
        action="store_true",
        help=f"run the shipped example specs ({', '.join(DEFAULT_EXAMPLES)})",
    )
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="output HTML path")
    args = parser.parse_args(argv)

    spec_paths: List[Path] = [Path(p) for p in (args.spec or [])]
    if args.all_examples or not spec_paths:
        spec_paths = [EXAMPLES_DIR / name for name in DEFAULT_EXAMPLES]

    runs: List[Dict[str, Any]] = []
    for path in spec_paths:
        if not path.exists():
            raise SystemExit(f"spec not found: {path}")
        spec = load_flux_spec(path)
        spec["_source"] = str(path)
        runs.append(flux_budget(spec))

    page = build_page(runs, commit=_commit())
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    print(f"wrote {out} ({len(page)} bytes, {len(runs)} flux run(s))", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
