"""
src/explain_compound.py -- ``maillard explain <compound>`` (Build Wave V1).

WHAT IT ANSWERS
---------------
"Where does this compound come from in your model, how good is each step's
evidence, and which papers am I actually trusting?"

It answers that from the LIVE registries -- the reaction tuples in
``network.py`` / ``sulfur.py`` / ``acrylamide.py``, the parameter registries
behind them, and the B4 threshold table -- so it moves the day the model moves.
NO NEW DATA is introduced here; every string printed comes out of a registry.

A REFUSAL IS AN ANSWER, HERE TOO
--------------------------------
Asking about HMF, DMHF, 2-pentylfuran or 1-hexanol does not produce an error.
It produces the engine's own declared reason for refusing them, verbatim, which
is more useful to a bench scientist than a stack trace: it says what would have
to be measured for the answer to exist.

EVIDENCE CLASSES, NORMALISED TO FOUR WORDS
------------------------------------------
The registries use their own vocabulary (``measured_rate``,
``measured_activation_energy``, ``derived_from_fit_data``,
``bounded_from_a_timescale_bracket``, ``structural_constant``,
``declared_assumption``). Downstream surfaces -- this command and the network
map -- collapse them to four, ONCE, here:

  measured -- a rate or an Ea printed in a paper
  fitted   -- estimated in this repository by least squares on FIT rows
  derived  -- computed from another constant by a stated relation
  pinned   -- held at a value nothing measured (including held at zero)

``evidence_class_of`` is the single place that mapping happens; the raw class
travels alongside it so nothing is lost.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

MEASURED = "measured"
FITTED = "fitted"
DERIVED = "derived"
PINNED = "pinned"

EVIDENCE_ORDER = (MEASURED, FITTED, DERIVED, PINNED)

EVIDENCE_MEANING: Mapping[str, str] = {
    MEASURED: "a rate constant or activation energy printed in a paper",
    FITTED: "estimated in this repository by least squares on declared FIT rows",
    DERIVED: "computed from another constant by a stated relation, not free",
    PINNED: "held at a value nothing measured -- including held at zero",
}


def evidence_class_of(
    raw_class: Optional[str], flags: Sequence[str] = ()
) -> str:
    """Collapse a registry evidence class to one of the four display classes."""
    flags = tuple(flags or ())
    if raw_class in ("measured_rate", "measured_activation_energy"):
        return MEASURED
    if raw_class == "derived_from_fit_data":
        return FITTED if "fitted_here" in flags else DERIVED
    if raw_class in ("bounded_from_a_timescale_bracket", "declared_assumption"):
        return PINNED
    if raw_class in ("structural_constant", "measured_ratio"):
        return DERIVED
    if raw_class is None:
        return DERIVED
    return PINNED


# ---------------------------------------------------------------------------
# Lane wiring
# ---------------------------------------------------------------------------


def _lane_reactions(lane: str) -> Tuple[Any, ...]:
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


def _lane_parameters(lane: str) -> Dict[str, Any]:
    from src.kinetic_core.engine import core_parameters

    if lane == "lipid":
        return {}
    return dict(core_parameters(lane))


def _species_label(key: str) -> str:
    """A human name for a species key, from whichever table owns it."""
    from src.kinetic_core.species import BY_KEY

    if key in BY_KEY:
        return BY_KEY[key].label
    from src.kinetic_core.species_sulfur import SULFUR_SPECIES

    for species in SULFUR_SPECIES:
        if species.key == key:
            return species.label
    from src.kinetic_core.species_acrylamide import ACRYLAMIDE_SPECIES

    for species in ACRYLAMIDE_SPECIES:
        if species.key == key:
            return species.label
    from src.kinetic_core.species_lipid import LIPID_SPECIES

    for species in LIPID_SPECIES:
        if species.key == key:
            return species.label
    return key


def _side(mapping: Mapping[str, int]) -> str:
    return " + ".join(
        (f"{coefficient} " if coefficient != 1 else "") + key
        for key, coefficient in mapping.items()
    ) or "-"


def _parameter_record(
    reaction: Any, parameters: Mapping[str, Any], lane: str
) -> Dict[str, Any]:
    """Everything defensible about the rate on one step."""
    key = reaction.parameter_key
    if key is None:
        return {
            "parameter_key": None,
            "evidence": DERIVED,
            "raw_evidence_class": None,
            "source_anchor": (
                "no parameter of its own: this step's constant is DERIVED from "
                "another step's by a pinned relation stated in the network module"
            ),
            "note": reaction.note or "",
            "flags": [],
            "temperature_range_c": None,
            "ph_of_measurement": None,
            "rate_transfer": "derived",
            "ea_kj_mol": None,
        }
    if lane == "sulfur" and key == "k_oligomer":
        from src.kinetic_core.parameters_sulfur import NO_MEASURED_RATE

        return {
            "parameter_key": key,
            "evidence": PINNED,
            "raw_evidence_class": "no_measured_rate",
            "source_anchor": NO_MEASURED_RATE[key],
            "note": "HELD AT 0.0 -- the channel is structurally present and switched off.",
            "flags": ["pinned_at_zero", "holdout_protected"],
            "temperature_range_c": None,
            "ph_of_measurement": None,
            "rate_transfer": "not_licensed",
            "ea_kj_mol": None,
        }
    parameter = parameters.get(key)
    if parameter is None:
        return {
            "parameter_key": key,
            "evidence": PINNED,
            "raw_evidence_class": None,
            "source_anchor": f"no parameter {key!r} in this lane's operative set",
            "note": "",
            "flags": [],
            "temperature_range_c": None,
            "ph_of_measurement": None,
            "rate_transfer": "unknown",
            "ea_kj_mol": None,
        }
    flags = list(getattr(parameter, "flags", ()) or ())
    raw = getattr(parameter, "evidence_class", None)
    return {
        "parameter_key": key,
        "evidence": evidence_class_of(raw, flags),
        "raw_evidence_class": raw,
        "transformation": getattr(parameter, "transformation", "") or "",
        "source_anchor": getattr(parameter, "source_anchor", "") or "",
        "dossier_anchor": getattr(parameter, "dossier_anchor", "") or "",
        "conditions": getattr(parameter, "conditions", "") or "",
        "note": getattr(parameter, "note", "") or "",
        "flags": flags,
        "temperature_range_c": list(getattr(parameter, "temperature_range_c", ()) or ()) or None,
        "ph_of_measurement": getattr(parameter, "ph_of_measurement", None),
        "rate_transfer": getattr(parameter, "rate_transfer", "unknown"),
        "ea_kj_mol": getattr(parameter, "ea_kj_mol", None),
        "k_ref": getattr(parameter, "k_ref", None),
        "unit": getattr(parameter, "unit", ""),
    }


# ---------------------------------------------------------------------------
# The lipid lane, which has routes but no mass-action reaction tuple
# ---------------------------------------------------------------------------


def _lipid_routes(species_key: str) -> List[Dict[str, Any]]:
    from src.kinetic_core.engine import core_lipid_model
    from src.kinetic_core.parameters_lipid import (
        K_LOOH_DECOMP_ANCHOR,
        Q10_ASSUMPTION,
    )
    from src.kinetic_core.species_lipid import (
        CLEAVAGE_MECHANISM,
        LOOH_POOLS,
        POSITION_PRODUCTS,
    )

    branch, composition = core_lipid_model()
    routes: List[Dict[str, Any]] = []
    for pool, (position, geometry) in LOOH_POOLS.items():
        if species_key not in POSITION_PRODUCTS.get(position, ()):
            continue
        share = float(branch.simplexes[(position, geometry)].get(species_key, 0.0))
        routes.append(
            {
                "step": f"{pool} -> {species_key}",
                "equation": f"{pool} -> {species_key} (+ co-products)",
                "branch_share_of_this_pool": share,
                "pool_fraction": float(composition.as_dict().get(pool, 0.0)),
                "evidence": FITTED,
                "raw_evidence_class": "derived_from_fit_data",
                "source_anchor": branch.provenance,
                "note": (
                    f"cleavage mechanism: {CLEAVAGE_MECHANISM.get(species_key, 'unassigned')}"
                    " (homolytic routes are suppressed by a hydrogen donor; "
                    "'both' routes are not). The DISTRIBUTION is fitted; the RATE "
                    "below is an assumption."
                ),
                "flags": ["branch_distribution_is_fit", "structural_zeros_enforced"],
                "temperature_range_c": None,
                "ph_of_measurement": None,
                "rate_transfer": "licensed",
            }
        )
    routes.append(
        {
            "step": "LOOH pool -> products (the RATE)",
            "equation": "first-order hydroperoxide decomposition",
            "evidence": PINNED,
            "raw_evidence_class": "declared_assumption",
            "source_anchor": K_LOOH_DECOMP_ANCHOR.source_anchor,
            "note": Q10_ASSUMPTION.warning,
            "flags": ["rate_is_an_assumption", "q10_band_2_to_3"],
            "temperature_range_c": [
                Q10_ASSUMPTION.licensed_span_c[0],
                Q10_ASSUMPTION.licensed_span_c[1],
            ],
            "ph_of_measurement": K_LOOH_DECOMP_ANCHOR.ph_of_measurement,
            "rate_transfer": K_LOOH_DECOMP_ANCHOR.rate_transfer,
        }
    )
    return routes


# ---------------------------------------------------------------------------
# The public entry point
# ---------------------------------------------------------------------------


def explain(compound: str) -> Dict[str, Any]:
    """Formation routes, evidence classes and literature anchors for ``compound``."""
    from src.kinetic_core.engine import (
        LANE_DEFAULT_TARGETS,
        TARGET_ALIASES,
        TRUNK,
        UNREPRESENTED_COMPOUNDS,
        _TARGET_LANE,
        _norm,
    )

    query = _norm(compound)
    payload: Dict[str, Any] = {
        "artifact": "maillard_explain",
        "query": str(compound),
        "known_targets": sorted({t for ts in LANE_DEFAULT_TARGETS.values() for t in ts}),
        "evidence_meaning": dict(EVIDENCE_MEANING),
    }

    if query in UNREPRESENTED_COMPOUNDS:
        payload.update(
            {
                "answered": False,
                "state": "refused",
                "reason": UNREPRESENTED_COMPOUNDS[query],
                "species_key": None,
                "lane": None,
                "routes": [],
            }
        )
        return payload

    species_key = TARGET_ALIASES.get(query)
    if species_key is None:
        payload.update(
            {
                "answered": False,
                "state": "no_vocabulary_entry",
                "reason": (
                    f"{compound!r} is not a species in any core lane, and it is "
                    "not on the named unrepresented-compound list either: the "
                    "engine has no vocabulary entry for it. That is a gap in the "
                    "vocabulary, not a claim about the chemistry."
                ),
                "species_key": None,
                "lane": None,
                "routes": [],
                "aliases_known": sorted(TARGET_ALIASES),
            }
        )
        return payload

    lane = _TARGET_LANE.get(species_key, TRUNK)
    payload.update({"species_key": species_key, "lane": lane, "answered": True,
                    "state": "answered", "label": _species_label(species_key)})

    if lane == "lipid":
        routes = _lipid_routes(species_key)
        payload["routes"] = routes
        payload["consumption"] = []
        payload["anchors"] = _rank_anchors(routes)
        payload["threshold"] = _threshold_block(species_key)
        return payload

    reactions = _lane_reactions(lane)
    parameters = _lane_parameters(lane)
    formation: List[Dict[str, Any]] = []
    consumption: List[Dict[str, Any]] = []
    for reaction in reactions:
        record = None
        if species_key in reaction.products:
            record = _parameter_record(reaction, parameters, lane)
            formation.append(
                {
                    "step": reaction.key,
                    "equation": f"{_side(reaction.reactants)} -> {_side(reaction.products)}",
                    "stoichiometry": int(reaction.products[species_key]),
                    **record,
                }
            )
        if species_key in reaction.reactants:
            record = record or _parameter_record(reaction, parameters, lane)
            consumption.append(
                {
                    "step": reaction.key,
                    "equation": f"{_side(reaction.reactants)} -> {_side(reaction.products)}",
                    **record,
                }
            )

    payload["routes"] = formation
    payload["consumption"] = consumption
    payload["anchors"] = _rank_anchors(formation + consumption)
    payload["threshold"] = _threshold_block(species_key)
    payload["evidence_census"] = {
        cls: sum(1 for r in formation if r["evidence"] == cls) for cls in EVIDENCE_ORDER
    }
    return payload


def _rank_anchors(routes: Sequence[Mapping[str, Any]], top_n: int = 6) -> List[Dict[str, Any]]:
    """
    The literature anchors behind these routes, best evidence first.

    Ranked by evidence class and then by how many steps each anchor carries, so
    the paper the model leans on hardest sorts to the top. NOTHING is summarised
    or paraphrased: the anchor string is the registry's own.
    """
    by_anchor: Dict[str, Dict[str, Any]] = {}
    for route in routes:
        anchor = str(route.get("source_anchor") or "").strip()
        if not anchor:
            continue
        entry = by_anchor.setdefault(
            anchor,
            {"source_anchor": anchor, "steps": [], "evidence": route.get("evidence", PINNED)},
        )
        entry["steps"].append(route.get("step"))
        if EVIDENCE_ORDER.index(route.get("evidence", PINNED)) < EVIDENCE_ORDER.index(
            entry["evidence"]
        ):
            entry["evidence"] = route.get("evidence")
    ranked = sorted(
        by_anchor.values(),
        key=lambda e: (EVIDENCE_ORDER.index(e["evidence"]), -len(e["steps"]), e["source_anchor"]),
    )
    return ranked[:top_n]


def _threshold_block(species_key: str) -> Dict[str, Any]:
    """The compound's water odour threshold, if the B4 table has one."""
    from src.kinetic_core.matrix_oav import NoMeasuredThreshold, select_threshold_verbose
    from src.kinetic_core.species_lipid import B4_COMPOUND_KEY, NO_B4_RECORD

    if species_key in NO_B4_RECORD:
        return {"available": False, "reason": NO_B4_RECORD[species_key]}
    b4 = B4_COMPOUND_KEY.get(species_key, species_key)
    record, diagnostics = select_threshold_verbose(b4, "water")
    if isinstance(record, NoMeasuredThreshold):
        return {"available": False, "reason": record.reason}
    return {
        "available": True,
        "value_ug_per_l": record.value_ug_per_l,
        "source": record.source,
        "criterion": record.criterion,
        "provenance_flag": record.provenance_flag,
        "notes": record.notes,
        "spread_x": diagnostics.get("spread_x"),
        "alternates_ug_per_l": diagnostics.get("alternates_ug_per_l"),
    }


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _wrap(text: str, width: int = 88, indent: str = "    ") -> str:
    words, lines, current = str(text).split(), [], ""
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(indent + current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(indent + current)
    return "\n".join(lines)


def render_explain_text(payload: Mapping[str, Any]) -> str:
    out: List[str] = []
    out.append("=" * 96)
    out.append(f"  EXPLAIN   {payload['query']}")
    out.append("=" * 96)
    out.append("")

    if not payload.get("answered"):
        state = payload.get("state")
        header = (
            "  REFUSED -- this compound is on the engine's named "
            "unrepresented-compound list."
            if state == "refused"
            else "  NO VOCABULARY ENTRY."
        )
        out.append(header)
        out.append("")
        out.append(_wrap(str(payload.get("reason", "")), indent="    "))
        out.append("")
        out.append("  A refusal is an answer: it says what would have to be measured.")
        out.append("")
        out.append("  Compounds this model CAN explain:")
        for target in payload.get("known_targets", []):
            out.append(f"    - {target}")
        out.append("")
        return "\n".join(out)

    out.append(
        f"  species key: {payload['species_key']}   lane: {payload['lane']}   "
        f"{payload.get('label', '')}"
    )
    census = payload.get("evidence_census") or {}
    if census:
        out.append(
            "  formation steps by evidence: "
            + ", ".join(f"{cls} {count}" for cls, count in census.items() if count)
        )
    out.append("")
    out.append("  EVIDENCE CLASSES")
    for cls in EVIDENCE_ORDER:
        out.append(f"    {cls:<9} {EVIDENCE_MEANING[cls]}")
    out.append("")

    routes = list(payload.get("routes") or [])
    out.append(f"  FORMATION ROUTES ({len(routes)})")
    out.append("  " + "-" * 94)
    if not routes:
        out.append("    (none -- this species is charged, never formed, in this lane)")
    for route in routes:
        out.append(
            f"    [{route.get('evidence', '?'):<8}] {str(route.get('step', '?')):<26} "
            f"{route.get('equation', '')}"
        )
        detail = []
        if route.get("parameter_key"):
            detail.append(f"k={route['parameter_key']}")
        if route.get("raw_evidence_class"):
            detail.append(f"class={route['raw_evidence_class']}")
        if route.get("temperature_range_c"):
            low, high = route["temperature_range_c"]
            detail.append(f"valid {low:g}-{high:g} C")
        if route.get("ph_of_measurement") is not None:
            detail.append(f"measured at pH {route['ph_of_measurement']:g}")
        if route.get("rate_transfer"):
            detail.append(f"transfer={route['rate_transfer']}")
        if route.get("branch_share_of_this_pool") is not None:
            detail.append(f"branch share {route['branch_share_of_this_pool']:.3g}")
        if detail:
            out.append("        " + "  ".join(detail))
        if route.get("source_anchor"):
            out.append(_wrap(f"src: {route['source_anchor']}", indent="        "))
        if route.get("note"):
            # The note is where a step says the thing that would otherwise be
            # invisible -- that its rate is an assumption, that its product was
            # never identified, that its channel is held at zero. Printing the
            # anchor without it would show the paper and hide the caveat.
            out.append(_wrap(f"note: {route['note']}", indent="        "))
        for flag in route.get("flags", []):
            out.append(f"        ! {flag}")
    out.append("")

    consumption = list(payload.get("consumption") or [])
    if consumption:
        out.append(f"  CONSUMPTION ROUTES ({len(consumption)})")
        out.append("  " + "-" * 94)
        for route in consumption:
            out.append(
                f"    [{route.get('evidence', '?'):<8}] {str(route.get('step', '?')):<26} "
                f"{route.get('equation', '')}"
            )
        out.append("")

    anchors = list(payload.get("anchors") or [])
    if anchors:
        out.append(f"  TOP LITERATURE ANCHORS ({len(anchors)}), best evidence first")
        out.append("  " + "-" * 94)
        for index, anchor in enumerate(anchors, start=1):
            out.append(
                f"    {index}. [{anchor['evidence']}] carries "
                f"{len(anchor['steps'])} step(s): {', '.join(str(s) for s in anchor['steps'][:6])}"
            )
            out.append(_wrap(anchor["source_anchor"], indent="       "))
        out.append("")

    threshold = dict(payload.get("threshold") or {})
    out.append("  ODOUR THRESHOLD (water) -- a model INPUT, never a prediction")
    if threshold.get("available"):
        out.append(f"    {threshold['value_ug_per_l']:g} ug/L")
        out.append(_wrap(f"src: {threshold['source']}", indent="       "))
        if threshold.get("provenance_flag"):
            out.append(f"       flag: {threshold['provenance_flag']}")
        spread = threshold.get("spread_x")
        if isinstance(spread, (int, float)) and spread > 1.0:
            out.append(
                f"       {len(threshold.get('alternates_ug_per_l') or [])} candidate "
                f"values spanning {spread:.3g}x -- carried, never averaged away"
            )
    else:
        out.append(_wrap(str(threshold.get("reason", "no measured threshold")), indent="    "))
    out.append("")
    out.append(
        "  Nothing on this page is new data: every line is read from a frozen "
        "registry at run time."
    )
    out.append("")
    return "\n".join(out)


__all__ = [
    "DERIVED",
    "EVIDENCE_MEANING",
    "EVIDENCE_ORDER",
    "FITTED",
    "MEASURED",
    "PINNED",
    "evidence_class_of",
    "explain",
    "render_explain_text",
]
