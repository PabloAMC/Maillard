"""The comparative front door: two formulations in, a ranked comparison out.

2026-09-03 (retirement step B5): this module drives ONE engine, the kinetic core
(``src/kinetic_core``). The FAST / "ordinal screening" lane that used to sit
beside it -- the SMIRKS enumeration path behind ``MaillardPipeline`` -- is
deleted, together with its ``--lane`` switch, its withheld-absolutes payload and
its directional-panel reliability tags (those were measured on the retired lane
and do not transfer; re-scoring the directional panel on the core is the next
retirement step).

What survives here is engine-neutral: spec loading and validation, the two-arm
document contract, ``predict_core`` / ``compare_core`` and their renderers, and
the value-of-information ``rank`` renderer over the core's own envelope.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from src import data_access, data_paths


RATIO_CAVEAT = (
    "Ratios are the reported quantity because comparisons cancel the systematic scale error. "
    "The interval is the CONSERVATIVE independent-error bound (A_p5/B_p95 .. A_p95/B_p5): it "
    "assumes the two arms' errors are unrelated, which the comparative argument says they are "
    "not. The correlation has never been measured, so no tighter interval is claimed. What HAS "
    "been measured is the direction, per axis, in the reliability column."
)


# ---------------------------------------------------------------------------------------
# Input specs
# ---------------------------------------------------------------------------------------

REQUIRED_SPEC_KEYS = ("precursors", "temp_C", "time_min", "ph", "aw")

SPEC_TEMPLATE = """\
# Maillard comparative-interface spec. Two arms, 'a' and 'b'.
# Precursor values are MILLIMOLAR (mM) -- the unit the projection budget consumes.
a:
  name: cysteine_ribose
  precursors:
    L-Cysteine: 10.0
    D-Ribose: 10.0
  temp_C: 140.0
  time_min: 30.0
  ph: 5.0
  aw: 0.98
  protein_type: free
b:
  name: cysteine_glucose
  precursors:
    L-Cysteine: 10.0
    D-Glucose: 10.0
  temp_C: 140.0
  time_min: 30.0
  ph: 5.0
  aw: 0.98
  protein_type: free
"""


class SpecError(ValueError):
    """A user-supplied spec is unusable. Never repaired silently."""


def load_spec_document(path: Path | str) -> Dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise SpecError(f"spec file not found: {path}")
    text = path.read_text(encoding="utf-8")
    try:
        data = yaml.safe_load(text)  # YAML is a JSON superset, so this reads .json too
    except yaml.YAMLError as exc:
        raise SpecError(f"{path}: could not parse as YAML/JSON -- {exc}") from exc
    if not isinstance(data, dict):
        raise SpecError(f"{path}: top level must be a mapping, got {type(data).__name__}")
    return data


def validate_spec(spec: Mapping[str, Any], *, label: str) -> Dict[str, Any]:
    if not isinstance(spec, Mapping):
        raise SpecError(f"{label}: must be a mapping")
    missing = [key for key in REQUIRED_SPEC_KEYS if spec.get(key) is None]
    if missing:
        raise SpecError(
            f"{label}: missing required key(s) {', '.join(missing)}. Every one of "
            f"{', '.join(REQUIRED_SPEC_KEYS)} must be given explicitly -- this interface does "
            "not default process conditions, because a defaulted temperature is a silent "
            "chemistry claim."
        )
    precursors = spec.get("precursors")
    if not isinstance(precursors, Mapping) or not precursors:
        raise SpecError(f"{label}: 'precursors' must be a non-empty mapping of name -> mM")
    resolved = dict(spec)
    resolved["name"] = str(spec.get("name") or label)
    resolved["precursors"] = {str(k): float(v) for k, v in precursors.items()}
    for key in ("temp_C", "time_min", "ph", "aw"):
        resolved[key] = float(spec[key])
    resolved["protein_type"] = str(spec.get("protein_type") or "free")
    return resolved


def split_comparison_document(
    document: Mapping[str, Any], *, source: str
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Pull the two arms out of a one-file comparison spec."""
    if "a" in document and "b" in document:
        return (
            validate_spec(document["a"], label=f"{source}:a"),
            validate_spec(document["b"], label=f"{source}:b"),
        )
    raise SpecError(
        f"{source}: a comparison spec needs top-level 'a:' and 'b:' arms (or pass two "
        "separate spec files). Write `maillard compare --template` for a starting point."
    )


def select_system(document: Mapping[str, Any], *, source: str, arm: Optional[str]) -> Dict[str, Any]:
    """Single-system selection, tolerating both flat specs and two-arm documents."""
    if arm:
        if arm not in document:
            raise SpecError(f"{source}: no system named '{arm}' in this spec")
        return validate_spec(document[arm], label=f"{source}:{arm}")
    if "a" in document and "b" in document:
        raise SpecError(
            f"{source}: this is a two-arm comparison spec; choose one with --system a|b, or "
            "use the `compare` verb."
        )
    return validate_spec(document, label=source)


# ---------------------------------------------------------------------------------------
# Running a spec through the existing pipeline
# ---------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------------------


def _ratio_bounds(
    a_env: Optional[Mapping[str, float]], b_env: Optional[Mapping[str, float]]
) -> Tuple[Optional[float], Optional[float]]:
    """Conservative independent-error bound on A/B. Not a correlation model."""
    if not a_env or not b_env:
        return (None, None)
    a_lo, a_hi = float(a_env["p5"]), float(a_env["p95"])
    b_lo, b_hi = float(b_env["p5"]), float(b_env["p95"])
    if b_hi <= 0 or b_lo <= 0:
        return (None, None)
    return (a_lo / b_hi, a_hi / b_lo)


# ---------------------------------------------------------------------------------------
# Single-system prediction
# ---------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------------------


def _fmt(value: Optional[float], places: int = 2) -> str:
    if value is None:
        return "-"
    if value != value:  # NaN
        return "-"
    if value == 0:
        return "0"
    if abs(value) >= 1e5 or abs(value) < 1e-3:
        return f"{value:.{places}e}"
    return f"{value:.{places}f}"


def _wrap(text: str, width: int = 92, indent: str = "  ") -> str:
    words, lines, current = text.split(), [], ""
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(indent + current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(indent + current)
    return "\n".join(lines)


def render_rank_text(payload: Mapping[str, Any]) -> str:
    out: List[str] = []
    out.append("=" * 96)
    out.append("  RANK-EXPERIMENTS   value-of-information over the uncertainty panel")
    out.append("=" * 96)
    out.append("")
    out.append(f"  source: {payload.get('source')}")
    out.append(
        f"  candidates: {payload.get('candidate_count')}   "
        f"outside the 90% CI: {payload.get('miss_count')}"
    )
    if payload.get("matrix_filter"):
        out.append(f"  matrix filter: {', '.join(payload['matrix_filter'])}")
    out.append("")
    out.append(
        f"  {'#':>3} {'compound':<28} {'VoI':>7} {'miss(log10)':>12}  {'template':<26} benchmark"
    )
    out.append("  " + "-" * 94)
    for candidate in payload.get("candidates", []):
        out.append(
            f"  {candidate['rank']:>3} {str(candidate['compound'])[:27]:<28} "
            f"{candidate['voi_score']:>7.3f} {candidate['envelope_miss_log10']:>12.3f}  "
            f"{str(candidate['suggested_doe_template'])[:25]:<26} {candidate['benchmark_id']}"
        )
    out.append("")
    out.append(
        _wrap(
            "This ranking is the model's honest product: every row is a place the model is "
            "measurably wrong, converted into a bookable measurement. It is computed from the "
            "cached Monte-Carlo envelope (results/validation/core_prediction_uncertainty.json), so it "
            "is only as current as the last trust-loop run."
        )
    )
    out.append("")
    return "\n".join(out)


def to_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, default=str)


# CORRECTED 2026-08-29 (Wave Q1). This caveat is printed to EVERY core-lane
# user, and three of its factual claims had been falsified by the waves that
# followed it: B6 added the lipid lane (so "no lipid-oxidation path" was wrong),
# B7 added HMF and DMHF as trunk targets (so "no HMF and no DMHF" was wrong),
# and B6's direct-sum co-integration means the lipid lane DOES compose with any
# one Maillard lane (so "its three lanes do not compose" was wrong, and there
# are four lanes). A caveat that overstates the model's limits is not the safe
# direction to be wrong in: it teaches users to distrust answers the model can
# actually support, and it goes stale invisibly because nothing tests prose.
def core_caveat() -> str:
    """The caveat printed with every core answer. Its one number is READ from the tracked
    scorecard (results/validation/core_panel_scores.json) rather than typed here: a
    hard-coded count was one wave behind twice (2026-09-03)."""
    try:
        summary = data_access.load_json(data_paths.CORE_PANEL_SCORES)["summary"]["out_of_sample"]
        headline = f"{summary['within_band']} of {summary['rows']} out-of-sample rows within 3x"
    except Exception:  # noqa: BLE001 - a missing scorecard is reported, not fatal
        headline = "scorecard not generated; run ./scripts/docker_maillard.sh core-scores"
    return (
        "KINETIC CORE. Absolute concentrations come from the mass-action network (frozen "
        "per-lane fit reports), and they are reported WITH their envelope declaration. The "
        "core refuses what it cannot name -- ask it for 1-hexanol, 2-pentylfuran or propanal and it "
        "will tell you why it will not answer, rather than answering. Its four lanes (trunk, sulfur, "
        "acrylamide, lipid) do not compose freely: the lipid lane co-integrates with ONE Maillard "
        "lane as a direct sum, and the Maillard lanes do not compose with each other. A refusal is "
        "an output, not a failure. Read the core's scorecard "
        f"(results/validation/core_panel_scores.md: {headline}) and its "
        "envelope (core_prediction_uncertainty.md) before trusting any number here; the "
        "directional panel (core_directional_scores.md) is what the reliability column reads."
    )


#: Evaluated once at import for callers that read the constant.
CORE_CAVEAT = core_caveat()

def _core_process(spec: Mapping[str, Any]):
    from src.kinetic_core.engine import ProcessSpec, ThermalProgram

    return ProcessSpec(
        thermal=ThermalProgram.isothermal(
            float(spec["temp_C"]), float(spec["time_min"])
        ),
        ph=float(spec["ph"]),
        water_activity=float(spec["aw"]) if spec.get("aw") is not None else None,
        matrix=str(spec.get("matrix") or spec.get("protein_type") or "water"),
    )


def spec_to_core(spec: Mapping[str, Any]):
    """Turn a validated CLI spec into the engine's FormulationSpec."""
    from src.kinetic_core.engine import FormulationSpec

    return FormulationSpec(
        name=str(spec.get("name") or "system"),
        precursors={str(k): float(v) for k, v in dict(spec["precursors"]).items()},
        process=_core_process(spec),
    )


def _core_targets(spec: Mapping[str, Any], targets: Optional[Sequence[str]]):
    """
    Which compounds to ask the core for.

    Build Wave V1 adds the middle branch: a spec may name ``targets:``
    explicitly. That is a usability fix with a governance benefit -- before it,
    a user could not ask "what does this model say about HMF here?", so the
    engine's most informative output (a NAMED REFUSAL with the reason) was
    unreachable from the front door. An explicitly named target that the core
    cannot represent now produces that refusal, which is the point of having
    written the refusals.
    """
    from src.kinetic_core.engine import default_targets_for

    if targets:
        return list(targets)
    declared = spec.get("targets")
    if declared:
        if isinstance(declared, str):
            return [declared]
        return [str(t) for t in declared]
    return list(default_targets_for(dict(spec["precursors"])))


def predict_core(
    spec: Mapping[str, Any], *, targets: Optional[Sequence[str]] = None
) -> Dict[str, Any]:
    """Single-formulation prediction through the kinetic core."""
    from src.kinetic_core.engine import engine_metadata, predict

    core_spec = spec_to_core(spec)
    requested = _core_targets(spec, targets)
    if not requested:
        return {
            "artifact": "maillard_predict_core",
            "lane": "core",
            "lane_label": "KINETIC CORE",
            "system": {"name": core_spec.name, "spec": dict(spec)},
            "answered": False,
            "declaration": {
                "state": "out_of_envelope",
                "lane": None,
                "reasons": [
                    "No precursor in this spec maps to a core species, so no lane "
                    "can be resolved and no target can be named. The core is a "
                    "named small-molecule network."
                ],
                "warnings": [],
                "unmapped_precursors": sorted(str(k) for k in spec["precursors"]),
                "unrepresented_targets": [],
                "mapped_precursors": {},
                "mapped_targets": {},
            },
            "rows": [],
            "caveats": {"core": CORE_CAVEAT},
            "engine": engine_metadata(),
        }

    run = predict(core_spec, requested)
    rows: List[Dict[str, Any]] = []
    oav_payload: Dict[str, Any] = {}
    if run.answered:
        oav_payload = dict(run.oav())
        # Q1: the rows come from the ENGINE, already carrying their interval and
        # their OAV entry, instead of being assembled here against a table keyed
        # in a different key-space. The hand-rolled version did
        # ``per_species[species_key]``, but ``per_species`` is keyed by B4 KEY --
        # identical for the trunk and sulfur lanes, DIFFERENT for the lipid one
        # ("HEXANAL" vs "hexanal", "DECADIENAL" vs "tt_2_4_decadienal"). So
        # every lipid compound's lookup missed and every report said "no
        # threshold" for two compounds that have a measured one. See
        # ``kinetic_core.keyspaces`` for why that failure is always silent.
        for row in run.interval_rows():
            rows.append(
                {
                    "compound": row["compound"],
                    "species_key": row["species_key"],
                    "predicted_ug_per_l": row["predicted_ug_per_l"],
                    # Q1: the interval travels WITH the point it belongs to. It
                    # used to be reachable only by joining the row back onto
                    # oav_table["per_species"][b4_key]["concentration"], which
                    # dropped it entirely for the four NO_B4_RECORD species --
                    # compounds with no odour threshold still have a perfectly
                    # well-defined concentration interval.
                    "interval_ug_per_l": row["interval_ug_per_l"],
                    "band_x": row["band_x"],
                    "oav": _oav_summary(row["oav"], no_b4_reason=row["no_b4_reason"]),
                    "lane": row["lane"],
                }
            )
    return {
        "artifact": "maillard_predict_core",
        "lane": "core",
        "lane_label": "KINETIC CORE",
        "system": {"name": core_spec.name, "spec": dict(spec)},
        "answered": run.answered,
        "declaration": run.declaration.as_dict(),
        "rows": rows,
        "oav_table": oav_payload,
        "run_metadata": dict(run.run_metadata),
        "caveats": {"core": CORE_CAVEAT},
        "engine": engine_metadata(),
    }


def _oav_summary(
    entry: Any, *, no_b4_reason: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Flatten one entry of the B4 ``oav_table``'s ``per_species`` block.

    A missing threshold is reported as ``available: False`` WITH its reason,
    never as an OAV of zero -- the B4 layer's whole design point is that the
    absence of a measured threshold is a state, not a value.

    Q1: ``no_b4_reason`` carries that same design point one step further out.
    A compound with NO B4 structural record has no entry in the table at all,
    and used to arrive here as a bare ``None`` -- so the report printed an
    empty reason, which reads as an unexplained blank rather than as the
    recorded, specific statement ``NO_B4_RECORD`` holds for each of them
    ("an alkane; the B4 registry has no hydrocarbon class at all").
    """
    if not isinstance(entry, Mapping):
        if no_b4_reason:
            return {
                "available": False,
                "threshold_state": "no_b4_structural_record",
                "reason": no_b4_reason,
            }
        return None
    interval = entry.get("OAV_interval") or [None, None]
    if entry.get("OAV_point") is None:
        return {
            "available": False,
            "threshold_state": entry.get("threshold_state"),
            "reason": (entry.get("diagnostics") or {}).get(
                "reason", "no measured threshold for this compound in this matrix"
            ),
        }
    return {
        "available": True,
        "oav": entry.get("OAV_point"),
        "low": interval[0],
        "high": interval[1],
        "threshold_ug_per_l": entry.get("threshold_ug_per_L"),
        "threshold_source": entry.get("threshold_source"),
        "threshold_state": entry.get("threshold_state"),
    }


def compare_core(
    spec_a: Mapping[str, Any],
    spec_b: Mapping[str, Any],
    *,
    targets: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Two-arm comparison through the kinetic core. Ratios lead, as in B4."""
    from src.kinetic_core.engine import compare as core_compare
    from src.kinetic_core.engine import engine_metadata

    core_a, core_b = spec_to_core(spec_a), spec_to_core(spec_b)
    requested = _core_targets(spec_a, targets) or _core_targets(spec_b, targets)
    payload = core_compare(core_a, core_b, requested) if requested else {
        "comparable": False,
        "reason": "no precursor in either arm maps to a core species",
        "declaration_a": {},
        "declaration_b": {},
    }
    return {
        "artifact": "maillard_compare_core",
        "lane": "core",
        "lane_label": "KINETIC CORE",
        "a": {"name": core_a.name, "spec": dict(spec_a)},
        "b": {"name": core_b.name, "spec": dict(spec_b)},
        "comparison": payload,
        "reliability": _reliability_block(spec_a, spec_b),
        "caveats": {"core": CORE_CAVEAT, "ratio": RATIO_CAVEAT},
        "engine": engine_metadata(),
    }


def _reliability_block(spec_a: Mapping[str, Any], spec_b: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Step B7 (2026-09-03): the per-axis reliability of THIS comparison, read from the
    core's own directional scorecard. The axes are the knobs the two arms differ on; the
    governing verdict is the weakest of them. When the scorecard is absent the block says
    so instead of inventing a tag.
    """
    from src.directional_reliability import describe_comparison

    try:
        described = describe_comparison(spec_a, spec_b)
    except (FileNotFoundError, ValueError) as exc:
        return {"available": False, "reason": str(exc), "axes": [], "per_axis": [], "governing": None}
    governing = described["governing"]
    return {
        "available": True,
        "source": "results/validation/core_directional_scores.json (strictly independent claims)",
        "axes": list(described["axes"]),
        "per_axis": [
            {"axis": r.axis, "agree": r.agree, "evaluable": r.evaluable, "verdict": r.verdict, "note": r.note}
            for r in described["per_axis"]
        ],
        "governing": None if governing is None else governing.render(),
        "governing_verdict": None if governing is None else governing.verdict,
        "no_axis_differs": bool(described["no_axis_differs"]),
    }


def _render_declaration(declaration: Mapping[str, Any], out: List[str]) -> None:
    state = declaration.get("state", "unknown")
    lanes = declaration.get("lanes") or (
        [declaration["lane"]] if declaration.get("lane") else []
    )
    out.append(
        f"  ENVELOPE: {state.upper()}   lane: "
        + (" + ".join(str(lane) for lane in lanes) if lanes else "NONE RESOLVED")
    )
    mapped = dict(declaration.get("mapped_precursors") or {})
    if mapped:
        out.append(
            "  mapped precursors: "
            + ", ".join(f"{key}={value:g} mM" for key, value in sorted(mapped.items()))
        )
    for reason in declaration.get("reasons") or []:
        out.append(_wrap(f"REFUSED -- {reason}", indent="    ! "))
    for warning in declaration.get("warnings") or []:
        out.append(_wrap(f"declared extrapolation -- {warning}", indent="    ~ "))


def envelope_error_text(declaration: Mapping[str, Any]) -> str:
    """
    The out-of-envelope message, naming the missing lane and the missing species.

    Build Wave V1. The engine already writes an excellent reason for every
    refusal; what it does not do is put the two questions a bench scientist
    actually asks -- "which lane could not be resolved?" and "which of MY inputs
    was the problem?" -- at the top. This composes both from the
    ``EnvelopeDeclaration`` itself, reusing its text verbatim and inventing
    nothing.
    """
    lines: List[str] = []
    lanes = declaration.get("lanes") or (
        [declaration["lane"]] if declaration.get("lane") else []
    )
    lines.append(
        "OUT OF ENVELOPE -- no number is emitted. "
        + (
            f"Lane resolved: {' + '.join(str(lane) for lane in lanes)}."
            if lanes
            else "NO LANE could be resolved for this request."
        )
    )
    unmapped = list(declaration.get("unmapped_precursors") or [])
    if unmapped:
        lines.append(
            "  missing species (precursors the core cannot name): "
            + ", ".join(str(u) for u in unmapped)
        )
    unrepresented = declaration.get("unrepresented_targets") or []
    if unrepresented:
        names = [
            (entry.get("compound") if isinstance(entry, Mapping) else entry[0])
            for entry in unrepresented
        ]
        lines.append(
            "  missing species (targets the core cannot name): "
            + ", ".join(str(n) for n in names)
        )
    for reason in declaration.get("reasons") or []:
        lines.append(_wrap(f"declared reason -- {reason}", indent="  "))
    lines.append(
        "  A refusal is an output, not a failure. Run "
        "`python scripts/maillard.py explain <compound>` to see what the core "
        "does carry for a compound, and why."
    )
    return "\n".join(lines)


def render_predict_core_text(payload: Mapping[str, Any]) -> str:
    out: List[str] = []
    out.append("=" * 96)
    out.append(f"  PREDICT [KINETIC CORE]   {payload['system']['name']}")
    out.append("=" * 96)
    out.append("")
    _render_declaration(payload.get("declaration") or {}, out)
    out.append("")
    if not payload.get("answered"):
        out.append("  NO NUMBER IS EMITTED. The core declined this request above.")
        out.append("")
        out.append(_wrap(payload["caveats"]["core"]))
        out.append("")
        return "\n".join(out)

    out.append(f"  {'compound':<38} {'ug/L (= ppb in water)':>24}  {'OAV':>14}")
    out.append("  " + "-" * 82)
    for row in payload["rows"]:
        oav = row.get("oav") or {}
        oav_text = (
            f"{oav['oav']:.3g}" if oav.get("available") and oav.get("oav") is not None
            else "no threshold"
        )
        out.append(
            f"  {row['compound'][:37]:<38} {_fmt(row['predicted_ug_per_l']):>24}  {oav_text:>14}"
        )
    out.append("")
    out.append(_wrap(payload["caveats"]["core"]))
    out.append("")
    return "\n".join(out)


def render_compare_core_text(payload: Mapping[str, Any]) -> str:
    out: List[str] = []
    a_name, b_name = payload["a"]["name"], payload["b"]["name"]
    out.append("=" * 96)
    out.append(f"  COMPARE [KINETIC CORE]   A = {a_name}   vs   B = {b_name}")
    out.append("=" * 96)
    out.append("")
    comparison = payload.get("comparison") or {}
    if not comparison.get("comparable"):
        out.append("  NOT COMPARABLE.")
        out.append(_wrap(str(comparison.get("reason", "")), indent="    ! "))
        for label, key in (("A", "declaration_a"), ("B", "declaration_b")):
            declaration = comparison.get(key) or {}
            if declaration:
                out.append(f"  arm {label}:")
                _render_declaration(declaration, out)
        out.append("")
        out.append(_wrap(payload["caveats"]["core"]))
        out.append("")
        return "\n".join(out)

    # Build Wave V1. The declarations used to be printed ONLY when the two arms
    # could not be compared, so a successful comparison silently dropped every
    # declared extrapolation the engine had attached to it -- the buffer that
    # was never declared, the lipid rate anchored 12 decades away. A warning
    # that appears only on failure is a warning nobody reads.
    for label, key in (("A", "declaration_a"), ("B", "declaration_b")):
        declaration = comparison.get(key) or {}
        if declaration:
            out.append(f"  arm {label}:")
            _render_declaration(declaration, out)
            out.append("")

    reliability = payload.get("reliability") or {}
    if reliability.get("available"):
        axes = reliability.get("axes") or []
        out.append("  Axes this comparison moves: " + (", ".join(axes) if axes else "none (the arms differ on no panel axis)"))
        out.append(f"  Governing reliability:      {reliability.get('governing') or 'n/a'}   [directional panel, independent claims]")
        for item in reliability.get("per_axis") or []:
            out.append(f"    - {item['axis']:<18} {item['verdict']:<12} panel {item['agree']}/{item['evaluable']}")
            if item.get("note"):
                out.append(_wrap(item["note"], indent="        "))
        out.append("")
    elif reliability:
        out.append(_wrap("reliability tags unavailable: " + str(reliability.get("reason")), indent="  ! "))
        out.append("")

    ratios = comparison.get("ratios") or {}
    rows = list(ratios.get("rows") or [])
    band = ratios.get("reliability_band_x")
    out.append(
        f"  {ratios.get('n_resolved', 0)} of {ratios.get('n_compared', 0)} ratios resolve "
        f"above the same-sample dispersion band"
        + (f" ({band:.3g}x)" if isinstance(band, (int, float)) else "")
    )
    # Q1: fixed at the source. ``compare_formulations`` no longer counts an
    # UNDEFINED ratio (one arm at exactly zero) as resolved, and publishes
    # ``n_undefined`` so this line does not have to re-derive the predicate --
    # the V1 workaround that lived here did, and it disagreed with the wider
    # predicate the per-row column below uses.
    undefined = int(
        ratios.get(
            "n_undefined",
            sum(1 for row in rows if str(row.get("direction")) == "undefined"),
        )
    )
    if undefined:
        out.append(
            f"  ...of which {undefined} are UNDEFINED (one arm at exactly zero, or on a route the "
            f"engine declares unidentified) and resolve nothing: see the 'resolved' column, not this count."
        )
    out.append("")
    out.append(f"  {'compound':<38} {'A/B':>12} {'direction':<14} resolved")
    out.append("  " + "-" * 88)
    for row in rows:
        ratio = row.get("ratio_a_over_b")
        if row.get("unidentified_arm"):
            ratio_text = f"{str(row['unidentified_arm']).upper()} unidentified"
        elif not isinstance(ratio, (int, float)):
            ratio_text = "-"
        elif ratio != ratio or ratio in (float("inf"), float("-inf")) or ratio == 0.0:
            # One arm is at zero: a ratio is undefined, not enormous.
            ratio_text = "A only" if ratio else "B only"
        else:
            ratio_text = f"{ratio:.3g}x"
        # 'within the band' means NOT resolved -- a ratio inside the dispersion
        # band is reported as unresolved rather than as a small effect. An
        # UNDEFINED ratio (one arm at zero) is neither: V1 stops labelling it
        # 'yes', which read as a resolved claim about a quantity that does not
        # exist.
        if str(row.get("direction")) == "undefined" or not isinstance(ratio, (int, float)) \
                or ratio != ratio or ratio in (float("inf"), float("-inf")) or ratio == 0.0:
            resolved = "n/a -- undefined"
        elif row.get("within_reliability_band"):
            resolved = "no -- inside band"
        else:
            resolved = "yes"
        out.append(
            f"  {str(row.get('compound', '?'))[:37]:<38} {ratio_text:>12} "
            f"{str(row.get('direction', '-')):<14} {resolved}"
        )
    out.append("")
    out.append(_wrap(payload["caveats"]["ratio"]))
    out.append("")
    out.append(_wrap(payload["caveats"]["core"]))
    out.append("")
    return "\n".join(out)
