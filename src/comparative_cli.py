"""The comparative interface: orchestration only, no science.

WHY THIS EXISTS
---------------
2026-08-28 (Wave S5). Every accuracy measurement this repository has made says the same
thing twice: the **absolute** numbers are bad and the **ordinal** ones are usable. The
free-precursor hold-out lands at a 6.04x median fold error; the matrix lane's three
observability modes land at 67-94x; CI coverage on genuine extrapolation is 1/5. Against
that, the directional panel scores 24/35 on strictly independent claims and 18/23 once pH
and water activity are set aside (2026-08-28, Wave W: was 21/29 and 17/19 before six
interlibrary-loan rows were added and scored 3/6 -- the ordinal claim got WEAKER, not
stronger, and the numbers here move with it).

So the front door leads with **ratios**, not with ppb. That is not a presentational
preference, it is the only reading of the model the evidence supports: a comparison of two
arms run through the same projection budget, the same marker yields and the same
observability factors cancels the systematic scale error those constants carry, and what
survives is the thing that was measured to work.

CORRECTED 2026-08-29 (Wave Q1). This line used to read "Absolutes are still available behind
``--absolute``, and they print their own caveat." That has not been true since the B5 cutover.
On the FAST lane ``--absolute`` now EXITS 2, and ``screening_payload`` strips every ppb field
from the payload before it can reach a renderer, unconditionally -- so no FAST absolute reaches
a user by any route. On the CORE lane absolutes ARE emitted, always with their envelope
declaration and always inside a reliability interval, and ``--absolute`` is a no-op note there
because the flag has nothing left to unlock.

WHAT THIS MODULE MAY AND MAY NOT DO
-----------------------------------
May: call ``MaillardPipeline``, reuse ``benchmark_to_formulation`` /
``benchmark_to_conditions`` for input handling, read the directional artifact through
``src.directional_reliability``, arrange the results in a table.

May NOT: introduce a constant, a correction, a correlation model, or any number that is not
either an input, a model output, or arithmetic on the two. In particular the ratio interval
below is the ordinary independent-error combination of two existing envelopes and is labelled
as the conservative bound it is; this module does not model the error correlation that the
comparative thesis rests on, because nobody has measured it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation
from src.config import DEFAULTS
from src.directional_reliability import (
    AxisReliability,
    compound_axes,
    describe_comparison,
    is_sulfur_compound,
    load_panel_counts,
    reliability_for_axis,
    weakest,
)
from src.pipeline import MaillardPipeline
from src.usability_reports import prepare_cli_confidence

#: The FAST lane's absolute caveat. CORRECTED 2026-08-29 (Wave Q1): this said
#: "Printed above every absolute number this interface emits, without
#: exception", which stopped being true at the B5 cutover. It is now printed
#: NOWHERE on the FAST lane -- ``screening_payload`` overwrites the caveat block
#: with the screening label and strips the ppb fields -- and it is not used on
#: the core lane at all, which has its own ``CORE_CAVEAT``. It is reachable only
#: from the library-level renderers, which the unit tests call directly on raw
#: payloads to pin that IF an absolute is ever printed it arrives with this text
#: in the same block. Keep it for that; do not read it as shipped behaviour.
ABSOLUTE_CAVEAT = (
    "ABSOLUTE ppb ARE NOT RELIABLE. Measured out-of-sample: median 6.0x fold error on the "
    "free-precursor hold-out (worst 52.6x), 67-94x on the matrix lane, and 1 of 5 genuine "
    "extrapolation rows inside the 90% CI. Read these as order-of-magnitude hypotheses that "
    "tell you which experiment to run, never as a specification."
)

RATIO_CAVEAT = (
    "Ratios are the reported quantity because comparisons cancel the systematic scale error. "
    "The interval is the CONSERVATIVE independent-error bound (A_p5/B_p95 .. A_p95/B_p5): it "
    "assumes the two arms' errors are unrelated, which the comparative argument says they are "
    "not. The correlation has never been measured, so no tighter interval is claimed. What HAS "
    "been measured is the direction, per axis, in the reliability column."
)

SULFUR_CAVEAT = (
    "SULFUR BRANCH: three absolute literature anchors, and the model fails all three. "
    "CORRECTED 2026-08-28 (Wave W); this caveat read 'zero absolute literature anchors' until "
    "the full text of 10.1021/jf9705983 was obtained. The anchors are the pH-5.0 aqueous "
    "isotope-dilution rows of Hofmann & Schieberle 1998 Table 1 "
    "(hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5); the model misses them by "
    "12.3x, 29.6x and 14.5x. The old benchmark cys_ribose_140C_Hofmann1998 remains retired -- "
    "both its values are marked no_verifiable_source and the paper confirms they appear nowhere "
    "in it. Its DIRECTION is not thereby worthless -- it predicted an unseen MFT measurement to "
    "1.52x -- but it has the temperature dependence backwards between 100 and 130 C, and it now "
    "also gets glucose-vs-fructose backwards by two orders of magnitude. Sulfur absolutes are "
    "anchored and wrong; sulfur temperature directions are wrong."
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


@dataclass
class CompoundRow:
    """One compound of one arm, with its three aliases collapsed into one record.

    ``Recommender`` deliberately keys ``predicted_ppb`` by canonical SMILES **and** species
    name **and** display label, so that downstream benchmark matching can find a compound
    under whichever spelling a paper used. That is right for matching and wrong for a table:
    reading it naively prints every compound three times. This record is the collapse, built
    from ``result.targets``, which is already one row per compound.
    """

    display: str
    species: Optional[str]
    canon: Optional[str]
    ppb: float
    envelope: Optional[Dict[str, float]] = None
    rate_limiting_family: Optional[str] = None
    route_count: int = 0


@dataclass
class SystemRun:
    """One evaluated arm."""

    spec: Dict[str, Any]
    name: str
    compounds: Dict[str, CompoundRow]
    warnings: List[str] = field(default_factory=list)
    confidence: Dict[str, Any] = field(default_factory=dict)


def _spec_to_bench(spec: Mapping[str, Any]) -> Dict[str, Any]:
    """Reuse the benchmark input contract rather than writing a second one.

    ``benchmark_to_formulation`` already owns the precursor->category routing (which names
    count as sugars, which as lipids, which as matrix cues) and ``benchmark_to_conditions``
    already owns the condition mapping. Re-implementing either here would create a second
    input path that could drift from the one every benchmark is scored through.
    """
    bench: Dict[str, Any] = {
        "benchmark_id": str(spec["name"]),
        "precursors": {
            str(name): {"concentration_mM": float(value)}
            for name, value in spec["precursors"].items()
        },
        "conditions": {
            "temp_C": float(spec["temp_C"]),
            "ph": float(spec["ph"]),
            "water_activity": float(spec["aw"]),
            "time_min": float(spec["time_min"]),
        },
        "protein_type": str(spec.get("protein_type", "free")),
    }
    for optional in ("denaturation_state", "moisture_regime", "sme_kj_per_kg"):
        if spec.get(optional) is not None:
            bench[optional] = spec[optional]
    return bench


def _rate_limiting_family(path: Sequence[Mapping[str, Any]]) -> Optional[str]:
    """The family of the highest-barrier step on a route -- the same rule the Wave S1 test uses."""
    steps = [step for step in path if step.get("barrier") is not None]
    if not steps:
        return None
    return str(max(steps, key=lambda step: float(step["barrier"])).get("family", "unknown"))


def _collapse_targets(result: Any, trace: Mapping[str, Any]) -> Dict[str, CompoundRow]:
    """One row per compound, with its route trace attached."""
    species_names = dict(trace.get("species_names", {}) or {})
    canon_by_species = {str(name): str(canon) for canon, name in species_names.items()}
    debug_paths = dict(trace.get("debug_paths", {}) or {})
    channel_flux = dict(trace.get("debug_channel_flux", {}) or {})

    envelopes = {
        str(compound): {
            "p5": float(envelope.predicted_p5),
            "p50": float(envelope.predicted_p50),
            "p95": float(envelope.predicted_p95),
            "ci_level_pct": int(envelope.ci_level_pct),
            "support_count": int(envelope.support_count),
            "envelope_source": str(envelope.envelope_source),
        }
        for compound, envelope in (result.uncertainty_envelopes or {}).items()
    }

    rows: Dict[str, CompoundRow] = {}
    for target in result.targets or []:
        display = str(target.get("name", "")).strip()
        if not display:
            continue
        ppb = float(target.get("concentration", 0.0) or 0.0)
        species = (target.get("projection") or {}).get("compound")
        species = None if species is None else str(species)
        canon = canon_by_species.get(species) if species else None

        envelope = None
        for alias in (display, species, canon):
            if alias and alias in envelopes:
                envelope = envelopes[alias]
                break

        path = debug_paths.get(canon, []) if canon else []
        rows[display] = CompoundRow(
            display=display,
            species=species,
            canon=canon,
            ppb=ppb,
            envelope=envelope,
            rate_limiting_family=_rate_limiting_family(path or []),
            route_count=len(channel_flux.get(canon, {}) or {}) if canon else 0,
        )
    return rows


def evaluate_system(
    spec: Mapping[str, Any],
    *,
    target_tag: str = DEFAULTS.default_target_tag,
    minimize_tag: str = DEFAULTS.default_minimize_tag,
) -> SystemRun:
    """Run one spec through the shipped forward path. No new science on this line or any other."""
    bench = _spec_to_bench(spec)
    formulation = benchmark_to_formulation(bench)
    formulation["time_minutes"] = float(spec["time_min"])
    conditions = benchmark_to_conditions(bench)

    designer = MaillardPipeline(target_tag=target_tag, minimize_tag=minimize_tag)
    try:
        result = designer.evaluate_single(formulation, conditions)
    except ValueError as exc:
        skipped = getattr(designer, "last_skipped_formulations", []) or []
        reasons = "; ".join(f"{row.get('name')}: {row.get('reason')}" for row in skipped)
        raise SpecError(
            f"{spec['name']}: the formulation could not be evaluated -- {exc}"
            + (f" ({reasons})" if reasons else "")
        ) from exc

    domain_warnings = prepare_cli_confidence(
        result,
        target_tag=target_tag,
        precursor_names=list(spec["precursors"].keys()),
        protein_type=str(formulation.get("protein_type", "free")),
        temp_c=float(spec["temp_C"]),
        ph=float(spec["ph"]),
        aw=float(spec["aw"]),
        formulation=formulation,
        baseline_conditions=conditions,
        designer=designer,
    )

    trace = (getattr(designer, "last_route_traces", {}) or {}).get(str(formulation.get("name")), {})

    return SystemRun(
        spec=dict(spec),
        name=str(spec["name"]),
        compounds=_collapse_targets(result, trace),
        warnings=[str(getattr(w, "description", w)) for w in domain_warnings],
        confidence=dict(result.confidence_metadata or {}),
    )


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


def compare_systems(
    run_a: SystemRun,
    run_b: SystemRun,
    *,
    counts: Optional[Mapping[str, tuple]] = None,
    top_n: Optional[int] = None,
) -> Dict[str, Any]:
    table = dict(counts if counts is not None else load_panel_counts())
    comparison = describe_comparison(run_a.spec, run_b.spec, counts=table)

    compounds = sorted(set(run_a.compounds) | set(run_b.compounds))
    rows: List[Dict[str, Any]] = []
    for compound in compounds:
        row_a = run_a.compounds.get(compound)
        row_b = run_b.compounds.get(compound)
        a_ppb = row_a.ppb if row_a else 0.0
        b_ppb = row_b.ppb if row_b else 0.0
        if a_ppb <= 0.0 and b_ppb <= 0.0:
            continue

        if b_ppb > 0.0 and a_ppb > 0.0:
            ratio: Optional[float] = a_ppb / b_ppb
            ratio_kind = "finite"
        elif a_ppb > 0.0:
            ratio, ratio_kind = None, "b_absent"
        else:
            ratio, ratio_kind = None, "a_absent"

        lo, hi = _ratio_bounds(
            row_a.envelope if row_a else None, row_b.envelope if row_b else None
        )

        axes = list(comparison["axes"]) + compound_axes(compound)
        reliabilities = [reliability_for_axis(axis, table) for axis in axes]
        governing = weakest(reliabilities)

        rows.append(
            {
                "compound": compound,
                "a_ppb": a_ppb,
                "b_ppb": b_ppb,
                "ratio": ratio,
                "ratio_kind": ratio_kind,
                "ratio_lo": lo,
                "ratio_hi": hi,
                "dominant_pathway_a": row_a.rate_limiting_family if row_a else None,
                "dominant_pathway_b": row_b.rate_limiting_family if row_b else None,
                "routes_a": row_a.route_count if row_a else 0,
                "axes": axes,
                "reliability": governing.render() if governing else "no-axis-differs",
                "reliability_verdict": governing.verdict if governing else "n/a",
                "sulfur": is_sulfur_compound(compound),
            }
        )

    rows.sort(key=lambda row: -max(row["a_ppb"], row["b_ppb"]))
    if top_n is not None:
        rows = rows[: max(int(top_n), 0)]

    return {
        "artifact": "maillard_compare",
        "a": {"name": run_a.name, "spec": run_a.spec},
        "b": {"name": run_b.name, "spec": run_b.spec},
        "axes_exercised": comparison["axes"],
        "governing_reliability": (
            comparison["governing"].render() if comparison["governing"] else "no-axis-differs"
        ),
        "per_axis": [
            {
                "axis": item.axis,
                "score": item.score,
                "verdict": item.verdict,
                "note": item.note,
            }
            for item in comparison["per_axis"]
        ],
        "rows": rows,
        "warnings_a": run_a.warnings,
        "warnings_b": run_b.warnings,
        "caveats": {
            "ratio": RATIO_CAVEAT,
            "absolute": ABSOLUTE_CAVEAT,
            "sulfur": SULFUR_CAVEAT if any(row["sulfur"] for row in rows) else None,
        },
    }


# ---------------------------------------------------------------------------------------
# Single-system prediction
# ---------------------------------------------------------------------------------------


def predict_system(
    run: SystemRun, *, counts: Optional[Mapping[str, tuple]] = None, top_n: Optional[int] = None
) -> Dict[str, Any]:
    table = dict(counts if counts is not None else load_panel_counts())
    rows: List[Dict[str, Any]] = []
    for compound, row in sorted(run.compounds.items(), key=lambda item: -item[1].ppb):
        if row.ppb <= 0.0:
            continue
        envelope = row.envelope
        lane = compound_axes(compound)
        lane_reliability = [reliability_for_axis(axis, table) for axis in lane]
        governing = weakest(lane_reliability)
        rows.append(
            {
                "compound": compound,
                "predicted_ppb": row.ppb,
                "range_p5": None if not envelope else envelope["p5"],
                "range_p95": None if not envelope else envelope["p95"],
                "ci_level_pct": None if not envelope else envelope["ci_level_pct"],
                "range_available": envelope is not None,
                "dominant_pathway": row.rate_limiting_family,
                "routes": row.route_count,
                "lane_reliability": governing.render() if governing else None,
                "sulfur": is_sulfur_compound(compound),
            }
        )
    if top_n is not None:
        rows = rows[: max(int(top_n), 0)]

    return {
        "artifact": "maillard_predict",
        "system": {"name": run.name, "spec": run.spec},
        "rows": rows,
        "warnings": run.warnings,
        "confidence_tier": run.confidence.get("tier"),
        "prediction_mode": run.confidence.get("prediction_mode"),
        "decision_mode": run.confidence.get("decision_mode"),
        "caveats": {
            "absolute": ABSOLUTE_CAVEAT,
            "sulfur": SULFUR_CAVEAT if any(row["sulfur"] for row in rows) else None,
            "no_range": (
                "A compound with no range has no counterpart in the Monte-Carlo uncertainty "
                "panel; it is a point prediction with NO measured interval, which is weaker "
                "evidence than a wide interval, not stronger."
            ),
        },
    }


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


def render_compare_text(payload: Mapping[str, Any], *, show_absolute: bool = False) -> str:
    """
    Render a FAST-lane compare payload as text.

    REACHABILITY, stated because it is not obvious and Q1 had to work it out.
    ``show_absolute`` CANNOT BE TRUE FROM THE CLI. Three gates stack: the only
    production caller (``scripts/maillard.py``) passes ``show_absolute=False``
    literally; ``--absolute`` exits 2 on this lane before reaching here; and the
    core lane renders through ``render_compare_core_text`` instead. So the
    absolute block below is dead from the front door.

    It is NOT dead from the library, and that is why it is kept rather than
    deleted: ``tests/unit/test_comparative_cli_2026_08.py`` calls it directly on
    a RAW payload to pin a policy that should outlive the current gating -- if
    absolutes are ever printed, the caveat ships in the SAME block as the
    numbers, not somewhere else on the page.

    THE CONSTRAINT THAT COMES WITH THAT: the absolute block reads ``a_ppb`` and
    ``b_ppb``, which ``screening_payload`` STRIPS. Passing a screened payload
    with ``show_absolute=True`` would raise ``KeyError``. Callers must pass the
    raw payload from ``compare_systems``, as the tests do.
    """
    out: List[str] = []
    a_name, b_name = payload["a"]["name"], payload["b"]["name"]
    banner = f" [{SCREENING_LABEL}]" if payload.get("absolutes_withheld") else ""
    out.append("=" * 96)
    out.append(f"  COMPARE{banner}   A = {a_name}   vs   B = {b_name}")
    out.append("=" * 96)
    out.append("")
    if payload.get("absolutes_withheld"):
        out.append(_wrap(payload["caveats"]["screening"], indent="  !! "))
        out.append("")
    axes = payload["axes_exercised"]
    out.append(
        f"  Axes this comparison moves: {', '.join(axes) if axes else '(none -- the two arms are identical)'}"
    )
    out.append(f"  Governing reliability:      {payload['governing_reliability']}")
    out.append("")
    for item in payload["per_axis"]:
        out.append(f"    - {item['axis']:<18} {item['verdict']:<12} panel {item['score']}")
        if item["note"]:
            out.append(_wrap(item["note"], indent="        "))
    out.append("")

    header = f"  {'compound':<34} {'A/B':>10} {'bound(lo..hi)':>22}  {'reliability':<18} pathway A"
    out.append(header)
    out.append("  " + "-" * 94)
    for row in payload["rows"]:
        if row["ratio"] is None:
            ratio_text = "A only" if row["ratio_kind"] == "b_absent" else "B only"
        else:
            ratio_text = f"{row['ratio']:.3g}x"
        if row["ratio_lo"] is None:
            bound = "no panel envelope"
        else:
            bound = f"{row['ratio_lo']:.3g} .. {row['ratio_hi']:.3g}"
        pathway = row["dominant_pathway_a"] or "-"
        out.append(
            f"  {row['compound'][:33]:<34} {ratio_text:>10} {bound:>22}  "
            f"{row['reliability']:<18} {pathway}"
        )
    out.append("")

    if show_absolute:
        out.append("  ABSOLUTE ppb (requested with --absolute)")
        out.append(f"  {'compound':<34} {'A ppb':>14} {'B ppb':>14}")
        out.append("  " + "-" * 64)
        for row in payload["rows"]:
            out.append(
                f"  {row['compound'][:33]:<34} {_fmt(row['a_ppb']):>14} {_fmt(row['b_ppb']):>14}"
            )
        out.append("")
        out.append("  !! " + "-" * 88)
        out.append(_wrap(payload["caveats"]["absolute"], indent="  !! "))
        out.append("  !! " + "-" * 88)
        out.append("")

    out.append("  HOW TO READ THIS")
    out.append(_wrap(payload["caveats"]["ratio"]))
    if payload["caveats"].get("sulfur"):
        out.append("")
        out.append(_wrap(payload["caveats"]["sulfur"]))
    for label, warnings in (("A", payload["warnings_a"]), ("B", payload["warnings_b"])):
        for warning in warnings:
            out.append(f"  [envelope warning, {label}] {warning}")
    out.append("")
    return "\n".join(out)


def render_predict_text(payload: Mapping[str, Any]) -> str:
    """
    Render a FAST-lane predict payload as text.

    REACHABILITY (Q1, same analysis as ``render_compare_text``): every payload
    that reaches this function FROM THE CLI has been through
    ``screening_payload``, which sets ``absolutes_withheld=True`` unconditionally
    and strips ``predicted_ppb`` / ``range_p5`` / ``range_p95``. So from the
    front door ``withheld`` is ALWAYS true, and the three ``else`` arms below --
    the "range (90% CI), ppb" column heading, the range and point-value bands,
    and the ``no_range`` caveat -- are unreachable, as is ``caveats["absolute"]``.

    They are kept because the library-level tests exercise this renderer on RAW
    payloads, and because those arms are the only thing that would render a FAST
    absolute correctly if the screening policy is ever revisited. Do not read
    them as live behaviour, and do not pass a screened payload to them: the
    fields they read have been deleted by then.
    """
    out: List[str] = []
    banner = f" [{SCREENING_LABEL}]" if payload.get("absolutes_withheld") else ""
    out.append("=" * 96)
    out.append(f"  PREDICT{banner}   {payload['system']['name']}")
    out.append("=" * 96)
    out.append("")
    out.append("  !! " + "-" * 88)
    if payload.get("absolutes_withheld"):
        out.append(f"  !! {SCREENING_LABEL}")
        out.append(_wrap(payload["caveats"]["screening"], indent="  !! "))
    else:
        out.append(_wrap(payload["caveats"]["absolute"], indent="  !! "))
    out.append("  !! " + "-" * 88)
    out.append("")
    out.append(
        f"  confidence tier: {payload.get('confidence_tier')}   "
        f"prediction mode: {payload.get('prediction_mode')}   "
        f"decision mode: {payload.get('decision_mode')}"
    )
    out.append(
        _wrap(
            "That tier is the run-level SCOPE vocabulary -- how close this formulation sits to "
            "systems the model has seen. It is not a validation claim and does not mean the "
            "numbers are right: 0 of 17 benchmarks in the panel are strict-ready, and a "
            "formulation can be 'high' on scope while carrying a 90% interval four decades "
            "wide, as the widths below will usually show."
        )
    )
    out.append("")
    withheld = bool(payload.get("absolutes_withheld"))
    column = "rank" if withheld else "range (90% CI), ppb"
    out.append(f"  {'compound':<34} {column:>28}  {'lane':<18} pathway")
    out.append("  " + "-" * 94)
    for index, row in enumerate(payload["rows"], start=1):
        if withheld:
            # ORDINAL SCREENING: the ordering survives, the magnitudes do not.
            band = f"#{index} (ppb withheld)"
        elif row.get("range_available"):
            band = f"{_fmt(row['range_p5'])} .. {_fmt(row['range_p95'])}"
        else:
            band = f"{_fmt(row.get('predicted_ppb'))} (point, no interval)"
        out.append(
            f"  {row['compound'][:33]:<34} {band:>28}  "
            f"{(row['lane_reliability'] or '-'):<18} {row['dominant_pathway'] or '-'}"
        )
    out.append("")
    if not withheld:
        out.append(_wrap(payload["caveats"]["no_range"]))
    if payload["caveats"].get("sulfur"):
        out.append("")
        out.append(_wrap(payload["caveats"]["sulfur"]))
    for warning in payload["warnings"]:
        out.append(f"  [envelope warning] {warning}")
    out.append("")
    out.append("  For a decision, use `maillard compare` instead: the ranking claims are the")
    out.append("  ones this model was measured to get right (24/35), and a comparison cancels")
    out.append("  the systematic scale error that makes the numbers above unreliable.")
    out.append("")
    return "\n".join(out)


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
            "cached Monte-Carlo panel (results/validation/prediction_uncertainty.json), so it "
            "is only as current as the last trust-loop run."
        )
    )
    out.append("")
    return "\n".join(out)


def to_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, default=str)


# ---------------------------------------------------------------------------------------
# Build Wave B5 -- THE PROPAGATOR CUTOVER
# ---------------------------------------------------------------------------------------
#
# From B5, ``compare`` and ``predict`` route through the KINETIC CORE
# (src/kinetic_core/engine.py). The FAST lane above is NOT deleted -- it is
# demoted to the ordinal-only front end its measured skill supports, and every
# surface it reaches is labelled and stripped of absolutes.
#
# The two rules this section enforces, mechanically rather than by convention:
#   1. NO ABSOLUTE ppb FROM THE FAST LANE REACHES A USER-FACING SURFACE.
#      ``screening_payload`` removes the fields; the CLI applies it to every
#      FAST payload before rendering or serialising, so there is no path from a
#      FAST ppb to a terminal.
#   2. Every FAST surface carries the ORDINAL SCREENING label.

SCREENING_LABEL = "ORDINAL SCREENING"

SCREENING_CAVEAT = (
    "ORDINAL SCREENING LANE. These are RANKINGS, not concentrations. The FAST lane's "
    "absolute ppb are withheld from every user-facing surface as of Wave B5, because its "
    "measured absolute skill does not support them (median 6.0x on the free-precursor "
    "hold-out, 67-94x on the matrix lane, 1 of 5 genuine extrapolation rows inside the 90% "
    "CI). What the FAST lane is measured to do is ORDER things: the directional panel scores "
    "24/35 on strictly independent claims. Use it to sort candidates; use `--lane core` for "
    "any quantity."
)

# CORRECTED 2026-08-29 (Wave Q1). This caveat is printed to EVERY core-lane
# user, and three of its factual claims had been falsified by the waves that
# followed it: B6 added the lipid lane (so "no lipid-oxidation path" was wrong),
# B7 added HMF and DMHF as trunk targets (so "no HMF and no DMHF" was wrong),
# and B6's direct-sum co-integration means the lipid lane DOES compose with any
# one Maillard lane (so "its three lanes do not compose" was wrong, and there
# are four lanes). A caveat that overstates the model's limits is not the safe
# direction to be wrong in: it teaches users to distrust answers the model can
# actually support, and it goes stale invisibly because nothing tests prose.
CORE_CAVEAT = (
    "KINETIC CORE. Absolute concentrations come from the mass-action network (frozen "
    "B1/B2.x/B3/B6/B7 parameters), and they are reported WITH their envelope declaration. The "
    "core refuses what it cannot name -- ask it for 1-hexanol, 2-pentylfuran or propanal and it "
    "will tell you why it will not answer, rather than answering. Its four lanes (trunk, sulfur, "
    "acrylamide, lipid) do not compose freely: the lipid lane co-integrates with ONE Maillard "
    "lane as a direct sum, and the Maillard lanes do not compose with each other. A refusal is "
    "an output, not a failure. Read the cutover final exam "
    "(results/validation/cutover_final_exam.md) for its measured out-of-sample accuracy before "
    "trusting any number here."
)

#: Fields removed from every FAST payload before it reaches a user.
#: Q1 added ``ci_level_pct`` and ``range_available``. Both are properties of an
#: absolute interval that is itself being withheld, so leaving them in described
#: a quantity the payload no longer carried -- ``range_available: true`` next to
#: no range. Nothing rendered them, so this strips two fields that were dead
#: weight rather than changing any output.
_FAST_ABSOLUTE_FIELDS = (
    "a_ppb", "b_ppb", "predicted_ppb", "range_p5", "range_p95",
    "ci_level_pct", "range_available",
)


def screening_payload(payload: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Strip FAST-lane absolutes and stamp the ORDINAL SCREENING label.

    Applied at the CLI boundary to every FAST payload. The ratio, the ordering
    and the reliability columns survive untouched -- those are the lane's
    measured product. The ppb columns do not.
    """
    out = dict(payload)
    out["lane"] = "fast"
    out["lane_label"] = SCREENING_LABEL
    rows = []
    for row in out.get("rows", []):
        clean = {k: v for k, v in dict(row).items() if k not in _FAST_ABSOLUTE_FIELDS}
        clean["absolute_ppb_withheld"] = True
        clean["lane_label"] = SCREENING_LABEL
        rows.append(clean)
    out["rows"] = rows
    caveats = dict(out.get("caveats") or {})
    caveats["screening"] = SCREENING_CAVEAT
    # The old absolute caveat described numbers this lane no longer emits.
    caveats["absolute"] = (
        "NOT APPLICABLE on the screening lane: absolute ppb are withheld here. "
        "Use `--lane core`."
    )
    out["caveats"] = caveats
    out["absolutes_withheld"] = True
    return out


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
        "caveats": {"core": CORE_CAVEAT, "ratio": RATIO_CAVEAT},
        "engine": engine_metadata(),
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
            f"  ...of which {undefined} are UNDEFINED (one arm at exactly zero) and "
            f"resolve nothing: see the 'resolved' column, not this count."
        )
    out.append("")
    out.append(f"  {'compound':<38} {'A/B':>12} {'direction':<14} resolved")
    out.append("  " + "-" * 88)
    for row in rows:
        ratio = row.get("ratio_a_over_b")
        if not isinstance(ratio, (int, float)):
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
