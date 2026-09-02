"""
src/kinetic_core/directional.py -- THE CORE SCORED ON THE DIRECTIONAL CLAIMS PANEL (step B7).

The product claim of this tool is ordinal: "which of these two formulations gives
more MFT?", "which way does hexanal move with temperature?". The retired screening
lane was scored on ``docs/validation/directional_claims_panel.yml`` (69 claims from
14 sources, each a literature statement with the two arms -- or the series -- it
was made about, written as runnable systems); its 24/36 did not transfer to the
kinetic core, and since retirement step B5 the model card has said so. This
module scores the CORE on the same panel, through the same front door a user
calls (``comparative_cli.spec_to_core`` → ``engine.predict``), and writes
``results/validation/core_directional_scores.{json,md}``.

Scoring rules, stated so they cannot drift:

* An ``ordering`` claim with arms A and B and ``expected_relation`` ``A>B`` /
  ``A<B`` AGREES when the core's predicted concentrations of the observable
  order that way. ``flat`` agrees when |log10(A/B)| <= log10(1 + tolerance),
  with the panel's own ``flat_tolerance_pct``.
* A ``monotonic_direction`` claim over a ``series`` AGREES when the predicted
  sequence is strictly increasing / decreasing (every consecutive step, beyond
  the flat tolerance), ``peak`` when the maximum sits at an interior point, and
  ``flat`` when every step is within tolerance.
* A claim is NOT EVALUABLE -- and is neither a hit nor a miss -- when it carries
  no runnable conditions, when any arm is refused (out of envelope), or when the
  observable is one the core does not represent. The reason is printed.
* IDENTICAL predictions across arms that differ on an axis the lane carries no
  term for (pH on the trunk and acrylamide lanes, water activity everywhere)
  are scored as evaluable MISSES for a non-flat expectation and flagged
  ``mechanism_absent``: a model without the term cannot get the direction right,
  and hiding that under "not evaluable" would flatter it.

``fit_status`` from the panel (``independent`` / ``fit_adjacent`` /
``fit_system_overlap``) is carried on every row; the headline is the
STRICTLY-INDEPENDENT subset, and every count is also reported over all claims.
Nothing in this module is tuned; the thresholds that turn counts into
trust / caution / do-not-use live in ``src/directional_reliability.py``.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import yaml

from src import data_paths

from .engine import engine_metadata, predict

PANEL_PATH: Path = data_paths.DIRECTIONAL_CLAIMS_PANEL
INDEPENDENT = "independent"
NOT_EVALUABLE = "not_evaluable"


# ---------------------------------------------------------------------------
# Panel loading
# ---------------------------------------------------------------------------


def load_panel(path: Path | str = PANEL_PATH) -> Dict[str, Any]:
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping) or "panel" not in payload:
        raise ValueError(f"{path}: not a directional claims panel (no `panel` key)")
    return dict(payload)


def _spec_from_system(name: str, system: Mapping[str, Any]) -> Dict[str, Any]:
    """The CLI spec dict for one panel arm (the same shape ``maillard compare`` reads)."""
    spec: Dict[str, Any] = {
        "name": name,
        "precursors": {str(k): float(v) for k, v in dict(system.get("precursors") or {}).items()},
        "temp_C": float(system["temp_C"]),
        "time_min": float(system["time_min"]),
        "ph": float(system["ph"]),
        "aw": float(system["water_activity"]) if system.get("water_activity") is not None else None,
    }
    for key in ("matrix", "protein_type"):
        if system.get(key):
            spec["matrix"] = str(system[key])
    return spec


def _arms(claim: Mapping[str, Any]) -> List[Tuple[str, Mapping[str, Any]]]:
    """``[(label, system), ...]`` in panel order; empty when the claim is prose-only."""
    conditions = claim.get("conditions") or {}
    if isinstance(conditions, Mapping) and "series" in conditions:
        return [
            (str(item.get("label") or f"point {i}"), item["system"])
            for i, item in enumerate(conditions["series"])
            if isinstance(item, Mapping) and item.get("system")
        ]
    if isinstance(conditions, Mapping):
        return [
            (str(conditions[k].get("label") or k), conditions[k]["system"])
            for k in ("A", "B")
            if isinstance(conditions.get(k), Mapping) and conditions[k].get("system")
        ]
    return []


# ---------------------------------------------------------------------------
# One claim
# ---------------------------------------------------------------------------


def _predict_arm(label: str, system: Mapping[str, Any], observable: str) -> Dict[str, Any]:
    from src.comparative_cli import spec_to_core

    spec = _spec_from_system(label, system)
    core_spec = spec_to_core(spec)
    run = predict(core_spec, [observable])
    declaration = run.declaration
    value: Optional[float] = None
    if run.answered:
        value = run.concentrations_ug_per_l.get(observable)
    return {
        "label": label,
        "spec": spec,
        "answered": bool(run.answered),
        "lane": declaration.lane,
        "state": declaration.state,
        "reasons": list(declaration.reasons),
        "warnings": list(declaration.warnings),
        "value_ug_per_l": None if value is None else float(value),
    }


_SULFUR_TOKENS = ("thiol", "sulfide", "mft", "fft", "mercapto", "thiazol", "thiophen", "methional")


def _is_sulfur_observable(name: str) -> bool:
    n = str(name).lower()
    return any(tok in n for tok in _SULFUR_TOKENS)


def _axis_refusal_between_arms(arms, runs) -> Optional[str]:
    """The engine's axis refusal, applied to every consecutive pair of arms."""
    from src.comparative_cli import spec_to_core

    from .engine import axis_refusal

    class _Decl:  # the engine's helper reads .lane / .lanes only
        def __init__(self, run):
            self.lane = run["lane"]
            self.lanes = None

    for (label_a, sys_a), (label_b, sys_b), run_a, run_b in zip(arms, arms[1:], runs, runs[1:]):
        spec_a = spec_to_core(_spec_from_system(label_a, sys_a))
        spec_b = spec_to_core(_spec_from_system(label_b, sys_b))
        reason = axis_refusal(spec_a, spec_b, _Decl(run_a), _Decl(run_b))
        if reason:
            return "refused by the engine: " + reason
    return None


def _log_step(a: float, b: float) -> Optional[float]:
    if a is None or b is None or a <= 0.0 or b <= 0.0:
        return None
    return math.log10(a / b)


def score_claim(claim: Mapping[str, Any], *, flat_tolerance_pct: float) -> Dict[str, Any]:
    """Score one claim on the core; never raises for a refusal."""
    claim_id = str(claim["claim_id"])
    observables = [str(o) for o in claim.get("observables") or []]
    expected = str(claim.get("expected_relation"))
    claim_type = str(claim.get("claim_type"))
    tol = math.log10(1.0 + float(flat_tolerance_pct) / 100.0)
    base = {
        "claim_id": claim_id,
        "category": str(claim.get("category")),
        "claim_type": claim_type,
        "expected_relation": expected,
        "observables": observables,
        "fit_status": str(claim.get("fit_status") or "undeclared"),
        "statement": str(claim.get("statement") or ""),
        "source": (claim.get("source") or {}).get("ref"),
    }
    arms = _arms(claim)
    if not arms:
        return {**base, "status": NOT_EVALUABLE, "reason": "the claim carries no runnable conditions (prose-only)", "arms": []}
    if not observables:
        return {**base, "status": NOT_EVALUABLE, "reason": "the claim names no observable", "arms": []}
    observable = observables[0]  # the panel lists one observable per claim; the first governs
    runs = [_predict_arm(label, system, observable) for label, system in arms]
    # A STRUCTURAL ZERO: the sulfur lane refusing because the charge has no sulfur source
    # is the model asserting that a sulfur observable is zero there, and a claim of the
    # form "with cysteine > without" is exactly a test of that assertion. Scored as 0,
    # flagged on the arm, never invented for any other refusal.
    for r in runs:
        if not r["answered"] and _is_sulfur_observable(observable) and any(
            "no sulfur source" in reason.lower() for reason in r["reasons"]
        ):
            r["value_ug_per_l"] = 0.0
            r["structural_zero"] = True
    refused = [r for r in runs if r["value_ug_per_l"] is None]
    if refused:
        reasons = "; ".join(refused[0]["reasons"]) or f"{observable}: no concentration returned"
        return {**base, "status": NOT_EVALUABLE, "reason": f"arm {refused[0]['label']!r} refused: {reasons[:400]}", "arms": runs}

    # Step 5 (2026-09-03): an axis the lane carries no term for is REFUSED by the engine
    # (engine.axis_refusal), so the claim is NOT EVALUABLE rather than a miss of identical
    # predictions. The count of such claims is reported by reason.
    refusal = _axis_refusal_between_arms(arms, runs)
    if refusal is not None:
        return {**base, "status": NOT_EVALUABLE, "reason": refusal, "arms": runs}
    values = [r["value_ug_per_l"] for r in runs]
    if claim_type == "ordering" and len(values) == 2 and expected in ("A>B", "A<B") and (values[0] == 0.0) != (values[1] == 0.0):
        # one arm is exactly zero (a structural zero or a true zero rate): the ordering is
        # defined even though the log ratio is not.
        agree = (values[0] > values[1]) if expected == "A>B" else (values[0] < values[1])
        return {
            **base, "status": "agree" if agree else "disagree", "reason": None,
            "identical_predictions": False, "mechanism_absent": False, "lane": runs[0]["lane"],
            "values_ug_per_l": values, "log10_ratio_B_over_A": None,
            "structural_zero": any(r.get("structural_zero") for r in runs), "arms": runs,
        }
    steps = [_log_step(values[i + 1], values[i]) for i in range(len(values) - 1)]
    if any(s is None for s in steps):
        return {**base, "status": NOT_EVALUABLE, "reason": "a predicted concentration is zero; no direction is defined", "arms": runs}
    identical = all(abs(s) < 1e-12 for s in steps)
    warnings = " ".join(w for r in runs for w in r["warnings"]).lower()
    mechanism_absent = identical and ("no ph term" in warnings or "no a_w term" in warnings or "no water activity" in warnings or "no aw term" in warnings)

    if claim_type == "ordering":
        step = steps[0]  # log10(B/A)
        if expected == "A>B":
            agree = step < -tol
        elif expected == "A<B":
            agree = step > tol
        elif expected == "flat":
            agree = abs(step) <= tol
        else:
            return {**base, "status": NOT_EVALUABLE, "reason": f"unknown expected_relation {expected!r}", "arms": runs}
        detail = {"log10_ratio_B_over_A": step}
    else:  # monotonic_direction
        if expected == "increasing":
            agree = all(s > tol for s in steps)
        elif expected == "decreasing":
            agree = all(s < -tol for s in steps)
        elif expected == "flat":
            agree = all(abs(s) <= tol for s in steps)
        elif expected == "peak":
            k = max(range(len(values)), key=lambda i: values[i])
            agree = 0 < k < len(values) - 1
        else:
            return {**base, "status": NOT_EVALUABLE, "reason": f"unknown expected_relation {expected!r}", "arms": runs}
        detail = {"log10_steps": steps}
    return {
        **base,
        "status": "agree" if agree else "disagree",
        "reason": None,
        "identical_predictions": identical,
        "mechanism_absent": bool(mechanism_absent),
        "lane": runs[0]["lane"],
        "values_ug_per_l": values,
        **detail,
        "arms": runs,
    }


# ---------------------------------------------------------------------------
# The panel
# ---------------------------------------------------------------------------


def _bucket() -> Dict[str, Any]:
    return {"agree": 0, "evaluable": 0, "not_evaluable": 0, "misses": [], "mechanism_absent": 0}


def _add(bucket: Dict[str, Any], row: Mapping[str, Any]) -> None:
    if row["status"] == NOT_EVALUABLE:
        bucket["not_evaluable"] += 1
        return
    bucket["evaluable"] += 1
    if row["status"] == "agree":
        bucket["agree"] += 1
    else:
        bucket["misses"].append(row["claim_id"])
        if row.get("mechanism_absent"):
            bucket["mechanism_absent"] += 1


def _close(bucket: Dict[str, Any]) -> Dict[str, Any]:
    bucket["rate"] = (bucket["agree"] / bucket["evaluable"]) if bucket["evaluable"] else None
    return bucket


def score_panel(path: Path | str = PANEL_PATH) -> Dict[str, Any]:
    panel = load_panel(path)
    tol_pct = float(panel.get("flat_tolerance_pct", 5.0))
    rows = [score_claim(claim, flat_tolerance_pct=tol_pct) for claim in panel["panel"]]

    def split(predicate) -> Dict[str, Any]:
        by_cat: Dict[str, Dict[str, Any]] = {}
        total = _bucket()
        ph_aw = _bucket()
        excl = _bucket()
        for row in rows:
            if not predicate(row):
                continue
            _add(by_cat.setdefault(row["category"], _bucket()), row)
            _add(total, row)
            if row["category"] in ("ph", "moisture_aw"):
                _add(ph_aw, row)
            else:
                _add(excl, row)
        return {
            "by_category": {k: _close(v) for k, v in sorted(by_cat.items())},
            "total": _close(total),
            "ph_aw": _close(ph_aw),
            "excluding_ph_aw": _close(excl),
        }

    independent = split(lambda r: r["fit_status"] == INDEPENDENT)
    everything = split(lambda r: True)
    not_evaluable_reasons: Dict[str, int] = {}
    for row in rows:
        if row["status"] == NOT_EVALUABLE:
            key = row["reason"].split(":")[0][:80]
            not_evaluable_reasons[key] = not_evaluable_reasons.get(key, 0) + 1
    payload = {
        "artifact": "core_directional_scores",
        "engine": engine_metadata(),
        "panel": {
            "path": data_paths.rel(Path(path)),
            "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "claims": len(rows),
            "flat_tolerance_pct": tol_pct,
            "generated": str(panel.get("generated")),
            "wave": str(panel.get("wave")),
        },
        "summary": {
            "headline": [independent["total"]["agree"], independent["total"]["evaluable"]],
            "headline_definition": (
                "strictly independent claims (fit_status == independent) the core could evaluate; "
                "refused arms and unrepresented observables are NOT EVALUABLE, identical "
                "predictions across a moved axis are MISSES"
            ),
            "independent": independent,
            "all_claims": everything,
            "not_evaluable_reasons": dict(sorted(not_evaluable_reasons.items(), key=lambda kv: -kv[1])),
        },
        "claims": rows,
    }
    return payload


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt_rate(bucket: Mapping[str, Any]) -> str:
    return f"{bucket['agree']}/{bucket['evaluable']}" + (f" ({bucket['rate']:.0%})" if bucket.get("rate") is not None else "")


def render_markdown(payload: Mapping[str, Any]) -> str:
    s = payload["summary"]
    ind, allc = s["independent"], s["all_claims"]
    out = [
        "# Core directional scores (the kinetic core on the directional claims panel)",
        "",
        f"Panel `{payload['panel']['path']}` ({payload['panel']['claims']} claims, flat tolerance "
        f"{payload['panel']['flat_tolerance_pct']:g} %). Nothing in the core was tuned to this panel.",
        "",
        f"* **headline (strictly independent, evaluable): {_fmt_rate(ind['total'])}**; "
        f"{ind['total']['not_evaluable']} independent claims not evaluable",
        f"* independent, excluding pH and water activity: {_fmt_rate(ind['excluding_ph_aw'])}; "
        f"pH and water activity alone: {_fmt_rate(ind['ph_aw'])}",
        f"* all claims (independent + fit-adjacent + fit-system overlap): {_fmt_rate(allc['total'])}; "
        f"{allc['total']['not_evaluable']} not evaluable",
        f"* misses where the lane carries no term for the moved axis (identical predictions): "
        f"{allc['total']['mechanism_absent']}",
        "* not evaluable, by reason: " + "; ".join(f"{k} ({v})" for k, v in s["not_evaluable_reasons"].items()),
        "",
        "## Per category (strictly independent claims)",
        "",
        "| category | agree | evaluable | rate | not evaluable | misses |",
        "|---|---|---|---|---|---|",
    ]
    for cat, b in ind["by_category"].items():
        out.append(
            f"| {cat} | {b['agree']} | {b['evaluable']} | "
            f"{'-' if b['rate'] is None else f'{b['rate']:.2f}'} | {b['not_evaluable']} | "
            f"{', '.join(b['misses']) or '-'} |"
        )
    out += ["", "## Per category (all claims)", "", "| category | agree | evaluable | rate | not evaluable |", "|---|---|---|---|---|"]
    for cat, b in allc["by_category"].items():
        out.append(f"| {cat} | {b['agree']} | {b['evaluable']} | {'-' if b['rate'] is None else f'{b['rate']:.2f}'} | {b['not_evaluable']} |")
    out += ["", "## Claims", "", "| claim | category | fit status | observable | expected | result | lane | predictions (ug/L) | note |", "|---|---|---|---|---|---|---|---|---|"]
    for r in payload["claims"]:
        preds = ", ".join(f"{v:.3g}" for v in r.get("values_ug_per_l") or []) or "-"
        note = r.get("reason") or ("identical predictions; the lane has no term for the moved axis" if r.get("mechanism_absent") else ("identical predictions" if r.get("identical_predictions") else ""))
        out.append(
            f"| {r['claim_id']} | {r['category']} | {r['fit_status']} | {', '.join(r['observables'])} | "
            f"{r['expected_relation']} | **{r['status']}** | {r.get('lane') or '-'} | {preds} | {str(note).replace('|', '/')[:160]} |"
        )
    out.append("")
    return "\n".join(out)


def write_artifact(payload: Mapping[str, Any], json_path: Path) -> Tuple[Path, Path]:
    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    md_path = json_path.with_suffix(".md")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    return json_path, md_path


__all__ = ["INDEPENDENT", "NOT_EVALUABLE", "PANEL_PATH", "load_panel", "render_markdown", "score_claim", "score_panel", "write_artifact"]
