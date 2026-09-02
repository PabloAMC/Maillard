"""
src/kinetic_core/scoring.py -- THE CORE'S OWN PANEL SCORECARD (retirement step B3).

Until B3 the kinetic core was scored by one script only, the pre-registered
cutover exam (``scripts/generators/generate_cutover_final_exam.py``), against a
single 3x band and always PAIRED with the legacy lane. Every headline the
repository publishes (panel size, evidence-role split, predictive passes,
strict-ready count, the hold-out's genuine-extrapolation hits) came from the
legacy harness ``src/benchmark_validation.py``, which is being retired.

This module scores the core on the SAME union panel the envelope uses
(``panel.panel_bundles``: trust loop + maillard_path hold-outs + external matrix
bundles), through the same bundle -> spec mapping, and reports in the legacy
harness's vocabulary so the B4 re-pin can put the two side by side:

* per (benchmark, compound): the point, the fold error, |log10| error, whether
  it lands inside the 3x pass band AND inside the bundle's OWN declared scale
  contract (``validation_contract.scale_thresholds``, falling back to the
  global ``BenchmarkThresholds`` exactly as the legacy resolver did), the core's
  B4 reliability interval and whether the measurement falls inside it;
* per benchmark: coverage, max ratio, mean |log10|, ``scale_status``,
  ``ranking_status`` (Pearson of log10 values over >= 3 answered ppb rows),
  ``overall_status`` in {pass, pass-no-ranking, scale-gap, ranking-gap,
  coverage-gap, refused}, ``strict_ready`` under the same tier / execution-path
  eligibility the legacy gate applied, ``blocking_issues``;
* evidence: the CORE's evidence role (``fit_targets.core_evidence_role``) next
  to the legacy one, the fit records of both worlds, and an ``in_core_fit`` flag
  on every row the sulfur fit read -- including the xylose pH-5 hold-out row.

Two coverage figures are printed, deliberately: ``honest_literature`` follows
the repository's leverage rule (per-row-recovery rows out, low-leverage fit
rows in, annotated) and ``out_of_sample`` additionally removes every row a core
fit saw. The second is the smaller number and the stronger claim.
"""
from __future__ import annotations

import hashlib
import math
import statistics
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from src import data_paths
from src.benchmark_metadata import (
    benchmark_signal_origin,
    get_benchmark_metadata,
    resolve_scale_thresholds,
)
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT

from .engine import engine_metadata, predict
from .fit_targets import core_evidence_role, core_fit_targets, fit_target_of, in_core_fit
from .panel import (
    RATIO_UNIT_FACTORS,
    SHARED_WITH_HOLDOUT_PANEL,
    bundle_targets,
    core_native_value,
    core_spec,
    limiting_precursor_molar,
    load_bundle,
    measured_value,
    panel_bundles,
    quantification_family,
)

#: The hold-out's pre-registered pass band (cutover exam, B2.x prereg): a point
#: is "within band" when predicted / measured lies in [1/3, 3].
PASS_BAND_LEVEL = 3.0

OVERALL_STATUSES = (
    "pass", "pass-no-ranking", "scale-gap", "ranking-gap", "coverage-gap", "refused",
)


# ---------------------------------------------------------------------------
# Row arithmetic
# ---------------------------------------------------------------------------


def fold_error(predicted: Optional[float], measured: Optional[float]) -> Optional[float]:
    """max(p/m, m/p) for two positive finite numbers, else None (same as the exam)."""
    if predicted is None or measured is None:
        return None
    if not (predicted > 0.0 and measured > 0.0):
        return None
    ratio = predicted / measured
    if not math.isfinite(ratio) or ratio <= 0.0:
        return None
    return max(ratio, 1.0 / ratio)


def summarise_folds(folds: Sequence[Optional[float]]) -> Dict[str, Any]:
    """n, median / worst / best fold, median |log10|, geometric-mean fold (the exam's summary)."""
    clean = [f for f in folds if f is not None and math.isfinite(f)]
    if not clean:
        return {
            "n": 0, "median_fold_error": None, "worst_fold_error": None,
            "best_fold_error": None, "median_abs_log10_error": None, "geometric_mean_fold": None,
        }
    logs = [abs(math.log10(f)) for f in clean]
    return {
        "n": len(clean),
        "median_fold_error": statistics.median(clean),
        "worst_fold_error": max(clean),
        "best_fold_error": min(clean),
        "median_abs_log10_error": statistics.median(logs),
        "geometric_mean_fold": 10.0 ** (sum(logs) / len(logs)),
    }


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy)


# ---------------------------------------------------------------------------
# One benchmark
# ---------------------------------------------------------------------------


def score_benchmark(
    path: Path | str, panel_tag: str, *, pass_band: float = PASS_BAND_LEVEL
) -> Dict[str, Any]:
    """Score every target of one bundle on the core; never raises for a refusal."""
    path = Path(path)
    bench = load_bundle(path)
    benchmark_id = str(bench.get("benchmark_id") or path.stem)
    metadata = get_benchmark_metadata(bench)
    protein_type = str(bench.get("protein_type", "free") or "free")
    thresholds = DEFAULT_VALIDATION_CONTRACT.thresholds
    contract = resolve_scale_thresholds(bench, protein_type=protein_type, thresholds=thresholds)
    declared = (bench.get("validation_contract") or {}).get("scale_thresholds") or {}
    contract_source = "bundle validation_contract" if declared else "global BenchmarkThresholds default"
    conditions = bench.get("conditions") or {}
    family, family_source = quantification_family(bench)
    spec = core_spec(bench, use_buffer=True)
    _, limiting_molar = limiting_precursor_molar(bench)

    rows: List[Dict[str, Any]] = []
    refused: List[Dict[str, Any]] = []
    for compound, target_spec in bundle_targets(bench).items():
        unit = str(target_spec.get("target_unit", "ppb"))
        measured = measured_value(bench, compound, target_spec)
        run = predict(spec, [compound])
        declaration = run.declaration
        reason = None
        if not run.answered:
            reason = " | ".join(declaration.reasons) or "out of envelope"
        elif unit != "ppb" and unit not in RATIO_UNIT_FACTORS:
            reason = f"target unit {unit!r} has no core conversion"
        elif measured is None:
            reason = "the bundle carries no measured value for this compound"
        if reason is not None:
            refused.append(
                {
                    "benchmark_id": benchmark_id, "panel": panel_tag, "compound": compound,
                    "target_unit": unit, "reason": reason,
                    "envelope_state": declaration.state, "lane": declaration.lane,
                }
            )
            continue
        predicted = core_native_value(run, compound, unit, limiting_molar)
        fold = fold_error(predicted, measured)
        interval = None
        within_interval = None
        if unit == "ppb":
            try:
                band = run.absolutes().get(compound)
            except Exception:  # noqa: BLE001 -- the interval is a bonus, never a refusal
                band = None
            if band is not None:
                interval = [float(band.lo_ug_per_l), float(band.hi_ug_per_l)]
                within_interval = bool(band.lo_ug_per_l <= float(measured) <= band.hi_ug_per_l)
        rows.append(
            {
                "compound": compound,
                "target_unit": unit,
                "measured": float(measured),
                "measured_ppb": float(measured) if unit == "ppb" else None,
                "predicted": None if predicted is None else float(predicted),
                "fold_error": fold,
                "abs_log10_error": None if fold is None else math.log10(fold),
                "within_band": None if fold is None else bool(fold <= pass_band),
                "within_contract_ratio": None if fold is None else bool(fold <= contract["max_ratio"]),
                "interval_ug_per_L": interval,
                "measured_within_interval": within_interval,
                "lane": declaration.lane,
                "envelope_state": declaration.state,
                "declaration_warnings": list(declaration.warnings),
                "shared_with": SHARED_WITH_HOLDOUT_PANEL.get((benchmark_id, compound)),
                "in_core_fit": in_core_fit(benchmark_id, compound),
            }
        )

    total = len(rows) + len(refused)
    coverage = (len(rows) / total) if total else 0.0
    folds = [r["fold_error"] for r in rows]
    clean = [f for f in folds if f is not None]
    max_ratio = max(clean) if clean else None
    mean_abs_log10 = statistics.fmean(math.log10(f) for f in clean) if clean else None

    ppb_pairs = [
        (math.log10(r["predicted"]), math.log10(r["measured"]))
        for r in rows
        if r["target_unit"] == "ppb" and r["predicted"] and r["measured"]
        and r["predicted"] > 0 and r["measured"] > 0
    ]
    pearson = (
        _pearson([p for p, _ in ppb_pairs], [m for _, m in ppb_pairs])
        if len(ppb_pairs) >= thresholds.min_matched_for_ranking else None
    )
    if len(ppb_pairs) >= thresholds.min_matched_for_ranking and pearson is not None:
        ranking_status = "pass" if pearson >= thresholds.ranking_threshold else "fail"
    elif rows:
        ranking_status = "n/a"
    else:
        ranking_status = "fail"

    if max_ratio is None:
        scale_status = "fail"
    elif max_ratio <= contract["max_ratio"] and (
        mean_abs_log10 is None or mean_abs_log10 <= contract["mean_abs_log10_error"]
    ):
        scale_status = "pass"
    else:
        scale_status = "fail"

    if not rows:
        overall = "refused"
    elif coverage < thresholds.full_coverage_threshold:
        overall = "coverage-gap"
    elif ranking_status == "fail":
        overall = "ranking-gap"
    elif scale_status == "fail":
        overall = "scale-gap"
    elif ranking_status == "n/a":
        overall = "pass-no-ranking"
    else:
        overall = "pass"

    blocking: List[str] = []
    if not rows:
        blocking.append("every target refused: " + "; ".join(r["reason"] for r in refused))
    if rows and coverage < thresholds.full_coverage_threshold:
        blocking.append(f"coverage {coverage:.1%} < {thresholds.full_coverage_threshold:.0%}")
    if ranking_status == "fail":
        blocking.append(f"ranking {'n/a' if pearson is None else f'{pearson:.3f}'} < {thresholds.ranking_threshold:.2f}")
    if rows and scale_status == "fail":
        blocking.append(f"max ratio {'n/a' if max_ratio is None else f'{max_ratio:.3f}'} > {contract['max_ratio']:.2f}")
        if mean_abs_log10 is not None and mean_abs_log10 > contract["mean_abs_log10_error"]:
            blocking.append(
                f"mean |log10 ratio| {mean_abs_log10:.3f} > {contract['mean_abs_log10_error']:.3f}"
            )
    eligible = DEFAULT_VALIDATION_CONTRACT.is_strict_gate_eligible(
        tier=metadata.tier, execution_path=metadata.execution_path
    )
    strict_ready = bool(
        rows and coverage >= thresholds.full_coverage_threshold
        and ranking_status != "fail" and scale_status == "pass" and eligible
    )
    if rows and not eligible and not blocking:
        blocking.append(
            f"tier {metadata.tier} / {metadata.execution_path} is outside the strict gate"
        )

    return {
        "benchmark_id": benchmark_id,
        "bench_file": data_paths.rel(path),
        "panel": panel_tag,
        "tier": metadata.tier,
        "family": metadata.family,
        "execution_path": metadata.execution_path,
        "protein_type": protein_type,
        "conditions": {
            k: conditions.get(k) for k in ("temp_C", "time_min", "ph", "water_activity")
        },
        "signal_origin": benchmark_signal_origin(path),
        "evidence_role": core_evidence_role(benchmark_id, path),
        "fit_target_of": fit_target_of(benchmark_id),
        "core_fit_targets": [t.as_dict() for t in core_fit_targets(benchmark_id)],
        "quantification_family": family,
        "quantification_source": family_source,
        "scale_thresholds": {**contract, "source": contract_source},
        "compounds": rows,
        "refused_compounds": refused,
        "matched_compounds": len(rows),
        "total_compounds": total,
        "coverage": coverage,
        "max_ratio": max_ratio,
        "mean_abs_log10_error": mean_abs_log10,
        "pearson_r_log10": pearson,
        "ranking_status": ranking_status,
        "scale_status": scale_status,
        "overall_status": overall,
        "strict_gate_eligible": eligible,
        "strict_ready": strict_ready,
        "blocking_issues": blocking,
        "within_band_count": sum(1 for r in rows if r["within_band"]),
        "lanes": sorted({str(r["lane"]) for r in rows if r["lane"]}),
    }


# ---------------------------------------------------------------------------
# The panel
# ---------------------------------------------------------------------------


def _bucket() -> Dict[str, Any]:
    return {"benchmarks": 0, "rows": 0, "within_band": 0, "contract_passes": 0,
            "strict_ready": 0, "folds": []}


def _close(bucket: Dict[str, Any]) -> Dict[str, Any]:
    folds = bucket.pop("folds")
    bucket["within_band_rate"] = (bucket["within_band"] / bucket["rows"]) if bucket["rows"] else None
    bucket["fold_summary"] = summarise_folds(folds)
    return bucket


def _report_sha(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score_panel(
    bundles: Optional[Sequence[Tuple[Path, str]]] = None, *, pass_band: float = PASS_BAND_LEVEL
) -> Dict[str, Any]:
    """Score the union panel (or ``bundles`` = ``[(path, panel_tag), ...]``) on the core."""
    if bundles is None:
        bundles, skipped = panel_bundles()
    else:
        bundles, skipped = list(bundles), []

    benches = [score_benchmark(path, tag, pass_band=pass_band) for path, tag in bundles]

    by_panel: Dict[str, Dict[str, Any]] = {}
    by_role: Dict[str, Dict[str, Any]] = {}
    by_lane: Dict[str, Dict[str, Any]] = {}
    honest = _bucket()
    out_of_sample = _bucket()
    in_fit = _bucket()
    all_rows = 0
    all_within = 0
    for b in benches:
        pb = by_panel.setdefault(b["panel"], _bucket())
        rb = by_role.setdefault(b["evidence_role"], _bucket())
        for bucket in (pb, rb):
            bucket["benchmarks"] += 1
            bucket["contract_passes"] += int(b["overall_status"] in {"pass", "pass-no-ranking"})
            bucket["strict_ready"] += int(b["strict_ready"])
        literature = b["signal_origin"] == "external_literature" and b["evidence_role"] != "fit_recovery"
        seen_h = seen_o = seen_f = False
        for r in b["compounds"]:
            all_rows += 1
            all_within += int(bool(r["within_band"]))
            lb = by_lane.setdefault(str(r["lane"]), _bucket())
            for bucket in (pb, rb, lb):
                bucket["rows"] += 1
                bucket["within_band"] += int(bool(r["within_band"]))
                bucket["folds"].append(r["fold_error"])
            if literature:
                honest["rows"] += 1
                honest["within_band"] += int(bool(r["within_band"]))
                honest["folds"].append(r["fold_error"])
                seen_h = True
                if r["in_core_fit"]:
                    in_fit["rows"] += 1
                    in_fit["within_band"] += int(bool(r["within_band"]))
                    in_fit["folds"].append(r["fold_error"])
                    seen_f = True
                else:
                    out_of_sample["rows"] += 1
                    out_of_sample["within_band"] += int(bool(r["within_band"]))
                    out_of_sample["folds"].append(r["fold_error"])
                    seen_o = True
        honest["benchmarks"] += int(seen_h)
        out_of_sample["benchmarks"] += int(seen_o)
        in_fit["benchmarks"] += int(seen_f)
        if literature and b["compounds"]:
            honest["contract_passes"] += int(b["overall_status"] in {"pass", "pass-no-ranking"})
            honest["strict_ready"] += int(b["strict_ready"])
    for lane_bucket in by_lane.values():
        lane_bucket["benchmarks"] = None  # a benchmark can span lanes; rows are the unit here

    predictive_passes = sorted(
        b["benchmark_id"] for b in benches
        if b["evidence_role"] == "predictive" and b["overall_status"] in {"pass", "pass-no-ranking"}
    )
    strict_ready = sorted(b["benchmark_id"] for b in benches if b["strict_ready"])
    holdout = by_panel.get("maillard_path_holdout")
    role_totals = {role: bucket["benchmarks"] for role, bucket in sorted(by_role.items())}
    reports = sorted(data_paths.VALIDATION_DIR.glob("kinetic_core_b*_fit_report.json"))

    payload: Dict[str, Any] = {
        "artifact": "core_panel_scores",
        "engine": engine_metadata(),
        "pass_band_level": float(pass_band),
        "parameter_sources": [
            {"report": data_paths.rel(p), "sha256": _report_sha(p)} for p in reports
        ],
        "summary": {
            "panel_benchmark_count": len(benches),
            "scored_benchmark_count": sum(1 for b in benches if b["compounds"]),
            "refused_benchmark_count": sum(1 for b in benches if not b["compounds"]),
            "matched_compound_count": all_rows,
            "refused_compound_count": sum(len(b["refused_compounds"]) for b in benches),
            "within_band_count": all_within,
            "within_band_rate": (all_within / all_rows) if all_rows else None,
            "fold_summary": summarise_folds(
                [r["fold_error"] for b in benches for r in b["compounds"]]
            ),
            "evidence_role_totals": role_totals,
            "predictive_passes": predictive_passes,
            "strict_ready": strict_ready,
            "by_panel": {k: _close(v) for k, v in sorted(by_panel.items())},
            "by_evidence_role": {k: _close(v) for k, v in sorted(by_role.items())},
            "by_lane": {k: _close(v) for k, v in sorted(by_lane.items())},
            "honest_literature": {
                **_close(honest),
                "definition": (
                    "external-literature rows with fit-recovery benchmarks removed, by the "
                    "repository's leverage rule (src/fit_target_index.py): rows a LOW-leverage "
                    "fit saw stay in, annotated in_core_fit"
                ),
            },
            "out_of_sample": {
                **_close(out_of_sample),
                "definition": (
                    "honest_literature rows MINUS every row a core fit read (fit_targets."
                    "kinetic_core_*_fit_targets.json): the strongest claim the panel supports"
                ),
            },
            "in_core_fit": {
                **_close(in_fit),
                "definition": "literature rows the sulfur fit read; reported, never pooled",
            },
            "holdout_within_band": (
                None if holdout is None else {"hits": holdout["within_band"], "total": holdout["rows"]}
            ),
            "shared_with_holdout_panel": sum(
                1 for b in benches for r in b["compounds"] if r.get("shared_with")
            ),
        },
        "benchmarks": benches,
        "refused_compounds": [r for b in benches for r in b["refused_compounds"]],
        "skipped_bundles": skipped,
    }
    return payload


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _fmt(value: Optional[float]) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if abs(value) >= 1000 or (abs(value) < 0.01 and value != 0):
        return f"{value:.3g}"
    return f"{value:.3f}"


def render_markdown(payload: Dict[str, Any]) -> str:
    s = payload["summary"]
    fs = s["fold_summary"]
    out = [
        "# Core panel scorecard (kinetic core on the union panel)",
        "",
        f"pass band = {payload['pass_band_level']:.1f}x; contracts from each bundle's "
        "`validation_contract.scale_thresholds`, else the global default.",
        "",
        f"* panel: **{s['panel_benchmark_count']}** benchmarks, {s['scored_benchmark_count']} scored, "
        f"{s['refused_benchmark_count']} fully refused; rows **{s['matched_compound_count']}**, "
        f"refused rows {s['refused_compound_count']}",
        f"* within {payload['pass_band_level']:.0f}x: {s['within_band_count']}/{s['matched_compound_count']} "
        f"({_fmt(s['within_band_rate'])}); median fold {_fmt(fs['median_fold_error'])}, "
        f"geometric mean {_fmt(fs['geometric_mean_fold'])}, worst {_fmt(fs['worst_fold_error'])}",
        f"* evidence roles (core): {s['evidence_role_totals']}",
        f"* predictive benchmarks passing their contract: {', '.join(s['predictive_passes']) or 'NONE'}; "
        f"strict-ready: {', '.join(s['strict_ready']) or 'NONE'}",
    ]
    for key, label in (
        ("honest_literature", "honest literature"),
        ("out_of_sample", "out-of-sample"),
        ("in_core_fit", "rows the sulfur fit read"),
    ):
        b = s[key]
        out.append(
            f"* **{label}: {b['within_band']}/{b['rows']} within band** "
            f"({_fmt(b['within_band_rate'])}), {b['benchmarks']} benchmarks, "
            f"median fold {_fmt(b['fold_summary']['median_fold_error'])}, "
            f"geometric mean {_fmt(b['fold_summary']['geometric_mean_fold'])}"
        )
    out += ["", "## Per panel / role / lane", "",
            "| split | key | benchmarks | rows | within band | rate | contract passes | strict-ready | median fold | geo-mean fold |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for split in ("by_panel", "by_evidence_role", "by_lane"):
        for key, b in s[split].items():
            out.append(
                f"| {split[3:]} | {key} | {b['benchmarks'] if b['benchmarks'] is not None else '-'} | "
                f"{b['rows']} | {b['within_band']} | {_fmt(b['within_band_rate'])} | "
                f"{b['contract_passes']} | {b['strict_ready']} | "
                f"{_fmt(b['fold_summary']['median_fold_error'])} | "
                f"{_fmt(b['fold_summary']['geometric_mean_fold'])} |"
            )
    out += ["", "## Benchmarks", "",
            "| benchmark | panel | tier | role | rows | coverage | max ratio | mean log10 | "
            "contract (ratio / log10) | status | strict | in core fit |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for b in payload["benchmarks"]:
        c = b["scale_thresholds"]
        fit_rows = sum(1 for r in b["compounds"] if r["in_core_fit"])
        out.append(
            f"| {b['benchmark_id']} | {b['panel']} | {b['tier']} | {b['evidence_role']} | "
            f"{b['matched_compounds']}/{b['total_compounds']} | "
            f"{_fmt(b['coverage'])} | {_fmt(b['max_ratio'])} | {_fmt(b['mean_abs_log10_error'])} | "
            f"{c['max_ratio']:.2f} / {c['mean_abs_log10_error']:.3f} | {b['overall_status']} | "
            f"{'yes' if b['strict_ready'] else 'no'} | {fit_rows or '-'} |"
        )
    out += ["", "## Rows", "",
            "| benchmark | compound | unit | measured | predicted | fold | within band | within contract | "
            "interval (ug/L) | measured inside | lane | in core fit |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for b in payload["benchmarks"]:
        for r in b["compounds"]:
            interval = "-" if r["interval_ug_per_L"] is None else (
                f"[{_fmt(r['interval_ug_per_L'][0])}, {_fmt(r['interval_ug_per_L'][1])}]"
            )
            out.append(
                f"| {b['benchmark_id']} | {r['compound']} | {r['target_unit']} | {_fmt(r['measured'])} | "
                f"{_fmt(r['predicted'])} | {_fmt(r['fold_error'])} | {_fmt(r['within_band'])} | "
                f"{_fmt(r['within_contract_ratio'])} | {interval} | {_fmt(r['measured_within_interval'])} | "
                f"{r['lane']} | {'yes' if r['in_core_fit'] else 'no'}"
                f"{' (shared: ' + r['shared_with'] + ')' if r.get('shared_with') else ''} |"
            )
    if payload["refused_compounds"]:
        out += ["", "## Refused rows", "", "| benchmark | panel | compound | reason |", "|---|---|---|---|"]
        for r in payload["refused_compounds"]:
            out.append(
                f"| {r['benchmark_id']} | {r['panel']} | {r['compound']} | "
                f"{str(r['reason']).replace('|', '/')[:220]} |"
            )
    if payload.get("skipped_bundles"):
        out += ["", "## Bundles kept off the scored panel", ""]
        for r in payload["skipped_bundles"]:
            out.append(f"- {r['benchmark_id']} ({r['panel']}): {r['reason']}")
    out.append("")
    return "\n".join(out)


def write_artifact(payload: Dict[str, Any], json_path: Path) -> Tuple[Path, Path]:
    import json

    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    md_path = json_path.with_suffix(".md")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    return json_path, md_path


__all__ = [
    "OVERALL_STATUSES",
    "PASS_BAND_LEVEL",
    "fold_error",
    "render_markdown",
    "score_benchmark",
    "score_panel",
    "summarise_folds",
    "write_artifact",
]
