from __future__ import annotations

import json
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional

from src.benchmark_validation import DEFAULT_TARGET_TAG, evaluate_benchmark, get_benchmark_files, get_benchmark_metadata, load_benchmark, summarize_evaluation

# 2026-08-27 (Wave I): `get_barrier()` resolves BARRIER_OFFSETS against the
# *normalised raw* reaction-family label (`_normalize_family_key`) BEFORE it
# canonicalises that label (`_canonical_fast_family`). We therefore have to
# canonicalise here, on the producing side, to know which offset keys can
# actually reach a step. These two helpers are private to `barrier_constants`
# but re-implementing them here would let the two copies drift silently, which
# is exactly the failure mode this fix is repairing.
from src.barrier_constants import _canonical_fast_family, _normalize_family_key


DEFAULT_FAMILY_OFFSET_KEYS: Dict[str, str] = {
    "schiff_condensation": "schiff",
    "amadori_rearrangement": "amadori",
    "1,2-enolisation": "enol",
    "strecker_degradation": "strecker",
    "cysteine_thermolysis": "cys",
    "thiol_addition": "thiol_addition",
    "thiol_oxidation": "thiol_oxidation",
    "aminoketone_condensation": "aminoketone_condensation",
    "retro_aldol": "retro_aldol",
}

# Raw `reaction_family` labels actually emitted by the template engine
# (`src/reaction_templates.py`) and the curated pathway layer
# (`src/curated_pathways.py`). BARRIER_OFFSETS is keyed on the NORMALISED form
# of these strings, not on the canonical family name, so the sensitivity sweep
# has to expand a canonical family into the labels that canonicalise to it.
#
# 2026-08-27 (Wave I): this list is asserted to be a superset of the labels
# found by AST-scanning `src/` in
# tests/unit/test_wave_i_tooling.py::test_engine_family_labels_cover_sources,
# so a newly added family cannot silently reintroduce a false zero.
ENGINE_FAMILY_LABELS: tuple[str, ...] = (
    "Additive_Thermal_Degradation",
    "Amadori_Rearrangement",
    "Aminoketone_Condensation",
    "Beta_Elimination",
    "Cysteine_Degradation",
    "DHA_Crosslinking",
    # Wave N 2026-08-27: corrected MFT route (Cerny & Davidek 2003/2004).
    "Deoxyosone_Reduction",
    "Enolisation",
    "Enolisation_1_2",
    "Enolisation_2_3",
    "Enolisation_2_3_Amadori",
    "Enolisation_Intermediate",
    "Fructofuranosyl_Dehydration",  # Wave P 2026-08-27: fructose HMF route correction
    "Furan_Ring_Aromatisation",
    "Furanone_Amino_Acid_Reduction",  # Wave P 2026-08-27: DMHF [HH] gate removed
    "Furanone_Cyclisation",
    "Furanone_Formation",
    "Furanone_Reductive_Opening",  # Wave P 2026-08-27: norfuraneol -> 2,3-pentanedione
    # Wave X 2026-08-28: norfuraneol + H2S + 2[H] -> MFT, re-added as a SLOW
    # PARALLEL channel constrained by Hofmann & Schieberle 1998 Table 4. NOT the
    # retired `Thiol_Addition_Norfuraneol`; new name, new barrier key.
    "Furanone_Reductive_Sulfhydrylation",
    "Generalized_Deamination",
    # Wave T4 2026-08-27: the ketose rearrangement. Emitted since before Wave I
    # (`reaction_templates.py:60`) and missing here the whole time, because the
    # guard test below AST-scans only string literals passed to
    # `ElementaryStep(...)` and this one is bound to a local variable first.
    "Heyns_Rearrangement",
    "Lipid_Homolysis",
    "Lipid_Schiff_Base",
    "Lipid_Strecker_Synergy",
    "Lipid_Thiazole_Condensation",
    # Wave P 2026-08-27: the C2+C3 recombination lane (Hofmann & Schieberle 1998).
    "Mercaptoketone_Aldol_Addition",
    "Mercaptoketone_Cyclodehydration",
    "Mercaptoketone_Formation",
    "Retro_Aldol_Fragmentation",
    "Safety_Risk_AGE",
    "Safety_Risk_Acrylamide",
    "Schiff_Base_Formation",
    "Strecker_Degradation",
    "Sugar_Dehydration",
    "Sugar_Ring_Opening",
    "Thiohemiacetal_Formation",
    "Thiol_Addition",
    "Thiol_Addition_H2",
    "Thiol_Addition_Hexose_Legacy_Shortcut",
    "Thiol_Addition_Legacy_Shortcut",
    "Thiol_Addition_Norfuraneol",  # retired Wave N 2026-08-27; kept: superset is harmless
    "Thiol_Addition_Pentodiulose",  # Wave N 2026-08-27: corrected MFT route
    "Thiol_Dehydration",
    "Thiol_Oxidation",
)


def resolve_offset_keys(reaction_family: str, offset_key: str) -> tuple[str, ...]:
    """Expand a canonical family + its short offset key into every BARRIER_OFFSETS key that bites.

    2026-08-27 (Wave I) -- FIX 6. The sweep used to publish a single short key
    (`schiff`, `enol`, `cys`, ...) into BARRIER_OFFSETS. `get_barrier()` matches
    those short keys through a hard-coded `offset_map` whose *values* must appear
    as a SUBSTRING of the normalised raw family label, and it does that match
    before canonicalisation. The result:

      * `schiff` -> needs "schiff_condensation" inside the label, but the engine
        emits `Schiff_Base_Formation` -> "schiff_base_formation". No match.
      * `enol` -> needs the literal "1,2-enolisation" inside the label, but
        normalisation has already turned every "-" into "_". Can never match.
      * `cys` -> needs "cysteine_thermolysis" inside the label, but the engine
        emits `Cysteine_Degradation`. No match.

    Those three offsets were therefore exact no-ops, and the sweep reported a
    sensitivity of 0.00 for them -- a false zero, not a finding. `retro_aldol`
    and `thiol_addition` were no-ops for the same reason via the exact-match
    branch (`Retro_Aldol_Fragmentation`, `Thiol_Addition_H2`, ...).

    The fix keeps the short key (harmless, and preserves the pre-existing
    behaviour for `amadori`/`strecker`, which did match) and adds the canonical
    family name plus every normalised engine label that canonicalises to it.
    """
    canonical_target = _canonical_fast_family(reaction_family) or _normalize_family_key(reaction_family)
    keys = {
        str(offset_key),
        str(reaction_family),
        _normalize_family_key(reaction_family),
        _normalize_family_key(canonical_target),
    }
    for label in ENGINE_FAMILY_LABELS:
        if _canonical_fast_family(label) == canonical_target:
            keys.add(_normalize_family_key(label))
    return tuple(sorted(key for key in keys if key))


def _status_score(status: str) -> int:
    return {
        "pass": 0,
        "pass-no-ranking": 1,
        "coverage-gap": 2,
        "ranking-gap": 3,
        "scale-gap": 4,
        "unsupported": 5,
    }.get(str(status), 4)


def _execution_weight(execution_path: str) -> float:
    normalized = str(execution_path).strip().lower()
    if normalized == "matrix_precursor_augmented":
        return 1.0
    if normalized == "free_precursor":
        return 0.65
    return 0.5


@contextmanager
def _barrier_offsets(offsets: Mapping[str, float]) -> Iterator[None]:
    previous = os.environ.get("BARRIER_OFFSETS")
    try:
        os.environ["BARRIER_OFFSETS"] = json.dumps(dict(offsets))
        yield
    finally:
        if previous is None:
            os.environ.pop("BARRIER_OFFSETS", None)
        else:
            os.environ["BARRIER_OFFSETS"] = previous


def build_family_sensitivity_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    delta_kcal: float = 3.0,
    family_offset_keys: Optional[Mapping[str, str]] = None,
) -> Dict[str, Any]:
    family_map = dict(family_offset_keys or DEFAULT_FAMILY_OFFSET_KEYS)
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    benchmark_contexts: List[Dict[str, Any]] = []

    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"free_precursor", "matrix_precursor_augmented"}:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        benchmark_contexts.append(
            {
                "bench_file": Path(bench_file),
                "benchmark_id": summary.benchmark_id,
                "execution_path": summary.execution_path,
                "protein_type": str(bench.get("protein_type", "free")),
                "baseline_summary": summary,
                "weight": _execution_weight(summary.execution_path),
            }
        )

    if not benchmark_contexts:
        return {
            "summary": {
                "evaluated_benchmark_count": 0,
                "family_count": 0,
                "delta_kcal": float(delta_kcal),
            },
            "families": [],
        }

    family_rows: List[Dict[str, Any]] = []
    for family, offset_key in family_map.items():
        # 2026-08-27 (Wave I) -- FIX 6: see `resolve_offset_keys`. Writing the
        # short key alone left schiff/enolisation/cysteine (and retro_aldol,
        # thiol_addition) as silent no-ops.
        offset_keys = resolve_offset_keys(family, offset_key)
        scenario_rows: List[Dict[str, Any]] = []
        max_weighted_impact = 0.0
        dominant_direction = "none"
        affected_benchmark_ids: set[str] = set()

        for direction_label, signed_delta in (("down", -abs(delta_kcal)), ("up", abs(delta_kcal))):
            scenario_details: List[Dict[str, Any]] = []
            status_change_count = 0
            ranking_contract_change_count = 0
            matrix_benchmark_change_count = 0
            mean_abs_relative_mae_shift = 0.0
            relative_mae_terms = 0
            weighted_impact_score = 0.0

            with _barrier_offsets({key: signed_delta for key in offset_keys}):
                for context in benchmark_contexts:
                    scenario_eval = evaluate_benchmark(context["bench_file"], target_tag=target_tag)
                    scenario_summary = summarize_evaluation(
                        scenario_eval,
                        protein_type=context["protein_type"],
                    )
                    baseline = context["baseline_summary"]
                    status_changed = scenario_summary.overall_status != baseline.overall_status
                    ranking_changed = scenario_summary.ranking_contract_status != baseline.ranking_contract_status
                    baseline_mae = float(baseline.mae_ppb or 0.0)
                    scenario_mae = float(scenario_summary.mae_ppb or 0.0)
                    mae_shift_ppb = scenario_mae - baseline_mae
                    relative_mae_shift = abs(mae_shift_ppb) / max(abs(baseline_mae), 1.0)
                    local_impact = (
                        abs(_status_score(scenario_summary.overall_status) - _status_score(baseline.overall_status))
                        + (1.0 if ranking_changed else 0.0)
                        + min(relative_mae_shift, 3.0)
                    )
                    weighted_impact_score += context["weight"] * local_impact
                    if status_changed:
                        status_change_count += 1
                    if ranking_changed:
                        ranking_contract_change_count += 1
                    if status_changed or ranking_changed or relative_mae_shift > 0.05:
                        affected_benchmark_ids.add(str(context["benchmark_id"]))
                    if context["execution_path"] == "matrix_precursor_augmented" and (status_changed or ranking_changed or relative_mae_shift > 0.05):
                        matrix_benchmark_change_count += 1
                    mean_abs_relative_mae_shift += relative_mae_shift
                    relative_mae_terms += 1
                    scenario_details.append(
                        {
                            "benchmark_id": str(context["benchmark_id"]),
                            "execution_path": str(context["execution_path"]),
                            "baseline_status": str(baseline.overall_status),
                            "scenario_status": str(scenario_summary.overall_status),
                            "baseline_ranking_contract_status": str(baseline.ranking_contract_status),
                            "scenario_ranking_contract_status": str(scenario_summary.ranking_contract_status),
                            "mae_shift_ppb": float(mae_shift_ppb),
                            "relative_mae_shift": float(relative_mae_shift),
                            "status_changed": bool(status_changed),
                            "ranking_contract_changed": bool(ranking_changed),
                        }
                    )

            mean_abs_relative_mae_shift = mean_abs_relative_mae_shift / max(relative_mae_terms, 1)
            scenario_row = {
                "direction": direction_label,
                "offset_kcal": float(signed_delta),
                "status_change_count": int(status_change_count),
                "ranking_contract_change_count": int(ranking_contract_change_count),
                "matrix_benchmark_change_count": int(matrix_benchmark_change_count),
                "mean_abs_relative_mae_shift": float(mean_abs_relative_mae_shift),
                "weighted_impact_score": float(weighted_impact_score),
                "benchmarks": scenario_details,
            }
            scenario_rows.append(scenario_row)
            if weighted_impact_score > max_weighted_impact:
                max_weighted_impact = float(weighted_impact_score)
                dominant_direction = direction_label

        family_rows.append(
            {
                "reaction_family": family,
                "offset_key": offset_key,
                # 2026-08-27 (Wave I): the keys actually written into
                # BARRIER_OFFSETS, so a reader can tell a real zero from a
                # routing no-op.
                "offset_keys": list(offset_keys),
                "dominant_direction": dominant_direction,
                "max_weighted_impact_score": float(max_weighted_impact),
                "affected_benchmark_ids": sorted(affected_benchmark_ids),
                "scenarios": scenario_rows,
            }
        )

    family_rows.sort(key=lambda row: (-float(row["max_weighted_impact_score"]), row["reaction_family"]))
    return {
        "summary": {
            "evaluated_benchmark_count": len(benchmark_contexts),
            "family_count": len(family_rows),
            "delta_kcal": float(delta_kcal),
            "benchmark_ids": [str(context["benchmark_id"]) for context in benchmark_contexts],
        },
        "families": family_rows,
    }


def render_family_sensitivity_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Sensitivity",
        "",
        "| Reaction Family | Dominant Direction | Max Weighted Impact | Affected Benchmarks | Up Score | Down Score |",
        "| --- | --- | ---: | --- | ---: | ---: |",
    ]
    for row in payload.get("families", []):
        scenarios = {str(item.get("direction")): item for item in row.get("scenarios", [])}
        up_score = float((scenarios.get("up") or {}).get("weighted_impact_score", 0.0))
        down_score = float((scenarios.get("down") or {}).get("weighted_impact_score", 0.0))
        lines.append(
            f"| {row['reaction_family']} | {row.get('dominant_direction', 'none')} | {float(row.get('max_weighted_impact_score', 0.0)):.2f} | {', '.join(row.get('affected_benchmark_ids', [])[:4]) or 'none'} | {up_score:.2f} | {down_score:.2f} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmarks evaluated: {int(summary.get('evaluated_benchmark_count', 0))}",
            f"Reaction families screened: {int(summary.get('family_count', 0))}",
            f"Offset magnitude screened: +/-{float(summary.get('delta_kcal', 0.0)):.1f} kcal/mol",
        ]
    )
    return "\n".join(lines) + "\n"