from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


ROOT = _repo_root()
DEFAULT_TARGET_REGISTRY_PATH = ROOT / "data" / "lit" / "computational_gap_closure_targets.json"
DEFAULT_VALIDATION_OUTPUT_DIR = ROOT / "results" / "validation"
DEFAULT_COMPUTATIONAL_GAP_REFINEMENT_OUTPUT_DIR = ROOT / "results" / "computational_gap_refinement"
DEFAULT_PRIORS_PATH = ROOT / "data" / "lit" / "computational_priors.json"
DEFAULT_DFT_METHOD_CHAIN = "wB97M-V/def2-tzvp // r2SCAN-3c/def2-svp + ddCOSMO(water)"

COMPUTATIONAL_GAP_PRIORITY_TARGET_IDS = (
    "hexanal_radical_quench",
    "lysinoalanine_crosslink",
    "aa_ring_open_dicarbonyl",
    "asparagine_sugar_explicit_water_cluster",
)

DEFAULT_WRITE_BACK_ARTIFACTS = (
    "data/lit/computational_priors.json",
    "results/validation/computational_gap_dft_ingestion_report.md",
    "results/validation/computational_gap_dft_ingestion_report.json",
    "results/validation/computational_gap_closure_scorecard.md",
    "results/validation/computational_gap_closure_scorecard.json",
    "results/validation/refinement_governance.md",
    "results/validation/refinement_governance.json",
)

DEFAULT_INGESTION_CEILING_BY_PROMOTION = {
    "ranking_only": "ranking_only_support",
    "bounded_calibration": "bounded_calibration",
    "uncertainty_narrowing_only": "uncertainty_narrowing_only",
}


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _to_repo_relative(path: Path) -> str:
    return path.resolve(strict=False).relative_to(ROOT).as_posix()


def _path_exists(relative_path: str) -> bool:
    return (ROOT / relative_path).exists()


def _xyz_atom_count(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                return int(stripped)
            except ValueError:
                return None
    return None


def load_computational_gap_targets(path: Optional[Path | str] = None) -> Dict[str, Any]:
    target_path = Path(path) if path is not None else DEFAULT_TARGET_REGISTRY_PATH
    return _load_json(target_path)


def _explicit_water_target() -> Dict[str, Any]:
    return {
        "id": "asparagine_sugar_explicit_water_cluster",
        "reaction_key": "asparagine_sugar_explicit_water_cluster",
        "slr_family": "13_safety",
        "gap": "Safety lane: proton-transfer-sensitive asparagine-sugar transition state under wet or HME-like conditions remains too broad with implicit-solvent-only support.",
        "decision_impact": "high",
        "cheap_first_method": "seed the cluster in xTB with a minimal explicit-water motif set and preserve the result as uncertainty narrowing only",
        "upgrade_method": "r2SCAN-3c refinement of a minimal explicit-water cluster before any single-point comparison",
        "acceptance_signal": "the explicit-water result narrows the moisture-sensitive acrylamide uncertainty band without claiming matrix closure",
        "method": "explicit-water cluster seed generation pending",
        "source_quality": "seed_required",
        "barrier_kcal_mol": None,
        "barrier_kj_mol": None,
        "uncertainty_kj": None,
        "promotion_ceiling": "uncertainty_narrowing_only",
        "dft_status": "seed_required",
        "active_arrhenius_key": "asparagine_sugar_explicit_water_cluster_pending",
    }


def _promotion_honest_label(promotion_ceiling: str) -> str:
    normalized = str(promotion_ceiling).strip()
    if normalized == "ranking_only":
        return "Selective DFT anchor, ranking support only until matrix benchmark closure"
    if normalized == "bounded_calibration":
        return "Selective DFT anchor, bounded calibration only until matrix benchmark closure"
    if normalized == "uncertainty_narrowing_only":
        return "Selective DFT anchor, uncertainty narrowing only until matrix benchmark closure"
    return "Selective DFT anchor, benchmark closure still withheld pending matrix observables"


def _build_surrogate_payload(source_row: Mapping[str, Any]) -> Dict[str, Any]:
    kind = str(source_row.get("surrogate_kind", "")).strip()
    short_label = str(source_row.get("surrogate_short_label", kind)).strip()
    reaction = str(source_row.get("surrogate_reaction", "")).strip()
    assumption = str(source_row.get("surrogate_assumption", "")).strip()
    note = str(source_row.get("surrogate_note", "")).strip()
    basis = [
        str(item).strip()
        for item in source_row.get("surrogate_basis", [])
        if str(item).strip()
    ]
    barrier_kcal = source_row.get("surrogate_barrier_kcal_mol", source_row.get("barrier_kcal_mol"))
    barrier_kj = source_row.get("surrogate_barrier_kj_mol", source_row.get("barrier_kj_mol"))
    available = any([kind, reaction, assumption, note, basis])
    return {
        "available": bool(available),
        "kind": kind if available else "",
        "short_label": short_label if available else "",
        "reaction": reaction if available else "",
        "assumption": assumption if available else "",
        "note": note if available else "",
        "basis": basis if available else [],
        "barrier_kcal_mol": barrier_kcal if available else None,
        "barrier_kj_mol": barrier_kj if available else None,
    }


def _build_target_row(source_row: Mapping[str, Any], *, priority_rank: int) -> Dict[str, Any]:
    reaction_key = str(source_row.get("reaction_key", source_row.get("id", ""))).strip()
    target_id = str(source_row.get("id", reaction_key)).strip()
    geometry_dir = ROOT / "data" / "geometries" / "xtb_inputs" / reaction_key
    geometry_dir_relative = _to_repo_relative(geometry_dir)
    reactant_relative = f"{geometry_dir_relative}/reactant.xyz"
    product_relative = f"{geometry_dir_relative}/product.xyz"
    ts_guess_relative = f"{geometry_dir_relative}/xtbpath_ts.xyz"
    run_script_relative = f"{geometry_dir_relative}/run_xtb.sh"
    is_explicit_water = target_id == "asparagine_sugar_explicit_water_cluster"
    reactant_atom_count = _xyz_atom_count(ROOT / reactant_relative)
    product_atom_count = _xyz_atom_count(ROOT / product_relative)
    atom_count_match = (
        reactant_atom_count is not None
        and product_atom_count is not None
        and reactant_atom_count == product_atom_count
    )

    xtb_ready = (
        not is_explicit_water
        and _path_exists(reactant_relative)
        and _path_exists(product_relative)
        and _path_exists(run_script_relative)
        and atom_count_match
    )
    if is_explicit_water:
        xtb_status = "seed_required"
        xtb_blocker = "Explicit-water cluster seed generation is still required before xTB path search."
    elif not _path_exists(reactant_relative) or not _path_exists(product_relative) or not _path_exists(run_script_relative):
        xtb_status = "seed_required"
        xtb_blocker = "Reactant/product geometries or the xTB runner script are still missing."
    elif not atom_count_match:
        xtb_status = "blocked_atom_count_mismatch"
        xtb_blocker = (
            f"Reactant and product XYZ files contain different atom counts "
            f"({reactant_atom_count} vs {product_atom_count}); xTB path search requires a balanced mapped pair."
        )
    else:
        xtb_status = "ready"
        xtb_blocker = None

    dft_ready = xtb_ready and _path_exists(ts_guess_relative)
    if is_explicit_water:
        dft_status = "seed_required"
    elif xtb_status == "blocked_atom_count_mismatch":
        dft_status = "blocked_atom_count_mismatch"
    elif dft_ready:
        dft_status = "ready_for_dft"
    else:
        dft_status = "blocked_missing_xtb_ts"

    promotion_ceiling = str(source_row.get("promotion_ceiling", "ranking_only")).strip() or "ranking_only"
    surrogate = _build_surrogate_payload(source_row)

    return {
        "id": target_id,
        "priority_rank": int(priority_rank),
        "reaction_key": reaction_key,
        "slr_family": str(source_row.get("slr_family", "unknown")),
        "gap": str(source_row.get("gap", "")),
        "decision_impact": str(source_row.get("decision_impact", "unknown")),
        "cheap_first_method": str(source_row.get("cheap_first_method", "")),
        "upgrade_method": str(source_row.get("upgrade_method", "")),
        "acceptance_signal": str(source_row.get("acceptance_signal", "")),
        "method": str(source_row.get("method", source_row.get("cheap_first_method", ""))),
        "source_quality": str(source_row.get("source_quality", "unknown")),
        "current_barrier_kcal_mol": source_row.get("barrier_kcal_mol"),
        "current_barrier_kj_mol": source_row.get("barrier_kj_mol"),
        "uncertainty_kj": source_row.get("uncertainty_kj"),
        "promotion_ceiling": promotion_ceiling,
        "maximum_claim_level": DEFAULT_INGESTION_CEILING_BY_PROMOTION.get(promotion_ceiling, promotion_ceiling),
        "active_arrhenius_key": str(source_row.get("active_arrhenius_key", "")),
        "surrogate": surrogate,
        "governance_status": "hold_observable_first",
        "governance_note": "Computational-gap refinement jobs can narrow priors and uncertainty only; they do not promote benchmark closure while observable-first governance remains active.",
        "write_back_artifacts": list(DEFAULT_WRITE_BACK_ARTIFACTS),
        "xtb": {
            "required": True,
            "status": xtb_status,
            "blocker_note": xtb_blocker,
            "surrogate_available": bool(surrogate.get("available", False)),
            "geometry_dir": geometry_dir_relative,
            "runner_script": None if is_explicit_water else run_script_relative,
            "required_inputs": [reactant_relative, product_relative],
            "reactant_atom_count": reactant_atom_count,
            "product_atom_count": product_atom_count,
            "atom_count_match": atom_count_match,
            "expected_outputs": [f"{geometry_dir_relative}/xtbpath.xyz", ts_guess_relative],
            "execute_with": f"python scripts/run_computational_gap_xtb.py --target {target_id} --execute",
        },
        "mlp": {
            "allowed_candidate_id": "mace_mp_small",
            "allowed_role": "geom_preopt",
            "authority_allowed": False,
            "blocked_candidates": ["mace_off_medium", "aimnet2_shortlist", "orbmol_shortlist"],
            "note": "MLP barrier authority remains blocked; only geometry preoptimization is allowed where benchmarked.",
        },
        "dft": {
            "required": True,
            "status": dft_status,
            "surrogate_available": bool(surrogate.get("available", False)),
            "runner_script": "scripts/run_computational_gap_dft.py",
            "reactant_path": reactant_relative,
            "ts_guess_path": ts_guess_relative,
            "solvent_name": "water",
            "temperature_k": 423.15,
            "charge": int(source_row.get("charge", 0) or 0),
            "spin": int(source_row.get("spin", 0) or 0),
            "use_explicit_solvent": bool(is_explicit_water),
            "n_water": 3 if is_explicit_water else 0,
            "expected_outputs": [
                "optimized_reactant.xyz",
                "optimized_ts.xyz",
                "barrier_summary.json",
                "quasi_harmonic_thermo.json",
                "surrogate_patch.json",
            ],
            "execute_with": f"python scripts/run_computational_gap_dft.py --target {target_id} --execute",
        },
    }


def build_computational_gap_refinement_plan_artifact(
    target_registry_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    registry = load_computational_gap_targets(target_registry_path)
    source_lookup = {
        str(row.get("id", "")).strip(): row
        for row in registry.get("targets", [])
        if str(row.get("id", "")).strip()
    }
    source_lookup["asparagine_sugar_explicit_water_cluster"] = _explicit_water_target()

    rows = [
        _build_target_row(source_lookup[target_id], priority_rank=index)
        for index, target_id in enumerate(COMPUTATIONAL_GAP_PRIORITY_TARGET_IDS, start=1)
    ]

    xtb_ready_count = sum(1 for row in rows if row["xtb"]["status"] == "ready")
    dft_ready_count = sum(1 for row in rows if row["dft"]["status"] == "ready_for_dft")
    explicit_solvent_count = sum(1 for row in rows if row["dft"]["use_explicit_solvent"])

    return {
        "generated_at": _utc_timestamp(),
        "summary": {
            "target_count": len(rows),
            "xtb_job_count": len(rows),
            "xtb_ready_count": xtb_ready_count,
            "dft_job_count": len(rows),
            "dft_ready_count": dft_ready_count,
            "explicit_solvent_job_count": explicit_solvent_count,
            "governance_status": "hold_observable_first",
            "maximum_promotion_without_wet_lab": "prior_narrowing_only",
        },
        "write_back_surfaces": list(DEFAULT_WRITE_BACK_ARTIFACTS),
        "source_registry_path": _to_repo_relative(Path(target_registry_path) if target_registry_path is not None else DEFAULT_TARGET_REGISTRY_PATH),
        "targets": rows,
    }


def build_computational_gap_xtb_job_manifest(
    plan_payload: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    payload = dict(plan_payload) if plan_payload is not None else build_computational_gap_refinement_plan_artifact()
    jobs = []
    for row in payload.get("targets", []):
        xtb = row.get("xtb", {})
        jobs.append(
            {
                "target_id": row.get("id"),
                "priority_rank": row.get("priority_rank"),
                "reaction_key": row.get("reaction_key"),
                "slr_family": row.get("slr_family"),
                "status": xtb.get("status"),
                "geometry_dir": xtb.get("geometry_dir"),
                "runner_script": xtb.get("runner_script"),
                "required_inputs": list(xtb.get("required_inputs", [])),
                "reactant_atom_count": xtb.get("reactant_atom_count"),
                "product_atom_count": xtb.get("product_atom_count"),
                "atom_count_match": xtb.get("atom_count_match"),
                "blocker_note": xtb.get("blocker_note"),
                "surrogate_available": bool(row.get("surrogate", {}).get("available", False)),
                "surrogate_kind": row.get("surrogate", {}).get("kind"),
                "surrogate_short_label": row.get("surrogate", {}).get("short_label"),
                "surrogate_barrier_kcal_mol": row.get("surrogate", {}).get("barrier_kcal_mol"),
                "surrogate_note": row.get("surrogate", {}).get("note"),
                "expected_outputs": list(xtb.get("expected_outputs", [])),
                "promotion_ceiling": row.get("promotion_ceiling"),
            }
        )
    return {
        "generated_at": _utc_timestamp(),
        "description": "Computational-gap xTB execution manifest for no-wet-lab narrowing.",
        "summary": {
            "job_count": len(jobs),
            "ready_count": sum(1 for job in jobs if job["status"] == "ready"),
            "seed_required_count": sum(1 for job in jobs if job["status"] == "seed_required"),
            "blocked_atom_count_mismatch_count": sum(1 for job in jobs if job["status"] == "blocked_atom_count_mismatch"),
        },
        "jobs": jobs,
    }


def build_computational_gap_dft_job_manifest(
    plan_payload: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    payload = dict(plan_payload) if plan_payload is not None else build_computational_gap_refinement_plan_artifact()
    jobs = []
    for row in payload.get("targets", []):
        dft = row.get("dft", {})
        jobs.append(
            {
                "target_id": row.get("id"),
                "priority_rank": row.get("priority_rank"),
                "reaction_key": row.get("reaction_key"),
                "slr_family": row.get("slr_family"),
                "status": dft.get("status"),
                "reactant_path": dft.get("reactant_path"),
                "ts_guess_path": dft.get("ts_guess_path"),
                "runner_script": dft.get("runner_script"),
                "solvent_name": dft.get("solvent_name"),
                "temperature_k": dft.get("temperature_k"),
                "charge": dft.get("charge"),
                "spin": dft.get("spin"),
                "use_explicit_solvent": dft.get("use_explicit_solvent"),
                "n_water": dft.get("n_water"),
                "blocker_note": row.get("xtb", {}).get("blocker_note"),
                "surrogate_available": bool(row.get("surrogate", {}).get("available", False)),
                "surrogate_kind": row.get("surrogate", {}).get("kind"),
                "surrogate_short_label": row.get("surrogate", {}).get("short_label"),
                "surrogate_barrier_kcal_mol": row.get("surrogate", {}).get("barrier_kcal_mol"),
                "surrogate_note": row.get("surrogate", {}).get("note"),
                "expected_outputs": list(dft.get("expected_outputs", [])),
                "promotion_ceiling": row.get("promotion_ceiling"),
                "write_back_artifacts": list(row.get("write_back_artifacts", [])),
            }
        )
    return {
        "generated_at": _utc_timestamp(),
        "description": "Computational-gap DFT execution manifest for no-wet-lab narrowing.",
        "summary": {
            "job_count": len(jobs),
            "ready_for_dft_count": sum(1 for job in jobs if job["status"] == "ready_for_dft"),
            "blocked_missing_xtb_ts_count": sum(1 for job in jobs if job["status"] == "blocked_missing_xtb_ts"),
            "blocked_atom_count_mismatch_count": sum(1 for job in jobs if job["status"] == "blocked_atom_count_mismatch"),
            "seed_required_count": sum(1 for job in jobs if job["status"] == "seed_required"),
        },
        "jobs": jobs,
    }


def build_computational_gap_dft_ingestion_artifact(
    *,
    manifest_payload: Optional[Mapping[str, Any]] = None,
    execution_payload: Optional[Mapping[str, Any]] = None,
    execution_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    manifest = dict(manifest_payload) if manifest_payload is not None else build_computational_gap_dft_job_manifest()
    execution: Dict[str, Any] = {}
    if execution_payload is not None:
        execution = dict(execution_payload)
    elif execution_path is not None and Path(execution_path).exists():
        execution = _load_json(Path(execution_path))

    execution_lookup = {
        str(row.get("target_id", "")).strip(): row
        for row in execution.get("jobs", [])
        if str(row.get("target_id", "")).strip()
    }
    rows: List[Dict[str, Any]] = []
    for job in manifest.get("jobs", []):
        target_id = str(job.get("target_id", ""))
        execution_row = execution_lookup.get(target_id, {})
        execution_status = str(execution_row.get("status", "not_run"))
        rows.append(
            {
                "target_id": target_id,
                "reaction_key": str(job.get("reaction_key", "")),
                "slr_family": str(job.get("slr_family", "unknown")),
                "manifest_status": str(job.get("status", "unknown")),
                "execution_status": execution_status,
                "available_result": execution_status == "completed",
                "barrier_kcal_mol": execution_row.get("barrier_kcal_mol"),
                "fast_mode": bool(execution_row.get("fast_mode", False)),
                "promotion_ready": bool(execution_row.get("promotion_ready", False)),
                "promotion_ceiling": str(job.get("promotion_ceiling", "unknown")),
                "write_back_artifacts": list(job.get("write_back_artifacts", [])),
            }
        )

    return {
        "generated_at": _utc_timestamp(),
        "summary": {
            "manifest_job_count": len(rows),
            "available_result_count": sum(1 for row in rows if row["available_result"]),
            "promotable_result_count": sum(1 for row in rows if row["promotion_ready"]),
            "ready_for_dft_count": sum(1 for row in rows if row["manifest_status"] == "ready_for_dft"),
            "execution_file_present": bool(execution),
        },
        "jobs": rows,
    }


def build_computational_gap_dft_promotion_payload(
    *,
    priors_payload: Mapping[str, Any],
    execution_payload: Mapping[str, Any],
) -> Dict[str, Any]:
    prior_entries = {
        str(row.get("reaction_key", "")).strip(): row
        for row in priors_payload.get("dft_kinetic_priors", {}).get("entries", [])
        if str(row.get("reaction_key", "")).strip()
    }
    promoted_rows: List[Dict[str, Any]] = []
    for row in execution_payload.get("jobs", []):
        reaction_key = str(row.get("reaction_key", "")).strip()
        if not reaction_key or reaction_key not in prior_entries:
            continue
        if str(row.get("status", "")) != "completed":
            continue
        if bool(row.get("fast_mode", False)):
            continue
        if not bool(row.get("promotion_ready", True)):
            continue
        barrier_kcal = row.get("barrier_kcal_mol")
        if barrier_kcal is None:
            continue
        prior_row = prior_entries[reaction_key]
        promoted_rows.append(
            {
                "reaction_key": reaction_key,
                "previous_tier": str(prior_row.get("current_tier", "unknown")),
                "new_tier": "selective_dft_anchor",
                "barrier_kcal_mol": float(barrier_kcal),
                "barrier_kj_mol": round(float(barrier_kcal) * 4.184, 2),
                "promotion_ceiling": str(prior_row.get("promotion_ceiling", row.get("promotion_ceiling", "unknown"))),
                "method_chain": str(row.get("method_chain", DEFAULT_DFT_METHOD_CHAIN)),
            }
        )
    return {
        "generated_at": _utc_timestamp(),
        "summary": {
            "candidate_result_count": len(list(execution_payload.get("jobs", []))),
            "promoted_count": len(promoted_rows),
        },
        "promotions": promoted_rows,
    }


def apply_computational_gap_dft_promotions(
    *,
    execution_payload: Optional[Mapping[str, Any]] = None,
    execution_path: Optional[Path | str] = None,
    priors_path: Optional[Path | str] = None,
    write_changes: bool = True,
) -> Dict[str, Any]:
    resolved_priors_path = Path(priors_path) if priors_path is not None else DEFAULT_PRIORS_PATH
    priors_payload = _load_json(resolved_priors_path)

    execution: Dict[str, Any] = {}
    if execution_payload is not None:
        execution = dict(execution_payload)
    elif execution_path is not None and Path(execution_path).exists():
        execution = _load_json(Path(execution_path))

    promotion_payload = build_computational_gap_dft_promotion_payload(
        priors_payload=priors_payload,
        execution_payload=execution,
    )
    promoted_lookup = {row["reaction_key"]: row for row in promotion_payload.get("promotions", [])}
    if promoted_lookup:
        for row in priors_payload.get("dft_kinetic_priors", {}).get("entries", []):
            reaction_key = str(row.get("reaction_key", "")).strip()
            promotion_row = promoted_lookup.get(reaction_key)
            if promotion_row is None:
                continue
            row["current_tier"] = promotion_row["new_tier"]
            row["target_tier"] = promotion_row["new_tier"]
            row["barrier_kcal_mol"] = promotion_row["barrier_kcal_mol"]
            row["barrier_kj_mol"] = promotion_row["barrier_kj_mol"]
            row["computational_method"] = promotion_row["method_chain"]
            row["honest_label"] = _promotion_honest_label(str(promotion_row["promotion_ceiling"]))
            row["promotion_note"] = "Promoted from completed computational-gap DFT execution; benchmark closure remains gated by matrix observables."
            row["last_refined_at"] = promotion_payload["generated_at"]
        if write_changes:
            _write_json(resolved_priors_path, priors_payload)
    promotion_payload["summary"]["wrote_priors"] = bool(promoted_lookup and write_changes)
    promotion_payload["summary"]["priors_path"] = _to_repo_relative(resolved_priors_path)
    return promotion_payload


def render_computational_gap_refinement_plan_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Computational Gap Refinement Plan",
        "",
        f"Generated: {payload.get('generated_at', 'unknown')}",
        f"Targets: {int(summary.get('target_count', 0))}",
        f"xTB-ready targets: {int(summary.get('xtb_ready_count', 0))}",
        f"DFT-ready targets: {int(summary.get('dft_ready_count', 0))}",
        f"Governance: {summary.get('governance_status', 'unknown')}",
        "",
        "| Priority | Target | Family | xTB | DFT | Surrogate | Ceiling | Write-back |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("targets", []):
        surrogate = row.get("surrogate", {})
        surrogate_text = surrogate.get("short_label", "-") if surrogate.get("available") else "-"
        lines.append(
            f"| {int(row.get('priority_rank', 0))} | {row.get('id', 'unknown')} | {row.get('slr_family', 'unknown')} | "
            f"{row.get('xtb', {}).get('status', 'unknown')} | {row.get('dft', {}).get('status', 'unknown')} | "
            f"{surrogate_text} | "
            f"{row.get('maximum_claim_level', row.get('promotion_ceiling', 'unknown'))} | "
            f"{', '.join(row.get('write_back_artifacts', [])[:2])} |"
        )
    return "\n".join(lines) + "\n"


def render_computational_gap_dft_ingestion_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Computational Gap DFT Ingestion Report",
        "",
    ]
    if int(summary.get("available_result_count", 0)) == 0:
        lines.append("No DFT results available yet.")
        lines.append("")
    lines.extend(
        [
            f"Manifest jobs: {int(summary.get('manifest_job_count', 0))}",
            f"Ready for DFT: {int(summary.get('ready_for_dft_count', 0))}",
            f"Available results: {int(summary.get('available_result_count', 0))}",
            f"Promotable results: {int(summary.get('promotable_result_count', 0))}",
            "",
            "| Target | Family | Manifest Status | Execution Status | Barrier (kcal/mol) | Promotable | Ceiling | Write-back |",
            "| --- | --- | --- | --- | ---: | --- | --- | --- |",
        ]
    )
    for row in payload.get("jobs", []):
        barrier = row.get("barrier_kcal_mol")
        barrier_text = "" if barrier is None else f"{float(barrier):.2f}"
        lines.append(
            f"| {row.get('target_id', 'unknown')} | {row.get('slr_family', 'unknown')} | {row.get('manifest_status', 'unknown')} | "
            f"{row.get('execution_status', 'not_run')} | {barrier_text} | {'yes' if row.get('promotion_ready') else '-'} | {row.get('promotion_ceiling', 'unknown')} | "
            f"{', '.join(row.get('write_back_artifacts', [])[:2])} |"
        )
    return "\n".join(lines) + "\n"


def render_computational_gap_dft_promotion_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Computational Gap DFT Promotion Report",
        "",
        f"Promotions applied: {int(summary.get('promoted_count', 0))}",
        f"Wrote priors: {'yes' if summary.get('wrote_priors') else 'no'}",
        f"Priors path: {summary.get('priors_path', 'unknown')}",
        "",
        "| Reaction Key | Previous Tier | New Tier | Barrier (kcal/mol) | Ceiling | Method Chain |",
        "| --- | --- | --- | ---: | --- | --- |",
    ]
    for row in payload.get("promotions", []):
        lines.append(
            f"| {row.get('reaction_key', 'unknown')} | {row.get('previous_tier', 'unknown')} | {row.get('new_tier', 'unknown')} | "
            f"{float(row.get('barrier_kcal_mol', 0.0)):.2f} | {row.get('promotion_ceiling', 'unknown')} | {row.get('method_chain', 'unknown')} |"
        )
    if not payload.get("promotions"):
        lines.extend(["", "No promotable DFT results were available."])
    return "\n".join(lines) + "\n"