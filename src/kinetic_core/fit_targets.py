"""
src/kinetic_core/fit_targets.py -- WHICH PANEL ROWS THE CORE'S OWN FITS READ.

``src/fit_target_index.py`` answers "was this benchmark a fit target?" from the fit
records under ``results/validation/``. The kinetic core's fit REPORTS key their rows
by fit-row id and a prose ``source_anchor``, so until 2026-09-03 this module carried
a hand-typed table mapping the sulfur fit's level rows to bundles. That table is
gone: the fit generators' row tables now carry ``benchmark_id`` /
``benchmark_compound`` next to each anchor, and
``scripts/generators/generate_core_fit_targets.py`` writes them as
``results/validation/kinetic_core_<lane>_fit_targets.json`` in the index's own
vocabulary (``fit_target_ids`` + ``fit_leverage``), which the index, the hold-out
guard and this module all read. One declaration, three consumers.

Leverage is READ from the record: 23 free parameters over 62 rows is 0.37 per row,
below ``fit_target_index.FIT_LEVERAGE_THRESHOLD`` (0.5), so by the repository's own
rule these rows are ``global_low_leverage`` -- they STAY in the coverage
denominator, annotated ``in_core_fit``, and are also reported as their own split.

The other lanes' fits name no panel bundle: B1 (trunk) fitted the Martins 2005
time series, B3 (acrylamide) Claeys / De Vleeschouwer / Knol rate constants, B6
(lipid) Frankel 1989 product shares, B7 (furanic) Blank 1997 yields.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from src import data_access, data_paths
from src.benchmark_metadata import benchmark_signal_origin
from src.fit_target_index import fit_leverage_class, fit_target_records_for

RECORD_GLOB = "kinetic_core_*_fit_targets.json"


@dataclass(frozen=True)
class CoreFitTarget:
    record: str            # fit-target record file name
    lane: str
    leverage_class: str    # per_row_recovery | global_low_leverage | undeclared
    free_parameters: Optional[int]
    fitted_rows: Optional[int]
    rows: Tuple[Tuple[str, str], ...]  # (fit-row id, compound) of THIS benchmark

    def as_dict(self) -> Dict[str, Any]:
        return {
            "record": self.record,
            "lane": self.lane,
            "leverage_class": self.leverage_class,
            "free_parameters": self.free_parameters,
            "fitted_rows": self.fitted_rows,
            "rows": [list(r) for r in self.rows],
        }


def record_paths() -> List[Path]:
    return sorted(data_paths.VALIDATION_DIR.glob(RECORD_GLOB))


@lru_cache(maxsize=1)
def _by_benchmark() -> Dict[str, Tuple[CoreFitTarget, ...]]:
    out: Dict[str, List[CoreFitTarget]] = {}
    for path in record_paths():
        payload = data_access.load_json(path)
        if payload.get("artifact") != "kinetic_core_fit_targets":
            continue
        leverage = dict(payload.get("fit_leverage") or {})
        klass = fit_leverage_class({"leverage": leverage})
        grouped: Dict[str, List[Tuple[str, str]]] = {}
        for row in payload.get("rows", []):
            grouped.setdefault(str(row["benchmark_id"]), []).append((str(row["id"]), str(row["compound"])))
        for benchmark_id, pairs in grouped.items():
            out.setdefault(benchmark_id, []).append(
                CoreFitTarget(
                    record=path.name, lane=str(payload.get("lane", "")), leverage_class=klass,
                    free_parameters=leverage.get("free_parameters"), fitted_rows=leverage.get("fitted_rows"),
                    rows=tuple(pairs),
                )
            )
    return {k: tuple(v) for k, v in out.items()}


def core_fit_targets(benchmark_id: str) -> Tuple[CoreFitTarget, ...]:
    """The core fit records that read a level row of this benchmark."""
    return _by_benchmark().get(str(benchmark_id), ())


def in_core_fit(benchmark_id: str, compound: str) -> bool:
    """True when THIS (benchmark, compound) level row was a row of a core fit."""
    return any(c == compound for target in core_fit_targets(benchmark_id) for _, c in target.rows)


#: The ONE evidence-role vocabulary (2026-09-03). What kind of claim a core score supports:
#:   fit_recovery       a core fit could reproduce the row (per-row leverage): algebra, not evidence
#:   internal_synthetic the comparator is frozen model output, not a measurement
#:   external_holdout   a measurement physically separated from every fit (external_validation/**)
#:   predictive         an external measurement on the trust loop the core was not fitted to
#: The bundle-level ``evidence_class`` label (``external_validation_only`` / ``diagnostic_only``)
#: is the DATA-SIDE marker the gates check; this function turns it and the fit records into the role.
EVIDENCE_ROLES = ("fit_recovery", "internal_synthetic", "external_holdout", "predictive")


def core_evidence_role(benchmark_id: str, bench_file: Path) -> str:
    """One of :data:`EVIDENCE_ROLES` for a bundle, from its fit records, signal origin and location."""
    path = Path(bench_file)
    if any(t.leverage_class in {"per_row_recovery", "undeclared"} for t in core_fit_targets(benchmark_id)):
        return "fit_recovery"
    if benchmark_signal_origin(path) == "internal_synthetic":
        return "internal_synthetic"
    try:
        if data_paths.EXTERNAL_VALIDATION_DIR in path.resolve().parents:
            return "external_holdout"
    except OSError:
        pass
    return "predictive"


def fit_target_of(benchmark_id: str) -> Dict[str, List[str]]:
    """``{"core": [...core fit-target records...], "index": [...every index record...]}``."""
    return {
        "core": [t.record for t in core_fit_targets(benchmark_id)],
        "index": list(fit_target_records_for(benchmark_id)),
    }


def clear_cache() -> None:
    _by_benchmark.cache_clear()


__all__ = [
    "EVIDENCE_ROLES",
    "CoreFitTarget",
    "RECORD_GLOB",
    "clear_cache",
    "core_evidence_role",
    "core_fit_targets",
    "fit_target_of",
    "in_core_fit",
    "record_paths",
]
