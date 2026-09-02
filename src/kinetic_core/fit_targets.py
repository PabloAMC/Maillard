"""
src/kinetic_core/fit_targets.py -- WHICH PANEL ROWS THE CORE'S OWN FITS SAW (retirement step B3).

``src/fit_target_index.py`` answers "was this benchmark a fit target?" for the
LEGACY engine: it reads ``*refit*.json`` / ``*rederivation*.json`` records, which
name their targets by benchmark file. The kinetic core's fit reports
(``results/validation/kinetic_core_b*_fit_report.json``) do not: their rows are
keyed by a fit-row id (``hofmann_ribose_FFT``) and a ``source_anchor`` sentence,
because most of them are ratios, conversions and shares that no benchmark bundle
carries. So the index could not see that the SULFUR lane's fit (B2 .. B8, 62 rows,
23 free parameters at B8) read the very Table 1 rows of Hofmann & Schieberle 1998
that the trust-loop panel scores as ``hofmann1998_*_cysteine_145C_20min_pH5`` --
nor that it read the XYLOSE pH-5 row that sits on the maillard_path HOLD-OUT
panel. Found 2026-09-03 while building the core scorecard; declared here.

The declaration is a table, checked by a unit test against the bundles it names
(the compound must be a target of that bundle, and the anchor must quote the
bundle's measured value). Leverage is READ from the B8 report, not typed: 23
free parameters over 62 rows is 0.37 per row, below
``fit_target_index.FIT_LEVERAGE_THRESHOLD`` (0.5), so by the repository's own
rule these rows are ``global_low_leverage``: they STAY in the coverage
denominator, annotated, and are additionally reported as an ``in_core_fit``
split so a reader can see the out-of-sample number on its own.

The other lanes' fits name no panel bundle: B1 (trunk) fitted the Martins 2005
time series, B3 (acrylamide) Claeys / De Vleeschouwer / Knol rate constants, B6
(lipid) Frankel 1989 product shares, B7 (furanic) Blank 1997 yields. None of
those is a scored bundle.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from src import data_access, data_paths
from src.benchmark_metadata import benchmark_signal_origin
from src.fit_target_index import FIT_LEVERAGE_THRESHOLD, fit_target_records_for

SULFUR_FIT_REPORT: Path = data_paths.VALIDATION_DIR / "kinetic_core_b8_fit_report.json"

MFT = "2-Methyl-3-furanthiol (MFT)"
FFT = "2-Furfurylthiol (FFT)"

#: fit-row id (as printed in the B2.x / B8 fit reports' ``rows[].id``) ->
#: (benchmark_id, compound as that bundle names it, the measured value the
#: anchor quotes, its unit). Level rows only: the ratio / share / conversion rows
#: of the same fit constrain nothing a bundle scores.
CORE_SULFUR_FIT_ROWS: Mapping[str, Tuple[str, str, float, str]] = {
    "hofmann_ribose_FFT": ("hofmann1998_ribose_cysteine_145C_20min_pH5", FFT, 121.0, "ppb"),
    "hofmann_ribose_MFT": ("hofmann1998_ribose_cysteine_145C_20min_pH5", MFT, 198.0, "ppb"),
    "hofmann_glucose_FFT": ("hofmann1998_glucose_cysteine_145C_20min_pH5", FFT, 28.0, "ppb"),
    "hofmann_glucose_MFT": ("hofmann1998_glucose_cysteine_145C_20min_pH5", MFT, 19.0, "ppb"),
    "hofmann_fructose_FFT": ("hofmann1998_fructose_cysteine_145C_20min_pH5", FFT, 32.0, "ppb"),
    "hofmann_fructose_MFT": ("hofmann1998_fructose_cysteine_145C_20min_pH5", MFT, 25.0, "ppb"),
    # THE HOLD-OUT ROW THE FIT SAW. Hofmann 1998 Table 1, xylose pH 5.0.
    "hofmann_xylose_FFT": ("mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5", FFT, 96.0, "ppb"),
    "hofmann_xylose_MFT": ("mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5", MFT, 143.0, "ppb"),
    # Fed-intermediate step-level rows (Tables 3, 4, 10): mol % of the fed precursor.
    "fed_furfural_h2s_FFT": ("hofmann1998_furan2aldehyde_h2s_145C_20min_pH5", FFT, 0.48, "mol_percent"),
    "fed_nf_h2s_MFT": ("hofmann1998_norfuraneol_h2s_145C_20min_pH5", MFT, 0.19, "mol_percent"),
    "fed_nf_cys_MFT": ("hofmann1998_norfuraneol_cysteine_145C_20min_pH5", MFT, 0.05, "mol_percent"),
    "fed_c2c3_MFT": ("hofmann1998_c2c3_recombination_145C_20min_pH5", MFT, 0.24, "mol_percent"),
    "fed_c2c3_MFT_pH3": ("hofmann1998_c2c3_recombination_145C_20min_pH3", MFT, None, "mol_percent"),
    "fed_c2c3_MFT_pH7": ("hofmann1998_c2c3_recombination_145C_20min_pH7", MFT, None, "mol_percent"),
}


@dataclass(frozen=True)
class CoreFitTarget:
    record: str            # fit report file name
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


@lru_cache(maxsize=1)
def sulfur_fit_leverage() -> Tuple[Optional[int], Optional[int], str]:
    """``(free_parameters, fitted_rows, leverage_class)`` read from the B8 report."""
    payload = data_access.load_json(SULFUR_FIT_REPORT, missing_ok=True)
    if not payload:
        return None, None, "undeclared"
    free = (payload.get("free_set") or {}).get("n_free")
    rows = max(
        (int(m.get("n_rows", 0)) for m in payload.get("members", []) if isinstance(m, Mapping)),
        default=0,
    )
    if not free or not rows:
        return None, None, "undeclared"
    klass = (
        "per_row_recovery" if (float(free) / float(rows)) >= FIT_LEVERAGE_THRESHOLD
        else "global_low_leverage"
    )
    return int(free), int(rows), klass


@lru_cache(maxsize=1)
def _by_benchmark() -> Dict[str, Tuple[CoreFitTarget, ...]]:
    free, rows, klass = sulfur_fit_leverage()
    grouped: Dict[str, List[Tuple[str, str]]] = {}
    for row_id, (benchmark_id, compound, _value, _unit) in CORE_SULFUR_FIT_ROWS.items():
        grouped.setdefault(benchmark_id, []).append((row_id, compound))
    return {
        benchmark_id: (
            CoreFitTarget(
                record=SULFUR_FIT_REPORT.name, lane="sulfur", leverage_class=klass,
                free_parameters=free, fitted_rows=rows, rows=tuple(pairs),
            ),
        )
        for benchmark_id, pairs in grouped.items()
    }


def core_fit_targets(benchmark_id: str) -> Tuple[CoreFitTarget, ...]:
    """The core fit records that read a level row of this benchmark."""
    return _by_benchmark().get(str(benchmark_id), ())


def in_core_fit(benchmark_id: str, compound: str) -> bool:
    """True when THIS (benchmark, compound) level row was a row of a core fit."""
    return any(
        c == compound for target in core_fit_targets(benchmark_id) for _, c in target.rows
    )


def core_evidence_role(benchmark_id: str, bench_file: Path) -> str:
    """
    What kind of claim a CORE score on this benchmark supports.

    ``fit_recovery`` only when a core fit could reproduce it row by row
    (``per_row_recovery``); ``internal_synthetic`` for a non-literature signal;
    else ``predictive`` -- annotated by :func:`core_fit_targets` when a
    low-leverage fit saw it, exactly as ``fit_target_index`` treats the legacy
    global fits. The legacy engine's own fit-recovery declarations (matrix
    observability factors solved from Pratap-Singh 2021, the projection budget)
    say nothing about the core and are NOT consulted here; the scorecard prints
    them beside this role as ``legacy_evidence_role`` so the two can be compared.
    """
    if any(t.leverage_class in {"per_row_recovery", "undeclared"} for t in core_fit_targets(benchmark_id)):
        return "fit_recovery"
    if benchmark_signal_origin(Path(bench_file)) == "internal_synthetic":
        return "internal_synthetic"
    return "predictive"


def fit_target_of(benchmark_id: str) -> Dict[str, List[str]]:
    """``{"core": [...core records...], "legacy": [...legacy index records...]}``."""
    return {
        "core": [t.record for t in core_fit_targets(benchmark_id)],
        "legacy": list(fit_target_records_for(benchmark_id)),
    }


__all__ = [
    "CORE_SULFUR_FIT_ROWS",
    "CoreFitTarget",
    "SULFUR_FIT_REPORT",
    "core_evidence_role",
    "core_fit_targets",
    "fit_target_of",
    "in_core_fit",
    "sulfur_fit_leverage",
]
