"""
src/kinetic_core/fit_rows.py -- FIT ROWS AS DATA, KEYED TO THE MEASUREMENT INDEX (2026-10-29, guard 1).

WHY THIS EXISTS
===============
Until now every wave generator typed its own rows: values, pots, pH, laboratory, copied from a dossier's summary. That
is how a pH 5.5 run entered fits as pH 7, a semi-quantitative ladder as calibrated data, one data set two or three
times, and a hold-out's pot with the wrong concentrations, buffer, temperature and time (conflicts register CF-001,
CF-004, CF-005, ...). Here a fit's data are a table in ``data/lit/fit_rows/<name>.yml``:

* ``pots``: each heated system, with its ledger id, the numeric conditions the engine integrates (temperature, time,
  initial mmol/L, buffer, label pH) and the sentence of the paper that states them;
* ``rows``: one per measurement (pot, species, time, value, unit, its error).

What a row may NOT state is read from the ledger entry instead: the laboratory (the leave-one-out group) and the
quantification class (which sets the measurement model's treatment of the series' calibration). A table may name a
method class only where the ledger says ``other`` or ``unknown``, and must then say why (``method_basis``). A row whose
ledger entry is not verified against its paper by two readers is refused (``source_verification_gate`` checks the same
in CI). A ledger entry marked as a copy of another (``same_data_as`` in its flags, or a confirmed duplicate in the
conflicts register) is dropped when its original is present (``fitting.count_once``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from src import data_access, data_paths

from . import fitting

FIT_ROWS_DIR: Path = data_paths.LIT_DIR / "fit_rows"
LEDGER_CLASS_TO_METHOD = {"sida": "sida", "calibrated": "calibrated", "semi_quantitative": "semi_quantitative"}
UNITS = {"mmol_L", "umol_L", "ug_L", "ug_kg", "pH", "mmol_per_L_naoh",
         "percent_of_unheated"}   # a thiol surviving heating as % of the unheated control (Kumazawa 2003): a ratio
ROW_FIELDS = ("pot", "species", "t_min", "value", "unit")
POT_FIELDS = ("ledger_id", "t_c", "minutes", "initial_mmol_L", "condition_quote")


@dataclass(frozen=True)
class FitRow:
    table: str
    pot: str
    ledger_id: str
    laboratory: str
    method: str
    species: str
    t_min: float
    value: float
    unit: str
    sigma: float          # in the residual's units: dex for a concentration, pH units for a pH
    series: str
    raw: Mapping[str, Any] = field(default_factory=dict)

    @property
    def log_scale(self) -> bool:
        return self.unit != "pH"


LABORATORIES: Path = data_paths.KEYS_DIR / "laboratories.yml"


def _laboratory_of() -> Dict[str, str]:
    """ledger lab string -> laboratory group key (data/keys/laboratories.yml)."""
    groups = (data_access.load_yaml(LABORATORIES, missing_ok=True) or {}).get("groups") or {}
    return {s: key for key, g in groups.items() for s in g.get("strings") or []}


def _ledger() -> Dict[str, dict]:
    out = {}
    for path in sorted(data_paths.MEASUREMENT_LEDGER_DIR.glob("*.yml")):
        for entry in (data_access.load_yaml(path) or {}).get("entries") or []:
            out[entry["id"]] = entry
    return out


def _verified() -> Dict[str, bool]:
    checks = ("conditions", "units", "quantification", "duplicates", "values")
    out = {}
    files = [data_paths.SOURCE_VERIFICATION, *sorted(data_paths.SOURCE_VERIFICATION_DIR.glob("*.yml"))]
    entries = [e for f in files for e in (data_access.load_yaml(f, missing_ok=True) or {}).get("entries") or []]
    for entry in entries:
        c = entry.get("checks") or {}
        second = str(entry.get("second_reader") or "").strip()
        out[entry["ledger_id"]] = all(c.get(k) in ("pass", "corrected") for k in checks) and \
            second not in ("", "pending", str(entry.get("read_by") or "").strip())
    return out


def _sigma(row: Mapping[str, Any]) -> float:
    """A row's error in residual units: the reported SD (as dex for a concentration) and the read error of a
    digitised point, added in quadrature. A row with neither is refused: an error is a property of the data."""
    unit, value = row["unit"], float(row["value"])
    parts = []
    if row.get("sd") is not None:
        sd = float(row["sd"])
        parts.append(sd if unit == "pH" else math.log10(math.e) * sd / value if value > 0 else float("inf"))
    if row.get("read_error") is not None:
        parts.append(float(row["read_error"]))
    if not parts:
        raise ValueError(f"row {row.get('pot')}/{row.get('species')}@{row.get('t_min')}: no sd and no read_error")
    return math.sqrt(sum(p * p for p in parts))


def load(name: str, require_verified: bool = True) -> Tuple[Dict[str, dict], Tuple[FitRow, ...]]:
    """(pots, rows) of ``data/lit/fit_rows/<name>.yml``; duplicates of an original that is present are dropped."""
    table = data_access.load_yaml(FIT_ROWS_DIR / f"{name}.yml")
    if require_verified and table.get("status", "draft") != "ready":
        raise PermissionError(f"{name}: the table is {table.get('status', 'draft')!r}, not 'ready': its sources are still being "
                              f"read; set status: ready once every source is verified by two readers")
    ledger, verified, labs = _ledger(), _verified(), _laboratory_of()
    pots: Dict[str, dict] = table.get("pots") or {}
    for pot_id, pot in pots.items():
        missing = [f for f in POT_FIELDS if f not in pot]
        if missing:
            raise ValueError(f"{name}: pot {pot_id} lacks {missing}")
        if pot["ledger_id"] not in ledger:
            raise KeyError(f"{name}: pot {pot_id} cites {pot['ledger_id']!r}, which is in no measurement ledger")
        if require_verified and not verified.get(pot["ledger_id"], False):
            raise PermissionError(f"{name}: pot {pot_id}'s source {pot['ledger_id']} is not verified against its paper "
                                  f"by two readers (data/lit/source_verification.yml)")
    rows: List[FitRow] = []
    for raw in table.get("rows") or []:
        missing = [f for f in ROW_FIELDS if f not in raw]
        if missing:
            raise ValueError(f"{name}: a row lacks {missing}: {raw}")
        if raw["unit"] not in UNITS:
            raise ValueError(f"{name}: unit {raw['unit']!r} not in {sorted(UNITS)}")
        pot = pots[raw["pot"]]
        entry = ledger[pot["ledger_id"]]
        ledger_class = str(entry.get("quant_class"))
        if raw["unit"] == "pH":
            method = "calibrated"     # a glass electrode reads pH on an absolute scale; its error is the row's sigma
        elif ledger_class in LEDGER_CLASS_TO_METHOD:
            if raw.get("method") not in (None, LEDGER_CLASS_TO_METHOD[ledger_class]):
                raise ValueError(f"{name}: row says method {raw['method']!r} but the ledger says {ledger_class!r}")
            method = LEDGER_CLASS_TO_METHOD[ledger_class]
        else:
            if raw.get("method") not in fitting.SCALE_TAU or not str(raw.get("method_basis") or "").strip():
                raise ValueError(f"{name}: the ledger's class is {ledger_class!r}, so the row must name a method class "
                                 f"from {sorted(fitting.SCALE_TAU)} and say why (method_basis)")
            method = raw["method"]
        lab_text = str(entry.get("lab"))
        if lab_text not in labs:
            raise KeyError(f"{name}: laboratory {lab_text!r} of {pot['ledger_id']} is not in data/keys/laboratories.yml; "
                           f"place it in a group (check the paper's authors) before fitting it")
        rows.append(FitRow(table=name, pot=raw["pot"], ledger_id=pot["ledger_id"], laboratory=labs[lab_text],
                           method=method, species=str(raw["species"]), t_min=float(raw["t_min"]),
                           value=float(raw["value"]), unit=raw["unit"], sigma=_sigma(raw),
                           series=f"{pot['ledger_id']}|{raw['pot']}|{raw['species']}", raw=raw))
    return pots, tuple(_drop_copies(rows, ledger))


def _drop_copies(rows: List[FitRow], ledger: Dict[str, dict]) -> List[FitRow]:
    """A row whose ledger entry is flagged a copy (``duplicate_of: <ledger id or paper stem>``, or the older
    ``same_data_as:``) is dropped when its original is among the rows: the original ledger id itself, or any entry
    of the original paper. The ledgers write ``duplicate_of``; reading only ``same_data_as`` (the first version)
    dropped nothing (first reader, sub-model 1, 2026-10-29)."""
    present_ids = {r.ledger_id for r in rows}
    present_papers = {str(ledger[i].get("dossier", "")).replace("_extraction.md", "") for i in present_ids}
    def original_present(entry_id: str) -> bool:
        for flag in ledger[entry_id].get("flags") or []:
            if isinstance(flag, str) and flag.split(":", 1)[0].strip() in ("duplicate_of", "same_data_as") and ":" in flag:
                target = flag.split(":", 1)[1].strip()
                if target != entry_id and (target in present_ids or target in present_papers - {
                        str(ledger[entry_id].get("dossier", "")).replace("_extraction.md", "")}):
                    return True
        return False
    return [r for r in rows if not original_present(r.ledger_id)]


def measurement_model() -> fitting.MeasurementModel:
    """The measurement model for FitRow rows: series and method from the table and the ledger, sigma from the row."""
    return fitting.MeasurementModel(series_of=lambda r: r.series, sigma_of=lambda r: r.sigma, method_of=lambda r: r.method)


def residual(row: FitRow, predicted: float) -> float:
    """The row's residual in its own units (before whitening): log10 ratio for a concentration, difference for pH."""
    if not row.log_scale:
        return float(predicted) - row.value
    return math.log10(max(float(predicted), 1e-30) / row.value) if row.value > 0 else float("nan")
