#!/usr/bin/env python
"""
pathway_atlas_gate -- the pathway-condition atlas is well-formed and every cell is grep-able in its sources (2026-11-01).

Checks every `docs/reference/pathway_atlas/F*.yml` against the shared vocabulary in `families.yml`:
  - `family` matches the filename and is an id in families.yml; the file has a non-empty `markers` list and `cells`;
  - every cell's `key` is exactly {ph, temperature, water}, each value a bin of families.yml (`not_stated` allowed);
  - the temperature key AGREES with `tags.T_C`: recomputed with low <= 60 < moderate <= 100 < high <= 140 < very_high;
    a range or list ("80-120", "100 / 130") may sit in any bin it spans; text in parentheses or after ';' is detail;
  - the water key agrees with `tags.buffer_M`: a buffer of >= 0.5 M cannot be keyed `dilute` or `not_stated`;
  - the tag axes (sugar, amine, oxygen, buffer, metals) hold a bin, or a list of bins when the cell spans arms of a
    contrast or paired pots (detail in `<axis>_note`); other tags (ph_basis, T_C, ...) are free-form, but their
    names must be plain identifiers;
  - `status` and `evidence` are from the vocabulary, and only the known cell fields appear;
  - `absent` carries a non-empty `detection_basis` (a detection limit, a positive control, or a label pattern that
    would have shown the route); a non-detect without one is `not_detected`;
  - the CAP: `dominant` only from labelling / fed_intermediate, or trapping with a mass-balance note on the cell;
  - every cited dossier exists in data/lit/extraction_dossiers/ with EXACT case (CI's Linux disk is case-sensitive);
  - every `quote` is <= 15 words and appears verbatim, after whitespace normalisation, in at least one cited dossier;
  - every CF id in `conflict` exists in docs/reference/conflicts_register.yml; a conflict not yet registered is
    written `unregistered: ...`; a `conflicting` cell names at least one registered CF id.
Structural and offline: it does not judge whether a cell's reading of its paper is right.

    python scripts/ci/pathway_atlas_gate.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "docs" / "reference" / "pathway_atlas"
FAMILIES = ATLAS / "families.yml"
DOSSIERS = ROOT / "data" / "lit" / "extraction_dossiers"
REGISTER = ROOT / "docs" / "reference" / "conflicts_register.yml"

KEY_AXES = ("ph", "temperature", "water")
TAG_AXES = ("sugar", "amine", "oxygen", "buffer", "metals")
STATUSES = {"dominant", "operating", "minor", "absent", "not_detected", "conflicting", "unattributed"}
EVIDENCE = {"labelling", "fed_intermediate", "trapping", "measured_intermediate", "mass_balance",
            "within_study_factor", "kinetic_fit", "qualitative"}
DOMINANT_EVIDENCE = {"labelling", "fed_intermediate"}
REQUIRED = ("key", "tags", "status", "evidence", "quote", "lab", "dossiers")
OPTIONAL = ("quantity", "where", "contrast", "note", "conflict", "detection_basis")
MAX_QUOTE_WORDS = 15
CF_ID = re.compile(r"CF-\d{3}")
MASS_BALANCE = re.compile(r"mass[- ]balance", re.IGNORECASE)
IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
T_ORDER = ("low", "moderate", "high", "very_high")
# a temperature number: not a time, pressure or count ("10 h", "25-40 MPa", "24-week")
T_NUMBER = re.compile(r"\d+(?:\.\d+)?(?!\s*(?:h\b|min\b|s\b|MPa|week|\d|\.\d))")
ROOM_T = re.compile(r"room temp|\bRT\b|ambient", re.IGNORECASE)


def _norm(text: str) -> str:
    return " ".join(str(text).split())


def _states_mass_balance(text: str) -> bool:
    """A mass balance is mentioned and not negated ('no mass balance on other fates' does not count)."""
    for match in MASS_BALANCE.finditer(text):
        before = text[max(0, match.start() - 12):match.start()].lower()
        if not re.search(r"\b(no|not|without)\b", before):
            return True
    return False


def _t_bin(t: float) -> str:
    return "low" if t <= 60 else "moderate" if t <= 100 else "high" if t <= 140 else "very_high"


def temperature_bins(t_c) -> set | None:
    """The temperature bins a `T_C` tag allows, or None when it states no temperature (null, 'not stated')."""
    if t_c is None or isinstance(t_c, bool):
        return None
    if isinstance(t_c, (int, float)):
        return {_t_bin(float(t_c))}
    main = re.split(r"[(;]", str(t_c))[0]
    numbers = [float(x) for x in T_NUMBER.findall(main)]
    if not numbers:
        return {"low"} if ROOM_T.search(main) else None
    lo, hi = T_ORDER.index(_t_bin(min(numbers))), T_ORDER.index(_t_bin(max(numbers)))
    return set(T_ORDER[lo:hi + 1])


def _vocabulary() -> tuple:
    spec = yaml.safe_load(FAMILIES.read_text(encoding="utf-8"))
    bins = {axis: set(values) for axis, values in spec["condition_bins"].items()}
    for axis in bins:
        bins[axis].add("not_stated")
    return {f["id"] for f in spec["families"]}, bins


def _register_ids() -> set:
    entries = (yaml.safe_load(REGISTER.read_text(encoding="utf-8")) or {}).get("entries") or []
    return {e["id"] for e in entries}


class _Dossiers:
    def __init__(self) -> None:
        self.names = {p.name for p in DOSSIERS.iterdir()}
        self._text = {}

    def text(self, name: str) -> str:
        if name not in self._text:
            self._text[name] = _norm((DOSSIERS / name).read_text(encoding="utf-8"))
        return self._text[name]


def _check_cell(where: str, cell: dict, bins: dict, register: set, dossiers: _Dossiers) -> list:
    problems = []
    missing = [f for f in REQUIRED if f not in cell]
    if missing:
        return [f"{where}: missing {missing}"]
    extra = sorted(set(cell) - set(REQUIRED) - set(OPTIONAL))
    if extra:
        problems.append(f"{where}: unknown cell field(s) {extra} (put them in tags or note)")

    key = cell["key"]
    if not isinstance(key, dict) or set(key) != set(KEY_AXES):
        problems.append(f"{where}: key must have exactly {list(KEY_AXES)}, got {key!r}")
    else:
        for axis in KEY_AXES:
            if key[axis] not in bins[axis]:
                problems.append(f"{where}: key.{axis} {key[axis]!r} not in {sorted(bins[axis])}")

    tags = cell["tags"] or {}
    if not isinstance(tags, dict):
        problems.append(f"{where}: tags is not a mapping")
        tags = {}
    for name, value in tags.items():
        if not IDENTIFIER.match(str(name)):
            problems.append(f"{where}: tag name {name!r} is not an identifier (an unquoted comma in a flow map?)")
        if name in TAG_AXES:
            values = value if isinstance(value, list) else [value]
            bad = [v for v in values if v not in bins[name]]
            if bad or not values:
                problems.append(f"{where}: tags.{name} {value!r} not a bin (or list of bins) of {sorted(bins[name])}")

    if isinstance(key, dict) and isinstance(tags, dict):
        allowed = temperature_bins(tags.get("T_C"))
        if allowed is not None and key.get("temperature") not in allowed:
            problems.append(f"{where}: key.temperature {key.get('temperature')!r} disagrees with T_C {tags.get('T_C')!r}"
                            f" (allowed: {sorted(allowed, key=T_ORDER.index)})")
        buffer_m = tags.get("buffer_M")
        if isinstance(buffer_m, (int, float)) and not isinstance(buffer_m, bool) and buffer_m >= 0.5 \
                and key.get("water") in ("dilute", "not_stated"):
            problems.append(f"{where}: key.water {key.get('water')!r} but buffer_M {buffer_m} (>= 0.5 M is concentrated)")

    status, evidence = cell["status"], cell["evidence"]
    if status not in STATUSES:
        problems.append(f"{where}: status {status!r} not in {sorted(STATUSES)}")
    if evidence not in EVIDENCE:
        problems.append(f"{where}: evidence {evidence!r} not in {sorted(EVIDENCE)}")
    if status == "dominant" and evidence not in DOMINANT_EVIDENCE:
        noted = " ".join(str(cell.get(f) or "") for f in ("note", "quantity", "contrast"))
        if not (evidence == "trapping" and _states_mass_balance(noted)):
            problems.append(f"{where}: CAP: `dominant` from evidence {evidence!r} (needs labelling / fed_intermediate,"
                            " or trapping with a mass-balance note)")

    if status == "absent" and not str(cell.get("detection_basis") or "").strip():
        problems.append(f"{where}: `absent` without `detection_basis` (a detection limit, a positive control or a"
                        " label pattern that would have shown the route); otherwise the status is `not_detected`")

    cited = cell["dossiers"] if isinstance(cell["dossiers"], list) else [cell["dossiers"]]
    present = []
    for name in cited:
        if name in dossiers.names:
            present.append(name)
        else:
            problems.append(f"{where}: dossier {name!r} not in data/lit/extraction_dossiers/ (names are case-sensitive)")
    quote = _norm(cell["quote"] or "")
    if not quote:
        problems.append(f"{where}: empty quote")
    else:
        if len(quote.split()) > MAX_QUOTE_WORDS:
            problems.append(f"{where}: quote has {len(quote.split())} words (> {MAX_QUOTE_WORDS})")
        if present and not any(quote in dossiers.text(name) for name in present):
            problems.append(f"{where}: quote not found verbatim in {present}: {quote!r}")

    conflict = cell.get("conflict")
    registered = []
    if conflict is not None:
        for item in conflict if isinstance(conflict, list) else [conflict]:
            ids = CF_ID.findall(str(item))
            if not ids and not str(item).startswith("unregistered:"):
                problems.append(f"{where}: conflict {item!r} names no CF id (write 'unregistered: ...' if none exists)")
            for cf in ids:
                if cf in register:
                    registered.append(cf)
                else:
                    problems.append(f"{where}: conflict {cf} not in docs/reference/conflicts_register.yml")
    if status == "conflicting" and not registered:
        problems.append(f"{where}: status `conflicting` but no registered CF id in `conflict`")
    return problems


def check() -> tuple:
    family_ids, bins = _vocabulary()
    register = _register_ids()
    dossiers = _Dossiers()
    problems, n_cells = [], 0
    paths = sorted(ATLAS.glob("F*.yml"))
    if not paths:
        problems.append(f"no F*.yml files in {ATLAS.relative_to(ROOT)}")
    for path in paths:
        try:
            doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError as exc:
            problems.append(f"{path.name}: YAML error: {exc}")
            continue
        family = doc.get("family")
        if family != path.stem:
            problems.append(f"{path.name}: family {family!r} does not match the filename")
        if family not in family_ids:
            problems.append(f"{path.name}: family {family!r} is not an id in families.yml")
        if not doc.get("markers"):
            problems.append(f"{path.name}: no `markers`")
        cells = doc.get("cells") or []
        if not cells:
            problems.append(f"{path.name}: no cells")
        for i, cell in enumerate(cells):
            n_cells += 1
            problems += _check_cell(f"{path.name} cell {i}", cell, bins, register, dossiers)
    return problems, len(paths), n_cells


def main() -> int:
    problems, n_files, n_cells = check()
    if problems:
        print(f"pathway_atlas_gate: FAIL ({len(problems)} problems, {n_files} files, {n_cells} cells)")
        for p in problems:
            print("  - " + p)
        return 1
    print(f"pathway_atlas_gate: PASS ({n_files} files, {n_cells} cells, docs/reference/pathway_atlas/)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
