"""The decisive-measurements record: which measurement would settle each structural question the literature on disk
cannot, loaded and checked.

2026-11-05. The data wishlist (src/data_wishlist.py) derives "what to measure next" from fit artifacts, so it only
sees fitted coordinates the data do not pin. Several of the plan's questions are structural instead (which rival
mechanism holds), and some have turned out undecidable on disk. ``docs/reference/decisive_measurements.yml`` is the
hand-kept record of those, one entry per question; this module loads it, refuses a malformed entry, and checks that
every record an entry cites exists: a CF- / EX- / CR- id in its register, anything else as a repository path. A record
that lives on an unmerged branch is accepted only when the file's ``pending_records`` names where it lives; once it
lands, the mark is reported as stale (a warning, never a failure, so a register merging first breaks nothing).

Nothing here derives a measurement: the record is curated, and the wishlist says so where it renders it.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from src import data_access, data_paths

RECORD_PATH = data_paths.DECISIVE_MEASUREMENTS
REGISTER_PATHS: Dict[str, Path] = {
    "CF": data_paths.CONFLICTS_REGISTER,
    "EX": data_paths.EXPERIMENT_REGISTER,
    "CR": data_paths.CHOICES_REGISTER,
}
STATUSES = ("open", "provisional", "decided")
PH_BASIS_TAGS = ("cooled", "label_25C", "in_situ", "not_stated")
NOT_WRITTEN = "design not yet written"
REQUIRED_FIELDS = ("id", "plan_step", "added", "question", "rivals_from", "rivals", "why_not_on_disk", "measurement",
                   "decides", "unlocks", "existing_data_route", "status")
WRITTEN_ARM_FIELDS = ("conditions", "species", "ph_values", "ph_basis", "temperatures_C", "o2", "replicates")
REGISTER_ID = re.compile(r"^(CF|EX|CR)-\d{3}$")
DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
#: The laptop-only PDF folder is private (the repository is public): a record cites the dossier, never the file.
_PRIVATE_PREFIX = data_paths.rel(data_paths.ARTICLES_DIR).rstrip("/") + "/"


class RecordError(ValueError):
    """The decisive-measurements record is malformed or cites a record that does not exist."""


def register_ids(paths: Mapping[str, Path] = REGISTER_PATHS) -> Dict[str, Set[str]]:
    """Every id in each register, keyed by prefix (CF, EX, CR)."""
    out: Dict[str, Set[str]] = {}
    for prefix, path in paths.items():
        payload = data_access.load_yaml(path)
        out[prefix] = {str(e["id"]) for e in (payload or {}).get("entries", []) if isinstance(e, dict) and "id" in e}
    return out


def _nonempty_str(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _check_ref(ref: Any, where: str, ids: Mapping[str, Set[str]], exists: Callable[[str], bool],
               pending: Mapping[str, str], problems: List[str], used_pending: Set[str]) -> None:
    if not _nonempty_str(ref):
        problems.append(f"{where}: a citation must be a non-empty string, got {ref!r}")
        return
    ref = ref.strip()
    if ref.startswith(_PRIVATE_PREFIX):
        problems.append(f"{where}: {ref} is a private, laptop-only file; cite its dossier instead")
        return
    m = REGISTER_ID.match(ref)
    present = (ref in ids.get(m.group(1), set())) if m else exists(ref)
    if present:
        return
    if ref in pending:
        used_pending.add(ref)
        return
    kind = "register id" if m else "path"
    problems.append(f"{where}: cites {ref}, a {kind} that is neither on this branch nor in pending_records")


def _check_arm(arm: Any, where: str, problems: List[str]) -> None:
    if not isinstance(arm, dict):
        problems.append(f"{where}: an arm must be a mapping")
        return
    if not _nonempty_str(arm.get("name")):
        problems.append(f"{where}: arm has no name")
    written = arm.get("design_written")
    if not isinstance(written, bool):
        problems.append(f"{where}: design_written must be true or false")
        return
    if not written:
        req = arm.get("requirement")
        if not (_nonempty_str(req) and req.strip().startswith(NOT_WRITTEN)):
            problems.append(f"{where}: an arm with no written design must say '{NOT_WRITTEN}' in its requirement")
        return
    for field in WRITTEN_ARM_FIELDS:
        value = arm.get(field)
        if field in ("species", "ph_basis", "temperatures_C"):
            if not (isinstance(value, list) and value):
                problems.append(f"{where}: {field} must be a non-empty list")
        elif not _nonempty_str(value):
            problems.append(f"{where}: {field} is missing (write 'not stated in the records' rather than leave it out)")
    for tag in arm.get("ph_basis") or []:
        if tag not in PH_BASIS_TAGS:
            problems.append(f"{where}: pH basis {tag!r} is not one of {', '.join(PH_BASIS_TAGS)}")


def validate(payload: Any, *, ids: Mapping[str, Set[str]], exists: Callable[[str], bool]) -> Tuple[List[str], List[str]]:
    """``(problems, warnings)`` for a parsed record. ``ids`` is :func:`register_ids`; ``exists`` answers whether a
    repository-relative path is on this branch. Problems make the record unusable; warnings name pending marks whose
    record has since landed, or that nothing cites."""
    problems: List[str] = []
    if not isinstance(payload, dict) or not isinstance(payload.get("entries"), list) or not payload["entries"]:
        return ["the record must be a mapping with a non-empty 'entries' list"], []
    pending = payload.get("pending_records") or {}
    if not isinstance(pending, dict) or not all(_nonempty_str(v) for v in pending.values()):
        problems.append("pending_records must map each record to where it lives (a non-empty string)")
        pending = {}
    used_pending: Set[str] = set()
    seen: Set[str] = set()

    def ref(r: Any, where: str) -> None:
        _check_ref(r, where, ids, exists, pending, problems, used_pending)

    for i, entry in enumerate(payload["entries"]):
        where = f"entries[{i}]"
        if not isinstance(entry, dict):
            problems.append(f"{where}: an entry must be a mapping")
            continue
        eid = entry.get("id")
        where = f"entry {eid or i}"
        missing = [f for f in REQUIRED_FIELDS if f not in entry]
        if missing:
            problems.append(f"{where}: missing {', '.join(missing)}")
            continue
        if eid in seen:
            problems.append(f"{where}: duplicate id")
        seen.add(eid)
        for field in ("id", "plan_step", "question"):
            if not _nonempty_str(entry[field]):
                problems.append(f"{where}: {field} must be a non-empty string")
        if not (isinstance(entry["added"], str) and DATE.match(entry["added"])):
            problems.append(f"{where}: added must be a YYYY-MM-DD string")
        if not (isinstance(entry["rivals"], list) and len(entry["rivals"]) >= 2 and all(_nonempty_str(r) for r in entry["rivals"])):
            problems.append(f"{where}: rivals must name at least two rivals")
        for r in entry["rivals_from"] if isinstance(entry["rivals_from"], list) else [None]:
            ref(r, f"{where}.rivals_from")
        blockers = entry["why_not_on_disk"]
        if not (isinstance(blockers, list) and blockers):
            problems.append(f"{where}: why_not_on_disk must list at least one blocker")
            blockers = []
        for j, b in enumerate(blockers):
            if not (isinstance(b, dict) and _nonempty_str(b.get("blocker")) and isinstance(b.get("cites"), list) and b["cites"]):
                problems.append(f"{where}.why_not_on_disk[{j}]: needs a blocker sentence and at least one citation")
                continue
            for r in b["cites"]:
                ref(r, f"{where}.why_not_on_disk[{j}]")
        m = entry["measurement"]
        if not (isinstance(m, dict) and _nonempty_str(m.get("source")) and _nonempty_str(m.get("summary"))
                and isinstance(m.get("arms"), list) and m["arms"]):
            problems.append(f"{where}.measurement: needs a source, a summary and at least one arm")
        else:
            for k, arm in enumerate(m["arms"]):
                _check_arm(arm, f"{where}.measurement.arms[{k}]", problems)
        decides = entry["decides"]
        if not (isinstance(decides, list) and decides and all(
                isinstance(d, dict) and _nonempty_str(d.get("outcome")) and _nonempty_str(d.get("means")) for d in decides)):
            problems.append(f"{where}: decides must list outcomes, each with what it means for the rivals")
        if not (isinstance(entry["unlocks"], list) and entry["unlocks"] and all(_nonempty_str(u) for u in entry["unlocks"])):
            problems.append(f"{where}: unlocks must list what the measurement could lift")
        if entry["existing_data_route"] is not None and not _nonempty_str(entry["existing_data_route"]):
            problems.append(f"{where}: existing_data_route must be a sentence or null")
        status = entry["status"]
        if status not in STATUSES:
            problems.append(f"{where}: status {status!r} is not one of {', '.join(STATUSES)}")
        if status == "provisional" and not _nonempty_str(entry.get("status_note")):
            problems.append(f"{where}: a provisional entry needs a status_note saying what it waits for")
        if status == "decided" and not _nonempty_str(entry.get("decided_by")):
            problems.append(f"{where}: a decided entry needs decided_by")
        for j, d in enumerate(entry.get("disagreements") or []):
            if not (isinstance(d, dict) and _nonempty_str(d.get("about")) and _nonempty_str(d.get("note"))
                    and isinstance(d.get("cites"), list) and d["cites"]):
                problems.append(f"{where}.disagreements[{j}]: needs about, note and at least one citation")
                continue
            for r in d["cites"]:
                ref(r, f"{where}.disagreements[{j}]")

    warnings: List[str] = []
    for key, where_it_lives in pending.items():
        m = REGISTER_ID.match(str(key))
        landed = (key in ids.get(m.group(1), set())) if m else exists(str(key))
        if landed:
            warnings.append(f"pending_records: {key} is now on this branch (was: {where_it_lives}); drop the mark")
        elif key not in used_pending:
            warnings.append(f"pending_records: {key} is cited by no entry")
    return problems, warnings


def _repo_path_exists(relative: str) -> bool:
    return (data_paths.REPO_ROOT / relative).exists()


def load(path: Path = RECORD_PATH, *, ids: Optional[Mapping[str, Set[str]]] = None,
         exists: Callable[[str], bool] = _repo_path_exists) -> Tuple[List[Dict[str, Any]], List[str]]:
    """The validated entries and the warnings. Raises :class:`RecordError` listing every problem."""
    payload = data_access.load_yaml(path)
    problems, warnings = validate(payload, ids=register_ids() if ids is None else ids, exists=exists)
    if problems:
        raise RecordError(f"{data_paths.rel(path)}: " + "; ".join(problems))
    return list(payload["entries"]), warnings


def _arm_text(arm: Mapping[str, Any]) -> str:
    if not arm.get("design_written"):
        return f"{arm['name']}: {arm['requirement'].strip()}"
    temps = ", ".join(f"{t} C" if isinstance(t, (int, float)) else str(t) for t in arm["temperatures_C"])
    return (f"{arm['name']}: {arm['conditions'].strip().rstrip('.')}; pH {arm['ph_values']} ({', '.join(arm['ph_basis'])}); "
            f"temperatures {temps}; O2 {arm['o2']}; replicates {arm['replicates']}")


def summarise(entries: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """One row per entry for the wishlist's JSON: the question, the step, the measurement, what it unlocks, the status."""
    rows = []
    for e in entries:
        m = e["measurement"]
        rows.append({
            "id": e["id"], "plan_step": e["plan_step"], "question": " ".join(e["question"].split()),
            "rivals": list(e["rivals"]),
            "measurement": " ".join(m["summary"].split()), "measurement_source": m["source"],
            "arms": [" ".join(_arm_text(a).split()) for a in m["arms"]],
            "design_written": all(a.get("design_written") for a in m["arms"]),
            "decides": [{"outcome": d["outcome"], "means": d["means"]} for d in e["decides"]],
            "unlocks": list(e["unlocks"]), "existing_data_route": e["existing_data_route"],
            "status": e["status"], "status_note": " ".join(str(e.get("status_note") or "").split()) or None,
            "disagreements": [{"about": d["about"], "note": " ".join(d["note"].split())} for d in e.get("disagreements") or []],
            "added": e["added"],
        })
    return rows


__all__ = ["RECORD_PATH", "RecordError", "load", "validate", "register_ids", "summarise", "PH_BASIS_TAGS", "STATUSES"]
