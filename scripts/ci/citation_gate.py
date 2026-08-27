#!/usr/bin/env python
"""Citation-contamination regression gate.

Context
-------
The 2026-08-26 citation sweep (``results/validation/citation_verification_ledger.md``)
found that **137 of 225 unique DOI anchors (61%) were defective**: 45 DEAD, 45
TOPIC-MISMATCH (live DOI, wrong paper), 47 METADATA-MISMATCH. Several of the dead
ones were not DOIs at all but fabricated strings with a recognisable confabulation
signature -- ``10.1021/acs.jafc.de_leyn_2019``, ``10.1021/acs.jafc.liardon_1991``,
``10.1016/j.foodchem.2015.00000``, ``10.1016/j.foodres.2025.001279``.

This gate exists so that class of damage cannot silently regrow. It is deliberately
**structural and offline**: it runs in seconds, needs no network, and cannot flake.
It does *not* try to replace the human/CrossRef verification pass -- a DOI that
resolves to the wrong paper is invisible to any structural check. Use ``--online``
(non-blocking, see ``ci.yml``) for liveness sampling.

Scan surfaces
-------------
* **Structured** (``SCAN_GLOBS``): fields NAMED for a DOI inside ``data/lit/*.json|yml``,
  ``data/benchmarks/**/*.json`` and ``data/species/*.yml``. All four checks apply.
* **Text** (``TEXT_SCAN_GLOBS``, added 2026-08-27, Wave I fix 4): every DOI-shaped
  token anywhere in ``src/**/*.py``, ``docs/**/*.md``, ``README.md`` and ``AUDIT.md``.
  Only checks 1 and 2 apply -- 3 and 4 are properties of a structured record.
  This surface exists because the flagship mechanism anchor Wave I had to repair
  (``10.1021/jf60200a038``, live but resolving to a gossypol/rat paper) was in
  Python comments, a docstring and a free-text calibration rationale, and the
  structured sweep could not see any of them.

Checks (offline, all blocking)
------------------------------
1. **DOI syntax** -- every DOI-bearing field, and every DOI-shaped token on the
   text surface, matches the DOI grammar.
2. **Confabulation signatures** -- author-name-shaped tails (``...de_leyn_2019``),
   all-zero article suffixes (``.00000``), and template placeholders
   (``XXXX``, ``TODO``, ``PLACEHOLDER``, ``10.xxxx/...``) are rejected outright,
   and are *not* waivable.
3. **Status coherence** -- no record may carry both a ``doi`` and
   ``source_status: no_verifiable_source``. Either the DOI verifies the source or
   there is no verifiable source; asserting both means one of the two is stale.
4. **Repair-record completeness** -- every ``doi_repair`` / ``citation_repair``
   record must carry ``old``, ``new``, ``date`` and ``basis``. A repair without a
   basis is indistinguishable from a fresh confabulation.
5. **Barrier-source disclosure** (``BARRIER_SOURCE_GLOBS``, added 2026-08-27,
   Wave J1) -- every entry in the ``data/qm/`` barrier benchmark files must name
   a source or carry ``source_status: no_verifiable_source``. This one checks for
   an ABSENT citation rather than a bad one: all nine literature windows and all
   eighteen claimed wB97M-V / revDSD / xTB values in those files shipped with no
   citation, no run record and no ledger entry, and silence was reading as
   provenance. The label satisfies the check by design -- it enforces disclosure,
   not quality -- and check 3 stops an entry from claiming the label and a DOI
   at once. Not waivable, because there is nothing to waive: adding the label is
   always possible and always free.

Known-debt baseline
-------------------
Checks 1 and 3 have a small explicit waiver list (``WAIVERS`` below) recording the
violations that already existed when this gate was written. Waivers are keyed by
exact (file, json-pointer, value), so any *new* violation fails the build; and a
waiver that stops matching is itself reported as stale, so the baseline can only
ratchet down. Check 2 has no waivers by design.

Usage
-----
    python scripts/ci/citation_gate.py               # offline, blocking (CI)
    python scripts/ci/citation_gate.py --online -n 25  # + liveness sample (advisory)
"""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
from pathlib import Path
from typing import Any, Iterator

ROOT = Path(__file__).resolve().parents[2]

# Globs swept for citation anchors. Kept in sync with the ledger's 8 file globs.
SCAN_GLOBS = (
    "data/lit/*.json",
    "data/lit/*.yml",
    "data/lit/*.yaml",
    "data/benchmarks/**/*.json",
    "data/species/*.yml",
    "data/species/*.yaml",
    # 2026-08-27 (Wave J1). data/qm/ was outside every sweep this repo has ever
    # run -- including the 2026-08-26 CrossRef ledger -- because `.gitignore`
    # line 33 (`data/*`) hid the directory from git entirely, so it was neither
    # reviewable nor scannable. Un-ignoring it (same commit) is what makes it
    # reachable here. See BARRIER_SOURCE_GLOBS below for the check it feeds.
    "data/qm/*.json",
)

# 2026-08-27 (Wave J1). Check 5's scope. These two files carry nine "literature"
# barrier windows and eighteen columns of claimed quantum-chemistry results, and
# between them they carried ZERO citations, run records or provenance of any kind.
# The check below makes that state impossible to re-enter silently: an entry must
# either name a source or admit that it has none.
BARRIER_SOURCE_GLOBS = (
    "data/qm/phase33_barrier_benchmarks.json",
    "data/qm/phase35_double_hybrid_benchmarks.json",
)

# Fields that count as naming a source for check 5. A DOI is the strong form; the
# typed identifier pair is the form used for genuinely DOI-less sources (theses,
# patents, journals that register no DOI) -- see the WAIVERS history note.
BARRIER_SOURCE_FIELDS = frozenset({"doi", "source_doi", "identifier", "citation"})

# 2026-08-27 (Wave I fix 4). The structured sweep above only ever saw fields
# NAMED for a DOI inside data files. The flagship mechanism anchor that Wave I
# had to repair -- 10.1021/jf60200a038, live but resolving to a gossypol/rat
# paper -- lived in none of them: it was sitting in Python comments, in a module
# docstring and inside a free-text calibration rationale string, i.e. in exactly
# the places a reader trusts most and the gate could not see. These globs close
# that hole. They are scanned as TEXT: any DOI-shaped token anywhere in the file
# is extracted and run through checks 1 and 2 (syntax, confabulation signature).
# Checks 3 and 4 (status coherence, repair-record completeness) are properties of
# a structured RECORD and do not apply to prose, so they are not run here.
TEXT_SCAN_GLOBS = (
    "src/**/*.py",
    "docs/**/*.md",
    "README.md",
    "AUDIT.md",
)

# A DOI-shaped token embedded in free text. Deliberately greedy on the suffix and
# then trimmed by `_trim_text_doi`, because prose ends DOIs with sentence
# punctuation, wraps them in brackets/backticks/quotes, and puts them in
# parentheses. Anything this misses is a false NEGATIVE (the gate stays silent),
# never a false positive that blocks a build on punctuation.
#
# `\` and `#` terminate the token, and neither is a loss: a DOI suffix cannot
# contain a backslash (in the archived research dumps under docs/research/ they
# are MARKDOWN ESCAPES, e.g. `10.1016/0891-5849\`), and `#` starts a URL fragment
# (`...1573830#:~:text=Generally%2C...`), which is part of the link, not the DOI.
# Both were producing pure-punctuation `doi-syntax` failures on first run.
TEXT_DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>()\[\]{},;\\#]+")

# Trailing characters that are prose, not part of the DOI. `.` is included
# because a DOI never legitimately ends in a full stop, and sentences do.
_TEXT_DOI_TRAILING = ".,;:`'\"*_-"

# Field names that are asserted to hold a DOI.
DOI_FIELDS = frozenset(
    {
        "doi",
        "source_doi",
        "old_doi",
        "new_doi",
        "literature_anchor_doi",
        "related_published_doi",
    }
)

REPAIR_FIELDS = frozenset({"doi_repair", "citation_repair"})
REQUIRED_REPAIR_KEYS = ("old", "new", "date", "basis")

# https://www.crossref.org/blog/dois-and-matching-regular-expressions/ (permissive form)
DOI_RE = re.compile(r"^10\.\d{4,9}/[-._;()/:A-Za-z0-9+<>\[\]]+$")

# --- Confabulation signatures. NOT waivable. -------------------------------
# Each entry is (compiled pattern, human explanation). These encode the exact
# shapes the 2026-08-26 sweep found among fabricated anchors.
CONFABULATION_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (
        re.compile(r"/(?:[a-z0-9.]*\.)?[a-z]{2,}(?:_[a-z]{2,})*_(?:19|20)\d{2}$"),
        "author-name-and-year embedded in the DOI suffix "
        "(e.g. 10.1021/acs.jafc.de_leyn_2019) - a model-invented identifier, "
        "not a registered DOI",
    ),
    (
        re.compile(r"\.0{4,}$"),
        "all-zero article suffix (e.g. 10.1016/j.foodchem.2015.00000) - "
        "a placeholder standing in for an unknown article number",
    ),
    (
        # Delimited-token match only. A bare substring match is unusable here:
        # "na" appears inside legitimate registrants such as 10.1073/pnas.*, and
        # "tbd"/"todo" can occur inside real article suffixes.
        re.compile(r"(?i)(?:^|[/._\-])(?:x{3,}|placeholder|todo|tbd|fixme|example|your[-_]?doi)(?:[/._\-]|$)"),
        "template placeholder text left inside a DOI",
    ),
    (
        re.compile(r"(?i)^10\.x"),
        "template registrant prefix (10.xxxx/...)",
    ),
    (
        # Elsevier "j.<journal>.<year>.<article-number>" DOIs use a 6-digit
        # article number allocated from 100000 upward; a zero-padded one
        # (10.1016/j.foodres.2025.001279) is an invented number, not a real
        # allocation. Verified against all 238 DOIs currently in the tree: no
        # legitimate anchor matches. Page-based legacy Elsevier DOIs
        # (10.1016/j.foodchem.2011.07.037) have a 2-digit issue segment and are
        # unaffected.
        re.compile(r"(?i)^10\.1016/j\.[a-z]+\.(?:19|20)\d{2}\.0\d{5}$"),
        "zero-padded Elsevier article number (e.g. 10.1016/j.foodres.2025.001279) - "
        "outside the allocated 6-digit range, an invented identifier",
    ),
)

# --- Known-debt baseline ---------------------------------------------------
# (check_id, repo-relative path, json pointer, exact value, why it is tolerated)
# Every entry must still be *fixed*; the waiver only stops the gate from failing
# the build today. Removing debt makes a waiver stale, which this gate reports.
WAIVERS: tuple[tuple[str, str, str, str, str], ...] = (
    # THE BASELINE IS EMPTY. Keep it that way: a new waiver is an admission of
    # new citation debt and needs a dated justification here, in the ledger.
    #
    # History of the ratchet:
    #  - The nine stale `no_verifiable_source` waivers that shipped with this
    #    gate were resolved on 2026-08-27 by clearing the stale flags in
    #    data/lit/deep_research_backlog.json.
    #  - The last four (2026-08-27) were the genuinely DOI-less sources whose
    #    identifiers were being stored in fields named `doi` / `source_doi`:
    #      * Huang (2022) Clemson TigerPrints MS thesis 3936 -> identifier
    #        https://open.clemson.edu/all_theses/3936/, identifier_scheme "url"
    #      * Fraser et al. (2018) US 9,943,096 B2 -> identifier "US9943096B2",
    #        identifier_scheme "patent"
    #      * Cantre et al. (2007) Philippine Agricultural Scientist 90(2):143-152
    #        -> identifier_scheme "journal_locator" (the journal registers no DOI)
    #      * Liu (2023) NC State thesis, the PPI off-note hold-out anchor ->
    #        identifier_scheme "citation" (no DOI or handle exists)
    #    Each now carries a typed `identifier` + `identifier_scheme` pair plus an
    #    `identifier_note` recording the retyping. `src.benchmark_validation.
    #    matrix_source_anchor` reads that pair alongside `source_doi`, so the
    #    external-data-status classification of every benchmark — the hold-out
    #    bundle's `external_validation_only` in particular — is unchanged
    #    (verified by a before/after classification diff over all 21 benchmark
    #    files; see tasks/audit_remediation.md).
)

# --- Text-surface carry-over: NOT a baseline, an UNPAID DEFECT --------------
# (check_id, repo-relative path, exact value, what is actually wrong)
#
# 2026-08-27 (Wave I fix 4). Widening the scan to code and prose immediately
# found a fabricated DOI in a document that predates the widening. This list is
# deliberately NOT part of WAIVERS and is NOT a known-debt baseline: WAIVERS is
# a ratchet for tolerated debt, whereas every entry here is a LIVE CONFABULATION
# that a named owner has to delete. It exists for one reason only -- so that
# turning the scan on does not red-line CI for work that is not the scan's fault
# -- and it is printed under its own alarmed heading on every run so it cannot
# become invisible.
#
# Rules for this list, which differ from WAIVERS on purpose:
#   * It applies ONLY to files matched by TEXT_SCAN_GLOBS. The structured
#     surface keeps its original rule that confabulation signatures are never
#     waivable; nothing here weakens that.
#   * Entries are keyed WITHOUT a pointer, because a line number in prose churns
#     for reasons that have nothing to do with the citation.
#   * Adding an entry requires naming the defect, not describing it as tolerable.
#   * An entry that stops matching is reported as stale and must be deleted, so
#     this list can only shrink.
# EMPTY, and it should stay that way.
#
# Ratchet history:
#   2026-08-27 (Wave I, opened): one entry — a confabulated Elsevier DOI in
#     docs/slr_benchmark_evaluation.md §7.3, surfaced the moment the scan was widened to
#     cover docs/**, README.md and AUDIT.md. It was carried for a matter of hours because
#     docs/** sat outside the editing scope of the sub-task that widened the scan.
#   2026-08-27 (Wave I, closed): FIXED AT SOURCE, not waived. The Hao et al. (2025)
#     reference was WITHDRAWN in the document: its score is struck, its "directly
#     calibrates the optimizer" claim is retracted (verified: no constant in the repo is
#     fitted to it — the identifier appears nowhere else), and the bad identifier is
#     retained in the withdrawal notice with spaces inserted so the record survives
#     without re-asserting the anchor.
#
# The gate reports a waiver that no longer matches as STALE and fails the build, so this
# list can only shrink. Adding to it requires a reason that would survive being read
# aloud.
TEXT_SURFACE_WAIVERS: tuple[tuple[str, str, str, str], ...] = ()


class Violation:
    __slots__ = ("check", "path", "pointer", "value", "detail")

    def __init__(self, check: str, path: str, pointer: str, value: str, detail: str) -> None:
        self.check = check
        self.path = path
        self.pointer = pointer
        self.value = value
        self.detail = detail

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.check, self.path, self.pointer, self.value)

    def __str__(self) -> str:
        return f"[{self.check}] {self.path}{self.pointer}\n    value: {self.value!r}\n    {self.detail}"


def _load(path: Path) -> Any:
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".json":
        return json.loads(text)
    import yaml  # imported lazily so the JSON-only path needs no PyYAML

    return yaml.safe_load(text)


def _iter_files() -> Iterator[Path]:
    seen: set[Path] = set()
    for pattern in SCAN_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if path.is_file() and path not in seen:
                seen.add(path)
                yield path


def _iter_text_files() -> Iterator[Path]:
    seen: set[Path] = set()
    for pattern in TEXT_SCAN_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if not path.is_file() or path in seen:
                continue
            # Build artefacts and caches are not sources of citations.
            parts = set(path.parts)
            if "__pycache__" in parts or ".git" in parts or any(
                part.endswith(".egg-info") for part in path.parts
            ):
                continue
            seen.add(path)
            yield path


def _trim_text_doi(raw: str) -> str:
    """Strip prose punctuation that the greedy text regex swallowed."""
    token = raw.strip()
    while token and token[-1] in _TEXT_DOI_TRAILING:
        token = token[:-1]
    return token


def _check_doi_token(raw: str, rel: str, pointer: str, out: list[Violation]) -> None:
    """Checks 1 and 2 on a single DOI string, from any source."""
    for pattern, why in CONFABULATION_PATTERNS:
        if pattern.search(raw):
            out.append(Violation("confabulation-signature", rel, pointer, raw, why))
            return
    if not DOI_RE.match(raw):
        out.append(
            Violation(
                "doi-syntax",
                rel,
                pointer,
                raw,
                "does not match the DOI grammar 10.<registrant>/<suffix>. "
                "If this is a non-DOI identifier (thesis, patent, standard) "
                "it does not belong in a field named for a DOI.",
            )
        )


def _scan_text(path: Path, rel: str, out: list[Violation], dois: list[str]) -> int:
    """Extract DOI-shaped tokens from a source/prose file and check them.

    Returns the number of tokens found. The pointer is `#L<line>` so a failure
    names the line, which is the only useful coordinate in an unstructured file.
    """
    found = 0
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        out.append(
            Violation("parse", rel, "", "", f"could not read: {type(exc).__name__}: {exc}")
        )
        return 0

    if "10." not in text:  # cheap bail-out; most source files have no DOI at all
        return 0

    for line_number, line in enumerate(text.splitlines(), start=1):
        for match in TEXT_DOI_RE.finditer(line):
            token = _trim_text_doi(match.group(0))
            if not token:
                continue
            found += 1
            dois.append(token)
            _check_doi_token(token, rel, f"#L{line_number}", out)
    return found


def _walk(node: Any, rel: str, pointer: str, out: list[Violation], dois: list[str]) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            child = f"{pointer}/{key}"

            if key in DOI_FIELDS and isinstance(value, str) and value.strip():
                raw = value.strip()
                dois.append(raw)
                _check_doi_token(raw, rel, child, out)

            if key in REPAIR_FIELDS and isinstance(value, dict):
                # `new` may be explicitly null: that is the deliberate
                # "DOI withdrawn, no defensible replacement found" record, paired
                # with doi=null and source_status=no_verifiable_source. The key
                # must still be present, so the withdrawal is stated rather than
                # implied by omission. old/date/basis must be non-empty.
                missing = [
                    k
                    for k in REQUIRED_REPAIR_KEYS
                    if k not in value
                    or (value[k] in (None, "") if k != "new" else value[k] == "")
                ]
                if missing:
                    out.append(
                        Violation(
                            "repair-record",
                            rel,
                            child,
                            ",".join(sorted(value)),
                            f"repair record is missing required key(s): {', '.join(missing)}. "
                            "Every repair must state what it replaced, what it replaced it with, "
                            "when, and on what evidence.",
                        )
                    )

            _walk(value, rel, child, out, dois)

        status = node.get("source_status")
        doi = node.get("doi")
        if (
            isinstance(status, str)
            and status.strip() == "no_verifiable_source"
            and isinstance(doi, str)
            and doi.strip()
        ):
            out.append(
                Violation(
                    "status-coherence",
                    rel,
                    pointer,
                    doi.strip(),
                    "record asserts source_status='no_verifiable_source' while also carrying a DOI. "
                    "One of the two is stale: either the DOI verifies the source (clear the flag) "
                    "or it does not (drop the DOI).",
                )
            )

    elif isinstance(node, list):
        for index, value in enumerate(node):
            _walk(value, rel, f"{pointer}[{index}]", out, dois)


def run_offline() -> tuple[list[Violation], list[Violation], list[tuple[str, ...]], list[str], int]:
    """Return (blocking, waived, stale_waivers, all_dois, files_scanned).

    Kept at five elements deliberately: it is the published shape, used by
    tests/unit/test_audit_remediation_carried_2026_08.py. The text-surface
    carry-over added on 2026-08-27 is reported through `run_offline_detailed`.
    """
    blocking, waived, stale, _carried, _stale_carried, dois, files_scanned = (
        run_offline_detailed()
    )
    return blocking, waived, stale, dois, files_scanned


def _check_barrier_sources(path: Path, rel: str, out: list[Violation]) -> None:
    """Check 5: every barrier benchmark entry names a source or admits it has none.

    Added 2026-08-27 (Wave J1). The failure this prevents is not a bad citation but
    the ABSENCE of one: nine literature windows and eighteen "computed" columns
    shipped with no DOI, no author-year, no run record and no ledger entry, and
    nothing in the repository could tell you that. Silence read as anchored.

    A `no_verifiable_source` label satisfies this check. That is deliberate: the
    check enforces DISCLOSURE, not quality. Check 3 (status coherence) still
    forbids claiming the label and a DOI at once, so an entry cannot satisfy both
    branches and pass as anchored.
    """
    try:
        payload = _load(path)
    except Exception:
        return  # the parse violation is already recorded by the caller
    if not isinstance(payload, dict):
        return
    benchmarks = payload.get("benchmarks")
    if not isinstance(benchmarks, list):
        return

    for index, entry in enumerate(benchmarks):
        if not isinstance(entry, dict):
            continue
        family = str(entry.get("family", f"index {index}"))
        status = entry.get("source_status")
        if isinstance(status, str) and status.strip() == "no_verifiable_source":
            continue
        if any(
            isinstance(entry.get(field), str) and entry[field].strip()
            for field in BARRIER_SOURCE_FIELDS
        ):
            continue
        out.append(
            Violation(
                "barrier-source",
                rel,
                f"/benchmarks[{index}]",
                family,
                "barrier benchmark entry names no source. Every window and every "
                "claimed computed value must carry one of "
                f"{sorted(BARRIER_SOURCE_FIELDS)} or the explicit admission "
                "source_status='no_verifiable_source'. An entry with neither is "
                "indistinguishable from a number someone typed.",
            )
        )


def run_offline_detailed() -> tuple[
    list[Violation],
    list[Violation],
    list[tuple[str, ...]],
    list[Violation],
    list[tuple[str, ...]],
    list[str],
    int,
]:
    """Return (blocking, waived, stale_waivers, text_carried, stale_text_carried,
    all_dois, files_scanned)."""
    violations: list[Violation] = []
    dois: list[str] = []
    files_scanned = 0

    for path in _iter_files():
        rel = path.relative_to(ROOT).as_posix()
        try:
            payload = _load(path)
        except Exception as exc:  # malformed data is itself a gate failure
            violations.append(
                Violation("parse", rel, "", "", f"could not parse: {type(exc).__name__}: {exc}")
            )
            continue
        files_scanned += 1
        _walk(payload, rel, "", violations, dois)

    # 2026-08-27 (Wave J1): check 5, barrier-source disclosure. Scoped to the two
    # QM benchmark files rather than to all of SCAN_GLOBS, because "must name a
    # source" is only well defined for a record that asserts a physical quantity.
    for pattern in BARRIER_SOURCE_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            if path.is_file():
                _check_barrier_sources(path, path.relative_to(ROOT).as_posix(), violations)

    # 2026-08-27 (Wave I fix 4): unstructured pass over code and prose.
    for path in _iter_text_files():
        rel = path.relative_to(ROOT).as_posix()
        if _scan_text(path, rel, violations, dois):
            files_scanned += 1

    waiver_index = {(c, p, ptr, v): why for c, p, ptr, v, why in WAIVERS}
    # Text-surface carry-over is keyed WITHOUT the pointer; see TEXT_SURFACE_WAIVERS.
    text_index = {(c, p, v) for c, p, v, _why in TEXT_SURFACE_WAIVERS}
    text_paths = {
        path.relative_to(ROOT).as_posix() for path in _iter_text_files()
    }

    blocking: list[Violation] = []
    waived: list[Violation] = []
    carried: list[Violation] = []
    matched: set[tuple[str, str, str, str]] = set()
    matched_text: set[tuple[str, str, str]] = set()

    for violation in violations:
        text_key = (violation.check, violation.path, violation.value)
        if violation.path in text_paths and text_key in text_index:
            matched_text.add(text_key)
            carried.append(violation)
        # Confabulation signatures are never waivable on the structured surface.
        elif violation.check != "confabulation-signature" and violation.key in waiver_index:
            matched.add(violation.key)
            waived.append(violation)
        else:
            blocking.append(violation)

    stale = [w for w in WAIVERS if (w[0], w[1], w[2], w[3]) not in matched]
    stale_text = [w for w in TEXT_SURFACE_WAIVERS if (w[0], w[1], w[2]) not in matched_text]
    return blocking, waived, stale, carried, stale_text, dois, files_scanned


def run_online(dois: list[str], sample: int, seed: int, timeout: float) -> list[tuple[str, str]]:
    """Resolve a random sample of DOIs against doi.org. Returns [(doi, reason), ...] for failures."""
    import urllib.error
    import urllib.request

    unique = sorted({d for d in dois if DOI_RE.match(d)})
    rng = random.Random(seed)
    picked = unique if sample >= len(unique) else rng.sample(unique, sample)
    dead: list[tuple[str, str]] = []

    for doi in picked:
        url = f"https://doi.org/{doi}"
        request = urllib.request.Request(
            url,
            method="HEAD",
            headers={"User-Agent": "maillard-citation-gate/1.0 (+https://github.com/) mailto:noreply@example.com"},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if response.status >= 400:
                    dead.append((doi, f"HTTP {response.status}"))
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 410):
                dead.append((doi, f"HTTP {exc.code}"))
            # 3xx/4xx from the *publisher* after a successful redirect is not a
            # dead DOI; only doi.org's own 404/410 proves non-registration.
        except Exception as exc:  # network hiccup: advisory mode, never fatal
            print(f"  ? {doi}: transport error ({type(exc).__name__}) - skipped", file=sys.stderr)

    print(f"\nOnline sample: {len(picked)} of {len(unique)} unique DOIs resolved.")
    return dead


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--online", action="store_true", help="also resolve a random DOI sample against doi.org")
    parser.add_argument("-n", "--sample", type=int, default=25, help="DOIs to sample in --online mode (default 25)")
    parser.add_argument("--seed", type=int, default=0, help="sample seed (0 = nondeterministic)")
    parser.add_argument("--timeout", type=float, default=15.0, help="per-DOI HTTP timeout, seconds")
    parser.add_argument(
        "--fail-on-online",
        action="store_true",
        help="exit non-zero on dead DOIs found in --online mode (default: advisory only)",
    )
    args = parser.parse_args(argv)

    (
        blocking,
        waived,
        stale,
        carried,
        stale_carried,
        dois,
        files_scanned,
    ) = run_offline_detailed()

    print(f"citation_gate: scanned {files_scanned} files, {len(dois)} DOI-bearing fields, "
          f"{len(set(dois))} unique DOIs.")

    if waived:
        print(f"\n{len(waived)} known-debt violation(s) waived by the baseline in this script:")
        for violation in waived:
            print(f"  - [{violation.check}] {violation.path}{violation.pointer} -> {violation.value}")

    if carried:
        print(
            f"\n!! {len(carried)} UNFIXED CITATION DEFECT(S) IN CODE/PROSE, carried by "
            f"TEXT_SURFACE_WAIVERS so the widened scan does not red-line CI.\n"
            f"!! These are NOT tolerated debt. Each one is a live defect awaiting an owner."
        )
        for violation in carried:
            why = next(
                (w[3] for w in TEXT_SURFACE_WAIVERS
                 if (w[0], w[1], w[2]) == (violation.check, violation.path, violation.value)),
                "",
            )
            print(f"  - [{violation.check}] {violation.path}{violation.pointer} -> {violation.value}")
            print(f"      {why}")

    exit_code = 0

    if stale:
        print(f"\nFAIL: {len(stale)} stale waiver(s) - the underlying violation is gone, "
              f"so the waiver must be deleted from WAIVERS in this file:")
        for check, path, pointer, value, _why in stale:
            print(f"  - [{check}] {path}{pointer} -> {value}")
        exit_code = 1

    if stale_carried:
        print(f"\nFAIL: {len(stale_carried)} stale text-surface carry-over(s) - the defect is "
              f"fixed, so the entry must be deleted from TEXT_SURFACE_WAIVERS in this file:")
        for check, path, value, _why in stale_carried:
            print(f"  - [{check}] {path} -> {value}")
        exit_code = 1

    if blocking:
        print(f"\nFAIL: {len(blocking)} citation violation(s):\n")
        for violation in blocking:
            print(str(violation) + "\n")
        exit_code = 1

    if args.online:
        seed = args.seed if args.seed else random.randrange(1 << 30)
        dead = run_online(dois, args.sample, seed, args.timeout)
        if dead:
            print(f"\n{len(dead)} DOI(s) did not resolve (seed={seed}):")
            for doi, reason in dead:
                print(f"  - {doi}: {reason}")
            if args.fail_on_online:
                exit_code = 1
        else:
            print("All sampled DOIs resolved.")

    if exit_code == 0:
        print("\ncitation_gate: PASS")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
