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

Checks (offline, all blocking)
------------------------------
1. **DOI syntax** -- every DOI-bearing field across ``data/lit/*.json|yml``,
   ``data/benchmarks/**/*.json`` and ``data/species/*.yml`` matches the DOI grammar.
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
)

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


def _walk(node: Any, rel: str, pointer: str, out: list[Violation], dois: list[str]) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            child = f"{pointer}/{key}"

            if key in DOI_FIELDS and isinstance(value, str) and value.strip():
                raw = value.strip()
                dois.append(raw)
                for pattern, why in CONFABULATION_PATTERNS:
                    if pattern.search(raw):
                        out.append(
                            Violation("confabulation-signature", rel, child, raw, why)
                        )
                        break
                else:
                    if not DOI_RE.match(raw):
                        out.append(
                            Violation(
                                "doi-syntax",
                                rel,
                                child,
                                raw,
                                "does not match the DOI grammar 10.<registrant>/<suffix>. "
                                "If this is a non-DOI identifier (thesis, patent, standard) "
                                "it does not belong in a field named for a DOI.",
                            )
                        )

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
    """Return (blocking, waived, stale_waivers, all_dois, files_scanned)."""
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

    waiver_index = {(c, p, ptr, v): why for c, p, ptr, v, why in WAIVERS}
    blocking: list[Violation] = []
    waived: list[Violation] = []
    matched: set[tuple[str, str, str, str]] = set()

    for violation in violations:
        # Confabulation signatures are never waivable.
        if violation.check != "confabulation-signature" and violation.key in waiver_index:
            matched.add(violation.key)
            waived.append(violation)
        else:
            blocking.append(violation)

    stale = [w for w in WAIVERS if (w[0], w[1], w[2], w[3]) not in matched]
    return blocking, waived, stale, dois, files_scanned


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

    blocking, waived, stale, dois, files_scanned = run_offline()

    print(f"citation_gate: scanned {files_scanned} files, {len(dois)} DOI-bearing fields, "
          f"{len(set(dois))} unique DOIs.")

    if waived:
        print(f"\n{len(waived)} known-debt violation(s) waived by the baseline in this script:")
        for violation in waived:
            print(f"  - [{violation.check}] {violation.path}{violation.pointer} -> {violation.value}")

    exit_code = 0

    if stale:
        print(f"\nFAIL: {len(stale)} stale waiver(s) - the underlying violation is gone, "
              f"so the waiver must be deleted from WAIVERS in this file:")
        for check, path, pointer, value, _why in stale:
            print(f"  - [{check}] {path}{pointer} -> {value}")
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
