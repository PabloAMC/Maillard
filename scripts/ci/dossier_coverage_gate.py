#!/usr/bin/env python
"""Dossier coverage gate.

On 2026-10-27 the owner asked why there were unread papers on disk. Of the PDFs in
``data/articles/`` thirteen had no dossier of their name: eleven had been read under another
name, two never had been, and one of those two (Cerny 2007 LWT) had been USED since wave B2 --
``k3_final_parameter_inventory.md`` cited a ``cerny2007b`` dossier that did not exist. Nothing
checked either direction. This gate checks both (tasks/data_restructure_plan.md, "PDFs on disk
without an extraction record").

Checks
------
(a) **Every dossier a file names exists.** Every ``<name>_extraction.md`` written anywhere in a
    tracked text file under src/, scripts/, data/, docs/, tests/ or tasks/ must be a file in the
    extraction-dossier directory. Names are compared lower-cased and reduced to letters and
    digits, so ``Zhou2023_extraction.md`` finds ``zhou2023_extraction.md``. A sentence that has
    to name a dossier that does not exist should name its stem without the suffix.
(b) **Every article on disk has a record.** Every ``.pdf``, ``.docx`` or ``.txt`` in the
    articles directory has a dossier of the same stem (same normalisation), or an entry in the
    hand-curated table ``data/keys/pdf_dossier_aliases.yml``: ``aliases`` (the dossier(s) that
    read it), ``triaged`` (why it needs none) or ``pending`` (downloaded, not yet read). The
    articles are local only (gitignored), so when the directory is absent -- CI, a worktree --
    check (b) is skipped with a note; the parts of the table that need no PDFs still run:
    every alias target exists, every triage line has a reason, no stem is listed twice, and no
    ``pending`` stem has meanwhile gained a dossier of its name (a stale queue line fails, so the
    queue can only shrink honestly).

Stdlib + PyYAML, offline, like the other gates in this tier.

Usage:
    python scripts/ci/dossier_coverage_gate.py
    python scripts/ci/dossier_coverage_gate.py --articles-dir /path/to/main/checkout/data/articles
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import yaml  # noqa: E402

from src import data_paths  # noqa: E402

SCAN_TREES = ("src", "scripts", "data", "docs", "tests", "tasks")
ARTICLE_SUFFIXES = (".pdf", ".docx", ".txt")
EXTRACTION_SUFFIX = "_extraction.md"
#: A dossier file name in running text. It must start with a letter (``zamora2007/2008/2015_extraction.md``
#: is shorthand, not a name) and must not follow a backslash (``\nyiltirak2026_extraction.md`` in a
#: string literal is a line break, then the name).
DOSSIER_REFERENCE = re.compile(r"(?<![A-Za-z0-9_.\-\\])([A-Za-z][A-Za-z0-9_.\-]*" + re.escape(EXTRACTION_SUFFIX) + r")")
TRIAGED_PREFIX = "triaged:"


def norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def dossier_stem(filename: str) -> str:
    if filename.endswith(EXTRACTION_SUFFIX):
        return filename[: -len(EXTRACTION_SUFFIX)]
    return filename[: -len(".md")] if filename.endswith(".md") else filename


def dossier_index(dossier_dir: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """(normalised stem of every ``*_extraction.md`` -> file name, exact stem of every ``*.md`` -> file name)."""
    files = sorted(p.name for p in dossier_dir.glob("*.md")) if dossier_dir.is_dir() else []
    extraction = {norm(dossier_stem(f)): f for f in files if f.endswith(EXTRACTION_SUFFIX)}
    every = {dossier_stem(f): f for f in files}
    return extraction, every


# --------------------------------------------------------------------------- check (a)


def tracked_text_files(root: Path, trees: Sequence[str] = SCAN_TREES) -> List[Path]:
    out = subprocess.run(["git", "ls-files", "-z", "--", *trees], cwd=root, capture_output=True, check=True).stdout
    return [root / p for p in out.decode("utf-8").split("\0") if p]


def missing_references(files: Iterable[Path], dossier_dir: Path, root: Path) -> Dict[str, List[str]]:
    """Dossier names cited in ``files`` that resolve to no ``*_extraction.md`` -> the files citing them."""
    extraction, _ = dossier_index(dossier_dir)
    missing: Dict[str, List[str]] = {}
    for path in files:
        try:
            raw = path.read_bytes()
        except OSError:  # tracked but deleted in the working tree
            continue
        if b"\0" in raw[:8192]:
            continue
        for name in set(DOSSIER_REFERENCE.findall(raw.decode("utf-8", errors="replace"))):
            if norm(dossier_stem(name)) not in extraction:
                try:
                    where = str(path.relative_to(root))
                except ValueError:
                    where = str(path)
                missing.setdefault(name, []).append(where)
    return {k: sorted(v) for k, v in sorted(missing.items())}


# --------------------------------------------------------------------------- check (b)


def load_alias_table(path: Path) -> Dict[str, Dict[str, object]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) if path.exists() else None
    payload = payload or {}
    return {section: dict(payload.get(section) or {}) for section in ("aliases", "triaged", "pending")}


def alias_table_problems(table: Dict[str, Dict[str, object]], dossier_dir: Path) -> List[str]:
    """Defects in the table itself; none of these needs the articles on disk."""
    extraction, every = dossier_index(dossier_dir)
    problems: List[str] = []
    seen: Dict[str, str] = {}
    for section, entries in table.items():
        for stem in entries:
            if stem in seen:
                problems.append(f"{stem}: listed under both `{seen[stem]}` and `{section}`")
            seen[stem] = section
    for stem, targets in table["aliases"].items():
        targets = [targets] if isinstance(targets, str) else list(targets or [])
        if not targets:
            problems.append(f"{stem}: alias with no dossier")
        for target in targets:
            if str(target) not in every:
                problems.append(f"{stem}: alias target `{target}` is not a dossier in {data_paths.rel(dossier_dir)}")
    for stem, reason in table["triaged"].items():
        reason = str(reason or "").strip()
        if reason.startswith(TRIAGED_PREFIX):
            reason = reason[len(TRIAGED_PREFIX):].strip()
        if not reason:
            problems.append(f"{stem}: triaged with no reason")
    for stem in table["pending"]:
        if norm(stem) in extraction:
            problems.append(
                f"{stem}: STALE pending line -- {extraction[norm(stem)]} now exists; delete the line"
            )
    return problems


def uncovered_articles(articles_dir: Path, dossier_dir: Path, table: Dict[str, Dict[str, object]]) -> Tuple[List[str], List[str]]:
    """(articles with no record, table stems with no file on this disk)."""
    extraction, _ = dossier_index(dossier_dir)
    listed = set(table["aliases"]) | set(table["triaged"]) | set(table["pending"])
    articles = sorted(p for p in articles_dir.iterdir() if p.is_file() and p.suffix.lower() in ARTICLE_SUFFIXES)
    uncovered = [p.name for p in articles if p.stem not in listed and norm(p.stem) not in extraction]
    on_disk = {p.stem for p in articles}
    absent = sorted(stem for stem in listed if stem not in on_disk)
    return uncovered, absent


# --------------------------------------------------------------------------- main


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--articles-dir",
        type=Path,
        default=data_paths.ARTICLES_DIR,
        help="where the local article files live (default: %(default)s); check (b) is skipped if absent",
    )
    args = parser.parse_args(argv)
    dossier_dir = data_paths.EXTRACTION_DOSSIERS_DIR
    failures: List[str] = []

    missing = missing_references(tracked_text_files(ROOT), dossier_dir, ROOT)
    for name, where in missing.items():
        failures.append(f"(a) {name} does not exist; cited in {', '.join(where[:4])}" + (" ..." if len(where) > 4 else ""))

    table = load_alias_table(data_paths.PDF_DOSSIER_ALIASES)
    failures += [f"(b) {p}" for p in alias_table_problems(table, dossier_dir)]

    articles_dir = args.articles_dir
    notes: List[str] = []
    if articles_dir.is_dir():
        uncovered, absent = uncovered_articles(articles_dir, dossier_dir, table)
        for name in uncovered:
            failures.append(
                f"(b) {name} has no dossier of its name and no line in {data_paths.rel(data_paths.PDF_DOSSIER_ALIASES)}"
            )
        notes.append(f"checked every article in {articles_dir}; {len(table['pending'])} pending reading")
        if absent:
            notes.append(f"{len(absent)} table stem(s) have no file in {articles_dir}: {', '.join(absent[:6])}")
    else:
        notes.append(f"{articles_dir} is absent; per-article check (b) skipped (the table itself was checked)")

    # Notes first: the last line is what the model card quotes, so it must depend neither on
    # whether this machine holds the articles nor on how long the reading queue is.
    notes.append(f"{sum(len(v) for v in table.values())} lines in {data_paths.rel(data_paths.PDF_DOSSIER_ALIASES)}")
    for note in notes:
        print(f"dossier_coverage_gate: note -- {note}")
    if not failures:
        print("dossier_coverage_gate: PASS (every cited dossier exists; the article table is coherent)")
        return 0
    print(f"dossier_coverage_gate: FAIL -- {len(failures)} problem(s)", file=sys.stderr)
    for line in failures:
        print(f"  {line}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
