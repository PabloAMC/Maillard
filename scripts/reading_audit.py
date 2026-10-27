#!/usr/bin/env python3
"""
reading_audit.py — which papers on this machine have been read, and which have not.

WHY THIS EXISTS (2026-09-09). The repository's rule is "no number without a dossier". The
inverse has never been checked: a PDF can sit in ``data/articles/`` for months with nothing
pointing at it, and nobody notices. On 9 September 2026 an audit by hand found fifty-five such
files, several of which answered gaps the repository names in its own artifacts (a second
laboratory's acrylamide constants; the amino-acid identity ratios two refused waves asked for;
the competing branch a third refused wave lacked). This command makes that gap countable
instead of accidental.

It is NOT a gate and NOT a tracked artifact. ``data/keys/papers.yml`` deliberately ignores the
local PDFs so that its output is the same on a machine that does not have them; this command is
the other half, and its answer depends on what is on the disk it runs on.

A PDF counts as READ when a dossier's file name matches it once both are reduced to letters and
digits, either way round (``charles-bernard2005.pdf`` matches ``charlesbernard2005_extraction.md``;
``Zhang2024.pdf`` matches ``Zhang2024_extraction.md`` and ``zhang2024b_extraction.md``). That is
a name match, not a claim that the dossier is about that paper: a dossier written from an
article whose PDF is not on disk, or a PDF read into a synthesis dossier under another name
(several of the older ``k1_``/``k2_``/``k3_`` inventories do this), is reported separately as
UNMATCHED rather than silently counted either way.

Usage:
    python scripts/reading_audit.py              # the report
    python scripts/reading_audit.py --json       # the same as JSON, for a script
    python scripts/reading_audit.py --quiet      # counts only
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_access, data_paths  # noqa: E402

DOSSIER_DIR = data_paths.EXTRACTION_DOSSIERS_DIR
ARTICLES_DIR = data_paths.ARTICLES_DIR
#: Dossiers that are syntheses over many papers, not one paper's reading; they never match a PDF.
SYNTHESIS_PREFIXES = ("k1_", "k2_", "k3_", "k4_", "k5_", "k6_", "research_", "wave_")
#: PDFs whose dossier is named after the paper rather than the file, and PDFs triaged without one,
#: live in data/keys/pdf_dossier_aliases.yml (2026-10-27; until then a map in this file). The same
#: table drives scripts/ci/dossier_coverage_gate.py, so the audit and the gate cannot disagree.
ALIAS_TABLE = data_paths.PDF_DOSSIER_ALIASES


def _alias_table() -> Dict[str, Dict[str, object]]:
    payload = data_access.load_yaml(ALIAS_TABLE) if ALIAS_TABLE.exists() else None
    payload = payload or {}
    return {section: dict(payload.get(section) or {}) for section in ("aliases", "triaged", "pending")}


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def audit() -> Dict[str, object]:
    pdfs = sorted(p.name for p in ARTICLES_DIR.glob("*.pdf")) if ARTICLES_DIR.exists() else []
    dossiers = sorted(p.name for p in DOSSIER_DIR.glob("*.md"))
    per_paper = [d for d in dossiers if not d.startswith(SYNTHESIS_PREFIXES)]
    stems = {d: _norm(d.replace("_extraction.md", "").replace(".md", "")) for d in per_paper}

    table = _alias_table()
    read: Dict[str, List[str]] = {}
    unread: List[str] = []
    triaged: Dict[str, str] = {}
    for pdf in pdfs:
        stem = Path(pdf).stem
        key = _norm(stem)
        if stem in table["aliases"]:
            targets = table["aliases"][stem]
            targets = [targets] if isinstance(targets, str) else list(targets)
            hits = [f for t in targets for f in (f"{t}_extraction.md", f"{t}.md") if f in dossiers]
            if hits:
                read[pdf] = sorted(hits)
                continue
        if stem in table["triaged"]:
            triaged[pdf] = str(table["triaged"][stem])
            continue
        hits = [d for d, stem in stems.items() if len(stem) > 4 and (key.startswith(stem) or stem.startswith(key))]
        if hits:
            read[pdf] = sorted(hits)
        else:
            unread.append(pdf)

    matched = {d for hits in read.values() for d in hits}
    unmatched = sorted(d for d in per_paper if d not in matched)
    return {
        "articles_dir": data_paths.rel(ARTICLES_DIR),
        "counts": {"pdfs": len(pdfs), "dossiers": len(dossiers), "per_paper_dossiers": len(per_paper),
                   "pdfs_read": len(read), "pdfs_triaged": len(triaged), "pdfs_unread": len(unread), "dossiers_without_a_pdf_here": len(unmatched)},
        "unread": unread,
        "triaged": triaged,
        "dossiers_without_a_pdf_here": unmatched,
        "note": ("a name match, not a verified identity; a PDF read into one of the synthesis inventories under "
                 "another name shows as unread, and a dossier written from an article this machine does not hold "
                 "shows as having no PDF"),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--json", action="store_true", help="print the report as JSON")
    parser.add_argument("--quiet", action="store_true", help="print the counts only")
    args = parser.parse_args(argv)
    report = audit()
    if args.json:
        print(json.dumps(report, indent=2))
        return 0
    c = report["counts"]
    print(f"{c['pdfs']} PDFs under {report['articles_dir']}; {c['per_paper_dossiers']} per-paper dossiers "
          f"({c['dossiers']} dossier files in all).")
    print(f"  read:   {c['pdfs_read']}")
    print(f"  triaged without a dossier: {c['pdfs_triaged']}")
    print(f"  UNREAD: {c['pdfs_unread']}")
    print(f"  dossiers whose paper is not on this disk: {c['dossiers_without_a_pdf_here']}")
    if args.quiet:
        return 0
    if report["unread"]:
        print("\nPDFs with no dossier of a matching name (read them, or record why not):")
        for name in report["unread"]:
            print(f"  {name}")
    print("\n" + report["note"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
