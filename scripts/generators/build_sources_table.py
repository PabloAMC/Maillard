#!/usr/bin/env python
"""
docs/guides/SOURCES.md -- every paper with an extraction dossier, alphabetically, with what the model takes
from it (measured rate constants, fit rows, benchmark pots, directional claims) and the dossier's own
one-line statement of what the repository uses. Counts come from `build_thiol_sink_figures.paper_usage_counts`.

    python scripts/generators/build_sources_table.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "scripts" / "generators"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from build_thiol_sink_figures import paper_usage_counts  # noqa: E402

OUT = ROOT / "docs" / "guides" / "SOURCES.md"
DOSSIERS = ROOT / "data" / "lit" / "extraction_dossiers"
TAKE_HEADING = re.compile(r"^#{2,3}\s*(?:§?\d+[a-z]?\.?\s*)?(?:what the repo(?:sitory)?.*(?:take|use|can use)|what (?:this paper|the paper|the sulfur module).*(?:suppl|takes)|which steps overlap)", re.I)
PAPER_STEM = re.compile(r"^([a-z\-]{2,})(\d{4})[a-z]*(?:_|$)", re.I)


def take_line(text: str) -> str:
    """The first sentence after the dossier's 'what the repo takes' heading, or the first bullet."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if TAKE_HEADING.match(line.strip()):
            body = []
            for l in lines[i + 1:]:
                if l.startswith("#"):
                    break
                s = l.strip().lstrip("-*> ").strip()
                if s.startswith("|") or s.startswith("```"):
                    continue          # tables and code carry the numbers; the sentence before them says what is taken
                if s:
                    body.append(s)
                if len(" ".join(body)) > 300:
                    break
            s = " ".join(body)
            s = re.sub(r"\*\*|`|\[|\]\([^)]*\)|[★⭐⚠️⇒]", "", s)
            s = re.sub(r"^\s*(TAKEN|SI: none(?: listed)?\.? Nothing further to retrieve\.?)[:\s]*", "", s, flags=re.I)
            s = re.sub(r"\s+", " ", s).strip()
            s = s.replace("|", "/")
            return (s[:200].rsplit(" ", 1)[0] + " …") if len(s) > 200 else s
    return ""


def norm(author: str) -> str:
    """'De Vleeschouwer', 'Vleeschouwer', 'Kocadağlı' and 'Kocadagli' all compare equal."""
    import unicodedata

    s = unicodedata.normalize("NFKD", author).encode("ascii", "ignore").decode().lower().replace("-", "")
    return re.sub(r"^(?:de|van|von|der)+", "", s)


def norm_key(key: str) -> str:
    a, y = key.rsplit(" ", 1)
    return f"{norm(a)} {y}"


def author_year(stem: str):
    """('Author', 'year') from a dossier file name such as 'hofmann1998_extraction'; None for a synthesis file."""
    m = PAPER_STEM.match(stem)
    return (m.group(1).capitalize(), m.group(2)) if m else None


def main() -> int:
    raw = paper_usage_counts()
    counts = {}
    for kind, c in raw.items():
        d = counts.setdefault(kind, {})
        for k, n in c.items():
            d[norm_key(k)] = d.get(norm_key(k), 0) + n
    reg = yaml.safe_load((ROOT / "data" / "keys" / "papers.yml").read_text(encoding="utf-8"))["papers"]
    by_dossier = {p["dossier"]: p for p in reg if p.get("dossier")}
    rows = []
    for f in sorted(DOSSIERS.glob("*.md")):
        ay = author_year(f.stem)
        if ay is None:
            continue              # wave syntheses and research-round notes are not papers
        author, year = ay
        key = norm_key(f"{author} {year}")
        entry = by_dossier.get(f.name) or by_dossier.get(str(f.relative_to(ROOT))) or {}
        citation = (entry.get("citation") or "").strip()
        doi = (entry.get("doi") or "").strip()
        text = f.read_text(encoding="utf-8")
        if not citation:
            title = next((l.lstrip("# ").strip() for l in text.splitlines() if l.startswith("# ")), f.stem)
            title = re.split(r"\s+[—–-]+\s+", title, maxsplit=1)[0]          # keep the author-year part of the heading
            title = re.sub(r"\s*\(`?10\.\d{4}.*$", "", title).strip()        # drop a DOI in brackets and what follows it
            citation = title[:110]
        if not doi:
            m = re.search(r"10\.\d{4,9}/[^\s`)\]]+", text[:1500])
            doi = m.group(0).rstrip(".,;") if m else ""
        doi = doi.strip("`* ").lower()
        citation = citation.replace("|", "/")
        if not re.search(r"(?:19|20)\d{2}", citation):
            citation = f"{citation} {year}"
        rows.append([key, citation, doi, take_line(text), [f.name]])
    # papers the model uses that have no dossier yet: rows from the registry, so the list is complete
    have = {r[0] for r in rows}
    used_keys = {k for c in counts.values() for k, n in c.items() if n}
    for key in sorted(used_keys - have):
        a, y = key.split(" ")
        hit = next((p for p in reg if p.get("citation") and re.search(rf"(?i){a}\b.*\b{y}\b", norm(p["citation"]) + " " + p["citation"])), None)
        if hit is None:
            hit = next((p for p in reg if norm(p["paper_id"]).startswith(a) and y in p["paper_id"]), None)
        if hit is None:
            continue
        rows.append([key, (hit.get("citation") or hit["paper_id"]).strip(), (hit.get("doi") or "").strip("`* ").lower(), "", []])
    # one row per paper: dossiers that share a DOI (a paper and its SI, or two readings) are merged
    merged = {}
    for key, citation, doi, take, names in rows:
        k = re.sub(r"[^0-9a-z]", "", doi) or names[0]
        if k in merged:
            m = merged[k]
            m[4] += names
            m[3] = m[3] or take
            m[1] = m[1] if len(m[1]) <= len(citation) else citation
        else:
            merged[k] = [key, citation, doi, take, names]
    # counts are keyed by author and year; a second paper by the same author in the same year cannot be told apart, so
    # the counts go to the first row of each author-year group and the others say so
    seen_keys = set()
    out_rows = []
    import unicodedata

    for key, citation, doi, take, names in sorted(merged.values(), key=lambda r: unicodedata.normalize("NFKD", r[1]).encode("ascii", "ignore").decode().lower()):
        if key in seen_keys:
            used = {k: 0 for k in counts}
            take = take or "(same author and year as the row above; any counts are listed there)"
        else:
            used = {k: c.get(key, 0) for k, c in counts.items()}
        seen_keys.add(key)
        out_rows.append((key, citation, doi, used, take, names))
    rows = out_rows
    out = ["# Sources: every paper the model uses, and what it takes from each", "",
           "*Generated by `scripts/generators/build_sources_table.py` from the paper registry, the extraction dossiers, "
           "the parameter registries, the frozen fit generators and the claims panel. One row per dossier, alphabetical. "
           "Counts: measured rate constants whose source is the paper; rows a calibration was tuned on; benchmark pots the "
           "model is scored against; directional claims it is scored on. A blank row means the paper was read and recorded "
           "but nothing quantitative was taken; the dossier says why. Papers with counts but no dossier were used through the paper registry "
           "and the frozen calibration scripts. Missing DOIs are papers with none printed.*", "",
           "| paper | DOI | constants | fit rows | pots | claims | what the repository takes | dossier |", "|---|---|---:|---:|---:|---:|---|---|"]
    for key, cit, doi, used, take, names in rows:
        doi_md = f"[{doi}](https://doi.org/{doi})" if doi else ""
        links = ", ".join(f"[{n}](../../data/lit/extraction_dossiers/{n})" for n in names) or "no dossier (registry entry only)"
        out.append(f"| {cit} | {doi_md} | {used['constants'] or ''} | {used['fit rows'] or ''} | {used['benchmark pots'] or ''} | "
                   f"{used['claims'] or ''} | {take} | {links} |")
    OUT.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(rows)} papers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
