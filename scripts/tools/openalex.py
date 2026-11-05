#!/usr/bin/env python
"""OpenAlex literature tooling (plan W1.4, 2026-11-05). Report-only: nothing here writes under ``data/``.

Three subcommands:

``verify``  Check every DOI in the paper registry (``data_paths.PAPER_REGISTRY``) against OpenAlex: does the DOI
            resolve, is the work retracted, does the registry's citation (first-author surname, year) match the
            record's authorships and publication year. Writes a dated report to ``results/validation`` (or ``--out``).
``sweep``   Forward (works that cite the seeds) and/or backward (works the seeds cite) citations of seed DOIs,
            de-duplicated against the registry, with title, year, venue, authors, institutions and abstract, to a
            JSON file for screening. A sweep is keyed to a case's decisive question (STATUS_AND_PLAN sec. 7 item 2).
``labs``    Institutions per DOI (authorships), as a FIRST PASS for laboratory assignment (H0 / H1); a person
            confirms group boundaries (a group that moved institution is still one laboratory).

Why: on 2026-11-05 two sweeps found Vauthey 2000 (a candidate second laboratory for absolute thiol levels) that a
topic search had filed under "microemulsion", and a title-only row led to an invented first author ("Zhou" for
Lopez 2019) -- ``verify`` and author-from-record discipline catch that class of error.

Etiquette: anonymous OpenAlex access, one request per ``--delay`` seconds (default 1.1), retries with back-off on
HTTP 429 / 5xx. No account or e-mail is sent. Network access only; the unit tests use recorded fixtures.

Usage (through the wrapper):
    ./scripts/docker_maillard.sh run "python scripts/tools/openalex.py verify"
    ./scripts/docker_maillard.sh run "python scripts/tools/openalex.py sweep --seed 10.1021/jf991254a --direction forward --out /tmp/x.json"
    ./scripts/docker_maillard.sh run "python scripts/tools/openalex.py labs --doi 10.1021/jf991254a"
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_access, data_paths  # noqa: E402

API = "https://api.openalex.org"
WORK_FIELDS = ("id,doi,display_name,publication_year,is_retracted,authorships,primary_location,open_access,"
               "abstract_inverted_index,referenced_works,cited_by_count")
REPORT_STEM = "literature_metadata_check"

Fetch = Callable[[str], Optional[Dict[str, Any]]]


# ----------------------------------------------------------------------------- network

def make_fetch(delay: float = 1.1, retries: int = 4, timeout: float = 60.0) -> Fetch:
    """A polite GET-JSON function: sleeps ``delay`` s before each request; returns None on 404 or after retries."""
    def fetch(url: str) -> Optional[Dict[str, Any]]:
        for attempt in range(retries):
            time.sleep(delay)
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "maillard-literature-tool"})
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    return json.load(resp)
            except urllib.error.HTTPError as err:
                if err.code == 404:
                    return None
                if err.code in (429, 500, 502, 503, 504):
                    time.sleep(5.0 * (attempt + 1))
                    continue
                raise
            except (urllib.error.URLError, TimeoutError):
                time.sleep(5.0 * (attempt + 1))
        return None
    return fetch


def work_by_doi(doi: str, fetch: Fetch) -> Optional[Dict[str, Any]]:
    return fetch(f"{API}/works/doi:{urllib.parse.quote(doi, safe='/')}?select={WORK_FIELDS}")


def citing_works(work_id: str, fetch: Fetch) -> List[Dict[str, Any]]:
    out, cursor = [], "*"
    short = work_id.rsplit("/", 1)[-1]
    while cursor:
        page = fetch(f"{API}/works?filter=cites:{short}&per-page=200&cursor={cursor}&select={WORK_FIELDS}")
        if page is None:
            raise RuntimeError(f"OpenAlex page failed for cites:{short}; the sweep is incomplete")
        out.extend(page.get("results", []))
        cursor = (page.get("meta") or {}).get("next_cursor")
    return out


def works_by_ids(ids: Sequence[str], fetch: Fetch) -> List[Dict[str, Any]]:
    out = []
    shorts = [i.rsplit("/", 1)[-1] for i in ids]
    for k in range(0, len(shorts), 50):
        page = fetch(f"{API}/works?filter=openalex:{'|'.join(shorts[k:k + 50])}&per-page=50&select={WORK_FIELDS}")
        if page is None:
            raise RuntimeError("OpenAlex page failed for referenced works; the sweep is incomplete")
        out.extend(page.get("results", []))
    return out


# ----------------------------------------------------------------------------- pure helpers

def clean_doi(doi: Optional[str]) -> str:
    d = (doi or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi.org/", "doi:"):
        if d.startswith(prefix):
            d = d[len(prefix):]
    return d.rstrip(".")


def _fold(text: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFKD", text) if not unicodedata.combining(c)).lower()


def parse_citation(citation: Optional[str]) -> Dict[str, Any]:
    """'van Boekel (2001)', 'Ma et al. (2024)', 'Hofmann & Schieberle 2000' -> {'surname': ..., 'year': ...}."""
    if not citation:
        return {"surname": None, "year": None}
    year = re.search(r"(19|20)\d{2}", citation)
    head = re.split(r"\(|,| et al| and | & |\d", citation, maxsplit=1)[0].strip()
    return {"surname": _fold(head) or None, "year": int(year.group(0)) if year else None}


def record_summary(work: Dict[str, Any]) -> Dict[str, Any]:
    authors = [((a.get("author") or {}).get("display_name") or "") for a in work.get("authorships", [])]
    insts = sorted({i.get("display_name") for a in work.get("authorships", []) for i in a.get("institutions", [])
                    if i.get("display_name")})
    ab = work.get("abstract_inverted_index") or {}
    abstract = " ".join(sorted(ab, key=lambda t: min(ab[t]))) if ab else ""
    venue = ((work.get("primary_location") or {}).get("source") or {}).get("display_name")
    return {"openalex_id": work.get("id"), "doi": clean_doi(work.get("doi")), "title": work.get("display_name"),
            "year": work.get("publication_year"), "venue": venue, "authors": authors, "institutions": insts,
            "is_retracted": bool(work.get("is_retracted")), "cited_by_count": work.get("cited_by_count"),
            "oa_url": (work.get("open_access") or {}).get("oa_url"), "abstract": abstract}


def check_entry(entry: Dict[str, Any], work: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Compare one registry row with its OpenAlex record. Flags are advisory: a person reads each."""
    doi = clean_doi(entry.get("doi"))
    cit = parse_citation(entry.get("citation"))
    row = {"paper_id": entry.get("paper_id"), "doi": doi, "citation": entry.get("citation"), "flags": []}
    if work is None:
        row["flags"].append("doi_not_found")
        return row
    s = record_summary(work)
    row.update({"openalex_title": s["title"], "openalex_year": s["year"],
                "openalex_first_author": s["authors"][0] if s["authors"] else None})
    if s["is_retracted"]:
        row["flags"].append("retracted")
    if cit["year"] and s["year"] and abs(int(cit["year"]) - int(s["year"])) > 1:
        row["flags"].append("year_mismatch")
    if cit["surname"] and s["authors"]:
        surnames = [_fold(a) for a in s["authors"]]
        if not any(cit["surname"] in a or a.endswith(cit["surname"].split()[-1]) for a in surnames):
            row["flags"].append("author_not_in_record")
    return row


def registry_entries(registry: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [p for p in registry.get("papers", []) if clean_doi(p.get("doi")).startswith("10.")]


def known_dois(registry: Dict[str, Any]) -> set:
    return {clean_doi(p.get("doi")) for p in registry.get("papers", []) if p.get("doi")}


# ----------------------------------------------------------------------------- subcommands

def run_verify(registry: Dict[str, Any], fetch: Fetch, limit: Optional[int] = None) -> Dict[str, Any]:
    entries = registry_entries(registry)[: limit or None]
    rows = [check_entry(e, work_by_doi(clean_doi(e["doi"]), fetch)) for e in entries]
    counts: Dict[str, int] = {}
    for r in rows:
        for f in r["flags"]:
            counts[f] = counts.get(f, 0) + 1
    return {"checked": len(rows), "flag_counts": counts, "flagged": [r for r in rows if r["flags"]],
            "clean": sum(1 for r in rows if not r["flags"])}


def run_sweep(seeds: Iterable[str], fetch: Fetch, known: set, direction: str = "forward") -> Dict[str, Any]:
    found: Dict[str, Dict[str, Any]] = {}
    seed_rows = []
    for doi in seeds:
        work = work_by_doi(clean_doi(doi), fetch)
        if work is None:
            seed_rows.append({"doi": clean_doi(doi), "status": "not_found"})
            continue
        seed_rows.append({"doi": clean_doi(doi), "openalex_id": work["id"], "cited_by_count": work.get("cited_by_count")})
        related = []
        if direction in ("forward", "both"):
            related += [("forward", w) for w in citing_works(work["id"], fetch)]
        if direction in ("backward", "both"):
            related += [("backward", w) for w in works_by_ids(work.get("referenced_works") or [], fetch)]
        for how, w in related:
            s = found.setdefault(w["id"], {**record_summary(w), "via": []})
            s["via"].append({"seed": clean_doi(doi), "direction": how})
    rows = sorted(found.values(), key=lambda r: (-(r["year"] or 0), r["title"] or ""))
    for r in rows:
        r["held"] = r["doi"] in known
    return {"seeds": seed_rows, "direction": direction, "works": rows,
            "n_works": len(rows), "n_held": sum(r["held"] for r in rows)}


def run_labs(dois: Iterable[str], fetch: Fetch) -> List[Dict[str, Any]]:
    out = []
    for doi in dois:
        w = work_by_doi(clean_doi(doi), fetch)
        out.append({"doi": clean_doi(doi), "found": w is not None,
                    "authorships": [] if w is None else [
                        {"author": (a.get("author") or {}).get("display_name"),
                         "institutions": [{"name": i.get("display_name"), "ror": i.get("ror"), "country": i.get("country_code")}
                                          for i in a.get("institutions", [])]}
                        for a in w.get("authorships", [])]})
    return out


def render_verify_markdown(payload: Dict[str, Any]) -> str:
    lines = [f"# Literature metadata check against OpenAlex ({payload['date']})", "",
             f"*Generated by `scripts/tools/openalex.py verify` from `{payload['registry']}` "
             f"(sha256 {payload['registry_sha256'][:12]}). Flags are advisory: a person reads each before any record "
             "changes; nothing under `data/` was written.*", "",
             f"Checked {payload['checked']} DOIs; {payload['clean']} without a flag.", "",
             "| flag | count |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in sorted(payload["flag_counts"].items())]
    lines += ["", "| paper_id | DOI | registry citation | OpenAlex first author / year | flags |", "|---|---|---|---|---|"]
    for r in payload["flagged"]:
        oa = f"{r.get('openalex_first_author') or '-'} / {r.get('openalex_year') or '-'}"
        lines.append(f"| {r['paper_id']} | {r['doi']} | {r.get('citation') or '-'} | {oa} | {', '.join(r['flags'])} |")
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("verify")
    v.add_argument("--limit", type=int)
    v.add_argument("--out", type=Path, help="directory for the report (default results/validation)")
    s = sub.add_parser("sweep")
    s.add_argument("--seed", action="append", required=True)
    s.add_argument("--direction", choices=("forward", "backward", "both"), default="forward")
    s.add_argument("--out", type=Path, required=True, help="JSON output path (use a scratch location)")
    lab = sub.add_parser("labs")
    lab.add_argument("--doi", action="append", required=True)
    lab.add_argument("--out", type=Path)
    for p in (v, s, lab):
        p.add_argument("--delay", type=float, default=1.1)
    a = ap.parse_args(argv)
    fetch = make_fetch(delay=a.delay)
    registry = data_access.load_yaml(data_paths.PAPER_REGISTRY)
    if a.cmd == "verify":
        out_dir = a.out or data_paths.VALIDATION_DIR
        payload = run_verify(registry, fetch, a.limit)
        payload.update({"date": _dt.date.today().isoformat(), "registry": data_paths.rel(data_paths.PAPER_REGISTRY),
                        "registry_sha256": hashlib.sha256(data_paths.PAPER_REGISTRY.read_bytes()).hexdigest(),
                        "generated_by": "scripts/tools/openalex.py verify", "source": "OpenAlex (api.openalex.org)"})
        stem = out_dir / f"{REPORT_STEM}_{payload['date']}"
        stem.with_suffix(".json").write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
        stem.with_suffix(".md").write_text(render_verify_markdown(payload))
        print(f"verify: {payload['checked']} checked, {payload['clean']} clean, flags {payload['flag_counts']} -> {stem}.{{json,md}}")
    elif a.cmd == "sweep":
        if str(a.out.resolve()).startswith(str(data_paths.DATA_ROOT.resolve())):
            ap.error("sweep output must not go under data/")
        payload = run_sweep(a.seed, fetch, known_dois(registry), a.direction)
        a.out.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
        print(f"sweep: {payload['n_works']} works ({payload['n_held']} already held) -> {a.out}")
    else:
        rows = run_labs(a.doi, fetch)
        text = json.dumps(rows, indent=1, ensure_ascii=False)
        if a.out:
            a.out.write_text(text + "\n")
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
