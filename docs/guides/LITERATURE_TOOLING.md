# Literature tooling: `scripts/tools/openalex.py` (plan W1.4, 2026-11-05)

A report-only OpenAlex client. Nothing it does writes under `data/`. It uses anonymous access at about one request per
second (`--delay`), backs off on HTTP 429 / 5xx, and sends no account or e-mail. The unit tests run on recorded
fixtures and need no network.

| subcommand | what it does | when to run it |
|---|---|---|
| `verify` | Every DOI in `data/keys/papers.yml` against OpenAlex: DOI resolves, `is_retracted`, registry citation (first-author surname, year) against the record. Writes `results/validation/literature_metadata_check_<date>.{json,md}`. Flags are advisory: a person reads each before any record changes. | once now; before each benchmark re-score; after a large dossier batch |
| `sweep` | Forward (citing) and/or backward (referenced) works of seed DOIs, de-duplicated against the registry, with authors, institutions and abstracts, to a JSON file **outside `data/`** (scratch) for screening. | before a case's items 1 - 5 freeze (STATUS_AND_PLAN sec. 7 item 2), keyed to that case's decisive criterion |
| `labs` | Authorship institutions (name, ROR, country) per DOI. | a first pass for laboratory ids (H0 / H1); a person confirms group boundaries, e.g. a group that moved institution is one laboratory |

**What earlier searches already did.** Broad forward chaining was run on 2026-10-28
([LITERATURE_SEARCH_2026-10_FORWARD.md](LITERATURE_SEARCH_2026-10_FORWARD.md): 2 432 citing works ranked against the
sulfur-lane gaps; "almost no modern kinetics the older rounds missed"), and Round 3 used OpenAlex records too. The value
of `sweep` now is narrower and sharper: run it against **one question** with its own decisive criterion. On
2026-11-05 a second-laboratory sweep (absolute MFT / FFT outside the DFA network) surfaced Vauthey 2000 (Nestec), which
Round 3 had read for mechanism and skipped.

**Discipline.** An author name or file key is taken from the record (authorships or the PDF's first page), never
inferred from a title row (lesson 2026-11-05). Screening is by title and abstract until a case's items 1 - 5 freeze;
abstracts are missing for many publisher records, so a title-only screen is recorded as such.

Usage, through the wrapper:

```bash
./scripts/docker_maillard.sh run "python scripts/tools/openalex.py verify"
./scripts/docker_maillard.sh run "python scripts/tools/openalex.py sweep --seed 10.1021/jf991254a --direction both --out /path/to/scratch/sweep.json"
./scripts/docker_maillard.sh run "python scripts/tools/openalex.py labs --doi 10.1021/jf991254a"
```
