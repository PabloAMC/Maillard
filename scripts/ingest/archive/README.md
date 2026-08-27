# Literature-curation script archive

**One-shot historical scripts. Do not re-run them.**

These files were written during the 2026 literature-curation campaigns to *land* deep-research
output into `data/lit/`. Every one of them is a write-once mutator: it loads a curated JSON
under `data/lit/`, splices in new entries, and writes the file back. Their outputs already live
in `data/lit/` (and are tracked in git), so the scripts have no remaining function — they are
kept only as provenance, so the curation decisions they encode are readable rather than
invisible.

They were moved here on 2026-08-26 as part of the audit remediation
(`tasks/audit_remediation.md`): they had been sitting untracked in the repo root and in
`scripts/`, which meant the literature-curation decisions in them existed nowhere in git
history.

## Caveats before you touch anything here

- **Re-running is destructive and mostly idempotency-free.** Most of these append entries
  unconditionally. Running one twice will duplicate registry rows in `data/lit/*.json`.
- **Their paths are now wrong.** They resolve the repo root relative to their own location:
  the `scratch_*.py` files use `Path(__file__).resolve().parent` (they lived in the repo root)
  and the `ingest_*.py` files use `.parents[1]` (they lived in `scripts/`). From
  `scripts/ingest/archive/` those resolve to the wrong directory. Any revival needs
  `.parents[3]` plus a review of the payloads.
- **Most of these files are not tracked in git.** They were moved with a plain `mv`. Whether to
  commit them is the repo owner's call; nothing here has been `git add`ed — *except*
  `scratch_check_citations.py` and `scratch_resolve_citations.py`, which were already tracked in
  the repo root and so were moved with `git mv` (2026-08-26) to keep their history.

## Catalogue

### Deep-research campaign ingestion

Landing scripts for the Gemini deep-research campaigns (`data/Gemini_Deep_Research/`). Each one
reads a campaign's curated references and writes them into the intake registry, the SLR
incorporation matrix, and the relevant prior/payload files.

| File | Purpose |
| --- | --- |
| `ingest_deep_research_campaigns.py` | Main campaign landing pass — intake registry, SLR matrix, computational priors, safety payloads, deep-research backlog. |
| `ingest_deep_research_campaigns_fgh.py` | Same, for campaigns **F/G/H** specifically. |
| `ingest_deep_research_campaign_i.py` | Same, for campaign **I**; also writes retention reference payloads. |
| `ingest_deep_research_priors.py` | Priors-only pass over deep-research output (computational priors + safety payloads + backlog). |

### Chemistry-family ingestion

Per-family batch landings. Each sets `status = "encoded"` in `benchmark_intake_registry.json`,
adds rows to `slr_incorporation_matrix.json`, and extends the prior/payload files its family
needs.

| File | Purpose |
| --- | --- |
| `ingest_family02.py` | Family 02 batch: intake registry, SLR matrix, computational priors, flavor payloads, backlog. |
| `ingest_family06_07.py` | Families 06 + 07 batch (adds safety payloads on top of the family-02 set). |
| `ingest_family06_07_backlog.py` | Families 06/07 backlog remainder; also writes process-state calibrations. |
| `ingest_family08_09_10_11.py` | Families 08–11 batch: priors, flavor payloads, process-state calibrations, backlog. |
| `ingest_batch_alternative_matrix.py` | `alternative_protein_matrix_scope` (27 refs) plus the remaining `amino_acid_sugar_core`, `glutathione_and_peptide_support`, `nucleotide_and_ribose_support`, `thiamine_fragmentation_support` and `fermentation_pretreatment` refs. Adds carbonyl-donor, sulfur-peptide, nucleotide, thiamine and shear-damage priors. |
| `ingest_batch_protein_damage_markers.py` | `protein_damage_markers` (25 refs) → CML / CEL / acrylamide / furosine / HCA / PhIP safety payloads and furosine-conversion + crosslink priors. |
| `ingest_batch_remaining_families.py` | Remaining families in one pass: melanoidin polymerization, ascorbic-acid Maillard, polyphenol amino capping, lipid-oxidation and lipid–Maillard crosstalk, phospholipid amine sink, off-note suppression. |
| `ingest_advanced_caps_damage_polymers.py` | Advanced capping / damage / polymer references across the full payload set (flavor, safety, process-state, priors, backlog). |
| `ingest_plant_protein_kinetics.py` | Plant-protein kinetics references → computational priors + safety payloads. |
| `ingest_volatile_sulfur_priors.py` | Volatile-sulfur priors → computational priors + retention reference payloads. |

### Citation and reference cleanup (one-off fixes)

| File | Purpose |
| --- | --- |
| `scratch_add_aliases.py` | Adds citation aliases (PMC/PMID/DOI variants) to registry, SLR matrix, priors and process-state calibrations so cross-references resolve. |
| `scratch_finalize_citations.py` | Rewrites a hand-curated list of `source` strings to their correct author-year form (e.g. `jafc_2019_ref21_...` → `Zha et al. (2019)`), recursively across the curated JSONs. |
| `scratch_fix_urugo.py` | Targeted text fix applied across `data/lit/*.json` **and** two of the ingest scripts here (a bad author token that had been propagated into the payloads). |
| `scratch_map_priors.py` | Read-only inspection: dumps `computational_priors.json` entry ids and their Ea/barrier/rate keys. Used to locate kinetic priors; writes nothing. |
| `scratch_parse_deep_research_i.py` | Read-only extraction: regex-pulls `"analyte"` JSON blocks out of `data/Gemini_Deep_Research/Gemini_deep_research_I.md` and test-parses them. Feeder for `ingest_deep_research_campaign_i.py`. |
| `scratch_check_citations.py` | Read-only regex validator for the citation strings in `benchmark_intake_registry.json`. **Tracked in git** — moved here with `git mv` on 2026-08-26, so its history (and the local modifications it carried) followed it. |
| `scratch_resolve_citations.py` | Resolves malformed citation strings in `benchmark_intake_registry.json` against Crossref over the network (`urllib` + rate-limit sleeps). **Tracked in git** — moved here with `git mv` on 2026-08-26. |

## Live tools left in place (NOT archived)

These stay in `scripts/` because they are re-runnable diagnostics, not one-shot mutators. They
are read-only over `data/lit/` and are still useful:

- `scripts/check_references.py` — cross-checks references in `benchmark_intake_registry.json`,
  `slr_incorporation_matrix.json` and `deep_research_backlog.json` against the deep-research
  source markdown in `data/Gemini_Deep_Research/`.
- `scripts/find_all_doi_gaps.py` — recursively scans the curated payload JSONs
  (intake registry, SLR matrix, flavor/safety/retention payloads, process-state calibrations)
  and reports references that lack a DOI. Backs the `[P]` "78% of backlog items lack DOIs"
  item in `tasks/audit_remediation.md`.

Also deliberately left in `scripts/`:

- `scripts/ingest_results.py` — tracked; the live GC-MS ingestion entry point wired into
  `scripts/docker_maillard.sh ingest`.
- `scripts/ingest_deep_research_markdown.py` — tracked; the live deep-research markdown reader.
- `scripts/sync_backlog.py` — untracked, but re-runnable: reconciles
  `deep_research_backlog.json` against `benchmark_intake_registry.json` rather than mutating
  curated payloads.
- `scripts/get_details_fgh.py` (+ `get_details_fgh_output.txt`) — untracked one-shot extraction
  helper for campaigns F/G/H whose output is checked in beside it; left in place because it is
  neither an `ingest_*` nor a `scratch_*` script. A candidate for this archive.
The two tracked `scratch_*.py` files that used to sit in the repo root
(`scratch_check_citations.py`, `scratch_resolve_citations.py`) were moved into this archive
with `git mv` on 2026-08-26, so their history followed them. They are catalogued above under
*Citation and reference cleanup*; the repo root is now clear of `scratch_*.py` scripts.
