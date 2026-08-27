# ⚠️ NOT PRIMARY LITERATURE — LLM-GENERATED RESEARCH OUTPUT

**Every file in this directory and its `raw/` subdirectory is machine-generated text produced by
a large language model ("Gemini Deep Research", plus a small number of Elicit exports). None of
it is a peer-reviewed source. None of it has been through editorial review. Nothing in it is a
measurement.**

These are *literature digests*: an LLM was asked to survey a topic and it returned prose,
tables of "kinetic parameters", and lists of citations. That output was then used as a starting
point for building this repository's literature registries. The audit of 2026-08-26/27 found
that a large fraction of the citations in those registries are defective, and this directory is
where most of them came from.

## The measured contamination

From `results/validation/citation_verification_ledger.md` — every DOI used as an *anchor*
anywhere in the repo's registries was resolved against CrossRef and compared with the claim
made at the anchor site:

| Class | Unique DOIs repo-wide | …of which also appear verbatim in **this directory** |
|---|---:|---:|
| `MATCH` (supports the claim) | 88 | 72 |
| `METADATA-MISMATCH` (right paper, wrong bibliographic claim) | 47 | 39 |
| `TOPIC-MISMATCH` (live DOI, **different subject**) | 45 | 39 |
| `DEAD` (does not resolve) | 45 | 39 |

**117 of the 137 defective anchors (85%) are present in these files.** That is not proof that
each one originated here, but this directory is the largest single reservoir of them in the
tree, and several carry the unmistakable confabulation signature of generated identifiers —
`10.1021/acs.jafc.de_leyn_2019`, `10.1016/j.foodchem.2015.00000`,
`10.1016/j.foodres.2025.001279` — DOIs that were never issued by any registrant.

There are **431 distinct DOI-shaped tokens** in this directory. Only **189** of them have ever
been checked against CrossRef (they are the ones that made it into a registry and so entered
the ledger). **The remaining 242 have never been verified by anything.** Their status is
unknown, and given the base rate above (117 of the 189 checked, 62%), most of them should be
expected to be defective too.

*Method, so these counts are reproducible rather than asserted:* tokens are extracted from
every `*.md` under this directory (including `raw/`) with `citation_gate.TEXT_DOI_RE`, trimmed
with `citation_gate._trim_text_doi`, and **lower-cased** before being matched against the
`doi` field of each entry in `results/validation/citation_verification_ledger.json`. The
case-folding matters: DOI suffixes are case-insensitive for resolution, and matching them
case-sensitively under-counts the overlap by one in three of the four classes above.
(2026-08-27 Wave J1 recount; an earlier draft of this file reported 430 / 188 / 116 / DEAD 38
from a case-sensitive comparison.)

## The rule

> **No number, rate constant, activation energy, threshold, ratio, or citation may enter the
> model from a file in this directory without independent verification against the primary
> source.** "The deep-research report says so" is not provenance. If the primary source cannot
> be obtained, the value does not ship; if it ships anyway, it is marked
> `source_status: no_verifiable_source` and its citation field is left empty.

`scripts/ci/citation_gate.py` enforces the mechanical half of this (DOI grammar, confabulation
signatures, status coherence, repair-record completeness) and runs blocking in CI. It cannot
detect a live DOI that resolves to the wrong paper — the largest defect class above — so the
gate is a floor, not a substitute for reading the paper.

## Why this is tracked at all

Because deleting it would destroy the evidence. This directory is the *provenance* of the
contamination the audit found, and a reviewer checking whether a given number is anchored has
to be able to follow the chain back to where it actually came from. It is committed for
forensics, not for citation.

## Known runtime coupling — an open owner item

Two shipped source modules name files in this directory as the authority for constants that are
live in the model. Neither opens the file at runtime; both cite it as the provenance a reader is
meant to trust:

* `src/matrix_correction.py` (module docstring) — the **soy** matrix correction factors are
  attributed to `maillard_plant_based.md` and `elicit_plant_based_cooking.md`.
* `src/barrier_constants.py` (module docstring) — `maillard_meat.md` and
  `maillard_plant_based.md` are listed among the starting-point sources for `FAST_BARRIERS`.

So the answer to "is an LLM digest currently load-bearing for runtime numbers?" is **yes, at
the level of provenance**: those constants have no other citation behind them. Re-anchoring
them on primary sources — or relabelling them as unanchored fits — is tracked as a `[P]` item
in `tasks/audit_remediation.md` (Wave J1).

Separately, the `source_slr_file` / scope-registry pointers in `data/lit/family_ingestion_plan.json`
and `data/lit/chemistry_family_scope_registry.json` reference **15 paths in this directory that do
not exist** (`01_amino_acid_sugar.md`, `11_maillard_lipid_crosstalk.md`, …). Fourteen are stale
names from before the `slr_family_NN_<family_id>.md` rename and resolve by hand; one,
`00_framework_scaffold.md`, has no counterpart in the tree at all. Also logged as `[P]`.

## Inventory

| Path | What it is |
|---|---|
| `raw/Gemini_deep_research_*.md` | The unedited model outputs, as generated. |
| `slr_family_NN_*.md` | Per-chemistry-family digests, generated by `scripts/generators/generate_slr_family_reports.py` **from** the raw dumps — machine output derived from machine output. |
| `maillard_*.md`, `plant_protein_thermal_processing.md`, `kinetic_thermodynamic_profiling.md`, `Deep_Research_*.md`, `*_Literature_Search.md` | Topic surveys. |
| `elicit_*.md` | Elicit exports (also LLM-summarised). |

---

*Added 2026-08-27, audit-remediation Wave J1. Before this commit the whole directory was
excluded by `.gitignore` line 33 (`data/*`) and invisible to anyone reading the repository —
including the reviewers auditing its citations.*
