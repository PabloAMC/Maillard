# ⛔ ARCHIVE OF MACHINE-GENERATED RESEARCH DUMPS — NOT A SOURCE

**Added 2026-08-27 (Wave T3 of the audit remediation; Wave T1 recommendation §7).**

## What is in this directory

Five files, all of them **LLM output**, none of them a paper:

| File | What it actually is |
|---|---|
| `Maillard_meat.md` | Gemini Deep Research digest |
| `Maillard_Plant_based.md` | Gemini Deep Research digest |
| `Kinetic and Thermodynamic Profiling of Aqueous Maillard and Advanced Sulfur Pathways_ An Exhaustive Analysis of Arrhenius Parameters.md` | Gemini Deep Research digest |
| `Elicit - Maillard Pathways in Plant-Based Cooking - Report.md` | Elicit report (LLM-summarised search results) |
| `Elicit - Maillard Reaction Pathways in Meat Cooking - Report.md` | Elicit report (LLM-summarised search results) |

They are **duplicates of the corpus under `data/Gemini_Deep_Research/`**, which carries a
strong warning README. These copies had none, and they sit under a `docs/` path, which reads
as authoritative. That is the defect this file closes.

## The rule — identical to `data/Gemini_Deep_Research/README.md`

> **No number, rate constant, activation energy, threshold, ratio, or citation may enter the
> model from a file in this directory without independent verification against the primary
> source.** "The deep-research report says so" is not provenance. If the primary source cannot
> be obtained, the value does not ship; if it ships anyway, it is marked
> `source_status: no_verifiable_source` and its citation field is left empty.

Additional cautions specific to these files:

* **Their citations have never been CrossRef-checked as a set.** The 2026-08-26 sweep of the
  live data layer found 61 % of DOI anchors defective (dead, topic-mismatched, or
  metadata-mismatched), several with a recognisable confabulation signature. Nothing has ever
  audited the identifiers *inside these documents*.
* **A number appearing here, near an author's name, is not evidence about that author's
  paper.** `scripts/trace_key_values.py` measures exactly that proximity and — until it was
  disarmed on 2026-08-27 — published the result as "Fully Verified". Read
  `results/validation/key_value_trace_report.md` and that script's docstring before treating
  any digest-adjacency as verification.

## Why these files are kept rather than deleted

They are the evidence trail. Several audit findings consist of showing that a shipped constant
traces to one of these documents and to nothing else; deleting them would destroy the
provenance record for the defect. Keep, warn, and never source from.

## Live dangling reference — `pathways.md` does not exist

`data/species/desirable_targets.yml` and `data/species/off_flavour_targets.yml` both declare
their `Source:` as, among these filenames, **`pathways.md`** — a file that exists nowhere in
this repository. Those two YAML files supply the odour thresholds that are the denominator of
every odour-activity value computed in `src/sensory.py`, so the header of a runtime-critical
data file names an LLM digest and a nonexistent document as its authority.

Wave T3 did **not** edit those two YAML headers: re-sourcing per-compound odour thresholds is a
retrieval job with a real risk of substituting plausible numbers for absent ones, which is the
defect this whole remediation exists to remove. Carried as **[P], owner decision** — see
`tasks/audit_remediation.md` § Wave T3 and Wave T1 finding T1-10. Until then, treat every
`odour_threshold_ug_per_kg` in those files as uncited unless it carries its own inline
`AUDIT`/source comment (hexanal and MFT do; most do not).

## Preferred long-term disposition

`git mv` these five files under `data/Gemini_Deep_Research/`, where the existing README already
covers them, and leave a pointer here. Not done in this wave because `data/species/*.yml` cite
them by bare filename and the reference-repair is a separate, owner-visible change.
