> **ARCHIVED 2026-09-03 (retirement step B5c).** Every command below drove the retired screening lane
> (`run_pipeline.py`, `optimize_formulation.py`, `run_campaign.py`, `ingest_results.py`) and no longer exists.
> The live quick start is `docs/guides/QUICKSTART.md`.

# Quick Start

This guide is for a scientist who wants to understand the tool quickly and get one trustworthy result out of it.

## Goal

In under 10 minutes, you should be able to:

- start the validated environment
- inspect the current benchmark envelope
- generate a graphical trust snapshot
- run one formulation prediction

## macOS Recommended Workflow

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validation-figures
```

Then open:

- [../../results/validation/benchmark_summary.md](../../results/validation/benchmark_summary.md)
- [../../results/validation/validation_overview.md](../../results/validation/validation_overview.md)

## Start here — `maillard compare`, not a single prediction

*Added 2026-08-28 (Wave S5).* Every out-of-sample measurement in this repository says the
absolute ppb numbers are wrong by 6x to 94x, and that the **ordinal** claims are right 21
times in 29. So the shortest useful first command is a comparison, not a prediction:

```bash
python scripts/maillard.py compare --template > my_comparison.yml   # edit the two arms
python scripts/maillard.py compare my_comparison.yml
```

It prints per-compound **A/B ratios**, the dominant (rate-limiting) pathway on each side, and
a reliability tag per comparison axis read live from the directional-accuracy panel — so a
pH sweep is labelled `do-not-use (4/7)` and a sugar swap `trust (8/8)`, on the page, at the
moment you ask. Three verbs:

| Verb | Question it answers |
| --- | --- |
| `compare` | "Which of these two formulations gives more X?" — ratios, CIs, per-axis reliability |
| `predict` | "What comes out of this one formulation?" — ranges not points, caveats inline |
| `rank-experiments` | "What should I measure next?" — the value-of-information queue |

`--absolute` adds the raw ppb columns to `compare`, with the caveat printed alongside them;
`--json` emits the machine-readable payload. The full designer below is unchanged and still
the right tool when you want a complete report or a grid screen.

## Two modes — pick the one that answers your question

`scripts/run_pipeline.py` has two distinct modes, and the difference matters:

| | **Forward mode** (no `--target`) | **Inverse design** (`--target TAG`) |
| --- | --- | --- |
| What it scores | *Exactly* the formulation you typed | The 15-entry grid in `data/formulation_grid.yml`, **plus** your formulation as one extra candidate if you supply precursors |
| Honours `--ratios`, `--time-minutes`, `--sugars`, … | Yes, all of them | Only for your own candidate — grid entries carry their own precursors and ratios |
| Answers | "What does *this* recipe do?" | "Which recipe should I try first?" |
| Report describes | Your formulation | The **winning candidate**, which is usually a grid entry, not yours |

If you want a report about the recipe you typed, use forward mode.

## First Prediction — forward mode (scores your formulation)

```bash
python scripts/run_pipeline.py \
  --sugars ribose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --time-minutes 45 \
  --protein-type pea_iso \
  --report \
  --output-dir results/quickstart_run
```

## Screening — inverse design (ranks candidates)

```bash
python scripts/run_pipeline.py \
  --sugars ribose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --time-minutes 45 \
  --protein-type pea_iso \
  --target meaty \
  --minimize beany \
  --report \
  --output-dir results/quickstart_screen
```

Your formulation is entered as the candidate `Your formulation (custom)` and ranked
against the grid. The CLI states which candidate the report describes; if a grid
entry wins, the report is about *that* entry, and §1 of `report.md` lists the
formulation that was actually evaluated (not your raw arguments).

> ### ⚠️ `--minimize beany` does **not** predict beaniness from your ingredients
>
> Added 2026-08-27 after a cold-start review found this command reads as off-note
> screening when it cannot do that. Be clear about what the flag does and does not do:
>
> - **`--minimize beany` only ranks off-note compounds that are already in the run.**
>   It re-weights the objective; it does not generate lipid off-notes.
> - **No fatty acid is registered as a precursor.** `Linoleic acid` is not a resolvable
>   precursor name — the resolver rejects it. `--lipids` accepts *the aldehydes
>   themselves* (`hexanal`, `nonanal`), so in forward mode you must supply the off-note
>   you wanted predicted. Supplying nothing gives **Off-Flavour Risk 0.00**, which means
>   "no off-note compound was in the system", not "this formulation is not beany".
> - **The lipid radical chain cannot bridge the gap either.** It enumerates to *zero*
>   steps from an unoxidised fatty acid plus O₂ (it needs a hydroperoxide seed), and the
>   network's `MAX_MW = 300 Da` prune removes 13-HPODE (312 Da) and the peroxyl radical
>   (311 Da), so even a seeded chain is cut before it can reach hexanal.
> - **The matrix lane's hexanal is fit-recovery, not prediction.** The matrix-only
>   hexanal number comes from observability factors that were back-solved from the very
>   benchmarks they are scored against (see the `fit-recovery` rows in
>   `results/validation/benchmark_summary.md`).
>
> Use the off-note lane to **rank formulations that already contain aldehydes**, never to
> ask "will my ingredients go beany". Registering a fatty-acid precursor is an open owner
> item; it would be dishonest to add one before the chain that consumes it works, because
> a registered precursor that yields nothing looks like a confident zero.

What you should look for:

- `predicted_ppb`: observable concentration estimate
- `confidence_metadata`: how strongly the result is benchmark-backed
- `projection_metadata`: matrix/headspace corrections and calibration provenance
- `provenance`: command, git revision, and scientific reference surface used to generate the artifact
- any validated-envelope warnings that say the run is speculative

## If You Only Want To Check Scientific Trust

Run:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
```

Interpretation:

- summary: exact benchmark-by-benchmark status
- validated-envelope: plain-language boundary of what is supported
- validation-overview figure: one-page view of reliability and current limitations

## If You Care About Matrix Systems Specifically

Run:

```bash
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-assertions
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh coverage-gaps
```

These commands tell you three different things:

- what matrix data exists
- whether the current matrix benchmarks pass their own ranking assertions
- what still blocks external scientific assessment

## If You Need A Shareable Package

For one formulation, use `--report` as shown above.

For a shareable multi-run campaign, use:

```bash
./scripts/docker_maillard.sh campaign \
  docs/examples/shareable_meaty_screen.yml \
  results/quickstart_campaign
```

This creates run-level reports plus campaign-level Markdown and JSON summaries.

## Ingesting Lab Results

`scripts/ingest_results.py` (or `./scripts/docker_maillard.sh ingest`) normalises a GC-MS
result file into the matrix intake contract. **Every flag below is required** — the values
come either from columns/metadata in your file or from the flag:

```bash
python scripts/ingest_results.py \
  --file path/to/my_results.csv \
  --protein-type pea_iso \
  --process-state extrusion_structured \
  --temp-c 105 \
  --ph 5.5 \
  --water-activity 0.85 \
  --time-min 45 \
  --source-reference "Internal run 2026-08, GC-MS/SIDA" \
  --precursor cysteine=15.0 \
  --precursor glucose=30.0
```

Required: `--file`, `--protein-type`, `--process-state`, `--temp-c`, `--ph`,
`--water-activity`, `--time-min`, `--source-reference`, and at least one `--precursor`.
Missing any of them fails fast with `Missing required ingest field: <name>`.
`--water-activity`, `--time-min` and `--source-reference` are the ones usually forgotten.
`--precursor` values are millimolar.

Without `--confirm` this is a **preview**: it writes preview and support-delta artifacts
into `results/ingest_previews/` and touches nothing under `results/validation/`. Adding
`--confirm` writes the canonical intake YAML into `results/validation/`. That is **one
YAML file** — it does not rebuild the benchmark panel and does not regenerate validation
artifacts; run the generators (`./scripts/docker_maillard.sh summary`) for that.

## Reproducibility

Runs are deterministic: the same inputs give the same numbers, and two in-process runs are
bit-identical. Two things legitimately differ between invocations and are excluded from the
scientific fingerprint:

- `timestamp` / `generated_at`, and the recorded `argv` / `output_directory`
- `provenance.input_fingerprint_sha256` is computed over the *scientific* inputs only, so
  the same run written to two different `--output-dir` values fingerprints identically

There is no `--seed` flag because no sampling happens in a single run; the Monte-Carlo
uncertainty generators take their own `--seed`. If you want to diff two report JSONs, ignore
the timestamp fields and compare the rest — they should match byte for byte.

## Before You Trust A Result

Read the **model card** in [README.md](../../README.md#model-card--the-validity-domain-generated-from-the-artifacts).
It is generated from the live validation artifacts and states, per claim type and system
class, whether the answer you are about to quote is `trust`, `caution` or `do-not-use` —
with the measurement behind each verdict. The one-line version: compare, do not quote
absolutes, and treat any pH or moisture direction as unsupported.

*(2026-08-28, Wave S5: this section used to point at `docs/architecture.md`, which described
a pre-audit trust surface — it opened with "High trust — use freely" at a time when 0 of 14
benchmarks are strict-ready. That file was folded into README's architecture section and
deleted rather than left to contradict the evidence.)*

---

## Appendix — Command reference

*Folded in 2026-08-28 (Wave S5) from `docs/reference/COMMAND_REFERENCE.md`, which was deleted.
It was accurate — all 32 subcommands verified against `scripts/docker_maillard.sh` — but it was
a third surface answering "how do I run this", alongside this file and README's Getting
Started. One fewer document, same content.*

## Environment Lifecycle

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh up` | Start or create the validated Docker container |
| `./scripts/docker_maillard.sh bootstrap` | Create the conda environment and apply required patches |
| `./scripts/docker_maillard.sh shell` | Open an interactive shell inside the validated environment |
| `./scripts/docker_maillard.sh status` | Check container and environment status |

## Core Validation Commands

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh summary` | Generate the benchmark summary artifact |
| `./scripts/docker_maillard.sh index` | Generate the benchmark index artifact |
| `./scripts/docker_maillard.sh validated-envelope` | Generate the plain-language domain-of-validity artifact |
| `./scripts/docker_maillard.sh validation-figures` | Generate graphical reliability and limitation figures |
| `./scripts/docker_maillard.sh thermo-gating` | Audit whether thermodynamic gating materially helps |

## Test Lanes

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh core` | Unit and integration correctness |
| `./scripts/docker_maillard.sh scientific-fast` | Fast scientific regression lane |
| `./scripts/docker_maillard.sh kinetics-validation` | Slower kinetics reference lane |
| `./scripts/docker_maillard.sh scientific` | Full scientific artifact generation plus regression lane |
| `./scripts/docker_maillard.sh qm-heavy` | QM and external-backend validation |

## Benchmark Inspection

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh targets data/benchmarks/<benchmark>.json` | Inspect target-level predictions for one benchmark |
| `./scripts/docker_maillard.sh targets-report` | Aggregate target-level report across supported benchmarks |
| `./scripts/docker_maillard.sh hofmann` | Diagnostic trace for the Hofmann sulfur benchmark |

## Matrix-Specific Commands

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh matrix-deltas` | Per-compound matrix benchmark deltas |
| `./scripts/docker_maillard.sh matrix-evidence` | External versus internal matrix evidence audit |
| `./scripts/docker_maillard.sh matrix-assertions` | Ranking assertions for current matrix benchmarks |
| `./scripts/docker_maillard.sh matrix-readiness` | Family-level readiness for matrix promotion |
| `./scripts/docker_maillard.sh matrix-branch-deltas main` | Compare current branch matrix outputs against main |
| `./scripts/docker_maillard.sh coverage-gaps` | Coverage-gap report for selective mechanistic expansion work |
| `./scripts/docker_maillard.sh literature-learning-loop` | Publish literature gaps, ready references, and recommended experiment requests in `results/validation/literature_learning_loop.{md,json}` plus `active_learning_requests.json` |
| `./scripts/docker_maillard.sh computational-gap-refinement-plan` | Generate the descriptive computational-gap refinement plan and manifests |
| `./scripts/docker_maillard.sh computational-gap-xtb` | Run xTB seed generation for the named computational-gap targets |
| `./scripts/docker_maillard.sh computational-gap-dft-preflight aa_ring_open_dicarbonyl` | Run the geometry-control preflight and write the target-specific DFT execution JSON without starting heavy QM |
| `./scripts/docker_maillard.sh computational-gap-dft aa_ring_open_dicarbonyl` | Run DFT refinement for one computational-gap target with preflight, phase tracking, and a target-specific execution JSON |
| `./scripts/docker_maillard.sh computational-gap-dft-ingest` | Build the DFT ingestion report for computational-gap refinement |
| `./scripts/docker_maillard.sh computational-gap-dft-promote` | Promote completed computational-gap DFT results into priors |
| `./scripts/docker_maillard.sh refinement-governance` | Generate the selective mechanistic refinement governance artifact |
| `./scripts/docker_maillard.sh campaign docs/examples/shareable_meaty_screen.yml` | Generate a review-ready campaign package with run-level and campaign-level artifacts |

## CLI Workflows Outside Docker

| Command | Purpose |
| --- | --- |
| `python scripts/run_pipeline.py --list-precursors` | List available precursors |
| `python scripts/run_pipeline.py --list-tags` | List sensory/optimization tags |
| `python scripts/run_pipeline.py ...` | Run one forward prediction |
| `python scripts/optimize_formulation.py ...` | Search a formulation space |
| `python scripts/run_campaign.py --names "A,B" --ph 5.5 --temp 105` | Generate a side-by-side comparison package via the campaign pipeline |
| `python scripts/run_campaign.py --spec docs/examples/shareable_meaty_screen.yml` | Generate a shareable multi-run campaign package |
| `python scripts/compare_sim_to_lit.py` | Compare the framework against literature benchmarks |

## Recommended Review Sequence

If you are reviewing the state of the repository rather than developing it, the shortest useful sequence is:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh literature-learning-loop
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh coverage-gaps
```

