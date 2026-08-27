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
  data/campaigns/shareable_meaty_screen.yml \
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

Read [architecture.md](../architecture.md).

That document explains the trust tiers: what the tool can support today, what remains directional only, and where wet-lab confirmation is still mandatory.
