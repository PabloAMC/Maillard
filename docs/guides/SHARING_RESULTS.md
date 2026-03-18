# Sharing Results

This guide is for the moment when you want to show Maillard outputs to alternative-protein scientists, reviewers, or collaborators without relying on live terminal context.

## What To Share

Use one of these three artifact levels:

- single-run report when you want to explain one formulation deeply
- comparison report when you want to compare a small set of named formulations
- campaign report when you want a review-ready package with multiple runs, a leaderboard, and reproducible provenance

## 1. Single-Run Report

Use this when you already know the candidate formulation.

```bash
python scripts/run_pipeline.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --ratios ribose:0.5,glucose:0.2,cysteine:0.2,leucine:0.1 \
  --ph 5.5 \
  --temp 105 \
  --time-minutes 45 \
  --protein-type pea_iso \
  --target meaty \
  --minimize beany \
  --report \
  --output-dir results/share/run_pea_meaty
```

The output directory contains:

- `report.md` for scientists and reviewers
- `report.json` for machine-readable downstream use
- a provenance block with branch, commit, dirty-state flag, command, and scientific reference surface

## 2. Small Comparison Report

Use this when the question is comparative.

```bash
python scripts/compare_formulations.py \
  --names "Cysteine Enrichment (Basic),Premium Meaty Mix,Soy-Specific Masking" \
  --ph 5.5 \
  --temp 105 \
  --target-tag meaty \
  --minimize-tag beany \
  --output-dir results/share/comparison_meaty
```

This creates `comparison.md` and `comparison.json` with the same provenance surface.

## 3. Review-Ready Campaign Package

Use this when you need a package you can hand to scientists without extra explanation.

```bash
./scripts/docker_maillard.sh campaign \
  data/campaigns/shareable_meaty_screen.yml \
  results/share/campaign_meaty
```

The campaign output contains:

- one report bundle per run in `runs/`
- `comparison.md` and `comparison.json`
- `campaign.md` and `campaign.json`
- provenance for the campaign plus the scientific reference files that define the current trust surface

## Minimum Sharing Checklist

Before sending artifacts externally, regenerate these repository-level trust surfaces:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-readiness
```

Share those alongside the run or campaign package whenever the claim depends on matrix behavior.

## What The Provenance Block Means

Every shareable artifact now records:

- the generating command
- git branch and commit
- whether the repository was dirty
- a stable hash of the input bundle
- the key scientific-reference files used to interpret the result

This is intended to make review and reproduction easier, not to imply that every result is externally benchmark-validated.

## Claim Discipline

When sharing with scientists:

- use free-precursor outputs for quantitative claims
- use pea and soy matrix outputs for directional prioritization unless external benchmark evidence exists
- include the validated-envelope artifact whenever you present matrix-heavy results