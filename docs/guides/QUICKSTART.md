# Quick Start

*Rewritten 2026-09-03 (retirement step B5c) for the one-engine tool. The previous quick start,
which drove the retired screening lane (`run_pipeline.py`, `optimize_formulation.py`,
`run_campaign.py`, `ingest_results.py`), is kept at
[`docs/history/QUICKSTART_legacy_lane_2026-09-03.md`](../history/QUICKSTART_legacy_lane_2026-09-03.md).*

## Goal

Decide, before a GC-MS run, which of two formulations or process settings is more likely to
give the aroma you want — and see exactly why the model refuses to answer when it cannot.

## Boot the environment

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
```

Every command below runs inside the container, either through
`./scripts/docker_maillard.sh run "<command>"` or from `./scripts/docker_maillard.sh shell`.

## Start here: `maillard compare`

```bash
python scripts/maillard.py compare --template > my_comparison.yml   # edit the two arms
python scripts/maillard.py compare my_comparison.yml
python scripts/maillard.py compare my_comparison.yml --report compare.html
```

A spec names two arms (`a`, `b`), each with `precursors` (name → mM), `temp_C`, `time_min`,
`ph`, `aw` and optionally `matrix` and `targets`. The output leads with the **ratio A/B per
compound** — the quantity the model's systematic scale error cancels out of — and prints each
arm's **envelope declaration**: which lane answered, which assumptions were extrapolated, and
any compound that was refused with its reason. `--json` gives the payload; `--report` writes a
self-contained HTML page with the OAV chart, refusal cards and every declared assumption.

## One formulation: `maillard predict`

```bash
python scripts/maillard.py predict my_comparison.yml --system a
```

Absolute µg/L (= ppb in water) per compound, **always with its reliability interval** and, where
an odour threshold exists, the odour-activity ratio. Read the interval before the point: an
absolute out of sample lands within 3x of the measurement on 5 of 48 panel rows.

## What does the model know about a compound? `maillard explain`

```bash
python scripts/maillard.py explain 2-methyl-3-furanthiol
python scripts/maillard.py explain hexanal          # the lipid lane
python scripts/maillard.py explain 2-pentylfuran    # a refusal, with the reason
```

## Which measurement would teach the model most? `maillard rank`

```bash
python scripts/maillard.py rank --top 10
python scripts/maillard.py rank --matrix pea_iso --json
```

Reads the core's Monte-Carlo envelope (`results/validation/core_prediction_uncertainty.json`)
and orders (benchmark, compound) rows by how badly and how uncertainly the model misses them.

## Your own measurements: `maillard score`

```bash
python scripts/maillard.py score --template > my_measurements.yml   # edit: precursors, process, measured ppb
python scripts/maillard.py score my_measurements.yml
```

Each measured compound is scored the way the panel scores a benchmark (fold error, the 3x band, the
reliability interval, or a named refusal), and a bundle-shaped record lands under `results/user/` with
your provenance. Nothing is refitted: calibration on new data is a new pre-registered fit wave.

## Before you trust a result

1. Read the **model card** in [README.md](../../README.md#how-well-calibrated-is-it): every
   absolute is *do-not-use* by rule; directional skill is not yet measured on the core.
2. Read the **envelope declaration** the verb printed: `IN_ENVELOPE` with declared
   extrapolations is the normal state; `OUT_OF_ENVELOPE` means no number was emitted.
3. Check whether your system resembles a panel bundle
   ([`core_panel_scores.md`](../../results/validation/core_panel_scores.md)): the sulfur lane at
   145 °C / 20 min is the best-covered region; long low-temperature holds and protein matrices
   are the worst.

## Regenerating the evidence

```bash
./scripts/docker_maillard.sh core-scores      # results/validation/core_panel_scores.*   (~15 s)
./scripts/docker_maillard.sh core-envelope    # results/validation/core_prediction_uncertainty.* (~40 min)
./scripts/docker_maillard.sh model-card       # re-splices the README model card
./scripts/docker_maillard.sh gates            # the six CI gates
./scripts/docker_maillard.sh pytest tests/unit tests/scripts
./scripts/docker_maillard.sh pytest tests/integration tests/scientific
```

## Command reference

| Command | Purpose |
| --- | --- |
| `./scripts/docker_maillard.sh up` / `bootstrap` / `shell` / `status` | Container lifecycle |
| `./scripts/docker_maillard.sh run "<cmd>"` | Any command inside the environment |
| `./scripts/docker_maillard.sh pytest [args]` | pytest inside the environment |
| `./scripts/docker_maillard.sh core-scores` | The core's panel scorecard |
| `./scripts/docker_maillard.sh core-envelope [--n-samples N --workers W]` | The core's Monte-Carlo envelope |
| `./scripts/docker_maillard.sh model-card` | Regenerate the README model card |
| `./scripts/docker_maillard.sh gates` | Run the six CI gates (incl. artifact freshness) |
| `./scripts/docker_maillard.sh data-readme` | Regenerate `data/README.md` |
| `./scripts/docker_maillard.sh experiment-value-ranking` | Rank experiments by value of information |
| `./scripts/docker_maillard.sh deep-research-audit` | Literature backlog audit |
| `python scripts/maillard.py {compare,predict,explain,rank,score}` | The front door |
