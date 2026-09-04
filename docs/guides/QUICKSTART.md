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

## Five minutes with the tool: what the output looks like

The template compares cysteine + ribose (arm A) with cysteine + glucose (arm B), both 10 mM at
140 °C for 30 min. Trimmed real output:

```
  COMPARE [KINETIC CORE]   A = cysteine_ribose   vs   B = cysteine_glucose

  arm B:
  ENVELOPE: IN_ENVELOPE_EXTRAPOLATED   lane: sulfur
    ~ declared extrapolation -- HEXOSE ENTRY UNIDENTIFIED (FFT, MFT): the only route from a hexose
    ~ to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary
    ~ measurement identifies ...

  Axes this comparison moves: sugar_identity
  Governing reliability:      do-not-use (4/8)   [directional panel, independent claims]

  compound                                        A/B direction      resolved
  2-acetylthiazole                            0.0251x higher_in_cysteine_glucose yes
  2-furfurylthiol (FFT)                  B unidentified undefined      n/a -- undefined
  2-methyl-3-furanthiol (MFT)            B unidentified undefined      n/a -- undefined
  furfural                                     0.923x higher_in_cysteine_glucose no -- inside band
```

How to read it: the **envelope declaration** says what the engine could and could not represent
(here, the glucose arm's thiol route is declared unidentified, so no ratio is claimed for the two
thiols); the **reliability line** says how well the model's direction sense is measured on the axis
the two arms differ along (sugar identity: 4 of 8 independent claims, "do-not-use"); the **table**
gives a ratio only where both arms are predictions and marks it resolved only when it clears the
same-sample dispersion band. A ratio you can act on is one with a resolved "yes" *and* a
reliability of "caution" or better.

## One formulation: `maillard predict`

```bash
python scripts/maillard.py predict my_comparison.yml --system a
```

Absolute µg/L (= ppb in water) per compound, **always with its reliability interval** and, where
an odour threshold exists, the odour-activity ratio. Read the interval before the point: an
absolute out of sample lands within 3x of the measurement on a small minority of panel rows (the
current count is the first line of `python scripts/maillard.py --help` and of the README model card).

## What does the model know about a compound? `maillard explain`

```bash
python scripts/maillard.py explain 2-methyl-3-furanthiol
python scripts/maillard.py explain hexanal          # the lipid lane
python scripts/maillard.py explain 2-pentylfuran    # a refusal, with the reason
```

## Which measurement would teach the model most? `maillard wishlist` and `maillard rank`

```bash
python scripts/maillard.py wishlist            # what to measure next, and what it would unlock
python scripts/maillard.py rank --top 10       # where the envelope misses most (alias of rank-experiments)
python scripts/maillard.py rank --matrix pea_iso --json
```

`wishlist` prints `results/validation/data_wishlist.md`, generated from the scorecard, the slice
profile, the directional scorecard and the value-of-information ranking: which fitted constants the
evidence does not pin and the one measurement that would identify each, which rows the engine
answers but declares not evaluable, what no lane represents, which directional axes are thin, and
what each measurement would let you predict. `rank` reads the core's Monte-Carlo envelope
(`results/validation/core_prediction_uncertainty.json`) and orders (benchmark, compound) rows by
how badly and how uncertainly the model misses them.

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
   absolute is *do-not-use* by rule; directional skill is measured per axis on the core's own
   claims panel and printed as the reliability column of every `compare`.
2. Read the **envelope declaration** the verb printed: `IN_ENVELOPE` with declared
   extrapolations is the normal state; `OUT_OF_ENVELOPE` means no number was emitted.
3. Check whether your system resembles a panel bundle
   ([`core_panel_scores.md`](../../results/validation/core_panel_scores.md)): the sulfur lane at
   145 °C / 20 min is the best-covered region; long low-temperature holds and protein matrices
   are the worst.

## Regenerating the evidence

```bash
./scripts/docker_maillard.sh core-scores      # results/validation/core_panel_scores.*   (~15 s)
./scripts/docker_maillard.sh core-envelope    # results/validation/core_prediction_uncertainty.* (~5-7 min natively)
./scripts/docker_maillard.sh model-card       # re-splices the README model card
./scripts/docker_maillard.sh gates            # the six CI gates
./scripts/docker_maillard.sh pytest tests/unit
./scripts/docker_maillard.sh pytest tests/scientific
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
| `./scripts/docker_maillard.sh wishlist` | Regenerate the data wishlist (what to measure next) |
| `./scripts/docker_maillard.sh deep-research-audit` | Literature backlog audit |
| `python scripts/maillard.py {compare,predict,explain,score,rank,wishlist}` | The front door |
