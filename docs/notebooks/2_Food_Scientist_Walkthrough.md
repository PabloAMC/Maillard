# Screening a Meaty Fix for a Pea-Protein Burger

Companion document for `2_Food_Scientist_Walkthrough.ipynb`.

This walkthrough is meant to feel like product usage, not framework internals. The notebook now calls the public API directly, keeps notebook-only logic to a minimum, and produces one main artifact package for review.

## Scenario

We compare two precursor strategies for a `pea_iso` burger prototype at `130 C`, `20 min`, `pH 6.5`, and `aw 0.95`:

- **Hypothesis A:** glucose + glycine + residual hexanal
- **Hypothesis B:** ribose + cysteine + residual hexanal

The question is not “which one is chemically interesting?” The question is “which one should go into the next wet-lab round?”

## Trust Boundary

The walkthrough starts with a domain-of-validity check. The current model warns that:

- `pea_iso` is still outside the strict free-precursor quantitative proof surface;
- `glycine` and `hexanal` do not carry PRIMARY quantitative validation in this exact benchmark sense.

That means the notebook should be read as a directional screening tool for formulation prioritization.

## What Changed In The Notebook

Three changes make the notebook easier to understand:

1. It no longer defines notebook-local helper functions just to interpret results.
2. It no longer surfaces two separate per-run `report.md` outputs as the main takeaway.
3. The comparison table now shows only columns that actually separate the candidates.

## Head-To-Head Result

The revised review table focuses on meaty signal generation rather than constant-value columns:

| Candidate | Meaty Score | Detected Meaty Targets | FFT proxy (ppb) | MFT proxy (ppb) | 2-pentyl-4-methylthiazole proxy (ppb) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Hypothesis A (Glucose+Glycine) | 0.00 | 0 | 0.000 | 0.000 | 0.000 |
| Hypothesis B (Ribose+Cysteine) | 32.74 | 4 | 0.390 | 0.162 | 0.023 |

This table is easier to trust because it reflects the real decision signal in the model.

## Why The Earlier Table Was Confusing

The previous version included columns such as `Beany Score`, residual off-note proxies, and trapping values that were constant in this specific screen. That created a visually noisy table full of `0.0` and `100.0`, which made the example feel unconvincing even though the core chemistry separation was real.

The cleaner interpretation is:

- A fails because it does not generate the sulfur-enabled meaty package.
- B succeeds because FFT, MFT, hydrogen sulfide, and related targets actually appear.
- In this particular example, the positive meaty signal is much more informative than the negative beany penalty.

## Optimization Follow-Up

The notebook still runs a seeded Bayesian optimization around the winning chemistry family so the user can see what the next formulation iteration might look like.

Current seeded result:

- `Best objective value: 49.14`
- `Best trial target score: 49.49`
- `Best trial beany score: 0.00`
- `Best trial safety score: 0.00`

Representative best parameters:

- `sugar_conc: 0.056`
- `aa_conc_sulfur: 0.797`
- `ph: 3.936`
- `temp: 105.808`
- `aw: 0.863`
- `time_minutes: 76.123`

## Main Artifacts

To keep the output easy to follow, the walkthrough now points users to one main artifact directory:

- `results/notebook_walkthrough/campaign_package/comparison.md`
- `results/notebook_walkthrough/campaign_package/comparison.json`
- `results/notebook_walkthrough/campaign_package/campaign.md`
- `results/notebook_walkthrough/campaign_package/campaign.json`

## Decision

Advance **Hypothesis B (Ribose+Cysteine)** to the wet lab.

The notebook now reads more like a real user workflow: define the screen, check trust limits, compare candidates, optimize around the winner, and export one clean handoff package.
