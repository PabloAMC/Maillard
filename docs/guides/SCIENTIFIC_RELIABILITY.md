# Scientific Reliability And Limits

This is the short, honest answer to the question: how much should a scientist trust Maillard today?

## Executive Summary

Maillard is currently strongest as a benchmark-aware formulation ranking tool inside a narrow validated envelope. It is not yet a universal quantitative predictor for all protein matrices or processing regimes.

## Current Validation Status

Based on the current Docker-generated artifacts:

| Benchmark family | Current status | What that means |
| --- | --- | --- |
| Free-precursor PRIMARY benchmarks | Strong | Four benchmarks are strict-ready and useful for quantitative regression inside the current envelope |
| Matrix off-flavour pea/soy benchmarks | Moderate | Executable and reproducible, but intentionally outside the strict gate |
| Matrix meaty-positive pea/soy candidates | Directional | Assertions pass, but they are internal reproducibility harnesses rather than external wet-lab evidence |
| New protein families or harsh process states | Low | Treat as exploratory until new benchmark data exists |

## Benchmarks You Can Trust Most Today

The strict-ready set currently includes:

- `acrylamide_asparagine_glucose_Parker2012`
- `cys_glucose_150C_Farmer1999`
- `cys_ribose_140C_Hofmann1998`
- `cys_ribose_150C_Mottram1994`

These are the most credible basis for:

- relative ranking across free-precursor formulations
- benchmark-aware model regression
- safety-aware formulation comparisons

## What Matrix Results Mean Today

For pea and soy, the repository now distinguishes three different things clearly:

- external off-flavour anchors
- internal meaty-positive reproducibility harnesses
- strict-gate-eligible free-precursor benchmarks

That distinction matters.

- `matrix_only` pea/soy benchmarks are useful executable intake/headspace checks.
- `matrix_precursor_augmented` pea/soy benchmarks are useful ranking harnesses.
- Neither one should be described as externally validated meaty-positive matrix prediction yet.

## Trust Tiers For Practical Use

### High Trust

Use the tool confidently for:

- free-precursor comparative formulation work
- identifying likely desirable versus adverse tradeoffs in the validated envelope
- checking whether a new intervention is directionally consistent with benchmark-backed chemistry

### Moderate Trust

Use the tool carefully for:

- pea/soy matrix comparisons when the output includes clear caveats
- process-state-aware ranking work where the question is comparative, not absolute
- deciding which experiment to run next rather than making a final quantitative claim

### Low Trust

Treat outputs as exploratory for:

- new protein families without nearby benchmarks
- extrusion-heavy or highly processed states without dedicated evidence
- peptide-bound and intact-protein systems
- absolute concentration claims in matrix systems without external measurements

## What The Tool Still Cannot Demonstrate

- external meaty-positive matrix assessment for pea or soy
- broad matrix strict-gate promotion
- robust quantitative transfer across all process states
- generalization to all alt-protein categories

## Can The Repo Generate Figures Showing Reliability And Limitations?

Yes.

Use:

```bash
./scripts/docker_maillard.sh validation-figures
```

This writes:

- [Versioned validation overview snapshot](../assets/validation_overview.png)
- [../../results/validation/validation_overview.md](../../results/validation/validation_overview.md)
- [../../results/validation/validation_overview.json](../../results/validation/validation_overview.json)

The figure summarizes:

- free-precursor benchmark ratio and error behavior
- benchmark status distribution
- matrix readiness by protein family
- current coverage gaps and missing scientific evidence

## Recommended Commands Before Presenting Results Externally

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-evidence
./scripts/docker_maillard.sh matrix-readiness
```

## Bottom Line

Maillard is already useful as a scientist-facing decision support tool when you keep the claims aligned with the validated envelope. It is strongest at ranking, prioritization, and transparent boundary-setting. It is not yet justified as a broad absolute predictor for external matrix meaty-positive performance.
