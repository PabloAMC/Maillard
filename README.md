# Maillard

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Docker Recommended](https://img.shields.io/badge/docker-recommended-blue.svg)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Maillard is a predictive computational framework designed to help alternative-protein scientists digitally explore and rank Maillard flavor chemistry *before* stepping into the wet lab.

By simulating the complex reactions between sugars and amino acids under specific conditions, it tells you what your current formulation will likely produce, what trade-offs exist (e.g., flavor vs. safety), and what you should try next.

---

## What This Repository Is For

Use Maillard when you want to answer questions like these before running a wet-lab experiment:

- Which precursor combination is most likely to generate meaty sulfur compounds under my process conditions?
- If I change pH, temperature, water activity, or reaction time, which flavour-active compounds move the most?
- Which candidate formulation improves desirable aroma while reducing beany or safety penalties?
- Which predictions are benchmark-backed, and which are only directional?

The library is most useful as a formulation-screening and prioritization system. It is not a replacement for final experimental confirmation.

---

## Does It Work?

Yes, within a clearly defined envelope. The repository now keeps two trust surfaces separate:

- a compact validation overview for the strongest quantitative benchmark evidence
- a validated-envelope figure that shows where that evidence boundary stops

Current in-repo validation summary:

- 8 supported benchmarks are tracked in Docker-validated artifacts
- 4 benchmarks are strict-ready today, all in the free-precursor envelope
- 9 matched compounds define the current authoritative quantitative proof surface
- the median matched-compound ratio in that proof surface is 1.118x

In practice, this means the software already reproduces a narrow but real set of literature systems closely enough to support quantitative screening inside that envelope.

The main validation figure now focuses only on the two panels that matter most for first-pass trust: parity against literature and benchmark-level quantitative error.

![Validation Overview](results/validation/validation_overview.png)

The benchmark-specific figure below shows what that looks like for an individual literature case: measured versus predicted concentrations on parity axes, the absolute ppb values per compound, and the benchmark summary used to assign status.

![Farmer Benchmark Comparison](results/validation/cys_glucose_150C_Farmer1999_comparison.png)

How to interpret trust:

- Free-precursor benchmarks are the quantitative proof surface. Use them when you need concentration-scale decisions.
- Pea and soy matrix paths are useful for directional prioritization, not yet for release-grade quantitative claims.
- Extrusion-heavy or intact-protein systems remain exploratory unless you add new benchmark evidence.

---

## Trust Levels

| System Type                                        | Trust Level        | What You Can Safely Use It For                                                                        |
| -------------------------------------------------- | ------------------ | ----------------------------------------------------------------------------------------------------- |
| **Free precursors**                          | **High**     | Quantitative ranking, concentration-scale comparison, candidate screening, safety-aware optimization. |
| **Pea / soy matrices**                       | **Moderate** | Directional comparison, off-flavour triage, hypothesis generation, deciding what to test next.        |
| **Intact protein / extrusion-heavy systems** | **Low**      | Exploratory use only; treat outputs as hypotheses until benchmarked.                                  |

Important caveat: matrix trust is lower because accessibility, retention, and pH-dependent headspace effects are not yet benchmark-closed across real plant matrices.

---

## What The Software Can Do

Within its supported envelope, Maillard can already do the following:

- **Rank likely volatile products for a formulation.**
  Given sugars, amino acids, additives, lipids, pH, temperature, water activity, and time, the pipeline predicts which compounds are most likely to dominate.
- **Estimate concentration-scale outputs in ppb inside the validated free-precursor envelope.**
  This is the main quantitative use case when your system resembles the strict-ready literature benchmarks.
- **Score formulations against sensory goals and penalties.**
  You can target tags such as meaty or roasted and simultaneously penalize beany or safety-related outcomes.
- **Generate scientist-facing reports.**
  The CLI can emit Markdown and JSON reports that include decision summaries, confidence metadata, and domain-of-validity warnings.
- **Optimize precursor mixtures.**
  The optimizer searches pH, temperature, water activity, time, and concentration choices to maximize a target sensory direction under penalties.
- **Explain why a prediction looks the way it does.**
  The reporting layer surfaces dominant compounds, penalties, likely bottlenecks, and why confidence is limited.

## What It Cannot Do Reliably Yet

- It cannot serve as a universal zero-shot predictor for complex plant matrices.
- It cannot replace wet-lab confirmation for extrusion-heavy or strongly peptide-bound systems.
- It cannot yet claim broad quantitative accuracy for matrix headspace release across process states.
- It should not be used as if every output were equally validated; some outputs are quantitative, others are directional.

---

## Recommended Workflow For Useful Results

For a new user, the most productive path is:

1. Check the validation surfaces first so you know whether your intended use case is benchmark-backed.
2. Run a forward prediction with your candidate precursors and process conditions.
3. Re-run with `--report` so you get Markdown and JSON outputs you can inspect or share.
4. If the prediction looks promising, use the optimizer to search nearby formulations instead of tuning by hand.
5. Treat matrix-heavy outputs as directional unless your system is visibly close to the validated envelope.

---

## Quick Start

Recommended setup: Docker for reproducibility. Local Python/conda setup is also documented in [Installation.md](Installation.md).

### 1. Start the environment and inspect the trust surface

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validation-figures
```

This gives you the two artifacts that should be read before interpreting predictions:

- [results/validation/benchmark_summary.md](results/validation/benchmark_summary.md) for the benchmark-by-benchmark contract
- [results/validation/validation_overview.md](results/validation/validation_overview.md) for the repository-level trust snapshot
- [results/validation/validated_envelope.md](results/validation/validated_envelope.md) for the current boundary of reliable use

## What To Run For Each Goal

| Goal                                              | Command                                                                             | What you get                                                                            |
| :------------------------------------------------ | :---------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------- |
| Check whether your use case is benchmark-backed   | `./scripts/docker_maillard.sh validation-figures`                                 | Overview PNGs plus Markdown/JSON summaries of trust and envelope boundaries.            |
| Inspect benchmark-by-benchmark status             | `./scripts/docker_maillard.sh summary`                                            | The current benchmark table with status, coverage, ratios, and notes.                   |
| Generate a prediction for a candidate formulation | `python scripts/run_pipeline.py ... --report`                                     | Compound predictions, decision summary, confidence warnings, and a saved report bundle. |
| Optimize a formulation around a target profile    | `python scripts/optimize_formulation.py ... --report`                             | The best trial, predicted tradeoffs, and an exportable report.                          |
| Inspect one literature benchmark directly         | `./scripts/docker_maillard.sh run python scripts/compare_sim_to_lit.py --lit ...` | A benchmark card with parity plot, absolute yields, and benchmark summary.              |

### 2. Run a forward prediction

Use [scripts/run_pipeline.py](scripts/run_pipeline.py) when you already have a candidate formulation and want to know what it is likely to produce.

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
  --output-dir results/first_run
```

This command is useful because it returns more than a list of compounds. It also gives you:

- a prediction summary for desirable compounds and penalties
- confidence and domain-of-validity warnings
- report artifacts you can compare between runs

Useful flags in this pipeline:

- `--list-precursors` to discover supported precursor names
- `--list-tags` to inspect available sensory tags
- `--aw` and `--time-minutes` to change process severity
- `--protein-type` and `--denaturation-state` to express matrix assumptions
- `--report` to persist a scientist-facing Markdown and JSON bundle

### 3. Optimize a formulation instead of guessing by hand

Use [scripts/optimize_formulation.py](scripts/optimize_formulation.py) when you want the software to search for a better operating point.

```bash
python scripts/optimize_formulation.py \
  --sugars ribose,glucose \
  --amino-acids cysteine,leucine \
  --lipids hexanal \
  --target-tag meaty \
  --minimize-tag beany \
  --protein-type pea_iso \
  --risk-aversion 1.5 \
  --n-iterations 50
```

The optimizer is useful when you want to balance several levers at once:

- maximize a sensory direction such as meaty
- penalize off-notes such as beany
- penalize safety risk through the objective
- search concentration ratios and process settings faster than manual trial-and-error

If you want a persistent report of the best trial, add `--report --output-dir ...`.

### 4. Generate benchmark-specific validation cards when you need them

If you need to inspect one literature benchmark directly, regenerate its comparison card:

```bash
./scripts/docker_maillard.sh run \
  python scripts/compare_sim_to_lit.py \
    --lit data/benchmarks/cys_glucose_150C_Farmer1999.json
```

This writes a PNG, Markdown summary, and JSON payload in [results/validation](results/validation).

---

## How To Decide Whether A Result Is Actionable

Use this rule of thumb:

- If your system is close to the free-precursor benchmarks, the output is suitable for quantitative screening and prioritization.
- If your system depends heavily on pea or soy matrix behavior, use the output to decide what to test next, not as a final concentration claim.
- If your system relies on extrusion history, intact proteins, or unsupported matrix physics, treat the result as exploratory.

The repository is designed to tell you not only what it predicts, but also when not to overclaim.

---

## Deeper Documentation

If you need to dig deeper into the science, validation mechanics, or architecture:

- **[docs/README.md](docs/README.md)** - The main index for all deep-dive reference and research documentation.
- **[docs/guides/SCIENTIFIC_RELIABILITY.md](docs/guides/SCIENTIFIC_RELIABILITY.md)** - Detailed breakdown of matrix predictability caps.
- **[docs/VALIDATION_GUIDE.md](docs/VALIDATION_GUIDE.md)** - Our strict validation methodology.
- **[docs/use_cases/README.md](docs/use_cases/README.md)** - Operational reports and the current pea/soy meaty benchmark candidate studies.
