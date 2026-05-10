# Command Reference

This reference focuses on the commands a scientist or reviewer is most likely to need.

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
| `./scripts/docker_maillard.sh campaign data/campaigns/shareable_meaty_screen.yml` | Generate a review-ready campaign package with run-level and campaign-level artifacts |

## CLI Workflows Outside Docker

| Command | Purpose |
| --- | --- |
| `python scripts/run_pipeline.py --list-precursors` | List available precursors |
| `python scripts/run_pipeline.py --list-tags` | List sensory/optimization tags |
| `python scripts/run_pipeline.py ...` | Run one forward prediction |
| `python scripts/optimize_formulation.py ...` | Search a formulation space |
| `python scripts/run_campaign.py --names "A,B" --ph 5.5 --temp 105` | Generate a side-by-side comparison package via the campaign pipeline |
| `python scripts/run_campaign.py --spec data/campaigns/shareable_meaty_screen.yml` | Generate a shareable multi-run campaign package |
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
