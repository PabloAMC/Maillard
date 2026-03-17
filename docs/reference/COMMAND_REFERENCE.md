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
| `./scripts/docker_maillard.sh coverage-gaps` | Coverage-gap report for P3-style expansion work |

## CLI Workflows Outside Docker

| Command | Purpose |
| --- | --- |
| `python scripts/run_pipeline.py --list-precursors` | List available precursors |
| `python scripts/run_pipeline.py --list-tags` | List sensory/optimization tags |
| `python scripts/run_pipeline.py ...` | Run one forward prediction |
| `python scripts/optimize_formulation.py ...` | Search a formulation space |
| `python scripts/compare_sim_to_lit.py` | Compare the framework against literature benchmarks |

## Recommended Review Sequence

If you are reviewing the state of the repository rather than developing it, the shortest useful sequence is:

```bash
./scripts/docker_maillard.sh summary
./scripts/docker_maillard.sh validated-envelope
./scripts/docker_maillard.sh validation-figures
./scripts/docker_maillard.sh matrix-readiness
./scripts/docker_maillard.sh coverage-gaps
```
