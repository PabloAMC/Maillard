# Project Structure

This is the shortest useful map of the repository.

## Top-Level Layout

| Path | Purpose |
| --- | --- |
| `src/` | Core chemistry, kinetics, projection, validation, and reporting code |
| `scripts/` | Reproducible entry points for workflows and artifact generation |
| `data/` | Benchmarks, species definitions, reactions, formulations, and process profiles |
| `tests/` | Unit, integration, scientific, and QM validation |
| `docs/` | User-facing guides, technical references, research notes, and development context |
| `results/validation/` | Generated validation artifacts used for review |
| `tasks/` | Backlog and lessons learned for ongoing repository work |

## Key Source Modules

| Path | Role |
| --- | --- |
| `src/recommend.py` | Main FAST observable prediction path |
| `src/benchmark_validation.py` | Benchmark execution, summaries, matrix evidence, assertions, and deltas |
| `src/validation_contract.py` | Single source of truth for benchmark policy |
| `src/headspace.py` | Observable release and matrix/headspace corrections |
| `src/matrix_correction.py` | Matrix accessibility and retention behavior |
| `src/reporting.py` | Human-readable and JSON report generation |
| `src/usability_reports.py` | Scientist-facing explainability payloads |

## Key Data Areas

| Path | Role |
| --- | --- |
| `data/benchmarks/` | Literature and internal benchmark definitions |
| `data/campaigns/` | Shareable screening campaign specifications |
| `data/protocols/` | Machine-readable primary benchmark protocol contracts |
| `data/species/` | Target sets for desirable, adverse, and toxic markers |
| `data/reactions/` | Curated pathways and reaction families |
| `data/formulation_grid.yml` | Named formulation definitions for explainability and comparisons |
| `data/temp_profiles/` | Temperature ramp profiles |

## Shareable Workflow Entrypoints

| Path | Role |
| --- | --- |
| `scripts/run_pipeline.py` | Single-run scientist-facing prediction and report generation |
| `scripts/run_campaign.py` | Campaign-level package generation with run artifacts, named comparisons, and provenance |

## Test Layout

| Path | Role |
| --- | --- |
| `tests/unit/` | Small, focused behavioral checks |
| `tests/integration/` | Cross-module behavior and workflow checks |
| `tests/scientific/` | Benchmark and scientific regression checks |
| `tests/qm/` | QM-heavy and backend-specific checks |

## Generated Artifacts To Inspect First

| Path | Role |
| --- | --- |
| `results/validation/benchmark_summary.md` | Exact benchmark state by case |
| `results/validation/validated_envelope.md` | Plain-language scientific boundary |
| `results/validation/validation_overview.png` | Visual summary of reliability and limits |
| `results/validation/matrix_benchmark_assertions.md` | Matrix ranking checks |
| `results/validation/benchmark_coverage_gaps.md` | Expansion and evidence gaps |

## Important Distinction

Not every document in `docs/` is meant for a first-time reader.

- use `docs/guides/` for onboarding
- use `docs/reference/` for commands and structure
- use `docs/protocols/` for review-ready experimental contracts
- use `docs/research/` for literature notes and scientific planning
- use `docs/development/` for logs, lessons, and backlog context
