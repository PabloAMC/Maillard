# Maillard Deep-Dive Documentation Index

If you are a scientist new to this repository, **please start with the root [../README.md](../README.md)**.

This directory contains the deeper technical, architectural, and scientific documentation for the framework. It is intended for those checking our math, extending the models, or debugging the core architecture.

---

## Validation & Scientific Limits

- [guides/SCIENTIFIC_RELIABILITY.md](guides/SCIENTIFIC_RELIABILITY.md) - Explains trust boundaries for different matrices in detail.
- [guides/SHARING_RESULTS.md](guides/SHARING_RESULTS.md) - How to generate shareable run, comparison, and campaign artifacts with provenance.
- [VALIDATION_GUIDE.md](VALIDATION_GUIDE.md) - Methodology for benchmarking primary literature against predictions.

## Technical Architecture & Pathways

- [architecture.md](architecture.md) - Core system architecture (Docker, Cantera, MACE/xTB).
- [pathways.md](pathways.md) - Overview of the modeled chemical pathways.
- [reference/SCIENTIFIC_REFERENCE.md](reference/SCIENTIFIC_REFERENCE.md) - Canonical scientific reference with validated articles, numeric anchors, comments, and pathway map.
- [SMIRKS_SYSTEM.md](SMIRKS_SYSTEM.md) - How structural transformations are managed in code.
- [xtb_limitations.md](xtb_limitations.md) - Specific physics constraints and workarounds regarding quantum calculations.
- [reference/PROJECT_STRUCTURE.md](reference/PROJECT_STRUCTURE.md) - Codebase organization list.
- [reference/COMMAND_REFERENCE.md](reference/COMMAND_REFERENCE.md) - Raw CLI options and developer commands.

## Protocols & Shareable Surfaces

- [protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md](protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md) - The review-ready internal protocol for the first primary pea/soy matrix benchmark package.

## Research & Developer Logs

- [research/README.md](research/README.md) - Literature notes and biological syntheses.
- [use_cases/README.md](use_cases/README.md) - Operational reports and candidate external benchmark studies.
- [development/README.md](development/README.md) - Previous implementation attempts, the backlog, and decision logs.
