# PPI/SPI Primary Benchmark Protocol

This document is the first review-ready internal protocol for generating the primary wet-lab benchmark package that still blocks broader matrix trust.

It is not a literature summary. It is the operational contract for the experiment that would let the repository promote pea and soy matrix systems from directional intake toward benchmark-backed target ranking.

## Why This Protocol Exists

The public literature now closes the search question but not the benchmark gap.

- pea isolate still lacks a quantitative meaty-positive matrix benchmark with ribose and cysteine
- soy isolate still lacks a quantitative SPI benchmark that measures meaty sulfur targets and adverse off-flavour markers in the same run
- Nishimura resolves soy chemistry qualitatively but not benchmark quantitation

That means the remaining blocker is primary data, not more reading.

## Experimental Goal

Generate a benchmark-grade aqueous isolate dataset for:

- pea protein isolate plus ribose plus cysteine
- soy protein isolate plus ribose plus cysteine

The dataset must support:

- ranked desirable sulfur targets
- ranked adverse lipid or matrix markers
- process-state calibration
- safety review for the high-temperature soy arm

## Matrix Arms

### Pea Arm

- matrix: pea protein isolate, 5 percent w/v buffered slurry
- pH: 5.5
- temperature: 95 C
- time points: 0, 30, 60, 120, 240 min
- purpose: establish the first quantitative PPI meaty-positive matrix benchmark under moderate aqueous heating

### Soy Arm

- matrix: soy protein isolate, 5 percent w/v buffered slurry
- pH: 5.8
- temperature: 120 C
- time points: 0, 30, 60, 120, 240 min
- purpose: establish the first quantitative SPI benchmark that captures meaty generation and off-flavour suppression under a realistic high-severity condition

## Shared Precursor Addition

- D-ribose: 1.0 mM exogenous
- L-cysteine: 1.0 mM exogenous

These levels are intentionally simple and benchmark-friendly. The goal is not commercial optimization in this protocol; it is a clean cross-benchmark contract.

## Required Analytical Panel

### Desirable Targets

- 2-methyl-3-furanthiol
- 2-furfurylthiol
- bis(2-methyl-3-furyl) disulfide
- 2,5-dimethylpyrazine

### Adverse Targets

- hexanal
- 2-pentylfuran
- furfural
- 1-octen-3-ol

### Safety

- acrylamide in the high-temperature soy arm

## Required Methods

- headspace: HS-SPME-GC-MS
- MFT and FFT: thiol-specific derivatization plus an internal-standard route suitable for absolute quantification
- aldehydes and furans: calibrated quantification with hexanal-d12 or equivalent internal standard
- acrylamide: LC-MS/MS in the soy 120 C arm

## Companion Process-State Assays

The benchmark should not be volatile-only.

Run these on the same lots:

- Ellman free SH
- OPA free amino groups
- DSC or equivalent denaturation proxy
- post-heating pH measurement

These assays are needed to close the gap between matrix observability and process-state calibration.

## Minimum Metadata Contract

Every run must report:

- matrix supplier and lot number
- protein content and solids basis
- exact pH before and after heating
- process vessel or heating format
- number of process replicates
- LOD and LOQ where relevant
- non-detect policy
- internal-standard details

## Promotion Criteria Inside The Repo

This protocol becomes benchmark-eligible only if it yields:

- absolute concentrations for the target panel
- at least three process replicates
- explicit uncertainty per compound or equivalent uncertainty statement
- a single-run meaty versus off-flavour tradeoff view
- traceable source metadata for every encoded numeric value

If those conditions are not met, the data can still support calibration, but not strict benchmark promotion.

## Machine-Readable Companion

The matching machine-readable contract is in [../../data/protocols/ppi_spi_primary_benchmark_contract.json](../../data/protocols/ppi_spi_primary_benchmark_contract.json).

That file should be treated as the future ingestion target for turning the experiment into benchmark and calibration artifacts.