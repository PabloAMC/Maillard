# Matrix Benchmark Intake Checklist

This checklist is the shortest path from a literature report to an honest go/no-go decision on benchmark encoding.

Use it before creating any new benchmark JSON entry for pea or soy matrix systems.

## Required For Benchmark JSON Promotion

Every candidate needs all of the following:

- a clearly specified matrix identity and concentration
- exogenous precursor concentrations or a defensible endogenous-only design
- explicit process conditions: pH, temperature, time, and post-process state notes
- a target panel that includes both desirable meaty sulfur compounds and adverse off-flavour markers when the benchmark claim depends on that tradeoff
- per-compound quantification semantics that are benchmark-usable: absolute concentration or a calibrated ratio with stated internal-standard method
- at least three process replicates or an equivalent uncertainty statement
- a non-detect policy and detection limits where relevant
- citation-level traceability for every encoded numeric value

If any item is missing, the document can still serve as a calibration or intake anchor, but not as a benchmark JSON source.

## Candidate Review Table

| Candidate | Matrix and process fully specified | Meaty targets present | Adverse markers present | Absolute quantification and IS | Replicates / uncertainty | Current verdict |
| --- | --- | --- | --- | --- | --- | --- |
| Pea meaty candidate | Partial | Partial | Partial | No | Partial | Intake design only |
| Soy meaty candidate | Partial | Yes, but qualitative-only proxy | Partial | No | Partial | Intake design only |

## Pea Meaty Candidate

Source report: [pea_matrix_meaty_benchmark.md](pea_matrix_meaty_benchmark.md)

Current usable anchors:

- Mottram and Hofmann support the free-precursor sulfur mechanism.
- Pratap-Singh supports the pea baseline off-flavour context.
- Asen 2022 and Malia 2025 support pea process-state calibration, not benchmark chemistry.

Current blockers:

- no external aqueous PPI study with ribose + cysteine measuring MFT or FFT quantitatively
- no single-run PPI tradeoff panel combining meaty sulfur and adverse lipid markers
- no benchmark-ready absolute MFT or FFT concentration in a pea matrix

Verdict:

- keep as a primary-data design spec
- do not encode as benchmark JSON

## Soy Meaty Candidate

Source report: [soy_matrix_meaty_benchmark.md](soy_matrix_meaty_benchmark.md)

Current usable anchors:

- Nishimura & Abe (2024) now reviewed in full text
- Pratap-Singh 2021 for SPI baseline off-flavour context
- Shu 2024 for adverse-marker reduction under severe soy treatment

What Nishimura actually gives us:

- soy workflow starts from a 75 mg/mL slurry and prepares MRPs from 62.5 mg/mL soy hydrolysate
- reaction condition is 16.5 mM cysteine + 16.5 mM ribose at 95 C for 90 min
- volatile analysis uses HS-SPME-GC/MS with retention-index and library matching
- output is relative peak area behavior with z-transformed clustering across n = 3 replicates

What Nishimura does not give us:

- absolute ppb concentrations for MFT or FFT
- an internal-standard quantitative method for the volatile sulfur panel
- a direct SPI-isolate benchmark condition matching the target unlock design

Verdict:

- promote to qualitative soy-hydrolysate intake anchor only
- do not encode as benchmark JSON
- still require primary SPI data for P3 unlock work

## Immediate P3 Implication

The checklist closes one uncertainty: soy no longer needs more literature triage to know whether the chemistry exists in a protein-matrix workflow.

The remaining blocker is now explicitly primary data:

- PPI and SPI benchmark runs must be generated with absolute compound quantification, internal standards, and a combined meaty/off-flavour panel.
