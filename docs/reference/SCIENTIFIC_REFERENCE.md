# Scientific Reference

This is the canonical scientific reference layer for the repository as of 2026-03-18.

It answers two different needs at once:

- a reader-facing overview of the Maillard pathways that matter for this codebase
- a validated reference table of the articles, numeric anchors, and scientific caveats currently trusted by the repository

## 1. Current State Of The Repo

The repository already had useful scientific material, but it was split across multiple files with different purposes.

The main pieces that existed before this document were:

- [../pathways.md](../pathways.md): broad pathway overview and target-compound framing
- [../Maillard_Plant_based.md](../Maillard_Plant_based.md): long-form review of plant-matrix chemistry and process tradeoffs
- [../slr_benchmark_evaluation.md](../slr_benchmark_evaluation.md): benchmark-oriented literature screening with explicit verdicts
- [../../data/lit/canonical_systems.json](../../data/lit/canonical_systems.json): a small machine-readable seed of canonical precursor systems
- [../../data/lit/benchmark_intake_registry.json](../../data/lit/benchmark_intake_registry.json): structured intake registry for the latest SLR triage

What did not exist was a single canonical file that combined:

- validated articles only
- the numeric values that matter operationally
- comments on what each article does and does not support
- an easy pathway map for new readers

This document fills that gap for human readers. The machine-readable companions are in [../../data/lit/process_state_calibrations.json](../../data/lit/process_state_calibrations.json), [../../data/lit/safety_reference_payloads.json](../../data/lit/safety_reference_payloads.json), and [../../data/lit/benchmark_intake_registry.json](../../data/lit/benchmark_intake_registry.json).

## 2. Easy Pathway Map

The codebase mainly cares about six reaction families.

### A. Core Carbonyl-Amine Entry

Reducing sugars react with amino groups to form a Schiff base and then an Amadori or Heyns product.

Operational meaning for the repo:

- this is the entry point for the flavor cascade
- sugar identity and pH change which downstream branch dominates
- ribose is much more valuable than glucose when the target is sulfur-driven meaty chemistry

### B. Strecker Degradation

Reactive dicarbonyls attack amino acids and generate Strecker aldehydes.

Operational meaning for the repo:

- methionine gives methional
- leucine gives 3-methylbutanal
- isoleucine gives 2-methylbutanal
- these compounds help roasted and savory character, but they are not enough by themselves to make a convincing meaty-positive plant matrix

### C. Sulfur Pathway

This is the decisive branch for meat-like aroma.

Operational meaning for the repo:

- ribose or another pentose plus cysteine is the main benchmark precursor family
- the key targets are 2-methyl-3-furanthiol and 2-furfurylthiol
- Mottram and Hofmann anchor the free-precursor chemistry
- the big current gap is not the free chemistry but the matrix translation in pea and soy

### D. Lipid-Maillard Crosstalk

Plant matrices bring lipid oxidation chemistry whether we want it or not.

Operational meaning for the repo:

- hexanal, 2-pentylfuran, nonanal, and 1-octen-3-ol are the main adverse markers
- these compounds are not just nuisance outputs; they define whether a matrix recommendation is usable
- the missing benchmark family is the one that measures desirable sulfur targets and adverse lipid markers in the same experiment

### E. Matrix Accessibility And Release

Even when the free chemistry is correct, real protein matrices cap what reaches headspace.

Operational meaning for the repo:

- `denaturation_state` controls how open the matrix is
- `cysteine_accessibility` controls how much sulfur chemistry is actually available
- `volatile_retention` controls how much of the generated signal escapes into observable headspace
- these are the main reasons matrix predictions remain directional instead of release-grade quantitative evidence

### F. Safety Branch

The main safety branch currently modeled is acrylamide.

Operational meaning for the repo:

- asparagine plus reducing sugar at high temperature is the key risk family
- Parker and Knol remain the kinetic anchors inside the current model
- Squeo 2023 now gives the repo an industrial endpoint range for plant-protein ingredients, but not a dynamic kinetic benchmark

## 3. Validated Scientific Anchors

The table below only includes references that the repository is currently willing to use as validated anchors, parameter references, or benchmark-intake signals. Every entry now includes a DOI string and a DOI URL so readers can verify the source directly.

| Reference | DOI | Verification URL | Role In Repo | Key Numeric Values | What It Supports | Comment |
| --- | --- | --- | --- | --- | --- | --- |
| Hofmann & Schieberle (1998) | 10.1021/jf9705983 | [https://doi.org/10.1021/jf9705983](https://doi.org/10.1021/jf9705983) | Free-precursor sulfur anchor | MFT and FFT formation in ribose plus cysteine model systems; article-level anchor validated | Confirms that pentose plus cysteine is the correct free-chemistry family for MFT and FFT | Strong chemistry anchor; not a plant-matrix benchmark |
| Mottram & Nobrega (2002) | 10.1021/jf0200826 | [https://doi.org/10.1021/jf0200826](https://doi.org/10.1021/jf0200826) | Free-precursor sulfur anchor | pH 5, 95 C, 4 h; ribose carbon skeleton retained in MFT and FFT; pages 4080-4086 validated | Mechanistic support for ribose plus cysteine sulfur pathway and benchmark candidate design | Strong mechanistic anchor; no direct PPI or SPI matrix data |
| Pratap-Singh et al. (2021) | 10.3390/molecules26134104 | [https://doi.org/10.3390/molecules26134104](https://doi.org/10.3390/molecules26134104) | Matrix headspace anchor | PPI 2-pentylfuran 638 +/- 49 ppb-equivalent; SPI 2-pentylfuran 2492 +/- 199 ppb-equivalent; compound-specific soy vs pea release ratios | Baseline off-flavour anchors for ambient pea and soy slurries | Valid for native matrix headspace, not meaty-positive induction |
| Shu et al. (2024) | 10.1016/j.ultsonch.2023.106675 | [https://doi.org/10.1016/j.ultsonch.2023.106675](https://doi.org/10.1016/j.ultsonch.2023.106675) | Conditional matrix tradeoff calibration | SPI hexanal reduction 70.60%; (E)-2-hexenal reduction 95.60%; 1-octen-3-ol reduction 61.23%; 2-pentylfuran not detected after treatment | Soy off-flavour attenuation under high-severity treatment | Useful for adverse-marker calibration only; no meaty sulfur panel |
| Asen et al. (2022) | 10.3389/fnut.2022.852225 | [https://doi.org/10.3389/fnut.2022.852225](https://doi.org/10.3389/fnut.2022.852225) | Pea process-state calibration | PPC 10% w/v; pH 3/5/7/9; 100 C, 30 min; base Td 74.45 C; heated fractions 124-206 C; triplicates | Best open denaturation-state anchor for pea thermal state versus pH | Parameter anchor, not benchmark chemistry |
| Li et al. (2025) | 10.1016/j.crfs.2025.101173 | [https://doi.org/10.1016/j.crfs.2025.101173](https://doi.org/10.1016/j.crfs.2025.101173) | Pea process-state calibration | Pea protein 3% w/w; Ellman DTNB with extinction coefficient 1.36e4; free SH in nmol/mg protein; triplicates | Best open free-SH accessibility anchor for pea heating response | Parameter anchor; still not the exact benchmark condition |
| Squeo et al. (2023) | 10.3390/foods12061331 | [https://doi.org/10.3390/foods12061331](https://doi.org/10.3390/foods12061331) | Safety reference anchor | Soy wet-extraction acrylamide 185-748 ug/kg; wet-extraction mean 451 ug/kg; 3 replicates; LC-MS/MS with d3-acrylamide; LOD 7 ng/mL; LOQ 24 ng/mL | Industrial endpoint reference for acrylamide in plant-protein ingredients | Honest safety reference, not a dynamic kinetic benchmark |
| Nishimura & Abe (2024) | 10.1016/j.foodchem.2024.141599 | [https://doi.org/10.1016/j.foodchem.2024.141599](https://doi.org/10.1016/j.foodchem.2024.141599) | Qualitative soy chemistry intake | Soy starting slurry 75 mg/mL; MRP mixture 62.5 mg/mL SPH + 16.5 mM cysteine + 16.5 mM ribose; 95 C, 90 min; HS-SPME-GC/MS with n = 3; volatile output reported as relative peak areas / z-transformed clustering | Confirms soy-hydrolysate protein-matrix sulfur chemistry and supports a soy benchmark-intake design | Full text now reviewed: useful as a qualitative intake anchor only, not as an absolute ppb or internal-standard benchmark |

## 4. Numeric Anchors By Module

### Headspace And Matrix Adverse Markers

- Pratap-Singh 2021:
  - pea 2-pentylfuran: 638 +/- 49 ppb-equivalent
  - soy 2-pentylfuran: 2492 +/- 199 ppb-equivalent
- Shu 2024:
  - soy hexanal post-treatment reduction: 70.60%
  - soy (E)-2-hexenal post-treatment reduction: 95.60%
  - soy 1-octen-3-ol post-treatment reduction: 61.23%
  - soy 2-pentylfuran: non-detected after the reported 120 C treatment
- Nishimura 2024:
  - soy starting slurry before hydrolysis: 75 mg/mL
  - soy hydrolysate used for MRP preparation: 62.5 mg/mL
  - cysteine: 16.5 mM
  - ribose: 16.5 mM
  - MRP condition: 95 C, 90 min
  - volatile method: HS-SPME-GC/MS, 200 uL sample, 10 min equilibration at 90 C, 15 min extraction at 90 C
  - output semantics: relative peak areas and z-transformed heatmaps, not absolute ppb concentrations

### Pea Process-State Calibration

- Asen 2022:
  - pea protein concentrate: 10% w/v
  - heating condition: 100 C, 30 min
  - pH sweep: 3, 5, 7, 9
  - base Td: 74.45 C
  - heated fraction Td range: 124-206 C
- Li 2025:
  - pea protein solution: 3% w/w
  - Ellman assay extinction coefficient: 1.36 x 10^4
  - free-SH units: nmol/mg protein
  - replication: triplicates

### Safety

- Squeo 2023:
  - soy wet-extraction acrylamide range: 185-748 ug/kg
  - wet-extraction mean: 451 ug/kg
  - replication: 3
  - LC-MS/MS calibration range: 1-230 ng/mL
  - calibration linearity: R^2 = 0.999
  - LOD: 7 ng/mL
  - LOQ: 24 ng/mL

## 5. Structural Gaps Still Open

These are not vague literature wishes. They are confirmed structural gaps in the current public evidence.

- No quantitative aqueous PPI benchmark with ribose plus cysteine measuring MFT or FFT and hexanal in the same run
- No quantitative aqueous SPI benchmark with ribose plus cysteine measuring MFT or FFT and adverse off-flavour markers in the same run; Nishimura closes the qualitative chemistry question for soy hydrolysate but not the quantitative benchmark gap
- No direct MFT or FFT retention measurement in PPI or SPI versus free-precursor systems
- No benchmark-eligible time series for PPI or SPI sulfur chemistry under the target conditions
- No single-run tradeoff benchmark combining meaty sulfur targets, adverse lipid markers, and safety readouts

## 6. Reader Guidance

If you want to understand the chemistry quickly, read the sections in this order:

1. sulfur pathway
2. lipid-Maillard crosstalk
3. matrix accessibility and release
4. safety branch

If you want the operational state of evidence, read:

1. [../slr_benchmark_evaluation.md](../slr_benchmark_evaluation.md)
2. [../../data/lit/benchmark_intake_registry.json](../../data/lit/benchmark_intake_registry.json)
3. [../../data/lit/process_state_calibrations.json](../../data/lit/process_state_calibrations.json)
4. [../../data/lit/safety_reference_payloads.json](../../data/lit/safety_reference_payloads.json)

If you want the benchmark candidates that translate this science into a missing-data package, read:

1. [../use_cases/pea_matrix_meaty_benchmark.md](../use_cases/pea_matrix_meaty_benchmark.md)
2. [../use_cases/soy_matrix_meaty_benchmark.md](../use_cases/soy_matrix_meaty_benchmark.md)
