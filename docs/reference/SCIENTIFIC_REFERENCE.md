# Scientific Reference

This is the canonical reader-facing scientific reference for the repository as of 2026-03-18.

It has one job: explain which parts of Maillard chemistry matter for this codebase, which papers are trusted for each part, and where the evidence is still structurally missing.

The machine-readable companions remain:

- [../../data/lit/process_state_calibrations.json](../../data/lit/process_state_calibrations.json)
- [../../data/lit/safety_reference_payloads.json](../../data/lit/safety_reference_payloads.json)
- [../../data/lit/benchmark_intake_registry.json](../../data/lit/benchmark_intake_registry.json)

## 1. How To Read This File

If you are new to the project, read the sections in this order:

1. What the repository is actually modeling
2. Reaction map for a meaty-positive system
3. What is well supported vs only approximated
4. Reference tables by scientific role
5. Structural gaps confirmed by the SLR

This file is intentionally not a full literature review. It is a map of the chemistry that is operationally important for the repository.

## 2. What The Repository Is Actually Modeling

The repository is not trying to simulate every elementary step in the full food matrix from first principles.

It is trying to predict a usable formulation signal for alternative-protein systems by combining:

- core free-precursor Maillard chemistry;
- sulfur-driven meaty target formation;
- adverse lipid-oxidation markers;
- matrix-dependent accessibility and release;
- a small safety layer centered on acrylamide.

In practical terms, the current question is not “can Maillard happen?” but “under which precursor, process, and matrix conditions do useful savory targets survive strongly enough to matter more than off-flavours and safety penalties?”

## 3. Reaction Map For A Meaty-Positive System

### 3.1 Entry Chemistry: Carbonyl Plus Amine

Reducing sugars react with amino groups to form a Schiff base and then Amadori or Heyns intermediates.

Why it matters here:

- this is the gateway to the whole flavor cascade;
- sugar identity changes the downstream branch balance;
- ribose is far more useful than glucose when the target is sulfur-driven meaty chemistry.

### 3.2 Dicarbonyl And Strecker Chemistry

Reactive carbonyl intermediates attack amino acids and generate Strecker aldehydes.

Why it matters here:

- methionine supports methional;
- leucine supports 3-methylbutanal;
- isoleucine supports 2-methylbutanal;
- these compounds help roasted and savory character, but they do not by themselves create a convincing meat-like plant system.

### 3.3 Sulfur Branch: The Main Meaty Pathway

This is the decisive branch for the repository's main savory targets.

Why it matters here:

- ribose plus cysteine is the canonical free-precursor family for meat-like sulfur chemistry;
- the main targets are 2-methyl-3-furanthiol (MFT) and 2-furfurylthiol (FFT);
- Hofmann and Mottram anchor the free chemistry;
- the main scientific gap is not whether the pathway exists, but how much of that chemistry survives translation into pea and soy matrices.

### 3.4 Lipid-Maillard Crosstalk: Why Matrices Fight Back

Plant matrices contribute their own oxidative volatile background.

Why it matters here:

- hexanal, 2-pentylfuran, nonanal, and 1-octen-3-ol are the main adverse markers in this repository;
- these compounds are part of the optimization target, not a side note;
- the missing benchmark family is the one that measures desirable sulfur targets and adverse lipid markers in the same experiment.

### 3.5 Matrix Accessibility And Release: Why Free Chemistry Is Not Enough

Even if the precursor chemistry is correct, protein matrices limit what can react and what can escape into headspace.

Why it matters here:

- `denaturation_state` controls how open the matrix is;
- `cysteine_accessibility` controls how much sulfur chemistry is actually available;
- `lysine_accessibility` affects how much amino reactivity remains available to the broader network;
- `volatile_retention` controls how much generated signal reaches observable headspace.

This is the main reason matrix predictions remain directional instead of fully quantitative.

### 3.6 Safety Branch

The main safety branch currently modeled is acrylamide.

Why it matters here:

- asparagine plus reducing sugar at elevated temperature is the core risk family;
- Parker and Knol remain the kinetic anchors inside the current model layer;
- Squeo 2023 provides an industrial endpoint range for plant-protein ingredients, but not a dynamic time-resolved kinetic benchmark.

## 4. What Is Well Supported Vs What Is Approximated

### Well supported in this repository

- free-precursor sulfur chemistry for ribose plus cysteine;
- off-flavour baselines for pea and soy native matrices;
- pea process-state calibration for denaturation and free sulfhydryl accessibility;
- industrial endpoint range for acrylamide in commercial plant-protein ingredients.

### Only approximated today

- translation from free chemistry to real pea and soy headspace intensity;
- direct retention of MFT and FFT inside protein matrices;
- matrix-specific tradeoff between desirable sulfur products and adverse lipid markers in one benchmark;
- process-history realism beyond the currently calibrated public anchors.

### Structurally missing in the public literature

- a quantitative PPI meaty-positive benchmark with ribose plus cysteine and simultaneous adverse-marker readout;
- a quantitative SPI meaty-positive benchmark with the same dual-panel readout;
- direct MFT or FFT retention measurements in pea or soy matrices;
- benchmark-eligible time series for matrix sulfur chemistry under the target conditions.

## 5. References By Scientific Role

The sections below separate papers by what they actually support. This is easier to read than a single mixed table because it prevents benchmark-ready evidence from being confused with calibration-only evidence.

### 5.1 Direct Chemistry Anchors

These papers support the core free-precursor chemistry that the repository treats as mechanistically reliable.

| Reference | DOI | Verification URL | Role In Repo | Key Numeric Values | What It Supports | Comment |
| --- | --- | --- | --- | --- | --- | --- |
| Hofmann & Schieberle (1998) | 10.1021/jf9705983 | [https://doi.org/10.1021/jf9705983](https://doi.org/10.1021/jf9705983) | Free-precursor sulfur anchor | MFT and FFT formation in ribose plus cysteine model systems; article-level anchor validated | Confirms that pentose plus cysteine is the correct free-chemistry family for MFT and FFT | Strong chemistry anchor; not a plant-matrix benchmark |
| Mottram & Nobrega (2002) | 10.1021/jf0200826 | [https://doi.org/10.1021/jf0200826](https://doi.org/10.1021/jf0200826) | Free-precursor sulfur anchor | pH 5, 95 C, 4 h; ribose carbon skeleton retained in MFT and FFT; pages 4080-4086 validated | Mechanistic support for ribose plus cysteine sulfur pathway and benchmark candidate design | Strong mechanistic anchor; no direct PPI or SPI matrix data |

### 5.2 Matrix And Headspace Calibration Anchors

These papers do not close the meaty-positive benchmark gap, but they do provide numerical anchors for how pea and soy systems behave before or after processing.

| Reference | DOI | Verification URL | Role In Repo | Key Numeric Values | What It Supports | Comment |
| --- | --- | --- | --- | --- | --- | --- |
| Pratap-Singh et al. (2021) | 10.3390/molecules26134104 | [https://doi.org/10.3390/molecules26134104](https://doi.org/10.3390/molecules26134104) | Matrix headspace anchor | PPI 2-pentylfuran 638 +/- 49 ppb-equivalent; SPI 2-pentylfuran 2492 +/- 199 ppb-equivalent; compound-specific soy vs pea release ratios | Baseline off-flavour anchors for ambient pea and soy slurries | Valid for native matrix headspace, not meaty-positive induction |
| Shu et al. (2024) | 10.1016/j.ultsonch.2023.106675 | [https://doi.org/10.1016/j.ultsonch.2023.106675](https://doi.org/10.1016/j.ultsonch.2023.106675) | Conditional matrix tradeoff calibration | SPI hexanal reduction 70.60%; (E)-2-hexenal reduction 95.60%; 1-octen-3-ol reduction 61.23%; 2-pentylfuran not detected after treatment | Soy off-flavour attenuation under high-severity treatment | Useful for adverse-marker calibration only; no meaty sulfur panel |
| Asen et al. (2022) | 10.3389/fnut.2022.852225 | [https://doi.org/10.3389/fnut.2022.852225](https://doi.org/10.3389/fnut.2022.852225) | Pea process-state calibration | PPC 10% w/v; pH 3/5/7/9; 100 C, 30 min; base Td 74.45 C; heated fractions 124-206 C; triplicates | Best open denaturation-state anchor for pea thermal state versus pH | Parameter anchor, not benchmark chemistry |
| Li et al. (2025) | 10.1016/j.crfs.2025.101173 | [https://doi.org/10.1016/j.crfs.2025.101173](https://doi.org/10.1016/j.crfs.2025.101173) | Pea process-state calibration | Pea protein 3% w/w; Ellman DTNB with extinction coefficient 1.36e4; free SH in nmol/mg protein; triplicates | Best open free-SH accessibility anchor for pea heating response | Parameter anchor; still not the exact benchmark condition |

### 5.3 Safety Anchor

This paper is the current reader-verifiable endpoint anchor for acrylamide in plant-protein ingredients.

| Reference | DOI | Verification URL | Role In Repo | Key Numeric Values | What It Supports | Comment |
| --- | --- | --- | --- | --- | --- | --- |
| Squeo et al. (2023) | 10.3390/foods12061331 | [https://doi.org/10.3390/foods12061331](https://doi.org/10.3390/foods12061331) | Safety reference anchor | Soy wet-extraction acrylamide 185-748 ug/kg; wet-extraction mean 451 ug/kg; 3 replicates; LC-MS/MS with d3-acrylamide; LOD 7 ng/mL; LOQ 24 ng/mL | Industrial endpoint reference for acrylamide in plant-protein ingredients | Honest safety reference, not a dynamic kinetic benchmark |

### 5.4 Qualitative Intake Anchors

These papers are useful because they show that the chemistry exists in a matrix-like system, but they are not benchmark-ready quantitative support.

| Reference | DOI | Verification URL | Role In Repo | Key Numeric Values | What It Supports | Comment |
| --- | --- | --- | --- | --- | --- | --- |
| Nishimura & Abe (2024) | 10.1016/j.foodchem.2024.141599 | [https://doi.org/10.1016/j.foodchem.2024.141599](https://doi.org/10.1016/j.foodchem.2024.141599) | Qualitative soy chemistry intake | Soy starting slurry 75 mg/mL; MRP mixture 62.5 mg/mL SPH + 16.5 mM cysteine + 16.5 mM ribose; 95 C, 90 min; HS-SPME-GC/MS with n = 3; volatile output reported as relative peak areas / z-transformed clustering | Confirms soy-hydrolysate protein-matrix sulfur chemistry and supports a soy benchmark-intake design | Full text reviewed: useful as a qualitative intake anchor only, not as an absolute ppb or internal-standard benchmark |

## 6. Numeric Anchors By Modeling Layer

### 6.1 Sulfur Chemistry

- Hofmann & Schieberle 1998:
  - validates MFT and FFT formation in ribose plus cysteine free systems.
- Mottram & Nobrega 2002:
  - pH 5;
  - 95 C;
  - 4 h;
  - ribose carbon skeleton retained in MFT and FFT.

### 6.2 Headspace And Adverse Markers

- Pratap-Singh 2021:
  - pea 2-pentylfuran: 638 +/- 49 ppb-equivalent;
  - soy 2-pentylfuran: 2492 +/- 199 ppb-equivalent.
- Shu 2024:
  - soy hexanal reduction after treatment: 70.60%;
  - soy (E)-2-hexenal reduction after treatment: 95.60%;
  - soy 1-octen-3-ol reduction after treatment: 61.23%;
  - soy 2-pentylfuran not detected after the reported 120 C treatment.
- Nishimura 2024:
  - soy starting slurry before hydrolysis: 75 mg/mL;
  - soy hydrolysate used for MRP preparation: 62.5 mg/mL;
  - cysteine: 16.5 mM;
  - ribose: 16.5 mM;
  - MRP condition: 95 C, 90 min;
  - volatile method: HS-SPME-GC/MS, 200 uL sample, 10 min equilibration at 90 C, 15 min extraction at 90 C;
  - output semantics: relative peak areas and z-transformed heatmaps, not absolute ppb concentrations.

### 6.3 Process-State Calibration

- Asen 2022:
  - pea protein concentrate: 10% w/v;
  - heating condition: 100 C, 30 min;
  - pH sweep: 3, 5, 7, 9;
  - base Td: 74.45 C;
  - heated fraction Td range: 124-206 C.
- Li 2025:
  - pea protein solution: 3% w/w;
  - Ellman assay extinction coefficient: 1.36 x 10^4;
  - free-SH units: nmol/mg protein;
  - replication: triplicates.

### 6.4 Safety

- Squeo 2023:
  - soy wet-extraction acrylamide range: 185-748 ug/kg;
  - wet-extraction mean: 451 ug/kg;
  - replication: 3;
  - LC-MS/MS calibration range: 1-230 ng/mL;
  - calibration linearity: R^2 = 0.999;
  - LOD: 7 ng/mL;
  - LOQ: 24 ng/mL.

## 7. Structural Gaps Confirmed By The SLR

These are confirmed public-evidence gaps, not just missing search effort.

- No quantitative aqueous PPI benchmark with ribose plus cysteine measuring MFT or FFT and hexanal in the same run.
- No quantitative aqueous SPI benchmark with ribose plus cysteine measuring MFT or FFT and adverse off-flavour markers in the same run.
- Nishimura closes the qualitative soy-hydrolysate chemistry question, but not the quantitative SPI benchmark gap.
- No direct MFT or FFT retention measurement in PPI or SPI versus free-precursor systems.
- No benchmark-eligible time series for PPI or SPI sulfur chemistry under the target conditions.
- No single-run tradeoff benchmark combining meaty sulfur targets, adverse lipid markers, and safety readouts.

## 8. What A Reader Should Conclude

The repository has a credible chemistry core and a meaningful matrix-calibration layer, but it does not yet have a public quantitative benchmark that closes the final translation from free sulfur chemistry to real pea or soy headspace.

That means the honest interpretation is:

- free-precursor sulfur direction is scientifically strong;
- matrix penalties and release effects are informed, not fully benchmark-closed;
- safety has an endpoint anchor but not a complete dynamic calibration;
- the missing data package is now clearly defined enough to guide either internal experiments or selective computational support.

## 9. Next Documents To Read

If you want the benchmark triage and evidence verdicts, read:

1. [../slr_benchmark_evaluation.md](../slr_benchmark_evaluation.md)
2. [../../data/lit/benchmark_intake_registry.json](../../data/lit/benchmark_intake_registry.json)

If you want the machine-readable calibration payloads, read:

1. [../../data/lit/process_state_calibrations.json](../../data/lit/process_state_calibrations.json)
2. [../../data/lit/safety_reference_payloads.json](../../data/lit/safety_reference_payloads.json)

If you want the missing-data packages that would close the main benchmark gaps, read:

1. [../use_cases/pea_matrix_meaty_benchmark.md](../use_cases/pea_matrix_meaty_benchmark.md)
2. [../use_cases/soy_matrix_meaty_benchmark.md](../use_cases/soy_matrix_meaty_benchmark.md)
