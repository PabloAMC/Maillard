# Vega & Brewer 1994 — EXTRACTION (six lipid-oxidation aldehydes dosed into a 3 % w/v gelatin gel, lipid-free, sniffed at 4 / 22 / 37 / 60 C by a 16-member panel, 3 replicates; 24 detectable odour thresholds in ppb from a linear 75 %-detection fit)

### THE ONLY MATRIX THRESHOLD LADDER IN THE CORPUS WITH A TEMPERATURE AXIS: six compounds x four temperatures = 24 measured thresholds in one matrix, one panel, one method, and the repository has been carrying all 24 at second hand through `k2_matrix_and_thresholds.md` with no dossier behind them until now.

**Source on disk:** `data/articles/vega1994.pdf` (17 pp., Journal of Food Lipids 1 (1994) 229-245).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/vega1994.txt`). The text layer is an OCR of a scanned typescript and is
**dirty in a specific, recoverable way**: the "±" glyph is rendered variously as `f`, `2`, `i`, `+`,
`.` or `_+`, and several compound names are mangled ("Heranal" for Hexanal, "Heptenal" for Heptanal,
"t-2-0ctena1" for t-2-octenal). **The threshold column itself (the DOT column of Table 1) is clean
and unambiguous** — every one of the 24 values is a plain two- or three-digit integer, and all 24
reproduce the values `k2_matrix_and_thresholds.md` sec. A.2 already carries. The percent-correct
cells are legible but their standard deviations must be read through the glyph substitution; each is
flagged where the reading is not certain. Figures 1 (sniffing apparatus), 2 (the regression
construction), 3 (DOT and viscosity vs temperature), 4 (DOT by compound and temperature), 5 (a
chromatogram) and 6 (detector response vs concentration) are images: **the viscosity values
themselves are figure-only** and are not typed as numbers here. There is no supplementary material.
Repo status before this dossier: Vega 1994 is cited by `src/kinetic_core/matrix_oav.py` (the
`_VEGA_GELATIN` ladder, the paraffin-oil rows and six of the water rows), by
`src/kinetic_core/parameters_matrix.py` (the `gelatin_3pct` `MatrixLoading` and the
`unsat_penalty_gelatin` observation) and by `k2_matrix_and_thresholds.md` sec. A.2 / A.3 / A.5 —
but has **no extraction dossier**; every number reaches the code through that synthesis file's
summary table.

## 0. Identity

| field | value |
|---|---|
| Title | "Detectable Odor Thresholds of Selected Lipid Oxidation Compounds at Various Temperatures in a Gelatin Model System" |
| Authors | J. D'Dios Vega and M. Susan Brewer (corresponding) — Division of Foods and Nutrition, University of Illinois, Urbana, IL 61801 |
| Venue | Journal of Food Lipids **1** (1994) 229-245. Received for publication 23 December 1993; accepted 20 March 1994. Copyright 1994 Food & Nutrition Press, Trumbull, Connecticut |
| DOI / article ID | none printed on the scan |
| Matrix | **3 % w/v gelatin dispersion in distilled water** (3 g Knox gelatin per 100 mL), microwave-melted, cooled to 22 C, dosed, then held 18 h at 4 C. **Lipid-free by design** — the authors chose gelatin precisely because it "is devoid of lipid components" |
| Quantity measured | **DOT = "detectable odor threshold"**, defined as the concentration at which **75 % of panellists detect a difference from a gelatin control**, read off a **linear least-squares regression** of percent-correct against concentration (Fig. 2). Not a BET, not a forced choice, **not corrected for chance** |
| Naming | "DOT" throughout; "CD" = coefficient of determination of that regression; "OD" = odor descriptor; "BP" = boiling point; "VP" = vapour pressure at 20 C |
| Companion paper | **Brewer & Vega 1995**, the same two authors, the same six compounds, in a cooked-beef model — `brewer1995_extraction.md` (this batch). The gelatin paper came first and the beef paper cites it |
| Second-hand content | Guadagni et al. 1963 and 1972 aqueous thresholds (Table 1 "Ref. ppb" column); Guadagni et al. 1972 paraffin-oil thresholds (p. 234); Wick et al. 1967 meat-slurry thresholds (p. 234). **None of these was measured here.** |

## 1. Why it matters

`src/kinetic_core/matrix_oav.py` is built on the rule "**THRESHOLDS ARE INPUTS, NEVER PREDICTIONS**"
(module docstring; declaration D.6 Module 7) and on the rule that **nothing is borrowed across
matrices**: a `(compound, matrix)` pair with no measurement returns a `NoMeasuredThreshold` object
carrying the reason, and that object is deliberately not a number so a caller dividing by it gets a
`TypeError` rather than a silent zero. The consequence is the live gap: `MATRIX_THRESHOLDS`
currently contains entries for exactly **two** matrices, and neither is a protein pot the engine can
cook.

**This paper is the source of one of those two.** `_VEGA_GELATIN` in `matrix_oav.py` lines 220-227
is this paper's Table 1 DOT column, and it is the table this dossier feeds:

| repository object | file | what this paper supplies |
|---|---|---|
| `_VEGA_GELATIN` -> `MATRIX_THRESHOLDS` (24 `ThresholdRecord`s, matrix `gelatin_3pct`) | `src/kinetic_core/matrix_oav.py` | **the whole ladder**: 6 compounds x 4 temperatures, Table 1 DOT column |
| the paraffin-oil rows of `MATRIX_THRESHOLDS` (3 records, matrix `paraffin_oil`) | same file | quoted second-hand from Guadagni 1972 on p. 234; **not measured here** |
| six `WATER_THRESHOLDS` records flagged `quoted_second_hand_by_vega1994` | same file | the Table 1 "Ref. ppb" column; **not measured here** |
| `MATRIX_LOADING["gelatin_3pct"]` (30 g protein/L, 0 lipid, `thermal_step_after_dosing = False`) | `src/kinetic_core/parameters_matrix.py` line 249 | the Methods paragraph: 3 g gelatin / 100 mL, no lipid, dosed at 22 C then held at 4 C |
| `ALPHA_BETA_UNSATURATION_OBSERVATIONS["unsat_penalty_gelatin"]` = 2.81x | `src/kinetic_core/parameters_matrix.py` line 495 | the t-2-hexenal / hexanal contrast at 22 C |

**Are these matrix thresholds the repository could carry for `gelatin_3pct`? Yes, and it already
does — but read what kind of number they are before leaning on them.** They are real, measured,
same-matrix, same-panel, same-method thresholds in a defined protein gel at four stated
temperatures. They are the cleanest matrix threshold set in the corpus for three reasons the paper
itself establishes: the matrix is lipid-free (nothing else in it can oxidise and contribute odour),
there is **no thermal step after dosing** (the compound present at perception is the compound
weighed in), and one panel of 16 measured every cell in three replicates, so within-table
comparisons are free of the cross-study, cross-method noise that wrecks every other matrix
comparison in the corpus. What they are **not** is comparable to an aqueous threshold: the criterion
is 75 % detection **uncorrected for chance** against a single control, which is a systematically
higher number than the 50 %-forced-choice values in the "Ref. ppb" column beside them. The sign of
that bias is known; its size is not measured anywhere in this paper.

**What it does NOT fix.** The refusal the odour-activity layer emits on protein pots is not lifted
by this paper. `gelatin_3pct` is a **collagen hydrolysate gel**, not one of the matrices
`data/species/protein_matrices.yml` charges (`blg`, `soy_isolate`, `pea_isolate`), and it has no site
densities on file — gelatin's amine and thiol counts are nowhere in this corpus, and gelatin is
famously cysteine-free, so the binding layer in `src/kinetic_core/matrix_sites.py` would charge it
nothing even if a loading were stated. Transferring a gelatin threshold to a soy or pea pot is
exactly the operation the house rules forbid and the module's `SEALED_OR_REFUSED_MATRICES` machinery
exists to prevent. **This paper turns a second-hand table into a first-hand one; it does not add a
matrix.**

## 2. Methods as they matter to a model

- **The gel.** Gelatin (Knox Inc., Englewood Cliffs, NJ), **3 g into 100 mL distilled water**.
  Microwaved 2 min at 100 % power (700 W), stirred, heated 1 more min, stirred, **cooled to 22 C**.
  So the protein loading is **30 g/L** and there is **no lipid, no sugar, no salt and no buffer**;
  **pH is never stated or measured** anywhere in the paper.
- **Dosing, and the fact that nothing is heated afterwards.** Aldehyde stock solutions were made up
  in distilled water at 22 C **at parts-per-million**, nitrogen-flushed, capped, stored at 4 C.
  **0.1 mL of stock was added to the cooled (22 C) gelatin** to give parts-per-billion, mixed 1 min,
  dispensed as 100 mL into 250 mL sniffing flasks, stoppered with ground glass and **stored at 4 C
  for 18 h**; samples were prepared 24 h before evaluation. **The dose is added AFTER the only heat
  step and is never cooked.** This is the single most important methodological fact in the paper and
  it is why `MATRIX_LOADING["gelatin_3pct"]` carries `thermal_step_after_dosing = False`. It is also
  the axis on which the companion beef paper differs (see `brewer1995_extraction.md`).
- **Concentration is nominal, never verified.** The concentrations in Table 1 are what was weighed
  in. **No headspace or liquid concentration was measured at the moment of sniffing.** The GC work
  (below) used separate 10 mL samples in a different vessel at a different temperature and reports
  only detector area counts. Losses to the headspace of a 250 mL flask holding 100 mL of gel, to the
  glass, and over 18 h at 4 C are not quantified. `concentration_verified = False` in the code is
  correct.
- **Panel.** **16 members**, trained in two 1-h sessions. Aged 21-32 (mean 23), **2 men and 14
  women, all non-smokers**. Trained on suprathreshold concentrations chosen "based on concentrations
  reported in meat systems": pentanal 50 ppm, hexanal 220 ppm, t-2-hexenal 120 ppm, heptanal
  100 ppm, t-2-octenal 150 ppm, t,t-2,4-decadienal 110 ppm.
- **Presentation.** A test flask was presented **alongside a matched flask of plain 3 % gelatin as
  the control** — one sample against one control, i.e. **not** a triangle test and **not** a 3-AFC.
  Each compound was presented **in ascending order of concentration**, to each panellist,
  **3 times**, at each of 4, 22, 37 and 60 C. The panellist squeezed a silicone bulb to push a
  stream of air over the sample surface and out through a sniff tube (Fig. 1: a 250 mL flask,
  1.2 x 12 cm silicone tubing, ground-glass stopper), sniffing as often as wanted with 10-20 s
  between sniffs, then scored intensity **0 = smells like control, 1 = different from control but
  cannot identify, 2 = point of recognition / very weak, 3 = weak, 4 = pronounced, 5 = strong,
  6 = very strong**.
- **Temperature control.** 37 C and 60 C flasks stood in a **37 ± 1 C or 60 ± 1 C water bath before
  and during** the sniff; 4 ± 1 C flasks stood in a tub of ice-chilled water. Flasks equilibrated
  **15 min at the selected temperature between panellists**. **The room was always 22 C, 50 % RH,
  positive air pressure** — so the sample temperature varies and the ambient does not.
- **How the threshold is computed, exactly.** "the percent of correct responses (judges detecting a
  difference compared to the control) was plotted versus concentration at a given temperature
  (Fig. 2). The **linear least squares regression** method was used to determine the best fit curve
  for data points; the concentration at which **75 % of the panellists detected a difference**
  compared to the control was considered to be the detectable odor threshold." Three consequences a
  modeller must carry: (i) it is an **interpolation on a straight line through five points**, not a
  dilution-series step, so the reported value can fall between the tested concentrations and does;
  (ii) **no correction for chance** is applied, although a same-different judgement against one
  control has a substantial guess rate; (iii) the **quality of that line is printed per row as the
  "CD" (coefficient of determination) and it is often poor** — see Flags 2.
- **GC, and what it is not.** Separate 10 mL samples of the same dosed gel went into 22 mL serum
  vials, Teflon-lined caps, **platen-heated at 60 C for 3 min**, 1 mL headspace withdrawn
  automatically. HP 5890 GC, Tekmar 7000/7050 equilibrium headspace autosampler; DB-5 fused silica
  60 m x 0.32 mm ID x 1.0 um film, 9 m x 0.32 mm deactivated guard as the interface; oven 35 C
  (5 min) then a ramp to 280 C (the ramp rate is printed as "SC/min" — an OCR failure, see Flags 6);
  platen 60 C, equilibrium 20 min, vial pressure 10 psi, valve and transfer line 120 C, injector
  210 C, FID 300 C; air 240, hydrogen 30, nitrogen 30 mL/min, helium carrier at 148 kPa. **The
  output is retention time and FID area counts only — there is no calibration to a headspace
  concentration anywhere in this paper**, so nothing here can be converted into a partition
  coefficient.
- **Viscosity.** Brookfield DV-I digital viscometer, **spindle 1, at 10 and 20 rpm**, at 4, 22, 37
  and 60 C, 10 mL samples. **The numbers appear only in Fig. 3** and are therefore figure-only. The
  prose says the viscosity "dropped dramatically between 4C and 22C then remained at a fairly
  constant low level to 60C".
- **Statistics.** Three replications of all analyses. Correlations by Statview 512+ v1.2; regressions
  plotted in Cricket Graph 1.2. **No ANOVA, no confidence interval on any threshold, and no
  significance test on any difference between temperatures or compounds is reported.**

## 3. Tables re-typed

### Table 1. "Sensory responses (percent correct), detectable odor threshold and descriptors for volatiles in a gelatin model system at selected concentrations"

Column headings as printed: `Percent Correct Responses` (footnote 1: "Percent of correct responses
± standard deviation") at five concentrations **A, B, C, D, E** (footnote 2: "A, B, C, D, and E are
selected organic concentrations (in parts per billion, ppb) of each compound" — the five values are
printed as a sub-header under each compound); then `CD` (footnote 3: "CD = coefficient of
determination"); then `DOT ppb`; then `Ref. ppb`; then `OD` (footnote 4: "OD = odor descriptor").
Footnote 5: "BP = boiling point". Footnote 6: "VP = vapor pressure, at 20 C". Footnote v: "Guadagni
et al. 1963." Footnote w: "Guadagni et al. 1972."

**All percent-correct cells are given as `value ± SD`. The "±" is an OCR reconstruction in every
cell** (the scan renders it as `f`, `2`, `i`, `+`, `.` or `_+`); the values either side of it are
legible. Cells where even that is uncertain are marked `[?]`.

**Pentanal** (printed `BP = 103 C, VP ~ 1.5 mm`) — concentrations A-E = **5, 10, 25, 35, 60 ppb**

| Temp, C | A (5) | B (10) | C (25) | D (35) | E (60) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 25.5 ± 15 | 25.5 ± 15 | 68.7 ± 15 | 68.7 ± 15 | 79.9 ± 10 | 0.80 | **47** | | pungent |
| 22 | 30.3 ± 12 | 42.5 ± 17 | 58.5 ± 13 | 92.5 ± 7 | 82.5 ± 15 | 0.72 | **41** | 12 (v) | oily |
| 37 | 50.0 ± 10 | 50.0 ± 10 | 75.5 ± 15 | 92.6 ± 8 | 82.5 ± 13 | 0.64 | **34** | | rancid |
| 60 | 50.0 ± 10 | 72.5 ± 15 | 86.7 ± 13 | 90.0 ± 6 | 92.5 ± 6 | 0.67 | **22** | | green |

**Hexanal** (printed "Heranal"; `BP = 131 C, VP ~ 1 mm`) — A-E = **20, 40, 70, 100, 150 ppb**

| Temp, C | A (20) | B (40) | C (70) | D (100) | E (150) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 30.1 ± 5 | 30.1 ± 8 | 81.3 ± 15 | 95.2 ± 3 | 97.3 ± 2 | 0.79 | **90** | | resinous |
| 22 | 45.2 ± 10 | 90.2 ± 7 | 82.4 ± 25 | 85.2 ± 10 | 97.3 ± 2 | 0.51 | **58** | 4.5 (w) | fresh green |
| 37 | 50.0 ± 12 | 92.3 ± 6 | 95.3 ± 4 | 99.1 ± 5 | 99.1 ± 5 | 0.49 | **34** | | fatty |
| 60 | 65.3 ± 25 | 75.5 ± 10 | 90.9 ± 5 | 90.9 ± 8 | 99.9 ± 1 | 0.87 | **38** | | rancid |

**t-2-Hexenal** (printed `BP = [?] 47 [?], VP ~ 10 mm` — the boiling point is garbled, see Flags 6)
— A-E = **5, 10, 30, 50, 100 ppb**

| Temp, C | A (5) | B (10) | C (30) | D (50) | E (100) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 9.2 ± 8 | 17.3 ± 6 | 26.6 ± 14 | 45.0 ± 12 | 43.8 ± 8 | 0.75 | **170** | | green |
| 22 | 13.9 ± 12 | 29.9 ± 9 | 43.4 ± 20 | 52.7 ± 28 | 66.0 ± 5 | 0.86 | **109** | 3 (v) | leafy |
| 37 | 24.9 ± 18 | 44.6 ± 12 | 56.4 ± 24 | 69.1 ± 26 | 79.9 ± 17 | 0.82 | **79** | | fragrant |
| 60 | 31.2 ± 18 | 53.3 ± 3 | 68.1 ± 8 | 83.3 ± 12 | 86.6 ± 6 | 0.72 | **60** | | sweet |

**Heptanal** (printed "Heptenal"; `BP = 153 C, VP = 0.5 mm`) — A-E = **5, 10, 20, 50, 100 ppb**

| Temp, C | A (5) | B (10) | C (20) | D (50) | E (100) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 12.3 ± 12 | 32.5 ± 15 | 32.5 ± 15 | 50.0 ± 13 | 68.4 ± 19 | 0.89 | **108** | | woody |
| 22 | 15.4 ± 15 | 43.2 ± 13 | 62.3 ± 15 | 75.4 ± 15 | 75.4 ± 15 | 0.57 | **79** | 3 (w) | nutty |
| 37 | 32.5 ± 15 | 52.5 ± 18 | 67.4 ± 12 | 72.6 ± 19 | 89.5 ± 10 | 0.77 | **62** | | grass |
| 60 | 45.5 ± 12 | 45.5 ± 12 | 65.8 ± 13 | 92.6 ± 7 | 92.6 ± 7 | 0.76 | **50** | | sharp sweet |

**t-2-Octenal** (printed `BP = 85 C, VP = 2 mm` — see Flags 6) — A-E = **10, 30, 50, 100, 150 ppb**

| Temp, C | A (10) | B (30) | C (50) | D (100) | E (150) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 10.2 ± 4 | 28.2 ± 4 | 35.9 ± 4 | 64.1 ± 4 | 74.4 ± 8 | 0.96 | **140** | | floral |
| 22 | 17.9 ± 8 | 38.4 ± 8 | 61.5 ± 8 | 79.4 ± 12 | 84.6 ± 15 | 0.85 | **109** | 3 (w) | fragrant |
| 37 | 30.7 ± 13 | 48.7 ± 9 | 64.1 ± 16 | 64.7 ± 19 | 94.9 ± 2 | 0.89 | **105** | | grass |
| 60 | 33.3 ± 12 | 56.4 ± 9 | 74.3 ± 4 | 87.2 ± 9 | 97.4 ± 2 | 0.86 | **81** | | herbal |

**t,t-2,4-Decadienal** (printed `BP = 115 C, VP = 8 mm` — see Flags 6) — A-E = **5, 10, 50, 75,
100 ppb**

| Temp, C | A (5) | B (10) | C (50) | D (75) | E (100) | CD | **DOT ppb** | Ref. ppb | OD |
|---:|---|---|---|---|---|---:|---:|---:|---|
| 4 | 10.2 ± 15 | 32.5 ± 15 | 32.5 ± 15 | 65.5 ± 10 | 65.5 ± 10 | 0.82 | **112** | | oily |
| 22 | 10.5 ± 13 | 62.5 ± 15 | 72.3 ± 13 | 88.7 ± 7 | 88.7 ± 7 | 0.67 | **64** | 0.07 (w) | fatty |
| 37 | 22.5 ± 17 | 54.8 ± 21 | 54.8 ± 21 | 75.6 ± 15 | 75.6 ± 15 | 0.72 | **89** | | painty |
| 60 | 55.4 ± 13 | 55.4 ± 13 | 68.5 ± 25 | 78.5 ± 13 | 88.9 ± 9 | 0.99 | **64** | | fragrant |

### Table 2. "Correlation coefficients between selected organics, matrix viscosity and temperature"

Footnote: "DOT values were correlated with viscosity at temperatures of 4 C, 22 C, 37 C, and 60 C."
So **each correlation is over four points**.

| Compound | Viscosity | Temperature |
|---|---:|---:|
| Pentanal | 0.69 | -0.99 |
| Hexanal | 0.92 | -0.87 |
| t-2-hexenal | 0.92 | -0.95 |
| Heptanal | 0.89 | -0.97 |
| t-2-octenal | 0.87 | -0.97 |
| t,t-2,4-decadienal | 0.85 | -0.69 |

### Table 3. "Retention time and peak responses for selected compounds at different concentrations"

Footnote 1: "Retention time ± standard deviation." Footnote *: "Detector response (area counts) ±
standard deviation: the symbol '*' indicates **failure to integrate peak response because organic
concentration was below the detection limit of the instrument**." The concentration headings are
printed in square brackets and are the same A-E ladders as Table 1; several bracket glyphs are
OCR-mangled (`~401` for `[40]`, `r 1001` for `[100]`, `"l51` for `[75]`) but the ladder is
recoverable from Table 1 and is written out below. **These are FID area counts, not
concentrations.**

| Organic | Retention time ± SD (min) | conc -> | | | | |
|---|---|---|---|---|---|---|
| Pentanal | 25.24 ± 0.03 | [5] * | [10] * | [25] * | [35] 6.9 ± 0.21 | [60] 7.5 ± 0.31 |
| Hexanal | 29.63 ± 0.23 | [20] * | [40] * | [70] 3.0 ± 0.5 | [100] 5.9 ± 0.19 | [150] 9.8 ± 0.4 [?] |
| t-2-hexenal | 30.54 ± 0.18 | [5] * | [10] * | [30] 2.6 ± 0.1 | [50] 6.1 ± 0.01 | [100] 8.8 ± 0.02 |
| Heptanal | 34.74 ± 0.41 | [5] 2.2 ± 0.03 | [10] 5.4 ± 0.02 | [20] 6.6 ± 0.04 | [50] 7.1 ± 0.02 | [100] 8.1 ± 0.01 |
| t-2-octenal | 38.72 ± 0.08 | [10] 1.8 ± 0.05 | [30] 3.8 ± 0.02 | [50] 4.0 ± 0.01 | [100] 4.9 ± 0.03 | [150] 5.0 ± 0.07 |
| t,t-2,4-decadienal | 48.31 ± 0.20 | [5] * | [10] * | [50] * | [75] 3.3 ± 0.01 | [100] 3.4 ± 0.08 |

The pentanal [5]/[10]/[25] cells: the scan shows the `*` symbol under [35] and [60] having values
and the first three columns blank-or-starred; the prose confirms the reading — "The minimum
concentration for GC detection of **pentanal was 35 ppb; hexanal was 70 ppb; t-2-hexenal was
30 ppb; and t,t-2,4-decadienal was 75 ppb**. GC detection occurred at the lowest concentration
tested for the remaining organics" (i.e. heptanal at 5 ppb and t-2-octenal at 10 ppb).

### Numbers printed in the running text (and NOT measured here)

| quantity | value | where | whose measurement |
|---|---|---|---|
| paraffin oil at 22 C: hexanal | 120 ppb | p. 234 | **Guadagni et al. 1972**, quoted |
| paraffin oil at 22 C: heptanal | 250 ppb | p. 234 | Guadagni et al. 1972, quoted |
| paraffin oil at 22 C: 2,4-decadienal | 135 ppb | p. 234 | Guadagni et al. 1972, quoted |
| meat slurries: methional | 6 100 ppb | p. 234 | **Wick et al. 1967**, quoted |
| meat slurries: phenylacetaldehyde | 940 ppb | p. 234 | Wick et al. 1967, quoted |
| meat slurries: 1-nonanal | 7 600 ppb | p. 234 | Wick et al. 1967, quoted |
| training suprathreshold doses | pentanal 50, hexanal 220, t-2-hexenal 120, heptanal 100, t-2-octenal 150, decadienal 110 **ppm** | Methods | this paper (not thresholds) |
| GC/detector-response linearity r | pentanal 0.985, hexanal 0.980, t-2-hexenal 0.946, heptanal 0.588, t-2-octenal 0.708, decadienal 0.780 | p. 241, describing Fig. 6 | this paper |
| DOT range over the whole study | "22 to 170 ppb" | Abstract | this paper |
| rank order of DOT | pentanal < hexanal < heptanal < t,t-2,4-decadienal < t-2-hexenal < t-2-octenal | Abstract | this paper |

**Viscosity: FIGURE-ONLY.** The Brookfield readings appear only as a dashed curve on Fig. 3 against
a right-hand axis in centipoise. Per house rule they are not typed as numbers. The same figure and
Fig. 4 are the only place the DOT-vs-temperature trend is drawn; the underlying values are Table 1
and are typed above.

### Arithmetic on the printed thresholds (all mine)

**1. The temperature effect within this matrix, per compound (DOT at 4 C over DOT at 60 C).**
pentanal 47/22 = **2.14x**; hexanal 90/38 = **2.37x**; t-2-hexenal 170/60 = **2.83x**; heptanal
108/50 = **2.16x**; t-2-octenal 140/81 = **1.73x**; decadienal 112/64 = **1.75x**. **Range 1.73x
to 2.83x over a 56 C span**, and **non-monotone in two of six**: hexanal rises from 34 (37 C) to 38
(60 C), and decadienal rises from 64 (22 C) to 89 (37 C) before falling back to 64. The abstract's
claim that DOT is "negatively correlated with temperature" is true on average and false in detail
for those two compounds. This is the number a modeller wants when asking how much a serving
temperature can move an odour-activity ratio inside one matrix: **a factor of about two, not a
decade.**

**2. The compound-to-compound spread at fixed temperature.** At 4 C: 47 to 170 ppb = **3.6x**. At
22 C: 41 to 109 = **2.7x**. At 37 C: 34 to 105 = **3.1x**. At 60 C: 22 to 81 = **3.7x**. So inside
this matrix the six compounds sit within a factor of four of each other **at every temperature** —
whereas their quoted aqueous values (Table 1 "Ref." column) span 0.07 to 12 ppb, a factor of
**171x**. The gelatin gel compresses the six compounds onto each other, and that compression is the
real finding of the paper.

**3. The gelatin/water ratios at 22 C (mine, and CROSS-METHOD in every row).** Dividing the 22 C DOT
by the Table 1 "Ref." value: pentanal 41/12 = **3.4x**; hexanal 58/4.5 = **12.9x**; heptanal
79/3 = **26.3x**; t-2-hexenal 109/3 = **36.3x**; t-2-octenal 109/3 = **36.3x**; t,t-2,4-decadienal
64/0.07 = **914x**. These reproduce `k2_matrix_and_thresholds.md` sec. A.2 exactly. **Every one of
these six ratios divides a 75 %-uncorrected sniffing threshold by a Guadagni forced-choice
threshold measured thirty years earlier in another laboratory.** They are not same-method pairs and
must never be reported as if they were. The 914x on decadienal in particular rests on a single
quoted `0.07 ppb` from a 1972 potato-chip paper.

**4. The unsaturation contrast the repository actually uses.** t-2-hexenal / hexanal at 22 C in the
gel = 109/58 = **1.88x**; the same pair in the quoted water column = 3/4.5 = **0.667x**. The ratio
of ratios = 1.88/0.667 = **2.82x** (the code carries **2.81**, the same number to rounding, computed
as 36.3/12.9). What makes it interesting is that it runs **against** hydrophobicity — t-2-hexenal is
the more water-soluble of the pair and should partition into the headspace *more*, not less — which
is why `parameters_matrix.py` reads it as adduction rather than partition. **But note what it rests
on**: both legs of the ratio of ratios carry a Guadagni water value, and the "cancelling" of the
75 %-criterion offset only works if that offset is the same multiplicative factor for both
compounds, which is assumed, not shown.

**5. Chain length among the saturated aldehydes.** C5 -> C6 -> C7 at 22 C: 41 -> 58 -> 79 ppb, i.e.
**+18 and +21 ppb per CH2, or 1.41x and 1.36x per CH2**. The trend holds at all four temperatures
(4 C: 47/90/108; 37 C: 34/34/62; 60 C: 22/38/50), though pentanal and hexanal are **equal** at 37 C.
Mean of the six per-CH2 ratios across the four temperatures = **1.35x per CH2 (mine)**. The paper's
own explanation is decreasing volatility with chain length.

**6. The sensory panel beat the instrument, and by how much.** Pentanal's DOT at 60 C is **22 ppb**
against a minimum GC-detectable **35 ppb** — the nose is **1.6x more sensitive** than the FID here.
Same comparison for hexanal at 37 C: DOT 34 vs GC floor 70 = **2.1x**. For t-2-hexenal at 60 C: DOT
60 vs GC floor 30 — here the **instrument wins by 2x**. Mixed, and the paper only reports the case
that favours the nose.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Of this paper's six compounds, **two are
keyed**: `hexanal` and `heptanal`. `e_2_octenal` is keyed and is the (E)-isomer, which is what
"t-2-octenal" means, so that is a third — but confirm the naming decision before binding them.
**Absent from the registry: pentanal, t-2-hexenal (trans-2-hexenal) and t,t-2,4-decadienal
((E,E)-2,4-decadienal).** Three of the six compounds in the cleanest matrix threshold set in the
corpus have no registry id. Every row below shares: **3 % w/v gelatin in distilled water (30 g
protein/L), lipid-free, no sugar/salt/buffer, pH not stated, dosed at 22 C from an aqueous stock and
held 18 h at 4 C with NO thermal step after dosing, 100 mL in a 250 mL stoppered flask, single
sample vs a gelatin control, ascending series, linear least-squares fit to 75 % detection
UNCORRECTED for chance, 16 panellists, 3 replicates, room at 22 C / 50 % RH.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| pentanal DOT in gelatin_3pct | 47 | ppb (µg/L) | 4 C, CD 0.80 | Table 1 p. 236 | **threshold** |
| pentanal DOT in gelatin_3pct | 41 | ppb | 22 C, CD 0.72 | Table 1 | **threshold** |
| pentanal DOT in gelatin_3pct | 34 | ppb | 37 C, CD 0.64 | Table 1 | **threshold** |
| pentanal DOT in gelatin_3pct | 22 | ppb | 60 C, CD 0.67 | Table 1 | **threshold** |
| hexanal DOT in gelatin_3pct | 90 | ppb | 4 C, CD 0.79 | Table 1 | **threshold** |
| hexanal DOT in gelatin_3pct | 58 | ppb | 22 C, CD 0.51 | Table 1 | **threshold** (weak fit, Flags 2) |
| hexanal DOT in gelatin_3pct | 34 | ppb | 37 C, CD 0.49 | Table 1 | **threshold** (weakest fit in the table, Flags 2) |
| hexanal DOT in gelatin_3pct | 38 | ppb | 60 C, CD 0.87 | Table 1 | **threshold** (non-monotone vs 37 C) |
| t-2-hexenal DOT in gelatin_3pct | 170 | ppb | 4 C, CD 0.75 | Table 1 | **threshold** (the highest in the paper) |
| t-2-hexenal DOT in gelatin_3pct | 109 | ppb | 22 C, CD 0.86 | Table 1 | **threshold** |
| t-2-hexenal DOT in gelatin_3pct | 79 | ppb | 37 C, CD 0.82 | Table 1 | **threshold** |
| t-2-hexenal DOT in gelatin_3pct | 60 | ppb | 60 C, CD 0.72 | Table 1 | **threshold** |
| heptanal DOT in gelatin_3pct | 108 | ppb | 4 C, CD 0.89 | Table 1 | **threshold** |
| heptanal DOT in gelatin_3pct | 79 | ppb | 22 C, CD 0.57 | Table 1 | **threshold** (weak fit) |
| heptanal DOT in gelatin_3pct | 62 | ppb | 37 C, CD 0.77 | Table 1 | **threshold** |
| heptanal DOT in gelatin_3pct | 50 | ppb | 60 C, CD 0.76 | Table 1 | **threshold** |
| t-2-octenal DOT in gelatin_3pct | 140 | ppb | 4 C, CD 0.96 | Table 1 p. 237 | **threshold** |
| t-2-octenal DOT in gelatin_3pct | 109 | ppb | 22 C, CD 0.85 | Table 1 | **threshold** |
| t-2-octenal DOT in gelatin_3pct | 105 | ppb | 37 C, CD 0.89 | Table 1 | **threshold** |
| t-2-octenal DOT in gelatin_3pct | 81 | ppb | 60 C, CD 0.86 | Table 1 | **threshold** |
| t,t-2,4-decadienal DOT in gelatin_3pct | 112 | ppb | 4 C, CD 0.82 | Table 1 | **threshold** |
| t,t-2,4-decadienal DOT in gelatin_3pct | 64 | ppb | 22 C, CD 0.67 | Table 1 | **threshold** |
| t,t-2,4-decadienal DOT in gelatin_3pct | 89 | ppb | 37 C, CD 0.72 | Table 1 | **threshold** (non-monotone: HIGHER than 22 C) |
| t,t-2,4-decadienal DOT in gelatin_3pct | 64 | ppb | 60 C, CD 0.99 | Table 1 | **threshold** |
| gelatin protein loading | 30 | g/L | 3 g gelatin in 100 mL distilled water | Methods p. 231 | level_only (a stated formulation, not an assay) |
| lipid in the matrix | 0 | % | "devoid of lipid components" | Methods p. 231 | level_only |
| panel size | 16 | panellists | 2 M / 14 F, 21-32 y, non-smokers, 3 replicates | Methods p. 231 | level_only |
| DOT temperature effect, 4 C / 60 C | 1.73 to 2.83 | x | within this matrix, per compound | derived from Table 1 (mine) | within_study_ratio |
| DOT chain-length effect, saturated C5-C7 | 1.35 | x per CH2 (mean of six pairs) | 4-60 C | derived from Table 1 (mine) | within_study_ratio |
| t-2-hexenal / hexanal DOT in the gel at 22 C | 1.88 | x | 22 C, same panel, same method | derived from Table 1 (mine) | within_study_ratio |
| unsaturation penalty (gel/water on t-2-hexenal over gel/water on hexanal) | 2.82 | x | 22 C | derived (mine); the code carries 2.81 | within_study_ratio (**cross-method in both legs**, Flags 3) |
| gelatin/water ratio at 22 C, six compounds | 3.4 / 12.9 / 26.3 / 36.3 / 36.3 / 914 | x | pentanal / hexanal / heptanal / t-2-hexenal / t-2-octenal / decadienal | derived (mine) from Table 1 DOT and Ref. columns | within_study_ratio (**cross-study AND cross-method**, Flags 3) |
| paraffin oil at 22 C: hexanal / heptanal / decadienal | 120 / 250 / 135 | ppb | quoted, method not restated | p. 234 (Guadagni 1972) | **threshold** — but `quoted_second_hand_by_vega1994`; not measured here |
| water: pentanal / hexanal / heptanal / t-2-hexenal / t-2-octenal / decadienal | 12 / 4.5 / 3 / 3 / 3 / 0.07 | ppb | quoted, method not restated | Table 1 "Ref." column (Guadagni 1963/1972) | **threshold** — `quoted_second_hand_by_vega1994`; not measured here |
| meat slurries: methional / phenylacetaldehyde / 1-nonanal | 6 100 / 940 / 7 600 | ppb | quoted, matrix described only as "meat slurries" | p. 234 (Wick 1967) | **threshold** — second-hand, and **the same three numbers are called "in beef" by Brewer 1995**, Flags 7 |
| correlation of DOT with viscosity | 0.69 / 0.92 / 0.92 / 0.89 / 0.87 / 0.85 | r, over 4 points | pentanal / hexanal / t-2-hexenal / heptanal / t-2-octenal / decadienal | Table 2 p. 238 | within_study_ratio (Flags 4) |
| correlation of DOT with temperature | -0.99 / -0.87 / -0.95 / -0.97 / -0.97 / -0.69 | r, over 4 points | same order | Table 2 | within_study_ratio (Flags 4) |
| FID area counts vs concentration | 1.8 to 9.8 | area counts | 60 C platen, 22 mL vial | Table 3 p. 240 | level_only (**no calibration to concentration exists in this paper**) |
| retention times on DB-5 | 25.24 / 29.63 / 30.54 / 34.74 / 38.72 / 48.31 ± SD | min | 60 m DB-5, 35 C ramp to 280 C | Table 3 | level_only |
| gel viscosity vs temperature | — | centipoise | 4 / 22 / 37 / 60 C, spindle 1 | Fig. 3 | **figure_only** |
| DOT vs temperature curves; DOT bar chart; detector response vs concentration | — | — | — | Figs. 3, 4, 6 | **figure_only** (the underlying DOTs are Table 1 and are typed above) |

### What these can and cannot be put next to

**(a) They are the `gelatin_3pct` rows and nothing else.** `matrix_oav.py` already carries all 24
with `criterion="75%_uncorrected"`, `thermal_step_after_dosing=False`,
`concentration_verified=False` and `cross_study_cross_method=False`. Every one of those four flags
is confirmed here against the primary text. **This dossier changes no value**; it moves the
provenance from "via `k2_matrix_and_thresholds.md` sec. A.2" to the printed table, and it adds the
CD column, which the code does not currently carry and should (Flags 2).

**(b) They cannot be moved to a protein pot, and the reason is stronger than the house rule.**
Gelatin is denatured collagen: essentially **no cysteine**, so no free thiol and no disulfide, and
its lysine content is low. The binding classes in `src/kinetic_core/matrix_sites.py`
(`saturated_aldehyde_amine`, `unsaturated_aldehyde_amine`, `hmf_thiol`, `hmf_amine`) all consume an
`amine` or `free_thiol` pool, and **`gelatin_3pct` has no entry in
`data/species/protein_matrices.yml`**, so it charges nothing. A gelatin threshold and a soy-isolate
threshold are not the same measurement made in two places; they are two different chemistries.

**(c) The one thing that DOES transport is the temperature slope, and only as a shape.** Within one
matrix, one panel and one method, four temperatures give a 1.7-2.8x fall from 4 C to 60 C for six
structurally similar aldehydes. That is a **within_study_ratio** with unusually good internal
control, and it is the corpus's only measurement of how much serving temperature moves a detection
threshold. It is still a gelatin-gel measurement, in which the fall is partly a melting gel (the
viscosity collapses between 4 and 22 C) and partly vapour pressure — the paper cannot separate the
two, and says so.

**(d) What CANNOT be transported at all.** The gelatin/water ratios (§3 arithmetic 3), because both
legs are cross-method; the paraffin-oil and meat-slurry rows, because Vega measured neither; and any
absolute headspace concentration, because none was ever measured.

## 5. Flags

1. **The criterion is 75 % detection, UNCORRECTED for chance, against a single control.** This is
   not a BET, not ASTM E679, not a 3-AFC and not a triangle test. A same-different judgement against
   one reference has a guess rate that is neither 1/2 nor 1/3 but depends on the panellist's
   criterion, and none of it is corrected. **Direction is known — these thresholds are
   systematically HIGHER than a chance-corrected 50 % forced-choice value would be — and the size is
   not measured anywhere in this paper.** Every ratio that divides one of these by a Guadagni value
   inherits that unmeasured offset.
2. **The regression quality is printed per row and the code does not carry it.** The CD column runs
   from **0.49 to 0.99**. Six of the 24 rows have CD below 0.70 (**hexanal 22 C 0.51 and 37 C 0.49;
   heptanal 22 C 0.57; pentanal 37 C 0.64 and 60 C 0.67; decadienal 22 C 0.67**). A CD of 0.49 means
   the straight line explains under half the variance of the five points it is interpolating, and
   the reported threshold is an interpolation on that line. **`ThresholdRecord` has a `notes` field
   and no CD field; adding the CD per record is the cheapest available improvement to this table.**
   Note which rows are affected: the two most-used compounds in the repository (hexanal at 22 and
   37 C) are the two worst-fitted cells in the paper.
3. **Every ratio to water in this paper is cross-study and cross-method, including the one the
   repository fits on.** The "Ref. ppb" column is Guadagni et al. 1963 and 1972, a different
   laboratory, a different decade, a different method and (for the 1972 rows) a potato-chip paper.
   `unsat_penalty_gelatin = 2.81` is a ratio of two such ratios; its cross-method offsets cancel
   only under the assumption that the 75 %-uncorrected criterion inflates hexanal and t-2-hexenal by
   the same factor. That assumption is not tested here and cannot be tested from this paper.
   **k2_matrix_and_thresholds.md already flags this; it is repeated because the number is FIT and
   the flag is not in the code's `provenance` beyond a string.**
4. **The correlations in Table 2 are over four points.** An r of 0.92 on n = 4 is not evidence of
   much. More seriously, **viscosity and temperature are not independent in this design** — the
   viscosity was varied only by changing the temperature, so the two columns of Table 2 are two
   readings of the same single experimental axis. The paper's conclusion that DOT is "positively
   correlated with viscosity, and negatively correlated with temperature" contains one fact, not two.
   The prose itself half-concedes this: viscosity is flat from 22 to 60 C while DOT keeps falling.
5. **Nominal concentrations, never verified, with 18 h of opportunity to be wrong.** Everything is
   weighed in; nothing is measured at the sniff. 0.1 mL of an aqueous ppm stock into 100 mL of gel,
   held 18 h at 4 C in a stoppered flask with 150 mL of headspace, then warmed to as much as 60 C
   and equilibrated 15 min between panellists. The paper's own discussion notes that long-chain
   aliphatics "can migrate to the solution surface or adhere to glassware". **The heavier compounds
   (decadienal, t-2-octenal) are the ones most at risk, and they are the ones with the largest
   apparent matrix effect.** A referee could reasonably read part of the 914x on decadienal as loss.
6. **Physical constants in the Table 1 sub-headers are garbled or wrong and should not be used.**
   Printed: pentanal BP 103 C (correct); hexanal BP 131 C (correct); heptanal BP 153 C (correct);
   **t-2-hexenal BP unreadable** (the scan gives `P 4 7 T`; the true value is ~146-147 C, so the
   glyphs are probably "= 147 C" but this is not asserted); **t-2-octenal BP printed as 85 C**,
   which cannot be right for a C8 enal (true ~177 C) — possibly a reduced-pressure value with the
   pressure dropped; **t,t-2,4-decadienal BP printed as 115 C**, likewise impossible at atmospheric
   pressure for a C10 dienal (true ~250 C). Vapour pressures are printed as 1.5, 1, 10, 0.5, 2 and
   8 mm "at 20 C" and are mutually inconsistent with those boiling points. Similarly the GC oven
   ramp prints as "SC/min", an OCR failure for a digit (5 C/min is the obvious reading and is **not
   asserted**). **None of these enters section 4.**
7. **Three thresholds in this paper belong to Wick et al. 1967 and are labelled differently by the
   two Vega/Brewer papers.** Vega calls them "meat slurries"; the companion Brewer 1995 calls the
   same three numbers "in beef". Attributing them to either paper would be laundering, and the
   matrix label itself is unresolved. Carry them as Wick 1967 or not at all.
8. **No pH anywhere.** The gel is unbuffered gelatin in distilled water and the pH is never stated
   or measured. `parameters_matrix.py` records `ph=None` for `gelatin_3pct`, which is correct and
   is a real gap: `PH_ADDUCT_GATE_BELOW` cannot be evaluated for this matrix.
9. **No replicate-level dispersion on the thresholds.** The ± values in Table 1 are the standard
   deviations of the **percent-correct** at each concentration, not of the threshold. There is no
   confidence interval, no inter-panellist range and no ANOVA on any DOT. (Contrast the companion
   beef paper, which does print individual-panellist ranges — see `brewer1995_extraction.md`.) So
   the 24 thresholds enter the repository as **points with no measured uncertainty**, and
   `dispersion_scale="not_stated"` in the code is the honest encoding.
10. **What this paper does not contain**: any protein other than gelatin; any lipid; any sugar,
    salt or buffer; any pH; any headspace or in-matrix concentration measurement; any
    partition coefficient; any binding constant; any same-method aqueous comparison; any
    chance-corrected threshold; any thermal treatment after dosing; any statistical test; any
    supplementary material.
11. **What to request from the authors or a follow-up**: (i) the raw percent-correct data behind the
    24 regressions, so the thresholds can be refitted with a chance-corrected model and given an
    interval; (ii) a same-panel water arm — six aqueous thresholds by this exact method would turn
    all six gelatin/water ratios from cross-method into same-method in one afternoon and is the
    single highest-value missing experiment in this dossier; (iii) the viscosity values behind
    Fig. 3; (iv) the pH of the gel.
12. **Registry gaps against `data/keys/compounds.yml`**: `hexanal` and `heptanal` are keyed;
    `e_2_octenal` is keyed and covers t-2-octenal if the naming is confirmed. **`pentanal`,
    `t_2_hexenal` and `tt_2_4_decadienal` have no registry id**, yet all three are used as keys in
    `matrix_oav.py`'s `_VEGA_GELATIN`, in its `WATER_THRESHOLDS` and in
    `parameters_matrix.py`'s `unsat_penalty_gelatin` — i.e. **the threshold layer is keyed on
    compound names that the compound registry does not know**. That mismatch is the most
    actionable registry finding in this dossier.
