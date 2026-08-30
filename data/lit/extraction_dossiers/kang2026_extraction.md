# Kang et al. 2026 — COMPLETE TRANSCRIPTION (Wave K4c)

**Full extraction of every number in `data/articles/Kang2026.pdf`.**
Extraction date 2026-08-28. Every value below was re-read off 500–900 dpi page rasters;
the PDF text layer was used only for prose and was **not** trusted for units (see §0b).

---

## 0. PAPER IDENTITY — MATCHES THE WAVE BRIEF'S DOI, BUT NOT ITS TOPIC DESCRIPTION

| field | value |
|---|---|
| Authors | **Di Kang**ᵇ, Lin Jiangᵃ, Songjin Zhengᵇ, Yuan Huᵃ, Haifeng Wangᵃ, Teng Liᵃ, **Yuying Fu**\*ᵃ, **Yun Zhai**\*ᵃ |
| Title | **"Unlocking the potential of cysteine–xylose Maillard reaction intermediates as natural flavor precursors: a comprehensive study on flavor regulation, storage stability, and antioxidant activity"** |
| Journal | **Sustainable Food Technology, 2026, 4, 3239–3252** |
| DOI | **10.1039/D5FB00932D** ✔ matches the brief |
| Type | RSC Paper, Open Access (CC BY-NC 3.0), published online **11 March 2026**; received 5 Dec 2025, accepted 10 Mar 2026 |
| Affiliations | ᵃ School of Food Science and Biotechnology, Zhejiang Gongshang University, Hangzhou 310018, PR China; ᵇ China Tobacco Hebei Industrial Co. Ltd, Shijiazhuang 050051, PR China |
| PDF | 14 pages, born-digital (Arbortext/Distiller), no OCR. Text layer complete but **unit-corrupting** (§0b) |

### 0a. **TOPIC MISMATCH vs the wave brief — READ THIS FIRST**

The brief expected *"the cysteine–xylose Maillard-reaction-intermediates study … 100/120/140 °C
series … MFT/FFT/volatile values"*. The DOI and journal are correct, and there **is** a
100/120/140 °C series, but two expectations are wrong in ways that matter:

1. **The subject intermediates are TTCA and ARPs, not free Cys–Xyl.**
   The paper studies **2-threityl-thiazolidine-4-carboxylic acid (TTCA)** and the
   **xylose–cysteine Amadori rearrangement product (ARP)** — two *isolated, purified*
   Maillard intermediates. All thermal-reaction experiments start from **purified TTCA at
   10 mmol L⁻¹**, *not* from a cysteine + xylose mixture. There is no Cys–Xyl kinetic run
   anywhere in the paper except the 90 °C/40 min synthesis step (§2).
2. **⚠ THE 100/120/140 °C VOLATILE NUMBERS ARE NOT IN THIS PDF.**
   The temperature ladder's numeric table is **Table S4 (Supplementary Information)**, and
   the pH ladder's is **Table S5 (SI)**. Neither S4 nor S5 is in `Kang2026.pdf`. The only
   representation of the temperature/pH ladders in the main paper is **Fig. 1a and 1b**,
   two 3-D bar charts of *compound-class totals* (7 classes), which I have digitised in §5.
   **Per-compound values (MFT = 2-methyl-3-furanthiol, FFT = 2-furfurylthiol) at
   100/120/140 °C exist only in Table S4 and are UNAVAILABLE from this file.** Only three
   individual numbers from the ladders are quoted in the running text (§5c).
   **Recommendation: fetch the SI (https://doi.org/10.1039/d5fb00932d) before using this
   paper as a temperature-ladder anchor.** (Note: the paper's own Data-availability
   statement cites a broken SI DOI, `10.1039/d500932d`.)
3. The **only complete per-compound volatile table in the main text is Table 1**, which is a
   **single-temperature (120 °C, pH 7.0, 120 min) exogenous-amino-acid panel** — 9 systems ×
   48 compounds. That is transcribed in full in §4.

### 0b. **CRITICAL UNITS FINDING — THE PDF TEXT LAYER RENDERS `μ` AS `m`**

`pdftotext` on this file silently converts every micro sign to a Latin `m`. Verified against
the page rasters:

| text layer says | **printed on the page** | anchor |
|---|---|---|
| "Characteristic flavor compounds (**mg L⁻¹**)" | **μg L⁻¹** | Table 1 caption, p. 3244 |
| "reaching 11.039 **mg L⁻¹** at 140 °C" | **11.039 μg L⁻¹** | p. 3242, §3.1 |
| "thiazoles increased sharply to 10.952 **mg L⁻¹**" | **10.952 μg L⁻¹** | p. 3242, §3.1 |
| "pyrazines … at pH 8 and 0.838 **mg L⁻¹**" | **0.838 μg L⁻¹** | p. 3242, §3.1 |
| "3 **mL** internal standard … 0.018 **mg mL⁻¹**" | **3 μL … 0.018 μg μL⁻¹** | p. 3242, §2.7 |
| "100 **mmol L⁻¹** DPPH ethanol solution" | **100 μmol L⁻¹** | p. 3241, §2.6.3 |
| "CAR/PDMS/DVB, 75 **mm**" fibre; "0.25 **mm** film"; "3.5 **mm**" HPLC particles | **75 μm; 0.25 μm; 3.5 μm** | §2.7, §2.2 |

**Genuinely `m`, verified on raster:** "0.5 to 3.0 **mg mL⁻¹**" (antioxidant concentration
range, abstract); "100 **mmol L⁻¹**" (storage solution, §2.5); "5 **mmol L⁻¹**"
FeCl₂–phenanthroline (§2.6.1); "10 **mmol L⁻¹**" TTCA reaction solution (§2.4);
"0.0827 **mol L⁻¹**" Cys/Xyl synthesis (§2.2).

**Consequence for the repo: the entire Kang volatile scale is μg L⁻¹ (ppb-equivalent),
i.e. 1000× smaller than any transcription taken from the text layer.** Anything already
ingested at "mg L⁻¹" from this paper is wrong by three orders of magnitude.
Reagent concentrations in §2.6.2/2.6.4/2.6.5/2.6.6 that I did **not** individually raster-check
are marked **[μ?]** below.

### 0c. Provenance codes used here

- **[M]** measured value printed by the authors (table or running text)
- **[D]** digitised by me off the figure raster (uncertainty stated per figure)
- **[Z]** derived by me from the paper's own numbers, never printed
- **[S]** exists only in the Supplementary Information — **not in this PDF**

Replication throughout: **n = 3**, mean ± SD, Duncan's multiple range test, p < 0.05
(§2.8, p. 3241). Superscript letters in Table 1 are Duncan groupings **within a row**.

---

## 1. THE ONE-PARAGRAPH ANSWER TO THE BRIEF'S CRITICAL CHECK

**The volatile values ARE absolute concentrations, not peak areas.** The internal standard is
**1,2-dichlorobenzene, 3 μL of a 0.018 μg μL⁻¹ methanolic solution spiked into 3 g of sample**
(= 0.054 μg IS per 3 g, i.e. 18 μg kg⁻¹). Two quantitation tiers coexist and **must be kept
distinct**:

- **Tier A — true external standard curves with IS calibration.** Compounds marked with an
  **asterisk (\*)** in Table 1 (footnote, p. 3246: *"Compounds marked with an asterisk (\*)
  were quantified using the standard curve method with internal standard calibration"*).
  **These are genuine absolute concentrations in μg L⁻¹.** 2-Methyl-3-furanthiol and
  2-furfurylthiol are both starred → **MFT and FFT are Tier A.**
- **Tier B — internal-standard-only semi-quantitation.** Compounds without an asterisk: *"for
  those without reference standards, quantification was performed just using the internal
  standard method"* (§2.7, p. 3242). These assume a response factor of 1 vs
  1,2-dichlorobenzene. **Ratio-comparable across conditions for the same compound, but the
  absolute magnitude carries an unquantified (likely 2–10×) response-factor error.**

Neither tier is peak-area-only, so **no value here is ordinal-only**. But every *subtotal*,
*total*, and *class* number mixes Tier A and Tier B and therefore inherits Tier B's error;
class totals and Fig. 1a/1b should be treated as semi-quantitative.
**No LOD, LOQ, recovery, calibration range, or R² for any standard curve is reported** —
the curves are relegated to Table S3 **[S]**, which is not in this PDF.

---

## 2. SYSTEM COMPOSITION AND CONDITIONS — applies to everything below

### 2a. Preparation of the intermediates (§2.2, p. 3240) — measured, from refs 4 and 26
| variable | value as printed |
|---|---|
| Substrates | equimolar **L-cysteine + D-xylose, 0.0827 mol L⁻¹** each, in deionised water |
| pH | **7.4 ± 0.01**, set with NaOH (6 mol L⁻¹) |
| Heat treatment | **90 °C for 40 min**, then immediate ice-bath quench |
| Primary separation | ion-exchange resin **H⁺ Dowex 50WX4, 200–400 mesh**, eluent **0.2 mol L⁻¹ NH₄OH** |
| Final purification | semi-preparative RP-HPLC, **Xbridge Amide 4.6 × 150 mm, 3.5 μm** (Waters) |
| Confirmation | UPLC-ESI-MS (Waters Synapt MALDI Q-TOF) + NMR (Bruker DRX 400 MHz); spectra in **Fig. S1 / Table S1** **[S]** |
| **Purity** | **NEVER STATED NUMERICALLY.** No % purity, no yield, no mass balance anywhere in the main text. |

### 2b. Thermal reaction model systems (§2.4, p. 3240)
| block | conditions |
|---|---|
| **Temperature ladder (Fig. 1a, Table S4 [S])** | **TTCA 10 mmol L⁻¹, pH 7.0**, heated **100 / 120 / 140 °C for 120 min** |
| **pH ladder (Fig. 1b, Table S5 [S])** | **TTCA 10 mmol L⁻¹**, pH adjusted to **5.5 / 7 / 8** with NaOH (6 mol L⁻¹), heated **120 °C for 120 min** |
| **Exogenous amino-acid panel (Table 1, Fig. 1c–e)** | **TTCA 10 mmol L⁻¹ + equimolar amino acid**, **pH 7.0**, **120 °C for 120 min**. Amino acids used: **Gly, Ala, Met, Leu, Ser, Thr, Pro, Val, Phe, His, Lys, Glu, Asp** (13) + a **Cys** system and a **no-addition control** |
| Vessel / heating | **high-temperature/high-pressure-resistant glass reaction vessels in an oil bath**; sealed (pressure-rated → superheated water above 100 °C) |
| Quench | immediate ice bath; samples then held at **4 °C** until analysis |
| **pH drift (measured, p. 3242)** | post-reaction pH fell to **5.1, 4.9 and 4.7** from initial 5.5, 7 and 8 respectively **[M]** — i.e. **all three "pH" runs finish acidic**; the pH label is an *initial* pH only |

### 2c. Volatile analysis (§2.7, p. 3242)
| variable | value as printed (μ corrected) |
|---|---|
| Method | **HS-SPME-GC-MS** |
| Sample | **3 g** of reaction solution |
| Internal standard | **3 μL of 1,2-dichlorobenzene, 0.018 μg μL⁻¹ in methanol** |
| Vial | **20 mL** headspace vial |
| Fibre | **CAR/PDMS/DVB, 75 μm** (divinylbenzene/carboxen/PDMS) |
| Extraction/desorption | *"following a previously published method (ref 15)"* — **times and temperatures NOT restated; extraction time, incubation temperature and equilibration time are UNAVAILABLE from this paper** |
| Desorption | **250 °C, 10 min**, in the injection port |
| Instrument | Agilent **7890B** GC + **5977B** MSD |
| Column | **DB-Wax, 30 m × 0.25 mm i.d. × 0.25 μm** film |
| Oven program / carrier / flow / split | **NOT GIVEN** (deferred to ref 15) |
| ID basis | NIST 17 + WILEY 07 spectral match, **Kovats RI vs a C7–C30 n-alkane series**, plus literature RI |
| Quantitation | Tier A (\*) standard curve + IS calibration (**Table S3 [S]**); Tier B IS-only |

### 2d. Storage-stability design (§2.5, p. 3241)
| variable | value |
|---|---|
| Forms | **solid lyophilised powder** (in self-sealing bags, or open Petri dish in a desiccator) **and solution** |
| Solution concentration | **100 mmol L⁻¹** in deionised water |
| Temperature arm | **4 °C, 25 °C, 40 °C** (solution) |
| pH arm | **pH 5.5, 7, 9** (solution, room temperature) |
| Water-activity arm | **solid** powder over saturated salt solutions; a_w **0.113, 0.432, 0.843** (salt identities in **Table S2 [S]**; the a_w–salt mapping is NOT in this PDF) |
| Sampling | **every 10 days**, total **60 days** |
| Metrics | residual concentration, **A₄₂₀**, **A₂₉₄**, colour **L\*, a\*, b\*, ΔE** |
| Browning index method (§2.3) | **A₄₂₀** on a UV-vis spectrophotometer, sample diluted so that **A₄₂₀ falls in 0.05–1.0** |

### 2e. Antioxidant assays (§2.6) — reagent concentrations, μ-status flagged
Fe²⁺ chelation: 1 mL sample + 1.85 mL water + 0.05 mL FeCl₂–phenanthroline **5 mmol L⁻¹**
(raster-verified mmol), 28 °C, 30 s then 10 min, 3000 rpm × 5 min, **A₅₆₂**.
Reducing power: 1 mL sample + 2.5 mL phosphate buffer pH 6.6 + 2.5 mL 1% K₃Fe(CN)₆, 50 °C
20 min, + 2.5 mL 10% TCA, 3000 rpm 10 min, 2.5 mL supernatant + 2.5 mL 0.1% FeCl₃ + 2.5 mL
water, 10 min, **A₇₀₀**.
DPPH: 3 mL **100 μmol L⁻¹** DPPH in ethanol (raster-verified μ) + 1 mL sample, dark, RT,
25 min, **A₅₁₇**.
•OH: 1 mL **2.5 mmol L⁻¹ [μ?]** 1,10-phenanthroline + 2 mL phosphate buffer pH 7.40 + 1 mL
**2.5 mmol L⁻¹ [μ?]** FeSO₄ + 1 mL water + 1 mL **20 mmol L⁻¹ [μ?]** H₂O₂, 37 °C 60 min, **A₅₃₆**.
O₂•⁻: 9 mL **50 mmol L⁻¹ [μ?]** Tris–HCl pH 8.2 + 0.5 mL sample, 25 °C 20 min, + 0.04 mL
1,2,4-benzenetriol (0.45 mmol in 10 mL of **10 mmol L⁻¹ [μ?]** HCl), 3 min, + 1 drop
**10 mmol L⁻¹ [μ?]** ascorbic acid, 5 min, **A₃₂₇.₂**.
Lipid peroxidation: 1 mL lecithin **10.0 mg mL⁻¹** in 0.01 mol L⁻¹ phosphate pH 7.4 +
1 mL **400 mmol L⁻¹ [μ?]** FeCl₃ + 1 mL **400 mmol L⁻¹ [μ?]** ascorbic acid + 1 mL sample,
37 °C 60 min; + 2 mL TBA reagent (15 g TCA + 0.37 g TBA + 2 mL conc. HCl per 100 mL),
95 °C 15 min, ice 10 min, centrifuge, **A₅₃₂**.
**All five clearance-rate formulae are printed and transcribed in §7a.**

---

## 3. WHAT KINETICS THIS PAPER DOES AND DOES NOT CONTAIN — read before using it

**There is not one rate constant, half-life, activation energy, reaction order, or Arrhenius
fit anywhere in this paper.** Every thermal experiment is a **single-endpoint** measurement
at **120 min**. The storage experiments are 7-point time courses (0–60 d) but the authors
fit nothing to them and report no k.

Consequences:
1. **The 100/120/140 °C series is a 3-point endpoint ladder, not a kinetic run.** With one
   time point per temperature you can extract an *apparent* temperature dependence of
   accumulated product, but **not** a rate constant, because formation and destruction are
   convolved (2-acetylthiazole and furfural are both formed and consumed at 140 °C).
2. **Any Eₐ derived from Fig. 1a is an Eₐ of *accumulated yield at 120 min*, not of an
   elementary step.** I compute one in §5b **[Z]** purely so the orchestrator can see the
   magnitude; it must not be shipped as a mechanistic Eₐ.
3. **The storage curves (§6) are the only genuine time-series in the paper** and are the only
   place a first-order k can be honestly fitted **[Z]** — see §6c.

---

## 4. TABLE 1 — COMPLETE TRANSCRIPTION (the paper's only per-compound volatile table)

**Anchor: Table 1, pp. 3244–3246 (PDF pages 6–8), printed landscape across three pages.**
Printed title: *"Characteristic flavor compounds (**μg L⁻¹**) derived from the reactions
between different amino acids and TTCA"*.
**Conditions for every column: TTCA 10 mmol L⁻¹ + equimolar amino acid, pH 7.0, 120 °C,
120 min** (§2.4). Columns are grouped by the added amino acid's character:
**Neutral** = Control (no addition), Cys, Gly, Ala, Met; **Acidic** = Glu, Asp;
**Basic** = His, Lys.
`RIᵇ` = RI calculated on DB-Wax vs C7–C30 n-alkanes; `RIᶜ` = literature RI from flavornet.org /
NIST WebBook; `—` = no literature RI available; **`ND` = not detected**.
**`*` = Tier A (standard curve + IS calibration); no asterisk = Tier B (IS-only semi-quant).**
All values **mean ± SD, n = 3**; superscript letters = Duncan groups **within the row**.
**Every cell below was read off the 200 dpi page raster. No cell was unreadable.**

### 4a. Sulfur-containing compounds — Thiols and sulfides (13 rows)

| Compound | RIᵇ | RIᶜ | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Dimethyl disulfide | 1085 | 1037 | ND | ND | ND | ND | **0.956 ± 0.013**ᵃ | ND | ND | ND | ND |
| 3-Mercapto-2-butanone\* | 1273 | 1283 | **1.798 ± 0.025**ᵃ | ND | ND | ND | ND | ND | ND | ND | ND |
| **2-Methyl-3-furanthiol\*** | 1302 | 1305 | **1.388 ± 0.207**ᶠ | **4.083 ± 0.092**ᵈ | **11.569 ± 0.165**ᵇ | **0.329 ± 0.008**ᵍ | — | **7.327 ± 0.222**ᶜ | **0.344 ± 0.007**ᵍ | **2.536 ± 0.149**ᵉ | **18.752 ± 1.061**ᵃ |
| 3-Mercapto-2-pentanone\* | 1352 | 1343 | **1.298 ± 0.021**ᵃ | **0.438 ± 0.01**ᵇ | ND | ND | ND | ND | ND | ND | ND |
| Dimethyl trisulfide | 1392 | 1370 | ND | ND | ND | ND | **2.031 ± 0.028**ᵃ | ND | ND | ND | ND |
| **2-Furfurylthiol\*** | 1426 | 1402 | **4.107 ± 0.152**ᵈ | **11.591 ± 0.26**ᶜ | **15.254 ± 0.217**ᵇ | **2.133 ± 0.003**ᵉ | ND | **12.107 ± 0.366**ᶜ | **3.176 ± 0.061**ᵈ | **0.38 ± 0.022**ᶠ | **29.295 ± 1.658**ᵃ |
| 3-Thiophenethiol\* | 1564 | 1530 | **0.178 ± 0.011**ᵉ | **12.892 ± 0.289**ᵃ | **0.74 ± 0.011**ᶜ | ND | ND | **0.405 ± 0.012**ᵉ | ND | ND | **2.31 ± 0.131**ᵇ |
| 2-Thiophenemethanethiol\* | 1702 | 1713 | **0.196 ± 0.012**ᵈ | **0.483 ± 0.011**ᵃ | **0.267 ± 0.004**ᶜ | ND | ND | **0.207 ± 0.006**ᵈ | ND | ND | **0.41 ± 0.023**ᵇ |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 1732 | — | **0.098 ± 0.006**ᵈ | ND | **2.804 ± 0.04**ᵇ | ND | ND | **0.574 ± 0.017**ᶜ | ND | ND | **4.428 ± 0.251**ᵃ |
| 1,2,3-Trithiolane | 1794 | — | **0.378 ± 0.026**ᵃ | ND | ND | ND | ND | ND | ND | ND | ND |
| 3,3′-Dithiobis[2-methyl-furan] | 2139 | — | ND | ND | ND | ND | ND | **3.946 ± 0.119**ᵇ | ND | ND | **6.862 ± 0.388**ᵃ |
| Furfuryl sulfide | 2274 | 2223 | ND | ND | **0.159 ± 0.002**ᵃ | ND | ND | ND | **0.025 ± 0**ᵇ | ND | ND |
| Bis(2-furfuryl)sulfide | 2419 | — | **0.042 ± 0.008**ᵇ | **0.972 ± 0.022**ᵃ | — | ND | ND | ND | ND | ND | ND |
| **Subtotal** | | | **9.483 ± 0.199**ᵈ | **30.46 ± 0.684**ᵇ | **30.794 ± 0.438**ᵇ | **2.461 ± 0.011**ᶠ | **2.987 ± 0.042**ᵉᶠ | **24.566 ± 0.743**ᶜ | **3.545 ± 0.068**ᵉ | **2.916 ± 0.171**ᵉᶠ | **62.057 ± 3.512**ᵃ |
| **Kinds** | | | **9** | **6** | **6** | **2** | **2** | **6** | **4** | **2** | **5** |

**Two printing/consistency anomalies in this block [Z]:**
- **MFT under Met is printed as an em-dash `—`, not `ND`.** Everywhere else in the table `—`
  means "no literature RI". A `—` in a data cell is undefined. Also
  **Bis(2-furfuryl)sulfide / Gly is `—`.** Treat both as **missing, not zero.**
- **The printed "Kinds" row does not match the non-ND cell count in three columns [Z].**
  Counting non-ND, non-`—` cells: Control 11 (printed 9), Cys 6 ✔, Gly 6 ✔ (with `—`
  excluded), Ala 2 ✔, Met 2 ✔, Glu 6 ✔, Asp 4 ✔, His 2 ✔, **Lys 7 (printed 5)**.
  The Control and Lys "Kinds" values are wrong. Do not use the Kinds row as a count.

### 4b. Sulfur-containing compounds — Thiophenes (15 rows)

| Compound | RIᵇ | RIᶜ | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Thiophene\* | 1016 | 1022 | ND | ND | ND | ND | ND | ND | ND | ND | ND |
| 3-Methyl-thiophene | 1119 | 1106 | ND | ND | **10.192 ± 0.145**ᵃ | ND | ND | ND | ND | ND | **0.083 ± 0.005**ᵇ |
| 2,3-Dihydro-5-methyl-thiophene | 1129 | 1156 | **0 ± 0**ᵉ | **0.058 ± 0.001**ᵈ | **0.167 ± 0.002**ᵇ | ND | ND | **0.126 ± 0.004**ᶜ | ND | ND | **0.204 ± 0.012**ᵃ |
| 2-Methyl-thiophene\* | 1135 | 1095 | **1.128 ± 0.159**ᶜ | **0.237 ± 0.005**ᵈ | **0.147 ± 0.002**ᵈ | **0.098 ± 0.002**ᵈ | ND | **5.389 ± 0.163**ᵇ | **0.141 ± 0.003**ᵈ | **0.48 ± 0.028**ᵈ | **12.21 ± 0.691**ᵃ |
| 2-Ethyl-thiophene | 1179 | 1167 | ND | ND | ND | ND | ND | **0.062 ± 0.002**ᶜ | ND | ND | **0.083 ± 0.005**ᵇ |
| 2,5-Dimethyl-thiophene | 1216 | 1202 | ND | ND | **0.464 ± 0.007**ᵇ | **0.76 ± 0.019**ᵃ | ND | ND | ND | ND | **0.461 ± 0.026**ᵇ |
| 2,3-Dimethyl-thiophene | 1232 | 1212 | ND | ND | **0.038 ± 0.001**ᵇ | ND | ND | **0.207 ± 0.006**ᵃ | ND | ND | ND |
| Dihydro-2-methyl-3(2H)-thiophenone | 1511 | 1506 | ND | ND | **1.167 ± 0.017**ᵇ | ND | ND | ND | ND | ND | **1.77 ± 0.1**ᵃ |
| 2-Thiophenecarboxaldehyde\* | 1674 | 1679 | **2.196 ± 0.076**ᵇ | **1.06 ± 0.024**ᶜ | **0.541 ± 0.008**ᵉ | ND | **3.109 ± 0.043**ᵃ | ND | **0.125 ± 0.002**ᶠ | **0.078 ± 0.005**ᶠ | **0.695 ± 0.039**ᵈ |
| 5-Methyl-2-thiophenecarboxaldehyde\* | 1701 | 1785 | **3.649 ± 0.075**ᵃ | **2.535 ± 0.057**ᶜ | ND | ND | ND | ND | ND | ND | **2.808 ± 0.159**ᵇ |
| 2-Acetyl-3-methylthiophene | 1754 | 1760 | ND | ND | **0.488 ± 0.007**ᵃ | ND | ND | ND | ND | ND | **0.388 ± 0.022**ᵇ |
| Thieno[3,2-b]thiophene | 1868 | 1843 | **1.302 ± 0.099**ᵇ | **7.967 ± 0.179**ᵃ | **0.659 ± 0.009**ᵈ | **0.038 ± 0.001**ᶠ | **0.063 ± 0.001**ᶠ | **0.39 ± 0.012**ᵉ | **0.192 ± 0.004**ᶠ | **0.105 ± 0.006**ᶠ | **0.964 ± 0.055**ᶜ |
| 1-(2-Thienyl)-1-propanone | 1872 | — | ND | ND | **0.46 ± 0.007**ᵇ | ND | ND | ND | ND | ND | **0.584 ± 0.033**ᵃ |
| 2,5-Thiophenedicarboxaldehyde\* | 1911 | — | **0.724 ± 0.055**ᵈ | **1.261 ± 0.028**ᵃ | **1.015 ± 0.014**ᶜ | ND | **1.079 ± 0.015**ᵇ | **0.292 ± 0.009**ᵉ | ND | ND | ND |
| 2-Methylthieno[2,3-b]thiophene | 1947 | — | **0.19 ± 0.023**ᶜ | **1.886 ± 0.042**ᵃ | ND | ND | **0.012 ± 0.000**ᵈ | ND | ND | ND | **0.24 ± 0.014**ᵇ |
| **Subtotal** | | | **9.189 ± 0.392**ᶜ | **15.003 ± 0.337**ᵇ | **15.338 ± 0.218**ᵇ | **0.895 ± 0.022**ᶠ | **4.263 ± 0.059**ᵉ | **6.466 ± 0.196**ᵈ | **0.457 ± 0.009**ᶠ | **0.663 ± 0.039**ᶠ | **20.49 ± 1.16**ᵃ |
| **Kinds** | | | **6** | **7** | **13** | **3** | **4** | **7** | **3** | **3** | **14** |

Note the Control cell for 2,3-dihydro-5-methyl-thiophene is printed **`0 ± 0`ᵉ**, i.e. a
*measured* zero, distinct from `ND`. Same convention recurs in §4c and §4d.

### 4c. Sulfur-containing compounds — Thiazoles (9 rows)

| Compound | RIᵇ | RIᶜ | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Thiazole\* | 1240 | 1262 | **1.796 ± 0.128**ᵃ | **0.461 ± 0.01**ᵈ | **0.556 ± 0.008**ᶜᵈ | **0.131 ± 0.003**ᶠ | — | **0.58 ± 0.018**ᶜ | **0.353 ± 0.007**ᵉ | **0.772 ± 0.045**ᵇ | **0.506 ± 0.029**ᶜᵈ |
| 2-Methyl-thiazole\* | 1272 | 1278 | **0.194 ± 0.017**ᵈ | **0.579 ± 0.013**ᵇ | ND | ND | ND | **0 ± 0**ᵉ | ND | **0.612 ± 0.036**ᵃ | **0.277 ± 0.016**ᶜ |
| 2,5-Dimethyl-thiazole\* | 1301 | 1326 | **2.203 ± 0.143**ᵃ | **0.127 ± 0.003**ᵇ | ND | ND | ND | **0.188 ± 0.006**ᵇ | ND | ND | ND |
| 2-Ethyl-thiazole\* | 1319 | 1304 | **0.273 ± 0.008**ᵃ | **0.152 ± 0.003**ᵇ | ND | ND | ND | **0 ± 0**ᶜ | ND | ND | ND |
| 2,4,5-Trimethyl-thiazole\* | 1373 | 1390 | **3.199 ± 0.127**ᵃ | **0.928 ± 0.021**ᵇ | **0.38 ± 0.005**ᶜ | ND | ND | **0.442 ± 0.013**ᶜ | ND | ND | **0.209 ± 0.012**ᵈ |
| 4,5-Dimethyl-thiazole\* | 1378 | 1843 | ND | ND | ND | ND | ND | ND | ND | ND | ND |
| 5-Ethyl-2,4-dimethyl-thiazole\* | 1437 | — | **0.587 ± 0.014**ᵃ | **0.174 ± 0.004**ᵈ | **0.339 ± 0.005**ᶜ | ND | ND | ND | ND | ND | **0.38 ± 0.022**ᵇ |
| 2-Ethyl-4-methylthiazole\* | 1449 | 1410 | **1.278 ± 0.117**ᵃ | **0.284 ± 0.006**ᵇ | **0.219 ± 0.003**ᵇ | **0.014 ± 0.000**ᶜ | ND | ND | **0.018 ± 0**ᶜ | ND | ND |
| **2-Acetylthiazole** | 1646 | 1643 | **8.795 ± 0.213**ᵃ | **0.538 ± 0.012**ᶜ | ND | ND | **2.837 ± 0.04**ᵇ | ND | ND | ND | ND |
| **Subtotal** | | | **18.325 ± 0.744**ᵃ | **3.255 ± 0.073**ᵇ | **1.494 ± 0.021**ᶜ | **0.145 ± 0.004**ᵈ | **2.837 ± 0.04**ᵇ | **1.211 ± 0.037**ᶜ | **0.371 ± 0.007**ᵈ | **1.385 ± 0.081**ᶜ | **1.371 ± 0.078**ᶜ |
| **Kinds** | | | **8** | **9** | **4** | **2** | **1** | **3** | **2** | **2** | **4** |

**RI anomaly [Z]: 4,5-dimethyl-thiazole is given RIᶜ = 1843, identical to thieno[3,2-b]thiophene's
literature RI.** For a compound eluting at RIᵇ 1378 on DB-Wax this is impossible; **the literature
RI for 4,5-dimethyl-thiazole is a copy-paste error.** The row is all-ND so it carries no data.
**Thiazole under Met is `—` (undefined) again.**
**The Cys "Kinds" = 9 exceeds the 8 non-ND cells in that column [Z]** — another Kinds error.

### 4d. Totals for sulfur-containing compounds

| | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|
| **Total contents of S-containing compounds (μg L⁻¹)** | **36.997 ± 1.1**ᶜ | **48.718 ± 1.094**ᵇ | **47.626 ± 0.678**ᵇ | **3.501 ± 0.037**ᵍ | **10.087 ± 0.14**ᵉ | **32.243 ± 0.975**ᵈ | **4.374 ± 0.084**ᶠ | **4.963 ± 0.291**ᶠ | **83.918 ± 4.749**ᵃ |
| **Total kinds of S-containing compounds** | **23** | **22** | **23** | **6** | **6** | **16** | **9** | **7** | **23** |

**Arithmetic check [Z]:** subtotals sum to 36.997 / 48.718 / 47.626 / 3.501 / 10.087 / 32.243 /
4.373 / 4.964 / 83.918. **All nine reconcile to ≤ 0.001 μg L⁻¹ except Asp (4.373 vs printed
4.374) and His (4.964 vs printed 4.963) — pure rounding. The S-block totals are internally sound.**

### 4e. Nitrogen-containing heterocycles (10 rows)

| Compound | RIᵇ | RIᶜ | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pyridine | 1189 | 1202 | ND | ND | ND | ND | ND | ND | ND | ND | **0.033 ± 0.002**ᵃ |
| Methylpyrazine\* | 1214 | 1263 | ND | **0.343 ± 0.008**ᵈ | **1.623 ± 0.023**ᵇ | ND | **0.423 ± 0.006**ᵈ | **1.512 ± 0.046**ᵇ | ND | **0.692 ± 0.041**ᶜ | **4.055 ± 0.229**ᵃ |
| 2-Pyridinecarboxaldehyde | 1186 | — | ND | ND | ND | ND | ND | ND | ND | ND | **1.239 ± 0.07**ᵃ |
| Pyrazine\* | 1216 | 1209 | ND | **0.063 ± 0.001**ᶜᵈ | **0.208 ± 0.003**ᶜᵈ | **0.039 ± 0.001**ᵈ | **0.12 ± 0.002**ᶜᵈ | **0.119 ± 0.004**ᶜᵈ | **0.294 ± 0.006**ᵇᶜ | **0.498 ± 0.029**ᵇ | **5.699 ± 0.323**ᵃ |
| 3-Methyl-pyridine\* | 1307 | 1346 | ND | **0.076 ± 0.002**ᵃ | ND | ND | ND | ND | ND | ND | ND |
| 2,5-Dimethyl-pyrazine\* | 1320 | 1328 | ND | **0.679 ± 0.015**ᵇ | **0.423 ± 0.006**ᶜ | ND | **0.069 ± 0.001**ᵉ | **0.402 ± 0.012**ᶜ | ND | **0.129 ± 0.008**ᵈ | **0.927 ± 0.052**ᵃ |
| 2-(n-Propyl)-pyrazine | 1476 | 1428 | ND | ND | ND | ND | ND | ND | ND | ND | **0.378 ± 0.021**ᵃ |
| Pyrrole | 1499 | 1518 | ND | ND | **0.228 ± 0.003**ᵃ | ND | ND | **0.12 ± 0.004**ᵇ | ND | ND | ND |
| 2-Methyl-1H-pyrrole | 1534 | 1551 | ND | ND | **0.371 ± 0.005**ᵇ | **0.036 ± 0.001**ᵈ | ND | **0.139 ± 0.004**ᶜ | ND | ND | **0.411 ± 0.023**ᵃ |
| 1H-Pyrrole-2-carboxaldehyde | 2002 | 2028 | ND | ND | ND | ND | **0.836 ± 0.012**ᵃ | ND | ND | ND | ND |
| **Subtotal** | | | **ND** | **1.161 ± 0.026**ᵈ | **2.853 ± 0.041**ᵇ | **0.076 ± 0.002**ᵉ | **1.448 ± 0.02**ᵈ | **2.292 ± 0.069**ᶜ | **0.294 ± 0.006**ᵉ | **1.319 ± 0.077**ᵈ | **12.741 ± 0.721**ᵃ |
| **Kinds** | | | **0** | **4** | **5** | **2** | **2** | **5** | **1** | **3** | **7** |

**Arithmetic check [Z]:** Cys 1.161 ✔, Gly 2.853 ✔, Ala 0.075 (printed 0.076), Met 1.448 ✔,
Glu 2.292 ✔, Asp 0.294 ✔, His 1.319 ✔, Lys 12.742 (printed 12.741). **All reconcile.**
**RI anomaly [Z]: 2-pyridinecarboxaldehyde is listed at RIᵇ 1186, i.e. *before* pyridine
(1189) and methylpyrazine (1214), yet it is printed after them — the table is otherwise
sorted by ascending RIᵇ.** Row order is not RI-monotone here; check before machine-parsing.

### 4f. Oxygen-containing heterocycles (9 rows)

| Compound | RIᵇ | RIᶜ | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Furan\* | 797 | 798 | **0.138 ± 0.01**ᵃ | ND | ND | ND | ND | **0 ± 0**ᵇ | **0 ± 0**ᵇ | ND | ND |
| 2-Methyl-furan\* | 851 | 829 | **0.019 ± 0.005**ᶠ | ND | **6.884 ± 0.098**ᵇ | **0.505 ± 0.012**ᵉ | ND | **4.741 ± 0.143**ᶜ | **0.276 ± 0.005**ᵉᶠ | **1.589 ± 0.093**ᵈ | **8.225 ± 0.465**ᵃ |
| 2-Ethyl-furan | 913 | 949 | ND | ND | **0.041 ± 0.001**ᵇ | ND | ND | **0.042 ± 0.001**ᵇ | **0 ± 0**ᶜ | ND | **0.089 ± 0.005**ᵃ |
| **Furfural\*** | 1457 | 1460 | **5.793 ± 0.155**ᵃ | ND | **1.05 ± 0.015**ᵉ | **0.242 ± 0.006**ᵍ | **1.968 ± 0.027**ᶜ | **1.019 ± 0.031**ᵉ | **3.589 ± 0.069**ᵇ | **1.369 ± 0.08**ᵈ | **0.664 ± 0.038**ᶠ |
| 2-Acetylfuran | 1491 | 1500 | ND | ND | ND | ND | ND | **1.189 ± 0.036**ᵃ | **0.378 ± 0.007**ᵇ | ND | ND |
| 1-(2-Furanyl)-ethanone\* | 1497 | 1501 | **0.136 ± 0.077**ᵃ | ND | ND | ND | ND | ND | ND | ND | ND |
| 5-Methylfurfural | 1542 | 1567 | ND | ND | **1.34 ± 0.019**ᵃ | **0.703 ± 0.017**ᵇ | ND | ND | ND | ND | ND |
| 2(5H)-Furanone\* | 1748 | 1767 | **0.072 ± 0.005**ᵃ | ND | ND | ND | ND | ND | ND | ND | ND |
| **4-Hydroxy-5-methyl-3(2H)-furanone\*** (norfuraneol) | 2108 | 2124 | ND | ND | ND | ND | ND | **0.19 ± 0.006**ᵃ | ND | ND | ND |
| **Subtotal** | | | **6.158 ± 0.231**ᶜ | **ND** | **9.314 ± 0.133**ᵃ | **1.45 ± 0.036**ᵍ | **1.968 ± 0.027**ᶠ | **7.182 ± 0.217**ᵇ | **4.243 ± 0.082**ᵈ | **2.958 ± 0.174**ᵉ | **8.978 ± 0.508**ᵃ |
| **Kinds** | | | **5** | **0** | **4** | **3** | **1** | **5** | **3** | **2** | **3** |

**Arithmetic check [Z]:** Control 6.158 ✔, Gly 9.315 (printed 9.314), Ala 1.450 ✔, Met 1.968 ✔,
Glu 7.181 (printed 7.182), Asp 4.243 ✔, His 2.958 ✔, Lys 8.978 ✔. **All reconcile.**
**Note for the K4c sibling extraction:** the same 4-hydroxy-5-methyl-3(2H)-furanone (HMFO /
norfuraneol) that is the whole subject of `nakamura2020` appears here **exactly once**, at
**0.19 ± 0.006 μg L⁻¹ in the TTCA+Glu system only**. Its RIᵇ 2108 / RIᶜ 2124 is a usable
DB-Wax retention anchor.

### 4g. Grand totals (p. 3246)

| | Control | Cys | Gly | Ala | Met | Glu | Asp | His | Lys |
|---|---|---|---|---|---|---|---|---|---|
| **Total (μg L⁻¹)** | **43.155 ± 1.223**ᵈ | **49.878 ± 1.12**ᶜ | **59.793 ± 0.851**ᵇ | **5.027 ± 0.075**ᶠ | **13.503 ± 0.188**ᵉ | **41.717 ± 1.262**ᵈ | **8.911 ± 0.172**ᵉ | **9.24 ± 0.542**ᵉ | **105.637 ± 5.978**ᵃ |
| **Total kinds** | **28** | **26** | **32** | **11** | **9** | **26** | **13** | **12** | **33** |

**Arithmetic check on the grand totals [Z]** (S-total + N-subtotal + O-subtotal):
Control 36.997 + 0 + 6.158 = **43.155** ✔ · Cys 48.718 + 1.161 + 0 = **49.879** (printed 49.878) ✔
Gly 47.626 + 2.853 + 9.314 = **59.793** ✔ · Ala 3.501 + 0.076 + 1.45 = **5.027** ✔
Met 10.087 + 1.448 + 1.968 = **13.503** ✔ · Glu 32.243 + 2.292 + 7.182 = **41.717** ✔
Asp 4.374 + 0.294 + 4.243 = **8.911** ✔ · His 4.963 + 1.319 + 2.958 = **9.240** ✔
Lys 83.918 + 12.741 + 8.978 = **105.637** ✔
**All nine grand totals reconcile exactly. Table 1's arithmetic is sound; only the per-block
"Kinds" counts are unreliable.**

### 4h. Derived quantities the paper does not print **[Z]**

Fold-change vs the no-addition Control, for the two flavour-critical thiols:

| system | MFT (μg L⁻¹) | ×Control | FFT (μg L⁻¹) | ×Control | FFT/MFT ratio |
|---|---|---|---|---|---|
| Control | 1.388 | 1.00 | 4.107 | 1.00 | **2.96** |
| + Cys | 4.083 | 2.94 | 11.591 | 2.82 | **2.84** |
| + Gly | 11.569 | 8.34 | 15.254 | 3.71 | **1.32** |
| + Ala | 0.329 | 0.24 | 2.133 | 0.52 | **6.48** |
| + Met | — (undefined) | — | ND | — | — |
| + Glu | 7.327 | 5.28 | 12.107 | 2.95 | **1.65** |
| + Asp | 0.344 | 0.25 | 3.176 | 0.77 | **9.23** |
| + His | 2.536 | 1.83 | 0.380 | 0.09 | **0.15** |
| + Lys | 18.752 | **13.51** | 29.295 | **7.13** | **1.56** |

**The FFT/MFT ratio spans 0.15–9.2 across amino-acid co-substrates at fixed T, pH and time.**
If the repo carries any fixed FFT:MFT branching ratio, this table is direct evidence that it
is not a constant — it is set by the nitrogen co-substrate.

---

## 5. FIGURE 1a AND 1b — THE TEMPERATURE AND pH LADDERS **[D]**

**Anchor: Fig. 1, p. 3243 (PDF page 5).** Caption as printed: *"(a) effect of temperature
(100 °C, 120 °C, and 140 °C) at pH 7.0 on volatile compounds; (b) effect of pH (5.5, 7, and 8)
at 120 °C on volatile compounds; (c) browning intensity of TTCA systems with different added
amino acids; (d) contents and (e) types of characteristic flavor compounds from TTCA with Gly,
Lys, or Glu; and (f) formation pathways…"*

Both panels are **3-D-effect clustered column charts**, y-axis printed
**"Concentration of flavor compound（μg/L）"**. Seven x-categories in both, in this printed
order: **Thiols/Sulfides · Thiophene · Thiazoles · Nitrogen-containing heterocycles ·
Oxygen-containing heterocycles · Sulfur-containing compounds · Total content**.

**Digitisation method and its error:** bars were located by colour mask and the value read at
the **rear top edge** of the 3-D column; the front (true) top sits ≈ 8.4 px lower, which is
**≈ 1.1 μg L⁻¹ on Fig. 1a's 0–75 axis and ≈ 0.36 μg L⁻¹ on Fig. 1b's 0–25 axis**. Both raw and
offset-corrected values are given. Residual digitisation error after correction: **± 0.3 μg L⁻¹
(Fig 1a), ± 0.1 μg L⁻¹ (Fig 1b)**. No error bars are drawn on either panel.

### 5a. Fig. 1a — TTCA 10 mmol L⁻¹, pH 7.0, 120 min, at 100 / 120 / 140 °C

| class | 100 °C raw | 120 °C raw | 140 °C raw | **100 °C corr.** | **120 °C corr.** | **140 °C corr.** |
|---|---|---|---|---|---|---|
| Thiols/Sulfides | 9.02 | 12.54 | 25.61 | **7.9** | **11.4** | **24.5** |
| Thiophene | 3.01 | 9.93 | 18.16 | **1.9** | **8.8** | **17.1** |
| Thiazoles | 6.14 | 17.12 | 20.25 | **5.0** | **16.0** | **19.2** |
| N-containing heterocycles | 1.57 | 1.44 | 1.57 | **0.5** | **0.3** | **0.5** |
| O-containing heterocycles | 4.70 | 7.58 | 14.37 | **3.6** | **6.5** | **13.3** |
| "Sulfur-containing compounds" | 16.33 | 15.42 | 20.51 | **15.2** | **14.3** | **19.4** ⚠ |
| **Total content** | 18.82 | 43.38 | 74.48 | **17.7** | **42.3** | **73.4** |

**Internal consistency [Z]:** the five class bars sum to **18.9 / 43.0 / 74.6** against printed
totals of **17.7 / 42.3 / 73.4** — agreement to **≤ 1.3 μg L⁻¹ (≤ 7 %)**, which is within the
digitisation offset. **The class bars and the Total bar are mutually consistent; the digitisation
is validated.**

**⚠ FIGURE ERROR — Fig. 1a's "Sulfur-containing compounds" bars are wrong at 120 and 140 °C.**
That category must equal Thiols + Thiophenes + Thiazoles = **14.8 / 36.2 / 60.8**. The plotted
bars are **15.2 / 14.3 / 19.4** — correct at 100 °C, then **2.5× and 3.1× too low**. Worse, those
three plotted values (15.2/14.3/19.4) are **numerically almost identical to the corresponding
group in Fig. 1b (15.0/14.1/19.3, §5b)**, i.e. **the pH-experiment series appears to have been
pasted into the temperature panel.** **DO NOT USE Fig. 1a's "Sulfur-containing compounds"
group.** Use the sum of the three sulfur class bars instead. Everything else in Fig. 1a is
self-consistent.

### 5b. Fig. 1b — TTCA 10 mmol L⁻¹, 120 °C, 120 min, at initial pH 5.5 / 7 / 8

| class | pH 5.5 raw | pH 7 raw | pH 8 raw | **pH 5.5 corr.** | **pH 7 corr.** | **pH 8 corr.** |
|---|---|---|---|---|---|---|
| Thiols/Sulfides | 12.05 | 7.86 | 7.01 | **11.7** | **7.5** | **6.7** |
| Thiophene | 1.92 | 1.84 | 2.09 | **1.6** | **1.5** | **1.7** |
| Thiazoles | 1.20 | 5.47 | 11.41 | **0.8** | **5.1** | **11.1** |
| N-containing heterocycles | *not resolved* | *not resolved* | 1.37 | **≈ 0** | **≈ 0** | **1.0** |
| O-containing heterocycles | 5.77 | 3.80 | 1.58 | **5.4** | **3.4** | **1.2** |
| Sulfur-containing compounds | 15.30 | 14.44 | 19.70 | **14.9** | **14.1** | **19.3** |
| **Total content** | 20.68 | 17.82 | 20.73 | **20.3** | **17.5** | **20.4** |

*Not resolved* = the pH 5.5 and pH 7 nitrogen bars are below the colour-mask detection floor
(< ≈ 0.3 μg L⁻¹); the text states pyrazines were **detected only at pH 8**, so ≈ 0 is right.

**Internal consistency [Z]:** here the "Sulfur-containing compounds" bar **does** equal
thiols + thiophenes + thiazoles: **14.1 vs 14.9** (pH 5.5), **14.1 vs 14.1** (pH 7),
**19.5 vs 19.3** (pH 8). **Fig. 1b is internally sound** — which is what convicts Fig. 1a.

### 5c. The only three ladder numbers printed as text — and their cross-validation **[M]**

From §3.1, p. 3242 (μ corrected):
1. *"The content of **furfural** (the dominant furan) increased significantly, reaching
   **11.039 μg L⁻¹ at 140 °C**, which was **3.3 times that at 100 °C**."*
   → **[Z] furfural at 100 °C = 11.039 / 3.3 = 3.345 μg L⁻¹** (and, by interpolation on the
   O-heterocycle class bar, ≈ 6 μg L⁻¹ at 120 °C — but that is inference, not measurement).
   **Cross-check:** my digitised O-heterocycle class totals are **3.6 / 6.5 / 13.3 μg L⁻¹** at
   100/120/140 °C. Furfural alone accounting for 3.35 of 3.6 and 11.04 of 13.3 is exactly what
   *"the dominant furan"* means. **Independent confirmation of both the digitisation and the
   μg L⁻¹ unit.**
2. *"As the pH increased from 5.5 to 8, thiol and furan contents decreased, while **thiazoles
   increased sharply to 10.952 μg L⁻¹**."*
   **Cross-check:** my digitised pH-8 thiazole bar = **11.1 μg L⁻¹**. Agreement to 1.4 %.
3. *"**2-Acetylthiazole** (nutty/popcorn odor) increased significantly with pH, comprising
   **27.01 % of total thiazoles at pH 8**."*
   → **[Z] 2-acetylthiazole at pH 8, 120 °C = 0.2701 × 10.952 = 2.958 μg L⁻¹.**
4. *"**Pyrazines** were detected **only at pH 8** and **0.838 μg L⁻¹**."*
   **Cross-check:** my digitised pH-8 N-heterocycle class bar = **1.0 μg L⁻¹**; pyrazines
   0.838 of that, remainder pyridines/pyrroles. Consistent.
5. *"Under the acidic condition (pH 5.5), thiols were most abundant, followed by thiophenes;
   **no nitrogen-containing heterocycles were detected**."* **[M]**

**These five statements are the ONLY per-compound numbers available from the temperature and
pH ladders without the SI.** Concretely: **MFT and FFT at 100/120/140 °C are stated only
qualitatively** — *"The concentrations of 2-methyl-3-furanthiol and 2-furfurylthiol also
increased accordingly"* (p. 3242). **No MFT or FFT ladder value exists in this PDF.**

### 5d. Apparent temperature dependence **[Z] — NOT a mechanistic activation energy**

Fitting ln(yield at 120 min) vs 1/T over three points (373.15 / 393.15 / 413.15 K),
R = 8.314 J mol⁻¹ K⁻¹, using the offset-corrected Fig. 1a values:

| class | apparent Eₐ (kJ mol⁻¹) | R² | note |
|---|---|---|---|
| Thiols/Sulfides | **36.0** | 0.948 | monotone |
| Thiophene | **70.8** | 0.962 | steepest; strongest T-gating |
| Thiazoles | **43.6** | 0.870 | saturating between 120 and 140 °C |
| O-heterocycles | **41.8** | 0.993 | furfural-dominated, near-perfect line |
| **Total content** | **45.7** | 0.990 | |
| N-heterocycles | *undefined* | — | non-monotone (0.5 → 0.3 → 0.5) |

**These are yield-accumulation exponents at a single 120 min endpoint, not elementary-step
barriers.** They are quoted here only so the orchestrator can judge magnitude; three points,
no replicate scatter on the figure, and formation/consumption convolved. **Do not ship as Eₐ.**
The furfural anchor gives an independent **single-compound, text-sourced [M]** check:
3.345 → 11.039 μg L⁻¹ over 100 → 140 °C is **Eₐ,apparent = 38.3 kJ mol⁻¹ [Z]**, against
**41.8 kJ mol⁻¹** for the digitised O-heterocycle class it dominates. **Agreement to 8 % —
a second independent validation of the Fig. 1a digitisation, this one on the temperature
dependence rather than the level.**

---

## 6. STORAGE STABILITY — the only genuine time-series in the paper

### 6a. Numbers printed in the text and abstract **[M]**

| stress | 60-day loss, TTCA | 60-day loss, ARP | anchor |
|---|---|---|---|
| **40 °C** (solution, 100 mmol L⁻¹) | **7.06 %** | **12.17 %** | p. 3247 §3.3; abstract |
| **pH 9** (solution, room temp.) | **11.19 %** | **21.25 %** | p. 3247 §3.3; abstract |
| **a_w 0.843** (solid) | **35.77 %** | **60.61 %** | abstract, p. 3239 |
| **a_w 0.113** (solid) | **13.5 %** | *(not stated)* | p. 3247 §3.3 |
| 4 °C, 25 °C | *"highly stable"*, no number | *"highly stable"*, no number | p. 3247 |
| pH 5.5, pH 7 | *"maintained stable concentrations"* | same | p. 3247 |

Also **[M]**: *"TTCA's degradation rate in high humidity (a_w 0.843) is approximately **60 %
slower** than that of ARP"* — explicitly labelled by the authors as **"a semi-quantitative
estimate"** (p. 3247). Conclusion, p. 3250: TTCA *"can retain **over 90 %** of its content after
60 days of storage under ambient, neutral, and dry conditions."*

### 6b. Fig. 2 digitisation — all six residual-content panels **[D]**

**Anchor: Fig. 2, p. 3247 (PDF page 9).** All panels: y-axis **"Concentration（mmol/L）", 0–120**,
x-axis storage day **0, 10, 20, 30, 40, 50, 60**. Starting concentration 100 mmol L⁻¹ in every
panel, so **the numbers below double as % remaining**. Error bars (SD, n = 3) are drawn but are
small; typical half-width ≈ **1–3 mmol L⁻¹** (largest ≈ 3 on the a_w panels).
Digitisation error **± 1 mmol L⁻¹**.

**Fig. 2a — TTCA in solution vs temperature**
| day | 4 °C | 25 °C | 40 °C |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 100 | 99.5 | 99 |
| 20 | 100 | 99.5 | 98.5 |
| 30 | 100 | 99 | 97.5 |
| 40 | 100 | 98 | 96.5 |
| 50 | 99.5 | 96 | 95 |
| 60 | 99 | 95.5 | **92.94** (pinned to the printed 7.06 % loss) |

**Fig. 2b — ARP in solution vs temperature**
| day | 4 °C | 25 °C | 40 °C |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 99 | 98 | 97 |
| 20 | 98 | 97 | 95.5 |
| 30 | 97.5 | 96.5 | 94.5 |
| 40 | 97 | 95 | 92 |
| 50 | 96 | 95.5 | 90.5 |
| 60 | 95.5 | 94 | **87.83** (pinned to the printed 12.17 % loss) |

**Fig. 2e — TTCA in solution vs pH (room temperature)**
| day | pH 5.5 | pH 7 | pH 9 |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 99.5 | 99 | 97.5 |
| 20 | 99 | 98.5 | 96 |
| 30 | 98.5 | 98 | 95.5 |
| 40 | 97.5 | 98 | 92.5 |
| 50 | 96.5 | 97.5 | 90 |
| 60 | 95.5 | 97 | **88.81** (pinned to the printed 11.19 % loss) |

**Fig. 2f — ARP in solution vs pH (room temperature)**
| day | pH 5.5 | pH 7 | pH 9 |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 98.5 | 98.5 | 95 |
| 20 | 97.5 | 97.5 | 92 |
| 30 | 97 | 94.5 | 86.5 |
| 40 | 95.5 | 93.5 | 82 |
| 50 | 94.5 | 91.5 | 81 |
| 60 | 94 | 90.5 | **78.75** (pinned to the printed 21.25 % loss) |

**Fig. 2i — solid TTCA vs water activity**
| day | a_w 0.113 | a_w 0.432 | a_w 0.843 |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 96.5 | 97 | 91.5 |
| 20 | 94.5 | 91.5 | 88 |
| 30 | 91.5 | 88 | 80.5 |
| 40 | 89 | 81 | 73 |
| 50 | 86 | 76.5 | 68.5 |
| 60 | **86.5** (= 13.5 % loss ✔ printed) | 72.5 | **64.23** (= 35.77 % loss ✔ printed) |

**Fig. 2j — solid ARP vs water activity**
| day | a_w 0.113 | a_w 0.432 | a_w 0.843 |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 10 | 96.5 | 92.5 | 87.5 |
| 20 | 92 | 86.5 | 80 |
| 30 | 89 | 80.5 | 70.5 |
| 40 | 86 | 72.5 | 59.5 |
| 50 | 82 | 64 | 50 |
| 60 | 78.5 | 57.5 | **39.39** (= 60.61 % loss ✔ printed) |

**Three of the six panels' day-60 endpoints are independently pinned by numbers printed in the
text, and all three digitisations landed within 0.5 mmol L⁻¹ before pinning. The Fig. 2
digitisation is validated.**

### 6c. First-order rate constants fitted by me **[Z] — the paper prints none**

OLS of ln(C/C₀) vs t over all seven points, 0–60 d; C₀ = 100 mmol L⁻¹.
Units **d⁻¹**; s⁻¹ given for direct repo use.

| system | condition | **k (d⁻¹)** | k (s⁻¹) | R² | t½ (d) |
|---|---|---|---|---|---|
| TTCA soln | 4 °C | 1.435 × 10⁻⁴ | 1.661 × 10⁻⁹ | **0.615** ⚠ | 4831 |
| TTCA soln | 25 °C | 8.034 × 10⁻⁴ | 9.298 × 10⁻⁹ | 0.887 | 863 |
| TTCA soln | 40 °C | **1.152 × 10⁻³** | 1.334 × 10⁻⁸ | 0.955 | 602 |
| ARP soln | 4 °C | 7.498 × 10⁻⁴ | 8.678 × 10⁻⁹ | 0.988 | 924 |
| ARP soln | 25 °C | 9.219 × 10⁻⁴ | 1.067 × 10⁻⁸ | 0.924 | 752 |
| ARP soln | 40 °C | **2.019 × 10⁻³** | 2.337 × 10⁻⁸ | 0.986 | 343 |
| TTCA soln | pH 5.5 | 7.665 × 10⁻⁴ | 8.872 × 10⁻⁹ | 0.969 | 904 |
| TTCA soln | pH 7 | 4.536 × 10⁻⁴ | 5.250 × 10⁻⁹ | 0.943 | 1528 |
| TTCA soln | pH 9 | **1.976 × 10⁻³** | 2.287 × 10⁻⁸ | 0.978 | 351 |
| ARP soln | pH 5.5 | 1.033 × 10⁻³ | 1.196 × 10⁻⁸ | 0.986 | 671 |
| ARP soln | pH 7 | 1.746 × 10⁻³ | 2.020 × 10⁻⁸ | 0.986 | 397 |
| ARP soln | pH 9 | **4.109 × 10⁻³** | 4.756 × 10⁻⁸ | 0.975 | 169 |
| TTCA solid | a_w 0.113 | 2.591 × 10⁻³ | 2.999 × 10⁻⁸ | 0.968 | 268 |
| TTCA solid | a_w 0.432 | 5.577 × 10⁻³ | 6.455 × 10⁻⁸ | 0.989 | 124 |
| TTCA solid | a_w 0.843 | **7.479 × 10⁻³** | 8.656 × 10⁻⁸ | 0.994 | 93 |
| ARP solid | a_w 0.113 | 3.998 × 10⁻³ | 4.627 × 10⁻⁸ | 0.998 | 173 |
| ARP solid | a_w 0.432 | 9.191 × 10⁻³ | 1.064 × 10⁻⁷ | 0.986 | 75 |
| ARP solid | a_w 0.843 | **1.504 × 10⁻²** | 1.740 × 10⁻⁷ | 0.978 | 46 |

⚠ **TTCA at 4 °C has R² = 0.615** — that series is five identical 100.0 readings followed by
99.5 and 99.0, i.e. it is at the digitisation floor. **Treat TTCA/4 °C as "no measurable loss
in 60 d", not as k = 1.4 × 10⁻⁴ d⁻¹.**

**Caveats that make these [Z] and not [M]:**
- The authors never claim first-order behaviour; R² is high in most rows but so would be a
  zero-order fit over ≤ 36 % conversion.
- Maximum conversion in the whole dataset is **60.6 %** (ARP, a_w 0.843); most curves reach
  ≤ 15 %, so these are effectively **initial-rate** constants.
- **The "room temperature" of the pH arm is never given a number.** It cannot be equated to
  the 25 °C of the temperature arm: TTCA at pH 7 gives k = 4.54 × 10⁻⁴ d⁻¹ while TTCA at
  25 °C gives 8.03 × 10⁻⁴ d⁻¹, a 1.8× discrepancy for what should be the same condition
  (the temperature arm's pH is also never stated). **The two solution arms are not
  cross-comparable.** ⚠
- Solid-state k's are for **lyophilised powder in a sealed desiccator**, a different phase
  entirely; do not mix them with the solution k's.

**[Z] Apparent Eₐ of storage degradation** from the three-temperature k's
(277.15 / 298.15 / 313.15 K): **TTCA 43.1 kJ mol⁻¹ (R² 0.947)**, **ARP 18.7 kJ mol⁻¹
(R² 0.808)**. Both are low — consistent with the authors' attribution to water autoionisation
rather than a covalent bond-breaking step — and both rest on three points with no error bars
on k. **The TTCA value inherits the bad 4 °C fit (R² 0.615) and should carry a wide prior.**

**[Z] a_w sensitivity**, ln k vs a_w over the three water activities:
**TTCA d(ln k)/d a_w = 1.42 (R² 0.898)**, **ARP d(ln k)/d a_w = 1.79 (R² 0.952)** — i.e.
**each 0.1 unit of a_w costs ≈ 15 % (TTCA) / 20 % (ARP) in shelf life**, over 0.113–0.843.
This is the single most repo-relevant derived quantity in the paper, and it is a **clean,
7-point, 3-level, 2-compound water-activity series** — directly comparable in *design* to
Bell 1995's PVP a_w series, though on intermediate *disappearance* rather than pigment
*formation*.

### 6d. Fig. 2c, 2d, 2g, 2h — browning during storage **[D], coarse**

Twin-axis panels: **left axis A₂₉₄ (dashed series), right axis A₄₂₀ (solid series)**.
⚠ **AXIS LABELLING ERROR:** the left axis in all four panels is printed
**"Concentration（mmol/L）"** although it plots **A₂₉₄, a dimensionless absorbance**
(range 0.2–0.6 in c/d, 0–1 in g/h). The right axis is correctly labelled A₄₂₀.
**Do not machine-read the left axis label.**

| panel | series | day-0 → day-60 range (approx.) |
|---|---|---|
| 2c (TTCA, temperature) | A₂₉₄: 4/25/40 °C | all start ≈ 0.40; 4 °C and 25 °C end ≈ 0.44; **40 °C ends ≈ 0.52** |
| 2c | A₄₂₀: 4/25/40 °C | all start ≈ 0.035; end **0.062 / 0.065 / 0.072** |
| 2d (ARP, temperature) | A₂₉₄ | start ≈ 0.39; end **0.45 / 0.52 / 0.58** |
| 2d | A₄₂₀ | start ≈ 0.036; end **0.062 / 0.068 / 0.078** |
| 2g (TTCA, pH) | A₂₉₄ | start ≈ 0.40; **pH 5.5 stays ≈ 0.40–0.42 (flat/slightly down)**; pH 7 ends ≈ 0.60; pH 9 ends ≈ 0.68 |
| 2g | A₄₂₀ | start ≈ 0.04; end ≈ **0.10 / 0.12 / 0.135** |
| 2h (ARP, pH) | A₂₉₄ | start ≈ 0.40; end ≈ **0.55 / 0.70 / 0.90** |
| 2h | A₄₂₀ | start ≈ 0.045; end ≈ **0.115 / 0.145 / 0.16** |

**[M] The authors' own reading of 2g (p. 3247):** at pH 5.5, TTCA concentration falls slightly
while **A₄₂₀ rises and A₂₉₄ *falls*** — opposite directions — which they interpret as
*"a transformation distinct from degradation, potentially due to acid-catalyzed hydrolysis or
reversion to N-xylosylamine."* **This is the paper's one genuinely mechanistic observation
that carries data behind it**, and it means **TTCA loss at acid pH is not all going to
browning.** Any model that equates intermediate disappearance with pigment formation is wrong
below pH ≈ 6 on this evidence.

---

## 7. TABLE 2 — SOLID-STATE COLOUR AFTER 60 DAYS **[M]**

**Anchor: Table 2, p. 3248 (PDF page 10).** Printed title: *"Variations in TTCA (a) and ARP (b)
contents in solid samples over 60 days of storage at different water activities"* — note the
title says *contents* but the table reports **colour coordinates**, not contents.
Values mean ± SD, n = 3; letters = Duncan groups within the row. **ΔE has no ± because it is
computed from the means. Control = the a_w-unexposed reference (never defined further).**

| MRI | Indicator | Control | a_w 0.113 | a_w 0.432 | a_w 0.843 |
|---|---|---|---|---|---|
| **TTCA** | L\* | **88.39 ± 0.03**ᵃ | **81.36 ± 0.02**ᵇ | **77.39 ± 0.04**ᶜ | **72.34 ± 0.02**ᵈ |
| | a\* | **5.19 ± 0.01**ᵈ | **7.71 ± 0.03**ᶜ | **9.19 ± 0.02**ᵇ | **11.63 ± 0.04**ᵃ |
| | b\* | **11.13 ± 0.02**ᵈ | **16.78 ± 0.04**ᶜ | **22.36 ± 0.03**ᵇ | **27.63 ± 0.01**ᵃ |
| | **ΔE** | — | **9.36**ᶜ | **16.22**ᵇ | **23.90**ᵃ |
| **ARP** | L\* | **87.68 ± 0.03**ᵃ | **79.36 ± 0.03**ᵇ | **62.18 ± 0.03**ᶜ | **51.03 ± 0.03**ᵈ |
| | a\* | **6.21 ± 0.03**ᵈ | **11.39 ± 0.03**ᶜ | **17.68 ± 0.03**ᵇ | **21.39 ± 0.03**ᵃ |
| | b\* | **13.19 ± 0.03**ᵈ | **20.38 ± 0.03**ᶜ | **29.63 ± 0.03**ᵇ | **37.68 ± 0.03**ᵃ |
| | **ΔE** | — | **12.16**ᶜ | **32.44**ᵇ | **46.62**ᵃ |

**ΔE recomputation [Z]** as √(ΔL\*² + Δa\*² + Δb\*²) vs each compound's own Control:
TTCA → **9.36 / 16.22 / 23.90**; ARP → **12.16 / 32.44 / 46.62**.
**All six printed ΔE values reproduce exactly to 2 decimals.** Table 2 is arithmetically airtight.

**[Z] The colour data are a second, independent browning probe on the same a_w ladder** and
they scale differently from concentration loss: at a_w 0.843, ARP loses **1.69×** more mass
than TTCA (60.61 % vs 35.77 %) but develops **1.95×** more colour (ΔE 46.62 vs 23.90). ARP's
degradation is therefore **more chromophore-productive per mole lost** — the kind of
per-pathway yield asymmetry a colour model needs and rarely gets.

### 7a. The five printed formulae (§2.6) **[M]** — transcribed exactly

- Fe²⁺ chelation rate = (A₀ − A₁)/A₀ × 100 %  (A₀ blank, A₁ sample, at 562 nm)
- DPPH• clearance rate = [1 − (Aᵢ − A_j)/A₀] × 100 %  (Aᵢ sample+DPPH, A_j sample+ethanol, A₀ ethanol+DPPH, at 517 nm)
- •OH clearance rate = (A_S − A₁)/(A₀ − A₁) × 100 %  (A₁ = phenanthroline+FeSO₄+H₂O₂ control, A₀ = FeSO₄+phenanthroline blank, at 536 nm)
- O₂•⁻ clearance rate = (A₀ − A₁)/A₀ × 100 %  (A₃₂₇ of blank vs sample)
- Lipid-peroxidation inhibition rate = (A₀ − A₁)/A₀ × 100 %  (at 532 nm)

Reducing power has **no formula** — the raw **A₇₀₀** is reported directly.

---

## 8. FIGURE 1c — BROWNING INTENSITY OF THE 15 AMINO-ACID SYSTEMS **[D]**

**Anchor: Fig. 1c, p. 3243.** y-axis **A₄₂₀, 0.0–1.0**; x-axis "Reaction models", 15 bars.
Conditions: **TTCA 10 mmol L⁻¹ + equimolar amino acid, pH 7.0, 120 °C, 120 min**.
**No error bars are drawn.** Duncan letters printed above each bar.
Digitisation error **± 0.005 A₄₂₀** (bar tops located from the baseline upward at ≥ 60 % bar
width, so the significance letters above the bars do not contaminate the read).

| # | system | **A₄₂₀ [D]** | Duncan letter [M] |
|---|---|---|---|
| 1 | TTCA (control) | **0.106** | l |
| 2 | Cys-TTCA | **0.072** | m |
| 3 | Gly-TTCA | **0.491** | d |
| 4 | Ala-TTCA | **0.414** | e |
| 5 | Met-TTCA | **0.254** | i |
| 6 | Leu-TTCA | **0.323** | g |
| 7 | Ser-TTCA | **0.334** | g |
| 8 | Thr-TTCA | **0.366** | f |
| 9 | Pro-TTCA | **0.528** | c |
| 10 | Val-TTCA | **0.318** | gh |
| 11 | Phe-TTCA | **0.428** | e |
| 12 | Glu-TTCA | **0.143** | k |
| 13 | Asp-TTCA | **0.195** | j |
| 14 | His-TTCA | **0.717** | b |
| 15 | Lys-TTCA | **0.774** | a |

**[M]** *"The browning intensity was significantly higher (p < 0.05) in the systems with added
amino acids than in the TTCA control (Fig. 1c)"* — **but Cys-TTCA (0.072) is BELOW the control
(0.106)**, and its letter `m` is the lowest group. **The sentence is falsified by the authors'
own figure for the Cys system.** Flagged. Mechanistically this is not surprising — added Cys
regenerates H₂S and diverts carbonyls into volatile sulfur compounds instead of melanoidin
(Table 1: Cys-TTCA has the **highest** thiol subtotal after Lys and **zero** oxygen
heterocycles) — but the claim as written is wrong.

**[Z] Browning vs volatile yield across the four systems with both measurements:**
| system | A₄₂₀ | total volatiles (μg L⁻¹) |
|---|---|---|
| TTCA | 0.106 | 43.155 |
| Gly-TTCA | 0.491 | 59.793 |
| Glu-TTCA | 0.143 | 41.717 |
| Lys-TTCA | 0.774 | 105.637 |
**Pearson r = 0.952 on four points** — suggestive of coupled pigment/volatile flux, but n = 4 and
Glu is the point that breaks any tidy story (near-control browning at near-control volatiles).

## 8a. FIGURES 1d and 1e — stacked class profiles **[D], redundant**
Stacked bars for TTCA / TTCA-Gly / TTCA-Glu / TTCA-Lys. Fig. 1d y-axis
**"Concentration（μg/L）" 0–100**, Fig. 1e y-axis **"Types" 0–50**. Segments: Thiols and
sulfides · Thiophenes · Thiazoles · Nitrogen-containing heterocycles · Oxygen-containing
heterocycles. Reading the stack tops: **1d ≈ 43 / 60 / 41 / 106**, **1e ≈ 28 / 32 / 26 / 33** —
which **exactly reproduce Table 1's "Total" and "Total kinds" rows** (43.155 / 59.793 / 41.717 /
105.637 and 28 / 32 / 26 / 33). **Fig. 1d/1e contain no information beyond Table 1.**

## 8b. FIGURE 1f — mechanism scheme, no numbers **[M], qualitative**
A large hand-drawn pathway map, TTCA → 1-DX / 3-DX → (2(5H)-furanone, 1,4-dideoxyosone,
4-hydroxy-5-methyl-3(2H)-furanone, MGO, glycolaldehyde, GO, furfural, 2-methylfuran,
2-deoxyaldotetrose) → (2-furfurylthiol, 2-methyl-3-furanthiol, 2-methylthiophene, thiophene,
furan, 2-acetylfuran, methylpyrazine, 2,5-dimethylpyrazine, pyrazine, pyridine,
2-pyridinecarboxaldehyde, pyrrole, 2-methyl-1H-pyrrole, 1H-pyrrole-2-carboxaldehyde,
5-aminopentanal, 1-iminopropan-2-one). **Every arrow is unlabelled — no rate, no branching
ratio, no yield.** Reagents shown on arrows: H₂O, H₂S, HCHO, NH₃/H₂N–, CO₂, [O], 2H₂O.
Useful only as a topology check against the repo's reaction network.

---

## 9. ANTIOXIDANT SECTION — the numbers that exist **[M]** and the ones that don't

Only **six numeric values** are printed for the whole antioxidant block (Figs 3 and 4, pp.
3248–3249); everything else is described in words with the data left in unlabelled figures.

| quantity | value | condition | anchor |
|---|---|---|---|
| Fe²⁺ chelating ability, TTCA | **71.35 %** | concentration not stated; presumably 3.0 mg mL⁻¹ | p. 3248 |
| Fe²⁺ chelating ability, ARP | **51.03 %** | ditto | p. 3248 |
| Lipid-peroxidation inhibition, MRPs | **2.18 % → 36.39 %** over the concentration sweep | 0.5 → 3.0 mg mL⁻¹ | p. 3248 |
| Lipid-peroxidation inhibition, ARP | **29.18 %** | 3.0 mg mL⁻¹ | p. 3248 |
| Lipid-peroxidation inhibition, TTCA | **31.30 %** | 3.0 mg mL⁻¹ | p. 3248 |
| Equivalence claim | **3.0 mg mL⁻¹ TTCA or ARP ≈ 0.5 mg mL⁻¹ ascorbic acid** for superoxide-anion scavenging; *"approximately **6 g of TTCA** could theoretically replace **1 g of ascorbic acid**"* | | p. 3249 |

Concentration range for all assays: **0.5–3.0 mg mL⁻¹** (genuinely mg, raster-verified).
Ascorbic acid reaches *"nearly 100 % scavenging across the same range"* for DPPH.
**"MRPs"** here means the complete Cys–Xyl Maillard reaction product mixture; **it is never
defined by heating time, temperature or concentration anywhere in the paper.** ⚠ That makes
every MRI-vs-MRP comparison in §3.4 unreproducible.
Fig. 4d is a proposed radical-scavenging mechanism drawing with **no numbers**.

---

## 10. USABILITY CAVEATS THAT APPLY TO EVERY NUMBER ABOVE

1. **Single endpoint (120 min) for every thermal experiment.** No time-resolution, therefore
   no rate constants and no way to separate formation from consumption.
2. **Sealed pressure vessels above 100 °C.** The 120 and 140 °C runs are superheated liquid
   water under autogenous pressure. Volatiles including H₂S cannot escape; the system is
   closed. This is **not** comparable to open-pan or retort geometry and it inflates
   sulfur-volatile retention relative to any vented process.
3. **The pH labels are initial pH only.** Measured final pH: **5.1, 4.9, 4.7** from initial
   5.5, 7.0, 8.0. All three runs end acidic; the pH-8 system spends its later half nearer
   pH 5. Do not model these as constant-pH runs.
4. **Tier-B semi-quantitation** (no asterisk) assumes an IS response factor of 1. Roughly
   **half** of Table 1's compounds — including every thiophene except five, every pyrrole,
   2-acetylthiazole, 2-acetylfuran and 5-methylfurfural — are Tier B.
5. **No LOD/LOQ, no recovery, no calibration range or R²** for any standard curve; Table S3
   is in the SI **[S]**.
6. **No purity figure for TTCA or ARP.** Starting-material quality is unquantified.
7. **SPME extraction time and temperature are not stated** (deferred to ref 15). Headspace
   partitioning coefficients cannot be back-computed. **This makes the μg L⁻¹ values
   "concentration in the analysed liquid as calibrated", not headspace-corrected liquid
   concentrations** — the distinction matters if the repo compares against liquid-phase models.
8. **Fig. 1a's "Sulfur-containing compounds" group is corrupted** (§5a). **Fig. 2c/d/g/h's
   left axis is mislabelled** (§6d). **Table 1's "Kinds" rows are miscounted** in at least
   three columns (§4a, §4c). **The Fig. 1c claim about added amino acids is falsified by its
   own Cys bar** (§8).
9. **"Room temperature" in the pH storage arm is never given a number**, and the temperature
   arm's pH is never given either — so the two solution storage arms cannot be joined (§6c).
10. **n = 3, and the SDs are implausibly small** — several Table 1 cells carry relative SDs
    below 1 % (e.g. 0.012 ± 0.000, 0.025 ± 0). For HS-SPME-GC-MS with n = 3 that is
    optimistic; treat the ± as within-run injection precision, not method precision.

---

## 11. VERDICT — what is usable from this PDF

### USABLE NOW

| block | count | status |
|---|---|---|
| **Table 1 — per-compound volatiles, 48 compounds × 9 systems** | **236 non-ND numeric cells + 24 subtotals/totals** | **FULLY USABLE.** Every cell legible; all 9 grand totals and all 26 block subtotals reconcile arithmetically. Units **μg L⁻¹**. Tier A (\*) vs Tier B must be preserved. |
| **Table 1 retention indices** | 48 × 2 | **USABLE** as DB-Wax RI anchors; two errors flagged (4,5-dimethyl-thiazole RIᶜ, 2-pyridinecarboxaldehyde ordering). |
| **Table 2 — colour after 60 d at 3 a_w** | 24 values + 6 ΔE | **FULLY USABLE and independently verified** (all 6 ΔE reproduce exactly). |
| **Fig. 1c — A₄₂₀ for 15 systems** | 15 | **USABLE [D]**, ± 0.005. |
| **Fig. 2a/b/e/f/i/j — 6 storage time-courses × 3 levels × 7 days** | **126 points** | **USABLE [D]**, ± 1 mmol L⁻¹, with 3 of 6 day-60 endpoints pinned by printed text. |
| **Fig. 1a/1b — class totals at 3 T and 3 pH** | 42 (minus the 6 corrupted) | **USABLE [D] with the Fig. 1a "S-containing" group EXCLUDED.** |
| **18 first-order storage k, 2 Eₐ, 2 a_w slopes** | 22 | **NEW — [Z], never printed by the authors.** §6c. |
| Printed storage-loss percentages, furfural/thiazole/pyrazine ladder values, antioxidant values | 6 + 5 + 6 | **[M], usable.** |

### NOT AVAILABLE FROM THIS FILE — the gap the orchestrator must decide about

1. **Table S4 [S]** — the full per-compound volatile table at **100/120/140 °C**. *This is the
   dataset the wave brief was actually after.* **Not in the PDF.**
2. **Table S5 [S]** — the full per-compound volatile table at **pH 5.5/7/8**. Not in the PDF.
3. **Table S3 [S]** — the standard curves (which compounds, what range, what R²). Not in the PDF.
4. **Table S2 [S]** — the saturated-salt ↔ a_w mapping. Not in the PDF.
5. **Table S1 / Fig. S1 [S]** — TTCA and ARP NMR/MS characterisation and tautomer forms.
6. **Fig. S2 / S3 [S]** — *"the thermal degradation of TTCA and free Cys"*, cited on p. 3242 as
   the evidence for increased H₂S release with temperature. **The H₂S-release evidence is in
   the SI, not here.**
7. **MFT and FFT at 100/120/140 °C** — qualitative statement only (§5c). **No numbers.**
8. **Any rate constant, half-life, Eₐ, or reaction order for the thermal chemistry** — the
   paper contains none, at any temperature.
9. **Any H₂S, NH₃, or α-dicarbonyl measurement.** All three are invoked constantly in the
   mechanism discussion; none is measured in the main text.

**Score against the brief: DOI ✔, temperature ladder present but figure-only ✔/✗,
MFT/FFT ladder ✗ (SI-only), absolute-vs-peak-area question ANSWERED — absolute, μg L⁻¹,
1,2-dichlorobenzene IS, two quantitation tiers.**

---

*Extraction performed 2026-08-28 (Wave K4c). Table 1, Table 2 and all captions re-read off
200 dpi page rasters; every unit re-read off 500–900 dpi crops; all figure values digitised
programmatically from 450–900 dpi renders with the method and residual error stated per
figure. No cell in Table 1 or Table 2 was unreadable.*
