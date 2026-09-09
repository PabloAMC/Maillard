# Zhang, Wang & Cao 2023 — EXTRACTION (DMDS and DMTS during 56-day storage at 50 °C of cooked methionine/cysteine/xylose and methionine/TTCA models; MeSH scavenging by MFT and FFT)
### The 2023 half of the BTBU pair. Same buffer, cook and instrument family as Zhang 2024, but here the time axis is storage at 50 °C, and the whole concentration table is printed in the main text.

**Source on disk:** `data/articles/Zhang2023.pdf` (owner's download, 2026-09-08). Read from the
`pdftotext` text layer in the scratchpad; clean. Table 1 re-typed row by row. Figures 2–5 (Met by
HPLC, verification series, antioxidant assays, correlations) are FIGURE-ONLY. The Supporting Information
(Table S1 labelling design, Table S2 standard curves, Table S3 isotopomer ratios, Fig. S1) is not on
disk. `Zhang2024_extraction.md` §1 (conditions) and §4 (the only rate constants) were read first; the two
systems are compared in §2 below.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Maillard Reaction Products Derived from Cysteine on the Formation of Dimethyl Disulfide and Dimethyl Trisulfide during Storage" |
| Authors | Zeyu Zhang, Bei Wang, Yanping Cao (Beijing Technology and Business University) |
| Venue | J. Agric. Food Chem. 2023, 71, 13043–13053; received 15 April 2023, published 25 August 2023 |
| DOI | 10.1021/acs.jafc.3c02450 |
| Systems | Met/Cys/Xyl ternary; Met/TTCA binary; Met alone; Met/Xyl — all pH 4.9 phosphate, cooked 115 °C / 60 min, then stored sealed at 50 °C for 1–56 days |
| Analytes | 14 volatiles by HS-SPME GC-MS with external standard curves (Table S2); Met by HPLC; ¹³C-Met and d3-Met isotopomer ratios; reducing power and DPPH |
| Companion | Zhang, Kuang, Wang & Cao 2024, Food Res. Int. 182, 114149 (the thiamine/MFT radical paper, dossier on disk) |

## 1. Why it matters

Programme 6 needs methional → MeSH → DMDS/DMTS on the sulfur lane's oxidant pool; B17 §6 says that
pool is what limits the disulfide channel. This paper gives the only printed concentration–time table
on disk for DMDS, DMTS and methional in a methionine pot: five points over 56 days at 50 °C in four
matrices, after a 115 °C cook. Its content for the model: (i) methional keeps forming for weeks at
50 °C (Met/Xyl reaches 49.9 µg/mL = 0.48 mM, 2.9 % of the methionine); (ii) DMDS climbs 30× and DMTS
9× in Met/Xyl over the same time, with air in a sealed vial as the only oxidant; (iii) cysteine-derived
material suppresses DMDS to 2 % of the Met/Xyl level while TTCA (which regenerates cysteine and H2S)
suppresses DMDS by half but produces 7× more DMTS — the trisulfide follows the H2S supply, as in Chin
1994; (iv) isotope labelling puts all of DMDS's and methional's methyl-sulfur on methionine and shows
DMTS in the cysteine pot carrying non-Met sulfur; (v) a verification series shows MFT and FFT consume
MeSH into mixed disulfides faster than MeSH dimerises, MFT faster than FFT (Fig. 3, no numbers). There
are no rate constants in the paper; the slopes in §3 are mine.

## 2. Methods as they matter to a model

- **Ternary model:** Met/Cys/Xyl **7.5:1:1 (w/w)** at **Met 2.5 mg/mL** → Met **16.8 mmol/L** (MW 149.21),
  Cys 0.333 mg/mL = **2.75 mmol/L** (121.16), Xyl 0.333 mg/mL = **2.22 mmol/L** (150.13). Molarities are
  my conversions.
- **Binary model:** Met/TTCA **1:1 (w/w)** → TTCA 2.5 mg/mL = **9.87 mmol/L** (2-threityl-thiazolidine-4-
  carboxylic acid, C8H15NO6S, MW 253.27; synthesised from Xyl + Cys 2:1 molar at pH 7.4, 115 °C, 60 min,
  purified to 98 %, NMR/MS-confirmed).
- **Controls:** Met alone; Met/Xyl (Xyl level not stated — presumably the ternary's 0.333 mg/mL).
- **Buffer:** phosphate, **pH 4.9 ± 0.1**; strength not stated (Zhang 2024 also does not state it).
- **Cook:** **115 °C, 60 min**, then ice bath. **Storage:** sealed screw-cap vials with PTFE/silicone
  septa, **50 °C**, sampled at **1, 14, 28, 42, 56 days**. Vial size and headspace during storage are not
  stated; the analysis vial is 40 mL with 10 mL liquid. **Oxidant: air in the sealed vial.** No metal
  control, no O2 measurement, no headspace flushing.
- **Quantification:** HS-SPME, 50/35 µm DVB/CAR/PDMS; **10 mL** sample + **1 µL IPDS at 1.886 ng/mL in
  methanol** as internal standard (as printed: 1.9 pg absolute; Zhang 2024 prints 0.0943 ng/mL — both
  implausibly small for an IS and best read as a typo of unknown direction); 45 °C 30 min equilibration,
  30 min extraction, desorption 250 °C 5 min. Agilent 7890B/5977A, DB-WAX 60 m × 0.25 mm × 0.25 µm,
  splitless, He 1.2 mL/min, 40 °C 2 min → 70 °C at 4 °C/min (1 min) → 220 °C at 3 °C/min (4 min) → 230 °C
  at 4 °C/min; EI 70 eV, scan 35–400. ID: NIST14 + RI vs C7–C40. **External standard curves (Table S2,
  absent)** for each compound. Triplicate; means ± SD. ⚠ **Methanethiol itself is not in Table 1** although
  a MeSH standard was bought; the paper's MeSH is inferred from its products.
- **Met:** Waters AccQ·Tag HPLC, mg/mL against a Met curve (Fig. 2a, figure-only).
- **Isotopes:** ¹³C-Met (99 atom %) and d3-Met (≥ 98 atom %), mixed 1:1 with unlabelled or used fully
  labelled, same cook and storage; isotopomer ratios [M]:[M+1]:[M+2] from normalised peak areas (Table S3,
  absent; ratios quoted in text carried below).
- **Verification (Fig. 3, figure-only):** (a) MeSH 1 µg/mL (20.8 µM) + FFT or MFT at 0.25, 0.5, 1, 2 µg/mL
  (2.2–17.5 µM), heated 115 °C 60 min (buffer not stated for this arm); (b) MeSH + FFT (or MFT) 1 µg/mL
  each in pH 4.9 phosphate, stored 50 °C, 1, 7, 14, 28, 42, 56 d, three vials per point. Reported: DMDS,
  DMTS, the mixed disulfide (methyl furfuryl disulfide or 2-methyl-3-(methyldisulfanyl)furan) and the
  thiol dimer.
- **Antioxidant assays:** reducing power (ferricyanide, A700, µM ascorbic-acid equivalents/mL; curve
  y = 0.0428x + 0.131); DPPH % scavenging.
- **Conversions used below (ng/mL → nmol/L):** DMDS ÷ 94.20; DMTS ÷ 126.26; methional ÷ 104.15; MFT and
  FFT ÷ 114.17; bis(2-methyl-3-furyl) disulfide ÷ 226.34; 2-methyl-3-(methyldisulfanyl)furan (MMFT) and
  methyl furfuryl disulfide ÷ 160.26; furfural ÷ 96.08; pyrazine ÷ 80.09.
- **Against Zhang 2024:** same lab, buffer pH, cook (115 °C / 60 min), fibre, GC column family and
  external-curve quantification; 2024's Fig. 2 system is Met:VB1:Xyl 1:3:3 at Met 15 mg/mL (6× this
  paper's Met, thiamine present, dose axis, no storage table); 2024's only rate constants (zero-order
  MMFT, 0.0028 / 0.0031, units unstated) are for its own thermal run, not for storage. This paper fits
  no rate constants at all.

## 3. Tables re-typed

### Table 1. "Changes in Concentrations of Volatile Compounds in Four Different Models ... Stored for 56 Days" — ng/mL, triplicate mean ± SD; nd = not detected

**Met-Cys-Xyl model**

| compound | CAS | 1 d | 14 d | 28 d | 42 d | 56 d |
|---|---|---:|---:|---:|---:|---:|
| thiophene | 110-02-1 | nd | nd | nd | 37 ± 2 | 72 ± 2 |
| **dimethyl disulfide** | 624-92-0 | **55 ± 0** | **68 ± 2** | **73 ± 4** | **90 ± 4** | **119 ± 14** |
| 2-methylthiophene | 554-14-3 | nd | 6 ± 0 | 9 ± 0 | 31 ± 0 | 33 ± 2 |
| 2-methyl-3-furanthiol | 28588-74-1 | nd | 118 ± 8 | 125 ± 3 | 202 ± 20 | 123 ± 7 |
| **dimethyl trisulfide** | 3658-80-8 | **nd** | **9 ± 0** | **9 ± 0** | **10 ± 0** | **11 ± 0** |
| 2-furfurylthiol | 98-02-2 | 72 ± 1 | 71 ± 0 | 68 ± 0 | 72 ± 1 | 65 ± 0 |
| **methional** | 3268-49-3 | **nd** | **507 ± 42** | **827 ± 103** | **2714 ± 172** | **9835 ± 64** |
| 2-methyl-3-(methyldisulfanyl)furan | 65505-17-1 | nd | nd | 61 ± 0 | 61 ± 0 | 64 ± 0 |
| 2-thiophenecarboxaldehyde | 98-03-3 | nd | nd | nd | 107 ± 1 | 108 ± 0 |
| methyl furfuryl disulfide | 57500-00-2 | nd | nd | 2 ± 0 | 2 ± 0 | 2 ± 0 |
| bis(2-methyl-3-furyl) disulfide | 28588-75-2 | nd | nd | 100 ± 0 | 101 ± 0 | 100.62 ± 0 |
| pyrazine | 290-37-9 | nd | 345 ± 24 | 311 ± 7 | 291 ± 7 | 257 ± 3 |
| furfural | 98-01-1 | 239 ± 13 | 700 ± 47 | 282 ± 24 | 495 ± 1 | 773 ± 1 |

**Met-TTCA model**

| compound | 1 d | 14 d | 28 d | 42 d | 56 d |
|---|---:|---:|---:|---:|---:|
| thiophene | nd | nd | nd | 11 ± 0 | 11 ± 0 |
| **dimethyl disulfide** | **61 ± 1** | **211 ± 5** | **672 ± 76** | **1918 ± 211** | **2475 ± 85** |
| 2-methylthiophene | nd | 5 ± 0 | 10 ± 1 | 6.91 ± 0 | 6 ± 0 |
| 2-methyl-3-furanthiol | nd | nd | nd | nd | 90 ± 0 |
| **dimethyl trisulfide** | **9 ± 0** | **10 ± 0** | **39 ± 2** | **318 ± 22** | **1143 ± 124** |
| 2-furfurylthiol | 88 ± 1 | nd | nd | nd | nd |
| **methional** | **nd** | **8766 ± 457** | **70173 ± 4447** | **24764 ± 1340** | **16012 ± 234** |
| 2-methyl-3-(methyldisulfanyl)furan | nd | nd | 61 ± 0 (a) | | |
| 2-thiophenecarboxaldehyde | nd | nd | 225 ± 11 | 275 ± 9 | 398 ± 0 |
| methyl furfuryl disulfide | nd | nd | nd | nd | 3 ± 0 |
| pyrazine | nd | 1213 ± 62 | 995 ± 51 | 804 ± 18 | 886 ± 32 |
| 2-methylpyrazine | nd | nd | nd | 140 ± 4 | 148 ± 3 |
| furfural | 292 ± 4 | 265 ± 7 | 268 ± 9 | 234 ± 6 | 351 ± 18 |

(a) This row is printed with only three cells ("nd nd 61 ± 0"); the text (§3.3) says MMFT was found in
the Met-TTCA model "only ... on the 56th day", so the 61 ± 0 most likely belongs under 56 d and the
42/28 d cells are nd; the print is ambiguous and is carried as printed.

**Met model**

| compound | 1 d | 14 d | 28 d | 42 d | 56 d |
|---|---:|---:|---:|---:|---:|
| **dimethyl disulfide** | **84 ± 0** | **58 ± 1** | **118 ± 6** | **165 ± 11** | **148 ± 12** |
| **dimethyl trisulfide** | **11 ± 0** | **9 ± 0** | **17 ± 0** | **17 ± 0** | **18 ± 1** |
| **methional** | **nd** | **nd** | **459 ± 20** | **321 ± 9** | **855 ± 0** |
| furfural | nd | nd | nd | nd | 148 ± 0 |

**Met-Xyl model**

| compound | 1 d | 14 d | 28 d | 42 d | 56 d |
|---|---:|---:|---:|---:|---:|
| **dimethyl disulfide** | **133 ± 7** | **454 ± 42** | **1323 ± 71** | **1948 ± 294** | **4086 ± 342** |
| **dimethyl trisulfide** | **19 ± 1** | **24 ± 2** | **50 ± 2** | **81 ± 5** | **163 ± 6** |
| **methional** | **nd** | **19039 ± 1189** | **27465 ± 1053** | **34731 ± 4425** | **49906 ± 1946** |
| furfural | nd | 3243 ± 150 | 4877 ± 434 | 6718 ± 342 | 11882 ± 570 |

### Derived from Table 1 (my arithmetic; storage at 50 °C after the cook)

Molar levels at day 56: DMDS 1.26 / 26.3 / 1.57 / 43.4 µmol/L (Met-Cys-Xyl / Met-TTCA / Met / Met-Xyl);
DMTS 0.087 / 9.05 / 0.143 / 1.29 µmol/L; methional 94 / 154 / 8.2 / 479 µmol/L.

| quantity | Met-Cys-Xyl | Met-TTCA | Met | Met-Xyl |
|---|---:|---:|---:|---:|
| DMDS, day 1 → day 56 (ng/mL) | 55 → 119 (×2.2) | 61 → 2475 (×41) | 84 → 148 (×1.8) | 133 → 4086 (×31) |
| DMDS mean slope, d1–d56 | 1.2 ng/mL/d = 12 nM/d | 44 ng/mL/d = 0.47 µM/d | 1.2 ng/mL/d | 72 ng/mL/d = 0.76 µM/d |
| DMDS steepest interval | — | d28–42: 89 ng/mL/d | — | d42–56: 153 ng/mL/d |
| DMTS, day 1 → day 56 | nd → 11 | 9 → 1143 (×127) | 11 → 18 | 19 → 163 (×8.6) |
| DMTS steepest interval | — | d42–56: 59 ng/mL/d = 0.47 µM/d | — | d42–56: 5.9 ng/mL/d |
| DMTS/DMDS, molar, day 56 | 0.069 | 0.34 | 0.091 | 0.030 |
| methional, peak | 9835 (d56, still rising) | 70173 (d28), then falls | 855 (d56) | 49906 (d56, rising ~linearly d14–56: 735 ng/mL/d = 7.1 µM/d) |
| methional as % of Met charged (16.8 mM) | 0.56 % | 0.92 % at d56 (4.0 % at d28) | 0.05 % | 2.9 % |
| DMDS as % of Met (2 S per DMDS) | 0.015 % | 0.31 % | 0.019 % | 0.52 % |

**Sulfur bookkeeping in Met-TTCA, d28 → d56:** methional falls by 54,161 ng/mL = 520 µM; DMDS rises by
1803 ng/mL = 19.1 µM (38 µM MeSH-eq); DMTS by 1104 ng/mL = 8.7 µM (17.5 µM CH3S-eq). The two disulfides
absorb ≈ 56 µM, i.e. **11 % of the methional lost**. Whatever removes methional at 50 °C in that pot,
it is mostly not the dimethyl polysulfides (candidates: further Strecker/aldol chemistry, thiazolidine or
thioacetal adducts with the regenerated cysteine, the 2-thiophenecarboxaldehyde that rises 225 → 398).

### Isotope ratios quoted in text (Table S3 absent)

| compound | model | observation | reading |
|---|---|---|---|
| DMDS, DMTS, methional, day 1 | both | [M]:[M+1]:[M+2] = 1:2:1 with ¹³C-Met : Met 1:1 | all methyl-sulfur from Met |
| DMDS, methional, all days | both | 1:2:1 maintained | entirely Met-derived throughout storage |
| DMTS | Met-TTCA | fully Met-derived at all times | |
| DMTS | Met-Cys-Xyl | ≈ 3:4:3 in the 1:1 mix; in fully labelled ¹³C-Met, [M+1]:[M+2] goes 0:10 → ~1:4 from d14 to d56; same with d3-Met | part of DMTS's sulfur comes from Cys (H2S) in the ternary pot |
| furfural, FFT, MFT | both | no labelled molecules before d42; MFT slightly labelled at d56 (ternary) | carbon/sulfur from Xyl and Cys, not Met |
| methyl furfuryl disulfide | both | 1:2:1; labelled in fully labelled pot | FFT + MeSH |
| 2-methyl-3-(methyldisulfanyl)furan, d56 | Met-Cys-Xyl | [M]:[M+1] = 10:90 (¹³C); [M]:[M+3] = 15:85 (d3) | 85–90 % of its methyl from Met; 10–15 % from elsewhere |

### Other numbers in text

Met-TTCA reducing power minimum 12.94 µM AE/mL ≈ 16 % of Met-Cys-Xyl at that time (→ ternary ≈ 81 µM
AE/mL, derived); reducing power equal in the two models for the first 14 d; Met-Xyl's antioxidant indices
rise with storage while the other three fall. Met-TTCA pyrazine ≈ 4× ternary. DMTS at d56: Met-TTCA
"63 times" Met alone (1143/18 = 63.5 ✓) and "7 times" Met-Xyl (1143/163 = 7.0 ✓); DMDS ternary "2 %" of
Met-Xyl (119/4086 = 2.9 %) and Met-TTCA "half" (2475/4086 = 61 %). Fig. 3 (verification): DMDS and
DMTS fall as FFT or MFT rises (thermal arm); during storage DMDS goes up then down, DMTS rises;
mixed-disulfide and dimer levels "far larger" than DMDS/DMTS, especially in the MeSH–MFT model; MeSH–MFT
production rates higher than MeSH–FFT — all figure-only, no values.

## 4. Numbers the repository can use

Registry keys: `methional`, `dimethyl_disulfide`, `dimethyl_trisulfide`, `2_methyl_3_furanthiol`,
`2_furfurylthiol`, `bis_2_methyl_3_furyl_disulfide`, `furfural`, `2_methylthiophene`, `methylpyrazine`,
`pyrazines` (class), `methanethiol` (not measured here). Not in registry: methionine, TTCA,
2-methyl-3-(methyldisulfanyl)furan, methyl furfuryl disulfide, pyrazine (parent), thiophene,
2-thiophenecarboxaldehyde.

| quantity | value | unit | conditions | source | evidence class | registry key |
|---|---|---|---|---|---|---|
| DMDS time series, Met/Xyl | 133, 454, 1323, 1948, 4086 at 1, 14, 28, 42, 56 d | ng/mL | 16.8 mM Met (+ Xyl), pH 4.9 phosphate, cooked 115 °C/60 min, stored 50 °C sealed under air | Table 1 | measured_rate (five-point series; no fit by authors) | dimethyl_disulfide |
| DMTS time series, Met/Xyl | 19, 24, 50, 81, 163 | ng/mL | same | Table 1 | measured_rate | dimethyl_trisulfide |
| methional time series, Met/Xyl | nd, 19039, 27465, 34731, 49906 | ng/mL | same | Table 1 | measured_rate | methional |
| DMDS / DMTS / methional series, Met/TTCA | see Table 1 | ng/mL | 16.8 mM Met + 9.87 mM TTCA, same | Table 1 | measured_rate | as above |
| DMDS / DMTS / methional series, Met/Cys/Xyl | see Table 1 | ng/mL | + 2.75 mM Cys, 2.22 mM Xyl | Table 1 | measured_rate | as above |
| DMDS / DMTS / methional series, Met alone | see Table 1 | ng/mL | 16.8 mM Met, same | Table 1 | measured_rate | as above |
| DMDS mean formation rate, Met/Xyl, 50 °C | 0.76 (mean d1–56); 1.6 (d42–56) | µmol/L/day | air-sealed vial, pH 4.9 | derived | within_study_ratio | dimethyl_disulfide |
| DMTS/DMDS molar ratio at 56 d | 0.030 (Met/Xyl); 0.069 (ternary); 0.091 (Met); 0.34 (Met/TTCA) | — | same | derived | within_study_ratio | both |
| DMDS suppression by cysteine MRPs | ternary 2.9 % of Met/Xyl; TTCA 61 % of Met/Xyl at 56 d | ratio | same | Table 1 | within_study_ratio | dimethyl_disulfide |
| DMTS enhancement by TTCA | 7.0 × Met/Xyl; 63 × Met at 56 d | ratio | same | Table 1, text | within_study_ratio | dimethyl_trisulfide |
| MFT dimer share, ternary pot at 28–56 d | bis-MFT 100 ng/mL (0.44 µM = 0.88 µM MFT-eq) vs free MFT 125–202 ng/mL (1.1–1.8 µM): dimer holds 33–45 % of MFT-eq | ratio | 50 °C storage, air | Table 1 | within_study_ratio | bis_2_methyl_3_furyl_disulfide |
| MMFT / bis-MFT in ternary pot | 61–64 vs 100–101 ng/mL → 0.38 / 0.44 = 0.86–0.90 nmol/nmol | ratio | same | Table 1 | within_study_ratio | not in registry / bis_2_methyl_3_furyl_disulfide |
| FFT decay, ternary | 72 → 65 ng/mL over 56 d (−10 %) | ng/mL | same | Table 1 | level_only | 2_furfurylthiol |
| FFT decay, Met/TTCA | 88 → nd by day 14 | ng/mL | same | Table 1 | level_only | 2_furfurylthiol |
| isotope: DMDS, methional methyl-S | 100 % from Met | — | both models | text | within_study_ratio | dimethyl_disulfide, methional |
| isotope: DMTS in ternary pot | part of the sulfur from Cys; fully-labelled [M+1]:[M+2] → ~1:4 by d56 | — | ternary | text | within_study_ratio | dimethyl_trisulfide |
| MFT vs FFT as MeSH scavengers | MFT faster (direction only) | qualitative | 115 °C / 60 min and 50 °C storage | Fig. 3 | figure_only | 2_methyl_3_furanthiol, 2_furfurylthiol |

## 5. Flags

1. **Methanethiol was not measured.** Every statement about MeSH is inferred from DMDS, DMTS and the
   mixed disulfides. The lane's MeSH balance cannot be closed with this paper.
2. **Methional is "nd" at day 1 in all four models** after a 115 °C / 60 min cook of 16.8 mM methionine
   with xylose, yet reaches 19 µg/mL by day 14 at 50 °C. Either the cook makes almost no methional at
   this pH 4.9 (possible — Zhang 2024's own cook reports MeSH-derived products but not methional) or the
   day-1 measurement is off. Treat day 1 as suspect for methional and use the d14–d56 slope.
3. **Several entries look like a calibration ceiling or floor**, not measurements: bis-MFT 100 ± 0,
   101 ± 0, 100.62 ± 0; MMFT 61 ± 0 at every time it appears in both models; 2-thiophenecarboxaldehyde
   107 ± 1 / 108 ± 0; DMTS 9 ± 0, 9 ± 0, 10 ± 0, 11 ± 0. "± 0" SDs on triplicates at four significant
   figures are not credible. The DMDS, DMTS (Met/TTCA and Met/Xyl) and methional rows carry real SDs and
   are the rows to use.
4. **Storage-vial headspace is not stated**, so the oxygen inventory per vial is unknown; the 10 mL in
   40 mL geometry is the analysis vial. Since B17's question is exactly how much oxidant a sealed pot
   holds, this is the missing number. If the storage vial were the 40 mL analysis vial with 10 mL liquid, its 30 mL
   air headspace holds ≈ 0.24 mmol O2 against 0.43 µmol DMDS made in the 10 mL (43.4 µM) — about
   500-fold; so oxidant supply was NOT limiting in Met/Xyl here, and the slow, accelerating DMDS rise
   reflects MeSH supply from methional, not oxidant exhaustion. (My estimate on an assumed geometry.)
5. **Methional's fate is not the disulfides** (§3 bookkeeping: 11 % in Met/TTCA). A model that routes
   all methional loss to MeSH → DMDS/DMTS will overpredict the polysulfides by ~10× in a cysteine-bearing
   pot.
6. **The isotope work is in Table S3 (absent)**; only the ratios quoted in text are carried.
7. **External standard curves (Table S2) absent**; HS-SPME with a single IS at a stated 1.9 pg — the
   absolute ng/mL are not commensurable with SIDA numbers elsewhere in the repo (same caution as the
   2024 dossier's §1.4). Ratios within a model and across models on the same day are the safer use.
8. **Two different Cys-derived pots do opposite things to DMTS**: free Cys + Xyl suppresses it (11 ng/mL)
   while TTCA drives it (1143 ng/mL). The paper's own explanation is antioxidant capacity plus MFT
   scavenging in the ternary pot; the isotope data add that TTCA's regenerated Cys supplies H2S, which
   Chin 1994 shows is required for DMTS. Any DMTS step in the model should take H2S as a co-reactant.
9. Buffer strength, storage vial volume, Xyl level in the Met/Xyl control, and the verification arm's
   buffer are all unstated.
