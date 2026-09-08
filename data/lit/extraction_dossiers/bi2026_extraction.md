# Bi et al. 2026 — EXTRACTION (pea milk grinding time course 0 s - 8 min: lipid substrate depletion by UPLC-MS/MS lipidomics; hexanal + hexanol by HS-SPME GC-MS SIM)
### The first time-resolved LOX-pathway study of pea milk on disk — but the product side (hexanal + hexanol) is printed only as a figure; what the tables carry is the SUBSTRATE side (free linoleic acid, phospholipids, triglycerides) at eight quench times.

**Source on disk:** `data/articles/Bi2026.pdf` (owner's download, 2026-09-08; Food Control 185 (2026)
112063, 10 pages, no supplementary material appended). Read from the scratchpad text layer
(`Bi2026.txt`); Table 1 (lipid subclasses x 8 times) and Table 2 (42 lipid species, 0 s vs 8 min)
are clean and re-typed below; page 3 was re-extracted with pypdf to confirm that Fig. 1 (the
hexanal + hexanol time course) leaves no numbers in the text layer — it is FIGURE-ONLY. Tables
A.1-A.4 and Figs A.1-A.2 are online-only supplementary and are NOT on disk. Repo status before this
dossier: backlog item LOX-01 ("the beany note before heat") holds only end-of-process levels
(Zhang 2020b: raw pea milk hexanal 164 ug/L, free linoleic acid ~14 300 ug/L; LOX-2 2160 U/mg) and
one isolate LOX residue (Gao 2020); no time course of anything.

## 0. Identity

| field | value |
|---|---|
| Title | "Quality-driven grinding control: Monitoring lipid oxidation markers for off-flavor mitigation in pea milk" |
| Authors | Shuang Bi, Lichang Yan, Xiaoying Xiong, Fan Yang, Ye Liu (corresponding) — Flavor Science Laboratory, Beijing Technology and Business University |
| Venue | Food Control 185 (2026) 112063; received 3 Nov 2025, accepted 16 Feb 2026, online 17 Feb 2026 |
| DOI | 10.1016/j.foodcont.2026.112063 (PII S0956-7135(26)00108-8) |
| Pea | Zhongwan No. 6 (ZW.6), Pisum sativum L., seeds stored at 4 C |
| Companions from the same group | Yan et al. 2024, Food Chem 445, 138696 (pea milk key odorants, the calibration method cited here); Bi et al. 2022, Food Chem 380, 132203 (pea seed milk flavour "via enzyme activity inhibition") — neither on disk; both are candidate acquisitions for LOX-01 (flag 9) |
| Data | "Data will be made available on request." |

## 1. Why it matters

The repository has no rate for the lipoxygenase step that makes hexanal / 1-hexanol before any heat
is applied. This paper is the first on disk that quenches pea milk at eight grinding times (0 s, 10
s, 20 s, 40 s, 60 s, 2, 4, 8 min) and measures something at each: a full lipidome (693 species,
UPLC-ESI-QTRAP-MS/MS, MRM, internal-standard quantification) and the sum hexanal + hexanol
(HS-SPME GC-MS SIM, five-point calibration). The product time course is only drawn (Fig. 1), so it
cannot be transcribed; the text does say what shape it has (fast rise in the first 2 min, slower
rise to 8 min, plateau after 8 min). The substrate time course IS printed (Table 1): total free
fatty acid falls 34.55 -> 13.53 ug/g, PE 82.96 -> 37.67 ug/g, TG 175.49 -> 69.66 ug/g, and Table 2
gives free linoleic acid 10.658 -> 0.309 ug/g between 0 s and 8 min. For a LOX module this is the
first measured depletion of the substrate pool against time, in the same 12.5 % (w/w) seed slurry
as Zhang 2020b, and it comes with the authors' own kinetic statements and lipid-vs-C6 correlation
coefficients (-0.83 to -0.99).

## 2. Methods as they matter to a model

- **Seeds and soak:** ~50 g pea seeds washed, soaked in 150 g distilled water at **4 C for 12 h**
  (seed:soak water 1:3 w/w). No blanching, no hot water anywhere.
- **Grinding:** "ground at a ratio of **1:7 (w/w) peas to distilled water** using a soymilk machine
  (P165, Joyoung)". Whether "peas" means dry seed or soaked seed is not stated; on a dry-seed
  basis 1:7 is a **12.5 % w/w slurry**, the same nominal strength as Zhang 2020b's 12.5 % seed
  slurry. **Grinding temperature is NOT stated**; water temperature is NOT stated (distilled water,
  presumably ambient; the seeds come out of a 4 C soak). ⚠ A household soymilk machine can heat
  during its cycle; the paper never says heating was disabled (flag 1). **pH during the reactive
  window is NOT stated** (natural pea slurry, no buffer).
- **Time design (the essential point):** "Pea milk was prepared at the initial grinding time of 0
  s, 10 s, 20 s, 40 s, 60 s, 2 min, 4 min, 8 min, and 12 min, respectively, following which citric
  acid was added to a final content of **5 %** to terminate the subsequent enzymatic reaction
  (Fischer et al., 2020)." After the quench "all samples were ground till 12 min" for texture.
  "Treat-X denotes the duration of active enzymatic reaction prior to acid quenching, despite
  equal total grinding time." For Treat-0 s the citric acid was added **prior to grinding**
  (§3.2). So each sample is: X min of enzyme time in unbuffered slurry, then (12 - X) min of
  grinding at ~5 % citric acid. The 12-min point exists in the preliminary experiment but is not
  in the tables (eight groups analysed: 0 s - 8 min).
- **Why 8 min:** "The results of the preliminary experiment showed that when the grinding time
  reached 8 min, the total content of hexanal and hexanol did not increase, indicating that 8 min
  was sufficient for enzymatic reactions."
- **Work-up:** slurry filtered through two layers of gauze, centrifuged **3000 rpm, 15 min, 4 C**;
  the **supernatant is the "pea milk"** that is analysed (volatiles and lipids). Snap-frozen in
  liquid N2, stored at -80 C. ⚠ Lipid values are therefore ug per g of centrifuged supernatant,
  and lipid partition between pellet and supernatant can move with grinding (the authors invoke
  "lipid release during grinding" to explain the 20 -> 40 s rise in Table 1).
- **Volatile isolation:** 5 mL pea milk + 0.5 g NaCl + **1 uL 2-methyl-3-heptanone internal
  standard** in a 20 mL headspace vial; equilibrated **50 C / 20 min**; DVB/CAR/PDMS 50/30 um
  fibre, **40 min** extraction, 450 rpm shaking; desorbed 250 C / 5 min.
- **GC-MS:** Agilent 7890B-5977A, DB-WAX 60 m x 0.25 mm x 0.25 um; 40 C (2 min) -> 75 C at 7 C/min
  -> 230 C at 5 C/min (2 min); splitless; He 1.8 mL/min; EI 70 eV; **SIM**, ions m/z 82, 56, 44
  for hexanal and 84, 56, 69 for hexanol; source 230 C, quadrupole 150 C.
- **Quantification of hexanal and hexanol, verbatim:** "The volatile compounds were quantified
  using a five-point standard calibration curve (Yan et al., 2024). A deodorized pea milk model
  solution was established using 2 % fat cow's milk as the matrix and was diluted to a protein
  content similar to that in pea milk (Zhang et al., 2020). The linear regression standard curves
  were plotted using the mass spectrum peak area ratio of the analyte and internal standard (1 uL
  of 2-methyl-3-heptanone at 27.2 ug/L) as the ordinate and the content ratio of the analyte and
  internal standard as the abscissa." §3.1 adds: "quantified using a standard calibration curve
  (Table A.1) in conjunction with selected ion monitoring (SIM) mode, yielding reliable results
  (R2 > 0.99)." So: **matrix-matched (cow's milk) internal-standard calibration, SIM, R2 > 0.99;
  the curve parameters are in Table A.1, which is not on disk.** The concentration unit of Fig. 1
  is not recoverable from the text layer (the same lab reports ug/L in Zhang 2020b-style work,
  but that is an inference, not a reading).
- **Lipid extraction:** 200 uL pea milk + 1 mL MTBE:MeOH 3:1 containing an internal-standard
  mixture; + 100 uL water; shake 2500 r/min 15 min; centrifuge 12 000 r/min 3 min 4 C; 500 uL
  upper phase evaporated at 20 C, reconstituted in 200 uL mobile phase B; 120 uL injected.
- **UPLC-ESI-QTRAP-MS/MS:** ExionLC + QTRAP 6500+; Thermo Accucore C30 2.6 um 2.1 x 100 mm; A =
  ACN/water 60/40, B = ACN/IPA 10/90, both 0.1 % formic acid + 10 mM ammonium formate; gradient
  80:20 (0 min) -> 70:30 (2) -> 40:60 (4) -> 15:85 (9) -> 10:90 (14) -> 5:95 (15.5-17.3) -> 80:20
  (17.5-20 min); 0.35 mL/min, 45 C, 2 uL. ESI 500 C, +5500 V / -4500 V; MRM, N2 collision gas 5
  psi. Identification against the MWDB database (Metware, Wuhan) by RT + ion pair;
  "quantitative analysis was performed by the internal standard method" (which standards, and
  whether one per class, is not stated). QC every 10 samples; > 75 % of species with CV < 0.3.
- **Replication / statistics:** all experiments in triplicate, mean +/- SD; ANOVA + Duncan (p <=
  0.05), letters printed in Tables 1-2; PCA on the eight groups (PC1 55.21 %, PC2 11.02 %);
  differential lipids by VIP > 1 and |log2 FC| >= 1.
- **NOT measured:** LOX activity (no assay of any kind); peroxide value; conjugated dienes /
  hydroperoxides; hexanal and hexanol as separate printed numbers; pH; temperature. The only
  "kinetic" quantities the paper prints are the substrate concentrations of Tables 1-2 and the
  correlation coefficients of §3.7.

## 3. Tables re-typed

### Fig. 1 — "Changes of total contents of hexanal and hexanol in pea milk under different enzymatic reaction time" — FIGURE-ONLY

No value, unit or axis range survives in the text layer. What the text states about it:
- "the sum of hexanal and hexanol increased significantly with prolonged enzymatic reaction time.
  Notably, no further increase in hexanal and hexanol was observed beyond 8 min of treatment."
- "the combined content of hexanal and hexanol increased rapidly within the first 2 min of
  grinding, followed by a slower yet continued increase up to 8 min."
- Contrast the authors draw: in soybeans "hexanal and hexanol formation plateaus after just 1 min
  (Trindler et al., 2022)"; attributed to "differences in LOX isoenzyme activity, substrate
  accessibility, or endogenous antioxidant systems".
- The initial (Treat-0 s) value, the plateau value and the separate hexanal and hexanol series are
  not given anywhere in the text. Hexanal and hexanol separately enter only Supplementary Fig. A2
  (correlations), not on disk.

### Table 1. "Changes in the content of each lipid subclass during the enzymatic reaction of pea milk."

Unit: **ug/g** (of centrifuged pea milk supernatant), mean +/- SD, n = 3; superscript letters =
Duncan groups within a row (p <= 0.05). Columns: Treat-0 s, 10 s, 20 s, 40 s, 60 s, 2 min, 4 min,
8 min. Footnote a spells the abbreviations (ADGGA acyl diacylglyceryl glucuronide; Cer ceramide;
Cert phytoceramide; CoQ coenzyme Q; DG diacylglycerol; DGDG digalactosyldiacylglycerol; DGGA
diacylglyceryl glucuronide; DGTS diacylglyceryl trimethylhomoserine; FFA free fatty acid; HexCer
hexosylceramide; LDGTS lyso-DGTS; LPA/LPC/LPE/LPG/LPI lyso-phosphatidic acid / -choline /
-ethanolamine / -glycerol / -inositol; MG monoglyceride; MGDG monogalactosyldiacylglycerol; PA
phosphatidic acid; PC phosphatidylcholine; PE phosphatidylethanolamine; PG phosphatidylglycerol;
PI phosphatidylinositol; PMeOH phosphatidylmethanol; PS phosphatidylserine; SPH sphingosine; SQDG
sulfoquinovosyldiacylglycerol; TG triglyceride).

| subclass | 0 s | 10 s | 20 s | 40 s | 60 s | 2 min | 4 min | 8 min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Total** | 525.12 +/- 24.88 a | 496.05 +/- 57.94 a | 494.89 +/- 29.10 a | 535.80 +/- 15.39 a | 500.06 +/- 45.95 a | 363.86 +/- 45.25 b | 293.43 +/- 42.36 bc | 271.67 +/- 60.35 c |
| ADGGA | 0.19 +/- 0.02 ab | 0.14 +/- 0.01 c | 0.17 +/- 0.04 abc | 0.20 +/- 0.01 a | 0.16 +/- 0.02 bc | 0.11 +/- 0.01 d | 0.07 +/- 0.02 de | 0.05 +/- 0.01 e |
| Cer | 0.06 +/- 0.01 ab | 0.06 +/- 0.01 ab | 0.06 +/- 0.01 ab | 0.07 +/- 0.00 ab | 0.07 +/- 0.01 a | 0.06 +/- 0.01 b | 0.05 +/- 0.01 b | 0.07 +/- 0.01 ab |
| Cert | 0.45 +/- 0.03 abc | 0.42 +/- 0.08 bc | 0.47 +/- 0.02 abc | 0.55 +/- 0.03 a | 0.52 +/- 0.04 ab | 0.42 +/- 0.04 bc | 0.38 +/- 0.03 c | 0.44 +/- 0.11 bc |
| CoQ | 0.50 +/- 0.04 a | 0.39 +/- 0.02 b | 0.40 +/- 0.04 b | 0.53 +/- 0.02 a | 0.48 +/- 0.03 a | 0.39 +/- 0.03 b | 0.27 +/- 0.05 c | 0.33 +/- 0.08 bc |
| DG | 7.06 +/- 0.59 ab | 6.34 +/- 0.50 abc | 5.97 +/- 0.60 bcd | 7.24 +/- 0.12 a | 6.46 +/- 0.59 abc | 5.58 +/- 0.70 cde | 4.79 +/- 0.50 e | 4.90 +/- 0.89 de |
| DGDG | 0.73 +/- 0.04 a | 0.67 +/- 0.04 ab | 0.62 +/- 0.04 b | 0.70 +/- 0.04 a | 0.67 +/- 0.05 ab | 0.50 +/- 0.03 c | 0.37 +/- 0.04 d | 0.32 +/- 0.03 d |
| DGGA | 0.28 +/- 0.06 ab | 0.26 +/- 0.01 bc | 0.26 +/- 0.07 bc | 0.33 +/- 0.02 a | 0.26 +/- 0.02 abc | 0.20 +/- 0.01 cd | 0.15 +/- 0.04 d | 0.13 +/- 0.03 d |
| DGTS | 1.88 +/- 0.04 a | 1.71 +/- 0.05 ab | 1.59 +/- 0.12 bc | 1.86 +/- 0.05 a | 1.79 +/- 0.13 a | 1.50 +/- 0.13 c | 1.24 +/- 0.04 d | 1.25 +/- 0.12 d |
| **FFA** | **34.55 +/- 3.71 a** | **34.12 +/- 6.30 a** | **33.46 +/- 3.52 a** | **29.63 +/- 5.55 ab** | **27.95 +/- 3.77 ab** | **23.15 +/- 4.59 bc** | **18.44 +/- 3.58 cd** | **13.53 +/- 3.98 d** |
| HexCer | 0.38 +/- 0.11 a | 0.30 +/- 0.04 a | 0.38 +/- 0.20 a | 0.33 +/- 0.05 a | 0.32 +/- 0.02 a | 0.32 +/- 0.06 a | 0.31 +/- 0.08 a | 0.32 +/- 0.06 a |
| LDGTS | 0.04 +/- 0.00 a | 0.01 +/- 0.00 b | 0.01 +/- 0.00 b | 0.02 +/- 0.00 b | 0.02 +/- 0.00 b | 0.02 +/- 0.00 b | 0.01 +/- 0.00 b | 0.02 +/- 0.01 b |
| LPA | 0.46 +/- 0.03 b | 0.48 +/- 0.03 b | 0.44 +/- 0.02 b | 0.42 +/- 0.01 b | 0.55 +/- 0.03 a | 0.35 +/- 0.03 c | 0.26 +/- 0.03 d | 0.26 +/- 0.06 d |
| LPC | 0.58 +/- 0.04 a | 0.66 +/- 0.05 a | 0.62 +/- 0.07 a | 0.47 +/- 0.01 b | 0.43 +/- 0.02 b | 0.41 +/- 0.08 b | 0.30 +/- 0.05 c | 0.41 +/- 0.09 b |
| LPE | 6.17 +/- 0.32 abc | 7.32 +/- 0.12 a | 6.45 +/- 0.87 abc | 4.76 +/- 0.19 cd | 5.36 +/- 0.28 bcd | 4.77 +/- 0.52 cd | 4.14 +/- 0.62 d | 6.67 +/- 2.31 ab |
| LPG | 0.09 +/- 0.01 ab | 0.10 +/- 0.01 a | 0.09 +/- 0.00 a | 0.08 +/- 0.00 ab | 0.09 +/- 0.00 ab | 0.06 +/- 0.01 cd | 0.05 +/- 0.01 d | 0.07 +/- 0.02 bc |
| LPI | 2.21 +/- 0.15 ab | 2.51 +/- 0.26 a | 2.34 +/- 0.03 ab | 2.18 +/- 0.08 ab | 2.35 +/- 0.14 ab | 1.58 +/- 0.36 c | 1.29 +/- 0.37 c | 1.80 +/- 0.61 bc |
| MG | 5.75 +/- 0.22 a | 5.85 +/- 0.50 a | 5.70 +/- 0.05 a | 5.85 +/- 0.47 a | 5.70 +/- 0.18 a | 6.00 +/- 0.19 a | 5.47 +/- 0.61 a | 5.51 +/- 0.30 a |
| MGDG | 0.12 +/- 0.01 a | 0.11 +/- 0.01 ab | 0.10 +/- 0.02 ab | 0.12 +/- 0.00 a | 0.12 +/- 0.02 a | 0.09 +/- 0.01 b | 0.07 +/- 0.01 c | 0.06 +/- 0.01 c |
| PA | 0.14 +/- 0.01 a | 0.14 +/- 0.01 ab | 0.13 +/- 0.00 a | 0.14 +/- 0.01 a | 0.16 +/- 0.00 a | 0.11 +/- 0.01 b | 0.09 +/- 0.01 c | 0.08 +/- 0.03 c |
| **PC** | 94.61 +/- 4.48 ab | 87.93 +/- 6.61 b | 84.72 +/- 3.78 b | 105.53 +/- 2.06 a | 95.47 +/- 2.92 ab | 70.50 +/- 7.86 c | 62.60 +/- 8.40 c | 61.00 +/- 15.91 c |
| **PE** | **82.96 +/- 11.68 a** | **79.12 +/- 14.05 a** | **71.40 +/- 8.67 a** | **86.72 +/- 1.12 a** | **82.72 +/- 5.25 a** | **55.75 +/- 5.94 b** | **44.29 +/- 5.86 bc** | **37.67 +/- 8.73 c** |
| PG | 14.22 +/- 0.86 a | 11.33 +/- 1.26 b | 11.42 +/- 0.76 b | 14.55 +/- 1.02 a | 13.89 +/- 0.15 a | 9.38 +/- 1.45 bc | 8.13 +/- 1.22 c | 8.18 +/- 2.02 c |
| PI | 51.14 +/- 3.08 a | 48.07 +/- 7.77 a | 46.88 +/- 0.77 a | 56.09 +/- 3.99 a | 49.52 +/- 8.21 a | 33.57 +/- 3.40 b | 27.94 +/- 3.61 b | 31.57 +/- 8.11 b |
| PMeOH | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a | 0.01 +/- 0.00 a |
| PS | 44.93 +/- 6.51 a | 48.87 +/- 4.25 a | 47.55 +/- 20.10 a | 50.67 +/- 7.82 a | 42.22 +/- 5.37 ab | 34.05 +/- 4.65 abc | 23.77 +/- 2.33 c | 27.28 +/- 5.48 bc |
| SPH | 0.02 +/- 0.00 a | 0.03 +/- 0.00 a | 0.02 +/- 0.00 a | 0.03 +/- 0.01 a | 0.03 +/- 0.01 a | 0.03 +/- 0.00 a | 0.03 +/- 0.01 a | 0.03 +/- 0.01 a |
| SQDG | 0.09 +/- 0.01 a | 0.08 +/- 0.01 ab | 0.08 +/- 0.01 ab | 0.09 +/- 0.00 a | 0.08 +/- 0.01 a | 0.06 +/- 0.00 bc | 0.06 +/- 0.01 c | 0.05 +/- 0.01 c |
| **TG** | 175.49 +/- 25.30 a | 159.01 +/- 26.47 a | 173.54 +/- 15.50 a | 166.62 +/- 17.92 a | 162.65 +/- 46.40 a | 114.90 +/- 20.96 b | 88.87 +/- 18.45 b | 69.66 +/- 13.94 b |

Shape of the table: nothing is significantly different from 0 s until 60 s (letters "a" across
0-60 s for total, FFA, PE, PI, TG, PS); the first significant drop is between 60 s and 2 min for
every major class; FFA alone declines monotonically from 20 s on. The authors' reading (§3.3.1):
0-20 s decrease; 20-40 s increase "lipid release during grinding temporarily exceeded the rate of
enzymatic consumption"; "a pronounced reduction in total lipids occurred between 1 and 2 min",
attributed to hydroperoxides formed in the first minute activating LOX.

### Table 2. "The important off-flavor precursor lipid molecules of pea milk screened on the basis of VIP > 1, |log2 FC| >= 1, and total lipid consumption > 1 ug/g."

Columns: no.; compound; formula; MW; Treat-0 s (ug/g); Treat-8 min (ug/g); lipid consumption
(ug/g) = difference; VIP; fold change; log2 FC; type (all "down"). Letters a/b mark 0 s vs 8 min
significance (all a vs b). Formulae omitted here except where useful; MW as printed.

| no. | compound | MW | 0 s (ug/g) | 8 min (ug/g) | consumed (ug/g) | VIP | FC | log2 FC |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | **PE (18:2_18:2)** | 739.515 | 29.088 +/- 4.569 | 9.145 +/- 2.220 | 19.942 | 1.188 | 0.329 | -1.605 |
| 2 | PS (18:0_17:1) | 775.536 | 14.671 +/- 2.264 | 4.019 +/- 0.965 | 10.652 | 1.240 | 0.269 | -1.896 |
| 3 | **FFA (18:2), free linoleic acid** (C18H32O2) | 280.240 | **10.658 +/- 0.809** | **0.309 +/- 0.056** | **10.349** | 1.266 | 0.034 | -4.873 |
| 4 | PC (18:2_18:2) | 781.562 | 17.728 +/- 1.106 | 7.594 +/- 1.399 | 10.134 | 1.217 | 0.453 | -1.144 |
| 5 | PE (18:2_16:0) | 715.515 | 16.123 +/- 1.906 | 6.177 +/- 1.161 | 9.946 | 1.194 | 0.403 | -1.312 |
| 6 | TG (18:2_18:2_18:3) | 876.721 | 14.457 +/- 2.756 | 4.559 +/- 0.805 | 9.898 | 1.217 | 0.327 | -1.613 |
| 7 | TG (16:0_18:2_18:2) | 854.736 | 16.266 +/- 2.401 | 6.428 +/- 0.859 | 9.838 | 1.225 | 0.406 | -1.300 |
| 8 | TG (18:1_18:2_18:3) | 878.736 | 13.612 +/- 2.147 | 5.431 +/- 0.822 | 8.181 | 1.217 | 0.411 | -1.282 |
| 9 | PE (18:1_18:2) | 741.531 | 13.224 +/- 2.274 | 5.235 +/- 1.040 | 7.989 | 1.175 | 0.396 | -1.338 |
| 10 | **FFA (18:1), free oleic acid** (C18H34O2) | 282.256 | 9.982 +/- 0.474 | 2.944 +/- 0.474 | 7.039 | 1.242 | 0.335 | -1.576 |
| 11 | TG (18:0_18:2_18:3) | 880.752 | 9.917 +/- 1.590 | 3.576 +/- 0.580 | 6.341 | 1.220 | 0.377 | -1.406 |
| 12 | PS (18:0_18:2) | 787.536 | 7.153 +/- 0.809 | 1.546 +/- 0.057 | 5.607 | 1.251 | 0.246 | -2.024 |
| 13 | TG (18:0_18:1_18:2) | 884.783 | 7.025 +/- 0.637 | 2.300 +/- 0.342 | 4.725 | 1.247 | 0.350 | -1.515 |
| 14 | TG (18:1_18:1_18:4) | 878.736 | 6.450 +/- 1.099 | 2.223 +/- 0.373 | 4.227 | 1.211 | 0.372 | -1.427 |
| 15 | TG (16:0_18:2_18:3) | 852.721 | 5.611 +/- 0.576 | 1.806 +/- 0.257 | 3.805 | 1.248 | 0.340 | -1.558 |
| 16 | PG (18:2_16:0) | 746.510 | 5.231 +/- 0.584 | 1.714 +/- 0.340 | 3.516 | 1.227 | 0.348 | -1.522 |
| 17 | TG (16:0_18:1_18:2) | 856.752 | 5.046 +/- 0.862 | 1.534 +/- 0.285 | 3.512 | 1.214 | 0.332 | -1.590 |
| 18 | TG (16:0_18:1_18:3) | 854.736 | 5.556 +/- 0.868 | 2.156 +/- 0.370 | 3.400 | 1.216 | 0.403 | -1.312 |
| 19 | PI (18:1_18:2) | 860.541 | 5.496 +/- 0.547 | 2.768 +/- 0.918 | 2.728 | 1.158 | 0.495 | -1.014 |
| 20 | TG (18:1_18:1_18:3) | 880.752 | 3.721 +/- 0.475 | 1.084 +/- 0.133 | 2.637 | 1.243 | 0.314 | -1.671 |
| 21 | PI (18:0_18:2) | 862.557 | 4.623 +/- 0.307 | 2.007 +/- 0.534 | 2.616 | 1.200 | 0.456 | -1.134 |
| 22 | TG (18:1_18:1_18:2) | 882.768 | 3.937 +/- 0.319 | 1.340 +/- 0.180 | 2.598 | 1.244 | 0.371 | -1.432 |
| 23 | TG (18:0_18:0_18:3) | 884.783 | 3.902 +/- 0.612 | 1.342 +/- 0.204 | 2.560 | 1.225 | 0.362 | -1.465 |
| 24 | TG (18:1_18:2_18:2) | 880.752 | 3.889 +/- 0.561 | 1.374 +/- 0.228 | 2.514 | 1.236 | 0.357 | -1.487 |
| 25 | TG (18:2_18:2_18:2) | 878.736 | 3.301 +/- 0.583 | 1.017 +/- 0.157 | 2.284 | 1.228 | 0.316 | -1.662 |
| 26 | PC (18:0_18:2) | 785.593 | 3.698 +/- 0.134 | 1.605 +/- 0.172 | 2.092 | 1.191 | 0.496 | -1.012 |
| 27 | TG (16:0_16:0_18:2) | 830.736 | 3.111 +/- 0.449 | 1.176 +/- 0.186 | 1.935 | 1.215 | 0.407 | -1.295 |
| 28 | **LPE (18:2)** | 477.286 | 2.018 +/- 0.189 | 0.113 +/- 0.030 | 1.905 | 1.261 | 0.068 | -3.872 |
| 29 | TG (18:2_18:3_18:3) | 874.705 | 2.493 +/- 0.425 | 0.592 +/- 0.088 | 1.901 | 1.237 | 0.254 | -1.975 |
| 30 | PS (18:2_18:1) | 785.521 | 2.770 +/- 0.643 | 0.927 +/- 0.156 | 1.843 | 1.185 | 0.352 | -1.506 |
| 31 | PS (18:2_16:0) | 759.505 | 3.137 +/- 0.479 | 1.303 +/- 0.270 | 1.834 | 1.197 | 0.395 | -1.341 |
| 32 | PC (18:2_18:3) | 779.547 | 2.749 +/- 0.133 | 0.929 +/- 0.159 | 1.820 | 1.232 | 0.378 | -1.405 |
| 33 | TG (16:0_18:0_18:2) | 858.768 | 2.665 +/- 0.466 | 0.856 +/- 0.146 | 1.809 | 1.225 | 0.330 | -1.599 |
| 34 | TG (16:0_16:2_18:2) | 826.705 | 2.183 +/- 0.211 | 0.563 +/- 0.077 | 1.620 | 1.249 | 0.287 | -1.799 |
| 35 | TG (18:0_18:3_20:0) | 912.815 | 2.089 +/- 0.377 | 0.557 +/- 0.065 | 1.532 | 1.221 | 0.301 | -1.734 |
| 36 | TG (18:0_18:1_18:3) | 882.768 | 2.083 +/- 0.256 | 0.594 +/- 0.067 | 1.489 | 1.248 | 0.303 | -1.725 |
| 37 | TG (18:1_18:3_20:0) | 910.799 | 1.838 +/- 0.257 | 0.467 +/- 0.080 | 1.371 | 1.244 | 0.270 | -1.891 |
| 38 | TG (18:1_18:2_20:0) | 912.815 | 1.840 +/- 0.343 | 0.524 +/- 0.098 | 1.316 | 1.225 | 0.295 | -1.763 |
| 39 | PC (18:3_18:1) | 781.562 | 2.024 +/- 0.090 | 0.745 +/- 0.145 | 1.280 | 1.242 | 0.387 | -1.369 |
| 40 | TG (18:1_18:2_18:4) | 876.721 | 1.792 +/- 0.310 | 0.525 +/- 0.129 | 1.267 | 1.228 | 0.296 | -1.754 |
| 41 | PE (18:0_18:2) | 743.547 | 1.785 +/- 0.222 | 0.633 +/- 0.141 | 1.152 | 1.182 | 0.392 | -1.351 |
| 42 | TG (14:0_18:2_18:2) | 826.705 | 1.437 +/- 0.244 | 0.381 +/- 0.044 | 1.056 | 1.231 | 0.291 | -1.783 |

Note: the "consumption" column is exactly (0 s) - (8 min); the printed FC column is NOT exactly (8
min)/(0 s) of the printed means (row 1: 9.145/29.088 = 0.314 vs printed 0.329; row 3: 0.029 vs
0.034) — FC was evidently computed on another basis (raw intensities or mean of per-replicate
ratios). Use the concentrations, not FC.

### Numbers printed only in the text

- 0 s lipidome: 693 species, 28 subclasses; GP 57.78 % of content (PC 18.02, PE 15.80, PI 9.74,
  PS 8.56, PG 2.71, LPE 1.18 %), TG 33.42 %, **FFA 6.58 %**, DG 1.34 %, MG 1.10 %. 31 FFAs (16 SFA,
  15 USFA); "linoleic acid (LA, C18:2), oleic acid (C18:1), palmitic acid (C16:0), and stearic
  acid (C18:0) were the most abundant". Highest GP species at 0 s: PE (18:2_18:2) 29.088, PI
  (18:2_16:0) 16.123, PC (18:1_18:2) 13.224 ug/g. Lyso-lipids at 0 s: LPE 6.176, LPI 2.211, LPC
  0.578, LPA 0.4652, LPG 0.085 ug/g; LPE (18:2) 2.018, LPE (18:1) 1.826, LPE (16:0) 1.234 ug/g.
  Top TGs: TG (16:0_18:2_18:2) 16.266, TG (18:2_18:2_18:3) 14.457, TG (18:1_18:2_18:3) 13.612 ug/g.
- 0 s vs 8 min: 319 differential lipids (16 FFA, 59 GP, 194 TG; 307 down, 12 up). "LA showed the
  most pronounced reduction, decreasing from 10.658 ug/g to 0.309 ug/g." LPE (18:2), LPI (18:3),
  LPC (18:2), LPA (18:2) fell by 1.91, 0.95, 0.15, 0.14 ug/g. PE "consumption rate (up to 54.59 %)"
  (= 1 - 37.67/82.96 from Table 1, reproduces). Key-precursor class totals consumed: GPs 64.27,
  TGs 34.26, FFAs 17.39 ug/g (17.39 = 10.349 + 7.039 reproduces).
- Interval counts of differential lipids: 0-10 s 8 (all DGTS/LDGTS, down); 20-40 s 7 (FFA 16:1
  and PE 22:6_16:0 down, LPI 18:3 appears); 40-60 s 10; 60 s-2 min 35 (15 down, 20 up); 2-4 min
  11 (all down, incl. FFA 18:2, LPE 18:2, LPA 18:2); 4-8 min 39 (12 down, 27 up; 20 TGs, 17 up,
  saturated/monounsaturated).
- §3.7 correlations of lipid content with hexanal, hexanol and their sum across the 0-8 min
  series: all 42 key lipids r = -0.83 to -0.99; PE (18:2_18:2) -0.94, PE (18:2_16:0) -0.93, PE
  (18:1_18:2) -0.95; TG (18:2_18:2_18:3) -0.97, TG (16:0_18:2_18:2) -0.98, TG (18:1_18:2_18:3)
  -0.97, TG (18:0_18:2_18:3) -0.97; **FFA (18:2) vs (hexanal + hexanol) r = -0.96**.
- Abstract: 12 key precursors by > 5 ug/g depletion (Table 2 rows 1-12); 8 by correlation.

## 4. Numbers the repository can use

Conversions below assume pea milk density ~1.0 g/mL so that ug/g of supernatant ~ ug/mL = mg/L;
free linoleic acid MW 280.45. Both are arithmetic on printed values, not readings.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **Hexanal + hexanol time course** | FIGURE-ONLY; shape: rapid rise 0-2 min, slower rise 2-8 min, no increase after 8 min | (unit not recoverable) | 12.5 % w/w pea slurry, unbuffered, T not stated, quench 5 % citric acid | Fig. 1; §2.2, §3.1 | figure_only (the shape statements are text) |
| Time to plateau of C6 products, pea milk | 8 (no further increase beyond 8 min; "rapidly within the first 2 min") | min of enzyme time | as above | §2.2, §3.1 | measured (qualitative time scale) |
| Same, soybean, for contrast | ~1 (plateau "after just 1 min") | min | Trindler 2022, cited | §3.1 | secondary |
| **Free linoleic acid, 0 s** | 10.658 +/- 0.809 (~10.7 mg/L ~ **38 uM**) | ug/g supernatant | 0 s = citric acid added before grinding | Table 2 row 3 | level_only (substrate pool at t = 0; compare Zhang 2020b ~14.3 mg/L = 51 uM in the same slurry strength) |
| **Free linoleic acid, 8 min** | 0.309 +/- 0.056 (~1.1 uM) | ug/g | end of enzyme window | Table 2 row 3 | level_only (end of process) |
| Free linoleic acid consumed 0 -> 8 min | 10.349 (~36.9 uM); 97.1 % of the pool | ug/g | | Table 2 | within_study_ratio |
| Apparent first-order constant for free LA, two-point | ln(10.658/0.309)/8 min = **0.44 min^-1** (t1/2 ~ 1.6 min) | min^-1 | endpoints only; LA is also being released by lipase/phospholipase during the window, so this is a lower bound on the LOX consumption constant | derived from Table 2 | measured_rate (crude, two-point; flag 4) |
| **Total FFA time course** | 34.55, 34.12, 33.46, 29.63, 27.95, 23.15, 18.44, 13.53 (SDs in Table 1) at 0, 10, 20, 40, 60 s, 2, 4, 8 min | ug/g | n = 3 each | Table 1 | measured_rate (fittable 8-point substrate series; net of release) |
| Apparent first-order constant, total FFA, 0-8 min | ln(34.55/13.53)/8 = **0.12 min^-1** (t1/2 ~ 6 min) | min^-1 | net series; non-LA FFAs dilute the LA-specific rate | derived from Table 1 | measured_rate (crude) |
| Free oleic acid 0 s -> 8 min | 9.982 -> 2.944 | ug/g | | Table 2 row 10 | level_only + ratio (0.295) |
| **PE time course** | 82.96, 79.12, 71.40, 86.72, 82.72, 55.75, 44.29, 37.67 | ug/g | same eight times | Table 1 | measured_rate (fittable; the phospholipid substrate the authors single out) |
| PE (18:2_18:2) 0 s -> 8 min | 29.088 -> 9.145 (-19.942; r = -0.94 with C6 sum) | ug/g | | Table 2 row 1; §3.7 | level_only + within_study_ratio (0.314) |
| PC, PI, PS, TG time courses | Table 1 rows | ug/g | | Table 1 | measured_rate (fittable) |
| Total lipid 0 s -> 8 min | 525.12 -> 271.67 (-48 %) | ug/g | | Table 1 | within_study_ratio |
| Class consumption 0 -> 8 min (key precursors only) | GP 64.27, TG 34.26, FFA 17.39 | ug/g | 12 lipids of > 5 ug/g depletion | §3.5 | within_study_ratio |
| PE consumption fraction | 54.59 % | % | 0 -> 8 min | abstract, §3.4; Table 1 | within_study_ratio |
| Ceiling on C6 yield from free LA alone | 36.9 uM hexanal-equivalent (= 3.7 mg/L hexanal) if every consumed free-LA molecule gave one C6 | uM | arithmetic | derived from Table 2 | derived bound (compare Zhang 2020b: 5.43 uM C6 measured = 10.6 % of free LA) |
| Correlation, FFA (18:2) vs hexanal + hexanol | -0.96 | r (Pearson, 8 points) | | §3.7 | measured statistic (product-side info that survives in text) |
| Correlations, all 42 key lipids vs C6 | -0.83 to -0.99 | r | | §3.7 | measured statistic |
| Grinding conditions | seed:soak water 1:3 w/w, 4 C, 12 h; grind 1:7 w/w peas:water; quench 5 % citric acid; 3000 rpm 15 min 4 C; supernatant | — | | §2.2 | protocol |
| LOX activity | NOT MEASURED | — | | — | — |
| Peroxide value / conjugated dienes / hydroperoxides | NOT MEASURED | — | | — | — |
| Hexanal and hexanol calibration | five-point, IS 2-methyl-3-heptanone, matrix 2 % cow's milk diluted to pea-milk protein, SIM, R2 > 0.99; parameters in Table A.1 (not on disk) | — | | §2.5, §3.1 | method |

What this changes for LOX-01: the repository can now carry (i) a measured substrate-depletion time
course in pea milk (Table 1 FFA / PE / TG rows, eight points, triplicates), (ii) a t = 0 free-LA
pool that agrees with Zhang 2020b to within 1.4x, (iii) an authors' time scale for the product
(rapid in 2 min, plateau by 8 min) and (iv) a correlation of -0.96 tying FFA (18:2) to the C6 sum
— but still NO printed product level at any time, so the product-side rate remains unfitted.

## 5. Flags

1. **Grinding temperature, water temperature and slurry pH are unstated.** A Joyoung P165 soymilk
   machine is a heating blender; the paper never says the heater was off, and "ground till 12
   min" is long enough for a household cycle to warm the slurry. Any rate taken from Table 1 has
   an unknown temperature attached (best guess ambient-to-warm, rising with time). Do not attach a
   temperature to these numbers in the repository.
2. **The product time course is FIGURE-ONLY** (Fig. 1) and even its unit is not in the text; hexanal
   and hexanol are never printed separately. Table A.1 (calibration) and Fig. A2 (correlations)
   are supplementary and not on disk. No value from Fig. 1 may enter the repository.
3. **"Treat-X" is enzyme time, not grinding time.** Every sample was ground 12 min in total; the
   quench (5 % citric acid, pH not stated but ~2 at that strength) is assumed complete and
   instantaneous; residual LOX activity at pH ~2 during the remaining grinding is not tested.
   The 0 s sample had acid added before grinding, so "0 s" is the seed's native lipidome, not a
   milk that ever saw active LOX.
4. **The substrate series is net of release.** Free fatty acids are simultaneously produced (lipase
   / phospholipase A1, A2 — the authors document LP rises in the first 20 s and a total-lipid
   rise at 40 s) and consumed (LOX). The two-point LA constant 0.44 min^-1 and the eight-point
   total-FFA constant 0.12 min^-1 are apparent net constants, not LOX turnover; the true
   LOX-on-free-LA constant is at least as large. No hydroperoxide or HPL product pool was
   measured, so the pathway cannot be closed.
5. **Lipids are ug per g of centrifuged supernatant**, after gauze filtration and 3000 rpm / 15 min /
   4 C. Pellet lipid is discarded; changing particle size with grinding time moves lipid between
   phases (the 20 -> 40 s bump). Absolute lipid levels are therefore supernatant levels, and the
   early-time (0-60 s) flatness partly reflects this partition, not zero reaction. Quantification
   is "internal standard method" against the vendor MWDB library with unstated standards —
   treat absolute ug/g as semi-quantitative, within-series ratios as sound.
6. **Table 2 FC does not reproduce from the printed means** (§3 note). Use concentrations.
7. **Large SDs on some rows** (PS at 20 s 47.55 +/- 20.10; TG at 60 s 162.65 +/- 46.40; LPE at 8 min
   6.67 +/- 2.31); Duncan letters show no class other than FFA moves significantly before 2 min.
   A fit should weight by the printed SDs.
8. **The "1:7 peas to water" basis is ambiguous** (dry vs soaked seed). On a dry basis it equals
   Zhang 2020b's 12.5 % slurry and the two t = 0 free-LA pools agree (10.7 vs 14.3 mg/L), which
   supports the dry-basis reading; on a soaked-seed basis the slurry would be ~1.8x more dilute.
9. **Acquisitions this paper points to for LOX-01:** Yan et al. 2024 (Food Chem 445, 138696) — same
   lab, pea milk hexanal / hexanol calibrated levels, the calibration this paper reuses; Bi et al.
   2022 (Food Chem 380, 132203) — pea seed milk flavour "via enzyme activity inhibition", likely to
   carry LOX activity and inhibition data; Trindler et al. 2022 (Food Chem 376, 131892) — review,
   source of the soybean 1-min plateau claim; Feng et al. 2021 (JAFC 70, 289) — hexanal formation
   mechanism in soybean processing at the subcellular level (may hold a soybean time course).
