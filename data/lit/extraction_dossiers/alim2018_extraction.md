# Alim, Song, Liu, Zou, Zhang & Zhang 2018 — EXTRACTION (file stem `alim2018`)
J. Sci. Food Agric. **98**(10):3774-3783, DOI 10.1002/jsfa.8891. Crossref 2026-11-03: "Flavour-active compounds in
thermally treated yeast extracts", authors Aygul Alim, Huanlu Song, Ye Liu, Tingting Zou, Yu Zhang, Songpei Zhang;
JSFA 98(10):3774-3783, online 2018-03, print 2018-08. **Matches the PDF and the brief.** [!] **The file is the Wiley
"Accepted Article" manuscript** (29 pp., unpaginated, tables and figures at the end), not the typeset version of
record; page numbers below are PDF pages. Not the same paper as `alim2019_extraction.md` (LWT 2019, CAMOLA).

### One commercial yeast-extract paste (Angel Yeast LA00) diluted 1 : 5 (w/v) in water, heated 1 h at 100 / 110 / 120 / 130 / 140 °C in a foil-sealed vial inside a Parr high-pressure reactor (25 °C = unheated). Precursors after heating: 15 free amino acids (g/kg, Table 1), 5'-nucleotides, Cys, ornithine, thiamine, fructose / glucose / sucrose, peptide MW distribution (Fig. 1, figure only). Volatiles by SPME (35 rows, Table 2) and SAFE (14 rows, Table 3) in "µg/L" against **one non-isotopic internal standard (2-methyl-3-heptanone, RF = 1)**. Triplicate heatings, ± and letters printed. **No pH (before or after), no moisture, no time series, no SIDA; MFT, HDMF, sotolon, 2-AP, phenylacetaldehyde not detected.**

**Read 2026-11-03 from `data/articles/alim2018.pdf`** (29 pp., accepted manuscript, text layer present). Prose from
`pdftotext -layout` (scratchpad `ingest/alim2018_work/text.txt`). **Table 1 (p. 26), Table 2 (pp. 27-28), Table 3
(p. 29) read by eye from 130 dpi renders** (`p-26.png` ... `p-29.png`); the text layer agrees with the renders cell
for cell (the split methional row of Table 3 was resolved on the render). **Fig. 1 (p. 22) digitised from a 150 dpi
render** (`p-22.png`), read error ±0.02 mg/g (a), ±0.05 (b), ±0.3 mg/g (c), ±2 % (d). Figs. 2-4 (p. 23-25: sensory
spider plots, PLSR loadings) not digitised (no model-relevant numbers). Unit conversions by hand [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Beijing Advanced Innovation Center for Food Nutrition and Human Health, Laboratory of Molecular
Sensory Science, **Beijing Technology and Business University**; corresponding author **Huanlu Song**
(songhl@th.btbu.edu.cn). **Existing id `btbu_song`** ("BTBU (Song H)") **fits** (P ~0.97): same corresponding author,
e-mail, laboratory, Parr-reactor + foil-sealed vial + 2-methyl-3-heptanone IS design, and the **same YE paste LA00
(Hubei Angel Yeast)** as `alim2019`. Same-lab dossiers: `alim2019_extraction.md`, `wang2012_extraction.md`,
`xu2013_extraction.md`, `raza2020_extraction.md`. Not independent of them for leave-one-lab-out.

**Same pots as alim2019? No** (confidence ~95 %). Same YE paste, but alim2019's YE pot is YE + xylose + Cys + thiamine
at 135 °C / 90 min (amounts unprinted) and its other pots are synthetic peptides; alim2018 heats the YE alone, 1 : 5 in
water, 100-140 °C / 60 min. No shared condition, no shared table. alim2018's composition data (Table 1, Fig. 1) are
the only printed composition of LA00 on file (alim2019 prints none).

## 1. What this paper answers for the model

| model need | what Alim 2018 gives | tag |
|---|---|---|
| MFT in heated YE | **not detected** (text, p. 13) | [M] / [!] |
| FFT, thiophenes, DMTS, disulfides vs T | Table 2 (SPME) / Table 3 (SAFE), µg/L, 25-140 °C, 1 h | [M] (semi-quant.) |
| Strecker aldehydes (2-/3-methylbutanal, methional, benzaldehyde) vs T | Table 2 / 3, µg/L | [M] (semi-quant.) |
| phenylacetaldehyde | **not in Tables 2 / 3** although the text implies it was detected | [!] |
| furans / furanones | furfural, 5-methylfurfural, 2-acetylfuran, 2-furanmethanol, 2,3-butanedione; **no HDMF, sotolon** | [M] |
| pyrazines vs T | 7 pyrazines (Table 2), several non-monotonic in T | [M] |
| amino-acid loss vs T | Table 1: 15 free AA, g/kg ± SD, 25-140 °C | [M] |
| sugar loss vs T | Fig. 1c: fructose, glucose, sucrose, mg/g, figure only | [F] |
| nucleotide loss vs T (ribose source) | Fig. 1a: 5'-CMP/GMP/IMP/AMP, mg/g; text gives 25 and 140 °C values | [F] / [M] |
| Cys, thiamine, ornithine vs T | Fig. 1b; thiamine axis labelled "ng/g" (unit doubtful) | [F] / [!] |
| peptide MW distribution vs T | Fig. 1d, % in four MW classes | [F] |
| pH, moisture / aw, heat-up | **none** | [!] |
| time series | **none** (one 1 h hold per T) | — |

**Answers.**
1. **Which odorants?** (confidence ~90 % on the list; identifications are MS-library + RI + odour, several RIs do
   not fit the named compound, sec. 5). Sulfur: FFT, thiophene, 2-methylthiophene, 2-methyl-3-thiolanone
   (= 2-methyltetrahydrothiophen-3-one), 2-acetylthiophene (twice), 2-thiophenecarboxaldehyde, 5-methyl-2-
   thiophenecarboxaldehyde, DMTS, "dimethyl sulfone", 4,5-dimethylthiazole, ethyl isopropyl disulfide, diisopropyl
   disulfide, dipropyl trisulfide, methional. Strecker: 2-/3-methylbutanal, methional, benzaldehyde. Pyrazines: methyl-,
   2,6-dimethyl-, trimethyl-, 2-ethyl-5-methyl-, 2-ethyl-3,6-dimethyl-, 2,6-diethyl-3-methyl-, 2,5-dimethyl-3-(3-
   methylbutyl)pyrazine. Furans: furfural, 5-methylfurfural, 2-acetylfuran, 2-furanmethanol, γ-nonalactone; dicarbonyl
   2,3-butanedione. **Not detected: MFT, 2-AP, HDMF, 3-mercapto-2-pentanone, 3-mercapto-2-butanone, 5-methylfurfurylthiol**
   [M, text p. 13].
2. **Units?** (confidence ~95 %). Volatiles "Sample Concn (µg/L)" from Ci = C_IS × A_i / A_IS with one IS
   (2-methyl-3-heptanone), **i.e. semi-quantitative, response factor 1 for every compound**, no calibration, no
   isotopologues. The basis (per L of the 1 : 5 reacted mixture, presumably) is not stated; the IS "500 ng/g" claim
   fits SAFE (20 µL × 0.816 µg/µL in 30 mL ≈ 0.54 mg/L) but not SPME (1 µL in 5 mL ≈ 0.16 mg/L) [Z] [!]. SPME and SAFE
   disagree by ×10-14 for the same pot (methional 120 °C: 28.94 vs 399.32; benzaldehyde 142.29 vs 1782.55) [M].
   Precursors: absolute (HPLC / amino-acid analyser with standards): g/kg (amino acids), mg/g (nucleotides, sugars,
   Cys, ornithine), basis (paste as is, DW or reacted solution) **not stated** [!].
3. **Heating?** (confidence ~90 %). 3 g YE + 15 mL water (precursors) or 6 g + 30 mL (volatiles), foil-sealed vial
   in a Parr high-pressure reactor, 100-140 °C, 1 h each; triplicate. Heat-up, vial volume, headspace, pH not given.
   Volatile isolation adds heat: SPME 55 °C / 20 + 40 min; SAFE preceded by 8 h shaking with pentane at 40 °C.
4. **Non-volatiles with units?** Table 1 free amino acids g/kg ± SD (6 T); Fig. 1 nucleotides mg/g (25 °C: CMP 0.089,
   GMP 0.745, IMP 0.733, AMP 0.692; 140 °C: 0.01 / 0.19 / 0.155 / 0.077, printed in text), Cys / ornithine mg/g,
   thiamine "ng/g", sugars mg/g, peptide MW %.
5. **A second lab's quantitative YE odorant data comparable with Münch 1998 (same odorant, absolute units)?
   Partly; not absolute** (confidence ~85 %). Same odorants with numbers: 3-methylbutanal, 2-methylbutanal, methional,
   FFT, 2,3-butanedione (and benzaldehyde, pyrazines). But: (a) single non-matched IS with RF = 1, so the thiols and
   aldehydes can be off by ×10 or more in either direction (P ~0.6 of > ×3 bias for FFT); (b) unit basis unstated;
   (c) different YE, 1 : 5 water vs Münch 1 + 3, 1 h at ≤ 140 °C vs 20 min at 145 °C; (d) SPME vs SAFE disagree ×10.
   Order-of-magnitude comparison, assuming µg per L of mixture and ~75 % paste solids, i.e. ×(5-8) to µg/kg DW [Z]:
   140 °C SPME 3-methylbutanal 103 µg/L → ~500-820 µg/kg DW (Münch CYE 1456); 2-methylbutanal 42.65 → ~210-340 (965);
   methional 16-29 (120-140 °C) → ~80-230 (1095); FFT 7.15 → ~36-57 (580); MFT nd (530). **Strecker aldehydes within
   ×2-5 of Münch, thiols ×10-16 lower** [Z]. P(usable as a second-lab check within ×2) ~15 %; within ×10 ~55 %. It is
   a **temperature-ordering** data set (rank of T effects within one lab, one method), not an absolute benchmark.
6. **Is anything simulable?** No (confidence ~90 %): food matrix with unknown pH and water content, composition
   basis unstated, no time course. The precursor T series (amino acids, sugars, nucleotides after 1 h) is the most
   usable part, as a matrix plausibility check, and is itself doubtful (Leu −56 % after 1 h at 100 °C, sec. 5).

## 2. Conditions

- **Material** [M]: YE paste LA00, Hubei Angel Yeast Co. (Yichang). Solids, salt, pH, composition before heating
  (other than the 25 °C column) not printed.
- **Precursor pots** [M]: 3 g YE + 15 mL ultrapure water, foil-sealed vial, Parr reactor, 100 / 110 / 120 / 130 /
  140 °C, 1 h; filtered 0.22 µm, stored 4 °C; "carried out three time" [V].
- **Volatile pots** [M]: 6 g YE + 30 mL water, same reactor and T, 1 h, triplicate.
- **pH / buffer**: none stated; unbuffered water [!]. **Headspace / O2**: not stated. **Heat-up / cooling**: not stated.
- **Sensory** [M]: 8 panellists, 10-point scale, taste (Fig. 2a) and aroma (Fig. 2b) profiles — [Q] only here.

## 3. Analytical methods

- **Amino acids, nucleotides, sugars, peptide MW** [C]: "Liu JB's modified method" (ref. 9; OPA/FMOC HPLC for AA);
  Cys (refs 10, 11), ornithine (ref. 12), thiamine (ref. 13) by cited methods; no LOD, no recovery printed.
- **SPME** [M]: 5 mL sample + 1 µL IS (0.816 µg/µL) in 20 mL vial, 55 °C 20 min equilibration, CAR/DVB/PDMS 40 min.
- **SAFE** [M]: 30 mL sample + 20 µL IS + pentane (80 mL), shaken 8 h at 40 °C, SAFE 1 h at 10⁻⁴ Torr, dried,
  concentrated to 0.5 mL.
- **GC-O-MS** [M]: Agilent 7890A-7000B, DB-WAX and DB-5MS 30 m; 3 GC-O sniffers; identification NIST 05 + RI (± 20)
  + odour; "MS,RI,O" printed for nearly every row although several RIs are far from literature (sec. 5).
- **Quantification** [M]: Ci = C_IS × A_i / A_IS (Eq. 2), no response factors. ± = SD of triplicates (P ~0.8;
  not stated); letters = ANOVA groups.

## 4. Data

### 4a. Table 1 (p. 26): free amino acids after 1 h at T, g/kg, mean ± (SD), ANOVA letters [M]

| AA | 25 °C | 100 °C | 110 °C | 120 °C | 130 °C | 140 °C |
|---|---|---|---|---|---|---|
| Asp | 3.03 ± 0.13 b | 2.76 ± 0.14 b | 1.86 ± 0.18 b | 2.22 ± 0.63 b | 1.73 ± 0.18 b | 8.99 ± 2.46 a |
| Glu | 2.66 ± 0.01 b | 2.57 ± 0.03 b | 2.09 ± 0.17 b | 2.31 ± 0.20 b | 10.66 ± 2.73 a | 10.67 ± 2.59 a |
| Ser | 2.78 ± 0.06 a | 2.64 ± 0.03 ab | 2.40 ± 0.81 ab | 1.40 ± 0.14 c | 1.34 ± 0.11 c | 1.71 ± 0.39 bc |
| His | 2.58 ± 0.07 b | 2.43 ± 0.16 b | 1.56 ± 0.03 c | 1.18 ± 0.23 c | 3.35 ± 0.31 a | 2.73 ± 0.36 b |
| Gly | 2.79 ± 0.20 b | 2.54 ± 0.05 b | 1.79 ± 0.13 c | 0.92 ± 0.04 d | 0.87 ± 0.03 d | 3.82 ± 0.06 a |
| Thr | 3.36 ± 0.078 b | 2.89 ± 0.18 c | 2.86 ± 0.03 c | 2.43 ± 0.23 d | 0.83 ± 0.04 e | 5.02 ± 0.26 a |
| Arg | 3.02 ± 0.24 c | 2.04 ± 0.07 c | 1.06 ± 0.04 c | 0.54 ± 0.04 c | 7.17 ± 0.16 b | 15.22 ± 3.02 a |
| Ala | 2.75 ± 0.21 b | 2.66 ± 0.05 b | 2.28 ± 0.42 c | 0.98 ± 0.27 c | 5.28 ± 0.331 a | 6.47 ± 1.26 a |
| Tyr | 4.4 ± 0.14 a | 3.98 ± 0.63 a | 3.13 ± 0.23 a | 2.22 ± 0.54 ab | 2.07 ± 2.01 ab | 0.27 ± 0.05 b |
| Val | 3.24 ± 0.23 c | 3.17 ± 0.68 c | 0.82 ± 0.14 d | 5.32 ± 0.45 b | 4.56 ± 0.42 b | 7.60 ± 0.25 a |
| Met | 8.67 ± 0.35 c | 7.83 ± 0.51 c | 5.38 ± 0.11 b | 2.61 ± 0.067 b | 2.36 ± 0.18 a | 1.19 ± 0.12 d |
| Phe | 13.04 ± 0.56 d | 9.98 ± 0.66 a | 9.02 ± 0.85 ab | 6.05 ± 0.51 b | 5.41 ± 0.79 b | 10.53 ± 1.01 a |
| Ile | 12.14 ± 4.98 a | 7.93 ± 0.47 ab | 5.14 ± 0.03 b | 4.92 ± 0.02 b | 4.42 ± 0.02 b | 3.74 ± 0.75 b |
| Leu | 21.69 ± 2.14 a | 9.56 ± 0.40 b | 8.16 ± 0.90 b | 7.06 ± 0.68 b | 4.41 ± 0.04 c | 2.53 ± 0.58 c |
| Lys | 2.85 ± 0.49 a | 2.90 ± 0.83 a | 1.57 ± 0.41 b | 0.91 ± 0.09 b | 0.76 ± 0.08 b | 0.59 ± 0.02 b |

### 4b. Fig. 1 (p. 22), precursors after 1 h at T [F] (25 and 140 °C nucleotides printed in text [M])

| quantity (unit as printed) | 25 | 100 | 110 | 120 | 130 | 140 °C |
|---|---|---|---|---|---|---|
| 5'-CMP (mg/g) | 0.089 [M] | 0.09 | 0.07 | 0.07 | 0.05 | 0.01 [M] |
| 5'-GMP (mg/g) | 0.745 [M] | 0.60 | 0.49 | 0.46 | 0.30 | 0.19 [M] |
| 5'-IMP (mg/g) | 0.733 [M] | 0.63 | 0.50 | 0.43 | 0.31 | 0.155 [M] |
| 5'-AMP (mg/g) | 0.692 [M] | 0.55 | 0.44 | 0.34 | 0.25 | 0.077 [M] |
| cysteine (mg/g) | 1.53 | 1.50 | 1.47 | 1.47 | 1.27 | 1.09 |
| ornithine (mg/g) | 0.47 | 0.24 | 0.26 | 0.14 | 0.12 | 0.11 |
| thiamine ("ng/g", sic) | 1.72 | 0.53 | 0.35 | 0.18 | 0 (nd) | 0 (nd) |
| fructose (mg/g) | 12.0 | 10.0 | 7.2 | 7.0 | 6.4 | 5.7 |
| glucose (mg/g) | 7.6 | 6.9 | 7.1 | 6.8 | 5.9 | 4.5 |
| sucrose (mg/g) | 1.4 | 0.9 | 0.7 | 0.6 | 0.3 | ~0.05 |
| peptides > 1000 Da (%) | 33 | 32 | 29 | 27 | 20 | 16 |
| 500-1000 Da (%) | 28 | 25 | 26 | 26 | 27 | 29 |
| 300-500 Da (%) | 7 | 7 | 8 | 9 | 14 | 15 |
| < 300 Da (%) | 32 | 36 | 37 | 38 | 39 | 40 |

Derived [Z], treating each 1 h hold as isothermal (it is not; heat-up unknown): apparent first-order loss of 5'-IMP
k = ln(C25/CT)/60 min = 0.0025 (100), 0.0064 (110), 0.0089 (120), 0.0143 (130), 0.0259 min⁻¹ (140); two-point
Arrhenius 100 → 140 °C Ea ≈ 75 kJ/mol. Fructose −52 %, glucose −41 %, Cys −29 % at 140 °C / 1 h. Illustrative only.

### 4c. Table 2 (pp. 27-28): SPME volatiles, "µg/L", mean ± (SD) letters, 25 / 100 / 110 / 120 / 130 / 140 °C [M]
Selected model-relevant rows (all 35 rows read; the others are acetic acid, octanal, nonanal, γ-nonalactone,
acetophenone, phenylethyl alcohol, ethyl phenylacetate, 2-acetylpyrrole and indole, given below the table):

| no. | compound (RI Wax / DB-5) | 25 | 100 | 110 | 120 | 130 | 140 °C |
|---|---|---|---|---|---|---|---|
| 1 | 2-methylbutanal (904 / –) | 2.12 ± 0.711 | 4.76 ± 0.98 | 5.21 ± 0.45 | 6.52 ± 3.12 | 33.65 ± 3.57 | 42.65 ± 28.06 |
| 2 | 3-methylbutanal (909 / –) | 5.28 ± 1.9 | 5.84 ± 1.61 | 7.09 ± 0.16 | 18.32 ± 5.85 | 66.19 ± 6.42 | 103.07 ± 55.93 |
| 3 | 2,3-butanedione (981 / –) | 0 | 0 | 0.37 ± 0.07 | 0.59 ± 0.03 | 0 | 0 |
| 4 | thiophene (1080 / –) | 0 | 0 | 0 | 0 | 1.93 ± 0.08 | 3.28 ± 0.03 |
| 5 | 2-methylthiophene (1083 / 776) | 0 | 0.42 ± 0.07 | 0.37 ± 0 | 0.6 ± 0.03 | 20.53 ± 2.42 | 13.03 ± 0.54 |
| 6 | "dimethyl sulfone" (1136 / –) | 0 | 0 | 0.06 ± 0.02 | 0.2 ± 0.17 | 1.3 ± 0.63 | 0.51 ± 0.72 |
| 7 | methylpyrazine (1262 / –) | 0 | 0.35 ± 0.07 | 0.61 ± 0.09 | 0.92 ± 0.09 | 2.63 ± 0.95 | 5.47 ± 0.99 |
| 8 | trimethylpyrazine (– / 950) | 17.17 ± 2.35 | 7.22 ± 3.41 | 6.73 ± 0.1 | 4.76 ± 0.36 | 0 | 0 |
| 10 | 2,6-dimethylpyrazine (1343 / 860) | 17.64 ± 0.92 | 0.53 ± 0.31 | 23.45 ± 6.58 | 2.77 ± 0.04 | 6.76 ± 3.11 | 47.21 ± 11.58 |
| 11 | 2-ethyl-5-methylpyrazine (– / 996) | 0 | 0 | 0.56 ± 0.29 | 2.52 ± 1.76 | 1.69 ± 0.87 | 6.41 ± 5.35 |
| 13 | 2-ethyl-3,6-dimethylpyrazine (1350 / 1542) | 0 | 0 | 0 | 1.53 ± 0.46 | 7.43 ± 2.25 | 4.78 ± 2.51 |
| 14 | 4,5-dimethylthiazole (1371 / –) | 0 | 0 | 0 | 0 | 0.52 ± 0.49 | 0 |
| 15 | 2,6-diethyl-3-methylpyrazine (1389 / 1585) | 17 ± 2.76 | 13.09 ± 4.3 | 35.45 ± 8.02 | 8.09 ± 0.85 | 7.21 ± 0.95 | 43.37 ± 22.59 |
| 16 | dimethyl trisulfide (1398 / –) | 0 | 0 | 0 | 0 | 0 | 2.64 ± 0.58 |
| 17 | 2-methyl-3-thiolanone (1415 / –) | 0 | 0 | 0 | 0 | 23.29 ± 4.11 | 34.7 ± 2.93 |
| 18 | **2-furfurylthiol** (1428 / –) | 0 | 0 | 0 | 1.64 ± 0.4 | 3.03 ± 0.3 | 7.15 ± 1.57 |
| 20 | **methional** (1449 / 902) | 5.05 ± 0.58 | 5.65 ± 1.06 | 9.5 ± 4.37 | 28.94 ± 1.86 | 26.37 ± 10.3 | 16 ± 4.65 |
| 22 | 2-acetylfuran (1516 / 913) | 0 | 0.31 ± 0.06 | 0.43 ± 0.05 | 1.2 ± 0.94 | 7.1 ± 0.6 | 14.23 ± 1.96 |
| 23 | benzaldehyde (1519 / 907) | 16.88 ± 6.41 | 46.57 ± 32 | 75.85 ± 33.02 | 142.29 ± 14.96 | 103.91 ± 2.02 | 194.19 ± 80.06 |
| 24 | 5-methylfurfural (– / 959) | 0 | 0 | 0.29 ± 0.02 | 4.02 ± 1.34 | 23.45 ± 1.18 | 57.62 ± 13.41 |
| 25 | furfural (1556 / 1275) | 0 | 0 | 8.6 ± 0.28 | 9.29 ± 1.62 | 9.68 ± 2.03 | 17.05 ± 4.07 |
| 27 | 2-furanmethanol (1652 / 859) | 0 | 1.82 ± 0.82 | 6.75 ± 0.57 | 6.15 ± 1.89 | 16.86 ± 6.2 | 0 |
| 28 | 2-acetylthiophene (– / 1014) | 0 | 0 | 0 | 7.97 ± 7.98 | 9.84 ± 1.47 | 15.31 ± 5.82 |
| 29 | 2,5-dimethyl-3-(3-methylbutyl)pyrazine (1683 / 1247) | 0 | 3.68 ± 1.39 | 9.66 ± 1.41 | 15.66 ± 1.39 | 19.72 ± 3.43 | 28.83 ± 9.39 |
| 30 | 2-thiophenecarboxaldehyde (1689 / 1055) | 0 | 0.64 ± 0.48 | 1.34 ± 0.28 | 6.55 ± 0.65 | 55.54 ± 1.49 | 18.7 ± 2.33 |
| 31 | "2-acetylthiophene" (1807 / 1054) | 0 | 2.33 ± 0.42 | 7.34 ± 3.02 | 76.49 ± 4.26 | 174.6 ± 9.31 | 196.55 ± 108.9 |

Other rows [M] (25 → 140 °C): acetic acid 0 / 3.04 / 18.24 / 22.14 / 25.25 / 0; octanal 0 / 0.51 / 0.78 / 0.38 /
1.67 / 4.82; nonanal 3.56 / 6.2 / 12.51 / 5.15 / 9.68 / 0; γ-nonalactone ("dihydro-5-pentyl-2(3H)-furanone") 14.6 /
12.09 / 6.88 / 2.59 / 0.34 / 0; acetophenone 0 / 9.21 / 10.08 / 16.62 / 18.1 / 21.17; phenylethyl alcohol 1.62 / 9.28 /
9.33 / 3.41 / 0 / 0; ethyl phenylacetate 0.96 / 18.18 / 19.77 / 35.99 / 47.88 / 156.39; 2-acetylpyrrole 3.46 / 3.5 /
4.23 / 7.88 / 19.07 / 0; indole 21.01 / 22.8 / 37.61 / 84.06 / 232.58 / 338.07. (SDs on the renders; letters omitted.)

### 4d. Table 3 (p. 29): SAFE volatiles, "µg/L", 25 / 100 / 110 / 120 / 130 / 140 °C [M]

| no. | compound | 25 | 100 | 110 | 120 | 130 | 140 °C |
|---|---|---|---|---|---|---|---|
| 6 | ethyl isopropyl disulfide (MS, O) | 355.43 ± 0.84 | 440.95 ± 3.11 | 713.15 ± 66.75 | 787.66 ± 57.42 | 2237.04 ± 3.44 | 1444.91 ± 3.24 |
| 7 | diisopropyl disulfide (MS, O) | 0 | 427.62 ± 46.13 | 2403.54 ± 10.78 | 2841.04 ± 9.1 | 5643.16 ± 48.08 | 1433.6 ± 36.77 |
| 8 | **methional** | 0 | 0 | 103.99 ± 9.06 | 399.32 ± 16.9 | 30.6 ± 22.63 | 0 |
| 9 | benzaldehyde | 0 | 49.52 ± 18.29 | 292.17 ± 24.65 | 1782.55 ± 13.58 | 552.33 ± 295.75 | 1204.6 ± 6.05 |
| 10 | dipropyl trisulfide (MS, O) | 151.15 ± 6.92 | 575.03 ± 16.55 | 4355 ± 59.4 | 8881.62 ± 1.84 | 3376.93 ± 12.02 | 1328.77 ± 49.71 |
| 11 | 5-methyl-2-thiophenecarboxaldehyde | 0 | 0 | 0.01 ± 0 | 824.86 ± 44.55 | 16.25 ± 6 | 0 |
| 13 | furfural | 0 | 0 | 0 | 88.69 ± 1.42 | 374.14 ± 12.05 | 0 |

Other rows: ethanol 0 / 126.1 / 572.43 / 878.52 / 2316.53 / 0; undecanal 0 / 145.19 / 399.33 / 0 / 0 / 0; 3-hexanone
474.67 / 51.42 / 24.15 / 22.11 / 14.69 / 0; octanal 0 / 0 / 4.7 / 12.56 / 9.83 / 0; limonene 2878.56 at 130 °C only;
benzyl alcohol 0 / 69.66 / 91.39 / 76.42 / 44 / 0; indole 94.12 / 163.58 / 259.96 / 877.47 / 1034.96 / 0.

## 5. Flags and internal consistency

1. [!] **Semi-quantitative**: one ketone IS, RF = 1, for thiols, aldehydes, pyrazines, disulfides alike; unit basis
   unstated; SPME vs SAFE ×10-14 apart for the same pot (methional, benzaldehyde).
2. [!] **IS concentration inconsistent**: "500 ng/g" vs 1 µL of 0.816 µg/µL into 5 mL (≈ 163 µg/L) for SPME [Z].
3. [!] **Identifications doubtful** (P ~0.6 at least three rows misassigned): "dimethyl sulfone" at Wax RI 1136 (lit.
   ~1900; dimethyl disulfide ~1070 is likelier); trimethylpyrazine DB-5 950 (lit. ~1000); 2-ethyl-3,6-dimethylpyrazine
   DB-5 1542 (lit. ~1080); methional Wax 1449 in Table 2 vs 1569 in Table 3; **2-acetylthiophene listed twice** (rows 28
   and 31), row 31's RIs 1807 / 1054 identical to Table 3's 5-methyl-2-thiophenecarboxaldehyde. Isopropyl / propyl
   disulfides and trisulfides (MS, O only) are unusual Maillard products in YE; first-time claims rest on library
   matches.
4. [!] **SAFE 140 °C column is 0 for ethanol, methional, furfural, indole, benzyl alcohol** while 130 °C is high;
   indole by SPME rises to 338 at 140 °C. A lost or failed 140 °C SAFE run is likelier (P ~0.6) than chemistry.
5. [!] Non-monotonic SPME rows beyond SD: 2,6-dimethylpyrazine 17.64 → 0.53 → 23.45 → 2.77 → 6.76 → 47.21;
   trimethylpyrazine highest unheated and 0 at 130-140 °C.
6. [!] **Amino-acid losses implausible as chemistry**: Leu 21.69 → 9.56 g/kg (−56 %) after 1 h at 100 °C, Ile −35 %,
   Phe −23 %, while 3-methylbutanal rises by ~0.6 µg/L; Arg ×5 and Glu ×4 at 130-140 °C (peptide hydrolysis per the
   authors [Q]); Met letters run c, c, b, b, a, d against monotonic falling means; Phe 25 °C "d" highest. Mass balance
   impossible to check without the basis. Glu 2.66 g/kg is low for a yeast extract (~50-80 g/kg DW typical [C, general
   knowledge]), suggesting a per-reaction-solution basis or a different extract (P ~0.5).
7. [!] **Thiamine axis "ng/g"**: 1.7 ng/g would be ~1000 × below Mahadevan 2006's 11-30 µg/g DW; µg/g is likelier
   (P ~0.6) but still ~10 × lower. Unit unresolved.
8. [!] **Fig. 1c: sucrose is labelled a "reducing sugar"**; fructose 12 mg/g is ~60-100 × Mahadevan's 0.11-0.21 mg/g
   DW in other pastes (different product; possible added sugar).
9. [!] Text contradictions: SPME found "more" compounds than SAFE, then SPME's total is "lower" (p. 13); peptides
   < 300 Da "rose gradually" (p. 11, Fig. 1d) then "decreased" (p. 13); text implies phenylacetaldehyde,
   2-acetyl-2-thiazoline and 2,3-diethyl-5-methylpyrazine were detected, but they are in neither table.
10. [!] No pH anywhere; unbuffered YE at 140 °C for 1 h will acidify (acetic acid row; sensory "sour" at 130-140 °C [Q]).
11. Accepted manuscript, not the version of record: values may differ from the typeset JSFA article (P ~0.1).

## 6. Existing records in the repo that cite this paper

grep 2026-11-03 for `jsfa.8891`, `alim2018`, "Flavour-active compounds in thermally", "Alim ... 2018":
- `data/keys/pdf_dossier_aliases.yml` l. 99 (pending line): Alim, Song et al. 2018 JSFA, DOI, 29 pp., "aroma-active vs
  non-volatile compounds of thermally reacted yeast extracts" — **match**.
- `data/lit/extraction_dossiers/alim2019_extraction.md` l. 6: "companion `alim2018` (JSFA 2018, thermally treated
  yeast extracts; queued, not read here)" — **match**.
- `data/research_corpus/raw/Gemini_deep_research_06.md` l. 256 (bibliography item 2, ResearchGate link) cited at l. 15
  for cysteine Strecker → H2S + ammonia → "highly potent meaty odorants like 2-methyl-3-furanthiol and various
  thiophenes" — **mismatch on MFT: the paper did not detect MFT** (p. 13); thiophenes yes; the H2S mechanism is not
  shown here [P]. Raw research corpus, not a curated record; no runtime effect.
- `docs/history/slr_benchmark_evaluation_2026-03.md` l. 109 cites Kang, Alim & Song 2018 (J Chromatogr B) — a
  different paper.
- No ledger, registry, benchmark, prior or code record cites it.

## 7. Proposed ledger entries (not written into the ledgers)

```yaml
# sulfur.yml
  # ---- alim2018 : Alim et al. 2018 (yeast-extract paste LA00 heated 1 h at 100-140 C; SPME / SAFE) (ingested 2026-11-03) ----
  - id: alim2018_tab2_tab3_sulfur_volatiles_temperature
    dossier: alim2018_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang & Zhang 2018 J Sci Food Agric 98:3774"
    lab: "BTBU (Song H)"
    system: "yeast-extract paste LA00 (Angel Yeast) 6 g + 30 mL water"
    medium: "unbuffered water, foil-sealed vial in Parr high-pressure reactor; food matrix"
    temperature_c: [25, 100, 110, 120, 130, 140]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: ["FFT", "thiophene", "2-methylthiophene", "2-methyl-3-thiolanone", "2-acetylthiophene", "2-thiophenecarboxaldehyde", "5-methyl-2-thiophenecarboxaldehyde", "dimethyl trisulfide", "4,5-dimethylthiazole", "ethyl isopropyl disulfide", "diisopropyl disulfide", "dipropyl trisulfide"]
    kind: endpoint
    location: "Table 2 (SPME) and Table 3 (SAFE) sulfur rows (dossier secs. 4c, 4d)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "let it react at different temperatures (100 C, 110 C, 120 C, 130 C, 140 C) for 1 h"
    flags: [matrix_food, ph_not_stated, single_is_no_response_factor, unit_basis_not_stated, "is_concentration_inconsistent: 500 ng/g vs 1 uL x 0.816 ug/uL in 5 mL", "spme_vs_safe_10x_disagreement", "identification_doubtful_rows", "safe_140C_column_zeros", accepted_manuscript_not_vor, "mft_not_detected", "complements: alim2018_tab2_tab3_strecker_pyrazine_temperature"]
    note: "T ordering within one lab: FFT first seen at 120 C (1.64 -> 7.15 at 140), thiophenes jump at 130 C. MFT not detected (Muench 1998 CYE 530 ug/kg DW). Not an absolute second-lab check (RF = 1, basis unstated)."
# aroma_other.yml
  - id: alim2018_tab2_tab3_strecker_pyrazine_temperature
    dossier: alim2018_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang & Zhang 2018 J Sci Food Agric 98:3774"
    lab: "BTBU (Song H)"
    system: "yeast-extract paste LA00 (Angel Yeast) 6 g + 30 mL water"
    medium: "unbuffered water, foil-sealed vial in Parr high-pressure reactor; food matrix"
    temperature_c: [25, 100, 110, 120, 130, 140]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: [2-methylbutanal, 3-methylbutanal, methional, benzaldehyde, methylpyrazine, "2,6-dimethylpyrazine", trimethylpyrazine, "2-ethyl-5-methylpyrazine", "2-ethyl-3,6-dimethylpyrazine", "2,6-diethyl-3-methylpyrazine", "2,5-dimethyl-3-(3-methylbutyl)pyrazine"]
    kind: endpoint
    location: "Table 2 (SPME) and Table 3 (SAFE) (dossier secs. 4c, 4d)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "let it react at different temperatures (100 C, 110 C, 120 C, 130 C, 140 C) for 1 h"
    flags: [matrix_food, ph_not_stated, single_is_no_response_factor, unit_basis_not_stated, "spme_vs_safe_10x_disagreement", "pyrazine_rows_non_monotonic", "identification_doubtful_rows", accepted_manuscript_not_vor, "phenylacetaldehyde_not_in_tables", "complements: alim2018_tab2_tab3_sulfur_volatiles_temperature", "also_relevant: trunk"]
    note: "Strecker aldehydes rise 5-20x from 25 to 140 C (3-methylbutanal 5.28 -> 103.07 ug/L SPME); methional peaks at 120 C. Order-of-magnitude vs Muench 1998 CYE within x2-5 after an assumed x5-8 basis conversion (dossier sec. 1 ans. 5)."
  - id: alim2018_tab1_fig1_yeast_extract_precursors_temperature
    dossier: alim2018_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang & Zhang 2018 J Sci Food Agric 98:3774"
    lab: "BTBU (Song H)"
    system: "yeast-extract paste LA00 (Angel Yeast) 3 g + 15 mL water"
    medium: "unbuffered water, foil-sealed vial in Parr high-pressure reactor; food matrix"
    temperature_c: [25, 100, 110, 120, 130, 140]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: [free amino acids (15), cysteine, ornithine, thiamine, "5'-CMP", "5'-GMP", "5'-IMP", "5'-AMP", peptide MW classes]
    kind: endpoint
    location: "Table 1 (amino acids, g/kg); Fig. 1a, 1b, 1d (digitised, dossier sec. 4b); nucleotide 25 / 140 C values in text"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "let it react in particular temperatures ... for 1 h"
    flags: [matrix_food, ph_not_stated, unit_basis_not_stated, "amino_acid_losses_implausible_leu_minus56pct_100C", "thiamine_unit_ng_per_g_doubtful", "anova_letters_inconsistent_met_phe", accepted_manuscript_not_vor, "also_relevant: sulfur"]
    note: "Only printed composition of LA00 (alim2019 prints none). 5'-IMP 0.733 -> 0.155 mg/g after 1 h at 140 C (apparent k 0.0025-0.026 /min, 100-140 C, Ea ~75 kJ/mol two-point, isothermal assumed) [dossier Z]. Cys -29 % at 140 C."
# trunk.yml
  - id: alim2018_fig1c_tab2_sugars_furans_temperature
    dossier: alim2018_extraction.md
    paper: "Alim, Song, Liu, Zou, Zhang & Zhang 2018 J Sci Food Agric 98:3774"
    lab: "BTBU (Song H)"
    system: "yeast-extract paste LA00 (Angel Yeast) in water, 1 : 5 (sugars: 3 g + 15 mL; furans: 6 g + 30 mL)"
    medium: "unbuffered water, foil-sealed vial in Parr high-pressure reactor; food matrix"
    temperature_c: [25, 100, 110, 120, 130, 140]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60]
    species: [fructose, glucose, sucrose, furfural, 5-methylfurfural, 2-acetylfuran, 2-furanmethanol, "2,3-butanedione", acetic acid]
    kind: endpoint
    location: "Fig. 1c (sugars, digitised, dossier sec. 4b); Table 2 / Table 3 furan and dicarbonyl rows (secs. 4c, 4d)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "let it react at different temperatures (100 C, 110 C, 120 C, 130 C, 140 C) for 1 h"
    flags: [matrix_food, ph_not_stated, unit_basis_not_stated, "sucrose_called_reducing_sugar", "sugars_calibrated_furans_semi_quantitative", "safe_140C_column_zeros", accepted_manuscript_not_vor]
    note: "Fructose 12.0 -> 5.7, glucose 7.6 -> 4.5 mg/g after 1 h at 140 C (read +-0.3). Furfural / 5-methylfurfural / 2-acetylfuran rise from 110-120 C (SPME, RF = 1). No HDMF, sotolon. Matrix plausibility only."
```

## 8. Supplementary material

None referenced (no "Supporting Information", "Table S" or appendix in the 29-page accepted manuscript; grep
2026-11-03). The version of record (JSFA 98:3774) might differ in layout; nothing indicates extra data. Not needed.
