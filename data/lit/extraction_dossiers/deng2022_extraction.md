# Deng 2022 — EXTRACTION (methionine + glucose 0.2 + 0.2 mol/L vs the Met-Glc Amadori compound 0.2 mol/L, initial pH 7.5 unbuffered, sealed vials, 100/120/130 C, 30-240 min; methional and eleven pyrazines by HS-SPME-GC-MS, semi-quantitative)
### The only paper on disk with a printed methional time course from methionine + a sugar at cooking temperature; the chain past methional (methanethiol, DMDS, DMTS) is invoked to explain the fall but never measured.

**Source on disk:** `data/articles/Deng2022.pdf` (10 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Deng2022.txt`); Tables 1 and 2 came through clean (value, SD and
significance letter in order for every cell) and are re-typed below; the column sums of Table 1 were
re-computed as a check (see Flags 4). Figures 1-4 (browning and pH; methional and total pyrazines at
three temperatures; GO, MGO, DA time courses; methionine) were not read and no value is taken from
them. The Supporting Information (Fig. S1 mechanism, S2-S3 ARP identity, S4 HPLC chromatogram, S5
sensory, Table S1 OAVs) is NOT on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Comparison of pyrazines formation in methionine/glucose and corresponding Amadori rearrangement product model" |
| Authors | Shibin Deng, Heping Cui, Khizar Hayat, Yun Zhai, Qiang Zhang, Xiaoming Zhang*, Chi-Tang Ho* (Jiangnan University, Wuxi; Putian University; Miami University; Anhui Qiangwang; Rutgers) |
| Venue | Food Chemistry 382 (2022) 132500. Received 18 October 2021, revised 14 February 2022, accepted 15 February 2022, online 19 February 2022 |
| DOI | 10.1016/j.foodchem.2022.132500 |
| Naming | Met/Glc = the binary mixture (the table headers write "Glc/Met", the figures "MG"); MG-ARP = N-(1-deoxy-D-fructos-1-yl)-methionine, the methionine-glucose Amadori compound (figures "MA"); GO = glyoxal, MGO = methylglyoxal, DA = diacetyl (2,3-butanedione); DL-methionine (racemic) throughout; "pyrazine" in the tables is the unsubstituted parent |
| Companions | Deng 2021 (JAFC 69, 5167; ARP preparation and the Met HPLC-ELSD method, not on disk); Zhou 2024 (same group, `zhou2024_extraction.md`, the fed-dicarbonyl pyrazine rates); Yu & Ho 1995 (`yu1995_extraction.md`, the methional -> methanethiol / DMDS / acrolein attribution cited here) |

## 1. Why it matters

Programme 6 (`tasks/roadmap_for_scientists.md` 5c) wants methionine as a distinguishable amino acid
on the sugar path with methional as its Strecker product and the chain methional -> methanethiol ->
DMDS / DMTS on the sulfur lane. The B19 reading log (methionine chain) has the chain documented but
no rates from a methionine + sugar pot. This paper prints, in Table 1, methional at five times (30 to
240 min) at 120 C from Met + glucose in water and from the Amadori compound, and prints in Tables 1
and 2 the pyrazine distribution at 100, 120 and 130 C. Methional rises and then falls in both pots
(Met/Glc: 87.77 µg/L at 180 min then 77.68 at 240; MG-ARP: 196.59 µg/L at 120 min then 85.60 at
240): the fall is what the chain to methanethiol would do, and the authors say so (citing Yu & Ho
1995), but methanethiol, DMDS and DMTS are not measured. So the paper gives (i) a methional level
series at one temperature and two precursor forms, (ii) a within-study ratio ARP vs free amino acid,
(iii) a pyrazine distribution from a methionine pot (which pyrazines a methionine system makes, and
that the parent pyrazine dominates in Met/Glc while 2,5-dimethylpyrazine dominates in the ARP pot),
and (iv) for B18, an order-of-magnitude of how much pyrazine a glucose pot's own dicarbonyl supply
yields against Zhou 2024's fed-glyoxal rate. Everything is semi-quantitative (response factor 1
against 1,2-dichlorobenzene), so levels validate at order of magnitude and ratios within the study
are the firmer quantities.

## 2. Methods as they matter to a model

- **Pots.** "Aqueous solutions of Met/Glc (each reagent at 0.2 mol/L) and MG-ARP model (0.2 mol/L)
  were prepared, and their pH values were adjusted to 7.5 ± 0.1 with NaOH (6 mol/L). The model
  solutions were heated in sealed vials at 100 ± 1, 120 ± 1, or 130 ± 1 C for different times (0,
  30, 60, 120, 180, and 240 min) with continuous stirring." Control: 0.2 mol/L Met alone, same
  conditions. Ice-water quench. So **[Met] = [Glc] = 200 mmol/L** (ten times Zhou 2024's 20 mmol/L),
  **[MG-ARP] = 200 mmol/L**, water, **no buffer** (NaOH only; the phosphate in 2.7 is the OPD
  reagent, not the pot). Vial volume, fill and headspace are not stated. The pH falls during heating
  (Fig. 1, figure-only; "much faster" in the ARP pot).
- **MG-ARP preparation.** Met + Glc 0.2 mol/L each, pH 7.5, 95 C water bath then vacuum (25 mbar)
  dehydration 5 + 20 min, dissolved, filtered (insoluble Met removed), Dowex 50WX4 H+ column,
  freeze-dried; **purity 96 %** (2.1); identity by UPLC-Q-TOF/MS and 1H/13C NMR (Fig. S2-S3).
- **Volatiles (methional, pyrazines).** HS-SPME-GC-MS: **2 mL** of product in a 20-mL vial (no salt
  mentioned) + 5 µL 1,2-dichlorobenzene 0.018 mg/mL in methanol = **0.09 µg internal standard per
  vial**; 50 C 10 min equilibration, DVB/CAR/PDMS 50/30 µm 50 C 30 min with stirring; desorption
  250 C 7 min, split 20:1; RTX-WAX 30 m x 0.25 mm x 0.25 µm; He 1.5 mL/min; EI 70 eV, m/z 35-500.
  Identification: NIST 17 + RI (C7-C30 alkanes) and, for eight compounds, authentic standards ("S":
  pyrazine, methylpyrazine, 2,5-, 2,6-, 2,3-dimethylpyrazine, trimethylpyrazine, acetylpyrazine,
  methional); ethylpyrazine, 2-ethyl-5-methylpyrazine, vinylpyrazine, 2-ethyl-3,5-dimethylpyrazine,
  2-vinyl-5-methylpyrazine and 2-ethyl-3-(methylthio)pyrazine by RI + MS only.
  **Quantification: semi-quantitative, response factor assumed 1** ("The compounds were
  semi-quantified by comparing their peak areas with that of the internal standard
  1,2-dichlorobenzene, assuming a response factor of 1 according to ... Kocadagli et al. (2021)");
  W_i (µg/L) = f' x A_i x m_s / (A_s x V) with f' = 1, m_s = 0.09 µg, V = 0.002 L, i.e. **W_i = 45
  µg/L x (A_i / A_s)**. No calibration curve, no LOD/LOQ for any volatile. All Table 1 and 2 values
  are therefore 1,2-dichlorobenzene-equivalents on this fibre and column.
- **alpha-Dicarbonyls (GO, MGO, DA).** 0.5 mL sample + 0.5 mL 0.1 mol/L phosphate pH 7.0 with 0.5 %
  OPD and 11 mmol/L DTPA, 4 h dark, room temperature; HPLC-PDA 315 nm, Sunfire C18, formic
  acid/methanol gradient; external calibration in water on the quinoxalines: GO y = 2.0954e6 x -
  6.8425e3 (R2 0.9999, LOD 0.0088, LOQ 0.0212 mmol/L); MGO y = 1.9055e6 x + 1.9028e3 (R2 0.9996, LOD
  0.0039, LOQ 0.0158 mmol/L); DA y = 3.1593e6 x - 3.0827e4 (R2 0.9996, LOD 0.0029, LOQ 0.0081
  mmol/L); x in mmol/L. Results only in Fig. 3 (figure-only).
- **Methionine.** HPLC-ELSD, XBridge Amide, water/acetonitrile with 0.1 % formic acid, external
  calibration log y = 1.4745 log x + 4.6812 (R2 0.9981), LOD 0.0897, LOQ 0.3018 mmol/L (sample
  diluted "to appropriate folds"). Results only in Fig. 4 (figure-only).
- **Browning** A420 after 70-fold dilution; **pH** by meter; both Fig. 1 (figure-only).
- **Replicates.** "All the results were the averages of three replicates"; mean ± SD; ANOVA p < 0.05
  (letters in the tables compare all ten columns of a row).
- **Unit conversions used below.** Methional M = 104.17 g/mol: 1 µg/L = 9.60e-3 µmol/L. Pyrazine
  80.09; methylpyrazine 94.12; dimethylpyrazines 108.14; trimethylpyrazine 122.17;
  2-ethyl-3,5-dimethylpyrazine 136.19 g/mol. Time in the tables is minutes.

## 3. Tables re-typed

### Table 1. "Overview of the methional and pyrazines concentrations (µg/L) for Met/Glc and MG-ARP (0.2 mol/L) model reactions at initial pH 7.5 and 120 C"

Mean ± SD, µg/L; letters a-g compare the ten cells of a row (p < 0.05); "–" = not detected. RI on
RTX-WAX; ID = identification (RI, MS, S = authentic standard).

| compound | RI | ID | Glc/Met 30 | 60 | 120 | 180 | 240 | MG-ARP 30 | 60 | 120 | 180 | 240 |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Methional** | 1472 | RI, MS, S | 18.52 ± 1.64 e | 25.27 ± 0.65 e | 75.76 ± 1.06 c | 87.77 ± 5.99 c | 77.68 ± 2.51 c | 25.42 ± 2.71 e | 54.14 ± 4.33 d | 196.59 ± 9.96 a | 153.85 ± 10.09 b | 85.60 ± 4.23 c |
| Pyrazine | 1215 | RI, MS, S | 4.96 ± 0.76 c | 6.31 ± 0.36 b | 7.16 ± 0.56 ab | 7.68 ± 0.57 a | 7.42 ± 0.36 a | – | 0.33 ± 0.05 d | 0.69 ± 0.06 d | 0.73 ± 0.07 d | 0.44 ± 0.07 d |
| Methylpyrazine | 1268 | RI, MS, S | 2.52 ± 0.13 d | 3.14 ± 0.13 c | 3.57 ± 0.27 b | 3.64 ± 0.11 b | 5.15 ± 0.12 a | 0.98 ± 0.10 e | 2.95 ± 0.25 d | 3.80 ± 0.20 b | 3.87 ± 0.13 b | 2.78 ± 0.09 d |
| 2,5-Dimethylpyrazine | 1329 | RI, MS, S | 2.42 ± 0.12 e | 3.89 ± 0.09 c | 3.83 ± 0.19 c | 2.98 ± 0.27 d | 2.28 ± 0.31 e | 4.60 ± 0.31 b | 5.34 ± 0.19 a | 5.70 ± 0.23 a | 4.26 ± 0.23 bc | 3.01 ± 0.14 d |
| 2,6-Dimethylpyrazine | 1337 | RI, MS, S | 0.28 ± 0.06 f | 0.57 ± 0.04 e | 0.65 ± 0.04 de | 0.59 ± 0.08 e | 0.81 ± 0.09 cd | 0.38 ± 0.08 f | 0.92 ± 0.08 c | 1.21 ± 0.07 b | 1.52 ± 0.10 a | 1.10 ± 0.08 b |
| Ethylpyrazine | 1341 | RI, MS | 0.29 ± 0.07 b | 0.67 ± 0.08 a | 0.34 ± 0.02 b | – | – | – | – | – | – | – |
| 2,3-Dimethylpyrazine | 1355 | RI, MS, S | – | 0.59 ± 0.07 a | 0.20 ± 0.03 b | 0.10 ± 0.03 c | 0.17 ± 0.02 bc | – | – | – | – | – |
| 2-Ethyl-5-methylpyrazine | 1403 | RI, MS | – | 1.03 ± 0.12 a | 0.51 ± 0.05 b | 0.41 ± 0.07 b | 0.68 ± 0.18 b | – | – | – | – | – |
| Trimethylpyrazine | 1417 | RI, MS, S | – | 1.03 ± 0.11 a | 0.80 ± 0.04 ab | 0.51 ± 0.04 b | 0.75 ± 0.23 ab | – | – | – | – | – |
| Vinylpyrazine | 1455 | RI, MS | – | 0.78 ± 0.06 b | 1.36 ± 0.16 a | 0.52 ± 0.18 cd | 0.61 ± 0.07 bc | – | – | – | 0.34 ± 0.03 de | 0.15 ± 0.03 e |
| 2-Ethyl-3,5-dimethylpyrazine | 1464 | RI, MS | 3.37 ± 0.31 a | 3.35 ± 0.07 a | 2.45 ± 0.20 c | 1.41 ± 0.22 d | – | – | – | – | – | – |
| Acetylpyrazine | 1666 | RI, MS, S | 0.48 ± 0.07 b | 0.58 ± 0.05 b | 1.16 ± 0.10 a | 1.05 ± 0.07 a | 1.17 ± 0.09 a | – | – | – | – | – |
| Total of pyrazines | | | 14.32 ± 1.51 c | 21.94 ± 1.15 a | 22.03 ± 1.14 a | 18.90 ± 1.30 b | 19.04 ± 0.40 b | 5.95 ± 0.49 g | 9.19 ± 0.57 ef | 11.40 ± 0.57 d | 10.71 ± 0.52 de | 7.48 ± 0.35 fg |
| Total pyrazines minus pyrazine | | | 9.36 ± 0.75 de | 15.63 ± 0.79 a | 14.87 ± 0.70 a | 11.22 ± 0.93 bc | 11.62 ± 0.59 b | 5.95 ± 0.49 f | 8.85 ± 0.52 e | 10.71 ± 0.51 bcd | 9.98 ± 0.45 cde | 6.04 ± 0.28 f |

**Arithmetic check (mine).** Column sums of the eleven pyrazine rows: Glc/Met 14.32 / 21.94 / 22.03 /
18.89 / 19.04 (printed totals 14.32 / 21.94 / 22.03 / 18.90 / 19.04: agree); MG-ARP 5.96 / **9.54** /
11.40 / 10.72 / 7.48 (printed 5.95 / **9.19** / 11.40 / 10.71 / 7.48). "Total minus pyrazine" MG-ARP
240 min: the eleven rows give **7.04**, printed **6.04**; the printed total 7.48 = 7.04 + 0.44
supports 7.04 (a typo in the last cell). The MG-ARP 60-min total (9.19 vs 9.54) does not close either
way; the rows are the values to keep.

**Methional in molar units (Table 1, 120 C).** Glc/Met: 0.178 / 0.243 / 0.727 / 0.843 / 0.746
µmol/L at 30 / 60 / 120 / 180 / 240 min. MG-ARP: 0.244 / 0.520 / 1.887 / 1.477 / 0.822 µmol/L. Against
200 mmol/L methionine (or ARP) the peak is 1e-5 of the precursor: the pot is never depleted by this
product.

### Table 2. "Overview of the methional and pyrazines concentrations (µg/L) for Met/Glc and MG-ARP (0.2 mol/L) model reactions at initial pH 7.5 at temperatures of 100 and 130 C"

Despite the title, **Table 2 has no methional row** (methional at 100 and 130 C is Fig. 2a only) and
no total rows; ethylpyrazine and acetylpyrazine rows are absent. Same columns and conventions as
Table 1.

Pyrazines at **100 C**:

| compound | RI | ID | Glc/Met 30 | 60 | 120 | 180 | 240 | MG-ARP 30 | 60 | 120 | 180 | 240 |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pyrazine | 1215 | RI, MS, S | – | 0.88 ± 0.08 c | 1.88 ± 0.09 b | 4.96 ± 0.24 a | 4.44 ± 0.38 a | – | – | – | – | – |
| Methylpyrazine | 1268 | RI, MS, S | – | – | – | 1.34 ± 0.11 a | 0.95 ± 0.10 bc | – | – | 1.04 ± 0.08 b | 0.82 ± 0.04 c | 0.97 ± 0.08 bc |
| 2,5-Dimethylpyrazine | 1329 | RI, MS, S | – | – | – | 1.12 ± 0.09 ab | 1.33 ± 0.09 a | – | – | 1.32 ± 0.11 a | 1.07 ± 0.14 b | 1.34 ± 0.06 a |

Pyrazines at **130 C**:

| compound | RI | ID | Glc/Met 30 | 60 | 120 | 180 | 240 | MG-ARP 30 | 60 | 120 | 180 | 240 |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pyrazine | 1215 | RI, MS, S | 7.42 ± 0.33 c | 7.87 ± 0.25 c | 14.43 ± 1.37 a | 11.34 ± 1.16 b | 8.43 ± 0.56 c | 0.09 ± 0.01 d | 0.35 ± 0.03 d | 0.61 ± 0.05 d | 0.92 ± 0.09 d | 1.23 ± 0.12 d |
| Methylpyrazine | 1268 | RI, MS, S | 3.39 ± 0.31 ef | 5.11 ± 0.16 c | 10.74 ± 1.70 a | 9.48 ± 0.36 a | 7.19 ± 0.18 b | 1.75 ± 0.13 g | 2.70 ± 0.24 fg | 4.47 ± 0.42 cde | 3.67 ± 0.20 def | 4.89 ± 0.53 cd |
| 2,5-Dimethylpyrazine | 1329 | RI, MS, S | 10.02 ± 0.82 c | 10.14 ± 0.72 c | 16.66 ± 0.71 a | 12.44 ± 0.44 b | 9.85 ± 0.68 c | 3.79 ± 0.17 e | 5.02 ± 0.33 d | 5.20 ± 0.23 d | 5.09 ± 0.46 d | 6.24 ± 0.56 d |
| 2,6-Dimethylpyrazine | 1337 | RI, MS, S | – | – | – | – | – | 0.91 ± 0.07 b | 1.38 ± 0.11 a | 1.47 ± 0.16 a | 1.43 ± 0.12 a | 1.67 ± 0.16 a |
| 2,3-Dimethylpyrazine | 1355 | RI, MS, S | – | – | – | – | – | – | – | 0.93 ± 0.06 a | 0.82 ± 0.08 a | 0.85 ± 0.05 a |
| 2-Ethyl-5-methylpyrazine | 1403 | RI, MS | 1.95 ± 0.16 b | 1.95 ± 0.17 b | 3.79 ± 0.36 a | 3.49 ± 0.40 a | 2.61 ± 0.29 b | – | – | – | – | – |
| Trimethylpyrazine | 1417 | RI, MS, S | – | – | 3.09 ± 0.20 a | 3.08 ± 0.18 a | 2.23 ± 0.11 b | 0.69 ± 0.03 d | 0.99 ± 0.09 cd | 1.16 ± 0.09 c | 1.17 ± 0.13 c | 1.99 ± 0.26 b |
| Vinylpyrazine | 1455 | RI, MS | 1.49 ± 0.19 d | 2.59 ± 0.18 c | 3.95 ± 0.34 a | 3.81 ± 0.35 a | 2.43 ± 0.12 c | – | – | – | – | 3.13 ± 0.16 b |
| 2-Ethyl-3,5-dimethylpyrazine | 1464 | RI, MS | 4.67 ± 0.17 b | 4.71 ± 0.25 b | 5.67 ± 0.47 a | 4.94 ± 0.44 b | 3.56 ± 0.23 c | 0.70 ± 0.06 d | 1.20 ± 0.10 d | 1.13 ± 0.07 d | 0.99 ± 0.10 d | 1.23 ± 0.06 d |
| 2-Vinyl-5-methylpyrazine | 1512 | RI, MS | – | – | – | 3.17 ± 0.11 a | 2.59 ± 0.17 b | – | – | 2.47 ± 0.25 b | 2.22 ± 0.19 b | 2.50 ± 0.18 b |
| 2-Ethyl-3-(methylthio)pyrazine | 1955 | RI, MS | – | – | – | – | – | – | – | – | 1.72 ± 0.13 | – |

### Numbers in the running text

- Methional from **Met alone** (0.2 mol/L, same conditions): "1-2 µg/L (Data not shown)" — the
  non-Maillard route is negligible.
- MG-ARP purity 96 %.
- Qualitative only (Figs. 2-4, all FIGURE-ONLY): total pyrazines and methional both rise with
  temperature in both pots; methional MG-ARP > Met/Glc at 100 and 120 C but Met/Glc > MG-ARP at 130 C
  ("less residual unreacted MG-ARP"); MGO is the predominant of the three dicarbonyls in both pots; GO
  and MGO rise then fall, faster at higher T; DA is low at 100 C and lower in MG-ARP than Met/Glc at
  120 and 130 C; Met in Met/Glc stays "considerably higher" than the Met released from MG-ARP; pH and
  browning move faster in the ARP pot (Fig. 1). No GO/MGO/DA/Met number is printed anywhere.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): methional -> `methional`; methylpyrazine ->
`methylpyrazine`; 2,5-/2,6-/2,3-dimethylpyrazine -> `2_5_dimethylpyrazine` / `2_6_dimethylpyrazine` /
`2_3_dimethylpyrazine`; ethylpyrazine -> `2_ethylpyrazine`; trimethylpyrazine ->
`trimethylpyrazine`; 2-ethyl-3,5-dimethylpyrazine -> `2_ethyl_3_5_dimethylpyrazine`; diacetyl ->
`2_3_butanedione` (alias "diacetyl"); **parent pyrazine -> not in registry as a molecule** (the class
row `pyrazines` carries the alias "pyrazine" meaning the family; same gap as Zhou 2024);
2-ethyl-5-methylpyrazine, vinylpyrazine, 2-vinyl-5-methylpyrazine, acetylpyrazine,
2-ethyl-3-(methylthio)pyrazine, methionine, glucose, glyoxal, methylglyoxal, MG-ARP -> not in registry.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| methional, Met + Glc, 30 / 60 / 120 / 180 / 240 min | 18.52 / 25.27 / 75.76 / 87.77 / 77.68 (= 0.178 / 0.243 / 0.727 / 0.843 / 0.746 µmol/L) | µg/L (dichlorobenzene-equivalent, RF = 1) | 200 + 200 mmol/L, water, initial pH 7.5 (NaOH), sealed vial, 120 C | none fitted | Table 1, p. 4 | level_only |
| methional, MG-ARP, same times | 25.42 / 54.14 / 196.59 / 153.85 / 85.60 (= 0.244 / 0.520 / 1.887 / 1.477 / 0.822 µmol/L) | µg/L (same) | 200 mmol/L ARP, same | none fitted | Table 1 | level_only |
| methional MG-ARP / Met+Glc at 30 / 60 / 120 / 180 / 240 min | 1.37 / 2.14 / 2.59 / 1.75 / 1.10 | — | 120 C | — | derived from Table 1 (same run, same method) | within_study_ratio |
| methional from Met alone | 1-2 | µg/L | 200 mmol/L Met, same conditions ("data not shown") | — | text 3.2, p. 5 | level_only |
| methional apparent formation rate, Met + Glc, 0-30 min (mine: 18.52/30) | 0.62 µg/L/min = 5.9e-3 µmol/L/min = 5.9e-6 mmol/L/min | µmol L-1 min-1 | 120 C, [Met] = [Glc] = 0.2 mol/L; assumes zero at t = 0 and no loss yet | zero-order secant | derived from Table 1 | derived_assumption |
| methional apparent formation rate, MG-ARP, 60-120 min (mine) | 2.37 µg/L/min = 0.0228 µmol/L/min | µmol L-1 min-1 | 120 C, 0.2 mol/L ARP; net of loss | zero-order secant | derived from Table 1 | derived_assumption |
| methional net loss, MG-ARP, 120-240 min (mine) | -0.925 µg/L/min; first-order equivalent ln(196.59/85.60)/120 = 6.9e-3 min-1 | min-1 | 120 C; a lower bound on the loss constant since formation continues | pseudo-first-order on the fall | derived from Table 1 | derived_assumption |
| parent pyrazine, Met + Glc, 120 C, 30-240 min | 4.96 -> 7.68 (0.062 -> 0.096 µmol/L) | µg/L | as above | none fitted | Table 1 | level_only |
| total pyrazines, Met + Glc / MG-ARP, 120 C | 14.32-22.03 / 5.95-11.40 | µg/L | as above | — | Table 1 | level_only |
| 2,5-dimethylpyrazine MG-ARP / Met+Glc, 120 C, 120 min | 5.70 / 3.83 = 1.49 | — | 120 C | — | Table 1 | within_study_ratio |
| parent pyrazine, Met + Glc, 130 C vs 120 C vs 100 C at 120 min | 14.43 / 7.16 / 1.88 | µg/L | same pot, three temperatures | — | Tables 1, 2 | level_only (a T-series of levels, usable as ratios 7.7 : 3.8 : 1) |
| 2,5-dimethylpyrazine, Met + Glc, 130 / 120 / 100 C at 120 min | 16.66 / 3.83 / n.d. | µg/L | same | — | Tables 1, 2 | level_only |
| methional at 100 and 130 C; total pyrazines vs time at 100/130 C | — | µg/L | — | — | Fig. 2 | figure_only |
| GO, MGO, DA time courses (both pots, three T) | — | mmol/L (OPD-HPLC, calibrated) | — | — | Fig. 3 | figure_only |
| Met consumption (Met/Glc) and Met release (MG-ARP) | — | mmol/L | — | — | Fig. 4 | figure_only |
| pH and A420 vs time, 120 C | — | — | — | — | Fig. 1 | figure_only |

**Within-corpus comparison for B18 (mine).** Zhou 2024 (same group, same GC method but external
calibration) fed [Ala] = [GO] = 20 mmol/L at pH 8, 120 C and got pyrazine at 0.1507 µmol/L/min, i.e.
about 18 µmol/L in 120 min. Here the whole Met + Glc pot at 200 + 200 mmol/L makes 0.09 µmol/L of
pyrazine in 120 min at 120 C, two hundred times less, with the glyoxal supplied by the sugar. This is
a response-factor-1 number against a calibrated one, so it is an order of magnitude, not a ratio; it
says the dicarbonyl supply in water, not the condensation, is what limits the pyrazine yield of a
sugar pot (the "what is missing" column of the roadmap's 2,5-dimethylpyrazine row).

## 5. Flags

1. **Semi-quantitative throughout.** Every volatile is peak-area against 0.09 µg of
   1,2-dichlorobenzene with response factor 1 on a DVB/CAR/PDMS fibre; no calibration, no LOD. The
   fibre's response to methional (polar, reactive) and to pyrazines against a chlorinated aromatic
   is unknown and can differ by several fold. Levels are order-of-magnitude; within-study ratios
   (ARP vs free Met, time ratios, temperature ratios) are the firmer content. Do not compare these
   µg/L with Pan 2025's or Zhou 2024's calibrated values as if they were the same quantity.
2. **Methional is printed only at 120 C** (Table 1). Table 2's title promises methional at 100 and
   130 C but the table has no methional row; those series exist only in Fig. 2a. Confirmed against
   the full text layer.
3. **The chain past methional is not measured.** Methanethiol, DMDS, DMTS and acrolein are named as
   the fate of methional (citing Yu & Ho 1995) but no sulfur volatile other than methional appears
   in the tables. The 6.9e-3 min-1 loss figure in section 4 is a net-decline number of mine, not a
   measured rate.
4. **Table 1 arithmetic.** MG-ARP 240 min "total minus pyrazine" printed 6.04, row sum 7.04 (the
   printed total 7.48 confirms 7.04). MG-ARP 60 min total printed 9.19, row sum 9.54. Keep the rows.
5. **Unbuffered, initial pH only.** pH 7.5 set with NaOH in water; the pH falls during heating
   (Fig. 1), faster in the ARP pot; the authors themselves attribute the ARP pot's low pyrazine yield
   to that fall. The working pH of every tabulated point is unknown. 0.2 mol/L reactants (ten times
   Zhou 2024).
6. **Sealed-vial volume, fill and headspace not stated**; the analysis takes 2 mL of product into a
   fresh 20-mL vial, so volatile loss on opening the reaction vial is not controlled and no salt is
   added for the SPME.
7. **Time zero is listed in the methods but absent from the tables**; the first tabulated point is
   30 min. The intercept needed for any rate fit is not available.
8. **Vessel heat-up time not stated** (relevant at 30 min).
9. **Minor inconsistency.** The text says "there was no methional-derived pyrazine observed" yet
   Table 2 lists 2-ethyl-3-(methylthio)pyrazine (1.72 µg/L, MG-ARP, 130 C, 180 min only, RI + MS
   identification).
10. **SI not on disk** (Table S1 OAVs, Fig. S5 sensory scores, ARP spectra).
11. **Registry gap**: parent pyrazine has no molecule row; five of the eleven pyrazines and all the
    precursors (methionine, glucose, GO, MGO, MG-ARP) are outside `compounds.yml`.
12. DL-methionine (racemic) was used; no consequence for Strecker chemistry, noted for the record.
