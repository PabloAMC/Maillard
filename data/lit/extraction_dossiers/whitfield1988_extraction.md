# Whitfield et al. 1988 — EXTRACTION (glycine, lysine or cysteine + ribose, with and without egg lecithin, 0.2 M phosphate pH 5.7 / 6.2, 140 C / 1 h; 68 heterocycles in ng per ampoule by Tenax headspace, n = 3)
### The first with/without-phospholipid Maillard pot: 2-pentylpyridine and four 2-alkylthiophenes appear only with lecithin, and MFT, FFT and the thiophenethiols fall 2-3.4 fold when it is present.

**Source on disk:** `data/articles/whitfield1988.pdf` (owner's download, 2026-09-08). Scanned
journal pages with an OCR text layer. The running text reads cleanly; Table 1 (pp. 266-267) is
scrambled in the plain text layer, so the two pages were re-extracted with `pypdf` in layout mode
and every cell was checked against a rendering of the pages; the re-typed table below is from that
check. No figures in the paper.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Phospholipid on the Formation of Volatile Heterocyclic Compounds in Heated Aqueous Solutions of Amino Acids and Ribose" |
| Authors | Frank B. Whitfield (CSIRO, on leave), Donald S. Mottram, Susan Brock, David J. Puckey, Linda J. Salter (AFRC Institute of Food Research, Bristol) |
| Venue | J. Sci. Food Agric. 42 (1988) 261-272; received 12 March 1987, accepted 25 June 1987 |
| DOI | not printed on the scan |
| Companions | Salter, Mottram, Whitfield 1988 (JSFA 46, 227-242: glycine + ribose + PC identifications); Farmer, Mottram, Whitfield 1989 (JSFA 49, 347-368: cysteine + ribose + PC identifications); Farmer & Mottram 1990 (`farmer1990_extraction.md`: the same cysteine pot with four lipids, relative areas) |
| Correction to carry | the row printed as "2-heptylthiophene" (622.2 ng) was re-identified as **2-pentylthiapyran** by Farmer & Mottram 1990 ("in our previous paper this compound was incorrectly identified as 2-heptylthiophene") |

## 1. Why it matters

The sulfur lane's reference pot is cysteine + ribose without lipid; the isolates carry 1-3 % lipid.
This is the only paper on disk that heats cysteine + ribose with and without a phospholipid and
prints absolute amounts (ng, external standards, n = 3) for MFT, FFT, the thiophenethiols, the
alkylthiazoles and the lipid–Maillard products side by side. The within-study ratios (no lecithin /
lecithin) are: MFT 3.0, FFT 2.0, 2-thiophenethiol 3.4, 3-thiophenethiol 2.1,
2-methyltetrahydrothiophen-3-one 2.1, 2-methylthiazole 3.0, 4,5-dimethylthiazole 1.7, methylpyrazine
2.1, the other pyrazines 1.0-1.5. The products that need lipid appear only in the lecithin pots:
2-pentylpyridine 194.7 ng with cysteine (9.0 with glycine, 1.9 with lysine); 2-butyl-, 2-pentyl-,
2-hexylthiophene 20.8 / 95.5 / 81.5 ng and the "2-heptylthiophene" (= 2-pentylthiapyran) 622.2 ng
only with cysteine + lecithin; the 2-alkylfurans in every lecithin pot, 2.6-6.4x more with an amino
acid + ribose present than in lecithin alone. The glycine and lysine columns are the same experiment
for the sugar lane (furfural 3.2x lower with lecithin in the glycine pot).

## 2. Methods as they matter to a model

- **Charge (per ampoule):** 2 mL of a stock of 5 g/L amino acid + 4.5 g/L D-ribose in 0.2 M sodium
  phosphate buffer; i.e. 10 mg amino acid + 9 mg ribose per ampoule. In molar terms (mine):
  cysteine 41.3 mmol/L, glycine 66.6 mmol/L, L-lysine 34.2 mmol/L if free base (the salt form is not
  stated), ribose 30.0 mmol/L; cysteine : ribose = 1.38 : 1. "Selected to approximate their relative
  compositions in mammalian muscle."
- **Lipid:** 300 µL of egg L-alpha-phosphatidylcholine (Sigma type III-E) in hexane at 100 mg/mL,
  evaporated under helium in the ampoule before the aqueous charge = **30 mg lecithin per 2 mL =
  15 g/L** (about 19 mmol/L PC, about 39 mmol/L acyl chains, my estimate from an average PC mass of
  770). Lecithin melted and sonicated 10 min to disperse. Fatty-acid composition of this PC is not
  given here; Farmer 1990 Table 1 gives it for the same Sigma product (16.2 % 18:2, 3.8 % 20:4).
- **Seven mixtures x 3 replicates = 21 ampoules:** each amino acid + ribose with and without
  lecithin, plus lecithin + buffer alone.
- **pH:** glycine and cysteine pots 5.7 before reaction, lysine 6.2; **fell by about one unit
  during heating**, the same with and without lecithin.
- **Vessel and heating:** 10 mL round-bottom Pyrex ampoules (flamed), purged with helium 2 min,
  flame-sealed (about 8 mL headspace, no oxygen added), laid horizontal in a Certoclav autoclave,
  **140 C, 0.28 MPa, 1 h**.
- **Volatile collection:** the cooled 2 mL transferred to a 250 mL flask with 18 mL phosphate
  buffer pH 5.7 (10x dilution), stirred at 60 C, purged with oxygen-free N2 at 60 mL/min for 30 min
  onto Tenax GC (trap at 20 C), then dried 5 min. Not exhaustive: the ng values are amounts trapped,
  not amounts in the liquid (flag 1).
- **GC-MS:** thermal desorption 250 C / 15 min onto a CPWAX 57CB 50 m x 0.32 mm column with
  cryofocusing; 60 C (5 min) -> 200 C at 4 C/min -> 200 C (10 min); Finnigan 4000, EI 40 eV, m/z
  33-400, 1 scan/s; LRI from C10-C20 alkanes added to the trap (5.5 mg/L, 1 µL, i.e. 5.5 ng each).
- **Quantification:** reverse-search against a library of the 68 target spectra in a narrow
  retention window; one characteristic ion per compound (the molecular ion except ethylpyrazine m/z
  107, pentylpyridine m/z 93, hexylthiophene m/z 97, heptylthiophene m/z 97); ion area compared with
  that of the authentic compound (40 ng) loaded on a Tenax trap and run the same way. **External
  standard placed on the trap, not in the liquid**, so recovery from the liquid is not corrected.
  Where no authentic compound existed, response estimated from an isomer or a close relative
  (rows marked c, tentative). Coefficients of variation 3-30 % across the triplicates; per-row SDs
  are not printed.

## 3. Tables re-typed

### Table 1. "Quantities (ng) of Selected Heterocyclic Components Formed in Maillard Reaction Mixtures in Both the Absence and Presence of Lecithin." Mean of three replicate reaction mixtures; ng trapped from one 2 mL reaction mixture. "—" = not found; blank = ratio not printed (one side not found). Ratio = without / with lecithin, as printed. LRI: reaction product / reference compound (CPWAX 57CB). c = tentative identification.

| compound | LRI rxn | LRI ref | Gly − | Gly + | Gly −/+ | Lys − | Lys + | Lys −/+ | Cys − | Cys + | Cys −/+ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Pyrazines** | | | | | | | | | | | |
| methylpyrazine | 1271 | 1272 | 23.8 | 7.2 | 3.3 | 534.7 | 317.9 | 1.7 | 27.3 | 12.9 | 2.1 |
| 2,5-dimethylpyrazine | 1329 | 1331 | 70.9 | 48.2 | 1.5 | 179.1 | 198.2 | 0.9 | 4.7 | 3.5 | 1.3 |
| ethylpyrazine | 1334 | 1336 | 1.8 | 1.3 | 1.4 | 12.8 | 23.0 | 0.6 | 9.0 | 8.1 | 1.1 |
| 2,6-dimethylpyrazine | 1334 | 1336 | 4.8 | 2.8 | 1.7 | 30.0 | 16.1 | 1.9 | 1.6 | 1.6 | 1.0 |
| 2,3-dimethylpyrazine | 1352 | 1354 | 3.3 | 1.9 | 1.7 | 10.5 | 7.1 | 1.5 | 6.7 | 5.3 | 1.3 |
| 2-ethyl-6-methylpyrazine | 1386 | 1387 | 1.6 | 0.9 | 1.8 | 13.0 | 8.4 | 1.5 | 0.9 | 0.9 | 1.0 |
| 2-ethyl-5-methylpyrazine | 1394 | 1393 | 41.8 | 23.7 | 1.8 | 107.5 | 98.0 | 1.1 | 2.1 | 1.9 | 1.1 |
| 2-ethyl-3-methylpyrazine | 1407 | 1407 | 0.2 | 0.1 | 2.0 | 1.2 | 1.7 | 0.7 | 2.3 | 2.0 | 1.1 |
| trimethylpyrazine | 1412 | 1412 | 17.8 | 13.9 | 1.3 | 10.8 | 9.3 | 1.2 | 1.1 | 0.8 | 1.4 |
| 2-ethyl-3,6-dimethylpyrazine | 1450 | 1449 | 0.9 | 0.8 | 1.1 | 6.1 | 5.5 | 1.1 | 0.4 | 0.2 | 2.0 |
| 2,5-diethylpyrazine | 1460 | 1463 | 3.3 | 2.7 | 1.2 | 7.8 | 7.2 | 1.1 | 0.3 | 0.2 | 1.5 |
| 2-ethyl-3,5-dimethylpyrazine | 1466 | 1466 | 9.3 | 10.0 | 0.9 | 4.0 | 4.0 | 1.0 | 0.5 | 0.5 | 1.0 |
| tetramethylpyrazine | 1485 | 1484 | 5.3 | 4.7 | 1.1 | 0.3 | 0.4 | 0.8 | 0.3 | — | |
| 3,5-diethyl-2-methylpyrazine c | 1496 | | 0.3 | 0.5 | 0.6 | 1.6 | 2.3 | 0.7 | 0.2 | 0.1 | 2.0 |
| **Pyridines** | | | | | | | | | | | |
| pyridine | 1209 | 1211 | 0.3 | — | | 13.1 | 6.5 | 2.0 | 0.1 | — | |
| 2-methylpyridine | 1242 | 1242 | 1.3 | — | | 1.4 | 0.8 | 1.8 | 0.9 | 0.9 | 1.0 |
| 4-methylpyridine | 1320 | 1322 | 0.5 | — | | 3.6 | 2.9 | 1.2 | — | — | |
| **2-pentylpyridine** | 1583 | 1584 | — | **9.0** | | — | **1.9** | | — | **194.7** | |
| 2-ethyl-4-methylpyridine c | 1587 | | — | — | | 154.7 | 344.3 | 0.4 | — | — | |
| an ethyldimethylpyridine c | 1662 | | — | — | | 1.3 | 3.6 | 0.4 | — | — | |
| an ethyldimethylpyridine c | 1732 | | — | — | | 28.5 | 113.6 | 0.2 | — | — | |
| **Furans** | | | | | | | | | | | |
| 2-butylfuran | 1119 | 1122 | 1.3 | 12.8 | 0.1 | 0.6 | 4.1 | 0.1 | 0.9 | 17.1 | 0.05 |
| 2-pentylfuran | 1219 | 1219 | 1.7 | 1425.2 | 0.001 | 1.2 | 674.1 | 0.002 | 1.2 | 704.7 | 0.002 |
| **2-methyl-3-furanthiol** | 1292 | 1295 | — | — | | — | — | | **2595.8** | **872.6** | **3.0** |
| 2-hexylfuran | 1312 | 1316 | 0.11 | 32.5 | 0.004 | 0.1 | 6.2 | 0.01 | — | 15.6 | |
| 2-methyl-3-(methylthio)furan | 1329 | 1333 | — | — | | — | — | | 1.7 | 2.4 | 0.7 |
| 2-heptylfuran | 1417 | 1416 | 0.17 | 82.1 | 0.002 | 0.1 | 9.9 | 0.01 | — | 67.6 | |
| **2-furanmethanethiol** | 1417 | 1421 | — | — | | — | — | | **2781.0** | **1384.7** | **2.0** |
| 2-furfural | 1454 | 1455 | 1063.2 | 328.7 | 3.2 | 7.4 | 5.4 | 1.4 | 41.3 | 24.2 | 1.7 |
| 3- or 4-methyl-2-furfural c | 1480 | | 7.0 | 2.4 | 2.9 | 1.3 | 0.5 | 2.6 | 0.2 | — | |
| benzofuran | 1488 | 1489 | 5.1 | 3.9 | 1.3 | 1.3 | — | | 0.2 | 0.2 | 1.0 |
| 2-acetylfuran | 1499 | 1499 | 33.0 | 25.7 | 1.3 | 6.8 | 8.9 | 0.8 | 4.3 | 4.6 | 0.9 |
| 2-octylfuran | 1520 | 1519 | 3.2 | 174.1 | 0.02 | 0.8 | 21.3 | 0.04 | — | 85.7 | |
| 3- or 4-methyl-2-furfural c | 1549 | | 1.7 | 1.6 | 1.1 | 0.6 | 1.0 | 0.6 | — | — | |
| 1-(2-furfuryl)pyrrole | 1823 | 1821 | 0.7 | 2.3 | 0.3 | 2.8 | 2.8 | 1.0 | 1.6 | 2.1 | 0.8 |
| **Thiophenes** (glycine and lysine pots: — throughout) | | | | | | | | | | | |
| 2-methylthiophene | 1092 | 1089 | — | — | | — | — | | 6.6 | 3.1 | 2.1 |
| 4,5-dihydro-2-methylthiophene c | 1142 | | — | — | | — | — | | 4.0 | 1.8 | 2.2 |
| 2,5-dimethylthiophene | 1149 | 1153 | — | — | | — | — | | 7.3 | 4.2 | 1.7 |
| 2-ethylthiophene | 1162 | 1161 | — | — | | — | — | | 2.5 | 1.4 | 1.8 |
| 2,4-dimethylthiophene c | 1181 | | — | — | | — | — | | 0.7 | 0.2 | 3.5 |
| 2,3-dimethylthiophene | 1203 | 1207 | — | — | | — | — | | 3.1 | 1.7 | 1.8 |
| 2-propylthiophene | 1238 | 1239 | — | — | | — | — | | 0.3 | 3.3 | 0.09 |
| trimethylthiophene c | 1265 | | — | — | | — | — | | 2.4 | 1.4 | 1.7 |
| **2-butylthiophene** | 1336 | 1339 | — | — | | — | — | | — | **20.8** | |
| **2-pentylthiophene** | 1437 | 1440 | — | — | | — | — | | — | **95.5** | |
| 2-methyltetrahydrothiophen-3-one | 1520 | 1518 | — | — | | — | — | | 574.5 | 272.1 | 2.1 |
| **2-hexylthiophene** | 1544 | 1545 | — | — | | — | — | | — | **81.5** | |
| **2-thiophenethiol** | 1561 | 1559 | — | — | | — | — | | **478.2** | **140.7** | **3.4** |
| 3-thiophenethiol c | 1581 | | — | — | | — | — | | 16.7 | 7.8 | 2.1 |
| **"2-heptylthiophene" (= 2-pentylthiapyran, Farmer 1990)** | 1660 | 1653 | — | — | | — | — | | — | **622.2** | |
| 2-thiophenemethanethiol | 1686 | 1689 | — | — | | — | — | | 1.6 | 0.9 | 1.8 |
| 2-formylthiophene | 1689 | 1688 | — | — | | — | — | | 31.8 | 26.1 | 1.2 |
| 2-formyl-4-methylthiophene c | 1709 | | — | — | | — | — | | 2.2 | 1.2 | 1.8 |
| 2-acetyl-3-methylthiophene | 1761 | 1760 | — | — | | — | — | | 10.9 | 10.1 | 1.1 |
| 3-acetylthiophene | 1771 | 1772 | — | — | | — | — | | 9.0 | 8.2 | 1.1 |
| 2-acetylthiophene | 1777 | 1777 | — | — | | — | — | | 4.6 | 5.4 | 0.9 |
| 2-formyl-5-methylthiophene | 1781 | 1780 | — | — | | — | — | | 3.9 | 3.1 | 1.3 |
| 2-formyl-3-methylthiophene | 1814 | 1813 | — | — | | — | — | | 85.5 | 60.1 | 1.4 |
| 2-propionylthiophene | 1840 | 1842 | — | — | | — | — | | 21.0 | 18.8 | 1.1 |
| thieno[3,2-b]thiophene | 1872 | 1876 | — | — | | — | — | | 527.1 | 344.9 | 1.5 |
| thieno[3,4-b]thiophene c | 1975 | | — | — | | — | — | | 24.8 | 22.2 | 1.1 |
| **Thiazoles** (glycine and lysine pots: — throughout) | | | | | | | | | | | |
| 2-methylthiazole c | 1243 | | — | — | | — | — | | 6.7 | 2.2 | 3.0 |
| 4,5-dihydro-2-methylthiazole c | 1312 | | — | — | | — | — | | 1.2 | 0.4 | 3.0 |
| 4,5-dimethylthiazole | 1377 | 1378 | — | — | | — | — | | 3.8 | 2.2 | 1.7 |
| trimethylthiazole | 1387 | 1385 | — | — | | — | — | | 17.2 | 16.5 | 1.1 |
| 5-ethyl-4-methylthiazole | 1440 | 1446 | — | — | | — | — | | 3.5 | 2.4 | 1.5 |
| 5-ethyl-2,4-dimethylthiazole | 1450 | 1455 | — | — | | — | — | | 10.4 | 7.8 | 1.3 |
| 2-acetylthiazole | 1640 | 1639 | — | — | | — | — | | 14.8 | 25.2 | 0.6 |

Footnote b (lecithin + buffer alone, no amino acid, no ribose; ng): 2-methylpyridine 0.5;
2-butylfuran 3.2; 2-pentylfuran 248; 2-hexylfuran 1.6; 2-heptylfuran 4.6; 2-octylfuran 15.3. All
other nominated compounds not detected in measurable amounts.

Notes on the re-typing: (a) the layout extraction misaligned the thiazole rows by one line; the
values above follow the rendering (4,5-dimethylthiazole 3.8 / 2.2 / 1.7, trimethylthiazole 17.2 /
16.5 / 1.1). (b) 2-octylfuran glycine "+" = 174.1 is legible only on the rendering. (c) Ratios are
as printed; they are the ratio of the means, rounded.

### Derived within-study numbers (mine, from Table 1)

| quantity | value |
|---|---|
| total 2-alkylfurans, cysteine + ribose + lecithin | 17.1 + 704.7 + 15.6 + 67.6 + 85.7 = 890.7 ng; lecithin alone 272.7 ng; ratio 3.27 (paper: 3.3; glycine 6.4, lysine 2.6) |
| 2-pentylpyridine, cysteine / glycine / lysine | 194.7 / 9.0 / 1.9 = 21.6 and 102 (paper: "some 20- and 100-fold") |
| "2-heptylthiophene" : 2-pentylthiophene (cysteine + lecithin) | 622.2 / 95.5 = 6.5 (paper: 7:1); 2-heptylfuran : 2-pentylfuran = 67.6 / 704.7 = 0.096 (paper: 1:10) |
| MFT : FFT without lecithin; with | 0.93; 0.63 |
| MFT in the liquid if the purge were exhaustive | 2595.8 ng / 2 mL = 1.3 mg/L (a lower bound; compare Hofmann 1998 Table 1: 19.8 µg / 100 mL = 0.2 mg/L at 33 mM Cys, 100 mM ribose, pH 5, 145 C, 20 min, SIDA) |

## 4. Routes and numbers the repository can use

| route or quantity | reactant -> product | mechanism as drawn | measured numbers, units, conditions | evidence class |
|---|---|---|---|---|
| W88-A lipid quench of MFT | cysteine + ribose (+ lecithin) -> MFT | none drawn; "competition for hydrogen sulphide ... by other reaction products" | 2595.8 -> 872.6 ng per 2 mL, ratio 3.0; 41 mM Cys, 30 mM ribose, 15 g/L PC, pH 5.7 -> about 4.7, 140 C, 1 h | within_study_ratio |
| W88-B lipid quench of FFT | -> FFT | as A | 2781.0 -> 1384.7 ng, ratio 2.0 | within_study_ratio |
| W88-C lipid quench of thiophenethiols | -> 2-thiophenethiol; 3-thiophenethiol | as A | 478.2 -> 140.7 (3.4); 16.7 -> 7.8 (2.1) | within_study_ratio |
| W88-D lipid quench of 2-methyltetrahydrothiophen-3-one | | as A | 574.5 -> 272.1 (2.1) | within_study_ratio |
| W88-E thiazoles under lipid | 2-methylthiazole; 4,5-dihydro-2-methylthiazole; 4,5-dimethyl; trimethyl; 5-ethyl-4-methyl; 5-ethyl-2,4-dimethyl; 2-acetylthiazole | as A ("aldehydes can also react with ammonia") | ratios 3.0; 3.0; 1.7; 1.1; 1.5; 1.3; 0.6 (2-acetylthiazole rises: 14.8 -> 25.2) | within_study_ratio |
| W88-F pyrazines under lipid (all three amino acids) | | competition of lipid aldehydes for NH3 / amino N; "a trend rather than a major effect" | glycine 0.6-3.3 (median about 1.5); lysine 0.6-2.0 (median about 1.1); cysteine 1.0-2.1; methylpyrazine the most reduced in all three (3.3, 1.7, 2.1) | within_study_ratio |
| W88-G furfural under lipid | ribose -> furfural | "condensation with lipid oxidation products ... to form non-volatile polymeric materials" (stated, not drawn) | glycine 1063.2 -> 328.7 (3.2); lysine 7.4 -> 5.4 (1.4); cysteine 41.3 -> 24.2 (1.7) | within_study_ratio |
| W88-H 2-pentylpyridine | 2,4-decadienal (from lecithin) + NH3 (from cysteine Strecker) -> 2-pentylpyridine | stated, citing Buttery 1977; drawn later in Farmer 1990 Fig 1 | only with lecithin: 194.7 ng (Cys), 9.0 (Gly), 1.9 (Lys); zero in lecithin alone | measured_yield (ng, external standard on trap) |
| W88-I 2-alkylthiophenes | lipid aldehydes + H2S -> 2-butyl-, 2-pentyl-, 2-hexylthiophene (and the row later re-assigned to 2-pentylthiapyran) | two candidates discussed: H2S on the 2-alkylfuran (rejected by the ratio argument: thiophene C7:C5 = 7:1 versus furan 1:10) and H2S on 4-oxoaldehydes or unsaturated aldehydes (favoured) | only with cysteine + lecithin: 20.8, 95.5, 81.5, ("heptyl") 622.2 ng; 2-propylthiophene 0.3 -> 3.3 | measured_yield + within_study_ratio |
| W88-J 2-alkylfurans | lecithin autoxidation -> 2-butyl- to 2-octylfuran | "known autoxidation product of linoleic acid" | lecithin alone 3.2 / 248 / 1.6 / 4.6 / 15.3 ng; with Cys + ribose 17.1 / 704.7 / 15.6 / 67.6 / 85.7; with Gly + ribose 12.8 / 1425.2 / 32.5 / 82.1 / 174.1; with Lys + ribose 4.1 / 674.1 / 6.2 / 9.9 / 21.3. The Maillard pot raises furan formation 2.6-6.4x ("unexpected, as Maillard reaction products are reported to act as antioxidants") | level_only |
| W88-K 2-methyl-3-(methylthio)furan | MFT + methanethiol source | not discussed | 1.7 -> 2.4 (0.7); the only MFT derivative that does not fall | level_only |

## 5. Rule sketches (reactant -> product in words; controls are mine, SMILES mine)

**S1. 2,4-decadienal + NH3 -> 2-pentylpyridine (W88-H; stated here, drawn in Farmer 1990 Fig 1:
carbinolamine at C1, N attacks C5, loss of water, oxidation; terminal).**
- positive: `CCCCC/C=C/C=C/C=O` (`DECADIENAL`) + `N` -> `CCCCCc1ccccn1`.
- negative: hexanal `CCCCCC=O` + `N` -> no pyridine (gives the aldimine / trialkylpyridine of Elmore 1997 instead); 2-pentylfuran `CCCCCc1ccco1` + `N` -> no fire.
- The NH3 comes from cysteine's Strecker branch; the paper's glycine and lysine numbers (9.0 and 1.9 ng versus 194.7) are the evidence that free NH3, not the amino acid, is the N donor. A rule keyed on `Cys`-derived NH3 would reproduce the ranking; one keyed on any amine would not.

**S2. lipid dienal + H2S -> 2-alkylthiophene (W88-I, the favoured alternative; see
`farmer1990_extraction.md` S2 and `mottram2002b_extraction.md` S2 for the drawn version).**
- positive: 2,4-nonadienal `CCCC/C=C/C=C/C=O` + `S` -> 2-pentylthiophene `CCCCCc1cccs1`.
- negative: 2-pentylfuran `CCCCCc1ccco1` + `S` -> the furan -> thiophene exchange must not be the rule (this paper's ratio argument); MFT + `S` -> no fire.

**S3. quench (a sink, not a structural rule).** The 2-3.4 fold fall in MFT, FFT and the
thiophenethiols with 15 g/L PC is the number a lipid-carrying isolate recipe should reproduce. The
paper assigns it to lipid carbonyls consuming H2S and NH3; a sink term on H2S proportional to the
unsaturated-aldehyde supply is the simplest encoding. Positive control for the sink: `DECADIENAL` +
`S` consumes H2S (products per S2 or the thiapyran rule); negative: methyl stearate + `S` -> no fire.

## 6. Flags

1. **External standard on the trap, non-exhaustive purge, no per-row SD.** The ng values are
   amounts trapped in 30 min from a 10x-diluted 2 mL sample; recovery from the liquid is not
   corrected, so absolute values are lower bounds. The −/+ ratios cancel this only if lecithin does
   not change the headspace partition; **lecithin at 15 g/L can bind hydrophobic volatiles**, which
   would lower the "+" column independently of chemistry. The paper does not test this (Farmer 1990
   raises the same point for triglyceride). Treat the thiol ratios as an upper bound on the chemical
   quench.
2. **The "2-heptylthiophene" row (622.2 ng) is 2-pentylthiapyran** (Farmer 1990). The quantification
   ion (m/z 97) and the response of "2-heptylthiophene" were used, so the 622.2 is a
   2-heptylthiophene-equivalent, not a calibrated thiapyran amount.
3. **pH drifted about one unit** (5.7 -> about 4.7; 6.2 -> about 5.2) in 0.2 M buffer; Farmer 1990
   moved to 0.5 M for this reason. The two studies are therefore not at the same pH history.
4. **Helium-purged sealed ampoules, yet lecithin oxidised to C4-C8 alkylfurans, dienals and
   alkanals**: the oxygen came from dissolved air, from the lipid's pre-existing hydroperoxides
   (Sigma egg PC in hexane) or from thermal (non-radical) routes; not resolved. The lipid lane's
   hydroperoxide charge for this pot is unknown.
5. **Trace 2-alkylfurans in the lipid-free pots** (1-2 ng) are "not explained"; the authors exclude
   cross-contamination.
6. **CV 3-30 %** is the only error statement; ratios near 1.0-1.5 (most pyrazines, most acylthiophenes,
   thiazoles other than 2-methyl) are within noise.
7. **Tentative rows (c)** were quantified with a relative's response; among the sulfur rows these
   are 4,5-dihydro-2-methylthiophene, 2,4-dimethylthiophene, trimethylthiophene, 3-thiophenethiol,
   2-formyl-4-methylthiophene, thieno[3,4-b]thiophene, 2-methylthiazole, 4,5-dihydro-2-methylthiazole.
8. **Lysine form not stated** (free base or hydrochloride), so its molarity is 34 or 27 mmol/L.
9. **Charge ratio differs from Hofmann 1998** (cysteine in excess of ribose here, 1.38:1; Hofmann
   1:3) and time is 3x longer at 5 C lower; the MFT lower bound of 1.3 mg/L versus Hofmann's 0.2 mg/L
   is consistent with excess cysteine and the longer time but cannot be used quantitatively.
10. **No lipid-only + cysteine control** (lipid + amino acid without ribose) and no time course.
