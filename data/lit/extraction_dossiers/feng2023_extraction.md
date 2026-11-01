# Feng, Yang, Xie, Yan, Ma, Hu & Yu 2024 — EXTRACTION (file stem `feng2023`)
DOI: 10.1002/jsfa.13408 (printed p. 5764; Crossref 2026-11-01: "The volatile organic compounds generated from the
Maillard reaction between L-ascorbic acid and L-cysteine in hot compressed water", J. Sci. Food Agric. 104(10),
5764-5775; online 2024-03-05, print 2024-08-15. **Matches the brief.**)
[!] **Year mismatch: the file is named `feng2023` but the paper is 2024** (received 8 Nov 2023, accepted 22 Feb 2024).
The stem is kept so the coverage gate matches. [!] **Possible wrong download:** the repo already cites a *different*
"Feng et al. (2023)" from the same group — Feng, Yang, Xie, Liu, Peng, Hu & Yu, "Kinetics of L-ascorbic acid
degradation and non-enzymatic browning development in hot-compressed water", Front. Nutr. 2023,
**10.3389/fnut.2022.1022254** (Crossref-checked 2026-11-01; `data/keys/papers.yml` paper_id
`frontiers_2022_hcw_aa_arrhenius_anchor`, dossier null; `data/lit/computational_priors.json` carries an
"Ea_hcw_ascorbic 82 kJ/mol" and Ea-by-pH from it, unverified). If `feng2023.pdf` was meant to be that paper, it is not;
that paper (open access, ASA loss kinetics vs T and pH) is the one with rate constants.

### L-ascorbic acid + L-cysteine, 2.45 mmol each (~70 mM [Z]) in 35 mL water, **unbuffered**, initial pH set with NaOH to 4.5 / 6.0 / 7.0 / 8.0 / 9.5, N2-purged sealed 316L autoclave (quartz liner), **170 and 210 °C × 30 min only**; HS-SPME GC-MS, 38 volatiles in µg/L against one internal standard (3-chloro-2,5-dimethylpyrazine, response factors not stated). **No time series, no final pH, no H2S, no residual Cys or ASA, no thiols.**

**Read 2026-11-01 from `data/articles/feng2023.pdf`** (12 pp., Wiley, born-digital). Prose and Table 1 from
`pdftotext -layout`; **Table 1 (rotated, pp. 5766-5770) checked against 130-200 dpi renders** (scratchpad
`ingest/feng2023_work/p3-03.png` .. `p7crop-07.png`); the text layer matches the renders cell for cell. Figs. 1-2
(p. 5771-5772) read from 170-200 dpi renders (`f8-08.png`, `f9-09.png`); Fig. 1 re-plots Table 1, Fig. 2 class totals
read to +-0.5 mg/L. Class sums of Table 1 computed by me (`feng2023_work/table1.tsv`) and tagged [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** School of Chemistry and Environmental Engineering, **Hubei Minzu University** (the former Hubei
University for Nationalities), Enshi; corresponding author **Yan Yang**, senior author **Ai-nong Yu**; Feng L also at
School of Light Industry, Beijing Technology and Business University. **Not in `data/keys/laboratories.yml`: NEW**,
proposed id `hubei_minzu_yu_an`. Evidence: same institution and senior author as `yu2016_extraction.md` (Yu & Tang
2016, ledger lab string "Hubei Univ. for Nationalities (Yu AN)", unmapped) and `yu2010_extraction.md` (Yu & Zhang
2010, Yu A-N lists both Central China Normal Univ. and Hubei Univ. for Nationalities; ledger string "Central China
Normal Univ. (Yu A-N / Zhang A-D)", unmapped). Same ASA + Cys pot family and the same GC programme as Yu 2010
(40 °C → 5 °C/min → 260 → 15 °C/min → 280 °C). Conservative choice for leave-one-laboratory-out: put all three under
one group. **Not** `btbu_xie` (the BTBU co-affiliation is Feng L at Light Industry; "Ya-ting Xie" is not Jianchun
Xie). **Not** the same group as `feng2021` (Feng T, Xia S, Jiangnan) or `feng2022` (Jiangnan Zhang X group, RG-ARP).

## 1. What this paper answers for the model

| model need | what Feng 2024 gives | tag |
|---|---|---|
| ascorbic-acid carbonyl lane | none in the engine; the paper gives no ASA loss, no dicarbonyl | — |
| H2S from cysteine vs T and pH | **not measured**; "H2S easier in acid" is the authors' inference | [P] |
| cysteine loss vs T / pH | **not measured** | — |
| polysulfide (H2S + aldehyde sink) vs pH | Σ6 polysulfides fall ×2.8 (170 °C) / ×2.0 (210 °C) from pH 4.5 to 9.5 | [Z] from [M] |
| thiophene formation vs pH, T | Σ16 thiophenes flat within ±20 % over pH 4.5-9.5; 210/170 ×1.0-1.3 | [Z] |
| pyrazine pH optimum | Σ5 pyrazines max pH 8.0 (170 °C) / pH 7.0 (210 °C), collapse at pH 9.5 (×6 drop at 210 °C) | [Z] |
| time course | **none** (30 min only, despite the abstract naming time) | [!] |
| absolute amounts | µg/L vs one IS, RF not stated → semi-quantitative | [M] / [!] |
| pH during / after heating | **not measured**; unbuffered pots at 170-210 °C | [!] |

**Answers.**
1. **Does anything transfer to the engine's sulfur lane? Very little** (confidence ~85 %). The engine has no ascorbic
   acid species; the carbon skeletons here come from ASA fragments, so no product row is simulable. What is left is
   directional: H2S-trapping products (trithiolanes, tetrathiocane, trithiolane) fall as the initial pH rises, in an
   unbuffered pot whose real pH is unknown.
2. **H2S vs pH: not answered** (~95 %). No sulfide was measured. The authors' "H2S more easily generated under acidic
   conditions" [Q/P] is an inference from polysulfide yields and **runs against** the engine's measured
   `k_cys_h2s` (Zheng & Ho 1994: neat cysteine thermolysis H2S rises x9-x17 from pH 5 to 9 at 80-110 °C, Table I as carried in `parameters_sulfur.ZHENG_TABLE_I_K_PER_MIN`). The
   polysulfide fall with pH has rival explanations (acetaldehyde / mercaptoacetaldehyde supply from ASA, H2S/HS-
   speciation and volatility, ring stability) that this paper cannot separate. Weight against the engine: negligible.
3. **Polysulfides vs pH: yes, directionally** (~70 % that the direction is real in this pot, given unknown in-situ pH).
   It agrees with Yu 2010 (same lab, 141 °C, phosphate): 3,5-dimethyl-1,2,4-trithiolanes and tetrathiocane fall from
   pH 6-7 to 8-9. It disagrees on 4,6-dimethyl-1,2,3-trithiane (Yu: rises pH 5 → 7; Feng: falls, nd at pH 8-9.5 at
   170 °C).
4. **Thiophenes from C4/C5 dicarbonyls: no rate, no precursor** (~95 %). Acyl thiophenes dominate (No. 10,
   5-methyl-2-acetylthiophene, up to 11.8 mg/L); the mechanism (dicarbonyl + H2S, Fig. 3B pathway 2) is drawn [P].
   The engine's SHET pool (norfuraneol + HS-) is a different precursor.
5. **Pyrazine pH optimum ~7-8 with a fall at 9.5** (~60 % real; two temperatures, unbuffered). Yu 2010 (phosphate)
   saw pyrazines rise monotonically to pH 9 — a disagreement at the alkaline end, plausibly buffer / in-situ pH.

## 2. Conditions

- **Reactants** [M]: "2.45 mmol ASA and Cys" in 35 mL double-distilled water [V] → 70 mM each if 2.45 mmol of each
  [Z; the 1:1 reading is inferred]. Analytical grade, Sigma-Aldrich.
- **Buffer / pH** [M]: none; pH adjusted to 4.5, 6.0, 7.0, 8.0, 9.5 with NaOH (Sartorius meter). Label = initial
  room-temperature pH. **No final pH, no pH during heating.** [!] ASA (pKa 4.1) + Cys at 70 mM would start near
  pH 2.5-3 [Z, not printed], so every arm carries a large, unstated NaOH load; at 170-210 °C acid formation will
  move the pH far from the label (cf. CF-054, CF-084).
- **Vessel** [M]: BE100-HP high-temperature high-pressure autoclave, 316L stainless steel with quartz lining
  (Shanghai LABE). Volume not printed ("BE100" suggests 100 mL, i.e. ~65 mL headspace [Z, unverified]).
- **Headspace / O2** [M]: N2 flowed through the inlet to expel air; reaction "oxygen-free" [Q].
- **Heating** [M]: aluminium-block furnace; 170 or 210 °C, 30 min, pressure = saturated vapour pressure
  (~0.8 / ~1.9 MPa [Z, steam tables]); quench in flowing cold water. **Heat-up time not reported** [!] — for a
  steel autoclave in a block it is likely a large share of 30 min.
- **Choice of T, t** [Q]: "based on preliminary experimental results" — not shown.
- **Replicates** [M]: "at least three times"; RSD "controlled within 10%" [V]; Table 1 ± values are RSD (footnote c),
  printed as fractions 0.014-0.065 (1.4-6.5 %).

## 3. Analytical method

- HS-SPME, DVB/CAR/PDMS 50/30 µm; extraction conditions **not printed** — "selected based on our previous paper"
  (ref. 28: Feng, Yang, Liu, Tan, Tan & Yu 2022, Food Chem 371:131155, ASA self-degradation volatiles in HCW; not on
  disk). NaCl in the materials list (salting presumably).
- GC-MS: 40 °C 4 min → 5 °C/min → 260 °C → 15 °C/min → 280 °C, 1 min; desorption 250 °C 2 min splitless.
  **Instrument and column not stated** [!].
- Identification: NIST 20 / Wiley 275, RI from n-C7-C22; co-injection (Co-GC) for 6 compounds (Nos. 1, 5, 6, 8, 16,
  20); 24 of 38 have no literature RI ("∼").
- **Quantification** [M]: from the concentration and peak area of the IS 3-chloro-2,5-dimethylpyrazine "and other
  standard substances reported in the literature" [V]. IS amount, addition point and response factors not stated →
  **semi-quantitative (IS, RF presumably 1)**. Unit µg/L (of reaction solution, presumably).
- No LOD; nd marked "–".
- Data availability statement points to an ORCID page, not a repository [!].

## 4. Data — Table 1 (pp. 5766-5770), µg/L, **170 °C / 210 °C**, mean of ≥ 3 (RSD 1.4-6.5 %) [M]

nd = not detected. Classes as in Fig. 1 (PZ pyrazines, TH thiophenes, PS polysulfides, OT others).

| No. | compound (RI calc/lit) | cl | pH 4.5 | pH 6.0 | pH 7.0 | pH 8.0 | pH 9.5 |
|---|---|---|---|---|---|---|---|
| 1 | 2-ethylpyrazine (915/914) | PZ | 759.09 / 597.19 | 1119.52 / 1305.60 | 1445.79 / 1724.78 | 1644.87 / 1807.03 | 368.48 / 123.13 |
| 2 | 2-methyl-5-ethylpyrazine (999/1000) | PZ | 470.77 / 422.05 | 618.09 / 1110.10 | 1076.45 / 1292.65 | 1206.77 / 1693.16 | 609.31 / 323.24 |
| 3 | 2,5-dimethyl-3-ethylpyrazine (1076/1074) | PZ | 1437.83 / 2183.10 | 840.33 / 1999.41 | 804.10 / 1591.82 | 787.41 / 568.36 | 754.42 / 462.23 |
| 4 | 6-methyl-2-ethylpyrazine (996/995) | PZ | 225.38 / 337.71 | 346.08 / 479.51 | 431.90 / 485.63 | 795.62 / 504.37 | 423.35 / nd |
| 5 | 2-methylpyrazine (825/826) | PZ | 5.70 / 21.67 | 18.31 / 466.31 | 326.20 / 876.86 | 373.58 / 920.82 | nd / nd |
| 6 | 2-acetyl-3-methylthiophene (1153/1161) | TH | 1402.98 / 2782.26 | 1135.30 / 2731.09 | 1006.20 / 2484.53 | 987.18 / 2361.62 | 817.73 / 1513.87 |
| 7 | 2-propanoylthiophene (1185/1188) | TH | 318.46 / 649.79 | 363.86 / 675.44 | 372.62 / 686.12 | 582.43 / 1164.55 | 1199.03 / 2111.39 |
| 8 | 3-acetylthiophene (1085/–) | TH | 1920.63 / 3368.52 | 1995.68 / 3582.15 | 2089.69 / 3709.41 | 3044.35 / 4323.99 | 5216.33 / 6009.72 |
| 9 | 2,5-thiophenedicarboxaldehyde (1194/1220) | TH | 2547.12 / 2030.55 | 1737.56 / 1686.38 | 1659.03 / 1519.81 | 1329.30 / 950.11 | 199.77 / 136.06 |
| 10 | 5-methyl-2-acetylthiophene (1208/–) | TH | 5938.77 / 11809.05 | 5488.27 / 10626.96 | 5255.54 / 8994.81 | 5147.60 / 6491.78 | 3793.99 / 3990.36 |
| 11 | 2,5-dimethyl-3-acetylthiophene (1280/–) | TH | 321.77 / 444.70 | 316.02 / 449.59 | 323.06 / 474.23 | 381.83 / 532.48 | 467.92 / 785.64 |
| 12 | 4-methyl-2-acetylthiophene (1188/–) | TH | 1891.98 / 2028.52 | 2424.16 / 2398.86 | 2461.67 / 2450.77 | 2595.07 / 3061.42 | 2719.87 / 3404.95 |
| 13 | 3-(methylsulfanyl)thiophene (1059/1082) | TH | 1146.78 / 410.67 | 1852.20 / 762.32 | 2041.72 / 899.74 | 3018.47 / 1381.03 | 2116.69 / 1350.34 |
| 14 | 2-(1-propenylthio)thiophene (1372/–) | TH | 547.85 / 2267.85 | 386.79 / 1614.43 | 352.58 / 1173.68 | 251.61 / 848.16 | 242.56 / 529.10 |
| 15 | 3-(vinylthio)thiophene (1309/–) | TH | 4433.26 / 2596.91 | 5576.85 / 3579.67 | 6306.05 / 3687.03 | 6601.28 / 4091.80 | 4133.03 / 2438.67 |
| 16 | 2-acetylthiophene (1090/1090) | TH | nd / 1005.97 | 1649.47 / 1446.29 | 1803.98 / 1511.50 | 2125.90 / 1860.68 | 2578.34 / 2374.62 |
| 17 | 2-methyl-5-propylthiophene (1218/–) | TH | 202.02 / nd | 242.92 / nd | 252.07 / 328.55 | 446.95 / 403.62 | 725.64 / 444.72 |
| 18 | 2-methylthieno[3,2-b]thiophene (1312/–) | TH | 1807.02 / 940.63 | 2170.71 / 1210.32 | 2392.28 / 1608.13 | 2601.56 / 1997.15 | 2783.04 / 2279.55 |
| 19 | thieno[2,3-b]thiophene (1249/1212) | TH | 4261.71 / 3679.81 | 4098.37 / 3459.97 | 3418.25 / 3316.48 | 3226.22 / 3009.05 | 3038.27 / 2680.92 |
| 33 | tetrahydro-3-methyl-2-propylthiophene (1124/–) | TH | 876.85 / 1526.04 | 818.49 / 1363.72 | 753.29 / 1164.29 | 517.18 / 772.03 | 494.76 / 703.48 |
| 38 | 2,3-dihydro-5-methylthiophene (905/–) | TH | nd / 428.95 | nd / 385.90 | nd / 381.52 | nd / 71.79 | nd / 52.55 |
| 22 | 3,5-dimethyl-1,2,4-trithiolane (1147/1148) | PS | 3937.48 / 5971.30 | 3280.99 / 5419.03 | 2880.60 / 5197.20 | 2395.15 / 4368.87 | 2032.21 / 4052.38 |
| 23 | 1,2,5-trithiepane (1253/–) | PS | 2314.49 / 1911.89 | 1952.63 / 1537.67 | 1875.25 / 1306.21 | 1811.83 / 1233.57 | 1633.40 / 1176.18 |
| 24 | 1,2,5,6-tetrathiacyclooctane (1394/–) | PS | 3611.83 / 8921.14 | 3584.09 / 7834.72 | 2782.35 / 7328.58 | 783.63 / 4545.95 | 333.92 / 3339.47 |
| 25 | 1,2,3-trithiolane (1534/–) | PS | 759.64 / 243.07 | 560.72 / 241.36 | 348.25 / 225.34 | 223.47 / 205.90 | 102.62 / 185.54 |
| 26 | 4-methyl-3,5-dithiaheptane (1320/–) | PS | 519.69 / 723.69 | 473.09 / 587.88 | 459.90 / 515.90 | 302.50 / 186.05 | 261.14 / 172.93 |
| 27 | 4,6-dimethyl-1,2,3-trithiane (1239/–) | PS | 890.22 / 1113.75 | 745.37 / 1051.54 | 694.74 / 1015.43 | nd / 981.00 | nd / 681.70 |
| 20 | 2-acetylthiazole (1019/1020) | OT | 574.71 / 1198.52 | 1141.33 / 1274.02 | 1368.80 / 1562.17 | 1754.52 / 1605.56 | 1323.64 / 917.73 |
| 21 | 2,4-dimethyl-4,5-dihydro-1,3-thiazole (953/–) | OT | 387.97 / 192.47 | 420.03 / 254.32 | 506.24 / 269.47 | 534.87 / 303.63 | 462.29 / 134.60 |
| 28 | 3,4-toluenedithiol (1364/–) | OT | 2618.50 / 954.98 | 2222.16 / 775.53 | 1943.40 / 635.92 | 1533.19 / 598.72 | 741.52 / 406.85 |
| 29 | 2-furancarbodithioic acid ethyl ester (1429/–) | OT | 2130.67 / 1546.09 | 1526.02 / 1296.04 | 1228.21 / 911.93 | 278.29 / 646.84 | 131.73 / 497.10 |
| 30 | 5-(methylthio)-2(1H)-pyridinone (1108/–) | OT | 418.31 / 322.12 | 1077.06 / 468.13 | 1381.28 / 473.84 | 2779.73 / 476.80 | 1567.41 / 458.52 |
| 31 | 3-methyl-4-(methylthio)phenol (1287/–) | OT | 763.83 / 601.26 | 788.50 / 719.01 | 854.96 / 865.61 | 870.73 / 890.72 | 953.01 / 923.98 |
| 32 | 2,3,5,6,7,8-hexahydro-1,4-benzodithiine (1441/–) | OT | 444.86 / 335.81 | 494.81 / 339.20 | 507.48 / 349.33 | 640.35 / 436.62 | 677.92 / 446.57 |
| 34 | 3-butyltetrahydro-2H-thiopyran (1232/–) | OT | 289.71 / 626.05 | 269.52 / 458.63 | 253.10 / 364.45 | 149.42 / 226.20 | 113.95 / 72.54 |
| 35 | 1-furfurylpyrrole (1181/1192) | OT | 3854.90 / 6122.46 | 1091.69 / 3457.60 | 886.95 / 2264.42 | 99.70 / 626.64 | nd / 515.68 |
| 36 | 4-(methylthio)thiophenol (1402/–) | OT | 1027.19 / 679.87 | 676.90 / 651.57 | 578.82 / 474.49 | 257.59 / 373.40 | nd / 335.89 |
| 37 | 2,4-dimethyl-3-ethyl-1H-pyrrole (1096/–) | OT | nd / 944.52 | nd / 890.97 | nd / 793.98 | nd / 784.65 | nd / 693.30 |

**Class sums of Table 1, µg/L** [Z] (170 / 210 °C):

| class (n) | pH 4.5 | pH 6.0 | pH 7.0 | pH 8.0 | pH 9.5 | pH 4.5 ÷ 9.5 |
|---|---|---|---|---|---|---|
| thiophenes (16) | 27 617 / 35 970 | 30 257 / 35 973 | 30 488 / 34 391 | 32 857 / 33 321 | 30 527 / 30 806 | 0.90 / 1.17 |
| polysulfides (6) | 12 033 / 18 885 | 10 597 / 16 672 | 9 041 / 15 589 | 5 517 / 11 521 | 4 363 / 9 608 | 2.76 / 1.97 |
| pyrazines (5) | 2 899 / 3 562 | 2 942 / 5 361 | 4 084 / 5 972 | 4 808 / 5 494 | 2 156 / 909 | 1.34 / 3.92 |
| others (11) | 12 511 / 13 524 | 9 708 / 10 585 | 9 509 / 8 966 | 8 898 / 6 970 | 5 971 / 5 403 | 2.10 / 2.50 |
| all 38 | 55 060 / 71 941 | 53 504 / 68 591 | 53 123 / 64 917 | 52 080 / 57 306 | 43 017 / 46 725 | — |

210/170 ratio [Z]: polysulfides 1.57 / 1.57 / 1.72 / 2.09 / 2.20 (pH 4.5 → 9.5); thiophenes 1.30 / 1.19 / 1.13 /
1.01 / 1.01; pyrazines 1.23 / 1.82 / 1.46 / 1.14 / 0.42.

**Single-compound pH contrasts (pH 4.5 ÷ 9.5)** [Z]: 3,5-dimethyl-1,2,4-trithiolane 1.94 (170) / 1.47 (210);
tetrathiocane 10.8 / 2.67; 1,2,3-trithiolane 7.4 / 1.31; 5-methyl-2-acetylthiophene 1.57 / 2.96; 2-acetylthiazole
0.43 / 1.31 (170 °C maximum at pH 8.0).

**Scale** [Z]: 3,5-dimethyl-1,2,4-trithiolane at 5 971 µg/L is 39 µM, i.e. ~5.6 × 10⁻⁴ mol per mol Cys (70 mM). All
38 volatiles together (≤ 72 mg/L) carry well under 1 % of the cysteine sulfur: these are trace sinks, not a mass
balance on H2S.

**Fig. 2** (p. 5772) [F]: class bars match the Table 1 sums for polysulfides and others (within ~0.2 mg/L) and
pyrazines (within ~0.9 mg/L); the thiophene bars run 3-10 % below the Table 1 sums (e.g. 30.3 vs 32.9 mg/L at pH 8.0 /
170 °C); the **"total" bars (39.3 → 31.8 mg/L at 170 °C, 53.4 → 37.3 at 210 °C) are ~15-17 mg/L below the sum of the
figure's own class bars** (e.g. 27.3 + 12.0 + 2.8 + 12.5 = 54.6 vs 39.3 at pH 4.5 / 170 °C). Fig. 1 re-plots Table 1.

## 5. Mechanism evidence and its kind

Within-study factor contrasts only (pH × 2 temperatures), endpoint, semi-quantitative. No labelling, no fed
intermediate, no omission arm (no ASA-alone or Cys-alone pot in this paper). Fig. 3 (ASA degradation to C5
polyhydroxy compounds and dicarbonyls; pyrazine, thiophene, polysulfide pathways) and Fig. 4 (precursor map) are [P].
The conclusion that Cys "did not participate in the reaction" except as N and S source [Q] rests on no labelling
experiment.

## 6. Bearing on open questions

- **Engine `k_cys_h2s` (Zheng & Ho 1994, H2S rising with pH 5 → 9)** — LEAVES UNCHANGED. The authors claim the
  opposite direction for H2S but measured no sulfide; their evidence is polysulfide yields, which also depend on
  aldehyde supply and in-situ pH. P(changes anything) ~2 %.
- **CF-012 (pH direction of MFT)** — NOT RELEVANT: no pentose, no MFT / FFT or any furanthiol reported (as in Yu
  2010). P ~0 %.
- **CF-054 / CF-084 (pH drift of unbuffered pots, engine charge balance)** — the paper is an example of the problem,
  not a datum: NaOH-adjusted, unbuffered, 170-210 °C, no final pH. P(helps) ~0 %.
- **SHET pool (norfuraneol + HS- → thiophenes, trithiolanes, lumped)** — qualitative support that, with H2S and C2-C5
  carbonyls, heterocyclic sulfur sinks dominate the volatile sulfur over thiols; no rate. P(changes) ~3 %.
- **Pyrazine pH optimum** — weak directional datum (max pH 7-8, fall at 9.5); disagrees with Yu 2010's monotone rise
  to pH 9 (phosphate). P(changes any pyrazine parameter) ~3 %.
- No other register entry on polysulfides, thiophenes, H2S or cysteine degradation vs pH was found
  (`grep -i "polysulf|trithiol|thiophene|ascorb"` of `docs/reference/conflicts_register.yml`).

## 7. Disagreements found

1. **vs Yu & Zhang 2010 (same lab family; ASA + Cys 100 mM, 0.2 M phosphate pH 5-9, 141 °C, 2 h, TIC areas)**:
   agree — dimethyltrithiolane and tetrathiocane fall to the alkaline side; acetylthiophenes (2- and 3-) rise with pH;
   2-acetylthiazole peaks near pH 8; thieno[2,3-b]thiophene falls with pH. Disagree — pyrazines (Yu: rise to pH 9;
   Feng: collapse at 9.5); 4,6-dimethyl-1,2,3-trithiane (Yu: rises 5 → 7; Feng: falls, nd at pH ≥ 8 / 170 °C);
   2-methylthieno[3,2-b]thiophene (Yu: falls; Feng: rises); thieno[3,2-b]thiophene, Yu's largest peak, is absent here.
   Candidate causes: buffer vs none (in-situ pH), 141 vs 170-210 °C, 2 h vs 30 min, identification (see flag 6).
2. **Authors' H2S-vs-pH claim vs Zheng & Ho 1994** (measured): see sec. 6; inference vs measurement.

## 8. Flags

1. [!] **Text contradicts Table 1 on pH optima.** Conclusion: Nos. 10, 22 and 24 "more easily generated" at pH 8.0;
   Table 1 has all three highest at pH 4.5 at both temperatures. Results: pyrazine No. 3 "readily generated" at pH 8.0;
   Table 1 has No. 3 highest at pH 4.5 (1437.83 / 2183.10) and near its lowest at pH 8.0.
2. [!] **"Only under neutral conditions"** is claimed for Nos. 4, 5, 16, 17, 27, 35, 36; Table 1 shows each detected
   at acidic and/or basic pH (e.g. No. 4 at all five pH at 170 °C).
3. [!] "Nos. 37 and 38 only identified at 201 °C" — Table 1 shows them at 210 °C only (typo). "Meldrum's reaction"
   for Maillard (p. 5765). Table title says "VPs".
4. [!] **Fig. 2 "total" ≠ sum of its class bars** by ~15-17 mg/L (sec. 4); thiophene bars 3-10 % below Table 1 sums.
5. [!] **Unbuffered, initial pH only**; large unstated NaOH load; no final pH; nominal pH 4.5-9.5 at 170-210 °C is not
   the reaction pH. Brief's "pH measured before/after": before only.
6. [!] Identifications: 32 of 38 by MS + RI only (Co-GC for 6), 24 with no literature RI; thieno[2,3-b]thiophene calc RI 1249 vs
   printed literature 1212 (Yu 2010 had it at 1192); 2,5-thiophenedicarboxaldehyde 1194 vs 1220. Yu 2010's
   "4-methyl-1,2-benzenedithiol" (RI 1358) is the same compound as No. 28 "3,4-toluenedithiol" (1364).
7. [!] Semi-quantitative: one IS (a chloropyrazine) for thiophenes, polysulfides and pyrroles; RF, IS amount and SPME
   conditions not printed; column and instrument not stated.
8. [!] **No time series** although the abstract names time as a studied parameter; T and t "chosen" from unshown
   preliminary runs; heat-up of the autoclave not reported.
9. [!] RSD "controlled within 10%" suggests repeats were selected to meet the target; ± values are RSD fractions, not
   SDs.
10. [!] Concentration basis inferred (2.45 mmol "ASA and Cys" read as each); autoclave volume not printed.

## 9. Ledger-ready measurements

- Table 1 sulfur volatiles (16 thiophenes, 6 polysulfides, 2 thiazoles, aromatic thiols / thioethers) × 5 pH × 2 T,
  endpoint 30 min — sulfur ledger, semi-quantitative, not simulable (no ascorbic-acid lane).
- Table 1 pyrazines (5) and pyrroles (2) × 5 pH × 2 T — aroma_other ledger, semi-quantitative, not simulable.
- Not keyable: any rate constant, Ea, H2S, cysteine or ASA loss, pH trajectory (none measured).
- Lead for the download queue: Feng et al. 2023 Front. Nutr. 10.3389/fnut.2022.1022254 (ASA degradation kinetics in
  HCW, same lab; the repo's "Feng et al. (2023)" anchor), and Feng et al. 2022 Food Chem 371:131155 (ASA
  self-degradation volatiles in HCW, the method reference).
