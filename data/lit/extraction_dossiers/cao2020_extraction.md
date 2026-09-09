# Cao et al. 2020 — EXTRACTION (olive and camellia oil heated in an open stirred fryer at 120 / 150 / 180 C for 24 h; DNPH-HPLC carbonyls; oleate 8-/9-/10-/11-OOH -> octanal / nonanal / decanal / 2-decenal / 2-undecenal assignment)
### The oleic-acid triglyceride aldehyde map at frying temperatures, with real concentrations (nmol/g) for the five oleate aldehydes.

**Source on disk:** `data/articles/cao2020.pdf` (owner's download, 2026-09-08; "Article in press" proof,
LWT 2020, article 108858). Read from the `pypdf` text layer; Tables 1 and 2 are clean and re-typed
below. Figures 1-5 are time-course plots: FIGURE-ONLY. Figure 6 is the drawn mechanism (page 8),
described in words from a rendering of the page.

## 0. Identity

| field | value |
|---|---|
| Title | "Oxidative stabilities of olive and camellia oils: Possible mechanism of aldehydes formation in oleic acid triglyceride at high temperature" |
| Authors | Jun Cao, Xin Jiang, Qianyuan Chen, Hao Zhang, Huihui Sun, Wei-Min Zhang, Chuan Li (Hainan Univ.) |
| Venue | LWT - Food Science and Technology 2020, 108858 (PII S0023643819312009) |
| DOI | 10.1016/j.lwt.2019.108858 |
| Naming | "8-ROOH, 9-ROOH, 10-ROOH, 11-ROOH" = oleate hydroperoxides with OOH at that carbon. "A-scission" / "B-scission" are the paper's own labels (defined in §3, Fig 6) and are NOT the repo's R18a/R18b side A/B. "2-decenal", "2-undecenal", "2,4-decadienal" are the trans / trans,trans standards (§2). |
| Companions | Cao et al. 2014a JAFC 62:12545 (camellia oil 62 C / 35 d: octanal, nonanal, decanal as indicators); Cao et al. 2014b Food Res Int 64:901 (HPLC-QqQ-MS parent-fatty-acid assignment) |

## 1. Why it matters

The repo has no dossier that says which oleate hydroperoxide isomer gives nonanal, octanal, decanal,
2-decenal and 2-undecenal, so oleate -> nonanal requests are refused. This paper (i) measures the five
aldehydes as DNPH derivatives with authentic-standard calibration in two ~78-79 % oleic oils at three
frying temperatures, (ii) draws the four hydroperoxide isomers and the two scissions of each (Fig 6),
and (iii) states the assignment in words: decanal and 2-undecenal from 8-OOH; nonanal from 10-OOH and
9-OOH; octanal from 11-OOH; 2-decenal from 9-OOH. It also gives the temperature trend: the saturated
aldehydes (nonanal, octanal, decanal — their "B-scission") dominate at 120 C and fall at 150/180 C,
while the 2-alkenals (their "A-scission") rise or plateau. The substrate is a real oil, not an isolated
hydroperoxide, so isomer assignment is inference from product pattern, not from isomer feeding.

## 2. Methods as they matter to a model

- **Oils:** commercial olive oil and camellia oil (Haikou market). Fatty acids (area %, Table 1):
  olive 77.70 % 18:1, 5.82 % 18:2, 1.08 % 18:3; camellia 79.08 % 18:1, 9.79 % 18:2, 18:3 nd.
  Tocopherols: olive alpha 342.67 µg/g + gamma 25.66 µg/g; camellia alpha 68.53 µg/g. Initial DPPH
  scavenging 23.4 % (olive) vs 3.19 % (camellia). Initial PV 3.44 / 1.81 mmol/kg; AV 0.598 / 0.137 mg
  KOH/g; p-AV 9.70 / 4.27. Triacylglycerol substrate, i.e. esterified oleate.
- **Heating:** 500 mL oil in an open stainless-steel pan (19 cm bottom, 22 cm top, 8 cm high) with an
  automatic stirrer "in order to contact enough oxygen", in a silicone-oil fryer at 120 +/- 5, 150 +/- 5
  and 180 +/- 5 C. 50 mL sampled at 1, 2, 4, 8, 12, 24 h; sealed under N2, stored -18 C. Air, stirred,
  open surface; no added water, no food.
- **Carbonyls:** ~1 g oil + 2 mL DNPH reagent (0.3 g/L DNPH in MeOH with 3 % 6 N HCl), 40 C / 1 h
  shaking; extracted MeOH/H2O 75:25 twice, then water + CH2Cl2; CH2Cl2 layer dried under N2, taken up in
  MeCN, 0.22 µm filtered. HPLC-PDA 360 nm, Zorbax Eclipse Plus C18 4.6 x 250 mm 5 µm, MeOH/water
  gradient 75 -> 100 % MeOH over 30 min, 0.8 mL/min, 10 µL. **Calibration:** DNPH-hydrazone standards
  synthesised in-house from the aldehyde standards (TCI): trans-2-pentenal, trans-2-hexenal,
  trans-2-heptenal, heptanal, trans-2-octenal, octanal, trans-2-nonenal, nonanal,
  trans,trans-2,4-decadienal, trans-2-decenal, decanal (plus C1-C6 alkanal hydrazones bought). Peak
  areas "converted to molarity using calibration curves of carbonyl-DNPH standards". Compounds marked
  (a) in Table 2 (hexanone, heptanal, nonanone, decanone, 2-undecenal) were "calculated from
  hexanal-equivalents using the molecular weight of DNPH-aldehydes/ketones", i.e. NOT calibrated with
  their own standard. Unit **nmol/g oil**. n = 3, mean +/- SD, Duncan letters.
- **Fatty acids:** FAME, GC-FID, area %. **Tocopherols:** HPLC-fluorescence. **PV, AV, p-AV:** AOCS.
- **Heating is of whole oil**: no hydroperoxide isomer analysis was done; isomers are inferred.

## 3. Tables re-typed

### Table 1. "Fatty acids composition (area %) of oil in the initial and after heated for 24 h." Mean +/- SD, n = 3 (letters omitted; nd = not detected)

| fatty acid | olive initial | olive 120 C | olive 150 C | olive 180 C | camellia initial | camellia 120 C | camellia 150 C | camellia 180 C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| C12:0 | 0.16 +/- 0.07 | 0.06 +/- 0.01 | 0.04 +/- 0.01 | 0.07 +/- 0.01 | 0.04 +/- 0.01 | 0.05 +/- 0.01 | 0.06 +/- 0.02 | 0.06 +/- 0.01 |
| C14:0 | 0.10 +/- 0.04 | 0.06 +/- 0.01 | 0.03 +/- 0.01 | nd | 0.02 +/- 0.01 | 0.05 +/- 0.01 | 0.05 +/- 0.01 | 0.04 +/- 0.01 |
| C16:0 | 9.97 +/- 0.03 | 9.89 +/- 0.04 | 10.41 +/- 0.12 | 10.90 +/- 0.01 | 7.96 +/- 0.01 | 8.23 +/- 0.10 | 8.67 +/- 0.20 | 9.30 +/- 0.20 |
| tC16:1 | 0.16 +/- 0.02 | 0.14 +/- 0.01 | 0.11 +/- 0.03 | 0.07 +/- 0.01 | 0.14 +/- 0.03 | 0.06 +/- 0.01 | 0.14 +/- 0.10 | 0.07 +/- 0.01 |
| cC16:1 | 0.78 +/- 0.06 | 0.69 +/- 0.03 | 0.66 +/- 0.06 | 0.60 +/- 0.01 | nd | nd | nd | nd |
| C17:0 | 0.11 +/- 0.01 | 0.11 +/- 0.01 | 0.12 +/- 0.01 | 0.12 +/- 0.01 | 0.07 +/- 0.02 | 0.06 +/- 0.01 | 0.07 +/- 0.02 | 0.07 +/- 0.01 |
| cC17:1 | 0.23 +/- 0.01 | 0.27 +/- 0.01 | 0.18 +/- 0.01 | 0.14 +/- 0.01 | 0.08 +/- 0.02 | 0.10 +/- 0.03 | 0.10 +/- 0.01 | 0.11 +/- 0.02 |
| C18:0 | 3.65 +/- 0.02 | 3.63 +/- 0.10 | 3.85 +/- 0.03 | 4.04 +/- 0.10 | 2.39 +/- 0.04 | 2.38 +/- 0.10 | 2.49 +/- 0.01 | 2.67 +/- 0.09 |
| **9cC18:1** | **77.70 +/- 0.38** | 78.24 +/- 0.29 | 78.93 +/- 0.42 | 78.98 +/- 0.04 | **79.08 +/- 0.07** | 80.37 +/- 0.10 | 80.49 +/- 0.37 | 80.04 +/- 0.30 |
| **9c12cC18:2** | **5.82 +/- 0.04** | 5.66 +/- 0.01 | 4.60 +/- 0.10 | 4.18 +/- 0.01 | **9.79 +/- 0.01** | 8.06 +/- 0.01 | 7.20 +/- 0.11 | 7.00 +/- 0.03 |
| C20:0 | nd | nd | nd | nd | 0.11 +/- 0.01 | 0.08 +/- 0.00 | 0.12 +/- 0.01 | 0.12 +/- 0.01 |
| 9c12c15cC18:3 | 1.08 +/- 0.01 | 0.99 +/- 0.07 | 0.82 +/- 0.01 | 0.66 +/- 0.03 | nd | nd | nd | nd |
| cC20:1 | 0.24 +/- 0.01 | 0.25 +/- 0.01 | 0.23 +/- 0.01 | 0.23 +/- 0.01 | 0.41 +/- 0.01 | 0.44 +/- 0.01 | 0.44 +/- 0.01 | 0.44 +/- 0.01 |
| sum SFA | 13.99 +/- 0.04 | 13.76 +/- 0.04 | 14.45 +/- 0.15 | 15.13 +/- 0.01 | 10.59 +/- 0.01 | 10.85 +/- 0.02 | 11.46 +/- 0.18 | 12.26 +/- 0.30 |
| sum UFA | 85.85 +/- 0.41 | 86.01 +/- 0.03 | 85.44 +/- 0.18 | 84.80 +/- 0.11 | 89.36 +/- 0.08 | 89.09 +/- 0.02 | 88.40 +/- 0.28 | 87.67 +/- 0.30 |
| sum cis MUFA | 78.95 +/- 0.35 | 79.45 +/- 0.11 | 80.02 +/- 0.29 | 79.96 +/- 0.11 | 79.57 +/- 0.07 | 80.81 +/- 0.11 | 80.93 +/- 0.35 | 80.48 +/- 0.30 |
| sum PUFA | 6.90 +/- 0.06 | 6.65 +/- 0.07 | 5.42 +/- 0.11 | 4.84 +/- 0.01 | 9.79 +/- 0.01 | 8.19 +/- 0.01 | 7.37 +/- 0.08 | 7.08 +/- 0.02 |
| sum n-3 PUFA | 1.08 +/- 0.01 | 0.99 +/- 0.07 | 0.82 +/- 0.01 | 0.66 +/- 0.03 | nd | nd | nd | nd |
| sum n-6 PUFA | 5.82 +/- 0.07 | 5.66 +/- 0.01 | 5.68 +/- 0.01 | 5.28 +/- 0.02 | 9.79 +/- 0.01 | 8.19 +/- 0.01 | 7.37 +/- 0.08 | 7.08 +/- 0.02 |

(Note: the printed "sum n-6 PUFA" row for olive at 150/180 C (5.68, 5.28) does not equal the 18:2 row
(4.60, 4.18); one of the two is a misprint. Use the 18:2 row.)

### Table 2. "Carbonyls composition (nmol/g) of oils in the initial and after heated for 24 h." Mean +/- SD, n = 3; nd = not detected. (a) = hexanal-equivalent calibration.

| carbonyl | olive initial | olive 120 C | olive 150 C | olive 180 C | camellia initial | camellia 120 C | camellia 150 C | camellia 180 C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| hexanal | 8.60 +/- 1.12 | 5.34 +/- 0.18 | 9.19 +/- 0.56 | 6.80 +/- 1.08 | 2.09 +/- 0.10 | 22.32 +/- 0.95 | 31.51 +/- 2.89 | 12.96 +/- 0.03 |
| hexanone (a) | nd | nd | 6.82 +/- 0.53 | 10.66 +/- 1.54 | nd | 12.75 +/- 0.49 | 14.12 +/- 8.53 | 12.42 +/- 0.21 |
| 2-heptenal | 8.55 +/- 0.20 | 7.20 +/- 0.17 | 8.57 +/- 1.38 | 9.64 +/- 1.31 | nd | 4.83 +/- 0.08 | 5.93 +/- 1.07 | 5.84 +/- 0.30 |
| heptanal (a) | nd | nd | nd | nd | nd | 7.62 +/- 0.21 | 8.57 +/- 0.10 | 9.63 +/- 0.18 |
| 2-octenal | nd | nd | 1.61 +/- 0.08 | 4.79 +/- 0.87 | nd | 6.40 +/- 0.08 | 5.13 +/- 0.23 | 4.00 +/- 0.21 |
| **octanal** | nd | 7.06 +/- 0.04 | 36.51 +/- 2.86 | 24.42 +/- 3.65 | nd | 33.57 +/- 3.03 | 22.54 +/- 3.57 | 15.06 +/- 0.93 |
| 2-nonenal | nd | 1.45 +/- 0.01 | 6.28 +/- 0.27 | 3.87 +/- 0.07 | nd | 6.36 +/- 0.20 | 8.18 +/- 0.75 | 9.10 +/- 1.00 |
| **nonanal** | 5.45 +/- 0.39 | 21.2 +/- 0.57 | 85.82 +/- 2.65 | 51.37 +/- 2.17 | nd | 67.37 +/- 5.22 | 43.71 +/- 7.76 | 34.21 +/- 0.11 |
| 2,4-decadienal | nd | 3.13 +/- 0.09 | 11.25 +/- 0.06 | 7.38 +/- 1.41 | nd | 16.54 +/- 0.59 | 21.09 +/- 1.82 | 16.77 +/- 1.28 |
| nonanone (a) | nd | nd | nd | nd | nd | 11.42 +/- 0.52 | 11.23 +/- 1.62 | 10.10 +/- 0.08 |
| **2-decenal** | nd | 1.29 +/- 0.05 | 22.87 +/- 0.12 | 24.22 +/- 3.52 | nd | 13.00 +/- 0.44 | 26.72 +/- 5.19 | 28.06 +/- 3.91 |
| **decanal** | nd | 7.36 +/- 0.22 | 25.91 +/- 0.37 | 13.11 +/- 1.85 | nd | 26.20 +/- 2.06 | 18.12 +/- 2.17 | 12.76 +/- 0.11 |
| decanone (a) | 2.07 +/- 0.20 | 4.28 +/- 0.16 | 7.15 +/- 0.45 | 7.45 +/- 1.30 | nd | 4.89 +/- 0.17 | 10.24 +/- 0.04 | 13.24 +/- 4.21 |
| **2-undecenal (a)** | nd | 1.79 +/- 0.10 | 28.47 +/- 0.18 | 25.35 +/- 2.84 | nd | 17.86 +/- 0.62 | 36.25 +/- 5.53 | 34.30 +/- 3.74 |

Numbers in the text not in the tables: nonanal 5.45 nmol/g in unheated olive oil; PV of camellia oil
150 C / 12 h 22.69 vs 120 C 23.09 mmol/kg; 150 C / 24 h 28.42 vs 120 C 60.89 mmol/kg; "both oils
demonstrated lower increments of the peroxide value at 180 C". Fatemi & Hammond 1980 relative
oxidation rates oleate:linoleate:linolenate 1:10.3:21.6 (cited). Time courses of AV, PV, p-AV, DPPH,
tocopherols and of the five aldehydes (Figs 2-5): FIGURE-ONLY.

### Figure 6 — "Formation pathways of hydroperoxide isomers (8-ROOH, 9-ROOH, 10-ROOH and 11-ROOH, I) and aldehydes (octanal, nonanal, decanal, 2-decenal and 2-undecenal, II-V) in the thermal oxidation of oleic acid triglyceride" (drawn; described)

Panel I (drawn on a triacylglycerol with the oleoyl chain numbered): (1) H lost at C8 -> + O2, H ->
**8-ROOH** (C9=C10 kept); (2) H lost at C8, "rearranged" (radical to C10, new C8=C9) -> **10-ROOH**;
(3) H lost at C11 -> **11-ROOH** (C9=C10 kept); (4) H lost at C11, rearranged (radical to C9, new
C10=C11) -> **9-ROOH**. Panels II-V: each ROOH loses •OH to the alkoxyl radical, then two C-C
scissions labelled A and B.

The paper's convention (read off the panels and the text): **A-scission** = the C-C bond on the side of
the alkoxyl carbon AWAY from the double bond (gives an alpha,beta-unsaturated aldehyde retaining the
C=C, plus a saturated alkyl radical); **B-scission** = the C-C bond BETWEEN the alkoxyl carbon and the
vinyl carbon (gives a saturated aldehyde plus a vinyl-type radical; the drawings add "+ •OH" on the B
branch to finish the second fragment).

| panel | isomer (C=C) | A-scission bond -> products | B-scission bond -> products |
|---|---|---|---|
| II | 8-ROOH (9=10) | C7-C8 -> **2-undecenal** (C8-C18) + heptanoyl-glyceride radical | C8-C9 -> **decanal** (C9-C18, after •OH) + 8-oxo-octanoyl core aldehyde |
| III | 9-ROOH (10=11) | C8-C9 -> **2-decenal** (C9-C18) + octanoyl-glyceride radical | C9-C10 -> **nonanal** (C10-C18) + 9-oxo-nonanoyl core aldehyde |
| IV | 10-ROOH (8=9) | C10-C11 -> 10-oxo-8-decenoyl core aldehyde + octyl radical (no volatile aldehyde drawn) | C9-C10 -> **nonanal** (C10-C18) + vinyl core radical |
| V | 11-ROOH (9=10) | C11-C12 -> 11-oxo-9-undecenoyl core aldehyde + heptyl radical | C10-C11 -> **octanal** (C11-C18) + vinyl core radical |

Authors' words (Section 3.5): "Homolytic A-scission of the carbon-carbon bond (C7-C8) produces
2-undecenal, or B-scission of the carbon-carbon bond (C8-C9) produces decanal." "Homolytic A-scission
of the carbon-carbon bond (C10-C11) produces a triglyceride core aldehyde ... and saturated alkyl
radical or B-scission of the carbon-carbon bond (C9-C10) produces nonanal." "A-scission of the
carbon-carbon bond (C11-C12) produces triglyceride core aldehyde and saturated alkyl radical, or
B-scission of the carbon-carbon bond (C10-C11) produces octanal." "Homolytic A-scission of the
carbon-carbon bond (C8-C9) produces 2-decenal, or B-scission of the carbon-carbon bond (C9-C10)
produces nonanal."

Temperature statement: "At 120 C, the illustrated aldehydes in olive and camellia oils mainly increased
in the order nonanal > octanal ~ decanal > 2-undecenal ~ 2-decenal. ... homolytic B-scission of the
carbon-carbon bond in the alkoxy radical was predominant at 120 C. The yields of nonanal, octanal and
decanal decreased at 150 C and 180 C, whereas the unsaturated aldehyde, 2-decenal and 2-undecenal
yields increased or tended to stabilize. This reflects that the higher the temperature, the easier it
is for the formation and homolytic A-scission of 8-ROOH and 9-ROOH to occur."

## 4. Routes and numbers the repository can use

All levels: nmol/g oil at 24 h, n = 3, DNPH-HPLC, whole oil (~78-79 % oleate, 6-10 % linoleate).

| route | reactant -> product | mechanism as drawn (figure) | measured numbers with units and conditions | evidence class |
|---|---|---|---|---|
| OL-10-B | oleate 10-OOH -> **nonanal** + 9-oxononanoyl core | alkoxyl C10, B-scission C9-C10 (Fig 6 IV) | nonanal 24 h: olive 21.2 (120 C), 85.82 (150 C), 51.37 (180 C); camellia 67.37 / 43.71 / 34.21 nmol/g | mechanism_drawn; measured_level |
| OL-9-B | oleate 9-OOH -> **nonanal** + 9-oxononanoyl core | alkoxyl C9, B-scission C9-C10 (Fig 6 III) | shares the nonanal numbers above; the paper does not split nonanal between 9- and 10-OOH | mechanism_drawn; measured_level (shared) |
| OL-9-A | oleate 9-OOH -> **2-decenal** + octanoyl-glyceride radical | A-scission C8-C9 (Fig 6 III) | 2-decenal: olive 1.29 / 22.87 / 24.22; camellia 13.00 / 26.72 / 28.06 nmol/g | mechanism_drawn; measured_level |
| OL-11-B | oleate 11-OOH -> **octanal** + vinyl core radical (-> 10-oxodecanoyl) | B-scission C10-C11 (Fig 6 V) | octanal: olive 7.06 / 36.51 / 24.42; camellia 33.57 / 22.54 / 15.06 nmol/g | mechanism_drawn; measured_level |
| OL-8-B | oleate 8-OOH -> **decanal** + 8-oxooctanoyl core | B-scission C8-C9 (Fig 6 II) | decanal: olive 7.36 / 25.91 / 13.11; camellia 26.20 / 18.12 / 12.76 nmol/g | mechanism_drawn; measured_level |
| OL-8-A | oleate 8-OOH -> **2-undecenal** + heptanoyl-glyceride radical | A-scission C7-C8 (Fig 6 II) | 2-undecenal (hexanal-equivalent calibration): olive 1.79 / 28.47 / 25.35; camellia 17.86 / 36.25 / 34.30 nmol/g | mechanism_drawn; measured_level (semi-quantitative) |
| OL-10-A / OL-11-A | 10-OOH / 11-OOH -> core aldehyde + octyl / heptyl radical | A-scission (Fig 6 IV, V); no volatile aldehyde drawn | none (heptanal (a) nd in olive, 7.6-9.6 in camellia) | mechanism_drawn |
| OL-T | temperature shift saturated -> unsaturated aldehydes | text §3.5 | within-oil ratios at 24 h, (2-decenal + 2-undecenal)/(octanal + nonanal + decanal): olive 0.087 (120 C), 0.35 (150 C), 0.56 (180 C); camellia 0.24 / 0.75 / 1.01 | measured_ratio |

Other numbers: linoleate-derived hexanal, 2-heptenal, 2-octenal, 2,4-decadienal, 2-nonenal are in
Table 2 too (2,4-decadienal up to 21.09 nmol/g in camellia 150 C), but the oils contain 6-10 %
linoleate, so these rows are not clean oleate products.

## 5. Rule sketches (repository suggestions, not the paper's)

Keys: `nonanal` exists in `data/keys/compounds.yml`; **octanal, decanal, (E)-2-decenal, 2-undecenal do
not.** `data/species/structures.yml` has `NONANAL` and the lump `LOOH_OL` ("oleate hydroperoxides, four
positional isomers lumped") but no per-isomer oleate hydroperoxide SMILES and no oleic acid / methyl
oleate / triolein entry. Because the paper assigns products per isomer, the lump must be split into
four molecules for these rules to have positive controls. Controls below use methyl esters to match the
repo's LOOH_9/13 convention; the free-acid analogues differ only in the ester end.

Proposed isomer SMILES (methyl esters, trans C=C as the autoxidation product):
- OL_8_OOH (C9=C10): `CCCCCCCC/C=C/C(OO)CCCCCCC(=O)OC`
- OL_9_OOH (C10=C11): `CCCCCCC/C=C/C(OO)CCCCCCCC(=O)OC`
- OL_10_OOH (C8=C9): `CCCCCCCCC(OO)/C=C/CCCCCCC(=O)OC`
- OL_11_OOH (C9=C10): `CCCCCCCC(OO)/C=C/CCCCCCCC(=O)OC`

**S1 (the paper's "B-scission", saturated aldehyde on the methyl side).** Pattern: `[#6:1][CH1:2]([O:3][OH])[CH1:4]=[CH1:5][#6:6]` where atom 1 is the METHYL-side chain. Change: bond C2-C4 breaks; C2 becomes CH=O; C4 (the vinyl carbon) ends as an aldehyde carbon of the core fragment (net, after •OH). This is the same pattern shape as R18b but with a mono-ene instead of a diene; the ordering of chain vs ester side must be enforced by mapping, since 8-/9-OOH have the C=C on the methyl side and 10-/11-OOH on the ester side.
- OL_10_OOH -> `CCCCCCCCC=O` (nonanal, key `nonanal`) + `O=CCCCCCCCC(=O)OC` (ME_9_OXONONANOATE)
- OL_11_OOH -> `CCCCCCCC=O` (octanal) + `O=CCCCCCCCCC(=O)OC` (methyl 10-oxodecanoate)
- OL_8_OOH -> `CCCCCCCCCC=O` (decanal) + `O=CCCCCCCC(=O)OC` (methyl 8-oxooctanoate)
- OL_9_OOH -> `CCCCCCCCC=O` (nonanal) + ME_9_OXONONANOATE
- negative controls: methyl stearate `CCCCCCCCCCCCCCCCCC(=O)OC` (no OOH, no C=C); hexanal `CCCCCC=O`; LOOH_13_ct should be handled by R18b, not by this rule (diene) — decide whether the mono-ene pattern should exclude the conjugated diene explicitly.

**S2 (the paper's "A-scission", 2-alkenal on the methyl side; only 8- and 9-OOH give a volatile).** Pattern: `[#6:1][CH2:7][CH1:2]([O:3][OH])[CH1:4]=[CH1:5][#6:6]` with atom 7 on the ESTER side. Change: bond C7-C2 breaks; C2 becomes CH=O conjugated to C4=C5 (2-alkenal, methyl side); C7 becomes a CH3 (alkane end of the core fragment, net after H transfer) — the shape of R18a.
- OL_9_OOH -> `CCCCCCC/C=C/C=O` ((E)-2-decenal) + `CCCCCCCC(=O)OC` (ME_OCTANOATE)
- OL_8_OOH -> `CCCCCCCC/C=C/C=O` ((E)-2-undecenal) + `CCCCCCC(=O)OC` (methyl heptanoate)
- OL_10_OOH / OL_11_OOH: this rule fires but the volatile is an alkane (octane, heptane), the aldehyde stays on the glyceride; the paper draws no volatile aldehyde for these and measured no heptanal in olive oil.
- negative: methyl stearate; nonanal.

**S3 (branch ratio, not a rule).** If a later wave wants a temperature-dependent B/A split for oleate,
the within-oil ratios in OL-T (0.087 -> 0.56 olive; 0.24 -> 1.01 camellia, 120 -> 180 C, 24 h) are the
only measured handle in this paper; they are end-point ratios of accumulated aldehydes after 24 h of
further reaction, not primary scission ratios.

## 6. Flags

1. **Isomer assignment is inference from an oil, not from fed isomers.** No hydroperoxide isomer was
   measured; the mapping 8/9/10/11-OOH -> aldehydes is the classic Frankel scheme restated. The
   measured content is the five aldehyde levels vs temperature and time.
2. **Cao's A/B labels are not the repo's R18a/R18b A/B**, and B-scission as drawn produces a vinyl
   radical (the step Miyazaki 2023 calls unfavourable for HpODE). The paper does not discuss radical
   stability; it labels bonds geometrically.
3. **Two calibrations:** octanal, nonanal, decanal, 2-decenal are calibrated with their own DNPH
   standards; 2-undecenal (and hexanone, heptanal, nonanone, decanone) are hexanal-equivalents by
   molecular-weight scaling — the 2-undecenal column is semi-quantitative.
4. **24 h open-pan heating with stirring at 120-180 C** is far past the point where aldehydes are
   themselves consumed (the paper says so: "at higher temperatures, the formation and decomposition of
   aldehydes were very complicated"); levels at 24 h are net accumulations, not yields. Time courses
   exist only as Figure 4 (FIGURE-ONLY); the 1-12 h numbers are not printed.
5. **Both oils contain linoleate (5.8 / 9.8 %) and olive contains linolenate (1.1 %)**; hexanal,
   2-heptenal, 2-octenal, 2,4-decadienal, 2-nonenal in Table 2 come from those. Camellia's higher
   linoleate also explains its higher 2,4-decadienal.
6. **"Article in press" proof:** page numbers "xxx"; the sum n-6 PUFA row in Table 1 is inconsistent
   with the 18:2 row for olive at 150/180 C (see note under Table 1).
7. **Triacylglycerol substrate**: the B-scission co-products are glyceride core aldehydes (Sjovall
   2002/2003 cited), not free oxo-acids; the volatile products are the same as for free oleic acid.
8. No DFT, no rate constants; nothing to mark inadmissible.
