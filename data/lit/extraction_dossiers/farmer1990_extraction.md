# Farmer & Mottram 1990 — EXTRACTION (cysteine + ribose with beef triglyceride, beef phospholipid, egg PC or egg PE, 0.5 M phosphate pH 5.6, 140 C / 1 h; 60 volatiles as relative peak areas against 1 ng 1,2-dichlorobenzene = 100, n = 4-5, with SDs)
### The four-lipid follow-up to Whitfield 1988: the drawn 2,4-decadienal routes to 2-pentylpyridine, 2-hexylthiophene and 2-pentylthiapyran, and the lipid-class dependence of the thiol and mercaptoketone suppression.

**Source on disk:** `data/articles/farmer1990.pdf` (owner's download, 2026-09-08). Scanned journal
pages with an OCR text layer. Text reads cleanly; Tables 1 and 2 (pp. 511-514) are partly scrambled
in the plain layer, so they were re-extracted with `pypdf` in layout mode and every cell was checked
against renderings of the four pages. Figure 1 (p. 517, the drawn mechanism) was rendered and is
described atom by atom in §4; no numbers were read from it. Table 3 (odour descriptors) from the
text layer.

## 0. Identity

| field | value |
|---|---|
| Title | "Interaction of Lipid in the Maillard Reaction between Cysteine and Ribose: the Effect of a Triglyceride and Three Phospholipids on the Volatile Products" |
| Authors | Linda J. Farmer (formerly Salter), Donald S. Mottram (AFRC Institute of Food Research, Bristol) |
| Venue | J. Sci. Food Agric. 53 (1990) 505-525; received 11 Dec 1989, accepted 27 May 1990 |
| DOI | not printed on the scan |
| Companions | Farmer, Mottram, Whitfield 1989 (JSFA 49, 347-368): identifications for cysteine + ribose + PC; Farmer & Mottram 1992 (JSFA 60, 489-497, cited here as "to be published"): the four lipids heated alone and with cysteine + ribose, the source of the 2,4-decadienal proportions quoted in §4; Whitfield 1988 (`whitfield1988_extraction.md`): the ng-scale predecessor |
| Lipid labels | BTG = beef subcutaneous triglyceride (Florisil-cleaned); BPL = beef M. semimembranosus phospholipid (silicic acid column; 1-3 % triglyceride carried, 3-10 % of the phospholipid lost); PC = Sigma egg phosphatidylcholine type III-E (c. 99 %); PE = Sigma egg phosphatidylethanolamine type III (c. 98 %) |

## 1. Why it matters

Same cysteine + ribose pot as Whitfield 1988 (41 mM cysteine, 30 mM ribose, 15 g/L lipid, 140 C,
1 h) but with four lipids and replicate SDs, and with the mechanism figure the earlier paper only
cited. For the sulfur lane it gives the lipid-class dependence of the thiol quench: MFT falls to 0.40
(BTG), 0.15 (BPL), 0.27 (PC), 0.24 (PE) of its lipid-free area; FFT to 0.62-0.72 with every lipid;
2-methyl-3-thiophenethiol to 0.076 (BTG) and 0.0025 (BPL); the four mercaptoketones to 0.6-1.0 with
triglyceride and about 0.5 with every phospholipid. For the lipid–Maillard products it gives the
lipid ranking PC >> BPL > PE >> BTG for 2-pentylpyridine (0 / 26 / 279 / 5210 / 429), 2-pentylthiapyran
(0 / 35 / 3150 / 34700 / 12500), 2-pentylthiophene (0 / 0 / 109 / 2700 / 1230) and the alkanethiols,
and it draws the three decadienal routes. Two structural facts for the rule writer: the triglyceride
(0.6 % 18:2, 0.4 % PUFA) makes almost none of the lipid–Maillard products, and PE's amino group did
not make it a better NH3 donor than PC for 2-pentylpyridine.

## 2. Methods as they matter to a model

- **Charge (per ampoule):** 30 mg lipid (from solution, solvent removed under N2 at 40-50 C) + 2 mL
  of cysteine 5 mg/mL + ribose 4.5 mg/mL in 0.5 M sodium phosphate pH 5.6. Molar (paper and mine):
  cysteine 41 mmol/L ("41 mM of cysteine amino groups"), ribose 30 mmol/L, lipid 15 g/L; PE amino
  groups c. 20 mmol/L, BPL amino groups c. 8 mmol/L (paper's values). Lipid-free control: the same
  2 mL without lipid. n = 4 per lipid, n = 5 for the lipid-free pot.
- **Vessel:** 2.5 mL Pyrex freeze-drying ampoules, flamed; N2 blown over 2 min; flame-sealed (about
  0.5 mL headspace); ultrasonic bath 40 C, 1 h to disperse the lipid.
- **Heating:** 140 C, 1 h, autoclave. Stored at −20 C for 1-2 days if not collected at once.
- **pH:** 5.6 initial; **fell by less than 0.5** during reaction (versus 1.1 at 0.2 M in Whitfield
  1988).
- **Volatile collection:** 2 mL diluted to 20 mL in phosphate buffer, 60 C, N2 50 mL/min onto Tenax
  GC (collection time not restated; 30 min in Whitfield 1988). Internal standard mixture on the
  trap before collection: C10-C24 alkanes 10 ng each + **1,2-dichlorobenzene 50 ng**.
- **GC-MS:** CPWAX 52CB 50 m x 0.32 mm, Unijector thermal desorption at 250 C, Finnigan 4000 / Incos
  2100; MS conditions as Whitfield 1988 (EI 40 eV).
- **Quantification: relative peak areas, not ng.** For each compound one characteristic ion (the
  molecular ion except no. 29 m/z 74, no. 45 m/z 57, no. 46 m/z 75, no. 52 m/z 93, no. 54 m/z 107); the
  ion area divided by that ion's relative abundance in the compound's spectrum (an approximation to
  the total ion area), then scaled so that 1 ng of 1,2-dichlorobenzene on the trap = 100. No response
  factors, no authentic-standard calibration ("the unavailability of authentic samples for many of
  the compounds quantified made it impossible to measure the absolute quantities"). Values above
  1000 are given to three significant figures. Areas are comparable across the five pots for one
  compound (the within-study ratio); across compounds only roughly.
- **Lipid analysis:** TLC purity (BTG: < 10 mg/g phospholipid, < 10 mg/g free acid / partial
  glycerides; BPL: 10-30 mg/g triglyceride); micro-Kjeldahl N (BPL 18.2 ± 0.89 mg/g vs 18-19 calculated;
  BTG < 0.5 mg/g); fatty acids by GC after saponification and diazomethane methylation (Table 1).
- **Odour:** sniffed at opening by three flavour chemists (Table 3); no GC-O, no thresholds.
- **Categories** assigned by the authors to each compound: A formed only with lipid; B reduced by
  lipid, with a = phospholipids suppress more than triglyceride, b = PE and BPL suppress more than BTG
  and PC, c = triglyceride suppresses more than phospholipids, d = all lipids similar; C largely
  unaffected; D none of the above.

## 3. Tables re-typed

### Table 1. "Fatty acid compositions for beef triglyceride (BTG), beef phospholipid (BPL), egg phosphatidylcholine (PC) and egg phosphatidylethanolamine (PE)." Percent of total fatty acids. "—" = below 0.1 %. Identification: I = co-injection with authentic ester, L = literature spectrum, t = tentative.

| fatty acid | id | BTG | BPL | PC | PE |
|---|---|---:|---:|---:|---:|
| 14:0 | I | 3.2 | — | 0.1 | — |
| 14:1 | | 0.7 | — | — | — |
| 15:0 | I | 0.5 | — | — | — |
| 16:0 | I | 27.9 | 15.3 | 32.8 | 17.7 |
| 16:1 (Δ9) | I | 3.1 | 1.5 | 1.0 | 0.5 |
| 17:0 | I | 1.5 | 0.7 | 0.2 | 0.2 |
| 17:1 | | 0.9 | 0.7 | — | — |
| 18:0 | I | 18.3 | 13.4 | 12.9 | 26.2 |
| 18:1 (Δ9) | I | 35.6 | 29.0 | 29.6 | 20.1 |
| 18:1 (Δ11) | I | 2.7 | 1.8 | 1.0 | 0.9 |
| 18:2 (ω6) | I | 0.6 | 8.0 | 16.2 | 13.8 |
| 18:3 (ω3) | I | 0.4 | 3.2 | — | — |
| 20:3 (ω9) | L | — | 1.1 | — | — |
| 20:3 (ω6) | L | — | 2.1 | 0.4 | 0.4 |
| 20:4 (ω6) | I | — | 8.2 | 3.8 | 13.7 |
| 20:4 (other isomer) | | — | 1.0 | — | — |
| 20:5 (ω3) | t | — | 5.4 | — | — |
| 22:4 (ω6) | t | — | 0.4 | 0.2 | 0.8 |
| 22:5 (ω3) | I | — | 5.3 | 1.1 | 3.5 |
| 22:6 (ω3) | I | — | 1.6 | 0.7 | 2.3 |
| minor fatty acids | | 4.4 | 1.5 | 0.9 | 1.1 |
| Σ saturated | | 51.4 | 29.4 | 46.0 | 44.1 |
| Σ monounsaturated | | 43.0 | 33.0 | 31.6 | 21.5 |
| Σ diunsaturated | | 0.6 | 8.0 | 16.2 | 13.8 |
| Σ 3 or more double bonds | | 0.4 | 28.3 | 6.2 | 20.7 |

(Row sums reproduce the printed Σ lines to 0.1 for all four lipids, which fixes the 17:0 / 17:1
placement the layout extraction blurred.)

### Table 2. "Relative peak areas for selected heterocyclic components formed from the Maillard reaction between cysteine and ribose in the absence and presence of BTG, BPL, PC and PE." Mean (SD); relative to 1 ng 1,2-dichlorobenzene = 100. Cat. = authors' category. Identification: c = MS + IR interpretation, d = literature spectrum (footnote markers as printed).

| no. | compound | MW | alone | + BTG | + BPL | + PC | + PE | cat. |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| | **Alkylthiophenes and related** | | | | | | | |
| 1 | 2-methylthiophene | 98 | 651 (341) | 420 (414) | 1280 (1060) | 217 (69) | 546 (204) | D |
| 2 | 4,5-dihydro-2-methylthiophene | 100 | 1220 (561) | 889 (421) | 1320 (646) | 494 (78) | 714 (161) | D |
| 3 | 2-ethylthiophene | 112 | 367 (132) | 236 (150) | 949 (449) | 267 (134) | 323 (26) | D |
| 4 | 2,5-dimethylthiophene | 112 | 2330 (875) | 1770 (944) | 4130 (1030) | 2100 (731) | 2240 (271) | D |
| 5 | 2,3-dimethylthiophene | 112 | 597 (166) | 337 (173) | 966 (285) | 586 (129) | 530 (53) | D |
| 6 | **2-butylthiophene** | 140 | 0 | 0 | 33 (66) | 176 (150) | 176 (204) | A |
| 7 | **2-pentylthiophene** | 154 | 0 | 0 | 109 (112) | 2700 (1970) | 1230 (916) | A |
| 8 | **2-hexylthiophene** | 168 | 0 | 0 | 184 (290) | 1220 (452) | 436 (315) | A |
| 9 | **2-pentylthiapyran** c | 168 | 0 | 35 (33) | 3150 (1880) | 34700 (14800) | 12500 (8250) | A |
| 10 | 2-(1-hexenyl)thiophene c (cis/trans) | 166 | 0 | 0 | 0 | 65 (35) | 42 (47) | A |
| 11 | 2-(1-hexenyl)thiophene c (cis/trans) | 166 | 0 | 0 | 0 | 374 (217) | 287 (308) | A |
| | **Acylthiophenes** | | | | | | | |
| 12 | 2-formylthiophene | 112 | 3450 (1280) | 2460 (647) | 1390 (233) | 2340 (564) | 1940 (411) | Bd |
| 13 | 3-acetylthiophene | 126 | 1270 (282) | 1070 (339) | 901 (71) | 1330 (222) | 1170 (319) | Bd |
| 14 | 2-acetylthiophene | 126 | 387 (121) | 224 (60) | 313 (71) | 318 (53) | 330 (119) | Bc |
| 15 | 2-propionylthiophene | 140 | 3600 (592) | 1950 (579) | 2950 (234) | 3170 (252) | 3000 (768) | Bc |
| 16 | 2-formyl-3-methylthiophene | 126 | 6940 (1790) | 5830 (994) | 2500 (1450) | 5810 (1920) | 4970 (1150) | Bd |
| 17 | 3-ethyl-2-formylthiophene c | 140 | 3120 (826) | 2030 (400) | 1070 (583) | 1980 (911) | 1950 (352) | Bd |
| 18 | a dimethylformylthiophene c | 140 | 880 (348) | 403 (85) | 237 (150) | 425 (144) | 691 (234) | Bd |
| 19 | a thienylethanal c | 126 | 1150 (651) | 1990 (829) | 1570 (1010) | 3420 (908) | 1650 (791) | D |
| | **Heterocyclic thiols** | | | | | | | |
| 20 | **2-furanmethanethiol (FFT)** | 114 | 50200 (7090) | 33400 (3560) | 31600 (7530) | 31000 (2490) | 35900 (7130) | Bd |
| 21 | **2-methyl-3-furanthiol (MFT)** | 114 | 30600 (3650) | 12300 (3800) | 4740 (1400) | 8240 (1110) | 7230 (2590) | Bb |
| 22 | **2-thiophenethiol** | 116 | 38300 (4390) | 12300 (3730) | 1230 (861) | 17600 (3890) | 5240 (3100) | Bb |
| 23 | 3-thiophenethiol c | 116 | 2560 (1260) | 1410 (458) | 57 (66) | 1220 (306) | 536 (109) | Bb |
| 24 | **2-methyl-3-thiophenethiol** d | 130 | 32600 (4410) | 2490 (2380) | 83 (147) | 6290 (1470) | 920 (1390) | Bb |
| | **Thiophenones** | | | | | | | |
| 25 | dihydro-3(2H)-thiophenone | 102 | 2280 (699) | 1660 (502) | 1520 (553) | 1220 (292) | 1300 (533) | Bd |
| 26 | dihydro-2-methyl-3(2H)-thiophenone | 116 | 38700 (12600) | 27400 (7450) | 36500 (15000) | 24400 (2480) | 26500 (5930) | Bd |
| 27 | trans-dihydro-2,(4/5)-dimethyl-3(2H)-thiophenone c | 130 | 4610 (1540) | 2900 (287) | 3050 (477) | 3270 (511) | 3390 (769) | Bd |
| 28 | cis-dihydro-2,(4/5)-dimethyl-3(2H)-thiophenone c | 130 | 5370 (1310) | 4090 (593) | 3830 (336) | 4280 (427) | 4340 (741) | Bd |
| 29 | dihydro-(2/5)-ethyl-3(2H)-thiophenone c | 130 | 5210 (1790) | 2830 (1220) | 3070 (1340) | 3320 (507) | 3180 (1130) | Bd |
| | **Dithianones and trithianes** | | | | | | | |
| 30 | 1,2-dithian-4-one c | 134 | 904 (194) | 937 (195) | 1100 (500) | 1030 (213) | 693 (111) | C |
| 31 | 3-methyl-1,2-dithian-4-one d | 148 | 559 (288) | 406 (243) | 690 (786) | 674 (109) | 507 (220) | C |
| 32 | trans-3,(5/6)-dimethyl-1,2-dithian-4-one c | 162 | 315 (90) | 199 (137) | 311 (111) | 588 (90) | 348 (125) | C |
| 33 | 3-methyl-1,2,4-trithiane | 152 | 825 (504) | 194 (135) | 2080 (1300) | 712 (198) | 709 (396) | Bc |
| | **Bicyclic compounds** | | | | | | | |
| 34 | 2,3-dihydro-6-methylthieno[2,3-c]furan (kahweofuran) | 140 | 23600 (5200) | 10400 (5100) | 21000 (2300) | 25000 (1810) | 20600 (4420) | Bc |
| 35 | thieno[2,3-b]thiophene | 140 | 7920 (2390) | 1980 (761) | 5410 (874) | 6350 (934) | 4680 (1720) | Bc |
| 36 | a methylthienothiophene c | 154 | 2210 (720) | 271 (126) | 981 (335) | 820 (115) | 741 (340) | Bc |
| 37 | a dihydrothienothiophene c | 142 | 57200 (29200) | 25900 (8590) | 56700 (5310) | 71000 (6300) | 52600 (12600) | Bc |
| 38 | a methyldihydrothienothiophene c | 156 | 15900 (2140) | 2440 (1030) | 7340 (4140) | 5710 (1400) | 5190 (2350) | Bc |
| | **Thiazoles** | | | | | | | |
| 39 | thiazole | 85 | 307 (198) | 218 (141) | 222 (125) | 172 (25) | 218 (140) | C |
| 40 | 2-methylthiazole d | 99 | 247 (201) | 242 (211) | 168 (110) | 136 (23) | 146 (95) | C |
| 41 | 4,5-dimethylthiazole | 113 | 502 (151) | 406 (24) | 285 (44) | 438 (66) | 389 (112) | C |
| 42 | trimethylthiazole | 127 | 1430 (419) | 866 (253) | 869 (328) | 1450 (319) | 1060 (469) | C |
| 43 | 2-acetylthiazole | 127 | 980 (354) | 772 (376) | 879 (479) | 1500 (528) | 1010 (534) | C |
| | **Mercaptocarbonyls** | | | | | | | |
| 44 | 2-mercapto-3-butanone c (= 3-mercapto-2-butanone) | 104 | 19600 (6850) | 19100 (8120) | 8800 (2300) | 10000 (1770) | 9330 (3780) | Ba |
| 45 | 2-mercapto-3-pentanone c (`MP3P`) | 118 | 27100 (2710) | 19600 (1850) | 13300 (4300) | 14400 (2800) | 14300 (2580) | Ba |
| 46 | 3-mercapto-2-pentanone d (`MP2P`) | 118 | 28200 (1490) | 21700 (1340) | 13200 (3350) | 14100 (2590) | 13800 (2890) | Ba |
| 47 | 1-mercapto-3-pentanone c | 118 | 1390 (646) | 841 (243) | 550 (473) | 690 (92) | 591 (223) | Ba |
| | **Alkanethiols** | | | | | | | |
| 48 | 1-heptanethiol | 132 | 0 | 0 | 0 | 456 (351) | 113 (126) | A |
| 49 | 1-octanethiol | 146 | 0 | 0 | 0 | 341 (247) | 104 (111) | A |
| | **Compounds without sulphur** | | | | | | | |
| 50 | 2-furfural | 96 | 1425 (481) | 1140 (265) | 1070 (394) | 705 (121) | 737 (159) | D |
| 51 | 1-(2-furyl)-2-propanone d | 124 | 2590 (694) | 2510 (489) | 2090 (122) | 2710 (150) | 2910 (593) | C |
| 52 | **2-pentylpyridine** | 149 | 0 | 26 (8) | 279 (90) | 5210 (2240) | 429 (239) | A |
| | Pyrazines | | | | | | | |
| 53 | methylpyrazine | 94 | 478 (243) | 555 (192) | 548 (221) | 300 (61) | 570 (74) | C |
| 54 | ethylpyrazine | 108 | 1009 (261) | 1140 (369) | 1040 (377) | 1280 (149) | 1140 (237) | C |
| 55 | 2,3-dimethylpyrazine | 108 | 362 (142) | 393 (155) | 369 (137) | 427 (58) | 392 (111) | C |
| 56 | 2-ethyl-5-methylpyrazine | 122 | 217 (43) | 187 (65) | 179 (74) | 226 (20) | 205 (86) | C |
| 57 | trimethylpyrazine | 122 | 94 (52) | 116 (55) | 64 (37) | 128 (15) | 117 (59) | C |
| | Oxazoles | | | | | | | |
| 58 | 5-ethyl-4-methyloxazole d | 111 | 325 (131) | 273 (39) | 280 (48) | 120 (17) | 200 (57) | D |
| 59 | 4-ethyl-5-methyloxazole | 111 | 2180 (536) | 1900 (96) | 2060 (295) | 997 (237) | 1470 (411) | D |
| 60 | trimethyloxazole | 111 | 231 (121) | 268 (104) | 229 (92) | 198 (47) | 235 (56) | C |

Footnote: compounds 23, 44 and 45 were not in the 1989 identification paper; 23 by comparison with
Heller & Milne 1978; 44 and 45 by interpretation of the spectra printed in the footnote (no. 44 base
peak m/z 43, M+ 104 at 20 %; no. 45 base peak 57, M+ 118 at 17 %).

### Table 3. Odour of the opened ampoules (three assessors, summary as printed)

| pot | descriptors |
|---|---|
| cysteine + ribose alone | strong sulphurous, rubber, H2S, slight meaty (ham) under |
| + BTG | strong sulphurous, H2S, some meaty (ham, roast, boiled) notes |
| + BPL | distinctly meaty (chicken, roasted) under sulphurous and H2S |
| + PC | predominantly sulphurous, H2S, rubber, meaty (ham, boiled) undertones |
| + PE | very distinct meaty (strong, roast, lamb, boiled) with some sulphurous, rubber; "the most meaty"; a roast-meat aroma perceptible at 2-3 m |

### Derived within-study ratios (mine, from Table 2; area with lipid / area without)

| compound | BTG | BPL | PC | PE | Whitfield 1988 (+PC / −), ng |
|---|---:|---:|---:|---:|---:|
| MFT (21) | 0.40 | 0.15 | 0.27 | 0.24 | 0.34 |
| FFT (20) | 0.67 | 0.63 | 0.62 | 0.72 | 0.50 |
| 2-thiophenethiol (22) | 0.32 | 0.032 | 0.46 | 0.14 | 0.29 |
| 3-thiophenethiol (23) | 0.55 | 0.022 | 0.48 | 0.21 | 0.47 |
| 2-methyl-3-thiophenethiol (24) | 0.076 | 0.0025 | 0.19 | 0.028 | not measured |
| 2-mercapto-3-butanone (44) | 0.97 | 0.45 | 0.51 | 0.48 | |
| 2-mercapto-3-pentanone (45) | 0.72 | 0.49 | 0.53 | 0.53 | |
| 3-mercapto-2-pentanone (46) | 0.77 | 0.47 | 0.50 | 0.49 | |
| 1-mercapto-3-pentanone (47) | 0.61 | 0.40 | 0.50 | 0.43 | |
| dihydro-2-methyl-3(2H)-thiophenone (26) | 0.71 | 0.94 | 0.63 | 0.68 | 0.47 (2-methyltetrahydrothiophen-3-one) |
| 2-furfural (50) | 0.80 | 0.75 | 0.49 | 0.52 | 0.59 |
| methylpyrazine (53) | 1.16 | 1.15 | 0.63 | 1.19 | 0.47 |
| 2-methylthiazole (40) | 0.98 | 0.68 | 0.55 | 0.59 | 0.33 |
| 4,5-dimethylthiazole (41) | 0.81 | 0.57 | 0.87 | 0.78 | 0.58 |
| 2-acetylthiazole (43) | 0.79 | 0.90 | 1.53 | 1.03 | 1.70 |

Lipid–Maillard products, PC : PE ratio: 2-pentylpyridine 12.1; 2-pentylthiapyran 2.8;
2-pentylthiophene 2.2; 2-hexylthiophene 2.8; 1-heptanethiol 4.0. Within the PC pot,
2-pentylthiapyran : 2-hexylthiophene = 28 (PE 29; BPL 17; BTG 35 : 0).

## 4. Routes and numbers the repository can use

| route or quantity | reactant -> product | mechanism as drawn (Figure 1) | measured numbers, units, conditions | evidence class |
|---|---|---|---|---|
| F90-A 2-pentylpyridine | 2,4-decadienal + NH3 -> 2-pentylpyridine | left branch: NH3 adds to the aldehyde carbon C1 (counting CHO as C1) giving the carbinolamine C5H11-CH=CH-CH=CH-CH(OH)-NH2; the N lone pair attacks C5 (the pentyl-bearing terminus of the diene) with the double bonds shifting and the C1 hydroxyl leaving, closing a six-membered N ring (2-pentyl-1,2-dihydropyridine as drawn, NH and two C=C); "[O]" gives 2-pentylpyridine, pentyl on C2 = the former C5 | relative areas 0 / 26 (8) / 279 (90) / 5210 (2240) / 429 (239) for alone / BTG / BPL / PC / PE; the lipids heated alone gave 2,4-decadienal in proportions BTG 1.4 : BPL 1.3 : PC 64 : PE 1 (Farmer & Mottram 1992, quoted) | mechanism_drawn + within_study_ratio |
| F90-B 2-hexylthiophene | 2,4-decadienal + H2S -> 2-hexylthiophene | middle branch: H2S adds S to C4 and H to C5 (1,4-addition) giving C5H11-CH2-CH(SH)-CH=CH-CHO; the S lone pair attacks C1, closing S-C1(OH)-C2=C3-C4(CH2C5H11), drawn as 2-hydroxy-5-hexyl-2,5-dihydrothiophene; −H2O gives 2-hexylthiophene (ring C1-C4, hexyl = C5-C10) | 0 / 0 / 184 (290) / 1220 (452) / 436 (315) | mechanism_drawn + within_study_ratio |
| F90-C 2-pentylthiapyran | 2,4-decadienal + H2S -> 2-pentyl-2H-thiapyran | right branch: H2S adds S to C5 and H to C4 (1,6-addition) giving C5H11-CH(SH)-CH2-CH=CH-CHO; S attacks C1, closing S-C1(OH)-C2=C3-C4H2-C5(C5H11), drawn as 2-hydroxy-6-pentyl-3,6-dihydro-2H-thiopyran; −H2O gives the 2H-thiapyran, pentyl on the sp3 carbon next to S | 0 / 35 (33) / 3150 (1880) / 34700 (14800) / 12500 (8250); the largest lipid-dependent product in every phospholipid pot; "the reaction of H2S with 2,4-decadienal may favour the formation of the thiapyran while 2,4-nonadienal and 2,4-octadienal tend to react to give the alkylthiophenes" (2-butyl and 2-propylthiapyran were detected below their isomeric thiophenes; not tabulated) | mechanism_drawn + within_study_ratio |
| F90-D shorter alkylthiophenes | 2,4-nonadienal, 2,4-octadienal + H2S -> 2-pentyl-, 2-butylthiophene | by analogy with B (stated) | 2-pentylthiophene 0 / 0 / 109 / 2700 / 1230; 2-butylthiophene 0 / 0 / 33 / 176 / 176; "the relative amounts of the three alkylthiophenes in the PC runs do not reflect the extreme preponderance of 2,4-decadienal" | within_study_ratio |
| F90-E alkenylthiophenes | lipid + H2S -> 2-(1-hexenyl)thiophene (two isomers) | not drawn | 0 / 0 / 0 / 65 + 374 / 42 + 287; egg phospholipids only | level_only |
| F90-F 1-alkanethiols | 1-heptanol, 1-octanol (lipid) + H2S -> 1-heptanethiol, 1-octanethiol | stated: "formed by the action of H2S on the corresponding alcohols rather than from the aldehydes" (PC made more heptanol / octanol than PE, PE more of the aldehydes) | 0 / 0 / 0 / 456 (351) / 113 (126); 0 / 0 / 0 / 341 (247) / 104 (111) | level_only (mechanism asserted) |
| F90-G thiol quench, class Bb | cysteine + ribose -> MFT, 2-thiophenethiol, 3-thiophenethiol, 2-methyl-3-thiophenethiol, suppressed most by PE and BPL | stated: PUFA (≥ 3 double bonds: BPL 28.3 %, PE 20.7 %, PC 6.2 %, BTG 0.4 %) and their carbonyls compete for H2S and NH3; ethanolamine amino groups (PE 20 mM, BPL 8 mM vs 41 mM cysteine) a secondary candidate | ratios in §3 derived table; MFT 0.40 / 0.15 / 0.27 / 0.24; 2-methyl-3-thiophenethiol 0.076 / 0.0025 / 0.19 / 0.028 | within_study_ratio |
| F90-H FFT quench, class Bd | -> FFT, "equally reduced by all four lipids" | as G | 0.67 / 0.63 / 0.62 / 0.72 | within_study_ratio |
| F90-I mercaptoketone quench, class Ba | -> the four mercaptoketones, reduced slightly by BTG, about halved by every phospholipid | stated: mercaptoketones come from H2S + dicarbonyl or alpha,beta-unsaturated carbonyl (Boelens 1975, Badings 1976, Takken 1976); "if the effect were simply caused by the removal of H2S ... one might have expected PC to exert a lesser effect" | 0.97 / 0.45 / 0.51 / 0.48 (44); 0.72 / 0.49 / 0.53 / 0.53 (45); 0.77 / 0.47 / 0.50 / 0.49 (46) | within_study_ratio |
| F90-J triglyceride-specific suppression, class Bc | bicyclic thiophenes, 3-methyl-1,2,4-trithiane, 2-acetyl- and 2-propionylthiophene fall most with BTG | authors suspect a physical cause: BTG formed a separate layer during headspace collection and these are among the least volatile products | e.g. dihydrothienothiophene (37) 57200 -> 25900 with BTG, unchanged with phospholipids | within_study_ratio (possibly artefact) |
| F90-K unaffected, class C | thiazoles, pyrazines, dithianones, 1-(2-furyl)-2-propanone | pyrazines from aminoketone condensation, "free NH3 ... may not be involved" | differences within SD | within_study_ratio (null) |
| F90-L furfural | ribose -> furfural, halved by PC and PE | polymerisation with lipid carbonyls (stated) | 0.80 / 0.75 / 0.49 / 0.52 | within_study_ratio |
| F90-M thienylethanal | rises with every lipid, most with PC | none | 1150 -> 1990 / 1570 / 3420 / 1650 | level_only |

## 5. Rule sketches (reactant -> product in words; controls are mine, SMILES mine)

**S1. 2,4-decadienal + NH3 -> 2-pentylpyridine + 2 H2O (F90-A; net of carbinolamine, 6-endo
cyclisation, dehydration, oxidation; terminal).** Required substructure: O=CH-CH=CH-CH=CH-C.
- positive: `CCCCC/C=C/C=C/C=O` (`DECADIENAL`) + `N` -> `CCCCCc1ccccn1`; 2,4-heptadienal `CC/C=C/C=C/C=O` + `N` -> 2-ethylpyridine `CCc1ccccn1`.
- negative: (E)-2-nonenal `CCCCCC/C=C/C=O` + `N` -> no fire (one C=C; cannot close a six-ring with N); hexanal + `N` -> no fire (Elmore 1997's trialkylpyridine needs three aldehydes and is a different rule); 2-pentylfuran + `N` -> no fire.
- The N donor is free NH3. Whitfield 1988's amino-acid ranking (Cys 194.7 >> Gly 9.0 > Lys 1.9 ng) is the control that a rule keyed on cysteine's amino group, or on lysine's epsilon-amine, must fail.

**S2. 2,4-alkadienal + H2S -> 2-(alkyl-CH2)-thiophene + H2O (F90-B; 1,4-addition regiochemistry;
terminal).**
- positive: `CCCCC/C=C/C=C/C=O` + `S` -> 2-hexylthiophene `CCCCCCc1cccs1`; 2,4-nonadienal `CCCC/C=C/C=C/C=O` + `S` -> 2-pentylthiophene `CCCCCc1cccs1`.
- negative: `CCCCCC/C=C/C=O` + `S` -> no fire; 2-pentylfuran `CCCCCc1ccco1` + `S` -> must not fire under this rule (the furan -> thiophene exchange is the route Whitfield 1988 argued against and this paper calls less probable).

**S3. 2,4-alkadienal + H2S -> 2-alkyl-2H-thiapyran + H2O (F90-C; 1,6-addition regiochemistry;
terminal). Write S2 and S3 as a pair on the same substructure.**
- positive: `CCCCC/C=C/C=C/C=O` + `S` -> 2-pentyl-2H-thiapyran `CCCCCC1C=CC=CS1`.
- negative: as S2. Branch note: in this pot S3 : S2 = 28 by area for decadienal (PC); `mottram2002b_extraction.md` flag 3 records the opposite in a methyl-linoleate pot; do not encode a branch ratio.

**S4. 1-alkanol + H2S -> 1-alkanethiol + H2O (F90-F; asserted, not shown; status proposed).**
- positive: 1-heptanol `CCCCCCCO` + `S` -> `CCCCCCCS`.
- negative: heptanal `CCCCCCC=O` + `S` -> no fire (the authors' argument is that the thiols track the alcohols, not the aldehydes); MFT + `S` -> no fire.

**S5. lipid quench (a sink, not a structural rule).** Suppression of MFT, the thiophenethiols and
the mercaptoketones scales with the lipid's polyunsaturation (Table 1 Σ ≥ 3 double bonds: BPL 28.3 %,
PE 20.7 %, PC 6.2 %, BTG 0.4 %) for the class Bb thiols, but not monotonically for MFT itself (PC
0.27 with 6.2 % PUFA; PE 0.24 with 20.7 %) or FFT (flat at 0.62-0.72). A single H2S sink
proportional to unsaturated-carbonyl supply reproduces the direction, not the per-compound ordering.

## 6. Flags

1. **Relative peak areas, not amounts.** No response factors; one ion per compound divided by its
   spectral abundance. The within-compound ratios across the five pots are the usable output; the
   ranking across compounds (e.g. FFT 50200 vs MFT 30600) is only approximate and disagrees with
   Whitfield 1988's calibrated ng (FFT 2781 vs MFT 2596, ratio 1.07 vs 1.64 here).
2. **Large SDs on the lipid-dependent products:** 2-pentylthiapyran PC 34700 ± 14800, 2-pentylthiophene
   PC 2700 ± 1970, 2-pentylpyridine PC 5210 ± 2240 (n = 4). The PC : PE ratios are within a factor of
   two of their uncertainty.
3. **Physical partitioning is a live alternative for two classes.** The authors themselves attribute
   the triglyceride-specific suppression (class Bc) to BTG forming a separate layer during the 60 C
   purge; phospholipid binding of hydrophobic thiols is not tested. The thiol ratios are an upper bound
   on the chemical quench, as in Whitfield 1988.
4. **PC made 64x the 2,4-decadienal of PE when heated alone, but only 12x the 2-pentylpyridine and
   2.8x the 2-pentylthiapyran in the Maillard pot**; the products are not linear in the dienal supply.
   Whatever rate the lipid lane assigns to decadienal + NH3 / H2S, the NH3 and H2S supply (from
   cysteine) or their competing sinks limit the product in the PC pot.
5. **The alkylthiophene chain distribution does not follow the dienal distribution** (stated), and the
   2-butyl / 2-propylthiapyrans were below their isomeric thiophenes while 2-pentylthiapyran was 28x
   2-hexylthiophene: the thiophene : thiapyran branch depends on the dienal chain length in a way the
   drawn mechanism does not explain. Keep S2 / S3 without a branch fraction.
6. **Correction of Whitfield 1988:** the 622.2 ng "2-heptylthiophene" is 2-pentylthiapyran (this
   paper, citing Mottram & Salter 1989). The corrected identification rests on MS similarity to
   2-hexylthiophene with a later retention time; no authentic standard.
7. **Tiny headspace (0.5 mL) under N2 in a 2.5 mL ampoule**, yet the phospholipids oxidised to
   dienals, alkanals and alkanols; the oxygen source (pre-formed hydroperoxides in the Sigma lipids,
   dissolved air, thermal routes) is not established. Egg PC in this study had 16.2 % 18:2 and 3.8 %
   20:4.
8. **Aroma paradox recorded by the authors:** the pots with the lowest MFT (BPL, PE) smelled most
   meaty; they suggest MFT above its pleasant range smells pungent and sulphurous. Not a modelling
   number, but it warns against using "more MFT = more meaty" as a target.
9. **pH 5.6 in 0.5 M phosphate with < 0.5 unit drift** is closer to the sulfur lane's reference pots
   than Whitfield 1988's 0.2 M buffer; the two studies' PC columns nevertheless agree (paper's
   statement; the MFT ratios 0.27 vs 0.34 and FFT 0.62 vs 0.50 bear it out).
10. **No time course, one temperature, one lipid level (15 g/L, about 10x an isolate's 1-3 % lipid
    on a dry basis but comparable to meat).** Whether the quench is linear in lipid content is not
    tested anywhere in the group's work on disk.
11. **Compound 44 is named "2-mercapto-3-butanone"** in the paper; it is the same molecule as
    3-mercapto-2-butanone, `CC(S)C(C)=O`, the product of Elmore 1997's acetoin + H2S step.
