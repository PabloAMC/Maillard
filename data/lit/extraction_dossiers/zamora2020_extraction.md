# Zamora, Lavado-Tena & Hidalgo 2020 — EXTRACTION (lipid-derived aldehydes + ammonia sources on silica gel, 180 C / 1 h; calibrated GC-MS yields of 2-, 3-, 2,5- and 2,6-alkylpyridines, D2O labelling)
### The only paper in the corpus that prints a calibrated yield for 2-pentylpyridine from 2,4-decadienal (14.68 µmol per mmol glutamine) and places it in a homologous series of 2,4-alkadienals.

**Source on disk:** `data/articles/zamora2020.pdf` (owner's download, 2026-09-08). Read from the `pypdf`
text layer, which is clean: Tables 1 and 2 came through row by row with their letters. Figures 1-4 are
drawn reaction schemes; the text describes every step they draw and those descriptions are what is
used below (nothing read off the drawings). Fig. S-1 (2-alkenal oligomerisation scheme) and Table S-1
(pyridines reported in foods) are in the Supplementary Material, which was not available.

## 0. Identity

| field | value |
|---|---|
| Title | "Oligomerization of reactive carbonyls in the presence of ammonia-producing compounds: A route for the production of pyridines in foods" |
| Authors | Rosario Zamora, Cristina M. Lavado-Tena, Francisco J. Hidalgo (Instituto de la Grasa, CSIC, Seville) |
| Venue | Food Chemistry 304 (2020) 125284; PII S0308-8146(19)31394-9 |
| DOI | 10.1016/j.foodchem.2019.125284 |
| Naming | Table 2 abbreviations: ACET acetaldehyde; ACR acrolein; BUT butanal; CROT crotonaldehyde (2-butenal); DD 2,4-decadienal; HDO 3,5-heptadien-2-one; HEP 2-heptenal; HEX 2-hexenal; HpD 2,4-heptadienal; HxD 2,4-hexadienal; MHDO 6-methyl-3,5-heptadien-2-one (the Materials section calls it 6-methyl-5-hepten-2-one; see flag 8); ND 2,4-nonadienal; OCT 2-octenal; OD 2,4-octadienal; PENT 2-pentenal; PROP propanal; Ala alanine. Geometry of the aldehydes is not stated (commercial, presumably E / E,E). |
| Companions | Zamora, Navarro, Aguilar & Hidalgo 2015 (Food Chem 174, 89: thermal breakdown of 2,4-alkadienals into 2-alkenals and alkanals; dossier `zamora2015_extraction.md`); Hidalgo, Delgado & Zamora 2009 (the silica-gel support method); Kim & Ho 1998 J. Food Lipids 5, 173-182 (the imine route they cite for 2,4-hexadienal; NOT the JAFC 46, 644 paper in this batch, see flag 6) |

## 1. Why it matters

Programme 7 (roadmap section 5d) needs a rule and, eventually, a rate for 2,4-decadienal + ammonia ->
2-pentylpyridine. This paper heats each of eleven lipid-derived aldehydes (six 2-alkenals, five
2,4-alkadienals) plus acetaldehyde, propanal, butanal, two dienones and alanine with an ammonia source
under one fixed protocol and quantifies the pyridines with authentic standards. It gives (i) the
2-pentylpyridine yield from 2,4-decadienal + glutamine (14.68 ± 0.58 µmol/mmol Gln) and the
3-pentylpyridine yield beside it (8.42 ± 0.53, MS-only), (ii) the same numbers for the C6-C9 dienals so
the chain-length trend is visible (34.07 -> 30.95 -> 27.15 -> 22.42 -> 14.68), (iii) the ranking of
ammonia sources on one aldehyde (ammonia ≈ urea > glutamine ≈ NH4Cl > creatinine; Table 1), (iv) the
demonstration that 2-alkenals alone give 2,5-dialkylpyridines at 10-34 % yield, an order of magnitude
above any 2-alkylpyridine, and (v) D2O labelling that fixes which ring protons are exchangeable. For
the model: glutamine (a species the engine has, `Gln`) is the ammonia source the authors chose as
"a common food component"; ammonia itself is not an engine species.

## 2. Methods as they matter to a model

- **Charge (section 2.2):** reactive carbonyl 50 µmol in 50 µL methanol (two carbonyls: 50 µmol of
  each); ammonia-producing compound 10 µmol in 45 µL water (alanine, when added, in that water);
  0.063-0.20 mm silica gel 300 mg as support; 30 µL of 0.3 mol/L sodium phosphate pH 6.5; 50 µL water.
  So the carbonyl : amine molar ratio is 5 : 1 (10 : 1 for two-carbonyl mixtures) and the "pH 6.5" is
  the nominal buffer pH of a 175 µL liquid phase on 300 mg silica, not a measured pH of the reaction.
- **Heating:** 180 C, closed test tubes, 1 h. Single time point, single temperature.
- **Work-up:** after cooling, 700 µL methanol + 30 µL internal standard (19 µmol methyl heptanoate per
  mL methanol, i.e. 0.57 µmol IS); stir 1 min; centrifuge 5 min at 2000 g; supernatant to GC-MS.
- **Menhaden oil run (section 2.2, last paragraph):** oil 1 g + glutamine 30 µmol in 150 µL water +
  25 µL methanol + 30 µL 0.3 M phosphate pH 6.5; 180 C / 1 h; 700 µL acetonitrile + same IS.
- **GC-MS:** Agilent 7820A / 5977 quadrupole; HP-5MS UI 30 m x 0.25 mm x 0.25 µm; He 1 mL/min;
  pulsed splitless 1 µL; injector 250 C; EI 70 eV; m/z 28-550; 40 C (3 min) -> 200 C at 20 C/min,
  hold 1 min. Fast oven; retention indices in Table 2 are on this column.
- **Quantification (section 2.5):** external standard curves for 2-methyl-, 3-methyl-, 2-ethyl-,
  3-ethyl-, 2-pentyl-, 2,5-dimethyl-, 2,6-dimethyl-, 5-ethyl-2-methyl-, 3-acetyl- and
  3-hydroxy-2-methylpyridine, seven levels 0-10 µmol, each standard added to 300 mg silica and taken
  through the same work-up without heating; content proportional to compound/IS area ratio, r > 0.99,
  RSD < 10 %. Pyridines without a standard (2-propyl, 2-butyl, 3-propyl, 3-butyl, 3-pentyl,
  2,5-diethyl, 2-ethyl-5-propyl, 5-butyl-2-propyl, 2-butyl-5-pentyl, 5-hexyl-2-pentyl) were estimated
  on the curve of "the compound with the closest chemical structure" (stated examples: 2-ethylpyridine
  curve for 2-propylpyridine, 2-pentylpyridine curve for 2-butylpyridine). **2-Pentylpyridine has its
  own authentic standard (Table 2: RI, MS, ST); 3-pentylpyridine is MS-only and borrowed-curve.**
- **Units:** µmol of pyridine per mmol of amine compound charged (Table 1) or per mmol glutamine
  (Table 2). With 10 µmol amine charged, 1 µmol/mmol = 0.01 µmol absolute = 0.1 % of the amine
  compound on a per-molecule basis. Yield on the carbonyl is never printed; it can be computed from
  the charge (below, marked as ours).
- **Replication:** mean ± SD of at least three independent experiments; ANOVA + Tukey, p < 0.05.
- **Labelling (section 2.3):** D2O in place of water; deuterium count in the product pyridine read by
  GC-MS. No 13C or 15N labelling.
- **Identification:** RI + MS + co-elution with standard ("ST") where a standard existed; MS-only
  identifications are declared "only tentative" by the authors (Table 2 footnote).

## 3. Tables re-typed

### Table 1. "Pyridines produced by crotonaldehyde oligomerization in the presence of ammonia and ammonia-producing compounds." µmol per mmol of amine compound; mean ± SD, n ≥ 3; letters = Tukey groups within a column.

Conditions: crotonaldehyde 50 µmol + amine compound 10 µmol, silica gel, pH 6.5 buffer, 180 C, 1 h.

| ammonia-producing compound | 2-methylpyridine | 3-methylpyridine | 2,5-dimethylpyridine | 5-ethyl-2-methylpyridine |
|---|---:|---:|---:|---:|
| ammonia | 22.89 ± 7.01 a | 1.10 ± 0.13 a | 4.44 ± 0.96 a,b | 463.3 ± 81.1 a |
| ammonium chloride | 9.56 ± 1.52 b,c | 4.41 ± 0.98 b | 2.04 ± 0.66 c,d | 258.1 ± 37.4 b |
| glutamine | 14.96 ± 0.33 c | 3.34 ± 0.66 b | 3.66 ± 0.70 b,d | 321.6 ± 11.9 b,d |
| creatinine | 7.59 ± 2.34 b,c | 1.11 ± 0.31 a | 2.09 ± 0.39 c,d | 119.8 ± 6.9 c |
| urea | 27.51 ± 4.10 a | 1.82 ± 0.28 a | 5.70 ± 1.10 a,b | 388.9 ± 59.0 a,d |

Text beside the table: the same four pyridines from every amine; no oligomerisation without an amine;
no part of the amine other than N enters the ring; with ammonia "almost 50 % of the initial nitrogen
was incorporated into the produced pyridines" (row sum 491.7 µmol/mmol, i.e. 49 % per mole of ammonia).

### Table 2. "Pyridines produced by reactive carbonyl oligomerization in the presence of glutamine." µmol per mmol glutamine; mean ± SD, n ≥ 3; letters compare rows within one pyridine.

Conditions: carbonyl(s) 50 µmol each + glutamine 10 µmol, silica gel, pH 6.5 buffer, 180 C, 1 h.
"Identification" is the paper's column: RI, MS, ST (co-elution with standard).

| pyridine | RI (HP-5MS) | identification | amount | carbonyl compound or precursor |
|---|---:|---|---:|---|
| 2-methyl | 814 | RI, MS, ST | 9.93 ± 1.35 a | CROT |
| | | | 34.07 ± 1.33 b | HxD |
| | | | 13.90 ± 2.79 a,c | ACET |
| | | | 16.58 ± 1.04 c | CROT/HxD |
| | | | 1.55 ± 0.26 d | ACET/ACR |
| | | | 50.21 ± 1.66 e | ACET/CROT |
| | | | 2.46 ± 0.59 d | PROP/CROT |
| | | | 1.59 ± 0.51 d | ACR/Ala |
| 3-methyl | 863 | RI, MS, ST | 17.46 ± 1.47 a | ACR |
| | | | 12.06 ± 1.40 a,b | HxD |
| | | | 2.37 ± 0.82 c,d | ACET |
| | | | 3.56 ± 0.46 c,d | PROP |
| | | | 0.71 ± 0.08 c,d | Ala |
| | | | 11.04 ± 3.81 a,b | ACR/CROT |
| | | | 4.92 ± 1.03 b,d | ACR/PENT |
| | | | 7.08 ± 1.05 b,d | ACR/HxD |
| | | | 5.49 ± 1.44 b,d | CROT/HxD |
| | | | 45.05 ± 0.82 e | ACET/ACR |
| | | | 1.40 ± 0.04 c,d | ACET/CROT |
| | | | 67.02 ± 4.18 f | PROP/ACR |
| | | | 1.70 ± 0.36 c,d | PROP/CROT |
| | | | 30.06 ± 0.65 g | ACR/Ala |
| | | | 5.23 ± 1.79 b,d | BUT/ACR |
| 2,6-dimethyl | 889 | RI, MS, ST | 26.93 ± 1.38 a | HDO |
| | | | 3.48 ± 0.47 b | MHDO |
| 2-ethyl | 909 | RI, MS, ST | 30.95 ± 0.34 a | HpD |
| | | | 1.85 ± 0.33 b | ACR/HpD |
| | | | 19.01 ± 0.50 c | CROT/HpD |
| | | | 8.81 ± 1.54 d | ACE/PENT (ACE = ACET) |
| 2,5-dimethyl | 938 | RI, MS, ST | 8.61 ± 1.43 a | ACET |
| | | | 0.75 ± 0.04 b | PROP |
| | | | 5.13 ± 0.56 c | ACR/CROT |
| | | | 3.66 ± 0.39 d | CROT/HpD |
| | | | 1.35 ± 0.47 b | ACET/ACR |
| | | | 7.22 ± 0.31 a | ACET/CROT |
| | | | 82.82 ± 0.75 e | PROP/CROT |
| | | | 0.85 ± 0.08 b | ACR/Ala |
| 3-ethyl | 965 | RI, MS, ST | 17.85 ± 0.59 a | HpD |
| | | | 15.93 ± 4.43 a,b | ACET |
| | | | 4.89 ± 1.69 c,d | ACR/CROT |
| | | | 8.60 ± 0.80 c,e | CROT/HpD |
| | | | 5.91 ± 0.66 c,d | ACET/ACR |
| | | | 1.92 ± 0.08 d | ACET/CROT |
| | | | 2.63 ± 0.31 d | ACR/Ala |
| | | | 10.99 ± 2.05 b,e | BUT/ACR |
| 2-propyl | 1003 | MS | 27.15 ± 1.82 a | OD |
| | | | 15.61 ± 1.04 b | ACE/HEX |
| 5-ethyl-2-methyl | 1034 | RI, MS, ST | 321.6 ± 11.9 a | CROT |
| | | | 32.63 ± 1.51 b,c | ACET |
| | | | 2.55 ± 0.96 d | ACR/CROT |
| | | | 24.33 ± 2.14 b | CROT/PENT |
| | | | 63.77 ± 3.05 e | CROT/HxD |
| | | | 43.58 ± 3.79 c | CROT/HpD |
| | | | 171.5 ± 12.5 f | ACET/CROT |
| | | | 33.32 ± 8.76 b,c | PROP/CROT |
| 3-propyl | 1060 | MS | 12.12 ± 1.38 | OD |
| 2-butyl | 1102 | MS | 22.42 ± 1.03 | ND |
| 2,5-diethyl | 1123 | MS | 19.09 ± 2.13 | CROT/PENT |
| 3-butyl | 1163 | MS | 12.18 ± 0.33 | ND |
| **2-pentyl** | 1205 | RI, MS, ST | **14.68 ± 0.58** | **DD** |
| 2-ethyl-5-propyl | 1211 | MS | 101.9 ± 12.5 a | PENT |
| | | | 1.83 ± 0.25 b | ACR/PENT |
| | | | 29.95 ± 8.44 c | CROT/PENT |
| | | | 1.34 ± 0.06 b | ACR/CROT/PENT |
| **3-pentyl** | 1263 | MS | **8.42 ± 0.53** | **DD** |
| 5-butyl-2-propyl | 1398 | MS | 156.9 ± 18.0 a | HEX |
| | | | 82.24 ± 10.94 b | ACE/HEX |
| 2-butyl-5-pentyl | 1599 | MS | 257.8 ± 9.3 | HEP |
| 5-hexyl-2-pentyl | 1806 | MS | 336.6 ± 43.9 | OCT |

Absences that the table implies (a carbonyl not listed under a pyridine was not reported to give it):
2,4-decadienal appears only under 2-pentyl and 3-pentyl; no 2,5-dialkylpyridine is listed for any
2,4-alkadienal alone; no 2-alkylpyridine is listed for any 2-alkenal alone. The paper does not say
whether "not listed" means not detected or below some limit; there is no LOD/LOQ statement.

### Menhaden oil + glutamine (text, section 3.6; no table)

1 g oil + 30 µmol glutamine, 180 C / 1 h: 3-methylpyridine 34.13 ± 6.61 µmol/mmol amine compound
("the main pyridine detected"), 2-ethylpyridine 4.34 ± 0.45 µmol/mmol. Hexanal appeared in parallel.
No other pyridine amounts printed for the oil.

### Yields on the carbonyl (ours, computed from the 50 µmol / 10 µmol charge; one carbonyl skeleton per 2-alkylpyridine, two per 2,5-dialkylpyridine)

| system | pyridine | µmol/mmol Gln | absolute µmol | % of Gln (per molecule) | % of carbonyl charged |
|---|---|---:|---:|---:|---:|
| DD + Gln | 2-pentylpyridine | 14.68 | 0.147 | 1.47 | 0.29 |
| DD + Gln | 3-pentylpyridine | 8.42 | 0.084 | 0.84 | 0.17 |
| ND + Gln | 2-butylpyridine | 22.42 | 0.224 | 2.24 | 0.45 |
| HxD + Gln | 2-methylpyridine | 34.07 | 0.341 | 3.41 | 0.68 |
| CROT + Gln | 5-ethyl-2-methylpyridine | 321.6 | 3.22 | 32.2 | 12.9 (2 CROT each) |
| OCT + Gln | 5-hexyl-2-pentylpyridine | 336.6 | 3.37 | 33.7 | 13.5 (2 OCT each) |

The authors' own summary of the last two rows: "the observed yields for this kind of pyridines were
10-34 %, and the best yield obtained for 2,5-dimethylpyridine was 8 %" (percent of glutamine).

## 4. Routes and numbers the repository can use

Conditions for every row unless stated: silica gel, nominal pH 6.5, 180 C, 1 h, carbonyl 50 µmol,
amine 10 µmol, n ≥ 3, calibrated GC-MS.

| route | reactant -> product | mechanism as drawn / stated | measured numbers (units, conditions) | evidence class |
|---|---|---|---|---|
| ZA-2ALK-DIENAL | **2,4-decadienal + NH3 (from Gln) -> 2-pentylpyridine** | Fig. 1, "alternative mechanism" for 2,4-alkadienals: aldimine of the dienal with ammonia, then "formation of the pyridine ring" (electrocyclisation of the aza-triene and oxidation are implied, not itemised in the text); attributed to Kim & Ho 1998 (J. Food Lipids) | 14.68 ± 0.58 µmol/mmol Gln (RI, MS, ST); homologues: HxD -> 2-methyl 34.07 ± 1.33; HpD -> 2-ethyl 30.95 ± 0.34; OD -> 2-propyl 27.15 ± 1.82 (MS); ND -> 2-butyl 22.42 ± 1.03 (MS) | measured_yield; mechanism_drawn (imine route) |
| ZA-3ALK-DIENAL | 2,4-decadienal + NH3 -> 3-pentylpyridine | not drawn for the dienal; stated as "a consequence of being the alkadienal the origin of both alkanals and acrolein" (dienal breaks to acrolein + alkanal, then the Fig. 2 route) | 8.42 ± 0.53 µmol/mmol Gln (MS-only, borrowed curve); ND -> 3-butyl 12.18 ± 0.33; OD -> 3-propyl 12.12 ± 1.38; HpD -> 3-ethyl 17.85 ± 0.59; HxD -> 3-methyl 12.06 ± 1.40 | measured_yield (tentative identity); mechanism proposed in words |
| ZA-2ALK-MAIN | crotonaldehyde + acetaldehyde + NH3 -> 2-methylpyridine (the main 2-alkylpyridine route) | Fig. 1: ammonia adds to crotonaldehyde (conjugate addition, 3-aminobutanal), the amine condenses with acetaldehyde to an aldimine, cyclisation, aromatisation; D2O gives a di-deuterated 2-methylpyridine ("deuteration of interchangeable protons during ring formation") | ACET/CROT 50.21 ± 1.66; CROT alone 9.93 ± 1.35; ACET alone 13.90 ± 2.79 µmol/mmol Gln; 2-pentenal/ACET -> 2-ethyl 8.81 ± 1.54; 2-hexenal/ACET -> 2-propyl 15.61 ± 1.04 | measured_yield; mechanism_drawn (D2O-supported) |
| ZA-3ALK-ACR | acrolein + alkanal + NH3 -> 3-alkylpyridine | Fig. 2: ammonia adds to acrolein (3-aminopropanal), imine with propanal, cyclisation; D2O gives mainly mono-deuterated 3-methylpyridine | PROP/ACR -> 3-methyl 67.02 ± 4.18; ACET/ACR 45.05 ± 0.82; ACR/Ala 30.06 ± 0.65; ACR alone 17.46 ± 1.47; BUT/ACR -> 3-ethyl 10.99 ± 2.05 | measured_yield; mechanism_drawn |
| ZA-25-MIX | crotonaldehyde + propanal + NH3 -> 2,5-dimethylpyridine | Fig. 3: ammonia adds to crotonaldehyde, then reacts with propanal, cyclisation, oxidation; D2O -> mono-deuterated | PROP/CROT 82.82 ± 0.75; CROT/PENT -> 2,5-diethyl 19.09 ± 2.13 | measured_yield; mechanism_drawn |
| ZA-25-SELF | **2 x 2-alkenal + NH3 -> 2-(C(n-3))-5-(C(n-2))-dialkylpyridine** (self-oligomerisation, no oxidation step) | Fig. S-1 (not available); text: as Fig. 3 but the extra double bond "facilitates the formation of the extended conjugated system before the cyclization" and "the oxidation step is not required" | CROT -> 5-ethyl-2-methyl 321.6 ± 11.9; PENT -> 2-ethyl-5-propyl 101.9 ± 12.5; HEX -> 5-butyl-2-propyl 156.9 ± 18.0; HEP -> 2-butyl-5-pentyl 257.8 ± 9.3; OCT -> 5-hexyl-2-pentyl 336.6 ± 43.9 µmol/mmol Gln (all but the first MS-only) | measured_yield; mechanism_drawn (in SI) |
| ZA-26-KETONE | 3,5-heptadien-2-one + NH3 -> 2,6-dimethylpyridine | Fig. 4: imine then cyclisation, or ammonia addition to the gamma,delta C=C then imine; oxidation of the cyclic intermediate; D2O -> mono-, di- and some tri-deuterated | HDO 26.93 ± 1.38; MHDO 3.48 ± 0.47 (with loss of the 6-methyl, "fate not investigated") | measured_yield; mechanism_drawn |
| ZA-NSOURCE | crotonaldehyde + {NH3, NH4Cl, Gln, urea, creatinine} -> same four pyridines | Table 1; N source ranking by total pyridine: NH3 491.7 ≈ urea 423.9 > Gln 343.6 > NH4Cl 274.1 > creatinine 130.6 µmol/mmol amine | Table 1 rows | measured_yield; within_study_ratio |
| ZA-OIL | menhaden oil + Gln -> 3-methylpyridine, 2-ethylpyridine | attributed to acrolein and 2,4-heptadienal from omega-3 oxidation | 34.13 ± 6.61 and 4.34 ± 0.45 µmol/mmol Gln (1 g oil, 30 µmol Gln) | measured_yield (level in oil; no carbonyl basis) |

Within-study ratios worth registering:
- Same protocol, same N source: 2-pentylpyridine : 3-pentylpyridine from 2,4-decadienal = 14.68 : 8.42
  ≈ 1.7 : 1 (the 3-isomer is MS-only).
- 2-alkylpyridine from Cn 2,4-alkadienal, n = 6, 7, 8, 9, 10: 34.07, 30.95, 27.15, 22.42, 14.68
  µmol/mmol Gln, a monotonic fall of 2.3x from C6 to C10 (authors: "a decrease ... when the chain
  length of the reactive carbonyl increased"). The C8 and C9 values are on borrowed calibration curves.
- Glutamine vs ammonia as N source on crotonaldehyde: 343.6 / 491.7 = 0.70 of the total pyridine;
  glutamine vs NH4Cl 1.25.
- 2,5-dialkylpyridine from a 2-alkenal alone vs 2-alkylpyridine from the 2,4-alkadienal of the same
  chain: 5-hexyl-2-pentyl (OCT) 336.6 vs 2-propyl (OD) 27.15: 12x on the µmol basis.

## 5. Rule sketches (repository suggestions, not the paper's)

Registry check (`data/keys/compounds.yml`): **2-pentylpyridine has no key; 2-pentylthiophene and
2-hexylthiophene have no key**; `2_pentyl_4_methylthiazole` and `2_hexyl_4_methylthiazole` exist
(not products of this paper); `acrolein`, `hexanal`, `nonanal`, `hydrogen_sulfide` exist. Species:
`DECADIENAL` (`CCCCC/C=C/C=C/C=O`, E,E), `Gln`, `Asn`, `Lys` exist in `data/species/structures.yml`;
**ammonia is not a species**. Every rule below that consumes NH3 therefore needs either an NH3
literature species or a two-step write (Gln -> NH3 + pyroglutamate; then NH3 + dienal), and the
engine would still make no rate.

**S1. 2,4-alkadienal + NH3 -> 2-alkylpyridine (net; ZA-2ALK-DIENAL; the Programme 7 target).**
Reactant: a 2,4-dienal R-CH=CH-CH=CH-CHO. Change: N condenses on C1 (aldimine), N bonds to C5, the
C1-C5 chain plus N becomes the ring, one H2 is lost (aromatisation) and one H2O (imine). Atom map
(from the Zhou 2000 5-13C label, `zhou2000_extraction.md`): dienal C5 -> pyridine C2 (bearing R),
dienal C1 -> pyridine C6.
- positive: `CCCCC/C=C/C=C/C=O` (DECADIENAL) + `N` -> `CCCCCc1ccccn1` (2-pentylpyridine) + `O` (+ H2)
- second positive: 2,4-nonadienal `CCCC/C=C/C=C/C=O` + `N` -> 2-butylpyridine `CCCCc1ccccn1` (Table 2, 22.42; Du 2023 system C 191 µg/L)
- negative: hexanal `CCCCCC=O` + `N` -> no fire (no diene; Table 2 lists no 2-alkylpyridine from any alkanal alone except via acetaldehyde/crotonaldehyde condensation); (E)-2-nonenal `CCCCCC/C=C/C=O` + `N` -> must not give 2-pentylpyridine by this rule (Du 2023 system B did report 2-pentylpyridine, 44 µg/L, but a C9 enal cannot give a C10 pyridine by S1: a different route).
- required substructure: `O=CH-CH=CH-CH=CH-C` with N-H2 on ammonia (not an amine: Table 1 shows the ring N comes from ammonia; an alpha-amino acid N does not enter the ring).

**S2. 2 x 2-alkenal + NH3 -> 2,5-dialkylpyridine (net; ZA-25-SELF).** Reactant: 2 x R-CH=CH-CHO.
Change: ammonia adds to C3 of one enal; the amine condenses with the CHO of the second; cyclisation,
loss of 2 H2O; no oxidation. Product from a Cn 2-alkenal: 2-[C(n-3) alkyl]-5-[C(n-2) alkyl]pyridine
(crotonaldehyde C4 -> 2-methyl-5-ethyl; 2-octenal C8 -> 2-pentyl-5-hexyl).
- positive: 2 x `CCCCC/C=C/C=O` (2-octenal) + `N` -> `CCCCCCc1ccc(CCCCC)nc1` (5-hexyl-2-pentylpyridine; 336.6 µmol/mmol Gln); 2 x `C/C=C/C=O` + `N` -> `CCc1ccc(C)nc1` (5-ethyl-2-methylpyridine; 321.6)
- proposed extension (no source): `DECENAL_2E` (`CCCCCCC/C=C/C=O`, literature_structures) -> 2-heptyl-5-octylpyridine `CCCCCCCCc1ccc(CCCCCCC)nc1`; mark `status: proposed`.
- negative: 2 x hexanal + `N` -> no fire (needs the alpha,beta C=C); DECADIENAL alone -> no 2,5-dialkylpyridine (Table 2 lists none for any dienal alone).

**S3. Acrolein + alkanal + NH3 -> 3-alkylpyridine (net; ZA-3ALK-ACR).**
- positive: `C=CC=O` + `CCC=O` + `N` -> `Cc1cccnc1` (3-methylpyridine; 67.02)
- negative: `CCC=O` + `N` alone -> 3.56 only (needs acrolein); hexanal + `N` -> no fire.
- The 3-pentylpyridine from 2,4-decadienal would be this rule after a dienal -> acrolein + hexanal (?) fragmentation the paper does not draw; leave 3-pentylpyridine as `proposed`.

**S4. Crotonaldehyde + acetaldehyde + NH3 -> 2-methylpyridine (net; ZA-2ALK-MAIN).**
- positive: `C/C=C/C=O` + `CC=O` + `N` -> `Cc1ccccn1` (50.21); negative: `CC=O` + `N` alone is allowed at 13.90 only because acetaldehyde self-condenses to crotonaldehyde; a SMIRKS on the pair should not fire on `CCC=O` + `CC=O`.

**S5. Glutamine -> NH3 (deamidation; the N source; Kim & Ho 1998 Fig. 2 draws intramolecular
cyclisation to pyroglutamic acid).** `NC(CCC(N)=O)C(=O)O` -> `O=C1CCC(N1)C(=O)O` + `N`. Negative: `Glu`
`NC(CCC(=O)O)C(=O)O` -> no fire (no amide). This is the step that lets the rule layer reach S1 from an
engine species.

## 6. Flags

1. **Silica-gel support, 175 µL liquid, 180 C, closed tube.** A dry, high-surface, low-water system;
   the authors' own earlier work (Hidalgo et al. 2009) chose it to mimic low-moisture frying/roasting.
   Not an aqueous buffer and not an oil; yields will not transfer to a 100 C soy slurry (Zhou 2000
   finds the reaction at room temperature in water, at far lower conversion).
2. **Carbonyl : amine = 5 : 1 and the yield is per mmol of amine.** The 0.29 % yield on the dienal
   (ours) is the number a rate fit against DECADIENAL would need; the paper's 14.68 is per Gln.
3. **Half the Table 2 identities are MS-only and on borrowed calibration curves** (every 2,5-dialkyl
   pyridine above C7, 3-pentylpyridine, 2-butylpyridine, 2-propylpyridine). The authors say to treat
   them as tentative. 2-Pentylpyridine itself is fully identified and calibrated.
4. **No time series, no temperature series, no LOD/LOQ, no blank without amine printed as a number**
   (the text says crotonaldehyde did not oligomerise without amine; no table row).
5. **Mechanisms are drawn, D2O-supported, not intermediate-isolated.** The D2O count fixes how many
   exchangeable positions there are, not the order of steps. Evidence class mechanism_drawn.
6. **Citation trap:** the "Kim & Ho, 1998" the paper credits with the dienal-imine route is J. Food
   Lipids 5, 173-182 (glutamine / glutamic acid with a mixture of alkadienals), not JAFC 46, 644-647
   (`kim1998_extraction.md`), which draws no pyridine mechanism.
7. **Nominal pH 6.5** is the buffer's pH; nothing was measured after heating. Glutamine at 10 µmol in
   175 µL is 57 mM; the buffer is 51 mM phosphate; the pH is not controlled in any meaningful sense.
8. **MHDO naming inconsistency:** Materials lists 6-methyl-5-hepten-2-one; the Table 2 footnote
   expands MHDO as 6-methyl-3,5-heptadien-2-one. Only the 2,6-dimethylpyridine row is affected.
9. **No thiazoles, no thiophenes, no H2S** in this paper; it covers the ammonia half of Programme 7
   only. The alkadienal geometry is not stated.
10. **Supplementary Fig. S-1 and Table S-1 not available**; the 2-alkenal oligomerisation scheme is
    described from the text only.
