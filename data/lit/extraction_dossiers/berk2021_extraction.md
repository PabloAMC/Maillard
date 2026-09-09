# Berk, Gürsul Aktağ & Gökmen 2021 — EXTRACTION (whole sesame seeds roasted at 180 / 200 / 220 C; sucrose, bound lysine, five α-dicarbonyls, fructoselysine, CML, CEL; nineteen-step multiresponse model with 95 % HPD and Arrhenius barriers)
### The Kocadağlı-Gökmen network applied to a real, dry, protein-bound-lysine food: CML from fructoselysine (not glyoxal), CEL from methylglyoxal + bound lysine, both with printed barriers.

**Source on disk:** `data/articles/berk2021.pdf` (owner's download, 2026-09-08; 14 pages, Eur. Food Res.
Technol., DOI 10.1007/s00217-021-03787-x, received 28 Mar 2021, accepted 22 May 2021; no volume/page in
the file). Read from the `pdftotext` text layer in the scratchpad; Table 1 is a rotated landscape table whose
text layer is clean and row-consistent (every row carries three k ± HPD pairs in order 180, 200, 220 C;
pypdf layout mode scrambles it and was not used). Table 2 and the Appendix ODEs are clean. Fig. 2 (all
concentration-time data, ten panels) and Fig. 3 (mass balance) are images: **FIGURE-ONLY**; only the axis
ranges and the maxima quoted in the text survive. Supplementary Fig. S1 (the discarded steps' fits) is not
on disk. Repo status before this dossier: not cited anywhere (grep of `data/lit`, `src/`, `results/`).
The same laboratory's glucose-glass constants are already on the trunk (`parameters_dicarbonyl.py`, B13,
from `kocadagli2016jafc_extraction.md` section 4); this paper is the same modelling machinery with a protein.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of α-dicarbonyl compounds and glycation products in sesame (Sesamum indicum L.) seeds during roasting: a multiresponse kinetic modelling approach" |
| Authors | Ecem Berk, Işıl Gürsul Aktağ, Vural Gökmen (corresponding; FoQuS group, Hacettepe University, Ankara) |
| Venue | European Food Research and Technology, 2021 (Springer); Original Paper |
| Abbreviations (theirs) | SUC sucrose; GLU glucose; INT intermediates = GLU + FFC; FFC fructofuranosyl cation; bLys protein-bound lysine; FL N-ε-fructoselysine (= furosine x 2.2); HP Heyns product (not measured); 3-DG 3-deoxyglucosone; 1-DG 1-deoxyglucosone; GO glyoxal; MG methylglyoxal; DA diacetyl; CML, CEL; P1..P3 unmeasured products |
| Lineage | Kocadağlı & Gökmen 2016 (Food Chem 211:892, glucose/wheat flour) and Göncüoğlu Taş & Gökmen 2017 (hazelnut roasting) networks; Athena Visual Studio 14.2; determinant criterion; HPD intervals |

## 1. Why it matters

Programme 7 needs protein-bound lysine as a slow reactant with CML and CEL as markers; the repository's
dicarbonyl trio came from an amine-free glucose glass and its docstring assumes glyoxal makes CML and
methylglyoxal makes CEL. This paper is the nearest thing on disk to a measured answer for a dry,
protein-bound system: sucrose + bound lysine in a ~20 % protein, ~50 % oil seed, 180-220 C, with all
five dicarbonyls, FL, CML and CEL measured and fitted together. Its results, in the authors' words: CEL
"originated from the reaction between methylglyoxal and bound lysine"; CML "predominantly by the oxidation
of N-ε-fructoselysine compared to the reaction of glyoxal with bound lysine" (the GO + bLys -> CML constant
fitted to 0 at all three temperatures); FL is the main 1-DG source; 1-DG fragmentation to MG and DA is
"the fastest one"; 3-DG and GO come mainly from glucose. Barriers are printed for FL -> CML (113 kJ/mol),
MG + bLys -> CEL (92), 1-DG -> MG (67) and sucrose -> intermediates (137), with R^2 of the Arrhenius line.
Nguyen 2016 (aqueous casein, 120-130 C) agrees on CML-from-AP and disagrees on CEL (AP route there).

Limits up front: a real food with unmeasured moisture and water activity; several constants indeterminate
or with HPD wider than the estimate; five negative activation energies; no CML or CEL sink steps; the
time courses are figure-only.

## 2. Methods as they matter to a model

- **Matrix:** Nigerian sesame seeds, whole, bought in Ankara. Composition from the introduction (generic):
  ~50 % oil, ~20 % protein, ~5 % sugar, unsaturated fatty acids 80-85 % of the oil, lignan antioxidants.
  Measured: **sucrose 22.5 ± 1.9 g/kg** = **65.7 mmol/kg** (342.30 g/mol; the only detectable sugar);
  **protein-bound lysine 10.4 ± 0.3 g/kg** = **71.1 mmol/kg** (146.19 g/mol) — i.e. **sucrose : lysine
  sites ≈ 0.92 : 1**, and per gram of protein (at 20 %) ≈ 0.36 mmol lysine/g. Raw seed already held 3-DG
  and MG (attributed to storage/transport); 1-DG, GO, DA absent in raw seed.
- **Moisture / water activity: NOT MEASURED, NOT PRINTED.** The text says only "owing to low moisture
  content, sesame seeds reach high temperatures immediately during roasting" and that sucrose "degrades
  into glucose and FFC rather than hydrolyzes". Raw sesame is typically 4-6 % moisture; treat a_w as
  unknown-low, falling during the roast.
- **Roasting:** 10 g seeds as a 3 mm layer in a petri dish, Memmert UN 55 oven at **180, 200, 220 C**,
  "time intervals ranging between 2.5 and 30 min"; the Fig. 2/3 axes give the grids: 180 C: 0, 5, 10, 20,
  30 min (Fig. 3a; Fig. 2 axes run to 30); 200 C: 0, 5, 10, 15, 20 min; 220 C: 0, 2.5, 5, 7.5, 10 min.
  Air oven, open dish, unstirred. Triplicate roasts; **individual replicate values used in the fit** (not
  averages). Cooled to ambient, ground, stored at −18 C.
- **Concentration unit:** µmol per kg of roasted seed (mass loss on roasting not corrected for or stated).
- **Sugars:** triple water extraction (10 + 5 + 5 mL per g), Carrez, Oasis HLB, HPLC-RI (Shodex KC-811,
  5 mM H2SO4); external calibration 0.001-1.0 %.
- **Protein-bound lysine:** 100 mg seed + 5 mL 8 N HCl, N2 headspace, 110 C 23 h; 50 µL hydrolysate dried,
  1 mL water; LC-MS/MS (Waters TQD, HILIC, ammonium formate / formic acid); theanine internal standard at
  0.5 mg/L; external calibration 0.25-10 mg/mL. Not isotope dilution.
- **α-Dicarbonyls (3-DG, 1-DG, MG, GO, DA):** 200 µL aqueous extract + 800 µL acetonitrile:water 5:3;
  500 µL supernatant + 150 µL 0.2 % o-phenylenediamine with 11 mM DETAPAC + 150 µL 0.5 M phosphate pH 7;
  2 h dark, room temperature; LC-MS/MS SIM of the quinoxalines (m/z 235.2 for 3-DG and 1-DG, 145 MG, 131
  GO, 159.2 DA); **5-methylquinoxaline internal standard**; calibration 0.02-2.0 mg/L of quinoxaline
  standards; glucosone and 3-DG working solutions derivatised alongside. Not isotope dilution; 1-DG is
  quantified against the 3-DG-derived response (no 1-DG standard is listed). Unit conversions used in
  section 3: 3-DG 162.14, 1-DG 162.14, MG 72.06, GO 58.04, DA 86.09 g/mol.
- **Furosine:** 200 µL hydrolysate, Oasis HLB, HPLC-DAD 280 nm (Atlantis HILIC, 1 % formic acid);
  external calibration 1-10 mg/L. **FL = furosine x 2.2** (Krause 2003). FL 308.33 g/mol.
- **CML and CEL:** 20 mg seed + 100 µL water + 450 µL 0.2 M borate pH 9.2 + 500 µL 1 M NaBH4 in 0.1 M
  NaOH, 4 h room temperature (**reduction: fructoselysine -> hexitol-lysine**, so AP cannot convert to CML
  in the hydrolysis); + 2 mL 8 N HCl, N2, 110 C 24 h; 20 µL hydrolysate dried, 1 mL water, Oasis HLB;
  LC-MS/MS MRM (Atlantis T3, 0.1 % formic acid / acetonitrile), m/z 205.10 -> 84 and 130 (CML), 219 -> 84
  and 130 (CEL). **Quantification by matrix-matched calibration** (0, 1, 2.5, 5, 10, 20 µg/mL of CML and
  CEL spiked into "blank" sesame — seeds roasted at 150 C — and carried through the whole procedure). **No
  internal standard is named for CML/CEL** — this is NOT isotope dilution. The blank matrix is itself
  heat-treated and will hold some CML/CEL; the calibration's zero point absorbs that.
- **Modelling:** Athena Visual Studio 14.2, numerical integration, non-linear regression on the determinant
  criterion (van Boekel 1996), per-temperature fits with 95 % HPD; Arrhenius line through the three k(T)
  gives Ea with an R^2 (Table 2) — **no HPD on Ea**, and the fit is on three points.

### The scheme in words (Fig. 1a, the proposed model; Fig. 1b is the comprehensive one; ODEs in the Appendix)

Nineteen steps; the ODEs make explicit what Table 1's labels blur (GLU and FFC are both "INT" in the
equations):
1. **SUC -> INT (k1)**, INT = glucose + fructofuranosyl cation, i.e. thermal (not hydrolytic) cleavage.
2. **INT -> GO (k2)** (Table label GLU -> GO); 4. **INT -> 3-DG (k4)** (GLU -> 3-DG).
3. **INT + bLys -> FL (k3)** (second order; GLU + bLys); 5. **INT + bLys -> HP (k5)** (FFC + bLys -> Heyns
   product, unmeasured).
6. **FL -> GO (k6)**; 7. **FL -> P1 (k7)**; 8. **FL -> CML (k8)**; 9. **FL -> 1-DG (k9)**. Steps 6, 7, 9
   return bLys to the pool (eq. 9); step 8 does not (CML keeps the lysine).
10. **GO + bLys -> CML (k10)**, second order — fitted 0 at all temperatures, kept in the model.
11. **HP -> 3-DG (k11)** (fixed at 0.10 x 10^-3 at 180 and 220 C); 13. **HP -> 1-DG (k13)**; 16. **HP -> CEL
    (k16)** (0 at all T); 19. **HP -> P3 (k19)** (0 at all T). Steps 11, 13, 19 return bLys.
12. **3-DG -> MG (k12)**; 14. **1-DG -> MG (k14)**; 17. **1-DG -> DA (k17)**; 18. **DA -> P2 (k18)**.
15. **MG + bLys -> CEL (k15)**, second order.
Absent by construction: any sink for 3-DG other than -> MG, any sink for GO other than -> CML (which is
zero, so GO only accumulates), any sink for MG other than -> CEL, any CML or CEL degradation ("omitted ...
because kinetic model did not produce good fits"), any glucose -> fructose, any FL <-> HP interconversion,
and the P4/P6/P7 dicarbonyl sinks and P1/P5 glucose/FFC sinks of Fig. 1b (all removed in discrimination).

## 3. Tables re-typed

### Table 1. "Reaction rate constants with their 95 % highest posterior density (HPD) intervals calculated for each elementary reaction step of the kinetic model shown in Fig. 1a"

Units as printed: first-order steps **min^-1 x 10^3**; second-order steps (3, 5, 10, 15) **kg µmol^-1 min^-1
x 10^3**. "ind*: Indeterminate, which means a large uncertainty in the estimated parameter within 95 %
confidence interval." "fixed" = held during estimation.

| # | step | unit | 180 C k | ± HPD | 200 C k | ± HPD | 220 C k | ± HPD |
|---|---|---|---:|---:|---:|---:|---:|---:|
| 1 | SUC -> INT | min^-1 x 10^3 | 9.29 | 1.12 | 61.6 | 4.7 | 174.8 | 11.7 |
| 2 | GLU -> GO | min^-1 x 10^3 | 0 | 0 | 0 | 0 | 857 | 226 |
| 3 | GLU + bLys -> FL | kg µmol^-1 min^-1 x 10^3 | 10.1 | ind | 0.48 | 0.06 | 0.75 | 0.33 |
| 4 | GLU -> 3-DG | min^-1 x 10^3 | 1842 | 882 | 518 | 663 | 125 | 244 |
| 5 | FFC + bLys -> HP | kg µmol^-1 min^-1 x 10^3 | 12.2 | 21.7 | 16.2 | ind | 14.0 | ind |
| 6 | FL -> GO | min^-1 x 10^3 | 0.20 | 0.08 | 17.3 | 3.10 | 0 | 0 |
| 7 | FL -> P1 | min^-1 x 10^3 | 721 | 762 | 141 | 94 | 1023 | 507 |
| 8 | FL -> CML | min^-1 x 10^3 | 5.54 | 0.88 | 29.0 | 0.90 | 62.2 | 7.84 |
| 9 | FL -> 1-DG | min^-1 x 10^3 | 58.5 | 9.88 | 107.1 | 87.8 | 102.6 | 18.5 |
| 10 | GO + bLys -> CML | kg µmol^-1 min^-1 x 10^3 | 0 | 0 | 0 | 0 | 0 | 0 |
| 11 | HP -> 3-DG | min^-1 x 10^3 | 0.10 | fixed | 0.17 | 0.43 | 0.10 | fixed |
| 12 | 3-DG -> MG | min^-1 x 10^3 | 0 | 0 | 222 | 622 | 131 | 107 |
| 13 | HP -> 1-DG | min^-1 x 10^3 | 0.64 | 0.68 | 0.45 | 1.11 | 0.87 | 0.11 |
| 14 | 1-DG -> MG | min^-1 x 10^3 | 3011 | 478 | 4795 | 3176 | 12939 | 1425 |
| 15 | MG + bLys -> CEL | kg µmol^-1 min^-1 x 10^3 | 0.0034 | 0.0008 | 0.01 | 0.001 | 0.024 | 0.003 |
| 16 | HP -> CEL | min^-1 x 10^3 | 0 | 0 | 0 | 0 | 0 | 0 |
| 17 | 1-DG -> DA | min^-1 x 10^3 | 3428 | ind | 2759 | 8093 | 342 | 77 |
| 18 | DA -> P2 | min^-1 x 10^3 | 1827 | 459 | 1137 | 3397 | 0 | 0 |
| 19 | HP -> P3 | min^-1 x 10^3 | 0 | 0 | 0 | 0 | 0 | 0 |

Rows whose HPD equals or exceeds the estimate (interval spans zero) at a given temperature: 4 (200, 220 C),
5 (180 C, plus ind at 200/220), 7 (180 C), 11 (200 C), 12 (200 C), 13 (180, 200 C), 17 (200 C, ind at
180), 18 (200 C). Rows determinate at all three temperatures: **1, 8, 9 (marginal at 200 C: 107.1 ± 87.8),
14, 15**.

Second-order unit conversion: 1 kg µmol^-1 min^-1 x 10^-3 = 1e-3 kg µmol^-1 min^-1 = 1e3 kg mol^-1 min^-1.
So k15 = 0.0034e-3 -> **3.4 kg mol^-1 min^-1** (180 C), 10 (200 C), 24 (220 C); k3 at 200 C = 0.48e-3 ->
480 kg mol^-1 min^-1. "kg" is kg of seed, not of water; the constant is not a molarity constant (flag 5).

### Table 2. "Activation energies (Ea) and coefficient of determination (R^2) calculated for each elementary reaction step of the kinetic model shown in Fig. 1a"

| # | step | Ea (kJ/mol) | R^2 |
|---|---|---:|---:|
| 1 | SUC -> INT | 137 | 0.980 |
| 2 | GLU -> GO | – | – |
| 3 | GLU + bLys -> FL | −122 | 0.648 |
| 4 | GLU -> 3-DG | −125 | 0.997 |
| 5 | FFC + bLys -> HP | 7 | 0.252 |
| 6 | FL -> GO | 399 | 1 |
| 7 | FL -> P1 | 14 | 0.020 |
| 8 | FL -> CML | 113 | 0.966 |
| 9 | FL -> 1-DG | 27 | 0.716 |
| 10 | GO + bLys -> CML | – | – |
| 11 | HP -> 3-DG | 0.7 | 0.0006 |
| 12 | 3-DG -> MG | −51 | 1 |
| 13 | HP -> 1-DG | 14 | 0.208 |
| 14 | 1-DG -> MG | 67 | 0.948 |
| 15 | MG + bLys -> CEL | 92 | 0.999 |
| 16 | HP -> CEL | – | – |
| 17 | 1-DG -> DA | −106 | 0.801 |
| 18 | DA -> P2 | −42 | 1 |
| 19 | HP -> P3 | – | – |

No uncertainty on any Ea. R^2 = 1 for rows 6, 12, 18 means a two-point line (one temperature had k = 0).
My two-point check (180 vs 220 C, R = 8.3145e-3, 1/453.15 − 1/493.15 = 1.790e-4 K^-1): k1 136.3, k8
112.3, k14 67.7, k15 90.8 kJ/mol — all agree with Table 2 to rounding, so the printed Ea are honest
Arrhenius lines through Table 1. Authors on the negative values: "due to the accumulation of intermediates,
these reactions acquire free energy, and in this way, they can progress without energy barrier" (citing
Vyazovkin 2016) — read as: those constants are compensating for missing steps, not as chemistry.

### Concentrations printed in the text (Fig. 2 is figure-only)

Per kg roasted seed; mg/kg as printed, µmol/kg by me.

| species | value | µmol/kg | condition | source |
|---|---|---:|---|---|
| sucrose, raw | 22.5 ± 1.9 g/kg | 65,700 | t = 0 | Results |
| sucrose loss at end of roast | 22 % / 80 % / 76 % | — | 180 C 30 min / 200 C 20 min / 220 C 10 min | Results |
| bound lysine, raw | 10.4 ± 0.3 g/kg | 71,100 | t = 0 | Results |
| 3-DG maximum | 6.9 ± 0.4 mg/kg | 42.6 | 180 C, 30 min ("quantitatively predominant") | Results |
| MG maximum | 5.7 ± 0.9 mg/kg | 79.1 | 220 C, 10 min | Results |
| GO maximum | 4.0 ± 0.2 mg/kg | 68.9 | 220 C, 10 min; GO "not quantified at the beginning", then "a sudden increase after prolonged roasting at all temperatures" | Results |
| 1-DG at end | 0.9 ± 0.1 mg/kg | 5.6 | 200 C, 20 min; "did not exhibit a particular change" | Results |
| DA at end | 1.2 ± 0.2 mg/kg | 13.9 | 200 C, 20 min; same | Results |
| FL maximum | 110.4 ± 32.6 / 64.9 ± 2.0 mg/kg | 358 / 210 | 180 C / 200 C, at 5 min | Results |
| FL maximum | 97.7 ± 2.0 mg/kg | 317 | 220 C, at 2.5 min | Results |
| CML "apparent maximum" | 28.2 ± 4.6 mg/kg | 138 | 220 C, 7.5 min (rising pattern at all T) | Results |
| CEL "apparent maximum" | 88.6 ± 10.9 mg/kg | 406 | 220 C, 7.5 min | Results |
| mass balance at 5 min | 92.1 / 83.2 / 71.3 % | — | 180 / 200 / 220 C | Results |
| mass balance at end | 82.1 / 59.2 / 61.3 % | — | 180 / 200 / 220 C | Results |

Fig. 2 axis ranges (µmol/kg): sucrose 0-70,000; bLys 0-100,000; FL 0-600; 3-DG 0-50; 1-DG 0-8; GO 0-80;
MG 0-80; DA 0-20; CML 0-160; CEL 0-500 — consistent with the maxima above. CEL exceeds CML about 3:1 at
220 C; CML 138 µmol/kg is 0.19 % of lysine sites; CEL 406 µmol/kg is 0.57 %; FL peak 358 µmol/kg is 0.5 %.
The dicarbonyls at their maxima are 40-80 µmol/kg, i.e. ~0.1 % of the sucrose carbon — small pools against
which k14 (3-13 min^-1) and k17 (0.3-3 min^-1) make 1-DG a transient of seconds.

## 4. Kinetic numbers the repository can use

Registry keys: CML -> `cml`; CEL -> `cel`; furosine -> `furosine`; diacetyl -> `2_3_butanedione`;
glyoxal, methylglyoxal, 3-deoxyglucosone, 1-deoxyglucosone, fructoselysine, protein-bound lysine, sucrose:
**not in registry**. Conditions for every row: whole sesame seed (~20 % protein, ~50 % oil), sucrose
65.7 mmol/kg, bLys 71.1 mmol/kg, open dish in an air oven, moisture and a_w unmeasured (low), no buffer,
no pH; 180 / 200 / 220 C; times to 30 / 20 / 10 min.

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| FL -> CML (oxidation of the Amadori product) | k8 | 5.54 ± 0.88 / 29.0 ± 0.90 / 62.2 ± 7.84 | 10^-3 min^-1 | 180 / 200 / 220 C | Table 1 | measured_rate (determinate at all T) |
| FL -> CML | Ea | 113 (R^2 0.966; my two-point 112) | kJ/mol | 180-220 C | Table 2 | measured_barrier (3 points, no interval) |
| MG + bLys -> CEL | k15 | 0.0034 ± 0.0008 / 0.01 ± 0.001 / 0.024 ± 0.003 (= 3.4 / 10 / 24 kg mol^-1 min^-1) | 10^-3 kg µmol^-1 min^-1 | 180 / 200 / 220 C | Table 1 | measured_rate (per kg seed, see flag 5) |
| MG + bLys -> CEL | Ea | 92 (R^2 0.999; my two-point 91) | kJ/mol | 180-220 C | Table 2 | measured_barrier |
| GO + bLys -> CML | k10 | 0 ± 0 at all T | kg µmol^-1 min^-1 | 180-220 C | Table 1 | within_study_ratio: GO route / FL route -> 0 (structural) |
| HP -> CEL | k16 | 0 ± 0 at all T | min^-1 | | Table 1 | structural zero |
| 1-DG -> MG | k14 | 3011 ± 478 / 4795 ± 3176 / 12939 ± 1425 | 10^-3 min^-1 | 180 / 200 / 220 C | Table 1 | measured_rate |
| 1-DG -> MG | Ea | 67 (R^2 0.948) | kJ/mol | | Table 2 | measured_barrier |
| 1-DG -> DA | k17 | ind / 2759 ± 8093 / 342 ± 77 | 10^-3 min^-1 | | Table 1 | measured_rate at 220 C only; Ea −106 unusable |
| DA -> P2 | k18 | 1827 ± 459 / 1137 ± 3397 / 0 | 10^-3 min^-1 | | Table 1 | measured_rate at 180 C only |
| FL -> 1-DG | k9 | 58.5 ± 9.88 / 107.1 ± 87.8 / 102.6 ± 18.5 | 10^-3 min^-1 | | Table 1 | measured_rate; Ea 27 (R^2 0.716) weak |
| FL -> GO | k6 | 0.20 ± 0.08 / 17.3 ± 3.10 / 0 | 10^-3 min^-1 | | Table 1 | unstable across T (flag 4); Ea 399 meaningless |
| GLU (INT) -> GO | k2 | 0 / 0 / 857 ± 226 | 10^-3 min^-1 | | Table 1 | 220 C only |
| GLU (INT) -> 3-DG | k4 | 1842 ± 882 / 518 ± 663 / 125 ± 244 | 10^-3 min^-1 | | Table 1 | 180 C only; Ea −125 unusable |
| 3-DG -> MG | k12 | 0 / 222 ± 622 / 131 ± 107 | 10^-3 min^-1 | | Table 1 | 220 C marginal; Ea −51 unusable |
| SUC -> INT | k1 | 9.29 ± 1.12 / 61.6 ± 4.7 / 174.8 ± 11.7 | 10^-3 min^-1 | | Table 1 | measured_rate |
| SUC -> INT | Ea | 137 (R^2 0.980) | kJ/mol | | Table 2 | measured_barrier |
| INT + bLys -> FL | k3 | ind / 0.48 ± 0.06 / 0.75 ± 0.33 | 10^-3 kg µmol^-1 min^-1 (= 480 / 750 kg mol^-1 min^-1) | | Table 1 | 200 C only reliable; Ea −122 unusable |
| FL -> P1 (other AP loss) | k7 | 721 ± 762 / 141 ± 94 / 1023 ± 507 | 10^-3 min^-1 | | Table 1 | 200 C only; non-monotone |
| share of FL loss that becomes CML | k8 / (k6 + k7 + k8 + k9) | 0.007 / 0.10 / 0.052 | — | 180 / 200 / 220 C (my arithmetic) | Table 1 | within_study_ratio (derived; k7 dominates and is ill-determined) |
| CEL : CML at 220 C, 7.5 min | 88.6 : 28.2 mg/kg = 406 : 138 µmol/kg | 2.9 | — | 220 C | Results | within_study_ratio |
| CML, CEL, FL, dicarbonyl, sucrose, bLys maxima | section 3 table | mg/kg and µmol/kg | | Results | level_only |
| all ten time courses at three temperatures | — | µmol/kg vs min | | Fig. 2 | figure_only |

Comparison with the trunk's Kocadağlı 2016 JAFC glucose-glass constants for the nominally same steps
(k_b at 180 C, 10^-3 min^-1; from `kocadagli2016jafc_extraction.md` section 4): 1-DG -> DA 12.2 there vs
3428 (ind) here; 3-DG -> MGO 304 there vs 0 here; Glc -> 3-DG 4.19 there vs 1842 here; GO sink 32.6 there
vs none here. Two to three decades apart in either direction for the same lab and machinery: these
"elementary" constants are matrix-specific lumps, and neither set transfers to a protein isolate in water
without a fit of its own. Recorded as a warning, not adjudicated.

## 5. Flags

1. **Moisture and water activity are not measured or printed.** For a repository that keys a_w, this paper
   is "dry roast, low and falling a_w, unknown". The 220 C series lasts 10 min; the 180 C series 30 min.
2. **Five negative activation energies** (k3, k4, k12, k17, k18) and two > 300 (k6 399). The authors'
   explanation is not chemistry; these are constants absorbing missing steps (no 3-DG, GO or MG sinks other
   than the AGE steps, no CML/CEL sinks, no glucose -> fructose, no FL <-> HP). Only k1, k8, k14, k15 have
   monotone k(T) with R^2 > 0.94 — those four barriers are the usable ones.
3. **k3 (glucose + bLys -> FL) is indeterminate at 180 C and falls 13-fold from 180 to 200 C.** The
   authors read the fall as consistent with CML formation from FL rising; more plausibly the FL pool at 5 min
   (the first point) is already past its maximum at 200 and 220 C, so the formation constant is set by one
   or two points. The Amadori-formation constant of this paper should not be adopted.
4. **The glyoxal route flips with temperature**: GO from FL (k6) at 180 and 200 C, from glucose (k2) at
   220 C, and each is zero where the other is active; the authors admit "the model of the formation of GO was
   not well fitted because of the limited data points in GO" (GO was below quantification early in every
   roast). GO + bLys -> CML = 0 is therefore a statement made with a poorly determined GO pool.
5. **Second-order constants are per kg of seed.** MG + bLys -> CEL = 3.4-24 kg mol^-1 min^-1 (180-220 C).
   To compare with an aqueous constant one would divide by the water fraction in which the reactants
   actually meet — unknown (a few per cent, falling). Read the Ea (92 kJ/mol) as transferable in principle;
   read the pre-factor as matrix-bound.
6. **No CML or CEL degradation steps**, although Nguyen 2016 needed them at 120-130 C in water and the
   sesame data end at "apparent maxima" at 220 C / 7.5 min. Extrapolating this model past 10 min at 220 C
   will overshoot.
7. **CML/CEL quantified by matrix-matched external calibration with no internal standard**, in a blank that
   was itself roasted at 150 C. Not isotope dilution; the absolute CML/CEL scale carries an unquantified
   recovery. The dicarbonyls carry a quinoxaline (non-isotopic) internal standard; 1-DG has no standard of
   its own. Bound lysine is external calibration against theanine.
8. **Mass balance falls to 59-61 % at 200 and 220 C** — 40 % of the counted moles are in unmeasured
   products (melanoidins, oil-derived carbonyls). Lipid-derived GO/MG is dismissed on the strength of the
   lignan antioxidants, not measured.
9. **FL = furosine x 2.2** (Nguyen: 3.1). FL levels and every FL-sourced constant depend on it.
10. **Individual replicate values were fitted** (a point in the paper's favour over Nguyen's averages); the
    HPDs already include roast-to-roast scatter.
11. **Table 1's 1-DG -> MG at 220 C (12.9 min^-1) and 1-DG -> DA (0.34 min^-1)** make 1-DG's lifetime ~5 s;
    with a 2.5-min sampling grid the 1-DG pool is at steady state and the two constants are only identified
    through the MG and DA pools. The ratio k14 / k17 (0.9 at 180 C if k17's ind value is taken, 1.7 at
    200 C, 38 at 220 C) is itself unstable.
12. **The CEL-from-MG finding contradicts Nguyen 2016's CEL-from-AP** (aqueous casein); the authors note
    "differences in physical state of reactants and the environment". The repository has one dry and one
    wet answer and should not choose without the isolate's own data.
