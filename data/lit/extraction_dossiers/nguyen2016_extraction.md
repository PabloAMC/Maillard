# Nguyen, van der Fels-Klerx & van Boekel 2016 — EXTRACTION (sodium caseinate + glucose or lactose, 0.1 M phosphate pH 6.8, 120 / 130 C, 0-30 min; CML and CEL by isotope dilution; eleven-step multiresponse model)
### The only aqueous protein-bound-lysine CML/CEL kinetic study in the corpus: Amadori product -> CML and -> CEL constants at sterilisation temperatures, with CML and CEL decay steps.

**Source on disk:** `data/articles/nguyen2016.pdf` (owner's download, 2026-09-08; 9 pages, Food Chemistry 192
(2016) 125-133). Read from the `pdftotext` text layer in the scratchpad; Table 1 is a rotated table whose
text layer is clean (pypdf reports "rotated text", `pdftotext -layout` on p. 7 re-produces every row and was
used to confirm the column order below). Figures 2 and 3 (all concentration-time data) are images: their
axis ranges survive in the text layer, their points do not. The supplementary material (Fig. S1 pH and mass
balance, Fig. S2 differential equations, Table S1 pairwise comparisons) was NOT on disk when this dossier
was written; it is now (`data/articles/nguyen2016_supplementary.docx`, read 2026-10-29: see the section
"Supplementary material" at the end, which supersedes flag 6 and the pH / mass-balance gaps). Repo status before
this dossier: not cited anywhere in `data/lit`, `src/`, `results/` (checked by grep); Quan 2020 (this
paper's methodological descendant) is already declared REFUSE in `k3_final_parameter_inventory.md` C.5.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of Nε-(carboxymethyl)lysine formation in aqueous model systems of sugars and casein" |
| Authors | Ha T. Nguyen (RIKILT Wageningen UR), H.J. van der Fels-Klerx (RIKILT, corresponding), Martinus A.J.S. van Boekel (Food Quality and Design, Wageningen University) |
| Venue | Food Chemistry 192 (2016) 125-133; received 17 Mar 2015, accepted 30 Jun 2015, online 2 Jul 2015 |
| DOI / PII | 10.1016/j.foodchem.2015.06.110 / S0308-8146(15)01002-X |
| Naming | "AP" = Amadori product (N-ε-fructosyl- or lactulosyl-lysine), reported as furosine x 3.1; "MRPs" = advanced Maillard reaction products (unmeasured lump); "LA" = Lobry de Bruyn-Alberda van Ekenstein isomerisation; X1..X4 = unmeasured sinks |
| Lineage | the casein-sugar network of Brands & van Boekel 2001 (JAFC 49:4667), 2002 (JAFC 50:6725) and 2003 (Food Chem 83:13) with CML and CEL bolted on; analytics from Troise et al. 2014 (Amino Acids 46:279) and Delatour et al. 2009 (J Chromatogr A 1216:2371) |

## 1. Why it matters

Programme 7 (`tasks/roadmap_for_scientists.md` section 5d) wants the isolate's protein-bound lysine as a
slow Maillard reactant with CML and CEL as the measured markers; the panel's CML and CEL rows are refused
today because the engine has neither species. This is the one paper on disk that fits protein-bound lysine
(casein, 30 g/L) + a reducing sugar in WATER at 120-130 C to a multiresponse model with explicit steps
sugar + lysine -> AP -> CML and AP -> CEL, and — unusually — with CML and CEL DECAY steps, because both
markers fell after 25 min. Its structural finding is the one the repository needs to decide on: the
glyoxal route to CML (sugar -> GO -> CML) fitted to ZERO; CML came from the Amadori product. That is the
opposite of the topology the B13 dicarbonyl trio was declared for ("glyoxal is the CML precursor",
`parameters_dicarbonyl.py`). Berk 2021 in a dry roast (see `berk2021_extraction.md`) reaches the same
conclusion for CML and the opposite one for CEL (methylglyoxal + lysine, not AP).

Limits stated up front: two temperatures only, no activation energies (the authors could not estimate
them; they give Q10 = k130/k120 instead), all concentration-time data in figures, several constants with
95 % HPD wider than the estimate.

## 2. Methods as they matter to a model

- **Buffer:** 0.1 M sodium phosphate (NaH2PO4/Na2HPO4), pH 6.8, water. pH fell during heating (more at
  130 C; M2 slightly below M1 at 130 C) — Fig. S1A, not on disk, no numbers in the main text.
- **M1 (glucose-casein):** 30 g sodium caseinate + 27 g D-glucose per litre buffer. Glucose 27 g / 180.16 =
  **149.9 mmol/L** (Fig. 2a axis runs 0-160).
- **M2 (lactose-casein):** 30 g sodium caseinate + 51.3 g lactose monohydrate per litre. 51.3 g / 360.31 =
  **142.4 mmol/L** (Fig. 3a axis 0-160; the authors note the measured t = 0 lactose was ~21 % below the
  weighed-in value).
- **Lysine residues:** NOT printed as a number. The design statement is "molar ratio of sugar to lysine
  residues of about 10:1" -> **~15 mmol/L**; casein at 8.2 g lysine / 100 g protein x 30 g/L / 146.19 =
  16.8 mmol/L, and the Fig. 2c / 3c axis runs 0-25 mmol/L. Take **15-17 mmol/L lysine residues, 30 g/L
  protein (sodium caseinate)**, i.e. **0.50-0.56 mmol lysine per g caseinate**.
- **Vessel / heating:** 6 mL aliquots in screw-capped glass tubes (1.2 x 10 cm) with a 1.2 x 3.5 cm
  headspace (≈ 4 mL air; sealed, air atmosphere, not stirred), aluminium heating blocks at 120 or 130 C
  for 30 min; heat-up 5 min (temperature stable afterwards); one tube removed every 5 min (0, 5, 10, 15,
  20, 25, 30) and cooled on ice 1 h; a little brown precipitate in M1 at 30 min (120 C) and 25 min
  (130 C) and in M2 at 25 and 30 min (130 C). Whole experiment duplicated; every analysis in duplicate.
  The model was fitted "including the heating-up time" (Fig. 2 caption), as if isothermal.
- **Sugars:** ethanol 1:1 protein precipitation; HPLC-ELSD, Alltech Prevail Carbohydrate ES, external
  calibration 25-250 µg/mL; LOD 12.5 µg/mL (galactose 250 µg/mL — galactose and tagatose were therefore
  NOT detected in M2 although Brands & van Boekel 2003 found ~8 and ~2.5 mmol/L at 30 min).
- **Lysine, furosine, CML, CEL (one hydrolysate, one injection):** 100 µL sample + 4 mL 6 N HCl, N2-saturated,
  110 C 24 h; PVDF filtration; 400 µL dried under N2; reconstituted in 370 µL water + 10 µL internal
  standards **d4-lysine, d2-CML, d4-CEL at 45.45 ng per mg protein**; Oasis HLB SPE; 10 µL onto
  Kinetex C18 2.6 µm, 5 mM perfluoropentanoic acid ion pairing, API 2000 triple quadrupole, ESI+, MRM.
  Quantifier transitions m/z 205 -> 84.1 (CML), 219.2 -> 84.0 (CEL), 255.1 -> 130.2 (furosine), 147.2 ->
  130.2 (lysine); CML against d2-CML (207 -> 144.1), CEL against d4-CEL (223 -> 134.1), furosine AND
  lysine against d4-lysine (151.2 -> 134.1). **Stable-isotope dilution for CML, CEL, lysine; furosine
  carries a surrogate standard.** Calibration 5 x 10^-3 to 1 µg/mL (lysine, CML, CEL) and 9 x 10^-3 to
  1 µg/mL (furosine); LOD 0.5, 3.0, 0.5, 1.0 x 10^-3 µg/mL for lysine, furosine, CML, CEL. **No NaBH4
  reduction** before hydrolysis (same as Troise 2015): fructosyl-lysine surviving in the sample can
  convert partly to CML during acid hydrolysis; the N2 saturation is the authors' guard.
- **Amadori product** = furosine x **3.1** (Brands & van Boekel 2001). Berk 2021 uses 2.2 (Krause 2003); a
  repository that compares AP levels across the two must undo the factors.
- **Reported units:** everything modelled in mmol/L. The CML level range is also printed per protein:
  M1 16.9-92.2, M2 16.3-103.7 mg CML / 100 g protein. Conversion at 30 g protein/L, CML 204.22 g/mol:
  92.2 mg/100 g -> 27.7 mg/L -> **0.135 mmol/L** (= 4.5 µmol per g protein); 103.7 -> **0.152 mmol/L**;
  16.9 -> 0.025 mmol/L. These match the Fig. 2e / 3e axes (0-0.16 mmol/L). CEL axis 0-0.12 mmol/L
  (CEL 218.25 g/mol: 0.12 mmol/L = 26 mg/L = 87 mg/100 g protein).
- **Time-zero background:** CEL at t = 0 was 12-13 % of the maximum reached, CML at t = 0 was 19 % of the
  maximum (both systems, both temperatures); "AP, CEL and CML were already present" in the caseinate.
- **Mass balance:** ~104-114 % at intermediate times (measurement scatter), falling to ~85 % (M1) and
  ~95 % (M2) at 30 min; the missing 15 % / 5 % is organic acids, melanoidins, and in M2 galactose/tagatose.
- **Fitting:** Athena Visual Studio 14.2, numerical integration, simultaneous non-linear regression on the
  averages of the duplicate analyses; 95 % HPD intervals. Two temperatures: "simultaneous estimation of
  rate constants and activation energies was not well possible. Estimates of Ea values were so low that it
  seemed as if there was no temperature dependence. Therefore the parameter estimation was performed as if
  it was an isothermal process." Q10 = k(130)/k(120) is the only temperature statement.

### The scheme in words (Fig. 5, the fitted model; Fig. 4 is the fuller candidate)

(a) M1 sugar branch: **glucose -> fructose (k1)**, fructose -> glucose (k4; fitted 0 in M1, dropped);
**fructose -> X1 (k2)**, X1 = organic acids (formic, acetic, lactic) and other unmeasured fructose products.
(b) M2 sugar branch: **lactose -> lactulose (k1)**, lactulose -> lactose (k4, minor, 0 at 120 C);
**lactulose -> X2 (k2)**, X2 = galactose, formic acid, C5/C6 products.
(c) Maillard branch, both systems: **sugar + lysine residue -> AP (k3, second order, L mmol^-1 min^-1)**;
**AP -> MRPs (k8)**; **AP -> CML (k7)**; **AP -> CEL (k9)** — the authors read this as AP degrading to
methylglyoxal (with lysine release) which then forms CEL, lumped into one first-order step; **CML -> X3
(k11)**; **CEL -> X4 (k10)**. Fig. 4's oxidative route (sugar -> GO under [O], k5; GO + lysine -> CML, k6)
was fitted: "the rate constant for the formation of CML via the oxidation steps of glucose or lactose with
GO formation as an intermediate was estimated as 0" and it is absent from Table 1 and Fig. 5. Likewise
"the formation of MGO via the oxidation of reducing sugars can be considered negligible".

Rate laws: the differential equations are in Fig. S2A/B (not on disk). All k except k3 are printed in
min^-1, so every other step is first order in its source species. Whether the AP -> CML / CEL / MRPs steps
return lysine to the pool is not recoverable from the main text (see flag 6).

## 3. Tables re-typed

### Table 1. "Rate constants of glucose-casein system (M1) and lactose-casein system (M2) heated at 120 C and 130 C."

Printed as value ± 95 % HPD half-width, with the relative HPD in parentheses. "—: Not determined."
Units as printed: k3 in min^-1 L mmol^-1; all others min^-1.

| constant (step) | unit | M1 120 C | M1 130 C | Q10 M1 | M2 120 C | M2 130 C | Q10 M2 |
|---|---|---|---|---:|---|---|---:|
| k1 (glucose -> fructose / lactose -> lactulose) | min^-1 | 7.4e-3 ± 1.1e-4 (15 %)* | 1.7e-2 ± 1.5e-3 (9 %) | 2.3 | 8.1e-3 ± 1.8e-3 (22 %) | 2.1e-2 ± 2.5e-3 (12 %) | 2.6 |
| k2 (fructose -> X1 / lactulose -> X2) | min^-1 | 1.3e-2 ± 1.1e-2 (90 %) | 2.3e-2 ± 6.3e-3 (28 %) | 1.8 | 2.0e-2 ± 1.6e-2 (80 %) | 1.7e-2 ± 9.4e-3 (55 %) | 0.9 |
| k3 (sugar + lysine -> AP) | min^-1 L mmol^-1 | 1.5e-4 ± 3.0e-5 (20 %) | 1.6e-4 ± 2.2e-5 (13 %) | 1.1 | 1.7e-4 ± 3.4e-5 (20 %) | 1.6e-4 ± 2.7e-5 (17 %) | 0.9 |
| k4 (ketose -> aldose, back-isomerisation) | min^-1 | 0 | 0 | — | 0 | 1.8e-2 ± 1.6e-2 (88 %) | — |
| k7 (AP -> CML) | min^-1 | 8.8e-3 ± 6.6e-3 (75 %) | 6.0e-3 ± 1.6e-3 (27 %) | 0.9 | 8.9e-3 ± 5.0e-3 (56 %) | 7.5e-3 ± 2.2e-3 (29 %) | 0.8 |
| k8 (AP -> MRPs) | min^-1 | 5.2e-2 ± 3.3e-2 (64 %) | 1.5e-1 ± 3.4e-2 (22 %) | 2.9 | 5.6e-2 ± 3.6e-2 (64 %) | 1.1e-1 ± 2.9e-2 (26 %) | 2.0 |
| k9 (AP -> CEL, via MGO) | min^-1 | 2.3e-3 ± 2.7e-3 (117 %) | 2.0e-3 ± 3.0e-4 (16 %) | 0.9 | 1.9e-3 ± 2.9e-3 (153 %) | 2.3e-3 ± 7.6e-4 (33 %) | 1.2 |
| k10 (CEL -> X4) | min^-1 | 1.9e-1 ± 2.7e-1 (143 %) | 0 | — | 1.5e-1 ± 2.9e-1 (196 %) | 2.8e-2 ± 2.8e-2 (99 %) | 0.2 |
| k11 (CML -> X3) | min^-1 | 2.9e-1 ± 2.7e-1 (93 %) | 7.7e-2 ± 3.3e-2 (43 %) | 0.3 | 2.4e-1 ± 2.0e-1 (82 %) | 8.1e-2 ± 3.7e-2 (46 %) | 0.3 |

\* printed "7.4 x 10^-3 ± 1.1 x 10^-4 (15 %)"; 1.1e-4 / 7.4e-3 = 1.5 %, not 15 %. Either the half-width
is 1.1e-3 (then 15 % holds) or the percentage is wrong. Every other cell's percentage agrees with its
value ± half-width to rounding (checked row by row). The (%) column was NOT re-derived for any other cell
beyond that check. k5 and k6 (the GO route of Fig. 4) are not in the table: fitted to zero and dropped.

### Concentration-time data

There is no numeric concentration table. Fig. 2 (M1) and Fig. 3 (M2) each carry six panels — sugar,
isomer (fructose or lactulose), lysine, AP, CML, CEL — at both temperatures with duplicate points and error
bars. **FIGURE-ONLY.** Table S1 (pairwise comparisons between systems and temperatures) is in the
supplement, not on disk. Numbers that ARE printed in the text:

| quantity | value | where |
|---|---|---|
| CML range in M1, all times and temperatures | 16.9-92.2 mg / 100 g protein (0.025-0.135 mmol/L) | Results 3.2 |
| CML range in M2 | 16.3-103.7 mg / 100 g protein (0.024-0.152 mmol/L) | Results 3.2 |
| CML(t = 0) / CML(max) | 19 % (both systems, both temperatures) | Results 3.2 |
| CEL(t = 0) / CEL(max) | 12 % (M1), 13 % (M2) | Results 3.2 |
| CML and CEL time course | rise "continuously until 25 min", then "decreased slightly" | Results 3.2 |
| AP time course | rises, falls after 20 min (120 C) / after 10 min (130 C); AP at 120 C higher than at 130 C; AP comparable between M1 and M2 | Results 3.2 |
| M2 vs M1 CML at 120 C | comparable to 15 min, then M2 higher by 0.4-17 % | Results 3.3 |
| Relative loss at 120 C, 30 min (M1): glucose, lysine; relative formation fructose, AP | "0.25, 0.45, 0.2 and 4.5" (stated as comparable to Brands & van Boekel 2001) | Results 3.2 |
| Relative loss at 120 C (M2): lactose, lysine; formation lactulose, AP | "0.3, 0.45, 0.2 and 4.0" | Results 3.2 |
| Mass balance at 30 min | ~85 % (M1), ~95 % (M2) | Results 3.2 |

The "relative" numbers above are quoted as printed; the paper does not define whether 0.45 for lysine is a
fraction lost (45 %) — read with the Brands 2001 convention (fraction of initial). With lysine ~15 mmol/L
that is ~6.8 mmol/L of lysine residues blocked in 30 min at 120 C; with AP "4.5" relative formation the
AP peak is a few mmol/L (Fig. 2d / 3d axes run 0-7 mmol/L), consistent.

## 4. Kinetic numbers the repository can use

Registry keys: CML -> `cml`; CEL -> `cel`; furosine -> `furosine`; glyoxal, methylglyoxal,
3-deoxyglucosone, fructosyl-lysine (AP), lysine residues: **not in registry** (`reactive_lysine` is a marker
set, not a species). All conditions: 0.1 M phosphate pH 6.8 (drifting down), 30 g/L sodium caseinate
(~15-17 mmol/L lysine residues), sealed tube with air headspace, unstirred, a_w ≈ 1.

| step | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| glucose + lysine residue -> AP | k3 | 1.5e-4 ± 3.0e-5 / 1.6e-4 ± 2.2e-5 | L mmol^-1 min^-1 (= 0.15 / 0.16 L mol^-1 min^-1) | 120 / 130 C, 150 mmol/L glucose | Table 1, M1 | measured_rate |
| lactose + lysine residue -> AP | k3 | 1.7e-4 ± 3.4e-5 / 1.6e-4 ± 2.7e-5 | L mmol^-1 min^-1 | 120 / 130 C, 142 mmol/L lactose | Table 1, M2 | measured_rate |
| AP -> CML | k7 | 8.8e-3 ± 6.6e-3 / 6.0e-3 ± 1.6e-3 | min^-1 | 120 / 130 C, M1 | Table 1 | measured_rate (120 C HPD 75 %) |
| AP -> CML | k7 | 8.9e-3 ± 5.0e-3 / 7.5e-3 ± 2.2e-3 | min^-1 | 120 / 130 C, M2 | Table 1 | measured_rate |
| AP -> CEL (via MGO, lumped) | k9 | 2.3e-3 ± 2.7e-3 / 2.0e-3 ± 3.0e-4 | min^-1 | 120 / 130 C, M1 | Table 1 | measured_rate at 130 C; 120 C interval spans zero |
| AP -> CEL | k9 | 1.9e-3 ± 2.9e-3 / 2.3e-3 ± 7.6e-4 | min^-1 | 120 / 130 C, M2 | Table 1 | as above |
| AP -> MRPs (all other AP loss) | k8 | 5.2e-2 ± 3.3e-2 / 1.5e-1 ± 3.4e-2 | min^-1 | 120 / 130 C, M1 | Table 1 | measured_rate |
| AP -> MRPs | k8 | 5.6e-2 ± 3.6e-2 / 1.1e-1 ± 2.9e-2 | min^-1 | 120 / 130 C, M2 | Table 1 | measured_rate |
| CML -> X3 (CML loss) | k11 | 2.9e-1 ± 2.7e-1 / 7.7e-2 ± 3.3e-2 | min^-1 | 120 / 130 C, M1 | Table 1 | measured_rate, wide; see flag 3 |
| CML -> X3 | k11 | 2.4e-1 ± 2.0e-1 / 8.1e-2 ± 3.7e-2 | min^-1 | 120 / 130 C, M2 | Table 1 | as above |
| CEL -> X4 (CEL loss) | k10 | 1.9e-1 ± 2.7e-1 / 0 | min^-1 | 120 / 130 C, M1 | Table 1 | interval spans zero / fitted zero |
| CEL -> X4 | k10 | 1.5e-1 ± 2.9e-1 / 2.8e-2 ± 2.8e-2 | min^-1 | 120 / 130 C, M2 | Table 1 | interval spans zero / touches zero |
| glucose -> fructose | k1 | 7.4e-3 / 1.7e-2 | min^-1 | 120 / 130 C | Table 1, M1 | measured_rate (see footnote on the 120 C HPD) |
| lactose -> lactulose | k1 | 8.1e-3 ± 1.8e-3 / 2.1e-2 ± 2.5e-3 | min^-1 | 120 / 130 C | Table 1, M2 | measured_rate |
| fructose -> X1 / lactulose -> X2 | k2 | 1.3e-2 / 2.3e-2 (M1); 2.0e-2 / 1.7e-2 (M2) | min^-1 | 120 / 130 C | Table 1 | measured_rate, HPD 28-90 % |
| sugar -> GO -> CML (oxidative route) | k5, k6 | 0 | — | both systems, both T | Results 3.3 | within_study_ratio: GO route / AP route -> 0 |
| CML-forming share of AP loss | k7 / (k7 + k8 + k9) | 0.14 / 0.038 (M1, 120 / 130 C); 0.13 / 0.063 (M2) | — | derived by me from Table 1 | Table 1 | within_study_ratio (derived) |
| CEL-forming share of AP loss | k9 / (k7 + k8 + k9) | 0.036 / 0.013 (M1); 0.028 / 0.019 (M2) | — | derived by me | Table 1 | within_study_ratio (derived) |
| k7 / k11 (sets the CML plateau: CML_ss ≈ (k7/k11) x AP) | 0.030 / 0.078 (M1); 0.037 / 0.093 (M2) | — | 120 / 130 C | derived by me | Table 1 | within_study_ratio (derived; see flag 3) |
| Q10 (k130 / k120) | k1 2.3-2.6; k8 2.0-2.9; k3 0.9-1.1; k7 0.8-0.9; k9 0.9-1.2; k11 0.3 | — | 120 -> 130 C | Table 1 | measured_barrier (two-point, weak): Q10 2.3 <=> Ea ≈ 100 kJ/mol for k1 |
| CML level | 0.025-0.135 (M1), 0.024-0.152 (M2) | mmol/L (16.9-92.2 / 16.3-103.7 mg per 100 g protein) | 0-30 min, 120-130 C | Results 3.2 | level_only |
| CML, CEL, AP, lysine, sugar time courses | — | mmol/L vs min | both systems, both T | Figs 2, 3 | figure_only |
| Temperature dependence of CML formation | none resolvable | — | 120 vs 130 C | Table 1, Results 3.2 ("CML comparable at 120 and 130 C") | level_only / structural |

Q10-to-Ea arithmetic (mine): Ea = R ln(Q10) / (1/393.15 − 1/403.15) = 131.8 kJ/mol x ln(Q10); Q10 2.3 ->
110 kJ/mol; Q10 2.9 -> 140; Q10 0.9 -> −14; Q10 0.3 -> −159. Only k1 and k8 are in the range chemistry
allows; the authors say so themselves ("the other reaction steps have much lower values of Q10 ... these
low rate constants are apparent rate constants comprising more than one reaction").

## 5. Flags

1. **Two temperatures, no barriers.** The authors could not estimate Ea and say the estimates "seemed as if
   there was no temperature dependence". The repository cannot take a barrier for AP -> CML or AP -> CEL
   from this paper; a Q10 of 0.8-1.2 on k3, k7, k9 across 120 -> 130 C is physically a fitting artefact
   (correlation with the sinks and with the heat-up), not a measurement of zero activation energy.
2. **All concentration-time data are figure-only.** Nothing but ranges and the t = 0 fractions is printed.
   If Programme 7 needs the curves, the supplement (Table S1) or a digitisation with a declared tolerance
   would be needed; the house rule forbids reading values off the figure.
3. **k7 and k11 are correlated, and k11 is enormous.** k11 = 0.29 min^-1 at 120 C would halve free CML in
   2.4 min, yet CML rose for 25 min; the fit holds because AP (~5 mmol/L) x k7 feeds it as fast as k11
   removes it (CML_ss ≈ k7 AP / k11 ≈ 0.15 mmol/L, which is the observed plateau). What the data pin is the
   ratio k7/k11 and the late-time decline, not the two constants separately (HPD 93 % and 75 % at 120 C).
   Any adoption should carry the pair, or the ratio, not k7 alone. The same holds for k9/k10 (CEL).
4. **The CML "instability" is a claim about total (protein-bound + free) CML after acid hydrolysis in a
   casein solution with browning precipitate at 25-30 min.** The authors themselves ask for a study "focusing
   on CML alone at high temperatures". The 30-min precipitate (M1 at 30 min / 120 C and 25 min / 130 C; M2
   at 25 and 30 min / 130 C) coincides with the decline; a recovery loss into an unhydrolysable pellet is
   not excluded by the text.
5. **No NaBH4 reduction before hydrolysis.** Fructosyl-lysine surviving to the hydrolysis step can convert
   to CML in hot HCl; the N2 saturation reduces but does not remove this. The t = 0 CML (19 % of max) and
   the AP-tracking shape of CML are both consistent with a small analytical contribution. Berk 2021 does
   reduce; Troise 2015 (the method source) explicitly chose not to. Comparisons across the three papers
   inherit this.
6. **Lysine regeneration is unrecoverable.** Fig. S2 (the ODEs) is not on disk. In the Brands network AP
   degradation returns the amino group; the Nguyen text says of the CEL step "MGO formation (and release of
   lysine)". Whether k7, k8, k9 return lysine changes the lysine balance by up to the AP flux. Do not
   transcribe the network without the supplement.
7. **Lysine-residue concentration is not printed**; 15-17 mmol/L is inferred (section 2). k3 is second
   order and does not need it, but any pseudo-first-order comparison with Berk 2021's kg µmol^-1 min^-1
   constants or with free-lysine pots does.
8. **The GO route fitted to zero is a model-discrimination result under these conditions** (sealed tube,
   ~4 mL air over 6 mL, no added oxidant, 120-130 C, pH 6.8). It says the AP route dominated, not that
   glyoxal + lysine is slow; glyoxal was not measured. Quan 2020 (same descent, free lysine, pH 7.0)
   measures glyoxal at 0.05-0.6 mmol/L and fits CML from it — the two papers are not directly comparable.
9. **Q10 for k10 (CEL loss) is 0.2 and k10 = 0 at 130 C in M1** — the CEL sink is not identified; the
   authors acknowledge "rather high" uncertainty. Treat CEL loss as unmeasured here.
10. **Table 1 k1 (M1, 120 C) footnote**: printed half-width 1.1e-4 is inconsistent with the printed 15 %;
    one of the two is a typo (section 3).
11. **Lactose vs glucose**: k3 equal within HPD in M1 and M2; the higher M2 CML late in the run is
    attributed by the authors to galactose released from lactulose (not detected — the column was blind to
    it), so the M2 network is missing a measured reactant. Within-study M2/M1 CML ratios (1.004-1.17 at
    120 C) are the safest carry-over.
12. **Mass balance falls to ~85 % (M1) at 30 min** and exceeds 100 % (to 114 %) in between; the model was
    fitted to duplicate averages only. Expect the constants to move if refitted to the individual replicates.

## Supplementary material (read 2026-10-29)

**Source on disk:** `data/articles/nguyen2016_supplementary.docx` (Word file; read only,
not modified). Read three ways in the scratchpad, never executing anything from it: `textutil -convert txt`
(captions and Table S1); the OOXML `word/document.xml` (the equations of Fig. S2 are Office Math objects that
`textutil` drops, linearised from the XML); and the figures. **There is no `word/charts/` part: the docx
embeds no chart data.** Every figure is a Windows Metafile (`word/media/image1-6.wmf`, each present twice:
the drawing and its compatibility fallback, byte-different but giving identical numbers) wrapping an
Enhanced Metafile written by Excel. The EMF draws every marker and every stacked-area vertex as vector
geometry, so the plotted values were **recovered from the vector coordinates** (axis calibrated on the drawn
tick marks and their text labels; scripts in the scratchpad), not digitised from pixels. Precision: pH one
drawing unit = 0.0005 pH, tick-calibration residual ≤ 0.0002 pH, and every recovered pH lands within 0.0005
of a multiple of 0.01, i.e. the values below are the plotted numbers to the meter's 0.01 resolution. Mass
balance: vertex quantum 0.006 percentage points, calibration residual ≤ 0.02 points. Times recovered within
0.006 min of 0, 5, ..., 30 min. Panel identity from the figure group's layout (image left = (a), right = (b);
the S1B grid (a) top-left, (b) top-right, (c) bottom-left, (d) bottom-right) and confirmed by content (M2
panels carry series labels 7 and 8; M2 at 130 C ends lower in pH than M1, as the main text says).

What the supplement holds (it has no methods text of its own; conditions and units come from the main paper):

| item | what it is | numbers in the file? |
|---|---|---|
| Fig. S1A | pH vs heating time, M1 (a) and M2 (b), 120 and 130 C, two points per time | yes, vector-recovered (below) |
| Fig. S1B | molar mass balance, stacked areas of each measured species as % of the t = 0 molar total, four panels | yes, vector-recovered (below) |
| Fig. S2A / S2B | NOT figures: the ODE systems for M1 and M2 as typed equations | yes, transcribed exactly (below) |
| Table S1 | paired t-tests of AP, CML, CEL between systems / temperatures | yes, typed table (below) |

### Conditions (quoted) and measurement basis

- Fig. S1A caption: "pH profile of M1 (a) and M2 (b) during heating at 120oC ([glyph]) and 130oC ([glyph])". The two
  legend glyphs are font symbols that do not survive as text (Wingdings 2 char 0xDA; Wingdings char 0x6C,
  the filled circle). In the drawing the
  120 C series is Excel's star marker (blue fill) and the 130 C series a filled black circle; the assignment
  rests on the Wingdings 0x6C = filled circle mapping and agrees with the main text (larger drop at 130 C).
- pH method (main paper 2.3): measured "at room temperature by a pH meter" after the one-hour ice cooling
  (Conductivity Proline Plus, QiS). So these are **cooled-sample, room-temperature readings**, not hot pH.
  Buffer 0.1 M sodium phosphate, set to pH 6.8; sealed 6 mL tubes; the 5 min heat-up is inside the time axis.
- Two markers per time point and no error bars. The paper says "The full experiment ... was replicated"
  and every analysis was duplicated; which of the two the pair of points is, is not stated. No SD is printed;
  the two values are given, and SD(n = 2) = |difference| / √2 if wanted.
- Fig. S1B caption: "Mass balance of reactants and reaction products in M1 heated at 120oC (a), M1 heated at
  130oC (b) and in M2 heated 120oC (c), M2 heated at 130oC (d)". Series (1) glucose, (2) fructose, (3) lysine
  residues, (4) Amadori product, (5) CML, (6) CEL, (7) lactose, (8) lactulose. The caption prints the CML
  band as 0.02-0.06 % (M1, both T), 0.02-0.09 % (M2 120 C), 0.02-0.11 % (M2 130 C) and the CEL band as
  0.01-0.04 % (M1 120), 0.01-0.06 % (M1 130), 0.02-0.05 % (M2 120), 0.01-0.06 % (M2 130).
- Mass-balance basis (main paper 2.6.2): "the percentage of total molar concentrations of reactants and
  reaction products ... at each time point to their concentrations at t = 0 min". Units: % (axis "Percentage
  (%)", 50-110, so the glucose/lactose band's base is clipped; the stacked tops are cumulative and exact).
  Sugars by HPLC-ELSD (external calibration), lysine / CML / CEL by LC-MS/MS stable-isotope dilution, AP =
  furosine × 3.1 (surrogate internal standard). One value per time (the paper fitted duplicate averages;
  the figure presumably shows them, not stated). AP is counted once per mole although it carries both a sugar
  and a lysine moiety, so this is a recovery of measured moles, not an element balance.

### Fig. S1A — pH vs time (vector-recovered; the two duplicate values, room-temperature readings)

| t (min) | M1 glucose 120 C | M1 glucose 130 C | M2 lactose 120 C | M2 lactose 130 C |
|---:|---|---|---|---|
| 0 | 6.82, 6.83 | 6.80, 6.81 | 6.84, 6.84 | 6.84, 6.84 |
| 5 | 6.82, 6.83 | 6.79, 6.79 | 6.83, 6.83 | 6.81, 6.81 |
| 10 | 6.79, 6.79 | 6.71, 6.71 | 6.80, 6.81 | 6.72, 6.73 |
| 15 | 6.75, 6.76 | 6.62, 6.62 | 6.77, 6.78 | 6.62, 6.62 |
| 20 | 6.72, 6.72 | 6.52, 6.54 | 6.73, 6.74 | 6.50, 6.51 |
| 25 | 6.68, 6.68 | 6.47, 6.47 | 6.69, 6.70 | 6.40, 6.42 |
| 30 | 6.65, 6.67 | 6.44, 6.50 | 6.65, 6.70 | 6.36, 6.39 |

Drop over 30 min (means): M1 120 C −0.17 (6.825 -> 6.66); M1 130 C −0.34 (6.805 -> 6.47); M2 120 C −0.17
(6.84 -> 6.675); M2 130 C −0.47 (6.84 -> 6.375). The drift is slow for the first 5 min (heat-up), then near
linear at 130 C (5 -> 20 min: about −0.017 pH/min in M1, −0.020 in M2; 120 C about −0.007). t = 0 is 6.80-6.84, not exactly 6.8.

### Fig. S1B — mass balance, % of the t = 0 molar total (vector-recovered from the stacked areas)

Bands = differences of consecutive stacked tops. Sugars and lysine to 0.1 point; AP, CML, CEL to 0.01 point
(quantum 0.006). "total" = top of the stack = the paper's recovery.

(a) M1 glucose-casein, 120 C

| t (min) | glucose | fructose | lysine | AP | CML | CEL | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 86.3 | 0.0 | 13.3 | 0.34 | 0.02 | 0.01 | 100.0 |
| 5 | 87.3 | 1.5 | 12.5 | 2.19 | 0.03 | 0.01 | 103.5 |
| 10 | 87.2 | 6.7 | 10.5 | 3.32 | 0.05 | 0.02 | 107.7 |
| 15 | 73.7 | 7.8 | 11.7 | 3.45 | 0.06 | 0.02 | 96.8 |
| 20 | 69.3 | 10.7 | 10.7 | 3.28 | 0.07 | 0.03 | 94.1 |
| 25 | 67.3 | 12.8 | 7.1 | 1.51 | 0.04 | 0.04 | 88.8 |
| 30 | 63.6 | 12.9 | 7.4 | 1.49 | 0.06 | 0.02 | 85.5 |

(b) M1 glucose-casein, 130 C

| t (min) | glucose | fructose | lysine | AP | CML | CEL | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 86.3 | 0.0 | 13.1 | 0.34 | 0.02 | 0.01 | 99.7 |
| 5 | 83.7 | 4.5 | 10.9 | 0.56 | 0.03 | 0.01 | 99.8 |
| 10 | 70.4 | 12.5 | 11.7 | 1.86 | 0.06 | 0.01 | 96.5 |
| 15 | 61.8 | 17.5 | 10.2 | 0.91 | 0.05 | 0.04 | 90.4 |
| 20 | 54.8 | 20.1 | 9.4 | 1.03 | 0.08 | 0.06 | 85.5 |
| 25 | 51.0 | 21.4 | 8.5 | 0.76 | 0.08 | 0.06 | 81.8 |
| 30 | 53.8 | 20.3 | 6.7 | 0.92 | 0.07 | 0.05 | 81.8 |

(c) M2 lactose-casein, 120 C

| t (min) | lactose | lactulose | lysine | AP | CML | CEL | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 84.6 | 0.0 | 15.0 | 0.46 | 0.02 | 0.02 | 100.0 |
| 5 | 80.5 | 1.6 | 13.7 | 2.48 | 0.04 | 0.03 | 98.4 |
| 10 | 78.4 | 5.3 | 13.1 | 3.65 | 0.06 | 0.02 | 100.5 |
| 15 | 86.4 | 9.3 | 11.0 | 3.76 | 0.08 | 0.02 | 110.6 |
| 20 | 80.9 | 12.5 | 13.3 | 4.17 | 0.09 | 0.03 | 111.0 |
| 25 | 76.1 | 14.2 | 8.5 | 1.38 | 0.09 | 0.05 | 100.4 |
| 30 | 72.5 | 13.8 | 7.8 | 1.55 | 0.09 | 0.01 | 95.8 |

(d) M2 lactose-casein, 130 C

| t (min) | lactose | lactulose | lysine | AP | CML | CEL | total |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 84.8 | 0.0 | 14.7 | 0.44 | 0.02 | 0.01 | 100.0 |
| 5 | 82.7 | 5.1 | 15.3 | 0.67 | 0.04 | 0.01 | 103.8 |
| 10 | 84.7 | 13.7 | 13.5 | 2.03 | 0.06 | 0.03 | 114.0 |
| 15 | 71.5 | 18.3 | 11.5 | 1.31 | 0.09 | 0.04 | 102.7 |
| 20 | 66.5 | 21.9 | 11.6 | 1.27 | 0.11 | 0.05 | 101.4 |
| 25 | 61.4 | 22.9 | 9.9 | 0.96 | 0.11 | 0.05 | 95.4 |
| 30 | 67.5 | 24.9 | 8.3 | 1.06 | 0.09 | 0.06 | 102.0 |

Cross-checks: the CML and CEL bands fall inside the caption's printed ranges except M1 130 C CML, which
reaches 0.08 % (caption 0.06 %), and M2 120 C CEL at 30 min, 0.01 % (caption 0.02 % floor); both are within
two vertex quanta. Recoveries match the main text (M2 120 C ≈ 110 % at 15-20 min; M2 130 C 101-114 %; M1
≈ 85 % and M2 ≈ 95 % at 30 min, 120 C) except **M1 120 C at 10 min: 107.7 % in the figure, "reached 104%"
in the text**. The main text's "relative" numbers for M1 120 C ("0.25, 0.45, 0.2 and 4.5") are reproduced
by this figure as fractional glucose loss (0.26), fractional lysine loss (0.44), fructose / initial glucose
(0.15) and AP(30) / AP(0) (4.4): the "4.5" is a ratio to the t = 0 AP background, not to initial lysine.

Inference (mine, not printed): on this molar basis, lysine / sugar at t = 0 is 0.154 (M1) and 0.176-0.177
(M2), i.e. sugar : lysine 6.5 : 1 and 5.7 : 1, not the design's "about 10:1". With the weighed-in glucose
150 mmol/L that would put lysine residues at ~23 mmol/L, above the 15-17 mmol/L inferred in section 2 and
above what 30 g/L caseinate can carry (~16.8 mmol/L); with the measured t = 0 sugar below the weighed value
(the paper says ~21 % below for lactose: 112 mmol/L -> lysine ~20 mmol/L) the gap narrows. The t = 0 sugar
and lysine concentrations remain unprinted; do not back out an absolute lysine from this ratio alone.

### Fig. S2A / S2B — the ODE systems (transcribed exactly from the typed equations)

Caption: "Reaction rates derived from the kinetic model for N-carboxymethyllysine (CML) formation in M1
system" (S2B: M2). Abbreviations as printed: Glu glucose, fru fructose, lacto lactose, lactu lactulose, lys
lysine, AP Amadori product, GO glyoxal, MRPs advanced Maillard reaction products, CEL. Concentrations in
mmol/L and time in min (main paper 2.6.2 and Table 1 units; the SI itself prints no units).

M1 (glucose-casein):

    d[glu]/dt  = -k1[glu] - k3[glu][lys] + k4[fru] - k5[glu]
    d[fru]/dt  =  k1[glu] - k2[fru] - k4[fru]
    d[lys]/dt  = -k3[glu][lys] - k6[lys][GO]
    d[AP]/dt   =  k3[glu][lys] - k7[AP] - k8[AP] - k9[AP]
    d[CEL]/dt  =  k9[AP] - k10[CEL]
    d[CML]/dt  =  k7[AP] - k11[CML] + k6[lys][GO]
    d[MRPs]/dt =  k8[AP]
    d[X1]/dt   =  k2[fru]
    d[X3]/dt   =  k11[CML]
    d[X4]/dt   =  k10[CEL]
    d[GO]/dt   =  k5[glu] - k6[lys][GO]

M2 (lactose-casein): identical with glu -> lacto, fru -> lactu and X1 -> X2:

    d[lacto]/dt = -k1[lacto] - k3[lacto][lys] + k4[lactu] - k5[lacto]
    d[lactu]/dt =  k1[lacto] - k2[lactu] - k4[lactu]
    d[lys]/dt   = -k3[lacto][lys] - k6[lys][GO]
    d[AP]/dt    =  k3[lacto][lys] - k7[AP] - k8[AP] - k9[AP]
    d[CEL]/dt   =  k9[AP] - k10[CEL]
    d[CML]/dt   =  k7[AP] - k11[CML] + k6[lys][GO]
    d[MRPs]/dt  =  k8[AP]
    d[X2]/dt    =  k2[lactu]
    d[X3]/dt    =  k11[CML]
    d[X4]/dt    =  k10[CEL]
    d[GO]/dt    =  k5[lacto] - k6[lys][GO]

Reaction scheme (the Fig. 4 candidate network these equations encode): sugar ⇌ isomer (k1 forward, k4 back);
isomer -> X1/X2 (k2); sugar + lys -> AP (k3, second order); sugar -> GO (k5); GO + lys -> CML (k6, second
order); AP -> CML (k7); AP -> MRPs (k8); AP -> CEL (k9); CEL -> X4 (k10); CML -> X3 (k11). Every other step
is first order in its one source. Table 1 then sets k4 = 0 (M1, and M2 at 120 C) and k5 = k6 = 0 (the GO
route, dropped from Fig. 5), so the fitted model is the system above without those terms.

What this settles (supersedes flag 6): **no step returns lysine.** Lysine is consumed only by k3 (and the
zeroed k6); AP -> CML, AP -> CEL and AP -> MRPs leave d[lys]/dt untouched, so the AP's lysine is counted as
lost into CML, CEL or MRPs. The main text's remark that the CEL route involves "release of lysine" is not in
the equations. GO has no sink other than k6. The isomer (fructose / lactulose) does not react with lysine.
Note for reuse: a network that writes AP -> 3-DG + lysine back (an Amadori decay that regenerates the amine)
and gives it Nguyen's k8 is NOT Nguyen's model; the constant was fitted with the lysine kept out.

### Table S1 — paired comparisons (typed table, transcribed)

Title: "Comparision of the concentrations of Amadori product (AP), CML, and CEL under certain tested conditions
using t-test with a significance value set to 0.05." Footnotes: "SME: Standard error mean"; "df: degrees of
freedom". **Unit of the paired difference not printed** (magnitudes consistent with mmol/L, the paper's
modelling unit; inference). df = 13 implies 14 pairs, i.e. 7 times × 2 replicates (inference).

| pair | compound | system | temp. (C) | paired difference, mean | SME | t | df | p |
|---:|---|---|---|---:|---:|---:|---:|---:|
| 1 | AP | M1 - M2 | 120 | 0.132 | 0.135 | 0.978 | 13 | 0.346 |
| 2 | AP | M1 - M2 | 130 | -0.065 | 0.063 | -1.038 | 13 | 0.318 |
| 3 | CML | M1 | 120 - 130 | -0.003×10^-1 (as printed, = -0.0003) | 0.003 | -0.126 | 13 | 0.902 |
| 4 | CML | M1 - M2 | 120 | -0.008 | 0.002 | -3.708 | 13 | 0.003 |
| 5 | CML | M1 - M2 | 130 | -0.015 | 0.004 | -3.398 | 13 | 0.005 |
| 6 | CEL | M1 - M2 | 120 | 0.004 | 0.001 | 2.633 | 13 | 0.021 |
| 7 | CEL | M1 - M2 | 130 | 0.009 | 0.004 | 2.260 | 13 | 0.042 |

Check: mean / SME reproduces t to the rounding of the printed SME (pair 1 0.98, pair 2 -1.03, pair 3 -0.10,
pair 7 2.25; pairs 4-6 drift because SME 0.001-0.002 is printed to one significant figure). Reading: AP does
not differ between the sugars; CML in M1 does not differ between 120 and 130 C; M2 has more CML and less CEL
than M1 at both temperatures.

### Supplement flags

S1. pH is a cooled-sample room-temperature reading; to compare with a model's hot-solution pH, correct for
    the phosphate pKa2 temperature coefficient (about −0.0028 /K; the hot pH is lower). Neither the SI nor the
    main text gives the hot value.
S2. The pH points are per replicate with no error bars; the 30 min duplicates differ by 0.06 pH (M1
    130 C) and 0.05 (M2 120 C), larger than the early-time spread. No acid concentrations were measured in this paper,
    so the drift is a forcing function / check for an acid-production sub-model, not a calibration of it.
S3. The mass balance mixes external-calibration sugars, SIDA lysine / CML / CEL and furosine-derived AP; its
    excursions to 107-114 % are measurement scatter by the authors' own reading. The bands are single values
    per time (no SD).
S4. Text-figure conflict: M1 120 C recovery at 10 min is 107.7 % in Fig. S1B, 104 % in the text.
S5. Authorship: the main paper's byline is Nguyen, van der Fels-Klerx & van Boekel (three authors); the SI
    file carries no byline. A four-author citation adding Peters does not match this paper's byline.
