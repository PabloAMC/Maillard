# Göncüoğlu Taş & Gökmen 2016 — EXTRACTION (hazelnut roasting, 26-step multiresponse model, 150/160/170 C)

**Source on disk:** `data/articles/Goncouglu2016.pdf` (54 pp., Elsevier accepted-manuscript PDF with line
numbers; text layer clean for body and Table 1; Appendix A differential equations are garbled in the text
layer — rate-constant subscripts are lost — and were NOT used; the scheme below is taken from Table 1's
step list and Fig. 2b). Supplementary Tables S1-S4 and Figures S1-S10 are **not in the PDF**; in
particular Table S4 (Ea and k_b at 160 C) is absent, so no activation energy from the authors is on disk.
Read-only extraction, 2026-09-07. Sibling: `goncuoglutas2017_extraction.md` (the 2017 paper by the same
authors); `kocadagli2016foodchem_extraction.md` (same lab, glucose/flour, same dicarbonyl set).

## 0. Identity

| field | value |
|---|---|
| Title | "Maillard Reaction and Caramelization during Hazelnut Roasting: A multiresponse kinetic study" |
| Authors | Neslihan Göncüoğlu Taş, Vural Gökmen (Hacettepe University, Ankara) |
| Venue | Food Chemistry, accepted manuscript (FOCH 20292); received 30 Aug 2016, revised 27 Nov 2016, accepted 30 Nov 2016 |
| DOI as printed | `http://dx.doi.org/10.1016/j.foodchem.2016.11.159` |

## 1. Why it matters to the model

Same laboratory and same analytical/modelling machinery as Kocadağlı & Gökmen 2016 (the source of the
repo's glucose-glass constants for glucosone, glyoxal, diacetyl), applied to a real low-moisture food at
150-170 C with a 26-step network that contains the trunk's steps by name: glucose -> 1,2-enediol ->
fructose; glucose + AA -> Amadori; Amadori -> 3-DG / 1-DG; 3-DG -> 3,4-DG -> HMF; 1-DG -> methylglyoxal /
dimethylglyoxal (diacetyl); glucose -> glyoxal; and dicarbonyl / HMF sinks. It is a second set of
fitted constants for those steps, in a **dry nut matrix (aw 0.40-0.55, moisture 2.5-5 %)**, at
temperatures 30-50 C above the trunk's range. Its main structural findings: HMF forms mainly from the
sucrose-derived fructofuranosyl cation, not via 3-DG; MGO and DMG come from 1-DG with the fastest constants
in the network; glyoxal comes from glucose directly; Strecker-type dicarbonyl + AA steps had to be removed.

## 2. Methods as they matter to a model

- **Material:** Tombul hazelnuts, 5 g unshelled per run, oven (Memmert UNE 400) at 150, 160, 170 C for
  15, 30, 60, 90, 120 min; quenched at -18 C; skins removed; ground. Replicates not stated for roasting;
  concentrations carry ± (SD, n unstated).
- **Initial composition (dry weight):** sucrose 5.5 ± 0.1 g/100 g (= 161 mmol/kg); fructose 0.4 ± 0.02
  g/100 g (22 mmol/kg); glucose 0.2 ± 0.06 g/100 g (11 mmol/kg); total free amino acids 2112 ± 49 mg/kg;
  protein-bound lysine 5401 ± 50 mg/kg; "total amino acids" used in the model 7513 ± 87 mg/kg
  (~51 mmol/kg, Fig. 1 axis); raw glyoxal 1.7 ± 0.6 mg/kg; glucosone absent at every point.
- **Moisture/pH:** decreased during roasting (Table S1, not on disk). Heat-transfer lag: authors cite
  4-6 min to reach oven temperature (Demir 2002) and choose NOT to model it.
- **Analytics:** water extraction (0.5 g in 10 mL, 3 steps); sugars by HPLC-RI (Shodex KC-811,
  0.1 % H3PO4); HMF by HPLC (Kocadağlı 2012); alpha-dicarbonyls as OPD quinoxalines by HPLC-ESI-MS SIM
  (m/z 251 glucosone, 235 1-/3-DG, 217 3,4-DG, 159 DMG, 145 MGO, 131 GO), quantified against quinoxaline,
  2-methylquinoxaline, 2,3-dimethylquinoxaline and 3-DG standards (1-DG and 3,4-DG therefore on the 3-DG
  response — semi-quantitative); free amino acids and acid-hydrolysed lysine by UPLC-ESI-MS/MS.
- **Quantification basis:** umol/kg hazelnut (dry weight) in Fig. 1; mg/kg dw in the text.
- **Modelling:** Athena Visual Studio 14.2, determinant criterion; 95 % HPD intervals; reparametrised
  Arrhenius with T_b = 160 C fitted to all three temperatures at once (Fig. S2, Table S4 — not on disk).

## 3. Tables and scheme

**Table 1** — "Reaction rate constants with 95% highest posterior density (HPD) intervals at different
temperatures according to the proposed kinetic model in Figure 2b for Maillard reaction and caramelization
during roasting of hazelnuts." Units as printed: first-order steps in min^-1 x 10^3; bimolecular steps
(5, 9, 10) in kg x umol^-1 x min^-1 x 10^3. "ind": "indeterminate, which means a large uncertainty in the
estimated parameter within 95% confidence interval."

| # | step | unit | k 150 C | HPD | k 160 C | HPD | k 170 C | HPD |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | SUC -> GLC + FFC | min^-1 x10^3 | 6.9 | ±0.8 | 15 | ±3.1 | 22 | ±2.0 |
| 2 | GLC -> 1,2-ED | min^-1 x10^3 | 141 | ±27.0 | 473 | ±131 | 698 | ±178 |
| 3 | 1,2-ED -> GLC | min^-1 x10^3 | 0 | ±0 | 8.5 | ±3.6 | 28 | ±8.3 |
| 4 | GLC -> 3-DG | min^-1 x10^3 | 0.03 | ±0.02 | 0 | ±0 | 0 | ±0 |
| 5 | GLC + AA -> AP | kg umol^-1 min^-1 x10^3 | 0.0009 | ±0.0007 | 0.003 | ±0.001 | 0.009 | ±0.001 |
| 6 | GLC -> GO | min^-1 x10^3 | 0.6 | ±0.2 | 2.5 | ±0.9 | 9.2 | ±1.0 |
| 7 | 1,2-ED -> FRU | min^-1 x10^3 | 1.3 | ±0.9 | 1.8 | ±0.7 | 4.2 | ±1.9 |
| 8 | FRU -> 1,2-ED | min^-1 x10^3 | 0 | ±0 | 0 | ±0 | 41 | ±14 |
| 9 | FRU + AA -> HP | kg umol^-1 min^-1 x10^3 | 0.00023 | ±0.00007 | 0.00062 | ±0.00015 | 0 | ±0 |
| 10 | FFC + AA -> HP | kg umol^-1 min^-1 x10^3 | 0.00094 | ±0.00028 | 0.00027 | ±0.00030 | 0.00004 | ±0.00002 |
| 11 | FFC -> HMF | min^-1 x10^3 | 0.58 | ±0.13 | 0.57 | ±0.19 | 2.02 | ±1.32 |
| 12 | HP -> 1-DG | min^-1 x10^3 | 0.23 | ±0.18 | 0.85 | ±0.58 | 267 | ±177 |
| 13 | HP -> 3-DG | min^-1 x10^3 | 0.009 | ±0.004 | 0.022 | ±0.030 | 12 | ind |
| 14 | AP -> 3-DG | min^-1 x10^3 | 0 | ±0 | 0.62 | ±0.03 | 0 | ±0 |
| 15 | AP -> 1-DG | min^-1 x10^3 | 3.47 | ±3.12 | 3.51 | ±1.33 | 0.56 | ±0.58 |
| 16 | 1-DG -> MGO | min^-1 x10^3 | 7012 | ±5510 | 33016 | ind | 47920 | ±33240 |
| 17 | 1-DG -> DMG | min^-1 x10^3 | 371 | ±241 | 895 | ±581 | 1073 | ±618 |
| 18 | 3-DG -> 3,4-DG | min^-1 x10^3 | 4.27 | ±0.57 | 29.4 | ±26.1 | 88.1 | ±22.7 |
| 19 | 3,4-DG -> HMF | min^-1 x10^3 | 0 | ±0 | 134 | ±127 | 390 | ±111 |
| 20 | HP -> P1 | min^-1 x10^3 | 11 | ±1.4 | 4.7 | ±5.6 | 59 | ind |
| 21 | AP -> P2 | min^-1 x10^3 | 140 | ±122 | 21.2 | ±34.4 | 5.24 | ±1.75 |
| 22 | 1-DG -> P3 | min^-1 x10^3 | 122 | ind | 0 | ind | 404 | ind |
| 23 | MGO -> P4 | min^-1 x10^3 | 126 | ±106 | 737 | ±113 | 918 | ±639 |
| 24 | DMG -> P5 | min^-1 x10^3 | 54 | ±41 | 130 | ±90 | 106 | ±65 |
| 25 | GO -> P6 | min^-1 x10^3 | 18 | ±8.3 | 61 | ±25 | 290 | ind |
| 26 | HMF -> P7 | min^-1 x10^3 | 12 | ±3.7 | 21 | ±11 | 103 | ±63.7 |

Abbreviations (footnote a): SUC sucrose; GLC glucose; FRU fructose; FFC fructofuranosyl cation; 1,2-ED
1,2-enediol; AP Amadori product; HP Heyns product; 1-DG 1-deoxyglucosone; 3-DG 3-deoxyglucosone; 3,4-DG
3,4-dideoxyglucosone; GO glyoxal; MGO methylglyoxal; DMG dimethylglyoxal; HMF 5-hydroxymethylfurfural;
AA total amino acids; P products.

**Scheme fitted (Fig. 2b):** the 26 steps above, exactly. Water is released at steps 5, 9, 10 (condensation),
11, 18, 19 (dehydration); AA is released at 12, 13, 14, 15 (AP/HP degradation). Steps present in the
comprehensive Fig. 2a but removed for Fig. 2b: 27 FRU -> FFC (fitted to zero), 28 3-DG -> MGO and 29
GLC -> MGO (removed after one-at-a-time exclusion; best fit with 1-DG -> MGO only), 30 3,4-DG -> P (and 3-DG
degradation; "their fits were better in that case"), and all dicarbonyl + AA (Strecker) steps (amino-acid
fit "not well compatible", Fig. 5). Alternative models 1-4 are in Figs. S3-S10 (not on disk).

**Activation energies:** text only — "The activation energies of elementary reaction steps were found to
range between 0-1174 kJ/mol with six zero and a few relatively high values" (Table S4, not on disk);
"The temperature dependence of the reactions was found to be more complicated than defined by the Arrhenius
equation." No per-step Ea is printed in the manuscript.

**Concentrations printed in the text (mg/kg dw unless stated):**

| quantity | 150 C | 160 C | 170 C |
|---|---:|---:|---:|
| sucrose loss at 120 min | 60 % | 70 % | 90 % |
| fructose / glucose after first significant drop (g/100 g) | 0.26 ± 0.02 / 0.08 ± 0.01 (30 min) | 0.21 ± 0.01 / 0.07 ± 0.01 (15 min) | 0.18 ± 0.01 / 0.06 ± 0.02 (15 min) |
| total amino acid loss at 120 min | 68 % | 81 % | 85 % |
| 3-DG at 120 min | 6.7 ± 0.1 | 6.1 ± 0.1 | max 5.4 ± 0.1 at 60 min, then falls |
| 3,4-DG | "approximately 5 times lower than 3-DG"; short lag | | |
| 1-DG maximum | 0.22 ± 0.01 | 0.31 ± 0.03 | 0.27 ± 0.01 |
| glyoxal | raw 1.7 ± 0.6; "increased up to 4 times after 15 min of roasting and did not change during prolonged roasting at all roasting temperatures" | | |
| methylglyoxal maximum | | 6.6 ± 0.5 (90 min) | |
| dimethylglyoxal | "no significant differences (p < 0.05)" between 150 and 160 C at 120 min | | |
| HMF at 120 min | 104 ± 0.5 | 238 ± 1.9 | 278 ± 0.7 |
| molar recovery of measured species | 90 % at 15 min, 39 % at 120 min | 71 % / 29 % | 62 % / 15 % |

Fig. 1 axis ranges (umol/kg): sucrose 0-200 000; fructose 0-25 000; glucose 0-12 000; 3-DG 0-50; 3,4-DG
0-12; HMF 0-2500; 1-DG 0-2.5; MGO 0-120; DMG 0-20; GO 0-150; total AA 0-60 000. Converting the text maxima:
3-DG 41 umol/kg, MGO 92, GO ~117 (4 x 1.7 mg/kg), 1-DG 1.9, DMG ~15 (axis), 3,4-DG ~8 — molar ordering in
roasted hazelnut **GO ~ MGO > 3-DG > DMG > 3,4-DG > 1-DG; glucosone absent**.

## 4. What the repo could take

Rate constants (FIT-eligible under the owner's rule, but see caveats — dry matrix, 3-DG-basis quantitation
of 1-DG/3,4-DG, HPD often > 50 % of k):

| candidate | numbers | note |
|---|---|---|
| k(GLC + AA -> AP), step 5 | 0.9 / 3 / 9 kg mol^-1 min^-1 at 150/160/170 C (from 0.0009e-3 kg umol^-1 min^-1 x 10^6) | trunk `k_schiff` (Martins X = 1.6e-5 L mmol^-1 min^-1 at 100 C, Ea 96.8) extrapolates to 0.64 / 1.2 / 2.2 L mol^-1 min^-1 at 150/160/170 C [D]: same order, hazelnut 1.4-4x higher. Implied Ea from the three k [D]: 180 kJ/mol. |
| k(GLC -> GO), step 6 | 0.6 / 2.5 / 9.2 e-3 /min | implied Ea [D] 213 kJ/mol; second lab for glucose -> glyoxal after Kocadağlı 2016 |
| k(1-DG -> DMG), step 17 | 0.371 / 0.895 / 1.073 /min | implied Ea [D] 83 kJ/mol; the only diacetyl-formation constant series outside Kocadağlı 2016 |
| k(1-DG -> MGO), step 16 | 7.0 / 33 (ind) / 48 /min | text: "7, 33 and 48 min^-1"; ~20x Kocadağlı 2016 flour (1.6 /min at 160 C); authors: lipid MGO "not expected to be more than 10 %" |
| k(3-DG -> 3,4-DG), step 18 | 4.27 / 29.4 / 88.1 e-3 /min | implied Ea [D] 236 kJ/mol; the HMF-route rate-determining step |
| k(3,4-DG -> HMF), step 19 | 0 / 134 / 390 e-3 /min | text: "almost 5 times higher than" step 18 (at 160/170 C); zero at 150 C — RATIO-ONLY |
| k(SUC -> GLC + FFC), step 1 | 6.9 / 15 / 22 e-3 /min | implied Ea [D] 90 kJ/mol; cleanest series (HPD <= 20 %) |
| k(MGO -> P), k(DMG -> P), k(GO -> P), k(HMF -> P), steps 23-26 | see table | sink order MGO > DMG > GO > HMF at all T (text); HMF sink 12 / 21 / 103 e-3 /min (Ea [D] 168) |
| k(GLC -> 1,2-ED), step 2 | 141 / 473 / 698 e-3 /min | Ea [D] 125 kJ/mol; but 1,2-ED -> FRU (step 7) is 1.3-4.2 e-3 only, so the enediol is a modelling reservoir, not a measured species |

Within-study ratios: k5 / k9 (Amadori vs Heyns formation) = 3.9 (150 C), 4.8 (160 C) — text "almost 5
times"; k16 / k17 (MGO vs DMG from 1-DG) = 19 / 37 / 45 — text "20, 40 and 45-fold"; k19 / k18 = 4.6 (160),
4.4 (170); HMF(120 min) 160/150 = 2.3, 170/150 = 2.7; 3-DG / 3,4-DG ~ 5; sucrose k 170/150 = 3.2.

Directional claims: (i) HMF from FFC dominates over the 3-DG route (Fig. 4: excluding step 11 leaves
predicted HMF "far below" measured) — in a sucrose-rich dry matrix; (ii) MGO and DMG originate from 1-DG,
not 3-DG or glucose (model discrimination); (iii) glyoxal from glucose directly, 4x the raw level by 15 min
then flat; (iv) glucosone not detected at 150-170 C in hazelnut (contrast Leitzen 2021 aqueous glucose
at 121 C, where glucosone > GO > MGO — matrix/temperature flips the ordering); (v) dicarbonyl + amino
acid steps are kinetically unnecessary here (Fig. 5); (vi) the 3-DG maximum moves earlier with T (120 min
at 150-160 C, 60 min at 170 C).

## 5. Caveats

- **Matrix:** dry roasted nut, aw 0.40-0.55, 56 % oil; concentrations per kg hazelnut; bimolecular units
  kg umol^-1 min^-1 — not transferable to mol/L water without an assumed reactive-phase volume.
- **Non-isothermal start:** 4-6 min to reach oven temperature, not modelled; sample is 5 g whole nuts.
- **Arrhenius fails by the authors' own statement** (Ea 0-1174 kJ/mol, six zeros); several constants are
  non-monotonic (steps 4, 9, 10, 14, 15, 20, 21, 24) or zero at one temperature (3, 8, 14, 19); five (13,
  16, 20, 22, 25) are "ind" at one or more temperatures. Only steps 1, 2, 6, 17, 18, 23 (partly), 26 are
  monotonic with HPD < 100 %.
- 1-DG, 3,4-DG quantified on the 3-DG calibration; glucosone/1-DG/3-DG share m/z 235/251 — same
  semi-quantitation flag as Kocadağlı 2016.
- Molar recovery of measured species falls to 15-39 % at 120 min: most of the carbon ends in unmeasured
  products, so the "P" sinks absorb everything.
- Table S4 (authors' Ea, k_b at 160 C), S1 (pH, moisture), S2 (individual amino acids), S3 (colour) and
  the alternative models are not in the PDF. Appendix A equations unreadable in the text layer.
- Lipid-derived GO/MGO cannot be excluded (authors estimate < 10 %).
