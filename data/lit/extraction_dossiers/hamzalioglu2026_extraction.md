# Hamzalıoğlu 2026 — EXTRACTION (pasteurised whole milk heated 110/120/130/140 C for 0.5-5 min; 24-step multiresponse model on lactose, total lysine, lactulosyllysine, 3-DG, 1-DG, glucosone, diacetyl, methylglyoxal, glyoxal, CML, CEL; rate constants at four temperatures with 95 % HPD and reparameterised Arrhenius barriers)
### A matrix comparator, not a trunk source: lactose + protein-bound lysine in milk, with every dicarbonyl the repository lacks printed as a fitted rate constant at 110-140 C — and a fitted glyoxal sink of zero at every temperature.

**Source on disk:** `data/articles/Hamzalioglu2026.pdf` (14 pp., open access CC-BY 4.0, owner's
download, 2026-09-08). Read from the text layer (`scratchpad/articles/Hamzalioglu2026.txt`); the
two kinetic tables came through with each mantissa, exponent and HPD on separate lines and were
re-assembled row by row below (checked against the in-text quotations of k1, k9, k10, k7, k5,
k24). Page 7 (journal p. 15782) was rasterised at 75 dpi to read the axis labels of Figure 3 only
(all y axes in "umol/kg dm", x "time (min)"). No value was read off any figure. The Supporting
Information (Figure S1 mass balances; Tables S1-S2 parameter correlation matrices) is NOT on disk.
The task brief gives PMC13220311; the PDF carries only the DOI.

## 0. Identity

| field | value |
|---|---|
| Title | "Multi-Response Kinetic Modeling of the Maillard Reaction in Milk During Heating at Ultra-High-Temperature Range" |
| Authors | Aytül Hamzalıoğlu, Işıl Aktağ, Vural Gökmen* (Hacettepe University, Ankara; Munzur University) |
| Venue | J. Agric. Food Chem. 2026, 74, 15776-15789. Received October 24 2025, revised May 1 2026, accepted May 5 2026, published May 14 2026 |
| DOI | 10.1021/acs.jafc.5c14296 |
| Naming | Lac = lactose; Lys = **total (bound + free) lysine**; LacLys = lactulosyllysine (the Amadori product, from furosine x 3.1); 3-DG, 1-DG = deoxyglucosones; G = glucosone; DA = diacetyl; MGO = methylglyoxal; GO = glyoxal; CML / CEL = N-epsilon-carboxymethyl / carboxyethyl-lysine; Int = an unquantified lactose isomer pool; P1-P10 = unquantified products |
| Companions | Kocadagli & Gokmen 2016 (same group, same Athena workflow; `kocadagli2016jafc_extraction.md` and `kocadagli2016foodchem_extraction.md`, the source of the trunk's B13 dicarbonyl constants); Brands & van Boekel 2001 / 2003 (disaccharide-casein, the furosine x 3.1 factor); Aktag 2019 (the group's UHT milk survey, ref 44); Kocadagli & Gokmen 2014 (the alpha-DC method, ref 33); Liu 2025 (glucose-lysine multiresponse, ref 54, not on disk) |

## 1. Why it matters

The trunk's glyoxal comes from a single dry-glass entry (glucose -> glucosone, Kocadagli 2016 at
160-200 C, `src/kinetic_core/parameters_dicarbonyl.py` `k_glc_g`), glucosone -> glyoxal is the same
paper's `k_g_go` (Ea 93.8), and the glyoxal sink `k_go_sink` is a 180 C rate with its barrier fixed
to zero. This paper is the same laboratory's multiresponse treatment of an **aqueous** food at
**110-140 C**, and it prints, at four temperatures with intervals:

- Amadori -> glucosone (`k4`), glucosone -> glyoxal (`k10`), Amadori -> glyoxal directly (`k7`);
- Amadori -> 1-DG (`k3`) and 1-DG -> methylglyoxal (`k9`), the route the fit prefers for
  methylglyoxal over any 3-DG route;
- sinks for 3-DG, 1-DG, glucosone, methylglyoxal, glyoxal (`k17`-`k21`): the **glyoxal sink is
  zero at all four temperatures** and the 3-DG sink is zero at 120-140 C;
- CML from glyoxal + lysine versus from Amadori oxidation, CEL from methylglyoxal + lysine versus
  from Amadori: the Amadori route wins for both.

What transfers to the trunk and what does not is set out in section 4 (bottom). The short version:
the barriers and the direction statements transfer as comparators; the constants do not, because
the sugar is lactose, the amine is casein lysine, the concentration basis is per kg of milk dry
matter, and the sampling (0.5-5 min) cannot resolve first-order constants above a few per minute.

## 2. Methods as they matter to a model

- **Matrix.** Pasteurised whole milk from a local market (already pasteurised, to mimic industrial
  practice). Lactose 41.77 +/- 1.25 mg/g milk = 939 150 umol/kg dry matter; total lysine 10.4 +/-
  0.3 g/100 g protein = 194 670 umol/kg dm; milk pH ~ 6.7, dropped by "only 0.1 after heating for
  5 min". **Dry-matter fraction implied by the two lactose figures: 41.77 / (0.93915 x 342.30) =
  0.130** (i.e. 13.0 % solids, the usual value for whole milk). **All concentrations in the model
  and the figures are umol/kg dm** (Figure 3 axes; the Methods sentence "umol/kg of milk" is loose).
- **Heating.** 1 mL of milk into pre-heated glass tubes in a heating metal block; temperature rise
  monitored with embedded thermocouples, lag "< 5 s". **110 C: 1, 2, 3, 4, 5 min; 120 C: 1-5 min;
  130 C: 1, 1.5, 2, 2.5, 3 min; 140 C: 0.5, 1, 1.5, 2, 2.5 min.** Ice quench, -80 C, freeze-dried,
  250 mg powder triple-extracted with water (1.25 + 0.625 + 0.625 mL). Triplicate heat treatments;
  individual replicate values used in the regression.
- **Lactose**: HPLC-RI, Shodex RSpak KC-811, 0.1 % H2SO4, external calibration 0.005-1 %. Glucose
  and galactose monitored, "a significant change ... could not be detected".
- **Lysine**: free lysine from the aqueous extract (acetonitrile-clarified); **total lysine from 8 N
  HCl hydrolysis at 110 C for 23 h**; LC-MS/MS (ZIC-HILIC), theanine internal standard, MRM,
  1-100 uM calibration. Total lysine loss: 15.9 % (110 C, 5 min), 27 % (140 C, 2.5 min).
- **Lactulosyllysine**: furosine by HILIC-DAD 280 nm (1-10 mg/L standards) on the acid
  hydrolysate, **LacLys = furosine x 3.1** (Brands & van Boekel). Furosine 1.04 g/kg protein
  (pasteurised) -> 3.08 g/kg protein (120 C, 5 min).
- **alpha-Dicarbonyls** (the method that matters here): aqueous extract + acetonitrile 1:1,
  15 000 g; **500 uL supernatant + 150 uL of 0.2 % OPD with 11 mM DTPA + 150 uL 0.5 M sodium
  phosphate pH 7; room temperature, dark, 2 h**; LC-MS/MS (Agilent 1260 II + QqQ, ESI+), Zorbax
  Eclipse XDB-C18 4.6 x 150 mm 5 um, 1 % formic acid in water / methanol 20 -> 60 % B in 8 min,
  1 mL/min, 40 C. SIM [M+H]+: glucosone-quinoxaline 251; 3-DG and 1-DG 235; DA (2,3-dimethyl-
  quinoxaline) 159; MGO (2-methylquinoxaline) 145; GO (quinoxaline) 131. Retention: G 3.0, 1-DG
  3.5, 3-DG 4.0, GO 5.6, MGO 6.6, DA 7.5, 5-methylquinoxaline (IS, 0.5 mg/L) 8.1 min. Calibration
  0.1-5 mg/L on derivatised G and 3-DG working solutions and on authentic quinoxaline and
  2-methylquinoxaline; **1-DG semi-quantified on the 3-DG curve**; recoveries 102.4 % (3-DG),
  103.6 % (GO), 99.5 % (MGO); LOD/LOQ at S/N 3 / 10 (values not printed here; "previously
  validated in infant formula", ref 33). Note the 2 h OPD incubation at neutral pH versus the 12-24 h
  incubations of the Jiangnan papers (`xia2022_extraction.md`, `yu2020_extraction.md`).
- **CML / CEL**: NaBH4 reduction (borate pH 9.2, 4 h) then 8 N HCl 110 C 24 h; LC-MS/MS (Atlantis
  T3), MRM 205.1 -> 84/130 (CML), 219.0 -> 84/130 (CEL); **matrix-matched calibration** in milk
  heated 110 C 1 min, 0-20 ug/mL.
- **Model.** Athena Visual Studio 14.2, determinant criterion, **each temperature fitted
  separately**, 95 % HPD intervals; model discrimination from a comprehensive network (Figure 1) to
  the 24-step network of Figure 2. Reparameterised Arrhenius k = k_ref exp(-(Ea/R)(1/T - 1/T_ref)),
  **T_ref = 120 C = 393.15 K**. Mass balance recovery (Figure S1): 89.8 / 86.1 / 82.5 / 80.6 % at
  110 / 120 / 130 / 140 C.
- **Concentration levels printed in the text** (the only ones; all else figure-only):

  | species | pasteurised milk | after heating |
  |---|---|---|
  | lactose | 41.77 mg/g (939 150 umol/kg dm) | 37.92 mg/g (836 600 umol/kg dm) at 110 C 5 min; -9.21 % (110 C 5 min), -16.10 % (140 C "3 min" as printed; the 140 C series ends at 2.5 min) |
  | total lysine | 194 670 umol/kg dm | -15.9 % (110 C series), -27 % (140 C 2.5 min) |
  | furosine | 1.04 g/kg protein | 3.08 g/kg protein (120 C 5 min) |
  | 3-DG | 0.31 mg/L (17.19 umol/kg dm) | 0.91 mg/L (51.17 umol/kg dm) at 110 C 5 min |
  | 1-DG | 0.04 mg/L (2.43 umol/kg dm) | 0.21 mg/L (11.79 umol/kg dm) at 110 C 5 min |
  | MGO | — | 0.08 mg/L (10.12 umol/kg dm) at 110 C 5 min |
  | GO | — | 0.1 mg/L (14.07 umol/kg dm) at 120 C 5 min |
  | CML | 0.11 mg/kg protein (0.14 umol/kg dm) | 1.51 mg/kg protein at 110 C 5 min |
  | CEL | 0.24 mg/kg protein (0.24 umol/kg dm) | 0.72 mg/kg protein at 110 C 5 min |

  Basis check: 0.31 mg/L / 162.14 = 1.91 umol/L of milk; 17.19 umol/kg dm x 0.130 = 2.23 umol/kg
  milk — consistent to 15 % (the mg/L figures are per litre of milk, the umol figures per kg dm).
  So in the aqueous phase of milk the dicarbonyls sit at **1-15 umol/L** after UHT-range heating,
  three orders of magnitude below the mmol/L axis of Xia 2022's glucose-glycine pot.

- **Unit conversions for a comparison with the trunk (mine).** First-order constants (min-1) need
  none. Second-order constants are in kg dm umol-1 min-1. To re-express per litre of milk water:
  1 umol/kg dm = 0.130 umol/kg milk = 0.130 / 0.87 = 0.149 umol/L water (density taken as 1), so
  k' [L umol-1 min-1] = k / 0.149 and k'' [L mmol-1 min-1] = k x 6.7e3. This is an assumption about
  where the reaction happens (in the serum, with casein lysine counted as dissolved). On a
  per-kg-of-whole-milk basis instead, k [kg dm umol-1 min-1] x (1/0.130) x 1000 = k x 7.7e3 in
  kg milk mmol-1 min-1. The two bases differ by 15 %; the dm-to-liquid step (x 6.7e3 to 7.7e3) is
  the large factor.

## 3. Tables re-typed

### Table 1. "Estimated Reaction Rate Constants (k) with 95 % Highest Posterior Density (HPD) Intervals at Different Temperatures According to the Proposed Kinetic Model in Figure 2"

"ind" = indeterminate ("a large uncertainty exists in the estimated parameter within the 95 %
confidence interval"). "0.0 +/- 0.0" is printed as such. Second-order steps (1, 11, 12) in
kg umol-1 min-1; all others min-1.

| # | step | unit | 110 C k | HPD | 120 C k | HPD | 130 C k | HPD | 140 C k | HPD |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Lac + Lys -> LacLys | kg umol-1 min-1 | 2.2e-9 | 2.6e-10 | 7.5e-9 | 1.3e-9 | 8.1e-9 | 1.8e-9 | 2.5e-8 | 4.5e-9 |
| 2 | LacLys -> 3-DG | min-1 | 4.4e-3 | 1.9e-3 | 5.9e-3 | 5.7e-4 | 6.8e-3 | 8.4e-4 | 7.7e-3 | 8.7e-4 |
| 3 | LacLys -> 1-DG | min-1 | 4.1e-2 | 2.6e-3 | 6.9e-2 | ind | 2.9e-3 | 2.3e-4 | 1.9e-1 | 9.2e-3 |
| 4 | LacLys -> G | min-1 | 1.9e-2 | 4.0e-3 | 2.1e-2 | 1.0e-3 | 2.5e-2 | 7.1e-2 | 1.5e-1 | 9.5e-3 |
| 5 | LacLys -> CML | min-1 | 1.7e-4 | 2.6e-3 | 3.6e-3 | 1.4e-2 | 1.4e-3 | 2.7e-3 | 1.1e-3 | 2.3e-3 |
| 6 | LacLys -> CEL | min-1 | 0.0 | 0.0 | 1.0e-3 | 2.8e-3 | 3.3e-3 | 2.2e-4 | 2.1e-3 | 1.3e-3 |
| 7 | LacLys -> GO | min-1 | 2.5e-4 | ind | 0.0 | 0.0 | 0.0 | 0.0 | 2.0e-3 | 1.8e-3 |
| 8 | 1-DG -> DA | min-1 | 1.1e-1 | 1.2e-2 | 6.5e-2 | 1.1e-2 | 2.7e-1 | 2.6e-2 | 2.0e-1 | 2.8e-2 |
| 9 | 1-DG -> MGO | min-1 | 10.8 | 7.6e-1 | 25.7 | 1.2 | 8.1e-1 | 9.5e-2 | 123.5 | 6.7 |
| 10 | G -> GO | min-1 | 7.4e-2 | 1.5 | 3.3e-1 | 3.9e-2 | 3.5e-1 | 3.7e-1 | 2.4e-3 | ind |
| 11 | GO + Lys -> CML | kg umol-1 min-1 | 0.0 | 0.0 | 8.2e-7 | 2.5e-6 | 1.9e-6 | 5.0e-6 | 2.4e-6 | 4.4e-6 |
| 12 | MGO + Lys -> CEL | kg umol-1 min-1 | 8.9e-7 | 4.5e-7 | 4.1e-6 | 4.4e-6 | 0.0 | 0.0 | 0.0 | 0.0 |
| 13 | Lac -> Int | min-1 | 1.6e-2 | 3.7e-3 | 4.0e-1 | ind | 9.0e-1 | ind | 8.2e-1 | 9.0e-1 |
| 14 | Int -> Lac | min-1 | 0.0 | 0.0 | 6.4 | 3.1 | 5.7 | 0.6 | 4.8 | 5.5 |
| 15 | Lac -> P1 | min-1 | 4.1e-6 | ind | 1.2e-2 | 9.1e-3 | 0.0 | 0.0 | 0.0 | 0.0 |
| 16 | LacLys -> P2 | min-1 | 0.0 | 0.0 | 3.2e-2 | 6.0e-2 | 0.0 | 0.0 | 2.9e-1 | 1.6e-1 |
| 17 | 3-DG -> P3 | min-1 | 9.2e-2 | 1.4e-1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 18 | 1-DG -> P4 | min-1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 9.5e-2 | ind |
| 19 | G -> P5 | min-1 | 6.5 | ind | 11.4 | ind | 7.4 | 22.4 | 45.9 | ind |
| 20 | MGO -> P6 | min-1 | 8.0 | ind | 16.1 | 1.0 | 0.0 | 0.0 | 42.8 | ind |
| 21 | GO -> P7 | min-1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 22 | CML -> P8 | min-1 | 2.6 | 2.8 | 6.0 | 24.1 | 2.2 | 2.6 | 2.5 | 1.2 |
| 23 | CEL -> P9 | min-1 | 2.0 | 1.1 | 14.5 | ind | 11.4 | ind | 5.5 | 3.5 |
| 24 | Lys -> P10 | min-1 | 3.7e-2 | 6.0e-3 | 5.8e-2 | 8.0e-3 | 1.2e-1 | 8.6e-3 | 1.9e-1 | 1.8e-2 |

The text says nine constants (k3, k7, k10, k13, k15, k18, k19, k20, k23) "were not determined in
the +/- 95 % HPD interval at all temperatures" and were kept because removing them did not improve
the fit.

### Table 2. "Activation Energies (Ea) According to the Proposed Kinetic Model in Figure 2" (reparameterised Arrhenius, T_ref = 120 C)

| # | step | Ea (kJ/mol) | HPD |
|---|---|---:|---:|
| 1 | Lac + Lys -> LacLys | 52.1 | 19.0 |
| 2 | LacLys -> 3-DG | 64.1 | 19.8 |
| 3 | LacLys -> 1-DG | 89.8 | ind |
| 4 | LacLys -> G | 75.9 | 21.1 |
| 5 | LacLys -> CML | 75.4 | 17.6 |
| 6 | LacLys -> CEL | 61.5 | 19.9 |
| 7 | LacLys -> GO | 43.8 | 34.7 |
| 8 | 1-DG -> DA | 18.9 | 7.9 |
| 9 | 1-DG -> MGO | 51.0 | 20.1 |
| 10 | G -> GO | 4.2 | 15.7 |
| 11 | GO + Lys -> CML | -22.4 | ind |
| 12 | MGO + Lys -> CEL | 32.6 | ind |
| 13 | Lac -> Int | 14.0 | ind |
| 14 | Int -> Lac | -15.4 | 15.4 |
| 15 | Lac -> P1 | -11.3 | 42.3 |
| 16 | LacLys -> P2 | 17.5 | ind |
| 17 | 3-DG -> P3 | -205.8 | ind |
| 18 | 1-DG -> P4 | -12.5 | 867.6 |
| 19 | G -> P5 | 5.4 | ind |
| 20 | MGO -> P6 | 43.9 | 19.7 |
| 21 | GO -> P7 | -0.3 | ind |
| 22 | CML -> P8 | 43.1 | 23.1 |
| 23 | CEL -> P9 | 9.7 | 17.7 |
| 24 | Lys -> P10 | 23.4 | 3.7 |

**Arithmetic (mine).** Two-point barriers from the Table 1 end members, 110 -> 140 C: k1 106.6
kJ/mol (Table 2: 52.1 — the Table 2 value is a four-point global fit dominated by the 120/130 C
pair, which are almost equal); k2 24.6 (Table 2: 64.1); k4 90.6 (75.9); k9 106.9 (51.0); k24 71.8
(23.4). The Table 2 barriers are not reproducible from Table 1 by a simple two-point rule; the
paper's own words: "The wide 95 % HPD intervals indeed point to limited parameter identifiability,
which may be attributed to the relatively narrow experimental temperature range." Store Table 2 as
printed, with its HPDs; do not refit.

**Flux arithmetic at 120 C (mine, from Table 1 and the printed levels):** LacLys formation k1 x
[Lac] x [Lys] = 7.5e-9 x 939 150 x 194 670 = 1.4e3 umol/kg dm/min; 3-DG formation k2 x [LacLys] ~
5.9e-3 x 3e3 (LacLys of a few thousand umol/kg dm per the Figure 3c axis) ~ 18 umol/kg dm/min; GO
from glucosone k10 x [G] ~ 0.33 x 10 ~ 3 umol/kg dm/min with no sink (k21 = 0); MGO: k9 x [1-DG]
~ 25.7 x 12 ~ 3e2 umol/kg dm/min against a sink k20 = 16.1 /min, giving a steady state of ~ 19
umol/kg dm. **k9 and k20 are identified only as a ratio** (k9/k20 = 1.35, 1.60, —, 2.89 at
110/120/130/140 C; at 130 C k20 is zero and k9 drops thirty-fold — a different solution branch, see
Flag 2). The same holds for glucosone (k19 = 6.5 to 45.9 /min against k4 supply): with 30 s to
1 min sampling, any first-order constant above ~ 5 /min is at quasi-steady state and unresolvable.

### Appendix A (the balances, transcribed from the broken text layer; consistent with Figure 2)

d[LacLys]/dt = k1[Lac][Lys] - (k2 + k3 + k4 + k5 + k6 + k7 + k16)[LacLys]
d[Lac]/dt = k14[Int] - (k13 + k15)[Lac] - k1[Lac][Lys]
d[G]/dt = k4[LacLys] - (k10 + k19)[G]
d[3-DG]/dt = k2[LacLys] - k17[3-DG]
d[GO]/dt = k10[G] + k7[LacLys] - k11[GO][Lys] - k21[GO]
d[CML]/dt = k5[LacLys] + k11[GO][Lys] - k22[CML]
d[CEL]/dt = k6[LacLys] + k12[MGO][Lys] - k23[CEL]
d[MGO]/dt = k9[1-DG] - k12[MGO][Lys] - k20[MGO]
d[1-DG]/dt = k3[LacLys] - (k8 + k9 + k18)[1-DG]
d[DA]/dt = k8[1-DG]
d[Int]/dt = k13[Lac] - k14[Int]
d[Lys]/dt = (k2 + k3 + k4 + k7)[LacLys] - k1[Lac][Lys] - k11[GO][Lys] - k12[MGO][Lys] - k24[Lys]
d[Pn]/dt = the matching k[precursor] for P1-P10.

Notes on the scheme: lysine is **regenerated** by the LacLys -> 3-DG / 1-DG / G / GO steps (as in
Martins 2005) but not by the CML/CEL steps; 3-DG has **no** route to MGO or GO (excluded by
discrimination: "MGO formation was incorporated only through the 1-DG pathway"; "GO from G instead
of from 3-DG"); glucose and galactose are absent ("not experimentally detected"; their inclusion
"resulted in poor model fits"); lactulose, lactosone, 1-deoxylactosone were tried as hypothetical
intermediates and dropped; no acid, no HMF, no furfural, no Strecker product; pH change (0.1 in 5
min) deemed negligible.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): CML -> `cml`; CEL -> `cel`; furosine -> `furosine`;
diacetyl -> `2_3_butanedione`; **glyoxal, methylglyoxal, 3-deoxyglucosone, 1-deoxyglucosone,
glucosone, lactulosyllysine, lactose -> not in registry**.

| step | quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|---|
| Lac + Lys(total) -> LacLys | k at 110/120/130/140 C | 2.2e-9 / 7.5e-9 / 8.1e-9 / 2.5e-8 | kg dm umol-1 min-1 | whole milk, pH ~ 6.7, 0.5-5 min | Table 1 step 1 | measured_rate (fitted, multiresponse), matrix = milk |
| Lac + Lys(total) -> LacLys | k re-expressed on a milk-water basis (mine) | 1.5e-5 / 5.0e-5 / 5.4e-5 / 1.7e-4 | L mmol-1 min-1 on a milk-water basis (x 6.7e3) | same | derived | derived_assumption; trunk `k_schiff` (glucose + glycine) extrapolates to 7.8e-5 at 120 C |
| Lac + Lys -> LacLys | Ea | 52.1 +/- 19.0 | kJ/mol | 110-140 C | Table 2 | measured_barrier (four temperatures, wide HPD); trunk carries 96.8 |
| LacLys -> 3-DG | k | 4.4e-3 / 5.9e-3 / 6.8e-3 / 7.7e-3 | min-1 | same | Table 1 step 2 | measured_rate (fitted); cf. trunk `k_ama_tdg` 1.1e-2 at 100 C |
| LacLys -> 3-DG | Ea | 64.1 +/- 19.8 | kJ/mol | same | Table 2 | measured_barrier; trunk 97.1 +/- 1.7; Yu 2020 (lysine) 83.4 +/- 5.9 |
| LacLys -> 1-DG | k; Ea | 4.1e-2 / 6.9e-2 (ind) / 2.9e-3 / 1.9e-1; 89.8 (ind) | min-1; kJ/mol | same | Table 1 step 3, Table 2 | measured_rate, HPD indeterminate at 120 C, 130 C value off-trend; barrier indeterminate; trunk `k_ama_odg` 1.6e-2 at 100 C, Ea 107.3 |
| LacLys -> glucosone | k; Ea | 1.9e-2 / 2.1e-2 / 2.5e-2 (+/- 7.1e-2) / 1.5e-1; 75.9 +/- 21.1 | min-1; kJ/mol | same | Table 1 step 4, Table 2 | measured_rate / measured_barrier — **the only aqueous Amadori -> glucosone entry on disk** (trunk `k_glc_g` is glucose -> glucosone in a 160-200 C glass, Ea 125.9) |
| glucosone -> glyoxal | k; Ea | 7.4e-2 (+/- 1.5) / 3.3e-1 / 3.5e-1 (+/- 3.7e-1) / 2.4e-3 (ind); 4.2 +/- 15.7 | min-1; kJ/mol | same | Table 1 step 10, Table 2 | measured_rate, poorly determined; barrier consistent with zero; trunk `k_g_go` 0.737 /min at 180 C, Ea 93.8 |
| LacLys -> glyoxal (direct) | k; Ea | 2.5e-4 (ind) / 0 / 0 / 2.0e-3; 43.8 +/- 34.7 | min-1; kJ/mol | same | Table 1 step 7, Table 2 | measured_rate, essentially unidentified |
| glucosone -> GO vs LacLys -> GO | k10 / k7 | ~ 100 (110 C) -> ~ 1 (140 C) | ratio of first-order constants (text: "7.4e-2 and 2.5e-4"; "2.0e-3 and 2.4e-3") | same | text 3.2.2 | within_study_ratio (of two ill-determined constants) |
| 1-DG -> methylglyoxal | k; Ea | 10.8 / 25.7 / 8.1e-1 / 123.5; 51.0 +/- 20.1 | min-1; kJ/mol | same | Table 1 step 9, Table 2 | measured_rate — but identified only jointly with k20 (Flag 2) |
| 1-DG -> diacetyl | k; Ea | 1.1e-1 / 6.5e-2 / 2.7e-1 / 2.0e-1; 18.9 +/- 7.9 | min-1; kJ/mol | same | Table 1 step 8, Table 2 | measured_rate / measured_barrier; trunk `k_odg_da` 12.2e-3 /min at 180 C, Ea 150.8 (dry glass) |
| methylglyoxal sink (-> P6) | k; Ea | 8.0 (ind) / 16.1 / 0 / 42.8 (ind); 43.9 +/- 19.7 | min-1; kJ/mol | same | Table 1 step 20, Table 2 | measured_rate, jointly with k9 |
| **glyoxal sink (-> P7)** | k | **0.0 +/- 0.0 at 110, 120, 130 and 140 C**; Ea -0.3 (ind) | min-1 | same | Table 1 step 21 | level_only (fitted zero): in milk at 110-140 C the data want **no glyoxal loss** beyond the GO + Lys -> CML step; the trunk's `k_go_sink` is 32.6e-3 /min at 180 C with Ea fixed to 0 |
| 3-DG sink (-> P3) | k | 9.2e-2 (+/- 1.4e-1) / 0 / 0 / 0; Ea -205.8 (ind) | min-1 | same | Table 1 step 17 | level_only (fitted zero at 120-140 C) |
| GO + Lys -> CML | k; Ea | 0 / 8.2e-7 / 1.9e-6 / 2.4e-6; -22.4 (ind) | kg dm umol-1 min-1 | same | Table 1 step 11 | measured_rate, HPD wider than value; "1000-10 000 times lower" than the LacLys route |
| MGO + Lys -> CEL | k | 8.9e-7 / 4.1e-6 / 0 / 0; Ea 32.6 (ind) | kg dm umol-1 min-1 | same | Table 1 step 12 | measured_rate at 110-120 C only |
| LacLys -> CML; LacLys -> CEL | k; Ea | 1.7e-4 / 3.6e-3 / 1.4e-3 / 1.1e-3 (75.4 +/- 17.6); 0 / 1.0e-3 / 3.3e-3 / 2.1e-3 (61.5 +/- 19.9) | min-1; kJ/mol | same | Table 1 steps 5, 6 | measured_rate / measured_barrier; the dominant AGE routes in this fit |
| Lys -> P10 (lysine loss not through LacLys) | k; Ea | 3.7e-2 / 5.8e-2 / 1.2e-1 / 1.9e-1; 23.4 +/- 3.7 | min-1; kJ/mol | same | Table 1 step 24 | measured_rate / measured_barrier (the best-determined barrier in the table) |
| levels | dicarbonyls in heated milk | 3-DG 51.17, 1-DG 11.79, MGO 10.12 umol/kg dm (110 C 5 min); GO 14.07 umol/kg dm (120 C 5 min) | umol/kg dm (= 0.91 / 0.21 / 0.08 / 0.1 mg/L milk) | as stated | text 3.1 | level_only (the only printed points; all time courses figure_only) |
| all | mass balance recovery | 89.8 / 86.1 / 82.5 / 80.6 % | % of initial moles | 110-140 C, end of run | text 3.2 | level_only |
| all | all species vs time at four temperatures | — | umol/kg dm (axis) | 0-5 min | Figure 3 | figure_only |

**What transfers to the trunk, and what does not.**

- Transfers as a comparator: (i) the Amadori-degradation barriers in an aqueous matrix at 110-140 C
  — LacLys -> 3-DG 64 +/- 20, -> glucosone 76 +/- 21, -> 1-DG 90 (ind) kJ/mol — against the trunk's
  97 / — / 107 and Yu 2020's 83 (3-DG); (ii) the finding that an aqueous multiresponse fit at these
  temperatures **wants no glyoxal sink** and no 3-DG sink above 110 C; (iii) the fit's preference
  for methylglyoxal from 1-DG rather than 3-DG, and glyoxal from glucosone rather than 3-DG, both
  of which agree with the routes the trunk already carries (`k_ama_mgo` is a DFG -> MG lump in
  Martins; `k_g_go` in the B13 lane); (iv) an aqueous Amadori -> glucosone constant at all, since the
  trunk's glucosone entry is a dry-glass extrapolation; (v) the absolute dicarbonyl levels in a
  neutral aqueous food (1-15 umol/L) as an order-of-magnitude anchor.
- Does not transfer: any constant as a number. The sugar is lactose (a 1,4-disaccharide whose
  Amadori product degrades through routes glucose lacks — the paper itself notes 4-DG shares m/z
  with 3-DG and cannot be excluded), the amine is casein-bound lysine at 195 mmol/kg dm with free
  lysine at ~ 50 umol/kg, the basis is per kg dry matter, the temperature range is 30 K with five
  points each, and the constants above ~ 5 /min (k9, k14, k19, k20, k22, k23) are quasi-steady-state
  artefacts of 30 s sampling. The second-order Amadori constant re-expressed on a water basis
  (5.0e-5 L mmol-1 min-1 at 120 C) lands within a factor 1.6 of the trunk's glucose + glycine value
  extrapolated to 120 C (7.8e-5) — a coincidence worth recording, not a calibration.

## 5. Flags

1. **Basis.** Concentrations are umol per kg of **dry matter** (Figure 3 axes; the lactose and
   lysine totals confirm it), not per kg of milk as one Methods sentence says. Second-order
   constants scale with the basis; first-order constants do not. My dm fraction (0.130) is derived
   from the paper's own lactose pair, not printed.
2. **Fast first-order constants are not resolved.** The 130 C column is a different solution
   branch: k9 (1-DG -> MGO) falls from 25.7 to 0.81 /min while k20 (MGO sink) goes to zero and k3
   (LacLys -> 1-DG) falls twenty-fold — the fit moved the bottleneck upstream. Likewise k19 (G sink,
   6.5-45.9 /min, ind) against k4; k22/k23 (CML/CEL sinks, 2-14 /min) with CML/CEL at umol/kg
   levels. Use ratios (k9/k20, k10/k19) or the barriers, not the constants.
3. **Barriers are weakly identified**: four temperatures over 30 K, fitted separately then joined;
   HPDs of 16-35 kJ/mol on the important steps; six negative Ea; the authors say so. Two-point
   end-member barriers (mine) disagree with Table 2 by 30-55 kJ/mol for k1, k2, k9, k24.
4. **Glyoxal sink fitted to exactly zero at all four temperatures** and 3-DG sink zero at 120-140 C:
   this is a statement about identifiability within 5 min as much as about chemistry (GO rises
   monotonically over 5 min at every temperature per Figure 3i, so no loss term is needed). It is
   nonetheless the only aqueous 110-140 C evidence on disk bearing on `k_go_sink`.
5. **1-DG semi-quantitative** (3-DG calibration); the authors note the largest fit deviations are for
   1-DG. 4-deoxyglucosone is isobaric with 3-DG and not separated.
6. **Lysine is total lysine** after 23 h acid hydrolysis, so "Lys" includes LacLys's lysine that
   hydrolysis releases as furosine + regenerated lysine; the paper does not say whether the
   furosine-derived fraction was subtracted. Free lysine (~ 7.8 ug/mL in milk per the introduction)
   is negligible against bound.
7. **LacLys = furosine x 3.1** (Brands & van Boekel's factor for lactulosyllysine hydrolysis yield);
   every LacLys-consuming constant inherits that factor.
8. **No glucose/galactose, lactulose, acids, HMF, furfural or Strecker products in the scheme**;
   pH change taken as 0.1 and ignored.
9. **Text-vs-table**: "16.10 % in the milk heated at 140 C for 3 min" but the 140 C series ends at
   2.5 min; the CML "roasting temperature" phrase is a copy from the group's nut papers.
10. **SI not on disk** (mass balances, correlation matrices); the paper says 120 C has the fewest
    strong parameter correlations.
11. Registry gap: glyoxal, methylglyoxal, 3-deoxyglucosone, 1-deoxyglucosone, glucosone and
    lactulosyllysine have no molecule row in `compounds.yml`; diacetyl is `2_3_butanedione`.
