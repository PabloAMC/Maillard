# De Vleeschouwer, Van der Plancken, Van Loey & Hendrickx 2008 — EXTRACTION (the a_w tables)
### Acrylamide formation/elimination in equimolar asparagine-glucose powders at four initial water activities, 120-200 C, multiresponse-fitted.

**Source on disk:** `data/articles/devleeschouwer2008.pdf` (owner's download, 2026-09-07). Read-only
extraction, wave B14. Text layer clean; Tables 1-3 re-typed from it below.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Investigation of the Influence of Different Moisture Levels on Acrylamide Formation/Elimination Reactions Using Multiresponse Analysis" |
| Authors | Kristel De Vleeschouwer, Iesel Van der Plancken, Ann Van Loey, Marc E. Hendrickx (KU Leuven, LFoRCe) |
| Venue | J. Agric. Food Chem. 2008, 56, 6460-6470 |
| DOI | 10.1021/jf8006294 |

## 1. Why the repo needed it

`parameters_acrylamide.py` carries `k_int1_acr` = 3.57e-3 +/- 1.38e-3 /min, Ea 159.2 +/- 29.5 kJ/mol
and `k_asn_glc` = 1.70 +/- 1.05 /M/min, Ea 117.5 +/- 25.2 at a_w 0.92, cited through De Vleeschouwer
2009 Part I (Food Chem 114:116, Table 3, glucose column). **Those are the a_w 0.92 column of this
paper's Table 2** (values identical to the printed digit), so the lane's constants come from the
equimolar Asn-Glc powder equilibrated at a_w 0.92, and this paper is the one that measures the same
constants at three other water activities.

## 2. Methods, as they matter to the lane

- Equimolar L-asparagine + D-glucose, freeze-dried, equilibrated at 4 C over saturated salts
  (Table 1: KCl 88 %, Sr(NO3)2 92 %, KNO3 96 %, K2SO4 99 % ERH); water contents 10.4 / 14.5 / ... %
  (basic system), 18.5 / 19.2 / ... % (potato-based, 10 % equimolar mixture in potato powder,
  glucose initially 5.27 g/kg dw).
- Heated (~50 mg) in hermetically closed custom Inox reactor tubes (8 x 100 mm) in a thermostated oil
  bath at 120 / 140 / 160 / 180 / 200 C; the registered temperature-time profile of each tube was
  integrated numerically in the fit (non-isothermal correction).
- Responses fitted together: acrylamide, glucose, asparagine, aspartic acid, browning (470 nm);
  Scheme 2 (the modified network): Asn + Glc -> INT (k_INT), INT -> acrylamide (k_F), acrylamide
  elimination (k_E), INT -> melanoidin (k_M), browning (k_B), k_C, Asn -> Asp (k_Asp), k_X.
  Determinant criterion, 95 % HPD intervals. T_ref 160 C.

## 3. Table 2 — basic Asn-Glc system, T_ref 160 C (value +/- 95 % HPD)

| parameter | a_w 0.88 | 0.92 | 0.96 | 0.99 |
|---|---:|---:|---:|---:|
| k_Fref (1e-3 /min) | 2.29 +/- 0.43 | **3.57 +/- 1.38** | 3.45 +/- 1.19 | 1.45 +/- 0.42 |
| k_Eref (/min) | 0.11 +/- 0.02 | 0.10 +/- 0.04 | 0.09 +/- 0.04 | 0.05 +/- 0.03 |
| k_INTref (/M/min) | 1.11 +/- 0.61 | **1.70 +/- 1.05** | 1.42 +/- 0.59 | 1.43 +/- 0.45 |
| k_Mref (/min) | 0.58 +/- 0.15 | 1.23 +/- 0.49 | 0.56 +/- 0.29 | 0.38 +/- 0.17 |
| k_Bref (/M/min) | 4.11 +/- 10.77 | 3.90 +/- 3.68 | 3.12 +/- 3.87 | 0.90 +/- 0.49 |
| k_Cref (/min) | indet. | indet. | 0.05 +/- 0.05 | indet. |
| k_Aspref (1e-3 /min) | 15.77 +/- 3.36 | 26.43 +/- 5.76 | 20.33 +/- 4.35 | 22.11 +/- 6.51 |
| k_Xref (1e-3 /min) | 5.64 +/- 4.06 | indet. | 0.70 +/- 1.97 | indet. |
| Ea_F (kJ/mol) | 145.5 +/- 15.9 | **159.2 +/- 29.5** | 135.8 +/- 27.3 | 84.6 +/- 17.9 |
| Ea_E | 102.5 +/- 15.4 | 113.2 +/- 32.3 | 101.8 +/- 34.5 | 143.4 +/- 34.3 |
| Ea_INT | 95.7 +/- 21.6 | **117.5 +/- 25.2** | 114.7 +/- 28.8 | 120.8 +/- 13.2 |
| Ea_M | 92.2 +/- 22.2 | 105.7 +/- 29.1 | 102.0 +/- 36.9 | 40.0 +/- 21.1 |
| Ea_B | 44.5 +/- 17.5 | 180.3 +/- 38.5 | 133.0 +/- 50.6 | 124.5 +/- 24.7 |
| Ea_C | -4.0 +/- 0.2 | -6.7 +/- 0.2 | 7.4 +/- 6.6 | -5.6 +/- 0.1 |
| Ea_Asp | 71.8 +/- 6.7 | 105.4 +/- 10.6 | 83.5 +/- 10.6 | 98.4 +/- 14.4 |
| Ea_X | 159.4 +/- 29.8 | 668.9 +/- 35.2 | 207.4 +/- 120.1 | indet. |

Bold = the lane's shipped constants. Authors' reading: formation and elimination parameters
"did not change significantly (based on a 95% confidence level) within the range of water
activities tested". The k_Fref point estimates relative to the 0.92 column: 0.64 / 1.00 / 0.97 / 0.41.

## 4. Table 3 — potato-based system, T_ref 160 C

| parameter | a_w 0.88 | 0.92 | 0.96 | 0.99 |
|---|---:|---:|---:|---:|
| k_Fref (1e-3 /min) | 2.51 +/- 1.70 | 1.96 +/- 0.50 | 0.83 +/- 1.16 | 2.40 +/- 0.94 |
| k_Eref (/min) | 0.07 +/- 0.04 | 0.05 +/- 0.01 | 0.08 +/- 0.04 | 0.04 +/- 0.02 |
| k_INTref (/M/min) | 3.94 +/- 3.46 | 3.09 +/- 0.52 | 2.82 +/- 0.69 | 3.50 +/- 0.90 |
| k_Mref (/min) | 0.96 +/- 0.70 | 0.59 +/- 0.17 | 0.17 +/- 0.03 | 0.58 +/- 0.24 |
| k_Bref (/M/min) | 0.61 +/- 0.33 | 3.00 +/- 1.46 | 4.37 +/- 2.54 | 1.76 +/- 0.80 |
| k_Aspref (1e-3 /min) | 28.27 +/- 8.00 | 15.71 +/- 4.11 | 19.45 +/- 4.63 | 36.45 +/- 8.21 |
| k_Xref (1e-3 /min) | 1.53 +/- 3.35 | 0.08 +/- 0.36 | indet. | 0.51 +/- 0.88 |
| Ea_F (kJ/mol) | 157.1 +/- 44.6 | 147.5 +/- 16.4 | 89.9 +/- 11.9 | 171.5 +/- 27.3 |
| Ea_E | 64.5 +/- 35.4 | 70.6 +/- 17.4 | 22.4 +/- 31.7 | 180.5 +/- 40.6 |
| Ea_INT | 109.2 +/- 36.4 | 105.3 +/- 7.5 | 112.3 +/- 10.7 | 125.7 +/- 11.3 |
| Ea_M | 117.0 +/- 43.3 | 121.0 +/- 19.4 | 121.1 +/- 8.0 | 118.2 +/- 26.0 |
| Ea_B | 77.2 +/- 33.2 | 114.2 +/- 29.5 | 197.0 +/- 44.0 | 129.5 +/- 30.7 |
| Ea_C | -23.8 +/- 1.2 | -21.7 +/- 0.33 | -22.2 +/- 0.6 | -12.0 +/- 0.51 |
| Ea_Asp | 93.2 +/- 15.5 | 99.0 +/- 14.3 | 87.0 +/- 11.6 | 104.0 +/- 11.6 |
| Ea_X | 357.9 +/- 160.7 | 426.0 +/- 242.5 | indet. | 475.8 +/- 149.3 |

Authors: the potato matrix changes neither the formation nor the elimination constants
significantly; Ea_E tends lower with the matrix.

## 5. What the repo takes (wave B14) and what it does not

- TAKEN: a DECLARED FLAT a_w term on `k_int1_acr` inside the measured window 0.88-0.99, band
  (0.41, 1.39) on the multiplier (`acrylamide_conditions.py`): the four point estimates relative to
  the shipped column united with that column's HPD. The window boundary is a refusal boundary.
- NOT taken: any a_w dependence below 0.88 (unmeasured here; the dry-side extrusion claims stay
  refused); the elimination constants (the lane's `k_acr_dp` is a declared band, not this source's);
  Table 3 (a validation set for a potato matrix the lane does not represent).
- Recorded caveat: the fit was non-isothermal (registered profiles), so the printed T_ref constants
  are model-integrated, not read off a plateau.
