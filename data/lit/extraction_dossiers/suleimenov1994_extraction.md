# Suleimenov & Krupp 1994 — EXTRACTION
### *"Solubility of hydrogen sulfide in pure water and in NaCl solutions, from 20 to 320°C and at saturation pressures"*
### Henry constants of H2S from P-V-T-X measurements in a 55 mL titanium bellows, 20-321 °C. log KH (bar per MOLE FRACTION) = 3.141 at 90 °C, 3.239 at 120 °C, 3.287 at 145 °C (Table 2, µ = 0). As a dimensionless gas/liquid concentration ratio: 0.86, 1.01, 1.09. Read 2026-10-29.

**Full read of `data/articles/suleimenov1994.pdf`** (12 pp., scan with OCR text layer; the OCR garbles
symbols, so the abstract, the correlation and Fig. 3's axis label were checked on 150-dpi renders of
pp. 2433 and 2437; Table 2 and Appendix A numbers were read from the layout text and spot-checked).

**Provenance:** **[M]** printed value · **[Z]** derived by me from printed values, never printed ·
**[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | O. M. Suleimenov and R. E. Krupp |
| Affiliation | Bundesanstalt für Geowissenschaften und Rohstoffe, Hannover |
| Venue | *Geochim. Cosmochim. Acta* **58**(11), 2433-2444 (1994), Pergamon |
| DOI | **not printed**; the page prints the code "0016-7037(94)00048-4" |
| Dates | received 25 Aug 1993, accepted in revised form 19 Jan 1994 |

---

## 1. THE ONE-PARAGRAPH ANSWER

H2S solubility in pure water falls with temperature to a minimum near 185 °C ("about 7.5 times less
than at 0°C") and then rises. The fitted Henry constants, converted to a gas/liquid concentration
ratio K = c_gas / c_aq(H2S, neutral), are **0.38 (25 °C), 0.86 (90 °C), 1.01 (120 °C), 1.09 (145 °C)**.
Above ~100 °C neutral H2S is therefore as concentrated in a headspace as in the liquid. For a 15 mL
vial with 5 mL liquid at 120 °C, **67 % of neutral H2S is in the gas; 31 % of total sulfide at pH 7**
(HS⁻ stays in the liquid; pK1 6.45 at 120 °C).

---

## 2. SYSTEM AND CONDITIONS [M], quoted

- **Vessel:** "a titanium pressure vessel ... subdivided by a titanium membrane bellows (volume 55 mL)
  containing the H2S-H2O vapor-liquid, 2-phase system". "All seals and gaskets wetted by H2S were made
  of gold." Total volumes per run 50.8-58.5 mL (Appendix A).
- **Temperature:** oil-immersed sheathed thermocouple; gradients "less than 2°C across the autoclave
  body"; accuracy "better than ±0.2°C up to about 200°C"; linear correction above 200 °C (max 1.0 °C
  at 320 °C).
- **Pressure:** piezo-resistive sensor, 0.1 % class, 250 bar; screw press holds the bellows relaxed
  (one turn = 0.15 mL → 0.01 bar, no hysteresis).
- **Water:** distilled, deionised (≤ 0.1 µS/cm), "degassed by boiling for a minimum of 30 min under a
  partial vacuum".
- **H2S:** "generated from analytical grade sodium sulfide crystals and concentrated orthophosphoric
  acid (in some cases Al2S3 + H2O) in an evacuated Kipp-type gas generator", dried (ice + CaCl2,
  Mg(ClO4)2), frozen at liquid-N2 temperature, non-condensables pumped off.
- **Filling:** evacuated bellows; "a weighed amount of degassed water was transferred under its own
  vapor pressure at 80-90°C ... followed by introduction of H2S gas"; quantities by weight difference.
- **Bulk composition:** ~1.0-1.67 mol H2O and 0.0072-0.0299 mol H2S per run (Appendix A).
- **pH / buffer:** none; pure water or NaCl (0-2.5 molal). H2S essentially un-ionised (pure water).
- **Equilibrium:** "Attainment of thermal and chemical equilibrium was demonstrated by means of a chart
  recorder" (times not given).

## 3. ANALYTICAL METHOD [M]

After the run the H2S was neutralised in situ by NaOH via an HPLC pump; after cooling, "the fully
neutralized H2S solution was immediately sampled and analysed for its H2S content by iodometric
titration ..., with typical uncertainties of 0.2%". Henry constants by "liquid-vapor flash
calculations", PRSV equation of state for the vapour (Carroll & Mather mixing rules; k_ij given),
a Pitzer molecule-molecule term for the liquid (temperature dependence from Edwards et al. 1975; the
OCR of its formula is garbled and it is not transcribed), Poynting correction. Uncertainty by Monte Carlo (> 10,000 draws per point): "typically ±0.01 log unit or less".
Pressure 2.2-138.6 bar.

---

## 4. RESULTS [M]

### 4a. Correlation (abstract and p. 2437)

log KH(T, P_sat,1) = A + B·T + C·T² + D/T + E·log T, 293 < T < 593 K, T in K:
A = +0.6342702616e+3, B = +0.2709284796e+0, C = −0.1113202904e−3, D = −0.1671907660e+5,
E = −0.2619219571e+3. **Units printed: "bar/molality".**

**[!] Unit check [Z].** The correlation evaluated at 298.15 K gives KH = 9.30 → log 0.968, which equals
Table 2's log KH(25 °C) = 2.713 plus log(0.0180153 kg/mol) = −1.744. So **the correlation is in bar
kg/mol (molality), while Table 2, Appendix A and Fig. 3 ("Log KH (bar/mol fr.)") are in bar per mole
fraction.** Appendix A confirms: 23.5 °C, P = 2.95 bar, x(H2S) = 0.005483 → (2.95 − 0.03)/0.005483 =
533 bar/mol fr. vs printed log KH 2.713 (516; the difference is the fugacity correction). Both are
correct as labelled; mixing them is a factor-55 error.

### 4b. Table 2 — skeleton log KH (bar/mol fr.) at saturation pressure, 0-200 °C, µ = ionic strength (NaCl) [M]

| T °C | µ = 0 | µ = 1 | µ = 2 | µ = 3 |
|---|---|---|---|---|
| 0 | 2.431 | 2.514 | 2.598 | 2.682 |
| 5 | 2.494 | 2.576 | 2.657 | 2.738 |
| 10 | 2.554 | 2.634 | 2.713 | 2.792 |
| 15 | 2.611 | 2.688 | 2.765 | 2.843 |
| 20 | 2.664 | 2.739 | 2.815 | 2.891 |
| 25 | 2.713 | 2.788 | 2.862 | 2.937 |
| 30 | 2.760 | 2.833 | 2.907 | 2.980 |
| 35 | 2.804 | 2.877 | 2.949 | 3.022 |
| 40 | 2.845 | 2.917 | 2.989 | 3.061 |
| 45 | 2.884 | 2.956 | 3.027 | 3.098 |
| 50 | 2.921 | 2.992 | 3.063 | 3.134 |
| 55 | 2.955 | 3.026 | 3.097 | 3.167 |
| 60 | 2.987 | 3.058 | 3.128 | 3.199 |
| 65 | 3.017 | 3.088 | 3.159 | 3.229 |
| 70 | 3.045 | 3.116 | 3.187 | 3.258 |
| 75 | 3.072 | 3.143 | 3.214 | 3.284 |
| 80 | 3.096 | 3.167 | 3.238 | 3.309 |
| 85 | 3.119 | 3.191 | 3.262 | 3.333 |
| **90** | **3.141** | 3.212 | 3.284 | 3.355 |
| 95 | 3.161 | 3.232 | 3.304 | 3.376 |
| 100 | 3.179 | 3.251 | 3.323 | 3.395 |
| 105 | 3.196 | 3.268 | 3.340 | 3.412 |
| 110 | 3.212 | 3.284 | 3.356 | 3.428 |
| 115 | 3.226 | 3.298 | 3.371 | 3.443 |
| **120** | **3.239** | 3.311 | 3.384 | 3.456 |
| 125 | 3.251 | 3.323 | 3.396 | 3.468 |
| 130 | 3.262 | 3.334 | 3.406 | 3.478 |
| 135 | 3.271 | 3.343 | 3.416 | 3.488 |
| 140 | 3.280 | 3.352 | 3.424 | 3.496 |
| **145** | **3.287** | 3.359 | 3.431 | 3.502 |
| 150 | 3.294 | 3.365 | 3.437 | 3.508 |
| 155 | 3.299 | 3.370 | 3.441 | 3.512 |
| 160 | 3.304 | 3.374 | 3.445 | 3.516 |
| 165 | 3.307 | 3.378 | 3.448 | 3.518 |
| 170 | 3.310 | 3.380 | 3.450 | 3.520 |
| 175 | 3.312 | 3.381 | 3.451 | 3.520 |
| 180 | 3.313 | 3.382 | 3.451 | 3.520 |
| 185 | 3.313 | 3.382 | 3.451 | 3.519 |
| 190 | 3.312 | 3.381 | 3.449 | 3.518 |
| 195 | 3.311 | 3.379 | 3.448 | 3.516 |
| 200 | 3.309 | 3.377 | 3.446 | 3.514 |

(205-365 °C printed in the same table; not transcribed — outside any food process.)

### 4c. Appendix A — primary data, pure water [M] (T °C, P bar, V mL, n H2O mol, n H2S mol, log KH ± σ, x H2S, y H2S)

| T | P | V | nH2O | nH2S | log KH | xH2S | yH2S |
|---|---|---|---|---|---|---|---|
| 23.5 | 2.95 | 57.9 | 1.6602 | 0.0125 | 2.713 ± .006 | 0.005483 | 0.990332 |
| 60.2 | 5.00 | 57.9 | 1.6602 | 0.0125 | 3.007 ± .005 | 0.004591 | 0.958000 |
| 100.8 | 7.44 | 57.9 | 1.6602 | 0.0125 | 3.176 ± .007 | 0.004145 | 0.854136 |
| 153.6 | 13.30 | 57.9 | 1.6602 | 0.0125 | 3.288 ± .009 | 0.003991 | 0.593201 |
| 206.2 | 27.07 | 57.9 | 1.6602 | 0.0125 | 3.324 ± .007 | 0.004115 | 0.328610 |
| 22.4 | 3.93 | 57.7 | 1.6657 | 0.0173 | 2.689 ± .004 | 0.007651 | 0.993125 |
| 151.0 | 16.03 | 57.7 | 1.6657 | 0.0173 | 3.287 ± .006 | 0.005501 | 0.679549 |
| 150.8 | 16.06 | 57.7 | 1.6657 | 0.0173 | 3.291 ± .006 | 0.005472 | 0.681716 |
| 35.6 | 5.56 | 53.3 | 1.663 | 0.0179 | 2.836 ± .004 | 0.007657 | 0.988772 |
| 234.3 | 44.98 | 53.3 | 1.663 | 0.0179 | 3.263 ± .004 | 0.007066 | 0.294868 |
| 167.1 | 20.45 | 53.3 | 1.663 | 0.0179 | 3.293 ± .005 | 0.006316 | 0.619299 |
| 156.2 | 18.04 | 53.3 | 1.663 | 0.0179 | 3.272 ± .005 | 0.006346 | 0.672040 |
| 142.6 | 15.82 | 53.3 | 1.663 | 0.0179 | 3.260 ± .005 | 0.006290 | 0.739948 |
| 140.8 | 15.91 | 53.3 | 1.663 | 0.0179 | 3.279 ± .005 | 0.006147 | 0.753517 |
| 20.8 | 4.25 | 50.8 | 1.248 | 0.0163 | 2.651 ± .005 | 0.008977 | 0.993857 |
| 232.4 | 43.24 | 50.8 | 1.248 | 0.0163 | 3.282 ± .007 | 0.006442 | 0.293170 |
| 254.8 | 56.89 | 52.8 | 1.248 | 0.0163 | 3.251 ± .012 | 0.006632 | 0.209280 |
| 22.3 | 4.96 | 52.8 | 1.050 | 0.01749 | 2.674 ± .004 | 0.009886 | 0.994173 |
| 242.6 | 49.46 | 52.8 | 1.050 | 0.01749 | 3.274 ± .007 | 0.006745 | 0.258729 |
| 273.0 | 73.33 | 52.8 | 1.050 | 0.01749 | 3.252 ± .009 | 0.007103 | 0.173177 |
| 242.5 | 49.29 | 52.8 | 1.050 | 0.01749 | 3.266 ± .007 | 0.006818 | 0.257846 |
| 234.7 | 44.81 | 52.8 | 1.050 | 0.01749 | 3.286 ± .005 | 0.006487 | 0.287884 |
| 221.1 | 37.42 | 52.8 | 1.050 | 0.01749 | 3.279 ± .006 | 0.006514 | 0.338558 |
| 198.2 | 28.45 | 52.8 | 1.050 | 0.01749 | 3.316 ± .008 | 0.006032 | 0.448500 |
| 85.8 | 8.66 | 52.8 | 1.230 | 0.01555 | 3.125 ± .005 | 0.005835 | 0.927107 |
| 155.1 | 16.36 | 52.8 | 1.230 | 0.01555 | 3.294 ± .016 | 0.005289 | 0.651042 |
| 22.8 | 5.41 | 53.4 | 1.0167 | 0.0185 | 2.693 ± .005 | 0.010297 | 0.994456 |
| 193.5 | 27.28 | 53.4 | 1.0167 | 0.0185 | 3.301 ± .005 (OCR ".OOS") | 0.006405 | 0.478679 |
| 229.6 | 42.59 | 53.4 | 1.0167 | 0.0185 | 3.292 ± .009 | 0.006691 | 0.315895 |
| 21.3 | 4.29 | 52.7 | 1.1716 | 0.0157 | 2.678 ± .006 | 0.008519 | 0.993719 |
| 291.8 | 89.38 | 52.7 | 1.1716 | 0.0157 | 3.112 ± .009 | 0.007747 | 0.109805 |
| 300.7 | 99.92 | 52.7 | 1.1716 | 0.0157 | 3.085 ± .009 | 0.008022 | 0.094903 |
| 310.9 | 112.86 | 52.7 | 1.1716 | 0.0157 | 3.030 ± .013 | 0.008561 | 0.077700 |
| 40.4 | 2.64 | 54.6 | 1.1426 | 0.00724 | 2.872 ± .005 | 0.000102 [!] | 0.970539 |
| 274.6 | 65.41 | 54.6 | 1.1426 | 0.00724 | 3.249 ± .031 | 0.002868 | 0.079152 |
| 40.0 | 4.53 | 54.4 | 1.2329 | 0.01313 | 2.851 ± .002 | 0.006065 | 0.982792 |
| 256.4 | 55.11 | 54.4 | 1.2329 | 0.01313 | 3.258 ± .006 | 0.005128 | 0.170784 |
| 266.0 | 62.63 | 54.4 | 1.2329 | 0.01313 | 3.235 ± .006 | 0.005319 | 0.147321 |
| 84.6 | 9.48 | 54.4 | 1.3698 | 0.01795 | 3.118 ± .003 | 0.006525 | 0.936083 |
| 196.2 | 27.99 | 54.4 | 1.3698 | 0.01795 | 3.312 ± .039 | 0.006176 | 0.462118 |
| 21.4 | 2.60 | 55.5 | 1.0580 | 0.0093 | 2.691 ± .010 | 0.005089 | 0.989826 |
| 155.4 | 11.59 | 55.5 | 1.0580 | 0.0093 | 3.289 ± .016 | 0.003028 | 0.516540 |
| 24.2 | 2.22 | 58.5 | 1.1206 | 0.00808 | 2.715 ± .004 | 0.004125 | 0.985968 |
| 199.8 | 21.12 | 58.5 | 1.1206 | 0.00808 | 3.332 ± .009 | 0.002475 | 0.257570 |
| 216.5 | 27.55 | 58.5 | 1.1206 | 0.00808 | 3.342 ± .012 | 0.002475 | 0.202746 |
| 22.8 | 7.90 | 57.1 | 1.1007 | 0.0299 | 2.667 ± .005 | 0.015546 | 0.996061 |
| 264.1 | 72.85 | 57.1 | 1.1007 | 0.0299 | 3.214 ± .007 | 0.011850 | 0.259991 |
| 279.1 | 86.35 | 57.1 | 1.1007 | 0.0299 | 3.178 ± .007 | 0.012592 | 0.211600 |
| 321.0 | 138.61 | 57.1 | 1.1007 | 0.0299 | 3.046 ± .009 | 0.015358 | 0.110385 |

[!] 40.4 °C, x = 0.000102 is inconsistent with its own P and log KH (≈ (2.64 − 0.07) / 745 ≈ 0.0034 expected);
a print or OCR error in x only; log KH 2.872 is in line with Table 2 (2.845 at 40 °C).

### 4d. Other printed results (not transcribed)

- Appendix B (NaCl solutions, primary data) and Table 1 (experimental activity coefficients of H2S
  in NaCl) — salting-out; Setchenov polynomial b(T, µ) coefficients printed on p. 2441
  (A = +0.8371062650e−01, B = −0.5135608863e−03, C = +0.6387039005e−05, D = −0.2217360319e−07,
  E = −0.5069412169e−10, F = +0.2827486651e−12), valid 25-320 °C. Salting-out "nearly independent of
  temperature up to about 250°C".
- Scaled-particle fit: σ₂ (Å) = 2.669 + 0.00296 T; ε₂/k (K) = 294.28 + 0.516417 T.
- ΔS of dissolution changes sign above 142 °C.

---

## 5. CONVERSION TO A HEADSPACE PARTITION COEFFICIENT [Z]

K_gl = c_gas / c_aq(H2S) = KH(bar/mol fr.) × 0.0180153 / (ρ_w · R · T), R = 0.0831446 L bar mol⁻¹ K⁻¹,
ρ_w in kg/L (ideal gas; Poynting and fugacity corrections neglected, < ~2 % at ≤ 145 °C and a few bar).

| T °C | log KH (bar/mol fr.) | KH | kH (bar kg/mol) | ρ_w | **K_gl** |
|---|---|---|---|---|---|
| 25 | 2.713 | 516 | 9.30 | 0.99705 | 0.376 |
| 32 (interp.) | — | — | — | — | ~0.43 |
| **90** | 3.141 | 1384 | 24.93 | 0.96534 | **0.855** |
| **120** | 3.239 | 1734 | 31.23 | 0.94306 | **1.013** |
| **145** | 3.287 | 1936 | 34.89 | 0.92141 | **1.089** |
| 160 | 3.304 | 2014 | 36.28 | 0.90744 | 1.110 |

Fraction of total sulfide in the gas, closed vial at equilibrium:
f = K_gl·V_g / (K_gl·V_g + V_l·(1 + 10^(pH − pK1))), pK1 from Suleimenov & Seward 1997
(`h2s_pk1_suleimenov_seward`: 6.489 at 90 °C, 6.454 at 120 °C, 6.486 at 145 °C).

| T °C | V_l / V_g (mL) | neutral (pH ≤ 5) | pH 5.5 | pH 7.0 |
|---|---|---|---|---|
| 90 | 5 / 10 | 0.63 | 0.61 | 0.29 |
| 120 | 5 / 10 | **0.67** | 0.65 | **0.31** |
| 145 | 5 / 10 | 0.69 | 0.66 | 0.34 |
| 120 | 1 / 9 | 0.90 | 0.89 | 0.67 |
| 120 | 1 / 19 | 0.95 | 0.94 | 0.81 |

The pH used is the in-situ pH at temperature; phosphate's own shift is not included.

---

## 6. WHAT THE PAPER DOES NOT STATE THAT A KINETIC MODEL NEEDS

1. **Mass-transfer rates** (only equilibrium): a sealed vial equilibrates its headspace on a minutes
   scale only if agitated; the paper says nothing about kinetics of exchange.
2. **Sugar, amino-acid or protein matrix effects** (pure water and NaCl only).
3. **pK1** (attempted, "extremely large uncertainties"; use Suleimenov & Seward 1997).
4. **Times to equilibrium** per run.

## 7. FLAGS

1. **[!] Two unit systems:** correlation in bar/molality, Table 2 / Appendix A / Fig. 3 in
   bar/mole fraction (factor 55.5).
2. **[!]** Appendix A row 40.4 °C, x = 0.000102 (print/OCR error in x).
3. Values above ~200 °C deviate from Drummond 1981 / Kozintseva 1964 (authors attribute to H2 / N2
   contamination in those studies) — irrelevant below 200 °C, where datasets agree.
4. Registry: paper absent from `data/keys/papers.yml`.

## Use in the model

- **The engine carries no H2S headspace.** The H2S species is the total aqueous sulfide pool
  (`species_sulfur.py`); the only headspace species is `OXR` (oxygen reservoir, B11); the only H2S
  loss is the fitted lumped `r_h2s_loss` ("volatilisation and every non-thiol sink",
  `sulfur.py`). Any sealed-vessel bundle with a recorded headspace could use K_gl above to split
  sulfide into gas and liquid; at 120-145 °C and pH ≤ 6, about half of the neutral H2S sits in a
  headspace equal in volume to the liquid (K_gl ≈ 1), unavailable to liquid-phase carbonyls, so the
  engine overstates liquid-phase H2S in such vessels by up to ~2× (more for larger headspaces).
- The headspace is a reversible reservoir in a sealed vessel, not a sink: it lowers the liquid
  concentration at temperature and returns on cooling.
- For **open / swept** vessels (Mizutani, Obata), K_gl sets the stripping driving force; the flow is
  still unknown.
</content>
</invoke>
