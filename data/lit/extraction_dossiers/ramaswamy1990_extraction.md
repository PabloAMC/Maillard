# Ramaswamy, Ghazala & van de Voort 1990 — EXTRACTION (first-order thermal loss of thiamine in water and in a glucose/glycine/ascorbate mixture, 110-150 °C)

**Source on disk:** `data/articles/ramaswamy1990.pdf` (the publisher's scan, 6 pages, with an OCR text layer).
Every number below was read from the page images and checked by eye on 2026-10-09, then cross-checked
against `pdftotext -layout`; the text layer agrees with the images on every value in Tables 1 and 2.
Page numbers are the journal's printed page numbers (PDF page = journal page − 124). Written for the
thiamine route of the sulfur lane (`src/kinetic_core/sulfur.py`, `r_thi_hmp` / `r_thi_mesh`), whose
activation energy is borrowed.

| field | value |
|---|---|
| Title | "Degradation Kinetics of Thiamine in Aqueous Systems at High Temperatures" |
| Authors | H. Ramaswamy, S. Ghazala, F. van de Voort (Macdonald College of McGill University, Ste Anne de Bellevue, PQ, Canada) |
| Venue | Canadian Institute of Food Science and Technology Journal 1990, 23(2/3), 125-130 |
| DOI | 10.1016/S0315-5463(90)70215-5 |

## 1. Methods (p. 126)

- **Thiamine form:** thiamine **hydrochloride** (Canlab Canada). Nowhere is the mononitrate used.
- **Two matrices**, both in double-distilled water (DDW):
  - **B1/DDW** ("thiamine in water"): 0.1 g thiamine per litre DDW.
  - **B1/MIX** ("thiamine in mixture"): 0.1 g thiamine, 1.0 g L-ascorbic acid, 10.0 g D-glucose and
    8.0 g glycine per litre DDW.
  - 0.1 g/L thiamine HCl = **0.297 mM** (derived here: 0.1 / 337.27 g/mol × 1000).
- **pH:** "at natural pH" (p. 126). **No pH value is printed** for either matrix, and there is no buffer.
- **Water activity:** not printed (dilute aqueous solutions).
- **Oxygen:** not controlled. Ampoules were 10-mL glass with 4-mL aliquots, tips sealed in an
  oxygen-natural gas flame, so an air headspace was present; capillaries 100 µL in sealed tubes.
- **Heating:** two techniques, (1) ampoule, (2) thin-walled capillary tube (90 mm long, 1.8 mm OD,
  0.15 mm wall), both in a circulating oil bath (± 0.2 °C), quenched in ice water. 110, 120, 130, 140,
  150 °C (p. 126, "Data gathering and analyses").
- **Assay:** HPLC (Waters; µ-Bondapak C18, 3.9 mm × 30 cm; methanol:water 25:75 with 20 % low-UV PIC B6
  hexane sulfonic acid; 1 mL/min; UV 254 nm; 15 µL injection; thiamine retention time 9.0 min). This
  measures the **intact parent thiamine peak**: the measured quantity is **total loss of thiamine by all
  routes**, not formation of any product. No product (HMP, MFT, H2S or otherwise) is measured.
- **Kinetic analysis:** retention (% of unheated) → linear regression of ln(retention) vs time; k = −slope
  (eq. 1), D = 2.303/k (eq. 2). Temperature dependence by Arrhenius (Ea = −slope × R on ln k vs 1/T,
  eq. 3) and by TDT (z = −1/slope on log D vs T, eq. 4). Reference temperature 121.1 °C.
- **Reaction order:** first order, "semi-logarithmic" plots linear (Figs 1, 2, 5, 6; p. 127, 129).
- **Replicates:** ampoules two per time point; capillaries two tubes per analysis, two analyses per
  time-temperature combination. Replicate variation "mostly less than 5%" (p. 128).

## 2. Findings that matter

### Table 1 (p. 127), "Kinetic parameters for thiamine degradation"; k in min⁻¹, D in min

| matrix | method | T (°C) | R² | k value (min⁻¹) | D value (min) |
|---|---|---|---|---|---|
| thiamine in water | ampoule | 110 | 0.9927 | 0.0015 | 1509 |
| | | 120 | 0.9803 | 0.0069 | 333.1 |
| | | 130 | 0.9654 | 0.0181 | 127.0 |
| | | 140 | 0.9944 | 0.0302 | 76.3 |
| | | 150 | 0.9955 | 0.0568 | 40.5 |
| thiamine in water | capillary | 110 | 0.9804 | 0.0020 | 1167 |
| | | 120 | 0.9807 | 0.0073 | 313.5 |
| | | 130 | 0.9799 | 0.0187 | 123.2 |
| | | 140 | 0.9978 | 0.0324 | 71.0 |
| | | 150 | 0.9840 | 0.0406 | 56.8 |
| thiamine in mixture | ampoule | 110 | 0.9958 | 0.0032 | 717.1 |
| | | 120 | 0.9951 | 0.0065 | 354.0 |
| | | 130 | 0.9942 | 0.0116 | 197.8 |
| | | 140 | 0.9960 | 0.0190 | 121.2 |
| | | 150 | 0.9910 | 0.0260 | 88.5 |
| thiamine in mixture | capillary | 110 | 0.9922 | 0.0034 | 682.7 |
| | | 120 | 0.9930 | 0.0069 | 334.4 |
| | | 130 | 0.9936 | 0.0121 | 190.0 |
| | | 140 | 0.9990 | 0.0194 | 119.0 |
| | | 150 | 0.9954 | 0.0302 | 75.8 |

No SE or CI is printed for the individual k values; R² is the only per-temperature statistic.
Check (derived here): 2.303/k reproduces the printed D to within rounding of k (e.g. 2.303/0.0015 = 1535
vs 1509 printed; 2.303/0.0069 = 333.8 vs 333.1).

### Table 2 (p. 128), "Arrhenius and thermal death time (TDT) parameters ... in the temperature range of 110-150°C"

Units as printed: Ea and its error in kJ/mole; z and its error in C°; D₀ in min; k₀ in min⁻¹ (both at
121.1 °C). Footnote: "ᵃR² Variance with the Arrhenius method; ᵇR² Variance with TDT method;
ᶜt₀.₉₅ = standard error." The column is headed t₀.₉₅ but the footnote calls it a standard error; which of
the two it actually is (an SE, or a 95 % t-interval half-width) is not stated beyond that footnote.

| matrix | method | Ea (kJ/mole) | R²ᵃ | t₀.₉₅ᶜ (kJ/mole) | z (C°) | R²ᵇ | t₀.₉₅ᶜ (C°) | D₀ (min) | k₀ (min⁻¹) |
|---|---|---|---|---|---|---|---|---|---|
| thiamine in water | ampoule | 118.0 | 0.96 | 13.7 | 26.4 | 0.95 | 3.4 | 394.3 | 0.0060 |
| thiamine in water | capillary | 102.6 | 0.94 | 15.2 | 30.6 | 0.92 | 5.1 | 349.5 | 0.0067 |
| thiamine in mixture | ampoule | 71.1 | 0.99 | 4.6 | 43.8 | 0.98 | 3.6 | 354.3 | 0.0066 |
| thiamine in mixture | capillary | 73.4 | 0.99 | 3.0 | 42.4 | 0.99 | 2.4 | 337.8 | 0.0069 |

Ea is printed in kJ/mole, so no unit conversion is needed. The abstract (p. 125) rounds the water range
to "103-118 kJ/mole".

Text statements (p. 128): the ampoule and capillary techniques differ by up to ± 7.2 % from the midpoint
in Ea and z for water, and a t-test found no significant difference (p > 0.05) between techniques in k,
D, Ea or z. For the mixture (p. 129) the techniques differ by ± 1.5 %. The authors note (p. 128) that the
Arrhenius plot of the water data (Fig. 3) is **"a bit convex upward"** and the TDT plot concave upward,
and that the large standard errors on the water Ea "were due to the linearization of the curve and not
due to experimental errors". One visible consequence (derived here): water/ampoule k₀ at 121.1 °C
(0.0060 min⁻¹) is below the measured k at 120 °C in Table 1 (0.0069 min⁻¹).

### Derived here: k at 100 °C and 140 °C from the paper's own Table 2 parameters

Arrhenius: k(T) = k₀ · exp[−(Ea/R)(1/T − 1/394.25 K)], R = 8.314 J mol⁻¹ K⁻¹, T in K.
TDT: D(T) = D₀ · 10^((121.1 − T)/z), k = 2.303/D. Both are the paper's parameters evaluated at a new T;
**100 °C is 10 °C below the measured range (110-150 °C), an extrapolation**; 140 °C is inside it.

| matrix | method | k(100 °C) Arrhenius | k(100 °C) TDT | k(140 °C) Arrhenius | k(140 °C) TDT |
|---|---|---|---|---|---|
| water | ampoule | 0.00078 min⁻¹ (0.047 h⁻¹) | 0.00093 min⁻¹ | 0.0311 min⁻¹ (1.87 h⁻¹) | 0.0304 min⁻¹ |
| water | capillary | 0.00114 min⁻¹ (0.069 h⁻¹) | 0.00135 min⁻¹ | 0.0281 min⁻¹ (1.68 h⁻¹) | 0.0273 min⁻¹ |
| mixture | ampoule | 0.00194 min⁻¹ (0.116 h⁻¹) | 0.00214 min⁻¹ | 0.0178 min⁻¹ (1.07 h⁻¹) | 0.0176 min⁻¹ |
| mixture | capillary | 0.00195 min⁻¹ (0.117 h⁻¹) | 0.00217 min⁻¹ | 0.0192 min⁻¹ (1.15 h⁻¹) | 0.0190 min⁻¹ |

Worked example (water/ampoule, 100 °C): 0.0060 × exp[−(118000/8.314)(1/373.15 − 1/394.25)] =
0.0060 × exp(−2.036) = 0.00078 min⁻¹. Same, TDT: D = 394.3 × 10^(21.1/26.4) = 2483.5 min; 2.303/2483.5 =
0.00093 min⁻¹. At 140 °C the derived values sit close to the measured Table 1 k (0.0302 and 0.0324
water; 0.0190 and 0.0194 mixture), as they should inside the range.

## 3. What it means for the model

- **Water, unbuffered, "natural pH", thiamine HCl: Ea 102.6-118.0 kJ/mol** (two techniques, SE-labelled
  errors 13.7-15.2 kJ/mol). This sits on the centre (~100) to upper half of the 55-145 kJ/mol band the
  sulfur lane searched (`thiol_assembly` route bounds). The live engine runs `k_thi_hmp` on the frozen
  lumped formation Ea, 64.08 kJ/mol (`kinetic_core_b9_fit_report.json`), rate anchored at 145 °C, so
  this water series is 39-54 kJ/mol above it.
- **With glucose + glycine + ascorbic acid present: Ea 71.1-73.4 kJ/mol** (errors 3.0-4.6), i.e. the
  composition of the pot moves the apparent Ea by ~35-45 kJ/mol, and **raises** k below ~121 °C (the
  derived k(100 °C) in the mixture is ~2× that in water) while **lowering** it above (at 140-150 °C the
  water k is higher). A Maillard pot is closer to B1/MIX than to B1/DDW, but the paper cannot say which
  component (ascorbate, glucose, glycine, or the pH shift they cause) does it.
- The quantity is **total thiamine loss** (HPLC parent peak). The model's thiamine sinks
  (`r_thi_hmp`, `r_thi_mesh`) are branches of that loss; their summed k can be capped by these values,
  and the branch Ea equals the total-loss Ea only if the branching ratio is temperature-independent,
  which this paper does not test.
- k₀ at 121.1 °C ≈ 0.006-0.007 min⁻¹ in both matrices; k at 140 °C (inside the range, measured) 0.019-0.032
  min⁻¹, i.e. half-lives of ~21-36 min (derived here: ln 2 / k).
- Measured crossover (Table 1): at 110 °C the mixture k is about twice the water k (0.0032-0.0034 vs
  0.0015-0.0020 min⁻¹), at 120 °C they are equal within rounding, at 140-150 °C water is faster.

## What it does not give

- No numeric pH (only "natural pH"), no buffer, no water activity, no oxygen control or headspace
  composition.
- No product measurement: nothing on HMP, MFT, H2S or any thiamine-derived volatile; no branching.
- No pH dependence, no thiamine mononitrate, no phosphorylated forms, no meat or tissue matrix.
- No data below 110 °C; 100 °C values above are extrapolations.
- No per-temperature SE/CI on k; the Table 2 error column's statistic is labelled ambiguously (t₀.₉₅ header,
  "standard error" footnote). No individual retention data tabulated (only Figs 1, 2, 5, 6).
- The authors themselves flag curvature in ln k vs 1/T for the water series, so a single Ea is a
  linearisation over 110-150 °C.
