# Zheng & Ho 1994 — COMPLETE KINETIC TRANSCRIPTION

**Full extraction of every kinetic number in `data/articles/zheng1994.pdf`.**
Supersedes the single-value transcription (Ea 123.0 kJ/mol at pH 9) that was previously
the only thing taken from this paper.

---

## 0. PAPER IDENTITY — MATCHES THE DOSSIER DECLARATION

The PDF is the expected paper. No mis-file.

| field | value |
|---|---|
| Authors | **Yan Zheng and Chi-Tang Ho** |
| Title | **"Kinetics of the Release of Hydrogen Sulfide from Cysteine and Glutathione During Thermal Treatment"** |
| Venue | *Sulfur Compounds in Foods*, Mussinan, C.; Keelan, M., Eds.; **ACS Symposium Series 564**, Chapter 11 |
| Publisher / date | American Chemical Society, Washington DC; **publication date July 29, 1994**; received March 23, 1994 |
| Pages | **138–146** (PDF pages 1–9) |
| DOI | **10.1021/bk-1994-0564.ch011** |
| Affiliation | Dept. of Food Science, Cook College, New Jersey Agricultural Experiment Station, Rutgers, New Brunswick NJ 08903 |
| Acknowledgement | NJ Agric. Exp. Stn. Publication No. D-10205-11-93 |
| PDF character | ABBYY FineReader OCR over a clean scan; text layer usable, **superscripts detached by OCR** — every value below was re-read off the 200-dpi page raster, not trusted to the text layer |

**One caveat on the declaration itself:** the repo's FIT declaration
(`docs/reference/FIT_HOLDOUT_DECLARATION.md:79` and
`data/lit/extraction_dossiers/k3_final_parameter_inventory.md:1285`) says
*"Zheng 1994 Tables I/III/V (36 k + 8 Ea)"*. **The true count is 32 rate constants
(16 in Table I + 16 in Table III) and 8 activation energies (Table V).** There is no
36th rate constant. What the declaration probably meant to also cover is **Tables II
and IV — 8 linear k-vs-log[OH⁻] regressions = 16 further fitted coefficients**, which
are transcribed here in §4. The "36" is a miscount and should be corrected to 32 (+16
regression coefficients if those are wanted).

**Provenance codes (same convention as `k1_kinetic_parameters.md`):**
- **[M]** measured / reported directly by the authors in a table
- **[F]** fitted by the authors (regression parameter they printed)
- **[Z]** derived by me from the paper's own numbers — never printed by the authors
- **±** none. **The paper reports no error bars, no SE, no CI, no replicate SD anywhere.**
  Duplicate tubes were withdrawn at each time point but replicate scatter is never shown.

---

## 1. THE ONE-PARAGRAPH ANSWER

All **32 first-order rate constants (Tables I and III)** and all **8 activation energies
(Table V)** are legible, complete, and now transcribed — **not one cell was unreadable.**
All 32 rate constants pass an internal consistency check against the paper's own printed
half-lives (t½ = ln2/k reproduces every printed T½ to **< 0.25 %**), so the transcription
is verified twice over. **All four cysteine activation energies reproduce from the paper's
own Table I data to ±0.1 kcal/mol with matching R², so the cysteine block is fully
self-consistent and directly usable.** Two of the four glutathione Ea values
(**pH 3.0 and pH 5.0**) **cannot be reproduced** from Table III by any 4-point Arrhenius
fit and are flagged. Two printing errors were found and are documented in §6. The single
largest usability caveat is **not** legibility but extent: **every rate constant is an
initial-rate constant fitted over ≤ 13 % conversion** — the steepest measured trace,
Figure 4 (110 °C, pH 9), reaches only ln(Ct/C₀) ≈ −0.13 at its last point (180 min) — so
nothing here constrains long-time, high-conversion or high-temperature behaviour.

---

## 2. SYSTEM COMPOSITION — applies to every number in this dossier

Stated on p. 140 ("Reactions" and "Determination of hydrogen sulfide").

| variable | value as printed |
|---|---|
| Substrates | **0.1 M cysteine** *or* **0.1 M glutathione** (separate experiments; never mixed) |
| Chemical source | Sigma Chemical Co. (all chemicals except as indicated) |
| Solvent | Milli-Q deionized water from distilled water, **degassed by sonication 15 min before use** |
| Buffer, pH 3.0 & 5.0 | **0.1 M citrate–phosphate** |
| Buffer, pH 7.0 & 9.0 | **0.1 M phosphate–sodium hydroxide** |
| pH values | **3.0, 5.0, 7.0, 9.0** (nominal, set at room temperature; **no hot-pH correction is mentioned**) |
| Temperatures | **80 °C, 90 °C, 100 °C, 110 °C** (glycerine bath) |
| Volume / vessel | **10 mL** in a Kimax glass test tube with PTFE-coated liner screw cap |
| Atmosphere | **headspace flushed with nitrogen to expel oxygen**, then capped and tightened → **anaerobic, closed, no headspace sweep** |
| Sampling | **two tubes withdrawn every 10 min for the first 30 min, then every 30 min thereafter**; quenched immediately in ice-water |
| Longest time point | **180 min** (from the Figure 1–4 x-axes) |
| Analysis | **Orion sulfide/silver ion-selective electrode** + Orion digital pH/mV meter, Nernstian E = E₀ + S·log(A) |
| Analysis matrix | sample diluted into **SAOB** = sulfur antioxidation buffer = **2 N NaOH + 1 N EDTA + 4 % ascorbic acid** |
| Calibration | standard curve from standardized Na₂S in SAOB |
| Sugar / carbonyl co-reactant | **NONE.** This is neat amino-acid/peptide thermolysis in buffer — **no Maillard partner present** |
| Replication | duplicate tubes per time point; **no error statistics reported** |

**What "first order" means here.** Figures 1–4 plot **ln(Ct/C₀) vs time in minutes**, and
the ordinate is **negative and falling**. So K is a **first-order decay constant of the
precursor**, back-computed from accumulated H₂S (i.e. Ct ≡ C₀ − [H₂S]t), not a rate
constant on H₂S accumulation directly. The paper never states this conversion explicitly —
it is the only reading consistent with a logarithm being defined at t = 0. **Rate law:
d[cys]/dt = −K·[cys], K in min⁻¹, pseudo-first-order in the amino compound at fixed pH.**

---

## 3. TABLE I — CYSTEINE, 16 first-order rate constants **[M]/[F]**

**Anchor: Table I, p. 143 (PDF page 6).** Title as printed: *"Table I. First-order Rate
Constants (K), the Corresponding Regression Coefficients (R²) and the Half-Life (T½) of the
Release of Hydrogen Sulfide from Cysteine during Thermal Processing at Different pH Value"*.
Column headers as printed: `Temp | pH | K (min⁻¹) | T½ (hrs) | R²`.
**Units quoted exactly as printed: K in `min⁻¹`, T½ in `hrs`.**
Conditions: 0.1 M cysteine, buffers and atmosphere per §2. All cells legible.

| Temp | pH | **K (min⁻¹)** *as printed* | **T½ (hrs)** *as printed* | **R²** *as printed* | K in s⁻¹ **[Z]** | t½ recomputed, hrs **[Z]** | dev. |
|---|---|---|---|---|---|---|---|
| 80 °C | pH 3.0 | **1.87 × 10⁻⁶** | 6171.2 | 0.999 | 3.117 × 10⁻⁸ | 6177.8 | +0.11 % |
| 80 °C | pH 5.0 | **2.02 × 10⁻⁶** | 5713.4 | 0.998 | 3.367 × 10⁻⁸ | 5719.0 | +0.10 % |
| 80 °C | pH 7.0 | **9.71 × 10⁻⁶** | 1188.5 | 0.998 | 1.618 × 10⁻⁷ | 1189.7 | +0.11 % |
| 80 °C | pH 9.0 | **2.95 × 10⁻⁵** | 391.3 | 0.991 | 4.917 × 10⁻⁷ | 391.6 | +0.08 % |
| 90 °C | pH 3.0 | **8.30 × 10⁻⁶** | 1395.2 | 0.995 | 1.383 × 10⁻⁷ | 1391.9 | −0.24 % |
| 90 °C | pH 5.0 | **1.24 × 10⁻⁵** | 930.1 | 0.997 | 2.067 × 10⁻⁷ | 931.6 | +0.17 % |
| 90 °C | pH 7.0 | **5.99 × 10⁻⁵** | 192.9 | 0.996 | 9.983 × 10⁻⁷ | 192.9 | −0.02 % |
| 90 °C | pH 9.0 | **1.33 × 10⁻⁴** | 86.7 | 0.985 | 2.217 × 10⁻⁶ | 86.9 | +0.19 % |
| 100 °C | pH 3.0 | **4.15 × 10⁻⁵** | 278.2 | 0.998 | 6.917 × 10⁻⁷ | 278.4 | +0.06 % |
| 100 °C | pH 5.0 | **2.65 × 10⁻⁵** | 435.6 | 0.998 | 4.417 × 10⁻⁷ | 435.9 | +0.08 % |
| 100 °C | pH 7.0 | **2.30 × 10⁻⁴** | 50.3 | 0.995 | 3.833 × 10⁻⁶ | 50.2 | −0.14 % |
| 100 °C | pH 9.0 | **4.53 × 10⁻⁴** | 25.5 | 0.996 | 7.550 × 10⁻⁶ | 25.5 | +0.01 % |
| 110 °C | pH 3.0 | **5.23 × 10⁻⁵** | 220.8 | 0.973 | 8.717 × 10⁻⁷ | 220.9 | +0.04 % |
| 110 °C | pH 5.0 | **7.97 × 10⁻⁵** | 144.8 | 0.985 | 1.328 × 10⁻⁶ | 144.9 | +0.10 % |
| 110 °C | pH 7.0 | **3.35 × 10⁻⁴** | 34.5 | 0.955 | 5.583 × 10⁻⁶ | 34.5 | −0.04 % |
| 110 °C | pH 9.0 | **7.45 × 10⁻⁴** | 15.5 | 0.995 | 1.242 × 10⁻⁵ | 15.5 | +0.04 % |

**Verification [Z]: t½ = ln2/K reproduces all 16 printed half-lives to < 0.25 %.** The
transcription is therefore self-checked; no cell is ambiguous.

**Real anomaly, not an OCR artefact — at 100 °C, pH 3 is FASTER than pH 5.**
K(100 °C, pH 3.0) = **4.15 × 10⁻⁵** > K(100 °C, pH 5.0) = **2.65 × 10⁻⁵**, breaking the
monotone pH ordering that holds at 80 °C, 90 °C and 110 °C. Confirmed against the printed
half-lives (278.2 h vs 435.6 h, i.e. the inversion is in the half-life column too) *and*
against Figure 3 (p. 142), where the filled-square pH 3.0 trace at 100 °C is visibly steeper
than the "+" pH 5.0 trace. The authors flag this themselves
(p. 140): *"The nonconformity of reactions at pH 3.0 may mean that a different reaction
mechanism is involved."* **Treat pH 3.0 as a separate mechanistic channel, not as the low
end of the base-catalysed one.**

---

## 4. TABLE III — GLUTATHIONE, 16 first-order rate constants **[M]/[F]**

**Anchor: Table III, p. 144 (PDF page 7).** Printed as "Table III." (OCR renders it
"Table ID"). Title: *"First-order Rate Constants (K), the Corresponding Regression
Coefficients (R²) and the Half-Life (T½) of the Release of Hydrogen Sulfide from
Glutathione during Thermal Treatment at Different pH Values"*.
Same headers and same printed units: `K (min⁻¹)`, `T½ (hrs)`.
Conditions: **0.1 M glutathione**, otherwise identical to §2. All cells legible.

| Temp | pH | **K (min⁻¹)** *as printed* | **T½ (hrs)** *as printed* | **R²** *as printed* | K in s⁻¹ **[Z]** | t½ recomputed, hrs **[Z]** | dev. |
|---|---|---|---|---|---|---|---|
| 80 °C | pH 3.0 | **4.04 × 10⁻⁶** | 2859.5 | 0.996 | 6.733 × 10⁻⁸ | 2859.5 | 0.00 % |
| 80 °C | pH 5.0 | **4.11 × 10⁻⁶** | 2810.8 | 0.996 | 6.850 × 10⁻⁸ | 2810.8 | 0.00 % |
| 80 °C | pH 7.0 | **3.25 × 10⁻⁵** | 355.5 | 0.996 | 5.417 × 10⁻⁷ | 355.5 | −0.01 % |
| 80 °C | pH 9.0 | **1.34 × 10⁻⁴** | 86.2 | 0.989 | 2.233 × 10⁻⁶ | 86.2 | +0.01 % |
| 90 °C | pH 3.0 | **1.17 × 10⁻⁵** | 987.4 | 0.987 | 1.950 × 10⁻⁷ | 987.4 | 0.00 % |
| 90 °C | pH 5.0 | **1.57 × 10⁻⁵** | 735.8 | 0.992 | 2.617 × 10⁻⁷ | 735.8 | 0.00 % |
| 90 °C | pH 7.0 | **7.69 × 10⁻⁵** | 150.2 | 0.997 | 1.282 × 10⁻⁶ | 150.2 | +0.02 % |
| 90 °C | pH 9.0 | **2.44 × 10⁻⁴** | 47.3 | 0.989 | 4.067 × 10⁻⁶ | 47.3 | +0.10 % |
| 100 °C | pH 3.0 | **2.28 × 10⁻⁵** | 506.7 | 0.965 | 3.800 × 10⁻⁷ | 506.7 | 0.00 % |
| 100 °C | pH 5.0 | **4.94 × 10⁻⁵** | 233.9 | 0.991 | 8.233 × 10⁻⁷ | 233.9 | −0.02 % |
| 100 °C | pH 7.0 | **2.26 × 10⁻⁴** | 51.1 | 0.993 | 3.767 × 10⁻⁶ | 51.1 | +0.03 % |
| 100 °C | pH 9.0 | **7.74 × 10⁻⁴** | 14.9 | 0.995 | 1.290 × 10⁻⁵ | 14.9 | +0.17 % |
| 110 °C | pH 3.0 | **4.57 × 10⁻⁵** | 252.8 | 0.979 | 7.617 × 10⁻⁷ | 252.8 | 0.00 % |
| 110 °C | pH 5.0 | **7.46 × 10⁻⁵** | 154.9 | 0.956 | 1.243 × 10⁻⁶ | 154.9 | −0.03 % |
| 110 °C | pH 7.0 | **3.89 × 10⁻⁴** | 29.7 | 0.992 | 6.483 × 10⁻⁶ | 29.7 | −0.01 % |
| 110 °C | pH 9.0 | **1.09 × 10⁻³** | 10.6 | 0.990 | 1.817 × 10⁻⁵ | 10.6 | −0.01 % |

**Verification [Z]: all 16 half-lives reproduce to < 0.20 %.**

**Glutathione vs cysteine, all 16 matched conditions — the ratio K(GSH)/K(Cys) [Z]:**

| | pH 3.0 | pH 5.0 | pH 7.0 | pH 9.0 |
|---|---|---|---|---|
| 80 °C | 2.16 | 2.03 | 3.35 | 4.54 |
| 90 °C | 1.41 | 1.27 | 1.28 | 1.83 |
| 100 °C | 0.55 | 1.86 | 0.98 | 1.71 |
| 110 °C | 0.87 | 0.94 | 1.16 | 1.46 |

Note this contradicts the authors' own blanket statement (p. 140, *"the rate constants
derived from first-order kinetics, however, were higher than those for cysteine under the
same reaction conditions"*) in **four of the sixteen** cells: 100 °C pH 3.0 (0.55),
100 °C pH 7.0 (0.98), 110 °C pH 3.0 (0.87) and 110 °C pH 5.0 (0.94). The generalisation
holds at 80–90 °C and at pH 9 at every temperature, but **not universally** — and the
advantage shrinks systematically as temperature rises (mean ratio 3.0 at 80 °C, 1.4 at
90 °C, 1.3 at 100 °C, 1.1 at 110 °C), which is the direct consequence of glutathione's
lower Ea. **[Z], flagged.**

---

## 5. TABLE V — 8 ACTIVATION ENERGIES **[F]**

**Anchor: Table V, p. 145 (PDF page 8).** Title: *"Activation Energy of Release of Hydrogen
Sulfide from Cysteine and Glutathione at Different pH Values"*. Column headers as printed:
`Ea (kcal/mol)` and `Linearity R²`, one pair per substrate.
Method (p. 144): Arrhenius `K = K₀ × Exp(−Ea/RT)`, with **K₀ = pre-exponential (absolute)
rate constant, Ea in kcal/mol, R = 1.987 cal/mol/°K, T in °K**. Fitted over the four
temperatures 80/90/100/110 °C. **Units quoted exactly as printed: `kcal/mol`.**
All 8 cells legible. **K₀ is never printed anywhere in the paper.**

| pH | Cysteine **Ea (kcal/mol)** | Cysteine Linearity R² | Glutathione **Ea (kcal/mol)** | Glutathione Linearity R² |
|---|---|---|---|---|
| pH 3.0 | **31.3** | 0.940 | **18.8** | 0.997 |
| pH 5.0 | **31.8** | 0.975 | **30.8** | 0.996 |
| pH 7.0 | **32.3** | 0.942 | **22.8** | 0.990 |
| pH 9.0 | **29.4** | 0.967 | **19.9** | 0.970 |

### 5a. SI conversion **[Z]** (× 4.184 kJ/kcal)

| pH | Cysteine Ea, kJ/mol | Glutathione Ea, kJ/mol |
|---|---|---|
| pH 3.0 | **130.96** | **78.66** |
| pH 5.0 | **133.05** | **128.87** |
| pH 7.0 | **135.14** | **95.40** |
| pH 9.0 | **123.01** | **83.26** |
| **mean pH 3–9** | **130.54** | **96.55** |

**This is where the repo's single transcribed value comes from: 29.4 kcal/mol × 4.184 =
123.01 kJ/mol = cysteine at pH 9.0.** Correct, but it is the *lowest* of the four cysteine
values and the *least* representative of neutral-to-acidic food pH.

### 5b. Independent Arrhenius refit from the paper's own Tables I and III **[Z]**

Ordinary least squares of ln K vs 1/T over all four temperatures, K in min⁻¹,
R = 1.987 × 10⁻³ kcal/mol/K. `A` converted to s⁻¹ (÷ 60) for direct repo use.

| substrate | pH | **Ea refit (kcal/mol)** | Ea printed | **Ea refit (kJ/mol)** | R² refit | R² printed | **A (s⁻¹)** | reconciled? |
|---|---|---|---|---|---|---|---|---|
| cysteine | 3.0 | 31.36 | 31.3 | 131.2 | 0.940 | 0.940 | **9.79 × 10¹¹** | **YES — exact** |
| cysteine | 5.0 | 31.78 | 31.8 | 133.0 | 0.975 | 0.975 | **1.93 × 10¹²** | **YES — exact** |
| cysteine | 7.0 | 32.37 | 32.3 | 135.5 | 0.942 | 0.942 | **2.36 × 10¹³** | **YES — exact** |
| cysteine | 9.0 | 29.47 | 29.4 | 123.3 | 0.967 | 0.967 | **1.04 × 10¹²** | **YES — exact** |
| glutathione | 3.0 | 21.40 | 18.8 | 89.6 | 0.991 | 0.997 | 1.30 × 10⁶ | **NO — see below** |
| glutathione | 5.0 | 26.59 | 30.8 | 111.2 | 0.964 | 0.996 | 2.31 × 10⁹ | **NO — unresolved** |
| glutathione | 7.0 | 22.96 | 22.8 | 96.0 | 0.989 | 0.990 | **8.88 × 10⁷** | **YES** |
| glutathione | 9.0 | 20.04 | 19.9 | 83.8 | 0.965 | 0.970 | **5.48 × 10⁶** | **YES** |

**All four cysteine Ea values and their R² reproduce exactly** — the cysteine block is
internally airtight, and the eight `A` values above are the **pre-exponentials the paper
never printed**, now recoverable.

**Glutathione pH 3.0** reconciles only if the **80 °C point is dropped**: a 3-point fit on
90/100/110 °C gives **Ea = 18.83 kcal/mol, R² = 0.9992**, matching the printed 18.8 / 0.997.
Likely the authors excluded it; not stated in the paper. **[Z], inferred.**

**Glutathione pH 5.0 = 30.8 kcal/mol cannot be reproduced by any subset.** 4-point fit
gives 26.59 (R² 0.964); the four 3-point leave-one-out fits give 21.63, 27.04, 25.43 and
**32.57** (dropping 110 °C, R² 0.9992) — none is 30.8 with R² 0.996. **This value is
unexplained and should not be used.** It is also the value that produces the paper's
headline claim that Ea peaks ~2 pH units above the isoelectric point (pI GSH = 2.83, so
the peak "should" be near pH 5) — **the discussion on p. 146 rests on the one number in
Table V that does not reconcile with Table III.**

### 5c. What the ABSTRACT says vs what Table V says — a real discrepancy

Abstract (p. 138): *"31.3, 31.8, **32.2** and 29.4 kcal/mol for cysteine and 18.8, 30.8,
22.8 and 19.9 kcal/mol for glutathione at pH 3.0, 5.0, 7.0 and 9.0"*.
**Table V (p. 145) prints 32.3, not 32.2, for cysteine pH 7.0.** My refit gives **32.37**,
so **Table V's 32.3 is the correct one and the abstract has the typo.**
Consequence for the repo: `data/lit/arrhenius_params.yml` (cysteine_thermolysis audit_flag)
quotes the abstract's 32.2 → 134.7 kJ/mol and a pH-mean of **130.4**. Using Table V's 32.3
the mean is **130.54 kJ/mol**. A 0.14 kJ/mol difference — cosmetic, but the *provenance*
should be re-pointed at Table V rather than the abstract.

---

## 6. TABLES II and IV — 8 pH regressions, 16 fitted coefficients **[F]**

Not named in the FIT declaration, but they are kinetic parameters and they are the paper's
only explicit pH law. Equations transcribed **exactly as printed**, including the
`log[OH⁻]` argument (note: **log of the hydroxide concentration, not pOH**, so the
independent variable runs −11, −9, −7, −5 across pH 3→9 at 25 °C ion product).

### 6a. Table II — cysteine. **Anchor: p. 143 (PDF page 6), below Table I.**
Title: *"Linear Relationship between the Rate Constants of the Release of Hydrogen Sulfide
from Cysteine and the pH"*. Headers: `Temp. | Relationship | R²`. K in min⁻¹.

| Temp. | Relationship *as printed* | R² *as printed* | reproduces from Table I? **[Z]** |
|---|---|---|---|
| 80 °C | **K = 6.87 × 10⁻⁶ × log[OH⁻] + 6.19 × 10⁻⁵** | 0.939 | **YES** (slope 6.87 × 10⁻⁶, int. 6.183 × 10⁻⁵, R² 0.944) |
| 90 °C | **K = 3.01 × 10⁻⁵ × log[OH⁻] + 2.78 × 10⁻⁴** | 0.985 | **YES** (3.015 × 10⁻⁵, 2.795 × 10⁻⁴) |
| 100 °C | **K = 1.07 × 10⁻⁴ × log[OH⁻] + 9.88 × 10⁻⁴** | 0.999 | **YES** (1.066 × 10⁻⁴, 9.829 × 10⁻⁴) |
| 110 °C | **K = 1.66 × 10⁻⁴ × log[OH⁻] + 1.55 × 10⁻⁴** | 0.982 | **slope YES; INTERCEPT EXPONENT IS A TYPO — see below** |

**PRINTING ERROR #1 (Table II, 110 °C).** The intercept is printed **1.55 × 10⁻⁴**. The
3-point least-squares fit (pH 5, 7, 9) of Table I at 110 °C gives slope 1.663 × 10⁻⁴ (matches
the printed 1.66 × 10⁻⁴ exactly) and **intercept 1.5509 × 10⁻³**. The mantissa 1.55 is right,
**the exponent is wrong by one decade.** The printed equation is also internally impossible:
at pH 9, 1.66 × 10⁻⁴ × (−5) + 1.55 × 10⁻⁴ = **−6.75 × 10⁻⁴**, a negative rate constant.
**Use 1.55 × 10⁻³.** [Z], high confidence.

### 6b. Table IV — glutathione. **Anchor: p. 145 (PDF page 8), above Table V.**
Title: *"Linear Relationship between the Rate Constants of the Release of Hydrogen Sulfide
from Glutathione and the pH"*.

| Temp. | Relationship *as printed* | R² *as printed* | reproduces from Table III? **[Z]** |
|---|---|---|---|
| 80 °C | **K = 3.25 × 10⁻⁵ × log[OH⁻] + 2.85 × 10⁻⁴** | 0.905 | **YES — exact** (3.247 × 10⁻⁵, 2.843 × 10⁻⁴, R² 0.9045) |
| 90 °C | **K = 5.71 × 10⁻⁵ × log[OH⁻] + 5.12 × 10⁻⁴** | 0.933 | **YES — exact** (5.708 × 10⁻⁵, 5.117 × 10⁻⁴) |
| 100 °C | **K = 1.81 × 10⁻⁴ × log[OH⁻] + 1.62 × 10⁻³** | 0.920 | **YES — exact** (1.812 × 10⁻⁴, 1.618 × 10⁻³) |
| 110 °C | **K = 2.54 × 10⁻⁴ × log[OH⁻] + 4.82 × 10⁻³** | 0.954 | **slope YES, R² YES, INTERCEPT NO — see below** |

**PRINTING ERROR #2 (Table IV, 110 °C).** Slope 2.54 × 10⁻⁴ and R² 0.954 both reproduce
*exactly* from the 3-point (pH 5/7/9) fit of Table III — so the fit is the right one — but
that same fit gives **intercept 2.29 × 10⁻³**, not the printed **4.82 × 10⁻³**. The printed
value over-predicts K(pH 9, 110 °C) by 3.3× (3.55 × 10⁻³ vs the measured 1.09 × 10⁻³).
**Use 2.29 × 10⁻³.** [Z], high confidence.

**Fitting window, stated by the authors:** these regressions use **pH 5.0–9.0 only**. p. 140:
*"A linear relationship … was observed except for those reactions at pH 3.0"*; p. 144: *"the
linear relationship … was also observed for the reactions of glutathione in the range from
pH 5.0 to pH 9.0."* My reproductions confirm this: 3-point (pH 5/7/9) fits match the printed
slopes and R² to 3 significant figures, 4-point fits do not. **Do not extrapolate these eight
equations to pH 3, and do not use them below pH 5.**

---

## 7. FIGURES 1–4 — the only raw kinetic data shown **[M]**

Cysteine only. **No glutathione figure exists in the paper.**
Ordinate on all four: **`Ln (Ct/Co)`**; abscissa: **`Time (minutes)`**, 0–210 plotted,
data to 180 min. Legends carry pH 3.0 / 5.0 / 7.0 / 9.0.

| Figure | page | Temp | printed y-axis range | max |ln(Ct/C₀)| at 180 min | **max conversion [Z]** |
|---|---|---|---|---|---|
| Figure 1 | 141 (PDF p. 4) | 80 °C | 0.001 to −0.006 | ≈ 0.0052 (pH 9.0) | **≈ 0.5 %** |
| Figure 2 | 141 (PDF p. 4) | 90 °C | 0.005 to −0.030 | ≈ 0.024 (pH 9.0) | **≈ 2.4 %** |
| Figure 3 | 142 (PDF p. 5) | 100 °C | 0.02 to −0.10 | ≈ 0.080 (pH 9.0) | **≈ 7.7 %** |
| Figure 4 | 142 (PDF p. 5) | 110 °C | 0.02 to −0.16 | ≈ 0.130 (pH 9.0) | **≈ 12.2 %** |

Legend ordering differs between figures (Figs 1, 2, 4 list pH 9.0 first; Fig 3 lists pH 3.0
first) — noted because it is an easy mis-read. Digitised slopes agree with Table I to within
the line width in every case checked (e.g. Fig 1 pH 9 → 2.9 × 10⁻⁵ vs table 2.95 × 10⁻⁵;
Fig 3 pH 7 → 2.28 × 10⁻⁴ vs table 2.30 × 10⁻⁴), so **the figures corroborate Table I but add
no new parameters.**

### 7a. Absolute H₂S yields recoverable from the rate constants **[Z] — INFERRED, use with care**
If (and only if) Ct ≡ C₀ − [H₂S]t with C₀ = 0.1 M, then [H₂S] at 180 min = 100 mM ×
(1 − e^(−K·180)). This is the reading required for the figures to be well-defined, but the
paper never states it. Under it:

| condition | cysteine [H₂S] @180 min, mM | glutathione [H₂S] @180 min, mM |
|---|---|---|
| 80 °C pH 3 / 5 / 7 / 9 | 0.03 / 0.04 / 0.17 / 0.53 | 0.07 / 0.07 / 0.58 / 2.38 |
| 90 °C pH 3 / 5 / 7 / 9 | 0.15 / 0.22 / 1.07 / 2.37 | 0.21 / 0.28 / 1.37 / 4.30 |
| 100 °C pH 3 / 5 / 7 / 9 | 0.74 / 0.48 / 4.06 / 7.83 | 0.41 / 0.89 / 3.99 / 13.01 |
| 110 °C pH 3 / 5 / 7 / 9 | 0.94 / 1.42 / 5.85 / 12.55 | 0.82 / 1.33 / 6.76 / 17.82 |

**This partially retires the repo's claim that Zheng 1994 carries "no absolute
concentrations"** (`k3_final_parameter_inventory.md:1285`) — absolute H₂S is recoverable,
but only through an unstated normalisation, so it is **[Z]-inferred, not [M]**, and should
not become a scored benchmark target on its own.

---

## 8. MECHANISM AS STATED BY THE AUTHORS — no numbers attached

- **β-elimination**, base-catalysed: hydroxide abstracts the **α-carbon proton** (p. 140,
  which the authors describe as *"the β position to the sulfhydryl group"*); **dehydroalanine**
  forms with loss of the sulfhydryl, which takes up a proton to give H₂S. Ref. 14
  (Tarbell & Harnish, *Chem. Rev.* **1951**, 49, 1–90) — **cited, not measured here.**
- **Cited, not measured:** β-elimination is favoured when the α-amino group of cysteine is
  **acetylated** or the carboxyl is **esterified** — ref. 15 (Schneider & Westley,
  *J. Biol. Chem.* **1969**, 244, 5735–5744). **No rate constants are carried over from
  either reference; both are qualitative citations.**
- **Cited, not measured:** glutathione evolves H₂S faster than cysteine — ref. 11 (Ohloff,
  Flament & Pickenhagen, *Food Rev. Intl.* **1985**, 1, 99–148); and H₂S release from chicken
  muscle rises with pH — ref. 10 (Mecchi, Pippen & Lineweaver, *J. Food Sci.* **1964**, 29,
  393–399). **No numbers imported from either.**
- **Structural rationalisation (p. 146), qualitative:** Ea peaks ≈ 2 pH units above the
  isoelectric point — **pI(cysteine) = 5.07, pI(glutathione) = 2.83** are the only other
  numeric constants in the paper. At pH 7 cysteine is mainly ⁺H₃N–CH(COO⁻)–CH₂SH and
  self-association is invoked to explain the Ea maximum. **This argument depends on the
  unreconcilable GSH pH 5.0 Ea (§5b) and should be treated as speculative.**

**There are no other kinetic tables in the paper.** Contents of the 9 PDF pages:
p. 138 abstract + intro, p. 139 intro, p. 140 Experimental + start of Results, p. 141
Figs 1–2, p. 142 Figs 3–4, p. 143 Tables I + II, p. 144 Table III + Arrhenius method,
p. 145 Tables IV + V, p. 146 discussion + acknowledgement + 15 references. **Complete.**

---

## 9. USABILITY CAVEATS THAT APPLY TO EVERY NUMBER ABOVE

1. **Initial-rate constants only.** Maximum observed conversion across the entire study is
   **≈ 12 %** (110 °C, pH 9, 180 min, Fig. 4). Nothing in this paper constrains behaviour
   beyond ~15 % consumption, and the first-order form is untested past that.
2. **Anaerobic, closed vessel.** N₂-flushed, sealed, no headspace sweep. H₂S accumulates in
   solution and headspace; **any reverse or secondary consumption is folded into K.**
3. **No Maillard co-reactant.** No sugar, no dicarbonyl, no lipid. These are **neat
   thermolysis** constants and are an **upper bound on the free-H₂S source term** in any real
   system, where H₂S is simultaneously consumed by carbonyls.
4. **Nominal room-temperature pH.** Phosphate buffers shift ~−0.3 to −0.5 pH units between
   25 °C and 110 °C; the paper applies no correction. **The effective hot pH at "pH 9" is
   probably nearer 8.4–8.6.**
5. **0.1 M substrate** — 10–100× above typical food free-cysteine levels. Any
   concentration-dependent bimolecular channel (e.g. cystine formation) is exaggerated here.
6. **No error bars anywhere.** R² is the only goodness figure. With 4-point Arrhenius fits,
   the Ea values have on the order of ±2–4 kcal/mol standing uncertainty that the paper
   never quotes (my refits recover the point estimates exactly, so the R² values are real,
   but R² = 0.940 on 4 points is a **weak** constraint).
7. **Electrode method.** Sulfide/silver ISE after dilution into 2 N NaOH SAOB. Measures
   **total sulfide (S²⁻ + HS⁻ + H₂S)**, not molecular H₂S, and not speciated.

---

## 10. VERDICT — what is now usable, what remains missing

### NOW USABLE (was missing before this transcription)

| block | count | status |
|---|---|---|
| **Table I — cysteine K (min⁻¹), 16 values** | 16 | **FULLY USABLE.** Every cell legible; all 16 verified against the paper's own half-lives to < 0.25 %. Converted to s⁻¹ in §3. |
| **Table III — glutathione K (min⁻¹), 16 values** | 16 | **FULLY USABLE.** Same verification, < 0.20 %. |
| **Table V — cysteine Ea, 4 values** | 4 | **FULLY USABLE and independently reproduced** from Table I to ±0.1 kcal/mol with matching R². 131.0 / 133.1 / **135.1** / 123.0 kJ/mol at pH 3/5/7/9. |
| **Table V — glutathione Ea, pH 7.0 and 9.0** | 2 | **USABLE** (22.8 → 95.4 and 19.9 → 83.3 kJ/mol; both reproduce). |
| **Table V — glutathione Ea, pH 3.0** | 1 | **USABLE WITH A CAVEAT** — 18.8 kcal/mol (78.7 kJ/mol) reproduces only if the 80 °C point is dropped, which the paper does not disclose. |
| **Pre-exponential factors A (s⁻¹), 8 values** | 8 | **NEW — [Z] derived, never printed by the authors.** Cysteine: 9.79 × 10¹¹ / 1.93 × 10¹² / 2.36 × 10¹³ / 1.04 × 10¹² s⁻¹ at pH 3/5/7/9. These are the A/Ea *pairs* the repo needs; see §10a. |
| **Tables II + IV — 8 pH regressions, 16 coefficients** | 16 | **USABLE for pH 5–9 only**, with the two intercept corrections in §6. |
| **Absolute H₂S yields, 32 conditions** | 32 | **[Z], inferred** from the ln(Ct/C₀) normalisation. Useful as a sanity bound; **not benchmark-grade.** |

**Score against the declaration: 32 of 32 rate constants and 7 of 8 activation energies are
now transcribed and usable.** (The declaration's "36" is a miscount — see §0.)

### STILL MISSING / NOT USABLE

1. **Glutathione Ea at pH 5.0 = 30.8 kcal/mol — REFUSE.** Irreconcilable with Table III by
   any 4-point or 3-point fit (best candidates 26.6 and 32.6). Use the refit **26.6 kcal/mol
   (111 kJ/mol, R² 0.964)** if a value is needed, and label it [Z]. **1 of 8 Ea unusable.**
2. **K₀ / pre-exponentials are never printed.** All eight A values in §5b are mine, not the
   paper's.
3. **No uncertainties of any kind.** Any σ attached to these constants in the repo is the
   repo's own assumption and must be labelled as such.
4. **No glutathione time-course figures** — Table III cannot be re-derived from raw data.
5. **Nothing above 110 °C.** Extrapolation into roasting/extrusion regimes (>140 °C) is
   unsupported by this paper, consistent with the existing `arrhenius_params.yml` note that
   a second dry/pyrolytic channel at ~190 kJ/mol should be a separate step.
6. **No H₂S consumption kinetics.** This paper gives the source term only. The sink is not
   here.

### 10a. Direct consequences for `data/lit/arrhenius_params.yml → cysteine_thermolysis`

The current entry carries **Ea = 130.4 kJ/mol with A = 1.0 × 10¹⁴ s⁻¹**, the A having been
re-fitted by the repo rather than taken from the source. Two corrections now available:

- **Provenance:** the mean 130.4 was computed from the **abstract's** 32.2 kcal/mol at pH 7.
  **Table V prints 32.3**, and my refit confirms 32.37. The pH 3–9 mean from Table V is
  **130.54 kJ/mol**. Re-point the citation at **Table V, p. 145**, not the abstract.
- **A is 2 orders of magnitude too high, and the pair is inconsistent.** The paper's own data
  give, per pH, the matched pairs
  (pH 3: Ea 131.2 kJ/mol, A 9.8 × 10¹¹ s⁻¹), (pH 5: 133.0, 1.9 × 10¹²),
  (pH 7: 135.5, **2.4 × 10¹³**), (pH 9: 123.3, 1.0 × 10¹²).
  **At pH 7 — the food-relevant one — the source-consistent pair is Ea ≈ 135.5 kJ/mol with
  A ≈ 2.4 × 10¹³ s⁻¹**, which returns k(100 °C) = 2.6 × 10⁻⁶ s⁻¹ against the measured
  3.8 × 10⁻⁶ s⁻¹ (within the fit's own residual — R² is only 0.942 on 4 points, so the
  Arrhenius line does not pass through every datum). The shipped pair (130.4 / 1.0 × 10¹⁴)
  gives k(100 °C) = **5.6 × 10⁻⁵ s⁻¹ — ≈ 15× faster than Table I's measured value at pH 7,
  and ≈ 7× faster than even the pH 9 value.** Recommend replacing the pH-averaged pair with a
  **matched (Ea, A) set at a single pH**; averaging Ea across pH while fitting A separately
  is what introduced the discrepancy. Best practice here: pin A directly to a measured K —
  e.g. A = K(100 °C, pH 7) · exp(Ea/RT) — rather than re-fitting it free.
- **A pH law is available and currently unused:** Table II, 100 °C,
  `K = 1.07 × 10⁻⁴ × log[OH⁻] + 9.88 × 10⁻⁴` min⁻¹, valid pH 5–9. This is a better
  representation than the current `pH: "3-9 (weak dependence)"` string — the dependence is
  **not** weak: K rises **17×** from pH 5 to pH 9 at 100 °C (14.6× at 80 °C, 10.7× at 90 °C,
  9.3× at 110 °C).

---

*Extraction performed 2026-08-28. Every value re-read off 200-dpi page rasters of PDF pages
4–8; the OCR text layer was used only for the prose. No cell in Tables I–V was illegible.*
