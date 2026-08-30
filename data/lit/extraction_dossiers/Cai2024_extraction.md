# Cai, Pan, Zhang, Yuan, Lao & Wu 2024 (10.1016/j.foodres.2024.114591) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/Cai2024.pdf` (8 pp.). Read method: **both** — born-digital text layer
(`pdftotext -layout`) for all body text and both tables, plus 600-dpi page rasters of PDF pp. 3–5 for
Fig. 2 (the kinetic scheme), Fig. 3A–E (the only source of the concentration–time data), and a
verification render of Table 1.

---

## 0. PAPER IDENTITY — **MATCHES** THE EXPECTED IDENTITY

| field | as printed |
|---|---|
| authors | Yanpei Cai, Xin Pan, Donghao Zhang, Lin Yuan, Fei Lao, Jihong Wu\* |
| title | "The kinetic study of 2-acetyl-1-pyrroline accumulation in the model system: An insight into enhancing rice flavor through the Maillard reaction" |
| venue | *Food Research International* |
| volume/pages/year | 191 (2024), article 114591, 8 pp. (running head "Food Research International 191 (2024) 114591") |
| DOI | https://doi.org/10.1016/j.foodres.2024.114591 |
| received/accepted | "Received 12 December 2023; Received in revised form 27 May 2024; Accepted 1 June 2024; Available online 13 June 2024" |
| copyright | "0963-9969/© 2024 Published by Elsevier Ltd." |
| affiliation | College of Food Science and Nutritional Engineering, China Agricultural University; National Engineering Research Center for Fruit & Vegetable Processing; Key Laboratory of Fruit & Vegetable Processing, Ministry of Agricultural and Rural Affairs; Beijing Key Laboratory for Food Non-thermal Processing, Beijing 100083, PR China |
| corresponding | Jihong Wu, wjh7268@cau.edu.cn, No.17 Qinghua East Road, Haidian District, 100083, Beijing |
| funding | National Key Research and Development Program of China, No. 2021YFD2100101 |
| PDF character | born-digital (Creator: Elsevier; Producer: Acrobat Distiller 8.1.0). Clean embedded text layer; no OCR. |

**Supplementary material is NOT in the repo.** The paper cites Table S1 (per-timepoint means with
ANOVA letters for all five responses), Table S2 (fitting performance), Figure S1 (mass balance) and
Figure S2 (2,3,5-trimethylpyrazine in cooked rice). None of these is available at
`data/articles/`; every number they contain is therefore recorded below only where the running text
quotes it, and marked **unavailable** otherwise.

---

## 1. ONE-PARAGRAPH VERDICT

This is the **highest-value multiresponse-kinetics paper in Wave K4a** and the only source in the
corpus that gives a *complete, closed, ten-step* ODE network for 2-acetyl-1-pyrroline (2-AP) in a
glucose/proline system with **all ten rate constants printed at three temperatures (100, 120, 140 °C)
and all ten activation energies printed**. It supplies both the 2-AP **formation** constants (k₆ for
MGO + PRO, k₇ for GO + PRO) and the 2-AP **depletion** constant (k₈, 2-AP → P₁), with separate Ea for
each — exactly the formation/depletion pair the flavour lane has never had for a pyrroline. It also
gives the α-dicarbonyl (methylglyoxal, glyoxal) time courses that feed those steps. **What it does
not give:** (i) **no uncertainties whatsoever on any rate constant or any Ea** — Table 1 is 40 bare
point estimates with not one standard error, confidence interval or highest-posterior-density
interval; (ii) **no RRMSE values**, despite RRMSE being defined in eq. (11); (iii) **no tabulated
concentration–time data** — every measured concentration in the main paper is figure-only (Fig. 3),
so all 105 concentration points below are `[fig]`; (iv) **no pH control and no pH measurement** — a
20 mmol/L phosphate buffer nominally at pH 9.0 has essentially zero buffer capacity there
(phosphate pK_a2 = 7.2, pK_a3 = 12.3) and the authors never report a post-heating pH; (v) **no
internal standard for 2-AP** — quantification is by 5-point *external* standard curve on HS-SPME/GC-MS,
which is the weakest possible calibration for a volatile in a changing matrix; (vi) no come-up-time
characterisation for a 40 mL sealed reactor whose first sampling point is 2 min. Two further defects
are load-bearing and are documented in §7: the fitted k₆/k₇ **over-predict the authors' own fed-
dicarbonyl experiment (Table 2) by ≈100×**, and three of ten Ea values are **negative**, one of them
(k₅) attached to a rate constant that is non-monotonic in temperature by seven orders of magnitude.
The paper is a rich **prior/structure** source and an excellent **hold-out** source; it is not a
source of absolute constants that can be trusted at face value.

---

## 2. SYSTEM DEFINITION — verbatim

### 2.1 Glucose/proline model system (the kinetics experiment, §2.2, p. 2)

Verbatim: *"Model systems were prepared by the study of Hofmann (1998) and Li et al. (2019) with
minor modifications. Proline (0.2 mol/L) and glucose (0.1 mol/L) were dissolved in the phosphate
buffer (20 mmol/L, pH = 9.0) to obtain the model system, which was thermally reacted in a high
temperature parallel reactor from YanZheng Instrument Co., ltd. (Shanghai, China). Referring to the
common processing of rice cooking (Leelayuthsoontorn & Thipayarat, 2006; Xu et al., 2019; Zhu et al.,
2020), the model system (40 mL) was sealed and heated at 100, 120, and 140 ◦C for 2, 5, 10, 15, 20,
25, and 30 min, respectively, followed by cooling in an ice-bath immediately."*

| variable | value as printed | provenance |
|---|---|---|
| L-proline, initial | 0.2 mol/L | [M] |
| D-glucose, initial | 0.1 mol/L | [M] |
| buffer | phosphate, 20 mmol/L | [M] |
| pH | 9.0 (nominal, set before heating) | [M] |
| hot-pH corrected? | **no** — no pH correction, no post-heating pH, no pH-vs-time data anywhere in the paper | [Z] |
| temperatures | 100, 120, 140 °C | [M] |
| times | 2, 5, 10, 15, 20, 25, 30 min | [M] |
| vessel | "high temperature parallel reactor", YanZheng Instrument Co. Ltd. (Shanghai) | [M] |
| volume | 40 mL | [M] |
| headspace/atmosphere | "sealed"; no gas specified, no purge, no headspace volume given. Authors later reason from *"the glucose/proline model system was heated in the sealed condition, the lack of oxygen…"* (p. 5), i.e. they treat it as O₂-limited but never measured O₂ | [M] + [Z] |
| agitation | **not stated** | [Z] |
| quench | *"cooling in an ice-bath immediately"* | [M] |
| come-up time | **not reported**; first sampling point is 2 min in a 40 mL sealed vessel | [Z] |
| replication | *"All experiments were conducted in triplicate."* (§2.9) | [M] |
| error bars | *"Values are given as means ± standard deviation (n = 3)"* (Fig. 3 caption, Tables 2 and 3 footnotes) | [M] |
| statistics | *"At 95 % confidence intervals, one-way ANOVA analysis was performed by Duncan's test in SPSS Statistics 26"* | [M] |

> **Internal inconsistency (initial proline).** Methods state 0.2 mol/L proline = 200 mmol/L. The
> t = 0 point of Fig. 3B reads **≈232 mmol/L** (§5.2, `[fig]`), a +16 % discrepancy. The glucose t = 0
> point reads 100 mmol/L and matches methods exactly. Both readings cannot be right; the model was
> presumably initialised at whatever value the fitted curves start from (all three predicted proline
> curves also start at ≈232). Flag: **any re-fit must decide which initial proline to use, and the
> choice moves k₃, k₆, k₇ and k₁₀ proportionally.**

### 2.2 α-Dicarbonyl/proline model system (§2.2, p. 2)

Verbatim: *"The methylglyoxal/proline and glyoxal/proline model systems reacted at 100 ◦C for 20 min,
initially containing 20 mmol/L proline and methylglyoxal or glyoxal with concentrations of 0, 10, 20,
or 40 mmol/L. Other reaction conditions were the same as the glucose/proline model system."*

So: 20 mmol/L PRO; MGO or GO at 0 / 10 / 20 / 40 mmol/L; 100 °C; 20 min; same 20 mmol/L phosphate
pH 9.0 buffer, same 40 mL sealed vessel, same ice-bath quench, n = 3.

### 2.3 Cooked-rice experiment (§2.3, p. 2)

Verbatim: *"In brief, 10 g rice was washed twice with 50 g water within 3 min and then mixed with
water with a ratio of 1:1.3. Before soaked for 40 min at 25 ◦C in a 100-mL round aluminum box with a
cover, the mixture was added with methylglyoxal or glyoxal of 0, 0.05, and 50 μmol/kg. The levels of
extra addition were determined by the preliminary study, both considering the color change in cooked
rice and potential health hazards caused by excessive addition of methylglyoxal (data not shown).
Finally, the aluminum box was heated in an electric steamer (100 ± 3 ◦C) for 15 min, followed by
insulating for 20 min."*

Rice: *"The same batch of milled fragrant rice (Arawana, Shanghai, China) from Wuchang China was
obtained from a local supermarket in Beijing, China. A total of 300 g of rice was fully mixed before
rice cooking."* Rice-matrix pH is stated only as a range at the end of the Discussion:
*"the pH difference between the real rice matrix (6.0–7.0) and model system in the study"* [M].

### 2.4 Analytical methods

**Proline (§2.4).** UPLC–MS/MS. Acquity UPLC + Xevo TQ-S triple quadrupole (Waters). Acquity UPLC BEH
C18, 1.7 μm, 2.1 × 100 mm, 40 °C. Injection 1.0 μL. Eluent A = 0.10 % (w/v) formic acid in water,
eluent B = acetonitrile. Gradient: 0.0–0.5 min 5 % B; 0.5–1.5 min 5→95 % B; 1.5–3.0 min 95 % B;
3.0–3.1 min 95→5 % B; 3.1–5.0 min 5 % B. Flow 0.4 mL/min. ESI+, MRM. *"The confirmation ion and
quantitation ion were m/z = 116.03 and 69.96. Additionally, the collision energy of the proline was
10 V. After optimization, the capillary voltage was 3.5 kV and the cone voltage was set at 20 V."*
Calibration: peak area vs standard concentration (**external**). SPE clean-up: PEP-SPE Cleanert resin
60 mg/3 mL (Bonna-Agela), 0.22 µm syringe filter.

**Glucose (§2.5).** HPLC–ELSD. Shimadzu LC-20A; polyamino HILIC column 5.0 μm, 250 × 4.6 mm (Dikma);
10 μL injection. A = acetonitrile, B = water. Gradient: 0–8 min 25 % B; 8–10 min 25→30 % B; 10–20 min
30 % B; 20–22 min 30→25 % B; 22–25 min 25 % B. Flow 1.0 mL/min, 40 °C. ELSD: nitrogen 350 kPa, drift
tube 70 °C. **External** calibration (peak area vs standard concentration).

**Methylglyoxal and glyoxal (§2.6).** *o*-Phenylenediamine (OPD) derivatisation then HPLC–DAD.
*"samples (400 μL) were reacted with 600 μL of 1.5 % (w/v) ortho-phenylenediamine dissolved in 50 %
methanol, in the dark at 37 ◦C for 24 h."* Agilent 1260 HPLC + DAD; Athena C18-WP 5 μm, 4.6 × 250 mm
(CNW), 35 °C; 10 μL injection. A = 0.075 % acetic acid in water; B = 80 % methanol + 20 % A. Flow
0.7 mL/min. Gradient: 0–23 min 10→50 % B; 23–28 min 50 % B; 28–31 min 50→70 % B; 31–40 min 70 % B;
40–44 min 70→10 % B; 44–50 min 10 % B. Detection 315 nm. **External** linear calibration curves.

**2-Acetyl-1-pyrroline (§2.7).** HS-SPME/GC–MS. *"Samples in model systems (2 mL) or cooked rice (3 g)
were added into 20 mL headspace bottles with the PTFE-silicone septum containing 0.8 g NaCl. Samples
were incubated at 70 ◦C for 20 min with agitation to reach the volatile equilibrium in sample
matrices and headspace. Then a 50/30 µm divinylbenzene/carboxen™/polydimethylsiloxane SPME fiber was
exposed to the headspace of samples at 70 ◦C for 30 min. Lastly, the SPME fiber was introduced into
the injector at 250 ◦C for 5 min."* GC: Agilent 7890 + 5975C MSD; HP-INNOWAX 30 m × 0.25 mm ×
0.25 µm; oven 40 °C (2 min) → 230 °C at 8 °C/min. Identification: NIST 17 library **and** the 2-AP
authentic standard. Quantification: *"The quantification of 2-acetyl-1-pyrroline was achieved through
5-point external standard curves. The target ion of m/z 111 operated in the selected ion monitoring
mode (SIM) was used as the quantitative ion (Rückriemen et al., 2015). Linear regression was employed
to plot the peak area of 2-acetyl-1-pyrroline standard (Y) against the concentration of 2-acetyl-1-
pyrroline standard (X) and expressed using the correlation coefficient (R²)."*

> **INTERNAL STANDARD FOR 2-AP: NONE.** Quantification is by **5-point external standard curve** only.
> No isotope-labelled 2-AP, no surrogate, no internal standard of any kind is named anywhere in the
> paper. The R² of the external calibration is said to be reported but **no numeric R² for the
> calibration is printed**. The 2-AP standard is *"2-acetyl-1-pyrroline (10 % w/w in toluene)…
> Toronto Research Chemicals (Toronto, Canada)"*, and the authors themselves note (Introduction)
> that *"Its purity was reported to degrade and turn red in 10 min (Fang & Cadwallader, 2014)."*
> **Absolute 2-AP levels from this paper carry an unbounded SPME-partition + standard-stability
> systematic.** Treat as `absolute_concentration: uncertain`.

### 2.5 Chemicals (§2.1)

L-Proline and D-glucose, HPLC grade, Solarbio (Beijing). HPLC-grade acetonitrile and methanol,
Thermo-Fisher. *o*-Phenylenediamine, TCI (Tokyo). Methylglyoxal, glyoxal, formic acid, acetic acid:
Macklin (Shanghai). 2-Acetyl-1-pyrroline standard (10 % w/w in toluene): Toronto Research Chemicals.
Other reagents analytical grade, Sinopharm.
**No concentration/assay purity is given for the commercial methylglyoxal or glyoxal** (both are
normally sold as ~40 % aqueous solutions and are self-oligomerising) — a hazard for §6 Table 2.

### 2.6 Kinetic modelling and fitting (§2.8)

*"According to Yu et al. (2020) with modifications, with MATLAB R2023a (MathWorks, USA), solved
differential equations were obtained through the "ode45" function, which were fitted to the
experimental data through the "lsqnonlin" function."*

Fit statistics defined: R² eq. (10) and RRMSE eq. (11):

- (10) `R² = 1 − Σᵢ₌₁ⁿ(ŷᵢ − yᵢ)² / Σᵢ₌₁ⁿ(yᵢ − ȳ)²`
- (11) `RRMSE = (1/ȳ) × sqrt( Σᵢ₌₁ⁿ(ŷᵢ − yᵢ)² / n )`

Arrhenius reparametrisation, eq. (12):

- (12) `k = k_av · exp( (Ea/R) × (1/T_av − 1/T) )`
  *"where R refers to the universal gas constant, T refers to the temperature concerned, k_av refers to
  the reaction rate constant at 120 ◦C."*

So **T_av = 120 °C = 393.15 K** and the reference constant is the 120 °C column of Table 1. `lsqnonlin`
is an unweighted (ordinary least-squares) non-linear solver with **no uncertainty output requested**,
which is why Table 1 has no error bars.

---

## 3. THE MULTIRESPONSE SCHEME — Fig. 2 and eqs. (1)–(9), COMPLETE

**Anchor: Fig. 2, p. 3 (PDF page 3), right column.** Caption as printed: *"Fig. 2. Kinetic model of
the glucose/proline system based on the proposed mechanism for the 2-acetyl-1-pyrroline thermal
generation. GLU: glucose; PRO: proline; MGO: methylglyoxal; GO: glyoxal; 2-AP: 2-acetyl-1-pyrroline;
INT: unquantifiable intermediates; P: products."* Read from a 600-dpi raster of PDF p. 3.

Every arrow in Fig. 2, transcribed with its number, exactly as drawn:

| step no. drawn on the arrow | reaction as drawn in Fig. 2 | co-substrate / co-product labelled beside the arrow | rate constant symbol (Table 1) |
|---|---|---|---|
| 1 | GLU → MGO | — | k₁ |
| 2 | GLU → GO | — | k₂ |
| 3 | GLU → INT | "PRO" entering on a curved arrow | k₃ |
| 4 | INT → MGO | "PRO" leaving on a curved arrow | k₄ |
| 5 | INT → GO | "PRO" leaving on a curved arrow | k₅ |
| 6 | MGO → 2-AP | "PRO" entering on a curved arrow | k₆ |
| 7 | GO → 2-AP | "PRO" entering on a curved arrow | k₇ |
| 8 | 2-AP → P₁ | — | k₈ |
| 9 | GLU → P₂ | — | k₉ |
| 10 | PRO → P₃ | — | k₁₀ |

The ODE system, verbatim from p. 4 (PDF page 4), eqs. (1)–(9):

```
(1)  r_glucose        = −(k1 + k2 + k9) × [GLU] − k3 × [GLU] × [PRO]
(2)  r_proline        = −k3 × [GLU] × [PRO] − k6 × [MGO] × [PRO] − k7 × [GO] × [PRO]
                        − k10 × [PRO] + (k4 + k5) × [INT]
(3)  r_methylglyoxal  = k1 × [GLU] + k4 × [INT] − k6 × [MGO] × [PRO]
(4)  r_glyoxal        = k2 × [GLU] + k5 × [INT] − k7 × [GO] × [PRO]
(5)  r_2-acetyl-1-pyrroline = k6 × [MGO] × [PRO] + k7 × [GO] × [PRO] − k8 × [2 − AP]
(6)  r_intermediates  = k3 × [GLU] × [PRO] − (k4 + k5) × [INT]
(7)  r_product1       = k8 × [2 − AP]
(8)  r_product2       = k9 × [GLU]
(9)  r_product3       = k10 × [PRO]
```

*"where r represents the reaction rate of each compound; k₁ − k₁₀ represents the reaction rate
constants for different reaction steps, respectively; [GLU], [PRO], [MGO],[GO], [2 − AP], and [INT]
are concentrations (mmol/L) of glucose, proline, methylglyoxal, glyoxal, 2-acetyl-1-pyrroline, and
unquantifiable intermediates, respectively."*

The scheme and the ODEs are **mutually consistent** — every arrow in Fig. 2 appears with the right
sign and the right order in eqs. (1)–(9), and vice versa. There is no scheme/equation mismatch.

> ### ⚠️ UNITS DEFECT (prominent — do not silently fix)
> §2.8 states concentrations are in **mmol/L**, but Table 1 prints the bimolecular constants k₃, k₆,
> k₇ in **L/(mol·min)**. Those two statements are mutually inconsistent by a factor of 10³ in every
> bimolecular term. A `[Z]` numerical check settles which is right: with k₇ = 0.38 L mol⁻¹ min⁻¹,
> [GO] = 1.58 × 10⁻⁴ mol/L and [PRO] = 0.198 mol/L (Fig. 3 at 100 °C, 20 min), plus k₆ = 0.17 and
> [MGO] = 4.8 × 10⁻⁵ mol/L, the 2-AP formation flux is 1.35 × 10⁻⁵ mol L⁻¹ min⁻¹; divided by
> k₈ = 1.77 min⁻¹ this gives a pseudo-steady-state **[2-AP] = 7.6 μmol/L**, against the observed
> 9.95 μmol/L. **The ODEs were therefore integrated in mol/L, not mmol/L.** The "(mmol/L)" in §2.8 is
> a printing error. Any downstream use of k₃, k₆, k₇ must use **mol/L**. [Z]

---

## 4. TABLE 1 — COMPLETE TRANSCRIPTION (the headline table)

**Anchor: Table 1, journal p. 5 (PDF page 5), spanning both columns.** Title as printed:
*"Kinetic parameters at different temperatures according to the proposed kinetic model in the
glucose/proline model system."*
Column headers as printed: `Elementary reaction step 1` | `Rate constant` | `100℃` | `120℃` | `140℃` |
`Activation energy (kJ/mol)`.
Footnote 1 as printed: *"GLU: glucose; PRO: proline; MGO: methylglyoxal; GO: glyoxal; 2-AP:
2-acetyl-1-pyrroline; INT: intermediates; P: products."*
Verified against a 600-dpi raster of PDF p. 5 — every cell below is legible and matches the text
layer. **No cell is unreadable. No cell carries an error bar, SE, CI or R².**

| # | Elementary reaction step | Rate constant (units as printed) | 100 ℃ | 120 ℃ | 140 ℃ | Activation energy (kJ/mol) |
|---|---|---|---|---|---|---|
| 1 | GLU → MGO | k₁ (1/min) | 3.81 × 10⁻⁵ | 2.02 × 10⁻⁴ | 0.0013 | 112.87 |
| 2 | GLU → GO | k₂ (1/min) | 5.07 × 10⁻⁸ | 0.0023 | 0.0028 | 355.28 |
| 3 | GLU + PRO → INT | k₃ (L/(mol·min)) | 8.15 × 10⁻⁴ | 0.026 | 0.49 | 205.32 |
| 4 | INT → MGO + PRO | k₄ (1/min) | 9.00 × 10⁻⁸ | 4.29 × 10⁻⁵ | 5.78 × 10⁻⁴ | 282.64 |
| 5 | INT → GO + PRO | k₅ (1/min) | 1.165 | 5.17 × 10⁻⁸ | 1.83 × 10⁻⁴ | − 294.05 |
| 6 | MGO + PRO → 2-AP | k₆ (L/(mol·min)) | 0.17 | 0.18 | 1.50 | 68.61 |
| 7 | GO + PRO → 2-AP | k₇ (L/(mol·min)) | 0.38 | 3.50 | 8.10 | 98.70 |
| 8 | 2-AP → P₁ | k₈ (1/min) | 1.77 | 40.75 | 239.71 | 157.89 |
| 9 | GLU → P₂ | k₉ (1/min) | 0.010 | 0.042 | 5.12 × 10⁻⁵ | − 164.48 |
| 10 | PRO → P₃ | k₁₀ (1/min) | 0.0060 | 0.0058 | 4.16 × 10⁻⁴ | − 84.02 |

All 30 rate constants and all 10 Ea are `[F]` — parameters fitted by the authors with `lsqnonlin` and
printed. None is `[M]`.

### 4.1 `[Z]` two-point Arrhenius reproduction of every printed Ea

Derived by me from the paper's own 100 °C and 140 °C columns via
`Ea = R·ln(k₁₄₀/k₁₀₀) / (1/373.15 − 1/413.15)`, with (1/373.15 − 1/413.15) = 2.5953 × 10⁻⁴ K⁻¹ and
R = 8.314 J mol⁻¹ K⁻¹. This is a *consistency check on the printed Ea*, not a new measurement.
Discrepancies are expected because the authors fitted eq. (12) to all three temperatures.

| k | k(100)→k(140) | Ea, `[Z]` two-point (kJ/mol) | Ea printed `[F]` (kJ/mol) | Δ |
|---|---|---|---|---|
| k₁ | 3.81e-5 → 0.0013 | +113.1 | 112.87 | +0.2 % |
| k₂ | 5.07e-8 → 0.0028 | +349.8 | 355.28 | −1.5 % |
| k₃ | 8.15e-4 → 0.49 | +205.1 | 205.32 | −0.1 % |
| k₄ | 9.00e-8 → 5.78e-4 | +280.9 | 282.64 | −0.6 % |
| k₅ | 1.165 → 1.83e-4 | **−280.6** | **− 294.05** | −4.6 % |
| k₆ | 0.17 → 1.50 | +69.7 | 68.61 | +1.6 % |
| k₇ | 0.38 → 8.10 | +98.0 | 98.70 | −0.7 % |
| k₈ | 1.77 → 239.71 | +157.3 | 157.89 | −0.4 % |
| k₉ | 0.010 → 5.12e-5 | **−168.9** | **− 164.48** | +2.7 % |
| k₁₀ | 0.0060 → 4.16e-4 | **−85.5** | **− 84.02** | +1.8 % |

**All ten printed Ea are internally reproducible from the printed rate constants to within 5 %.**
Table 1 is arithmetically self-consistent. The problem is not arithmetic; it is that three of the ten
Ea are negative and one (k₅) is attached to a rate constant that falls seven orders of magnitude from
100 → 120 °C and then rises 3500× from 120 → 140 °C.

### 4.2 `[Z]` derived ratios that the model actually needs

| quantity | 100 ℃ | 120 ℃ | 140 ℃ | note |
|---|---|---|---|---|
| k₇/k₆ (GO vs MGO reactivity with proline) | 2.24 | 19.4 | 5.40 | glyoxal is the more reactive Strecker partner at **every** temperature |
| k₂/k₅ (GO from glucose vs GO from INT) | 4.35 × 10⁻⁸ | 4.45 × 10⁴ | 15.3 | the authors' "GO comes mostly from glucose at 120/140 °C" claim; at 100 °C the ordering **inverts** |
| k₁/k₄ (MGO from glucose vs from INT) | 423 | 4.71 | 2.25 | MGO always predominantly from direct glucose fragmentation, but the margin collapses with T |
| k₈/(k₆+k₇) — depletion vs formation capacity | 3.21 min·mol/L | 11.1 | 24.9 | 2-AP depletion out-scales formation ~8× faster than formation across the 40 °C span |
| k₈(140)/k₈(100) | — | — | 135× | the whole "maximum at 100 °C" result rests on this one ratio |

---

## 5. FIGURE 3 — CONCENTRATION–TIME DATA, DIGITISED

**Anchor: Fig. 3, journal p. 4 (PDF page 4), full page width.** Caption as printed: *"Fig. 3. Kinetic
model fit (lines) to the experimental data (symbols) of reactants and products in the glucose/proline
model system. Values are given as means ± standard deviation (n = 3). R² means coefficient of
determination for the fitting performance."*

Panels: (A) Glucose (mmol/L); (B) Proline (mmol/L); (C) Methylglyoxal (μmol/L); (D) Glyoxal (μmol/L);
(E) 2-Acetyl-1-pyrroline (μmol/L). All five share `Time (min)`, 0–30, and three series: 100 °C blue
squares, 120 °C red triangles, 140 °C black diamonds, each with an "observed" and a "predicted" entry.

**All values in §5.1–5.5 are `[fig]`** — digitised by me from 600-dpi crops of PDF page 4
(`pdftoppm -r 600 -f 4`, one crop per panel). **No numeric version of these data is printed anywhere
in the main paper**; Table S1 (which holds the means and ANOVA letters) is not in the repo. Estimated
digitisation uncertainty is given per panel. Error bars are drawn on only a minority of the symbols
and are not individually resolvable for most points; where a bar is clearly resolvable I give it.

### 5.1 Fig. 3A — Glucose (mmol/L). Digitisation uncertainty ±2 mmol/L.

| t (min) | 100 ℃ | 120 ℃ | 140 ℃ |
|---|---|---|---|
| 0 | 100 | 100 | 100 |
| 2 | 100 (visible SD bar ≈ ±10) | 80 | 47 |
| 5 | 97 (SD bar ≈ ±3) | 68 | 35 |
| 10 | 92 | 57 | 28 |
| 15 | 89 | 49 | 24 |
| 20 | 84 | 44 | 21 |
| 25 | 84 (SD bar ≈ ±2) | 42 | 20 |
| 30 | 80 | 39 | 19 |

Printed fit statistics on the panel: **R² = 0.80 (100 °C), R² = 0.80 (120 °C), R² = 0.62 (140 °C)** `[F]`.

### 5.2 Fig. 3B — Proline (mmol/L). Digitisation uncertainty ±3 mmol/L.

| t (min) | 100 ℃ | 120 ℃ | 140 ℃ |
|---|---|---|---|
| 0 | **232** (see §2.1 inconsistency) | 232 | 232 |
| 2 | 224 (SD ≈ ±10) | 215 | 202 |
| 5 | 212 | 211 (SD ≈ ±15) | 190 |
| 10 | 202 | 196 (SD ≈ ±10) | 181 |
| 15 | 207 (SD ≈ ±7) | 199 (SD ≈ ±13) | 178 (SD ≈ ±10) |
| 20 | 198 (SD ≈ ±8) | 198 (SD ≈ ±14) | 170 (SD ≈ ±5) |
| 25 | 205 (SD ≈ ±4) | 198 (SD ≈ ±10) | 181 (SD ≈ ±4) |
| 30 | 188 | 195 | 163 (SD ≈ ±9) |

Printed fit statistics: **R² = 0.48 (100 °C), R² = 0.32 (120 °C), R² = 0.42 (140 °C)** `[F]`.
Note the 100 °C and 120 °C proline series are **non-monotonic** (rise at 15 and 25 min); the fits are
correspondingly poor, and the authors concede this: *"parallel and consecutive reaction steps were
bypassed, possibly causing free energy barriers and relatively low fitting performance of proline and
glucose in Table S2."*

### 5.3 Fig. 3C — Methylglyoxal (μmol/L). Digitisation uncertainty ±5 μmol/L (100/120 °C), ±10 (140 °C).

| t (min) | 100 ℃ | 120 ℃ | 140 ℃ |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 2 | 10 | 46 (SD ≈ ±8) | 200 (SD ≈ ±38) |
| 5 | 19 | 96 | 273 |
| 10 | 26 | 125 (SD ≈ ±12) | 296 |
| 15 | 38 | 130 | 337 |
| 20 | 48 | 163 (SD ≈ ±13) | 306 |
| 25 | 56 (SD ≈ ±8) | 184 | 317 |
| 30 | 62 (SD ≈ ±8) | 207 (SD ≈ ±14) | 272 (SD ≈ ±30) |

Printed fit statistics: **R² = 0.98 (100 °C), R² = 0.93 (120 °C), R² = 0.97 (140 °C)** `[F]`.
Cross-check against the running text (§3.2): *"The concentrations of methylglyoxal and glyoxal,
ranging from 10-340 μmol/L and 30–320 μmol/L respectively"* — my digitised MGO range is
**10 → 337 μmol/L**, matching the printed 10–340 to within digitisation error. This validates the
Fig. 3C calibration.

### 5.4 Fig. 3D — Glyoxal (μmol/L). Digitisation uncertainty ±6 μmol/L.

| t (min) | 100 ℃ | 120 ℃ | 140 ℃ |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 2 | 52 (SD ≈ ±6) | 283 (SD ≈ ±37) | 257 |
| 5 | 85 | 231 | 123 |
| 10 | 125 (SD ≈ ±19) | 202 (SD ≈ ±13) | 74 |
| 15 | 163 (SD ≈ ±7) | 169 (SD ≈ ±23) | 60 |
| 20 | 158 | 130 (SD ≈ ±6) | 51 |
| 25 | 137 (SD ≈ ±8) | 100 (SD ≈ ±10) | 47 |
| 30 | 127 (SD ≈ ±20) | 86 (SD ≈ ±5) | 35 |

Printed fit statistics: **R² = 0.83 (100 °C), R² = 0.92 (120 °C), R² = 0.64 (140 °C)** `[F]`.

> **Discrepancy between text and figure.** The text states glyoxal ranges **"30–320 μmol/L"**. My
> digitised maximum is **283 μmol/L** (120 °C, 2 min), whose SD bar reaches ≈320. The digitised
> minimum is 35, against the stated 30. Either the text quotes mean+SD as the upper bound, or the
> Fig. 3D 120 °C/2 min point sits a little higher than I read it. **Recorded, not resolved.** The
> lower bound agrees.

### 5.5 Fig. 3E — 2-Acetyl-1-pyrroline (μmol/L). Digitisation uncertainty ±0.15 μmol/L.

| t (min) | 100 ℃ | 120 ℃ | 140 ℃ |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 2 | 0.7 (SD ≈ ±0.4) | 4.05 (SD ≈ ±0.7) | 2.45 (SD ≈ ±0.45) |
| 5 | 2.6 (SD ≈ ±0.6) | 6.2 (SD ≈ ±0.7) | 0.80 |
| 10 | 6.6 (SD ≈ ±0.6) | 5.1 (SD ≈ ±0.8) | 0.73 |
| 15 | 6.65 (SD ≈ ±1.05) | 4.1 (SD ≈ ±0.2) | 0.89 (SD ≈ ±0.55) |
| 20 | **9.95** (SD ≈ ±1.3) | 2.65 (SD ≈ ±1.2) | 0.62 |
| 25 | 8.4 (SD ≈ ±1.4) | 1.6 (SD ≈ ±0.4) | 0.46 |
| 30 | 7.4 (SD ≈ ±2.3) | 1.15 (SD ≈ ±0.2) | 0.38 |

Printed fit statistics: **R² = 0.93 (100 °C), R² = 0.82 (120 °C), R² = 0.54 (140 °C)** `[F]`.

**Anchor for the one 2-AP value the authors print in text (journal p. 4, §3.3):** *"The highest yield
of 2-acetyl-1-pyrroline was observed in 100 ◦C, up to 9.95 μmol/L after 20 min of reaction."*
→ 2-AP(100 °C, 20 min) = **9.95 μmol/L** `[M]`. My independent digitisation of that point returned
9.95, which calibrates the whole of Fig. 3E.

### 5.6 Fit-statistic summary (all `[F]`, all read off the panels of Fig. 3)

| response | R² 100 ℃ | R² 120 ℃ | R² 140 ℃ |
|---|---|---|---|
| glucose | 0.80 | 0.80 | 0.62 |
| proline | 0.48 | 0.32 | 0.42 |
| methylglyoxal | 0.98 | 0.93 | 0.97 |
| glyoxal | 0.83 | 0.92 | 0.64 |
| 2-AP | 0.93 | 0.82 | 0.54 |

**RRMSE: defined in eq. (11), never reported.** Table S2 (the "fitting performance" table the text
points to) is not in the repo — **unavailable**.

---

## 6. TABLE 2 — fed-α-dicarbonyl / proline model system

**Anchor: Table 2, journal p. 6 (PDF page 6), left column.** Title as printed: *"Quantification of
2-acetyl-1-pyrroline in the α-dicarbonyl/proline model system with different additions of
α-dicarbonyls."*
Column headers as printed: `α-Dicarbonyl addition (mmol/L) 1` | `2-Acetyl-1-pyrroline concentration
(μmol/L) 2`.
Footnote 1: *"Values are given as means ± standard deviation (n = 3). n.d. means not detected due to
the concentration of 2-acetyl-1-pyrroline being below the detection limit."*
Footnote 2: *"The superscripts of a, b, and c indicate significant difference at p < 0.05."*
Conditions (from §2.2): 100 °C, 20 min, 20 mmol/L proline, same buffer/vessel/quench.

| α-dicarbonyl | addition (mmol/L) | 2-AP (μmol/L) as printed | dicarbonyl : proline molar ratio `[Z]` |
|---|---|---|---|
| Methylglyoxal | 0 | n.d. | 0 |
| Methylglyoxal | 10 | 0.19 ± 0.012 ᵇ | 1 : 2 |
| Methylglyoxal | 20 | 0.12 ± 0.043 ᵇᶜ | 1 : 1 |
| Methylglyoxal | 40 | n.d. | 2 : 1 |
| Glyoxal | 0 | n.d. | 0 |
| Glyoxal | 10 | **0.43 ± 0.088 ᵃ** | 1 : 2 |
| Glyoxal | 20 | 0.079 ± 0.0079 ᶜ | 1 : 1 |
| Glyoxal | 40 | n.d. | 2 : 1 |

All eight cells legible; none unreadable. All values `[M]`.

Authors' reading of the table (journal p. 4, §3.3): *"The highest productivity was observed when the
molar ratio of glyoxal and proline was 1:2 at 100 ◦C for 20 min. The 2-acetyl-1-pyrroline
concentration was below the detection limit when 2 equivalents of methylglyoxal or glyoxal were
added, which could be explained by competitive reactions occurring increasingly with a higher molar
ratio of α-dicarbonyl and proline in the alkaline aqueous model system (Adams et al., 2004; Wang &
Ho, 2012)."*

> ### ⚠️ THE PAPER'S OWN MODEL FAILS ITS OWN TABLE 2 BY ≈100× — `[Z]`
> Table 2 was run at 100 °C for 20 min in the same buffer, so it is directly predictable from Table 1.
> Pseudo-steady state of eq. (5) gives `[2-AP] ≈ (k₆[MGO] + k₇[GO])·[PRO] / k₈` with k₈(100 °C) =
> 1.77 min⁻¹ (integration over 20 min at k₈ = 1.77 min⁻¹ is fully relaxed; τ = 0.57 min):
>
> | fed system | k·[dicarbonyl]·[PRO] (mol L⁻¹ min⁻¹) | predicted 2-AP (μmol/L) | **observed** (μmol/L) | over-prediction |
> |---|---|---|---|---|
> | GO 10 mmol/L + PRO 20 mmol/L | 0.38 × 0.010 × 0.020 = 7.6 × 10⁻⁵ | **43** | 0.43 | **100×** |
> | MGO 10 mmol/L + PRO 20 mmol/L | 0.17 × 0.010 × 0.020 = 3.4 × 10⁻⁵ | **19.2** | 0.19 | **101×** |
>
> The same calculation applied to the *in-situ* glucose/proline system at 100 °C/20 min
> ([MGO] = 48 μmol/L, [GO] = 158 μmol/L, [PRO] = 198 mmol/L) predicts 7.6 μmol/L against an observed
> 9.95 — i.e. **the model is self-consistent on the data it was fitted to and wrong by two orders of
> magnitude on the one independent experiment in the same paper.** The gap is the same (100× / 101×)
> for both dicarbonyls, which points at a systematic cause — most plausibly that the commercial
> methylglyoxal/glyoxal solutions are largely present as hydrates/oligomers so that the nominal
> 10 mmol/L is not 10 mmol/L of free α-dicarbonyl, and/or that 2-AP recovery collapses in the fed
> system. **This is the single most important finding in this dossier for the model: k₆ and k₇ must
> NOT be propagated as absolute bimolecular rate constants.** Their *ratio* (k₇/k₆) and their
> *temperature dependence* (Ea) survive; their magnitude does not.

---

## 7. TABLE 3 — cooked rice with added α-dicarbonyl

**Anchor: Table 3, journal p. 6 (PDF page 6), right column.** Title as printed: *"Effects of the
addition of methylglyoxal or glyoxal on the 2-acetyl-1-pyrroline concentration in the cooked rice."*
Column headers as printed: `α-Dicarbonyl addition (μmol/kg) 1` | `2-Acetyl-1-pyrroline
concentration(μmol/kg) 2`.
Footnote 1: *"Values are given as means ± standard deviation (n = 3)."*
Footnote 2: *"The superscripts of a, b, and c indicate significant difference at p < 0.05."*

| α-dicarbonyl | addition (μmol/kg) | 2-AP (μmol/kg) as printed | Δ vs control `[Z]` |
|---|---|---|---|
| Methylglyoxal | 0 | 6.75 ± 0.14 ᶜ | — |
| Methylglyoxal | 0.05 | **9.01 ± 0.54 ᵃ** | **+33.5 %** |
| Methylglyoxal | 50 | 7.92 ± 0.094 ᵇ | +17.3 % |
| Glyoxal | 0 | 6.75 ± 0.14 ᶜ | — (same control row, reprinted) |
| Glyoxal | 0.05 | 7.18 ± 0.89 ᵇᶜ | +6.4 % |
| Glyoxal | 50 | 7.79 ± 0.31 ᵇ | +15.4 % |

All six cells legible. All values `[M]`; the Δ column is `[Z]`.
Authors' text confirms the +33.5 %: *"especially with the 33.5 % concentration promoted after adding
0.05 μmol/kg methylglyoxal."* Note the **0 μmol/kg row is one and the same control** printed twice —
the table has 5 distinct samples, not 6.

> **Non-monotonic dose response, both compounds, in opposite directions.** MGO: 0.05 > 50 > 0.
> GO: 50 > 0.05 > 0. The authors' explanation for the MGO inversion is *"the acceleration of the
> extra methylglyoxal to the advanced stages of the Maillard reaction, including the 2-acetyl-1-
> pyrroline transition to brown polymers"*, supported by the observation (Figure S2, **unavailable**)
> that 50 μmol/kg MGO rice contained more 2,3,5-trimethylpyrazine. The GO ordering is not explained.
> A 1000× dose change moving the response by ±15 % in either direction is not a dose–response curve
> the model can be scored on; treat Table 3 as directional.

---

## 8. THE GLYOXAL FINDING — IN FULL, VERBATIM

This is the paper's headline mechanistic claim. Four separate strands, each recorded uncompressed.

### 8.1 The claim: glyoxal *formation* is rate-determining for the whole 2-AP network

**Abstract, p. 1, verbatim:** *"Using the multi-response kinetic modeling to derive kinetic
parameters, the formation of glyoxal, as the reactive intermediate, was rate-determining for the
overall generation rate of 2-acetyl-1-pyrroline."*

**§3.2, journal p. 5, verbatim (the full passage, uncut):** *"By estimating rate constants, glyoxal
might predominantly originate from glucose degradation (k2) compared with the cleavage of
Maillard-derived intermediates (k5) at higher reaction temperatures of 120 ◦C and 140 ◦C, consistent
with previous studies (Kocadağlı & Gökmen, 2016; Berk et al., 2021). Besides, concerning both the
lowest rate constants and highest Ea in the overall reaction of 2-acetyl-1-pyrroline generation,
glyoxal formation was identified as the rate-determining step in the kinetic scheme of
2-acetyl-1-pyrroline. However, showing the negative Ea, the pathway of glyoxal formation through the
degradation of intermediates was discrepant with the Arrhenius equation, possibly due to the
complexity of reaction pathways and the involvement of oxygen. Glyoxal was reported to be generated
via the Namiki pathway from the Schiff base and the cleavage of glucosone enolized from glucose, both
dependent on the presence of the oxidative mechanism (Gobert & Glomb, 2009). Since the glucose/proline
model system was heated in the sealed condition, the lack of oxygen contributed to complex changes of
k5 with temperature."*

**Conclusion, journal p. 6, verbatim:** *"During the Maillard reaction, the formation of the
intermediate glyoxal, through the degradation of glucose and proline was identified as the
rate-determining step for the overall generation of 2-acetyl-1-pyrroline."*

**The numbers behind the claim** (Table 1, all `[F]`):

| supporting number | value | where |
|---|---|---|
| k₂ (GLU → GO) at 100 °C | **5.07 × 10⁻⁸ min⁻¹** — the smallest number in the entire table | Table 1 row 2 |
| Ea(k₂) | **355.28 kJ/mol** — the largest Ea in the entire table | Table 1 row 2 |
| k₂ at 120 / 140 °C | 0.0023 / 0.0028 min⁻¹ | Table 1 row 2 |
| k₅ (INT → GO + PRO) | 1.165 / 5.17 × 10⁻⁸ / 1.83 × 10⁻⁴ min⁻¹ at 100/120/140 °C | Table 1 row 5 |
| Ea(k₅) | **− 294.05 kJ/mol** (negative) | Table 1 row 5 |
| k₂ > k₅ at 120 °C | 0.0023 vs 5.17 × 10⁻⁸ → **k₂/k₅ = 4.45 × 10⁴** `[Z]` | derived |
| k₂ > k₅ at 140 °C | 0.0028 vs 1.83 × 10⁻⁴ → **k₂/k₅ = 15.3** `[Z]` | derived |
| k₂ ≪ k₅ at 100 °C | 5.07 × 10⁻⁸ vs 1.165 → **k₂/k₅ = 4.35 × 10⁻⁸** `[Z]` — the ordering **inverts** at the one temperature where 2-AP is maximal | derived |

### 8.2 The second glyoxal strand: GO is the *more reactive* Strecker partner than MGO

**§3.3, journal p. 6, verbatim:** *"As shown in Table 1, with the higher rate constant, glyoxal (k7)
was more reactive to combine with proline to yield 2-acetyl-1-pyrroline than methylglyoxal (k6), which
also contributed to the earlier descent of glyoxal concentrations than methylglyoxal in Fig. 3C and
Fig. 3D."* And the mechanistic hypothesis, verbatim: *"Hence, glyoxal was proposed to experience a
similar reaction mechanism with proline for the final 2-acetyl-1-pyrroline generation, by previously
yielding 1-pyrroline or 2-methyl-1-pyrroline."*

Numbers: k₇ / k₆ = 0.38 / 0.17 (100 °C), 3.50 / 0.18 (120 °C), 8.10 / 1.50 (140 °C) `[F]` →
**k₇/k₆ = 2.24, 19.4, 5.40** `[Z]`. Ea(k₇) = 98.70 vs Ea(k₆) = 68.61 kJ/mol `[F]`.

### 8.3 The third strand: fed glyoxal beats fed methylglyoxal in the model system…

Table 2: at 10 mmol/L (1 : 2 dicarbonyl : proline), glyoxal gives **0.43 ± 0.088 μmol/L** 2-AP versus
methylglyoxal's **0.19 ± 0.012 μmol/L** — a **2.3× advantage** for glyoxal `[Z]`, in the same
direction as k₇/k₆ = 2.24 at 100 °C. That agreement (2.3 vs 2.24) is the one place in the paper where
Table 1 and Table 2 agree, *in ratio*, even though they disagree 100× in magnitude (§6).

### 8.4 …but glyoxal LOSES in real cooked rice, and the authors do not resolve it

Table 3: 0.05 μmol/kg glyoxal gives 7.18 ± 0.89 (letter "bc" — **not** significantly above the
6.75 ± 0.14 control at p < 0.05), while 0.05 μmol/kg methylglyoxal gives 9.01 ± 0.54 (letter "a").
Authors, verbatim: *"Lower 2-acetyl-1-pyrroline enhancement was shown in rice added with glyoxal than
methylglyoxal, indicating possibly changed kinetic parameters in its generation from glyoxal and
proline, considering the complex effect of physical properties in rice (Vyazovkin, 2016)."*
They also note the regulatory asymmetry, verbatim: *"Methylglyoxal has been allowed in China
(GB2760-2014) as a food flavoring substance… generally recognized as safe (GRAS) within maximum use
levels of 0.05–3 mg/kg in baked food by Flavor and Extract Manufacturers Associations (FEMA)
(Hall, 1965). However, glyoxal has not been listed as a food additive or flavoring by relative
standards."*

### 8.5 Objections to the glyoxal claim that the model must carry forward — `[Z]`

1. **"Lowest rate constants" is true only at 100 °C, and only for k₂.** At 120 °C, k₂ = 0.0023 is
   *larger* than k₁ (2.02 × 10⁻⁴), k₄ (4.29 × 10⁻⁵) and k₅ (5.17 × 10⁻⁸). At 140 °C, k₂ = 0.0028 is
   larger than k₁ (0.0013) and k₄ (5.78 × 10⁻⁴). The rate-determining-step claim is therefore
   supported by the 100 °C column only. The text's plural, *"the lowest rate constants"*, is not
   supported by Table 1 at two of the three temperatures.
2. **Magnitude comparison across mixed units is not meaningful.** k₁, k₂, k₄, k₅, k₈, k₉, k₁₀ are
   min⁻¹; k₃, k₆, k₇ are L mol⁻¹ min⁻¹. Comparing "lowest rate constant" across the two families is
   dimensionally invalid; only fluxes are comparable.
3. **The Fig. 3D data do not obviously support a glyoxal bottleneck.** Glyoxal *peaks earlier and
   higher than* methylglyoxal at 120 °C (283 vs 46 μmol/L at 2 min) and is comparable at 100 °C
   (163 vs 38 at 15 min). An intermediate whose formation is rate-determining should be the scarce
   one; here it is the abundant one at short times.
4. **The Conclusion contradicts the scheme.** The Conclusion says glyoxal forms *"through the
   degradation of glucose and proline"*, but Fig. 2 step 2 (k₂) is **GLU → GO with no proline**;
   proline enters glyoxal formation only via k₅ (INT → GO + PRO), whose Ea is negative and which the
   authors themselves discount. Text and scheme disagree; both are quoted above.
5. **k₅'s negative Ea is an artefact the authors half-acknowledge.** k₅ falls 2.3 × 10⁷-fold from
   100 → 120 °C and then rises 3540-fold from 120 → 140 °C. This is not an Arrhenius process; it is a
   parameter absorbing model error, in the exact same lump (INT) that they also admit is
   unidentifiable: *"Unquantifiable intermediates in Fig. 2 included the Schiff base, Amadori
   rearrangement product, deoxyglucosones, etc., which could not be measured quantitively due to
   experimental restrictions."* Their oxygen explanation is untested — no O₂ was measured, and the
   vessel is described only as "sealed".

---

## 9. OTHER NUMERIC CLAIMS IN THE RUNNING TEXT

| claim | value | anchor | provenance |
|---|---|---|---|
| average glucose depletion rate vs proline | glucose **1.7×** faster | §3.1, journal p. 4 | [M] |
| MGO concentration range across all runs | **10–340 μmol/L** | §3.2, journal p. 4 | [M] |
| GO concentration range across all runs | **30–320 μmol/L** | §3.2, journal p. 4 | [M] (conflicts with Fig. 3D, §5.4) |
| 2-AP maximum | **9.95 μmol/L**, 100 °C, 20 min | §3.3, journal p. 4 | [M] |
| Ea, glucose degradation (quoted in text as the k₁ route) | **112.87 kJ/mol** | §3.2, journal p. 4 | [F] |
| Ea, "Maillard-related pathways" (= k₄) | **282.64 kJ/mol** | §3.2, journal p. 4 | [F] |
| Ea, 2-AP generation, both routes | **68.61 and 98.70 kJ/mol** | §3.3, journal p. 4 | [F] |
| Ea, 2-AP depletion | **157.89 kJ/mol** | §3.3, journal p. 4 | [F] |
| mass-balance recovery at 30 min | **80.41 % (100 °C), 69.50 % (120 °C), 54.40 % (140 °C)** | §3.4, journal p. 4; underlying plot = Figure S1, **unavailable** | [M] |
| 2-AP odour threshold in starch | **0.0073 μg/kg** | Introduction, p. 1 | [C] (no primary source given at that sentence; context is Verma & Srivastav 2022) |
| MGO FEMA GRAS use level in baked food | **0.05–3 mg/kg** | §3.5, journal p. 6 | [C] (Hall & Oser 1965) |
| rice matrix pH | **6.0–7.0** | §3.5, journal p. 6 | [M/C] |
| global rice production 2023 | 523.5 million tons | Introduction, p. 1 | [C] (FAO 2023) |
| 2-AP standard purity degrades / turns red | in **10 min** at ambient | Introduction, p. 1 | [C] (Fang & Cadwallader 2014) |

Authors' own account of the mass-balance loss, verbatim: *"The decline of the mass balance can be
explained by the importance of melanoidins formation (Martins & Van Boekel, 2005), which was regarded
as advanced and final stages in the Maillard and caramelization reaction."* — a direct hand-off to the
melanoidin ε lane (see the `brands2002b` and `martins2003c` dossiers).

Authors' own account of the negative Ea for k₉ and k₁₀, verbatim: *"However, as shown in Table 1, the
Ea of the decomposition of glucose and proline was negative, which could be explained by the
simplified kinetic model. In order to minimize unknown parameters, parallel and consecutive reaction
steps were bypassed, possibly causing free energy barriers and relatively low fitting performance of
proline and glucose in Table S2. Similar results were also reported by Göncüoğlu Taş & Gökmen (2017)
and Berk et al. (2021). Lower kinetic rate constants (k9, k10) were observed at the higher
temperature, which might also imply that the self-degradation of glucose and proline was not as
crucial as other reaction pathways for 2-acetyl-1-pyrroline generation."*

Contradiction with prior literature that the authors flag themselves, verbatim: *"Inconsistent with
our findings, Chan & Reineccius found that 2-acetyl-1-pyrroline formation followed the zero-order
reaction kinetics (2005), and it could be explained by the dependence of the kinetics of
2-acetyl-1-pyrroline accumulation on different reaction conditions."*

---

## 10. DEFECT REGISTER (consolidated, for the orchestrator)

| # | defect | severity | consequence |
|---|---|---|---|
| D1 | **No uncertainty on any of the 30 k or 10 Ea.** `lsqnonlin` point estimates only. | high | any prior built from Table 1 must carry an *assumed* σ; nothing in the paper bounds it |
| D2 | **Fitted k₆/k₇ over-predict the paper's own Table 2 by ≈100× (§6).** | **critical** | magnitudes of k₆, k₇ are not transferable; only ratio and Ea are |
| D3 | Three negative Ea (k₅ −294.05, k₉ −164.48, k₁₀ −84.02); k₅ non-monotonic over 7 orders | high | k₅, k₉, k₁₀ are error sinks, not barriers. Do not ingest as barriers |
| D4 | pH 9.0 set in 20 mmol/L phosphate — **essentially no buffer capacity at pH 9** (pK_a2 7.2, pK_a3 12.3); no hot-pH correction; no post-heat pH | high | the true reaction pH is unknown and drifting; every constant carries an unquantified pH systematic |
| D5 | 2-AP quantified by **external** 5-point curve, HS-SPME, no internal standard, standard known to degrade in 10 min | high | absolute 2-AP levels uncertain; flag `absolute_concentration: uncertain` |
| D6 | Initial proline: methods 200 mmol/L vs Fig. 3B t=0 ≈232 mmol/L | medium | ±16 % on every bimolecular constant |
| D7 | Units: §2.8 says mmol/L, Table 1 bimolecular constants are L/(mol·min); `[Z]` check shows mol/L was used | medium | resolved above, but must be recorded |
| D8 | Text says GO range 30–320 μmol/L; Fig. 3D max digitises to 283 | low | recorded, unresolved |
| D9 | RRMSE defined, never reported; Table S1, S2, Fig. S1, S2 unavailable | medium | all concentration data are `[fig]`; no per-point means/SD available |
| D10 | No come-up time for a 40 mL sealed vessel with a 2-min first point at 140 °C | medium | short-time points at 140 °C are not isothermal |
| D11 | Conclusion ("glyoxal from degradation of glucose **and proline**") contradicts Fig. 2 step 2 (GLU → GO) | low | text/scheme mismatch, both quoted |
| D12 | "Lowest rate constants" for k₂ holds only at 100 °C; cross-unit magnitude comparison invalid | medium | weakens the rate-determining-step claim (§8.5) |
| D13 | No purity/assay for commercial MGO and GO (usually ~40 % aqueous, oligomerised) | high | the most likely root cause of D2 |

---

## NEW-PARAMETER TABLE (consolidated)

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| k₁ (GLU → MGO) | 3.81 × 10⁻⁵ | 1/min | 100 °C, Glc 0.1 M / Pro 0.2 M, 20 mM phosphate pH 9.0, sealed | Table 1 row 1, journal p. 5 (PDF p. 5) | [F] |
| k₁ | 2.02 × 10⁻⁴ | 1/min | 120 °C, ditto | Table 1 row 1, p. 5 | [F] |
| k₁ | 0.0013 | 1/min | 140 °C, ditto | Table 1 row 1, p. 5 | [F] |
| Ea(k₁) | 112.87 | kJ/mol | 100–140 °C, T_av = 120 °C | Table 1 row 1, p. 5 | [F] |
| k₂ (GLU → GO) | 5.07 × 10⁻⁸ | 1/min | 100 °C | Table 1 row 2, p. 5 | [F] |
| k₂ | 0.0023 | 1/min | 120 °C | Table 1 row 2, p. 5 | [F] |
| k₂ | 0.0028 | 1/min | 140 °C | Table 1 row 2, p. 5 | [F] |
| Ea(k₂) | 355.28 | kJ/mol | 100–140 °C | Table 1 row 2, p. 5 | [F] |
| k₃ (GLU + PRO → INT) | 8.15 × 10⁻⁴ | L/(mol·min) | 100 °C | Table 1 row 3, p. 5 | [F] |
| k₃ | 0.026 | L/(mol·min) | 120 °C | Table 1 row 3, p. 5 | [F] |
| k₃ | 0.49 | L/(mol·min) | 140 °C | Table 1 row 3, p. 5 | [F] |
| Ea(k₃) | 205.32 | kJ/mol | 100–140 °C | Table 1 row 3, p. 5 | [F] |
| k₄ (INT → MGO + PRO) | 9.00 × 10⁻⁸ | 1/min | 100 °C | Table 1 row 4, p. 5 | [F] |
| k₄ | 4.29 × 10⁻⁵ | 1/min | 120 °C | Table 1 row 4, p. 5 | [F] |
| k₄ | 5.78 × 10⁻⁴ | 1/min | 140 °C | Table 1 row 4, p. 5 | [F] |
| Ea(k₄) | 282.64 | kJ/mol | 100–140 °C | Table 1 row 4, p. 5 | [F] |
| k₅ (INT → GO + PRO) | 1.165 | 1/min | 100 °C | Table 1 row 5, p. 5 | [F] ⚠️ D3 |
| k₅ | 5.17 × 10⁻⁸ | 1/min | 120 °C | Table 1 row 5, p. 5 | [F] ⚠️ D3 |
| k₅ | 1.83 × 10⁻⁴ | 1/min | 140 °C | Table 1 row 5, p. 5 | [F] ⚠️ D3 |
| Ea(k₅) | − 294.05 | kJ/mol | 100–140 °C | Table 1 row 5, p. 5 | [F] ⚠️ negative |
| **k₆ (MGO + PRO → 2-AP)** | **0.17** | L/(mol·min) | 100 °C | Table 1 row 6, p. 5 | [F] ⚠️ D2 |
| k₆ | 0.18 | L/(mol·min) | 120 °C | Table 1 row 6, p. 5 | [F] ⚠️ D2 |
| k₆ | 1.50 | L/(mol·min) | 140 °C | Table 1 row 6, p. 5 | [F] ⚠️ D2 |
| **Ea(k₆), 2-AP formation via MGO** | **68.61** | kJ/mol | 100–140 °C | Table 1 row 6, p. 5 | [F] |
| **k₇ (GO + PRO → 2-AP)** | **0.38** | L/(mol·min) | 100 °C | Table 1 row 7, p. 5 | [F] ⚠️ D2 |
| k₇ | 3.50 | L/(mol·min) | 120 °C | Table 1 row 7, p. 5 | [F] ⚠️ D2 |
| k₇ | 8.10 | L/(mol·min) | 140 °C | Table 1 row 7, p. 5 | [F] ⚠️ D2 |
| **Ea(k₇), 2-AP formation via GO** | **98.70** | kJ/mol | 100–140 °C | Table 1 row 7, p. 5 | [F] |
| **k₈ (2-AP → P₁), 2-AP depletion** | **1.77** | 1/min | 100 °C | Table 1 row 8, p. 5 | [F] |
| **k₈** | **40.75** | 1/min | 120 °C | Table 1 row 8, p. 5 | [F] |
| **k₈** | **239.71** | 1/min | 140 °C | Table 1 row 8, p. 5 | [F] |
| **Ea(k₈), 2-AP depletion** | **157.89** | kJ/mol | 100–140 °C | Table 1 row 8, p. 5 | [F] |
| k₉ (GLU → P₂) | 0.010 | 1/min | 100 °C | Table 1 row 9, p. 5 | [F] |
| k₉ | 0.042 | 1/min | 120 °C | Table 1 row 9, p. 5 | [F] |
| k₉ | 5.12 × 10⁻⁵ | 1/min | 140 °C | Table 1 row 9, p. 5 | [F] |
| Ea(k₉) | − 164.48 | kJ/mol | 100–140 °C | Table 1 row 9, p. 5 | [F] ⚠️ negative |
| k₁₀ (PRO → P₃) | 0.0060 | 1/min | 100 °C | Table 1 row 10, p. 5 | [F] |
| k₁₀ | 0.0058 | 1/min | 120 °C | Table 1 row 10, p. 5 | [F] |
| k₁₀ | 4.16 × 10⁻⁴ | 1/min | 140 °C | Table 1 row 10, p. 5 | [F] |
| Ea(k₁₀) | − 84.02 | kJ/mol | 100–140 °C | Table 1 row 10, p. 5 | [F] ⚠️ negative |
| k₇/k₆ ratio | 2.24 / 19.4 / 5.40 | dimensionless | 100 / 120 / 140 °C | derived from Table 1 rows 6–7 | [Z] |
| k₈/(k₆+k₇) | 3.21 / 11.1 / 24.9 | min·mol/L | 100 / 120 / 140 °C | derived from Table 1 rows 6–8 | [Z] |
| Ea check, all ten k | reproduce printed Ea to ≤ 5 % | kJ/mol | 100 ↔ 140 °C two-point | §4.1 | [Z] |
| 2-AP max yield | 9.95 | μmol/L | 100 °C, 20 min, glucose/proline | text §3.3, journal p. 4 | [M] |
| 2-AP time course, 100 °C | 0.7 / 2.6 / 6.6 / 6.65 / 9.95 / 8.4 / 7.4 | μmol/L | t = 2/5/10/15/20/25/30 min | Fig. 3E, journal p. 4 | [fig] |
| 2-AP time course, 120 °C | 4.05 / 6.2 / 5.1 / 4.1 / 2.65 / 1.6 / 1.15 | μmol/L | t = 2…30 min | Fig. 3E, p. 4 | [fig] |
| 2-AP time course, 140 °C | 2.45 / 0.80 / 0.73 / 0.89 / 0.62 / 0.46 / 0.38 | μmol/L | t = 2…30 min | Fig. 3E, p. 4 | [fig] |
| glucose time course, 100/120/140 °C | see §5.1 (21 points) | mmol/L | t = 2…30 min | Fig. 3A, p. 4 | [fig] |
| proline time course, 100/120/140 °C | see §5.2 (21 points) | mmol/L | t = 2…30 min | Fig. 3B, p. 4 | [fig] |
| MGO time course, 100/120/140 °C | see §5.3 (21 points) | μmol/L | t = 2…30 min | Fig. 3C, p. 4 | [fig] |
| GO time course, 100/120/140 °C | see §5.4 (21 points) | μmol/L | t = 2…30 min | Fig. 3D, p. 4 | [fig] |
| MGO range | 10–340 | μmol/L | all runs | text §3.2, p. 4 | [M] |
| GO range | 30–320 | μmol/L | all runs | text §3.2, p. 4 | [M] |
| R², glucose | 0.80 / 0.80 / 0.62 | — | 100 / 120 / 140 °C | Fig. 3A panel, p. 4 | [F] |
| R², proline | 0.48 / 0.32 / 0.42 | — | 100 / 120 / 140 °C | Fig. 3B panel, p. 4 | [F] |
| R², methylglyoxal | 0.98 / 0.93 / 0.97 | — | 100 / 120 / 140 °C | Fig. 3C panel, p. 4 | [F] |
| R², glyoxal | 0.83 / 0.92 / 0.64 | — | 100 / 120 / 140 °C | Fig. 3D panel, p. 4 | [F] |
| R², 2-AP | 0.93 / 0.82 / 0.54 | — | 100 / 120 / 140 °C | Fig. 3E panel, p. 4 | [F] |
| RRMSE, all responses | **not reported** (eq. 11 defined only) | — | — | §2.8, p. 4 | — |
| mass-balance recovery, 30 min | 80.41 / 69.50 / 54.40 | % | 100 / 120 / 140 °C | text §3.4, p. 4 (Fig. S1 unavailable) | [M] |
| 2-AP, fed MGO 10 mmol/L | 0.19 ± 0.012 | μmol/L | 100 °C, 20 min, PRO 20 mmol/L | Table 2, journal p. 6 | [M] |
| 2-AP, fed MGO 20 mmol/L | 0.12 ± 0.043 | μmol/L | ditto | Table 2, p. 6 | [M] |
| 2-AP, fed MGO 40 mmol/L | n.d. (< LOD) | μmol/L | ditto | Table 2, p. 6 | [M] |
| **2-AP, fed GO 10 mmol/L** | **0.43 ± 0.088** | μmol/L | ditto | Table 2, p. 6 | [M] |
| 2-AP, fed GO 20 mmol/L | 0.079 ± 0.0079 | μmol/L | ditto | Table 2, p. 6 | [M] |
| 2-AP, fed GO 40 mmol/L | n.d. (< LOD) | μmol/L | ditto | Table 2, p. 6 | [M] |
| model over-prediction of Table 2 | ≈100× (GO) and ≈101× (MGO) | dimensionless | 100 °C, 20 min | §6, derived from Table 1 + Table 2 | [Z] ⚠️ D2 |
| 2-AP, cooked rice control | 6.75 ± 0.14 | μmol/kg | steamed 100 ± 3 °C, 15 min + 20 min hold | Table 3, journal p. 6 | [M] |
| 2-AP, rice + 0.05 μmol/kg MGO | 9.01 ± 0.54 | μmol/kg | ditto | Table 3, p. 6 | [M] |
| 2-AP, rice + 50 μmol/kg MGO | 7.92 ± 0.094 | μmol/kg | ditto | Table 3, p. 6 | [M] |
| 2-AP, rice + 0.05 μmol/kg GO | 7.18 ± 0.89 | μmol/kg | ditto | Table 3, p. 6 | [M] |
| 2-AP, rice + 50 μmol/kg GO | 7.79 ± 0.31 | μmol/kg | ditto | Table 3, p. 6 | [M] |
| MGO 0.05 μmol/kg enhancement | +33.5 % | % | cooked rice | text §3.5, p. 6 (and Table 3) | [M] |
| GO 0.05 / GO 50 / MGO 50 enhancement | +6.4 / +15.4 / +17.3 | % | cooked rice | derived from Table 3 | [Z] |
| 2-AP odour threshold in starch | 0.0073 | μg/kg | starch matrix | Introduction, p. 1 | [C] |
| MGO FEMA GRAS level | 0.05–3 | mg/kg | baked food | §3.5, p. 6 | [C] (Hall & Oser 1965) |
| rice matrix pH | 6.0–7.0 | pH | cooked rice | §3.5, p. 6 | [M/C] |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

The natural cut axis for this paper is **temperature** (100 / 120 / 140 °C), and there is a second,
much sharper cut available along **system** (in-situ dicarbonyl generation vs fed dicarbonyl), and a
third along **matrix** (aqueous model vs cooked rice).

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 1, k₆ / k₇ / k₈ at 100 °C and 140 °C, plus their three Ea (68.61, 98.70, 157.89 kJ/mol)** | **FIT — as Ea priors and as a k₇/k₆ ratio prior ONLY; NOT as absolute magnitudes** | — | These are the corpus's only 2-AP formation *and* depletion barriers. The Ea and the ratio survive D2; the magnitudes do not. Ingest with `absolute_rate: false`, `prior_only: true`, and an assumed σ (nothing in the paper bounds it, D1) |
| **Table 1, all ten k at 120 °C** | **★ HOLD-OUT** | **temperature (interpolation)** | 120 °C is the Arrhenius reference T_av of eq. (12), so it is the *most* strongly determined column — which makes it the least informative to fit and the best interpolation test. A model given the 100 and 140 °C ends must reproduce the middle. Two rows are pathological (k₅ = 5.17 × 10⁻⁸, k₉ = 0.042) and a model that misses those is telling us something true |
| **Fig. 3E 2-AP time course, 100 °C (7 points)** | **FIT** | temperature | the maximum-yield arm, the one with the best fit (R² = 0.93), and the arm anchored by the one text-printed value (9.95 μmol/L). Ingest as `absolute_concentration: uncertain` (D5) |
| **★ Fig. 3E 2-AP time courses, 120 °C and 140 °C (14 points)** | **★ HOLD-OUT** | **temperature (extrapolation upward)** | This is the sharpest test the paper offers. The entire scientific claim is that 2-AP *accumulation* peaks at the LOWEST temperature because k₈ outruns k₆ + k₇ — a **sign reversal in the temperature response of the observable**. A model that only knows the 100 °C arm and the two formation Ea must predict the 120 °C peak-then-crash (6.2 μmol/L at 5 min → 1.15 at 30 min) and the 140 °C near-total suppression (≤ 2.45 μmol/L throughout). No fitted continuous term can fake a maximum-versus-temperature inversion |
| **Fig. 3C + 3D MGO and GO time courses, 100 °C (14 points)** | **FIT** | temperature | the α-dicarbonyl pool the 2-AP node divides by; R² 0.98 / 0.83 at this temperature. Needed to make the 120/140 °C 2-AP hold-out a *prediction* rather than a guess |
| **★ Fig. 3D glyoxal 120 °C early rise (283 μmol/L at 2 min, falling to 86 at 30 min)** | **★ HOLD-OUT** | temperature × time shape | The glyoxal peak *moves earlier and higher* from 100 → 120 °C and then *collapses* at 140 °C. It is a build-up/decay shape with a moving maximum — the classic multiresponse test. ⚠️ conflicts with the text's stated 320 μmol/L upper bound (D8); score with a ±15 % tolerance |
| **★ Table 2 in full (8 cells, fed MGO/GO at 0/10/20/40 mmol/L, 100 °C, 20 min)** | **★★ HOLD-OUT — the highest-value row in this dossier** | **system: in-situ vs fed dicarbonyl** | Same temperature, same time, same buffer, same vessel as the fitted arm — but the dicarbonyl is *fed* rather than generated. §6 shows the authors' own model over-predicts it by **≈100× for both compounds**. A model that fits the in-situ arm and then predicts Table 2 correctly has discovered something the original authors did not; a model that reproduces the 100× error has honestly inherited the defect. Either outcome is informative, which is exactly what a hold-out is for. **⚠️ Circularity risk if any wave uses Table 2 to re-scale k₆/k₇ — that would convert the test into a fit. Forbid.** |
| **Table 3 (cooked rice, 5 distinct samples)** | **HOLD-OUT — directional only, not scored on level** | **matrix: aqueous pH 9 model → rice pH 6–7** | A 1000× dose change moves the response ±15 % non-monotonically, in *opposite* directions for MGO and GO. There is no level here worth scoring, but the **sign** of the MGO 0.05 μmol/kg effect (+33.5 %, letter "a", p < 0.05) is a real matrix-transfer prediction. Score sign + significance only. Also note the pH jump (9.0 → 6.0–7.0) makes this the paper's only pH datapoint, and it is confounded with matrix |
| **Table 1, k₅ / k₉ / k₁₀ and their negative Ea** | **neither** | — | Three negative Ea, one of them attached to a rate constant that is non-monotonic by seven orders of magnitude (D3). These are error sinks for a lumped INT pool the authors admit is unmeasurable. Record as an `audit_flag`; do not fit, do not score |
| **Fig. 3B proline (all 21 points)** | **neither** | — | Non-monotonic in two of three arms, R² = 0.32–0.48, and the initial condition is internally contradictory (D6, 200 vs 232 mmol/L). Not usable as a level. Record the **1.7× glucose-vs-proline depletion-rate ratio** (§9) as a shape constraint only |
| Mass-balance recovery 80.41 / 69.50 / 54.40 % at 30 min | **FIT (as a one-sided ceiling on the melanoidin sink)** | temperature | The authors attribute the whole deficit to melanoidin formation and cite Martins & Van Boekel 2005 for it. It bounds the trunk's browning flux from above at three temperatures — a constraint the Module-4 browning lane has nothing else for. Cross-reference the `brands2002b` and `martins2003c` ε dossiers |

**Coverage note.** This paper opens a node the corpus does not currently have — an
**α-acetyl-N-heterocycle (2-AP) formation/depletion pair with separate Ea** — and it does so in the
same glucose + amino-acid chemistry that Module 4's trunk already covers. Recommend it be filed
against **Module 4 (trunk/browning)** for the mass-balance and dicarbonyl rows, and against a
**2-AP / proline-specific flavour node** for k₆, k₇, k₈. Nothing in the current declaration overlaps
it, so no disjointness violation is created.

**Circularity risks, explicit.**
1. Table 1 is *derived from* Fig. 3 by the authors' own fit. **Fitting Table 1 constants AND scoring
   against Fig. 3 data is circular.** Choose one: either ingest Table 1 as priors and hold out Fig. 3,
   or fit Fig. 3 directly and use Table 1 only as a sanity check. The proposal above takes the second
   route for the 100 °C arm and the first route for the Ea.
2. Ea(k₆), Ea(k₇), Ea(k₈) are reproducible to ≤ 2 % from the printed 100 and 140 °C constants (§4.1).
   Fitting the 100 and 140 °C constants and *also* fitting the Ea double-counts the same three
   numbers. Fit one, derive the other.
3. Table 2 must never be used to re-scale k₆/k₇ (see the ★★ row above).
