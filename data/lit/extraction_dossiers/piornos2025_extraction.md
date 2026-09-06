# Piornos, Balagiannis, Koussissi, Bekkers, Vissenaekens, Brouwer & Parker 2025 (10.1016/j.foodchem.2024.141532) — per-paper extraction 2026-09-07

**Source PDF:** `data/articles/Piornos2025.pdf` (11 pp., 1,449,284 bytes, SHA-256 `f1462091bdda2d0590ed887cf0d8be5137fe822f313b99f93d6cfd0e23bdcc80`). Born-digital Elsevier PDF (Acrobat Distiller 8.1.0), clean text layer.
Read method: **both** — full text layer (`pdftotext -layout`) read end to end, **plus** 200 dpi rasters of pp. 4, 6, 7, 8, 9 (`pp-04.png` … `pp-09.png`) for Fig. 1 (concentration–time curves), Fig. 2 (kinetic scheme), Table 1 (re-verified cell by cell), Fig. 3 (Arrhenius plot), Fig. 4 (the fifteen R² values, which live only in the raster) and Fig. 5 (chemical mechanism). **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ HEADLINE: the only multi-response Maillard fit in the corpus that runs INSIDE a real cereal matrix at 65–90 °C. Eight (k′ at 70 °C, Ea) pairs, five branching fractions. **All Ea are FIXED, not estimated — no confidence interval exists for any of them.** Ea(ARP → dicarbonyl pool) = 315 kJ/mol and Ea(Strecker step) = 255 kJ/mol are 2–3× any aqueous value in the corpus (Martins 2005: 97–125 kJ/mol) and are almost certainly moisture-confounded (§8).

---

## §0. IDENTITY

| item | value | how verified |
|---|---|---|
| **title** | ***"Multi-response kinetic modelling of the formation of five Strecker aldehydes during kilning of barley malt"*** | p.1 text + `pdfinfo` Title |
| **authors** | José A. Piornos^a,1, Dimitris P. Balagiannis^a,* (corresponding), Elisabeth Koussissi^b,2, August Bekkers^b, Johan Vissenaekens^c, Eric Brouwer^b, Jane K. Parker^a | p.1 |
| affiliations | ^a Dept. Food and Nutritional Sciences, University of Reading, RG6 6DZ, UK; ^b Heineken Supply Chain BV, Global Innovation & Research, Zoeterwoude, NL; ^c Mouterij Albert, Ruisbroek-Sint-Amands, BE. Current addresses: Piornos → CSGA, CNRS/INRAE, Dijon; Koussissi → Univ. of West Attica | p.1 footnotes |
| venue | *Food Chemistry* **464** (2025) 141532 | header every page; `pdfinfo` Subject |
| **DOI** | **10.1016/j.foodchem.2024.141532** | p.1 |
| dates | Received 12 July 2024; revised 28 September 2024; accepted 2 October 2024; online 5 October 2024 | p.1 |
| licence | open access, **CC BY** | p.1 |
| **funding** | "This study has been fully funded by Heineken Supply Chain BV." Writing retreat funded by the Agriculture, Food and Health research theme, University of Reading | p.10 Acknowledgments |
| **conflicts** | declared: Piornos (support, equipment from Heineken and Mouterij Albert); Koussissi, Bekkers, Brouwer (Heineken employment); Vissenaekens (Mouterij Albert employment) | p.10 |
| data availability | "Data will be made available on request." | p.10 |
| **supplementary** | Table S1 (green-malt moisture), Table S2 (MRM settings), Table S3 (ANOVA), Fig. S1 (full concentration data, panels a–s), Fig. S2 (predicted Int1) — **[NEG] NOT on disk** (`data/articles/` holds no Piornos SI) | `ls data/articles` |
| software | Athena Visual Studio v14.2 (AthenaVISUAL Inc.); SPSS 22 | p.3, p.4 |
| PDF character | born-digital; Table 1 and all equations are in the text layer; Fig. 3/4 numerics only in raster | — |

**Correct file for its expected identity.** Text-layer Table 1 vs raster: all 18 rows identical (`[F]` cells below are double-checked).

---

## §1. ONE-PARAGRAPH VERDICT — READ THIS BEFORE USING ANY NUMBER HERE

This is a Balagiannis/Parker-school multi-response fit (same Athena software, determinant criterion and reparametrised Arrhenius form as Martins & Van Boekel 2005, which the repository already uses for its trunk) applied to **dried green malt cured isothermally at 65, 78 and 90 °C for 8.4 h** in a pilot micro-malting kiln. It quantifies free amino acids, glucose, fructose, seven Amadori compounds (FruLeu/Ile/Val/Phe/Ala/Gly/Pro) and five Strecker aldehydes over time, and fits a **glucose + AA → ARP → Int1 (a lumped short-chain-dicarbonyl pool) → Strecker aldehyde** network with fructose feeding Int1 directly. What it delivers is **eight rate constants at T_ref = 70 °C with 95 % HPD intervals** (Table 1) and **five branching fractions F_AA**. What it does **not** deliver is any estimated activation energy: **all four Ea (175, 315, 255, 121 kJ/mol) were fixed in the final run** and carry no interval; they were obtained in alternating runs (fix k′, fit Ea; fix Ea, fit k′) that the paper does not report. Three further limits: (i) the system is a **low-moisture, drying solid** — moisture of the malt during curing is never stated in the main text (only pre-kiln green-malt moisture, in the missing Table S1), so the Ea values fold together chemistry and water-activity loss; (ii) **Int1 is never measured** — no dicarbonyl analytics of any kind — so k′_M has a 92 % error and Ea_M = 255 is the temperature dependence of an unobserved pool; (iii) **no sulfur chemistry beyond methional** — methionine is pooled under k_i1 and FruMet is unmeasured (model-estimated). The paper's own ranking of temperature sensitivity is the one usable qualitative result: **ARP cleavage (Ea₂ = 315) > Strecker step on the dicarbonyl pool (Ea_M = 255) > glucose–AA condensation (Ea₁ = 175) > fructose degradation (Ea_F = 121)**. Use the k′ values only as 65–90 °C solid-matrix anchors and the Ea only as *ordinal* priors (§8).

---

## §2. SYSTEM — verbatim `[M]`

### 2.1 Material and kiln (§2.1, p.2)

> "Barley was received at Mouterij Albert (Ruisbroek-Sint-Amands, Belgium), steeped in water for one day and germinated for five days at industrial scale. Two different varieties of barley were used in this study: the two-row spring variety 'RGT Planet' and six-row winter variety 'Etincel'. The germinated grains, i.e., the green malt, were kilned using pilot-scale micro-malting equipment from Nordon & Cie. (Nancy, France). Moisture of the green malt samples was measured before the kilning process (Supplementary Table S1). The green malt was placed in cubic shape stainless-steel baskets (150 mm side), split diagonally into two parts by a piece of steel (650 ± 2 g in each side). The micro-malting equipment was provided with eight baskets with grilled bottom to allow hot air to circulate throughout."

### 2.2 Drying and curing programme (§2.1, p.2)

> "The kilning programmes had an initial drying process, starting at 25 °C and reaching 55 °C in 10 min, then increasing to 64 °C in 45 min, kept constant for 4 h and 50 min and raised to 65 °C in 3.25 h. After the drying process, the temperature was increased to the curing temperature (65 °C, 78 °C or 90 °C) in 10 min and kept constant for 8.4 h. Sampling was done only during the curing stage, since no formation of aroma compounds was observed during the drying stage. The baskets were taken out randomly from the different positions in the oven and the empty space replaced by an empty basket with a lid. The total duration of the malting experiments was 16 h and 52 min, and the samples were taken every 72 min only during the curing stage. After kilning, the rootlets were removed by manual rubbing, separated by sieving through a 1.8 × 23 mm mesh and stored in a freezer at −30 ± 1 °C to limit thermal reactions. The kilning experiments were performed in duplicate from two different batches of barley for each variety in different days, except for the experiments at 65 °C, where the duplicates were from the same batch due to availability at the industry on the day of collection."

| variable | value | tag |
|---|---|---|
| Curing temperatures | **65, 78, 90 °C** (isothermal, reached in 10 min from 65 °C) | [F] |
| Curing duration | **8.4 h** | [F] |
| Sampling interval | **every 72 min** → t = 0, 1.2, 2.4, 3.6, 4.8, 6.0, 7.2, 8.4 h = **8 points per run** (Fig. 1 x-axis ticks confirm) | [F]/[D] |
| Drying stage preceding curing | 25 → 55 °C (10 min) → 64 °C (45 min), hold 4 h 50 min, → 65 °C over 3.25 h; total programme 16 h 52 min | [F] |
| Runs | 2 varieties × 3 temperatures × duplicate = **12 kilning runs** | [D] |
| Vessel | 150 mm stainless baskets, 650 ± 2 g green malt per half-basket, 8 baskets, hot-air through-flow | [F] |
| Moisture | **only pre-kiln green-malt moisture, in Table S1 (not on disk)**. Moisture, a_w or temperature-in-bed during curing: **[NEG] not reported anywhere in the main text** | [NEG] |
| pH | **[NEG] not reported** | [NEG] |
| Temperature homogeneity | authors admit "non-homogeneous temperature of the air and flowrates through the boxes" (p.5) as a variability source | [F] |

---

## §3. ANALYTES AND QUANTIFICATION — verbatim `[M]`

### 3.1 Non-volatile extraction (§2.3, p.2–3)

> "The ground samples (1.0 g) were extracted using 10 mL of ultrapure water (18.2 MΩ) containing 1.25 mM of L-norvaline and 15 μM of trehalose as internal standards for amino acids and sugars, respectively. … centrifugation at 5500 ×g for 15 min at 4 °C … the pellet was reextracted twice (5 mL × 2). … The extractions were performed in duplicate."

| analyte class | method | internal standard | calibration | tag |
|---|---|---|---|---|
| Free amino acids (18) | HILIC LC-MS/MS (Synchronis HILIC 150 × 4.6 mm, dynamic MRM, ESI+) | **L-norvaline** | 0–2.5 mM standard of 18 amino acids | [M] |
| Glucose, fructose | LC-MS/MS (LUNA Omega SUGAR, ESI−, MRM) | **trehalose** | 0–1 mM | [M] |
| ARP: FruVal, FruLeu, FruIle, FruPhe, FruAla, FruGly, FruPro | LC-MS/MS (Discovery HS F5-3, 55 °C, ESI+, MRM); extracts diluted 50× | **"no internal standard was used for ARP"** | 0–1000 μg/L in ultrapure water (authentic standards, Toronto Research Chemicals, 95–97 %) | [M] |
| **FruMet** | **[NEG] not quantified** — "FruMet was not quantified instrumentally and for this reason was estimated by the model" (p.6) | — | — | [NEG] |
| **Dicarbonyls / Int1** | **[NEG] not measured at all** — "These intermediates were not quantified in this study, but we assumed that they corresponded to a pool of SCDC, such as glyoxal and methylglyoxal" (p.6) | — | — | [NEG] |

### 3.2 Strecker aldehydes (§2.4, p.3)

> "Ground malt samples (1.0 g for experiments at 65 and 78 °C; 0.5 g for 90 °C) were weighed in 20-mL screw-capped SPME vials, together with 5 mL of saturated NaCl aqueous solution and 5 μL of internal standard solution (100 mg/L of 2-methylpentanal and 100 mg/L of 2-methylbenzaldehyde in absolute ethanol). 2-Methylpentanal was used as internal standard for 2-methylpropanal, 2-methylbutanal, 3-methylbutanal, and methional; and 2-methylbenzaldehyde for phenylacetaldehyde. The samples were incubated at 50 °C for 10 min and then a PDMS/DVB/Carboxen® SPME fibre was exposed to the headspace of the vial for 20 min. … The standards for calibration (0–1000 μg/L) were spiked in freeze-dried green malt in order to account for the matrix effects on the release of the volatiles to the headspace of the samples. The analyses were performed in duplicate."

Quantifier ions (p.3): m/z 41 (2-methylpropanal), 41 (2-/3-methylbutanal), 58 (2-methylpentanal IS), 48 (methional), 91 (phenylacetaldehyde and 2-methylbenzaldehyde IS). GC: ZB-5MSi 30 m × 0.25 mm × 1 μm, 50 °C (2 min) → 300 °C at 6 °C/min, SIM.

**Absolute basis:** every concentration and every rate constant is in **mmol per kg of ground kilned malt** (Fig. 1, Fig. 4 caption "all in mmol/kg"; second-order k in kg mmol⁻¹ h⁻¹). **[NEG] Whether that kg is fresh or dry weight is never stated.** Matrix-matched calibration for volatiles (spiked freeze-dried green malt) is a genuine strength; the ARP calibration is in plain water with no IS.

### 3.3 Magnitudes on the raster (Fig. 1, p.4; Fig. 4 axes, p.8) `[M]`

| species | range observed | source |
|---|---|---|
| Glucose | ≈ 5–45 mmol/kg (Fig. 4a axis 0–60) | Fig. 4 |
| Leu / Ile / Val / Phe / Met | ≈ 4.5–9 / 5–9 / 12–27 / 10–21 / 5.3–8.2 mmol/kg | Fig. 4b–f |
| FruLeu, FruIle, FruVal, FruPhe | start ≈ 0.10–0.20 mmol/kg; 78 °C rises monotonically to ≈ 0.5–0.8 by 8.4 h; **90 °C peaks at ≈ 3.6–4.8 h (FruLeu ≈ 0.78, FruPhe ≈ 0.77, FruVal ≈ 0.51, FruIle ≈ 0.47) then falls** | Fig. 1a–d |
| 3-MB / 2-MB / 2-MP / PhAc / methional at 90 °C, 8.4 h (run 90S1) | ≈ 0.46 / 0.22 / 0.20 / 0.18 / 0.12 mmol/kg; at 65 °C all aldehydes stay ≈ 0 | Fig. 1e–i |
| Int1 (model only) | "maximal predicted level of around 0.17 mmol/kg (Supplementary Fig. S2)" | p.6 [F] |

---

## §4. THE REACTION SCHEME — every step `[F]`

Backbone stated in prose (p.5): "Glucose + Free Amino Acids —k1→ ARP —k2 (− Free Amino Acids)→ Int1 —k3 (+ Strecker Amino Acids)→ Strecker Aldehydes". Final scheme is Fig. 2 (p.6) and the ODE system (p.7). Underlined (measured) species in Fig. 2: Glucose, Fructose, Leu, Ile, Val, Phe, Met, AAi, FruLeu, FruIle, FruVal, FruPhe, FruAAi, and the five aldehydes. Not measured: FruMet, AAj, FruAAj, Int1, MRP, SDP.

| # | step | rate law (p.7) | constant | tag |
|---|---|---|---|---|
| 1a | Glucose + Leu → FruLeu | k₂₁[Gluc][Leu] | k₂₁ | [F] |
| 1b | Glucose + Ile → FruIle | k₃₁[Gluc][Ile] | k₃₁ | [F] |
| 1c | Glucose + Val → FruVal | k₄₁[Gluc][Val] | k₄₁ | [F] |
| 1d | Glucose + Phe → FruPhe | k₅₁[Gluc][Phe] | k₅₁ | [F] |
| 1e | Glucose + Met → FruMet; Glucose + AAi → FruAAi; Glucose + AAj → FruAAj | k_i1[Gluc][X] (one shared constant) | k_i1 | [F] |
| 2 | each ARP → Int1 + regenerated amino acid (7 ARPs, one shared constant: k₂₂ = k₃₂ = k₄₂ = k₅₂ = k_i2 ≡ k₂) | k₂[ARP] | k₂ | [F] |
| F | Fructose → Int1 (single step, first order) | k_F[Fruc] | k_F | [F] |
| 3a | Int1 + Leu → 3-methylbutanal + SDP | k_M F_Leu [Int1][Leu] | k_M, F_Leu | [F] |
| 3b | Int1 + Ile → 2-methylbutanal + SDP | k_M F_Ile [Int1][Ile] | k_M, F_Ile | [F] |
| 3c | Int1 + Val → 2-methylpropanal + SDP | k_M F_Val [Int1][Val] | k_M, F_Val | [F] |
| 3d | Int1 + Phe → phenylacetaldehyde + SDP | k_M F_Phe [Int1][Phe] | k_M, F_Phe | [F] |
| 3e | Int1 + Met → methional + SDP | k_M F_Met [Int1][Met] | k_M, F_Met | [F] |
| 4 | Int1 + any amino acid → MRP (the (1 − F_AA) share of Leu/Ile/Val/Phe/Met, and 100 % of AAi, AAj) | k_M (1 − F_AA)[Int1][AA]; k_M[Int1]([AAi]+[AAj]) | k_M | [F] |

Full ODEs as printed (p.7), e.g. d[Gluc]/dt = −[Gluc]·(k₂₁[Leu] + k₃₁[Ile] + k₄₁[Val] + k₅₁[Phe] + k_i1([Met]+[AAi]+[AAj])); d[Fruc]/dt = −k_F[Fruc]; d[Leu]/dt = −k₂₁[Gluc][Leu] + k₂₂[FruLeu] − k_M[Int1][Leu]; d[FruLeu]/dt = k₂₁[Gluc][Leu] − k₂₂[FruLeu]; d[Int1]/dt = k_F[Fruc] + k₂₂[FruLeu] + k₃₂[FruIle] + k₄₂[FruVal] + k₅₂[FruPhe] + k_i2·([FruMet]+[FruAAi]+[FruAAj]) − k_M[Int1]·([Leu]+[Ile]+[Val]+[Phe]+[Met]+[AAi]+[AAj]); d[3MB]/dt = k_M F_Leu[Int1][Leu]; d[MRP]/dt = k_M[Int1]·{(1−F_Leu)[Leu] + (1−F_Ile)[Ile] + (1−F_Val)[Val] + (1−F_Phe)[Phe] + (1−F_Met)[Met] + [AAi] + [AAj]}; d[SDP]/dt = k_M[Int1]·(F_Leu[Leu] + F_Ile[Ile] + F_Val[Val] + F_Phe[Phe] + F_Met[Met]). Note the ODEs write k₂₂, k₃₂ … k_i2 separately but the text (p.7) states "For all ARP, the kinetic rate constants for the degradation reaction were the same (k₂)". **Glucose is consumed only via ARP; there is no direct glucose → dicarbonyl step, no sugar isomerisation, no reversibility, no loss term for aldehydes.**

Pooling (p.6): AAi = Ala + Gly + Pro (their ARPs measured, FruAAi = FruAla + FruGly + FruPro); AAj = all remaining amino acids (ARP unmeasured, model-estimated).

**Routes tested and rejected (p.5)** `[F]`: direct ARP → Strecker aldehyde + SDP without Int1 (Cremer et al. 2000 / Yaylayan 2003 route) — "AIC = 672" vs "AIC = 639" for the Int1 model; degradation of ARP, Int1 or aldehydes into unquantified species; a second intermediate before Int1; MRP from amino acids alone; individual (unshared) constants — all "either did not produce any relevant improvement … or the parameters related to them were null". Adding the one-step fructose → Int1 route "improved the fit of the model in terms of a lower sum of squared residuals".

**Chemical interpretation (Fig. 5, p.9)** `[F]`: glucose + AA → Schiff base (−H₂O) → ARP → (releases amino acid) → **deoxyglucosones** (3-DG and 1-DG drawn) → **SCDC** (R₁-CO-CO-R₂) ← fructose (direct arrow); SCDC + amino acid → Strecker aldehyde + amino ketone (−CO₂), and SCDC + AA → "colour and other aroma compounds". Text (p.10): "the sugar backbone forms a vicinal dicarbonyl, such as 3-deoxyglucosone from glucose … These dicarbonyls, in turn, break down further into a pool of reactive SCDC, such as glyoxal, methylglyoxal, 2,3-butanedione". So **Int1 lumps deoxyosones AND their C2–C4 fragments into one pool**.

---

## §5. FITTED PARAMETERS — Table 1 re-typed exactly (p.6) `[F]`

Arrhenius form (eq. 1, p.3): k = k′·exp[(Ea/R_g)(1/T_ref − 1/T)], **T_ref = 343.15 K (70.00 °C)**, R_g = 8.314 J mol⁻¹ K⁻¹, Ea in J mol⁻¹. Estimation: Bayesian multi-response, diagonal covariance, determinant criterion (Box & Draper 1965), Athena Visual Studio 14.2; discrimination by SS and AIC = n ln(SS/n) + 2(p+1).

| Parameter | used in | Optimal estimate ± 95 % HPD (% error) | tag |
|---|---|---|---|
| k′₂₁ (kg mmol⁻¹ h⁻¹) | k₂₁ | **1.05·10⁻⁴ ± 7·10⁻⁶ (7 %)** | [F] |
| k′₃₁ (kg mmol⁻¹ h⁻¹) | k₃₁ | **6.05·10⁻⁵ ± 4·10⁻⁶ (7 %)** | [F] |
| k′₄₁ (kg mmol⁻¹ h⁻¹) | k₄₁ | **2.57·10⁻⁵ ± 1·10⁻⁶ (6 %)** | [F] |
| k′₅₁ (kg mmol⁻¹ h⁻¹) | k₅₁ | **5.08·10⁻⁵ ± 4·10⁻⁶ (7 %)** | [F] |
| k′_i1 (kg mmol⁻¹ h⁻¹) | k_i1 | **1.16·10⁻⁵ ± 7·10⁻⁷ (6 %)** | [F] |
| **Ea₁ (kJ mol⁻¹)** | k₂₁, k₃₁, k₄₁, k₅₁, k_i1 | **175 (fixed)** | [F] |
| k′₂ (h⁻¹) | k₂ | **7.30·10⁻⁴ ± 3·10⁻⁵ (5 %)** | [F] |
| **Ea₂ (kJ mol⁻¹)** | k₂ | **315 (fixed)** | [F] |
| F_Leu | — | **1.000 (upper bound)** | [F] |
| F_Ile | — | **0.518 ± 0.03 (6 %)** | [F] |
| F_Val | — | **0.172 ± 0.01 (6 %)** | [F] |
| F_Phe | — | **0.207 ± 0.01 (5 %)** | [F] |
| F_Met | — | **0.269 ± 0.02 (9 %)** | [F] |
| k′_M (kg mmol⁻¹ h⁻¹) | k_M | **5.19·10⁻⁴ ± 5·10⁻⁴ (92 %)** | [F] |
| **Ea_M (kJ mol⁻¹)** | k_M | **255 (fixed)** | [F] |
| k′_F (h⁻¹) | k_F | **1.82·10⁻² ± 3·10⁻³ (18 %)** | [F] |
| **Ea_F (kJ mol⁻¹)** | k_F | **121 (fixed)** | [F] |

\* HPD = highest posterior density. 13 estimated parameters (5 k′₁-type + k′₂ + 5 F + k′_M + k′_F); "≤10 % of the actual estimated value for 10 out of 13 parameters" (p.9) — consistent: F_Leu is at its bound, k′_F 18 %, k′_M 92 %.

**On the Ea (p.9, verbatim):** "The activation energies, Ea, were kept fixed at the last calculation run in order to reduce the system's HPD intervals, hence no confidence intervals are reported for them. Activation energies were found to be in the range between 121 and 315 kJ mol⁻¹ (Table 1), generally higher than similar kinetic studies. In our study, the activation energy related to the formation of Strecker aldehydes, EaM, was 255 kJ mol⁻¹, but Ea2, for the cleavage of ARP, was even higher (315 kJ mol⁻¹)." And: "The estimation was performed in several runs, alternatively keeping one group of parameters fixed (k′ or Ea) and estimating the rest." **[NEG] The Ea-estimating runs, their intervals and their starting values are nowhere reported.**

**On k′_M (p.9):** "The high uncertainty of k′M was due to the fact that it was associated with the degradation of Int1 whose concentration was not determined analytically."

**HPD rounding check** `[Z]`: printed half-widths are one significant figure, so the % error column is the more precise datum — k′₄₁: 1·10⁻⁶/2.57·10⁻⁵ = 3.9 % vs printed 6 %; k′₅₁: 7.9 % vs 7 %; k′_M: 96 % vs 92 %; k′_F: 16.5 % vs 18 %; F_Met: 7.4 % vs 9 %. All others agree within rounding. Use the printed % when reconstructing an interval.

### 5.1 Temperature sensitivity — the paper's ranking `[F]` and Fig. 3 reconstruction `[Z]`

Verbatim (p.9): "since the intermediate steps (k2, kM) showed higher Eas than the first steps (k1, kF), and as processing temperature increased, the degradation rate of the intermediate compounds increased to a larger extent than that of the precursors. Therefore, for the same amount of degraded initial precursors, a higher amount of Strecker aldehydes was formed in relation to lower temperatures. … also due to the considerable reduction of k2′ in relation to the other reaction rate constants since Ea2 has the highest value." Abstract: "degradation of Amadori rearrangement products and short-chain dicarbonyls was more sensitive to temperature change due to their higher activation energies compared to other kinetic steps." Conclusions: "lower temperatures restrain the degradation of the Amadori products in particular".

k(T) from Table 1 via eq. 1 `[Z]` (these reproduce the Fig. 3 tick positions on the raster — e.g. k₂ ≈ 1.4·10⁻⁴ at 65 °C and ≈ 3·10⁻¹ at 90 °C; k_F ≈ 1·10⁻² → 1.9·10⁻¹):

| constant | Ea | k(65 °C) | k(78 °C) | k(90 °C) | **k(100 °C) extrapolated** | k(90)/k(65) | Q₁₀ (65→75 °C) |
|---|---|---|---|---|---|---|---|
| k₂₁ (kg mmol⁻¹ h⁻¹) | 175 | 4.24·10⁻⁵ | 4.25·10⁻⁴ | 3.08·10⁻³ | 1.46·10⁻² | 73 | 6.0 |
| k₃₁ | 175 | 2.44·10⁻⁵ | 2.45·10⁻⁴ | 1.77·10⁻³ | 8.38·10⁻³ | 73 | 6.0 |
| k₄₁ | 175 | 1.04·10⁻⁵ | 1.04·10⁻⁴ | 7.53·10⁻⁴ | 3.56·10⁻³ | 73 | 6.0 |
| k₅₁ | 175 | 2.05·10⁻⁵ | 2.05·10⁻⁴ | 1.49·10⁻³ | 7.04·10⁻³ | 73 | 6.0 |
| k_i1 | 175 | 4.68·10⁻⁶ | 4.69·10⁻⁵ | 3.40·10⁻⁴ | 1.61·10⁻³ | 73 | 6.0 |
| **k₂ (h⁻¹), ARP → Int1** | **315** | 1.43·10⁻⁴ | 9.03·10⁻³ | 3.19·10⁻¹ | **5.23** | **2240** | **25** |
| **k_M (kg mmol⁻¹ h⁻¹), Int1 + AA** | **255** | 1.38·10⁻⁴ | 3.98·10⁻³ | 7.13·10⁻² | **0.685** | **515** | **13.5** |
| k_F (h⁻¹), Fru → Int1 | 121 | 9.72·10⁻³ | 4.78·10⁻² | 1.88·10⁻¹ | 0.551 | 19 | 3.4 |

ARP half-life 1/k₂ · ln 2 `[Z]`: **4 860 h at 65 °C, 77 h at 78 °C, 2.2 h at 90 °C, 8 min at 100 °C**. A Q₁₀ of 25 for a single elementary step is physically implausible for pure chemistry; it is the signature of a *fixed* Ea absorbing the moisture drop of the bed as T rises (§8). Note that at 65 °C the model has ARP essentially inert (t½ ≈ 200 days), which is what makes the 65 °C aldehyde curves flat (Fig. 1).

---

## §6. GOODNESS OF FIT `[F]` (Fig. 4, p.8 raster; text p.9)

| panel | species | linear fit predicted = a·observed + b | R² |
|---|---|---|---|
| a | Glucose | y = 0.303x + 15.699 | **0.2195** |
| b | Leucine | y = 0.4654x + 3.1494 | 0.3697 |
| c | Isoleucine | y = 0.4872x + 3.1995 | 0.4493 |
| d | Valine | y = 0.4466x + 9.4952 | 0.4331 |
| e | Phenylalanine | y = 0.3971x + 8.2695 | 0.378 |
| f | Methionine | y = 0.4647x + 3.2789 | 0.3908 |
| g | FruLeu | y = 0.8128x + 0.076 | 0.8667 |
| h | FruIle | y = 0.7856x + 0.0555 | 0.8349 |
| i | FruVal | y = 0.8512x + 0.0448 | 0.8914 |
| j | FruPhe | y = 0.7827x + 0.0987 | 0.8322 |
| k | 3-Methylbutanal | y = 0.9663x − 0.0011 | 0.9501 |
| l | 2-Methylbutanal | y = 0.9354x + 0.0043 | 0.9465 |
| m | 2-Methylpropanal | y = 0.9327x + 0.0046 | 0.9458 |
| n | Phenylacetaldehyde | y = 0.9789x − 0.0012 | 0.9558 |
| o | Methional | y = 0.7996x + 0.0077 | 0.8507 |

Text (p.9): "R² close to or higher than 0.80 for ARP and around 0.95 for all Strecker aldehydes apart from methional (0.81)" — **⚠ the Fig. 4o panel prints R² = 0.8507 for methional, not 0.81** [D]; minor internal inconsistency. "The residual plots (not shown) were checked, and they were randomly scattered." AIC of the accepted model **639** (vs 672 for the direct ARP → aldehyde alternative). **[NEG] No SS, no n, no per-response variance, no fructose predicted-vs-observed panel** (fructose is absent from Fig. 4 although it is a fitted response). The sugar and amino-acid fits are poor by the paper's own admission — "sugars and amino acids had the least good fit" — because raw green malt from industrial germination boxes is heterogeneous (p.5). **The precursor pools are effectively unconstrained by the data; the fit is carried by ARP and aldehyde curves.** Fructose panel absence + glucose R² = 0.22 means the k₁-type constants are anchored by ARP formation curves, not by glucose depletion.

Also noted (p.7–9): at 90 °C aldehydes fell between the last two sampling points in some runs — attributed to "evaporation and/or thermal degradation" — and no loss term was fitted "since it would have led to overfitting". Two-way ANOVA: no significant variety effect on any aldehyde (p > 0.05, Table S3).

---

## §7. WHAT THE REPOSITORY CAN USE

Trunk lane = mass-action network (Amadori compounds, 1-/3-deoxyosones, methylglyoxal as species; fitted (k at 100 °C, Ea) pairs). Sulfur lane = Strecker degradation of cysteine (H₂S release), free-cysteine depletion Ea = 55.1 kJ/mol (`kang2026_SI_extraction.md` §5b, provenance re-pointed to `zhai2023foodchem` per `k6a_sulfur_ladders_synthesis.md`).

| this paper's parameter | repo step it could inform | role | verdict |
|---|---|---|---|
| **k′₂ = 7.30·10⁻⁴ h⁻¹ at 70 °C, Ea₂ = 315 (fixed)** — ARP → dicarbonyl pool + AA | trunk `k_ama_tdg` (DFG → 3-DG, Martins Ea 97.1 ± 1.7), `k_ama_odg` (→ 1-DG, 107.3 ± 7.3), `k_ama_mgo` (→ MG, 124.5 ± 4.7) | **validation only, ordinal.** Piornos's total ARP-cleavage rate at 100 °C extrapolates to 5.2 h⁻¹ = 0.087 min⁻¹ `[Z]` vs Martins's summed DFG sinks 1.1·10⁻² + 1.6·10⁻² + 7.1·10⁻³ = 3.4·10⁻² min⁻¹ — **2.6× faster, from a 315 kJ/mol extrapolation over 30 K with no interval.** The Ea itself (315) is 2.5–3× the aqueous values and must not replace them. | ordinal prior: "ARP cleavage is the most T-sensitive trunk step" — consistent in sign with Martins (ARP-sink Ea 97–125 > 3-DG → acid 29.6) |
| **k′_M = 5.19·10⁻⁴ kg mmol⁻¹ h⁻¹ (92 % error), Ea_M = 255 (fixed)** — dicarbonyl pool + AA → Strecker aldehyde | trunk Strecker consumption of MG / 3-DG by amino acids; sulfur lane cysteine Strecker step | **do not use numerically.** 92 % HPD; Int1 unmeasured; Int1 lumps 3-DG, 1-DG, MG, glyoxal, diacetyl. At 100 °C `[Z]` k_M = 0.685 kg mmol⁻¹ h⁻¹ = 1.1·10⁻² L mmol⁻¹ min⁻¹ if 1 kg ≈ 1 L — same order as Martins's bimolecular 3-DG + Gly → Mel (8.1·10⁻⁴) only within a decade. | Ea_M = 255 vs sulfur-lane 55.1 kJ/mol for cysteine depletion: **a 4.6× discrepancy**; the two are not the same step (cysteine depletion is the aqueous free-AA sink; Ea_M is a fixed, solid-matrix, unobserved-pool value). Piornos gives **no support** for raising the 55.1; it shows only that solid-matrix apparent Ea for Strecker steps can be inflated by moisture loss |
| **F_Met = 0.269 ± 0.02** (share of Met consumed by Int1 that becomes methional; rest → MRP) | sulfur lane: branching of methionine Strecker vs other sinks | usable as a **branching prior** for methionine at 65–90 °C in a solid matrix (methional yield per Met consumed ≈ 27 %); F_Leu = 1.0 (bound), F_Ile 0.52, F_Val 0.17, F_Phe 0.21 for the non-sulfur Strecker aldehydes | **USE (prior, wide)** — but F is confounded with the shared k_M and with volatile losses at 90 °C |
| **Ea₁ = 175 (fixed)**, k′₂₁…k′_i1 at 70 °C | trunk `k_schiff` (Glu + Gly → E1, Martins 96.8 ± 2.8) | k₂₁ at 100 °C `[Z]` 1.46·10⁻² kg mmol⁻¹ h⁻¹ = 2.4·10⁻⁴ per mmol·min vs Martins 1.6·10⁻⁵ L mmol⁻¹ min⁻¹ — 15× faster, in a matrix where water is scarce (condensation favoured). Ea 175 vs 97: again ~1.8× | validation of *ordering* only: Leu > Ile > Phe > Val > (Met, pooled) for ARP formation rate — a per-amino-acid reactivity ladder the repo lacks |
| **Ea_F = 121, k′_F = 1.82·10⁻² h⁻¹** — fructose → dicarbonyl pool in one step | trunk fructose branch (`k_fru_glc` 93.4; Martins has no direct Fru → dicarbonyl step) | the repo's trunk routes fructose via isomerisation; Piornos found a direct one-step Fru → Int1 improved SS (Kato 1969, Mundt & Wedzicha 2003 precedent). Ea_F = 121 is the **least inflated** of the four and closest to Martins's `k_glc_fru` 122.6 ± 5.2 | **candidate structural prior** (a direct fructose → 1-DG/MG-type sink) rather than a numeric one |
| Ea ordering Ea₂ (315) > Ea_M (255) > Ea₁ (175) > Ea_F (121) | whole trunk | the one robust qualitative result: **downstream (intermediate-consuming) steps are steeper than upstream (precursor-consuming) steps**, so aldehyde yield per unit precursor loss rises with T | USE as an ordinal sanity check on the trunk's Ea set |

**Transfer limits.** (1) **Window 65–90 °C** — every repository trunk k is at 100 °C; using Piornos's k′ requires extrapolating a fixed 175–315 kJ/mol slope 10–30 K beyond the data, which is where a moisture-confounded Ea does the most damage (k₂ changes 2240× across 25 K in-window; 16× more over the next 10 K). (2) **Matrix** — a drying cereal bed with unreported moisture, unreported pH, unreported bed temperature; units are per kg malt, not per litre; bimolecular k′ cannot be converted to L mmol⁻¹ without a water content that is not given. (3) **Species** — no measured dicarbonyl, no cysteine, no H₂S, no thiol; "Int1" is not any one of the repo's species. (4) **No Ea intervals** — the paper cannot supply a prior width; any prior built on these Ea must be declared as "point value, fixed by the authors, unknown width".

**Net recommendation:** cite as an **ordinal/structural** source (Ea ranking; direct fructose sink; per-amino-acid ARP-formation ladder; F_Met ≈ 0.27 methional branching); do **not** ingest any Ea numerically into the trunk or sulfur lane; the k′ values may serve a 65–90 °C solid-matrix *hold-out* if the repo ever runs a low-moisture scenario, provided the moisture confound is flagged.

---

## §8. CAVEATS

1. **All four activation energies are fixed point values with no interval.** The paper's alternating k′/Ea runs are not reported. Any use of 175/315/255/121 must carry the flag "author-fixed, no HPD".
2. **Moisture confound.** Curing is a drying step; the bed's water content falls with time and is lower at higher T. A fixed Arrhenius slope fitted across 65/78/90 °C therefore absorbs the a_w trend; Q₁₀ ≈ 25 (k₂) and ≈ 13.5 (k_M) `[Z]` are not chemical Q₁₀ values. Main text gives **no moisture during curing, no a_w, no pH**; Table S1 (pre-kiln moisture) is not on disk.
3. **Int1 is unmeasured** — a modelling construct lumping deoxyosones and C2–C4 dicarbonyls; k′_M has a 92 % HPD, and Ea_M is the temperature dependence of a species nobody observed. FruMet, AAj and FruAAj are also model-estimated, not measured.
4. **Precursor fits are poor** (glucose R² = 0.22; amino acids 0.37–0.45; slopes 0.30–0.49). The paper attributes this to raw-material heterogeneity and kiln non-uniformity. The k₁-type constants are constrained by ARP curves, not by precursor depletion; there is **no fructose predicted-vs-observed panel** although fructose is a fitted response.
5. **Shared constants by design**: one k₂ for seven ARPs, one k_M for all Int1 + AA reactions, one k_i1 for Met + AAi + AAj. Individual constants "gained in quality" only when merged — i.e. the data cannot resolve per-ARP degradation rates.
6. **Minor internal inconsistency**: text says methional R² = 0.81; Fig. 4o prints 0.8507. HPD half-widths are printed to one significant figure and disagree with the % column by up to 2 points (k′₄₁: 3.9 % vs 6 %).
7. **Volatile losses at 90 °C** in the last interval (evaporation / degradation) were deliberately not modelled; the F factors and k_M therefore carry a small downward bias at the top temperature.
8. **No loss/degradation term for aldehydes, no reversibility, no direct glucose → dicarbonyl route** — the network is minimal by construction and rejects the ARP → aldehyde direct route only on AIC (639 vs 672), not on mechanism.
9. **Sulfur relevance is thin**: methionine only (k_i1 pooled; FruMet unmeasured; F_Met = 0.269 ± 0.02); no cysteine, H₂S, thiols or sulfur volatiles. The Ea_M = 255 kJ/mol has no bearing on the sulfur lane's 55.1 kJ/mol cysteine depletion except as a warning that solid-matrix apparent Ea are inflated.
10. **Funding/COI**: fully Heineken-funded, four industry co-authors — no reason to doubt the numbers, but the modelling goal (kiln control for alcohol-free beer) drove the choice of window and species.
11. **Supplementary data absent from the repo** (Tables S1–S3, Figs S1–S2): the full 12-run dataset, the Int1 prediction and the moisture table cannot be checked here.
