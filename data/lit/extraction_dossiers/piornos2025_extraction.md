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

---

## Supplementary material (read 2026-10-29)

Source: `data/articles/piornos2025_supplementary.docx` (owner's download). Tables S1–S3 read with `textutil`. **Figs. S1 and S2
are Word-native charts whose data are stored as numeric caches inside the .docx; the values below were read from those caches
(exact as stored — not digitised from an image)**, by unzipping the chart XML to a scratch directory and parsing it. This closes
caveat 11 (§8) and the [NEG] of §0 ("NOT on disk").

### Table S1 — *"Moisture content of the green malts used for the kilning experiments."* (pre-kiln, %w, mean ± SD)

| sample | variety | replicate | curing T (°C) | moisture (%w) |
|---|---|---|---|---|
| 65S0 | S | 1 | 65 | 38.6 ± 2.73 |
| 65W0 | W | 1 | 65 | 39.7 ± 4.60 |
| 78S0.1 | S | 1 | 78 | 44.0 ± 1.29 |
| 78S0.2 | S | 2 | 78 | 38.3 ± 0.08 |
| 78W0.1 | W | 1 | 78 | 41.6 ± 0.98 |
| 78W0.2 | W | 2 | 78 | 40.5 ± 0.16 |
| 90S0.1 | S | 1 | 90 | 40.3 ± 0.03 |
| 90S0.2 | S | 2 | 90 | 40.9 ± 0.06 |
| 90W0.1 | W | 1 | 90 | 42.3 ± 0.26 |
| 90W0.2 | W | 2 | 90 | 41.4 ± 0.11 |

Footnote: *"Moisture expressed as mean value ± standard deviation. Varieties of barley: S for spring (Etincel), W for winter
(Planet)."* ⚠️ **This contradicts the main text (§2.1), which calls "RGT Planet" the two-row spring variety and "Etincel" the
six-row winter variety.** The S/W labels of the data are therefore ambiguous as to cultivar; the ANOVA finds no variety effect,
so nothing numeric depends on it. Only one moisture per 65 °C condition (the 65 °C duplicates came from one batch, §2.2).
**Still no moisture during curing** — Table S1 is the green malt before the 16 h kiln programme (≈ 38–44 %w).

### Table S2 — MRM conditions (LC-MS/MS) for 19 amino acids + norvaline (IS), glucose, fructose, trehalose and seven ARPs
(FruVal, FruIle, FruLeu, FruPhe, FruAla, FruGly, FruPro). Method parameters only; no LOD/LOQ.

### Table S3 — two-way ANOVA (temperature × variety) on the five Strecker aldehydes, α = 0.05

Temperature significant for all five (F = 99.8 3MB, 84.8 2MB, 94.0 2MP, 209.7 PhAc, 25.3 methional; p ≤ 0.001); variety not
significant (p = 0.243, 0.625, 0.473, 0.103, 0.152); interaction not significant (p = 0.107–0.421). Error DF 6, total DF 12
(i.e. one value per run — presumably the 8.4 h endpoint; not stated). Full sums of squares are in the source table.

### Fig. S1 — *"Concentration of amino acids (a-g), sugars (h-i), ARP (j-n) and Strecker aldehydes (o-s) as a function of time
during the curing stage of kilning at different temperatures (65, 78 or 90 °C) for two different varieties of barley (S for spring,
W for winter) and two replicates."*

Each chart stores an OBSERVED and a PREDICTED (model) series; only OBSERVED is transcribed. Units **mmol/kg** (kilned malt; fresh
or dry basis not stated, §3). 12 runs × 8 times (curing time 0–8.4 h, every 1.2 h). Single values per run (no SD). "AAi" and "AAj"
are the paper's pooled amino-acid groups; "FruAAi" the pooled ARPs of AAi. Values shown to 4 significant figures.

⚠️ **Data-integrity findings in the stored data `[D]`:**
1. **65 °C Strecker aldehydes: replicate 2 is a copy of replicate 1.** For all five aldehydes, row 65S2 is identical to 65S1 and
   65W2 to 65W1 in every one of the 8 points. Amino acids, sugars and ARPs at 65 °C differ between replicates. Either the aldehyde
   duplicates were never measured separately or the data were copied; the 65 °C aldehyde "replicates" carry no independent information.
2. **Amino acids (all seven panels a–g), run 90W2: the 7.2 h and 8.4 h values repeat the 4.8 h and 6.0 h values exactly**
   (e.g. Leu 5.613, 5.327, 5.613, 5.327). The last two 90W2 amino-acid points are not measurements.
3. **Glucose and fructose are erratic and bimodal** within runs (e.g. glucose 65S2: 14.1, 34.9, 41.6, 40.2, 12.9, 44.6, 26.0,
   48.7 mmol/kg; values cluster near ~9–14 and ~30–48 regardless of time or temperature). The ARPs and aldehydes are smooth. Sugar
   loss cannot be read from these series; the fit's glucose and fructose constants rest on noise of ±50 %.

**Panel (chart 1): Leucine, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 6.498 | 7.513 | 8.275 | 7.306 | 7.221 | 7.035 | 7.39 | 6.489 |
| 65S2 | 6.595 | 7.044 | 6.985 | 7.487 | 6.708 | 6.218 | 5.379 | 5.615 |
| 65W1 | 5.809 | 5.17 | 5.753 | 5.185 | 5.61 | 4.913 | 5.931 | 5.814 |
| 65W2 | 5.14 | 5.961 | 4.94 | 5.803 | 5.254 | 5.417 | 7.317 | 7.814 |
| 78S1 | 6.888 | 6.378 | 7.176 | 6.43 | 7.595 | 7.048 | 8.398 | 8.753 |
| 78S2 | 7.371 | 8.3 | 8.536 | 7.253 | 8.193 | 8.106 | 7.13 | 6.464 |
| 78W1 | 6.805 | 6.611 | 6.568 | 7.026 | 6.206 | 5.927 | 6.656 | 6.586 |
| 78W2 | 7.434 | 7.598 | 7.282 | 6.343 | 5.983 | 7.661 | 8.232 | 5.831 |
| 90S1 | 6.75 | 7.02 | 6.534 | 6.314 | 6.3 | 6.674 | 9.245 | 7.386 |
| 90S2 | 7.686 | 7.687 | 7.242 | 7.134 | 6.932 | 6.148 | 6.798 | 6.23 |
| 90W1 | 6.147 | 5.795 | 5.902 | 5.689 | 4.878 | 5.454 | 7.638 | 6.551 |
| 90W2 | 6.476 | 5.763 | 5.608 | 5.897 | 5.613 | 5.327 | 5.613 | 5.327 |

**Panel (chart 2): Isoleucine, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 6.451 | 7.518 | 8.706 | 7.294 | 7.312 | 7.204 | 7.606 | 6.551 |
| 65S2 | 6.758 | 7.126 | 7.222 | 7.486 | 7.048 | 6.483 | 5.467 | 5.976 |
| 65W1 | 5.627 | 5.289 | 5.831 | 5.185 | 5.564 | 5 | 6.13 | 5.812 |
| 65W2 | 5.148 | 5.981 | 5.126 | 5.849 | 5.386 | 5.463 | 7.411 | 7.958 |
| 78S1 | 6.789 | 6.326 | 7.365 | 6.979 | 7.604 | 7.086 | 8.341 | 8.708 |
| 78S2 | 7.541 | 8.219 | 8.506 | 7.173 | 8.125 | 8.135 | 6.96 | 6.358 |
| 78W1 | 6.732 | 6.995 | 6.912 | 7.003 | 6.39 | 6.343 | 6.734 | 6.595 |
| 78W2 | 7.447 | 7.335 | 7.553 | 6.144 | 5.937 | 7.636 | 8.51 | 5.901 |
| 90S1 | 7.139 | 7.04 | 6.726 | 6.569 | 6.398 | 6.476 | 9.171 | 7.451 |
| 90S2 | 7.498 | 8.197 | 7.481 | 7.307 | 7.045 | 6.304 | 6.902 | 6.174 |
| 90W1 | 6.187 | 6.104 | 6.059 | 5.828 | 5.002 | 5.634 | 7.473 | 6.331 |
| 90W2 | 6.271 | 5.758 | 5.643 | 5.788 | 5.777 | 5.366 | 5.777 | 5.366 |

**Panel (chart 3): Valine, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 18.01 | 19.23 | 21.17 | 18.83 | 17.46 | 17.32 | 25.46 | 17.64 |
| 65S2 | 19.59 | 18.28 | 22.88 | 18.16 | 18.31 | 18.29 | 13.79 | 13.56 |
| 65W1 | 13.27 | 13.54 | 14.49 | 16.03 | 13.78 | 13.35 | 14.26 | 14.37 |
| 65W2 | 13.89 | 14.35 | 15.91 | 16.14 | 12.88 | 13.26 | 22.27 | 19.88 |
| 78S1 | 18.41 | 20.28 | 16.95 | 17.85 | 21.84 | 15.38 | 22.41 | 21.65 |
| 78S2 | 19.25 | 20.52 | 21.52 | 18.59 | 19.41 | 21.61 | 25.16 | 19.87 |
| 78W1 | 23.34 | 18.89 | 16.97 | 19.27 | 20.22 | 18.38 | 19.88 | 19.11 |
| 78W2 | 18.86 | 18.6 | 18.73 | 16.23 | 15.51 | 26.31 | 20.77 | 15.16 |
| 90S1 | 16.92 | 16.37 | 18.58 | 16 | 14.67 | 15.4 | 24.1 | 19.83 |
| 90S2 | 19.72 | 21.68 | 23.04 | 22.18 | 18.34 | 18.89 | 19.19 | 15.5 |
| 90W1 | 18.96 | 15.02 | 14.28 | 13.49 | 14.23 | 16.88 | 19.5 | 15.38 |
| 90W2 | 15.26 | 15.29 | 13.67 | 14.71 | 14.83 | 12.12 | 14.83 | 12.12 |

**Panel (chart 4): Phenylalanine, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 12.78 | 14.42 | 17.82 | 14.23 | 13.58 | 15.05 | 16.21 | 12.59 |
| 65S2 | 13.92 | 13.78 | 15.26 | 15.7 | 14.16 | 13.92 | 11.35 | 13.27 |
| 65W1 | 13.04 | 11.07 | 13 | 12.15 | 12.22 | 11.75 | 13.78 | 12.1 |
| 65W2 | 11.13 | 14.22 | 11.63 | 13.78 | 11.1 | 11.22 | 15.63 | 14.77 |
| 78S1 | 13.49 | 13.25 | 13.58 | 13.44 | 15.2 | 12.76 | 17.06 | 17.21 |
| 78S2 | 14.98 | 16.52 | 17.12 | 14.24 | 17.62 | 16.69 | 20.12 | 15.63 |
| 78W1 | 17.7 | 16.47 | 15.68 | 18.43 | 15.04 | 14.54 | 16.19 | 16.03 |
| 78W2 | 17.3 | 19.25 | 19.43 | 14.11 | 13.66 | 19.62 | 18.18 | 11.26 |
| 90S1 | 13.81 | 13.05 | 12.69 | 12.51 | 11.85 | 11.68 | 19.75 | 15.94 |
| 90S2 | 15.48 | 18.38 | 15.09 | 14.65 | 13.06 | 13 | 17.89 | 13.96 |
| 90W1 | 14.61 | 13.11 | 12.83 | 12.83 | 11.06 | 12.31 | 17.46 | 14.84 |
| 90W2 | 13.91 | 12.05 | 11.69 | 13.32 | 12.57 | 10.59 | 12.57 | 10.59 |

**Panel (chart 5): Methionine, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 6.481 | 7.107 | 7.377 | 6.713 | 6.499 | 6.669 | 8.131 | 6.373 |
| 65S2 | 6.896 | 6.482 | 7.615 | 6.822 | 6.893 | 6.558 | 5.49 | 5.759 |
| 65W1 | 5.823 | 5.446 | 5.93 | 5.977 | 5.769 | 5.333 | 6.021 | 5.713 |
| 65W2 | 5.38 | 5.861 | 5.819 | 6.426 | 5.293 | 5.298 | 7.449 | 6.933 |
| 78S1 | 6.705 | 6.651 | 6.765 | 6.619 | 7.426 | 6.122 | 7.157 | 7.385 |
| 78S2 | 6.721 | 6.941 | 7.113 | 6.455 | 6.77 | 7.122 | 7.105 | 6.24 |
| 78W1 | 6.987 | 6.435 | 6.144 | 6.722 | 6.44 | 6.029 | 6.15 | 5.92 |
| 78W2 | 6.289 | 6.346 | 6.429 | 5.543 | 5.582 | 7.485 | 7.544 | 5.992 |
| 90S1 | 6.457 | 6.054 | 6.501 | 6.107 | 5.94 | 5.755 | 7.556 | 6.613 |
| 90S2 | 6.525 | 7.044 | 6.979 | 6.66 | 5.88 | 6.139 | 6.416 | 5.721 |
| 90W1 | 6.399 | 5.615 | 5.422 | 5.354 | 5.36 | 5.857 | 6.597 | 5.777 |
| 90W2 | 5.947 | 5.726 | 5.342 | 5.597 | 5.876 | 5.349 | 5.876 | 5.349 |

**Panel (chart 6): AAi, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 159.1 | 190.4 | 218.3 | 196 | 177.8 | 171.1 | 149.3 | 158.9 |
| 65S2 | 156.1 | 179.5 | 159.9 | 192.6 | 168.6 | 147.8 | 116 | 119.4 |
| 65W1 | 116.8 | 116.4 | 118.9 | 110.9 | 118.2 | 104.7 | 121.9 | 123.4 |
| 65W2 | 115.7 | 124.4 | 108.4 | 134.8 | 116.2 | 116 | 167.5 | 180.4 |
| 78S1 | 147.1 | 143.1 | 163.2 | 152.7 | 162.2 | 153.9 | 180 | 184.7 |
| 78S2 | 153.4 | 172.4 | 180.7 | 157.2 | 180.8 | 168.4 | 137.3 | 133.9 |
| 78W1 | 151.1 | 149.5 | 156 | 163.3 | 140.2 | 132.8 | 143.9 | 144 |
| 78W2 | 160.1 | 154.5 | 156.4 | 133.9 | 129.2 | 139.2 | 192.6 | 130.8 |
| 90S1 | 172.2 | 175.9 | 145.5 | 167.4 | 156.8 | 164.8 | 200.6 | 172.3 |
| 90S2 | 172.5 | 184.9 | 126.3 | 136.1 | 150.5 | 120.4 | 145.3 | 138.5 |
| 90W1 | 121.2 | 137.8 | 136.6 | 129.8 | 99.58 | 112.8 | 165.5 | 142.2 |
| 90W2 | 142 | 118.4 | 133.4 | 133 | 136.1 | 118.8 | 136.1 | 118.8 |

**Panel (chart 7): AAj, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 116.9 | 105 | 134.5 | 98.44 | 97.68 | 93 | 141.7 | 114.6 |
| 65S2 | 117.8 | 96.91 | 134.2 | 99.57 | 111.8 | 109.4 | 98.75 | 90.62 |
| 65W1 | 89.97 | 99.15 | 101.1 | 113.9 | 97.18 | 96.67 | 92.26 | 92.36 |
| 65W2 | 102.9 | 92.01 | 115 | 104.1 | 91.77 | 89.43 | 134.3 | 103.8 |
| 78S1 | 114.8 | 120.5 | 102.9 | 103.8 | 112.8 | 90.59 | 126.3 | 111.4 |
| 78S2 | 114.1 | 103.5 | 107.1 | 111.7 | 96.1 | 120.3 | 130.2 | 123.5 |
| 78W1 | 130.6 | 115.3 | 99.05 | 102.2 | 124.4 | 116 | 124.8 | 115 |
| 78W2 | 97.28 | 98.3 | 98.27 | 98.44 | 95.4 | 130.2 | 104.2 | 103.2 |
| 90S1 | 95.84 | 90.64 | 110.4 | 89.69 | 87.56 | 84.39 | 118.2 | 111.7 |
| 90S2 | 112.1 | 98.76 | 116.3 | 117 | 102.8 | 105.6 | 117.8 | 101.8 |
| 90W1 | 121.6 | 90.27 | 85.14 | 87.81 | 103.8 | 115.3 | 105.9 | 90.49 |
| 90W2 | 93.68 | 112 | 90.28 | 89.23 | 94.6 | 90.03 | 94.6 | 90.03 |

**Panel (chart 8): Glucose, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 44.13 | 43.3 | 23.6 | 37.39 | 40.49 | 37.27 | 33.99 | 41.68 |
| 65S2 | 14.05 | 34.86 | 41.63 | 40.19 | 12.89 | 44.57 | 26 | 48.71 |
| 65W1 | 27.66 | 40.39 | 43.2 | 30.92 | 12.63 | 41.6 | 9.134 | 41.56 |
| 65W2 | 41.72 | 30.15 | 27.21 | 32.28 | 41.76 | 9.755 | 28.92 | 29.49 |
| 78S1 | 44.92 | 37.17 | 14.31 | 43.78 | 24.03 | 12.39 | 10.93 | 14.43 |
| 78S2 | 33.45 | 36.02 | 24.9 | 38.15 | 36.13 | 30.29 | 39.98 | 9.16 |
| 78W1 | 39.96 | 43.09 | 37.67 | 5.884 | 9.677 | 9.599 | 43.08 | 9.541 |
| 78W2 | 20.66 | 30.53 | 30.29 | 39.52 | 41.84 | 36.39 | 37.26 | 8.934 |
| 90S1 | 45.68 | 16.9 | 12.61 | 33.21 | 13.21 | 12.18 | 12.2 | 39.09 |
| 90S2 | 35.31 | 33.17 | 31.25 | 11.17 | 11.04 | 9.931 | 9.195 | 35.53 |
| 90W1 | 46.87 | 9.722 | 10.98 | 29.19 | 28.71 | 11.14 | 35.05 | 9.479 |
| 90W2 | 36.26 | 39.43 | 29.97 | 11.67 | 10.24 | 8.736 | 8.478 | 10.82 |

**Panel (chart 9): Fructose, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 15.7 | 15.91 | 10.2 | 14.63 | 16.97 | 14.11 | 11.22 | 17.17 |
| 65S2 | 9.128 | 14.66 | 16.22 | 15.49 | 9.111 | 16.44 | 11.75 | 17.97 |
| 65W1 | 12.84 | 20 | 20.38 | 14.24 | 9.009 | 22.17 | 9.24 | 18.86 |
| 65W2 | 20.67 | 13.27 | 12.4 | 22.97 | 20.73 | 9.451 | 13.49 | 11.37 |
| 78S1 | 20.29 | 14.74 | 7.089 | 20.17 | 10.02 | 6.517 | 9.042 | 7.721 |
| 78S2 | 15.02 | 15.35 | 11.04 | 17.11 | 15.62 | 13.37 | 17.15 | 12.45 |
| 78W1 | 22.13 | 26.36 | 26 | 13.45 | 7.972 | 10.95 | 26.62 | 8.285 |
| 78W2 | 14.13 | 17.4 | 14.65 | 22.31 | 24.64 | 16.28 | 18.54 | 13.18 |
| 90S1 | 18.72 | 8.748 | 8.615 | 11.68 | 10.09 | 8.82 | 8.136 | 16.83 |
| 90S2 | 14.61 | 13.66 | 13.36 | 9.522 | 9.734 | 9.329 | 6.345 | 16.27 |
| 90W1 | 22.53 | 9.336 | 9.882 | 13.65 | 14.03 | 7.92 | 18.53 | 10.18 |
| 90W2 | 20.23 | 19.24 | 13.95 | 9.196 | 8.904 | 9.858 | 10.67 | 8 |

**Panel (chart 10): FruLeu, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.1406 | 0.1211 | 0.1216 | 0.1515 | 0.1727 | 0.2234 | 0.1563 | 0.1815 |
| 65S2 | 0.1431 | 0.1235 | 0.1447 | 0.1329 | 0.1604 | 0.2015 | 0.1816 | 0.1512 |
| 65W1 | 0.0559 | 0.0857 | 0.0433 | 0.0616 | 0.0682 | 0.0798 | 0.085 | 0.0865 |
| 65W2 | 0.0776 | 0.0535 | 0.0815 | 0.1068 | 0.0695 | 0.0742 | 0.0698 | 0.0736 |
| 78S1 | 0.1795 | 0.2365 | 0.4205 | 0.3527 | 0.3744 | 0.4723 | 0.4719 | 0.4504 |
| 78S2 | 0.1728 | 0.3864 | 0.3328 | 0.5737 | 0.6524 | 0.8189 | 0.7265 | 0.77 |
| 78W1 | 0.2375 | 0.3591 | 0.5649 | 0.554 | 0.7399 | 0.6012 | 0.7702 | 0.7332 |
| 78W2 | 0.2712 | 0.5596 | 0.621 | 0.7971 | 1.007 | 0.675 | 0.9353 | 1.17 |
| 90S1 | 0.1818 | 0.5121 | 0.6423 | 0.6267 | 0.7816 | 0.6848 | 0.6914 | 0.4668 |
| 90S2 | 0.308 | 0.6196 | 0.6134 | 0.8611 | 0.9612 | 0.9165 | 0.788 | 0.6554 |
| 90W1 | 0.1849 | 0.4512 | 0.6104 | 0.4707 | 0.5084 | 0.5043 | 0.4248 | 0.497 |
| 90W2 | 0.2647 | 0.547 | 0.6987 | 0.755 | 0.654 | 0.5924 | 0.6267 | 0.585 |

**Panel (chart 11): FruIle, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.1003 | 0.0839 | 0.0803 | 0.1115 | 0.1245 | 0.1593 | 0.0989 | 0.1271 |
| 65S2 | 0.1065 | 0.0829 | 0.1041 | 0.0927 | 0.1069 | 0.1435 | 0.1258 | 0.1053 |
| 65W1 | 0.0459 | 0.0625 | 0.032 | 0.041 | 0.0516 | 0.0599 | 0.0578 | 0.0623 |
| 65W2 | 0.0588 | 0.0398 | 0.0519 | 0.0781 | 0.0535 | 0.056 | 0.0454 | 0.0608 |
| 78S1 | 0.1299 | 0.1605 | 0.2662 | 0.1853 | 0.2379 | 0.3076 | 0.2854 | 0.2604 |
| 78S2 | 0.1445 | 0.2517 | 0.2678 | 0.3636 | 0.4196 | 0.4576 | 0.4608 | 0.4593 |
| 78W1 | 0.1673 | 0.244 | 0.3851 | 0.3362 | 0.478 | 0.3767 | 0.5087 | 0.4208 |
| 78W2 | 0.1792 | 0.3511 | 0.4101 | 0.5194 | 0.6597 | 0.4059 | 0.5473 | 0.6985 |
| 90S1 | 0.1304 | 0.3018 | 0.3915 | 0.3626 | 0.4675 | 0.3874 | 0.3874 | 0.2715 |
| 90S2 | 0.2215 | 0.3853 | 0.3748 | 0.5132 | 0.5956 | 0.5778 | 0.4645 | 0.3936 |
| 90W1 | 0.1372 | 0.2899 | 0.375 | 0.2915 | 0.3185 | 0.3048 | 0.2648 | 0.3194 |
| 90W2 | 0.1818 | 0.3389 | 0.44 | 0.446 | 0.4061 | 0.3624 | 0.3826 | 0.3511 |

**Panel (chart 12): FruVal, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.111 | 0.0956 | 0.1133 | 0.1234 | 0.1368 | 0.1769 | 0.1428 | 0.1424 |
| 65S2 | 0.1173 | 0.1166 | 0.1148 | 0.1041 | 0.1283 | 0.1626 | 0.1448 | 0.1164 |
| 65W1 | 0.0694 | 0.0807 | 0.0444 | 0.0642 | 0.0663 | 0.0777 | 0.077 | 0.0776 |
| 65W2 | 0.0747 | 0.0542 | 0.0801 | 0.0963 | 0.0692 | 0.0717 | 0.0757 | 0.0977 |
| 78S1 | 0.1382 | 0.1789 | 0.3008 | 0.2483 | 0.2761 | 0.3513 | 0.3544 | 0.3756 |
| 78S2 | 0.1807 | 0.2659 | 0.342 | 0.3709 | 0.4211 | 0.5996 | 0.4653 | 0.5263 |
| 78W1 | 0.1771 | 0.2581 | 0.3907 | 0.4034 | 0.4915 | 0.3904 | 0.5245 | 0.5119 |
| 78W2 | 0.2433 | 0.4634 | 0.4433 | 0.5299 | 0.6611 | 0.5485 | 0.7251 | 0.7429 |
| 90S1 | 0.1366 | 0.402 | 0.4361 | 0.4991 | 0.5087 | 0.4472 | 0.436 | 0.3233 |
| 90S2 | 0.2191 | 0.4895 | 0.4973 | 0.5977 | 0.6 | 0.5729 | 0.5757 | 0.4282 |
| 90W1 | 0.1453 | 0.3126 | 0.4442 | 0.3224 | 0.3502 | 0.3331 | 0.3007 | 0.3417 |
| 90W2 | 0.1964 | 0.3614 | 0.4403 | 0.5487 | 0.4303 | 0.3798 | 0.3827 | 0.3784 |

**Panel (chart 13): FruPhe, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.1475 | 0.1188 | 0.126 | 0.1617 | 0.1869 | 0.2328 | 0.1536 | 0.1949 |
| 65S2 | 0.1548 | 0.1212 | 0.1518 | 0.1438 | 0.1203 | 0.2129 | 0.1866 | 0.1553 |
| 65W1 | 0.0952 | 0.111 | 0.0575 | 0.0745 | 0.0944 | 0.106 | 0.1052 | 0.1096 |
| 65W2 | 0.1012 | 0.0659 | 0.0991 | 0.14 | 0.0908 | 0.1041 | 0.0845 | 0.1058 |
| 78S1 | 0.1776 | 0.225 | 0.4006 | 0.3321 | 0.3539 | 0.4476 | 0.4317 | 0.3981 |
| 78S2 | 0.2149 | 0.3642 | 0.405 | 0.5398 | 0.6239 | 0.6868 | 0.7058 | 0.7663 |
| 78W1 | 0.2751 | 0.4146 | 0.6613 | 0.6611 | 0.8357 | 0.6837 | 0.8703 | 0.9866 |
| 78W2 | 0.2966 | 0.5869 | 0.7062 | 0.9187 | 1.136 | 0.6933 | 1.01 | 1.253 |
| 90S1 | 0.195 | 0.4635 | 0.6132 | 0.5745 | 0.77 | 0.6172 | 0.6466 | 0.4234 |
| 90S2 | 0.3191 | 0.593 | 0.5911 | 0.7824 | 0.953 | 0.9204 | 0.888 | 0.6248 |
| 90W1 | 0.2351 | 0.567 | 0.8126 | 0.538 | 0.6059 | 0.6058 | 0.5048 | 0.6233 |
| 90W2 | 0.3086 | 0.6209 | 0.8368 | 1.031 | 0.7599 | 0.6922 | 0.7456 | 0.7092 |

**Panel (chart 14): FruAAi, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.2845 | 0.2427 | 0.2359 | 0.3492 | 0.3676 | 0.4925 | 0.3234 | 0.4137 |
| 65S2 | 0.3214 | 0.2551 | 0.298 | 0.2898 | 0.3526 | 0.4584 | 0.4126 | 0.3206 |
| 65W1 | 0.1859 | 0.2938 | 0.1724 | 0.1914 | 0.2322 | 0.2965 | 0.2975 | 0.2714 |
| 65W2 | 0.2756 | 0.2099 | 0.228 | 0.343 | 0.2486 | 0.2885 | 0.2366 | 0.2882 |
| 78S1 | 0.3435 | 0.4963 | 0.942 | 0.7489 | 0.9612 | 1.267 | 1.175 | 1.095 |
| 78S2 | 0.331 | 0.6948 | 0.6894 | 1.035 | 1.156 | 1.588 | 1.37 | 1.434 |
| 78W1 | 0.5121 | 0.8344 | 1.286 | 1.209 | 1.738 | 1.429 | 1.929 | 1.637 |
| 78W2 | 0.6178 | 1.283 | 1.463 | 1.691 | 2.088 | 1.824 | 2.295 | 2.587 |
| 90S1 | 0.3145 | 1.116 | 1.595 | 1.48 | 2.102 | 1.756 | 1.938 | 1.332 |
| 90S2 | 0.506 | 1.263 | 1.314 | 1.662 | 2.201 | 2.158 | 1.849 | 1.692 |
| 90W1 | 0.3813 | 1.033 | 1.481 | 1.194 | 1.283 | 1.285 | 1.245 | 1.48 |
| 90W2 | 0.5779 | 1.196 | 1.613 | 1.722 | 1.541 | 1.655 | 1.738 | 1.527 |

**Panel (chart 15): 3-Methylbutanal, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.00294 | 0.004 | 0.0036 | 0.00402 | 0.00555 | 0.00799 | 0.00724 | 0.00713 |
| 65S2 | 0.00294 | 0.004 | 0.0036 | 0.00402 | 0.00555 | 0.00799 | 0.00724 | 0.00713 |
| 65W1 | 0.00322 | 0.00315 | 0.00308 | 0.0031 | 0.00303 | 0.00529 | 0.00512 | 0.0061 |
| 65W2 | 0.00322 | 0.00315 | 0.00308 | 0.0031 | 0.00303 | 0.00529 | 0.00512 | 0.0061 |
| 78S1 | 0.00361 | 0.01219 | 0.02428 | 0.02626 | 0.05846 | 0.06401 | 0.1197 | 0.08976 |
| 78S2 | 0.00636 | 0.03456 | 0.03022 | 0.05105 | 0.05751 | 0.1115 | 0.1059 | 0.1545 |
| 78W1 | 0.01214 | 0.02644 | 0.03882 | 0.07179 | 0.1188 | 0.1132 | 0.1729 | 0.1686 |
| 78W2 | 0.0156 | 0.03801 | 0.06882 | 0.09035 | 0.086 | 0.1378 | 0.1814 | 0.1311 |
| 90S1 | 0.00571 | 0.1441 | 0.1846 | 0.2189 | 0.3554 | 0.4397 | 0.4821 | 0.4582 |
| 90S2 | 0.01388 | 0.09345 | 0.1556 | 0.2327 | 0.3901 | 0.4778 | 0.5601 | 0.587 |
| 90W1 | 0.00868 | 0.09408 | 0.1998 | 0.2035 | 0.2949 | 0.2989 | 0.408 | 0.3613 |
| 90W2 | 0.01942 | 0.1026 | 0.1381 | 0.3118 | 0.3039 | 0.4936 | 0.5051 | 0.422 |

**Panel (chart 16): 2-Methylbutanal, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.00093 | 0.00108 | 0.00098 | 0.00114 | 0.00141 | 0.00222 | 0.00189 | 0.00185 |
| 65S2 | 0.00093 | 0.00108 | 0.00098 | 0.00114 | 0.00141 | 0.00222 | 0.00189 | 0.00185 |
| 65W1 | 0.00119 | 0.00119 | 0.00099 | 0.00096 | 0.00099 | 0.00157 | 0.00168 | 0.00188 |
| 65W2 | 0.00119 | 0.00119 | 0.00099 | 0.00096 | 0.00099 | 0.00157 | 0.00168 | 0.00188 |
| 78S1 | 0.00128 | 0.00386 | 0.00801 | 0.00822 | 0.0184 | 0.0207 | 0.03701 | 0.02722 |
| 78S2 | 0.00209 | 0.01184 | 0.01006 | 0.01706 | 0.01951 | 0.04014 | 0.03758 | 0.05576 |
| 78W1 | 0.00379 | 0.00951 | 0.01362 | 0.02449 | 0.04513 | 0.04174 | 0.06837 | 0.06567 |
| 78W2 | 0.00296 | 0.01173 | 0.01574 | 0.03065 | 0.02591 | 0.05118 | 0.072 | 0.04837 |
| 90S1 | 0.00204 | 0.05231 | 0.0707 | 0.08495 | 0.1551 | 0.2034 | 0.2288 | 0.215 |
| 90S2 | 0.00253 | 0.02906 | 0.05451 | 0.08141 | 0.1651 | 0.2246 | 0.2777 | 0.2945 |
| 90W1 | 0.00298 | 0.03708 | 0.08495 | 0.08935 | 0.1412 | 0.144 | 0.2261 | 0.1983 |
| 90W2 | 0.00608 | 0.03822 | 0.05153 | 0.1368 | 0.1338 | 0.271 | 0.2734 | 0.2323 |

**Panel (chart 17): 2-Methylpropanal, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.00019 | 0.00241 | 0.00217 | 0.00238 | 0.00284 | 0.00374 | 0.00331 | 0.00311 |
| 65S2 | 0.00019 | 0.00241 | 0.00217 | 0.00238 | 0.00284 | 0.00374 | 0.00331 | 0.00311 |
| 65W1 | 0.00298 | 0.00265 | 0.0033 | 0.00378 | 0.00295 | 0.0037 | 0.00398 | 0.00455 |
| 65W2 | 0.00298 | 0.00265 | 0.0033 | 0.00378 | 0.00295 | 0.0037 | 0.00398 | 0.00455 |
| 78S1 | 0.0016 | 0.00419 | 0.00851 | 0.0086 | 0.0178 | 0.02062 | 0.03312 | 0.02606 |
| 78S2 | 0.00235 | 0.00968 | 0.00767 | 0.01392 | 0.01509 | 0.03534 | 0.03159 | 0.0471 |
| 78W1 | 0.00682 | 0.00956 | 0.01194 | 0.02386 | 0.04282 | 0.03878 | 0.06754 | 0.05624 |
| 78W2 | 0.00883 | 0.01358 | 0.02048 | 0.02952 | 0.02606 | 0.04579 | 0.06386 | 0.04622 |
| 90S1 | 0.0022 | 0.04159 | 0.06037 | 0.07179 | 0.1351 | 0.1744 | 0.2143 | 0.2033 |
| 90S2 | 0.00622 | 0.02427 | 0.04156 | 0.06846 | 0.1459 | 0.2017 | 0.2425 | 0.2806 |
| 90W1 | 0.00579 | 0.03344 | 0.07491 | 0.08022 | 0.1277 | 0.1298 | 0.2135 | 0.1889 |
| 90W2 | 0.00794 | 0.03207 | 0.04358 | 0.117 | 0.1169 | 0.2386 | 0.2432 | 0.2003 |

**Panel (chart 18): Phenylacetaldehyde, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.00239 | 0.0027 | 0.00277 | 0.00285 | 0.00365 | 0.00452 | 0.00396 | 0.00382 |
| 65S2 | 0.00239 | 0.0027 | 0.00277 | 0.00285 | 0.00365 | 0.00452 | 0.00396 | 0.00382 |
| 65W1 | 0.00298 | 0.00259 | 0.00278 | 0.00314 | 0.00249 | 0.00295 | 0.00319 | 0.0034 |
| 65W2 | 0.00298 | 0.00259 | 0.00278 | 0.00314 | 0.00249 | 0.00295 | 0.00319 | 0.0034 |
| 78S1 | 0.00247 | 0.00497 | 0.01183 | 0.0124 | 0.0279 | 0.0319 | 0.04283 | 0.03228 |
| 78S2 | 0.00372 | 0.01242 | 0.01283 | 0.01933 | 0.02488 | 0.04771 | 0.04636 | 0.06168 |
| 78W1 | 0.00463 | 0.00994 | 0.01538 | 0.02603 | 0.05127 | 0.05093 | 0.08106 | 0.08243 |
| 78W2 | 0.00616 | 0.01459 | 0.02372 | 0.03504 | 0.03491 | 0.05872 | 0.0778 | 0.05348 |
| 90S1 | 0.00263 | 0.05338 | 0.06617 | 0.09201 | 0.1503 | 0.1478 | 0.1693 | 0.1833 |
| 90S2 | 0.00591 | 0.03627 | 0.06014 | 0.09585 | 0.1536 | 0.1536 | 0.1747 | 0.188 |
| 90W1 | 0.00394 | 0.0379 | 0.08543 | 0.1008 | 0.1603 | 0.1706 | 0.1952 | 0.2235 |
| 90W2 | 0.00595 | 0.04067 | 0.05888 | 0.1189 | 0.125 | 0.2011 | 0.227 | 0.1991 |

**Panel (chart 19): Methional, mmol/kg**

| run | 0 h | 1.2 h | 2.4 h | 3.6 h | 4.8 h | 6 h | 7.2 h | 8.4 h |
|---|---|---|---|---|---|---|---|---|
| 65S1 | 0.00036 | 0.00046 | 0.0005 | 0.00051 | 0.00067 | 0.00074 | 0.00073 | 0.00079 |
| 65S2 | 0.00036 | 0.00046 | 0.0005 | 0.00051 | 0.00067 | 0.00074 | 0.00073 | 0.00079 |
| 65W1 | 0.00028 | 0.0003 | 0.00029 | 0.00041 | 0.00032 | 0.0005 | 0.00048 | 0.00056 |
| 65W2 | 0.00028 | 0.0003 | 0.00029 | 0.00041 | 0.00032 | 0.0005 | 0.00048 | 0.00056 |
| 78S1 | 0.00021 | 0.00058 | 0.00132 | 0.00137 | 0.00388 | 0.00441 | 0.00906 | 0.00576 |
| 78S2 | 0.00035 | 0.00133 | 0.0013 | 0.00202 | 0.00307 | 0.00906 | 0.00852 | 0.01345 |
| 78W1 | 0.00059 | 0.00106 | 0.00168 | 0.0032 | 0.00908 | 0.00853 | 0.01844 | 0.01865 |
| 78W2 | 0.00121 | 0.00154 | 0.00308 | 0.00493 | 0.00424 | 0.01076 | 0.02244 | 0.01133 |
| 90S1 | 0.00028 | 0.01484 | 0.02435 | 0.03516 | 0.0831 | 0.1139 | 0.1219 | 0.1198 |
| 90S2 | 0.00119 | 0.00553 | 0.01467 | 0.03944 | 0.1148 | 0.1593 | 0.1787 | 0.2141 |
| 90W1 | 0.00053 | 0.00843 | 0.03044 | 0.03338 | 0.0577 | 0.0638 | 0.09318 | 0.08511 |
| 90W2 | 0.00076 | 0.00677 | 0.01134 | 0.04855 | 0.04821 | 0.09956 | 0.09736 | 0.08264 |

### Fig. S2 — *"Predicted concentration of the compounds non-quantified analytically"* (FruMet, FruAAj, Int1, MRP, SDP), model
output only, same 12 runs. Not transcribed (not data).

**Model use.** **Yes, with the integrity caveats.** Fig. S1 gives the full per-run data behind the fit: ARPs (FruLeu/Ile/Val/Phe/
AAi) and the five Strecker aldehydes are usable 12-run × 8-point series (but the 65 °C aldehyde replicates are copies, so 65 °C has
2, not 4, independent aldehyde runs); the amino acids are usable except 90W2 at 7.2/8.4 h; glucose and fructose are too erratic to
constrain sugar loss. Table S1 adds pre-kiln moisture only.
