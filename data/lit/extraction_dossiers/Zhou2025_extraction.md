# Zhou, Diao, Zhang, Yu, Wang, Gu, Ren, Li, Dong & Yi 2025 (10.1016/j.lwt.2025.117469) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/Zhou2025.pdf` (10 pp.). Born-digital Elsevier PDF, clean text layer.
Read method: **both** — full text layer read directly, **plus** 400 dpi page-8 raster
(`z_p8-08.png`, crops `z_fig5ABC.png`, `z_fig5DEF.png`) for Figure 5, which carries the only
concentration–time data and the decisive axis label.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY

| field | value |
|---|---|
| Authors | Xiong Zhou, Yuhua Diao, Jiawen Zhang, Xinxin Yu, Gang Wang, Yin Gu, Dabing Ren, Siyu Li\*\*, Wenjiang Dong\*\*\*, Lunzhao Yi\* |
| Title | **"Sulfury/roasty fading indicators in roasted coffees: Their contribution and applicability in coffee freshness perception and prediction"** |
| Venue | *LWT – Food Science and Technology* **218** (2025), Article **117469** |
| DOI | **10.1016/j.lwt.2025.117469** |
| Dates | Received 10 May 2024; revised 26 January 2025; accepted 30 January 2025; online 31 January 2025 |
| Licence | Open access, CC BY |
| Affiliations | Kunming Univ. of Science & Technology; Spice and Beverage Research Institute, CATAS (Wanning, Hainan); Kunming Institute for Food and Drug Control; Kunming Gemmy Food Co. |
| PDF character | born-digital; text layer complete and reliable; **all figures are raster with no numeric text layer** |

**Correct file for its expected identity.** This is indeed the only candidate in the corpus for
an **FFT first-order-fade + Arrhenius model**.

---

## 1. ONE-PARAGRAPH VERDICT — READ THIS BEFORE USING ANY NUMBER HERE

The paper does deliver what the brief hoped for in form: **nine first-order decay constants
(3 compounds × 3 temperatures), three Arrhenius regressions, and three extrapolated 25 °C
rate constants**, all for whole roasted coffee beans under accelerated storage. But the
sulfur-relevant part of it **does not survive an internal-consistency check, and the failure is
specifically in the 2-furfurylthiol row**. Three findings, in descending order of importance:

1. **★ THE PRINTED FFT ARRHENIUS FIT IS NOT A FIT TO THE FFT RATE CONSTANTS.** Table 2 gives
   FFT k = 0.035 / 0.039 / 0.042 day⁻¹ at 40 / 50 / 60 °C, i.e. **Ln k = −3.352 / −3.244 /
   −3.170**. Figure 5D, the panel captioned "D: 2-furfurylthiol", has a **Ln k axis spanning
   −3.8 to −4.2** — the FFT values lie **entirely off that axis** — and its three plotted points
   digitise to **k = 0.0210 / 0.0190 / 0.0150 day⁻¹**, which is the **methanethiol** magnitude
   range. Regressing Table 2's own FFT constants gives **slope −0.953, Ea = 7.92 kJ/mol**;
   the paper prints **slope −1.762, Ea = 14.65 kJ/mol**. **The paper's FFT activation energy is
   ~1.85× too large relative to its own FFT data**, and its extrapolated 25 °C FFT rate
   (0.0102 day⁻¹, t½ = 68 d) is **~3× slower** than the honest extrapolation from its own
   FFT constants (0.0302 day⁻¹, t½ = 23 d). See §6 for the full forensic trail.
2. **★ EVEN THE HONEST FFT Ea IS NOT A CHEMICAL ACTIVATION ENERGY.** k rises only from 0.035 to
   0.042 day⁻¹ across a 20 K span — **Q₁₀ = 1.11 (40→50 °C) and 1.08 (50→60 °C)**, against
   Q₁₀ ≈ 2.4 for a genuine 80 kJ/mol chemical process. Ea ≈ 8 kJ/mol is in the
   diffusion/partition regime, not the bond-breaking regime. **Whatever limits FFT loss from
   whole beans over 30 days at 40–60 °C, it is almost certainly not the thiol–quinone /
   thiol–melanoidin chemistry the reaction network models.** Treat this as a *matrix-transport*
   dataset, not a *reaction-kinetics* dataset.
3. **The Table 1 OAV block is irreconcilable with the Table 2 / Figure 5 concentrations** by
   three different, non-constant factors (18.5× for FFT, 102.7× for methanethiol, 10.0× for
   2-ethyl-3-methyl-pyrazine). §4.1. **Table 1 must not be used as a concentration source.**

**There are no elementary rate constants, no mechanism, no pH axis, no water-activity axis,
and no reactant concentrations.** Every k here is a lumped, empirical, whole-bean
disappearance constant at one moisture and one packaging condition. Quantitation is
**semi-quantitative** (single-point internal-standard peak-area ratio, no authentic-standard
calibration curve), and there is **no SI on disk** (Tables S1–S3, which hold the 86-compound
concentration matrix, are not in the repo).

---

## 2. SYSTEM DEFINITION — verbatim

### 2.1 Material and roast (§2.2, p. 2)

| variable | value as printed |
|---|---|
| Botanical material | **Robusta *coffea canephora*** fruits, experimental station of the Spice and Beverage Research Institute, CATAS, Hainan, China |
| Harvest | **March 2023**; "commercially mature", "no visible defects and physical damage" |
| Post-harvest | twigs/leaves removed manually; fruits **mechanically peeled and degummed**; wet beans held at **4 °C for a maximum of 24 h** before roasting |
| Pre-roast drying | **hot air drying until moisture < 12 %** |
| Replication of roast | **triplicate, 500 g each** |
| Roaster | PRE 1 Z drum Coffee Roaster (Emmerich am Rhein, Germany) |
| Roast profile | **180 °C for 7.5 min** — a single, light-to-medium profile; **no other roast level was studied** |
| Post-roast | "quickly cooled down to room temperature" |

### 2.2 Accelerated storage — the kinetic experiment (§2.2, p. 2)

> "were then submitted to three different accelerated storage experiments (**40 °C, 50 °C,
> 60 °C**) in individual electro-thermal constant-temperature dry boxes (DHG-9030A, Shanghai
> Yiheng Laboratory Instrument Co., Ltd., Shanghai, China) for **30 days** to mimic the
> long-time storage of coffee under normal room temperature, and were **collected every five
> days**. After collection, samples were put into food-grade packaging bags, sealed, degassed
> and stored at **− 20 °C**."

| variable | value |
|---|---|
| Temperatures | **40, 50, 60 °C** (= 313.15, 323.15, 333.15 K, as the paper itself converts them, p. 6) |
| Duration | **30 days** |
| Sampling | **every 5 days** → t = 0, 5, 10, 15, 20, 25, 30 d (**7 points**) |
| Storage state | **whole roasted beans** in an electro-thermal dry box |
| Humidity | **NOT CONTROLLED AND NOT REPORTED** — "dry box" is all that is said. **No a_w, no RH, no bean moisture during storage.** This is a first-order gap: an FFT fade rate without a water activity is not transferable |
| Oxygen / atmosphere during storage | **NOT REPORTED.** Nothing states whether beans were in air, sealed, or degassed *during* the 40–60 °C hold; the sealing/degassing described applies **after** collection, for the −20 °C archive |
| Packaging during storage | not stated (only the post-collection food-grade bag is described) |
| Light | not stated |

### 2.3 Brewing — what was actually analysed (§2.3, p. 2)

Every concentration in this paper is measured **in the brew, not in the bean.**

> "all samples were crushed by an ultra-fine pulverizer (VTA-6S3 MAHLKONIG, Germany) with a
> **grinding degree of 3** and passed through a **60-mesh sieve** … the French filter press
> kettle was preheated with deionized water at **92 °C for 10 s**, subsequently, **10 g of
> coffee powders** were put into a French filter press pot, and **200 mL of 92 °C deionized
> water** was poured in … then extracted for **5 min**, pressed the handle of the press kettle,
> and gently shook the press kettle for 5 s"

| variable | value |
|---|---|
| Brew ratio | **10 g / 200 mL = 1:20 (w/v), 50 g/L** |
| Brew water | deionized water, **92 °C** |
| Contact time | **5 min**, French press |
| Brew pH | **NOT MEASURED AND NOT REPORTED** — a significant omission for a thiol paper |

### 2.4 Analysis (§2.4, p. 2)

| variable | value as printed |
|---|---|
| Method | **HS-SPME-GC-MS**, per Dong et al. 2019 |
| Internal standard | **3-heptanone**, purity > 98.0 %, Sigma-Aldrich (Poole, UK); prepared as **10 µL 3-heptanone in 100 mL methanol**; **30 µL of that solution added to each 10 mL sample**, held **1 h** before analysis |
| Sample | **10 mL brew in a 20 mL headspace vial** |
| Fibre | 30/50 µm **CAR/PDMS/DVB**, aged at 250 °C for 1 h |
| Equilibration / extraction | stabilise at **40 °C for 20 min**, adsorb **40 min**, desorb at **250 °C for 5 min** |
| GC-MS | Agilent 7890A GC / 5975C MSD, quadrupole |
| Column | **DB-WAX, 30 m × 0.25 mm × 0.25 µm** |
| Oven | 40 °C hold 5 min → 1.5 °C/min to 130 °C → 8 °C/min to 200 °C → 10 °C/min to 240 °C, hold 5 min |
| Carrier | helium, **1.0 mL/min**; splitless injection at 250 °C |
| MS | EI; ion source 250 °C, quadrupole 150 °C, transfer line 250 °C; full scan **30–350 m/z at 3.06 scan/s** |
| Identification | NIST14.0, **match similarity > 80 %**, plus LRI from C8–C40 alkanes |
| **Quantitation** | > "The relative content of each compound was calculated by **comparing its peak area with that of the internal standard**." — **SEMI-QUANTITATIVE. No authentic-standard calibration curve, no response factors, no recovery, no LOD/LOQ.** The µg/L axis labels in Fig. 5 are therefore 3-heptanone-equivalents, not true FFT mass concentrations |
| Replication | **"GC-MS analysis were completed in triplicate"** (n = 3) |
| Error convention | §2.7: "all data are expressed as the average standard deviation (SD) of three tests"; ANOVA + Duncan's new multiple range test, two-tailed, **p < 0.05**. **No error bars are printed on any kinetic figure or in Table 2** |

### 2.5 Sensory (§2.5, p. 3) — context for the S-R index

Quantitative descriptive analysis (QDA). Panel: **10 graduate volunteers** (5 F, 5 M, aged
20–25) of the Spice and Beverage Institute, trained with **Le Nez du Café**. Five descriptors
generated by the panel: **sulfury/roasty, smoky, nutty, caramelized, overall aroma intensity**.
Serving: **20 mL per aliquot**, identical cups, room at **26 °C**, unsalted biscuits + purified
water as palate cleanser, **10-min rest between samples**. **No scale range is stated and no
sensory numbers are tabulated anywhere in the paper** — Fig. 2D–F are radar plots only.

### 2.6 The two model equations, verbatim (§2.6, p. 3)

> "Ct = C0 e⁻ᵏᵗ
> Ct represents the concentration of a specific volatile compounds at a specific storage time t.
> C0 represents the initial concentration of this specific volatile compound. **k represents the
> reaction rate constant at a specific temperature.** t represents the storage time."

> "k = k0 exp(− Ea/RT)
> k0 refers to the pre-factor; **Ea is the activation energy (J/mol)**; **R is the molar gas
> constant of 8.314 (KJ/mol)**; T is the temperature expressed as Kelvin degree (K)."

⚠️ **UNIT DEFECT #1 (methods).** *"R is the molar gas constant of 8.314 (KJ/mol)"* is wrong
twice over: R = 8.314 **J** mol⁻¹ **K⁻¹**, not 8.314 kJ/mol, and the printed form omits K⁻¹
entirely. The numeric value 8.314 is right; the unit string is not.

⚠️ **UNIT DEFECT #2 (§3.3.2, p. 6).** > "the Arrhenius formula was modified into Ln k =
− Ea/RT + Ln k0, **Lnk and 1/T were set as Y and X variables**". Taken literally this makes
the printed slopes (−1.762, −1.233, −2.094) equal to −Ea/R with x in K⁻¹, giving
**Ea = 14.6 / 10.3 / 17.4 J mol⁻¹** — absurd. **The x-axis of Fig. 5D–F is printed
`1000/T`** (read off the 400 dpi raster; the label is unambiguous). **The self-consistent
reading is x = 1000/T, hence Ea = |slope| × 8.314 kJ mol⁻¹.** All Ea values below use that
reading and are tagged **[Z]** because the authors never print an Ea value anywhere in the
paper.

**Units of k.** The x-axis of Fig. 5A–C is **"Storage time (day)"**, so every k in Table 2 is in
**day⁻¹**. The paper never states the unit of k. Tagged **[Z]** where it matters.

---

## 3. TABLE 2 — THE COMPLETE KINETIC RESULT (the whole quantitative core of the paper)

**Anchor: Table 2, p. 8 (PDF page 8).** Title as printed: *"Formulas obtained from the
mathematical models in Fig. 5."*
Upper block column headers as printed: `Compounds | Temperature (°C) | First-order kinetic
reaction formula | R² | Corresponding figure`.
Lower block column headers as printed: `Compounds | Arrhenius formula | R² | Corresponding figure`.
Transcribed verbatim from the 400 dpi raster; **all cells legible, none unreadable.**
Note the paper's own inconsistent capitalisation of the dependent variable ("Y" in the first
row, "y" in all others) — reproduced as printed.

### 3.1 Upper block — the nine first-order models **[F]**

| Compound | Temp (°C) | First-order kinetic reaction formula *as printed* | R² *as printed* | Fig. | **C₀, µg/L [Z]** | **k, day⁻¹ [Z]** | **t½ = ln2/k, day [Z]** |
|---|---|---|---|---|---|---|---|
| 2-Furfurylthiol | 40 | **Y = 12.621 exp (− 0.035x)** | **0.987** | Fig. 5A | 12.621 | **0.035** | 19.80 |
| 2-Furfurylthiol | 50 | **y = 12.500 exp (− 0.039x)** | **0.986** | Fig. 5A | 12.500 | **0.039** | 17.77 |
| 2-Furfurylthiol | 60 | **y = 12.348 exp (− 0.042x)** | **0.985** | Fig. 5A | 12.348 | **0.042** | 16.50 |
| Methanethiol | 40 | **y = 6.202 exp (− 0.015x)** | **0.984** | Fig. 5B | 6.202 | **0.015** | 46.21 |
| Methanethiol | 50 | **y = 5.976 exp (− 0.018x)** | **0.964** | Fig. 5B | 5.976 | **0.018** | 38.51 |
| Methanethiol | 60 | **y = 6.03 exp (− 0.021x)** | **0.970** | Fig. 5B | 6.03 | **0.021** | 33.01 |
| 2-Ethyl-3-methyl-pyrazine | 40 | **y = 159.330 exp (− 0.018x)** | **0.992** | Fig. 5C | 159.330 | **0.018** | 38.51 |
| 2-Ethyl-3-methyl-pyrazine | 50 | **y = 157.205 exp (− 0.0211x)** | **0.953** | Fig. 5C | 157.205 | **0.0211** | 32.85 |
| 2-Ethyl-3-methyl-pyrazine | 60 | **y = 153.556 exp (− 0.025x)** | **0.979** | Fig. 5C | 153.556 | **0.025** | 27.73 |

**C₀ units [Z]:** the paper never labels C₀. The Fig. 5A–C ordinates are printed
**"Concentration (µg/L)"** and the fitted intercepts (≈12.6, ≈6.2, ≈159) sit exactly on those
axes, so **C₀ is in µg/L of brew, in 3-heptanone-equivalents.**

**Note on the C₀ trend [Z]:** for all three compounds C₀ *decreases* monotonically with storage
temperature (FFT 12.621 → 12.500 → 12.348). Since t = 0 is the *same* freshly roasted bean for
all three arms, the three C₀ values are three fits to what should be one number. Their spread
(**±1.1 % for FFT, ±1.9 % for methanethiol, ±1.8 % for 2E3MP**) is a usable, if crude,
**empirical estimate of the fit/measurement precision**, and it is the only precision estimate
the paper offers on the kinetics at all.

### 3.2 Lower block — the three Arrhenius regressions **[F]**

| Compound | Arrhenius formula *as printed* | R² *as printed* | Fig. | **Ea, kJ/mol [Z]** (= \|slope\|×8.314, x = 1000/T) | **k₀, day⁻¹ [Z]** (= e^intercept) |
|---|---|---|---|---|---|
| 2-Furfurylthiol | **y = − 1.762x + 1.447** | **0.956** | Fig. 5D | **14.65** | **4.250** |
| Methanethiol | **y = − 1.233x − 0.717** | **0.997** | Fig. 5E | **10.25** | **0.4882** |
| 2-Ethyl-3-methyl-pyrazine | **y = − 2.094x + 1.814** | **0.992** | Fig. 5F | **17.41** | **6.135** |

⚠️ **Every row of this block fails at least one internal-consistency check. See §6.**

### 3.3 The extrapolated 25 °C rate constants **[F]** (body text, p. 6, §3.3.2)

> "it can be further calculated that the theoretical reacting rates of 2-furfurylthiol,
> methanethiol and 2-ethyl-3-methyl-pyrazine under actual room temperature storage condition
> (**25 °C**) were **0.0102, 0.0076, 0.0055**, respectively."

| Compound | **k(25 °C) as printed** | unit **[Z]** | **t½ at 25 °C, day [Z]** | **k(25 °C) recomputed from the printed Arrhenius formula [Z]** | deviation |
|---|---|---|---|---|---|
| 2-Furfurylthiol | **0.0102** | day⁻¹ | **68.0** | 0.01153 | **+13.0 %** |
| Methanethiol | **0.0076** | day⁻¹ | **91.2** | 0.00781 | **+2.7 %** |
| 2-Ethyl-3-methyl-pyrazine | **0.0055** | day⁻¹ | **126.0** | 0.00547 | **−0.6 %** |

The methanethiol and pyrazine values round-trip; **the FFT value does not** (13 % off its own
formula) — the third independent symptom of the FFT row being corrupt.

### 3.4 The final predictive form **[F]** (p. 6)

> "it can be speculated that the evolution formula of 2-furfurylthiol, methanethiol and
> 2-ethyl-3-methyl-pyrazine under actual storage condition should be **Y = a\*exp(kX)**, of which
> Y and X represent the actual concentrations of the compounds and the storage time,
> respectively, a and k represent the initial concentrations of the compounds and the
> theoretical reacting rates, respectively."

⚠️ **Sign defect.** As printed, `Y = a*exp(kX)` with a *positive* k is a growth law. The
intended form is `Y = a·exp(−kX)`, consistent with Table 2's `exp(−0.035x)`. Note also the
sentence swaps the definitions ("Y and X represent the actual concentrations … and the storage
time" — correct; then "a and k represent the initial concentrations … and the theoretical
reacting rates" — correct), so only the sign is wrong.

---

## 4. TABLE 1 — the OAV / odour-threshold block (COMPLETE, but see the warning)

**Anchor: Table 1, p. 5 (PDF page 5).** Title as printed: *"Theoretical aroma profiles of coffee
brews produced from freshly roasted coffee beans."*
Column headers as printed: `Volatile compound | Threshold (µg/L) ᵃ | Odor description ᵇ | OAV values`.
Footnote **a**: *"The threshold values were referenced from the works previously reported"* —
**so every threshold is [C], cited, not measured here**; cited sources: Amanpour & Selli 2016;
Belitz, Grosch & Schieberle 2009; Czerny, Christlbauer, Christlbauer et al. 2008; Giri, Osako,
Okamoto & Ohshima 2010; Miyazato, Nakamura, Hashimoto & Hayashi 2013; Nishimura & Mihara 1990;
Piccino, Boulanger, Descroix & Sing 2014; Puvipirom & Chaiseri 2012; Steinhaus & Schieberle 2007.
Footnote **b**: odour descriptions likewise **[C]**.
The printed table groups compounds under one merged odour-description cell per block; the
description is reproduced against every member of its block below.

### 4.0 Full transcription

| Volatile compound | Threshold (µg/L) **[C]** | Odor description **[C]** | OAV values **[F]** |
|---|---|---|---|
| Methanethiol | **0.20** | Sulfury/roasty (roasted coffee, roasted barley, roasted nuts, sulfur, alliaceous, onion) | **3184.35** |
| 5-Methyl-2-furanmethanethiol | **0.05** | Sulfury/roasty (as above) | **18411.00** |
| 2-[(Methyldithio)methyl]-furan | **0.17** | Sulfury/roasty (as above) | **1605.12** |
| **2-Furfurylthiol** | **0.01** | Sulfury/roasty (as above) | **23297.83** |
| 2-Methyl-3-furanthiol *(printed "2-Methy1-3-furanthiol", OCR/typesetting of the l as 1)* | **0.05** | Sulfury/roasty (as above) | **525.00** |
| Dimethyl disulfide | **0.82** | Sulfury/roasty (as above) | **119.21** |
| Dimethyl trisulfide | **0.001** | Sulfury/roasty (as above) | **37620.00** |
| 2-[(Methylthio)methyl]-furan | **0.40** | Sulfury/roasty (as above) | **4197.95** |
| 2,4-Dimethyl-thiazole | **18.00** | Sulfury/roasty (as above) | **21.61** |
| 2,6-Dimethyl-pyrazine | **1.72** | Sulfury/roasty (as above) | **875.92** |
| **Total sulfury/roasty OAV** | — | — | **88857.99** |
| 2,5-Dimethyl-pyrazine | **1.82** | Nutty (baked corn, cocoa, baked hazelnut) | **893.37** |
| **2-Ethyl-3-methyl-pyrazine** | **130.00** | Nutty (as above) | **12.25** |
| 2,3-Dimethyl-pyrazine | **0.88** | Nutty (as above) | **584.45** |
| 2-Ethyl-6-methyl-pyrazine | **40.00** | Nutty (as above) | **42.53** |
| 2,5-Diethyl-pyrazine | **20.00** | Nutty (as above) | **37.86** |
| 3-Ethyl-2,5-dimethyl-pyrazine | **0.40** | Nutty (as above) | **3991.23** |
| 2,6-Diethyl-pyrazine | **6.00** | Nutty (as above) | **64.61** |
| **Total nutty OAV** | — | — | **5626.3** |
| 2-Ethyl-3,5-dimethyl-pyrazine | **0.04** | Smoky (green clove, curry, phenolic) | **13749.75** |
| 1-Methyl-1H-pyrrole-2-carboxaldehyde | **37.00** | Smoky (as above) | **41.34** |
| 1-(2-Furanylmethyl)-1H-pyrrole | **100.00** | Smoky (as above) | **17.69** |
| Ethyl acetate | **1.30** | Smoky (as above) | **348.74** |
| 2-Methoxy-phenol | **2.50** | Smoky (as above) | **791.56** |
| 4-Ethyl-2-methoxy-phenol | **25.00** | Smoky (as above) | **56.76** |
| **Total smoky OAV** | — | — | **15005.84** |
| 5-Methyl-2-furancarboxaldehyde | **20.00** | Caramelized (acidic almond, burnt sugar, creamy, fruity, sweet) | **119.97** |
| 2,3-Pentanedione | **5.00** | Caramelized (as above) | **45.96** |
| 3-Methyl-4-heptanone | **0.05** | Caramelized (as above) | **12323.40** |
| 2-Methyl-butanal | **1.30** | Caramelized (as above) | **2348.58** |
| 2-Methyl-benzaldehyde | **1.50** | Caramelized (as above) | **338.94** |
| Ethyl acetate | **1.30** | Caramelized (as above) | **348.74** |
| **Total caramelized OAV** | — | — | **15525.59** |
| Methyl-pyrazine | **105.00** | pungent | **14.66** |
| Pyridine | **77.00** | rancid | **16.97** |
| 2-Ethyl-5-methyl-pyrazine | **16.00** | fruity, pungent, sweet | **84.19** |
| 2-Methyl-5-propyl-pyrazine | **2.00** | earthy | **173.07** |
| 3-Methylbutylester-2-butenoic acid | **0.50** | paper | **309.72** |
| 2-Ethyl-phenol | **0.260** | phenolic | **629.19** |

**Ethyl acetate is printed twice** (once under smoky, once under caramelized) with the same
threshold **1.30** and the same OAV **348.74** — it is double-counted into both group totals.
**[Z]** verification: the "Total sulfury/roasty OAV" 88857.99 is **not** the sum of its ten
member rows (which sum to **89 837.99**); the difference is exactly **980.00**. The nutty block
sums to **5626.30** ✓ matches. Smoky sums to **15 005.84** ✓ matches. Caramelized sums to
**15 525.59** ✓ matches. **So only the sulfury/roasty total is arithmetically wrong**, by 980.

### 4.1 ⚠️ TABLE 1 CANNOT BE RECONCILED WITH TABLE 2 / FIGURE 5 — DO NOT USE IT FOR CONCENTRATIONS

OAV ≡ concentration / threshold. Inverting Table 1 gives an implied day-0 brew concentration;
comparing it with the day-0 concentration the same paper prints in Table 2:

| Compound | Threshold, µg/L | OAV (T1) | **implied [C] = OAV × threshold, µg/L [Z]** | **C₀ from Table 2, µg/L** | **ratio [Z]** |
|---|---|---|---|---|---|
| 2-Furfurylthiol | 0.01 | 23297.83 | **232.98** | **12.621** | **18.46×** |
| Methanethiol | 0.20 | 3184.35 | **636.87** | **6.202** | **102.7×** |
| 2-Ethyl-3-methyl-pyrazine | 130.00 | 12.25 | **1592.50** | **159.330** | **10.00×** |

**Three compounds, three different discrepancy factors, none of them 1.** There is no single
unit error, dilution factor, or bean-vs-brew basis that reconciles them. Conclusion:
**Table 1's OAVs and Table 2/Figure 5's concentrations are on mutually inconsistent scales, and
at most one of them can be right.** The kinetics (Table 2 / Fig. 5) are internally coherent with
each other and with the figure axes, so **Table 2/Fig. 5 is the survivable block; Table 1's OAV
column is not usable as a quantitative input.** The *thresholds* column remains usable — but it
is **[C]**, cited from nine other papers, and should be taken from those primary sources.

---

## 5. FIGURE 5A–C — THE CONCENTRATION–TIME DATA **[fig]**

**Anchor: Fig. 5A–C, p. 8 (PDF page 8).** Caption as printed: *"First-order kinetic reaction
models (A: 2-furfurylthiol; B: methanethiol; C: 2-ethyl-3-methyl-pyrazine) and Arrhenius
formulas (D: 2-furfurylthiol; E: methanethiol; F: 2-ethyl-3-methyl-pyrazine) elaborating the
correlations between concentrations and storage time, storage temperature (T) and reacting
rates (k) of sulfury/roasty correlated volatile compounds, respectively."*
Ordinate as printed: **"Concentration (µg/L)"**. Abscissa: **"Storage time (day)"**.
Series: 40 °C = black solid line + filled squares; 50 °C = red dashed + filled circles;
60 °C = blue dash-dot + filled triangles.

**These are the only measured data points in the paper.** Digitised from the 400 dpi raster
(`z_fig5ABC.png`). **No error bars are drawn on any point**, despite n = 3 and the stated SD
convention. Digitisation uncertainty stated per panel.

### 5.1 Fig. 5A — 2-furfurylthiol, µg/L **[fig]** (±0.15 µg/L)

| Storage day | **40 °C** | **50 °C** | **60 °C** |
|---|---|---|---|
| 0 | 12.8 | 13.0 | 12.7 |
| 5 | 10.0 | 9.6 | 9.5 |
| 10 | 8.2 | 7.8 | 7.55 |
| 15 | 7.3 | 6.85 | 6.8 |
| 20 | 6.1 | 5.7 | 5.6 |
| 25 | 5.1 | 5.0 | 4.85 |
| 30 | 4.6 | 3.8 | 3.6 |

Axis ticks printed: 3, 6, 9, 12, 15 µg/L. **Day-0 points of the three arms are visually
coincident at ≈12.7–13.0 µg/L**, consistent with one shared fresh-bean starting material.
**Overall 30-day loss [Z]: 64 % (40 °C), 71 % (50 °C), 72 % (60 °C).** The three arms are
separated by less than the day-to-day scatter until about day 20 — a visual restatement of the
Q₁₀ ≈ 1.1 finding.

### 5.2 Fig. 5B — methanethiol, µg/L **[fig]** (±0.06 µg/L)

| Storage day | **40 °C** | **50 °C** | **60 °C** |
|---|---|---|---|
| 0 | 6.35 | 6.35 | 6.35 |
| 5 | 5.60 | 5.29 | 5.09 |
| 10 | 5.28 | 4.90 | 4.79 |
| 15 | 4.99 | 4.67 | 4.47 |
| 20 | 4.72 | 4.33 | 4.14 |
| 25 | 4.25 | 3.92 | 3.77 |
| 30 | 4.01 | 3.65 | 3.32 |

Axis ticks printed: 3, 4, 5, 6, 7 µg/L. The three day-0 points are exactly superimposed at
**≈6.35 µg/L**, notably **above** all three fitted C₀ values (6.202 / 5.976 / 6.03) — the fits
systematically undershoot t = 0. **Overall 30-day loss [Z]: 37 % / 43 % / 48 %.**

### 5.3 Fig. 5C — 2-ethyl-3-methyl-pyrazine, µg/L **[fig]** (±2 µg/L)

| Storage day | **40 °C** | **50 °C** | **60 °C** |
|---|---|---|---|
| 0 | 160 | 158 | 154 |
| 5 | 145 | 140 | 131 |
| 10 | 133 | 126 | 116 |
| 15 | 120 | 111 | 103 |
| 20 | 111 | 102 | 93 |
| 25 | 102 | 93 | 84 |
| 30 | 96 | 86 | 75 |

Axis ticks printed: 80, 120, 160 µg/L. **Overall 30-day loss [Z]: 40 % / 46 % / 51 %.**

---

## 6. ★ THE FFT ARRHENIUS FORENSICS — the single most important section of this dossier

### 6.1 The axis label settles the unit question

Read off the 400 dpi raster (`z_fig5DEF.png`), the common abscissa label under Fig. 5F is
printed **`1000/T`**. The methods text (§2.6) says the X variable is `1/T`. **The figure is
right and the text is wrong.** Consequence: **Ea = |slope| × 8.314 kJ mol⁻¹**, not J mol⁻¹.
With this reading, the printed 25 °C extrapolations for methanethiol and 2E3MP round-trip to
within 2.7 % and 0.6 % (§3.3), which independently confirms x = 1000/T.

### 6.2 Fig. 5D does not plot 2-furfurylthiol

**Panel D's Ln k axis is printed −3.8 / −3.9 / −4.0 / −4.1 / −4.2.**
The FFT rate constants from Table 2 are 0.035 / 0.039 / 0.042 day⁻¹, i.e.
**Ln k = −3.352 / −3.244 / −3.170 — every one of them above the top of panel D's axis.**
It is not possible for panel D to be a plot of the FFT constants.

Digitising panel D's three plotted squares (x = 1000/T = 3.002, 3.094, 3.193, i.e.
T = 60, 50, 40 °C) gives:

| panel | digitised Ln k at 60 / 50 / 40 °C **[fig]** | implied k, day⁻¹ **[Z]** | matches which published triple? |
|---|---|---|---|
| **D (captioned "2-furfurylthiol")** | −3.862 / −3.962 / −4.200 | **0.0210 / 0.0190 / 0.0150** | the **methanethiol** magnitude range — *not* FFT |
| **E (methanethiol)** | −3.962 / −4.075 / −4.200 | **0.0190 / 0.0170 / 0.0150** | **exactly** the body-text methanethiol triple (p. 6): "0.015 under 40 °C, 0.017 under 50 °C, 0.019 under 60 °C" ✓ |
| **F (2-ethyl-3-methyl-pyrazine)** | −3.670 / −3.857 / −4.020 | **0.0255 / 0.0211 / 0.0180** | **exactly** the Table 2 2E3MP triple (0.025 / 0.0211 / 0.018) ✓ |

**Panels D and E plot k values in the same 0.015–0.021 day⁻¹ band.** Two panels for two
different compounds cannot both be right when the compounds' Table-2 constants differ by 2.3×.

### 6.3 Regression check — what each printed formula actually is a fit to

Ordinary least squares of Ln k on 1000/T (T = 313.15 / 323.15 / 333.15 K):

| candidate k triple (40/50/60 °C) | fitted slope | fitted intercept | R² | Ea, kJ/mol |
|---|---|---|---|---|
| **FFT, Table 2** (0.035, 0.039, 0.042) | **−0.9526** | −0.3059 | 0.9920 | **7.92** |
| FFT, body text (0.035, 0.039, 0.041) | −0.8282 | −0.6990 | 0.9638 | 6.89 |
| **MeSH, Table 2** (0.015, 0.018, 0.021) | **−1.7561** | **+1.4110** | 0.9991 | 14.60 |
| MeSH, body text (0.015, 0.017, 0.019) | **−1.2334** | −0.2599 | 0.9997 | 10.26 |
| 2E3MP, Table 2 (0.018, 0.0211, 0.025) | −1.7125 | +1.4477 | 0.9987 | 14.24 |
| 2E3MP, body text (0.017, 0.021, 0.025) | −2.0131 | +2.3580 | 0.9986 | 16.74 |
| *digitised panel D* | −1.778 | +1.496 | 0.9564 | 14.78 |
| *digitised panel E* | −1.246 | −0.220 | 0.9999 | 10.36 |
| *digitised panel F* | −1.830 | +1.818 | 0.9963 | 15.22 |

Now compare against what is printed:

| printed row | printed slope / intercept / R² | nearest regression | verdict |
|---|---|---|---|
| **2-Furfurylthiol** | **−1.762 / +1.447 / 0.956** | **MeSH Table 2** (−1.7561 / +1.4110), and digitised panel D (−1.778 / +1.496, R² 0.9564 — R² matches to 3 dp) | **★ The printed "FFT" Arrhenius fit is a fit to the METHANETHIOL rate constants.** It reproduces its own figure perfectly and its own compound's data not at all |
| **Methanethiol** | **−1.233 / −0.717 / 0.997** | slope matches MeSH body text (−1.2334) and digitised panel E (−1.246) **exactly**; **intercept does not** (should be ≈ −0.22 to −0.26) | **slope correct, intercept wrong.** With intercept −0.717 the line at x = 3.2 would sit at Ln k = −4.66, **below the bottom of panel E's own axis** (−4.2). Printing error |
| **2-Ethyl-3-methyl-pyrazine** | **−2.094 / +1.814 / 0.992** | intercept matches digitised panel F (+1.818); **slope does not** (panel F gives −1.830, 2E3MP body text gives −2.013) | **intercept correct, slope suspect.** With slope −2.094 and intercept +1.814 the line at x = 3.193 sits at Ln k = −4.87, again **below panel F's axis floor** (−4.1). Printing error |

**None of the three printed Arrhenius rows is fully self-consistent.** The FFT row is the worst
case because it is not merely mistyped — it is a fit to the wrong compound's data.

### 6.4 The honest FFT numbers, and the size of the error

| quantity | **paper's printed value** | **honest value recomputed from the paper's own Table 2 FFT constants [Z]** | error factor |
|---|---|---|---|
| Ea (FFT) | 14.65 kJ/mol | **7.92 kJ/mol** | **1.85× too high** |
| k₀ (FFT) | 4.250 day⁻¹ | **0.7365 day⁻¹** | 5.8× too high |
| Arrhenius R² | 0.956 | **0.9920** | — |
| **k(25 °C), FFT** | **0.0102 day⁻¹** | **0.0302 day⁻¹** | **2.96× too slow** |
| **t½(25 °C), FFT** | **68.0 day** | **23.0 day** | **2.96× too long** |

**If any wave ingests this paper's headline FFT shelf-life number (68 days at 25 °C), it will be
ingesting a ~3× error that the paper's own Table 2 refutes.**

### 6.5 Even the honest Ea is not a chemical activation energy

| quantity **[Z]** | value |
|---|---|
| Q₁₀ (40 → 50 °C), FFT | **1.114** |
| Q₁₀ (50 → 60 °C), FFT | **1.077** |
| Q₁₀ implied by a typical 80 kJ/mol chemical step at 50 °C | ≈ 2.4 |
| Ea (FFT, honest) | **7.92 kJ/mol** |
| Ea (methanethiol, from Table 2 triple) | 14.60 kJ/mol |
| Ea (2E3MP, from Table 2 triple) | 14.24 kJ/mol |

An Ea of 8–15 kJ/mol is the range of **viscous flow, diffusion, and partition**, not of covalent
bond formation or cleavage. For reference, the corpus's own measured sulfur chemistry sits an
order of magnitude higher (Zheng & Ho 1994 H₂S release from cysteine: **~123 kJ/mol at pH 9**).
**Structural conclusion: over 30 days at 40–60 °C in whole beans, FFT disappearance is
transport-limited or driven by an already-present oxidant pool, not by a thermally activated
Maillard/quinone step.** This makes the dataset valuable as a **matrix-retention / storage-loss**
constraint and close to useless as a *reaction-network* constraint — and it is a direct,
quantitative caution against any repo term that extrapolates a thiol-consumption Arrhenius law
from roasting temperatures down to storage temperatures.

---

## 7. QUALITATIVE / STRUCTURAL FINDINGS (no numbers attached, but they constrain the model)

**Anchors: §3.2.1–3.2.2, pp. 4–5; §3.3.1, p. 6; Fig. 3, Fig. 4; Conclusion, p. 7.**

| # | finding | provenance |
|---|---|---|
| 7.1 | **86 compounds identified**: furans 16, pyrazines 14, pyrroles 10, sulfurous 8, thiazoles 8, aldehydes 7, phenols 5, ketones 5, alcohols 4, pyridines 3, esters 2, others 4 (styrene, trimethyl-oxazole, 2,4,5-trimethyl-1,3-dioxolane, toluene) | [M] §3.1.1, p. 3 |
| 7.2 | Class shares of total volatiles, ranges across all conditions: furans **30.85–37.16 %**, pyrazines **17.72–21.95 %**, pyrroles **8.54–10.16 %**, ketones **7.82–9.57 %**, phenols **7.54–9.01 %**, aldehydes **5.30–6.95 %**, sulfurous **3.80–5.32 %**, thiazoles **0.93–1.68 %** | [M] §3.1.1, p. 3, Fig. 1D–F. ⚠️ the sentence lists eight class names against eight ranges but the name order ("furans, pyrazines, pyrroles, ketones, phenols aldehydes, sulfurous compounds, and thiazoles") is missing a comma and the pairing of "phenols aldehydes" to two ranges is ambiguous as printed |
| 7.3 | **35 compounds have OAV > 1** | [M] §3.1.1, p. 3 |
| 7.4 | PLS-R correlation coefficients (R² of blue lines) between volatiles and sulfury/roasty score: **0.99 (40 °C), 0.97 (50 °C), 0.90 (60 °C)** | [F] §3.2.1, p. 4, Fig. 3 |
| 7.5 | **Positive** S-R correlates at 40/50 °C: 2-methyl-3-furanthiol, methanethiol, 2-[(methylthio)methyl]-furan, 2-[(methyldithio)methyl]-furan, **2-furfurylthiol**, 2,4-dimethyl-thiazole, 2,6-dimethyl-pyrazine, 2-ethyl-3-methyl-pyrazine | [F] §3.2.1, p. 4 |
| 7.6 | At **60 °C** the set changes: 2-[(methyldithio)methyl]-furan and 2,4-dimethyl-thiazole **drop out**, replaced by **5-methyl-2-furanmethanethiol**. A genuine temperature-dependent change in the correlate set | [F] §3.2.1, p. 4 |
| 7.7 | **Negative** S-R correlates at all temperatures: **dimethyl disulfide, methyl-pyrazine, pyridine, 2-methyl-butanal**; these "all increased significantly during the first 15–25 days of the storage, then remained stable" | [F]/[M] §3.2.1 p. 4 and §3.2.2 p. 5, Fig. 4 |
| 7.8 | ★ **Dimethyl disulfide is treated as a THIOL SINK MARKER**: > "dimethyl disulfide may be generated by the oxidation reaction or melanoidin-mediated non-covalent interactions that uses other low-threshold corresponding thiol as the basic precursors … therefore, the presence of dimethyl disulfide indicates a decreased sulfury/roasty OAV" — an explicit thiol → disulfide oxidative channel operating at 40–60 °C in the *dry bean* | [M] (mechanistic claim, [C] to Gigl/Hofmann/Frank 2021; Hofmann & Schieberle 2002; Hofmann, Schieberle & Grosch 1996) §3.2.2, p. 5 |
| 7.9 | Proposed FFT loss channels, all cited: interaction with **melanoidin polymers**, **hydroxyhydroquinone** and its oxidation product **quinone**; reaction with **chlorogenic acids** and their thermal degradation products such as **(4-ethyl)catechol** | [C] §3.2.2, p. 5 (Gigl, Hofmann & Frank 2021; Müller & Hofmann 2007; Gigl, Frank, Irmer & Hofmann 2022) |
| 7.10 | The authors state the loss is partly a **release** effect, not only a destruction effect: > "These interactions and reactions may not only lead to a decreased free 2-furfurylthiol in coffee brews, but also reduce its quantity, thus greatly **inhibit the release** of 2-furfurylthiol from the coffee brews into the air" — i.e. the measured headspace decline conflates **binding** and **degradation**. Directly relevant to the Sun 2019 reversible/irreversible split | [M] §3.2.2, p. 5 |
| 7.11 | 2,4-dimethylthiazole "decreased significantly during the **first 20 days**" at 40 and 50 °C | [M] §3.2.2, p. 5, Fig. 4A–B |
| 7.12 | Of the 13 correlated volatiles, **only 3** (FFT, methanethiol, 2E3MP) obeyed the first-order model with R² 0.97–0.99 — the other 10 did not, and **the paper never reports their fits or says why they failed** | [M] §3.3.1, p. 6 |

---

## 8. WHAT IS NOT IN THIS PAPER (declared gaps)

| missing | consequence |
|---|---|
| **Supplementary Tables S1–S3** (the 86-compound × condition concentration matrix) | **not on disk**; only the 3 modelled compounds have retrievable time courses. Obtain from https://doi.org/10.1016/j.lwt.2025.117469 if the full matrix is wanted |
| **Any pH measurement** | brew pH never measured; no pH axis at all. Cannot support any pH term |
| **Any water activity / RH / bean moisture during storage** | only "moisture < 12 %" *before roasting*. **A storage fade rate without a_w is not transferable to another matrix** |
| **Any oxygen specification during storage** | unknown whether the 40–60 °C hold was aerobic. Given §7.8's disulfide channel, this is the single most important missing variable |
| **Error bars on any kinetic quantity** | n = 3 stated, SD convention stated, **nothing printed**. No SE on any k, C₀ or Ea |
| **Sensory numbers** | radar plots only; no tabulated S-R scores, no scale range |
| **Absolute quantitation** | semi-quantitative peak-area ratio to 3-heptanone; the µg/L axes are 3-heptanone-equivalents |
| **More than one roast level** | single 180 °C / 7.5 min Robusta profile |
| **Ea printed by the authors** | the paper never writes down a single Ea value; only the regression slopes |

---

## NEW-PARAMETER TABLE (consolidated)

k units are **day⁻¹** [Z] (from the "Storage time (day)" abscissa; the paper never states them).
C₀ units are **µg/L of brew, 3-heptanone-equivalents** [Z] (from the Fig. 5A–C ordinate).
Ea is computed as |slope| × 8.314 kJ mol⁻¹ with x = 1000/T [Z] (the authors print no Ea at all).

| parameter | value | units (as printed / [Z] where inferred) | conditions | anchor | provenance |
|---|---|---|---|---|---|
| k, 2-furfurylthiol | **0.035** | day⁻¹ [Z] | 40 °C, whole Robusta beans, 30 d, dry box | Table 2, p. 8 | [F] |
| k, 2-furfurylthiol | **0.039** | day⁻¹ [Z] | 50 °C, ditto | Table 2, p. 8 | [F] |
| k, 2-furfurylthiol | **0.042** | day⁻¹ [Z] | 60 °C, ditto | Table 2, p. 8 | [F] |
| k, 2-furfurylthiol (body text variant) | **0.041** | day⁻¹ [Z] | 60 °C — **conflicts with Table 2's 0.042** | p. 6, §3.3.2 | [F] |
| C₀, 2-furfurylthiol | **12.621 / 12.500 / 12.348** | µg/L [Z] | 40 / 50 / 60 °C arms | Table 2, p. 8 | [F] |
| R², FFT first-order | **0.987 / 0.986 / 0.985** | — | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| t½, 2-furfurylthiol | **19.80 / 17.77 / 16.50** | day | 40 / 50 / 60 °C | from Table 2 k | [Z] |
| **Ea, 2-furfurylthiol — AS PRINTED** | **14.65** | kJ/mol [Z] (slope −1.762) | 40–60 °C | Table 2 lower block, p. 8 | [F] — ⚠️ **CORRUPT, is a fit to methanethiol data, see §6** |
| **Ea, 2-furfurylthiol — HONEST** | **7.92** | kJ/mol | 40–60 °C, refit of Table 2's own FFT k | §6.3/§6.4 | **[Z]** |
| k₀, 2-furfurylthiol — as printed | **4.250** | day⁻¹ [Z] (intercept +1.447) | — | Table 2, p. 8 | [F] — ⚠️ corrupt |
| k₀, 2-furfurylthiol — honest | **0.7365** | day⁻¹ | refit of Table 2's own FFT k (R² 0.9920) | §6.3 | **[Z]** |
| **k(25 °C), FFT — as printed** | **0.0102** | day⁻¹ [Z] (t½ 68.0 d) | extrapolated | p. 6, §3.3.2 | [F] — ⚠️ **~3× too slow** |
| **k(25 °C), FFT — honest** | **0.0302** | day⁻¹ (t½ 23.0 d) | extrapolated from the honest fit | §6.4 | **[Z]** |
| Q₁₀, FFT | **1.114 (40→50), 1.077 (50→60)** | — | — | from Table 2 k | **[Z]** |
| k, methanethiol | **0.015 / 0.018 / 0.021** | day⁻¹ [Z] | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| k, methanethiol (body text variant) | **0.015 / 0.017 / 0.019** | day⁻¹ [Z] | 40 / 50 / 60 °C — **conflicts with Table 2 at 50 and 60 °C** | p. 6, §3.3.2 | [F] |
| C₀, methanethiol | **6.202 / 5.976 / 6.03** | µg/L [Z] | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| R², methanethiol first-order | **0.984 / 0.964 / 0.970** | — | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| t½, methanethiol | **46.21 / 38.51 / 33.01** | day | 40 / 50 / 60 °C | from Table 2 k | [Z] |
| Ea, methanethiol — as printed | **10.25** | kJ/mol [Z] (slope −1.233) | 40–60 °C | Table 2, p. 8 | [F] — slope OK, **intercept −0.717 is a printing error, see §6.3** |
| Ea, methanethiol — from Table 2 triple | **14.60** | kJ/mol (k₀ 4.100 day⁻¹, R² 0.9991) | 40–60 °C | §6.3 | **[Z]** |
| k(25 °C), methanethiol | **0.0076** | day⁻¹ [Z] (t½ 91.2 d) | extrapolated; round-trips to +2.7 % | p. 6 | [F] |
| k, 2-ethyl-3-methyl-pyrazine | **0.018 / 0.0211 / 0.025** | day⁻¹ [Z] | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| k, 2E3MP (body text variant) | **0.017 / 0.021 / 0.025** | day⁻¹ [Z] | 40 / 50 / 60 °C — **conflicts with Table 2 at 40 and 50 °C** | p. 6, §3.3.2 | [F] |
| C₀, 2E3MP | **159.330 / 157.205 / 153.556** | µg/L [Z] | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| R², 2E3MP first-order | **0.992 / 0.953 / 0.979** | — | 40 / 50 / 60 °C | Table 2, p. 8 | [F] |
| t½, 2E3MP | **38.51 / 32.85 / 27.73** | day | 40 / 50 / 60 °C | from Table 2 k | [Z] |
| Ea, 2E3MP — as printed | **17.41** | kJ/mol [Z] (slope −2.094) | 40–60 °C | Table 2, p. 8 | [F] — **slope inconsistent with its own Fig. 5F, see §6.3** |
| Ea, 2E3MP — from Table 2 triple | **14.24** | kJ/mol (k₀ 4.253 day⁻¹, R² 0.9987) | 40–60 °C | §6.3 | **[Z]** |
| k(25 °C), 2E3MP | **0.0055** | day⁻¹ [Z] (t½ 126.0 d) | extrapolated; round-trips to −0.6 % | p. 6 | [F] |
| FFT 30-day loss | **64 / 71 / 72** | % | 40 / 50 / 60 °C | Fig. 5A | **[fig]/[Z]** |
| Methanethiol 30-day loss | **37 / 43 / 48** | % | 40 / 50 / 60 °C | Fig. 5B | **[fig]/[Z]** |
| 2E3MP 30-day loss | **40 / 46 / 51** | % | 40 / 50 / 60 °C | Fig. 5C | **[fig]/[Z]** |
| FFT concentration–time series | 21 points, §5.1 | µg/L | 40/50/60 °C × 7 days | Fig. 5A, p. 8 | **[fig]** ±0.15 |
| Methanethiol concentration–time series | 21 points, §5.2 | µg/L | 40/50/60 °C × 7 days | Fig. 5B, p. 8 | **[fig]** ±0.06 |
| 2E3MP concentration–time series | 21 points, §5.3 | µg/L | 40/50/60 °C × 7 days | Fig. 5C, p. 8 | **[fig]** ±2 |
| Odour threshold, FFT | **0.01** | µg/L | water | Table 1, p. 5 | **[C]** — cited, not measured here |
| Odour threshold, methanethiol | **0.20** | µg/L | water | Table 1, p. 5 | **[C]** |
| Odour threshold, 2-methyl-3-furanthiol | **0.05** | µg/L | water | Table 1, p. 5 | **[C]** |
| Odour threshold, dimethyl trisulfide | **0.001** | µg/L | water | Table 1, p. 5 | **[C]** |
| Odour threshold, dimethyl disulfide | **0.82** | µg/L | water | Table 1, p. 5 | **[C]** |
| Odour threshold, 5-methyl-2-furanmethanethiol | **0.05** | µg/L | water | Table 1, p. 5 | **[C]** |
| (remaining 30 thresholds) | see §4.0 | µg/L | water | Table 1, p. 5 | **[C]** |
| All 36 OAV values | see §4.0 | — | — | Table 1, p. 5 | [F] — ⚠️ **irreconcilable with Table 2 concentrations, §4.1; do not use** |
| PLS-R S-R correlation R² | **0.99 / 0.97 / 0.90** | — | 40 / 50 / 60 °C | Fig. 3B/D/F | [F] |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> Zhou 2025 is **not** in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. **A declaration amendment
> is required before any wave may fit it.** This section is a proposal only; the declaration was
> not edited.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 2 FFT first-order rows, 40 °C and 60 °C** (k = 0.035, 0.042 day⁻¹; C₀ 12.621, 12.348 µg/L) | **FIT** | temperature (end points) | The corpus has **no** thiol-loss constant anywhere near storage temperature; the nearest is Hofmann 2002's 80 °C brew. These two rows pin the low-T end of the thiol-consumption lane. Fit the **k values**, never the paper's Arrhenius row |
| **Table 2 FFT 50 °C row** (k = 0.039 day⁻¹) | **★ HOLD-OUT** | temperature (interpolation) | The only interior point on the FFT T-axis; scoring it tests whether the model's T-dependence has the right *curvature* between two fitted end points. Cheap and genuinely blind |
| **Fig. 5A FFT concentration–time series, all 21 points** | **HOLD-OUT** | — | `digitised_from_figure`, ±0.15 µg/L, and it is the *same information* as the Table 2 fits. Fitting both double-counts. Keep as a shape check on whether the decay is truly single-exponential (the 40 °C arm visibly flattens after day 20 — a first-order model should not) |
| **Table 2 methanethiol rows (all 3 T)** | **FIT** | — | An independent low-molecular-weight thiol on the same matrix and clock; constrains the shared oxidative-sink term without touching the FFT scoring axis |
| **Table 2 2-ethyl-3-methyl-pyrazine rows (all 3 T)** | **neither** | — | Pyrazine storage loss is outside the sulfur branch's scope and the repo has no pyrazine-persistence node. Ingest only as the qualitative observation that a non-thiol also fades with Ea ≈ 14 kJ/mol — i.e. **the fade is not thiol-specific**, which is itself evidence for the transport-limited reading in §6.5 |
| **The three printed Arrhenius formulas (Table 2 lower block)** | **★ neither — DO NOT INGEST AS PUBLISHED** | — | **The FFT row is a fit to methanethiol data (§6.2–6.3); the methanethiol intercept and the 2E3MP slope are each inconsistent with their own figures.** If an FFT Arrhenius parameter is wanted, use the **[Z] refit: Ea = 7.92 kJ/mol, k₀ = 0.7365 day⁻¹, R² = 0.9920**, and label it as a repo-side derivation from Zhou 2025 Table 2, not as a Zhou 2025 result |
| **The 25 °C extrapolations (k = 0.0102 / 0.0076 / 0.0055 day⁻¹)** | **neither** | — | The FFT value carries the full 3× error of §6.4. The methanethiol and 2E3MP values are sound but are extrapolations, not measurements — they add no information beyond the k values already proposed for FIT |
| **Table 1 OAV column** | **★ neither — REJECT** | — | Irreconcilable with the paper's own concentrations by 18.5× / 102.7× / 10.0× (§4.1), and the sulfury/roasty group total is arithmetically wrong by 980. **Contamination hazard if ingested** |
| **Table 1 threshold column** | **neither (route to primaries)** | — | All **[C]**. The repo already holds Czerny 2008 and Belitz 2009 thresholds; take them from the primaries rather than from this secondary transcription |
| **§7.8 dimethyl disulfide as a thiol-oxidation marker** | **qualitative ingest only** | — | Supports a thiol → disulfide channel active at 40–60 °C in the dry bean. No number attached; the DMDS time course is in Fig. 4 only, which is not digitised here |
| **§7.10 binding-vs-destruction conflation** | **qualitative ingest, HIGH PRIORITY** | — | The authors state explicitly that measured headspace decline conflates reduced release with reduced quantity. **This is the same distinction Sun 2019 operationalises**, and it is the reason the k values above must be modelled as *apparent* disappearance, not as chemical consumption |

### ⚠️ Circularity and honesty warnings the orchestrator must carry forward

1. **Do not let this dataset certify the sulfur branch's temperature dependence.** Ea ≈ 8 kJ/mol
   over 40–60 °C is not commensurable with the ~120 kJ/mol regime the roasting-temperature
   anchors (Zheng 1994) live in. A model that reproduces both is either genuinely
   two-mechanism or is being fit to a transport artefact. **If the repo's thiol-consumption lane
   is currently a single Arrhenius law, these rows will break it — and that break is the
   informative result, not a failure to be tuned away.**
2. **The 40 °C and 60 °C FIT rows and the 50 °C HOLD-OUT row come from one experiment on one
   roast of one Robusta lot with n = 3 analytical (not biological) replicates and no printed
   error bars.** Any sigma assigned to them is a judgement call, not a measurement. Suggest
   propagating the **C₀ spread (±1.1 % for FFT)** as the floor and inflating it substantially for
   the unreported a_w and O₂.
3. **The paper is a shelf-life-prediction paper, not a mechanism paper.** Its k values are
   defensible as *what happens*; none of its Arrhenius parameters are defensible as *why*.
