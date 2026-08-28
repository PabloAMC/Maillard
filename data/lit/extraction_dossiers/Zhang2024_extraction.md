# Zhang, Kuang, Wang & Cao 2024 (10.1016/j.foodres.2024.114149) — Wave Z3 extraction 2026-08-28

Source: `data/articles/Zhang2024.pdf`, read in full (10 pp). *Food Research International*
**182 (2024) 114149**. "Metabolomics reveals factors affecting the radical reaction of
sulfides during thermal processing for meaty aroma." Beijing Technology and Business
University (BTBU). Received 30 Nov 2023; accepted 17 Feb 2024.

All figure panels re-read from **300 dpi page renders** (`z3img/zh-0{4,5,7}.png` and crops
`f2a..f2g`, `f4a`, `f4c`, `zh04_fig1a`). Table 1 matched the text layer exactly.

---

## ⚠ READ THIS BEFORE THE BRIEF'S VERDICT SECTION — WHAT THIS PAPER IS AND IS NOT

The brief nominates this as **THE priority test** of the hypothesis that Wave U's backwards
sulfur temperature dependence is explained by **MFT consumption at processing temperature**.

**The paper supplies the MECHANISM and its magnitudes. It does NOT supply a temperature
series for MFT, and it contains no activation energy, no Arrhenius fit, and no rate constant
for any MFT-consuming step.** Every kinetic number in it is at **115 °C**, one temperature.

Worse for direct use: **every concentration–time profile in this paper lives in the
Supplementary Information** (`Fig. S1`, `Fig. S2`, `Fig. S3`, `Table S1`–`Table S7`), and the
SI is **NOT in the retrieved PDF**. The dynamic profiles the brief asks for (Fig. S1),
the single-factor temperature sweep (Fig. S2), the Box–Behnken design and its responses
(Tables S2/S3), and the zero-order regression table (Tables S6/S7) are all absent.
**Only the summary numbers quoted in the running text survive.** They are transcribed below,
each with its page anchor, and they are the whole of the paper's kinetics.

**What IS here, and it is not nothing:** a measured, quantified, three-channel MFT sink at
processing temperature, plus a *direction* statement about temperature from the paper's own
optimisation. See §7.

---

## 1. CONDITIONS — verbatim, all four systems

### 1.1 "Simple models" (§2.2 Preparation of MFT, p. 2) — the Fig. 1 system
> "The binary mixture of **VB1 at 15.0 mg/mL (w/v) and Xyl in a ratio 1:1** was solubilized in
> **pH 4.9 ± 0.1 phosphate buffered solution**. The mixture was made to react for **60 min**,
> then immediately cooled in ice water. … The treatment groups were samples prepared at the
> same thermal procedure with adding **Cys, GSH and GCys of equal weight as VB1 and Xyl**,
> respectively."

⚠ **The temperature of this step is NOT STATED in §2.2.** Everywhere else in the paper the
"same thermal process" is 115 °C / 60 min (§2.3, §2.5, §2.6). **Treat Fig. 1 as 115 °C only
with that caveat attached; the paper does not say so.**

VB1 = thiamine, Xyl = xylose, Cys = cysteine, GSH = glutathione, GCys = cystine.
15.0 mg/mL thiamine hydrochloride (MW 337.27) ≈ **44.5 mM**; xylose (MW 150.13) at equal
weight ≈ **99.9 mM**; Cys (121.16) ≈ 123.8 mM; GSH (307.32) ≈ 48.8 mM; GCys (240.30) ≈ 62.4 mM.
*(These molarities are MY conversion from the paper's w/v figures, not the paper's.)*

### 1.2 "Multi-component models" (§2.3, p. 2) — the Fig. 2 system
> "**Met (15.0 mg/mL, w/v)** as the source of MeSH mixed with **VB1 and Xyl as a ratio of
> 1:1:1**, which dissolved in phosphate buffered solution (**pH 4.9 ± 0.1**). For treatment
> groups, Cys (15.0 mg/mL, w/v), GSH (15.0 mg/mL, w/v), and GCys (15.0 mg/mL, w/v) were mixed
> into pH 4.9 phosphate buffered solution, respectively. The mixture solutions were incubated
> for up to **60 min at 115 °C**. Subsequently, mixtures were cooled with ice water."

⚠ **BUT §2.6 (p. 3) uses a DIFFERENT ratio for the Fig. 2 dose–response:**
> "The **ternary mixture of Met, VB1 and Xyl (1:3:3, w/w)** mixed with a series of Cys (or GSH)
> solutions with different concentrations (**0.5, 2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 62.5
> and 75.0 mg/mL, w/v**), respectively. The mixture was incubated at **115 °C for 60 min**…"

and §2.5 states the 1:3:3 ratio is "**approximately equal to the ratio in real thermal process
flavorings sample**". **So Fig. 1 = 1:1 VB1:Xyl; Fig. 2 = 1:3:3 Met:VB1:Xyl. Two different
stoichiometries. Do not merge them.**

### 1.3 Thermal process flavourings (§2.4, p. 2) — the Fig. 4 system, a REAL FOOD
> "**Longissimus thoracic muscles** (trimmed of subcutaneous fat) mixed with pre-weighted
> ingredients … The mixture was heated with a **115 ℃ conventional thermal process for
> 60 min**. … In the treatment group … With the addition of **Cys (2.5 mg/mL) and GSH
> (6.0 mg/mL)**, the mixture was prepared at **115 ℃ for 50 min and then at 95 ℃ for 10 min**."

CHPBF = control; THPBF = treatment. **Pork muscle matrix, not a model solution.**

### 1.4 Quantification (§2.10, p. 3) — read this before converting any ppb
> "Samples of **10 mL** mixed with **1 μL IPDS (0.0943 ng/mL, w/v)** were sealed in a 40 mL
> headspace glass … 30 min incubation at 45 ℃ … absorbed with a manual SPME holder with a
> 50/35 μm DVB/CAR/PDMS fiber (2 cm) at **45 °C for 30 min** … GC–MS/SIM … splitless."
> "The **concentrations (ng/mL) of identified each above-mentioned individual volatile were
> through construction of standard curves (Table S5)**. The corresponding **external standard
> curves** were obtained by mass ratio and peak area ratios of eight different concentrations
> of standard solutions with the internal standard ion peak. Furthermore, the
> **semiquantitative data** of volatile compounds **for the verification test** were estimated
> by comparing the peak area of each chemical with that of the internal standard."

⚠ Two quality tiers in one paper: Figs. 1/2 are **external-standard-curve quantified (ng/mL)**;
Fig. 4c is explicitly **semiquantitative**. Fig. 4a's DMTS is in **µg/mL**, Fig. 4c's in
**ng/mL** — a 1000× unit change between two panels of the same figure. Triplicate; Duncan's
test, p < 0.05.

⚠ **HS-SPME with an external standard curve is not isotope dilution.** These ng/mL values are
headspace-partition-dependent and matrix-dependent. They are **not commensurable with
Hofmann's SIDA ppb**, and any repo row that puts them side by side must say so.

---

## 2. FIGURE 1 (p. 4) — the simple VB1/Xyl model. **300 dpi read-off, flagged.**

"Impact of Cys, GSH and GCys on **concentrations (a)** and **percentage (b)** of MFT, MMFT and
MFT-MFT in the simple models."

### 2a — concentration, ng/mL (error bars drawn on every bar)
| group | MMFT | **MFT** | MFT-MFT |
|---|---:|---:|---:|
| VB1-Xyl (control) | 0.035 | **0.07** | 0.04 |
| VB1-Xyl-**Cys** | 0.04 | **1.34** | 0.115 |
| VB1-Xyl-**GSH** | 0.035 | **1.05** | 0.07 |
| VB1-Xyl-**GCys** | 0.04 | **2.01** | 1.09 |

### 2b — percentage of the three-species sum, %
| group | MMFT | MFT | MFT-MFT |
|---|---:|---:|---:|
| VB1-Xyl | 24 | 49 | 28 |
| VB1-Xyl-Cys | 3 | 90 | 7.5 |
| VB1-Xyl-GSH | 3.5 | 91 | 6.5 |
| VB1-Xyl-GCys | 1.5 | 64 | 34.5 |

**MY READ-OFF IS SELF-VALIDATING.** Panel (b) is the normalisation of panel (a), and the two
were read independently. Recomputing (b) from my (a) numbers:
control 0.07/0.145 = **48.3 %** (read 49); Cys 1.34/1.495 = **89.6 %** (read 90);
GSH 1.05/1.155 = **90.9 %** (read 91); GCys 2.01/3.14 = **64.0 %** (read 64).
**All four within 1 percentage point.** The panel-(a) values are therefore good to ~±3 %.

### The number that matters
**MFT-MFT (the oxidative dimer) is 54.2 % of the MFT pool by mass in the GCys system
(1.09 vs 2.01) and only 8.6 % in the Cys system (0.115 vs 1.34) — a 6.3× swing driven by
nothing but which sulfur additive is present, at ONE temperature and ONE time.**
The paper's own words (p. 4): *"the GCys group exhibited the lowest MFT percentage (Fig. 1),
which was attributed to **MFT being readily oxidized by free radicals, forming MFT-MFT, its
dimer**."*

---

## 3. FIGURE 2 (p. 5) — dose–response in the multi-component model, 115 °C / 60 min, pH 4.9

Ten additive levels: **0.5, 2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 62.5, 75.0 mg/mL**.
All values ng/mL. **300 dpi read-off from individual panel crops; ±3–5 % of full scale.**

### (a) MeSH
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 7.3 | 23.2 | 36.6 | 43.4 | 48.8 | 52.9 | **71.5** | 61.5 | 56.5 | 47.3 |
| GSH | 12.0 | 18.9 | 18.7 | 25.5 | 46.5 | 54.0 | **83.5** | 68.7 | 58.2 | 60.0 |

### (b) DMDS — the only monotone panel
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 0.0755 | 0.0648 | 0.0578 | 0.0522 | 0.0511 | 0.0473 | 0.0454 | 0.0424 | 0.0421 | 0.0421 |
| GSH | 0.0700 | 0.0596 | 0.0637 | 0.0523 | 0.0494 | 0.0456 | 0.0440 | 0.0439 | 0.0412 | 0.0410 |

### (c) DMTS
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 0.163 | 0.286 | 0.333 | 0.369 | 0.427 | 0.428 | 0.470 | **0.472** | 0.254 | 0.208 |
| GSH | 0.112 | 0.166 | 0.148 | 0.161 | **0.241** | 0.205 | 0.222 | 0.233 | 0.211 | 0.121 |

### (d) **MFT** — the panel the brief needs
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 0.4 | 0.75 | 1.35 | 1.6 | 2.4 | 2.5 | 4.0 | 5.7 | **11.0** | 11.0 |
| GSH | 0.2 | 0.25 | 0.9 | 1.4 | 4.35 | 6.4 | 12.2 | **23.6** | 19.6 | 18.2 |

**MFT in the GSH arm is NON-MONOTONIC: it peaks at 50 mg/mL (23.6 ng/mL) and falls 1.30× by
75 mg/mL.** The Cys arm plateaus flat at 11.0 from 62.5 mg/mL. Neither is a saturating
Michaelis curve; the GSH arm turns over.

### (e) MMFT (= MFT + MeSH radical product)
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 0.050 | 0.085 | 0.090 | 0.089 | 0.108 | 0.102 | **0.195** | 0.194 | 0.174 | 0.142 |
| GSH | 0.039 | 0.042 | 0.048 | 0.056 | 0.121 | 0.164 | 0.171 | **0.234** | 0.204 | 0.202 |

### (f) MFT-MFT (= MFT oxidative dimer)
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 0.02 | 0.03 | 0.04 | 0.04 | 0.05 | 0.10 | 0.55 | **1.00** | 0.68 | 0.55 |
| GSH | 0.02 | 0.05 | 0.05 | 0.06 | 0.40 | 1.28 | 2.07 | **5.75** | 4.65 | 4.85 |

### (g) MMFT / MFT-MFT ratio — the paper's own response variable R2
| mg/mL | 0.5 | 2.5 | 5 | 10 | 20 | 30 | 40 | 50 | 62.5 | 75 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cys | 1.80 | **1.97** | 1.66 | 1.42 | 1.90 | 1.02 | 0.32 | 0.18 | 0.24 | 0.23 |
| GSH | 1.29 | **1.40** | 1.25 | 0.83 | 0.32 | 0.13 | 0.08 | 0.04 | 0.045 | 0.04 |

**This ratio falls ~11× (Cys) and ~35× (GSH) across the additive range.** The paper (p. 5):
*"the concentration ratios (Fig. 2 g) decreased and the values in the Cys group were
significantly larger than those in the GSH group, which indicated that MFT in the Cys group
could readily react with MeSH, leading to the generation of MMFT rather than self-oxidation
and the formation of disulfides."*

### THE FLUX ARITHMETIC — MFT'S TWO MEASURED SINKS, SIDE BY SIDE (my derivation from d/e/f)
At each additive level, dimer-sink : radical-sink : surviving free MFT, by mass:

| mg/mL | GSH: MFT-MFT | MMFT | MFT | **(MMFT+2·dimer_MFT_equiv) / MFT** |
|---|---:|---:|---:|---:|
| 20 | 0.40 | 0.121 | 4.35 | ~0.21 |
| 30 | 1.28 | 0.164 | 6.4 | ~0.44 |
| 40 | 2.07 | 0.171 | 12.2 | ~0.36 |
| 50 | 5.75 | 0.234 | 23.6 | ~0.51 |
| 75 | 4.85 | 0.202 | 18.2 | ~0.55 |

*(dimer carries 2 MFT units; MW MFT-MFT 226.32, MFT 114.17, so 1 ng dimer ≈ 1.009 ng MFT.
`(MMFT_as_MFT + 2·dimer_as_MFT)/MFT` with MMFT MW 160.25 → ×0.7125.)*
**In the GSH arm, MFT that has already been consumed into its two measured products is
21–55 % of the MFT still free, at 115 °C, in 60 min.** In the Cys arm the same ratio runs
~0.03–0.19. **These are lower bounds on total MFT consumption** — the paper names a third
sink (melanoidin/quinone ionic reactions, §Introduction) that it never measures at all.

---

## 4. THE ONLY KINETIC CONSTANTS IN THE PAPER (§3.1, p. 4, verbatim)

> "The dynamic analysis results (**Table S6**) showed that the **zero-order kinetic equations
> of MMFT levels** in the Cys and GSH treatment groups had better-fitted data. The two groups'
> respective **reaction rate constants were 0.0028 and 0.0031**, which indicated there was not
> much of a difference between either. Subsequently, **linear regression (shown in Table S7)**
> was applied to analyze the effects of Cys and GSH on the rate of changes (i.e., the slope of
> the regression line) in MMFT levels, as the two treatments and control group showed
> interesting **two-stage trends before and after 90 min** of thermal processing. Results
> showed that during the **first 90 min** of thermal processing, the **slope of MMFT level in
> Cys and GSH groups reached 0.0043 and 0.0033**, which are about **3 times and 2.5 times the
> slope in the control group**, respectively. After that, the **slope in the Cys group dropped
> to a comparable level to that in the control group**. The **slope of MMFT in the GSH
> treatment group is 0.0032, which is 4 times that of the control group**."

| quantity | Cys | GSH | control (derived) |
|---|---:|---:|---:|
| zero-order k, MMFT, whole run | **0.0028** | **0.0031** | not given |
| slope, MMFT, **t < 90 min** | **0.0043** | **0.0033** | 0.00143 / 0.00132 → ≈ **0.0014** |
| slope, MMFT, **t > 90 min** | ≈ control | **0.0032** | 0.0032/4 = **0.0008** |

⚠ **UNITS ARE NEVER STATED.** MMFT is in ng/mL and time in min throughout, so the natural
reading is **ng·mL⁻¹·min⁻¹**, but the paper does not say so and Table S6/S7 are not in the
retrieved PDF. **Do not ingest these as numbers without the unit resolved.**
⚠ The control slope is MY back-calculation from "about 3 times" and "2.5 times" — the two
routes give 0.00143 and 0.00132, consistent, but it is a division of a rounded ratio.
⚠ **The thermal run in §2.3 is "up to 60 min". These slopes are quoted over 90+ min.** The
kinetic run behind Table S6/S7 is therefore a LONGER experiment than §2.3 describes, and the
PDF never states its duration or its temperature. Flag this before any use.

**The one structural fact that survives all of that: the MMFT slope BREAKS at 90 min and, for
Cys, falls back to the control rate.** The paper: *"After that, time may take over, making the
impact of additives less dominant."* A one-stage rate law cannot express this.

---

## 5. FIGURE 4 (p. 7) — verification in the pork thermal process flavouring

### 4a — response variables. **Twin-axis figure; read at 300 dpi; the right axis is 0–3.**
| group | DMTS (µg/mL) | MMFT / MFT-MFT |
|---|---:|---:|
| **CHPBF** (control, 115 °C / 60 min) | **3.60 ± 0.53** | **≈ 1.09** |
| **THPBF** (treatment, 115 °C/50 min → 95 °C/10 min, + Cys 2.5 + GSH 6.0 mg/mL) | **1.70 ± 0.15** | **≈ 2.71** |

⚠ The MMFT/MFT-MFT values are read against the *secondary* axis (0–3, with "3" aligned to
"4.45" on the primary axis). **±0.15 at best.** Direction is safe: **DMTS ×0.47, ratio ×2.5.**

### 4c — DMTS during 112 days of storage, **semiquantitative ng/mL**
| day | 1 | 28 | 56 | 84 | 112 |
|---|---:|---:|---:|---:|---:|
| CHPBF | 2.2 | 6.3 | 16.1 | 10.7 | 17.6 |
| THPBF | *(no bar drawn)* | *(no bar drawn)* | 2.35 | 6.5 | 6.9 |

⚠ **No THPBF bar is drawn at days 1 and 28.** The paper never says why (below LOD? not
measured?). Do not read those as zeros.
⚠ CHPBF is **non-monotonic** (16.1 at d56 → 10.7 at d84 → 17.6 at d112) — inconsistent with a
single accumulation process and unremarked by the authors.

### 4b — aroma profile, radar, 1–9 scale (read: meaty CHPBF ≈ 5.0, THPBF ≈ 6.1; treatment
lower on burnt and spicy). **Ordinal only.** Text (p. 7): *"treatment group with the addition
of Cys and GSH has higher score of meaty intensity and lower burnt and spicy intensity."*

---

## 6. THE OTHER QUANTITATIVE STATEMENTS IN THE RUNNING TEXT

**RSM (§3.4.1, p. 6):**
- R² = **0.905** for DMTS (R1) and **0.873** for MMFT/MFT-MFT (R2); model p-values **0.0004**
  and **0.0019**; lack of fit not significant. 27 experimental points, Box–Behnken, 4 factors
  × 3 levels, 3 centre replicates.
- **Optimum: Cys:GSH = 3:7, initial reaction temperature 115 ℃, second reaction time 10 min,
  second stage temperature 95 ℃.**
- *"reaction time and reaction temperature had a relatively dominant effect on flavor
  compounds compared to the concentration ratio of Cys and GSH."*

**Verification of the optimum in the MODEL (§3.4.2, p. 6):**
- measured MMFT/MFT-MFT **1.422** vs predicted **1.294** (9.9 % high);
- DMTS: *"the difference between the actual measured DMTS concentrations and the predicted
  value predicted is **7.0 %**."*

**Sensory threshold (Table 1, p. 5).** Forced-choice ascending, 12 trained panelists, ratios of
Cys : Met concentration = 0.1 / 0.5 / 1 / 2 / 4 / 6 / 8 / 10 / 12.5 / 15 / 17.5 / 20.
Full 24-row × 12-column ±-matrix transcribed below.

| group | panel | 0.1 | 0.5 | 1 | 2 | 4 | 6 | 8 | 10 | 12.5 | 15 | 17.5 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Model | 1 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 2 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 3 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 4 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 5 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 6 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 7 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 8 | – | – | – | – | + | + | + | + | + | + | + | + |
| | 9 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 10 | – | – | – | – | + | + | + | + | + | + | + | + |
| | 11 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 12 | – | – | – | – | – | + | + | + | + | + | + | + |
| Thermal process flavorings | 1 | – | – | – | – | – | – | – | + | + | + | + | + |
| | 2 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 3 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 4 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 5 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 6 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 7 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 8 | – | – | – | – | – | + | + | + | + | + | + | + |
| | 9 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 10 | – | – | – | – | – | – | + | + | + | + | + | + |
| | 11 | – | – | – | – | – | – | – | + | + | + | + | + |
| | 12 | – | – | – | – | – | – | + | + | + | + | + | + |

("+" = difference detectable by the panelist; "–" = not detectable.)
⇒ *"the differences among all samples could be distinguished up until the concentration ratio
of Cys and Met … reached **6** in the ternary mixture and **8** in thermal process flavorings.
… the Cys/Met concentration ratio should be **smaller than 6** (corresponds to the optimal Cys
used level of **30.0 mg/mL** in the model)."* **Sensory, not chemical — do not ingest as a
concentration constraint.**

**Metabolomics (§3.5, pp. 7–8).** OPLS-DA 84 % prediction accuracy, 2 PCs; **232 significant
metabolites** (VIP > 1.0 and p < 0.05) assigned to **86 pathways**. Down-regulation of thiamin
monophosphate; up-regulation of Cys. **Nothing here is a rate or a concentration.**

---

## 7. ★ THE PAPER'S THREE TEMPERATURE STATEMENTS — the entire temperature content

All three are in §3.4.1 (p. 6), all refer to **Fig. S2, which is not in the retrieved PDF**:

1. > "A **lower initial reaction temperature** could **significantly reduce the concentrations
   > of DMTS** and **promote an increase in the MMFT/MFT-MFT concentration ratio**."
2. > "DMTS had a lower concentration when the reaction time was **less than 20 min** or the
   > **second stage temperature was less than 105 ℃**, meanwhile, the concentration ratio of
   > MMFT/MFT-MFT could reach the highest value."
3. The RSM optimum itself: **the second stage is 95 °C, i.e. 20 °C BELOW the first stage.**
   A four-factor optimisation over a designed space chose to **drop** the temperature for the
   final 10 min in order to maximise MMFT/MFT-MFT and minimise DMTS.

**These are direction-of-temperature statements about DMTS and about a RATIO. Neither is a
statement about the absolute MFT concentration, and no MFT-vs-temperature number exists
anywhere in this paper.**

---

## ★★ THE VERDICT THE BRIEF REQUIRES

> **Does Zhang 2024 explain a NET DECREASE of MFT with rising temperature, and with what
> numbers?**
>
> **NO — not as a measurement. It cannot: it is a single-temperature (115 °C) study whose
> temperature sweep is in an SI that is not in the retrieved file. There is no number in this
> paper that can be put on a temperature axis for MFT.**
>
> **What it DOES establish, with numbers, is the PRECONDITION for a net decrease: that MFT is
> consumed, fast, at processing temperature, by at least two measured channels and one
> named-but-unmeasured third.**
> - MFT → **MFT-MFT** (oxidative dimerisation): up to **5.75 ng/mL of dimer against 23.6 ng/mL
>   of free MFT** (GSH, 50 mg/mL) — the dimer alone carries **~49 % as much MFT as remains
>   free** (Fig. 2d/2f).
> - MFT + MeSH → **MMFT** (radical coupling): the paper's whole thesis, and the only channel it
>   fits a rate to (**zero order, k = 0.0028 / 0.0031**).
> - MFT + **melanoidins/quinones** (ionic): named in the Introduction (p. 2, citing Tang 2013
>   and Charles-Bernard 2005) and **never measured**.
>
> **And it establishes the DIRECTION, from the paper's own optimisation rather than from a
> kinetic series: lowering the temperature of the final processing stage from 115 °C to 95 °C
> IMPROVES the meaty-thiol outcome** — DMTS ×0.47 and MMFT/MFT-MFT ×2.5 in a real pork matrix
> (Fig. 4a). That is a designed-experiment result, it is in the right direction for Wave U, and
> it is a *ratio and a co-product*, not MFT itself.
>
> **The honest one-line answer: this paper supplies the mechanism and the sign, not the
> magnitude. It cannot close Wave U's ×0.249, and any wave that claims it does is
> over-reading the SI it does not have.**

---

## WHAT THIS UNLOCKS

### 1. It converts "the network has no MFT sink" from an inference into a two-source finding
Wave Z2 established from van Seeventer 2001 that MFT is destroyed at **59 %/day at 50 °C** by
a **non-oxidative** electrophilic oligomerisation. This paper adds, at **115 °C**, an
**oxidative** dimerisation channel that is **large enough to see in the mass balance** —
and the two papers do not conflict, because van Seeventer's system had residual cysteine and
an anti-oxidative matrix (air ≈ argon), while Zhang's GCys arm is deliberately oxidising.
**Together: three named MFT sinks (oligomerisation, dimerisation, melanoidin adduction), of
which the repository models zero.** That is the defect, and it is structural, not a barrier.

### 2. A measured branching ratio between two MFT sinks, at ten additive levels
`MMFT / MFT-MFT` is exactly a **branch-ratio observable** of the kind the repo's step-level
panel is short of, and it is a *ratio within one chromatogram*, so it is immune to the SPME
response-factor caveat that disqualifies the absolute ng/mL. Ten points × two additives = 20
ingestable ratio rows, plus one in a real matrix (Fig. 4a).

### 3. The 90-minute rate break
A single Arrhenius rate cannot produce a slope that is 3× control for 90 min and then equals
control. Whatever the units turn out to be, **the shape is a two-regime observable** and it is
the repo's only evidence of one on the sulfur branch.

### 4. The thiamine lane gets its first quantitative context
Hofmann 1998 Table 10 ranks thiamin **30× weaker than the C₂+C₃ route** for MFT (8.2 µg vs
268.1 µg). Zhang's whole system is thiamine/xylose and it makes MFT at 0.07–23.6 ng/mL. These
are not commensurable (SIDA µg per 50 mL vs SPME ng/mL), but the paper does confirm the
mechanism the repo would need: thiamine → **5-hydroxy-3-mercapto-2-pentanone** → MFT, citing
**Cerny & Guntz-Dubini 2008 (10.1021/jf801762c)** — a paper the repo does not hold and which is
the direct primary source for that intermediate. **[RETRIEVE]**

### 5. A negative result worth recording: additive effects are NON-MONOTONIC on every panel
MeSH peaks at 40 mg/mL, DMTS peaks at 50 (Cys) / 20 (GSH), MFT peaks at 50 (GSH), MMFT peaks
at 40–50, MFT-MFT peaks at 50. **Six of seven panels turn over.** Only DMDS is monotone. Any
repo term that makes a sulfur output monotone in cysteine loading is contradicted by six
independent curves in one figure — and this corroborates van Seeventer's sign-crossing
cysteine minimum (Z2 §Figure 1) from a different lab, different system, different temperature.

---

## INGESTION RECOMMENDATION

| item | class | note |
|---|---|---|
| `MMFT/MFT-MFT` vs additive concentration, 10 pts × 2 additives (Fig. 2g) | **benchmark set (branch ratio)** | ratio within one run; the safest thing in the paper |
| `MFT-MFT / MFT` and `MMFT / MFT` at each additive level (Figs. 2d/e/f) | **benchmark set (sink fraction)** | derived ratios; state the derivation |
| Fig. 1 percentage panel (24/49/28 → 3/90/7.5 → 3.5/91/6.5 → 1.5/64/34.5) | **benchmark (composition)** | already normalised by the authors; ⚠ temperature not stated in §2.2 |
| CHPBF→THPBF: DMTS ×0.47, MMFT/MFT-MFT ×2.5 (Fig. 4a) | **directional claim, real matrix** | pork muscle; not a model system |
| "lower initial reaction temperature reduces DMTS and raises MMFT/MFT-MFT" | **directional claim (temperature)** | ⚠ underlying Fig. S2 not in the retrieved PDF |
| RSM optimum = second stage at **95 °C**, 20 °C below stage 1 | **directional claim** | the paper's strongest temperature statement |
| Six-of-seven non-monotone dose responses | **structural claim, sign-crossing** | falsifies any monotone cysteine term |
| MMFT slope breaks at 90 min; Cys falls back to control | **structural claim (two-regime)** | |
| Zero-order k = 0.0028 (Cys), 0.0031 (GSH); slopes 0.0043 / 0.0033 / 0.0032 | **DO NOT INGEST AS NUMBERS** | units unstated, Table S6/S7 absent, run duration inconsistent with §2.3 |
| Absolute ng/mL from Figs. 1, 2 | **DO NOT use as ppb anchors** | HS-SPME + external curve; not commensurable with SIDA |
| Fig. 4c storage series | **semiquantitative, non-monotone, two missing bars** | record, do not gate |
| Any Ea, any temperature dependence of MFT | **NOT INGESTABLE — single temperature (115 °C)** | |
| Cerny & Guntz-Dubini 2008, 10.1021/jf801762c | **[RETRIEVE]** | primary source for 5-hydroxy-3-mercapto-2-pentanone, the thiamine→MFT intermediate |
| The paper's own **Supplementary Information** | **[RETRIEVE]** | Fig. S2 (temperature single-factor), Tables S6/S7 (the kinetics), Table S5 (standard curves) — **the SI holds everything the brief actually asked for** |

---
---

# ★ WAVE K3 ADDENDUM (2026-08-28) — completion of the interrupted Z3 write, and the verdict K3 asks for

The Z3 predecessor was terminated mid-write. Everything above §"INGESTION RECOMMENDATION" was
already complete and has been **verified against the source text** (`z3_zhang2024.txt`,
78 322 bytes, the full converted PDF) and the page renders in `z3img/`. Nothing above is
changed. This addendum adds the four items the K3 brief names and that the Z3 write did not
reach, plus one cross-paper link that only became available this wave.

## K3.1 — CHECKLIST AGAINST THE K3 BRIEF

| K3 brief item | where it is | status |
|---|---|---|
| every concentration–time profile | §3 (Fig. 2 dose–response, 10 levels × 7 panels × 2 additives) | ⚠️ **PARTIAL, AND THE LIMIT IS THE PAPER'S, NOT THE ANALYST'S.** See K3.2 |
| the zero-order rate constants 0.0028 / 0.0031, with UNITS and exact conditions | §4 | **units NOT RECOVERABLE — see K3.3** |
| the two-stage temperature protocol (115 °C then 95 °C) | §1.3, §5, §6, §7 | **complete — see K3.4 for the consolidated statement** |
| cysteine / glutathione / **cystine** effects | §2 (all three, Fig. 1), §3 (Cys and GSH only, Fig. 2) | **complete — consolidated in K3.5** |
| explicit verdict: does it quantify MFT formation **and** loss? | — | **K3.6, new** |
| what it licenses that K1's 30 °C rate constant cannot | — | **K3.7, new — this is the load-bearing section** |

## K3.2 — THE CONCENTRATION–TIME PROFILES: WHAT EXISTS AND WHAT DOES NOT

**There is no concentration–time profile anywhere in the retrieved PDF.** Every axis in every
figure of the main article is a **dose** axis (mg/mL of additive), a **group** axis (control /
Cys / GSH / GCys), or a **storage-day** axis. The word "dynamic" appears in the paper only as a
pointer out of it:

> §3.1, p. 4: "The **dynamic analysis results (Table S6)** showed that the zero-order kinetic
> equations of MMFT levels in the Cys and GSH treatment groups had better-fitted data."

Enumerated, the paper's own **time-resolved** content is:

| what | where the paper puts it | in the retrieved PDF? |
|---|---|---|
| dynamic concentration profiles of MeSH/DMDS/DMTS/MFT/MMFT/MFT-MFT vs time | **Fig. S1** | **NO** |
| single-factor temperature sweep | **Fig. S2** | **NO** |
| Box–Behnken design matrix and its 27 responses | **Tables S2 / S3** | **NO** |
| zero-order regression table (the source of 0.0028 / 0.0031) | **Table S6** | **NO** |
| two-stage linear regression (the source of 0.0043 / 0.0033 / 0.0032) | **Table S7** | **NO** |
| standard curves | **Table S5** | **NO** |

**The ONLY time axis in the retrieved article is Fig. 4c, the 112-day storage series** (§5),
which is **semiquantitative** by the paper's own §2.10 and has **two missing bars**. It is a
storage-stability series at ambient, not a thermal-processing profile.

⇒ **The brief's request for "every concentration-time profile" is satisfied by reporting that
the paper contains none in the retrieved file, and by naming exactly which SI object holds each
one.** The dose–response tables in §3 above are the paper's entire quantitative surface, and
they are a **dose** surface. Requesting the SI (`10.1016/j.foodres.2024.114149`, Appendix A) is
the single highest-value retrieval attached to this paper and would convert **0 time series into
6 time series × 3 groups**, plus the temperature sweep the whole sulfur-temperature question
turns on.

## K3.3 — ★ THE ZERO-ORDER CONSTANTS: THE UNITS CANNOT BE RECOVERED, AND HERE IS THE PROOF

The Z3 write flagged the units as unstated. **This wave attempted the recovery and it fails.
Recording the failure so no later wave re-attempts it.**

The verbatim sentence is quoted in full at §4. The four candidate readings and why each fails:

| candidate unit | test | verdict |
|---|---|---|
| **ng·mL⁻¹·min⁻¹** | MMFT in Fig. 2e spans 0.039 → 0.234 ng/mL. A zero-order rate of 0.0028 ng·mL⁻¹·min⁻¹ over the §2.3 run of 60 min gives **0.168 ng/mL** — the right order. | **plausible, and the natural reading**, but see the duration objection below |
| **µg·mL⁻¹·min⁻¹** | would give 168 ng/mL, **~700× above every MMFT value in the paper** | **refuted** |
| **fraction·min⁻¹** (normalised) | would give 0.168 = 17 % conversion; nothing in the paper is expressed as a fraction | **unsupported** |
| **ng·mL⁻¹·h⁻¹** | 0.0028 × 1.5 h = 0.0042 ng/mL, **~50× below the observed MMFT range** | **refuted** |

**But the plausible reading does not close either**, and this is the part that blocks ingestion:
the slopes are quoted over a **two-stage trend "before and after 90 min"** while §2.3 describes
the multi-component model as incubated **"for up to 60 min at 115 °C"**. A regression that has a
break at 90 min cannot come from a 60-minute experiment. **The run behind Table S6/S7 is a
different, longer experiment whose duration, temperature and sampling grid the article never
states.** Fitting units to a run of unknown length is not recovery; it is invention.

> **VERDICT ON 0.0028 / 0.0031: DO NOT INGEST AS NUMBERS.** The most probable unit is
> **ng·mL⁻¹·min⁻¹** and it should be recorded as *probable, unverified*, but the underlying
> experiment's duration and temperature are unstated and its regression table is absent.
> **What IS ingestable is the RATIO 0.0031/0.0028 = 1.107** — "there was not much of a
> difference between either" — and the **two-stage structure**, both of which are
> unit-independent. The same applies to 0.0043 / 0.0033 / 0.0032: ingest
> **Cys/control ≈ 3×, GSH/control ≈ 2.5× before 90 min; GSH/control ≈ 4× after**, and
> **Cys → control after 90 min**, as ratios.

## K3.4 — THE TWO-STAGE TEMPERATURE PROTOCOL, CONSOLIDATED

Assembled from §1.3, §5, §6 and §7 into the single statement the model needs:

**CHPBF (control, "conventional"):** pork *longissimus thoracis*, trimmed of subcutaneous fat,
mixed with pre-weighed ingredients, **115 °C for 60 min, one stage**.

**THPBF (treatment, the RSM optimum):** the same matrix **+ Cys 2.5 mg/mL + GSH 6.0 mg/mL**
(= the optimised **Cys:GSH = 3:7** ratio at the optimised total loading),
**115 °C for 50 min, then 95 °C for 10 min**.

The optimum was **chosen by a 4-factor × 3-level Box–Behnken RSM over 27 runs** with two
responses: R1 = DMTS concentration (minimise), R2 = MMFT/MFT-MFT ratio (maximise).
Model quality: **R² = 0.905 (R1)** and **0.873 (R2)**; model p = **0.0004** and **0.0019**;
lack of fit not significant. Verification error: **9.9 %** on R2 and **7.0 %** on R1.

**The four optimised factors and their chosen levels:**
| factor | optimum |
|---|---|
| Cys : GSH concentration ratio | **3 : 7** |
| initial (first-stage) reaction temperature | **115 °C** |
| **second-stage temperature** | **95 °C** (i.e. **−20 °C**) |
| second-stage reaction time | **10 min** |

**Measured effect of the two-stage protocol in the real pork matrix (Fig. 4a):**
**DMTS 3.60 ± 0.53 → 1.70 ± 0.15 µg/mL (×0.47)** and **MMFT/MFT-MFT ≈1.09 → ≈2.71 (×2.5)**.
Storage (Fig. 4c, semiquantitative): THPBF DMTS stays **2.4–2.6× below** CHPBF at days 56–112.

**★ Read as a temperature statement, this is the paper's strongest and it is a DESIGNED one:**
a four-factor optimisation, free to put the second stage anywhere in its design space, **chose to
drop the temperature by 20 °C for the last 10 min** in order to maximise the desirable thiol
ratio and minimise the polysulfide. Corroborated by the two single-factor statements in §7
("a lower initial reaction temperature could significantly reduce the concentrations of DMTS and
promote an increase in the MMFT/MFT-MFT concentration ratio"; "DMTS had a lower concentration
when … the second stage temperature was less than 105 ℃").

⚠️ **What it is NOT.** It is a statement about **DMTS** and about a **RATIO**. Neither is a
statement about **[MFT]** itself, and the paper contains no MFT-vs-temperature number anywhere.
A model that reproduces "lower final-stage T ⇒ higher MMFT/MFT-MFT and lower DMTS" has matched
this paper. A model that infers "therefore [MFT] rises when T falls" has gone beyond it.

## K3.5 — Cys vs GSH vs CYSTINE (GCys), CONSOLIDATED

**All three additives appear only in the simple VB1/Xyl model (Fig. 1, §2). The dose–response
(Fig. 2) and the kinetics (§4) drop GCys entirely and carry only Cys and GSH.** GCys is
therefore a **single-point** comparison, at equal weight (15.0 mg/mL) and one condition.

| additive | MW | **MFT ng/mL** | **MFT-MFT ng/mL** | **dimer/monomer (mass)** | MFT % of the 3-species sum |
|---|---:|---:|---:|---:|---:|
| none (VB1-Xyl control) | — | 0.07 | 0.04 | 0.571 | 49 % |
| **Cys** (cysteine) | 121.16 | **1.34** | 0.115 | **0.086** | **90 %** |
| **GSH** (glutathione) | 307.32 | 1.05 | 0.07 | **0.067** | **91 %** |
| **GCys** (cystine) | 240.30 | **2.01** | **1.09** | **0.542** | **64 %** |

**The three statements the paper makes, verbatim, and what each licenses:**

1. > "the **GCys** treatment group consistently had **higher concentrations of** … **DMDS and
   > DMTS** … with higher values in the case of GCys addition with 1 % (w/w) compared to the Cys
   > group. Moreover, the GCys group exhibited the **lowest MFT percentage** (Fig. 1), which was
   > attributed to **MFT being readily oxidized by free radicals, forming MFT-MFT, its dimer**."
2. > "an … in MFT and corresponding disulfides indicated that **Cys and GSH were more capable of
   > controlling MFT to form dimers than GCys**."
3. > "**Poor control ability may be explained that GCys is an oxidized form** of [Cys]" (p. 4);
   and the abstract: "**no discernible difference between Cys and GSH in dynamic profiles**."

**★ THE STRUCTURAL RESULT, in one line:** *the sulfur additive's REDOX STATE, not its
concentration, sets the MFT dimerisation branch.* GCys (the **disulfide** form) makes the **most**
MFT (2.01, 1.5× Cys) and simultaneously loses the **most** of it to the dimer (54.2 % vs 8.6 %) —
**a 6.3× swing in the branch fraction driven by nothing but which sulfur species is present, at
one temperature and one time.** Cys and GSH, both **free thiols**, agree with each other to
within 1.3× on the branch fraction (0.086 vs 0.067) while differing 1.28× in MFT level.

⇒ **A model that carries only a "cysteine concentration" term cannot express this.** The
falsifiable structural test is: **a reduced-sulfur additive and an oxidised-sulfur additive at
equal molar sulfur must give branch fractions differing by ≳5×, in the direction
oxidised ⇒ more dimer.** ⚠️ Note the equal-**weight** dosing means the molar sulfur is not
matched: Cys 15 mg/mL = 123.8 mM S; GSH = 48.8 mM S; GCys = 62.4 mM × 2 S = **124.8 mM S**.
**Cys and GCys are, by luck, nearly sulfur-matched (123.8 vs 124.8 mM S — within 0.8 %), while
GSH is 2.5× lower.** That makes the **Cys-vs-GCys contrast a genuinely sulfur-controlled
comparison** and strengthens the redox-state reading considerably. *(This molar conversion is
mine, not the paper's.)*

## K3.6 — ★★ VERDICT: DOES THIS PAPER QUANTIFY MFT FORMATION **AND** LOSS?

> **FORMATION: YES, but only as a dose response at one temperature and one time.**
> Fig. 1 gives four absolute MFT levels (0.07 / 1.34 / 1.05 / 2.01 ng/mL) across four sulfur
> additives; Fig. 2d gives twenty more across ten additive levels × two additives
> (0.2 → 23.6 ng/mL, a **118× span**). All at **115 °C / 60 min / pH 4.9 phosphate**.
> These are formation *endpoints*, not formation *rates*. **No rate of MFT formation is
> reported or derivable** — there is no time axis.
>
> **LOSS: YES — and this is the paper's real contribution. Two channels, both quantified as
> co-product mass, at processing temperature.**
> - **Oxidative dimerisation, MFT → MFT-MFT.** Fig. 2f: up to **5.75 ng/mL dimer against
>   23.6 ng/mL free MFT** (GSH, 50 mg/mL). Fig. 1: **54.2 %** of the MFT pool by mass sits in
>   the dimer in the GCys system.
> - **Radical coupling with MeSH, MFT + MeSH → MMFT.** Fig. 2e, ten levels × two additives; the
>   only channel the paper fits any rate law to (zero order, k = 0.0028 / 0.0031, **units
>   unresolved**, K3.3).
> - **A named third sink that is never measured**: "melanoidin/quinone ionic reactions"
>   (Introduction, p. 2, citing Tang 2013 and **Charles-Bernard 2005**). ⚠️ **That is the same
>   Charles-Bernard 2005 that supplies K1's thiol-consumption constants** — the two papers are
>   describing the same sink, one at 25 °C and one at 115 °C, and neither measures it at the
>   other's temperature.
>
> **THE AGGREGATE (my derivation, §3's flux arithmetic):** in the GSH arm at 115 °C / 60 min,
> **MFT already consumed into its two MEASURED products is 21–55 % of the MFT still free**;
> in the Cys arm, 3–19 %. These are **lower bounds** on total MFT consumption, because the third
> sink is unmeasured.
>
> **THE HONEST ONE-LINE ANSWER: the paper quantifies MFT LOSS as a MASS FRACTION at a single
> processing temperature, and MFT FORMATION as a DOSE RESPONSE at the same single temperature.
> It quantifies neither as a RATE with resolvable units, and it puts nothing on a temperature
> axis.**

## K3.7 — ★★ WHAT THIS LICENSES THAT K1's 30 °C CONSTANT CANNOT

**K1's declared gap, verbatim** (`k1_kinetic_parameters.md` §0 and §1d):

> "**No activation energy for thiol consumption exists anywhere in this basket.** Both thiol
> papers are single-temperature. Any T-extrapolation of the FFT/MFT sink is the repo's own
> assumption and must be labelled as such."
> … and, on the only Ea in reach: "**DO NOT INGEST.** Two labs, two matrices, two analytical
> bases, two points. The same treatment applied to FFT gives a **negative** Ea … **The
> literature reviewed supports NO activation energy for thiol consumption.**"

K1's usable thiol-sink constant is **k(FFT sink) ≈ 9 × 10⁻⁴ s⁻¹ at 30 °C**, in a ~10 g/L
coffee-solids matrix at pH 5–6 (Hofmann & Schieberle 2002 Fig. 6 [K1-fitted] 9.4 × 10⁻⁴,
Table 2 9.8 × 10⁻⁴, bounded below by Charles-Bernard's > 7.7 × 10⁻⁴ at 25 °C).

### What Zhang 2024 DOES license

| # | licence | why K1's constant cannot give it | strength |
|---:|---|---|---|
| **1** | **A thiol sink is OPERATING at 115 °C, with measured magnitude.** 21–55 % of free MFT already consumed into measured products in 60 min. | K1's constant is a **30 °C** number in a **coffee** matrix. Using it at 115 °C requires an Ea that **K1 explicitly refuses to supply.** Zhang gives a **direct observation at processing temperature** that needs no Ea at all. | **strong — this is the licence** |
| **2** | **The 115 °C sink is at least TWO channels, and they are DIFFERENT channels from K1's.** K1's is **covalent thioether addition to melanoidins/electrophiles** (Hofmann's < 1.5 % disulfide branch **rules dimerisation out at 30 °C**). Zhang's are **oxidative dimerisation** and **radical coupling to MeSH** — both essentially absent from K1's 30 °C system. | K1's 9 × 10⁻⁴ s⁻¹ is a **one-channel** constant, and Hofmann's own control says the disulfide route carries **< 1.5 %** of it. At 115 °C the dimer carries **up to 49 %** of the MFT pool. **The channel mix INVERTS between the two temperatures.** | **★ strongest structural finding** |
| **3** | **⇒ A SINGLE Arrhenius term on "thiol consumption" is therefore WRONG, and now demonstrably so.** Two different reactions dominate at the two temperatures; one Ea cannot describe the switch. | K1 declared "no Ea exists". Zhang upgrades that from *"we don't have one"* to **"one wouldn't be right if we did"**. | **★ this is the real gift** |
| **4** | **A branch-fraction target at 115 °C**: dimer/monomer 0.086 (Cys) to 0.542 (GCys); MMFT/MFT-MFT 1.97 → 0.04 across the additive dose range. | K1 has no branch information at any temperature above 30 °C. | strong, ratio-based |
| **5** | **A DIRECTION on temperature, from a designed experiment**: dropping the final stage 115 → 95 °C gives DMTS ×0.47 and MMFT/MFT-MFT ×2.5 in a real pork matrix. | K1 has **no** temperature direction of any kind for thiols. | moderate — about DMTS and a ratio, not [MFT] |
| **6** | **A two-regime rate structure**: the MMFT slope **breaks at 90 min** and, for Cys, falls back to the control rate. | K1's constants are all **single-exponential first order** with no regime change. | moderate — units unresolved, but the *break* is unit-free |
| **7** | **A redox-state axis on the additive** (K3.5): oxidised sulfur ⇒ 6.3× more dimerisation at near-matched molar sulfur. | K1's basket has no additive-redox axis at all. | moderate — one point, one condition |

### What Zhang 2024 DOES NOT license — stated as flatly as K1 stated its gap

> **IT DOES NOT SUPPLY AN ACTIVATION ENERGY FOR THIOL CONSUMPTION, AND IT CANNOT BE MADE TO.**
> Every kinetic number in it is at **115 °C**. The temperature sweep is **Fig. S2**, which is not
> in the retrieved file. **Pairing Zhang's 115 °C rate with K1's 30 °C rate to extract an Ea would
> repeat exactly the two-point, two-lab, two-matrix, two-method error that K1 named and refused
> (§1d) — and it would be worse, because the two rates are not even the same reaction (licence
> #2).** Record this as a **prohibited derivation**.
>
> It also does not supply: any absolute MFT concentration commensurable with a SIDA anchor
> (HS-SPME + external curve, §1.4); any MFT-vs-temperature number; any rate constant with
> resolved units; any dose–response for the *precursor* rather than the *additive*.

### ⇒ THE CONSOLIDATED POSITION FOR THE CONSUMPTION MODULE

**Three sources, three temperatures, three channels, and they do not compose into one Arrhenius
term.** This is now a three-point picture and it should be recorded as one:

| T | matrix | dominant measured MFT/thiol sink | magnitude | source |
|---:|---|---|---|---|
| **30 °C** | coffee melanoidins, pH 6, aqueous | **covalent thioether addition** (disulfide branch **< 1.5 %**) | k ≈ 9 × 10⁻⁴ s⁻¹, t½ ≈ 12 min | Hofmann & Schieberle 2002 + Charles-Bernard 2005, via **K1** |
| **50 °C** | ribose/cysteine process flavour, pH 5.0, 0.5 M phosphate | **acid-catalysed electrophilic oligomerisation at C-5** — *not* oxidative; air ≈ argon; mass balance fails to close as thiol + disulfide | **59 % of initial per day, ZERO ORDER** | van Seeventer 2001, via **Z2 finding 13** |
| **115 °C** | thiamine/xylose/methionine, pH 4.9 phosphate | **oxidative dimerisation** + **radical coupling to MeSH** | **21–55 % of free MFT consumed in 60 min** (lower bound) | **Zhang 2024, this paper** |
| **120 °C** | Amadori + cysteine, unbuffered water | **dimerisation, 6.5–9.6 % of MFT-equivalents, pH-INVARIANT over pH 6–8** | see `zhou2023_extraction.md` §2.1 | **Zhou 2023, Wave K3** |

**Four temperatures, four papers, and the dominant channel is different at each.** The 30 °C
result explicitly excludes the 115 °C channel (< 1.5 % disulfide); the 50 °C result explicitly
excludes oxidation (air ≈ argon); the 115 °C result is *predominantly* the oxidative channel that
the other two exclude.

> **RECOMMENDATION TO THE BUILD WAVES: implement thiol consumption as a SET OF NAMED CHANNELS
> with a declared temperature range of validity each, and an explicit `no_Ea_available` state —
> NOT as one lumped Arrhenius sink. Zhang 2024's contribution is that this is no longer a
> preference; it is what four papers at four temperatures measure.**

### K3.8 — ONE CROSS-PAPER LINK THAT ONLY EXISTS AFTER THIS WAVE

**Zhou 2023 (`zhou2023_extraction.md`, extracted this wave) independently measures Zhang's
dimerisation channel** — different lab (Jiangnan/Rutgers vs BTBU), different feedstock (purified
Amadori vs thiamine/xylose), different temperature (120 vs 115 °C), different pH (6–8 unbuffered
vs 4.9 buffered), different quantification (HS-SPME/1,2-DCB vs HS-SPME/IPDS).

| | Zhang 2024, 115 °C, pH 4.9 | Zhou 2023, 120 °C, pH 6–8 |
|---|---|---|
| dimer/monomer, **mass** | 0.067 (GSH) · 0.086 (Cys) · **0.542 (GCys, oxidised)** | 0.065 (pH 7) · 0.086 (pH 6) · 0.095 (pH 8) |
| dimer/monomer, **MFT-equivalents (molar ×2)** | **6.7 – 8.7 %** (the two free-thiol additives) | **6.5 – 9.6 %** (across pH 6–8) |
| **agreement** | — | **the two ranges OVERLAP; their midpoints are 7.7 % and 8.1 %, i.e. within 1.05×** |

**Two labs, two feedstocks, two pH regimes, ~5 °C apart, agree that ≈7–10 % of MFT sits as its
own disulfide dimer under free-thiol conditions.** That is the tightest cross-validation the
sulfur consumption lane has ever had, and it appears only when these two papers are read
together. **Zhou additionally shows the branch fraction is near-invariant in pH while [MFT]
swings 3.0× — so the dimerisation fraction is a much more stable quantity to target than the
level.**

⚠️ **And Zhou adds the qualification that changes what the sink MEANS**: its SI Table S2 gives
bis(2-methyl-3-furyl) disulfide an odour threshold of **0.00032 µg/L** against MFT's **0.005**,
i.e. **the dimer is 15.6× MORE potent than the monomer**, and carries a marginally *higher* OAV
than the MFT it came from. **Mass lost to the dimer is not aroma lost.** Any objective that
scores the dimerisation channel as a pure loss is wrong by roughly the ratio of the thresholds.
