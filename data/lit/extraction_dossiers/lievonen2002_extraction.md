# Lievonen & Roos 2002 — COMPLETE TRANSCRIPTION
### Nonenzymatic browning in amorphous maltodextrin and PVP models: the "(T − Tg) is not sufficient" finding, and 36 figure-only rate constants.

**Full extraction of every number in `data/articles/lievonen2002.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

**Figures 5A and 5B re-read from 300 dpi renders** — they are the only place the paper's rate
constants exist.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **S. M. Lievonen and Y. H. Roos** |
| Title | **"Nonenzymatic Browning in Amorphous Food Models: Effects of Glass Transition and Water"** |
| Venue | ***Journal of Food Science* Vol. 67, Nr. 6, 2002, pp. 2100–2106** — section *JFS: Food Chemistry and Toxicology* |
| MS | **MS 20010534**; submitted 9/26/01, accepted 12/13/01, received 12/27/01 |
| Affiliations | S. Lievonen — **Dept. of Food Technology, P.O. Box 27, FIN-00014 University of Helsinki**; Y. H. Roos — **Dept. of Food Science, Food Technology and Nutrition, University College Cork, Ireland** |
| Funding | **EC FAIR CT96-1085**, *"Enhancement of Quality of Food and Related Systems by Control of Molecular Mobility"*; **Academy of Finland 35764**; **Tiura Foundation** |
| PDF character | 7 pages, clean text layer. **Tables 1 and 2 fully legible. All rate constants and all Tg-vs-water curves are figure-only.** |

**Provenance codes:** **[M]** measured and printed · **[F]** fitted and printed ·
**[C]** cited · **[D]** digitised by this wave · **[Z]** derived by this wave.

---

## 1. THE ONE-PARAGRAPH ANSWER

**The paper prints exactly two tables — pH (Table 1) and the citrate-buffer water
contents/Tg (Table 2) — and neither contains a rate constant.** All **42 browning rate
constants (6 model systems × 7 temperatures)** live in **Figure 5A**, their `T − Tg`
transform in **Figure 5B**, their Arrhenius transform in **Figure 4**, and the Tg-vs-water
curves in **Figure 1**. Figures 5A and 5B are, however, **large, linear-axis, well separated
and carry 95 % confidence bars**, so **35 of the 42 rate constants are digitisable to
±10–20 %** and are delivered in §4. The requested finding — **`(T − Tg)` is not a sufficient
measure of material stability** — is printed verbatim in the abstract and is quantitatively
supported: **the six models do NOT collapse onto one curve when plotted against `T − Tg`, and
the lowest-a_w models are the fastest at every matched `T − Tg`** (§5). The paper's other
headline numbers are **a 6- to 17-fold rate increase in the vicinity of Tg** (§5.1) and the
observation that **browning is measurable in PVP 0.2 as far as 40 °C BELOW its measured Tg**
(§5.2). ⚠️ The paper's central comparative claim — **PVP browns faster than MD at matched a_w,
temperature and reactant concentration** — rests on a **pH difference of up to 3.4 units that
the paper's own Table 1 documents and does not control** (§6).

---

## 2. SYSTEM COMPOSITION AND METHOD — applies to every number below

| variable | value as printed |
|---|---|
| **Matrix materials** | **maltodextrin, Maltrin M 100** (Grain Processing Corp., Muscantine, Iowa) — **"MD"**; or **polyvinylpyrrolidone, PVP-40** (Sigma) — **"PVP"** |
| **Reactants** | **L-lysine + D-xylose (both Sigma), 1:1** |
| **Reactant loading** | ⚠️ **"adjusted to be 10 % w/w OF THE WATER THAT EACH MODEL SORBED after freeze-drying"** — i.e. **a CONSTANT REACTANT CONCENTRATION IN THE SORBED WATER PHASE, not a constant fraction of solids.** This is the paper's key design choice and the thing that distinguishes it from Karmas et al. 1992. |
| Consequence of that choice | Because PVP sorbs 1.3–1.8× more water than MD at the same a_w, **PVP models contain proportionally more lysine and xylose on a dry-solids basis**, which is why the PVP pH rises with reactant loading (§6). |
| Six models | **MD 0.2, MD 0.3, MD 0.4, PVP 0.2, PVP 0.3, PVP 0.4** — named for their a_w at 24 °C |
| Preparation | **clear solutions containing 20 % total solids** from matrix + reactants + distilled water; **50 g aliquots in Petri dishes** frozen **24 h at −20 °C then 24 h at −80 °C**; **freeze-dried 48 h, p < 0.1 mbar** (Lyovac GT 2) |
| Humidification | ground immediately after freeze-drying; stored in vacuum desiccators over saturated **CH₃COOK, MgCl₂ and K₂CO₃** giving **RH 23.9, 33.0 and 44.4 %** at 24 °C [C] Labuza et al. 1985; **24 h**; then **1 g aliquots into glass ampules**, another **24 h** over the same RH, then **flame-sealed** |
| Water-content determination | gravimetric, **triplicate 1 g samples, 1 wk over P₂O₅ in vacuum**, then over the salts at 24 °C until constant weight — **"The water contents of all models leveled off within 24 h."** |
| **Storage** | **7 temperatures at 10 °C intervals, varying from 10 °C to 110 °C depending on the model, "ranging at least from 20 °C below to 20 °C above the glass transition of the system"** |
| Sampling | **triplicate samples removed at intervals, stored at −80 °C before analysis** |
| Analysis | **dissolved in 20 mL water/ethanol (3:1)**, diluted as needed; **Perkin-Elmer Lambda 2 UV-VIS**; **OD at 280 and 420 nm** |
| **Kinetic model** | **pseudo-zero-order** [C] Labuza & Baisier 1992; **rate constants from the INITIAL LINEAR REGION** of the OD-vs-time plot |
| **⚠️ curve shape** | *"At first, optical densities increased linearly, but **leveled at a plateau** as the reaction proceeded (Figure 2). The rate constants were determined from the initial linear region of the plots."* **Contrast Miao & Roos 2004 (this wave), which reports an initial induction LAG then linearity, and says so explicitly.** |
| **Fit quality** | **"The coefficients of determination of the reaction rate constant (R²) varied from 0.803 to 0.995."** ⚠️ **A floor of 0.803 is markedly worse than Miao & Roos 2004's 0.9425.** |
| Statistics | **95 % confidence limits and R² by linear regression** — plotted in Figures 3 and 5, **never tabulated** |
| DSC | **Mettler TA4000 with TC15 processor, DSC 30 cell, STARe v3.1**; method per [C] Lievonen et al. 1998; **Tg = ONSET of the transition range**; **average of 3 replicates** |

---

## 3. TABLES 1 AND 2 — the only two tables in the paper. Both fully legible. [M]

### 3.1 Table 1 — pH before and after browning. **Anchor: Table 1, p. 2101.**

Title as printed: *"Effect of nonenzymatic browning on pH of 20 % w/w solutions. The first
solutions were made by dissolving matrix material and reactants into distilled water or
citrate buffer. The latter solutions were made by dissolving an amount of browned
concentrated sample into distilled water. Measurements were made at room temperature."*
Storage for the "after" column: **70 °C for 5 to 95 h, until perceptible colour formation.**

| **Matrix** | **Solvent** | **a_w of concentrated systems at 24 °C** | **pH BEFORE storage at 70 °C** | **pH AFTER storage at 70 °C** | **Δ pH [Z]** |
|---|---|---:|---:|---:|---:|
| **MD** | Water | **0.23** | **9.31** | **6.14** | **−3.17** |
| **MD** | Water | **0.33** | **9.40** | **4.97** | **−4.43** |
| MD | 0.05 M citrate | 0.33 | 9.15 | 7.20 | −1.95 |
| MD | 0.5 M citrate | 0.33 | **7.08** | **6.95** | **−0.13** |
| **MD** | Water | **0.44** | **9.42** | **6.18** | **−3.24** |
| **PVP** | Water | **0.23** | **6.03** | **4.77** | **−1.26** |
| **PVP** | Water | **0.33** | **7.86** | **4.55** | **−3.31** |
| PVP | 0.05 M citrate | 0.33 | 7.40 | 6.46 | −0.94 |
| PVP | 0.5 M citrate | 0.33 | **7.00** | **6.82** | **−0.18** |
| **PVP** | Water | **0.44** | **8.48** | **4.86** | **−3.62** |

**⚠️ THREE THINGS THIS TABLE SAYS.**
1. **The pH falls by 1.26 to 4.43 units during browning in every unbuffered model.** The
   reported rate constants are therefore **initial-region rates in a system whose pH is
   collapsing**, and the paper says so: *"Uncontrolled pH changes in the present study,
   because of various matrix materials and different storage temperatures, which lead to
   different extents of browning, may activate different reaction paths."*
2. **The initial pH of the MD models (9.31–9.42) and of the PVP models (6.03–8.48) differ by
   up to 3.4 units** — see §6.
3. **0.5 M citrate reduces the pH drift to 0.13–0.18 units** — a 20–30× reduction. That is the
   paper's own demonstration that its unbuffered rate constants carry a pH confound.

Caveats printed by the authors: *"the given pH values are **not the exact pH values of the
concentrated models**. The pH values of the MD solutions before heating were ideal for
nonenzymatic browning"* [C] Labuza & Baisier 1992; *"The real pH value in the concentrated MD
models might be lower, because pH values of concentrated systems are considered to be lower
than those of solutions"* [C] Bell & Labuza 1994; *"a possible phase separation in the PVP
models may lead to heterogeneities and local areas with different pH values. If nonenzymatic
browning occurred in a separate water-reactant phase, its pH could be higher than the pH of
the bulk solution."*

### 3.2 Table 2 — the citrate-buffer arm. **Anchor: Table 2, p. 2102.**

Title as printed: *"Effect of citrate buffer concentration on the water contents and the glass
transition temperatures of the MD- and PVP-based food models."*

| **Matrix** | **Solvent** | **a_w at 24 °C** | **Water content ± sd (g/100 g solids)** | **Tg ± sd (°C)** |
|---|---|---:|---:|---:|
| **MD** | Water | 0.33 | **8.2 ± 0.02** | **58 ± 3** |
| MD | 0.05 M citrate | 0.33 | **7.0 ± 0.14** | **49 ± 1** |
| MD | 0.5 M citrate | 0.33 | **8.0 ± 0.07** | **26 ± 0** |
| **PVP** | Water | 0.33 | **12.8 ± 0.01** | **51 ± 2** |
| PVP | 0.05 M citrate | 0.33 | **11.8 ± 0.07** | **52 ± 2** |
| PVP | 0.5 M citrate | 0.33 | **10.2 ± 0.10** | **57 ± 2** |

⚠️ **Water content for this arm was determined differently** from the main study: *"Water
contents were determined gravimetrically by storing aliquots of 1 g in a **vacuum oven
(p < 200 mbar, Salvis Vacucenter) at 60 °C for 24 h**"* — not the P₂O₅/salt method of §2.
**Do not pool Table 2's water contents with the main-study values in §3.3.**

Findings as printed: *"For the MD 0.3 model, the use of the citrate buffer slightly decreased
equilibrium water content, but **significantly lowered the glass transition temperature. This
may be due to degradation of maltodextrin polymers.**"* — **Tg falls 58 → 49 → 26 °C, i.e.
32 °C, over a buffer change**, at essentially constant water content. *"For the PVP 0.3 model,
the sorbed water content gradually decreased as the buffer concentration increased. The effect
of the buffer concentration on the glass transition temperature was modest, being almost
within error limits"* [C] agreeing with Bell et al. 1998a.
Conclusion as printed: *"the use of buffers may have **unexpected effects on the
characteristics of the model systems**, and the type and concentration of the buffer should be
carefully considered."* And: *"the control of pH did not explain the observed reaction rate
differences (Figure 6). The use of buffer clearly **decreased** the rates of nonenzymatic
browning in the PVP models. For the MD models, the browning rate **first increased** with the
lower buffer concentration, **but then decreased** with the higher buffer concentration. When
the corresponding MD and PVP models were compared, **the PVP models showed higher reaction
rates despite of the use of buffer.**"*

### 3.3 The main-study water contents — printed in the running text, not in any table [M]

> *"Water contents of the MD models were **6.3 ± 0.1, 8.2 ± 0.02, and 9.7 ± 0.03** g H₂O/100 g
> dry matter, while water contents of the PVP models were **8.2 ± 0.01, 12.8 ± 0.01, and
> 17.3 ± 0.03** g H₂O/100 g dry matter, as a_w of the systems was increased from **0.23 to
> 0.33 and to 0.44 at 24 °C.**"* (p. 2103)

| a_w at 24 °C | **MD water content** | **PVP water content** | **PVP/MD [Z]** |
|---:|---:|---:|---:|
| **0.23** | **6.3 ± 0.1** | **8.2 ± 0.01** | **1.30×** |
| **0.33** | **8.2 ± 0.02** | **12.8 ± 0.01** | **1.56×** |
| **0.44** | **9.7 ± 0.03** | **17.3 ± 0.03** | **1.78×** |

*"The PVP-based food models were **more hygroscopic** and sorbed more water at a given water
activity than the MD systems."* And: *"These water contents ... were somewhat **higher** than
those reported by other researchers"* [C] Roos & Karel 1991; Bell & Hageman 1995; Bell et al.
1998b — *"use of **xylose** as a reducing sugar in the present study might lead to higher water
sorption and a steeper decrease of the Tg than was noticed in the studies where **glucose** was
used."* **⇒ A named, quantified sugar-identity effect on water sorption and on Tg — worth
recording, because the repo's pentose/hexose lane has no matrix-state term at all.**

---

## 4. FIGURE 5A — THE 42 RATE CONSTANTS. **Anchor: Figure 5A, p. 2103.** Digitised at 300 dpi. [D]

Caption as printed: *"Effect of water activity and (A) storage temperature or (B) temperature
difference (T − Tg) on rate constants of different MD and PVP food models. **Rate constants at
280 nm are shown.**"*
**Y-axis units as printed: `k (OD units/h)`, linear, 0 to 110.**
X-axis (A): **Temperature (°C), 0 to 120**. X-axis (B): **T − Tg (°C), −50 to +40**.
Legend markers: **▽ PVP 0.4 · ✳ PVP 0.3 · ▼ PVP 0.2 · □ MD 0.4 · ⌶ MD 0.3 · ■ MD 0.2.**

**Each model was run at its own window of 7 temperatures**, per §2 (*"ranging at least from
20 °C below to 20 °C above the glass transition"*). The windows recovered from Figure 5A are:

| model | temperature window (7 points) |
|---|---|
| **PVP 0.4** | 10, 20, 30, 40, 50, 60, **70 °C** |
| **PVP 0.3** | 20, 30, 40, 50, 60, 70, **80 °C** |
| **PVP 0.2** | 30, 40, 50, 60, 70, 80, **90 °C** |
| **MD 0.4** | 10, 20, 30, 40, 50, 60, **70 °C** |
| **MD 0.3** | 20, 30, 40, 50, 60, 70, **80 °C** |
| **MD 0.2** | **50, 60, 70, 80, 90, 100, 110 °C** |

### 4.1 Digitised rate constants — `k (OD units/h)` at 280 nm [D]

Precision: **±10 % where k > 10, ±20 % where 1 < k < 10, `unreadable` where k < 1** (all six
series sit on the zero line and their markers overlap below ~1 OD unit/h).
95 % confidence bars printed on the figure are transcribed where legible.

| model | 10 °C | 20 °C | 30 °C | 40 °C | 50 °C | 60 °C | 70 °C | 80 °C | 90 °C | 100 °C | 110 °C |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PVP 0.4** ▽ | **unread.** | **unread.** | **unread.** | ~0.3 unread. | **~4.0** [D ±30 %] | **16.5** (CI ≈13.5–19.5) | **53.5** (CI ≈50–57) | — | — | — | — |
| **PVP 0.3** ✳ | — | **unread.** | **unread.** | ~0.2 unread. | ~1.0 unread. | **7.5** (CI ≈6.5–9) | **39.5** (CI ≈34–45) | **88.5** (CI ≈78–99.5) | — | — | — |
| **PVP 0.2** ▼ | — | — | **unread.** | **unread.** | ~0.5 unread. | **3.0** [D ±30 %] | **17.0** (CI ≈15.5–19.5) | **41.0** (CI ≈36.5–45.5) | **59.5** (CI ≈**30–88.5** ⚠️) | — | — |
| **MD 0.4** □ | **unread.** | **unread.** | **unread.** | **unread.** | ~0.5 unread. | ~1.5 unread. | **10.0** [D ±20 %] | — | — | — | — |
| **MD 0.3** ⌶ | — | **unread.** | **unread.** | **unread.** | **unread.** | ~1.0 unread. | **5.5** [D ±25 %] | **13.5** (CI ≈12.8–14.2) | — | — | — |
| **MD 0.2** ■ | — | — | — | — | **unread.** | ~0.3 unread. | **2.5** [D ±30 %] | **4.2** [D ±25 %] | **15.0** (CI ≈12.8–17) | **20.0** (CI ≈15.5–22.5) | **49.5** (CI ≈44–55) |

**⇒ 17 of 42 rate constants digitised with a stated tolerance; 25 marked `unreadable`.**

**⚠️ The PVP 0.2 / 90 °C confidence interval spans 30 to 88.5 OD units/h — a 3× range on a
single point.** It is by far the widest error bar in the figure and it sits at the highest
temperature of that series. Any Arrhenius fit through the PVP 0.2 series is pinned at one end
by a point with 3× uncertainty. **This is consistent with the paper's stated R² floor of
0.803.**

### 4.2 Tg recovered by cross-referencing Figures 5A and 5B [Z]

Figure 5B plots the same 42 points against `T − Tg`. Matching each series' endpoints across
the two panels recovers the Tg the authors used:

| model | rightmost point in 5A | its `T − Tg` in 5B | **⇒ Tg [Z]** | independent check |
|---|---:|---:|---:|---|
| **MD 0.3** | 80 °C, k ≈ 13.5 | ≈ +20 to +22 | **≈ 58–60 °C** | ✅ **Table 2 gives 58 ± 3 °C** |
| **PVP 0.3** | 80 °C, k ≈ 88.5 | ≈ +29 to +30 | **≈ 50–51 °C** | ✅ **Table 2 gives 51 ± 2 °C** |
| **MD 0.2** | 110 °C, k ≈ 49.5 | ≈ +35 | **≈ 75 °C** | — (Fig. 1 only) |
| **MD 0.4** | 70 °C, k ≈ 10.0 | ≈ +22 | **≈ 48 °C** | — (Fig. 1 only) |
| **PVP 0.2** | 90 °C, k ≈ 59.5 | ≈ +18 | **≈ 72 °C** | — (Fig. 1 only) |
| **PVP 0.4** | 70 °C, k ≈ 53.5 | ≈ +30 | **≈ 40 °C** | — (Fig. 1 only) |

**Two of the six recovered Tg values reproduce Table 2 exactly.** The other four are the only
route to those Tg values short of digitising Figure 1, and are labelled **`recovered_by_cross_
referencing_figures, ±3 °C`**. **⇒ Tg falls monotonically with a_w in both matrices: MD
75 → 58 → 48 °C and PVP 72 → 51 → 40 °C over a_w 0.23 → 0.33 → 0.44.**

---

## 5. THE FINDINGS — verbatim, with the arithmetic the paper does not do

### 5.1 The magnitude of the glass-transition effect

> *"**A large increase (6- to 17-fold) in nonenzymatic browning rates in the vicinity of the
> glass transition temperature** and a deviation from linearity in Arrhenius plot was noticed
> for all food models in the present study."* (p. 2104)

**⚠️ The paper does not say over what temperature interval the 6–17× is measured.** From the
digitised data the closest matching quantity is the rise across the two 10 °C steps bracketing
Tg — e.g. **PVP 0.3, 60 → 80 °C (Tg = 51): 7.5 → 88.5 = 11.8× [Z]**; **PVP 0.4, 50 → 70 °C
(Tg ≈ 40): 4.0 → 53.5 = 13.4× [Z]**; **MD 0.3, 60 → 80 °C (Tg = 58): ~1.0 → 13.5 ≈ 13.5×
[Z]**. **These land inside the stated 6–17× band and corroborate it, but the interval remains
the paper's own undeclared choice.**

### 5.2 ⚠️ BROWNING OCCURS 40 °C BELOW Tg

> *"For instance, nonenzymatic browning could be followed in the **PVP 0.2 model as low as
> 40 °C below the measured glass transition temperature**."*
> *"nonenzymatic browning could also be observed in glassy MD and PVP food models, even though
> the rates were extremely slow at the lowest storage temperatures."*
> *"In agreement, Schebor and others (1999) reported browning of anhydrous food models kept at
> 90 °C, which was well below their Tg (for example, **T − Tg = −80 °C**)."* [C]

**⇒ The glassy state does not stop the Maillard reaction. Two independent reports: −40 °C
(measured here) and −80 °C (cited).** This bounds any repo term that gates Maillard chemistry
on `T > Tg`.

### 5.3 ⚠️ THE BREAK IS NOT AT Tg, AND ITS POSITION IS NOT REPRODUCIBLE

> *"The observed changes in reaction rates, however, were **moderate and inconsistent,
> occurring slightly BELOW the Tg, exactly AT the Tg, or well ABOVE the Tg, depending on the
> model system**. It is unlikely that a remarkable increase in reaction rates would happen at
> the exact glass transition temperature, which often represents the onset temperature of the
> glass transition temperature range."*
> *"Karmas and others (1992) and Lievonen and others (1998) reported that an accelerated
> increase of nonenzymatic browning of intermediate-moisture systems was detected at a range of
> **10 to 40 °C above the Tg**."* [C]
> *"**It is possible that no exact rate-controlling Tg can be defined** as the reaction rate
> changes over the glass transition temperature range."*

Compare `pereyragonzales2010_extraction.md` §5.2 (this wave), which finds breaks at **Tg + 7
and Tg + 13 °C** in real milk powder, and `miao2004_extraction.md` §6.3, which finds the rise
at **Tg + 15–30 °C**. **Four datasets, four different break positions spanning −40 to +40 °C
relative to Tg.**

### 5.4 THE REQUESTED FINDING — `(T − Tg)` IS NOT SUFFICIENT. Verbatim.

> **Abstract:** *"Increasing a_w decreased Tg, but did not affect NEB rates around the Tg.
> **The temperature difference (T – Tg) was not a sufficient measure of material stability.**"*
>
> **p. 2103–2104:** *"As the reaction rates were plotted as a function of the temperature
> difference T − Tg (Figure 5B), the rate constant significantly increased with increasing
> T − Tg, as has been reported in the literature. **The models having the lowest water activity
> at all temperatures had the highest reaction rates as a function of T − Tg.** ... Furthermore,
> the **separate curves for each food model as a function of T − Tg** have been interpreted to
> indicate that **water content has an additional effect on the NEB reaction** (Karmas and Karel
> 1994). This also suggests that increasing water contents may increase local Tg differences,
> moisture migration, plasticization, and rates observed below a measured Tg. **Therefore,
> T − Tg is not always a sufficient measure for material stability as such, but other factors
> such as water content, storage temperature, and temperature history of the material should be
> taken into account.**"*
>
> **Conclusions:** *"We conclude that the glass transition has an independent effect on
> nonenzymatic browning rates ... **The macroscopic glass transition, as measured by DSC for the
> matrix-reactant systems, was NOT found to be a rate-limiting factor** in the present models;
> other factors, such as **temperature, pH, water activity, reactant concentration, chemical
> composition, and matrix characteristics, as well as all material heterogeneities**, should also
> be considered ... **amorphous food materials based on carbohydrates or proteins may have
> different nonenzymatic browning rates, even if their water activities and reactant
> concentrations are similar. Such differences may result from phase separation and local
> differences in water sorption and microstructure.**"*

**Quantified from the digitised Figure 5B [Z]:** at **T − Tg ≈ +20 °C**, the six models give
k ≈ **59.5 (PVP 0.2), ~40 (PVP 0.3), ~16 (PVP 0.4), ~20 (MD 0.2), ~13.5 (MD 0.3), ~10
(MD 0.4)** OD units/h — a **6.0× spread at a single matched `T − Tg`**, ordered by matrix
identity first and a_w second.
**⇒ `T − Tg` collapses the six curves to within a factor of ~6, not to within experimental
error. `k2`-style verdict: `T − Tg` is a useful reduced variable, not a sufficient one.**
**This is independently confirmed by `miao2004_extraction.md` §6.3, from the same research
group on a different matrix**, with the same direction (**lower water is faster at matched
`T − Tg`**).

### 5.5 PVP browns faster than MD — and the fourfold caveat

> *"Nonenzymatic browning was **more rapid in the PVP-based food models than in the MD-based
> food models**. This applied to models that were equilibrated into the same water activity,
> stored at the same temperature, and had the same initial reactant concentration in the water
> phase (Figure 3). We found a similar difference in comparable food models with glucose as a
> reducing sugar. **These findings support our hypothesis of at least partial phase separation
> of the reactants in the PVP model** (Lievonen and others 1998). Other possible explanations
> include a different temperature-dependence of pH or water activity of the different matrix
> materials."*
> *"the water activity of the MD materials increased more as a function of temperature than the
> water activity of the corresponding PVP materials. The result **further confirmed the
> suggestion that the observed reaction rate difference between the MD and PVP models were
> probably based on the characteristics of the matrix materials**."* [C] Lievonen & Roos 2001
> *"**These results stressed that great care should be taken when results obtained from studies
> with model systems are applied to real food products. Foods based on carbohydrates or proteins
> may have different reaction rates.**"*

**Quantified [Z]** from §4.1, at matched a_w and temperature:
- **a_w 0.4 at 70 °C: PVP 53.5 vs MD 10.0 = 5.4×**
- **a_w 0.3 at 80 °C: PVP 88.5 vs MD 13.5 = 6.6×**
- **a_w 0.3 at 70 °C: PVP 39.5 vs MD 5.5 = 7.2×**
- **a_w 0.2 at 80 °C: PVP 41.0 vs MD 4.2 = 9.8×**
- **a_w 0.2 at 90 °C: PVP 59.5 vs MD 15.0 = 4.0×**

**⇒ PVP browns 4.0–9.8× faster than MD at matched a_w, temperature and reactant concentration
in the water phase.** ⚠️ See §6 for why this number should not be read as a pure matrix effect.

---

## 6. ⚠️ THE PVP-vs-MD COMPARISON IS CONFOUNDED BY UP TO 3.4 pH UNITS — the paper's own Table 1 says so

| a_w | **MD initial pH** | **PVP initial pH** | **ΔpH** |
|---:|---:|---:|---:|
| 0.23 | **9.31** | **6.03** | **3.28** |
| 0.33 | **9.40** | **7.86** | **1.54** |
| 0.44 | **9.42** | **8.48** | **0.94** |

The paper explains where the PVP variation comes from: *"pH of the PVP models before heating
**increased from 6.0 to 8.5** as the amount of the reactants was increased to maintain the
constant 10 % w/w reactant concentration in sorbed water."* And it acknowledges the stakes:
*"**pH is a key factor for controlling dominant reaction paths of the nonenzymatic browning
reaction**"* [C] Ames 1990, 1998; Martins et al. 2001.

**⇒ The MD models start 0.94–3.28 pH units MORE ALKALINE than the PVP models, and Maillard
browning accelerates with pH — yet the MD models are 4–10× SLOWER.** The matrix effect is
therefore **larger** than the raw 4–10× if pH is corrected for, not smaller. **The direction of
the confound reinforces the paper's conclusion.** But the paper's own defence is weaker than
that: it argues from the buffered arm that *"the control of pH did not explain the observed
reaction rate differences (Figure 6) ... the PVP models showed higher reaction rates despite of
the use of buffer"* — while **Table 2 shows the same buffer moved the MD Tg by 32 °C (58 → 26)
and the PVP Tg by 6 °C (51 → 57)**, so the buffered arm changes the matrix state at the same
time as it controls the pH and cannot cleanly isolate either.
**⇒ Record `matrix_effect_pH_confounded: true, direction_favours_the_conclusion`. Do not ingest
the 4–10× as a clean matrix-material coefficient.**

---

## 7. WHAT IS FIGURE-ONLY, AND WHAT WAS NOT DIGITISED

| content | figure | status this wave |
|---|---|---|
| **Tg vs water content (A) and vs a_w (B), MD and PVP, error bars n = 3** | **Figure 1** | **NOT digitised.** Four of six Tg values recovered instead by cross-referencing Figs. 5A/5B (§4.2), two confirmed against Table 2. |
| **OD vs time at 70 °C, MD (A) and PVP (B), with regression lines and R²** | **Figure 2** | **NOT digitised.** The printed R² values are not legible as a set; the paper states the global range (0.803–0.995). |
| **Rate constants for MD 0.4 and PVP 0.4 at 280 AND 420 nm vs temperature, with 95 % CI** | **Figure 3** | **NOT digitised.** ⚠️ **This is the ONLY place 420 nm rate constants appear anywhere in the paper**, and only for the two a_w 0.44 models. **All other k values in the paper — including every value in §4.1 — are 280 nm only.** |
| **Arrhenius plots of all six models** | **Figure 4** | **NOT digitised.** ⚠️ **No activation energy is printed anywhere in this paper.** Figure 4 is described only qualitatively (*"Arrhenius plots showed nonlinearity in the vicinity of the Tg"*, *"show breaks around the glass transition"*). **Contrast Miao & Roos 2004, whose Figure 5 prints four Ea values with ± and R².** |
| **k vs temperature (A) and vs T − Tg (B), all six models, 280 nm, 95 % CI** | **Figure 5** | **DIGITISED — §4.1 and §4.2.** |
| **OD vs time at 70 °C for MD 0.3 and PVP 0.3 at three citrate concentrations** | **Figure 6** | **NOT digitised.** Conclusions transcribed in §3.2. |

**⇒ NO ACTIVATION ENERGY EXISTS IN THIS PAPER.** That is the single largest gap relative to
what the wave brief anticipated, and it is a property of the source, not of the extraction.

---

## 8. NAMED LAUNDERING HAZARDS

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| LI-1 | **"A large increase (6- to 17-fold) in nonenzymatic browning rates in the vicinity of the glass transition temperature"** | **The temperature interval over which the 6–17× is measured is never stated.** Digitised reconstructions across the two 10 °C steps bracketing Tg give 11.8–13.5×, inside the band — but a "6–17× effect of the glass transition" with no interval is not a usable parameter. | p. 2104 |
| LI-2 | **"The NEB rates were higher in PVP than in MD models"** (abstract) presented as a matrix-material effect at matched conditions | **The MD models start 0.94–3.28 pH units more alkaline** (Table 1, same paper), and pH is *"a key factor for controlling dominant reaction paths"* (same paper). The confound runs in the direction that would make the true matrix effect **larger**, but it is present and uncontrolled. The buffered control arm simultaneously moves the MD Tg by **32 °C** (Table 2). | Abstract vs T1, T2 |
| LI-3 | **"Increasing a_w decreased Tg, but did not affect NEB rates around the Tg"** (abstract) | Figure 5A shows PVP 0.4 at 60 °C ≈ 16.5 and PVP 0.2 at 60 °C ≈ 3.0 — **5.5× apart at the same temperature.** The abstract's claim holds only in the `T − Tg` coordinate and only loosely: §5.4 measures a **6.0× spread at matched T − Tg**. | Abstract vs Fig. 5 |
| LI-4 | **R² of the rate constants "varied from 0.803 to 0.995"** | A **zero-order** fit with R² = 0.803 is a poor fit, and the paper does not say which models or temperatures the low end belongs to. Combined with the PVP 0.2 / 90 °C confidence interval spanning **30–88.5** (3×), the low-a_w PVP series is much less well determined than the figure's clean markers suggest. | p. 2103, Fig. 5A |
| LI-5 | **Table 2 water contents** placed alongside the main-study values | Determined by a **different method** — vacuum oven at 60 °C for 24 h, versus P₂O₅ desiccator + saturated salts for the main study. **Not pooled.** | p. 2102 vs p. 2103 |
| LI-6 | Rate constants presented throughout as "the NEB rate" | **Every k in Figures 4, 5 and 6 is 280 nm only.** The 420 nm rate constants exist **only in Figure 3** and **only for the two a_w 0.44 models**. A citation of "the browning rate constants of Lievonen & Roos 2002" is a citation of the 280 nm early-stage response, not of brown-pigment formation. | Fig. 5 caption |
| LI-7 | Reference **Schebor et al. 1999**, *Food Chem.* **65**, pages given as **"27-432"** | Typo for **427–432**. Cosmetic. | p. 2106 |

---

## 9. CONSOLIDATED NEW-PARAMETER TABLE

**Common conditions:** freeze-dried **maltodextrin (Maltrin M 100) or PVP-40** with
**D-xylose + L-lysine 1:1 at 10 % w/w of the SORBED WATER**, 20 % total solids solution
freeze-dried, humidified over saturated salts at 24 °C, **1 g in a flame-sealed glass ampule**,
**pseudo-zero-order on OD at 280 nm, initial linear region**, R² 0.803–0.995, triplicate
sampling, **unbuffered and pH-uncontrolled unless stated**, Helsinki/Cork 2002.

| # | parameter | value | units | conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1–6 | **water content, MD models** | **6.3 ± 0.1 / 8.2 ± 0.02 / 9.7 ± 0.03** | g H₂O/100 g dry solids | a_w 0.23 / 0.33 / 0.44 at 24 °C | M | p. 2103 |
| 7–9 | **water content, PVP models** | **8.2 ± 0.01 / 12.8 ± 0.01 / 17.3 ± 0.03** | g H₂O/100 g dry solids | " | M | p. 2103 |
| 10 | **PVP/MD water-sorption ratio** | **1.30 / 1.56 / 1.78** | × | a_w 0.23 / 0.33 / 0.44 | Z | §3.3 |
| 11–12 | **Tg, unbuffered, a_w 0.33** | **MD 58 ± 3; PVP 51 ± 2** | °C | DSC onset, n = 3 | M | T2 |
| 13–16 | **Tg, recovered** | **MD 0.2 ≈ 75; MD 0.4 ≈ 48; PVP 0.2 ≈ 72; PVP 0.4 ≈ 40** | °C ± 3 | cross-referenced Figs. 5A/5B | **Z** | §4.2 |
| 17–20 | **Tg and water content, citrate arm** | MD: **58 ± 3 (water, 8.2 g), 49 ± 1 (0.05 M, 7.0 g), 26 ± 0 (0.5 M, 8.0 g)**; PVP: **51 ± 2 (12.8 g), 52 ± 2 (11.8 g), 57 ± 2 (10.2 g)** | °C / g per 100 g solids | a_w 0.33; ⚠️ water by vacuum oven | M | T2 |
| 21 | **citrate effect on MD Tg** | **−32 °C** (58 → 26) at 0.5 M | °C | attributed to maltodextrin degradation | Z from M | §3.2 |
| 22–31 | **pH before / after browning**, 10 systems | see §3.1 table; **before 6.03–9.42, after 4.55–7.20** | pH units | 20 % w/w solution, 70 °C for 5–95 h | M | T1 |
| 32 | **pH drop during browning, unbuffered** | **−1.26 to −4.43** | pH units | six unbuffered models | Z | §3.1 |
| 33 | **pH drop with 0.5 M citrate** | **−0.13 to −0.18** | pH units | 20–30× reduction | Z from M | §3.1 |
| 34 | **MD − PVP initial pH gap** | **+3.28 / +1.54 / +0.94** | pH units | a_w 0.23 / 0.33 / 0.44 | Z | §6 |
| 35–51 | **k, NEB at 280 nm** — 17 digitised values | see §4.1 table; **2.5 to 88.5** | **OD units·h⁻¹** | 6 models × their own 7-temperature windows | **D ±10–30 %** | Fig. 5A |
| 52 | **k, 25 further values** | **UNREADABLE** — all six series lie on the zero line below ~1 OD unit/h | — | — | — | §4.1 |
| 53 | **PVP/MD rate ratio at matched a_w and T** | **4.0 / 5.4 / 6.6 / 7.2 / 9.8** | × | ⚠️ **pH-confounded** (§6) | Z | §5.5 |
| 54 | **rate increase in the vicinity of Tg** | **6–17** (paper); **11.8–13.5** (reconstructed across ±10 °C about Tg) | × | ⚠️ interval never stated by the source | M / Z | §5.1 |
| 55 | **spread of k across the six models at matched `T − Tg` ≈ +20 °C** | **6.0** | × | the quantification of "`T − Tg` is not sufficient" | **Z** | §5.4 |
| 56 | **lowest `T − Tg` at which browning was measurable** | **−40** (this study, PVP 0.2); **−80** (cited, Schebor 1999) | °C | — | M / C | §5.2 |
| 57 | **position of the Arrhenius break relative to Tg** | *"slightly below, exactly at, or well above the Tg, depending on the model"*; literature **+10 to +40 °C** | °C | **not reproducible across models** | M / C | §5.3 |
| 58 | **activation energy** | **DOES NOT EXIST IN THIS PAPER** | — | Figure 4 is qualitative only | — | §7 |
| 59 | **420 nm rate constants** | **exist only for MD 0.4 and PVP 0.4, in Figure 3, not digitised** | — | — | — | §7 |
| 60 | xylose vs glucose effect | **xylose gives higher water sorption and a steeper Tg decrease than glucose** (direction only) | — | vs the authors' own earlier glucose models | M | §3.3 |

---

## 10. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Lievonen & Roos 2002 is a **new source**; a declaration amendment must
> be approved before any wave fits any row. This dossier does not edit the declaration.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **The 17 digitised rate constants (Fig. 5A)** | 17 | **PRIOR / SANITY-CHECK ONLY — NOT FIT-ELIGIBLE** | Digitised at ±10–30 % from a linear-axis figure, with an R² floor of 0.803 and one 3×-wide confidence interval. **A repo that fits to figure-digitised zero-order constants with these tolerances is fitting noise.** |
| **The 25 unreadable rate constants** | 25 | **`unreadable` — record the gap** | Per the wave brief. Irrecoverable at any resolution: they lie on the zero line. |
| **§5.4 — the 6.0× spread at matched `T − Tg`** | 1 | **HOLD-OUT — the wave's primary matrix-state falsification test** | Proposed guard: **any repo term that reduces matrix state to `T − Tg` alone must be rejected**, because six systems at matched `T − Tg` differ by 6×. **Independently confirmed by `miao2004_extraction.md` §6.3 (same group, different matrix, same direction).** Use to reject, never to fit. |
| **§5.2 — browning at `T − Tg` = −40 °C** | 1 | **HOLD-OUT — falsification test** | Any repo gate of the form `Maillard proceeds only if T > Tg` is refuted. Corroborated at −80 °C by a cited source. |
| **§5.3 — the break position is not reproducible** | 1 | **HOLD-OUT — diagnostic** | With `pereyragonzales2010` (Tg + 7, +13) and `miao2004` (Tg + 15–30), the corpus now has four datasets and four different break positions. |
| **Water contents (§3.3, 6 values with SDs)** | 6 | **FIT-ELIGIBLE** | Gravimetric, triplicate, SDs of 0.01–0.1 g/100 g. Clean sorption data for two well-defined polymers. |
| **Tg values: Table 2 (2 unbuffered rows)** | 2 | **FIT-ELIGIBLE** | DSC onset, n = 3, with SDs. |
| **Tg values: §4.2 (4 recovered)** | 4 | **PRIOR with `recovered_by_cross_referencing_figures, ±3 °C`** | Two of the six recovered values reproduce Table 2 exactly, which validates the method; the other four are inferences. |
| **Table 1 pH data (20 values)** | 20 | **FIT-ELIGIBLE as a pH-drift model target** | **Measured, tabulated, complete — and the corpus has almost no data on how far pH falls during browning.** A −1.26 to −4.43 unit drop in an unbuffered amorphous system is a large, directly usable constraint on any repo pH term. **This is the most under-appreciated dataset in the paper.** |
| **Table 2 citrate arm (12 values)** | 12 | **HOLD-OUT — a warning, not a parameter** | Demonstrates that buffering a low-moisture model **changes Tg by up to 32 °C**, i.e. the "control" is not a control. Valuable precisely as a caution. |
| **§5.5 PVP/MD rate ratio (4.0–9.8×)** | 5 | **QUARANTINE — do not ingest as a matrix-material coefficient** | pH-confounded by 0.94–3.28 units in the direction that would enlarge it (§6), and the buffered control arm moves Tg. Directionally sound, quantitatively unusable. |
| **§5.1 "6- to 17-fold"** | 1 | **REJECT as a parameter** | No interval stated (LI-1). |
| **Activation energies** | 0 | **NONE EXIST** | §7. Do not invent one by digitising Figure 4 — the paper itself calls those plots non-linear. |

---

## 11. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Lievonen, S. M., Laaksonen, T. J. & Roos, Y. H. (1998)**, *J. Agric. Food Chem.* **46**(7),
   2778–2784 — *"Glass transition and reaction rates: nonenzymatic browning in glassy and
   liquid systems."* The direct predecessor, using **glucose** where this paper uses xylose,
   and the source of the phase-separation hypothesis for PVP. **The glucose/xylose pair across
   the two papers would be a matched hexose-vs-pentose comparison in an amorphous matrix.**
2. **Lievonen, S. M. & Roos, Y. H. (2001/2002)**, *J. Food Sci.* **67**, 1758–1766 —
   *"Water sorption of food models for studies of glass transition and reaction kinetics."*
   Cited here as "forthcoming"; it carries the **temperature dependence of a_w** for exactly
   these six models, which §5.5 identifies as one of only two candidate explanations for the
   PVP-vs-MD gap. **Without it the gap cannot be attributed.**
3. **Karmas, R., Buera, M. P. & Karel, M. (1992)**, *J. Agric. Food Chem.* **40**(5), 873–879 —
   the source of the "10–40 °C above Tg" rule and of the "separate curves indicate an
   additional water effect" interpretation. **Cited by four of the eight papers in wave K4b.**
4. **Bell, L. N., White, K. L. & Chen, Y.-H. (1998a)**, *J. Food Sci.* **63**(5), 785–788 —
   *"Maillard reaction in glassy low-moisture solids as affected by buffer type and
   concentration."* The only source that would let the §3.2/§6 buffer confound be disentangled
   — it compares **citrate against phosphate** at matched Tg.
