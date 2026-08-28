# Martins & Van Boekel 2003 (10.1016/S0008-6215(03)00174-5) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/martins2003.pdf` (14 pp.). Read method: both — `pdftotext -layout`
text layer for the body and all tables, plus page rasters read with the Read tool for Scheme 3
(PDF p. 7 / journal p. 1671) and Fig. 7 (PDF p. 11 / journal p. 1675), which the text layer
cannot render.

---

## 0. PAPER IDENTITY — **DOES NOT MATCH THE ORCHESTRATOR'S EXPECTED IDENTITY**

> **IDENTITY CORRECTION, STATED IN BOLD AT THE TOP AS REQUIRED.**
> **`data/articles/martins2003.pdf` IS PART II (Kinetic analysis), Carbohydr. Res. 338 (2003)
> 1665–1678.** The orchestrator's wave brief expected Part II to be `martins2003c.pdf`.
> **That is wrong.** `martins2003c.pdf` is a completely different paper (melanoidin extinction
> coefficients, Food Chemistry 83 (2003) 135–142). **This file — `martins2003.pdf` — is where the
> rate-constant table lives.**
>
> A second, independent naming hazard: **the authors' own lettered self-citations are permuted
> relative to our repo keys.** In Martins & Van Boekel (2005) Food Chem. 92, the reference list
> reads:
> - "Martins & Van Boekel (**2003a**)" = Part II = our `martins2003.pdf`
> - "Martins & Van Boekel (**2003b**)" = extinction coefficient = our `martins2003c.pdf`
> - "Martins, Marcelis & Van Boekel (**2003c**)" = Part I = our `martins2003b.pdf`
>
> In Martins & Van Boekel (2005) Food Chem. 90 the letters are permuted *again* ("2003a" = Part I,
> "2003b" = Part II, "2003c" = extinction coefficient). **Never resolve a Martins 2003 citation by
> its letter suffix. Resolve it by journal + volume + page.**

| field | value |
|---|---|
| authors | Sara I.F.S. Martins, Martinus A.J.S. Van Boekel |
| title | "Kinetic modelling of Amadori N-(1-deoxy-D-fructos-1-yl)-glycine degradation pathways. Part II — Kinetic analysis" |
| venue | Carbohydrate Research |
| volume/pages/year | 338 (2003) 1665–1678 |
| DOI | 10.1016/S0008-6215(03)00174-5 |
| received / accepted | Received 11 December 2002; accepted 11 April 2003 |
| affiliation | Department of Agrotechnology and Food Sciences, Product Design and Quality Management Group, Wageningen University, P.O. Box 8129, 6700 EV Wageningen, The Netherlands |
| PDF character | born-digital (Elsevier / Acrobat Distiller 4.05 for Windows); clean embedded text layer; the minus sign and the "±" glyph both degrade to `/` or `9/` in `pdftotext` output — every number below was re-verified against the page raster |
| corresponding author | tiny.vanboekel@wur.nl |
| funding | Fundação para a Ciência e Tecnologia (Portugal) |

### PRE-EXISTING-COVERAGE FLAG — DOUBLE-COUNTING HAZARD

`docs/reference/FIT_HOLDOUT_DECLARATION.md` (read-only; Amendment 1, 2026-08-28, lines ~209–215)
declares roles for "Martins 2003 Wageningen thesis (edepot.wur.nl/121418)" Tables 4.2.3, 4.1.1,
3.3.1 and Ch. 6. **These four journal papers are the thesis chapters**, so those declarations
already cover this material in some form.

**My mapping proposal for this paper (proposal, NOT a verified fact — I have not opened the
thesis PDF):**

| declared thesis item | my proposed mapping | confidence | verdict |
|---|---|---|---|
| Table **4.2.3** "reverse rates (Amadori → parent sugar), 100/120 °C × pH 5.5/6.8" | **Candidate match on conditions only.** This paper's Table 3 is the *only* Martins table on the 100/120 °C × pH 5.5/6.8 grid, so it is the most likely referent. **BUT this paper's Table 3 contains NO reverse rate constants of any kind.** All 16 estimated constants are irreversible forward steps of Model 2 (Scheme 3). | medium on the referent, **high on the negative finding** | **The declared "reverse rates (Amadori → parent sugar)" DO NOT EXIST in this paper.** See §5 and the DEFECT box below. |
| Table **5.1** "model discrimination (SB ↔ DFG reversible, ΔAICc 276 over irreversible)" | **Maps to a DIFFERENT paper, not this one:** Martins & Van Boekel (2005) *Food Chem.* **90**, 257–269, `data/articles/martins2005.pdf`, its Table 1 (ΔAICc: A = 0, B = 276.00, C = 287.46). This paper's own Tables 1 and 2 are M1-vs-M2 discrimination with ΔAICc values 0/10.96/11.50/123.62/127.86/173.32/423.88 — **276 appears nowhere here.** | high | not a double-count against this paper |
| Table 4.1.1 glycine-release yields | maps to Part I (`martins2003b.pdf`) Table 1, not here | high | n/a |
| Table 3.3.1 melanoidin ε | maps to `martins2003c.pdf` Table 1, not here | high | n/a |
| Ch. 6 pH ladder | maps to `martins2005b.pdf` Tables 2 and 3, not here | high | n/a |

**Are the numbers the same measurements re-published?** Within this paper's own family, yes and
importantly so: **Figs. 4–7 of this paper plot the SAME experimental dataset as Part I
(`martins2003b.pdf`) Figs. 1–8**, with model curves overlaid. Part I and Part II are one
experiment reported twice. Any wave that fits Part I's time-courses and also fits Part II's
Figs. 4–7 would be double-counting a single dataset. Part II's Table 3, by contrast, is *derived
from* that dataset (regression output), not an independent measurement — so **Table 3 and the
Part I time-courses are not independent evidence either.**

---

## 1. ONE-PARAGRAPH VERDICT

This paper gives the model a complete, printed, 16-parameter multiresponse rate-constant table
(Table 3) for thermal degradation of the isolated Amadori compound DFG, on a 2 × 2 grid of
100/120 °C × pH 5.5/6.8, each constant with a 95% highest-posterior-density (HPD) interval, plus
a fully specified reaction scheme (Scheme 3) and the differential equations for its predecessor
(Scheme 2). It also gives model-discrimination statistics (SS, n, AICc, ΔAICc, PPB) for two
competing topologies and two fructose-origin hypotheses. **It reports NO activation energies —
there is no Arrhenius analysis anywhere in this paper**; the two temperatures are fitted as
separate independent systems, so any Ea must be derived by the reader ([Z]) from two points only.
**It reports NO reverse rate constants** — every one of the 16 steps in Scheme 3 is a
unidirectional arrow, and the one "reverse-looking" constant (k12, mannose → glucose) is
explicitly the surviving half of a pair whose opposite direction (Glc → Man) the authors
*deliberately deleted* on literature grounds. **The starting material is DFG, not glucose +
glycine**, so this paper contains no Maillard initiation step and no glucose/glycine
second-order rate constant. Concentrations are all mmol/l (mM); every rate constant is
first-order (min⁻¹) except none — there is no second-order constant except k14 (Gly + Cn → Mel),
whose printed unit is *also* min⁻¹, which is a **unit defect** (see DEFECT box).

---

## 2. SYSTEM DEFINITION — verbatim

Part II's own §2 is only "2.1. Kinetic modelling". **The experimental system is defined in Part I**
(`martins2003b.pdf` §2), whose conditions Part II inherits verbatim. Reproducing the operative
facts, with the Part I anchor given:

| variable | value as printed | anchor |
|---|---|---|
| reactant | N-(1-deoxy-D-fructos-1-yl)-glycine (DFG), compound **1**, synthesised in-house (19% yield on glycine) plus material donated by Nestlé Research Center Lausanne | Part I §2.2, p. 1653 |
| initial concentration | "Compound 1 (0.237 g, 10 mmol) was dissolved in phosphate buffer (100 mL …)" → **10 mmol/l** | Part I §2.4, p. 1654 |
| buffer | "phosphate buffer (100 mL, **0.1 M** K₂HPO₄/KH₂PO₄, pH 6.8 and 5.5)" | Part I §2.4, p. 1654 |
| pH | initial 5.5 and 6.8, **set at room temperature, NOT hot-pH corrected, NOT held constant.** Measured pH drop at initial pH 6.8: 0.19 units after 180 min at 100 °C; 0.21 units after 120 min at 120 °C. At initial pH 5.5 the drop was **0.8 and 0.9** respectively | Part I §3.1.6 + Fig. 7, p. 1658–1659 |
| temperature | 100 °C (max 4 h) and 120 °C (max 3 h), oil bath | Part I §2.4 |
| vessel | screw-capped glass tubes (Schott, 16 × 160 mm) | Part I §2.4 |
| headspace/atmosphere | **not stated** (sealed screw-cap; no inerting reported) | — |
| volume | 100 mL prepared; sub-samples withdrawn | Part I §2.4 |
| agitation | **not stated** | — |
| quench | "samples were taken and immediately cooled in ice water, prior to analyses" | Part I §2.4 |
| replication | "Each reaction mixture was prepared, heated and analysed in **at least duplicate**" | Part I §2.4 |
| error definition (data) | Part I figures: error bars are drawn but their definition is not stated in Part I; Part II says only "the scatter in replicates was not very high which allowed a goodness-of-fit test" (§3.4, p. 1670) | — |
| error definition (parameters) | "**95% highest posterior density (HPD)**" — Table 3 caption, p. 1676. This is the Bayesian analogue of a confidence interval, not a standard error. |

Analytical methods, all inherited from Part I §2.5–2.10:
- **DFG and glycine:** Biochrom 20 Amino Acid Analyser (Pharmacia), lithium-citrate step gradient,
  post-column ninhydrin, 570 nm. Internal standard **glutamic acid, 25 nmol/µL (printed "25 mmol/mL")**.
  Quantification by external standards. t_R: DFG 9.5 min, glycine 38.8 min.
- **3-DG, 1-DG, methylglyoxal:** derivatised with 1,2-diaminobenzene (1 mol/L methanolic, overnight
  at 25 °C) to the quinoxalines; RP-HPLC, acetonitrile/ammonium acetate (pH 3.5, 20 mmol/L)
  gradient 10→30% MeCN in 50 min, 0.4 mL/min, detection 320 nm. External standards
  2-methyl-3-(1,2,3-trihydroxypropyl)quinoxaline (5′) and 2-methylquinoxaline (10′);
  "Both quinoxalines have the same extinction coefficient."
- **Glc, Man, formic acid, acetic acid:** HPLC ION-300 ion-exchange, 2.5 mmol/L H₂SO₄, 0.4 mL/min,
  85 °C; RI for sugars (Glc 15.1 min, Man 16.4 min), 210 nm for acids (FA 21.6 min, AcOH 23.7 min).
  **Man and Fru co-elute on this column.**
- **Glc / Man / Fru separated:** Dionex CarboPac PA100, NaOH/NaOAc gradient, 1 mL/min, ED-40
  electrochemical detection (Glc 12.95, Man 13.9, Fru 16.59 min).
- **Melanoidins:** A₄₇₀ on a Pharmacia Biotech spectrophotometer, converted to concentration by
  Lambert–Beer with **ε = 0.64 ± 0.03 l mmol⁻¹ cm⁻¹ at 470 nm**, cited to the extinction-coefficient
  paper (= our `martins2003c.pdf`). Concentration is therefore expressed as **moles of sugar
  incorporated into the brown polymer**, not moles of melanoidin.

### Fitting / statistics methods (Part II §2.1, p. 1666)

> "Based on the established reaction network a kinetic model was proposed and translated into a
> mathematical model by setting-up differential equations for each reaction step. The software
> package **Athena Visual Workbench** was used for numerical integration as well as for parameter
> estimation. The model parameters, the rate constants, were estimated by non-linear regression
> using the **determinant criterion**, that is to minimize the determinant of the matrix of
> cross-products of the various responses, so called dispersion matrix."

Integration subroutine: **DDAPLUS** (§3.3, p. 1668). Discrimination tools: multivariate
goodness-of-fit test built into Athena; **Posterior Probability (PPB)**; **Akaike Criterion
(AIC)**, corrected form used.

Printed equations, verbatim (p. 1671):

- AIC = n ln(σ²) + 2(p + 1)   … (15)
- σ² = SS / n   … (16)
- AIC_c = n ln(σ²) + 2(p + 1) · ( n / (n − p) )
- ΔAIC = AIC − AIC_min   … (17)

Rule of thumb, verbatim (p. 1672):
> "models with ΔAIC ≤ 2–3 are worthwhile to consider, values of ΔAIC between 4 and 7 indicate that
> models are less supported, and values higher than 10 indicate that models may be discarded."

---

## 3. SIMPLE-KINETICS RESULT (§3.1, p. 1666) — for the record

**Anchor: §3.1 body text, p. 1666 (PDF p. 2). No table.**

- Model form as printed: `C = C0 exp(−kt)`, C in mmol l⁻¹, t in min.
- DFG disappearance "followed indeed a first-order reaction model" [M].
- At **100 °C, pH 5.5**, when reaction order was itself estimated by minimising SS:
  **order = 1.6, 95% CI 0.83–2.4** [F]. Anchor: p. 1666, §3.1; illustrated in Fig. 1.
- Authors' own conclusion: "These results show that the use of simple kinetics is very limited."
- **No k value is printed for the simple-kinetics fit. Unreadable/absent — there is no number to take.**

---

## 4. THE REACTION SCHEMES — step-by-step, with the arrow directions

### 4.1 Scheme 1 (p. 1667) — "established reaction network"

Qualitative, not fitted. Primary routes solid, secondary dashed. Species: DFG, Gly, MG, Glc, Man,
Fru, AA, FA, 3-DG, 1-DG, Cₙ (n ≤ 6, "unidentified carbohydrate fragments"), STD (Strecker
degradation products), E1, E2. Caption verbatim:

> "E1 and E2 are unidentified key compounds involved in rate-determining steps that can be the
> Schiff's base, the cation form of the Schiff's base, the 1,2 enaminol or the 2,3-enaminol,
> respectively."

### 4.2 Scheme 2 = **Model 1 (M1)** (p. 1668) — via its printed differential equations

The differential equations printed on pp. 1667–1668 (eqs. 1–14) define M1 unambiguously. Note the
paper **misnumbers eq. (10) as "(19)"** — a printing error (see DEFECT box). Transcribed verbatim,
then read as arrows:

| eq. | printed equation | implied step(s) |
|---|---|---|
| (1) | d[DFG]/dt = −k1[DFG] − k2[DFG] − k3[DFG] | DFG consumed by 1, 2, 3 |
| (2) | d[E1]/dt = k1[DFG] − k4[E1] − k10[E1] − k11[E1] | **1: DFG → E1**; E1 consumed by 4, 10, 11 |
| (3) | d[E2]/dt = k2[DFG] − k7[E2] − k16[E2] | **2: DFG → E2**; E2 consumed by 7, 16 |
| (4) | d[MG]/dt = k3[DFG] | **3: DFG → MG** (M1: direct, no intermediate) |
| (5) | d[3DG]/dt = k4[E1] − k5[3DG] − k6[3DG] | **4: E1 → 3-DG (+Gly)**; **5: 3-DG → Cₙ**; **6: 3-DG → FA** |
| (6) | d[1DG]/dt = k7[E2] − k8[1DG] − k9[1DG] | **7: E2 → 1-DG (+Gly)**; **8: 1-DG → Cₙ**; **9: 1-DG → AA** |
| (7) | d[Man]/dt = k10[E1] − k12[Man] | **10: E1 → Man (+Gly)**; **12: Man → Glc** |
| (8) | d[Glc]/dt = k11[E1] + k12[Man] − k13[Glc] | **11: E1 → Glc (+Gly)**; **13: Glc → Cₙ** (M1) |
| (9) | d[Fru]/dt = k16[E2] | **16: E2 → Fru** (hypothesis B; hypothesis A is 16: E1 → Fru) |
| **(19)** ← misnumbered, should be (10) | d[FA]/dt = k6[3DG] | 6: 3-DG → FA |
| (11) | d[AA]/dt = k9[1DG] | 9: 1-DG → AA |
| (12) | d[Gly]/dt = k3[DFG] + k4[E1] + k10[E1] + k11[E1] + k7[E2] + k16[E2] − k14[Gly][Cₙ] | glycine released by 3, 4, 10, 11, 7, 16; consumed by **14** |
| (13) | d[Cₙ]/dt = k3[DFG] + k5[3DG] + k8[1DG] + k13[Glc] − k14[Gly][Cₙ] | Cₙ from 3, 5, 8, 13; consumed by 14 |
| (14) | d[Mel]/dt = k14[Gly][Cₙ] | **14: Gly + Cₙ → Mel** (the only bimolecular step) |

M1 also carried a **step 15: MG → FA + AA**, added during iteration. Verbatim result (p. 1668):
> "it appeared that the estimation of this parameter lead to an **indeterminate result for all the
> studied systems**, which means that this step is not important for the model. Also it suggests
> that the organic acids formation does not result from MG degradation."

### 4.3 Scheme 3 = **Model 2 (M2)** (p. 1671) — the scheme Table 3 belongs to

**Anchor: Scheme 3, p. 1671 (PDF p. 7).** Read off the page raster. Caption verbatim:

> "Second proposed kinetic model, Model 2 (M2), based on the established reaction network. E₁ and
> E₂ are unidentified key compounds involved in rate-determining steps that can be the Schiff's
> base, the cation form of the Schiff's base, the 1,2 enaminol or the 2,3-enaminol, respectively.
> Glycine (Gly); methylglyoxal (MG); Glc; Man; fructose (Fru); acetic acid (AA); formic acid (FA);
> 3-deoxyosone (3-DG); 1-deoxyosone (1-DG); unidentified carbohydrate fragments (Cₙ); melanoidins
> (Mel). **A** (secondary route) and **B** (primary route) hypothesis for fructose formation."

**M2 written out as arrows — THIS IS THE AUTHORITATIVE MAPPING FOR TABLE 3:**

```
                                       k12
                              Man  ─────────────►  Glu
                              ▲                     │
                        k10  /                      │ k13
                    (−Gly)  /            k11        │
                           /   ┌────────────────────┘
                          /    ▼        (−Gly)
       k16              E1 ──────────────►  3-DG  ──── k6 ────►  FA
  Fru ◄──── E1           ▲       k4 (−Gly)    │
   (hypothesis A)        │                    │ k5
                         │ k1                 ▼
                        DFG  ──── k3 ────►   Cₙ  ──── k15 ────►  MG + FA + AA
                         │                    ▲
                         │ k2 (−Gly)          │ k8
                         ▼                    │
       k16              E2  ──── k7 ────►  1-DG  ──── k9 ────►  AA
  Fru ◄──── E2                   (−Gly)
   (hypothesis B)

                     Cₙ + Gly  ──── k14 ────►  Mel
```

**Step-by-step, M2 (each row is one arrow; ALL ARROWS ARE UNIDIRECTIONAL):**

| step | reaction | note |
|---|---|---|
| 1 | DFG → E₁ | 1,2-enolisation branch entry |
| 2 | DFG → E₂ (+ Gly) | 2,3-enolisation branch entry |
| 3 | DFG → Cₙ | retro-aldolisation; **changed from M1**, where step 3 went straight to MG |
| 4 | E₁ → 3-DG (+ Gly) | |
| 5 | 3-DG → Cₙ | colour-precursor route |
| 6 | 3-DG → FA | formic acid |
| 7 | E₂ → 1-DG (+ Gly) | |
| 8 | 1-DG → Cₙ | |
| 9 | 1-DG → AA | acetic acid |
| 10 | E₁ → Man (+ Gly) | parent-sugar formation |
| 11 | E₁ → Glc (+ Gly) | parent-sugar formation |
| 12 | Man → Glc | **sugar isomerisation, ONE DIRECTION ONLY** |
| 13 | Glc → 3-DG | **changed from M1**, where step 13 was Glc → Cₙ |
| 14 | Cₙ + Gly → Mel | **only second-order step** |
| 15 | Cₙ → MG + FA + AA | **changed from M1**, where step 15 was MG → FA + AA and was indeterminate |
| 16 | **hypothesis A:** E₁ → Fru; **hypothesis B:** E₂ → Fru | **Table 3 is hypothesis B ⇒ k16 = E₂ → Fru** |

### 4.4 REVERSE RATE CONSTANTS — the explicit negative finding

> **THERE ARE NO REVERSE RATE CONSTANTS IN THIS PAPER.**

Every arrow in Scheme 3 (and Scheme 2) is single-headed. There is no k₋₁, no k₋₄, no
DFG → Glc + Gly step, no Schiff-base back-reaction, and no equilibrium constant anywhere in the
paper. Concretely:

1. **DFG is the starting reactant**, so the Amadori-forming step does not appear at all and
   therefore neither does its reverse.
2. **Parent-sugar formation (Glc, Man, Fru) is modelled as FORWARD steps 10, 11 and 16 out of the
   enolisation intermediates E₁/E₂** — not as a reversal of Amadori rearrangement. This is a
   deliberate mechanistic choice, argued in Part I §3.2.4.
3. **The closest thing to a reverse constant is k12 (Man → Glc)**, and the authors state in print
   why the opposite direction was deleted rather than fitted (p. 1667):
   > "From a previous study²⁰ … it was reported that Man isomerised into Glc at a rate of
   > **18 × 10⁻³ min⁻¹** while the reverse was approximately half (**7 × 10⁻³ min⁻¹**), as well as
   > the degradation of Man into Cₙ (n ≤ 6) compounds. … Concerning Glc the rate of its degradation
   > into Cₙ (n ≤ 6) compounds was found to be higher than its rate of isomerization into Man.
   > **As a result the degradation step of Glc into Man was neglected.**"
   Those two numbers (18 × 10⁻³ and 7 × 10⁻³ min⁻¹) are **[C] cited**, from De Bruin, J. M.,
   Ph.D. Thesis, University of Technology, Delft, 1986 (ref. 20), for **alkaline conditions** —
   they are not measurements of this system and must not be used as if they were.
4. The paper closes by explicitly *deferring* the reversibility question (p. 1676):
   > "This result in combination with the fact that 3-DG is the main precursor of carbohydrate
   > fragments involved in colour formation raises an important question about the importance of
   > ARP reversibility in Maillard reaction. **This question will be addressed in a following
   > paper.**"
   That following paper is **Martins & Van Boekel (2005) Food Chem. 90, 257–269**
   (`data/articles/martins2005.pdf`), which resolves it by model discrimination only and
   **also prints no numeric reverse rate constant** — its Scheme 4 collapses E₁ out by a
   steady-state assumption. I checked that file directly to close the question.

**Consequence for the declaration:** the FIT_HOLDOUT_DECLARATION line promising "Table 4.2.3
reverse rates (Amadori → parent sugar, 100/120 °C × pH 5.5/6.8)" as a **FIT** dataset that
"closes the declared reverse-rate structural gap with measured values" is **not supported by any
of the four journal papers in this wave, nor by martins2005.pdf.** Either the thesis contains a
table that never made it into any of the papers, or the declaration mis-describes what is there.
**Flagged for the orchestrator as a blocking item; it must be checked against the thesis PDF
before anything is fitted to it.**

---

## 5. TABLE 3 — THE RATE-CONSTANT TABLE (complete transcription)

**Anchor: Table 3, p. 1676 (PDF p. 12).** Title as printed:

> "Rate constants (10⁻²) ± 95% highest posterior density (HPD) as found by kinetic modelling for
> hypothesis B of Model 2 (M2)"

Column headers as printed: `Rate constant (min⁻¹)` | `A (100 °C, pH 5.5)` | `B (120 °C, pH 5.5)` |
`C (100 °C, pH 6.8)` | `D (120 °C, pH 6.8)`

**Units, exactly as printed:** the row label carries `min⁻¹`; the table title carries the scale
factor `(10⁻²)`. So **every tabulated number must be multiplied by 10⁻² to obtain min⁻¹.**
Columns 1–4 below are transcribed **exactly as printed** (i.e. in units of 10⁻² min⁻¹); the
`[Z]` columns give my conversion to bare min⁻¹ and are labelled as mine.

All values are **[F] — fitted by the authors** (multiresponse non-linear regression, determinant
criterion). None is a direct measurement.

### 5.1 As printed (units of 10⁻² min⁻¹; error is ±95% HPD)

| Rate constant (min⁻¹) | A (100 °C, pH 5.5) | B (120 °C, pH 5.5) | C (100 °C, pH 6.8) | D (120 °C, pH 6.8) |
|---|---|---|---|---|
| k1  | 0.19 ± 0.02 | 1.11 ± 0.07 | 0.57 ± 0.03 | 8.89 ± 0.83 |
| k2  | 0.10 ± 0.04 | 0.86 ± 0.40 | 1.56 ± 0.09 | 6.29 ± 0.66 |
| k3  | 0.18 ± 0.01 | 0.88 ± 0.28 | 1.55 ± 0.11 | 8.62 ± 0.36 |
| k4  | 20.14 ± 3.37 | 31.13 ± 9.97 | 7.94 ± 1.51 | 215.79 ± 50.08 |
| k5  | 1.38 ± 0.39 | 2.23 ± 1.17 | 9.07 ± 3.45 | 506.69 ± 62.87 |
| k6  | 0.19 ± 0.04 | 4.30 ± 0.69 | 2.74 ± 2.82 | 30.41 ± 2.17 |
| k7  | 60.17 ± 9.11 | 76.79 ± 10.65 | 21.25 ± 9.18 | 55.93 ± 2.14 |
| k8  | 5.90 ± 4.06 | 15.41 ± 1.71 | 0.00 ± 0.00 | 0.00 ± 0.00 |
| k9  | 3.93 ± 0.21 | 22.50 ± 1.81 | 190.85 ± 22.85 | 653.55 ± 94.08 |
| k10 | 11.31 ± 1.94 | 19.42 ± 6.84 | 7.07 ± 1.11 | 25.07 ± 5.28 |
| k11 | 6.42 ± 1.07 | 10.96 ± 4.25 | 11.31 ± 1.81 | 50.55 ± 10.63 |
| k12 | 0.39 ± 0.13 | 1.27 ± 0.24 | 0.08 ± 0.05 | 1.06 ± 0.17 |
| k13 | 0.73 ± 0.28 | 1.41 ± 0.26 | 0.22 ± 0.05 | 2.03 ± 0.23 |
| k14 | 0.12 ± 0.03 | 70.68 ± 3.93 | 0.34 ± 0.06 | 2.47 ± 0.99 |
| k15 | 1.45 ± 0.42 | 5.15 ± 0.99 | 1.59 ± 0.22 | 16.82 ± 5.91 |
| k16 | *(cell empty — not estimated)* | *(cell empty — not estimated)* | 1.34 ± 0.59 | 4.51 ± 1.72 |

The k16 cells at pH 5.5 are **blank in the original, not zero.** The reason is given in §3.4,
p. 1673: "this test is only relevant for the systems at pH 6.8, since at lower pH Fru was not
detected."

### 5.2 Same table with the step meaning attached, converted to bare min⁻¹ [Z]

Conversion by me: printed value × 10⁻². **[Z] on every value in the four numeric columns; the
step column is [M] (read from Scheme 3).**

| k | step (M2, hypothesis B) | A: 100 °C pH 5.5 (min⁻¹) [Z] | B: 120 °C pH 5.5 (min⁻¹) [Z] | C: 100 °C pH 6.8 (min⁻¹) [Z] | D: 120 °C pH 6.8 (min⁻¹) [Z] |
|---|---|---|---|---|---|
| k1  | DFG → E₁ | 1.9 × 10⁻³ | 1.11 × 10⁻² | 5.7 × 10⁻³ | 8.89 × 10⁻² |
| k2  | DFG → E₂ (+Gly) | 1.0 × 10⁻³ | 8.6 × 10⁻³ | 1.56 × 10⁻² | 6.29 × 10⁻² |
| k3  | DFG → Cₙ | 1.8 × 10⁻³ | 8.8 × 10⁻³ | 1.55 × 10⁻² | 8.62 × 10⁻² |
| k4  | E₁ → 3-DG (+Gly) | 2.014 × 10⁻¹ | 3.113 × 10⁻¹ | 7.94 × 10⁻² | 2.1579 |
| k5  | 3-DG → Cₙ | 1.38 × 10⁻² | 2.23 × 10⁻² | 9.07 × 10⁻² | 5.0669 |
| k6  | 3-DG → FA | 1.9 × 10⁻³ | 4.30 × 10⁻² | 2.74 × 10⁻² | 3.041 × 10⁻¹ |
| k7  | E₂ → 1-DG (+Gly) | 6.017 × 10⁻¹ | 7.679 × 10⁻¹ | 2.125 × 10⁻¹ | 5.593 × 10⁻¹ |
| k8  | 1-DG → Cₙ | 5.90 × 10⁻² | 1.541 × 10⁻¹ | 0.00 | 0.00 |
| k9  | 1-DG → AA | 3.93 × 10⁻² | 2.250 × 10⁻¹ | 1.9085 | 6.5355 |
| k10 | E₁ → Man (+Gly) | 1.131 × 10⁻¹ | 1.942 × 10⁻¹ | 7.07 × 10⁻² | 2.507 × 10⁻¹ |
| k11 | E₁ → Glc (+Gly) | 6.42 × 10⁻² | 1.096 × 10⁻¹ | 1.131 × 10⁻¹ | 5.055 × 10⁻¹ |
| k12 | Man → Glc | 3.9 × 10⁻³ | 1.27 × 10⁻² | 8 × 10⁻⁴ | 1.06 × 10⁻² |
| k13 | Glc → 3-DG | 7.3 × 10⁻³ | 1.41 × 10⁻² | 2.2 × 10⁻³ | 2.03 × 10⁻² |
| k14 | Cₙ + Gly → Mel | 1.2 × 10⁻³ | 7.068 × 10⁻¹ | 3.4 × 10⁻³ | 2.47 × 10⁻² |
| k15 | Cₙ → MG + FA + AA | 1.45 × 10⁻² | 5.15 × 10⁻² | 1.59 × 10⁻² | 1.682 × 10⁻¹ |
| k16 | E₂ → Fru | not estimated | not estimated | 1.34 × 10⁻² | 4.51 × 10⁻² |

### 5.3 DEFECT BOX — internal problems in Table 3 that a fitting wave must not paper over

1. **k14 unit defect.** k14 is the **only bimolecular** step (d[Mel]/dt = k14·[Gly]·[Cₙ], eq. 14),
   so its unit must be **l mmol⁻¹ min⁻¹** (or equivalent), not min⁻¹. Table 3's row label prints
   min⁻¹ for **all** constants including k14. **The printed unit for k14 is wrong.** Any transfer
   of k14 into the trunk must re-derive the unit from eq. (14) and must know the concentration
   basis (mmol/l) used in the fit. *(In the later Martins & Van Boekel 2005 Food Chem. 92 paper,
   the analogous bimolecular constant k1 IS correctly labelled `l mol⁻¹ min⁻¹` — confirming that
   the min⁻¹ label here is an error, and confirming that unit basis is mol vs. mmol
   inconsistent across the Martins series.)*
2. **k14 is wildly non-monotone in temperature at pH 5.5.** 0.12 (100 °C) → **70.68** (120 °C) is a
   **589-fold** jump for a 20 °C rise, against 0.34 → 2.47 (7.3-fold) at pH 6.8. That is a
   ~180 kJ mol⁻¹ apparent Ea at one pH and ~48 kJ mol⁻¹ at the other. This is almost certainly a
   fitting artefact (the pH 5.5 / 120 °C melanoidin response is small and the parameter poorly
   identified) or a typographic error. **Do not fit or Arrhenius-extrapolate k14 from this table.**
3. **k8 is estimated as exactly 0.00 ± 0.00 at both pH 6.8 conditions.** This is a *deleted*
   parameter, not a measured zero. It also explains the parameter counts in Tables 1 and 2:
   16 constants − k8 = 15 parameters at pH 6.8; at pH 5.5, k16 absent gives 15. A zero with a zero
   interval carries no information and must never be entered into a model as a measured rate.
4. **k5 at 120 °C / pH 6.8 = 506.69 × 10⁻² = 5.07 min⁻¹**, and k9 = 653.55 × 10⁻² = 6.54 min⁻¹.
   These are ~10³ times k1 in the same column. On a 5-min sampling grid these are far past the
   point where the data can constrain them; the printed HPD intervals (±12% and ±14%) are
   almost certainly over-confident. The authors themselves say only that "the rate constant of
   step 9 is the highest of the system."
5. **Eq. (10) is printed as "(19)"** on p. 1668 — a numbering error in the original.
6. **No Ea anywhere.** Two temperatures per pH is the *minimum* for a two-point Arrhenius slope,
   and the authors deliberately did not compute one. Any Ea derived from this table is [Z], has
   zero degrees of freedom, and inherits both endpoints' HPD width.
7. **The pH label is a room-temperature initial pH that then drifts**, by 0.19–0.21 units at
   pH 6.8 and by **0.8–0.9 units at pH 5.5**. The "pH 5.5" columns are therefore really a pH
   *trajectory* from 5.5 down to ~4.6, not a pH-5.5 measurement. This matters enormously for any
   pH-axis hold-out design.

### 5.4 The authors' own reading of Table 3 (verbatim, pp. 1675–1676) — [M] interpretive claims

- "at pH 5.5 step 1 prevailed to steps 2 and 3, especially when the temperature was increased. At
  higher pH (6.8) on the other hand step 2 gained importance, which is evident at lower
  temperature, suggesting that 2,3-enolization becomes more relevant by increasing the pH."
- "independently of the reaction conditions, deoxyosones formation prevail to sugars formation.
  The rate constant for step 4 is always higher than for step 10 and for step 11."
- "Fru formation is only a minor step from the 2,3-enolization pathway. The rate constant for step
  7 is always higher than for step 16."
- "At lower pH Man formation was favored towards Glc and no other sugar was detected, whereas at
  higher pH not only Fru was formed but also Glc was formed in higher amounts than Man."
- "steps 6 and 9 are the most significant in formic and acetic acid formation, respectively."
- "Not only the rate constant for the degradation of 1-DG into acetic acid (step 9) prevailed to
  the degradation of 1-DG into carbohydrate fragments (Cₙ) (step 8), in particular at pH 6.8, but
  also under these conditions the rate constant of step 9 is the highest of the system."
- "As the pH increased 1-DG is no longer involved in carbohydrate fragments responsible for colour
  formation, but mainly in acetic acid formation. Under these conditions 3-DG becomes the main
  precursor in colour formation, through step 5."

---

## 6. TABLE 1 — model discrimination at pH 5.5

**Anchor: Table 1, p. 1675 (PDF p. 11).** Title as printed: "Model discrimination tests for the
systems studied at pH 5.5". Footnote a as printed: "Model 1 (M1) presented in Scheme 2; Model 2
(M2) presented in Scheme 3."

Column headers as printed: `System` | `Model ᵃ` | `Parameters` | `SS` | `n` | `AIC_c` | `Δ_AICc` | `PPB`
(no units are printed on any column; SS is in the multiresponse determinant-criterion sense,
PPB is a log posterior probability on Athena's scale).

| System | Model | Parameters | SS | n | AIC_c | Δ_AICc | PPB |
|---|---|---|---|---|---|---|---|
| A (100 °C) | M1 | 14 | 5.14 | 100 | −261.92 | 0 | 34.11 |
| A (100 °C) | M2 | 15 | 5.58 | 100 | −250.96 | 10.96 | 32.28 |
| B (120 °C) | M1 | 13 | 17.01 | 90 | −117.22 | 11.50 | 14.78 |
| B (120 °C) | M2 | 15 | 14.05 | 90 | −128.72 | 0 | 17.65 |

All [F]. Authors' reading (p. 1672): "At 100 °C, M2 is less supported whereas at 120 °C M1 can be
discarded, which was also supported by the PPB values." Their caveat, verbatim:
> "It should be taken into account that under these conditions (pH 5.5), for the studied heating
> period, the DFG degradation rate is quite small as well as the amount of products formed, which
> gives the models more flexibility to fit the experimental data."

**Note the direct contradiction between the two rows of this table** — at 100 °C the *preferred*
model is M1 and at 120 °C it is M2, on the same pH, with the same chemistry. The authors do not
resolve it; they rely on the pH 6.8 result (Table 2) instead.

---

## 7. TABLE 2 — model discrimination at pH 6.8 (M1 vs M2 × hypothesis A vs B)

**Anchor: Table 2, p. 1676 (PDF p. 12).** Title as printed: "Model discrimination tests for the
systems studied at pH 6.8". Footnotes as printed: "* Model 1 (M1) presented in Scheme 2; Model 2
(M2) presented in Scheme 3." and "ᵃ Indeterminate (low trust region)."

Column headers as printed: `System` | `Hypotheses` | `Model *` | `Parameters` | `SS` | `n` |
`AIC_c` | `Δ_AICc` | `PPB`

| System | Hypotheses | Model | Parameters | SS | n | AIC_c | Δ_AICc | PPB |
|---|---|---|---|---|---|---|---|---|
| C (100 °C) | A | M1 | 14 | 50.05 | 99 | −32.60 | 173.32 | 16.20 |
| C (100 °C) | A | M2 | 14 | 8.69 | 99 | −205.92 | 0 | 19.80 |
| D (120 °C) | A | M1 | 13 | 49.85 | 99 | −35.69 | 127.86 | 15.78 |
| D (120 °C) | A | M2 | 15 | 12.96 | 99 | −163.55 | 0 | 21.25 |
| C (100 °C) | B | M1 | 12 | 754.67 | 99 | 230.67 | 423.88 | **Indt.ᵃ** |
| C (100 °C) | B | M2 | 15 | 9.61 | 99 | −193.21 | 0 | 29.09 |
| D (120 °C) | B | M1 | 13 | 49.08 | 99 | −37.22 | 123.62 | 14.59 |
| D (120 °C) | B | M2 | 15 | 13.32 | 99 | −160.85 | 0 | 21.77 |

All [F]. "Indt." = indeterminate (low trust region) — **unreadable as a number by design; there is
no value.**

Authors' reading (p. 1673), verbatim:
> "According to the Akaike criterion, independently of the chosen hypothesis A or B, M2 always
> performed better than M1. The results for M1 indicated that this model could be discarded
> (ΔAIC ≥ 10). These findings are confirmed by the obtained PPB values that were always higher in
> M2. … When comparing hypotheses A and B for M2 the PPB values as well as the AIC values of
> hypothesis B are higher than those of hypothesis A. These results suggest that **fructose was
> mainly formed through the intermediate E₂**."

**Internal inconsistency worth flagging:** the authors say "the AIC values of hypothesis B are
higher than those of hypothesis A" and treat that as support for B. But for AIC, **lower is
better**, and M2's AIC_c is −205.92 under A vs −193.21 under B at 100 °C, and −163.55 under A vs
−160.85 under B at 120 °C — i.e. by AIC, **hypothesis A fits better at both temperatures.** The
PPB values do favour B (29.09 > 19.80 at 100 °C; 21.77 > 21.25 at 120 °C). **Their conclusion is
carried by PPB alone, and their sentence about AIC has the sign of the criterion backwards.**
Table 3 is nonetheless the hypothesis-B fit. This is a real, quotable defect and directly weakens
any downstream claim that "E₂ → Fru is the established topology."

---

## 8. FIGURES — inventory and digitised values

**Digitisation policy for this dossier:** Part II's Figs. 2 and 4–7 replot the **same experimental
points as Part I** (see the double-counting note in §0). Per-point digitisation of those response
panels is therefore done in the **Part I dossier** (`martins2003b_extraction.md` §5–§11), where
the data are primary. Here I record only what is unique to Part II: panel inventory, axis ranges,
and the mass-balance figure, which is Part II-only.

### 8.1 Fig. 1 — simple-kinetics comparison
**Anchor: Fig. 1, p. 1666 (PDF p. 2).** Caption verbatim: "Simple kinetics analysis of
N-(1-deoxy-D-fructos-1-yl)-glycine (DFG) thermal degradation at 100 °C, pH 5.5. Comparison of a
first-order (—), estimated order of 1.6 (---) and second-order (– –) plot." No numeric values are
recoverable beyond the order = 1.6 already transcribed in §3. **No axis tick labels legible at
page-raster resolution: unreadable.**

### 8.2 Fig. 2 — M1 fit at 120 °C, pH 5.5 (A) and 6.8 (B)
**Anchor: Fig. 2, pp. 1669–1670 (PDF pp. 5–6).** Caption verbatim: "Model 1 fit (lines) to
experimental data (dots) of DFG thermal degradation at 120 °C: pH 5.5 (A) and 6.8 (B). Glycine
(Gly); methylglyoxal (MG); Glc; Man; fructose (Fru); acetic acid (AA); formic acid (FA);
3-deoxyosone (3-DG); 1-deoxyosone (1-DG)." **Superseded by M2 (Figs. 4–7) — no values extracted;
this is a rejected model's fit.** The authors' verbatim residual diagnosis is the useful content:
- pH 5.5: "the model fitted the data reasonably well … Glycine was slightly overestimated, while
  formic acid and melanoidins were underestimated at the beginning and overestimated at the end of
  heating period."
- pH 6.8: "there was clearly a miss fit for the organic acids formation, namely for formic acid,
  as well as for Glc. Both compounds were underestimated. Also a lack of fit was observed for MG
  formation. It was overestimated at the beginning and underestimated as the reaction proceeded."

### 8.3 Fig. 3 — MASS BALANCE (Part II-only content)
**Anchor: Fig. 3, p. 1671 (PDF p. 7).** Caption verbatim: "Mass balance: evolution of each
intermediate towards the reactant initial concentration in heated
N-(1-deoxy-D-fructos-1-yl)-glycine (DFG) at 100 °C, pH 6.8 (A) and 120 °C, pH 6.8 (B). Glycine
(Gly); Deoxyosones (DG); Sugars (Sug); organic acids (OA); methylglyoxal (MG); melanoidins (Mel)."

Stacked bar chart, y-axis "%" 0–120 with a dotted reference line at 100%. Stack order in the
legend: DFG, Gly, DG, Sug, OA, MG, Mel. **Only the stack TOTAL is legible at page-raster
resolution; the individual segment boundaries below ~5% are not.** Per-segment values: **unreadable**.

Totals digitised from the page raster, all **[fig]** (read from PDF p. 7 raster; my read precision
is ±3 percentage points, limited by bar-top thickness):

| panel | time (min) | total mass balance (%) [fig] |
|---|---|---|
| A (100 °C, pH 6.8) | 0 | 100 |
| A | 15 | ~81 |
| A | 30 | ~78 |
| A | 45 | ~80 |
| A | 60 | ~85 |
| A | 90 | ~89 |
| A | 120 | ~93 |
| A | 150 | ~96 |
| A | 180 | ~100 |
| B (120 °C, pH 6.8) | 0 | 100 |
| B | 5 | ~90 |
| B | 10 | ~72 |
| B | 15 | ~78 |
| B | 30 | ~85 |
| B | 45 | ~92 |
| B | 60 | ~95 |
| B | 90 | ~98 |
| B | 120 | ~100 |

Authors' verbatim statement, which is the citable form of this result [M]:
> "at the initial stage of the reaction the products identified and quantified in the present study
> do not count for the total DFG degradation. However, as the degradation reaction proceeded,
> **within experimental error, 100% was reached**. … A possible explanation for the observed gap is
> that, as the reaction proceeds besides DFG other compounds formed during its degradation become
> reactants themselves with the formation of the same end products."

**No carbon-balance (as opposed to mole-balance) figure is given anywhere in this paper.** The
authors' "mass balance" is explicitly "evolution of each intermediate towards the reactant initial
concentration", i.e. a **molar** balance against 10 mmol/l DFG, with no carbon-number weighting.
Requests for a carbon balance: **not present in this paper.**

### 8.4 Figs. 4–7 — M2 fits (the accepted model)
Captions verbatim:
- **Fig. 4**, p. 1672: "Model 2 (M2) fit (lines) to experimental data (dots) of
  N-(1-deoxy-D-fructos-1-yl)-glycine (DFG) thermal degradation at **100 °C and pH 5.5**. Glycine
  (Gly); methylglyoxal (MG); Glc; Man; acetic acid (AA); formic acid (FA); 3-deoxyosone (3-DG);
  1-deoxyosone (1-DG)."
- **Fig. 5**, p. 1673: same, "at **120 °C and pH 5.5**", legend lists "glucose (Glu); mannose (Man)".
- **Fig. 6**, p. 1674: same, "at **100 °C and pH 6.8**".
- **Fig. 7**, p. 1675: same, "at **120 °C and pH 6.8**".

**Fig. 7 panel inventory and axis ranges** (read from PDF p. 11 raster; this is the only one I
rastered in full):

| panel | series (legend as printed) | y-axis label & range | x-axis range | terminal / peak value [fig] |
|---|---|---|---|---|
| "DFG" | DFG (○), Gly (▲), MG (×) | mM, 0–12 (ticks every 2) | Time (min) 0–150 | DFG → ~0 by 30 min; Gly plateau ~7.6 mM; MG plateau ~2.9 mM |
| "Sugars" | Glc (◆), Man (△), Fru (●) | mM, 0–1 (ticks 0.2) | 0–150 | Glc peak ~0.85 mM at ~10 min, → ~0.20 at 120 min; Man ~0.40 peak → ~0.17; Fru rises to ~0.30 |
| "Organic Acids" | AA (▲), FA (◇) | mM, 0–8 (ticks 2) | 0–150 | AA → ~5.5 mM; FA → ~3.0 mM |
| "Deoxyosones" | 3-DG (◆), 1-DG (○) | mM, 0–0.08 (ticks 0.02) | 0–150 | 3-DG peak ~0.045 mM at ~5 min; 1-DG peak ~0.060 mM at ~2 min; both ≈0 by 40 min |
| "Melanoidins" | (×, unlabelled series) | mM, 0–3 (ticks 1) | 0–150 | plateau ~2.4 mM by ~45 min |

All **[fig]**, from the PDF p. 11 raster; read precision ~±5% of full scale. These are the same
measurements as Part I Figs. 1–8 at 120 °C / pH 6.8.

Authors' summary of the M2 fit quality (p. 1669): "A major improvement in the organic acids fit was
observed, as well as in the sugars and MG formation. Independently of the reaction conditions the
model seems to fit the experimental data quite well."

---

## NEW-PARAMETER TABLE (consolidated)

Every rate constant below is **[F]** (author-fitted). "as printed" = the Table 3 cell; "min⁻¹" =
my ×10⁻² conversion, tagged [Z]. Anchor for every rate-constant row: **Table 3, p. 1676**.

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| k1 (DFG → E₁) | 0.19 ± 0.02 → 1.9 × 10⁻³ [Z] | 10⁻² min⁻¹ ("min⁻¹" with title factor 10⁻²) | 100 °C, pH 5.5 init., 0.1 M phosphate, [DFG]₀ 10 mmol/l | Table 3, p. 1676 | [F] |
| k1 | 1.11 ± 0.07 → 1.11 × 10⁻² [Z] | 10⁻² min⁻¹ | 120 °C, pH 5.5 | Table 3, p. 1676 | [F] |
| k1 | 0.57 ± 0.03 → 5.7 × 10⁻³ [Z] | 10⁻² min⁻¹ | 100 °C, pH 6.8 | Table 3, p. 1676 | [F] |
| k1 | 8.89 ± 0.83 → 8.89 × 10⁻² [Z] | 10⁻² min⁻¹ | 120 °C, pH 6.8 | Table 3, p. 1676 | [F] |
| k2 (DFG → E₂ + Gly) | 0.10 ± 0.04 / 0.86 ± 0.40 / 1.56 ± 0.09 / 6.29 ± 0.66 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k3 (DFG → Cₙ) | 0.18 ± 0.01 / 0.88 ± 0.28 / 1.55 ± 0.11 / 8.62 ± 0.36 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k4 (E₁ → 3-DG + Gly) | 20.14 ± 3.37 / 31.13 ± 9.97 / 7.94 ± 1.51 / 215.79 ± 50.08 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k5 (3-DG → Cₙ) | 1.38 ± 0.39 / 2.23 ± 1.17 / 9.07 ± 3.45 / 506.69 ± 62.87 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k6 (3-DG → FA) | 0.19 ± 0.04 / 4.30 ± 0.69 / 2.74 ± 2.82 / 30.41 ± 2.17 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k7 (E₂ → 1-DG + Gly) | 60.17 ± 9.11 / 76.79 ± 10.65 / 21.25 ± 9.18 / 55.93 ± 2.14 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k8 (1-DG → Cₙ) | 5.90 ± 4.06 / 15.41 ± 1.71 / **0.00 ± 0.00** / **0.00 ± 0.00** | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] — the zeros are deleted parameters, not measurements |
| k9 (1-DG → AA) | 3.93 ± 0.21 / 22.50 ± 1.81 / 190.85 ± 22.85 / 653.55 ± 94.08 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k10 (E₁ → Man + Gly) | 11.31 ± 1.94 / 19.42 ± 6.84 / 7.07 ± 1.11 / 25.07 ± 5.28 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k11 (E₁ → Glc + Gly) | 6.42 ± 1.07 / 10.96 ± 4.25 / 11.31 ± 1.81 / 50.55 ± 10.63 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k12 (Man → Glc) | 0.39 ± 0.13 / 1.27 ± 0.24 / 0.08 ± 0.05 / 1.06 ± 0.17 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k13 (Glc → 3-DG) | 0.73 ± 0.28 / 1.41 ± 0.26 / 0.22 ± 0.05 / 2.03 ± 0.23 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k14 (Cₙ + Gly → Mel) | 0.12 ± 0.03 / **70.68 ± 3.93** / 0.34 ± 0.06 / 2.47 ± 0.99 | 10⁻² min⁻¹ **(unit WRONG — bimolecular)** | A / B / C / D | Table 3, p. 1676 | [F] — see DEFECT box items 1–2 |
| k15 (Cₙ → MG + FA + AA) | 1.45 ± 0.42 / 5.15 ± 0.99 / 1.59 ± 0.22 / 16.82 ± 5.91 | 10⁻² min⁻¹ | A / B / C / D | Table 3, p. 1676 | [F] |
| k16 (E₂ → Fru, hyp. B) | not estimated / not estimated / 1.34 ± 0.59 / 4.51 ± 1.72 | 10⁻² min⁻¹ | C / D only (pH 6.8) | Table 3, p. 1676 | [F] |
| DFG apparent reaction order | 1.6, 95% CI 0.83–2.4 | dimensionless | 100 °C, pH 5.5 | §3.1, p. 1666 | [F] |
| Man → Glc isomerisation rate | 18 × 10⁻³ | min⁻¹ | alkaline conditions, monosaccharide isomerisation study | §3.3, p. 1667 (cites De Bruin 1986, ref. 20) | **[C]** — do not use as this system's value |
| Glc → Man isomerisation rate | 7 × 10⁻³ ("approximately half") | min⁻¹ | as above | §3.3, p. 1667 (De Bruin 1986) | **[C]** |
| relative enolisation rates Glc : Man | 1.0 and 0.5 | h⁻¹ | Isbell et al. 1971 (ref. 22) | §3.3, p. 1668 | **[C]** |
| melanoidin ε at 470 nm (used to convert A₄₇₀ → mM) | 0.64 ± 0.03 | l mmol⁻¹ cm⁻¹ | glucose/glycine | Part I §2.10; cites `martins2003c.pdf` | **[C]** within this paper |
| 3-DG / 1-DG ratio | 3.2 | dimensionless | 100 °C, pH 6.8 | Part I §3.2.1, p. 1659 | [M] (Part I) |
| Model discrimination, pH 5.5 | ΔAICc: M1 0 / M2 10.96 (100 °C); M1 11.50 / M2 0 (120 °C) | dimensionless | pH 5.5 | Table 1, p. 1675 | [F] |
| Model discrimination, pH 6.8, hyp. A | ΔAICc: M1 173.32 / M2 0 (100 °C); M1 127.86 / M2 0 (120 °C) | dimensionless | pH 6.8 | Table 2, p. 1676 | [F] |
| Model discrimination, pH 6.8, hyp. B | ΔAICc: M1 423.88 / M2 0 (100 °C); M1 123.62 / M2 0 (120 °C) | dimensionless | pH 6.8 | Table 2, p. 1676 | [F] |
| PPB, M2 hyp. B | 29.09 (100 °C) / 21.77 (120 °C) | dimensionless (Athena log-posterior scale) | pH 6.8 | Table 2, p. 1676 | [F] |
| Mass balance closure | ~72–81% minimum at 10–30 min, recovering to ~100% at the end | % of [DFG]₀ (molar) | 100 and 120 °C, pH 6.8 | Fig. 3, p. 1671 | **[fig]** (totals only; segments unreadable) |
| **Reverse rate constants** | **NONE PRESENT** | — | — | Schemes 2 & 3, Table 3 | **negative finding, high confidence** |
| **Activation energies (Ea)** | **NONE PRESENT** | — | — | whole paper | **negative finding, high confidence** |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Blocking issue first.** Amendment 1 of the declaration already assigns **FIT** to "Martins 2003
thesis Table 4.2.3 reverse rates (Amadori → parent sugar, 100/120 °C × pH 5.5/6.8)" on the stated
rationale that it "closes the declared reverse-rate structural gap with measured values."
**I could not find any such reverse rate constant in this paper (the journal version of that
chapter), nor in Part I, nor in the 2005 Food Chem. 90 follow-up that explicitly took up the
reversibility question.** The declaration also attributes ΔAICc = 276 to a "Table 5.1", which
actually belongs to `martins2005.pdf` Table 1, where the winning hypothesis A beats the *fully
irreversible* hypothesis C by **287.46**, not 276 (276 is the A-vs-B gap, B being a *differently*
reversible topology). **Recommendation: freeze that FIT row until someone opens
edepot.wur.nl/121418 and reads Table 4.2.3 directly.** If the thesis table turns out to be this
paper's Table 3 under a different number, then the "reverse rates" description is simply wrong and
the row should be re-scoped or dropped.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| Table 3, **columns C + D (pH 6.8, 100 and 120 °C)**, steps k1, k2, k3, k4, k6, k7, k9, k10, k11 | **FIT** | pH | pH 6.8 is the trunk's fixed pH, so these are in-condition. These nine steps are the well-identified ones (HPD ≤ ~30% of the estimate, monotone in T, chemically load-bearing: both enolisation branches, retro-aldol, both deoxyosone→acid steps, both parent-sugar steps). |
| Table 3, **columns A + B (pH 5.5, 100 and 120 °C)**, same nine steps | **HOLD-OUT** | pH | Genuine off-pH extrapolation for a pH-fixed trunk, mirroring the Hofmann convention and the declaration's existing Ch. 6 cut. **Circularity risk: LOW but non-zero** — same lab, same apparatus, same regression code, and (critically) the same 2 × 2 experiment already reported as Part I. If Part I's time-courses are FIT anywhere, this arm is contaminated. |
| Table 3, **k5, k8, k12, k13, k15** (all columns) | **quarantine — do not fit, do not score** | — | k8 is an explicit deleted-zero at pH 6.8; k5 at 120 °C/pH 6.8 is 5.07 min⁻¹ with a sampling grid that cannot constrain it; k12/k13 are the surviving halves of pairs whose partners were deleted on cited (not measured) grounds; k15 lumps three products into one arrow. |
| Table 3, **k14 (all columns)** | **quarantine — do not fit, do not score, do NOT unit-convert silently** | — | Printed unit (min⁻¹) contradicts the printed rate law (eq. 14, bimolecular); the pH-5.5 temperature ratio is 589× vs 7.3× at pH 6.8. Any use requires a unit re-derivation and an explicit decision about the mmol/l vs mol/l basis. |
| Table 3, **k16 (E₂ → Fru), columns C + D only** | **diagnostic_only** | — | Its topological justification (hypothesis B over A) rests on a PPB comparison the authors themselves mis-described in AIC terms (see §7). Structurally interesting, evidentially weak. |
| **Any Ea derived from Table 3 columns A→B or C→D** | **[Z] prior at most; never a scored target** | temperature | Two points, zero residual degrees of freedom, and the authors deliberately did not publish one. If the trunk needs a Martins Ea, take it from `martins2005.pdf` Table 2 (5-temperature Arrhenius fit with printed ± intervals) instead — that is a strictly better source and is a different paper. |
| **Tables 1 and 2 (discrimination statistics)** | **structure-selection evidence only; never a fit target and never a scored hold-out** | — | These are model-comparison numbers, not physical observables. Usable to justify a topology choice (M2 over M1: ΔAICc 123.6–423.9 at pH 6.8, unambiguous). **Do not** use them to justify "E₂ → Fru" without carrying the §7 caveat. |
| **Fig. 3 mass-balance totals** | **diagnostic_only** | — | [fig]-only, segments unreadable, ±3 pp read precision, and it is a molar not a carbon balance. Useful as a sanity check that a simulated DFG run closes to ~100% at long times; useless as a scored target. |
| **Figs. 4–7 response time-courses** | **DO NOT ADD — already covered as Part I** | — | Same experiment, published twice. Assign the time-course role once, in the Part I dossier, and cross-reference it here. Assigning it in both places is a rule-1 disjointness violation of exactly the kind Amendment 2 had to correct for Cerny 2007. |

**Cut-axis summary:** the clean cut is **pH (6.8 FIT / 5.5 HOLD-OUT)**, not temperature. The
temperature axis has only two levels per pH and both are needed to identify any of these
constants, so cutting on T would leave one arm unidentifiable. The pH cut is also the one axis
where the paper contains a real physical surprise for a pH-fixed trunk — the 2,3-/1,2-enolisation
branching ratio k2/k1 flips from 0.53 (100 °C, pH 5.5) to 2.74 (100 °C, pH 6.8) [Z, from Table 3].

**One further honesty flag on the pH cut:** the "pH 5.5" arm drifts by **0.8–0.9 pH units** during
the run (Part I Fig. 7 and §3.1.6), versus 0.19–0.21 units at pH 6.8. A hold-out scored at
"pH 5.5" is really scored at a moving pH averaging nearer 5.1. If the trunk has any pH dependence
at all, this must be modelled or the hold-out will fail for the wrong reason.
