# Martins & Van Boekel 2005 (10.1016/j.foodchem.2004.08.013) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/martins2005b.pdf` (12 pp.). Read method: both — `pdftotext -layout`
text layer for the body, plus page rasters of PDF pp. 6, 7 and 9 (journal pp. 442, 443, 445) read
with the Read tool to recover **every minus sign and exponent sign in Tables 1–4**, which the text
layer silently drops (see the DEFECT note in §4).

---

## 0. PAPER IDENTITY — **MATCHES**, but the sibling keys around it DO NOT

> **IDENTITY, STATED IN BOLD AT THE TOP AS REQUIRED.**
> **`data/articles/martins2005b.pdf` IS Martins & Van Boekel, "Kinetics of the glucose/glycine
> Maillard reaction pathways: influences of pH and reactant initial concentrations", Food Chemistry
> 92 (2005) 437–448.** This is correct as expected. **This is the pH-ladder paper.**
>
> **BUT THE WAVE BRIEF'S SIBLING ASSIGNMENTS WERE WRONG AND MUST BE RESTATED HERE:**
> - **`martins2003.pdf` is Part II (Kinetic analysis), Carbohydr. Res. 338 (2003) 1665–1678 — this
>   is where the DFG rate-constant table lives.** It was *expected* to be `martins2003c.pdf`. It is not.
> - **`martins2003c.pdf` is a DIFFERENT paper entirely**: "Melanoidins extinction coefficient in the
>   glucose/glycine Maillard reaction", Food Chemistry 83 (2003) 135–142.
> - **`martins2003b.pdf` is Part I** (Reaction mechanism), Carbohydr. Res. 338 (2003) 1651–1663.
>
> **This paper's own reference list is the cleanest disambiguator in the whole set** and is quoted
> verbatim here, because it settles the letter-suffix confusion for anyone reading later:
> > "Martins, S. I. F. S., & Van Boekel, M. A. J. S. (**2003a**). Kinetic modelling … **Part II –
> > Kinetic analysis. Carbohydrate Research, 338, 1665–1678.**"
> > "Martins, S. I. F. S., & Van Boekel, M. A. J. S. (**2003b**). **Extinction coefficient of
> > melanoidins** … *Food Chemistry, 83, 135–142.*"
> > "Martins, S. I. F. S., Marcelis, A. T. M., & Van Boekel, M. A. J. S. (**2003c**). … **Part I –
> > Reaction mechanism. Carbohydrate Research, 338, 1651–1663.**"
>
> **These letters are the OPPOSITE of our repo keys, and they are permuted differently again in
> Martins & Van Boekel (2005) Food Chem. 90 (`martins2005.pdf`), where 2003a = Part I,
> 2003b = Part II, 2003c = extinction coefficient. Never resolve a Martins 2003 citation by its
> letter suffix. Resolve it by journal + volume + page.**

| field | value |
|---|---|
| authors | Sara I.F.S. Martins, Martinus A.J.S. Van Boekel* |
| title | "Kinetics of the glucose/glycine Maillard reaction pathways: influences of pH and reactant initial concentrations" |
| venue | Food Chemistry |
| volume/pages/year | 92 (2005) 437–448 |
| DOI | 10.1016/j.foodchem.2004.08.013 |
| received / accepted | Received 17 May 2004; accepted 2 August 2004 |
| affiliation | Product Design and Quality Management, Dept. of Agrotechnology and Food Sciences, Wageningen University, P.O. Box 8129, 6700 EV Wageningen, The Netherlands |
| PDF character | born-digital (Elsevier / Acrobat Distiller 5.0.5 for Windows); text layer present but **systematically drops the minus sign in exponents and in negative numbers** — e.g. "1.1 × 10⁻⁵ ± 4 × 10⁻⁶" renders as "1.1 · 105 ± 4 · 106", and "−8.43 ± 0.91" renders as "8.43 ± 0.91". **Every sign in Tables 1–4 was re-verified against the page raster.** |
| corresponding author | tiny.vanboekel@wur.nl |
| funding / acknowledgement | Fundação para a Ciência e Tecnologia (Portugal); "The authors also thank **Ruud van Gorp** for performing the heating experiment at constant pH." |

### PRE-EXISTING-COVERAGE FLAG — DOUBLE-COUNTING HAZARD

`docs/reference/FIT_HOLDOUT_DECLARATION.md` (read-only; Amendment 1, 2026-08-28, lines ~209–215)
declares roles for "Martins 2003 Wageningen thesis (edepot.wur.nl/121418)" Tables 4.2.3, 4.1.1,
3.3.1 and Ch. 6. **These four journal papers are the thesis chapters.**

**My mapping proposal for this paper (proposal, NOT a verified fact — I have not opened the thesis
PDF):**

| declared thesis item | my proposed mapping | confidence | verdict |
|---|---|---|---|
| **Ch. 6 pH ladder** — "k at pH **4.8/5.5/6.0/6.8/7.5** + **fitted pH exponents**", role **"HOLD-OUT except pH 5.5/6.8 columns (which are the already-FIT conditions)"**, cut axis pH | **THIS PAPER.** Its Table 2 is the k table at exactly pH **4.8, 5.5, 6.0, 6.8, 7.5**, and its Table 3 is the **fitted pH exponents (pD) with 95% CIs**. The declared pH list is a verbatim match; no other Martins paper contains a five-pH ladder. | **VERY HIGH** | **DIRECT DOUBLE-COUNT. Same measurements, re-published.** Ch. 6 = this paper. Do not re-declare. |
| Table 4.2.3 "reverse rates (Amadori → parent sugar), 100/120 °C × pH 5.5/6.8" | not in this paper — this paper is **single-temperature (100 °C only)** and has **no reverse rate constant**. Its step 3 (D-Fru → D-Glu) is the reverse *isomerisation*, not an Amadori reversal. | high | see the blocking note below |
| Table 4.1.1 glycine-release yields | maps to `martins2003b.pdf` (Part I) Table 1 | high | n/a |
| Table 3.3.1 melanoidin ε | maps to `martins2003c.pdf` Table 1; **this paper only *cites* it** ("0.64 ± 0.03 l mmol⁻¹ cm⁻¹", §2.2) | high | n/a |

**Two precision notes on the existing Ch. 6 declaration row, which I recommend amending:**
1. The row's carve-out — "**except pH 5.5/6.8 columns (which are the already-FIT conditions)**" —
   needs a source pointer. The pH 6.8 column of Table 2 is the **same experiment** as Table 1's
   "pH drop ≤ 1" column and as Table 4's "1:1" column: **all three print identical values for all
   ten constants** (verified cell-by-cell in §7 below). It is *one* measurement printed *three
   times in one paper*. Whichever declaration row claims it, the other two must be excluded.
2. The row says "cut axis: pH". Correct — **but the pH 4.8 column contains four zeros
   (k3 = k5 = k10 = 0.00) which are DELETED parameters, not measured zeros.** A hold-out scored on
   that column will score four cells that carry no information. See DEFECT box item 3.

**Blocking issue inherited from this wave:** Amendment 1 also assigns **FIT** to "thesis
Table 4.2.3 reverse rates (Amadori → parent sugar)" on the rationale that it "closes the declared
reverse-rate structural gap with **measured** values." **I could not find any such reverse rate
constant in this paper, in Part I, in Part II, or in Martins & Van Boekel (2005) Food Chem. 90
(`martins2005.pdf`), which is the paper that explicitly took up the reversibility question and
resolved it by model discrimination *without ever printing a numeric reverse constant*.** Details
in `martins2003_extraction.md` §4.4. **Recommend freezing that FIT row until the thesis PDF is
read directly.**

---

## 1. ONE-PARAGRAPH VERDICT

This is the highest-value paper of the four for a pH-aware trunk. It gives a **complete
10-parameter rate-constant table at five initial pH values (4.8, 5.5, 6.0, 6.8, 7.5)** on a single
system (glucose + glycine, 0.2 M equimolar, 0.1 M phosphate, 100 °C), each constant with a
**95% HPD interval**; a **second 10-parameter table at three glucose:glycine molar ratios**
(1:1, 2:1, 1:2); a **third table comparing a free-drifting-pH run against a pH-stat-controlled
run**; and — the unique deliverable — **fitted pH-dependence exponents (pD) and intercepts
(log kₑ) for all ten steps, each with a 95% confidence interval**, together with the explicit
power-law equation they were regressed in. **It is SINGLE-TEMPERATURE (100 °C only) and therefore
contains NO activation energies** — the Ea's for this same 10-step model live in the sister paper
`martins2005.pdf` (Food Chem. 90, 257–269) Table 2. It contains **no reverse rate constants**;
the only backward-looking step is k3 (fructose → glucose), one leg of the sugar-isomerisation
pair. A sensitivity analysis (Fig. 10) is presented for k1, k3, k6 and k8 but **only as unlabelled
plots — no numeric sensitivities are printed.**

---

## 2. SYSTEM DEFINITION — verbatim

| variable | value as printed | anchor |
|---|---|---|
| reactants | "**Equimolar solutions of glucose and glycine (0.2 mol/l)**" | §2.1, p. 439 |
| ratio experiments | "the concentrations of glucose and glycine were changed to give **molar ratios of 2:1 and 1:2**", with the equimolar system as control | §3.3, p. 445 |
| buffer | "prepared in **phosphate buffer (10 ml, 0.1 mol/l)**" | §2.1, p. 439 |
| pH | "set to have an initial pH **4.8, 5.5, 6.0, 6.8 and 7.5 at 20 °C**" | §2.1, p. 439 |
| **hot-pH correction** | **EXPLICITLY NOT APPLIED, and the authors say so.** Verbatim: "It is common knowledge that the pH decreases with temperature because of increased water dissociation. **It is therefore stressed that the initial pH values reported are as measured at room temperature.**" | §2.1, p. 439 |
| pH control | **Normally uncontrolled and drifting.** One control experiment used a **pH-stat**: "The pH stat controller monitored the pH and, when it dropped below its set point, an auto burette (**ABU901 Radiometer Copenhagen**) added a solution of **NaOH (1 mol/l)** to the reaction mix. … When the pH was kept constant, this was done at **the pH the solution had at 100 °C**; the corresponding initial pH was 6.8 at room temperature." | §2.1, p. 439 |
| **pH-drop data filter — CRITICAL** | "to avoid an interfering affect of pH drop on the reaction rate inhibition, **only data points were taken for which the pH drop was not higher than 1 unit**" | Table 2 footnote, p. 442; §3.2, p. 441 |
| filtration | "The solutions were **filtered (0.2 µm, Schleicher & Schuell)**" | §2.1, p. 439 |
| temperature | **100 °C only**, oil bath | §2.1, p. 439 |
| vessel | "screw-capped glass tubes (**Schott, 16 × 160 mm**)"; the pH-stat run used a septum-sealed vessel sampled with "a **Hamilton SampleLock (10 ml)**" | §2.1, p. 439 |
| headspace / atmosphere | **not stated** | — |
| volume | 10 ml | §2.1, p. 439 |
| agitation | **not stated** | — |
| quench | "samples were taken and **immediately cooled in ice**, prior to analyses"; pH-stat run: "into a glass tube and **immediately cooled in ice**" | §2.1, p. 439 |
| duration | "heated at 100 °C **for up to 4 h**" (pH study); ratio study runs extend to ~180–200 min | §3.1, p. 440 |
| temperature monitoring | "The temperature was monitored by a **Consort R305**." | §2.1, p. 439 |
| replication (n) | "Each reaction mixture was prepared, heated and analysed, **at least in duplicate**." | §2.1, p. 439 |
| error definition, parameters | "**± 95% HPDᵃ** … ᵃ **Highest posterior density**" (Tables 1 and 4); "**Estimates ± 95% confidence intervals**" (Table 3) | Tables 1, 3, 4 |
| error definition, data | **not stated** for the response figures | — |

### Analytes and analytical methods

> "The following compounds were identified and quantified: **glycine, D-glucose, D-fructose,
> N-(1-deoxy-D-fructos-1-yl)glycine, 1-deoxy-2,3-hexodiulose, 3-deoxy-2-hexosulose, formic acid,
> acetic acid, methylglyoxal and 5-(hydroxymethylfurfural).** The used methodology has been
> described previously (Martins & Van Boekel, 2004)." (§2.2, p. 439)

**So the analytical methods are inherited by citation** from `martins2005.pdf` (Food Chem. 90),
which in turn inherits them from Part I. See `martins2003b_extraction.md` §2 for the full
chromatographic detail (Biochrom 20 amino-acid analyser; 1,2-diaminobenzene quinoxaline
derivatisation for the α-dicarbonyls; ION-300 for sugars and acids; Dionex CarboPac PA100 for
Glc/Man/Fru separation).

**Melanoidins, verbatim (§2.2, p. 439):**
> "The melanoidins concentration was calculated from the measured absorbance at **470 nm**, using
> the extinction coefficient **0.64 ± 0.03 l mmol⁻¹ cm⁻¹** (Martins & Van Boekel, 2003b), and the
> resulting concentration reflects **the amount of glucose molecules incorporated in the
> melanoidins**."
**[C] within this paper — the primary is `martins2003c.pdf`.** Note the population caveat recorded
in `martins2003c_extraction.md` §3.4 item 1: that ε was measured on the >3500 Da fraction, which
carries only 12–18% of the colour, yet it is applied here to the total A₄₇₀.

**HMF was measured and then discarded, verbatim (§3.1, p. 440):**
> "In contrast, the formation of 5-hydroxymethylfurfural (HMF) **increased with decreasing pH**;
> however, its amounts were **one order of magnitude lower (i.e. in the µmolar range)** than the
> other detected reaction products, and **this compound was therefore neglected in the further
> analysis.**"
**No HMF numbers are published anywhere in the paper. Unquantified.**

### Fitting method, verbatim (§2.3, p. 439)

> "The proposed kinetic model (Scheme 1) was translated into a mathematical model by deriving
> differential equations for each reaction step, where the rate constants are the parameters to be
> estimated. The software package **Athena Visual Workbench** (www.athenavisual.com) was used for
> numerical integration of the differential equations, as well as for parameter estimation. The
> model parameters were estimated by non-linear regression, using the **determinant criterion**,
> that is to minimise the determinant of the matrix of cross-products of the various responses, the
> so called dispersion matrix. **The determinant criterion replaces the familiar least-squares
> objective, because the multiresponse method requires a different statistical approach** (Stewart,
> Caracotsios, & Sørensen, 1992). The goodness of fit test is installed in the used software
> package. It gives a sampling probability by which the adequacy of the model can be judged and was
> based on replicate experiments (Stewart, Shon, & Box, 1998)."

**No goodness-of-fit statistic, no SS, no AIC and no PPB value is printed anywhere in this paper.**
Unlike Part II, this paper reports **no model-discrimination table at all** — model selection was
done in the sister paper (`martins2005.pdf` Table 1). **Explicit negative finding.**

---

## 3. SCHEME 1 — THE 10-STEP MODEL, STEP BY STEP

**Anchor: Scheme 1, p. 439 (PDF p. 3).** Caption verbatim:
> "Kinetic model for the glucose/glycine Maillard reaction. D-Glucose (D-Glu); glycine (Gly);
> D-fructose (D-Fru); N-(1-deoxy-D-fructos-1-yl)–glycine (DFG); 3-deoxglucosone (3-DG);
> 1-deoxyglucosone (1-DG); methylglyoxal (MG); acetic acid (AA); formic acid (FA); melanoidins
> (Mel)."

**This model is imported wholesale from `martins2005.pdf` Scheme 4** — §1, p. 438: "Recently, the
influence of temperature on the glucose/glycine Maillard reaction was studied and a comprehensive
kinetic model was proposed, **which is reproduced here as Scheme 1** (Martins & Van Boekel, 2004)."
**It is NOT re-derived here. This paper only re-estimates its parameters at new pH values.**

**Scheme 1 written out as arrows:**

```
                                   k10
                       D-Glu ───────────────►  FA + AA
                         │
              k2  ┌──────┴──────┐  k3
        D-Fru ◄───┤             ├───► D-Glu        (reversible isomerisation pair)
                  └─────────────┘
                                   k1
                       D-Glu + Gly ─────────►  DFG

                       DFG ──── k4 ────►  3-DG + Gly   (1,2-enolisation)
                       DFG ──── k7 ────►  1-DG + Gly   (2,3-enolisation)
                       DFG ──── k6 ────►  MG   (+ Gly) (retro-aldolisation)

                      3-DG ──── k5 ────►  FA
                      1-DG ──── k8 ────►  AA

               3-DG + Gly ──── k9 ────►  Mel
```

| step | reaction | printed unit | authors' own description |
|---|---|---|---|
| **1** | D-Glu + Gly → DFG | **l mol⁻¹ min⁻¹** | "the initial Maillard reaction step 1" (§4, p. 447) |
| **2** | D-Glu → D-Fru | min⁻¹ | "sugar isomerisation" (§3.2, p. 443); "the one that describes sugar isomerisation (k2)" |
| **3** | D-Fru → D-Glu | min⁻¹ | "the isomerisation step of fructose into glucose (k3)" (§3.4, p. 447) — **the reverse leg of the isomerisation pair** |
| **4** | DFG → 3-DG + Gly | min⁻¹ | "the **1,2-enaminol route**, with the formation of 3-deoxyglucosone (step 4)" (§4, p. 447) |
| **5** | 3-DG → FA | min⁻¹ | "3-DG also showed low pH-dependence in **formic acid** formation (step 5)" (§4, p. 447) |
| **6** | DFG → MG (retro-aldol) | min⁻¹ | "the **DFG retro-aldolisation** step (k6)"; "the formation of methylglyoxal from the Amadori compound (step 6)" (§3.4, §4) |
| **7** | DFG → 1-DG + Gly | min⁻¹ | "the **2,3-enaminol** route (step 7)"; "k7 (1-DG formation from DFG)" (§3.4, §4) |
| **8** | 1-DG → AA | min⁻¹ | "the formation of **acetic acid** from 1-deoxyglucosone (k8)" (§3.4, p. 447) |
| **9** | 3-DG + Gly → Mel | min⁻¹ (**unit defect — bimolecular**) | "3-DG degraded (preferably into **melanoidin** formation) by reaction with glycine (step 9)" (§4, p. 447) |
| **10** | D-Glu → FA + AA | min⁻¹ | "**glucose degradation into the organic acids**" — described as "steps 2, 3 and 10" being sugar isomerisation and degradation (§4, p. 447) |

**Only ONE constant is printed with a second-order unit (k1, `l mol⁻¹ min⁻¹`). k9 is also
bimolecular by its own description ("3-DG + Gly → Mel") but is printed as `min⁻¹`.** See DEFECT
box item 1.

**Note the k-numbering is COMPLETELY DIFFERENT from Part II's 16-step DFG model.** In Part II, k9
is 1-DG → AA and k1 is DFG → E₁. Here k9 is 3-DG + Gly → Mel and k1 is Glu + Gly → DFG. **A
cross-paper "k9" comparison is meaningless. Always carry the model with the constant.**

---

## 4. TABLE 2 — THE pH-LADDER RATE-CONSTANT TABLE (complete transcription)

**Anchor: Table 2, p. 442 (PDF p. 6).** Title as printed:
> "Rate constants (k) estimation by the proposed glucose/glycine kinetic model (Scheme 1)"

Column headers as printed: `K` | `pH 4.8` | `pH 5.5` | `pH 6.0` | `pH 6.8` | `pH 7.5`
Row labels as printed: `1 (l mol⁻¹ min⁻¹)`, then `2 (min⁻¹)` … `10 (min⁻¹)`.
Footnote as printed:
> "Samples heated at 100 °C, in phosphate buffer (0.1 M); to avoid an interfering affect of pH drop
> on the reaction rate inhibition, only data points were taken for which the pH drop was not higher
> than 1 unit."

**Note: Table 2 carries NO error-type footnote.** Tables 1 and 4 both say "± 95% HPDᵃ / ᵃ Highest
posterior density"; Table 2 prints "±" values with no label. **By continuity with Tables 1 and 4
(and because the pH 6.8 column of Table 2 is numerically identical to Table 1's "pH drop ≤ 1"
column and Table 4's "1:1" column, errors included), these are certainly also 95% HPD intervals —
but the table itself does not say so. Recorded as a printing omission.**

**All values [F] — fitted by the authors.** Read from the PDF p. 6 raster; the text layer drops
every exponent minus sign.

| k | pH 4.8 | pH 5.5 | pH 6.0 | pH 6.8 | pH 7.5 |
|---|---|---|---|---|---|
| **1** (l mol⁻¹ min⁻¹) | 7.9 × 10⁻⁷ ± 2 × 10⁻⁸ | 2.3 × 10⁻⁶ ± 5 × 10⁻⁸ | 4.9 × 10⁻⁶ ± 1 × 10⁻⁷ | 1.1 × 10⁻⁵ ± 4 × 10⁻⁶ | 1.8 × 10⁻⁵ ± 8 × 10⁻⁷ |
| **2** (min⁻¹) | 8.8 × 10⁻⁶ ± 1 × 10⁻⁶ | 5.5 × 10⁻⁵ ± 4 × 10⁻⁶ | 2.3 × 10⁻⁴ ± **3 × 10⁻³** | 8.9 × 10⁻⁴ ± 4 × 10⁻⁵ | 2.1 × 10⁻³ ± 7 × 10⁻⁵ |
| **3** (min⁻¹) | **0.00** | **0.00** | 2.8 × 10⁻³ ± 2 × 10⁻³ | 5.9 × 10⁻³ ± 8 × 10⁻⁴ | 1.4 × 10⁻² ± 9 × 10⁻⁴ |
| **4** (min⁻¹) | 1.5 × 10⁻³ ± 7 × 10⁻⁵ | 3.1 × 10⁻³ ± 1 × 10⁻⁴ | 4.5 × 10⁻³ ± 3 × 10⁻⁴ | 6.5 × 10⁻³ ± 3 × 10⁻⁴ | 1.5 × 10⁻² ± 1 × 10⁻³ |
| **5** (min⁻¹) | **0.00** | 6.3 × 10⁻³ ± 2 × 10⁻³ | 2.0 × 10⁻² ± 2 × 10⁻³ | 4.1 × 10⁻³ ± 5 × 10⁻⁴ | 1.9 × 10⁻² ± 1 × 10⁻² |
| **6** (min⁻¹) | 5.0 × 10⁻⁴ ± 8 × 10⁻⁵ | 9.1 × 10⁻⁴ ± 1 × 10⁻⁴ | 2.7 × 10⁻³ ± 1 × 10⁻⁴ | 5.2 × 10⁻³ ± 2 × 10⁻⁴ | 1.1 × 10⁻² ± 9 × 10⁻⁴ |
| **7** (min⁻¹) | 7.7 × 10⁻⁴ ± 7 × 10⁻⁵ | 4.0 × 10⁻⁴ ± 1 × 10⁻⁴ | 3.1 × 10⁻³ ± 3 × 10⁻⁴ | 1.3 × 10⁻² ± 4 × 10⁻⁴ | 1.5 × 10⁻² ± 1 × 10⁻³ |
| **8** (min⁻¹) | 9.9 × 10⁻² ± 9 × 10⁻³ | 6.2 × 10⁻² ± 2 × 10⁻² | 1.0 × 10⁻¹ ± 1 × 10⁻² | 1.3 × 10⁺⁰ ± 2 × 10⁻¹ | 4.5 × 10⁻¹ ± 4 × 10⁻² |
| **9** (min⁻¹) | 5.0 × 10⁻⁵ ± 4 × 10⁻⁶ | 2.1 × 10⁻⁴ ± 1 × 10⁻⁵ | 2.4 × 10⁻⁴ ± 2 × 10⁻⁵ | 9.1 × 10⁻⁴ ± 2 × 10⁻⁵ | 5.3 × 10⁻⁴ ± 5 × 10⁻⁵ |
| **10** (min⁻¹) | **0.00** | 2.1 × 10⁻⁵ ± 3 × 10⁻⁶ | 6.3 × 10⁻⁵ ± **6 × 10⁻⁴** | 1.0 × 10⁻⁴ ± 5 × 10⁻⁵ | 3.1 × 10⁻⁴ ± 6 × 10⁻⁵ |

### 4.1 DEFECT BOX — problems in Table 2 that a fitting wave must not paper over

1. **k9 unit defect.** Step 9 is "3-DG + Gly → Mel" by the authors' own description (§4, p. 447:
   "3-DG degraded … by reaction with glycine (step 9)"), i.e. **bimolecular**, so its unit must be
   a second-order one (l mol⁻¹ min⁻¹ or l mmol⁻¹ min⁻¹). **Table 2 prints `min⁻¹`.** Only k1 gets
   a second-order label. Any transfer of k9 must re-derive the unit from the rate law and must fix
   the concentration basis. *(Compare `martins2003.pdf` Table 3, where the analogous melanoidin
   constant k14 has exactly the same defect.)*
2. **Concentration-basis inconsistency across the Martins series.** k1's unit here is
   **l mol⁻¹ min⁻¹** (mol basis) while every concentration in the figures is plotted in **mmol/l**.
   Part II uses mmol/l throughout with no second-order label at all. **Anyone porting these
   constants must resolve mol vs mmol explicitly — a factor of 1000 on every bimolecular step.**
3. **THREE CELLS PRINTED AS "0.00" ARE DELETED PARAMETERS, NOT MEASURED ZEROS:** k3 at pH 4.8 and
   5.5; k5 at pH 4.8; k10 at pH 4.8. They carry **no error bar at all** (not even ±0), and the
   corresponding species were not detected at those pHs. A hold-out scored on the pH 4.8 column
   would be scoring 3 of 10 cells that contain no information. **Chemically these are consistent
   (no fructose → no k3; negligible browning/acids below pH 6) — but they are absences, not zeros.**
4. **Two error bars are LARGER THAN THEIR OWN ESTIMATE:**
   - k2 at pH 6.0: **2.3 × 10⁻⁴ ± 3 × 10⁻³** — the interval is **13× the estimate** and spans zero
     and beyond.
   - k10 at pH 6.0: **6.3 × 10⁻⁵ ± 6 × 10⁻⁴** — the interval is **9.5× the estimate**.
   Both are visible on the page raster and are not text-layer artefacts. Given that k2 at pH 5.5 is
   5.5 × 10⁻⁵ ± 4 × 10⁻⁶ and at pH 6.8 is 8.9 × 10⁻⁴ ± 4 × 10⁻⁵ (both ~5% relative), the pH 6.0
   entries look like **typographic errors — most plausibly 3 × 10⁻⁵ and 6 × 10⁻⁶.** I have NOT
   corrected them; they are transcribed as printed and flagged here. **Neither cell should be used
   with its printed uncertainty.**
5. **NON-MONOTONE CELLS that contradict the paper's own "overall, the rate constants increased with
   increasing initial pH" (§3.2, p. 441):**
   - **k5:** 0.00 → 6.3 × 10⁻³ → **2.0 × 10⁻²** → **4.1 × 10⁻³** → 1.9 × 10⁻². The pH 6.0 value is
     **4.9× the pH 6.8 value.** This is the single worst offender and it is exactly why Table 3
     gives k5 the near-zero, hugely uncertain exponent pD = 0.09 ± 1.10.
   - **k7:** 7.7 × 10⁻⁴ → **4.0 × 10⁻⁴** → 3.1 × 10⁻³ → … The pH 5.5 value is **below** the
     pH 4.8 value.
   - **k8:** 9.9 × 10⁻² → 6.2 × 10⁻² → 1.0 × 10⁻¹ → **1.3 × 10⁺⁰** → **4.5 × 10⁻¹**. The pH 6.8
     value is **2.9× the pH 7.5 value**, and pH 5.5 dips below pH 4.8.
   - **k9:** 9.1 × 10⁻⁴ (pH 6.8) → **5.3 × 10⁻⁴** (pH 7.5) — a 42% **decrease** at the highest pH.
   The authors' blanket "overall, the rate constants increased with increasing initial pH" is true
   on average and false in detail for k5, k7, k8 and k9. **Quote both.**
6. **The pH label is a room-temperature INITIAL pH that then drifts downward during the run.** The
   authors mitigate by discarding data beyond a 1-unit drop, but that means **each pH column is a
   trajectory, and the columns are not equally trimmed** — the higher-pH runs drop faster (Fig. 4)
   and so lose more late data. **The effective mean pH of the "pH 7.5" column is materially below
   7.5.**
7. **k8 at pH 6.8 is printed `1.3 × 10⁺⁰`** — the only positive exponent in the table, ~2 orders
   above its neighbours in the same row's pH 6.0 cell. It is reproduced identically in Tables 1 and
   4 for the same experiment, so it is at least internally consistent, but it makes k8 the fastest
   step in the model by ~100× and is the reason the paper says "the limiting factor on the amount
   of acetic acid seems to be the formation of 1-DG" rather than its consumption.

---

## 5. TABLE 3 — THE FITTED pH EXPONENTS (the priority deliverable)

**Anchor: Table 3, p. 443 (PDF p. 7).** Title as printed:
> "Estimates ± 95% confidence intervals of parameters **kₑ** and **pD** from pH rate profiles"

Column headers as printed: `Rate constant k_obs` | `pD` | `Log(kₑ)`

### 5.1 THE EXACT FUNCTIONAL FORM THE EXPONENTS WERE FITTED IN — quoted verbatim

This is the part the wave brief asked for verbatim. Reproduced exactly, in the order the paper
presents it (§3.2, pp. 442–443):

Step 1, the regression that was actually run — a straight line in log-space:
> "As a next step, the pH-dependence was analysed quantitatively by making a so-called **pH rate
> profile**, i.e. a **plot of log(k_obs) vs. pH**. A general equation for a pH rate profile is:
>
> **log(k_obs) = log a + b × pH,**
>
> where **a** and **b** are regression coefficients. The coefficient **a** reflects one or more
> elementary rate constants and this coefficient is therefore called **kₑ**, while **b** reflects
> the pH-dependence, which will be called **pD**."

Step 2, the same relationship re-expressed as the power law that Table 3's parameters plug into:
> "The pH-dependence of observed rate constants can then be expressed as:
>
> **k_obs = kₑ (10^{pH · pD}).**
>
> This equation basically fulfils the same function as the Arrhenius equation does for the
> temperature-dependence of rate constants."

**⚠ TYPOGRAPHIC DISCREPANCY BETWEEN THE ABSTRACT AND THE BODY, FLAGGED AS THE SPEC REQUIRES.**
The Abstract (p. 437) prints the same equation with the exponent's factors transposed and the
multiplication dot placed differently:
> Abstract, verbatim: "a quantitative relationship was derived: **k_obs = kₑ(10^{pD · pH})**, in
> which k_obs is the estimated rate constants from experiments, **kₑ** an expression for the
> elementary reaction, and **pD** the parameter describing the pH-dependence."
> Body §3.2, p. 442, verbatim: "**k_obs = kₑ(10^{pH · pD})**."
These are **mathematically identical** (pD·pH = pH·pD); the difference is cosmetic. **I flag it
only because a careless reader could mistake "10^{pD}·pH" for a different functional form.** The
authoritative reading is unambiguous from the log-linear regression it was derived from:
**log₁₀(k_obs) = log₁₀(kₑ) + pD × pH.**

### 5.2 What the parameters mean, and the reference values — verbatim (§3.2, pp. 441–443)

The diagnostic framework, quoted in full because it is what makes pD interpretable:
> "For specific acid–base catalysis, it can be derived that the observed rate constant, k_obs,
> should have the following relationship (Loudon, 1991):
> **k_obs = k₁[H⁺]  for acid catalysis;**
> **k_obs = k₂ K_w/[H⁺]  for base catalysis;**
> where k₁ and k₂ are elementary rate constants for the reaction under study, and K_w is the water
> dissociation constant, so a plot of log(k_obs) vs. pH should have a **slope of −1 for acid
> catalysis and +1 for base-catalysis**."

> "For specific acid catalysis, kₑ reflects the elementary rate constant for the acid-catalysed step
> and b, or pD **should be −1** in that case. For pure specific base-catalysed reactions, kₑ
> reflects the **product of the water dissociation constant K_w and the elementary rate constant**
> for the base-catalysed step, while b or pD **should equal +1**. If there is **no pH-dependence
> (pD 0)**, then the observed rate constant reflects the elementary rate constant. The two
> parameters kₑ and pD **need to be determined by experiment**, as we have done here."

Loudon's five candidate forms, quoted verbatim from §1, p. 438, because they define the space the
authors chose *not* to use:
> "1. The observed rate constant is independent of pH: **k_obs = c1**.
> 2. The observed rate constant is directly proportional to H⁺ ions: **k_obs = c2[H⁺]**.
> 3. The observed rate constant is inversely proportional to H⁺ ions (or equivalently directly
>    proportional to [OH⁻]): **k_obs = c3/[H⁺]**.
> 4. The observed rate constant is independent of pH at low pH and directly proportional to pH at
>    high pH: **k_obs = c4[H⁺]/(c5 + [H⁺])**.
> 5. The observed rate constant is inversely proportional to pH at low pH and independent at high
>    pH: **k_obs = c6/(c7 + [H⁺])**."
**The authors fitted form 2/3 in its log-linear generalisation (a free slope pD), not the saturating
forms 4 or 5.** That is a modelling choice worth carrying: it forbids any pH plateau by construction.

### 5.3 THE TABLE — complete transcription

**All values [F] — fitted by the authors, by linear regression of log₁₀(k_obs) on pH across the
five pH columns of Table 2.**

**Units, exactly as printed and as they must be read:**
- **pD** — the table prints **no unit**. It is dimensionless: it is the slope of
  log₁₀(k) vs. pH, i.e. **decades of rate per pH unit**. *(Note this is the negative of the
  conventional "n" in k ∝ [H⁺]ⁿ: a pD of +0.50 means k ∝ [H⁺]^−0.50.)*
- **Log(kₑ)** — the table prints **no unit**, and this is a genuine defect: kₑ carries the unit of
  the constant it belongs to, so **row 1's log(kₑ) = −8.43 is log₁₀ of a quantity in
  l mol⁻¹ min⁻¹, while rows 2–10 are log₁₀ of quantities in min⁻¹.** The column mixes two
  dimensions. See DEFECT box item 2 below.
- Errors are **"± 95% confidence intervals"** per the table title — note this differs from
  Tables 1, 2 and 4, which use **95% HPD**. Two different interval types in one paper.

| Rate constant k_obs | step (Scheme 1) | **pD** (dimensionless; decades per pH unit) | **Log(kₑ)** (log₁₀; **unit = that of the step, l mol⁻¹ min⁻¹ for row 1, min⁻¹ otherwise**) |
|---|---|---|---|
| **1** | D-Glu + Gly → DFG | **0.50 ± 0.14** | **−8.43 ± 0.91** |
| **2** | D-Glu → D-Fru | **0.88 ± 0.27** | **−9.14 ± 1.68** |
| **3** | D-Fru → D-Glu | **0.46 ± 0.46** | **−5.36 ± 3.26** |
| **4** | DFG → 3-DG + Gly | **0.35 ± 0.10** | **−4.46 ± 0.60** |
| **5** | 3-DG → FA | **0.09 ± 1.10** | **−2.59 ± 7.54** |
| **6** | DFG → MG | **0.51 ± 0.14** | **−5.75 ± 0.85** |
| **7** | DFG → 1-DG + Gly | **0.60 ± 0.52** | **−6.24 ± 3.25** |
| **8** | 1-DG → AA | **0.40 ± 0.60** | **−3.12 ± 2.40** |
| **9** | 3-DG + Gly → Mel | **0.40 ± 0.38** | **−6.03 ± 2.35** |
| **10** | D-Glu → FA + AA | **0.54 ± 0.36** | **−7.56 ± 2.35** |

**Sign convention check [Z]:** every pD is **positive**, i.e. every rate increases with pH, i.e.
every step is net base-catalysed to some degree. On the authors' own scale (−1 = pure acid
catalysis, +1 = pure base catalysis), all ten sit between 0 and +1 — "a mix between pure acid and
pure base catalysis."

**Reconstruction check [Z], to verify I have the equation the right way round.** Taking k1 and
predicting k_obs at pH 6.8 from Table 3:
log₁₀(k) = −8.43 + 0.50 × 6.8 = −5.03 → k = 9.4 × 10⁻⁶ l mol⁻¹ min⁻¹.
Table 2 prints **1.1 × 10⁻⁵** at pH 6.8. **Agreement to 15%** — the equation and the parameter
signs are correct as transcribed. Same check on k2: −9.14 + 0.88 × 6.8 = −3.16 →
6.9 × 10⁻⁴ min⁻¹, against Table 2's **8.9 × 10⁻⁴**. **Agreement to 22%.** Both well inside the
CIs. **The transcription is verified.**

### 5.4 DEFECT BOX — problems in Table 3

1. **Five of the ten exponents are not significantly different from zero at 95%:** k3
   (0.46 ± 0.46, bound touches zero), **k5 (0.09 ± 1.10, wildly indeterminate)**, k7 (0.60 ± 0.52),
   k8 (0.40 ± 0.60), k9 (0.40 ± 0.38, bound touches zero at 0.02). The authors acknowledge some of
   this: "**One rate constant does not really depend on pH (k5), while the slopes for k3 and k8 do
   not significantly differ from zero (or are too imprecise to be estimated well).**" **They do not
   flag k7 and k9, whose intervals are equally or more marginal.** Only **k1, k2, k4 and k6** have
   exponents that are convincingly non-zero.
2. **The Log(kₑ) column mixes units** (l mol⁻¹ min⁻¹ for row 1, min⁻¹ for the rest) with **no unit
   printed at all**. Anyone reconstructing k from kₑ must know which row they are on.
3. **k5's parameters are effectively meaningless as printed:** pD = 0.09 ± 1.10 and
   log(kₑ) = −2.59 ± 7.54. The kₑ interval spans **fifteen orders of magnitude**. This is a direct
   consequence of the non-monotone k5 column in Table 2 (DEFECT §4.1 item 5). **Do not use k5's
   pH parameters for anything.**
4. **n = 5 points per regression, and for k3, k5 and k10 the zero-valued cells cannot be
   log-transformed at all** — so those three exponents are fitted on **4 or even 3 points**. The
   paper never states the per-row n. **This is why their CIs are the widest.**
5. **The regression is unweighted in log space** as far as one can tell, i.e. it ignores the very
   different HPD widths in Table 2. Not stated either way.
6. **Table 3 is derived from Table 2. They are NOT independent evidence.** Fitting a trunk to
   Table 2 and then validating it against Table 3 (or vice versa) is circular.

### 5.5 Fig. 6 — the pH–rate profiles behind Table 3

**Anchor: Fig. 6, p. 443 (PDF p. 7).** Caption verbatim: "pH–rate profiles for the estimated rate
constants describing reactions depicted in Scheme 1. k₁ (◆), k₂ (■), k₃ (●), k₄ (◇), k₅ (□),
k₆ (▲), k₇ (○), k₈ (×), k₉ (△), k₁₀ (+)."

Three panels, all with y = `log(k_obs)` from **−7 to 0** (panel c: −7 to **+1**) and x = `pH` from
**4 to 8**:

| panel | series plotted | shape [fig] |
|---|---|---|
| (a) | k₁, k₂, k₃ | three near-parallel ascending lines; k₃ highest (log ≈ −2.8 → −1.9), k₂ middle (≈ −5.1 → −2.6), k₁ lowest (≈ −6.1 → −4.8) |
| (b) | k₄, k₅, k₆ | k₅ (□) nearly **flat** at log ≈ −1.9 to −1.7 — visibly the flat one; k₄ ascending ≈ −2.8 → −1.9; k₆ ascending ≈ −3.9 → −2.0 |
| (c) | k₇, k₈, k₉, k₁₀ | k₈ (×) highest, scattered around log ≈ −1.0 to +0.1 with a **visible outlier high point near pH 6.8**; k₇ ≈ −3.2 → −1.8; k₉ ≈ −4.3 → −3.3; k₁₀ ≈ −4.7 → −3.5 |

All **[fig]**, read from the PDF p. 7 raster; read precision ≈ ±0.15 log units. **These are simply
log₁₀ of the Table 2 cells — no new information, but they confirm visually that k₅ is the flat one
and that k₈'s pH 6.8 point is the outlier driving its wide CI.**

### 5.6 Literature comparison slopes — all [C], all derived by these authors from others' data

**Anchor: §3.2, pp. 443–444.** Verbatim, with provenance marked:
> "Dworschák and Örsi (1977) **determined** a slope of **0.38** for the pH rate profile describing
> the loss of **methionine** and **0.13** for **tryptophan** in the Maillard reaction with glucose.
> From results published by Pílková, Pokorný, and Davídek (1990) for browning of Heyns compounds,
> likewise, a slope of **0.27 can be derived**. From results of Ashoor and Zent (1984) a slope of
> **0.64 can be calculated** for browning. From results of Nicoli et al. (1993), a slope of
> **0.37 can be calculated** from browning data, a slope of **around 0.1** for the loss of glucose
> and glycine, and a slope of **0.1** for production of CO₂ in the Strecker reaction. Even though
> these literature results are not directly comparable to our results because of the different
> approach taken, **it is clear that the pH-dependence of most rate constants is less than that for
> pure base-catalysed reactions.**"

| slope | quantity | source | provenance |
|---|---|---|---|
| 0.38 | methionine loss | Dworschák & Örsi (1977) | **[C]** — determined by the original authors |
| 0.13 | tryptophan loss | Dworschák & Örsi (1977) | **[C]** |
| 0.27 | browning of Heyns compounds | Pílková et al. (1990) | **[C]-derived** — "can be derived", i.e. computed by Martins from others' data |
| 0.64 | browning | Ashoor & Zent (1984) | **[C]-derived** — "can be calculated" |
| 0.37 | browning | Nicoli et al. (1993) | **[C]-derived** |
| ~0.1 | glucose and glycine loss | Nicoli et al. (1993) | **[C]-derived** |
| 0.1 | CO₂ from Strecker reaction | Nicoli et al. (1993) | **[C]-derived** |

**None of these is a measurement from this paper. The four "can be derived / can be calculated"
values are Martins's own recomputations from published figures and inherit unknown digitisation
error.**

---

## 6. TABLE 1 — pH-STAT CONTROL vs. FREE-DRIFTING pH

**Anchor: Table 1, p. 442 (PDF p. 6).** Title as printed:
> "Rate constants (k) estimation ± 95% HPDᵃ intervals: system with observed pH drop vs. pH kept
> constant by NaOH addition"

Column headers as printed: (blank) | `pH drop ⩽ 1` | `pH constant`
Footnotes as printed: "Samples heated at 100 °C, pH 6.8; **n.a. not analysed (formation of
methylglyoxal)**." and "ᵃ Highest posterior density."

**All values [F].** Read from the PDF p. 6 raster.

| k (unit as printed) | pH drop ⩽ 1 | pH constant |
|---|---|---|
| **1** (l mol⁻¹ min⁻¹) | 1.1 × 10⁻⁵ ± 4 × 10⁻⁶ | 1.8 × 10⁻⁵ ± 1 × 10⁻⁶ |
| **2** (min⁻¹) | 8.9 × 10⁻⁴ ± 4 × 10⁻⁵ | 1.3 × 10⁻³ ± 1 × 10⁻⁴ |
| **3** (min⁻¹) | 5.9 × 10⁻³ ± 8 × 10⁻⁴ | 2.5 × 10⁻³ ± 1 × 10⁻³ |
| **4** (min⁻¹) | 6.5 × 10⁻³ ± 3 × 10⁻⁴ | 1.3 × 10⁻² ± 9 × 10⁻⁴ |
| **5** (min⁻¹) | 4.1 × 10⁻³ ± 5 × 10⁻⁴ | **1.3 × 10⁻¹ ± 1 × 10⁻²** |
| **6** (min⁻¹) | 5.2 × 10⁻³ ± 2 × 10⁻⁴ | **n.a.** |
| **7** (min⁻¹) | 1.3 × 10⁻² ± 4 × 10⁻⁴ | 2.0 × 10⁻² ± 2 × 10⁻³ |
| **8** (min⁻¹) | 1.3 × 10⁺⁰ ± 2 × 10⁻¹ | 8.7 × 10⁻¹ ± 5 × 10⁻² |
| **9** (min⁻¹) | 9.1 × 10⁻⁴ ± 2 × 10⁻⁵ | 7.8 × 10⁻⁴ ± 3 × 10⁻⁵ |
| **10** (min⁻¹) | 1.0 × 10⁻⁴ ± 5 × 10⁻⁵ | 5.7 × 10⁻⁴ ± 4 × 10⁻⁵ |

**"n.a." for k6 means methylglyoxal was not analysed in the pH-stat run — a genuine missing cell,
not a zero.**

Authors' conclusion, verbatim (§3.2, p. 441):
> "The estimated rate constants **did not change drastically** whether or not the system had a
> constant pH. This made us conclude that it was **justified – within limits –** to use non-constant
> pH studies to study pH effects on rate constants under the studied conditions."

**INTERNAL INCONSISTENCY, FLAGGED — the "did not change drastically" claim does not survive a
cell-by-cell check.** Ratios (pH constant / pH drop ⩽ 1), computed by me [Z]:

| k | ratio [Z] | verdict |
|---|---|---|
| 1 | 1.64× | outside both HPD intervals |
| 2 | 1.46× | outside both HPD intervals |
| 3 | 0.42× | 2.4-fold **decrease** |
| 4 | **2.00×** | doubled |
| 5 | **31.7×** | **thirty-two-fold** — and both intervals are tight (±12% and ±8%), so this is far outside error |
| 6 | — | not analysed |
| 7 | 1.54× | outside both HPD intervals |
| 8 | 0.67× | 1.5-fold decrease |
| 9 | 0.86× | ~within combined error |
| 10 | **5.70×** | nearly six-fold |

**Only k9 is genuinely unchanged. Three constants (k5, k10, k4) change by 2–32×, and k5's
32-fold shift is the largest single discrepancy in the paper.** The authors' own §3.4 later leans
on k5's behaviour ("3-DG also showed low pH-dependence in formic acid formation (step 5)") without
reconciling it with this. **Quote both:**
- §3.2, p. 441: "The estimated rate constants did not change drastically whether or not the system
  had a constant pH."
- §4, p. 447: "the observed **consistency for most of the estimated parameters** when the pH was
  kept constant gives an additional indication that the deduced pH-dependence is quite accurate."
- Table 1, p. 442: k5 changes by **31.7×**; k10 by **5.7×**; k4 by **2.0×**.

**This is the single most important honesty flag in this dossier.** The entire justification for
using drifting-pH data to build a pH-rate profile rests on this comparison, and the comparison is
weaker than stated. **Any hold-out built on Table 2's pH ladder inherits this.**

---

## 7. TABLE 4 — MOLAR-RATIO EXPERIMENTS (1:1 vs 2:1 vs 1:2)

**Anchor: Table 4, p. 445 (PDF p. 9).** Title as printed:
> "Estimated rate constants (k) ± 95% HPDᵃ interval at different initial concentrations of the
> reactants, expressed as **molar ratio glucose/glycine**"

Column headers as printed: `k` | `1:1` | `2:1` | `1:2`
Footnotes as printed: "Samples heated at **100 °C and reaction initial pH 6.8**." and "ᵃ Highest
posterior density."

**All values [F].** Read from the PDF p. 9 raster.

| k (unit as printed) | 1:1 | 2:1 | 1:2 |
|---|---|---|---|
| **1** (l mol⁻¹ min⁻¹) | 1.1 × 10⁻⁵ ± 4 × 10⁻⁶ | 1.5 × 10⁻⁵ ± 1 × 10⁻⁶ | 1.4 × 10⁻⁵ ± 8 × 10⁻⁷ |
| **2** (min⁻¹) | 8.9 × 10⁻⁴ ± 4 × 10⁻⁵ | 1.2 × 10⁻³ ± 2 × 10⁻⁴ | 1.1 × 10⁻³ ± 2 × 10⁻⁴ |
| **3** (min⁻¹) | 5.9 × 10⁻³ ± 8 × 10⁻⁴ | 5.1 × 10⁻³ ± 3 × 10⁻³ | 3.5 × 10⁻³ ± 2 × 10⁻³ |
| **4** (min⁻¹) | 6.5 × 10⁻³ ± 3 × 10⁻⁴ | 8.0 × 10⁻³ ± 9 × 10⁻⁴ | 8.5 × 10⁻³ ± 9 × 10⁻⁴ |
| **5** (min⁻¹) | 4.1 × 10⁻³ ± 5 × 10⁻⁴ | 3.1 × 10⁻³ ± 1 × 10⁻⁴ | 2.9 × 10⁻³ ± 2 × 10⁻⁴ |
| **6** (min⁻¹) | 5.2 × 10⁻³ ± 2 × 10⁻⁴ | 5.3 × 10⁻³ ± 6 × 10⁻⁴ | 8.5 × 10⁻³ ± 6 × 10⁻⁴ |
| **7** (min⁻¹) | 1.3 × 10⁻² ± 4 × 10⁻⁴ | 1.1 × 10⁻² ± 1 × 10⁻³ | 1.3 × 10⁻² ± 3 × 10⁻³ |
| **8** (min⁻¹) | 1.3 × 10⁺⁰ ± 2 × 10⁻¹ | 9.3 × 10⁻¹ ± 9 × 10⁻² | 9.6 × 10⁻¹ ± 9 × 10⁻² |
| **9** (min⁻¹) | 9.1 × 10⁻⁴ ± 2 × 10⁻⁵ | 1.4 × 10⁻³ ± 3 × 10⁻⁴ | 1.0 × 10⁻³ ± 1 × 10⁻⁴ |
| **10** (min⁻¹) | 1.0 × 10⁻⁴ ± 5 × 10⁻⁵ | 1.0 × 10⁻⁴ ± 6 × 10⁻⁵ | 2.1 × 10⁻⁴ ± 5 × 10⁻⁵ |

**The exact molar concentrations behind "2:1" and "1:2" are NEVER STATED.** §3.3 says only "the
concentrations of glucose and glycine were changed to give molar ratios of 2:1 and 1:2". Fig. 8's
panel (a) shows glucose starting near **200 mmol/l** and glycine near **100 mmol/l** for the 2:1
run [fig], implying glucose was held at 0.2 M and glycine halved rather than glucose doubled.
Fig. 9 (1:2) likewise shows glucose ~200 and the second species ~... **the marker assignment in
Fig. 9(a) is ambiguous at page-raster resolution: I could not confirm which of the two curves is
glycine. Recorded as unreadable.** **This is a real gap: the second-order constant k1 cannot be
independently checked without the actual concentrations.**

Authors' conclusion, verbatim (§3.3, p. 445):
> "As observed in Table 4, the **estimated parameters show no significant difference**. Some
> variation (as is to be expected because of experimental error) is observed, but on the whole it
> could be concluded that the **variation in the values remained within the 95% confidence
> intervals**."

**PARTIAL CONTRADICTION, FLAGGED.** Checking the claim cell by cell [Z] against the printed HPDs:
- **k1**: 1.1 ± 0.4 vs 1.5 ± 0.1 vs 1.4 ± 0.08 (×10⁻⁵). The 2:1 and 1:2 intervals
  (1.4–1.6 and 1.32–1.48) **do not overlap** the 1:1 point estimate's own tight range in the other
  two columns, though 1:1's wide ±4 × 10⁻⁶ does cover them. Marginal.
- **k3**: 5.9 ± 0.8 vs 5.1 ± 3.0 vs **3.5 ± 2.0** (×10⁻³). The 1:2 interval (1.5–5.5) **barely
  touches** the 1:1 interval (5.1–6.7). A 1.7-fold shift.
- **k6**: 5.2 ± 0.2 vs 5.3 ± 0.6 vs **8.5 ± 0.6** (×10⁻³). The 1:2 interval (7.9–9.1) and the 1:1
  interval (5.0–5.4) **do not overlap at all.** A clean **1.6-fold increase** in the retro-aldol
  step when glycine is doubled — which is chemically *interesting* (more amino acid, more
  amino-catalysed retro-aldol) but directly contradicts "no significant difference."
- **k4**: 6.5 ± 0.3 vs 8.0 ± 0.9 vs 8.5 ± 0.9 (×10⁻³). 1:1 (6.2–6.8) vs 1:2 (7.6–9.4):
  **no overlap.** A 1.3-fold increase.
- **k9**: 9.1 ± 0.2 vs **14 ± 3** vs 10 ± 1 (×10⁻⁴). 1:1 (8.9–9.3) vs 2:1 (11–17): **no overlap.**
- **k5**: 4.1 ± 0.5 vs 3.1 ± 0.1 vs 2.9 ± 0.2 (×10⁻³). 1:1 (3.6–4.6) vs 1:2 (2.7–3.1):
  **no overlap.**
So **six of the ten constants have non-overlapping 95% HPD intervals across the ratio arms.** Quote
both:
- §3.3, p. 445: "the estimated parameters show no significant difference … the variation in the
  values remained within the 95% confidence intervals."
- Table 4, p. 445: k4, k5, k6 and k9 have **non-overlapping** 95% HPD intervals between the 1:1 and
  at least one other arm; k6 differs by 1.6× with no overlap at all.
**The pattern is not random — k4, k6 and k9 (the three steps that consume or release glycine) all
move in the direction expected if glycine concentration matters.** The authors' framing ("the model
is robust to change in initial concentrations") is the conclusion they wanted; the table is more
equivocal than they say.

### 7.1 The three-way identity of the pH 6.8 / 1:1 column — DISJOINTNESS HAZARD

**Verified cell by cell:** Table 1's `pH drop ⩽ 1` column, Table 2's `pH 6.8` column and Table 4's
`1:1` column print **identical values and identical error bars for all ten constants**. They are
**one experiment printed three times in this single paper.**

**Consequence:** a wave that fits "Table 2 pH 6.8" and separately holds out "Table 4 1:1" (or
"Table 1 pH drop") would be holding out data it had already fit. **This is precisely the
Cerny-2007 rule-1 violation that Amendment 2 of the declaration had to fix retroactively.**
**Recommendation: declare this column exactly once, in exactly one place, and record the other two
appearances as aliases.**

---

## 8. FIGURES — inventory and digitised values

### 8.1 Fig. 1 — glucose degradation across the pH ladder
**Anchor: Fig. 1, p. 440 (PDF p. 4).** Caption verbatim: "Glucose degradation during heating of
glucose/glycine reaction mixture in phosphate buffer (0.1 M) at 100 °C and initial pH 4.8 (◆);
pH 5.5 (□); pH 6.0 (▲); pH 6.8 (✳) and pH 7.5 (●)."
Axes as printed: y = `mmol/l`, 0–250, ticks 50; x = `Time (min)`, 0–250, ticks 50.

Digitised from the text-layer axis scaffold plus the PDF p. 4 layout. **[fig]**, read precision
≈ ±8 mmol/l. **The five series overlap heavily in the upper half of the plot and I could not
separate pH 4.8 from pH 5.5 reliably below 150 min: those cells are marked unreadable.**

| time (min) | pH 4.8 | pH 5.5 | pH 6.0 | pH 6.8 | pH 7.5 |
|---|---|---|---|---|---|
| 0 | ~200 | ~200 | ~200 | ~200 | ~200 |
| 60 | unreadable | unreadable | ~185 | ~175 | ~155 |
| 120 | unreadable | ~190 | ~178 | ~160 | ~130 |
| 180 | ~190 | ~185 | ~170 | ~148 | ~112 |
| 240 | ~188 | ~180 | ~163 | ~138 | ~100 |

**Prefer the authors' own numbers [M] (§3.1, p. 440):**
> "As the initial pH increased, glucose degradation rate also increased (Fig. 1). **A maximum of
> 50% degradation was observed at pH 7.5** whereas, for the same heating time, **only 6% of glucose
> was degraded at pH 4.8**."

| pH | glucose degraded at end of run | provenance |
|---|---|---|
| 4.8 | **6%** | [M] |
| 7.5 | **50%** | [M] |
| 5.5, 6.0, 6.8 | **not stated numerically** | — |

### 8.2 Fig. 2 — mass balance (two panels only)
**Anchor: Fig. 2, p. 440 (PDF p. 4).** Caption verbatim: "Mass balance of reactants and reaction
products in heated glucose–glycine reaction mixture in phosphate buffer (0.1 M) at 100 °C and
**(a) pH 6.0, and (b) pH 7.5**. Glucose (Glu); fructose (Fru); Amadori compound (DFG); sum of
1-deoxy and 3-deoxyglucosone (DG); sum of formic and acetic acid (OA); methylglyoxal (MG);
melanoidins (Mel)."

**Only pH 6.0 and pH 7.5 are shown — no mass balance is published for pH 4.8, 5.5 or 6.8.**
Per-segment values are **unreadable** at page-raster resolution. The operative content is textual
[M] (§3.1, p. 440):
> "The results of the mass balance calculations showed that, **for pH ⩽ 6.0, the main quantified
> compound was the Amadori compound (DFG)** (Fig. 2(a)). The observed **gap (mass balance <100%)
> between 45 and 180 min** suggested that other products, not identified in the present study, were
> also formed. However, the fact that we calculated **100% recovery, eventually**, indicates that
> the acids formed are **stable end-products of scission reactions, leading to C1–C5 reaction
> products**. Moreover, under these reaction conditions **the predominance of 3-DG over 1-DG was
> also observed at pH ⩽ 6.0**."
> "At **pH 7.5, 20% of the nearly 50% of degraded glucose represented both formic and acetic acid,
> of which 14% was due to acetic acid** (Fig. 2(b))."

| quantity | value | conditions | provenance |
|---|---|---|---|
| organic acids as share of degraded glucose | **20%** | pH 7.5, 100 °C | [M] |
| acetic acid alone as share of degraded glucose | **14%** | pH 7.5, 100 °C | [M] |
| formic acid alone, by difference | **6%** | pH 7.5, 100 °C | **[Z]** |
| mass-balance gap window | 45–180 min | pH 6.0 | [M] |

**As in Part II, this is a MOLAR balance against initial reactant, not a carbon balance.** No
carbon-number weighting is applied anywhere.

### 8.3 Fig. 3 — organic acids at pH 5.5
**Anchor: Fig. 3, p. 441 (PDF p. 5).** Caption verbatim: "Organic acids formation during heating of
glucose/glycine reaction mixture in phosphate buffer (0.1 M) at 100 °C and pH 5.5. Formic acid (■);
acetic acid (□)."
Axes as printed: y = `Acids mmol/l`, 0.0–1.5, ticks 0.5; x = `Time (min)`, 0–250, ticks 50.
**Series values are unreadable at page-raster resolution** (the two curves are close together below
1.5 mmol/l). Operative content [M] (§3.1, p. 440): "at lower pH values, a **slightly higher
concentration of formic acid than acetic acid** was observed at the beginning of the heating
period (Fig. 3). The opposite was observed when the pH was increased."

### 8.4 Fig. 4 — pH drop across the ladder (a CONDITION, not an outcome)
**Anchor: Fig. 4, p. 441 (PDF p. 5).** Caption verbatim: "pH drop during heating of glucose/glycine
reaction mixture in phosphate buffer (0.1 M) at 100 °C and pH 4.8 (◆); pH 5.5 (□); pH 6.0 (▲);
pH 6.8 (✳) and pH 7.5 (●)."
Axes as printed: y = `pH`, **4 to 8, ticks every 1**; x = `Time (min)`, 0–250, ticks 50.

Digitised from the PDF p. 5 layout. **[fig]**, read precision ≈ ±0.15 pH unit (the y-axis is
compressed over 4 units, so this figure digitises **poorly** — treat these as indicative only).

| time (min) | pH 4.8 | pH 5.5 | pH 6.0 | pH 6.8 | pH 7.5 |
|---|---|---|---|---|---|
| 0 | ~4.8 | ~5.5 | ~6.0 | ~6.8 | ~7.5 |
| 60 | ~4.75 | ~5.4 | ~5.8 | ~6.4 | ~6.7 |
| 120 | ~4.7 | ~5.35 | ~5.7 | ~6.2 | ~6.4 |
| 180 | ~4.7 | ~5.3 | ~5.6 | ~6.0 | ~6.2 |
| 240 | ~4.65 | ~5.25 | ~5.55 | ~5.9 | ~6.1 |

Authors' statement [M] (§3.1, p. 440):
> "The importance of carboxylic acid formation is that it is related to the observed pH drop, **even
> though a buffer system (0.1 M) was used** (Fig. 4). **This phenomenon was particularly important at
> higher pH values.**"

**Key consequence, and the reason the Table 2 footnote exists:** the pH 7.5 run loses ~1.3 pH units
over 240 min, so the 1-unit data filter cuts it earliest. **The "pH 7.5" column of Table 2 is
fitted on the SHORTEST time window of the five, and its effective mean pH is nearer ~7.0.** This
is a structural asymmetry in the pH ladder that the paper does not quantify. **Flagged.**

### 8.5 Fig. 5 — the pH-stat run (10 responses, 5 panels)
**Anchor: Fig. 5, p. 442 (PDF p. 6).** Caption verbatim: "Glucose–glycine Maillard reaction at
100 °C with pH kept constant by addition of 1 N NaOH (initial pH set at 6.8 at room temperature).
Model fit (lines, calculated according to Scheme 1) to experimental data (markers). (a) Glucose (◇),
glycine (■); (b) fructose (▲), N-(1-deoxy-D-fructos-1-yl)–glycine (○); (c) 1-deoxyglucosone (●),
3-deoxyglucosone (△); (d) formic acid (◆), acetic acid (□); (e) melanoidins (✳)."

Panel inventory and terminal values, read from the PDF p. 6 raster. **All [fig]**, read precision
as noted per panel.

| panel | y-axis label & range | x range | terminal values at ~120 min [fig] |
|---|---|---|---|
| (a) | `mmol/l`, 0–250, ticks 50 | 0–150 | glycine (■) ~155; glucose (◇) ~100 (±8 mmol/l) |
| (b) | `mmol/l`, 0–25, ticks 5 | 0–150 | fructose (▲) ~20 rising; DFG (○) peak ~13 at ~50 min then ~10 (±1 mmol/l) |
| (c) | `DG mmol/l`, 0–0.7, ticks 0.1 | 0–150 | 3-DG (△) plateau ~0.60; 1-DG (●) plateau ~0.30 (±0.03 mmol/l) |
| (d) | `Acids mmol/l`, 0–32, ticks 4 | 0–150 | acetic (□) ~28 rising steeply; formic (◆) ~8 (±1 mmol/l) |
| (e) | `Mel mmol/l`, 0–10, ticks 2 | 0–150 | melanoidins (✳) ~8.5, essentially linear (±0.4 mmol/l) |

**Note the striking result visible in panel (c): under pH-stat control, 3-DG (~0.60 mmol/l) is
2× 1-DG (~0.30 mmol/l) and both PLATEAU rather than peaking-and-decaying.** Compare the same
species in the drifting-pH 2:1 run (Fig. 8c), where the ratio is ~5:1. **This is the visual
counterpart of the 31.7× shift in k5 between Tables 1's two columns (§6).**

### 8.6 Fig. 7 — the simultaneous fit across all five pHs (10 panels)
**Anchor: Fig. 7, p. 444 (PDF p. 8).** Caption verbatim: "Model fit (lines, calculated according to
Scheme 1) to experimental data (markers) of glucose/glycine aqueous system heated at 100 °C and
initial pH 4.8 (◆); pH 5.5 (✳); pH 6.0 (○); pH 6.8 (+) and pH 7.5 (▲). (a) Glucose (Glu);
(b) glycine (Gly); (c) N-(1-deoxy-D-fructos-1-yl)–glycine (DFG); (d) fructose (Fru);
(e) **3-deoxyglycosone (1-DG)**; (f) 1-deoxyglucosone (1-DG). (g) formic acid (FA); (h) acetic acid
(AA); (i) methylglyoxal (MG); (j) melanoidins (Mel). The initial pH values are reported as set at
room temperature."

**PRINTING ERROR IN THE CAPTION, FLAGGED:** panel (e) is labelled "**3-deoxyglycosone (1-DG)**" —
both the species name is misspelled and the abbreviation is wrong. From the axis label on the page
raster, panel (e)'s y-axis reads `3-DG mmol/l` and panel (f)'s reads `1-DG mmol/l`, so **(e) is
3-deoxyglucosone (3-DG) and (f) is 1-deoxyglucosone (1-DG)**. The caption's "(1-DG)" for panel (e)
is a typo.

Panel inventory with axis ranges, read from the PDF p. 8 raster. **This is the paper's central data
figure — all five pH columns of Table 2 were fitted to these ten responses simultaneously.**

| panel | y-axis label & range | x range | note |
|---|---|---|---|
| (a) | `Glu mmol/l`, 100–220, ticks 20 (**suppressed zero**) | 0–250 | see §8.1 for values |
| (b) | `Gly mmol/l`, 100–220, ticks 20 (**suppressed zero**) | 0–250 | |
| (c) | `DFG mmol/l`, 0–14, ticks 2 | 0–250 | |
| (d) | `Fru mmol/l`, 0–20, ticks 5 | 0–250 | |
| (e) | `3-DG mmol/l`, 0–1.4, ticks 0.2 | 0–250 | caption says "3-deoxyglycosone (1-DG)" — typo |
| (f) | `1-DG mmol/l`, 0–0.5, ticks 0.1 | 0–250 | |
| (g) | `FA mmol/l`, 0–10, ticks 2 | 0–250 | |
| (h) | `AA mmol/l`, 0–30, ticks 5 | 0–250 | |
| (i) | `MG mmol/l`, 0–14, ticks 2 | 0–250 | |
| (j) | `Mel mmol/l`, 0–12, ticks 2 | 0–250 | |

**Per-point digitisation of Fig. 7 was NOT attempted.** Ten panels × five overlapping pH series
with markers of similar size makes reliable per-series separation impossible at page-raster
resolution, and a wrong series assignment here would be worse than no number. **The axis ranges
above are the honest deliverable; the fitted summary of these data is Table 2, which is
transcribed completely in §4.** Panels (a) and (b) additionally use **suppressed zeros** (y starts
at 100), which exaggerates the apparent degradation — worth knowing before reading the figure by eye.

### 8.7 Figs. 8 and 9 — the ratio experiments
**Anchor: Fig. 8, p. 445 (PDF p. 9); Fig. 9, p. 446 (PDF p. 10).** Captions verbatim:
- Fig. 8: "**Glucose/glycine molar ratio 2:1.** Model fit (lines, calculated according to Scheme 1)
  to experimental data (markers) for the glucose/glycine aqueous system heated at 100 °C and
  pH 6.8. (a) Glucose (◇), glycine (■); (b) fructose (▲), N-(1-deoxy-D-fructos-1-yl)–glycine (○),
  methylglyoxal (+); (c) 1-deoxyglucosone (●), 3-deoxyglucosone (△); (d) formic acid (◆), acetic
  acid (□); (e) melanoidins (✳)."
- Fig. 9: "**Glucose/glycine molar ratio 1:2.**" — otherwise identical structure and legend.

Fig. 8 axis ranges and terminal values, from the PDF p. 9 raster. **All [fig]**:

| panel | y-axis & range | x range | terminal values at ~180 min [fig] |
|---|---|---|---|
| (a) | `mmol/l`, 0–250, ticks 50 | 0–200 | glucose (◇) ~140 from ~200; glycine (■) ~62 from ~100 |
| (b) | `mmol/l`, 0–25, ticks 5 | 0–200 | fructose (▲) ~21; DFG (○) ~7; MG (+) ~4 |
| (c) | `DG mmol/l`, 0–0.6, ticks 0.2 | 0–200 | 3-DG (△) plateau ~0.55; 1-DG (●) ~0.10 |
| (d) | `Acids mmol/l`, 0–16, ticks 4 | 0–200 | acetic (□) ~14; formic (◆) ~3 |
| (e) | `Mel mmol/l`, 0–10, ticks 2 | 0–200 | melanoidins (✳) ~7.5, near-linear |

**Fig. 8(a) is the only direct evidence of the actual concentrations behind "2:1": glucose starts
at ~200 mmol/l and glycine at ~100 mmol/l**, i.e. **glycine was halved, glucose held at 0.2 M**.
[fig], ±8 mmol/l. **Fig. 9(a)'s two curves could not be separated reliably: the actual 1:2
concentrations are unreadable from the figure and are not stated in the text.**

**Note that the 2:1 3-DG:1-DG ratio here is ~5.5:1, against ~2:1 in the pH-stat run (Fig. 5c) and
Table 1's 31.7× k5 discrepancy. The deoxyosone partition is the least stable thing in this paper.**

### 8.8 Fig. 10 — sensitivity analysis
**Anchor: Fig. 10, p. 446 (PDF p. 10).** Caption verbatim: "Sensitivity analysis of the responses
depicted in the proposed kinetic model (Scheme 1) for rate constants **k₁, k₃, k₆ and k₈**. Within
each graph, **the sensitivity was zero for the responses not shown**. Glucose (Glu); glycine (Gly);
fructose (Fru); N-(1-deoxy-D-fructos-1-yl)–glycine (DFG); acetic acid (AA); formic acid (FA);
1-deoxyglucosone (1-DG); 3-deoxyglucosone (3-DG) and melanoidins (Mel)."

Four panels, y = ∂Y/∂kᵢ, x = Time (min) 0–120. Axis ranges as printed:

| panel | y-axis label | y range as printed | responses shown |
|---|---|---|---|
| (a) | `dY/dk1` | −3.E+03 … 2.E+06 (mixed labelling: −3.E+06 to 2.E+06 on the ticks) | DFG (large +), Fru (+), Gly (large −), Glu (large −) |
| (b) | `dY/dk6` | −1000 … 1500 | AA (large +), MG (+), Gly (+), FA/Mel (≈0), DFG (−) |
| (c) | `dY/dk3` | −2.E+03 … 2.E+03 | Glu (+), Fru (−) |
| (d) | `dY/dk8` | −0.4 … 0.4 | AA (+), 1-DG (−) |

**NO NUMERIC SENSITIVITY VALUES ARE PRINTED.** The curves are unlabelled beyond species names, and
the y-axis tick labelling in panel (a) is internally inconsistent between the text layer and the
raster. **Per-curve digitisation would be unreliable; recorded as unreadable.** The operative
content is the authors' qualitative conclusions [M] (§3.4, pp. 446–447):
> "As expected, **k₁ has a strong influence in the main products** formed throughout the Maillard
> reaction. Moreover, together with k₄ and k₇ (DFG enolisation step), **the DFG retro-aldolisation
> step (k₆) appeared not to be a redundant parameter.** Besides methylglyoxal formation, it has a
> **positive influence in glycine regeneration**. Also, the isomerisation step of fructose into
> glucose (k₃) becomes more evident for longer heating periods, which may explain its dependence on
> the pH drop. In acetic acid formation, besides k₁, **k₇ (1-DG formation from DFG) also shows a
> strong positive influence (results not shown)**, which explains the **low sensitivity of acetic
> acid to k₈**. This result is well in line with the observed low pH-dependence of k₈. **The
> limiting factor on the amount of acetic acid seems to be the formation of 1-DG.**"

Note the "(results not shown)" for the k₇ sensitivity — **another unpublished result**.
Conclusion, verbatim (Abstract, p. 437): "**Again, the model performed well; all steps were
important and the model was consistent with the established reaction mechanism.**"

---

## NEW-PARAMETER TABLE (consolidated)

**This paper is SINGLE-TEMPERATURE (100 °C). It contains NO activation energies and NO reverse rate
constants.**

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| **k1** (Glu + Gly → DFG) | 7.9×10⁻⁷ ± 2×10⁻⁸ / 2.3×10⁻⁶ ± 5×10⁻⁸ / 4.9×10⁻⁶ ± 1×10⁻⁷ / 1.1×10⁻⁵ ± 4×10⁻⁶ / 1.8×10⁻⁵ ± 8×10⁻⁷ | **l mol⁻¹ min⁻¹** | 100 °C, 0.2 M equimolar, 0.1 M phosphate; pH 4.8 / 5.5 / 6.0 / 6.8 / 7.5 | Table 2, p. 442 | [F] |
| **k2** (Glu → Fru) | 8.8×10⁻⁶ ± 1×10⁻⁶ / 5.5×10⁻⁵ ± 4×10⁻⁶ / 2.3×10⁻⁴ ± **3×10⁻³** / 8.9×10⁻⁴ ± 4×10⁻⁵ / 2.1×10⁻³ ± 7×10⁻⁵ | min⁻¹ | as above | Table 2, p. 442 | [F] — pH 6.0 error defective |
| **k3** (Fru → Glu) | **0.00** / **0.00** / 2.8×10⁻³ ± 2×10⁻³ / 5.9×10⁻³ ± 8×10⁻⁴ / 1.4×10⁻² ± 9×10⁻⁴ | min⁻¹ | as above | Table 2, p. 442 | [F] — two deleted-zero cells |
| **k4** (DFG → 3-DG + Gly) | 1.5×10⁻³ ± 7×10⁻⁵ / 3.1×10⁻³ ± 1×10⁻⁴ / 4.5×10⁻³ ± 3×10⁻⁴ / 6.5×10⁻³ ± 3×10⁻⁴ / 1.5×10⁻² ± 1×10⁻³ | min⁻¹ | as above | Table 2, p. 442 | [F] — **the only fully monotone row** |
| **k5** (3-DG → FA) | **0.00** / 6.3×10⁻³ ± 2×10⁻³ / 2.0×10⁻² ± 2×10⁻³ / 4.1×10⁻³ ± 5×10⁻⁴ / 1.9×10⁻² ± 1×10⁻² | min⁻¹ | as above | Table 2, p. 442 | [F] — **badly non-monotone** |
| **k6** (DFG → MG) | 5.0×10⁻⁴ ± 8×10⁻⁵ / 9.1×10⁻⁴ ± 1×10⁻⁴ / 2.7×10⁻³ ± 1×10⁻⁴ / 5.2×10⁻³ ± 2×10⁻⁴ / 1.1×10⁻² ± 9×10⁻⁴ | min⁻¹ | as above | Table 2, p. 442 | [F] — **monotone** |
| **k7** (DFG → 1-DG + Gly) | 7.7×10⁻⁴ ± 7×10⁻⁵ / 4.0×10⁻⁴ ± 1×10⁻⁴ / 3.1×10⁻³ ± 3×10⁻⁴ / 1.3×10⁻² ± 4×10⁻⁴ / 1.5×10⁻² ± 1×10⁻³ | min⁻¹ | as above | Table 2, p. 442 | [F] — pH 5.5 dips below pH 4.8 |
| **k8** (1-DG → AA) | 9.9×10⁻² ± 9×10⁻³ / 6.2×10⁻² ± 2×10⁻² / 1.0×10⁻¹ ± 1×10⁻² / **1.3×10⁺⁰ ± 2×10⁻¹** / 4.5×10⁻¹ ± 4×10⁻² | min⁻¹ | as above | Table 2, p. 442 | [F] — non-monotone; pH 6.8 outlier |
| **k9** (3-DG + Gly → Mel) | 5.0×10⁻⁵ ± 4×10⁻⁶ / 2.1×10⁻⁴ ± 1×10⁻⁵ / 2.4×10⁻⁴ ± 2×10⁻⁵ / 9.1×10⁻⁴ ± 2×10⁻⁵ / 5.3×10⁻⁴ ± 5×10⁻⁵ | min⁻¹ **(unit WRONG — bimolecular)** | as above | Table 2, p. 442 | [F] — decreases at pH 7.5 |
| **k10** (Glu → FA + AA) | **0.00** / 2.1×10⁻⁵ ± 3×10⁻⁶ / 6.3×10⁻⁵ ± **6×10⁻⁴** / 1.0×10⁻⁴ ± 5×10⁻⁵ / 3.1×10⁻⁴ ± 6×10⁻⁵ | min⁻¹ | as above | Table 2, p. 442 | [F] — one deleted zero; pH 6.0 error defective |
| **pD, k1** | **0.50 ± 0.14** | dimensionless (decades per pH unit) | fitted across pH 4.8–7.5, 100 °C | **Table 3, p. 443** | **[F]** |
| **pD, k2** | **0.88 ± 0.27** | dimensionless | as above | Table 3, p. 443 | [F] — closest to pure base catalysis |
| **pD, k3** | **0.46 ± 0.46** | dimensionless | as above | Table 3, p. 443 | [F] — not significant |
| **pD, k4** | **0.35 ± 0.10** | dimensionless | as above | Table 3, p. 443 | [F] |
| **pD, k5** | **0.09 ± 1.10** | dimensionless | as above | Table 3, p. 443 | [F] — **indeterminate** |
| **pD, k6** | **0.51 ± 0.14** | dimensionless | as above | Table 3, p. 443 | [F] |
| **pD, k7** | **0.60 ± 0.52** | dimensionless | as above | Table 3, p. 443 | [F] — marginal |
| **pD, k8** | **0.40 ± 0.60** | dimensionless | as above | Table 3, p. 443 | [F] — not significant |
| **pD, k9** | **0.40 ± 0.38** | dimensionless | as above | Table 3, p. 443 | [F] — marginal |
| **pD, k10** | **0.54 ± 0.36** | dimensionless | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k1** | **−8.43 ± 0.91** | log₁₀ of a value in **l mol⁻¹ min⁻¹** (**unit not printed**) | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k2** | **−9.14 ± 1.68** | log₁₀ of min⁻¹ (**unit not printed**) | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k3** | **−5.36 ± 3.26** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k4** | **−4.46 ± 0.60** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k5** | **−2.59 ± 7.54** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] — **spans 15 orders of magnitude** |
| **Log(kₑ), k6** | **−5.75 ± 0.85** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k7** | **−6.24 ± 3.25** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k8** | **−3.12 ± 2.40** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k9** | **−6.03 ± 2.35** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **Log(kₑ), k10** | **−7.56 ± 2.35** | log₁₀ of min⁻¹ | as above | Table 3, p. 443 | [F] |
| **pH-dependence functional form** | **k_obs = kₑ(10^{pH · pD})**, equivalently **log(k_obs) = log kₑ + pD × pH** | — | 100 °C, pH 4.8–7.5 | §3.2, pp. 442–443 | [F] — **abstract prints the exponent as `pD · pH`; identical, see §5.1** |
| pD reference values | **−1 = pure specific acid catalysis; +1 = pure specific base catalysis; 0 = no pH dependence** | dimensionless | — | §3.2, p. 443 | [C] (Loudon 1991) |
| Rate constants, pH-stat vs drifting pH | full 10-row comparison | as printed | 100 °C, pH 6.8 | Table 1, p. 442 | [F] — **k5 differs 31.7×, k10 5.7×; see §6** |
| Rate constants, molar ratios | full 10-row comparison | as printed | 100 °C, pH 6.8; 1:1, 2:1, 1:2 | Table 4, p. 445 | [F] — **6/10 rows have non-overlapping HPDs; see §7** |
| glucose degraded | **6** | % of initial | pH 4.8, 100 °C, end of run | §3.1, p. 440 | [M] |
| glucose degraded | **50** | % of initial | pH 7.5, 100 °C, end of run | §3.1, p. 440 | [M] |
| organic acids | **20** | % of degraded glucose | pH 7.5, 100 °C | §3.1, p. 440 | [M] |
| acetic acid | **14** | % of degraded glucose | pH 7.5, 100 °C | §3.1, p. 440 | [M] |
| formic acid | **6** | % of degraded glucose | pH 7.5, 100 °C | Fig. 2b / §3.1 | **[Z]** (20 − 14) |
| melanoidin ε at 470 nm | 0.64 ± 0.03 | l mmol⁻¹ cm⁻¹ | glucose/glycine | §2.2, p. 439 | **[C]** — primary is `martins2003c.pdf` |
| HMF | "**one order of magnitude lower … µmolar range**"; increases with **decreasing** pH | µmol/l | all pHs | §3.1, p. 440 | [M] — **no numeric value published** |
| literature pH slopes | 0.38 / 0.13 / 0.27 / 0.64 / 0.37 / ~0.1 / 0.1 | dimensionless | various Maillard systems | §3.2, pp. 443–444 | **[C]** and **[C]-derived** |
| Sensitivity analysis | qualitative only; **no numeric values printed** | — | k1, k3, k6, k8 | Fig. 10, p. 446 | [M] — **unreadable as numbers** |
| goodness-of-fit statistic | **NONE PRINTED** | — | — | whole paper | negative finding |
| **activation energies** | **NONE — single temperature (100 °C)** | — | — | whole paper | negative finding, high confidence |
| **reverse rate constants** | **NONE (k3 is the reverse isomerisation leg only, not an Amadori reversal)** | — | — | Scheme 1, p. 439 | negative finding, high confidence |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Blocking issue first: this paper is, on my reading, ALREADY DECLARED as "Martins 2003 thesis
Ch. 6 pH ladder", with role "HOLD-OUT except pH 5.5/6.8 columns (which are the already-FIT
conditions)".** The declared pH list (4.8/5.5/6.0/6.8/7.5) and the declared "fitted pH exponents"
are verbatim matches to Tables 2 and 3 here. **Do not create a second declaration row.**
**Recommendation: amend the existing Ch. 6 row to (i) name the journal anchor
`data/articles/martins2005b.pdf`, Tables 2 and 3, pp. 442–443; (ii) record the three-way column
identity of §7.1; (iii) record that the pH 4.8 column has three information-free cells; (iv) record
the Table 1 pH-stat discrepancy of §6, which weakens the whole pH-ladder's foundation.**

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 2, pH 6.8 column, all 10 constants** | **FIT** (matches the declared carve-out) | pH | The trunk's own pH. **Declare this column ONCE. It is printed three times in this paper** (Table 1 `pH drop ⩽ 1`, Table 2 `pH 6.8`, Table 4 `1:1`) with identical values — record the other two as aliases so no later wave can hold out data already fit. |
| **Table 2, pH 5.5 column, all 10 constants** | **FIT** (matches the declared carve-out) | pH | pH 5.5 is a Martins-standard condition already used elsewhere in the trunk. Note k3 = 0.00 here is a deleted parameter. |
| **Table 2, pH 4.8, 6.0 and 7.5 columns** | **HOLD-OUT** | pH | The genuine extrapolation arm, exactly as the declaration reasons. **Exclusions within the arm:** (a) k3, k5, k10 at pH 4.8 — deleted zeros, no information to score; (b) k2 and k10 at pH 6.0 — printed error bars exceed their own estimates (probable typos); (c) k5 everywhere — the row is non-monotone by 4.9× and cannot be scored meaningfully. **That leaves ~21 of 30 scorable cells, which is still a strong test.** |
| **Table 3, pD and Log(kₑ) for k1, k2, k4, k6** | **FIT as a functional-form prior, NOT as a scored target** | — | These four are the only exponents convincingly different from zero. They give the trunk a defensible pH-scaling law for the initiation step, sugar isomerisation, 1,2-enolisation and retro-aldolisation. **They are derived FROM Table 2 and are not independent of it — so they must never be scored against a model that was fit to Table 2.** |
| **Table 3, pD for k3, k5, k7, k8, k9, k10** | **diagnostic_only** | — | Not significantly different from zero at 95% (k5 catastrophically so: 0.09 ± 1.10, with log kₑ spanning 15 orders of magnitude). Recording them as "pH-insensitive within error" is the honest use; treating them as measured exponents is not. |
| **Table 1 (pH-stat vs drifting pH)** | **HOLD-OUT — and a high-value one** | pH control | This is the cleanest *methodological* hold-out in the whole Martins corpus: same nominal condition, same lab, one variable changed (pH held vs. allowed to drift). If the trunk models in-run pH at all, it should predict which constants move and which don't. **It currently would fail: the authors say "no drastic change" but k5 moves 31.7×, k10 5.7×, k4 2.0× (§6).** Scoring this honestly is more informative than scoring another concentration curve. |
| **Table 4, 2:1 and 1:2 columns** | **HOLD-OUT** | reactant ratio | Orthogonal to the pH axis and to everything else in the trunk. A rate constant that is genuinely a rate constant must not move with initial concentration. **The paper claims it doesn't; six of ten rows have non-overlapping 95% HPDs (§7), and the three glycine-coupled steps (k4, k6, k9) move in the chemically expected direction.** This is a real, sharp test of whether the trunk's rate law is correctly specified. **Caveat: the actual 2:1 and 1:2 molar concentrations are never stated (only inferable from Fig. 8a as ~200/~100 mmol/l); the 1:2 concentrations are unreadable. A scoring wave must resolve this or the second-order step k1 cannot be checked.** |
| **Table 4, 1:1 column** | **excluded — alias of the FIT pH 6.8 column** | — | Identical values. See §7.1. |
| **§8.1 glucose degradation (6% at pH 4.8, 50% at pH 7.5)** | **HOLD-OUT** | pH | Model-free, cheap, and the two extreme pHs are exactly the extrapolation arm. Two numbers, but they bracket an 8-fold span in reactant conversion. |
| **§8.2 acid partition at pH 7.5 (20% acids, 14% acetic)** | **HOLD-OUT** | pH | Model-free selectivity constraint at the far end of the pH arm. Tests the FA/AA branching ratio, which is the k5-vs-k8 competition — precisely the pair the model handles worst. |
| **Fig. 7 (the 10-response × 5-pH fit)** | **diagnostic_only unless someone re-digitises it properly** | pH | Five overlapping series per panel; I judged per-point separation unreliable at page-raster resolution and did not fake it. If the trunk needs the raw time-courses rather than the fitted constants, this figure must be re-digitised from a higher-resolution source. **Note panels (a) and (b) use suppressed zeros.** |
| **Fig. 10 sensitivity analysis** | **not usable — no numbers printed** | — | Qualitative curves only, plus a "(results not shown)" for the k7 sensitivity. |
| **HMF** | **not usable — no numbers published** | — | Measured, judged negligible ("µmolar range"), and discarded without publication. Only the *direction* (increases with decreasing pH) is quotable. |

**Cut-axis summary.** This paper supports **two orthogonal cuts**, which is rare and valuable:
1. **pH** — FIT at 5.5 and 6.8, HOLD-OUT at 4.8, 6.0 and 7.5. This is the declaration's existing
   convention and should stand.
2. **Reactant ratio** — FIT at 1:1, HOLD-OUT at 2:1 and 1:2. **This axis is currently undeclared
   and is the cheapest genuine extrapolation test available anywhere in the Martins corpus**,
   because a correctly specified rate law *must* be ratio-invariant. Recommend adding it.
A third, methodological cut (**pH-stat vs. drifting**, Table 1) is available and, in my judgement,
the most diagnostic of the three.

**Circularity risks, stated plainly:**
1. **Highest — the three-way column identity (§7.1).** Table 1's `pH drop ⩽ 1`, Table 2's `pH 6.8`
   and Table 4's `1:1` are one experiment printed three times. Declaring them separately would
   both inflate the evidence base and create a fit/hold-out overlap.
2. **Structural — Table 3 is derived from Table 2.** The pH exponents are a regression on the pH
   ladder. Fitting the ladder and then "validating" against the exponents (or vice versa) is
   circular. Pick one.
3. **Cross-paper — the model itself is imported.** Scheme 1 is `martins2005.pdf` Scheme 4,
   unchanged. The **topology** was selected there (by the ΔAICc = 0 / 276 / 287.46 discrimination),
   and the **Ea's** are there too. If a wave fits this paper's constants *and* that paper's Ea's,
   it is fitting one model's parameters from two papers — legitimate, but it must not then be
   described as two independent validations.
4. **Foundational — the pH-drift caveat (§6 and §8.4).** Every column of Table 2 is fitted to a
   *drifting* pH trimmed at 1 unit, and the higher-pH runs drift fastest, so the "pH 7.5" column is
   fitted on the shortest window and its effective mean pH is nearer ~7.0. The one control that
   was meant to justify this (Table 1) shows a 31.7× shift in k5. **A hold-out that fails at
   pH 7.5 may be failing on the pH trajectory rather than on the chemistry, and the report must
   say so rather than scoring it as a chemistry miss.**
