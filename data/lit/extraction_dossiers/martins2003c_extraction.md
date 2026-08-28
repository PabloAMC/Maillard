# Martins & van Boekel 2003 (10.1016/S0308-8146(03)00219-X) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/martins2003c.pdf` (8 pp.). Read method: both — `pdftotext -layout`
text layer for the body and Tables 2–4, plus the page raster of PDF p. 5 (journal p. 139) read with
the Read tool to recover **Table 1's minus signs**, which the text layer drops entirely (see the
DEFECT note in §3).

---

## 0. PAPER IDENTITY — **DOES NOT MATCH THE ORCHESTRATOR'S EXPECTED IDENTITY**

> **IDENTITY CORRECTION, STATED IN BOLD AT THE TOP AS REQUIRED.**
> **`data/articles/martins2003c.pdf` IS NOT PART I AND IS NOT PART II of the Amadori-degradation
> pair.** It is a **different paper in a different journal**: Martins & van Boekel, "Melanoidins
> extinction coefficient in the glucose/glycine Maillard reaction", **Food Chemistry 83 (2003)
> 135–142**. It contains **no rate constants, no reaction scheme, and no kinetic model.** It is an
> analytical-methods paper about ¹⁴C-radiolabelling and extinction coefficients.
>
> **The wave brief expected `martins2003c.pdf` to be Part II. That is wrong. Part II is
> `martins2003.pdf`** (Carbohydr. Res. 338 (2003) 1665–1678), which is where the rate-constant
> table lives. Part I is `martins2003b.pdf` (Carbohydr. Res. 338 (2003) 1651–1663).
>
> A second, independent hazard: **the authors' own lettered self-citations are permuted relative to
> our repo keys, and permuted differently in different papers.** In Martins & Van Boekel (2005)
> Food Chem. **92** this paper is cited as "**2003b**"; in Martins & Van Boekel (2005) Food Chem.
> **90** the same paper is "**2003c**"; and in Part I (`martins2003b.pdf`) it is simply "ref. 17".
> **Never resolve a Martins 2003 citation by its letter suffix. Resolve it by journal + volume +
> page.**

| field | value |
|---|---|
| authors | Sara I.F.S. Martins, Martinus A.J.S. van Boekel* |
| title | "Melanoidins extinction coefficient in the glucose/glycine Maillard reaction" |
| section | "Analytical, Nutritional and Clinical Methods" |
| venue | Food Chemistry |
| volume/pages/year | 83 (2003) 135–142 |
| DOI | 10.1016/S0308-8146(03)00219-X |
| received / revised / accepted | Received 2 January 2003; received in revised form 28 April 2003; accepted 28 April 2003 |
| affiliation | Product Design and Quality Management Group, Dept. of Agrotechnology and Food Science, Wageningen University, PO Box 8129, 6700 EV Wageningen, The Netherlands |
| PDF character | born-digital (Elsevier / Acrobat Distiller 4.05 for Windows); text layer present but **systematically drops the minus sign and the "±" glyph** — Table 1 is unusable from the text layer alone and was re-read from the page raster |
| corresponding author | tiny.vanboekel@wur.nl |
| funding / collaboration | European Commission Short-Term Scientific Mission (STSM) of **COST Action 919, "Melanoidins in Food and Health"**; Portuguese FCT. Work partly done at the Procter Department, University of Leeds, UK, with **Prof. B.L. Wedzicha**. Microanalysis by Hennie Halm, Wageningen. |

### PRE-EXISTING-COVERAGE FLAG — DOUBLE-COUNTING HAZARD

`docs/reference/FIT_HOLDOUT_DECLARATION.md` (read-only; Amendment 1, 2026-08-28, lines ~209–215)
declares roles for "Martins 2003 Wageningen thesis (edepot.wur.nl/121418)" Tables 4.2.3, 4.1.1,
3.3.1 and Ch. 6. **These four journal papers are the thesis chapters.**

**My mapping proposal for this paper (proposal, NOT a verified fact — I have not opened the thesis
PDF):**

| declared thesis item | my proposed mapping | confidence | verdict |
|---|---|---|---|
| Table **3.3.1** "melanoidin ε (**1.00 ± 0.03 at 420 nm, 0.65 ± 0.02 at 470 nm**; T- and DP-invariant)", role **FIT (measured input, not a scored target)** | **This paper's Table 1, the `120 °C, pH 6.8` block** (p. 139). The declaration's two quoted values are **verbatim cell matches**: ε₄₂₀ = 1.00 ± 0.03 and ε₄₇₀ = 0.65 ± 0.02 are exactly the 120 °C / pH 6.8 row. The "T- and DP-invariant" descriptor matches this paper's §3.1 conclusion. Thesis chapter numbering "3.x" also matches: this is the melanoidin/colour chapter, which precedes the DFG chapters (4.x) and the pH chapter (6). | **VERY HIGH** | **DIRECT DOUBLE-COUNT. Same measurements, re-published.** If Table 3.3.1 is already FIT, this paper's Table 1 is already FIT and must not be re-declared as new. |
| Table 4.2.3 "reverse rates (Amadori → parent sugar)" | not in this paper — **there are no rate constants here at all** | high | n/a |
| Table 4.1.1 glycine-release yields | maps to `martins2003b.pdf` (Part I) Table 1 | high | n/a |
| Ch. 6 pH ladder | maps to `martins2005b.pdf` Tables 2 and 3 | high | n/a |

**Two precision notes on the existing declaration row, which I recommend amending:**
1. The declaration quotes only the **120 °C, pH 6.8** values but describes ε as "T- and
   DP-invariant." This paper actually measures **three** condition sets and reports ε₄₇₀ = 0.65,
   0.66 and **0.62** — invariance is the authors' *conclusion from a non-significance test*, not a
   measured identity. The pH 5.5 block also carries a **2.5× larger error bar** (±0.05 vs ±0.02).
   The declaration should carry all three, not one.
2. The value actually **used downstream** by the Martins kinetic papers is neither of these: it is
   the **three-system average, ε₄₇₀ = 0.64 ± 0.03 l mmol⁻¹ cm⁻¹** (§3.1, p. 138). Part I
   (`martins2003b.pdf` §2.10) and the 2005 Food Chem. 92 paper (`martins2005b.pdf` §2.2) both cite
   0.64 ± 0.03 (the 2005 paper prints it as "0.64 ± 0.03 l mmol⁻¹ cm⁻¹"). **If the trunk's
   melanoidin observable model is meant to reproduce the Martins papers, 0.64 ± 0.03 is the number
   to carry, not 0.65 ± 0.02.**

---

## 1. ONE-PARAGRAPH VERDICT

This paper gives the model exactly one thing, but gives it well: a **directly measured, radiotracer-
calibrated extinction coefficient** that converts a browning absorbance into a molar concentration
of **glucose molecules incorporated into non-dialysable melanoidins** — at four wavelengths
(420, 450, 470, 490 nm) × three condition sets (120 °C/pH 6.8; 100 °C/pH 6.8; 100 °C/pH 5.5), each
with an error bar and a 95% CI on the regression intercept. **There are NO rate constants, NO
activation energies, NO reaction scheme, and NO kinetic model in this paper.** It also supplies
three secondary results that matter for any browning observable: (a) the **non-dialysable (>3500 Da)
fraction carries only 12–18% of the total colour** — i.e. the ε measured here describes a minority
of the chromophore population; (b) ε is **constant over the whole heating period** in all three
systems, so it is degree-of-polymerisation-independent; and (c) elemental analysis giving a
melanoidin C/N ratio (11–19) and stoichiometry (≈1.2 sugar per amino acid, ≈3 water eliminated per
sugar). **The critical caveat, stated by the authors themselves, is that ε here is defined on the
>3500 Da fraction only, whereas the Martins kinetic papers apply it to the TOTAL A₄₇₀** — a
population mismatch that the papers never reconcile. See the DEFECT box in §3.4.

---

## 2. SYSTEM DEFINITION — verbatim

| variable | value as printed | anchor |
|---|---|---|
| reactants | glucose + glycine, **equimolar, final concentration 0.2 M each** | §2.2, p. 136 |
| radiotracer | "**1 MBq of D-[U-¹⁴C]-glucose** was spiked into the solution (**the concentration of added labelled sugar was negligible**)". Supplied by Amersham Life Science Ltd (UK); "Specific Activity **111 MBq/mmol and 11.4 GBq/mmol**" | §2.1–2.2, p. 136 |
| buffer | **phosphate buffer 0.1 M**, pH 6.8 or 5.5 | §2.2, p. 136 |
| pH | 6.8 or 5.5. **Not stated whether hot-pH corrected; not held constant; no pH-drift data given in this paper.** (Part I, the companion, measures the drift on the same buffer: 0.19–0.21 units at pH 6.8, 0.8–0.9 units at pH 5.5.) | §2.2, p. 136 |
| temperature | 120 °C and 100 °C (pH 6.8), and 100 °C (pH 5.5) — **oil bath**, "the proper safety measures taken" | §2.2, p. 136 |
| system labels | **A** = 120 °C, pH 6.8; **B** = 100 °C, pH 6.8; **C** = 100 °C, pH 5.5 | §1, p. 136 |
| volume | reaction mixtures **100 ml**, distributed over tubes "each containing a minimum of **10 ml**" | §2.2, p. 136 |
| vessel | "glass, screw-capped, Schott tubes (**16 × 160 mm**)" | §2.2, p. 136 |
| headspace / atmosphere | **not stated** | — |
| agitation | **not stated** | — |
| quench | "At timed intervals, the samples were withdrawn, **immediately cooled in ice** and then dialysed." | §2.2, p. 136 |
| replication (n) | "Samples were heated in **at least duplicate**." For ε itself: "The extinction coefficient was calculated taking **the average of the repetitions** carried out under the same conditions." Fig. 3 shows **two regression lines per system**, i.e. n = 2 repetitions. | §2.2, p. 136; §3.1, p. 138; Fig. 3 caption |
| error definition, ε | **not explicitly defined.** ε is the slope of A₄₇₀ vs incorporated-sugar concentration; the "±" is most plausibly the standard error of that regression slope, but the paper never says. Fig. 2's caption does define its own bars: "**The error bars represent the standard deviation for each observation.**" | Table 1; Fig. 2 caption, p. 138 |

### Melanoidin definition and separation — the operative methodological choice

> "Melanoidins were, **rather arbitrarily**, defined as high-molecular-weight by a **lower limit of
> 3500 Da**, which was the nominal cut-off value in the dialysis system used." (§1, p. 136)

**Dialysis — §2.3, p. 136.** "Approximately **2 ml** of the reaction mixture were injected into
dialysis cassettes (Mr > 3500) (Slide-A-Lyzer Dialysis Cassette, **3.5 K MWCO**, Pierce, USA) and
dialyzed against distilled water. The optimum dialysis time was established by carrying out a test,
where the retentate of the same sample was counted after different dialysis days. Each day
corresponds to **two water changes, 2 l each**. The ¹⁴C-activity reached equilibrium after
**7 days (14 water replacements)** (Fig. 1). Only then the contents were removed and the volume
adjusted to **10 ml** with distilled water."

**Scintillation counting — §2.4, pp. 136–137.** 1 ml of diluted dialysed retentate + **10 ml
scintillation liquid** (Emulsifier Scintillator Plus, Packard) in a plastic vial; counted
**10 min** on a Liquid Scintillation Analyzer 1600TR (Packard). "The count due to ¹⁴C was
**corrected for quenching** … **Quench correction was done by the internal standard method**
(Price, 1973). After counting the sample was spiked with an accurately known amount of isotope in
an unquenched form and recounted. The specific activity of ¹⁴C-glucose in the reaction mixture was
calculated from the counts obtained from 1 ml of a **100-fold diluted unheated reaction mixture**
and was expressed as number of **disintegration per minute (dpm) per mol of glucose**. Once the
quench-corrected number of counts for a certain sample was known the concentration of ¹⁴C-sugar
incorporated into the high molecular weight fraction could be calculated by dividing the number of
counts per minute by the specific activity of the sugar."

**Spectrometry — §2.5, p. 137.** "Browning was measured spectrophotometrically as the absorbance at
**420, 450, 470 and 490 nm. Optical pathlength of the cuvette 1 cm. Spectrophotometer Cary 50 Bio,
Varian.**"

**Microanalysis — §2.6, p. 137.** "The reaction mixtures (**without addition of radio labelled
sugar**) were dialyzed in the same way as described above and the retentates **freeze-dried**.
Microanalysis was carried out using a **CE Instruments Element Analyzer, Type EA 1110 CHN**."

### How ε is defined, verbatim (§3.1, p. 138)

> "the extinction coefficient was calculated based on the Lambert–Beer equation (A = ε · c · l).
> Since the factor l is constant, there is a direct linear relation between absorbance (A) and
> concentration (c), through the extinction coefficient (ε). **The concentration of nondialysable
> melanoidins was expressed as the concentration of U-¹⁴C glucose incorporated. The extinction
> coefficient was then the slope of the plot of the absorbance of the retentate after dialysis
> versus the concentration of radiolabelled glucose incorporated into those melanoidins.**"

Why the intercepts are non-zero, verbatim:
> "The finding that **the regression line not always passed through the origin** was ascribed to the
> background measured in the system. **Radioactive glucose, not incorporated into the melanoidins
> was possibly still retained inside of the retentate** to some extent. A possible reason is the
> **fouling of the membrane** by the melanoidins, slowing down the diffusion, or a **complex
> formation of the melanoidins with the sugar**, which would also delay the diffusion."

---

## 3. TABLE 1 — THE EXTINCTION COEFFICIENTS (complete transcription)

**Anchor: Table 1, p. 139 (PDF p. 5).** Title as printed:
> "Extinction coefficient (ε) of glucose/glycine nondialysable melanoidins measured for different
> wavelengths under different reaction conditions"

Column headers as printed: a spanning header `Wavelength (nm)` over four columns
`420` | `450` | `470` | `490`.
Row labels as printed: `ε (l mmol⁻¹ cm⁻¹)`, `Intcp.ᵃ Lower 95%`, `Intcp.ᵃ Upper 95%`, under three
italic condition headings `120 °C, pH 6.8`, `100 °C, pH 6.8`, `100 °C, pH 5.5`.
Footnote as printed: "ᵃ Intercept with the origin (95% confidence interval)."

**All ε values are [F] — regression slopes fitted by the authors** (least-squares slope of A vs.
incorporated-sugar concentration). The intercept bounds are likewise [F]. The units printed are
**l mmol⁻¹ cm⁻¹** for ε; the intercept rows carry **no unit** in the original (they are in
absorbance units, dimensionless, by construction — **the missing unit is a printing omission**).

**Read from the PDF p. 5 raster.** The `pdftotext` text layer of this table is **unusable**: it
drops every minus sign, turning e.g. "−0.23" into "0.23". Every sign below was verified visually.

### System A — 120 °C, pH 6.8

| row | 420 nm | 450 nm | 470 nm | 490 nm |
|---|---|---|---|---|
| ε (l mmol⁻¹ cm⁻¹) | **1.00 ± 0.03** | **0.77 ± 0.02** | **0.65 ± 0.02** | **0.53 ± 0.01** |
| Intcp. Lower 95% | −0.23 | −0.18 | −0.15 | −0.12 |
| Intcp. Upper 95% | −0.02 | −0.002 | −0.01 | −0.002 |

### System B — 100 °C, pH 6.8

| row | 420 nm | 450 nm | 470 nm | 490 nm |
|---|---|---|---|---|
| ε (l mmol⁻¹ cm⁻¹) | **1.01 ± 0.02** | **0.79 ± 0.02** | **0.66 ± 0.02** | **0.53 ± 0.01** |
| Intcp. Lower 95% | −0.12 | −0.09 | −0.09 | −0.057 |
| Intcp. Upper 95% | −0.02 | −0.01 | −0.004 | **0.01** |

### System C — 100 °C, pH 5.5

| row | 420 nm | 450 nm | 470 nm | 490 nm |
|---|---|---|---|---|
| ε (l mmol⁻¹ cm⁻¹) | **0.97 ± 0.07** | **0.71 ± 0.05** | **0.62 ± 0.05** | **0.50 ± 0.04** |
| Intcp. Lower 95% | −0.20 | −0.16 | −0.14 | −0.11 |
| Intcp. Upper 95% | **0.07** | **0.05** | **0.05** | **0.05** |

### 3.1 Cross-condition summary [M] and derived spread [Z]

| wavelength | A (120/6.8) | B (100/6.8) | C (100/5.5) | range [Z] | max relative error bar [Z] |
|---|---|---|---|---|---|
| 420 nm | 1.00 ± 0.03 | 1.01 ± 0.02 | 0.97 ± 0.07 | 0.04 (4.0%) | 7.2% (system C) |
| 450 nm | 0.77 ± 0.02 | 0.79 ± 0.02 | 0.71 ± 0.05 | 0.08 (10.7%) | 7.0% (system C) |
| 470 nm | 0.65 ± 0.02 | 0.66 ± 0.02 | 0.62 ± 0.05 | 0.04 (6.3%) | 8.1% (system C) |
| 490 nm | 0.53 ± 0.01 | 0.53 ± 0.01 | 0.50 ± 0.04 | 0.03 (5.8%) | 8.0% (system C) |

**Authors' headline, verbatim (Abstract, p. 135) [M]:**
> "At 470 nm, values of **0.65 (±0.02) l mmol⁻¹ cm⁻¹**; **0.66 (±0.02) l mmol⁻¹ cm⁻¹** and
> **0.62 (±0.05) l mmol⁻¹ cm⁻¹**, were obtained at **120 °C pH 6.8, 100 °C pH 6.8 and 100 °C
> pH 5.5**, respectively. **The difference is not significant.** The extinction coefficient appeared
> to not to vary within the pH and temperature range studied."

**The three-system average, which is the value the kinetic papers actually use (§3.1, p. 138) [F/Z
by the authors]:**
> "In the present study the **molar extinction coefficient averaged over the three systems**, for
> glucose/glycine HMW melanoidins was calculated to be **0.64 ± 0.03 l mmol⁻¹ cm⁻¹ at 470 nm**."

*(Note: the arithmetic mean of 0.65, 0.66, 0.62 is 0.6433 → 0.64. Self-consistent.)*

### 3.2 Constancy over the heating period [M]

> "In Fig. 3 the results show that **over the observation period the extinction coefficient of the
> nondialysable melanoidins remained constant**, in the three systems studied."
> "**If the degree of polymerization of melanoidins increases with the reaction time, the fact that
> the extinction coefficient does not seem to change with time implies that it does not depend on
> the degree of polymerization.**"

This is the "DP-invariant" claim the declaration cites. **It is an inference from the linearity of
Fig. 3, not an independent measurement.** Its strength is the linearity of the plots, which is
excellent by eye (see §5).

### 3.3 Literature comparisons — all [C], NOT measurements of this system

| ε at 470 nm | system | source, as cited |
|---|---|---|
| 0.94 l mmol⁻¹ cm⁻¹ | glucose/glycine, HMW and LMW the same | Leong (1999) PhD thesis, Leeds |
| 0.34 l mmol⁻¹ cm⁻¹ | **alanine** (lower bound of Leong's amino-acid range) | Leong (1999) |
| 1.0 l mmol⁻¹ cm⁻¹ | glucose/glycine, pH 5.5, 55 °C — **estimated as a kinetic model parameter, not radiolabelled** | Leong & Wedzicha (2000) |
| 0.3 l mmol⁻¹ cm⁻¹ (recalculated to 470 nm) | sugar/**casein**, 120 °C, pH 6.8 | Brands et al. (2002) |
| 0.48 and 0.53 l mmol⁻¹ cm⁻¹ **at 420 nm** | glucose/casein and fructose/casein | Brands et al. (2002) |
| 0.41 l mmol⁻¹ cm⁻¹ **at 450 nm**, 90 °C | glucose/glycine, radiolabelled | Wedzicha & Kaputo (1992) |
| 0.37 l mmol⁻¹/cm⁻¹ **at 450 nm**, 55 °C | glucose/glycine, radiolabelled | Wedzicha & Kaputo (1992) — note the printed unit "l mmol⁻¹/cm⁻¹" is a **typo** for l mmol⁻¹ cm⁻¹ |

Authors' comparison, verbatim: "In the present study, under the same conditions the obtained
extinction coefficient was **1.00 ± 0.03 l mmol⁻¹ cm⁻¹ at 420 nm, about twice as high than the ε of
sugar/casein melanoidins.** This means that in sugar/amino acid reaction mixtures, **less glucose
molecules (or glucose fragments) have to be incorporated into the melanoidins** than in the
sugar/casein systems, to increase the absorbance by one unit."

### 3.4 DEFECT BOX — problems a downstream user must not paper over

1. **POPULATION MISMATCH — the most consequential defect.** ε here is measured on the **>3500 Da
   retentate only**, which carries only **12–18%** of the total colour (§4). But **Part I
   (`martins2003b.pdf` §2.10) and Martins & Van Boekel 2005 Food Chem. 92 (`martins2005b.pdf` §2.2)
   apply this ε to the TOTAL, undialysed A₄₇₀** to obtain "melanoidin concentration". **These are
   different populations.** The kinetic papers' melanoidin concentrations are therefore
   *A₄₇₀,total / ε_HMW*, a quantity with no clean physical interpretation unless the LMW and HMW
   chromophores happen to share an ε. Leong (1999) is cited as evidence that they do
   ("at 470 nm the extinction coefficient for HMW and LMW glucose/glycine melanoidins was the
   same") — **but that is a [C] result from a different lab under different conditions (55 °C,
   acetate buffer, pH 5.5), and this paper does not verify it.** Any trunk that scores browning
   against a Martins melanoidin curve inherits this unverified assumption.
2. **Non-zero, systematically NEGATIVE intercepts.** In 10 of the 12 (system × wavelength) cells,
   **both** 95% bounds on the intercept are negative, i.e. the intercept is significantly below
   zero. The authors explain this as retained un-incorporated radioglucose in the retentate, which
   would bias the *x*-axis (concentration) **upward** and hence bias **ε downward**. **The reported
   ε values are therefore plausibly lower bounds.** The authors do not correct for this.
3. **Error bars are undefined for ε.** Table 1's "±" is never identified as SD, SE or CI. Fig. 2's
   caption defines *its* bars as SD, but that is a different quantity. **Treat the ε "±" as an
   unlabelled dispersion measure.**
4. **The intercept rows carry no unit** in the original. They are absorbance units by construction.
5. **Only n = 2 repetitions per system** (two lines per panel in Fig. 3). A "not significant"
   verdict on a three-way comparison at n = 2 is weak. The pH 5.5 system in particular has error
   bars 2.5× the others, which is what one would expect from its much smaller absorbance range
   (see §5) rather than from any real invariance.
6. **The pH 5.5 system was run to 480 min while pH 6.8 systems ran to 60–300 min** (Tables 3–4).
   The three systems are matched on *absorbance reached*, not on time. This is deliberate and
   sensible, but it means the "temperature/pH invariance" conclusion is really an
   "invariance-with-extent-of-browning" conclusion.

---

## 4. THE ¹⁴C / COLOUR-PARTITION RESULTS — % of glucose and % of colour

**Anchor: §3.1, pp. 137–138, and Abstract, p. 135.** No table — these are body-text numbers and
Fig. 2 / Fig. 4 readings.

### 4.1 Fraction of colour retained above 3500 Da [M]

> "After dialysis the absorbance results show that independently of the reaction conditions,
> **more than 80% of browning, measured spectrophotometrically at 470 nm, passed into the
> dialysate, namely 88, 88 and 82% in system A, B and C, respectively.** The majority of coloured
> compounds were thus not retained in the HMW fraction (Mr > 3500) but below it."

| system | conditions | % colour passing into dialysate (LMW) [M] | % colour retained (>3500 Da, HMW) [Z] |
|---|---|---|---|
| A | 120 °C, pH 6.8 | 88 | 12 |
| B | 100 °C, pH 6.8 | 88 | 12 |
| C | 100 °C, pH 5.5 | 82 | 18 |

Abstract, verbatim [M]: "These melanoidins only represented **12% of the total colour formed**."
**Note the abstract quotes 12% and silently omits system C's 18%.** Flagged.

### 4.2 HMW-vs-LMW proportionality [F]

**Anchor: Fig. 4, p. 140.** Caption verbatim: "The absorbance at 470 nm of high molecular weight
(HMW) fraction plotted against the low molecular weight (LMW) fraction. **A—Temperature 120 °C
pH 6.8 (HMW = 0.12 LMW; R² = 0.99); B—Temperature 100 °C pH 6.8 (HMW = 0.14 LMW; R² = 0.99);
C—Temperature 100 °C pH 5.5 (HMW = 0.24 LMW; R² = 0.96)**."

| system | HMW/LMW slope [F] | R² [F] | HMW as % of total, from the slope [Z] |
|---|---|---|---|
| A (120 °C, pH 6.8) | 0.12 | 0.99 | 10.7% |
| B (100 °C, pH 6.8) | 0.14 | 0.99 | 12.3% |
| C (100 °C, pH 5.5) | 0.24 | 0.96 | 19.4% |

**Internal consistency check [Z]:** the slope-derived HMW shares (10.7 / 12.3 / 19.4%) agree with
the §4.1 direct measurements (12 / 12 / 18%) to within ~1.5 percentage points. **Consistent.**
Note that the *lowest*-severity system (100 °C, pH 5.5) has the **largest** HMW share — i.e. the
polymer/oligomer partition does depend on condition, even though ε does not. The authors' own
sentence — "The fact that independently of the studied reaction conditions, the same % of high
molecular weight fraction was obtained, suggests that the formation of nondialysable (HMW)
melanoidins does not appear to be sensitive to temperature and pH" — **overstates this**: 18–19%
is 1.5× 12%, which is not "the same %". **Flagged as an author overstatement contradicted by their
own Fig. 4 slopes.**

### 4.3 % of glucose incorporated into melanoidins

> **THIS PAPER DOES NOT REPORT A "% OF GLUCOSE INCORPORATED INTO MELANOIDINS" AS A PERCENTAGE OF
> THE STARTING GLUCOSE.** The wave brief asked for it; I searched the whole text and it is not
> there.

What the paper *does* give is the **absolute** incorporated-sugar concentration, as the x-axis of
Fig. 3 — and only as a figure. The maximum x-values reached, digitised from the PDF p. 5 raster
(**[fig]**, read precision ≈ ±0.05 mmol/l):

| system | x-axis label | x range shown | max incorporated sugar reached [fig] | as % of [glucose]₀ = 200 mmol/l [Z] |
|---|---|---|---|---|
| A (120 °C, pH 6.8) | mmol/l | 0.0–3.0 | ~2.85 mmol/l | ~1.4% |
| B (100 °C, pH 6.8) | mmol/l | 0–2 | ~1.65 mmol/l | ~0.8% |
| C (100 °C, pH 5.5) | mmol/l | 0–1.5 | ~1.25 mmol/l | ~0.6% |

**So: of 200 mmol/l starting glucose, roughly 0.6–1.4% ends up in the >3500 Da melanoidin
fraction under these conditions.** [Z], derived by me from a [fig] read against the stated 0.2 M
initial concentration. **This number is never printed by the authors and must carry both tags.**
If total colour is ~6–8× the HMW colour (§4.1–4.2), a *total*-melanoidin glucose incorporation
would be of order **4–11%** [Z] — but that extrapolation assumes the LMW/HMW ε equality that
this paper does not verify (DEFECT box item 1), so it should be treated as an order-of-magnitude
statement only.

### 4.4 Browning induction times [M]

**Anchor: §3.1, p. 137, describing Fig. 2.**
> "An absorbance of **three units** was reached after **0.5, 1.8 and 8 h** for systems A (120 °C,
> pH 6.8), B (100 °C, pH 6.8) and C (100 °C, pH 5.5), respectively."

| system | time to A₄₇₀ = 3 (before dialysis) [M] |
|---|---|
| A (120 °C, pH 6.8) | 0.5 h (30 min) |
| B (100 °C, pH 6.8) | 1.8 h (108 min) |
| C (100 °C, pH 5.5) | 8 h (480 min) |

[Z] ratios: A:B = 3.6× for +20 °C at pH 6.8; B:C = 4.4× for +1.3 pH units at 100 °C. **Compare
Part I's finding for DFG degradation that −1.3 pH units ≡ +20 °C. The same equivalence holds here
for browning, to within 20%.** A useful independent corroboration.

### 4.5 Oligomer size [M/C]

> "the result that HMW colourants with molecular weight up to several thousand daltons could not be
> observed (Hofmann, 1998b) suggests that the built up of LMW into HMW structures **only reaches
> the oligomer size, approximately 13 molecules of glucose and glycine**. Once the melanoidin
> chromophore is formed its concentration increases proportionally with time."

**This "13" is an inference the authors draw from a cited result, not a measurement.** Tag as
[C]-derived. (Sanity: 13 × ~180 Da ≈ 2340 Da, just below the 3500 Da cut-off — the reasoning is
that anything retained is barely over the cut-off.)

---

## 5. FIGURES 1–5 — inventory and digitised values

### 5.1 Fig. 1 — dialysis equilibration
**Anchor: Fig. 1, p. 137 (PDF p. 3).** Caption verbatim: "Change in ¹⁴C-activity in the retentate as
a function of dialysis time (two water changes per day) at room temperature. Counts per minute
(cpm)." Axes: y = cpm, x = dialysis time (days). **Axis tick values not legible at page-raster
resolution: unreadable.** The operative result is in §2.3 text: **equilibrium after 7 days /
14 water replacements** [M].

### 5.2 Fig. 2 — browning before and after dialysis
**Anchor: Fig. 2, p. 138 (PDF p. 4).** Caption verbatim: "Browning [measured as absorbance (Abs) at
470 nm] before (closed markers) and after (open markers) dialysis. A—Temperature 120 °C pH 6.8;
B—Temperature 100 °C pH 6.8; C—Temperature 100 °C pH 5.5. **The error bars represent the standard
deviation for each observation.**" Three stacked panels, y = Abs 470 nm, x = time. Per-point
digitisation not performed (the operative numbers — 88/88/82% and the 0.5/1.8/8 h induction times —
are printed in the text and transcribed in §4.1 and §4.4 above). **This is the only figure in the
paper whose error bars are explicitly defined.**

### 5.3 Fig. 3 — the ε regressions (the primary data)
**Anchor: Fig. 3, p. 139 (PDF p. 5).** Caption verbatim: "Browning [measured as absorbance (Abs) at
470 nm] as function of the melanoidin concentration (measured as incorporated sugar);
A—Temperature 120 °C pH 6.8; B—Temperature 100 °C pH 6.8; C—Temperature 100 °C pH 5.5. **Each line
corresponds to a repetition of the experiment. The slope of the lines determines the extinction
coefficient value.**"

Axes as printed, read from the PDF p. 5 raster:

| panel | y-axis label & range | x-axis label & range | n lines | shape [fig] |
|---|---|---|---|---|
| A (120 °C, pH 6.8) | `Abs 470 nm`, 0–2, ticks 0.4 | `mmol/l`, 0.0–3.0, ticks 0.5 | 2 (nearly superimposed) | linear above ~0.5 mmol/l; **visible upward curvature / offset below ~0.4 mmol/l** |
| B (100 °C, pH 6.8) | `Abs 470 nm`, 0–1.2, ticks 0.2 | `mmol/l`, 0–2, ticks 0.5 | 2 (nearly superimposed) | linear throughout |
| C (100 °C, pH 5.5) | `Abs 470 nm`, 0–1, ticks 0.2 | `mmol/l`, 0–1.5, ticks 0.5 | 2 (nearly superimposed) | linear; **noticeably more scatter**, consistent with the larger ± in Table 1 |

**Terminal points [fig]** (PDF p. 5 raster, ±0.05 mmol/l on x, ±0.05 on y): A reaches
(~2.85 mmol/l, ~1.75 Abs); B reaches (~1.65 mmol/l, ~1.02 Abs); C reaches (~1.25 mmol/l,
~0.83 Abs).

**Slope sanity check [Z]:** A: 1.75/2.85 = 0.61; B: 1.02/1.65 = 0.62; C: 0.83/1.25 = 0.66. Against
the printed ε₄₇₀ of 0.65, 0.66, 0.62. **Agreement to within ~7%, i.e. within my read precision.
The Table 1 values are corroborated by the figure.** (My through-origin ratios do not account for
the negative intercepts, which is why they run slightly low for A and B — exactly the direction
DEFECT box item 2 predicts.)

**FLAG — panel A's low-concentration curvature.** Below ~0.4 mmol/l the panel-A points sit visibly
above a straight line through the rest. If real, ε is **not** constant at the very start of
browning in the hottest system, which would qualify the "constant over the observation period"
claim. The authors do not mention it. **[fig] observation, low confidence, but it is there.**

### 5.4 Fig. 5 — schematic only
**Anchor: Fig. 5, p. 140.** Caption verbatim: "Formation of high molecular weight components through
the combination of coloured low molecular weight subunits (**adapted from Leong, 1999**). R and R′
may carry a chromophore which does not change extensively when combined with the high molecular
weight molecule." **A cartoon. No numbers. [C].**

---

## 6. TABLE 2 — literature C/N ratios (all [C], none measured here)

**Anchor: Table 2, p. 141 (PDF p. 7).** Title as printed: "Microanalysis results reported in
literature under different reaction conditions". Column headers as printed: `Reference` |
`Reaction conditions` | `C/N`.

| Reference | Reaction conditions (verbatim) | C/N |
|---|---|---|
| Cämmerer and Kroh (1995) | H₂O; 60 °C; 160 h; pH 5; [gly] = [glu] = 0.1 M | 7 |
| Cämmerer and Kroh (1995) | H₂O; 100 °C; 10 h; pH 5; [gly] = [glu] = 0.1 M | 9 |
| Wedzicha and Kaputo (1992) | Acetate Buffer 0.2 M; 90 °C; 22 h; pH 5.5; [gly] = 1.0 M; [glu] = 1.0 M | 8 |
| Leong (1999) | Acetate Buffer 0.2 M; 55 °C; 90 h; pH 5.5; [gly] = 0.5 M; [glu] = 1.0 M | 8 |
| Feather and Nelson (1984) | 100 °C; 8 h; pH 3.5; [gly] = 1.0 M; [glu] = 0.2 M | 10 |
| Bobbio et al. (1981) | Citrate Buffer 0.05 M; [gly] = 0.66 M; [glu] = 1.25 M; 70 °C; 415 h; pH 3.0 | 12 |
| Bobbio et al. (1981) | (same) 70 °C; 415 h; pH 6.0 | 11 |
| Bobbio et al. (1981) | (same) 80 °C; 80 h; pH 3.0 | 13 |
| Bobbio et al. (1981) | (same) 80 °C; 80 h; pH 6.0 | 10 |
| Hayashi and Namiki (1986) | [ala] = [glu] = 2 M; 100 °C; 190 min; pH 2.3 | 12 |
| Hayashi and Namiki (1986) | (same) 100 °C; 53 min; pH 6.5 | 9 |
| Hayashi and Namiki (1986) | (same) 100 °C; 25 min; pH 9.2 | 7.5 |
| Olsson, Pernemalm, and Theander (1978) | H₂O; 100 °C; 120 h; pH 5; [gly] = 0.5 M; [glu] = 0.75 M | 12 |

**Every row is [C].** No error bars, no units beyond those shown. The authors' own comment:
"the C/N values reported in literature are **not consistent, either with pH or temperature**."

---

## 7. TABLE 3 — microanalysis, system A (120 °C, pH 6.8)

**Anchor: Table 3, p. 141 (PDF p. 7).** Title as printed: "Microanalysis results from system Aᵃ".
Footnote as printed: "ᵃ Calculated number of molecules of sugar (a) per molecule of amino acid (b)
and calculated number of molecules of water (y) per molecule of sugar (a)."
Column headers as printed: `Time (min)` | (spanning) `120 °C pH 6.8` over `C/N` | `a/b` | `y/a` |
`Abs 470 nm`.

| Time (min) | C/N | a/b | y/a | Abs 470 nm |
|---|---|---|---|---|
| 15 | 11 | 1.2 | 1.8 | 3.6 |
| 30 | 11 | 1.2 | 2.6 | 6.3 |
| 45 | 11 | 1.2 | 3.1 | 8.1 |
| 60 | 11 | 1.2 | 3.2 | 9.2 |

C/N and Abs are **[M]**; a/b and y/a are **[F]** — fitted by solving the stoichiometric model below.
No error column is printed.

**The stoichiometric model, verbatim (§3.2, p. 140):**
> "The overall reaction for the formation of melanoidin is a combination of **a** molecules of sugar
> consisting of l, m and n atoms of C, H and O, respectively and **b** molecules of amino acid
> consisting of p, q, r and s atoms of C, H, N and O, to give a melanoidin formula, where **y** is
> the number of water molecules: **C_{la+pb} H_{ma+qb−2y} N_{rb} O_{na+sb−y}**. In the
> glucose/glycine system **l = 6, m = 12, n = 6, p = 2, q = 5, r = 1 and s = 2. Assuming b = 1**,
> the unknowns a and y can be found by solving the following equations: **C = 6a + 2** and
> **H = 12a + 5 − 2y**."
> "**The number of carbon dioxide molecules was not calculated.**"

Authors' reading [M]: "the number of incorporated mol of sugar (or its corresponding degradation
product) per amino acid **remains constant, around 1.2**. This is consistent with the fact that
**almost 80% of glycine was recovered** after the reaction heating period."
And: "the estimated number [of water molecules eliminated per mol of sugar] was **3**."

**INTERNAL INCONSISTENCY, FLAGGED:** the text says y/a "was 3", but Table 3's y/a column is
**1.8 → 2.6 → 3.1 → 3.2**, i.e. it rises by 78% over the heating period and only reaches ~3 at the
end. **The "constant, independent of reaction conditions" claim for y/a is not supported by the
authors' own table.** Quote both, as the spec requires:
- Text (p. 140): "the number of molecules of water eliminated per mol of sugar or corresponding
  degradation product incorporated **seemed to be constant independently of the reaction
  conditions**. In the present study the estimated number was **3**."
- Table 3 (p. 141): y/a = **1.8, 2.6, 3.1, 3.2** at 15, 30, 45, 60 min.

Literature comparisons for a/b, all [C]: **0.91** (Yaylayan & Kaminsky 1998, H₂O/MeOH, 65 °C, 7 h)
and **2.19** (Cämmerer & Kroh 1995, solvent-free, 170 °C, 20 min).

---

## 8. TABLE 4 — microanalysis, systems B and C

**Anchor: Table 4, p. 141 (PDF p. 7).** Title as printed: "Microanalysis results from systems B and
C". Column headers as printed: `Time (min)` | (spanning) `100 °C pH 6.8` over `C/N` | `Abs 470 nm` |
(spanning) `100 °C pH 5.5` over `C/N` | `Abs 470 nm`.

| Time (min) | 100 °C pH 6.8: C/N | 100 °C pH 6.8: Abs 470 nm | 100 °C pH 5.5: C/N | 100 °C pH 5.5: Abs 470 nm |
|---|---|---|---|---|
| 30 | 15 | 0.5 | *(blank)* | *(blank)* |
| 60 | 14 | 2.0 | 19 | 0.1 |
| 120 | 11 | 4.6 | *(blank)* | *(blank)* |
| 180 | 11 | 6.3 | 16 | 0.6 |
| 300 | 11 | 10.0 | 12 | 1.3 |
| 420 | *(blank)* | *(blank)* | 12 | 2.2 |
| 480 | *(blank)* | *(blank)* | 12 | 2.6 |

All **[M]**. Blank cells are blank in the original (those time points were not sampled in that
system) — **not zeros**. No a/b or y/a columns are given for systems B and C; **only system A has
the full stoichiometric analysis.** No error columns.

Authors' reading [M] (§3.2, pp. 140–141):
> "at higher pH and temperature, in a way the most favourable reaction conditions for formation of
> melanoidins, **the C/N ratio remained constant throughout the heating period**. If we decreased
> the temperature and/or the pH, we observed in the initial period of the reaction **a higher value
> for C/N followed by a decrease till a constant value was reached**. It seems that the ratio is
> indeed **lower at higher pH** and therefore seems to depend on the reaction conditions. These
> results are consistent with the fact that at **lower pH around 90% of glycine was recovered**, in
> contrast with the **80% at higher pH**, independently of the temperature. However, **we should be
> careful with interpretation because elemental analysis could be sensitive to experimental
> errors.**"

Abstract [M]: "A trend was observed in the melanoidins C/N ratio: **it decreased with increasing
reaction pH** as well as it **changed to a lower level, of about 8, as the extent of browning
increased**."

**INTERNAL INCONSISTENCY, FLAGGED:** the abstract says C/N "changed to a lower level, **of about
8**, as the extent of browning increased." **No cell in Tables 3 or 4 is 8, or even close: the
minimum measured value anywhere is 11.** The value 8 appears only in **Table 2**, i.e. in the
*literature* rows (Wedzicha & Kaputo 1992; Leong 1999). **The abstract imports a literature number
into a sentence about "the melanoidins C/N ratio" without saying so.** Quote both:
- Abstract (p. 135): "it changed to a lower level, of **about 8**, as the extent of browning
  increased."
- Tables 3 and 4 (p. 141): the measured C/N values are **11, 11, 11, 11 / 15, 14, 11, 11, 11 /
  19, 16, 12, 12, 12**. Minimum = **11**.
- Body text (p. 141) is careful and correct: "there is a trend that the C/N ratio decreases with
  increasing pH (Bobbio…; Hayashi & Namiki) **as well as it changes to a lower level, of about 8,
  as the extent of browning increases**" — the citations make clear this half-sentence is about the
  literature. **The abstract drops the attribution.**

---

## NEW-PARAMETER TABLE (consolidated)

**This paper contains NO rate constants, NO equilibrium constants and NO activation energies.**

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| ε₄₂₀ | 1.00 ± 0.03 | l mmol⁻¹ cm⁻¹ | 120 °C, pH 6.8, glu/gly 0.2 M equimolar, 0.1 M phosphate, >3500 Da fraction | Table 1, p. 139 | [F] |
| ε₄₅₀ | 0.77 ± 0.02 | l mmol⁻¹ cm⁻¹ | 120 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₇₀ | 0.65 ± 0.02 | l mmol⁻¹ cm⁻¹ | 120 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₉₀ | 0.53 ± 0.01 | l mmol⁻¹ cm⁻¹ | 120 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₂₀ | 1.01 ± 0.02 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₅₀ | 0.79 ± 0.02 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₇₀ | 0.66 ± 0.02 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₉₀ | 0.53 ± 0.01 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 6.8 | Table 1, p. 139 | [F] |
| ε₄₂₀ | 0.97 ± 0.07 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 5.5 | Table 1, p. 139 | [F] |
| ε₄₅₀ | 0.71 ± 0.05 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 5.5 | Table 1, p. 139 | [F] |
| ε₄₇₀ | 0.62 ± 0.05 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 5.5 | Table 1, p. 139 | [F] |
| ε₄₉₀ | 0.50 ± 0.04 | l mmol⁻¹ cm⁻¹ | 100 °C, pH 5.5 | Table 1, p. 139 | [F] |
| **ε₄₇₀, three-system average — the value the kinetic papers use** | **0.64 ± 0.03** | l mmol⁻¹ cm⁻¹ | average over A, B, C | §3.1, p. 138 | [F] (author-computed) |
| Intercept 95% CI, 120 °C pH 6.8 | 420: [−0.23, −0.02]; 450: [−0.18, −0.002]; 470: [−0.15, −0.01]; 490: [−0.12, −0.002] | (absorbance; **unit omitted in original**) | 120 °C, pH 6.8 | Table 1, p. 139 | [F] |
| Intercept 95% CI, 100 °C pH 6.8 | 420: [−0.12, −0.02]; 450: [−0.09, −0.01]; 470: [−0.09, −0.004]; 490: [−0.057, +0.01] | (absorbance) | 100 °C, pH 6.8 | Table 1, p. 139 | [F] |
| Intercept 95% CI, 100 °C pH 5.5 | 420: [−0.20, +0.07]; 450: [−0.16, +0.05]; 470: [−0.14, +0.05]; 490: [−0.11, +0.05] | (absorbance) | 100 °C, pH 5.5 | Table 1, p. 139 | [F] |
| ε constancy with time | "remained constant" in all three systems | — | all | §3.1, p. 138 + Fig. 3, p. 139 | [M]/[fig] inference |
| ε independence of degree of polymerisation | asserted | — | all | §3.1, p. 138 | **inference, not a measurement** |
| colour passing into dialysate (<3500 Da) | 88 / 88 / 82 | % of A₄₇₀ | A / B / C | §3.1, p. 137 | [M] |
| colour retained above 3500 Da | 12 / 12 / 18 | % of A₄₇₀ | A / B / C | §3.1, p. 137 | [Z] (100 − above) |
| HMW/LMW absorbance slope | 0.12 (R² 0.99) / 0.14 (R² 0.99) / 0.24 (R² 0.96) | dimensionless | A / B / C | Fig. 4 caption, p. 140 | [F] |
| time to A₄₇₀ = 3 (undialysed) | 0.5 / 1.8 / 8 | h | A / B / C | §3.1, p. 137 | [M] |
| max sugar incorporated into >3500 Da fraction | ~2.85 / ~1.65 / ~1.25 | mmol/l | A / B / C | Fig. 3, p. 139 | **[fig]** (±0.05 mmol/l) |
| **% of starting glucose incorporated into >3500 Da melanoidins** | **~1.4 / ~0.8 / ~0.6** | % of 200 mmol/l | A / B / C | Fig. 3, p. 139 + §2.2 | **[Z] from [fig] — NEVER PRINTED BY THE AUTHORS** |
| glucose per glycine in melanoidin (a/b) | 1.2 (constant, 15–60 min) | mol/mol | 120 °C, pH 6.8 | Table 3, p. 141 | [F] |
| water eliminated per sugar (y/a) | 1.8 / 2.6 / 3.1 / 3.2 at 15/30/45/60 min; **text claims "3", constant** | mol/mol | 120 °C, pH 6.8 | Table 3, p. 141 vs §3.2 p. 140 | [F] — **table contradicts text, see §7** |
| C/N, system A | 11 at all four times | mol/mol | 120 °C, pH 6.8, 15–60 min | Table 3, p. 141 | [M] |
| C/N, system B | 15 → 14 → 11 → 11 → 11 | mol/mol | 100 °C, pH 6.8, 30–300 min | Table 4, p. 141 | [M] |
| C/N, system C | 19 → 16 → 12 → 12 → 12 | mol/mol | 100 °C, pH 5.5, 60–480 min | Table 4, p. 141 | [M] |
| C/N "about 8" | 8 | mol/mol | **NOT MEASURED HERE** — appears only in Table 2 literature rows | Abstract p. 135 vs Table 2 p. 141 | **[C] mis-presented as [M] in the abstract; see §8** |
| glycine recovery | ~90 (lower pH) / ~80 (higher pH) | % | as stated | §3.2, p. 141 | [M] |
| oligomer size | ~13 molecules of glucose + glycine | molecules | inference from Hofmann (1998b) | §3.1, p. 139 | **[C]-derived inference** |
| dialysis equilibration | 7 days / 14 water replacements | days | rt | §2.3, p. 137 + Fig. 1 | [M] |
| dialysis cut-off | 3500 (3.5 K MWCO) | Da | Slide-A-Lyzer, Pierce | §2.3, p. 136 | [M] |
| ¹⁴C specific activity | 111 MBq/mmol and 11.4 GBq/mmol | — | Amersham D-[U-¹⁴C]-glucose | §2.1, p. 136 | [M] |
| literature ε₄₇₀, glu/gly | 0.94 | l mmol⁻¹ cm⁻¹ | Leong (1999), 55 °C, pH 5.5, acetate | §3.1, p. 138 | **[C]** |
| literature ε₄₇₀, glu/gly | 1.0 | l mmol⁻¹ cm⁻¹ | Leong & Wedzicha (2000), **kinetic-model estimate, no radiolabel** | §3.1, p. 138 | **[C]** |
| literature ε₄₇₀, sugar/casein | 0.3 (recalc.) | l mmol⁻¹ cm⁻¹ | Brands et al. (2002), 120 °C, pH 6.8 | §3.1, p. 138 | **[C]** |
| literature ε₄₂₀, sugar/casein | 0.48 (glucose) / 0.53 (fructose) | l mmol⁻¹ cm⁻¹ | Brands et al. (2002) | §3.1, pp. 138–139 | **[C]** |
| literature ε₄₅₀, glu/gly | 0.41 (90 °C) / 0.37 (55 °C) | l mmol⁻¹ cm⁻¹ | Wedzicha & Kaputo (1992) | §3.1, p. 138 | **[C]** |
| literature C/N table | 7 to 13 across 13 condition sets | mol/mol | various | Table 2, p. 141 | **[C]** |
| **rate constants / Ea** | **NONE PRESENT** | — | — | whole paper | negative finding, high confidence |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Blocking issue first: Table 1 of this paper is, on my reading, ALREADY DECLARED.** Amendment 1
assigns "**FIT (measured input, not a scored target)**" to "Martins 2003 thesis Table 3.3.1
melanoidin ε (1.00 ± 0.03 at 420 nm, 0.65 ± 0.02 at 470 nm; T- and DP-invariant)". Those two values
are **verbatim cell matches** to this paper's Table 1, 120 °C / pH 6.8 block. Assessed probability
that these are the same measurements re-published: **very high.**
**Recommendation: do NOT create a second declaration row. Amend the existing one to (i) name the
journal anchor `data/articles/martins2003c.pdf` Table 1, p. 139; (ii) carry all three condition
sets and all four wavelengths, not one cell each; (iii) record that the value the Martins kinetic
papers actually apply is the three-system average 0.64 ± 0.03 at 470 nm, not 0.65 ± 0.02;
(iv) record the >3500 Da population caveat (DEFECT box item 1).**

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 1, ε₄₇₀ all three systems (0.65 / 0.66 / 0.62), and the 0.64 ± 0.03 average** | **FIT — measured input, NOT a scored target** (matches the existing declared role) | — | ε is an observable-model constant that converts the trunk's melanoidin state variable into a predicted A₄₇₀. It is a *calibration*, not a prediction. Scoring the model against it would be circular. This preserves the declaration's existing, correct reasoning: "the Martins step-9 browning HOLD-OUT is unaffected and stays held out." |
| **Table 1, ε at 420 / 450 / 490 nm** | **FIT — measured input**, same status | wavelength | Needed only if the trunk ever predicts browning at a wavelength other than 470 nm. Carrying all four costs nothing and prevents a future wave from re-extracting them as "new". |
| **Table 1 intercept 95% CIs** | **diagnostic_only** | — | They are the quantitative evidence for the retained-radioglucose bias (DEFECT item 2), i.e. evidence that ε is a lower bound. Useful for setting an asymmetric prior on ε; not a target. |
| **§4.1 colour partition (12 / 12 / 18% above 3500 Da) and Fig. 4 slopes (0.12 / 0.14 / 0.24)** | **diagnostic_only** | — | This is the single most important number for anyone who wants to *interpret* a melanoidin concentration: the ε being applied describes ~12–18% of the chromophore population. It should be recorded as a standing caveat on every browning comparison, not fitted. |
| **§4.4 induction times (0.5 / 1.8 / 8 h to A₄₇₀ = 3)** | **HOLD-OUT candidate** | pH and temperature jointly | This *is* a genuine kinetic observable — a browning-rate measurement on the glucose/glycine system at 0.2 M — and it is on a **different concentration basis (0.2 M) from the Martins DFG papers (10 mM)** and from the 2005 pH-ladder paper (0.2 M, matching). It is a clean, cheap, three-point test of whether the trunk gets browning *timescales* right across a 20 °C and 1.3-pH span. **Circularity risk: MEDIUM.** Same lab, and the 0.2 M / pH 6.8 / 100 °C condition is the same nominal system as `martins2005b.pdf`'s 1:1 baseline, whose k9 may be fit. If k9 is fit anywhere, the 1.8 h point is contaminated; the 0.5 h (120 °C) and 8 h (pH 5.5) points remain genuine extrapolations. **Recommend: hold out the 120 °C and pH 5.5 points only.** |
| **Tables 3 and 4 C/N ratios** | **diagnostic_only** | — | Elemental composition is not something the trunk models. The authors themselves warn "elemental analysis could be sensitive to experimental errors", and the abstract mis-attributes a literature value (§8). |
| **Table 3 a/b = 1.2 and y/a** | **diagnostic_only** | — | a/b = 1.2 is a genuinely useful stoichiometric sanity constraint (≈1 glucose per glycine in the polymer, consistent with ~80–90% glycine recovery). y/a is internally contradictory (§7) and must not be used. |
| **Table 2 (literature C/N)** | **not usable** | — | Entirely [C]. Belongs to other papers under other conditions. |
| **The "13-molecule oligomer" figure** | **not usable as a parameter** | — | An inference from a cited absence-of-evidence result. |
| **Fig. 3 x-axis-derived "% glucose into melanoidins" (~0.6–1.4%)** | **diagnostic_only** | — | [Z] from [fig], never printed by the authors, and it describes only the >3500 Da fraction. Useful as an order-of-magnitude reality check on carbon flux into browning; far too soft to score. |

**Cut-axis summary:** this paper does not support a clean fit/hold-out cut of its own, because its
principal deliverable (ε) is a **calibration constant, not a prediction target**. Its correct role
is the one the declaration already assigned: a measured input that the observable model consumes.
The one exception worth harvesting is §4.4's induction times, which *are* predictions and which cut
cleanly on both pH and temperature.

**Circularity risks, stated plainly:**
1. **Highest — the thesis Table 3.3.1 identity.** Re-declaring this paper's Table 1 as a new source
   would double-count a calibration constant, inflating the apparent evidence base.
2. **Structural, and larger than it looks — the population mismatch.** ε is measured on the
   >3500 Da fraction; the Martins kinetic papers apply it to total A₄₇₀. Every Martins "melanoidin
   mmol/l" curve in `martins2003b.pdf` Fig. 8, `martins2003.pdf` Figs. 4–7, `martins2005b.pdf`
   Figs. 5/7/8/9 and `martins2005.pdf` is therefore **A₄₇₀,total divided by ε_HMW**. If the trunk
   fits any of those curves *and* also treats ε as an independently-known input, the two uses are
   not independent, and the ~6–8× population factor is silently absorbed into whatever rate
   constant produces melanoidins. **This should be recorded as a standing caveat on the melanoidin
   observable, whatever roles are eventually assigned.**
3. **Minor — the 0.64 vs 0.65 discrepancy.** The declaration quotes 0.65 ± 0.02 (one system); the
   kinetic papers use 0.64 ± 0.03 (three-system average). The 1.5% difference is negligible
   numerically but matters for traceability: a reader trying to reproduce a Martins melanoidin
   curve with 0.65 will be off by 1.5% and will not know why.
