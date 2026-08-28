# Baek, Linforth, Blake & Taylor 1999 — COMPLETE TRANSCRIPTION
### Model gelatine gels: perceived intensity tracks the RATE of in-nose volatile release, not the maximum concentration — and protein binding is experimentally excluded.

**Full extraction of every number in `data/articles/baek1999.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

**Figures 5 and 7 re-read from a 300 dpi render.** The PDF is encrypted (`copy:no`) but the
text layer and the raster are both fully usable.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **Inger Baek, Rob S. T. Linforth, Anthony Blake¹, Andrew J. Taylor** |
| Title | **"Sensory Perception is Related to the Rate of Change of Volatile Concentration In-nose during Eating of Model Gels"** |
| Venue | ***Chemical Senses* 24: 155–160, 1999** — Oxford University Press |
| Accepted | **9 December 1998** |
| Affiliations | **Samworth Flavour Laboratory, Division of Food Sciences, University of Nottingham, Sutton Bonington Campus, Loughborough LE12 5RD, UK**; ¹ **Corporate R&D, Firmenich SA, Route des Jeunes 1, Geneva 8, Switzerland** |
| Funding | **UK MAFF** and **BBSRC LINK scheme**, industrial partners **Firmenich, Micromass and Stable Micro Systems** |
| PDF character | 6 pages, encrypted (RC4, print-only), clean text layer. ⚠️ **The paper contains exactly ONE table (Table 1, a dosing matrix). Every result is a figure.** |

**Provenance codes:** **[M]** measured and printed · **[F]** fitted and printed ·
**[C]** cited · **[D]** digitised by this wave · **[Z]** derived by this wave.

---

## 1. THE ONE-PARAGRAPH ANSWER

**This paper prints no results table.** Its entire quantitative output is seven figures plus
five p-values, one correlation coefficient and one equation in the running text. What it does
supply is **the single cleanest experimental exclusion of protein binding in the corpus**
(§4) and **a direct, paired, within-subject demonstration that perceived aroma intensity is
NOT proportional to the in-nose volatile concentration** (§5): across five gelatine
concentrations, the **maximum in-nose concentration falls by only 1.5× and is NOT
statistically significant**, while **perceived intensity falls by 3.3× at P < 0.001**. The
quantity that does track perception is the **rate of release** (`Imax/2(t75 − t25)`), which
falls 3.0× at P < 0.01. Digitising Figure 7 and computing the correlations — **which the
paper states qualitatively and never quantifies** — gives **r = 0.968 between the release
gradient and TI Imax versus r = 0.860 for the maximum concentration**, and **r = 0.886 versus
r = 0.851** against the sensory scaling score (§5.2). ⚠️ **The binding exclusion is specific
to gelatine and furfuryl acetate and does not generalise** (§4.3). For the repo, the
load-bearing consequence is that **a model which predicts a headspace or in-nose concentration
and then maps it monotonically to perceived intensity is refuted by a within-subject
experiment**, and the refuting variable is a **time derivative** the repo does not currently
compute.

---

## 2. SYSTEM COMPOSITION AND METHOD

| variable | value as printed |
|---|---|
| **Gel matrix** | **gelatine, 250 bloom strength (DGF STOESS, Germany)** at **2.0, 3.5, 5.0, 6.5, 8.0 % w/w** |
| **Sugars** | **20.5 % w/w sucrose + 22.5 % w/w glucose syrup (Cerestar, UK)** — constant in every gel |
| **Volatile** | **furfuryl acetate (FFA)**, Firmenich SA, **dispersed in propylene glycol** |
| Rationale for FFA | *"a single volatile ... which had an easily recognizable, characteristic flavour so that sensory panellists could readily recognize the compound. Another advantage of FFA was that it shows **relatively low persistence in-mouth** and therefore minimizes problems of carry over."* |
| Preparation | gelatine dissolved at **60 °C**; sucrose–glucose mix **boiled and allowed to cool**; **FFA added to the COOLED gelatine–sugar mixture BEFORE gelation** |
| ⚠️ dosing order | **The flavour is added after the thermal step and before gelation.** Unlike Brewer 1995 (`k2_matrix_and_thresholds.md` §D.2 ii), **there is no thermal step after dosing** — `thermal_step_after_dosing: false`. |
| Portion | cut into cubes of **6.0 ± 0.1 g**, **stored overnight at 4 °C**, then **equilibrated at room temperature** before presentation |
| **Instrument** | **atmospheric pressure ionisation mass spectrometry (API-MS)** — **Platform II quadrupole** (Micromass, Altrincham) with a modified API source [C] Linforth & Taylor 1998, **European Patent EP 0819 937 A2** |
| Ionisation | **4 kV corona discharge, positive ion mode, cone voltage 20 V** |
| **Ion monitored** | **m/z 80.8, a fragment ion of FFA** |
| **Breath sampling** | **air drawn from ONE nostril at 25–30 mL/min over a 2 min period**; *"no specific eating instructions were given"* |
| Signal processing | *"The raw breath by breath traces were converted into TR curves by **smoothing the peak height data** and converting peak height into **actual concentration in air (nanolitres of volatile per litre of air; ppbv)** after **calibration of the API-MS interface with a series of FFA standards**"* |
| **Detection level** | **10 to 100 ppbv in the gas phase** [C] Linforth & Taylor 1998 |
| **TR parameters extracted** | **Imax** (maximum concentration), **Tmax** (time to maximum), and **the rate of release, defined as the gradient over the up-slope between t25 and t75** |
| **Equation (1) as printed** | **`Gradient = Imax / 2(t75 − t25)`** where **t25 = time to 25 % of Imax** and **t75 = time to 75 % of Imax** |
| Rationale for t25–t75 | *"The points t75 and t25 were used because, in most cases, **this region of the TR curve was essentially linear** whereas the situation was more complex above, and below, these values."* |
| **Sensory scaling** | **16 trained panellists**, relative maximum FFA intensity, **two reference samples anchoring the scale: an UNFLAVOURED 5 % w/w gel = 0 and a FLAVOURED 2 % w/w gel = 10**; *"trained to focus on the flavour and **ignore differences in sweetness and texture**"*; **≥2 min between samples**; odourless water and dry crackers available; **each panellist assessed one set of gel samples** |
| **Time–intensity (TI)** | **11 trained panellists**; **TI performed SIMULTANEOUSLY with TR data collection**, the TI signal fed into an **analogue channel of the mass spectrometer**, **one data point per second**; **4 min between samples**; **breath monitored by API-MS before the next gel and no significant carry over observed** |
| TI averaging | *"first normalized in the intensity direction to the numeric average of Imax, followed by averaging in the time dimension"* [C] Overbosch et al. 1986 |
| **Model (in-vitro) system** | **a cube of gel dissolving in 100 mL of distilled water at 37 °C in a 125 mL flask with stirring**; headspace drawn into the API-MS at **12 mL/min**; **monitored in DUPLICATE for up to 45 min** |

### 2.1 Table 1 — the only table. **Anchor: Table 1, p. 156.** [M]

Title as printed: *"Concentrations of gelatine and FFA in the gels"*.

| **Gelatine concentration (% w/w)** | **FFA concentration in gel (mg/kg)** |
|---:|---|
| **2.0** | **100, 300, 600, 1000, 1500** |
| **3.5** | **300** |
| **5.0** | **0, 100, 300, 600, 1000, 1500** |
| **6.5** | **300** |
| **8.0** | **100, 300, 600, 1000, 1500** |

**⇒ 22 gels. The 5 % series is the only one with a true zero-flavour control** (which doubles
as the sensory scale's zero anchor). **All the comparative results use the 300 mg/kg row**,
i.e. the one FFA level present in all five gelatine concentrations.

---

## 3. THE STATISTICAL RESULTS — every p-value and coefficient the paper prints [M]/[F]

**The paper prints exactly six inferential statistics. All are transcribed here; there are no
others.**

| # | statistic as printed | value | test | anchor |
|---:|---|---|---|---|
| 1 | sensory scaling, **2 % vs 8 % w/w gel** significantly different | **P < 0.001** | ANOVA | p. 157 |
| 2 | **TI Imax decreasing** with gelatine concentration | **P < 0.001** | ANOVA | p. 157 |
| 3 | **TI Tmax increasing** with gelatine concentration | **P < 0.001** | ANOVA | p. 157 |
| 4 | **TR Imax across the three gel concentrations** | **NO significant difference** | ANOVA | p. 158 |
| 5 | **instrumental gradient decreasing** | **P < 0.01** | ANOVA | p. 159 |
| 6 | **line of best fit, TI:TR Tmax ratio vs TR Tmax** | **R = 0.421** | linear regression | p. 159 |

**⚠️ NO CORRELATION COEFFICIENT IS PRINTED FOR THE PAPER'S CENTRAL CLAIM.** The claim —
*"A plot of the relative values against gel concentration (Figure 7) suggested that **the
gradient correlated better with the sensory values than TR Imax**"* — is supported only by
*"To confirm this trend, the actual values were analysed statistically"*, followed by the
four ANOVA results above. **The paper never computes a correlation between an instrumental
and a sensory variable.** §5.2 supplies them.

---

## 4. THE PROTEIN-BINDING EXCLUSION — the evidence, in full

This is the requested item and it rests on **two independent experiments**.

### 4.1 The in-vitro dissolution experiment (Figure 3, p. 158) [M]

> *"Gels containing the same concentration of FFA (300 mg/kg) but different gelatine
> concentrations (**2, 5 and 8 % w/w**) were placed in water and stirred, and the FFA release
> was measured by monitoring the HS until the gel had completely dissolved. The results
> (Figure 3) demonstrated that the three gels released FFA at different rates, with the
> **softest gel (2 % w/w gelatine) becoming fully dissolved after 10 min, whereas it took
> 40 min for the harder gel (8 % w/w gelatine) to dissolve completely**. However,
> **irrespective of gelatine concentration, ALL GELS REACHED THE SAME CONCENTRATION OF FFA IN
> THE HS WHEN FULLY DISSOLVED, showing there was NO BINDING of FFA to gelatine under these
> conditions** and that slower release was the most likely explanation for the difference in
> TI performance of the gels."*

**⇒ THE EXCLUSION LOGIC, stated cleanly: at equilibrium after complete dissolution, the
headspace concentration is independent of the protein loading over a 4× range (2 → 8 % w/w,
i.e. 20 → 80 g/L gelatine). If gelatine bound FFA, the 8 % gel would equilibrate at a lower
headspace concentration. It does not.**

**⚠️ The equality is asserted in prose and shown in Figure 3; no numerical tolerance is
given.** The digitisable content of Figure 3 is two duplicate release curves per gel; the
plateau heights are the load-bearing quantity and the paper does not print them.
**Record `binding_exclusion_tolerance: not_stated`.** A generous reading of the figure puts
the three plateaux within ~10 % of each other; **not asserted**.

**Corroboration [C]:** *"This agrees with work by Harrison and Hills (1996), who found that
increasing the gelatine and sucrose concentrations of a gel system resulted in reduced release
rate for a **water-soluble dye** which they used as a marker for gel dissolution."* — i.e. a
**non-volatile, non-binding tracer shows the same rate effect**, which is exactly what a
dissolution-limited (not binding-limited) mechanism predicts.

### 4.2 The in-vivo linearity experiment (Figure 4, p. 158) [M]

> *"there was a possibility that FFA could be binding to **membranes in the mouth and nose** so
> that only a portion of the FFA release was actually transported to the olfactory receptors.
> To test this hypothesis, the in-nose concentration of FFA during eating was monitored for a
> series of gels with different volatile and gel concentrations **using ONE panellist to
> minimize person-to-person variation**. **If no significant binding of the volatile was
> occurring, plotting the maximum concentration of FFA in-nose against the gel FFA
> concentration should produce a LINEAR plot.** Figure 4 shows the data obtained from **three
> replicates of three different gel concentrations with five FFA concentrations**. Despite the
> variable nature of the eating event, **the lines are linear, suggesting that no significant
> binding of FFA occurs between the release event in-mouth and perception in-nose.** The above
> experiment was repeated using **three extra panellists**. The relationships between the gel
> volatile concentration and breath volatile concentration were **still linear for each
> panellist** (data not shown), but there were **substantial differences between panellists in
> the amount of FFA released from the same gels**, presumably because of the different chewing
> patterns of the panellists"* [C] Brown et al. 1995.

**⇒ Nine dose–response lines (3 gel concentrations × 3 replicates) plus four panellists, all
linear over a 15× FFA range (100 → 1 500 mg/kg). Mucosal binding is excluded to the precision
of the linearity.**

**⚠️ This is a stronger and more relevant result than the paper makes of it.** It is a
**four-panellist replication of a dose–response linearity test** at concentrations spanning
15×, and it directly addresses the mucosal-adsorption term that
`k2_matrix_and_thresholds.md` §B.6 flagged as **stranded** for 3-sulfanylhexan-1-ol in saliva.
**For a neutral ester in gelatine, mucosal binding is negligible. For a thiol in saliva,
Starkenmann measured a 6 × 10⁴ analytical quench.** ⇒ **Mucosal/salivary interaction is
compound-class-specific, not a general term.** That is a new, defensible constraint.

### 4.3 ⚠️ THE LIMITS OF THE EXCLUSION — three, all material

1. **One volatile.** FFA is a **neutral, moderately polar ester** with no electrophilic
   carbonyl available for Schiff-base or Michael chemistry. **It is precisely the class of
   compound that `k2_matrix_and_thresholds.md` §B.3 identifies as binding weakly and
   reversibly.** The exclusion says nothing about aldehydes — and Meynier 2002 (this wave)
   measures **6.88× headspace suppression of t-2-hexenal by dairy protein at 34 g/L**, on the
   same class of experiment.
2. **One protein, and an unusually inert one.** **Gelatine is a collagen hydrolysate with
   essentially no hydrophobic binding pockets** — the same point `k2_matrix_and_thresholds.md`
   §D.2(iii) makes about Vega 1994. A null binding result in gelatine does not transfer to
   β-lactoglobulin, casein or soy.
3. **The concentrations are enormous.** **100–1 500 mg/kg FFA** is **10³–10⁶×** above a typical
   odour threshold. **Any saturable binding site is saturated**, so the experiment tests the
   **non-saturable** part of the isotherm only. **A binding site present at, say, 10 µM would
   be invisible at 1 500 mg/kg and dominant at threshold.**

**⇒ Cite this as "binding excluded for a neutral ester in gelatine at 100–1 500 mg/kg", never
as "binding excluded".**

---

## 5. THE CENTRAL RESULT — perception tracks the rate, not the maximum

### 5.1 Figure 5 — the paired TR and TI curves. **Anchor: Figure 5, p. 159.** Digitised. [D]

Caption as printed: *"Volatile release profiles (TR; solid markers) and TI sensory data (open
markers) from ●, 2 %, ▲, 5 % and ■, 8 % w/w gelatine gels flavoured with 300 mg/kg FFA. Each
curve is the mean value obtained from 11 panellists each eating one sample of each gel."*
**Left axis: `[ppbv]`, 0–500. Right axis: `Sensory perception`, 0–12. X axis: `Time [min]`,
0–1.6.**

| gel | **TR Imax (ppbv)** [D ±5 %] | **TR Tmax (min)** [D] | **TI Imax (sensory units)** [D ±5 %] | **TI Tmax (min)** [D] |
|---|---:|---:|---:|---:|
| **2 % w/w** | **≈475** | **≈0.47** | **≈10.2** | **≈0.55** |
| **5 % w/w** | **≈332** | **≈0.57** | **≈6.7** | **≈0.60** |
| **8 % w/w** | **≈313** | **≈0.80–0.85** | **≈3.4** | **≈0.85** |
| **8 %/2 % ratio [Z]** | **0.659 — a 1.52× fall, NOT SIGNIFICANT** | **1.8× later** | **0.333 — a 3.00× fall, P < 0.001** | **1.5× later** |

**⇒ THE HEADLINE CONTRAST, in one line: a 3.0× drop in perceived intensity is accompanied by
a 1.5× drop in in-nose concentration that does not reach significance.**

### 5.2 Figure 7 — the normalised comparison, digitised, and the correlations the paper never computes

Caption as printed: *"Relationship between sensory and instrumental data as gel concentration
changes: **■, TR Imax, ◇, sensory score, ○, TI Imax and ▲, TR gradient**. For a comparison,
all data have been normalized so that the values for the 2 % w/w gel are 100 % w/w. Each point
is the mean value from 11 panellist eating five different gels flavoured with 300 mg/kg FFA."*
**Y axis: `Normalised TI and TR parameters [%]`, 0–120. X axis: `Gelatine concentration [%]`,
0–10.**

| gelatine % | **TR Imax ■** | **sensory score ◇** | **TI Imax ○** | **TR gradient ▲** |
|---:|---:|---:|---:|---:|
| **2.0** | **100** | **100** | **100** | **100** |
| **3.5** | **85** | **96.5** | **62** | **62** |
| **5.0** | **71** | **58.5** | **58.5** | **44.5** |
| **6.5** | **81.5** | **53** | **41** | **41** |
| **8.0** | **67** | **44** | **30.5** | **33.5** |
| **total fall [Z]** | **1.49×** | **2.27×** | **3.28×** | **2.99×** |
| **significance [M]** | **NOT significant** | **P < 0.001** | **P < 0.001** | **P < 0.01** |

*(all values [D], ±3 percentage points; the 2 % row is 100 by construction. Cross-check: the
Figure 5 TR Imax ratio 313/475 = 65.9 % against Figure 7's 67 % ✅, and TI Imax 3.4/10.2 =
33.3 % against 30.5 % ✅ — the two figures are mutually consistent.)*

**⚠️ TR Imax is also NON-MONOTONE: 100 → 85 → 71 → 81.5 → 67.** It rises from 5 % to 6.5 %
gelatine. **The three sensory/rate variables are monotone; the concentration variable is not.**

**THE CORRELATIONS, computed by this wave [Z] — the paper prints none:**

| instrumental variable | **r vs sensory scaling score** | **r vs TI Imax** |
|---|---:|---:|
| **TR Imax (maximum in-nose concentration)** | **+0.851** | **+0.860** |
| **TR gradient (rate of release)** | **+0.886** | **+0.968** |

**⇒ The release gradient beats the maximum concentration on both sensory measures, and
decisively on the like-for-like comparison (TI Imax, both collected simultaneously from the
same 11 panellists): r = 0.968 vs 0.860.** ⚠️ n = 5 points, so neither r is individually
well determined; **the load-bearing evidence is not the correlation but the significance
contrast in §5.1 — the concentration variable does not reach significance at all while the
rate variable does.**

**The paper's own conclusion, verbatim:** *"These findings suggest that, in this system,
**perception is not directly related to the maximum volatile concentration in-nose as might be
predicted by the basic Power Law**. This supports the necessity to modify the Power Law to
take account of temporal changes."* And: *"These experiments demonstrate that **the temporal
aspects of volatile release are related to aroma perception**."*

### 5.3 Figure 6 — the adaptation test, and its weak result [F]

Caption as printed: *"Relationship between the ratio of the sensory Tmax and Instrumental
Tmax plotted against instrumental Tmax."*

> *"If adaptation was occurring, then the receptors should become progressively less receptive
> to an increase in TR so that the TI peak should occur before the TR peak and the ratio of
> TI:TR should be <1. Figure 6 shows that the data are scattered but **the line of best fit
> (R = 0.421) reveals that the TI:TR ratio is in fact >1 until the TR Tmax reaches 0.6 min,
> after which TI:TR is <1**. Although **these results need to be treated with caution, due to
> the scatter and the fact that they represent only one volatile**, they suggest that there
> are **two processes taking place, an initial lag phase where the Tmax for perception occurs
> LATER than the Tmax of the stimulus and then an adaptive phase where the Tmax for perception
> occurs EARLIER than the Tmax of the stimulus.**"*

**⚠️ R = 0.421 ⇒ R² = 0.177 — the fit explains 18 % of the variance.** The crossover at
**TR Tmax = 0.6 min** is read off a line that accounts for less than a fifth of the scatter,
and the authors say so. **Record the 0.6 min crossover as `low_confidence, R2=0.18`; do not
ingest it as a parameter.**

---

## 6. WHAT IS NOT IN THE PAPER

| content | figure | status |
|---|---|---|
| **Absolute sensory scaling means ± SE for the five gels** | **Figure 1** | **NOT digitised this wave.** The normalised values are in Figure 7 (§5.2) and the scale anchors are in the methods (unflavoured 5 % = 0, flavoured 2 % = 10). |
| **Averaged TI profiles for all five gelatine concentrations** | **Figure 2** | **NOT digitised.** Figure 5 gives three of the five as paired TR/TI curves. |
| **In-vitro release curves, 2/5/8 % gels, duplicate, to 45 min** | **Figure 3** | **NOT digitised** — the load-bearing content is the *equality of the plateaux*, which the paper asserts in prose (§4.1) without a tolerance. |
| **TR Imax vs gel FFA concentration, 3 gels × 5 FFA levels × 3 replicates, one panellist** | **Figure 4** | **NOT digitised.** The load-bearing content is *linearity*, asserted in prose (§4.2). No slopes or intercepts are printed. |
| **Individual TI:TR Tmax ratios** | **Figure 6** | **NOT digitised**; only the fitted R = 0.421 and the 0.6 min crossover are printed. |
| **Any absolute t25, t75 or gradient value in physical units** | **nowhere** | ⚠️ **Equation (1) is defined but never evaluated numerically anywhere in the paper.** All gradient values are normalised to the 2 % gel. **No gradient in ppbv/min exists in the source.** |
| **Gel texture, melting point or breakdown measurements** | **nowhere** | ⚠️ The mechanism — *"different rates of gel breakdown in-mouth"* — is inferred, never measured. The paper offers three candidate causes and tests none: *"an increase in melting point"* [C] Wilson & Brown 1997, *"a different rate of dissolution"*, *"a different rate of breakdown due to chewing efficiency"* [C] Brown et al. 1995. |
| **Any FFA odour threshold** | **nowhere** | Not measured, not cited. The gels are dosed 10³–10⁶× above any plausible threshold. |

---

## 7. NAMED LAUNDERING HAZARDS

| # | claim, as printed | reality | anchor |
|---:|---|---|---|
| B-1 | Abstract: **"there was NO BINDING of volatile to protein in the gel, nor to mucous membranes"** | True for **one neutral ester (furfuryl acetate)** in **one nearly-pocketless protein (gelatine)** at **100–1 500 mg/kg**, i.e. 10³–10⁶× above threshold, where any saturable site is saturated (§4.3). **Not a general exclusion of flavour–protein binding**, and Meynier 2002 (this wave) measures 6.88× suppression of an α,β-unsaturated aldehyde by dairy protein on the same class of experiment. | Abstract vs §4.3 |
| B-2 | **"the gradient correlated better with the sensory values than TR Imax"** | The paper computes **no correlation at all** between any instrumental and any sensory variable. The claim rests on visual inspection of Figure 7 plus four separate ANOVAs. (Computed here: **r = 0.968 vs 0.860** against TI Imax — the claim is **correct**, but it was **unsupported as published**.) | p. 159 |
| B-3 | **"the TI:TR ratio is in fact >1 until the TR Tmax reaches 0.6 min"** | Read off a line of best fit with **R = 0.421, i.e. R² = 0.18**. The authors themselves say the results *"need to be treated with caution, due to the scatter"*. **The 0.6 min crossover is not a parameter.** | p. 159 |
| B-4 | **"There were no significant differences in the maximum in-nose volatile concentrations"** used as evidence that concentration does not matter | The comparison is over **three gel concentrations** (2/5/8 %) and the observed fall is **1.5×**. **A non-significant 1.5× is a null result with modest power, not a demonstration of no effect.** The strong evidence is the *contrast* with the simultaneously-measured 3.0× significant fall in TI Imax, not the null on its own. | Abstract, p. 158 |
| B-5 | **"were due to different rates of gel breakdown in-mouth"** | **No texture, melting-point or breakdown measurement was made.** Three candidate mechanisms are named and none is tested; the in-vitro experiment measures **dissolution in stirred water at 37 °C**, not in-mouth breakdown. | Abstract vs §6 |
| B-6 | **Equation (1), `Gradient = Imax/2(t75 − t25)`**, presented as a defined instrumental parameter | **It is never evaluated in physical units anywhere in the paper.** Only 2 %-normalised values exist. A repo that implements eq. (1) has **no numerical value to validate against**. | p. 159 vs §6 |
| B-7 | Reference **"Linforth, R.S.T., Baek, I. and Taylor A.J. (1998) ... Food Chem., in press"** | **Cited three times as support for the two-phase lag/adaptation interpretation, and it was unpublished at the time.** Any downstream chain that rests on that interpretation rests on an in-press citation. | p. 160 |

---

## 8. CONSOLIDATED NEW-PARAMETER TABLE

**Common conditions:** gelatine (250 bloom) gels at 2.0–8.0 % w/w with **20.5 % sucrose +
22.5 % glucose syrup**, flavoured with **furfuryl acetate at 300 mg/kg** (unless stated),
6.0 ± 0.1 g cubes, overnight at 4 °C then equilibrated to room temperature; **API-MS at
m/z 80.8, one nostril, 25–30 mL/min, 2 min**; Nottingham 1998.

| # | parameter | value | units | conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1–3 | **TR Imax (max in-nose FFA)** | **≈475 / ≈332 / ≈313** | **ppbv** | 2 / 5 / 8 % w/w gelatine, 300 mg/kg FFA, 11 panellists | **D ±5 %** | Fig. 5 |
| 4–6 | **TR Tmax** | **≈0.47 / ≈0.57 / ≈0.80–0.85** | min | " | **D** | Fig. 5 |
| 7–9 | **TI Imax (perceived intensity)** | **≈10.2 / ≈6.7 / ≈3.4** | sensory units (0 = unflavoured 5 % gel, 10 = flavoured 2 % gel) | " | **D ±5 %** | Fig. 5 |
| 10–12 | **TI Tmax** | **≈0.55 / ≈0.60 / ≈0.85** | min | " | **D** | Fig. 5 |
| 13 | **TR Imax fall, 2 → 8 % gelatine** | **1.49–1.52× — NOT SIGNIFICANT** | × | ANOVA over 3 gels | M + Z | §5.1, §5.2 |
| 14 | **TI Imax fall, 2 → 8 %** | **3.00–3.28×, P < 0.001** | × | ANOVA | M + Z | §5.1, §5.2 |
| 15 | **sensory scaling score fall, 2 → 8 %** | **2.27×, P < 0.001** | × | 16 panellists | M + Z | §5.2 |
| 16 | **TR gradient fall, 2 → 8 %** | **2.99×, P < 0.01** | × | ANOVA | M + Z | §5.2 |
| 17–36 | **normalised Fig. 7 grid** | 4 variables × 5 gelatine concentrations; see §5.2 table | % of the 2 % gel | 300 mg/kg FFA, 11 panellists | **D ±3 pp** | Fig. 7 |
| 37 | **r(TR gradient, TI Imax)** | **+0.968** | r, n = 5 | paired, simultaneous | **Z** | §5.2 |
| 38 | **r(TR Imax, TI Imax)** | **+0.860** | r, n = 5 | " | **Z** | §5.2 |
| 39 | **r(TR gradient, sensory score)** | **+0.886** | r, n = 5 | different panels (11 vs 16) | **Z** | §5.2 |
| 40 | **r(TR Imax, sensory score)** | **+0.851** | r, n = 5 | " | **Z** | §5.2 |
| 41 | **rate-of-release definition** | **`Gradient = Imax / 2(t75 − t25)`** | — | ⚠️ never evaluated in physical units | M | eq. (1), p. 159 |
| 42 | **in-vitro dissolution time, 2 % gel** | **10** | min | 6 g cube in 100 mL stirred water at 37 °C | M | Fig. 3 / p. 158 |
| 43 | **in-vitro dissolution time, 8 % gel** | **40** | min | " | M | Fig. 3 / p. 158 |
| 44 | **plateau headspace FFA at full dissolution** | **EQUAL for 2, 5 and 8 % gelatine ⇒ NO BINDING** ⚠️ tolerance not stated | — | 20–80 g/L gelatine, 300 mg/kg FFA, 37 °C | M | §4.1 |
| 45 | **in-nose Imax vs gel FFA concentration** | **LINEAR** over 100–1 500 mg/kg, 3 gel concentrations × 3 replicates, replicated in 4 panellists ⇒ **no mucosal binding** ⚠️ no slope or intercept printed | — | one panellist + 3 replicates | M | §4.2 |
| 46 | **inter-panellist variation in released amount** | *"substantial differences between panellists"* — **magnitude not stated** | — | same gels, 4 panellists | M | p. 158 |
| 47 | **TI:TR Tmax crossover** | **0.6** | min | ⚠️ **R = 0.421, R² = 0.18 — low confidence** | F | Fig. 6 / §5.3 |
| 48 | **API-MS in-breath detection range** | **10–100** | ppbv | cited | C | p. 156 |
| 49 | gel FFA dosing matrix | 22 gels; 5 gelatine levels × up to 6 FFA levels | mg/kg | see Table 1 | M | T1 |

---

## 9. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Baek et al. 1999 is a **new source**; a declaration amendment must be
> approved before any wave fits any row. This dossier does not edit the declaration.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **§5.1/§5.2 — the significance contrast (TR Imax n.s. vs TI Imax P < 0.001)** | 2 | **HOLD-OUT — the wave's cleanest falsification test for concentration→intensity mapping** | Proposed guard: **any repo path that maps a predicted headspace/in-nose concentration monotonically onto a perceived intensity is refuted by a within-subject experiment in which the two variables move by 1.5× (n.s.) and 3.0× (P < 0.001) respectively.** This is a stronger form of `k2_matrix_and_thresholds.md` §B.5, because TR and TI were recorded **simultaneously, on the same panellists, on the same mouthful**, rather than compared across experiments. **Use to reject; never fit.** |
| **§4.1 — the binding exclusion (equal plateaux over 20–80 g/L gelatine)** | 1 | **HOLD-OUT — a scoped null** | Ingest **only** as `binding_excluded: {compound_class: neutral_ester, protein: gelatine, concentration: 100-1500 mg/kg}`. §4.3 lists three reasons it does not generalise. It is the *matched control* for Meynier 2002's 6.88× aldehyde suppression, not a contradiction of it. |
| **§4.2 — the mucosal-binding exclusion (linear over 15×, 4 panellists)** | 1 | **HOLD-OUT, and a new constraint on the saliva lane** | Together with `k2_matrix_and_thresholds.md` §B.6 (Starkenmann's 6 × 10⁴ salivary thiol quench), this establishes that **mucosal/salivary interaction is compound-class-specific, not a general multiplicative term.** A repo that applies a single saliva factor to all odourants is refuted by the pair. |
| **The 20 normalised Figure 7 values** | 20 | **PRIOR / SHAPE ONLY — NOT FIT-ELIGIBLE** | Digitised at ±3 pp from a figure, and **normalised** — they carry no absolute scale. The shape (monotone for three variables, **non-monotone for TR Imax**) is the usable content. |
| **The 12 Figure 5 Imax/Tmax values** | 12 | **PRIOR** | Digitised, ±5 %. Cross-validated against Figure 7 to within 3 pp, which is reassuring but does not make them measurements. |
| **The 4 correlation coefficients (§5.2)** | 4 | **DERIVED — report, do not fit** | n = 5 on digitised data. Their value is that they **quantify a claim the paper made without support** (B-2), not that they are precise. |
| **Equation (1), the gradient definition** | 1 | **ADOPT AS A MODEL FORM, with no calibration value** | The definition is clean and the paper shows the derived quantity outperforms Imax. But **no absolute gradient value exists in the source** (B-6), so a repo implementation cannot be numerically validated against this paper. |
| **The 0.6 min TI:TR crossover** | 1 | **REJECT** | R² = 0.18, authors' own caution (B-3). |
| **In-vitro dissolution times (10 and 40 min)** | 2 | **METADATA** | Establishes the 4× rate range that drives the whole result; not a repo parameter (stirred water at 37 °C is not a mouth). |
| **The five printed p-values** | 5 | **METADATA / acceptance gates** | Usable as pass/fail criteria on any repo model of this system. |

---

## 10. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Linforth, R. S. T., Baek, I. & Taylor, A. J. (1998/1999)**, *Food Chem.* —
   *"Simultaneous instrumental and sensory analysis of volatile release from gelatine and
   pectin/gelatine gels."* Cited three times, **in press at the time** (B-7), by the same
   authors, on the same system, **with a second gelling agent (pectin)**. It is the direct
   companion paper and would (a) provide absolute gradient values that this paper lacks,
   (b) test whether the concentration/perception dissociation survives a change of matrix
   polymer. **Highest-value retrieval from this paper.**
2. **Harrison, M. & Hills, B. P. (1996)**, *Int. J. Food Sci. Tech.* **31**, 167–176 —
   *"A mathematical model to describe flavour release from gelatine gels."* **A published
   release-rate model for exactly this matrix**, and the source of the water-soluble-dye
   control that corroborates the dissolution-limited mechanism (§4.1). ⚠️ Note
   `k2_matrix_and_thresholds.md` §"three retrievals" already asks for **Harrison & Hills
   1997** (*JAFC* 45, 1883–1890) for a different reason — the same authors, adjacent years.
   **Requesting both together is efficient.**
3. **Overbosch, P. (1986)**, *Chem. Senses* **11**, 315–329 and **Overbosch & de Jong (1989)**,
   *Physiol. Behav.* **45**, 607–613 — the adaptation model this paper tests and finds
   inadequate. If the repo ever needs a temporal perception model, these are the reference
   objects and this paper is the evidence that they need modification.
4. **Brown, W. E., Dauchel, C. & Wakeling, I. (1995)**, *J. Text. Studies* **27**, 433–450 —
   *"Influences of chewing efficiency on texture and flavour perception of food."* The source
   of the inter-panellist variation this paper observes but does not quantify (row 46).
