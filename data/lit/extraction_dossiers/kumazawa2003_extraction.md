# Kumazawa & Masuda 2003 (10.1021/jf021025f) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/kumazawa2003.pdf` (5 pp.). Read method: **both** — born-digital text
layer (`pdftotext -layout`, clean) for all body text and Tables 1–3; page raster at 900 dpi
(PDF page 4 = printed p. 2677) for Figure 3, which is the **only** place the full seven-point pH grid
of residual ratios appears.

## 0. PAPER IDENTITY — **MATCHES** THE EXPECTED IDENTITY

| field | value |
|---|---|
| authors | Kenji Kumazawa,* Hideki Masuda |
| title | "Investigation of the Change in the Flavor of a Coffee Drink during Heat Processing" |
| venue | *Journal of Agricultural and Food Chemistry* |
| volume/pages/year | Vol. 51, No. 9, pp. 2674–2678, 2003 |
| DOI | 10.1021/jf021025f (article ID printed as `JF021025F`) |
| received / revised / accepted | Received for review **October 7, 2002**; Revised manuscript received **January 11, 2003**; Accepted **January 12, 2003**; Published on Web **March 21, 2003** |
| affiliation | Material Research and Development Laboratories, Ogawa and Company, Ltd., 15-7 Chidori Urayasushi, Chiba 279-0032, Japan |
| PDF character | **Born-digital** (Creator: Parlance Publisher 5.0 / Xyvision PostScript Formatter; Producer: Acrobat Distiller 3.01). Text layer complete; figures are embedded rasters with no numeric data layer. |

Expected identity was "coffee drink heat processing / retort sterilisation at 121 °C, with a pH 3–7
grid" — **confirmed**, with two corrections to the expectation:
1. The pH grid is **3.0, 4.0, 5.0, 5.4, 6.0, 6.4, 7.0** — seven levels, not integer steps; there is an
   extra pair at 5.4 and 6.4. A separate two-row table (Table 2) adds **pH 6.8**, extending the range
   to 3.0–6.8/7.0.
2. The compound set in the pH grid is **only four** compounds (2-furfurylthiol and three of its
   degradation products), not the eight potent odorants of Table 1.

## 1. ONE-PARAGRAPH VERDICT

This paper delivers the requested **pH grid for 2-furfurylthiol thermal survival under retort
conditions**, in two complementary forms: **Figure 3** (residual ratio %, seven pH levels, 121 °C /
10 min — figure-only, no numbers printed in text or table) and **Table 3** (relative GC peak-area
ratios for 2-furfurylthiol, difurfuryl disulfide, furfural and furfuryl alcohol across the same seven
pH levels). **Table 2** adds a small 2 × 2 temperature/time × pH block (pH 6.0 and 6.8 × [121 °C,
10 min] and [123 °C, 20 min]) which is the paper's only temperature/time contrast. **There are NO rate
constants, NO half-lives, NO activation energies, NO D or z values and NO regression of any kind.**
There is also **NO error bar, NO standard deviation, NO replicate count and NO statistical test
anywhere in Tables 2 or 3 or Figure 3** — the entire pH dataset is reported as bare single numbers
with no stated n. Critically for the task: **the retort schedule is given only as "121 °C for 10 min"
(and "123 °C for 20 min"); no F₀ value, no come-up time, no come-down time and no lethality
calculation appear anywhere in the paper.** Likewise **there is NO storage test** — the only mention
of storage is a literature aside about a Thermos flask. The pH-grid experiments are **aqueous
citric/Na₂HPO₄ buffer model solutions at 1 ppm 2-furfurylthiol, not coffee**; the only coffee-drink
number is a single GC-O-derived residual ratio of "∼20 %".

## 2. SYSTEM DEFINITION — verbatim

### 2.1 The coffee drink (used for the sensory / GCO-H work, Table 1 and Figures 1–2)

> "Coffee Drinks. The Arabica coffee beans were medium-roasted (luminosity value = 22.9), sealed under
> vacuum, and stored at −20 °C until used. Deionized hot water (∼90 °C, 500 mL) was poured on the
> ground coffee powders (50 g) in a filter. The filtrate (∼450 mL) was immediately cooled to ∼20 °C in
> tap water, and then NaHCO3 was added to a pH of 6.0; the mixture was canned **without the
> deoxidization process**. The canned coffee drink was sterilized at **121 °C for 10 min** followed by
> immediate cooling to ∼10 °C in tap water. Sterilization was performed in a retort pasteurizer (model
> PRS-02-II-VC; supplied from Nissen, Japan). **The pH value after heat processing of the coffee drink
> was 5.2.**" (p. 2674)

Key facts recorded:
- Brew ratio [Z]: 50 g ground coffee / 500 mL water = **100 g L⁻¹**; recovered filtrate ~450 mL.
- Roast: Arabica, medium, **luminosity value = 22.9** (units/scale not defined by the authors).
- pH adjusted to **6.0 with NaHCO₃ before canning**; **pH after heat processing = 5.2** — a **drop of
  0.8 pH units caused by the retort process itself**. This is a measured, printed number and is one of
  the more useful facts in the paper.
- **Headspace/atmosphere: canned "without the deoxidization process"**, i.e. **air/oxygen was left in
  the can**. Explicitly stated. No inert gas, no vacuum, no deaeration.
- Cooling: "immediate cooling to ∼10 °C in tap water."
- Vessel: can (size/fill volume not stated). Agitation during retort: not stated.

### 2.2 THE RETORT SCHEDULE — everything the paper says

The complete set of retort/heat-processing statements in the paper, verbatim:

> "The canned coffee drink was sterilized at **121 °C for 10 min** followed by immediate cooling to
> ∼10 °C in tap water. Sterilization was performed in a retort pasteurizer (**model PRS-02-II-VC;
> supplied from Nissen, Japan**)." (p. 2674, coffee drink)

> "Half of each solution was canned without the deoxidization process and treated by heat processing
> (**121 °C, 10 min**; the instrument was described above), followed by immediate cooling to ∼10 °C in
> tap water, and the other half of the solution was not heated." (p. 2675, model solutions)

> "...the different temperatures (**121 °C for 10 min and 123 °C for 20 min**) and pH (pH 6.0 and 6.8)
> conditions." (p. 2677)

> Table 2 column headers: "**121 °C, 10 min**" and "**123 °C, 20 min**".

> Figure 3 caption: "Residual ratios of 2-furfurylthiol after model solutions (**121 °C, 10 min**) were
> heated at various pH values."

**NOT PRESENT ANYWHERE IN THE PAPER:**
- **F₀ value** — the string "F0"/"F₀" does not occur. No lethality figure of any kind.
- **Come-up time** — not stated.
- **Come-down / cooling ramp time** — only "immediate cooling to ∼10 °C in tap water".
- **Retort pressure, overpressure, medium (steam vs water-immersion vs steam-air)** — not stated.
- **Whether "10 min" is hold time at 121 °C or total process time** — not disambiguated.
- **z value, D value, or any thermal-death-time parameter** — absent.
- **Storage test** — absent. The only storage reference in the paper is a literature aside:
  "Recent investigations indicated that the melanoidins in the coffee brew are involved in the loss of
  2-furfurylthiol when brewed coffee is kept warm in a Thermos flask" (p. 2674, citing refs 9–11).
  **No storage experiment was run by these authors.**

**Consequence:** the thermal history of these samples is **underdetermined**. A model that needs an
integrated time–temperature exposure cannot get one from this paper; only the nominal
"121 °C / 10 min" and "123 °C / 20 min" labels are available. This must be carried as a systematic.

The 123 °C / 20 min condition is **not an experiment reported in this paper's own Results** — it comes
from the authors' earlier study, ref. 8 (Kumazawa, Masuda, Nishimura & Hiraishi, *Nippon Shokuhin
Kagaku Kogaku Kaishi* 1998, 45, 108–113), which is cited for the 3.6 % residual ratio (see §5.2).
Table 2's right-hand column is therefore **partly [C]** — see the discussion under Table 2.

### 2.3 The aqueous model solutions — the pH grid system (Table 3, Figure 3)

> "Model Experiments on the Thermal Stability of 2-Furfurylthiol in Aqueous Solution. At a
> concentration of **1 ppm**, 2-furfurylthiol was dissolved in **citric/Na2HPO4 buffer solutions of
> pH 3.0, 4.0, 5.0, 5.4, 6.0, 6.4, and 7.0 (mixed 1/10 M citric acid and 1/5 M Na2HPO4)**. Half of each
> solution was canned without the deoxidization process and treated by heat processing (121 °C, 10 min;
> the instrument was described above), followed by immediate cooling to ∼10 °C in tap water, and the
> other half of the solution was not heated. 2-Furfurylthol in the treated and untreated solutions
> (**100 mL**) were isolated by extraction with methylene chloride (50 mL × 2). After drying over
> anhydrous sodium sulfate, the solvent was evaporated to ∼5 mL in volume. A further concentration was
> conducted with a nitrogen stream to ∼100 µL. The internal standard solution (20 µL) prepared from
> methyl octanoate (**53.3 mg**) in methylene chloride (**10 mL**) was added to the concentrate before
> the solvent was removed by the evaporator. Quantification of 2-furfurylthiol was performed by GC, and
> the **residual ratios (percent) were calculated from the measured values of 2-furfurylthiol in the
> thermally treated and nontreated solutions**." (p. 2675)

Key facts:
- **Concentration: 1 ppm 2-furfurylthiol** = 1 mg L⁻¹ = **8.76 µM** [Z] (MW 114.17).
- **Buffer: McIlvaine-type citrate/phosphate**, made by "mixing 1/10 M citric acid and 1/5 M Na₂HPO₄"
  (i.e. **0.1 M citric acid and 0.2 M disodium hydrogen phosphate**). The mixing ratio per pH level is
  **not printed** — it is the standard McIlvaine table, but the authors do not give it, so the ionic
  strength varies across the grid and is **not determinable from the paper**. Flagged.
- **pH is measured at room temperature before heating.** The paper says nothing about hot-pH
  correction, and citrate/phosphate buffers shift with temperature. **The nominal pH labels are
  cold-pH labels; the pH during the 121 °C hold is unknown.** This is a material caveat for any model
  that uses pH as a mechanistic variable. Flagged.
- **Atmosphere: canned "without the deoxidization process"** — air-saturated, oxygen present.
- **Residual ratio definition, verbatim:** "the residual ratios (percent) were calculated from the
  measured values of 2-furfurylthiol in the thermally treated and nontreated solutions." So
  **residual ratio = [heated] / [unheated] × 100**, with the unheated half of the *same* solution as
  the denominator. This is a clean paired design.
- Extraction: CH₂Cl₂, 50 mL × 2 from 100 mL; dried over Na₂SO₄; to ~5 mL; N₂ stream to ~100 µL.
- Internal standard: **methyl octanoate**, 53.3 mg / 10 mL CH₂Cl₂ stock, 20 µL added
  ⇒ **106.6 µg of methyl octanoate per sample** [Z].
- **Replication (n): NOT STATED.** No error bars, no SD, no statistical test for Table 2, Table 3 or
  Figure 3.

### 2.4 Sensory evaluation (Figure 1)

> "The coffee drinks before and after heat processing were placed in glass beakers (∼20 °C). The
> samples were sniffed by **15 panelists (12 males and 3 females, between the ages of 27 and 38)**. All
> panelists had previously received extensive training in descriptive sensory analysis of the coffee
> drinks and had experience in the sensory profiling of various food samples. The intensities of the
> attributes of the heated coffee drink was scored on a **category scale of 1, 2, 3, ..., 7**. For the
> unheated coffee drink, all of the category scales were defined as 4. The results obtained by the
> panelists were then averaged. The resulting data were evaluated by an **analysis of variance**."
> (p. 2675)

Figure 1 caption as printed: "Odor attributes of after-heating coffee drink. For the nonheated coffee
drink, all of the category scales were defined as 4. Significant difference from the control
(nonheated coffee drink): **, p < 0.01; ***, p < 0.001."

**Figure 1 is a spider/bar plot of averaged category scores. The individual attribute scores are not
printed in text or table. The attribute names recoverable from the body text are: pungent, putrid,
heavy (tallowy), sour, caramel-like (all increased) and roasty (significantly decreased).** The
numeric scores themselves are **unreadable** from the text layer; I did not digitise Figure 1 because
category-scale sensory scores are not usable model parameters and the figure is on PDF page 2 which
was not required for the numeric grid.

### 2.5 GCO-H (Table 1, Figure 2)

> "GCO-H was performed with a Hewlett-Packard (HP) model 5890 series II gas chromatograph connected to
> the purge and trap system TCT/PTI (Chrompack). The empty glass tube in the desorption heating block
> of the purge and trap facility was deactivated by treatment with 5% 1,1,1,3,3,3-hexamethyldisilazane,
> which was dissolved in toluene. The gas chromatograph was equipped with a DB-Wax fused silica
> capillary column (30 m × 0.25 mm i.d.; film thickness = 0.25 µm; J&W Scientific) and a thermal
> conductivity detector (TCD). The column temperature was programmed from 40 to 210 °C at the rate of
> 3 °C/min in all runs. The flow rate of the helium carrier gas was 1.5 mL/min. A glass sniffing port
> was connected to the outlet of the TCD and was heated by a ribbon heater. Moist air was pumped into
> the sniffing port. **The coffee drink (10 mL) was placed in a vessel (volume = 100 mL), sealed with a
> septum, and then held in a water bath at 40 °C for 10 min. The headspace volumes [20−0.1 mL
> corresponding to flavor dilution (FD) factors 1 to 200] were drawn by a gastight syringe** and then
> injected into the purge system, which operated in the desorption mode for 5 min at a temperature of
> 150 °C. The carrier gas helium (flow = 10 mL/min) swept the headspace sample into the trap, which was
> cooled with liquid nitrogen at −150 °C. To start the GC run, the cold trap was heated rapidly to
> 200 °C, and this temperature was held for 5 min." (p. 2675)

[Z] The FD scale: headspace volume 20 mL → FD 1; 0.1 mL → FD 200. So **FD factor = 20 mL / (injected
volume in mL)**, a 2-fold (or here, sometimes 5-fold) dilution ladder over 1–200.

### 2.6 Enrichment for identification, GC, GC-MS, headspace GC-MS

> "The coffee drinks (∼1 L, before and after the heat processing) were distilled under reduced pressure
> (**40 °C, 20 mmHg**). The steam distillate (∼250 mL) was passed through a column packed with 10 g of
> Porapak Q (Waters). The adsorbed compounds were eluted with methylene chloride (100 mL). The eluate
> was dried over anhydrous sodium sulfate, and the solvent was removed using a rotary evaporator to
> ∼5 mL in volume. A further concentration was conducted with a nitrogen stream to ∼100 µL." (p. 2675)

> "Gas Chromatography. An HP 5890 series II gas chromatograph equipped with a flame ionization detector
> (FID) was used. A fused silica column (**60 m × 0.25 mm i.d.; coated with a 0.25 µm film of DB-Wax**;
> J&W Scientific) was used. The column temperature was programmed from **80 to 210 °C at the rate of
> 3 °C/min** for all runs. The injector and detector temperatures were both **250 °C**. The flow rate of
> the nitrogen carrier gas was **1.2 mL/min**, and the split ratio was **1:28**." (p. 2675)
> — **This is the instrument that produced Table 3 and Figure 3.** Note it is FID with nitrogen carrier
> and a 1:28 split; no sulfur-selective detector was used for the pH grid.

> "Gas Chromatography−Mass Spectrometry. An HP 5890 series II gas chromatograph coupled to an HP model
> 5972 series mass spectrometer was used. The column was a 60 m × 0.25 mm i.d. DB-Wax fused silica
> capillary column (J&W Scientific) with a film thickness of 0.25 µm. The column temperature was
> programmed from 40 to 210 °C at the rate of 3 °C/min. The injector temperature was 250 °C. The flow
> rate of the helium carrier gas was 1 mL/min, and samples (0.2 µL) were applied using the splitless
> injection technique. The mass spectrometer was used under the following conditions: ionization
> voltage, 70 eV (EI); ion source temperature, 150 °C." (p. 2675)

> "Identification of Components. The identification of the components was made by comparison of their
> Kovats GC retention indices, mass spectra, and odor quality to those of authentic compounds."
> (p. 2675)

### 2.7 Chemicals

> "Methanethiol, 2-furfuryl methyl disulfide, 4-hydroxy-2,5-dimethyl-3(2H)-furanone, and difurfuryl
> disulfide were obtained from Sigma-Aldrich (Tokyo, Japan). 2-Furfurylthiol, methional, acetic acid,
> 3-methylbutanoic acid, furfural, furfuryl alcohol, and methyl octanoate were obtained from Tokyo Kasei
> Kogyo (Tokyo, Japan). 3-Mercapto-3-methylbutyl formate was synthesized according to the literature
> (2)." (pp. 2674–2675) — **no purities are stated for any chemical.**

---

## 3. TABLE 1 — COMPLETE TRANSCRIPTION

**Anchor: Table 1, p. 2676 (PDF page 3).** Title as printed: "Potent Odorants of Coffee Drinks Showing
Differences in Their FD Factors before and after Heat Processing".
Column headers as printed: `no.^a` | `RI^b` | `RI^c` | `compound` | `odor quality^f` | `FD factor^g`
split into `before heating` and `after heating`.
Footnotes as printed: "^a Numbering refers to Figure 2. ^b Retention index on DB-Wax column (30 m ×
0.25 mm i.d.; coated with a 0.25 µm film) observed for GCO-H. ^c Retention index on DB-Wax column
(60 m × 0.25 mm i.d.; coated with a 0.25 µm film) observed for GC-MS. ^d Identified by
headspace-GC-MS. ^e Identified by GC-MS using enrichment of odorants samples. ^f Odor quality assigned
during GCO-H. ^g Relationship between FD factor and headspace volume is given under Materials and
Methods."
Conditions: coffee drink, before vs after **121 °C / 10 min** retort; GCO-H at 40 °C, 10 min
equilibration; `nd` = not detected.

| no. | RI (30 m, GCO-H) | RI (60 m, GC-MS) | compound | odor quality | FD before heating | FD after heating |
|---|---|---|---|---|---|---|
| 1a | 696 | 690 | methanethiol^d | putrid | 1 [M] | 4 [M] |
| 1b | 696 | (blank) | unknown | pungent | 1 [M] | 4 [M] |
| 4 | 887 | (blank) | unknown | roasty | 1 [M] | nd [M] |
| 5 | 908 | (blank) | unknown | pungent | 20 [M] | 40 [M] |
| 10 | 1116 | (blank) | unknown | roasty | 4 [M] | 2 [M] |
| 15 | 1434 | 1435 | 2-furfurylthiol^e | roasty | **200** [M] | **40** [M] |
| 16 | 1453 | 1448 | methional^e | potato-like | 4 [M] | 1 [M] |
| 17 | 1454 | 1447 | acetic acid^e | sour | nd [M] | 2 [M] |
| 21 | 1523 | 1521 | 3-mercapto-3-methylbutyl formate^e | roasty | 20 [M] | nd [M] |
| 27 | 1684 | 1676 | 3-methylbutanoic acid^e | sour | nd [M] | 4 [M] |
| 28 | 1806 | 1813 | 2-furfuryl methyl disulfide^e | meaty | 4 [M] | 20 [M] |
| 34 | 2043 | 2040 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone^e | caramel-like | 40 [M] | **200** [M] |
| 36 | 2384 | (blank) | unknown | sour | nd [M] | 1 [M] |

Blank cells in the `RI^c` column are for the five unidentified peaks, for which no GC-MS RI is given.

**FD factors are ordinal dilution factors, not concentrations.** They must not be treated as
proportional to amount. The authors say so themselves:
> "However, it is known that FD factors in general contain some errors (13), and the static headspace
> sampling technique has the problem that the high water solubility and high boiling point compounds
> suffer losses. Therefore, to reveal the change in the absolute amount of these compounds, which were
> contained in the coffee drinks consisting of complex mixtures of food matrices, a further study would
> be required." (p. 2676)

The GCO-H run detected **34 odor-active peaks before heating and 35 after**, "with flavor dilution (FD)
factors in the range of 1−200" (p. 2676) [M]. Twelve peaks changed; eight of the twelve were
identified.

**Figure 2** (p. 2676, PDF page 3) caption as printed: "Flavor dilution chromatograms (DB-Wax
stationary phase column) of coffee drinks before (top) and after (bottom) heat processing. Numbers
indicate the positions at which an odor was perceived at the sniffing port." — **it is a stick plot of
the Table 1 FD factors plus the 22–23 unchanged peaks; the unchanged peaks' FD values are not
tabulated and are not readable to better than a factor of 2 from the render. Recorded as unreadable;
Table 1 supersedes it for the twelve changed peaks.**

---

## 4. TABLE 2 — COMPLETE TRANSCRIPTION (the temperature/time × pH block)

**Anchor: Table 2, p. 2677 (PDF page 4).** Title as printed: "Residual Ratios^a of 2-Furfurylthiol
after Heating of Model Solutions at Different Temperatures and pH Values".
Column headers as printed: `pH` | `121 °C, 10 min` | `123 °C, 20 min`.
Footnote as printed: "^a Residual ratios are denoted by percentage."
Conditions: aqueous citrate/phosphate model solution, 1 ppm 2-furfurylthiol, canned without
deoxidization. **No error bars, no SD, no n, no significance test.**

| pH | 121 °C, 10 min (%) | 123 °C, 20 min (%) |
|---|---|---|
| 6.8 | 5.5 [M] | 5.3 [M] |
| 6.0 | 48.3 [M] | 43.4 [M] |

**Only four cells.** This is the entire temperature/time axis of the paper.

### 4.1 What Table 2 does and does not support

The authors' reading, verbatim (p. 2677):
> "The results shown in Table 2 indicate that **more pronounced differences in the residual ratios
> occurred at the different pH conditions compared with the different time and temperature
> conditions.**"

[Z] Quantifying that claim from the four cells:
- pH effect at 121 °C/10 min: 48.3 → 5.5 %, a **factor of 8.8** for 0.8 pH units.
- pH effect at 123 °C/20 min: 43.4 → 5.3 %, a **factor of 8.2**.
- Time/temperature effect at pH 6.0: 48.3 → 43.4 %, a **factor of 1.11** for +2 °C and +10 min.
- Time/temperature effect at pH 6.8: 5.5 → 5.3 %, a **factor of 1.04**.

**The authors' claim is strongly supported: the 0.8-pH-unit step is ~8× more consequential than
doubling the hold time and adding 2 °C** [Z].

[Z] Apparent pseudo-first-order constants (the authors compute none):

| condition | residual | apparent k = −ln(residual)/t | apparent t½ |
|---|---|---|---|
| pH 6.0, 121 °C, 10 min | 0.483 | 0.0727 min⁻¹ | 9.5 min |
| pH 6.0, 123 °C, 20 min | 0.434 | 0.0417 min⁻¹ | 16.6 min |
| pH 6.8, 121 °C, 10 min | 0.055 | 0.290 min⁻¹ | 2.4 min |
| pH 6.8, 123 °C, 20 min | 0.053 | 0.147 min⁻¹ | 4.7 min |

**These are internally inconsistent by a factor of ~1.7–2:** doubling the hold time at a *higher*
temperature yields a *lower* apparent rate constant at both pH values. Either the loss is **not
first-order** (most likely — it saturates once the accessible oxidant/oxygen in the can is consumed),
or the nominal "10 min"/"20 min" labels do not represent the true thermal exposure (see §2.2 — no
come-up time is given, so the 10 min process may carry a proportionally larger ramp contribution).
**Do not extract a rate constant or a z value from Table 2.** Recorded as [Z] with that warning.

### 4.2 Provenance caveat on the 123 °C / 20 min column

The paper's Results text ties the 123 °C / 20 min, pH 6.8 condition to the authors' **earlier** study:
> "However, the residual ratio of 2-furfurylthiol (∼20%), which was calculated from the measured values
> of the GC-O experiment before and after heat processing of the coffee drinks, was much higher than
> that of the model experiment (**3.6%: pH 6.8, 123 °C for 20 min**) **in our previous study (8)**."
> (p. 2677)

**Table 2's pH 6.8 / 123 °C-20 min cell reads 5.3 %, while the text attributes 3.6 % to the same
nominal condition from ref. 8.** The paper never reconciles the two. The most consistent reading is
that Table 2 reports **new** measurements at all four cells (which the Materials & Methods section
"the following model experiment was done to clarify the residual ratios of 2-furfurylthiol affected by
the different conditions" supports), and that 3.6 % is a **[C]** value from ref. 8 for a nominally
identical condition. **Both numbers recorded; the 5.3 % is [M] and the 3.6 % is [C]; the ~1.5×
disagreement between two runs of the same nominal condition is a direct measure of this paper's
irreproducibility and should inform σ.**

---

## 5. TABLE 3 — COMPLETE TRANSCRIPTION (the full pH grid, all four compounds)

**Anchor: Table 3, p. 2677 (PDF page 4).** Title as printed: "Quantitative Data (Relative Amounts)^a
for 2-Furfurylthiol and Degradation Products after Heating of Model Solutions at Different pH Values".
Column headers as printed: `compound` under a spanning header `pH` with columns
`3.0  4.0  5.0  5.4  6.0  6.4  7.0`.
Footnote as printed: "^a Relative amounts of compounds are shown as **the ratio of the peak area of
each compound to the peak area of an internal standard (methyl octanoate) on the GC**."
Conditions: **after heating at 121 °C, 10 min**; 1 ppm initial 2-furfurylthiol; citrate/phosphate
buffer; canned without deoxidization; GC-FID, 60 m DB-Wax. **No error bars, no SD, no n.**
Units: **dimensionless peak-area ratio to methyl octanoate**, NOT concentration and NOT µg/L.

| compound | pH 3.0 | pH 4.0 | pH 5.0 | pH 5.4 | pH 6.0 | pH 6.4 | pH 7.0 |
|---|---|---|---|---|---|---|---|
| 2-furfurylthiol | 0.46 [M] | 0.44 [M] | 0.41 [M] | 0.37 [M] | 0.21 [M] | 0.05 [M] | <0.01 [M] |
| difurfuryl disulfide | <0.01 [M] | 0.01 [M] | 0.03 [M] | 0.07 [M] | 0.13 [M] | 0.14 [M] | 0.17 [M] |
| furfural | <0.01 [M] | 0.01 [M] | 0.01 [M] | 0.01 [M] | 0.02 [M] | 0.02 [M] | 0.03 [M] |
| furfuryl alcohol | <0.01 [M] | <0.01 [M] | <0.01 [M] | <0.01 [M] | <0.01 [M] | <0.01 [M] | 0.01 [M] |

All `<0.01` entries are **left-censored at the reporting resolution of 0.01**, not zeros.

### 5.1 [Z] Derived quantities from Table 3

**(a) Product/precursor ratios and the disulfide branching.** Two 2-furfurylthiol molecules make one
difurfuryl disulfide. On a *peak-area-ratio* basis (which is not a molar basis — see the warning
below):

| pH | FFT remaining | disulfide formed | disulfide / FFT [Z] | disulfide as fraction of (FFT + 2×disulfide) [Z] |
|---|---|---|---|---|
| 3.0 | 0.46 | <0.01 | <0.022 | <4.2 % |
| 4.0 | 0.44 | 0.01 | 0.023 | 4.3 % |
| 5.0 | 0.41 | 0.03 | 0.073 | 12.8 % |
| 5.4 | 0.37 | 0.07 | 0.189 | 27.5 % |
| 6.0 | 0.21 | 0.13 | 0.619 | 55.3 % |
| 6.4 | 0.05 | 0.14 | 2.80 | 84.8 % |
| 7.0 | <0.01 | 0.17 | >17 | >97.1 % |

**Warning:** the "peak area ratio to methyl octanoate" is a **detector-response-uncorrected** quantity.
The authors give **no FID response factors** for any of the four compounds. A carbon-count (effective
carbon number) correction would change these ratios by tens of percent. **Treat the column as a
within-compound pH series, never as a cross-compound mass balance.** [Z] with that caveat.

**(b) Mass-balance shortfall — the authors' own observation, verbatim (p. 2677):**
> "However, **the total amount of the volatile degradation products of 2-furfurylthiol detected by GC
> analysis was less than the amount of 2-furfurylthiol lost during heat processing.** It can be assumed
> that the volatile and nonvolatile degradation products were produced at the same time according to a
> Fenton-type reaction (11)."

[Z] Checking that with a naive carbon-blind area balance, taking the pH 3.0 row as the closest thing to
an unreacted reference (FFT = 0.46): at pH 7.0, FFT lost = 0.46 − <0.01 ≈ 0.45 area units, while
products formed = 0.17 (disulfide) + 0.03 (furfural) + 0.01 (furfuryl alcohol) = 0.21 area units.
Even without response correction, and even crediting the disulfide with two thiol equivalents
(0.34 + 0.03 + 0.01 = 0.38), **the accounted fraction is at most ~84 % and is likely well under that
once response factors are applied.** The authors' qualitative statement is supported; a nonvolatile
sink exists. Ref. 11 is Blank et al. 2002 (the Fenton paper extracted in this same wave) — the two
papers cross-reference each other on exactly this point.

**(c) Reconstructing an initial value.** Table 3 reports only **heated** solutions. The unheated
half-samples were measured (they are the denominator of the residual ratios in Fig. 3 / Table 2) but
their peak-area ratios are **never printed**. Combining Table 3 with Fig. 3:
at pH 3.0, residual ratio ≈ 99 % [fig] with FFT area ratio 0.46 [M] ⇒ **unheated FFT area ratio ≈ 0.465**
[Z]; at pH 6.0, residual 45 % [fig] with area 0.21 ⇒ unheated ≈ 0.467 [Z]; at pH 6.4, residual 11 %
[fig] with area 0.05 ⇒ unheated ≈ 0.45 [Z]. **The three independent reconstructions agree to within
4 %, which is a strong internal-consistency check across Table 3 and Figure 3** ✓ and implies a common
unheated area ratio of **≈0.46 ± 0.01** [Z] for 1 ppm 2-furfurylthiol against 106.6 µg methyl octanoate.

At pH 5.4 the check is looser: residual 79.5 % [fig] × 0.465 = 0.370, table 0.37 ✓ exact.
At pH 5.0: 89.1 % × 0.465 = 0.414, table 0.41 ✓.
At pH 4.0: 96.2 % × 0.465 = 0.447, table 0.44 ✓.
At pH 7.0: 0.1 % × 0.465 = 0.0005, table <0.01 ✓ (consistent with censoring).
**Figure 3 and Table 3 are mutually consistent at all seven pH levels** — the strongest quality signal
in this paper.

### 5.2 The authors' text on the pH grid, verbatim (p. 2677)

> "Figure 3 shows the residual ratios of 2-furfurylthiol in aqueous model solutions at various pH
> values. The data clearly indicate that the stability of 2-furfurylthiol during heat processing has a
> close correlation with the pH of the aqueous solutions. The residual ratios of 2-furfurylthiol in the
> after-heating model solutions decreased with increasing pH. **Especially, in the pH range of 5.0−7.0,
> the residual ratios sharply decreased with an increase in the pH.** This finding can account for the
> residual ratios of the 2-furfurylthiol difference between pH 6.0 and 6.8. In the after-heating model
> solutions, difurfuryl disulfide, furfural, and furfuryl alcohol were detected by GC as the volatile
> degradation products of 2-furfurylthiol. **Among these compounds, difurfuryl disulfide was the major
> degradation product, and its amount in the after-heating model solutions increased with increasing
> pH (Table 3).** It is well-known that 2-furfurylthiol easily oxidizes to disulfide (14, 15). The
> oxidation mechanism of thiols may be either polar, radical, or both (16)."

> "In general, the canned coffee drinks are sterilized after the pH range has been adjusted to between
> 5 and 7 to prevent an increase in the sour taste and the cohesion of the coffee and milk components
> (12). However, as for the concentration of 2-furfurylthiol in the after-heating model solutions, it
> was found that the significant difference were brought about by the condition of pH (in the range of
> 5−7) of aqueous solutions. **Therefore, the sulfury−roasty flavor of canned coffee drinks would be
> expected to change with the difference in only a few pH values during heat processing.**"

> "The concentration of 2-furfurylthiol in the coffee drinks was much lower than in the model solutions
> in these experiments; moreover, the coffee drinks consist of complex mixtures of food matrices.
> Therefore, it is thought that **the change in the flavor of the coffee drink during heat processing is
> associated with the various reactions [for example, thiols binding to coffee melanoidins (9, 10) and
> the Fenton-type reaction (11)], together with pH-dependent degradation.**" (p. 2678)

The last quote is the paper's explicit acknowledgement that the **model-solution pH curve is not
directly transferable to coffee** — the coffee matrix adds melanoidin binding (refs 9, 10 = Hofmann &
Schieberle) and Fenton chemistry (ref 11 = Blank et al. 2002). Important for hold-out design.

**The 1 ppm model / coffee mismatch, verbatim (p. 2677 and 2678):**
> "the residual ratio of 2-furfurylthiol (**∼20%**), which was calculated from the measured values of
> the GC-O experiment before and after heat processing of the coffee drinks, was much higher than that
> of the model experiment (3.6%: pH 6.8, 123 °C for 20 min) in our previous study (8)."
> "The concentration of 2-furfurylthiol in the coffee drinks was much lower than in the model solutions
> in these experiments."

**Note the "∼20 %" for the real coffee drink is derived by the authors from the FD-factor drop
200 → 40 in Table 1 (40/200 = 20 %).** It is therefore an **ordinal-dilution-based estimate, not a
quantitative measurement** — [Z]-by-authors at best, and the authors themselves warn (p. 2676) that
"FD factors in general contain some errors". **Do not use ∼20 % as a measured residual ratio.**

---

## 6. FIGURE 3 — DIGITISED TRANSCRIPTION (the full seven-point pH grid)

**Anchor: Figure 3, p. 2677 (PDF page 4).** Caption as printed: "Residual ratios of 2-furfurylthiol
after model solutions (121 °C, 10 min) were heated at various pH values."
Axes as printed: y = `Residual ratio (%)`, ticks 0, 10, 20, …, 100; x = `pH`, ticks 2.5, 3.0, 3.5, 4.0,
4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5 (linear).
Symbols: filled circles joined by straight line segments. **No error bars are drawn.** No n is stated.
Read method: 900 dpi raster of PDF p. 4, cropped to the figure (`k4a/kuma_fig3-4.png`), linear
calibration on the printed y ticks (100 % at y = 83 px, 0 % at y = 1592 px in the displayed render).
Digitisation uncertainty **±1.5 percentage points**.

| pH | residual ratio /% | provenance |
|---|---|---|
| 3.0 | 99.5 | [fig] |
| 4.0 | 96.2 | [fig] |
| 5.0 | 89.1 | [fig] |
| 5.4 | 79.5 | [fig] |
| 6.0 | 45.1 | [fig] |
| 6.4 | 11.3 | [fig] |
| 7.0 | 0.1 | [fig] |

The seven plotted pH values are exactly the seven buffer levels named in the Methods (3.0, 4.0, 5.0,
5.4, 6.0, 6.4, 7.0) ✓.

### 6.1 CONFLICT between Figure 3 and Table 2 at pH 6.0 (flagged, both quoted)

Both describe **2-furfurylthiol residual ratio at pH 6.0 after 121 °C / 10 min in the same model
solution**:
- **Table 2 prints 48.3 %** [M].
- **Figure 3 digitises to 45.1 %** [fig], using a calibration that reproduces the pH 3.0 point at
  99.5 % (i.e. the top of the axis) and the pH 7.0 point at ~0.

The gap is **3.2 percentage points, ~7 % relative** — larger than my ±1.5 pp digitisation uncertainty.
I cannot resolve it: either the figure marker is plotted slightly low, or the two are separate runs.
**Recommendation: prefer Table 2's 48.3 % [M] at pH 6.0 and use the Figure 3 values only for the six
pH levels Table 2 does not cover.** Do not average them. The two numbers are also reconcilable through
Table 3 (§5.1c) only to ±4 %, which is consistent with the gap being real run-to-run scatter.

Note also that **pH 6.8 appears in Table 2 (5.5 %) but not in Figure 3**, and **pH 3.0, 4.0, 5.0, 5.4,
6.4 and 7.0 appear in Figure 3 but not in Table 2**. The two are complementary, overlapping in exactly
one cell — the one that disagrees.

### 6.2 [Z] Shape of the pH dependence

| pH interval | Δ(residual)/ΔpH | provenance |
|---|---|---|
| 3.0 → 4.0 | −3.3 %/pH unit | [Z] from Fig. 3 |
| 4.0 → 5.0 | −7.1 %/pH unit | [Z] |
| 5.0 → 5.4 | −24 %/pH unit | [Z] |
| 5.4 → 6.0 | **−57 %/pH unit** | [Z] |
| 6.0 → 6.4 | **−85 %/pH unit** | [Z] |
| 6.4 → 7.0 | −19 %/pH unit (floor reached) | [Z] |

[Z] Apparent pseudo-first-order rate constant across the grid (10 min nominal hold):

| pH | residual | k = −ln(residual)/10 min /min⁻¹ [Z] | log₁₀ k [Z] |
|---|---|---|---|
| 3.0 | 0.995 | 5.0 × 10⁻⁴ | −3.30 |
| 4.0 | 0.962 | 3.9 × 10⁻³ | −2.41 |
| 5.0 | 0.891 | 1.15 × 10⁻² | −1.94 |
| 5.4 | 0.795 | 2.29 × 10⁻² | −1.64 |
| 6.0 | 0.451 (0.483 per Table 2) | 7.96 × 10⁻² (7.27 × 10⁻²) | −1.10 (−1.14) |
| 6.4 | 0.113 | 2.18 × 10⁻¹ | −0.66 |
| 6.8 (Table 2) | 0.055 | 2.90 × 10⁻¹ | −0.54 |
| 7.0 | 0.001 | ≥ 6.9 × 10⁻¹ (censored) | ≥ −0.16 |

**A log₁₀k-vs-pH regression over pH 3.0–6.8 gives a slope of ≈0.71 decades per pH unit** [Z]
(least squares on the seven uncensored points, excluding the censored pH 7.0). That is **sub-unit
slope**, i.e. **not first-order in [OH⁻] and not a clean thiolate-fraction dependence either** (a pure
thiolate mechanism with pKa(2-furfurylthiol) ≈ 9–10 would give a slope of ~1.0 well below the pKa).
The sub-unit slope is the single most model-relevant [Z] result available from this paper. **It is
mine, not the authors' — they print no regression of any kind.** Treat with care: it rests on a
pseudo-first-order assumption that §4.1 already shows to be internally inconsistent across time.

---

## 7. WHAT THE PAPER DOES NOT CONTAIN (explicit negative findings for the orchestrator)

| requested / expected item | status |
|---|---|
| F₀ value | **ABSENT** — no lethality value of any kind is printed |
| come-up time | **ABSENT** |
| come-down / cooling ramp time | **ABSENT** (only "immediate cooling to ∼10 °C in tap water") |
| hold time | present, but only as the nominal "10 min" (and "20 min"); not disambiguated from total process time |
| retort medium / pressure / rotation | **ABSENT** |
| storage test / shelf-life data | **ABSENT** — no storage experiment was performed |
| replicate count (n) for Tables 2, 3 or Fig. 3 | **ABSENT** |
| error bars / SD / CV on any pH-grid number | **ABSENT** |
| any statistical test on Tables 2 or 3 | **ABSENT** (ANOVA is used only for the Fig. 1 sensory data) |
| rate constants, half-lives, D or z values, Ea | **ABSENT** |
| FID response factors for the four Table 3 compounds | **ABSENT** |
| McIlvaine mixing ratios / ionic strength per pH level | **ABSENT** |
| hot-pH correction | **ABSENT** — nominal pH values are cold-measured |
| absolute concentrations for Table 3 | **ABSENT** — only peak-area ratios |
| unheated (t = 0) peak-area ratios | **ABSENT** — measured but never printed |

---

## NEW-PARAMETER TABLE (consolidated)

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| 2-FFT residual ratio, pH 6.8 | 5.5 | % | 121 °C, 10 min, citrate/phosphate, 1 ppm | Table 2, p. 2677 | [M] |
| 2-FFT residual ratio, pH 6.0 | 48.3 | % | 121 °C, 10 min | Table 2, p. 2677 | [M] |
| 2-FFT residual ratio, pH 6.8 | 5.3 | % | 123 °C, 20 min | Table 2, p. 2677 | [M] |
| 2-FFT residual ratio, pH 6.0 | 43.4 | % | 123 °C, 20 min | Table 2, p. 2677 | [M] |
| 2-FFT residual ratio, pH 3.0 | 99.5 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] (900 dpi crop, PDF p. 4) |
| 2-FFT residual ratio, pH 4.0 | 96.2 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] |
| 2-FFT residual ratio, pH 5.0 | 89.1 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] |
| 2-FFT residual ratio, pH 5.4 | 79.5 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] |
| 2-FFT residual ratio, pH 6.0 | 45.1 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] — **conflicts with Table 2's 48.3** |
| 2-FFT residual ratio, pH 6.4 | 11.3 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] |
| 2-FFT residual ratio, pH 7.0 | 0.1 | % | 121 °C, 10 min | Fig. 3, p. 2677 | [fig] |
| 2-furfurylthiol relative amount, pH 3.0–7.0 | 0.46 / 0.44 / 0.41 / 0.37 / 0.21 / 0.05 / <0.01 | peak-area ratio to methyl octanoate | after 121 °C, 10 min | Table 3, p. 2677 | [M] |
| difurfuryl disulfide relative amount, pH 3.0–7.0 | <0.01 / 0.01 / 0.03 / 0.07 / 0.13 / 0.14 / 0.17 | peak-area ratio to methyl octanoate | as above | Table 3, p. 2677 | [M] |
| furfural relative amount, pH 3.0–7.0 | <0.01 / 0.01 / 0.01 / 0.01 / 0.02 / 0.02 / 0.03 | peak-area ratio to methyl octanoate | as above | Table 3, p. 2677 | [M] |
| furfuryl alcohol relative amount, pH 3.0–7.0 | <0.01 ×6 / 0.01 at pH 7.0 | peak-area ratio to methyl octanoate | as above | Table 3, p. 2677 | [M] |
| unheated 2-FFT peak-area ratio (reconstructed) | 0.46 ± 0.01 | peak-area ratio | 1 ppm 2-FFT vs 106.6 µg methyl octanoate | derived from Table 3 + Fig. 3 | **[Z]** |
| apparent pseudo-1st-order k vs pH (10 min) | 5.0e−4 / 3.9e−3 / 1.15e−2 / 2.29e−2 / 7.96e−2 / 2.18e−1 / 2.90e−1 at pH 3.0/4.0/5.0/5.4/6.0/6.4/6.8 | min⁻¹ | 121 °C | derived from Fig. 3 + Table 2 | **[Z]** |
| slope of log₁₀ k vs pH, pH 3.0–6.8 | ≈ **0.71** | decades per pH unit | 121 °C, 10 min | derived from Fig. 3 + Table 2 | **[Z]** — sub-unit slope |
| pH sensitivity vs time/temperature sensitivity | 8.8× (pH 6.0→6.8) vs 1.11× (10→20 min, 121→123 °C) | ratio of residual ratios | model solution | derived from Table 2 | **[Z]** |
| disulfide / 2-FFT area ratio vs pH | <0.022 / 0.023 / 0.073 / 0.189 / 0.619 / 2.80 / >17 | dimensionless | 121 °C, 10 min | derived from Table 3 | **[Z]** |
| pH drop caused by retort processing (coffee drink) | 6.0 → **5.2** (Δ = −0.8) | pH units | canned Arabica brew, 121 °C, 10 min | p. 2674 | [M] |
| model-solution 2-FFT concentration | 1 (8.76 µM) | ppm | citrate/phosphate buffer | p. 2675 | [M] (µM is [Z]) |
| buffer composition | 1/10 M citric acid + 1/5 M Na₂HPO₄ (mixing ratio per pH **not printed**) | M | pH 3.0–7.0 | p. 2675 | [M] |
| brew ratio, coffee drink | 100 (50 g / 500 mL) | g L⁻¹ | ~90 °C pour-over, filter | p. 2674 | [Z] |
| roast luminosity value | 22.9 | (scale not defined) | Arabica, medium roast | p. 2674 | [M] |
| retort schedule | 121 °C, 10 min (and 123 °C, 20 min) — **no F₀, no come-up time** | °C, min | Nissen PRS-02-II-VC retort pasteurizer | pp. 2674–2675 | [M] |
| headspace of cans | air — canned "without the deoxidization process" | — | both coffee drink and model solutions | pp. 2674–2675 | [M] |
| FD factor, 2-furfurylthiol, before → after heating | 200 → 40 | FD factor | coffee drink, 121 °C 10 min, GCO-H at 40 °C | Table 1, p. 2676 | [M] |
| FD factor, methional, before → after | 4 → 1 | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, 3-mercapto-3-methylbutyl formate | 20 → nd | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, methanethiol | 1 → 4 | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, 2-furfuryl methyl disulfide | 4 → 20 | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 40 → 200 | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, acetic acid | nd → 2 | FD factor | as above | Table 1, p. 2676 | [M] |
| FD factor, 3-methylbutanoic acid | nd → 4 | FD factor | as above | Table 1, p. 2676 | [M] |
| retention indices (DB-Wax) for the eight identified odorants | see Table 1 §3 | RI | DB-Wax 30 m and 60 m | Table 1, p. 2676 | [M] |
| odor-active peaks detected, before / after | 34 / 35 | count | GCO-H, FD 1–200 | p. 2676 | [M] |
| 2-FFT residual ratio in the real coffee drink | ∼20 | % | 121 °C, 10 min, canned coffee | p. 2677 | **[Z]-by-authors** — computed from the FD 200→40 drop, not measured |
| 2-FFT residual ratio, previous study | 3.6 | % | pH 6.8, 123 °C, 20 min | p. 2677, ref. 8 | **[C]** — conflicts with Table 2's 5.3 % for the same nominal condition |
| sensory panel | 15 panelists (12 M, 3 F, ages 27–38), 7-point category scale, control fixed at 4 | — | 20 °C, glass beakers | p. 2675 | [M] |

**NO rate constants, NO half-lives, NO Ea, NO D/z values, NO R², NO fitted model, NO error bars and NO
replicate count are printed anywhere in this paper.** Every k, slope and ratio above is tagged [Z] and
is mine.

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Figure 3, acidic arm: pH 3.0, 4.0, 5.0, 5.4** (p. 2677) | **FIT** | pH | Four points on the shallow side of the curve, 121 °C / 10 min, one buffer family, one concentration. Enough to pin the low-pH baseline and the onset of the transition. |
| **Figure 3, alkaline arm: pH 6.0, 6.4, 7.0 + Table 2 pH 6.8** (p. 2677) | **HOLD-OUT** | pH | The steep side, where the residual ratio falls 79.5 % → 0.1 % over 1.6 pH units. A model fitted on pH ≤ 5.4 that can predict the *position and steepness of the collapse* has learned the mechanism; one that cannot will fail unmistakably. This is a clean, well-separated pH cut and the best hold-out the paper offers. **Use Table 2's 48.3 % at pH 6.0, not Fig. 3's 45.1 % (see §6.1); pH 7.0 is left-censored at ~0 and should be scored as ≤1 %.** |
| **Table 3, 2-furfurylthiol row, all seven pH** (p. 2677) | **EXCLUDE from fitting — use as an internal cross-check on Fig. 3 only** | — | It is arithmetically the same measurement as Fig. 3 (§5.1c shows they reconstruct to a common unheated ratio of 0.46 ± 0.01). Fitting both double-counts one experiment. |
| **Table 3, difurfuryl disulfide row, all seven pH** (p. 2677) | **FIT** (disulfide-branching term only) | pH | An orthogonal *product* observable to the *precursor* observable of Fig. 3, from the same runs. It constrains branching, not loss rate. **Uncalibrated FID area ratios — fit the shape, never the absolute level, and only if the model can be reduced to a relative branching fraction.** |
| **Table 3, furfural and furfuryl alcohol rows** (p. 2677) | **EXCLUDE** | — | Furfural spans 0.01–0.03 and furfuryl alcohol is `<0.01` at six of seven pH values. Both are at or near the 0.01 reporting resolution across the whole grid; there is no signal to fit. |
| **Table 2, the 123 °C / 20 min column (pH 6.0 and 6.8)** (p. 2677) | **HOLD-OUT** | time × temperature | The only time/temperature contrast in the paper: two cells, orthogonal to the whole pH-only fit. **Severe caveat: §4.1 shows the apparent first-order k *falls* when time and temperature both increase, so this arm will falsify any first-order model by construction. That is arguably the point — it is diagnostic of saturation/oxidant limitation. Score it as a qualitative direction test (residual ratio must decrease, and only slightly), not as a quantitative target.** |
| **Table 1, FD factors** (p. 2676) | **EXCLUDE from any fit** | — | Ordinal 2- and 5-fold dilution factors on a 1–200 ladder, in a matrix (real coffee) that the authors themselves say is not comparable to the model solutions. The paper's own "∼20 %" coffee residual ratio is derived from these and is not a measurement. |
| **The coffee-drink pH drop 6.0 → 5.2** (p. 2674) | **HOLD-OUT** (if the model predicts brew pH evolution) | — | A single [M] number, but a real and cleanly stated one: the retort process itself acidifies the brew by 0.8 pH units. If the model has an acid-formation term (acetic and 3-methylbutanoic acid both rose in Table 1), this is a direct, otherwise-unconstrained check. Otherwise EXCLUDE. |

**Circularity and structural flags for the orchestrator:**
1. **Figure 3 and the 2-furfurylthiol row of Table 3 are the same experiment.** Fit at most one.
2. **Do not fit both arms of the Fig. 3 curve.** The whole value of this source is that the pH axis
   splits cleanly into a flat arm (≤5.4) and a collapse arm (≥6.0). Fitting both leaves no pH
   extrapolation anywhere in the source.
3. **The pH 6.0 cell is contested** (Table 2 = 48.3 % [M] vs Fig. 3 = 45.1 % [fig]). Pick Table 2, say
   so, and do not average.
4. **The 123 °C / 20 min pH-6.8 cell has a competing published value** (5.3 % here [M] vs 3.6 % in
   ref. 8 [C]). That ~1.5× spread between two runs of one nominal condition is the only direct
   estimate of this paper's reproducibility and should set a floor on σ for all its numbers — there
   are no error bars anywhere to use instead.
5. **These are 1 ppm aqueous buffer solutions, not coffee.** The authors state explicitly (p. 2678)
   that the coffee result "is associated with the various reactions [thiols binding to coffee
   melanoidins and the Fenton-type reaction], together with pH-dependent degradation". Do not use this
   pH curve to predict coffee-matrix thiol survival without a matrix term; if it is used that way,
   record the matrix gap as a systematic.
6. **The pH labels are cold-measured and uncorrected for the 121 °C hold.** Citrate/phosphate buffers
   shift substantially with temperature. If the model uses pH mechanistically, the *effective* hot pH
   is an unquantified offset from the axis labels used here.
7. **The thermal history is underdetermined** (no F₀, no come-up time). Any integrated-lethality
   treatment must invent the missing ramp, and that invention is a systematic, not a σ.
8. **n is unstated and no error bars exist** for any pH-grid number. Do not assume n = 3; do not assign
   a σ derived from anything in this paper other than the reproducibility gap in flag 4.
