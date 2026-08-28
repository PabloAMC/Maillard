# Stack, Conrad & Mahmud 2018 (10.1021/acs.chemrestox.7b00239) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/stack2018.pdf` (7 pp.). Born-digital ACS PDF, clean text layer.
Read method: **both** — text layer for the running text, **plus** 400 dpi page-5 raster
(`s_p5-5.png`, crops `s_fig4.png`, `s_fig4CD.png`, `s_fig4D.png`) for **Figure 4C and 4D, which
are images with no text layer and which carry every rate constant and thermodynamic parameter
in the paper**. Nothing quantitative here is in the text layer except the abstract summary.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY

| field | value |
|---|---|
| Authors | **Douglas E. Stack\*, John A. Conrad, Bejan Mahmud** |
| Title | **"Structural Identification and Kinetic Analysis of the in Vitro Products Formed by Reaction of Bisphenol A-3,4-quinone with N-Acetylcysteine and Glutathione"** |
| Venue | *Chemical Research in Toxicology* **2018**, 31 (1), **81−87** |
| DOI | **10.1021/acs.chemrestox.7b00239** |
| Dates | Received August 25, 2017; **Published December 27, 2017** (© 2017 ACS) |
| Licence | **Open access, CC-BY** |
| Affiliation | Department of Chemistry, University of Nebraska at Omaha, 6001 Dodge Street, Durham Science Center, Omaha, Nebraska 68182 |
| Funding | Funds for Undergraduate Scholarly Experiences (FUSE), University of Nebraska at Omaha |
| SI | Supporting Information (PDF) exists — chromatography conditions, 1-D/2-D NMR, **and the Arrhenius plots used to calculate the free energies of activation (Figures S3, S4)**. **NOT ON DISK.** |

**Correct file for its expected identity.**

---

## 1. ONE-PARAGRAPH VERDICT

This is a **stopped-flow, temperature-resolved, fully reversible thiol + o-quinone Michael
addition** measured at pH 6.2 — exactly the reaction class the repo's thiol-consumption /
quinone-sink lane needs, and the **only dataset in the corpus that supplies a reverse rate
constant for a thiol–quinone adduct**. It gives **14 primary rate constants** (Figure 4C: NAC
forward and reverse at four temperatures, GSH forward and reverse at three) plus a complete set
of derived activation parameters (Figure 4D). **The rate constants are sound and should be
used. The activation parameters should not.** §6 shows that the paper's printed ΔH‡ values are
**systematically 2.303× larger** than an Arrhenius refit of the paper's own Figure 4C
constants — the signature of a natural-log Arrhenius plot processed as if it were a base-10
plot — and that the printed ΔG‡ values regenerate the measured rate constants **126× to 5870×
too fast** under the Eyring equation, so they are neither Arrhenius activation energies nor
Gibbs free energies of activation. **This resolves the ΔG-versus-k₂ convention mismatch the
brief flagged: the two quantities were never on the same footing, and the k values are the
survivable half.** Caveats that limit transfer to a food matrix: the solvent is **50 %
methanol**, the temperature range is **4.6–19.4 °C** (nowhere near any thermal-process
temperature), the electrophile is **bisphenol A-3,4-quinone**, not a coffee/food quinone, and
the thiols are **NAC and GSH**, not 2-furfurylthiol.

---

## 2. SYSTEM DEFINITION — verbatim

### 2.1 Stopped-flow kinetics — the block that matters (Experimental Procedures, p. 83)

> "**Stopped-flow Kinetic Analysis.** A **125 µM BPAQ solution (final concentration) in 50 %
> methanol** was mixed with various concentrations of N-acetylcysteine (**0.050, 0.10, and
> 0.125 M**; final concentration) in **0.050 M phosphate buffer pH 6.2** and **50 % methanol**
> using a stopped-flow spectrophotometer (TgK Scientific). Concentration dependence was carried
> out at **5 different temperatures: 4.6 °C, 9.6 °C, 14.5 °C, and 19.4 °C**. Reactions were
> monitored at **400 nm**. Reaction traces were fit to a **sum of three exponentials** using
> KaleidaGraph (Synergy Software). GSH reactions under the same conditions with the following
> changes: **BPAQ concentration 65 µM; GSH concentrations 40, 20, and 10 mM**."

⚠️ **DEFECT: "5 different temperatures" is followed by a list of FOUR.** The Results section
(p. 84) and Figure 4B both say **four** ("Each NAC concentration was carried out at four
temperatures"; Fig. 4B caption: "run at four different temperatures"), and Figure 4C tabulates
exactly four NAC rows. **Read "5" as a typo for 4.** No fifth temperature exists anywhere in
the paper.

| variable | NAC arm | GSH arm |
|---|---|---|
| Electrophile | **BPAQ, 125 µM final** | **BPAQ, 65 µM final** |
| Nucleophile concentrations | **0.050, 0.10, 0.125 M** (50, 100, 125 mM) | **10, 20, 40 mM** |
| Buffer | **0.050 M phosphate, pH 6.2** | same |
| Solvent | **50 % methanol** (aqueous/methanol 1:1) | same; > "more dilute solutions due to the **reduced solubility of GSH in 50 % methanol** when compared to NAC" |
| Temperatures | **4.6, 9.6, 14.5, 19.4 °C** | **9.6, 14.5, 19.4 °C** — > "Only three different temperatures were conducted with GSH as **insolubility limited lower temperatures**" |
| Probe | **absorbance at 400 nm** (BPAQ disappearance; BPAQ is deep red) | same |
| Instrument | TgK Scientific **KinetAsyst** stopped-flow | same |
| Kinetic regime | **pseudo-first order**, large excess of thiol | same |
| Trace model | **sum of three exponentials**, fit in KaleidaGraph | same |
| pH correction for methanol | **NOT MENTIONED.** "pH 6.2" is the aqueous buffer's nominal pH; the apparent pH of a 50 % methanol mixture is not the same quantity and no correction is reported | — |
| Ionic strength | not reported beyond "0.050 M phosphate" | — |
| Replication | **not stated**; no n, no replicate scatter | — |

### 2.2 The three-phase decomposition and which phase the constants come from (Results, p. 84)

> "Fitting software indicates the decay is the **sum of three different exponentials**,
> indicating three separate phases. The **first phase, which represents 86 % of absorbance
> decrease**, is both NAC concentration and temperature dependent. The **second phase is NAC
> concentration independent and temperature dependent**. The **third phase is both temperature
> and NAC concentration independent, which suggests it is not part of the reaction process but
> an artifact of how components were mixed.** Reaction with GSH displayed the same three
> exponential decays."

> "**Analysis of the first phase** was carried out by plotting the observed rate constants as a
> function of NAC or GSH concentration generating liner [sic] plots, **the slope of which
> corresponds to the forward rate constants, k₁, of the bimolecular conjugate addition**
> (Figure 4B). **Nonzero intercepts indicates a reversible process, and the value of these
> nonzero intercept represents the reverse rate constant, k₋₁.**"

**★ This is the rate-law definition and it must be carried with every number below.**
The fitted model is
**k_obs = k₁·[RSH] + k₋₁**, applied to the **first phase only (86 % of the amplitude)**.
Therefore:
- **k₁ is second-order, units M⁻¹ s⁻¹**, and is defined **on total thiol [RSH], not on
  thiolate [RS⁻]** — the concentrations plotted in Fig. 4B are "[N-acetylcysteine] (M)", the
  analytical concentration. **Any pH transfer of k₁ must re-apply the thiolate fraction
  explicitly (§7.3).**
- **k₋₁ is first-order, units s⁻¹**, obtained as a **y-intercept extrapolated to zero thiol**,
  which is the least well-determined quantity in the paper (the authors put 7–10 % error on
  intercepts vs 4–6 % on slopes, §2.3).
- **The second phase — assigned to "the proton transfers involved in the tautomerization to the
  final catechol product" — is NOT included in k₁ or k₋₁.** So k₁/k₋₁ describes the reversible
  Michael addition to the *enolate/cyclohexadienone* intermediate, **not** formation of the
  final aromatic catechol adduct (Fig. 4D scheme: BPAQ + RS⁻ ⇌ enolate intermediate → *proton
  transfers* → catechol product). **The overall adduct formation is therefore effectively
  irreversible once tautomerization occurs; the reversibility measured here is on the first
  step only.** This is the single most important structural caveat for the model.

### 2.3 Error statement, verbatim (Results, p. 84)

> "**Error analysis for the line slopes was between 4 and 6 %, while the y-intercepts were
> between 7 and 10 %; the final ΔG‡ values reported are ca. ± 0.5 kcal/mol.**"

So: **k₁ carries ±4–6 %, k₋₁ carries ±7–10 %, and the (unusable) ΔG‡ carry ±0.5 kcal/mol.**
No error bars are printed on any individual cell of Figure 4C.

### 2.4 Preparative chemistry conditions (Experimental, pp. 82–83) — different from the kinetics

**Do not confuse these with the kinetic conditions.** The preparative runs are at room
temperature and much higher concentration.

| variable | NAC preparation | GSH preparation |
|---|---|---|
| BPA start | **0.250 g, 1.10 mmol** in **50 mL 3:2 CHCl₃:CH₃OH** | same, then BPAQ redissolved in **50 mL 50:50 methanol:water** |
| Oxidant | **IBX 0.337 g (1.20 mmol)**, **−15 °C, 12 h**; solution turns deep red after ~15 min | same |
| BPAQ workup | chloroform layer washed **6 × 50 mL with 0.1 M phosphate buffer pH 6.0**; HPLC showed **> 95 % BPAQ**, used without further purification | same |
| Thiol | **NAC 235 mg, 1.43 mmol** in **5 mL 0.1 M phosphate buffer pH 6.2** | **GSH 572 mg, 1.43 mmol** in **5 mL 0.1 M phosphate buffer pH 6.2** |
| Stoichiometry | **1.3 equiv thiol** | **1.3 equiv thiol** |
| Observation | > "The mixture **decolorized instantly** from bright red to a light yellow color" | > "stirred at rt for **15 min**" |
| Reaction time | stirred **10 min at rt** | **15 min at rt** |

⚠️ IBX safety note reproduced as printed: > "*Caution: IBX is explosive under impact if heated
to more than 200 °C.*"

---

## 3. FIGURE 4C — THE RATE CONSTANT TABLE (the quantitative core) **[M]/[F]**

**Anchor: Figure 4C, p. 85 (PDF page 5).** Caption as printed: *"(C) Forward and reverse rate
constants for the reaction BPAQ with NAC and GSH at various temperatures."*
Header block as printed: `pH =6.2` / `NAC` / `GSH`; column headers as printed:
`Temp (°C) | k₁ (M⁻¹s⁻¹) | k₋₁ (s⁻¹) | k₁ (M⁻¹s⁻¹) | k₋₁ (s⁻¹)`.
**Units quoted exactly as printed.** Transcribed from the 400 dpi raster; **all cells legible,
none unreadable.** The two GSH cells at 4.6 °C are printed as **em-dashes (—)**, i.e. not
measured (GSH insolubility), not unreadable.

| Temp (°C) | **NAC k₁ (M⁻¹s⁻¹)** | **NAC k₋₁ (s⁻¹)** | **GSH k₁ (M⁻¹s⁻¹)** | **GSH k₋₁ (s⁻¹)** |
|---|---|---|---|---|
| **19.4** | **496** | **88** | **1547** | **84** |
| **14.5** | **478** | **63** | **1403** | **57** |
| **9.6** | **446** | **48** | **1165** | **33** |
| **4.6** | **391** | **37** | **—** (not measured) | **—** (not measured) |

Provenance: **[F]** — every value is a regression parameter (slope or intercept) of the
Fig. 4B `k_obs` vs `[RSH]` lines, not a directly observed quantity.

### 3.1 Derived equilibrium constants **[Z]** (never printed by the authors)

K = k₁/k₋₁ for the **first (reversible Michael addition) step only**, units **M⁻¹**.

| Temp (°C) | **K_NAC (M⁻¹)** | **K_GSH (M⁻¹)** |
|---|---|---|
| 19.4 | **5.64** | **18.42** |
| 14.5 | **7.59** | **24.61** |
| 9.6 | **9.29** | **35.30** |
| 4.6 | **10.57** | — |

**K decreases with temperature for both thiols** ⇒ the addition is **exothermic**, and the
adduct is *less* favoured as temperature rises. **This is a directly food-relevant sign:
heating pushes the thiol–quinone Michael equilibrium back toward free thiol.**

### 3.2 van 't Hoff analysis of K **[Z]**

| quantity | **NAC** | **GSH** |
|---|---|---|
| ΔH°_rxn (first step) | **−6.81 kcal/mol** (−28.5 kJ/mol) | **−10.92 kcal/mol** (−45.7 kJ/mol) |
| ΔS°_rxn (first step) | **−19.8 cal mol⁻¹ K⁻¹** | **−31.6 cal mol⁻¹ K⁻¹** |
| R² of the van 't Hoff line | 0.9610 (4 pts) | 0.9972 (3 pts) |

**Internal consistency check [Z]:** ΔH°_rxn must equal Ea(forward) − Ea(reverse). From the §6.2
refits: NAC 2.58 − 9.39 = **−6.81** ✓ exact; GSH 4.76 − 15.69 = **−10.93** vs −10.92 ✓ exact.
The refit Arrhenius parameters and the van 't Hoff analysis are mutually consistent to the last
digit, which independently validates the Figure 4C transcription.

---

## 4. FIGURE 4D — THE PRINTED ACTIVATION PARAMETERS (transcribed, then refuted in §6)

**Anchor: Figure 4D, p. 85 (PDF page 5).** Caption as printed: *"(D) Free energies of
activation, **derived from Arrhenius plots**, for both forward and reverse conjugate
additions."* Panel D also carries the mechanism scheme
`BPAQ + RS⁻ ⇌(k₁/k₋₁) enolate intermediate —proton transfers→ catechol adduct`.
Transcribed verbatim from the 400 dpi raster; all four lines legible.

| system | line as printed |
|---|---|
| **NAC** | **ΔG‡_k1 = 9.2 kcal/mol (ΔH‡ = 5.9 kcal/mol, ΔS‡ = −12 eu)** |
| **NAC** | **ΔG‡_k−1 = 11.7 kcal/mol (ΔH‡ = 20.0 kcal/mol, ΔS‡ = +28 eu)** |
| **GSH** | **ΔG‡_k1 = 7.8 kcal/mol (ΔH‡ = 11.0 kcal/mol, ΔS‡ = +11 eu)** |
| **GSH** | **ΔG‡_k−1 = 11.2 kcal/mol (ΔH‡ = 35.6 kcal/mol, ΔS‡ = +82 eu)** |

Units as printed: **kcal/mol** for ΔG‡ and ΔH‡; **eu** (entropy units = cal mol⁻¹ K⁻¹) for ΔS‡.
Stated precision (p. 84): **ca. ± 0.5 kcal/mol** on the ΔG‡ values.

**Corresponding abstract text**, verbatim (p. 81):

> "Stopped-flow kinetic analysis reveals the 1,6-conjugate additions to be reversible with a
> **forward free energy of activation of 9.2 and 7.8 kcal/mol for the NAC and GSH reactions**,
> respectively. The bimolecular forward rate constant at 19.4 °C was approximately three time
> faster for GSH compared to NAC, **1547 vs 496 M⁻¹ s⁻¹**. The free energy of activation for the
> reverse reactions were similar, **11.7 and 11.2 kcal/mol** for NAC and GSH, respectively."

**Cited comparison [C]:** > "These ΔG‡ values are similar to computational calculations,
**10.0 kcal/mol**, that model the reaction of GSH with BPAQ **even though the modeled reaction
assumed a 1,4-conjugate addition**" — Kolšek, Sollner Dolenc & Mavri 2013, *Chem. Res. Toxicol.*
**26**, 106−111 (ref 24). **This 10.0 kcal/mol is a DFT/computational value from another paper,
not a measurement, and per the repo's no-DFT policy it should not be ingested as an anchor.**
Note also that the agreement the authors claim rests on the ΔG‡ values that §6 shows are
mis-derived — so the "similar to computational calculations" corroboration is spurious.

---

## 5. FIGURE 4A / 4B — the raw traces **[fig]**

**Anchor: Figure 4A/4B, p. 85.**
- **4A**: *"The decrease in absorption of a 125 µM solution of BPAQ in 50 % methanol when mixed
  with various, excess concentrations of NAC."* Ordinate **"Absorbance (400 nm)"**, ticks
  0 / 0.2 / 0.4 / 0.6; abscissa **"Time (s)", logarithmic**, ticks 0.001 / 0.01 / 0.1 / 1 / 10.
  Annotated **"increasing concentration of N-acetylcysteine"** with a downward arrow, and the
  three phases bracketed as **phase 1** (≈0.001–0.06 s), **phase 2** (≈0.1–1 s), **phase 3**
  (≈2–20 s). Starting absorbance ≈ **0.55**; plateau after phase 1 ≈ **0.09**; final ≈ **0.03**.
  **[fig]** — consistent with the stated 86 % of the absorbance decrease occurring in phase 1
  ([Z] check: (0.55 − 0.09)/(0.55 − 0.03) = **88 %** ✓).
- **4B**: *"The observed rate constants as a function of NAC concentration run at four different
  temperatures."* Ordinate **"k_obs (s⁻¹)"**, ticks 0 / 50 / 100 / 150; abscissa
  **"[N-acetylcysteine] (M)"**, ticks 0 / 0.05 / 0.1. Four labelled lines, top to bottom:
  **19.4 °C, 14.5 °C, 9.6 °C, 4.6 °C**. Three data points per line at
  **[NAC] = 0.050, 0.075, 0.10 M**.

⚠️ **DISCREPANCY [Z]:** the Experimental section states NAC concentrations of **"0.050, 0.10,
and 0.125 M"**, but Figure 4B's plotted points sit at approximately **0.050, 0.075 and 0.100 M**
and the x-axis stops at 0.1 — **no point is plotted at 0.125 M**. Either the middle
concentration was 0.075 M (not 0.10 M) and the highest was 0.10 M (not 0.125 M), or a point is
omitted from the figure. Not resolvable from the paper; **flagged, not fixed.**
**[fig]** digitised intercepts and slopes are not transcribed here because Figure 4C already
reports them as fitted values; the figure is recorded for provenance only.

The **GSH** equivalents of Fig. 4B, and **all Arrhenius plots for both forward and reverse
reactions (Figures S3 and S4)**, are in the Supporting Information, which is **not on disk**.
That is why §6's refits had to be done from the Figure 4C constants rather than checked against
the authors' own plots.

---

## 6. ★ RESOLUTION OF THE ΔG-VERSUS-k CONVENTION MISMATCH

The brief flagged a "ΔG-vs-k₂ convention mismatch". It is real, it is the paper's fault, and it
resolves cleanly. Three independent tests, all run against the paper's own Figure 4C numbers.

### 6.1 Test 1 — are the printed triples internally consistent? **YES.**

ΔG‡ = ΔH‡ − TΔS‡ at T = 292.55 K (19.4 °C):

| system | ΔH‡ − TΔS‡ **[Z]** | printed ΔG‡ | verdict |
|---|---|---|---|
| NAC k₁ | 5.9 − (292.55)(−0.012) = **9.41** | 9.2 | consistent (within the stated ±0.5) |
| NAC k₋₁ | 20.0 − (292.55)(+0.028) = **11.81** | 11.7 | consistent |
| GSH k₁ | 11.0 − (292.55)(+0.011) = **7.78** | 7.8 | consistent |
| GSH k₋₁ | 35.6 − (292.55)(+0.082) = **11.61** | 11.2 | consistent |

So the four triples are a coherent set *on their own terms*. The problem is not arithmetic
within Figure 4D.

### 6.2 Test 2 — do the printed ΔH‡ match an Arrhenius refit of Figure 4C? **NO — every one is 2.303× too large.**

Ordinary least squares of ln k on 1/T, using the Figure 4C constants and nothing else:

| series | **refit Ea, kcal/mol [Z]** | refit A **[Z]** | R² | **printed ΔH‡** | **ratio printed/refit** | **2.303 × refit** | dev. |
|---|---|---|---|---|---|---|---|
| NAC k₁ (4 T) | **2.58** | 4.267 × 10⁴ M⁻¹s⁻¹ | 0.9405 | 5.9 | **2.287** | **5.94** | **+0.7 %** |
| NAC k₋₁ (4 T) | **9.39** | 8.862 × 10⁸ s⁻¹ | 0.9938 | 20.0 | **2.130** | 21.63 | +8.1 % |
| GSH k₁ (3 T) | **4.76** | 5.690 × 10⁶ M⁻¹s⁻¹ | 0.9721 | 11.0 | **2.311** | **10.96** | **−0.3 %** |
| GSH k₋₁ (3 T) | **15.69** | 4.489 × 10¹³ s⁻¹ | 0.9923 | 35.6 | **2.269** | **36.13** | **+1.5 %** |

**Four independent series, four ratios clustered at 2.13–2.31, mean 2.25.** The natural
constant in that band is **ln(10) = 2.303**, and multiplying the refit by exactly 2.303
reproduces **three of the four printed values to within 1.5 %**.

**★ DIAGNOSIS: the authors applied a spurious factor of ln(10) to their Arrhenius slopes — a
natural-log plot (ln k vs 1/T) processed with the base-10 formula Ea = −2.303·R·slope.** Every
printed ΔH‡ in Figure 4D is therefore **2.303× too large**. The ΔS‡ values were then evidently
chosen to bring ΔG‡ = ΔH‡ − TΔS‡ back to a plausible magnitude, which is why Test 1 passes and
Test 3 fails.

*(The NAC k₋₁ series is the one 8 % outlier. Refitting it on only the upper three temperatures,
to match the GSH protocol, moves it further away (Ea 10.16, ×2.303 = 23.40, +17 %), so the
four-point fit reported above is the best available and the outlier is not an artefact of point
selection. Given the stated 7–10 % error on intercepts, an 8 % residual on the reverse series is
within the paper's own error budget.)*

### 6.3 Test 3 — do the printed ΔG‡ regenerate the measured k under Eyring? **NO — by 2 to 4 orders of magnitude.**

Eyring, k = (k_B T/h)·exp(−ΔG‡/RT), at T = 292.55 K where k_B T/h = **6.096 × 10¹² s⁻¹**:

| series | measured k (19.4 °C) | **TRUE Eyring ΔG‡ from that k, kcal/mol [Z]** | printed ΔG‡ | **k regenerated from the printed ΔG‡** | error |
|---|---|---|---|---|---|
| NAC k₁ | 496 M⁻¹s⁻¹ | **13.51** | 9.2 | 8.17 × 10⁵ | **1647× too fast** |
| NAC k₋₁ | 88 s⁻¹ | **14.51** | 11.7 | 1.11 × 10⁴ | **126× too fast** |
| GSH k₁ | 1547 M⁻¹s⁻¹ | **12.84** | 7.8 | 9.08 × 10⁶ | **5870× too fast** |
| GSH k₋₁ | 84 s⁻¹ | **14.54** | 11.2 | 2.62 × 10⁴ | **312× too fast** |

**The printed ΔG‡ values are not Gibbs free energies of activation for these reactions.** The
true Eyring barriers all cluster at **12.8–14.5 kcal/mol** — a narrow, chemically sensible band
for a fast bimolecular Michael addition — against printed values spread over 7.8–11.7.

### 6.4 The honest replacement set **[Z]**

If any activation parameter is wanted from this paper, use these, derived from Figure 4C alone
and labelled as a repo-side derivation, **never** the Figure 4D values:

| quantity | **NAC forward** | **NAC reverse** | **GSH forward** | **GSH reverse** |
|---|---|---|---|---|
| **Ea (Arrhenius)** | **2.58 kcal/mol** = **10.8 kJ/mol** | **9.39 kcal/mol** = **39.3 kJ/mol** | **4.76 kcal/mol** = **19.9 kJ/mol** | **15.69 kcal/mol** = **65.6 kJ/mol** |
| **A (pre-exponential)** | 4.267 × 10⁴ M⁻¹s⁻¹ | 8.862 × 10⁸ s⁻¹ | 5.690 × 10⁶ M⁻¹s⁻¹ | 4.489 × 10¹³ s⁻¹ |
| **R² of the refit** | 0.9405 | 0.9938 | 0.9721 | 0.9923 |
| **ΔG‡ (Eyring, 19.4 °C)** | **13.51 kcal/mol** | **14.51 kcal/mol** | **12.84 kcal/mol** | **14.54 kcal/mol** |
| **ΔH‡ = Ea − RT (19.4 °C)** | 2.00 kcal/mol | 8.81 kcal/mol | 4.18 kcal/mol | 15.10 kcal/mol |
| **ΔS‡ = (ΔH‡ − ΔG‡)/T** | **−39.3 eu** | **−19.5 eu** | **−29.6 eu** | **+1.9 eu** |

⚠️ **Honesty note on the refit itself.** The forward series are fit over a **14.8 K span
(NAC) or 9.8 K span (GSH)** with **four or three points** and R² as low as **0.94**. The NAC
forward k changes by only **1.27×** across the whole range. **These Ea values are weakly
determined and their error bars are large — plausibly ±30–50 % on the forward channels.** They
are better than the published ones, not good in absolute terms.

### 6.5 Note on the "Stack vs Liu 2023" cross-check

The brief said to use Stack's own numbers where they conflict with a Liu 2023. **The only
"Liu 2023" in the repo is a matrix-path hold-out benchmark**
(`docs/reference/FIT_HOLDOUT_DECLARATION.md:117`, `VALIDATION_CONTRACT.md:327` — a
hexanal/matrix dataset), **not a thiol–quinone kinetics paper, and no thiol–quinone Liu 2023
PDF exists in `data/articles/`.** The instruction is therefore vacuous as stated: **there is no
conflict on disk to resolve, and Stack's Figure 4C values stand unopposed.** Flagging this so
the orchestrator can correct the reference if a different Liu 2023 was meant.

---

## 7. PRODUCT DISTRIBUTION, REGIOCHEMISTRY AND MECHANISM

### 7.1 Product yields **[M]**

**Anchor: Results, pp. 83–84, and Experimental, pp. 82–83.** Two different yield bases are
reported and must not be mixed.

| system | product | **HPLC relative % (crude, 1.3 equiv thiol)** | **isolated %, large scale** | mass isolated |
|---|---|---|---|---|
| **BPAQ + NAC** | **5-NAC** (3-hydroxy-5-NAC-BPA) | **89 %** (major) | **81 %** | **319 mg** |
| **BPAQ + NAC** | **2-NAC** (3-hydroxy-2-NAC-BPA) | **11 %** (minor) | **7 %** | **27 mg** |
| **BPAQ + NAC** | 3-OHBPA (reduction byproduct) | "trace amounts" | not isolated | — |
| **BPAQ + NAC** | any bis-adduct | **none observed** at 1.3, 2.0 or 5.0 equiv NAC | — | — |
| **BPAQ + GSH** | **2,5-diGSH** (bis-adduct) | present at 1.3 equiv | **17 %** | **151 mg** |
| **BPAQ + GSH** | **5-GSH** | major mono-adduct | **36 %** | **212 mg** |
| **BPAQ + GSH** | **2-GSH** | minor mono-adduct | **7 %** | **43 mg** |
| **BPAQ + GSH** | **3-OHBPA** (reduction byproduct) | > "observed in **much larger quantities** when compared to the NAC reaction" | **32 %** | **81 mg** |

⚠️ The paper prints "43 mg of 3-hydoxy-2-GSH (7%)" — the compound name is missing "-BPA" and
"hydroxy" is misspelled; reproduced as printed.

### 7.2 Regiochemistry — a clean structural result **[M]**

| finding | evidence | anchor |
|---|---|---|
| **Both NAC and GSH add exclusively 1,6**; the major (5-SR) and minor (2-SR) isomers are *both* 1,6-addition products | HMBC three-bond correlations: **5-NAC shows three cross peaks** between the quaternary sp³ carbon and arene protons H2, H6, H2′, whereas a 6-NAC isomer would show only two | Fig. 2A, p. 84; Scheme 2, p. 82 |
| **No 1,4-addition product (6-SR) was detected**, mono- or di-adducted | > "No products, mono- or diadducted, were the result of 1,4-addition conjugate addition." | p. 84 |
| The **2,5-diGSH** bis-adduct is **the result of two separate 1,6-conjugate additions** | two resolved α-CH₂ sets (one downfield, one upfield of GSH); NOE shows only the downfield methylene correlates to the lone remaining catechol arene proton | p. 84 |
| **GSH excess suppresses the bis-adduct**: > "At 5-fold excess of GSH, **little of the di-2,5-GSH adduct or the reduced 3-OHBPA catechol** were observed" | Fig. 3 (0.5 mM BPAQ scale, 1.3 / 2.0 / 5.0 equiv GSH) | p. 84 |
| **Proposed bis-adduct mechanism**: oxidation of the newly formed catechol mono-adduct **by unreacted BPAQ**, which also generates the 3-OHBPA reduction product — explaining why 3-OHBPA and 2,5-diGSH rise together | Scheme 3, p. 85 | p. 84 |
| GSH adducts have **lower oxidation potentials than NAC adducts (30–65 mV)** ⇒ GSH is the stronger electron donor ⇒ larger k₂ in Scheme 3 ⇒ bis-adduct forms with GSH but not NAC | **[C]** Macedo et al. 2007, *J. Health Sci.* **53**, 31−42 (ref 28), measured on α-methyldopamine and N-methyl-α-dopamine o-quinones, **not on BPAQ** | p. 85 |

⚠️ **Naming collision the orchestrator must not trip over:** Scheme 3 uses **k₁ and k₂** for the
*bis-adduct* pathway (mono-adduct formation vs oxidation of the mono-adduct), while Figure 4C/4D
use **k₁ and k₋₁** for the *stopped-flow* forward/reverse Michael addition. **The k₁ of Scheme 3
and the k₁ of Figure 4C are not the same constant.** The Scheme 3 k₂ is **never measured** —
the only statement about it is the qualitative > "If k₁ is greater than k₂, increasing GSH
concentration relative to BPAQ would decrease the rate of 5-GSH-BPAQ formation and subsequent
formation of 2,5-diGSH."

### 7.3 ★ The thiolate explanation — the paper's claim, and where it over-reaches **[Z]**

The authors' explanation for GSH being ~3× faster than NAC, verbatim (p. 86):

> "**GSH, pKa 8.8, is more acidic than NAC, pKa 9.5.**²⁹ Thus, at a pH of 6.2, the effective
> concentration of thiolate anion would be much larger with GSH solutions compared to NAC."

pKa values are **[C]**, cited to Trujillo & Radi 2002, *Arch. Biochem. Biophys.* **397**, 91−98
(ref 29). Testing the claim quantitatively against the paper's own k₁:

| quantity **[Z]** | NAC | GSH |
|---|---|---|
| pKa **[C]** | 9.5 | 8.8 |
| thiolate fraction at pH 6.2 | **5.01 × 10⁻⁴** | **2.51 × 10⁻³** |
| observed k₁ (19.4 °C) on **total thiol** | 496 M⁻¹s⁻¹ | 1547 M⁻¹s⁻¹ |
| **k₁ per thiolate anion** | **9.90 × 10⁵ M⁻¹s⁻¹** | **6.17 × 10⁵ M⁻¹s⁻¹** |

- **Thiolate availability alone predicts a 5.00× GSH advantage. The observed advantage is
  3.12×.** The authors' explanation therefore **over-predicts by 1.6×**.
- **Per thiolate anion, NAC is the FASTER nucleophile (1.60× faster than GSH).** The GSH
  advantage is entirely a speciation effect and is in fact *partly offset* by lower intrinsic
  anion reactivity — plausibly steric, given GSH's tripeptide bulk.
- **Modelling consequence:** a repo term that transfers this k₁ to another pH **must** apply the
  thiolate fraction explicitly, using **k₁_intrinsic ≈ 6–10 × 10⁵ M⁻¹ s⁻¹ on [RS⁻]**, not the
  pH-6.2 apparent constant. Applying 496 or 1547 M⁻¹s⁻¹ at a different pH would be wrong by the
  ratio of thiolate fractions — e.g. **~16× at pH 7.4**.
- ⚠️ The pKa values are aqueous literature values applied to a **50 % methanol** medium, where
  both the thiol pKa and the apparent pH shift substantially. **The absolute thiolate fractions
  above are indicative only.**

**Cited supporting precedent [C]:** Freeman et al. (Baker et al. 2007, *J. Biol. Chem.* **282**,
31085−31093, ref 21) observed "much faster conjugate addition of GSH relative to NAC with
α,β-unsaturated nitro compounds derived from fatty acids and postulated this was a function of
thiol acidity."

### 7.4 Reversibility precedents cited by the authors **[C]** (no numbers attached)

| system | observation | ref |
|---|---|---|
| 4-hydroxyequilenin o-quinone + GSH | adducts **liberated 4-hydroxyequilenin and GSH when exposed to NADPH reduction** | 18 (Peng et al. 2010) |
| quercetin o-quinone + GSH | isomeric adducts **interconvert when allowed to equilibrate in solution** | 30 (Boersma et al. 2000) |
| allyl / benzyl isothiocyanate + GSH | conjugates **release the isothiocyanate electrophile in vitro** | 31 (Bruggeman et al. 1986) |

**These are the qualitative literature basis for treating thiol–electrophile capture as
reversible.** Useful as structural support; no rate constants.

---

## 8. WHAT IS NOT IN THIS PAPER (declared gaps)

| missing | consequence |
|---|---|
| **Supporting Information** (Figures S3/S4, the actual Arrhenius plots; GSH concentration-dependence plots; all NMR) | **not on disk.** §6's diagnosis had to be reconstructed from Figure 4C rather than checked against the authors' own plots. **Fetching the SI would let the 2.303 diagnosis be confirmed directly from their plot axes** |
| Any temperature above **19.4 °C** | the entire dataset sits at 4.6–19.4 °C. **No thermal-process temperature is covered.** Extrapolating k₁ to 100 °C+ spans ~80 K on an Ea determined over ~15 K with R² 0.94 |
| Any pH other than **6.2** | no pH axis; the thiolate correction in §7.3 is a repo-side inference, not a measurement |
| Any aqueous-only measurement | **everything is in 50 % methanol.** No water-only control, and no comment on the methanol effect on either rate or pKa |
| Ionic strength control / activity corrections | not reported |
| **n, replicates, or per-cell error bars** | only the global "4–6 % on slopes, 7–10 % on intercepts" statement |
| Rate constants for the **second phase** (tautomerization) | explicitly identified as temperature-dependent and concentration-independent, but **never quantified**. So the paper measures the reversible step and *not* the step that commits the adduct |
| **k₂ of Scheme 3** (oxidation of mono-adduct by BPAQ) | never measured; only a qualitative inequality |
| Any food-relevant thiol (FFT, MFT, methanethiol, cysteine) | **the nucleophiles are NAC and GSH only** |
| Any food-relevant quinone (caffeoylquinone, 4-ethylcatechol quinone, o-benzoquinone) | **the electrophile is BPAQ only** |

---

## NEW-PARAMETER TABLE (consolidated)

All conditions: **0.050 M phosphate buffer pH 6.2 in 50 % methanol**, stopped-flow at 400 nm,
pseudo-first-order in excess thiol, first-phase (86 % amplitude) analysis,
**k_obs = k₁[RSH] + k₋₁**. k₁ is defined on **total thiol**, not thiolate.

| parameter | value | units (as printed) | conditions | anchor | provenance |
|---|---|---|---|---|---|
| k₁, BPAQ + NAC | **496** | M⁻¹ s⁻¹ | 19.4 °C, pH 6.2, 50 % MeOH | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₁, BPAQ + NAC | **478** | M⁻¹ s⁻¹ | 14.5 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₁, BPAQ + NAC | **446** | M⁻¹ s⁻¹ | 9.6 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₁, BPAQ + NAC | **391** | M⁻¹ s⁻¹ | 4.6 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₋₁, BPAQ–NAC adduct | **88** | s⁻¹ | 19.4 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₋₁, BPAQ–NAC adduct | **63** | s⁻¹ | 14.5 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₋₁, BPAQ–NAC adduct | **48** | s⁻¹ | 9.6 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₋₁, BPAQ–NAC adduct | **37** | s⁻¹ | 4.6 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₁, BPAQ + GSH | **1547** | M⁻¹ s⁻¹ | 19.4 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₁, BPAQ + GSH | **1403** | M⁻¹ s⁻¹ | 14.5 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₁, BPAQ + GSH | **1165** | M⁻¹ s⁻¹ | 9.6 °C | Fig. 4C, p. 85 | [F] ±4–6 % |
| k₋₁, BPAQ–GSH adduct | **84** | s⁻¹ | 19.4 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₋₁, BPAQ–GSH adduct | **57** | s⁻¹ | 14.5 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| k₋₁, BPAQ–GSH adduct | **33** | s⁻¹ | 9.6 °C | Fig. 4C, p. 85 | [F] ±7–10 % |
| K = k₁/k₋₁, NAC | **5.64 / 7.59 / 9.29 / 10.57** | M⁻¹ | 19.4 / 14.5 / 9.6 / 4.6 °C | from Fig. 4C | **[Z]** |
| K = k₁/k₋₁, GSH | **18.42 / 24.61 / 35.30** | M⁻¹ | 19.4 / 14.5 / 9.6 °C | from Fig. 4C | **[Z]** |
| ΔH°_rxn (first step), NAC | **−6.81** | kcal/mol (−28.5 kJ/mol) | van 't Hoff, R² 0.9610 | §3.2 | **[Z]** |
| ΔH°_rxn (first step), GSH | **−10.92** | kcal/mol (−45.7 kJ/mol) | van 't Hoff, R² 0.9972 | §3.2 | **[Z]** |
| ΔS°_rxn (first step), NAC | **−19.8** | cal mol⁻¹ K⁻¹ | van 't Hoff | §3.2 | **[Z]** |
| ΔS°_rxn (first step), GSH | **−31.6** | cal mol⁻¹ K⁻¹ | van 't Hoff | §3.2 | **[Z]** |
| **Ea (Arrhenius), NAC forward** | **2.58** | kcal/mol (10.8 kJ/mol), A = 4.267 × 10⁴ M⁻¹s⁻¹ | 4.6–19.4 °C, R² 0.9405 | §6.2/§6.4 | **[Z]** |
| **Ea (Arrhenius), NAC reverse** | **9.39** | kcal/mol (39.3 kJ/mol), A = 8.862 × 10⁸ s⁻¹ | 4.6–19.4 °C, R² 0.9938 | §6.2/§6.4 | **[Z]** |
| **Ea (Arrhenius), GSH forward** | **4.76** | kcal/mol (19.9 kJ/mol), A = 5.690 × 10⁶ M⁻¹s⁻¹ | 9.6–19.4 °C, R² 0.9721 | §6.2/§6.4 | **[Z]** |
| **Ea (Arrhenius), GSH reverse** | **15.69** | kcal/mol (65.6 kJ/mol), A = 4.489 × 10¹³ s⁻¹ | 9.6–19.4 °C, R² 0.9923 | §6.2/§6.4 | **[Z]** |
| ΔG‡ (Eyring, true), NAC fwd / rev | **13.51 / 14.51** | kcal/mol | 19.4 °C, from measured k | §6.3 | **[Z]** |
| ΔG‡ (Eyring, true), GSH fwd / rev | **12.84 / 14.54** | kcal/mol | 19.4 °C, from measured k | §6.3 | **[Z]** |
| ΔG‡ forward, NAC — **AS PRINTED** | **9.2** | kcal/mol | ±0.5 | Fig. 4D, p. 85 | [F] — ⚠️ **CORRUPT, see §6** |
| ΔG‡ reverse, NAC — **AS PRINTED** | **11.7** | kcal/mol | ±0.5 | Fig. 4D, p. 85 | [F] — ⚠️ **CORRUPT** |
| ΔG‡ forward, GSH — **AS PRINTED** | **7.8** | kcal/mol | ±0.5 | Fig. 4D, p. 85 | [F] — ⚠️ **CORRUPT** |
| ΔG‡ reverse, GSH — **AS PRINTED** | **11.2** | kcal/mol | ±0.5 | Fig. 4D, p. 85 | [F] — ⚠️ **CORRUPT** |
| ΔH‡ printed, NAC fwd / rev | **5.9 / 20.0** | kcal/mol | — | Fig. 4D, p. 85 | [F] — ⚠️ **2.303× too large** |
| ΔH‡ printed, GSH fwd / rev | **11.0 / 35.6** | kcal/mol | — | Fig. 4D, p. 85 | [F] — ⚠️ **2.303× too large** |
| ΔS‡ printed, NAC fwd / rev | **−12 / +28** | eu | — | Fig. 4D, p. 85 | [F] — ⚠️ back-computed from the corrupt ΔH‡ |
| ΔS‡ printed, GSH fwd / rev | **+11 / +82** | eu | — | Fig. 4D, p. 85 | [F] — ⚠️ back-computed from the corrupt ΔH‡ |
| **k₁ per thiolate anion, NAC** | **9.90 × 10⁵** | M⁻¹ s⁻¹ on [RS⁻] | 19.4 °C, pKa 9.5 | §7.3 | **[Z]** (pKa is [C]) |
| **k₁ per thiolate anion, GSH** | **6.17 × 10⁵** | M⁻¹ s⁻¹ on [RS⁻] | 19.4 °C, pKa 8.8 | §7.3 | **[Z]** (pKa is [C]) |
| thiol pKa, NAC / GSH | **9.5 / 8.8** | — | aqueous | p. 86, ref 29 | **[C]** Trujillo & Radi 2002 |
| GSH/NAC k₁ ratio | **3.12×** | — | 19.4 °C, pH 6.2 | Fig. 4C | **[Z]** ("approximately three time faster", abstract, [M]) |
| product split, BPAQ+NAC | **89 % 5-NAC / 11 % 2-NAC** | HPLC % | 1.3 equiv NAC, rt, 10 min | p. 84 | [M] |
| isolated yields, BPAQ+NAC | **81 % 5-NAC (319 mg) / 7 % 2-NAC (27 mg)** | % | large scale, 1.10 mmol | p. 83 | [M] |
| isolated yields, BPAQ+GSH | **17 % 2,5-diGSH (151 mg) / 36 % 5-GSH (212 mg) / 7 % 2-GSH (43 mg) / 32 % 3-OHBPA (81 mg)** | % | 1.3 equiv GSH | p. 84 | [M] |
| regiochemistry | **exclusively 1,6-addition; no 1,4-adduct detected** | — | both thiols, all stoichiometries | Fig. 2A, p. 84 | [M] |
| phase-1 amplitude | **86 %** of the absorbance decrease | % | both thiols | p. 84 | [M] |
| adduct oxidation potential, GSH vs NAC | **30–65 mV lower** for GSH adducts | mV | α-methyldopamine / N-methyl-α-dopamine o-quinones — **not BPAQ** | p. 85, ref 28 | **[C]** Macedo 2007 |
| computational ΔG‡, GSH + BPAQ | **10.0** | kcal/mol | assumed **1,4**-addition | p. 84, ref 24 | **[C]** Kolšek 2013 — ⚠️ **DFT; excluded by the repo's no-DFT policy** |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> Stack 2018 is **not** in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. **A declaration
> amendment is required before any wave may fit it.** This section is a proposal only; the
> declaration was not edited.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Fig. 4C NAC k₁ at 4.6, 9.6, 14.5 °C** (391, 446, 478 M⁻¹s⁻¹) | **FIT** | temperature | Pins the absolute magnitude of thiol + o-quinone Michael addition. The corpus currently has **no measured second-order thiol–quinone rate constant at all**; this is the anchor for the whole quinone-sink lane |
| **Fig. 4C NAC k₁ at 19.4 °C** (496 M⁻¹s⁻¹) | **★ HOLD-OUT** | temperature (top of range) | The only genuine T-extrapolation available within the NAC arm. Weak (a 1.27× total span) but honest, and free |
| **Fig. 4C NAC k₋₁, all four temperatures** (88, 63, 48, 37 s⁻¹) | **FIT** | — | **The corpus's only reverse rate constant for a thiol–electrophile adduct.** Without it the repo can only model thiol capture as irreversible, which §3.1 shows is wrong at the first step. Fit all four: the reverse Ea (9.39 kcal/mol) is the better-determined of the two channels (R² 0.9938) |
| **Fig. 4C GSH k₁ and k₋₁, all rows** | **★ HOLD-OUT (whole arm)** | **nucleophile identity** | This is the strongest cut available. Fit the NAC arm, predict the GSH arm. It tests exactly the thing the model needs to get right — **whether a thiol's reactivity can be predicted from its pKa and sterics** — and §7.3 shows the naive pKa answer is wrong by 1.6×, so it is a genuinely discriminating test rather than a formality |
| **Derived K = k₁/k₋₁ and the van 't Hoff ΔH°_rxn (§3.1–3.2)** | **neither — derived, would double-count** | — | K is algebraically determined by the k values already assigned above. Ingesting it as an independent constraint double-counts. **Use it as a diagnostic**: if the model reproduces both k's it must reproduce K, and the negative ΔH° (heating releases thiol) is a qualitative behaviour worth asserting in a regression test |
| **Fig. 4D ΔG‡, ΔH‡, ΔS‡ — all twelve printed values** | **★ neither — DO NOT INGEST** | — | **ΔH‡ is 2.303× too large across all four series (§6.2); ΔG‡ regenerates the measured k 126–5870× too fast (§6.3).** If activation parameters are wanted, use the **[Z] refit in §6.4** and label it a repo-side derivation from Stack Fig. 4C, never as a Stack result |
| **The §6.4 refit Ea values** | **FIT only as loose priors, never as scored targets** | — | Determined over 10–15 K with R² down to 0.94; plausibly ±30–50 %. Adequate to say "the forward barrier is small and the reverse barrier is 3–4× larger"; **not adequate to extrapolate to process temperatures** |
| **The §7.3 per-thiolate intrinsic k₁** (9.90 × 10⁵ and 6.17 × 10⁵ M⁻¹s⁻¹) | **FIT (as the pH-transferable form)** | — | If the repo needs a thiol–quinone rate at any pH other than 6.2, **this is the form to use**, with the 50 %-methanol pKa caveat attached. Ingesting the pH-6.2 apparent constant instead would be wrong by ~16× at pH 7.4 |
| **Product regiochemistry (100 % 1,6-addition), yields, phase-1 amplitude** | **qualitative ingest** | — | Structural facts, not kinetics. The 1,6-exclusivity is a useful constraint on any adduct-structure enumeration |
| **Kolšek 2013 computational ΔG‡ = 10.0 kcal/mol** | **★ REJECT** | — | DFT, from a third paper, on an assumed-wrong (1,4) regiochemistry, and the "agreement" the authors claim rests on the corrupt ΔG‡ values. **Excluded by the repo's no-DFT policy on all three counts** |

### ⚠️ Transferability warnings the orchestrator must carry forward

1. **Solvent.** Every number is in **50 % methanol**. Both the thiol pKa and the quinone's
   electrophilicity shift in mixed solvent, and the paper applies aqueous pKa values without
   comment. **Any use in an aqueous food matrix carries an unquantified solvent correction.**
2. **Temperature range.** 4.6–19.4 °C. **Nothing here constrains behaviour at brewing (92 °C),
   pasteurisation (121 °C) or roasting (180 °C+) temperatures.** The forward k changes by only
   1.27× across the whole measured range, so the Ea is barely identified.
3. **Chemical analogy is the real assumption.** BPAQ is a *para*-substituted bisphenol
   o-quinone; the food-relevant electrophiles (chlorogenic-acid quinones, 4-ethylcatechol
   quinone, o-benzoquinone) differ in ring substitution and redox potential. **Treat k₁ ≈ 5 ×
   10² M⁻¹ s⁻¹ (total thiol, pH 6.2) as an order-of-magnitude anchor for the class, not as a
   transferable constant for any specific food quinone.**
4. **The reversibility is on the first step only.** k₁/k₋₁ describes addition to the enolate
   intermediate; the subsequent proton-transfer/tautomerization step (phase 2) is **not
   quantified** and commits the adduct. **A model that treats the measured K as the equilibrium
   constant for net adduct formation will over-estimate thiol recovery.** This is the single
   most likely way to misuse this paper.
5. **Both thiols here are far better nucleophiles than 2-furfurylthiol is likely to be**, and
   neither is present in a food matrix at the concentrations used (10–125 mM). The rate
   constants are transferable; the *conditions* are not.
