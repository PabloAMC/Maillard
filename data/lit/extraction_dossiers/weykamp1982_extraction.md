# Weykamp & Penders 1982 (Clin. Chim. Acta 125, 341–350; PII 0009-8981(82)90265-0) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/weykamp1982.pdf` (10 pp.). Read method: both — OCR text layer
(Acrobat 3.0 Capture plug-in; heavily corrupted, subscripts lost, "10⁻³" rendered as "10-j", "10e4",
"lo-'") plus 200-dpi **page rasters of PDF pp. 1, 3, 5, 6, 7**, which are the pages carrying every
number in this dossier. **Every value below was read off the raster, not the OCR layer.** The OCR
layer is unsafe for this file and must not be used by any downstream wave.

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY

| field | value |
|---|---|
| authors | Cas W. Weykamp* and Theo J. Penders |
| title | "Mechanism and speed of reactions between haemoglobin and glucose; consequences for the measurement of glycosylated haemoglobins in patient material" |
| venue | *Clinica Chimica Acta*, Elsevier Biomedical Press; article code CCA 2312 |
| volume/pages/year | **125 (1982) 341–350** |
| DOI | none printed; PII 0009-8981(82)90265-0. Article code `0009-8981/82/0000-0000/$02.75` |
| received/revision | "(Received April 20th; revision July 12th, 1982)" |
| affiliation | Department of Clinical Chemistry, Winterswijkse Ziekenhuizen, Postbus 9005, 7100 GG Winterswijk (The Netherlands) |
| PDF character | **OCR scan** (Acrobat 3.0 Capture), no page-count metadata, image-backed with a poor text layer |

Expected identity was "Schiff base formation kinetics at 37 °C, forward k1, reverse k−1, equilibrium
K; haemoglobin-glycation / glucose-amine". **Confirmed on all counts.**

## 1. ONE-PARAGRAPH VERDICT — AND THE UNIT ANSWER UP FRONT

This paper gives exactly four kinetic constants, all at a single temperature (37 °C), for the
glucose + haemoglobin system in intact human erythrocytes: the Schiff-base (aldimine) **formation**
constant k₁, the Schiff-base **degradation** constant k₂ (which is the "k₋₁" of Maillard notation),
the Amadori/HbA₁c **formation** constant k₃, and the Schiff-base **equilibrium** constant K₁. **There
is no temperature series and therefore no Ea** — the paper's only temperature remark is qualitative:
"In an experiment, which we have not described further here, the reaction took 2½ times longer at
room temperature". There is also **no stated pH anywhere in the paper** (see §2 — this is the single
most important caveat). The forward constant k₁ is *derived*, not measured: k₂, k₃ and K₁ were
measured and k₁ = K₁·k₂. The paper additionally gives a six-species algebraic model (Table I) and one
validation table (Table II). Three internal-consistency defects are documented in §6.

### THE UNIT QUESTION — RESOLVED. k₁ is in **(mmol/l)⁻¹ h⁻¹**, i.e. **l mmol⁻¹ h⁻¹**.

**It is NOT M⁻¹ s⁻¹. It is NOT M⁻¹ h⁻¹. It is NOT l mol⁻¹ min⁻¹.**

The units are established by one sentence, in the Summary, p. 341 (PDF page 1). **Verbatim,
raster-read:**
> "Reaction speed constants were found for the formation of HbA₁c (4.2 · 10⁻³) and the formation and
> degradation of the Schiff Base (9.22 · 10⁻⁴ and 0.435, respectively) **for concentrations in
> mmol/l and time in hours** at 37°C. An equilibrium constant of 2.12 · 10⁻³ was found for the
> reversible formation of the Schiff Base."

That clause — "for concentrations in mmol/l and time in hours" — is the **only** place in the paper
where units are assigned to any of the four constants. **No units appear next to any constant in the
body text, and there is no table of constants.** The assignment is therefore:

| constant | numeric value | units implied by the Summary clause | dimensional reasoning |
|---|---|---|---|
| k₁ (Schiff formation, 2nd order) | 9.22 · 10⁻⁴ | **l mmol⁻¹ h⁻¹** = (mmol/l)⁻¹ h⁻¹ | Reaction 1 is Hb + G → HbG, bimolecular; concentrations in mmol/l, time in h |
| k₂ (Schiff degradation, 1st order) | 0.435 | **h⁻¹** | Reaction 2 is HbG → Hb + G, unimolecular; "t, time in hours" is stated again under the k₂ derivation |
| k₃ (HbA₁c formation, 1st order) | 4.20 · 10⁻³ | **h⁻¹** | Reaction 3 is HbG → HbA₁c, unimolecular |
| K₁ (Schiff equilibrium) | 2.12 · 10⁻³ | **l mmol⁻¹** = (mmol/l)⁻¹ | K₁ = HbG/(Hb·G) with all three in mmol/l |

**Is the assignment ambiguous?** No — but it is *implicit*, and one specific ambiguity must be
flagged: the Summary clause is written once and applies to all four numbers collectively. The
individual reaction orders (which determine whether a given constant carries the concentration unit
at all) come from the mechanism in Fig. 1 and from the first-order integrated law printed on p. 344,
not from the units clause. **The internal 4-digit consistency check confirms the reading:**
K₁ · k₂ = 2.12 × 10⁻³ × 0.435 = **9.222 × 10⁻⁴** = the printed k₁ [Z]. This closes only if K₁ is
(mmol/l)⁻¹ and k₂ is h⁻¹, giving k₁ in (mmol/l)⁻¹ h⁻¹. **Resolution is unambiguous.**

Conversions, **all [Z], never printed by the authors:**

| constant | as printed | × 10³ (mmol → mol) | in M⁻¹ h⁻¹ / h⁻¹ | in M⁻¹ s⁻¹ / s⁻¹ |
|---|---|---|---|---|
| k₁ | 9.22 × 10⁻⁴ l mmol⁻¹ h⁻¹ | 0.922 l mol⁻¹ h⁻¹ | **0.922 M⁻¹ h⁻¹** | **2.56 × 10⁻⁴ M⁻¹ s⁻¹** |
| k₂ | 0.435 h⁻¹ | — | **0.435 h⁻¹** | **1.208 × 10⁻⁴ s⁻¹** |
| k₃ | 4.20 × 10⁻³ h⁻¹ | — | **4.20 × 10⁻³ h⁻¹** | **1.167 × 10⁻⁶ s⁻¹** |
| K₁ | 2.12 × 10⁻³ l mmol⁻¹ | 2.12 l mol⁻¹ | **2.12 M⁻¹** | **2.12 M⁻¹** |

Also [Z]: Schiff-base half-life at 37 °C = ln2 / 0.435 = **1.594 h**; the paper's own equivalent
statement (Discussion, p. 348) is "after 5, 10 and 15 h, 10⁻¹, 10⁻² and 10⁻³ parts, respectively, of
the original HbG quantity remain", which corresponds to a decadic decay time of ~5 h, i.e.
k₂ = ln10/5 = 0.461 h⁻¹ — consistent with 0.435 to ~6%.

## 2. SYSTEM DEFINITION — verbatim

**Reaction scheme (Fig. 1, p. 342 / PDF page 2).** Caption verbatim (OCR-corrupted; the substantive
definitions are legible): "Reactions of haemoglobin and glucose with indication of rate of reaction,
specific rate constants and equilibrium constants. Hb, free haemoglobin; G, free glucose; HbG,
haemoglobin with one β-chain occupied by glucose in aldimine form; HbA₁c, haemoglobin with one
β-chain occupied by glucose in ketoamine form; HbGG, haemoglobin with both β-chains occupied by
glucose in aldimine form; HbA₁c₁c, haemoglobin with both β-chains occupied by glucose in ketoamine
form; HbA₁cG, haemoglobin with both β-chains occupied, one by glucose in aldimine form and one by
glucose in ketoamine form; Kₙ, equilibrium constant for equilibrium reaction n; Sₙ, rate of reaction
for reaction n; kₙ, specific rate constant for reaction n."

Numbering, therefore: **Reaction 1 = Hb + G → HbG (Schiff formation, k₁); Reaction 2 = HbG → Hb + G
(Schiff reversion, k₂); Reaction 3 = HbG → HbA₁c (Amadori, k₃).** In Maillard notation the paper's
**k₂ is k₋₁**.

**Symmetry relations, p. 343 (PDF page 3), verbatim, raster-read:**
> "sites, HbA₁c only one, so that Reaction 4 can proceed half as fast
> k₄ = ½k₁.
> In this way the following may also be calculated: k₅ = k₂, k₆ = k₃, k₇ = ½k₁, k₈ = 2k₂, k₉ = 2k₃,
> K₂ = ½K₁, K₃ = ¼K₁. The number of variables is hereby reduced to four (k₁, k₂, k₃ and K₁). Because
> K₁ = k₁/k₂ in equilibrium, it follows than whenever three of the four values can be measured all
> other constants can be calculated. **In the investigation k₂, k₃ and K₁ were measured
> experimentally and k₁ was derived from K₁ and k₂.**"

That last sentence fixes the provenance: **k₂, k₃, K₁ are [M]; k₁ is [Z]/[F] — derived, not
measured.**

**Methods, "Measurement of reaction constants", p. 343 (PDF page 3), verbatim, raster-read:**
> "Fresh EDTA blood from a non-diabetic laboratory worker was washed 3 times in 10 times the volume
> of physiological saline. 4 ml sterile physiological saline containing 100 mmol/l glucose and 72
> mmol/l NaHCO₃ was mixed with 4 ml of the washed erythrocytes. The glucose contained ¹⁴C-labelled
> glucose (D-[2-¹⁴C]-glucose, New England Nuclear, 549 Albany Street, Boston, MA, USA). The specific
> activity was 0.0855 nCi/nmol glucose. The mixture was incubated at 37°C on a rotor. After 8 h four
> washes were carried out with 20 times the cell volume of physiological saline at 37°C. The washed
> erythrocytes were then incubated at 37°C with one cell volume of physiological saline. After 3 h a
> single washing with 10 times the cell volume of physiological saline was carried out after which
> the erythrocytes were incubated again with 1 cell volume of physiological saline at 37°C.
> 500-µl samples were taken at 0, ½, 1, 2, 4, 8 and 12 h. The sample-cells were washed once in 20
> times the cell volume of physiological saline. 200 µl erythrocytes were haemolysed with 250 µl
> water. The haemoglobin concentration was determined from this haemolysate. 400 µl haemolysate were
> mixed with 400 µl oxalic acid, 0.3 mol/l, and warmed to 100°C to break down the
> glucose-haemoglobin band [7]. 1 ml 20% trichloroacetic acid was added to this reaction mixture to
> sediment the protein. After centrifugation the upper fluid layer with free ¹⁴C-glucose was
> pipetted into a test tube. The sediment was resuspended again 3 times with 1 ml water and
> centrifuged. The pooled supernatants were evaporated down to 300 µl. The concentrated fluid was
> put into counting vials, mixed with 10 ml counting fluid, and measured in a scintillation counter
> (Delta 300, Searle Analytic, Searle Nederland BV, Wiegerbruinlaan 75, Uithoorn, The Netherlands).
> Following this 29.5 nCi ¹⁴C was added as an internal standard and counted."

**Model-checking experiment, p. 343–344 (PDF pages 3–4), verbatim:**
> "100 µl erythrocytes, which had been washed 3 times in 10 cell volumes physiological saline, were
> haemolysed with 100 µl water and mixed with 1 ml glucose solution with a concentration varying
> from 0 to 3720 mmol/l. The haemolysates were incubated at 37°C for 10 h. The components were
> separated by isoelectric focusing [5] (Power Supply 2003, Multiphor 2117, Focusing medium:
> Ampholine-PAG 1804-131, LKB, LKB-Produkter AB, S-161, 26 Bromma, Sweden). Scanning was carried out
> with a densitometer after fixation. (Beckman R112 with R-115 integration-unit, Beckman Instruments
> Incorporated, 200 South Kraemer Boulevard, Brea, CA 92621, USA)."

**Condition table**

| condition variable | value as printed |
|---|---|
| temperature | **37 °C** everywhere (incubation, washes, model-check). One qualitative room-temperature remark, no number for T |
| **pH** | **NOT STATED ANYWHERE IN THE PAPER.** The word "pH" does not appear. The medium is "physiological saline" (composition not given beyond the name); the *loading* step alone contains 72 mmol/l NaHCO₃. The reaction compartment for k₂, K₁ and k₃ is the **intact erythrocyte cytosol**, whose pH the authors never measure, state, or assume |
| buffer identity/molarity | physiological saline (unbuffered, composition unstated) for all kinetic incubations; 72 mmol/l NaHCO₃ present only in the 8 h ¹⁴C-loading mixture |
| hot-pH correction | not applicable; no pH at all |
| glucose | 100 mmol/l in the loading mixture; **47.4 mmol/l at the equilibrium point used for K₁**; 0–3720 mmol/l in the model check; 1250 mmol/l for Table II; physiological reference value 5.4 mmol/l (cited [10]) |
| haemoglobin | total 8.70 mmol/l at the equilibrium point; free Hb 7.42 mmol/l; physiological reference 9.3 mmol/l (cited [10]) |
| tracer | D-[2-¹⁴C]-glucose, specific activity 0.0855 nCi/nmol glucose |
| internal standard | **29.5 nCi ¹⁴C added after the first count** |
| vessel / volume | 4 ml saline + 4 ml packed erythrocytes; 500 µl samples; 200 µl erythrocytes + 250 µl water for haemolysis |
| headspace / atmosphere | not stated |
| agitation | "incubated at 37°C **on a rotor**" |
| quench | 400 µl haemolysate + 400 µl oxalic acid 0.3 mol/l, **warmed to 100 °C** to break the glucose–haemoglobin (Schiff) bond [7], then 1 ml 20% TCA to precipitate protein |
| sampling times | 0, ½, 1, 2, 4, 8 and 12 h (7 points) |
| analytical method | liquid scintillation counting of released ¹⁴C-glucose (for HbG); isoelectric focusing + densitometry (for Hb / HbG / HbA₁c / HbGG separation) |
| replication (n) | **n = 1 donor.** "Fresh EDTA blood from **a** non-diabetic laboratory worker". No replicate count anywhere |
| error bars | **NONE.** No ±, no SD, no SE anywhere in the paper. The only statistics are one correlation coefficient (0.9986) and one p-value (p < 0.001) |

**HbG quantification formula, p. 344 (PDF page 4), as printed:**
`P = (x·u) / (y·v·w·(y − x)) × 100%` — the OCR renders this as
"P = X.U / ,_.U.W(J-x) x 100%" and the raster is legible but the variable subscripting is degraded.
Definitions as printed: "P, percentage HbG; u, nCi ¹⁴C in the internal standard; v, µl haemolysate
used in the counting; w, haemoglobin concentration in the haemolysate in mmol/l; x, dpm in the
sample; y, dpm in the sample plus internal standard; z, specific activity." **The exact algebraic
arrangement of this formula is partially unreadable** (the denominator's grouping cannot be resolved
with confidence from the scan, and the definitions list a variable `z` that does not visibly appear
in the printed formula). Recorded as **unreadable** rather than reconstructed.

## 3. DETERMINATION OF k₂ (= the Maillard k₋₁) — the headline number

**Anchor: Results, "Determination of k₂", pp. 344–345 (PDF pages 4–5).**

Rate law as printed, p. 344:
`ln (HbG₀ / HbGᵢ) = k₂ · t`
"HbG₀, HbG concentration at time 0; HbGᵢ, HbG concentration at time i; k₂, reaction speed constant
for Reaction 2; **t, time in hours**."

Justification as printed, p. 344: "In a glucose-free milieu only Reaction 2 will run (Fig. 1) as
Reaction 1 will not occur and S₃ ≪ S₂ (this hypothesis will be proven later). For Reaction 2 a first
order reaction is most probable".

**Result, p. 345 (PDF page 5), verbatim, raster-read:**
> "We found a correlation coefficient of 0.9986, which demonstrates the linear connection, and thus
> the order of the reaction (p < 0.001).
> The reaction speed constant k₂, in agreement with the slope of the line, is **0.435**."

| parameter | value | units | conditions | provenance |
|---|---|---|---|---|
| k₂ (Schiff-base degradation = Maillard k₋₁) | **0.435** | **h⁻¹** (from the Summary units clause) | 37 °C, intact erythrocytes in glucose-free physiological saline, pH not stated, n = 1 donor, 7 time points 0–12 h | **[M]** |
| correlation coefficient of the ln-linear fit | 0.9986 | — | as above | [F] |
| significance | p < 0.001 | — | as above | [F] |

Supporting figure: **Figure 2, p. 344 (PDF page 4)**, caption verbatim: "Decrease of HbG as a
function of time." Axes: percentage HbG vs "Time in hours", x-axis ticks 2, 4, 6, 8, 10. I have not
digitised it; the fitted k₂ is printed in the text so digitisation is unnecessary.

## 4. DETERMINATION OF K₁, AND THE DERIVATION OF k₁

**Anchor: Results, "Determination of K₁" and "Calculation of k₁", p. 345 (PDF page 5).** Verbatim,
raster-read:
> "Before sampling at t = 0, the erythrocytes were washed in physiological saline. During this time
> (0.717 h) conversion of HbG into Hb and G has taken place. The HbG concentration at time −0.717,
> the equilibrium concentration of HbG, can be calculated from the following equation [9]:
> log HbGᵢ = (−k₂ · t)/2.303 + log HbG₀.
> For t = −0.717 a HbG of 8.56% was found. The total haemoglobin concentration at equilibrium was
> 8.70 mmol/l. Of this 6.2% was HbA₁c (determined in the original blood) and 8.56% was HbG. The free
> haemoglobin concentration was then:
> ((100 − 6.2 − 8.56)/100) × 8.7 = 7.42 mmol/l.
> The glucose level at equilibrium was 47.4 mmol/l.
> At equilibrium the following equation applies:
> K₁ = HbG / (Hb · G).
> Substitution of the values determined results in **K₁ = 2.12 · 10⁻³**."

> "**Calculation of k₁.** The reaction speed constant k₁ is derived from
> k₁ = K₁ · k₂ = **9.22 · 10⁻⁴**."

| parameter | value | units | conditions | provenance |
|---|---|---|---|---|
| back-extrapolation time | 0.717 | h | wash duration before t = 0 | [M] |
| HbG at equilibrium | 8.56 | % of total Hb | 37 °C | [M] (back-extrapolated with k₂ → strictly [Z] from k₂) |
| total haemoglobin at equilibrium | 8.70 | mmol/l | 37 °C | [M] |
| HbA₁c in the original blood | 6.2 | % | donor baseline | [M] |
| free Hb at equilibrium | 7.42 | mmol/l | 37 °C | [F] (computed, arithmetic printed) |
| glucose at equilibrium | 47.4 | mmol/l | 37 °C | [M] |
| **K₁** | **2.12 × 10⁻³** | **l mmol⁻¹** (= 2.12 M⁻¹ [Z]) | 37 °C, pH not stated | **[M]** (authors call it measured; strictly it is computed from five measured quantities, arithmetic printed) |
| **k₁** | **9.22 × 10⁻⁴** | **l mmol⁻¹ h⁻¹** (= 0.922 M⁻¹ h⁻¹, = 2.56 × 10⁻⁴ M⁻¹ s⁻¹ [Z]) | 37 °C, pH not stated | **[F]/[Z] — derived as K₁·k₂, NOT measured** |

**My reproduction [Z], all steps:**
- Free Hb: (100 − 6.2 − 8.56)/100 × 8.70 = 0.8524 × 8.70 = **7.416 mmol/l** → printed 7.42 ✓
- HbG: 0.0856 × 8.70 = **0.7447 mmol/l** (this intermediate is not printed)
- K₁ = 0.7447 / (7.42 × 47.4) = 0.7447 / 351.71 = **2.118 × 10⁻³** → printed 2.12 × 10⁻³ ✓
- k₁ = 2.12 × 10⁻³ × 0.435 = **9.222 × 10⁻⁴** → printed 9.22 × 10⁻⁴ ✓ (4-digit agreement)

**The K₁ / k₁ / k₂ triad is internally consistent to 4 significant figures.**

## 5. DETERMINATION OF k₃ (Amadori / HbA₁c formation)

**Anchor: Results, "Determination of k₃", pp. 345–346 (PDF pages 5–6).** Verbatim, raster-read:
> "K₃ could be determined by measuring the increase in HbA₁c with time. The problem with this is that
> t cannot be made large enough to give a reasonably measurable rise in HbA₁c. Erythrocytes convert
> glucose to lactate, resulting in total cell haemolysis within 2 days. The results of measurement
> within 48 h proved too inaccurate, and we therefore chose another approach.
> S₃ = ΔHbA₁c/Δt = k₃ · HbG.
> HbG = K₁ · Hb · G,
> from which follows: k₃ = ΔHbA₁c / (Δt · K₁ · Hb · G)"

(**Printing error, reported not fixed:** the section heading and first sentence say "K₃", capital K,
where the quantity being determined is the rate constant k₃. The equations that follow all use lower
case k₃.)

> "In this equation we have filled in the values for non-diabetics: (Hb = 9.3 mmol/l [10]), G = 5.4
> mmol/l [10], Δt = 1380 h, average age of erythrocytes [10], K₁ = 2.12 · 10⁻³ and ΔHbA₁c = 0.577
> (median value determined according to [5].
> A k₃ value of **4.20 · 10⁻³** was found."

| parameter | value | units as printed | conditions | provenance |
|---|---|---|---|---|
| Hb (non-diabetic reference) | 9.3 | mmol/l | — | **[C]** — Geigy, *Wissenschaftliche Tabellen*, Documenta Geigy 1960 (ref. 10) |
| G (non-diabetic reference) | 5.4 | mmol/l | — | **[C]** — same source |
| Δt (mean erythrocyte age) | 1380 | h | — | **[C]** — same source |
| ΔHbA₁c | 0.577 | (units not printed; context implies mmol/l) | median, "determined according to [5]" (Mortensen 1980, IEF method) | [M] / [C] |
| K₁ | 2.12 × 10⁻³ | l mmol⁻¹ | 37 °C | [M] (§4) |
| **k₃** | **4.20 × 10⁻³** | **h⁻¹** | 37 °C, pH not stated | **[F]** — computed, and **it does not reproduce (see §6.2)** |

## 6. INTERNAL-CONSISTENCY CHECKS AND DEFECTS — three findings

### 6.1 The k₁ / k₂ / K₁ triad closes exactly. ✓
K₁ · k₂ = 9.222 × 10⁻⁴ vs printed k₁ = 9.22 × 10⁻⁴ [Z]. **4-digit agreement.** This is the check that
also fixes the units (§1).

### 6.2 The printed k₃ does NOT reproduce from the printed inputs — 7% high.
My arithmetic [Z], using the authors' own five inputs:
```
k3 = ΔHbA1c / (Δt · K1 · Hb · G)
   = 0.577 / (1380 × 2.12e-3 × 9.3 × 5.4)
   = 0.577 / (1380 × 2.12e-3 × 50.22)
   = 0.577 / (1380 × 0.106466)
   = 0.577 / 146.92
   = 3.928e-3
```
**Printed: 4.20 × 10⁻³. Recomputed: 3.93 × 10⁻³. The printed value is 6.9% higher than the printed
inputs support.** Reported as a defect, not corrected. Note that Table I (§7) *is* built on
4.20 × 10⁻³ (see §7 check), so the model and the abstract are self-consistent with each other and
only the k₃ derivation line is off.

### 6.3 The printed ratio k₃/k₂ = 1.03 · 10⁻² is wrong by a factor of ~1.07; the reciprocal is what was computed.
**Anchor: "Testing hypothesis", p. 346 (PDF page 6).** As printed:
`S₃/S₂ = k₃/k₂ = 1.03 · 10⁻². From this follows that S₃ ≪ S₂.`

My arithmetic [Z]: k₃/k₂ = 4.20 × 10⁻³ / 0.435 = **9.66 × 10⁻³**, not 1.03 × 10⁻².
Conversely k₂/k₃ = 0.435 / 4.20 × 10⁻³ = **103.6**.
The printed "1.03 · 10⁻²" is the mantissa of **103.6** with a 10⁻² exponent attached — i.e. the
authors computed the ratio the other way round (103.6) and then wrote it as if it were the
reciprocal. The Discussion confirms the intended meaning, verbatim, p. 348: "we found formation of
HbA₁c to be **100 times slower** than that of HbG", which matches 103.6.
**The conclusion (S₃ ≪ S₂) is unaffected; only the printed number is wrong.** Reported, not fixed.

## 7. TABLE I — the six-species algebraic model

**Anchor: TABLE I, p. 346 (PDF page 6).** Title as printed:
"MATHEMATICAL MODEL FOR CALCULATION OF HAEMOGLOBIN AND GLYCOSYLATED HAEMOGLOBINS"

Raster-read from PDF p. 6 (the OCR layer mangles every exponent here and must not be used).

**Input, as printed**

| symbol | definition as printed |
|---|---|
| w | Total haemoglobin in mmol/l = w |
| v | Glucose in mmol/l = v |
| u | Time in h = u |

**Calculation, as printed**

```
              w
Hb = ──────────────────────────────────────────────────────────────────────────────────
                       8.9 × 10⁻⁶ vu
     1 + 2.12 × 10⁻³ v + ───────────── + 3.96 × 10⁻¹¹ v²u² + 9.4 × 10⁻⁹ v²u + 1.12 × 10⁻⁶ v²
                       1 + 4.4 × 10⁻⁶ vu

Hb        = x
HbG       = 2.12 × 10⁻³ xv
                8.9 × 10⁻⁶ xvu
HbA1c     = ──────────────────
             1 + 4.4 × 10⁻⁶ vu
HbA1c1c   = 3.96 × 10⁻¹¹ xv²u²
HbA1cG    = 9.4 × 10⁻⁹ xv²u
HbGG      = 1.12 × 10⁻⁶ xv²
```

**Output, as printed:** "Concentrations of free Hb, HbG, HbA₁c, HbA₁c₁c, HbA₁cG and HbGG in mmol/l"

| model coefficient | value as printed | provenance | my [Z] reconstruction from the four constants |
|---|---|---|---|
| HbG prefactor | 2.12 × 10⁻³ | [F] | = K₁ ✓ exact |
| HbA₁c numerator | 8.9 × 10⁻⁶ | [F] | = k₃ · K₁ = 4.20 × 10⁻³ × 2.12 × 10⁻³ = 8.904 × 10⁻⁶ ✓ **confirms k₃ = 4.20 × 10⁻³ is the value actually used** |
| HbA₁c denominator | 4.4 × 10⁻⁶ | [F] | not reconstructible from the four constants alone; **unreadable as to origin** |
| HbA₁c₁c prefactor | 3.96 × 10⁻¹¹ | [F] | origin not stated; **unreadable as to origin** |
| HbA₁cG prefactor | 9.4 × 10⁻⁹ | [F] | origin not stated; **unreadable as to origin** |
| HbGG prefactor | 1.12 × 10⁻⁶ | [F] | = K₁ · (¼K₁) = 2.12 × 10⁻³ × 5.30 × 10⁻⁴ = 1.1236 × 10⁻⁶ ✓ **consistent with the stated K₃ = ¼K₁ relation (p. 343)** |

## 8. TABLE II — model validation

**Anchor: TABLE II, p. 347 (PDF page 7).** Title as printed:
"OBSERVED AND CALCULATED PERCENTAGES OF FREE HAEMOGLOBIN AND GLYCOSYLATED HAEMOGLOBINS IN AN
HAEMOLYSATE WITH A GLUCOSE CONCENTRATION OF 1250 mmol/l AFTER 10 h INCUBATION AT 37°C"

Column headers as printed: `Component` | `Observed` | `Calculated`. **No units column; the title
says percentages. No error bars.**

| Component | Observed | Calculated | provenance |
|---|---|---|---|
| Hb | 19 | 17 | Observed [M], Calculated [F] |
| HbA₁c | 12 | 8 | [M] / [F] |
| HbG | 46 | 43 | [M] / [F] |
| HbGG | 23 | 28 | [M] / [F] |
| HbA₁c₁c | – | *(blank)* | — |
| HbA₁cG | – | 2 | [F] |

**Transcription note:** the HbA₁c₁c "Calculated" cell is **blank in the print** (raster-verified on
PDF p. 7); the HbA₁c₁c "Observed" cell holds an em-dash. The HbA₁cG "Observed" cell holds an em-dash
and the Calculated cell holds 2. Recorded exactly as printed rather than filled in.

[Z] check: observed column sums to 19 + 12 + 46 + 23 = **100**; calculated column sums to
17 + 8 + 43 + 28 + 2 = **98**. The authors' verdict, verbatim (p. 348): "the good agreement between
the values confirms the correctness of the model." The largest single discrepancy is HbA₁c
(12 observed vs 8 calculated, a 50% relative error), which the authors do not comment on.

## 9. OTHER PRINTED QUANTITIES

**Anchor: Discussion, p. 348 (PDF page 8), verbatim:**
> "The assumed [1] slow formation of HbA₁c and quick formation and breakdown of HbG is confirmed by
> our investigation: we found formation of HbA₁c to be **100 times slower** than that of HbG.
> Furthermore, the demonstration of HbGG (the 'double Schiff Base') is a strong argument in favour of
> the supposition that glucose only binds to the β-chains.
> From the reaction-speed comparison for the dissociation of HbG in a glucose-free milieu it follows
> that after **5, 10 and 15 h, 10⁻¹, 10⁻² and 10⁻³ parts**, respectively, of the original HbG
> quantity remain. We use a **10-h incubation** time to remove HbG from patient samples.
> In an experiment, which we have not described further here, **the reaction took 2½ times longer at
> room temperature**, so that **25 h incubation** would be necessary to remove 99% of the original
> amount of HbG."

**This last sentence is the paper's ONLY temperature information beyond 37 °C, and it is
unquantified** ("room temperature", no number). [Z]: if "2½ times longer" means k₂ is 2.5× smaller,
then k₂(RT) ≈ 0.174 h⁻¹. With "room temperature" taken as 20 °C (an assumption the paper does not
make), a two-point Arrhenius estimate would give Ea ≈ R·ln(2.5)/(1/293.15 − 1/310.15) = 8.314 ×
0.9163 / 1.869 × 10⁻⁴ = **4.1 × 10⁴ J/mol ≈ 41 kJ/mol**. **This number is [Z] built on an
un-stated temperature and must NOT be treated as an Ea from this paper.** It is recorded only so
that no downstream wave invents it independently and then treats it as measured.

**Clinical simulation outputs, Figs. 4–6, pp. 348–349 (PDF pages 8–9)**, all [fig], none digitised:
- Fig. 4 caption: "HbG, HbA₁c and HbG/HbA₁c + HbG (the fraction HbG in fast haemoglobins) as a
  function of a constant glucose concentration." Text result [M]: "the HbG forms a constant **1/7**
  part of the fast haemoglobins, independent of the glucose concentration."
- Fig. 5 caption: "Measured fast haemoglobins (HbG + HbA₁c) as a function of the glucose
  concentration at the moment of blood sampling." Text result [M]: "for glucose concentrations of
  **2 and 36 mmol/l**, glycosylated haemoglobin values of **12.5 and 16.5**, respectively, are found.
  In both cases the HbA₁c level is **12%**."
- Fig. 6 caption: "Concentration of HbG, HbA₁c, and fast haemoglobins during a 4-day period as the
  result of changing glucose concentration."
- Fig. 3 caption, p. 347 (PDF page 7): "Isoelectric focusing patterns of reaction products of
  haemoglobin and glucose after 10 h incubation at 37 °C with a concentration of glucose varying
  from 0 to **3100 mmol/l** (one pattern shown in detail). At a glucose concentration of 0 only Hb
  (**92%**) and HbA₁c (**8%**) are present. With increasing glucose concentration HbG appears and at
  the highest glucose concentrations even the double Schiff Base HbGG. HbA₁c remains constant and Hb
  decreases to **5%** at the highest glucose level."
  (**Inconsistency, reported not fixed:** the Methods say the model-check glucose range was
  "0 to 3720 mmol/l"; the Fig. 3 caption says "0 to 3100 mmol/l".)

## 10. DEFECTS AND CAVEATS TO CARRY FORWARD

1. **NO pH ANYWHERE.** The word does not occur in the paper. Any use of these constants must carry an
   explicit "pH unstated; intraerythrocytic, physiological saline" flag. Do not silently assume 7.4.
2. **NO Ea; single temperature (37 °C).** The room-temperature remark is unquantified.
3. **k₁ is derived, not measured** (k₁ = K₁·k₂). Only k₂, k₃ and K₁ are measurements.
4. **n = 1 donor. No error bars anywhere.**
5. **k₃ does not reproduce from the printed inputs** (3.93 × 10⁻³ recomputed vs 4.20 × 10⁻³ printed;
   §6.2). Table I nevertheless uses 4.20 × 10⁻³.
6. **The printed k₃/k₂ = 1.03 × 10⁻² is wrong**; the correct ratio from the printed constants is
   9.66 × 10⁻³ (§6.3).
7. **Table I coefficients 4.4 × 10⁻⁶, 3.96 × 10⁻¹¹ and 9.4 × 10⁻⁹ have no stated derivation** and
   cannot be reconstructed from k₁/k₂/k₃/K₁ alone.
8. **Table II's HbA₁c₁c "Calculated" cell is blank in the print**; observed column sums to 100, the
   calculated column to 98.
9. **Fig. 3 caption (3100 mmol/l) contradicts the Methods (3720 mmol/l).**
10. **The HbG quantification formula on p. 344 is partially unreadable** in the scan (§2).
11. **The OCR text layer for this file is unusable for numbers.** Every exponent is corrupted. Any
    re-extraction must go through page rasters.
12. **This is a biomedical, intracellular system**, not a food matrix: intact erythrocytes,
    haemoglobin β-chain N-terminal valine, physiological saline, 37 °C. Transfer to food kinetics is
    an extrapolation of ~60–140 °C in temperature and of the entire matrix.

## NEW-PARAMETER TABLE (consolidated)

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| k₂ — Schiff-base (aldimine) degradation, **= Maillard k₋₁** | 0.435 | h⁻¹ (units from the Summary clause "concentrations in mmol/l and time in hours") | 37 °C, intact human erythrocytes, glucose-free physiological saline, **pH not stated**, n = 1 donor, 7 points 0–12 h | Results "Determination of k₂", p. 345 / PDF p. 5; value also in Summary, p. 341 / PDF p. 1 | **[M]** |
| k₂ converted | 1.208 × 10⁻⁴ | s⁻¹ | as above | — | **[Z]** |
| Schiff-base half-life | 1.594 | h | as above | — | **[Z]** |
| K₁ — Schiff-base equilibrium constant | 2.12 × 10⁻³ | l mmol⁻¹ (= (mmol/l)⁻¹) | 37 °C, Hb_total 8.70 mmol/l, free Hb 7.42 mmol/l, G 47.4 mmol/l, **pH not stated** | Results "Determination of K₁", p. 345 / PDF p. 5; Summary p. 341 / PDF p. 1 | **[M]** |
| K₁ converted | 2.12 | M⁻¹ (l mol⁻¹) | as above | — | **[Z]** |
| k₁ — Schiff-base formation, 2nd order | 9.22 × 10⁻⁴ | **l mmol⁻¹ h⁻¹** (= (mmol/l)⁻¹ h⁻¹) | 37 °C, **pH not stated** | Results "Calculation of k₁", p. 345 / PDF p. 5; Summary p. 341 / PDF p. 1 | **[F]/[Z] — derived as K₁·k₂, NOT measured** |
| k₁ converted | 0.922 | M⁻¹ h⁻¹ | as above | — | **[Z]** |
| k₁ converted | 2.56 × 10⁻⁴ | M⁻¹ s⁻¹ | as above | — | **[Z]** |
| k₃ — HbA₁c (Amadori) formation, 1st order | 4.20 × 10⁻³ | h⁻¹ | 37 °C, computed from Hb 9.3 mmol/l, G 5.4 mmol/l, Δt 1380 h, ΔHbA₁c 0.577, K₁ 2.12 × 10⁻³ | Results "Determination of k₃", pp. 345–346 / PDF pp. 5–6; Summary p. 341 / PDF p. 1 | **[F]** — and does not reproduce (§6.2) |
| k₃ converted | 1.167 × 10⁻⁶ | s⁻¹ | as above | — | **[Z]** |
| k₃ recomputed from the paper's own inputs | 3.93 × 10⁻³ | h⁻¹ | as above | §6.2 | **[Z]** — flags a 6.9% discrepancy |
| ln-linear fit correlation coefficient for k₂ | 0.9986 | — | 37 °C, 7 points | p. 345 / PDF p. 5 | [F] |
| significance of the k₂ fit | p < 0.001 | — | as above | p. 345 / PDF p. 5 | [F] |
| k₄ = ½k₁; k₅ = k₂; k₆ = k₃; k₇ = ½k₁; k₈ = 2k₂; k₉ = 2k₃; K₂ = ½K₁; K₃ = ¼K₁ | (statistical-factor relations) | — | two equivalent β-chains | p. 343 / PDF p. 3 | [F] (assumed, not measured) |
| Table I: HbG prefactor | 2.12 × 10⁻³ | (mmol/l)⁻¹ | 37 °C | Table I, p. 346 / PDF p. 6 | [F] (= K₁) |
| Table I: HbA₁c numerator | 8.9 × 10⁻⁶ | (mmol/l)⁻¹ h⁻¹ | 37 °C | Table I, p. 346 / PDF p. 6 | [F] (= k₃·K₁ ✓ [Z]) |
| Table I: HbA₁c denominator coefficient | 4.4 × 10⁻⁶ | (mmol/l · h)⁻¹ | 37 °C | Table I, p. 346 / PDF p. 6 | [F] — **origin not stated** |
| Table I: HbA₁c₁c prefactor | 3.96 × 10⁻¹¹ | (mmol/l)⁻² h⁻² | 37 °C | Table I, p. 346 / PDF p. 6 | [F] — **origin not stated** |
| Table I: HbA₁cG prefactor | 9.4 × 10⁻⁹ | (mmol/l)⁻² h⁻¹ | 37 °C | Table I, p. 346 / PDF p. 6 | [F] — **origin not stated** |
| Table I: HbGG prefactor | 1.12 × 10⁻⁶ | (mmol/l)⁻² | 37 °C | Table I, p. 346 / PDF p. 6 | [F] (= K₁·¼K₁ ✓ [Z]) |
| Table II observed / calculated, glucose 1250 mmol/l, 10 h, 37 °C | Hb 19/17; HbA₁c 12/8; HbG 46/43; HbGG 23/28; HbA₁cG –/2; HbA₁c₁c –/*(blank)* | % | 37 °C, haemolysate | Table II, p. 347 / PDF p. 7 | [M] / [F] |
| HbG fraction of "fast haemoglobins" | 1/7 | (dimensionless) | model output, glucose-independent | Discussion p. 348 / PDF p. 8 (Fig. 4) | [F] |
| HbG removal by incubation | 10⁻¹, 10⁻², 10⁻³ remaining after 5, 10, 15 h | (fraction) | 37 °C, glucose-free saline | Discussion p. 348 / PDF p. 8 | [F] |
| room-temperature slowdown | "2½ times longer"; 25 h to remove 99% | (dimensionless) | **"room temperature" — no number given** | Discussion p. 348 / PDF p. 8 | [M], qualitative |
| implied Ea if RT = 20 °C is assumed | ~41 | kJ mol⁻¹ | **BUILT ON AN ASSUMPTION THE PAPER DOES NOT MAKE** | §9 | **[Z] — DO NOT USE AS AN Ea** |
| **pH** | **NOT STATED ANYWHERE IN THE PAPER** | — | — | — | — |
| **Activation energy Ea** | **NOT PRESENT IN THIS PAPER** | — | — | — | — |

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Standing status note.** `docs/reference/FIT_HOLDOUT_DECLARATION.md` Amendment 1 (2026-08-28)
already contains the row: *"Weykamp & Penders 1982 Schiff k−1 = 0.435 /h at 37 C (biomedical) —
**FIT (as a prior, rate_transfer: not_licensed to food T without the Ea)** — only measured Schiff
reverse rate in any literature; internal k1/k−1/K consistency verified to 4 digits."* This extraction
**confirms** that row: the 4-digit consistency is verified here in §6.1, and the absence of an Ea is
confirmed in §9. The proposals below are additive to that standing row and do not move it.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **k₂ = 0.435 h⁻¹** (Schiff reversion at 37 °C) | **FIT — as a prior only**, already declared; keep the existing `rate_transfer: not_licensed` guard | temperature (this is the low-T anchor; food data are the high-T arm) | It is the only directly measured Schiff-base reverse rate in the literature at any temperature, and it is the one number in this paper that is (a) measured, (b) from a 7-point ln-linear fit with r = 0.9986, and (c) unit-resolved. Its role as a *prior* rather than a scored target is right, because it is a single-donor biomedical measurement with no error bar. |
| **K₁ = 2.12 × 10⁻³ l mmol⁻¹ (= 2.12 M⁻¹)** | **FIT — as a prior, same guard** | temperature | Independent of k₂ in its derivation (it comes from an equilibrium concentration balance, not from a rate fit), so it is a genuine second constraint on the Schiff branch, not a restatement of k₂. |
| **k₁ = 9.22 × 10⁻⁴ l mmol⁻¹ h⁻¹** | **DO NOT FIT SEPARATELY — DERIVED** | — | k₁ = K₁·k₂ exactly. Fitting k₁ *and* k₂ *and* K₁ would triple-count one pair of measurements. Use k₂ + K₁ and let k₁ fall out, or use k₁ alone; never all three. **This is the single most important circularity in this paper.** |
| **k₃ = 4.20 × 10⁻³ h⁻¹** (Amadori formation) | **HOLD-OUT, with a defect flag** | reaction step (Schiff arm FIT → Amadori arm held out) | Once k₁/k₂/K₁ set the Schiff branch, predicting the Amadori rate at 37 °C is a genuine orthogonal test of the next mechanistic step. **But it must be scored against a band, not a point:** the printed 4.20 × 10⁻³ does not reproduce from the paper's own inputs (§6.2 gives 3.93 × 10⁻³), so the honest target is 3.9–4.2 × 10⁻³ h⁻¹. Also note k₃ is itself computed *using* K₁, which weakens its independence — score it as a weak hold-out. |
| **k₂/k₃ ≈ 100** (Amadori 100× slower than Schiff reversion) | **HOLD-OUT (ordinal)** | — | The robust, defect-free version of the k₃ datum: the *ratio* survives both arithmetic defects. Cheap and discriminating. |
| **Table II observed column** (Hb 19, HbA₁c 12, HbG 46, HbGG 23 at 1250 mmol/l glucose, 10 h, 37 °C) | **HOLD-OUT** | glucose concentration (1250 mmol/l is 230× physiological) | This is the only *speciation* measurement in the paper and it was made by an orthogonal method (isoelectric focusing + densitometry) rather than by ¹⁴C counting. Predicting the four-way split at extreme glucose is a real test of whether the model's double-adduct (HbGG) branch is right. **Do not fit to the Calculated column** — that is the authors' own model output, not data. |
| **Table I coefficients (2.12 × 10⁻³, 8.9 × 10⁻⁶, 1.12 × 10⁻⁶)** | **DO NOT USE — REDUNDANT** | — | Reconstructed here as K₁, k₃·K₁ and ¼K₁² respectively. They are the four constants re-expressed, not new information. |
| **Table I coefficients (4.4 × 10⁻⁶, 3.96 × 10⁻¹¹, 9.4 × 10⁻⁹)** | **DO NOT USE — UNTRACEABLE** | — | No stated derivation and not reconstructible from k₁/k₂/k₃/K₁. Using them would import an unexplained saturation/branching assumption. |
| **The "2½ times longer at room temperature" remark** | **NOT A DATUM** | — | No temperature is given. The ~41 kJ/mol figure in §9 is [Z] built on an invented 20 °C and must never enter a parameter, a bound, or an initialisation. It is recorded in this dossier solely so that a later wave does not independently invent it and mistake it for evidence. |

**Circularity risks flagged.**
(i) **k₁ = K₁·k₂ exactly** — the three constants are two independent numbers wearing three hats. Any
fit that treats them as three data points inflates the effective evidence by 50%.
(ii) **k₃ is computed using K₁**, so a "hold-out" on k₃ after fitting K₁ is only partially orthogonal.
Weight it accordingly.
(iii) **Table II's "Calculated" column is the authors' model**, built from the same four constants.
Scoring against it would be scoring against a model, not data.
(iv) **No pH.** Any pH-dependent term in the trunk that is calibrated against this source is
calibrated against an unstated pH. This must be declared, not assumed.
(v) **No Ea, single temperature.** The declaration's existing `rate_transfer: not_licensed to food T
without the Ea` guard is exactly correct and this extraction supplies no basis for lifting it.
(vi) **n = 1 donor, no error bars.** The uncertainty on every number in this paper is unknown; a
prior built on it should be given a deliberately wide sigma rather than a nominal one.
