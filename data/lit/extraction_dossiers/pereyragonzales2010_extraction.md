# Pereyra Gonzales, Naranjo, Leiva & Malec 2010 — COMPLETE TRANSCRIPTION
### Maillard kinetics in real skim milk powder: 18 first-order rate constants and 5 activation energies across a_w 0.33–0.98 at 37/50/60 °C.

**Full extraction of every number in `data/articles/pereyragonzales2010.pdf`.**
Wave K4b, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified.

Table 1 re-read from a **400 dpi render** to resolve the unit string, which the text layer
mangles to `k  103 (hs1)`.

---

## 0. PAPER IDENTITY — MATCHES THE EXPECTED IDENTITY. No mis-file.

| field | value as printed |
|---|---|
| Authors | **A. S. Pereyra Gonzales, G. B. Naranjo, G. E. Leiva, L. S. Malec** (corresponding) |
| Title | **"Maillard reaction kinetics in milk powder: Effect of water activity at mild temperatures"** |
| Venue | ***International Dairy Journal* 20 (2010) 40–45** |
| DOI | **10.1016/j.idairyj.2009.07.007** |
| Dates | Received 20 April 2009; revised 9 July 2009; **accepted 20 July 2009**; PDF created Nov 2009; **volume year 2010** |
| Affiliation | Departamento de Química Orgánica, Área Química y Microbiología de Alimentos, **Facultad de Ciencias Exactas y Naturales, Universidad de Buenos Aires**, 1428 Buenos Aires, Argentina |
| Funding | **UBACYT X021** and **UBACYT X062**, Secretaría de Ciencia y Técnica, Universidad de Buenos Aires |
| PDF character | 6 pages, clean digital text layer. **Table 1 fully legible; no cell unreadable.** |

⚠️ **Filename/venue note.** The repo file is `pereyragonzales2010.pdf`; the paper is
`Int. Dairy J.` **20**, 40–45, **2010** (accepted 2009). The expected identity in the wave
brief ("the milk-powder a_w Maillard study, idairyj 2009/2010") is **confirmed exactly.**

**Provenance codes:** **[M]** measured and printed · **[F]** fitted by the authors ·
**[C]** cited · **[Z]** derived by this wave, never printed.

---

## 1. THE ONE-PARAGRAPH ANSWER

**This is the cleanest tabulated kinetic dataset in wave K4b and the only one in a REAL
food.** Table 1 prints **18 first-order rate constants with 95 % confidence limits**
(6 water activities × 3 temperatures) and **5 activation energies with 95 % confidence
limits**, all for the loss of available lysine from commercial skim milk powder. **Every one
of the five printed Ea values reproduces from the paper's own k values by a three-point
Arrhenius regression to within 2.5 %** (§5.1), and the unit string — mangled in the text
layer — is confirmed as **h⁻¹** by an independent arithmetic check against a prose statement
elsewhere in the paper (§3.1). **Nothing is unreadable.** The headline scientific content is
threefold: (i) at 50 and 60 °C the rate is **flat across a_w 0.31–0.71** and falls **4-fold
at a_w 0.98** — i.e. **there is no intermediate-a_w maximum in real milk powder**, unlike the
lactose–casein model system from the same laboratory; (ii) at 37 °C the a_w 0.33 rate
collapses by **56-fold** relative to a_w 0.43, because 37 °C sits **4 °C above the reported
Tg of 33 °C**; and (iii) the a_w 0.33 Arrhenius plot **breaks into three linear zones** with
activation energies of **65.2, 498.0 and 79.9 kJ/mol** — a middle-zone value **3.7× larger
than any physically ordinary Maillard Ea** and confined to a **10 °C window starting 7 °C
above Tg**. ⚠️ **The paper labels its six samples with three mutually inconsistent sets of
water activities** (§3.2).

---

## 2. SYSTEM COMPOSITION AND METHOD — applies to every number below

| variable | value as printed |
|---|---|
| **Substrate** | **Commercial skim milk powder**, local market, **51.7 % w/w lactose and 34.2 % w/w protein, dry basis** |
| Preservative | **0.06 % w/w potassium sorbate** (antimicrobial), Mallinckrodt |
| Pretreatment | **20 % w/w dispersion freeze-dried** — Stokes model 21, **condenser −40 °C, chamber < 100 µm Hg, 48 h** |
| Equilibration | **10 g portions at 25 °C to constant weight in vacuum desiccators** over saturated salts: **MgCl₂ (0.33), K₂CO₃ (0.43), Mg(NO₃)₂ (0.52), KI (0.69), KCl (0.85), K₂SO₄ (0.98)** |
| Sample vessel | **~200 mg sealed into glass flasks immediately after equilibration**, to prevent water exchange; **a_w checked before heat treatment, "no changes were detected due to the transfer"** |
| Temperatures | **37, 50 and 60 °C**; the a_w 0.33 system **additionally at 30, 40 and 45 °C** |
| Sampling | **two flasks removed at regular intervals**; **t = 0 controls analysed in triplicate** |
| Longest run | **35 months at a_w 0.33 / 37 °C** (only the first 11 months are plotted in Fig. 1) |
| **Response variable** | **available lysine**, by the **o-phthaldialdehyde / N-acetyl-L-cysteine (OPA/NAC) spectrophotometric method** [C] Medina Hernández & García Alvarez-Coque 1992 |
| OPA/NAC reagent | **25 mL of 0.05 M ethanolic OPA + 25 mL of 0.05 M aqueous NAC + 200 mL of 0.02 M boric acid–borate buffer, pH 9.5, made to 1 L**; **2.0 mL sample into 10 mL reagent, diluted to 25 mL** |
| Sample dissolution | **5 % w/v sodium dodecylsulphate** |
| Read | **absorbance at 335 nm**, HP 8453 |
| Standard curve | **casein**, per Malec et al. 2002 [C] |
| Replication | **3 replicates per sample × 2 flasks per time point = the mean of six measurements is reported**; **assay CV < 3 %** |
| Nitrogen | Kjeldahl in duplicate, Büchi 430 digestor + 320 distillation unit |
| **Kinetic model** | **pseudo first-order**, `L_t = L₀·e^{kt}`; *"the loss of available lysine content ... followed better linear relationships when ln(100·L_t/L₀) was plotted against time than when 100·L_t/L₀ was plotted"* |
| Statistics | linear regression; **95 % CI by Student t**; **ANOVA** to test whether rate constants differ; **F-test for lack of fit**; **Ea with 95 % CI by the point-by-point analysis of Labuza (1984)** [C] |
| **a_w reporting convention** | *"As a_w may change during storage due to crystallization and/or Maillard reaction, it was reported as **initial a_w**."* |

**⚠️ NO BUFFER, NO pH CONTROL.** The paper contrasts itself with its own earlier lactose–
casein model (Malec et al. 2002), which **did** use phosphate buffer: *"The differences could
be mostly attributed to the use of phosphate buffer in the model system; this avoids the
decrease in pH due to Maillard reaction, and thus the reduction of the reaction rate
throughout the experiment."* **⇒ Every rate constant below is an unbuffered, drifting-pH
rate constant in a real matrix, and that is exactly what makes it valuable and exactly what
makes it non-transferable to a buffered model.**

---

## 3. TABLE 1 — THE 18 RATE CONSTANTS AND 5 ACTIVATION ENERGIES. **Anchor: Table 1, p. 42.** [M]/[F]

Title as printed: *"First-order rate constants (k) and activation energies (E_a) with 95 %
confidence limits for available lysine loss from skim milk powder stored at 37, 50 and 60 °C
and equilibrated at various water activities (a_w)."*
Column headers as printed at 400 dpi: **`a_w`** · **`k × 10³ (hs⁻¹)`** with sub-headers
**`37 °C`, `50 °C`, `60 °C`** · **`E_a (kJ mol⁻¹)`**.

| **a_w** | **k × 10³ at 37 °C** | **k × 10³ at 50 °C** | **k × 10³ at 60 °C** | **E_a (kJ mol⁻¹)** |
|---:|---:|---:|---:|---:|
| **0.33** | **0.0198 ± 0.0020** | **7.14 ± 0.51** | **17.46 ± 0.93** | ⚠️ **blank — Arrhenius non-linear** |
| **0.43** | **1.115 ± 0.075** | **7.24 ± 0.78** | **27.4 ± 1.9** | **121.1 ± 2.7** |
| **0.52** | **0.762 ± 0.095** | **6.35 ± 0.32** | **26.0 ± 1.3** | **129.2 ± 3.7** |
| **0.69** | **0.802 ± 0.046** | **7.08 ± 0.53** | **27.2 ± 1.9** | **131.4 ± 3.1** |
| **0.85** | **0.538 ± 0.051** | **4.56 ± 0.34** | **20.8 ± 1.4** | **135.4 ± 2.9** |
| **0.98** | **0.245 ± 0.020** | **1.89 ± 0.13** | **7.54 ± 0.58** | **128.5 ± 4.1** |

**⇒ Actual k values in h⁻¹ [Z] (divide the printed entries by 1 000):**

| a_w | **k, 37 °C (h⁻¹)** | **k, 50 °C (h⁻¹)** | **k, 60 °C (h⁻¹)** |
|---:|---:|---:|---:|
| 0.33 | **1.98 × 10⁻⁵** | **7.14 × 10⁻³** | **1.746 × 10⁻²** |
| 0.43 | **1.115 × 10⁻³** | **7.24 × 10⁻³** | **2.74 × 10⁻²** |
| 0.52 | **7.62 × 10⁻⁴** | **6.35 × 10⁻³** | **2.60 × 10⁻²** |
| 0.69 | **8.02 × 10⁻⁴** | **7.08 × 10⁻³** | **2.72 × 10⁻²** |
| 0.85 | **5.38 × 10⁻⁴** | **4.56 × 10⁻³** | **2.08 × 10⁻²** |
| 0.98 | **2.45 × 10⁻⁴** | **1.89 × 10⁻³** | **7.54 × 10⁻³** |

### 3.1 ⚠️ RESOLVING THE UNIT — `hs⁻¹` is `h⁻¹`, verified arithmetically

The header prints **`hs⁻¹`**, which is not a unit. Two candidate readings: **h⁻¹** (hours,
with a stray `s`) or **s⁻¹**. **Resolved by the paper's own prose, p. 42:** *"when milk powder
was stored at an initial a_w of 0.47, there was a **20 % loss of lysine after 7 days**."*
The a_w 0.43 row's 37 °C constant is 1.115 × 10⁻³. Over **7 days = 168 h**:
`ln(L/L₀) = −1.115 × 10⁻³ × 168 = −0.187 ⇒ 17 % loss` **[Z]** — against a stated 20 %.
On the `s⁻¹` reading the same 7 days would give complete destruction (`k·t = 605`).
**⇒ The unit is h⁻¹.** Independently corroborated by Fig. 2, whose `ln k` axis runs −12 to
−1: `ln(1.98 × 10⁻⁵) = −10.83` and `ln(2.74 × 10⁻²) = −3.60`, both inside the axis.

### 3.2 ⚠️ THE SIX SAMPLES CARRY THREE MUTUALLY INCONSISTENT SETS OF WATER ACTIVITIES

| source in the paper | the six values |
|---|---|
| **§2.2 (salt list) and Table 1** | **0.33, 0.43, 0.52, 0.69, 0.85, 0.98** |
| **Fig. 1 caption** (*"initial water activities of..."*) | **0.32, 0.47, 0.56, 0.70, 0.85, 0.98** |
| **Running text, p. 42** | *"initial a_w range **0.31–0.71**"*, *"an initial a_w of **0.47**"*, *"at a_w of **0.97** their values were approximately four-fold lower"*, *"the initial a_w of **0.32**"* |

**Every one of the six samples is labelled with two or three different a_w values, differing
by up to 0.05.** The Table 1 column is the nominal saturated-salt value at 25 °C; the Fig. 1
caption is presumably the measured initial a_w. **Neither is identified as such and the
measured values are never tabulated.** ⚠️ **Ingest Table 1's a_w with
`a_w_uncertainty: ±0.05, three_conflicting_labellings_in_source`.** This matters: an a_w
error of 0.05 near the Tg boundary is the difference between glassy and rubbery.

---

## 4. THE a_w DEPENDENCE — the paper's findings, complete

### 4.1 At 50 and 60 °C: **flat to a_w 0.71, then a fall** (p. 42)

> *"The curves at 50 and 60 °C showed similar trends: the rate constants of samples stored in
> the initial a_w range **0.31–0.71 were not significantly different (P > 0.05)**. Therefore,
> at 50 and 60 °C, **no maximum values were observed at intermediate a_w**. At an initial a_w
> of 0.85, the rate constants decreased slightly, and at a_w of 0.97 their values were
> approximately **four-fold lower** than at the initial a_w of 0.70."*

**Check [Z]:** 60 °C, a_w 0.69 → 0.98: 27.2/7.54 = **3.61×**. 50 °C: 7.08/1.89 = **3.75×**.
**Both "approximately four-fold". ✅**

Mechanism offered, both [C]: *"the dilution of reactants and the **inhibitory effect of
water, a product of the reaction**"* — Eichner & Karel 1972; van Boekel 2001.

**⚠️ This CONTRADICTS the classical picture the paper itself cites in its introduction:**
*"The maximal rate of NEB reaction was reported to occur in the a_w range of **0.5–0.75**"*
[C] Kaanane & Labuza 1989; Labuza & Saltmarch 1981. And it contradicts the same laboratory's
own model system: *"in the model system at 37 and 50 °C, the largest rate constant value was
observed at an initial a_w of **0.52**, and at 60 °C the rate was maximal at a wider range of
a_w. Instead, **no maximum was observed in milk powder** stored at intermediate a_w."*
**⇒ A real dairy matrix and a compositionally matched lactose–casein model system, from the
same lab under the same conditions, disagree on whether an intermediate-a_w rate maximum
exists at all.** That is a first-class finding for any repo term of the form
`rate = f(a_w)` with a hump: **the hump is a model-system artefact, not a milk-powder
property.** The paper attributes the difference to the model system's **phosphate buffer**
(§2) and to *"the presence of salts and whey proteins"*.

### 4.2 At 37 °C: the glass-transition collapse

> *"At 37 °C, the behaviour at intermediate and high a_w was the same. However, at the
> initial a_w of 0.32, the Maillard reaction occurred at a **much slower rate** than at a_w of
> 0.47–0.70. This substantial reduction in the rate constant at 37 °C **seems not to be
> caused only by the decline in a_w**, since at 50 and 60 °C **no lessening was observed at
> a_w of 0.32**. It should be noted that **37 °C is very close to the Tg reported by Jouppila
> and Roos (1994) for skim milk stored at a_w 0.33 (33 °C)**."*

**Quantified [Z]:** at 37 °C, k(0.33)/k(0.43) = 0.0198/1.115 = **1/56.3 — a 56-fold
collapse.** At 50 °C the same pair is 7.14/7.24 = **0.986 — no effect at all.** At 60 °C,
17.46/27.4 = **0.637**.
**⇒ The a_w 0.33 sample is 56× slower than its neighbour when stored 4 °C above Tg, and
indistinguishable from it when stored 17 °C above Tg.** That is a clean, quantified,
temperature-gated matrix-state effect on a Maillard rate in a real food.

**Cited Tg values [C] Jouppila & Roos 1994:** **Tg = 33 °C at a_w 0.33**; **Tg = 9 °C at
a_w 0.43**. Both for skim milk. *"At a_w 0.43 ... a Tg of 9 °C, much lower than 37 °C. Hence,
at this temperature, the increase in the mobility of the reactants at a_w values above 0.33,
as shown by the dramatic rise in reaction rate, could be mainly attributed to the
**plasticizing effect of water, which lessened the glass transition temperature, rather than
to the influence of water activity**."*
Corroborating [C]: *"Bell (1996) ... concluded that the rate of pigment formation was **more
significantly influenced by the Tg of a material than by its a_w**."*

### 4.3 Lactose crystallisation — measured, and shown NOT to matter (Fig. 3, p. 43) [M]

Loss of adsorbed water was tracked gravimetrically in duplicate:

| condition | RH | result as printed |
|---|---:|---|
| a_w 0.33 at **50 °C** | 32 % | *"the system did not show decreasing water content **over 200 h** of storage. Thus, **no crystallization of lactose occurred** during the experiment."* |
| a_w 0.33 at **60 °C** | 31 % | *"a loss of adsorbed water was detected **after a delay of 40 h**; hence, under these conditions, crystallization of lactose **did not affect the reaction rate**, as it occurred **when more than 50 % lysine was lost**."* |
| a_w 0.43 at **37 °C** | 43 % | *"adsorbed water was lost because of lactose crystallization in a rather short time (**within 23 h**)."* |

Conclusion as printed: *"Apparently, **in this study crystallization of lactose had no effect
on the reactivity of lysine**."* Supported by *"there were no significant differences between
the rate constant values of samples equilibrated at a_w 0.33 and 0.43 at 50 and 60 °C"* and
by [C] Morgan et al. 2005 and Vuataz 1988: *"crystallization immobilizes lactose and thus
less dissolved lactose is available for reaction with proteins."*

**⚠️ A latent inconsistency the paper does not resolve.** Crystallisation happened at
**37 °C / a_w 0.43 within 23 h** but not at **50 °C / a_w 0.33 in 200 h**. The paper says
*"at 50 and 60 °C and at the same RH, faster crystallization could be expected"* — but the
50 and 60 °C crystallisation runs were at **RH 32/31 %**, not 43 %. **The a_w 0.43 system was
never tested for crystallisation at 50 or 60 °C.** The conclusion "crystallisation had no
effect" is therefore established at a_w 0.33 and asserted at a_w 0.43.

Also noted: *"a_w could also rise during storage, as water is a product of Maillard reaction.
However, in this study the early stage of this reaction, which only produces a small amount
of water, was predominant, and its contribution to the overall increase of a_w was not
relevant"* [C] Vuataz 2002.

---

## 5. TEMPERATURE DEPENDENCE — the five Arrhenius fits and the three-zone anomaly

### 5.1 ⚠️ ALL FIVE PRINTED Ea VALUES REPRODUCE FROM TABLE 1'S OWN k VALUES [Z]

Three-point ordinary least squares on `ln k` vs `1/T` (T = 310.15, 323.15, 333.15 K):

| a_w | **Ea recomputed [Z]** | **Ea printed [F]** | deviation |
|---:|---:|---:|---:|
| 0.43 | **119.7 kJ/mol** | 121.1 ± 2.7 | **−1.2 %** |
| 0.52 | **132.1** | 129.2 ± 3.7 | +2.2 % |
| 0.69 | **132.2** | 131.4 ± 3.1 | +0.6 % |
| 0.85 | **136.5** | 135.4 ± 2.9 | +0.8 % |
| 0.98 | **128.2** | 128.5 ± 4.1 | −0.2 % |

**Every value reproduces within 2.5 %.** The small residual is expected: the paper used the
**point-by-point method of Labuza (1984)**, not OLS. **⇒ Table 1 is internally self-consistent
and verified twice over. This is the highest-confidence kinetic block in wave K4b.**

Conclusion as printed: *"The values varied from **121 to 135 kJ mol⁻¹** and were within the
range found by other authors for loss of lysine in milk-related systems"* [C] Kessler & Fink
1986; Malec et al. 2002; Morales et al. 1995; Naranjo et al. 1998.
And: *"in this study, the **water activity did not affect the temperature-dependence** of the
rate of loss of lysine, as Ea values for the systems equilibrated at a_w values of 0.43–0.98
were **not significantly different (P > 0.05)**."*

**⚠️ THIS DIRECTLY CONTRADICTS Miao & Roos 2004** (this wave, `miao2004_extraction.md`),
which the paper cites by name: *"Some authors (Malec et al., 2002; **Miao & Roos, 2004**)
found a higher temperature-dependence of Maillard reaction with decreasing moisture content.
However, in this study, the water activity did not affect the temperature-dependence."*
**Two papers in this wave reach opposite conclusions on the same question.** Miao's Ea range
is **117.8 → 183.6 kJ/mol** over a water-content change of 3.78 → 6.86 g/100 g solids;
Pereyra Gonzales's is **121.1 → 135.4 kJ/mol, non-significant**, over a_w 0.43 → 0.98.
**Both are correctly reported and both are probably right — of different systems.** The
discriminating variables are named in §7.

### 5.2 THE a_w 0.33 THREE-ZONE ARRHENIUS PLOT (Fig. 4, p. 43) — the anomaly

> *"Regarding the Arrhenius plot of the system equilibrated at a_w 0.33, **no linear
> relationship was observed**. In order to analyze the temperature-dependence ... other
> temperatures near Tg were studied (**30, 40 and 45 °C**). The equation was fitted to the
> experimental data **in three parts, with different slopes** (Fig. 4) ... Three zones were
> observed and the three straight lines corresponded to temperature ranges of **30–40 °C,
> 40–50 °C and 50–60 °C**. The activation energies for lysine loss were calculated separately
> for each temperature range and were **65.2, 498.0 and 79.9 kJ mol⁻¹**, respectively."*

| zone | temperature range | **Ea (kJ mol⁻¹)** [F] | relation to Tg = 33 °C |
|---|---|---:|---|
| glassy / near-Tg | **30–40 °C** | **65.2** | spans Tg; **below the usual 85–166 range** |
| **transition** | **40–50 °C** | **498.0** | **break at 7 °C above Tg** |
| rubbery | **50–60 °C** | **79.9** | **second break at 13 °C above Tg** |

**Verification [Z]:** the 50–60 °C zone can be checked against Table 1 directly.
`Ea = R·ln(17.46/7.14)/(1/323.15 − 1/333.15) = 8.314 × 0.8942/9.289 × 10⁻⁵ =` **80.0
kJ/mol** — against the printed **79.9. ✅ Exact to 0.1 %.** The 30–40 and 40–50 zones cannot
be checked because **the k values at 30, 40 and 45 °C are in Figure 4 only and are never
tabulated** — see §6.

Interpretation as printed: *"E_a in the intermediate region (40–50 °C) was considerably
higher than the two others, and **also much higher than typical values for the Maillard
reaction**. This higher temperature-dependence of the reaction rates was explained by the
**abrupt change of the viscosity of the system above the glass transition temperature**. In
the present study, the break was observed at **7 °C above Tg**, which was in agreement with
previous studies of Karmas et al. (1992), Lievonen and Roos (2002) and Roos, Jouppila and
Zielasko (1996), who reported that a large increase in NEB occurred at a range of
**2 °C–40 °C above Tg**. At 30 °C, where the system was in the glassy state, the Maillard
reaction still occurred, although at a very low rate, and the E_a corresponding to the
temperature range in the vicinity of the Tg (30–40 °C) was **lower than the usual values
reported for the loss of lysine (85–166 kJ mol⁻¹)** [C] Kessler & Fink 1986; Labuza &
Saltmarch 1981; Malec et al. 2002; Naranjo et al. 1998. **The second break occurred at 13 °C
above the Tg** and the E_a at the higher temperature zone (50–60 °C) was closer to the range
... indicating that **at temperatures well above Tg, the reaction became diffusion-controlled
and followed Arrhenius kinetics**."*

**⚠️ The 498.0 kJ/mol figure carries no confidence interval, is fitted over a 10 °C window,
and — because the k values at 40 and 45 °C are figure-only — cannot be verified from any
printed number.** It is **3.7×** the largest ordinary value in the paper (135.4) and
**7.6×** the smallest (65.2). **It should never be ingested as a Maillard activation energy.
It is a description of a state transition, not a chemical barrier**, and the paper says as
much ("the abrupt change of the viscosity"). If the repo carries an Ea prior, this value
would corrupt it.

---

## 6. WHAT IS **NOT** IN THE PAPER — the figure-only content

| content | where it lives | recoverable? |
|---|---|---|
| **k at 30, 40 and 45 °C for a_w 0.33** (3 values) | **Figure 4 only** (Arrhenius plot, `ln k` vs `1/T`, axis 0.00295–0.00335, `ln k` −12 to 0) | Digitisable in principle; **NOT digitised in this wave** — the plot is small, six of the nine points overlap the three fitted segments, and the two Ea values that would be checked (65.2, 498.0) come from those points. **Marked `unreadable_this_wave`, not `unreadable`.** |
| **Lysine-loss time courses** (Fig. 1, 6 curves at 37 °C, semilog, 0–350 days) | Figure 1 only | The k values are in Table 1; the raw time courses are not needed. |
| **`ln k` vs a_w curves** for this study **and** the Malec 2002 model system (Fig. 2, 6 series) | Figure 2 only | **The three dotted lines are the Malec et al. 2002 lactose–casein model system** — the comparison dataset for §4.1's contradiction. **Not digitised here; the primary source (Malec et al. 2002) should be retrieved instead** (§9). |
| Water-sorption kinetics (Fig. 3, 3 curves, water content vs time, 0–100 h) | Figure 3 only | Qualitative conclusions transcribed in §4.3. |
| **Measured (as opposed to nominal) a_w of the six samples** | **nowhere** | ⚠️ Never tabulated (§3.2). |
| **pH of the milk powder at any point** | **nowhere** | ⚠️ Never measured. The paper's central explanation for its difference from Malec 2002 is a pH-drift argument, made without a single pH measurement. |
| **Moisture content (g water/100 g solids) of the six samples** | **nowhere** | ⚠️ Only a_w is given. Comparison with Miao & Roos 2004 and Lievonen & Roos 2002, both of which report **water content**, therefore requires an external sorption isotherm. |

---

## 7. RECONCILING THIS PAPER WITH MIAO & ROOS 2004 AND LIEVONEN & ROOS 2002 [Z]

All three are in wave K4b. All three study Maillard kinetics vs water and temperature. They
disagree, and the disagreements are informative rather than contradictory:

| | **Pereyra Gonzales 2010** | **Miao & Roos 2004** | **Lievonen & Roos 2002** |
|---|---|---|---|
| matrix | **real skim milk powder** | lactose/trehalose 1:1 model | maltodextrin or PVP model |
| reactants | **the food's own lactose + protein lysine** | **added D-xylose + L-lysine, 5 % of solids** | **added D-xylose + L-lysine, 10 % of sorbed water** |
| response | **available lysine loss (reactant)** | **optical density 280/420 nm (product)** | **optical density 280/420 nm (product)** |
| order | **pseudo first-order** | **zero-order** | **pseudo-zero-order** |
| buffer | **none** | none | none (a citrate arm tested separately) |
| T range | **30–60 °C** | **40–90 °C** | **10–110 °C** |
| **does Ea depend on water?** | **NO** (121–135, N.S., over a_w 0.43–0.98) | **YES** (117.8 → 183.6 over 3.78 → 6.86 g/100 g) | not computed |
| **is there an intermediate-a_w maximum?** | **NO** in milk powder; **YES** in the matched model system | rate rises monotonically with water | rate rises monotonically with a_w |

**Three candidate reconcilers, none of which the papers themselves state:**
1. **Different response variables.** Pereyra Gonzales measures **loss of a reactant** in the
   *initial* stage; Miao and Lievonen measure **appearance of a product** (colour). These are
   different steps of the cascade and there is no reason their activation energies should
   agree. ⚠️ **The repo must not pool lysine-loss Ea with browning Ea.**
2. **Different water axes.** Pereyra Gonzales varies **a_w at fixed composition**; Miao varies
   **water content at fixed a_w-to-Tg relationship**. Miao's Ea trend is with **water content**,
   Pereyra Gonzales's null is over **a_w**. In a milk powder the two are not interchangeable
   and this paper never reports water content (§6).
3. **Different distance from Tg.** Miao's three water contents have Tg = 55.0, 34.8, 26.4 °C
   against a 40–90 °C range, i.e. **all rubbery**. Pereyra Gonzales's a_w 0.43–0.98 samples
   have Tg ≤ 9 °C against 37–60 °C, i.e. **all well rubbery** — which is precisely the regime
   in which this paper finds **no** Ea dependence, and precisely the regime in which its
   a_w 0.33 sample (Tg 33 °C) shows the **largest** effect in the whole paper. **The two
   results are compatible if the Ea dependence lives in the transition region, not in the
   rubbery bulk.**

**⇒ Proposed corpus rule, for the orchestrator's consideration: `Ea_depends_on_water` should
be scoped to `T − Tg < ~20 °C`, and treated as absent well above Tg.** Both papers' data are
consistent with that; neither states it.

---

## 8. CONSOLIDATED NEW-PARAMETER TABLE

**Common conditions:** commercial skim milk powder (51.7 % lactose, 34.2 % protein, dry
basis) + 0.06 % potassium sorbate, 20 % dispersion freeze-dried, equilibrated at 25 °C over
saturated salts, ~200 mg sealed in glass flasks, **unbuffered, pH uncontrolled and never
measured**, available lysine by OPA/NAC at 335 nm, **pseudo-first-order**, n = 6 measurements
per point, assay CV < 3 %.

| # | parameter | value | units | conditions | class | anchor |
|---:|---|---:|---|---|:--:|---|
| 1–6 | **k, lysine loss, 37 °C** at a_w 0.33 / 0.43 / 0.52 / 0.69 / 0.85 / 0.98 | **1.98e-5 ± 2.0e-6 / 1.115e-3 ± 7.5e-5 / 7.62e-4 ± 9.5e-5 / 8.02e-4 ± 4.6e-5 / 5.38e-4 ± 5.1e-5 / 2.45e-4 ± 2.0e-5** | **h⁻¹** (95 % CI) | 37 °C, skim milk powder | M | T1 p.42 |
| 7–12 | **k, lysine loss, 50 °C**, same a_w order | **7.14e-3 ± 5.1e-4 / 7.24e-3 ± 7.8e-4 / 6.35e-3 ± 3.2e-4 / 7.08e-3 ± 5.3e-4 / 4.56e-3 ± 3.4e-4 / 1.89e-3 ± 1.3e-4** | **h⁻¹** | 50 °C | M | T1 |
| 13–18 | **k, lysine loss, 60 °C**, same a_w order | **1.746e-2 ± 9.3e-4 / 2.74e-2 ± 1.9e-3 / 2.60e-2 ± 1.3e-3 / 2.72e-2 ± 1.9e-3 / 2.08e-2 ± 1.4e-3 / 7.54e-3 ± 5.8e-4** | **h⁻¹** | 60 °C | M | T1 |
| 19–23 | **E_a, lysine loss** at a_w 0.43 / 0.52 / 0.69 / 0.85 / 0.98 | **121.1 ± 2.7 / 129.2 ± 3.7 / 131.4 ± 3.1 / 135.4 ± 2.9 / 128.5 ± 4.1** | kJ mol⁻¹ (95 % CI) | 37–60 °C, 3-point Arrhenius | F | T1 |
| 24 | **E_a is independent of a_w over 0.43–0.98** | **P > 0.05** | — | 37–60 °C | M | p.43 |
| 25 | **E_a, a_w 0.33, zone 30–40 °C** | **65.2** | kJ mol⁻¹ (**no CI printed**) | spans Tg = 33 °C | F | p.44 |
| 26 | ⚠️ **E_a, a_w 0.33, zone 40–50 °C** | **498.0** | kJ mol⁻¹ (**no CI, unverifiable, DO NOT INGEST as a chemical Ea**) | break at Tg + 7 °C | F | p.44 |
| 27 | **E_a, a_w 0.33, zone 50–60 °C** | **79.9** | kJ mol⁻¹ | break at Tg + 13 °C; **reproduces to 0.1 % from T1** | F+Z | p.44 |
| 28 | **rate ratio a_w 0.43/0.33 at 37 °C** | **56.3** | × | 4 °C above Tg | Z | §4.2 |
| 29 | rate ratio a_w 0.43/0.33 at 50 °C | **0.99 (no effect)** | × | 17 °C above Tg | Z | §4.2 |
| 30 | rate ratio a_w 0.43/0.33 at 60 °C | **0.64** | × | 27 °C above Tg | Z | §4.2 |
| 31 | **dilution effect, a_w 0.69 → 0.98** | **3.61× (60 °C), 3.75× (50 °C)** rate reduction | × | independent of temperature | Z | §4.1 |
| 32 | **a_w range over which rate is flat** | **0.31–0.71, P > 0.05** | — | 50 and 60 °C | M | p.42 |
| 33 | **no intermediate-a_w rate maximum in real milk powder** | qualitative, but the contrast with the matched model system is the finding | — | 37–60 °C | M | p.43 |
| 34 | Tg of skim milk at a_w 0.33 | **33** | °C | cited | C | Jouppila & Roos 1994 via p.42 |
| 35 | Tg of skim milk at a_w 0.43 | **9** | °C | cited | C | " |
| 36 | crystallisation onset, a_w 0.43, 37 °C | **within 23 h** | h | RH 43 % | M | Fig.3 / p.43 |
| 37 | crystallisation onset, a_w 0.33, 60 °C | **40 h** (after >50 % lysine loss) | h | RH 31 % | M | " |
| 38 | crystallisation, a_w 0.33, 50 °C | **none in 200 h** | — | RH 32 % | M | " |
| 39 | literature E_a range for lysine loss | **85–166** | kJ mol⁻¹ | cited | C | p.44 |
| 40 | literature a_w of maximal NEB | **0.5–0.75** | a_w | cited — **NOT observed here** | C | p.40 |
| 41 | literature a_w range of maximal browning, milk powder 70–130 °C | **0.44–0.85** | a_w | cited | C | Acevedo 2006, Franzen 1990 via p.42 |
| 42 | break in NEB rate above Tg | **2–40 °C above Tg** (lit.); **7 and 13 °C above Tg** (this study) | °C | — | C + F | p.44 |
| 43 | milk powder composition | **51.7 % lactose, 34.2 % protein, dry basis** | % w/w | commercial skim milk powder | M | p.41 |

---

## 9. PROPOSED FIT / HOLD-OUT ROLE — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **Proposal only.** Pereyra Gonzales et al. 2010 is a **new source**; a declaration
> amendment must be approved before any wave fits any row. This dossier does not edit the
> declaration.

| dataset | rows | **proposed role** | rationale |
|---|---:|---|---|
| **Table 1, 15 rate constants at a_w 0.43–0.98** | 15 | **FIT-ELIGIBLE — the best-quality candidate in wave K4b** | Tabulated, with 95 % CIs, internally verified twice (§3.1 unit check, §5.1 Ea reproduction), in a **real food matrix**, with a **reactant-loss** response variable rather than a colour proxy. If the repo's Maillard lane wants one new fit target this wave, this is it. |
| **Table 1, 3 rate constants at a_w 0.33** | 3 | **HOLD-OUT** | These three are the glass-transition probe (§4.2, 56× collapse at 37 °C). They are the only rows in the paper that test a matrix-state effect and should adjudicate, not calibrate. |
| **Table 1, 5 activation energies** | 5 | **FIT-ELIGIBLE, but redundant** | They are derived from rows 1–18 (§5.1 shows they reproduce). **Fitting both k and Ea double-counts.** Recommend: fit the k values, use the Ea values as a consistency check. |
| **The "E_a independent of a_w" null (P > 0.05)** | 1 | **HOLD-OUT — falsification test** | Directly contradicts Miao & Roos 2004 (this wave). Any repo term making Ea a function of water must reproduce **both** results, which §7 argues requires scoping by `T − Tg`. |
| **§5.2 three-zone Arrhenius, a_w 0.33** | 3 | **HOLD-OUT for the 65.2 and 79.9 values; REJECT the 498.0** | 79.9 is verified to 0.1 %; 65.2 and 498.0 rest on figure-only k values with no CI. **498.0 kJ/mol must never enter an Ea prior** (§5.2). |
| **§4.1 flat-then-fall a_w profile** | 1 | **HOLD-OUT — acceptance criterion** | Proposed gate: a repo `f(a_w)` term must be **flat within error over a_w 0.31–0.71 and 3.6–3.8× lower at a_w 0.98**, at 50 and 60 °C, in a milk-powder matrix. It must **not** produce a hump — the hump exists only in the buffered model system. |
| **§4.3 crystallisation timings** | 3 | **METADATA** | Establishes that crystallisation is not confounding rows 1–18. Also a bound on how long a lactose-bearing powder stays amorphous. |
| **a_w labels** | 6 | **INGEST with `a_w_uncertainty: ±0.05`** | Three conflicting labellings in the source (§3.2). |
| **Figure-only k at 30, 40, 45 °C** | 3 | **NOT EXTRACTED — flagged `unreadable_this_wave`** | Would verify the 65.2 and 498.0 values. Digitisation deferred (§6). |
| **Fig. 2 dotted lines (Malec 2002 model system)** | ~18 | **DO NOT DIGITISE — retrieve the primary instead** | See §10.1. |

---

## 10. RETRIEVALS THIS PAPER MAKES WORTH REQUESTING

1. **Malec, L. S., Pereyra Gonzales, A. S., Naranjo, G. B. & Vigo, M. S. (2002)**,
   *Food Res. Int.* **35**, 849–853 — *"Influence of water activity and storage temperature
   on lysine availability of a milk like system."* **Same laboratory, same temperatures, same
   a_w values, same response variable, same kinetic order — but a phosphate-buffered
   lactose–casein model instead of a real milk powder.** This is a **near-perfect matched
   pair** and the source of the §4.1 contradiction (model system shows an intermediate-a_w
   maximum, real milk powder does not). **Retrieving it would give the corpus a controlled
   model-vs-real-food comparison on Maillard kinetics that nothing else in the corpus
   provides.** Highest-value retrieval from this paper by a wide margin.
2. **Jouppila, K. & Roos, Y. H. (1994)**, *J. Dairy Sci.* **77**, 2907–2915 — *"Glass
   transition and crystallization in milk powders."* Source of the Tg = 33 °C (a_w 0.33) and
   Tg = 9 °C (a_w 0.43) values on which §4.2's entire interpretation rests, **and** the source
   of a full Tg-vs-a_w curve for skim milk powder, which the repo would need to apply any
   `T − Tg` scoping rule (§7) to a dairy matrix.
3. **Vuataz, G. (2002)**, *Lait* **82**, 485–500 — *"The phase diagram of milk: a new tool for
   optimising."* A published **state diagram for milk**, i.e. Tg, crystallisation and
   solubility boundaries in one place. If the repo ever needs a dairy matrix-state model, this
   is the reference object.
4. **Kessler, H.-G. & Fink, R. (1986)**, *J. Food Sci.* **51**, 1105–1111, 1155 — the classical
   source of the 85–166 kJ/mol lysine-loss Ea range that §5.2 uses to judge the 65.2 and
   498.0 outliers.
