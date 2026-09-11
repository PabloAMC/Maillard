# Gökmen, Kocadağlı, Göncüoğlu & Ataç Mogol 2012 — EXTRACTION (asparagine + glucose and asparagine + HMF, 0.01 mmol each on 500 mg silica gel, sealed tube, oil bath at 90/120/150/180 C for 5-60 min; acrylamide, 3-APA and HMF by LC-MS/MS, plus Orbitrap confirmation of eleven intermediates)

### HMF LOSS IN AN AMINE POT, MEASURED AT LAST — **"Approximately 67 % of the HMF was lost within 10 min"** at 180 C with equimolar asparagine, which is a pseudo-first-order **k ≈ 0.111 /min and a half-life of ≈ 6.3 min (mine)** against the repository's shipped `k_hmf_self` = 8.97e-7 /min with Ea zero by declaration — a gap of about **1.2e5-fold at 180 C (mine)**. But it is **ONE temperature for the HMF measurement, so NO activation energy for HMF degradation is derivable from this paper**, and the loss is the AMINE sink, not self-degradation.

**Source on disk:** `data/articles/Gokmen2012.pdf` (7 pp., *Food Chemistry* **132** (2012) 168-174).
The `pdftotext -layout` text layer is a clean digital layer and the running prose came through intact.
**THIS PAPER CONTAINS NO NUMBERED TABLES AT ALL.** Its entire quantitative content is (i) numbers
stated in the running text and (ii) six figures. Figures 1a/1b (HMF formation and degradation),
2a/2b (acrylamide kinetics and temperature effect), 3 (Arrhenius plots), 5 (Orbitrap spectra) and 6
(3-APA) are **FIGURE-ONLY**; Figure 4 is a mechanism scheme. **No value below is read off any
plot.** Figure 1b was rastered (`scratchpad/img/gk12f1b.png`) only to confirm its axes and point
spacing, which are recorded qualitatively in section 3 and contribute no number.

## 0. Identity

| field | value |
|---|---|
| Title | "Model studies on the role of 5-hydroxymethyl-2-furfural in acrylamide formation from asparagine" |
| Authors | **Vural Gökmen** (corresponding, vgokmen@hacettepe.edu.tr), **Tolgahan Kocadağlı**, **Neslihan Göncüoğlu**, **Burçe Ataç Mogol** — Department of Food Engineering and Food Research Center, Hacettepe University, 06800 Beytepe, Ankara, Turkey |
| Venue | *Food Chemistry* **132** (2012) 168-174. Received 20 June 2011; revised 28 July 2011; accepted 12 October 2011; available online 19 October 2011 |
| **DOI, exactly as printed in the PDF** | **`doi:10.1016/j.foodchem.2011.10.048`** (printed at the foot of page 168, in that lower-case `doi:` form) |
| Abbreviations as the paper uses them | **ASN** = L-asparagine; **GLC** = D-glucose; **HMF** = 5-hydroxymethyl-2-furfural; **3-APA** = 3-aminopropionamide |
| Companions already on disk | **`hamzalioglu2018_extraction.md`** — the same laboratory (Gökmen), the source of the shipped `k_hmf_self`; **`goncuoglu2016_extraction.md`** — the roasted-hazelnut second laboratory whose step k26 gives the 23 000x disagreement; **`kocadagli2016*_extraction.md`**, **`goncuoglutas2016/2017`**, **`goncuoglu2026`** — the same Hacettepe group. **Three of this paper's four authors appear elsewhere in the corpus as first authors.** See flag 9 |
| Repo status before this dossier | **Not cited anywhere in `src/kinetic_core/`.** No extraction dossier existed |

## 1. Why it matters

**The repository's furanic channel prints its own alarm about a hole, and this paper is a
measurement inside that hole.** `src/kinetic_core/parameters_furanic.py` carries

```
_HMF_SELF_DEGRADATION_FRACTION_7D_5C = 0.009
_HMF_SELF_DEGRADATION_K_PER_MIN = -ln(1 - 0.009) / (7*24*60)      # = 8.97e-7 /min
```

with `flags=("single_temperature_no_ea_licensed", "ea_zero_by_declaration",
"negligible_at_cooking_temperature")` and a note that reads, in the module's own capitals, "THE
MODEL HAS NO EFFECTIVE HMF SINK AT COOKING TEMPERATURE and must therefore be expected to
OVER-PREDICT HMF … THE 50-150 C WINDOW IS EMPTY." `parameters_dicarbonyl.py` sizes the
disagreement against the roasted-hazelnut laboratory at "**23 000x APART** … Half-life at 160 C:
about 33 minutes measured against about 1.5 years shipped."

**This paper puts a number inside the empty window's upper edge.** It measures HMF *disappearance*
in a pot that also contains an amino acid, at **180 C**, and prints the loss: **"Approximately 67 %
of the HMF was lost within 10 min"**. Treating that as pseudo-first-order — which the paper's own
sentence "HMF content decreased exponentially" licenses — gives **k ≈ 0.111 /min and t½ ≈ 6.3 min
(both mine)**. The shipped constant, with Ea zero by declaration, predicts **8.97e-7 /min at every
temperature, i.e. a half-life of about 1.5 years**. The ratio is **≈ 1.24e5 (mine)** — five orders
of magnitude, an order of magnitude worse than the hazelnut laboratory's 23 000x, and at a
temperature 20 C higher.

**And it lands next to the hazelnut number in a way that is worth reporting.**
`parameters_dicarbonyl.py` records the second laboratory at "12 / 21 / 103 x 1e-3 /min at 150 / 160 /
170 C". This paper's 180 C value, **0.111 /min (mine)**, sits directly above that laboratory's
170 C value of **0.103 /min**, in the right order and by a plausible step. **Two independent
laboratories, two entirely different matrices (a roasting hazelnut and a silica-supported
asparagine melt), and the HMF sink comes out within ~8 % across a 10 C gap.** That is the strongest
convergence anywhere in the corpus on this constant, and it makes the shipped value's isolation
much harder to defend. It does not, by itself, license a refit — see flags 2, 3 and 6.

**But three things this paper is NOT, stated up front.**

1. **It is NOT `k_hmf_self`.** The repository's constant is HMF *alone*, no amino acid, aqueous
   pH 3.5, from Hamzalıoğlu & Gökmen 2018's model-free control. **This paper's pot contains
   equimolar asparagine**, and the paper's whole argument is that the loss *is* the reaction with
   asparagine ("indicating its reaction with ASN to form certain Maillard reaction products"). So
   the number belongs to the **amine sink**, which `parameters_furanic.py`'s own note already warns
   is not independent of the self sink. Installing 0.111 /min as `k_hmf_self` would be a category
   error; it is an upper bound on the self sink and a direct measurement of the combined sink.
2. **It gives NO activation energy for HMF degradation.** Methods, verbatim: "Degradation of HMF was
   also determined in ASN–HMF model system", and Fig. 1's caption fixes the condition — "during
   heating at **180 C**". The four-temperature series (90, 120, 150, 180 C) was run for **acrylamide
   formation only**. **One temperature for HMF. No Ea. The brief's question is answered: no.**
3. **The one Ea it does print is for a different quantity.** Acrylamide formation in the ASN–HMF
   system has **Ea = 138.78 kJ/mol** over 90-180 C, and in ASN–GLC **83.94 kJ/mol**. Those are
   product-appearance activation energies, not HMF-disappearance ones, and equating them requires
   assuming the same rate-limiting step. Section 4 carries the 138.78 as a `derived_assumption`
   with that assumption written on it, and flag 5 argues against using it.

**A second thing this paper supplies that the corpus is short of: a fed-intermediate yield for
HMF → acrylamide.** 0.01 mmol HMF plus 0.01 mmol ASN gives **7.31e-4 mmol acrylamide in 60 min at
180 C**, i.e. **7.31 mol % on the HMF (or ASN) charged (mine)** — against **1.62 mol % (mine)** from
the equimolar glucose control and **0.061 mol % (mine)** from asparagine alone. That is a clean
`fed_intermediate_yield` triple, from one pot, one instrument, one heating protocol, with the
no-carbonyl blank included.

## 2. Methods as they matter to a model

- **The pot — a silica-supported low-moisture melt, not a solution.** "A portion of (0.1 ml) of the
  solution containing **0.01 mmol of ASN alone**, a binary mixture of **GLC and ASN (0.01 mmol
  each)**, or a binary mixture of **HMF and ASN (0.01 mmol each)** was transferred to a glass tube
  containing **50 mg of silica gel**. Then, **450 mg of silica gel was added to cover the reaction
  mixture**, and the tube was **tightly closed with a screw cap**."
- **Heating.** "The reactions were performed in an **oil bath at 90, 120, 150, and 180 C for 5, 10,
  20, 30, and 60 min** in order to obtain kinetic and thermodynamic data for acrylamide formation."
- **Which measurement got which temperatures — the crux.** "In addition, **HMF formation was
  determined in both GLC and ASN–GLC model systems** … **during heating at 180 C. Degradation of
  HMF was also determined in ASN–HMF model system.**" Figure 1's caption: "Formation and degradation
  of HMF in the model systems **during heating at 180 C**". **So: acrylamide at four temperatures;
  HMF formation at one; HMF degradation at one.**
- **Replication.** "All reactions were performed **in triplicate**, and **mean values were
  reported**." Error bars appear on the figures. **No standard deviation is printed in the text for
  any value.**
- **Atmosphere, pH, water activity.** A tightly capped tube; **no headspace gas is specified, no pH
  is measured or controlled, and no water activity or moisture content is reported.** The paper
  calls the regime "low moisture conditions" and quantifies it nowhere. 0.1 mL of aqueous solution
  on 500 mg silica is the only handle, and the water is not removed before heating — it is heated in
  a sealed tube, so the tube contains its own steam at 180 C. **This is not a defined water
  activity.**
- **Effective concentration.** 0.01 mmol in 0.1 mL of the applied solution is **0.1 M** at the
  moment of application; after distribution over 500 mg of silica in a sealed tube at 180 C, no
  concentration is definable. **Every rate below is a rate of a supported melt, not of a solution**,
  and cannot be converted to a bimolecular constant.
- **Analysis of HMF.** Extracted with 10 mL of 10 mM formic acid, vortex 2 min, centrifuge 11 180 g
  for 5 min, 0.45 µm nylon filter; Shimadzu UFLC with DAD, Atlantis dC18 250 x 4.6 mm 5 µm, isocratic
  10 mM aqueous formic acid : acetonitrile **90:10 v/v**, 1.0 mL/min, 25 C, **detection at 285 nm**.
  "The concentration of HMF was calculated by means of a **calibration curve built in the range
  between 0.5 and 10 µg/ml**."
- **Analysis of acrylamide.** Same extraction, then Oasis MCX SPE clean-up (first eight drops
  discarded); Waters Acquity H Class UPLC + TQ detector, ESI positive, HSS T3 column, 10 mM formic
  acid with 0.5 % methanol, 0.3 mL/min, 40 C; **MRM 72 → 55 (CE 9 V) and 72 → 44 (CE 12 V)**, dwell
  0.2 s. Calibration **1.0-100 ng/mL**.
- **Analysis of 3-APA.** Agilent 1200 HPLC + 6130 MS, ESI positive, SIM; Atlantis T3, 10 mM formic
  acid : methanol **70:30 v/v**, 0.8 mL/min, 40 C; **m/z 89 [M+H]+ for quantification, m/z 72 for
  confirmation**; retention time **3.1 min**; calibration **10-100 ng/mL**.
- **Structure confirmation.** Thermo Accela LC + **Exactive Orbitrap**, APCI positive (chosen over
  ESI for sensitivity), Atlantis T3, 0.05 % formic acid : methanol 70:30, 0.5 mL/min, 30 C, **full
  scan m/z 50-300 at R = 100 000**, AGC 5e5, max injection 100 ms. **APCI, not ESI — a detail that
  matters for the paper's own caveat about response factors (flag 8).**
- **What the paper measured but did not print.** "The change of the concentrations of ASN and GLC
  were monitored … (**data not shown**)". So the asparagine and glucose time courses exist and are
  unavailable; only the qualitative statement survives: "**No remarkable GLC and ASN remained in the
  model system of ASN–GLC after 5 min at 180 C**", and "The rate of the change of ASN was slightly
  lower in the model system of ASN–HMF than that of ASN–GLC. **The ASN–HMF model system still
  contained available reactants after 10 min of heating at 180 C.**"

## 3. Tables re-typed

**THERE ARE NO TABLES IN THIS PAPER.** What follows is every printed number, transcribed from the
running text with its page anchor, in the paper's own order. Marks: `[M]` measured, `[C]` cited,
`[F]` fitted.

### 3.1 HMF formation and degradation (p. 170, right column; Fig. 1)

| statement, verbatim | value | mark |
|---|---|---|
| "HMF content of the model system containing 0.01 mmol of GLC **linearly increased to 3.29 x 10^-5 mmol within a reaction time of 30 min at 180 C**" | 3.29e-5 mmol | `[M]` |
| "HMF content increased to **2.01 x 10^-5 mmol within 5 min at 180 C**, and remained relatively stable afterward (Fig. 1a)" (ASN–GLC system) | 2.01e-5 mmol | `[M]` |
| "An increase of the **initial rate of HMF formation (5.21 times)** when GLC was heated with ASN" | 5.21x | `[F]` (a ratio of two initial rates, neither printed) |
| "apparent maximum level of HMF was found to be **1.5 times higher in the model system of GLC** than that of ASN–GLC at 180 C" | 1.5x | `[F]` |
| **"The results revealed that HMF content decreased exponentially in the model system of ASN–HMF during heating at 180 C (Fig. 1b). Approximately 67 % of the HMF was lost within 10 min indicating its reaction with ASN to form certain Maillard reaction products."** | **67 % lost in 10 min at 180 C** | **`[M]` — THE HEADLINE NUMBER** |

**Figure 1b, described without extracting values.** A single monotonically decreasing curve with
symbols at **0, 10, 30 and 60 min** and vertical error bars; ordinate labelled "HMF (mmol)" with
printed ticks at 0.0E+00 through 1.2E-02; abscissa "Reaction time (min)", 0 to 60. The 0-min point
sits at the 1.0E-02 tick, consistent with the stated 0.01 mmol charge. **The curve does not reach
zero by 60 min.** No data labels are printed and **no value is read from it.** The 10-min point is
the one the 67 % sentence quantifies.

### 3.2 Acrylamide formation (p. 171, left column; Fig. 2a)

| statement, verbatim | value | mark |
|---|---|---|
| "acrylamide content reached **6.14 x 10^-6 mmol within 60 min of heating at 180 C** in the model system containing **ASN alone**" | 6.14e-6 mmol | `[M]` |
| "In ASN–GLC model system, acrylamide content rapidly reached to an **apparent maximum of 1.62 x 10^-4 mmol within 10 min at 180 C**. It remained relatively stable over the heating period of 10-30 min, and slightly decreased afterward" | 1.62e-4 mmol | `[M]` |
| "acrylamide content exponentially increased to **7.31 x 10^-4 mmol in the ASN–HMF model system within 60 min at 180 C**" | 7.31e-4 mmol | `[M]` |
| "The amount of acrylamide formed in the ASN–HMF model system was **2.27, 4.14, and 5.06 times higher** than in the ASN–GLC model system after **10, 30, and 60 min**, respectively" | 2.27 / 4.14 / 5.06 | `[F]` |

### 3.3 Temperature dependence of acrylamide formation (p. 171, right column; Figs. 2b, 3)

| statement, verbatim | value | mark |
|---|---|---|
| "There was a **linear increase of acrylamide with temperature** in the model system of **ASN–GLC**, while acrylamide content **exponentially increased with temperature** in the model system of **ASN–HMF** (Fig. 2b)" | — | `[M]` |
| "the model system of **ASN–GLC generated more acrylamide** during heating at temperatures **equal to 120 C or lower**. However, the model system of **ASN–HMF generated larger amounts** of acrylamide at temperatures **exceeding 120 C**" | crossover at **≈ 120 C** | `[M]` |
| "acrylamide formation obey the Arrhenius law with **very high correlation coefficients** in a temperature range of **90-180 C**" | (no R² is printed) | `[F]` |
| **"The activation energy of acrylamide formation was found to be 83.94 kJ/mol, and 138.78 kJ/mol for the model systems of ASN–GLC and ASN–HMF, respectively."** | **83.94 and 138.78 kJ/mol** | **`[F]`** |

Fig. 2b's caption: "Effect of temperature on acrylamide formation in ASN–GLC and ASN–HMF model
systems for **heating time of 60 min**." Fig. 3's caption: "Arrhenius plots of acrylamide formation
in the ASN–GLC and ASN–HMF model systems." **Both are figure-only; the two Ea values are printed in
the text and are transcribed as printed.** No pre-exponential factor, no k value at any temperature,
and no R² is printed anywhere.

### 3.4 3-APA (p. 172-173; Fig. 6)

| statement, verbatim | value | mark |
|---|---|---|
| "The amount of 3-APA formed in the model system of **ASN–HMF was 5.30 and 6.37 times higher** than that of GLC–ASN after heating at 180 C for **5 and 30 min**, respectively" | 5.30 / 6.37 | `[F]` |
| "Increasing the heating time from **5 to 30 min at 180 C significantly decreased the amount of 3-APA** in the model systems of ASN–HMF and ASN–GLC, while the amount of acrylamide increased in both" | — | `[M]` |
| **No absolute 3-APA amount is printed anywhere in the paper.** Fig. 6 is figure-only | — | — |

### 3.5 The intermediates confirmed by Orbitrap (Fig. 4 scheme, Fig. 5 spectra; p. 171-172)

Molecular weights as printed in Scheme Fig. 4, with the compound numbers the paper assigns.
**These are structural identifications, not quantities** — the paper prints no concentration for any
of them.

| # | identity as the scheme labels it | MW printed |
|---|---|---:|
| — | glucose | 180 |
| 1 | first glucose dehydration product (−H2O) | 162 |
| 2 | second dehydration product (−2 H2O) | 144 |
| 3 | third dehydration product (−3 H2O) = **HMF** | 126 |
| 4 | Schiff base of glucose with ASN | 294 |
| 5 | Schiff base of compound 1 with ASN | 276 |
| 6 | Schiff base of compound 2 with ASN | 258 |
| 7 | Schiff base of compound 3 (HMF) with ASN | 240 |
| 8 | decarboxylated Schiff base (azomethine ylide) from 4 | 250 |
| 9 | azomethine ylide from 5 | 232 |
| 10 | azomethine ylide from 6 | 214 |
| 11 | azomethine ylide from 7 | 196 |
| 12 | **3-aminopropionamide (3-APA)** | 88 |
| — | **acrylamide** | 71 |

Printed mass accuracies: Schiff bases "**Δ < 1.5 ppm**"; azomethine ylides "**Δ < 0.5 ppm**";
the β-elimination ions at theoretical [M+H]+ **180.08665** and **126.05496**, "**Δ < 0.5 ppm**".
Fig. 5's caption fixes the condition: "reaction intermediates and products formed in the ASN–GLC
model system heated at **180 C for 3 min**" — **a 3-min time point that appears nowhere in the
kinetic series (5, 10, 20, 30, 60 min).**

### 3.6 Derived numbers (mine, arithmetic on the printed values — NOT the paper's)

**The pseudo-first-order HMF loss constant.** The paper prints 67 % lost in 10 min and states the
decay is exponential. Assuming first order in HMF:

- **k = −ln(1 − 0.67) / 10 min = 0.1109 /min ≈ 0.111 /min** (mine)
- **t½ = ln 2 / k = 6.25 min** (mine)

**Against the shipped constant** (`_HMF_SELF_DEGRADATION_K_PER_MIN` = 8.97e-7 /min, Ea = 0 by
declaration, therefore the same at 180 C):

- **ratio = 0.1109 / 8.97e-7 = 1.24e5** (mine)
- **shipped half-life = ln 2 / 8.97e-7 min = 7.73e5 min = 537 days ≈ 1.47 years** (mine), against
  6.25 min measured here.

**Against the second laboratory** (`parameters_dicarbonyl.py`: 12 / 21 / 103 x 1e-3 /min at 150 /
160 / 170 C, roasted hazelnut, `goncuoglu2016_extraction.md` step k26):

| source | T (C) | k (/min) |
|---|---:|---:|
| hazelnut, k26 | 150 | 0.012 `[C]` |
| hazelnut, k26 | 160 | 0.021 `[C]` |
| hazelnut, k26 | 170 | 0.103 `[C]` |
| **this paper (mine, from the printed 67 %/10 min)** | **180** | **0.111** |
| shipped `k_hmf_self` (Ea = 0) | any | 8.97e-7 |

**The 180 C value sits 1.08x above the hazelnut 170 C value (mine)** — the right order and a
plausible step, from a different laboratory in a different matrix. It is also, on the hazelnut
laboratory's own 150→170 C trend, *lower* than a naive extrapolation would give, which is what one
would expect if the hazelnut trend (an 8.6x rise over 20 C) is partly matrix-driven.

**Molar yields (mine, all on the 0.01 mmol charge, 180 C):**

| system | product | printed amount | yield on charge (mine) |
|---|---|---|---:|
| GLC alone, 30 min | HMF | 3.29e-5 mmol | **0.329 mol %** |
| ASN–GLC, 5 min | HMF | 2.01e-5 mmol | **0.201 mol %** |
| ASN alone, 60 min | acrylamide | 6.14e-6 mmol | **0.0614 mol %** |
| ASN–GLC, 10 min (apparent max) | acrylamide | 1.62e-4 mmol | **1.62 mol %** |
| **ASN–HMF, 60 min** | **acrylamide** | **7.31e-4 mmol** | **7.31 mol %** |

**Acrylamide per HMF actually consumed (mine, chaining two printed numbers).** At 10 min the paper
gives ASN–GLC = 1.62e-4 mmol and the ASN–HMF/ASN–GLC ratio = 2.27, so **ASN–HMF at 10 min =
3.68e-4 mmol acrylamide**; and 67 % of 0.01 mmol HMF = **6.7e-3 mmol HMF consumed**. So
**acrylamide / HMF consumed = 3.68e-4 / 6.7e-3 = 5.5 mol % at 10 min, 180 C (mine)**. **About 95 %
of the HMF that disappears goes somewhere other than acrylamide**, and the paper identifies none of
it quantitatively.

**Internal consistency check (mine).** At 60 min the printed ASN–HMF value is 7.31e-4 mmol and the
printed ratio is 5.06, implying **ASN–GLC at 60 min = 1.445e-4 mmol** — slightly below its 10-min
apparent maximum of 1.62e-4, exactly matching the prose "slightly decreased afterward". **The
printed numbers are mutually consistent.**

## 4. Numbers the repository can use

**All rows: 0.01 mmol of each reactant applied in 0.1 mL water onto 500 mg silica gel in a
screw-capped glass tube, oil bath, no controlled pH, no measured water activity, triplicate with
means reported and no printed SD.**

### The HMF sink — what the brief asked for

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **HMF remaining after 10 min** | **≈ 33 %** (i.e. "approximately 67 % lost") | % of charge | **180 C**, equimolar ASN present, silica-supported | p. 170 right col. | **measured_rate** (as a fractional loss) |
| **pseudo-first-order HMF loss constant** | **0.111** (mine) | /min | 180 C, +ASN | derived from the 67 % sentence, assuming first order (the paper says "decreased exponentially") | **measured_rate (mine)** — the arithmetic is mine, the 67 % and the exponential form are the paper's |
| **HMF half-life** | **6.25** (mine) | min | 180 C, +ASN | as above | **measured_rate (mine)** |
| ratio to shipped `k_hmf_self` | **1.24e5** (mine) | — | 180 C vs a constant with Ea = 0 | derived | measured_ratio (mine) |
| **HMF loss at any other temperature** | **NOT MEASURED** | — | — | — | — |
| **activation energy for HMF degradation** | **NOT DERIVABLE.** One temperature only (180 C) | — | — | Methods p. 169; Fig. 1 caption | — |
| HMF remaining at 30 and 60 min | **FIGURE-ONLY** (Fig. 1b has points there; no values printed) | — | 180 C, +ASN | Fig. 1b | figure_only |

**Explicitly, because the distinction decides whether this can be installed:** this is the
**HMF + asparagine** sink. `parameters_furanic.k_hmf_self` is the **HMF-alone** sink from
Hamzalıoğlu & Gökmen 2018 (pH 3.5, 5 C, no amino acid). **They are different reactions.** The
honest uses of 0.111 /min are (a) as a **`measured_bound`**: the total HMF sink at 180 C in an
amine-containing pot is at least this fast, so `k_hmf_self` + the amine sink must together be at
least 0.111 /min there; and (b) as the **amine sink's own constant**, if the channel grows one.
Installing it as `k_hmf_self` would silently attribute an amine-driven loss to self-degradation.

### Temperature dependence — for acrylamide, not for HMF

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **Ea, acrylamide formation, ASN–GLC** | **83.94** | kJ/mol | 90-180 C, 60 min, silica-supported | p. 171 right col.; Fig. 3 | **measured_rate** (an Arrhenius fit; the underlying k values are not printed) |
| **Ea, acrylamide formation, ASN–HMF** | **138.78** | kJ/mol | 90-180 C, 60 min | p. 171 right col.; Fig. 3 | **measured_rate** |
| Ea of the HMF + ASN channel, taken as the ASN–HMF acrylamide Ea | 138.78 | kJ/mol | 90-180 C | — | **derived_assumption** — valid only if acrylamide appearance and HMF disappearance share a rate-limiting step. **Flag 5 argues they do not.** Do not install silently |
| crossover temperature, ASN–GLC vs ASN–HMF acrylamide yield | ≈ **120 C** | C | 60 min | p. 171 right col.; Fig. 2b | **measured_bound** — "equal to 120 C or lower" favours GLC, "exceeding 120 C" favours HMF |
| pre-exponential factor, or k at any single temperature | **NOT PRINTED** for either system | — | — | — | — |
| R² of either Arrhenius fit | **NOT PRINTED** ("very high correlation coefficients") | — | — | — | — |

### Fed-intermediate yields — HMF fed as the carbonyl

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **acrylamide from ASN + HMF** | **7.31e-4 mmol = 7.31 mol % on charge (mine)** | mmol / mol % | 0.01 mmol each, **180 C, 60 min** | p. 171 left col. | **fed_intermediate_yield** |
| **acrylamide from ASN + GLC** | 1.62e-4 mmol = 1.62 mol % (mine), apparent maximum | mmol / mol % | 0.01 mmol each, 180 C, **10 min** | p. 171 left col. | **fed_intermediate_yield** — note the different time |
| **acrylamide from ASN alone (blank)** | 6.14e-6 mmol = 0.0614 mol % (mine) | mmol / mol % | 0.01 mmol, 180 C, 60 min | p. 171 left col. | **fed_intermediate_yield** — the no-carbonyl control |
| acrylamide, ASN–HMF / ASN–GLC | **2.27 / 4.14 / 5.06** at 10 / 30 / 60 min | — | 180 C | p. 171 left col. | **within_study_ratio** |
| acrylamide, ASN–HMF / ASN-alone at 60 min | **119x** (mine, 7.31e-4 / 6.14e-6) | — | 180 C | derived | within_study_ratio (mine) |
| **acrylamide per HMF consumed** | **5.5 mol %** at 10 min (mine, chaining 1.62e-4 x 2.27 against 67 % of 0.01 mmol) | mol % | 180 C | derived from p. 170 and p. 171 | **derived_assumption (mine)** — a carbon-accounting bound: ~95 % of consumed HMF is unaccounted |

### HMF formation (the channel's other end)

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| HMF from GLC alone | 3.29e-5 mmol = **0.329 mol % (mine)**, "linearly increased" | mmol | 0.01 mmol GLC, 180 C, 30 min | p. 170 right col. | **fed_intermediate_yield** |
| HMF from ASN + GLC | 2.01e-5 mmol = **0.201 mol % (mine)**, then "relatively stable" | mmol | 0.01 mmol each, 180 C, 5 min | p. 170 right col. | **fed_intermediate_yield** |
| initial-rate enhancement of HMF formation by ASN | **5.21x** | — | 180 C | p. 170 right col. | **within_study_ratio** — the two initial rates themselves are not printed |
| apparent-maximum HMF, GLC / ASN–GLC | **1.5x** | — | 180 C | p. 170 right col. | within_study_ratio |

### 3-APA

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| 3-APA, ASN–HMF / ASN–GLC | **5.30** at 5 min, **6.37** at 30 min | — | 180 C | p. 172-173 | **within_study_ratio** |
| absolute 3-APA at any time | **NOT PRINTED** (Fig. 6 only) | — | — | Fig. 6 | figure_only |
| 3-APA time behaviour | rises to a maximum by 5 min, then **decreases** from 5 to 30 min in both systems while acrylamide rises | — | 180 C | p. 172-173, abstract | **level_only** — the qualitative shape is the paper's central kinetic claim |

### Reactant depletion

| quantity | value | conditions | anchor | evidence class |
|---|---|---|---|---|
| GLC and ASN remaining in ASN–GLC | "**No remarkable GLC and ASN remained** … after 5 min at 180 C" | 180 C | p. 170 right col. | **measured_bound**, qualitative; the time course is "data not shown" |
| ASN in ASN–HMF | "still contained available reactants after **10 min**"; rate of ASN change "slightly lower" than in ASN–GLC | 180 C | p. 170 right col. | level_only; "data not shown" |

### Absent quantities

| quantity | status |
|---|---|
| any HMF rate constant at 90, 120 or 150 C | **NOT PRESENT** |
| any activation energy for HMF degradation or formation | **NOT PRESENT** |
| any pH, water activity, or moisture content | **NOT PRESENT** |
| any standard deviation, printed as a number | **NOT PRESENT** (error bars are figure-only) |
| any concentration of any Fig. 4 intermediate | **NOT PRESENT** — all eleven are identifications only |
| a mass balance on HMF | **NOT PRESENT** — the ~95 % of consumed HMF that is not acrylamide is unaccounted |
| any HMF measurement in an amine-free pot above 25 C | **NOT PRESENT** — the GLC-alone system measures HMF *formation*, not its decay, and no HMF-alone control is heated |

## 5. Flags

1. **The pot is a silica-supported melt in a sealed tube, not a solution and not a food.** 0.01 mmol
   in 0.1 mL applied to 500 mg of silica gel, capped, and dropped into an oil bath at 180 C. The
   water is not removed, so the tube generates its own steam; the paper calls this "low moisture
   conditions" and never measures the moisture. **No concentration, no pH, no water activity and no
   headspace composition is definable.** Every rate here is a rate of that specific support at that
   specific loading, and silica gel is a surface with acidic silanols that can catalyse both
   dehydration and carbonyl-amine condensation. **Transferring 0.111 /min to a dough or a nut is a
   matrix jump the paper does nothing to license.**
2. **The headline number is one printed word: "Approximately".** "Approximately 67 % of the HMF was
   lost within 10 min." There is no SD, no n stated at that point (triplicate is stated globally),
   and no digit beyond two significant figures. **My k = 0.111 /min inherits that: if the true loss
   is 60 % or 75 %, k is 0.092 or 0.139 /min — a ±25 % band from the rounding alone (mine).** The
   1.24e5 comparison survives that band untouched; a refit would not.
3. **Two points, not a curve, define the constant.** The 67 % is a single interval (0 → 10 min). The
   30- and 60-min points exist only in Figure 1b and are not printed, so **first-order behaviour is
   asserted by the word "exponentially" and is not verifiable from the printed record.** Figure 1b's
   visible shape — a steep drop to 10 min then a long shallow tail that does not reach zero — is
   **not** what a single first-order decay looks like; it is what a decay to a non-zero plateau, or
   a two-pool decay, looks like. **If a fraction of the HMF is unreactive or re-formed from the ASN
   adducts, a single-exponential k fitted to the first interval OVERSTATES the long-time sink.**
   This is the strongest argument against installing the number as a rate constant rather than as a
   bound.
4. **It is the amine sink, not the self sink.** Repeated here because it is the single easiest error
   to make with this paper. See section 4.
5. **Do not import 138.78 kJ/mol as the HMF-degradation Ea.** It is the Ea of **acrylamide
   appearance** in the ASN–HMF system, fitted over 90-180 C at a fixed 60 min. Three reasons it is
   not the HMF-loss Ea: (i) only ~5 % of consumed HMF becomes acrylamide (mine, section 3.6), so
   acrylamide tracks a minor branch; (ii) the paper's own kinetics show acrylamide still rising at
   60 min while HMF is largely gone by 10 min, i.e. the two are not rate-coupled at 180 C; (iii)
   3-APA sits between them and peaks at 5 min then falls, so the path is at least three steps and
   the observed Ea is a composite. **The value is real and belongs in an acrylamide row, not a
   furanic one.**
6. **No k values underlie the Arrhenius fits in printed form.** The paper reports two activation
   energies and no pre-exponential factor, no rate constant at any temperature, and no R² ("very
   high correlation coefficients"). **The fits cannot be audited or re-derived**, and Figure 3 is
   figure-only. A four-point Arrhenius fit across 90 C of range in a system whose mechanism the same
   paper says changes across that range (the ASN–GLC/ASN–HMF crossover at ~120 C) is fragile.
7. **The ASN–GLC comparison is not like-for-like on time.** ASN–GLC's 1.62e-4 mmol is an
   **apparent maximum at 10 min**; ASN–HMF's 7.31e-4 mmol is the **60-min value on a still-rising
   curve**. The paper's own printed ratios (2.27 / 4.14 / 5.06 at matched times) are the correct
   comparison; the two headline absolute numbers are not at the same time point.
8. **The paper flags its own response-factor problem, and it applies to the Orbitrap panel.** "it is
   not easy to correlate lower response of potential β-elimination products with their lower
   occurrence in the reaction mixture, because **these forms may have lesser tendencies to ionise
   under the stated APCI condition**." So the Fig. 4/Fig. 5 intermediates are **identifications with
   no quantitative meaning**, including the paper's argument that path II (via 3-APA) dominates path
   I (direct β-elimination) — that argument rests on a signal comparison the authors themselves say
   is not interpretable as an abundance comparison.
9. **Author overlap with the corpus is heavy and creates an independence problem.** Vural Gökmen is
   also the corresponding author of Hamzalıoğlu & Gökmen 2018, from which the shipped `k_hmf_self`
   is taken; Neslihan Göncüoğlu is the Göncüoğlu of `goncuoglu2016_extraction.md`, the roasted-hazelnut
   "second laboratory" whose k26 gives the 23 000x figure; Tolgahan Kocadağlı appears in
   `kocadagli2016*`. **The "two independent laboratories" framing on this constant is weaker than it
   looks — the Hacettepe group is on both sides of it.** The convergence noted in section 1 between
   this paper at 180 C and the hazelnut at 170 C is therefore a *within-group* convergence, not a
   cross-group replication. It still contradicts the shipped value, but it should not be scored as
   independent confirmation.
10. **A 3-min time point appears only in Figure 5's caption.** All Orbitrap identifications were made
    on ASN–GLC heated "at 180 C for **3 min**", a time not in the kinetic series (5, 10, 20, 30, 60).
    So the structural evidence and the kinetic evidence are from different samples.
11. **Nothing here constrains anything below 90 C**, and the HMF degradation measurement constrains
    only 180 C. **The repository's 50-150 C window remains empty after this paper.** What this paper
    does is put a hard number on the far side of that window and make it impossible to argue that
    Ea = 0 is harmless: a constant that is right at 5 C and right at 180 C cannot have Ea = 0 across
    a 1.24e5 gap.
12. **The unaccounted 95 %.** By my own arithmetic, ~95 % of the consumed HMF at 10 min goes to
    something other than acrylamide, and the paper identifies none of it quantitatively — the Schiff
    base (MW 240) and azomethine ylide (MW 196) of HMF with ASN are detected and never quantified.
    Whatever the furanic channel does with this number, **the sink it represents is overwhelmingly a
    sink to unassigned products**, which is what a `FRAG_C`-style route is for. The repository's
    existing `r_hmf_self` edge (`{"HMF": 1} → {"FRAG_C": 6}` in `network.py`) already has the right
    shape for it.
13. **What to request.** (i) **The same ASN–HMF degradation series at 90, 120 and 150 C** — the
    apparatus already ran those temperatures for acrylamide, so the samples may exist; that single
    addition would fill the repository's empty window and license a real Ea. Worth writing to
    Gökmen. (ii) **An HMF-alone control heated at 180 C on the same silica support** — that, and
    only that, would separate `k_hmf_self` from the amine sink at cooking temperature. (iii) The
    "data not shown" ASN and GLC time courses, which would close the mass balance the paper leaves
    open.
