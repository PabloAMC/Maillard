# Wang & Ho 2008 (`10.1021/jf8012025`) — single-paper extraction, Wave K5b, 2026-08-29

> Wang, Y.; Ho, C.-T.*
> "**Formation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone through Methylglyoxal: A Maillard Reaction
> Intermediate**"
> *J. Agric. Food Chem.* **2008**, 56 (16), 7405−7409.
> Department of Food Science, Rutgers University, New Brunswick, New Jersey 08901-8520.
> Corresponding author: C.-T. Ho, 65 Dudley Rd., New Brunswick, NJ; ho@aesop.rutgers.edu.
> Received for review April 16, 2008. Revised May 29, 2008. Accepted June 1, 2008.
> Published on Web **July 2, 2008**. Article ID `JF8012025`.

**Source PDF:** `data/articles/wang2008.pdf` (5 pp., 242,757 bytes, on disk 2026-08-29 09:41).
**Read method:** full text layer, `pdftotext -layout`. **⚠️ THE PAPER HAS NO TABLES.** Its entire
quantitative content is **Figure 1, a bar chart with no text layer** (verified with
`pdftotext -bbox`: page 2 carries the running head at y = 22 pt and then nothing until the figure
caption at y = 193 pt — the plot area is vector line-art with no glyphs). Figure 1 was therefore
**digitised from a 300 dpi render by pixel measurement**, with the axis calibration verified against
the tick marks (§3.1). Figures 2−6 (a mechanism scheme and four EI mass spectra) carry no numbers I
used.

---

## 0. IDENTITY VERDICT — ★ CORRECT FILE, BUT THE "QUANTIFIED YIELDS" ARE IN A BAR CHART, NOT A TABLE

| Brief expected | Found |
|---|---|
| `10.1021/jf8012025`, Wang & Ho, DMHF via methylglyoxal | **EXACT MATCH.** Title, authors, journal, volume, issue, pages and article ID all confirmed from the paper's own header/footer. |
| "the sign-reversing pH ladder (MG+Cys rises with pH, MG+Gly falls)" | **★ CONFIRMED AND QUANTIFIED** — §3. And it is a **three-way** ladder, not two: MG alone, MG+Cys **and** MG+Gly, with MG alone the steepest riser. |
| "quantified yields" | **⚠️ QUALIFIED YES.** Real HPLC-UV quantification against an external standard, with SD and n = 3 — but published **only as a bar chart with no data table and no numeric annotations**. Every number in §3 of this dossier is a **digitisation [D-fig]**, not a transcription. |
| "the ¹³C₆ intact-skeleton test" | **★ CONFIRMED**, and there are **two** CAMOLA experiments, not one (§4). The second one — `[¹³C₆]glucose + [¹²C₃]MG` — is the more informative and is **not** mentioned in the abstract. |
| "a SECOND DMHF formation edge (from MGO)" | **★ CONFIRMED, and it is mechanistically DISTINCT**: the MG edge runs **without acetylformoin** (§4.3), which is a hard structural separation from the intact-skeleton edge. |

**No wrong-file problem.**

---

## 1. SYSTEM DEFINITION — verbatim (Materials and Methods, p. 7406)

### 1.1 Materials
> "2,5-Dimethyl-4-hydroxy-3(2H)-furanone, methylglyoxal (**40 wt % in water**), L-cysteine,
> L-glycine, sodium hydroxide and **[¹³C₆] or [¹²C₆] D-glucose** were purchased from Sigma Chemical
> Co. (St. Louis, MO). HPLC grade water, acetonitrile, dichloromethane and methanol were purchased
> from Fisher Scientific (Springfield, NJ)."

### 1.2 The reaction — verbatim
> "**Preparation of Model Systems.** MG, glucose, cysteine and glycine were dissolved in the
> **phosphate buffer (0.5 M, pH 3.0, 5.0, and 8.0)**, separately. **The pH was adjusted with 1 N
> sodium hydroxide.** The concentrations were **1.4 M, 1.4 and 1 M for MG, glucose and amino acids**,
> respectively. Two model systems were set up. The first one contained **an aliquot (1 mL) of MG and
> phosphate buffer solution**. In the second group, an aliquot of MG was mixed with either glycine
> or cysteine. All these samples were prepared in **sealed glass tubes** and heated at **120 °C for
> 60 min**. All reacted samples were cooled by an ice bath and **centrifuged at 14 × 1000 rpm
> (16000g) for 5 min** before HPLC analysis."

**Derived, arithmetic only:** MG : amino acid = **1.4 : 1.0 mol/mol**. Phosphate **0.5 M** —
i.e. **2.5× the phosphate of Blank 1997** (0.2 M), which matters given that Blank 1997 Table 3
measures a 3−9× specific phosphate catalysis on the *other* DMHF edge.

### 1.3 ⚠️ pH: SET-POINT ONLY, AND THE ADJUSTMENT IS ONE-SIDED
- The pH values 3.0 / 5.0 / 8.0 are **buffer set-points adjusted before heating with 1 N NaOH**.
- **No final pH is reported for any run.** No statement that pH was held during the 60 min.
- The adjustment reagent is **NaOH only** — no acid is mentioned — so the buffers were titrated
  **up** to their set-points from whatever the dissolved-precursor pH was. MG at 1.4 M is acidic
  and generates **pyruvic acid** during the run (§5), so a downward drift over 60 min at 120 °C is
  expected and unmeasured.
- **Carry as: "pH 3/5/8, 0.5 M phosphate, set-point, adjusted with NaOH before heating; hold not
  stated; no final pH reported."** Same defect class as `k3` §C.11 (Zhou 2023) and
  `blank1997_extraction.md` §1.1.

### 1.4 CAMOLA labelling experiments — verbatim
> "The carbon module labeling (CAMOLA) technique (19) was used to verify and evaluate the formation
> pathways of DMHF from MG and glucose. **A 1:1 mixture of 1.4 M [¹³C₆]glucose and 1.4 M
> [¹²C₆]glucose (containing natural abundant ¹³C) were reacted with 1 M glycine or cysteine at
> 120 °C and pH 5 for 60 min**, and another similar model, **a 1:1 mixture of 1.4 M [¹³C₆]glucose
> and 1.4 M [¹²C₃]methylglyoxal (containing natural abundant ¹³C) were reacted with 1 M glycine or
> cysteine at 120 °C and pH 5 for 60 min**. Then, the reaction mixture was extracted three times
> with 10 mL methylene chloride. The organic phase was separated, dried over anhydrous sodium
> sulfate, and concentrated under nitrogen gas for the GC-MS analysis."

*(⚠️ Note the CAMOLA reference number: the Methods text attributes CAMOLA to "(19)", which in the
reference list is **Namiki & Hayashi 1982, "Formation of novel free radicals products at the early
stage of Maillard reaction"** — a free-radical paper that has nothing to do with CAMOLA. The
Introduction attributes CAMOLA correctly to **"(18) Schieberle, P. The carbon module labeling
(CAMOLA) technique … Ann. N.Y. Acad. Sci. 2005, 1043, 236−248."* **The Methods citation is a
mis-numbered reference; the correct source is ref 18, Schieberle 2005.** Recorded so nobody chases
Namiki 1982 for a CAMOLA method.)*

### 1.5 Analytics — verbatim

**HPLC (this is what produced Figure 1):**
> "Dionex UltiMate 3000 LC Modules … **UV−vis detector (model VWD-3400)** … A **Luna C18
> (Phenomenex) column (150 × 4.6 mm i.d., 3 µm particle size)** was used for DMHF analysis. The
> column temperature was maintained at **25 °C** … The mobile phase … consisted of HPLC grade water
> with **0.15 % acetic acid (v/v; solvent A)**, acetonitrile (solvent B) and methanol (solvent C)
> with a constant flow rate set at **0.5 mL/min**. HPLC gradient programs were modified according to
> the method reported by Haleva-Toledo et al. (20) for DMHF analysis as follows: **6 % solvent B and
> 6 % solvent C, and they increased together to 15 % over 16 min, then decreased to 6 % over 4 min.
> The whole program ran for 20 min. DMHF was detected with a UV wavelength at 290 nm** and the
> injection volumes were 30 µL. **The external standard quantification method was applied in this
> study. Every single peak area for the quantification was laid in the linear range of standard
> curve.**"

**⚠️ EXTERNAL STANDARD, NOT INTERNAL.** No internal standard, no isotope dilution, no stated
recovery. The samples went straight from the sealed tube through a centrifugation into the HPLC
with no extraction step, which limits the exposure — but **any matrix suppression at 290 nm is
uncorrected, and no recovery figure is given.** This is a materially weaker quantification than
Blank 1997's SIDA, and the two must not be treated as equivalent-quality anchors.

**GC-MS (isotopomers only, not quantification):**
> "HP6890 … Agilent 6890 Series … mass spectral detector (**EI, 70 eV**). The column was **HP-1701
> (14 % (cyanopropyl-phenyl) methylpolysiloxane capillary, 60 m × 0.25 mm id, film thickness =
> 0.25 µm)**. The hydrogen, air, and makeup gas (helium) flow rate was at 30.0, 300.0, and
> 5.0 mL/min … **The injector was in 1:1 split mode.** The 1.0 mL/min constant carrier gas (helium)
> flow rate was set. The GC oven temperature was programmed as follows: the initial oven temperature
> **40 °C** was set and increased to **280 °C at a rate of 5 °C/min held for 12 min. The total run
> time was 60 min.** The injector temperature was **250 °C** and detector temperature was **250 °C**."

**Statistics — verbatim, and ⚠️ IT IS WRITTEN BACKWARDS:**
> "**Statistical Analysis.** Data were expressed as means ± standard deviation (SD) and represent
> **three independent analyses**. Statistical significance was examined using **Student's t test**
> comparison between the means. **A p value of >0.05 was considered statistically significant.**"

**⚠️ `p > 0.05` = significant is inverted; the convention is `p < 0.05`.** The identical sentence
also appears in the Figure 1 caption. **Every "significantly" in this paper's text is therefore
unanchored** — the reader cannot tell which comparisons the authors considered significant. **Use
the SD bars (which are real and readable, §3) and ignore the significance claims.** Recorded as a
paper defect, not corrected.

---

## 2. FIGURE 1 CAPTION — verbatim (p. 7406)

> **Figure 1.** Effect of pH on 2,5-dimethyl-4-hydroxy-3(2H)-furanone formation in **1.4 M MG
> phosphate buffer solution** incubated with or without **1.0 M glycine** or **1.0 M cysteine**:
> **MG only (open bars); MG−cysteine (slashed bars); MG−glycine (dotted bars).** Values for DMHF are
> the **means ± standard deviation (SD), each analyzed three independent times.** Statistical
> significance was examined using Student's *t* test comparison between the means. A *p* value of
> >0.05 was considered statistically significant.

**Y-axis label, read from the render: `DMHF content (mg/ mole MG)`. X-axis: `pH`, categories 3, 5, 8.
Y range 0−250 with major ticks every 50.**

---

## 3. ★★ FIGURE 1, DIGITISED — the paper's entire quantitative content

### 3.1 Digitisation method and its accuracy — stated so the numbers can be trusted or rejected

Page 2 rendered at 300 dpi, cropped to the plot, upsampled 2×. The y-axis was located as the longest
contiguous dark column (x = 1379 px, spanning y = 177 → 1154) and the baseline as the longest
contiguous dark row (y = 1136, x = 1363 → 3764). **Tick marks were then located independently** and
their positions back-solved through the linear calibration `value = (1136 − y) × 250/(1136 − 177)`:

| tick found at y (px) | value implied by the calibration | nominal |
|---|---|---|
| 177 | **250.0** | 250 |
| 369 | **199.9** | 200 |
| 561 | **149.9** | 150 |
| 749 | **100.9** | 100 |
| 943 | **50.3** | 50 |
| 1135 | **0.3** | 0 |

**Maximum deviation 0.9 units on a 250-unit axis — the calibration is linear to better than 0.4 %.**
Bar tops and error-bar caps were then read as the topmost dark pixel per column, grouped into
contiguous runs (bar tops span ≈ 175 px, error-bar caps ≈ 37 px, so the two are unambiguously
separable).

**Quoted precision: ±2 units on the bar heights, ±3 on the SDs.** This is a **[D-fig]**
digitisation, **not a transcription**, and every consumer of these numbers must carry that.

### 3.2 ★★ THE TABLE THE PAPER SHOULD HAVE PRINTED

**Units: mg DMHF per mole MG.** n = 3 independent analyses. Conditions: **1.4 M MG, 0.5 M phosphate,
120 °C, 60 min, sealed tube; amino acid, where present, 1.0 M.**

| pH | **MG only** (open) | **MG + cysteine** (slashed) | **MG + glycine** (dotted) |
|---|---|---|---|
| **3** | **4 ± 2** | **20 ± 2** | **72 ± 21** |
| **5** | **13 ± 1** | **34** (SD < 1, not resolvable) | **18 ± 3** |
| **8** | **199 ± 19** | **70 ± 16** | **8** (SD < 1, not resolvable) |

Raw digitised values before rounding, for the record:
MG only 4.1 / 13.3 / 199.3 (caps 6.3 / 14.3 / 218.7);
MG+Cys 20.5 / 34.2 / 69.6 (caps 22.7 / — / 85.5);
MG+Gly 71.5 / 18.2 / 7.8 (caps 92.8 / 21.4 / —).

### 3.3 ★★ THE THREE pH TRENDS — this is the paper's headline and it is a three-way pattern

| system | pH 3 → pH 8 | direction | fold |
|---|---|---|---|
| **MG alone** | 4 → 199 | **RISES, steeply** | **≈ 49×** |
| **MG + cysteine** | 20 → 70 | **RISES** | **≈ 3.5×** |
| **MG + glycine** | 72 → 8 | **★ FALLS** | **≈ 9× down** |

**Abstract, verbatim:** *"MG alone or MG with cysteine could produce **increased level of DMHF with
pH increased**, whereas **MG with glycine had contrary trend**."* — **the digitisation reproduces the
abstract's three-way claim exactly, in both sign and rank order, which is the best available
validation that the digitisation is right.**

### 3.4 ★ THE AMINO-ACID CROSSOVER — a second, independent structural fact in the same figure

Reading Figure 1 **down the columns** instead of across gives an ordering that reverses with pH:

| pH | ordering | best system | worst |
|---|---|---|---|
| **3** | **Gly (72) > Cys (20) > none (4)** | glycine | MG alone |
| **5** | **Cys (34) > Gly (18) ≈ none (13)** | cysteine | MG alone |
| **8** | **none (199) > Cys (70) > Gly (8)** | **MG alone** | glycine |

**★ At pH 3, adding glycine multiplies DMHF by 18×; at pH 8, adding glycine DIVIDES it by 26×.**
That is a **460× swing in the effect of a single amino acid, driven by pH alone**. It is a ratio, so
it is immune to the external-standard scale error. **This is the strongest hold-out in the paper.**

### 3.5 [DERIVED] mol %, so the two DMHF edges can be compared on one scale

MW(DMHF, C₆H₈O₃) = 128.13. `mol % on MG = mg/mol × 0.1 / 128.13 = mg/mol × 7.805 × 10⁻⁴`.

| pH | MG only | MG + Cys | MG + Gly |
|---|---|---|---|
| 3 | 0.0032 % | 0.0160 % | **0.0558 %** |
| 5 | 0.0104 % | 0.0267 % | 0.0142 % |
| 8 | **0.1556 %** | 0.0543 % | 0.0061 % |

**⚠️ THE BASIS IS PER MOLE OF MG, AND DMHF NEEDS TWO C3 UNITS.** On a carbon-conserving basis
(2 MG → 1 DMHF) the maximum conversion of MG carbon is **2 × 0.156 = 0.31 mol %**.

**Cross-edge comparison, with its caveat stated first:**
Blank 1997's best pentose cell is **0.0141 mol % on sugar** (EHMF, xylose/Ala 1:4) and its best HDMF
cell is **0.0040 mol %** (arabinose/Gly). Wang & Ho's best is **0.156 mol % on MG** — **≈ 39× higher
than Blank's best HDMF cell.** **BUT the temperatures differ (120 °C vs 90 °C) and the bases differ
(a C3 fragment vs a C5 sugar), so this is an order-of-magnitude orientation and NOT a
commensurable comparison.** What it does establish is that **both edges live in the
10⁻³ − 10⁻¹ mol % band**, i.e. the MGO edge is not negligible relative to the pentose edge and
cannot be dropped on magnitude grounds alone.

---

## 4. ★★ THE TWO CAMOLA EXPERIMENTS — the structural core of the paper

### 4.1 Experiment A: `[¹³C₆]glucose + [¹²C₆]glucose (1:1)` + glycine **or** cysteine, 120 °C, pH 5, 60 min

**Design logic, verbatim (p. 7407):**
> "**If the glucose carbon skeleton keeps intact in DMHF formation, equal molar of [¹³C₆] labeled
> DMHF and [¹²C₆] unlabeled DMHF should be obtained. However, if the fragmentation of glucose occurs
> before DMHF formation, up to seven isotopomers with different numbers of labeled carbons may be
> formed.**"

**Result, verbatim:**
> "The results demonstrated that **five isotopomers [¹³C₁] to [¹³C₅] from glucose fragmentation, were
> not observed in the presence of glycine and cysteine, and a 1:1 mixture of [¹³C₆]DMHF and
> [¹²C₆]DMHF was obtained** (Figure 3), suggesting **no breakdown of glucose during DMHF formation.**
> Generally, fragmentation degree is related to the temperature. **When reaction temperature was
> increased to 165 °C, still no fragmentation of glucose occurred in the DMHF formation (data not
> shown).**"

**And the same test applied to the intermediate:**
> "Previous studies suggested that acetylformoin was an important precursor which could be reduced
> to DMHF (17). In the current experiments, **molecular ion of m/z 144 representing the
> [¹²C₆]acetylformoin was present in an equal intensity to m/z 150 ([¹³C₆]acetylformoin)** indicating
> that **acetylformoin as one of DMHF intermediates was also generated from intact glucose**
> (Figure 4). It is therefore concluded that **in the presence of glycine or cysteine, DMHF can only
> be formed through the intact glucose.**"

**★ Two hard numbers here, both isotopomer ratios and therefore response-factor-immune:**

| # | measurement | value | meaning |
|---|---|---|---|
| A1 | `[¹³C₆]DMHF : [¹²C₆]DMHF` from a 1:1 labelled/unlabelled glucose feed | **1 : 1**, with **[¹³C₁]−[¹³C₅] ABSENT** | **Hexose → DMHF is 100 % intact-C6. Fragment recombination contributes ~0 %.** |
| A2 | `[¹²C₆]acetylformoin (m/z 144) : [¹³C₆]acetylformoin (m/z 150)` | **equal intensity, 1 : 1** | **Acetylformoin is itself made from the intact hexose, so the intact-C6 constraint holds at the INTERMEDIATE, not only at the product.** |

**★ The 165 °C null result is a temperature-robustness statement**, and it is one of the very few
temperature statements anywhere in the DMHF literature — but it is **"data not shown"**, so it
carries no number and cannot be scored.

### 4.2 ★★ Experiment B: `[¹³C₆]glucose + [¹²C₃]MG (1:1, both 1.4 M)` + glycine or cysteine, 120 °C, pH 5, 60 min

**This experiment is NOT in the abstract and it is the one that licenses a second edge.**

**Result, verbatim (p. 7408):**
> "As an important intermediate during thermal degradation of glucose, MG itself could generate DMHF
> with or without amino acids. **If a 1:1 mixture of [¹³C₆]glucose and [¹²C₃]MG reacted with glycine
> or cysteine at 120 °C, a 4:1 mixture of [¹³C₆]DMHF and [¹²C₆]DMHF was obtained. Because some of
> the MG involved in the Strecker degradation, only 20 % of DMHF was formed from MG, and the rest
> 80 % was formed from glucose** (Figure 5). **However, no [¹²C₆]acetylformoin was observed
> suggesting that acetylformoin was not a precursor during the DMHF formation from MG** (Figure 6).
> In all of DMHF and acetylformoin isotopomers, **only [¹³C₆] labeled and [¹²C₆] unlabeled, were
> observed. Glucose kept carbon skeleton intact during DMHF formation even its fragment MG was
> present, which indicated that MG and cyclization of intact glucose pathways were parallel since
> the precursors of these two pathways were different.**"

**★ Three hard results:**

| # | measurement | value | meaning |
|---|---|---|---|
| B1 | `[¹³C₆]DMHF : [¹²C₆]DMHF` with an equimolar `[¹³C₆]glucose / [¹²C₃]MG` feed | **4 : 1** | **★ A SPIKED MG POOL AT 1.4 M SUPPLIES 20 % OF THE DMHF; the hexose supplies 80 %.** The `[¹²C₆]`DMHF can only have come from **two ¹²C₃ MG units recombining** — so **the MGO edge IS a C3+C3 recombination and it demonstrably runs.** |
| B2 | acetylformoin isotopomers in the same run | **`[¹²C₆]`acetylformoin ABSENT**; only `[¹³C₆]` seen | **★★ THE MG EDGE BYPASSES ACETYLFORMOIN ENTIRELY.** This is the structural separation of the two edges — they do not share an intermediate. |
| B3 | absence of `[¹³C₃]`/mixed isotopomers in either experiment | **only [¹³C₆] and [¹²C₆] ever observed** | The two edges are **PARALLEL and non-mixing**: no DMHF is assembled from one glucose-derived C3 plus one MG-derived C3. |

### 4.3 ★ RECONCILING A1 AND B1 — they look contradictory and they are not

A1 says a glucose-only system gives **zero** fragment-recombination DMHF. B1 says a
glucose + spiked-MG system gives **20 %** from MG, i.e. from C3+C3 recombination. **Both are true
and the resolution is the MG concentration:**

- In experiment A there is no added MG. Any MG generated *in situ* is `¹³C₃` or `¹²C₃` in proportion
  to its parent glucose, and recombining one of each would give a **`[¹³C₃]`DMHF (m/z 131)** — one of
  the very isotopomers the authors looked for and **did not find**. So at the *in-situ* MG level in a
  1.4 M glucose system at 120 °C, **the C3+C3 edge carries a share small enough to be undetectable.**
- In experiment B the MG is **spiked to 1.4 M**, i.e. to parity with the glucose. At that
  concentration the same edge carries **20 %**.

**→ THE CORRECT READING: the MGO → DMHF edge is REAL and its rate is CONCENTRATION-CONTROLLED, and
at the MGO levels a hexose system generates on its own it is a MINOR-to-NEGLIGIBLE contributor.**
A model that gives the MGO edge a large share in a plain sugar+amine system contradicts A1; a model
that omits the edge entirely contradicts B1. **Both are falsifiable and they bracket the edge.**

**⚠️ Note also that experiment B's 20 % is NOT a clean branch fraction**: the authors themselves
attribute the shortfall to MG being consumed by Strecker degradation (*"Because some of the MG
involved in the Strecker degradation"*), so **20 % is a lower bound on what the MG pool would give
if it were not competed for.** No measurement of that competition is reported.

---

## 5. THE MECHANISM AS THE AUTHORS STATE IT

### 5.1 Figure 2 legend, verbatim
> **Figure 2.** Reaction pathway leading methylglyoxal to 2,5-dimethyl-4-hydroxy-3(2H)-furanone
> **via Cannizzaro reaction**.

### 5.2 The route in the authors' words (p. 7406−7407)
> "Generally, **MG and 1-hydroxy-2-propanone, the two major degradation products in Maillard
> Reaction, could react together to form DMHF via 2,5-dioxo-3,4-dihydroxyhexane** (18). However,
> when MG was heated alone at 120 °C, the formation of DMHF was observed (Figure 1), and the [DMHF]
> level was significantly increased as pH of the reaction increased. **MG, one of the dicarbonyl
> compounds, may transform into 1-hydroxy-2-propanone and pyruvic acid through the Cannizzaro
> reaction (Figure 2), and subsequently lead to DMHF by reacting MG with 1-hydroxy-2-propanone (18).
> DMHF formation from MG was pH-dependent because the Cannizzaro reaction is a base preferential
> reaction.**"

*(⚠️ **Typo in the source**: the sentence reads *"the **MG level** was significantly increased as pH
of the reaction increased"*. Figure 1 plots DMHF, not MG, and the abstract says DMHF. **Read it as
DMHF.** Recorded, not corrected.)*

**As a node-and-edge list:**

| # | step | note |
|---|---|---|
| 1 | **2 MG → MG + 1-hydroxy-2-propanone (acetol) + pyruvic acid** | **Cannizzaro disproportionation of MG. BASE-CATALYSED** — this is the authors' explanation for the pH rise. |
| 2 | **MG (C3) + 1-hydroxy-2-propanone (C3) → 2,5-dioxo-3,4-dihydroxyhexane (C6)** | aldol-type C3+C3 condensation, cited to Schieberle 2005 (ref 18) |
| 3 | **2,5-dioxo-3,4-dihydroxyhexane → DMHF** | cyclisation/dehydration |
| — | **acetylformoin is NOT on this route** | §4.2 B2, measured |

**★ Step 1 is the pH lever.** The authors' mechanism makes the MGO edge base-favoured *because the
reductant (acetol) is produced by a base-catalysed disproportionation*. That is a testable
mechanistic claim, and the `MG alone` column (4 → 13 → 199, **49×**) is exactly the shape a
base-catalysed rate-limiting first step predicts.

### 5.3 ★ The cysteine story — three DIFFERENT roles at three pHs, in the authors' own words (p. 7407−7408)
> "**At pH 8, Strecker degradation was the major reaction which consumed most of MG in the presence
> of amino acids. However, cysteine could be degraded into hydrogen sulfide which can be used as a
> reducing agent to produce 1-hydroxy-2-propanone from MG.** **At pH 5, Cannizzaro reaction and
> Strecker degradation became weaker, and reduction activity of cysteine was the main effect on DMHF
> formation, consequently cysteine reacting with MG generated a high level of DMHF.** **Cysteine may
> change its role from a reductant to an inhibitor at pH 3.** The inhibitory effect of thiol group in
> cysteine on DMHF formation has been observed particularly at pH 3 (21−23). Haleva-Toledo et al.
> (20) showed that **cysteine and N-acetylcysteine inhibited DMHF formation at pH 3 by a nucleophilic
> attack of thiol group to the open carbonyl form of DMHF** (20). **It was very interesting to
> observe that the generation of DMHF from MG and glycine increased as pH decreased. Further studies
> are required to explain this phenomenon.**"

**Three claims, three status levels:**

| claim | status |
|---|---|
| H₂S from cysteine reduces MG to acetol → **cysteine is a REDUCTANT feeding the edge** | **[P] proposed**, consistent with the measured pH 5 maximum of the MG+Cys/MG-alone ratio (34/13 = 2.6×) |
| At pH 3 cysteine's thiol **attacks DMHF's open carbonyl form → cysteine is a SINK** | **[C] cited** to Haleva-Toledo 1999 and Friedman/Molnar-Perl 1990 ×3, **not measured here** |
| **The MG+glycine falling trend is UNEXPLAINED** — verbatim: *"Further studies are required to explain this phenomenon."* | **★ THE AUTHORS DECLINE TO EXPLAIN THEIR OWN HEADLINE RESULT.** Any model that reproduces it does so without mechanistic guidance from this paper. |

**★★ THE CYSTEINE DOUBLE ROLE IS THE MOST MODEL-RELEVANT SENTENCE IN THE PAPER.** Cysteine is
simultaneously (a) a **source** of the reductant that drives the MGO edge and (b) a **sink** that
consumes DMHF by thiol addition. Which one wins is **pH-determined**, with the crossover between
pH 3 and pH 5. **The repo currently has neither: DMHF is only ever a product, never a reactant
(`grep -n dmhf src/reaction_templates.py` → two hits, both `products=[dmhf, ...]`).** See
`shu1988_extraction.md` for the fed-DMHF measurement of the sink side.

### 5.4 Conclusion, verbatim (p. 7408)
> "In conclusion, the results of this study indicate that **MG, depending on the pH differently
> affected DMHF generation in the presence or absence of amino acids. DMHF level increased as pH
> increased when cysteine reacted with MG, whereas the trend was reversed in the presence of
> glycine. When glucose reacted with glycine or cysteine, glucose skeleton kept intact in the formed
> DMHF as well as its precursor acetylformoin. Acetylformoin was not formed in the reaction between
> MG and either glycine or cysteine.**"

---

## 6. THE NAMING TRAP — what "HMF"/"HDMF"/"DMHF" denotes IN THIS PAPER

**There is no abbreviation table.** Abbreviations are defined inline at first use. I read every
occurrence:

| token | denotes, in this paper | defined where |
|---|---|---|
| **DMHF** | **2,5-Dimethyl-4-hydroxy-3(2H)-furanone**, *"known as Furaneol"* — **the repo's `DMHF`, i.e. furaneol.** | Title, abstract line 1, and Introduction sentence 1 (p. 7405) |
| **MG** | **Methylglyoxal**, *"also known as 2-oxopropanal or pyruvaldehyde"* — the repo's `MGO`. | Abstract; Introduction p. 7405 |
| **CAMOLA** | carbon module labeling | Introduction |
| **HDMF** | **NEVER USED.** Zero occurrences. | — |
| **"HMF" standing alone** | **NEVER USED in the body.** | — |
| **⚠️ `5-(hydroxymethyl)furfural`** | **APPEARS ONCE, in reference 20's TITLE** (Haleva-Toledo et al. 1999: *"Effects of L-cysteine and N-acetyl-L-cysteine on 4-hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol), **5-(hydroxymethyl)furfural**, and 5-methylfurfural formation and browning …"*). | Reference list, p. 7409 |
| **norfuraneol / 4-hydroxy-5-methyl-3(2H)-furanone** | **NOT MENTIONED ANYWHERE.** | — |

**→ THE TRAP IN THIS PAPER IS THE OPPOSITE ONE.** `DMHF` unambiguously means furaneol throughout the
body, but **the reference list contains the OTHER "HMF" (5-hydroxymethylfurfural) inside a title.**
Any automated harvest of "HMF"-like strings from this PDF will pick up a 5-HMF reference from a
DMHF paper. Also note the ordering of the letters: this paper writes **DMHF**, Blank 1996/1997 write
**HDMF**, Poisson 2019 writes **HDMF** — three orderings, one molecule.

**Norfuraneol is entirely absent from this paper**, so it offers **no** evidence for or against the
`norfuraneol ≫ DMHF` ordering.

---

## 7. KINETICS — **NONE**

No rate constant, no activation energy, no time course, no half-life. **One temperature (120 °C, with
a single unquantified "data not shown" excursion to 165 °C) and one time (60 min) for every
measurement in the paper.** Three pH points and three amino-acid conditions at that one (T, t).
**Nothing here can be turned into a rate parameter by any transformation.**

---

## 8. CITED, NOT MEASURED

| quantity / claim | source, as this paper cites it | ref # |
|---|---|---|
| DMHF from hexose via **acetylformoin reduction**, by **disproportionation OR a Strecker reaction with amino acids** | **Hofmann & Schieberle 2001**, *"Acetylformoin — an important progenitor of 4-hydroxy-2,5-dimethyl-3(2H)-furanone and 2-acetyltetrahydropyridine during thermal food processing"*, Flavour 2000 Proc. 6th Wartburg Aroma Symp., pp. 311−322 | 17 |
| **hexose or pentose can be cleaved into MG + 1-hydroxy-2-propanone to generate DMHF**, "**a major pathway when glucose is reacted with proline in aqueous solution**" | **Schieberle 2005**, CAMOLA, *Ann. N.Y. Acad. Sci.* **1043**, 236−248 | 18 |
| **6-deoxysugars (rhamnose) are MORE effective than hexoses at forming DMHF, via 2,3-dioxo-4,5-dihydroxyhexane, "which is not easily formed from hexoses and cannot be generated from pentoses"** | **Hofmann & Schieberle 1997**, *JAFC* **45**, 898−906 | 13 |
| DMHF formation studied in 6-deoxysugar/hexose/pentose systems | **Schieberle 1992**, ACS Symp. Ser. 490, 164−175 | 14 |
| DMHF from pentoses by **Strecker-aldehyde chain elongation with acetylformoin as intermediate** | **Blank & Fay 1996**, JAFC 44, 531−536 (`10.1021/jf950439o`) | 15 |
| DMHF reacts with **cysteine or H₂S to give meat-like aroma at low pH** | **Shu & Ho 1988**, JAFC **36**, 801−803 (`data/articles/shu1988.pdf`); **Zheng, Brown, Leding, Mussinan & Ho 1997**, JAFC 45, 894−897 | 10, 11 |
| DMHF + **phenylalanine** gives alkylpyrazines | **Jutta & Werner 1990**, *Z. Lebensm. Unters. Forsch.* **190**, 14−16 | 12 |
| **cysteine and N-acetylcysteine INHIBIT DMHF formation at pH 3 by nucleophilic thiol attack on the open carbonyl form** | **Haleva-Toledo, Naim, Zehavi & Rouseff 1999**, JAFC **47**, 4140−4145 | 20 |
| sulfur amino acids inhibit browning | **Friedman & Molnar-Perl 1990** ×3, JAFC 38, 1642 / 1648 / 1652 | 21−23 |
| 3(2H)-furanones from 2,3-enolisation via 1-deoxyosones | **Hodge, Mills & Fisher 1972**, Cereal Sci. Today 17, 34−40 | 16 |

**★ Two of these cited papers are already on disk**: `shu1988.pdf` (ref 10, dossiered in this wave)
and the whole Hofmann/Schieberle cluster. **Ref 13, Hofmann & Schieberle 1997 JAFC 45:898−906, is
the source of the claim that DMHF CANNOT be generated from pentoses via the 2,3-dioxo-4,5-
dihydroxyhexane route** — which, read literally, sits in tension with Blank & Fay's pentose route
(they use a *different* route, chain elongation, so there is no actual contradiction, but the
sentence would read as one if quoted alone). Flagged.

---

## 9. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured · **[M-fig]** measured but published only as a bar chart · **[D-fig]** my
digitisation of that chart · **[C]** cited · **[D]** derived arithmetic · **[P]** proposed by authors.
Conditions for all **[M]**/**[D-fig]** rows: **120 °C, 60 min, sealed glass tube, 0.5 M phosphate,
1.4 M MG (and/or 1.4 M glucose), 1.0 M amino acid, n = 3.**

| # | quantity | value | class | verdict |
|---|---|---|---|---|
| 1 | **DMHF from MG alone, pH 3 / 5 / 8** | **4 ± 2 / 13 ± 1 / 199 ± 19 mg per mol MG** | **[D-fig]** | ★ USE, with the digitisation caveat. |
| 2 | **DMHF from MG + cysteine, pH 3 / 5 / 8** | **20 ± 2 / 34 / 70 ± 16 mg per mol MG** | **[D-fig]** | ★ USE. |
| 3 | **DMHF from MG + glycine, pH 3 / 5 / 8** | **72 ± 21 / 18 ± 3 / 8 mg per mol MG** | **[D-fig]** | ★ USE. |
| 4 | **★ pH trend, MG alone** | **RISES ≈ 49× from pH 3 to 8** | **[D-fig]** / **[M]** (abstract states the sign) | ★★ USE — the **sign** is transcribed, the magnitude digitised. |
| 5 | **★ pH trend, MG + cysteine** | **RISES ≈ 3.5×** | same | ★★ USE. |
| 6 | **★★ pH trend, MG + glycine** | **FALLS ≈ 9×** | same | ★★ USE. **The sign-reversing result.** |
| 7 | **★★ Amino-acid effect swing across pH** | **glycine: ×18 at pH 3, ÷26 at pH 8 → a 460× swing** | **[D]** from #1/#3 | ★★ USE. Pure ratio, scale-free. **Best hold-out here.** |
| 8 | **mol % on MG** | **0.0032 − 0.156 mol %**; carbon basis (2 MG/DMHF) up to **0.31 mol %** | **[D]** | Use for cross-edge orientation only; see §3.5's caveat. |
| 9 | **★★ Hexose → DMHF isotopomer distribution** | **[¹³C₆] : [¹²C₆] = 1 : 1**, with **[¹³C₁]−[¹³C₅] ABSENT** | **[M]** | ★★ USE. **Hexose DMHF is 100 % intact-C6.** Response-factor-immune. |
| 10 | **★★ Acetylformoin isotopomer distribution, glucose system** | **m/z 144 : 150 = 1 : 1 (equal intensity)** | **[M]** | ★★ USE. The intact-C6 constraint holds at the **intermediate**. |
| 11 | **★★ MG spike branch fraction** | **[¹³C₆]DMHF : [¹²C₆]DMHF = 4 : 1 → 20 % from MG, 80 % from glucose** (1.4 M glucose + 1.4 M MG) | **[M]** | ★★ USE — but as a **LOWER BOUND** on the MG edge (§4.2, Strecker competition). |
| 12 | **★★ No acetylformoin on the MG route** | **[¹²C₆]acetylformoin ABSENT in the MG-spiked run** | **[M]** | ★★ USE. **The structural separation of the two edges.** |
| 13 | **The two edges do not mix** | only [¹³C₆] and [¹²C₆] isotopomers ever observed; **no [¹³C₃]DMHF** | **[M]** | ★ USE. No cross-assembly from one hexose-C3 + one MG-C3. |
| 14 | **In-situ MG edge share in a plain hexose system** | **below detection** (no [¹³C₃]DMHF in Experiment A) | **[M]**, upper bound | ★★ USE as a **CEILING**. Pairs with #11 as a bracket. |
| 15 | No glucose fragmentation at **165 °C** either | qualitative | **[M]**, **"data not shown"** | Note only. Not scorable. |
| 16 | Uncertainty | **SD from n = 3 independent analyses**; SDs on the bars: 2, 2, 21 / 1, <1, 3 / 19, 16, <1 | **[M-fig]**/**[D-fig]** | ★ USE as σ. |
| 17 | Quantification method | **HPLC-UV 290 nm, EXTERNAL standard, no internal standard, no recovery stated** | **[M]** | **⚠️ Materially weaker than Blank 1997's SIDA. Weight accordingly.** |
| 18 | Phosphate | **0.5 M**, one level | **[M]** | **No phosphate series.** 2.5× Blank 1997's 0.2 M. |
| 19 | pH conditioning | **set-points 3/5/8, adjusted with 1 N NaOH before heating; hold not stated; no final pH** | **[M]** | **⚠️ §1.3. Carry the caveat.** |
| 20 | **MG route: 2 MG → acetol + pyruvic acid (Cannizzaro), then MG + acetol → 2,5-dioxo-3,4-dihydroxyhexane → DMHF** | route | **[P]** authors, **[C]** Schieberle 2005 | ★ USE as topology for the MGO edge. |
| 21 | **Cannizzaro is base-preferential → the MG edge is base-favoured** | mechanism | **[P]** | ★ USE as the *reason* for #4; consistent with the data. |
| 22 | **Cysteine → H₂S → reductant producing acetol from MG** | mechanism | **[P]** | USE as a hypothesis only. Not measured. |
| 23 | **Cysteine as an INHIBITOR at pH 3, by thiol attack on DMHF's open carbonyl** | mechanism | **[C]** Haleva-Toledo 1999 | **Secondary. Not measured here.** Pairs with `shu1988`. |
| 24 | **The MG+glycine falling trend is UNEXPLAINED by the authors** | — | **[M]** admission | ★ Record. No mechanism is offered. |
| 25 | Statistical criterion as printed | **"p > 0.05 was considered statistically significant"** | **[M]**, **inverted** | **⚠️ Ignore all significance claims; use the SDs.** |
| 26 | rate constant, Ea, time course, temperature series | **NONE** | — | **ABSENT. §7.** |
| 27 | norfuraneol | **NOT MENTIONED** | — | **ABSENT.** No bearing on the norfuraneol ≫ DMHF ordering. |
| 28 | DMHF consumption rate / sink kinetics | **NONE** | — | **ABSENT.** The sink is cited (#23), never measured. |

---

## 10. CROSS-CHECKS AGAINST THE REST OF THE CLUSTER

| constraint | this paper's evidence | verdict |
|---|---|---|
| **`blank1996` #9: norfuraneol ≫ DMHF** | norfuraneol not mentioned | **SILENT.** No evidence either way. |
| **`blank1996` #7: HEMF requires alanine** | HEMF/EHMF not measured | **SILENT.** |
| **`blank1996` #3/#5: hexose DMHF is intact-C6** | **items #9, #10, #13** | **★★ CORROBORATED, at both the product and the acetylformoin intermediate.** This is the CAMOLA result the repo already pins in `tests/unit/test_chemistry_soundness.py` and in `src/reaction_templates.py:440−442`; **it is now transcribed from the primary source rather than from an abstract.** |
| **`blank1996` §5.3 step 7: acetylformoin → DMHF by reduction** | **item #10** puts acetylformoin on the intact-hexose route with the right isotopomer signature; **item #12** removes it from the MG route | **★ CORROBORATED for the hexose edge, EXCLUDED for the MGO edge.** Acetylformoin is edge-specific, not universal. |
| **`blank1997` §3.1: "methylglyoxal, acetol and dihydroxyacetone do generate HDMF (data not shown)"** — Blank's lab, 1997 | items #1−#3, #11 measure exactly that | **★★ CORROBORATED. Two independent labs (Nestlé Lausanne and Rutgers), eleven years apart, agree the MGO → DMHF edge exists.** Blank never published his; Wang & Ho did. |
| **`blank1997` Table 4: xylose/glycine HDMF RISES weakly with pH (1.19× over pH 5→7)** | **item #6: MG/glycine DMHF FALLS 9× over pH 3→8** | **⚠️ OPPOSITE SIGNS IN THE SAME AMINO-ACID SYSTEM.** Not a contradiction — different precursors, different pH windows, and Wang & Ho attribute their fall to Strecker consumption of MG rather than to the furanone step. **But it means the two edges CANNOT share one pH term.** Full treatment in `k5b_dmhf_synthesis.md` §4. |
| **`poisson2019` Table 10: HDMF from sucrose in coffee is 7.0−10.6 % C3/C3 recombination at mid-roast** | item #14 (in-situ C3/C3 below detection at 120 °C aqueous) and item #11 (20 % with a 1.4 M MG spike) | **★ CONSISTENT AND COMPLEMENTARY.** Poisson's 7−10.6 % in a real, hotter, drier matrix sits **between** Wang & Ho's aqueous in-situ null and their spiked 20 %. **Three independent points on the same edge, all small.** |

---

## 11. VERDICT

**A real second formation edge, established structurally and quantified only in a bar chart.**

| requirement | met? | detail |
|---|---|---|
| **A second DMHF formation edge, from a species the model already tracks** | **★★ YES.** `MGO` exists as a first-class species in `src/kinetic_core/species.py:77` and is produced in all three lanes (`r_dpo_c2c3`, `r_glc_c2c3`, `r_ama_mgo`). | The edge is `2 × C3 → DMHF`, acetylformoin-free. |
| **A quantitative anchor** | **⚠️ YES, but [D-fig].** Nine cells with SDs, digitised from a chart, external-standard HPLC. | **Not equivalent in quality to Blank 1997's SIDA.** |
| **A pH dependence** | **★★ YES, and sign-reversing.** | Three systems × three pHs. |
| **A structural constraint** | **★★ YES, two of them** (intact-C6; acetylformoin absent from the MG route). | Isotopomer ratios, response-factor-immune. |
| **A branch fraction** | **⚠️ PARTIAL.** 20 % at a 1:1 glucose:MG spike; **below detection** at in-situ MG. | A bracket, not a fraction. |
| **Kinetics** | **❌ NO.** One T, one t. | §7. |
| **The DMHF sink** | **❌ NO.** Cited to Haleva-Toledo 1999 and asserted for pH 3; never measured. | See `shu1988_extraction.md`. |

**One-line honest answer:** *Wang & Ho establish a second, acetylformoin-free, C3+C3 DMHF formation
edge from methylglyoxal and give it a sign-reversing pH ladder and an in-situ ceiling — but they
publish its only numbers as an un-annotated bar chart quantified against an external standard, so
the edge is STRUCTURALLY well-founded and NUMERICALLY the weakest anchor in this cluster.*

---

## 12. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION.**

### 12.1 Proposed **HOLD-OUT** — this paper's numbers belong here, not in a fit set

| candidate | proposed role | why |
|---|---|---|
| **★★ Item #9/#10 — hexose DMHF and its acetylformoin are 100 % intact-C6 ([¹³C₁]−[¹³C₅] absent)** | **HOLD-OUT, structural.** The repo already pins this; **re-declare it against the PRIMARY source** rather than the abstract-derived version in `src/reaction_templates.py:440−442`. | Isotopomer ratio; scale-free; falsifiable in one run; already implemented, so it is a free regression. |
| **★★ Item #12 — no acetylformoin on the MG route** | **HOLD-OUT, structural.** If an MGO edge is ever added, it must **not** be routed through acetylformoin. | A negative structural constraint is cheap to test and expensive to discover late. |
| **★★ Item #7 — the 460× glycine-effect swing across pH** | **HOLD-OUT, ratio.** | Pure ratio; unaffected by the external-standard weakness; a two-sided directional test a model can fail in either direction. |
| **★ Item #14 + #11 — the MGO edge bracket** (below detection in-situ; 20 % at a 1.4 M spike) | **HOLD-OUT, bracket.** | Bounds the edge from both sides without needing an absolute yield. |
| **Items #1−#3 — the nine absolute bar values** | **HOLD-OUT ONLY, at loose tolerance (suggest ≥3× band).** | **They are a digitisation of a chart, quantified by external standard with no recovery.** Suitable for "is the model in the right decade and does it have the right shape", not for a tight contract. |

### 12.2 Proposed **FIT** — ⚠️ recommend **NONE**

| candidate | recommendation |
|---|---|
| Any barrier or rate constant | **PROHIBITED.** One temperature (§7). |
| The absolute levels of items #1−#3 | **DO NOT FIT.** A fit target that is (a) digitised from a chart by an agent, (b) quantified against an external standard with no recovery, and (c) drawn from a system whose pH hold is not stated, is **three independent provenance defects deep**. Fitting a constant to it would be a slower-motion repeat of `thiol_addition_pentodiulose` (`src/barrier_constants.py:307`), where a barrier was fitted to a number that turned out not to be a measurement. **The difference here — and it is a real one — is that these ARE measurements; the defect is transmission, not fabrication. That still does not make them fit-grade.** |
| The pH **shapes** (#4/#5/#6) | **Fit-eligible ONLY as directional/sign constraints**, never as levels, and only if the model has independent edges to carry them. |

### 12.3 Prerequisite

**If anyone wants the nine Figure 1 values as anything tighter than a 3× band, the numeric data must
be obtained from the authors or from a supporting-information file — there is none in this PDF.**
Otherwise carry them permanently as `absolute_concentration: false`-class, the same treatment
`k4c_wave_consolidated.md` §220 and `k3` §1243 give to Shu 1988 and Hemmler 2018.

---

## 13. DECLARED GAPS — verbatim, for `k3` §C

> **"This paper has no tables. Its entire quantitative content is Figure 1, a nine-bar chart with no
> numeric annotations and no text layer in the plot area.** All nine values and their SDs in this
> dossier are a pixel digitisation performed during this extraction, calibrated against the axis
> tick marks to better than 0.4 % and validated against the abstract's three-way directional claim.
> **They are `[D-fig]`, not transcriptions, and must never be presented as quoted values.**"

> **"Quantification was by HPLC-UV at 290 nm against an EXTERNAL standard.** No internal standard,
> no isotope dilution, no recovery correction and no recovery figure are reported. This is a
> materially weaker quantification than the stable-isotope-dilution assay of Blank et al. 1997, and
> the two sources' absolute numbers must not be given equal weight."

> **"The statistical criterion is printed inverted: 'A p value of >0.05 was considered statistically
> significant', in both the Methods and the Figure 1 caption.** Every use of the word
> 'significantly' in this paper is therefore uninterpretable. Use the reported SDs (n = 3) and
> ignore the significance claims."

> **"pH 3 / 5 / 8 are set-points adjusted with 1 N NaOH before heating.** No final pH is reported
> for any run and no statement is made about whether pH was held over the 60 min at 120 °C. The
> system generates pyruvic acid by the authors' own proposed Cannizzaro step, so downward drift is
> expected and unmeasured."

> **"One temperature and one time.** 120 °C / 60 min for every quantified point. The single
> excursion to 165 °C is reported as 'data not shown' with no number. **No rate constant, no
> activation energy and no time course exist in this paper.**"

> **"The 20 % MG contribution (4:1 isotopomer ratio) is a LOWER BOUND, not a branch fraction.** The
> authors attribute the shortfall to MG being consumed by Strecker degradation and do not measure
> that competition."

> **"The authors do not explain their own headline result.** On the MG + glycine falling pH trend,
> p. 7408 verbatim: 'It was very interesting to observe that the generation of DMHF from MG and
> glycine increased as pH decreased. Further studies are required to explain this phenomenon.'"

> **"The Methods section mis-cites CAMOLA to reference (19), Namiki & Hayashi 1982, a free-radical
> paper.** The correct source, cited correctly in the Introduction, is reference (18), Schieberle
> 2005, Ann. N.Y. Acad. Sci. 1043:236−248."

> **"The DMHF SINK is cited, never measured here.** The claim that cysteine inhibits DMHF at pH 3 by
> nucleophilic thiol attack on the open carbonyl form is attributed to Haleva-Toledo et al. 1999,
> JAFC 47:4140−4145, and to Friedman & Molnar-Perl 1990. No consumption rate, no fed-DMHF
> experiment and no DMHF mass balance appears in this paper."

---

## 14. WHAT TO FETCH NEXT — ranked

| # | paper | why | confidence |
|---|---|---|---|
| 1 | **Haleva-Toledo, E.; Naim, M.; Zehavi, U.; Rouseff, R. L. 1999**, *"Effects of L-cysteine and N-acetyl-L-cysteine on 4-hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol), 5-(hydroxymethyl)furfural, and 5-methylfurfural formation and browning in buffer solutions containing either rhamnose or glucose and arginine"*, **JAFC 47, 4140−4145** | **★ THE DMHF SINK, WITH A pH AXIS AND A CYSTEINE AXIS, IN A BUFFER.** It is the source of this paper's only sink claim, it also carries **5-HMF and 5-methylfurfural in the same buffer** — so it would serve the HMF channel (`research_round3_channels.md` §A) at the same time. **Highest-value single fetch identified by this wave.** | Full citation transcribed verbatim from the reference list. Very high. |
| 2 | **Hofmann, T.; Schieberle, P. 2001**, *"Acetylformoin — an important progenitor of 4-hydroxy-2,5-dimethyl-3(2H)-furanone and 2-acetyltetrahydropyridine during thermal food processing"*, Flavour 2000 Proceedings, 6th Wartburg Aroma Symposium, pp. 311−322 | The primary source for **acetylformoin → DMHF by disproportionation or Strecker reaction with amino acids** — i.e. for the reductant identity that `src/barrier_constants.py:325` currently mis-attributes to Blank & Fay 1996 (defect **D1** of `blank1996_extraction.md` §0c). **Fetching this could RESOLVE D1.** | Conference proceedings, no DOI. Identity high; availability doubtful. |
| 3 | **Schieberle, P. 2005**, *"The carbon module labeling (CAMOLA) technique"*, **Ann. N.Y. Acad. Sci. 1043, 236−248** | Carries the claim that **hexose/pentose cleavage into MG + 1-hydroxy-2-propanone is "a major pathway when glucose is reacted with proline in aqueous solution"** — a *quantified-sounding* claim about the MGO edge's share that this paper only paraphrases. Also cited by Poisson 2019 (ref 8). | High. |
| 4 | **Hofmann, T.; Schieberle, P. 1997**, JAFC **45**, 898−906 | Source of the rhamnose ≫ hexose ≫ pentose DMHF reactivity ordering and of the claim that the 2,3-dioxo-4,5-dihydroxyhexane route **cannot** operate from pentoses. Would give the 6-deoxysugar arm the repo has no data on. | High. |
| 5 | **Zheng, Brown, Leding, Mussinan & Ho 1997**, JAFC **45**, 894−897, *"Formation of sulphur-containing flavour compounds from reactions of furaneol and cysteine, glutathione, hydrogen sulphide, and alanine/hydrogen sulfide"* | A **second fed-DMHF + cysteine study**, same group as Shu & Ho 1988, nine years later and with more sulfur donors. Would pair directly with `shu1988_extraction.md` on the sink side. ⚠️ Note `data/lit/extraction_dossiers/zheng1994_extraction.md` exists — **check whether that is this same group's earlier paper before ordering.** | High. |

---

## 15. SOURCES USED IN THIS EXTRACTION BEYOND THE PDF

**None.** Everything above is from `data/articles/wang2008.pdf` (text layer + a 300 dpi render of
page 2 measured in-place), from the sibling dossiers in this directory, and from read-only
inspection of `src/kinetic_core/species.py`, `src/kinetic_core/sulfur.py`,
`src/reaction_templates.py` and `src/barrier_constants.py` in the working tree. No web lookup was
performed. Nothing in `src/`, `tests/`, `results/` or `data/benchmarks/` was modified.

*End of dossier. Nothing outside this file was created or modified.*
