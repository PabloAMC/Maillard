# Shu & Ho 1988 (JAFC 36:801−803) — single-paper extraction, Wave K5b, 2026-08-29

> Shu, C.-K.*,†; Ho, C.-T.
> "**Effect of pH on the Volatile Formation from the Reaction between Cysteine and
> 2,5-Dimethyl-4-hydroxy-3(2H)-furanone**"
> *J. Agric. Food Chem.* **1988**, 36 (4), 801−803.
> Department of Food Science, Cook College, New Jersey Agricultural Experiment Station, Rutgers,
> The State University of New Jersey, New Brunswick, NJ 08903.
> † Present address: R. J. Reynolds Tobacco Co., Bowman Gray Technical Center, Winston-Salem,
> NC 27102.
> Received for review **July 18, 1986**. Revised **February 5, 1988**. Accepted **March 11, 1988**.
> CCC line: `0021-8561/88/1436-0801$01.50/0`. NJ Agric. Exp. Stn. Publication No. D-10205-7-86.
> **⚠️ NO DOI IS PRINTED ON THE PAPER** (1988 predates ACS DOI stamping). The DOI was **not**
> looked up and is **not** asserted here.

**Source PDF:** `data/articles/shu1988.pdf` (3 pp. of this article + the start of an unrelated
article, 407,790 bytes, on disk since 2026-08-28 11:45). Producer metadata: *Acrobat Capture 3.0 /
Adobe PDF Library 4.0, created 2001-08-29* — **it is an OCR'd scan of a 1988 print page.**
**Read method:** `pdftotext -layout` for prose. **⚠️ THE TEXT LAYER SCRAMBLES TABLE I's COLUMNS
BADLY** — it collapses the three pH columns wherever a row has empty cells, so that e.g.
`3,5-dimethyl-1,2,4-trithiolane` reads as `3.1 / 3.2` with the 24.7 orphaned on the previous line.
**Every cell of Table I below was therefore re-read visually from a 300 dpi render of page 2.**
Figure 1 (the yield-vs-pH plot) has **no numeric text layer** and was **digitised by pixel
measurement** (§4.1). **Do not trust a naive text-layer harvest of this PDF.**

**Status:** `research_round3_channels.md` §B.3 and §E identified this file as *"already in the repo,
no dossier … the DMHF analogue of the Whitfield fed-norfuraneol runs."* **This is that dossier.**

---

## 0. VERDICT UP FRONT

### 0a. IDENTITY — ★ CORRECT FILE

| Brief expected | Found |
|---|---|
| Shu & Ho 1988, JAFC 36:801−803 | **EXACT MATCH** — title, both authors, journal, volume, issue and page range all confirmed from the paper's own header and footer. |
| "fed-DMHF + cysteine at pH 2.2 / 5.1 / 7.1" | **★ CONFIRMED**, and the pH values are chosen for a stated reason: *"pH values of 2.2, 5.1, and 7.1 represented points **below, around, and above the isoelectric points of cysteine**."* |
| "the DMHF analogue of the Whitfield fed-norfuraneol runs" | **★ CONFIRMED and it is exactly that** — a stoichiometric 1:1 feed of the furanone and cysteine, heated, and the product spectrum read out. Whitfield & Mottram 1999 do the same for **norfuraneol** at pH 4.5. |
| "constrains DMHF CONSUMPTION by cysteine, the sink side of the channel" | **⚠️ PARTIALLY.** It constrains the **product spectrum** of that consumption completely, and the **total volatile yield** semi-quantitatively (Figure 1, digitised). **It does NOT measure how much DMHF was consumed** — there is no DMHF mass balance, no residual-DMHF measurement and no conversion figure anywhere. See §0b. |

**No wrong-file problem.**

### 0b. ★★ WHAT THIS PAPER IS AND IS NOT — read before using it

**IT IS:** the only fed-DMHF + cysteine experiment in the corpus; a three-point pH series; a complete
57-compound product spectrum with GC area percentages; **a total-volatile yield curve vs pH for
three systems (DMHF+Cys, DMHF alone, cysteine alone)**; and a set of clean **on/off** product
switches.

**IT IS NOT:** a quantification of DMHF consumption. There is **no residual DMHF measurement, no
conversion, no mass balance, no rate**, and **Table I is GC area % with no internal standard and no
response factors** — the same class `k3` §1243 already assigns it (*"GC area %, no internal
standard; pH-shape source only"*). **This dossier does not upgrade that class. It sharpens what the
shape is.**

### 0c. ★ THE ONE STRUCTURAL FACT THE MODEL COULD USE TODAY

**`2,5-dimethyl-4-hydroxy-3(2H)-THIOPHENONE` — the direct O→S ring swap of DMHF — is the marker of
DMHF + cysteine coupling, and it is 6.0 % / 5.8 % / NOT DETECTED at pH 2.2 / 5.1 / 7.1.**
**The DMHF→thiophenone sink switches OFF between pH 5.1 and pH 7.1.** Together with the identical
on/off behaviour of **2,4-hexanedione** (4.7 / 6.7 / not detected), that is a **two-marker,
same-direction pH gate on the DMHF sink**, and it is the most usable thing in the paper. See §5.2.

---

## 1. SYSTEM DEFINITION — verbatim (Experimental Section, p. 801)

> "Each solution was prepared by dissolving **0.05 mol of L-cysteine hydrochloride** (Ajinomoto Co.,
> Tokyo, Japan) and **0.05 mol of 2,5-dimethyl-4-hydroxy-3(2H)-furanone (99 %)** (International
> Flavors and Fragrances, Union Beach, NJ) in approximately 450 mL of distilled water, **titrating
> the solution to the appropriate pH value (2.2, 5.1, 7.1) with 10 % Na₂CO₃**, and diluting the
> sample to a **final volume of 500 mL** with distilled water. The solution was then placed in a
> **2-L Parr bomb** (Parr Instrument Co., Moline, IL) and **heated at 160 °C for 0.5 h**. The reaction
> mass was subjected to a **vacuum steam distillation, extraction, and concentration**, and the
> concentrate was analyzed by **gas chromatography−mass spectrometry (GC-MS) on fused silica columns
> (OV-1 and Carbowax 20M)** as described previously (Shu et al., 1985a). **Aroma properties of the
> volatile mixture from each reaction were evaluated by Dr. Manfred Vock, IFF Chief Flavorist.**"

**Derived, arithmetic only: 0.10 mol/L DMHF and 0.10 mol/L cysteine·HCl, exactly 1:1, in 500 mL.**

### 1.1 ⚠️ FOUR CONDITION CAVEATS, ALL LOAD-BEARING

1. **★★ THERE IS NO BUFFER. AT ALL.** The pH was set by **titrating with 10 % Na₂CO₃** and the
   solution was then sealed and heated. **No buffer salt appears anywhere in the paper.** In an
   unbuffered 0.1 M cysteine·HCl solution held at 160 °C for 30 min — a reaction that generates
   H₂S, CO₂, organic acids and NH₃ — **the pH cannot possibly have held.** **All three pH labels
   are INITIAL, ROOM-TEMPERATURE values on an unbuffered system.** This is the **most severe**
   instance of the pH-labelling defect in the whole DMHF cluster: Blank 1996/1997 and Wang & Ho at
   least had 0.2−0.5 M phosphate. **Carry every "pH 2.2 / 5.1 / 7.1" in this paper as
   `initial pH, unbuffered, final pH unknown`.**
2. **Sodium carbonate is a reactant, not an inert base.** At pH 7.1 the titration consumed
   substantially more Na₂CO₃ than at pH 2.2, so **the three arms differ in sodium and carbonate
   loading by an unreported amount**, and CO₂ is liberated into a sealed bomb. Same confounder class
   as `apriyantono1993_extraction.md` §1.3 item 3, but larger here because the pH span is wider.
3. **160 °C in a sealed Parr bomb** — well above the boiling point of water, so the system is at
   autogenous pressure. **This is the hottest system in the DMHF cluster** (compare 90 °C Blank,
   120 °C Wang & Ho, 171−220 °C Poisson's dry roast, ≈100 °C reflux Apriyantono).
4. **Isolation is vacuum steam distillation + extraction + concentration**, i.e. another thermally
   aggressive, volatility-selective workup with **no recovery correction and no internal standard**.

### 1.2 ⚠️ NO REPLICATION IS STATED ANYWHERE

I searched the full text for `n =`, `replicate`, `duplicate`, `triplicate`, `standard deviation`,
`SD`, `mean`, `±`. **There is nothing.** No table footnote, no figure error bar, no methods
sentence. **Every value in Table I and Figure 1 is, on the face of the paper, a single
determination of unknown precision.** *(Compare Blank 1997: "maximum SD ≤ 10 %, n ≥ 2 assays × 2
injections"; Wang & Ho 2008: SD from n = 3; Apriyantono 1993: ±16 %, means of two runs; Poisson
2019: also nothing.)*

---

## 2. ABSTRACT — verbatim (p. 801)

> "The effect of pH on the volatile components formed by the reaction between cysteine and
> 2,5-dimethyl-4-hydroxy-3(2H)-furanone was studied in a **Parr bomb model system**. pH values of
> **2.2, 5.1, and 7.1 represented points below, around, and above the isoelectric points of
> cysteine**. **At pH 2.2, the major components formed were 3-methyl-2-(2-oxopropyl)thiophene,
> 1,2,3-trithia-5-cycloheptene, 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone, and 2,4-hexanedione. At
> pH 5.1, the major components formed were 3,5-dimethyl-1,2,4-trithiolane isomers and
> 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone. At pH 7.1, various classes of heterocyclic compounds,
> including thiazoles, thiazolines, and pyrazines, were identified.** All three of the reactions
> produced roasted and meaty aromas, but **the reaction at pH 2.2 generated the best meat flavor** in
> terms of the roasted and meaty notes evaluated."

---

## 3. ★★ TABLE I — verbatim, re-read from the page image, all 57 compounds

> **Table I. Volatile Components Identified from the Reaction of Cysteine and DMHF at Different pH
> Values**
> Column head: **GC area, %** · pH 2.2 · pH 5.1 · pH 7.1
> Footnote *a*: **"T = trace, less than 1 %."**

*(Blank cell = compound not reported for that pH. The paper does not distinguish "not detected" from
"not reported"; there is no footnote for the blanks.)*

| component | **pH 2.2** | **pH 5.1** | **pH 7.1** |
|---|---|---|---|
| ***Acyclic Carbonyls*** | | | |
| acetaldehyde | **T**^a | **T** | **T** |
| ethyl acetate | | **4.3** | |
| acetol acetate | **T** | **T** | |
| acetoin acetate | **T** | **T** | |
| **acetone** | **★ 13.3** | **1.5** | **2.4** |
| methyl ethyl ketone | **T** | **4.8** | **5.4** |
| 2-pentanone | | **2.3** | **T** |
| 3-pentanone | | **2.4** | |
| **2,4-hexanedione** | **★ 4.7** | **★ 6.7** | **★ —** |
| 2,3-pentanedione | **T** | **T** | |
| acetol | | | **T** |
| 3-hydroxy-2-butanone | | | **T** |
| 3-hydroxy-2-pentanone | **1.1** | **3.6** | **1.4** |
| 2-hydroxy-3-pentanone | **T** | **1.2** | **1.0** |
| ***Cyclic Carbonyls*** | | | |
| 2,5-dimethyl-4,5-dihydro-3(2H)-furanone | | | **T** |
| **2,5-dimethyl-3(2H)-furanone** | **1.4** | **2.4** | **3.6** |
| 2,4,5-trimethyl-3(2H)-furanone | **T** | **T** | **1.3** |
| 2,4-dimethyl-5-ethyl-3(2H)-furanone | **T** | **T** | |
| 2-ethyl-2-cyclopenten-1-one | | | **T** |
| 5-methyl-2-ethyl-2-cyclopenten-1-one | | | **1.8** |
| 2-methyl-2-pentenolide | | | **T** |
| 2-ethyl-2-pentenolide | | | **1.2** |
| ***Sulfur Compounds*** | | | |
| **3,5-dimethyl-1,2,4-trithiolane** | **3.1** | **★ 24.7** | **3.2** |
| 4,6-dimethyl-1,2,3,5-tetrathiane | | **T** | |
| 3,6-dimethyl-1,2,4,5-tetrathiane | | **T** | |
| 3-methyl-1,2,4-trithiane | | | **T** |
| 2-thiophenethiol | **2.2** | | |
| 1-mercapto-2-propanone | | **T** | |
| 2-acetylthiophene | | | **1.6** |
| 3-methyl-2-acetylthiophene | | | **T** |
| 5-methyl-2-acetylthiophene | | | **1.4** |
| 3,5-dimethylthiophene-3-carboxaldehyde | | | **T** |
| **3-methyl-2-(2-oxopropyl)thiophene** | **★ 23.0** | **T** | |
| 2-methyl-3-propionylthiophene | **2.0** | | |
| tetrahydrothiophene | **T** | | |
| thienothiophene | **T** | **T** | |
| **★ 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** | **★ 6.0** | **★ 5.8** | **★ —** |
| 2,5-dimethyl-2,4-dihydroxy-3(2H)-thiophenone | **T** | **T** | |
| 2,5-dimethyl-2-hydroxy-3(2H)-thiophenone | **T** | **T** | |
| 3,4-dimethylisothiazole | | **T** | |
| 2,4-dimethyl-3-thiazoline | | | **1.1** |
| 2,4,5-trimethyl-3-thiazoline | | | **T** |
| 4,5-dimethyl-2-ethyl-3-thiazoline | | | **T** |
| thiazole | | **T** | |
| 2-methylthiazole | **T** | **T** | |
| **2-acetylthiazole** | **T** | **T** | **★ 7.8** |
| 2,5-dimethylthiazole | | **T** | |
| 2-methyl-5-ethylthiazole | | **T** | |
| 4,5-dimethyl-2-acetylthiazole | | | **T** |
| 2,4,5-trimethylthiazole | | **T** | **T** |
| 2-thiazolyl ethyl ketone | | | **4.6** |
| **1,2,3-trithia-5-cycloheptene** | **★ 10.6** | | |
| ***Nitrogen Compounds*** | | | |
| pyrazine | | | **T** |
| 2,5-dimethylpyrazine | | | **4.4** |
| 2,6-dimethylpyrazine | | | **1.9** |
| 2-methyl-5-ethylpyrazine | | | **T** |
| 2-methyl-6-ethylpyrazine | | | **T** |
| 2,3,5-trimethylpyrazine | | | **2.0** |
| 2,5-dimethyl-3-ethylpyrazine | | | **1.4** |
| 2,6-dimethyl-3-ethylpyrazine | | | **1.4** |
| 2,6-diethyl-3-methylpyrazine | | | **1.1** |
| 2,4,5-trimethyloxazole | | | **T** |
| 2,3-dimethylpiperidine | | | **T** |
| 2,4,6-trimethyl-1,3,5-dithiazine | | | **T** |

**⚠️ The columns do NOT sum to 100 %** (pH 2.2 ≈ 67 %, pH 5.1 ≈ 60 %, pH 7.1 ≈ 44 %, counting `T` as
0). **The paper reports GC area percentages of an unreported total, with an unreported number of
unidentified peaks. These are shares of the chromatogram, not of the volatile mass, and not of the
DMHF consumed.**

---

## 4. ★★ FIGURE 1, DIGITISED — the paper's only yield data, and it is a three-curve pH series

### 4.1 Digitisation method and accuracy

Page 2 rendered at 300 dpi, the figure region cropped and upsampled 2×. The y-axis was located as
the longest contiguous dark column (**x = 327 px, spanning y = 363 → 1254**) and the baseline as the
longest contiguous dark row (**y = 1252**). **Y-axis ticks were located independently** and
back-solved through the linear calibration `value = (1248 − y) × 100 / 203.5`:

| tick found at y (px) | value implied | nominal |
|---|---|---|
| 434 | **400.0** | 400 |
| 636 | **300.7** | 300 |
| 840 | **200.5** | 200 |
| 1044 | **100.2** | 100 |
| 1248 | **0.0** | 0 |

**Deviation ≤ 0.7 units on a 400-unit axis — linear to better than 0.2 %.** X-axis ticks were found
at **x = 503, 778, 1046** (pH 2, 4, 6), giving **135.7 px per pH unit**.

**⚠️ ONE HONEST CAVEAT ON THE X-AXIS.** The three curves' left endpoints all begin at
**x ≈ 548−556 px, i.e. pH ≈ 2.33−2.39**, not at pH 2.2, and their right endpoints at
**x ≈ 1188−1200, i.e. pH ≈ 7.05−7.14**. The plotted points are therefore **drawn a little off their
nominal pH values** — this is a hand-drafted 1988 line figure, not a plotted dataset.
**The y-values are reliable to ±5 mg; the x-positions are nominal.**

### 4.2 Figure 1 caption — verbatim (p. 802)
> **Figure 1.** pH effect on the **yield of volatiles** formed from the reaction of **cysteine and
> DMHF**, **DMHF degradation**, and **cysteine degradation**.

Y-axis label: **`MG`** (milligrams). X-axis: **`pH`**, 2 → 8. Curve labels, printed to the right of
each line: **`CYSTEINE & DMHF REACTION`** (top), **`DMHF DEGRADATION`** (middle),
**`CYSTEINE DEGRADATION`** (bottom).

**Source of the two comparison curves, verbatim (p. 801−802):**
> "**Figure 1 shows the yields obtained from the reaction between cysteine and DMHF as well as from
> the degradation of each reactant at different pH's (Shu et al., 1985a,b).**"

**★★ THE TWO DEGRADATION CURVES ARE CITED, NOT MEASURED HERE.** They come from **Shu, Hagedorn,
Mookherjee & Ho 1985a, JAFC 33, 442−446** (cysteine) and **Shu, Mookherjee & Ho 1985b, JAFC 33,
446−448** (DMHF). Only the top curve is this paper's own data.

### 4.3 ★★ THE TABLE THE PAPER SHOULD HAVE PRINTED

**Units: mg of total volatiles** (basis: 0.05 mol of each reactant in 500 mL; see §4.4).

| system | **pH 2.2** | **pH 5.1** | **pH 7.1** | shape | provenance |
|---|---|---|---|---|---|
| **★ Cysteine + DMHF reaction** | **≈ 178** | **≈ 321** | **≈ 247** | **★ PEAKED at pH 5.1** | **[M]** this paper |
| **DMHF degradation alone** | **≈ 160** | **≈ 122** | **≈ 93** | **monotone DECREASING, near-linear** | **[C]** Shu 1985b |
| **Cysteine degradation alone** | **≈ 57** | **≈ 202** | **≈ 25** | **★ SHARPLY PEAKED at pH 5.1** | **[C]** Shu 1985a |

Raw digitised values before rounding: reaction **177.9 / 320.6 / 246.9**; DMHF degradation
**160.4 / 121.9 / 93.1**; cysteine degradation **≈ 57 / 201.5 / ≈ 25**.
*(The DMHF-degradation curve is drawn as a straight line; its pH 5.1 value is an interpolation the
authors themselves drew, not a plotted point.)*

### 4.4 [DERIVED] the yield on a molar basis, so it can be compared

Basis: 0.05 mol DMHF (MW 128.13 → 6.41 g) + 0.05 mol cysteine·HCl (MW 157.62 → 7.88 g) = **14.29 g
of reactants**.

| system | max yield | **as mg total volatiles per mol DMHF fed** | **as % of reactant mass** |
|---|---|---|---|
| Cysteine + DMHF, pH 5.1 | 321 mg | **6420 mg/mol** | **2.25 %** |
| Cysteine + DMHF, pH 2.2 | 178 mg | 3560 mg/mol | 1.25 % |
| DMHF alone, pH 2.2 | 160 mg | 3200 mg/mol | 2.50 % *(of 6.41 g DMHF)* |

**★ Roughly 1−2.5 % of the fed mass ends up as recovered volatiles.** That is a **recovered-volatile
yield, not a conversion**: the other 97.5−99 % is unaccounted for (melanoidins, non-volatiles,
losses in the steam distillation) and the paper does not account for it.

### 4.5 ★★ THE THREE SHAPES, AND WHAT EACH ONE CONSTRAINS

| # | observation | constraint it places |
|---|---|---|
| **S1** | **The DMHF + cysteine reaction PEAKS at pH ≈ 5.1**, at 1.8× its pH 2.2 value and 1.3× its pH 7.1 value | The **coupling** of the two reactants is maximal near cysteine's isoelectric point. The authors' reading: *"[3,5-dimethyl-1,2,4-trithiolane] is apparently formed easily **when pH = pI and cysteine is in its neutral form**."* |
| **S2** | **DMHF ALONE degrades MONOTONICALLY MORE at LOW pH** — 160 → 122 → 93 mg from pH 2.2 → 7.1, a **1.7× fall** | **★ The DMHF SELF-degradation sink is ACID-favoured.** ⚠️ **[C]**, from Shu 1985b, not measured here. |
| **S3** | **Cysteine ALONE degrades in a SHARP peak at pH ≈ 5.1** — 57 → 202 → 25 mg, an **8× swing** | The H₂S/reductant supply peaks at the pI. **[C]**, from Shu 1985a. |
| **S4** | **The reaction curve's shape resembles S3 (peaked), not S2 (monotone)** | The authors' own conclusion, verbatim (p. 802): *"**From these data, the reaction between cysteine and DMHF appears to be more similar to cysteine degradation than to DMHF degradation. It is therefore likely that cysteine has a stronger influence on the reaction than DMHF.**"* **★ i.e. the RATE-CONTROLLING partner in the DMHF sink is the CYSTEINE, not the furanone.** |

**★ S2 is worth pausing on, because it appears to conflict with a citation in Blank 1997.**
Blank et al. 1997 (p. 2644) adjust to **pH 4.0** during workup *"which is the pH optimum for the
stability of HDMF in aqueous solutions (Hirvi et al., 1980)"* — i.e. DMHF is **most stable near
pH 4**. Shu's S2 curve says **more volatiles are produced from DMHF alone as pH falls toward 2.2**.
**These are compatible**: a stability optimum at pH 4 means degradation increases both below and
above it, and S2 measures only the *volatile* degradation products, which need not track total
decomposition (the alkaline branch may go preferentially to non-volatile browning). **But the two
should not be quoted side by side as if they agreed. Both are [C] from unread sources.**

---

## 5. ★★ THE PRODUCT-SPECTRUM RESULT — the sink's chemistry, resolved by pH

### 5.1 TABLE II — verbatim (p. 802). The authors' own summary.

> **Table II. Summary of pH Effect on the Volatile Formation from the Reaction of Cysteine and DMHF**

| pH | class of compounds formed |
|---|---|
| **2.2** | (1) **thiophenes and cyclic trisulfide dominating**; (2) **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone and 2,4-hexanedione next dominating**; (3) 3,5-dimethyl-1,2,4-trithiolanes, **low level**; (4) **no pyrazines formed** |
| **5.1** | (1) **3,5-dimethyl-1,2,4-trithiolanes dominating**; (2) **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone next dominating**; (3) **no pyrazines formed** |
| **7.1** | (1) **pyrazines, thiazoles, thiazolines formed**; (2) 3,5-dimethyl-1,2,4-trithiolanes, **low level**; (3) **no 2,4-hexanedione formed** |

### 5.2 ★★ THE FOUR CLEAN pH SWITCHES — the usable content

| # | marker | pH 2.2 | pH 5.1 | pH 7.1 | reading |
|---|---|---|---|---|---|
| **W1** | **2,5-dimethyl-4-hydroxy-3(2H)-THIOPHENONE** — the direct O→S ring swap of DMHF | **6.0** | **5.8** | **NOT DETECTED** | **★★ THE DMHF+CYSTEINE COUPLING MARKER. Flat across pH 2.2−5.1, then OFF at 7.1.** |
| **W2** | **2,4-hexanedione** | **4.7** | **6.7** | **NOT DETECTED** | **★ Same on/off gate as W1**, and Table II states it explicitly. Two independent markers, same switch. |
| **W3** | **3,5-dimethyl-1,2,4-trithiolane** | **3.1** | **★ 24.7** | **3.2** | **★★ PEAKED at pI, 8× over both neighbours.** The single largest cell in the table. |
| **W4** | **pyrazines** (9 species) | **NONE** | **NONE** | **★ ALL NINE PRESENT**, 2,5-dimethyl at 4.4 % | **★★ HARD ON/OFF at pH 7.1.** Authors: *"Pyrazines were identified only at this pH, suggesting that **the amino group of cysteine was more reactive at pH 7.1** than at pH 2.2 or 5.1."* |

**★★ W1 IS THE MOST DIRECTLY MODEL-RELEVANT NUMBER IN THE PAPER.** `2,5-dimethyl-4-hydroxy-3(2H)-
thiophenone` is DMHF with the ring oxygen replaced by sulfur — it can only come from a
thiol/H₂S attack on DMHF followed by ring closure. **It is therefore an unambiguous, structurally
specific tracer of DMHF consumption by cysteine**, and its behaviour is: **present at ≈ 6 % of the
chromatogram at pH 2.2 AND 5.1, absent at 7.1.** Two further members of the same family confirm the
same gate at trace level: `2,5-dimethyl-2,4-dihydroxy-3(2H)-thiophenone` (**T / T / —**) and
`2,5-dimethyl-2-hydroxy-3(2H)-thiophenone` (**T / T / —**).

### 5.3 The authors' compound-by-compound attributions — verbatim (p. 801)
> "At pH 2.2 the major components from this reaction were 3-methyl-2-(2-oxopropyl)thiophene,
> 1,2,3-trithia-5-cycloheptene, acetone, 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone, 2,4-hexanedione,
> and 3,5-dimethyl-1,2,4-trithiolane. **3-Methyl-2-(2-oxopropyl)thiophene,
> 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone, and 2,4-hexanedione appear to be formed by the
> interaction of cysteine and DMHF** (Shu et al., 1985c). **1,2,3-Trithia-5-cycloheptene is a
> cysteine degradation product** (Shu et al., 1985a), and **acetone can be derived from DMHF** (Shu
> et al., 1985b). **Cysteine degradation can also yield 3,5-dimethyl-1,2,4-trithiolane** (Shu et al.,
> 1985a), but this product may alternatively be formed by **the reaction of H₂S, another cysteine
> degradation product, with acetaldehyde, which can be derived from either cysteine or DMHF**."

**★ THE AUTHORS PARTITION THEIR OWN PRODUCT LIST BY ORIGIN. Transcribed as a table because it is a
route assignment, not an observation:**

| product | attributed origin | source of the attribution |
|---|---|---|
| **3-methyl-2-(2-oxopropyl)thiophene** (23.0 % at pH 2.2) | **★ cysteine × DMHF interaction** | **[C]** Shu 1985c |
| **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** | **★ cysteine × DMHF interaction** | **[C]** Shu 1985c |
| **2,4-hexanedione** | **★ cysteine × DMHF interaction** | **[C]** Shu 1985c |
| 1,2,3-trithia-5-cycloheptene (10.6 % at pH 2.2) | **cysteine degradation alone** | **[C]** Shu 1985a |
| acetone (13.3 % at pH 2.2) | **DMHF degradation alone** | **[C]** Shu 1985b |
| 3,5-dimethyl-1,2,4-trithiolane (24.7 % at pH 5.1) | **ambiguous** — cysteine degradation **OR** H₂S + acetaldehyde (acetaldehyde from either reactant) | **[C]** Shu 1985a + **[P]** authors |

**★ So of the six major products, THREE are declared to be genuine cross-products of the DMHF sink,
ONE is pure cysteine chemistry, ONE is pure DMHF chemistry, and ONE is ambiguous.** That is a
usable partition — and **every one of those attributions is CITED to the 1985 papers, none is
measured here.**

### 5.4 The pH 7.1 secondary-reaction argument — verbatim (p. 801)
> "At pH 7.1, several classes of compounds were identified, including **thiazoles, thiazolines,
> thiolanes, furanones, cyclopentenones, pyrazines, oxazoles, and dithiazine.** … **Several factors
> indicated that secondary reactions occurred readily at this pH. First, a large number of compounds
> were identified in the GC profile of the volatile fraction, but none could be designed as a major
> reaction product. Second, 2,4-hexanedione and the thiophenones were not detected as products, and
> 3,5-dimethyl-1,2,4-trithiolane accounted for only 3 % of the volatiles formed.**"

**★ The authors' own reading of the pH 7.1 arm is that the DMHF+cysteine PRIMARY products were made
and then CONSUMED**, not that they were never formed. **That matters for how W1/W2 should be
scored**: they are a gate on the *survival* of the primary products, and only conditionally a gate
on the *formation* of the coupling itself.

### 5.5 ★ TABLE III — verbatim (p. 802). Sensory, one expert, no panel.

> **Table III. Aroma Properties of the Volatiles Formed from the Reaction of Cysteine and DMHF at
> Different pH's**

| pH | aroma properties |
|---|---|
| **2.2** | **roasted, meaty, bread crust** |
| **5.1** | **oniony, roasted, biting** |
| **7.1** | **solvent, roasted, burnt** |

> "Comparison of the aroma properties of each volatile mixture indicates that **the best meat flavor
> was produced by the reaction at pH 2.2**, although all three of the reactions produced similar
> roasted/meaty characters."

**⚠️ THIS IS ONE PERSON'S DESCRIPTION.** Methods, verbatim: *"Aroma properties of the volatile
mixture from each reaction were evaluated by **Dr. Manfred Vock, IFF Chief Flavorist**."* No panel,
no replication, no scale, no blinding, no reference standards. **Descriptive only. Not a sensory
dataset.**

**★ Note the tension the authors leave unresolved: the BEST meat flavour is at pH 2.2, while the
HIGHEST volatile YIELD is at pH 5.1 (Figure 1).** Yield and quality anti-correlate here. That is a
genuinely useful qualitative fact for any flavour-optimisation objective and it is stated nowhere
else in the cluster.

---

## 6. THE NAMING TRAP — what "HMF"/"HDMF"/"DMHF" denotes IN THIS PAPER

**No abbreviation table.** One abbreviation is defined, in the first sentence of the body:

| token | occurrences | denotes, in this paper | defined where |
|---|---|---|---|
| **DMHF** | throughout | **2,5-Dimethyl-4-hydroxy-3(2H)-furanone** = **furaneol** = the repo's `DMHF`. | p. 801, body sentence 1: *"2,5-Dimethyl-4-hydroxy-3(2H)-furanone (DMHF) is an important α-dicarbonyl…"* |
| **`HDMF`** | **ZERO** | — | — |
| **`HMF` standing alone** | **ZERO** | — | — |
| 5-hydroxymethylfurfural | **ZERO** — not in the body, not in any reference title | — | — |
| **norfuraneol / 4-hydroxy-5-methyl-3(2H)-furanone** | **1**, in a reference title only | Tonsbeek, Plancken & van den Werkhof 1968, *"…Isolation of **4-Hydroxy-5-methyl-3(2H)-furanone** and Its 2,5-Homologue from Beef Broth"* | Reference list |

**→ TRAP-FREE. `DMHF` here means furaneol and nothing else.** Note the letter order: this paper and
Wang & Ho 2008 write **`DMHF`**; Blank 1996/1997 and Poisson 2019 write **`HDMF`**. Same molecule.

### 6.1 ⚠️ A CHEMICAL ERROR IN THE PAPER'S FIRST SENTENCE
> "**2,5-Dimethyl-4-hydroxy-3(2H)-furanone (DMHF) is an important α-dicarbonyl**, which has been
> identified in many food sources."

**DMHF is not an α-dicarbonyl.** It is an enolic 3(2H)-furanone: a ring ketone at C-3 with an enol
at C-4, i.e. an **α-hydroxy enone / reductone**, not a 1,2-diketone. *(Its tautomers and its
ring-opened form do present a 1,2-dicarbonyl-like reactivity, which is presumably what the authors
meant, but the statement as written is wrong.)* **Recorded because the sentence is quotable and
would propagate a species mis-classification if lifted.**

---

## 7. KINETICS — **NONE**

No rate constant, no activation energy, no time course, no half-life, no conversion.
**One temperature (160 °C) and one time (0.5 h) for every measurement in the paper.** Three pH
points at that single (T, t). **There is no second time point and no second temperature anywhere,
so not even a two-point rate can be built, and there is no measurement of how much DMHF was
consumed for such a rate to be expressed in.**

---

## 8. ★★ THE SINK QUESTION — what this paper does and does not settle

The brief asked for this paper because **"this constrains DMHF CONSUMPTION by cysteine, the sink
side of the channel."** Precisely what it delivers:

| sink question | answered? | evidence |
|---|---|---|
| **Does cysteine consume DMHF?** | **★★ YES, structurally and unambiguously.** | The thiophenone series (W1) is DMHF with S for O — it cannot arise any other way. Three members detected. |
| **What are the products?** | **★★ YES, completely — 57 compounds, three pH states.** | Table I, §3. |
| **Does the product spectrum depend on pH?** | **★★ YES, dramatically, with four clean switches.** | §5.2, Table II. |
| **Is the total volatile output pH-dependent?** | **★ YES — peaked at pH ≈ 5.1, 1.8× over pH 2.2.** | Figure 1, digitised, §4.3. |
| **Which partner controls the rate?** | **★ YES, qualitatively: CYSTEINE.** | §4.5 S4, the authors' own comparison of curve shapes. |
| **HOW MUCH DMHF is consumed?** | **❌ NO.** | No residual DMHF, no conversion, no mass balance, no yield of any single product on a molar basis. |
| **How FAST?** | **❌ NO.** | One T, one t. §7. |
| **At what pH is the SINK strongest?** | **⚠️ AMBIGUOUS, and the paper says so.** | Total volatiles peak at pH 5.1; the DMHF-specific coupling markers (W1/W2) are flat 2.2−5.1 and vanish at 7.1; but §5.4 argues the pH 7.1 disappearance is **secondary consumption of the primary products**, not absence of the coupling. **The paper cannot separate "the sink is off" from "the sink ran and its products were eaten."** |

**★ THE HONEST BOTTOM LINE FOR THE SINK: this paper establishes that the DMHF + cysteine sink EXISTS,
gives its complete product spectrum at three pH states, shows the total volatile output peaks near
cysteine's pI, and identifies cysteine as the controlling partner — and it supplies NO magnitude and
NO rate. A sink edge added on this paper alone would have a correct, evidenced structure and a
completely free rate constant.** That is the same position `blank1996_extraction.md` §9 described
for the *formation* edge before `blank1997` arrived — **and for the sink, there is no companion
paper. The nearest candidates are Haleva-Toledo 1999 and Hirvi 1980 (§12).**

---

## 9. CROSS-CHECKS AGAINST THE CLUSTER

| constraint | this paper's evidence | verdict |
|---|---|---|
| **`wang2008` §5.3: "Cysteine may change its role from a reductant to an inhibitor at pH 3"**, cited to Haleva-Toledo 1999 | **W1/W2**: the DMHF→thiophenone and 2,4-hexanedione coupling markers are **present at pH 2.2** at 6.0 % and 4.7 % | **★★ CORROBORATED, and this paper is the direct evidence.** Wang & Ho assert a *consumption* of DMHF by cysteine's thiol at low pH; Shu & Ho **show the products of exactly that consumption at pH 2.2**. **Two papers, twenty years apart, same Rutgers group — Wang & Ho even cite this paper as their ref 10.** |
| **`wang2008` Figure 1: MG + cysteine DMHF yield RISES with pH (20 → 34 → 70 mg/mol over pH 3 → 8)** | **W1/W2**: the DMHF **sink** markers are strongest at low pH and vanish at 7.1 | **★★ MUTUALLY CONSISTENT AND MECHANISTICALLY COHERENT.** If the cysteine **sink** is acid-favoured and switches off by pH 7−8, then **net** DMHF in a cysteine system must rise with pH — which is exactly Wang & Ho's measured curve. **Two independent experiments, one from the formation side and one from the sink side, produce the same pH-dependence with the same sign.** ★ This is the strongest formation/sink cross-validation available in the cluster. |
| **`apriyantono1993` §4.4 Introduction**, which cites *this* paper for "pyrazines only at pH 7.1" | **W4** — nine pyrazines at pH 7.1, none at 2.2 or 5.1 | **★ TRANSCRIPTION CONFIRMED.** Apriyantono's paraphrase of Shu & Ho is accurate. |
| **`apriyantono1993` items #10/#11: pyrazines 414× down on acidification (8940 → 21.6 nmol/mol)** | **W4** — pyrazines on/off between pH 5.1 and 7.1 | **★★ CORROBORATED ACROSS SYSTEMS.** Two completely different systems (xylose/lysine at ≈100 °C; DMHF/cysteine at 160 °C), same direction, same approximate threshold. **Pyrazine formation being strongly base-favoured is now a two-source result.** |
| **`k3` §B2.2, which already records this paper as a shape source**: *"three shapes in ONE experiment: monotone-acid thiophenes/thiols; peaked-at-pI trithiolanes; monotone-alkaline pyrazines"* | Table I, §5.2 | **★ CONFIRMED AND SHARPENED.** That register entry is correct. **This dossier adds the numbers (W1−W4) and adds Figure 1's three yield curves, which `k3` §B2.2 does not mention.** |
| **`blank1997` §11: "pH 4.0 is the HDMF stability optimum" ([C] Hirvi 1980)** | **S2**: DMHF-alone volatile output falls monotonically 160 → 93 mg from pH 2.2 → 7.1 | **⚠️ NOT DIRECTLY COMPARABLE.** §4.5. Volatile output ≠ total decomposition. **Both are [C] from unread sources; do not quote them as agreeing.** |
| **`poisson2019`** | cysteine-free system | **NO CONTACT.** |
| **The repo's `norfuraneol + cysteine/H₂S` chemistry (`furanone_reductive_sulfhydrylation`, `furanone_reductive_opening`)** | This is the **DMHF** analogue: **DMHF + cysteine → 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** (W1) is the exact structural counterpart of **norfuraneol + H₂S → 4-mercapto-5-methyl-3(2H)-furanone → MFT** in `src/barrier_constants.py:324` | **★★ A MISSING SIBLING EDGE, IDENTIFIED.** The repo implements the C5 furanone's thiol chemistry and has **no** C6 counterpart — DMHF is only ever a product (`grep dmhf src/reaction_templates.py` → two hits, both `products=[...]`). **This paper is the evidence that the C6 edge exists.** |

---

## 10. CONSOLIDATED PARAMETER TABLE

Legend: **[M]** measured here · **[M-fig]** measured but published only as a line figure ·
**[D-fig]** my digitisation · **[C]** cited · **[D]** derived · **[P]** proposed · **[S]** sensory.
Conditions for all **[M]** rows: **0.10 M DMHF + 0.10 M L-cysteine·HCl in 500 mL water, UNBUFFERED,
titrated to pH 2.2 / 5.1 / 7.1 with 10 % Na₂CO₃ before heating, 2-L Parr bomb, 160 °C, 0.5 h,
vacuum steam distillation + extraction + concentration, GC-MS on OV-1 and Carbowax 20M,
NO replication stated.**

| # | quantity | value | class | verdict |
|---|---|---|---|---|
| 1 | **★★ 2,5-dimethyl-4-hydroxy-3(2H)-THIOPHENONE (the DMHF O→S swap)** | **6.0 % / 5.8 % / NOT DETECTED** (GC area) at pH 2.2 / 5.1 / 7.1 | **[M]** | ★★ USE **as an on/off structural marker of the DMHF sink**. Not as a level. |
| 2 | **★ 2,4-hexanedione** | **4.7 / 6.7 / NOT DETECTED** | **[M]** | ★ USE as a **second, independent marker of the same gate**. |
| 3 | **The two minor thiophenones** (2,4-dihydroxy and 2-hydroxy) | **T / T / not detected** | **[M]** | ★ USE — **confirms the W1 gate at trace level, three members of one family.** |
| 4 | **★★ 3,5-dimethyl-1,2,4-trithiolane** | **3.1 / 24.7 / 3.2 %** — **PEAKED at pI, 8×** | **[M]** | ★★ USE as a shape. Largest cell in the table. |
| 5 | **★★ Pyrazines (9 species)** | **ALL ABSENT at pH 2.2 and 5.1; ALL PRESENT at pH 7.1** (2,5-dimethyl 4.4 %, 2,3,5-trimethyl 2.0 %, 2,6-dimethyl 1.9 %, three at 1.1−1.4 %, three trace) | **[M]** | ★★ USE as a **hard on/off**. Corroborated by `apriyantono1993`. |
| 6 | **3-methyl-2-(2-oxopropyl)thiophene** | **23.0 / T / not detected** | **[M]** | ★ USE as an acid-favoured on/off; **attributed to the cysteine × DMHF interaction** ([C] Shu 1985c). |
| 7 | **1,2,3-trithia-5-cycloheptene** | **10.6 / — / —** | **[M]** | USE; **pure cysteine-degradation product** ([C] Shu 1985a). |
| 8 | **acetone** | **13.3 / 1.5 / 2.4 %** | **[M]** | USE; **DMHF-derived** ([C] Shu 1985b). Acid-favoured, 9×. |
| 9 | **2-acetylthiazole** | **T / T / 7.8 %** | **[M]** | ★ USE — a base-favoured on-switch, ≥8×. |
| 10 | Thiazolines (3 species) | **absent / absent / 1.1 % + 2×T** | **[M]** | USE as an on/off at pH 7.1. |
| 11 | 2,5-dimethyl-3(2H)-furanone | **1.4 / 2.4 / 3.6 %** | **[M]** | USE — the only **monotone base-favoured** carbonyl. |
| 12 | 3-hydroxy-2-pentanone | 1.1 / 3.6 / 1.4 % | **[M]** | USE. Peaked. |
| 13 | Full 57-compound product spectrum, 3 pH states | Table I, §3 | **[M]** | ★ USE as a **structural** spectrum. **NOT as levels.** |
| 14 | **★★ Total volatiles, cysteine + DMHF** | **≈ 178 / 321 / 247 mg** at pH 2.2 / 5.1 / 7.1 — **PEAKED at 5.1** | **[D-fig]** | ★★ USE as a **shape**, with the digitisation caveat. |
| 15 | **Total volatiles, DMHF alone** | **≈ 160 / 122 / 93 mg** — **monotone decreasing, 1.7×** | **[C]** Shu 1985b, **[D-fig]** | ★ USE as a shape. **CITED, not measured here.** |
| 16 | **Total volatiles, cysteine alone** | **≈ 57 / 202 / 25 mg** — **sharply peaked at pI, 8×** | **[C]** Shu 1985a, **[D-fig]** | ★ USE as a shape. **CITED, not measured here.** |
| 17 | **★ The reaction's pH shape resembles CYSTEINE degradation, not DMHF degradation** | qualitative | **[M]** interpretation | ★★ USE. *"cysteine has a stronger influence on the reaction than DMHF."* |
| 18 | **[D] recovered-volatile yield** | **1.25−2.25 % of reactant mass**; **3560−6420 mg volatiles per mol DMHF fed** | **[D]** from #14 | Orientation only. **NOT a conversion.** |
| 19 | **DMHF consumed / conversion / mass balance** | **NONE** | — | **⚠️ ABSENT. §8.** |
| 20 | **Aroma at pH 2.2 / 5.1 / 7.1** | **"roasted, meaty, bread crust" / "oniony, roasted, biting" / "solvent, roasted, burnt"**; best meat flavour at **pH 2.2** | **[S]**, **one flavourist, no panel** | ★ Record as qualitative. **Not a sensory dataset.** |
| 21 | **★ Yield and quality ANTI-CORRELATE** | max yield at pH 5.1; best flavour at pH 2.2 | **[D]** from #14 + #20 | ★ Note — useful for any flavour objective, stated nowhere else in the cluster. |
| 22 | **Route attributions of the six major products** | 3 are cysteine × DMHF cross-products; 1 pure cysteine; 1 pure DMHF; 1 ambiguous | **[C]** Shu 1985a,b,c + **[P]** | ★ USE as a partition. **All attributions are cited, none measured here.** |
| 23 | **pH 7.1 = secondary reactions consumed the primary products** | qualitative, 3 stated indicators | **[P]** authors | ★ **CRITICAL for interpreting #1/#2.** §5.4. |
| 24 | Trace convention | **`T` = less than 1 %** | **[M]** | ★★ Defines the floor. |
| 25 | **Column sums** | **≈ 67 / 60 / 44 %** — **the table does NOT sum to 100** | **[D]** | ⚠️ GC area % of an unreported total with unreported unidentified peaks. |
| 26 | **Quantification basis** | **GC area %, NO internal standard, NO response factors** | **[M]** | **⚠️ `absolute_concentration: false`. Unchanged from `k3` §1243.** |
| 27 | **⚠️ Buffer** | **NONE. pH set by Na₂CO₃ titration before heating; final pH never measured.** | **[M]** | **⚠️⚠️ THE MOST SEVERE pH-LABEL DEFECT IN THE CLUSTER. §1.1.** |
| 28 | **Replication, SD, error bars** | **NONE STATED ANYWHERE** | — | **⚠️ ABSENT. §1.2. Do not assign a σ.** |
| 29 | Conditions | **0.1 M : 0.1 M, 500 mL, 160 °C, 0.5 h, Parr bomb** | **[M]** | Context. **The hottest system in the cluster.** |
| 30 | rate constant, Ea, time course, temperature series | **NONE** | — | **ABSENT. §7.** |
| 31 | norfuraneol | **not measured** (appears in one reference title only) | — | **ABSENT.** No bearing on the `norfuraneol ≫ DMHF` ordering. |
| 32 | **⚠️ "DMHF is an important α-dicarbonyl"** | — | **[M]**, **CHEMICALLY WRONG** | **Do not propagate. §6.1.** |

---

## 11. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR. NOT A DECLARATION.**

### 11.1 Proposed **HOLD-OUT**

| candidate | proposed role | why |
|---|---|---|
| **★★ Items #1 + #2 — the DMHF→thiophenone and 2,4-hexanedione coupling markers, present at pH 2.2 AND 5.1, absent at 7.1** | **HOLD-OUT, structural on/off, PAIRED.** Score the two markers together, and **only if a DMHF + cysteine sink edge is ever implemented.** | Two independent structural markers of the same gate. Free of any calibration. ⚠️ **Must carry §5.4's caveat**: the pH 7.1 absence is, on the authors' own reading, secondary consumption — so the hold-out tests *net survival*, not *coupling rate*. |
| **★★ Item #5 — pyrazines ALL absent at pH 2.2 and 5.1, ALL present at pH 7.1** | **HOLD-OUT, structural on/off. Score JOINTLY with `apriyantono1993` items #10/#11.** | Two systems, two labs, two sugars/precursors, same direction and threshold. **A joint hold-out across two papers is far harder to pass by accident than either alone.** |
| **★ Item #4 — trithiolane peaks 8× at cysteine's pI** | **HOLD-OUT, shape.** Only if the model represents cysteine's ionisation state. | `src/conditions.py` already applies a pH-ionisation correction by family-name substring; **this is a direct test of whether that correction has the right shape for a thiol/amine amphoteric species.** |
| **★ Item #14 — the reaction's total-volatile yield peaks at pH 5.1 (≈1.8× over pH 2.2)** | **HOLD-OUT, directional, LOOSE.** | **[D-fig]** digitisation of a hand-drafted 1988 figure; the x-positions are nominal (§4.1). Direction and rough magnitude only. |
| **Item #17 — the reaction tracks cysteine degradation, not DMHF degradation** | **HOLD-OUT, qualitative.** | A model in which the DMHF sink's pH shape is dominated by the furanone rather than the thiol contradicts it. |
| **Items #15/#16** | **NEITHER — they are [C] from Shu 1985a,b**, digitised at second hand. | Record as provenance for #17 only. |

### 11.2 Proposed **FIT** — **NONE. Recommend a hard prohibition, and say why in the register.**

| candidate | recommendation |
|---|---|
| **Anything in this paper** | **PROHIBITED AS A FIT TARGET.** Five independent reasons, any one sufficient: **(1)** **GC area %, no internal standard, no response factors** — the values are not concentrations; **(2)** **no replication or SD is stated anywhere** (§1.2), so no honest likelihood exists; **(3)** **the system is UNBUFFERED and the pH labels are initial-only** (§1.1) — the *independent variable* is not what it says it is; **(4)** **one temperature, one time** (§7) — nothing identifies a rate; **(5)** **no DMHF conversion or mass balance** (§8) — there is no sink magnitude to fit to. |
| **★ A specific warning for the register** | This paper is a **near-perfect setup for the failure mode `src/barrier_constants.py:307` records**: a sulfur-lane barrier with **zero absolute literature anchors**, and a fed-precursor sulfur paper that *looks* like it supplies one. **It does not.** `thiol_addition_pentodiulose` was fitted to `cys_ribose_140C_Hofmann1998`'s 342/200 ppb, which Wave S2b proved was a repo-internal derivation, and had to be reverted. **A `thiol_addition_dmhf` or similar fitted to Shu & Ho's 6.0 % GC area would repeat that exactly, with the added defect that 6.0 % is not even a concentration.** Record as a **prohibited derivation**, by name, before anyone builds the edge. |

### 11.3 ★ The edge this paper licenses — proposal only

**A DMHF sink edge is structurally warranted and numerically unconstrained:**

```
DMHF + H2S  ->  2,5-dimethyl-4-hydroxy-3(2H)-thiophenone + H2O
    C6H8O3 + H2S -> C6H8O2S + H2O                          [balances exactly]
```

**This is the exact C6 counterpart of the C5 edge the repo already carries as
`furanone_reductive_sulfhydrylation` (`src/barrier_constants.py:324`, norfuraneol + H₂S + 2[H] →
2-methyl-3-furanthiol + 2 H₂O).** The repo has the C5 thiol chemistry and no C6 counterpart at all.

**If it is added, on this evidence:**
- its **barrier must be an ESTIMATE, explicitly labelled UNCONSTRAINED**, in the same register as
  `thiol_addition_pentodiulose` and `deoxyosone_reduction`;
- the natural neutral starting value is the **`thiol_addition` class value 28.60**, i.e. the same
  un-fitted class value Wave N gave its siblings — **not** a value tuned to make item #1 come out at
  6 %;
- items #1/#2/#5 become **hold-outs on the edge**, never fit targets;
- **it must be built only if the orchestrator wants it.** Adding a sink with a free rate to a
  channel whose formation edges are only now becoming calibratable would put two unconstrained
  parameters in series. **`blank1997` calibrates formation; nothing calibrates this.** Sequencing
  matters more than completeness here.

---

## 12. DECLARED GAPS — verbatim, for `k3` §C

> **"THE SYSTEM IS UNBUFFERED.** pH 2.2, 5.1 and 7.1 were set by titrating an aqueous
> 0.1 M cysteine·HCl + 0.1 M DMHF solution with 10 % Na₂CO₃ **before** sealing it in a Parr bomb at
> 160 °C for 30 min. No buffer salt appears anywhere in the paper and no final pH is reported. A
> reaction generating H₂S, CO₂, NH₃ and organic acids in an unbuffered solution at 160 °C cannot
> have held its pH. **All three pH labels are initial, room-temperature values, and this is the most
> severe instance of the initial-pH defect in the DMHF cluster.** The three arms also differ in
> sodium and carbonate loading by an unreported amount, because holding a higher pH required more
> Na₂CO₃."

> **"Table I is GC area percentage with no internal standard and no response factors.** The columns
> sum to approximately 67 %, 60 % and 44 %, so an unreported number of unidentified peaks is
> excluded. These are shares of a chromatogram, not concentrations, not molar yields, and not shares
> of the DMHF consumed. `absolute_concentration: false` — unchanged from `k3` §1243."

> **"No replicate count, no standard deviation and no error bar appears anywhere in this paper.**
> Every value in Table I and Figure 1 must be treated as a single determination of unknown
> precision. No σ may be assigned."

> **"There is NO measurement of how much DMHF was consumed.** No residual DMHF, no conversion, no
> mass balance, no molar yield of any product. The paper establishes that the DMHF + cysteine sink
> EXISTS and gives its complete product spectrum and its pH shape; it supplies no magnitude and no
> rate. **A sink edge built on this paper alone would have a correct structure and a completely free
> rate constant.**"

> **"One temperature and one time.** 160 °C / 0.5 h for every measurement. No rate constant, no
> activation energy, no time course and no temperature series exist in this paper."

> **"Figure 1's two comparison curves are CITED, not measured here.** The DMHF-degradation and
> cysteine-degradation curves come from Shu, Mookherjee & Ho 1985b (JAFC 33, 446−448) and Shu,
> Hagedorn, Mookherjee & Ho 1985a (JAFC 33, 442−446) respectively. Only the top curve
> (cysteine + DMHF) is this paper's own data. All three were digitised at second hand from a
> hand-drafted line figure whose plotted points sit a little off their nominal pH values."

> **"The route attributions of the major products are all CITED, none measured here.** That
> 3-methyl-2-(2-oxopropyl)thiophene, 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone and 2,4-hexanedione
> arise from the cysteine × DMHF interaction is attributed to Shu et al. 1985c (JAFC 33, 638−641);
> that 1,2,3-trithia-5-cycloheptene is a cysteine-degradation product to Shu et al. 1985a; that
> acetone derives from DMHF to Shu et al. 1985b. **None of these origin assignments is demonstrated
> by an experiment in this paper.**"

> **"The pH 7.1 disappearance of the coupling markers is ambiguous by the authors' own reading.**
> p. 801: 'Several factors indicated that secondary reactions occurred readily at this pH … 2,4-
> hexanedione and the thiophenones were not detected as products.' The paper cannot distinguish
> 'the DMHF+cysteine coupling did not run at pH 7.1' from 'it ran and its products were consumed'."

> **"The aroma data are one flavourist's free description, not a sensory dataset.** Table III's nine
> descriptors come from a single named expert (Dr Manfred Vock, IFF) with no panel, no scale, no
> replication, no blinding and no reference standards."

> **"The paper's opening sentence mis-classifies DMHF as 'an important α-dicarbonyl'.** It is an
> enolic 3(2H)-furanone / reductone, not a 1,2-diketone. Do not propagate the classification."

> **"No DOI is printed on this 1988 paper and none was looked up.** Cite it as J. Agric. Food Chem.
> 1988, 36 (4), 801−803 until a DOI is Crossref-verified."

---

## 13. WHAT TO FETCH NEXT — ranked

| # | paper | why | confidence |
|---|---|---|---|
| 1 | **★ Shu, C.-K.; Hagedorn, M. L.; Mookherjee, B. D.; Ho, C.-T. 1985c**, *"Two Novel 2-Hydroxy-3(2H)-thiophenones from the Reaction between Cystine and 2,5-Dimethyl-4-hydroxy-3(2H)-furanone"*, **JAFC 33, 638−641** | **The source of every "formed by the interaction of cysteine and DMHF" attribution in this paper** (§5.3), and it is a **cystine** (disulfide) rather than cysteine feed — a second sulfur donor on the same sink. Without it the route assignments here are unsupported. | Full citation transcribed from the reference list. High. |
| 2 | **Shu, C.-K.; Mookherjee, B. D.; Ho, C.-T. 1985b**, *"Volatile Components of the Thermal Degradation of 2,5-Dimethyl-4-hydroxy-3(2H)-furanone"*, **JAFC 33, 446−448** | **★ THE DMHF SELF-DEGRADATION PAPER — a pure DMHF sink with a pH axis** (it is the source of Figure 1's middle curve). This is the closest thing in the literature to an uncomplicated DMHF stability/decomposition dataset. | High. |
| 3 | **Shu, C.-K.; Hagedorn, M. L.; Mookherjee, B. D.; Ho, C.-T. 1985a**, *"pH Effect on the Volatile Components in the Thermal Degradation of Cysteine"*, **JAFC 33, 442−446** | Figure 1's bottom curve, and the analytical method this paper says it follows *"as described previously."* **The repo's sulfur lane has no cysteine-alone pH-degradation source.** | High. |
| 4 | **★★ Haleva-Toledo, E.; Naim, M.; Zehavi, U.; Rouseff, R. L. 1999**, JAFC **47**, 4140−4145 | **Already ranked #1 in `wang2008_extraction.md` §14 — this dossier RAISES it further.** It is the only identified source that quantifies DMHF *inhibition* by cysteine against pH **in a buffer**, i.e. it could supply the sink magnitude that Shu & Ho structurally establishes and cannot measure. **Two of the five papers in this wave point at it independently.** | Very high. |
| 5 | **Hirvi, T.; Honkanen, E.; Pyysalo, T. 1980**, *"Stability of 2,5-dimethyl-4-hydroxy-3(2H)-furanone and 2,5-dimethyl-4-methoxy-3(2H)-furanone in aqueous buffer solutions"*, **LWT 13, 324−325** | Cited by `blank1997` for the pH-4 stability optimum. **A DMHF stability-vs-pH study IN BUFFER** — the buffered counterpart of Figure 1's middle curve, and the thing that would settle the §4.5 tension. Short paper. | High (transcribed from `blank1997`'s reference list). |
| 6 | **Shu, C.-K.; Hagedorn, M. L.; Ho, C.-T. 1986**, *"Two Novel Thiophenes Identified from the Reaction between Cysteine and 2,5-Dimethyl-4-hydroxy-3(2H)-furanone"*, **JAFC 34, 344−346** | The paper this one says it is a follow-up to (*"We have previously reported the compounds identified in this reaction"*). Identification support for Table I. | High; lower value. |

---

## 14. SOURCES USED IN THIS EXTRACTION BEYOND THE PDF

**None.** Every number and quotation is from `data/articles/shu1988.pdf` — the text layer for prose,
and a **300 dpi render of page 2 read visually for every cell of Table I** (the text layer scrambles
the columns) plus **pixel-measured for Figure 1** (no numeric text layer). Comparisons are to
sibling dossiers in this directory and to read-only inspection of `src/barrier_constants.py`,
`src/reaction_templates.py` and `src/conditions.py` in the working tree. **No DOI lookup and no web
search were performed.** Nothing in `src/`, `tests/`, `results/` or `data/benchmarks/` was modified.

*End of dossier. Nothing outside this file was created or modified.*
