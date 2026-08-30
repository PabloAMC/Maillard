# Göncüoğlu Taş & Gökmen 2017 — Maillard reaction and caramelization during hazelnut roasting

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### Wave brief: *"extract STRUCTURE, flag every T-dependence as non-transferable."* Done — and the paper turns out to contain **no activation energies at all** in the file we hold.

---

## §0. IDENTITY — ⚠️ THE FILENAME YEAR IS WRONG

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/Goncouglu2016.pdf` (729,635 bytes) | `ls` |
| **DOI** | **`10.1016/j.foodchem.2016.11.159`** ✔ matches the round-3 expectation | Elsevier cover page, line 9 |
| PII | `S0308-8146(16)32008-8` · Reference `FOCH 20292` | cover page |
| title | *"Maillard Reaction and Caramelization during Hazelnut Roasting: A multiresponse kinetic study"* | cover page |
| authors | **Neslihan Göncüoğlu Taş**, Vural Gökmen | cover page |
| **correct citation year** | ★ **2017**, *Food Chemistry* **221**, 1911–1922 | Accepted 30 Nov 2016, **published 2017**; volume/pages from `research_round3_channels.md` §A.2 row 3 (Crossref-verified there), **not** printed in this PDF |
| dates | Received 30 Aug 2016 · Revised 27 Nov 2016 · Accepted 30 Nov 2016 | cover page |
| version | **Elsevier "Accepted Manuscript"** — not the version of record | cover page |
| funding | TÜBİTAK IntenC German–Turkish programme, Project 113O178 | Acknowledgements |

**⚠️ NAMING WARNINGS, REPORT UPSTREAM.**
1. **The file is named `Goncouglu2016` but the paper is cited as 2017.** Anything keying on the
   filename year will mis-cite it. I have named this dossier `goncuoglutas2017_extraction.md` to
   match the citation, and flag the file stem here.
2. ⚠️ **There is a SECOND, UNRELATED file with a near-identical stem on disk:**
   `data/articles/Goncouglu2026.pdf` = *"Effect of extrusion parameters and feed composition on
   physical characteristics, aroma profile and acrylamide content in pea protein-enriched corn
   extrudates"*, **Food Chemistry: X 36 (2026) 104019**, also by Neslihan Göncüoğlu Taş. It is an
   **extrusion** paper and belongs to the matrix/extrusion lane, **not** to the HMF cluster.
   `Goncouglu2016` and `Goncouglu2026` differ by one digit and must not be confused.
3. The author's surname is **Göncüoğlu Taş** (two words); the file stem drops the `Taş`. The JAFC
   companion paper's reference list spells it `Göncüoğlu`; the Şen 2022 paper cites her as
   `Taş and Gökmen (2017)`. All three refer to this paper.

---

## §1. ★★ THE PRE-REGISTERED ARRHENIUS DISCLAIMER — and what is actually in the file

The abstract's last sentence, verbatim (line 26–27):

> "**The temperature dependence of the reactions was complicated and could not be explained by the
> Arrhenius equation.**"

The Conclusion repeats it, verbatim (line 616–617):

> "**The temperature dependence of the reactions was found to be more complicated than defined by
> the Arrhenius equation.**"

And §3.2 gives the numbers behind it — **this is the passage the wave brief was after**, verbatim
(lines 405–421):

> "The temperature dependence of elementary reactions during hazelnut roasting was determined by the
> activation energies (Ea) and reaction rate constants (kb) at reference temperature of 160 °C
> (**Table S4 in Supplementary Material**). **The activation energies of elementary reaction steps
> were found to range between 0‒1174 kJ/mol with six zero and a few relatively high values.**
> However, it was reported that the activation energies of most of the chemical reactions were in
> the order of 120 kJ/mol (Van Boekel, 1998) … **The inconsistent activation energy values of
> elementary steps of Maillard reaction and caramelization during hazelnut roasting might be
> explained with the effect of different complex mechanisms rather than the temperature dependence
> defined by the Arrhenius equation that made an over simplification** (Martins et al., 2001). The
> model fits at 150, 160 and 170 °C were obtained all at once during parameter estimation by using
> the reparametrized Arrhenius equation (Fig. S2) … **it should be expressed that these findings
> represented a temperature range of 150‒170 °C and it might be better to study a wider temperature
> range, by including lower temperatures (100‒120 °C), to explain the Arrhenius behavior of the
> reactions.**"

### ★ THREE THINGS THIS ESTABLISHES

**(1) `Ea = 1174 kJ/mol` is a chemically absurd number, and it is the authors' own upper bound.**
For scale: the C–C bond dissociation energy is ~350 kJ/mol; van Boekel's food-reaction norm is
~120 kJ/mol; the largest Ea anywhere in `k3` §A is 131.8. **A 1174 kJ/mol "barrier" is not a
barrier; it is a fitting artefact.** Any Ea set whose range includes it is disqualified as a whole,
not row by row.

**(2) SIX of the ~26 activation energies were FIXED TO ZERO** — i.e. the authors declared six
elementary steps to be temperature-independent over 150–170 °C. (The JAFC companion uses exactly
this device, with the footnote *"Zero activation energy indicates that the reaction rate constant
does not follow Arrhenius equation and the Ea was fixed to zero during parameter estimation."*)

**(3) ⚠️ THE Ea VALUES THEMSELVES ARE NOT IN THIS DOCUMENT.** They live in **Table S4 of the
Supplementary Material, which is NOT in the PDF on disk** (the file ends at the Highlights page,
p.53). **I am not reporting numbers I did not read.** ⇒ **The corpus holds ZERO hazelnut activation
energies. `research_round3_channels.md` §A.2 row 3 correctly ranked this paper low for Ea; it turns
out to be worse than "no usable Ea" — it is "no Ea at all in the file we hold."**

⇒ **RECORD AS A PROHIBITED DERIVATION, in the register of `k3` §C.1 and round-3 §A.3.4: no
activation energy may be attributed to `10.1016/j.foodchem.2016.11.159` from this file, and if
Table S4 is ever obtained, its values are disqualified by the authors' own disclaimer and by the
1174 kJ/mol bound.**

---

## §2. SYSTEM AND CONDITIONS `[M]`

| item | value | anchor |
|---|---|---|
| matrix | **hazelnuts** (*Corylus avellana* L.), raw | §2 |
| **roasting** | **150, 160 and 170 °C** for **15, 30, 60, 90 and 120 min** | §3.1 |
| model units | mg/kg dw and g/100 g dw; rate constants in `min⁻¹×10³` and `kg·µmol⁻¹·min⁻¹×10³` | Table 1 |
| software | Athena Visual Studio v14.2, non-linear regression, determinant criterion (Van Boekel 1996), HPD intervals | §2.6 |
| Arrhenius form | reparameterised, **T_b = 160 °C**, `R = 8.3145×10⁻³ kJ/mol K`, *"the data at 150, 160, and 170 °C were fitted to experimental data all at once"* | §2.6, verbatim |
| **pH** | ⚠️ **measured but reported only in Table S1 (not in the PDF).** All that is stated in the body: *"the moisture content and pH of the hazelnuts **decreased** during roasting (Table S1)"* | §3.1 |
| **moisture / a_w** | same — **Table S1 only, not held** | §3.1 |

### 2.1 Initial reactant loads and absolute product yields — the strongest part of the paper `[M]`

| quantity | value | anchor |
|---|---|---|
| **initial sucrose** | **5.5 ± 0.1 g/100 g dw** (*"80–90 % of total sugars"* in hazelnut, `[C]` Alasalvar 2010) | §3.1 |
| initial fructose | **0.4 ± 0.02 g/100 g dw** | §3.1 |
| initial glucose | **0.2 ± 0.06 g/100 g dw** | §3.1 |
| sucrose loss (120 min) | **60 %, 70 %, 90 %** at 150 / 160 / 170 °C | §3.1 |
| initial total free amino acids | **2112 ± 49 mg/kg dw** | §3.1 |
| initial protein-bound lysine | **5401 ± 50 mg/kg dw** | §3.1 |
| **initial total amino acids ("AA")** | **7513 ± 87 mg/kg dw** | §3.1 ⇒ ⚠️ **same free + bound lump as the wheat-flour paper** |
| total AA loss (120 min) | **68 %, 81 %, 85 %** at 150 / 160 / 170 °C | §3.1 |
| ★ **HMF at 120 min** | ★ **104 ± 0.5, 238 ± 1.9, 278 ± 0.7 mg/kg dw** at 150 / 160 / 170 °C | §3.1 |
| 3-DG at 120 min | **6.7 ± 0.1 (150 °C), 6.1 ± 0.1 (160 °C)** mg/kg dw; at **170 °C an apparent maximum of 5.4 ± 0.1 at 60 min, then decreasing** | §3.1 |
| **3,4-DG** | *"approximately **5 times lower** than those of 3-deoxyglucosone"*, with **a short lag phase** | §3.1 |
| 1-DG max | **0.22 ± 0.01, 0.31 ± 0.03, 0.27 ± 0.01 mg/kg dw** at 150 / 160 / 170 °C — *"the least abundant"* deoxyhexulose | §3.1 |
| **glyoxal, raw hazelnut** | ★ **1.7 ± 0.6 mg/kg dw — the only α-dicarbonyl present before roasting**; rises ~4× in 15 min then plateaus at all temperatures | §3.1 |
| MGO max | **6.6 ± 0.5 mg/kg dw** at 160 °C × 90 min | §3.1 |
| DMG | forms at all temperatures; **no significant difference (p < 0.05) between 150 and 160 °C at 120 min** | §3.1 |
| **glucosone** | ★ **NOT PRESENT at any time–temperature** — *"Glucosone was not present in hazelnuts roasted at selected time-temperature combinations."* | §3.1 |
| regulatory context | EFSA dietary HMF intake estimate **1.6 mg/person/day**; a 30 g hazelnut portion can exceed it at long roasting | §3.1 `[C]` |
| industry reference | *"the generally used thermal treatment by industry of **145 °C for 15 min**, HMF will not be a health risk"* | §3.1 |

★ **`HMF = 104 / 238 / 278 mg/kg dw` at 120 min across three temperatures, with authentic-standard
HPLC quantification, in a named real matrix — a clean absolute anchor and a strong hold-out
candidate.** (Şen 2022 §3.2 independently cites this paper's 278.0 ± 0.7 mg/kg dw figure, so it is
also the one number in this cluster that is corroborated across two papers of the same group.)

⚠️ **Comparison caveat carried by the authors themselves**, §3.1: methylglyoxal and glyoxal in
hazelnut are *"almost 10 times"* higher than in olive oil heated at 200 °C for 1 h — *"Therefore,
most of the methylglyoxal and glyoxal formed in hazelnuts during roasting was considered to
originate from fractions other than lipids."* ★ A useful **lipid-vs-carbohydrate source
attribution** for the dicarbonyl lane.

---

## §3. THE REACTION NETWORK — 26 steps (Figure 2b), and the HMF ODE

Abbreviations verbatim (Table 1 footnote a): `SUC` sucrose; `GLC` glucose; `FRU` fructose;
**`FFC` fructofuranosyl cation**; `1,2-ED` 1,2-enediol; **`AP` Amadori product**;
**`HP` Heyns product**; `1-DG`; `3-DG`; `3,4-DG` 3,4-dideoxyglucosone; `GO` glyoxal;
`MGO` methylglyoxal; `DMG` dimethylglyoxal; `HMF`; `AA` total amino acids; `P` products.

```
 1  SUC → GLC + FFC        14  AP → 3-DG
 2  GLC → 1,2-ED           15  AP → 1-DG
 3  1,2-ED → GLC           16  1-DG → MGO
 4  GLC → 3-DG             17  1-DG → DMG
 5* GLC + AA → AP          18  3-DG → 3,4-DG
 6  GLC → GO               19  3,4-DG → HMF   ★ terminal, 3-DG LIMB
 7  1,2-ED → FRU           20  HP → P1
 8  FRU → 1,2-ED           21  AP → P2
 9* FRU + AA → HP          22  1-DG → P3
10* FFC + AA → HP          23  MGO → P4
11  FFC → HMF   ★ terminal, FRUCTOSE/CATION LIMB     24  DMG → P5
12  HP → 1-DG              25  GO → P6
13  HP → 3-DG              26  HMF → P7   ★ THE HMF SINK
```
(`*` = bimolecular, units `kg×µmol⁻¹×min⁻¹×10³`.) A 27th step, **FRU → FFC**, was tested and
**removed**: *"Contribution of free fructose to the fructofuranosyl cation formation was not
included in the model (Fig. 2b) as the rate constants (k₂₇) were found to equal to zero in the
comprehensive model (Fig. 2a), which means it was kinetically less significant"* (§3.5).

**★ THE HMF NODE ODE (Appendix, decoded from the subscript-stripped text layer against the Fig. 2b
edge list and Table 1 numbering):**

> **`d[HMF]/dt = k₁₁·[FFC] + k₁₉·[3,4-DG] − k₂₆·[HMF]`**

**Structurally identical to both Kocadağlı papers**: two parallel first-order sources, one
first-order sink. **No HMF + amine step** (unlike Şen 2022).

⚠️ **`FFC`, `AP`, `HP` and `1,2-ED` are ALL UNMEASURED.** The authors name this as the cause of
their worst uncertainties, verbatim §3.2:

> "reaction rate constants of some steps (**k₁₃, k₁₆, k₂₀, k₂₂, k₂₅**) showed a higher uncertainty
> and could not be estimated in this interval. The reason for that might be the **involvement of
> compounds like Heyns product or degradation products of dicarbonyl compounds, which could not be
> quantified analytically**, in the model. However, exclusion of these compounds from the model …
> **did not give better solutions.**"

---

## §4. TABLE 1 — the complete 26-step × 3-temperature constant set

**Transcribed verbatim from p.51–52 of the accepted manuscript.** Caption verbatim: *"Reaction rate
constants with 95 % highest posterior density (HPD) intervals at different temperatures according to
the proposed kinetic model in Figure 2b for Maillard reaction and caramelization during roasting of
hazelnuts."* Footnote `ᵇ ind`: *"indeterminate, which means a large uncertainty in the estimated
parameter within 95 % confidence interval."* Provenance **`[F]`**.

| # | elementary step | units | k(150) | HPD | k(160) | HPD | k(170) | HPD |
|---|---|---|---|---|---|---|---|---|
| 1 | SUC → GLC + FFC | min⁻¹×10³ | 6.9 | ±0.8 | 15 | ±3.1 | 22 | ±2.0 |
| 2 | GLC → 1,2-ED | min⁻¹×10³ | 141 | ±27.0 | 473 | ±131 | 698 | ±178 |
| 3 | 1,2-ED → GLC | min⁻¹×10³ | **0** | ±0 | 8.5 | ±3.6 | 28 | ±8.3 |
| 4 | GLC → 3-DG | min⁻¹×10³ | 0.03 | ±0.02 | **0** | ±0 | **0** | ±0 |
| 5\* | GLC + AA → AP | kg·µmol⁻¹·min⁻¹×10³ | 0.0009 | ±0.0007 | 0.003 | ±0.001 | 0.009 | ±0.001 |
| 6 | GLC → GO | min⁻¹×10³ | 0.6 | ±0.2 | 2.5 | ±0.9 | 9.2 | ±1.0 |
| 7 | 1,2-ED → FRU | min⁻¹×10³ | 1.3 | ±0.9 | 1.8 | ±0.7 | 4.2 | ±1.9 |
| 8 | FRU → 1,2-ED | min⁻¹×10³ | **0** | ±0 | **0** | ±0 | 41 | ±14 |
| 9\* | FRU + AA → HP | kg·µmol⁻¹·min⁻¹×10³ | 0.00023 | ±0.00007 | 0.00062 | ±0.00015 | **0** | ±0 |
| 10\* | FFC + AA → HP | kg·µmol⁻¹·min⁻¹×10³ | 0.00094 | ±0.00028 | 0.00027 | ±0.00030 | 0.00004 | ±0.00002 |
| **11** | **FFC → HMF** | min⁻¹×10³ | **0.58** | **±0.13** | **0.57** | **±0.19** | **2.02** | **±1.32** |
| 12 | HP → 1-DG | min⁻¹×10³ | 0.23 | ±0.18 | 0.85 | ±0.58 | 267 | ±177 |
| 13 | HP → 3-DG | min⁻¹×10³ | 0.009 | ±0.004 | 0.022 | ±0.030 | 12 | **indᵇ** |
| 14 | AP → 3-DG | min⁻¹×10³ | **0** | ±0 | 0.62 | ±0.03 | **0** | ±0 |
| 15 | AP → 1-DG | min⁻¹×10³ | 3.47 | ±3.12 | 3.51 | ±1.33 | 0.56 | ±0.58 |
| 16 | 1-DG → MGO | min⁻¹×10³ | 7012 | ±5510 | 33016 | **indᵇ** | 47920 | ±33240 |
| 17 | 1-DG → DMG | min⁻¹×10³ | 371 | ±241 | 895 | ±581 | 1073 | ±618 |
| 18 | 3-DG → 3,4-DG | min⁻¹×10³ | 4.27 | ±0.57 | 29.4 | ±26.1 | 88.1 | ±22.7 |
| **19** | **3,4-DG → HMF** | min⁻¹×10³ | **0** | **±0** | **134** | **±127** | **390** | **±111** |
| 20 | HP → P1 | min⁻¹×10³ | 11 | ±1.4 | 4.7 | ±5.6 | 59 | **indᵇ** |
| 21 | AP → P2 | min⁻¹×10³ | 140 | ±122 | 21.2 | ±34.4 | 5.24 | ±1.75 |
| 22 | 1-DG → P3 | min⁻¹×10³ | 122 | **indᵇ** | **0** | **indᵇ** | 404 | **indᵇ** |
| 23 | MGO → P4 | min⁻¹×10³ | 126 | ±106 | 737 | ±113 | 918 | ±639 |
| 24 | DMG → P5 | min⁻¹×10³ | 54 | ±41 | 130 | ±90 | 106 | ±65 |
| 25 | GO → P6 | min⁻¹×10³ | 18 | ±8.3 | 61 | ±25 | 290 | **indᵇ** |
| **26** | **HMF → P7** | min⁻¹×10³ | **12** | **±3.7** | **21** | **±11** | **103** | **±63.7** |

**Cross-check of the sucrose row against the running text `[D]`:** §3.1 states *"The degradation
rates of sucrose to glucose and fructofuranosyl cation were **0.0069, 0.015 and 0.022 min⁻¹**"* —
Table 1 step 1 gives 6.9 / 15 / 22 ×10⁻³ min⁻¹ ✔ **exact match; the units convention is confirmed.**

### 4.1 Cells that carry no information

**5 `ind` cells** (13@170, 16@160, 20@170, 22 all three, 25@170).
**HPD ≥ estimate:** step 3@150 (0±0), 4@160/170 (0), 8@150/160 (0), 9@170 (0), 13@160 (0.022 ±
0.030), 14@150/170 (0), 15@150 (3.47 ± 3.12) and @170 (0.56 ± 0.58), **16@150 (7012 ± 5510, 79 %)
and @170 (47920 ± 33240, 69 %)**, 18@160 (29.4 ± 26.1, 89 %), **19@160 (134 ± 127, 95 %)**,
20@160 (4.7 ± 5.6), 21@150 (140 ± 122) and @160 (21.2 ± 34.4), 23@150 (126 ± 106) and @170
(918 ± 639), 24@150 (54 ± 41) and @160 (130 ± 90) and @170 (106 ± 65), **26@170 (103 ± 63.7, 62 %)**,
11@170 (2.02 ± 1.32, 65 %).
**⇒ ~28 of 78 cells (36 %) are unusable — the highest defect rate in the K5a cluster.** The
contamination reaches **both** HMF-forming steps and the sink.

---

## §5. ★★ THE HMF BRANCH — where the abstract and the table disagree, and how the authors resolve it

### 5.1 The abstract and Highlights

Abstract, verbatim (lines 20–22): *"**5-hydroxymethylfurfural formation mainly proceeded via
fructofuranosyl cation dehydration rather than 3-deoxglucosone**"* (`deoxglucosone` sic).
Highlight #2, verbatim: *"**5-Hydroxymethylfurfural formation proceeds via fructofuranosyl cation
dehydration.**"*

### 5.2 ★ What Table 1 actually says about the two terminal steps `[D]`

| | k₁₁ FFC→HMF | k₁₉ 3,4-DG→HMF | **k₁₉ / k₁₁** |
|---|---|---|---|
| 150 °C | 0.58 ± 0.13 | **0 ± 0** | 0 |
| 160 °C | 0.57 ± 0.19 | **134 ± 127** | **235×** |
| 170 °C | 2.02 ± 1.32 | **390 ± 111** | **193×** |

**At 160 and 170 °C the 3-DG limb's terminal rate constant is 193–235× LARGER than the fructose
limb's.** ★ **This is the same apparent contradiction found in Gürsul Aktağ 2020 (§1.2 of that
dossier) — and it recurs in the same lab's Şen 2022, where the 3-DG constant is again the larger.
In three of the five Gökmen multiresponse papers, the fructose-limb *headline* is accompanied by a
smaller fructose-limb *rate constant*.**

### 5.3 ★★ AND HERE — UNLIKE IN GÜRSUL AKTAĞ — THE AUTHORS SEE IT AND EXPLAIN IT

§3.5, verbatim, in full because this is the decisive methodological passage of the whole cluster:

> "The role of fructofuranosyl cation in 5-hydroxymethylfurfural formation during hazelnut roasting
> was tested by **excluding this reaction step from the proposed model. The predicted values of
> 5-hydroxymethylfurfural were found to be far below the experimental values and indicated only the
> contribution of the 3-deoxyglucosone pathway (Fig. 4).** Therefore, contribution of
> 5-hydroxymethylfurfural formation through the dehydration of fructofuranosyl cation was found to
> be crucial compared to formation through the 3-deoxyglucosone pathway during roasting of hazelnut.
> The reason for a smaller contribution of 3-deoxyglucosone to 5-hydroxymethylfurfural formation was
> due to lower reaction rate constants of 3-deoxyglucosone formation from glucose, Heyns product and
> Amadori product (**k₄, k₁₃, and k₁₄**) (Table 1). **The rate constants of 3,4-dideoxyglucosone
> formation from 3-deoxyglucosone (k₁₈) and 5-hydroxymethylfurfural formation from
> 3,4-dideoxyglucosone (k₁₉) were higher than the rate constants of 5-hydroxymethylfurfural
> formation through fructofuranosyl cation (k₁₁).** ★ **The lower reaction rate constants of
> 5-hydroxymethylfurfural formation through fructofuranosyl cation could be attributed to the fact
> that concentration of fructofuranosyl cation could not be measured because of experimental
> restrictions and the reaction steps of 5-hydroxymethylfurfural formation through fructofuranosyl
> cation were reduced to one dehydration step, in order to simplify the proposed model.**"

**★ THIS IS THE SINGLE MOST IMPORTANT PARAGRAPH THE K5a WAVE FOUND, AND IT SHOULD GOVERN THE HMF
NODE DESIGN.** In plain terms, the authors state:

1. **The fructose-limb claim rests on MODEL DISCRIMINATION (delete the edge → the fit collapses),
   NOT on the magnitude of the rate constant.**
2. **The 3-DG limb's rate constants are the larger ones, and they say so explicitly.**
3. **`k₁₁` is small BECAUSE `[FFC]` is unmeasured** — only the product `k₁₁·[FFC]` is constrained,
   and the fitted `k₁₁` absorbs whatever scale the unmeasured pool takes. **A rate constant on an
   unmeasured node is not commensurable with one on a measured node.**
4. `k₁₁` additionally **lumps a multi-step dehydration into one edge**, so it is not an elementary
   constant at all.

**⇒ THE "FRUCTOSE LIMB vs 3-DG LIMB" QUESTION IS NOT A CONTEST BETWEEN RATE CONSTANTS. It is a
contest between FLUXES, and the fluxes are set by the precursor pools. Every "the fructose limb
dominates" statement in this literature is a flux statement about a matrix; every rate-constant
comparison that appears to contradict it is comparing a measured node with an unmeasured one.**

### 5.4 The 3-DG limb's internal rate-determining step `[F]`

§3.5, verbatim: *"The reaction rate constant of 5-hydroxymethylfurfural formation from
3,4-dideoxyglucosone was found to be **almost 5 times higher** than the rate constant of
3,4-dideoxyglucosone formation from 3-deoxyglucosone. In the 5-hydroxymethylfurfural formation
through 3-deoxyglucosone pathway, **5-hydroxymethylfurfural formation from 3,4-dideoxyglucosone was
found to be the fast step and 3,4-dideoxyglucosone formation from 3-deoxyglucosone was the rate
determining step.**"

Check `[D]`, k₁₉/k₁₈: **150 °C → 0/4.27 = 0 (the claim fails); 160 °C → 134/29.4 = 4.56 ✔;
170 °C → 390/88.1 = 4.43 ✔.** The "almost 5 times" holds at 160 and 170 °C only.
★ **Independent corroboration: the wheat-flour paper gives k₁₃/k₁₂ = 2.17 / 3.50 / 2.18 at
160/180/200 °C — same ordering, same lab, different matrix, factor 2–5 rather than ~5.**
⇒ **`3-DG → 3,4-DG` is the rate-determining step of the 3-DG limb in TWO independent matrices.**
This is a genuinely transferable structural constraint.

**⚠️ AND IT DIRECTLY OPPOSES Gürsul Aktağ 2020 §3.2(ii)** for the *fructose* limb: there,
*"formation of FFC from fructose was found to be the fast step and the rate determining step was the
HMF formation from FFC."* **The two limbs have their rate-determining step in opposite positions.**

---

## §6. ★ WAVE K5a CROSS-CHECK — a 3-point Arrhenius refit of Table 1 `[D]`

Since Table S4 is not held, this refit is the only temperature-dependence information available.
OLS `ln k` vs `1/T` over 423.15 / 433.15 / 443.15 K.

| step | Ea refit (kJ/mol) | R² | flag |
|---|---|---|---|
| 1 SUC→GLC+FFC | 90.6 | 0.968 | ✔ |
| 2 GLC→1,2-ED | 125.2 | 0.926 | ✔ |
| 6 GLC→GO | **212.8** | 1.000 | high but log-linear |
| **11 FFC→HMF** | **96.5** | **0.728** | ⚠️ **non-monotonic (0.58 → 0.57 → 2.02)** |
| 16 1-DG→MGO | 150.5 | 0.897 | ⚠️ but 2 of 3 cells are `ind`/HPD-dominated |
| 17 1-DG→DMG | 83.2 | 0.882 | HPD-dominated at 150 |
| 18 3-DG→3,4-DG | **236.4** | 0.979 | ⚠️ far above van Boekel's 120 kJ/mol norm |
| **19 3,4-DG→HMF** | **NOT COMPUTABLE** | — | ★ **k = 0 exactly at 150 °C ⇒ ln k = −∞** |
| 21 AP→P2 | **−256.3** | 0.995 | ⚠️ **monotone DECREASING (140→21.2→5.24)** |
| 23 MGO→P4 | 155.7 | 0.842 | |
| 24 DMG→P5 | 53.2 | 0.551 | non-monotonic |
| 25 GO→P6 | 216.4 | 0.993 | |
| **26 HMF→P7** | **166.9** | **0.922** | ⚠️ HPD is 62 % of the estimate at 170 °C |

**Refit range: −256 to +236 kJ/mol** — and **five** further steps (3, 4, 8, 9, 14, 19) contain an
exact zero and admit **no** Arrhenius line at all. That is at minimum **six** underivable Ea, which
matches the authors' *"six zero"* count in §3.2 — **so the six zeros in Table S4 are almost
certainly these steps, set to zero because ln(0) is undefined, not because the chemistry is
athermal.** ⇒ **The paper's own "0–1174 kJ/mol" range is bounded below by a numerical artefact and
above by an absurdity. Neither end is chemistry.**

**★ Note the direction of the failure vs the JAFC companion.** In the amine-free glucose melt,
`Fru→Int` and `Int→HMF` were the *only* steps whose Arrhenius behaviour reproduced (R² = 1.000).
Here, in a real nut matrix at 150–170 °C, **`FFC→HMF` is non-monotonic (R² = 0.728) and
`3,4-DG→HMF` is undefined.** Together with the wheat-flour result (§4 of that dossier: `Int→HMF`
R² = 0.322, Ea = −97.6), **all three real-matrix systems destroy the temperature dependence that
the amine-free model system supported.**

---

## §7. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL` — the paper's real contribution

| # | constraint | anchor |
|---|---|---|
| S1 | ★ **Deleting the FFC→HMF edge makes the model predict HMF "far below the experimental values"** | §3.5, Fig. 4 — the same test the JAFC paper (Fig. S5) and the wheat-flour paper (Fig. 5b) ran, with the same outcome. **Three independent matrices, one conclusion: a 3-DG-only HMF node cannot reproduce observed HMF.** |
| S2 | ★ **`3-DG → 3,4-DG` is the rate-determining step of the 3-DG limb; `3,4-DG → HMF` is fast** — factor ~5 here, factor 2–3.5 in wheat flour | §3.5, §5.4 |
| S3 | ★ **`FRU → FFC` is kinetically negligible under dry roasting (k₂₇ = 0 in the comprehensive model)** — FFC comes from **sucrose cleavage**, not from free fructose | §3.5 ⇒ **In a sucrose matrix the "fructose limb" is really a SUCROSE limb. This is why Şen 2022's network has no fructose node at all.** |
| S4 | **1,2-enolisation is a rate-determining step**; omitting it makes glucose and fructose *"continuously decrease"*, contradicting the observed plateau | §3.3 — agrees with the wheat-flour paper, opposes the JAFC amine-free paper |
| S5 | **No mannose was detected**; epimerisation is negligible | §3.3 |
| S6 | **Glucose contributes more than fructose and FFC to the early stage of the MR** | abstract |
| S7 | **MGO and DMG both come from 1-DG** | abstract, §3.6 |
| S8 | **GO comes from glucose degradation** (k₆), not from glucosone — **glucosone does not exist in this matrix** | Highlight #4, §3.1 |
| S9 | **Degradation-rate ordering: MGO (k₂₃) > DMG (k₂₄) > GO (k₂₅) > HMF (k₂₆)**, and *"they all increased with roasting temperature"* | §3.6, verbatim. ★ **HMF is the LEAST reactive of the four toward whatever consumes them** — check `[D]` at 170 °C: 918 > 106 (fails: DMG < GO 290) — ⚠️ **the ordering holds at 150 and 160 °C but DMG and GO swap at 170 °C.** Quote with temperature. |
| S10 | Rate constants of glucose/fructose degradation (k₂, k₄, k₆, k₈) **exceed** those of amino-acid regeneration from AP/HP (k₁₂–k₁₅) | §3.1 |
| S11 | **HMF is a stable end product that accumulates**, with a first-order sink only | ODE, §3.1 |
| S12 | Most MGO and GO in roasted hazelnut originate from **non-lipid** fractions, despite ~60 % oil | §3.1 (the olive-oil comparison) |

---

## §8. VERIFIED NEGATIVES `[NEG]`

- **★ NO ACTIVATION ENERGY IS IN THIS FILE.** Table S4 (Supplementary) is not held. The authors
  disclaim Arrhenius twice (abstract + conclusion) and report a range including 1174 kJ/mol.
- **No pH value, no moisture value, no water activity value** — all in Table S1, not held.
- **No colour / browning / melanoidin measurement** (colour is in Table S3, not held).
- **No glucosone** (not detected).
- **No furfural.**
- **No acrylamide.**
- **Supplementary Material — Tables S1–S4 and Figures S2–S10 — is entirely absent from the PDF.**
  This includes the model-discrimination figures S3–S10 that underpin the exclusion decisions.
- **AP, HP, FFC and 1,2-ED are never quantified** (author-declared, §3.2).

---

## §9. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| **any activation energy for this DOI** | — | ★ **REFUSE. Not in the file; author-disclaimed; the reported range (0–1174 kJ/mol) is disqualified at both ends. RECORD AS A PROHIBITED DERIVATION.** |
| Table 1, all 78 k ± HPD cells | `[F]` | **USE-Q** — `min⁻¹×10³` (or `kg·µmol⁻¹·min⁻¹×10³` for steps 5, 9, 10), per-temperature, hazelnut only; **drop the ~28 cells of §4.1 (36 %)** |
| **k₁₉ (3,4-DG→HMF)** | `[F]` | **USE-Q at 170 °C only** (390 ± 111). **REFUSE at 150 °C (exactly 0) and 160 °C (134 ± 127, HPD 95 % of estimate).** |
| **k₁₁ (FFC→HMF)** | `[F]` | **RATIO-ONLY — `[FFC]` unmeasured and the step is a lumped multi-dehydration; author-declared (§5.3).** Non-monotonic. |
| **k₂₆ (HMF sink)** | `[F]` | **USE-Q at 150 and 160 °C** (12 ± 3.7; 21 ± 11); **REFUSE at 170 °C** (103 ± 63.7). ⚠️ **This is the only HMF-sink constant in the whole cluster with a usable HPD at more than one temperature** — see the synthesis §4. |
| k₁ (sucrose) 6.9/15/22 ×10⁻³ min⁻¹ | `[F]` | ★ **USE** — tight HPD at all three T, cross-verified against the running text |
| k₁₈ (3-DG→3,4-DG) | `[F]` | **USE-Q** — 150 and 170 °C only (HPD 89 % at 160 °C) |
| Ea refits of §6 | `[D]` | **PRIOR-ONLY, and only for steps 1, 2, 6, 25** (monotone, R² ≥ 0.93). **Never for 11, 19, 21, 24, 26.** |
| **HMF = 104 ± 0.5 / 238 ± 1.9 / 278 ± 0.7 mg/kg dw @ 120 min, 150/160/170 °C** | `[M]` | ★ **USE — proposed HOLD-OUT.** Independently re-cited by Şen 2022. |
| initial sucrose/fructose/glucose, AA loads and losses, 3-DG / 3,4-DG / 1-DG / MGO / GO / DMG levels | `[M]` | **USE** |
| **3,4-DG ≈ 3-DG / 5**, and **glyoxal 1.7 ± 0.6 mg/kg in RAW hazelnut** | `[M]` | ★ **USE — a rare pre-heating baseline** |
| S1–S12 | `[M]`/`[F]` | **STRUCTURAL** — **S1, S2 and S3 are the highest-value rows in the whole K5a wave** |
| §5.3 (the unmeasured-node admission) | `[M]` | ★ **STRUCTURAL — governs the HMF node design. Quote it into the declaration if the node ships.** |
