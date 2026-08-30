# Gürsul Aktağ & Gökmen 2020 — α-dicarbonyls and HMF in fruit juices during storage

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

---

## §0. IDENTITY

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/Gursul2020.pdf` (585,691 bytes) | `ls` |
| **DOI** | **`10.1016/j.foodchem.2020.126620`** ✔ matches the round-3 expectation | printed on p.1 of the PDF |
| title | *"Multiresponse kinetic modelling of α-dicarbonyl compounds formation in fruit juices during storage"* | p.1 |
| authors | **Işıl Gürsul Aktağ**, Vural Gökmen | p.1 |
| journal | *Food Chemistry* **320** (2020) 126620 | running head |
| dates | Received 26 Oct 2019 · Revised 17 Feb 2020 · Accepted 15 Mar 2020 · online 17 Mar 2020 | p.1 |
| version | **published version of record** (Elsevier typeset, `0308-8146/© 2020 Elsevier Ltd.`) | p.1 |

Filename/DOI mapping is **correct** for this file.

---

## §1. ★★ THE SINGLE MOST IMPORTANT FINDING OF WAVE K5a: THE ROUND-3 PREMISE FOR THIS PAPER IS WRONG

`research_round3_channels.md` §A.2 row 5 records this paper as:

> "**The only acidic, food-pH, low-temperature member of the set** … **the only one where the
> fructose-vs-3-DG partition was tested for statistical significance**"

The first half is correct. **The second half is not, and the abstract sentence it rests on does not
survive contact with the paper's own tables.** Three independent defects, all verifiable from the
PDF:

### 1.1 The Methods describe **no** significance test on any rate constant

The paper's entire statistics section, verbatim (§2.8, p.4):

> "**Free amino acid data** were subjected to analysis of variance (one-way ANOVA). The SPSS 17.0
> statistical package was used for the evaluation of statistical significance of the differences
> between mean values by Tukey test. P < 0.05 was considered to be statistically significant for
> the results."

**The ANOVA is declared for the free amino acid data only.** No test — ANOVA, *t*, bootstrap,
posterior overlap, likelihood ratio — is described anywhere for the rate constants. The only
uncertainty attached to any *k* in this paper is the Bayesian 95 % HPD interval from Athena Visual
Studio. Elsewhere the paper uses "significantly" loosely for HPD non-overlap, e.g. §3.2(iii):
*"The rate of glucosone formation through 1,2-enediol intermediate (k7) increased (**± 95 % HPD**)
with increasing storage temperature"*. **⇒ The abstract's "(p < 0.05)" is not backed by a described
procedure.**

### 1.2 At the terminal HMF-forming step, the table says the **opposite** of the abstract

Abstract, verbatim: *"The contribution of fructose dehydration through fructofuranosyl cation on
the formation of 5-hydroxymethylfurfural was significantly higher (p < 0.05) than 3-deoxyglucosone
dehydration."*

The two steps that make HMF are **k₈ (3-DG → HMF)** and **k₁₁ (FFC → HMF)**. From the paper's own
Table 1 (×10⁻³ week⁻¹), at 37 °C — the only temperature where HMF accumulates at all:

| juice | k₈ 3-DG→HMF | k₁₁ FFC→HMF | which is larger |
|---|---|---|---|
| Apple | 0.7 ± **1.98** | 5.3 **ind\*** | k₁₁, but **k₈'s HPD is 2.8× its estimate and k₁₁ is author-declared indeterminate** |
| **Orange** | **2.0 ± 0.68** | **0.0 ± 0.02** | ★ **k₈ — the 3-DG limb, by ∞** |
| **Peach** | **11.8 ± 2.79** | **0.4 ± 0.22** | ★ **k₈ — the 3-DG limb, by 29.5×** |

**In two of the three matrices the 3-DG→HMF rate constant is larger than the FFC→HMF rate constant,
and in one of those the FFC step is exactly zero.** In the third (apple) the only number favouring
the fructose limb is the one the authors themselves mark `ind` — *"indeterminate, which means a
large uncertainty in the estimated parameter within 95 % confidence interval"*.

### 1.3 The "7, 4 and 2 times higher" claim in the body **mixes two different reaction steps**

Body, §3.2(ii), verbatim: *"Comparing to 3-DG, HMF formation from fructose were found almost **7, 4
and 2 times higher** in apple juice, orange juice and peach nectar, respectively."*

I reproduced these three ratios from Table 1 `[D]`:

| juice | claimed | k₁₁/k₈ | k₁₀/k₈ | which one matches |
|---|---|---|---|---|
| Apple | ~7 | 5.3/0.7 = **7.6** ✔ | 0.2/0.7 = 0.29 ✘ | **k₁₁/k₈** (the HMF-forming step) |
| Orange | ~4 | 0.0/2.0 = **0** ✘ | 7.1/2.0 = **3.6** ✔ | **k₁₀/k₈** (FRU→FFC — **not** an HMF-forming step) |
| Peach | ~2 | 0.4/11.8 = 0.034 ✘ | 22.6/11.8 = **1.9** ✔ | **k₁₀/k₈** |

**⇒ The headline ratio set is assembled from `k₁₁/k₈` for apple and `k₁₀/k₈` for orange and peach.**
`k₁₀` is `FRU → FFC`, the *entry* to the fructose limb, not the step that makes HMF. And the one
cell that does use the HMF-forming step (apple `k₁₁ = 5.3`) is the cell the authors flagged as
indeterminate. Consistently using `k₁₁/k₈` gives **7.6 / 0 / 0.034** — i.e. the fructose limb wins
in one matrix out of three.

### 1.4 What the paper *does* legitimately establish about the branch

The defensible, in-paper statement is the narrower one in §3.2(ii), verbatim:

> "The rate constants of fructose dehydration to FFC (k₁₀) were found significantly higher than
> that of FFC dehydration to HMF (k₁₁) in orange juice and peach nectar (Table 1). For apple juice,
> the rate constant of HMF formation from FFC (k₁₁) showed a large uncertainty in the estimated
> parameter within 95 % confidence interval. … **In HMF formation through fructose pathway,
> formation of FFC from fructose was found to be the fast step and the rate determining step was
> the HMF formation from FFC.**"

✔ Verified from Table 1 `[D]`: orange 7.1 > 0.0; peach 22.6 > 0.4. **That is a rate-determining-step
claim about the fructose limb's internal structure, not a partition claim against 3-DG.**

And the genuinely strong structural argument the paper makes is a *pool-size* argument, verbatim:

> "It should be noted here that **fructose as a precursor leading to HMF was the dominant form of
> sugars in these samples.** Initial concentration of fructose was **2.3 times and 9.3 times higher
> in apple juice than orange juice and peach nectar** used in this study."

⇒ **In apple juice the fructose limb wins on POOL SIZE, not on rate constant.** That is a
matrix-composition statement and it is transferable as such.

**★ RECOMMENDED CORRECTION TO `research_round3_channels.md` §0a and §A.2 row 5:** delete "the only
one where the fructose-vs-3-DG partition was tested for statistical significance"; replace with
"the abstract asserts a significant fructose-limb advantage, but the paper describes no test on
rate constants, its terminal HMF-forming constant favours the 3-DG limb in 2 of 3 juices, and the
supporting ratio set mixes two different reaction steps (K5a §1)."

---

## §2. SYSTEM AND CONDITIONS `[M]`

| item | value | anchor |
|---|---|---|
| matrices | **apple juice** (clear), **orange juice** (cloudy), **peach nectar** (added-sugar) | §2.2 |
| production | laboratory-made from Golden Delicious apples, Washington oranges, Bursa peaches (Fig. 1 flowcharts). Apple: heat 85 °C×5 min → press → pectinase 1 mL/L + amylase 0.2 mL/L at 50 °C → gelatine-bentonite clarification → filter. Orange: pulper only. Peach: heat 85 °C×5 min → pulper → nectar (**66 % sugar syrup**, citric acid, water, pulp) | Fig. 1 |
| pasteurisation | **90 °C × 10 min** in a water bath, in falcons, for all three | Fig. 1 |
| **storage temperatures** | **4, 27 and 37 °C** | §2.2 |
| **storage duration** | **24 weeks**, sub-sampled **3 parallel every 2 weeks**, frozen −18 °C until analysis | §2.2 |
| **pH** | ★ **3.4**, and **it does not change during storage** — *"no change was observed in the pH values of samples during storage in our study (Table S3)"* | §3.2(i) and §3.1 |
| °Brix (initial) | apple **15°**, orange **12°**, peach **16°** | §3.1 |
| model units | **mmol/L** | §2.8 |
| ★ **kinetics restricted to 27 and 37 °C** | *"there were no significant changes in concentrations of reactants and reaction products in the samples stored at 4 °C (Table S1). Therefore, kinetic analysis was limited with the data observed for 27 and 37 °C."* | §3.1, verbatim |

### 2.1 ★ THE STRUCTURAL GIFT: THE MAILLARD REACTION IS ABSENT, MEASURED `[M]`

Verbatim, §3.1:

> "the concentrations of **free amino acids remained relatively stable in all samples during storage
> at all temperatures** (Table S2) … **Our findings do not support that Maillard reaction takes
> place in apple juice, orange juice and peach nectar during storage.** Comparing to acid
> conditions, Maillard reaction occurs more quickly under neutral conditions … As acidic and
> high-water activity systems, fruit juices are not considered suitable for the Maillard reaction to
> take place during storage. **Therefore, formation of α-dicarbonyl compounds in juices and nectars
> during storage would be due to sugar dehydration reactions.**"

★ **This is a measured amine-off control at food pH and food temperature** — the counterpart of the
Kocadağlı JAFC amine-free melt at 160–200 °C. **The two together bracket the caramelization-only
HMF route from 4 °C to 200 °C.** The claimed absence rests on the ANOVA of §1.1, which *is* the test
the Methods describe; this is the one significance claim in the paper that is properly supported.

### 2.2 Analytical quality `[M]`

| analyte | method | calibration | LOD / LOQ | verdict |
|---|---|---|---|---|
| sucrose, glucose, fructose | HPLC-RID, Shodex SH-1011, 5 mM H₂SO₄, 50 °C | external, 0.25–2.5 g/L | 1.0 / 3.0 mg/L | absolute |
| **HMF** | Shimadzu UFLC, Atlantis dC18, 10 mM formic acid/MeCN 90:10, DAD **285 nm** | external, **1–10 mg/L** | ⚠️ **"10 mg/L and 30 µg/L"** | **absolute standard, but see the unit defect below** |
| 3-DG, glucosone | LC-ESI-MS SIM (m/z 235, 251), quinoxalines | external, authentic 3-DG (75 %) and glucosone (≥98 %), 0.1–1 mg/L | 2.5–15 / 8.3–50 µg/L | absolute |
| MGO, GO | SIM m/z 145, 131 | 2-methyl- / quinoxaline, 0.1–1.0 mg/L | as above | absolute |
| **threosone** | SIM m/z 191 | ⚠️ **semi-quantitated on the glucosone curve** — *"the calibration curve of glucosone was used for semi-quantitation of threosone derivatives, since both have same proton-accepting groups"* | — | **`absolute_concentration: false`** |
| free amino acids | Waters Acquity UPLC-TQD, Atlantis HILIC | external, 0.1–5.0 mg/L | 0.2–5 / 0.7–16.7 µg/L | absolute |

⚠️ **PRINTED UNIT DEFECT, §2.4:** *"The LOD and LOQ values for HMF were **10 mg/L and 30 µg/L**"* —
an LOD three orders of magnitude *above* its own LOQ, and above the entire calibration range's
lower bound (1 mg/L). Almost certainly "10 µg/L and 30 µg/L". **Do not ingest either number.**

---

## §3. THE REACTION NETWORK — 15 steps (Fig. 3), and the HMF ODE

Abbreviations verbatim (Fig. 2/3 captions): `SUC` sucrose; `GLU` glucose; `FRU` fructose;
**`FFC` fructofuranosyl cation**; `3-DG`; `G` glucosone; `MGO`; `GO`; `T` threosone; `HMF`;
`P` products; `1,2-ED` 1,2-enediol.

```
 1  SUC → FRU + GLU          9  3-DG → MGO
 2  GLU → 1,2-ED            10  FRU → FFC
 3  1,2-ED → GLU            11  FFC → HMF     ★ terminal, FRUCTOSE LIMB
 4  1,2-ED → FRU            12  G → GO
 5  FRU → 1,2-ED            13  G → T
 6  1,2-ED → 3-DG           14  G → P1
 7  1,2-ED → G              15  3-DG → P2
 8  3-DG → HMF   ★ terminal, 3-DG LIMB
```

**★ THE HMF NODE ODE — Appendix A, verbatim from the typeset PDF:**

> **`d[HMF]/dt = k₁₁·[FFC] + k₈·[3-DG]`**

⚠️⚠️ **HMF HAS NO SINK IN THIS MODEL AT ALL.** Not a first-order decay (as in both Kocadağlı
papers and Göncüoğlu Taş), not a bimolecular amine reaction (as in Şen 2022). **HMF is a terminal,
strictly accumulating species here.** That is defensible at 27–37 °C over 24 weeks *and it is a
falsifiable modelling choice*, but it means **this paper contributes zero information about HMF
degradation** and must never be used to argue the sink is absent at cooking temperature.

Full ODE set, verbatim from Appendix A:
```
d[SUC]/dt     = −k1[SUC]
d[GLU]/dt     =  k1[SUC] + k3[1,2−enediol] − k2[GLU]
d[FRU]/dt     =  k1[SUC] + k4[1,2−enediol] − (k5+k10)[FRU]
d[HMF]/dt     =  k11[FFC] + k8[3−DG]
d[3−DG]/dt    =  k6[1,2−enediol] − (k8+k9+k15)[3−DG]
d[G]/dt       =  k7[1,2−enediol] − (k12+k13+k14)[G]
d[GO]/dt      =  k12[G]
d[T]/dt       =  k13[G]
d[MGO]/dt     =  k9[3−DG]
d[1,2−enediol]/dt = k2[GLU] + k5[FRU] − (k3+k4+k6+k7)[1,2−enediol]
d[FFC]/dt     =  k10[FRU] − k11[FFC]
d[P1]/dt      =  k14[G]
d[P2]/dt      =  k15[3−DG]
```

⚠️ **BOTH `1,2-ED` AND `FFC` ARE UNMEASURED.** The authors say so for the enediol, verbatim
§3.2(i): *"A possible explanation for this observation is that **1,2-enediol intermediate could not
be quantified. Hence, an unquantified compound causes the highest intervals or indetermination of
parameters in a mechanistic model.**"* The same holds for FFC (never analysed; no standard). ⇒
**k₁₀ and k₁₁ are identified only up to a common scale on the FFC pool; k₂–k₇ share a second
unidentified scale on the enediol pool. Nine of the fifteen constants are scale-free-only.**

---

## §4. TABLE 1 — rate constants, ± 95 % HPD

**Transcribed verbatim from p.8.** Caption verbatim: *"Estimated reaction rate constants
(k, week⁻¹ × 10³) with 95 % highest posterior density (HPD) intervals according to the proposed
kinetic model shown in Fig. 3 for sugar degradation reactions in apple juice, orange juice and
peach nectar during storage at different temperatures."* Provenance **`[F]`**.

⚠️ **PRINTED HEADER DEFECT — the column header row lists only FIVE temperature labels
("37 °C | 37 °C | 27 °C | 37 °C | 27 °C") for SIX k/HPD column pairs. The first "27 °C" is missing.**
I recovered the assignment from the paper's own running text, §3.2(i), verbatim: *"The hydrolysis
rates of sucrose to glucose and fructose (k₁) were **0.04, 0.05 and 0.05 week⁻¹ at 27 °C; 0.1, 0.2
and 0.2 week⁻¹ at 37 °C** for apple juice, orange juice and peach nectar, respectively."*
Check `[D]`: 27 °C column values 36.3 / 50.7 / 44.8 ×10⁻³ = 0.036 / 0.051 / 0.045 ✔;
37 °C column values 123.7 / 147.4 / 147.2 ×10⁻³ = 0.12 / 0.15 / 0.15 ✔.
**Confirmed column order: Apple 37, Apple 27, Orange 37, Orange 27, Peach 37, Peach 27.**

| # | step | Apple 37 | HPD | Apple 27 | HPD | Orange 37 | HPD | Orange 27 | HPD | Peach 37 | HPD | Peach 27 | HPD |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SUC → FUR + GLU | 123.7 | ±10.03 | 36.3 | ±3.06 | 147.4 | ±9.68 | 50.7 | ±3.29 | 147.2 | ±13.1 | 44.8 | ±3.79 |
| 2 | GLU → 1,2-ED | 598.4 | **ind\*** | 33.1 | ±22.02 | 58.3 | ±29.84 | 37.5 | ±8.48 | 3690.7 | **ind\*** | 595.5 | ±556.3 |
| 3 | 1,2-ED → GLU | 781.8 | ±111.70 | 66.6 | ±104.20 | 495.4 | ±277.1 | 99.3 | ±29.90 | 2047.3 | ±161.5 | 441.6 | ±531.4 |
| 4 | 1,2-ED → FRU | 544.3 | ±390.00 | 90.5 | ±83.33 | 28.5 | ±30.53 | 39.2 | ±34.08 | 407.8 | ±282 | 65.0 | ±68.81 |
| 5 | FRU → 1,2-ED | 143.1 | ±98.33 | 0.0 | ±0.0 | 0.9 | **fixed** | 17.8 | ±9.34 | 684.2 | ±463.6 | 105.0 | ±66.6 |
| 6 | 1,2-ED → 3-DG | 0.7 | ±0.08 | 2.1 | ±1.22 | 6.7 | ±2.14 | 0.4 | ±0.05 | 0.3 | ±0.037 | 36.4 | ±33.03 |
| 7 | 1,2-ED → G | 2.0 | ±0.27 | 9.6 | ±5.73 | 20.4 | ±6.59 | 1.2 | ±0.14 | 0.0001 | ±0.002 | 0.00002 | ±0.009 |
| **8** | **3-DG → HMF** | **0.7** | **±1.98** | **0.0** | **fixed** | **2.0** | **±0.68** | **0.0** | **fixed** | **11.8** | **±2.79** | **0.0** | **fixed** |
| 9 | 3-DG → MGO | 4.6 | ±0.41 | 3.5 | ±0.43 | 6.0 | ±0.60 | 11.7 | ±2.51 | 0.0 | **fixed** | 0.0 | **fixed** |
| 10 | FRU → FFC | 0.2 | ±0.04 | 0.0 | ±0 | 7.1 | ±2.26 | 0.0 | ±0.0 | 22.6 | ±12.69 | 0.0 | **fixed** |
| **11** | **FFC → HMF** | **5.3** | **ind\*** | **0.0** | **fixed** | **0.0** | **±0.02** | **0.0** | **fixed** | **0.4** | **±0.22** | **0.0** | **fixed** |
| 12 | G → GO | 10.0 | ±1.02 | 5.1 | ±1.36 | 9.2 | ±0.74 | 8.7 | ±1.49 | 23.2 | ±1.57 | 4.8 | ±1.25 |
| 13 | G → T | 0.0 | **fixed** | 0.0 | **fixed** | 0.0 | **fixed** | 0.0 | **fixed** | 22.1 | ±1.19 | 29.4 | ±1.29 |
| 14 | G → P1 | 56.1 | **fixed** | 562.0 | **fixed** | 453.2 | **fixed** | 21.1 | **fixed** | 0.0 | **fixed** | 29.7 | ±17.9 |
| 15 | 3-DG → P2 | 107.6 | **fixed** | 123.3 | **fixed** | 228.2 | **fixed** | 71.7 | **fixed** | 427.2 | **fixed** | 94365.5 | ±58130 |

`*ind` footnote verbatim: *"indeterminate, which means a large uncertainty in the estimated
parameter within 95 % confidence interval."*

### 4.1 ⚠️ HOW MUCH OF THIS TABLE IS ACTUALLY ESTIMATED

Of **90 numeric cells**: **21 are marked `fixed`** (not estimated — set by hand), **2 are `ind`**,
**19 more have HPD ≥ estimate**. **⇒ 42 of 90 cells (47 %) carry no usable information.**
On the two HMF-forming steps specifically, **10 of 12 cells are either `fixed` at zero, `ind`, or
have an HPD wider than the estimate.** The two survivors are orange 37 °C k₈ = 2.0 ± 0.68 and
peach 37 °C k₈ = 11.8 ± 2.79 — **both on the 3-DG limb.**

⚠️ Note also **step 15 peach 27 °C: 94 365.5 ± 58 130**, i.e. 94.4 week⁻¹, which is **220× larger
than the same step at 37 °C (427.2 ×10⁻³)** and ~10⁴× every other constant in the table. That is
the stiff-unidentified-pair signature the repo already refused in Han 2025 (round-3 §A.3.1 (iii)).

---

## §5. TABLE 2 — activation energies — ★ **REFUSE ENTIRELY. NOT ONE VALUE REPRODUCES.**

**Transcribed verbatim from p.9.** Caption verbatim: *"Estimated activation energies (Ea, kJ/mol)
according to the proposed kinetic model shown in Fig. 3 for sugar degradation reactions in apple
juice, orange juice and peach nectar during storage at different temperatures."*

| # | step | Apple Ea | R² | Orange Ea | R² | Peach Ea | R² |
|---|---|---|---|---|---|---|---|
| 1 | SUC → FUR + GLU | 54.16 | 0.96 | 55.51 | 0.98 | 48.23 | 0.97 |
| 2 | GLU → 1,2-ED | 11.10 | 0.03 | 48.19 | 0.99 | 34.95 | 0.43 |
| 3 | 1,2-ED → GLU | −5.94 | 0.01 | 56.29 | 0.90 | 79.16 | 0.86 |
| 4 | 1,2-ED → FRU | 2.41 | 0.00 | −81.09 | 0.96 | −18.61 | 0.51 |
| 5 | FRU → 1,2-ED | 65.09 | 0.53 | −145.77 | 0.97 | 19.99 | 0.18 |
| 6 | 1,2-ED → 3-DG | −20.02 | 0.52 | 112.67 | 0.94 | −141.35 | 0.90 |
| 7 | 1,2-ED → G | 10.12 | 0.07 | 94.62 | 0.88 | 23.43 | 0.72 |
| **8** | **3-DG → HMF** | **−127.02** | **0.51** | **−107.71** | **0.51** | **−77.04** | **0.51** |
| 9 | 3-DG → MGO | 18.43 | 1.00 | −24.21 | 0.90 | – | – |
| 10 | FRU → FFC | −30.71 | 0.03 | 65.16 | 0.21 | −65.77 | 0.51 |
| **11** | **FFC → HMF** | **−90.94** | **0.51** | **−230.31** | **0.51** | **−188.99** | **0.51** |
| 12 | G → GO | −9.82 | 0.72 | −0.15 | 0.02 | −97.29 | 0.71 |
| 13 | G → T | – | – | – | – | 26.54 | 0.79 |
| 14 | G → P1 | −36.53 | 0.46 | 15.13 | 0.06 | 54.27 | 0.34 |
| 15 | 3-DG → P2 | −26.52 | 0.97 | −41.38 | 0.56 | 121.55 | 0.41 |

The authors' own comment, verbatim §3.2: *"the activation energies for each reaction step was found
to be in the range of **−230 and 122 kJ/mol** … The calculated **negative activation energy values
might indicate that no energy barriers were presented in these reaction steps** due to the
accumulation of intermediate compounds (Vyazovkin, 2016)."*

### 5.1 ★ WAVE K5a CROSS-CHECK — I recomputed every Ea from the paper's own Table 1 `[D]`

**Method.** Kinetics exist at exactly **two** temperatures (27 and 37 °C; 4 °C was excluded, §2).
A two-point Arrhenius is therefore fully determined:
`Ea = −R·(ln k₃₇ − ln k₂₇) / (1/T₃₇ − 1/T₂₇)`, R = 8.3145×10⁻³ kJ mol⁻¹ K⁻¹.

**Result: NOT ONE of the 43 published Ea values reproduces. Not one. Zero out of forty-three.**
Full comparison:

| step | juice | k₂₇ | k₃₇ | **Ea recomputed** | Ea published |
|---|---|---|---|---|---|
| SUC→FRU+GLU | Apple | 36.3 | 123.7 | **94.9** | 54.16 |
| SUC→FRU+GLU | Orange | 50.7 | 147.4 | **82.6** | 55.51 |
| SUC→FRU+GLU | Peach | 44.8 | 147.2 | **92.1** | 48.23 |
| GLU→1,2-ED | Apple | 33.1 | 598.4 | **224.1** | 11.10 |
| GLU→1,2-ED | Orange | 37.5 | 58.3 | **34.2** | 48.19 |
| GLU→1,2-ED | Peach | 595.5 | 3690.7 | **141.2** | 34.95 |
| 1,2-ED→GLU | Apple | 66.6 | 781.8 | **+190.6** | **−5.94** ⚠️ sign flip |
| 1,2-ED→GLU | Orange | 99.3 | 495.4 | **124.4** | 56.29 |
| 1,2-ED→GLU | Peach | 441.6 | 2047.3 | **118.7** | 79.16 |
| 1,2-ED→FRU | Apple | 90.5 | 544.3 | **+138.9** | +2.41 |
| 1,2-ED→FRU | Orange | 39.2 | 28.5 | −24.7 | −81.09 |
| 1,2-ED→FRU | Peach | 65.0 | 407.8 | **+142.1** | **−18.61** ⚠️ sign flip |
| FRU→1,2-ED | Apple | **0.0** | 143.1 | **NOT COMPUTABLE** | 65.09 |
| FRU→1,2-ED | Orange | 17.8 | 0.9 | −231.0 | −145.77 |
| FRU→1,2-ED | Peach | 105.0 | 684.2 | **+145.1** | +19.99 |
| 1,2-ED→3-DG | Apple | 2.1 | 0.7 | −85.0 | −20.02 |
| 1,2-ED→3-DG | Orange | 0.4 | 6.7 | **218.1** | 112.67 |
| 1,2-ED→3-DG | Peach | 36.4 | 0.3 | −371.4 | −141.35 |
| 1,2-ED→G | Apple | 9.6 | 2.0 | **−121.4** | **+10.12** ⚠️ sign flip |
| 1,2-ED→G | Orange | 1.2 | 20.4 | **219.3** | 94.62 |
| 1,2-ED→G | Peach | 0.00002 | 0.0001 | 124.6 | 23.43 |
| **3-DG→HMF** | **Apple** | **0.0 (fixed)** | 0.7 | ★ **NOT COMPUTABLE** | **−127.02** |
| **3-DG→HMF** | **Orange** | **0.0 (fixed)** | 2.0 | ★ **NOT COMPUTABLE** | **−107.71** |
| **3-DG→HMF** | **Peach** | **0.0 (fixed)** | 11.8 | ★ **NOT COMPUTABLE** | **−77.04** |
| 3-DG→MGO | Apple | 3.5 | 4.6 | 21.2 | 18.43 |
| 3-DG→MGO | Orange | 11.7 | 6.0 | −51.7 | −24.21 |
| **FRU→FFC** | Apple / Orange / Peach | **0.0** in all three | — | ★ **NOT COMPUTABLE** | −30.71 / 65.16 / −65.77 |
| **FFC→HMF** | **Apple** | **0.0 (fixed)** | 5.3 | ★ **NOT COMPUTABLE** | **−90.94** |
| **FFC→HMF** | **Orange** | **0.0 (fixed)** | **0.0** | ★ **NOT COMPUTABLE** | **−230.31** |
| **FFC→HMF** | **Peach** | **0.0 (fixed)** | 0.4 | ★ **NOT COMPUTABLE** | **−188.99** |
| G→GO | Apple | 5.1 | 10.0 | **+52.1** | **−9.82** ⚠️ sign flip |
| G→GO | Orange | 8.7 | 9.2 | +4.3 | −0.15 |
| G→GO | Peach | 4.8 | 23.2 | **+121.9** | **−97.29** ⚠️ sign flip |
| G→T | Peach | 29.4 | 22.1 | **−22.1** | **+26.54** ⚠️ sign flip |
| G→P1 | Apple | 562.0 | 56.1 | −178.4 | −36.53 |
| G→P1 | Orange | 21.1 | 453.2 | **237.4** | 15.13 |
| G→P1 | Peach | 29.7 | **0.0** | NOT COMPUTABLE | 54.27 |
| 3-DG→P2 | Apple | 123.3 | 107.6 | −10.5 | −26.52 |
| 3-DG→P2 | Orange | 71.7 | 228.2 | **+89.6** | **−41.38** ⚠️ sign flip |
| 3-DG→P2 | Peach | 94365.5 | 427.2 | −417.8 | +121.55 |

### 5.2 Four independent, individually disqualifying defects

**(i) ★ THE SIX HMF ACTIVATION ENERGIES ARE MATHEMATICALLY UNDERIVABLE FROM THE PAPER'S OWN DATA.**
Both HMF-forming steps have `k = 0.0`, **explicitly marked `fixed` by the authors**, at 27 °C in
**all three juices**. `ln(0) = −∞`. **No finite activation energy of any value can be obtained.**
Yet Table 2 prints six specific negative numbers to two decimal places. **This is not a
misestimate; the numbers have no derivation.** They belong in the same register as the three
fabrication-class citation defects of `k3` §0.3 — with the difference that here it is the *source
paper*, not the repo, that produced them.

**(ii) THE R² COLUMN IS IMPOSSIBLE.** With exactly two temperatures, a straight line through two
points has **R² = 1.000 identically**, always. Table 2 prints R² values of 0.00, 0.01, 0.02, 0.03,
0.06, 0.07 … 1.00. **Whatever generated this column, it was not a regression on the Table 1 k
values.** Note further that **`R² = 0.51` appears on six different rows** with six different data
(steps 8×3, 11×3 and 10-peach) — the exact rows whose k is fixed at zero. A repeated constant on
precisely the underivable rows is a placeholder, not a statistic.

**(iii) SEVEN OUTRIGHT SIGN FLIPS** between the recomputation and the publication, on steps that
*are* computable (1,2-ED→GLU apple, 1,2-ED→FRU peach, 1,2-ED→G apple, G→GO apple and peach, G→T
peach, 3-DG→P2 orange). A sign flip is not a rounding difference.

**(iv) EVEN THE BEST-BEHAVED STEP DISAGREES BY ~1.8×.** Sucrose hydrolysis (step 1) is the one
process in the paper with clean, monotone, tight-HPD k in all six cells. Recomputed Ea:
**94.9 / 82.6 / 92.1**. Published: **54.16 / 55.51 / 48.23**. The paper's own text quotes the
published range approvingly — *"especially hydrolysis of sucrose into glucose and fructose was
fairly temperature dependent in juices (Ea; 48–56 kJ/mol)"*. **The correct value from their own
table is 83–95 kJ/mol.**

**⇒ VERDICT: Gürsul Aktağ & Gökmen 2020 Table 2 is REFUSED IN ITS ENTIRETY. Record the reason so a
later wave does not re-ingest it.** Table 1 is unaffected by this and remains usable under §4.1.

---

## §6. ABSOLUTE MEASUREMENTS `[M]` — the part of the paper that is sound

| quantity | apple | orange | peach nectar | anchor |
|---|---|---|---|---|
| **max HMF (37 °C, 24 wk)** | **16.2 ± 0.7 mg/L** | **3.8 ± 0.2 mg/L** | **12.2 ± 0.5 mg/L** | §3.1 |
| HMF at 27 °C | ★ **no accumulation at all** | " | " | §3.1: *"formation of HMF followed a typical kinetic pattern in juices stored at 37 °C, while **there was no accumulation of HMF at 27 °C**"* |
| max glucosone (37 °C) | 3.5 mmol/L @ 16 wk | 1.7 mmol/L @ 14 wk | (lower) | §3.1 |
| max 3-DG (37 °C) | — | — | 0.4 mmol/L @ 14 wk | §3.1 |
| max GO (37 °C, 24 wk) | 0.4 mmol/L | 0.2 mmol/L | 0.08 mmol/L | §3.1 |
| **dominant α-dicarbonyl** | **glucosone** | **glucosone** | **3-DG** | §3.1 |
| MGO | detected | detected | **not detected** (k₉ fixed 0) | §3.1 |
| threosone | not detected | not detected | **detected only here** | §3.1 |
| **1-DG** | ★ **not detected in any sample** | " | " | §3.2(iii): *"1-DG is formed by 2,3-enolization of fructose under alkaline conditions. Expectedly, **1-DG was not detected in the acidic samples**"* |
| sucrose loss, 24 wk, 27 °C | 50 % | 62 % | 54 % | §3.1 |
| sucrose loss, 24 wk, 37 °C | 93 % | 89 % | 96 % | §3.1 |
| initial fructose ratio | **2.3× orange, 9.3× peach** | — | — | §3.2(ii) |
| regulatory ceiling | **AIJN max 10 mg/L HMF for fruit juices** | | | §3.1 `[C]` |

★ **`max HMF = 16.2 / 3.8 / 12.2 mg/L` at 37 °C after 24 weeks, with zero accumulation at 27 °C, in
three named matrices at pH 3.4 with a measured amine-off control — this is a clean, absolutely
quantified, authentic-standard hold-out dataset.** It is the *only* food-pH, storage-temperature
absolute HMF anchor in the whole K5a cluster.

---

## §7. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL`

| # | constraint | anchor |
|---|---|---|
| S1 | **At pH 3.4 and 4–37 °C the Maillard reaction does not run; sugar dehydration alone accounts for the α-dicarbonyls and HMF** (amino acids ANOVA-flat) | §3.1 |
| S2 | **pH does not drift over 24 weeks** — unlike the initial-pH problem of `k3` §C.11 (Zhou 2023) | §3.1, Table S3 |
| S3 | **Acid sucrose hydrolysis is the entry point**, not thermal cleavage | §3.2(i) |
| S4 | **1-DG is absent at acidic pH** (2,3-enolisation needs alkali) | §3.2(iii) — ★ a clean pH sign constraint for `k3` §B.2 |
| S5 | **1,2-enolisation is a required step in acidic aqueous juice** — omitting it breaks the sucrose-hydrolysis estimates | §3.2 |
| S6 | ★ **Same lab, opposite verdict on the same step:** *"Contrarily, Kocadağlı and Gökmen (2016) indicated that 1,2-enediol formation in interconversion of glucose-fructose was unnecessary … This suggests that **the physical form of the reactants in food systems is a strong determinant in whether the enolisation step was important or not.**"* | §3.2(i), verbatim |
| S7 | **Within the fructose limb, FRU→FFC is fast and FFC→HMF is rate-determining** (orange, peach) | §3.2(ii) — ★ **the exact opposite of Göncüoğlu Taş 2017's hazelnut finding for the 3-DG limb, where the terminal step is fast** |
| S8 | **3-DG → MGO out-runs 3-DG → HMF** | §3.2(iii): *"degradation of 3-DG to MGO (k₉) was kinetically more important than its degradation to HMF (k₈)"* — check `[D]`: apple 4.6 vs 0.7 (6.6×), orange 6.0 vs 2.0 (3.0×); ⚠️ **fails in peach**, where k₉ = 0 (fixed) and k₈ = 11.8 |
| S9 | Under acid, **dehydration to furaldehydes is favoured over isomerisation** | §3.2 `[C]` (Simpson 2012) |
| S10 | **GO and threosone come from retro-aldolisation of glucosone**, not from 3-DG | §3.2(iii) |
| S11 | Glucosone dominates in apple/orange; 3-DG dominates in peach nectar — **the dominant dicarbonyl is matrix-set, at identical pH** | §3.1 ★ another branch-fraction-is-not-a-constant datum for `k3` §0 |

---

## §8. VERIFIED NEGATIVES `[NEG]`

- **No HMF sink of any kind in the model.** (§3)
- **No 3,4-DG.** The 3-DG→HMF step is written as a single lumped dehydration; the intermediate that
  both Kocadağlı papers and the hazelnut paper measure is **absent from this network**.
  ⇒ **This paper's `k₈` is NOT commensurable with their `k₅`/`k₁₃`/`k₁₉`; it lumps two steps.**
- **No 1-DG** (not detected), **no diacetyl**, **no furfural**, **no dimethylglyoxal**.
- **No melanoidin / browning measurement.**
- **No amino-acid consumption** (measured to be flat — that is the point).
- **No water activity.**
- Supplementary Tables S1–S4 and Figures S1–S4 are **not in the PDF on disk**; the comprehensive
  (pre-simplification) model's rate constants live in Table S4 and are **not held**.

---

## §9. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| **Table 2 — all 43 activation energies** | `[F]` | ★ **REFUSE ENTIRELY.** Six are mathematically underivable (k = 0 fixed at 27 °C); the R² column is impossible for a 2-point fit; seven sign flips; the best-behaved step is off by 1.8×. **Record the refusal.** |
| Table 1, the 48 estimated cells | `[F]` | **USE-Q** — `week⁻¹×10³`, per-temperature, per-juice; **drop all 21 `fixed`, both `ind`, and the 19 HPD ≥ estimate cells** |
| **k₈ (3-DG→HMF), orange 37 °C = 2.0 ± 0.68; peach 37 °C = 11.8 ± 2.79** | `[F]` | **USE-Q — the only two well-determined HMF-forming constants in the paper, and both are on the 3-DG limb.** ⚠️ lumps 3-DG→3,4-DG→HMF into one step |
| **k₁₀, k₁₁ (fructose limb)** | `[F]` | **RATIO-ONLY** — `[FFC]` unmeasured; and 10 of 12 cells are `fixed`/`ind`/HPD-dominated |
| **the abstract's "significantly higher (p < 0.05)" fructose-limb claim** | — | ★ **REFUSE. Not supported by the described statistics, contradicted by Table 1 in 2 of 3 juices, and assembled from two different reaction steps (§1).** |
| max HMF 16.2 / 3.8 / 12.2 mg/L @ 37 °C, 24 wk; **no accumulation at 27 °C** | `[M]` | ★ **USE — proposed HOLD-OUT** (see synthesis §7) |
| glucosone / 3-DG / GO maxima, sucrose losses, initial fructose ratios | `[M]` | **USE** |
| pH 3.4, stable over 24 weeks | `[M]` | **USE** |
| HMF LOD/LOQ "10 mg/L and 30 µg/L" | `[M]` | **REFUSE — internally impossible** |
| threosone concentrations | `[M]` | **`absolute_concentration: false`** |
| S1–S11 | `[M]`/`[F]` | **STRUCTURAL** — S1, S4, S6 and S11 are the highest-value |
