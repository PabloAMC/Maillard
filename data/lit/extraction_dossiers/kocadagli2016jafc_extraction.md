# Kocadağlı & Gökmen 2016 (JAFC) — glucose caramelization ± NaCl

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

---

## §0. IDENTITY — ⚠️ FILENAME/DOI MAPPING IS THE **OPPOSITE** OF THE NAÏVE GUESS

| item | value | how verified |
|---|---|---|
| **file on disk** | `data/articles/Kocada2016.pdf` (810,015 bytes) | `ls` |
| **DOI** | **`10.1021/acs.jafc.6b01862`** | printed on the ACS cover page of the PDF text layer, line 8 |
| title | *"Effect of sodium chloride on α-dicarbonyl compounds and 5-hydroxymethyl-2-furfural formations from glucose under caramelization conditions – A multiresponse kinetic modelling approach"* | PDF p.1 (cover) and p.1 of the manuscript body |
| authors | Tolgahan Kocadağlı, Vural Gökmen | manuscript p.1 |
| journal | *J. Agric. Food Chem.*, **Just Accepted Manuscript**, Publication Date (Web): **30 Jul 2016** | ACS cover page |
| corresponding | Prof. Dr. Vural Gökmen, `vgokmen@hacettepe.edu.tr` | manuscript lines 8–10 |
| pages of PDF | 38 (ACS "Page N of 38" running head) | throughout |

**⚠️ WRONG-FILE WARNING, REPORT UPSTREAM.** The task brief and the round-3 dossier both implied
the longer filename would be the JAFC paper. **It is not.** The mapping is:

| filename | is actually | DOI |
|---|---|---|
| **`Kocada2016.pdf`** | ★ **the JAFC glucose + NaCl caramelization paper — THE HMF PICK** | `10.1021/acs.jafc.6b01862` |
| **`Kocadagli2016.pdf`** | the *Food Chemistry* glucose/wheat-flour paper | `10.1016/j.foodchem.2016.05.150` |

Anything that ingested by filename stem alone has them swapped.

**⚠️ SECOND WARNING — THIS IS A "JUST ACCEPTED" MANUSCRIPT, NOT THE VERSION OF RECORD.**
Verbatim from the ACS cover page (PDF p.1):

> "'Just Accepted' manuscripts have been peer-reviewed and accepted for publication. They are
> posted online prior to technical editing, formatting for publication and author proofing …
> **should not be considered the official version of record** … Note that technical editing may
> introduce minor changes to the manuscript text and/or graphics which could affect content"

The final version is *J. Agric. Food Chem.* **64**, 6333–6342 (the citation carried in the
round-3 dossier). **Every number below is from the Just Accepted PDF and may differ from the
version of record.** Do not cite page numbers from the printed article against this file.

---

## §1. SYSTEM AND CONDITIONS — all verbatim, `[M]`

| item | value | anchor |
|---|---|---|
| system A | **glucose 0.1 M in water**, 0.5 mL per glass tube, frozen −80 °C, freeze-dried | lines 116–125 |
| system B | **glucose–NaCl (0.1 M each)** in water, same handling | lines 116–125 |
| amine | **NONE. There is no amino acid and no protein anywhere in this system.** | §Preparation, and abstract |
| physical state | freeze-dried solid → *"two different types of amorphous states, which was not characterized in this study"* | lines 122–125, verbatim |
| temperatures | **160, 180 and 200 °C** | line 126 |
| times | **up to 30 min**, in **duplicates** | lines 126–127 |
| heating | oil bath (Memmert), PTFE-sealed screw caps (**closed system, headspace present**) | lines 125–128 |
| extraction | 2.5 mL water, vortex/shake 1 min | lines 128–129 |
| units of the model | **µmol** (absolute amount per tube, not a concentration) | line 194: *"The amount of reactants and products were expressed as μmol"* |
| initial glucose, system A | **47.1 ± 0.66 µmol** | line 215 |
| initial glucose, system B | **56.4 ± 0.77 µmol** | line 216 |

⚠️ **The two systems are NOT iso-loaded: 56.4 vs 47.1 µmol = a 19.7 % higher glucose charge in the
NaCl system**, despite the authors' stated intent (*"all tubes contained same amount of glucose by
pipetting rather than weighting sugar crystals"*, lines 117–118). Per `k3` §0, a 2× precursor
change moves a branch fraction by 3×; a 1.2× change is smaller but is **not** nothing, and it is
in the same direction as the reported NaCl effect. The mole-conversion numbers (§2.3) are
normalised and are therefore the safer comparison; the raw µmol are not.

### 1.1 Analytical methods — quantification quality varies by analyte `[M]`

| analyte | method | calibration | verdict |
|---|---|---|---|
| glucose, fructose | HPLC-RID, Shodex Sugar SH-1011, 5 mM H₂SO₄, 50 °C | external, **0.005–1 %** | absolute |
| **HMF** | Shimadzu UFLC, Atlantis dC18 250×4.6, 10 mM formic acid/MeCN 90:10, DAD **285 nm** | external, **0.1–10 mg/L**, authentic HMF (98 %, Acros) | **absolute, authentic standard** |
| 3-DG | LC-ESI-MS SIM, quinoxaline derivative, m/z 235 | external, **0.1–5 mg/L**, authentic 3-DG (75 %) | absolute |
| glyoxal / methylglyoxal / diacetyl | SIM m/z 131 / 145 / 159 | external, quinoxaline / 2-methyl- / 2,3-dimethylquinoxaline, 0.1–2 mg/L | absolute |
| **glucosone, 1-DG, 3,4-DG** | SIM m/z 251 / 235 / 217 | ⚠️ **SEMI-QUANTITATED against the 3-DG curve** — *"this calibration curve was used for semi-quantitation of glucosone, 1–deoxyglucosone and 3,4-dideoxyglucosone quinoxaline derivatives since both have same proton accepting groups"* (lines 163–165) | **`absolute_concentration: false`** |

⚠️ **3,4-DG is semi-quantitated.** 3,4-DG is the *immediate precursor of HMF on the 3-DG limb*
(step 5). Its concentration scale is therefore borrowed from 3-DG's response factor, which means
**k₅ (3,4-DG→HMF) carries an unknown multiplicative response-factor error.** This is the same
class flagged in `k3` §B.9 for semi-quant HS-SPME rows. Note that 1-DG and 3-DG share the same
SIM ion (m/z 235) and are separated only chromatographically.

---

## §2. THE REACTION NETWORK — 18 elementary steps, Figure 2, transcribed from the ASCII figure and Appendix A

**Edge list (Figure 2, PDF p.34).** Abbreviations verbatim from the paper (lines 483–486):
`Glc` glucose; `Fru` fructose; `1–DG` 1–deoxyglucosone; `3–DG` 3–deoxyglucosone;
`3,4–DG` 3,4–dideoxyglucosone; `G` glucosone; `GO` glyoxal; `MG` methylglyoxal; `DA` diacetyl;
`HMF` 5–hydroxymethyl–2–furfural; `Int` intermediate; `P` products.

```
 1  Glc  → Fru            (isomerisation, forward)
 2  Fru  → Glc            (isomerisation, reverse)
 3  Glc  → 3-DG           (dehydration, −H2O)
 4  3-DG → 3,4-DG         (dehydration, −H2O)
 5  3,4-DG → HMF          ★ terminal step, 3-DG LIMB
 6  Fru  → Int            (dehydration over the intact fructofuranose ring, −H2O)
 7  Int  → HMF            ★ terminal step, FRUCTOSE LIMB
 8  Fru  → 1-DG           (2,3-enolisation)
 9  Glc  → G              (oxidation)
10  G    → GO
11  3-DG → MGO
12  1-DG → DA
13  3-DG → P1             (sink)
14  1-DG → P2             (sink)
15  GO   → P3             (sink)
16  MGO  → P4             (sink)
17  DA   → P5             (sink)
18  HMF  → P6             ★ THE HMF SINK, first order in HMF
```

**★ THE HMF NODE ODE, from Appendix A (PDF p.37).** The text layer strips the subscripts; the
equation is recovered unambiguously by cross-reading Appendix A against the Figure 2 edge list and
Table 1 step numbering:

> **`d[HMF]/dt = k₅·[3,4-DG] + k₇·[Int] − k₁₈·[HMF]`**

i.e. **two parallel first-order sources and exactly one first-order sink.** ⚠️ Note the
architectural difference from Han 2025 (`research_round3_channels.md` §A.3.1), whose sink is
`HMF → melanoidin`: here the sink is to an unnamed lumped `P6` and **melanoidins are not modelled
at all** — the paper's own conclusion asks for them to be added.

⚠️ **`Int` IS AN UNMEASURED, UNQUANTIFIED SPECIES.** Verbatim, lines 255–260:

> "Dehydration of fructose to form HMF via cyclic intermediates was reduced to two steps,
> comprising an **undetermined intermediate (Int)** for simplification. These undetermined
> intermediates **could be** the enol form of 2,5-anhydro-D-mannose … and subsequently this enol
> intermediate dehydrate from C3 to form a 2,3-dihydrofuran and further dehydration produce HMF"

**Consequence, and it is the single most important caveat in this dossier:** because `[Int]` is
never measured, only the *product* `k₇·[Int]` is constrained by the data. `k₆` and `k₇` are
therefore identifiable **only up to a common scale factor** on the `Int` pool. **No absolute
transfer of k₆ or k₇ into another model is licensed.** Ratios of k₆ (or k₇) *between the two
systems in this same paper* are safe, because the scale cancels. See §5.

### 2.1 Model discrimination — the one test that is genuinely decisive `[M]`

Verbatim, lines 424–432:

> "A kinetic model constructed by **omitting the formation of HMF from fructose** through an
> undetermined intermediate **did not fit to the experimentally observed data of either or both
> 3-deoxyglucosone, 3,4-dideoxyglucosone and HMF by no means** (Figure S5, see Supporting
> Information). In this test, other responses not much affected (not shown). **This indicates that
> HMF formation through only dehydration of 3-deoxyglucosone does not correspond to amounts of HMF
> observed quantitatively.**"

⚠️ Figure S5 is in the Supporting Information, which is **NOT in the PDF on disk** (the PDF ends at
the TOC graphic, p.38). The claim is quotable; the supporting figure is not held.

Other discrimination outcomes, all verbatim `[M]`:
- *"The rate constant of 3-deoxyglucosone formation **from fructose** was estimated to be **zero in
  every case**"* (lines 344–345) → **the Fru→3-DG edge is refused by the data.**
- *"If degradation of glucosone to other undetermined products included to model, its rate constant
  was estimated to be **zero**"* (lines 338–340) → no G→P sink.
- *"The model tended to estimate the rate constant of methylglyoxal formation **from
  1-deoxyglucosone** zero in most cases"* (lines 390–391) → MGO comes from 3-DG here.
- *"no mannose formation was observed"* (line 280) → epimerisation edge refused.
- *"the enediol intermediate was not included as the interconversion was obviously fast"*
  (lines 280–282) → **1,2-ED is NOT in this model** (contrast: it IS in the wheat-flour paper and in
  Gürsul Aktağ 2020 — see the synthesis).

---

## §3. TABLE 1 — reaction rate constants at each temperature, ± 95 % HPD

**Transcribed verbatim from PDF p.31.** Caption verbatim: *"Reaction rate constants with 95 %
highest posterior density (HPD) intervals at different temperatures according to the proposed
kinetic model (Figure 2) for caramelization of glucose and glucose-NaCl mixture."*

**Units: `min⁻¹ × 10³` for every row** (i.e. divide the printed number by 1000 to get min⁻¹).
All steps are first order; there is no bimolecular step in this paper (no amine).
`*ind` = *"indeterminate, which means a large uncertainty in the estimated parameter within 95 %
confidence interval."* (footnote, PDF p.31). Provenance for the whole table: **`[F]`**.

### GLUCOSE (no NaCl)

| # | step | k(160) | HPD | k(180) | HPD | k(200) | HPD |
|---|---|---|---|---|---|---|---|
| 1 | Glc→Fru | 237 | ±123 | 1804 | ±81 | 3845 | ±305 |
| 2 | Fru→Glc | 1284 | ±737 | 10409 | **ind\*** | 17657 | **ind\*** |
| 3 | Glc→3-DG | 0.91 | ±0.19 | 4.14 | ±1.71 | 3.60 | ±1.26 |
| 4 | 3-DG→3,4-DG | 23.1 | ±4.03 | 30.5 | ±4.71 | 49.3 | ±10.1 |
| **5** | **3,4-DG→HMF** | **160** | **±35.0** | **110** | **±28.2** | **137** | **±46.1** |
| 6 | Fru→Int | 100 | ±8.6 | 344 | ±26.0 | 1058 | ±96.6 |
| **7** | **Int→HMF** | **0.31** | **±0.07** | **1.87** | **±0.15** | **9.31** | **±1.74** |
| 8 | Fru→1-DG | 0.61 | ±0.15 | 2.47 | ±0.73 | 5.89 | ±1.93 |
| 9 | Glc→G | 0.023 | ±0.002 | 0.054 | ±0.006 | 0.294 | ±0.017 |
| 10 | G→GO | 361 | ±34.9 | 594 | ±80.8 | 2129 | ±149 |
| 11 | 3-DG→MGO | 96.0 | ±20.8 | 338 | ±29.4 | 863 | ±98.1 |
| 12 | 1-DG→DA | 2.71 | ±0.58 | 14.3 | ±1.79 | 68.8 | ±4.99 |
| 13 | 3-DG→P1 | 555 | ±153 | 2241 | ±1117 | 841 | ±678 |
| 14 | 1-DG→P2 | 347 | ±94.5 | 925 | ±293 | 2035 | ±719 |
| 15 | GO→P3 | 66.1 | ±11.1 | 6.49 | ±16.9 | 33.1 | ±10.1 |
| 16 | MGO→P4 | 23.7 | ±19.9 | 85.0 | ±13.4 | 65.6 | ±13.2 |
| 17 | DA→P5 | 5.53 | ±18.1 | 28.8 | ±14.4 | **0** | *(blank)* |
| **18** | **HMF→P6** | **20.6** | **±11.9** | **36.7** | **±7.76** | **203** | **±56.4** |

### GLUCOSE–NaCl

| # | step | k(160) | HPD | k(180) | HPD | k(200) | HPD |
|---|---|---|---|---|---|---|---|
| 1 | Glc→Fru | 212 | ±79 | 4712 | ±599 | 9489 | ±1388 |
| 2 | Fru→Glc | 1000 | ±543 | 24962 | **ind\*** | 52998 | **ind\*** |
| 3 | Glc→3-DG | 0.39 | ±0.07 | 0.99 | ±0.31 | 4.06 | ±2.34 |
| 4 | 3-DG→3,4-DG | 10.6 | ±2.26 | 43.3 | ±9.58 | 101 | ±28.0 |
| **5** | **3,4-DG→HMF** | **46.0** | **±29.5** | **163** | **±57.0** | **418** | **±120** |
| 6 | Fru→Int | 391 | ±60.5 | 1335 | ±184 | 4297 | ±622 |
| **7** | **Int→HMF** | **1.15** | **±0.15** | **10.0** | **±2.62** | **41.1** | **±13.1** |
| 8 | Fru→1-DG | 0.50 | ±0.13 | 1.34 | ±0.19 | 1.96 | ±0.52 |
| 9 | Glc→G | 0.020 | ±0.002 | 0.054 | ±0.009 | 0.27 | ±0.03 |
| 10 | G→GO | 770 | ±159 | 787 | ±171 | 3129 | ±611 |
| 11 | 3-DG→MGO | 167 | ±24.0 | 257 | ±65.3 | 890 | ±217 |
| 12 | 1-DG→DA | 10.7 | ±1.45 | 11.6 | ±0.81 | 88.4 | ±17.4 |
| 13 | 3-DG→P1 | 202 | ±105 | 169 | ±222 | 827 | ±1157 |
| 14 | 1-DG→P2 | 398 | ±126 | 516 | ±76.0 | 445 | ±99.6 |
| 15 | GO→P3 | 83.2 | ±20.7 | 35.1 | ±30.0 | 5.69 | ±27.4 |
| 16 | MGO→P4 | 91.6 | ±24.4 | 84.2 | ±49.3 | 35.6 | ±44.7 |
| 17 | DA→P5 | 2.15 | ±14.7 | **0** | *(blank)* | **0** | *(blank)* |
| **18** | **HMF→P6** | **0.0** | **±0.0** | **263** | **±86.0** | **955** | **±313** |

### 3.1 Rows where the HPD equals or exceeds the estimate (⇒ interval spans zero ⇒ **REFUSE**)

Glucose: step 15 @180 (6.49 ± 16.9), step 16 @160 (23.7 ± 19.9), step 17 @160 (5.53 ± 18.1),
step 13 @200 (841 ± 678, 81 %), step 18 @160 (20.6 ± 11.9, 58 %).
Glucose-NaCl: step 13 @180 (169 ± 222), step 13 @200 (827 ± 1157), step 15 @180 (35.1 ± 30.0) and
@200 (5.69 ± 27.4), step 16 @200 (35.6 ± 44.7), step 17 @160 (2.15 ± 14.7),
**step 5 @160 (46.0 ± 29.5, 64 %)**.
Plus the four `ind` cells on step 2. **13 of 102 numeric cells (12.7 %) are unusable as printed.**

---

## §4. TABLE 2 — reparameterised Arrhenius, ± 95 % HPD

**Transcribed verbatim from PDF p.32.** Reparameterisation, verbatim from lines 200–208:

> "Temperature dependence of the rate constants were evaluated by means of activation energies Ea
> (kJ/mol) defined by Arrhenius equation, which is **reparameterized to the average base
> temperature studied (Tb = 180 °C)** for statistical reasons … where **kb is reparameterized
> pre-exponential factor (equals to the rate constant at Tb)**, R is the universal gas constant
> (8.3145×10⁻³ kJ·mol⁻¹·K⁻¹) … **the data for 160, 180 and 200 °C was simultaneously fitted during
> parameter estimation.**"

⇒ `k(T) = k_b · exp[ −(Ea/R)·(1/T − 1/T_b) ]`, **T_b = 453.15 K**, `k_b` in `min⁻¹ × 10³`.
Provenance **`[F]`**. Footnote `δ` verbatim: *"Zero activation energy (Ea) indicates that the
reaction rate constant (k) of the elementary step **does not follow Arrhenius equation** and the Ea
was **fixed to zero** during parameter estimation."*

| # | step | GLUCOSE k_b | HPD | GLUCOSE Ea | HPD | NaCl k_b | HPD | NaCl Ea | HPD |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Glc→Fru | 2039 | ±83 | 151.5 | ±99.8 | 5279 | ±484 | 263.6 | ±23.1 |
| 2 | Fru→Glc | 10942 | **ind\*** | 146.4 | ±104.2 | 27705 | **ind\*** | 280.8 | ±29.8 |
| 3 | Glc→3-DG | 4.19 | ±2.44 | 107.2 | ±52.7 | 1.11 | ±0.20 | 85.2 | ±11.2 |
| 4 | 3-DG→3,4-DG | 30.5 | ±3.39 | 36.9 | ±6.3 | 36.6 | ±4.18 | 117.7 | ±11.1 |
| **5** | **3,4-DG→HMF** | **119** | **±19.8** | **0** ᵟ | — | **103** | **±26.2** | **169.8** | **±22.0** |
| 6 | Fru→Int | 330 | ±22.8 | 100.4 | ±6.6 | 1402 | ±122 | 94.9 | ±8.1 |
| **7** | **Int→HMF** | **1.84** | **±0.70** | **151.4** | **±34.3** | **8.79** | **±1.67** | **149.8** | **±19.0** |
| 8 | Fru→1-DG | 2.11 | ±0.40 | 99.3 | ±21.8 | 1.36 | ±0.23 | 93.9 | ±18.8 |
| 9 | Glc→G | 0.069 | ±0.005 | 125.9 | ±4.9 | 0.059 | ±0.007 | 131.8 | ±9.2 |
| 10 | G→GO | 737 | ±58.9 | 93.8 | ±6.4 | 1033 | ±170 | 95.2 | ±15.1 |
| 11 | 3-DG→MGO | 304 | ±33.5 | 84.8 | ±6.9 | 309 | ±52.8 | 95.3 | ±13.9 |
| 12 | 1-DG→DA | 12.2 | ±1.12 | 150.8 | ±8.8 | 16.1 | ±3.00 | 139.9 | ±18.3 |
| 13 | 3-DG→P1 | 2239 | ±1504 | 90.9 | ±59.5 | 231 | ±131 | **0** ᵟ | — |
| 14 | 1-DG→P2 | 873 | ±178 | 77.9 | ±23.9 | 574 | ±107 | 55.3 | ±21.1 |
| 15 | GO→P3 | 32.6 | ±8.83 | **0** ᵟ | — | 28.5 | ±21.0 | **0** ᵟ | — |
| 16 | MGO→P4 | 55.4 | ±13.2 | **0** ᵟ | — | 73.4 | ±30.1 | **0** ᵟ | — |
| 17 | DA→P5 | **0** | — | *(blank)* | — | **0** | — | *(blank)* | — |
| **18** | **HMF→P6** | **36.9** | **±64.2** | **152.8** | **±154.4** | **227** | **±63.0** | **137.8** | **±24.6** |

---

## §5. ★ WAVE K5a CROSS-CHECK — I refitted a 3-point Arrhenius to the paper's OWN Table 1

**Method (`[D]` — derived by this analyst, never printed by the authors).** For each step I fitted
`ln k` vs `1/T` by ordinary least squares over the three Table-1 temperatures, and evaluated
`k_b` at T_b = 453.15 K, to compare with the published Table 2. **These two numbers are not
required to agree** — Table 2 comes from a *simultaneous global fit of the Arrhenius-substituted
ODE system to all three temperatures*, not from the per-temperature k. But **the size of the
disagreement, and the R² of the naïve refit, are a direct identifiability diagnostic**: a step
whose per-temperature k already trace a clean Arrhenius line is well determined; a step where the
two disagree by ≥1.5× is one where the global fit is being driven by the other responses, not by
the compound's own curve.

### GLUCOSE (no NaCl)

| step | Ea refit `[D]` | Ea published | k_b refit `[D]` | k_b published | R² of refit | flag |
|---|---|---|---|---|---|---|
| Glc→Fru | 119.4 | 151.5 | 1230 | 2039 | 0.947 | mismatch 1.27× |
| Fru→Glc | 112.6 | 146.4 | 6425 | 10942 | 0.909 | mismatch 1.30× |
| Glc→3-DG | 59.6 | 107.2 | 2.43 | 4.19 | **0.698** | **non-monotonic (0.91→4.14→3.60)** |
| 3-DG→3,4-DG | 32.1 | 36.9 | 33.0 | 30.5 | 0.969 | ✔ agrees |
| **3,4-DG→HMF** | **−7.0** | **0** ᵟ | 133.8 | 119 | **0.189** | ⚠️ **non-monotonic (160→110→137)** |
| **Fru→Int** | **100.5** | **100.4** | **343.1** | **330** | **1.000** | ★ **agrees to 0.1 %** |
| **Int→HMF** | **145.0** | **151.4** | **1.844** | **1.84** | **1.000** | ★ **agrees to 4 %** |
| Fru→1-DG | 96.9 | 99.3 | 2.14 | 2.11 | 0.988 | ✔ agrees |
| Glc→G | 108.0 | 125.9 | 0.074 | 0.069 | 0.955 | mismatch 1.17× |
| G→GO | 75.1 | 93.8 | 790 | 737 | 0.927 | mismatch 1.25× |
| 3-DG→MGO | 93.7 | 84.8 | 314 | 304 | 0.997 | ✔ agrees |
| 1-DG→DA | 137.7 | 150.8 | 14.5 | 12.2 | 1.000 | ✔ agrees |
| 3-DG→P1 | 19.2 | 90.9 | 1022 | 2239 | **0.099** | ⚠️ non-monotonic, no signal |
| 1-DG→P2 | 75.4 | 77.9 | 891 | 873 | 0.999 | ✔ agrees |
| GO→P3 | −31.9 | **0** ᵟ | 24.0 | 32.6 | 0.099 | non-monotonic (authors set 0) |
| MGO→P4 | 44.3 | **0** ᵟ | 51.7 | 55.4 | 0.594 | non-monotonic (authors set 0) |
| **HMF→P6** | **96.7** | **152.8** | 55.4 | 36.9 | 0.910 | ⚠️ **mismatch 1.58×; published HPD ±154.4 > estimate** |

### GLUCOSE–NaCl

| step | Ea refit `[D]` | Ea published | k_b refit `[D]` | k_b published | R² | flag |
|---|---|---|---|---|---|---|
| Glc→Fru | 163.3 | 263.6 | 2239 | 5279 | 0.899 | ⚠️ mismatch **1.61×** |
| Fru→Glc | 170.6 | 280.8 | 11644 | 27705 | 0.902 | ⚠️ mismatch **1.65×** |
| Glc→3-DG | 99.4 | 85.2 | 1.20 | 1.11 | 0.980 | ✔ close |
| 3-DG→3,4-DG | 96.3 | 117.7 | 37.1 | 36.6 | 0.986 | ✔ close |
| **3,4-DG→HMF** | **94.2** | **169.8** | 151 | 103 | **0.997** | ⚠️ **mismatch 1.80× on a step whose own k trace R²=0.997** |
| **Fru→Int** | **102.1** | **94.9** | **1356** | **1402** | **1.000** | ★ **agrees to 3 %** |
| **Int→HMF** | **152.7** | **149.8** | **8.21** | **8.79** | **0.991** | ★ **agrees to 2 %** |
| Fru→1-DG | 58.5 | 93.9 | 1.12 | 1.36 | 0.950 | mismatch 1.61× |
| Glc→G | 110.4 | 131.8 | 0.069 | 0.059 | 0.974 | mismatch 1.19× |
| G→GO | 58.8 | 95.2 | 1263 | 1033 | 0.740 | mismatch 1.62× |
| 3-DG→MGO | 70.7 | 95.3 | 345 | 309 | 0.914 | mismatch 1.35× |
| 1-DG→DA | 88.7 | 139.9 | 22.9 | 16.1 | 0.757 | mismatch 1.58× |
| 3-DG→P1 | 58.9 | **0** ᵟ | 311 | 231 | 0.632 | non-monotonic |
| 1-DG→P2 | 5.0 | 55.3 | 451 | 574 | 0.204 | non-monotonic |
| GO→P3 | −113.6 | **0** ᵟ | 24.5 | 28.5 | 0.949 | monotonic **decreasing** |
| MGO→P4 | −39.7 | **0** ᵟ | 64.1 | 73.4 | 0.796 | non-monotonic |
| **HMF→P6** | **NOT COMPUTABLE** | 137.8 ± 24.6 | — | 227 | — | ⚠️ **k(160 °C) = 0.0 exactly ⇒ ln k = −∞. No Arrhenius line passes through this data. The published Ea = 137.8 cannot be derived from Table 1.** |

### 5.1 What the cross-check establishes

**★ FINDING 1 — the fructose limb is the only part of this paper whose temperature dependence is
internally reproducible.** `Fru→Int` and `Int→HMF` are the **only two steps that agree between the
naïve refit and the published global fit in BOTH systems**, and both have refit R² ≥ 0.991. Four
independent numbers (100.4/100.5 and 94.9/102.1 for k₆; 151.4/145.0 and 149.8/152.7 for k₇) land
within 7 %. This is a real result and it is the strongest single piece of HMF kinetics in the
entire K5a cluster.

**★ FINDING 2 — in the amine-free glucose system, the 3-DG limb's terminal step HAS NO
TEMPERATURE DEPENDENCE AT ALL, AND THE AUTHORS SAY SO.** `3,4-DG→HMF` runs **160 → 110 → 137**
(×10⁻³ min⁻¹) across 160→180→200 °C. My refit: Ea = **−7.0 kJ/mol, R² = 0.189**. The authors' own
Table 2 sets **Ea = 0 with the δ footnote "does not follow Arrhenius equation"**. Adding NaCl
restores monotonicity (46.0→163→418) but the published Ea (169.8) is **1.80× my refit (94.2)**
even though that k trace is nearly perfectly log-linear (R² = 0.997) — meaning the global fit is
pulling this parameter away from what its own compound's data support. **⇒ No Ea for the 3-DG→HMF
edge is defensible from this paper in either system.**

**★ FINDING 3 — the HMF SINK Ea is not usable in either system.**
Glucose: published `152.8 ± 154.4` — **the 95 % HPD half-width exceeds the estimate**, the same
`SD ≥ estimate` pathology the repo already refused in Knol (`k3` §C.6) and in Han 2025 k8/k15.
NaCl: `k(160 °C) = 0.0 ± 0.0` exactly, so **no Arrhenius fit is possible from the published k**,
yet Table 2 prints `137.8 ± 24.6`. **⇒ REFUSE both. The corpus still has no usable HMF-sink Ea from
this paper** (see `hamzalioglu2018_extraction.md` for the one that is usable, at 5–50 °C).

**★ FINDING 4 — the glucose↔fructose isomerisation Ea are the paper's own admitted artefact.**
Published: `151.5 ± 99.8` and `146.4 ± 104.2` (glucose); `263.6` and `280.8` (NaCl). HPD is 66 %
and 71 % of the estimate in the glucose system, and both `k_b` are `ind`. The authors state,
verbatim (lines 312–330):

> "The temperature dependence of the interconversion was higher in the presence of NaCl. This
> seemed interesting because the rate constants were higher at 180 and 200 °C. But it was due to
> the **limitation of Arrhenius equation, which does not take the physical conditions of the system
> into consideration** … NaCl **decreased the matrix molecular mobility** … **Therefore, the
> activation energy estimated here should not be considered as measure of an energy barrier for the
> reaction.** The limitations of Arrhenius equation in food systems and complex reactions were well
> discussed by Peleg, Normand and Corradini (2012)."

**This is an author-declared prohibition on reading their own Ea as barriers.** It must travel with
any ingestion, in the same register as `k3` §C.1 and the Lee 2024 prohibition in round-3 §A.3.4.

---

## §6. THE NaCl EFFECT — the paper's headline, re-derived as ratios `[D]` from Table 1

Ratios cancel the `Int` scale factor (§2) **and** the semi-quantitation response factors (§1.1),
so they are the transferable form. Computed as k(NaCl)/k(glucose) at each temperature:

| step | 160 °C | 180 °C | 200 °C | verdict |
|---|---|---|---|---|
| 1 Glc→Fru | 0.89 | **2.61** | **2.47** | matches the abstract's *"2.5 times faster … at 180 and 200 °C"* ✔; **no effect at 160 °C** |
| 2 Fru→Glc | 0.78 | 2.40 | 3.00 | same pattern (but k is `ind` at 180/200) |
| 3 Glc→3-DG | **0.43** | **0.24** | 1.13 | **decrease at 160/180, REVERSES at 200 °C** |
| 4 3-DG→3,4-DG | 0.46 | 1.42 | 2.05 | sign flip between 160 and 180 °C |
| 5 3,4-DG→HMF | 0.29 | 1.48 | 3.05 | sign flip |
| **6 Fru→Int** | **3.91** | **3.88** | **4.06** | ★ **3.9–4.1×, flat across 40 °C — the most stable ratio in the paper** |
| **7 Int→HMF** | **3.71** | **5.35** | **4.41** | ★ **3.7–5.4×** |
| 8 Fru→1-DG | 0.82 | 0.54 | 0.33 | decrease, monotone |
| 9 Glc→G | 0.87 | 1.00 | 0.92 | **no effect** (≈1 at all T) |
| 11 3-DG→MGO | 1.74 | 0.76 | 1.03 | no consistent effect |
| 18 HMF→P6 | 0.0 | 7.17 | 4.70 | not interpretable (k=0 at 160) |

**★ The abstract's "the rate constants increase 4 fold in the presence of NaCl" is CONFIRMED and is
temperature-stable for k₆ (3.91 / 3.88 / 4.06).** This is a genuinely well-behaved measured effect
and — being a ratio within one lab, one method, one pair of tubes — it is **immune to the `Int`
identifiability problem of §2**. It is the highest-quality NaCl number in the cluster.

**⚠️ The abstract's "A decrease in rate constants of 3-deoxyglucosone and 1-deoxyglucosone
formations by the presence of NaCl" is NOT temperature-general.** k₃ falls 2.3× at 160 °C and
4.2× at 180 °C but **rises 1.13× at 200 °C.** The authors concede it: *"There was no significant
difference at 200 °C, due to the higher uncertainty for the parameter in glucose-NaCl system"*
(lines 351–352). Quote the ratio with its temperature or not at all.

### 6.1 Absolute HMF yields — mole conversions, `[M]`, verbatim lines 235–242

> "The amount of HMF formed in the presence of NaCl was **almost 10 fold higher** for the same
> time-temperature treatment. The **maximum mole conversions of glucose to HMF** were **0.4 % at
> 160 °C (30 min), 1.6 % at 180 °C (20 min) and 3.5 % at 200 °C (15 min)** during caramelization of
> glucose only. These maximum conversions in the presence of NaCl altered to **1.4 % at 160 °C
> (20 min), 3.1 % at 180 °C (10 and 15 min) and 3.7 % at 200 °C (3 min)**."

★ **These six numbers are the single most valuable ingestible quantity in this paper**: they are
**dimensionless mole yields**, so they are immune to the µmol-per-tube normalisation, the initial-load
mismatch of §1, and every response-factor caveat. NaCl/glucose conversion ratio: **3.5× (160 °C),
1.9× (180 °C), 1.06× (200 °C)** — i.e. **the NaCl enhancement of the HMF YIELD collapses to nothing
by 200 °C**, even though the rate constant ratio stays at ~4×. The difference is that at 200 °C the
NaCl system also destroys HMF ~4.7× faster (k₁₈). **A model that carries only the formation
enhancement and not the sink enhancement will over-predict NaCl-HMF at high temperature by ~3×.**

### 6.2 Other measured ratios `[M]`, verbatim lines 228–231

> "In the glucose system, the average ratio of the amount of **3-deoxyglucosone to
> 1-deoxyglucosone** was **4.6 ± 0.5, 4.0 ± 1.2 and 3.9 ± 1.0** at 160, 180 and 200 °C … The ratio
> was slightly increased to **5.6 ± 0.3, 5.2 ± 1.5 and 5.3 ± 1.2** in the glucose-NaCl system"

(Compare the wheat-flour paper's 10.6 ± 2.8 / 7.8 ± 2.7 / 6.9 ± 1.5 — same lab, same method,
**2.3× higher** in the amine-containing flour matrix. See the synthesis §5.)

---

## §7. DIRECTIONAL / STRUCTURAL CONSTRAINTS — `STRUCTURAL`, no number to fit

| # | constraint | verbatim anchor |
|---|---|---|
| S1 | **HMF is primarily formed from fructose dehydration in an amine-free glucose melt** | *"The estimated rate constants clearly indicated that HMF is primarily formed from fructose dehydration."* (line 433) |
| S2 | **A 3-DG-only HMF model cannot reproduce the observed HMF** (Fig. S5 test) | lines 424–432, §2.1 |
| S3 | **Fru→3-DG is zero.** 3-DG comes only from glucose in this system. | lines 344–345 |
| S4 | **1-DG comes only from fructose** (2,3-enolisation), *"which does not occur from glucose"* | lines 371–374 |
| S5 | **MGO comes from 3-DG here, but from 1-DG in the wheat-flour system** — the authors explicitly reconcile the two | lines 396–406: *"the main source of methylglyoxal may quantitatively depend on the amount of precursor α-dicarbonyl compound formed"* ★ **this is a same-lab, same-method demonstration that a branch source is matrix-determined, not a constant** — the exact claim of `k3` §0 |
| S6 | **Isomerisation is 5× faster in reverse** (Fru→Glc / Glc→Fru ≈ 5) in **both** systems at **all** temperatures | line 284–286; check from Table 1: 5.4/5.8/4.6 (glucose), 4.7/5.3/5.6 (NaCl) `[D]` ✔ |
| S7 | **NaCl switches glucose degradation from the open-chain α-dicarbonyl route to the cyclic route** | lines 463–467: *"degradation of glucose switch to cyclic intermediates in the presence of sodium chloride, as evident from decrease in the rate of α-dicarbonyl compounds formation and elevation of fructose degradation to HMF through cyclic intermediates"* |
| S8 | **Volatiles leave the system.** GO, MGO and DA degradation rates *decrease* with temperature; the authors attribute it to headspace loss, not chemistry | lines 384–388: *"attempted to higher volatility of the compound, which diffuse to the headspace of the tube during heating"* — ⚠️ **this is a transport artefact contaminating k₁₅, k₁₆, k₁₇, exactly the Lee 2024 defect class of round-3 §A.3.4. Those three Ea are set to 0 by the authors for this reason.** |
| S9 | Oxidation is quantitatively minor vs dehydration | line 336 |

---

## §8. WHAT THIS PAPER DOES **NOT** CONTAIN — verified negatives `[NEG]`

- **No amino acid, no Amadori product, no melanoidin, no browning measurement.** (By design — it is
  a pure caramelization study. This is its greatest strength for the repo and its main limit.)
- **No pH.** The system is a freeze-dried solid; no pH is measured or reported anywhere.
- **No water activity, no moisture content.** The amorphous state is explicitly *not characterised*.
- **No 1,2-enediol node.**
- **No furfural.** Only HMF among the furanics.
- **No absolute HMF concentration table.** HMF vs time exists only as Figure 3/4 rasters
  (y-axis max 1.8 µmol for glucose, 2.5 µmol for glucose-NaCl — readable from the axis labels in
  the text layer, but the curves themselves are not in the text layer). **The mole conversions of
  §6.1 are the only extractable absolute HMF numbers.**
- **Supporting Information (Figures S1–S5) is NOT in the PDF on disk.** The model-discrimination
  figure S5 that underpins constraint S2 cannot be inspected.

---

## §9. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| Table 1, all 102 k ± HPD cells | `[F]` | **USE-Q** — as *per-temperature* constants of **this** system only; `min⁻¹×10³`; drop the 13 cells of §3.1 |
| **k₆, k₇ (fructose limb) absolute values** | `[F]` | ⚠️ **RATIO-ONLY** — `[Int]` is unmeasured, so only `k₆·[Int]` and `k₇·[Int]` are identified |
| **k₆(NaCl)/k₆(glucose) = 3.91 / 3.88 / 4.06** | `[D]` | ★ **USE** — temperature-stable, scale-free, the paper's best number |
| **k₇(NaCl)/k₇(glucose) = 3.71 / 5.35 / 4.41** | `[D]` | **USE-Q** — 1.4× spread over 40 °C |
| **k₃ ratio (3-DG formation)** | `[D]` | **USE-Q — quote with temperature. Sign flips at 200 °C.** |
| Table 2 Ea for steps 6 and 7 (100.4 / 94.9 and 151.4 / 149.8 kJ/mol) | `[F]`, reproduced `[D]` | ★ **USE-Q** — *"apparent Ea for a lumped fructose-dehydration-to-HMF step in an amine-free amorphous glucose melt, 160–200 °C, 3 points, over an unmeasured intermediate pool"*. Never "the HMF barrier". |
| **Table 2 Ea for step 5 (3,4-DG→HMF)** | `[F]` | **REFUSE.** Glucose: authors set it to 0 with the "does not follow Arrhenius" footnote. NaCl: published 169.8 vs refit 94.2 = **1.80× discrepancy**. |
| **Table 2 Ea for step 18 (HMF sink)** | `[F]` | **REFUSE.** Glucose: HPD (±154.4) > estimate (152.8). NaCl: k(160) = 0 ⇒ **not derivable from the paper's own table.** |
| **Table 2 Ea for steps 1 and 2 (isomerisation)** | `[F]` | **REFUSE — author-declared not a barrier** (§5 Finding 4, verbatim). |
| Table 2 Ea for steps 13/15/16/17 | `[F]` | **REFUSE** — either fixed to 0 by the authors, or contaminated by headspace loss (S8). |
| Table 2 Ea for steps 3, 4, 8–12, 14 | `[F]` | **USE-Q** — not HMF-lane; carry the §5 refit column as a companion uncertainty. Steps 4, 8, 11, 12, 14 agree between refit and publication and are the sound ones. |
| **Mole conversions (0.4/1.6/3.5 % and 1.4/3.1/3.7 %)** | `[M]` | ★ **USE — the best absolute HMF quantity in the cluster.** |
| 3-DG:1-DG ratios (4.6/4.0/3.9 and 5.6/5.2/5.3 ± SD) | `[M]` | **USE** |
| S1–S9 | `[M]`/`[F]` | **STRUCTURAL** |
| 3,4-DG, 1-DG, glucosone concentrations | `[M]` | **`absolute_concentration: false`** — semi-quantitated on the 3-DG curve |
