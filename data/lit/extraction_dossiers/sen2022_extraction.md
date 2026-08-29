# Şen & Gökmen 2022 — Maillard + caramelization in sucrose-rich, low-moisture nuts and seeds

### Wave K5a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### Purpose in the wave brief: **BOUND the matrix dependence of the HMF branch.** It does.

---

## §0. IDENTITY

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/Sen2022.pdf` (1,477,351 bytes) | `ls` |
| **DOI** | **`10.1016/j.foodchem.2022.133583`** ✔ matches the round-3 expectation | p.1 |
| title | *"Kinetic modeling of Maillard and caramelization reactions in sucrose-rich and low moisture foods applied for roasted nuts and seeds"* | p.1 |
| authors | **Dilara Şen**, Vural Gökmen | p.1 |
| journal | *Food Chemistry* **395** (2022) 133583 | running head |
| dates | Received 13 Mar 2022 · Revised 20 Jun 2022 · Accepted 24 Jun 2022 · online 27 Jun 2022 | p.1 |
| version | **published version of record** | p.1 |

Filename/DOI mapping is **correct** for this file.

---

## §1. SYSTEM AND CONDITIONS `[M]`

| item | value | anchor |
|---|---|---|
| matrices | **sunflower seed** (*Helianthus annuus*), **pumpkin seed** (*Cucurbita moschata*), **flaxseed** (*Linum usitatissimum*), **peanut** (*Arachis hypogaea*), **almond** (*Prunus dulcis*) — bought from a local market, Ankara | §2.2 |
| roasting | **30 g** on an aluminium plate, conventional oven (Memmert UN 55), **160 and 180 °C**, **5 to 60 min** | §2.2 |
| replication | **triplicate** for each nut × temperature × time | §2.2 |
| model units | **µmol/kg** | §2.11 |
| statistics | one-way ANOVA + **Duncan's test (p < 0.05)**, SPSS 16.0; PCA in XLSTAT | §2.12 |

### 1.1 ★ THE MATRIX ENVELOPE — Table 1, transcribed verbatim (p.4)

Caption verbatim: *"Proximate composition, pH, and water activities of raw samples."*
Superscript letters: *"Different letters on the same line indicate statistically significant
differences."*

| | Sunflower Seed | Pumpkin Seed | Flaxseed | Peanut | Almond |
|---|---|---|---|---|---|
| Moisture (g/100 g) | 5.26 ± 0.2ᵃ | **6.7 ± 0.07ᵇ** | 3.86 ± 0.15ᶜ | 4.62 ± 0ᵈ | 4.01 ± 0.02ᶜ |
| Oil (g/100 g) | 42.24 ± 14.07ᵃ | 49.08 ± 1.33ᵃ | 42.94 ± 3.04ᵃ | 43.77 ± 0.19ᵃ | 48.02 ± 0.24ᵃ |
| Protein (g/100 g) | 20.11 ± 4.3ᵇ | 17.36 ± 1.56ᵃ,ᵇ | 12.76 ± 1.62ᵃ | 12.61 ± 0.42ᵃ | 20.48 ± 2.17ᵇ |
| Carbohydrate (g/100 g) | 28.32 ± 18.07ᵃ | 22.61 ± 2.82ᵃ | 36.65 ± 4.85ᵃ | 36.39 ± 0.62ᵃ | 24.31 ± 2.43ᵃ |
| Ash (g/100 g) | 4.08 ± 0.16ᵃ | 4.24 ± 0.01ᵃ | 3.79 ± 0.08ᵇ | 2.6 ± 0.02ᶜ | 3.18 ± 0ᵈ |
| **Water activity (25 °C)** | 0.48 ± 0.02ᵇ,ᶜ | **0.62 ± 0ᵈ** | 0.43 ± 0.02ᵃ | 0.46 ± 0.01ᵃ,ᵇ | 0.5 ± 0.01ᶜ |
| **pH** | **7.1 ± 0.03ᶜ** | **7.14 ± 0.02ᶜ** | **6.62 ± 0.03ᵃ** | **7.01 ± 0.04ᵇ** | **6.62 ± 0.01ᵃ** |

★ **This is the corpus's cleanest same-method pH + a_w + composition envelope for the HMF lane:
five matrices, pH 6.62–7.14, a_w 0.43–0.62, oil 42–49 %, at two temperatures.** Compare the
matched conditions elsewhere in the cluster: Gürsul Aktağ pH **3.4** at 27–37 °C, Kocadağlı flour
7 % moisture / pH unmeasured at 160–200 °C, Kocadağlı JAFC amine-free freeze-dried melt.

### 1.2 Initial reactant loads `[M]`

| quantity | value | anchor |
|---|---|---|
| sucrose (raw) | sunflower **3.2 ± 0.0**, pumpkin **1.3 ± 0.0**, flaxseed **2.3 ± 0.0**, peanut **4.0 ± 0.2**, almond **3.4 ± 0.1** g/100 g | §3.2 |
| max sucrose loss (180 °C, 40 min) | **84, 50, 58, 49, 46 %** (same order) | §3.2 |
| total free amino acids (raw) | **2131.7 ± 601.1 to 4217.7 ± 493.1 mg/kg** | §3.2 |
| protein-bound lysine (raw) | highest **peanut 11719.4 ± 835.0**, **pumpkin 11771.0 ± 1313.6**; lowest **almond 6933.7 ± 189.1 mg/kg** | §3.2 |
| total amino-acid loss on roasting | flaxseed **72 %**, sunflower **71 %**, peanut **44 %**, almond **35 %**, pumpkin **28 %** | §3.2 |
| furosine (raw) | sunflower **45.6 ± 0.1**, pumpkin **195.6 ± 49.9**, flaxseed **15.6 ± 1.0**, peanut **47.1 ± 4.0** mg/kg; **almond: not detectable** | §3.2 |
| **FL conversion factor** | ⚠️ **`[FL] = furosine × 2.2`** — *"Since furosine are generated through degradation of N-ε-fructoselysine (FL) to furoyl derivatives with acid hydrolysis, FL content was calculated via multiplying furosine concentration by 2.2 (Berk, Aktağ, & Gökmen, 2021)"* | §2.8 `[C]` |

⚠️ **`[FL]` — one of the network's central species — is not measured; it is a cited ×2.2 stoichiometric
inflation of furosine.** Any k involving FL inherits that borrowed factor.

### 1.3 ★ ABSOLUTE HMF, and the α-dicarbonyl envelope `[M]`

| quantity | value | anchor |
|---|---|---|
| **max HMF** | ★ **247.0 ± 8.3 mg/kg in sunflower seed, 30 min at 180 °C** (declared *"highest (p < 0.05)"*) | §3.2 |
| **min HMF** | **42.7 ± 0.4 mg/kg in pumpkin seed** — *"due to the least amount of sucrose among other samples"* | §3.2 |
| trend | *"an overall increasing trend in all samples in relation to the decrease in sugar content"* | §3.2 |
| 3-DG | **> 10 mg/kg in all roasted samples except peanut** (slightly lower) | §3.2 |
| glucosone | **the least abundant** dicarbonyl in every roasted sample | §3.2 |
| GO / DMG max | **pumpkin seed, 10.3 and 5.9 mg/kg** | §3.2 |
| max CML | **43.0 ± 3.1 mg/kg, sunflower**, generally at 180 °C within 30–40 min | §3.2 |
| max CEL | **164.1 ± 9.8 mg/kg, pumpkin @ 180 °C × 40 min**; sunflower **132.8 ± 5.0**, flaxseed **132.6 ± 8.6**, peanut **150.0 ± 5.9**, almond **116.2 ± 11.3 mg/kg** | §3.2 |
| ordering | **CEL > CML** in every matrix; **furosine > CML** in every matrix | §3.2 |

★ HMF measured by HPLC-DAD at 285 nm against **authentic HMF (98 %, Acros), 1–10 mg/L** — an
absolute, standard-calibrated number.
⚠️ α-dicarbonyls, CML and CEL were quantified against authentic standards too (3-DG 75 %, MGO 40 %,
quinoxalines), **but CML/CEL used matrix-matched calibration with raw samples as blank** (§2.10),
i.e. the raw-sample background is subtracted, not zero.

---

## §2. THE REACTION NETWORK — 24 constants, TWO models (nuts/seeds vs pumpkin seed)

Abbreviations verbatim (Table 2 footnote): `SUC` sucrose; `GLU`/`GLC` glucose;
**`FFC` fructofuranosyl cation**; `HMF`; `bLYS` protein-bound lysine; **`FL` N-ε-fructoselysine
(the Amadori product)**; **`HP` Heyns product**; `3-DG`; `GO`; `G` glucosone; `1-DG`; `MGO`;
`DMG` dimethylglyoxal; `CML` N-ε-carboxymethyllysine; `CEL` N-ε-carboxyethyllysine; `AA` total
amino acids; `P` products.

**★ Two distinct networks were fitted.** PCA (Fig. 1) separated **pumpkin seed** from the other
four on water activity, sucrose and free amino acids, and it gets its own model — verbatim §3.3:
*"it is envisaged that all samples can be explained with a single model, whereas probably different
reaction pathways will be dominant in pumpkin seed samples."*

**★ THE HMF NODE ODE — Appendix A, Eq. (5) for nuts/seeds and Eq. (28) for pumpkin seed, verbatim
and IDENTICAL in both:**

> **`d[HMF]/dt = k₂·[FFC] + k₁₀·[3-DG] − k₁₉·[HMF]·[AA]`**

⚠️⚠️ **THE HMF SINK HERE IS BIMOLECULAR IN HMF AND TOTAL AMINO ACIDS.** This is the **only**
multiresponse paper in the K5a cluster with an amine-dependent HMF sink — and it is **the exact
step that the same lab's Kocadağlı & Gökmen 2016 wheat-flour paper TESTED AND REJECTED**
(*"HMF did not show a good fit when its degradation included amino acids in a bimolecular reaction.
Therefore, this step was excluded"*). Same lab, opposite architecture. See the synthesis §4.
It is also the step that **Hamzalıoğlu & Gökmen 2018 measures directly** — see that dossier.

Full nuts/seeds ODE set, verbatim from Appendix A:
```
(1)  d[SUC]/dt  = −k1[SUC]
(2)  d[Lys]/dt  = −(k8+k21)[HP] + (k9+k18)[FL] − k3[FFC][Lys] − k4[GLC][Lys] − k14[GO][Lys] − k16[MGO][Lys]
(3)  d[AA]/dt   = −k17[GLC][AA] − k19[HMF][AA] − k20[FFC][AA] − k22[3−DG][AA]
(4)  d[FL]/dt   = −k4[GLC][Lys] − (k9+k13+k18)[FL]           ⚠️ see §2.1
(5)  d[HMF]/dt  =  k2[FFC] + k10[3−DG] − k19[HMF][AA]
(6)  d[3−DG]/dt =  k5[GLC] + k8[HP] − k10[3−DG] − k22[3−DG][AA]
(7)  d[G]/dt    =  k7[GLC] − k24[G]
(8)  d[DMG]/dt  =  k12[1−DG] − k23[DMG]
(9)  d[MGO]/dt  =  k11[1−DG] − k16[MGO][Lys]
(10) d[GO]/dt   =  k6[GLC] + k24[G] − k14[GO][Lys]
(11) d[CML]/dt  =  k13[FL] + k14[GO][Lys]
(12) d[CEL]/dt  =  k15[HP] + k16[MGO][Lys]
(13) d[FFC]/dt  =  k1[SUC] − k2[FFC] − k3[FFC][Lys] − k20[FFC][AA]
(14) d[GLC]/dt  =  k1[SUC] − k4[GLC][Lys] − (k5+k6+k7)[GLC] − k17[GLC][AA]
(15) d[1−DG]/dt =  k9[FL] − (k11+k12)[1−DG]
(16) d[HP]/dt   =  k3[FFC][Lys] − (k8+k15+k21)[HP]
(17–23) d[P1..P7]/dt = k17[GLC][AA], k18[FL], k19[HMF][AA], k20[FFC][AA], k21[HP], k22[3−DG][AA], k23[DMG]
```
Pumpkin-seed differences (Eqs. 24–37): `k7` is redefined as **FL → GO + bLYS** (not GLC → G);
`k24` is redefined as **HP → 1-DG + bLYS** (not G → GO); `k21` becomes **1-DG → P5b**;
`k20` becomes a first-order FFC drain rather than `FFC + AA → P4`; the G and GO-from-glucosone
nodes are dropped. **The glucosone node does not exist in the pumpkin-seed model.**

### 2.1 ⚠️ PRINTED SIGN ERROR IN Eq. (4) AND Eq. (27)

Eq. (4) reads `d[FL]/dt = − k4[GLC][Lys] − (k9+k13+k18)[FL]`. **The `k4` term must be `+`** — `k4`
is the *formation* of FL from glucose and bound lysine (Table 2 names it "GLC + bLYS → FL"), and it
appears with a minus sign in Eq. (2) for `[Lys]`, which is only consistent if it is a source for FL.
As printed, FL has **no source at all** and can only decay from an initial value. The same error
appears in Eq. (27) (pumpkin). **Transcribed verbatim above; flagged as a typesetting error, not a
model choice.**

### 2.2 What model discrimination removed `[M]`

- *"the reaction of all α-dicarbonyl compounds with free amino acids was also included … However,
  when all the product formation steps were included in the model, the predicted data were not fit
  well … Finally, the **self-reaction of DMG (k23) and the reaction of 3-DG with amino acids (k22)**
  were found to be important."* (§3.4.3)
- MGO formation: only the **1-DG** route survives — *"the predicted concentration data were well
  fitted with the experimental data **when only the formation of MGO through 1-DG pathway (k11) was
  included**"* (§3.4.3).
- **3-DG formation through the FL step was excluded** *"as it was estimated in lower bound"*.
- **GO formation from FL was removed for all but pumpkin seed.**
- **HP → FL interconversion was excluded** *"Since the increase in the number of unmeasurable
  reactants and products in the model, such as HP, reduces the possibility of accurate estimates"*.

⚠️ **`FFC` AND `HP` ARE BOTH UNMEASURED**, and the authors name this as the cause of the paper's
worst uncertainties, verbatim §3.4:

> "However, **some rate constants could not be estimated within the specified confidence interval.
> This can be explained by the fact that the relevant rate constants are in the steps which include
> unquantified compounds, such as reactions involving HP, GLC, or FFC.**"

⇒ **k₂ (FFC→HMF) is identified only up to the scale of the unmeasured FFC pool.** k₁₀ (3-DG→HMF)
is **not** subject to this: 3-DG is measured against an authentic standard. **The two arms of the
HMF branch in this paper are not epistemically symmetric.**

---

## §3. TABLE 2 — the complete constant set, ± 95 % HPD

**Transcribed verbatim from p.7.** Caption verbatim: *"Estimated reaction rate constants
(k, min-1x10³) with 95 % highest posterior density (HPD) intervals at different temperatures
according to the proposed model for Maillard reaction and caramelization during roasting of
samples."* Footnote `*`: *"Reaction rate constants in kg·µmol⁻¹·min⁻¹ × 10³."* (applies to
k17, k19, k20, k22). Footnote `***ind`: *"Indeterminate, which means a large uncertainty in the
estimated parameter within 95 % confidence interval."* Provenance **`[F]`**.

`–` = the step does not exist in that matrix's model.

| k | elementary step | Sunf 160 | Sunf 180 | Flax 160 | Flax 180 | Peanut 160 | Peanut 180 | Almond 160 | Almond 180 | Pumpkin 160 | Pumpkin 180 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| k1 | SUC → FFC + GLC | 11.24 ± 1.31 | 29.43 ± 3.3 | 11.72 ± 1.16 | 33.37 ± 30.19 | 3.72 ± 0.5 | 14.62 ± 3.94 | 9.09 ± 0.65 | 23.07 ± 2.96 | 5.34 ± 0.4 | 14.78 ± 2.15 |
| **k2** | **FFC → HMF** | **0 ± 0** | **89.49 ± 44.3** | **0 ± 0** | **0 ± 0** | **13.2 ± 15.83** | **0.85 ± 1.25** | **0 ± 0** | **2.77 ± 6.3** | **0 ± 0** | **0 ± 0** |
| k3 | FFC + bLYS → HP | 96.83 ± 44.55 | 61.50 **ind** | 27.3 ± 12.54 | 5.51 ± 30.36 | 22.04 **ind** | 1.94 ± 2.75 | 0.78 ± 0.3 | 1.02 ± 0.29 | 1.98 ± 12.45 | 29.3 ± 43.3 |
| k4 | GLC + bLYS → FL | 5.43 **ind** | 8.98 ± 8.22 | 0.40 **ind** | 2.72 **ind** | 11.91 ± 3.85 | 212.65 ± 115.3 | 50.3 **ind** | 8.01 ± 9.78 | 1.87 ± 0.43 | 0.13 ± 0.15 |
| k5 | GLC → 3-DG | 38.23 ± 80.83 | 5.10 **ind** | 34.94 ± 21.88 | 208.64 ± 107.2 | 15.09 ± 22.47 | 0 ± 0 | 544.68 ± 535.8 | 2.53 **ind** | 10.91 ± 1.4 | 8.12 ± 10.36 |
| k6 | GLC → GO | 0.17 ± 2.4 | 2.12 ± 1.36 | 1.21 ± 0.28 | 8.51 ± 6.9 | 4.35 ± 2.1 | 24.97 **ind** | 3.13 ± 2.58 | 0.27 ± 0.52 | 4.3 ± 2.34 | 0.49 ± 0.57 |
| k7a | GLC → G | 202.14 ± 418.3 | 298.11 ± 524.9 | 885.53 ± 341.1 | 2295.25 ± 557.4 | 1.64 ± 0.61 | 5.8 ± 2.95 | 0.22 ± 0.07 | 0.03 ± 0.03 | – | – |
| k7b | FL → GO + bLYS | – | – | – | – | – | – | – | – | 0 ± 0 | 5.36 ± 1.83 |
| k8 | HP → 3-DG + bLYS | 0.57 ± 3.76 | 0 ± 0 | 2.67 ± 1.8 | 0 ± 0 | 0.5 ± 0.36 | 2.23 ± 1.16 | 0 ± 0 | 0.96 ± 1.77 | 0 ± 0 | 0 ± 0 |
| k9 | FL → 1-DG + bLYS | 364.36 ± 365.5 | 96.26 ± 93 | 18.17 ± 7.83 | 111.85 ± 22.7 | 28.22 ± 3.29 | 27.83 ± 13.51 | 17.11 ± 2.19 | 46.49 ± 11.22 | 180.61 ± 10.78 | 156.17 ± 31.28 |
| **k10** | **3-DG → HMF** | **568.79 ± 143.8** | **109.98 ± 90.19** | **1024.28 ± 612.3** | **918.19 ± 431.4** | **279.32 ± 180** | **807.82 ± 291.8** | **2528.66 ± 1625** | **601.45 ± 822** | **37.53 ± 1.19** | **552.83 ± 278.7** |
| k11 | 1-DG → MGO | 5.36 ± 5.63 | 20.73 ± 26.49 | 217.84 ± 185.6 | 232.42 ± 172.4 | 1266.25 ± 656.4 | 342 ± 290.5 | 124.06 ± 26.3 | 218.46 ± 164 | 18.02 ± 3.06 | 14.4 **ind** |
| k12 | 1-DG → DMG | 94 ± 34.71 | 5.27 ± 4.89 | 69.49 ± 57.57 | 18.6 ± 10.93 | 199.59 **ind** | 29.65 ± 26.09 | 11.82 ± 5.48 | 7.7 ± 5.95 | 6.36 ± 1.75 | 5.36 ± 1.11 |
| k13 | FL → CML | 0.8 ± 1.7 | 0 ± 0 | 6.87 ± 0.78 | 11.94 ± 8.12 | 3.43 ± 1.53 | 7.4 ± 2.05 | 1.07 ± 1.77 | 4.12 ± 7.52 | 3.44 ± 1.44 | 0.98 ± 2.86 |
| k14 | GO + bLYS → CML | 0.56 ± 0.35 | 2.7 ± 0.56 | 0 ± 0 | 1.96 ± 2.34 | 0.53 ± 0.35 | 0.41 ± 0.51 | 1.05 ± 1.37 | 1.68 ± 4.69 | 0.04 ± 0.04 | 0.23 ± 0.09 |
| k15 | HP → CEL | 0 ± 0 | 0 ± 0 | 0.55 ± 0.22 | 0 ± 0 | 0 ± 0 | 0.96 ± 0.49 | 0 ± 0 | 0 ± 0 | 1.67 ± 2.49 | 3.03 ± 1.36 |
| k16 | MGO + bLYS → CEL | 2.16 ± 0.51 | 5.12 ± 2.22 | 0.86 ± 0.82 | 9.87 ± 2.1 | 1.8 ± 0.15 | 1.31 ± 0.82 | 2.56 ± 0.33 | 5.37 ± 1.18 | 0.08 ± 0.07 | 0.13 ± 0.09 |
| k17\* | GLC + AA → P1 | 34.37 ± 32.36 | 54.6 ± 45.15 | 52.15 ± 14.04 | 327.74 ± 504.1 | 61.84 **ind** | 463.95 ± 416.1 | 0 ± 0 | 0 ± 0 | 50.8 **ind** | 9.5 ± 13.06 |
| k18 | FL → P2 + bLYS | 0 ± 0 | 974.55 ± 1172 | 0 ± 0 | 0 ± 0 | 525.68 ± 83.19 | 2243.87 ± 314.4 | 1275.86 ± 324.3 | 4742.99 ± 924.9 | 3.3 **ind** | 37.46 ± 33.29 |
| **k19\*** | **HMF + AA → P3** | **9.85 ± 5.59** | **0 ± 12.9** | **26.48 ± 29.5** | **13.35 ± 27.76** | **0 ± 0** | **0 ± 0** | **29.29 ± 25.55** | **24.79 ± 62.13** | **0 ± 0** | **43.66 ± 28.72** |
| k20\* | FFC + AA → P4 | 317.56 **ind** | 183.93 ± 133.2 | 72.74 **ind** | 2.68 ± 39.94 | 51.53 ± 40.31 | 3.44 ± 8.73 | 0.71 ± 0.99 | 2.27 ± 1.27 | 0 ± 0 | 98.3 **ind** |
| k21a | HP → P5a + bLYS | 6.75 ± 7.4 | 16.52 ± 14.19 | 2.85 ± 5.94 | 8.3 ± 21.12 | 0.2 ± 13.55 | 17.82 ± 24.66 | 7.46 ± 8.81 | 0 ± 0 | – | – |
| k21b | 1-DG → P5b | – | – | – | – | – | – | – | – | 475.6 **ind** | 179.69 ± 78.25 |
| k22\* | 3-DG + AA → P6 | 9.89 ± 61.77 | 3.13 ± 3.7 | 0 ± 0 | 85.75 ± 58.98 | 24.67 ± 45.31 | 42.3 ± 109.4 | 369.63 ± 300.9 | 36.06 ± 44.36 | 0 ± 0 | 0 ± 0 |
| k23 | DMG → P7 | 18591.15 ± 19730 | 229.67 ± 237.2 | 41.88 ± 18.15 | 23.34 ± 19.55 | 90.85 ± 52.51 | 32.02 ± 25.85 | 45.79 ± 29.39 | 16.8 ± 12.37 | 8.7 ± 5.12 | 16.28 ± 6.97 |
| k24a | G → GO | 0.22 ± 0.4 | 0.02 ± 0.16 | 0 ± 0 | 0 ± 0 | 40.39 ± 20.9 | 0 ± 0 | 0 ± 0 | 0 ± 0 | – | – |
| k24b | HP → 1-DG + bLYS | – | – | – | – | – | – | – | – | 1.32 ± 1.92 | 0 ± 0 |

### 3.1 ⚠️ HOW MUCH OF TABLE 2 IS USABLE

**12 cells are `ind`.** Counting cells whose HPD half-width **equals or exceeds** the estimate
(⇒ the 95 % interval spans zero): I count **≥ 40 of the ~230 numeric cells**, including
**k10 almond @180 (601.45 ± 822 — the interval spans zero on the paper's headline step)**,
k10 sunflower @180 (109.98 ± 90.19, 82 %), k10 almond @160 (2528.66 ± 1625, 64 %),
k2 peanut @160 (13.2 ± 15.83), k2 almond @180 (2.77 ± 6.3), k19 sunflower @180 (0 ± 12.9),
k19 flaxseed both, k19 almond @180 (24.79 ± 62.13), k5 sunflower @160 (38.23 ± 80.83),
k5 almond @160 (544.68 ± 535.8), k22 sunflower @160, peanut both, almond @180,
k23 sunflower @160 (18591.15 ± 19730). **⇒ roughly a fifth of the table carries no signal, and the
contamination reaches the paper's own headline conclusion.**

⚠️ Note **k23 sunflower @160 °C = 18 591 ×10⁻³ min⁻¹ = 18.6 min⁻¹**, which is **~340× the same step
at 180 °C (0.230 min⁻¹)** and 7–2000× every other constant. Stiff-unidentified-pair signature.

---

## §4. ★ THE DISSENT, QUANTIFIED — and correctly bounded

### 4.1 What the paper claims

Abstract, verbatim: *"Accordingly, **3-deoxyglucosone formation via sugar degradation;
5-hydroxymethylfurfural formation from 3-deoxyglucosone** and only in pumpkin seeds the conversion
of N-ε-fructoselysine to glyoxal and Heyns product to 1-deoxyglucosone **were found to be
quantitatively important**."*

§3.4.1, verbatim:

> "Sucrose, glucose, and FFC are the direct precursors of HMF in nuts and seeds during roasting.
> During heat treatment in dry systems, HMF can be formed directly from FFC or glucose via 3-DG.
> **It was determined that the 3-DG pathway was predominant in HMF formation for all samples
> according to reaction rate constants (Table 2).** The results were in agreement with the findings
> of Taş and Gökmen (2017) who stated that **the contribution of FFC to HMF formation is < 3-DG.**
> **However, they also emphasized that FFC was an important source in HMF formation, and the model
> was not suitable when this step was removed** due to the fact that FFC is an intermediate that
> cannot be detected experimentally and the dehydration steps in the formation of HMF from FFC are
> reduced to one step."

Conclusion, verbatim: *"the formation of HMF through 3-DG was found to be quantitatively important;
**sugar degradation was more predominant than MR in the formation of 3-DG**"*.

### 4.2 ★ WAVE K5a VERIFICATION — the claim IS supported by the table, unlike Gürsul Aktağ's

k₁₀ (3-DG→HMF) vs k₂ (FFC→HMF), directly comparable (both `min⁻¹×10³`, both terminal HMF-forming
steps in the same fitted model):

| matrix | k₁₀ 160 | k₂ 160 | ratio 160 | k₁₀ 180 | k₂ 180 | ratio 180 |
|---|---|---|---|---|---|---|
| Sunflower | 568.79 | **0** | **∞** | 109.98 | 89.49 | **1.2** |
| Flaxseed | 1024.28 | **0** | **∞** | 918.19 | **0** | **∞** |
| Peanut | 279.32 | 13.20 | **21.2** | 807.82 | 0.85 | **950** |
| Almond | 2528.66 | **0** | **∞** | 601.45 | 2.77 | **217** |
| Pumpkin | 37.53 | **0** | **∞** | 552.83 | **0** | **∞** |

**k₂ is exactly zero in 6 of 10 cells. k₁₀ > k₂ in all 10.** The narrowest margin is 1.2× (sunflower
at 180 °C). **The dissent is real and it is large.** `[D]`

### 4.3 ★★ BUT — THE THREE REASONS IT DOES **NOT** OVERTURN THE FRUCTOSE-LIMB FINDING

**(i) The comparison is not between equals.** `[3-DG]` is measured against an authentic standard;
`[FFC]` is never measured. In a fitted ODE only the *product* `k₂·[FFC]` is constrained. If the
model's implicit FFC pool is small, `k₂` inflates; if large, `k₂` collapses to zero — which is
exactly what it does in 6 of 10 cells. **A rate constant on an unmeasured node cannot be compared
in magnitude with one on a measured node. The authors say as much in the passage quoted in §4.1.**

**(ii) FFC is drained by two large competing sinks before it can dehydrate.** In this sucrose-rich
matrix `FFC` is consumed by `k₃` (FFC + bLYS → HP) and `k₂₀` (FFC + AA → P4), and both are
**orders of magnitude larger than k₂**: e.g. sunflower 160 °C, `k₃ = 96.83`, `k₂₀ = 317.56`,
`k₂ = 0`. So `k₂ ≈ 0` is a **flux-competition** result specific to a matrix carrying
7–12 g/kg of bound lysine — **not** evidence that the FFC → HMF elementary step is intrinsically
slow. The paper's own §3.4.2 confirms the drain: *"HP formation from FFC (k₃) appears to be more
important than FL formation from glucose (k₄) in samples indicating that **FFC was more important in
the early stages of MR**."*

**(iii) It is a sucrose system, not a hexose system.** The only entry to FFC here is
`SUC → FFC + GLC` (k₁). There is **no free-fructose node at all** — fructose is not in the network,
and glucose enters the 3-DG limb directly (`k₅`). The Kocadağlı systems, by contrast, generate a
large *free fructose* pool by isomerisation and route it through `Fru → Int → HMF`.
**⇒ Şen 2022 and Kocadağlı 2016 are not testing the same partition.** Şen tests
*sucrose-derived cation vs glucose-derived 3-DG*; Kocadağlı tests *isomerisation-derived fructose vs
glucose-derived 3-DG*.

**⇒ CORRECT FRAMING FOR THE ORCHESTRATOR: the branch fraction is not "fructose-limb wins" or
"3-DG-limb wins"; it is set by which precursor pool the matrix supplies and by what competes for
that pool. Şen 2022 bounds the 3-DG-dominant end of that range.** See the synthesis §5.

### 4.4 The 3-DG source is caramelization, not Maillard `[F]`

§3.4.3, verbatim: *"it was found that **the formation pathway of 3-DG through glucose dehydration
(k₅) was predominant in all samples**, while in peanut samples, the formation over the MR (k₈)
became predominant as high temperature was influential."*
Check `[D]`, k₅ vs k₈: sunflower 38.23 vs 0.57 (67×); flaxseed 34.94 vs 2.67 (13×); almond 544.68
vs 0 (∞); pumpkin 10.91 vs 0 (∞); peanut 160 °C 15.09 vs 0.5 (30×) but **peanut 180 °C 0 vs 2.23 —
the claimed reversal ✔**. ★ **Another same-lab, same-method demonstration that a branch source
flips with matrix and temperature.**

---

## §5. ★ WAVE K5a CROSS-CHECK — implied 2-point activation energies `[D]`

**THERE IS NO ARRHENIUS FIT AND NO ACTIVATION ENERGY IN THIS PAPER** (verified negative, §6). With
exactly two temperatures I computed the implied `Ea = −R·Δln k / Δ(1/T)` for the HMF branch, purely
as a transferability diagnostic:

| matrix | k₁₀(160) | k₁₀(180) | **implied Ea(k₁₀)** | k₂(160) | k₂(180) | implied Ea(k₂) |
|---|---|---|---|---|---|---|
| Sunflower | 568.79 | 109.98 | **−134.1 kJ/mol** | 0 | 89.49 | not computable |
| Flaxseed | 1024.28 | 918.19 | **−8.9** | 0 | 0 | not computable |
| Peanut | 279.32 | 807.82 | **+86.7** | 13.20 | 0.85 | −223.8 |
| Almond | 2528.66 | 601.45 | **−117.2** | 0 | 2.77 | not computable |
| Pumpkin | 37.53 | 552.83 | **+219.5** | 0 | 0 | not computable |

**★ The implied Ea for the SAME elementary step, by the SAME lab, with the SAME method, at the SAME
two temperatures, spans −134 to +220 kJ/mol across five matrices — and flips sign three times.**
`k₁₀` *decreases* with temperature in 3 of 5 matrices. **⇒ There is no transferable temperature
dependence for the 3-DG → HMF edge, and this is the cleanest demonstration of that fact anywhere in
the cluster** (five independent matrices, one lab, one method).

The same non-monotonicity infects the rest of the table: k₅ falls 38.23→5.10 (sunflower) and
544.68→2.53 (almond, a **215× drop**); k₁₂ falls in 4 of 5 matrices — the authors note it, verbatim
§3.4.3: *"The rate constant of DMG formation through 1-DG (k₁₂) **decreased when the roasting
temperature increased**."*

---

## §6. VERIFIED NEGATIVES `[NEG]`

- **No activation energy, no Arrhenius fit, no pre-exponential, no `k_b`, no Table 3.** Grepped
  the full text for `Arrhenius`, `activation`, `Ea`: **zero hits.** Each temperature is fitted
  separately.
- **No fructose node** — fructose is measured as a raw-material sugar (*"glucose and fructose were
  rapidly degraded upon roasting"*) but does not appear in either network.
- **No 3,4-dideoxyglucosone.** ⇒ **k₁₀ lumps `3-DG → 3,4-DG → HMF` into one step and is NOT
  commensurable with Kocadağlı's k₅/k₁₃ or Göncüoğlu Taş's k₁₉.**
- **No melanoidin / browning / colour measurement.**
- **No furfural.**
- **No lipid-oxidation node**, despite 42–49 % oil. The authors invoke lipid oxidation only
  qualitatively, for sunflower's GO: *"sunflower seeds have a matrix rich in polyunsaturated fatty
  acids which might be a source of GO due to lipid oxidation (Fujioka & Shibamoto, 2004)"* (§3.4.4).
- **Water is omitted from the model** — verbatim §3.4: *"Although the water content should be
  included in the model as it affects the reaction, **the heat and mass transfer coefficient of
  water were omitted** due to the restricted amount of water in nuts and seeds."*
- **HP was never measured** (no standard); **FL is furosine × 2.2**, not a direct measurement.
- Supplementary Tables S1–S5 and Figs. S1–S5 (the per-sample concentration time courses) are
  **not in the PDF on disk**.

---

## §7. DIRECTIONAL / STRUCTURAL CONSTRAINTS `STRUCTURAL`

| # | constraint | anchor |
|---|---|---|
| S1 | **In a sucrose-rich, reducing-sugar-poor, neutral-pH, low-a_w matrix the 3-DG limb carries HMF; the FFC limb is drained to Heyns product first** | §3.4.1, §3.4.2, §4.2–4.3 |
| S2 | **3-DG comes from glucose dehydration (caramelization), not from Amadori/Heyns degradation** — except peanut at 180 °C | §3.4.3, verified §4.4 |
| S3 | **k₃ > k₄ (HP from FFC beats FL from glucose) in every matrix except peanut@180 and almond** | §3.4.2 |
| S4 | **MGO comes only from 1-DG** here (the 3-DG route was rejected) — matching the wheat-flour paper, opposing the JAFC NaCl paper | §3.4.3 |
| S5 | **CEL > CML in every matrix**, because `k₁₁ > k₆, k₂₄` (MGO forms faster than GO) **and** `k₁₆ > k₁₄` (MGO reacts with lysine faster than GO) | §3.4.4, verbatim |
| S6 | **CML comes mostly from FL oxidation (k₁₃)**, except sunflower @180 °C where GO + bLYS (k₁₄) takes over — attributed to sunflower's PUFA-derived GO | §3.4.4 |
| S7 | **Furosine > CML in all samples ⇒ the early stage of the MR dominates** | §3.2 `[C]`+`[M]` |
| S8 | **Protein-bound lysine is the reactive amine; bound arginine and histidine do not change significantly (p > 0.05) in any sample** | §3.4.2 ★ a clean amine-selectivity constraint |
| S9 | **Glucosone is the least abundant dicarbonyl everywhere**, *"because glucosone is susceptible to oxidation and can be easily degraded due to its reductone structure"* | §3.2 |
| S10 | **Pumpkin seed needs a different network** (PCA-separated on a_w 0.62, low sucrose, low free AA): GO from FL, 1-DG from HP, no glucosone node | §3.3, §3.4 ★ **a measured demonstration that network TOPOLOGY, not just parameters, is matrix-dependent** |
| S11 | Amino acids decrease more slowly than sucrose **because AP/HP breakdown regenerates them** | §3.2 |
| S12 | CML and CEL **decrease** at long times / high temperature; the authors declined to model CEL degradation | §3.2, §3.4.4 |

---

## §8. USABILITY VERDICTS

| item | provenance | verdict |
|---|---|---|
| Table 2, all ~230 cells | `[F]` | **USE-Q** — `min⁻¹×10³` (or `kg·µmol⁻¹·min⁻¹×10³` for k17/k19/k20/k22), **per-matrix, per-temperature only**; drop the 12 `ind` and the ~40 HPD ≥ estimate cells (§3.1) |
| **k₁₀ (3-DG→HMF)** | `[F]` | **USE-Q** — well-determined in 6 of 10 cells; ⚠️ **lumps 3-DG→3,4-DG→HMF; not commensurable with Kocadağlı k₅/k₁₃ or Göncüoğlu k₁₉** |
| **k₂ (FFC→HMF)** | `[F]` | **RATIO-ONLY, and even the ratio is weak** — `[FFC]` unmeasured, zero in 6 of 10 cells |
| **k₁₀/k₂ ratios (∞, ∞, ∞, ∞, ∞, ∞, 1.2, 21.2, 217, 950)** | `[D]` | ★ **USE as a BOUND, not as a branch fraction.** Reads: *"in a sucrose-rich, low-moisture, reducing-sugar-poor, neutral-pH nut/seed matrix, the fitted terminal 3-DG→HMF constant exceeds the fitted terminal FFC→HMF constant by ≥1.2× and usually by ∞."* |
| **k₁₉ (HMF + AA → P3), the bimolecular sink** | `[F]` | **USE-Q with a conflict flag** — 6 of 10 cells are zero or HPD-dominated; **and the same lab's wheat-flour paper rejected this exact step.** The two well-determined cells are sunflower @160 (9.85 ± 5.59) and pumpkin @180 (43.66 ± 28.72), both `kg·µmol⁻¹·min⁻¹×10³` |
| **implied Ea range −134 to +220 kJ/mol for k₁₀** | `[D]` | ★ **STRUCTURAL — USE as the falsifier of any fixed 3-DG→HMF barrier.** Never as a value. |
| **any published activation energy** | — | ★ **DOES NOT EXIST IN THIS PAPER.** |
| Table 1 (pH, a_w, proximate) | `[M]` | ★ **USE — the cluster's best matrix envelope** |
| max HMF 247.0 ± 8.3 mg/kg (sunflower, 180 °C, 30 min) and 42.7 ± 0.4 (pumpkin) | `[M]` | ★ **USE — proposed HOLD-OUT** |
| CML / CEL maxima, furosine, sucrose losses, amino-acid losses | `[M]` | **USE** |
| `[FL] = furosine × 2.2` | `[C]` | **USE-Q — a borrowed stoichiometric factor, not measured here** |
| Eq. (4)/(27) sign on the `k₄` term | — | **printed error — read as `+`** |
| S1–S12 | `[M]`/`[F]` | **STRUCTURAL** — S1, S2, S8 and S10 are the highest-value |
