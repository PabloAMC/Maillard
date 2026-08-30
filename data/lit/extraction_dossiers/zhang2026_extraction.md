# Zhang et al. 2026 — glucose–glutathione multiresponse kinetics: 19 rate constants and 19 activation energies, audited

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★★ THE ONLY PAPER IN WAVE K6a THAT PUBLISHES RATE CONSTANTS AND ACTIVATION ENERGIES. **ALL 19 Eₐ REPRODUCE FROM THE AUTHORS' OWN k VALUES (mean |deviation| 2.1 %) — the Arrhenius arithmetic is sound. THE PROBLEM IS ELSEWHERE, AND IT IS WORSE: six of the constants are pinned to a node that was never measured, six more sit on a saturating top leg, and the whole table has no confidence intervals.**

---

## §0. IDENTITY `[M]` — **the file is what the brief expected**

| item | value | anchor |
|---|---|---|
| file | `data/articles/Zhang2026.pdf`, 6 061 008 B, 13 pages, Elsevier/Distiller, born-digital | `pdfinfo` |
| title | **"Early reaction pathways in a glucose-glutathione Maillard system to form meat flavor"** | p. 1 |
| authors | **Chenping Zhang**, Dawei Zhen, Jinxin He, Kexin Li, Yu Liang, Qingqing Hu, **Jianchun Xie\*** | p. 1 |
| affiliations | School of Food Science and Health, **Beijing Technology and Business University (BTBU)**; School of Chinese Materia Medica, Beijing University of Chinese Medicine | p. 1 |
| **DOI** | ★ **`10.1016/j.foodchem.2026.148681`** ✔ exactly as the brief expected | p. 1 |
| journal | ***Food Chemistry* 510 (2026) 148681** | running head |
| funding | **NSFC 32372462** | p. 12 |
| ethics | BTBU Ethical Permission **No. 234, 20 Nov 2024** (sensory/GC-O panel) | p. 12 |
| ★ unusual declaration | *"The authors declare in this article **no AI TOOLS are used at all**."* | p. 12 |
| competing interests | none | p. 12 |
| SI | **NOT ON DISK.** Tables **S1–S5** (¹H/¹³C NMR of the five prepared intermediates) and ★ **Tables S6–S8 — the volatile flavour compounds at 130, 140 and 150 °C**, which are the paper's only volatile data | pp. 5, 8 |

**Lineage note:** this is the Xie group at BTBU, methodologically downstream of **Kocadağlı & Gökmen 2016** and **Martins & Van Boekel 2003** (both cited in the Introduction as the multiresponse template). It is therefore a **sibling of the five Gökmen-school papers audited in wave K5a**, and the same audit classes apply.

---

## §1. SYSTEM AND CONDITIONS `[M]`

| parameter | value |
|---|---|
| charge | **GSH 0.3 mmol (0.092 g) + D-glucose 0.3 mmol (0.054 g)** in **5 mL** buffer → **60 mmol/L each** |
| ★ buffer | ★ **sodium dihydrogen phosphate, 0.2 mol/L**, adjusted to **pH 6.5** with 1 M NaOH — **BUFFERED, and at a pH the corpus rarely has** |
| vessel | **15 mL pressure glass vial**, HEL **Parallel Synthesis Poly-Block 4** |
| **temperatures** | ★ **130, 140, 150 °C** — three rungs, and the **highest** aqueous sulfur ladder in the K6a corpus |
| times | **10, 20, 30, 40, 50, 60, 90, 120 min** — ★ **eight time points per temperature** |
| control | **GSH alone (no glucose)**, same protocol — used to establish that Cys-Gly and pGlu come from GSH hydrolysis, not from the Maillard route |
| quench | ice-water |
| replication | **n = 3**, mean ± SD, one-way ANOVA + Duncan, SPSS 24 |
| ⚠ measured pH after heating | **NOT REPORTED** `[NEG]` — but the pot is buffered at 0.2 M, which is ~3× the substrate load, so drift should be modest |

### 1.1 The prepared intermediates
Cys-Amadori: cysteine 3.5 g + glucose 3.3 g + Na₂HPO₄ 7.2 g + 2.5 mL water, **50 °C / 120 min**.
Gly-Amadori: glycine 5 g + glucose 10 g + sodium bisulfite 8 g + 8 mL water, **90 °C / 20 min**.
Both purified on **Bio-Rad AG 50W-X4 (H⁺), 35 × 1.6 cm**, UV 220 nm. Purity by HPLC-ELSD.
⚠️ **No purity figure is printed** `[NEG]`.

---

## §2. QUANTIFICATION BASIS — ★ **the non-volatiles are genuinely absolute; the two model targets are not**

### 2.1 Non-volatiles: **ABSOLUTE, external standard, LC-MS/MS in MRM mode** — Table 1, transcribed in full

| compound | quantitative ions (m/z) | cone V | collision energy (eV) | standard range (mmol/L) | calibration curve | R² |
|---|---|---:|---:|---|---|---:|
| Gly | 76 | 65 | 10 | 0.1–15 | y = 13 535 x − 1023 | 0.9992 |
| Cys | 122 > 76 | 65 | 14 | 0.1–15 | y = 14 511 x − 989 | 0.9980 |
| pGlu | 130 > 84 | 105 | 21 | 0.1–30 | y = 69 353 x − 891 | 0.9931 |
| Cys-Gly | 179 > 161 | 75 | 6 | 0.1–30 | y = 87 837 x − 1242 | 0.999 |
| Gly-Amadori | 238 > 202 | 80 | 13 | 0.1–15 | y = 15 284 x − 912 | 0.9958 |
| Cys-Amadori | 284 > 122 | 110 | 13 | 0.1–15 | y = 18 583 x − 8821 | 0.9975 |
| **pGlu-Amadori** | 314 > 130 | 130 | 20 | ★ **–** | ★ **–** | – |
| GSH | 308 > 179 | 100 | 12 | 0.5–80 | y = 118 488 x + 1876 | 0.9939 |
| (Cys-Gly)-TTCA | 341 > 179 | 95 | 10 | 0.1–20 | y = 60 891 x + 1011 | 0.9967 |
| (Cys-Gly)-Amadori | 341 > 305 | 95 | 14 | 0.1–20 | y = 8834 x + 1021 | 0.9987 |
| **di-(Cys-Gly)** | 355 > 177 | 115 | 21 | ★ **–** | ★ **–** | – |
| GSH-Amadori | 470 > 386 | 115 | 18 | 0.1–20 | y = 31 573 x + 1103 | 0.998 |
| GSSG | 613 > 355 | 160 | 23 | 0.1–20 | y = 85 929 x + 1981 | 0.9982 |
| glucose | 179 > 81 | 145 | 18 | 0.5–80 | y = 9205 x − 698 | 0.9987 |

> *"pGlu-Amadori and dimer of Cys-Gly were not quantified by the external standard method as their
> standards were unavailable."* (p. 3)

**⇒ Twelve of fourteen non-volatile species are absolutely quantified with R² ≥ 0.993. This is the
best-instrumented non-volatile data set in the K6a corpus.** ⚠️ **But two model species —
pGlu-Amadori and di-(Cys-Gly) — are UNQUANTIFIED, and both carry rate constants
(k5, k8 and k17 respectively).**

⚠️ **A calibration-curve red flag:** Cys-Amadori's intercept is **−8821** against a slope of 18 583,
i.e. the curve crosses zero at **x = 0.475 mmol/L**, which is **above the lowest standard (0.1)**.
Every Cys-Amadori reading below ~0.5 mmol/L is an extrapolation into negative concentration.
Every other curve has an intercept ≤ 2 % of the slope; this one is **47 %**. `[!]`

### 2.2 ★★ THE TWO MODEL TARGETS THAT ARE NOT SPECIES — the paper's own words

> *"concentrations of **α-DC (α-dicarbonyl compounds) and melanoidins were respectively represented
> by 294 and 420 nm UV absorbance values**, and concentration of volatile flavor compounds was
> indicated by their **total amount at the reaction time of 120 min** by GC–MS analysis."* (p. 10)

and the authors' own limitation statement:

> *"the reaction pathways and kinetic modelling results are limited by the simplifying assumptions
> using 294 nm UV-absorbance for α-dicarbonyls and 420 nm UV-absorbance for melanoidins, and
> **rather than specific chemical entities**."* (p. 9)

★★ **THIS IS THE K5a §2.3 UNIDENTIFIABILITY CLASS, EXACTLY.** In a fitted ODE an unmeasured node's
concentration scale is not identified; only the product `k·[node]` is. Here **three nodes are on
arbitrary scales** — α-DC (absorbance), melanoidin (absorbance), and volatiles (a GC total in
ng/g). **Every rate constant that touches one of them — k4, k8, k11, k13, k14, k16 (production of
α-DC), k18 (→melanoidin) and k19 (→volatiles), plus the (k18+k19) sink terms that appear in the
GSH, Cys-Gly, Cys and Gly balances — is not commensurable in magnitude with a rate constant on a
measured node.** That is **8 of 19**, and they include the entire "first group" the authors
identify as having the largest k values.

### 2.3 Volatiles: **semi-quant, reported in ng/g, and only in the SI**
SPME (**DVB/CAR/PDMS, 2 cm, 50/30 μm**), 1 μL of **0.06 μg/μL 1,2-dichlorobenzene** IS, equilibrate
50 °C / 10 min, GC-MS + GC-O. 19 authentic chemicals listed for identification. **No calibration
curve for any volatile.** The only volatile numbers in the main text (p. 8, 150 °C / 120 min):
**3-methyl-2-thiophenecarboxaldehyde 48.10 ± 7.88; dihydro-2-methyl-3(2H)-thiophenone 22.41 ± 3.98;
2-acetylthiazole 24.41 ± 5.12 ng/g.**
**⇒ `RATIO-ONLY` on the volatile side; the ng/g values are IS-normalised areas.**

---

## §3. THE REACTION NETWORK (Scheme 3) — the fitted model

**19 rate constants** over **14 measured species + 3 unmeasured indices.** Verbatim rate laws, as
printed (equations 1–16, pp. 10–11; `α-DC` = UV₂₉₄, melanoidin = UV₄₂₀):

```
r(GSH)              = −k1[GSH][glucose] − k2[GSH] + k4[GSH-Amadori] − (k18+k19)[α-DC][GSH]
r(glucose)          = −k1[GSH][glucose] − k7[Cys-Gly][glucose] − k15[Cys][glucose] − k13[Gly][glucose]
r(GSH-Amadori)      =  k1[GSH][glucose] − (k4+k6)[GSH-Amadori]
r(Cys-Gly)          = −(k9+k17)[Cys-Gly] − k7[Cys-Gly][glucose] + k11[(Cys-Gly)-Amadori]
                       + k6[GSH-Amadori] − (k18+k19)[α-DC][Cys-Gly]
r((Cys-Gly)-TTCA)   =  k7[Cys-Gly][glucose] − k10[(Cys-Gly)-TTCA]
r((Cys-Gly)-Amadori)=  k10[(Cys-Gly)-TTCA] − (k11+k12)[(Cys-Gly)-Amadori]
r(Cys)              =  k9[Cys-Gly] − k15[Cys][glucose] + k16[Cys-Amadori] − (k18+k19)[α-DC][Cys]
r(Gly)              =  k9[Cys-Gly] + k12[(Cys-Gly)-Amadori] − k13[Gly][glucose] + k14[Gly-Amadori]
                       − (k18+k19)[α-DC][Gly]
r(Cys-Amadori)      =  k12[(Cys-Gly)-Amadori] + k15[Cys][glucose] − k16[Cys-Amadori]
r(Gly-Amadori)      =  k13[Gly][glucose] − k14[Gly-Amadori]
r(α-DC)             =  k16[Cys-Amadori] + k14[Gly-Amadori] + k8[pGlu-Amadori] + k4[GSH-Amadori]
                       + k11[(Cys-Gly)-Amadori] − (k18+k19)[α-DC]
r(melanoidin)       =  k18 × [α-DC] × [Gly] × [Cys] × [Cys-Gly] × [GSH]                    (14)
r(volatile cpds)    =  k19 × [α-DC] × [Gly] × [Cys]                                        (15)
r(di-(Cys-Gly))     =  k17 × [Cys-Gly]                                                     (16)
```

### 3.1 ★★ TWO STRUCTURAL DEFECTS IN THE RATE LAWS, BOTH SERIOUS

**(a) Equation 14 is FIFTH ORDER.** `r(melanoidin) = k18[α-DC][Gly][Cys][Cys-Gly][GSH]` — a product
of five concentrations. **No elementary or lumped step in Maillard chemistry is fifth order**, and
`k18` therefore carries units of **L⁴ mol⁻⁴ min⁻¹** (or absorbance-corrected equivalents). Equation
15 is third order for the same reason. **A fifth-order lump means k18's magnitude is a pure
artefact of the concentration scale**, and its Eₐ (118.0 kJ mol⁻¹) is the temperature dependence of
a five-way product, not of a barrier.

**(b) The α-DC sink is written two different ways in the same model.** In the α-DC balance the sink
is `−(k18+k19)[α-DC]` — **first order in α-DC alone**. In equations 14 and 15 the same two channels
are **fifth and third order**. **The mass balance and the product equations are mutually
inconsistent.** ⚠️ Either the α-DC balance is missing the co-reactant factors, or 14–15 have
extra ones. As printed, the model does not conserve.

⚠️ **(c) `k1–k22` are declared in the text; only k1–k19 appear in Table 2.** `k20, k21, k22` are
never defined and never reported. `[NEG]`

---

## §4. ★★ TABLE 2 — the complete kinetic table, transcribed `[M]` `[F]`

**Caption:** *"The relative RMSE and R² values in multiresponse kinetic fitting the experimental data
of glutathione-glucose system heated at 140 °C by Scheme 1, Scheme 2, and Scheme 3, respectively;
and the k values and activation energies calculated for glutathione-glucose system heated at 130,
140 and 150 °C, respectively, based on the modelling results by Scheme 3."*

### 4a. Model discrimination (all at 140 °C) — **the fit-quality half**

| response | Scheme 1 rRMSE (%) | Scheme 1 R² | Scheme 2 rRMSE (%) | Scheme 2 R² | ★ Scheme 3 rRMSE (%) | ★ Scheme 3 R² |
|---|---:|---:|---:|---:|---:|---:|
| glucose | 0.8335 | 0.4911 | 0.6930 | 0.5548 | **0.0277** | **0.9233** |
| GSH | 0.2214 | 0.7703 | 0.1911 | 0.8188 | **0.0090** | **0.9633** |
| GSH-Amadori | 0.3016 | 0.6927 | 0.2571 | 0.7331 | **0.0110** | **0.9521** |
| Cys-Gly | 0.5437 | 0.5925 | 0.2580 | 0.7243 | **0.0224** | **0.9367** |
| (Cys-Gly)-TTCA | 0.3245 | 0.6718 | 0.2772 | 0.7031 | **0.0114** | **0.9501** |
| (Cys-Gly)-Amadori | 0.4634 | 0.6136 | 0.3931 | 0.6446 | **0.0177** | **0.9440** |
| pGlu | 0.3333 | 0.6671 | 0.2853 | 0.7002 | **0.0125** | **0.9496** |
| Cys-Amadori | 0.3644 | 0.6518 | 0.3112 | 0.6841 | **0.0141** | **0.9441** |
| Gly-Amadori | 0.6935 | 0.5548 | 0.5830 | 0.5803 | **0.0246** | **0.9346** |
| pGlu-Amadori | 0.8942 | 0.4811 | 0.7320 | 0.5319 | **0.0354** | **0.9102** |
| Cys | 1.1045 | 0.3274 | 0.9021 | 0.4733 | **0.0366** | **0.9086** |
| Gly | 0.8514 | 0.4904 | 0.7022 | 0.5431 | **0.0312** | **0.9164** |
| di-(Cys-Gly) | 0.8252 | 0.5055 | 0.6852 | 0.5672 | **0.0297** | **0.9231** |
| GSSG | 0.9746 | 0.4656 | 0.7951 | 0.5115 | **0.0371** | **0.9066** |
| UV 294 nm | 0.5914 | 0.5745 | 0.5022 | 0.6007 | **0.0198** | **0.9437** |
| UV 420 nm | 0.9745 | 0.4656 | 0.7935 | 0.5146 | **0.0392** | **0.9041** |
| volatile compounds | 0.8926 | 0.4811 | 0.7331 | 0.5208 | **0.0338** | **0.9132** |

★ **The model discrimination is clean and unambiguous: Scheme 3 beats Schemes 1 and 2 on all
17 responses, by a factor of 20–30 on rRMSE.** That is a `STRUCTURAL` result the repo can use
independent of any rate constant.
⚠️ **But note what the comparison is:** Scheme 3 is the **union** of Schemes 1 and 2 (the authors say
so). A superset model beating its own submodels is nearly guaranteed and is **not** evidence that
Scheme 3 is right — only that neither hydrolysis route alone suffices. **The correct reading is
"both the GSH-hydrolysis route and the GSH-Amadori route are needed", which is exactly the claim
the abstract makes, and it is well supported.**

### 4b. ★ The rate constants and activation energies `[F]`

⚠️ **Table 2 prints NO UNITS for any k.** The rate laws are mixed-order (first, second, third and
fifth), so **the 19 constants do not share a unit.** From the Arrhenius section's `k₀ (s⁻¹)` the
implied time unit is seconds, but an 8-point time course over 120 minutes with k ≈ 0.09 would be
finished in ~1 minute at that reading. **The values are almost certainly min⁻¹-based.
`[!]` UNIT UNRESOLVED — do not ingest any k magnitude until it is settled.**

| # | step (authors' description) | **k(130 °C)** | **k(140 °C)** | **k(150 °C)** | **Eₐ published (kJ/mol)** | ★ **Eₐ refit `[D]`** | **R² `[D]`** | **dev.** | ★ leg 130→140 | ★ leg 140→150 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| k1 | GSH + glucose → GSH-Amadori | 0.0010 | 0.0029 | 0.0045 | 109.1 | **107.0** | 0.952 | −1.9 % | 147.4 | 63.9 |
| k2 | GSH → GSSG (oxidation) | 0.0001 | 0.0004 | 0.0008 | 138.9 | **147.8** | 0.969 | **+6.4 %** | 192.0 | 100.8 |
| k3 | GSH hydrolysis | 0.0004 | 0.0011 | 0.0022 | 123.5 | **121.1** | 0.991 | −2.0 % | 140.1 | 100.8 |
| **k4** | **GSH-Amadori → α-DC** | 0.0300 | 0.0643 | 0.0748 | 68.1 | **65.1** | ⚠ **0.879** | −4.4 % | 105.6 | ⚠ **22.0** |
| k5 | pGlu → pGlu-Amadori | 0.00002 | 0.00010 | 0.00026 | 184.3 | **182.2** | 0.983 | −1.1 % | 222.9 | 138.9 |
| k6 | GSH-Amadori → pGlu-Amadori + Cys-Gly | 0.0010 | 0.0019 | 0.0049 | 113.2 | **112.5** | 0.985 | −0.6 % | 88.9 | 137.7 |
| k7 | Cys-Gly + glucose → (Cys-Gly)-TTCA | 0.0017 | 0.0041 | 0.0077 | 105.4 | **107.3** | 0.993 | +1.8 % | 121.9 | 91.6 |
| **k8** | **pGlu-Amadori → α-DC** | 0.0431 | 0.0882 | 0.0911 | 53.5 | **53.5** | ⚠ **0.794** | −0.1 % | 99.2 | ⚠ **4.7** |
| k9 | Cys-Gly hydrolysis → Cys + Gly | 0.0002 | 0.0009 | 0.0014 | 132.8 | **138.6** | 0.917 | +4.4 % | 208.3 | 64.2 |
| k10 | (Cys-Gly)-TTCA → (Cys-Gly)-Amadori | 0.0017 | 0.0024 | 0.0074 | 110.6 | **103.8** | 0.906 | **−6.1 %** | ⚠ 47.8 | ⚠ **163.7** |
| **k11** | **(Cys-Gly)-Amadori → α-DC** | 0.0377 | 0.0734 | 0.0937 | 66.25 | **64.8** | 0.940 | −2.2 % | 92.3 | 35.5 |
| k12 | (Cys-Gly)-Amadori → Cys-Amadori + Gly | 0.0046 | 0.0083 | 0.0199 | 103.2 | **103.7** | 0.984 | +0.5 % | 81.7 | 127.1 |
| **k13** | **Gly + glucose → Gly-Amadori** | 0.0351 | 0.0529 | 0.0997 | 71.4 | **73.9** | 0.981 | +3.5 % | 56.8 | 92.1 |
| **k14** | **Gly-Amadori → α-DC** | 0.0453 | 0.0912 | 0.0959 | 53.4 | **53.6** | ⚠ **0.811** | +0.3 % | 96.9 | ⚠ **7.3** |
| k15 | Cys + glucose → Cys-Amadori | 0.0047 | 0.0096 | 0.0195 | 100.7 | **100.9** | **1.000** | +0.2 % | 98.9 | 103.0 |
| **k16** | **Cys-Amadori → α-DC** | 0.0263 | 0.0492 | 0.0881 | 87.9 | **85.7** | **1.000** | −2.5 % | 86.7 | 84.7 |
| k17 | Cys-Gly → di-(Cys-Gly) (oxidation) | 0.0008 | 0.0015 | 0.0045 | 122.9 | **122.2** | 0.971 | −0.6 % | 87.1 | 159.7 |
| **k18** | **α-DC → melanoidin** (5th order!) | 0.0010 | 0.0019 | 0.0053 | 118.0 | **118.0** | 0.979 | **0.0 %** | 88.9 | 149.1 |
| **k19** | **α-DC → volatile compounds** (3rd order) | 0.0023 | 0.0065 | 0.0098 | 104.9 | **103.1** | 0.947 | −1.7 % | 143.9 | 59.7 |

---

## §5. ★★ THE AUDIT — five findings, in decreasing order of consequence

### 5.1 ✅ **THE ARITHMETIC IS SOUND. All 19 Eₐ reproduce.**
Refitting `ln k` vs `1/T` on the authors' own three k values:
**mean |deviation| = 2.1 %; median 1.7 %; worst case k2 at +6.4 %; 12 of 19 within 2 %.**
**⇒ Unlike the Ma-2022 and Han-2025 defect classes and unlike four of the six Hamzalıoğlu
pre-exponentials (`hamzalioglu2018_extraction.md` §6), there is NO transcription or arithmetic
error in this Ea column.** The deviations are consistent with the authors having fitted to
unrounded k values and printed rounded ones.

### 5.2 ⚠️ **BUT SEVEN OF NINETEEN ARE BADLY CURVED, AND THE AUTHORS SAY SO THEMSELVES**
Three-point Arrhenius R² `[D]`: **7 of 19 fall below 0.95**, and three below 0.90:
**k8 (0.794), k14 (0.811), k4 (0.879)**, with k10 at 0.906 and k9 at 0.917.

The leg-by-leg breakdown shows why, and it is systematic:

| class | steps | leg 130→140 | leg 140→150 | ratio |
|---|---|---|---|---|
| ★ **Amadori → α-DC degradations** | k4, k8, k11, k14 | 92–106 | ⚠ **4.7 – 35.5** | **0.05 – 0.38** |
| **formation / hydrolysis / oxidation** | k6, k10, k12, k13, k17, k18 | 47–89 | **92–164** | **1.55 – 3.43** |
| clean | k15, k16 | 87–99 | 85–103 | **0.98 – 1.04** |

★ **The authors observed this and reported it in words** (p. 11):
> *"Regarding the reaction steps with smaller Eₐ, most of which involved Amadori compounds'
> degradation reaction …, **k values increased rapidly from 130 to 140 °C but slowly from 140 to
> 150 °C**, including k4, k8, k11, k13, k14, and k16. In contrast, k values kept increasing
> prominently from 140 to 150 °C, for the reaction steps with greater Eₐ."*

**⇒ They named the curvature and then published a single Arrhenius Eₐ for each step anyway.**
For k8 the published 53.5 kJ mol⁻¹ is the average of a **99.2** low leg and a **4.7** high leg —
a 21-fold spread. For k14 it is 96.9 and 7.3. **These are not barriers; they are averages across a
saturation.** ★★ **THIS IS THE MA-2022 PLATEAU-ARTEFACT CLASS that wave K5a §6.3 defined, and it is
the cleanest example of it in the corpus, because here the same table contains steps that do NOT
show it (k15, k16 are flat to within 4 %).**

**⇒ RULE FOR INGESTION: the four Amadori-degradation constants (k4, k8, k11, k14) must NOT be
ingested as barriers. Their low-leg values (92–106 kJ mol⁻¹) are the better estimate of the barrier;
their published values (53–68) are contaminated by a top-leg plateau.**

### 5.3 ⚠️ **THE SMALL CONSTANTS ARE ROUNDING-LIMITED, AND ONE IS UNINTERPRETABLE**
k values are printed to a fixed 4–5 decimals, so the smallest carry **one significant figure**:

| step | printed k(130) | implied rounding interval | ★ **Eₐ range from that one cell alone** | published |
|---|---|---|---|---|
| **k2** | **0.0001** | 0.00005 – 0.000149 | ★ **119.2 – 196.6 kJ mol⁻¹** | 138.9 |
| **k5** | **0.00002** | 0.000015 – 0.0000249 | ★ **166.4 – 202.3 kJ mol⁻¹** | 184.3 |
| k9 | 0.0002 | — | ± ~25 | 132.8 |
| k3, k6, k17, k18 | 0.0004–0.0010 | — | ± ~10 | — |

**⇒ k2's activation energy is uncertain by ±40 kJ mol⁻¹ from typesetting alone. `REFUSE` k2 and k5
as quantitative Eₐ.**

### 5.4 ⚠️ **NO CONFIDENCE INTERVALS, ANYWHERE — so identifiability is untested**
`[NEG]` **Table 2 prints no standard error, no HPD interval, no confidence bound on any of the 19
constants.** The Gökmen-school papers audited in K5a all printed at least an approximate interval,
which is what let that wave detect the **Han-2025 class (HPD ≥ estimate)**. **Here that test cannot
be run at all.** With 19 free constants, 8 time points and 3 temperatures, and **three of the
17 responses being arbitrary-scale indices**, the a-priori expectation is that several constants are
poorly identified — and there is no way to find out from the paper.

### 5.5 ★★ **THE DEEPEST PROBLEM: 8 OF 19 CONSTANTS ARE PINNED TO AN UNMEASURED NODE**
Per §2.2, α-DC and melanoidin are **UV absorbances**, and volatiles are a **GC total**. Every
constant producing or consuming those — **k4, k8, k11, k13(→Gly-Amadori→α-DC), k14, k16, k18, k19**
— has a magnitude that depends on the arbitrary absorbance scale.
★ **The Eₐ SURVIVE this** (a constant scale factor cancels in `ln k₂/k₁`), **but the k magnitudes
do not, and neither does any comparison between a k on a measured node and a k on an absorbance
node.** The authors' central magnitude claim — that "the first group with the greatest k values"
consists of **k4, k8, k11, k13, k14, k16**, i.e. **five of six are α-DC-producing steps on the
absorbance scale** — is therefore **not established**: those constants are large partly because
their product node is measured in absorbance units of order 0.1–1 rather than mmol/L of order 10.

**⇒ ★ PROPOSED CORRECTION FOR THE CORPUS:** the sentence *"most of [the highest rate constants]
were the degradation of Amadori compounds into α-dicarbonyls"* is a **scale artefact plus a
plateau artefact**, and should not be carried forward as a chemistry claim. **What survives is the
Eₐ ordering, not the k ordering** — the same distinction wave K5a §2.3 drew for the fructose limb.

---

## §6. WHICH STEPS OVERLAP THE REPO'S NETWORK

| Zhang step | repo counterpart | ingestible? |
|---|---|---|
| **k1** GSH + glucose → GSH-Amadori, **Eₐ 107.0** | the trunk's amine + sugar → Amadori edge | ★ **`USE-Q`** — measured node both sides, R² 0.952, no plateau. **The single most transferable number in this paper.** |
| **k15** Cys + glucose → Cys-Amadori, **Eₐ 100.9, R² 1.000** | ★ **the sulfur lane's own Amadori-formation edge** | ★★ **`USE-Q` — the best-conditioned constant in the table** (both legs 98.9 / 103.0, a 4 % spread) |
| **k16** Cys-Amadori → α-DC, **Eₐ 85.7, R² 1.000** | Amadori degradation on the sulfur limb | ★ **`USE-Q`** — flat legs (86.7 / 84.7) despite the absorbance node, because Eₐ is scale-free |
| **k13** Gly + glucose → Gly-Amadori, **Eₐ 73.9** | the glycine/generic amine limb | `USE-Q`, mild curvature |
| k4, k8, k11, k14 Amadori → α-DC | Amadori degradation | ⚠ **`PRIOR-ONLY` at the LOW-LEG value (92–106), never at the published 53–68** |
| **k3** GSH hydrolysis **123.5**, **k9** Cys-Gly hydrolysis **138.6** | ★ **nothing in the repo** — peptide hydrolysis is not modelled | ★ `STRUCTURAL` + `PRIOR-ONLY`: **a peptide-bound cysteine releases free Cys with a ~125–140 kJ mol⁻¹ barrier**, which is far above the 54.7 kJ mol⁻¹ of free-Cys *degradation* (`zhai2023foodchem_extraction.md` §9.4). ⇒ **in a peptide system the rate-limiting step for H₂S supply is hydrolysis, not Cys degradation** |
| **k2** GSH → GSSG, **k17** Cys-Gly → dimer | ★ **thiol oxidation to disulfide — the repo's `zhou_dimer_share` channel** | ⚠ k2 `REFUSE` (rounding); **k17 = 122.2 kJ mol⁻¹, R² 0.971, `USE-Q`** — ★ **the corpus's only measured Eₐ for thiol→disulfide oxidation** |
| **k19** α-DC → volatiles, **Eₐ 103.1** | the volatile-formation lump | `PRIOR-ONLY` — 3rd-order lump on two arbitrary scales |
| **k18** α-DC → melanoidin, **Eₐ 118.0** | the browning lump | ⚠ `REFUSE` as a barrier — **5th-order rate law** (§3.1a). ★ **Note it is 118 vs the ~28–31 kJ mol⁻¹ that two independent A₄₂₀ ladders give for browning** (`zhai2023jafc` §5.1, `wang2026` §3.2) — **a 4× disagreement that is fully explained by the 5th-order form** |

★★ **THE HEADLINE CROSS-CHECK: k18's Eₐ (118.0) is four times the directly measured browning Eₐ
(28.3 and 30.6 kJ mol⁻¹ from two other labs).** A fifth-order rate law absorbs the temperature
dependence of four co-reactant concentrations into its own "barrier". **This is a worked example of
why a lumped multi-order sink's Eₐ is not a barrier — and it is directly relevant to the repo's
`Ea_decay_thiol_sink` landing at 248 kJ mol⁻¹ against a 250 ceiling** (`kinetic_core_b2_2_diagnosis.md`
F-5): **a lumped sink whose true rate law has co-reactants will always inflate its fitted Eₐ.**

---

## §7. THE Kang SWITCH-ON CROSS-CHECK

**Not testable here in the usual sense**: this paper's ladder is **130/140/150 °C**, entirely *above*
the 120–140 °C window, and it reports **no per-compound volatile ladder in the main text** (Tables
S6–S8 are the SI and are not on disk). What it does supply is the **kinetic-constant analogue**:

★ **Every one of the four Amadori→α-DC steps SATURATES on the 140→150 °C leg** (Eₐ falls to
4.7–35.5 kJ mol⁻¹), while every formation/hydrolysis/oxidation step **accelerates** on the same leg
(Eₐ rises to 92–164). **⇒ Above 140 °C the precursor-supply chemistry is already saturating while
the consumption chemistry is still accelerating.**

**⇒ This is the mechanism that would produce a peak-and-decline in free thiols above ~140 °C, and
it is exactly what Zhai's Fig. 3 measures directly** (FFT at 140 °C peaks at 80 min and falls 2.5×,
`zhai2023jafc_extraction.md` §8.5) **and what Wang measures over 105→125 °C in a different pot**
(`wang2026_extraction.md` §4). **Three independent lines, three labs, one conclusion: above roughly
120–140 °C, thiol consumption outruns thiol formation.** That is the opposite of a switch-ON.

⚠️ **Honest limit:** Zhang's constants are for a **glucose–glutathione** system at **pH 6.5,
buffered**, and its α-DC node is an absorbance. The agreement is at the level of *shape*, not of
number.

---

## §8. OTHER MEASURED CONTENT `[M]`

| finding | value | anchor |
|---|---|---|
| volatiles at 130 °C / 120 min | **19 compounds**: 15 sulfur, 3 oxygen heterocycles, **1 nitrogen heterocycle** | p. 8 |
| effect of temperature on the profile | *"the kinds of volatile flavor compounds almost fixed whereas their amounts increased"* | p. 8 |
| top volatiles at **150 °C / 120 min** | 3-methyl-2-thiophenecarboxaldehyde **48.10 ± 7.88**; dihydro-2-methyl-3(2H)-thiophenone **22.41 ± 3.98**; 2-acetylthiazole **24.41 ± 5.12** ng/g | p. 8 |
| ★ **MFT** | *"2-methyl-3-furanthiol (meaty) … exhibited **higher FD values by GC-O** under each of the three reaction temperatures at 120 min"* — **FD values themselves are in Tables S6–S8, not on disk** | p. 8 |
| UV₂₉₄ and UV₄₂₀ | both increase monotonically with time and with temperature | Fig. 3 |
| Cys, Gly, pGlu | *"always increased over reaction times within 120 min for the three reaction temperatures"* — **no turnover in the peptide hydrolysates** | p. 8 |
| di-(Cys-Gly), GSSG | *"always ascended with reaction times"* — **oxidation is monotone** | p. 8 |
| fructose | ★ **excluded from the network because *"fructose was almost absent in the reaction solutions during the LC-MS analysis"*** | p. 9 |

★★ **THE FRUCTOSE EXCLUSION IS A DIRECT, MEASURED CONTRADICTION OF THE K5a HMF CLUSTER.**
Wave K5a §2 established, from three Gökmen matrices, that **deleting the fructose/cation → HMF edge
makes the model badly under-predict**, and treated the fructose limb as a hard structural
constraint. **Zhang's LC-MS finds essentially no fructose in a buffered pH-6.5 glucose pot at
130–150 °C and therefore omits the isomerisation entirely — and still fits 17 responses at
R² 0.90–0.96.** The authors flag their own choice: *"The disregardance of fructose might induce
potential deviations in the calculated kinetic parameters."*
**⇒ `[!]` Recorded as an open contradiction for the orchestrator: the fructose limb may be
pH- and matrix-conditional rather than universal. K5a's evidence is at 160–200 °C in low-moisture
and melt systems; this is 130–150 °C in 0.2 M phosphate at pH 6.5.** Not adjudicated here.

---

## §9. VERIFIED NEGATIVES `[NEG]`

| # | sought | verdict |
|---|---|---|
| 1 | Confidence intervals / SEs on any rate constant | ★ **ABSENT** (§5.4) |
| 2 | Units on any rate constant | ★ **ABSENT** (§4b) |
| 3 | Pre-exponential factors (k₀) | **ABSENT** — the Arrhenius form is printed with k₀ but no value is reported for any step |
| 4 | k20, k21, k22 | **DECLARED IN TEXT, NEVER REPORTED** (§3.1c) |
| 5 | α-dicarbonyl identification or quantification as chemical species | ★ **ABSENT — UV₂₉₄ only** |
| 6 | Melanoidin characterisation | **ABSENT — UV₄₂₀ only** |
| 7 | Per-compound volatile ladder in the main text | **ABSENT — Tables S6–S8, not on disk** |
| 8 | Measured pH after heating | **ABSENT** |
| 9 | MFT, FFT, or any thiol concentration | ★ **ABSENT from the main text** — MFT appears only as a GC-O FD-value claim |
| 10 | H₂S | **ABSENT** |
| 11 | Any temperature below 130 °C | **ABSENT** |
| 12 | Sensory scores | **ABSENT from the main text** despite the ethics approval for a panel |
| 13 | Purity of the prepared intermediates | **ABSENT** |

---

## §10. CONSOLIDATED PARAMETER TABLE

**Common conditions:** GSH 60 mmol/L + glucose 60 mmol/L in **0.2 M NaH₂PO₄ buffer, pH 6.5**,
15 mL pressure vial, 130/140/150 °C, 8 time points to 120 min, n = 3.
**Rate-constant units UNRESOLVED (§4b); Eₐ in kJ mol⁻¹.**

| # | parameter | value | condition | prov | verdict |
|---:|---|---|---|---|---|
| 1 | ★ **Eₐ, Cys + glucose → Cys-Amadori (k15)** | ★ **100.9** (published 100.7), R² **1.000**, legs 98.9/103.0 | 130–150 °C, pH 6.5 buffered | **F/D** | ★ `USE-Q` |
| 2 | ★ **Eₐ, Cys-Amadori → α-DC (k16)** | ★ **85.7** (published 87.9), R² **1.000**, legs 86.7/84.7 | " | **F/D** | ★ `USE-Q` |
| 3 | **Eₐ, GSH + glucose → GSH-Amadori (k1)** | **107.0** (published 109.1), R² 0.952 | " | **F/D** | `USE-Q` |
| 4 | **Eₐ, Gly + glucose → Gly-Amadori (k13)** | **73.9** (published 71.4), R² 0.981 | " | **F/D** | `USE-Q` |
| 5 | ★ **Eₐ, Cys-Gly → di-(Cys-Gly) (thiol oxidation, k17)** | ★ **122.2** (published 122.9), R² 0.971 | " | **F/D** | ★ `USE-Q` — **the corpus's only thiol→disulfide Eₐ** |
| 6 | **Eₐ, GSH hydrolysis (k3)** | **121.1** (published 123.5), R² 0.991 | " | **F/D** | `USE-Q` |
| 7 | **Eₐ, Cys-Gly hydrolysis (k9)** | **138.6** (published 132.8), R² 0.917 | " | **F/D** | `PRIOR-ONLY` |
| 8 | **Eₐ, (Cys-Gly)-TTCA formation (k7)** | 107.3 (published 105.4), R² 0.993 | " | **F/D** | `USE-Q` |
| 9 | **Eₐ, α-DC → volatiles (k19)** | 103.1 (published 104.9) | " | **F/D** | `PRIOR-ONLY` (3rd-order lump) |
| 10 | ⚠ **Eₐ, α-DC → melanoidin (k18)** | 118.0 | " | **F** | ★ `REFUSE` — 5th-order rate law; 4× the measured browning Eₐ |
| 11 | ⚠ **Eₐ, Amadori → α-DC (k4/k8/k11/k14)** | published **68.1 / 53.5 / 66.25 / 53.4**; ★ **low-leg values 105.6 / 99.2 / 92.3 / 96.9** | " | **F/D** | ★ `PRIOR-ONLY` **at the low-leg value only** (§5.2) |
| 12 | ⚠ **Eₐ, k2 (GSH→GSSG) and k5** | 138.9 and 184.3 | " | **F** | ★ `REFUSE` — rounding-limited to ±40 and ±18 (§5.3) |
| 13 | **All 57 k values** | see §4b | " | **F** | ⚠ `REFUSE` **until units are established** |
| 14 | ★ **Model discrimination** | Scheme 3 ≫ Schemes 1, 2 on all 17 responses (rRMSE 20–30× lower) | 140 °C | **F** | ★ `STRUCTURAL` |
| 15 | ★ **Fructose is absent** in a buffered pH-6.5 glucose pot at 130–150 °C | qualitative, LC-MS | " | **M** | ★ `STRUCTURAL` ⚠ contradicts K5a |
| 16 | Table 1 calibration curves (14 species) | see §2.1 | — | **M** | `USE` |
| 17 | Volatiles at 150 °C/120 min | 48.10 ± 7.88 / 22.41 ± 3.98 / 24.41 ± 5.12 ng/g | — | **M** | `RATIO-ONLY` |
| 18 | Compound count | 19 volatiles at 130 °C (15 S, 3 O, 1 N) | — | **M** | `STRUCTURAL` |

---

## §11. USABILITY VERDICTS — summary

- **`USE`** — Table 1's calibration curves; the model-discrimination result.
- **`USE-Q`** — the six well-conditioned Eₐ: **k15 (100.9), k16 (85.7), k1 (107.0), k17 (122.2),
  k3 (121.1), k7 (107.3)**, each qualified by: no confidence interval published, buffered pH 6.5,
  a glutathione (peptide) system rather than free cysteine, and a 130–150 °C range that is **above**
  most of the repo's cooking window so any use below 130 °C is an extrapolation.
- **`PRIOR-ONLY`** — k4/k8/k11/k14 at their **low-leg** values; k9; k19.
- **`RATIO-ONLY`** — all volatile ng/g values.
- **`STRUCTURAL`** — Scheme 3 vs 1/2; the fructose absence; "Amadori→α-DC saturates above 140 °C
  while formation/oxidation steps accelerate"; "peptide hydrolysis (121–139 kJ/mol) is a higher
  barrier than free-Cys degradation (55 kJ/mol)".
- **★ `REFUSE`** — every k **magnitude** (units unresolved, and 8 of 19 sit on arbitrary-scale
  nodes); **k18** as a barrier (5th order); **k2 and k5** as quantitative Eₐ (rounding);
  the published **53–68 kJ mol⁻¹** for the Amadori-degradation steps (plateau artefact);
  and the authors' claim that Amadori degradation has the largest rate constants (scale artefact).

---

## §12. DRAFT FIT / HOLD-OUT ROLES — **FOR THE ORCHESTRATOR, NOT A DECLARATION EDIT**

### 12a. Recommended **FIT** candidates
| candidate | why |
|---|---|
| ★ **Eₐ(Cys + glucose → Cys-Amadori) = 100.9 kJ mol⁻¹** | R² 1.000, both legs within 4 %, measured node on both sides, directly on the sulfur limb. **The strongest single kinetic number this wave found.** |
| ★ **Eₐ(Cys-Amadori → α-DC) = 85.7 kJ mol⁻¹** | R² 1.000, flat legs; scale-free despite the absorbance node |
| **Eₐ(thiol → disulfide) = 122.2 kJ mol⁻¹** | the repo's dimer channel has no measured barrier at all |

### 12b. Recommended **HOLD-OUT**
| candidate | what it tests |
|---|---|
| ★ **"Amadori→α-DC steps saturate on the 140→150 °C leg while formation steps accelerate"** — as a sign test on the model's own effective barriers | whether the core reproduces the *crossover* in temperature sensitivity between supply and consumption |
| **"Peptide hydrolysis is slower than free-Cys degradation by ~70 kJ mol⁻¹ of barrier"** | the repo has no peptide route; this bounds how a protein-hydrolysate matrix would differ |
| **"Fructose is undetectable at pH 6.5, 130–150 °C"** | ⚠ **conflicts with K5a's structural constraint; adjudicate before scoring** |

### 12c. **DO NOT USE**
Every k magnitude; k18 as a barrier; the published Amadori-degradation Eₐ; k2 and k5.

---

## §13. DECLARED GAPS

| # | gap | closure |
|---|---|---|
| G1 | ★ **Tables S6–S8 — the volatile compounds and GC-O FD values at 130/140/150 °C** — are the paper's only per-compound volatile data and are not on disk. **They would give a fourth independent MFT ladder, at 130–150 °C, in a peptide system.** | order the SI of `10.1016/j.foodchem.2026.148681` — ★ **second-highest-value SI request of this wave, after Wang's Table S2** |
| G2 | ★ **The units of the 19 rate constants** | ask the authors, or re-derive by reconstructing the ODE from Fig. 3's digitised curves (large job) |
| G3 | **No confidence intervals** — identifiability untestable | not recoverable |
| G4 | **α-DC and melanoidin are absorbances, not species** | not recoverable from this paper; the Gökmen-school papers do measure 3-DG and are the better source for that node |
| G5 | ⚠ **The α-DC mass balance is inconsistent with equations 14–15** (§3.1b) | would need the authors' MATLAB code |
| G6 | **k20–k22 undefined** | — |
| G7 | ⚠ **The fructose contradiction with K5a** | needs adjudication at orchestrator level |
