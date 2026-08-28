# Wave K4b — paired matrix/water thresholds, air-matrix partition, and amorphous browning kinetics
### Eight-paper synthesis, 2026-08-28. Read-only wave: **no repo file outside `data/lit/extraction_dossiers/` was written, staged or modified. Nothing was committed.**

Eight per-paper dossiers accompany this file:
`hong2020_`, `Meynier2002_`, `leksrisompong2010_`, `baek1999_`, `lievonen2002_`,
`miao2004_`, `pereyragonzales2010_`, `anantharamkrishnan2020_extraction.md`.
**All eight PDFs are the expected papers. No mis-file.** Two naming caveats only:
`Meynier2002.pdf` is published as ***Lait* 83 (2003) 223–235** (2002 = acceptance year), and
`anantharamkrishnan2020.pdf` is a **JAFC ASAP article with no volume or page numbers assigned**.

Builds on and extends `k2_matrix_and_thresholds.md`. **Its central conclusion — no general
matrix correction factor — is confirmed, not overturned.** Two of its subsidiary claims are
now sharpened by NEW measured data and one is flagged as a partial contradiction (§C).

---

# THE FIVE-LINE ANSWER

1. **The corpus now has TWO same-panel, same-method paired matrix/water threshold sets** —
   Hong 2020 (10 compounds, water vs soybean paste) and Leksrisompong 2010 (3 compounds,
   8 matrices, water vs caseinate vs 10 %/20 % emulsion vs oil). **Every ratio in
   `k2_matrix_and_thresholds.md` §A.8 was `cross_study_cross_method`. These are not.** §A.
2. **The K2 verdict survives.** Hong's same-panel ratios still span **29× to 2 035× (70×
   max/min)** and **one of ten inverts**; Leksrisompong's span **0.54× to 23.5×** with
   **three of twenty-one below 1**. A single correction factor is still refuted — but the
   **1-σ band collapses from K2's 27–41× to 4.7×**, so most of K2's *dispersion* was method
   noise while its *tails* were real. §A.2.
3. **Hydrophobicity does not predict the matrix shift, and this is now measured three ways.**
   Hong: **r = −0.05** between log P and the log matrix shift over ten compounds. Leksrisompong:
   the partition-based prediction of the oil/water threshold shift errs by **30.6×, 7.2× and
   2.4×**, and **gets the sign wrong for two of three compounds**. Meynier: a two-phase
   partition model of a real emulsion **fails in both directions on two structural isomers**.
   **`k2_matrix_and_thresholds.md` §D.4.5 ("stop deriving matrix odour activity from `f_free`")
   is now demonstrated on matched samples rather than inferred.** §B.
4. **A new, quantified, previously unrecorded contributor to the hexanal over-prediction:
   the air/water partition constant itself is uncertain by ~10×.** Meynier measures
   **K_aw(hexanal) = 2.5 × 10⁻³ at 30 °C**, which is **4.6× BELOW the EPI-Suite Henry value
   printed in its own Table I** and **7.2× below Hall & Anderson's**. Two independent
   static-headspace papers in this wave (Meynier, Leksrisompong) land **6–17× below** their
   own cited reference values, for the same methodological reason. §D.
5. **The browning-kinetics trio delivers one fully-tabulated real-food dataset and two
   figure-bound model-system datasets — and they contradict each other on whether Ea depends
   on water.** Pereyra Gonzales (real skim milk powder, 18 k + 5 Ea, all tabulated, all
   verified) says **no**; Miao & Roos (lactose/trehalose model, 4 Ea recovered from a figure)
   says **yes, −12 to −18 kJ/mol per g water/100 g solids**. §E proposes the reconciliation:
   **scope `Ea_depends_on_water` to `T − Tg < ~20 °C`.**

---

# (A) THE PAIRED THRESHOLD TABLE — the corpus's first same-method matrix/water ratios

## A.1 Hong & Kim 2020 — water vs soybean paste, one panel, n = 20, 3-AFC ASTM E679-91

Matrix A: **distilled water, 12 mL, dosed and served immediately.**
Matrix B: **soaked 20 h + autoclaved 121 °C/40 min + ground 24 000 rpm paste, 12 g, stomached
30 s, vacuum-sealed, 24 h at 4 °C.** Matrix-matched blanks. Group BET = geometric mean of
individual BETs. **`cross_study_cross_method: FALSE`.**

| compound | **water BET µg/kg** | **soybean BET µg/kg** | **ratio** | Δlog₁₀ | log P |
|---|---:|---:|---:|---:|---:|
| butyric acid | 0.9 | 1,831.2 | **2,035×** | +3.31 | 1.00 |
| dimethyl disulfide | 2.7 | 4,737.9 | **1,755×** | +3.24 | 1.22 |
| 2-pentyl furan | 8.6 | 11,690.2 | **1,359×** | +3.13 | 3.70 |
| 3-methyl butanal | 0.5 | 131.6 | **263×** | +2.42 | 1.79 |
| 2-methyl butanal | 2.1 | 548.9 | **261×** | +2.42 | 1.76 |
| 4-vinyl phenol | 166.0 | 25,139.9 | **151×** | +2.18 | 2.40 |
| **hexanal** | **12.6** | **1,669.6** | **132.5×** | **+2.12** | 1.77 |
| 2,3-dimethyl pyrazine | 401.0 | 19,600.8 | **48.9×** | +1.69 | 0.54 |
| 4-ethyl phenol | 23.2 | 675.1 | **29.1×** | +1.46 | 2.58 |
| **ethyl-4-methylpentanoate** | **9.7** | **1.5** | **0.155× (6.5× EASIER)** | **−0.81** | 2.10 |

⚠️ The printed **Z column is not a test statistic** — every entry equals `BET2 − BET1` to
within 1.4 %. ⚠️ The **± column's log base is not stated**. Both detailed in
`hong2020_extraction.md` §4.

## A.2 The dispersion, against K2's cross-study estimates

| statistic | **Hong 2020 (same panel), 9 elevated** | K2 beef/water | K2 gelatin/water |
|---|---:|---:|---:|
| geometric mean | **277×** | 906× | 33× |
| max / min | **70×** | 87× | 269× |
| **log₁₀ SD** | **0.668 dec** | 0.71 dec | 0.80 dec |
| **1-σ band width** | **4.7×** | **27×** | **41×** |

**⇒ The band narrows 6–9× when the method confound is removed, and the max/min spread does
not. `k2_matrix_and_thresholds.md` §D.4.1 (ship a lookup table, not a formula) stands.**
Including the inverted ester, the 10-compound log₁₀ SD is **1.206 decades (16× band)**.

## A.3 Leksrisompong 2010 — one panel, n = 25, 8 matrices, 3 compounds, thresholds AND partition coefficients on the SAME samples

BET (ppb, 22 °C) with SE; K_HS/matrix (dimensionless, 40 °C) with SD. Full tables in
`leksrisompong2010_extraction.md` §3, §4.

| matrix | diacetyl BET | δ-decalactone BET | furaneol BET | diacetyl K | δ-dec. K | furaneol K |
|---|---:|---:|---:|---:|---:|---:|
| **Water** | **6.0** | **66.0** | **22.3** | 2.60e-5 | 1.70e-5 | 3.1e-7 |
| 0 % fat pH 7.0 | 40.8 | 43.5 | 90.8 | 1.50e-5 | 2.20e-5 | 8.2e-7 |
| 0 % fat pH 5.5 | 44.9 | 35.8 | 46.4 | 2.70e-5 | 1.60e-5 | 2.2e-7 |
| 10 % fat pH 7.0 | 5.6 | 546 | 56.4 | 1.50e-5 | 2.3e-6 | 2.5e-7 |
| 10 % fat pH 5.5 | 9.2 | 294 | 67.8 | 1.60e-5 | 2.6e-6 | 9.0e-7 |
| 20 % fat pH 7.0 | 21.8 | 113 | 66.9 | 1.70e-5 | 1.5e-6 | 2.9e-7 |
| 20 % fat pH 5.5 | 8.6 | 329 | 134 | 1.60e-5 | 1.0e-6 | 7.2e-7 |
| **Oil** | **99.5** | **1,550** | **27.4** | 4.80e-5 | 1.0e-7 | 6.0e-7 |
| **oil/water BET ratio** | **16.6×** | **23.5×** | **1.23×** | — | — | — |

**⚠️ pH (5.5 vs 7.0) affects the PARTITION coefficient in 4 of 12 comparisons and the
THRESHOLD in 0 of 12.** A statistically detectable change in headspace concentration produced
no detectable change in perception — a third independent confirmation of
`k2_matrix_and_thresholds.md` §B.5.

## A.4 The hexanal ladder, complete

| matrix | protein | hexanal matrix threshold | aqueous ref. | **ratio** | pairing quality |
|---|---|---:|---:|---:|---|
| 3 % gelatin, 22 °C (Vega 1994) | 30 g/L collagen | 58 ppb | 4.5 ppb (Guadagni) | **12.9×** | cross-study, cross-method |
| **soybean paste, ~22 °C (Hong 2020)** | **not stated** | **1,669.6 µg/kg** | **12.6 µg/kg, same panel** | **132.5×** | **SAME STUDY ✅** |
| cooked beef 1:1, 45 °C (Brewer 1995) | ~100 g/L | 5.87 ppm | 4.5 ppb (Guadagni) | **1,304×** | cross-study + `dose_added_pre_cook` |

**Hong's 132.5× lands exactly where `k2_matrix_and_thresholds.md` §D.2(ii) predicted a
properly-controlled number would fall once Brewer's pre-cook thermal loss was removed.
That is independent corroboration of K2 owner-decision item #2.**

---

# (B) THE MECHANISM VERDICT — partition and binding are refuted THREE independent ways

| evidence | result | source |
|---|---|---|
| **corr(log P, log matrix shift), 10 compounds, same panel** | **r = −0.052**; slope **−0.070 dec per log P unit**. Largest shift is a log P 1.00 acid (2 035×); smallest elevated is a log P 2.58 phenol (29.1×) | Hong 2020 §6.2 |
| **corr(log K_HS/matrix, log BET) over 8 matched matrices** | **δ-decalactone (log P 3.4) r = −0.907 ✅ correct sign; diacetyl (log P −2.0) r = +0.600 ❌; furaneol (log P 1.4) r = +0.427 ❌** | Leksrisompong 2010 §6.1 |
| **partition-model prediction of the oil/water threshold shift** | errs by **30.6× (diacetyl, under), 7.2× (δ-decalactone, over), 2.4× (furaneol, under)**; over all 21 matrix cells the error spans **0.14× to 30.6×** | Leksrisompong 2010 §6.2 |
| **two-phase emulsion partition model vs measurement** | fails **in both directions on two structural isomers** (amyl acetate over-releases, isoamyl acetate under-releases; both log P ≈ 2.3); ethyl pentanoate alone is N.S. across 30–80 °C | Meynier 2002 §8.1 |
| **max headspace suppression by protein, four independent systems** | **1.07–1.39× (34 g/L dairy protein, Meynier); 0.79–1.69× (10–30 g/L caseinate, Leksrisompong — note >1 means ENHANCEMENT); 1.25–3.7× (40 g/L β-lg, Andriot via K2 §B.4)** | this wave + K2 |
| **perceived intensity vs in-nose concentration, simultaneous, within-subject** | in-nose Imax falls **1.5×, NOT significant**; perceived TI Imax falls **3.0×, P < 0.001**, over the same mouthfuls of the same 11 panellists | Baek 1999 §5.1 |
| **the variable that DOES track perception** | the **rate of release** `Imax/2(t75−t25)`: **r = 0.968 vs TI Imax**, against **r = 0.860** for the maximum concentration | Baek 1999 §5.2 [derived; the paper prints no correlation] |

**⇒ The mechanism is not partition and not reversible binding. Three candidates remain, and
all three now have direct support in this wave:**
1. **Background-odour masking** — the only mechanism that is compound-specific in a way that
   ignores hydrophobicity. Hong's own proposal; predicts the ester inversion (the one compound
   with no soybean counterpart); consistent with r = −0.05.
2. **Irreversible carbonyl chemistry** — Meynier's t-2-hexenal is retained **6.88×** by dairy
   protein against hexanal's **1.39×**, on a browning system; Anantharamkrishnan shows
   benzaldehyde forms a **+88 Da Schiff base at pH 7–8 and NOTHING at pH 3**, and citral
   **cross-links the protein** at 45 °C; Leksrisompong shows diacetyl binds caseinate at
   **pH 7 (P<0.05)** and **not at pH 5.5 (N.S.)**. **Four independent lines, one direction.**
3. **Delivery kinetics** — Baek shows perception tracks a **time derivative** of concentration,
   which no equilibrium partition model computes at all.

---

# (C) WHERE THIS WAVE EXTENDS, AND WHERE IT PARTIALLY CONTRADICTS, K2

| K2 claim | status after K4b |
|---|---|
| **§D.1 "No general matrix correction factor"** | ✅ **CONFIRMED on same-method data.** Hong's same-panel ratios still span 70× max/min with one inversion. |
| **§D.2(i) "No same-method matrix-vs-water pair exists anywhere"** | ⚠️ **NOW OBSOLETE — and this is a NEW-DATA correction, flagged loudly.** Two such sets now exist (Hong 2020, Leksrisompong 2010). K2 was true of its ten papers; it is no longer true of the corpus. |
| **§A.8's 1-σ bands of 27–41×** | ⚠️ **PARTIALLY CONTRADICTED by new measured data.** On a same-panel set the 1-σ band is **4.7×**, not 27–41×. **K2's bands were inflated by cross-study method noise.** The *tails* (70× max/min) are real; the *width of the bulk* was not. **This does not change K2's recommendation — a 4.7× band is still far too wide for a single factor, and the inversions are fatal to one regardless.** |
| **§D.2(ii) "Reclassify Brewer 1995 as `dose_added_pre_cook`"** | ✅ **STRENGTHENED.** Hong's properly-controlled hexanal ratio (132.5×) sits an order of magnitude below Brewer's 1 304×, exactly as predicted. |
| **§B.4 "Reversible binding is a single-digit-to-low-tens factor"** | ✅ **CONFIRMED twice more.** Meynier: **1.07–1.39×** at 34 g/L dairy protein. Leksrisompong: **0.79–1.69×** at 10–30 g/L caseinate. |
| **§B.5 "Headspace suppression ≠ perceived suppression"** | ✅ **CONFIRMED twice more, and hardened.** Leksrisompong: pH changes K in 4/12 comparisons and BET in **0/12**. Baek: a **simultaneous, within-subject** demonstration — 1.5× n.s. in concentration against 3.0× P<0.001 in perception. |
| **§D.4.5 "Stop deriving matrix odour activity from `f_free`"** | ✅ **CONFIRMED, and now measured rather than inferred.** Leksrisompong §6.2 quantifies the error at up to **30.6×** with **two of three signs wrong**. |
| **§D.3 "The α,β-unsaturation penalty of ≈2–3× is the only defensible parametric term"** | ✅ **STRENGTHENED to four independent lines** (K2's two matrices + Meynier's 6.88× vs 1.39× + Anantharamkrishnan's citral cross-linking). ⚠️ But Meynier's alkenal/alkanal ratio is **4.9×**, above K2's 2–3×; the term's magnitude is matrix-dependent. |
| **§B.3 "`method` must be a first-class field on every aldehyde binding record"** | ✅ **CONFIRMED and extended to a NEW field: `pH`.** Anantharamkrishnan shows carbonyl–lysine adduct formation is **abolished at pH 3**; Leksrisompong shows the same for diacetyl–caseinate at pH 5.5. **Propose adding `pH` alongside `method` on aldehyde binding records.** |
| **§B.6 "Saliva thiol binding is STRANDED"** | ⚠️ **NOW SCOPED.** Baek 1999 excludes **mucosal** binding for a neutral ester over a 15× dose range in four panellists. **⇒ Mucosal/salivary interaction is compound-class-specific, not a general multiplicative term.** A repo that applies one saliva factor to all odourants is refuted by the Baek/Starkenmann pair. |

---

# (D) THE PARTITION-CONSTANT PROBLEM — a new, quantified contributor to the hexanal over-prediction

## D.1 Meynier's own measured air/water K is 6.24× below the Henry's constants in its own Table I

| compound | Table I Henry ⇒ K_aw(25 °C) | Meynier measured ⇒ K_aw(25 °C) | ratio |
|---|---:|---:|---:|
| isoamyl acetate | 2.399e-2 | 3.549e-3 | 6.76× |
| amyl acetate | 1.586e-2 | 2.281e-3 | 6.95× |
| ethyl pentanoate | 1.520e-2 | 2.141e-3 | 7.10× |
| **hexanal** | **8.706e-3** | **1.911e-3** | **4.56×** |
| t-2-hexenal | 1.999e-3 | 3.214e-4 | 6.22× |
| | | **geometric mean** | **6.24× (4.56–7.10)** |

**Five compounds, two chemical classes, one direction, a 1.6× spread. That is a systematic
calibration offset, not chemistry.** Both papers calibrate a **gas-phase headspace injection
against liquid-phase standards** (Meynier: aroma in cyclohexane; Leksrisompong: methanol or
diethyl ether). Leksrisompong lands **6.1–17.3×** below its own cited comparisons for the same
reason. **Two labs, two decades, same direction, same method flaw.**

## D.2 The literature range on hexanal air/water K

| source | K_aw | provenance |
|---|---:|---|
| EPI-Suite Henry, via Meynier Table I | **8.7 × 10⁻³** @25 °C | estimator |
| **Meynier measured** | **1.9 × 10⁻³** @25 °C (2.5 × 10⁻³ @30 °C) | static HS, liquid-calibrated FID |
| Hall & Anderson 1983, via Meynier | **1.8 × 10⁻²** @30 °C | + ammonium sulphate (salting-out) |
| **total spread** | **9.5×** | — |

**⇒ If the repo computes matrix headspace ppb by multiplying a liquid concentration by an
estimator-class Henry constant (≈8.7 × 10⁻³), that alone sits 4.6× above the only direct
measurement in this wave — before any binding or matrix term is applied.**
⚠️ **This is NOT a licence to swap the constant.** The 6.24× is systematic across five
compounds and is at least as likely to be a Meynier method bias. **The finding is that the
uncertainty on air/water K is ±0.5 decades, and the repo should carry that band.**

## D.3 The full validated K(T) surface, reconstructed

`Meynier2002_extraction.md` §7.2 delivers a **5 compound × 4 media × 6 temperature (30–80 °C)
grid**, reconstructed as `K(T) = K(303.15)·exp(−ΔH_aff/R·(1/T − 1/303.15))` from Table III +
Table IV, and **validated against four of the paper's own printed prose checkpoints to within
2–7 %**. Three cells are excluded by the source's own ΔH footnotes (ethyl pentanoate water
>50 °C and skim milk >60 °C; t-2-hexenal skim milk >60 °C).

**Hexanal, the row the repo needs** (dimensionless air/matrix):

| medium | 30 °C | 40 °C | 50 °C | 60 °C | 70 °C | 80 °C |
|---|---:|---:|---:|---:|---:|---:|
| **water** | **2.50e-3** [M] | 4.17e-3 | 6.74e-3 | 1.06e-2 | 1.62e-2 | 2.42e-2 |
| **skim milk** | **1.80e-3** [M] | 2.84e-3 | 4.36e-3 | 6.51e-3 | 9.51e-3 | 1.36e-2 |
| **anhydrous milk fat** | **8.20e-5** [M] | 1.33e-4 | 2.08e-4 | 3.18e-4 | 4.73e-4 | 6.89e-4 |
| **full-fat cream (30 % fat)** | **1.60e-4** [M] | 2.55e-4 | 3.95e-4 | 5.96e-4 | 8.78e-4 | 1.26e-3 |
| skim-milk retention vs water | 28.0 % | 31.9 % | 35.4 % | 38.5 % | 41.3 % | **43.8 %** (paper prints 42 % ✅) |

---

# (E) THE BROWNING-KINETICS TRIO — one real food, two model systems, one contradiction

| | **Pereyra Gonzales 2010** | **Miao & Roos 2004** | **Lievonen & Roos 2002** |
|---|---|---|---|
| matrix | **real skim milk powder** (51.7 % lactose, 34.2 % protein) | lactose/trehalose 1:1 | maltodextrin M100 **or** PVP-40 |
| reactants | the food's own lactose + protein lysine | **D-xylose + L-lysine, 5 % of solids** | **D-xylose + L-lysine, 10 % of SORBED WATER** |
| response | **available lysine loss (a reactant)** | OD 280/420 nm (a product) | OD 280 nm (a product) |
| order | **pseudo first-order** | zero-order | pseudo-zero-order |
| T range | 30–60 °C | 40–90 °C | 10–110 °C |
| **data availability** | **ALL 18 k + 5 Ea TABULATED with 95 % CI** ✅ | **0 of 36 k tabulated; 4 Ea in Figure 5 only** | **0 of 42 k tabulated; NO Ea EXISTS** |
| **k recovered this wave** | **18 of 18** | **11 digitised ±10–30 %, 25 unreadable** | **17 digitised ±10–30 %, 25 unreadable** |
| R² of the rate fits | not stated per row | **0.9425–0.9990** | ⚠️ **0.803–0.995** |
| **does Ea depend on water?** | **NO** (121.1–135.4, P > 0.05, over a_w 0.43–0.98) | **YES** (280 nm: 155.8 → 117.8; 420 nm: 183.6 → 128.6, over 3.78 → 6.86 g/100 g) | **no Ea computed** |
| **is there an intermediate-a_w maximum?** | **NO in real milk powder; YES in the same lab's matched buffered model system** | rate rises monotonically with water | rate rises monotonically with a_w |
| **is `T − Tg` sufficient?** | break at **Tg + 7 and Tg + 13 °C** | rise at **Tg + 15–30 °C** | **NO — 6.0× spread at matched `T − Tg`**; browning observed at **Tg − 40 °C** |

## E.1 The Ea contradiction, and the proposed reconciliation

**Miao's Ea trend is with WATER CONTENT at conditions that are all rubbery. Pereyra Gonzales's
null is over a_w at conditions that are all well rubbery (Tg ≤ 9 °C against 37–60 °C) — and
the one sample of hers that IS near Tg (a_w 0.33, Tg 33 °C) shows the largest effect in her
whole paper: a 56× rate collapse at 37 °C that vanishes entirely at 50 °C.**

**⇒ Proposed corpus rule for the orchestrator: scope `Ea_depends_on_water` to
`T − Tg < ~20 °C` and treat it as absent well above Tg.** Both datasets are consistent with
that; neither paper states it. Two further reconcilers are also live and non-exclusive:
**(a)** the response variables differ (reactant loss vs product colour) and **must not be
pooled**; **(b)** the water axes differ (a_w vs water content) and are not interchangeable in
a milk powder, which reports only a_w.

## E.2 `T − Tg` is NOT a sufficient reduced variable — two independent papers, one direction

Lievonen & Roos state it in their abstract and this wave quantifies it: **at matched
`T − Tg` ≈ +20 °C, the six models span 6.0× in rate**, ordered by matrix identity first and
a_w second, with **the lowest-a_w model fastest**. Miao & Roos, same research group, different
matrix, reach the same conclusion with the same direction (**lower water is faster at matched
`T − Tg`**). And **browning proceeds 40 °C BELOW Tg** (measured) and **80 °C below** (cited).
**⇒ Any repo term that gates or collapses matrix state onto `T − Tg` alone is refuted by two
papers.**

## E.3 Two under-appreciated datasets

- **Lievonen Table 1 is a complete, tabulated pH-drift dataset**: pH falls **1.26 to 4.43 units**
  during browning in six unbuffered amorphous models, and **0.5 M citrate cuts the drift to
  0.13–0.18 units (a 20–30× reduction)**. The corpus has almost nothing on how far pH falls
  during browning. ⚠️ But the same buffer **moves the maltodextrin Tg by 32 °C (58 → 26)**, so
  the "control" is not a control.
- **Pereyra Gonzales's crystallisation timings** establish that lactose crystallisation is
  **not** confounding her rate constants (no crystallisation at a_w 0.33/50 °C in 200 h; onset
  at 40 h at 60 °C, after >50 % lysine loss).

---

# (F) NAMED LAUNDERING HAZARDS — wave K4b's contribution to the "342/200 list"

Ranked by how likely a future wave is to import them. Full detail in the per-paper dossiers.

| # | claim, as printed | reality | source |
|---:|---|---|---|
| **1** | **"A high partition coefficient value ideally correlates to a low threshold value ... and we observed these trends in our study"** | Observed for **one of three compounds**. r = **−0.91** (δ-decalactone) but **+0.60** (diacetyl) and **+0.43** (furaneol). **Diacetyl in oil has the highest K and the highest threshold simultaneously.** | Leksrisompong 2010 p. 364 |
| **2** | Hong's **Z-value column** (up to 24,942), presented as `Z = (BET2−BET1)/√(SD1²+SD2²)` | Every entry equals **BET2 − BET1** to within 1.4 %; the denominator is 1.00. It is a difference in µg/kg. **Diagnosed: Hong inherited a correct statistic from Leksrisompong 2010 and broke it by log-scaling only the denominator.** | Hong 2020 Table 2 |
| **3** | **"BET in the soybean-based model was 1.5 µg/kg (p < .05)"** — the one inverted compound | **No p-value is computed anywhere in the paper.** On the SD convention Hong says it used, the correct Z is **−0.79, not significant**; it reaches significance only on the SE convention Hong says it did NOT use. | Hong 2020 p. 4 |
| **4** | **Leksrisompong Table 5's `K_HS/matrix E⁵` header** | Must be read as **K = printed × 10⁻⁵**. Taken at face value it gives partition coefficients of order 10⁵ — **wrong by 10 orders of magnitude**. | Leksrisompong T5 |
| **5** | **"Our results agreed with previous published values: Salvador et al. (1994) also found diacetyl in sunflower oil to have higher K (6.3E-4) compared with water (4.5E-4)"** | Only the **ordering** agrees. The magnitudes differ by **13–17×**. | Leksrisompong p. 361 |
| **6** | **"The activation energies of different model systems varied from 117.8 to 183.6 kJ/mol"** | There are exactly **four** Ea values, not six — the 5.85 g/100 g model was **never Arrhenius-fitted**. **183.6 is never attributed** to a water content or wavelength, and **128.6 never appears in the running text at all.** Both are recoverable only from a 300 dpi render of Figure 5. | Miao & Roos p. 5255 |
| **7** | **"The rate (and, therefore, extent) of covalent bond formation increased with pH, temperature and water activity"** | **No rate was measured.** No rate constant, yield, calibration or internal standard exists anywhere in the paper; all conclusions read relative peak heights in MaxEnt-deconvoluted spectra. For **three of four compounds the 45 °C spectra are at noise level** because the products cross-linked out of the mass window. | Anantharamkrishnan 2020 |
| **8** | **E_a = 498.0 kJ mol⁻¹** for lysine loss at a_w 0.33, 40–50 °C | **No confidence interval, fitted over a 10 °C window, and unverifiable** — the 40 and 45 °C rate constants exist only in Figure 4. It is **3.7×** the largest ordinary value in the same paper. It describes a viscosity transition, not a chemical barrier. **Must never enter an Ea prior.** | Pereyra Gonzales p. 44 |
| **9** | Hong's **"± SD of individual threshold value"** with no stated scale | Cannot be linear (implied RSD ranges **0.0024 % to 120 %** across rows). Must be a **log SD**, base not stated. Ingesting as linear µg/kg would understate threshold uncertainty by **four orders of magnitude** on the soybean rows. | Hong 2020 Table 2 note |
| **10** | Baek's **"there was no binding of volatile to protein in the gel, nor to mucous membranes"** | True for **one neutral ester in gelatine at 100–1 500 mg/kg** — i.e. 10³–10⁶× above threshold, where any saturable site is saturated, in a collagen hydrolysate with essentially no hydrophobic pockets. **Not a general exclusion**; Meynier measures **6.88×** suppression of an alkenal by dairy protein on the same class of experiment. | Baek 1999 abstract |
| **11** | **"the gradient correlated better with the sensory values than TR Imax"** | The paper computes **no correlation at all** between any instrumental and any sensory variable. (Computed here: **r = 0.968 vs 0.860** — the claim is correct but was **unsupported as published**.) | Baek 1999 p. 159 |
| **12** | Meynier abstract: **"nearly 90 %"** t-2-hexenal retention; **"6 % for isoamyl acetate to 40 % for hexanal"** | Results say **"approximately 85 %"** and Table III gives **85.5 %**. And the 6 % / 40 % range mixes **30 °C and ~80 °C** values as if comparable. | Meynier abstract |
| **13** | Meynier **Table IV** — 20 enthalpies presented uniformly | **Three are fitted over restricted ranges** per the table's own 6 pt superscript footnotes (ethyl pentanoate water **30–50 °C**, ethyl pentanoate skim milk and t-2-hexenal skim milk **30–60 °C**). Extrapolating them to 80 °C is wrong by up to **1.7×**. | Meynier T4 footnotes |
| **14** | Lievonen abstract: **"The NEB rates were higher in PVP than in MD models"** at matched conditions | The **MD models start 0.94–3.28 pH units MORE ALKALINE** (the paper's own Table 1), and the buffered control arm **moves the MD Tg by 32 °C**. The confound runs in the direction that would **enlarge** the true effect, but it is uncontrolled. | Lievonen abstract vs T1/T2 |
| **15** | **Leksrisompong's 0/10/20 % fat series** presented as a fat effect | Aqueous protein **falls** as fat rises (constant 1 % total emulsifier) and **droplet size rises 2.7×** from 10 % to 20 %. Three variables move together; the authors concede the effect "is unknown". | Leksrisompong pp. 352, 358 |
| **16** | **Meynier / Leksrisompong absolute partition coefficients** | Both calibrate a **gas-phase injection against liquid-phase standards**. Meynier is **6.24×** below its own Table I Henry values; Leksrisompong is **6.1–17.3×** below its own cited comparisons. **Ratios are sound; absolutes carry a factor-of-10 downward bias risk.** | §D.1 |
| **17** | Anantharamkrishnan **ref. 2**: *"Gremli, H. A.; Dimitrova, N.; ... Interaction of Flavor Compounds with Soy Protein. Chem. Res. Toxicol. 2010, 51, 1050−1059"* | The correct citation is **Gremli, H. A. (1974), *J. Amer. Oil Chem. Soc.* 51, 95A–97A**, single-authored. Wrong journal, wrong year, wrong pages, fabricated co-authors. ⚠️ Meynier 2002 cites the same paper **correctly**. | Anantharamkrishnan ref. 2 |
| **18** | **Miao Table 4's 5.85 g/100 g row, T_ref = 50 °C** | Violates the paper's own stated rule (*"the experimental temperature closest to the predicted Tg"*) — Tg is 34.8 °C, so T_ref should be **40**. Unexplained. | Miao T4 |
| **19** | Hong's Table 1 dosing ranges for **butyric acid/water (81×)**, **2,3-DMP/soybean (81×)**, **DMDS/soybean (72.9×)** | Inconsistent with the stated **7-step** series (should be 729×). Two are 5-step; the DMDS lower bound (400) is the only value in the table that is not an exact member of its own geometric series — **360** would fit. | Hong T1 |
| **20** | **"Guadagni, Buttery & Turnbaugh (1972) ... *Journal of Agricultural and Food Chemistry*, 23, 1435–1444"** | Correct venue is ***J. Sci. Food Agric.*** 23, 1435–1444. **Second corpus paper to mis-cite the primary source of four of K2 §A.8's six aqueous anchors.** | Hong 2020 ref. list |

**Plus a first-class calibration fact, not a hazard:** two independent papers in this wave
establish that **cross-study aqueous orthonasal thresholds are reproducible to no better than
~3 orders of magnitude.** Hong quotes two water values for **butyric acid** differing by
**5,540×** and two for **3-methylbutanal** differing by **800×**, and calls the spread "within
the range of variability"; Leksrisompong quotes a published **furaneol** range of
**1 to 1,700 ppb (1,700×)**. **This is the quantitative justification for
`k2_matrix_and_thresholds.md` §D.2(i) and it is why the two same-panel sets in §A matter.**

---

# (G) CONSOLIDATED NEW-PARAMETER TABLE — wave level

Per-paper tables with full conditions and anchors are in the eight dossiers. This is the
wave-level roll-up of what is newly available.

**Class: M = measured and printed · F = fitted and printed by the authors · D = digitised by
this wave (tolerance stated) · Z = derived by this wave, never printed · C = cited by the
source · Q = qualitative, no number in the source.**

| # | parameter | value | units | conditions | class | source |
|---:|---|---:|---|---|:--:|---|
| **THRESHOLDS — same-panel paired sets (44 new values)** | | | | | | |
| 1–20 | **BET, 10 compounds × 2 matrices** | 0.5 to 25,139.9 (full table §A.1) | µg/kg | 3-AFC ASTM E679-91, n = 20, ~22 °C, water vs soybean paste, matrix-matched blanks | M | Hong T2 |
| 21–44 | **BET, 3 compounds × 8 matrices** | 5.6 to 1,550, each with a linear SE (23.6–45.2 % RSE) | ppb | 3-AFC ASTM E679-91, n = 25, **22 °C**, duplicate on different days | M | Leksrisompong T3 |
| 45–54 | **matrix/water shift, soybean paste** | **0.155× to 2,035×**; geomean **277×** (9 elevated) | × | same panel, `cross_study_cross_method: FALSE` | Z | Hong |
| 55 | **log₁₀ SD of the soybean shift** | **0.668** (9 compounds) / **1.206** (all 10) | decades ⇒ **4.7×** / **16.1×** band | " | Z | Hong |
| 56 | **corr(log P, log matrix shift)** | **−0.052**; slope **−0.070 dec/log P** | r, n = 10 | " | Z | Hong |
| 57–77 | **matrix/water shift, 7 matrices × 3 compounds** | **0.54× to 23.5×**, **3 of 21 below 1** | × | same panel | Z | Leksrisompong |
| 78 | **per-compound panel dispersion** | **0.4–0.9, mean 0.65** ⚠️ **base assumed log₁₀, NOT STATED** | decades | n = 20 individual BETs | M (scale Z) | Hong T2 |
| **PARTITION — 44 new coefficients + a validated K(T) surface** | | | | | | |
| 79–98 | **K_air/medium at 30 °C**, 5 compounds × 4 media, each ± SD | **1.0e-5 to 4.5e-3** (full table `Meynier2002` §5) | dimensionless | static HS-GC-FID, 3 g in 22.4 mL, 30 min (60 min skim milk) | M | Meynier T3 |
| 99–118 | **ΔH_affinity**, 5 compounds × 4 media, each ± SD | **32.0 to 56.5** ⚠️ **3 cells restricted to 30–50/30–60 °C** | kJ·mol⁻¹ | van 't Hoff slope, 30–80 °C | F | Meynier T4 |
| 119–238 | **K(T) grid** | 5 × 4 × 6 = **120 values**, validated against 4 printed checkpoints to 2–7 % | dimensionless | 30/40/50/60/70/80 °C | **Z** | `Meynier2002` §7.2 |
| 239–262 | **K_HS/matrix**, 3 compounds × 8 matrices, each ± SD | **1.0e-7 to 4.8e-5** (⚠️ printed as `× E⁵`) | dimensionless | static HS-GC-MS SIM, 10 g in 40 mL, **40 °C** | M | Leksrisompong T5 |
| 263 | **estimator-vs-measurement offset on K_aw** | **6.24× (range 4.56–7.10)** | × over-statement | 5 compounds, 25 °C | Z | `Meynier2002` §7.3 |
| 264 | **absolute-scale offset vs cited K values** | **6.1–17.3× low** | × | vs Salvador 1994, Guyot 1996 | Z | `leksrisompong2010` §4.1 |
| 265 | **hexanal K_aw literature spread** | **1.9e-3 to 1.8e-2 = 9.5×** | dimensionless | 25–30 °C | Z + C | §D.2 |
| 266–271 | **skim-milk retention of hexanal** | **28.0 / 31.9 / 35.4 / 38.5 / 41.3 / 43.8** | % vs water | 30→80 °C; 80 °C validated against the printed 42 % | Z | §D.3 |
| 272–277 | **skim-milk retention of t-2-hexenal** | **85.5 → 83.7 — FLAT** | % vs water | 30→80 °C ⚠️ 70–80 °C outside ΔH validity; measured on a **browning** system at up to 5.8 g/L | Z | `Meynier2002` §7.2 |
| 278 | **water/skim-milk headspace suppression** | hexanal **1.39×**; **t-2-hexenal 6.88×**; esters **1.07–1.33×** | × | 30 °C, ~34 g/L dairy protein | Z | `Meynier2002` §5.1 |
| 279 | **max caseinate effect on headspace, 10→30 g/L** | **0.787× (diacetyl, pH 7)** to **1.694× (furaneol, pH 7 — ENHANCEMENT)**; **all N.S. at pH 5.5** | × | 40 °C | Z from M | `leksrisompong2010` §7 |
| 280 | **partition-model error on the oil/water threshold shift** | **30.6× / 7.2× / 2.4×**; over 21 cells **0.14× to 30.6×** | × | paired, same samples | **Z** | `leksrisompong2010` §6.2 |
| 281 | **corr(log K, log BET) over 8 matched matrices** | **+0.600 / −0.907 / +0.427** | r, n = 8 | diacetyl / δ-decalactone / furaneol | **Z** | `leksrisompong2010` §6.1 |
| 282 | Gremli 1974 soy retention ⇒ per-gram | hexanal **1.33e-2**, t-2-hexenal **4.90e-2** | L/g | 50 g/L soy protein, vacuum extraction | Z from C | `Meynier2002` §9 |
| 283 | emulsion droplet size | d[4,3] **0.59 µm** @10 % fat, **1.58 µm** @20 %; d[0.9] **1.10 / 1.64 µm** | µm | 80 MPa, pH-independent | M | Leksrisompong |
| **BROWNING KINETICS — 18 tabulated k, 5 + 4 Ea, 28 digitised k** | | | | | | |
| 284–301 | **k, available lysine loss, real skim milk powder** | **1.98e-5 to 2.74e-2**, each with a 95 % CI | **h⁻¹** | 6 a_w × 3 T (37/50/60 °C), pseudo-first-order, unbuffered | **M** | Pereyra Gonzales T1 |
| 302–306 | **E_a, lysine loss** | **121.1 ± 2.7 / 129.2 ± 3.7 / 131.4 ± 3.1 / 135.4 ± 2.9 / 128.5 ± 4.1** — **all five reproduce from the paper's own k to within 2.5 %** | kJ mol⁻¹ | a_w 0.43 / 0.52 / 0.69 / 0.85 / 0.98 | **F + Z** | Pereyra Gonzales T1 |
| 307 | **E_a is independent of a_w over 0.43–0.98** | **P > 0.05** | — | 37–60 °C, real food | M | Pereyra Gonzales |
| 308–310 | **E_a, a_w 0.33, three zones** | **65.2 / ⚠️498.0 / 79.9** (the last reproduces to 0.1 %) | kJ mol⁻¹ | 30–40 / 40–50 / 50–60 °C; breaks at **Tg + 7** and **Tg + 13 °C** | F | Pereyra Gonzales p. 44 |
| 311 | **rate collapse at a_w 0.33, 37 °C** | **56.3×** vs a_w 0.43; **no effect at 50 °C (0.99×)** | × | 4 °C above Tg vs 17 °C above | Z | Pereyra Gonzales |
| 312 | **dilution effect a_w 0.69 → 0.98** | **3.61× (60 °C), 3.75× (50 °C)** | × | temperature-independent | Z | Pereyra Gonzales |
| 313–316 | **E_a, NEB in lactose/trehalose model** | **117.8 ± 6.7 (R² 0.9896) · 155.8 ± 12.1 (0.9939) · 128.6 ± 11.3 (0.9983) · 183.6 ± 4.9 (0.9959)** — **two recovered from a 300 dpi render, absent or unattributed in the text** | kJ mol⁻¹ | 280 nm 6.86 g / 280 nm 3.78 g / **420 nm 6.86 g** / **420 nm 3.78 g** per 100 g solids | **F** | Miao Fig. 5 |
| 317 | **dEa/d(water content)** | **−12.3 (280 nm), −17.9 (420 nm)** | kJ mol⁻¹ per g water/100 g solids | 3.78 → 6.86 g/100 g | Z | Miao |
| 318 | **Ea(420 nm) − Ea(280 nm)** | **+27.8 (3.78 g), +10.8 (6.86 g)** | kJ mol⁻¹ | brown pigment vs early yellow | Z | Miao |
| 319 | **E_a at 5.85 g/100 g** | **DOES NOT EXIST IN THE SOURCE** | — | — | — | Miao §5.1 |
| 320–330 | **k, NEB, lactose/trehalose** | **0.58 to 24.0** at 70/80/90 °C | **OD units·g solid⁻¹·min⁻¹** | 3 water contents × 2 wavelengths | **D ±10–30 %** | Miao Fig. 4A |
| 331 | **25 further Miao k values** | **UNREADABLE** — on the zero line of a linear axis | — | 40/50/60 °C | — | Miao §6.1 |
| 332–348 | **k, NEB, MD and PVP models** | **2.5 to 88.5** | **OD units·h⁻¹** | 6 models × their own 7-temperature windows, 280 nm only | **D ±10–30 %** | Lievonen Fig. 5A |
| 349 | **25 further Lievonen k values** | **UNREADABLE** | — | — | — | Lievonen §4.1 |
| 350 | **activation energy in Lievonen & Roos 2002** | **NONE EXISTS** — Figure 4 is qualitative | — | — | — | Lievonen §7 |
| 351 | **spread of k at matched `T − Tg` ≈ +20 °C** | **6.0×** across six models | × | the quantification of "`T − Tg` is not sufficient" | **Z** | Lievonen §5.4 |
| 352 | **lowest `T − Tg` at which browning was measured** | **−40** (measured, PVP 0.2); **−80** (cited, Schebor 1999) | °C | — | M / C | Lievonen §5.2 |
| 353 | **position of the Arrhenius break relative to Tg** | **+7, +13** (real milk powder) · **+15 to +30** (lactose/trehalose) · *"slightly below, exactly at, or well above"* (MD/PVP) · **+10 to +40** (literature) | °C | **four datasets, four answers** | F / M / C | §E |
| 354–373 | **pH before and after browning, 10 amorphous systems** | before **6.03–9.42**, after **4.55–7.20**; **drift −1.26 to −4.43**, cut to **−0.13 to −0.18** by 0.5 M citrate | pH units | 20 % w/w solution, 70 °C, 5–95 h | **M** | Lievonen T1 |
| 374 | **citrate effect on maltodextrin Tg** | **−32 °C** (58 → 26) at 0.5 M, at constant water content | °C | ⚠️ the pH "control" changes the matrix state | Z from M | Lievonen T2 |
| 375–380 | **Tg, MD and PVP models** | MD **≈75 / 58 ± 3 / ≈48**; PVP **≈72 / 51 ± 2 / ≈40** | °C | a_w 0.23 / 0.33 / 0.44; two measured, four **recovered by cross-referencing Figs. 5A/5B ±3 °C** | M + **Z** | Lievonen T2, §4.2 |
| 381–395 | **Tg, DSC onset, lactose/trehalose systems** | **0.17 ± 0.36 to 59.78 ± 0.25** | °C | 5 RVP × 3 systems, triplicate | M | Miao T2 |
| 396 | **Gordon–Taylor k** | **7.3 ± 0.7** (model), **7.7 ± 0.8** (lactose) | — | with Tg(water) = −135 °C | F | Miao |
| 397–409 | **BET and GAB isotherm constants** | BET m_m **6.89**, GAB m_m **3.95** (+ 11 shape constants) | g H₂O/100 g solids | 25 °C, **a_w 0.114–0.441 only, n = 4** | F | Miao T1 |
| 410–415 | **water content vs a_w, MD and PVP** | MD **6.3 / 8.2 / 9.7**; PVP **8.2 / 12.8 / 17.3** (± 0.01–0.1) | g H₂O/100 g dry solids | a_w 0.23 / 0.33 / 0.44 at 24 °C | M | Lievonen |
| 416–421 | **WLF C1, C2 (both bases), 280 nm** | **8.1/87.3, 9.7/149.3, 16.6/229.1** and **5.3/49.6, 4.3/56.6, 8.0/95.4** | — | 3 water contents; **C1 doubles, C2 nearly triples** ⇒ not transferable | F | Miao T4 |
| 422–427 | **crystallisation onset** | lactose/trehalose **300 / 56 / 32 h** vs pure lactose **11 / 4 / 3 h** | h | RVP 54.5 / 65.6 / 76.1 %, 25 °C | M | Miao |
| 428–430 | **crystallisation in milk powder** | **none in 200 h** (a_w 0.33, 50 °C); **40 h** (a_w 0.33, 60 °C, after >50 % lysine loss); **23 h** (a_w 0.43, 37 °C) | h | ⇒ crystallisation does not confound the rate constants | M | Pereyra Gonzales |
| **PERCEPTION AND MECHANISM** | | | | | | |
| 431 | **in-nose Imax vs perceived intensity, same mouthfuls** | concentration falls **1.49–1.52×, NOT significant**; perception falls **3.00–3.28×, P < 0.001** | × | 2 → 8 % w/w gelatine gels, 300 mg/kg FFA, 11 panellists, simultaneous TR + TI | **M + Z** | Baek §5 |
| 432 | **rate of release falls** | **2.99×, P < 0.01** | × | `Gradient = Imax/2(t75−t25)` | M + Z | Baek §5.2 |
| 433–436 | **correlations the paper never computed** | gradient vs TI Imax **r = 0.968**; Imax vs TI Imax **r = 0.860**; gradient vs sensory score **0.886**; Imax vs sensory score **0.851** | r, n = 5 | digitised Fig. 7 ±3 pp | **Z** | Baek §5.2 |
| 437 | **binding exclusion, scoped** | equal headspace plateau at **20, 50 and 80 g/L gelatine** ⚠️ tolerance not stated | — | furfuryl acetate at **300 mg/kg**, 37 °C stirred water | M | Baek §4.1 |
| 438 | **mucosal-binding exclusion** | in-nose Imax **LINEAR** in gel concentration over **100–1 500 mg/kg**, 3 gels × 3 replicates, replicated in 4 panellists | — | ⚠️ no slope or intercept printed | M | Baek §4.2 |
| 439 | **in-vitro dissolution time** | **10 min** (2 % gel) vs **40 min** (8 % gel) | min | 6 g cube, 100 mL stirred water, 37 °C | M | Baek §4.1 |
| 440–453 | **covalent adduct mass shifts on β-lactoglobulin** | **+88, +99, +152, +198, +268, +297, +324, +402, +554, +46, +78, +47, +32 increments, −18** | Da, **± 1 Da** | 1 wt % BLG variant A (18 362 Da), **12 ppth flavour (saturating)**, UPLC-ESI-qTOF | M | Anantharamkrishnan §3 |
| 454 | **carbonyl–lysine adduct formation at pH 3** | **NONE** (benzaldehyde, citral, DMTS); **AITC forms adducts more slowly** | directional gate | 24 h, ~20 °C | **Q** | Anantharamkrishnan §4.1 |
| 455 | **adduct formation vs pH and temperature** | **increases pH 3 → 7 → 8**; **increases 4 → 20 °C** ⚠️ **the authors' own claim stops at 20 °C** | directional gate | ⚠️ 45 °C spectra are at noise level for 3 of 4 compounds | **Q** | Anantharamkrishnan §4.2, §5 |
| 456 | **water activity gates Schiff-base formation** | benzaldehyde and DMDS: **nothing from a_w 0.11 to 0.64, a trace at 0.75**; AITC: monotone rise; citral: **no effect in 48 h** | directional gate | freeze-dried BLG, 6-week equilibration | **Q** | Anantharamkrishnan §4.3 |
| 457 | **α,β-unsaturated aldehyde cross-links protein** | citral, at 45 °C — monomer peak vanishes into noise | mechanism | pH 7, 12 ppth, 24 h | **Q** | Anantharamkrishnan §4.2 |
| 458 | **rate constants, yields, half-lives or equilibrium constants for covalent adduct formation** | **NONE EXIST IN ANANTHARAMKRISHNAN 2020** | — | — | — | §1 of that dossier |
| **CALIBRATION FACTS** | | | | | | |
| 459 | **cross-study aqueous threshold reproducibility, butyric acid** | **5,540×** between two cited water values | × | — | Z | Hong §5 |
| 460 | cross-study aqueous reproducibility, 3-methylbutanal | **800×** | × | — | Z | Hong §5 |
| 461 | published furaneol water-threshold range | **1 to 1,700 ppb = 1,700×** | × | — | C | Leksrisompong §8 |

---

# (H) PROPOSED FIT / HOLD-OUT ROLE PER DATASET — **DRAFT FOR ORCHESTRATOR**

> ⚠️ **ALL EIGHT PAPERS ARE NEW SOURCES.** Under the repo's
> `docs/reference/FIT_HOLDOUT_DECLARATION.md` regime, **a declaration amendment must be
> approved before any wave fits any row below.** This wave did not edit the declaration, did
> not touch `src/`, `tests/`, `results/` or `data/benchmarks/`, and did not commit. **No
> section below is authorisation.**

## H.1 Recommended FIT-ELIGIBLE (subject to amendment) — in priority order

| # | dataset | rows | why it is the strongest candidate |
|---:|---|---:|---|
| **1** | **Pereyra Gonzales 2010, Table 1 — 15 rate constants at a_w 0.43–0.98** | 15 | **The best-quality dataset in the wave.** Tabulated with 95 % CIs, in a **real food matrix**, on a **reactant-loss** response variable, and **verified twice internally** (the unit resolves by an independent arithmetic check; all five Ea reproduce from the k values to within 2.5 %). If the orchestrator amends the declaration for one thing this wave, this is it. |
| **2** | **Meynier 2002, Tables III + IV — matrix/water RATIOS and 17 of 20 ΔH_aff** | 37 | The 6.24× absolute-scale offset (§D.1) **cancels exactly in a within-study ratio**, leaving the ratios sound. **Exclude the 3 ΔH cells the source's own footnotes restrict.** |
| **3** | **Miao & Roos 2004, Figure 5 — 4 activation energies** | 4 | Printed with ± and R²; two are confirmed by the running text. ⚠️ Ingest with `Ea_biased_low_by_uncorrected_come_up_time`. |
| **4** | **Leksrisompong 2010, Table 5 — K_water/K_matrix ratios** | 21 | Same ruling as #2: ratios survive the 6–17× absolute offset. |
| **5** | **Lievonen & Roos 2002, Table 1 — 20 pH values before/after browning** | 20 | **The most under-appreciated dataset in the wave.** The corpus has almost nothing on how far pH falls during browning; this gives −1.26 to −4.43 units in six unbuffered amorphous systems, plus a 20–30× buffered control. |
| **6** | **Miao Tables 1–3 + Lievonen §3.3 — sorption isotherms, Tg and Gordon–Taylor** | 40 | Clean DSC and gravimetric data with SDs. Scope the isotherms to a_w 0.114–0.441. |
| **7** | Pereyra Gonzales, 5 Ea | 5 | ⚠️ **Redundant** with #1 — they are derived from the same k values. **Fit the k, use the Ea as a check; do not fit both.** |

## H.2 Recommended HOLD-OUT — the acceptance and falsification set

| # | dataset | rows | proposed use |
|---:|---|---:|---|
| **1** | **Hong 2020, all 20 BETs and the 10 paired ratios** | 30 | **The wave's flagship hold-out.** Suggested acceptance criterion for any matrix-correction work: **reproduce ≥7/10 ratios within 5× (the measured 1-σ band) AND reproduce the sign on all 10, including the ester inversion.** Nothing in the repo today can do the last part. |
| **2** | **Leksrisompong 2010, all 24 BETs** | 24 | The second same-panel matrix threshold set, and the only one with an **oil** arm. |
| **3** | **Leksrisompong §6 — the 3 correlations and the 21-cell error table** | 24 | **The wave's primary falsification test for `f_free` / partition→perception.** Guard: a model that cannot produce a **positive** log K / log BET correlation for a hydrophilic carbonyl (**r = +0.60, diacetyl**) is refuted. |
| **4** | **Hong §6.2 — corr(log P, matrix shift) = −0.05** | 1 | Guard: **any shipped matrix term that is a monotone function of log P is refuted.** |
| **5** | **Baek 1999 §5.1 — 1.5× n.s. concentration vs 3.0× P<0.001 perception, simultaneous** | 2 | **The cleanest falsification of concentration→intensity mapping in the corpus**, because TR and TI were recorded on the same mouthfuls of the same panellists. |
| **6** | **Lievonen §5.4 — 6.0× spread at matched `T − Tg`**, and **§5.2 — browning at `T − Tg` = −40 °C** | 2 | Guards against any repo term that reduces matrix state to `T − Tg` or gates Maillard on `T > Tg`. **Independently confirmed by Miao §6.3.** |
| **7** | **Pereyra Gonzales §4.1 — the flat-then-fall a_w profile** | 1 | Guard: a repo `f(a_w)` term must be **flat within error over a_w 0.31–0.71 and 3.6–3.8× lower at a_w 0.98** in a milk-powder matrix, and **must NOT produce a hump** — the hump exists only in the buffered model system. |
| **8** | **The Ea-vs-water contradiction (Miao yes / Pereyra Gonzales no)** | 3 | **Neither should be fit until the orchestrator rules on scope.** §E.1 proposes scoping `Ea_depends_on_water` to `T − Tg < ~20 °C`. |
| **9** | **Meynier §8.1 — two-phase emulsion model fails in both directions** | 5 | Guard on any composition-only partition model. |
| **10** | **Anantharamkrishnan C1 — no carbonyl–lysine adduct at pH 3** | 1 | Guard: **any repo aldehyde-loss-to-protein term must go to zero at acid pH.** Corroborated in-wave by Leksrisompong's diacetyl/caseinate pH result. |
| **11** | **Meynier hexanal water column (2.5e-3 @30 °C → 2.42e-2 @80 °C)** | 6 | Cheap decisive guard: if the repo's hexanal partition term at 30 °C in water falls outside **1.9e-3 to 1.8e-2**, it is outside the entire measured literature. |

## H.3 Recommended QUARANTINE, PRIOR-ONLY, or REJECT

| dataset | rows | ruling |
|---|---:|---|
| **Hong's printed Z column** | 10 | **REJECT** — arithmetically not the stated statistic (hazard 2). |
| **Hong's ± column** | 20 | **BLOCKED pending clarification of the log base**; if confirmed log₁₀, propose role = PRIOR on sensory-threshold uncertainty (mean 0.65 decades). |
| **Pereyra Gonzales E_a = 498.0 kJ/mol** | 1 | **REJECT** — no CI, unverifiable, 3.7× the largest ordinary value, describes a viscosity transition (hazard 8). |
| **All 28 digitised browning rate constants (Miao Fig. 4A, Lievonen Fig. 5A)** | 28 | **PRIOR / SANITY-CHECK ONLY.** ±10–30 % from linear-axis figures, an R² floor of 0.803 in one paper, and a 1.4–1.7× disagreement between Miao's own Figure 4A and Figure 5. **Fitting these is fitting noise.** |
| **The 50 unreadable rate constants** | 50 | **Recorded as `unreadable`** per the wave brief. ~32 are irrecoverable at any resolution (they lie on the zero line of a linear axis). |
| **All absolute partition coefficients (Meynier T3, Leksrisompong T5)** | 44 | **INGEST only with `absolute_scale_suspect: true, offset_vs_literature: 6-17x_low`** (hazard 16). |
| **Meynier's t-2-hexenal / skim-milk rows** | 8 | **QUARANTINE** — measured on a visibly browning system at up to 5 800 mg/L. |
| **Leksrisompong's 12 emulsion rows** | 12 | **QUARANTINE from any pure-fat-dependence fit** — fat, aqueous protein and droplet size are confounded by design and the authors concede it. |
| **Lievonen's PVP/MD rate ratio (4.0–9.8×)** | 5 | **QUARANTINE** — pH-confounded by 0.94–3.28 units; the buffered "control" moves Tg by 32 °C. |
| **Miao's WLF coefficients** | 18 | **REJECT as transferable parameters; RETAIN as evidence against universal WLF.** C1 doubles and C2 nearly triples over 3.1 g/100 g of water. |
| **Meynier Table I (EPI-Suite/UNIFAC)** | 35 | **REJECT as measurements; ingest only as `cited_from_EPI_Suite`.** |
| **Anantharamkrishnan 2020, entire paper** | — | **NO DECLARATION AMENDMENT NEEDED.** There is nothing fittable — no rate, yield, calibration or internal standard exists. **Record as a `mechanism_reference` in the citation graph, not a data source.** The 14 mass shifts are structural metadata; C1–C9 are directional gates. |
| **Baek's 0.6 min TI:TR crossover** | 1 | **REJECT** — R = 0.421, i.e. R² = 0.18; the authors themselves counsel caution. |
| **Baek's binding exclusion** | 1 | **Ingest ONLY as** `binding_excluded: {compound_class: neutral_ester, protein: gelatine, concentration: 100-1500 mg/kg}` — three reasons it does not generalise. |

## H.4 New provenance fields this wave requires

Extending `k2_matrix_and_thresholds.md` §D.4.2:
- **`cross_study_cross_method: false`** — the first records in the corpus that can carry it
  (Hong's 20 BETs, Leksrisompong's 24 BETs).
- **`dispersion_scale`** — `log10 | linear | not_stated`. Hong and Leksrisompong report
  dispersion on **different scales** for the same kind of quantity, and Hong does not say which.
- **`absolute_scale_suspect` / `offset_vs_literature`** on any static-headspace partition
  coefficient (§D.1).
- **`pH` as a first-class field on aldehyde binding records**, alongside the `method` field K2
  §B.3 already required. Two papers in this wave show carbonyl–protein binding is abolished at
  acid pH.
- **`response_variable`** on every browning rate constant — `reactant_loss` vs `product_colour`.
  §E.1 argues these must not be pooled.
- **`equilibration_asymmetric`** on paired matrix/water thresholds (Hong doses water at 0 h and
  soybean at 24 h/4 °C).
- **`dose_saturating`** on any covalent-adduct record (Anantharamkrishnan's 12 ppth is 10⁶–10⁷×
  a food-relevant loading).

---

# (I) THE SIX RETRIEVALS WORTH REQUESTING NEXT — wave-level, ranked

1. **Anantharamkrishnan, V., Hoye, T. & Reineccius, G. (2020)**, *J. Agric. Food Chem.* **68**,
   6395–6402 — *"Covalent Adduct Formation Between Flavor Compounds of Various Functional Group
   Classes and the Model Protein β-Lactoglobulin."* The parent survey covering **"various
   functional group classes"** — **the only plausible place an n-alkanal, and possibly hexanal
   itself, appears.** Directly targets the repo's hexanal-over-prediction question, which the
   2020 paper in hand cannot answer because it studies benzaldehyde, citral, AITC and DMTS/DMDS
   and **no straight-chain aldehyde**.
2. **Kim, M. K., Drake, S. L. & Drake, M. A. (2011)**, *J. Sensory Studies* **26**(4), 278–290 —
   orthonasal thresholds in **three matrices (water, pH 5.5 buffer, oil)** by the same protocol
   lineage as Hong and Leksrisompong. **A third same-panel multi-matrix paired set**, with an
   oil arm directly comparable to Leksrisompong's and to the Guadagni paraffin-oil values held
   second-hand in `k2_matrix_and_thresholds.md` §A.3.
3. **Malec, L. S., Pereyra Gonzales, A. S., Naranjo, G. B. & Vigo, M. S. (2002)**,
   *Food Res. Int.* **35**, 849–853 — the **phosphate-buffered lactose–casein model system** run
   by the same laboratory at the same temperatures and water activities as the real milk powder.
   **A near-perfect matched model-vs-real-food pair on Maillard kinetics**, and the source of
   §E's intermediate-a_w-maximum contradiction. Nothing else in the corpus provides that control.
4. **Hall, G. & Anderson, J. (1983)**, *Lebensm. Wiss. Technol.* **16**, 362–366 — *"Volatile fat
   oxidation products. II. Influence of temperature on volatility of saturated, mono- and
   di-unsaturated aldehydes in liquid media."* The **only other direct measurement of hexanal
   and t-2-hexenal air/water K vs temperature**, disagreeing with Meynier by 7.2–8.6× (§D.2).
   Also likely carries a **t,t-2,4-decadienal** partition value — the compound that is the
   extreme outlier in all three media of K2 §A.8 and for which the corpus has **no** partition
   constant at all.
5. **Lievonen, S. M., Laaksonen, T. J. & Roos, Y. H. (2002)**, *J. Agric. Food Chem.* **50**,
   7034–7041 — *"Nonenzymatic browning in food models in the vicinity of glass transition:
   effects of fructose, glucose, and xylose as a reducing sugar."* The third paper of the
   Roos-group trilogy and the only one that **varies the reducing sugar**. Would supply a
   **pentose-vs-hexose comparison in an amorphous matrix** — directly relevant to the pentose ≫
   hexose ordering that `k2_matrix_and_thresholds.md` §C.2 found **inverted** in a real extruder.
6. **Guyot, C. et al. (1996)** — the source of both the log P values and a δ-decalactone water
   partition coefficient 6.1× above Leksrisompong's, and (per Leksrisompong p. 349) itself a
   **paired partition/sensory-intensity study**. Would adjudicate the §D.1 absolute-scale
   question and add a third paired instrumental/sensory dataset.

*Runners-up:* **Fares et al. (1998)** — quantitative **covalent** diacetyl–caseinate binding
with two site classes and a glycosylated-caseinate control, i.e. the quantitative covalent
flavour–protein paper Anantharamkrishnan is not; **Gremli (1974)**, *JAOCS* **51**, 95A–97A —
soy-protein hexanal (37–44 %) and t-2-hexenal (68–75 %) retention at 50 g/L, on the repo's
plant-protein path; **Karmas, Buera & Karel (1992)** — cited by four of this wave's eight
papers and the source of the "10–40 °C above Tg" rule; and **Linforth, Baek & Taylor (1998/99)**,
*Food Chem.* — Baek's in-press companion, the only route to **absolute** release-gradient values.
