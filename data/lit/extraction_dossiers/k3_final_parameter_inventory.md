# K3 — FINAL PARAMETER INVENTORY
### The single consolidated document for the build waves. 2026-08-28.

**Scope.** This file consolidates every extraction dossier produced by waves **Z0, Z1, Z2, Z3,
K1, K2 and K3** into four things a build wave can act on:

- **(A) the MASTER PARAMETER TABLE** — every usable rate constant, activation energy,
  stoichiometry, binding constant, threshold and threshold ratio, with value, units,
  conditions, source anchor and a usability verdict;
- **(B) the DIRECTIONAL / ORDINAL CONSTRAINT LIST** — the structural facts the model must
  satisfy whether or not any parameter is fitted;
- **(C) the DECLARED GAPS** — carried **verbatim** from the analysts who declared them;
- **(D) a proposed FIT vs HOLD-OUT SPLIT** — **DRAFT FOR ORCHESTRATOR**, not applied.

**Nothing in the repository was written, staged or modified by wave K3.**

## HOW TO READ THIS FILE

**Provenance codes**, used on every numeric row:

| code | meaning |
|---|---|
| **[M]** | measured / reported directly by the paper's authors |
| **[C]** | cited by the paper from another source — the paper did not measure it |
| **[F]** | fitted by the paper's authors (a regression parameter) |
| **[D]** | derived by an extraction analyst from the paper's own data; **never printed by the authors** |

**Usability verdicts:** `USE` (ingest as printed) · `USE-Q` (ingest with the stated
qualification) · `RATIO-ONLY` (the absolute value is not commensurable; the ratio is)
· `STRUCTURAL` (not a number to fit, a shape to satisfy) · `REFUSE` (do not ingest, reason
given).

**Source dossiers** (all in this scratchpad): the per-paper `*_extraction.md` files, plus the
wave indices `z0_index.md`, `z1_disjointness_verdict.md`, `z2_index.md`, `z3_index.md`,
`k1_kinetic_parameters.md`, `k2_matrix_and_thresholds.md`, and this wave's
`zhou2023_extraction.md`, `Xin2026b_extraction.md`, `Zhang2024_extraction.md`
(+ its **Wave K3 addendum**), `zhu2023_kafirin_binding_extraction.md`.

---

# §0. THE EXECUTIVE SUMMARY — eleven things the build waves must not get wrong

**0. ★ BRANCH FRACTIONS ARE NOT CONSTANTS, AND THIS IS NOW MEASURED RATHER THAN ARGUED.**
A **2× change in precursor loading moves the xylose share of MFT from 15 % to 46 %** — one lab,
one method, one pH, one temperature (Cerny 2007 Table 5; §A.3.6 ii, §B10.1). Wave Z1 §E deduced
that the Hofmann-vs-Cerny disjointness is a *concentration-handling* defect and had to estimate
the effect size from an assumed rate law; **this measures it.** Any fixed thiamine : sugar MFT
split — or any fixed branch fraction anywhere — is falsified by one row pair. Frankel's
hydroperoxide-isomer effect (§A.5) and Cerny's temperature effect (§B3.6) are the same finding on
two other lanes.

**0b. The corpus finally has a norfuraneol activation energy** — **104.9 ± 8.9 / 121.1 ± 8.1 /
122.3 ± 19.5 kJ/mol**, verified three independent ways (§A.3.6 i) — **and it is not a constant**:
across three studies it spans **64–122 kJ/mol** as a function of precursor loading and matrix.
It is an *apparent, lumped* approach-to-plateau rate in an **alkaline (pH 8.4–9.5)** food gel,
80–100 °C, 3 points. **Never call it "the norfuraneol barrier".**

**0c. An alkaline-pH wall now qualifies most of the new trunk and browning constants.** Ajandouz
(pH 8.0 / 9.7) and both Bornhorst papers (pH 8.4–9.5) supply nearly every new Ea, and all three
dossiers refuse rate transfer to pH 5–7 — with Ea transfer licensed **only** for glucose-loss and
amino-loss (±10 %), **explicitly not for browning** (§C.13).


**1. There is no activation energy for thiol consumption, and the corpus now shows there
should not be ONE.** K1 declared the gap; Z2 and K3 turned it from "we lack a number" into
"a single Arrhenius term would be the wrong model". Four papers at four temperatures measure
four *different* dominant sinks (§A.4, §B.7).

**2. Two constants in `arrhenius_params.yml` follow the same forensic signature — a confirmed
Ea bolted to an invented prefactor.** `cysteine_thermolysis` (Z0 finding 2) and
`schiff_condensation` (Z2 finding 3). **Two of two audited.** The remaining `A_value` entries
have not been checked and should be.

**3. Three fabrication-class citation defects are now closed as un-re-pointable.**
`Ea = 129 kJ/mol` for acrylamide (absent from all three Knol papers, Z0 #1 + Z2 #2); the
`52.1 / 72.9 kJ/mol` pair (true values 94.4 ± 11 / 85.1 ± 14, Z2 #4); and the retired
`342 / 200 ppb` sulfur anchor (true values 198 / 121 ppb, `hofmann1998_extraction.md`).

**4. No family-level pH term can pass.** Four independent sign-crossings are now on record
(§B.2). The newest — **Zhou 2023 Figure 3** — has the C0/C1 pyrazines and the C2 pyrazines
moving in **opposite** pH directions **in the same pot, from the same α-dicarbonyl**.

**5. The matrix threshold question is settled and the answer is "no general factor".**
Ratios span **3.4× to 6 714×** with a 1-σ band 27–41× wide, and a large part of the beef
numbers is thermal loss, not perception (K2 §D.2). Ship a lookup table with an explicit
`no_factor_available` state.

**6. Reversible binding cannot deliver matrix odour activity, and this is quantified.**
Ceiling **1.3–3.7× at 4 % protein**, **≤7.6× extrapolated to 10 %** (K2 §B.4), against
threshold shifts of two to four orders. And perception responds *less* than headspace does
(K2 §B.5), so a headspace-calibrated `f_free` *over*-predicts.

**7. The pentose ≫ hexose ordering is matrix-scoped, not universal.** Aqueous 145 °C:
ribose/glucose = 10.4× on MFT (Hofmann 1998 T1). Real 65 %-moisture extruder: **fructose beats
ribose 4.8× and xylose falls below the no-sugar control** (Xin 2026, corroborated by colour).

**8. ★ Single-internal-standard HS-SPME µg/kg are not absolute concentrations, and the corpus can
now prove it from inside the literature.** The same physical extrudate, reported by the same
group in two 2026 papers, differs by **10–23× on three compounds with non-overlapping ranges**
while agreeing to **3 %** on a fourth (§E.2.2, §B9.1). Flag `absolute_concentration: false` on
every semi-quant HS-SPME row — Conti ×2, Xin ×2, Zhang 2024, Zhou 2023.

**9. The corpus's biggest remaining hole is not a parameter, it is a PAIRING.** There is no
same-method matrix-vs-water threshold pair anywhere (K2 §D.2 i); no thiol Ea (K1 §1d); no
MFT-vs-temperature number (Zhang 2024 §7); and no precursor dose–response in any of the three
extrusion papers (K2 §C.4). Four missing *comparisons*, not four missing constants.

---

# §A. MASTER PARAMETER TABLE

## A.1 — TRUNK: condensation, Amadori, browning, organic acids

| parameter | value | units | step | conditions | source anchor | prov. | verdict |
|---|---|---|---|---|---|---|---|
| X (step 1) | **1.6 × 10⁻⁵** | **L mmol⁻¹ min⁻¹** (= 2.67 × 10⁻⁴ L mol⁻¹ s⁻¹) | Glc + Gly → DFG (condensation) | aqueous, T_ref 100 °C | Martins 2005 T2 | [F] | **USE** — ⚠️ the PDF text layer **strips the minus sign from every exponent** in this table; read from the 300 dpi render (Z2) |
| Ea (step 1) | **96.8 ± 2.8** | kJ/mol | " | " | Martins 2005 T2 | [F] | **USE** — this is the repo's `schiff_condensation` Ea = 97.0, **confirmed** |
| A (`schiff_condensation`, shipped) | 1.5 × 10¹¹ | L mol⁻¹ s⁻¹ | " | — | `arrhenius_params.yml` | — | **REFUSE — over-predicts its own cited source by 14.8×.** Implied correct **A ≈ 1.0 × 10¹⁰** (Z2 #3) |
| X (step 9) | **8.1 × 10⁻⁴ ± 1.7 × 10⁻⁵** | L mmol⁻¹ min⁻¹ | **3-DG + Gly → melanoidins** | aqueous, T_ref 100 °C | Martins 2005 T2 | [F] | **USE — the first browning parameter this lane has ever had** (Z2 #15) |
| Ea (step 9) | **95.2 ± 2.3** | kJ/mol | " | " | Martins 2005 T2 | [F] | **USE** |
| ε (melanoidin, Glc/Gly) | **0.64 ± 0.03** | L mmol⁻¹ cm⁻¹ (= **640** L mol⁻¹ cm⁻¹) | browning readout at 470 nm | aqueous | Martins 2005 | [M] | **USE-Q — amine-specific** |
| ε (melanoidin, Glc/**Asn**) | **282** | L mol⁻¹ cm⁻¹ | " | aqueous | Knol 2005 | [M] | **USE-Q — 2.3× from the Gly value. Any single ε across amino acids is wrong by that factor** (Z2 #18) |
| Martins Table 2, all 10 steps | X ± 95 % HPD **and** Ea ± 95 % HPD | as printed | incl. organic acids, sugar isomerisation, methylglyoxal | aqueous, T_ref 100 °C | Martins 2005 T2 | [F] | **USE — the repo's existing audit_flag reproduces the ten rounded Ea and carries NO uncertainties. Z2 supplies every ±** |
| k1 (Glc + Asn → Schiff base) | Ea **57.6 ± 8.0** | kJ/mol | condensation | 120–200 °C | Knol 2005 T1 | [F] | **USE-Q — ⚠️ CONFLICT with Martins step 1 (96.8 ± 2.8), intervals non-overlapping.** Owner question before either is called "the" condensation barrier (Z2 #17) |
| isomerisation | **61 ± 8** | kJ/mol | sugar isomerisation | 120–200 °C | Knol 2010 T2 | [F] | **USE** |
| acetic acid formation | **75 ± 10** | kJ/mol | organic acid | " | Knol 2010 T2 | [F] | **USE** |
| formic acid formation | **84 ± 14** | kJ/mol | organic acid | " | Knol 2010 T2 | [F] | **USE** |
| Knol 2010 Table 2, complete | 7 steps × 5 T + Ea + HPD | as printed | trunk | 120–200 °C | Knol 2010 T2 | [F] | **USE** |
| Schiff/Amadori split | — | — | — | — | Martins 2005 T1 | [M] | **STRUCTURAL — REFUSE THE SPLIT.** The author tested E1 three ways, then removing it entirely "fitted the data equally well"; "E1 is not a rate-determining step". **Carry ONE composite** (Z2 #5) |
| Ea (`beta_elimination_dha`, shipped) | 79.5 | kJ/mol | β-elimination | — | repo | — | **REFUSE — 43.5 kJ/mol below the measured aqueous barrier (123.0 at pH 9, Zheng 1994), and its prefactor is a silent NaN substitution** (Z0 #21) |

### A.1.1 — ★ NEW VIA Z3: the Ajandouz 2008 activation-energy set — 24 own-measurement Ea

Shared conditions for **every** row: **glucose 0.2 M; BSA or casein 5 mg/mL; 0.2 M sodium
phosphate at pH 8.0 or 0.2 M sodium borate at pH 9.7; 60 / 70 / 80 / 100 °C; heating times
pH 8.0 = 240 / 180 / 120 / 40 min, pH 9.7 = 100 / 100 / 100 / 10 min; screw-cap tubes,
ice-cooled; each kinetic point ≥ duplicate; rate constants from INITIAL SLOPES ONLY.**
Paper's own error bar: **mean deviation ≤ 12 %** (§2.2, p. 1245). All **[M]**, Ajandouz 2008
Tables 1–3 and §3.4.

| observable | glucose alone | + casein | + BSA | pH | anchor |
|---|---:|---:|---:|---|---|
| **Ea, glucose disappearance** (kJ/mol) | **93** | **85** | **87** | 8.0 | T1 p. 1246 |
| " | 96 | 93 | 92 | 9.7 | T1 |
| **Ea, free amino-group loss** (kJ/mol, TNBS) | **143** (casein alone) / **116** (BSA alone) | **106** | **102** | 8.0 | T2 p. 1248 |
| " | — | 105 | 92 | 9.7 | T2 |
| **Ea, browning A₄₂₀** (kJ/mol) | **164** | **120** | **130** | 8.0 | T3 p. 1250 |
| " | 126 | 92 | 95 | 9.7 | T3 |
| **Ea, UV A₂₉₄** (kJ/mol) | **152** | **128** | **129** | 8.0 | **§3.4 running text — NOT TABULATED** |
| " | 123 | 107 | 101 | 9.7 | §3.4 text |

**Verdicts.** Glucose-disappearance rows: **USE as a referee for `k_glc_other` (amine-free) and
`k_schiff` (carbonyl side)**. Amino-group-loss rows: **USE as the AMINE-side referee for
`k_schiff`** — the corpus has never had one. Browning rows: **USE — the first A₄₂₀ Ea set in the
corpus**, and a third-lab cross-check on Martins' melanoidin Ea 95.2.

**★ AND THE PARTITION BENCHMARK, which the corpus has nothing else like** (§3.4, p. 1250, **[M]**):
| quantity | value | conditions |
|---|---:|---|
| **caramelisation share of A₂₉₄** | **25–80 %** of the glucose–protein system's UV absorbance | pH 8.0 and 9.7, 60–100 °C; **"percentage increases with temperature"** |
| **caramelisation share of A₄₂₀** | **7–55 %** of the browning | " |
| caramelisation share of A₂₉₄ | 40–62 % | **[C]** — 100 °C, **pH 4.0–7.0**, glucose + lysine (Ajandouz & Puigserver 1999). ⚠️ second-hand, **but this is the row at the REPO's own pH range** |
| caramelisation share of A₄₂₀ | 10–36 % | **[C]**, same |

⇒ **the amine-independent lane carries a quarter to four-fifths of the UV signal and up to half
the colour. This sizes Wave U's `structural_zero` directly.**

⚠️ **THE ALKALINE-pH WALL applies to this whole block — see §C.13.** Ajandouz licenses Ea
transfer to pH 5–7 **only for glucose-loss and amino-loss (±10 %)**, and **explicitly NOT for
browning**, whose Ea fall **15–29 %** between pH 8.0 and 9.7.

⚠️ Ajandouz Tables 1–3 also carry **~14 second-hand [C] rows** (Warmbier 1976, Morales & van
Boekel 1998, Lea & Hannan 1949, Jokinen 1976, Thompson 1976, Carpenter 1962, Malec 2002,
Ben-Gara 1972, Fabriani 1972, Tsao 1978, Brands & van Boekel 2001/2002, Bluestein 1975). Two are
**recomputed by Ajandouz from the cited authors' data** (footnote a) — a double-laundering risk.
**The one that matters most: Ea(amine loss, ribose + BSA) = 130 kJ/mol, 85–145 °C** (Carpenter
1962) — **the closest system to the repo's ribose lane in the whole table**, and its a_w/%-H₂O
cell prints an ambiguous "14.1". **[RETRIEVE] before use.**

## A.2 — ACRYLAMIDE / SAFETY

The best-parameterised block in the corpus.

| parameter | value | units | step | conditions | source anchor | prov. | verdict |
|---|---|---|---|---|---|---|---|
| k_Fref (control) | **0.451 ± 0.023 × 10⁻³** | min⁻¹ | Reactant → acrylamide | **0.01 M Asn + 0.01 M Glc**, 0.05 M citrate **pH 6**, closed inox, 140–200 °C, **T_ref 160 °C**, non-isothermal corrected, n = 35 | Claeys 2005 T2 p. 1528 | [F], ± = SE | **USE** — self-validates against its own Fig. 2 (plateau 2886 vs ~2800 ppb; t_max 49.8 vs 40–50 min) |
| k_Eref (control) | **111.1 ± 8.9 × 10⁻³** | min⁻¹ | acrylamide → DP | " | Claeys 2005 T2 | [F] | **USE-Q — k_F and k_E are 0.94–0.98 CORRELATED (authors state it). Ingest the RATIO or carry the correlation** |
| Ea_F / Ea_E (control) | **168.25 ± 3.80 / 167.21 ± 4.30** | kJ/mol | " | " | Claeys 2005 T2 | [F] | **USE** |
| Claeys competitor set | Gln, Cys, Lys, Ala full rows | as printed | " | + 0.01 M competitor | Claeys 2005 T2 | [F] | **USE** — k_E/k_F ratio Gln 167 < Ala 222 < control 246 < Lys 477 < **Cys 536** [D] |
| **k_E2ref** | **49.36 ± 1.18** | **M⁻¹ min⁻¹** | **acrylamide + cysteine → adduct, SECOND ORDER** | a_w 0.92 freeze-dried powder, ~3 M, T_ref 160 °C | De Vleeschouwer 2009 II T3 p. 542 | [F], ± = 95 % HPD | **★ USE — the tightest parameter in the corpus (2.4 % RSE). Generalises to any thiol/Michael-acceptor pair** |
| **Ea_E2** | **51.3 ± 1.5** | kJ/mol | " | " | De Vleeschouwer 2009 II T3 | [F] | **★ USE — HALF of Ea_E (113.2). The two channels CROSS OVER in temperature** |
| k_INTgref | **1.70 ± 1.05** | **M⁻¹ min⁻¹** | Asn + Glc → Int1, **SECOND ORDER** | a_w 0.92, ~3 M, T_ref 160 °C | De Vleeschouwer 2009 I T3 p. 124 | [F] | **USE — the only genuinely bimolecular Maillard-initiation constant in the corpus** |
| k_Fref (glucose) | **3.57 ± 1.38 × 10⁻³** | min⁻¹ | Int1 → acrylamide | " | De Vleeschouwer 2009 I T3 | [F] | **USE** — ⚠️ printed unit "(10⁻³ mm⁻¹)" is a typo for min⁻¹ |
| k_Eref (glucose) | **0.10 ± 0.04** | min⁻¹ | acrylamide → DP | " | De Vleeschouwer 2009 I T3 | [F] | **USE** |
| k_Aspref (glucose) | **26.43 ± 5.76 × 10⁻³** | min⁻¹ | Asn → Asp deamidation | " | De Vleeschouwer 2009 I T3 | [F] | **USE** |
| Ea_F / Ea_E / Ea_INTg / Ea_Asp (glucose) | **159.2 ± 29.5 / 113.2 ± 32.3 / 117.5 ± 25.2 / 105.4 ± 10.6** | kJ/mol | " | " | De Vleeschouwer 2009 I T3 | [F] | **USE** — Ea_Asp is remarkably sugar-independent (105–109 across all three sugars) |
| k_HYref / k_Iref (sucrose) | 0.47 ± 0.09 / 0.50 ± 0.20 × 10⁻³ | min⁻¹ | sucrose hydrolysis / Glc→Frc | " | De Vleeschouwer 2009 I T3 | [F] | **USE** |
| Ea_hy / Ea_I | 140.6 ± 8.5 / 113.2 ± 13.8 | kJ/mol | " | " | De Vleeschouwer 2009 I T3 | [F] | **USE** |
| k_Yref / Ea_Y | 0.35 ± 0.01 min⁻¹ / **110.5 ± 8.5** kJ/mol | — | **unspecified** Cys sink | " | De Vleeschouwer 2009 II T3 | [F] | **USE-Q — the product is not identified** |
| *k_M, k_B, k_C, k_X and their Ea* | *various* | — | — | — | De Vleeschouwer 2009 I T3, **italicised** | [F] | **REFUSE — the authors mark them "NO PHYSICAL MEANING". Includes Ea_Cg = −6.7 (negative) and Ea_X = 668.9 (absurd)** |
| **k6** | **7.96 / 28.1 / 88.1 / 250 / 650 × 10⁻³** | min⁻¹ at 120/140/160/180/200 °C | **Acrylamide → Product X, FIRST ORDER** | aqueous Asn/Glc | **Knol 2005 T1** | [F] | **★ USE — this is the repo's `− ke·[A]` term at `src/reaction_templates.py:2029`. Re-point the citation from "Knol 2009"/"Knol 2010" to Knol 2005** (Z2 #1) |
| Ea (k6) | **85.1 ± 14** | kJ/mol | " | " | Knol 2005 T1 | [F] | **USE-Q — carry the author's caveat: "the model was NOT restrained by experimental data for the products formed in the degradation reaction"** |
| Ea (k4, formation) | **94.4 ± 11** | kJ/mol | acrylamide formation | " | Knol 2005 T1 | [F] | **USE** — this and 85.1 are the true pair that `safety_reference_payloads.json` entries[27] mis-states as 52.1 / 72.9 |
| Knol 2005 full set | **30 rate constants + 6 Ea, ±95 % HPD, 120–200 °C** | as printed | trunk + acrylamide | aqueous | Knol 2005 T1 | [F] | **USE** — concentration–time data are **figures only**; flag `digitised_from_figure` |
| **`Ea = 129 kJ/mol`** (shipped) | 129 | kJ/mol | acrylamide | — | `src/barrier_constants.py:274, :493`, `src/safety.py:790–795` | — | **★ REFUSE — FABRICATION-CLASS, TRIPLE-CONFIRMED AND UN-RE-POINTABLE.** Knol 2005: max 102 (upper-95 % 116); Knol 2009: **no Ea at all**; Knol 2010: max 93 ± 12 |
| `A_f = 1.6 × 10¹³` (shipped) | — | — | acrylamide | — | repo | — | **REFUSE — fit-circularity signature: `A_f·exp(−129000/RT)` at 160 °C = 4.4 × 10⁻³ vs Knol's k2 = 4.5 × 10⁻³, i.e. back-solved to Knol's own T_av with the wrong Ea. The ~2000× separating them in shipped code is `_acrylamide_ph_factor(5.5)`** (Z0 #13) |
| acrylamide **degradation** parameters generally | — | — | — | — | Knol 2005 + 2009 + 2010 | — | **REFUSE AS CONSTANTS — UNIDENTIFIABLE, three times over.** Knol 2005: unconstrained by data. Knol 2009: **every** degradation parameter has SD ≥ estimate (`tc2 = 0.60 ± 52`, `k2 = 3.5 ± 8`, `τ = 22 ± 15`). Knol 2010: the Ea went negative and was deleted (Z2 #12) |
| Quan 2020 Table 1, all 88 constants | — | **NO UNITS PRINTED** | 11 steps | 0.1 M phosphate pH 7.0, 100 & 130 °C, 7 points 0–21 min, n = 3, SIDA | Quan 2020 T1 | [F] | **REFUSE — four disqualifying defects: no units; no rate laws (SI Fig. 1 absent, so every ORDER is unknown); mislabelled continuation headers; no ± despite a footnote claiming them. Authors: "specific activation energies cannot be estimated"** |
| acrylamide yield, real system | ~3 % of initial asparagine ⇒ **~2.1 × 10⁵ ppb** at 0.1 M, **~4 × 10⁵** scaled to 0.2 M | ppb | end-to-end | — | Knol 2010 | [M]+[D] | **USE as the arbiter of the repo's 480× two-lane contradiction: `predict_acrylamide` 544 870 ppb is within ~1.3×; the network's 1143 ppb is ~350–380× LOW** (Z0 #3) |

## A.3 — SULFUR: THIOL AND HETEROCYCLE FORMATION

### A.3.1 The absolute aqueous anchors (SIDA — the highest-quality tier)

Shared conditions for every Hofmann 1998 row: **cysteine 3.3 mmol + carbohydrate 10.0 mmol in
100 mL, 0.5 mol/L phosphate, pH 3.0 / 5.0 / 7.0, 20 min at 145 °C** (= cys 33 mM, sugar 100 mM,
**1:3 cys:sugar**). SIDA with [²H₂]-FFT and [²H₃]-MFT; triplicates "differed by not more than
10 %". `mol% = mol product / mol limiting precursor × 100`, basis verified.

| system | FFT µg/100 mL | MFT µg/100 mL | ⇒ FFT ppb | ⇒ MFT ppb | prov. | verdict |
|---|---:|---:|---:|---:|---|---|
| **ribose, pH 5.0** | **12.1** | **19.8** | **121** | **198** | [M] | **★ USE — the sulfur branch's primary absolute anchor** |
| xylose, pH 5.0 | 9.6 | 14.3 | 96 | 143 | [M] | **USE** |
| fructose, pH 5.0 | 3.2 | 2.5 | 32 | 25 | [M] | **USE** |
| glucose, pH 5.0 | 2.8 | 1.9 | 28 | 19 | [M] | **USE** |
| glucose-6-P / rhamnose / maltose, pH 5.0 | 0.9 / 0.8 / 0.6 | 0.6 / 0.8 / 0.3 | 9 / 8 / 6 | 6 / 8 / 3 | [M] | **USE** |
| ribose, pH 3.0 / 7.0 | 22.9 / 1.2 | 55.3 / 2.5 | 229 / 12 | 553 / 25 | [M] | **USE** |
| glucose, pH 3.0 / 7.0 | 0.7 / 0.6 | 0.3 / 0.4 | 7 / 6 | 3 / 4 | [M] | **USE** |
| rhamnose, pH 3.0 / 7.0 | 0.2 / 0.1 | 0.1 / 0.1 | 2 / 1 | 1 / 1 | [M] | **USE** |
| ribose / glucose / rhamnose, **dry 180 °C / 6 min** on silica | 97.2 / 1.4 / 0.4 | 25.1 / 4.2 / 3.1 | — | — | [M] | **USE-Q — DRY system (3.0 g silica + 300 µL of 2 M buffer). NOT comparable to the aqueous rows** |
| **`342 / 200 ppb`** (retired) | — | — | — | — | — | **REFUSE — appears NOWHERE in the paper. Fabricated/real = 1.73× and 1.65× high.** The repo's conditions were also wrong on all four of T, t, stoichiometry and buffer |

### A.3.2 Step-level (fed-precursor) yields — the first step-level constraints the repo ever had

All Hofmann 1998, **50 mL, 1 mmol precursor + 1 mmol co-reactant, pH 5.0, 145 °C / 20 min**
unless stated.

| system | product | µg | mol % | prov. | verdict |
|---|---|---:|---:|---|---|
| ribose + H₂S | FFT | 9.2 | 0.008 | [M] | USE |
| 3-deoxyribosulose + H₂S | FFT | 78.6 | 0.08 | [M] | USE |
| **furan-2-aldehyde + H₂S** | FFT | **550.8** | **0.48** | [M] | **USE — furfural is 60× ribose and 7× the 3-deoxyosone as an FFT source** |
| ribose + H₂S | MFT | 15.1 | 0.01 | [M] | USE |
| **norfuraneol + H₂S** | MFT | **211.2** | **0.19** | [M] | **USE — 14× ribose** |
| norfuraneol + cysteine | MFT | 50.8 | 0.05 | [M] | USE |
| **hydroxyacetaldehyde + mercapto-2-propanone (C2+C3)** | MFT | **268.1** | **0.24** | [M] | **★ USE — the single most effective MFT route measured** |
| thiamin | MFT | 8.2 | 0.01 | [M] | **USE — thiamin is 30× weaker than the C2+C3 route** |
| 2-oxopropanal + H₂S **1:1** | mercapto-2-propanone | 1650 | 1.8 | [M] | **USE** |
| 2-oxopropanal + H₂S **1:2** | mercapto-2-propanone | 3600 | **4.0** | [M] | **★ USE — 2× H₂S gives 2.2× product. H₂S is rate-limiting and SUPER-LINEAR** |
| C2+C3, pH 3.0 / 5.0 / 7.0 | FFT | 26.1 / 40.5 / 58.2 | 0.02 / 0.04 / 0.05 | [M] | USE |
| C2+C3, pH 3.0 / 5.0 / 7.0 | MFT | 15.5 / 268.1 / 311.5 | 0.01 / 0.23 / 0.27 | [M] | USE |
| C2+C3, **dry 180 °C / 6 min** | MFT | **1553.9** | **1.39** | [M] | **USE-Q — this is the abstract's headline "1.4 mol%" and it is a DRY system** |
| C2+C3, pH 3/5/7 | FA / NF | 101.2, 153.1 / 364.5, 885.1 / 1443.1, 3610.5 | 0.11, 0.13 / 0.40, 0.78 / 1.60, 3.18 | [M] | USE |
| glucose + cysteine (100 mL), pH 3/5/7 | mercapto-2-propanone | 11.4 / 59.6 / 26.5 | — | [M] | **USE — peaked at pH 5** |
| ribose (100 mL, pH 5.0) | FA / FFT / **NF** / MFT | 67.5 / 12.1 / **54 530.0** / 19.8 | — | [M] | **★ USE — ribose makes 54 530 µg NF against 19.8 µg MFT, a factor of ~2 750, while glucose's ratio is ~10. The authors' own conclusion: "MFT formation might not run exclusively via NF as the key intermediate"** |
| **NF + cysteine, pH 4.5, 140 °C / 60 min** | MFT | 15, 15 µg/10 mg NF | **0.150** | [M], DHS RF = 1 | **USE-Q — semi-quantitative (response factors ASSUMED 1)** — Whitfield 1999 T1 |
| **NF + H₂S (1:2), pH 4.5, 140 °C / 60 min** | MFT | — | **0.120** | [M], DHS | **USE-Q — two labs, two methods, agreement within ~1.6× on the H₂S channel** |
| **NF + cysteine, pH 6.5, 140 °C / 60 min** | MFT | **nd (< 0.1 = LOD)** | **< 0.0010** | [M] | **★ USE — ≥150× collapse over two pH units within one lab** — Whitfield 2001 T1 |

### A.3.3 The fed-Amadori system — NEW THIS WAVE (Zhou 2023)

Conditions: **purified N-(1-deoxy-D-xylulos-1-yl)-alanine 20 mmol/L ± L-cysteine 20 mmol/L
(1:1), DEIONIZED WATER — UNBUFFERED, initial pH 6 / 7 / 8, 120 °C, 60 min**, HS-SPME +
external calibration vs 1,2-dichlorobenzene, n = 3. Full 22-compound table in
`zhou2023_extraction.md` §2.

| quantity | pH 6 | pH 7 | pH 8 | prov. | verdict |
|---|---:|---:|---:|---|---|
| **MFT, µg/L** (= ppb) | **696.99 ± 17.44** | **1588.57 ± 21.24** | **525.62 ± 16.88** | [M] | **USE-Q — NOT SIDA; not commensurable with A.3.1 in absolute terms** |
| **FFT, µg/L** | 813.65 ± 13.37 | 757.965 ± 13.03 | 325.22 ± 14.61 | [M] | **USE-Q** |
| bis(2-methyl-3-furyl) disulfide, µg/L | 59.70 ± 2.75 | 102.59 ± 4.17 | 50.07 ± 2.15 | [M] | **USE-Q — ⚠️ its calibration curve zeroes at x = 34.3 µg/L (SI T1 intercept −0.7126), so 50.07 is only 1.5× above a pseudo-LOD** |
| bis(2-furfuryl) disulfide, µg/L | 35.31 ± 1.63 | 17.68 ± 0.76 | **ND** | [M] | USE-Q — "ND" may mean "< 4.1 µg/L" |
| 2-acetylthiazole, µg/L | 2.43 ± 0.08 | 11.70 ± 2.14 | **582.34 ± 14.57** | [M] | **USE** |
| thiazole, µg/L | 5.27 ± 0.22 | 10.90 ± 0.99 | 107.38 ± 8.50 | [M] | USE |
| pyrazine / methylpyrazine, µg/L (ARP+Cys) | ND / ND | 6.06 / 1.03 | **177.05 / 99.00** | [M] | USE |
| 2-furfural, µg/L (ARP **alone**) | 1606.92 ± 41.62 | 1339.37 ± 83.04 | 436.63 ± 19.93 | [M] | USE |
| 2,3-butanedione, µg/L (ARP alone → ARP+Cys) | 8.93 → **ND** | 33.63 → **ND** | 126.18 → 20.00 | [M] | **USE — Cys suppresses it completely at pH 6–7 and by 84 % at pH 8** |
| **MFT mol% of ARP** | 0.0305 % | **0.0696 %** | 0.0230 % | [D] | **USE-Q — same order as Hofmann's fed-precursor mol% (0.01–0.24 %)** |
| **MFT / FFT** (mass = molar) | **0.857** | **2.096** | **1.616** | [D] | **★ RATIO-ONLY, USE — cross-validates against Hofmann 1998 T1 pH 5 (1.636). Two labs, two methods, two feedstocks** |
| **dimer / MFT, thiol-equivalents** | **8.6 %** | **6.5 %** | **9.6 %** | [D] | **★ USE — near-invariant in pH while [MFT] swings 3.0×** |
| dimer / FFT, thiol-equivalents | 4.4 % | 2.4 % | 0 (< LOD) | [D] | USE-Q |
| **final pH** (ARP alone) | **3.99** | **4.13** | **4.99** | [M] **[fig]** | **★ USE — MUST be ingested alongside the table or the table is misleading** |
| **final pH** (ARP + Cys) | **3.22** | **3.42** | **5.07** | [M] **[fig]** | **★ USE** |
| 2,5-DMP : 2,6-DMP at 90 min (Cys+MGO) | 8.4 | 8.4 | 9.2 | [D] | **★ RATIO-ONLY, USE — a pH-insensitive regiochemical branch ratio from a single α-dicarbonyl** |
| 2-acetylthiazole at 60 min: ARP-fed vs MGO-fed, pH 8 | — | — | **582 vs 665 µg/L (14 %)** | [M]+[M] | **★ USE — a step-level, two-system constraint on the ARP → MGO flux** |

**Cys + MGO fed system** (both **20 mmol/L**, 1:1, unbuffered, pH 6/7/8, **120 °C, 0/30/60/90
min**): the full 8-panel × 3-pH × 4-time grid is transcribed in `zhou2023_extraction.md` §4.
All values **[M] [fig, 300 dpi]** with per-panel uncertainties stated. **This is the only
α-dicarbonyl-fed, time-resolved, pH-resolved dataset in the corpus.**

### A.3.4 Precursor conversion — a constraint class the repo has none of

| quantity | value | conditions | source | prov. | verdict |
|---|---|---|---|---|---|
| **cysteine consumed** | **33 → 15 mM (55 %)** | Hofmann process flavour, **130 °C / 20 min** | van Seeventer 2001 | [M] | **★ USE — scores REACTANT conversion, independent of every downstream branching ratio** (Z2 #14) |
| **ribose consumed** | **100 → 25 mM (75 %)** | " | van Seeventer 2001 | [M] | **★ USE** |

### A.3.5 Sulfur budget

| quantity | value | source | prov. | verdict |
|---|---|---|---|---|
| H₂S supplied by fed-precursor experiments | **1–2 mmol** | Hofmann 1998 T3/T4/T7 protocols | [M] | — |
| H₂S from cysteine thermolysis in situ | **≈ 0.17 mmol** at 145 °C / 20 min | Zheng 1994 + Hofmann 1998 | [D] | **STRUCTURAL — a ≈6–12× SURPLUS in the fed experiments. This is the quantitative mechanism behind the Hofmann/Cerny disjointness** (Z0 #8) |
| NF concentration, **fed** vs **in situ** | 20 mM vs **4.78 mM** (= 4.2×) | Hofmann 1998 T4 vs T5 | [M]+[D] | **STRUCTURAL** (Z1 §E) |

### A.3.6 — ★ NEW VIA Z3: the norfuraneol Ea, the thiamine route, and a furfural sink Ea

#### (i) The corpus's FIRST norfuraneol activation energy — Bornhorst 2017b, LWT `S0023643817302384`

Conditions: **mashed-potato model gel** (15 g instant potato flakes + 0.5 g low-acyl gellan gum
+ 0.13 g CaCl₂ + 0–2 g D-ribose + 0–2 g L-lysine + 80.37–84.37 g DDI water per 100 g), firm gel,
aluminium test cells, **80 / 90 / 100 °C** (the 90 °C row is **imported** from the 2017 companion,
not re-measured), **come-up time 1.75 min EXCLUDED from every reported time**, triplicate.
Model: `C = C∞ − (C∞ − C₀)·exp(−k·t)` — **an approach-to-plateau FORMATION law with NO
destruction term**, M-2₀ fixed at zero. **M-2 = 4-hydroxy-5-methyl-3(2H)-furanone = NORFURANEOL.**
**★ pH 8.4–9.5 throughout; CaCl₂ present in every row.**

| formula (ribose g / lysine g per 100 g) | **Ea, M-2 formation** | **z-value** | k at 80 / 90* / 100 °C (10⁻³ min⁻¹) | M-2∞ at 80 / 90* / 100 °C (mg/g) | prov. |
|---|---:|---:|---|---|---|
| **1_R 0.5_L** | **121.1 ± 8.1 kJ/mol** (28.9 kcal/mol) | 20.8 ± 1.4 °C | 1.5 ± 0.5 / 5.1 ± 0.9 / 13.2 ± 1.1 | 0.57 ± 0.17 / 0.47 ± 0.06 / 0.43 ± 0.02 | [F] |
| **1_R 1_L** | **122.3 ± 19.5 kJ/mol** (29.2 kcal/mol) | 20.6 ± 3.4 °C | 3.2 ± 1.8 / 7.4 ± 1.8 / 30.0 ± 3.7 | 0.18 ± 0.08 / 0.24 ± 0.03 / 0.16 ± 0.01 | [F] |
| **2_R 2_L** | **104.9 ± 8.9 kJ/mol** (25.1 kcal/mol) | 24.0 ± 2.0 °C | 6.0 ± 1.8 / 18.3 ± 4.1 / 40.3 ± 3.7 | 0.15 ± 0.03 / 0.10 ± 0.01 / 0.14 ± 0.01 | [F] |

**Verdict: ★ USE — this is the repo's FIRST norfuraneol Ea, and it is verified three independent
ways** (the Z3 analyst's Arrhenius refit reproduces 119.4 / 123.0 / 104.7; the `D = 2.303/k`
checks reproduce every printed D-value to ≤ 3 %; the z-values reproduce from
`z = 2.303·R·T_ref²/Ea` at a **T_ref = 363.15 K that the analyst had to recover from the
arithmetic — the paper never states it**).

⚠️ **FOUR MANDATORY QUALIFICATIONS.**
1. **It is an APPARENT, LUMPED rate of approach to a plateau of NET accumulation — not an
   elementary step.** The model has no destruction term at all.
2. **pH 8.4–9.5.** See §C.13. *(The un-supplemented matrices are pH 5.2 / 6.0 / 6.1; adding
   ribose + lysine drives them to 8.4–9.9. Every row is alkaline **because of the precursors**.)*
3. **Only 3 temperatures per Arrhenius fit, and one of them is imported from another paper.**
4. **CaCl₂ is present in every row**, and the M-2 data for the gellan matrix were **withheld by
   the authors** ("excluded … due to low M-2 concentrations", Table 1 header).

**The 90 °C companion (Bornhorst 2017, `S0023643816305825`) supplies the matrix comparison**
(all **[F]**, 90 °C): M-2∞ **egg white 0.54 / 0.28 / 0.14** vs **mashed potato 0.47 / 0.24 /
0.10 mg per g sample** for 1_R0.5_L / 1_R1_L / 2_R2_L; k **egg white 19.9 ± 0.7 / 21.7 ± 1.8 /
22.0 ± 2.9** vs **mashed potato 5.1 ± 0.9 / 7.4 ± 1.8 / 18.3 ± 4.1 × 10⁻³ min⁻¹**.
**★ A structural-zero benchmark comes free: M-2 = 0 in all three matrices with no added
precursors** (§3.1). Full CIELAB L\*/a\* blocks (with their own Ea and z) are in `z3_index.md` §1C/1D.

#### (ii) ★ THE THIAMINE-vs-SUGAR MFT BRANCHING RATIO — measured, with BOTH single-route controls

Cerny & Briffod 2007 (`10.1021/jf062874w`) and Cerny 2007b (`10.1016/j.lwt.2006.09.008`).
Conditions: **cysteine : thiamine : xylose = 1 : 1 : 3**, K-phosphate **0.5 mol/L**, **145 °C /
20 min**, Teflon vials. Cerny 2007b system "A" = **0.0999 M cysteine / 0.0999 M thiamine /
0.2997 M xylose, pH 5.00**, [¹³C₅]xylose at 98 % enrichment.

| system | MFT isotopomer split (% unlabelled = **thiamine** : % labelled = **xylose**) | anchor | prov. | verdict |
|---|---|---|---|---|
| **full ternary** (cys + thiamine + xylose) | **54 : 46** | Cerny 2007b Table 3 | [M] | **★ USE — a MEASURED two-route branching ratio, not a guess** |
| **NO CYSTEINE** (xylose + thiamine) | **> 99 : < 1** | Cerny 2007b Table 5 | [M] | **★ USE — single-route control. MFT is >99 % thiamine-derived when cysteine is absent** |
| **NO THIAMINE** (xylose + cysteine) | **< 5 : > 95** | Cerny 2007b Table 6 | [M] | **★ USE — the other single-route control. >95 % xylose-derived** |
| ternary at **2× concentration** (xylose 0.3 / cys 0.1 / thiamine 0.1 M) | **54 : 46** | Cerny 2007 Table 5, "A" | [M] | **★★ USE** |
| ternary at **1× concentration** (0.15 / 0.05 / 0.05 M) | **85 : 15** | Cerny 2007 Table 5, "B" | [M] | **★★ USE** |

> **★★ THE CONCENTRATION → ROUTE-MIX COEFFICIENT. This is the measurement Wave Z1 §E asked for
> and could not find.** A **2× change in precursor concentration moves the xylose share of MFT
> from 15 % to 46 % — a 3.1× change in the branch fraction.** Z1's §E(i) argued that the
> Hofmann/Cerny disjointness is explained by concentration handling rather than by a barrier, and
> estimated a 13–42× effect from a 4.2× concentration ratio *on the assumption* that the step is
> first order in each reagent. **Here the effect is measured directly, in one lab, one method,
> one pH, one temperature.** The branch fractions are **not constants**; they are functions of
> concentration. **Any model that carries a fixed thiamine:sugar MFT split is falsified by this
> one row pair.**

⚠️ **INGEST THE TABLES, NOT THE PROSE.** Both Cerny papers' running text says "54 % labelled"
where the tables say **46 % labelled / 54 % unlabelled**. The prose error repeats in **both**
papers. See §F #15.

**The pH series that comes with it** (Cerny 2007 Table 2, pH 4.0 / 5.0 / 5.5 / 6.0 / 7.0,
145 °C / 20 min, **GC-TIC integrated peak areas ×10⁶ — ARBITRARY UNITS, explicitly NOT corrected
for MS response factors or partition coefficients; NEVER a concentration**):

| compound | pH 4.0 | 5.0 | 5.5 | 6.0 | 7.0 | shape |
|---|---:|---:|---:|---:|---:|---|
| **MFT** | 415 | 336 | 371 | **539** | 391 | **★ MAXIMUM AT pH 6.0** |
| **FFT** | **431** | 368 | 364 | 185 | **0** | monotone down; **falls to exact zero** |
| **2-furaldehyde** | **208** | 158 | 165 | **0** | **0** | **falls to exact zero at pH ≥ 6** |
| 3-mercapto-2-pentanone | 145 | 141 | 147 | 62 | 12 | monotone down |
| 4,5-dihydro-2-methyl-3(2H)-thiophenone | 0 | 0 | 0 | **372** | **470** | **THRESHOLD SWITCH — on only above pH 5.5** |
| 3-acetyl-1,2-dithiolane | 0 | 0 | 0 | 117 | 240 | threshold switch |
| 5-(2-hydroxyethyl)-4-methylthiazole | **525** | 54 | 43 | 161 | 80 | strongly non-monotone |

**And the isotope splits across that pH ladder (Cerny 2007 Table 4, [M], HARD — fractions are
response-factor-immune):** MFT is **75 / 85 / 90 / 80 % thiamine-derived** at pH 4 / 5 / 6 / 7,
while **FFT is > 99 % XYLOSE-derived at every pH**, and 2-furaldehyde is **> 99 % xylose**.
**2-mercapto-3-pentanone is 94–>95 % XYLOSE while its isomer 3-mercapto-2-pentanone is
77–90 % THIAMINE** — an isomer-split diagnostic.

**[C], and the only absolute thiamine-MFT numbers available** (Cerny 2007 p. 1556, citing
Zeiler ref 26, a thesis): **thiamine-derived MFT 45 µg/L vs ribose/cysteine MFT 3.4 µg/L
(13×)**, both at 120 mmol/L, pH 5.7; and **~1500× in favour of thiamine "when the reagent
concentrations are as low as in meat"** — ⚠️ **that concentration is never numerically
specified; directional only.**

#### (iii) ★ NEW VIA Z3: a real, fitted Ea for a FURFURAL sink — Yaghmur 2005

Conditions (the **water** arm, the only one comparable to our systems): **cysteine 6 mmol +
furfural 3 mmol in 10 g phosphate buffer 0.5 M, pH 5.0** (⇒ cys : furfural = **2 : 1**;
≈ 600 / 300 mmol per kg — *analyst conversion; the paper never states molarity*), **40 / 50 /
65 / 70 °C**, GC-FID, naphthalene internal standard, three replicates, ±1.0 % relative error.

| parameter | value | units | anchor | prov. | verdict |
|---|---:|---|---|---|---|
| k_obs, **water** | **3.2 / 6.7 / 12.6 / 16.0** | **h⁻¹ × 10²** (0.032–0.160 h⁻¹) at 40 / 50 / 65 / 70 °C | Table 1 p. 233 | [F] | **USE** — r² 0.97 / 0.96 / 0.94 / 0.95 |
| **Ea, water** | **46.50 ± 1.0** | **kJ mol⁻¹** (11.11 kcal/mol) | §3.7 p. 232 | [F] | **★ USE — but ONLY under the step name "disappearance of furfural in the presence of cysteine". NEVER as "furfural → FFT"** |
| Arrhenius A, and the Arrhenius r² | **NOT REPORTED** | — | §3.7 | [M] | **UNAVAILABLE** |
| **FFT conversion, water** | **< 1** | **% of furfural charged** at 65 °C / 15 h | §3.1 p. 226 | [M] | **USE — an upper bound** |
| furfural conversion, water, 15 h | 84.9 | % | `1 − exp(−0.126 × 15)` | [D] | USE if labelled derived |
| **⇒ FFT share of the consumed furfural** | **≲ 1.2** | **% of the furfural flux** | derivation | [D] | **★ USE — THIS IS WHY the 46.50 Ea is a lump: the FFT branch is a MINOR channel off a large unidentified sink** |
| t½ of furfural, water, 65 °C | 5.5 | h | k_obs 0.126 h⁻¹ | [D] | USE if labelled derived |
| order in cysteine | **linear, r² = 0.98** over **0–12 mmol** cysteine | — | Fig. 13 | [F] | **USE (structural: FIRST ORDER in cysteine)** ⚠️ **microemulsion arm only** |

⚠️ **Verification: the Z3 analyst refit Table 1 and got 46.47 kJ/mol against the printed 46.50 —
the table and the printed Ea agree to better than 1 kJ/mol.** ✔
⚠️ **REFUSE the other two arms as comparators:** the O/W microemulsion (Ea 32.3) is a surfactant
interface, and the arm labelled "aqueous solution" is **water / propylene glycol 1 : 1 w/w — a
50 % w/w polyol, not water** (Ea 56.72; **and its column has NO r² in the printed table**).
⚠️ The microemulsion components are added to the aqueous reference **after** heating, so the
water arm is not a bare-water control end to end.

#### (iv) 2-acetyl-2-thiazoline, Sotolon, and a required oxidant — Hofmann & Schieberle 1996b

| parameter | value | units | conditions | anchor | prov. | verdict |
|---|---:|---|---|---|---|---|
| **Sotolon** from butane-2,3-dione + hydroxyacetaldehyde | **13.5 / 764.7 / 273.1** | µg | 1.0 mmol each in **50 mL 0.5 M Na/K phosphate, autoclave 145 °C / 20 min**, pH **3.0 / 5.0 / 7.0** | Table 5 expt 1 p. 186 | [M] | **★ USE — the first Sotolon anchor in the corpus.** Peaked at pH 5 |
| Sotolon, **dry-heat** | **131.5** | µg | silica gel 3.0 g + 20.4 mg K₂HPO₄ + 0.3 mL water, **5 min at 180 °C**, pH 5.0 only | Table 5 expt 2 | [M] | USE-Q |
| Sotolon mol % | 0.011 / 0.597 / 0.213 / 0.103 | mol % | pH 3 / 5 / 7 aqueous, then dry pH 5 | derivation | [D] | USE if labelled derived |
| AT (2-acetyl-2-thiazoline) yield, **air** | **9.8** | µg from 190 µg HDT | 20 min in 10 mL water at **80 °C** | Table 4 expt 1 | [M] | USE |
| AT yield, **+ Cu²⁺ 0.01 µmol** | **20.1** | µg | as above | Table 4 expt 2 | [M] | USE |
| AT yield, **+ higher Cu²⁺** | **38.4** | µg | ⚠️ dose printed as "0.01" for both rows in the text layer; the render is ambiguous — **ingest as "a higher Cu²⁺ level"** | Table 4 expt 3 | [M] | **VALUE usable, DOSE UNCONFIRMED** |
| **AT yield, ARGON** | **1.7** | µg | degassed + argon-flushed | Table 4 expt 4 | [M] | **★ USE — the corpus's FIRST required-oxidant measurement. Removing O₂ costs 5.8× of the yield; Cu²⁺ buys 2.1–3.9×** |
| ATD / HDT / AT after charging **ATD** | 18.1 / 76.4 / **51.5** | µg of 150 µg charged | refluxed 20 min in tap water (10 mL), 100 °C | Table 3 p. 185 | [M] | **USE — mass balance closes at 97.3 %** |
| ATD / HDT / AT after charging **HDT** | 22.6 / 120.1 / **10.9** | µg | as above | Table 3 | [M] | **USE — closes at 102.4 %** |
| **Ea, AT formation** | **≳ 97** (a LOWER BOUND) | kJ/mol (23.2 kcal/mol) | `ln Y(t = 1 min)` vs 1/T, 5 points, **50–100 °C**; the paper prints **no rate constant** | derivation from Fig. 2 | **[D]** | **USE ONLY labelled "Wave Z3 derivation from a figure read-off"; biased low** |
| **Ea, AT net decay** | **≳ 166** (a LOWER BOUND) | kJ/mol (39.7 kcal/mol) | **2 points only**; 100 °C 16.1→1.75 µg over 10→60 min; 90 °C 10.3→7.6 over 30→60 min | derivation | **[D]** | **order-of-magnitude only** |
| **Ea(destruction) − Ea(formation)** | **≳ 70** | kJ/mol | both bias directions widen the gap | derivation | [D] | **★ USE AS AN INEQUALITY, NOT A PARAMETER — see §B10.4** |

The full **AT time × temperature surface** (6 times × 5 temperatures, µg from 150 µg HDT charged,
**[D]** read off Figure 2 at 300 dpi, ±0.3 µg) is transcribed in `z3_index.md` §1G. It
self-validates: 16.1 µg at 10 min / 100 °C = **10.9 mol %** against the paper's own stated
"nearly 10 % yield" ✔, and Table 3's 10.9 µg at 20 min / 100 °C agrees with the figure's ≈10.1 µg
to **8 %** — a consistency check the paper never makes, and it passes.

⚠️ **DOUBLE-COUNTING HAZARD: Hofmann & Schieberle 1996b Tables 1 and 2 are VERBATIM
RE-PUBLICATIONS of Hofmann & Schieberle 1998 Tables 2, 8 and 10.** The MFT rows
(ribose 55.3 / 19.8 / 2.5; glucose 0.3 / 1.9 / 0.4; rhamnose 0.1 / 0.8 / 0.1 µg) and the C₂+C₃
rows (268.1 µg aqueous, 1553.9 µg dry) are **not second measurements**. **DO NOT INGEST THEM AS
NEW.** See §F #16.

## A.4 — THIOL CONSUMPTION: FOUR TEMPERATURES, FOUR DIFFERENT SINKS

**★ This is the section a build wave is most likely to get wrong. Read §B.7 before using it.**

| T | matrix / conditions | channel | parameter | units | source | prov. | verdict |
|---:|---|---|---|---|---|---|---|
| **25 °C** | reconstituted coffee, **1 % total solids, 0.01 M acetate pH 5.2, aerobic**, thiols 3.2–9.7 µmol/L, SPME vs coffee-free blank, duplicates, **no CIs** | covalent addition to matrix electrophiles | **k(FFT) > 7.70 × 10⁻⁴** | **s⁻¹** | Charles-Bernard 2005 T2 p. 4430 | [F] | **USE — lower bound only.** ⚠️ **UNIT ERRATUM: the table header prints "[mol⁻¹ s⁻¹]"; the values are s⁻¹**, proven against five of the paper's own half-lives |
| 25 °C | " | " | k(PentSH) 1.83 × 10⁻⁴ · k(BuSH) 1.42 × 10⁻⁴ · k(PropSH) 1.32 × 10⁻⁴ · k(EtSH) 1.19 × 10⁻⁴ · k(2BT) 8.11 × 10⁻⁵ · k(MMBF) 2.12 × 10⁻⁵ · **k(2M2P) 1.02 × 10⁻⁶** | s⁻¹ | " | Charles-Bernard 2005 T2 | [F] | **USE — a ten-constant ladder spanning 755×** |
| 25 °C | " | — | t½: FFT/BnSH/PhSH **< 15.0 min**; PentSH 63; BuSH 81; PropSH 87; EtSH 97 min; 2BT 2.37 h; MMBF 9.1 h; 2M2P **189 h** | — | " | [D] | USE |
| **30 °C** | **pH 6.0, 0.1 M phosphate, 12.5 g/L coffee melanoidin (MW > 3000), FFT 438 µM** | FFT + melanoidin → **thioether** | **k = 5.6 × 10⁻² min⁻¹ = 9.4 × 10⁻⁴ s⁻¹**, t½ **12.4 min** | min⁻¹ / s⁻¹ | Hofmann & Schieberle **2002** (file `hofmann2001.pdf`), Fig. 6 | **[D]** — the paper publishes **no** rate constant | **★ USE — first order CONFIRMED against the data (k constant ±18 % over 80× extent)** |
| 30 °C | pH 6.0, lys/glycolaldehyde CROSSPY model, 30 min | FFT + pyrazinium dication | 9.8 × 10⁻⁴ | s⁻¹ | Hofmann 2002 T2 | [D] | USE — single point |
| 30 °C | pH 6.0, albumin/glycolaldehyde | — | 6.5 × 10⁻⁴ | s⁻¹ | Hofmann 2002 T2 | [D] | USE — single point |
| 30 °C | pH 6.0, albumin/glucose | — | 3.0 × 10⁻⁴ | s⁻¹ | Hofmann 2002 T2 | [D] | USE — single point |
| 30 °C | pH 6.0, (heated) chlorogenic acid | — | 8.4 × 10⁻⁵ / 4.6 × 10⁻⁵ | s⁻¹ | Hofmann 2002 T2 | [D] | **USE — essentially inert** |
| **80 °C** | real coffee brew, thermos, 50 g powder/L | FFT loss | 0.023 min⁻¹ = 3.8 × 10⁻⁴ s⁻¹, t½ ~30 min | min⁻¹ | Hofmann 2002 Fig. 1 | [D] | **USE-Q — ⚠️ SLOWER than the 30 °C model systems, because the brew's electrophile pool was partly consumed during extraction. Do NOT pair with a 30 °C value to get an Ea** |
| 80 °C | " | MMBF loss | 0.0084 min⁻¹ = 1.4 × 10⁻⁴ s⁻¹, t½ ~83 min | min⁻¹ | Hofmann 2002 Fig. 1 | [D] | USE — clean first order (0.0071–0.0093 over 7× extent) |
| **50 °C** | ribose/cysteine process flavour, **pH 5.0, 0.5 M phosphate**; air ≈ argon | **acid-catalysed electrophilic oligomerisation at C-5** — NOT oxidative | **MFT 59 % of initial per day, ZERO ORDER**; FFT 28 %/day | %/day | van Seeventer 2001 T1 | [M] | **★ USE-Q — ZERO ORDER; DO NOT convert to a first-order half-life** (Z2). Mechanism evidence: **85 % H–D exchange at C-5 vs 10 % at C-4**; MFT mass balance **fails to close** as thiol + disulfide while MB, FFT and DMFT all close |
| **115 °C** | thiamine/xylose/methionine, **pH 4.9 ± 0.1 phosphate**, 60 min | **oxidative dimerisation** + **radical coupling to MeSH** | **21–55 % of free MFT already consumed into measured products** (GSH arm); 3–19 % (Cys arm) | fraction | Zhang 2024 Figs. 2d/e/f | [D] | **★ USE-Q — LOWER BOUNDS; a third sink (melanoidin/quinone) is named and never measured** |
| 115 °C | " | MMFT (zero order) | **0.0028 (Cys), 0.0031 (GSH)** | **UNITS NOT STATED** | Zhang 2024 §3.1 p. 4, from absent Table S6 | [F] | **REFUSE AS NUMBERS — units unrecoverable (K3 addendum §K3.3); run duration inconsistent with §2.3. INGEST THE RATIO 1.107 and the two-stage structure only** |
| 115 °C | " | MMFT slope, t < 90 min | Cys **≈3× control**; GSH **≈2.5× control** | ratio | Zhang 2024 §3.1 | [F] | **RATIO-ONLY, USE** |
| 115 °C | " | MMFT slope, t > 90 min | Cys **→ control**; GSH **≈4× control** | ratio | Zhang 2024 §3.1 | [F] | **★ RATIO-ONLY, USE — the slope BREAKS at 90 min. A one-stage rate law cannot express this** |
| **120 °C** | Amadori + cysteine, unbuffered, pH 6–8 | dimerisation | **6.5–9.6 % of MFT-equivalents**, pH-invariant | fraction | Zhou 2023 T1 | [D] | **★ USE — independently corroborates Zhang's dimer channel: the free-thiol arms of the two papers agree within 1.3×** |
| **30 °C** | pH 6.0 | 2 RSH → RSSR | **< 1.5 % of the thiol flux** (< 6 µg disulfide from 400 µg FFT consumed) | fraction | Hofmann 2002 | [M] | **★ USE — KILLS the disulfide sink at 30 °C, and fixes the stoichiometry as 1:1 thioether addition, FIRST order in thiol** |
| **6 °C** | **diethyl ether / DCM / n-pentane**, 1–10 days, purified thiols, no sugar/amino acid/buffer/H₂S | O₂-driven oxidative dimerisation | MFT **≈13.5×** more oxidation-labile than FFT (10-day loss **53.3 % vs 5.5 %**); mass balance closes to 0.7 % | ratio | Hofmann 1996 | [M] | **RATIO-ONLY, USE-Q — ⚠️ the paper's OWN framing is "to check this source of ARTIFACT formation". ⚠️ Z2 transplant warning: the ORDERING holds in aqueous but the MAGNITUDE is 2.1×, not 13.5× (van Seeventer, 50 °C: MFT 59 vs FFT 28 %/day). Ingest the split; do NOT carry 13.5× into aqueous systems** |
| — | — | matrix electrophile site density | **8–10 mmol/g dry coffee solids** (= 0.08–0.10 M at 1 % t.s.) | mmol g⁻¹ | Charles-Bernard 2005 p. 4428, hydroxylamine titration | [M] | **★ USE — the only published site density. Converts the pseudo-first-order sink into a bimolecular one that scales with matrix loading** |
| — | — | covalent capacity | FFT bound **≥ 0.028 mmol per g melanoidin** | mmol g⁻¹ | Hofmann 2002 | [D] | **USE — LOWER BOUND (FFT was limiting)** |
| — | — | bimolecular recast | k₂(EtSH) ≈ **1.3 × 10⁻³**; k₂(FFT) ≳ **8.6 × 10⁻³** | **M⁻¹ s⁻¹** | Charles-Bernard T2 × site density | [D] | **USE-Q — order-of-magnitude only** |
| **Ea(thiol consumption)** | — | — | **DOES NOT EXIST** | — | — | — | **★ REFUSE ANY VALUE. See §C.1 for K1's verbatim declaration and §B.7 for why one term is the wrong model** |

**Rate-attenuation table** (Charles-Bernard 2005 T3, dimensionless `1/k_rel`, same conditions;
additive in mmol per g dry coffee solids) — **[M], STRUCTURAL, directly usable as branch-fraction
constraints**:

| additive (mmol/g) | EtSH | PropSH | BuSH | PentSH | **FFT** |
|---|---:|---:|---:|---:|---:|
| aerobic (ref) | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| **anaerobic** | 17.9 | 19.8 | 15.9 | 16.4 | **45.1** |
| hydroxylamine (10) | 10.5 | 11.4 | 10.8 | 13.4 | **4.0** |
| anaerobic + hydroxylamine (21) vs aerobic | 64.8 | 82.5 | 101.3 | 236.1 | 92.6 |
| Na₂SO₃ (1) | 51.5 | 61.7 | 18.4 | 62.7 | 67.9 |
| ascorbic acid (1) | 9.0 | 7.8 | 4.4 | 3.7 | **23.6** |
| **caffeic acid (0.28)** | **0.3** | **0.3** | **0.3** | **0.3** | 1.0 |
| DTPA (0.5) | 0.5 | 0.7 | 0.7 | 1.0 | 1.0 |

⇒ aliphatic thiols have **one** composite O₂-dependent electrophilic-addition channel;
**FFT retains ~25 % of its rate with all electrophiles blocked**, so benzylic/furfuryl thiols
need a **second, non-electrophilic, O₂/radical channel**. Caffeic acid at 0.3 is a **pro-oxidant**
(3× acceleration); DTPA ≈ 1 means **Fenton chemistry is not the main route**.

**★ THE CROSS-VALIDATION** — two labs, two methods, two matrices, agreement within 25 %:
Hofmann Fig. 6 **9.4 × 10⁻⁴** (30 °C, SIDA) · Hofmann T2 **9.8 × 10⁻⁴** (30 °C, SIDA) ·
Charles-Bernard **> 7.7 × 10⁻⁴** (25 °C, headspace ratio).
⇒ **adopt k(FFT sink) ≈ 9 × 10⁻⁴ s⁻¹ at 30 °C, ~10 g/L coffee solids, pH 5–6.**

## A.5 — LIPID OXIDATION

| parameter | value | units | conditions | source | prov. | verdict |
|---|---|---|---|---|---|---|
| hexanal share of the aldehyde distribution | **11–26 %**, across **26 columns at ONE temperature (180 °C)** | fraction | methyl linoleate hydroperoxides | Frankel 1989 | [M] | **★ USE-Q — this REFUTES the shipped `hexanal 0.37` in `lipid_oxidation_calibration.json` (used at `src/lipid_oxidation.py:446, :512`): 0.37 sits ABOVE the paper's entire measured range** |
| hexanal share vs hydroperoxide isomer | **11 / 13 / 20 %** with zero additive | fraction | " | Frankel 1989 | [M] | **STRUCTURAL — the branch fraction is not a constant** |
| hexanal share vs hydrogen-donor loading | **2.0×** | ratio | " | Frankel 1989 | [M] | **STRUCTURAL — and the two effects INTERACT (absent for pure *cis,trans*-13)** |
| α-tocopherol effect | **total volatiles ÷ 2.5–6.7×** *while* **shifting the distribution TOWARD hexanal** | ratio | " | Frankel 1989 | [M] | **★ STRUCTURAL — a fixed split gets both effects wrong, in OPPOSITE directions for total vs hexanal-specific output. Real protein matrices carry tocopherols** (Z2 #19) |
| **nonanal from linoleate** | **ABSENT — appears in no table, figure or sentence** | — | " | Frankel 1989 | [M] | **★ USE AS A NEGATIVE — SUPPORTS Wave P's oleate-nonanal routing. `hexanal_path: nonanal 0.15` is unsupported if fed linoleate hydroperoxide** |
| 9-hydroperoxide's aldehyde partner | **2-nonenal**; in the methyl ester, **methyl 9-oxononanoate** (2nd most abundant, **up to 47 %**) | fraction | " | Frankel 1989 | [M] | **STRUCTURAL — the repo is likely missing 2-nonenal and 9-oxononanoate entirely** |
| any absolute yield; any Ea | **NONE EXIST** | — | **one temperature only (180 °C)** | Frankel 1989 | — | **REFUSE** |

## A.6 — PROTEIN BINDING

`per-gram` in the basis column means **no protein molar mass is needed** — the strongest
provenance available.

| ligand | protein | constant | **basis** | conditions | source | prov. | verdict |
|---|---|---:|---|---|---|---|---|
| 2-heptanone | β-lactoglobulin | K_b **330 M⁻¹** ⇒ **8.97 × 10⁻³ L/g** | **36 800 g/mol — RECOVERED BY ARITHMETIC (dimer), NOT stated by source** | pH 3, 50 mM NaCl, 30 °C, 50 µL/L aroma, 0–4 % protein | Andriot 2000 | [F]+[D] | **USE-Q — label `recovered_by_arithmetic`, NOT `stated_by_source`** |
| 2-octanone | β-lactoglobulin | K_b **950 M⁻¹** ⇒ **2.58 × 10⁻² L/g** | " | " | Andriot 2000 | [F]+[D] | USE-Q |
| **2-nonanone** | β-lactoglobulin | K_b **2 440 M⁻¹** ⇒ **6.63 × 10⁻² L/g** | " | " | Andriot 2000 | [F]+[D] | USE-Q |
| hexanal | lupin protein isolate (91 % protein) | **5.15 × 10⁻² L/g** (34 % at 10 g/L) | **per-gram** | pH 7.0, 30 °C, 3 h, 125 rpm, 5 mg/L, static HS-GC-MS | Barallat-Pérez 2024 Fig. 4 | [D] | **USE-Q — headspace method; see the split rule below** |
| nonanal | lupin protein isolate | **3.17 × 10⁻¹ L/g** (76 %) | **per-gram** | " | Barallat-Pérez 2024 | [D] | USE-Q |
| 2-nonanone | lupin protein isolate | **4.71 × 10⁻² L/g** (32 %) | **per-gram** | " | Barallat-Pérez 2024 | [D] | **USE** |
| hexanal | pig gastric mucin | **6.4 × 10⁻¹ L/g** (6 % at 0.1 g/L) | **per-gram** | " | Barallat-Pérez 2024 | [D] | USE-Q |
| nonanal | pig gastric mucin | **2.35 L/g** (19 %) | **per-gram** | " | Barallat-Pérez 2024 | [D] | USE-Q |
| 2-nonanone | pig gastric mucin | **1.36 L/g** (12 %) | **per-gram** | " | Barallat-Pérez 2024 | [D] | USE |
| 2-heptanone / 2-octanone / 2-nonanone / nonanal | soy | 4.40 × 10⁻³ / 1.24 × 10⁻² / **3.72 × 10⁻²** / 4.38 × 10⁻² L/g | **100 000 g/mol — STATED BY SOURCE, three times** | pH 8.0, 30 mM Tris, **10 mM 2-mercaptoethanol**, 25 °C, equilibrium dialysis | Damodaran 1981 | [M]+[D] | **★ USE — the stated molar basis is the datum S4 named as blocking** |
| hexanal / 1-hexanol | soy, partially denatured | 1.47 × 10⁻³ / 6.99 × 10⁻⁴ L/g; saturation **0.847 and 0.889 mg/g protein**; K **173.4 / 80.3 M⁻¹** | **per-gram** | gel filtration | Arai 1970 **via** Damodaran | [C] | **USE-Q — second-hand. The 1-hexanol row is the FIRST anchor of any kind on that lane; it implies f_free ≈ 0.93 against the shipped fitted 0.063 — 15× off** |
| **MFT (MFO)** | **kafirin** | **K = 1092.66 (298 K) / 848.23 (318 K) / 404.24 (338 K) M⁻¹**; ΔH **−20.58 kJ/mol**, ΔS **−10.17 J/mol/K** | 24 652 Da, stated | **53 % aqueous ethanol**, kafirin 200 mg/L (8.1 µM) | Zhu 2023 (`zhu2023_kafirin_binding_extraction.md`) | [F] | **USE-Q — the ONLY MFT binding constant in the corpus; van't Hoff checks exactly. ⚠️ 53 % ETHANOL — NOT transferable to aqueous or low-moisture systems** |
| BMFD / DT / MMFD | kafirin | K(298 K) 1934.33 / 327.42 / 345.36 M⁻¹ | " | " | Zhu 2023 | [F] | **USE-Q** — ⚠️ **DT's printed ΔH (14.22) disagrees with a van't Hoff fit to its own K's (+18.3)** |
| soy K, native → **denatured** | soy | **930 → 1240 M⁻¹ (+33 %)**, sites unchanged | 100 kDa | after **90 °C / 1 h** | Damodaran 1981 | [M] | **★ USE — fills a declared gap where the repo's file records "qualitative only, NOT MODELLED"** (Z0 #18) |
| soy K vs temperature | **flat 25 → 45 °C** | — | " | — | Damodaran 1981 | [M] | **USE — a modelling licence** |
| 3-sulfanylhexan-1-ol | whole human saliva | **STRANDED** | — | pH 7.6, 0.59 g/L protein, 22 °C, 45 min | Starkenmann 2008 | — | **REFUSE — no basis, no site stoichiometry, mechanism unresolved. Apparent per-gram avidity 10³–10⁴ L/g is 3–5 orders above every hydrophobic constant. A 1:1 stoichiometry would need an effective MW < ≈8 000 g/mol, smaller than any major salivary protein. Authors: "Whether this absorption results from physicochemical interactions, chemical transformations, or covalent linkage to glycoproteins remains unclear."** |
| chain-length slope | **2.72×/CH₂** (Andriot β-lg) · **2.9×/CH₂** (Damodaran soy) | ratio | two proteins, two labs | — | K2 §B.6 | [D] | **★ USE — they agree to 6 %, against the shipped `a_p·Pow` implied ≈3.6×/CH₂** |
| **binding ceiling** | **1.25× / 2.4× / 3.7×** at 4 % β-lg (C7/C8/C9); **1.9× / 3.6× / 7.6×** extrapolated to 100 g/L | ratio | Andriot 2000 Fig. 1, eq 3 | — | K2 §B.4 | [M]+[D] | **★ STRUCTURAL — the number that closes the "can binding explain matrix thresholds?" question. NO.** |
| gas–liquid partition, K_ga | 7.2 × 10⁻³ / 1.03 × 10⁻² / 2.4 × 10⁻² (C7/C8/C9) | — | pH 3, 30 °C | Andriot 2000 | [D] **[fig]** | USE-Q, ±5–8 % |
| mass-transfer h_D | ≈ 2.0 / 2.5 / 5.0 × 10⁻⁵ | m/s | " | Andriot 2000 | [D] | **REFUSE FOR TRANSFER — apparatus-specific; the authors attribute the fit deviation to unstirred layers** |

## A.7 — ODOUR THRESHOLDS

### A.7.1 In water (µg/L) — the largest single list, NEW THIS WAVE (Zhou 2023 SI Table S2)

⚠️ **PROVENANCE: the SI gives NO citation for any threshold** — no reference column, no
footnote. The **basis is declared** ("Threshold in water") but the numbers are un-anchored.
**Ingest as `basis_declared: true, provenance: unknown`. Do NOT record as measured.**

| compound | threshold, µg/L in water | | compound | threshold, µg/L in water |
|---|---:|---|---|---:|
| **bis(2-methyl-3-furyl) disulfide** | **0.00032** | | thiazole | 38 |
| **2-methyl-3-furanthiol (MFT)** | **0.005** | | methylpyrazine | 60 |
| **2-furfurylthiol (FFT)** | **0.006** | | 2,5-dimethylpyrazine | 80 |
| **2-acetylthiazole** | **10** | | thiophene | 84 |
| 2,4-dimethylthiazole | 18 | | 2,6-dimethylpyrazine | 400 |
| — | — | | furfural | 3000 |
| 2-thiophenecarboxaldehyde | 5000 | | furan | 4500 |
| 2-acetylfuran | 10000 | | pyrazine | 75000 |

⚠️ **Absent from Table S2** (7 of the 22 quantified compounds): 2,3-butanedione,
2-methylthiophene, 5-methyl-2-thiophenecarboxaldehyde, thieno[3,2-b]thiophene, 2-methylthiazole,
di-2-furfuryl sulfide, bis(2-furfuryl) disulfide.

### A.7.2 Other aqueous anchors

| compound | matrix | threshold | method | source | prov. | verdict |
|---|---|---:|---|---|---|---|
| **(R/S)-3-sulfanylhexan-1-ol** | mineral water | **22 ng/L** | 3-AFC, ascending 4–500 ng/L, 30 trained panellists, geometric mean of 2 sessions, **RETRONASAL (30 mL held in mouth 5 s)** | Starkenmann 2008 p. 9577 | [M] | **★ USE-Q — ALREADY AN IN-MOUTH NUMBER. A saliva correction on top would DOUBLE-COUNT** |
| same | water | 1–60 ng/L | not stated | Vermeulen 2005 via Starkenmann | [C] | USE-Q |
| hexanal | water | **5.0 µg/kg** | not stated | **recovered by K2 from Xin 2026's own OAV arithmetic** | [D] | USE-Q |
| nonanal | water | **≈1.1 µg/kg** | not stated | " | [D] | USE-Q |
| 2-pentylfuran | water | **5.800 µg/kg** | not stated | " (exact to 4 s.f. at both ends of the printed range) | [D] | **USE — ★ now recovered TWICE, from FOUR independent numbers across TWO papers (K2 from Xin 2026; K3 from Xin 2026b). Firmly established as this group's value** |
| **1-octen-3-ol** | water | **1.5 µg/kg** | not stated | **printed VERBATIM in Xin 2026b p. 9 — no citation given** | **[C]** | **★ USE — the only threshold in either Xin paper that is actually printed rather than inferred** |
| 1-octen-3-ol | water | **1.500 µg/kg** | not stated | back-solved from Xin 2026b's own OAV range | [D] | **★ USE — it MATCHES the printed 1.5, which independently VALIDATES K2's OAV-inversion method** |
| **2-heptanone** | water | **140.0 µg/kg** | not stated | back-solved from Xin 2026b's own OAV range | [D] | USE-Q |
| 14 further thresholds | — | — | — | Hofmann 1995 | [M]/[C] | USE-Q — ordinal-grade (Z0) |

⚠️ **Both Xin papers declare the OAV denominator an "odor threshold in AIR (µg/kg)" while the
values recovered are unambiguously the classical WATER thresholds.** The identical defective
sentence appears in **two different journals**. It is this group's standing practice, not a
one-off slip — which strengthens K2's recommendation to fix the OAV denominator (§C.2, §F #9).
| **2-acetyl-1-pyrroline** | — | **NO THRESHOLD EXISTS IN EITHER 2-AP PAPER** | — | Schieberle 1989 + 2005 | — | **GAP — 2-AP's whole claim to being a top odorant is an OAV argument the repo cannot score. [RETRIEVE] Buttery, Ling & Juliano 1982** |

### A.7.3 ★ The threshold RATIOS — the form that survives unknown provenance

| ratio | value | consequence |
|---|---:|---|
| **MFT : FFT** | **0.833** | MFT is 1.2× more potent than FFT, in water |
| **bis(2-methyl-3-furyl) disulfide : MFT** | **0.064** | **★ the dimer is 15.6× MORE potent than its own monomer** |
| MFT : 2-acetylthiazole | 1 : 2000 | |
| methylpyrazine : pyrazine | 1 : 1250 | methyl substitution buys 3 orders in the pyrazine family |

**★ THE CONSEQUENCE THAT CHANGES A MODULE.** Only **6.5–9.6 %** of MFT-equivalents sit in the
disulfide dimer (§A.3.3), but the dimer is **15.6× more potent**, so the dimer's OAV **matches or
exceeds** the monomer's — Zhou SI Table S2 prints, at pH 7, **MFT OAV 3.18 × 10⁵ vs disulfide
3.21 × 10⁵**. ⇒ **Mass lost to MFT dimerisation is NOT aroma lost.** Any objective that scores
the Zhang 2024 / Zhou 2023 dimerisation channel as a pure loss is wrong by roughly the threshold
ratio. **This qualifies the entire thiol-consumption module, §A.4.**

### A.7.4 In matrices — and the reclassification that must travel with them

| compound | **cooked beef** (ppm) | **3 % gelatin @22 °C** (ppb) | paraffin oil (ppb) | milk (`?/kg`) | aqueous (ppb) |
|---|---:|---:|---:|---:|---:|
| pentanal | 2.67 | 41 | — | — | 12 |
| **hexanal** | **5.87** | **58** | 120 | — | 4.5 |
| heptanal | 0.23 | 79 | 250 | — | 3 |
| t-2-hexenal | 7.87 | 109 | — | — | 3 |
| t-2-octenal | 4.20 | 109 | — | — | 3 |
| t,t-2,4-decadienal | 0.47 | 64 | 135 | — | 0.07 |
| octanal | — | — | — | 160 | — |
| nonanal | — | — | — | 1 600 | ≈1.1 µg/kg |
| 2-nonanone | — | — | — | 52 000 | — |
| propanoic / butanoic / octanoic acid | — | — | — | 51 200 / 7 500 / 25 600 | — |
| ethyl hexanoate | — | — | — | 1 024 | — |

Gelatin at other temperatures (ppb, 4 / 22 / 37 / 60 °C): pentanal 47/41/34/22 · hexanal
90/58/34/38 · heptanal 108/79/62/50 · t-2-hexenal 170/109/79/60 · t-2-octenal 140/109/105/81 ·
decadienal 112/64/89/64. **Temperature effect 1.7–2.8× over 4→60 °C, non-monotone in 2 of 6.**

⚠️ **MANDATORY QUALIFICATIONS, all from K2 §D.2:**
- **Brewer 1995's beef values are NOMINAL DOSES ADDED TO RAW BEEF BEFORE A 70 °C COOK**, never
  analytically verified at the moment of sniffing. **RECLASSIFY as `dose_added_pre_cook`, not
  `threshold_in_matrix`.**
- **Tian 2020's unit cell prints a literal `?`**, verified at 900 dpi. Contextually µg/kg;
  **not asserted**. All seven milk rows carry a **factor-of-1000 basis risk**.
- **No same-method matrix-vs-water pair exists anywhere.** Criteria differ: Brewer 50 %
  forced-choice; **Vega 75 % detection UNCORRECTED for chance**; Tian 50 % 3-AFC;
  Guadagni not stated.
- Wick et al. 1967's methional 6.1 / phenylacetaldehyde 0.94 / **nonanal 7.6 ppm** are cited by
  **both** Brewer and Vega, **measured by neither**, and labelled "in beef" by one and "in meat
  slurries" by the other. **Attributing them to Brewer or Vega would be laundering.**

### A.7.5 Matrix / water ratios — and the verdict on a general factor

| compound | beef/water | gelatin/water | oil/water | beef/gelatin |
|---|---:|---:|---:|---:|
| pentanal | 223× | 3.4× | — | 65× |
| hexanal | **1 304×** | **12.9×** | 26.7× | **101×** |
| heptanal | 77× | 26.3× | 83× | **2.9×** |
| t-2-hexenal | 2 623× | 36.3× | — | 72× |
| t-2-octenal | 1 400× | 36.3× | — | 38.5× |
| t,t-2,4-decadienal | **6 714×** | **914×** | 1 929× | 7.3× |
| **geometric mean** | **906×** | **33×** | 244× | **27×** |
| **max / min** | 87× | **269×** | 72× | 35× |
| **1-σ band (log10 SD)** | **0.71 dec ⇒ 27× wide** | **0.80 dec ⇒ 41× wide** | — | 0.63 dec ⇒ 18× wide |

**★ VERDICT: NO GENERAL MATRIX CORRECTION FACTOR.** Applying the gelatin/water geometric mean
(33×) uniformly would make pentanal **9.8× too high** and decadienal **27.7× too low** — a 10×
and 28× error in opposite directions. **Ship a LOOKUP TABLE with an explicit
`no_factor_available` state**, carrying the fields `criterion`, `thermal_step_after_dosing`,
`concentration_verified`, `aqueous_reference_source`, `cross_study_cross_method: true`
(which is **every** ratio above).

**The ONE parametric term that is defensible: the α,β-unsaturation penalty ≈ 2–3×** —
gelatin 2.8×, beef 2.0×, two matrices, two panels, and a mechanism (2-alkenals are Michael
acceptors; t-2-hexenal has a **lower** logP than hexanal and a 2.8× **larger** penalty).
**The chain-length term does NOT transfer** — monotone in gelatin (3.4 → 12.9 → 26.3, i.e.
2.78×/CH₂) and collapsed in beef (223 → 1 304 → 77; heptanal breaks it by 17×).

### A.7.6 Saliva

| quantity | value | modality | verdict |
|---|---:|---|---|
| analytical free-thiol LOD, **water** | 0.001 mg/L | HPLC-fluorescence, Acrylodan | — |
| analytical free-thiol LOD, **crude saliva** | **60 ± 10 mg/L** | same method, same instrument | — |
| **⇒ same-method saliva/water quench** | **6 × 10⁴** | **analytical** | **USE — [D]** |
| the paper's printed "3 × 10⁶" | — | **analytical LOD ÷ SENSORY threshold** | **REFUSE — cross-modality** |
| sensory saliva quench | **NO NUMBER EXISTS** | qualitative only | **GAP** |
| salivary protein | **0.59 ± 0.2 mg/mL, pH 7.6 ± 0.1**, pooled from 4 donors | [M] | USE |
| the paper's quoted literature range "0.6–1.2 mg/mL **or 230–250 mg/mL**" | — | — | **REFUSE the second half — a transcription error in the source** |

## A.8 — EXTRUSION / MATRIX PATH

| study | base | dosed precursor | process | key quantified outputs | verdict |
|---|---|---|---|---|---|
| **Guo 2020** *Food Hydrocolloids* 105, 105752 | SPI 94.2 % + wheat gluten at 0/10/20/30/40 % of SPI | **NOT a Maillard precursor** — encapsulated flavour powder (β-cyclodextrin, 11 % fat) at **1 %** | twin-screw FT-36, 8 zones **20/50/80/150/140/100/80/60 °C**, L/D 26:1, compression 4.6:1, die 10.0 mm, moisture **50/60/70/80 %**, feed **6 kg/h *or* 30 g/min (contradictory, 3.3×)** | **retention %** of 16 terpenes: total **52.46 → 44.07 → 32.24 → 23.04 %** over 50→80 % moisture (**2.28×**, all steps significant, 14/14 compounds agree); wheat gluten **35.46 → 42.04 → 44.07 → 43.78 → 33.39 %** (optimum 20–30 %) | **USE — RETENTION only** |
| **Conti 2025 I + II** *Food Res. Int.* 208, 116169 / 218, 116938 | soy protein **concentrate** (Arcon SM, ≥70 g/100 g protein dry) | **THIAMINE·HCl > 99 %, 1.5 % w/w**, added 2 h before extrusion. **ONE dose level** | **single-screw** RXPQ Labor 24, 5 zones, z1–3 ~40/60/80 °C, z4 = z5 − 15, **z5 = 180/160/140 °C**, moisture **30/34/38 % DRY basis (= 23–28 % wet — LOW-moisture TVP, NOT HME)**, L/D 15.5:1, die 3.6 mm, feed 170 g/min, **216 rpm**. ⚠️ **moisture confounded with temperature by design (diagonal, not factorial)** | 98 and 71 volatiles, **µg/g SEMI-QUANT vs hexanal-D12**: **pyrazines 49.74 vs 3.38 µg/g (14.7×)** · all sulfur **44.81 vs 19.78 (2.27×)** · lipid-oxidation aldehydes **3.80 vs 7.65 (0.50× — MILD wins)** · **hexanal 1.95 vs 4.70 µg/g** (the ONLY absolute values; the IS is hexanal-D12) · largest single volatile **4-methyl-5-(2-hydroxyethyl)thiazole 20.25 µg/g** | **USE-Q — NET (formation − retention)** |
| **Xin 2026** *Food Hydrocolloids* 182, 113124 | NPI : SPI : WG = **20 : 70 : 10** | **EIGHT carbohydrates at 6 % w/w dry, ADDITIVE, + a true no-carbohydrate control.** **ONE dose level** | **twin-screw**, screw 704 mm × 22 mm, L/D 32:1, **cooling die 70 °C**, **240 rpm**, dry feed 2 kg/h, water 3.7 kg/h, **final moisture ≈65 % — TRUE HME**, zones **60/80/110/135/140/140/140 °C**, **3 independent runs per formulation** | **total pyrazines µg/kg: FR 6 621.64 > SU 2 200.79 > RI 1 389.14 > GL 1 149.61 > {BG, WS, MA, control} ≈700–780 > XY ≈600** · 2-ethyl-3,6-dimethylpyrazine in FR **5 638.91** · 2-pentylfuran **BG 100 259.10 / WS 97 794.41 / control ≈97 000 / RI ≈52 000** · aldehyde totals GL 3 506.49 / FR 3 483.59 / control 2 433.19 | **★ USE — the only dosed-carbohydrate factorial in a real HME** |

**Reconciliation (K2 §C.1):** Guo measures **retention**; Conti measures the **net**; they do not
contradict. Conti's own text: *"the greater number and quantity of volatile compounds in the
FTSPs obtained at 30 % M/180 °C may be attributed to thiamine degradation, which results in the
generation of more compounds."* ⇒ **FORMATION BEATS RETENTION under precursor dosing. The
extrusion lane needs BOTH terms and must not fit one to data generated by the other.**

## A.9 — 2-ACETYL-1-PYRROLINE (the repo's #1 named completeness gap)

| parameter | value | conditions | source | prov. | verdict |
|---|---|---|---|---|---|
| **1-pyrroline (2 mM) + pyruvaldehyde (0.1 mM)** | **1140 µg Acp, 72 % of the volatile fraction** | 100 mL, **pH 7.0, 2 h**, SIDA | Schieberle 1989 | [M] | **★ USE — the best 2-AP anchor available; replaces an `[ABS]` citation with a full-text primary** |
| synthesised 2-(1-hydroxy-2-oxo-1-propyl)-pyrrolidine → ATHP | **23.4 % at pH 7** vs **0.9 %** from the free precursors under identical pH/T/t | — | Schieberle 2005 | [M] | **★ STRUCTURAL — the ADDITION step is rate-limiting, not the ring enlargement** |
| 1-pyrroline + pyruvaldehyde vs proline + pyruvaldehyde | **240×** | — | Schieberle 1989/2005 | [M] | **★ STRUCTURAL — the STRECKER step is the 2-AP bottleneck, not the acylation** |
| proline + glucose / fructose / sucrose | **< 0.1 µg Acp** each | dilute aqueous | Schieberle 2005 | [M] | **★ USE AS A NEGATIVE** |
| proline + **DHAP** | **13.6 µg** | " | Schieberle 2005 | [M] | **★ ⇒ >136× separation. The carbonyl must come from a FRAGMENTATION/glycolytic product, not the intact sugar** |
| 2-AP aqueous vs dry | **58.5 vs 1.1 µg (53× higher AQUEOUS)** | — | Schieberle 2005 | [M] | **USE-Q — opposite to the usual low-moisture heuristic. ⚠️ confounded by T and t** |
| donor selectivity | **hydroxy-2-propanone → ONLY ATHP; 2-oxopropanal → ONLY AP** | — | Schieberle 2005 | [M] | **★ STRUCTURAL — MUTUALLY EXCLUSIVE. Corrects `isotope_topology_evidence.md` §25 BEFORE implementation** |
| Schieberle 2005 pH ladder | 5 rows incl. **pH 9.0: 38.4 µg, 3.1 mol %** | — | Schieberle 2005 T2 | [M] | **USE — ⚠️ the PDF text layer SILENTLY DROPPED this entire row; read from the 300 dpi render** |
| odour threshold for 2-AP | **DOES NOT EXIST IN EITHER PAPER** | — | — | — | **GAP** |

## A.10 — ISOTOPE / TOPOLOGY / STOICHIOMETRY

| constraint | value | conditions | source | prov. | verdict |
|---|---|---|---|---|---|
| MFT skeleton origin | **intact-C₅: 49 % unlabelled / 46 % ¹³C₅**, no fragment pattern; *"pathways via ribose fragmentation were not relevant"* | ribose + cysteine, **95 °C / 4 h, pH 5.0** | Cerny 2003 T2 | [M] | **★ STRUCTURAL — ⇒ ~93 % of MFT must come from the intact-skeleton (1,4-dideoxyosone) route** |
| **NF share of MFT** | **≤ 7 %** — and this is an **UPPER** bound (NF was spiked in at 1.5× the cysteine) | 95 °C / 4 h / pH 5.00 | Cerny 2003 T3 | [M] | **★ STRUCTURAL — but the gate must be evaluated AT CERNY'S CONDITIONS, not at 145 °C** |
| 2-mercapto-3-pentanone (the NF-diagnostic isomer) | **96 % unlabelled** (i.e. from the NF spike) | " | Cerny 2003 T3 | [M] | **USE — the sharpest single-species NF-route marker in the corpus** |
| same isomer, temperature switch | **NOT DETECTED at 95 °C**; **77.5 µg/10 mg NF at 140 °C** | Cerny 2003 vs Whitfield 1999 T1 | [M] | **★ STRUCTURAL — the NF route is SWITCHED ON by temperature** |
| butane-2,3-dione dual source | **54 : 46** (loss of ribose C-1 vs C-5) | in situ, isotope-derived, no fed-precursor artefact | Cerny 2004 | [M] | **★ USE — a scoreable step-level branching ratio** |
| thiazole / methylthio splits | 65 : 35 / 87 : 13 | " | Cerny 2004 | [M] | USE |
| C₂+C₃ route to MFT at 95 °C | *"was not relevant"* | 95 °C / 4 h / pH 5 | Cerny 2004 | [M] | **★ STRUCTURAL — a direct challenge to the Wave P C₂+C₃ lane, which Hofmann 1998 T10 shows is the STRONGEST route at 145 °C. Owner call: the lane is temperature-scoped** |
| **1,4-dideoxyosone formation** | `1-deoxyosone + α-amino acid → 1,4-dideoxyosone + RCHO + CO₂ + NH₃`, balance **verified C₇H₁₃NO₆ both sides** for glycine | — | Nedvidek 1992 Scheme 3 | [M] | **★ STRUCTURAL — REPLACES the repo's `[HH]`-token `Deoxyosone_Reduction`, which is the paper's MINORITY channel. One oxygen of the osone leaves in the Strecker aldehyde; the amino acid is a REAGENT, not a spectator** |
| norfuraneol's own dominant fate | **mercaptoketones : MFT = 16.3 : 1**; MFT is **2.6 %** of everything fed NF produces; dithiolanones alone are **33×** the MFT in the H₂S system | 140 °C, pH 4.5, fed NF + cysteine | Whitfield 1999 | [M] | **★ STRUCTURAL — four falsifiable ratio tests, immune to the response-factor caveat** |
| α-dicarbonyl + H₂S → mercaptoketone | proceeds to major products **between −15 °C and room temperature in ~1–2 h**, ethanol, ~100× H₂S | — | Mottram 1995 | [M] | **★ STRUCTURAL — convert to a barrier UPPER BOUND, labelling the Eyring conversion as the repo's inference (the paper gives no rate)** |
| glycine's own N-compound yield | **1.8× UP with lysine present; 0.35× with arginine** | 180 °C, 1 h, pH 7, ¹⁵N | Hwang 1995b | [M] | **★ STRUCTURAL — a value > 1 is UNREACHABLE by any shared-pool competition model; it forces a lysine-catalysed sugar-fragmentation channel** |
| Strecker aldehyde selectivity in a 6-AA mixture | **leucine 1.24× UP, isoleucine 0.41×, phenylalanine 0.37×, methionine 0.19×** | 2 min | Martin & Ames 2001 | [M] | **USE — ordinal** |
| total pyrazines, mixture vs sum of singles | **38 %** (K1 recompute: 36.5 %) | — | Martin & Ames 2001 | [M]+[D] | **USE-Q — ⚠️ the mix is SUGAR-LIMITED at 1.75 mol AA per mol glucose, so it is partly stoichiometric; the authors prove it by adding fructose** |
| Cys : sugar = **1:3** | chosen by an **unpublished SENSORY sweep** for "roasty, meatlike" character, **not yield**. 1:1 → H₂S-like/rubbery; 1:10 → burnt/caramel | — | Hofmann 1995 | [M] | **PROVENANCE NOTE — the repo's benchmark stoichiometry is a sensory optimum** |

---

# §B. DIRECTIONAL / ORDINAL CONSTRAINT LIST

These are **structural obligations**. Most are ratios within a single experiment and are
therefore immune to the response-factor and absolute-calibration caveats that qualify §A.

## B.1 — Sugar reactivity ordering (matrix-scoped)

| # | constraint | measured value | conditions | source |
|---:|---|---|---|---|
| B1.1 | **pentose ≫ hexose** on sulfur volatiles | ribose/glucose **MFT 10.42×, FFT 4.32×** | aqueous, 0.5 M phosphate pH 5.0, **145 °C / 20 min**, 1:3 cys:sugar | Hofmann 1998 T1 |
| B1.2 | **★ THE ORDERING INVERTS IN A REAL EXTRUDER** | **FR 6621.64 > SU 2200.79 > RI 1389.14 > GL 1149.61 > {BG,WS,MA,control} ≈700–780 > XY ≈600 µg/kg**; fructose beats ribose **4.8×**; **xylose falls BELOW the no-sugar control** | **HME, 140 °C, 65 % moisture, sheared, plant protein**, 6 % w/w each | Xin 2026 |
| B1.3 | B1.2 corroborated **independently by colour** | xylose L\* **58.49**, ribose **58.80** vs control **57.60**; ΔE **42.19 / 41.96** vs control **43.41** — **both pentoses brown LESS than adding no sugar at all** | " | Xin 2026 T1 |
| B1.4 | **sucrose (NON-reducing) is 2nd on pyrazines AND the only addition that LOWERS hexanal (−39 %)** — unremarked by the paper | — | " | Xin 2026 |
| B1.5 | fructose > glucose at a_w 0.92 but **glucose > fructose in dilute solution** | — | De Vleeschouwer I vs Claeys | K1 §5.10 |

⇒ **The pentose ≫ hexose directional row must be RE-SCOPED BY MATRIX** (K2 owner item 3),
exactly as Z0 re-scoped PH-05/PH-07. ⚠️ Xin's own explanation (ribose diverted to furans) is
**contradicted by its own data — RI has the LOWEST 2-pentylfuran and 2-butylfuran of all nine
samples.** The mechanism is open.

## B.2 — pH: FIVE independent sign-crossings. **NO family-level pH term can pass.**

| # | constraint | measured | conditions | source |
|---:|---|---|---|---|
| B2.1 | **pentose monotonic, hexose peaked** | ribose MFT rises monotonically as pH falls 7→3 (**2.5 → 19.8 → 55.3**); glucose and rhamnose go through a **MAXIMUM at pH 5** | 0.5 M phosphate, 145 °C / 20 min | Hofmann 1998 T2 |
| B2.2 | **three shapes in ONE experiment** | monotone-acid thiophenes/thiols; **peaked-at-pI trithiolanes**; monotone-alkaline pyrazines | — | Shu 1988 |
| B2.3 | **★ NEW — the C0/C1 vs C2 pyrazine inversion** | in the SAME Cys+MGO pot from the SAME α-dicarbonyl: **2,5-DMP and 2,6-DMP INCREASE with pH** (90 min: 66/106/237 and 7.9/12.6/25.9 µg/L) while **pyrazine and methylpyrazine DECREASE with pH** (47.0/20.8/15.4 and 160/97/117) | 20 mM Cys + 20 mM MGO, unbuffered, 120 °C, pH 6/7/8 | **Zhou 2023 Fig. 3a,b,g,h** |
| B2.4 | **★ NEW — MFT is NON-MONOTONE with a MAXIMUM at pH 7** while FFT is monotone-decreasing | MFT **696.99 → 1588.57 → 525.62**; FFT **813.65 → 757.97 → 325.22** µg/L | fed Amadori + Cys, **unbuffered**, 120 °C / 60 min | **Zhou 2023 T1** |
| B2.5 | ⚠️ **B2.4 CONFLICTS IN SIGN with B2.1 over their overlap** | Hofmann: MFT falls 19.8 → 2.5 from pH 5 → 7. Zhou: MFT rises 697 → 1589 from pH 6 → 7. Absolute levels **64× apart at pH 7** | see §A.3.1 vs §A.3.3 | — |
| B2.6 | MFT from fed NF collapses **≥150×** from pH 4.5 to 6.5 while **total volatiles fall only 1.7×** | 0.150 mol % → < 0.0010 mol %; class total 294 → 7 µg/10 mg NF (42×) | 140 °C, one lab, one instrument | Whitfield 1999/2001 |
| B2.7 | furfural falls **15–49×** over pH 4.5→6.5 in **aqueous** model systems, but is **unaffected pH 4→9** in roasted seed at **a_w ≈ 0.3** | — | two matrices | Meynier 1995 vs Laemont (PH-07) |
| B2.8 | **H₂S availability RISES with pH (measured) while thiols FALL with pH (measured)** | — | — | Zheng + Meynier |
| B2.9 | the **formic : acetic acid ratio INVERTS** when the amine is removed | 25 % AA / 5 % FA with glycine → **1.2 % FA / 0.65 % AA** without | — | Martins 2005 |
| B2.10 | MFT falls **> 152×**, FFT **6.1×**, furfural **15–49×** over pH 4.5 → 6.5 | constant-pH, 4-amino-acid, single-sugar factorial, ~59 series | Meynier 1995 |

**★ B2.5 IS THE MOST CONSEQUENTIAL NEW CONFLICT AND MUST NOT BE AVERAGED AWAY.** Five candidate
reconciliations are set out in `zhou2023_extraction.md` §6. The strongest is that Zhou's system is
**unbuffered** and its pH-7 run **ends at final pH 3.42** — only 0.2 units above its pH-6 run's
3.22 — so the "1 pH unit" between the two columns is really a difference in *how long the system
spent above pH 4*. **RECOMMENDATION: ingest the two as SEPARATE constraint records qualified
`buffered / pH-clamped / free-sugar-fed` and `unbuffered / pH-trajectory / Amadori-fed`. Do NOT
merge them into one pH response curve.**

## B.3 — Route topology

| # | constraint | measured | source |
|---:|---|---|---|
| B3.1 | ribose makes **54 530 µg NF against 19.8 µg MFT (~2 750×)** while glucose's ratio is only ~10 ⇒ *"MFT formation might not run exclusively via NF"* (the authors) | — | Hofmann 1998 T5 |
| B3.2 | **mercaptoketones : MFT = 16.3 : 1** from fed NF; MFT is **2.6 %** of the total | 140 °C, pH 4.5 | Whitfield 1999 |
| B3.3 | both mercaptopentanone isomers emitted, ratio **≈ 1 : 1** (74.5 : 77.5) — **the repo emits neither** | " | Whitfield 1999 |
| B3.4 | **MFT / total_volatiles ≤ 0.05** (measured 0.026) | " | Whitfield 1999 |
| B3.5 | **NF share of ribose MFT flux ≤ 7 %**, evaluated **at 95 °C / 4 h / pH 5.00** | — | Cerny 2003 T3 |
| B3.6 | the C₂+C₃ route is the **STRONGEST** MFT route at 145 °C (0.24 mol %) and *"not relevant"* at 95 °C | — | Hofmann 1998 T10 vs Cerny 2004 |
| B3.7 | **1-deoxypentosone partitions** between norfuraneol and **2,3,4-pentanetrione** (verified by two negative controls) — funnelling the whole pool into NF over-supplies everything downstream | — | Nedvidek 1992 Scheme 2 |
| B3.8 | furfural is **60× ribose** and **7× the 3-deoxyosone** as an FFT source | pH 5.0, 145 °C | Hofmann 1998 T3 |
| B3.9 | **★ NEW — the furfural → FFT carbon balance closes at ~70 % at pH 6–7 and FAILS at pH 8 (129 %, over-accounted)** | 62 % / 65 % / 48 % of furfural removed by adding Cys; FFT gained accounts for 69 / 74 / **129** % | **Zhou 2023 T1** |

## B.4 — Competition and selectivity

| # | constraint | measured | source |
|---:|---|---|---|
| B4.1 | **Cys SUPPRESSES 2,3-butanedione completely at pH 6 and 7, and by 84 % at pH 8**, while simultaneously RAISING pyrazines | 8.93 → ND · 33.63 → ND · 126.18 → 20.00 µg/L | **Zhou 2023 T1** |
| B4.2 | the cysteine acrylamide-scavenging channel equals the basic elimination channel at **[Cys] = k_E/k_E2 = 2.0 mM**; at the ~1–2 M loading used it outruns it **~1000×** | — | De Vleeschouwer II [D] |
| B4.3 | glycine's own N-compound yield **1.8× UP with lysine, 0.35× with arginine** | 180 °C, 1 h, pH 7 | Hwang 1995b |
| B4.4 | reactivity order **lysine > cysteine > isoleucine ≈ glycine** (formula count, 10 h, 100 °C) but the **browning order is DIFFERENT: Lys > Ile > Gly ≫ Cys** (~90× below lysine) | — | Hemmler 2018 |
| B4.5 | the Lys/Gly formula-count ratio is **time-dependent: ≈17× at 2 h, ≈2× at 10 h** | — | Hemmler 2018 |
| B4.6 | **sugars drive rates, amino acids drive product identity** (45–95 % formula overlap across six sugars) | — | Hemmler 2018 |
| B4.7 | **★ NEW — the sulfur additive's REDOX STATE, not its concentration, sets the MFT dimerisation branch**: cystine (oxidised) gives **54.2 %** dimer, cysteine (reduced) **8.6 %**, at near-matched molar sulfur (123.8 vs 124.8 mM S) | — | **Zhang 2024 Fig. 1 + K3 addendum §K3.5** |

## B.5 — Shape constraints (non-monotone / turning points)

| # | constraint | measured | source |
|---:|---|---|---|
| B5.1 | MFT's degradation rate has a **MINIMUM at ~50 mM cysteine** — cysteine both stabilises and destabilises it | — | van Seeventer 2001 Fig. 1 |
| B5.2 | **six of seven** additive dose responses are **non-monotone**; MFT in the GSH arm peaks at 50 mg/mL (23.6 ng/mL) and falls 1.30× by 75 | 10 levels, 115 °C / 60 min | Zhang 2024 Fig. 2 |
| B5.3 | the MMFT slope **BREAKS at 90 min**; for Cys it falls back to the control rate | — | Zhang 2024 §3.1 |
| B5.4 | glutamine's k_E Arrhenius **crosses the control's at ~180 °C** (stated explicitly) | — | Claeys 2005 |
| B5.5 | glutamine's acrylamide promotion **GROWS with T in liquid (155 %→322 %) but SHRINKS with T at a_w 0.92 (267 %→120 %)** — same lab, same amino acid, **neither paper remarks on it** | — | Claeys vs De Vleeschouwer |
| B5.6 | furan-ring ligands (BMFD, MFT) bind protein **WEAKER on heating**; polysulfides (DT, MMFD) bind **STRONGER** | ΔH −9.79 / −20.58 vs +14.22 / +19.40 kJ/mol | Zhu 2023 |
| B5.7 | **★ NEW — thiophene and 2-methylthiazole show an INDUCTION PERIOD**: ~nothing at 30 min, then a jump. **Zero-order extrapolation through the origin is wrong for these** | thiophene 0 → 0 → 11.4 → 17.6 µg/L | **Zhou 2023 Fig. 3c,f** |
| B5.8 | **★ NEW — methylpyrazine is non-monotone in pH**: the pH 7 and pH 8 series **cross** between 30 and 60 min | — | **Zhou 2023 Fig. 3h** |

## B.6 — Matrix / threshold structure

| # | constraint | measured | source |
|---:|---|---|---|
| B6.1 | **α,β-unsaturation penalty ≈ 2–3×**, reproducing in two independent matrices with two independent panels, **with a mechanism** | gelatin **2.8×**, beef **2.0×** | Vega 1994 + Brewer 1995 |
| B6.2 | the **chain-length** term does NOT transfer: monotone in gelatin (2.78×/CH₂), collapsed in beef (heptanal breaks it by **17×**) | — | " |
| B6.3 | polyunsaturated **t,t-2,4-decadienal is the extreme outlier in ALL THREE media** | 914× / 1 929× / 6 714× | Vega, Guadagni, Brewer |
| B6.4 | **perception is LESS sensitive to matrix binding than headspace is** — so a headspace-calibrated `f_free` OVER-predicts | lupin 1 wv%: headspace 31.4/72.4/36.8 % vs perceived 19.4/15.5/11.4 % (1.6×/4.7×/3.2×) | Barallat-Pérez 2024 |
| B6.5 | Andriot's ANOVA finds **A × BLG NOT significant** — the protein effect on perceived intensity is an **additive offset**, not a multiplicative factor. *"this effect is not well correlated with the retention of the aromas for the protein"* | — | Andriot 2000 p. 4250 |
| B6.6 | **protein + mucin is strongly SUPER-ADDITIVE**: +40 pp (hexanal), +35 pp (2-nonanone) above independent-Langmuir; ≈0 for nonanal (saturated). **Never compose saliva and matrix binding as independent terms** | — | Barallat-Pérez 2024 |
| B6.7 | **aldehyde binding constants are METHOD-DEPENDENT by 7–35×** where the ketone gap is 1.3×, for a mechanism the source states (cysteine–aldehyde condensation / Schiff base at pH 6–10, counted by headspace depletion, suppressed by 2-mercaptoethanol and invisible to dialysis) | hexanal 35×, nonanal 7.2×, 2-nonanone 1.3× | Barallat-Pérez vs Damodaran |
| B6.8 | **three proteins, three methods, three labs, 43 years agree within 1.8× on 2-nonanone** ⇒ the **ketone** per-gram constants ARE transferable | 3.72 / 4.71 / 6.63 × 10⁻² L/g | Damodaran, Barallat-Pérez, Andriot |
| B6.9 | **7S carries essentially all of soy's carbonyl binding; 11S almost none** — a composition confound that explains the 3–8× Snel/Damodaran spread without either being wrong | — | Damodaran 1981b (abstract only) |

## B.7 — ★ THIOL CONSUMPTION: THE STRUCTURAL FINDING

**Four papers, four temperatures, four dominant channels, and each excludes the others'.**

| T | dominant sink | what EXCLUDES the neighbouring channel |
|---:|---|---|
| **30 °C** | covalent thioether addition to melanoidins/electrophiles | **the disulfide branch is < 1.5 % of the flux** (Hofmann 2002, measured) |
| **50 °C** | acid-catalysed electrophilic **oligomerisation at C-5** | **air ≈ argon**, and the mass balance **fails to close** as thiol + disulfide (van Seeventer 2001) |
| **115 °C** | **oxidative dimerisation + radical coupling to MeSH** | dimer carries **up to 49 %** of the MFT pool — the channel the 30 °C system rules out at < 1.5 % |
| **120 °C** | dimerisation, **6.5–9.6 %** of MFT-equivalents, **pH-invariant over pH 6–8** | independent corroboration of the 115 °C channel, different lab/feedstock/pH/buffer |

> **⇒ IMPLEMENT THIOL CONSUMPTION AS A SET OF NAMED CHANNELS, EACH WITH A DECLARED
> TEMPERATURE RANGE OF VALIDITY AND AN EXPLICIT `no_Ea_available` STATE — NOT AS ONE LUMPED
> ARRHENIUS SINK. Pairing any two of the four rates to extract an Ea is a PROHIBITED DERIVATION
> (§C.1), and it would be worse than the two-point error K1 already refused, because the two
> rates are not even the same reaction.**

**And the aroma qualification (§A.7.3): the dimer is 15.6× more potent than the monomer, so
its OAV matches the monomer's. Mass lost to dimerisation is not aroma lost.**

## B.8 — Model-diagnostic constraints (about the model, not the chemistry)

| # | constraint | source |
|---:|---|---|
| B8.1 | **the disjointness error is TWO-SIDED and therefore CANNOT be a barrier error**: **2.76× TOO LOW** on fed norfuraneol and **≥17.7× TOO HIGH** on in-situ norfuraneol — a **≈49× two-sided spread** | Z1 §D |
| B8.2 | the model's response is **SATURATED**: a 2.30 kcal/mol barrier change should move the rate **15.9×** at 418 K; the Table-4 observable moved **1.22×** and the NF channel ≈**2.65×** — **13× and 6× weaker than Eyring** | Z1 §E(iii) |
| B8.3 | ⇒ **chase the saturation BEFORE touching any barrier. Until then every fit against Table 4 is fitting a knob that barely moves the dial** | Z1 recommendation 2 |
| B8.4 | the concentration difference between the two systems is the right size to explain the spread: NF fed:in-situ **4.2×**, H₂S plausibly another **3–10×** ⇒ **13–42×** if first order in each | Z1 §E(i) |
| B8.5 | the source author **REFUSED to transfer** his own model to real food: *"these mechanistic studies in model systems are not easily applied to real food systems … the parameters are only applicable for specific experimental conditions such as time–temperature profile of frying, potato genotype, slice thickness and initial concentration of precursors."* | Knol 2009, Z2 #11 |
| B8.6 | **composition moves 2-pentylfuran; process severity does NOT** — two papers, one clean dissociation | Conti 2025b (0.26 vs 0.24 µg/g, n.s.) vs Xin 2026 (RI −46 %, GL −30 %, FR −21 %) |
| B8.7 | the **ramp-vs-hold** question resolves toward **HOLD**: an independent industrial lab reproducing the Hofmann protocol writes *"heated in an autoclave from room temperature to 130 °C in 10 min **and was kept at 130 °C for 20 min**"*. ⇒ **the nine live `hofmann1998_*` benchmarks are NOT 8× over-exposed** | van Seeventer 2001, Z2 #6 |

## B.9 — ★ NEW THIS WAVE: quantification validity — protein-blend and method constraints

| # | constraint | measured | source |
|---:|---|---|---|
| B9.1 | **★ single-IS HS-SPME µg/kg are NOT absolute concentrations, proved INTERNALLY to the literature.** The same physical sample, reported by the same group in two 2026 papers, differs by **≈10× (hexanal), ≈17× (2-pentylfuran), ≈23× (nonanal) with NON-OVERLAPPING ranges** — while agreeing to **3 %** on 2-heptanone. The second paper never states its SPME fibre | see §E.2.2 | **Xin 2026 vs Xin 2026b** |
| B9.2 | protein source moves total aldehydes **12.4×** at constant total protein and **zero added carbohydrate** | ERPI 2398.22 vs ENPI 193.95 µg/kg | Xin 2026b |
| B9.3 | but it moves **2-heptanone only 1.10×** across the same six formulations | — | Xin 2026b |
| B9.4 | ⚠️ **β-sheet content and fibrous degree RANK DIFFERENTLY** — EMPI 32 % is top on β-sheet, ENPI 1.59 is top on fibrousness — **contradicting the paper's own stated mechanism and its abstract** | β-sheet EMPI 32 > EPPI 31 > ENPI 27 > ERPI/EYP 22 > ESPI 20 %; fibrous ENPI 1.59 > EMPI 1.27 > EPPI 1.24 | Xin 2026b |
| B9.5 | total free amino acids span **1613.74 → 3347.27 µg/g** across six protein blends — **the amine precursor pool varies 2.1× at constant total protein** | validated HILIC-MS/MS | Xin 2026b |

**★ B9.1 IS A COLLECTION-LEVEL FINDING, NOT A PAPER-LEVEL ONE.** Apply it as a flag —
`absolute_concentration: false` — to **every** semi-quantitative HS-SPME row in the corpus:
Conti 2025 I & II, both Xin papers, Zhang 2024, Zhou 2023. None may be scored against a SIDA
anchor or against each other on absolute level. **Within-paper ratios are unaffected**, which is
why every "RATIO-ONLY" verdict in §A now has a citable, literature-internal basis rather than
resting on the corpus's own assertion.

## B.10 — ★ NEW VIA Z3: route mix, branch fractions and the formation/destruction split

| # | constraint | measured | source |
|---:|---|---|---|
| **B10.1** | **★★ THE BRANCH FRACTION IS A FUNCTION OF CONCENTRATION, MEASURED.** A **2× change in precursor loading moves the xylose share of MFT from 15 % to 46 %** — a 3.1× change in the branch fraction, one lab, one method, one pH, one temperature | 85 : 15 (1×) vs 54 : 46 (2×) % thiamine : xylose | Cerny 2007 Table 5, pH 5.00, 145 °C / 20 min |
| **B10.2** | **and BOTH single-route controls are run**: MFT is **> 99 % thiamine-derived with no cysteine** and **> 95 % xylose-derived with no thiamine**, against **54 : 46** in the full ternary | — | Cerny 2007b Tables 3 / 5 / 6 |
| B10.3 | **FFT is > 99 % XYLOSE-derived at every pH and in every control**, including the no-cysteine control — **⇒ thiamine can serve as the sulfur donor for FFT** | — | Cerny 2007 T4, Cerny 2007b T5/T6 |
| **B10.4** | **★ Ea(destruction) exceeds Ea(formation) by ≳ 70 kJ/mol** for a fed heterocyclic intermediate — a **measured sign reversal on the same 75 → 100 °C step: ×4.2 at 10 min but ×0.177 at 60 min** | Ea_form ≳ 97, Ea_dest ≳ 166 kJ/mol | Hofmann & Schieberle 1996b Fig. 2 **[D]**, 2-acetyl-2-thiazoline |
| B10.5 | ⚠️ **B10.4 is on 2-acetyl-2-thiazoline, NOT on MFT**, and it is a fed-intermediate experiment. *"It does not measure MFT's temperature response and must never be cited as if it did. What it licenses is the MODEL CHANGE …, not a PARAMETER TRANSFER."* | — | `hofmann1996b_extraction.md`, verbatim |
| B10.6 | **★ oxygen is REQUIRED**: removing it costs **5.8×** of the 2-acetyl-2-thiazoline yield (9.8 → 1.7 µg); Cu²⁺ buys **2.1–3.9×** | argon 1.7 vs air 9.8 vs Cu²⁺ 20.1 / 38.4 µg | Hofmann 1996b Table 4 |
| B10.7 | **the norfuraneol Ea is not a constant**: across three studies it spans **64–122 kJ/mol** and is a function of **precursor loading and matrix** | Lau 2003 64.0–122.3; Pandit 2006 81.4–96.1; Bornhorst 104.9–122.3 | `bornhorst2017b_extraction.md` |
| B10.8 | **norfuraneol plateau FALLS as precursor loading RISES** — 0.47 → 0.24 → 0.10 mg/g for 1_R0.5_L → 1_R1_L → 2_R2_L at 90 °C, while **k rises** 5.1 → 7.4 → 18.3 × 10⁻³ min⁻¹ | — | Bornhorst 2017 Table 1 |
| B10.9 | ⚠️ **B10.8's temperature counterpart is INSIDE ONE STANDARD ERROR** — *"the 1_R,0.5_L trend is 0.57 ± 0.17 → 0.43 ± 0.02, so the fall is inside one standard error at 80 °C; the other two formulas are flat within noise. **This is suggestive corroboration, not a measurement of a sink. State it that way or not at all.**"* | — | `bornhorst2017b_extraction.md`, verbatim |
| **B10.10** | **★ MFT peaks at pH 6.0 in a THIAMINE-containing system**, while FFT and 2-furaldehyde fall **monotonically to exact zero** and three H₂S-consuming heterocycles **switch on only above pH 5.5** | MFT 415/336/371/**539**/391; FFT 431/368/364/185/**0**; furfural 208/158/165/**0**/**0**; thiophenone **0**/0/0/372/470 (peak areas ×10⁶) | Cerny 2007 Table 2, pH 4.0/5.0/5.5/6.0/7.0 |
| B10.11 | ⚠️ **B10.10 is a THIRD, DIFFERENT MFT-vs-pH shape** — Hofmann 1998 is monotone-decreasing over pH 3→7; Zhou 2023 peaks at pH 7; Cerny 2007 peaks at pH 6.0. **Three papers, three shapes, three feedstocks.** *"Any repo pH term applied to 'MFT' without knowing which route supplies it is fitting two different functions at once."* | — | §B2.1, §B2.4, and Cerny 2007 |
| B10.12 | **the caramelisation (amine-independent) lane carries 25–80 % of A₂₉₄ and 7–55 % of A₄₂₀**, and **the share RISES with temperature** | — | Ajandouz 2008 §3.4 |
| B10.13 | **2,5-dimethylpyrazine's BOTH methyl carbons come from the SUGAR** — 0 % ¹³C label from [3-¹³C]Ala or [2-¹³C]Gly in all four systems — while **2,6-dimethylpyrazine picks up 25–30 % from the amino acid** | 0/0/0/0 vs 0/0/30/25 % labelled | Amrani-Hemaimi 1995 Table 2 |
| B10.14 | **3-ethyl-2,5-dimethylpyrazine is 100 % ¹³C-labelled from [3-¹³C]alanine and ND from glycine** ⇒ *"one single reaction route exists"*, and the compound is **present with alanine (20 / 19 %) and absent with glycine (0 / 0 %)** — an **on/off amino-acid switch** | — | Amrani-Hemaimi 1995 Tables 1 & 2 |
| B10.15 | **the FFT branch is ≲ 1.2 % of the furfural flux** in water at 65 °C — so a "furfural sink" Ea is overwhelmingly an Ea for something else | < 1 % FFT conversion against 84.9 % furfural conversion | Yaghmur 2005 §3.1 **[M]** + **[D]** |
| B10.16 | **the furfural sink is FIRST ORDER in cysteine**, linear over 0–12 mmol, r² = 0.98 ⚠️ **microemulsion arm only** | — | Yaghmur 2005 Fig. 13 |
| B10.17 | **an interfacial matrix accelerates the furfural sink 3.7–6.3×** and **shifts its Ea from 46.50 to 32.3 kJ/mol**, while a 50 % w/w polyol *slows* it and raises the Ea to 56.72 | k_obs water 3.2–16.0 vs O/W 20.1–59.8 vs PG 1.7–12.2 (h⁻¹ ×10²) | Yaghmur 2005 Table 1 |
| B10.18 | **a structural zero**: M-2 (norfuraneol) = **0** in all three food matrices with no added precursors | — | Bornhorst 2017 §3.1 |

**★ B10.1 IS THE MOST CONSEQUENTIAL SINGLE ROW Z3 ADDS.** Wave Z1's §E argued that the
Hofmann-vs-Cerny disjointness is a *concentration-handling* defect rather than a barrier error,
and had to estimate the size of the effect (13–42×) from an assumed first-order rate law.
**Cerny 2007 Table 5 measures it directly.** Any model carrying a fixed thiamine : sugar MFT
split — or, by extension, any fixed branch fraction anywhere — is falsified by one row pair.


---

# §C. DECLARED GAPS — carried VERBATIM

**These are quoted, not paraphrased. Do not soften them.**

## C.1 — Thiol consumption has no activation energy — K1, verbatim

> "**No activation energy for thiol consumption exists anywhere in this basket.** Both thiol
> papers are single-temperature. Any T-extrapolation of the FFT/MFT sink is the repo's own
> assumption and must be labelled as such."
> — `k1_kinetic_parameters.md` §0

> "Ea(MMBF loss) ≈ **30 kJ/mol** … **[K1] — DO NOT INGEST.** Two labs, two matrices, two
> analytical bases, two points. The same treatment applied to FFT gives a **negative** Ea (the
> 80 °C brew is *slower* than the 30 °C models, because the real brew's electrophile pool was
> partly consumed during extraction). **The literature reviewed supports NO activation energy
> for thiol consumption.**"
> — `k1_kinetic_parameters.md` §1d

**K3 extends this** (`Zhang2024_extraction.md` §K3.7): *"**IT DOES NOT SUPPLY AN ACTIVATION
ENERGY FOR THIOL CONSUMPTION, AND IT CANNOT BE MADE TO.** … Pairing Zhang's 115 °C rate with
K1's 30 °C rate to extract an Ea would repeat exactly the two-point, two-lab, two-matrix,
two-method error that K1 named and refused — and it would be worse, because the two rates are
not even the same reaction. Record this as a **prohibited derivation**."*

## C.2 — No same-method matrix-vs-water threshold pair — K2, verbatim

> "**NO SAME-METHOD MATRIX-VS-WATER PAIR EXISTS ANYWHERE IN THESE TEN PAPERS.** Every ratio in
> §A.8 crosses labs, decades and criteria… Vega asserts its water references were obtained
> '**using similar sensory techniques**'; that is **not verifiable from the paper** and is
> contradicted by its own methods."
> — `k2_matrix_and_thresholds.md` §D.2(i)

> "**They are not even the same kind of quantity as the aqueous values they are compared to.**
> … **Brewer's beef 'thresholds' are nominal doses added BEFORE a 70 °C cook** — a large part of
> the 100–6 700× is thermal loss, not perception."
> — `k2_matrix_and_thresholds.md`, four-line answer, item 2

> "**The repo cannot obtain matrix-relevant odour activity from `f_free`, and the gap is not
> closable by tuning `a_p`.**"
> — `k2_matrix_and_thresholds.md` §B.4

## C.3 — No temperature dependence for MFT — Zhang 2024, verbatim

> "**The paper supplies the MECHANISM and its magnitudes. It does NOT supply a temperature
> series for MFT, and it contains no activation energy, no Arrhenius fit, and no rate constant
> for any MFT-consuming step.** Every kinetic number in it is at **115 °C**, one temperature."

> "**The honest one-line answer: this paper supplies the mechanism and the sign, not the
> magnitude. It cannot close Wave U's ×0.249, and any wave that claims it does is over-reading
> the SI it does not have.**"
> — `Zhang2024_extraction.md`

## C.4 — No precursor dose–response anywhere in the extrusion lane — K2, verbatim

> "**No dose–response for any precursor.** Conti uses one thiamine level (1.5 %); Xin uses one
> carbohydrate level (6 %); Guo uses one flavour level (1 %). **No sensitivity of yield to
> precursor loading is recoverable from any of the three.**"

> "**No zero-precursor control in Conti.** Every 'thiamine-derived' assignment is
> literature-attributed (Dreher et al. 2003, a model *orange juice* system), not demonstrated.
> And the largest class in Conti's table, the pyrazines, is **not** a thiamine product at all."

> "**No absolute quantification, except Conti's hexanal.** … Xin's 2-pentylfuran at
> **50–100 ppm** is a DVB/CAR/PDMS partition artefact, not a concentration."
> — `k2_matrix_and_thresholds.md` §C.4

## C.5 — Quan 2020: 88 constants, all unusable — K1, verbatim

> "**Four disqualifying defects: (1) NO UNITS on any k; (2) NO RATE LAWS — Supplementary
> Figure 1 is absent from the PDF, so every reaction ORDER is unknown; (3) the 'Table 1
> Continued' headers are MISLABELLED; (4) no ± SD printed despite a footnote claiming them.**
> Authors state explicitly: '**specific activation energies cannot be estimated**'."
> — `k1_kinetic_parameters.md` §2d

## C.6 — Acrylamide degradation is unidentifiable, not mis-valued — Z2, verbatim

> "Knol 2005: k6 fitted but '**the model was not restrained by experimental data** for the
> products formed in the degradation reaction'. Knol 2009: **every** degradation parameter has
> SD ≥ estimate (`tc2 = 0.60 ± 52`, `k2 = 3.5 ± 8`, `τ = 22 ± 15`). Knol 2010 (Z0): the Ea
> **went negative and was deleted**."
> — `z2_index.md` finding 12

## C.7 — No number to fit a barrier to, for the 1,4-dideoxyosone route — Z1, verbatim

> "**Nedvidek 1992 states no reaction temperature** for any of its aqueous runs, and no numeric
> measurement of any 1,4-dideoxyosone exists in it. The α-amino-acid dependence is stated
> qualitatively only. **There is no number here to fit a barrier to.**"
> — `z1_disjointness_verdict.md`, residual uncertainty 6

## C.8 — Meynier and Cerny: no absolute quantification of any kind — Z0, verbatim

> "**Directional/ratio benchmarks: YES (~20 rows, with fold ratios). Absolute ppb: NEVER** —
> 'response factors not determined', µg/10 mg ribose is semi-quantitative."
> — `z0_index.md`, meynier1995 verdict

> "**Directional/topology only (8 rows), plus one step-level branching ratio (54:46). No
> quantification of any kind exists in this paper.**"
> — `z0_index.md`, cerny2004 verdict

## C.9 — Frankel: the brief's premise was wrong — Z2, verbatim

> "**NOT a yield source — the brief's premise is wrong; no absolute yield and no Ea exist in it
> (one temperature, 180 °C).**"
> — `z2_index.md`, frankel1989 verdict

## C.10 — van Seeventer: zero order, single temperature — Z2, verbatim

> "**Benchmark: the precursor-conversion pair (55 % / 75 %) and Table 1's degradation rates —
> but they are ZERO-ORDER, % of initial per day; do not convert to first-order half-lives.**
> **Single temperature (50 °C): no T-dependence is ingestable.**"
> — `z2_index.md`, vanseeventer2001 verdict

## C.11 — ★ NEW THIS WAVE: Zhou 2023's pH labels are initial, not held

> "⚠️ **THIS IS THE MOST IMPORTANT METHODOLOGICAL NUMBER IN THE PAPER.** The 'pH 6 / 7 / 8'
> labels on Table 1 are **initial** pH of an **unbuffered** system. By the end of the 60-min
> hold… The pH 6 and pH 7 runs **converge to within 0.2 units of each other**. … **If the model
> has no pH-trajectory state, Zhou's Table 1 must be ingested as an ORDINAL/DIRECTIONAL
> constraint, not an absolute-pH one.**"
> — `zhou2023_extraction.md` §3

## C.12 — ★ NEW THIS WAVE: Zhang 2024's rate-constant units cannot be recovered

> "**VERDICT ON 0.0028 / 0.0031: DO NOT INGEST AS NUMBERS.** … the underlying experiment's
> duration and temperature are unstated and its regression table is absent. **What IS ingestable
> is the RATIO 0.0031/0.0028 = 1.107 … and the two-stage structure**, both of which are
> unit-independent."
> — `Zhang2024_extraction.md` §K3.3

## C.13 — ★ THE ALKALINE-pH WALL — Z3, verbatim. **The single biggest caveat on §A.1.1 and §A.3.6(i).**

**Ajandouz (pH 8.0 / 9.7) and both Bornhorst papers (pH 8.4–9.5) carry nearly every new constant
Z3 contributes, and all three dossiers refuse rate transfer to the repo's pH 5–7.**

> "⚠ **pH 8.0 and 9.7 are ALKALINE.** The repo's systems are pH 5–7; Wave S3's fitting corpus is
> at **pH 6.8**. This is the single biggest caveat on everything below."
> "⚠ **Proteins (BSA, casein), not free amino acids.** The amine is a protein-bound lysine
> ε-NH₂ … not free glycine."
> "**Do not transfer any rate.**" / "**the rate transfer is not defensible at all.**"
> "**The 'Ea is pH-insensitive' licence applies to the glucose-loss and amino-loss steps ONLY,
> not to browning.**"
> "Any RATE transfer to pH 5–7 | **NOT LICENSED — 3.5–50× pH effects measured**"
> "**Two of sixteen parameters get a referee; fourteen do not.** Say that plainly rather than
> implying the fit has been validated."
> — `ajandouz2008_extraction.md`

> "⇒ **Adding free lysine raises the mashed potato from pH 5.2 to 8.4–9.5.** … **Neither the
> rates nor the M-2∞ values transfer to a pH-5 sulfur benchmark.**"
> "★ **THIS IS AN APPROACH-TO-PLATEAU FORMATION LAW, NOT A DEGRADATION LAW.** … Anyone reading
> these `k` values as norfuraneol *degradation* rates would be inverting the paper."
> "`k` read as a norfuraneol DEGRADATION rate | **FORBIDDEN**"
> — `bornhorst2017_extraction.md`

> "★ **Across three studies the norfuraneol Ea spans 64–122 kJ/mol and is a function of precursor
> loading and matrix, not a constant.** A repo that reads any single value as 'the norfuraneol
> barrier' is over-fitting a lump."
> "**MUST carry 'approach-to-plateau accumulation in an alkaline mashed-potato gel', never 'the
> norfuraneol barrier'**"
> "Transfer of any RATE to pH 5–7 | **NOT LICENSED — pH 8.4–9.5**"
> "Extrapolation of these Ea above 100 °C | **NOT LICENSED — 80–100 °C, 3 points**"
> — `bornhorst2017b_extraction.md`

## C.14 — The Yaghmur Ea is a lump over ≳98.8 % non-FFT flux — Z3, verbatim

> "**IT IS NOT AN Ea FOR FURFURAL → FFT. It is an Ea for the LUMPED DISAPPEARANCE OF FURFURAL in
> the presence of cysteine, and in the aqueous arm FFT accounts for LESS THAN 1.2 % of the
> furfural consumed.** … **The lumped Ea of 46.50 kJ/mol therefore describes a channel the paper
> never identifies.**"
> "⚠ **40–70 °C.** Every processing temperature in this repository (115–180 °C) is **outside**
> this range by 45–110 °C. **Nothing here licenses extrapolation to a Maillard processing
> temperature.**"
> "⚠ **THIS IS NOT A CONSTANT REFUTATION AND MUST NOT BE FILED AS ONE.** Three reasons, all real:
> (i) the measured quantity is a *lump* over ≳98.8 % non-FFT flux, not either elementary step;
> (ii) an Arrhenius Ea and an Eyring ΔG‡ are different quantities and the repo's barriers are the
> latter; (iii) the measurement is at 40–70 °C, far outside the repo's operating window."
> — `yaghmur2005_extraction.md`

## C.15 — Three Z3 papers have NO absolute quantification at all — verbatim

> "★★ **THE ANALYTICAL LIMIT, STATED FIRST: there is NO absolute quantification anywhere in this
> paper.** … **No internal standard. No response factors. No µg, no ppb, no mol %. One experiment
> per condition — no replicates.** ⇒ **Directional/ratio and isotope source ONLY.**"
> — `amrani-hemaimi1995_extraction.md`

> "⚠ **THE QUANTIFICATION CAVEAT, IN THE AUTHORS' OWN WORDS (p. 1554):** 'The data in Table 2
> represent **integrated peak areas** and are **not corrected by MS response factors** or taking
> into account the **different partition coefficients** … Also, **for certain compounds, the pH
> might have an effect on the partition coefficient**.'"
> "Absolute concentrations of anything | **NOT AVAILABLE — peak areas only**"
> "**Any repo pH term applied to 'MFT' without knowing which route supplies it is fitting two
> different functions at once.**"
> — `cerny2007_extraction.md`

> "Any concentration or yield | **NOT AVAILABLE — isotope ratios only, no absolute data of any
> kind**" · "the control reactions C–F are **UNREPLICATED**"
> — `cerny2007b_extraction.md`

**And the general hazard, flagged across three dossiers:** *"MFT/pyrazine 'amounts' in
Amrani-Hemaimi and Cerny & Briffod are peak areas, not concentrations. … Any downstream row that
carries them in µg, ppb or mol % has been fabricated in transit."*

## C.16 — van Boekel 2005 is a dead end, and it is NOT a repo defect — verbatim

> "⚠ **A NEGATIVE RESULT: THE PDF IS A ONE-PAGE POSTER ABSTRACT WITH ZERO PARAMETERS** …
> **Zero rate constants. Zero activation energies. Zero enthalpies or entropies. Zero tables.
> Zero figures. Zero concentrations. Zero melanoidin parameters.** - The word **'melanoidin' does
> not appear.** … The parameter table the brief hoped for is described in the **FUTURE TENSE** …
> **This is a statement of intended work, not a result.**"
> "**Not a repo citation defect.** … **This is a dead-end RETRIEVAL, not a defect.**"
> "⚠ **SELF-CORRECTION** … **Any claim that 'no paper in the corpus models melanoidins' is wrong;
> cite Martins 2005.**"
> — `vanboekel2005_extraction.md`

## C.17 — The van Seeventer 130 °C sink CANNOT be derived — Z3, verbatim

> "**⇒ All five degradation rates (MB 14, MFT 59, FFT 28, HDF < 10, sotolone n.d. %/day) are
> STORAGE rates at 50 °C. No processing-temperature rate exists in this paper.**
> **⇒ No activation energy can be derived. One temperature, no Arrhenius, no z-value, no Ea.**
> That is not a limitation of the extraction; it is a limitation of the experiment."
> "⚠ **EVERYTHING IN THIS SECTION IS MY ARITHMETIC, NOT THE PAPER'S.** The paper licenses none of
> it. I am running it precisely to establish that it **cannot** be run, which is a result."
> "**van Seeventer's storage rate is NOT the sink that explains Wave U, and no choice of Ea can
> make it one.**"
> "**Do not bolt van Seeventer's 59 %/day onto the network as the MFT sink.**"
> — `vanseeventer2001_z3_addendum.md`

**What the exercise DOES yield, as a bounded target rather than a parameter:** solving
`−0.5·k₁₃₀ + 4·k₁₀₀ = ln(0.0547)` for Wave U's `100 °C/4 h → 130 °C/0.5 h` hold-out gives
**Ea ≈ 266 / 175 / 124 / 116 kJ/mol for an assumed k₁₀₀ of 0.01 / 0.10 / 0.50 / 0.70 h⁻¹** —
**a SENSITIVITY, not a determination** (k₁₀₀ is a free assumption). Back-extrapolating those to
50 °C gives **~1e-6 to ~5 %/day, i.e. 12× to ~10⁷× slower than the measured 59 %/day** ⇒ **the
processing sink and the storage sink are DIFFERENT CHANNELS.** That qualitative conclusion is
the robust result; **no single row of the table is a constant.**

## C.18 — Corpus-level gaps still open

| gap | consequence | status |
|---|---|---|
| Xin 2026 **Tables S2 and S3** (all 45 volatiles × 9 formulations, + the OAVs **with the threshold list**) | would turn ~12 quoted compounds into **405 quantified cells** | ⚠️ **partially addressed this wave — see §A.7.1 and the `Xin2026b_extraction.md` jackpot check in §E** |
| **Zhang 2024 Supplementary Information** (Fig. S1 dynamic profiles, **Fig. S2 the temperature sweep**, Tables S6/S7 the kinetics, Table S5 the standard curves) | **holds every time series and the entire temperature axis** the sulfur-T question turns on | **OPEN — the single highest-value retrieval in the corpus** |
| Tian et al. 2019 (`10.3168/jds.2019-16796`) | the only way to settle the `?/kg` unit; without it **seven milk thresholds carry a factor-of-1000 basis risk** | OPEN |
| Quan 2020 **Supplementary Figure 1** + version of record | would unlock 88 fitted constants | OPEN |
| Hemmler 2018 **SI Figs S1–S5** | the six-sugar comparison | OPEN |
| Damodaran & Kinsella **1981b** (only p. 1 held) | would close the **7S/11S composition confound** — the most plausible non-error explanation of the Snel/Damodaran 3–8× gap | OPEN |
| Buttery, Ling & Juliano 1982 | **the only route to a 2-AP odour threshold** | OPEN |
| Hwang's pyrazine twin (JAFC 1995, 43, 179−184) | — | OPEN |
| **an odour threshold for 2-methylthiophene, 2,3-butanedione, thieno[3,2-b]thiophene, 2-methylthiazole, the furfuryl sulfides** | 7 of Zhou's 22 compounds cannot be scored | OPEN |

---

# §D. PROPOSED FIT vs HOLD-OUT SPLIT — **DRAFT FOR ORCHESTRATOR**

> **STATUS: DRAFT. NOT APPLIED. No repo file was touched. This section is a proposal for the
> orchestrator to accept, amend or reject.**

## D.0 — The three rules this split obeys

1. **No dataset appears in both columns.** Verified by construction — each source paper is
   assigned exactly one role per module, and where a paper spans modules the *rows* are
   partitioned by module and the partition is stated.
2. **The existing frozen hold-outs are untouched.** Enumerated in D.1; none is moved, renamed,
   re-scoped or promoted to fit.
3. **Every module holds out at least one dataset.** Checked in D.7.

## D.1 — THE FROZEN HOLD-OUTS — do not touch, do not re-scope, do not promote

**21 bundles, all `evidence_class: external_validation_only`.**

**Matrix path — `data/benchmarks/external_validation/*.json` (4):**
`external_validation_bi_2020_raw_pea_hexanal` · `external_validation_bi_2020_roasted_pea_hexanal`
· `external_validation_li_2026_spi_wg_hme_control` ·
`external_validation_liu_2023_ppi_offnote_baseline`

**Maillard path — `data/benchmarks/external_validation/maillard_path/*.json` (17):**
`mp_holdout_fructose_asparagine_180C_Lin2022` ·
`mp_holdout_glucose_alanine_130C_2h_pH50_Schibilsky2019` ·
`mp_holdout_glucose_alanine_130C_2h_pH80_Schibilsky2019` ·
`mp_holdout_glucose_asparagine_180C_10min_Chang2021` ·
`mp_holdout_glucose_asparagine_180C_30min_Chang2021` ·
`mp_holdout_glucose_asparagine_180C_30min_water_Chang2021` ·
`mp_holdout_glucose_asparagine_180C_Ye2024` ·
`mp_holdout_glucose_only_autoclave_121C_Steinhagen2021` ·
`mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH3` ·
`mp_holdout_hofmann1998_glucose_cysteine_145C_20min_pH7` ·
`mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3` ·
`mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7` ·
`mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5` ·
`mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026` ·
`mp_holdout_ribose_cysteine_buffer_110C_2h_Yiltirak2026` ·
`mp_holdout_ribose_cysteine_buffer_120C_1h_Yiltirak2026` ·
`mp_holdout_ribose_cysteine_buffer_130C_30min_Yiltirak2026`

⚠️ **Note the existing Hofmann 1998 partition and PRESERVE it:** the **pH 5.0** rows are in the
FIT panel (`data/benchmarks/hofmann1998_*.json`); the **pH 3 and pH 7** rows and the **xylose**
row are HELD OUT. **Any new Hofmann-1998-derived row must respect that same axis-based cut:
pH 5.0 → fit, off-pH → hold-out.** This is the cleanest existing convention in the repo and the
split below extends it rather than inventing a new one.

## D.2 — THE EXISTING FIT PANEL (for reference; unchanged by this proposal)

23 files in `data/benchmarks/`: `acrylamide_spi_extrusion_130C_ACSRef3` ·
`cml_cel_commercial_pbma_Foods2023` · `cys_ribose_140C_Hofmann1998` ·
`furosine_extrusion_crossover_140C_RamirezJimenez2000` ·
`hofmann1998_c2c3_recombination_145C_20min_pH3 / pH5 / pH7` ·
`hofmann1998_fructose_cysteine_145C_20min_pH5` ·
`hofmann1998_furan2aldehyde_h2s_145C_20min_pH5` ·
`hofmann1998_glucose_cysteine_145C_20min_pH5` ·
`hofmann1998_norfuraneol_cysteine_145C_20min_pH5` ·
`hofmann1998_norfuraneol_h2s_145C_20min_pH5` ·
`hofmann1998_ribose_cysteine_145C_20min_pH5` · `pea_isolate_40C_PratapSingh2021` ·
`pea_isolate_ribose_cysteine_100C_45min_Internal2026` ·
`pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026` *(diagnostic_only)* ·
`pea_isolate_uht_140C_Trikusuma2019` · `resconi_2023_pbma_beef_identity_benchmark` ·
`soy_isolate_40C_PratapSingh2021` · `soy_isolate_ribose_cysteine_100C_45min_Internal2026` ·
`soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026` *(diagnostic_only)* ·
`thiamine_cys_glucose_120C_Bolton1994` · `thiamine_cys_xylose_145C_Cerny2008`

## D.3 — MODULE 1: SULFUR FORMATION

| dataset | **role** | one-line reasoning |
|---|---|---|
| Hofmann 1998 T1/T2 **pH 5.0** rows (ribose, glucose, fructose, xylose already split) | **FIT** *(already)* | the corpus's only SIDA absolute anchors at the repo's reference pH; unchanged |
| Hofmann 1998 T1/T2 **pH 3 / pH 7** rows | **HOLD-OUT** *(already frozen)* | the pH axis is the model's weakest (6/10); keeping it out is the only honest test of the pH term |
| Hofmann 1998 T3/T4/T6/T7/T8/T10 (step-level fed-precursor) | **FIT** | step-level constraints are what the network is *for*; they identify individual barriers rather than lumps |
| Hofmann 1998 T2 **dry 180 °C / 6 min** rows | **HOLD-OUT (new)** | a genuinely different regime (silica, 6 min, 180 °C); the only a_w extrapolation test the sulfur branch has |
| **Zhou 2023 Table 1, pH 7 column (ARP and ARP+Cys)** | **FIT (new)** | the fed-Amadori initial condition the corpus has never had; pH 7 is its midpoint and its MFT maximum |
| **★ Zhou 2023 Table 1, pH 6 and pH 8 columns** | **★ HOLD-OUT (new)** | **the strongest new hold-out available.** A model fitted at the pH-7 maximum must *predict* the fall on both sides — and it can only do so if the two-factor law (§B2.4 / §A.3.3) is structurally present rather than fitted |
| **Zhou 2023 Figure 3, Cys+MGO time × pH grid** | **HOLD-OUT (new)** | a completely different feedstock (α-dicarbonyl, not Amadori) and the only time axis; nothing in the fit panel constrains it |
| **Zhou 2023 §A.3.3 cross-system pair (582 vs 665 µg/L)** | **HOLD-OUT (new)** | a *derived consistency* test, not a level; it scores the ARP → MGO flux with no free downstream parameter |
| Whitfield 1999 (fed NF, pH 4.5, 140 °C) | **FIT** | the only independent replication of the Hofmann NF chemistry; needed to pin the NF channel |
| Whitfield 2001 (fed NF, **pH 6.5**) | **HOLD-OUT** | the ≥150× pH collapse is the sharpest single pH prediction the model can be asked to make ⚠️ its H₂S column should NOT carry a standalone row (Methods omits HMF — likely a typo) |
| Cerny 2003 T2/T3 (isotope, 95 °C / 4 h) | **HOLD-OUT** | 95 °C is 50 °C below the fit panel and the route mix is measured to change there; the NF ≤ 7 % ceiling must be **evaluated at Cerny's conditions** |
| Cerny 2004 (in-situ branching 54:46 etc.) | **FIT** | pure ratios, free of the fed-precursor artefact; they identify topology, not magnitude |
| van Seeventer 2001 precursor conversion (55 % / 75 %) | **FIT** | a *reactant*-side constraint; it cannot be traded against any product-side residual, so it adds information without competing |
| Meynier 1995 (~20 directional rows) | **HOLD-OUT** | directional only, and the pH axis again; excellent as a shape test, useless as a level |
| Shu 1988 | **neither** | GC area %, no internal standard; pH-shape source only. Record as a shape constraint (§B2.2), not a scored dataset |

## D.4 — MODULE 2: THIOL CONSUMPTION

| dataset | **role** | one-line reasoning |
|---|---|---|
| Charles-Bernard 2005 T2 ladder, **25 °C** | **FIT** | ten first-order constants spanning 755× at one well-specified condition — the module's calibration backbone |
| Charles-Bernard 2005 T3 attenuation matrix | **FIT** | branch weights, not levels; they identify the channel decomposition the ladder alone cannot |
| Hofmann 2002 Fig. 6 + T2, **30 °C** | **FIT** | the cross-validating partner (agreement within 25 %); fitting one without the other wastes the corroboration |
| **Hofmann 2002 Fig. 1, real coffee brew, 80 °C** | **★ HOLD-OUT** | the only temperature extrapolation available — **and it is the one the model will get WRONG in the informative direction** (the brew is *slower* than the 30 °C models). A pass here would be strong evidence; a fail localises the electrophile-pool depletion term |
| **van Seeventer 2001 Table 1, 50 °C zero-order MFT/FFT** | **★ HOLD-OUT** | a third temperature, a third mechanism, and **zero order** — it tests the *functional form*, not just the magnitude |
| **Zhang 2024 Figs. 2d/e/f dimer & MMFT fractions, 115 °C** | **HOLD-OUT** | the highest temperature available and a channel the 30 °C data explicitly exclude; the sharpest test of the "named channels, not one Arrhenius" architecture |
| Zhang 2024 Fig. 1 (four sulfur additives, redox-state axis) | **FIT** | one condition, four levels — usable to identify the redox term that nothing else in the corpus constrains |
| **Zhou 2023 dimer/monomer fractions, 120 °C, pH 6–8** | **HOLD-OUT** | independent lab; **its agreement with Zhang to 1.3× is the score to beat**, and holding both out keeps that agreement as an out-of-sample fact |
| Hofmann 1996 (6 °C, organic solvent) | **neither** | an artefact-control study in ether/DCM/pentane; ingest only the MFT ≫ FFT **ordering** (§A.4), with Z2's magnitude warning |
| Zhu 2023 kafirin binding constants | **neither** | 53 % ethanol; not a Maillard matrix. Use to bound the reversible-partition share of an observed headspace loss, not as a scored row |

## D.5 — MODULE 3: ACRYLAMIDE / SAFETY

| dataset | **role** | one-line reasoning |
|---|---|---|
| Claeys 2005 T2 **control** system | **FIT** | the cleanest complete two-step set with SEs and a non-isothermal correction, and it self-validates |
| Claeys 2005 T2 **competitor** systems (Gln, Cys, Lys, Ala) | **FIT** | they identify the competition term; without them the control row alone is unidentifiable in a mixture |
| De Vleeschouwer I, glucose non-italic subset | **FIT** | the only bimolecular initiation constant; needed to give the trunk a real second-order form |
| **De Vleeschouwer I, FRUCTOSE and SUCROSE systems** | **★ HOLD-OUT** | a different sugar at the same a_w — the sugar-transfer test — and fructose's HPDs span zero, which the model should reproduce as *low confidence*, not as a fitted value |
| De Vleeschouwer II, **cysteine** system (k_E2, Ea_E2) | **FIT** | the tightest parameters in the corpus; they anchor the thiol/Michael-acceptor family that the sulfur module also needs |
| **De Vleeschouwer II, GLUTAMINE system** | **HOLD-OUT** | it carries the B5.5 sign-crossing (promotion grows with T in liquid, shrinks at a_w 0.92) — a shape no fitted term should be allowed to see first |
| Knol 2005 T1 (30 k + 6 Ea) | **FIT** | the largest self-consistent Arrhenius set for the Asn/Glc trunk |
| **Knol 2010 T2 (7 steps × 5 T)** | **★ HOLD-OUT** | a *third* lab on the same trunk; the only genuine cross-lab extrapolation test the acrylamide module can have |
| Knol 2009 real-food band (9.3 × 10³–2.6 × 10⁴ µg/kg dm) | **HOLD-OUT** | real food, not a model system, and the author himself refused to transfer (B8.5) — exactly what a hold-out is for |
| Knol 2009 **degradation** parameters | **neither** | SD ≥ estimate on every one; unidentifiable (C.6) |
| Quan 2020 | **neither** | no units, no orders (C.5) |
| existing `acrylamide_spi_extrusion_130C_ACSRef3` | **FIT** *(already)* | unchanged |

## D.6 — MODULES 4–8

### Module 4 — trunk / browning
| dataset | role | reasoning |
|---|---|---|
| Martins 2005 T2, all 10 steps + HPDs | **FIT** | the repo already fits this corpus; now with uncertainties |
| **Martins 2005 browning (step 9, ε 0.64)** | **★ HOLD-OUT** | the browning lane has **never** had a parameter; holding it out is the only way to learn whether the trunk predicts colour or merely accommodates it |
| Knol 2005 ε = 282 L mol⁻¹ cm⁻¹ (Glc/Asn) | **HOLD-OUT** | tests the amine-specificity of ε (§A.1) rather than fitting around it |
| Zheng 1994 Tables I/III/V (36 k + 8 Ea) | **FIT** | kinetic reference for cysteine thermolysis and β-elimination; **not a benchmark source** (no absolute concentrations) |
| Hemmler 2018 | **neither** | ESI intensity ≠ concentration; pH uncontrolled and drifting 2–3 units. Ordinal only |

### Module 5 — lipid oxidation
| dataset | role | reasoning |
|---|---|---|
| Frankel 1989, the 26-column distribution at zero additive | **FIT** | the only measured branch distribution; it must replace the shipped `hexanal 0.37` |
| **Frankel 1989, the α-tocopherol arms** | **★ HOLD-OUT** | tocopherol moves total and hexanal-share in **opposite** directions — a two-sided test a fitted split cannot fake |
| **Frankel 1989, the nonanal ABSENCE** | **HOLD-OUT (negative test)** | the model must predict ~zero nonanal from fed linoleate hydroperoxide; a fitted split would never be asked this otherwise |
| existing matrix-path hexanal hold-outs (Bi 2020 ×2) | **HOLD-OUT** *(already frozen)* | unchanged |

### Module 6 — protein binding
| dataset | role | reasoning |
|---|---|---|
| Damodaran 1981 soy K's (stated 100 kDa basis) | **FIT** | the strongest provenance in the batch; the basis is printed, not inferred |
| Andriot 2000 β-lactoglobulin K's | **FIT** | needed for the chain-length slope; label `recovered_by_arithmetic` |
| **Barallat-Pérez 2024 lupin + mucin constants** | **★ HOLD-OUT** | a third protein, a third method, and it carries the **method-dependence** finding (B6.7); fitting it would hide the very effect it proves |
| **Andriot 2000 sensory-intensity arm** | **HOLD-OUT** | the model predicts headspace; whether that maps to perceived intensity is precisely the open question (B6.4/B6.5) |
| Starkenmann 2008 saliva binding | **neither** | stranded, no basis, mechanism unresolved (§A.6) |

### Module 7 — thresholds / matrix correction
| dataset | role | reasoning |
|---|---|---|
| **Zhou 2023 SI Table S2 thresholds (15 compounds, water)** | **reference table, NOT a scored dataset** | thresholds are model *inputs*, not predictions; scoring them would be a category error |
| **Zhou 2023 SI Table S2 OAVs** | **HOLD-OUT (arithmetic check)** | already verified exact on 4 of 15 spot-checks; useful as a regression test on the OAV code path |
| Vega 1994 gelatin ladder (6 compounds × 4 T) | **FIT** *(as lookup-table entries)* | the cleanest matrix threshold set — no thermal step after dosing |
| **Brewer 1995 beef set** | **HOLD-OUT, and RECLASSIFIED `dose_added_pre_cook`** | it is not a threshold in the repo's sense (C.2); holding it out avoids fitting thermal loss into a perception term |
| Tian 2020 milk set | **neither, until the unit is settled** | a literal `?` in the units cell; factor-of-1000 basis risk |
| Guadagni 1963/72 aqueous references | **reference only** | held second-hand through Vega; every ratio in §A.7.5 divides by them |

### Module 8 — extrusion / matrix path
| dataset | role | reasoning |
|---|---|---|
| **Xin 2026 (Food Hydrocolloids) 9-formulation carbohydrate factorial** | **FIT** | the only dosed-precursor factorial in a true HME; it is what the extrusion lane exists to reproduce |
| **Xin 2026 xylose and ribose arms specifically** | **★ HOLD-OUT (carved out of the above)** | these are the two arms that **invert** the pentose ≫ hexose ordering (B1.2). Fitting them would let the model absorb the inversion as a parameter instead of explaining it. ⚠️ **This is the one place the split cuts *within* a paper — the cut is by ARM, and no arm appears twice** |
| Conti 2025 I + II (thiamine, 3 severities) | **FIT** | the only dosed-thiamine extrusion data; needed for the formation term |
| **Conti 2025 hexanal absolute pair (1.95 vs 4.70 µg/g)** | **HOLD-OUT** | the only absolute values in either Conti paper (the IS is hexanal-D12, so hexanal alone is properly quantified) — the strongest single number in the lane deserves to be out-of-sample |
| Guo 2020 retention series | **FIT** | pure retention; it identifies the loss term that Conti's net cannot separate |
| existing matrix-path hold-outs (Li 2026, Liu 2023) | **HOLD-OUT** *(already frozen)* | unchanged |
| **Xin 2026b six-formulation PROTEIN-substitution series** (Food Res. Int. 233, 119010) | **★ HOLD-OUT (new)** | it varies **protein composition at zero added carbohydrate** — precisely the axis the fitted Xin 2026 companion holds constant. Fit the carbohydrate axis, **predict** the protein axis. Same extruder, same settings, same group ⇒ a genuinely controlled cross |
| **Xin 2026b total free amino acids (HILIC-MS/MS, 1613.74 → 3347.27 µg/g)** | **FIT** | the only *validated absolute* measurement in either Xin paper; it constrains the amine precursor pool that every Maillard module divides by, and nothing else in the corpus measures it in an extrudate |
| **Xin 2026 vs Xin 2026b same-sample 10–23× discrepancy** | **⚠️ NOT A DATASET — a CALIBRATION FACT** (§E.2.2, §B9.1) | it must be applied as an `absolute_concentration: false` flag on other rows. **Scoring it would be a category error** |

### Module 9 — ★ NEW: the Z3 block (norfuraneol, thiamine route, furfural sink, browning Ea)

| dataset | **role** | one-line reasoning |
|---|---|---|
| **Bornhorst 2017b Ea(norfuraneol) = 104.9 / 121.1 / 122.3 kJ/mol + z** | **FIT** | the corpus's only norfuraneol Ea; it must go in as a *prior*, not a target, and it must carry the label "approach-to-plateau accumulation in an alkaline mashed-potato gel" |
| **Bornhorst 2017 90 °C matrix pair (egg white vs mashed potato, M-2∞ and k)** | **★ HOLD-OUT** | a *matrix*-transfer test at fixed T and fixed formulation — the model should predict a 3.9× k difference between two food gels it was not fitted to |
| Bornhorst 2017 structural zero (M-2 = 0 with no precursors, all three matrices) | **FIT** | a free, unambiguous zero; fitting it costs nothing and catches sign errors |
| **Bornhorst 2017b `2_R,2_L` a\* row at 100 °C, and the Ea 245.6 it drives** | **neither** | **D = 0.5 min < the 1.75 min come-up time.** Not a kinetic measurement (§F #24) |
| **Cerny 2007b full ternary MFT split (54 : 46)** | **FIT** | the branching ratio the thiamine lane has never had a number for |
| **★ Cerny 2007b the two single-route controls (no-cysteine > 99 : < 1; no-thiamine < 5 : > 95)** | **★ HOLD-OUT** | **the sharpest structural test in the whole split.** A model fitted on the ternary must *predict* both limiting cases; getting 54 : 46 right while missing either control means the routes are wrong and the ratio was fitted |
| **★ Cerny 2007 Table 5 concentration pair (85 : 15 at 1× vs 54 : 46 at 2×)** | **★ HOLD-OUT** | **the single highest-value hold-out row Z3 adds** (B10.1). It scores whether the model's branch fractions *respond to concentration at all* — the exact defect Z1 §E diagnosed. A model with fixed branch fractions fails it by construction, which is the point |
| Cerny 2007 Table 2 pH ladder (65 peak-area cells, 5 pH) | **HOLD-OUT (directional only)** | peak areas, never concentrations; but six sign-level switches and a third MFT-vs-pH shape (B10.10, B10.11) — a shape test the pH term should not see first |
| Cerny 2007 Table 4 isotope splits across pH | **FIT** | fractions, response-factor-immune; they identify which route supplies MFT at each pH, which is a prerequisite for any pH term at all |
| **Yaghmur 2005 water arm (k_obs ×4, Ea 46.50 ± 1.0)** | **neither — `audit_flag` only** | 40–70 °C, and a lump over ≳98.8 % non-FFT flux. **Not comparable to an Eyring barrier** (§C.14, §F #29). Record it; do not score against it |
| Yaghmur FFT-share bound (≲ 1.2 % of the furfural flux) | **FIT (as a ceiling)** | a one-sided constraint on the furfural → FFT branch, and the corpus has no other |
| **Ajandouz 2008 own-measurement Ea set (24 values)** | **FIT (as priors on Ea only)** | glucose-loss and amino-loss Ea transfer to pH 5–7 at ±10 % by the paper's own licence; **no rate transfers, and the browning Ea explicitly do not transfer** (§C.13) |
| **★ Ajandouz caramelisation partition (25–80 % of A₂₉₄; 7–55 % of A₄₂₀)** | **★ HOLD-OUT** | it sizes the amine-independent lane — a quantity the model computes but has never been scored on. Holding it out is the only way to learn whether `structural_zero` is right |
| Hofmann 1996b Sotolon anchors (13.5 / 764.7 / 273.1 µg, pH 3/5/7) | **FIT** | the corpus's first Sotolon numbers; nothing else constrains that node |
| **Hofmann 1996b oxidant series (argon 1.7 / air 9.8 / Cu²⁺ 20.1, 38.4 µg)** | **★ HOLD-OUT** | the only required-oxidant measurement in the corpus, and a 5.8× effect the model currently cannot express at all. ⚠️ ingest the second Cu²⁺ dose as "a higher level" (§F #19) |
| Hofmann 1996b AT time × temperature surface (Fig. 2, 30 points) | **HOLD-OUT** | `digitised_from_figure`, ±0.3 µg — but it is the corpus's **only** measured T × t surface with a **sign reversal**, and B10.4's Ea inequality is derived from it, so scoring it would be circular if it were also fitted |
| **Hofmann 1996b Tables 1 and 2** | **neither — DUPLICATES** | verbatim re-publications of Hofmann 1998 Tables 2, 8, 10 (§F #16). **Double-counting hazard** |
| Amrani-Hemaimi 1995 Table 2 isotope fractions (40 cells) | **FIT** | fractions, response-factor-immune, and the only pyrazine carbon-origin bookkeeping the corpus has |
| **Amrani-Hemaimi 1995 the alanine-vs-glycine ethyl-pyrazine ON/OFF switch** | **★ HOLD-OUT** | 20 / 19 % with alanine, **0 / 0 %** with glycine, and 100 % ¹³C-labelled ⇒ *"one single reaction route exists"*. An on/off amino-acid switch is unfakeable by a fitted continuous term |
| Amrani-Hemaimi Table 1 row 6+7 (co-elution) | **neither** | not a single-species value, and its labelling figure is an average of two compounds (§F #28) |
| **van Boekel 2005** | **neither** | one-page poster abstract, zero parameters (§C.16) |
| van Seeventer Z3 addendum §2/§3 tables | **neither** | arithmetic exercises run to prove an extrapolation cannot be made (§C.17, §F #30) |

## D.7 — COVERAGE CHECK

| module | ≥1 hold-out? | which |
|---|---|---|
| 1. sulfur formation | ✅ | Zhou pH 6/8 · Zhou Fig. 3 · Hofmann dry-180 · Whitfield 2001 · Cerny 2003 · Meynier (+ 5 frozen Hofmann rows) |
| 2. thiol consumption | ✅ | Hofmann 80 °C brew · van Seeventer 50 °C · Zhang 115 °C · Zhou 120 °C |
| 3. acrylamide / safety | ✅ | De Vleeschouwer fructose+sucrose · De Vleeschouwer glutamine · Knol 2010 · Knol 2009 real food (+ 12 frozen mp_holdouts) |
| 4. trunk / browning | ✅ | Martins step 9 + ε · Knol ε |
| 5. lipid oxidation | ✅ | Frankel tocopherol arms · Frankel nonanal absence (+ 2 frozen Bi 2020) |
| 6. protein binding | ✅ | Barallat-Pérez lupin + mucin · Andriot sensory arm |
| 7. thresholds / matrix | ✅ | Brewer beef (reclassified) · Zhou OAV arithmetic |
| 8. extrusion / matrix path | ✅ | Xin xylose + ribose arms · **Xin 2026b six-formulation protein series** · Conti hexanal pair (+ 2 frozen) |
| **9. Z3 block (norfuraneol / thiamine route / browning Ea)** | ✅ | **Cerny 2007 concentration pair · Cerny 2007b two single-route controls · Bornhorst matrix pair · Ajandouz caramelisation partition · Hofmann 1996b oxidant series + Fig. 2 surface · Amrani-Hemaimi on/off switch** |

**Disjointness check: no dataset appears in both columns.** The five places the split cuts inside
a single paper are stated explicitly and each cuts along a declared axis:

| # | paper | axis of the cut | fit side | hold-out side |
|---:|---|---|---|---|
| i | **Hofmann 1998** | **pH** | pH 5.0 rows | pH 3 / pH 7 / xylose (already frozen) |
| ii | **Zhou 2023** | **pH column** | pH 7 | pH 6 and pH 8 |
| iii | **Xin 2026** | **carbohydrate arm** | 7 sugars + control | **xylose, ribose** |
| iv | **Cerny 2007 / 2007b** | **system composition** | full ternary + isotope splits | **the two single-route controls, and the 1×/2× concentration pair** |
| v | **Bornhorst 2017 / 2017b** | **what is scored** | the Ea (as a prior) and the structural zero | **the 90 °C matrix pair (egg white vs mashed potato)** |

Convention (i) is the repo's existing one; (ii) mirrors it deliberately; (iii)–(v) each hold out
the *limiting case* or the *transfer axis* while fitting the interior, which is the only
arrangement under which a pass is evidence of structure rather than of fitting.

## D.8 — FOUR THINGS THE ORCHESTRATOR MUST DECIDE BEFORE THIS IS APPLIED

1. **Is cutting Zhou 2023 by pH column legitimate, given §C.11?** The pH labels are *initial*,
   and the pH-6 and pH-7 runs converge to within 0.2 units. If the model has no pH-trajectory
   state, the pH-6 hold-out may be **unfairly hard** — it would be asked to predict a 2.3×
   difference between two runs that spent most of their time at the same pH. **My
   recommendation: keep the cut, but score the pH-6 row as DIAGNOSTIC rather than gating until
   a pH-trajectory state exists.**
2. **Does holding out Knol 2010 orphan the trunk?** Knol 2010 T2 is also Z0's source for three
   "missing chemistry" lanes (isomerisation 61 ± 8, acetic 75 ± 10, formic 84 ± 14). If those
   are needed as *constants* they cannot also be a hold-out. **Recommendation: split by step —
   hold out the acrylamide steps, fit the organic-acid and isomerisation steps** (they belong to
   Module 4, not Module 3, so this does not violate rule 1).
3. **★ Can the alkaline block be scored at all?** Every Z3 constant in §A.1.1 and §A.3.6(i) is
   at **pH 8.0–9.7**, and all three source dossiers refuse rate transfer to pH 5–7 (§C.13). The
   split above puts the Ea in as *priors* and holds out the matrix pair and the caramelisation
   partition — but **if the model has no pH-dependent Ea, even the prior is a smuggled
   assumption.** **Recommendation: ingest the Ea with a mandatory `pH_of_measurement` field and a
   `rate_transfer: not_licensed` flag, and score the two hold-outs as DIAGNOSTIC until a
   pH-dependent barrier exists.** Note the sub-question the dossier raises and does not settle:
   *"Two of sixteen parameters get a referee; fourteen do not. Say that plainly rather than
   implying the fit has been validated."*
4. **Should the two `diagnostic_only` internal pilots be re-roled?**
   `pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026` and its soy twin are currently in
   the fit directory but flagged `diagnostic_only`. They are the only 100 °C protein-matrix
   sulfur data. **Recommendation: leave them exactly as they are** — this wave found no
   evidence bearing on them, and moving a diagnostic row into a scored role is precisely the
   kind of quiet promotion the audit trail exists to prevent.

---

# §E. WAVE K3 NEW EXTRACTIONS

Three dossiers were produced or completed this wave: `zhou2023_extraction.md` (new, the correct
paper at last), `Xin2026b_extraction.md` (new), and the **Wave K3 addendum** to
`Zhang2024_extraction.md` (completing an interrupted write). Their parameters are already folded
into §A–§D above. This section records what each one *changed* about the corpus.

## E.1 — Zhou et al. 2023, JAFC 71, 2472−2481 — **the corpus gap K1 declared is now CLOSED**

K1 §6 gap 1 read: *"`Zhou2023.pdf` is the wrong paper … **Zhou et al. 2023 JAFC — purified
deoxyxylulosyl-alanine Amadori compound + cysteine, pH 6/7/8, MFT rise-then-fall — is NOT in the
corpus.** It was the only pH-resolved MFT source and the only purified-Amadori source assigned to
me. **Re-fetch, and rename the existing file so no wave cites it as Zhou.**"

**Both halves are now done.** `data/articles/Zhou2023b.pdf` + `zhou2023b_supporting.pdf` are the
correct article and its SI. The mis-filed paper (Zhu et al., kafirin × Baijiu, *Food Chem.* 400
(2023) 133854) has had its dossier **preserved unchanged** as
`zhu2023_kafirin_binding_extraction.md` — it holds the corpus's only MFT–protein binding
constants and must not be lost. ⚠️ **`data/articles/Zhou2023.pdf` is still the wrong file for
DOI 10.1021/acs.jafc.2c08360 and should be renamed.**

**What it delivered:** a purified fed-Amadori initial condition (§A.3.3); a complete pH 6/7/8 ×
22-compound table; the first α-dicarbonyl-fed time × pH grid in the corpus; the author-stated
two-factor product law behind the MFT pH-7 maximum (§B2.4); a pH-invariant dimerisation branch
fraction that independently corroborates Zhang 2024 to 1.3×; a 15-compound threshold list
(§A.7.1); and a **new sign-crossing** that falsifies any family-level pyrazine pH term (§B2.3).

**And one new conflict the build waves must not average away: §B2.5.**

## E.2 — Xin et al. 2026b, *Food Res. Int.* 233 (2026) 119010 — **the jackpot was PARTIAL, and the by-catch is bigger**

Dossier: `Xin2026b_extraction.md`. 15 pp.; Figs. 1/3/4 re-read at 400–600 dpi.
**Multiprotein blends-regulated high-moisture extruded seafood analogs**, Dalian Polytechnic —
the same group, same extruder and same settings as the Xin 2026 companion.

**System:** six high-moisture extrudates at **70 SPI + 20 X + 10 WG** (X = SPI / **NPI** peanut /
**RPI** rice / **YP** yeast / **MPI** mung bean / **PPI** pea), **substitutive, constant total
protein, NO added carbohydrate**. Extruder: zones **60/80/110/135/140/140/140 °C**, **240 rpm**,
dry feed **2 kg/h**, **≈65 % moisture**, **70 °C cooling die**, **L/D 32:1**. Seven commercial
grilled seafoods are texture/colour comparators only — **no seafood volatile data**.

### E.2.1 The threshold jackpot — PARTIAL, but better evidenced than a table

**The threshold/OAV table is Table S5, online, NOT in the PDF.** However:

| threshold | value | basis | how obtained | provenance |
|---|---:|---|---|---|
| **1-octen-3-ol** | **1.5 µg/kg** | water | **printed verbatim in the body, p. 9** | **[C], NO CITATION GIVEN** |
| **hexanal** | **5.000 µg/kg** | water | **back-solved exactly from the paper's own printed OAV range** | [D] |
| **1-octen-3-ol** | **1.500 µg/kg** | water | " — **agrees with the printed 1.5** | [D] |
| **2-pentylfuran** | **5.800 µg/kg** | water | " | [D] |
| **2-heptanone** | **140.0 µg/kg** | water | " | [D] |

Three of the four back-solve **exact to 4–5 s.f. at BOTH ends of the printed range**.

**★ TWO INDEPENDENT CONFIRMATIONS FALL OUT:**
1. **The printed 1.5 µg/kg matches the back-solved 1.500** — an internal check that K2's
   OAV-inversion method is sound, using a value the paper actually prints.
2. **2-pentylfuran = 5.800 µg/kg now agrees across TWO papers from FOUR independent numbers**
   (K2 recovered it from the Food Hydrocolloids companion; this wave recovers it again here).
   **That threshold is now firmly established as the value this group uses.**

⚠️ **AND THE DEFECT REPRODUCES IN A SECOND JOURNAL.** §2.11 declares the OAV denominator an
**"air threshold value (μg/kg)"** while the values are unambiguously the classical **water**
thresholds — the identical defective sentence as the companion paper. **No source citation is
given for any threshold anywhere in either paper.** K2's laundering hazard #1 is therefore not a
one-off slip in one manuscript; it is this group's standing practice across two journals.
**It strengthens, rather than weakens, K2's §D.4.6 recommendation to fix the OAV denominator.**

### E.2.2 ★ THE BIGGEST FIND — a same-lab internal proof that HS-SPME µg/kg are not absolute

**This paper's `ENPI` formulation is bit-for-bit the companion paper's no-carbohydrate control**
— same protein blend, same extruder, same settings, same internal standard, same authors.
The two papers report:

| compound | Xin 2026 (Food Hydrocolloids), control | Xin 2026b (Food Res. Int.), ENPI | **ratio** |
|---|---:|---:|---:|
| **hexanal** | — | **74.86 µg/kg** | **≈10×** |
| **2-pentylfuran** | ≈97 000 µg/kg | **5 742.21 µg/kg** | **≈17×, RANGES NON-OVERLAPPING** |
| **nonanal** | — | **12.44 µg/kg** | **≈23×, RANGES NON-OVERLAPPING** |
| 2-heptanone | — | **2 778.1 µg/kg** | **1.03× — agrees to 3 %** |

**Likely cause, named in the dossier: this paper never states the SPME fibre or the extraction
conditions.** The companion uses DVB/CAR/PDMS.

> **★ THIS IS THE STRONGEST AVAILABLE EVIDENCE FOR A POSITION THE CORPUS HAS SO FAR ONLY
> ASSERTED: single-internal-standard HS-SPME µg/kg are NOT absolute concentrations.** Two papers
> by the same group, on the same physical sample, in the same year, differ by **10–23× on three
> compounds with non-overlapping ranges** — while agreeing to **3 %** on a fourth. That is not
> lab-to-lab variation and it is not a transcription error; it is fibre/partition dependence,
> and it is **internal to the literature**, requiring no assumption on the repo's part.
>
> **CONSEQUENCE: every semi-quantitative HS-SPME µg/kg row in the corpus — Conti 2025 I & II,
> both Xin papers, Zhang 2024, and Zhou 2023 — must carry an explicit `absolute_concentration:
> false` flag, and none of them may be scored against a SIDA anchor or against each other on
> absolute level. Within-paper ratios remain sound.** This retro-actively justifies every
> "RATIO-ONLY" verdict in §A and gives it a citable basis.

### E.2.3 Quantified content — usability by class

| class | numbers | verdict |
|---|---|---|
| **total free amino acids** | **1613.74 → 3347.27 µg/g**, by **validated HILIC-MS/MS** | **★ USE — the only trustworthy ABSOLUTE numbers in the paper** |
| secondary structure (β-sheet %) | **EMPI 32 % > EPPI 31 % > ENPI 27 % > ERPI/EYP 22 % > ESPI 20 %** | **USE** |
| fibrous degree | **ENPI 1.59 > EMPI 1.27 > EPPI 1.24** | **USE** |
| colour a\* | as printed | USE |
| aldehyde totals | **ERPI 2398.22 vs ENPI 193.95 µg/kg (12.4×)**; hexanal **74.86 (ENPI) → 1264.62 (EPPI)**; nonanal **12.44 → 195.04**; 2-pentylfuran **5742.21 → 17 033.09**; 2-heptanone spans only **1.10×** | **RATIO-ONLY / DIRECTIONAL — see E.2.2** |
| hardness | ≈17 000 (ENPI) → ≈35 500 g (EYP) **[fig]** | DIRECTIONAL |
| WAC / LF-NMR / E-nose | — | DIRECTIONAL |
| sensory | **ZERO sensory numbers are printed** | — |
| **every OAV** | — | **REFUSE — the denominator is mislabelled and uncited (E.2.1)** |
| **ΔE as a browning index** | — | **REFUSE — measured against a WHITE PLATE** |
| **"EPPI alcohols 1871.78"** | — | **REFUSE — a copy-paste of its own aldehyde total** |
| **"3 % / 5 %" water loss** | — | **REFUSE — the figure says 6.7–7.8 %** |
| **the β-sheet → fibrousness claim** | — | **REFUSE — the two rankings DISAGREE (EMPI is top on β-sheet, ENPI on fibrous degree), contradicting the paper's own mechanism AND its abstract** |

⚠️ **24 documented internal defects**, including **five malformed or misassigned references
(~8 % of the reference list)** and **dimethyl phthalate — a plasticiser — among the VIP
discriminants (VIP 1.38)**. Full register in `Xin2026b_extraction.md`.

### E.2.4 Its role in the split

| dataset | **role** | reasoning |
|---|---|---|
| **Xin 2026b six-formulation protein-substitution series** | **HOLD-OUT (new)** | it varies **protein composition at constant carbohydrate (zero)**, which is exactly the axis the fitted Xin 2026 companion holds constant. The two together are a clean composition-factorial cross: fit the carbohydrate axis, predict the protein axis |
| Xin 2026b **total free amino acids (HILIC-MS/MS)** | **FIT** | the only validated absolute measurement; it constrains the amine precursor pool that every Maillard module divides by, and nothing else in the corpus measures it in an extrudate |
| Xin 2026b **recovered thresholds** | **reference table** (`used_by_source`) | inputs, not predictions — same treatment as §A.7.1 |
| Xin 2026b **vs Xin 2026 same-sample discrepancy (E.2.2)** | **⚠️ NOT A DATASET — a CALIBRATION FACT** | it must be applied as a flag on other rows, never scored |

## E.3 — Zhang et al. 2024 — the interrupted write, completed

`Zhang2024_extraction.md` was complete through its ingestion table; the **Wave K3 addendum**
(§K3.1–K3.8) adds: the checklist against the K3 brief; the finding that **no concentration–time
profile exists in the retrieved PDF at all** and exactly which SI object holds each one; the
**failed unit recovery** for 0.0028 / 0.0031, recorded so no later wave re-attempts it; the
consolidated two-stage protocol; the consolidated Cys/GSH/**cystine** comparison with the
molar-sulfur correction that makes Cys-vs-GCys a genuinely sulfur-matched contrast; the explicit
formation-and-loss verdict; **§K3.7, the seven things it licenses that K1's 30 °C constant
cannot, and the one prohibited derivation**; and **§K3.8, the Zhou cross-validation to 1.3 %**.

## E.4 — `z3_index.md` — Wave Z3's missing consolidation, written this wave

Wave Z3 produced ten per-paper dossiers and was **interrupted before it wrote an index**, so its
parameters had never been consolidated the way Z0, Z1, Z2, K1 and K2 had. **`z3_index.md` now
exists** (≈250 parameter rows across 10 blocks, plus a directional list, verbatim gaps, an errata
register, per-paper verdicts and a [RETRIEVE] list). Everything in it is folded into §A.1.1,
§A.3.6, §B.10, §C.13–C.17, §F #15–30, §G and §D Module 9 above.

**Papers it covers, and what each turned out to be worth:**

| paper | what it is |
|---|---|
| **Ajandouz 2008** | **constant source (Ea only)** — 24 own-measurement Ea + the caramelisation partition; **alkaline, protein-bound amines, no rate transfer** |
| **Bornhorst 2017b** | **the corpus's first norfuraneol Ea**, verified three ways; apparent/lumped, alkaline, 3 points |
| **Bornhorst 2017** | benchmark source (norfuraneol formation at 90 °C in two food gels) + a structural zero |
| **Cerny 2007b** | **hard isotope constraints, no concentrations** — the only paper with BOTH two-component controls run |
| **Cerny 2007** | directional + the **concentration → route-mix coefficient** (§B10.1); peak areas only |
| **Hofmann & Schieberle 1996b** | the only measured **T × t surface with a sign reversal**, the only **required-oxidant** measurement, the first Sotolon anchors — ⚠️ **but Tables 1–2 are duplicates of Hofmann 1998** |
| **Yaghmur 2005** | one lumped furfural-sink Ea (46.50 ± 1.0 kJ/mol), 40–70 °C, over ≳98.8 % non-FFT flux |
| **Amrani-Hemaimi 1995** | directional + 40 hard isotope fractions; **zero absolute quantification** |
| **van Boekel 2005** | **nothing** — a one-page poster abstract. A dead-end retrieval, **not a repo defect** |
| **van Seeventer 2001 (Z3 addendum)** | a verification pass (Z2's dossier confirmed, nothing to correct) **plus a negative structural result** |

⚠️ **The Z3 consolidation opened no PDFs and no figures** — every value in `z3_index.md` is
transcribed from dossier text, and figure-derived values are marked **[D]** with the original
dossiers' read-off uncertainties. **One item in it is a consolidator observation rather than a
dossier claim and is labelled as such** (the Ajandouz "14.1" moisture-basis ambiguity, §F #23).

---

# §F. ERRATA AND LAUNDERING REGISTER — consolidated pointers

The full registers live in `k1_kinetic_parameters.md` §7 (7 entries) and
`k2_matrix_and_thresholds.md` "NAMED LAUNDERING HAZARDS" (14 entries + 4 unresolved internal
contradictions). The corpus-level headline entries:

| # | claim as printed | reality | source |
|---:|---|---|---|
| 1 | sulfur anchor **342 / 200 ppb** | **198 / 121 ppb**; 342/200 appear nowhere in Hofmann 1998 | repo (retired) |
| 2 | acrylamide **Ea = 129 kJ/mol** "from Knol" | absent from all three Knol papers; maxima 93, 102, none | `barrier_constants.py:274, :493`, `safety.py:790–795` |
| 3 | **Ea 52.1 / 72.9** "from Knol 2005" | true pair **94.4 ± 11 / 85.1 ± 14** | `safety_reference_payloads.json` entries[27] |
| 4 | `cysteine_thermolysis` **A = 1.0e14 s⁻¹** | refuted by its own source by **7–150×**; correct prefactors 9.8e11/1.9e12/2.4e13/1.0e12 at pH 3/5/7/9. **Ea 130.4 CONFIRMED exactly** | `arrhenius_params.yml` |
| 5 | `schiff_condensation` **A = 1.5e11** | over-predicts its cited source by **14.8×**; implied **≈1.0e10**. **Ea 97.0 CONFIRMED** | `arrhenius_params.yml` |
| 6 | Charles-Bernard T2 units "**mol⁻¹ s⁻¹**" | the values are **s⁻¹** | Charles-Bernard 2005 |
| 7 | De Vleeschouwer I k_Fref unit "10⁻³ **mm**⁻¹" | means **min⁻¹** | De Vleeschouwer 2009 I |
| 8 | Hwang 1995b: "asparagine had the highest oxazole contribution" | **asparagine produced NO oxazoles**; and the unit header "mg/g of glucose" should be **µg/g** | Hwang 1995b |
| 9 | Xin 2026: OAV "divided by its odor threshold **in air**" | the thresholds used are the classical **AQUEOUS** ones (recovered exactly) | Xin 2026 p. 3 |
| 10 | Brewer 1995: "**three orders of magnitude**" vs gelatin | true: **65×, 101×, 72×, 38.5×, 2.9×, 7.3×** — heptanal wrong by 35× | Brewer 1995 p. 593 |
| 11 | Barallat-Pérez: "mucin increased aroma binding **four to 12-fold**" | that is the ratio to **mucin alone**; adding mucin to the protein system gives **1.07–2.2×** | Barallat-Pérez abstract |
| 12 | Starkenmann: "**3 × 10⁶** above its odor threshold" | analytical LOD ÷ **sensory** threshold; same-method figure is **6 × 10⁴** | Starkenmann p. 9578 |
| 13 | Zhu 2023: DT's **ΔH = 14.22** kJ/mol | van't Hoff on its own K's gives **+18.3** | Zhu 2023 |
| 14 | **★ NEW** — Zhou 2023: bis(2-methyl-3-furyl) disulfide calibration `y = 0.0208x − 0.7126` | **zeroes at x = 34.3 µg/L**, so the pH-8 value of 50.07 is only 1.5× above a pseudo-LOD; "ND" for bis(2-furfuryl) disulfide may mean "< 4.1 µg/L" | Zhou 2023 SI T1 |
| **15** | **★ VIA Z3** — Cerny's prose: *"a higher precursor concentration (A) affords a **higher proportion (54 and 60 %, respectively) of ¹³C-labeled**"* | **The table's header reads `unlabeled : labeled`, so A gives 46 % and 40 % LABELLED, not 54 and 60.** The **same error repeats in BOTH papers**, and Cerny 2007b additionally **cites the wrong table** (says Table 5; the numbers are in Table 3). **INGEST THE TABLES (46 % / 40 %), NOT THE PROSE.** The *direction* the sentence asserts is correct | Cerny 2007 p. 1556 **and** Cerny 2007b p. 1314 |
| **16** | **★ VIA Z3** — Hofmann & Schieberle **1996b** Tables 1 and 2 | **VERBATIM RE-PUBLICATIONS of Hofmann & Schieberle 1998 Tables 2, 8 and 10.** *"This is a re-publication, not an independent replicate — the repo must not count it as a second measurement."* **DOUBLE-COUNTING HAZARD** | Hofmann 1996b |
| **17** | **★ VIA Z3** — the Hofmann DRY-heat protocol, described twice | **1996b: 5 min, 20.4 mg KH₂PO₄ + 0.3 mL water (= 0.50 M). 1998: 6 min, 300 µL of 2 mol/L phosphate.** *"Two published descriptions of one experiment differ by 20 % in time and 4× in buffer molarity while reporting the same number to four significant figures (1553.9 µg, 1.39 mol %)."* **A second, independent instance of Z0 finding #11 — on the DRY series** | Hofmann 1996b vs 1998 |
| 18 | **★ VIA Z3** — Hofmann 1996b names the silica phosphate salt two ways | Table 2 fn b says **KH₂PO₄**; Table 5 fn b says **K₂HPO₄** for the same protocol | Hofmann 1996b |
| 19 | **★ VIA Z3** — Hofmann 1996b Table 4's second Cu²⁺ dose | the text layer gives **"0.01" for BOTH rows** and the 300 dpi render is unresolvable. **The AT values 20.1 and 38.4 µg are certain; the 0.1 µmol dose IS NOT CONFIRMED.** Ingest as "0.01 and a higher Cu²⁺ level" | Hofmann 1996b T4 |
| 20 | **★ VIA Z3** — Ajandouz Table 2's NFDM (46, 39) and pasta (31) kJ/mol rows | **recomputed by Ajandouz from other people's published data** (footnote a). *"Citing these as 'Ben-Gara & Zimmerman' would launder an Ajandouz recomputation as a primary measurement."* | Ajandouz 2008 T2 |
| 21 | **★ VIA Z3** — Ajandouz's UV A₂₉₄ Ea values (152/128/129; 123/107/101) | exist **only in the running text of §3.4 and are NOT TABULATED.** There is no table anchor to point at | Ajandouz 2008 §3.4 |
| 22 | **★ VIA Z3** — Ajandouz's a_w–Ea direction | **cited and then contradicted in the same paragraph** (Malec 2002 says Ea falls with a_w; Thompson 1976 and Jokinen 1976 saw no effect). *"Quoting either half alone would launder an open question as a settled direction."* | Ajandouz 2008 |
| 23 | **★ VIA Z3 (consolidator observation, not a dossier claim)** — Ajandouz Table 2's Ribose + BSA row prints **"14.1" WITHOUT parentheses** in a column headed *"a_w or (% H₂O)"* | under the table's own convention that reads as **a_w = 14.1 — physically impossible.** Presumably 14.1 % H₂O with the parentheses lost. **This is the most repo-relevant row in the table (ribose, 85–145 °C): treat its moisture basis as UNITS NOT STATED and resolve it from Carpenter 1962 before use** | Ajandouz 2008 T2 |
| 24 | **★ VIA Z3** — Bornhorst 2017b Table 3, `2_R,2_L` a\* at 100 °C | **D = 0.5 min against a come-up time of 1.75 min** — *"the reaction is essentially complete during come-up"* — and it drives the paper's largest Ea (**245.6 ± 72.1**). **Ingesting that Ea would launder a come-up artefact.** The authors publish the row without comment | Bornhorst 2017b T3 |
| 25 | **★ VIA Z3** — Bornhorst 2017b never states the z-value reference temperature | **T_ref = 363.15 K recovered by the Z3 analyst from the arithmetic.** Any downstream use of these z-values without stating T_ref is unanchored | Bornhorst 2017b |
| 26 | **★ VIA Z3** — Bornhorst 2017 §3.2 pairs the *low* k with the *low* D | cosmetic ordering only; the range is right (`D = 2.303/k` → 451.6 and 125.8 min), and 2017b's Table 1 confirms | Bornhorst 2017 §3.2 |
| 27 | **★ VIA Z3** — van Boekel 2005 misnames its own Amadori compound, **twice** | the system is glucose + **glycine** (⇒ fructosyl-**glycine**) but the abstract twice says *"the Amadori product, **fructosyl-lysine**"*. *"Either the abstract is sloppy or the experiment is not the one its title claims."* ⚠️ **No repo file cites this DOI — a dead-end retrieval, not a defect** | van Boekel 2005 |
| 28 | **★ VIA Z3** — Amrani-Hemaimi Table 1 row 6+7 | an **unresolved co-elution the authors declined to resolve** (*"we did not try to verify this point"*). **The 35 % and 42 % glycine-column values are NOT single-species**, and Table 2's 50/40/80/100 labelling figures for that row are an **AVERAGE OF TWO COMPOUNDS** | Amrani-Hemaimi 1995 |
| 29 | **★ VIA Z3** — Yaghmur's Ea filed as a refutation of a repo barrier | **category error.** *"(i) the measured quantity is a lump over ≳98.8 % non-FFT flux; (ii) an Arrhenius Ea and an Eyring ΔG‡ are different quantities and the repo's barriers are the latter; (iii) the measurement is at 40–70 °C."* **File as an `audit_flag`, value unchanged** | Yaghmur 2005 |
| 30 | **★ VIA Z3** — the van Seeventer addendum's §2 and §3 tables | **ARITHMETIC EXERCISES, not parameters.** §2 is *"run precisely to establish that it cannot be run"*; §3 is *"a sensitivity, not a determination."* **Extracting any single row as a constant would invert the dossier's purpose** | `vanseeventer2001_z3_addendum.md` |

**Two PDF-extraction hazards that must be re-checked on every future ingest:**
`schieberle2005.pdf`'s text layer **silently drops an entire row of Table 2** (pH 9.0, 38.4 µg,
3.1 mol %) and garbles Table 3 (`*9` for `>99`); `martins2005.pdf`'s layer **strips the minus
sign from every exponent** in Table 2, turning `1.6 × 10⁻⁵` into `1.6 × 10⁵`.
**Both were only caught by re-rendering at 300–400 dpi.**

---

# §G. RETRIEVAL QUEUE — ranked

| # | paper | why | requested by |
|---:|---|---|---|
| **1** | **Zhang et al. 2024 Supplementary Information** (`10.1016/j.foodres.2024.114149`, Appendix A) | **Fig. S2 is the entire temperature axis of the sulfur-consumption question**; Fig. S1 is every time series; Tables S6/S7 are the kinetics whose units cannot otherwise be recovered | Z3 / K3 |
| 2 | **Xin et al. 2026 SI** (`10.1016/j.foodhyd.2026.113124`) Tables S2/S3 | ~12 quoted compounds → **405 quantified cells**, plus the threshold list | K2 |
| **2b** | **Xin et al. 2026b SI, Table S5** (`10.1016/j.foodres.2026.119010`) | the **20 remaining thresholds**, and — the real question — **whether any of them carries a citation**. Four are already recovered by arithmetic (§A.7.2); S5 would settle the provenance of the whole list. Then S4 and S1 | **K3, new** |
| 3 | **Tian et al. 2019** `10.3168/jds.2019-16796` | the only way to settle the `?/kg` unit; 7 milk thresholds are unusable without it | K2 |
| 4 | **Zhou et al., JAFC 2022, 70, 15202−15212** — "Competitive formation of 2,3-butanedione and pyrazines through intervention of added cysteine during thermal processing of alanine-xylose Amadori compound" | the **dedicated** study of the Cys/α-dicarbonyl competition that Zhou 2023 measures only in passing (B4.1) | **K3, new** |
| 5 | **Liu et al., JAFC 2022, 70, 11643−11651** — "Formation priority of pyrazines and 2-acetylthiazole dependent on the added cysteine and fragments of deoxyosones…" | the **MGO-vs-GO selectivity** study — the missing half of the B2.3 paradox | **K3, new** |
| **5b** | **★ Brands & van Boekel 2002** — lactose + casein, **90–130 °C, 15-step multiresponse INCLUDING MELANOIDINS, Ea 71–159 kJ/mol** | *"**THIS is the paper the brief was actually after**, and Wave S3 already names it as its independent comparator"* — serves two purposes at once. ⚠️ **`brands2001.pdf` IS in the corpus; brands2002 is NOT.** No DOI given; also cited as "Brands (2002), thesis Chapter 4" | **Z3 — named its highest-value retrieval** |
| **5c** | **★ Carpenter, Morgan, Lea & Parr 1962** — **ribose + BSA, 85–145 °C**, FDNB, Ea 130 kJ/mol | *"Ribose, and a temperature window that actually contains the repo's processing conditions. **No other row in the 38-paper corpus has both.**"* Also resolves the "14.1" moisture-basis ambiguity (§F #23). No DOI given (second-hand via Ajandouz Table 2) | **Z3** |
| 6 | **Buttery, Ling & Juliano 1982**, *Chem. Ind.* 958–959 | the **only** route to a 2-AP odour threshold; 2-AP's whole claim is an OAV argument | Z2 |
| **6b** | **Lau et al. 2003** (whey protein gel, **116–131 °C**, M-2 Ea 64.0–122.3) and **Pandit et al. 2006** (mashed potato, **121 °C**, M-2 Ea 81.4–96.1, **zero added lysine**) | *"between them they would extend the norfuraneol Ea record into the repo's actual processing window"* — the Bornhorst Ea are 80–100 °C and alkaline. Lau is also the source of the `C = C∞ − (C∞−C₀)e^{−kt}` equation. No DOIs given | **Z3** |
| **6c** | **Zeiler** (ref 26 of Cerny 2007, a thesis) — thiamine vs ribose/cysteine MFT, **45 vs 3.4 µg/L** at 120 mmol/L pH 5.7, rising to **1500×** at meat concentrations | the only absolute thiamine-MFT concentration pair in the corpus, and the largest concentration→route-mix claim — currently held only second-hand, with the "meat concentration" never numerically specified. No DOI (thesis) | **Z3, medium** |
| **6d** | **Güntert et al. 1990** (thiamine pH series, 130 °C / 6 h) · **Dwivedi & Arnold** ([³⁵S]thiamine, pH 3.5–8.0) · **Grosch et al. 1993** · **Hincelin et al. 1992** | the four counterparties Cerny 2007/2007b names but the corpus does not hold — including the two whose contradiction Cerny reconciles (>100× concentration difference) and the norfuraneol + H₂S → MFT hypothesis Cerny refutes. No DOIs given | **Z3, medium** |
| **6e** | **Zhang et al. 2014** (egg white, 75–100 °C, **zero-order** M-2) · **Cerny & Grosch 1994**, *Z. Lebensm. Unters. Forsch.* **198**:210–214 | the reaction-order counterparty to Bornhorst's first-order fit; and the replicate of the alanine ethyl-pyrazine on/off switch (B10.14) | **Z3, low** |
| 7 | **Damodaran & Kinsella 1981b**, *JAFC* 29, 1253–1257 | closes the **7S/11S confound** — the most plausible non-error explanation of the 3–8× soy spread | Z0 |
| 8 | **Cerny & Guntz-Dubini 2008**, `10.1021/jf801762c` | primary source for 5-hydroxy-3-mercapto-2-pentanone, the thiamine → MFT intermediate | Z3 |
| 9 | **Harrison & Hills 1997**, *JAFC* 45, 1883–1890 | primary source of Andriot's eqs 1–6; would upgrade the β-lg basis from `recovered_by_arithmetic` to `stated_by_source` | K2 |
| 10 | **Quan 2020 version of record** + its Supplementary Figure 1 | would unlock 88 fitted constants currently unusable | K1 |
| 11 | **Laroque et al. 2008**, *Food Chem* 111:1032 | cited as the source of the pentose > hexose **rate** claim; an actual kinetic paper on sugar reactivity | K1 |
| 12 | **Milani, Menis-Henrique & Conti 2022** | the only plausible source of a **thiamine dose–response**, which neither Conti paper provides | K2 |
| 13 | **Guadagni, Buttery & Turnbaugh 1972**, *J. Sci. Food Agric.* 23, 1435–1444 | primary source of **four of the six aqueous thresholds** every ratio in §A.7.5 divides by | K2 |
| 14 | **Knol et al. 2008**, *Mol. Nutr. Food Res.* 52:313–321 | the fourth genotype underpinning Knol 2009 Eq. 6's R² = 0.98 sugar regression | Z2 |
| 15 | **Hwang, Hartman & Ho 1995**, *JAFC* 43, 179−184 (the pyrazine twin) | the pyrazine data currently reach the corpus only as an aggregate in a figure | K1 |
