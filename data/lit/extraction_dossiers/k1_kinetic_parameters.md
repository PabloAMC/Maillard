# Wave K1 — CONSOLIDATED KINETIC PARAMETERS 2026-08-28

Every rate constant, activation energy, half-life, binding constant and rate-attenuation
factor found in the 10 papers assigned to Wave K1. Per-paper dossiers:
`charles-bernard2005_extraction.md`, `hofmann2001_extraction.md`, `claeys2005_extraction.md`,
`devleeschouwer2009b_extraction.md` (Part I), `devleeschouwer2009_extraction.md` (Part II),
`martin2001_extraction.md`, `hwang1995b_extraction.md`, `hemmler2018_extraction.md`,
`quan2020_extraction.md`, `Zhou2023_extraction.md`.

**Provenance codes used throughout:**
- **[M]** measured/reported directly by the authors
- **[F]** FITTED by the authors (regression parameter)
- **[K1]** derived by me from the paper's own data — never reported by the authors
- **±** meanings differ by paper and are stated per block (SE vs 95% HPD vs SD vs none)

---

## 0. THE ONE-PARAGRAPH ANSWER

Of the ten papers, **three carry parameters that a mass-action ODE core can use today**
(Claeys 2005; De Vleeschouwer 2009 I & II), **two carry the thiol-consumption constants
the repo's biggest structural defect needs** (Charles-Bernard 2005; Hofmann & Schieberle
2002), **one carries equilibrium binding constants for MFT** (the mis-filed `Zhou2023.pdf`),
**one has 88 rate constants that are unusable as printed** (Quan 2020 — no units, no rate
laws, mislabelled headers), and **three contain no kinetic parameters at all**
(Martin & Ames 2001; Hwang 1995b; Hemmler 2018 — all ordinal/structural only).

**Two hard negatives worth stating up front:**
1. **No activation energy for thiol consumption exists anywhere in this basket.** Both
   thiol papers are single-temperature. Any T-extrapolation of the FFT/MFT sink is the
   repo's own assumption and must be labelled as such.
2. **Two requested papers are effectively missing**: `Zhou2023.pdf` is a different paper
   entirely (Zhu et al., kafirin/Baijiu), so the Amadori/cysteine pH 6/7/8 MFT study is
   **not in the corpus**; and Hwang's pyrazine twin (JAFC 43:179−184) is not there either.

---

## 1. THIOL CONSUMPTION — the parameters that fix the repo's missing sink

### 1a. Charles-Bernard, Roberts & Kraehenbuehl 2005, Table 2 (p. 4430) — **[F]**
**Shared conditions for this whole block:** reconstituted coffee brew, **1% total solids
in 0.01 M acetate buffer, pH 5.2, 25 °C, aerobic**, thiols spiked at 3.2–9.7 µmol/L,
headspace SPME-GC-MS vs a coffee-free blank, duplicates, **no CIs reported**.
**Rate law stated by the authors: pseudo-FIRST-ORDER in thiol**, ln C vs t[s], slope = −k.
**UNIT CORRECTION: the table header prints "[mol⁻¹ s⁻¹]" — it is wrong. The values are
s⁻¹.** Proven in the dossier against five of the paper's own Figure-2 half-lives.

| parameter | value | units | reaction/step | usability verdict |
|---|---|---|---|---|
| k(PhSH) | **> 7.70 × 10⁻⁴** | s⁻¹ | thiophenol + coffee matrix → adduct | directly usable, lower bound only |
| **k(FFT)** | **> 7.70 × 10⁻⁴** | s⁻¹ | **2-furfurylthiol + coffee matrix → adduct** | **directly usable, lower bound only** |
| k(BnSH) | > 7.70 × 10⁻⁴ | s⁻¹ | benzylthiol + matrix | directly usable, lower bound |
| k(PentSH) | 1.83 × 10⁻⁴ | s⁻¹ | pentanethiol + matrix | directly usable |
| k(BuSH) | 1.42 × 10⁻⁴ | s⁻¹ | butanethiol + matrix | directly usable |
| k(PropSH) | 1.32 × 10⁻⁴ | s⁻¹ | propanethiol + matrix | directly usable |
| k(EtSH) | 1.19 × 10⁻⁴ | s⁻¹ | ethanethiol + matrix | directly usable |
| k(2BT) | 8.11 × 10⁻⁵ | s⁻¹ | 2-butanethiol + matrix | directly usable |
| k(MMBF) | 2.12 × 10⁻⁵ | s⁻¹ | 3-mercapto-3-methylbutyl formate + matrix | directly usable |
| k(2M2P) | 1.02 × 10⁻⁶ | s⁻¹ | 2-methyl-2-propanethiol + matrix | directly usable |

**Derived half-lives [K1]** (t½ = ln2/k, same conditions):
FFT/BnSH/PhSH **< 15.0 min**; PentSH 63 min; BuSH 81 min; PropSH 87 min; EtSH 97 min;
2BT 2.37 h; MMBF 9.1 h; 2M2P 189 h. Span 2M2P→FFT **> 755×**.

**Second-order recast [K1]** — combining Table 2 with the p. 4428 hydroxylamine titration
of **8–10 mmol electrophilic sites per g dry coffee** (= 0.08–0.10 M at 1% t.s.):
| parameter | value | units | step | verdict |
|---|---|---|---|---|
| k₂(EtSH + matrix electrophile) | ≈ **1.3 × 10⁻³** | M⁻¹ s⁻¹ | bimolecular thiol + electrophile | **[K1] derived, order-of-magnitude only** — makes the sink scale with matrix loading |
| k₂(FFT + matrix electrophile) | ≳ **8.6 × 10⁻³** | M⁻¹ s⁻¹ | bimolecular | **[K1] derived, order-of-magnitude only** |
| [E] matrix electrophile capacity | **8–10 mmol/g dry coffee solids** | mmol g⁻¹ | titratable sink pool | **[M] directly usable** — the only published site density |

### 1b. Charles-Bernard 2005, Table 3 — rate attenuation 1/k_rel, **[M]**, dimensionless
Same conditions; additive amounts in **mmol per g dry coffee solids**.
| additive (amount mmol/g) | EtSH | PropSH | BuSH | PentSH | **FFT** |
|---|---|---|---|---|---|
| aerobic (ref) | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| anaerobic | 17.9 | 19.8 | 15.9 | 16.4 | **45.1** |
| hydroxylamine (10) | 10.5 | 11.4 | 10.8 | 13.4 | **4.0** |
| anaerobic + hydroxylamine (21), vs aerobic | 64.8 | 82.5 | 101.3 | 236.1 | 92.6 |
| anaerobic + hydroxylamine (21), vs anaerobic | 3.6 | 4.2 | 6.4 | 14.4 | 4.1 |
| Na₂SO₃ (1) | 51.5 | 61.7 | 18.4 | 62.7 | 67.9 |
| ascorbic acid (1) | 9.0 | 7.8 | 4.4 | 3.7 | **23.6** |
| caffeic acid (0.28) | **0.3** | **0.3** | **0.3** | **0.3** | 1.0 |
| DTPA (0.5) | 0.5 | 0.7 | 0.7 | 1.0 | 1.0 |
**Verdict: structural, directly usable as branch-fraction constraints.** Reading these as
channel weights: aliphatic thiols have **one** composite O₂-dependent electrophilic-addition
channel; **FFT retains ~25% of its rate with all electrophiles blocked**, so benzylic
thiols need a **second, non-electrophilic, O₂/radical channel**. Note caffeic acid at 0.3
is a **pro-oxidant** (accelerates 3×) and DTPA ≈ 1 means **Fenton chemistry is not the
main route**.

### 1c. Hofmann & Schieberle 2002 (file `hofmann2001.pdf`) — **NO parameters published**
*JAFC* **2002**, 50, 319−326 (note: filename year is wrong). The paper reports no rate
constant, Ea or half-life. Everything below is **[K1]** — my first-order regressions on
figures digitized at 400–450 dpi. SIDA quantification ([²H₂]-FFT, [²H₆]-MMBF etc.).

| parameter | value | units | reaction/step | conditions | verdict |
|---|---|---|---|---|---|
| **k(FFT + melanoidins)** | **5.6 × 10⁻² min⁻¹ = 9.4 × 10⁻⁴ s⁻¹** (t½ **12.4 min**) | min⁻¹ / s⁻¹ | FFT + coffee melanoidin → thioether | **30 °C, pH 6.0, 0.1 M phosphate, 12.5 g/L melanoidin (MW>3000), FFT 438 µM** | **[K1] fit, Fig. 6; first order confirmed (k constant ±18% over 80× extent)** |
| k(FFT + lys/glycolaldehyde CROSSPY model) | 9.8 × 10⁻⁴ s⁻¹ | s⁻¹ | FFT + pyrazinium dication | 30 °C, pH 6.0, 30 min, dry-heated 5 min/230 °C model | **[K1] from Table 2 single point** |
| k(FFT, real brew) | 0.023 min⁻¹ = 3.8 × 10⁻⁴ s⁻¹ (t½ ~30 min) | min⁻¹ / s⁻¹ | FFT loss in coffee brew | **80 °C**, thermos, 50 g powder/L | **[K1] fit, Fig. 1** |
| k(MMBF, real brew) | 0.0084 min⁻¹ = 1.4 × 10⁻⁴ s⁻¹ (t½ ~83 min) | min⁻¹ / s⁻¹ | MMBF loss in coffee brew | 80 °C, thermos | **[K1] fit, Fig. 1 — clean first order (k varies only 0.0071–0.0093 over 7× extent)** |
| k(FFT + albumin/glycolaldehyde) | 6.5 × 10⁻⁴ | s⁻¹ | — | 30 °C, pH 6.0, 30 min | [K1] single point |
| k(FFT + albumin/glucose) | 3.0 × 10⁻⁴ | s⁻¹ | — | 30 °C, pH 6.0, 30 min | [K1] single point |
| k(FFT + heated chlorogenic acid) | 8.4 × 10⁻⁵ | s⁻¹ | — | 30 °C, pH 6.0, 30 min | [K1] — essentially inert |
| k(FFT + chlorogenic acid) | 4.6 × 10⁻⁵ | s⁻¹ | — | 30 °C, pH 6.0, 30 min | [K1] — essentially inert |
| FFT bound per g melanoidin | ≥ **0.028 mmol/g** | mmol g⁻¹ | covalent capacity | 30 °C, 90 min | **[K1] LOWER BOUND — FFT was limiting** |
| disulfide branch | **< 1.5%** of thiol flux | — | 2 RSH → RSSR | 30 °C, pH 6.0 | **[M] — kills the oxidation-to-disulfide sink implementation** |

**CROSS-VALIDATION — the strongest result in the whole basket:**
| system | T | matrix loading | k(FFT) s⁻¹ | basis |
|---|---|---|---|---|
| Hofmann Fig. 6 | 30 °C | 12.5 g/L melanoidin | **9.4 × 10⁻⁴** | SIDA, absolute |
| Hofmann Table 2 | 30 °C | lys/glycolaldehyde model | **9.8 × 10⁻⁴** | SIDA, absolute |
| Charles-Bernard Table 2 | 25 °C | 10 g/L coffee solids | **> 7.7 × 10⁻⁴** | headspace ratio |
**Two labs, two methods, two matrices, agreement within 25%. Adopt k(FFT sink) ≈ 9 × 10⁻⁴
s⁻¹ at 30 °C in a ~10 g/L coffee-solids matrix, pH 5–6.**

### 1d. The only Ea in reach for a thiol — and I recommend REFUSING it
| parameter | value | derivation | verdict |
|---|---|---|---|
| Ea(MMBF loss) | ≈ **30 kJ/mol** | two points: Charles-Bernard 2.12 × 10⁻⁵ s⁻¹ (25 °C) vs Hofmann 1.4 × 10⁻⁴ s⁻¹ (80 °C) | **[K1] — DO NOT INGEST.** Two labs, two matrices, two analytical bases, two points. The same treatment applied to FFT gives a **negative** Ea (the 80 °C brew is *slower* than the 30 °C models, because the real brew's electrophile pool was partly consumed during extraction). **The literature reviewed supports NO activation energy for thiol consumption.** |

---

## 2. ACRYLAMIDE FORMATION / ELIMINATION — the well-parameterised block

### 2a. Claeys, De Vleeschouwer & Hendrickx 2005, Table 2 (p. 1528) — **[F]**, ± = **SE**
**Shared conditions:** **0.01 M L-asparagine + 0.01 M D-glucose + 0.01 M competitor**, in
**0.05 M citrate buffer pH 6**, closed inox tubes, **140/160/180/200 °C**, non-isothermal
correction by Euler integration of a 2 s thermocouple log. **T_ref = 160 °C.**
**Rate law, stated: two consecutive FIRST-ORDER reactions**, dC_R/dt = −k_F·C_R;
dC_AA/dt = k_F·C_R − k_E·C_AA. Arrhenius in reference form
k = k_ref·exp[(Ea/R)(1/T_ref − 1/T)], R = 8.314. n = 35 points per system.

| system | k_Fref (×10⁻³ min⁻¹) | k_Eref (×10⁻³ min⁻¹) | Ea_F (kJ/mol) | Ea_E (kJ/mol) | pseudo-R² | Pr<W |
|---|---|---|---|---|---|---|
| **control** | **0.451 ± 0.023** ᵃ | **111.1 ± 8.9** ᵃ | **168.25 ± 3.80** ᵃ | **167.21 ± 4.30** ᵃ | 0.975 | 0.934 |
| glutamine | **1.640 ± 0.416** ᵇ | 274.1 ± 81.6 ᵃᵇ | 166.8 ± 14.4 ᵃᵇ | **103.9 ± 17.5** ᵇ | 0.912 | 0.620 |
| cysteine | 0.501 ± 0.116 ᵃ | 268.7 ± 81.6 ᵃᵇ | **206.3 ± 13.5** ᵇ | 180.0 ± 17.0 ᵃᶜ | **0.765** | 0.599 |
| lysine | 0.587 ± 0.074 ᵃ | **280.2 ± 43.4** ᵇ | 179.3 ± 7.7 ᵃᵇ | **140.0 ± 9.0** ᵇᶜ | 0.950 | 0.572 |
| alanine | 0.465 ± 0.034 ᵃ | 103.1 ± 12.0 ᵃ | 173.3 ± 5.4 ᵃᵇ | 169.7 ± 6.1 ᵃᶜ | 0.956 | 0.096 |
Letters = significantly different on 95% asymptotic CIs.

**Derived [K1]:** k_E/k_F ratio — glutamine 167 < alanine 222 < control 246 < lysine 477 <
**cysteine 536**. AA half-life from k_E at 160 °C: control 6.24 min, glutamine 2.53,
cysteine 2.58, lysine 2.47, alanine 6.72 min.

**Verdict: directly usable in a mass-action ODE core — no unit conversion needed beyond
the ×10⁻³ prefix.** TWO CONDITIONS: (i) **k_F and k_E are 0.94–0.98 correlated** (authors
state this explicitly) — ingest the **ratio** as the well-determined quantity or carry the
correlation, otherwise the repo will overstate its information; (ii) the constants are
**lumped and matrix-specific** to 10 mM/10 mM/10 mM, pH 6, closed aqueous — they will not
transfer to a different competitor loading or to low a_w.

**Internal validation [K1] (both pass):** quasi-plateau (k_F/k_E)·C_R0 = 4.06 × 10⁻⁵ M =
**2886 ppb** vs Fig. 2a's ~2800 ppb plateau; t_max = ln(k_E/k_F)/(k_E−k_F) = **49.8 min**
vs Fig. 2a's ~40–50 min peak at 160 °C.

### 2b. De Vleeschouwer et al. 2009 **Part I** (`devleeschouwer2009b.pdf`), Table 3 (p. 124) — **[F]**, ± = **95% HPD interval**
**Shared conditions:** equimolar asparagine–sugar, **freeze-dried powder equilibrated to
a_w 0.92 at 4 °C** (moisture 12.6–17.8%), **~3 M reactant concentrations**, hermetically
sealed inox tubes, **120/140/160/180/200 °C**, 4 s thermocouple log, **multiresponse fit
by the determinant criterion** (Athena Visual Studio v11.0). **T_ref = 160 °C.**
**Rate laws (eqs 2–10): k_INT·[Asn]·[Sugar] SECOND ORDER; k_C·[Sugar] first; k_F·[Int1]
first; k_M·[Int1] first; k_B·[Int2]·[Asn] SECOND ORDER; k_E·[AA] first; k_Asp·[Asn] first;
k_X·[Asp] first.**
**⚠️ Parameters in *italic* below are marked by the authors as having "NO PHYSICAL
MEANING" — DO NOT INGEST THEM.** **⚠️ k_Fref's printed unit "(10⁻³ mm⁻¹)" is a typo for
min⁻¹.**

| parameter | Glucose (Scheme 4) | Fructose (Scheme 4) | Sucrose (Scheme 3) | units | step | verdict |
|---|---|---|---|---|---|---|
| k_Fref | **3.57 ± 1.38** | 7.40 ± 9.48 | 3.56 ± 0.86 | ×10⁻³ min⁻¹ | Int1 → acrylamide | usable (glucose); fructose HPD spans 0 |
| k_Eref | **0.10 ± 0.04** | 0.09 ± 0.02 | **0.74 ± 0.17** | min⁻¹ | acrylamide → DP | usable |
| k_INTgref | **1.70 ± 1.05** | – | 1.28 ± 1.09 | **M⁻¹ min⁻¹** | Asn + Glc → Int1 | usable |
| k_INTfref | – | 0.22 ± 0.38 | 1.73 ± 0.80 | M⁻¹ min⁻¹ | Asn + Frc → Int1 | HPD spans 0 for fructose |
| *k_Mref* | *1.23 ± 0.49* | *0.58 ± 0.16* | *0.04 ± 0.01* | min⁻¹ | Int1 → Int2 | **REFUSE — no physical meaning** |
| *k_Bref* | *3.90 ± 3.68* | *0.63 ± 30.64* | *4.49* ᶜ | M⁻¹ min⁻¹ | browning | **REFUSE** (fructose HPD 49× the estimate) |
| *k_Cgref* | *Indeterminate* | – | *0.48 ± 0.23* | ×10⁻³ min⁻¹ | caramelisation | **REFUSE** |
| *k_Cfref* | – | *1.03 ± 0.86* | *716.18 ± 122.10* | ×10⁻³ min⁻¹ | caramelisation | **REFUSE** |
| k_Aspref | **26.43 ± 5.76** | 13.62 ± 4.08 | 7.29 ± 2.60 | ×10⁻³ min⁻¹ | Asn → Asp (deamidation) | usable |
| *k_Xref* | *Indeterminate* | *2.12 ± 3.25* | *1.28 ± 5.11* | ×10⁻³ min⁻¹ | Asp consumption | **REFUSE** |
| k_HYref | – | – | 0.47 ± 0.09 | ×10⁻³ min⁻¹ | sucrose hydrolysis | usable |
| k_Iref | – | – | 0.50 ± 0.20 | ×10⁻³ min⁻¹ | Glc → Frc isomerisation | usable |
| **Ea_F** | **159.2 ± 29.5** | 122.0 ± 66.5 | 96.2 ± 22.3 | kJ/mol | acrylamide formation | usable |
| **Ea_E** | **113.2 ± 32.3** | 95.4 ± 18.3 | 108.9 ± 21.8 | kJ/mol | acrylamide elimination | usable |
| Ea_INTg | **117.5 ± 25.2** | – | 328.2 ± 61.3 | kJ/mol | initial Maillard | usable (glucose) |
| Ea_INTf | – | 149.1 ± 87.7 | 180.7 ± 21.1 | kJ/mol | initial Maillard | wide HPD |
| *Ea_M* | *105.7 ± 29.1* | *90.1 ± 27.8* | *64.7 ± 17.8* | kJ/mol | — | **REFUSE** |
| *Ea_B* | *180.3 ± 38.5* | *119.6 ± 44.3* | *34.7* ᶜ | kJ/mol | — | **REFUSE** |
| *Ea_Cg* | ***−6.7 ± 0.2*** | – | *Indeterminate* | kJ/mol | — | **REFUSE — NEGATIVE Ea** |
| *Ea_Cf* | – | *152.7 ± 44.1* | *124.3 ± 10.5* | kJ/mol | — | **REFUSE** |
| Ea_Asp | 105.4 ± 10.6 | 108.3 ± 11.4 | 109.4 ± 16.1 | kJ/mol | deamidation | usable; remarkably sugar-independent |
| *Ea_X* | ***668.9 ± 35.2*** | *167.6 ± 68.8* | *322.2 ± 296.0* | kJ/mol | — | **REFUSE — physically absurd** |
| Ea_hy | – | – | 140.6 ± 8.5 | kJ/mol | sucrose hydrolysis | usable |
| Ea_I | – | – | 113.2 ± 13.8 | kJ/mol | isomerisation | usable |
ᶜ = fixed, not estimated.
**Authors' own caveat to carry: "a linear Arrhenius curve was assumed"; deviations are
worst at 120–140 °C and "could possibly be a consequence of a break in the Arrhenius
curve."**

### 2c. De Vleeschouwer et al. 2009 **Part II** (`devleeschouwer2009.pdf`), Table 3 (p. 542) — **[F]**, ± = **95% HPD**
**Shared conditions:** as Part I, **plus an equimolar third amino acid**; a_w 0.92,
120–200 °C, T_ref = 160 °C. Moisture differs by system (Gln 9.58%, control 14.53%,
Cys 19.65%). **[Cys] is measured as cysteine + cystine.**
**Added rate laws: k_INT2·[AA2]·[Glc] SECOND ORDER; k_E2·[Cys]·[AA] SECOND ORDER
(order assumed, never tested); k_Y·[Cys] first; k_Glu·[Gln] first.**
⚠️ Rows marked ᶜ are **fixed at Part I's control values** — no independent information.
⚠️ Part II does not repeat Part I's "no physical meaning" italics; they still apply.
⚠️ **Three typos in the Scheme-1 (glutamine) equations, verified at 400 dpi:** eq 2
`k_IN1`→`k_INT`; eq 7 `k_Int2·[Glc]·[Glc]`→`k_INT2·[Gln]·[Glc]`; eq 8
`k_B·[Int1]`→`k_F·[Int1]`.

| parameter | Gln (Scheme 1) | Cys (Scheme 2) | units | step | verdict |
|---|---|---|---|---|---|
| k_Fref | **8.05 ± 0.90** | 3.57 ᶜ | ×10⁻³ min⁻¹ | Int1 → acrylamide | usable (Gln) |
| k_Eref | **0.36 ± 0.05** | 0.10 ᶜ | min⁻¹ | acrylamide → DP | usable (Gln) |
| **k_E2ref** | – | **49.36 ± 1.18** | **M⁻¹ min⁻¹** | **acrylamide + cysteine → adduct** | **★ directly usable, tightest parameter in the basket (2.4% RSE)** |
| k_INTref | 1.70 ᶜ | 1.70 ᶜ | M⁻¹ min⁻¹ | Asn + Glc | fixed |
| k_INT2ref | 0.01 ± 0.00 | 0.26 ± 0.02 | M⁻¹ min⁻¹ | AA2 + Glc → melanoidin (LUMPED) | **apparent only — authors say NOT comparable to k_INT** |
| k_Mref / k_Bref / k_Cref / k_Aspref | 1.23 ᶜ / 3.90 ᶜ / Indet. ᶜ / 26.43 ᶜ | same | — | — | fixed; several have no physical meaning |
| k_Xref | Indeterminate ᶜ | 0.04 ± 0.01 | min⁻¹ | Asp consumption | Cys value estimated |
| k_Gluref | 0.62 ± 0.33 | – | min⁻¹ | Gln → glutamic acid | usable, wide |
| k_Yref | – | 0.35 ± 0.01 | min⁻¹ | unspecified Cys sink | usable but unspecified product |
| Ea_F | **124.1 ± 9.3** | 159.2 * | kJ/mol | acrylamide formation | Gln usable; Cys marker undefined (treat as fixed) |
| Ea_E | **92.4 ± 12.0** | 113.2 ᶜ | kJ/mol | acrylamide elimination | usable (Gln) |
| **Ea_E2** | – | **51.3 ± 1.5** | kJ/mol | **acrylamide + cysteine** | **★ directly usable — HALF of Ea_E** |
| Ea_INT2 | 13.2 ± 4.3 | 30.3 ± 1.6 | kJ/mol | lumped AA2+Glc | apparent only |
| Ea_C | **−6.7** ᶜ | **−6.7** ᶜ | kJ/mol | caramelisation | **REFUSE — negative** |
| Ea_X | 668.9 ᶜ | 97.2 ± 8.3 | kJ/mol | Asp consumption | **REFUSE the Gln value** |
| Ea_Glu | 35.9 ± 6.8 | – | kJ/mol | Gln deamidation | usable |
| Ea_Y | – | 110.5 ± 8.5 | kJ/mol | unspecified Cys sink | usable |

**Derived [K1]:** the cysteine channel equals the basic elimination channel when
[Cys] = k_E/k_E2 = 0.10/49.36 = **2.0 mM**. At the ~1–2 M cysteine loading used, the
cysteine channel outruns the basic one by **~1000×** — which is the ">99% reduction" of
the abstract.

### 2d. Quan et al. 2020, Table 1 — **[F]**, **88 constants, ALL UNUSABLE AS PRINTED**
**Conditions:** aqueous, **0.1 M phosphate pH 7.0**, 10 mL sealed, **100 and 130 °C only**,
7 timepoints 0–21 min, n = 3, SIDA quantification (d₄-CML, d₄-CEL, ¹³C₃-acrylamide),
ordinary non-linear least squares (NOT the determinant criterion).
**Four disqualifying defects: (1) NO UNITS on any k; (2) NO RATE LAWS — Supplementary
Figure 1 is absent from the PDF, so every reaction ORDER is unknown; (3) the "Table 1
Continued" headers are MISLABELLED (K₁–K₅ repeated instead of K₆–K₁₁); (4) no ± SD
printed despite a footnote claiming them.** Authors state explicitly: "**specific
activation energies cannot be estimated**".
**Column mapping recovered by me from six text cross-references (all pass) — see dossier.**

| Group / T | k₁ GO ×10⁻³ | k₂ MGO ×10⁻² | k₃ CML ×10³ | k₄ CEL ×10² | k₅ AA(GO) ×10² | k₆ acetald ×10⁻³ | k₇ acrolein ×10⁻² | k₈ AA(acro) ×10³ | k₉ harmane ×10² | k₁₀ norharm ×10² | k₁₁ melanoidin ×10⁻³ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Lys-Glu/100 | 8.34 | 0.66 | 2.94 | 0.25 | – | 0.27 | 1.50 | – | – | – | 1.41 |
| Trp-Glu/100 | 0.02 | 0.20 | – | – | – | 0.16 | 0.94 | – | – | – | 0.11 |
| Asn-Glu/100 | 1.80 | 0.16 | – | – | – | 0.15 | 0.63 | – | – | – | 0.59 |
| Mix-Glu/100 | 5.58 | 0.68 | 1.19 | 0.19 | – | 0.19 | 0.79 | – | – | – | 1.76 |
| Lys-Glu/130 | 8.15 | 3.70 | 6.11 | 1.25 | – | 1.58 | 11.4 | – | – | – | 1.58 |
| Trp-Glu/130 | 3.15 | 1.55 | – | – | – | 1.09 | 9.58 | – | 4.47 | 8.35 | 0.65 |
| Asn-Glu/130 | 4.43 | 3.16 | – | – | 6.05 | 5.77 | 445 | 4.41 | – | – | 1.21 |
| Mix-Glu/130 | 4.67 | 4.05 | 1.91 | 0.77 | 1.09 | 5.03 | 405 | 4.02 | 4.79 | 7.84 | 1.66 |
R² across all 88 fits: **0.75–0.99** (the paper says 0.80–0.99; eight fits are below 0.80).
**Verdict for every cell: NOT USABLE — no units, no order. Re-fetch the version of record.**
Additional confound: **lysine is 30 mmol in Lys/Glu but 100 mmol in Mix/Glu (3.3×)**, so
the competition comparison is not controlled. Two of the paper's six quantitative claims
are contradicted by its own table (k₁ and k₂ are *not* lower in Mix).

---

## 3. EQUILIBRIUM BINDING (not kinetics) — from the mis-filed `Zhou2023.pdf`
Zhu et al. 2023, *Food Chem* 400:133854 (kafirin/Baijiu). **[F]** from UV double-reciprocal
plots. **Conditions: 53% aqueous ethanol, kafirin 200 mg/L (8.1 µM, MW 24652 Da),
25/45/65 °C.**

| ligand | K(298 K) | K(318 K) | K(338 K) | ΔH kJ/mol | ΔS J/mol/K | units | verdict |
|---|---|---|---|---|---|---|---|
| BMFD | 1934.33 | 1448.00 | 1214.00 | −9.79 | +29.96 | M⁻¹ | usable as equilibrium, 53% EtOH only |
| DT | 327.42 | 515.21 | 806.88 | +14.22 | +97.95 | M⁻¹ | **ΔH inconsistent with its own K's (van't Hoff gives +18.3)** |
| **MFO (= MFT)** | **1092.66** | **848.23** | **404.24** | **−20.58** | **−10.17** | M⁻¹ | **usable — the only MFT binding constant in the basket; van't Hoff checks exactly** |
| MMFD | 345.36 | 707.76 | 863.90 | +19.40 | +114.28 | M⁻¹ | usable |
ΔG (kJ/mol) at 298/318/338 K: BMFD −18.72/−19.32/−19.92; DT −14.96/−16.92/−18.88;
MFO −17.55/−17.34/−17.14; MMFD −14.66/−16.94/−19.23. R² 0.948–0.998.
**Sign-crossing: furan-ring ligands (BMFD, MFO) bind WEAKER on heating; polysulfides (DT,
MMFD) bind STRONGER.**
Associated **[M]** headspace suppression by 200 mg/L kafirin, 35 °C: BMFD −64.28%,
**MFO −54.74%**, MMFD −37.31%, DT −22.05% (all p<0.01).
Associated **[M]** odour thresholds in 53% aq. ethanol: BMFD 1.141 ± 0.105 → 2.137 ± 0.086;
DT 0.710 ± 0.066 → 1.461 ± 0.094; **MFO 0.657 ± 0.058 → 2.185 ± 0.131 µg/L (+232.57%)**;
MMFD 0.901 ± 0.073 → 2.091 ± 0.177.

---

## 4. PAPERS WITH NO KINETIC PARAMETERS — stated plainly

| paper | what it has instead | usability |
|---|---|---|
| **Martin & Ames 2001** | Relative GC peak areas at ONE time (2 min) and a falling 180→160 °C profile. | **Ordinal only.** Best datum: total pyrazines in the six-amino-acid mix = **38%** of the sum of the six single systems (my recompute 36.5%) — **but the mix is sugar-limited at 1.75 mol AA per mol glucose**, so it is partly stoichiometric, and the authors prove it (adding fructose raises the yield). The part that survives: **Strecker aldehyde selectivity in the mixture — leucine 1.24× UP, isoleucine 0.41×, phenylalanine 0.37×, methionine 0.19×.** |
| **Hwang, Hartman & Ho 1995b** | ¹⁵N isotope partitioning + yields in µg/g glucose at ONE time (1 h), ONE T (180 °C), ONE pH (7). **No replication, no error bars, no statistics anywhere.** | **Ordinal, but with hard stoichiometric constraints.** Best datum: **glycine's OWN N-compound yield rises 1.8× with lysine present and falls to 0.35× with arginine present** — a value >1 that no shared-pool competition model can reach; it requires a lysine-catalysed sugar-fragmentation channel. |
| **Hemmler et al. 2018** | Counts of molecular formulae by FT-ICR-MS and relative ESI(−) intensities. pH **uncontrolled and drifts 2–3 units**. | **Ordinal over a proxy (formula count).** Reactivity order **lysine > cysteine > isoleucine ≈ glycine** (10 h, 100 °C) — but it is time-dependent (Lys/Gly ≈ 17× at 2 h, ≈ 2× at 10 h) and the **browning order is different** (Lys > Ile > Gly >> Cys). Best structural datum: **sugars drive rates, amino acids drive product identity** (45–95% formula overlap across six sugars). One [K1] number — ARP decay ~0.135 h⁻¹, t½ ≈ 5.1 h — **which I recommend refusing** (ESI intensity ≠ concentration; 2 points; pH drifting). |

---

## 5. TOP 10 BY VALUE TO A MASS-ACTION KINETIC CORE

**1. k(FFT sink) ≈ 9 × 10⁻⁴ s⁻¹ at 30 °C, ~10 g/L coffee solids, pH 5–6.**
Hofmann Fig. 6 [K1] 9.4 × 10⁻⁴ and Table 2 [K1] 9.8 × 10⁻⁴, corroborated by
Charles-Bernard's independent bound > 7.70 × 10⁻⁴ at 25 °C. *This single number converts
the repo's "no thiol consumption at all" defect into a bounded one.* First order in thiol,
confirmed against the data, not assumed. **Highest priority in the whole basket.**

**2. k_E2 = 49.36 ± 1.18 M⁻¹ min⁻¹ at 160 °C, Ea_E2 = 51.3 ± 1.5 kJ/mol** — acrylamide +
cysteine, second order. De Vleeschouwer Part II. Already bimolecular, already in
M⁻¹ min⁻¹, already Arrhenius-parameterised at a stated T_ref, and the **tightest
parameters in the basket** (2.4% and 2.9% relative HPD). Generalises beyond acrylamide to
any thiol/Michael-acceptor pair.

**3. Claeys 2005 control set: k_F = 0.451 ± 0.023 × 10⁻³ min⁻¹, k_E = 111.1 ± 8.9 × 10⁻³
min⁻¹, Ea_F = 168.25 ± 3.80, Ea_E = 167.21 ± 4.30 kJ/mol (T_ref 160 °C).**
The cleanest complete two-step first-order set with SEs, fit diagnostics and a
non-isothermal correction — and it **self-validates** against its own Figure 2 (plateau
2886 vs ~2800 ppb; t_max 49.8 vs ~40–50 min). Ingest the **ratio** (correlation 0.94–0.98).

**4. The structural claim that Ea(thiol scavenging) is HALF Ea(generic elimination)**
— 51.3 vs 113.2 kJ/mol (De Vleeschouwer II). The two channels cross over in temperature,
so a core that lumps all product loss into one Arrhenius term gets the T-dependence wrong
in any protein-bearing matrix. Tightly determined, and mechanistically portable.

**5. Charles-Bernard's full thiol k ladder + the two-channel decomposition.**
Ten first-order constants spanning **755×** (2M2P 1.02 × 10⁻⁶ → FFT > 7.70 × 10⁻⁴ s⁻¹),
plus Table 3 showing benzylic thiols retain ~25% of their rate with all electrophiles
blocked while aliphatics do not. Gives both magnitudes and the branch structure.

**6. De Vleeschouwer Part I glucose control (non-italic subset): k_INT = 1.70 ± 1.05
M⁻¹ min⁻¹ (SECOND ORDER, Asn + Glc), k_F = 3.57 ± 1.38 × 10⁻³ min⁻¹, k_E = 0.10 ± 0.04
min⁻¹, k_Asp = 26.43 ± 5.76 × 10⁻³ min⁻¹; Ea 117.5 / 159.2 / 113.2 / 105.4 kJ/mol.**
The only genuinely **bimolecular** Maillard-initiation constant in the basket, from a
proper multiresponse determinant-criterion fit. Matrix-specific to a_w 0.92 and ~3 M.

**7. The matrix electrophile capacity: 8–10 mmol sites per g dry coffee solids**
(Charles-Bernard, hydroxylamine titration) — with Hofmann's covalent lower bound of
≥ 0.028 mmol FFT bound per g melanoidin. The only published site density; it converts the
pseudo-first-order thiol sink into a bimolecular one that scales with matrix loading
([K1] k₂ ≈ 1.3 × 10⁻³ M⁻¹ s⁻¹ for EtSH, ≳ 8.6 × 10⁻³ for FFT).

**8. Hofmann's disulfide control: < 1.5% of the thiol flux goes to RSSR** (< 6 µg
disulfide from 400 µg FFT consumed), plus the DTE-irreversibility result and the ²H-NMR
line broadening. **Kills an entire family of wrong sink implementations** and fixes the
stoichiometry as 1:1 thioether addition — first order in thiol, not second.

**9. Hwang's lysine synergy: glycine's own N-compound yield rises 1.8× with lysine and
falls to 0.35× with arginine.** A value **> 1** is unreachable by any shared-substrate
competition scheme; it forces a **lysine-catalysed sugar-fragmentation channel** into the
network. Structural, not numerical, but it constrains architecture rather than parameters.

**10. The sign-crossing validation rows.** Five of them, all from this basket, all
falsifiable, none reproducible by a smooth family-level factor:
(i) **fructose > glucose at a_w 0.92 but glucose > fructose in dilute solution**
(De Vleeschouwer I vs Claeys); (ii) **glutamine's promotion of acrylamide grows with T in
liquid (155%→322%) but shrinks with T at a_w 0.92 (267%→120%)** — same lab, same amino
acid, neither paper remarks on it; (iii) **glutamine's k_E Arrhenius crosses the control's
at ~180 °C** (stated explicitly by Claeys); (iv) **cysteine ranks 2nd by MRP count but
last by browning** (Hemmler, ~90× below lysine); (v) **furanthiol–protein binding weakens
with T while polysulfide binding strengthens** (Zhu 2023).

---

## 6. CORPUS GAPS THIS WAVE EXPOSES

1. **`Zhou2023.pdf` is the wrong paper.** It contains Zhu et al. 2023 *Food Chem*
   400:133854 (kafirin/Baijiu). **Zhou et al. 2023 JAFC — purified deoxyxylulosyl-alanine
   Amadori compound + cysteine, pH 6/7/8, MFT rise-then-fall — is NOT in the corpus.**
   It was the only pH-resolved MFT source and the only purified-Amadori source assigned to
   me. **Re-fetch, and rename the existing file so no wave cites it as Zhou.**
2. **Hwang's pyrazine twin (JAFC 1995, 43, 179−184)** is absent. `hwang1995b.pdf` is the
   pyridine/pyrrole/oxazole paper; the pyrazine data reach me only as an aggregate in its
   Figure 4.
3. **Quan 2020's Supplementary Figure 1** (the rate equations) and Figure 2 are absent,
   which is why 88 fitted constants cannot be used. Also, the file is a **pre-proof** —
   obtain the version of record, where the mislabelled Table 1 headers may be fixed.
4. **Hemmler 2018's SI (Figs S1–S5)** is absent; the six-sugar comparison and the
   compositional-space counts reach me only through the main text.
5. **Laroque et al. 2008, *Food Chem* 111:1032**, "Kinetic study on the Maillard reaction.
   Consideration of sugar reactivity" — cited by Hemmler as the source of the pentose >
   hexose rate claim. An actual kinetic paper on sugar reactivity, not in the corpus.

## 7. ERRATA FOUND — record these so no later wave re-imports them

| paper | erratum |
|---|---|
| Charles-Bernard 2005 | Table 2 / Fig. 11 units printed "mol⁻¹ s⁻¹"; **the values are s⁻¹** (proven against 5 of the paper's own half-lives). |
| Hofmann & Schieberle | Filename says 2001; **published 2002**, JAFC 50:319−326. Table 1's incubation temperature is stated three different ways (30 / 40 / 45 °C). |
| De Vleeschouwer Part I | k_Fref unit printed "(10⁻³ **mm**⁻¹)" — means min⁻¹. |
| De Vleeschouwer Part II | **Three equation typos** (eqs 2, 7, 8 of the glutamine scheme), verified at 400 dpi. Part I's "no physical meaning" italics are not repeated in Part II's table though they still apply. Cys Ea_F carries an undefined footnote marker. |
| Hwang 1995b | Abstract and body claim **asparagine had the highest oxazole contribution — asparagine produced NO oxazoles** (Table 1 and Fig. 3). Table 1 unit header "mg/g of glucose" should be **µg/g** (proven: figure bars equal the column sums exactly). |
| Quan 2020 | Mislabelled continuation headers; no units; no SD despite the footnote; the k₁/k₂ competition claim is contradicted by the table; "melanoidin has the highest k" contradicted by the printed prefixes; R² range understated. |
| Zhu 2023 (`Zhou2023.pdf`) | DT's ΔH (14.22 kJ/mol) disagrees with a van't Hoff fit to its own K values (+18.3). |
