# Wave K4c — CONSOLIDATED PARAMETER TABLE + PROPOSED FIT/HOLD-OUT ROLES

**DRAFT FOR THE ORCHESTRATOR. This file does NOT amend
`docs/reference/FIT_HOLDOUT_DECLARATION.md` and must not be treated as a declaration.**
Extraction date 2026-08-28. Companion dossiers:
`kang2026_extraction.md`, `bell1995_extraction.md`, `nakamura2020_extraction.md`.

---

## 0. PAPER-IDENTITY CHECK — all three PDFs are the expected papers

| file | expected (brief) | **actual** | verdict |
|---|---|---|---|
| `Kang2026.pdf` | 10.1039/D5FB00932D, Sustainable Food Technol., Cys–Xyl MRI, 100/120/140 °C | **DOI, journal, volume and page range all match** — Kang et al., *Sustain. Food Technol.* 2026, **4**, 3239–3252 | ✔ **correct paper**, ⚠ **topic partly mis-described in the brief — see §1** |
| `bell1995.pdf` | 10.1016/0963-9969(95)00052-6, Bell, non-enzymatic browning in amorphous solids | **DOI matches.** Bell, L. N., *Food Research International* **28**(6), 591–597 | ✔ **correct paper**, ⚠ **printed year is 1996, not 1995** |
| `nakamura2020.pdf` | 10.1093/bbb/zbaa019, xylose + 4-hydroxy-5-methyl-3(2H)-furanone browning | **DOI matches.** Nakamura, Mikami, Noda & Murata, *Biosci. Biotechnol. Biochem.* **85**(2), 401–410 | ✔ **correct paper**, ⚠ **issue of record is 2021** |

**No mis-filed PDF. Two filename-year mismatches (bell1995 → 1996, nakamura2020 → 2021) that
should be carried as aliases in the bibliography.**

---

## 1. THE THREE FINDINGS THE ORCHESTRATOR NEEDS BEFORE ANYTHING ELSE

### 1.1 ⚠⚠ Kang 2026 — units: the PDF text layer renders **μ as m**, corrupting every concentration by 1000×

Verified against 500–900 dpi rasters. **Table 1's printed unit is μg L⁻¹, not mg L⁻¹.**
Likewise the internal standard is **3 μL of 0.018 μg μL⁻¹**, not 3 mL of 0.018 mg mL⁻¹;
DPPH is **100 μmol L⁻¹**, not 100 mmol L⁻¹; the SPME fibre is **75 μm**; the GC film **0.25 μm**.
Genuinely `m`: the 0.5–3.0 **mg** mL⁻¹ antioxidant range, the 100 **mmol** L⁻¹ storage solution,
the 10 **mmol** L⁻¹ TTCA reaction solution, the 5 **mmol** L⁻¹ chelation reagent.

> **ACTION: if anything from Kang 2026 has already been ingested at "mg L⁻¹", it is wrong by
> three orders of magnitude.** Any other paper in the corpus processed through the same
> pdftotext path with an RSC/Arbortext PDF should be re-checked for the same corruption.

### 1.2 ⚠⚠ Kang 2026 — the 100/120/140 °C per-compound data are **NOT IN THE PDF**

The brief anticipated *"the sulfur branch's first temperature ladder — extract every
time-temperature point."* The ladder exists as an experiment but its numeric table is
**Table S4 (Supplementary Information)**, and the pH ladder is **Table S5 (SI)**. Neither is in
`Kang2026.pdf`. The main text carries only **Fig. 1a/1b (compound-class totals, digitised in the
dossier)** plus **five individual numbers** quoted in prose. **MFT and FFT at 100/120/140 °C are
mentioned only qualitatively — no numbers.** The standard curves (Table S3), the a_w/salt map
(Table S2), and the H₂S-release evidence (Figs S2/S3) are all SI-only too.

> **ACTION: fetch the SI from https://doi.org/10.1039/d5fb00932d before assigning Kang 2026 any
> temperature-ladder role.** (The paper's own Data-availability line cites a broken SI DOI,
> `10.1039/d500932d`.) **Until then, Kang 2026 cannot serve as the sulfur branch's temperature
> ladder — it can only serve the roles proposed in §4.**

### 1.3 ✔ Kang 2026 — the brief's CRITICAL CHECK, answered: **absolute, not peak areas**

Quantitation is by **1,2-dichlorobenzene internal standard**, in two tiers that must be kept
distinct in the repo:
- **Tier A** — compounds marked `*` in Table 1: **external standard curve + IS calibration**.
  **Both MFT and FFT are Tier A.** These are genuine absolute μg L⁻¹.
- **Tier B** — unstarred compounds: **IS-only semi-quantitation**, response factor assumed = 1.
  Ratio-comparable within a compound across conditions; absolute magnitude carries an
  unquantified (plausibly 2–10×) error.
**Nothing in Kang 2026 is peak-area-ordinal.** But every class subtotal and grand total mixes
the two tiers and inherits Tier B's error.

---

## 2. CONSOLIDATED PARAMETER TABLE — everything the three papers yield

Codes: **[M]** printed by the authors · **[D]** digitised by me (per-panel error in the dossier)
· **[Z]** derived by me, never printed · **[S]** exists only in Supplementary Information.

### 2a. Kang et al. 2026 — TTCA/ARP, Cys–Xyl intermediates

| # | quantity | value | conditions | code | anchor |
|---|---|---|---|---|---|
| K1 | **Per-compound volatile table**, 48 compounds × 9 amino-acid systems | **236 non-ND cells, μg L⁻¹, mean ± SD, n = 3** | TTCA 10 mM + equimolar AA, **pH 7.0, 120 °C, 120 min**, sealed | [M] | Table 1, pp. 3244–3246 |
| K2 | MFT (2-methyl-3-furanthiol) | 1.388 / 4.083 / **11.569** / 0.329 / — / 7.327 / 0.344 / 2.536 / **18.752** μg L⁻¹ (Control/Cys/Gly/Ala/Met/Glu/Asp/His/Lys) | as K1 | [M] **Tier A** | Table 1 |
| K3 | FFT (2-furfurylthiol) | 4.107 / 11.591 / 15.254 / 2.133 / ND / 12.107 / 3.176 / 0.380 / **29.295** μg L⁻¹ | as K1 | [M] **Tier A** | Table 1 |
| K4 | **FFT/MFT ratio spans 0.15–9.23** across nine nitrogen co-substrates at fixed T, pH, t | 2.96 (control) → 1.32 (Gly) → 0.15 (His) → 1.56 (Lys) | as K1 | [Z] | dossier §4h |
| K5 | Total volatiles per system | 43.155 / 49.878 / 59.793 / 5.027 / 13.503 / 41.717 / 8.911 / 9.240 / **105.637** μg L⁻¹ | as K1 | [M] | Table 1, p. 3246 |
| K6 | Total compound kinds per system | 28 / 26 / 32 / 11 / 9 / 26 / 13 / 12 / 33 | as K1 | [M] | Table 1 |
| K7 | **Volatile class totals vs temperature** (thiols, thiophenes, thiazoles, N-het, O-het, total) | 6 classes × 3 T; total **17.7 / 42.3 / 73.4 μg L⁻¹** at 100/120/140 °C | TTCA 10 mM, pH 7.0, 120 min | [D] ± 0.3 | Fig. 1a, p. 3243 |
| K8 | **Volatile class totals vs pH** | 6 classes × 3 pH; total **20.3 / 17.5 / 20.4 μg L⁻¹** at pH 5.5/7/8 | TTCA 10 mM, 120 °C, 120 min | [D] ± 0.1 | Fig. 1b |
| K9 | Furfural ladder | **3.345 (100 °C) → 11.039 (140 °C) μg L⁻¹**, ratio 3.3 | pH 7.0, 120 min | [M] + [Z] | p. 3242 |
| K10 | Thiazoles at pH 8 | **10.952 μg L⁻¹**; 2-acetylthiazole = **27.01 %** of it (= **2.958 μg L⁻¹** [Z]) | 120 °C, 120 min | [M] | p. 3242 |
| K11 | Pyrazines detected **only at pH 8** | **0.838 μg L⁻¹** | 120 °C, 120 min | [M] | p. 3242 |
| K12 | **Measured pH drift during reaction** | initial 5.5 / 7.0 / 8.0 → final **5.1 / 4.9 / 4.7** | 120 °C, 120 min | [M] | p. 3242 |
| K13 | Browning A₄₂₀, 15 amino-acid systems | 0.072 (Cys) – **0.774 (Lys)**; control 0.106 | pH 7.0, 120 °C, 120 min | [D] ± 0.005 | Fig. 1c |
| K14 | **Storage residual, 6 designs × 3 levels × 7 time points** | **126 points**, mmol L⁻¹ from 100 | 60 d; T 4/25/40 °C; pH 5.5/7/9; a_w 0.113/0.432/0.843 | [D] ± 1 | Fig. 2a/b/e/f/i/j |
| K15 | 60-day losses | TTCA 7.06 % (40 °C), 11.19 % (pH 9), 13.5 % (a_w 0.113), **35.77 % (a_w 0.843)**; ARP 12.17 / 21.25 / — / **60.61 %** | as K14 | [M] | pp. 3239, 3247 |
| K16 | **18 apparent first-order storage rate constants** | **1.44 × 10⁻⁴ – 1.50 × 10⁻² d⁻¹** (1.7 × 10⁻⁹ – 1.74 × 10⁻⁷ s⁻¹) | as K14 | **[Z]** | dossier §6c |
| K17 | Apparent Eₐ of storage degradation | **TTCA 43.1 kJ mol⁻¹ (R² 0.947); ARP 18.7 kJ mol⁻¹ (R² 0.808)** | 4–40 °C solution | **[Z]** | dossier §6c |
| K18 | **a_w sensitivity of degradation** | **d(ln k)/d a_w = 1.42 (TTCA), 1.79 (ARP)**; ≈ 15 % / 20 % per 0.1 a_w | solid, 25 °C, 60 d | **[Z]** | dossier §6c |
| K19 | Solid-state colour after 60 d | L\*, a\*, b\*, ΔE at 3 a_w × 2 compounds; **ΔE 9.36 → 23.90 (TTCA), 12.16 → 46.62 (ARP)** | as K14 | [M], all 6 ΔE reproduce exactly [Z] | Table 2, p. 3248 |
| K20 | **ARP browns more per mole lost than TTCA** | 1.69× the mass loss but 1.95× the ΔE at a_w 0.843 | 60 d | **[Z]** | dossier §7 |
| K21 | Apparent Eₐ of volatile accumulation (yield-at-120-min exponents) | thiols 36.0, thiophenes 70.8, thiazoles 43.6, O-het 41.8, **total 45.7 kJ mol⁻¹**; furfural alone **38.3** | 100–140 °C | **[Z] — NOT mechanistic** | dossier §5d |
| K22 | Antioxidant single values | Fe²⁺ chelation TTCA **71.35 %**, ARP **51.03 %**; lipid-peroxidation inhibition TTCA 31.30 %, ARP 29.18 %, MRPs 36.39 % at 3.0 mg mL⁻¹ | assay conditions in dossier §2e | [M] | pp. 3248–3249 |
| K23 | DB-Wax retention indices | 48 compounds, RI_calc and RI_lit | DB-Wax 30 m, C7–C30 | [M] | Table 1 |
| — | **Table S4/S5 per-compound ladders, Table S3 standard curves, Table S2 a_w map, Figs S2/S3 H₂S** | — | — | **[S] NOT IN PDF** | §1.2 |

### 2b. Bell 1996 (`bell1995`) — non-enzymatic browning in amorphous PVP

| # | quantity | value | conditions | code | anchor |
|---|---|---|---|---|---|
| B1 | **13 pseudo-zero-order browning rate constants** with **95 % CL** and **R²** | **0.0053 ± 0.0016 → 0.035 ± 0.004 OD/g/d** | **25 °C, pH 7, 0.1 molal phosphate, 1 molal Glc + 1 molal Gly**, PVP matrix | [M]/[F] | Table 2, p. 594 |
| B2 | **13 glass transition temperatures** (DSC **onset**, 5 °C/min) | **+55 → −76 °C** | as B1 | [M] | Table 2 |
| B3 | T − T_g for all 13 systems | **−30 → +101 °C** | as B1 | [Z] (matches Fig. 2 abscissa) | dossier §4 |
| B4 | **4 rate constants at non-1-molal reactant** | **0.0034 ± 0.0013 → 0.068 ± 0.011 OD/g/d** at **0.28 → 4.9 molal** | 25 °C, pH 7 | [M]/[F] | Table 3, p. 595 |
| B5 | Full composition of 12 systems | a_w, % moisture db, PVP g, sorbed water g, added solution μL, buffer molal, reactant molal | — | [M], self-consistent [Z] | Table 1, p. 593 |
| B6 | **Glass-transition effect, FULLY DECONFOUNDED** (a_w 0.54, moisture 18.3–18.5 % db, 1.1 molal all fixed) | **k = 0.014 (glassy) → 0.024 → 0.034 (rubbery) = 2.4×** | 25 °C | **[Z]** | dossier §4a, §6 |
| B7 | Glass-transition effect, author's headline | *"7-fold"* = **k_max/k_min = 6.6×**; **mean glassy → mean rubbery = 3.0×** | 25 °C | [M] + [Z] | dossier §4b |
| B8 | **No resolvable Tg effect INSIDE the glass** | Tg 38 / 43 / 55 °C at a_w 0.33 → k 0.0053 / 0.0073 / 0.0054, **p > 0.05** | 25 °C | [M] | p. 594 |
| B9 | **a_w effect at fixed 1 molal: rise then PLATEAU, no maximum** | 0.0053 (a_w 0.33) → 0.034 (0.54) then **flat 0.023–0.035 to a_w 0.96** | 25 °C, PVP-LMW | [M] + [Z] | Table 2, Fig. 3 |
| B10 | **Apparent reaction order in reactant molality** | **n = 1.40** (K15, a_w 0.33, 1.0 → 4.9 molal) and **n = 1.05** (K30, a_w 0.76, 0.28 → 1.0 molal) — **not 2** | 25 °C | **[Z]** | dossier §5a |
| B11 | Effect-size ranking | **concentration 20× ≫ Tg (deconfounded) 2.4×**; a_w 0.33→0.54 6.4×, a_w > 0.54 flat | 25 °C | **[Z]** | dossier §6 |
| B12 | **Non-zero t = 0 pigment** | intercept **0.31–0.50 OD/g** = 9–19 % of the 67-day absorbance | after 1 week equilibration at 25 °C | **[Z]** | Fig. 1, dossier §7 |

### 2c. Nakamura et al. 2021 (`nakamura2020`) — xylose / HMFO browning

| # | quantity | value | conditions | code | anchor |
|---|---|---|---|---|---|
| N1 | **HMFO formation + browning, 2 pH × 5 times** | A400 0 → **4.02** (pH 6.5) / **3.47** (pH 8.0); HMFO peak **1.91 mM at 3 h** (pH 6.5) / **0.59 mM at 0.5 h** (pH 8.0) | **13 mM Xyl + 27 mM Lys, 0.2 M Na/K phosphate, 100 °C** | [D] ± 0.02, **4 text cross-checks pass** | Fig. 2, p. 404 |
| N2 | **HMFO molar yield from xylose** | **14.7 % (pH 6.5), 4.5 % (pH 8.0)** | as N1 | **[Z]** | dossier §4a |
| N3 | **pH effect on browning INVERTS with time** | pH 8.0 is **1.89×** faster over 0–0.5 h, **0.62×** as fast over 1–6 h; curves cross at ≈ 3 h | as N1 | **[Z]** | dossier §4a |
| N4 | **HMFO first-order decay, 5 conditions** | **0.094 h⁻¹** (Xyl-Lys pH 6.5) · **0.328 h⁻¹** (Xyl-Lys pH 8.0, 0.5–3 h) · **1.433 h⁻¹** (pure HMFO +OPD) · **2.062 h⁻¹** (pure HMFO −OPD) = **2.6 × 10⁻⁵ – 5.7 × 10⁻⁴ s⁻¹** | 100 °C, pH 6.5/8.0 | **[Z]** | dossier §4a, §5c |
| N5 | **MGO + DA from pure HMFO, 2 trapping modes × 5 times** | cumulative at 6 h: MGO **0.235 mM**, DA **0.244 mM** from 1.727 mM HMFO = **27.7 % molar** | **1.8 mM HMFO, pH 8.0, 100 °C** | [D], **3 text cross-checks pass** | Fig. 3 |
| N6 | **Cumulative ≫ standing pool** | MGO **≈ 120×**; 1-DX **≥ 500×** | pH 8.0, 100 °C, 6 h | **[Z]** | dossier §5c, §7c |
| N7 | **Seven-dicarbonyl inventory × 5 times** (HMFO, pentosone, xylosone, 1-DX, 3-DX, tetrosone, MGO, DA) | **1-DX 9.84 mM** and **MGO 4.16 mM** at 6 h; all others ≤ 1.5 mM | 13 mM Xyl + 27 mM Lys + 28 mM OPD, pH 8.0, 100 °C | [D], 1-DX cross-check 0.4 % | Fig. 9a |
| N8 | **Molar branching from 13 mM xylose** | **1-DX = 76 %**, **MGO = 32 %** (C₃ basis), everything else ≤ 8 % | as N7 | **[Z]** | dossier §7a |
| N9 | Standing-pool dicarbonyls × 5 times | HMFO peak **0.378 mM**; MGO **0.150 mM**; DA **0.097 mM**; 1-DX ≤ 0.02 mM | as N7 but OPD after heating | [D] ± 0.01 | Fig. 9b |
| N10 | **Intrinsic browning ranking at matched ≈ 13 mM** | initial dA400/dt: **MGO+Lys 18.7 ≫ DA 10.6 ≈ DA+Lys 9.1 > HMFO 6.5 > MGO 5.8 > HMFO+Lys 4.2 ≫ Xyl+Lys 1.54 ≫ Xyl 0.32 h⁻¹** | **100 °C, pH 8.0, ± 27.4 mM Lys** | **[Z]** from [D] | Fig. 10, dossier §8c |
| N11 | **MGO+Lys browns 12.2× faster than Xyl+Lys per mole** | ratio 12.2 | as N10 | **[Z]** | dossier §8c |
| N12 | **Lys reduces HMFO browning to exactly half** | 4.09 → 2.04 A400 at 60 min (0.499×) | as N10 | [M] + [D], exact | Fig. 10b |
| N13 | **8 browning curves × 6–7 time points** | 50 points | as N10 | [D] ± 0.03 | Fig. 10 |
| N14 | **HMFO decomposition is O₂- and metal-dependent** | strongly repressed by **40 mM DTPA** and by **N₂ degassing 5 min** | 1.8 mM HMFO, pH 8.0, 100 °C, 3 h | [M], **qualitative — no % given** | Fig. 7a |
| N15 | New dicarbonyl intermediates identified | **2-hydroxy-3,4-dioxopentanal (C₅H₆O₄)** and **2,3-dioxopent-4-enal (C₅H₄O₃)**, exact masses 203.0824 / 185.0709 | RT, pH 8.0 | [M] | Figs 4–6 |
| N16 | Melanoidin MW window | **≈ 2000–5000**, eluting just before tannic acid (MW 1701) on GP-HPLC | pH 8.0 | [M] | Fig. 11 |
| N17 | **Fluorescence ≠ pigment** | *"fluorescent substances and pigments do not seem to be identical"*; fluorescence requires Lys, pigment does not | pH 8.0 | [M] | Fig. 11 |

---

## 3. WHAT THE THREE PAPERS COLLECTIVELY DO AND DO NOT ADD

**They add, jointly:**
1. **The corpus's cleanest deconfounded water-activity / glass-transition evidence** (Bell), plus
   a second, independent a_w series on a different observable — intermediate *disappearance*
   rather than pigment *formation* (Kang K18). **Two papers, two matrices, two observables, same
   axis.**
2. **A quantitative branching constraint on pentose degradation** (Nakamura N8: 1-DX 76 %,
   MGO 32 % from xylose at pH 8/100 °C) and **the flux-vs-pool correction that goes with it**
   (N6: pools understate flux by 10²–10³).
3. **An intrinsic browning-power ranking at matched molarity** (N10/N11), which is exactly the
   kind of relative constraint that survives the absence of an extinction coefficient.
4. **A large, arithmetically airtight per-compound volatile matrix** for a Cys–Xyl intermediate
   at one condition (Kang K1–K6), including **Tier-A MFT and FFT**.
5. **A demonstration that FFT/MFT branching is NOT a constant** (K4: 0.15–9.23 across nitrogen
   co-substrates at fixed T/pH/t).

**They do NOT add:**
- **Any activation energy.** Bell has one temperature (25 °C); Nakamura has one (100 °C); Kang's
  thermal work is single-endpoint at 120 min, so its "Eₐ" (K21) is a yield exponent, not a
  barrier. **Kang K17 (43.1 / 18.7 kJ mol⁻¹) is a *storage* Eₐ over 4–40 °C and belongs to a
  different regime entirely.**
- **Any molar pigment rate.** Bell's k is OD/g/d in an extract; Nakamura's A400 is
  dilution-multiplied. Both are relative.
- **Any measured H₂S, NH₃, or free-amine concentration.**

---

## 4. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT ONLY**

Vocabulary and reasoning style follow `docs/reference/FIT_HOLDOUT_DECLARATION.md`
(**FIT** / **HOLD-OUT** / **neither**). Module numbers refer to that document.
**★ marks the rows I judge highest-value.** No dataset below appears in two columns, and where
a cut falls inside a paper the cut axis is named.

### 4a. Module 4 (trunk / browning) and Module 7 (matrix correction) — **Bell 1996**

| dataset | **proposed role** | one-line reasoning |
|---|---|---|
| **★ Bell Table 2, the a_w 0.54 triple (k 0.014 / 0.024 / 0.034 at Tg 28 / 15 / 8 °C)** | **★ HOLD-OUT** | **the only fully deconfounded glass-transition measurement in the literature**: a_w, moisture and reactant molality are all identical across the three rows, so a model that gets the 2.4× span right has genuinely predicted a mobility term rather than absorbed it. Fitting it would destroy the one test worth having |
| **Bell Table 2, the a_w 0.33 triple (k 0.0053 / 0.0073 / 0.0054, p > 0.05)** | **FIT** | a **null** result — no resolvable Tg effect inside the glass across a 17 °C Tg span. Fitting a null costs no degrees of freedom and catches sign errors, exactly like the Bornhorst structural zero in Module 9 |
| **Bell Table 2, the seven PVP-LMW rows (a_w 0.33 → 0.96 at fixed 1 molal)** | **FIT** | the a_w axis at constant internal molality. This is the row set that should replace any bell-shaped a_w term the repo carries — **the maximum at a_w 0.55–0.7 is absent here by construction** |
| **★ Bell Table 3 + its 1-molal partners (apparent order n = 1.05 and 1.40)** | **★ HOLD-OUT** | a **functional-form** test, not a level: browning in a reduced-moisture matrix is **sub-bimolecular**. A model that assumes k ∝ [Glc][Gly] predicts n = 2 and fails by construction — which is the point. Cut axis: **reactant molality** (the 1-molal rows are in the FIT column, the 0.28 / 0.77 / 2.31 / 4.9-molal rows are out) |
| Bell Table 1 (composition), Table 2 Tg column | **reference table, NOT a scored dataset** | Tg and moisture are model *inputs*; scoring them is a category error (same treatment as the Zhou threshold table in Module 7) |
| Bell Fig. 1 intercepts (0.31–0.50 OD/g at t = 0) | **neither — an `initial_condition` flag** | it says the rate constants are slopes on an already 9–19 % browned system. Apply as a non-zero-initial-pigment note on B1/B4, do not score |

⚠ **Structural limitation the declaration should record:** Bell's design separates a_w from Tg
but **cannot separate Tg from polymer molecular weight** — they are perfectly collinear across
PVP-LMW / K15 / K30. PVP-K30 is the slowest system at every a_w, which is equally consistent
with an MW effect. **Any Tg term fitted or scored on Bell inherits this confound.**

### 4b. Module 4 (trunk / browning) and Module 9 (the Z3 norfuraneol block) — **Nakamura 2021**

Nakamura is the natural companion to the existing Module 9 norfuraneol entries
(Bornhorst 2017b Ea, Cerny 2007 branching) because **HMFO = 4-hydroxy-5-methyl-3(2H)-furanone
is the same compound**, here measured on the *pentose* rather than the potato/egg matrix.

| dataset | **proposed role** | one-line reasoning |
|---|---|---|
| **★ Nakamura Fig. 10b initial browning rates (MGO / DA / HMFO ± Lys, 6 curves)** | **★ HOLD-OUT** | the **intrinsic browning power of three named intermediates at matched molarity, temperature and pH**. A trunk fitted on sugar–amine browning must *predict* that MGO+Lys is 12× faster than Xyl+Lys and that Lys *halves* HMFO browning. Both are unfakeable by a fitted continuous term, and the second has a sign the model has no reason to get right |
| **Nakamura Fig. 10a (Xyl and Xyl+Lys, 7 time points each)** | **FIT** | the ordinary sugar–amine browning curve at 100 °C / pH 8.0, in the same experiment and on the same instrument as the hold-out above. Cut axis: **substrate identity** (xylose in, its three degradation products out) |
| **★ Nakamura Fig. 9a molar branching (1-DX = 76 %, MGO = 32 % of 13 mM xylose)** | **★ HOLD-OUT** | a **flux-partition** constraint the corpus has never had for a pentose. It scores whether the network routes xylose predominantly through the 1-deoxyosone, which is a topology question, not a level question |
| **Nakamura Fig. 9b standing-pool time course (5 species × 5 times)** | **FIT** | it identifies the *consumption* side that Fig. 9a alone cannot separate — the same logic that puts Guo 2020 retention in the FIT column of Module 8 |
| **★ Nakamura Fig. 3a/3b cumulative-vs-standing contrast (MGO ≈ 120×, 1-DX ≥ 500×)** | **★ HOLD-OUT (derived consistency test)** | not a level: it scores whether the model's **steady-state pool** for a fast intermediate is 2–3 orders below its **cumulative flux**. A model with a slow sink fails it. Directly analogous to the existing Zhou §A.3.3 cross-system pair |
| **Nakamura Fig. 2 HMFO decay constants (0.094 h⁻¹ at pH 6.5 vs 0.328 h⁻¹ at pH 8.0)** | **FIT (as priors on the HMFO sink)** | the corpus's only measured norfuraneol *disappearance* rates in a pentose system; they pair with Bornhorst 2017b's *formation* Ea, which is already FIT |
| **★ Nakamura Fig. 2 pH-inversion (pH 8 faster to 1 h, slower after; crossover at 3 h)** | **★ HOLD-OUT (shape test)** | the pH axis is the model's weakest (per the declaration's own §D.3 reasoning) and this is a **sign reversal in time**, which no monotone pH factor can produce. Same class as the Cerny 2007 Table 2 directional hold-out |
| **Nakamura Fig. 7 DTPA / degassing result** | **neither — a `requires_oxidant` gate, not a dataset** | it establishes that HMFO → MGO/DA is **oxidative and metal-catalysed** but quantifies nothing (chromatograms only). Treat exactly like the Hofmann 1996b oxidant finding: record the mechanism, do not score. ⚠ **It also means every Nakamura rate constant is an aerobic, unchelated rate** and must carry that flag |
| Nakamura Fig. 11 pigment characterisation, Figs 4–6/8 identifications | **reference only** | MW window 2000–5000, the 340 nm unknown, the colourless 280 nm polymer, the seven quinoxaline assignments. Qualitative; no scoring possible |
| Nakamura pH 6.5 dicarbonyl series | **unavailable** | *"data not shown"* (p. 407) |

### 4c. Module 1 (sulfur formation) — **Kang 2026**, with a gating condition

⚠ **RECOMMENDATION: assign Kang 2026 NO Module-1 role until the SI is in hand** (§1.2). What
follows is conditional.

| dataset | **proposed role** | one-line reasoning |
|---|---|---|
| **★ Kang Table 1, the MFT and FFT rows across 9 nitrogen co-substrates (Tier A)** | **★ HOLD-OUT** | **the sharpest new test the sulfur branch can be given.** The model must predict that FFT/MFT ranges over **0.15–9.23** as the nitrogen co-substrate changes at fixed T, pH and time. A model with a fixed MFT:FFT branch fails by construction — the identical logic that puts the Cerny 2007 Table 5 concentration pair in the declaration's highest-value slot |
| **Kang Table 1, the remaining Tier-A (\*) compounds, 9 systems** | **FIT** | absolute μg L⁻¹ with standard curves; a genuine amino-acid-panel constraint on the sulfur/thiazole/thiophene network at a single well-specified condition. Cut axis: **compound identity** (MFT and FFT out, the rest in) |
| **Kang Table 1, all Tier-B (unstarred) compounds** | **neither — record as `absolute_concentration: false`** | IS-only semi-quantitation with an assumed response factor of 1. Ratio-comparable within a compound; **not scorable as levels.** Same treatment as Shu 1988 and Hemmler 2018 in the existing declaration |
| **Kang Table 1 subtotals, totals and "kinds" rows** | **neither** | every subtotal mixes Tier A and Tier B; and the **"Kinds" counts are demonstrably miscounted** in at least three columns (dossier §4a, §4c) |
| **Kang Fig. 1a/1b class totals (temperature and pH ladders)** | **neither — pending the SI** | figure-digitised class sums of mixed-tier data, and **Fig. 1a's "Sulfur-containing compounds" group is corrupted** (2.5–3.1× too low at 120/140 °C, apparently pasted from the pH experiment). If Tables S4/S5 arrive, use those instead and retire the digitisation |
| Kang K12 measured pH drift (5.5/7/8 → 5.1/4.9/4.7) | **neither — a `condition_correction` fact** | it means none of the "pH" arms is a constant-pH run. Apply as a flag on any Kang row used, do not score |

### 4d. Module 7 (thresholds / matrix correction) and a possible new stability lane — **Kang storage block**

| dataset | **proposed role** | one-line reasoning |
|---|---|---|
| **★ Kang Fig. 2i/2j, solid TTCA and ARP vs a_w (2 compounds × 3 a_w × 7 times = 42 points)** | **★ HOLD-OUT** | **the corpus's only water-activity series on intermediate *stability*.** Bell gives a_w on pigment *formation*; this gives a_w on precursor *survival*, in a different matrix, from a different lab. Holding it out keeps the **agreement or disagreement with Bell's a_w term an out-of-sample fact** — which is exactly the argument the declaration uses for the Zhou/Zhang cross-lab pair in Module 2 |
| **Kang Fig. 2a/2b/2e/2f, solution storage vs T and pH (84 points)** | **FIT** | the temperature and pH axes of the same stability lane, in solution. Cut axis: **phase and stressor** (solution T/pH in, solid a_w out) |
| **Kang K17 storage Eₐ (43.1 / 18.7 kJ mol⁻¹)** | **FIT (as a wide prior only)** | 4–40 °C, three points, no error bars on k, and the TTCA value inherits a bad 4 °C fit (R² 0.615). Usable as a magnitude prior; **not as a target** |
| **Kang Table 2 colour (ΔE 9.36 → 46.62 across a_w)** | **HOLD-OUT (second observable)** | an independent browning probe on the *same* a_w ladder as Fig. 2i/2j, and it scales **differently** from mass loss (K20: ARP loses 1.69× more mass but develops 1.95× more colour). It scores per-pathway chromophore yield, which nothing else in the corpus does |
| **★ Kang Fig. 2g, the pH 5.5 TTCA divergence (A₄₂₀ ↑ while A₂₉₄ ↓)** | **★ HOLD-OUT (negative/sign test)** | the authors read it as **acid-catalysed reversion to N-xylosylamine, not degradation to pigment**. A model that equates intermediate loss with browning predicts the wrong sign on A₂₉₄. One of the very few sign tests available anywhere in the corpus |
| Kang Fig. 1c browning A₄₂₀ across 15 amino acids | **FIT** | 15 levels at one condition, and it is the only amino-acid-identity browning ladder in the corpus. ⚠ **Note that the paper's own claim that added amino acids always raise browning is falsified by its Cys bar (0.072 < control 0.106)** — ingest the numbers, not the sentence |
| Kang §3.4 antioxidant block (K22) | **neither** | "MRPs" is never defined by any heating time, temperature or concentration, so no comparison in that section is reproducible |

### 4e. Cross-paper consistency checks worth registering as scored derived tests

| proposed test | what it scores | why it is a good test |
|---|---|---|
| **★ Bell B9 (a_w plateau above 0.54 at fixed molality) vs Kang K18 (d ln k/d a_w = 1.4–1.8, monotone to a_w 0.843)** | whether the repo's a_w term can be **flat for pigment formation and monotone for precursor loss at the same time** | two papers, two matrices, two observables, one axis. If the repo carries a single a_w multiplier applied to both lanes, this pair falsifies it |
| **★ Nakamura N11 (MGO+Lys : Xyl+Lys = 12.2×) vs the trunk's own MGO node** | whether the trunk's α-dicarbonyl branch carries enough browning flux | a pure ratio, response-factor-immune, and the corpus has no other measurement of relative browning power at matched molarity |
| Nakamura N8 (1-DX 76 % of xylose) vs the repo's pentose deoxyosone split | topology, not level | a fitted split cannot be asked this unless it is held out |
| Kang K4 (FFT/MFT = 0.15–9.23) vs any fixed MFT:FFT ratio in the sulfur branch | whether the branch fractions respond to the amine at all | structurally identical to the Cerny 2007 concentration-pair test already in Module 9 |

---

## 5. THINGS THE ORCHESTRATOR SHOULD DECIDE

1. **Fetch the Kang SI** (Tables S2–S5, Figs S1–S3) or formally record that Kang 2026 cannot
   serve as the sulfur temperature ladder. **This is the single largest open item from this
   wave.**
2. **Audit the corpus for the μ → m text-layer corruption** (§1.1). Any RSC/Arbortext PDF
   ingested via `pdftotext` is at risk of a 1000× unit error.
3. **Decide whether Kang's storage block opens a new lane** ("intermediate stability under
   T / pH / a_w stress"). It does not fit cleanly into Modules 1–9 as written; §4d proposes
   Module 7 by analogy but a Module 10 may be cleaner.
4. **Rule on the Bell Tg/MW confound** (§4a note). If a Tg term is fitted or scored on Bell, the
   declaration should say so explicitly.
5. **Rule on whether figure-digitised data may carry ★ HOLD-OUT status.** Nakamura is
   *entirely* figure-based (no tables at all), and four of my highest-value proposals are
   digitisations. Precedent exists — the declaration already grants HOLD-OUT to the
   `digitised_from_figure` Hofmann 1996b T × t surface — but Nakamura would be the first paper
   where *every* row is digitised.
6. **Note two filename/year mismatches** for the bibliography: `bell1995` → 1996;
   `nakamura2020` → 2021 (issue of record).

---

*Drafted 2026-08-28 by Wave K4c. Values in §2 are transcribed or derived in the three companion
dossiers, where every anchor, calibration and error estimate is given per table and per figure
panel. §4 is a proposal for the orchestrator and amends nothing.*
