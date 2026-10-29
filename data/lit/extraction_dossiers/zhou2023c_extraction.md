# Zhou T., Xia, Cui, Hayat, Zhang & Ho 2023 (10.1021/acs.jafc.3c04424) — EXTRACTION
### *"Promotion or Inhibition Effects of Exogenous Glutathione-Degraded Amino Acids on the Formation of 2,3-Butanedione and Pyrazines via Varied Pathways of Interaction with α-Dicarbonyl Compounds Derived from N-(1-Deoxy-D-xylulos-1-yl)-alanine"*
### The only paper in the corpus with a measured glutathione time course and its degradation products. GSH 20 mmol/L heated alone or with 20 mmol/L Ala-xylose Amadori compound (ARP), water, initial pH 7.5, 120 °C, 0–120 min: GSH, pGlu, Glu, Cys-Gly, Cys, Gly, ARP, Ala (UPLC-MS/MS, external standards); 3-/1-deoxypentosone, GO, MGO (OPD/HPLC); 16 volatiles incl. **MFT and FFT in µg/L at 30/60/90/120 min** (standard-calibrated), with an ARP + Cys control. Read 2026-10-27.

**Source on disk:** `data/articles/zhou2023c.pdf` (publisher's PDF, 10 pp.; stamped "Downloaded … by
DOKUZ EYLUL UNIV user on 07 September 2026"). Text from `pdftotext -layout`; Figs. 1, 3, 4, 5 read from
300 dpi renders (pp. 3, 5, 7). **Fig. 1 is a heat map with every cell printed as a number** — transcribed
[M], not digitised. Supporting Information (Table S1 MRM parameters; Fig. S1 GSH–3-deoxypentosone adduct
MS/MS) **not on disk**.

**Duplicate check (requested):** `zhou2023_extraction.md` is **Zhou et al. 2023, 10.1021/acs.jafc.2c08360**
(JAFC 71(5):2472, "Cysteine-induced pH-dependent formation of thiols … N-(1-deoxy-D-xylulos-1-yl)-
alanine", source `data/articles/Zhou2023b.pdf` + `zhou2023b_supporting.pdf`). There is **no
`zhou2023b_extraction.md`**. No dossier carries 10.1021/acs.jafc.3c04424 (grep of
`data/lit/extraction_dossiers/`); this is a new paper, not a duplicate. The 2c08360 paper is this one's
ref. 19 (same ARP preparation) and the ARP + Cys sister study is its ref. 15 (Zhou et al. 2022, JAFC
70:15202).

**Provenance:** **[M]** printed number · **[Q]** qualitative statement · **[F]** read off a figure by me
(read error stated) · **[Z]** derived by me, never printed · **[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Tong Zhou, Xue Xia, Heping Cui, Khizar Hayat, Xiaoming Zhang\*, Chi-Tang Ho\* |
| Title | as above |
| Venue | *J. Agric. Food Chem.* **2023**, **71** (39), **14312–14321** |
| DOI | 10.1021/acs.jafc.3c04424 (printed) |
| Affiliation | Jiangnan University (Wuxi); Miami University; Rutgers |
| Dates | received 29 Jun 2023; revised 30 Aug 2023; accepted 13 Sep 2023; published 22 Sep 2023 |

---

## 1. THE ONE-PARAGRAPH ANSWER

**Glutathione's fate.** Heated alone (20 mM, water, pH 7.5 start, 120 °C), GSH falls 20 → 7.45 mM in
2 h (**62.75 % degraded**, printed), close to first order at **k ≈ 8.2 × 10⁻³ min⁻¹** [Z]. The products
measured are **pyroglutamic acid 4.7 mM (the largest), Cys-Gly 2.3 mM, Glu 0.8 mM, Gly 0.85 mM and free
Cys ≈ 0.05 mM** — γ-glutamyl cleavage with Glu cyclising to pGlu. Together they recover **44 % of the lost
γ-Glu and 19 % of the lost Cys** [Z]; the rest (incl. any GSSG, H₂S and Cys-Gly cyclisation) is not
measured. With the ARP present, GSH falls **twice as fast (20 → 2.9 mM, k ≈ 1.6 × 10⁻² min⁻¹)**, pGlu
still accumulates (3.6 mM), while Cys-Gly and Cys are almost absent (consumed) and the ARP's own
breakdown is slowed (k 0.022 → 0.014 min⁻¹). Both GSH loss rates are **7–14 × the Zheng & Ho H₂S-release
constant the model uses** (extrapolated to 120 °C, pH 7.5). **Thiols.** ARP + GSH gives **MFT 59, 64,
57, 8 µg/L and FFT 124, 383, 270, 116 µg/L at 30/60/90/120 min** — peak ≈ 2.8 × 10⁻⁵ (MFT) and
1.7 × 10⁻⁴ (FFT) mol/mol GSH [Z]; the ARP alone gives none; **ARP + free Cys at 60 min gives MFT
1150 µg/L (18 × GSH) and FFT 511 µg/L (1.3 × GSH)**. GSH lowers 3-deoxypentosone (1.14 → 0.41 mM at 2 h)
and glyoxal (1.41 → 0.55 mM) but not 1-deoxypentosone at 30 min.

---

## 2. METHODS [M] (pp. 14313–14314)

| item | as printed |
|---|---|
| ARP | N-(1-deoxy-D-xylulos-1-yl)-alanine, prepared: xylose 0.4 mol/L + Ala 0.2 mol/L in 100 mL water, pH 7.5 ± 0.1 (6 M NaOH), 80 °C 60 min, then rotary vacuum dehydration 80 °C 20 min; Dowex 50WX8 (H⁺), NH₄OH 0.2 mol/L; HPLC-ELSD fractionation; freeze-dried; NMR as in ref. 19 |
| ARP ± GSH | **ARP 20 mmol/L ± GSH 20 mmol/L**, aqueous (no buffer), **initial pH 7.5 with 3 mol/L NaOH**; **pressure-resistant bottles, oil bath 120 °C, 30, 60, 90, 120 min**; ice-water quench |
| ARP + Cys control | ARP + Cys (1 : 1, 20 mmol/L each implied), 120 °C, 60 min, pH 7.5 |
| GSH alone | **20 mmol/L GSH, pH 7.5, 120 °C, 30/60/90/120 min** |
| "Degraded GSH" series | GSH pre-heated alone 0/30/60/90/120 min (120 °C, pH 7.5); then **+ ARP 20 mmol/L**, initial pH 7.5, **120 °C, 60 min** |
| ARP + single amino acid | ARP 20 mmol/L + Cys, Gly, pGlu or Glu (1 : 1), pH 7.5, 120 °C, 60 min |
| Headspace / atmosphere | not stated (bottle volume, fill, air/N₂) [!] |
| pH during heating | not measured; text: ARP degradation releases acids, "rapid reduction of the pH" [Q] [!] |
| α-Dicarbonyls | OPD derivatisation → quinoxalines; HPLC-DAD 315 nm, SunFire C18; 3-DP, 1-DP, GO, MGO (method of ref. 15) |
| Amino acids / peptides | ACQUITY UPLC-TQD, BEH Amide, ESI⁺ MRM; **external standards ARP, Ala, GSH, Gly, pGlu, Cys, Glu, "Gly−Cys"**. Eight calibration lines printed (in the stated order): y = 4.05×10⁶x − 5164.14 (R² 0.9996, LOD 0.15, LOQ 0.59 µmol/L); y = 1.55×10⁷x + 2.66×10⁶x [sic] (0.9996; 0.38; 1.12); y = 5.97×10⁶x + 2.40×10⁴ (0.9983; 0.03; 0.07); y = 3.96×10⁶x + 2396.22 (0.9984; 1.33; 5.62); y = 2.78×10⁷x + 1.87×10⁵ (0.9916; 0.05; 0.19); y = 7.47×10⁶x − 2296.61 (0.9983; 0.18; 0.64); y = 2.54×10⁷x + 2.17×10⁴ (0.9986; 0.02; 0.07); y = 3.83×10⁶x − 1.48×10⁴ (0.9928; 0.04; 0.17); x in mg/mL, y peak area |
| Volatiles | HS-SPME: 3 mL sample + 2 mL saturated NaCl + IS 1,2-dichlorobenzene (5 µL of 0.0018 µg/µL); DVB/CAR/PDMS 60 °C 30 min; 7890B/5977B, DB-WAX; **quantified with calibration curves of authentic standards in water** (MFT, FFT, 2,3-butanedione, the pyrazines, thiophenes and thiazoles listed as purchased standards) |
| Replicates | triplicate, mean ± SD; Duncan p < 0.05 |

[!] Naming: "Gly−Cys" in Methods and the Fig. 3c legend, "Cys−Gly" in Fig. 3b and the text; the
γ-glutamyl-cleavage product is Cys-Gly. The second calibration line prints "+ 2.66 × 10⁶x" (a typo).

---

## 3. RESULTS

### 3.1 Fig. 1 (p. 14314): volatiles heat map, all cells printed [M]

Units are not printed on the figure; they are **µg/L** [Z: the ARP+GSH-60 and ARP+Cys-60 cells reproduce
Table 1 (µg/L) exactly — 6.87×10¹ = 68.7, 8.33 = 8.33, 9.25×10⁻¹ = 0.925, 4.51 = 4.51]. 0 = not detected.

| compound | ARP 30 | ARP+GSH 30 | ARP 60 | ARP+GSH 60 | ARP+Cys 60 | ARP 90 | ARP+GSH 90 | ARP 120 | ARP+GSH 120 |
|---|---|---|---|---|---|---|---|---|---|
| 2,3-butanedione | 69.2 | 105 | 23.0 | 68.7 | 0.925 | 16.7 | 94.5 | 22.4 | 79.5 |
| furfural | 341 | 113 | 619 | 397 | 230 | 874 | 575 | 1130 | 716 |
| 2-propionylfuran | 1.08 | 0.0230 | 1.87 | 1.24 | 0 | 4.46 | 4.63 | 6.93 | 10.8 |
| 2-acetylfuran | 0.439 | 0 | 15.8 | 5.02 | 0.232 | 34.7 | 14.9 | 61.7 | 27.3 |
| pyrazine | 0 | 8.92 | 0 | 8.33 | 11.0 | 0 | 7.25 | 0 | 12.3 |
| methylpyrazine | 0 | 8.26 | 0 | 13.1 | 4.51 | 0 | 8.54 | 0 | 10.5 |
| 2,5-dimethylpyrazine | 0 | 0.226 | 0 | 3.75 | 0 | 0 | 5.39 | 0 | 7.97 |
| 2-methylthiophene | 0 | 2.07 | 0 | 6.40 | 33.4 | 0 | 6.73 | 0 | 4.40 |
| 2-thiophenecarboxaldehyde | 0 | 0.265 | 0 | 1.33 | 51.9 | 0 | 4.25 | 0 | 7.70 |
| 5-methyl-2-thiophenecarboxaldehyde | 0 | 0 | 0 | 0 | 8.54 | 0 | 0 | 0 | 0 |
| thieno[3,2-b]thiophene | 0 | 0 | 0 | 0 | 3.46 | 0 | 0 | 0 | 0 |
| 2-methylthiazole | 0 | 0 | 0 | 0 | 0.924 | 0 | 0 | 0 | 0 |
| thiazole | 0 | 0 | 0 | 0 | 23.2 | 0 | 0 | 0 | 3.23 |
| 2-acetylthiazole | 0 | 0 | 0 | 9.61 | 18.6 | 0 | 11.8 | 0 | 11.0 |
| **2-methyl-3-furanthiol** | 0 | **59.1** | 0 | **64.1** | **1150** | 0 | **57.1** | 0 | **8.10** |
| **2-furfurylthiol** | 0 | **124** | 0 | **383** | **511** | 0 | **270** | 0 | **116** |

No SD printed in the heat map [!]. **Molar [Z]** (MW 114.17; 20 mmol/L GSH):

| | 30 min | 60 min | 90 min | 120 min | ARP+Cys 60 |
|---|---|---|---|---|---|
| MFT µmol/L | 0.518 | 0.561 | 0.500 | 0.071 | 10.07 |
| MFT mol/mol GSH (or Cys) | 2.6 × 10⁻⁵ | **2.8 × 10⁻⁵** | 2.5 × 10⁻⁵ | 3.5 × 10⁻⁶ | 5.0 × 10⁻⁴ |
| FFT µmol/L | 1.09 | 3.35 | 2.37 | 1.02 | 4.48 |
| FFT mol/mol GSH (or Cys) | 5.4 × 10⁻⁵ | **1.7 × 10⁻⁴** | 1.2 × 10⁻⁴ | 5.1 × 10⁻⁵ | 2.2 × 10⁻⁴ |

MFT plateaus 30–90 min then **falls 7-fold between 90 and 120 min**; FFT peaks at 60 min and falls 3.3-fold
by 120 min — both are net formation/loss in a closed bottle.

### 3.2 Table 1 (p. 14318): ARP + exogenous GSH or single amino acid, 120 °C, 60 min, pH 7.5 start, µg/L, mean ± SD (n = 3) [M]

| compound | RI | KI | ID | ARP + GSH | ARP + Glu | ARP + Cys | ARP + Gly | ARP + pGlu |
|---|---|---|---|---|---|---|---|---|
| 2,3-butanedione | 979 | 977 | MS, RI, S | 68.7 ± 2.09b | 74.1 ± 2.30b | 0.925 ± 0.739a | 112 ± 8.38d | 85.5 ± 4.66c |
| pyrazine | 1205 | 1209 | MS, RI, S | 8.33 ± 1.06b | 8.33 ± 1.00b | 11.0 ± 1.28c | 5.97 ± 1.25a | ND |
| methylpyrazine | 1265 | 1264 | MS, RI, S | 13.1 ± 2.72b | 29.5 ± 2.30c | 4.51 ± 0.494a | ND | ND |
| 2,5-dimethylpyrazine | 1321 | 1323 | MS, RI, S | 3.75 ± 0.738a | 15.1 ± 2.48b | ND | 1.16 ± 0.164a | ND |
| 2,3-dimethylpyrazine | 1343 | 1346 | MS, RI, S | ND | 7.63 ± 1.07 | ND | ND | ND |
| trimethylpyrazine | 1400 | 1408 | MS, RI, S | ND | 12.5 ± 3.30 | ND | ND | ND |
| total pyrazines | | | | 25.2 ± 4.50c | 72.8 ± 7.26d | 15.5 ± 0.803b | 7.13 ± 1.40a | ND |

### 3.3 Fig. 3 (p. 14316): GSH, its products, ARP and Ala, mmol/L [F]

**(b) GSH alone, 20 mmol/L, pH 7.5, 120 °C.** GSH on a broken axis (8–20); read error ±0.5 mM for GSH,
±0.1 mM for the rest. 120 min GSH from the printed 62.75 % [M].

| t (min) | GSH | pGlu ("PCA") | Cys-Gly | Glu | Gly | Cys |
|---|---|---|---|---|---|---|
| 0 | 20 | 0 | 0 | 0 | 0 | 0 |
| 30 | 15.8 | 1.95 | 0.75 | 0.65 | 1.12 | ~0.05 |
| 60 | 12.5 | 3.45 | 1.50 | 0.78 | 0.98 | ~0.06 |
| 90 | 9.2 | 4.42 | 1.85 | 0.75 | 0.88 | ~0.07 |
| 120 | **7.45** [Z from 62.75 % M] | 4.72 | 2.32 | 0.83 | 0.85 | ~0.05 |

**(c, d) ARP + GSH, 20 + 20 mmol/L, pH 7.5, 120 °C.** Read error ±0.3 mM (GSH), ±0.05 mM (rest).

| t (min) | GSH (d) | pGlu | Glu | Gly | Cys | Cys-Gly |
|---|---|---|---|---|---|---|
| 0 | 20 | 0 | 0 | 0 | 0 | 0 |
| 30 | 14.1 | 1.50 | 0.77 | 0.52 | ~0.03 | ~0.02 |
| 60 | 8.6 | 2.57 | 0.78 | 0.35 | ~0.03 | ~0.02 |
| 90 | 5.7 | 3.35 | 0.76 | 0.17 | ~0.03 | ~0.02 |
| 120 | 2.9 | 3.60 | 0.72 | 0.07 | ~0.03 | ~0.02 |

**(a) ARP and released Ala, mmol/L** (±0.3 mM)

| t (min) | ARP (ARP alone) | ARP (ARP+GSH) | Ala (ARP alone) | Ala (ARP+GSH) |
|---|---|---|---|---|
| 0 | 20 | 20 | 0 | 0 |
| 30 | 10.0 | 12.9 | 7.2 | 6.8 |
| 60 | 5.4 | 9.25 | 7.8 | 8.0 |
| 90 | 3.2 | 5.65 | 7.9 | 8.4 |
| 120 | 1.45 | 3.6 | 8.1 | 7.95 |

### 3.4 Fig. 4 (p. 14316): α-dicarbonyls, mmol/L, ARP 20 mM ± GSH 20 mM, 120 °C [F]

Read error ±0.02 mM (3-DP, 1-DP, GO), ±0.003 mM (MGO).

| t (min) | 3-DP −GSH | 3-DP +GSH | 1-DP −GSH | 1-DP +GSH | GO −GSH | GO +GSH | MGO −GSH | MGO +GSH |
|---|---|---|---|---|---|---|---|---|
| 30 | 0.65 | 0.60 | 0.43 | 0.43 | 0.70 | 0.25 | 0.068 | 0.071 |
| 60 | 0.92 | 0.59 | 0.31 | 0.28 | 1.07 | 0.32 | 0.083 | 0.076 |
| 90 | 1.06 | 0.48 | 0.23 | 0.175 | 1.31 | 0.40 | 0.101 | 0.073 |
| 120 | 1.14 | 0.41 | 0.19 | 0.125 | 1.41 | 0.55 | 0.111 | 0.062 |

### 3.5 Fig. 5 (p. 14318) [F]

**(a) ARP + pre-degraded GSH, 120 °C 60 min, µg/L** (±1.5 µg/L):

| GSH pre-heating (min) | 0 | 30 | 60 | 90 | 120 |
|---|---|---|---|---|---|
| 2,3-butanedione | 68.5 | 77 | 89.5 | 91.5 | 100.5 |
| total pyrazines | 25 | 34 | 42.5 | 49.5 | 58 |

**(b) GO and MGO after ARP + X, 120 °C 60 min, mmol/L** (±0.015):

| X | Glu | Gly | GSH | pGlu | Cys |
|---|---|---|---|---|---|
| GO | 0.50 | 0.725 | 0.30 | 0.94 | 0.195 |
| MGO | 0.18 | 0.24 | 0.075 | 0.28 | ~0.005 |

The 0-min bar in 5a (68.5) and the GSH bar in 5b (GO 0.30) agree with Fig. 1 (68.7) and Fig. 4c (0.32 at
60 min) [Z].

### 3.6 Derived quantities [Z]

**GSH apparent first-order loss, 120 °C, pH 7.5 start:**

| system | 0–30 | 30–60 | 60–90 | 90–120 | 0–120 min |
|---|---|---|---|---|---|
| GSH alone (min⁻¹) | 7.9 × 10⁻³ | 7.8 × 10⁻³ | 1.0 × 10⁻² | 7.0 × 10⁻³ | **8.2 × 10⁻³** |
| GSH + ARP (min⁻¹) | 1.2 × 10⁻² | 1.6 × 10⁻² | 1.4 × 10⁻² | 2.3 × 10⁻² | **1.6 × 10⁻²** |

ARP loss: alone 2.2 × 10⁻² min⁻¹ (0–120 min), with GSH 1.4 × 10⁻² min⁻¹.

**Against the model's constant.** An unweighted 4-point Arrhenius fit of Zheng & Ho 1994 Table III
(`zheng1994_extraction.md` §4; H₂S release, 0.1 M GSH, 80–110 °C) gives Ea 96 kJ/mol at pH 7 and 84 kJ/mol
at pH 9, and at 120 °C k = 9.2 × 10⁻⁴ (pH 7), 2.4 × 10⁻³ (pH 9), **1.2 × 10⁻³ min⁻¹ at pH 7.5**
(log-interpolated). Measured GSH loss is **~7 × that (alone) and ~14 × (with ARP)** — a 10 K
extrapolation, smaller than du2025's 35 K. (The engine's own refit may differ slightly; this is my
re-derivation, not a model run.)

**Mass balance, GSH alone at 120 min** (12.55 mM lost):
- γ-Glu side: pGlu 4.72 + Glu 0.83 = 5.55 mM = **44 %** of lost GSH.
- Cys side: Cys-Gly 2.32 + Cys 0.05 = 2.37 mM = **19 %**. Free Gly 0.85 mM implies ≥ 0.85 mM Cys-Gly
  was hydrolysed, so Cys released ≥ 0.85 mM but free Cys is ~0.05 mM: **≥ 0.8 mM Cys consumed** (to H₂S
  per the authors; GSSG/disulfides and cyclo(Cys-Gly) not measured).
- Unaccounted: ~56 % of the γ-Glu and ~81 % of the Cys moiety. A Zheng & Ho H₂S channel alone
  (1.2 × 10⁻³ min⁻¹ × ∫GSH dt ≈ 1550 mM·min) would release only ~1.8 mM H₂S; **most of the missing
  sulfur is in a pool the paper does not measure (likely GSSG/mixed disulfides in an air headspace)** [!].

**With ARP at 120 min** (17.1 mM lost): pGlu + Glu = 4.32 mM (25 %); Cys-Gly + Cys ≈ 0.05 mM (0.3 %).
The ARP adds ~4.5 mM extra GSH loss over GSH alone, and Cys-Gly/Cys are scavenged (deoxypentosone
adducts; an MS/MS-only GSH–3-DP adduct in Fig. S1).

---

## 4. WHAT THIS CAN AND CANNOT ANCHOR IN THE MODEL

**Model context.** `GSH` → H₂S + FRAG_C(10) + FRAG_N(3) is its only reaction (`r_gsh_h2s`,
`k_gsh_h2s` = Zheng & Ho 1994, pH-resolved, flagged `h2s_release_only`; `parameters_sulfur.py`
l. 1934, `sulfur.py` l. 246, 1175). Pentose ARPs, 3-/1-deoxyosones, GO/MGO, MFT and FFT are in the
sulfur/trunk lanes for free cysteine. Nothing here was run against the model.

**(1) GSH total loss rate, 120 °C, pH 7.5 start — usable as a measured constant with caveats.**
**k_GSH,loss ≈ 8.2 × 10⁻³ min⁻¹ (alone, near first order, 4 time points, ±10 % from the read)**, and
≈ 1.6 × 10⁻² min⁻¹ with 20 mM pentose ARP. This is the first direct measurement in the corpus of GSH
disappearance, with external-standard UPLC-MS/MS. It **contradicts the model's implicit assumption that
GSH disappears only via H₂S release**: the H₂S constant is ~7–14 × smaller. Caveats: one temperature, one
pH (start only, unbuffered), no atmosphere control, n = 3 but SD visible only as bars.

**(2) GSH's product branching (usable, semi-quantitative).** The measured fate is **γ-glutamyl cleavage
→ pGlu (dominant identified product) + Cys-Gly, then Cys-Gly → Cys + Gly**, with free Cys held near
zero. A minimal, evidence-backed GSH extension: GSH → Cys-Gly + pGlu (k ≈ 3–4 × 10⁻³ min⁻¹ for the
identified part, i.e. 44 % of 8.2 × 10⁻³), Cys-Gly → Cys + Gly (slower; Gly 0.85 mM vs Cys-Gly 2.3 mM at
2 h), and Cys entering the existing cysteine chemistry. The unidentified ~50 % (probably oxidation to
GSSG) must stay an explicit sink, not be routed to H₂S. Glycation is a small branch here (the GSH–3-DP
adduct is shown by MS only; compare `wang2020_extraction.md`, `sun2019b_extraction.md`: GSH-ARP ≤ 0.3 %
at 140 °C, ≤ 9 % at 90–100 °C).

**(3) MFT/FFT from GSH vs Cys — usable direction/ratio test.** Same ARP, same pot, 60 min:
**MFT(Cys)/MFT(GSH) = 18; FFT(Cys)/FFT(GSH) = 1.3**. GSH is a far weaker MFT source than free cysteine
but almost as good an FFT source ([Z, my interpretation, not the authors'] GSH's slow Cys release starves the H₂S + 1-deoxyosone/norfuraneol → MFT
route more than the furfural → FFT route). Absolute: MFT ≈ 2.8 × 10⁻⁵, FFT ≈ 1.7 × 10⁻⁴ mol/mol GSH at
the 60-min peak, both decaying afterwards (closed bottle). A model with GSH → Cys-Gly → Cys at the rates
in (1)–(2) should reproduce the 18× and ~1× ratios with the existing Cys chemistry; a model with only
`r_gsh_h2s` routes all GSH sulfur through H₂S and has no way to make that asymmetry except through
H₂S supply. Note the ARP here is an **Ala-xylose** ARP (a pentose-ARP precursor, not free xylose).

**(4) Deoxyosone and dicarbonyl suppression by GSH (usable as direction + magnitude).** At 120 min GSH
lowers 3-DP 2.8 ×, GO 2.6 ×, MGO 1.8 ×, 1-DP 1.5 ×; at 30 min only GO differs. Consistent with thiol
trapping of 3-DP (MS adduct) and with the model's Cys + deoxyosone chemistry, but the paper offers no
rate constants.

**(5) No lipid in this paper.** No lipid × sulfur numbers.

**(6) Not usable.** pH-resolved anything (pH 7.5 start only, drift unmeasured); the 2,3-butanedione /
pyrazine mechanism claims (Gly + MGO → 2,3-butanedione) as rates; the heat map's lack of SD.

**Caveats, collected.**
1. Unbuffered water, NaOH to pH 7.5 at start; ARP breakdown acidifies (authors).
2. Bottle size, fill and atmosphere unstated; GSSG not measured → sulfur balance open.
3. GSH at 30/60/90 min on a compressed broken axis (±0.5 mM); 120 min from the printed 62.75 %.
4. Fig. 1 units inferred (µg/L) from agreement with Table 1.
5. "Gly-Cys" vs "Cys-Gly" naming inconsistency; one calibration line typo.
