# Liu, Yu, Zhou, Xu, Hayat, Zhang & Ho 2022 — EXTRACTION
### *"Formation Priority of Pyrazines and 2-Acetylthiazole Dependent on the Added Cysteine and Fragments of Deoxyosones during the Thermal Process of the Glycine–Ribose Amadori Compound"*
### 20 mM glycine–ribose Amadori compound (GR-ARP) ± 10/20/40 mM cysteine, and fed 20 mM methylglyoxal, glyoxal or both ± 10–40 mM cysteine, pH 7.5 (initial), 120 °C, 15–60 min. Time courses of cysteine, MGO and GO in mmol/L, and of 2-acetylthiazole, pyrazine and methylpyrazine in µg/L. The dicarbonyls and cysteine vanish in ≤ 15–30 min while 2-acetylthiazole plateaus at ~4.5 µM: the product captures ~0.02 % of the cysteine + MGO consumed. Read 2026-10-27.

**Full read of `data/articles/liu2022.pdf`** (9 pp., born-digital ACS PDF). Prose and Tables 1–2 from
the `pdftotext -layout` text layer. **All kinetic data are in Figures 2–5 and were read from 200-dpi
renders of journal pp. 11646, 11648, 11649, 11650**; read error ≈ ±2 % of each axis span (±20 µg/L on
the 0–1000 axes, ±0.4 mM on the 0–20 mM axes, ±0.05 on the 0–5 µg/L axis). No supplementary
information. **No existing `liu2022_extraction.md` in the repository** (checked: `liu2021`, `liu2023b`,
`liu2025` are other papers), so this dossier takes the plain name. Sister paper, same Jiangnan/Rutgers
group and design: `zhou2023_extraction.md` (xylose–alanine ARP + cysteine, pH 5–8).

**Provenance:** **[M]** printed number · **[F]** read off a figure (±2 % of axis) · **[Q]**
qualitative · **[Z]** derived by me · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Meichen Liu, Jingyang Yu, Tong Zhou, Huaneng Xu, Khizar Hayat, Xiaoming Zhang,* Chi-Tang Ho* |
| Title | as above |
| Venue | *J. Agric. Food Chem.* **2022, 70, 11643–11651** (Crossref: issue 37) |
| DOI | **10.1021/acs.jafc.2c04874** (printed on every page footer) |
| Affiliations | State Key Lab of Food Science and Technology, Jiangnan University, Wuxi (all but two); Rutgers, New Brunswick (Ho); Miami University, Oxford OH (Hayat) |
| Dates | received 9 July 2022; revised 28 Aug 2022; accepted 30 Aug 2022; published 7 Sep 2022 |
| Funding | NSFC 32172330; JUFSTR20180204 |

---

## 1. THE ONE-PARAGRAPH ANSWER

Purified N-(1-deoxy-D-ribulos-1-yl)-glycine (GR-ARP, 20 mM), heated alone at pH 7.5 / 120 °C / 60 min,
gives only furans (furfural 1503 µg/L) and releases **glyoxal (GO) to 1.57 mM and methylglyoxal (MGO)
to 0.44 mM** by 60 min [F]. Adding cysteine (10/20/40 mM) produces **2-acetylthiazole ≈ 850 / 300 /
120 µg/L** (falling with dose) and **pyrazine ≈ 1.5 / 6.3 / 8.7 µg/L** (rising with dose) [F]. In fed
dicarbonyl models (20 mM MGO or GO ± 10/20/40 mM Cys, same pH/T), **cysteine is gone within 15 min in
most arms** (40 mM Cys + MGO: within 30 min), MGO falls to ≤ 3 mM by 15 min at ≥ 20 mM Cys, GO to
< 0.5 mM, while **2-acetylthiazole plateaus by 15 min at 350–580 µg/L (2.8–4.6 µM)** and pyrazine
stays < 4 µg/L. The authors read this as Cys reacting "preferentially" with MGO to give
2-acetylthiazole, then with GO via carbonylcysteimine → cysteamine + recovered GO. **No H2S, no
thiols other than MFT (in one table) and no pH after heating are reported.**

---

## 2. METHODS [M] (p. 11645)

| item | as printed |
|---|---|
| GR-ARP synthesis | equimolar Gly + ribose (**0.2 mol/L**), pH 7.5 (6 M NaOH), 100 mL, rotary vacuum dehydration **80 °C / 20 min**, < 25 mbar within 5 min; purified on Dowex 50WX4 (H+), eluted with water then 0.18 M NH4OH; HPLC-ELSD (BEH Amide) fractions freeze-dried (purity not stated) [!] |
| ARP ± Cys | Cys **10, 20, 40 mmol/L** added to **20 mmol/L GR-ARP**; "pH levels ... adjusted to **7.5** before heat treatment"; **oil bath 120 °C for 15, 30, 45, 60 min** "separately", in "temperature- and pressure-resistant bottles"; ice-water quench. Volume, headspace, buffer: **not stated** (no buffer named) [!] |
| Cys + dicarbonyl | Cys **10, 20, 40 mmol/L** added to MGO, GO, or MGO/GO aqueous solutions (**20 mmol/L** each); pH 7.5; 120 °C; 15–60 min; same bottles |
| α-Dicarbonyls | o-phenylenediamine derivatisation (20 g/L OPD + 11 mM DTPA in pH 7.4 phosphate), 500 µL + 500 µL, 4 h dark, 3,4-hexanedione internal standard (50 µL of 30 mM); HPLC-DAD 315 nm, SunFire C18 |
| Cys | HPLC-ELSD, external standard, **y = −40.196x² + 97379x + 58661, R² = 0.9938** (x mmol/L, y peak area) |
| Volatiles | HS-SPME: 3 mL sample + 1.2 g NaCl in 20 mL vial, 55 °C 10 min equilibration, DVB/CAR/PDMS 75 µm 30 min; GC-MS DB-WAX 30 m × 0.25 mm × 0.25 µm, 40 °C (3 min) → 80 (5 °C/min) → 160 (10 °C/min, 0.5 min) → 175 (2 °C/min) → 230 °C (10 °C/min, 7 min); EI 70 eV, 35–450 m/z |
| Quantification | **external standard** for pyrazine, methylpyrazine, 2-acetylthiazole: "ln y = 2.09 × 10⁻⁶x (R² 0.9993)", "ln y = 7.31 × 10⁻⁷x (R² 0.994)", "ln y = 2.19 × 10⁻⁸x (R² 0.9903)", x peak area, y concentration [!] |
| Replication | "All experiments ... repeated three times"; Duncan's test, p < 0.05 |

---

## 3. RESULTS

### 3.1 Table 1 — GR-ARP alone, 120 °C (µg/L) [M] (p. 11646)

| compound | KI^B | RI^C | ID^D | concentration |
|---|---|---|---|---|
| 2-methylfuran | 898 | 911 | MS, RI, S | 3.1 ± 0.8 d |
| 2,2′-bifuran | 1614 | 1615 | MS, RI, S | 4.7 ± 0.6 d |
| 2-furfural | 1444 | 1460 | MS, RI, S | **1503.0 ± 23.0 a** |
| 3-acetyl-2,5-dimethylfuran | 1597 | 1567 | MS, RI, S | 10.4 ± 2.3 c |
| 2(5H)-furanone | 1763 | 1767 | MS, RI, S | 2.0 ± 0.4 e |
| 4-hydroxy-5-methyl-3(2H)-furanone | 2079 | 2124 | MS, RI, S | 16.4 ± 1.2 b |

^A "All results were determined by the SPME method and quantitatively measured by an external
standard. ... means ± standard deviation. ... Duncan's multiple comparison test (n = 3)." ^B LRI on
DB-WAX. ^C literature RI (flavornet, NIST webbook). ^D identification methods. Time not stated in the
table; Figure 2 caption says 60 min.

### 3.2 Table 2 — GR-ARP + Cys, 120 °C (µg/L) [M] (p. 11647)

Text: "At the same concentration of 20 mmol/L, GR−ARP reacted with Cys at 120 °C for 1 h, and the
initial pH was 7.5."

| compound | KI | RI | ID | concentration |
|---|---|---|---|---|
| 2-methyl-3-furanthiol | 1307 | 1304 | MS, RI, S | **78.4 ± 2.5 b** |
| 2-thiophenecarboxaldehyde | 1695 | 1678 | MS, RI, S | 40.1 ± 2.0 c |
| 2-acetylthiazole | 1657 | 1653 | MS, RI, S | **302.4 ± 6.4 a** |
| pyrazine | 1209 | 1213 | MS, RI, S | 1.5 ± 0.2 e |
| methylpyrazine | 1256 | 1244 | MS, RI, S | 3.0 ± 0.4 d |

### 3.3 Figure 2 — GR-ARP ± Cys, pH 7.5, 120 °C, 60 min [F] (p. 11646)

(a) total volatiles: **20 GR-ARP ≈ 970 µg/L; "20 Gly + 20 Rib" ≈ 325 µg/L** (caption says "the model
of Gly/Rib/Cys ... all 20 mmol/L"; x-axis label says Gly + Rib only [!]).

(b), (c) at 60 min:

| arm | pyrazine (µg/L) | methylpyrazine (µg/L) | 2-acetylthiazole (µg/L) |
|---|---|---|---|
| 20 GR-ARP | 0 | 0 | 0 |
| + 10 Cys | 1.5 | 3.0 | **850** |
| + 20 Cys | 6.3 | 2.7 | **300** |
| + 40 Cys | 8.7 | 2.2 | **120** |

(d) α-dicarbonyls from GR-ARP (one series each; with or without Cys not stated — read as without) [!]:

| min | 0 | 15 | 30 | 45 | 60 |
|---|---|---|---|---|---|
| GO (mmol/L) | 0 | 0.18 | 0.55 | 0.82 | **1.57** |
| MGO (mmol/L) | 0 | 0.11 | 0.13 | 0.15 | **0.44** |

### 3.4 Figure 3 — 20 mM MGO + Cys [F] (p. 11648)

| quantity | Cys dose | 0 | 15 | 30 | 45 | 60 min |
|---|---|---|---|---|---|---|
| 2-acetylthiazole (µg/L) | 10 mM | 0 | 352 | 390 | 398 | 366 |
| | 20 mM | 0 | **566** | 561 | 583 | 577 |
| | 40 mM | 0 | 492 | 480 | 506 | 489 |
| MGO (mmol/L) | 10 mM | 20 | 7.1 | 3.0 | 1.3 | 0.4 |
| | 20 mM | 20 | 2.9 | 0.5 | 0.1 | 0.1 |
| | 40 mM | 20 | 1.6 | 0.2 | 0.1 | 0.1 |
| Cys (mmol/L) | 10 mM | 10 | 0 | 0 | 0 | 0 |
| | 20 mM | 20 | 0 | 0 | 0 | 0 |
| | 40 mM | 40 | **26.7** | **0** | 0 | 0 |

### 3.5 Figure 4 — 20 mM GO + Cys [F] (p. 11649)

| quantity | Cys dose | 0 | 15 | 30 | 45 | 60 min |
|---|---|---|---|---|---|---|
| pyrazine (µg/L) | 10 mM | 0 | 0.30 | 0.35 | 0.75 | 0.43 |
| | 20 mM | 0 | 0.80 | 0.85 | 1.03 | 1.15 |
| | 40 mM | 0 | 2.62 | 2.35 | 3.58 | 3.35 |
| GO (mmol/L) | 10 mM | 20 | 4.7 | 2.8 | 1.9 | 1.5 |
| | 20 mM | 20 | ~0.3 | ~0.4 | ~0.4 | ~0.4 |
| | 40 mM | 20 | ~0.3 | ~0.4 | ~0.4 | ~0.4 |
| Cys (mmol/L) | 10 mM | 10 | 0 | 0 | 0 | 0 |
| | 20 mM | 20 | **1.9** | 1.7 | 1.3 | 1.1 |
| | 40 mM | 40 | 0 | 0 | 0 | 0 |

(20 and 40 mM GO curves overlap within read error after 15 min.)

### 3.6 Figure 5 — 20 mM MGO + 20 mM GO + Cys [F] (p. 11650)

y-axes of (a)–(c) are labelled **"mmol/L"**; read as µg/L (Flag 3). Panels (a)–(c) carry **six** doses
(10, 20, 25, 30, 35, 40 mM); Methods name three [!].

2-acetylthiazole (a):

| Cys | 15 | 30 | 45 | 60 min |
|---|---|---|---|---|
| 10 | 240 | 209 | 250 | 246 |
| 20 | 196 | 205 | 216 | 221 |
| 25 | 128 | 112 | 116 | 124 |
| 30 | 105 | 119 | 116 | 122 |
| 35 | 82 | 98 | 108 | 115 |
| 40 | 80 | 84 | 99 | 108 |

pyrazine (b) and methylpyrazine (c):

| Cys | pyrazine 15 / 30 / 45 / 60 | methylpyrazine 15 / 30 / 45 / 60 |
|---|---|---|
| 10 | 0 / 0 / 0 / 0 | 0.52 / 0.24 / 0.87 / 0.97 |
| 20 | 0 / 0 / 0.20 / 0.23 | 0.65 / 0.91 / 1.21 / 1.75 |
| 25 | 0 / 0 / 0.23 / 0.33 | 0.84 / 1.03 / 1.45 / 2.01 |
| 30 | 0.05 / 0.19 / 0.33 / 0.39 | 1.05 / 1.87 / 2.43 / 3.57 |
| 35 | 0.16 / 0.43 / 0.44 / 0.58 | 1.35 / 1.93 / 2.78 / 3.27 |
| 40 | 0.13 / 0.36 / 0.47 / 0.65 | 1.63 / 2.48 / 3.12 / 3.95 |

MGO (d), GO (e), Cys (f), mmol/L (doses 10/20/40 only):

| quantity | Cys | 0 | 15 | 30 | 45 | 60 |
|---|---|---|---|---|---|---|
| MGO | 10 | 20 | 12.1 | 10.2 | 9.0 | 8.5 |
| | 20 | 20 | 5.7 | 3.6 | 3.0 | 1.8 |
| | 40 | 20 | 2.1 | 0.5 | 0.3 | 0.2 |
| GO | 10 | 20 | 7.4 | 5.8 | 4.8 | 4.7 |
| | 20 | 20 | 1.6 | 0.7 | 0.5 | 0.4 |
| | 40 | 20 | 0.5 | 0.6 | 0.4 | 0.3 |
| Cys | 10, 20, 40 | dose | **0** | 0 | 0 | 0 |

### 3.7 Mechanistic claims [Q]
Cys + α-dicarbonyl → Schiff base → (i) thiazolidine → −HCOOH → thiazole (2-acetylthiazole from MGO;
2-thiazolecarboxaldehyde from GO, **not detected**); (ii) −CO2 → carbonylcysteimine (keto/enol);
keto form hydrolyses to **cysteamine + recovered dicarbonyl** ("predominant" for GO); enol form
hydrolyses to Strecker aldehyde + α-aminocarbonyl → pyrazines. "Cys preferentially reacted with MGO
first rather than GO to exclusively generate 2-acetylthiazole at a high yield." None of the
intermediates (Schiff bases, thiazolidines, carbonylcysteimines, cysteamine, aminoacetaldehyde) is
measured.

---

## 4. DERIVED QUANTITIES [Z]

**Molar yields** (MW: 2-acetylthiazole 127.16, furfural 96.08, norfuraneol 114.1, pyrazine 80.09,
methylpyrazine 94.11):

| quantity | µg/L | µM | per reactant |
|---|---|---|---|
| 2-acetylthiazole, 20 Cys + 20 MGO, plateau | 570 | **4.5** | **0.022 % of MGO**; MGO consumed in 15 min ≈ 17 mM → **~0.03 % of the MGO consumed** |
| 2-acetylthiazole, 10 Cys + 20 MGO | 390 | 3.1 | 0.031 % of the Cys |
| 2-acetylthiazole, 20 ARP + 10 Cys, 60 min | 850 | 6.7 | 0.067 % of the Cys; ~1.5 % of the 0.44 mM MGO the ARP releases alone |
| furfural, 20 ARP, 60 min | 1503 | 15.6 | 0.078 mol % of the ARP |
| norfuraneol, 20 ARP | 16.4 | 0.14 | 0.0007 mol % |
| MFT, 20 ARP + 20 Cys, 60 min | 78.4 | 0.69 | 0.0034 % of the Cys |
| pyrazine, best arm (ARP + 40 Cys) | 8.7 | 0.11 | — |
| GO from ARP, 60 min | — | 1570 | **7.9 mol % of the ARP** |
| MGO from ARP, 60 min | — | 440 | **2.2 mol % of the ARP** |

**The "high yield" of the abstract is ~0.02–0.07 % of the cysteine or MGO consumed.** > 99.9 % of the
consumed Cys and MGO goes to unmeasured products.

**Cysteine + MGO consumption rate, lower bound.** 20 + 20 mM, MGO 20 → 2.9 mM in 900 s; if the loss
were second order in equal partners, k ≈ (1/2.9 − 1/20) mM⁻¹ / 900 s ≈ **0.33 M⁻¹ s⁻¹ at 120 °C,
initial pH 7.5**. Not clean: MGO also disappears by other routes (10 Cys + 20 MGO loses 12.9 mM MGO
for 10 mM Cys in 15 min), and no MGO-only or Cys-only control is shown.

**Cysteine disappears faster than any partner can explain.** In 40 mM Cys + 20 mM MGO, **26.7 mM Cys
vanishes between 15 and 30 min while only ~1.4 mM MGO is consumed**. Apparent first-order loss
≥ ln(26.7/0.6)/900 s ≈ **4 × 10⁻³ s⁻¹** (0.6 mM = the calibration intercept, Flag 5). For comparison,
the model's measured thermolysis set (`ZHENG_CYSTEINE_THERMOLYSIS`, Zheng & Ho 1994, interpolated to
pH 7.5) gives **2.7 × 10⁻⁵ s⁻¹ at 120 °C** (2.4 % in 15 min): **≥ 150× slower**. Likewise in 20 Cys +
20 GO, Cys stops at ~1.9 mM once GO is gone, but in 10 and 40 Cys + 20 GO it reaches 0 by 15 min.
Candidate explanations, none tested in the paper: oxidation to cystine (bottle headspace and O2 not
stated), reaction with dicarbonyl-derived intermediates or products (thiazolidines, cysteamine
adducts), or the Cys assay's floor.

**Mass balance on dicarbonyls.** 10 Cys + 20 MGO + 20 GO consumes 7.9 mM MGO + 12.6 mM GO = 20.5 mM
of dicarbonyl for 10 mM Cys in 15 min: at least half of the dicarbonyl loss is not via cysteine.

---

## 5. WHAT THE MODEL CAN TAKE FROM IT

The model's only cysteine + dicarbonyl step is `r_cys_actz`: **Cys + MGO → 2-acetylthiazole (ACTZ) +
CBX**, stoichiometric, `k_cys_actz` (thiol_assembly family barrier). There is no Cys + GO step.

- **(a) The strongest content: Cys + MGO is mostly NOT a 2-acetylthiazole channel.** The measured
  2-acetylthiazole is ~0.03 % of the Cys + MGO consumed. A single stoichiometric step cannot fit both
  observables:
  - if `k_cys_actz` is set to reproduce the 2-acetylthiazole level (as the Zhou 2023 hold-out does:
    582 / 665 µg/L at pH 8), it **under-states cysteine consumption by dicarbonyls by ~10³–10⁴**;
  - if it were set to the measured Cys/MGO decay (≥ ~0.3 M⁻¹ s⁻¹ apparent, 120 °C, pH 7.5), it would
    **over-predict 2-acetylthiazole by ~3000×**.
  The structural fix the data point to is a **Cys + MGO (and Cys + GO) adduct sink** with
  2-acetylthiazole as a ~10⁻⁴ branch. That sink consumes cysteine directly, i.e. it removes thiol and
  H2S precursor, which bears on the open "unidentified thiol sink" question. It is a cysteine sink, not
  a thiolate-of-product sink, and it is measured only at pH 7.5. **Usable as a fit target / hold-out
  for the rate of Cys loss with fed MGO or GO at 120 °C, pH 7.5 (initial), under the caveats below.**
- **(b) ARP → GO and MGO fragmentation, 120 °C, pH 7.5:** GO 1.57 mM and MGO 0.44 mM at 60 min from
  20 mM GR-ARP (7.9 and 2.2 mol %), with time points at 15/30/45 min. A direct test of the model's
  Amadori → dicarbonyl fragmentation if the lane carries a pentose-glycine ARP (the model's MGO source
  is `r_ama_mgo`); GO from a pentose ARP is not in the network today.
- **(c) Product branching vs Cys dose in the ARP pot** (60 min): 2-acetylthiazole **falls** 850 → 300
  → 120 µg/L and pyrazine **rises** 1.5 → 6.3 → 8.7 µg/L as Cys goes 10 → 20 → 40 mM. A falling
  2-acetylthiazole with more Cys is impossible for a pure Cys + MGO → ACTZ step with MGO limiting;
  it needs Cys to divert MGO or its precursors (consistent with (a)).
- **(d) Shape constraint for the hold-out already in the model:** Liu's Cys + MGO plateau (566 µg/L at
  20 + 20 mM, pH 7.5) is within 15 % of Zhou 2023's MGO-fed value at pH 8 (665 µg/L) and Zhou's ARP
  value (582 µg/L). Same laboratory, same SPME method: a reproducibility check, not independent
  evidence.
- **(e) pH: no axis.** One initial pH (7.5), unbuffered as far as the text says, not re-measured.
  Cannot test any pH law; cannot test cysteine thermolysis.
- **(f) MFT 78.4 µg/L (20 ARP + 20 Cys, 1 h, pH 7.5)** is a single level from a glycine ARP pot;
  usable only as a qualitative "MFT forms from a pentose ARP + Cys at pH 7.5" check.

---

## 6. FLAGS

1. **Table 2 does not match Figure 2.** Table 2 (stated: 20 ARP + 20 Cys, 1 h) prints 2-acetylthiazole
   302.4, pyrazine 1.5, methylpyrazine 3.0. Fig. 2c at 20 Cys gives 2-acetylthiazole ≈ 300 (matches),
   but Fig. 2b at 20 Cys gives pyrazine ≈ 6.3, methylpyrazine ≈ 2.7; **Table 2's pyrazine values are
   Fig. 2b's 10 mM Cys bars**. One of the two is mislabelled.
2. **Furfural (1503 µg/L, Table 1) exceeds Fig. 2a's "total volatile flavor compounds" for GR-ARP
   (≈ 970 µg/L).** The total cannot be smaller than one of its components; Fig. 2a's basis is unclear.
3. **Figure 5a–c y-axes say "mmol/L".** 250 mmol/L of 2-acetylthiazole from 20 mM MGO is impossible;
   the values match the µg/L scale of Figs. 3a/4a. Read as µg/L.
4. **Figure 5 has doses (25, 30, 35 mM Cys) not in the Methods**, and its (d)–(f) panels carry only
   10/20/40.
5. **Cysteine assay floor.** ELSD calibration y = −40.196x² + 97379x + 58661 has an intercept
   equivalent to **~0.6 mM** (58661/97379); "0" on Figs. 3c, 4c, 5f means "below ~0.6 mM", not zero.
   The quadratic coefficient is negligible over 0–40 mM.
6. **Volatile calibrations "ln y = a·x"** pass through y = 1 µg/L at zero area and are exponential in
   peak area — an unusual form; taken as printed. Levels depend entirely on it.
7. **No buffer, no post-heating pH, no volume/headspace, no atmosphere stated** for the 120 °C runs.
   Cysteine + dicarbonyl chemistry releases CO2 and formic acid (their own Fig. 1); the pH almost
   certainly drifts below 7.5.
8. **Fig. 2d caption says "with and without Cys addition" but shows one series per dicarbonyl.**
9. **Fig. 3c legend markers are swapped relative to 3a/3b** (10 mM = circle in 3c, square in 3a/b);
   the curves are identified by their t = 0 values, so the reading is unaffected.
10. **"Cys preferentially reacted with MGO"** is inferred from product yields that are ~10⁻⁴ of the
    reactants; the measured consumption (Figs. 4b, 5d,e) shows **GO consumed faster than MGO**, which
    the authors also state. The priority claim is about the 2-acetylthiazole branch, not about where
    cysteine goes.
11. **GR-ARP purity and identity** (α-ribulosyl glycine) are asserted, not shown; trace free ribose or
    glycine would change Fig. 2.
12. **Registry (`data/keys/compounds.yml`)**: `2_methyl_3_furanthiol`, `furfural`, `norfuraneol`,
    `methylpyrazine` keyed; **2-acetylthiazole not keyed**, although the engine carries it as `ACTZ`.
    **Unsubstituted pyrazine is an alias of the class id `pyrazines`**, so a "pyrazine" row would
    resolve to the whole class; key it separately before scoring anything from Figs. 2b/4a/5b.
