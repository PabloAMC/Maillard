# Brands & van Boekel 2002c — EXTRACTION
### *"Kinetic Modeling of Reactions in Heated Monosaccharide−Casein Systems"*
### Glucose/fructose (and galactose/tagatose) + 3 % caseinate, pH 6.8 phosphate, 90–130 °C: four successive multiresponse models, rate constants at five temperatures, activation energies. The only sugar → C-fragment step is **fructose → 2 trioses (step 5)**, fitted through acetic acid; there is no glycolaldehyde anywhere. Read 2026-10-27.

**Full read of `data/articles/brands2002c.pdf`** (text layer, 15 pages; Schemes 1–5, Figure 1 and
equations 6.1–6.11 checked on 150-dpi rasters). Companions: `brands2001_extraction.md` (the network
paper, figures only, no constants) and `brands2002b_extraction.md`. The trunk lane's backbone is
Martins 2005 (`martins2005_extraction.md`), not this paper.

**Provenance:** **[M]** measured · **[F]** fitted by the authors (multiresponse, determinant criterion) ·
**[Q]** qualitative · **[D]** derived by this extraction from printed numbers · **[Z]** extrapolated by
this extraction outside the paper's range · **[!]** flag.

---

## 0. IDENTITY — matches

| field | value as printed |
|---|---|
| Authors | Carline M. J. Brands, Martinus A. J. S. van Boekel* |
| Venue | *J. Agric. Food Chem.* 2002, 50 (23), 6725–6739 |
| DOI | 10.1021/jf011164h (printed p. 6725 footer) |
| Affiliation | Product Design and Quality Management Group, Wageningen University |
| Dates | received 28 Aug 2001; revised 18 Jul 2002; accepted 1 Aug 2002; web 11 Oct 2002 |
| Funding | EC FAIR CT96-1080 ("Optimization of the Maillard Reaction") |
| PDF character | 15 pp.; Figure 1 (network), Schemes 1–5 (kinetic models), Tables 1–7 (rate constants), Figures 2–13 (fits, not digitised) |

---

## 1. THE ONE-PARAGRAPH ANSWER

A monosaccharide–casein network (Brands 2001) is cut down to a fittable kinetic model and refined
four times. The final model (**Scheme 4**, 13 parameters, log posterior probability −31.233) has
glucose ⇄ fructose isomerisation (k1, k2), sugar → formic acid + C5 (k3, k4), **fructose → 2 trioses
(k5)**, trioses → acetic acid (k6), glucose + lysine → Amadori (k7), Amadori → lysine + 2 trioses
(k8), Amadori → lysine + Cn (k9, fitted to zero), fructose → C6 (k12), Cn + lysine → melanoidins
(k13), **glucose → C6 (k14)** and Cn + arginine → melanoidins (k15). Activation energies over 90–130
°C are mostly 100–160 kJ/mol (Table 6). The headline chemistry: **ketoses do not react with lysine
directly; only their degradation products do** (k10 → 0 in Scheme 3, dropped in Scheme 4). The only
step that breaks a hexose into short fragments is **k5, fructose → 2 C3 (glyceraldehyde/
dihydroxyacetone, "both fast reacting toward methylglyoxal")**; it is fitted through acetic acid
(the trioses are never measured), its Ea is 113 ± 4 kJ/mol in Scheme 4, and its value at 120 °C is
7.1e-4 /min. Glucose has no direct C3 or C2 step: it reaches the trioses only through fructose (or
the Amadori compound, k8). **No glycolaldehyde, no C2 fragment other than acetic acid, and no
methylglyoxal measurement anywhere in the paper.**

---

## 2. METHODS [M]

| item | as printed |
|---|---|
| Standard system | 150 mM sugar (glucose, fructose, galactose or tagatose) + sodium caseinate 3 % w/w (DMV, 90 % protein) in **0.1 M phosphate, pH 6.8**; sugar : lysine residues ≈ 10 : 1 |
| Variants | 75 mM sugar (Fig. 11); 1.5 % caseinate (Fig. 12); initial pH 5.9 at 120 °C (Fig. 5, Table 3); pH-stat at constant pH 6.8, 100 °C (Fig. 7) |
| Heating | oil bath, screw-capped Schott tubes 16 × 160 mm; **90, 100, 110, 120, 130 °C**; times include 2–3 min heat-up; cooled in ice water; at least duplicate |
| pH | initial only (uncontrolled except the pH-stat run); falls with acid formation, handled in Schemes 3–4 by multiplying [lys] by 10^−ΔpH |
| Analyses | sugars and acids by HPLC after Sephadex G25 separation; available lysine by OPA fluorescence (ex 340 / em 430 nm, SDS); Amadori via furosine (HPLC, conversion factor **3.1**); melanoidins from A420 of the protein fraction with ε = **500 L mol⁻¹ cm⁻¹** (expressed as sugar units incorporated) |
| Fitting | numerical integration + determinant criterion (Gregpak/Athena Visual Workbench); model discrimination by posterior probability; least-squares comparison with DynaFit; temperature by reparametrised Arrhenius k = X exp(−Y·Ea), X = k0 exp(−Ea/RT_av), Y = (1/R)(1/T − 1/T_av) |
| Units | first-order k in **min⁻¹**; bimolecular k7, k10, k13, k15 in **L mmol⁻¹ min⁻¹** |

---

## 3. THE MODELS (step numbers as printed)

| step | Scheme 1 | Scheme 2 | Scheme 3 | Scheme 4 (final) | Scheme 5 (aldose/ketose general) |
|---|---|---|---|---|---|
| 1 | Glc → Fru | same | same | same | aldose → ketose |
| 2 | Fru → Glc | same | same | same | ketose → aldose |
| 3 | Glc → formic + C5 | same | same | same | aldose → formic + C5 |
| 4 | Fru → formic + C5 | same | same | same | ketose → formic + C5 |
| **5** | **Fru → 2 trioses** (eq 1.5: d[triose]/dt = 2k5[fru] − …) | same | same | same (eq 6.5) | ketose → trioses |
| 6 | trioses → acetic + Cn | same | trioses → acetic + Cn | trioses → acetic (eq 6.4) | same |
| 7 | Glc + Lys-R → Amadori | same | same | same (× 10^−ΔpH) | aldose + Lys-R → Amadori |
| 8 | Amadori → Lys-R + acetic | same | Amadori → Lys-R + trioses | **Amadori → Lys-R + 2 trioses** (eq 6.5) | same |
| 9 | Amadori → AMP | same | same | Amadori → Lys-R + Cn (**k9 = 0**) | — |
| 10 | Fru + Lys-R → AMP | — | same | dropped | — |
| 11 | AMP → melanoidins | same | same | dropped | — |
| 12 | — | Fru → Cn | same | Fru → C6 | ketose → C6 |
| 13 | — | Cn + Lys-R → AMP | same | Cn + Lys-R → melanoidins | same |
| 14 | — | — | — | **Glc → C6** | aldose → C6 |
| 15 | — | — | — | Cn + Arg-R → melanoidins | same |

Figure 1 (the network, [Q]): aldose ⇄ 1,2-enediol ⇄ ketose ⇄ 2,3-enediol; 1,2-enediol → Cn and →
formic acid + C5; ketose → trioses → acetic acid; 2,3-enediol → **acetic acid + C4** and → Cn; the
same three branches for the Amadori/Heyns enaminols with lysine-R released. Text p. 6736 identifies
C5 as "possibly 2-deoxyribose", C3 as trioses (glyceraldehyde, 1,3-dihydroxyacetone, "both fast
reacting toward methylglyoxal"), C6 as deoxyosones [Q].

---

## 4. EVERY PRINTED NUMBER

### 4.1 Table 1 — Scheme 1, 120 °C, initial pH 6.8, ±95 % CI [F] (p. 6728)

| k | determinant | least squares |
|---|---|---|
| k1 (min⁻¹) | 0.010 28 ± 0.000 46 | 0.010 20 ± 0.000 21 |
| k2 (min⁻¹) | 0.005 09 ± 0.000 23 | 0.005 08 ± 0.000 22 |
| k3 (min⁻¹) | 0.000 47 ± 0.000 04 | 0.000 66 ± 0.000 17 |
| k4 (min⁻¹) | 0.001 10 ± 0.000 06 | 0.001 05 ± 0.000 18 |
| **k5 (min⁻¹)** | **0.007 12 ± 0.002 14** | 0.004 61 ± 0.000 46 |
| k6 (min⁻¹) | 0.004 39 ± 0.001 43 | 0.006 73 ± 0.001 52 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 18 ± 0.000 01 | 0.000 24 ± 0.000 02 |
| k8 (min⁻¹) | 0.111 34 ± 0.011 46 | 0.112 24 ± 0.046 83 |
| k9 (min⁻¹) | 0.143 59 ± 0.017 10 | 0.168 31 ± 0.075 95 |
| k10 (L mmol⁻¹ min⁻¹) | 0.000 15 ± 0.000 01 | 0.000 17 ± 0.000 03 |
| k11 (min⁻¹) | 0.125 14 ± 0.030 48 | 0.073 02 ± 0.026 21 |

Log posterior probability −31.992. Text: the CIs of k5, k6 (ketose → trioses → acetic acid) and k11
"were somewhat wide … the data obtained do not contain enough information to estimate them" [Q].

### 4.2 Table 2 — Scheme 1, 90–130 °C, Arrhenius [F] (p. 6730)

| k | 90 °C | 100 °C | 110 °C | 120 °C | 130 °C | Ea (kJ/mol) |
|---|---|---|---|---|---|---|
| k1 (min⁻¹) | 0.000 39 | 0.001 18 | 0.003 42 | 0.009 34 | 0.024 31 | 126 ± 2 |
| k2 (min⁻¹) | 0.000 26 | 0.000 72 | 0.001 93 | 0.004 91 | 0.011 94 | 117 ± 2 |
| k3 (min⁻¹) | 0.000 01 | 0.000 04 | 0.000 12 | 0.000 36 | 0.001 03 | 137 ± 7 |
| k4 (min⁻¹) | 0.000 03 | 0.000 11 | 0.000 42 | 0.001 49 | 0.004 98 | 159 ± 5 |
| **k5 (min⁻¹)** | 0.000 12 | 0.000 37 | 0.001 09 | **0.003 06** | 0.008 16 | **129 ± 17** |
| k6 (min⁻¹) | 0.000 71 | 0.001 87 | 0.004 66 | 0.011 10 | 0.025 31 | 109 ± 19 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 01 | 0.000 03 | 0.000 08 | 0.000 20 | 0.000 46 | 114 ± 2 |
| k8 (min⁻¹) | 0.005 10 | 0.015 28 | 0.043 23 | 0.116 01 | 0.296 47 | 124 ± 4 |
| k9 (min⁻¹) | 0.006 16 | 0.019 14 | 0.056 07 | 0.155 54 | 0.410 17 | 128 ± 4 |
| k10 (L mmol⁻¹ min⁻¹) | 0.000 00 | 0.000 01 | 0.000 04 | 0.000 12 | 0.000 34 | 138 ± 4 |
| k11 (min⁻¹) | 0.016 31 | 0.031 72 | 0.059 58 | 0.108 38 | 0.191 38 | 75 ± 11 |

The ± on Ea is not labelled; by the paper's convention elsewhere, read as 95 % interval [!].
Fit quality: glucose–casein "very well"; fructose–casein poor for formic acid, lysine loss and
browning [Q].

### 4.3 Scheme 2 [F] (p. 6730)

n (mean carbon number of Cn) = **5.8 ± 0.4**. Step 11's activation energy "approaches zero" (not
shown). No rate-constant table printed for Scheme 2.

### 4.4 Table 3 — Scheme 1, 120 °C, **initial pH 5.9**, ±95 % CI [F] (p. 6732)

| k | value |
|---|---|
| k1 (min⁻¹) | 0.002 57 ± 0.000 31 |
| k2 (min⁻¹) | 0.000 87 ± 0.000 03 |
| k3 (min⁻¹) | 0.000 07 ± 0.000 01 |
| k4 (min⁻¹) | 0.000 19 ± 0.000 02 |
| **k5 (min⁻¹)** | **0.001 49 ± 0.000 69** |
| k6 (min⁻¹) | 0.002 14 ± 0.001 05 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 08 ± 0.000 00 |
| k8 (min⁻¹) | 0.025 28 ± 0.002 25 |
| k9 (min⁻¹) | 0.073 34 ± 0.006 09 |
| k10 (L mmol⁻¹ min⁻¹) | 0.000 03 ± 0.000 00 |
| k11 (min⁻¹) | 0.120 67 ± 0.030 18 |

"The rate of all reactions decreased with decreasing pH" [Q]. k5(pH 5.9)/k5(pH 6.8) = 0.21 [D].
HMF rises at lower pH but stays "an order of magnitude lower than those of other reaction products"
[Q, Fig. 6].

### 4.5 Table 4 — Scheme 3, 120 °C, pH 6.8, with 10^−ΔpH on lysine, ±95 % CI [F] (p. 6734)

| k | value |
|---|---|
| k1 (min⁻¹) | 0.010 39 ± 0.000 47 |
| k2 (min⁻¹) | 0.005 03 ± 0.000 24 |
| k3 (min⁻¹) | 0.000 47 ± 0.000 04 |
| k4 (min⁻¹) | 0.001 09 ± 0.000 06 |
| **k5 (min⁻¹)** | **0.001 04 ± 0.000 20** |
| k6 (min⁻¹) | 0.053 82 ± 0.018 40 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 32 ± 0.000 03 |
| k8 (min⁻¹) | 0.154 28 ± 0.037 65 |
| k9 (min⁻¹) | 0.161 65 ± 0.035 38 |
| k10 (L mmol⁻¹ min⁻¹) | 0.000 06 ± 0.000 04 |
| k11 (min⁻¹) | 0.089 85 ± 0.027 29 |
| k12 (min⁻¹) | 0.008 38 ± 0.002 25 |
| k13 (L mmol⁻¹ min⁻¹) | 0.001 75 ± 0.000 44 |

n = **10.5 ± 2.0** (physically impossible, > 6; "To predict enough browning, this was apparently
necessary") [!]. Log posterior probability −39.871 (Scheme 1 + pH effect) → −34.882 (Scheme 3).

### 4.6 Table 5 — Scheme 4, 120 °C, pH 6.8, ±95 % CI [F] (p. 6735)

| k | value |
|---|---|
| k1 (min⁻¹) | 0.009 93 ± 0.000 38 |
| k2 (min⁻¹) | 0.005 12 ± 0.000 23 |
| k3 (min⁻¹) | 0.000 47 ± 0.000 04 |
| k4 (min⁻¹) | 0.001 09 ± 0.000 06 |
| **k5 (min⁻¹)** | **0.000 86 ± 0.000 16** |
| k6 (min⁻¹) | 0.083 80 ± 0.034 49 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 15 ± 0.000 02 |
| k8 (min⁻¹) | 0.131 15 ± 0.027 08 |
| k9 (min⁻¹) | 0 |
| k12 (min⁻¹) | 0.006 23 ± 0.001 33 |
| k13 (L mmol⁻¹ min⁻¹) | 0.002 92 ± 0.000 60 |
| **k14 (min⁻¹)** | **0.002 57 ± 0.000 58** |
| k15 (L mmol⁻¹ min⁻¹) | 0.000 18 ± 0.000 24 |

Log posterior probability −31.233 ("the best model obtained so far"). Arginine residues ≈ 6 mmol/L.

### 4.7 Table 6 — Scheme 4, 90–130 °C, Arrhenius [F] (p. 6735)

| k | 90 °C | 100 °C | 110 °C | 120 °C | 130 °C | Ea (kJ/mol) |
|---|---|---|---|---|---|---|
| k1 (min⁻¹) | 0.000 39 | 0.001 17 | 0.003 36 | 0.009 22 | 0.025 28 | 126 ± 2 |
| k2 (min⁻¹) | 0.000 27 | 0.000 72 | 0.001 92 | 0.004 86 | 0.012 32 | 116 ± 2 |
| k3 (min⁻¹) | 0.000 01 | 0.000 04 | 0.000 13 | 0.000 38 | 0.001 16 | 139 ± 7 |
| k4 (min⁻¹) | 0.000 03 | 0.000 10 | 0.000 39 | 0.001 38 | 0.004 94 | 159 ± 6 |
| **k5 (min⁻¹)** | **0.000 04** | **0.000 11** | **0.000 29** | **0.000 71** | **0.001 75** | **113 ± 4** |
| k6 (min⁻¹) | 0.003 61 | 0.012 65 | 0.043 12 | 0.138 60 | 0.445 57 | 146 ± 11 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 01 | 0.000 03 | 0.000 06 | 0.000 13 | 0.000 28 | 97 ± 3 |
| k8 (min⁻¹) | 0.007 82 | 0.019 22 | 0.046 26 | 0.106 81 | 0.246 59 | 105 ± 5 |
| k9 (min⁻¹) | 0 | 0 | 0 | 0 | 0 | 0 |
| k12 (min⁻¹) | 0.000 10 | 0.000 34 | 0.001 09 | 0.003 37 | 0.010 37 | 141 ± 12 |
| k13 (L mmol⁻¹ min⁻¹) | 0.000 25 | 0.000 61 | 0.001 47 | 0.003 42 | 0.007 94 | 105 ± 10 |
| **k14 (min⁻¹)** | 0.000 13 | 0.000 35 | 0.000 94 | 0.002 43 | 0.006 28 | 119 ± 11 |
| k15 (L mmol⁻¹ min⁻¹) | 0.000 02 | 0.000 04 | 0.000 08 | 0.000 14 | 0.000 25 | **71 ± 103** |

Check [D]: an ordinary least-squares Arrhenius refit of the five printed k5 values gives Ea = 114.7
kJ/mol, consistent with the printed 113 ± 4 (the printed values are rounded to 1 significant figure
at 90 °C).

### 4.8 Table 7 — Scheme 5, galactose/tagatose–casein, 120 °C, pH 6.8, ±95 % CI [F] (p. 6738)

| k | value |
|---|---|
| k1 (min⁻¹) | 0.009 92 ± 0.001 15 |
| k2 (min⁻¹) | 0.015 59 ± 0.001 45 |
| k3 (min⁻¹) | 0.000 70 ± 0.000 08 |
| k4 (min⁻¹) | 0.002 24 ± 0.000 19 |
| **k5 (min⁻¹)** (tagatose → trioses) | **0.002 19 ± 0.000 28** |
| k6 (min⁻¹) | 0.141 47 ± 0.055 28 |
| k7 (L mmol⁻¹ min⁻¹) | 0.000 18 ± 0.000 03 |
| k8 (min⁻¹) | 0.375 25 ± 0.060 70 |
| k12 (min⁻¹) | 0.015 17 ± 0.003 85 |
| k13 (L mmol⁻¹ min⁻¹) | 0.005 09 ± 0.001 30 |
| k14 (min⁻¹) | 0.003 22 ± 0.000 86 |
| k15 (L mmol⁻¹ min⁻¹) | 0 |

### 4.9 Other printed numbers [M]/[Q] (literature values the paper quotes are marked "cited")

| quantity | value | location |
|---|---|---|
| Acyclic fructose, neutral pH | 0.7 % at RT, 13.1 % at 80 °C (cited, Yaylayan 1993) | p. 6728 |
| Acyclic glucose, neutral pH, RT | 0.002 % (Hayward & Angyal) to 0.022 % (Cantor & Peniston) (cited) | p. 6728 |
| Acyclic glucose vs pH | ≈ 3-fold increase pH 6.5 → 7.5 at 25 °C (cited); fructose ×1.1 over pH 2–9 at 80 °C (cited) | p. 6732 |
| Lysine ε-NH2 pKa | 10.53 at RT (cited); unprotonated 0.0186 % at pH 6.8, 0.0023 % at pH 5.9 (25 °C) | p. 6732 |
| C2/C3 fragments vs pH | "formation of C2 and C3 sugar fragments increased with increasing pH" (cited, Hayashi & Namiki 1986); fragments "much more reactive toward amino groups than the original sugars" | p. 6732 [Q] |
| Activation energies "in the order of 120 kJ/mol" | as expected for chemical reactions | p. 6729 [Q] |
| Amadori maximum | same maximum concentration at all temperatures | p. 6729 [Q] |
| Galactose acyclic form | ≈ 10× glucose (cited); tagatose acyclic 15 % below fructose (cited) | pp. 6737–6738 |
| pH-stat run, 100 °C, pH 6.8 | Table 2 100 °C constants predict glucose and fructose "very well", acids, browning and lysine loss under-predicted | p. 6732 [Q, Fig. 7] |
| 75 mM sugar; 1.5 % caseinate | Table 5 constants predict sugars and acids "perfect", lysine "reasonable", fructose browning over-predicted | p. 6736 [Q, Figs. 11–12] |

---

## 5. FLAGS [!]

1. **k5 is not stable across models.** At 120 °C, pH 6.8: 0.00712 (Table 1, Scheme 1, one
   temperature), 0.00306 (Table 2, Scheme 1, five temperatures), 0.00104 (Table 4, Scheme 3), 0.00086
   (Table 5, Scheme 4, one temperature), 0.00071 (Table 6, Scheme 4, five temperatures): a factor of
   **10** depending on model structure. Its Ea moves from 129 ± 17 to 113 ± 4 kJ/mol. The authors
   say the data "do not contain enough information" for k5/k6 in Scheme 1.
2. **The trioses are never measured.** k5 is identified only through acetic acid (eq 6.4) and, in
   Schemes 3–4, through lysine loss and browning via k13. Martins 2005 attributes acetic acid to
   1-deoxyglucosone (step 8), a route this paper lumps into k5/k8 → triose → acetic acid; so k5 also
   carries any non-triose acetic acid.
3. Eq 6.11 (p. 6736) writes the arginine terms with **k14** ("3/6 k14[triose][arg] + 6/6 k14[C6][arg]")
   where Scheme 4 and eqs 6.5, 6.9, 6.10 make it k15; the lysine terms in 6.11 also lack the 10^−ΔpH
   factor used in 6.5–6.9. Typo or model-code difference; not resolvable from the paper.
4. Figure 7 caption says "based on the reaction rate constants of Table 1"; the text (p. 6732) says
   the 100 °C values of Table 2 were used.
5. Figure 12 caption says "initial sugar concentration of 75 mM"; the text says 150 mM sugar with
   1.5 % caseinate. Caption copied from Figure 11.
6. Table 2 k10 at 90 °C printed 0.000 00 and Table 6 k1–k3 and k7 at 90 °C carry 1 significant
   figure; refitting Ea from the table loses precision there.
7. Table 6 k15 Ea 71 ± 103 kJ/mol: unidentified. Table 6 k9 = 0 at every temperature.
8. Scheme 3's n = 10.5 ± 2.0 exceeds its own definition (n ≤ 6).
9. Planck's constant printed with units "J s⁻¹" (p. 6729); typographic.

---

## 6. WHAT THE MODEL CAN TAKE FROM IT

**The need:** `k_glc_ha`, the sulfur lane's `r_glc_c2c3`, Glc → HA + MGO + FRAG_C, first order,
1/min, reference 145 °C. Frozen B9 value log10 k = −8.59 (2.6e-9 /min); B47 refit −4.95 (1.1e-5
/min); unidentified (flat profile).

- **No step in this paper produces hydroxyacetaldehyde (glycolaldehyde).** The only C2 the network
  names is acetic acid (from the 2,3-enediol, Fig. 1, and from trioses, k6). So nothing here measures
  the HA half of `r_glc_c2c3`.
- **The one hexose → C-fragment step is k5, fructose → 2 trioses** (C3 + C3; glyceraldehyde/DHA, which
  the authors say go fast to methylglyoxal). Stoichiometry 1 Fru → 2 C3 (eq 6.5), sugar-only (no amine
  in the step, though casein is present in every pot; Brands 2001 found ketose acid formation about
  equal with and without casein). It bears on the **MGO half only**, and only via fructose.
- **Glucose has no direct C3/C2 step.** Glucose's own losses are k1 (→ Fru), k3 (→ formic + C5), k7
  (+ Lys) and k14 (→ C6 deoxyosones). Glucose reaches trioses through Fru (k1 then k5) or through the
  Amadori compound (k8, Amadori → Lys + 2 trioses; amine-dependent).
- **Numbers at 145 °C [Z]** (15 °C above the paper's 130 °C ceiling; pH 6.8 phosphate, falling during
  heating; casein present):
  - k5 (Table 6, Ea 113 ± 4, from 0.000 71 /min at 120 °C): **5.6e-3 /min** (5.2–6.0e-3 over the Ea
    interval), log10 = −2.25. From Table 2 (Scheme 1, Ea 129): 3.2e-2 /min.
  - Glucose-referenced effective rate to 2 trioses through a quasi-steady fructose pool,
    k1·k5/(k2 + k4 + k5 + k12), Table 6 constants: 6.3e-4 /min at 120 °C, **4.5e-3 /min at 145 °C**
    [Z/D]; only ≈ 5 % of the fructose flux goes to trioses.
  - pH: k5 at initial pH 5.9 is 0.21× the pH 6.8 value (Scheme 1, 120 °C) [D].
- **Verdict:** k5 is an author-fitted [F] constant with an Ea, at the right temperature and pH, but it
  is (a) a fructose step, (b) C3 + C3 with no C2, (c) identified only through acetic acid with an
  unmeasured intermediate, and (d) model-dependent by a factor of 10 at 120 °C. It **cannot anchor
  `k_glc_ha`** as a measured value. It can serve as an **order-of-magnitude soft ceiling** on
  hexose → C3 fragmentation without amine catalysis in pH 6.8 phosphate: ~10^−2.3 /min at 145 °C.
  Both B9 (−8.59) and B47 (−4.95) sit below that ceiling, so it does not discriminate between them.
  Caveat on the ceiling: amino compounds catalyse sugar fragmentation (Hayashi & Namiki, cited here),
  so a cysteine pot could exceed the sugar-only rate.
- **Cross-check for the trunk lane:** k1/k2 (126/116 kJ/mol) agree with Martins 2005's glucose ⇄
  fructose steps (122.6/93.4 kJ/mol) on k1 and not on k2; different system (casein vs glycine,
  150 vs 200 mM). Not a write-back candidate; recorded only.
