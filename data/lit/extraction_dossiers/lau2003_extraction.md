# Lau, Tang, Taub, Yang, Edwards & Mao 2003 — EXTRACTION (20 % whey-protein-concentrate gel + 1 % D-ribose, 79 % water, sealed glass capillaries in an oil bath at 116 / 121 / 126 / 131 °C, up to ~40 min; M-2 = norfuraneol by anion-exclusion HPLC-DAD 285 nm in peak area per g; pseudo-first-order formation to a temperature-independent plateau; Ea by three regression methods; one 915 MHz microwave validation)

### NORFURANEOL (M-2) FORMATION FROM RIBOSE + PROTEIN AMINES IS FIRST ORDER WITH k = 0.110-0.281 min⁻¹ OVER 116-131 °C AND **Ea 81.5 ± 3.2 kJ/mol** (authors' preferred one-step non-linear fit; 64.0 ± 8.7 by two-step, 122 ± 30 by multi-linear) — i.e. close to the engine's lumped 64 kJ/mol, NOT to the 120-147 kJ/mol of furfural. The plateau M∞ is the same at all four temperatures (10 500-11 500 area/g) and **no M-2 loss is visible**; a decay fit of the 116 °C curve bounds k_decay ≲ 6 × 10⁻³ min⁻¹ (~2σ). pH is NOT reported.

**Source on disk:** `data/articles/Lau2003.pdf` (9 pp., Elsevier publisher PDF, clean text layer; file name capitalised).
Prose and Tables 1-3 from `pdftotext -layout`, checked on 130/300-dpi renders (`scratch/lit8/lau-6.png`, `lau-7.png`).
Fig. 3 digitised [D] from the 130-dpi render (read ±150 area/g). Arithmetic: `scratch/lit8/lau_fit.py` (maillard env).

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of chemical marker formation in whey protein gels for studying microwave sterilization" |
| Authors | M. H. Lau, **J. Tang** (corresponding, jtang@mail.wsu.edu), I. A. Taub, T. C. S. Yang, C. G. Edwards, R. Mao — Washington State University; US Army Natick Soldier Center |
| Venue | Journal of Food Engineering **60** (2003) 397-405. Received 14 Oct 2002, accepted 8 Feb 2003 |
| DOI as printed | **10.1016/S0260-8774(03)00063-3** ✔ matches the brief. **Identity verified.** |
| Naming | M-2 = 4-hydroxy-5-methyl-3(2H)-furanone = norfuraneol = engine `NF`. M-1 (cited) = DDMP; M-3 = HMF |

## 1. Why it matters

The engine's NF formation rides the lumped 64 kJ/mol barrier; the furfural barrier swap under discussion would move that
lumped value. B61 (`results/validation/kinetic_core_b61_prereg.md`) adds a hydroxide-catalysed NF decay, measured at one
point (Mikami 2017: 0.032 min⁻¹, pH 8, 100 °C) with a DECLARED barrier of 103.7 kJ/mol. Lau is a candidate independent
Ea for ribose + amine → NF, and a candidate test of NF stability at 116-131 °C.

## 2. Methods as they matter to a model (pp. 399-401)

- **Matrix.** 20 g WPC (Alacen 882, 80 % protein dry basis) + 1 g D-ribose + 79 g deionised water, dispersed at 20 °C;
  gel set at 80 °C / 40 min (water bath) before kinetics; "minimum marker yield (< 8 %)" from the set step.
  Ribose 1 % w/w = **66.6 mmol/kg gel ≈ 84 mmol/kg water** (mine; WPC moisture not stated). Protein ~16 % w/w. **Free
  amine / lysine content NOT reported** (only "lysine, arginine, histidine and methionine" named as reactive).
- **pH: NOT REPORTED** (the only pH in the methods is a literature WPC denaturation range "in 0.7 M phosphate buffer at
  pH 6.0"). The paper's framing is "low acid conditions (pH > 5)". No buffer added.
- **Heating.** ~2 mL into glass capillaries ("inner diameter 1.5 mm; outer diameter 1.55 mm; length 15 mm" — see
  Flags), flame-sealed, oil bath at **116, 121, 126, 131 °C**, "different time intervals"; ice-water quench. Duplicates.
  Come-up time not stated.
- **Analysis.** Gel homogenised in 4 mL 10 mM H2SO4, centrifuged, 0.45 µm filtered; Bio-Rad anion-exclusion column,
  10 mM H2SO4 1 mL/min, DAD 285 nm, M-2 at 5.6 min. **Result = peak area per g gel; NO calibration to mass or mol.**
- **Model.** dM/dt = k(M∞ − M)ⁿ; n by non-linear regression (Table 1) and by r² of linearised plots (Table 2); then first
  order: M = M∞ − (M∞ − M0)e^(−kt), Arrhenius about Tref = 123.5 °C (396.8 K), R = 1.987 cal/mol/K.

## 3. Numbers

### Table 1 (p. 402), non-linear fit of order, M0, M∞

| T (°C) | n | M0 (area/g) | M∞ (area/g) | r² |
|---|---|---|---|---|
| 116 | 0.95 ± 0.13 | 817 ± 171 | 10 700 ± 300 | 0.998 |
| 121 | 0.84 ± 0.10 | 812 ± 178 | 10 500 ± 400 | 0.997 |
| 126 | 1.04 ± 0.20 | 786 ± 132 | 11 500 ± 700 | 0.999 |
| 131 | 0.77 ± 0.19 | 873 ± 119 | 10 900 ± 200 | 0.999 |

Measured M0 = 836 ± 40 area/g. ANOVA: M∞ not different between temperatures (p > 0.05).

### Table 2 (p. 402), r² of linearised plots: first order best at all T (0.977-0.998); half order 0.962-0.994.

### Table 3 (p. 403), first-order k (min⁻¹) and Ea

| | two-step | multi-linear | **one-step non-linear** | M-1, broccoli (Kim & Taub 1993) |
|---|---|---|---|---|
| k 116 °C | 0.112 ± 0.005 | 0.089 ± 0.008 | **0.110 ± 0.005** | 0.0099 |
| k 121 °C | 0.166 ± 0.011 | 0.094 ± 0.008 | **0.152 ± 0.006** | 0.0146 |
| k 126 °C | 0.190 ± 0.008 | 0.238 ± 0.018 | **0.207 ± 0.007** | 0.0213 |
| k 131 °C | 0.243 ± 0.011 | 0.311 ± 0.021 | **0.281 ± 0.013** | 0.0309 |
| k_ref (123.5 °C) | 0.172 ± 0.006 | 0.159 ± 0.020 | **0.178 ± 0.005** | 0.0198 |
| Ea (kcal/mol) | 15.29 ± 2.09 | 29.22 ± 7.16 | **19.48 ± 0.76** | 23.7 |
| Ea (kJ/mol, mine) | 64.0 ± 8.7 | 122.3 ± 30.0 | **81.5 ± 3.2** | 99.2 |
| r² | 0.964 | 0.893 | 0.995 | |

Refit of the printed k's (mine): two-step 64.4, multi-linear 122.3, non-linear 81.7 kJ/mol ✓.
Microwave validation (Table 4, 915 MHz, F0 3.2 and 4.9 min): measured M/M∞ 0.58 ± 0.05 and 0.63 ± 0.03; the non-linear
parameters predict 0.54 / 0.64, two-step 0.57 / 0.66, multi-linear 0.46 / 0.56 (rejected).

### Fig. 3 (p. 402), one M-2 time course, temperature not labelled [D]

| t (min) | 0 | 2 | 4 | 6 | 8 | 10 | 15 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| M-2 (area/g) | 900 | 2 800 | 4 350 | 5 600 | 6 550 | 7 900 | 8 900 | 9 600 | 10 550 | 10 600 |

A first-order fit gives M0 836, M∞ 10 790, k = 0.111 min⁻¹ — **this is the 116 °C series** (Table 1/3 match; inferred).
Fig. 4 (ln(M∞ − M) vs t) shows data to ~20 min at 116 °C, ~11 min at 121 °C, ~8 min at 126 °C and ~9 min at 131 °C, so
**the plateau is directly observed only at 116 °C** (Fig. 3); Fig. 5 (3-D surface) carries points to ~40 min at the low
temperatures but cannot be read per temperature.

### Arithmetic (mine): an NF decay bound from Fig. 3 (116 °C)

Model: ribose-limited first-order formation (k_f) with first-order M-2 loss (k_d); M0, M∞, k_f free, k_d fixed; σ ≈ 147
area/g (rms of the k_d = 0 fit).

| k_d (min⁻¹) | 0 | 1e-3 | 2e-3 | 3e-3 | 4.3e-3 | 6e-3 | 1e-2 |
|---|---|---|---|---|---|---|---|
| SSE (1e3) | 216 | 209 | 213 | 226 | 252 | 298 | 429 |
| Δχ² vs best | 0.3 | 0 | 0.2 | 0.8 | 2.0 | 4.1 | 10 |

So **k_d(116 °C) ≲ 6 × 10⁻³ min⁻¹ at ~2σ, ≲ 1 × 10⁻² at ~3σ**; k_d = 0 is fully consistent. One curve, figure-read.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| ribose + protein amines → NF (overall) | Ea | **81.5 ± 3.2** (preferred); 64.0 ± 8.7; 122 ± 30 | kJ/mol | WPC gel 20 %, ribose 1 %, 79 % water, pH n.r., 116-131 °C | Table 3 | measured (method-dependent) |
| same | first-order k | 0.110 / 0.152 / 0.207 / 0.281 | min⁻¹ | 116 / 121 / 126 / 131 °C | Table 3 (non-linear) | measured |
| same | k_ref | 0.178 ± 0.005 | min⁻¹ | 123.5 °C | Table 3 | measured |
| NF plateau | M∞ T-independence | 10 500-11 500 (n.s.) | area/g | 116-131 °C | Table 1 | measured (relative) |
| NF loss | k_d upper bound | **≲ 6e-3 (2σ)** | min⁻¹ | 116 °C, gel pH n.r. | Fig. 3 [D], arithmetic | derived bound |
| NF molar yield | — | **NOT REPORTED** (no calibration) | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** an apparent Ea band for ribose + amine → NF at sterilisation temperatures (64-82 kJ/mol on the two accepted
methods); a weak upper bound on NF loss at 116 °C; the fact that the NF yield ceiling does not change between 116 and
131 °C.

**Cannot:** no pH, no amine concentrations, no absolute NF; not a buffered aqueous pot (protein gel); no NF decay data at
131 °C beyond ~9 min; no sulfur; the overall k lumps ribose → Amadori → 2,3-enolisation → DPO → NF.

## 6. Flags

1. Capillary "inner diameter 1.5 mm; outer diameter 1.55 mm; length 15 mm" cannot hold "about 2 ml" (26 µL); Pandit 2006
   (same lab) uses 1.75 mm × 15 **cm** (~0.36 mL). Probably 15 cm; the fill volume is still inconsistent.
2. Ea depends on the regression method by a factor of ~2 (64 vs 122 kJ/mol); the authors' validation rejects the 122.
3. "kcal mol⁻¹ K⁻¹" in Table 3 is a unit misprint for kcal/mol.
4. M-1 (DDMP) comparison values are from broccoli extract (Kim & Taub 1993), not this paper.
5. Table 3 "non-linear" k's are back-calculated from k_ref and Ea, so they are Arrhenius-smooth by construction.

## 7. Use in the model

- **(a) The furfural barrier swap.** Lau's preferred NF-formation Ea (81.5 kJ/mol; 64.0 two-step) sits near the engine's
  lumped 64 kJ/mol and far below furfural's 120-147. **If the swap raises the SHARED lumped barrier, NF formation would
  leave Lau's band; the swap should be furfural-specific (split `k_tdp_fur`, and `k_glc_fur`, from the lumped value), with
  the NF route kept near 64-95 kJ/mol** (Pandit 2006 gives 93). This also supports `SULFUR_FLUX_MAP.md` sec. 5's reading
  that the furfural and norfuraneol branches need different temperature dependences.
  Transfer caveat: Lau's Ea is the overall ribose → NF rate in a protein gel; which engine step carries it depends on
  which step is rate-limiting in the engine (pentose → DPO or DPO → NF).
- **(b) B61's declared 103.7 kJ/mol NF-decay barrier.** No NF loss is observed. B61's law (k_OH(100 °C, pH 5) = 10^−4.495
  min⁻¹, 10^(pH−5), 103.7 kJ/mol, plus the neutral term) predicts at 116 °C: **1.6e-3 (pH 6.0), 2.8e-3 (6.3), 4.3e-3 (6.5),
  1.3e-2 min⁻¹ (pH 7.0)**. Against the Fig.-3 bound (≲ 6e-3 at 2σ) this is **consistent if the gel pH ≤ ~6.6 and
  disfavoured at pH ≥ 7**. Whey-protein gel pH is not in the paper (WPC dispersions are typically ~6.3-6.8 — outside
  knowledge, not a measurement here). **Lau neither tests nor constrains the 103.7 kJ/mol barrier itself**: there is no
  decay measured at two temperatures, and the T-independent plateau is insensitive (with B61 at pH 6.5 the plateau would
  fall only ~3 % more at 131 °C than at 116 °C).
- **Hold-out suitability:** poor for the engine as it stands (protein matrix, unknown pH and amine pool, relative units).
  Usable as a **shape** check only: k(131)/k(116) = 2.55 (non-linear) for the NF-appearance rate.
