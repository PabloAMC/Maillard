# Pandit, Tang, Mikhaylenko & Liu 2006 — EXTRACTION (mashed potato from flakes, 83.12 % water, + 1.5 % D-ribose, sealed capillaries in an oil bath at 116 / 121 / 126 / 131 °C up to 60 min; M-2 = norfuraneol by anion-exclusion HPLC-DAD 285 nm, CALIBRATED against a commercial standard in mg/g; first-order formation to a T-independent plateau; Ea by three regression methods; ribose 1.5-4.5 %, added lysine 0-1 %, salt 0-1 %, five potato sources at 121 °C)

### SECOND, INDEPENDENT NORFURANEOL-FORMATION Ea: **22.23 ± 1.54 kcal/mol = 93.0 ± 6.4 kJ/mol** (authors' preferred one-step non-linear; 96.1 two-step, 81.4 multi-linear), k = 0.030-0.089 min⁻¹ over 116-131 °C, ~3.5× slower than in whey gel (fewer amines). Plateau **C∞ ≈ 0.27 mg/g (≈ 2.4 µmol/g, 2.4 mol % of the ribose)**, independent of T. **No M-2 loss is resolvable at any temperature up to 60 min**; the flat 35→60 min segment at 121 °C with 4.5 % ribose is consistent with either amine exhaustion or a slow loss (≤ ~6 × 10⁻³ min⁻¹). pH is NOT reported.

**Source on disk:** `data/articles/pandit2006.pdf` (9 pp., Elsevier publisher PDF, clean text layer). Prose and Tables 1-3
from `pdftotext -layout`; Figs. 2-5 from 200-dpi renders (`scratch/lit8/pan-6.png`, `pan-7.png`), tagged [D] (read
±0.005 mg/g). Fig. 3 is a 3-D scatter that cannot be read per point. Arithmetic: `scratch/lit8/pandit_fit.py`.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of chemical marker M-2 formation in mashed potato—a tool to locate cold spots under microwave sterilization" |
| Authors | R. B. Pandit, **J. Tang** (corresponding), G. Mikhaylenko, F. Liu — Washington State University, Pullman |
| Venue | Journal of Food Engineering **76** (2006) 353-361. Received 27 Nov 2004, accepted 19 May 2005, online 9 Aug 2005 |
| DOI as printed | **10.1016/j.jfoodeng.2005.05.032** ✔ matches the brief. **Identity verified.** |

## 1. Why it matters

As for `lau2003_extraction.md`: an independent Ea for ribose + amine → NF (the engine's NF formation rides the lumped
64 kJ/mol), and any evidence of NF loss at 116-131 °C bearing on B61's declared 103.7 kJ/mol hydroxide-decay barrier.

## 2. Methods as they matter to a model (pp. 354-356)

- **Matrix.** 1.5 g D-ribose dissolved in 83.12 g distilled water, mixed with 15.38 g dry potato flakes (Washington
  Potato Co.); 83.12 % moisture. Ribose **99.9 µmol/g sample = 120 mmol/kg water** (mine). Protein 7.3 % (new batch).
  Free amino acids (pre-extracted, 80 % EtOH, Beckman 6300): **Met 1.41, Lys 1.70, His 1.33, Arg 3.70 µmol/g** (only these
  four reported; total free amino acids NOT reported). Unheated potato: no M-2 detected.
- **pH: NOT REPORTED.** Only the framing "weak acidic (pH > 5) environment" (citing Prakash 1997). No buffer.
- **Heating.** ~2 mL paste into glass capillaries (ID 1.75 mm, length 15 cm), flame-sealed; oil bath at **116, 121, 126,
  131 °C**, "upto 60 min"; come-up time **12-14 s**; ice quench. Two replicates per T.
- **Analysis.** 0.16-0.17 g ground in 2 mL 10 mM H2SO4 + 5 mM citric acid; frozen overnight at −20 °C; centrifuged ×2;
  0.45 µm PTFE; Agilent 1100, Bio-Rad fast-acid column, 285 nm. **Calibrated** with commercial M-2 (Givaudan): peak area =
  57 348 × (mg/mL), R² = 0.9983 (Fig. 2, 0-1.2 mg/mL). Yield in mg M-2 per g sample (eq. 2).
- **Model.** dC/dt = k(C∞ − C)ⁿ, C0 = 0; SAS NLIN; Arrhenius about T0 = 396.7 K; R = 1.987 cal/mol/K.

## 3. Numbers

### Table 1 (p. 357), non-linear order and plateau (C0 fixed at 0)

| T (°C) | n | C∞ (mg/g) | R² |
|---|---|---|---|
| 116 | 1.050 ± 0.44 | 0.267 ± 0.04 | 0.989 |
| 121 | 1.0271 ± 0.23 | 0.284 ± 0.018 | 0.990 |
| 126 | 1.062 ± 0.18 | 0.273 ± 0.012 | 0.991 |
| 131 | 1.066 ± 0.40 | 0.268 ± 0.015 | 0.981 |

ANOVA: C∞ not different between temperatures (95 %). Table 2 (r² of linearised plots): first order best at all T
(0.959-0.987).

### Table 3 (p. 359), first-order k (min⁻¹) and Ea

| | modified two-step | two-step multi-linear | **one-step non-linear** | whey gel (Lau 2003) |
|---|---|---|---|---|
| k116 | 0.031 ± 0.014 | 0.029 ± 0.004 | **0.030 ± 0.001** | 0.110 |
| k121 | 0.040 ± 0.007 | 0.0451 ± 0.004 | **0.044 ± 0.002** | 0.152 |
| k126 | 0.056 ± 0.007 | 0.053 ± 0.005 | **0.063 ± 0.003** | 0.207 |
| k131 | 0.096 ± 0.024 | 0.072 ± 0.006 | **0.089 ± 0.005** | 0.281 |
| k0 (123.5 °C) | 0.050 ± 0.048 | 0.048 ± 0.066 | **0.052 ± 0.002** | 0.178 |
| Ea (kcal/mol) | 22.96 ± 2.73 | 19.45 ± 3.79 | **22.23 ± 1.54** | 19.48 |
| Ea (kJ/mol, mine) | 96.1 ± 11.4 | 81.4 ± 15.9 | **93.0 ± 6.4** | 81.5 |
| r² | 0.978 | 0.914 | 0.983 | |

Refit of printed k's (mine): 97.3 / 75.7 / 94.7 kJ/mol (the multi-linear k's give 75.7, not 81.4).

### Fig. 4 (p. 359), ribose level, 121 °C, means of two [D]

| t (min) | 5 | 10 | 20 | 35 | 60 |
|---|---|---|---|---|---|
| 1.5 % ribose (mg/g) | 0.042 | 0.091 | 0.170 | 0.214 | 0.260 |
| 3 % ribose | 0.068 | 0.113 | 0.185 | 0.236 | 0.279 |
| 4.5 % ribose | 0.079 | 0.157 | 0.230 | **0.251** | **0.251** |

### Fig. 5 (p. 359), added L-lysine, 121 °C [D]

| t (min) | 5 | 10 | 20 | 35 | 60 |
|---|---|---|---|---|---|
| original | 0.045 | 0.090 | 0.165 | 0.218 | 0.265 |
| + 0.5 % Lys | 0.130 | 0.207 | 0.320 | 0.408 | 0.455 |
| + 1 % Lys | 0.135 | 0.212 | 0.265 | 0.300 | 0.315 |

Fig. 6: five potato sources at 121 °C, yields not different (P = 0.989). Fig. 7: 0 vs 1 % salt, no effect (P > 0.9988).

### Arithmetic (mine)

1. **Molar yield.** C∞ 0.27 mg/g ÷ 114.1 g/mol = **2.37 µmol/g = 2.4 mol % of the 1.5 % ribose**, 29 % of the sum of the
   four reported amino acids, 1.4× the free lysine. The ribose-independence of the ceiling (Fig. 4: 0.25-0.28 mg/g at 1.5,
   3, 4.5 % ribose) and its rise with added lysine support the authors' "amino acid was the limiting element". The
   first-order "k" is then not a ribose-consumption constant; it is the rate of approach to an amine-set ceiling.
2. **Decay fits (Fig. 4, 121 °C; σ = 0.006 mg/g assumed).** Ribose-limited formation + first-order loss k_d:
   1.5 % ribose χ² 7.2 (k_d = 0) → 6.7 (k_d = 1e-2): **insensitive**. 4.5 % ribose χ² 13.1 (0) → 6.8 (4-6e-3) → 8.8 (1e-2):
   the flat 35→60 min segment is fitted *better* with k_d ≈ 4-6e-3 min⁻¹, but the run is amine-exhausted (formation
   stops), so a flat plateau is also what k_d = 0 + exhaustion gives. **Conclusion: no M-2 loss resolvable; k_d(121 °C)
   ≲ 1e-2 min⁻¹ is all the figure supports.**
3. **Lysine 1 % below lysine 0.5 % at 60 min** (0.315 vs 0.455 mg/g), plateauing from 35 min: either a different ceiling
   (browning competition) or NF consumption by lysine — compare Mikami 2017 ("Lys seemed to repress the colorization";
   NF + Lys). Not resolvable here.
4. **131 °C:** the 3-D Fig. 3 shows the 131 °C series rising monotonically to its last point (highest yield at the longest
   time); no visible decline. Not digitisable per point.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| ribose + potato amines → NF (approach to ceiling) | Ea | **93.0 ± 6.4** (preferred); 96.1 ± 11.4; 81.4 ± 15.9 | kJ/mol | mashed potato, 83 % water, 1.5 % ribose, pH n.r., 116-131 °C | Table 3 | measured |
| same | first-order k | 0.030 / 0.044 / 0.063 / 0.089 | min⁻¹ | 116 / 121 / 126 / 131 °C | Table 3 | measured |
| NF ceiling | C∞ | 0.267-0.284 (n.s. vs T) | mg/g | 116-131 °C | Table 1 | measured, calibrated |
| NF ceiling | molar yield | 2.4 | mol % of ribose | 1.5 % ribose | arithmetic 1 | derived |
| NF ceiling vs ribose | Fig. 4 table | 0.25-0.28 at 60 min | mg/g | 121 °C, 1.5-4.5 % ribose | Fig. 4 [D] | digitised |
| NF vs added Lys | Fig. 5 table | up to 0.455 | mg/g | 121 °C, +0.5 % Lys | Fig. 5 [D] | digitised |
| NF loss | k_d | not resolvable; ≲ 1e-2 | min⁻¹ | 121 °C | arithmetic 2 | weak bound |

## 5. What this can and cannot anchor

**Can:** a second NF-formation Ea (81-96 kJ/mol), with a calibrated, amine-limited NF ceiling of ~2.4 µmol/g; the
T-independence of that ceiling; lysine raising it.
**Cannot:** no pH; free amine pool only partly characterised (four amino acids); food matrix (starch, 83 % water), not a
buffered pot; no NF decay measurement; no sulfur; time resolution 5 min at the start.

## 6. Flags

1. ~2 mL "injected" into a 1.75 mm × 15 cm capillary (max ~0.36 mL) — volume inconsistent, as in Lau 2003.
2. Multi-linear Ea printed 19.45 kcal/mol (81.4 kJ) but its own k's give 75.7 kJ/mol.
3. Table 3 "kcal mol⁻¹ K⁻¹" is a unit misprint.
4. The first-order constant describes approach to an amine-limited ceiling (arithmetic 1), so it is not directly the engine's
   pentose → DPO or DPO → NF constant.
5. Extracts were frozen overnight before analysis (stability not tested here; Kim 1996 says frozen extracts are stable).

## 7. Use in the model

- **(a) The furfural barrier swap.** Pandit's NF Ea (93 ± 6 kJ/mol preferred) and Lau's (81.5 ± 3) bracket the NF route at
  **~80-95 kJ/mol** — above the engine's lumped 64 but well below the 120-147 kJ/mol of furfural formation (Garrett 1969,
  Danon 2014, Schirle-Keller 1992). Two matrices, two labs' worth of amine pools, the same answer: **the NF branch must
  NOT be moved with a furfural-specific barrier increase.** A swap of the shared lumped barrier to ~120-147 would put NF
  formation 30-55 kJ/mol above both measurements. If anything, NF's own route could move modestly upward (64 → ~85).
- **(b) B61's declared 103.7 kJ/mol NF decay.** No measurable M-2 loss over 60 min at 116-131 °C. B61's law at 121 °C
  predicts 2.4e-3 (pH 6.0), 4.3e-3 (6.3), 6.5e-3 (6.5), 2.0e-2 min⁻¹ (pH 7.0); over the 35→60 min window that is 6 %, 10 %,
  15 %, 39 % loss of an already-formed plateau. The 4.5 % ribose plateau (0.251 → 0.251, ±~3 %) is therefore **consistent
  with B61 at pH ≲ 6.3 (or with amine exhaustion + slow loss), and inconsistent with pH ≥ 7**; mashed-potato pH is not
  reported. Like Lau, **Pandit does not test the barrier value**: no loss is measured at two temperatures.
- **Hold-out suitability:** low (matrix, unknown pH, incompletely known amine pool). Its most portable fact is the ratio of
  NF-appearance rates k(131)/k(116) = 2.97 and a T-independent amine-limited ceiling.
