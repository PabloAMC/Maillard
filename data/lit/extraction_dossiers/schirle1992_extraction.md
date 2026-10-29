# Schirle-Keller & Reineccius 1992 — EXTRACTION (glucose 0.5 M + L-cysteine 0.1 M in distilled water, unbuffered, 400 mL stirred Parr reactor, 80 / 100 / 120 / 150 °C, 2-10 h; GC-AED of a CH2Cl2 extract: furfural, furfuryl alcohol, 2-acetylfuran, 5-methylfurfural, DDMP, HMF; zero-order "formation" rates in internal-standard-normalised area units; Arrhenius from three temperatures)

### THE ONLY MULTI-TEMPERATURE FURFURAL TIME COURSE IN A SUGAR + CYSTEINE POT IN THE CORPUS: net furfural accumulation rises **~30× from 100 to 120 °C and ~14-40× from 120 to 150 °C**; the printed apparent Ea is **35.18 kcal/mol = 147 kJ/mol** (my refit of Table 1: 156; 171-189 if the 150 °C rate is per 20-min sampling interval, as Fig. 5 suggests). This is the direction of Kang 2026 and far steeper than the engine's lumped 64 kJ/mol furfural formation (which gives ×2.9 per 100→120 °C). **No furfural loss is visible at any temperature**, but the series are lag-shaped, so the paper does NOT yield a furfural decay constant. **No absolute concentrations, no pH values, no sulfur volatiles.** Hexose, not pentose.

**Source on disk:** `data/articles/schirle1992.pdf` (15 pp., ACS Symp. Ser. publisher scan with OCR text layer). Prose
from `pdftotext -layout`; Table 1 (printed rotated) and Table 2 re-read from a 220-dpi render (`scratch/lit8/schT1.png`);
Figs. 3-6 read from 200-dpi renders (`scratch/lit8/f3c.png`, `f4c.png`, `f5c.png`) and tagged **[D]** with read error.
Arithmetic in `scratch/lit8/schirle_arr.py`, `schirle_arr2.py` (run in the maillard env).

## 0. Identity

| field | value |
|---|---|
| Title | "Reaction Kinetics for the Formation of Oxygen-Containing Heterocyclic Compounds in Model Systems" |
| Authors | J. P. Schirle-Keller, **G. A. Reineccius** — Dept. of Food Science and Nutrition, University of Minnesota, St. Paul |
| Venue | ACS Symposium Series **490** (Flavor Precursors; Teranishi, Takeoka, Güntert eds.), ch. 19, pp. 244-258, 1992. Received 28 Jan 1991 |
| DOI as printed | **10.1021/bk-1992-0490.ch019** ✔ matches the brief. **Identity verified.** |

## 1. Why it matters

`docs/guides/SULFUR_FLUX_MAP.md` sec. 5: furfural formation (`r_tdp_fur`, and the hexose entry `r_glc_fur`, both on the
lumped 64.1 kJ/mol) sits under a decay barrier of 143.5 kJ/mol, so the model's furfural falls 2.5× from 100 to 120 °C
where Kang 2026 measures a 1.7× rise. A barrier swap is being considered; it needs an independent multi-temperature
furfural series in a sugar + thiol pot, and any decay information.

## 2. Methods as they matter to a model (p. 245-246)

- **Pot.** "Glucose (0.5 mole - 36.03 g) and cysteine (0.1 mole - 4.84 g) ... dissolved in 400 ml of distilled water".
  **36.03 g glucose is 0.200 mol and 4.84 g cysteine is 0.040 mol**; in 400 mL these are **0.50 M glucose and 0.10 M
  cysteine** — the "mole" in the text is molarity. Glc:Cys = 5:1. **No buffer.**
- **pH.** "Temperature and pH were recorded prior to reaction"; after cooling each sample "the pH was measured".
  **No pH value is reported anywhere** (initial or final).
- **Vessel.** Parr 4560 series 400 mL stirred pressure reactor (stirrer 60 %), sealed; come-up time 8 min (80 °C) to
  18 min (150 °C); t = 0 when set temperature reached. **Headspace volume not stated.** Each 50 mL sample withdrawn
  through a port and replaced by N2 to restore pressure (so liquid volume falls 400 → ~150 mL over a run and the
  headspace grows; volatile partitioning is uncorrected).
- **Sampling.** 100 and 120 °C: every 2 h to 10 h. 150 °C: every 20 min to 140 min (Fig. 5 shows points only to
  ~1.67 h). 80 °C data were collected but **none were usable** (no MS confirmation).
- **Analysis.** Sample cooled to 25 °C, extracted 3 × 5 mL CH2Cl2 containing octane 10 ppm (internal standard), dried,
  concentrated 300× (150× at 150 °C); GC-AED carbon channel (193 nm) for quantitation, oxygen channel (777 nm) for
  identity; GC-MS for identity. **"Normalized area" = AED carbon area normalised to the internal standard; no response
  factors, no calibration, no recovery correction.** Units "UNA" are arbitrary, and the conversion factor is not given.
- **Kinetics.** Labuza & Kamman program; zero- and first-order fits at each T (Table 1); Arrhenius on zero-order k
  (Table 2, Fig. 6).
- **Replicates.** None stated; single runs.

## 3. Numbers

### Table 1 (p. 255), zero-order formation k in UNA h⁻¹ (verbatim, re-read from the render)

| compound | 100 °C k (r²) | 120 °C k (r²) | 150 °C k (r²) |
|---|---|---|---|
| furfural | 1.692e3 (0.91) | 5.030e4 (0.99) | 6.829e5 (0.89) |
| furfuryl alcohol ("furfurol") | — (blank) | 4.217e4 (0.84) | 6.159e5 (0.91) |
| 2-acetylfuran | 4.151e3 (0.93) | 1.945e5 (0.97) | 1.165e6 (0.96) |
| 5-methylfurfural | 4.866e4 (0.80) | 2.537e5 (0.87) | 3.878e6 (0.98) |
| DDMP ("di(H)-di(OH)-6-Me pyranone") | 1.377e3 (0.73) | 4.238e4 (0.95) | 3.258e5 (0.97) |
| HMF | 4.224e4 (0.89) | 2.922e5 (0.72) | 3.674e6 (0.76) |

First-order fits in the same table (k ≈ 1.8-2.5 h⁻¹ at 100-120 °C, 11-12 h⁻¹ at 150 °C, r² 0.45-0.92) are nearly
identical for every compound, which is the signature of fitting ln(area) of a growing series; the authors reject them
for 120/150 °C ("much better fit with zero order"). **Do not use the first-order k.**

### Table 2 (p. 256), Arrhenius of zero-order k, three temperatures (100, 120, 150 °C)

| compound | slope (K) | intercept | r² | Ea printed (kcal/mol) | Ea (kJ/mol, mine) |
|---|---|---|---|---|---|
| furfural | −1.77e4 | 54.9 | 0.94 | **35.18** | **147.2** |
| 2-acetylfuran | −1.82e4 | 56.8 | 0.87 | 36.20 | 151.3 |
| 5-methylfurfural | −1.86e4 | 59.3 | 0.86 | 37.03 | 154.6 |
| DDMP | −1.54e4 | 49.1 | 0.85 | 30.73 | 128.0 |
| HMF | −1.41e4 | 48.2 | 0.83 | 28.14 | 117.2 |

**Furfuryl alcohol has no Ea** (no 100 °C rate). The text (p. 257) says Ea "ranged from 28 to 37 kcal/mole"; the
abstract says "28 to 33" — the abstract is wrong against its own Table 2.

### Figures, digitised [D]

Fig. 3 (100 °C, carbon channel; read ±0.1e4 UNA). Fig. 4 (120 °C, carbon; ±0.1e5). Fig. 5 (150 °C, carbon; ±0.3e6).

| t (h) | 0 | 2.2 | 4.2 | 6.2 | 8.2 | 10.2 |
|---|---|---|---|---|---|---|
| furfural 100 °C (1e4 UNA) | 0.1 | 0.1 | 0.1 | 0.7 | 0.85 | 1.38 |
| HMF 100 °C (1e4) | 0.1 | 0.2 | 0.7 | 2.0 | 2.25 | 4.5 |
| furfural 120 °C (1e5) | 0.1 | 0.65 | 1.25 | 2.55 | 3.57 | 4.5 |
| furfuryl alcohol 120 °C (1e5) | 0.1 | 0.57 | 0.95 | 2.08 | 2.97 | 4.35 |
| HMF 120 °C (1e5) | 0.1 | 1.53 | 2.28 | 3.12 | **3.67** | **2.70** |

| t (h) | 0 | 0.33 | 0.67 | 1.0 | 1.33 | 1.67 |
|---|---|---|---|---|---|---|
| furfural 150 °C (1e6) | 0.5 | 0.9 | 1.6 | 2.3 | 2.85 | 4.2 |
| HMF 150 °C (1e6) | 0.5 | 3.2 | 5.4 | 6.1 | 10.8 | 22.3 |

### Arithmetic (mine)

1. **Table 1 → Table 2 does not reproduce for every compound.** Least-squares on the three Table-1 values: furfural
   155.6 kJ/mol (printed 147.2), 2-acetylfuran 144.6 (151.3), **5-methylfurfural 115.5 (printed 154.6)**, DDMP 140.9
   (128.0), HMF 117.2 (117.2 ✓). The paper's slopes were fitted to something slightly different from Table 1.
2. **Table 1 vs the figures.** At 100 and 120 °C, furfural, furfuryl alcohol and DDMP Table-1 rates match the figure
   slopes (furfural 100 °C: my fit of Fig. 3 points 2.2-10.2 h gives 1.66e3 UNA/h, intercept −0.40e4 vs printed 1.692e3
   / −4.700e3; 120 °C: 5.01e4 / −6.0e4 vs printed 5.030e4 / **−5.942e3** — the intercept exponent is a misprint).
   2-Acetylfuran, 5-methylfurfural and HMF Table-1 rates at 100 and 120 °C are **~10× the slopes in Figs. 3-4**; at
   150 °C all six are ~1/3 of the Fig.-5 slope in per-hour units (furfural: Fig. 5 slope 2.15-2.35e6 UNA/h vs printed
   6.829e5) — consistent with a regression on **sampling index (20-min steps) instead of hours** at 150 °C. Unresolved;
   it changes the furfural Ea:

   | 150 °C furfural rate used | 3-point Ea (kJ/mol) | 120→150 °C pairwise Ea |
   |---|---|---|
   | 6.829e5 (Table 1 as printed) | 156 | 120 |
   | 2.05e6 (×3, per hour) | 186 | 171 |
   | 2.35e6 (Fig. 5 slope, t > 0) | 189 | 177 |

   100→120 °C pairwise Ea for furfural is 207 kJ/mol, inflated by the 4-h lag at 100 °C (zero-order intercept negative).
3. **Temperature ratios of net furfural accumulation (the shape observable):** k(120)/k(100) = **29.7**; k(150)/k(120) =
   **13.6** (as printed) or **~41** (×3 reading). For comparison: Ea 64.1 kJ/mol gives 2.9 and 4.0; 133 gives 8.9 and 17.9;
   147 gives 11.2 and 24.4.
4. **Amounts at a fixed time** (furfural, normalised area at ~10 h): 1.4e4 (100 °C) vs 4.5e5 (120 °C) → **×32**.
5. **Furfuryl alcohol / furfural** (AED carbon response is nominally per carbon atom, both C5, so area ratio ≈ molar ratio
   in the extract): 0.84 at 120 °C, 0.90 at 150 °C (Table 1 rates); 0.97 at 120 °C / 10 h (Fig. 4). Extraction recovery is
   not corrected (furfuryl alcohol is more water-soluble), so this is likely a **lower bound** on the in-pot ratio.
6. **Furfural decay.** No maximum or plateau in any furfural series. But the 100 °C series has a 4-h lag and the
   120 °C series is near-linear after 2 h, which a consecutive scheme (glucose → intermediate → furfural) produces with
   or without a fast furfural sink: with fast decay, furfural tracks the intermediate (quasi-steady state); without it,
   furfural grows faster than linearly. The engine's `k_fur_decay` at 120 °C (143.5 kJ/mol from 0.38 min⁻¹ at 145 °C:
   0.028 min⁻¹, t½ 25 min) is therefore **neither confirmed nor excluded** by these curves. **No decay estimate.**
7. **HMF loss at 120 °C:** 3.67e5 → 2.70e5 from 8.2 to 10.2 h (−26 %), seen in both carbon and oxygen channels
   (O: 2.5e3 → 2.0e3). If first order and formation stopped, k ≥ ln(3.67/2.70)/120 min = **≥ 2.6e-3 min⁻¹**; single point,
   no replicate. At 150 °C the oxygen channel shows an HMF dip (0.7 → 1.35 h) that the carbon channel does not.

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| net furfural accumulation | apparent Ea (zero order) | **147.2** (printed); 156 (refit); 171-189 (if 150 °C is per 20 min) | kJ/mol | Glc 0.5 M + Cys 0.1 M, water, unbuffered, pH unreported, 100-150 °C | Table 2; arithmetic 1-2 | measured (3 T, single runs) |
| net furfural accumulation | k(120)/k(100) | **29.7** | — | same | Table 1 | within-study ratio (lag-affected at 100 °C) |
| net furfural accumulation | k(150)/k(120) | 13.6 (printed) / ~41 (×3) | — | same | Table 1 / Fig. 5 [D] | within-study ratio, ambiguous |
| HMF | apparent Ea | 117.2 | kJ/mol | same | Table 2 (reproduces) | measured |
| 2-acetylfuran / 5-Me-furfural / DDMP | apparent Ea | 151 / 155 / 128 (printed) | kJ/mol | same | Table 2 (5-Me-furfural does not reproduce) | measured, flagged |
| furfuryl alcohol / furfural | area ratio | 0.84-0.97 | — | 120-150 °C | Table 1, Fig. 4 [D] | ratio, uncorrected recovery |
| HMF loss | k floor | ≥ 2.6e-3 | min⁻¹ | 120 °C, 8-10 h | Fig. 4 [D] | one-point bound |
| absolute concentrations | — | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:**
1. A **unitless temperature-shape test**: the ratio of furfural accumulation between 100, 120 and 150 °C in a sugar +
   cysteine pot. It says the net furfural flux rises steeply with T (apparent Ea ≥ ~120 kJ/mol on any reading of the
   120→150 °C pair), which no network with a 64 kJ/mol formation step under a 143.5 kJ/mol sink can produce.
2. Furfuryl alcohol ≈ furfural in amount at 120-150 °C in a cysteine pot: a large furfural → furfuryl alcohol reduction
   flux exists (engine has no furfuryl alcohol species).
3. HMF apparent Ea 117 kJ/mol (glucose + Cys).

**Cannot:**
- No absolute amounts (arbitrary normalised area), no response factors, no calibration → no mM, no ppb, no yields.
- No pH (initial or final) → the engine's pH factors cannot be set from the paper.
- No sulfur volatiles: thiophenones / "thiofurans" were seen but deferred to "a separate publication" (p. 247).
- No furfural decay constant (arithmetic 6); no furfuryl alcohol Ea.
- Three temperatures, single runs; the 100 °C points are lag-dominated; the 150 °C rate scale is ambiguous.

## 6. Flags

1. **"0.5 mole / 0.1 mole"** are 0.5 M / 0.1 M (the masses are 0.2 / 0.04 mol). `docs/guides/LITERATURE_SEARCH_2026-10_MECHANISMS.md`
   (row for Schirle-Keller) repeats "0.5 mol + 0.1 mol in 400 mL"; read it as molarity.
2. **Abstract Ea "28 to 33 kcal/mole"** contradicts Table 2 (28.14-37.03). The same MECHANISMS doc quotes "Ea 117-138 kJ/mol"
   from the abstract; the furfural value is **147** (Table 2), and there is **no** furfuryl-alcohol Ea.
3. Text p. 247: "data were not available at 100 °C so only three temperatures were used" — means 80 °C.
4. Table 1 exponent misprints and a Table 1 ↔ figure ~10× mismatch for 2-acetylfuran, 5-Me-furfural and HMF at 100/120 °C;
   a ×3 time-unit ambiguity at 150 °C (arithmetic 2). Furfural at 100/120 °C is internally consistent.
5. Sample withdrawal (50 mL each) shrinks the liquid and grows the headspace during a run; uncorrected.
6. "Furfurol" = furfuryl alcohol here (the abstract lists "furfural, furfuryl alcohol"); in older literature "furfurol"
   can mean furfural itself.

## 7. Use in the model

- **(a) The furfural barrier swap.** Supports raising furfural FORMATION well above 64 kJ/mol: the net accumulation ratio
  100→120 °C is ~30× (Ea_app ≥ 147 on the printed table; 120-177 on the 120→150 °C pair alone). This is a net
  (formation − loss) observable, so it cannot apportion between a higher formation barrier and a lower decay barrier, but
  any candidate must make furfural *rise* steeply with T in a sugar + Cys pot; the shipped engine makes it fall. It gives
  **no independent decay barrier**.
- **Hold-out suitability: weak, shape only, and it tests the hexose entry.** Glucose reaches furfural in the engine only via
  `r_glc_fur` (`k_glc_fur`, frozen log10 k(145 °C) = −4.92, B55 Ea 64.08 = the lumped value, UNIDENTIFIED since B2). The
  `HEXOSE ENTRY UNIDENTIFIED` declaration (`src/kinetic_core/engine.py`) covers only MFT/FFT and their products, not
  furfural, but the same argument applies: an absolute furfural number from glucose is a band-floor artefact. What can be
  pre-registered is the **ratio k(120)/k(100) = 29.7** (and, with the 150 °C ambiguity declared, k(150)/k(120) ∈ [13.6, 41])
  of furfural accumulation, with initial pH declared as an assumption (a 0.1 M cysteine solution is near its pI; the
  paper gives no number). If the swap is applied to the *shared* lumped barrier (which `k_glc_fur` also carries), this
  ratio does test it; it does not test `r_tdp_fur` (pentose route).
- **Furfuryl alcohol** ≈ furfural (0.84-0.97 by area) in a cysteine pot at 120-150 °C is a reductive furfural sink of the
  same order as furfural itself — a candidate identity for part of the "large unidentified" `k_fur_decay` sink (see
  `rizzi2010_extraction.md` for the reductant pool). Not a rate.
- **HMF:** Ea_app 117 kJ/mol and a ≥ 2.6e-3 min⁻¹ late loss at 120 °C, both in a cysteine pot — a cross-check for the
  trunk's HMF branch, not an anchor.
