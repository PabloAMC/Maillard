# Cui, Xin, Ji, Wang, Sun, Yu, Wei, Wang, Zhang & Wang 2026 — EXTRACTION
DOI: 10.1016/j.lwt.2026.120020 (printed on p. 1; Crossref 2026-11-01: "Formation of volatile compounds and glycation products in thermally processed chestnut protein–sugar systems: Effects of water activity and temperature", Cui, Xin, Ji, Wang, Sun, Yu, Wei, Wang, Zhang, Wang, LWT 258, article 120020, 2026-10. **Matches the brief.**)
### Chestnut protein isolate (CPI) alone, + maltose (CPI-M) or + sucrose (CPI-S), protein : sugar 1 : 1.5 (w/w), freeze-dried, equilibrated over saturated salts to a_w 0.53 / 0.75 / 0.97 (AquaLab-verified BEFORE heating), then ONE 20-min heating in a forced-convection oven at 100 / 150 / 200 °C. End points: colour, SDS-PAGE / PAS, FTIR, soluble protein, 3-DG / GO / MGO and "protein-bound" CML / CEL (LC-MS/MS), 36 volatiles (HS-SPME-GC-MS, semi-quantitative). **The only a_w x T factorial on disk, but it is a 3 x 3 grid of single 20-min end points: no time course, no rate, no pH, a_w not held or measured during heating.**

**Read 2026-11-01 from `data/articles/cui2026.pdf`** (14 pp., born-digital, clean text layer). Prose and Tables 1-3
from `pdftotext -layout`; **Fig. 4 (the 3 x 5 grid of 3-DG / GO / MGO / CML / CEL) digitised from 300 dpi crops of
p. 9** and Fig. 5 (volatile heatmap, PCA, Venn) checked against a 200 dpi render of p. 10 (scratchpad
`ing5/cuir*-09.png`, `cuiq*-09.png`, `cui-10.png`). Read error on Fig. 4: +-0.3 ug/g (GO), +-1.5 (MGO), +-3-5 (CML),
+-4-8 (CEL), +-100 (3-DG above the axis break); 3-DG values BELOW the axis break (~0-500) are not readable. The
supplement (Table S1 volatile concentrations, Table S2, Fig. S1) is not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Tianjin Key Laboratory of Food Science and Health, School of Medicine, Nankai University (Wang
Shuo; Wang Ruican), with Hebei Agricultural University (Wang Xianghong). **Not in `lab_ids.yml`: NEW** (proposed id
`nankai_wang`). [!] Not the Cui of the Ho network (Cui Heping, Jiangnan; `cui2022_extraction.md`): different first
author (Cui Haoxin), different institute; do not merge on the surname.

## 1. What this paper answers for the model

| model need | what Cui 2026 gives | tag |
|---|---|---|
| a_w x T dependence of dicarbonyls and CML / CEL on a PLANT protein | 3 a_w x 3 T x 3 matrices, one 20-min end point each (Fig. 4); a_w set by salts and verified by AquaLab before heating | [M] / [F] |
| a RATE vs a_w (is there a maximum near a_w 0.7?) | **No rate.** One heating time; the grid is end-point yield vs a_w and T. Fig. 4 is captioned "kinetic profiles" but its x-axis is temperature [!] | [M] |
| shape vs a_w | **Depends on the marker and on T** (sec. 3): dicarbonyls rise monotonically with a_w at 150 °C (no maximum); CML and CEL peak at a_w 0.75 at 150 °C for both sugars; at 100 °C the maltose AGEs peak at 0.75, the sucrose ones at 0.97 | [F] / [Z] |
| a_w during heating | **Not measured, not held.** Equilibrated in open Petri dishes; heating vessel and closure not stated; no mass loss reported | [!] |
| pH | **Not measured** anywhere in the heating study | — |
| O2 | forced-convection oven, air; not varied | [M] |
| protein-bound amine identity | chestnut protein isolate, 92.62 % purity, with 4.31 % residual carbohydrate and 141.05 mg GAE/100 g DM phenolics; no lysine or free-amine measurement | [M] |
| volatiles vs a_w and T | **Only qualitatively in the main text** (z-scored heatmap, Fig. 5A) plus three printed ranges; the per-cell numbers are in Table S1 (not on disk) | [M] / [F] |

**Answers.**
1. **Dicarbonyls vs a_w (CPI-M / CPI-S, ug/g protein, 20 min):** at 150 °C, MGO 44 / 47 / 57 and 32 / 31 / 54; GO
   10.6 / 11.3 / 19.8 and 6.6 / 7.2 / 19.5 at a_w 0.53 / 0.75 / 0.97 — **monotone rise, no intermediate maximum** [F].
   3-DG is highest at 100 °C and falls with T; at 100 °C CPI-M 2750 / 4000 / 5500, CPI-S below the break / below the
   break / 8300 [F]. 3-DG, MGO and GO all fall from 150 to 200 °C (consumed).
2. **CML / CEL vs a_w:** at 150 °C, CML 55 / 63 / 40 (M) and 94 / **137.83** / 66 (S); CEL 81 / 98 / 72 (M) and 110 /
   **238.48** / 130 (S) — **a maximum at a_w 0.75** (x1.2-1.6 for M, x1.5-2.2 for S over the flanks) [F; bold printed].
   At 100 °C: CML 27 / 50 / 44 (M), 20 / 17 / 47 (S); CEL 58 / 84 / 61 (M), 24 / 11 / 109 (S). At 200 °C everything
   with sugar falls back (CML 18-33, CEL 35-67). "CEL/CML ratio of approximately 1.7:1" [V].
3. **Is this a rate maximum near a_w 0.7?** **Not demonstrable** (confidence ~85 %). An end-point maximum of an
   intermediate-like product (CML / CEL fall from 150 to 200 °C, so they are consumed) can arise from faster
   consumption at high a_w as easily as from slower formation. The engine-relevant reading is narrower: at 150 °C
   the MGO pool rises with a_w while CML / CEL peak at 0.75 [Z].
4. **Protein-bound-amine relevance:** a plant-protein glycation end-point set with CML / CEL on the protein (hydrolysis
   and any NaBH4 reduction NOT described; method cited to Wang et al. 2024) [!]. The sugar-free CPI makes
   **81-158 ug/g CEL at 200 °C**, MORE than either sugar system at 200 °C [F] — a non-sugar CEL source (authors: protein
   oxidation or residual sugars) [P] or an analytical artefact; either way CEL here is not a clean sugar-glycation
   marker at 200 °C.
5. **Volatiles:** 100 °C "volatile generation was negligible across all matrices" [V]; 150 °C CPI-M dominated by maltol
   (2156-5463 ug/kg) and furfurals [M]; 150 °C CPI-S at a_w 0.53-0.75 gives the alkylpyrazines (heatmap hot cells for
   2,5-dimethyl-, 2,3,5-trimethyl-, 3-ethyl-2,5-dimethyl-, 3,5-diethyl-2-methylpyrazine in the CPI-S-0.53-150 and
   -0.75-150 columns) [F, z-score only]; 200 °C CPI-M 5-methylfurfural 6396-10,223 and maltol 750-6295 ug/kg [M];
   200 °C CPI-S "near-total volatile attenuation" [V]; sugar-free CPI at 200 °C gives C6-C9 aldehydes and
   benzaldehyde. **No volatile number can be placed on an a_w x T cell without Table S1.**

## 2. Conditions

- **Materials** [M]: CPI by alkaline extraction (pH 8, 4 °C, 2 h) / isoelectric precipitation (pH 4.32), neutralised
  to pH 7, dialysed 72 h; yield 30.20 +- 2.47 %. Maltose and sucrose >= 99 %.
- **Mixtures** [M]: CPI-M and CPI-S solutions 12.0 mg/mL, protein : sugar 1 : 1.5 (w/w), stirred 1 h at 25 °C,
  hydrated 14 h at 4 °C, freeze-dried. No buffer.
- **a_w** [M]: ~500 mg in open Petri dishes over saturated Mg(NO3)2·6H2O (0.53), NaCl (0.75), K2SO4 (0.97), 7-10 d;
  verified with an AquaLab-TLE meter. Moisture contents not printed.
- **Heating** [M]: forced-convection oven, 100 / 150 / 200 °C, 20 min; cooled, sealed, -20 °C. Vessel, sample depth
  and come-up time not stated [!]. Open dishes at 150-200 °C for 20 min would lose most of their water: the a_w
  label is the starting a_w, not the reaction a_w [Z].
- **Analysis** [M]: 3-DG / GO / MGO by LC-MS/MS (MRM 235/199, 131/77, 145/77 — quinoxaline masses, so OPD
  derivatisation is implied though not described [Z]); CML 205.2/84.02, CEL 218.9/84.0 on HILIC; recoveries
  85.5-100.4 %; volatiles HS-SPME (DVB/CAR/PDMS) GC-MS, NIST + RI, semi-quantified against 2-methyl-3-heptanone (no
  response factors). n = 3; three-way ANOVA.

## 3. Quantitative results

**Table 1** (whole chestnuts, g/100 g DM, 20 min) [M]: unheated sucrose 11.33; 100 °C sucrose 12.47, maltose 4.77;
150 °C sucrose 14.67, maltose 2.23; 200 °C sucrose 12.27, maltose < 0.5 (LOD). Glucose, fructose, lactose never
detected.

**Fig. 4** (ug/g protein, 20 min; M = CPI-M, S = CPI-S, P = CPI alone; 100 / 150 / 200 °C) [F unless printed]:

| a_w | marker | M | S | P |
|---|---|---|---|---|
| 0.53 | 3-DG | 2750 / 1950 / < 500 | below break x3 | ~0 |
| 0.53 | GO | 4.6 / 10.6 / 4.4 | 3.7 / 6.6 / 2.6 | 2.3 / 1.8 / 1.6 |
| 0.53 | MGO | 13.7 / 44.2 / 20.2 | 8.6 / 31.9 / 11.6 | 4.4 / 5.3 / 2.9 |
| 0.53 | CML | 27 / 55 / 31 | 20 / 94 / 32 | 7 / 12 / 17 |
| 0.53 | CEL | 58 / 81 / 60 | 24 / 110 / 64 | 6 / 17 / 87 |
| 0.75 | 3-DG | 4000 / 2100 / < 500 | below break x3 | ~0 |
| 0.75 | GO | 4.3 / 11.3 / 4.6 | 2.7 / 7.2 / 3.0 | 2.0 / 2.2 / 1.5 |
| 0.75 | MGO | 15.8 / 46.9 / 21.9 | 6.9 / 30.6 / 14.4 | 4.4 / 5.2 / 3.1 |
| 0.75 | CML | 50 / 63 / 18 | 17 / **137.83** / 29 | 4 / 18 / 29 |
| 0.75 | CEL | 84 / 98 / 35 | 11 / **238.48** / 50 | 5 / 27 / 158 |
| 0.97 | 3-DG | 5500 / 2300 / < 500 | 8300 / 3750 / < 500 | ~0 |
| 0.97 | GO | 7.8 / 19.8 / 3.5 | 10.2 / 19.5 / 1.9 | 4.7 / 2.1 / 1.5 |
| 0.97 | MGO | 20.5 / 57 / 18 | 31.3 / 53.6 / 6.9 | 4.7 / 5.4 / 4.0 |
| 0.97 | CML | 44 / 40 / 26 | 47 / 66 / 33 | 9 / 14 / 16 |
| 0.97 | CEL | 61 / 72 / 45 | 109 / 130 / 67 | 10 / 17 / 81 |

Unheated CPI-M already holds **209.23 ug/g 3-DG**, "approximately 47-fold higher than CPI-S" [M] (formed during
co-lyophilisation). [Z] 3-DG yield: 2750 ug/g protein at a_w 0.53 / 100 °C = 0.39 mol % of the maltose; 8300 ug/g
(CPI-S, a_w 0.97, 100 °C) = 1.2 mol % of the sucrose.

**Table 2** (colour, Delta E vs unheated) [M]: 100 °C: M 10.40 / 12.27 / 30.66, S 12.54 / 12.53 / 39.11, P 13.33 /
11.97 / 8.89 at a_w 0.53 / 0.75 / 0.97; 150 °C: M 39.04 / 39.25 / 41.34, S 37.18 / 35.85 / 45.37; 200 °C: M 51.18 /
55.13 / 53.51, S 53.72 / 53.01 / 53.48, P 30.81 / 30.94 / 33.91.

**Table 3** (soluble protein, % of total) [M]: unheated P 44.34, M 59.55, S 51.43; 100 °C M 49.23 / 44.72 / 30.86,
S 36.22 / 33.36 / 27.00, P 31.44 / 28.57 / 22.76; 150 °C 1.40-6.68; 200 °C 0.01-0.62.

Correlations (Spearman, Fig. 6A) [M]: MGO-CML r = 0.828; Delta E-SPC r = -0.821; MGO with furans / pyrroles r =
0.609-0.737.

## 4. Mechanism evidence and its kind

Within-study factors only (sugar type x a_w x T, one time point), plus correlations and PCA. No labelling, no fed
intermediate, no mass balance, no kinetics. "Sucrose inversion at high a_w" is [P]: glucose and fructose were never
detected (Table 1 in chestnuts; not measured in the model mixtures). "MGO as the link between glycation and aroma"
is a correlation [Q]. Laboratory: Nankai (Wang Shuo), NEW.

## 5. Bearing on open questions

- **Atlas sec. 1.4, "intermediate a_w with a measured a_w"** — **PARTLY filled** for F05 / F08-F10 end points: a_w is
  measured (before heating) at 0.53 / 0.75 / 0.97 and crossed with 100 / 150 / 200 °C. Not a rate, and the reaction
  a_w is unknown (open oven, no mass loss). Confidence it can anchor an a_w TERM in the engine: low (~0.2).
- **Atlas sec. 1.4, "plant-protein-bound amines"** — **PARTLY filled**: CML / CEL on a plant protein (chestnut) with
  two disaccharides; no lysine, no furosine, no rate. LEAVES the F05 unknown "pea / soy glycation KINETICS" unchanged.
- **F05 unknown "high T (100-140 C) at intermediate a_w with a measured a_w"** — the 100 °C row and the 150 °C row
  are the nearest cells now on disk (with `wolf1977_extraction.md`, 125 °C, a_w 0.65, sealed).
- **F05 engine mismatch (d), no a_w term on glycation** — SUPPORTS that water matters for protein glycation end
  points (x1.2-2.2 between a_w levels at 150 °C), but gives no functional form; the direction differs by marker and
  temperature, so it cannot choose between "maximum near 0.5-0.75" (Warmbier) and "flat then falls at 0.98"
  (Pereyra).
- **F05 CEL route (MGO + bound Lys vs Amadori -> CEL)** — LEAVES UNCHANGED. MGO and CEL co-vary across cells
  (Fig. 6A), consistent with either route; CEL from sugar-free CPI at 200 °C argues that some CEL here is not from the
  added sugar at all.
- **CF-072 (lysine return)** — LEAVES UNCHANGED (no lysine measured).
- **F30 (lipid crossover)** — LEAVES UNCHANGED: the sugar-free CPI's C6-C9 aldehydes at 200 °C come from defatted
  protein with unstated residual lipid; no lipid-present-vs-absent contrast.

## 6. Disagreements found

1. **Water dependence at the wet end (NEW).** Cui: at 100 °C browning, 3-DG and the sucrose AGEs are HIGHEST at a_w
   0.97 (Delta E 30.66 / 39.11 vs ~10-12.5 at 0.53-0.75). Pereyra (skim milk powder, available lysine, 37-60 °C):
   ~4x SLOWER at a_w 0.98 than at 0.31-0.71; Warmbier (casein-glucose-glycerol): maximum at a_w ~0.5. Candidate causes:
   (i) sucrose needs water to invert (authors' reading) and maltose's 3-DG also rises with a_w; (ii) the reaction a_w
   is not the label: an open sample at 100 °C dries during 20 min, so "0.97" passes through intermediate a_w while
   "0.53" goes glassy-dry; (iii) end-point vs initial rate; (iv) 100 °C vs 37-60 °C (above Tg for every cell);
   (v) different markers (3-DG / colour vs available lysine).
2. **Internal inconsistency on chestnut a_w** [!]: sec. 2.2.4 says heating whole chestnuts at 100 / 150 / 200 °C
   gives a_w ~0.91 / 0.84 / 0.55; sec. 3.1 says 0.97 (raw) -> 0.75 and 0.53. Not both can be the same measurement.
3. **CEL in protein without sugar exceeds CEL with sugar at 200 °C** [!] (CPI 87 / 158 / 81 vs CPI-M 60 / 35 / 45,
   CPI-S 64 / 50 / 67 ug/g protein). No other dossier on disk reports sugar-free CEL above the sugar arm; causes:
   carbonyls from protein / phenolic oxidation in a near-insoluble 200 °C residue (SPC 0.57-0.62 %), the residual
   4.31 % carbohydrate, or an extraction / hydrolysis artefact on the insoluble sugar-containing residues (SPC
   0.01-0.06 %) that under-recovers their CEL.

## 7. Ledger-ready measurements

- Fig. 4 MGO, GO, CML, CEL for CPI-M and CPI-S at the 9 a_w x T cells (two readers; anchors 137.83 / 238.48 ug/g
  printed) — the only a_w x T grid of glycation end points on a plant protein. Key as end points at 20 min with the
  starting a_w, NOT as rates, and flag the reaction-a_w caveat.
- Fig. 4 3-DG above the axis break (CPI-M all a_w at 100 / 150 °C; CPI-S a_w 0.97) with the printed 209.23 ug/g
  unheated anchor.
- Table 2 Delta E (printed) as a browning end point per cell.
- **Request Table S1** (volatile concentrations per cell) before any volatile is keyed.
