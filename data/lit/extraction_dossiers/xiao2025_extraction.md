# Xiao, Flory, Alavi & Li 2025 — EXTRACTION (six soy / pea / gluten blends before and after low-moisture twin-screw extrusion: free NH2 by TNBS, free SH by Ellman's, H0, IVPD, functionality)
### Blends, not isolates: the only paper on disk with a paired before / after-extrusion free-amino AND free-thiol measurement on the same protein material, on a TNBS scale that is again 10-18x the lysine content.

**Source on disk:** `data/articles/Xiao2025.pdf` (owner's download, 2026-09-08). Read from the pdftotext
layer in the scratchpad (`articles/Xiao2025.txt`, 1757 lines); Tables 2, 4 and 5 have a clean text
layer (values and significance letters interleaved line by line) and are re-typed in full below.
Fig. 1 (solubility in six buffers), Fig. 2 (foaming) and Fig. 3 (PCA) are figure-only. The
**supplementary material (Tables S1-S4, Fig. S1) is NOT on disk**; it holds the per-ingredient
values (SPI, two SPCs, soy flour, PPI, gluten) that would have been the isolate rows. Repo status
before this dossier: soy and pea isolates carry native free-thiol and disulfide densities in
`protein_matrices.yml`; no amine density; no after-heat or after-extrusion factor.

## 0. Identity

| field | value |
|---|---|
| Title | "Physicochemical and functional properties of plant proteins before and after extrusion texturization" |
| Authors | Ruoshi Xiao, Jenna Flory, Sajid Alavi, Yonghui Li (Kansas State University, Grain Science and Industry) |
| Venue | Food Hydrocolloids 163 (2025) 111119; received 30 Sep 2024, revised 16 Jan 2025, accepted 20 Jan 2025, online 21 Jan 2025; PII S0268005X25000797 |
| DOI | 10.1016/j.foodhyd.2025.111119 |
| Naming | CS = cold-swelling proteins (SPI, SPC "Arcon S", PPI); HS = heat-swelling (SPC "Arcon F", wheat gluten); "x % CS" = a raw-material blend whose CS share is x %; (R) = raw blend, (T) = texturized extrudate (TVP). Table 1 is a literature table; Table 2 = FTIR; Table 3 = which buffer breaks which bond; Table 4 = free NH2, free SH, H0, IVPD; Table 5 = WHC, OHC, EAI, ESI, LGC. §3.2 refers to "the test of free -SH group determination (Table 2)" — that is Table 4. |
| Companions | Flory, Xiao, Li, Dogan, Talavera & Alavi 2023, Foods 12, 3232 (the same six blends: formulations, extrusion conditions, texture) — NOT on disk. Method parents: Großmann 2021 (TNBS), Hao 2022 (SH), Tang 2021 (SDS-binding H0), Tinus 2012 (IVPD). Same lab and same two assays as Xiao 2024 (book chapter, `xiao2024_extraction.md`) and Shen 2022. |

## 1. Why it matters

The matrix layer has native soy and pea thiol / disulfide densities and nothing for the amine
pool, and nothing that says how any density moves when the isolate is cooked or extruded. This
paper measures free amino groups (TNBS) and free sulfhydryl (Ellman's after 8 M urea) on the same
six protein blends before and after pilot-scale low-moisture extrusion. Two things carry over:
(i) the **within-blend after/before ratios** — free SH falls to 0.41-0.92x (mean 0.68x), the TNBS
free-amino number to 0.34-0.71x, SDS-binding hydrophobicity to 0.58-0.86x — which are the first
extrusion factors on disk for soy/pea material; (ii) another confirmation that this lab's TNBS scale
(3.5-6.9 mmol/g protein for blends of ~71 % protein, 3.07 for gluten) is an order of magnitude
above any lysine content and cannot feed the amine pool. What does NOT carry over: no row is a pure
isolate (the per-ingredient data are in the supplement, not on disk), the blend recipes are only
described qualitatively, and disulfide / total SH were not reported for the blends.

## 2. Methods as they matter to a model

- **Materials (§2.1):** "six treatments with varying cold-swelling protein ratios: 0% CS, 30% CS,
  40% CS, 50% CS, 60% CS, and 90% CS. These samples consisted of blends of soybean protein isolate
  (SPI), two variants of soybean protein concentrate (Arcon F and Arcon S SPC), soybean flour, pea
  protein isolate (PPI), and wheat gluten (Gluten)." "the overall protein content across all samples
  remained relatively consistent (around 71.48%)." Exact recipes are in Flory et al. 2023 (not on
  disk). What this paper's text lets one reconstruct: **0 % and 40 % CS are the only blends with
  wheat gluten** (§3.2, §3.3); **50 % and 90 % CS are the only blends with PPI** (§3.3, §3.5);
  **30 % and 60 % CS contain only soybean protein and soy flour** (§3.7.2). Ingredients, blends and
  extrudates were "ground into fine powder using a coffee grinder" and kept at 4 °C.
- **Extrusion:** "pilot-scale twin-screw extruder" (Abstract), low-moisture TVP route (Introduction:
  30-40 % moisture); temperature, screw speed, moisture, die and drying are NOT in this paper
  (Flory 2023).
- **Free amino groups (§2.4, TNBS):** "0.5 mL of a sample solution, with a concentration of 4 mg/mL
  and dissolved in a 1% (w/v) SDS solution, was combined with 4 mL of 0.2125 M phosphate buffer (pH
  8.2 ...) and 4 mL of 0.1% (v/v) 2,4,6-trinitrobenzene sulfonic acid (TNBS). This suspension was
  then incubated in a 50 °C water bath at 200 rpm for 1 h in darkness. Following the hour-long
  reaction, 8 mL of 0.1 N HCl was introduced ... stored at room temperature for 30 min ...
  centrifugation at 8000 g for 10 min, and its absorbance was measured ... at 340 nm. Additionally,
  L-leucine solutions were prepared following the same procedure, but with varying concentrations
  (ranging from 0 to 2.4 mM), to establish a standard curve." The conversion to "mmol/g protein" is
  not printed. **Calibration ceiling:** 2.4 mM in a 0.5 mL aliquot = 1.2 µmol NH2; with 2 mg of
  sample in the aliquot the curve tops out at 0.6 mmol/g sample = 0.84 mmol/g protein at 71.48 %
  protein. Table 4's raw blends (3.47-6.90 mmol/g protein = 2.5-4.9 mmol/g sample) are 4-8x above
  the top standard.
- **Free sulfhydryl (§2.5, Hao 2022):** "A 75 mg sample was dispersed in 10 mL of Tris-Gly-urea
  buffer (containing 0.086 mol/L Tris, 0.09 mol/L glycine, 0.004 mol/L EDTA, and 8 mol/L urea) and
  shaken overnight. For free sulfhydryl content determination, 1 mL sample was mixed with 4 mL
  Tris-Gly buffer and 0.05 mL Ellman's reagent (4 mg/mL), shaken for 15 min in the dark, centrifuged
  (8000 × g, 8 min), and absorbance measured at 412 nm." Total SH (β-mercaptoethanol, 12 % TCA) and
  S-S = (total − free)/2 are described, with SH (µmol/g) = 73.53 × A412 × D / C, C in mg/mL, D = 5
  (free) or 10 (total) — identical to the book chapter — **but no total-SH or S-S value for any
  blend is reported in the paper**; Table 4 has free SH only. C is the powder concentration
  (7.5 mg/mL), so Eq. 1 yields µmol per g of powder; Table 4 says "µmol/g protein" and whether a
  division by 0.7148 was applied is not stated.
- **Surface hydrophobicity (§2.5, Tang 2021):** SDS-binding method: 10 mg protein sample + 40 mL
  0.1 mmol/L SDS, 1 h; dialysis 48 h (3.5 kDa); methylene blue / chloroform partition, A655;
  reported as µg SDS bound per mg protein. Not an ANS slope; not comparable in scale to Shen 2022
  or Chen 2022.
- **IVPD (§2.6):** pH-drop, trypsin / chymotrypsin / protease, 6.25 mg protein/mL, pH 8.00, 37 °C,
  10 min; IVPD % = 65.66 + 18.10 × ΔpH(10 min) (Tinus 2012).
- **Solubility (§2.3):** 200 mg in 10 mL of six extractants (IEF = 8 M urea + 50 mM DTT + 2 % SDS +
  2 M thiourea + 2 % CHAPS in phosphate; and the drop-one variants; PB = 100 mM phosphate pH 7.5),
  2 h, 8000 g 15 min, Bradford; extrudates only (Fig. 1, figure-only; §3.2 prints some numbers).
- **Functional (§2.7):** WHC/OHC 0.25 g in 7.5 mL, 4500 g 15 min; EAI/ESI turbidimetric at 500 nm;
  FC/FS per Shen 2021; LGC per Shen & Li 2021.
- **Statistics (§2.8):** "at least duplicate", mean ± SD; Table notes say n = 2; Tukey P < 0.05.

## 3. Tables re-typed

### Table 2. "Protein secondary structures from FTIR" (%; mean ± SD, n = 2)

| sample | α-helix | β-sheet | β-turn | random coil |
|---|---:|---:|---:|---:|
| 0 % CS (R) | 33.20 ± 3.09 | 60.67 ± 3.32 | 6.13 ± 0.23 | 0.00 ± 0.00 |
| 30 % CS (R) | 27.65 ± 0.00 | 57.51 ± 0.00 | 3.24 ± 0.00 | 11.60 ± 0.00 |
| 40 % CS (R) | 39.86 ± 2.56 | 57.13 ± 1.38 | 3.01 ± 1.19 | 0.00 ± 0.00 |
| 50 % CS (R) | 22.45 ± 0.49 | 63.58 ± 0.50 | 3.26 ± 0.23 | 10.72 ± 0.21 |
| 60 % CS (R) | 23.05 ± 0.92 | 51.23 ± 0.05 | 6.79 ± 0.31 | 18.93 ± 1.18 |
| 90 % CS (R) | 25.98 ± 4.82 | 57.74 ± 0.90 | 5.86 ± 3.72 | 10.42 ± 0.20 |
| 0 % CS (T) | 22.95 ± 1.03 | 71.84 ± 0.11 | 5.21 ± 0.91 | 0.00 ± 0.00 |
| 30 % CS (T) | 20.52 ± 1.27 | 71.86 ± 1.15 | 7.62 ± 0.11 | 0.00 ± 0.00 |
| 40 % CS (T) | 11.40 ± 0.20 | 57.91 ± 0.90 | 5.43 ± 0.24 | 25.26 ± 0.94 |
| 50 % CS (T) | 25.81 ± 3.20 | 55.18 ± 0.39 | 3.64 ± 0.09 | 15.37 ± 2.90 |
| 60 % CS (T) | 23.87 ± 1.72 | 62.14 ± 3.06 | 13.98 ± 1.33 | 0.00 ± 0.00 |
| 90 % CS (T) | 30.26 ± 3.34 | 65.56 ± 2.21 | 4.19 ± 1.12 | 0.00 ± 0.00 |

(Rows sum to ~100 %. Text (§3.1): wheat gluten alone had 46.82 % α-helix, Table S1.)

### Table 4. "Protein free amino content, free sulfhydryl content, surface hydrophobicity, and in vitro protein digestibility (IVPD)" (mean ± SD, n = 2; letters = Tukey groups within a column)

| sample | free amino (mmol/g protein) | free SH (µmol/g protein) | H0 (µg SDS/mg protein) | IVPD (%) |
|---|---:|---:|---:|---:|
| 0 % CS (R) | 3.47 ± 0.05 g | 2.21 ± 0.10 de | 59.73 ± 1.96 c | 85.21 ± 0.26 g |
| 30 % CS (R) | 4.68 ± 0.02 e | 2.66 ± 0.12 bcd | 64.14 ± 1.20 b | 87.38 ± 0.26 e |
| 40 % CS (R) | 4.88 ± 0.02 d | 2.61 ± 0.02 cd | 64.28 ± 0.35 b | 87.47 ± 0.13 e |
| 50 % CS (R) | 5.53 ± 0.03 b | 3.12 ± 0.24 ab | 54.83 ± 0.97 d | 89.37 ± 0.26 c |
| 60 % CS (R) | 5.29 ± 0.05 c | 2.83 ± 0.03 abc | 69.05 ± 1.67 a | 88.83 ± 0.26 d |
| 90 % CS (R) | 6.90 ± 0.05 a | 3.24 ± 0.04 a | 57.27 ± 0.55 cd | 91.18 ± 0.51 a |
| 0 % CS (T) | 1.75 ± 0.02 i | 2.03 ± 0.02 e | 35.80 ± 1.72 g | 86.11 ± 0.00 f |
| 30 % CS (T) | 3.11 ± 0.05 h | 1.99 ± 0.03 e | 38.41 ± 1.35 g | 89.64 ± 0.13 c |
| 40 % CS (T) | 1.67 ± 0.03 i | 2.06 ± 0.06 e | 37.14 ± 1.44 g | 86.20 ± 0.13 f |
| 50 % CS (T) | 3.90 ± 0.05 f | 2.03 ± 0.14 e | 42.22 ± 0.83 f | 90.37 ± 0.13 b |
| 60 % CS (T) | 3.19 ± 0.03 h | 1.96 ± 0.03 e | 42.00 ± 0.14 f | 89.46 ± 0.13 c |
| 90 % CS (T) | 3.85 ± 0.06 f | 1.32 ± 0.04 f | 49.10 ± 1.29 e | 89.82 ± 0.13 c |

Text-only numbers that belong with this table (from the supplement, quoted in §3.3-3.5): wheat
gluten free amino **3.07 mmol/g protein** ("the lowest ... among the ingredients"); PPI "had the
highest free amino content among all the raw ingredients" (value not printed); PPI H0 **50.69 µg
SDS/mg protein** (lowest of the ingredients); the 60 % CS ingredients' H0 ranged 53.28-77.65.
Means quoted in the text: free SH 2.78 (R) → 1.90 (T) µmol/g protein; H0 61.55 (R) → 40.78 (T).

### Table 5. "Water holding capacity (WHC), oil holding capacity (OHC), emulsifying properties, and least gelation concentration (LGC)" (n = 2)

| sample | WHC (g water/g material) | OHC (g oil/g material) | EAI (m²/g) | ESI (min) | LGC (%) |
|---|---:|---:|---:|---:|---:|
| 0 % CS (R) | 1.62 ± 0.01 | 0.82 ± 0.00 | 8.44 ± 0.13 | 15.99 ± 0.43 | 18 |
| 30 % CS (R) | 2.77 ± 0.03 | 0.85 ± 0.00 | 10.82 ± 0.26 | 19.10 ± 0.94 | 20 |
| 40 % CS (R) | 3.23 ± 0.12 | 0.94 ± 0.01 | 10.05 ± 0.26 | 20.47 ± 0.18 | 20 |
| 50 % CS (R) | 2.32 ± 0.03 | 0.59 ± 0.02 | 11.41 ± 0.05 | 21.77 ± 0.88 | 17 |
| 60 % CS (R) | 3.65 ± 0.03 | 0.96 ± 0.01 | 9.72 ± 0.02 | 22.84 ± 0.80 | 17 |
| 90 % CS (R) | 4.69 ± 0.02 | 0.93 ± 0.05 | 10.01 ± 0.47 | 21.34 ± 0.64 | 18 |
| 0 % CS (T) | 2.55 ± 0.03 | 0.89 ± 0.04 | 1.59 ± 0.02 | 58.81 ± 2.27 | 19 |
| 30 % CS (T) | 3.01 ± 0.01 | 1.04 ± 0.00 | 6.50 ± 0.10 | 14.78 ± 0.58 | 20 |
| 40 % CS (T) | 2.76 ± 0.01 | 0.84 ± 0.00 | 1.84 ± 0.20 | 40.60 ± 0.11 | 20 |
| 50 % CS (T) | 1.98 ± 0.02 | 1.03 ± 0.00 | 3.73 ± 0.18 | 32.92 ± 0.91 | > 20 |
| 60 % CS (T) | 2.98 ± 0.02 | 0.89 ± 0.00 | 4.54 ± 0.45 | 21.19 ± 2.00 | > 20 |
| 90 % CS (T) | 2.97 ± 0.02 | 1.19 ± 0.00 | 2.75 ± 0.10 | 36.18 ± 2.07 | 19 |

(Column heads say "g/g material"; Eqs. 4-5 say "g/g protein". Same numbers.)

### Solubility numbers printed in §3.2 (Fig. 1 itself is figure-only)

Raw blends in water 20.04-46.20 % (Flory 2023). Extrudates: PB 0.64-1.74 % (mean 1.05 %); IEF
70.72-88.95 % (0 % CS highest, 88.95 %; mean 80.93 %); removing urea cost 12.22-31.92 % (0 % CS
13.96, 40 % CS 12.22, 50 % CS 31.92); removing DTT cost 60.08 % on average (0 % CS −72.24 %, 40 % CS
−72.88 %, others ~−50 %); removing both: 0 % CS −81.29 %; removing SDS + thiourea + CHAPS −33.52 %.
Authors' ranking of what holds the extrudate together: disulfide > hydrophobic > hydrogen bonding.

## 4. Site densities the repository can use

No row is an isolate. Basis assumption: "per g protein" as printed; if Eq. 1 was applied to the
powder without dividing by the protein fraction, multiply by 1/0.7148 = 1.40 (the blends' stated
protein content). Native = raw powder blend, no heat in the assay beyond 50 °C (TNBS).

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN (arithmetic; basis) | conditions | source | evidence |
|---|---|---:|---|---|---|---|---|
| all-soy blend, 30 % CS (SPI + SPC + soy flour; ~71.5 % protein) | free SH | 2.66 ± 0.12 | µmol/g protein | 0.00266 as printed (0.0037 if per g powder) | native (R), 8 M urea overnight, DTNB | Table 4 | measured, blend, basis as printed |
| all-soy blend, 60 % CS | free SH | 2.83 ± 0.03 | µmol/g protein | 0.00283 (0.0040) | native (R) | Table 4 | measured, blend |
| PPI-containing blends, 50 % / 90 % CS | free SH | 3.12 ± 0.24 / 3.24 ± 0.04 | µmol/g protein | 0.00312 / 0.00324 (0.0044 / 0.0045) | native (R) | Table 4 | measured, blend |
| gluten-containing blends, 0 % / 40 % CS | free SH | 2.21 ± 0.10 / 2.61 ± 0.02 | µmol/g protein | 0.00221 / 0.00261 | native (R) | Table 4 | measured, blend |
| same six, after extrusion (T) | free SH | 2.03 / 1.99 / 2.06 / 2.03 / 1.96 / 1.32 (0 / 30 / 40 / 50 / 60 / 90 % CS) | µmol/g protein | 0.00203 / 0.00199 / 0.00206 / 0.00203 / 0.00196 / 0.00132 | low-moisture TVP, conditions in Flory 2023 | Table 4 | measured, blend |
| **free-SH extrusion factor (T/R), within blend** | ratio | 0.92 / 0.75 / 0.79 / 0.65 / 0.69 / 0.41 (mean of blends 1.90/2.78 = **0.68**) | — | basis-independent | as above | derived from Table 4 | within-study ratio |
| any blend | S-S, total SH | not reported (method described, no data) | — | — | — | — | absent |
| any blend | free amino (TNBS) | 3.47-6.90 (R); 1.67-3.90 (T) | mmol/g protein | **do not use as a density** (flag 1): 10-18x the lysine content of the constituent proteins | | Table 4 | measured, scale implausible |
| **TNBS free-amino extrusion factor (T/R), within blend** | ratio | 0.50 / 0.66 / 0.34 / 0.71 / 0.60 / 0.56 (0 / 30 / 40 / 50 / 60 / 90 % CS) | — | only meaningful if the scale error is multiplicative (flag 2) | as above | derived from Table 4 | within-study ratio, conditional |
| wheat gluten (ingredient) | free amino (TNBS) | 3.07 | mmol/g protein | do not use; 16x the lysine of whole-wheat protein and more for gluten | native | §3.3 text (Table S2) | measured, scale implausible |
| pea protein isolate (ingredient) | free amino | "highest among the raw ingredients" — value not printed | — | — | | §3.3 text (Table S2) | supplement, not on disk |
| pea protein isolate (ingredient) | H0, SDS-binding | 50.69 | µg SDS/mg protein | n/a (not a site density; = 0.176 µmol SDS/mg protein) | native | §3.5 text (Table S2) | measured, scale method-specific |
| blends | H0 extrusion factor (T/R) | 0.60 / 0.60 / 0.58 / 0.77 / 0.61 / 0.86 | — | — | | derived from Table 4 | within-study ratio |
| blends | protein content | ~71.48 | % | f_p = 0.7148 for the per-powder reading | | §2.1 | stated, method not given |

Comparison with the repo table (native, mmol per g protein): repo `soy_isolate.free_thiol` =
0.0078; the all-soy blends here read 0.0027-0.0028 as printed (0.34-0.36x) or 0.0037-0.0040 per
powder (0.47-0.51x); the same lab's SPI in `xiao2024_extraction.md` reads 0.0047. The blend value
is diluted by SPC and soy flour (non-protein solids, heat-treated concentrate), so the isolate row
should stay Ruan / Shimada; this paper adds only the extrusion factor.

## 5. Flags

1. **The free-amino column is not a lysine density.** Constituent proteins carry 0.19 (wheat),
   0.38 (soy) and 0.52 (pea) mmol lysine per g protein (USDA values via the book's Chapter 1,
   Table 2; Lys 146.19 g/mol); the blends should sit at ≤ 0.4-0.6 mmol/g protein plus a few
   hundredths of N-terminal α-NH2. Table 4 prints 3.47-6.90 for raw blends, 10-18x too high, and
   4-8x above the top of the paper's own leucine calibration (§2). Same finding as Shen 2022 (8.44)
   and the book chapter (9.48 SPI, 9.90 PPI, 5.01 gluten): a systematic scale of this lab's TNBS
   protocol. The soy amine pool must not be taken from here either. The all-soy blends (4.68, 5.29)
   are half the same lab's SPI (9.48) on the same protocol, so the scale is not even internally
   stable between papers.
2. **Extrusion factors for free NH2 are conditional.** If the artefact is a constant multiplier
   (e.g. a volume-basis slip), the T/R ratios 0.34-0.71 are real losses of TNBS-reactive amines
   during low-moisture extrusion (Maillard with the soy flour's sugars, isopeptide cross-links).
   If the artefact is an offset or a calibration non-linearity above the top standard (where all
   these readings sit), the ratios are meaningless. Osen 2015 (cited by the authors) found no free
   amino change on extruding PPI. Record as directional only.
3. **No isolate row.** Every measured sample is a blend; the SPI / SPC / soy flour / PPI / gluten
   values are in Tables S2-S3 of the supplement, which is not on disk (only gluten's 3.07 mmol/g
   protein and PPI's H0 50.69 are quoted in text). Retrieving the supplement would give the
   per-isolate free SH, S-S (Table S3, described as showing gluten "rich in disulfide bonds") and
   TNBS numbers under the book-chapter protocol.
4. **Blend recipes are qualitative** (which blends contain gluten / PPI / only soy) and the
   extrusion conditions are absent; both are in Flory 2023, not on disk. The "extrusion factor" is
   therefore for an unspecified low-moisture TVP process.
5. **Protein basis ambiguous**: Eq. 1 gives per g powder (C = 7.5 mg powder/mL), Table 4 says per
   g protein; the blends' 71.48 % protein (method unstated) is the only number available to
   convert, and whether it was applied is unknown (1.40x difference).
6. **Total SH and S-S not reported for any blend** although the method is described; the
   disulfide-driven insolubility of the extrudates is argued from the DTT solubility test (−60 %
   without DTT), not from an S-S assay.
7. **Free-SH heating direction is not universal**: the authors note that SPI, PPI and oat
   concentrate extruded alone showed INCREASED free SH in other studies (Li 2023; Pöri 2022;
   Zhang, Zhao 2022) while blends with gluten or starch decreased (Gao 2023; Zhang & Ryu 2023);
   this paper's 0.41-0.92x is one process on six blends.
8. **H0 is an SDS-binding number** (µg SDS / mg protein), not an ANS slope; it cannot be compared
   with Shen 2022 (ANS, 10^5 scale) or Chen 2022 (ANS, 10^2 scale) except by direction.
9. n = 2 throughout; several SDs of 0.00 in Table 2 (identical duplicates) suggest a single
   spectrum fitted twice for some rows.
10. Table 5 heads say "g/g material" while Eqs. 4-5 say "g/g protein"; the 0.25 g weighed is
    powder, so "material" is right.
