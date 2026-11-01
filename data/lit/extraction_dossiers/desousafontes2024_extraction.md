# de Sousa Fontes, Galvao, de Carvalho, Guedes, Lima, Bezerra & Madruga 2024 — EXTRACTION
DOI: 10.1016/j.foodchem.2024.138398 (printed p. 1; Crossref 2026-11-01: "Thiamine, cysteine and xylose added to the Maillard reaction of goat protein hydrolysate potentiates the formation of meat flavoring compounds", Food Chemistry 445, article 138398, print July 2024; seven authors, de Sousa Fontes first, Madruga last. **Matches the brief.**)
### Goat-viscera (liver : heart : lung 1:1:1) Alcalase hydrolysate, 20 mL + xylose 1.0 g + cysteine 0.1 g + thiamine 0.5 g, made up to 50 mL, initial pH 4.0 (F4) or 6.0 (F6) set with 1 M HCl, **unbuffered**, autoclave 121 °C / 60 min in 100 mL screw-cap flasks. HS-SPME GC-MS, two columns, **semi-quantified against one internal standard (1,2-dichlorobenzene), response factor 1, µg/L**. **Two pH points, one temperature, one time; final pH not reported; no thiamine-free, cysteine-free, xylose-free or lipid-free arm with numbers** (the supplement-only model systems M4 / M6 are "data not shown").

**Read 2026-11-01 from `data/articles/desousa-fontes2024.pdf`** (12 pp., born-digital). Prose from
`pdftotext -layout`; **Table 1 (p. 3) and Table 2 (pp. 4-6) checked cell by cell against 170 dpi renders of pp. 3-5**
(scratchpad `desousa-fontes2024_work/p3-03.png`, `p4-04.png`, `p5-05.png`); the renders also show which Table 2 values
are italic (= quantified on the polar HP-INNOWAX column, footnote p. 6), which the text layer loses. Figs. 1 (PCA biplot)
and 2 (heat map) carry no numbers and were not digitised. **Supplementary Tables S1 (fatty acids) and S2 (volatiles per
column) are not on disk.**

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Universidade Federal da Paraíba (UFPB), Department of Food Engineering, Technology Center, João Pessoa
(Madruga corresponding, msmadruga@uol.com.br; Bezerra, Galvão, de Carvalho, Guedes), with dos Santos Lima (IF Sertão
Pernambucano, Petrolina; amino-acid HPLC). **Not a group in `data/keys/laboratories.yml`: NEW**, proposed id
`ufpb_madruga`, ledger string "UFPB Joao Pessoa (Madruga)" — the string the sulfur ledger already uses for
`madruga1997` (Madruga sole author, DTQA/UFPB). [!] Do **not** map it to `reading_mottram`: `madruga1995` is Madruga's
Reading work with Mottram (Reading affiliation, UFPB only as "present address"); this 2024 paper is UFPB's own
method set (DVB/CAR/PDMS SPME, 1,2-dichlorobenzene IS, Agilent 7890B), sharing only the dichlorobenzene-IS habit.
Overlap in the corpus: `madruga1995_extraction.md` (beef pH 5.6 → 4.0 series, 140 °C) and `madruga1997_extraction.md`
(beef + IMP / thiamine / cysteine) are the same senior author's earlier pH and precursor work; no other goat or
offal-hydrolysate paper is keyed (the group's Cordeiro 2020 / 2022, Queiroz 2017 and Cunha 2023 hydrolysate papers
are cited, not on disk). Closest systems keyed: `gao2020b` (fish hydrolysate + Cys + Xyl), `cerny2007` (Cys + Xyl +
thiamine pH ladder, 145 °C).

## 1. What this paper answers for the model

| model need | what de Sousa Fontes 2024 gives | tag |
|---|---|---|
| pH direction of MFT / FFT (CF-012) | MFT 0.93, FFT 0.77, MFT-MFT 1.45 µg/L at initial pH 4; all **nd at pH 6** (no LOD printed) | [M] |
| pH direction of thiazoles / pyrazines / Strecker aldehydes | all higher at pH 6: thiazoles x8.5, pyrazines x100, Strecker aldehydes x26 [Z] | [M] / [Z] |
| thiamine with vs without (CF-089, CF-090, CF-007) | **no arm**: thiamine, cysteine and xylose are in every pot with numbers | — |
| lipid present (CF-031) | yes: ~0.5 g lipid / 100 g flavouring [Z] (organ lipid); **no lipid-free arm with numbers** | [M] / [Z] |
| absolute concentrations | **no**: one IS (1,2-dichlorobenzene), RF 1, no calibration, no LOD for volatiles | [M] |
| pH after heating | **not reported** (pH measurement is in the methods; no pH row in Table 1) | [!] |
| T, t | 121 °C, 60 min, autoclave (come-up time not stated) | [M] |
| replicates | n = 3, "triplicate sampling"; whether three pots or three analyses of one pot is not stated | [M] / [!] |
| precursor consumption | Table 1 free amino acids and sugars before (H) and after (F4, F6), external standards; basis ambiguous (sec. 6) | [M] / [!] |

**Answers.**
1. **pH series of MFT / FFT / thiazoles / pyrazines?** Two points only (initial pH 4.0 and 6.0), unbuffered, initial pH
   only (confidence ~97 %). MFT, FFT, bis(2-methyl-3-furyl) disulfide and 2-furfuryl 2-methyl-3-furyl disulfide are
   reported only at pH 4; thiazoles (dominated by 4-methyl-5-(2-hydroxyethyl)thiazole, HET), pyrazines, Strecker
   aldehydes, furfural, thiophenes and lipid volatiles are all higher at pH 6 (Table 2) [M].
2. **CF-012 direction:** same sign as the buffered laboratories (MFT falls as pH rises), here at 121 °C with thiamine
   present and in a protein-hydrolysate matrix. Weight is low: unbuffered pots, no final pH, nd without an LOD, and the
   FFT pool is **not** zero at pH 6 — 2-furfuryl methyl disulfide 2.05 µg/L appears only there (≈ 12.8 nM FFT-equivalents
   [Z], more than F4's free FFT + mixed disulfide, ≈ 6.9 nM [Z]). So at least the FFT "pH 4 only" claim is a
   speciation shift (free thiol → methyl disulfide), not absence. Probability this datum changes CF-012's status:
   ~10 %.
3. **Thiamine arms (CF-089, CF-090, CF-007):** none (confidence ~99 %). The model systems M4 / M6 (xylose + cysteine +
   thiamine in water, pH set with NaOH) contain all three; their 25 volatiles are "data not shown" [V], only in the
   Fig. 2 heat map. Nothing bounds a cysteine or xylose effect on thiamine routes.
4. **Lipid (CF-031):** the hydrolysate carries lipid (H 15.25 g/100 g dry, 89.85 % moisture → ≈ 1.5 g/100 g wet; the
   flavourings 9.07 / 8.71 g/100 g dry → ≈ 0.53 g/100 g wet [Z]); organ lipid (liver, heart, lung), so
   phospholipid-rich [Q from the matrix, not measured here; Table S1 missing]. No lipid ± arm, so no CF-031 contrast.
5. **Absolute or relative?** Semi-quantitative (confidence ~99 %): "semi-quantification" [V] by the IS-area ratio, IS
   50 ng in 2 mL; no response factors, no calibration, no LOD for volatiles; compounds read on the polar column (italics)
   are a different detector-response regime from the non-polar ones. Only within-compound F4 / F6 ratios are usable.
6. **pH measured after heating?** No (confidence ~90 %; methods say pH was measured on all samples but no value is
   printed anywhere, including Table 1).
7. **T and t:** 121 °C / 60 min in an autoclave, flasks cooled in ice [M]; come-up and cool-down not stated.

## 2. Conditions

- **Hydrolysate (H)** [M]: Saanen goat viscera (liver : heart : lung 1:1:1 by mass), minced, frozen -15 °C ≤ 120 days;
  endogenous enzymes killed 90 °C / 15 min; water 1:2 (w:v); Alcalase 0.8 % (E:S), 50 °C, pH-stat 8.0 with 0.5 M NaOH,
  240 min, 3600 rpm; enzyme killed 90 °C / 15 min. Degree of hydrolysis 26.49 % [M]. pH of H not printed (pH-stat at
  8.0, so presumably ~8 [Z]).
- **Flavourings F4 / F6** [M]: 20 mL H + 0.1 g cysteine + 1.0 g xylose + 0.5 g thiamine in a 100 mL screw-cap flask;
  pH adjusted to 4.0 or 6.0 "using 1 M HCl" [V]; made up to 50 mL with water; 121 °C / 60 min autoclave; ice bath;
  stored frozen ≤ 8 days for the non-volatile analyses. Formulation after Guo et al. 2010 [C].
  - Charged concentrations [Z]: cysteine 2.0 g/L = 16.5 mM (+ ~1.1 mM free Cys from H, if Table 1 is wet basis);
    xylose 20 g/L = 133 mM; thiamine 10 g/L = 29.6 mM as thiamine HCl (37.7 mM as the free cation) — **salt form not
    stated** [!]. Molar xylose : cysteine : thiamine ≈ 8 : 1 : 1.8 [Z] — thiamine exceeds cysteine, unlike Cerny's
    0.3 : 0.1 : 0.1 M ternary (CF-007).
  - H at 40 % v/v of the pot [Z]; flavouring dry matter 5.87 % (F4), 6.13 % (F6) [M, from moisture].
  - Buffer: **none** (hydrolysate peptides only). Headspace: ~50 mL air over 50 mL liquid in a 100 mL flask [Z];
    O2 not controlled.
- **Model systems M4 / M6** [M]: same supplements in 20 mL water, pH 4.0 / 6.0 "using 1 M NaOH" [V], made up to 50 mL,
  121 °C / 60 min. [!] F4 / F6 adjusted down with HCl, M4 / M6 up with NaOH: the ionic backgrounds differ. Volatiles
  not tabulated.
- **Replicates** [M]: "triplicate sampling" [V]; Table 1 and Table 2 give mean ± SD, n = 3; Tukey (Table 2), t-test
  (Table 1).

## 3. Quantitative results

### 3.1 Table 2 (pp. 4-6): Maillard / sulfur rows, µg/L (IS-equivalent, RF 1), mean ± SD, n = 3 [M]

Italic in the PDF (polar column) marked (p). nd = not detected (no LOD). GV and H columns are nd for every row below
unless given.

| compound | LRI | H | F4 (pH 4) | F6 (pH 6) | F6 / F4 [Z] |
|---|---|---|---|---|---|
| 2-methyl-3-furanthiol (MFT) # | 867 | nd | 0.93 ± 0.15 | nd | — |
| 2-furfurylthiol (FFT) # | 910 | nd | 0.77 ± 0.08 | nd | — |
| bis(2-methyl-3-furyl) disulfide | 1537 | nd | 1.45 ± 0.05 | nd | — |
| 2-furfuryl 2-methyl-3-furyl disulfide | 1641 | nd | 0.04 ± 0.00 | nd | — |
| 2-furfuryl methyl disulfide | 1213 | nd | nd | 2.05 ± 0.29 | — |
| methyl disulfide (p) | 1080 | nd | nd | 0.77 ± 0.10 | — |
| methional | 907 | nd | 0.77 ± 0.06 | 2.01 ± 0.01 | 2.6 |
| thiols and sulfides, total | | | 3.96 ± 0.34 | 4.83 ± 0.40 | 1.2 |
| thiazole (p) | 1262 | nd | 0.04 ± 0.00 | 4.31 ± 0.13 | 108 |
| 4-methylthiazole | 819 | nd | nd | 0.41 ± 0.02 | — |
| 2-acetylthiazole | 1017 | nd | 1.06 ± 0.00 | 3.46 ± 0.04 | 3.3 |
| 4-methyl-5-ethylthiazole (p) | 1026 | nd | 0.53 ± 0.01 | 23.00 ± 3.45 | 43 |
| 4-methyl-5-hydroxyethylthiazole (HET) | 1285 | nd | 47.14 ± 2.81 | 337.05 ± 6.73 | 7.2 |
| 5-ethenyl-4-methylthiazole (p) | 1542 | nd | nd | 47.48 ± 12.73 | — |
| thiazoles, total | | | 48.77 ± 2.82 | 415.71 ± 23.10 | 8.5 |
| 3-methylthiophene # | <800 | nd | 0.88 ± 0.08 | 1.25 ± 0.05 | 1.4 |
| 2-acetyl-5-methylthiophene | 1155 | nd | 1.01 ± 0.12 | nd | — |
| 2-formyl-3-methylthiophene | 1122 | nd | 2.97 ± 0.33 | 6.15 ± 0.03 | 2.1 |
| 2-formyl-2,3-dihydrothiophene | 1079 | nd | 1.37 ± 0.32 | 4.96 ± 0.35 | 3.6 |
| 3-thiophenecarboxaldehyde | 1000 | nd | 5.68 ± 0.00 | 9.26 ± 1.26 | 1.6 |
| thiophenes, total | | | 11.91 ± 0.85 | 21.62 ± 1.69 | 1.8 |
| 2-methylpropanal (p) | 820 | nd | 0.29 ± 0.14 | 13.38 ± 2.62 | 46 |
| 2-methylbutanal (p) | 923 | nd | 0.45 ± 0.10 | 23.29 ± 1.92 | 52 |
| 3-methylbutanal (p) | 925 | 0.03 ± 0.00 | 9.25 ± 1.84 | 372.38 ± 34.08 | 40 |
| 2-methylpentanal | <800 | nd | 6.08 ± 1.57 | 2.98 ± 0.19 | 0.49 |
| 2-methylpyrazine (p) | 1277 | 0.03 ± 0.01 | nd | 12.65 ± 0.41 | — |
| 2-propylpyrazine (p) | 1431 | nd | nd | 19.97 ± 1.81 | — |
| 2,5-dimethylpyrazine (p) | 913 | 2.34 ± 0.13 | 0.64 ± 0.26 | 32.22 ± 3.78 | 50 |
| 2-methyl-5-isopropylpyrazine | 1055 | 1.23 ± 0.01 | nd | nd | — |
| 1-furfurylpyrrole | 1184 | nd | nd | 5.41 ± 0.06 | — |

\# = "thermal degradation product of thiamine" per the authors' footnote [Q]; FFT carries the mark too.

Molar [Z] (MW 114.17 thiols, 226.3 disulfides, 160.26 2-furfuryl methyl disulfide, 143.21 HET): F4 MFT 8.1 nM free
+ 12.8 nM in MFT-MFT + 0.18 nM in the mixed disulfide ≈ 21 nM MFT-equivalents; FFT 6.7 + 0.18 ≈ 6.9 nM;
MFT / FFT 1.2 free, ≈ 3.1 in equivalents. F6 FFT-equivalents in 2-furfuryl methyl disulfide ≈ 12.8 nM. HET 0.33 µM
(F4), 2.35 µM (F6). MFT-equivalents per charged thiamine ≈ 7 × 10⁻⁷ (F4).

### 3.2 Table 2: sugar-derived furans, µg/L [M]

| compound | LRI | F4 | F6 | F6 / F4 [Z] |
|---|---|---|---|---|
| furfural | 845 | 1891.84 ± 166.08 | 3188.44 ± 36.86 | 1.69 |
| 2-(hydroxymethyl)furan | 859 | 2.50 ± 0.12 | 11.55 ± 0.94 | 4.6 |
| 2-methylfuran # (p) | 876 | 0.52 ± 0.34 | 15.14 ± 0.35 | 29 |
| 2-acetylfuran (p) | 1522 | nd | 40.50 ± 0.69 | — |
| 2-furyl methyl ketone | 913 | nd | 2.44 ± 0.00 | — |
| 2,2'-bifuran | 1039 | 2.87 ± 0.00 | 6.46 ± 0.40 | 2.3 |
| furan (p) | 808 | 0.13 ± 0.02 | 4.20 ± 0.97 | 32 |

Furfural ≈ 19.7 µM (F4) and 33.2 µM (F6), ≈ 1.5-2.5 × 10⁻⁴ of the charged xylose [Z]. Absent from H.

### 3.3 Table 2: lipid-type volatiles (selected), µg/L [M]

| compound | LRI | GV (µg/kg) | H | F4 | F6 | F6 / F4 [Z] |
|---|---|---|---|---|---|---|
| pentanal (p) | 979 | nd | 0.03 ± 0.00 | 0.39 ± 0.11 | 11.44 ± 0.21 | 29 |
| hexanal | 800 | 0.69 ± 0.03 | 0.48 ± 0.03 | 5.85 ± 0.62 | 8.40 ± 0.81 | 1.4 |
| heptanal | 902 | 0.21 ± 0.03 | 0.46 ± 0.03 | 1.65 ± 0.20 | 3.21 ± 0.26 | 1.9 |
| octanal (p) | 1294 | nd | 0.04 ± 0.01 | 0.50 ± 0.03 | 20.50 ± 0.76 | 41 |
| nonanal | 1104 | 1.83 ± 0.37 | 3.40 ± 0.18 | 20.88 ± 2.42 | 37.25 ± 2.73 | 1.8 |
| decanal | 1206 | nd | 0.14 ± 0.03 | 2.53 ± 0.26 | 2.92 ± 0.29 | 1.2 |
| hexadecanal | 1817 | 1.16 ± 0.09 | nd | 33.01 ± 3.68 | 58.67 ± 4.45 | 1.8 |
| (E)-2-octenal | 1060 | nd | nd | 1.82 ± 0.16 | 5.68 ± 0.28 | 3.1 |
| (Z)-2-nonenal | 1161 | nd | nd | 2.54 ± 0.17 | 4.46 ± 0.30 | 1.8 |
| (E)-2-decenal | 1267 | nd | nd | 0.86 ± 0.09 | 1.57 ± 0.07 | 1.8 |
| (E,E)-2,4-decadienal | 1319 | nd | nd | 0.50 ± 0.03 | nd | — |
| benzaldehyde | 960 | nd | 29.35 ± 9.70 | 69.80 ± 2.85 | 68.88 ± 3.49 | 0.99 |
| 1-octen-3-ol | 981 | 1.27 ± 0.10 | 0.53 ± 0.29 | 2.38 ± 0.27 | 8.96 ± 0.86 | 3.8 |
| 2-pentylfuran | 992 | nd | 3.43 ± 0.12 | 6.57 ± 0.97 | 62.72 ± 7.38 | 9.5 |
| (E)-2-(2-pentenyl)furan (p) | 1303 | nd | nd | 1.18 ± 0.62 | 23.82 ± 1.69 | 20 |
| 2-octanone (p) | 1289 | nd | 0.23 ± 0.07 | 0.41 ± 0.15 | 19.90 ± 3.81 | 49 |
| 2-nonanone | 1093 | nd | 1.38 ± 0.23 | 0.83 ± 0.20 | 3.91 ± 0.34 | 4.7 |
| 2-decanone | 1193 | nd | 3.61 ± 0.02 | 2.66 ± 0.11 | 7.80 ± 0.18 | 2.9 |
| caproic (hexanoic) acid (p) | 1850 | nd | nd | 12.41 ± 1.02 | 152.78 ± 14.59 | 12 |
| aldehydes, total (25) | | 6.42 ± 0.61 | 34.76 ± 10.17 | 177.34 ± 16.33 | 327.42 ± 22.52 | 1.8 |
| ketones, total (19) | | 5.58 ± 0.52 | 7.91 ± 0.86 | 9.89 ± 1.03 | 53.60 ± 6.74 | 5.4 |

Grand totals (text p. 7): F4 2252.53, F6 4943.98 µg/L [M]; the class totals of Table 2 sum to 2252.53 (F4, exact) and
4941.89 (F6) [Z].

### 3.4 Table 1 (p. 3): composition and precursors, mean ± SD, n = 3 [M]

| parameter | GV | H | F4 | F6 |
|---|---|---|---|---|
| protein, g/100 g dry | 79.93 ± 0.40 | 72.01 ± 0.69 | 65.14 ± 0.34 | 61.92 ± 0.16 |
| moisture, g/100 g wet | 75.08 ± 0.51 | 89.85 ± 0.03 | 94.13 ± 0.02 | 93.87 ± 0.01 |
| lipids, g/100 g dry | 22.09 ± 0.39 | 15.25 ± 0.72 | 9.07 ± 0.65 | 8.71 ± 0.25 |
| carbohydrates (by difference), g/100 g dry | 2.40 ± 1.73 | 2.32 ± 0.37 | 18.03 ± 1.63 | 20.08 ± 3.68 |
| ash, g/100 g dry | 4.43 ± 0.13 | 10.55 ± 0.23 | 7.16 ± 0.32 | 7.45 ± 0.41 |
| glucose, mg/100 g | 1691.95 ± 37.31 | 465.74 ± 5.09 | 165.05 ± 4.21 | 180.82 ± 2.82 |
| fructose | 240.02 ± 0.42 | 115.28 ± 4.85 | 80.27 ± 9.21 | 131.33 ± 7.77 |
| ribose | 19.87 ± 0.83 | 9.21 ± 1.77 | 11.14 ± 0.64 | 19.09 ± 0.93 |
| xylose (added) | na | na | 2258.23 ± 70.68 | 1998.66 ± 13.68 |
| cysteine (free), mg/100 g | na | 33.57 ± 0.76 | 16.76 ± 12.03 | 5.61 ± 0.00 |
| lysine | na | 145.21 ± 4.66 | 26.58 ± 0.03 | 17.59 ± 10.23 |
| glycine | na | 46.64 ± 3.86 | 13.01 ± 0.09 | 9.44 ± 3.08 |
| alanine | na | 84.42 ± 4.18 | 19.47 ± 0.35 | 13.28 ± 7.15 |
| leucine | na | 20.72 ± 0.00 | < LOD | < LOD |
| valine + methionine | na | 40.36 ± 0.45 | 5.82 ± 0.10 | 4.73 ± 1.53 |

(Other amino acids, mannose and maltose are in Table 1; not transcribed.) Sugars by HPLC-RID with external standards
(R² > 0.9992); amino acids by OPA-HPLC-DAD with external standards (R² > 0.995). **The basis (wet / dry) of the sugar
and amino-acid rows is not stated** (sec. 6).

## 4. Analytical method

- **Volatiles** [M]: 2 mL liquid (or 2 g GV) + 1 µL 1,2-dichlorobenzene (50 µg/mL in MeOH, i.e. 50 ng) in a 20 mL vial;
  equilibrate 50 °C / 10 min, DVB/CAR/PDMS 50/30 µm fibre 30 min with agitation; Agilent 7890B / 5977B, VF-5MS
  (30 m) and HP-INNOWAX (60 m), splitless 250 °C, 40 °C (3 min) → 5 °C/min → 120 °C → 10 °C/min → 230 °C (5 min);
  EI 70 eV, m/z 50-400. Identification by NIST 2014 + LRI (C6-C20 alkanes); **no authentic standards named**.
  Quantity = area ratio to the IS × IS amount (RF 1) [M]; LOD / LOQ for volatiles not given.
- **Fatty acids** [M]: direct transesterification (Hartmann & Lago), GC-FID, C21 IS; results in Table S1 (missing).
- **Proximate** [M]: AOAC 926.07B, 923.03, 930.25 (N × 6.25); lipid by Folch. pH meter calibrated at pH 4 and 7
  (values never reported).

## 5. Bearing on open questions

- **CF-012 (pH direction of MFT)** — weakly SUPPORTS "MFT falls as pH rises" at 121 °C with thiamine present, against
  the reading that a near-neutral optimum takes over at 90-120 °C. Weight low (unbuffered, no final pH, nd with no
  LOD, semi-quantitative). P(changes CF-012) ≈ 10 %.
- **CF-089 / CF-090 / CF-007** — LEAVES UNCHANGED: no thiamine, cysteine or xylose omission arm with numbers
  (P(changes) < 2 %).
- **CF-031 (lipid effect on MFT / FFT)** — LEAVES UNCHANGED: lipid is present (~0.5 % w/w) but never varied
  (P < 2 %).
- **CF-033 (lipid volatiles in lipid-free pots)** — LEAVES UNCHANGED: the lipid-free M4 / M6 volatiles are not shown.
- **CF-088 (fall of lipid aldehydes with Cys + pentose)** — LEAVES UNCHANGED: no ± Cys / pentose contrast; aldehydes
  rise H → F because the pot is heated, not because reactants are removed.
- **CF-013 / CF-015** — LEAVES UNCHANGED (one time point; MFT / FFT only at one pH).
- **Thiamine route pH (engine: `k_thi_hmp` first order, no pH term; `parameters_sulfur.thiamine_h2s_share` zero below
  pH 6.25)** — HET, the thiazole half of thiamine, is x7.2 higher at initial pH 6 than pH 4, consistent with faster
  thiamine cleavage at higher pH; but `cerny2007` Table 2 (buffered, 145 °C) has HET **falling** pH 4 → 6 (525 → 161
  peak-area units). The authors misreport Cerny & Briffod as finding HET not pH-dependent. A within-pair
  disagreement, not a datum to fit (SPME of polar HET with RF 1).
- **Furfural pH direction** — furfural x1.7 higher at initial pH 6; `cerny2007` Table 2 has furfural 208 at pH 4 and
  0 at pH 6. Unbuffered pots drift acid on heating, so F6's in-situ pH may be well below 6 (not measured).

## 6. Flags

1. [!] **Final pH not reported** although measured (methods 2.3.1); both pots unbuffered, so the label pH is the
   start only (ph_start_only, unbuffered).
2. [!] **No LOD for volatiles**: the pH-6 "nd" for MFT / FFT / MFT-MFT cannot be bounded.
3. [!] **FFT not absent at pH 6**: 2-furfuryl methyl disulfide (F6 only, 2.05 µg/L) carries ≈ 12.8 nM FFT-equivalents
   [Z]; the "FFT only at pH 4" conclusion confuses free-thiol speciation with formation.
4. [!] **Sugar / amino-acid basis ambiguous** [Z]: read wet, xylose 2258 mg/100 g (F4) exceeds the charged 2000 mg/100 mL
   (113 % recovery) and exceeds the total carbohydrate by difference (18.03 % of 5.87 % dry matter ≈ 1.06 g/100 g wet);
   read dry, xylose would be ~93 % consumed in 1 h at 121 °C, implausible beside furfural at 2 × 10⁻⁴ of the xylose.
   GV glucose (1.69 g/100 g) fits only a dry basis (carbohydrate 2.40 % dry). Do not use Table 1 for xylose or
   cysteine consumption. Free cysteine by OPA (poor for cysteine) with SD 12.03 on a mean of 16.76 adds doubt.
5. [!] Amino-acid LOD printed as "LD < 0.289 mg/L" (methods) but "< LOD (leucine) = 1,358 mg/L" (Table 1 footnote).
6. [!] **Thiamine salt form not stated**; 29.6 vs 37.7 mM [Z].
7. [!] **pH adjusted with HCl in F4 / F6 but NaOH in M4 / M6**; H's own pH not printed.
8. [!] **Replicate unit unclear** ("triplicate sampling"); several SDs print as 0.00 (2-acetylthiazole F4,
   3-thiophenecarboxaldehyde F4, 2,2'-bifuran F4), suggesting analytical rather than pot replicates.
9. [!] Significance letters inverted for nonanal and (Z)-2-nonenal (higher F6 value lettered B; elsewhere A = higher).
10. [!] Identification by library + LRI only; the single LRI column mixes polar and non-polar values without saying
    which (octanal 1294 and pentanal 979 are polar values and are italic, but 2,5-dimethylpyrazine, also italic, is
    listed at 913, a non-polar value); some identities are doubtful (3-methoxythiophenol 6.36 µg/L in F4 only).
11. [!] Text vs table: hexanal F6 8.5 (text) vs 8.40 (Table 2); thiophenes F4 11.92 (text) vs 11.91.
12. [!] Supplement (Tables S1, S2) not on disk.

## 7. Ledger-ready measurements

- Table 2 sulfur / thiazole / thiophene rows, F4 vs F6 (sulfur lane) — two-point, initial-pH, semi-quantitative.
- Table 2 Strecker aldehydes, pyrazines, pyrrole (aroma_other) and sugar-derived furans (trunk), F4 vs F6, plus H.
- Table 2 lipid-type rows, GV → H → F4 / F6 (lipid lane) — a lipid-containing Maillard pot at two pH.
- Table 1 precursor composition (aroma_other) — external standards, but basis ambiguous; composition context only.
- Not keyable: any rate constant, any omission contrast (thiamine, cysteine, xylose, lipid), any final pH, M4 / M6.
