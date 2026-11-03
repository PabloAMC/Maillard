# Jacob, Hutzler & Methner 2018/2019 — EXTRACTION (file stem `jacob2018`)
**Eur. Food Res. Technol. 245(1):95-109**, DOI **10.1007/s00217-018-3143-z** (printed p. 1; the PDF is the
online-first version without page numbers). Crossref 2026-11-03: "Comparison of various industrially applicable
disruption methods to produce yeast extract using spent yeast from top-fermenting beer production: influence on amino
acid and protein content", Friedrich Felix Jacob, Mathias Hutzler, Frank-Jürgen Methner, published online
**2018-08-23**, EFRT 245(1), 95-109 (print issue 2019). **Matches the PDF and the brief.**
[!] Naming: stem says 2018 (online year); `jacob2019_extraction.md` and `tome2021_extraction.md` cite it as
"Jacob, Hutzler & Methner 2019" (print year). Same paper. **Companion of** `jacob2019_extraction.md` (Jacob,
Striegel, Rychlik, Hutzler & Methner 2019, EFRT 245:1169, same yeast batch design, same three disruption methods,
sums and sugars / vitamins / GSH).

### Spent top-fermenting brewer's yeast (S. cerevisiae TUM 68, 50 L pilot CCT, washed ×3, 7 % dry matter, 200 mL) disrupted by Vibrogen cell mill (0.5 mm glass beads, 15 min, 7 °C), ultrasonic sonotrode (400 W, 30 min, 7 °C) or autolysis (24 h, 50 °C, NaCl 0.086 M + ethyl acetate 0.051 M); centrifuged, freeze-dried; plus one industrial autolysed yeast extract of unknown origin. **Individual free amino acids (20 proteinogenic incl. Cys, Gln, Asn; plus GABA, ornithine) and individual total amino acids after acid hydrolysis, mg/g dry mass of extract, mean ± 95 % CI of triplicates (Figs. 3-7, bar charts only)**; protein (Bradford), free- and total-amino-acid sums (Fig. 8); essential amino acids per g protein and chemical scores (Table 2). **No sugars, no nucleotides, no glutathione, no pH, no heating** (those are in `jacob2019`, except pH and heating, which neither has).

**Read 2026-11-03 from `data/articles/jacob2018.pdf`** (15 pp., Springer, born-digital, all bar charts vector).
Prose from `pdftotext`; **pp. 5-10 rendered at 110 dpi** (scratchpad `ingest/jacob2018_work/p-0N.png`) to identify
figures and axes; **Figs. 3-8 digitised from the PDF vector paths** (`pdftocairo -svg`; bar rectangles and error-bar
strokes converted with the gridline positions; scripts `scratchpad/scripts/svg_bars.py`, `bars_to_values.awk`;
outputs `fig3_bars.txt`, `p7-p10_bars.txt`). Read error ≈ ±0.05 mg/g on Figs. 3-5, ±0.1 on Figs. 6-7, ±1 on Fig. 8
(bar edges are exact vector coordinates; the error is the gridline calibration). **Table 2 (p. 10) read from the
render and checked cell by cell against the text layer.** Figs. 2, 9-11 seen, not digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** **TU Berlin, Fachgebiet Brauwesen** (Institut für Lebensmitteltechnologie und Lebensmittelchemie);
corresponding and first author **F. F. Jacob**, senior author **F.-J. Methner**; Hutzler at Forschungszentrum
Weihenstephan (TU München; strain, QC PCR). **Existing id `tu_berlin_methner`** ("TU Berlin Brewing (Jacob /
Methner)", `data/keys/laboratories.yml` l. 84, evidence already names the 2018 paper). Same lab:
`jacob2019_extraction.md`.

## 1. What this paper answers for the model

| model need (yeast extract as a forced precursor input) | what Jacob 2018 gives | tag |
|---|---|---|
| individual free amino acids | Figs. 3-4: 20 amino acids × 4 extracts, mg/g dry mass ± 95 % CI | [F] |
| free cysteine | Fig. 3: 4.0 / 6.1 / 10.0 / 13.0 mg/g (cell mill / sonotrode / autolysis / industrial); method VDLUFA 4.11.1, redox state not stated | [F] / [!] |
| free methionine | Fig. 4: 3.0 / 5.2 / 9.1 / 11.3 mg/g | [F] |
| free Asn, Gln (acrylamide, Maillard amines) | Asn 4.2 / 6.4 / 14.0 / 18.8; Gln 6.3 / 9.7 / 13.4 / 14.6 mg/g | [F] |
| ornithine, GABA | Fig. 5 | [F] |
| total (hydrolysed) amino acids | Figs. 6-7 individual; Fig. 8 sums (440-633 mg/g) | [F] |
| glutathione | **not measured here** (Ellman thiol is in `jacob2019` Fig. 4) | — |
| reducing sugars, ribose, nucleotides | **not measured here** | [!] |
| heating / kinetics | **none**; extracts analysed as made within 12 h | [!] |
| basis | mg per g dry mass of extract (powder had ~10 % residual water) | [M] |

**Answers.**
1. **Per-amino-acid amounts with units and basis: yes, figure-only** (confidence ~90 % in the digitised values,
   which are vector-exact to ±0.05-0.1 mg/g). mg/g dry mass of freeze-dried extract; ± = 95 % CI from Student's t on
   triplicates [M]. The industrial product carries a constant nominal bar (±0.10 on Figs. 3-5, ±0.2 on Figs. 6-7)
   — read as **no replicate error** for it (P ≈ 80 %).
2. **Cysteine: free and total given; reduced cysteine not established** (confidence ~70 % that the free "Cys" is
   mostly cystine in the autolysates). Free Cys and Pro were measured by a separate feed-analysis method (VDLUFA 4.11.1)
   [M]; that method number is, to my knowledge, the performic-acid-oxidation amino-acid method, which reads cysteine
   + cystine together as cysteic acid [C, not stated in the paper]. Cross-check with the companion [Z]: lab autolysate
   free Cys 10.0 mg/g = 82 µmol/g, but `jacob2019` Ellman total free thiol (as GSH) is only 1.2 mg/g = 3.9 µmol/g;
   reduced Cys at 82 µmol/g would have read ≥ 20× higher on DTNB. So in the 24 h / 50 °C autolysate the free
   cyst(e)ine is essentially all disulfide, or one of the two assays is wrong. For mechanical extracts the two agree in
   magnitude (free Cys 33 / 50 µmol/g vs DTNB 42 / 43 µmol/g GSH-eq), so their "GSH" may be partly free Cys.
3. **Methionine** (confidence ~90 %): free 2.97 ± 0.62 / 5.17 ± 1.82 / 9.13 ± 1.11 / 11.32 mg/g; total 18.3-28.3 mg/g.
   **Glutathione: not here.** **Sugars, nucleotides: not here** (see `jacob2019`: glucose 1.7-3.2, fructose
   0.5-1.0 mg/g, RNA 49-55 mg/g, no 5'-nucleotides).
4. **Does disruption method matter? Strongly for free amino acids, not for totals** (~90 %): autolysis frees ×2-5
   more of each amino acid than the cell mill (except Glu, Ser); totals do not differ (ANOVA p > 0.05) except Glu,
   lower after autolysis [M]. Heat-inactivated yeast (80 °C, 15 min) then milled gave < 1 mg/g total free amino acids
   [M]: the free pool is made by the yeast's own proteases after disruption, not pre-existing.

## 2. Conditions

- **Yeast** [M]: S. cerevisiae LeoBavaricus TUM 68 (top-fermenting), propagated 40 mL → 400 mL → 4 L, 20 °C; 50 L
  aerated wort (12.4 °P, pH 5.19, Döhler all-malt concentrate, Table 1) in a CCT, 5 d at 18 °C unpressurised; ~1 L
  core yeast cropped day 6 (70 % live); washed ×3 with distilled water at 10 % DM, then 7 % DM.
- **Disruption** [M]: cell mill 75 Hz, beads:suspension 1:2 w/w, 15 min, jacket-cooled to 7 °C; sonotrode HD 3400,
  25 mm tip, 20 kHz, 400 W continuous, 30 min, 7 °C; autolysis 24 h, 50 °C, 100 rpm, NaCl 0.086 mol/L + ethyl acetate
  0.051 mol/L. Nitrogen yield 98 % (autolysis) vs 80 % (mechanical) [M] (total N 554 → 543 / 443 mg/100 mL).
- **Work-up** [M]: 10 000 g, 20 min, 4 °C; supernatant freeze-dried (−50 °C, 0.1 mbar) to ~10 % water; analysed
  within 12 h, resuspended in distilled water. Industrial autolysed extract: origin yeast unknown [M].
- pH of extracts, aw, O2: **not reported**. Replicates: triplicate disruptions [M].

## 3. Analytical methods

- **Free amino acids** [M]: HPLC, MEBAK 2.6.4.1 (derivatisation not described), except **Pro and Cys: VDLUFA method
  4.11.1** (feedstuff analysis). **Total amino acids**: acid hydrolysis, same method 4.11.1. No LOD, no recovery,
  no standards listed [!]. Trp not quantified in totals (acid-labile; "Nq" in Table 2).
- **Protein**: Bradford, 595 nm (> 3 kDa only); Kjeldahl total N (MEBAK 2.6.1.1) for disruption yield. **FAN**:
  ninhydrin (MEBAK 2.6.4.1). Protein size: Agilent Bioanalyzer Protein 80.
- **Statistics** [M]: means of three; ± = 95 % confidence interval (Student's t); paired t test between methods;
  one-way ANOVA; significance p < 0.05.

## 4. Data (mg/g dry mass of yeast extract; CM cell mill, S sonotrode, A lab autolysis, IP industrial autolysed product; ± 95 % CI) [F, vector]

### 4a. Figs. 3-4, free proteinogenic amino acids

| aa | CM | S | A | IP |
|---|---|---|---|---|
| Gln | 6.25 ± 0.69 | 9.65 ± 1.65 | 13.43 ± 0.99 | 14.57 |
| His | 4.38 ± 0.63 | 4.43 ± 0.81 | 9.08 ± 0.37 | 8.39 |
| Val | 6.89 ± 1.07 | 8.74 ± 0.52 | 21.39 ± 1.11 | 27.19 |
| Ile | 4.15 ± 0.35 | 6.71 ± 0.88 | 18.43 ± 2.48 | 22.99 |
| Leu | 5.70 ± 0.57 | 10.32 ± 0.74 | 31.50 ± 6.45 | 41.26 |
| Pro | 5.44 ± 0.40 | 10.60 ± 2.73 | 16.47 ± 0.16 | 16.12 |
| **Cys** | **4.01 ± 0.50** | **6.09 ± 0.55** | **9.98 ± 1.25** | **13.04** |
| Asp | 3.70 ± 0.21 | 3.34 ± 1.90 | 13.42 ± 0.46 | 15.41 |
| Glu | 18.38 ± 2.91 | 19.85 ± 2.58 | 23.26 ± 3.02 | 40.00 |
| Asn | 4.22 ± 0.16 | 6.35 ± 0.66 | 14.00 ± 0.51 | 18.83 |
| Ser | 8.86 ± 0.16 | 7.28 ± 1.36 | 9.36 ± 0.81 | 16.96 |
| Gly | 3.77 ± 1.05 | 6.35 ± 0.81 | 7.88 ± 0.51 | 7.99 |
| Thr | 4.09 ± 0.36 | 4.99 ± 0.66 | 12.97 ± 1.08 | 17.39 |
| Ala | 9.54 ± 0.66 | 10.37 ± 0.39 | 21.59 ± 1.28 | 39.61 |
| Arg | 13.53 ± 0.99 | 11.32 ± 1.90 | 22.04 ± 4.22 | 19.96 |
| Tyr | 4.44 ± 0.21 | 7.25 ± 1.09 | 16.92 ± 2.69 | 9.22 |
| **Met** | **2.97 ± 0.62** | **5.17 ± 1.82** | **9.13 ± 1.11** | **11.32** |
| Trp | 2.83 ± 0.07 | 6.27 ± 0.47 | 7.85 ± 1.26 | 10.30 |
| Phe | 4.49 ± 0.42 | 8.03 ± 1.24 | 19.41 ± 2.10 | 26.25 |
| Lys | 5.30 ± 0.76 | 8.31 ± 0.50 | 18.27 ± 0.75 | 22.76 |
| sum of 20 [Z] | 122.9 | 161.4 | 316.4 | 399.6 |
| printed sum (Fig. 8; abstract) | 115.7 (115) | 155.3 (155) | 306.0 (307) | 433.1 |

Molar [Z] (MW 121.16 / 149.21 / 132.12): free Cys 33 / 50 / 82 / 108 µmol/g; Met 20 / 35 / 61 / 76; Asn 32 / 48 / 106 /
143. Text [M]: industrial product has "up to 20 mg more" Ala and Glu than the lab autolysate (Fig.: +18.0, +16.7 ✓).

### 4b. Fig. 5, free non-proteinogenic amino acids

| | CM | S | A | IP |
|---|---|---|---|---|
| GABA | 1.89 ± 0.21 | 3.28 ± 0.70 | 9.87 ± 0.54 | 3.27 |
| ornithine | 9.03 ± 2.04 | 8.13 ± 1.02 | 9.80 ± 1.23 | 16.64 |

ANOVA: GABA p < 0.01, ornithine p > 0.05 [M]. Text "at 10 mg/g" GABA after autolysis ✓.

### 4c. Figs. 6-7, total amino acids after acid hydrolysis

| aa | CM | S | A | IP |
|---|---|---|---|---|
| His | 11.73 ± 1.34 | 12.35 ± 1.73 | 11.37 ± 2.17 | 14.07 |
| Thr | 16.66 ± 1.16 | 18.90 ± 1.18 | 20.96 ± 2.24 | 27.01 |
| Val | 21.71 ± 1.34 | 23.02 ± 0.99 | 24.99 ± 3.74 | 35.88 |
| Met | 18.29 ± 1.39 | 20.30 ± 1.50 | 20.90 ± 1.48 | 28.30 |
| Ile | 17.16 ± 1.15 | 19.39 ± 0.92 | 19.10 ± 1.74 | 32.46 |
| Phe | 19.01 ± 0.43 | 20.84 ± 1.26 | 21.79 ± 2.10 | 31.27 |
| Leu | 29.45 ± 1.18 | 32.59 ± 2.07 | 34.09 ± 4.38 | 51.07 |
| Lys | 30.17 ± 0.94 | 29.94 ± 3.32 | 29.56 ± 1.30 | 40.80 |
| Asp | 35.39 ± 1.23 | 36.58 ± 2.82 | 36.83 ± 0.93 | 47.57 |
| Glu | 61.82 ± 0.96 | 65.89 ± 5.39 | 53.42 ± 4.75 | 114.61 |
| Asn | 16.63 ± 2.15 | 18.59 ± 2.57 | 17.99 ± 4.71 | 28.39 |
| Ser | 14.13 ± 1.30 | 16.55 ± 1.46 | 16.76 ± 0.64 | 24.14 |
| Gln | 15.40 ± 2.04 | 16.46 ± 1.61 | 17.48 ± 1.00 | 20.37 |
| Gly | 16.37 ± 0.81 | 17.44 ± 0.96 | 18.96 ± 1.19 | 22.91 |
| Ala | 22.78 ± 0.93 | 24.77 ± 1.78 | 25.46 ± 0.57 | 44.80 |
| Arg | 23.51 ± 1.13 | 24.52 ± 2.82 | 26.35 ± 2.25 | 34.11 |
| Tyr | 14.72 ± 0.95 | 16.51 ± 0.95 | 19.98 ± 4.71 | 12.81 |
| Pro | 15.70 ± 2.46 | 18.37 ± 1.46 | 17.78 ± 6.90 | 26.90 |
| **Cys** | **16.89 ± 1.35** | **16.00 ± 1.23** | **16.51 ± 1.00** | **18.08** |
| sum of 19 [Z] | 417.5 | 449.0 | 450.3 | 655.6 |
| printed total (Fig. 8) | 440.1 | 449.3 | 447.2 | 633.3 |

### 4d. Fig. 8, protein (Bradford), free amino acids and total amino acids (stacked bars)

| | CM | S | A | IP |
|---|---|---|---|---|
| protein (Bradford) | 321.6 ± 9.8 | 285.6 ± 8.9 | 53.3 ± 12.6 | 102.2 ± 2.1 |
| free amino acids | 115.7 ± 1.2 | 155.3 ± 1.2 | 306.0 ± 2.7 | 433.1 |
| total amino acids | 440.1 ± 5.0 | 449.3 ± 5.0 | 447.2 ± 5.0 | 633.3 ± 4.9 |

Text [M]: protein 30 / 27 / 5 % dm (Fig.: 32.2 / 28.6 / 5.3); total amino acids 43 % for all three (Fig.: 44.0 /
44.9 / 44.7); **FAN 16.5 / 25.8 / 45.1 / 58.1 mg/g** (CM / S / A / IP).

### 4e. Table 2 (p. 10), essential amino acids, mg/g protein [M]

| | CM | S | A | IP | egg | FAO/WHO |
|---|---|---|---|---|---|---|
| His | 26.49 | 27.50 | 25.46 | 22.23 | 25.00 | 15.00 |
| Thr | 37.57 | 42.08 | 46.91 | 42.66 | 47.00 | 11.00 |
| Val | 48.99 | 51.27 | 55.93 | 56.67 | 66.00 | 15.00 |
| Leu | 65.08 | 72.58 | 76.28 | 80.66 | 86.00 | 21.00 |
| Ile | 38.72 | 43.21 | 42.73 | 51.25 | 54.00 | 15.00 |
| Met + Cys | 41.30 | 45.20 | 46.75 | 44.69 | 57.00 | 20.00 |
| Lys | 68.11 | 66.68 | 66.13 | 64.44 | 70.00 | 18.00 |
| Phe + Tyr | 42.90 | 46.38 | 48.77 | 49.38 | 93.00 | 21.00 |
| Trp | Nq | Nq | Nq | Nq | 25.00 | 15.00 |
| Total | 369.17 | 394.89 | 408.95 | 411.98 | 498.00 | 136.00 |

Chemical scores (CS egg 46-110; CS FAO/WHO 148-426) printed; nutritional, not transcribed. "Protein" here = total
amino acids: Fig. 6 His / Table 2 His = 0.443 / 0.449 / 0.447 / 0.633 g/g, i.e. the Fig. 8 totals [Z].

### 4f. Not digitised
Fig. 2 (N-based disruption effectiveness vs time, rising to 80 % mechanical / 98 % autolysis), Fig. 9 (variance vs
mean of released amino acids, autolysis R² = 0.7558 [M]), Fig. 10 (free / total per amino acid, derivable from 4a /
4c), Fig. 11 (Bioanalyzer gel: autolysates < 4 kDa only, mechanical 3.5-63 kDa [Q]).

## 5. Flags and internal consistency

1. [!] **"Cys" redox state unknown**; free Cys in the lab autolysate (82 µmol/g) is ≥ 20× the companion's DTNB thiol
   (3.9 µmol/g GSH-eq) — the free "Cys" must be mostly cystine, or an assay is wrong (sec. 1, ans. 2). Do not enter it
   as reduced cysteine.
2. [!] **Table 2 "Met + Cys" and "Phe + Tyr" are Met and Phe alone** [Z]: CM Met 18.29 / 0.4401 = 41.6 vs 41.30 printed;
   Phe 19.01 / 0.4401 = 43.2 vs 42.90; same for S and A. Adding Cys (16-18 mg/g) and Tyr would raise them by ~38 and
   ~33 mg/g protein; the sulfur and aromatic chemical scores are understated.
3. [!] **Total Gln and Asn after acid hydrolysis** (Fig. 7, 15-28 mg/g) are chemically implausible: acid hydrolysis
   deamidates them to Glu and Asp. Either a non-standard hydrolysis or the bars are not what they say (P ≈ 75 % the
   totals for Gln / Asn are artefacts). Total Cys 16-18 mg/g (3.8 % of amino acids) is also high for yeast protein
   (typically ~1-1.5 %) [C, general knowledge].
4. [!] **Sums do not close**: individual free amino acids sum to 122.9 / 161.4 / 316.4 / 399.6 mg/g vs printed 115.7 /
   155.3 / 306.0 / 433.1 (lab extracts +4-6 %, industrial −8 %); totals 417.5 / 449.0 / 450.3 / 655.6 vs 440.1 / 449.3 /
   447.2 / 633.3 (CM −5 %, IP +4 %) [Z].
5. [!] Industrial product error bars are constant (±0.10 / ±0.2) and Fig. 8 total error bars are ±5.0 for all four
   extracts: nominal, not computed (P ≈ 80 %). The companion's Table 2 prints the same totals with ±1.2-2.5:
   446.13 / 449.16 / 447.02 / 633.31 — S, A, IP match Fig. 8 to ±0.2, but **CM is 440.1 here vs 446.13 there**
   (6 mg/g, beyond the ±1 read error).
6. [!] Mechanical-disruption temperature: methods say 7 °C, discussion says 10 °C.
7. [!] Free amino acids are 300× the cold-methanol intracellular pool [Q, ref. 11]; they are generated by autoproteolysis
   after disruption, so they depend on hold time and temperature before freeze-drying (analysed within 12 h; hold
   conditions before freezing not stated).

## 6. Existing records in the repo that cite this paper

`grep -rIn` 2026-11-03 for `00217-018-3143`, `jacob2018`, "Comparison of various industrially", "Jacob … 2018"
across `data/ docs/ src/ scripts/ tasks/`: no code, no ledger entry of its own; citations below.

| record | field / claim | repo value | paper | verdict |
|---|---|---|---|---|
| `jacob2019_extraction.md` ll. 6-8, 215-216 | citation | Jacob, Hutzler & Methner 2019, EFRT 245:95-109, 10.1007/s00217-018-3143-z | same (online 2018-08-23) | **match** |
| same | content "individual free and total amino acids", same yeast, same three methods | Figs. 3-7 | **match** |
| same, sec. 1 table | free amino acid sums 115.7-433.2 mg/g dw (its Table 2: 115.68 / 155.38 / 306.36 / 433.21) | Fig. 8: 115.7 / 155.3 / 306.0 / 433.1 | **match** (to the read error) |
| same, sec. 4 | FAN 16.50 / 25.80 / 45.10 / 58.1 | text 16.5 / 25.8 / 45.1 / 58.1 | **match** |
| same, sec. 4 (its Table 2) | total amino acids 446.13 / 449.16 / 447.02 / 633.31 | Fig. 8: 440.1 / 449.3 / 447.2 / 633.3 | **match except cell mill** (440.1 vs 446.13; an inter-paper discrepancy, not a repo error) |
| same, sec. 1 ans. 3 | autolysate "GSH" bar "may be partly Cys" | free Cys 82 µmol/g vs DTNB 3.9 µmol/g | **unlikely for the autolysates** (free "Cys" is mostly cystine or DTNB-silent); plausible for mechanical extracts |
| `data/lit/measurement_ledger/aroma_other.yml` l. 3093 (`jacob2019_tab2_tab3_spent_yeast_extract_composition` note) | "Individual amino acids in the companion paper 10.1007/s00217-018-3143-z (not on disk)" | DOI matches | **content match; "not on disk" now stale** |
| `tome2021_extraction.md` ll. 50-52 (ref. 31) | title, EFRT 245(1):95-109, DOI | same | **match** |
| `tome2021_extraction.md` l. 123 | "industrial autolysate free 433 / total 633; lab autolysis 306 / 447" (attributed to Jacob 2019 Table 2) | Fig. 8: 433.1 / 633.3; 306.0 / 447.2 | **match** (also printed in `jacob2019` Table 2) |
| `data/keys/laboratories.yml` l. 84 (`tu_berlin_methner`) | "Jacob et al. 2018 / 2019: Chair of Brewing Science, TU Berlin (Methner)" | Fachgebiet Brauwesen, TU Berlin; Methner senior | **match** |

## 7. Proposed ledger entries (`data/lit/measurement_ledger/aroma_other.yml`; not written)

```yaml
  # ---- jacob2018 (ingested 2026-11-03) ----
  - id: jacob2018_figs3_5_free_amino_acids_by_disruption_method
    dossier: jacob2018_extraction.md
    paper: "Jacob, Hutzler & Methner 2018 Eur Food Res Technol 245:95 (online 2018, issue 2019)"
    lab: "TU Berlin Brewing (Jacob / Methner)"
    system: "freeze-dried extracts of spent brewer's yeast TUM 68 (cell mill, sonotrode, autolysis 24 h 50 C) + one industrial autolysed yeast extract"
    medium: "dry powders (~10 % water), resuspended in water, analysed unheated within 12 h"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [Gln, His, Val, Ile, Leu, Pro, Cys, Asp, Glu, Asn, Ser, Gly, Thr, Ala, Arg, Tyr, Met, Trp, Phe, Lys, GABA, ornithine]
    kind: other
    location: "Figs. 3-5 bar charts (mg/g dry mass, mean +- 95 % CI, n = 3), digitised from vector paths in dossier secs. 4a-4b"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "The free amino acids (except proline and cysteine) were quantified using ... HPLC"
    flags: [matrix_precursor_composition, matrix_food, dry_weight_basis, figure_read, cys_redox_state_unknown, industrial_product_no_replicate_error, individual_sum_differs_from_printed_sum, "also_relevant: sulfur"]
    note: "Free Cys 4.0 / 6.1 / 10.0 / 13.0 and Met 3.0 / 5.2 / 9.1 / 11.3 mg/g dw (CM / S / A / IP); autolysis frees x2-5 more of most amino acids. Lab-autolysate free Cys (82 umol/g) is >20x the companion's DTNB thiol, so treat it as cyst(e)ine total, not reduced Cys."
  - id: jacob2018_figs6_8_tab2_total_amino_acids_and_protein
    dossier: jacob2018_extraction.md
    paper: "Jacob, Hutzler & Methner 2018 Eur Food Res Technol 245:95 (online 2018, issue 2019)"
    lab: "TU Berlin Brewing (Jacob / Methner)"
    system: "as jacob2018_figs3_5_free_amino_acids_by_disruption_method"
    medium: "dry powders, acid-hydrolysed for total amino acids; Bradford protein"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [total amino acids (19 individual, no Trp), protein (Bradford), free amino acids (sum), total amino acids (sum)]
    kind: other
    location: "Figs. 6-7 (individual totals), Fig. 8 (protein / free / total sums), Table 2 (essential amino acids per g total amino acid), dossier secs. 4c-4e"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "The resuspended yeast extracts underwent acid hydrolysis to determine all free and protein-bound amino acids"
    flags: [matrix_precursor_composition, matrix_food, dry_weight_basis, figure_read, total_gln_asn_after_acid_hydrolysis_implausible, tab2_met_cys_is_met_only, nominal_error_bars]
    note: "Totals do not differ by disruption method (ANOVA p > 0.05) except Glu; total amino acids 44-45 % dm (lab extracts), 63 % (industrial). Bradford protein 32 / 29 / 5 / 10 % dm misses < 3 kDa peptides in autolysates."
```

## 8. Supplementary material
None referenced (no ESM, no "Online Resource"). All model-relevant numbers are in Figs. 3-8 and Table 2 of the
main text. Sugars, RNA, thiamine and the Ellman thiol for the same extracts are in the companion `jacob2019`.
