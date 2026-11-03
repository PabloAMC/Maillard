# Jacob, Striegel, Rychlik, Hutzler & Methner 2019 — EXTRACTION (file stem `jacob2019`)
**Eur. Food Res. Technol. 245:1169-1182**, DOI **10.1007/s00217-019-03237-9** (printed p. 1). Crossref 2026-11-03:
"Yeast extract production using spent yeast from beer manufacture: influence of industrially applicable disruption
methods on selected substance groups with biotechnological relevance", Jacob, Striegel, Rychlik, Hutzler, Methner,
issued 2019-01-24, EFRT 245, 1169-1182. **Matches the PDF and the brief.**
**Not to be confused with** the companion Jacob, Hutzler & Methner 2019, EFRT 245:95-109,
10.1007/s00217-018-3143-z (same yeast, same three disruption methods, **individual free and total amino acids**);
not on disk (cited here as ref. 15 and by Tomé 2021 as ref. 31).

### Spent top-fermenting brewer's yeast (S. cerevisiae TUM 68, 50 L pilot CCT, washed ×3) disrupted three ways at 200 mL, 7 % dry matter: Vibrogen cell mill (0.5 mm glass beads, 15 min, 7 °C), ultrasonic sonotrode (400 W, 30 min, 7 °C), autolysis (24 h, 50 °C, NaCl 0.086 M + ethyl acetate 0.051 M); centrifuged, freeze-dried; plus one **industrial autolysed yeast extract** of unknown starting yeast. Composition per g dry weight: protein (4 ways), total / free amino acids (sums only), FAN, glycogen, β-glucan, trehalose, glucose, fructose, DNA, RNA, ash, fat (Table 2); eight B vitamins incl. thiamine (Table 3); minerals (Table 4); glutathione and polyphenols (Fig. 4, figure only). **No heating of the extracts, no individual amino acids, no cysteine, no ribose, no 5'-nucleotides, no pH.**

**Read 2026-11-03 from `data/articles/jacob2019.pdf`** (14 pp., Springer, born-digital). Prose from
`pdftotext -layout`; **Table 2 (p. 1174) and Table 3 + Fig. 2 (p. 1176) read from 170 dpi renders**
(scratchpad `ingest/jacob2019_work/p-06.png`, `p-08.png`); **Fig. 4 (p. 1180) read from a 170 dpi render and a 400 dpi
crop** (`p-12.png`, `fig4-12.png`), read error ±0.15 mg/g on the left axis (2 mg/g gridlines ≈ 82 px at 400 dpi).
Table 4 (minerals) seen in the text layer only; not transcribed.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** **TU Berlin, Fachgebiet Brauwesen** (Institut für Lebensmitteltechnologie und Lebensmittelchemie);
corresponding author **F. F. Jacob**, senior author **F.-J. Methner**; folate vitamers by Striegel & Rychlik (TU
München, Analytische Lebensmittelchemie), yeast strain / PCR by Hutzler (TUM Weihenstephan). **Not in
`data/keys/laboratories.yml`: NEW**, proposed id **`tu_berlin_methner`**, string "TU Berlin Brewing (Jacob / Methner)".
No same-lab dossier (grep "Methner", "Jacob", "Berlin" in `data/lit/extraction_dossiers/`: none).

## 1. What this paper answers for the model

| model need (yeast extract as a forced input) | what Jacob 2019 gives | tag |
|---|---|---|
| free amino acids | **sum only** (115.7-433.2 mg/g dw) and FAN; individual amino acids in the companion paper (not on disk) | [M] / [!] |
| cysteine | measured by a separate method [Q] but **not reported here** | [!] |
| glutathione | Fig. 4 bars, 1.2-13.2 mg/g dw (Ellman's reagent: total free thiol, expressed as GSH) | [F] |
| reducing sugars | glucose 1.73-3.19, fructose 0.47-0.98 mg/g dw (HPLC); trehalose 25-44 (non-reducing) | [M] |
| ribose, IMP / GMP | **not measured**; only total RNA 48.8-55.2 mg/g dw (A260) | [!] |
| thiamine | 5.18-7.46 mg/100 g dw (HPLC, DIN EN 14122) | [M] |
| Fe / Cu (thiol-oxidation catalysts) | Table 4 exists (ICP-OES); not transcribed | — |
| basis | dry weight; moisture 6.8-11.3 % printed | [M] |
| heating | **none**; extracts characterised as made | [!] |

**Answers.**
1. **Composition of a real (autolysed) yeast extract as model input? Partly** (confidence ~75 % the printed values
   are sound; ± are 95 % confidence limits of triplicates). For the industrial autolysed product: free amino acids
   433 mg/g dw, glucose 1.73 + fructose 0.47 mg/g (≈ 9.6 + 2.6 µmol/g [Z]), thiamine 7.46 mg/100 g (≈ 0.28 µmol/g as
   thiamine cation [Z]), GSH-equivalent thiol ≈ 1.9 mg/g (≈ 6 µmol/g [Z]). Ribose and nucleotides are absent, so the
   pentose / IMP side of the forced input stays unfilled.
2. **Does the disruption method matter for the Maillard precursors? Yes, strongly for GSH and free amino acids**
   (confidence ~85 %). Autolysis (24 h, 50 °C) leaves ~1.2 mg/g GSH-equivalent vs ~13 mg/g after mechanical
   disruption (×11 [F]) and raises free amino acids ×2-2.6 (306 vs 116-155 mg/g dw [M]); reducing sugars do not
   differ (ANOVA p > 0.05 [M]). Commercial (autolysed) YE therefore delivers Maillard amines mostly as free amino
   acids and little GSH; a mechanically made YE delivers GSH but fewer free amino acids.
3. **Is the "glutathione" GSH? Only nominally** (confidence ~70 % that it is total free thiol). The Ellman / DTNB
   assay reads every free SH (GSH, γ-Glu-Cys, Cys, protein thiols) at 412 nm; the paper reports it as glutathione.
   For autolysed extracts, where free Cys rises [Q, ref. 15], the "GSH" bar may be partly Cys.
4. **Reducing sugars are low** (confidence ~90 %): ~0.3 % glucose and 0.09 % fructose on average [M text], consistent
   with spent yeast at the end of fermentation. A YE alone is sugar-poor; reducing sugar in a YE-based flavour comes
   from added ribose / xylose / glucose or from RNA-derived ribose (not measured here).

## 2. Conditions

- **Yeast** [M]: S. cerevisiae LeoBavaricus TUM 68 (top-fermenting), propagated in standardised hopped wort
  (12.4 °P, pH 5.19, Table 1), fermented 5 d at 18 °C in a 50 L CCT, cropped day 6 (70 % live), washed three times
  (10 % dry matter, 1000 g), diluted to **7 % dry matter**.
- **Disruption** [M]: cell mill (Vibrogen VI 6, 75 Hz, 200 mL, 0.5 mm glass beads 1:2 w/w, 15 min, jacket-cooled to
  7 °C); sonotrode (Bandelin HD 3400, 25 mm, 20 kHz, 400 W continuous, 30 min, 7 °C); autolysis (200 mL, 24 h, 50 °C,
  100 rpm, NaCl 0.086 mol/L + ethyl acetate 0.051 mol/L).
- **Extract** [M]: centrifuged 20 min, 10,000 g, 4 °C; supernatant freeze-dried (−50 °C, 0.1 mbar) to a powder with
  ~10 % water [Q] (Table 2: 9.4-11.3 %).
- **Industrial product** [M]: a commercial autolysed YE; starting yeast and medium unknown [Q].
- **Replicates** [M]: "All experiments were performed in triplicate" [V]; ± = 95 % confidence limits (Student t);
  Dixon outlier test.
- **No heating, pH, water activity or reaction step**: a composition paper.

## 3. Analytical methods

- Protein: Bradford (595 nm, > 3 kDa); Kjeldahl N × 6.25 / × 5.5 (MEBAK 2.6.1.1). Free amino acids: HPLC (MEBAK
  2.6.4.1) except Pro and Cys (separate method, VDLUFA 4.11.1); total amino acids after acid hydrolysis. FAN: ninhydrin
  570 nm [M]. **Sums only are printed.**
- Sugars: glucose, fructose, maltose, sucrose, maltotriose by HPLC (MEBAK 2.7.1); glycogen by amyloglucosidase +
  HPLC glucose; trehalose by anthrone (620 nm); β-glucan enzymatic [M]. Calibrated, absolute.
- Thiamine: HPLC (DIN EN 14122); other B vitamins by DIN / AOAC / microbiological assays; folates by LC-MS/MS [M].
- DNA / RNA: A260 (NanoDrop) after kit purification [M] — total RNA, no nucleotide speciation.
- **Glutathione: Ellman's reagent (DTNB), 412 nm**, after the Ortiz-Julien patent method [M]; non-specific for GSH [!].
- Polyphenols: Fe(III) complex 600 nm (MEBAK 2.16.1); antioxidant potential: Sigma kit (Trolox eq.); reduction
  potential: DCPIP 520 nm [M].

## 4. Data

### 4a. Table 2 (p. 1174): general composition, mg/g dw, mean ± 95 % CL (n = 3) [M]

| quantity | cell mill | sonotrode | autolysis | industrial (autolysed) |
|---|---|---|---|---|
| protein N × 6.25 | 485.24 ± 8.90 | 487.75 ± 10.49 | 482.01 ± 2.06 | 743.10 ± 0.5 |
| protein N × 5.50 | 427.02 ± 9.53 | 429.22 ± 9.21 | 424.24 ± 2.98 | 653.93 ± 0.5 |
| protein Bradford | 321.56 ± 9.80 | 285.40 ± 8.90 | 52.90 ± 12.50 | 102.00 ± 1.00 |
| total amino acids | 446.13 ± 1.20 | 449.16 ± 2.00 | 447.02 ± 2.51 | 633.31 ± 2.00 |
| **free amino acids** | 115.68 ± 1.23 | 155.38 ± 1.26 | 306.36 ± 2.65 | 433.21 ± 0.10 |
| FAN | 16.50 ± 0.21 | 25.80 ± 0.32 | 45.10 ± 0.28 | 58.1 ± 0.10 |
| glycogen | 263.56 ± 10.21 | 266.87 ± 12.92 | 258.38 ± 8.63 | 108.48 ± 2.55 |
| β-glucan | 1.45 ± 0.80 | 1.38 ± 0.01 | 1.36 ± 0.01 | 1.40 ± 0.01 |
| trehalose | 44.00 ± 1.90 | 42.10 ± 2.10 | 25.00 ± 3.60 | 35.20 ± 0.10 |
| **glucose** | 2.57 ± 0.23 | 3.19 ± 0.45 | 2.82 ± 0.15 | 1.73 ± 0.01 |
| **fructose** | 0.98 ± 0.16 | 0.75 ± 0.21 | 0.87 ± 0.17 | 0.47 ± 0.01 |
| DNA | 20.56 ± 0.21 | 19.67 ± 0.47 | 2.30 ± 0.54 | 3.49 ± 0.10 |
| **RNA** | 50.01 ± 0.21 | 50.07 ± 0.37 | 48.84 ± 0.87 | 55.18 ± 0.10 |
| ash | 130.94 ± 2.56 | 133.37 ± 2.47 | 131.86 ± 2.49 | 135.87 ± 2.30 |
| fat | 10.3 ± 0.5 | 11.4 ± 0.6 | 5.1 ± 0.3 | 6.7 ± 0.01 |
| moisture (%) | 9.50 ± 0.12 | 9.40 ± 0.11 | 11.3 ± 0.30 | 6.80 ± 0.10 |

Maltose, sucrose, maltotriose: methods given, values not printed [!].

### 4b. Table 3 (p. 1176): B vitamins, mg/100 g dw, mean ± 95 % CL [M] (model-relevant row first)

| vitamin | cell mill | sonotrode | autolysis | industrial |
|---|---|---|---|---|
| **thiamine (B1)** | 6.88 ± 0.31 | 7.05 ± 0.43 | 5.18 ± 0.20 | 7.46 ± 0.5 |
| riboflavin (B2) | 2.16 ± 0.25 | 2.41 ± 0.33 | 1.16 ± 0.19 | 10.55 ± 0.5 |
| niacin (B3) | 94.19 ± 1.21 | 103.62 ± 2.15 | 68.34 ± 0.98 | 78.58 ± 2.0 |
| pantothenic acid (B5) | 20.36 ± 0.54 | 18.56 ± 0.78 | 15.73 ± 0.28 | 43.18 ± 2.0 |
| pyridoxine (B6) | 4.86 ± 0.28 | 5.10 ± 0.51 | 3.09 ± 0.19 | 5.90 ± 0.5 |
| biotin (B7) | 113.92 ± 2.67 | 127.94 ± 3.69 | 138.55 ± 1.58 | 618.65 ± 2.0 |
| total folate (B9) | 4.52 ± 0.26 | 4.94 ± 0.42 | 1.35 ± 0.18 | 5.29 ± 0.1 |
| cobalamin (B12) | 0.18 ± 0.02 | 0.12 ± 0.04 | 1.12 ± 0.05 | 0.16 ± 0.1 |

[!] Biotin and B12 at "mg/100 g" are implausibly high for yeast (biotin is normally µg/100 g); the unit applies to
the whole table as printed. Thiamine is in the expected range. Fig. 2 (folate vitamers) not digitised (not model-relevant).

### 4c. Fig. 4 (p. 1180): glutathione (Ellman) and polyphenols, mg/g dw, stacked bars [F, ±0.15]

| | cell mill | sonotrode | autolysis | industrial |
|---|---|---|---|---|
| "glutathione" (lower bar) | 12.9 | 13.2 | 1.2 | 1.9 |
| bar top (GSH + polyphenols) | 16.0 | 16.6 | 2.7 | 3.0 |
| polyphenols (difference) [Z] | 3.1 | 3.4 | 1.5 | 1.1 |
| GSH-equivalent µmol/g [Z, MW 307.3] | 42 | 43 | 3.9 | 6.2 |

Error bars on GSH ≈ ±0.4 (cell mill), ±0.6 (sonotrode) mg/g [F]. Antioxidant potential (right axis, mg TE/g dw) ≈ 23 / 25 / 1 / 2
and reduction potential (right axis, %) ≈ 91 / 92 / 4 / 4 [F, ±2]; ANOVA p < 0.01 for all [M]. (The diamonds sit near
4 on the left axis but belong to the right axis: the text's "around ten times" above Vieira's 1.3-2.6 mg TE/g fits ~23.)

## 5. Flags and internal consistency

1. [!] **"Glutathione" by DTNB is total free thiol** (sec. 3); no GSH-specific assay; no GSSG.
2. Text averages check out [Z]: glucose mean of three lab extracts 2.86 mg/g = 0.29 % (text "0.3 %"); fructose 0.87
   mg/g = 0.087 % ("0.09 %"); trehalose gain 17.1-19 mg/g ("more than 17"); glycogen ~26 %; ash 13.2-13.3 %.
3. [!] **Autolysis mass balance**: free (306) + Bradford (53) = 359 ≪ total amino acids 447 mg/g; the authors put the
   gap in peptides < 3 kDa [Q]. Mechanical extracts close within ~2 % (116 + 322 = 437 vs 446).
4. [!] **Industrial-product ± are tiny** (±0.001 on minerals, ±0.01 on sugars): probably analytical replicates of one
   lot, not process replicates; not comparable with the lab extracts' CL.
5. [!] Moisture "10 %" in the methods vs 6.8-11.3 % in Table 2; values are on a dry-weight basis, so harmless.
6. [!] Biotin / B12 units (sec. 4b).
7. [!] Cys and Pro were analysed (methods) but are not reported anywhere in this paper.

## 6. Existing records in the repo that cite this paper

`grep -rIn` 2026-11-03 for `jacob2019`, `00217-019-03237`, `00217-018-3143`, "Jacob FF", "Jacob, F", "Methner",
"TUM 68" across `data/ docs/ src/ scripts/ tasks/`: **only `data/keys/pdf_dossier_aliases.yml`** (pending triage line:
"Jacob et al. 2019 Eur Food Res Technol 245:1169 (10.1007/s00217-019-03237-9): yeast extract from spent brewer's yeast vs
cell-disruption method (composition)": **matches**). No ledger, registry or code record cites it.

**What it adds to existing yeast-extract composition.** `munch1997_tab3_7_9_yeast_free_amino_acids` and
`munch1998_tab6_7_*` give individual free amino acids and free Cys for baker's-yeast extracts, **but no sugars, RNA,
thiamine or GSH**; `raza2020_tab5_fig5_fig7_ye_precursors` has nucleotides / ribose for a heated GSH-YE (figure only).
Jacob adds, for the first time in the repo: **reducing sugars (glucose, fructose) of yeast extracts, thiamine, total RNA,
a GSH-equivalent thiol, and the free / total amino-acid split** for brewer's-yeast extracts and one commercial
autolysate, plus the disruption-method effect (autolysis destroys ~90 % of GSH and doubles free amino acids). It does
not give ribose, IMP / GMP or Cys.

## 7. Proposed ledger entries (aroma_other.yml; not written)

```yaml
  # ---- jacob2019 : Jacob, Striegel, Rychlik, Hutzler & Methner 2019 (spent brewer's yeast extract vs disruption method) (ingested 2026-11-03) ----
  - id: jacob2019_tab2_tab3_spent_yeast_extract_composition
    dossier: jacob2019_extraction.md
    paper: "Jacob, Striegel, Rychlik, Hutzler & Methner 2019 Eur Food Res Technol 245:1169"
    lab: "TU Berlin Brewing (Jacob / Methner)"
    system: "freeze-dried extracts of spent brewer's yeast TUM 68 (cell mill, sonotrode, autolysis 24 h 50 C) + one industrial autolysed yeast extract"
    medium: "dry powders, analysed unheated (moisture 6.8-11.3 %)"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [free amino acids (sum), total amino acids (sum), free alpha-amino nitrogen, glucose, fructose, trehalose, glycogen, RNA, thiamine]
    kind: other
    location: "Table 2 (mg/g dw) and Table 3 (mg/100 g dw), mean +- 95 % CL, n = 3 (transcribed in dossier secs. 4a, 4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "a yeast extract powder was produced ... that had a residual water content of 10%"
    flags: [matrix_precursor_composition, matrix_food, dry_weight_basis, amino_acids_sums_only, no_ribose_no_nucleotides, industrial_product_origin_unknown, "biotin_b12_units_implausible"]
    note: "Reducing sugars low (glucose 1.7-3.2, fructose 0.5-1.0 mg/g dw); thiamine 5.2-7.5 mg/100 g dw; RNA 49-55 mg/g dw total (no 5'-nucleotides). Autolysis doubles free amino acids (306 vs 116-155 mg/g). Individual amino acids in the companion paper 10.1007/s00217-018-3143-z (not on disk)."
  - id: jacob2019_fig4_glutathione_by_disruption_method
    dossier: jacob2019_extraction.md
    paper: "Jacob, Striegel, Rychlik, Hutzler & Methner 2019 Eur Food Res Technol 245:1169"
    lab: "TU Berlin Brewing (Jacob / Methner)"
    system: "as jacob2019_tab2_tab3_spent_yeast_extract_composition"
    medium: "dry powders, analysed unheated"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [glutathione (Ellman total thiol)]
    kind: other
    location: "Fig. 4 stacked bars (digitised, +-0.15 mg/g, dossier sec. 4c)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "after reacting with Ellman's reagent ... the absorption was measured spectrophotometrically at 412 nm"
    flags: [matrix_precursor_composition, matrix_food, dry_weight_basis, nonspecific_thiol_assay, figure_read, "also_relevant: sulfur"]
    note: "GSH-equivalent thiol 12.9 / 13.2 mg/g dw (cell mill / sonotrode) vs 1.2 (lab autolysis) and 1.9 (industrial autolysate): autolysis loses ~90 %. DTNB reads all free SH (GSH, gamma-Glu-Cys, Cys), so an upper bound on reduced GSH."
```

## 8. Supplementary material

None referenced (no ESM, no Supplementary Table). The long-chain fatty-acid data are "not shown in detail" [Q]. The
model-relevant missing item is not SI but the **companion paper** (individual free amino acids incl. Cys by disruption
method): Jacob, Hutzler & Methner 2019, EFRT 245:95-109, DOI 10.1007/s00217-018-3143-z.
