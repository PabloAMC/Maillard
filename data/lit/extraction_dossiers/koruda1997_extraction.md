# Kuroda, Nagaba, Tsubuku, Kawajiri & Ueda 1997 — EXTRACTION (file stem `koruda1997`)
**Food Sci. Technol. Int. Tokyo 3(3):239-240** (Note, 2 pp.). DOI **10.3136/fsti9596t9798.3.239** (not printed on
the PDF; found by Crossref bibliographic query 2026-11-03: "Simultaneous Determination of Glutathione,
γ-Glutamylcysteine and Cysteine in Commercial Yeast Extract by HPLC with Fluorimetric Detection", Kuroda, Nagaba,
Tsubuku, Kawajiri, Ueda, 1997, vol. 3, pp. 239-240). **Matches the PDF** (title, authors, journal, pages).
[!] **The file stem misspells the first author: `koruda` for Kuroda.** The stem is kept so the coverage gate matches;
any ledger `paper:` string should say Kuroda.

### Four commercial glutathione-enriched yeast extracts (GSH-YE, brands A-D, Japan) analysed for reduced GSH, γ-Glu-Cys (γ-GC) and free Cys by aminopropyl-silica HPLC with post-column N-acridinylmaleimide (NAM) fluorescence. One table of contents (g/100 g product), one recovery table, LODs, CVs. **No heating, no sugars, no free amino acids, no nucleotides, no oxidised forms (GSSG, cystine), no dry-matter basis, no SD per brand.**

**Read 2026-11-03 from `data/articles/koruda1997.pdf`** (2 pp., scanned image with an OCR layer that garbles
numbers: "lO", "O.I O"). p. 1 prose from `pdftotext -layout`; **p. 2 (Tables 1-2, Fig. 1) read from a 200 dpi render**
(scratchpad `ingest/koruda1997_work/p-2.png`). All numbers below are from the render [M]; unit conversions [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Food Research Laboratories, **Ajinomoto Co., Ltd.**, Kawasaki; no corresponding author marked;
first author Motonaka Kuroda, senior (last) author Yoichi Ueda (Ueda also on the garlic / onion sulfur-taste papers
cited, and co-inventor with Kuroda of the cited GSH-YE seasoning patent). **Not in `data/keys/laboratories.yml`: NEW**,
proposed id **`ajinomoto_kuroda_ueda`**, string "Ajinomoto Food Research Labs (Kuroda / Ueda)". No same-lab dossier
in `data/lit/extraction_dossiers/` (grep "Ajinomoto", "Kuroda", "Ueda": none).

## 1. What this paper answers for the model

| model need (yeast extract as a forced input) | what Kuroda 1997 gives | tag |
|---|---|---|
| free cysteine in a YE | 4 brands, 0.020-1.71 g/100 g (×86 spread) | [M] |
| glutathione (reduced) | 4 brands, 0.46-5.02 g/100 g | [M] |
| γ-Glu-Cys (a second thiol peptide pool) | 4 brands, 0.057-0.50 g/100 g | [M] |
| GSSG, cystine (oxidised pools) | **not measured** (NAM labels free SH only) | [!] |
| free amino acids, reducing sugars, ribose, IMP/GMP, thiamine | **none** | [!] |
| basis (dry weight / as sold / moisture) | "g/100 g" of product, moisture not given | [!] |
| product type | **GSH-enriched** YE only (a speciality seasoning), not ordinary autolysed YE | [!] |
| method quality | LOD 2.3 / 4.5 / 1.1 pmol (GSH / γ-GC / Cys); recovery 95.6-98.8 %; CV 1.1-2.0 % | [M] |

**Answers.**
1. **Does it give cysteine and GSH in commercial YE? Yes, for four GSH-enriched products** (confidence ~95 % that the
   numbers are as printed and method-sound; recoveries and CVs are good). Molar (as sold) [Z]: GSH 15-163 µmol/g,
   γ-GC 2.3-20 µmol/g, Cys 1.7-141 µmol/g. GSH is 73-85 % of the molar reduced thiol in brands B-D, but only 50 %
   in brand A, where free Cys (141 µmol/g) nearly equals GSH (163 µmol/g) [Z].
2. **Is it a valid input for an ordinary (non-enriched) yeast extract? No** (confidence ~80 %). All four are sold as
   GSH-YE; brand D (0.46 % GSH) may be a weakly enriched or ordinary extract, but the paper does not say. Ordinary
   autolysed YE has ~0.1-0.2 % GSH-equivalent thiol (`jacob2019_extraction.md`, Fig. 4) and ~0.2 % free Cys
   (`munch1997`, CYE 2.13 mg/g dw).
3. **Can the brand spread be used as a prior band? Yes, as a product-to-product band for GSH-YE** (confidence ~70 %):
   Cys/GSH ratio 0.02-0.34, γ-GC/GSH 0.08-0.16 [Z]. The authors attribute the γ-GC and Cys variation to fermentation
   conditions [Q].
4. **Reduced-only is a bias of unknown size** (confidence ~60 % that GSSG + cystine are a material fraction in a dried
   powder after storage). The model would need total (reduced + oxidised) pools as charged; this paper bounds the
   reduced part only.

## 2. Conditions

- **Samples** [M]: commercial GSH-YE "purchased from four manufacturers (A-D)" [V]; Japanese market, 1996-97; brand
  names, form (powder / paste), moisture, salt and lot not given.
- **Preparation** [M]: ~10 mg dissolved in 10 mL water; 3.2 mL + 0.4 mL 50 % (w/v) metaphosphoric acid + 0.4 mL
  0.2 M EDTA, 10 min room T, 0.45 µm filter; 2.5 mL filtrate + 1.5 mL 0.8 N KOH; 20 µL injected.
- **No heating, no pH of the product, no temperature, no replicates per brand stated** (Table 2 gives single values).

## 3. Analytical methods

- **HPLC** [M]: Finesil-NH2 10 µm, 250 × 4.6 mm (JASCO), 25 °C; isocratic 0.05 M K phosphate pH 3.5 / acetonitrile
  90:10, 1.0 mL/min. Post-column 0.045 M KOH (0.3 mL/min) then NAM 2 µg/mL in acetonitrile (0.3 mL/min), 10 m × 0.25 mm
  coil at 50 °C; fluorescence ex 365 / em 435 nm. Retention: Cys 6.3, GSH 9.1, γ-GC 10.7 min.
- **Quantification** [M]: external standards (10 µM each in 20 mM EDTA; γ-GC from Kojin Co.); linear r = 0.998 for all
  three, 20 pmol-2 nmol. Absolute, calibrated (quant_class `calibrated`).
- **LOD** [M]: GSH 2.3, γ-GC 4.5, Cys 1.1 pmol. **Precision**: CV 1.1 / 1.5 / 2.0 % (GSH / γ-GC / Cys), "n=95" [!]
  (probably n = 5 misprinted, or 95 repeat injections; not resolvable).
- **What the ± means**: Table 2 carries no ±. Table 1 column 4 header reads "Standard of recovery (%)" [V] —
  truncated; likely the RSD of the recovered amount.

## 4. Data

### 4a. Table 2 (p. 240): contents in commercial GSH-YE, g/100 g [M]; molar [Z]

| brand | GSH | γ-GC | Cys | GSH µmol/g [Z] | γ-GC µmol/g [Z] | Cys µmol/g [Z] |
|---|---|---|---|---|---|---|
| A | 5.02 | 0.50 | 1.71 | 163.3 | 20.0 | 141.1 |
| B | 5.00 | 0.48 | 0.12 | 162.7 | 19.2 | 9.9 |
| C | 1.52 | 0.25 | 0.10 | 49.5 | 10.0 | 8.3 |
| D | 0.46 | 0.057 | 0.020 | 15.0 | 2.28 | 1.65 |

[Z] MW 307.32 (GSH), 250.27 (γ-GC), 121.16 (Cys); basis "g/100 g" of the product as sold (dry-matter basis not
stated). Text range "0.057-0.50" and "0.020-1.71 g/100 g" matches the table.

### 4b. Table 1 (p. 240): recovery in yeast-extract samples, n = 5 per level [M]

| analyte | added (pmol) | mean recovered (pmol) | "standard of recovery" (%) |
|---|---|---|---|
| GSH | 40 / 100 / 200 | 39 / 98 / 199 | 1.1 / 0.9 / 1.2 |
| γ-GC | 40 / 100 / 200 | 38 / 98 / 192 | 1.3 / 1.0 / 1.1 |
| Cys | 40 / 100 / 200 | 39 / 95 / 189 | 1.7 / 1.3 / 1.5 |

Text: average recoveries 98.8 / 96.8 / 95.6 % [M].

### 4c. Fig. 1: chromatograms (standards 10 µM each; brand C sample). Qualitative only; not digitised.

## 5. Flags and internal consistency

1. [!] **Recovery averages do not quite reproduce from Table 1** [Z]: GSH 97.5 / 98 / 99.5 → mean 98.3 % (printed
   98.8); γ-GC 95 / 98 / 96 → 96.3 % (96.8); Cys 97.5 / 95 / 94.5 → 95.7 % (95.6). Within the rounding of the
   integer pmol values; harmless.
2. [!] **No replicate count or ± for Table 2**; brand values are single numbers.
3. [!] **Basis**: "g/100 g" with no moisture; if the products are ~5 % water (dry powder) the dry-basis values are
   ~5 % higher; pastes would differ more.
4. [!] **Reduced thiols only**: NAM is a maleimide (labels free SH); GSSG, mixed disulfides and cystine are invisible.
   Total cysteine-equivalent sulfur in these products is therefore ≥ the sum shown.
5. [!] "n=95" in the precision sentence (sec. 3).
6. [!] Stem misspelling (header).

## 6. Existing records in the repo that cite this paper

`grep -rIn` 2026-11-03 for `koruda`, `kuroda`, `fsti9596`, "Kawajiri", "Nagaba" across `data/ docs/ src/ scripts/
tasks/`: **only `data/keys/pdf_dossier_aliases.yml`** (the pending triage line, which says "Kuroda, Nagaba, Tsubuku,
Kawajiri & Ueda 1997 Food Sci Technol Int Tokyo 3:239 (2 pp.)": **matches**, and correctly notes the misspelt stem).
No ledger, registry or code record cites it.

**What it adds to existing yeast-extract composition in the repo.** Existing: `munch1997_tab3_7_9_yeast_free_amino_acids`
(15 free amino acids + free Cys in one commercial baker's-YE and lab fractions; ninhydrin analyser, Cys as
pyridylethyl-Cys), `munch1998_tab6_7_yeast_amino_acids_fermentation` (6 amino acids), `raza2020_tab5_fig5_fig7_ye_precursors`
(a GSH-YE "YE-GSH-10%" whose GSH and free Cys are **never measured**; nucleotides / ribose figure-only). **Kuroda is
the first record in the repo of GSH and γ-GC contents of any yeast extract**, and the first free-Cys numbers for
GSH-enriched products: it supplies the missing composition of the class of product Raza 2020 heated (a "10 % GSH"
label vs 0.46-5.02 % measured here across brands, a 2-10× gap worth noting [Z]).

## 7. Proposed ledger entries (aroma_other.yml; not written)

```yaml
  # ---- koruda1997 : Kuroda, Nagaba, Tsubuku, Kawajiri & Ueda 1997 (thiols in commercial GSH-enriched yeast extracts) (ingested 2026-11-03) ----
  - id: koruda1997_tab2_gsh_yeast_extract_thiols
    dossier: koruda1997_extraction.md
    paper: "Kuroda, Nagaba, Tsubuku, Kawajiri & Ueda 1997 Food Sci Technol Int Tokyo 3:239 (file stem koruda1997)"
    lab: "Ajinomoto Food Research Labs (Kuroda / Ueda)"
    system: "four commercial glutathione-enriched yeast extracts (brands A-D, Japan), unheated"
    medium: "product dissolved in water (about 1 mg/mL), metaphosphoric acid / EDTA, neutralised with KOH"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [glutathione, gamma-glutamylcysteine, cysteine]
    kind: other
    location: "Table 2 (g/100 g; transcribed with molar conversions in dossier sec. 4a); Table 1 recoveries (sec. 4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "commercial GSH-YE purchased from four manufacturers (A-D)"
    flags: [matrix_precursor_composition, matrix_food, reduced_thiols_only, basis_as_sold_moisture_not_reported, no_sd_per_brand, gsh_enriched_products_only, "also_relevant: sulfur"]
    note: "GSH 0.46-5.02, gamma-Glu-Cys 0.057-0.50, free Cys 0.020-1.71 g/100 g product (HPLC, post-column NAM fluorescence, external standards, recovery 95.6-98.8 %). Reduced forms only (no GSSG, cystine). Product-to-product band for GSH-enriched YE (the class Raza 2020 heated), not for ordinary autolysed YE."
```

## 8. Supplementary material

None referenced (a two-page Note). The cited earlier GSH-YE analyses (Aoyagi & Miya 1991, Gekkan Food Chemical 78:70,
in Japanese; Hashimoto et al. 1992, J Jpn Soc Nutr Food Sci 45:363) report GSH only and are not on disk; not needed.
