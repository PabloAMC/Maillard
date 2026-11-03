# Tomé 2021 — EXTRACTION (file stem `tome2021`; review)
**ACS Food Sci. Technol. 1(4):487-494**, DOI **10.1021/acsfoodscitech.0c00131** (printed p. 487). Crossref 2026-11-03:
"Yeast Extracts: Nutritional and Flavoring Food Ingredients", Tomé, issued 2021-04-28, ACS Food Science & Technology
1, 487-494. **Matches the PDF and the brief.** Narrative review; received 5 Dec 2020, published 28 Apr 2021.

### Narrative review of yeast extracts (production by autolysis / plasmolysis / hydrolysis; nutrition; flavour; glutamate). Two composition tables compiled from secondary and grey sources: Table 1 (proximate, B vitamins, minerals for dried yeast, dried yeast extract, yeast-extract paste) and Table 2 (18 free and total amino acids, g/100 g DM, dried YE). Table 3 (umami-peptide taste thresholds) and Table 4 (aroma-active compounds, descriptors only) are qualitative for the model. **No glutathione, no 5'-nucleotides (IMP / GMP), no ribose, no free cysteine number; no new data.**

**Read 2026-11-03 from `data/articles/tome2021.pdf`** (8 pp., ACS, born-digital). Prose and reference list from
`pdftotext -layout`; **Tables 1 (p. 488) and 2 (p. 489) read from 170 dpi renders** (scratchpad
`ingest/tome2021_work/p-2.png`, `p-3.png`). Tables 3-4 seen in the text layer only (not model-relevant). No deep
extraction, per the brief.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag. In a review every number is [C] in substance; [M] below means "printed in this review".

**Laboratory.** Sole and corresponding author **Daniel Tomé**, UMR PNCA, AgroParisTech / INRAE, Université
Paris-Saclay. **Not in `data/keys/laboratories.yml`: NEW** if keyed, proposed id **`agroparistech_tome`** (review;
no laboratory data of its own, so it should never be a leave-one-laboratory-out unit — better to key the entry to the
cited primaries' labs when those are ingested). No same-author dossier (grep "Tom", "AgroParisTech": none).

## 1. What this paper answers for the model

| model need (yeast extract as a forced input) | what Tomé 2021 tabulates | traced to (cited primary, by title) | tag |
|---|---|---|---|
| free amino acids (18) | Table 2, free and total, g/100 g DM, one "typical" dried YE | refs 27, 28, 30, 31 (table-level, not per value; below) | [M]/[C] |
| free cysteine | **only "cystine" 0.3 free / 0.4 total g/100 g DM** (Table 2) | same | [M]/[!] |
| glutathione | **no number**; named as a component only | refs 27, 32-35 (named) | [Q] |
| sugars | Table 1: "sugar g" 8.0 (dried yeast), blank (dried YE), 0 (YE paste); carbohydrate 7-13 g (dried YE) | refs 5, 8, 10, 26-29 | [M]/[C] |
| ribose, IMP / GMP | **none** (nucleotides named, never quantified) | — | [!] |
| thiamine | Table 1: 40 (dried yeast), 10-12 (dried YE), 9.7 (YE paste) mg/100 g; text: "about 3 mg" per 100 g "standard yeast extract" | Table 1 refs 5, 8, 10, 26-29; text ref 27 | [M]/[C] |
| glutamate | free ~5 % and total ~10 % w/w (abstract, text); Table 2 Glu 5.2 free / 8.1 total | refs 27-31 | [M]/[C] |

**Cited primaries behind the tables (titles from the reference list) [C]:**
- ref 5: Sommer 1998, "Yeast extracts: production, properties and components", Food Australia 50(4):181-183.
- ref 8: Podpora & Świderski 2010, "Preparaty żywnościowe otrzymywane z odpadowych drożdży pofermentacyjnych jako
  przykład innowacji", Postępy Techniki Przetwórstwa Spożywczego 20:94-99 (Polish).
- ref 10: Podpora & Waszkiewicz-Robak 2007, "Nowoczesne metody otrzymywania ekstraktów drożdżowych jako naturalnych
  składników smakowych", Postępy Techniki Przetwórstwa Spożywczego, 67-71 (Polish).
- ref 26: Podpora & Świderski 2010, "Możliwość wykorzystania drożdży pofermentacyjnych w produkcji autolizatów
  drożdżowych", Przemysł Fermentacyjny i Owocowo-Warzywny 3(54):14 (Polish).
- ref 27: Podpora, Świderski, Sadowska, Piotrowska & Rakowska 2015, "Spent brewer's yeast autolysates as a new and
  valuable component of functional food and dietary supplements", J Food Process Technol 6(12).
- ref 28: Podpora, Świderski, Sadowska, Rakowska & Wasiak-Zys 2016, "Spent brewer's yeast extracts as a new component
  of functional food", Czech J Food Sci 34:554-563.
- ref 29: CIQUAL, "Composition nutritionnelle des aliments", 2012 (French food-composition table; the YE-paste column,
  with sodium 3.6 g/100 g, looks like a Marmite-type spread [P]).
- ref 30: Rozmierska, Mrówka, Komorowska & Stecka 2001, "Wartość żywieniowa białek preparatów pochodzenia drożdżowego",
  Żywność Nauka Technologia Jakość suppl. 2(08):144 (Polish).
- ref 31: Jacob, Hutzler & Methner 2019, "Comparison of various industrially applicable disruption methods to produce
  yeast extract using spent yeast from top-fermenting beer production: influence on amino acid and protein content",
  EFRT 245(1):95-109 (10.1007/s00217-018-3143-z; companion of `jacob2019`, not on disk).

**Answers.**
1. **Usable as a forced-input composition? Only as a plausibility band** (confidence ~80 %). Values are "typical",
   compiled across sources without saying which source gives which number, with no n, no ±, and no production method
   per column. A primary (Münch 1997 for baker's YE; Jacob 2019 for brewer's YE; the Podpora papers if ingested) should
   be the input; Tomé's Table 2 is a cross-check.
2. **Does it fill the ribose / nucleotide / GSH gap? No** (confidence ~99 %). Nothing quantitative on any of them.
3. **Is Table 2 internally sound? Not quite** (confidence ~90 %): rows do not add to the printed totals, and the text's
   free-fraction claims contradict the table (sec. 5).

## 2. Conditions

Not applicable (review). The tables give no production method, strain, lot, moisture (other than Table 1's "water g"
row: 4 / 4 / 37 g/100 g) or analytical method.

## 3. Analytical methods

None stated. Basis: Table 1 "for 100 g" of product as sold; Table 2 "g/100 g DM". Quant class of the underlying
data unknown.

## 4. Data

### 4a. Table 2 (p. 489): "typical" free and total amino acids of dried S. cerevisiae yeast extract, g/100 g DM [M]

| amino acid | free | total | | amino acid | free | total |
|---|---|---|---|---|---|---|
| alanine | 3.7 | 4.4 | | lysine | 2.2 | 4.0 |
| arginine | 1.8 | 2.5 | | methionine | 0.7 | 0.9 |
| aspartic acid | 1.9 | 4.9 | | phenylalanine | 1.8 | 2.3 |
| cystine | 0.3 | 0.4 | | proline | 1.0 | 2.0 |
| glutamic acid | 5.2 | 8.1 | | serine | 1.8 | 2.3 |
| glycine | 1.2 | 2.4 | | threonine | 1.7 | 2.1 |
| histidine | 0.9 | 1.0 | | tyrosine | 1.2 | 1.4 |
| isoleucine | 2.2 | 2.7 | | tryptophan | 0.6 | 0.6 |
| leucine | 3.5 | 3.8 | | valine | 2.6 | 2.9 |
| **total (printed)** | **34.6** | **48.0** | | | | |

No Asn, Gln, ornithine, GABA rows (acid hydrolysis would merge Asn / Gln into Asp / Glu; their absence from the free
column is a loss for the model) [!].

### 4b. Table 1 (p. 488): per 100 g product, model-relevant rows only [M]

| row | dried yeast | dried yeast extract | yeast extract paste |
|---|---|---|---|
| water g | 4 | 4 | 37 |
| protein (N × 6.25) g | 48 | 69 | 27.8 |
| carbohydrate g | 23.7 | 7-13 | 11.8 |
| – sugar g | 8.0 | (blank) | 0 |
| fat g | 5.4 | < 0.5 | 0.0 |
| B1 (thiamine) mg | 40 | 10-12 | 9.7 |
| iron mg | 5 | (blank) | 3.7 |
| copper mg | 5.3 | (blank) | 0.3 |
| sodium mg | 40 | (blank) | 3600 |

Other rows (kcal, fibre, fatty-acid classes, B2-B12, Ca, P, Mg, K, Zn, Mn, Se) exist; not transcribed.

## 5. Flags and internal consistency

1. [!] **Table 2 rows do not sum to the printed totals** [Z]: free 34.3 (printed 34.6), total 48.7 (printed 48.0).
2. [!] **Free fraction: text vs table.** Text: free amino acids are "35-40% of the total protein" [V]; Table 2 gives
   34.6 / 48.0 = 72 % of total amino acids [Z]. The 35-40 % may be of N × 6.25 protein (e.g. 69 g in Table 1:
   34.6 / 69 = 50 %) — still not 35-40 %.
3. [!] **Per-residue free fractions contradict the table** [Z]: text says ≥ 85 % of Leu, Ala, Met, Phe are free —
   table gives 92, 84, 78, 78 %; text says only 14-37 % of Asp, Gly, Arg are free — table gives 39, 50, 72 %.
   The text and the table come from different sources.
4. [!] **Thiamine**: Table 1 dried YE 10-12 mg/100 g vs text "about 3 mg" per 100 g standard YE (ref 27), a ×3-4
   disagreement inside the review. Jacob 2019 measures 5.2-7.5 mg/100 g dw — between the two.
5. [!] **Total glutamate**: abstract "10 % w/w total" vs Table 2 total Glu 8.1 g/100 g DM.
6. [!] **"Cystine" only** in Table 2: free cysteine (the Maillard-relevant thiol) is not given; a free "cystine" of 0.3
   g/100 g DM (≈ 12.5 µmol/g cystine, 25 µmol/g Cys-equivalent [Z, MW 240.3]) is comparable with Münch 1997's CYE free Cys 2.13 mg/g (17.6 µmol/g).
7. [!] Cross-check with Jacob 2019 (Table 2 sums, mg/g dw) [Z]: industrial autolysate free 433 / total 633; lab
   autolysis 306 / 447. Tomé's 346 / 480 sits between, consistent with ref 31 being one source, but unverified.

## 6. Existing records in the repo that cite this paper

`grep -rIn` 2026-11-03 for `tome2021`, `acsfoodscitech.0c00131`, "Tomé", "Tome 2021" across `data/ docs/ src/ scripts/
tasks/`: **only `data/keys/pdf_dossier_aliases.yml`** (pending line "Tome 2021 ACS Food Sci Technol
(10.1021/acsfoodscitech.0c00131): review of yeast extracts (composition, flavouring); background only": **matches**;
"background only" is a fair verdict). No ledger, registry or code record cites it.

**What it adds to existing yeast-extract composition.** Against `munch1997_tab3_7_9_yeast_free_amino_acids` (15 free
amino acids, baker's YE) and the new `jacob2019` entries (sums only), Tomé adds an **18-residue free + total amino-acid
profile** (incl. Glu, Asp, Lys, Trp, Tyr — useful for the amine pool and the Strecker-aldehyde precursors) and a
thiamine range, both secondary. It adds nothing for GSH, ribose or nucleotides. The Podpora papers (refs 27, 28) are
the primaries worth fetching if a Polish-spent-yeast autolysate column is wanted.

## 7. Proposed ledger entries (aroma_other.yml; not written)

One entry, for the printed Table 2 (and the thiamine / sugar rows of Table 1), flagged as a review compilation.

```yaml
  # ---- tome2021 : Tome 2021 (review of yeast extracts; typical composition tables) (ingested 2026-11-03) ----
  - id: tome2021_tab1_tab2_yeast_extract_typical_composition
    dossier: tome2021_extraction.md
    paper: "Tome 2021 ACS Food Sci Technol 1:487 (review)"
    lab: "AgroParisTech / INRAE (Tome D.), review"
    system: "'typical' dried Saccharomyces cerevisiae yeast extract (Table 2); dried yeast, dried yeast extract, yeast-extract paste (Table 1)"
    medium: "dry product / paste as sold; no production method per column"
    temperature_c: []
    ph: {initial: null, measured_during: false, final: null}
    time_min: []
    species: [free amino acids (18), total amino acids (18), cystine, sugar, thiamine]
    kind: other
    location: "Table 2 (g/100 g DM; transcribed in dossier sec. 4a); Table 1 rows (per 100 g; sec. 4b)"
    status: digitised_in_dossier
    quant_class: unknown
    simulable: false
    condition_quote: "Typical Total and Free Amino Acid Content of Dried Saccharomyces cerevisiae Yeast Extract"
    flags: [matrix_precursor_composition, matrix_food, review_compilation, sources_not_per_value, row_sums_do_not_match_totals, text_contradicts_table, no_free_cysteine_only_cystine, no_ribose_no_nucleotides_no_gsh]
    note: "Secondary: Table 2 from refs 27, 28, 30, 31 (Podpora 2015/2016, Rozmierska 2001, Jacob, Hutzler & Methner 2019 EFRT 245:95) without per-value attribution. Free sum 34.6 (rows 34.3), total 48.0 (rows 48.7) g/100 g DM. Plausibility band only; key the primaries instead when ingested."
```

## 8. Supplementary material

None referenced (no Supporting Information). Nothing model-relevant is held back; the gaps are in the cited
primaries.
