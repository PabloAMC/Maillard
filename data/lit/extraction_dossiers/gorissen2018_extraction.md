# Gorissen et al. 2018 — EXTRACTION (amino-acid composition of 35 commercial protein isolates by one UPLC-MS/MS method: soy (7 products), pea (3 products) and 13 other sources, as g per 100 g raw powder and as % of N x 6.25 protein)
### One laboratory, one hydrolysis, one instrument across every commercial isolate — the cleanest cross-source lysine / arginine / methionine / cysteine comparison on disk, but supplier-pooled, without free amino acids, free sugars, moisture or ash.

**Source on disk:** `data/articles/gorissen2018.pdf` (owner's download, 2026-09-08). Read from the
scratchpad text layer (`gorissen2018.txt`, clean; Table 1 and Table 2 extract row-aligned and
were checked against their printed column sums). Two tables (Table 1, Table 2), four figures
(Figs 1-4: protein content, EAA, and per-amino-acid bars as % of protein, mean +/- SEM,
FIGURE-ONLY except where the text prints the value). Repo status before this dossier: the paper is
cited in `data/benchmarks/maillard_validation_benchmarks.md` §2.3 with "approximate" key values
(pea Lys ~7.2, Cys ~0.9, Met ~0.9; soy Lys ~6.4, Cys ~1.1, Met ~1.3 g/100 g protein; an Asn+Asp
column) and in `docs/protocols/{pea,soy}_matrix_meaty_benchmark.md` for sulfur amino acids —
none of those numbers is in this paper (flag 1). `data/species/protein_matrices.yml` carries the
amine pools from flour compositions (pea 0.524, soy 0.379 mmol/g protein).

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Protein content and amino acid composition of commercially available plant-based protein isolates" |
| Authors | Stefan H. M. Gorissen, Julie J. R. Crombag, Joan M. G. Senden, W. A. Huub Waterval, Jörgen Bierau, Lex B. Verdijk, Luc J. C. van Loon (corresponding) — NUTRIM, Maastricht University Medical Centre+; Clinical Genetics, MUMC+ |
| Venue | Amino Acids (2018), original article; received 19 Mar 2018, accepted 24 Aug 2018; open access CC BY 4.0. (The benchmarks file cites it as 50(12):1685-1695.) |
| DOI | 10.1007/s00726-018-2640-5 |
| Products | 35 commercial powders "presently commercially available as isolated protein powder suitable for application in human nutrition or animal feeds": oat 1, lupin 1, wheat 7, hemp 1, microalgae 1, **soy 7**, brown rice 1, **pea 3**, corn 3, potato 2, milk 1, whey 3, caseinate 1, casein 2, egg 1; plus human m. vastus lateralis from 10 volunteers. Suppliers listed as a pool (Agri Nutrition, Agridient, Avebe, Cargill, Chamtor, Cosucra, FrieslandCampina DMV / Domo, L.I. Frank, MRM, Roquette, Selecta, Tate & Lyle, Tereos, Volac, Vitablend, Wulro) — **no product is tied to a supplier or a value** (flag 2). Samples obtained and analysed Dec 2014 - Jun 2018. |
| Nitrogen factor | **6.25** for every source, deliberately ("to enable direct comparisons") |
| Naming | "raw material" = the powder as received; "% of total protein" = g amino acid per 100 g of N x 6.25 protein; ΣEAA = His + Ile + Leu + Lys + Met + Phe + Thr + Val (Trp not measured) |
| Compound registry | none (amino acids only; no volatiles) |

## 1. Why it matters

For Programme 7 the number wanted is the lysine (and arginine, cysteine, methionine) content of a
commercial pea or soy isolate per gram of protein, to replace the flour-derived amine pool. This
paper gives it for **the mean of 3 pea and 7 soy commercial isolates measured identically**:
lysine 4.7 and 3.4 g per 100 g powder, 5.9 and 4.6 % of N x 6.25 protein, i.e. **0.40 (pea) and
0.31 (soy) mmol lysine per g protein** — 23 % and 17 % BELOW the flour values the repo carries
and 25 % below Jaeger 2023's single-product values (0.539 / 0.410). Arginine 5.9 / 4.8 g/100 g
powder (0.42 / 0.37 mmol/g protein); methionine 0.3 / 0.3 and cysteine 0.2 / 0.2 g/100 g powder
(0.025-0.027 and 0.021-0.022 mmol/g protein, lower bounds because of the hydrolysis). The paper
measures nothing else the programme needs — no free amino acids, no sugars, no moisture, no ash
— and its absolute values run low for reasons it states itself (12-h acid hydrolysis without
oxidation; N x 6.25 overstating legume protein). Its worth to the repo is the spread: the pea
lysine per g protein from three commercial products and the soy from seven are pooled to one mean
each (no SD in the table; SEM only in the figures), which is the closest thing on disk to a
market band for the amine pool.

## 2. Methods as they matter to a model

- **Protein content.** ~10 mg powder in duplicate, Dumas combustion (vario MAX cube CN), N x
  6.25. Ranges by source (text): soy 61-91 %, pea 77-81 %, wheat 74-88 %, corn 58-75 %, potato
  77-83 %, whey 72-84 %, casein 67-78 %; means printed in the text for pea (80 %), potato (80 %),
  wheat (81 %), brown rice (79 %), hemp (51 %), lupin (61 %), oat (64 %), corn (65 %), egg
  (51 %), caseinate (86 %), muscle (84 %). The **soy mean is FIGURE-ONLY** (Fig. 1); it is
  recovered below as 74 % from lysine (3.4 g/100 g powder = 4.6 % of protein) and methionine
  (0.3 = 0.4 %) — both give 74-75 % (inferred; flag 3).
- **Hydrolysis (verbatim core):** "Approximately 6 mg of protein powder ... was hydrolyzed in 3 mL
  6 M HCl for 12 h at 110 C. After hydrolysis, samples were cooled down to 4 C ... HCl was
  evaporated under nitrogen stream and the dried amino acids were reconstituted in 5 mL water."
  "The acid hydrolysis was performed in the absence of oxygen and the hydrolyzation process was
  terminated after 12 h of incubation to minimize the reduction of cysteine and methionine.
  Although the acid hydrolysis is not optimal for all amino acids, we used this procedure for
  all protein samples to enable direct comparisons." No performic-acid oxidation, so cysteine
  and methionine are partly destroyed (flag 4); 12 h is half the conventional 24 h, so the
  hydrophobic Ile / Val / Leu bonds are incompletely cleaved (their values are low bounds too).
- **Quantification.** UPLC-MS/MS, underivatised amino acids (Waterval et al. 2009): 10 uL
  hydrolysate + 1500 uL 0.5 mM TDFHA + 10 uL stable-isotope-labelled internal standards; Acquity
  BEH C18 1.7 um 2.1 x 100 mm, 30 C, TDFHA / acetonitrile gradient, 650 uL/min; Quattro Premier
  XE, ESI+, MRM; six-point standards 31.25-500 uM. **Asparagine and glutamine converted to Asp
  and Glu, and tryptophan destroyed**, in the hydrolysis; **aspartic acid was not reported
  either** (Table 1 footnote: "Tryptophan, aspartic acid, asparagine, and glutamine were not
  measured"), so the amino-acid sums cover about 70-75 % of the protein (flag 5).
- **Replicates.** Protein content in duplicate per powder; amino acids: one hydrolysate per
  powder implied (not stated). Table 1 prints the mean over the products of a source with no
  dispersion; Figs 1-4 print mean +/- SEM across products (FIGURE-ONLY). For n = 1 sources the
  bar is a single powder.
- **Free amino acids, sugars, moisture, ash, lipid:** NOT MEASURED.

## 3. Tables re-typed

### Table 1. "Amino acid content of various dietary protein sources and human skeletal muscle" — "Values are presented in g per 100 g raw material. Tryptophan, aspartic acid, asparagine, and glutamine were not measured"

Soy and pea columns row by row; the other 14 columns summarised after the table.

| amino acid | Soy (n = 7) | Pea (n = 3) |
|---|---|---|
| Threonine | 2.3 | 2.5 |
| **Methionine** | **0.3** | **0.3** |
| Phenylalanine | 3.2 | 3.7 |
| Histidine | 1.5 | 1.6 |
| **Lysine** | **3.4** | **4.7** |
| Valine | 2.2 | 2.7 |
| Isoleucine | 1.9 | 2.3 |
| Leucine | 5.0 | 5.7 |
| ΣEAA | 19.9 | 23.6 |
| Serine | 3.4 | 3.6 |
| Glycine | 2.7 | 2.8 |
| Glutamic acid (incl. Gln) | 12.4 | 12.9 |
| Proline | 3.3 | 3.1 |
| **Cysteine** | **0.2** | **0.2** |
| Alanine | 2.8 | 3.2 |
| Tyrosine | 2.2 | 2.6 |
| **Arginine** | **4.8** | **5.9** |
| ΣNEAA | 31.9 | 34.4 |

Checks. Soy: EAA rows sum to 19.8 (printed 19.9), NEAA to 31.8 (31.9); total 51.7 g amino acids
per 100 g powder = 70 % of a 74 % protein content. Pea: 23.5 (23.6) and 34.3 (34.4); total 58.0 =
72.5 % of 80 %. The shortfall is the unmeasured Asp/Asn (~10-12 % of legume protein), Trp, and the
water of hydrolysis convention. One-decimal rounding makes cysteine 0.2 mean 0.15-0.25 (+/- 25 %).

Other columns (g/100 g raw material), lysine / methionine / cysteine / arginine only: oat 1.3 /
0.1 / 0.4 / 3.1; lupin 2.1 / 0.2 / 0.2 / 5.5; wheat 1.1 / 0.7 / 0.7 / 2.4; hemp 1.4 / 1.0 / 0.2 /
5.3; microalgae 3.6 / 0.0 / 0.1 / 3.4; brown rice 1.9 / 2.0 / 0.6 / 5.4; corn 1.0 / 1.1 / 0.3 /
1.7; potato 4.8 / 1.3 / 0.3 / 3.3; whey 7.1 / 1.8 / 0.8 / 1.7; milk 5.9 / 2.1 / 0.2 / 2.6;
caseinate 5.9 / 2.2 / 0.1 / 2.9; casein 4.6 / 1.6 / 0.1 / 2.1; egg 2.7 / 1.4 / 0.4 / 2.6; human
muscle 6.6 / 1.7 / 0.0 / 4.4. (Whey cysteine 0.8 g/100 g powder at 72-84 % protein = ~1.0 g/100 g
protein, against the 2.4-3.0 g/100 g protein of the literature — a direct display of how much
cysteine this hydrolysis loses; flag 4.)

### Values printed in the text as % of total protein (Figs 2-4 otherwise FIGURE-ONLY)

| item | Soy | Pea | where |
|---|---|---|---|
| Protein content, % of raw material | FIGURE-ONLY (range 61-91 %); 74 % inferred | 80 % (range 77-81 %) | Results "Protein content"; Fig. 1 |
| ΣEAA, % of protein | 27 | 30 | Results "Essential amino acid content"; Fig. 2 |
| **Lysine, % of protein** | **4.6** | **5.9** | Results "Amino acid profiles", Discussion; Fig. 4a |
| **Methionine, % of protein** | **0.4** | **0.4** | same; Fig. 4b |
| Leucine, % of protein | not printed (plant mean 7.1 +/- 0.8) | not printed | Fig. 3a |
| Plant-source means | lysine 3.6 +/- 0.6 %, methionine 1.0 +/- 0.3 % (n = 10 sources; +/- SEM) | | Abstract, Discussion |
| Animal-source means | lysine 7.0 +/- 0.6 %, methionine 2.5 +/- 0.1 % | | Abstract |
| WHO/FAO/UNU 2007 requirements quoted | lysine 4.5 %, methionine 1.6 %, leucine 5.9 % of protein | | Discussion |

Consistency: pea 4.7 g/100 g powder / 80 % = 5.9 % of protein (printed 5.9); soy 3.4 / 4.6 % =
73.9 % protein, 0.3 / 0.4 % = 75 % — the inferred soy protein content is 74-75 %.

### Table 2. "Representative amount of protein" (g of protein / g of raw material to supply 2.7 g leucine or 10.9 g EAA, the amounts in 25 g whey)

| source | matched for leucine: protein g / raw g | matched for ΣEAA: protein g / raw g |
|---|---|---|
| Soy | 40 / 55 | 40 / 55 |
| Pea | 38 / 48 | 37 / 46 |

(The ratio raw / protein gives the protein content used: soy 40/55 = 72.7 %, pea 38/48 = 79 % —
a second route to the FIGURE-ONLY soy protein content; 73-75 %.) Other rows: oat 47/73, lupin
52/86, wheat 45/55, hemp 54/105, microalgae 48/69, brown rice 37/47, corn 20/31, potato 33/41,
whey 25/32, milk 31/39, caseinate 30/35, casein 34/47, egg 39/77.

## 4. Numbers the repository can use

Molar masses: Lys 146.19, Arg 174.20, Met 149.21, Cys 121.16 g/mol. "Powder" = raw material as
received (moisture not measured, so not a dry-matter basis; flag 6). Per g protein = per g powder
/ protein fraction (pea 0.80 printed; soy 0.74 inferred), or directly from the printed % of
protein where the paper gives it. Values are means over 3 (pea) and 7 (soy) commercial products
with no printed SD.

| product | quantity | value +/- sd | unit as printed | mmol per g protein (arithmetic) | method | source | evidence class |
|---|---|---|---|---|---|---|---|
| Pea, mean of 3 commercial isolates | protein content | 80 (range 77-81) | % of raw material | — (N x 6.25) | Dumas, duplicate | text; Fig. 1 | measured (mean printed; SEM figure-only) |
| Soy, mean of 7 commercial isolates | protein content | FIGURE-ONLY (range 61-91) | % of raw material | — | Dumas | Fig. 1 | figure_only; **74 % inferred** from 3.4 g Lys / 4.6 % and 0.3 g Met / 0.4 % (and Table 2's 40/55) |
| **Pea** | **total lysine** | **4.7** (no sd) | g/100 g raw material | 47 mg/g powder / 146.19 = 0.322 mmol/g powder; / 0.80 = 58.75 mg/g protein, / 146.19 = **0.402 mmol/g protein** (paper's own 5.9 % -> 0.404) | UPLC-MS/MS after 6 M HCl 12 h | Table 1; Fig. 4a; text | measured (pooled mean of 3 products) |
| **Soy** | **total lysine** | **3.4** (no sd) | g/100 g raw material | 34 / 146.19 = 0.233 mmol/g powder; paper's 4.6 % of protein = 46 mg/g protein, / 146.19 = **0.315 mmol/g protein** | as above | Table 1; Fig. 4a; text | measured (pooled mean of 7 products; protein basis from the printed %) |
| Pea | total arginine | 5.9 | g/100 g raw material | 59 / 174.20 = 0.339 mmol/g powder; / 0.80 = 73.75 mg/g protein, / 174.20 = **0.423 mmol/g protein** | as above | Table 1 | measured |
| Soy | total arginine | 4.8 | g/100 g raw material | 48 / 174.20 = 0.276 mmol/g powder; / 0.74 = 64.9 mg/g protein, / 174.20 = **0.372 mmol/g protein** | as above | Table 1 | measured (protein basis inferred) |
| Pea | total methionine | 0.3 (0.4 % of protein) | g/100 g raw material | 3 / 149.21 = 0.020 mmol/g powder; 4 mg/g protein / 149.21 = **0.027 mmol/g protein** (0.025 via 0.3 / 0.80) | as above; no oxidation step | Table 1; Fig. 4b; text | measured, LOWER BOUND (flag 4); one-decimal rounding +/- 17 % |
| Soy | total methionine | 0.3 (0.4 % of protein) | g/100 g raw material | 0.020 mmol/g powder; **0.027 mmol/g protein** | as above | Table 1; Fig. 4b; text | measured, lower bound |
| Pea | total cysteine | 0.2 | g/100 g raw material | 2 / 121.16 = 0.0165 mmol/g powder; / 0.80 = 2.5 mg/g protein, / 121.16 = **0.021 mmol/g protein** | as above; no performic oxidation | Table 1 | measured, LOWER BOUND (flag 4); rounding +/- 25 %; compare Gao 2020's 0.061-0.074 mmol half-cystine / g protein by Ellman's for a pea isolate |
| Soy | total cysteine | 0.2 | g/100 g raw material | 0.0165 mmol/g powder; / 0.74 = 2.7 mg/g protein, / 121.16 = **0.022 mmol/g protein** | as above | Table 1 | measured, lower bound; compare the 0.100-0.114 mmol half-cystine / g protein of Ruan 2014 / Shimada 1988 in `protein_matrices.yml` |
| Pea / Soy | free amino acids | NOT MEASURED | — | — | — | — | — |
| Pea / Soy | free sugars, moisture, ash, lipid | NOT MEASURED | — | — | — | — | — |
| Pea / Soy | ΣEAA | 30 / 27 | % of protein | — | derived from Table 1 | text; Fig. 2 | measured |
| Plant isolates (10 sources) | lysine; methionine | 3.6 +/- 0.6; 1.0 +/- 0.3 | % of protein (mean +/- SEM over sources) | 0.25 +/- 0.04; 0.067 +/- 0.02 mmol/g protein | as above | Abstract | measured (cross-source mean) |

Comparison the matrix layer can print (arithmetic, not a fit): amine pool, mmol lysine per g
protein — repo (flour, USDA) pea 0.524 / soy 0.379; Jaeger 2023 (one commercial isolate each,
ion chromatography) 0.539 / 0.410; this paper (3 and 7 commercial isolates pooled, UPLC-MS/MS
after 12-h hydrolysis) 0.40 / 0.31. The two isolate sources disagree by 25 %, more than either's
stated uncertainty; §5 flag 7 gives the reasons and a band of **0.40-0.54 (pea) and 0.31-0.41
(soy) mmol/g protein** is what the evidence supports for a commercial isolate.

## 5. Flags

1. **The repo's existing citation of this paper is wrong in every number.**
   `data/benchmarks/maillard_validation_benchmarks.md` §2.3 attributes to Gorissen 2018 "pea
   isolate Lys ~7.2, Cys ~0.9, Met ~0.9; soy isolate Lys ~6.4, Cys ~1.1, Met ~1.3 g/100 g
   protein; Asn+Asp ~4.5 / ~5.0", and `docs/protocols/soy_matrix_meaty_benchmark.md` §3.2 /
   `pea_matrix_meaty_benchmark.md` cite "cysteine + methionine ~1.3 g/100 g protein (soy)" and
   "total cysteine 0.8-1.1 g/100 g protein (pea)" from it. The paper prints pea Lys 5.9 % / Cys
   0.25 % / Met 0.4 % and soy Lys 4.6 % / Cys 0.27 % / Met 0.4 % of protein, and did not measure
   aspartic acid or asparagine at all. Those files are outside this dossier's remit and were not
   edited; the benchmark row should be re-sourced (the values it quotes look like USDA / FAO
   flour tables, not this paper).
2. **No product is identifiable.** Seventeen suppliers are listed for 35 powders as a pool; the
   pea (n = 3) and soy (n = 7) values are means across unnamed products, so this paper cannot
   supply a "product | brand" row — only a market mean and (figure-only) SEM. The wide soy
   protein range (61-91 %) says the seven "soy isolates" include concentrates.
3. **Soy protein content is FIGURE-ONLY** and inferred here as 74-75 % from two printed ratios
   (lysine, methionine) and Table 2's 40 g protein / 55 g raw; the per-g-protein soy numbers for
   arginine and cysteine in §4 depend on that inference (lysine and methionine do not — the
   paper prints them as % of protein).
4. **Cysteine and methionine are lower bounds by method.** 6 M HCl hydrolysis without performic
   oxidation destroys a large, variable share of cysteine (as cystine / cysteic acid losses) and
   some methionine; the paper shortened hydrolysis to 12 h to limit this, at the cost of
   incomplete release of Ile / Val / Leu. The whey column (cysteine 0.8 g/100 g powder, ~1 g/100 g
   protein against a literature 2.4-3.0) shows a loss of roughly 60 %. Pea / soy cysteine here
   (0.021-0.022 mmol/g protein) is 3-5x below the Ellman half-cystine values already in
   `protein_matrices.yml` (Gao 2020 pea 0.061-0.074; Ruan / Shimada soy 0.10-0.11); do not
   replace those with these.
5. **Aspartic acid, asparagine, glutamine and tryptophan were not measured**; the amino-acid sums
   are 70-72 % of the N x 6.25 protein. Glutamic acid (12.4 / 12.9 g/100 g) includes glutamine.
   Nothing here for the asparagine pool (acrylamide route).
6. **No moisture, ash, lipid or sugar.** "g per 100 g raw material" is as-received powder, so the
   per-g-powder numbers are on a wetter basis than Jaeger 2023's dry-matter numbers (by the
   unmeasured ~5-8 % moisture); the per-g-protein numbers are basis-free.
7. **Absolute lysine runs 20-25 % below Jaeger 2023 and the flour tables.** Candidate causes:
   (i) different products; (ii) the one-decimal rounding (4.7 -> 4.65-4.75, 1 %); (iii) the 12-h
   hydrolysis (lysine is normally fully released, so a small effect); (iv) N x 6.25 inflating the
   protein denominator equally for both papers, so not a differential cause; (v) an
   underivatised-UPLC-MS/MS calibration bias with labelled internal standards in a hydrolysate
   matrix. None can be assigned from the paper. Carry both as the band 0.40-0.54 (pea) and
   0.31-0.41 (soy) mmol/g protein rather than choosing one.
8. **Acid hydrolysis regenerates part of any glycated lysine** (Amadori-lysine gives back roughly
   half its lysine, the rest as furosine / pyridosine), so "total lysine" by this method counts
   some already-blocked epsilon-amines as available; a minor upward bias for unheated isolates,
   larger for any heat-treated ingredient. The `amine_available_band` [0.4, 1.0] in the matrix
   table already covers it.
9. **Human muscle protein (84 %) is listed with proline 0.0 and cysteine 0.0** — plainly method
   artefacts for muscle, a reminder that the zeros in Table 1 (microalgae methionine 0.0) are
   detection failures, not absences.
