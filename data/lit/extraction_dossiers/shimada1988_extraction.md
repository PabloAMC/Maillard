# Shimada & Cheftel 1988 — EXTRACTION (commercial SPI, free SH 8 µmol/g and half-cystine ~100 µmol/g protein, gels 80-130 C)
### The reference soy SH / S-S balance by Ellman's DTNB in 6 M urea + 0.5 % SDS, with a DTT + gel-filtration half-cystine assay cross-checked by cysteic acid.

**Source on disk:** `data/articles/shimada1988.pdf` (owner's download, 2026-09-08). Read from the OCR text
layer in the scratchpad (`articles/shimada1988.txt`). The PDF's first page begins mid-references of the
preceding article (Fennema, Kramer, Labuza ...) and the paper itself starts on that page; the last page
carries the start of the next article (Zhuge, Posner & Deyoe, cottonseed gossypol). OCR renders "µmol"
as "pmol" and "1988" as "1Q08"/"1900" in running heads; both are silently corrected below. The paper has
NO numeric tables: all SH / half-cystine / solubility / firmness data are in Figures 1-8, and only the
values the authors state in the text are re-typed here. Repo status before this dossier: no soy site
density on file.

## 0. Identity

| field | value |
|---|---|
| Title | "Determination of Sulfhydryl Groups and Disulfide Bonds in Heat-Induced Gels of Soy Protein Isolate" |
| Authors | Kazuko Shimada (present address Yamaguchi Women's University) and Jean Claude Cheftel (Laboratoire de Biochimie et Technologie Alimentaires, Université des Sciences et Techniques du Languedoc, Montpellier) |
| Venue | J. Agric. Food Chem. 1988, 36, 147-153; received 17 Feb 1987, accepted 13 Jul 1987 |
| DOI | 10.1021/jf00079a038 (not printed in this 1988 scan; from the brief) |
| Naming | "free SH groups" = SH of protein soluble in the non-denaturing standard buffer; "total SH groups" = SH in 6 M urea + 0.5 % SDS (no reductant); "half-cystine" = total SH + 2 x S-S, measured after 10 mM DTT reduction and DTT removal. The paper never prints an S-S number; S-S = (half-cystine - total SH)/2 is left to the reader. |
| Cited by | Ruan et al. 2014 (`ruan2014_extraction.md`) for the NEM effect on soy gels |

## 1. Why it matters

This is the canonical soy-isolate SH / S-S balance: a commercial SPI with 8.0 µmol total SH per g
protein and ~100 µmol half-cystine per g protein, i.e. only 8 % of the cysteine pool is free thiol, the
rest is disulfide. It is the number the matrix layer needs for a `soy_isolate` entry and the natural
independent check on Ruan 2014's ~114 µmol/g total. It also states how the free-SH pool responds to
heating at gelling concentrations (down by up to 40 % at 80 C, down further to 115-120 C, back up at
130 C when S-S bonds start to break with 15-20 % half-cystine loss) and quantifies the accessibility
question: only 5.0 of the 8.0 µmol/g is seen in native buffer, and only 0.4-1.4 of 1.6-5.6 mol SH per
mol 11S is "surface" in the literature it reviews.

## 2. Methods as they matter to a model

- **Protein:** "Soy protein isolate (SPI) (Purina Protein 500E) was purchased from SIO, Boulogne sur
  Seine (France). It contained 90% protein (N x 6.25, db) and 4.6% moisture. The nitrogen solubility
  index (NSI) estimated according to AOCS (1970) was 68." Commercial, "partly heat denatured" (authors'
  own words in the Discussion). Protein in solutions by Lowry (Bensadoun & Weinstein modification).
- **Heating:** aqueous SPI dispersions, native pH 6.7 (or adjusted with 10 N NaOH to pH 7-10), in
  capped glass bottles, 30 min in a water bath (80-100 C) or autoclave (105-130 C; 130 C also 1 h), then
  tap-water cooled and held 4 C for 15-20 h before analysis. Protein concentrations: 0.9-18.2 %
  (concentration series, 80 C); 13 % (pH series and NEM series); 16 % (temperature series, 13 %
  also for firmness). No added salt; ionic strength is whatever the isolate brings.
- **NEM / succinimide:** 10 mM (pH series) or 13 mM (temperature series), 1 h room temperature before
  heating; "about 10 times the free SH concentration"; residual total SH after NEM "below 1 µmol/g of
  protein".
- **Solubilisation for SH assays:** samples brought to 0.2 % protein (0.1 g protein / 50 mL) in one of
  three buffers, Ultra-Turrax below 25 C 3 min, 20 000 g 15 min, supernatant assayed: (i) standard
  buffer = 0.086 M Tris - 0.09 M glycine - 4 mM Na2EDTA pH 8.0 -> "free SH"; (ii) standard buffer + 6 M
  urea + 0.5 % (17.3 mM) SDS -> "total SH"; (iii) as (ii) + 10 mM DTT -> half-cystine. Solubility of
  unheated SPI in (i)/(ii)/(iii): ~40 / 97 / 99 %; of heated gels: lower in (i) (down to ~15-40 %),
  ~90 % in (ii), ~100 % in (iii). ⚠ "free SH" is therefore expressed per g of TOTAL protein although
  only ~40 % (unheated) or less (gels) of the protein is in the assayed supernatant.
- **SH method, verbatim:** "SH groups were determined with use of 5,5'-dithiobis(2-nitrobenzoic acid)
  (DTNB) according to Ellman (1959) with some modifications. To a 3-mL aliquot of the protein
  supernatant in the standard buffer with or without denaturants was added 0.03 mL of Ellman's reagent
  solution (4 mg of DTNB/mL of standard buffer). After the solution was rapidly mixed and allowed to
  stand at room temperature for 15 min, absorbance was read at 412 nm. [...] A molar extinction
  coefficient of 1.36 x 10^4 M^-1 cm^-1 was used for calculating micromoles of SH/gram of protein."
  Reagent and protein blanks both used. DTNB charge: 0.12 mg = 0.30 µmol per 3 mL (100 µM) against
  ≤ 6 mg protein x 8 µmol/g = 0.048 µmol SH: reagent in ≥ 6-fold excess for total SH, ~2.5-fold for
  half-cystine (0.6 µmol/3 mL at 100 µmol/g), adequate.
- **Half-cystine by DTNB, verbatim:** "Protein gels homogenized in the standard buffer with urea, SDS,
  and DTT were centrifuged. The supernatants were incubated at room temperature for 6-8 h and then
  subjected to gel filtration on a Sephadex G-25 column (2.0 cm (i.d.) x 8.0 cm) equilibrated with the
  standard buffer containing denaturants, in order to remove DTT. [...] The protein fraction was
  collected, and the half-cystine content was determined with DTNB as described above." Control: BSA
  gave 0.8 mol free SH and 35.1 mol half-cystine per mol (literature: ~0.7 and 35), and reduced BSA
  held its SH for 20 h in air at 20 C.
- **Half-cystine by amino-acid analysis, verbatim:** "Half-cystine was determined as cysteic acid
  according to Moore (1963). SPI (5 mg) was treated with 2 mL of performic acid [...] at 0 C for 4 h.
  [...] protein hydrolysis was carried out in 2 mL of 6 N HCl at 110 C for 18 h. [...] The half-cystine
  content is the mean of two independent determinations."
- **Replicates:** SH and solubility in duplicate (range bars), unheated SPI in triplicate; firmness
  n = 5 (sd bars); amino-acid analysis n = 2.
- **Units in the repo:** all values are per g protein (Lowry); 1 µmol/g = 0.001 mmol/g. Powder basis if
  ever needed: x 0.90 (protein, db) x (1 - 0.046) (moisture) = x 0.859 g protein per g powder as sold.

## 3. Tables re-typed

No tables in the paper. Text-stated numbers:

### 3.1 Unheated SPI (Results, "Influence of Protein Concentration"; Fig. 2 legend values quoted in text)

| quantity | value | unit | where |
|---|---:|---|---|
| total SH (6 M urea + 0.5 % SDS) | **8.0** | µmol/g total protein | text; abstract says "initial 8 µmol/g of protein" |
| free SH (standard-buffer-soluble protein) | **5.0** | µmol/g total protein | text |
| half-cystine, Ellman after DTT | **"close to 100"** | µmol/g protein | text, Fig. 2 |
| half-cystine, amino-acid analysis | **104.0 ± 6.6** | µmol/g unheated soy protein | text (n = 2) |
| total SH / half-cystine | **~8 %** | — | text |
| protein solubility, standard / +urea,SDS / +DTT | ~40 / 97 / 99 | % | text |
| residual total SH after 10 mM NEM | < 1 | µmol/g protein | text |

### 3.2 Literature values for 11S globulin quoted by the authors (not their measurements)

| quantity | value | source cited |
|---|---|---|
| free SH, 11S | ~2 mol/mol = 6.3 µmol/g (MW 320 000) | Draper & Catsimpoolas 1978; Nakamura et al. 1984b |
| free SH, 11S (higher value) | 5.6 mol/mol = 15.7 µmol/g (MW 356 000) | Simard & Boulet 1978 |
| S-S bonds, 11S | 18-20 per mol | Draper & Catsimpoolas 1978; Kim & Kinsella 1986 |
| surface SH, 11S (phosphate pH 7.6) | 0.4-1.4 mol/mol = 1.1-3.9 µmol/g | Simard & Boulet 1978; Nakamura et al. 1984b |
| total SH, 11S (denaturant) | 1.6-5.6 mol/mol = 4.5-15.7 µmol/g | same |
| soy protein, dilute, phosphate pH 7.6, before / after 80 C 10 min | 6.2 / 5.5 µmol/g protein | Hashizume & Watanabe 1979 |
| BLG A, 1 % pH 7.0, 2 min: half-cystine loss | 3 % at 125 C; 10 % at 145 C | Watanabe & Klostermeyer 1976 |

### 3.3 Heated SPI: what the text states (all curves themselves are FIGURE-ONLY)

| series | statement in text |
|---|---|
| 80 C / 30 min, 0.9 % protein | total SH reduced, free SH (standard buffer) not significantly changed vs unheated |
| 80 C / 30 min, 0.9-18.2 % | total SH decreases with protein concentration, "40% maximum decrease"; free SH decreases "by 50% at maximum", tracking a 60 % maximum solubility loss; half-cystine constant ~100 |
| 80 C / 30 min, 13 %, pH 7-10 | total SH and free SH both fall with pH; half-cystine constant pH 7-10 |
| 16 %, pH 6.7, 80-120 C / 30 min | total SH and free SH fall with temperature up to 115 C; half-cystine "close to 100" |
| 16 %, 130 C, 30 min / 1 h | total SH and free SH INCREASE vs 115-120 C; half-cystine losses "about 15% after 30 min and 20% after 1 h" |
| gel firmness | no change 80-105 C; marked increase 115-120 C; drastic softening 130 C; NEM cuts firmness at all T except 130 C; succinimide has no effect; minimum gelling concentration 12-13 % at 80 C / 30 min |

## 4. Site densities the repository can use

Protein basis: per g protein (Lowry) as printed; no conversion assumption.

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN | conditions | source | evidence |
|---|---|---:|---|---:|---|---|---|
| soy protein isolate, commercial (Purina Protein 500E, 90 % protein N x 6.25 db, NSI 68) | total free SH (6 M urea + 0.5 % SDS, DTNB) | 8.0 (triplicate; range not printed) | µmol/g of total protein | **0.0080** | native, 0.2 % protein, pH 8.0 Tris-glycine-EDTA | Results §1; abstract | measured |
| same | native-buffer-accessible free SH | 5.0 | µmol/g of total protein | 0.0050 (lower bound on accessible SH; assayed on the ~40 % soluble fraction but divided by total protein) | native, standard buffer pH 8.0 | Results §1 | measured |
| same | half-cystine (Ellman after 10 mM DTT + G-25) | ~100 | µmol/g protein | **~0.100** | native and 80-120 C gels | text, Fig. 2/4/8 | measured (text gives "close to 100" only) |
| same | half-cystine (cysteic acid) | 104.0 ± 6.6 | µmol/g unheated soy protein | 0.104 | native | text | measured (n = 2) |
| same | S-S, native | (100 - 8.0)/2 = **46**; with 104.0: 48 | not printed | **0.046 (0.048)** | native | derived from the two rows above | inferred |
| same | total SH after 80 C / 30 min at 18.2 % | 8.0 x (1 - 0.40) = ~4.8 | not printed | ~0.0048 | 80 C, 30 min, pH 6.7 | "40% maximum decrease" | inferred from a stated percentage; exact point FIGURE-ONLY |
| same | S-S after 80 C / 30 min at 18.2 % | 46 + (8.0 - 4.8)/2 = ~47.6 | not printed | ~0.048 | same | arithmetic | inferred |
| same | half-cystine after 130 C | 100 x 0.85 = ~85 (30 min); ~80 (1 h) | not printed | ~0.085 / ~0.080 | 16 %, pH 6.7 | "losses about 15% ... 20%" | inferred |
| same | total SH after 115-120 C, and after 130 C | — | — | — | 16 % | Fig. 8 | FIGURE-ONLY |
| same | total SH vs pH 7-10 | — | — | — | 13 %, 80 C | Fig. 4 | FIGURE-ONLY |

Accessible fraction: 5.0 / 8.0 = 62 % of the total free SH is titratable without denaturant (on a
total-protein basis); on the literature 11S numbers the surface fraction is 25-100 % of the total
(0.4-1.4 of 1.6-5.6 mol/mol). Free SH as a fraction of half-cystine: 8 % (paper's own figure).

## 5. Flags

1. **No table, no printed sd for the central 8.0 / 5.0 / 100 values.** The 8.0 and 5.0 are means of
   triplicates whose range is drawn in Fig. 2 but not printed; the "100" is a reading the authors give
   as "close to 100". Only the cysteic-acid value carries a printed uncertainty (104.0 ± 6.6).
2. **S-S is never printed.** 46 µmol/g is this dossier's arithmetic, using the half-cystine definition
   the paper itself gives, [SH + (2 x S-S)].
3. **"Free SH" basis mismatch.** The 5.0 µmol/g is per g of TOTAL protein but was measured on the
   standard-buffer supernatant containing ~40 % of the protein; per g of soluble protein it would be
   ~12.5 µmol/g, which exceeds the total-SH value and shows the soluble fraction is SH-rich (or that
   the basis is as stated and the number is an amount, not a density). Use 5.0 only as "amount of SH
   titratable in native buffer per g of isolate protein".
4. **Commercial, partly heat-denatured isolate (NSI 68).** Native-state numbers are for this ingredient,
   not for undenatured soy globulins; Ruan 2014's lab isolate gives 7.5 µmol/g by a different reagent,
   so the two agree within method offsets.
5. **All heated-state values are FIGURE-ONLY** except the percentages quoted; the inferred rows in §4
   carry that percentage's rounding.
6. **Ionic strength unstated** for every heating series (aqueous dispersions of the isolate, no buffer).
7. **130 C reversal.** The free-SH pool RISES and half-cystine FALLS at 130 C (S-S breakdown with H2S
   release is the authors' reading). A model that only oxidises SH with temperature will be wrong above
   ~120 C for soy; this is a directional hold-out shape for the sulfur lane.
8. **DOI** is not printed in the scan; taken from the brief.
