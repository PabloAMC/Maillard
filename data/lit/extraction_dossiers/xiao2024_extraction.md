# Xiao, Hu, Kumar & Li 2024 — EXTRACTION (book chapter: TNBS free-amino and DTNB free / total SH protocols, with example values for SPI, PPI and wheat gluten)
### The Li-lab protocol behind Shen 2022 and Xiao 2025: one clean SH / S-S row each for a soy and a pea isolate, and a free-amino scale that is 19-25x the lysine content.

**Source on disk:** `data/articles/Xiao2024.pdf` (owner's download, 2026-09-08). The PDF is the WHOLE book
*Plant-Based Proteins: Production, Physicochemical, Functional, and Sensory Properties* (Methods and
Protocols in Food Science, ed. Yonghui Li, Springer; 442 pages). Read from the pdftotext layer in the
scratchpad (`articles/Xiao2024.txt`, 17 040 lines); the chapter is lines 4684-5144 of that file =
**book pages 123-132** (Chapter 10). The text layer of Tables 1-3 is clean (the "±" renders as `/C6`,
the degree sign as `/C14`); all three re-typed below and Table 3's S-S column re-derived from its own
two SH columns. Fig. 1 (leucine standard curve) is figure-only. Repo status before this dossier:
`protein_matrices.yml` carries soy and pea free-thiol / disulfide densities from Ruan 2014, Shimada
1988, Gao 2020, Shen 2022, Chihi 2016; **no amine density for either isolate**.

## 0. Identity

| field | value |
|---|---|
| Title | "Free Amino and Sulfhydryl Group Content" (Chapter 10) |
| Authors | Ruoshi Xiao, Ruijia Hu, Nandan Kumar, Yonghui Li (Kansas State University, Grain Science and Industry) |
| Venue | in Y. Li (ed.), *Plant-Based Proteins: Production, Physicochemical, Functional, and Sensory Properties*, Methods and Protocols in Food Science, Springer Science+Business Media, pp. 123-132; the chapter footer says "© The Author(s) ... 2025" although the volume is catalogued as 2024 |
| DOI | 10.1007/978-1-0716-4272-6_10 |
| Naming | the Abstract calls TNBS "2,4-dinitrobenzene-1-sulfonic acid (DNBS)"; the reagent actually listed and used (§2.1, §3.1.1) is 2,4,6-trinitrobenzenesulfonic acid (TNBS). "Free SH", "Total SH" and "SS bond" as in Table 3; the "Gluten—Lab" row is a wheat gluten prepared in the authors' laboratory, "Gluten—Commercial" a purchased vital wheat gluten |
| Kin | Shen 2022 (Food Chem 385, 132687) and Xiao 2025 (Food Hydrocolloids 163, 111119) use the same two assays; Shen 2022's method parent (Shen & Li 2021) is not on disk, so this chapter is the closest thing to a printed method for the lab's earlier free-NH2 and free-SH numbers |

## 1. Why it matters

The matrix layer (`src/kinetic_core/matrix_sites.py`, `data/species/protein_matrices.yml`,
`results/validation/matrix_sites_prereg.md` §5) charges `soy_isolate` and `pea_isolate` with measured
free-thiol and disulfide densities but has no amine density for either, so neither isolate binds an
aldehyde or HMF. This chapter is the first paper on disk that prints, for the SAME soy and pea
isolates under the SAME protocol, a free-amino-group number (TNBS) AND free SH, total SH and S-S
(Ellman's in 8 M urea, total after β-mercaptoethanol reduction and TCA precipitation). The SH / S-S
rows are usable and sit inside or near the repo's bands. The free-amino rows are NOT usable as
lysine densities: 9.48 (SPI) and 9.90 (PPI) mmol/g protein are 25x and 19x the lysine content of
those proteins and exceed the total number of peptide bonds per gram of soy protein (7.8 meq/g,
this same book's Chapter 12, Table 1, p. 151). The chapter therefore (a) confirms that Shen 2022's
"8.44 mmol/g" is a lab-wide scale artefact of this TNBS protocol, not a one-off typo, and (b) still
leaves the amine pool unfilled.

## 2. Methods as they matter to a model

### 2.1 Free amino groups (TNBS)

- **Reagents (§2.1, §3.1.1):** 1 % (w/v) SDS; 0.2125 M phosphate buffer pH 8.20 ± 0.02 ("Adding
  4.5 mL of 0.2125 M ... NaH2PO4 to 100 mL of 0.2125 M ... Na2HPO4"); 0.1 % (v/v) TNBS made fresh
  and in the dark from the 5 % (v/v) methanolic stock; 0.1 N HCl; 0-2.4 mM L-leucine in 1 % SDS.
- **Sample (§3.1.2):** "Add 2 mg of sample (see Note 4) and 5 mL of 1% (w/v) SDS solution to a 50
  mL centrifuge tube to obtain a sample solution at a concentration of 4 mg/mL (see Note 5)."
  ⚠ 2 mg in 5 mL is 0.4 mg/mL, not 4 mg/mL; 4 mg/mL needs 20 mg. Xiao 2025 §2.4 states 4 mg/mL.
  Which mass was actually weighed is not recoverable from the text. Note 5: "The sample
  concentration could be adjusted according to the properties of the user's samples."
- **Standard (§3.1.3, Table 1):** 2.4 mM L-leucine = 15.74 mg in 50 mL 1 % SDS (15.74 mg / 131.17
  g/mol / 0.050 L = 2.40 mM, checks); dilutions to 2.0, 1.6, 1.2, 0.8, 0.4, 0 mM. So the
  calibration is in **leucine α-amino equivalents**, concentration expressed in the 0.5 mL aliquot.
- **Reaction (§3.1.4):** 0.5 mL sample solution + 4 mL phosphate buffer + 4 mL 0.1 % TNBS; 50 °C,
  200 rpm, 1 h, dark; + 8 mL 0.1 N HCl; 30 min dark at room temperature; centrifuge 10 min at
  8000 g, 20 °C; A340 of the supernatant. Total volume 16.5 mL. No heating step beyond 50 °C, no
  hydrolysis: the number is meant to be the primary amines (ε-NH2 of lysine + α-NH2 of chain
  termini) accessible to TNBS in an SDS dispersion of the intact protein. The chapter's own
  Introduction says so: "the determination of free amino groups in native proteins essentially
  requires quantification of the number of amino groups in lysine."
- **Calculation:** NOT printed. No equation converts A340 to "mmol/g protein"; the protein content
  of the three samples is not stated anywhere in the chapter. Whether "per g protein" is a
  per-powder number divided by an unstated protein fraction, or a loose label, cannot be told.
- **Range of the calibration, worked out:** the top standard (2.4 mM in a 0.5 mL aliquot) carries
  1.2 µmol of NH2. With 2 mg of sample in the aliquot (4 mg/mL) the curve tops out at **0.6 mmol/g
  sample**; with 0.2 mg (0.4 mg/mL) at 6 mmol/g. The reported 9.5-9.9 mmol/g protein lies above the
  calibrated range under either reading of the sample mass.

### 2.2 Free and total sulfhydryl, disulfide (Ellman's)

- **Buffers (§2.2):** Tris-Gly = 0.086 M Tris, 0.09 M glycine, 0.004 M EDTA, pH 8.0; the same with
  **8 M urea** for the sample dispersion. Ellman's reagent 4 mg/mL DTNB in Tris-Gly, fresh (Note 8).
  β-mercaptoethanol; 12 % (w/v) TCA. Method "combined the method of Beveridge [18] and Tang [19]".
- **Sample (§3.2.1):** 75 mg powder in 10 mL Tris-Gly-8 M urea (**7.5 mg powder / mL**), vortex,
  shake overnight at 200 rpm. Fine powder (Note 11); homogeneous suspension without clumps (Note 13).
- **Free SH (§3.2.2):** 1 mL suspension + 4 mL Tris-Gly + 0.05 mL Ellman's; 15 min dark; centrifuge
  8000 g, 4 °C, 10 min; A412. So the assay itself runs in ~1.6 M urea (1 mL of 8 M into 5.05 mL);
  the overnight 8 M urea soak is the unfolding step. Dilution factor D = 5.
- **Total SH (§3.2.3):** 1 mL suspension + 4 mL Tris-Gly + 0.05 mL β-mercaptoethanol, 1 h; + 10 mL
  12 % TCA, 1 h; centrifuge; wash pellet twice with 5 mL 12 % TCA; redissolve in 10 mL Tris-Gly; 4 mL
  + 0.04 mL Ellman's, 15 min dark; centrifuge; A412. D = 10. Note the redissolution is in Tris-Gly
  WITHOUT urea after TCA denaturation.
- **Equations (§3.2.4):** µmol SH/g = 73.53 × A412 × D / C, with C the sample concentration in
  mg/mL and 73.53 = 10^6 / 13 600 (ε412 of TNB^2- = 1.36 × 10^4 M^-1 cm^-1); "10^6 is for
  conversions from the molar basis to the µmol/mL basis and from mg solids to g solids". S-S
  (µmol/g) = (total SH − free SH) / 2. Blank correction (Note 15): A412 = A1 − A2 − A3 (sample,
  reagent blank, protein blank).
- **Basis:** the equation yields **µmol per gram of solids (powder)** because C is the powder
  concentration (7.5 mg/mL). Table 3 is headed "µmol/g protein". Either the authors divided by a
  protein fraction they do not print, or the header is loose. Treated below as "per g protein as
  printed; if actually per g powder, multiply by 1/f_p (about 1.1-1.2 for 85-90 % isolates)".
- **Replicates:** not stated in the chapter; ± values in Tables 2-3 are presumably SD of at least
  duplicates (the lab's papers say "at least duplicate"; Xiao 2025 prints n = 2).
- **Materials:** the SPI, PPI and commercial gluten are not identified (supplier, protein content,
  process); the lab gluten's preparation is not described.

## 3. Tables re-typed

### Table 1 (p. 127). "Preparation of standards for determination of free amino groups"

| tube | 2.4 mM L-leucine (mL) | 1 % SDS (mL) | resulting L-leucine (mM) |
|---|---|---|---|
| 1 | 2.5 | 0.5 | 2.0 |
| 2 | 2.0 | 1.0 | 1.6 |
| 3 | 1.5 | 1.5 | 1.2 |
| 4 | 1.0 | 2.0 | 0.8 |
| 5 | 0.5 | 2.5 | 0.4 |
| 6 | 0.0 | 3.0 | 0.0 |

(The 2.4 mM stock itself is the seventh point; §3.1.4 step 6 says "0-2.4 mM".)

### Table 2 (p. 128). "Free amino group content (mmol/g protein) of soy protein isolate (SPI), pea protein isolate (PPI), and wheat gluten"

| sample | free amino group (mmol/g protein) |
|---|---:|
| SPI | 9.48 ± 0.04 |
| PPI | 9.90 ± 0.05 |
| Gluten | 5.01 ± 0.03 |

### Table 3 (p. 130). "Free SH group (µmol/g protein), total SH group (µmol/g protein), and SS bond contents (µmol/g protein) of soy protein isolate (SPI), pea protein isolate (PPI), and wheat gluten"

| sample | free SH (µmol/g protein) | total SH (µmol/g protein) | SS bond (µmol/g protein) | (total − free)/2, re-derived |
|---|---:|---:|---:|---:|
| SPI | 4.74 ± 0.15 | 80.88 ± 1.09 | 38.07 ± 0.51 | 38.07 |
| PPI | 6.26 ± 0.17 | 63.92 ± 2.32 | 28.83 ± 1.24 | 28.83 |
| Gluten—Commercial | 2.86 ± 0.23 | 94.71 ± 2.41 | 45.92 ± 1.09 | 45.93 |
| Gluten—Lab | 2.83 ± 0.28 | 100.75 ± 2.17 | 48.96 ± 1.11 | 48.96 |

The S-S column is exactly Eq. 2 applied to the two SH columns (no independent S-S assay). Free SH as
a share of total SH: SPI 5.9 %, PPI 9.8 %, gluten 3.0 % / 2.8 %.

## 4. Site densities the repository can use

Basis assumption for every row: the printed "per g protein" is taken at face value; the
per-powder alternative (Eq. 1 with C = 7.5 mg powder/mL) would raise each number by 1/f_p with f_p
unstated (about 1.1-1.2x for a typical 85-90 % isolate). Native = dry commercial-type powder,
no heat step in the assay beyond 50 °C (TNBS) or room temperature (DTNB).

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN (arithmetic; basis) | conditions | source | evidence |
|---|---|---:|---|---|---|---|---|
| soy protein isolate (unidentified commercial-type SPI) | free SH | 4.74 ± 0.15 | µmol/g protein | **0.00474 ± 0.00015** (÷1000; per g protein as printed; 0.0052-0.0056 if per powder at f_p 0.85-0.90) | native, overnight 8 M urea, DTNB pH 8.0 | Table 3, p. 130 | measured, basis as printed |
| same | total SH (after β-ME reduction) | 80.88 ± 1.09 | µmol/g protein | 0.0809 ± 0.0011 (= half-cystine; × 121.16 mg/mmol = 9.8 mg Cys / g protein = 0.98 g / 100 g protein) | as above | Table 3 | measured |
| same | S-S | 38.07 ± 0.51 | µmol/g protein | **0.0381 ± 0.0005** (derived by the paper as (total − free)/2, same as the repo's Shimada / Ruan rows) | as above | Table 3 | derived from two measured SH values |
| same | free amino groups | 9.48 ± 0.04 | mmol/g protein | **do not use** (flag 1); leucine-equivalent TNBS scale, 25x the lysine content | 1 % SDS, TNBS pH 8.2, 50 °C 1 h | Table 2, p. 128 | measured, scale implausible |
| pea protein isolate (unidentified commercial-type PPI) | free SH | 6.26 ± 0.17 | µmol/g protein | **0.00626 ± 0.00017** (per g protein as printed; 0.0070-0.0074 if per powder) | native, as above | Table 3 | measured, basis as printed |
| same | total SH | 63.92 ± 2.32 | µmol/g protein | 0.0639 ± 0.0023 (7.7 mg Cys / g protein = 0.77 g / 100 g) | as above | Table 3 | measured |
| same | S-S | 28.83 ± 1.24 | µmol/g protein | **0.0288 ± 0.0012** | as above | Table 3 | derived |
| same | free amino groups | 9.90 ± 0.05 | mmol/g protein | **do not use** (flag 1); 19x the lysine content | as above | Table 2 | measured, scale implausible |
| wheat gluten, commercial | free SH / total SH / S-S | 2.86 ± 0.23 / 94.71 ± 2.41 / 45.92 ± 1.09 | µmol/g protein | 0.00286 / 0.0947 / **0.0459** | native | Table 3 | measured / measured / derived |
| wheat gluten, lab-prepared | free SH / total SH / S-S | 2.83 ± 0.28 / 100.75 ± 2.17 / 48.96 ± 1.11 | µmol/g protein | 0.00283 / 0.1008 / **0.0490** | native | Table 3 | measured / measured / derived |
| wheat gluten (both) | free amino groups | 5.01 ± 0.03 | mmol/g protein | do not use; 27x the lysine content of whole-wheat protein (0.186 mmol/g) and more for gluten | | Table 2 | measured, scale implausible |
| any | amine after heating | not measured | — | — | no heated sample in the chapter | — | absent |

Comparison with the repo table (`protein_matrices.yml`, native isolates, mmol per g protein):

| matrix | quantity | repo centre (band) | this chapter | ratio |
|---|---|---|---|---|
| soy isolate | free thiol | 0.0078 (0.0075-0.0080) | 0.0047 | 0.61x — below the band |
| soy isolate | disulfide | 0.050 (0.046-0.053) | 0.0381 | 0.76x — below the band |
| soy isolate | half-cystine (total SH) | ~0.10-0.11 (Shimada, Ruan; free thiol is "7-8 % of half-cystine") | 0.081 | 0.75x |
| pea isolate | free thiol | 0.0159 (0.0021-0.0174) | 0.0063 | 0.39x — inside the band, between Chihi's globulins and Gao's whole isolate |
| pea isolate | disulfide | 0.0257 (0.0042-0.0297) | 0.0288 | 1.12x — inside the band |
| pea isolate | half-cystine | 0.062-0.073 (Gao 2020: 52-62 µmol/g powder at 83-85 % protein) | 0.064 | agrees |

Both isolates' S-S rows here rest, like the repo's, on half-cystine minus free thiol over two; no
paper on disk yet prints a direct S-S assay of a native isolate.

## 5. Flags

1. **The free-amino numbers cannot be lysine densities.** Lysine content of the parent proteins,
   from this same book's Chapter 1, Table 2 (USDA / FAO, mg per g protein): pea (split, USDA-16085)
   76.6 → **0.524 mmol/g protein**; soybean flour (USDA-16115) 55.4 → **0.379 mmol/g**; whole-wheat
   flour 27.2 → 0.186 mmol/g (Lys 146.19 g/mol). Add ~0.02-0.05 mmol/g of N-terminal α-amines. A
   TNBS number on an intact protein should therefore be ≤ 0.4-0.6 mmol/g protein, and in practice
   lower because buried lysines under-react. Table 2 prints 9.48, 9.90 and 5.01: 25x, 19x and 27x
   the lysine content, and larger than the total peptide-bond count of soy protein (h_tot = 7.8
   meq/g, Chapter 12 Table 1, p. 151). The chapter's own Introduction defines free amino groups of
   a native protein as essentially the lysine ε-NH2, so the table contradicts its own premise.
   Together with Shen 2022 (8.44 mmol/g for a lab PPI) and Xiao 2025 (3.1-6.9 mmol/g protein for
   blends, 3.07 for gluten) this is a **systematic scale of the Li-lab TNBS protocol**, not a typo.
   The repo's amine pool for soy and pea must NOT be taken from this chapter.
2. **A candidate mechanism, conjecture only:** the calibration is in mM within the 0.5 mL aliquot;
   if the equivalent concentration were multiplied by the full 16.5 mL reaction volume instead of
   the 0.5 mL aliquot, every value would be inflated 33x, which would put SPI at 0.29, PPI at 0.30
   and gluten at 0.15 mmol/g protein — all ordinary TNBS values for intact proteins (0.6-0.8x
   lysine). The text does not print the calculation, so this is not a correction and no repo
   number may be derived from it. The relative order (PPI ≥ SPI > gluten) is at least consistent
   with lysine, but the PPI/SPI ratio (1.04) does not reproduce the lysine ratio (1.38).
3. **Sample-mass sentence is internally inconsistent** (2 mg in 5 mL ≠ 4 mg/mL, §3.1.2); the
   reported values lie above the top of the leucine curve under either reading (§2.1).
4. **"Per g protein" is not demonstrated.** Eq. 1 produces µmol per g of powder (C = 7.5 mg
   powder/mL; "from mg solids to g solids"), and no protein content, nitrogen factor or protein
   assay is given for the SPI, PPI or glutens. §4 therefore carries both readings; the per-powder
   reading raises the SH / S-S densities by ~10-20 %.
5. **Isolates are not identified** (no supplier, lot, process, protein content, moisture). The
   repo's soy row is built from alcohol-washed lab and commercial isolates that agree (0.0075-0.0080);
   this chapter's 0.0047 is 40 % lower, which is within what supplier and drying history do to
   free SH (Gao 2020 shows 20 % movement from extraction pH alone).
6. **Total SH protocol detail:** the TCA-precipitated, reduced protein is redissolved in Tris-Gly
   WITHOUT urea before the second Ellman's step; incomplete redissolution would under-read total
   SH and hence S-S. The soy half-cystine here (0.081 mmol/g) is 0.75x Shimada's, consistent with
   some under-recovery or a different isolate; cannot be separated.
7. **No heated sample, no S-S after heating, no amine after heating.** The chapter is a protocol;
   it gives the native baseline only.
8. **Replicate number not stated**; ± presumably SD of duplicates as in the lab's papers.
9. Free SH was assayed after an overnight 8 M urea soak but in ~1.6 M urea during the 15 min
   colour step. This is the same protocol as Xiao 2025 and (by the lab's reference chain) Shen 2022,
   so those three free-SH numbers are on one footing; Gao 2020's (8 M urea throughout) and
   Chihi 2016's (1.5 M GdnCl) are not identical footings.

## 6. Other chapters of the book with pea / soy isolate numbers the repository could use later (not extracted)

- **Chapter 1, Du, Rajpurohit, Kumar & Li, "Overview of Plant-Based Proteins", pp. 3-20** — Table 1
  (pp. 5-6) proximate composition of seeds/flours (pea 23.1 % protein, soybean full-fat flour
  38.6 %); Table 2 (pp. 7-11) full amino-acid composition in mg per g protein for pea (Lys 76.6,
  Cys 11.8, Met 8.4; USDA-16085, p. 9) and soybean flour (Lys 55.4, Cys 19.8, Met 12.6; USDA-16115,
  p. 10), plus wheat, rice, other pulses and oilseeds; Table 3 (p. 13) Osborne-class shares.
  Seed/flour values, not isolates, but the only lysine numbers in the book.
- **Chapter 9, Liu, Zhao, Gupta & Carrillo, "Quantification of Crude and Soluble Protein Content",
  pp. 107-122** — Table 1 (pp. 113-114) nitrogen-to-protein conversion factors: pea 5.44 / 5.40,
  soybean 5.71 / 5.44 / 5.52, wheat 5.83 / 5.49 / 5.33 / 5.75 (vs the 6.25 used by Gao 2020 and
  Chen 2022). No isolate protein contents printed.
- **Chapter 12, Huang, Zhang, Zhou, Zhang & Sui, "Protein Digestibility Through In Vitro
  Gastrointestinal Digestion", pp. 143-154** — Table 1 (p. 151) h_tot = 7.8 meq peptide bonds per g
  protein for soy (8.3 gluten), the α/β constants for OPA degree-of-hydrolysis; Fig. 4 SPI DH curves
  are figure-only.
- **Chapter 21, Zhang & Xiao, "Surface Hydrophobicity", pp. 279-286** — SPI at 0.1 % w/v in 0.01 M
  phosphate pH 7.0: ANS index S0 = 779.91 native, 1708.62 after 95 °C / 10 min; bound BPB 23.76 µg
  native, 56.60 µg heated (p. 283). Instrument-specific scale, but a native-vs-heated pair for SPI.
- **Chapter 26, Hong, Kwon, Xu & Li, "Emulsifying Properties", pp. 323-334** — Table 1 (p. 329):
  PPI EA 71.8 %, ES 62.0 %, EAI 566 m²/g, ESI 105 min, and PPI:guar blends; Table 2 (p. 332)
  emulsion capacity 802.5 mL oil / g protein. Functional, not compositional.
- **Chapter 27, Rajpurohit, Chen & Li, "Water and Oil Holding Capacity", pp. 335-344** — Tables 1-3
  (pp. 337-339): commercial SPI, PPI and vital wheat gluten WHC (SPI 4.17-4.29, PPI 3.95-4.16,
  VWG 1.50-1.63 g/g; pH 5.5 and AACC variants), with powder moisture contents (SPI 5.99 %, PPI
  7.26 %, VWG 5.93 %, Table 3, p. 339). Functional; the moisture numbers are the only isolate
  composition data there.
- **Chapter 11, Chang & Chen, "Amino Acid Composition", pp. 133-142** — HPLC (AccQ-Tag) protocol
  only; no example pea / soy composition table found in the text layer.
- Chapters 5-6 (wet extraction; Osborne fractionation, pp. 61-74) describe procedures without
  isolate composition tables in the text layer; Chapter 32 (twin-screw extrusion, pp. 405-432) is
  process-only.
